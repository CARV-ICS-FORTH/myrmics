/*****************************************************************************

 Copyright 2006 Sandia Corporation. Under the terms of Contract
 DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 retains certain rights in this software.

 This program is free software; you can redistribute it and/or
 modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation; either version 2
 of the License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, 
 Boston, MA  02110-1301, USA.

 name:		mpi_overhead.c

 purpose:	A benchmark to calculate the host overhead associated with 
 		sending and receiving MPI messages.

*****************************************************************************/
#include <kernel_toolset.h>
#include <arch.h>
#include <fmpi.h>

#define AUTO_ITERATIONS       0
#define MAX_ITERATIONS        1000
#define MIN_ITERATIONS        10
#define MIN_WORK              1      /* # of work loop iterations */
#define WORK_FACTOR_1         2.0
#define TIMING_RESOLUTION     10 // clock cycles

#define OVRCV		      0
#define OVSND		      1
#define OVNOP		      2
#define OVDONE		      3


float smb_overhead_wtime(unsigned int *tmr) {
#ifdef ARCH_MB
  unsigned int new;
  unsigned int diff;

  new = ar_glob_timer_read();

  if (new > *tmr) {
    diff = new - *tmr;
  }
  else {
    diff = 0xFFFFFFFF - (*tmr - new);
  }

  *tmr = new;

  return (float) (diff / TIMING_RESOLUTION);
#else
  return 0.0;
#endif
}

void smb_overhead_show_partials(int label, float iter_t, float base_t) {

  kt_printf( "work      iter_t       base_t\r\n");

  kt_printf( "%-8d %-13.3f %-13.3f\r\n", label, iter_t, base_t);
}

void smb_overhead_show_results(int verbose, int nohdr, int data_size, 
                               int iterations, int work, float iter_t, 
                               float work_t, float overhead, float base_t) {
  float availability = 100.0 * (1.0 - overhead / base_t);

  if (!nohdr) {
      kt_printf(
          "msgsize iterations   iter_t  work_t  overhead   base_t avail(%)\r\n");
  }

  kt_printf( "\
%-8d\
%-12d \
%-12.3f \
%-12.3f \
%-12.3f \
%-12.3f \
%-12.1f\n\
", data_size, iterations, iter_t, work_t, overhead, base_t, 
availability);
}


// direction: 0 = send, 1 = recv        default: 0
// iterations                           default: 0 (= auto)
// data_size: bytes                     default: 8
// threshold: ratio(?)                  default: 1.5
// base_threshold: ratio(?)             default: 1.02
// nohdr: 0 = print headers, 1 = don't  default: 0
// verbose: 0 = don't print partials    default: 0
int smb_overhead_mpi(int direction, int iterations, int data_size, 
                     float threshold, float base_threshold, int nohdr, 
                     int verbose) {

  MPI_Status rstatus, sstatus;
  MPI_Request rrequest, srequest;
  void *msg = NULL;
  unsigned char *tmsg;
  float time1, base_time, iter_time, work_time, accum_time, overhead, 
    work_factor;
  int count, size, rank, dest_node, iter,
    work, new_work, last_work, accum_count, continue_flag;
  struct {
      int command;
      int iterations;
  } *message;
  volatile int x;
  volatile float y, a = 1.0, b = 1.0;
  unsigned int tmr;


  /**************************************************************
    Initialize MPI
  **************************************************************/
  MPI_Comm_size( MPI_COMM_WORLD, &size );
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  if ( size % 2 ) {
      if ( rank == 0 ) {
          kt_printf( "ERROR: This program requires # processors be a \
multiple of 2\n");
      }
      return 1;
  }


  /**************************************************************
    Allocate messaging resources
  **************************************************************/
  if ((message = kt_malloc(2 * sizeof(int))) == NULL) {
    kt_printf( "Unable to allocate memory\n");
    return 1;
  }
  if (msg == NULL) {
    if ((msg = kt_malloc(data_size)) == NULL) {
      kt_printf( "Unable to allocate memory\n");
      return 1;
    }
    for (tmsg = (unsigned char *)msg, count = 0; count < data_size; 
    ++tmsg, ++count)
      *tmsg = count;
  }

  if (!rank && verbose) {
    if (direction == 0) // send
      kt_printf( "Calculating send overhead\n");
    else
      kt_printf( "Calculating receive overhead\n");
    if (iterations != AUTO_ITERATIONS)
      kt_printf( "Using %d iterations per test\n", iterations);
    else
      kt_printf( "Iterations are being calculated automatically\n");
    kt_printf( "Message size is %d bytes\n", data_size);
    kt_printf( "Overhead comparison threshold is %f (%.1f%%)\n", 
            threshold, 100.0 * (threshold - 1.0));
    kt_printf( "Base time comparison threshold is %f (%.1f%%)\n", 
            base_threshold, 100.0 * (base_threshold - 1.0));
    kt_printf( "Timing resolution is %d cc\n", TIMING_RESOLUTION);
  }

  /* autocalculate the number of iterations per test */
  if (iterations == AUTO_ITERATIONS) {
    if      (data_size < 64*1024)     iterations = MAX_ITERATIONS;
    else if (data_size < 8*1024*1024) iterations = MAX_ITERATIONS / 10;
    else                              iterations = MAX_ITERATIONS / 100;

    if (iterations < MIN_ITERATIONS) iterations = MIN_ITERATIONS;

    if (!rank && verbose)
      kt_printf( "Using %d iterations per work value\n", iterations);
  }

  MPI_Barrier( MPI_COMM_WORLD );


  /**************************************************************
    Timing Node
  **************************************************************/
  if ( rank < (size / 2) ) {     /* lower half of the rank */
    dest_node = rank + size / 2;
    work = MIN_WORK;
    work_factor = WORK_FACTOR_1;
    accum_time = 0.0;
    base_time = 0.0;
    accum_count = 0;
    do {
      message->iterations = iterations;
      if (direction == 0) // send
        message->command = OVRCV;
      else
	message->command = OVSND;
      if ( MPI_Send( (int *)message,
	             2,
		     MPI_INT,
		     dest_node,
		     1,
		     MPI_COMM_WORLD) != MPI_SUCCESS ) {
	kt_printf( "command send failed\n");
        return 1;
      }
      MPI_Barrier(MPI_COMM_WORLD);

      tmr = 0;
      smb_overhead_wtime(&tmr);
      if (direction == 0) { // send
        for (iter = 0; iter < iterations; ++iter) {

	  /* This barrier ensures the transfer from the previous 
	     iteration is complete and the recieve side is ready.
	     It will be subtracted from the total time as a 
	     part of the work time calculation. */
	  MPI_Barrier(MPI_COMM_WORLD);
	  
	  /* Send the data to the slave */
          if ( MPI_Isend( msg,
		          data_size,
		          MPI_BYTE,
		          dest_node,
		          2,
		          MPI_COMM_WORLD,
			  &srequest ) != MPI_SUCCESS ) {
	    kt_printf( "Unable to send data\n");
            return 1;
          }

	  /* Work */
	  for (x = 0; x < work; ++x)
	    y = a * (float)x + b;

	  /* wait for the send to complete */
	  MPI_Wait(&srequest, &sstatus);
        }
      }

      else { /* direction == recv */
        for (iter = 0; iter < iterations; ++iter) {

          /* Slave sends data */
          if ( MPI_Irecv( msg,
                          data_size,
                          MPI_BYTE,
	                  dest_node,
	                  2,
	                  MPI_COMM_WORLD,
	                  &rrequest ) != MPI_SUCCESS ) {
	    kt_printf( "Error receiving data\n");
            return 1;
          }

	  MPI_Barrier(MPI_COMM_WORLD);

	  /* Work */
	  for (x = 0; x < work; ++x)
	    y = a * (float)x + b;

	  /* wait for the send to complete */
	  MPI_Wait(&rrequest, &rstatus);
	}
      }
      time1 = smb_overhead_wtime(&tmr) / iterations;

      if (work == MIN_WORK) 
	base_time = time1;

      last_work = work;
      iter_time = time1;

      /* check to see if we're past the knee of the curve, 
         i.e. time1 is rising */
      if (time1 > (base_time * threshold)) {
        if (verbose)
          smb_overhead_show_partials(work, time1, base_time);
	break;
      }

      /* low pass filter the flat part of curve to determine 
         the base time */
      if ( time1 < (base_time * base_threshold) ) {
        accum_time += time1;
        base_time = accum_time / ++accum_count;
      }

      if (verbose)
        smb_overhead_show_partials(work, time1, base_time);

      if (work == 0) 
        new_work = 2;
      else
        new_work = work * work_factor;

      if (new_work == work)
        break;
      else
        work = new_work;

    } while(1);

    /* Work includes everything extra required for this benchmark
       i.e. the timer logic and the barrier */
    message->iterations = iterations;
    message->command = OVNOP;
    if ( MPI_Send( (int *)message,
	           2,
		   MPI_INT,
		   dest_node,
		   1,
		   MPI_COMM_WORLD) != MPI_SUCCESS ) {
      kt_printf( "command send failed\n");
      return 1;
    } 
    MPI_Barrier(MPI_COMM_WORLD);

    tmr = 0;
    smb_overhead_wtime(&tmr);
    for (iter = 0; iter < iterations; ++iter) {
      MPI_Barrier(MPI_COMM_WORLD);
      for (x = 0; x < last_work; ++x)
        y = a * (float)x + b;
    }
    work_time = smb_overhead_wtime(&tmr) / iterations;

    /* Overhead is the iteration time minus the work time */
    overhead = iter_time - work_time;

    if (!rank)
      smb_overhead_show_results(verbose, nohdr, data_size, iterations, 
                                last_work, iter_time, work_time, overhead, 
                                base_time);

    /* terminate the session. */
    message->command = OVDONE;
    if ( MPI_Send( (int *)message,
                   2,
                   MPI_INT,
	           dest_node,
	           1,
	           MPI_COMM_WORLD) != MPI_SUCCESS ) {
      kt_printf( "command send failed\n");
      return 1;
    }
    MPI_Barrier(MPI_COMM_WORLD);

  } /* end of timing node logic */

  /**************************************************************
    Slave Node
  **************************************************************/
  else {
    dest_node = rank - size / 2;

    /* receive data until told otherwise */
    continue_flag = 1;
    do {


      /* get the command */
      if ( MPI_Recv(  (int *)message,
                      2,
                      MPI_INT,
	              dest_node,
	              1,
	              MPI_COMM_WORLD,
	              &rstatus ) != MPI_SUCCESS ) {
	kt_printf( "Error receiving command\n");
        return 1;
      }
      MPI_Barrier(MPI_COMM_WORLD);

      switch (message->command) {
        case OVRCV: /* receive the message */
	  for (iter = 0; iter < message->iterations; ++iter) {
	    /* prepost the receive to avoid an unexpected message */
            if ( MPI_Irecv( msg,
                            data_size,
                            MPI_BYTE,
	                    dest_node,
	                    2,
	                    MPI_COMM_WORLD,
	                    &rrequest ) != MPI_SUCCESS ) {
	      kt_printf( "Error receiving data from the client\n");
              return 1;
            }
	    MPI_Barrier(MPI_COMM_WORLD);
            MPI_Wait(&rrequest, &rstatus);
	  }
	  break;

        case OVSND: /* send the message */
	  for (iter = 0; iter < message->iterations; ++iter) {
	    MPI_Barrier(MPI_COMM_WORLD);
            if ( MPI_Send( msg,
                           data_size,
                           MPI_BYTE,
	                   dest_node,
	                   2,
	                   MPI_COMM_WORLD) != MPI_SUCCESS ) {
	      kt_printf( "Error sending data to the client\n");
              return 1;
            }
	  }
	  break;
  
        case OVNOP: /* just the barrier in the loop */
	  for (iter = 0; iter < message->iterations; ++iter) {
	    MPI_Barrier(MPI_COMM_WORLD);
	  }
	  break;

        case OVDONE: /* we're all done */
          continue_flag = 0;
	  break;
      }

    } while (continue_flag);

  } /* end of slave node */

  kt_free(msg);
  kt_free(message);

  return 0;
}
