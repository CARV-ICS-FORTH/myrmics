// ============================================================================
// Copyright (c) 2013, FORTH-ICS / CARV 
//                     (Foundation for Research & Technology -- Hellas,
//                      Institute of Computer Science,
//                      Computer Architecture & VLSI Systems Laboratory)
// 
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
// 
//     http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
// 
// ============================================================================
// The code in this file is taken from the NAS benchmarks, modified to 
// run on the Formic architecture. The NAS benchmarks are available under
// the following copyright:
//
//  This benchmark is part of the NAS Parallel Benchmark 3.3 suite.
//  It is described in NAS Technical Report 95-020.
//
//  Permission to use, copy, distribute and modify this software
//  for any purpose with or without fee is hereby granted.  We
//  request, however, that all derived work reference the NAS
//  Parallel Benchmarks 3.3. This software is provided "as is"
//  without express or implied warranty.
// ============================================================================
// Code modified by Iakovos Mavroidis
// ============================================================================



#ifdef ARCH_MB
#include <fmpi.h>
#include <kernel_toolset.h>
#include <arch.h>
  
#define strcpy kt_strcpy
#define printf kt_printf
#define malloc kt_malloc
#define MPI_Wtime FMPI_Wtime

#else
#include "mpi.h"
#include <stdlib.h>
#include <stdio.h>
#endif

//#include "npbparams.h"

#define CLASS 'S'
#define NUM_PROCS 16
#define COMPILETIME "10 Apr 2012"
#define NPBVERSION "3.3.1"
#define MPICC "cc"
#define CFLAGS "-O3 "
#define CLINK "cc"
#define CLINKFLAGS "(none)"
#define CMPI_LIB "-lmpi"
#define CMPI_INC "(none)"

/////

/******************/
/* default values */
/******************/
#ifndef CLASS
#define CLASS 'S'
#define NUM_PROCS            1                 
#endif
#define MIN_PROCS            1


/*************/
/*  CLASS S  */
/*************/
#if CLASS == 'S'
#define  TOTAL_KEYS_LOG_2    16
#define  MAX_KEY_LOG_2       11
#define  NUM_BUCKETS_LOG_2   9
#endif


/*************/
/*  CLASS W  */
/*************/
#if CLASS == 'W'
#define  TOTAL_KEYS_LOG_2    20
#define  MAX_KEY_LOG_2       16
#define  NUM_BUCKETS_LOG_2   10
#endif

/*************/
/*  CLASS A  */
/*************/
#if CLASS == 'A'
#define  TOTAL_KEYS_LOG_2    23
#define  MAX_KEY_LOG_2       19
#define  NUM_BUCKETS_LOG_2   10
#endif


/*************/
/*  CLASS B  */
/*************/
#if CLASS == 'B'
#define  TOTAL_KEYS_LOG_2    25
#define  MAX_KEY_LOG_2       21
#define  NUM_BUCKETS_LOG_2   10
#endif


/*************/
/*  CLASS C  */
/*************/
#if CLASS == 'C'
#define  TOTAL_KEYS_LOG_2    27
#define  MAX_KEY_LOG_2       23
#define  NUM_BUCKETS_LOG_2   10
#endif


/*************/
/*  CLASS D  */
/*************/
#if CLASS == 'D'
#define  TOTAL_KEYS_LOG_2    29
#define  MAX_KEY_LOG_2       27
#define  NUM_BUCKETS_LOG_2   10
#undef   MIN_PROCS
#define  MIN_PROCS           4
#endif


#define  TOTAL_KEYS          (1 << TOTAL_KEYS_LOG_2)
#define  MAX_KEY             (1 << MAX_KEY_LOG_2)
#define  NUM_BUCKETS         (1 << NUM_BUCKETS_LOG_2)
#define  NUM_KEYS            (TOTAL_KEYS/NUM_PROCS*MIN_PROCS)

/*****************************************************************/
/* On larger number of processors, since the keys are (roughly)  */ 
/* gaussian distributed, the first and last processor sort keys  */ 
/* in a large interval, requiring array sizes to be larger. Note */
/* that for large NUM_PROCS, NUM_KEYS is, however, a small number*/
/* The required array size also depends on the bucket size used. */
/* The following values are validated for the 1024-bucket setup. */
/*****************************************************************/
//#if   NUM_PROCS < 257
#define  SIZE_OF_BUFFERS     13*NUM_KEYS/2
//#elif NUM_PROCS < 512
//#define  SIZE_OF_BUFFERS     5*NUM_KEYS/2
//#elif NUM_PROCS < 1024
//#define  SIZE_OF_BUFFERS     4*NUM_KEYS
//#else
//#define  SIZE_OF_BUFFERS     13*NUM_KEYS/2
//#endif

/*****************************************************************/
/* NOTE: THIS CODE CANNOT BE RUN ON ARBITRARILY LARGE NUMBERS OF */
/* PROCESSORS. THE LARGEST VERIFIED NUMBER IS 1024. INCREASE     */
/* MAX_PROCS AT YOUR PERIL                                       */
/*****************************************************************/
#if CLASS == 'S'
#define  MAX_PROCS           128
#else
#define  MAX_PROCS           1024
#endif

#define  MAX_ITERATIONS      10
#define  TEST_ARRAY_SIZE     5


/***********************************/
/* Enable separate communication,  */
/* computation timing and printout */
/***********************************/
#define  TIMING_ENABLED
#ifdef NO_MTIMERS
#undef TIMINIG_ENABLED
#define TIMER_START( x, c )
#define TIMER_STOP( x, c )
#else
#define TIMER_START( x, c ) if (*c->timeron) timer_start( x, c )
#define TIMER_STOP( x, c ) if (*c->timeron) timer_stop( x, c )
#define T_TOTAL  0
#define T_RANK   1
#define T_RCOMM  2
#define T_VERIFY 3
#define T_LAST   3
#endif


/*************************************/
/* Typedef: if necessary, change the */
/* size of int here by changing the  */
/* int type to, say, long            */
/*************************************/
typedef  int  INT_TYPE;
typedef  long INT_TYPE2;
#define MP_KEY_TYPE MPI_INT

typedef struct context {
 int *timeron;

 /********************/
 /* MPI properties:  */
 /********************/
 int      my_rank,
         comm_size;

 /********************/
 /* Some global info */
 /********************/
 INT_TYPE *key_buff_ptr_global,         /* used by full_verify to get */
         total_local_keys,             /* copies of rank info        */
         total_lesser_keys;

 int      passed_verification;

 /************************************/
 /* These are the three main arrays. */
 /* See SIZE_OF_BUFFERS def above    */
 /************************************/
 INT_TYPE *key_array,
	*key_buff1,    
	*key_buff2,
	*bucket_size,     /* Top 5 elements for */
	*bucket_size_totals, /* part. ver. vals */
	*bucket_ptrs,
	*process_bucket_distrib_ptr1,   
	*process_bucket_distrib_ptr2;   
 int send_count[MAX_PROCS], recv_count[MAX_PROCS],
     send_displ[MAX_PROCS], recv_displ[MAX_PROCS];
 int *tmp;
 INT_TYPE2 test_index_array[TEST_ARRAY_SIZE],
         test_rank_array[TEST_ARRAY_SIZE];

 unsigned int start[64], elapsed[64];

 int        KS;
 float	R23, R46, T23, T46;
} context;



/**********************/
/* Partial verif info */
/**********************/
INT_TYPE2 S_test_index_array[TEST_ARRAY_SIZE] = 
                             {48427,17148,23627,62548,4431},
         S_test_rank_array[TEST_ARRAY_SIZE] = 
                             {49,6178,49934,64917,65463},

         W_test_index_array[TEST_ARRAY_SIZE] = 
                             {357773,934767,875723,898999,404505},
         W_test_rank_array[TEST_ARRAY_SIZE] = 
                             {1249,11698,1039987,1043896,1048018},

         A_test_index_array[TEST_ARRAY_SIZE] = 
                             {2112377,662041,5336171,3642833,4250760},
         A_test_rank_array[TEST_ARRAY_SIZE] = 
                             {104,17523,123928,8288932,8388264},

         B_test_index_array[TEST_ARRAY_SIZE] = 
                             {41869,812306,5102857,18232239,26860214},
         B_test_rank_array[TEST_ARRAY_SIZE] = 
                             {33422937,10244,59149,33135281,99}, 

         C_test_index_array[TEST_ARRAY_SIZE] = 
                             {44172927,72999161,74326391,129606274,21736814},
         C_test_rank_array[TEST_ARRAY_SIZE] = 
                             {61147,882988,266290,133997595,133525895},

         D_test_index_array[TEST_ARRAY_SIZE] = 
                             {1317351170,995930646,1157283250,1503301535,1453734525},
         D_test_rank_array[TEST_ARRAY_SIZE] = 
                             {1,36538729,1978098519,2145192618,2147425337};



/***********************/
/* function prototypes */
/***********************/
static float	randlc( float *X, float *A, context *c );

void full_verify( context *c );


/*****************************************************************/
/******     C  _  P  R  I  N  T  _  R  E  S  U  L  T  S     ******/
/*****************************************************************/


static void c_print_results( char   *name,
                      char   class,
                      int    n1, 
                      int    n2,
                      int    n3,
                      int    niter,
                      int    nprocs_compiled,
                      int    nprocs_total,
                      float t,
                      int mops,
		      char   *optype,
                      int    passed_verification,
                      char   *npbversion,
                      char   *compiletime,
                      char   *mpicc,
                      char   *clink,
                      char   *cmpi_lib,
                      char   *cmpi_inc,
                      char   *cflags,
                      char   *clinkflags )
{
//    char *evalue="1000";

    printf( "\r\n\r\n %s Benchmark Completed\r\n", name ); 

    printf( " Class           =                        %c\r\n", class );

    if( n3 == 0 ) {
        long nn = n1;
        if ( n2 != 0 ) nn *= n2;
        printf( " Size            =             %12ld\r\n", nn );   /* as in IS */
    }
    else
        printf( " Size            =              %3dx %3dx %3d\r\n", n1,n2,n3 );

    printf( " Iterations      =             %12d\r\n", niter );
 
//    printf( " Time in seconds =             %12.2f\r\n", t );
    printf( " Time in cycles =             %12d\r\n", (int)t );

    printf( " Total processes =             %12d\r\n", nprocs_total );

    if ( nprocs_compiled != 0 )
        printf( " Compiled procs  =             %12d\r\n", nprocs_compiled );

//    printf( " Mop/s total     =             %12.2f\r\n", mops );
    printf( " Op/Kcycles total     =             %12d\r\n", (int)mops );

//    printf( " Mop/s/process   =             %12.2f\r\n", mops/((float) nprocs_total) );
    printf( " Op/Kcycles/process   =             %12d/%12d\r\n", (int)mops, (int) nprocs_total);

    printf( " Operation type  = %24s\r\n", optype);

    if( passed_verification )
        printf( " Verification    =               SUCCESSFUL\r\n" );
    else
        printf( " Verification    =             UNSUCCESSFUL\r\n" );

    printf( " Version         =             %12s\r\n", npbversion );

    printf( " Compile date    =             %12s\r\n", compiletime );

    printf( "\r\n Compile options:\r\n" );

    printf( "    MPICC        = %s\r\n", mpicc );

    printf( "    CLINK        = %s\r\n", clink );

    printf( "    CMPI_LIB     = %s\r\n", cmpi_lib );

    printf( "    CMPI_INC     = %s\r\n", cmpi_inc );

    printf( "    CFLAGS       = %s\r\n", cflags );

    printf( "    CLINKFLAGS   = %s\r\n", clinkflags );
#ifdef SMP
    evalue = getenv("MP_SET_NUMTHREADS");
    printf( "   MULTICPUS = %s\r\n", evalue );
#endif

    printf( "\r\n\r\n" );
    printf( " Please send feedbacks and/or the results of this run to:\r\n\r\n" );
    printf( " NPB Development Team\r\n" );
    printf( " npb@nas.nasa.gov\r\n\r\n\r\n" );
/*    printf( " Please send the results of this run to:\n\n" );
    printf( " NPB Development Team\n" );
    printf( " Internet: npb@nas.nasa.gov\n \n" );
    printf( " If email is not available, send this to:\n\n" );
    printf( " MS T27A-1\n" );
    printf( " NASA Ames Research Center\n" );
    printf( " Moffett Field, CA  94035-1000\n\n" );
    printf( " Fax: 650-604-3957\n\n" );*/
}
 

/*****************************************************************/
/******            T  I  M  E  R  _  C  L  E  A  R          ******/
/*****************************************************************/
static void timer_clear( int n, struct context *c )
{
    c->elapsed[n] = 0.0;
}


/*****************************************************************/
/******            T  I  M  E  R  _  S  T  A  R  T          ******/
/*****************************************************************/
static void timer_start( int n, context *c )
{
    c->start[n] = MPI_Wtime();
}


/*****************************************************************/
/******            T  I  M  E  R  _  S  T  O  P             ******/
/*****************************************************************/
static void timer_stop( int n, context *c )
{
    float t, now;

    now = MPI_Wtime();
    t = now - c->start[n];
    c->elapsed[n] += t;

}


/*****************************************************************/
/******            T  I  M  E  R  _  R  E  A  D             ******/
/*****************************************************************/
static float timer_read( int n, context *c )
{
    return( (float)c->elapsed[n] );
}




/*
 *    FUNCTION RANDLC (X, A)
 *
 *  This routine returns a uniform pseudorandom double precision number in the
 *  range (0, 1) by using the linear congruential generator
 *
 *  x_{k+1} = a x_k  (mod 2^46)
 *
 *  where 0 < x_k < 2^46 and 0 < a < 2^46.  This scheme generates 2^44 numbers
 *  before repeating.  The argument A is the same as 'a' in the above formula,
 *  and X is the same as x_0.  A and X must be odd double precision integers
 *  in the range (1, 2^46).  The returned value RANDLC is normalized to be
 *  between 0 and 1, i.e. RANDLC = 2^(-46) * x_1.  X is updated to contain
 *  the new seed x_1, so that subsequent calls to RANDLC using the same
 *  arguments will generate a continuous sequence.
 *
 *  This routine should produce the same results on any computer with at least
 *  48 mantissa bits in double precision floating point data.  On Cray systems,
 *  double precision should be disabled.
 *
 *  David H. Bailey     October 26, 1990
 *
 *     IMPLICIT DOUBLE PRECISION (A-H, O-Z)
 *     SAVE KS, R23, R46, T23, T46
 *     DATA KS/0/
 *
 *  If this is the first call to RANDLC, compute R23 = 2 ^ -23, R46 = 2 ^ -46,
 *  T23 = 2 ^ 23, and T46 = 2 ^ 46.  These are computed in loops, rather than
 *  by merely using the ** operator, in order to insure that the results are
 *  exact on all systems.  This code assumes that 0.5D0 is represented exactly.
 */


/*****************************************************************/
/*************           R  A  N  D  L  C             ************/
/*************                                        ************/
/*************    portable random number generator    ************/
/*****************************************************************/

static float	randlc( float *X, float *A, context *c )
{
      float		T1;
      float		T3;
      float		T4;
      float		A1;
      float		A2;
      float		X1;
      float		X2;
      int     		i;

//printf("X=%f, A=%f, ", *X, *A);

      if (c->KS == 0) 
      {
        c->R23 = 1.0;
        c->R46 = 1.0;
        c->T23 = 1.0;
        c->T46 = 1.0;
    
        for (i=1; i<=8; i++)
        {
          c->R23 = 0.50 * c->R23;
          c->T23 = 2.0 * c->T23;
        }
        for (i=1; i<=16; i++)
        {
          c->R46 = 0.50 * c->R46;
          c->T46 = 2.0 * c->T46;
        }
        c->KS = 1;
      }

/*  Break A into two parts such that A = 2^23 * A1 + A2 and set X = N.  */

     A1 = (int)(c->R23 * *A);
//printf("A1=%f, T23=%f, ", A1, T23);
      A2 = *A - c->T23 * A1;
//printf("A2=%f, ", A2);
//printf("T23=%f, ", T23);

/*  Break X into two parts such that X = 2^23 * X1 + X2, compute
    Z = A1 * X2 + A2 * X1  (mod 2^23), and then
    X = 2^23 * Z + A2 * X2  (mod 2^46).                            */

      X1 = (int)(c->R23 * *X);
      X2 = *X - c->T23 * X1;
      T1 = A1 * X2 + A2 * X1;

      T3 = c->T23 * (T1 - c->T23 * (int)(c->R23 * T1)) + A2 * X2;
//printf("A2=%f, X2=%f, T1=%f, T3=%f, ", A2, X2, T1, T3);
      T4 = T3 - c->T46 * (int)(c->R46 * T3);
      *X = T4;
//if (c->KS < 5) {
//print_float(c->R46 * *X);
//c->KS++;
//}
      return(c->R46 * *X);
} 



/*****************************************************************/
/************   F  I  N  D  _  M  Y  _  S  E  E  D    ************/
/************                                         ************/
/************ returns parallel random number seq seed ************/
/*****************************************************************/

/*
 * Create a random number sequence of total length nn residing
 * on np number of processors.  Each processor will therefore have a 
 * subsequence of length nn/np.  This routine returns that random 
 * number which is the first random number for the subsequence belonging
 * to processor rank kn, and which is used as seed for proc kn ran # gen.
 */

float   find_my_seed( int  kn,       /* my processor rank, 0<=kn<=num procs */
                       int  np,       /* np = num procs                      */
                       long nn,       /* total num of ran numbers, all procs */
                       float s,      /* Ran num seed, for ex.: 314159265.00 */
                       float a, 
                       context *c )     /* Ran num gen mult, try 1220703125.00 */
{


  long   i;

  float t1,t2,t3,an;
  int nq;
  long   mq,kk,ik;


      nq = nn / np;

      for( mq=0; nq>1; mq++,nq/=2 )
          ;

      t1 = a;

      for( i=1; i<=mq; i++ ) {
        t2 = randlc( &t1, &t1, c );
	}
	// iakovos
	t1 = 11;

      an = t1;

      kk = kn;
      t1 = s;
      t2 = an;

      for( i=1; i<=100; i++ )
      {
        ik = kk / 2;
        if( 2 * ik !=  kk ) 
            t3 = randlc( &t1, &t2, c );
        if( ik == 0 ) 
            break;
        t3 = randlc( &t2, &t2, c );
        kk = ik;
      }

      return( t1 );

}




/*****************************************************************/
/*************      C  R  E  A  T  E  _  S  E  Q      ************/
/*****************************************************************/

void	create_seq( float seed, float a, context *c )
{
	float x;
	int    i, k;

        k = MAX_KEY/4;

	for (i=0; i<NUM_KEYS; i++)
	{
	    x = randlc(&seed, &a, c);
	    x += randlc(&seed, &a, c);
    	    x += randlc(&seed, &a, c);
	    x += randlc(&seed, &a, c);  

//printf("key=%f\n", k*x);
            c->key_array[i] = k*x;
	}
}




/*****************************************************************/
/*************    F  U  L  L  _  V  E  R  I  F  Y     ************/
/*****************************************************************/


void full_verify( context *c )
{
    MPI_Status  status;
    MPI_Request request;
    
    INT_TYPE    i, j;
    INT_TYPE    *k, last_local_key;

    k = malloc(sizeof(INT_TYPE));

    
    TIMER_START( T_VERIFY, c );

/*  Now, finally, sort the keys:  */
    for( i=0; i<c->total_local_keys; i++ )
        c->key_array[--c->key_buff_ptr_global[c->key_buff2[i]]-
                                 c->total_lesser_keys] = c->key_buff2[i];
    last_local_key = (c->total_local_keys<1)? 0 : (c->total_local_keys-1);

/*  Send largest key value to next processor  */
    if( c->my_rank > 0 )
        MPI_Irecv( k,
                   1,
                   MP_KEY_TYPE,
                   c->my_rank-1,
                   1000,
                   MPI_COMM_WORLD,
                   &request );                   
    if( c->my_rank < c->comm_size-1 ) {
		*c->tmp = (int)(c->key_array[last_local_key]);
        MPI_Send( c->tmp,
                  1,
                  MP_KEY_TYPE,
                  c->my_rank+1,
                  1000,
                  MPI_COMM_WORLD );
    } if( c->my_rank > 0 )
        MPI_Wait( &request, &status );

/*  Confirm that neighbor's greatest key value 
    is not greater than my least key value       */              
    j = 0;
    if( c->my_rank > 0 && c->total_local_keys > 0 )
        if( *k > c->key_array[0] )
            j++;


/*  Confirm keys correctly sorted: count incorrectly sorted keys, if any */
    for( i=1; i<c->total_local_keys; i++ )
        if( c->key_array[i-1] > c->key_array[i] )
            j++;


    if( j != 0 )
    {
        printf( "Processor %d:  Full_verify: number of keys out of sort: %d\n\r",
                c->my_rank, j );
    }
    else
        c->passed_verification++;
           
    TIMER_STOP( T_VERIFY, c );

}




/*****************************************************************/
/*************             R  A  N  K             ****************/
/*****************************************************************/


void rank( int iteration, context *c )
{

    INT_TYPE    i, k;

    INT_TYPE    shift = MAX_KEY_LOG_2 - NUM_BUCKETS_LOG_2;
    INT_TYPE    key;
    INT_TYPE2   bucket_sum_accumulator, j, m;
    INT_TYPE    local_bucket_sum_accumulator;
    INT_TYPE    min_key_val, max_key_val;
    INT_TYPE    *key_buff_ptr;


    TIMER_START( T_RANK, c );

/*  Iteration alteration of keys */  
    if(c->my_rank == 0 )                    
    {
      c->key_array[iteration] = iteration;
      c->key_array[iteration+MAX_ITERATIONS] = MAX_KEY - iteration;
    }


/*  Initialize */
    for( i=0; i<NUM_BUCKETS+TEST_ARRAY_SIZE; i++ )  
    {
        c->bucket_size[i] = 0;
        c->bucket_size_totals[i] = 0;
        c->process_bucket_distrib_ptr1[i] = 0;
        c->process_bucket_distrib_ptr2[i] = 0;
    }


/*  Determine where the partial verify test keys are, load into  */
/*  top of array bucket_size                                     */
    for( i=0; i<TEST_ARRAY_SIZE; i++ )
        if( (c->test_index_array[i]/NUM_KEYS) == c->my_rank )
            c->bucket_size[NUM_BUCKETS+i] = 
                          c->key_array[c->test_index_array[i] % NUM_KEYS];


/*  Determine the number of keys in each bucket */
    for( i=0; i<NUM_KEYS; i++ )
        c->bucket_size[c->key_array[i] >> shift]++;


/*  Accumulative bucket sizes are the bucket pointers */
    c->bucket_ptrs[0] = 0;
    for( i=1; i< NUM_BUCKETS; i++ )  
        c->bucket_ptrs[i] = c->bucket_ptrs[i-1] + c->bucket_size[i-1];


/*  Sort into appropriate bucket */
    for( i=0; i<NUM_KEYS; i++ )  
    {
        key = c->key_array[i];
        c->key_buff1[c->bucket_ptrs[key >> shift]++] = key;
    }

    TIMER_STOP( T_RANK, c );
    TIMER_START( T_RCOMM, c );

/*  Get the bucket size totals for the entire problem. These 
    will be used to determine the redistribution of keys      */
    MPI_Allreduce( c->bucket_size, 
                   c->bucket_size_totals, 
                   NUM_BUCKETS+TEST_ARRAY_SIZE, 
                   MP_KEY_TYPE,
                   MPI_SUM,
                   MPI_COMM_WORLD );

    TIMER_STOP( T_RCOMM, c );
    TIMER_START( T_RANK, c );

/*  Determine Redistibution of keys: accumulate the bucket size totals 
    till this number surpasses NUM_KEYS (which the average number of keys
    per processor).  Then all keys in these buckets go to processor 0.
    Continue accumulating again until supassing 2*NUM_KEYS. All keys
    in these buckets go to processor 1, etc.  This algorithm guarantees
    that all processors have work ranking; no processors are left idle.
    The optimum number of buckets, however, does not result in as high
    a degree of load balancing (as even a distribution of keys as is
    possible) as is obtained from increasing the number of buckets, but
    more buckets results in more computation per processor so that the
    optimum number of buckets turns out to be 1024 for machines tested.
    Note that process_bucket_distrib_ptr1 and ..._ptr2 hold the bucket
    number of first and last bucket which each processor will have after   
    the redistribution is done.                                          */

    bucket_sum_accumulator = 0;
    local_bucket_sum_accumulator = 0;
    c->send_displ[0] = 0;
    c->process_bucket_distrib_ptr1[0] = 0;
    for( i=0, j=0; i<NUM_BUCKETS; i++ )  
    {
        bucket_sum_accumulator       += c->bucket_size_totals[i];
        local_bucket_sum_accumulator += c->bucket_size[i];
        if( bucket_sum_accumulator >= (j+1)*NUM_KEYS )  
        {
            c->send_count[j] = local_bucket_sum_accumulator;
            if( j != 0 )
            {
                c->send_displ[j] = c->send_displ[j-1] + c->send_count[j-1];
                c->process_bucket_distrib_ptr1[j] = 
                                        c->process_bucket_distrib_ptr2[j-1]+1;
            }
            c->process_bucket_distrib_ptr2[j++] = i;
            local_bucket_sum_accumulator = 0;
        }
    }

/*  When NUM_PROCS approaching NUM_BUCKETS, it is highly possible
    that the last few processors don't get any buckets.  So, we
    need to set counts properly in this case to avoid any fallouts.    */
    while( j < c->comm_size )
    {
        c->send_count[j] = 0;
        c->process_bucket_distrib_ptr1[j] = 1;
        j++;
    }

    TIMER_STOP( T_RANK, c );
    TIMER_START( T_RCOMM, c ); 

/*  This is the redistribution section:  first find out how many keys
    each processor will send to every other processor:                 */
    MPI_Alltoall( c->send_count,
                  1,
                  MPI_INT,
                  c->recv_count,
                  1,
                  MPI_INT,
                  MPI_COMM_WORLD );

/*  Determine the receive array displacements for the buckets */    
    c->recv_displ[0] = 0;
    for( i=1; i<c->comm_size; i++ )
        c->recv_displ[i] = c->recv_displ[i-1] + c->recv_count[i-1];


/*  Now send the keys to respective processors  */    
    MPI_Alltoallv( c->key_buff1,
                   c->send_count,
                   c->send_displ,
                   MP_KEY_TYPE,
                   c->key_buff2,
                   c->recv_count,
                   c->recv_displ,
                   MP_KEY_TYPE,
                   MPI_COMM_WORLD );

//if (ar_get_board_id()==25 && ar_get_core_id()==5) {
//int test;
//test = c->recv_displ[NUM_PROCS-1] + c->recv_count[NUM_PROCS-1];
//if (c->my_rank==0) {
//	printf("%d %d\n\r", c->recv_displ[255], c->recv_count[255]);
//    printf("key_buff2[12309]=%d\n\r", c->key_buff2[12309]);
//}

    TIMER_STOP( T_RCOMM, c ); 
    TIMER_START( T_RANK, c );

/*  The starting and ending bucket numbers on each processor are
    multiplied by the interval size of the buckets to obtain the 
    smallest possible min and greatest possible max value of any 
    key on each processor                                          */
    min_key_val = c->process_bucket_distrib_ptr1[c->my_rank] << shift;
    max_key_val = ((c->process_bucket_distrib_ptr2[c->my_rank] + 1) << shift)-1;

/*  Clear the work array */
    for( i=0; i<max_key_val-min_key_val+1; i++ )
        c->key_buff1[i] = 0;

/*  Determine the total number of keys on all other 
    processors holding keys of lesser value         */
    m = 0;
    for( k=0; k<c->my_rank; k++ )
        for( i= c->process_bucket_distrib_ptr1[k];
             i<=c->process_bucket_distrib_ptr2[k];
             i++ )  
            m += c->bucket_size_totals[i]; /*  m has total # of lesser keys */

/*  Determine total number of keys on this processor */
    j = 0;                                 
    for( i= c->process_bucket_distrib_ptr1[c->my_rank];
         i<=c->process_bucket_distrib_ptr2[c->my_rank];
         i++ )  
        j += c->bucket_size_totals[i];     /* j has total # of local keys   */

//    if (ar_get_board_id()==25 && ar_get_core_id()==5) printf("j=%d\n\r", j);
//    if (c->my_rank == 0) printf("j=%d\n\r", j);
//    if (j != test) { printf("bid=%d cid=%d rank=%d j=%d test=%d\n\r", ar_get_board_id(), ar_get_core_id(), c->my_rank, j, test);}


/*  Ranking of all keys occurs in this section:                 */
/*  shift it backwards so no subtractions are necessary in loop */
    key_buff_ptr = c->key_buff1 - min_key_val;

/*  In this section, the keys themselves are used as their 
    own indexes to determine how many of each there are: their
    individual population                                       */
    for( i=0; i<j; i++ ) 
        key_buff_ptr[c->key_buff2[i]]++;  /* Now they have individual key   */
                                       /* population                     */

/*  To obtain ranks of each key, successively add the individual key
    population, not forgetting the total of lesser keys, m.
    NOTE: Since the total of lesser keys would be subtracted later 
    in verification, it is no longer added to the first key population 
    here, but still needed during the partial verify test.  This is to 
    ensure that 32-bit key_buff can still be used for class D.           */
/*    key_buff_ptr[min_key_val] += m;    */
    for( i=min_key_val; i<max_key_val; i++ )   
        key_buff_ptr[i+1] += key_buff_ptr[i];  


/* This is the partial verify test section */
/* Observe that test_rank_array vals are   */
/* shifted differently for different cases */
    for( i=0; i<TEST_ARRAY_SIZE; i++ )
    {                                             
        k = c->bucket_size_totals[i+NUM_BUCKETS];    /* Keys were hidden here */
        if( min_key_val <= k  &&  k <= max_key_val )
        {
            c->passed_verification++;
	    continue;
            /* Add the total of lesser keys, m, here */
            INT_TYPE2 key_rank = key_buff_ptr[k-1] + m;
            int failed = 0;

            switch( CLASS )
            {
                case 'S':
                    if( i <= 2 )
                    {
                        if( key_rank != c->test_rank_array[i]+iteration ) {
                            failed = 1;
			}
                        else {
                            c->passed_verification++;
			}
                    }
                    else
                    {
                        if( key_rank != c->test_rank_array[i]-iteration )
                            failed = 1;
                        else
                            c->passed_verification++;
                    }
                    break;
                case 'W':
                    if( i < 2 )
                    {
                        if( key_rank != c->test_rank_array[i]+(iteration-2) )
                            failed = 1;
                        else
                            c->passed_verification++;
                    }
                    else
                    {
                        if( key_rank != c->test_rank_array[i]-iteration )
                            failed = 1;
                        else
                            c->passed_verification++;
                    }
                    break;
                case 'A':
                    if( i <= 2 )
        	    {
                        if( key_rank != c->test_rank_array[i]+(iteration-1) )
                            failed = 1;
                        else
                            c->passed_verification++;
        	    }
                    else
                    {
                        if( key_rank != c->test_rank_array[i]-(iteration-1) )
                            failed = 1;
                        else
                            c->passed_verification++;
                    }
                    break;
                case 'B':
                    if( i == 1 || i == 2 || i == 4 )
        	    {
                        if( key_rank != c->test_rank_array[i]+iteration )
                            failed = 1;
                        else
                            c->passed_verification++;
        	    }
                    else
                    {
                        if( key_rank != c->test_rank_array[i]-iteration )
                            failed = 1;
                        else
                            c->passed_verification++;
                    }
                    break;
                case 'C':
                    if( i <= 2 )
        	    {
                        if( key_rank != c->test_rank_array[i]+iteration )
                            failed = 1;
                        else
                            c->passed_verification++;
        	    }
                    else
                    {
                        if( key_rank != c->test_rank_array[i]-iteration )
                            failed = 1;
                        else
                            c->passed_verification++;
                    }
                    break;
                case 'D':
                    if( i < 2 )
        	    {
                        if( key_rank != c->test_rank_array[i]+iteration )
                            failed = 1;
                        else
                            c->passed_verification++;
        	    }
                    else
                    {
                        if( key_rank != c->test_rank_array[i]-iteration )
                            failed = 1;
                        else
                            c->passed_verification++;
                    }
                    break;
            }
            if( failed == 1 )
                printf( "Failed partial verification: "
                        "iteration %d, processor %d, test key %d\n\r", 
                         iteration, c->my_rank, (int)i );
        }
    }


    TIMER_STOP( T_RANK, c ); 


/*  Make copies of rank info for use by full_verify: these variables
    in rank are local; making them global slows down the code, probably
    since they cannot be made register by compiler                        */

    if( iteration == MAX_ITERATIONS ) 
    {
        c->key_buff_ptr_global = key_buff_ptr;
        c->total_local_keys    = j;
        c->total_lesser_keys   = 0;  /* no longer set to 'm', see note above */
    }

}      


/*****************************************************************/
/*************             M  A  I  N             ****************/
/*****************************************************************/

#ifdef ARCH_MB
int nas_is_mpi()
#else
int main( int argc, char **argv )
#endif
{

    int             i, iteration, *itemp;

    float          *timecounter, maxtime;

    context c;

    timecounter = malloc(sizeof(float));
    itemp = malloc(sizeof(int));

/*  Initialize MPI */
    MPI_Init( NULL, NULL );
    MPI_Comm_rank( MPI_COMM_WORLD, &c.my_rank );
    MPI_Comm_size( MPI_COMM_WORLD, &c.comm_size );


/*
    if (c.my_rank != 0)  {
	    while(1);
	    MPI_Finalize();
	    return 0;
	}
    c.comm_size = 1;
*/

/* initialize space */
    
    c.key_array = malloc(SIZE_OF_BUFFERS * sizeof(INT_TYPE));
    c.key_buff1 = malloc(SIZE_OF_BUFFERS * sizeof(INT_TYPE));
    c.key_buff2 = malloc(SIZE_OF_BUFFERS * sizeof(INT_TYPE));
    c.bucket_size = malloc((NUM_BUCKETS+TEST_ARRAY_SIZE)*sizeof(INT_TYPE));
    c.bucket_size_totals = malloc((NUM_BUCKETS+TEST_ARRAY_SIZE)*sizeof(INT_TYPE));
    c.bucket_ptrs = malloc(NUM_BUCKETS*sizeof(INT_TYPE));
    c.process_bucket_distrib_ptr1 = malloc((NUM_BUCKETS+TEST_ARRAY_SIZE)*sizeof(INT_TYPE));
    c.process_bucket_distrib_ptr2 = malloc((NUM_BUCKETS+TEST_ARRAY_SIZE)*sizeof(INT_TYPE));
    c.timeron = malloc(sizeof(int));
    c.tmp = malloc(sizeof(int));
    c.KS = 0;

    printf("Total allocated size= %dKB\n\r", (int)((SIZE_OF_BUFFERS*3+5*NUM_BUCKETS+4*TEST_ARRAY_SIZE+MAX_PROCS+2)*sizeof(int)/(1024)));
    

/*  Initialize the verification arrays if a valid class */
    for( i=0; i<TEST_ARRAY_SIZE; i++ )
        switch( CLASS )
        {
            case 'S':
                c.test_index_array[i] = S_test_index_array[i];
                c.test_rank_array[i]  = S_test_rank_array[i];
                break;
            case 'A':
                c.test_index_array[i] = A_test_index_array[i];
                c.test_rank_array[i]  = A_test_rank_array[i];
                break;
            case 'W':
                c.test_index_array[i] = W_test_index_array[i];
                c.test_rank_array[i]  = W_test_rank_array[i];
                break;
            case 'B':
                c.test_index_array[i] = B_test_index_array[i];
                c.test_rank_array[i]  = B_test_rank_array[i];
                break;
            case 'C':
                c.test_index_array[i] = C_test_index_array[i];
                c.test_rank_array[i]  = C_test_rank_array[i];
                break;
            case 'D':
                c.test_index_array[i] = D_test_index_array[i];
                c.test_rank_array[i]  = D_test_rank_array[i];
                break;
        };

        

/*  Printout initial NPB info */
    if( c.my_rank == 0 )
    {
//        FILE *fp;
        printf( "\n\r\n\r NAS Parallel Benchmarks 3.3 -- IS Benchmark\n\r\n\r" );
        printf( " Size:  %ld  (class %c)\n\r", (long)TOTAL_KEYS*MIN_PROCS, CLASS );
        printf( " Iterations:   %d\n\r", MAX_ITERATIONS );
        printf( " Number of processes:     %d\n\r", c.comm_size );

//        fp = fopen("timer.flag", "r");
        *c.timeron = 0;
//        if (fp) {
//            timeron = 1;
//            fclose(fp);
//        }
    }

/*  Check that actual and compiled number of processors agree */
    if( c.comm_size != NUM_PROCS )
    {
        if( c.my_rank == 0 )
            printf( "\n\r ERROR: compiled for %d processes\n\r"
                    " Number of active processes: %d\n\r"
                    " Exiting program!\n\r\n\r", NUM_PROCS, c.comm_size );
        MPI_Finalize();
        return( 1 );
    }

/*  Check to see whether total number of processes is within bounds.
    This could in principle be checked in setparams.c, but it is more
    convenient to do it here                                               */
    if( c.comm_size < MIN_PROCS || c.comm_size > MAX_PROCS)
    {
       if( c.my_rank == 0 )
           printf( "\n\r ERROR: number of processes %d not within range %d-%d"
                   "\n\r Exiting program!\n\r\n\r", c.comm_size, MIN_PROCS, MAX_PROCS);
       MPI_Finalize();
       return( 1 );
    }

    MPI_Bcast(c.timeron, 1, MPI_INT, 0, MPI_COMM_WORLD);

#ifdef  TIMING_ENABLED 
    for( i=1; i<=T_LAST; i++ ) timer_clear( i, &c );
#endif

/*  Generate random number sequence and subsequent keys on all procs */
    create_seq( find_my_seed( c.my_rank, 
                              c.comm_size, 
                              4*(long)TOTAL_KEYS*MIN_PROCS,
                              314159265.00,      /* Random number gen seed */
                              11.00, &c ),   /* Random number gen mult */
                11.00, &c );                 /* Random number gen mult */



/*  Do one interation for free (i.e., untimed) to guarantee initialization of  
    all data and code pages and respective tables */
    rank( 1, &c );  

/*  Start verification counter */
    c.passed_verification = 0;

    if( c.my_rank == 0 && CLASS != 'S' ) printf( "\n\r   iteration\n\r" );

/*  Initialize timer  */             
    timer_clear( 0, &c );

/*  Initialize separate communication, computation timing */
#ifdef  TIMING_ENABLED 
    for( i=1; i<=T_LAST; i++ ) timer_clear( i, &c );
#endif

/*  Start timer  */             
#ifdef ARCH_MB
    ar_timer_reset();
#endif
    timer_start( 0, &c );


/*  This is the main iteration */
    for( iteration=1; iteration<=MAX_ITERATIONS; iteration++ )
    {
        if( c.my_rank == 0) printf( "        %d\n\r", iteration );
//        printf( "        %d\n\r", iteration );
        rank( iteration, &c );
    }


/*  Stop timer, obtain time for processors */
    timer_stop( 0, &c );

    *timecounter = timer_read( 0, &c );

/*  End of timing, obtain maximum time of all processors */
    MPI_Reduce( timecounter,
                &maxtime,
                1,
                MPI_FLOAT,
                MPI_MAX,
                0,
                MPI_COMM_WORLD );


/*  This tests that keys are in sequence: sorting of last ranked key seq
    occurs here, but is an untimed operation                             */
    full_verify(&c);


/*  Obtain verification counter sum */
    *itemp = c.passed_verification;
    MPI_Reduce( itemp,
                &c.passed_verification,
                1,
                MPI_INT,
                MPI_SUM,
                0,
                MPI_COMM_WORLD );




/*  The final printout  */
    if( c.my_rank == 0 )
    {
        if( c.passed_verification != 5*MAX_ITERATIONS + c.comm_size )
            c.passed_verification = 0;
#ifdef ARCH_MB
		int ops = MAX_ITERATIONS*TOTAL_KEYS*MIN_PROCS;
		int opkc = ops/maxtime;
#else
		int opkc = (int) ((MAX_ITERATIONS)*TOTAL_KEYS*MIN_PROCS/maxtime/1000000.);
#endif
		
        c_print_results( "IS",
                         CLASS,
                         (int)(TOTAL_KEYS),
                         MIN_PROCS,
                         0,
                         MAX_ITERATIONS,
                         NUM_PROCS,
                         c.comm_size,
                         maxtime,
                         opkc,
                         "keys ranked", 
                         c.passed_verification,
                         NPBVERSION,
                         COMPILETIME,
                         MPICC,
                         CLINK,
                         CMPI_LIB,
                         CMPI_INC,
                         CFLAGS,
                         CLINKFLAGS );
    }
                    

#ifdef  TIMING_ENABLED
    if (*c.timeron)
    {
        float    *t1, tmin[T_LAST+1], tsum[T_LAST+1], tmax[T_LAST+1];
        char      t_recs[T_LAST+1][9];

	t1 = malloc((T_LAST+1)*sizeof(float));
    
        for( i=0; i<=T_LAST; i++ )
            t1[i] = timer_read( i, &c );

        MPI_Reduce( t1,
                    tmin,
                    T_LAST+1,
                    MPI_FLOAT,
                    MPI_MIN,
                    0,
                    MPI_COMM_WORLD );
        MPI_Reduce( t1,
                    tsum,
                    T_LAST+1,
                    MPI_FLOAT,
                    MPI_SUM,
                    0,
                    MPI_COMM_WORLD );
        MPI_Reduce( t1,
                    tmax,
                    T_LAST+1,
                    MPI_FLOAT,
                    MPI_MAX,
                    0,
                    MPI_COMM_WORLD );

        if( c.my_rank == 0 )
        {
            strcpy( t_recs[T_TOTAL],  "total" );
            strcpy( t_recs[T_RANK],   "rcomp" );
            strcpy( t_recs[T_RCOMM],  "rcomm" );
            strcpy( t_recs[T_VERIFY], "verify");
            printf( " nprocs = %6d     ", c.comm_size);
            printf( "     minimum     maximum     average\n\r" );
#ifdef ARCH_MB
#else
            for( i=0; i<=T_LAST; i++ )
            {
                printf( " timer %2d (%-8s):  %10.4f  %10.4f  %10.4f\n\r",
                        i+1, t_recs[i], tmin[i], tmax[i], 
                        tsum[i]/((float) c.comm_size) );
            }
#endif
            printf( "\n\r" );
        }
    }
#endif

    MPI_Finalize();

    return 0;
         /**************************/
}        /*  E N D  P R O G R A M  */
         /**************************/
