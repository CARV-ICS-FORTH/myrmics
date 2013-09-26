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
#define free   kt_free
#define sqrtf  ar_float_sqrt

#define MPI_Wtime FMPI_Wtime
  
#else
#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#endif

//npbparams.h
#define CLASS 'S'
#define NPROCS 16

#if CLASS == 'S'
#define NA 1400
#define NONZER 7
#define SHIFT 10.
#endif

#if CLASS == 'W'
#define NA 7000
#define NONZER 8
#define SHIFT 12.
#endif

#if CLASS == 'A'
#define NA 14000
#define NONZER 11
#define SHIFT 20.
#endif

#if NPROCS == 1
#define NUM_PROC_COLS 1
#define NUM_PROC_ROWS 1
#endif

#if NPROCS == 2
#define NUM_PROC_COLS 2
#define NUM_PROC_ROWS 1
#endif

#if NPROCS == 4
#define NUM_PROC_COLS 2
#define NUM_PROC_ROWS 2
#endif

#if NPROCS == 8
#define NUM_PROC_COLS 4
#define NUM_PROC_ROWS 2
#endif

#if NPROCS == 16
#define NUM_PROC_COLS 4
#define NUM_PROC_ROWS 4
#endif

#if NPROCS == 32
#define NUM_PROC_COLS 8
#define NUM_PROC_ROWS 4
#endif

#if NPROCS == 64
#define NUM_PROC_COLS 8
#define NUM_PROC_ROWS 8
#endif

#if NPROCS == 128
#define NUM_PROC_COLS 16
#define NUM_PROC_ROWS 8
#endif

#if NPROCS == 256
#define NUM_PROC_COLS 16
#define NUM_PROC_ROWS 16
#endif

#if NPROCS == 512
#define NUM_PROC_COLS 32
#define NUM_PROC_ROWS 16
#endif

#define TMP_SIZE ((NA/1400)*350*3+50)

#define NITER 15
#define RCOND 0.1
///

#define NUM_PROCS NUM_PROC_COLS * NUM_PROC_ROWS
#define NZ NA*(NONZER+1)/NUM_PROCS*(NONZER+1)+NONZER+NA*(NONZER+2+NUM_PROCS/256)/NUM_PROC_COLS
//#define NZ 50000 // NPROCS=2 CLASS=S
//#define NZ 100000 // NPROCS=1 CLASS=S
//#define NZ  700000 // NPROCS=1 CLASS=W
//#define NZ 1050000 // NPROCS=1 CLASS=A


// f2c
#define TRUE_ (1)
#define FALSE_ (0)
typedef int ftnlen;
#define max(a,b) ((a) >= (b) ? (a) : (b))
//
#define abs(a) ((a) >= 0 ? (a) : (-a))

float  randlc_( float *X, float *A );

void full_verify( void );

int print_results__(char *name__, char *class__, int *n1,
         int *n2, int *n3, int *niter, int *nprocs_compiled__,
         int *nprocs_total__, float *t, float *mops, char *
        optype, int *verified, char *npbversion, char *compiletime, char *
        cs1, char *cs2, char *cs3, char *cs4, char *cs5, char *cs6, char *cs7,
         ftnlen name_len, ftnlen class_len, ftnlen optype_len, ftnlen
        npbversion_len, ftnlen compiletime_len, ftnlen cs1_len, ftnlen
        cs2_len, ftnlen cs3_len, ftnlen cs4_len, ftnlen cs5_len, ftnlen
        cs6_len, ftnlen cs7_len);

#ifndef ARCH_MB


#include <math.h>

float ar_float_pow (float x, float y)
  {
    float result = 1;
    int result_int;
    int inverse;
    int i;
  
    if (y < 0) {
      inverse = 1;
      y = -y;
    }
  
    for (i = 0; i < 16; i++) {
      x = sqrtf(x);
      y *= 2;
    }
  
    x = sqrtf(x);
    while (y > 0.5) {
      result *= x;
      y -= 0.5;
    }
  
    if (inverse) result = 1.0 / result;
  
    return result;

  }

#endif

typedef struct context {
 struct {
    float *v, *aelt, *a, *x, *z__, *p, *q, *r__, *w;
 } main_flt_mem__;

 struct {
    int *colidx, *rowstr, *iv, *arow, *acol;
 } main_int_mem__;

 struct {
    int naa, nzz, npcols, nprows, proc_col__, proc_row__, firstrow, 
	    lastrow, firstcol, lastcol, exch_proc__, exch_recv_length__, 
	    send_start__, send_len__;
 } partit_size__;

 struct {
    MPI_Datatype dp_type__;
    int me, nprocs, root, dummy_dp_type__;
 } mpistuff_;

 struct {
    int timeron;
 } timers_;

 struct {
    float amult, tran;
 } urando_;

 unsigned int start[4], elapsed[4];

} context;

/* Subroutine */ int setup_proc_info__(int *num_procs__, int *
	num_proc_rows__, int *num_proc_cols__, context *c);

/* Subroutine */ int initialize_mpi__(context *c);

/* Subroutine */ int setup_proc_info__(int *num_procs__, int *
	num_proc_rows__, int *num_proc_cols__, context *c);

/* Subroutine */ int setup_submatrix_info__(int *l2npcols, int *
	reduce_exch_proc__, int *reduce_send_starts__, int *
	reduce_send_lengths__, int *reduce_recv_starts__, int *
	reduce_recv_lengths__, context *c);

/* Subroutine */ int conj_grad__(int *colidx, int *rowstr, float 
	*x, float *z__, float *a, float *p, float *q, 
	float *r__, float *w, float *rnorm, int *l2npcols, 
	int *reduce_exch_proc__, int *reduce_send_starts__, int *
	reduce_send_lengths__, int *reduce_recv_starts__, int *
	reduce_recv_lengths__, context *c);

/* Subroutine */ int makea_(int *n, int *nz, float *a, int *
	colidx, int *rowstr, int *nonzer, int *firstrow, int *
	lastrow, int *firstcol, int *lastcol, float *rcond, 
	int *arow, int *acol, float *aelt, float *v, 
	int *iv, float *shift, context *c);

/* Subroutine */ int sparse_(float *a, int *colidx, int *rowstr, 
	int *n, int *arow, int *acol, float *aelt, int *
	firstrow, int *lastrow, float *x, int *mark, int *
	nzloc, int *nnza);

/* Subroutine */ int sprnvc_(int *n, int *nz, float *v, int *
	iv, int *nzloc, int *mark, context *c);

/* Subroutine */ int vecset_(int *n, float *v, int *iv, int *
	nzv, int *i__, float *val);

int icnvrt_(float *x, int *ipwr2);

/* Subroutine */ static int s_stop(char *a, ftnlen b) {
//    printf("Done\n\r");
	while (1);
}


/* randdp.f -- translated by f2c (version 20090411).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#define R23 (0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5)
#define T23 (2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0)
#define R46 (0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5)
#define T46 (2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0)

/* --------------------------------------------------------------------- */
/* --------------------------------------------------------------------- */
/*<       double precision function randlc (x, a) >*/
float randlc_(float *x, float *a)
{
    /* System generated locals */
    float ret_val;

    /* Local variables */
    float z__, a1, a2, t1, t2, t3, t4, x1, x2;

/* --------------------------------------------------------------------- */
/* --------------------------------------------------------------------- */
/* --------------------------------------------------------------------- */

/*   This routine returns a uniform pseudorandom double precision number in the */
/*   range (0, 1) by using the linear congruential generator */

/*   x_{k+1} = a x_k  (mod 2^46) */

/*   where 0 < x_k < 2^46 and 0 < a < 2^46.  This scheme generates 2^44 numbers */
/*   before repeating.  The argument A is the same as 'a' in the above formula, */
/*   and X is the same as x_0.  A and X must be odd double precision integers */
/*   in the range (1, 2^46).  The returned value RANDLC is normalized to be */
/*   between 0 and 1, i.e. RANDLC = 2^(-46) * x_1.  X is updated to contain */
/*   the new seed x_1, so that subsequent calls to RANDLC using the same */
/*   arguments will generate a continuous sequence. */

/*   This routine should produce the same results on any computer with at least */
/*   48 mantissa bits in double precision floating point data.  On 64 bit */
/*   systems, double precision should be disabled. */

/*   David H. Bailey     October 26, 1990 */

/* --------------------------------------------------------------------- */
/*<       implicit none >*/
/*<       double precision r23,r46,t23,t46,a,x,t1,t2,t3,t4,a1,a2,x1,x2,z >*/
/*<        >*/
/* --------------------------------------------------------------------- */
/*   Break A into two parts such that A = 2^23 * A1 + A2. */
/* --------------------------------------------------------------------- */
/*<       t1 = r23 * a >*/
    t1 = *a * R23;
/*<       a1 = int (t1) >*/
    a1 = (float) ((int) t1);
/*<       a2 = a - t23 * a1 >*/
    a2 = *a - a1 * T23;
/* --------------------------------------------------------------------- */
/*   Break X into two parts such that X = 2^23 * X1 + X2, compute */
/*   Z = A1 * X2 + A2 * X1  (mod 2^23), and then */
/*   X = 2^23 * Z + A2 * X2  (mod 2^46). */
/* --------------------------------------------------------------------- */
/*<       t1 = r23 * x >*/
    t1 = *x * R23;
/*<       x1 = int (t1) >*/
    x1 = (float) ((int) t1);
/*<       x2 = x - t23 * x1 >*/
    x2 = *x - x1 * T23;
/*<       t1 = a1 * x2 + a2 * x1 >*/
    t1 = a1 * x2 + a2 * x1;
/*<       t2 = int (r23 * t1) >*/
    t2 = (float) ((int) (t1 * R23));
/*<       z = t1 - t23 * t2 >*/
    z__ = t1 - t2 * T23;
/*<       t3 = t23 * z + a2 * x2 >*/
    t3 = z__ * T23 + a2 * x2;
/*<       t4 = int (r46 * t3) >*/
    t4 = (float) ((int) (t3 * R46));
/*<       x = t3 - t46 * t4 >*/
    *x = t3 - t4 * T46;
/*<       randlc = r46 * x >*/
    ret_val = *x * R46;
/*<       return >*/
    return ret_val;
/*<       end >*/
} /* randlc_ */


/* Subroutine */ int print_results__(char *name__, char *class__, int *n1,
         int *n2, int *n3, int *niter, int *nprocs_compiled__,
         int *nprocs_total__, float *t, float *mops, char * 
        optype, int *verified, char *npbversion, char *compiletime, char *
        cs1, char *cs2, char *cs3, char *cs4, char *cs5, char *cs6, char *cs7,
         ftnlen name_len, ftnlen class_len, ftnlen optype_len, ftnlen
        npbversion_len, ftnlen compiletime_len, ftnlen cs1_len, ftnlen
        cs2_len, ftnlen cs3_len, ftnlen cs4_len, ftnlen cs5_len, ftnlen
        cs6_len, ftnlen cs7_len)
{

    printf("Benchmark Completed.\n");
    printf(" Class = %c\n", *class__);
    printf(" Iterations = %d\n", *niter);
    printf(" Time in cycles = %f\n", *t);
    printf(" Total processes = %d\n", *nprocs_total__);
    printf(" Compiled procs  = %d\n", *nprocs_compiled__);
//    printf(" Mop/s total     = %f\n", *mops);
//    printf(" Mop/s/process   = %f\n", (float)*mops / (float)*nprocs_total__);
    printf(" Operation type  = %s\n", optype);
    if (*verified) {
        printf(" Verification     = SUCCESSFUL\n");
    } else {
        printf(" Verification     = UNSUCCESSFUL\n");
    }
    return 0;
} /* print_results__ */


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
  
#ifdef ARCH_MB
int nas_cg_mpi()
#else
int main(void)
#endif
{

    int i__1, i__2;
    float d__1;

    int l2npcols;
    int verified;
    int reduce_exch_proc__[NUM_PROC_COLS], i__, j, k;
    float t;
    float t1[6];
    int it;
    float *norm_temp1__;
    float *norm_temp2__;
    volatile float zeta_verify_value__;
    int reduce_send_starts__[NUM_PROC_COLS], reduce_recv_starts__[NUM_PROC_COLS];
    float err;
    int reduce_send_lengths__[NUM_PROC_COLS], reduce_recv_lengths__[NUM_PROC_COLS];
    float zeta, tmax, tsum[6];
    char class__[1];
    float tming[6], tmaxg[6], rnorm;
    float mflops;
    MPI_Status status;
    float epsilon;
    MPI_Request request;
    
    context c;

    c.main_flt_mem__.v = malloc((NA+1)*sizeof(float));
    c.main_flt_mem__.aelt = malloc(NZ*sizeof(float));
    c.main_flt_mem__.a = malloc(NZ*sizeof(float));
    c.main_flt_mem__.x = malloc((NA/NUM_PROC_ROWS+2)*sizeof(float));
    c.main_flt_mem__.z__ = malloc((NA/NUM_PROC_ROWS+2)*sizeof(float));
    c.main_flt_mem__.p = malloc((NA/NUM_PROC_ROWS+2)*sizeof(float));
    c.main_flt_mem__.q = malloc((NA/NUM_PROC_ROWS+2)*sizeof(float));
    c.main_flt_mem__.r__ = malloc((NA/NUM_PROC_ROWS+2)*sizeof(float));
    c.main_flt_mem__.w = malloc((NA/NUM_PROC_ROWS+2)*sizeof(float));
    c.main_int_mem__.colidx = malloc(NZ*sizeof(int));
    c.main_int_mem__.rowstr = malloc((NA+1)*sizeof(int)); 
    c.main_int_mem__.iv = malloc((2*NA+1)*sizeof(int));
    c.main_int_mem__.arow = malloc(NZ*sizeof(int));
    c.main_int_mem__.acol = malloc(NZ*sizeof(int));

    norm_temp2__ = malloc(2*sizeof(float));
    norm_temp1__ = malloc(2*sizeof(float));

    initialize_mpi__(&c);
    if (NA == 1400 && NONZER == 7 && NITER == 15 && SHIFT == 10) {
	*(unsigned char *)class__ = 'S';
	zeta_verify_value__ = 1.5047763586;
    } else if (NA == 7000 && NONZER == 8 && NITER == 15 && SHIFT == 12) {
	*(unsigned char *)class__ = 'W';
	zeta_verify_value__ = 12.762039;
	zeta_verify_value__ = 12.7812662125;
    } else if (NA == 14000 && NONZER == 11 && NITER == 15 && SHIFT == 20) {
	*(unsigned char *)class__ = 'A';
	zeta_verify_value__ = 17.130235054029;
	zeta_verify_value__ = 0.1117730737;
    } else if (NA == 75000 && NONZER == 13 && NITER == 75 && SHIFT == 60) {
	*(unsigned char *)class__ = 'B';
	zeta_verify_value__ = 22.712745482631;
    } else if (NA == 150000 && NONZER == 15 && NITER == 75 && SHIFT == 110) {
	*(unsigned char *)class__ = 'C';
	zeta_verify_value__ = 28.973605592845;
    } else if (NA == 1500000 && NONZER == 21 && NITER == 100 && SHIFT == 500) {
	*(unsigned char *)class__ = 'D';
	zeta_verify_value__ = 52.514532105794;
    } else if (NA == 9000000 && NONZER == 26 && NITER == 100) {
	*(unsigned char *)class__ = 'E';
	zeta_verify_value__ = 77.522164599383;
    } else {
	*(unsigned char *)class__ = 'U';
    }
    if (c.mpistuff_.me == c.mpistuff_.root) {
    }
//    if (TRUE_) {
	c.mpistuff_.dp_type__ = MPI_FLOAT;
//    } else {
//	c.mpistuff_.dp_type__ = MPI_REAL;
//    }
    c.partit_size__.naa = NA;
    c.partit_size__.nzz = NZ;
    int c__num_procs = NUM_PROCS;
    int c__num_proc_rows = NUM_PROC_ROWS;
    int c__num_proc_cols = NUM_PROC_COLS;
    setup_proc_info__(&c__num_procs, &c__num_proc_rows, &c__num_proc_cols, &c);
    setup_submatrix_info__(&l2npcols, reduce_exch_proc__, 
	    reduce_send_starts__, reduce_send_lengths__, reduce_recv_starts__,
	     reduce_recv_lengths__, &c);
    for (i__ = 1; i__ <= 4; ++i__) {
	timer_clear(i__, &c);
    }
    c.urando_.tran = 314159265.;
    c.urando_.amult = 11.;
    zeta = randlc_(&c.urando_.tran, &c.urando_.amult);
    int c__7 = NONZER;
    float c_b28 = SHIFT;
    float c_b34 = .1;
    makea_(&c.partit_size__.naa, &c.partit_size__.nzz, c.main_flt_mem__.a, 
	    c.main_int_mem__.colidx, c.main_int_mem__.rowstr, &c__7, &
	    c.partit_size__.firstrow, &c.partit_size__.lastrow, &
	    c.partit_size__.firstcol, &c.partit_size__.lastcol, &c_b34, 
	    c.main_int_mem__.arow, c.main_int_mem__.acol, c.main_flt_mem__.aelt, 
	    c.main_flt_mem__.v, c.main_int_mem__.iv, &c_b28, &c);
    i__1 = c.partit_size__.lastrow - c.partit_size__.firstrow + 1;

    for (j = 1; j <= i__1; ++j) {
	i__2 = c.main_int_mem__.rowstr[j] - 1;
	for (k = c.main_int_mem__.rowstr[j - 1]; k <= i__2; ++k) {
	    c.main_int_mem__.colidx[k - 1] = c.main_int_mem__.colidx[k - 1] - 
		    c.partit_size__.firstcol + 1;
	}
    }

    for (i__ = 1; i__ <= NA/NUM_PROC_ROWS+1; ++i__) {
	c.main_flt_mem__.x[i__ - 1] = 1.;
    }

    zeta = 0.;

    for (it = 1; it <= 1; ++it) {
	conj_grad__(c.main_int_mem__.colidx, c.main_int_mem__.rowstr, 
		c.main_flt_mem__.x, c.main_flt_mem__.z__, c.main_flt_mem__.a, 
		c.main_flt_mem__.p, c.main_flt_mem__.q, c.main_flt_mem__.r__, 
		c.main_flt_mem__.w, &rnorm, &l2npcols, reduce_exch_proc__, 
		reduce_send_starts__, reduce_send_lengths__, 
		reduce_recv_starts__, reduce_recv_lengths__, &c);
	norm_temp1__[0] = 0.;
	norm_temp1__[1] = 0.;
	i__1 = c.partit_size__.lastcol - c.partit_size__.firstcol + 1;
	for (j = 1; j <= i__1; ++j) {
	    norm_temp1__[0] += c.main_flt_mem__.x[j - 1] * c.main_flt_mem__.z__[
		    j - 1];
	    norm_temp1__[1] += c.main_flt_mem__.z__[j - 1] * 
		    c.main_flt_mem__.z__[j - 1];
	}
	i__1 = l2npcols;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (c.timers_.timeron) {
		timer_start(4, &c);
	    }
	    MPI_Irecv(norm_temp2__, 2, c.mpistuff_.dp_type__, reduce_exch_proc__[i__ - 1], i__, MPI_COMM_WORLD, &request)		    ;
	    MPI_Send(norm_temp1__, 2, c.mpistuff_.dp_type__, reduce_exch_proc__[i__ - 1], i__, MPI_COMM_WORLD);
	    MPI_Wait(&request, &status);
	    if (c.timers_.timeron) {
		timer_stop(4, &c);
	    }
	    norm_temp1__[0] += norm_temp2__[0];
	    norm_temp1__[1] += norm_temp2__[1];
	}
	norm_temp1__[1] = 1. / sqrtf(norm_temp1__[1]);
	i__1 = c.partit_size__.lastcol - c.partit_size__.firstcol + 1;
	for (j = 1; j <= i__1; ++j) {
	    c.main_flt_mem__.x[j - 1] = norm_temp1__[1] * c.main_flt_mem__.z__[
		    j - 1];
	}
    }


    for (i__ = 1; i__ <= NA/NUM_PROC_ROWS+1; ++i__) {
	c.main_flt_mem__.x[i__ - 1] = 1.;
    }
    zeta = 0.;
    for (i__ = 1; i__ <= 4; ++i__) {
	timer_clear(i__, &c);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    timer_clear(1, &c);
    timer_start(1, &c);
    for (it = 1; it <= NITER; ++it) {
	conj_grad__(c.main_int_mem__.colidx, c.main_int_mem__.rowstr, 
		c.main_flt_mem__.x, c.main_flt_mem__.z__, c.main_flt_mem__.a, 
		c.main_flt_mem__.p, c.main_flt_mem__.q, c.main_flt_mem__.r__, 
		c.main_flt_mem__.w, &rnorm, &l2npcols, reduce_exch_proc__, 
		reduce_send_starts__, reduce_send_lengths__, 
		reduce_recv_starts__, reduce_recv_lengths__, &c);
	norm_temp1__[0] = 0.;
	norm_temp1__[1] = 0.;
	i__1 = c.partit_size__.lastcol - c.partit_size__.firstcol + 1;
	for (j = 1; j <= i__1; ++j) {
	    norm_temp1__[0] += c.main_flt_mem__.x[j - 1] * c.main_flt_mem__.z__[
		    j - 1];
	    norm_temp1__[1] += c.main_flt_mem__.z__[j - 1] * 
		    c.main_flt_mem__.z__[j - 1];
	}
	i__1 = l2npcols;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (c.timers_.timeron) {
		timer_start(4, &c);
	    }
	    MPI_Irecv(norm_temp2__, 2, c.mpistuff_.dp_type__, reduce_exch_proc__[i__ - 1], i__, MPI_COMM_WORLD, &request)		    ;
	    MPI_Send(norm_temp1__, 2, c.mpistuff_.dp_type__, reduce_exch_proc__[i__ - 1], i__, MPI_COMM_WORLD);
	    MPI_Wait(&request, &status);
	    if (c.timers_.timeron) {
		timer_stop(4, &c);
	    }
	    norm_temp1__[0] += norm_temp2__[0];
	    norm_temp1__[1] += norm_temp2__[1];
	}

	norm_temp1__[1] = 1. / sqrtf(norm_temp1__[1]);
	if (c.mpistuff_.me == c.mpistuff_.root) {
	    zeta = 1. / norm_temp1__[0] + SHIFT;
	    if (it == 1) {
	    }
            printf("\t%d\t%f\t%f\n\r", it, rnorm, zeta);
	}
	i__1 = c.partit_size__.lastcol - c.partit_size__.firstcol + 1;
	for (j = 1; j <= i__1; ++j) {
	    c.main_flt_mem__.x[j - 1] = norm_temp1__[1] * c.main_flt_mem__.z__[
		    j - 1];
	}
    }
    timer_stop(1, &c);
    t = timer_read(1, &c);
    MPI_Reduce(&t, &tmax, 1, c.mpistuff_.dp_type__, MPI_MAX, c.mpistuff_.root, MPI_COMM_WORLD);
    if (c.mpistuff_.me == c.mpistuff_.root) {
        printf("Benchmark completed\n");
	epsilon = 1e-4;
	if (*(unsigned char *)class__ != 'U') {
	    err = (d__1 = zeta - zeta_verify_value__, abs(d__1)) / 
		    zeta_verify_value__;
printf("err=%f, epsilon=%f, zeta=%f, zeta_ver=%f\n\r", err, epsilon, zeta, zeta_verify_value__);
	    if (err <= epsilon) {
		verified = TRUE_;
	    } else {
		verified = FALSE_;
	    }
	} else {
	    verified = FALSE_;
	}
	if (tmax != 0.f) {
	    mflops = ( 2*NITER*NA ) * ( 3.+(float)( NONZER*(NONZER+1) ) + 25.*(5.+(float)( NONZER*(NONZER+1) )) + 3. ) / tmax / 1e6f;
	} else {
	    mflops = 0.f;
	}
        int c__1400 = NA;
        int c__15 = NITER;
        int c__1 = 1;
        int c__0 = 0;
	print_results__("CG", class__, &c__1400, &c__0, &c__0, &c__15, &c__1,
		 &c.mpistuff_.nprocs, &tmax, &mflops, "          floating poi"
		"nt", &verified, "3.3.1", "04 May 2012", "f77 ", "f77", "-lmpi"
		, "(none)", "-O3", "(none)", "randi8", (ftnlen)2, (ftnlen)1, (
		ftnlen)24, (ftnlen)5, (ftnlen)11, (ftnlen)4, (ftnlen)3, (
		ftnlen)5, (ftnlen)6, (ftnlen)3, (ftnlen)6, (ftnlen)6);
    }
    if (! c.timers_.timeron) {
	goto L999;
    }
    for (i__ = 1; i__ <= 4; ++i__) {
	t1[i__ - 1] = timer_read(i__, &c);
    }
    t1[1] -= t1[2];
    t1[5] = t1[2] + t1[3];
    t1[4] = t1[0] - t1[5];
    MPI_Reduce(t1, &tsum, 6, c.mpistuff_.dp_type__, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(t1, &tming, 6, c.mpistuff_.dp_type__, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(t1, &tmaxg, 6, c.mpistuff_.dp_type__, MPI_MAX, 0, MPI_COMM_WORLD);
    if (c.mpistuff_.me == 0) {
	for (i__ = 1; i__ <= 6; ++i__) {
	    tsum[i__ - 1] /= c.mpistuff_.nprocs;
	}
    }
L999:
    MPI_Finalize();
    return 0;
} /* MAIN__ */

/* Subroutine */ int initialize_mpi__(context *c)
{

    int fstatus;

    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &c->mpistuff_.me);
    MPI_Comm_size(MPI_COMM_WORLD, &c->mpistuff_.nprocs);
    c->mpistuff_.root = 0;
    if (c->mpistuff_.me == c->mpistuff_.root) {
	fstatus = -1;
	c->timers_.timeron = FALSE_;
	if (fstatus == 0) {
	    c->timers_.timeron = TRUE_;
	}
    }
    int *tmp;
    tmp = malloc(sizeof(int));
    *tmp = c->timers_.timeron;
    MPI_Bcast(tmp, 1, MPI_INT, 0, MPI_COMM_WORLD);
    c->timers_.timeron = *tmp;
	free(tmp);
    return 0;
} /* initialize_mpi__ */

/* Subroutine */ int setup_proc_info__(int *num_procs__, int *
	num_proc_rows__, int *num_proc_cols__, context *c)
{


    int i__, log2nprocs;

    if (c->mpistuff_.nprocs != *num_procs__) {
	if (c->mpistuff_.me == c->mpistuff_.root) {
       printf("Error with number of processors %d != %d\n\r", c->mpistuff_.nprocs, *num_procs__);
	}
	MPI_Finalize();
	s_stop("", (ftnlen)0);
    }
    i__ = *num_proc_cols__;
L100:
    if (i__ != 1 && i__ / 2 << 1 != i__) {
	if (c->mpistuff_.me == c->mpistuff_.root) {
            printf("Error with num_proc_cols\n");
	}
	MPI_Finalize();
	s_stop("", (ftnlen)0);
    }
    i__ /= 2;
    if (i__ != 0) {
	goto L100;
    }
    i__ = *num_proc_rows__;
L200:
    if (i__ != 1 && i__ / 2 << 1 != i__) {
	if (c->mpistuff_.me == c->mpistuff_.root) {
            printf("Error with num_proc_rows\n");
	}
	MPI_Finalize();
	s_stop("", (ftnlen)0);
    }
    i__ /= 2;
    if (i__ != 0) {
	goto L200;
    }
    log2nprocs = 0;
    i__ = c->mpistuff_.nprocs;
L300:
    if (i__ != 1 && i__ / 2 << 1 != i__) {
        printf("Error with nprocs\n");
	MPI_Finalize();
	s_stop("", (ftnlen)0);
    }
    i__ /= 2;
    if (i__ != 0) {
	++log2nprocs;
	goto L300;
    }
    c->partit_size__.npcols = *num_proc_cols__;
    c->partit_size__.nprows = *num_proc_rows__;
    return 0;
} /* setup_proc_info__ */


int setup_submatrix_info__(int *l2npcols, int *
	reduce_exch_proc__, int *reduce_send_starts__, int *
	reduce_send_lengths__, int *reduce_recv_starts__, int *
	reduce_recv_lengths__, context *c)
{
    int i__1;

    int col_size__, row_size__, i__, j, div_factor__;

    --reduce_recv_lengths__;
    --reduce_recv_starts__;
    --reduce_send_lengths__;
    --reduce_send_starts__;
    --reduce_exch_proc__;

    c->partit_size__.proc_row__ = c->mpistuff_.me / c->partit_size__.npcols;
    c->partit_size__.proc_col__ = c->mpistuff_.me - c->partit_size__.proc_row__ * 
	    c->partit_size__.npcols;
    if (c->partit_size__.naa / c->partit_size__.npcols * c->partit_size__.npcols == 
	    c->partit_size__.naa) {
	col_size__ = c->partit_size__.naa / c->partit_size__.npcols;
	c->partit_size__.firstcol = c->partit_size__.proc_col__ * col_size__ + 1;
	c->partit_size__.lastcol = c->partit_size__.firstcol - 1 + col_size__;
	row_size__ = c->partit_size__.naa / c->partit_size__.nprows;
	c->partit_size__.firstrow = c->partit_size__.proc_row__ * row_size__ + 1;
	c->partit_size__.lastrow = c->partit_size__.firstrow - 1 + row_size__;
    } else {
	if (c->partit_size__.proc_row__ < c->partit_size__.naa - 
		c->partit_size__.naa / c->partit_size__.nprows * 
		c->partit_size__.nprows) {
	    row_size__ = c->partit_size__.naa / c->partit_size__.nprows + 1;
	    c->partit_size__.firstrow = c->partit_size__.proc_row__ * row_size__ 
		    + 1;
	    c->partit_size__.lastrow = c->partit_size__.firstrow - 1 + row_size__;
	} else {
	    row_size__ = c->partit_size__.naa / c->partit_size__.nprows;
	    c->partit_size__.firstrow = (c->partit_size__.naa - 
		    c->partit_size__.naa / c->partit_size__.nprows * 
		    c->partit_size__.nprows) * (row_size__ + 1) + (
		    c->partit_size__.proc_row__ - (c->partit_size__.naa - 
		    c->partit_size__.naa / c->partit_size__.nprows * 
		    c->partit_size__.nprows)) * row_size__ + 1;
	    c->partit_size__.lastrow = c->partit_size__.firstrow - 1 + row_size__;
	}
	if (c->partit_size__.npcols == c->partit_size__.nprows) {
	    if (c->partit_size__.proc_col__ < c->partit_size__.naa - 
		    c->partit_size__.naa / c->partit_size__.npcols * 
		    c->partit_size__.npcols) {
		col_size__ = c->partit_size__.naa / c->partit_size__.npcols + 1;
		c->partit_size__.firstcol = c->partit_size__.proc_col__ * 
			col_size__ + 1;
		c->partit_size__.lastcol = c->partit_size__.firstcol - 1 + 
			col_size__;
	    } else {
		col_size__ = c->partit_size__.naa / c->partit_size__.npcols;
		c->partit_size__.firstcol = (c->partit_size__.naa - 
			c->partit_size__.naa / c->partit_size__.npcols * 
			c->partit_size__.npcols) * (col_size__ + 1) + (
			c->partit_size__.proc_col__ - (c->partit_size__.naa - 
			c->partit_size__.naa / c->partit_size__.npcols * 
			c->partit_size__.npcols)) * col_size__ + 1;
		c->partit_size__.lastcol = c->partit_size__.firstcol - 1 + 
			col_size__;
	    }
	} else {
	    if (c->partit_size__.proc_col__ / 2 < c->partit_size__.naa - 
		    c->partit_size__.naa / (c->partit_size__.npcols / 2) * (
		    c->partit_size__.npcols / 2)) {
		col_size__ = c->partit_size__.naa / (c->partit_size__.npcols / 2) 
			+ 1;
		c->partit_size__.firstcol = c->partit_size__.proc_col__ / 2 * 
			col_size__ + 1;
		c->partit_size__.lastcol = c->partit_size__.firstcol - 1 + 
			col_size__;
	    } else {
		col_size__ = c->partit_size__.naa / (c->partit_size__.npcols / 2);
		c->partit_size__.firstcol = (c->partit_size__.naa - 
			c->partit_size__.naa / (c->partit_size__.npcols / 2) * (
			c->partit_size__.npcols / 2)) * (col_size__ + 1) + (
			c->partit_size__.proc_col__ / 2 - (c->partit_size__.naa - 
			c->partit_size__.naa / (c->partit_size__.npcols / 2) * (
			c->partit_size__.npcols / 2))) * col_size__ + 1;
		c->partit_size__.lastcol = c->partit_size__.firstcol - 1 + 
			col_size__;
	    }
	    if (c->mpistuff_.me % 2 == 0) {
		c->partit_size__.lastcol = c->partit_size__.firstcol - 1 + (
			col_size__ - 1) / 2 + 1;
	    } else {
		c->partit_size__.firstcol = c->partit_size__.firstcol + (
			col_size__ - 1) / 2 + 1;
		c->partit_size__.lastcol = c->partit_size__.firstcol - 1 + 
			col_size__ / 2;
	    }
	}
    }
    if (c->partit_size__.npcols == c->partit_size__.nprows) {
	c->partit_size__.send_start__ = 1;
	c->partit_size__.send_len__ = c->partit_size__.lastrow - 
		c->partit_size__.firstrow + 1;
    } else {
	if (c->mpistuff_.me % 2 == 0) {
	    c->partit_size__.send_start__ = 1;
	    c->partit_size__.send_len__ = (c->partit_size__.lastrow + 1 - 
		    c->partit_size__.firstrow + 1) / 2;
	} else {
	    c->partit_size__.send_start__ = (c->partit_size__.lastrow + 1 - 
		    c->partit_size__.firstrow + 1) / 2 + 1;
	    c->partit_size__.send_len__ = (c->partit_size__.lastrow - 
		    c->partit_size__.firstrow + 1) / 2;
	}
    }
    if (c->partit_size__.npcols == c->partit_size__.nprows) {
	c->partit_size__.exch_proc__ = c->mpistuff_.me % c->partit_size__.nprows * 
		c->partit_size__.nprows + c->mpistuff_.me / c->partit_size__.nprows;
    } else {
	c->partit_size__.exch_proc__ = ((c->mpistuff_.me / 2 % 
		c->partit_size__.nprows * c->partit_size__.nprows + c->mpistuff_.me 
		/ 2 / c->partit_size__.nprows) << 1) + c->mpistuff_.me % 2;
    }
    i__ = c->partit_size__.npcols / 2;
    *l2npcols = 0;
    while(i__ > 0) {
	++(*l2npcols);
	i__ /= 2;
    }
    div_factor__ = c->partit_size__.npcols;
    i__1 = *l2npcols;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j = (c->partit_size__.proc_col__ + div_factor__ / 2) % div_factor__ + 
		c->partit_size__.proc_col__ / div_factor__ * div_factor__;
	reduce_exch_proc__[i__] = c->partit_size__.proc_row__ * 
		c->partit_size__.npcols + j;
	div_factor__ /= 2;
    }
    for (i__ = *l2npcols; i__ >= 1; --i__) {
	if (c->partit_size__.nprows == c->partit_size__.npcols) {
	    reduce_send_starts__[i__] = c->partit_size__.send_start__;
	    reduce_send_lengths__[i__] = c->partit_size__.send_len__;
	    reduce_recv_lengths__[i__] = c->partit_size__.lastrow - 
		    c->partit_size__.firstrow + 1;
	} else {
	    reduce_recv_lengths__[i__] = c->partit_size__.send_len__;
	    if (i__ == *l2npcols) {
		reduce_send_lengths__[i__] = c->partit_size__.lastrow - 
			c->partit_size__.firstrow + 1 - 
			c->partit_size__.send_len__;
		if (c->mpistuff_.me / 2 << 1 == c->mpistuff_.me) {
		    reduce_send_starts__[i__] = c->partit_size__.send_start__ + 
			    c->partit_size__.send_len__;
		} else {
		    reduce_send_starts__[i__] = 1;
		}
	    } else {
		reduce_send_lengths__[i__] = c->partit_size__.send_len__;
		reduce_send_starts__[i__] = c->partit_size__.send_start__;
	    }
	}
	reduce_recv_starts__[i__] = c->partit_size__.send_start__;
    }
    c->partit_size__.exch_recv_length__ = c->partit_size__.lastcol - 
	    c->partit_size__.firstcol + 1;
    return 0;
} 

inline void float_cp(float *a, float *b, int size) {
	int i;
	for (i = 0; i < size; i++) { a[i] = b[i]; }
}

/* Subroutine */ int conj_grad__(int *colidx, int *rowstr, float 
	*x, float *z__, float *a, float *p, float *q, 
	float *r__, float *w, float *rnorm, int *l2npcols, 
	int *reduce_exch_proc__, int *reduce_send_starts__, int *
	reduce_send_lengths__, int *reduce_recv_starts__, int *
	reduce_recv_lengths__, context *c)
{

    int cgitmax = 25;

    int i__1, i__2, i__3, i__4;

    float d__;
    int i__, j, k;
    float rho, sum, rho0;
    float beta;
    int cgit;
    float alpha;
    MPI_Request request;
    MPI_Status status;

    --colidx;
    --rowstr;
    --x;
    --z__;
    --a;
    --p;
    --q;
    --r__;
    --w;
    --reduce_recv_lengths__;
    --reduce_recv_starts__;
    --reduce_send_lengths__;
    --reduce_send_starts__;
    --reduce_exch_proc__;
    float *tmp0, *tmp1;

    tmp0 = malloc(TMP_SIZE*sizeof(float));
    tmp1 = malloc(TMP_SIZE*sizeof(float));

    if (c->timers_.timeron) {
	    timer_start(2, c);
    }
    i__1 = c->partit_size__.naa / c->partit_size__.nprows + 1;
    for (j = 1; j <= i__1; ++j) {
	q[j] = 0.;
	z__[j] = 0.;
	r__[j] = x[j];
	p[j] = r__[j];
	w[j] = 0.;
    }
    sum = 0.;
    i__1 = c->partit_size__.lastcol - c->partit_size__.firstcol + 1;
    for (j = 1; j <= i__1; ++j) {
	sum += r__[j] * r__[j];
    }
    i__1 = *l2npcols;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (c->timers_.timeron) {
	    timer_start(3, c);
	}
	MPI_Irecv(tmp0, 1, c->mpistuff_.dp_type__, reduce_exch_proc__[i__], i__, MPI_COMM_WORLD, &request);
        tmp1[0] = sum;
	MPI_Send(tmp1, 1, c->mpistuff_.dp_type__, reduce_exch_proc__[i__], i__, MPI_COMM_WORLD);
	MPI_Wait(&request, &status);
    rho = tmp0[0];
	if (c->timers_.timeron) {
	    timer_stop(3, c);
	}
	sum += rho;
    }
    rho = sum;
    i__1 = cgitmax;
    for (cgit = 1; cgit <= i__1; ++cgit) {
	i__2 = c->partit_size__.lastrow - c->partit_size__.firstrow + 1;
	for (j = 1; j <= i__2; ++j) {
	    sum = 0.;
	    i__3 = rowstr[j + 1] - 1;
	    for (k = rowstr[j]; k <= i__3; ++k) {
		sum += a[k] * p[colidx[k]];
	    }
	    w[j] = sum;
	}
	for (i__ = *l2npcols; i__ >= 1; --i__) {
	    if (c->timers_.timeron) {
		timer_start(3, c);
	    }
	    MPI_Irecv(tmp0, reduce_recv_lengths__[i__], c->mpistuff_.dp_type__, reduce_exch_proc__[i__], 		    i__, MPI_COMM_WORLD, &request);
            float_cp(tmp1, &w[reduce_send_starts__[i__]], reduce_send_lengths__[i__]);
	    MPI_Send(tmp1, reduce_send_lengths__[i__], c->mpistuff_.dp_type__, reduce_exch_proc__[i__], 		    i__, MPI_COMM_WORLD);
	    MPI_Wait(&request, &status);
            float_cp(&q[reduce_recv_starts__[i__]], tmp0, reduce_recv_lengths__[i__]);
	    if (c->timers_.timeron) {
		timer_stop(3, c);
	    }
	    i__2 = c->partit_size__.send_start__ + reduce_recv_lengths__[i__] - 
		    1;
	    for (j = c->partit_size__.send_start__; j <= i__2; ++j) {
		w[j] += q[j];
	    }
	}
	if (*l2npcols != 0) {
	    if (c->timers_.timeron) {
		timer_start(3, c);
	    }
            
	    MPI_Irecv(tmp0, c->partit_size__.exch_recv_length__,c->mpistuff_.dp_type__, c->partit_size__.exch_proc__, 1,MPI_COMM_WORLD, &request);
            float_cp(tmp1, &w[c->partit_size__.send_start__],c->partit_size__.send_len__);
	    MPI_Send(tmp1,c->partit_size__.send_len__, c->mpistuff_.dp_type__,c->partit_size__.exch_proc__, 1, MPI_COMM_WORLD);
	    MPI_Wait(&request, &status);
            float_cp(&q[1], tmp0, c->partit_size__.exch_recv_length__);
	    if (c->timers_.timeron) {
		timer_stop(3, c);
	    }
	} else {
	    i__2 = c->partit_size__.exch_recv_length__;
	    for (j = 1; j <= i__2; ++j) {
		q[j] = w[j];
	    }
	}
	i__3 = c->partit_size__.lastrow - c->partit_size__.firstrow + 1, i__4 = 
		c->partit_size__.lastcol - c->partit_size__.firstcol + 1;
	i__2 = max(i__3,i__4);
	for (j = 1; j <= i__2; ++j) {
	    w[j] = 0.;
	}
	sum = 0.;
	i__2 = c->partit_size__.lastcol - c->partit_size__.firstcol + 1;
	for (j = 1; j <= i__2; ++j) {
	    sum += p[j] * q[j];
	}
	i__2 = *l2npcols;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    if (c->timers_.timeron) {
		timer_start(3, c);
	    }
	    MPI_Irecv(tmp0, 1, c->mpistuff_.dp_type__, reduce_exch_proc__[i__], i__, MPI_COMM_WORLD, &request);
            tmp1[0] = sum;
	    MPI_Send(tmp1, 1, c->mpistuff_.dp_type__ ,reduce_exch_proc__[i__], i__, MPI_COMM_WORLD);
	    MPI_Wait(&request, &status);
            d__ = tmp0[0];
	    if (c->timers_.timeron) {
		timer_stop(3, c);
	    }
	    sum += d__;
	}
	d__ = sum;
	alpha = rho / d__;
	rho0 = rho;
	i__2 = c->partit_size__.lastcol - c->partit_size__.firstcol + 1;
	for (j = 1; j <= i__2; ++j) {
	    z__[j] += alpha * p[j];
	    r__[j] -= alpha * q[j];
	}
	sum = 0.;
	i__2 = c->partit_size__.lastcol - c->partit_size__.firstcol + 1;
	for (j = 1; j <= i__2; ++j) {
	    sum += r__[j] * r__[j];
	}
	i__2 = *l2npcols;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    if (c->timers_.timeron) {
		timer_start(3, c);
	    }
	    MPI_Irecv(tmp0, 1, c->mpistuff_.dp_type__, reduce_exch_proc__[i__], i__, MPI_COMM_WORLD, &request);
            tmp1[0] = sum;
	    MPI_Send(tmp1, 1, c->mpistuff_.dp_type__, reduce_exch_proc__[i__], i__, MPI_COMM_WORLD);
	    MPI_Wait(&request, &status);
            rho = tmp0[0];
	    if (c->timers_.timeron) {
		timer_stop(3, c);
	    }
	    sum += rho;
	}
	rho = sum;
	beta = rho / rho0;
	i__2 = c->partit_size__.lastcol - c->partit_size__.firstcol + 1;
	for (j = 1; j <= i__2; ++j) {
	    p[j] = r__[j] + beta * p[j];
	}
    }
    i__1 = c->partit_size__.lastrow - c->partit_size__.firstrow + 1;
    for (j = 1; j <= i__1; ++j) {
	sum = 0.;
	i__2 = rowstr[j + 1] - 1;
	for (k = rowstr[j]; k <= i__2; ++k) {
	    sum += a[k] * z__[colidx[k]];
	}
	w[j] = sum;
    }
    for (i__ = *l2npcols; i__ >= 1; --i__) {
	if (c->timers_.timeron) {
	    timer_start(3, c);
	}
	MPI_Irecv(tmp0, reduce_recv_lengths__[i__], c->mpistuff_.dp_type__, reduce_exch_proc__[i__], i__, MPI_COMM_WORLD, &request);
        float_cp(tmp1, &w[reduce_send_starts__[i__]], reduce_send_lengths__[i__]);
	MPI_Send(tmp1, reduce_send_lengths__[i__],c->mpistuff_.dp_type__, reduce_exch_proc__[i__], i__, MPI_COMM_WORLD);
	MPI_Wait(&request, &status);
        float_cp(&r__[reduce_recv_starts__[i__]], tmp0, reduce_recv_lengths__[i__]);
	if (c->timers_.timeron) {
	    timer_stop(3, c);
	}
	i__1 = c->partit_size__.send_start__ + reduce_recv_lengths__[i__] - 1;
	for (j = c->partit_size__.send_start__; j <= i__1; ++j) {
	    w[j] += r__[j];
	}
    }
    if (*l2npcols != 0) {
	if (c->timers_.timeron) {
	    timer_start(3, c);
	}
	MPI_Irecv(tmp0, c->partit_size__.exch_recv_length__, c->mpistuff_.dp_type__, c->partit_size__.exch_proc__, 1, MPI_COMM_WORLD, &request);
        float_cp(tmp1, &w[c->partit_size__.send_start__], c->partit_size__.send_len__);
	MPI_Send(tmp1, c->partit_size__.send_len__, c->mpistuff_.dp_type__, c->partit_size__.exch_proc__, 1, MPI_COMM_WORLD);
	MPI_Wait(&request, &status);
        float_cp(&r__[1], tmp0, c->partit_size__.exch_recv_length__);
	if (c->timers_.timeron) {
	    timer_stop(3, c);
	}
    } else {
	i__1 = c->partit_size__.exch_recv_length__;
	for (j = 1; j <= i__1; ++j) {
	    r__[j] = w[j];
	}
    }
    sum = 0.;
    i__1 = c->partit_size__.lastcol - c->partit_size__.firstcol + 1;
    for (j = 1; j <= i__1; ++j) {
	d__ = x[j] - r__[j];
	sum += d__ * d__;
    }
    i__1 = *l2npcols;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (c->timers_.timeron) {
	    timer_start(3, c);
	}
	MPI_Irecv(tmp0, 1, c->mpistuff_.dp_type__, reduce_exch_proc__[i__], i__, MPI_COMM_WORLD, &request);
        tmp1[0] = sum;
	MPI_Send(tmp1, 1, c->mpistuff_.dp_type__, reduce_exch_proc__[i__], i__, MPI_COMM_WORLD);
	MPI_Wait(&request, &status);
        d__ = tmp0[0];
	if (c->timers_.timeron) {
	    timer_stop(3, c);
	}
	sum += d__;
    }
    d__ = sum;
    if (c->mpistuff_.me == c->mpistuff_.root) {
	*rnorm = sqrtf(d__);
    }
    if (c->timers_.timeron) {
	timer_stop(2, c);
    }
    free(tmp0);
    free(tmp1);
    return 0;
} /* conj_grad__ */

/* Subroutine */ int makea_(int *n, int *nz, float *a, int *
	colidx, int *rowstr, int *nonzer, int *firstrow, int *
	lastrow, int *firstcol, int *lastcol, float *rcond, 
	int *arow, int *acol, float *aelt, float *v, 
	int *iv, float *shift, context *c)
{
    int i__1, i__2, i__3;
    float d__1;


    int i__, nzv, jcol, nnza;
    float size;
    int irow;
    float scale, ratio;
    int ivelt, ivelt1;
    int iouter;

//    cilist io___90 = { 0, 6, 0, 0, 0 };
//    cilist io___91 = { 0, 6, 0, 0, 0 };
//    cilist io___92 = { 0, 6, 0, 0, 0 };

    --iv;
    --v;
    --rowstr;
    --aelt;
    --acol;
    --arow;
    --colidx;
    --a;

    size = 1.;
    d__1 = 1. / (float) (*n);
    ratio = ar_float_pow(*rcond, d__1);
    nnza = 0;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	iv[*n + i__] = 0;
    }
    i__1 = *n;
    for (iouter = 1; iouter <= i__1; ++iouter) {
	nzv = *nonzer;
	sprnvc_(n, &nzv, &v[1], &colidx[1], &iv[1], &iv[*n + 1], c);
        float c_b220 = .5;
	vecset_(n, &v[1], &colidx[1], &nzv, &iouter, &c_b220);
	i__2 = nzv;
	for (ivelt = 1; ivelt <= i__2; ++ivelt) {
	    jcol = colidx[ivelt];
	    if (jcol >= *firstcol && jcol <= *lastcol) {
		scale = size * v[ivelt];
		i__3 = nzv;
		for (ivelt1 = 1; ivelt1 <= i__3; ++ivelt1) {
		    irow = colidx[ivelt1];
		    if (irow >= *firstrow && irow <= *lastrow) {
			++nnza;
			if (nnza > *nz) {
			    goto L9999;
			}
			acol[nnza] = jcol;
			arow[nnza] = irow;
			aelt[nnza] = v[ivelt1] * scale;
		    }
		}
	    }
	}
	size *= ratio;
    }
    i__1 = *lastrow;
    for (i__ = *firstrow; i__ <= i__1; ++i__) {
	if (i__ >= *firstcol && i__ <= *lastcol) {
	    iouter = *n + i__;
	    ++nnza;
	    if (nnza > *nz) {
		goto L9999;
	    }
	    acol[nnza] = i__;
	    arow[nnza] = i__;
	    aelt[nnza] = *rcond - *shift;
	}
    }
    sparse_(&a[1], &colidx[1], &rowstr[1], n, &arow[1], &acol[1], &aelt[1], 
	    firstrow, lastrow, &v[1], &iv[1], &iv[*n + 1], &nnza);
    return 0;
L9999:
    printf("Error with space for matrix elements in makea\n");
    return 0;
} /* makea_ */

/* Subroutine */ int sparse_(float *a, int *colidx, int *rowstr, 
	int *n, int *arow, int *acol, float *aelt, int *
	firstrow, int *lastrow, float *x, int *mark, int *
	nzloc, int *nnza)
{
    int i__1, i__2;

    int i__, j, k;
    float xi;
    int nza, jajp1, nrows, nzrow;

    --a;
    --colidx;
    --rowstr;
    --nzloc;
    --mark;
    --x;
    --arow;
    --acol;
    --aelt;

    nrows = *lastrow - *firstrow + 1;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	rowstr[j] = 0;
	mark[j] = FALSE_;
    }
    rowstr[*n + 1] = 0;
    i__1 = *nnza;
    for (nza = 1; nza <= i__1; ++nza) {
	j = arow[nza] - *firstrow + 2;
	++rowstr[j];
    }
    rowstr[1] = 1;
    i__1 = nrows + 1;
    for (j = 2; j <= i__1; ++j) {
	rowstr[j] += rowstr[j - 1];
    }
    i__1 = *nnza;
    for (nza = 1; nza <= i__1; ++nza) {
	j = arow[nza] - *firstrow + 1;
	k = rowstr[j];
	a[k] = aelt[nza];
	colidx[k] = acol[nza];
	++rowstr[j];
    }
    for (j = nrows; j >= 1; --j) {
	rowstr[j + 1] = rowstr[j];
    }
    rowstr[1] = 1;
    nza = 0;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x[i__] = 0.f;
	mark[i__] = FALSE_;
    }
    jajp1 = rowstr[1];
    i__1 = nrows;
    for (j = 1; j <= i__1; ++j) {
	nzrow = 0;
	i__2 = rowstr[j + 1] - 1;
	for (k = jajp1; k <= i__2; ++k) {
	    i__ = colidx[k];
	    x[i__] += a[k];
	    if (! mark[i__] && x[i__] != 0.) {
		mark[i__] = TRUE_;
		++nzrow;
		nzloc[nzrow] = i__;
	    }
	}
	i__2 = nzrow;
	for (k = 1; k <= i__2; ++k) {
	    i__ = nzloc[k];
	    mark[i__] = FALSE_;
	    xi = x[i__];
	    x[i__] = 0.;
	    if (xi != 0.) {
		++nza;
		a[nza] = xi;
		colidx[nza] = i__;
	    }
	}
	jajp1 = rowstr[j + 1];
	rowstr[j + 1] = nza + rowstr[1];
    }
    return 0;
} /* sparse_ */

/* Subroutine */ int sprnvc_(int *n, int *nz, float *v, int *
	iv, int *nzloc, int *mark, context *c)
{
    int i__1;

    int i__, ii, nn1, nzv, nzrow;
    float vecloc, vecelt;

    --mark;
    --nzloc;
    --v;
    --iv;

    nzv = 0;
    nzrow = 0;
    nn1 = 1;
L50:
    nn1 <<= 1;
    if (nn1 < *n) {
	goto L50;
    }
L100:
    if (nzv >= *nz) {
	goto L110;
    }
    vecelt = randlc_(&c->urando_.tran, &c->urando_.amult);
    vecloc = randlc_(&c->urando_.tran, &c->urando_.amult);
    i__ = icnvrt_(&vecloc, &nn1) + 1;
    if (i__ > *n) {
	goto L100;
    }
    if (mark[i__] == 0) {
	mark[i__] = 1;
	++nzrow;
	nzloc[nzrow] = i__;
	++nzv;
	v[nzv] = vecelt;
	iv[nzv] = i__;
    }
    goto L100;
L110:
    i__1 = nzrow;
    for (ii = 1; ii <= i__1; ++ii) {
	i__ = nzloc[ii];
	mark[i__] = 0;
    }
    return 0;
} /* sprnvc_ */

int icnvrt_(float *x, int *ipwr2)
{
    int ret_val;

    ret_val = (int) (*ipwr2 * *x);
    return ret_val;
} /* icnvrt_ */

/* Subroutine */ int vecset_(int *n, float *v, int *iv, int *
	nzv, int *i__, float *val)
{
    int i__1;

    int k;
    int set;

    --iv;
    --v;

    set = FALSE_;
    i__1 = *nzv;
    for (k = 1; k <= i__1; ++k) {
	if (iv[k] == *i__) {
	    v[k] = *val;
	    set = TRUE_;
	}
    }
    if (! set) {
	++(*nzv);
	v[*nzv] = *val;
	iv[*nzv] = *i__;
    }
    return 0;
} /* vecset_ */

