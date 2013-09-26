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
#define log    ar_float_log

#define MPI_Wtime FMPI_Wtime

#else
#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#endif

float ar_float_log(float x)
{
if (x <= 0.005) return -5.298317;
if (x <= 0.015) return -4.199705;
if (x <= 0.030) return -3.506558;
if (x <= 0.270) return -1.309333;
if (x <= 0.300) return -1.203973;
if (x <= 0.510) return -0.673345;
if (x <= 0.570) return -0.562119;
if (x <= 0.810) return -0.210721;
if (x <= 0.990) return -0.010051;
if (x <= 0.995) return -0.005013;
if (x <= 1.000) return -0.000001;

printf("Error: value too large\n");
while(1);
return 0;
}



// f2c
#define TRUE_ (1)
#define FALSE_ (0)
typedef int ftnlen;
typedef int integer;
typedef int logical;
#define max(a,b) ((a) >= (b) ? (a) : (b))
//
#define abs(a) ((a) >= 0 ? (a) : (-a))

#define CLS 'S'
#define NPROCS 16
#define NPM NPROCS

#if CLS == 'S'
#define CLASS "S"
#define M 24
#endif

#if CLS == 'W'
#define CLASS "W"
#define M 25
#endif

#if CLS == 'A'
#define CLASS "A"
#define M 28
#endif

void full_verify( void );

typedef struct context {
 struct {
     MPI_Datatype dp_type__;
     integer me, nprocs, root, dummy_dp_type__;
 } mpistuff_;

 struct {
//     float x[131072], q[10], qq[10000];
     float *x;
     float *q;
 } storage_;

 unsigned int start[5], elapsed[5];
} context;

static void timer_clear( int n, struct context *c )
{
    c->elapsed[n] = 0.0;
}

static void timer_start( int n, context *c )
{
    c->start[n] = MPI_Wtime();
}

static void timer_stop( int n, context *c )
{
    float t, now;

    now = MPI_Wtime();
    t = now - c->start[n];
    c->elapsed[n] += t;

}

static float timer_read( int n, context *c )
{
    return( (float)c->elapsed[n] );
}

/* Subroutine */ static int print_results__(char *name__, char *class__, int *n1,
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
    printf(" Operation type  = %s\n", optype);
    if (*verified) {
        printf(" Verification     = SUCCESSFUL\n");
    } else {
        printf(" Verification     = UNSUCCESSFUL\n");
    }
    return 0;
} /* print_results__ */


#define R23 (0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5)
#define T23 (2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0)
#define R46 (0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5)
#define T46 (2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0)

static float randlc_(float *x, float *a)
{
    float ret_val;

    float z__, a1, a2, t1, t2, t3, t4, x1, x2;

    t1 = *a * R23;
    a1 = (float) ((integer) t1);
    a2 = *a - a1 * T23;
    t1 = *x * R23;
    x1 = (float) ((integer) t1);
    x2 = *x - x1 * T23;
    t1 = a1 * x2 + a2 * x1;
    t2 = (float) ((integer) (t1 * R23));
    z__ = t1 - t2 * T23;
    t3 = z__ * T23 + a2 * x2;
    t4 = (float) ((integer) (t3 * R46));
    *x = t3 - t4 * T46;
    ret_val = *x * R46;
    return ret_val;
} /* randlc_ */

/* Subroutine */ int vranlc_(integer *n, float *x, float *a, 
	float *y)
{
    integer i__1;

    integer i__;
    float z__, a1, a2, t1, t2, t3, t4, x1, x2;

    --y;

    t1 = *a * R23;
    a1 = (float) ((integer) t1);
    a2 = *a - a1 * T23;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	t1 = *x * R23;
	x1 = (float) ((integer) t1);
	x2 = *x - x1 * T23;
	t1 = a1 * x2 + a2 * x1;
	t2 = (float) ((integer) (t1 * R23));
	z__ = t1 - t2 * T23;
	t3 = z__ * T23 + a2 * x2;
	t4 = (float) ((integer) (t3 * R46));
	*x = t3 - t4 * T46;
	y[i__] = *x * R46;
    }
    return 0;
} /* vranlc_ */


/* Subroutine */ static int s_stop(char *a, ftnlen b) {
//    printf("Done\n\r");
        while (1);
}

#ifdef ARCH_MB
int nas_ep_mpi()
#else
int main(void)
#endif
{

    float dum[3] = { 1.,1.,1. };
    integer i__1;
    float d__1, d__2;

    float sx_verify_value__, sy_verify_value__;
    integer ierrcode;
    logical verified;
    integer k_offset__, no_nodes__;
    integer i__, k, l;
    float t1, t2, t3, t4, x1, x2, gc, an;
    integer ik, kk, np;
    float tm;
    float tt;
    float sx;
    float sy;
    float t1m[6];
    integer nit;
    integer node;
    float mops, tsum[6];
    float tming[6], tmaxg[6];
    integer np_add__;
    float sx_err__, sy_err__;
    logical timers_enabled__;
    integer no_large_nodes__;

    context c;
    integer c__0 = 0;
    integer c_b39 = 131072;
    float c_b32 = 11.;
 
    int   *tmp;
    tmp = malloc(sizeof(int));

//     float x[131072], q[10], qq[10000];
    c.storage_.x = malloc(131072*sizeof(float));
    c.storage_.q = malloc(10*sizeof(float));
    printf("Total allocated size = %dKB\n", (int)(131082*sizeof(float))/1024);

     
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &node);
    MPI_Comm_size(MPI_COMM_WORLD, &no_nodes__);
    c.mpistuff_.root = 0;
    if (TRUE_) {
	c.mpistuff_.dp_type__ = MPI_FLOAT;
    } 
    if (node == c.mpistuff_.root) {
        printf(" NAS Parallel Benchmarks 3.3 -- EP Benchmark\n\n");
        printf(" Number of random numbers generated:        %d\n", 1<<M);
        printf(" Number of active processes:                      %d\n", no_nodes__);

        int fstatus = -1; 
	timers_enabled__ = FALSE_;
	if (fstatus == 0) {
	    timers_enabled__ = TRUE_;
	}
    }
    *tmp = timers_enabled__;
    MPI_Bcast(tmp, 1, MPI_INT, c.mpistuff_.root, MPI_COMM_WORLD);
    verified = FALSE_;
    np = (1<<(M - 16)) / no_nodes__;
    no_large_nodes__ = (1<<(M - 16)) % no_nodes__;
    if (node < no_large_nodes__) {
	np_add__ = 1;
    } else {
	np_add__ = 0;
    }
    np += np_add__;
    if (np == 0) {
        printf("Too many nodes:%d\n", no_nodes__);
	ierrcode = 1;
//	MPI_Abort(MPI_COMM_WORLD, ierrcode);
	s_stop("", (ftnlen)0);
    }
    vranlc_(&c__0, dum, &dum[1], &dum[2]);
    dum[0] = randlc_(&dum[1], &dum[2]);

//    for (i__ = 1; i__ <= 131072; ++i__) {
//	c.storage_.x[i__ - 1] = -1e99;
//    }
    mops = 0; //log(sqrtf(1.));
    for (i__ = 1; i__ <= 4; ++i__) {
	timer_clear(i__, &c);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    timer_start(1, &c);
    vranlc_(&c__0, &t1, &c_b32, &c.storage_.x[0]);
    for (i__ = 1; i__ <= 17; ++i__) {
	t2 = randlc_(&t1, &t1);
    }

    t1 = 11;

    an = t1;
    tt = 271828183.;
    gc = 0.;
    sx = 0.;
    sy = 0.;
    for (i__ = 0; i__ <= 9; ++i__) {
	c.storage_.q[i__] = 0.;
    }
    if (np_add__ == 1) {
	k_offset__ = node * np - 1;
    } else {
	k_offset__ = no_large_nodes__ * (np + 1) + (node - no_large_nodes__) *
		 np - 1;
    }
    i__1 = np;
    for (k = 1; k <= i__1; ++k) {
	kk = k_offset__ + k;
	t1 = 271828183.;
	t2 = an;
	for (i__ = 1; i__ <= 100; ++i__) {
	    ik = kk / 2;
	    if (ik << 1 != kk) {
		t3 = randlc_(&t1, &t2);
	    }
	    if (ik == 0) {
		goto L130;
	    }
	    t3 = randlc_(&t2, &t2);
	    kk = ik;
	}
L130:
	if (timers_enabled__) {
	    timer_start(3, &c);
	}
	vranlc_(&c_b39, &t1, &c_b32, &c.storage_.x[0]);
//printf("x=%f %f %f\n", c.storage_.x[0], c.storage_.x[1], c.storage_.x[2]);

	if (timers_enabled__) {
	    timer_stop(3, &c);
	}
	if (timers_enabled__) {
	    timer_start(2, &c);
	}
	for (i__ = 1; i__ <= 65536; ++i__) {
	    x1 = c.storage_.x[(i__ << 1) - 2] * 2. - 1.;
	    x2 = c.storage_.x[(i__ << 1) - 1] * 2. - 1.;
	    d__1 = x1;
	    d__2 = x2;
	    t1 = d__1 * d__1 + d__2 * d__2;
	    if (t1 <= 1.) {
		t2 = sqrtf(ar_float_log(t1) * -2. / t1);
		t3 = x1 * t2;
		t4 = x2 * t2;
		d__1 = abs(t3), d__2 = abs(t4);
		l = (integer) max(d__1,d__2);
		c.storage_.q[l] += 1.;
		sx += t3;
		sy += t4;
	    }
	}
    FMPI_Wtime();
    if (node == c.mpistuff_.root) printf("k=%d->%d\n", k, i__1);
	if (timers_enabled__) {
	    timer_stop(2, &c);
	}
    }
    if (timers_enabled__) {
	timer_start(4, &c);
    }
    MPI_Allreduce(&sx, &c.storage_.x[0], 1, c.mpistuff_.dp_type__, MPI_SUM, MPI_COMM_WORLD);
    sx = c.storage_.x[0];
    MPI_Allreduce(&sy, &c.storage_.x[0], 1, c.mpistuff_.dp_type__, MPI_SUM, MPI_COMM_WORLD);
    sy = c.storage_.x[0];
    MPI_Allreduce(c.storage_.q, &c.storage_.x[0], 10, c.mpistuff_.dp_type__, MPI_SUM, MPI_COMM_WORLD);
    if (timers_enabled__) {
	timer_stop(4, &c);
    }
    for (i__ = 1; i__ <= 10; ++i__) {
	c.storage_.q[i__ - 1] = c.storage_.x[i__ - 1];
    }
    for (i__ = 0; i__ <= 9; ++i__) {
	gc += c.storage_.q[i__];
    }
    timer_stop(1, &c);
    tm = timer_read(1, &c);
    MPI_Allreduce(&tm, &c.storage_.x[0], 1, c.mpistuff_.dp_type__, MPI_MAX, MPI_COMM_WORLD);
    tm = c.storage_.x[0];
    if (node == c.mpistuff_.root) {
	nit = 0;
	verified = TRUE_;
	if (M == 24) {
	    sx_verify_value__ = -3247.83465203474;
	    sy_verify_value__ = -6958.407078382297;
	} else if (M == 25) {
	    sx_verify_value__ = -2863.319731645753;
	    sy_verify_value__ = -6320.053679109499;
	} else if (M == 28) {
	    sx_verify_value__ = -4295.875165629892;
	    sy_verify_value__ = -15807.32573678431;
	} else if (M == 30) {
	    sx_verify_value__ = 40338.15542441498;
	    sy_verify_value__ = -26606.69192809235;
	} else if (M == 32) {
	    sx_verify_value__ = 47643.67927995374;
	    sy_verify_value__ = -80840.72988043731;
	} else if (M == 36) {
	    sx_verify_value__ = 198248.1200946593;
	    sy_verify_value__ = -102059.6636361769;
	} else if (M == 40) {
	    sx_verify_value__ = -531971.744153;
	    sy_verify_value__ = -368883.4557731;
	} else {
	    verified = FALSE_;
	}
	if (verified) {
	    sx_err__ = (d__1 = (sx - sx_verify_value__) / sx_verify_value__, 
		    abs(d__1));
	    sy_err__ = (d__1 = (sy - sy_verify_value__) / sy_verify_value__, 
		    abs(d__1));
	    verified = sx_err__ <= 1e-8 && sy_err__ <= 1e-8;
	}
	mops = 33554432. / tm / 1e6;
       
    int time_ovfl = mm_get_context(ar_get_core_id())->fmpi->time_ovfl;
    float time_cycles = (float)tm+(float)(0xFFFFFFFF)*(float)(time_ovfl);
    printf("EP Benchmark Results:\n\nCPU Cycles/Overflow = %f/%d\n\nCycles=%f\n N = 2^   %d\nNo. Gaussian Pairs =    %f\nSums =     %0.20f    %0.20f\nCounts:\n", tm, time_ovfl, time_cycles, M, gc, sx, sy);
	for (i__ = 0; i__ <= 9; ++i__) {
		printf("%d %10.f\n", i__, c.storage_.q[i__]);
	}

        integer c__25 = M+1;
        integer c__16 = NPM;
	print_results__("EP", CLASS, &c__25, &c__0, &c__0, &nit, &c__16, &
		no_nodes__, &tm, &mops, "Random numbers generated", &verified,
		 "3.3.1", "19 May 2012", "f77 ", "f77", "-lmpi", "(none)", 
		"-O3", "(none)", "randi8", (ftnlen)2, (ftnlen)1, (ftnlen)24, (
		ftnlen)5, (ftnlen)11, (ftnlen)4, (ftnlen)3, (ftnlen)5, (
		ftnlen)6, (ftnlen)3, (ftnlen)6, (ftnlen)6);
    }
    if (! timers_enabled__) {
	goto L999;
    }
    for (i__ = 1; i__ <= 4; ++i__) {
	t1m[i__ - 1] = timer_read(i__, &c);
    }
    t1m[5] = t1m[3];
    t1m[4] = t1m[0] - t1m[5];
    MPI_Reduce(t1m, &tsum, 6, c.mpistuff_.dp_type__, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(t1m, &tming, 6, c.mpistuff_.dp_type__, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(t1m, &tmaxg, 6, c.mpistuff_.dp_type__, MPI_MAX, 0, MPI_COMM_WORLD);
    if (node == 0) {
    }
L999:
    MPI_Finalize();
    return 0;
} /* MAIN__ */

