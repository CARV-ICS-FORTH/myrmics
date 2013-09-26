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
#include <math.h>
#endif

#define TRUE_ (1)
#define FALSE_ (0)
typedef int ftnlen;
typedef int integer;
typedef int logical;
typedef float real;
#define max(a,b) ((a) >= (b) ? (a) : (b))
#define abs(a) ((a) >= 0 ? (a) : (-a))

#define CLASS 'S'
#define NPROCS 16

#define NNODES_COMPILED NPROCS
#if NPROCS==1
#define NNODES_XDIM 1
//#define ISIZ1 12
//#define ISIZ2 12
#endif
#if NPROCS==2
#define NNODES_XDIM 2
//#define ISIZ1 6
//#define ISIZ2 12
#endif
#if NPROCS==4
#define NNODES_XDIM 2
//#define ISIZ1 6
//#define ISIZ2 6
#endif
#if NPROCS==8
#define NNODES_XDIM 4
//#define ISIZ1 6
//#define ISIZ2 3
#endif
#if NPROCS==16
#define NNODES_XDIM 4
//#define ISIZ1 3
//#define ISIZ2 3
#endif
#if NPROCS==32
#define NNODES_XDIM 8
//#define ISIZ1 2
//#define ISIZ2 3
#endif
#if NPROCS==64
#define NNODES_XDIM 8
//#define ISIZ1 2
//#define ISIZ2 2
#endif
#if NPROCS==128
#define NNODES_XDIM 16
//#define ISIZ1 1
//#define ISIZ2 2
#endif
#if NPROCS==256
#define NNODES_XDIM 16
//#define ISIZ1 1
//#define ISIZ2 1
#endif
#if NPROCS==512
#define NNODES_XDIM 32
//#define ISIZ1 1
//#define ISIZ2 1
#endif


#if CLASS == 'S'
#define ISIZ01 12
#define ISIZ02 12
#define ISIZ03 12
//#define ISIZ1 6
#define ISIZ1 12
#define ISIZ2 12
#define ISIZ3 ISIZ03
#define ITMAX_DEFAULT 50
#define INORM_DEFAULT 50
#define DT_DEFAULT 0.5
#endif

#if CLASS == 'W'
#define ISIZ01 33
#define ISIZ02 33
#define ISIZ03 33
#define ISIZ1 17
#define ISIZ2 33
#define ISIZ3 ISIZ03
#define ITMAX_DEFAULT 300
#define INORM_DEFAULT 300
#define DT_DEFAULT 0.0015
#endif

#if CLASS == 'A'
#define ISIZ01 64
#define ISIZ02 64
#define ISIZ03 64
#define ISIZ1 32
#define ISIZ2 64
#define ISIZ3 ISIZ03
#define ITMAX_DEFAULT 250
#define INORM_DEFAULT 250
#define DT_DEFAULT 2.0
#endif


typedef struct context {
 struct {
    MPI_Datatype dp_type__;
    integer me, nprocs, root;
 } mpistuff_;

 struct {
    float dxi, deta, dzeta, tx1, tx2, tx3, ty1, ty2, ty3, tz1, tz2, tz3;
    integer nx, ny, nz, nx0, ny0, nz0, ipt, ist, iend, jpt, jst, jend, ii1, 
	    ii2, ji1, ji2, ki1, ki2;
 } cgcon_;

 struct {
    float dx1, dx2, dx3, dx4, dx5, dy1, dy2, dy3, dy4, dy5, dz1, dz2, 
	    dz3, dz4, dz5, dssp;
 } disp_;

 struct {
//    float u[5*(ISIZ1+4)*(ISIZ2+4)*ISIZ3]/* was [5][10][16][12] */, rsd[5*(ISIZ1+4)*(ISIZ2+4)*ISIZ3]	/* was [5][10]
//	    [16][12] */, frct[5*(ISIZ1+4)*(ISIZ2+4)*ISIZ3]	/* was [5][10][16][12] */, flux[5*(ISIZ1+2)*(ISIZ2+2)*ISIZ3]	
    float *u, *rsd, *frct, *flux;
 } cvar_;

 struct {
    integer ipr, inorm;
 } cprcon_;

 struct {
    float dt, omega, *tolrsd, rsdnm[5], *errnm, frc, ttotal;
    integer itmax, invert;
 } ctscon_;

 struct {
//    float a[5*5*ISIZ1*ISIZ2]	/* was [5][5][6][12] */, b[5*5*ISIZ1*ISIZ2]	/* was [5][5][
//	    6][12] */, c__[5*5*ISIZ1*ISIZ2]	/* was [5][5][6][12] */, d__[5*5*ISIZ1*ISIZ2]	
//	    /* was [5][5][6][12] */;
    float *a, *b, *c__, *d__;
 } cjac_;

 struct {
//    float ce[65]	/* was [5][13] */;
    float *ce;
 } cexact_;

 struct {
    integer id, ndim, num, xdim, ydim, row, col;
 } dim_;

 struct {
    integer north, south, east, west;
 } neigh_;

 struct {
//    float buf[5*2*ISIZ2*ISIZ3]	/* was [5][288] */, buf1[5*2*ISIZ2*ISIZ3]	/* was [5][
//	    288] */;
    float *buf, *buf1;
 } comm_;

 struct {
    float maxtime;
    logical timeron;
 } timer_;

 unsigned int start[11], elapsed[11];

 float *ftmp;

} context;


#ifdef ARCH_MB
void MPI_Abort(MPI_Comm comm, int error) {
   printf("MPI Abort!\n");
   while(1);
}
#endif

int rhs_(context *c);
int l2norm_(integer *ldx, integer *ldy, integer *ldz, 
	integer *nx0, integer *ny0, integer *nz0, integer *ist, integer *iend,
	 integer *jst, integer *jend, float *v, float *sum, context *c);
int jacld_(integer *k, context *c);
int blts_(integer *ldmx, integer *ldmy, integer *ldmz, 
	integer *nx, integer *ny, integer *nz, integer *k, float *omega, 
	float *v, float *ldz, float *ldy, float *ldx, 
	float *d__, integer *ist, integer *iend, integer *jst, integer *
	jend, integer *nx0, integer *ny0, integer *ipt, integer *jpt, context *c);
int jacu_(integer *k, context *c);
int buts_(integer *ldmx, integer *ldmy, integer *ldmz, 
	integer *nx, integer *ny, integer *nz, integer *k, float *omega, 
	float *v, float *tv, float *d__, float *udx, 
	float *udy, float *udz, integer *ist, integer *iend, 
	integer *jst, integer *jend, integer *nx0, integer *ny0, integer *ipt,
	 integer *jpt, context *c);
int exact_(integer *i__, integer *j, integer *k, float *
	u000ijk, context *c);
int exchange_1__(float *g, integer *k, integer *iex, context *c);
int exchange_3__(float *g, integer *iex, context *c);
int exchange_4__(float *g, float *h__, integer *
	ibeg, integer *ifin1, integer *jbeg, integer *jfin1, context *c);
int exchange_5__(float *g, integer *ibeg, integer *
	ifin1, context *c);
int exchange_6__(float *g, integer *jbeg, integer *
	jfin1, context *c);
int bcast_inputs__(context *c);
int verify_(float *xcr, float *xce, float *
	xci, char *class__, logical *verified, ftnlen class_len, context *c);

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

/* Subroutine */ int verify_(float *xcr, float *xce, float *
	xci, char *class__, logical *verified, ftnlen class_len, context *c)
{

    float d__1;

    integer m;
    float dtref, xcedif[5], xcidif, xceref[5], xcrdif[5], xciref, 
	    xcrref[5], epsilon;

    dtref = 0.0;
    --xce;
    --xcr;

    epsilon = 1e-8;
    *(unsigned char *)class__ = 'U';
    *verified = TRUE_;
    for (m = 1; m <= 5; ++m) {
	xcrref[m - 1] = 1.f;
	xceref[m - 1] = 1.f;
    }
    xciref = 1.f;
    if (c->cgcon_.nx0 == 12 && c->cgcon_.ny0 == 12 && c->cgcon_.nz0 == 12 && 
	    c->ctscon_.itmax == 50) {
	*(unsigned char *)class__ = 'S';
	dtref = .5;
	xcrref[0] = .016196343210976702;
	xcrref[1] = .0021976745164821318;
	xcrref[2] = .0015179927653399185;
	xcrref[3] = .0015029584435994323;
	xcrref[4] = .034264073155896461;
	xceref[0] = 6.4223319957960924e-4;
	xceref[1] = 8.4144342047347926e-5;
	xceref[2] = 5.8588269616485186e-5;
	xceref[3] = 5.847422259515735e-5;
	xceref[4] = .0013103347914111294;
	xciref = 7.8418928865937083;
    } else if (c->cgcon_.nx0 == 33 && c->cgcon_.ny0 == 33 && c->cgcon_.nz0 == 33 && 
	    c->ctscon_.itmax == 300) {
	*(unsigned char *)class__ = 'W';
	dtref = .0015;
	xcrref[0] = 12.36511638192;
	xcrref[1] = 1.317228477799;
	xcrref[2] = 2.550120713095;
	xcrref[3] = 2.326187750252;
	xcrref[4] = 28.26799444189;
	xceref[0] = .4867877144216;
	xceref[1] = .05064652880982;
	xceref[2] = .0928181810196;
	xceref[3] = .08570126542733;
	xceref[4] = 1.084277417792;
	xciref = 11.61399311023;
    } else if (c->cgcon_.nx0 == 64 && c->cgcon_.ny0 == 64 && c->cgcon_.nz0 == 64 && 
	    c->ctscon_.itmax == 250) {
	*(unsigned char *)class__ = 'A';
	dtref = 2.;
	xcrref[0] = 779.02107606689367;
	xcrref[1] = 63.40276525969287;
	xcrref[2] = 194.99249727292479;
	xcrref[3] = 178.45301160418537;
	xcrref[4] = 1838.4760349464247;
	xceref[0] = 29.964085685471943;
	xceref[1] = 2.8194576365003349;
	xceref[2] = 7.3473412698774742;
	xceref[3] = 6.7139225687777051;
	xceref[4] = 70.715315688392578;
	xciref = 26.030925604886277;
    } else if (c->cgcon_.nx0 == 102 && c->cgcon_.ny0 == 102 && c->cgcon_.nz0 == 102 
	    && c->ctscon_.itmax == 250) {
	*(unsigned char *)class__ = 'B';
	dtref = 2.;
	xcrref[0] = 3553.2672969982736;
	xcrref[1] = 262.14750795310692;
	xcrref[2] = 883.3372185095219;
	xcrref[3] = 778.12774739425265;
	xcrref[4] = 7308.7969592545314;
	xceref[0] = 114.01176380212709;
	xceref[1] = 8.1098963655421574;
	xceref[2] = 28.480597317698308;
	xceref[3] = 25.905394567832939;
	xceref[4] = 260.54907504857413;
	xciref = 47.887162703308227;
    } else if (c->cgcon_.nx0 == 162 && c->cgcon_.ny0 == 162 && c->cgcon_.nz0 == 162 
	    && c->ctscon_.itmax == 250) {
	*(unsigned char *)class__ = 'C';
	dtref = 2.;
	xcrref[0] = 10376.6980323537846;
	xcrref[1] = 892.212458801008552;
	xcrref[2] = 2562.38814582660871;
	xcrref[3] = 2191.94343857831427;
	xcrref[4] = 17807.8057261061185;
	xceref[0] = 215.986399716949279;
	xceref[1] = 15.57895592398636;
	xceref[2] = 54.1318863077207766;
	xceref[3] = 48.2262643154045421;
	xceref[4] = 455.902910043250358;
	xciref = 66.64045535721813;
    } else if (c->cgcon_.nx0 == 408 && c->cgcon_.ny0 == 408 && c->cgcon_.nz0 == 408 
	    && c->ctscon_.itmax == 300) {
	*(unsigned char *)class__ = 'D';
	dtref = 1.;
	xcrref[0] = 48684.17937025;
	xcrref[1] = 4696.371050071;
	xcrref[2] = 12181.14549776;
	xcrref[3] = 10338.01493461;
	xcrref[4] = 71423.98413817;
	xceref[0] = 375.2393004482;
	xceref[1] = 30.84128893659;
	xceref[2] = 94.34276905469;
	xceref[3] = 82.30686681928;
	xceref[4] = 700.262063621;
	xciref = 83.34101392503;
    } else if (c->cgcon_.nx0 == 1020 && c->cgcon_.ny0 == 1020 && c->cgcon_.nz0 == 
	    1020 && c->ctscon_.itmax == 300) {
	*(unsigned char *)class__ = 'E';
	dtref = .5;
	xcrref[0] = 209964.1687874;
	xcrref[1] = 21304.03143165;
	xcrref[2] = 53192.28789371;
	xcrref[3] = 45097.61639833;
	xcrref[4] = 293236.000659;
	xceref[0] = 480.0572578333;
	xceref[1] = 42.21993400184;
	xceref[2] = 121.0851906824;
	xceref[3] = 104.788898677;
	xceref[4] = 836.3028257389;
	xciref = 95.12163272273;
    } else {
	*verified = FALSE_;
    }
    for (m = 1; m <= 5; ++m) {
	xcrdif[m - 1] = (d__1 = (xcr[m] - xcrref[m - 1]) / xcrref[m - 1], abs(
		d__1));
	xcedif[m - 1] = (d__1 = (xce[m] - xceref[m - 1]) / xceref[m - 1], abs(
		d__1));
    }
    xcidif = (d__1 = (*xci - xciref) / xciref, abs(d__1));
    if (*(unsigned char *)class__ != 'U') {
	*verified = (d__1 = c->ctscon_.dt - dtref, abs(d__1)) <= epsilon;
	if (! (*verified)) {
	    *(unsigned char *)class__ = 'U';
	}
    } else {
    }
    if (*(unsigned char *)class__ != 'U') {
      printf("Comparison of RMS-norms of residual\n");
    } else {
    }
    for (m = 1; m <= 5; ++m) {
	if (*(unsigned char *)class__ == 'U') {
	} else if (xcrdif[m - 1] <= epsilon) {
	} else {
	    *verified = FALSE_;
            printf("FAILURE:  %d %f %f %f\n", m, xcr[m], xcrref[m-1], xcrdif[m-1]);
	}
    }
    if (*(unsigned char *)class__ != 'U') {
        printf("Comparison of RMS-norms of solution error\n");
    } else {
    }
    for (m = 1; m <= 5; ++m) {
	if (*(unsigned char *)class__ == 'U') {
	} else if (xcedif[m - 1] <= epsilon) {
	} else {
	    *verified = FALSE_;
            printf("FAILURE: %d %f %f %f\n", m, xce[m], xceref[m-1], xcedif[m-1]);
	}
    }
    if (*(unsigned char *)class__ != 'U') {
    } else {
    }
    if (*(unsigned char *)class__ == 'U') {
    } else if (xcidif <= epsilon) {
    } else {
	*verified = FALSE_;
        printf("FAILURE: %f %f %f\n", *xci, xciref, xcidif);
    }
    if (*(unsigned char *)class__ == 'U') {
    } else if (*verified) {
    } else {
    }
    return 0;
} /* verify_ */

/* Subroutine */ int subdomain_(context *c)
{

    integer errorcode, mm;

    mm = c->cgcon_.nx0 % c->dim_.xdim;
    if (c->dim_.row <= mm) {
	c->cgcon_.nx = c->cgcon_.nx0 / c->dim_.xdim + 1;
	c->cgcon_.ipt = (c->dim_.row - 1) * c->cgcon_.nx;
    } else {
	c->cgcon_.nx = c->cgcon_.nx0 / c->dim_.xdim;
	c->cgcon_.ipt = (c->dim_.row - 1) * c->cgcon_.nx + mm;
    }
    mm = c->cgcon_.ny0 % c->dim_.ydim;
    if (c->dim_.col <= mm) {
	c->cgcon_.ny = c->cgcon_.ny0 / c->dim_.ydim + 1;
	c->cgcon_.jpt = (c->dim_.col - 1) * c->cgcon_.ny;
    } else {
	c->cgcon_.ny = c->cgcon_.ny0 / c->dim_.ydim;
	c->cgcon_.jpt = (c->dim_.col - 1) * c->cgcon_.ny + mm;
    }
    c->cgcon_.nz = c->cgcon_.nz0;
    if (c->cgcon_.nx < 3 || c->cgcon_.ny < 3 || c->cgcon_.nz < 3) {
	errorcode = 1;
	MPI_Abort(MPI_COMM_WORLD, errorcode);
    }
    if (c->cgcon_.nx > ISIZ1 || c->cgcon_.ny > ISIZ2 || c->cgcon_.nz > ISIZ3) {
	errorcode = 1;
	MPI_Abort(MPI_COMM_WORLD, errorcode);
    }
    c->cgcon_.ist = 1;
    c->cgcon_.iend = c->cgcon_.nx;
    if (c->neigh_.north == -1) {
	c->cgcon_.ist = 2;
    }
    if (c->neigh_.south == -1) {
	c->cgcon_.iend = c->cgcon_.nx - 1;
    }
    c->cgcon_.jst = 1;
    c->cgcon_.jend = c->cgcon_.ny;
    if (c->neigh_.west == -1) {
	c->cgcon_.jst = 2;
    }
    if (c->neigh_.east == -1) {
	c->cgcon_.jend = c->cgcon_.ny - 1;
    }
    return 0;
} /* subdomain_ */

/* Subroutine */ int ssor_(integer *niter, context *c)
{

    integer i__1, i__2, i__3, i__4;

    integer i__, j, k, m;
    float tv[5*ISIZ1*ISIZ2]	/* was [5][6][12] */;
    float tmp;
    integer istep;
    float wtime;
    float delunm[5];

    c->mpistuff_.root = 0;
    tmp = 1. / (c->ctscon_.omega * (2. - c->ctscon_.omega));
    for (m = 1; m <= ISIZ2; ++m) {
	for (k = 1; k <= ISIZ1; ++k) {
	    for (j = 1; j <= 5; ++j) {
		for (i__ = 1; i__ <= 5; ++i__) {
		    c->cjac_.a[i__ + (j + (k + m * 6) * 5) * 5 - 181] = 0.;
		    c->cjac_.b[i__ + (j + (k + m * 6) * 5) * 5 - 181] = 0.;
		    c->cjac_.c__[i__ + (j + (k + m * 6) * 5) * 5 - 181] = 0.;
		    c->cjac_.d__[i__ + (j + (k + m * 6) * 5) * 5 - 181] = 0.;
		}
	    }
	}
    }
    rhs_(c);
    int c__6 = ISIZ1;
    int c__12 = ISIZ2;
    int isiz3 = ISIZ3;
    l2norm_(&c__6, &c__12, &isiz3, &c->cgcon_.nx0, &c->cgcon_.ny0, &c->cgcon_.nz0, &
	    c->cgcon_.ist, &c->cgcon_.iend, &c->cgcon_.jst, &c->cgcon_.jend, 
	    c->cvar_.rsd, c->ctscon_.rsdnm, c);
    for (i__ = 1; i__ <= 10; ++i__) {
	timer_clear(i__, c);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    timer_clear(1, c);
    timer_start(1, c);
    i__1 = *niter;
    for (istep = 1; istep <= i__1; ++istep) {
	if (c->dim_.id == 0) {
	    if (istep % 20 == 0 || istep == c->ctscon_.itmax || istep == 1) {
		if (*niter > 1) {
		}
	    }
	}
	i__2 = c->cgcon_.nz - 1;
	for (k = 2; k <= i__2; ++k) {
	    i__3 = c->cgcon_.jend;
	    for (j = c->cgcon_.jst; j <= i__3; ++j) {
		i__4 = c->cgcon_.iend;
		for (i__ = c->cgcon_.ist; i__ <= i__4; ++i__) {
		    for (m = 1; m <= 5; ++m) {
			c->cvar_.rsd[m + (i__ + (j + (k << 4)) * 10) * 5 - 746] 
				= c->ctscon_.dt * c->cvar_.rsd[m + (i__ + (j + (k 
				<< 4)) * 10) * 5 - 746];
		    }
		}
	    }
	}
	i__2 = c->cgcon_.nz - 1;
	for (k = 2; k <= i__2; ++k) {
	    jacld_(&k, c);
	    blts_(&c__6, &c__12, &isiz3, &c->cgcon_.nx, &c->cgcon_.ny, &
		    c->cgcon_.nz, &k, &c->ctscon_.omega, c->cvar_.rsd, c->cjac_.a, 
		    c->cjac_.b, c->cjac_.c__, c->cjac_.d__, &c->cgcon_.ist, &
		    c->cgcon_.iend, &c->cgcon_.jst, &c->cgcon_.jend, &c->cgcon_.nx0, &
		    c->cgcon_.ny0, &c->cgcon_.ipt, &c->cgcon_.jpt, c);
	}
	for (k = c->cgcon_.nz - 1; k >= 2; --k) {
	    jacu_(&k, c);
	    buts_(&c__6, &c__12, &isiz3, &c->cgcon_.nx, &c->cgcon_.ny, &
		    c->cgcon_.nz, &k, &c->ctscon_.omega, c->cvar_.rsd, tv, 
		    c->cjac_.d__, c->cjac_.a, c->cjac_.b, c->cjac_.c__, &c->cgcon_.ist, 
		    &c->cgcon_.iend, &c->cgcon_.jst, &c->cgcon_.jend, &c->cgcon_.nx0, 
		    &c->cgcon_.ny0, &c->cgcon_.ipt, &c->cgcon_.jpt, c);
	}
	i__2 = c->cgcon_.nz - 1;
	for (k = 2; k <= i__2; ++k) {
	    i__3 = c->cgcon_.jend;
	    for (j = c->cgcon_.jst; j <= i__3; ++j) {
		i__4 = c->cgcon_.iend;
		for (i__ = c->cgcon_.ist; i__ <= i__4; ++i__) {
		    for (m = 1; m <= 5; ++m) {
			c->cvar_.u[m + (i__ + (j + (k << 4)) * 10) * 5 - 746] +=
				 tmp * c->cvar_.rsd[m + (i__ + (j + (k << 4)) * 
				10) * 5 - 746];
		    }
		}
	    }
	}
	if (istep % c->cprcon_.inorm == 0) {
	    l2norm_(&c__6, &c__12, &isiz3, &c->cgcon_.nx0, &c->cgcon_.ny0, &
		    c->cgcon_.nz0, &c->cgcon_.ist, &c->cgcon_.iend, &c->cgcon_.jst, &
		    c->cgcon_.jend, c->cvar_.rsd, delunm, c);
	}
	rhs_(c);
	if (istep % c->cprcon_.inorm == 0 || istep == c->ctscon_.itmax) {
	    l2norm_(&c__6, &c__12, &isiz3, &c->cgcon_.nx0, &c->cgcon_.ny0, &
		    c->cgcon_.nz0, &c->cgcon_.ist, &c->cgcon_.iend, &c->cgcon_.jst, &
		    c->cgcon_.jend, c->cvar_.rsd, c->ctscon_.rsdnm, c);
	}
	if (c->ctscon_.rsdnm[0] < c->ctscon_.tolrsd[0] && c->ctscon_.rsdnm[1] < 
		c->ctscon_.tolrsd[1] && c->ctscon_.rsdnm[2] < c->ctscon_.tolrsd[2] 
		&& c->ctscon_.rsdnm[3] < c->ctscon_.tolrsd[3] && c->ctscon_.rsdnm[4]
		 < c->ctscon_.tolrsd[4]) {
	    if (c->dim_.id == 0) {
	    }
	    goto L900;
	}
    MPI_Wtime();
    }
L900:
    timer_stop(1, c);
    wtime = timer_read(1, c);
    c->ftmp[0] = c->timer_.maxtime;
    MPI_Allreduce(&wtime, &c->ftmp[0],  1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
    c->timer_.maxtime = c->ftmp[0];
#ifdef ARCH_MB
    printf("Time ovfl=%d\n", mm_get_context(ar_get_core_id())->fmpi->time_ovfl);
#endif
//    MPI_Allreduce(&wtime, &c->timer_.maxtime, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
    return 0;
} /* ssor_ */

/* Subroutine */ int setiv_(context *c)
{
    integer i__1, i__2, i__3;

    float ue_nx0jk__[5], ue_iny0k__[5];
    integer i__, j, k, m;
    float xi, eta, pxi, peta, zeta;
    integer iglob, jglob;
    float pzeta, ue_ij1__[5], ue_i1k__[5], ue_1jk__[5], ue_ijnz__[ 5];

    integer c__1 = 1;
    i__1 = c->cgcon_.nz - 1;
    for (k = 2; k <= i__1; ++k) {
	zeta = (float) (k - 1) / (c->cgcon_.nz - 1);
	i__2 = c->cgcon_.ny;
	for (j = 1; j <= i__2; ++j) {
	    jglob = c->cgcon_.jpt + j;
	    if (jglob != 1 && jglob != c->cgcon_.ny0) {
		eta = (float) (jglob - 1) / (c->cgcon_.ny0 - 1);
		i__3 = c->cgcon_.nx;
		for (i__ = 1; i__ <= i__3; ++i__) {
		    iglob = c->cgcon_.ipt + i__;
		    if (iglob != 1 && iglob != c->cgcon_.nx0) {
			xi = (float) (iglob - 1) / (c->cgcon_.nx0 - 1);
			exact_(&c__1, &jglob, &k, ue_1jk__, c);
			exact_(&c->cgcon_.nx0, &jglob, &k, ue_nx0jk__, c);
			exact_(&iglob, &c__1, &k, ue_i1k__, c);
			exact_(&iglob, &c->cgcon_.ny0, &k, ue_iny0k__, c);
			exact_(&iglob, &jglob, &c__1, ue_ij1__, c);
			exact_(&iglob, &jglob, &c->cgcon_.nz, ue_ijnz__, c);
			for (m = 1; m <= 5; ++m) {
			    pxi = (1. - xi) * ue_1jk__[m - 1] + xi * 
				    ue_nx0jk__[m - 1];
			    peta = (1. - eta) * ue_i1k__[m - 1] + eta * 
				    ue_iny0k__[m - 1];
			    pzeta = (1. - zeta) * ue_ij1__[m - 1] + zeta * 
				    ue_ijnz__[m - 1];
			    c->cvar_.u[m + (i__ + (j + (k << 4)) * 10) * 5 - 
				    746] = pxi + peta + pzeta - pxi * peta - 
				    peta * pzeta - pzeta * pxi + pxi * peta * 
				    pzeta;
			}
		    }
		}
	    }
	}
    }
    return 0;
} /* setiv_ */

/* Subroutine */ int setcoeff_(context *c)
{
    float d__1;

    c->cgcon_.dxi = 1. / (c->cgcon_.nx0 - 1);
    c->cgcon_.deta = 1. / (c->cgcon_.ny0 - 1);
    c->cgcon_.dzeta = 1. / (c->cgcon_.nz0 - 1);
    c->cgcon_.tx1 = 1. / (c->cgcon_.dxi * c->cgcon_.dxi);
    c->cgcon_.tx2 = 1. / (c->cgcon_.dxi * 2.);
    c->cgcon_.tx3 = 1. / c->cgcon_.dxi;
    c->cgcon_.ty1 = 1. / (c->cgcon_.deta * c->cgcon_.deta);
    c->cgcon_.ty2 = 1. / (c->cgcon_.deta * 2.);
    c->cgcon_.ty3 = 1. / c->cgcon_.deta;
    c->cgcon_.tz1 = 1. / (c->cgcon_.dzeta * c->cgcon_.dzeta);
    c->cgcon_.tz2 = 1. / (c->cgcon_.dzeta * 2.);
    c->cgcon_.tz3 = 1. / c->cgcon_.dzeta;
    c->cgcon_.ii1 = 2;
    c->cgcon_.ii2 = c->cgcon_.nx0 - 1;
    c->cgcon_.ji1 = 2;
    c->cgcon_.ji2 = c->cgcon_.ny0 - 2;
    c->cgcon_.ki1 = 3;
    c->cgcon_.ki2 = c->cgcon_.nz0 - 1;
    c->disp_.dx1 = .75;
    c->disp_.dx2 = c->disp_.dx1;
    c->disp_.dx3 = c->disp_.dx1;
    c->disp_.dx4 = c->disp_.dx1;
    c->disp_.dx5 = c->disp_.dx1;
    c->disp_.dy1 = .75;
    c->disp_.dy2 = c->disp_.dy1;
    c->disp_.dy3 = c->disp_.dy1;
    c->disp_.dy4 = c->disp_.dy1;
    c->disp_.dy5 = c->disp_.dy1;
    c->disp_.dz1 = 1.;
    c->disp_.dz2 = c->disp_.dz1;
    c->disp_.dz3 = c->disp_.dz1;
    c->disp_.dz4 = c->disp_.dz1;
    c->disp_.dz5 = c->disp_.dz1;
    d__1 = max(c->disp_.dx1,c->disp_.dy1);
    c->disp_.dssp = max(d__1,c->disp_.dz1) / 4.;
    c->cexact_.ce[0] = 2.;
    c->cexact_.ce[5] = 0.;
    c->cexact_.ce[10] = 0.;
    c->cexact_.ce[15] = 4.;
    c->cexact_.ce[20] = 5.;
    c->cexact_.ce[25] = 3.;
    c->cexact_.ce[30] = .5;
    c->cexact_.ce[35] = .02;
    c->cexact_.ce[40] = .01;
    c->cexact_.ce[45] = .03;
    c->cexact_.ce[50] = .5;
    c->cexact_.ce[55] = .4;
    c->cexact_.ce[60] = .3;
    c->cexact_.ce[1] = 1.;
    c->cexact_.ce[6] = 0.;
    c->cexact_.ce[11] = 0.;
    c->cexact_.ce[16] = 0.;
    c->cexact_.ce[21] = 1.;
    c->cexact_.ce[26] = 2.;
    c->cexact_.ce[31] = 3.;
    c->cexact_.ce[36] = .01;
    c->cexact_.ce[41] = .03;
    c->cexact_.ce[46] = .02;
    c->cexact_.ce[51] = .4;
    c->cexact_.ce[56] = .3;
    c->cexact_.ce[61] = .5;
    c->cexact_.ce[2] = 2.;
    c->cexact_.ce[7] = 2.;
    c->cexact_.ce[12] = 0.;
    c->cexact_.ce[17] = 0.;
    c->cexact_.ce[22] = 0.;
    c->cexact_.ce[27] = 2.;
    c->cexact_.ce[32] = 3.;
    c->cexact_.ce[37] = .04;
    c->cexact_.ce[42] = .03;
    c->cexact_.ce[47] = .05;
    c->cexact_.ce[52] = .3;
    c->cexact_.ce[57] = .5;
    c->cexact_.ce[62] = .4;
    c->cexact_.ce[3] = 2.;
    c->cexact_.ce[8] = 2.;
    c->cexact_.ce[13] = 0.;
    c->cexact_.ce[18] = 0.;
    c->cexact_.ce[23] = 0.;
    c->cexact_.ce[28] = 2.;
    c->cexact_.ce[33] = 3.;
    c->cexact_.ce[38] = .03;
    c->cexact_.ce[43] = .05;
    c->cexact_.ce[48] = .04;
    c->cexact_.ce[53] = .2;
    c->cexact_.ce[58] = .1;
    c->cexact_.ce[63] = .3;
    c->cexact_.ce[4] = 5.;
    c->cexact_.ce[9] = 4.;
    c->cexact_.ce[14] = 3.;
    c->cexact_.ce[19] = 2.;
    c->cexact_.ce[24] = .1;
    c->cexact_.ce[29] = .4;
    c->cexact_.ce[34] = .3;
    c->cexact_.ce[39] = .05;
    c->cexact_.ce[44] = .04;
    c->cexact_.ce[49] = .03;
    c->cexact_.ce[54] = .1;
    c->cexact_.ce[59] = .3;
    c->cexact_.ce[64] = .2;
    return 0;
} /* setcoeff_ */

/* Subroutine */ int setbv_(context *c)
{
    integer i__1, i__2;

    integer i__, j, k, iglob, jglob;

    integer c__1 = 1;
    i__1 = c->cgcon_.ny;
    for (j = 1; j <= i__1; ++j) {
	jglob = c->cgcon_.jpt + j;
	i__2 = c->cgcon_.nx;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    iglob = c->cgcon_.ipt + i__;
	    exact_(&iglob, &jglob, &c__1, &c->cvar_.u[(i__ + (j + 16) * 10) * 5 
		    - 745], c);
	    exact_(&iglob, &jglob, &c->cgcon_.nz, &c->cvar_.u[(i__ + (j + (
		    c->cgcon_.nz << 4)) * 10) * 5 - 745], c);
	}
    }
    if (c->neigh_.west == -1) {
	i__1 = c->cgcon_.nz;
	for (k = 1; k <= i__1; ++k) {
	    i__2 = c->cgcon_.nx;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		iglob = c->cgcon_.ipt + i__;
		exact_(&iglob, &c__1, &k, &c->cvar_.u[(i__ + ((k << 4) + 1) * 
			10) * 5 - 745], c);
	    }
	}
    }
    if (c->neigh_.east == -1) {
	i__1 = c->cgcon_.nz;
	for (k = 1; k <= i__1; ++k) {
	    i__2 = c->cgcon_.nx;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		iglob = c->cgcon_.ipt + i__;
		exact_(&iglob, &c->cgcon_.ny0, &k, &c->cvar_.u[(i__ + (c->cgcon_.ny 
			+ (k << 4)) * 10) * 5 - 745], c);
	    }
	}
    }
    if (c->neigh_.north == -1) {
	i__1 = c->cgcon_.nz;
	for (k = 1; k <= i__1; ++k) {
	    i__2 = c->cgcon_.ny;
	    for (j = 1; j <= i__2; ++j) {
		jglob = c->cgcon_.jpt + j;
		exact_(&c__1, &jglob, &k, &c->cvar_.u[((j + (k << 4)) * 10 + 1) 
			* 5 - 745], c);
	    }
	}
    }
    if (c->neigh_.south == -1) {
	i__1 = c->cgcon_.nz;
	for (k = 1; k <= i__1; ++k) {
	    i__2 = c->cgcon_.ny;
	    for (j = 1; j <= i__2; ++j) {
		jglob = c->cgcon_.jpt + j;
		exact_(&c->cgcon_.nx0, &jglob, &k, &c->cvar_.u[(c->cgcon_.nx + (j + 
			(k << 4)) * 10) * 5 - 745], c);
	    }
	}
    }
    return 0;
} /* setbv_ */

/* Subroutine */ int rhs_(context *c)
{
    integer i__1, i__2, i__3;
    float d__1, d__2, d__3, d__4, d__5, d__6, d__7, d__8;

    integer i__, j, k, m;
    float q;
    integer l1, l2;
    float u21, u31, u41;
    float u21i, u31i, u41i, u51i;
    integer iex;
    float u21j, u31j, u41j, u51j, tmp, u21k, u31k, u41k, u51k;
    integer ist1, jst1;
    integer iend1, jend1;
    float u21im1, u31im1, u41im1, u51im1, u21jm1, u31jm1, u41jm1, 
	    u51jm1, u21km1, u31km1, u41km1, u51km1;

    if (c->timer_.timeron) {
	timer_start(2, c);
    }
    i__1 = c->cgcon_.nz;
    for (k = 1; k <= i__1; ++k) {
	i__2 = c->cgcon_.ny;
	for (j = 1; j <= i__2; ++j) {
	    i__3 = c->cgcon_.nx;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		for (m = 1; m <= 5; ++m) {
		    c->cvar_.rsd[m + (i__ + (j + (k << 4)) * 10) * 5 - 746] = 
			    -c->cvar_.frct[m + (i__ + (j + (k << 4)) * 10) * 5 
			    - 746];
		}
	    }
	}
    }
    iex = 0;
    if (c->timer_.timeron) {
	timer_start(7, c);
    }
    exchange_3__(c->cvar_.u, &iex, c);
    if (c->timer_.timeron) {
	timer_stop(7, c);
    }
    l1 = 0;
    if (c->neigh_.north == -1) {
	l1 = 1;
    }
    l2 = c->cgcon_.nx + 1;
    if (c->neigh_.south == -1) {
	l2 = c->cgcon_.nx;
    }
    ist1 = 1;
    iend1 = c->cgcon_.nx;
    if (c->neigh_.north == -1) {
	ist1 = 4;
    }
    if (c->neigh_.south == -1) {
	iend1 = c->cgcon_.nx - 3;
    }
    i__1 = c->cgcon_.nz - 1;
    for (k = 2; k <= i__1; ++k) {
	i__2 = c->cgcon_.jend;
	for (j = c->cgcon_.jst; j <= i__2; ++j) {
	    i__3 = l2;
	    for (i__ = l1; i__ <= i__3; ++i__) {
		c->cvar_.flux[(i__ + ((j + k * 14) << 3)) * 5 - 560] = c->cvar_.u[(
			i__ + (j + (k << 4)) * 10) * 5 - 744];
		u21 = c->cvar_.u[(i__ + (j + (k << 4)) * 10) * 5 - 744] / 
			c->cvar_.u[(i__ + (j + (k << 4)) * 10) * 5 - 745];
		q = (c->cvar_.u[(i__ + (j + (k << 4)) * 10) * 5 - 744] * 
			c->cvar_.u[(i__ + (j + (k << 4)) * 10) * 5 - 744] + 
			c->cvar_.u[(i__ + (j + (k << 4)) * 10) * 5 - 743] * 
			c->cvar_.u[(i__ + (j + (k << 4)) * 10) * 5 - 743] + 
			c->cvar_.u[(i__ + (j + (k << 4)) * 10) * 5 - 742] * 
			c->cvar_.u[(i__ + (j + (k << 4)) * 10) * 5 - 742]) * .5 
			/ c->cvar_.u[(i__ + (j + (k << 4)) * 10) * 5 - 745];
		c->cvar_.flux[(i__ + ((j + k * 14) << 3)) * 5 - 559] = c->cvar_.u[(
			i__ + (j + (k << 4)) * 10) * 5 - 744] * u21 + (
			c->cvar_.u[(i__ + (j + (k << 4)) * 10) * 5 - 741] - q) *
			 .4;
		c->cvar_.flux[(i__ + ((j + k * 14) << 3)) * 5 - 558] = c->cvar_.u[(
			i__ + (j + (k << 4)) * 10) * 5 - 743] * u21;
		c->cvar_.flux[(i__ + ((j + k * 14) << 3)) * 5 - 557] = c->cvar_.u[(
			i__ + (j + (k << 4)) * 10) * 5 - 742] * u21;
		c->cvar_.flux[(i__ + ((j + k * 14) << 3)) * 5 - 556] = (c->cvar_.u[(
			i__ + (j + (k << 4)) * 10) * 5 - 741] * 1.4 - q * .4) 
			* u21;
	    }
	    i__3 = c->cgcon_.iend;
	    for (i__ = c->cgcon_.ist; i__ <= i__3; ++i__) {
		for (m = 1; m <= 5; ++m) {
		    c->cvar_.rsd[m + (i__ + (j + (k << 4)) * 10) * 5 - 746] -= 
			    c->cgcon_.tx2 * (c->cvar_.flux[m + (i__ + 1 + ((j + k *
			     14) << 3)) * 5 - 561] - c->cvar_.flux[m + (i__ - 1 
			    + ((j + k * 14) << 3)) * 5 - 561]);
		}
	    }
	    i__3 = l2;
	    for (i__ = c->cgcon_.ist; i__ <= i__3; ++i__) {
		tmp = 1. / c->cvar_.u[(i__ + (j + (k << 4)) * 10) * 5 - 745];
		u21i = tmp * c->cvar_.u[(i__ + (j + (k << 4)) * 10) * 5 - 744];
		u31i = tmp * c->cvar_.u[(i__ + (j + (k << 4)) * 10) * 5 - 743];
		u41i = tmp * c->cvar_.u[(i__ + (j + (k << 4)) * 10) * 5 - 742];
		u51i = tmp * c->cvar_.u[(i__ + (j + (k << 4)) * 10) * 5 - 741];
		tmp = 1. / c->cvar_.u[(i__ - 1 + (j + (k << 4)) * 10) * 5 - 745]
			;
		u21im1 = tmp * c->cvar_.u[(i__ - 1 + (j + (k << 4)) * 10) * 5 - 
			744];
		u31im1 = tmp * c->cvar_.u[(i__ - 1 + (j + (k << 4)) * 10) * 5 - 
			743];
		u41im1 = tmp * c->cvar_.u[(i__ - 1 + (j + (k << 4)) * 10) * 5 - 
			742];
		u51im1 = tmp * c->cvar_.u[(i__ - 1 + (j + (k << 4)) * 10) * 5 - 
			741];
		c->cvar_.flux[(i__ + ((j + k * 14) << 3)) * 5 - 559] = 
			c->cgcon_.tx3 * 1.3333333333333333 * (u21i - u21im1);
		c->cvar_.flux[(i__ + ((j + k * 14) << 3)) * 5 - 558] = 
			c->cgcon_.tx3 * (u31i - u31im1);
		c->cvar_.flux[(i__ + ((j + k * 14) << 3)) * 5 - 557] = 
			c->cgcon_.tx3 * (u41i - u41im1);
		d__1 = u21i;
		d__2 = u31i;
		d__3 = u41i;
		d__4 = u21im1;
		d__5 = u31im1;
		d__6 = u41im1;
		d__7 = u21i;
		d__8 = u21im1;
		c->cvar_.flux[(i__ + ((j + k * 14) << 3)) * 5 - 556] = 
			c->cgcon_.tx3 * -.47999999999999987 * (d__1 * d__1 + 
			d__2 * d__2 + d__3 * d__3 - (d__4 * d__4 + d__5 * 
			d__5 + d__6 * d__6)) + c->cgcon_.tx3 * 
			.16666666666666666 * (d__7 * d__7 - d__8 * d__8) + 
			c->cgcon_.tx3 * 1.9599999999999997 * (u51i - u51im1);
	    }
	    i__3 = c->cgcon_.iend;
	    for (i__ = c->cgcon_.ist; i__ <= i__3; ++i__) {
		c->cvar_.rsd[(i__ + (j + (k << 4)) * 10) * 5 - 745] += 
			c->disp_.dx1 * c->cgcon_.tx1 * (c->cvar_.u[(i__ - 1 + (j + (
			k << 4)) * 10) * 5 - 745] - c->cvar_.u[(i__ + (j + (k <<
			 4)) * 10) * 5 - 745] * 2. + c->cvar_.u[(i__ + 1 + (j + 
			(k << 4)) * 10) * 5 - 745]);
		c->cvar_.rsd[(i__ + (j + (k << 4)) * 10) * 5 - 744] = 
			c->cvar_.rsd[(i__ + (j + (k << 4)) * 10) * 5 - 744] + 
			c->cgcon_.tx3 * .1 * 1. * (c->cvar_.flux[(i__ + 1 + ((j + 
			k * 14) << 3)) * 5 - 559] - c->cvar_.flux[(i__ + ((j + k *
			 14) << 3)) * 5 - 559]) + c->disp_.dx2 * c->cgcon_.tx1 * (
			c->cvar_.u[(i__ - 1 + (j + (k << 4)) * 10) * 5 - 744] - 
			c->cvar_.u[(i__ + (j + (k << 4)) * 10) * 5 - 744] * 2. 
			+ c->cvar_.u[(i__ + 1 + (j + (k << 4)) * 10) * 5 - 744])
			;
		c->cvar_.rsd[(i__ + (j + (k << 4)) * 10) * 5 - 743] = 
			c->cvar_.rsd[(i__ + (j + (k << 4)) * 10) * 5 - 743] + 
			c->cgcon_.tx3 * .1 * 1. * (c->cvar_.flux[(i__ + 1 + ((j + 
			k * 14) << 3)) * 5 - 558] - c->cvar_.flux[(i__ + ((j + k *
			 14) << 3)) * 5 - 558]) + c->disp_.dx3 * c->cgcon_.tx1 * (
			c->cvar_.u[(i__ - 1 + (j + (k << 4)) * 10) * 5 - 743] - 
			c->cvar_.u[(i__ + (j + (k << 4)) * 10) * 5 - 743] * 2. 
			+ c->cvar_.u[(i__ + 1 + (j + (k << 4)) * 10) * 5 - 743])
			;
		c->cvar_.rsd[(i__ + (j + (k << 4)) * 10) * 5 - 742] = 
			c->cvar_.rsd[(i__ + (j + (k << 4)) * 10) * 5 - 742] + 
			c->cgcon_.tx3 * .1 * 1. * (c->cvar_.flux[(i__ + 1 + ((j + 
			k * 14) << 3)) * 5 - 557] - c->cvar_.flux[(i__ + ((j + k *
			 14) << 3)) * 5 - 557]) + c->disp_.dx4 * c->cgcon_.tx1 * (
			c->cvar_.u[(i__ - 1 + (j + (k << 4)) * 10) * 5 - 742] - 
			c->cvar_.u[(i__ + (j + (k << 4)) * 10) * 5 - 742] * 2. 
			+ c->cvar_.u[(i__ + 1 + (j + (k << 4)) * 10) * 5 - 742])
			;
		c->cvar_.rsd[(i__ + (j + (k << 4)) * 10) * 5 - 741] = 
			c->cvar_.rsd[(i__ + (j + (k << 4)) * 10) * 5 - 741] + 
			c->cgcon_.tx3 * .1 * 1. * (c->cvar_.flux[(i__ + 1 + ((j + 
			k * 14) << 3)) * 5 - 556] - c->cvar_.flux[(i__ + ((j + k *
			 14) << 3)) * 5 - 556]) + c->disp_.dx5 * c->cgcon_.tx1 * (
			c->cvar_.u[(i__ - 1 + (j + (k << 4)) * 10) * 5 - 741] - 
			c->cvar_.u[(i__ + (j + (k << 4)) * 10) * 5 - 741] * 2. 
			+ c->cvar_.u[(i__ + 1 + (j + (k << 4)) * 10) * 5 - 741])
			;
	    }
	    if (c->neigh_.north == -1) {
		for (m = 1; m <= 5; ++m) {
		    c->cvar_.rsd[m + ((j + (k << 4)) * 10 + 2) * 5 - 746] -= 
			    c->disp_.dssp * (c->cvar_.u[m + ((j + (k << 4)) * 10 
			    + 2) * 5 - 746] * 5. - c->cvar_.u[m + ((j + (k << 4)
			    ) * 10 + 3) * 5 - 746] * 4. + c->cvar_.u[m + ((j + (
			    k << 4)) * 10 + 4) * 5 - 746]);
		    c->cvar_.rsd[m + ((j + (k << 4)) * 10 + 3) * 5 - 746] -= 
			    c->disp_.dssp * (c->cvar_.u[m + ((j + (k << 4)) * 10 
			    + 2) * 5 - 746] * -4. + c->cvar_.u[m + ((j + (k << 
			    4)) * 10 + 3) * 5 - 746] * 6. - c->cvar_.u[m + ((j 
			    + (k << 4)) * 10 + 4) * 5 - 746] * 4. + c->cvar_.u[
			    m + ((j + (k << 4)) * 10 + 5) * 5 - 746]);
		}
	    }
	    i__3 = iend1;
	    for (i__ = ist1; i__ <= i__3; ++i__) {
		for (m = 1; m <= 5; ++m) {
		    c->cvar_.rsd[m + (i__ + (j + (k << 4)) * 10) * 5 - 746] -= 
			    c->disp_.dssp * (c->cvar_.u[m + (i__ - 2 + (j + (k << 
			    4)) * 10) * 5 - 746] - c->cvar_.u[m + (i__ - 1 + (j 
			    + (k << 4)) * 10) * 5 - 746] * 4. + c->cvar_.u[m + (
			    i__ + (j + (k << 4)) * 10) * 5 - 746] * 6. - 
			    c->cvar_.u[m + (i__ + 1 + (j + (k << 4)) * 10) * 5 
			    - 746] * 4. + c->cvar_.u[m + (i__ + 2 + (j + (k << 
			    4)) * 10) * 5 - 746]);
		}
	    }
	    if (c->neigh_.south == -1) {
		for (m = 1; m <= 5; ++m) {
		    c->cvar_.rsd[m + (c->cgcon_.nx - 2 + (j + (k << 4)) * 10) * 5 
			    - 746] -= c->disp_.dssp * (c->cvar_.u[m + (c->cgcon_.nx 
			    - 4 + (j + (k << 4)) * 10) * 5 - 746] - c->cvar_.u[
			    m + (c->cgcon_.nx - 3 + (j + (k << 4)) * 10) * 5 - 
			    746] * 4. + c->cvar_.u[m + (c->cgcon_.nx - 2 + (j + (
			    k << 4)) * 10) * 5 - 746] * 6. - c->cvar_.u[m + (
			    c->cgcon_.nx - 1 + (j + (k << 4)) * 10) * 5 - 746] *
			     4.);
		    c->cvar_.rsd[m + (c->cgcon_.nx - 1 + (j + (k << 4)) * 10) * 5 
			    - 746] -= c->disp_.dssp * (c->cvar_.u[m + (c->cgcon_.nx 
			    - 3 + (j + (k << 4)) * 10) * 5 - 746] - c->cvar_.u[
			    m + (c->cgcon_.nx - 2 + (j + (k << 4)) * 10) * 5 - 
			    746] * 4. + c->cvar_.u[m + (c->cgcon_.nx - 1 + (j + (
			    k << 4)) * 10) * 5 - 746] * 5.);
		}
	    }
	}
    }
    iex = 1;
    if (c->timer_.timeron) {
	timer_start(7, c);
    }
    exchange_3__(c->cvar_.u, &iex, c);
    if (c->timer_.timeron) {
	timer_stop(7, c);
    }
    l1 = 0;
    if (c->neigh_.west == -1) {
	l1 = 1;
    }
    l2 = c->cgcon_.ny + 1;
    if (c->neigh_.east == -1) {
	l2 = c->cgcon_.ny;
    }
    jst1 = 1;
    jend1 = c->cgcon_.ny;
    if (c->neigh_.west == -1) {
	jst1 = 4;
    }
    if (c->neigh_.east == -1) {
	jend1 = c->cgcon_.ny - 3;
    }
    i__1 = c->cgcon_.nz - 1;
    for (k = 2; k <= i__1; ++k) {
	i__2 = l2;
	for (j = l1; j <= i__2; ++j) {
	    i__3 = c->cgcon_.iend;
	    for (i__ = c->cgcon_.ist; i__ <= i__3; ++i__) {
		c->cvar_.flux[(i__ + ((j + k * 14) << 3)) * 5 - 560] = c->cvar_.u[(
			i__ + (j + (k << 4)) * 10) * 5 - 743];
		u31 = c->cvar_.u[(i__ + (j + (k << 4)) * 10) * 5 - 743] / 
			c->cvar_.u[(i__ + (j + (k << 4)) * 10) * 5 - 745];
		q = (c->cvar_.u[(i__ + (j + (k << 4)) * 10) * 5 - 744] * 
			c->cvar_.u[(i__ + (j + (k << 4)) * 10) * 5 - 744] + 
			c->cvar_.u[(i__ + (j + (k << 4)) * 10) * 5 - 743] * 
			c->cvar_.u[(i__ + (j + (k << 4)) * 10) * 5 - 743] + 
			c->cvar_.u[(i__ + (j + (k << 4)) * 10) * 5 - 742] * 
			c->cvar_.u[(i__ + (j + (k << 4)) * 10) * 5 - 742]) * .5 
			/ c->cvar_.u[(i__ + (j + (k << 4)) * 10) * 5 - 745];
		c->cvar_.flux[(i__ + ((j + k * 14) << 3)) * 5 - 559] = c->cvar_.u[(
			i__ + (j + (k << 4)) * 10) * 5 - 744] * u31;
		c->cvar_.flux[(i__ + ((j + k * 14) << 3)) * 5 - 558] = c->cvar_.u[(
			i__ + (j + (k << 4)) * 10) * 5 - 743] * u31 + (
			c->cvar_.u[(i__ + (j + (k << 4)) * 10) * 5 - 741] - q) *
			 .4;
		c->cvar_.flux[(i__ + ((j + k * 14) << 3)) * 5 - 557] = c->cvar_.u[(
			i__ + (j + (k << 4)) * 10) * 5 - 742] * u31;
		c->cvar_.flux[(i__ + ((j + k * 14) << 3)) * 5 - 556] = (c->cvar_.u[(
			i__ + (j + (k << 4)) * 10) * 5 - 741] * 1.4 - q * .4) 
			* u31;
	    }
	}
	i__2 = c->cgcon_.jend;
	for (j = c->cgcon_.jst; j <= i__2; ++j) {
	    i__3 = c->cgcon_.iend;
	    for (i__ = c->cgcon_.ist; i__ <= i__3; ++i__) {
		for (m = 1; m <= 5; ++m) {
		    c->cvar_.rsd[m + (i__ + (j + (k << 4)) * 10) * 5 - 746] -= 
			    c->cgcon_.ty2 * (c->cvar_.flux[m + (i__ + ((j + 1 + k *
			     14) << 3)) * 5 - 561] - c->cvar_.flux[m + (i__ + ((j 
			    - 1 + k * 14) << 3)) * 5 - 561]);
		}
	    }
	}
	i__2 = l2;
	for (j = c->cgcon_.jst; j <= i__2; ++j) {
	    i__3 = c->cgcon_.iend;
	    for (i__ = c->cgcon_.ist; i__ <= i__3; ++i__) {
		tmp = 1. / c->cvar_.u[(i__ + (j + (k << 4)) * 10) * 5 - 745];
		u21j = tmp * c->cvar_.u[(i__ + (j + (k << 4)) * 10) * 5 - 744];
		u31j = tmp * c->cvar_.u[(i__ + (j + (k << 4)) * 10) * 5 - 743];
		u41j = tmp * c->cvar_.u[(i__ + (j + (k << 4)) * 10) * 5 - 742];
		u51j = tmp * c->cvar_.u[(i__ + (j + (k << 4)) * 10) * 5 - 741];
		tmp = 1. / c->cvar_.u[(i__ + (j - 1 + (k << 4)) * 10) * 5 - 745]
			;
		u21jm1 = tmp * c->cvar_.u[(i__ + (j - 1 + (k << 4)) * 10) * 5 - 
			744];
		u31jm1 = tmp * c->cvar_.u[(i__ + (j - 1 + (k << 4)) * 10) * 5 - 
			743];
		u41jm1 = tmp * c->cvar_.u[(i__ + (j - 1 + (k << 4)) * 10) * 5 - 
			742];
		u51jm1 = tmp * c->cvar_.u[(i__ + (j - 1 + (k << 4)) * 10) * 5 - 
			741];
		c->cvar_.flux[(i__ + ((j + k * 14) << 3)) * 5 - 559] = 
			c->cgcon_.ty3 * (u21j - u21jm1);
		c->cvar_.flux[(i__ + ((j + k * 14) << 3)) * 5 - 558] = 
			c->cgcon_.ty3 * 1.3333333333333333 * (u31j - u31jm1);
		c->cvar_.flux[(i__ + ((j + k * 14) << 3)) * 5 - 557] = 
			c->cgcon_.ty3 * (u41j - u41jm1);
		d__1 = u21j;
		d__2 = u31j;
		d__3 = u41j;
		d__4 = u21jm1;
		d__5 = u31jm1;
		d__6 = u41jm1;
		d__7 = u31j;
		d__8 = u31jm1;
		c->cvar_.flux[(i__ + ((j + k * 14) << 3)) * 5 - 556] = 
			c->cgcon_.ty3 * -.47999999999999987 * (d__1 * d__1 + 
			d__2 * d__2 + d__3 * d__3 - (d__4 * d__4 + d__5 * 
			d__5 + d__6 * d__6)) + c->cgcon_.ty3 * 
			.16666666666666666 * (d__7 * d__7 - d__8 * d__8) + 
			c->cgcon_.ty3 * 1.9599999999999997 * (u51j - u51jm1);
	    }
	}
	i__2 = c->cgcon_.jend;
	for (j = c->cgcon_.jst; j <= i__2; ++j) {
	    i__3 = c->cgcon_.iend;
	    for (i__ = c->cgcon_.ist; i__ <= i__3; ++i__) {
		c->cvar_.rsd[(i__ + (j + (k << 4)) * 10) * 5 - 745] += 
			c->disp_.dy1 * c->cgcon_.ty1 * (c->cvar_.u[(i__ + (j - 1 + (
			k << 4)) * 10) * 5 - 745] - c->cvar_.u[(i__ + (j + (k <<
			 4)) * 10) * 5 - 745] * 2. + c->cvar_.u[(i__ + (j + 1 + 
			(k << 4)) * 10) * 5 - 745]);
		c->cvar_.rsd[(i__ + (j + (k << 4)) * 10) * 5 - 744] = 
			c->cvar_.rsd[(i__ + (j + (k << 4)) * 10) * 5 - 744] + 
			c->cgcon_.ty3 * .1 * 1. * (c->cvar_.flux[(i__ + ((j + 1 + 
			k * 14) << 3)) * 5 - 559] - c->cvar_.flux[(i__ + ((j + k *
			 14) << 3)) * 5 - 559]) + c->disp_.dy2 * c->cgcon_.ty1 * (
			c->cvar_.u[(i__ + (j - 1 + (k << 4)) * 10) * 5 - 744] - 
			c->cvar_.u[(i__ + (j + (k << 4)) * 10) * 5 - 744] * 2. 
			+ c->cvar_.u[(i__ + (j + 1 + (k << 4)) * 10) * 5 - 744])
			;
		c->cvar_.rsd[(i__ + (j + (k << 4)) * 10) * 5 - 743] = 
			c->cvar_.rsd[(i__ + (j + (k << 4)) * 10) * 5 - 743] + 
			c->cgcon_.ty3 * .1 * 1. * (c->cvar_.flux[(i__ + ((j + 1 + 
			k * 14) << 3)) * 5 - 558] - c->cvar_.flux[(i__ + ((j + k *
			 14) << 3)) * 5 - 558]) + c->disp_.dy3 * c->cgcon_.ty1 * (
			c->cvar_.u[(i__ + (j - 1 + (k << 4)) * 10) * 5 - 743] - 
			c->cvar_.u[(i__ + (j + (k << 4)) * 10) * 5 - 743] * 2. 
			+ c->cvar_.u[(i__ + (j + 1 + (k << 4)) * 10) * 5 - 743])
			;
		c->cvar_.rsd[(i__ + (j + (k << 4)) * 10) * 5 - 742] = 
			c->cvar_.rsd[(i__ + (j + (k << 4)) * 10) * 5 - 742] + 
			c->cgcon_.ty3 * .1 * 1. * (c->cvar_.flux[(i__ + ((j + 1 + 
			k * 14) << 3)) * 5 - 557] - c->cvar_.flux[(i__ + ((j + k *
			 14) << 3)) * 5 - 557]) + c->disp_.dy4 * c->cgcon_.ty1 * (
			c->cvar_.u[(i__ + (j - 1 + (k << 4)) * 10) * 5 - 742] - 
			c->cvar_.u[(i__ + (j + (k << 4)) * 10) * 5 - 742] * 2. 
			+ c->cvar_.u[(i__ + (j + 1 + (k << 4)) * 10) * 5 - 742])
			;
		c->cvar_.rsd[(i__ + (j + (k << 4)) * 10) * 5 - 741] = 
			c->cvar_.rsd[(i__ + (j + (k << 4)) * 10) * 5 - 741] + 
			c->cgcon_.ty3 * .1 * 1. * (c->cvar_.flux[(i__ + ((j + 1 + 
			k * 14) << 3)) * 5 - 556] - c->cvar_.flux[(i__ + ((j + k *
			 14) << 3)) * 5 - 556]) + c->disp_.dy5 * c->cgcon_.ty1 * (
			c->cvar_.u[(i__ + (j - 1 + (k << 4)) * 10) * 5 - 741] - 
			c->cvar_.u[(i__ + (j + (k << 4)) * 10) * 5 - 741] * 2. 
			+ c->cvar_.u[(i__ + (j + 1 + (k << 4)) * 10) * 5 - 741])
			;
	    }
	}
	if (c->neigh_.west == -1) {
	    i__2 = c->cgcon_.iend;
	    for (i__ = c->cgcon_.ist; i__ <= i__2; ++i__) {
		for (m = 1; m <= 5; ++m) {
		    c->cvar_.rsd[m + (i__ + ((k << 4) + 2) * 10) * 5 - 746] -= 
			    c->disp_.dssp * (c->cvar_.u[m + (i__ + ((k << 4) + 2) 
			    * 10) * 5 - 746] * 5. - c->cvar_.u[m + (i__ + ((k <<
			     4) + 3) * 10) * 5 - 746] * 4. + c->cvar_.u[m + (
			    i__ + ((k << 4) + 4) * 10) * 5 - 746]);
		    c->cvar_.rsd[m + (i__ + ((k << 4) + 3) * 10) * 5 - 746] -= 
			    c->disp_.dssp * (c->cvar_.u[m + (i__ + ((k << 4) + 2) 
			    * 10) * 5 - 746] * -4. + c->cvar_.u[m + (i__ + ((k 
			    << 4) + 3) * 10) * 5 - 746] * 6. - c->cvar_.u[m + (
			    i__ + ((k << 4) + 4) * 10) * 5 - 746] * 4. + 
			    c->cvar_.u[m + (i__ + ((k << 4) + 5) * 10) * 5 - 
			    746]);
		}
	    }
	}
	i__2 = jend1;
	for (j = jst1; j <= i__2; ++j) {
	    i__3 = c->cgcon_.iend;
	    for (i__ = c->cgcon_.ist; i__ <= i__3; ++i__) {
		for (m = 1; m <= 5; ++m) {
		    c->cvar_.rsd[m + (i__ + (j + (k << 4)) * 10) * 5 - 746] -= 
			    c->disp_.dssp * (c->cvar_.u[m + (i__ + (j - 2 + (k << 
			    4)) * 10) * 5 - 746] - c->cvar_.u[m + (i__ + (j - 1 
			    + (k << 4)) * 10) * 5 - 746] * 4. + c->cvar_.u[m + (
			    i__ + (j + (k << 4)) * 10) * 5 - 746] * 6. - 
			    c->cvar_.u[m + (i__ + (j + 1 + (k << 4)) * 10) * 5 
			    - 746] * 4. + c->cvar_.u[m + (i__ + (j + 2 + (k << 
			    4)) * 10) * 5 - 746]);
		}
	    }
	}
	if (c->neigh_.east == -1) {
	    i__2 = c->cgcon_.iend;
	    for (i__ = c->cgcon_.ist; i__ <= i__2; ++i__) {
		for (m = 1; m <= 5; ++m) {
		    c->cvar_.rsd[m + (i__ + (c->cgcon_.ny - 2 + (k << 4)) * 10) * 
			    5 - 746] -= c->disp_.dssp * (c->cvar_.u[m + (i__ + (
			    c->cgcon_.ny - 4 + (k << 4)) * 10) * 5 - 746] - 
			    c->cvar_.u[m + (i__ + (c->cgcon_.ny - 3 + (k << 4)) * 
			    10) * 5 - 746] * 4. + c->cvar_.u[m + (i__ + (
			    c->cgcon_.ny - 2 + (k << 4)) * 10) * 5 - 746] * 6. 
			    - c->cvar_.u[m + (i__ + (c->cgcon_.ny - 1 + (k << 4)) 
			    * 10) * 5 - 746] * 4.);
		    c->cvar_.rsd[m + (i__ + (c->cgcon_.ny - 1 + (k << 4)) * 10) * 
			    5 - 746] -= c->disp_.dssp * (c->cvar_.u[m + (i__ + (
			    c->cgcon_.ny - 3 + (k << 4)) * 10) * 5 - 746] - 
			    c->cvar_.u[m + (i__ + (c->cgcon_.ny - 2 + (k << 4)) * 
			    10) * 5 - 746] * 4. + c->cvar_.u[m + (i__ + (
			    c->cgcon_.ny - 1 + (k << 4)) * 10) * 5 - 746] * 5.);
		}
	    }
	}
    }
    i__1 = c->cgcon_.nz;
    for (k = 1; k <= i__1; ++k) {
	i__2 = c->cgcon_.jend;
	for (j = c->cgcon_.jst; j <= i__2; ++j) {
	    i__3 = c->cgcon_.iend;
	    for (i__ = c->cgcon_.ist; i__ <= i__3; ++i__) {
		c->cvar_.flux[(i__ + ((j + k * 14) << 3)) * 5 - 560] = c->cvar_.u[(
			i__ + (j + (k << 4)) * 10) * 5 - 742];
		u41 = c->cvar_.u[(i__ + (j + (k << 4)) * 10) * 5 - 742] / 
			c->cvar_.u[(i__ + (j + (k << 4)) * 10) * 5 - 745];
		q = (c->cvar_.u[(i__ + (j + (k << 4)) * 10) * 5 - 744] * 
			c->cvar_.u[(i__ + (j + (k << 4)) * 10) * 5 - 744] + 
			c->cvar_.u[(i__ + (j + (k << 4)) * 10) * 5 - 743] * 
			c->cvar_.u[(i__ + (j + (k << 4)) * 10) * 5 - 743] + 
			c->cvar_.u[(i__ + (j + (k << 4)) * 10) * 5 - 742] * 
			c->cvar_.u[(i__ + (j + (k << 4)) * 10) * 5 - 742]) * .5 
			/ c->cvar_.u[(i__ + (j + (k << 4)) * 10) * 5 - 745];
		c->cvar_.flux[(i__ + ((j + k * 14) << 3)) * 5 - 559] = c->cvar_.u[(
			i__ + (j + (k << 4)) * 10) * 5 - 744] * u41;
		c->cvar_.flux[(i__ + ((j + k * 14) << 3)) * 5 - 558] = c->cvar_.u[(
			i__ + (j + (k << 4)) * 10) * 5 - 743] * u41;
		c->cvar_.flux[(i__ + ((j + k * 14) << 3)) * 5 - 557] = c->cvar_.u[(
			i__ + (j + (k << 4)) * 10) * 5 - 742] * u41 + (
			c->cvar_.u[(i__ + (j + (k << 4)) * 10) * 5 - 741] - q) *
			 .4;
		c->cvar_.flux[(i__ + ((j + k * 14) << 3)) * 5 - 556] = (c->cvar_.u[(
			i__ + (j + (k << 4)) * 10) * 5 - 741] * 1.4 - q * .4) 
			* u41;
	    }
	}
    }
    i__1 = c->cgcon_.nz - 1;
    for (k = 2; k <= i__1; ++k) {
	i__2 = c->cgcon_.jend;
	for (j = c->cgcon_.jst; j <= i__2; ++j) {
	    i__3 = c->cgcon_.iend;
	    for (i__ = c->cgcon_.ist; i__ <= i__3; ++i__) {
		for (m = 1; m <= 5; ++m) {
		    c->cvar_.rsd[m + (i__ + (j + (k << 4)) * 10) * 5 - 746] -= 
			    c->cgcon_.tz2 * (c->cvar_.flux[m + (i__ + ((j + (k + 1)
			     * 14) << 3)) * 5 - 561] - c->cvar_.flux[m + (i__ + (
			    (j + (k - 1) * 14) << 3)) * 5 - 561]);
		}
	    }
	}
    }
    i__1 = c->cgcon_.nz;
    for (k = 2; k <= i__1; ++k) {
	i__2 = c->cgcon_.jend;
	for (j = c->cgcon_.jst; j <= i__2; ++j) {
	    i__3 = c->cgcon_.iend;
	    for (i__ = c->cgcon_.ist; i__ <= i__3; ++i__) {
		tmp = 1. / c->cvar_.u[(i__ + (j + (k << 4)) * 10) * 5 - 745];
		u21k = tmp * c->cvar_.u[(i__ + (j + (k << 4)) * 10) * 5 - 744];
		u31k = tmp * c->cvar_.u[(i__ + (j + (k << 4)) * 10) * 5 - 743];
		u41k = tmp * c->cvar_.u[(i__ + (j + (k << 4)) * 10) * 5 - 742];
		u51k = tmp * c->cvar_.u[(i__ + (j + (k << 4)) * 10) * 5 - 741];
		tmp = 1. / c->cvar_.u[(i__ + (j + ((k - 1) << 4)) * 10) * 5 - 745]
			;
		u21km1 = tmp * c->cvar_.u[(i__ + (j + ((k - 1) << 4)) * 10) * 5 - 
			744];
		u31km1 = tmp * c->cvar_.u[(i__ + (j + ((k - 1) << 4)) * 10) * 5 - 
			743];
		u41km1 = tmp * c->cvar_.u[(i__ + (j + ((k - 1) << 4)) * 10) * 5 - 
			742];
		u51km1 = tmp * c->cvar_.u[(i__ + (j + ((k - 1) << 4)) * 10) * 5 - 
			741];
		c->cvar_.flux[(i__ + ((j + k * 14) << 3)) * 5 - 559] = 
			c->cgcon_.tz3 * (u21k - u21km1);
		c->cvar_.flux[(i__ + ((j + k * 14) << 3)) * 5 - 558] = 
			c->cgcon_.tz3 * (u31k - u31km1);
		c->cvar_.flux[(i__ + ((j + k * 14) << 3)) * 5 - 557] = 
			c->cgcon_.tz3 * 1.3333333333333333 * (u41k - u41km1);
		d__1 = u21k;
		d__2 = u31k;
		d__3 = u41k;
		d__4 = u21km1;
		d__5 = u31km1;
		d__6 = u41km1;
		d__7 = u41k;
		d__8 = u41km1;
		c->cvar_.flux[(i__ + ((j + k * 14) << 3)) * 5 - 556] = 
			c->cgcon_.tz3 * -.47999999999999987 * (d__1 * d__1 + 
			d__2 * d__2 + d__3 * d__3 - (d__4 * d__4 + d__5 * 
			d__5 + d__6 * d__6)) + c->cgcon_.tz3 * 
			.16666666666666666 * (d__7 * d__7 - d__8 * d__8) + 
			c->cgcon_.tz3 * 1.9599999999999997 * (u51k - u51km1);
	    }
	}
    }
    i__1 = c->cgcon_.nz - 1;
    for (k = 2; k <= i__1; ++k) {
	i__2 = c->cgcon_.jend;
	for (j = c->cgcon_.jst; j <= i__2; ++j) {
	    i__3 = c->cgcon_.iend;
	    for (i__ = c->cgcon_.ist; i__ <= i__3; ++i__) {
		c->cvar_.rsd[(i__ + (j + (k << 4)) * 10) * 5 - 745] += 
			c->disp_.dz1 * c->cgcon_.tz1 * (c->cvar_.u[(i__ + (j + ((k - 
			1) << 4)) * 10) * 5 - 745] - c->cvar_.u[(i__ + (j + (k <<
			 4)) * 10) * 5 - 745] * 2. + c->cvar_.u[(i__ + (j + ((k 
			+ 1) << 4)) * 10) * 5 - 745]);
		c->cvar_.rsd[(i__ + (j + (k << 4)) * 10) * 5 - 744] = 
			c->cvar_.rsd[(i__ + (j + (k << 4)) * 10) * 5 - 744] + 
			c->cgcon_.tz3 * .1 * 1. * (c->cvar_.flux[(i__ + 
                              ((j + (k + 1) * 14) << 3)) * 5 - 559] - c->cvar_.flux[(i__ + ((j + k 
			* 14) << 3)) * 5 - 559]) + c->disp_.dz2 * c->cgcon_.tz1 * (
			c->cvar_.u[(i__ + (j + ((k - 1) << 4)) * 10) * 5 - 744] - 
			c->cvar_.u[(i__ + (j + (k << 4)) * 10) * 5 - 744] * 2. 
			+ c->cvar_.u[(i__ + ((j + (k + 1)) << 4) * 10) * 5 - 744])
			;
		c->cvar_.rsd[(i__ + (j + (k << 4)) * 10) * 5 - 743] = 
			c->cvar_.rsd[(i__ + (j + (k << 4)) * 10) * 5 - 743] + 
			c->cgcon_.tz3 * .1 * 1. * (c->cvar_.flux[(i__ + 
                              ((j + (k + 1) * 14) << 3)) * 5 - 558] - c->cvar_.flux[(i__ + 
                              ((j + k * 14) << 3)) * 5 - 558]) + c->disp_.dz3 * c->cgcon_.tz1 * (
			c->cvar_.u[(i__ + (j + ((k - 1) << 4)) * 10) * 5 - 743] - 
			c->cvar_.u[(i__ + (j + (k << 4)) * 10) * 5 - 743] * 2. 
			+ c->cvar_.u[(i__ + (j + ((k + 1) << 4)) * 10) * 5 - 743])
			;
		c->cvar_.rsd[(i__ + (j + (k << 4)) * 10) * 5 - 742] = 
			c->cvar_.rsd[(i__ + (j + (k << 4)) * 10) * 5 - 742] + 
			c->cgcon_.tz3 * .1 * 1. * (c->cvar_.flux[(i__ + 
                              ((j + (k + 1) * 14) << 3)) * 5 - 557] - c->cvar_.flux[(i__ + 
                              ((j + k * 14) << 3)) * 5 - 557]) + c->disp_.dz4 * c->cgcon_.tz1 * (
			c->cvar_.u[(i__ + (j + ((k - 1) << 4)) * 10) * 5 - 742] - 
			c->cvar_.u[(i__ + (j + (k << 4)) * 10) * 5 - 742] * 2. 
			+ c->cvar_.u[(i__ + (j + ((k + 1) << 4)) * 10) * 5 - 742])
			;
		c->cvar_.rsd[(i__ + (j + (k << 4)) * 10) * 5 - 741] = 
			c->cvar_.rsd[(i__ + (j + (k << 4)) * 10) * 5 - 741] + 
			c->cgcon_.tz3 * .1 * 1. * (c->cvar_.flux[(i__ + 
                              ((j + (k + 1) * 14) << 3)) * 5 - 556] - c->cvar_.flux[(i__ + 
                              ((j + k * 14) << 3)) * 5 - 556]) + c->disp_.dz5 * c->cgcon_.tz1 * (
			c->cvar_.u[(i__ + (j + ((k - 1) << 4)) * 10) * 5 - 741] - 
			c->cvar_.u[(i__ + (j + (k << 4)) * 10) * 5 - 741] * 2. 
			+ c->cvar_.u[(i__ + (j + ((k + 1) << 4)) * 10) * 5 - 741])
			;
	    }
	}
    }
    i__1 = c->cgcon_.jend;
    for (j = c->cgcon_.jst; j <= i__1; ++j) {
	i__2 = c->cgcon_.iend;
	for (i__ = c->cgcon_.ist; i__ <= i__2; ++i__) {
	    for (m = 1; m <= 5; ++m) {
		c->cvar_.rsd[m + (i__ + (j + 32) * 10) * 5 - 746] -= 
			c->disp_.dssp * (c->cvar_.u[m + (i__ + (j + 32) * 10) * 5 
			- 746] * 5. - c->cvar_.u[m + (i__ + (j + 48) * 10) * 5 
			- 746] * 4. + c->cvar_.u[m + (i__ + (j + 64) * 10) * 5 
			- 746]);
		c->cvar_.rsd[m + (i__ + (j + 48) * 10) * 5 - 746] -= 
			c->disp_.dssp * (c->cvar_.u[m + (i__ + (j + 32) * 10) * 5 
			- 746] * -4. + c->cvar_.u[m + (i__ + (j + 48) * 10) * 5 
			- 746] * 6. - c->cvar_.u[m + (i__ + (j + 64) * 10) * 5 
			- 746] * 4. + c->cvar_.u[m + (i__ + (j + 80) * 10) * 5 
			- 746]);
	    }
	}
    }
    i__1 = c->cgcon_.nz - 3;
    for (k = 4; k <= i__1; ++k) {
	i__2 = c->cgcon_.jend;
	for (j = c->cgcon_.jst; j <= i__2; ++j) {
	    i__3 = c->cgcon_.iend;
	    for (i__ = c->cgcon_.ist; i__ <= i__3; ++i__) {
		for (m = 1; m <= 5; ++m) {
		    c->cvar_.rsd[m + (i__ + (j + (k << 4)) * 10) * 5 - 746] -= 
			    c->disp_.dssp * (c->cvar_.u[m + (i__ + (j + ((k - 2) << 
			    4)) * 10) * 5 - 746] - c->cvar_.u[m + (i__ + (j + (
			    (k - 1) << 4)) * 10) * 5 - 746] * 4. + c->cvar_.u[m + 
			    (i__ + (j + (k << 4)) * 10) * 5 - 746] * 6. - 
			    c->cvar_.u[m + (i__ + (j + ((k + 1) << 4)) * 10) * 5 
			    - 746] * 4. + c->cvar_.u[m + (i__ + (j + ((k + 2) << 
			    4)) * 10) * 5 - 746]);
		}
	    }
	}
    }
    i__1 = c->cgcon_.jend;
    for (j = c->cgcon_.jst; j <= i__1; ++j) {
	i__2 = c->cgcon_.iend;
	for (i__ = c->cgcon_.ist; i__ <= i__2; ++i__) {
	    for (m = 1; m <= 5; ++m) {
		c->cvar_.rsd[m + (i__ + (j + ((c->cgcon_.nz - 2) << 4)) * 10) * 5 - 
			746] -= c->disp_.dssp * (c->cvar_.u[m + (i__ + (j + (
			(c->cgcon_.nz - 4) << 4)) * 10) * 5 - 746] - c->cvar_.u[m + 
			(i__ + (j + ((c->cgcon_.nz - 3) << 4)) * 10) * 5 - 746] * 
			4. + c->cvar_.u[m + (i__ + (j + ((c->cgcon_.nz - 2) << 4)) *
			 10) * 5 - 746] * 6. - c->cvar_.u[m + (i__ + (j + (
			(c->cgcon_.nz - 1) << 4)) * 10) * 5 - 746] * 4.);
		c->cvar_.rsd[m + (i__ + (j + ((c->cgcon_.nz - 1) << 4)) * 10) * 5 - 
			746] -= c->disp_.dssp * (c->cvar_.u[m + (i__ + (j + (
			(c->cgcon_.nz - 3) << 4)) * 10) * 5 - 746] - c->cvar_.u[m + 
			(i__ + (j + ((c->cgcon_.nz - 2) << 4)) * 10) * 5 - 746] * 
			4. + c->cvar_.u[m + (i__ + (j + ((c->cgcon_.nz - 1) << 4)) *
			 10) * 5 - 746] * 5.);
	    }
	}
    }
    if (c->timer_.timeron) {
	timer_stop(2, c);
    }
    return 0;
} /* rhs_ */

/* Subroutine */ int proc_grid__(context *c)
{

	integer     pow_ii(integer *, integer *);

    integer xdim0, ydim0;

    xdim0 = NNODES_XDIM;
    ydim0 = NNODES_COMPILED / xdim0;
    c->dim_.ydim = (integer) (sqrtf((float) c->dim_.num) + .001);
    c->dim_.xdim = c->dim_.num / c->dim_.ydim;
    while(c->dim_.ydim >= ydim0 && c->dim_.xdim * c->dim_.ydim != c->dim_.num) {
	--c->dim_.ydim;
	c->dim_.xdim = c->dim_.num / c->dim_.ydim;
    }
    if (c->dim_.xdim < xdim0 || c->dim_.ydim < ydim0 || c->dim_.xdim * c->dim_.ydim !=
	     c->dim_.num) {
	if (c->dim_.id == 0) {
	}
	MPI_Abort(MPI_COMM_WORLD, 16);
    }
    if (c->dim_.id == 0 && c->dim_.num != 1 << c->dim_.ndim) {
    }
    c->dim_.row = c->dim_.id % c->dim_.xdim + 1;
    c->dim_.col = c->dim_.id / c->dim_.xdim + 1;
    return 0;
} /* proc_grid__ */

/* Subroutine */ int pintgr_(context *c)
{
    integer i__1, i__2;
    float d__1, d__2, d__3;

    integer i__, j, k;
    integer ind1, ind2;
    float frc1, phi1[(ISIZ2+2)*(ISIZ3+2)]	/* was [14][14] */, phi2[(ISIZ2+2)*(ISIZ3+2)]	/* 
	    was [14][14] */, frc2, frc3;
    integer ibeg, jbeg, ifin, jfin, ifin1, jfin1, iglob, jglob;
    float dummy;
    integer iglob1, iglob2, jglob1, jglob2;

    ibeg = c->cgcon_.nx + 1;
    ifin = 0;
    iglob1 = c->cgcon_.ipt + 1;
    iglob2 = c->cgcon_.ipt + c->cgcon_.nx;
    if (iglob1 >= c->cgcon_.ii1 && iglob2 < c->cgcon_.ii2 + c->cgcon_.nx) {
	ibeg = 1;
    }
    if (iglob1 > c->cgcon_.ii1 - c->cgcon_.nx && iglob2 <= c->cgcon_.ii2) {
	ifin = c->cgcon_.nx;
    }
    if (c->cgcon_.ii1 >= iglob1 && c->cgcon_.ii1 <= iglob2) {
	ibeg = c->cgcon_.ii1 - c->cgcon_.ipt;
    }
    if (c->cgcon_.ii2 >= iglob1 && c->cgcon_.ii2 <= iglob2) {
	ifin = c->cgcon_.ii2 - c->cgcon_.ipt;
    }
    jbeg = c->cgcon_.ny + 1;
    jfin = 0;
    jglob1 = c->cgcon_.jpt + 1;
    jglob2 = c->cgcon_.jpt + c->cgcon_.ny;
    if (jglob1 >= c->cgcon_.ji1 && jglob2 < c->cgcon_.ji2 + c->cgcon_.ny) {
	jbeg = 1;
    }
    if (jglob1 > c->cgcon_.ji1 - c->cgcon_.ny && jglob2 <= c->cgcon_.ji2) {
	jfin = c->cgcon_.ny;
    }
    if (c->cgcon_.ji1 >= jglob1 && c->cgcon_.ji1 <= jglob2) {
	jbeg = c->cgcon_.ji1 - c->cgcon_.jpt;
    }
    if (c->cgcon_.ji2 >= jglob1 && c->cgcon_.ji2 <= jglob2) {
	jfin = c->cgcon_.ji2 - c->cgcon_.jpt;
    }
    ifin1 = ifin;
    jfin1 = jfin;
    if (c->cgcon_.ipt + ifin1 == c->cgcon_.ii2) {
	ifin1 = ifin - 1;
    }
    if (c->cgcon_.jpt + jfin1 == c->cgcon_.ji2) {
	jfin1 = jfin - 1;
    }
    for (i__ = 0; i__ <= 13; ++i__) {
	for (k = 0; k <= 13; ++k) {
	    phi1[i__ + k * 14] = 0.f;
	    phi2[i__ + k * 14] = 0.f;
	}
    }
    i__1 = jfin;
    for (j = jbeg; j <= i__1; ++j) {
	jglob = c->cgcon_.jpt + j;
	i__2 = ifin;
	for (i__ = ibeg; i__ <= i__2; ++i__) {
	    iglob = c->cgcon_.ipt + i__;
	    k = c->cgcon_.ki1;
	    d__1 = c->cvar_.u[(i__ + (j + (k << 4)) * 10) * 5 - 744];
	    d__2 = c->cvar_.u[(i__ + (j + (k << 4)) * 10) * 5 - 743];
	    d__3 = c->cvar_.u[(i__ + (j + (k << 4)) * 10) * 5 - 742];
	    phi1[i__ + j * 14] = (c->cvar_.u[(i__ + (j + (k << 4)) * 10) * 5 - 
		    741] - (d__1 * d__1 + d__2 * d__2 + d__3 * d__3) * .5 / 
		    c->cvar_.u[(i__ + (j + (k << 4)) * 10) * 5 - 745]) * .4;
	    k = c->cgcon_.ki2;
	    d__1 = c->cvar_.u[(i__ + (j + (k << 4)) * 10) * 5 - 744];
	    d__2 = c->cvar_.u[(i__ + (j + (k << 4)) * 10) * 5 - 743];
	    d__3 = c->cvar_.u[(i__ + (j + (k << 4)) * 10) * 5 - 742];
	    phi2[i__ + j * 14] = (c->cvar_.u[(i__ + (j + (k << 4)) * 10) * 5 - 
		    741] - (d__1 * d__1 + d__2 * d__2 + d__3 * d__3) * .5 / 
		    c->cvar_.u[(i__ + (j + (k << 4)) * 10) * 5 - 745]) * .4;
	}
    }
    exchange_4__(phi1, phi2, &ibeg, &ifin1, &jbeg, &jfin1, c);
    frc1 = 0.;
    i__1 = jfin1;
    for (j = jbeg; j <= i__1; ++j) {
	i__2 = ifin1;
	for (i__ = ibeg; i__ <= i__2; ++i__) {
	    frc1 += phi1[i__ + j * 14] + phi1[i__ + 1 + j * 14] + phi1[i__ + (
		    j + 1) * 14] + phi1[i__ + 1 + (j + 1) * 14] + phi2[i__ + 
		    j * 14] + phi2[i__ + 1 + j * 14] + phi2[i__ + (j + 1) * 
		    14] + phi2[i__ + 1 + (j + 1) * 14];
	}
    }
    dummy = frc1;
    c->ftmp[0] = frc1;
    MPI_Allreduce(&dummy, &c->ftmp[0], 1, c->mpistuff_.dp_type__, MPI_SUM, MPI_COMM_WORLD);
    frc1 = c->cgcon_.dxi * c->cgcon_.deta * c->ftmp[0];
    for (i__ = 0; i__ <= 13; ++i__) {
	for (k = 0; k <= 13; ++k) {
	    phi1[i__ + k * 14] = 0.f;
	    phi2[i__ + k * 14] = 0.f;
	}
    }
    jglob = c->cgcon_.jpt + jbeg;
    ind1 = 0;
    if (jglob == c->cgcon_.ji1) {
	ind1 = 1;
	i__1 = c->cgcon_.ki2;
	for (k = c->cgcon_.ki1; k <= i__1; ++k) {
	    i__2 = ifin;
	    for (i__ = ibeg; i__ <= i__2; ++i__) {
		iglob = c->cgcon_.ipt + i__;
		d__1 = c->cvar_.u[(i__ + (jbeg + (k << 4)) * 10) * 5 - 744];
		d__2 = c->cvar_.u[(i__ + (jbeg + (k << 4)) * 10) * 5 - 743];
		d__3 = c->cvar_.u[(i__ + (jbeg + (k << 4)) * 10) * 5 - 742];
		phi1[i__ + k * 14] = (c->cvar_.u[(i__ + (jbeg + (k << 4)) * 10) 
			* 5 - 741] - (d__1 * d__1 + d__2 * d__2 + d__3 * d__3)
			 * .5 / c->cvar_.u[(i__ + (jbeg + (k << 4)) * 10) * 5 - 
			745]) * .4;
	    }
	}
    }
    jglob = c->cgcon_.jpt + jfin;
    ind2 = 0;
    if (jglob == c->cgcon_.ji2) {
	ind2 = 1;
	i__1 = c->cgcon_.ki2;
	for (k = c->cgcon_.ki1; k <= i__1; ++k) {
	    i__2 = ifin;
	    for (i__ = ibeg; i__ <= i__2; ++i__) {
		iglob = c->cgcon_.ipt + i__;
		d__1 = c->cvar_.u[(i__ + (jfin + (k << 4)) * 10) * 5 - 744];
		d__2 = c->cvar_.u[(i__ + (jfin + (k << 4)) * 10) * 5 - 743];
		d__3 = c->cvar_.u[(i__ + (jfin + (k << 4)) * 10) * 5 - 742];
		phi2[i__ + k * 14] = (c->cvar_.u[(i__ + (jfin + (k << 4)) * 10) 
			* 5 - 741] - (d__1 * d__1 + d__2 * d__2 + d__3 * d__3)
			 * .5 / c->cvar_.u[(i__ + (jfin + (k << 4)) * 10) * 5 - 
			745]) * .4;
	    }
	}
    }
    if (ind1 == 1) {
	exchange_5__(phi1, &ibeg, &ifin1, c);
    }
    if (ind2 == 1) {
	exchange_5__(phi2, &ibeg, &ifin1, c);
    }
    frc2 = 0.;
    i__1 = c->cgcon_.ki2 - 1;
    for (k = c->cgcon_.ki1; k <= i__1; ++k) {
	i__2 = ifin1;
	for (i__ = ibeg; i__ <= i__2; ++i__) {
	    frc2 += phi1[i__ + k * 14] + phi1[i__ + 1 + k * 14] + phi1[i__ + (
		    k + 1) * 14] + phi1[i__ + 1 + (k + 1) * 14] + phi2[i__ + 
		    k * 14] + phi2[i__ + 1 + k * 14] + phi2[i__ + (k + 1) * 
		    14] + phi2[i__ + 1 + (k + 1) * 14];
	}
    }
    dummy = frc2;
    c->ftmp[0] = frc2;
    MPI_Allreduce(&dummy, &c->ftmp[0], 1, c->mpistuff_.dp_type__, MPI_SUM, MPI_COMM_WORLD);
    frc2 = c->cgcon_.dxi * c->cgcon_.dzeta * c->ftmp[0];
    for (i__ = 0; i__ <= 13; ++i__) {
	for (k = 0; k <= 13; ++k) {
	    phi1[i__ + k * 14] = 0.f;
	    phi2[i__ + k * 14] = 0.f;
	}
    }
    iglob = c->cgcon_.ipt + ibeg;
    ind1 = 0;
    if (iglob == c->cgcon_.ii1) {
	ind1 = 1;
	i__1 = c->cgcon_.ki2;
	for (k = c->cgcon_.ki1; k <= i__1; ++k) {
	    i__2 = jfin;
	    for (j = jbeg; j <= i__2; ++j) {
		jglob = c->cgcon_.jpt + j;
		d__1 = c->cvar_.u[(ibeg + (j + (k << 4)) * 10) * 5 - 744];
		d__2 = c->cvar_.u[(ibeg + (j + (k << 4)) * 10) * 5 - 743];
		d__3 = c->cvar_.u[(ibeg + (j + (k << 4)) * 10) * 5 - 742];
		phi1[j + k * 14] = (c->cvar_.u[(ibeg + (j + (k << 4)) * 10) * 5 
			- 741] - (d__1 * d__1 + d__2 * d__2 + d__3 * d__3) * 
			.5 / c->cvar_.u[(ibeg + (j + (k << 4)) * 10) * 5 - 745])
			 * .4;
	    }
	}
    }
    iglob = c->cgcon_.ipt + ifin;
    ind2 = 0;
    if (iglob == c->cgcon_.ii2) {
	ind2 = 1;
	i__1 = c->cgcon_.ki2;
	for (k = c->cgcon_.ki1; k <= i__1; ++k) {
	    i__2 = jfin;
	    for (j = jbeg; j <= i__2; ++j) {
		jglob = c->cgcon_.jpt + j;
		d__1 = c->cvar_.u[(ifin + (j + (k << 4)) * 10) * 5 - 744];
		d__2 = c->cvar_.u[(ifin + (j + (k << 4)) * 10) * 5 - 743];
		d__3 = c->cvar_.u[(ifin + (j + (k << 4)) * 10) * 5 - 742];
		phi2[j + k * 14] = (c->cvar_.u[(ifin + (j + (k << 4)) * 10) * 5 
			- 741] - (d__1 * d__1 + d__2 * d__2 + d__3 * d__3) * 
			.5 / c->cvar_.u[(ifin + (j + (k << 4)) * 10) * 5 - 745])
			 * .4;
	    }
	}
    }
    if (ind1 == 1) {
	exchange_6__(phi1, &jbeg, &jfin1, c);
    }
    if (ind2 == 1) {
	exchange_6__(phi2, &jbeg, &jfin1, c);
    }
    frc3 = 0.;
    i__1 = c->cgcon_.ki2 - 1;
    for (k = c->cgcon_.ki1; k <= i__1; ++k) {
	i__2 = jfin1;
	for (j = jbeg; j <= i__2; ++j) {
	    frc3 += phi1[j + k * 14] + phi1[j + 1 + k * 14] + phi1[j + (k + 1)
		     * 14] + phi1[j + 1 + (k + 1) * 14] + phi2[j + k * 14] + 
		    phi2[j + 1 + k * 14] + phi2[j + (k + 1) * 14] + phi2[j + 
		    1 + (k + 1) * 14];
	}
    }
    dummy = frc3;
    c->ftmp[0] = frc3;
    MPI_Allreduce(&dummy, &c->ftmp[0], 1, c->mpistuff_.dp_type__, MPI_SUM, MPI_COMM_WORLD);
    frc3 = c->cgcon_.deta * c->cgcon_.dzeta * c->ftmp[0];
    c->ctscon_.frc = (frc1 + frc2 + frc3) * .25;
    return 0;
} /* pintgr_ */

integer nodedim_(integer *num)
{
    integer ret_val;

    float fnum;

    fnum = (float) (*num);
//    ret_val = (integer) (log(fnum) / log(2.) + 1e-5f);
switch (*num) {
    case 1: ret_val = 0; break;
    case 2: ret_val = 1; break;
    case 4: ret_val = 2; break;
    case 8: ret_val = 3; break;
    case 16: ret_val = 4; break;
    case 32: ret_val = 5; break;
    case 64: ret_val = 6; break;
    case 128: ret_val = 7; break;
    case 256: ret_val = 8; break;
    case 512: ret_val = 9; break;
        default: printf("Unknown num fnum=%f\n", fnum); while(1);
}
    return ret_val;
} /* nodec->dim_ */

/* Subroutine */ int neighbors_(context *c)
{

    c->neigh_.south = -1;
    c->neigh_.east = -1;
    c->neigh_.north = -1;
    c->neigh_.west = -1;
    if (c->dim_.row > 1) {
	c->neigh_.north = c->dim_.id - 1;
    } else {
	c->neigh_.north = -1;
    }
    if (c->dim_.row < c->dim_.xdim) {
	c->neigh_.south = c->dim_.id + 1;
    } else {
	c->neigh_.south = -1;
    }
    if (c->dim_.col > 1) {
	c->neigh_.west = c->dim_.id - c->dim_.xdim;
    } else {
	c->neigh_.west = -1;
    }
    if (c->dim_.col < c->dim_.ydim) {
	c->neigh_.east = c->dim_.id + c->dim_.xdim;
    } else {
	c->neigh_.east = -1;
    }
    return 0;
} /* neighbors_ */

/* Subroutine */ int l2norm_(integer *ldx, integer *ldy, integer *ldz, 
	integer *nx0, integer *ny0, integer *nz0, integer *ist, integer *iend,
	 integer *jst, integer *jend, float *v, float *sum, context *c)
{
    integer v_dim2, v_dim3, v_offset, i__1, i__2, i__3;

    integer i__, j, k, m;
    float dummy[5];

    v_dim2 = *ldx + 2 + 1 + 1;
    v_dim3 = *ldy + 2 + 1 + 1;
    v_offset = 1 + 5 * (-1 + v_dim2 * (-1 + v_dim3));
    v -= v_offset;
    --sum;

    for (m = 1; m <= 5; ++m) {
	dummy[m - 1] = 0.;
    }
    i__1 = *nz0 - 1;
    for (k = 2; k <= i__1; ++k) {
	i__2 = *jend;
	for (j = *jst; j <= i__2; ++j) {
	    i__3 = *iend;
	    for (i__ = *ist; i__ <= i__3; ++i__) {
		for (m = 1; m <= 5; ++m) {
		    dummy[m - 1] += v[m + (i__ + (j + k * v_dim3) * v_dim2) * 
			    5] * v[m + (i__ + (j + k * v_dim3) * v_dim2) * 5];
		}
	    }
	}
    }
    if (c->timer_.timeron) {
	timer_start(10, c);
    }

    int i;
    for (i = 0; i < 5; i++) c->ftmp[i] = sum[1+i];
    MPI_Allreduce(dummy, c->ftmp, 5, c->mpistuff_.dp_type__, MPI_SUM, MPI_COMM_WORLD);
    for (i = 0; i < 5; i++) sum[1+i] = c->ftmp[i];

    if (c->timer_.timeron) {
	timer_stop(10, c);
    }
    for (m = 1; m <= 5; ++m) {
	sum[m] = sqrtf(sum[m] / ((*nx0 - 2) * (*ny0 - 2) * (*nz0 - 2)));
    }
    return 0;
} /* l2norm_ */

/* Subroutine */ int jacu_(integer *k, context *c)
{
    integer i__1, i__2;
    float d__1, d__2, d__3, d__4, d__5, d__6, d__7, d__8, d__9;

    integer i__, j;
    float c34, r43;
    float c1345, tmp1, tmp2, tmp3;

    if (c->timer_.timeron) {
	timer_start(6, c);
    }
    r43 = 1.3333333333333333;
    c1345 = .19599999999999998;
    c34 = .10000000000000001;
    i__1 = c->cgcon_.jend;
    for (j = c->cgcon_.jst; j <= i__1; ++j) {
	i__2 = c->cgcon_.iend;
	for (i__ = c->cgcon_.ist; i__ <= i__2; ++i__) {
	    tmp1 = 1. / c->cvar_.u[(i__ + (j + (*k << 4)) * 10) * 5 - 745];
	    tmp2 = tmp1 * tmp1;
	    tmp3 = tmp1 * tmp2;
	    c->cjac_.d__[((i__ + j * 6) * 5 + 1) * 5 - 180] = c->ctscon_.dt * 2. *
		     (c->cgcon_.tx1 * c->disp_.dx1 + c->cgcon_.ty1 * c->disp_.dy1 + 
		    c->cgcon_.tz1 * c->disp_.dz1) + 1.;
	    c->cjac_.d__[((i__ + j * 6) * 5 + 2) * 5 - 180] = 0.;
	    c->cjac_.d__[((i__ + j * 6) * 5 + 3) * 5 - 180] = 0.;
	    c->cjac_.d__[((i__ + j * 6) * 5 + 4) * 5 - 180] = 0.;
	    c->cjac_.d__[((i__ + j * 6) * 5 + 5) * 5 - 180] = 0.;
	    c->cjac_.d__[((i__ + j * 6) * 5 + 1) * 5 - 179] = c->ctscon_.dt * 2. *
		     (c->cgcon_.tx1 * (-r43 * c34 * tmp2 * c->cvar_.u[(i__ + (j + 
		    (*k << 4)) * 10) * 5 - 744]) + c->cgcon_.ty1 * (-c34 * tmp2 
		    * c->cvar_.u[(i__ + (j + (*k << 4)) * 10) * 5 - 744]) + 
		    c->cgcon_.tz1 * (-c34 * tmp2 * c->cvar_.u[(i__ + (j + (*k << 
		    4)) * 10) * 5 - 744]));
	    c->cjac_.d__[((i__ + j * 6) * 5 + 2) * 5 - 179] = c->ctscon_.dt * 2. *
		     (c->cgcon_.tx1 * r43 * c34 * tmp1 + c->cgcon_.ty1 * c34 * 
		    tmp1 + c->cgcon_.tz1 * c34 * tmp1) + 1. + c->ctscon_.dt * 2. *
		     (c->cgcon_.tx1 * c->disp_.dx2 + c->cgcon_.ty1 * c->disp_.dy2 + 
		    c->cgcon_.tz1 * c->disp_.dz2);
	    c->cjac_.d__[((i__ + j * 6) * 5 + 3) * 5 - 179] = 0.;
	    c->cjac_.d__[((i__ + j * 6) * 5 + 4) * 5 - 179] = 0.;
	    c->cjac_.d__[((i__ + j * 6) * 5 + 5) * 5 - 179] = 0.;
	    c->cjac_.d__[((i__ + j * 6) * 5 + 1) * 5 - 178] = c->ctscon_.dt * 2. *
		     (c->cgcon_.tx1 * (-c34 * tmp2 * c->cvar_.u[(i__ + (j + (*k <<
		     4)) * 10) * 5 - 743]) + c->cgcon_.ty1 * (-r43 * c34 * tmp2 
		    * c->cvar_.u[(i__ + (j + (*k << 4)) * 10) * 5 - 743]) + 
		    c->cgcon_.tz1 * (-c34 * tmp2 * c->cvar_.u[(i__ + (j + (*k << 
		    4)) * 10) * 5 - 743]));
	    c->cjac_.d__[((i__ + j * 6) * 5 + 2) * 5 - 178] = 0.;
	    c->cjac_.d__[((i__ + j * 6) * 5 + 3) * 5 - 178] = c->ctscon_.dt * 2. *
		     (c->cgcon_.tx1 * c34 * tmp1 + c->cgcon_.ty1 * r43 * c34 * 
		    tmp1 + c->cgcon_.tz1 * c34 * tmp1) + 1. + c->ctscon_.dt * 2. *
		     (c->cgcon_.tx1 * c->disp_.dx3 + c->cgcon_.ty1 * c->disp_.dy3 + 
		    c->cgcon_.tz1 * c->disp_.dz3);
	    c->cjac_.d__[((i__ + j * 6) * 5 + 4) * 5 - 178] = 0.;
	    c->cjac_.d__[((i__ + j * 6) * 5 + 5) * 5 - 178] = 0.;
	    c->cjac_.d__[((i__ + j * 6) * 5 + 1) * 5 - 177] = c->ctscon_.dt * 2. *
		     (c->cgcon_.tx1 * (-c34 * tmp2 * c->cvar_.u[(i__ + (j + (*k <<
		     4)) * 10) * 5 - 742]) + c->cgcon_.ty1 * (-c34 * tmp2 * 
		    c->cvar_.u[(i__ + (j + (*k << 4)) * 10) * 5 - 742]) + 
		    c->cgcon_.tz1 * (-r43 * c34 * tmp2 * c->cvar_.u[(i__ + (j + (*
		    k << 4)) * 10) * 5 - 742]));
	    c->cjac_.d__[((i__ + j * 6) * 5 + 2) * 5 - 177] = 0.;
	    c->cjac_.d__[((i__ + j * 6) * 5 + 3) * 5 - 177] = 0.;
	    c->cjac_.d__[((i__ + j * 6) * 5 + 4) * 5 - 177] = c->ctscon_.dt * 2. *
		     (c->cgcon_.tx1 * c34 * tmp1 + c->cgcon_.ty1 * c34 * tmp1 + 
		    c->cgcon_.tz1 * r43 * c34 * tmp1) + 1. + c->ctscon_.dt * 2. * 
		    (c->cgcon_.tx1 * c->disp_.dx4 + c->cgcon_.ty1 * c->disp_.dy4 + 
		    c->cgcon_.tz1 * c->disp_.dz4);
	    c->cjac_.d__[((i__ + j * 6) * 5 + 5) * 5 - 177] = 0.;
	    d__1 = c->cvar_.u[(i__ + (j + (*k << 4)) * 10) * 5 - 744];
	    d__2 = c->cvar_.u[(i__ + (j + (*k << 4)) * 10) * 5 - 743];
	    d__3 = c->cvar_.u[(i__ + (j + (*k << 4)) * 10) * 5 - 742];
	    d__4 = c->cvar_.u[(i__ + (j + (*k << 4)) * 10) * 5 - 744];
	    d__5 = c->cvar_.u[(i__ + (j + (*k << 4)) * 10) * 5 - 743];
	    d__6 = c->cvar_.u[(i__ + (j + (*k << 4)) * 10) * 5 - 742];
	    d__7 = c->cvar_.u[(i__ + (j + (*k << 4)) * 10) * 5 - 744];
	    d__8 = c->cvar_.u[(i__ + (j + (*k << 4)) * 10) * 5 - 743];
	    d__9 = c->cvar_.u[(i__ + (j + (*k << 4)) * 10) * 5 - 742];
	    c->cjac_.d__[((i__ + j * 6) * 5 + 1) * 5 - 176] = c->ctscon_.dt * 2. *
		     (c->cgcon_.tx1 * (-(r43 * c34 - c1345) * tmp3 * (d__1 * 
		    d__1) - (c34 - c1345) * tmp3 * (d__2 * d__2) - (c34 - 
		    c1345) * tmp3 * (d__3 * d__3) - c1345 * tmp2 * c->cvar_.u[(
		    i__ + (j + (*k << 4)) * 10) * 5 - 741]) + c->cgcon_.ty1 * (
		    -(c34 - c1345) * tmp3 * (d__4 * d__4) - (r43 * c34 - 
		    c1345) * tmp3 * (d__5 * d__5) - (c34 - c1345) * tmp3 * (
		    d__6 * d__6) - c1345 * tmp2 * c->cvar_.u[(i__ + (j + (*k << 
		    4)) * 10) * 5 - 741]) + c->cgcon_.tz1 * (-(c34 - c1345) * 
		    tmp3 * (d__7 * d__7) - (c34 - c1345) * tmp3 * (d__8 * 
		    d__8) - (r43 * c34 - c1345) * tmp3 * (d__9 * d__9) - 
		    c1345 * tmp2 * c->cvar_.u[(i__ + (j + (*k << 4)) * 10) * 5 
		    - 741]));
	    c->cjac_.d__[((i__ + j * 6) * 5 + 2) * 5 - 176] = c->ctscon_.dt * 2. *
		     (c->cgcon_.tx1 * (r43 * c34 - c1345) * tmp2 * c->cvar_.u[(
		    i__ + (j + (*k << 4)) * 10) * 5 - 744] + c->cgcon_.ty1 * (
		    c34 - c1345) * tmp2 * c->cvar_.u[(i__ + (j + (*k << 4)) * 
		    10) * 5 - 744] + c->cgcon_.tz1 * (c34 - c1345) * tmp2 * 
		    c->cvar_.u[(i__ + (j + (*k << 4)) * 10) * 5 - 744]);
	    c->cjac_.d__[((i__ + j * 6) * 5 + 3) * 5 - 176] = c->ctscon_.dt * 2. *
		     (c->cgcon_.tx1 * (c34 - c1345) * tmp2 * c->cvar_.u[(i__ + (j 
		    + (*k << 4)) * 10) * 5 - 743] + c->cgcon_.ty1 * (r43 * c34 
		    - c1345) * tmp2 * c->cvar_.u[(i__ + (j + (*k << 4)) * 10) * 
		    5 - 743] + c->cgcon_.tz1 * (c34 - c1345) * tmp2 * c->cvar_.u[(
		    i__ + (j + (*k << 4)) * 10) * 5 - 743]);
	    c->cjac_.d__[((i__ + j * 6) * 5 + 4) * 5 - 176] = c->ctscon_.dt * 2. *
		     (c->cgcon_.tx1 * (c34 - c1345) * tmp2 * c->cvar_.u[(i__ + (j 
		    + (*k << 4)) * 10) * 5 - 742] + c->cgcon_.ty1 * (c34 - 
		    c1345) * tmp2 * c->cvar_.u[(i__ + (j + (*k << 4)) * 10) * 5 
		    - 742] + c->cgcon_.tz1 * (r43 * c34 - c1345) * tmp2 * 
		    c->cvar_.u[(i__ + (j + (*k << 4)) * 10) * 5 - 742]);
	    c->cjac_.d__[((i__ + j * 6) * 5 + 5) * 5 - 176] = c->ctscon_.dt * 2. *
		     (c->cgcon_.tx1 * c1345 * tmp1 + c->cgcon_.ty1 * c1345 * tmp1 
		    + c->cgcon_.tz1 * c1345 * tmp1) + 1. + c->ctscon_.dt * 2. * (
		    c->cgcon_.tx1 * c->disp_.dx5 + c->cgcon_.ty1 * c->disp_.dy5 + 
		    c->cgcon_.tz1 * c->disp_.dz5);
	    tmp1 = 1. / c->cvar_.u[(i__ + 1 + (j + (*k << 4)) * 10) * 5 - 745];
	    tmp2 = tmp1 * tmp1;
	    tmp3 = tmp1 * tmp2;
	    c->cjac_.a[((i__ + j * 6) * 5 + 1) * 5 - 180] = -c->ctscon_.dt * 
		    c->cgcon_.tx1 * c->disp_.dx1;
	    c->cjac_.a[((i__ + j * 6) * 5 + 2) * 5 - 180] = c->ctscon_.dt * 
		    c->cgcon_.tx2;
	    c->cjac_.a[((i__ + j * 6) * 5 + 3) * 5 - 180] = 0.;
	    c->cjac_.a[((i__ + j * 6) * 5 + 4) * 5 - 180] = 0.;
	    c->cjac_.a[((i__ + j * 6) * 5 + 5) * 5 - 180] = 0.;
	    d__1 = c->cvar_.u[(i__ + 1 + (j + (*k << 4)) * 10) * 5 - 744] * 
		    tmp1;
	    c->cjac_.a[((i__ + j * 6) * 5 + 1) * 5 - 179] = c->ctscon_.dt * 
		    c->cgcon_.tx2 * (-(d__1 * d__1) + (c->cvar_.u[(i__ + 1 + (j + 
		    (*k << 4)) * 10) * 5 - 744] * c->cvar_.u[(i__ + 1 + (j + (*
		    k << 4)) * 10) * 5 - 744] + c->cvar_.u[(i__ + 1 + (j + (*k 
		    << 4)) * 10) * 5 - 743] * c->cvar_.u[(i__ + 1 + (j + (*k << 
		    4)) * 10) * 5 - 743] + c->cvar_.u[(i__ + 1 + (j + (*k << 4))
		     * 10) * 5 - 742] * c->cvar_.u[(i__ + 1 + (j + (*k << 4)) * 
		    10) * 5 - 742]) * .20000000000000001 * tmp2) - 
		    c->ctscon_.dt * c->cgcon_.tx1 * (-r43 * c34 * tmp2 * c->cvar_.u[
		    (i__ + 1 + (j + (*k << 4)) * 10) * 5 - 744]);
	    c->cjac_.a[((i__ + j * 6) * 5 + 2) * 5 - 179] = c->ctscon_.dt * 
		    c->cgcon_.tx2 * (c->cvar_.u[(i__ + 1 + (j + (*k << 4)) * 10) *
		     5 - 744] * tmp1 * 1.6000000000000001) - c->ctscon_.dt * 
		    c->cgcon_.tx1 * (r43 * c34 * tmp1) - c->ctscon_.dt * 
		    c->cgcon_.tx1 * c->disp_.dx2;
	    c->cjac_.a[((i__ + j * 6) * 5 + 3) * 5 - 179] = c->ctscon_.dt * 
		    c->cgcon_.tx2 * (c->cvar_.u[(i__ + 1 + (j + (*k << 4)) * 10) *
		     5 - 743] * tmp1 * -.4);
	    c->cjac_.a[((i__ + j * 6) * 5 + 4) * 5 - 179] = c->ctscon_.dt * 
		    c->cgcon_.tx2 * (c->cvar_.u[(i__ + 1 + (j + (*k << 4)) * 10) *
		     5 - 742] * tmp1 * -.4);
	    c->cjac_.a[((i__ + j * 6) * 5 + 5) * 5 - 179] = c->ctscon_.dt * 
		    c->cgcon_.tx2 * .4;
	    c->cjac_.a[((i__ + j * 6) * 5 + 1) * 5 - 178] = c->ctscon_.dt * 
		    c->cgcon_.tx2 * (-(c->cvar_.u[(i__ + 1 + (j + (*k << 4)) * 10)
		     * 5 - 744] * c->cvar_.u[(i__ + 1 + (j + (*k << 4)) * 10) * 
		    5 - 743]) * tmp2) - c->ctscon_.dt * c->cgcon_.tx1 * (-c34 * 
		    tmp2 * c->cvar_.u[(i__ + 1 + (j + (*k << 4)) * 10) * 5 - 
		    743]);
	    c->cjac_.a[((i__ + j * 6) * 5 + 2) * 5 - 178] = c->ctscon_.dt * 
		    c->cgcon_.tx2 * (c->cvar_.u[(i__ + 1 + (j + (*k << 4)) * 10) *
		     5 - 743] * tmp1);
	    c->cjac_.a[((i__ + j * 6) * 5 + 3) * 5 - 178] = c->ctscon_.dt * 
		    c->cgcon_.tx2 * (c->cvar_.u[(i__ + 1 + (j + (*k << 4)) * 10) *
		     5 - 744] * tmp1) - c->ctscon_.dt * c->cgcon_.tx1 * (c34 * 
		    tmp1) - c->ctscon_.dt * c->cgcon_.tx1 * c->disp_.dx3;
	    c->cjac_.a[((i__ + j * 6) * 5 + 4) * 5 - 178] = 0.;
	    c->cjac_.a[((i__ + j * 6) * 5 + 5) * 5 - 178] = 0.;
	    c->cjac_.a[((i__ + j * 6) * 5 + 1) * 5 - 177] = c->ctscon_.dt * 
		    c->cgcon_.tx2 * (-(c->cvar_.u[(i__ + 1 + (j + (*k << 4)) * 10)
		     * 5 - 744] * c->cvar_.u[(i__ + 1 + (j + (*k << 4)) * 10) * 
		    5 - 742]) * tmp2) - c->ctscon_.dt * c->cgcon_.tx1 * (-c34 * 
		    tmp2 * c->cvar_.u[(i__ + 1 + (j + (*k << 4)) * 10) * 5 - 
		    742]);
	    c->cjac_.a[((i__ + j * 6) * 5 + 2) * 5 - 177] = c->ctscon_.dt * 
		    c->cgcon_.tx2 * (c->cvar_.u[(i__ + 1 + (j + (*k << 4)) * 10) *
		     5 - 742] * tmp1);
	    c->cjac_.a[((i__ + j * 6) * 5 + 3) * 5 - 177] = 0.;
	    c->cjac_.a[((i__ + j * 6) * 5 + 4) * 5 - 177] = c->ctscon_.dt * 
		    c->cgcon_.tx2 * (c->cvar_.u[(i__ + 1 + (j + (*k << 4)) * 10) *
		     5 - 744] * tmp1) - c->ctscon_.dt * c->cgcon_.tx1 * (c34 * 
		    tmp1) - c->ctscon_.dt * c->cgcon_.tx1 * c->disp_.dx4;
	    c->cjac_.a[((i__ + j * 6) * 5 + 5) * 5 - 177] = 0.;
	    d__1 = c->cvar_.u[(i__ + 1 + (j + (*k << 4)) * 10) * 5 - 744];
	    d__2 = c->cvar_.u[(i__ + 1 + (j + (*k << 4)) * 10) * 5 - 743];
	    d__3 = c->cvar_.u[(i__ + 1 + (j + (*k << 4)) * 10) * 5 - 742];
	    c->cjac_.a[((i__ + j * 6) * 5 + 1) * 5 - 176] = c->ctscon_.dt * 
		    c->cgcon_.tx2 * (((c->cvar_.u[(i__ + 1 + (j + (*k << 4)) * 10)
		     * 5 - 744] * c->cvar_.u[(i__ + 1 + (j + (*k << 4)) * 10) * 
		    5 - 744] + c->cvar_.u[(i__ + 1 + (j + (*k << 4)) * 10) * 5 
		    - 743] * c->cvar_.u[(i__ + 1 + (j + (*k << 4)) * 10) * 5 - 
		    743] + c->cvar_.u[(i__ + 1 + (j + (*k << 4)) * 10) * 5 - 
		    742] * c->cvar_.u[(i__ + 1 + (j + (*k << 4)) * 10) * 5 - 
		    742]) * .4 * tmp2 - c->cvar_.u[(i__ + 1 + (j + (*k << 4)) * 
		    10) * 5 - 741] * tmp1 * 1.4) * (c->cvar_.u[(i__ + 1 + (j + (
		    *k << 4)) * 10) * 5 - 744] * tmp1)) - c->ctscon_.dt * 
		    c->cgcon_.tx1 * (-(r43 * c34 - c1345) * tmp3 * (d__1 * d__1)
		     - (c34 - c1345) * tmp3 * (d__2 * d__2) - (c34 - c1345) * 
		    tmp3 * (d__3 * d__3) - c1345 * tmp2 * c->cvar_.u[(i__ + 1 + 
		    (j + (*k << 4)) * 10) * 5 - 741]);
	    c->cjac_.a[((i__ + j * 6) * 5 + 2) * 5 - 176] = c->ctscon_.dt * 
		    c->cgcon_.tx2 * (c->cvar_.u[(i__ + 1 + (j + (*k << 4)) * 10) *
		     5 - 741] * tmp1 * 1.4 - (c->cvar_.u[(i__ + 1 + (j + (*k << 
		    4)) * 10) * 5 - 744] * 3. * c->cvar_.u[(i__ + 1 + (j + (*k 
		    << 4)) * 10) * 5 - 744] + c->cvar_.u[(i__ + 1 + (j + (*k << 
		    4)) * 10) * 5 - 743] * c->cvar_.u[(i__ + 1 + (j + (*k << 4))
		     * 10) * 5 - 743] + c->cvar_.u[(i__ + 1 + (j + (*k << 4)) * 
		    10) * 5 - 742] * c->cvar_.u[(i__ + 1 + (j + (*k << 4)) * 10)
		     * 5 - 742]) * tmp2 * .20000000000000001) - c->ctscon_.dt * 
		    c->cgcon_.tx1 * (r43 * c34 - c1345) * tmp2 * c->cvar_.u[(i__ 
		    + 1 + (j + (*k << 4)) * 10) * 5 - 744];
	    c->cjac_.a[((i__ + j * 6) * 5 + 3) * 5 - 176] = c->ctscon_.dt * 
		    c->cgcon_.tx2 * (c->cvar_.u[(i__ + 1 + (j + (*k << 4)) * 10) *
		     5 - 743] * c->cvar_.u[(i__ + 1 + (j + (*k << 4)) * 10) * 5 
		    - 744] * -.4 * tmp2) - c->ctscon_.dt * c->cgcon_.tx1 * (c34 - 
		    c1345) * tmp2 * c->cvar_.u[(i__ + 1 + (j + (*k << 4)) * 10) 
		    * 5 - 743];
	    c->cjac_.a[((i__ + j * 6) * 5 + 4) * 5 - 176] = c->ctscon_.dt * 
		    c->cgcon_.tx2 * (c->cvar_.u[(i__ + 1 + (j + (*k << 4)) * 10) *
		     5 - 742] * c->cvar_.u[(i__ + 1 + (j + (*k << 4)) * 10) * 5 
		    - 744] * -.4 * tmp2) - c->ctscon_.dt * c->cgcon_.tx1 * (c34 - 
		    c1345) * tmp2 * c->cvar_.u[(i__ + 1 + (j + (*k << 4)) * 10) 
		    * 5 - 742];
	    c->cjac_.a[((i__ + j * 6) * 5 + 5) * 5 - 176] = c->ctscon_.dt * 
		    c->cgcon_.tx2 * (c->cvar_.u[(i__ + 1 + (j + (*k << 4)) * 10) *
		     5 - 744] * tmp1 * 1.4) - c->ctscon_.dt * c->cgcon_.tx1 * 
		    c1345 * tmp1 - c->ctscon_.dt * c->cgcon_.tx1 * c->disp_.dx5;
	    tmp1 = 1. / c->cvar_.u[(i__ + (j + 1 + (*k << 4)) * 10) * 5 - 745];
	    tmp2 = tmp1 * tmp1;
	    tmp3 = tmp1 * tmp2;
	    c->cjac_.b[((i__ + j * 6) * 5 + 1) * 5 - 180] = -c->ctscon_.dt * 
		    c->cgcon_.ty1 * c->disp_.dy1;
	    c->cjac_.b[((i__ + j * 6) * 5 + 2) * 5 - 180] = 0.;
	    c->cjac_.b[((i__ + j * 6) * 5 + 3) * 5 - 180] = c->ctscon_.dt * 
		    c->cgcon_.ty2;
	    c->cjac_.b[((i__ + j * 6) * 5 + 4) * 5 - 180] = 0.;
	    c->cjac_.b[((i__ + j * 6) * 5 + 5) * 5 - 180] = 0.;
	    c->cjac_.b[((i__ + j * 6) * 5 + 1) * 5 - 179] = c->ctscon_.dt * 
		    c->cgcon_.ty2 * (-(c->cvar_.u[(i__ + (j + 1 + (*k << 4)) * 10)
		     * 5 - 744] * c->cvar_.u[(i__ + (j + 1 + (*k << 4)) * 10) * 
		    5 - 743]) * tmp2) - c->ctscon_.dt * c->cgcon_.ty1 * (-c34 * 
		    tmp2 * c->cvar_.u[(i__ + (j + 1 + (*k << 4)) * 10) * 5 - 
		    744]);
	    c->cjac_.b[((i__ + j * 6) * 5 + 2) * 5 - 179] = c->ctscon_.dt * 
		    c->cgcon_.ty2 * (c->cvar_.u[(i__ + (j + 1 + (*k << 4)) * 10) *
		     5 - 743] * tmp1) - c->ctscon_.dt * c->cgcon_.ty1 * (c34 * 
		    tmp1) - c->ctscon_.dt * c->cgcon_.ty1 * c->disp_.dy2;
	    c->cjac_.b[((i__ + j * 6) * 5 + 3) * 5 - 179] = c->ctscon_.dt * 
		    c->cgcon_.ty2 * (c->cvar_.u[(i__ + (j + 1 + (*k << 4)) * 10) *
		     5 - 744] * tmp1);
	    c->cjac_.b[((i__ + j * 6) * 5 + 4) * 5 - 179] = 0.;
	    c->cjac_.b[((i__ + j * 6) * 5 + 5) * 5 - 179] = 0.;
	    d__1 = c->cvar_.u[(i__ + (j + 1 + (*k << 4)) * 10) * 5 - 743] * 
		    tmp1;
	    c->cjac_.b[((i__ + j * 6) * 5 + 1) * 5 - 178] = c->ctscon_.dt * 
		    c->cgcon_.ty2 * (-(d__1 * d__1) + (c->cvar_.u[(i__ + (j + 1 + 
		    (*k << 4)) * 10) * 5 - 744] * c->cvar_.u[(i__ + (j + 1 + (*
		    k << 4)) * 10) * 5 - 744] + c->cvar_.u[(i__ + (j + 1 + (*k 
		    << 4)) * 10) * 5 - 743] * c->cvar_.u[(i__ + (j + 1 + (*k << 
		    4)) * 10) * 5 - 743] + c->cvar_.u[(i__ + (j + 1 + (*k << 4))
		     * 10) * 5 - 742] * c->cvar_.u[(i__ + (j + 1 + (*k << 4)) * 
		    10) * 5 - 742]) * tmp2 * .20000000000000001) - 
		    c->ctscon_.dt * c->cgcon_.ty1 * (-r43 * c34 * tmp2 * c->cvar_.u[
		    (i__ + (j + 1 + (*k << 4)) * 10) * 5 - 743]);
	    c->cjac_.b[((i__ + j * 6) * 5 + 2) * 5 - 178] = c->ctscon_.dt * 
		    c->cgcon_.ty2 * (c->cvar_.u[(i__ + (j + 1 + (*k << 4)) * 10) *
		     5 - 744] * tmp1 * -.4);
	    c->cjac_.b[((i__ + j * 6) * 5 + 3) * 5 - 178] = c->ctscon_.dt * 
		    c->cgcon_.ty2 * (c->cvar_.u[(i__ + (j + 1 + (*k << 4)) * 10) *
		     5 - 743] * tmp1 * 1.6000000000000001) - c->ctscon_.dt * 
		    c->cgcon_.ty1 * (r43 * c34 * tmp1) - c->ctscon_.dt * 
		    c->cgcon_.ty1 * c->disp_.dy3;
	    c->cjac_.b[((i__ + j * 6) * 5 + 4) * 5 - 178] = c->ctscon_.dt * 
		    c->cgcon_.ty2 * (c->cvar_.u[(i__ + (j + 1 + (*k << 4)) * 10) *
		     5 - 742] * tmp1 * -.4);
	    c->cjac_.b[((i__ + j * 6) * 5 + 5) * 5 - 178] = c->ctscon_.dt * 
		    c->cgcon_.ty2 * .4;
	    c->cjac_.b[((i__ + j * 6) * 5 + 1) * 5 - 177] = c->ctscon_.dt * 
		    c->cgcon_.ty2 * (-(c->cvar_.u[(i__ + (j + 1 + (*k << 4)) * 10)
		     * 5 - 743] * c->cvar_.u[(i__ + (j + 1 + (*k << 4)) * 10) * 
		    5 - 742]) * tmp2) - c->ctscon_.dt * c->cgcon_.ty1 * (-c34 * 
		    tmp2 * c->cvar_.u[(i__ + (j + 1 + (*k << 4)) * 10) * 5 - 
		    742]);
	    c->cjac_.b[((i__ + j * 6) * 5 + 2) * 5 - 177] = 0.;
	    c->cjac_.b[((i__ + j * 6) * 5 + 3) * 5 - 177] = c->ctscon_.dt * 
		    c->cgcon_.ty2 * (c->cvar_.u[(i__ + (j + 1 + (*k << 4)) * 10) *
		     5 - 742] * tmp1);
	    c->cjac_.b[((i__ + j * 6) * 5 + 4) * 5 - 177] = c->ctscon_.dt * 
		    c->cgcon_.ty2 * (c->cvar_.u[(i__ + (j + 1 + (*k << 4)) * 10) *
		     5 - 743] * tmp1) - c->ctscon_.dt * c->cgcon_.ty1 * (c34 * 
		    tmp1) - c->ctscon_.dt * c->cgcon_.ty1 * c->disp_.dy4;
	    c->cjac_.b[((i__ + j * 6) * 5 + 5) * 5 - 177] = 0.;
	    d__1 = c->cvar_.u[(i__ + (j + 1 + (*k << 4)) * 10) * 5 - 744];
	    d__2 = c->cvar_.u[(i__ + (j + 1 + (*k << 4)) * 10) * 5 - 743];
	    d__3 = c->cvar_.u[(i__ + (j + 1 + (*k << 4)) * 10) * 5 - 742];
	    c->cjac_.b[((i__ + j * 6) * 5 + 1) * 5 - 176] = c->ctscon_.dt * 
		    c->cgcon_.ty2 * (((c->cvar_.u[(i__ + (j + 1 + (*k << 4)) * 10)
		     * 5 - 744] * c->cvar_.u[(i__ + (j + 1 + (*k << 4)) * 10) * 
		    5 - 744] + c->cvar_.u[(i__ + (j + 1 + (*k << 4)) * 10) * 5 
		    - 743] * c->cvar_.u[(i__ + (j + 1 + (*k << 4)) * 10) * 5 - 
		    743] + c->cvar_.u[(i__ + (j + 1 + (*k << 4)) * 10) * 5 - 
		    742] * c->cvar_.u[(i__ + (j + 1 + (*k << 4)) * 10) * 5 - 
		    742]) * .4 * tmp2 - c->cvar_.u[(i__ + (j + 1 + (*k << 4)) * 
		    10) * 5 - 741] * tmp1 * 1.4) * (c->cvar_.u[(i__ + (j + 1 + (
		    *k << 4)) * 10) * 5 - 743] * tmp1)) - c->ctscon_.dt * 
		    c->cgcon_.ty1 * (-(c34 - c1345) * tmp3 * (d__1 * d__1) - (
		    r43 * c34 - c1345) * tmp3 * (d__2 * d__2) - (c34 - c1345) 
		    * tmp3 * (d__3 * d__3) - c1345 * tmp2 * c->cvar_.u[(i__ + (
		    j + 1 + (*k << 4)) * 10) * 5 - 741]);
	    c->cjac_.b[((i__ + j * 6) * 5 + 2) * 5 - 176] = c->ctscon_.dt * 
		    c->cgcon_.ty2 * (c->cvar_.u[(i__ + (j + 1 + (*k << 4)) * 10) *
		     5 - 744] * c->cvar_.u[(i__ + (j + 1 + (*k << 4)) * 10) * 5 
		    - 743] * -.4 * tmp2) - c->ctscon_.dt * c->cgcon_.ty1 * (c34 - 
		    c1345) * tmp2 * c->cvar_.u[(i__ + (j + 1 + (*k << 4)) * 10) 
		    * 5 - 744];
	    c->cjac_.b[((i__ + j * 6) * 5 + 3) * 5 - 176] = c->ctscon_.dt * 
		    c->cgcon_.ty2 * (c->cvar_.u[(i__ + (j + 1 + (*k << 4)) * 10) *
		     5 - 741] * tmp1 * 1.4 - (c->cvar_.u[(i__ + (j + 1 + (*k << 
		    4)) * 10) * 5 - 744] * c->cvar_.u[(i__ + (j + 1 + (*k << 4))
		     * 10) * 5 - 744] + c->cvar_.u[(i__ + (j + 1 + (*k << 4)) * 
		    10) * 5 - 743] * 3. * c->cvar_.u[(i__ + (j + 1 + (*k << 4)) 
		    * 10) * 5 - 743] + c->cvar_.u[(i__ + (j + 1 + (*k << 4)) * 
		    10) * 5 - 742] * c->cvar_.u[(i__ + (j + 1 + (*k << 4)) * 10)
		     * 5 - 742]) * tmp2 * .20000000000000001) - c->ctscon_.dt * 
		    c->cgcon_.ty1 * (r43 * c34 - c1345) * tmp2 * c->cvar_.u[(i__ 
		    + (j + 1 + (*k << 4)) * 10) * 5 - 743];
	    c->cjac_.b[((i__ + j * 6) * 5 + 4) * 5 - 176] = c->ctscon_.dt * 
		    c->cgcon_.ty2 * (c->cvar_.u[(i__ + (j + 1 + (*k << 4)) * 10) *
		     5 - 743] * c->cvar_.u[(i__ + (j + 1 + (*k << 4)) * 10) * 5 
		    - 742] * -.4 * tmp2) - c->ctscon_.dt * c->cgcon_.ty1 * (c34 - 
		    c1345) * tmp2 * c->cvar_.u[(i__ + (j + 1 + (*k << 4)) * 10) 
		    * 5 - 742];
	    c->cjac_.b[((i__ + j * 6) * 5 + 5) * 5 - 176] = c->ctscon_.dt * 
		    c->cgcon_.ty2 * (c->cvar_.u[(i__ + (j + 1 + (*k << 4)) * 10) *
		     5 - 743] * tmp1 * 1.4) - c->ctscon_.dt * c->cgcon_.ty1 * 
		    c1345 * tmp1 - c->ctscon_.dt * c->cgcon_.ty1 * c->disp_.dy5;
	    tmp1 = 1. / c->cvar_.u[(i__ + (j + ((*k + 1) << 4)) * 10) * 5 - 745];
	    tmp2 = tmp1 * tmp1;
	    tmp3 = tmp1 * tmp2;
	    c->cjac_.c__[((i__ + j * 6) * 5 + 1) * 5 - 180] = -c->ctscon_.dt * 
		    c->cgcon_.tz1 * c->disp_.dz1;
	    c->cjac_.c__[((i__ + j * 6) * 5 + 2) * 5 - 180] = 0.;
	    c->cjac_.c__[((i__ + j * 6) * 5 + 3) * 5 - 180] = 0.;
	    c->cjac_.c__[((i__ + j * 6) * 5 + 4) * 5 - 180] = c->ctscon_.dt * 
		    c->cgcon_.tz2;
	    c->cjac_.c__[((i__ + j * 6) * 5 + 5) * 5 - 180] = 0.;
	    c->cjac_.c__[((i__ + j * 6) * 5 + 1) * 5 - 179] = c->ctscon_.dt * 
		    c->cgcon_.tz2 * (-(c->cvar_.u[(i__ + (j + ((*k + 1) << 4)) * 10)
		     * 5 - 744] * c->cvar_.u[(i__ + (j + ((*k + 1) << 4)) * 10) * 
		    5 - 742]) * tmp2) - c->ctscon_.dt * c->cgcon_.tz1 * (-c34 * 
		    tmp2 * c->cvar_.u[(i__ + (j + ((*k + 1) << 4)) * 10) * 5 - 
		    744]);
	    c->cjac_.c__[((i__ + j * 6) * 5 + 2) * 5 - 179] = c->ctscon_.dt * 
		    c->cgcon_.tz2 * (c->cvar_.u[(i__ + (j + ((*k + 1) << 4)) * 10) *
		     5 - 742] * tmp1) - c->ctscon_.dt * c->cgcon_.tz1 * c34 * 
		    tmp1 - c->ctscon_.dt * c->cgcon_.tz1 * c->disp_.dz2;
	    c->cjac_.c__[((i__ + j * 6) * 5 + 3) * 5 - 179] = 0.;
	    c->cjac_.c__[((i__ + j * 6) * 5 + 4) * 5 - 179] = c->ctscon_.dt * 
		    c->cgcon_.tz2 * (c->cvar_.u[(i__ + (j + ((*k + 1) << 4)) * 10) *
		     5 - 744] * tmp1);
	    c->cjac_.c__[((i__ + j * 6) * 5 + 5) * 5 - 179] = 0.;
	    c->cjac_.c__[((i__ + j * 6) * 5 + 1) * 5 - 178] = c->ctscon_.dt * 
		    c->cgcon_.tz2 * (-(c->cvar_.u[(i__ + (j + ((*k + 1) << 4)) * 10)
		     * 5 - 743] * c->cvar_.u[(i__ + (j + ((*k + 1) << 4)) * 10) * 
		    5 - 742]) * tmp2) - c->ctscon_.dt * c->cgcon_.tz1 * (-c34 * 
		    tmp2 * c->cvar_.u[(i__ + (j + ((*k + 1) << 4)) * 10) * 5 - 
		    743]);
	    c->cjac_.c__[((i__ + j * 6) * 5 + 2) * 5 - 178] = 0.;
	    c->cjac_.c__[((i__ + j * 6) * 5 + 3) * 5 - 178] = c->ctscon_.dt * 
		    c->cgcon_.tz2 * (c->cvar_.u[(i__ + (j + ((*k + 1) << 4)) * 10) *
		     5 - 742] * tmp1) - c->ctscon_.dt * c->cgcon_.tz1 * (c34 * 
		    tmp1) - c->ctscon_.dt * c->cgcon_.tz1 * c->disp_.dz3;
	    c->cjac_.c__[((i__ + j * 6) * 5 + 4) * 5 - 178] = c->ctscon_.dt * 
		    c->cgcon_.tz2 * (c->cvar_.u[(i__ + (j + ((*k + 1) << 4)) * 10) *
		     5 - 743] * tmp1);
	    c->cjac_.c__[((i__ + j * 6) * 5 + 5) * 5 - 178] = 0.;
	    d__1 = c->cvar_.u[(i__ + (j + ((*k + 1) << 4)) * 10) * 5 - 742] * 
		    tmp1;
	    c->cjac_.c__[((i__ + j * 6) * 5 + 1) * 5 - 177] = c->ctscon_.dt * 
		    c->cgcon_.tz2 * (-(d__1 * d__1) + (c->cvar_.u[(i__ + (j + ((*k 
		    + 1) << 4)) * 10) * 5 - 744] * c->cvar_.u[(i__ + (j + ((*k + 
		    1) << 4)) * 10) * 5 - 744] + c->cvar_.u[(i__ + (j + ((*k + 1)
		    << 4)) * 10) * 5 - 743] * c->cvar_.u[(i__ + (j + ((*k + 1) << 
		    4)) * 10) * 5 - 743] + c->cvar_.u[(i__ + (j + ((*k + 1) << 4))
		     * 10) * 5 - 742] * c->cvar_.u[(i__ + (j + ((*k + 1) << 4)) * 
		    10) * 5 - 742]) * tmp2 * .20000000000000001) - 
		    c->ctscon_.dt * c->cgcon_.tz1 * (-r43 * c34 * tmp2 * c->cvar_.u[
		    (i__ + (j + ((*k + 1) << 4)) * 10) * 5 - 742]);
	    c->cjac_.c__[((i__ + j * 6) * 5 + 2) * 5 - 177] = c->ctscon_.dt * 
		    c->cgcon_.tz2 * (c->cvar_.u[(i__ + (j + ((*k + 1) << 4)) * 10) *
		     5 - 744] * tmp1 * -.4);
	    c->cjac_.c__[((i__ + j * 6) * 5 + 3) * 5 - 177] = c->ctscon_.dt * 
		    c->cgcon_.tz2 * (c->cvar_.u[(i__ + (j + ((*k + 1) << 4)) * 10) *
		     5 - 743] * tmp1 * -.4);
	    c->cjac_.c__[((i__ + j * 6) * 5 + 4) * 5 - 177] = c->ctscon_.dt * 
		    c->cgcon_.tz2 * 1.6000000000000001 * (c->cvar_.u[(i__ + (j + (
		    (*k + 1) << 4)) * 10) * 5 - 742] * tmp1) - c->ctscon_.dt * 
		    c->cgcon_.tz1 * (r43 * c34 * tmp1) - c->ctscon_.dt * 
		    c->cgcon_.tz1 * c->disp_.dz4;
	    c->cjac_.c__[((i__ + j * 6) * 5 + 5) * 5 - 177] = c->ctscon_.dt * 
		    c->cgcon_.tz2 * .4;
	    d__1 = c->cvar_.u[(i__ + (j + ((*k + 1) << 4)) * 10) * 5 - 744];
	    d__2 = c->cvar_.u[(i__ + (j + ((*k + 1) << 4)) * 10) * 5 - 743];
	    d__3 = c->cvar_.u[(i__ + (j + ((*k + 1) << 4)) * 10) * 5 - 742];
	    c->cjac_.c__[((i__ + j * 6) * 5 + 1) * 5 - 176] = c->ctscon_.dt * 
		    c->cgcon_.tz2 * (((c->cvar_.u[(i__ + (j + ((*k + 1) << 4)) * 10)
		     * 5 - 744] * c->cvar_.u[(i__ + (j + ((*k + 1) << 4)) * 10) * 
		    5 - 744] + c->cvar_.u[(i__ + (j + ((*k + 1) << 4)) * 10) * 5 
		    - 743] * c->cvar_.u[(i__ + (j + ((*k + 1) << 4)) * 10) * 5 - 
		    743] + c->cvar_.u[(i__ + (j + ((*k + 1) << 4)) * 10) * 5 - 
		    742] * c->cvar_.u[(i__ + (j + ((*k + 1) << 4)) * 10) * 5 - 
		    742]) * .4 * tmp2 - c->cvar_.u[(i__ + (j + ((*k + 1) << 4)) * 
		    10) * 5 - 741] * tmp1 * 1.4) * (c->cvar_.u[(i__ + (j + ((*k 
		    + 1) << 4)) * 10) * 5 - 742] * tmp1)) - c->ctscon_.dt * 
		    c->cgcon_.tz1 * (-(c34 - c1345) * tmp3 * (d__1 * d__1) - (
		    c34 - c1345) * tmp3 * (d__2 * d__2) - (r43 * c34 - c1345) 
		    * tmp3 * (d__3 * d__3) - c1345 * tmp2 * c->cvar_.u[(i__ + (
		    j + ((*k + 1) << 4)) * 10) * 5 - 741]);
	    c->cjac_.c__[((i__ + j * 6) * 5 + 2) * 5 - 176] = c->ctscon_.dt * 
		    c->cgcon_.tz2 * (c->cvar_.u[(i__ + (j + ((*k + 1) << 4)) * 10) *
		     5 - 744] * c->cvar_.u[(i__ + (j + ((*k + 1) << 4)) * 10) * 5 
		    - 742] * -.4 * tmp2) - c->ctscon_.dt * c->cgcon_.tz1 * (c34 - 
		    c1345) * tmp2 * c->cvar_.u[(i__ + (j + ((*k + 1) << 4)) * 10) 
		    * 5 - 744];
	    c->cjac_.c__[((i__ + j * 6) * 5 + 3) * 5 - 176] = c->ctscon_.dt * 
		    c->cgcon_.tz2 * (c->cvar_.u[(i__ + (j + ((*k + 1) << 4)) * 10) *
		     5 - 743] * c->cvar_.u[(i__ + (j + ((*k + 1) << 4)) * 10) * 5 
		    - 742] * -.4 * tmp2) - c->ctscon_.dt * c->cgcon_.tz1 * (c34 - 
		    c1345) * tmp2 * c->cvar_.u[(i__ + (j + ((*k + 1) << 4)) * 10) 
		    * 5 - 743];
	    c->cjac_.c__[((i__ + j * 6) * 5 + 4) * 5 - 176] = c->ctscon_.dt * 
		    c->cgcon_.tz2 * (c->cvar_.u[(i__ + (j + ((*k + 1) << 4)) * 10) *
		     5 - 741] * tmp1 * 1.4 - (c->cvar_.u[(i__ + (j + ((*k + 1) << 
		    4)) * 10) * 5 - 744] * c->cvar_.u[(i__ + (j + ((*k + 1) << 4))
		     * 10) * 5 - 744] + c->cvar_.u[(i__ + (j + ((*k + 1) << 4)) * 
		    10) * 5 - 743] * c->cvar_.u[(i__ + (j + ((*k + 1) << 4)) * 10)
		     * 5 - 743] + c->cvar_.u[(i__ + (j + ((*k + 1) << 4)) * 10) * 
		    5 - 742] * 3. * c->cvar_.u[(i__ + (j + ((*k + 1) << 4)) * 10) 
		    * 5 - 742]) * tmp2 * .20000000000000001) - c->ctscon_.dt * 
		    c->cgcon_.tz1 * (r43 * c34 - c1345) * tmp2 * c->cvar_.u[(i__ 
		    + (j + ((*k + 1) << 4)) * 10) * 5 - 742];
	    c->cjac_.c__[((i__ + j * 6) * 5 + 5) * 5 - 176] = c->ctscon_.dt * 
		    c->cgcon_.tz2 * (c->cvar_.u[(i__ + (j + ((*k + 1) << 4)) * 10) *
		     5 - 742] * tmp1 * 1.4) - c->ctscon_.dt * c->cgcon_.tz1 * 
		    c1345 * tmp1 - c->ctscon_.dt * c->cgcon_.tz1 * c->disp_.dz5;
	}
    }
    if (c->timer_.timeron) {
	timer_stop(6, c);
    }
    return 0;
} /* jacu_ */

/* Subroutine */ int jacld_(integer *k, context *c)
{
    integer i__1, i__2;
    float d__1, d__2, d__3, d__4, d__5, d__6, d__7, d__8, d__9;

    integer i__, j;
    float c34, r43;
    float c1345, tmp1, tmp2, tmp3;

    if (c->timer_.timeron) {
	timer_start(5, c);
    }
    r43 = 1.3333333333333333;
    c1345 = .19599999999999998;
    c34 = .10000000000000001;
    i__1 = c->cgcon_.jend;
    for (j = c->cgcon_.jst; j <= i__1; ++j) {
	i__2 = c->cgcon_.iend;
	for (i__ = c->cgcon_.ist; i__ <= i__2; ++i__) {
	    tmp1 = 1. / c->cvar_.u[(i__ + (j + (*k << 4)) * 10) * 5 - 745];
	    tmp2 = tmp1 * tmp1;
	    tmp3 = tmp1 * tmp2;
	    c->cjac_.d__[((i__ + j * 6) * 5 + 1) * 5 - 180] = c->ctscon_.dt * 2. *
		     (c->cgcon_.tx1 * c->disp_.dx1 + c->cgcon_.ty1 * c->disp_.dy1 + 
		    c->cgcon_.tz1 * c->disp_.dz1) + 1.;
	    c->cjac_.d__[((i__ + j * 6) * 5 + 2) * 5 - 180] = 0.;
	    c->cjac_.d__[((i__ + j * 6) * 5 + 3) * 5 - 180] = 0.;
	    c->cjac_.d__[((i__ + j * 6) * 5 + 4) * 5 - 180] = 0.;
	    c->cjac_.d__[((i__ + j * 6) * 5 + 5) * 5 - 180] = 0.;
	    c->cjac_.d__[((i__ + j * 6) * 5 + 1) * 5 - 179] = c->ctscon_.dt * 2. *
		     (c->cgcon_.tx1 * (-r43 * c34 * tmp2 * c->cvar_.u[(i__ + (j + 
		    (*k << 4)) * 10) * 5 - 744]) + c->cgcon_.ty1 * (-c34 * tmp2 
		    * c->cvar_.u[(i__ + (j + (*k << 4)) * 10) * 5 - 744]) + 
		    c->cgcon_.tz1 * (-c34 * tmp2 * c->cvar_.u[(i__ + (j + (*k << 
		    4)) * 10) * 5 - 744]));
	    c->cjac_.d__[((i__ + j * 6) * 5 + 2) * 5 - 179] = c->ctscon_.dt * 2. *
		     (c->cgcon_.tx1 * r43 * c34 * tmp1 + c->cgcon_.ty1 * c34 * 
		    tmp1 + c->cgcon_.tz1 * c34 * tmp1) + 1. + c->ctscon_.dt * 2. *
		     (c->cgcon_.tx1 * c->disp_.dx2 + c->cgcon_.ty1 * c->disp_.dy2 + 
		    c->cgcon_.tz1 * c->disp_.dz2);
	    c->cjac_.d__[((i__ + j * 6) * 5 + 3) * 5 - 179] = 0.;
	    c->cjac_.d__[((i__ + j * 6) * 5 + 4) * 5 - 179] = 0.;
	    c->cjac_.d__[((i__ + j * 6) * 5 + 5) * 5 - 179] = 0.;
	    c->cjac_.d__[((i__ + j * 6) * 5 + 1) * 5 - 178] = c->ctscon_.dt * 2. *
		     (c->cgcon_.tx1 * (-c34 * tmp2 * c->cvar_.u[(i__ + (j + (*k <<
		     4)) * 10) * 5 - 743]) + c->cgcon_.ty1 * (-r43 * c34 * tmp2 
		    * c->cvar_.u[(i__ + (j + (*k << 4)) * 10) * 5 - 743]) + 
		    c->cgcon_.tz1 * (-c34 * tmp2 * c->cvar_.u[(i__ + (j + (*k << 
		    4)) * 10) * 5 - 743]));
	    c->cjac_.d__[((i__ + j * 6) * 5 + 2) * 5 - 178] = 0.;
	    c->cjac_.d__[((i__ + j * 6) * 5 + 3) * 5 - 178] = c->ctscon_.dt * 2. *
		     (c->cgcon_.tx1 * c34 * tmp1 + c->cgcon_.ty1 * r43 * c34 * 
		    tmp1 + c->cgcon_.tz1 * c34 * tmp1) + 1. + c->ctscon_.dt * 2. *
		     (c->cgcon_.tx1 * c->disp_.dx3 + c->cgcon_.ty1 * c->disp_.dy3 + 
		    c->cgcon_.tz1 * c->disp_.dz3);
	    c->cjac_.d__[((i__ + j * 6) * 5 + 4) * 5 - 178] = 0.;
	    c->cjac_.d__[((i__ + j * 6) * 5 + 5) * 5 - 178] = 0.;
	    c->cjac_.d__[((i__ + j * 6) * 5 + 1) * 5 - 177] = c->ctscon_.dt * 2. *
		     (c->cgcon_.tx1 * (-c34 * tmp2 * c->cvar_.u[(i__ + (j + (*k <<
		     4)) * 10) * 5 - 742]) + c->cgcon_.ty1 * (-c34 * tmp2 * 
		    c->cvar_.u[(i__ + (j + (*k << 4)) * 10) * 5 - 742]) + 
		    c->cgcon_.tz1 * (-r43 * c34 * tmp2 * c->cvar_.u[(i__ + (j + (*
		    k << 4)) * 10) * 5 - 742]));
	    c->cjac_.d__[((i__ + j * 6) * 5 + 2) * 5 - 177] = 0.;
	    c->cjac_.d__[((i__ + j * 6) * 5 + 3) * 5 - 177] = 0.;
	    c->cjac_.d__[((i__ + j * 6) * 5 + 4) * 5 - 177] = c->ctscon_.dt * 2. *
		     (c->cgcon_.tx1 * c34 * tmp1 + c->cgcon_.ty1 * c34 * tmp1 + 
		    c->cgcon_.tz1 * r43 * c34 * tmp1) + 1. + c->ctscon_.dt * 2. * 
		    (c->cgcon_.tx1 * c->disp_.dx4 + c->cgcon_.ty1 * c->disp_.dy4 + 
		    c->cgcon_.tz1 * c->disp_.dz4);
	    c->cjac_.d__[((i__ + j * 6) * 5 + 5) * 5 - 177] = 0.;
	    d__1 = c->cvar_.u[(i__ + (j + (*k << 4)) * 10) * 5 - 744];
	    d__2 = c->cvar_.u[(i__ + (j + (*k << 4)) * 10) * 5 - 743];
	    d__3 = c->cvar_.u[(i__ + (j + (*k << 4)) * 10) * 5 - 742];
	    d__4 = c->cvar_.u[(i__ + (j + (*k << 4)) * 10) * 5 - 744];
	    d__5 = c->cvar_.u[(i__ + (j + (*k << 4)) * 10) * 5 - 743];
	    d__6 = c->cvar_.u[(i__ + (j + (*k << 4)) * 10) * 5 - 742];
	    d__7 = c->cvar_.u[(i__ + (j + (*k << 4)) * 10) * 5 - 744];
	    d__8 = c->cvar_.u[(i__ + (j + (*k << 4)) * 10) * 5 - 743];
	    d__9 = c->cvar_.u[(i__ + (j + (*k << 4)) * 10) * 5 - 742];
	    c->cjac_.d__[((i__ + j * 6) * 5 + 1) * 5 - 176] = c->ctscon_.dt * 2. *
		     (c->cgcon_.tx1 * (-(r43 * c34 - c1345) * tmp3 * (d__1 * 
		    d__1) - (c34 - c1345) * tmp3 * (d__2 * d__2) - (c34 - 
		    c1345) * tmp3 * (d__3 * d__3) - c1345 * tmp2 * c->cvar_.u[(
		    i__ + (j + (*k << 4)) * 10) * 5 - 741]) + c->cgcon_.ty1 * (
		    -(c34 - c1345) * tmp3 * (d__4 * d__4) - (r43 * c34 - 
		    c1345) * tmp3 * (d__5 * d__5) - (c34 - c1345) * tmp3 * (
		    d__6 * d__6) - c1345 * tmp2 * c->cvar_.u[(i__ + (j + (*k << 
		    4)) * 10) * 5 - 741]) + c->cgcon_.tz1 * (-(c34 - c1345) * 
		    tmp3 * (d__7 * d__7) - (c34 - c1345) * tmp3 * (d__8 * 
		    d__8) - (r43 * c34 - c1345) * tmp3 * (d__9 * d__9) - 
		    c1345 * tmp2 * c->cvar_.u[(i__ + (j + (*k << 4)) * 10) * 5 
		    - 741]));
	    c->cjac_.d__[((i__ + j * 6) * 5 + 2) * 5 - 176] = c->ctscon_.dt * 2. *
		     (c->cgcon_.tx1 * (r43 * c34 - c1345) * tmp2 * c->cvar_.u[(
		    i__ + (j + (*k << 4)) * 10) * 5 - 744] + c->cgcon_.ty1 * (
		    c34 - c1345) * tmp2 * c->cvar_.u[(i__ + (j + (*k << 4)) * 
		    10) * 5 - 744] + c->cgcon_.tz1 * (c34 - c1345) * tmp2 * 
		    c->cvar_.u[(i__ + (j + (*k << 4)) * 10) * 5 - 744]);
	    c->cjac_.d__[((i__ + j * 6) * 5 + 3) * 5 - 176] = c->ctscon_.dt * 2. *
		     (c->cgcon_.tx1 * (c34 - c1345) * tmp2 * c->cvar_.u[(i__ + (j 
		    + (*k << 4)) * 10) * 5 - 743] + c->cgcon_.ty1 * (r43 * c34 
		    - c1345) * tmp2 * c->cvar_.u[(i__ + (j + (*k << 4)) * 10) * 
		    5 - 743] + c->cgcon_.tz1 * (c34 - c1345) * tmp2 * c->cvar_.u[(
		    i__ + (j + (*k << 4)) * 10) * 5 - 743]);
	    c->cjac_.d__[((i__ + j * 6) * 5 + 4) * 5 - 176] = c->ctscon_.dt * 2. *
		     (c->cgcon_.tx1 * (c34 - c1345) * tmp2 * c->cvar_.u[(i__ + (j 
		    + (*k << 4)) * 10) * 5 - 742] + c->cgcon_.ty1 * (c34 - 
		    c1345) * tmp2 * c->cvar_.u[(i__ + (j + (*k << 4)) * 10) * 5 
		    - 742] + c->cgcon_.tz1 * (r43 * c34 - c1345) * tmp2 * 
		    c->cvar_.u[(i__ + (j + (*k << 4)) * 10) * 5 - 742]);
	    c->cjac_.d__[((i__ + j * 6) * 5 + 5) * 5 - 176] = c->ctscon_.dt * 2. *
		     (c->cgcon_.tx1 * c1345 * tmp1 + c->cgcon_.ty1 * c1345 * tmp1 
		    + c->cgcon_.tz1 * c1345 * tmp1) + 1. + c->ctscon_.dt * 2. * (
		    c->cgcon_.tx1 * c->disp_.dx5 + c->cgcon_.ty1 * c->disp_.dy5 + 
		    c->cgcon_.tz1 * c->disp_.dz5);
	    tmp1 = 1. / c->cvar_.u[(i__ + (j + ((*k - 1) << 4)) * 10) * 5 - 745];
	    tmp2 = tmp1 * tmp1;
	    tmp3 = tmp1 * tmp2;
	    c->cjac_.a[((i__ + j * 6) * 5 + 1) * 5 - 180] = -c->ctscon_.dt * 
		    c->cgcon_.tz1 * c->disp_.dz1;
	    c->cjac_.a[((i__ + j * 6) * 5 + 2) * 5 - 180] = 0.;
	    c->cjac_.a[((i__ + j * 6) * 5 + 3) * 5 - 180] = 0.;
	    c->cjac_.a[((i__ + j * 6) * 5 + 4) * 5 - 180] = -c->ctscon_.dt * 
		    c->cgcon_.tz2;
	    c->cjac_.a[((i__ + j * 6) * 5 + 5) * 5 - 180] = 0.;
	    c->cjac_.a[((i__ + j * 6) * 5 + 1) * 5 - 179] = -c->ctscon_.dt * 
		    c->cgcon_.tz2 * (-(c->cvar_.u[(i__ + (j + ((*k - 1) << 4)) * 10)
		     * 5 - 744] * c->cvar_.u[(i__ + (j + ((*k - 1) << 4)) * 10) * 
		    5 - 742]) * tmp2) - c->ctscon_.dt * c->cgcon_.tz1 * (-c34 * 
		    tmp2 * c->cvar_.u[(i__ + (j + ((*k - 1) << 4)) * 10) * 5 - 
		    744]);
	    c->cjac_.a[((i__ + j * 6) * 5 + 2) * 5 - 179] = -c->ctscon_.dt * 
		    c->cgcon_.tz2 * (c->cvar_.u[(i__ + (j + ((*k - 1) << 4)) * 10) *
		     5 - 742] * tmp1) - c->ctscon_.dt * c->cgcon_.tz1 * c34 * 
		    tmp1 - c->ctscon_.dt * c->cgcon_.tz1 * c->disp_.dz2;
	    c->cjac_.a[((i__ + j * 6) * 5 + 3) * 5 - 179] = 0.;
	    c->cjac_.a[((i__ + j * 6) * 5 + 4) * 5 - 179] = -c->ctscon_.dt * 
		    c->cgcon_.tz2 * (c->cvar_.u[(i__ + (j + ((*k - 1) << 4)) * 10) *
		     5 - 744] * tmp1);
	    c->cjac_.a[((i__ + j * 6) * 5 + 5) * 5 - 179] = 0.;
	    c->cjac_.a[((i__ + j * 6) * 5 + 1) * 5 - 178] = -c->ctscon_.dt * 
		    c->cgcon_.tz2 * (-(c->cvar_.u[(i__ + (j + ((*k - 1) << 4)) * 10)
		     * 5 - 743] * c->cvar_.u[(i__ + (j + ((*k - 1) << 4)) * 10) * 
		    5 - 742]) * tmp2) - c->ctscon_.dt * c->cgcon_.tz1 * (-c34 * 
		    tmp2 * c->cvar_.u[(i__ + (j + ((*k - 1) << 4)) * 10) * 5 - 
		    743]);
	    c->cjac_.a[((i__ + j * 6) * 5 + 2) * 5 - 178] = 0.;
	    c->cjac_.a[((i__ + j * 6) * 5 + 3) * 5 - 178] = -c->ctscon_.dt * 
		    c->cgcon_.tz2 * (c->cvar_.u[(i__ + (j + ((*k - 1) << 4)) * 10) *
		     5 - 742] * tmp1) - c->ctscon_.dt * c->cgcon_.tz1 * (c34 * 
		    tmp1) - c->ctscon_.dt * c->cgcon_.tz1 * c->disp_.dz3;
	    c->cjac_.a[((i__ + j * 6) * 5 + 4) * 5 - 178] = -c->ctscon_.dt * 
		    c->cgcon_.tz2 * (c->cvar_.u[(i__ + (j + ((*k - 1) << 4)) * 10) *
		     5 - 743] * tmp1);
	    c->cjac_.a[((i__ + j * 6) * 5 + 5) * 5 - 178] = 0.;
	    d__1 = c->cvar_.u[(i__ + (j + ((*k - 1) << 4)) * 10) * 5 - 742] * 
		    tmp1;
	    c->cjac_.a[((i__ + j * 6) * 5 + 1) * 5 - 177] = -c->ctscon_.dt * 
		    c->cgcon_.tz2 * (-(d__1 * d__1) + (c->cvar_.u[(i__ + (j + ((*k 
		    - 1) << 4)) * 10) * 5 - 744] * c->cvar_.u[(i__ + (j + ((*k - 
		    1) << 4)) * 10) * 5 - 744] + c->cvar_.u[(i__ + (j + ((*k - 1)
		    << 4)) * 10) * 5 - 743] * c->cvar_.u[(i__ + (j + ((*k - 1) << 
		    4)) * 10) * 5 - 743] + c->cvar_.u[(i__ + (j + ((*k - 1) << 4))
		     * 10) * 5 - 742] * c->cvar_.u[(i__ + (j + ((*k - 1) << 4)) * 
		    10) * 5 - 742]) * tmp2 * .20000000000000001) - 
		    c->ctscon_.dt * c->cgcon_.tz1 * (-r43 * c34 * tmp2 * c->cvar_.u[
		    (i__ + (j + ((*k - 1) << 4)) * 10) * 5 - 742]);
	    c->cjac_.a[((i__ + j * 6) * 5 + 2) * 5 - 177] = -c->ctscon_.dt * 
		    c->cgcon_.tz2 * (c->cvar_.u[(i__ + (j + ((*k - 1) << 4)) * 10) *
		     5 - 744] * tmp1 * -.4);
	    c->cjac_.a[((i__ + j * 6) * 5 + 3) * 5 - 177] = -c->ctscon_.dt * 
		    c->cgcon_.tz2 * (c->cvar_.u[(i__ + (j + ((*k - 1) << 4)) * 10) *
		     5 - 743] * tmp1 * -.4);
	    c->cjac_.a[((i__ + j * 6) * 5 + 4) * 5 - 177] = -c->ctscon_.dt * 
		    c->cgcon_.tz2 * 1.6000000000000001 * (c->cvar_.u[(i__ + (j + (
		    (*k - 1) << 4)) * 10) * 5 - 742] * tmp1) - c->ctscon_.dt * 
		    c->cgcon_.tz1 * (r43 * c34 * tmp1) - c->ctscon_.dt * 
		    c->cgcon_.tz1 * c->disp_.dz4;
	    c->cjac_.a[((i__ + j * 6) * 5 + 5) * 5 - 177] = -c->ctscon_.dt * 
		    c->cgcon_.tz2 * .4;
	    d__1 = c->cvar_.u[(i__ + (j + ((*k - 1) << 4)) * 10) * 5 - 744];
	    d__2 = c->cvar_.u[(i__ + (j + ((*k - 1) << 4)) * 10) * 5 - 743];
	    d__3 = c->cvar_.u[(i__ + (j + ((*k - 1) << 4)) * 10) * 5 - 742];
	    c->cjac_.a[((i__ + j * 6) * 5 + 1) * 5 - 176] = -c->ctscon_.dt * 
		    c->cgcon_.tz2 * (((c->cvar_.u[(i__ + (j + ((*k - 1) << 4)) * 10)
		     * 5 - 744] * c->cvar_.u[(i__ + (j + ((*k - 1) << 4)) * 10) * 
		    5 - 744] + c->cvar_.u[(i__ + (j + ((*k - 1) << 4)) * 10) * 5 
		    - 743] * c->cvar_.u[(i__ + (j + ((*k - 1) << 4)) * 10) * 5 - 
		    743] + c->cvar_.u[(i__ + (j + ((*k - 1) << 4)) * 10) * 5 - 
		    742] * c->cvar_.u[(i__ + (j + ((*k - 1) << 4)) * 10) * 5 - 
		    742]) * .4 * tmp2 - c->cvar_.u[(i__ + (j + ((*k - 1) << 4)) * 
		    10) * 5 - 741] * tmp1 * 1.4) * (c->cvar_.u[(i__ + (j + ((*k 
		    - 1) << 4)) * 10) * 5 - 742] * tmp1)) - c->ctscon_.dt * 
		    c->cgcon_.tz1 * (-(c34 - c1345) * tmp3 * (d__1 * d__1) - (
		    c34 - c1345) * tmp3 * (d__2 * d__2) - (r43 * c34 - c1345) 
		    * tmp3 * (d__3 * d__3) - c1345 * tmp2 * c->cvar_.u[(i__ + (
		    j + ((*k - 1) << 4)) * 10) * 5 - 741]);
	    c->cjac_.a[((i__ + j * 6) * 5 + 2) * 5 - 176] = -c->ctscon_.dt * 
		    c->cgcon_.tz2 * (c->cvar_.u[(i__ + (j + ((*k - 1) << 4)) * 10) *
		     5 - 744] * c->cvar_.u[(i__ + (j + ((*k - 1) << 4)) * 10) * 5 
		    - 742] * -.4 * tmp2) - c->ctscon_.dt * c->cgcon_.tz1 * (c34 - 
		    c1345) * tmp2 * c->cvar_.u[(i__ + (j + ((*k - 1) << 4)) * 10) 
		    * 5 - 744];
	    c->cjac_.a[((i__ + j * 6) * 5 + 3) * 5 - 176] = -c->ctscon_.dt * 
		    c->cgcon_.tz2 * (c->cvar_.u[(i__ + (j + ((*k - 1) << 4)) * 10) *
		     5 - 743] * c->cvar_.u[(i__ + (j + ((*k - 1) << 4)) * 10) * 5 
		    - 742] * -.4 * tmp2) - c->ctscon_.dt * c->cgcon_.tz1 * (c34 - 
		    c1345) * tmp2 * c->cvar_.u[(i__ + (j + ((*k - 1) << 4)) * 10) 
		    * 5 - 743];
	    c->cjac_.a[((i__ + j * 6) * 5 + 4) * 5 - 176] = -c->ctscon_.dt * 
		    c->cgcon_.tz2 * (c->cvar_.u[(i__ + (j + ((*k - 1) << 4)) * 10) *
		     5 - 741] * tmp1 * 1.4 - (c->cvar_.u[(i__ + (j + ((*k - 1) << 
		    4)) * 10) * 5 - 744] * c->cvar_.u[(i__ + (j + ((*k - 1) << 4))
		     * 10) * 5 - 744] + c->cvar_.u[(i__ + (j + ((*k - 1) << 4)) * 
		    10) * 5 - 743] * c->cvar_.u[(i__ + (j + ((*k - 1) << 4)) * 10)
		     * 5 - 743] + c->cvar_.u[(i__ + (j + ((*k - 1) << 4)) * 10) * 
		    5 - 742] * 3. * c->cvar_.u[(i__ + (j + ((*k - 1) << 4)) * 10) 
		    * 5 - 742]) * tmp2 * .20000000000000001) - c->ctscon_.dt * 
		    c->cgcon_.tz1 * (r43 * c34 - c1345) * tmp2 * c->cvar_.u[(i__ 
		    + (j + ((*k - 1) << 4)) * 10) * 5 - 742];
	    c->cjac_.a[((i__ + j * 6) * 5 + 5) * 5 - 176] = -c->ctscon_.dt * 
		    c->cgcon_.tz2 * (c->cvar_.u[(i__ + (j + ((*k - 1) << 4)) * 10) *
		     5 - 742] * tmp1 * 1.4) - c->ctscon_.dt * c->cgcon_.tz1 * 
		    c1345 * tmp1 - c->ctscon_.dt * c->cgcon_.tz1 * c->disp_.dz5;
	    tmp1 = 1. / c->cvar_.u[(i__ + (j - 1 + (*k << 4)) * 10) * 5 - 745];
	    tmp2 = tmp1 * tmp1;
	    tmp3 = tmp1 * tmp2;
	    c->cjac_.b[((i__ + j * 6) * 5 + 1) * 5 - 180] = -c->ctscon_.dt * 
		    c->cgcon_.ty1 * c->disp_.dy1;
	    c->cjac_.b[((i__ + j * 6) * 5 + 2) * 5 - 180] = 0.;
	    c->cjac_.b[((i__ + j * 6) * 5 + 3) * 5 - 180] = -c->ctscon_.dt * 
		    c->cgcon_.ty2;
	    c->cjac_.b[((i__ + j * 6) * 5 + 4) * 5 - 180] = 0.;
	    c->cjac_.b[((i__ + j * 6) * 5 + 5) * 5 - 180] = 0.;
	    c->cjac_.b[((i__ + j * 6) * 5 + 1) * 5 - 179] = -c->ctscon_.dt * 
		    c->cgcon_.ty2 * (-(c->cvar_.u[(i__ + (j - 1 + (*k << 4)) * 10)
		     * 5 - 744] * c->cvar_.u[(i__ + (j - 1 + (*k << 4)) * 10) * 
		    5 - 743]) * tmp2) - c->ctscon_.dt * c->cgcon_.ty1 * (-c34 * 
		    tmp2 * c->cvar_.u[(i__ + (j - 1 + (*k << 4)) * 10) * 5 - 
		    744]);
	    c->cjac_.b[((i__ + j * 6) * 5 + 2) * 5 - 179] = -c->ctscon_.dt * 
		    c->cgcon_.ty2 * (c->cvar_.u[(i__ + (j - 1 + (*k << 4)) * 10) *
		     5 - 743] * tmp1) - c->ctscon_.dt * c->cgcon_.ty1 * (c34 * 
		    tmp1) - c->ctscon_.dt * c->cgcon_.ty1 * c->disp_.dy2;
	    c->cjac_.b[((i__ + j * 6) * 5 + 3) * 5 - 179] = -c->ctscon_.dt * 
		    c->cgcon_.ty2 * (c->cvar_.u[(i__ + (j - 1 + (*k << 4)) * 10) *
		     5 - 744] * tmp1);
	    c->cjac_.b[((i__ + j * 6) * 5 + 4) * 5 - 179] = 0.;
	    c->cjac_.b[((i__ + j * 6) * 5 + 5) * 5 - 179] = 0.;
	    d__1 = c->cvar_.u[(i__ + (j - 1 + (*k << 4)) * 10) * 5 - 743] * 
		    tmp1;
	    c->cjac_.b[((i__ + j * 6) * 5 + 1) * 5 - 178] = -c->ctscon_.dt * 
		    c->cgcon_.ty2 * (-(d__1 * d__1) + (c->cvar_.u[(i__ + (j - 1 + 
		    (*k << 4)) * 10) * 5 - 744] * c->cvar_.u[(i__ + (j - 1 + (*
		    k << 4)) * 10) * 5 - 744] + c->cvar_.u[(i__ + (j - 1 + (*k 
		    << 4)) * 10) * 5 - 743] * c->cvar_.u[(i__ + (j - 1 + (*k << 
		    4)) * 10) * 5 - 743] + c->cvar_.u[(i__ + (j - 1 + (*k << 4))
		     * 10) * 5 - 742] * c->cvar_.u[(i__ + (j - 1 + (*k << 4)) * 
		    10) * 5 - 742]) * tmp2 * .20000000000000001) - 
		    c->ctscon_.dt * c->cgcon_.ty1 * (-r43 * c34 * tmp2 * c->cvar_.u[
		    (i__ + (j - 1 + (*k << 4)) * 10) * 5 - 743]);
	    c->cjac_.b[((i__ + j * 6) * 5 + 2) * 5 - 178] = -c->ctscon_.dt * 
		    c->cgcon_.ty2 * (c->cvar_.u[(i__ + (j - 1 + (*k << 4)) * 10) *
		     5 - 744] * tmp1 * -.4);
	    c->cjac_.b[((i__ + j * 6) * 5 + 3) * 5 - 178] = -c->ctscon_.dt * 
		    c->cgcon_.ty2 * (c->cvar_.u[(i__ + (j - 1 + (*k << 4)) * 10) *
		     5 - 743] * tmp1 * 1.6000000000000001) - c->ctscon_.dt * 
		    c->cgcon_.ty1 * (r43 * c34 * tmp1) - c->ctscon_.dt * 
		    c->cgcon_.ty1 * c->disp_.dy3;
	    c->cjac_.b[((i__ + j * 6) * 5 + 4) * 5 - 178] = -c->ctscon_.dt * 
		    c->cgcon_.ty2 * (c->cvar_.u[(i__ + (j - 1 + (*k << 4)) * 10) *
		     5 - 742] * tmp1 * -.4);
	    c->cjac_.b[((i__ + j * 6) * 5 + 5) * 5 - 178] = -c->ctscon_.dt * 
		    c->cgcon_.ty2 * .4;
	    c->cjac_.b[((i__ + j * 6) * 5 + 1) * 5 - 177] = -c->ctscon_.dt * 
		    c->cgcon_.ty2 * (-(c->cvar_.u[(i__ + (j - 1 + (*k << 4)) * 10)
		     * 5 - 743] * c->cvar_.u[(i__ + (j - 1 + (*k << 4)) * 10) * 
		    5 - 742]) * tmp2) - c->ctscon_.dt * c->cgcon_.ty1 * (-c34 * 
		    tmp2 * c->cvar_.u[(i__ + (j - 1 + (*k << 4)) * 10) * 5 - 
		    742]);
	    c->cjac_.b[((i__ + j * 6) * 5 + 2) * 5 - 177] = 0.;
	    c->cjac_.b[((i__ + j * 6) * 5 + 3) * 5 - 177] = -c->ctscon_.dt * 
		    c->cgcon_.ty2 * (c->cvar_.u[(i__ + (j - 1 + (*k << 4)) * 10) *
		     5 - 742] * tmp1);
	    c->cjac_.b[((i__ + j * 6) * 5 + 4) * 5 - 177] = -c->ctscon_.dt * 
		    c->cgcon_.ty2 * (c->cvar_.u[(i__ + (j - 1 + (*k << 4)) * 10) *
		     5 - 743] * tmp1) - c->ctscon_.dt * c->cgcon_.ty1 * (c34 * 
		    tmp1) - c->ctscon_.dt * c->cgcon_.ty1 * c->disp_.dy4;
	    c->cjac_.b[((i__ + j * 6) * 5 + 5) * 5 - 177] = 0.;
	    d__1 = c->cvar_.u[(i__ + (j - 1 + (*k << 4)) * 10) * 5 - 744];
	    d__2 = c->cvar_.u[(i__ + (j - 1 + (*k << 4)) * 10) * 5 - 743];
	    d__3 = c->cvar_.u[(i__ + (j - 1 + (*k << 4)) * 10) * 5 - 742];
	    c->cjac_.b[((i__ + j * 6) * 5 + 1) * 5 - 176] = -c->ctscon_.dt * 
		    c->cgcon_.ty2 * (((c->cvar_.u[(i__ + (j - 1 + (*k << 4)) * 10)
		     * 5 - 744] * c->cvar_.u[(i__ + (j - 1 + (*k << 4)) * 10) * 
		    5 - 744] + c->cvar_.u[(i__ + (j - 1 + (*k << 4)) * 10) * 5 
		    - 743] * c->cvar_.u[(i__ + (j - 1 + (*k << 4)) * 10) * 5 - 
		    743] + c->cvar_.u[(i__ + (j - 1 + (*k << 4)) * 10) * 5 - 
		    742] * c->cvar_.u[(i__ + (j - 1 + (*k << 4)) * 10) * 5 - 
		    742]) * .4 * tmp2 - c->cvar_.u[(i__ + (j - 1 + (*k << 4)) * 
		    10) * 5 - 741] * tmp1 * 1.4) * (c->cvar_.u[(i__ + (j - 1 + (
		    *k << 4)) * 10) * 5 - 743] * tmp1)) - c->ctscon_.dt * 
		    c->cgcon_.ty1 * (-(c34 - c1345) * tmp3 * (d__1 * d__1) - (
		    r43 * c34 - c1345) * tmp3 * (d__2 * d__2) - (c34 - c1345) 
		    * tmp3 * (d__3 * d__3) - c1345 * tmp2 * c->cvar_.u[(i__ + (
		    j - 1 + (*k << 4)) * 10) * 5 - 741]);
	    c->cjac_.b[((i__ + j * 6) * 5 + 2) * 5 - 176] = -c->ctscon_.dt * 
		    c->cgcon_.ty2 * (c->cvar_.u[(i__ + (j - 1 + (*k << 4)) * 10) *
		     5 - 744] * c->cvar_.u[(i__ + (j - 1 + (*k << 4)) * 10) * 5 
		    - 743] * -.4 * tmp2) - c->ctscon_.dt * c->cgcon_.ty1 * (c34 - 
		    c1345) * tmp2 * c->cvar_.u[(i__ + (j - 1 + (*k << 4)) * 10) 
		    * 5 - 744];
	    c->cjac_.b[((i__ + j * 6) * 5 + 3) * 5 - 176] = -c->ctscon_.dt * 
		    c->cgcon_.ty2 * (c->cvar_.u[(i__ + (j - 1 + (*k << 4)) * 10) *
		     5 - 741] * tmp1 * 1.4 - (c->cvar_.u[(i__ + (j - 1 + (*k << 
		    4)) * 10) * 5 - 744] * c->cvar_.u[(i__ + (j - 1 + (*k << 4))
		     * 10) * 5 - 744] + c->cvar_.u[(i__ + (j - 1 + (*k << 4)) * 
		    10) * 5 - 743] * 3. * c->cvar_.u[(i__ + (j - 1 + (*k << 4)) 
		    * 10) * 5 - 743] + c->cvar_.u[(i__ + (j - 1 + (*k << 4)) * 
		    10) * 5 - 742] * c->cvar_.u[(i__ + (j - 1 + (*k << 4)) * 10)
		     * 5 - 742]) * tmp2 * .20000000000000001) - c->ctscon_.dt * 
		    c->cgcon_.ty1 * (r43 * c34 - c1345) * tmp2 * c->cvar_.u[(i__ 
		    + (j - 1 + (*k << 4)) * 10) * 5 - 743];
	    c->cjac_.b[((i__ + j * 6) * 5 + 4) * 5 - 176] = -c->ctscon_.dt * 
		    c->cgcon_.ty2 * (c->cvar_.u[(i__ + (j - 1 + (*k << 4)) * 10) *
		     5 - 743] * c->cvar_.u[(i__ + (j - 1 + (*k << 4)) * 10) * 5 
		    - 742] * -.4 * tmp2) - c->ctscon_.dt * c->cgcon_.ty1 * (c34 - 
		    c1345) * tmp2 * c->cvar_.u[(i__ + (j - 1 + (*k << 4)) * 10) 
		    * 5 - 742];
	    c->cjac_.b[((i__ + j * 6) * 5 + 5) * 5 - 176] = -c->ctscon_.dt * 
		    c->cgcon_.ty2 * (c->cvar_.u[(i__ + (j - 1 + (*k << 4)) * 10) *
		     5 - 743] * tmp1 * 1.4) - c->ctscon_.dt * c->cgcon_.ty1 * 
		    c1345 * tmp1 - c->ctscon_.dt * c->cgcon_.ty1 * c->disp_.dy5;
	    tmp1 = 1. / c->cvar_.u[(i__ - 1 + (j + (*k << 4)) * 10) * 5 - 745];
	    tmp2 = tmp1 * tmp1;
	    tmp3 = tmp1 * tmp2;
	    c->cjac_.c__[((i__ + j * 6) * 5 + 1) * 5 - 180] = -c->ctscon_.dt * 
		    c->cgcon_.tx1 * c->disp_.dx1;
	    c->cjac_.c__[((i__ + j * 6) * 5 + 2) * 5 - 180] = -c->ctscon_.dt * 
		    c->cgcon_.tx2;
	    c->cjac_.c__[((i__ + j * 6) * 5 + 3) * 5 - 180] = 0.;
	    c->cjac_.c__[((i__ + j * 6) * 5 + 4) * 5 - 180] = 0.;
	    c->cjac_.c__[((i__ + j * 6) * 5 + 5) * 5 - 180] = 0.;
	    d__1 = c->cvar_.u[(i__ - 1 + (j + (*k << 4)) * 10) * 5 - 744] * 
		    tmp1;
	    c->cjac_.c__[((i__ + j * 6) * 5 + 1) * 5 - 179] = -c->ctscon_.dt * 
		    c->cgcon_.tx2 * (-(d__1 * d__1) + (c->cvar_.u[(i__ - 1 + (j + 
		    (*k << 4)) * 10) * 5 - 744] * c->cvar_.u[(i__ - 1 + (j + (*
		    k << 4)) * 10) * 5 - 744] + c->cvar_.u[(i__ - 1 + (j + (*k 
		    << 4)) * 10) * 5 - 743] * c->cvar_.u[(i__ - 1 + (j + (*k << 
		    4)) * 10) * 5 - 743] + c->cvar_.u[(i__ - 1 + (j + (*k << 4))
		     * 10) * 5 - 742] * c->cvar_.u[(i__ - 1 + (j + (*k << 4)) * 
		    10) * 5 - 742]) * .20000000000000001 * tmp2) - 
		    c->ctscon_.dt * c->cgcon_.tx1 * (-r43 * c34 * tmp2 * c->cvar_.u[
		    (i__ - 1 + (j + (*k << 4)) * 10) * 5 - 744]);
	    c->cjac_.c__[((i__ + j * 6) * 5 + 2) * 5 - 179] = -c->ctscon_.dt * 
		    c->cgcon_.tx2 * (c->cvar_.u[(i__ - 1 + (j + (*k << 4)) * 10) *
		     5 - 744] * tmp1 * 1.6000000000000001) - c->ctscon_.dt * 
		    c->cgcon_.tx1 * (r43 * c34 * tmp1) - c->ctscon_.dt * 
		    c->cgcon_.tx1 * c->disp_.dx2;
	    c->cjac_.c__[((i__ + j * 6) * 5 + 3) * 5 - 179] = -c->ctscon_.dt * 
		    c->cgcon_.tx2 * (c->cvar_.u[(i__ - 1 + (j + (*k << 4)) * 10) *
		     5 - 743] * tmp1 * -.4);
	    c->cjac_.c__[((i__ + j * 6) * 5 + 4) * 5 - 179] = -c->ctscon_.dt * 
		    c->cgcon_.tx2 * (c->cvar_.u[(i__ - 1 + (j + (*k << 4)) * 10) *
		     5 - 742] * tmp1 * -.4);
	    c->cjac_.c__[((i__ + j * 6) * 5 + 5) * 5 - 179] = -c->ctscon_.dt * 
		    c->cgcon_.tx2 * .4;
	    c->cjac_.c__[((i__ + j * 6) * 5 + 1) * 5 - 178] = -c->ctscon_.dt * 
		    c->cgcon_.tx2 * (-(c->cvar_.u[(i__ - 1 + (j + (*k << 4)) * 10)
		     * 5 - 744] * c->cvar_.u[(i__ - 1 + (j + (*k << 4)) * 10) * 
		    5 - 743]) * tmp2) - c->ctscon_.dt * c->cgcon_.tx1 * (-c34 * 
		    tmp2 * c->cvar_.u[(i__ - 1 + (j + (*k << 4)) * 10) * 5 - 
		    743]);
	    c->cjac_.c__[((i__ + j * 6) * 5 + 2) * 5 - 178] = -c->ctscon_.dt * 
		    c->cgcon_.tx2 * (c->cvar_.u[(i__ - 1 + (j + (*k << 4)) * 10) *
		     5 - 743] * tmp1);
	    c->cjac_.c__[((i__ + j * 6) * 5 + 3) * 5 - 178] = -c->ctscon_.dt * 
		    c->cgcon_.tx2 * (c->cvar_.u[(i__ - 1 + (j + (*k << 4)) * 10) *
		     5 - 744] * tmp1) - c->ctscon_.dt * c->cgcon_.tx1 * (c34 * 
		    tmp1) - c->ctscon_.dt * c->cgcon_.tx1 * c->disp_.dx3;
	    c->cjac_.c__[((i__ + j * 6) * 5 + 4) * 5 - 178] = 0.;
	    c->cjac_.c__[((i__ + j * 6) * 5 + 5) * 5 - 178] = 0.;
	    c->cjac_.c__[((i__ + j * 6) * 5 + 1) * 5 - 177] = -c->ctscon_.dt * 
		    c->cgcon_.tx2 * (-(c->cvar_.u[(i__ - 1 + (j + (*k << 4)) * 10)
		     * 5 - 744] * c->cvar_.u[(i__ - 1 + (j + (*k << 4)) * 10) * 
		    5 - 742]) * tmp2) - c->ctscon_.dt * c->cgcon_.tx1 * (-c34 * 
		    tmp2 * c->cvar_.u[(i__ - 1 + (j + (*k << 4)) * 10) * 5 - 
		    742]);
	    c->cjac_.c__[((i__ + j * 6) * 5 + 2) * 5 - 177] = -c->ctscon_.dt * 
		    c->cgcon_.tx2 * (c->cvar_.u[(i__ - 1 + (j + (*k << 4)) * 10) *
		     5 - 742] * tmp1);
	    c->cjac_.c__[((i__ + j * 6) * 5 + 3) * 5 - 177] = 0.;
	    c->cjac_.c__[((i__ + j * 6) * 5 + 4) * 5 - 177] = -c->ctscon_.dt * 
		    c->cgcon_.tx2 * (c->cvar_.u[(i__ - 1 + (j + (*k << 4)) * 10) *
		     5 - 744] * tmp1) - c->ctscon_.dt * c->cgcon_.tx1 * (c34 * 
		    tmp1) - c->ctscon_.dt * c->cgcon_.tx1 * c->disp_.dx4;
	    c->cjac_.c__[((i__ + j * 6) * 5 + 5) * 5 - 177] = 0.;
	    d__1 = c->cvar_.u[(i__ - 1 + (j + (*k << 4)) * 10) * 5 - 744];
	    d__2 = c->cvar_.u[(i__ - 1 + (j + (*k << 4)) * 10) * 5 - 743];
	    d__3 = c->cvar_.u[(i__ - 1 + (j + (*k << 4)) * 10) * 5 - 742];
	    c->cjac_.c__[((i__ + j * 6) * 5 + 1) * 5 - 176] = -c->ctscon_.dt * 
		    c->cgcon_.tx2 * (((c->cvar_.u[(i__ - 1 + (j + (*k << 4)) * 10)
		     * 5 - 744] * c->cvar_.u[(i__ - 1 + (j + (*k << 4)) * 10) * 
		    5 - 744] + c->cvar_.u[(i__ - 1 + (j + (*k << 4)) * 10) * 5 
		    - 743] * c->cvar_.u[(i__ - 1 + (j + (*k << 4)) * 10) * 5 - 
		    743] + c->cvar_.u[(i__ - 1 + (j + (*k << 4)) * 10) * 5 - 
		    742] * c->cvar_.u[(i__ - 1 + (j + (*k << 4)) * 10) * 5 - 
		    742]) * .4 * tmp2 - c->cvar_.u[(i__ - 1 + (j + (*k << 4)) * 
		    10) * 5 - 741] * tmp1 * 1.4) * (c->cvar_.u[(i__ - 1 + (j + (
		    *k << 4)) * 10) * 5 - 744] * tmp1)) - c->ctscon_.dt * 
		    c->cgcon_.tx1 * (-(r43 * c34 - c1345) * tmp3 * (d__1 * d__1)
		     - (c34 - c1345) * tmp3 * (d__2 * d__2) - (c34 - c1345) * 
		    tmp3 * (d__3 * d__3) - c1345 * tmp2 * c->cvar_.u[(i__ - 1 + 
		    (j + (*k << 4)) * 10) * 5 - 741]);
	    c->cjac_.c__[((i__ + j * 6) * 5 + 2) * 5 - 176] = -c->ctscon_.dt * 
		    c->cgcon_.tx2 * (c->cvar_.u[(i__ - 1 + (j + (*k << 4)) * 10) *
		     5 - 741] * tmp1 * 1.4 - (c->cvar_.u[(i__ - 1 + (j + (*k << 
		    4)) * 10) * 5 - 744] * 3. * c->cvar_.u[(i__ - 1 + (j + (*k 
		    << 4)) * 10) * 5 - 744] + c->cvar_.u[(i__ - 1 + (j + (*k << 
		    4)) * 10) * 5 - 743] * c->cvar_.u[(i__ - 1 + (j + (*k << 4))
		     * 10) * 5 - 743] + c->cvar_.u[(i__ - 1 + (j + (*k << 4)) * 
		    10) * 5 - 742] * c->cvar_.u[(i__ - 1 + (j + (*k << 4)) * 10)
		     * 5 - 742]) * tmp2 * .20000000000000001) - c->ctscon_.dt * 
		    c->cgcon_.tx1 * (r43 * c34 - c1345) * tmp2 * c->cvar_.u[(i__ 
		    - 1 + (j + (*k << 4)) * 10) * 5 - 744];
	    c->cjac_.c__[((i__ + j * 6) * 5 + 3) * 5 - 176] = -c->ctscon_.dt * 
		    c->cgcon_.tx2 * (c->cvar_.u[(i__ - 1 + (j + (*k << 4)) * 10) *
		     5 - 743] * c->cvar_.u[(i__ - 1 + (j + (*k << 4)) * 10) * 5 
		    - 744] * -.4 * tmp2) - c->ctscon_.dt * c->cgcon_.tx1 * (c34 - 
		    c1345) * tmp2 * c->cvar_.u[(i__ - 1 + (j + (*k << 4)) * 10) 
		    * 5 - 743];
	    c->cjac_.c__[((i__ + j * 6) * 5 + 4) * 5 - 176] = -c->ctscon_.dt * 
		    c->cgcon_.tx2 * (c->cvar_.u[(i__ - 1 + (j + (*k << 4)) * 10) *
		     5 - 742] * c->cvar_.u[(i__ - 1 + (j + (*k << 4)) * 10) * 5 
		    - 744] * -.4 * tmp2) - c->ctscon_.dt * c->cgcon_.tx1 * (c34 - 
		    c1345) * tmp2 * c->cvar_.u[(i__ - 1 + (j + (*k << 4)) * 10) 
		    * 5 - 742];
	    c->cjac_.c__[((i__ + j * 6) * 5 + 5) * 5 - 176] = -c->ctscon_.dt * 
		    c->cgcon_.tx2 * (c->cvar_.u[(i__ - 1 + (j + (*k << 4)) * 10) *
		     5 - 744] * tmp1 * 1.4) - c->ctscon_.dt * c->cgcon_.tx1 * 
		    c1345 * tmp1 - c->ctscon_.dt * c->cgcon_.tx1 * c->disp_.dx5;
	}
    }
    if (c->timer_.timeron) {
	timer_stop(5, c);
    }
    return 0;
} /* jacld_ */

/* Subroutine */ int init_comm__(context *c)
{


    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &c->dim_.id);
    MPI_Comm_size(MPI_COMM_WORLD, &c->dim_.num);
    if (c->dim_.num < NNODES_COMPILED) {
	if (c->dim_.id == 0) {
	}
	MPI_Abort(MPI_COMM_WORLD, 16);
    }
    c->dim_.ndim = nodedim_(&c->dim_.num);
	c->mpistuff_.dp_type__ = MPI_FLOAT;
    return 0;
} /* init_comm__ */

/* Subroutine */ int exchange_6__(float *g, integer *jbeg, integer *
	jfin1, context *c)
{
    integer i__1;

    integer k;
    //float dum[1024];
    MPI_Status status;
    MPI_Request msgid3;

    if (*jfin1 == c->cgcon_.ny) {
	MPI_Irecv(c->ftmp, c->cgcon_.nz, c->mpistuff_.dp_type__, c->neigh_.east, 3, MPI_COMM_WORLD, &msgid3);
	MPI_Wait(&msgid3, &status);
	i__1 = c->cgcon_.nz;
	for (k = 1; k <= i__1; ++k) {
	    g[c->cgcon_.ny + 1 + k * 14] = c->ftmp[k - 1];
	}
    }
    if (*jbeg == 1) {
	i__1 = c->cgcon_.nz;
	for (k = 1; k <= i__1; ++k) {
	    c->ftmp[k - 1] = g[k * 14 + 1];
	}
	MPI_Send(c->ftmp, c->cgcon_.nz, c->mpistuff_.dp_type__, c->neigh_.west, 3, MPI_COMM_WORLD);
    }
    return 0;
} /* exchange_6__ */

/* Subroutine */ int exchange_5__(float *g, integer *ibeg, integer *
	ifin1, context *c)
{
    integer i__1;

    integer k;
//    float dum[1024];
    MPI_Status status;
    MPI_Request msgid1;

    if (*ifin1 == c->cgcon_.nx) {
	MPI_Irecv(c->ftmp, c->cgcon_.nz, c->mpistuff_.dp_type__, c->neigh_.south, 1, MPI_COMM_WORLD, &msgid1);
	MPI_Wait(&msgid1, &status);
	i__1 = c->cgcon_.nz;
	for (k = 1; k <= i__1; ++k) {
	    g[c->cgcon_.nx + 1 + k * 14] = c->ftmp[k - 1];
	}
    }
    if (*ibeg == 1) {
	i__1 = c->cgcon_.nz;
	for (k = 1; k <= i__1; ++k) {
	    c->ftmp[k - 1] = g[k * 14 + 1];
	}
	MPI_Send(c->ftmp, c->cgcon_.nz, c->mpistuff_.dp_type__, c->neigh_.north, 1, MPI_COMM_WORLD);
    }
    return 0;
} /* exchange_5__ */

/* Subroutine */ int exchange_4__(float *g, float *h__, integer *
	ibeg, integer *ifin1, integer *jbeg, integer *jfin1, context *c)
{
    integer i__1;

    integer i__, j;
    integer ny2;
//    float dum[1024];
    MPI_Status status;
    MPI_Request msgid1;
    MPI_Request msgid3;

    ny2 = c->cgcon_.ny + 2;
    if (*jfin1 == c->cgcon_.ny) {
	i__1 = c->cgcon_.nx << 1;
	MPI_Irecv(c->ftmp, i__1, c->mpistuff_.dp_type__, c->neigh_.east, 3, MPI_COMM_WORLD, &msgid3);
	MPI_Wait(&msgid3, &status);
	i__1 = c->cgcon_.nx;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    g[i__ + (c->cgcon_.ny + 1) * 14] = c->ftmp[i__ - 1];
	    h__[i__ + (c->cgcon_.ny + 1) * 14] = c->ftmp[i__ + c->cgcon_.nx - 1];
	}
    }
    if (*jbeg == 1) {
	i__1 = c->cgcon_.nx;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    c->ftmp[i__ - 1] = g[i__ + 14];
	    c->ftmp[i__ + c->cgcon_.nx - 1] = h__[i__ + 14];
	}
	i__1 = c->cgcon_.nx << 1;
	MPI_Send(c->ftmp, i__1, c->mpistuff_.dp_type__, c->neigh_.west, 3, MPI_COMM_WORLD);
    }
    if (*ifin1 == c->cgcon_.nx) {
	i__1 = ny2 << 1;
	MPI_Irecv(c->ftmp, i__1, c->mpistuff_.dp_type__, c->neigh_.south, 1, 		MPI_COMM_WORLD, &msgid1);
	MPI_Wait(&msgid1, &status);
	i__1 = c->cgcon_.ny + 1;
	for (j = 0; j <= i__1; ++j) {
	    g[c->cgcon_.nx + 1 + j * 14] = c->ftmp[j];
	    h__[c->cgcon_.nx + 1 + j * 14] = c->ftmp[j + ny2];
	}
    }
    if (*ibeg == 1) {
	i__1 = c->cgcon_.ny + 1;
	for (j = 0; j <= i__1; ++j) {
	    c->ftmp[j] = g[j * 14 + 1];
	    c->ftmp[j + ny2] = h__[j * 14 + 1];
	}
	i__1 = ny2 << 1;
	MPI_Send(c->ftmp, i__1, c->mpistuff_.dp_type__, c->neigh_.north, 1, MPI_COMM_WORLD);
    }
    return 0;
} /* exchange_4__ */

/* Subroutine */ int exchange_3__(float *g, integer *iex, context *c)
{
    integer i__1, i__2;

    integer i__, j, k;
    MPI_Request mid;
    MPI_Status status;
    integer ipos1, ipos2;

    g -= 746;

    if (*iex == 0) {
	if (c->neigh_.north != -1) {
	    i__1 = c->cgcon_.ny * 10 * c->cgcon_.nz;
	    MPI_Irecv(c->comm_.buf1, i__1, c->mpistuff_.dp_type__, 		    c->neigh_.north, 2, MPI_COMM_WORLD, &mid);
	}
	if (c->neigh_.south != -1) {
	    i__1 = c->cgcon_.nz;
	    for (k = 1; k <= i__1; ++k) {
		i__2 = c->cgcon_.ny;
		for (j = 1; j <= i__2; ++j) {
		    ipos1 = (k - 1) * c->cgcon_.ny + j;
		    ipos2 = ipos1 + c->cgcon_.ny * c->cgcon_.nz;
		    c->comm_.buf[ipos1 * 5 - 5] = g[(c->cgcon_.nx - 1 + (j + (k <<
			     4)) * 10) * 5 + 1];
		    c->comm_.buf[ipos1 * 5 - 4] = g[(c->cgcon_.nx - 1 + (j + (k <<
			     4)) * 10) * 5 + 2];
		    c->comm_.buf[ipos1 * 5 - 3] = g[(c->cgcon_.nx - 1 + (j + (k <<
			     4)) * 10) * 5 + 3];
		    c->comm_.buf[ipos1 * 5 - 2] = g[(c->cgcon_.nx - 1 + (j + (k <<
			     4)) * 10) * 5 + 4];
		    c->comm_.buf[ipos1 * 5 - 1] = g[(c->cgcon_.nx - 1 + (j + (k <<
			     4)) * 10) * 5 + 5];
		    c->comm_.buf[ipos2 * 5 - 5] = g[(c->cgcon_.nx + (j + (k << 4))
			     * 10) * 5 + 1];
		    c->comm_.buf[ipos2 * 5 - 4] = g[(c->cgcon_.nx + (j + (k << 4))
			     * 10) * 5 + 2];
		    c->comm_.buf[ipos2 * 5 - 3] = g[(c->cgcon_.nx + (j + (k << 4))
			     * 10) * 5 + 3];
		    c->comm_.buf[ipos2 * 5 - 2] = g[(c->cgcon_.nx + (j + (k << 4))
			     * 10) * 5 + 4];
		    c->comm_.buf[ipos2 * 5 - 1] = g[(c->cgcon_.nx + (j + (k << 4))
			     * 10) * 5 + 5];
		}
	    }
	    i__1 = c->cgcon_.ny * 10 * c->cgcon_.nz;
	    MPI_Send(c->comm_.buf, i__1, c->mpistuff_.dp_type__, 		    c->neigh_.south, 2, MPI_COMM_WORLD);
	}
	if (c->neigh_.north != -1) {
	    MPI_Wait(&mid, &status);
	    i__1 = c->cgcon_.nz;
	    for (k = 1; k <= i__1; ++k) {
		i__2 = c->cgcon_.ny;
		for (j = 1; j <= i__2; ++j) {
		    ipos1 = (k - 1) * c->cgcon_.ny + j;
		    ipos2 = ipos1 + c->cgcon_.ny * c->cgcon_.nz;
		    g[((j + (k << 4)) * 10 - 1) * 5 + 1] = c->comm_.buf1[ipos1 *
			     5 - 5];
		    g[((j + (k << 4)) * 10 - 1) * 5 + 2] = c->comm_.buf1[ipos1 *
			     5 - 4];
		    g[((j + (k << 4)) * 10 - 1) * 5 + 3] = c->comm_.buf1[ipos1 *
			     5 - 3];
		    g[((j + (k << 4)) * 10 - 1) * 5 + 4] = c->comm_.buf1[ipos1 *
			     5 - 2];
		    g[((j + (k << 4)) * 10 - 1) * 5 + 5] = c->comm_.buf1[ipos1 *
			     5 - 1];
		    g[(j + (k << 4)) * 50 + 1] = c->comm_.buf1[ipos2 * 5 - 5];
		    g[(j + (k << 4)) * 50 + 2] = c->comm_.buf1[ipos2 * 5 - 4];
		    g[(j + (k << 4)) * 50 + 3] = c->comm_.buf1[ipos2 * 5 - 3];
		    g[(j + (k << 4)) * 50 + 4] = c->comm_.buf1[ipos2 * 5 - 2];
		    g[(j + (k << 4)) * 50 + 5] = c->comm_.buf1[ipos2 * 5 - 1];
		}
	    }
	}
	if (c->neigh_.south != -1) {
	    i__1 = c->cgcon_.ny * 10 * c->cgcon_.nz;
	    MPI_Irecv(c->comm_.buf1, i__1, c->mpistuff_.dp_type__, 		    c->neigh_.south, 1, MPI_COMM_WORLD, &mid);
	}
	if (c->neigh_.north != -1) {
	    i__1 = c->cgcon_.nz;
	    for (k = 1; k <= i__1; ++k) {
		i__2 = c->cgcon_.ny;
		for (j = 1; j <= i__2; ++j) {
		    ipos1 = (k - 1) * c->cgcon_.ny + j;
		    ipos2 = ipos1 + c->cgcon_.ny * c->cgcon_.nz;
		    c->comm_.buf[ipos1 * 5 - 5] = g[((j + (k << 4)) * 10 + 2) * 
			    5 + 1];
		    c->comm_.buf[ipos1 * 5 - 4] = g[((j + (k << 4)) * 10 + 2) * 
			    5 + 2];
		    c->comm_.buf[ipos1 * 5 - 3] = g[((j + (k << 4)) * 10 + 2) * 
			    5 + 3];
		    c->comm_.buf[ipos1 * 5 - 2] = g[((j + (k << 4)) * 10 + 2) * 
			    5 + 4];
		    c->comm_.buf[ipos1 * 5 - 1] = g[((j + (k << 4)) * 10 + 2) * 
			    5 + 5];
		    c->comm_.buf[ipos2 * 5 - 5] = g[((j + (k << 4)) * 10 + 1) * 
			    5 + 1];
		    c->comm_.buf[ipos2 * 5 - 4] = g[((j + (k << 4)) * 10 + 1) * 
			    5 + 2];
		    c->comm_.buf[ipos2 * 5 - 3] = g[((j + (k << 4)) * 10 + 1) * 
			    5 + 3];
		    c->comm_.buf[ipos2 * 5 - 2] = g[((j + (k << 4)) * 10 + 1) * 
			    5 + 4];
		    c->comm_.buf[ipos2 * 5 - 1] = g[((j + (k << 4)) * 10 + 1) * 
			    5 + 5];
		}
	    }
	    i__1 = c->cgcon_.ny * 10 * c->cgcon_.nz;
	    MPI_Send(c->comm_.buf, i__1, c->mpistuff_.dp_type__, 		    c->neigh_.north, 1, MPI_COMM_WORLD);
	}
	if (c->neigh_.south != -1) {
	    MPI_Wait(&mid, &status);
	    i__1 = c->cgcon_.nz;
	    for (k = 1; k <= i__1; ++k) {
		i__2 = c->cgcon_.ny;
		for (j = 1; j <= i__2; ++j) {
		    ipos1 = (k - 1) * c->cgcon_.ny + j;
		    ipos2 = ipos1 + c->cgcon_.ny * c->cgcon_.nz;
		    g[(c->cgcon_.nx + 2 + (j + (k << 4)) * 10) * 5 + 1] = 
			    c->comm_.buf1[ipos1 * 5 - 5];
		    g[(c->cgcon_.nx + 2 + (j + (k << 4)) * 10) * 5 + 2] = 
			    c->comm_.buf1[ipos1 * 5 - 4];
		    g[(c->cgcon_.nx + 2 + (j + (k << 4)) * 10) * 5 + 3] = 
			    c->comm_.buf1[ipos1 * 5 - 3];
		    g[(c->cgcon_.nx + 2 + (j + (k << 4)) * 10) * 5 + 4] = 
			    c->comm_.buf1[ipos1 * 5 - 2];
		    g[(c->cgcon_.nx + 2 + (j + (k << 4)) * 10) * 5 + 5] = 
			    c->comm_.buf1[ipos1 * 5 - 1];
		    g[(c->cgcon_.nx + 1 + (j + (k << 4)) * 10) * 5 + 1] = 
			    c->comm_.buf1[ipos2 * 5 - 5];
		    g[(c->cgcon_.nx + 1 + (j + (k << 4)) * 10) * 5 + 2] = 
			    c->comm_.buf1[ipos2 * 5 - 4];
		    g[(c->cgcon_.nx + 1 + (j + (k << 4)) * 10) * 5 + 3] = 
			    c->comm_.buf1[ipos2 * 5 - 3];
		    g[(c->cgcon_.nx + 1 + (j + (k << 4)) * 10) * 5 + 4] = 
			    c->comm_.buf1[ipos2 * 5 - 2];
		    g[(c->cgcon_.nx + 1 + (j + (k << 4)) * 10) * 5 + 5] = 
			    c->comm_.buf1[ipos2 * 5 - 1];
		}
	    }
	}
    } else {
	if (c->neigh_.west != -1) {
	    i__1 = c->cgcon_.nx * 10 * c->cgcon_.nz;
	    MPI_Irecv(c->comm_.buf1, i__1, c->mpistuff_.dp_type__, 		    c->neigh_.west, 4, MPI_COMM_WORLD, &mid);
	}
	if (c->neigh_.east != -1) {
	    i__1 = c->cgcon_.nz;
	    for (k = 1; k <= i__1; ++k) {
		i__2 = c->cgcon_.nx;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    ipos1 = (k - 1) * c->cgcon_.nx + i__;
		    ipos2 = ipos1 + c->cgcon_.nx * c->cgcon_.nz;
		    c->comm_.buf[ipos1 * 5 - 5] = g[(i__ + (c->cgcon_.ny - 1 + (k 
			    << 4)) * 10) * 5 + 1];
		    c->comm_.buf[ipos1 * 5 - 4] = g[(i__ + (c->cgcon_.ny - 1 + (k 
			    << 4)) * 10) * 5 + 2];
		    c->comm_.buf[ipos1 * 5 - 3] = g[(i__ + (c->cgcon_.ny - 1 + (k 
			    << 4)) * 10) * 5 + 3];
		    c->comm_.buf[ipos1 * 5 - 2] = g[(i__ + (c->cgcon_.ny - 1 + (k 
			    << 4)) * 10) * 5 + 4];
		    c->comm_.buf[ipos1 * 5 - 1] = g[(i__ + (c->cgcon_.ny - 1 + (k 
			    << 4)) * 10) * 5 + 5];
		    c->comm_.buf[ipos2 * 5 - 5] = g[(i__ + (c->cgcon_.ny + (k << 
			    4)) * 10) * 5 + 1];
		    c->comm_.buf[ipos2 * 5 - 4] = g[(i__ + (c->cgcon_.ny + (k << 
			    4)) * 10) * 5 + 2];
		    c->comm_.buf[ipos2 * 5 - 3] = g[(i__ + (c->cgcon_.ny + (k << 
			    4)) * 10) * 5 + 3];
		    c->comm_.buf[ipos2 * 5 - 2] = g[(i__ + (c->cgcon_.ny + (k << 
			    4)) * 10) * 5 + 4];
		    c->comm_.buf[ipos2 * 5 - 1] = g[(i__ + (c->cgcon_.ny + (k << 
			    4)) * 10) * 5 + 5];
		}
	    }
	    i__1 = c->cgcon_.nx * 10 * c->cgcon_.nz;
	    MPI_Send(c->comm_.buf, i__1, c->mpistuff_.dp_type__, 		    c->neigh_.east, 4, MPI_COMM_WORLD);
	}
	if (c->neigh_.west != -1) {
	    MPI_Wait(&mid, &status);
	    i__1 = c->cgcon_.nz;
	    for (k = 1; k <= i__1; ++k) {
		i__2 = c->cgcon_.nx;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    ipos1 = (k - 1) * c->cgcon_.nx + i__;
		    ipos2 = ipos1 + c->cgcon_.nx * c->cgcon_.nz;
		    g[(i__ + ((k << 4) - 1) * 10) * 5 + 1] = c->comm_.buf1[
			    ipos1 * 5 - 5];
		    g[(i__ + ((k << 4) - 1) * 10) * 5 + 2] = c->comm_.buf1[
			    ipos1 * 5 - 4];
		    g[(i__ + ((k << 4) - 1) * 10) * 5 + 3] = c->comm_.buf1[
			    ipos1 * 5 - 3];
		    g[(i__ + ((k << 4) - 1) * 10) * 5 + 4] = c->comm_.buf1[
			    ipos1 * 5 - 2];
		    g[(i__ + ((k << 4) - 1) * 10) * 5 + 5] = c->comm_.buf1[
			    ipos1 * 5 - 1];
		    g[(i__ + k * 160) * 5 + 1] = c->comm_.buf1[ipos2 * 5 - 5];
		    g[(i__ + k * 160) * 5 + 2] = c->comm_.buf1[ipos2 * 5 - 4];
		    g[(i__ + k * 160) * 5 + 3] = c->comm_.buf1[ipos2 * 5 - 3];
		    g[(i__ + k * 160) * 5 + 4] = c->comm_.buf1[ipos2 * 5 - 2];
		    g[(i__ + k * 160) * 5 + 5] = c->comm_.buf1[ipos2 * 5 - 1];
		}
	    }
	}
	if (c->neigh_.east != -1) {
	    i__1 = c->cgcon_.nx * 10 * c->cgcon_.nz;
	    MPI_Irecv(c->comm_.buf1, i__1, c->mpistuff_.dp_type__, 		    c->neigh_.east, 3, MPI_COMM_WORLD, &mid);
	}
	if (c->neigh_.west != -1) {
	    i__1 = c->cgcon_.nz;
	    for (k = 1; k <= i__1; ++k) {
		i__2 = c->cgcon_.nx;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    ipos1 = (k - 1) * c->cgcon_.nx + i__;
		    ipos2 = ipos1 + c->cgcon_.nx * c->cgcon_.nz;
		    c->comm_.buf[ipos1 * 5 - 5] = g[(i__ + ((k << 4) + 2) * 10) 
			    * 5 + 1];
		    c->comm_.buf[ipos1 * 5 - 4] = g[(i__ + ((k << 4) + 2) * 10) 
			    * 5 + 2];
		    c->comm_.buf[ipos1 * 5 - 3] = g[(i__ + ((k << 4) + 2) * 10) 
			    * 5 + 3];
		    c->comm_.buf[ipos1 * 5 - 2] = g[(i__ + ((k << 4) + 2) * 10) 
			    * 5 + 4];
		    c->comm_.buf[ipos1 * 5 - 1] = g[(i__ + ((k << 4) + 2) * 10) 
			    * 5 + 5];
		    c->comm_.buf[ipos2 * 5 - 5] = g[(i__ + ((k << 4) + 1) * 10) 
			    * 5 + 1];
		    c->comm_.buf[ipos2 * 5 - 4] = g[(i__ + ((k << 4) + 1) * 10) 
			    * 5 + 2];
		    c->comm_.buf[ipos2 * 5 - 3] = g[(i__ + ((k << 4) + 1) * 10) 
			    * 5 + 3];
		    c->comm_.buf[ipos2 * 5 - 2] = g[(i__ + ((k << 4) + 1) * 10) 
			    * 5 + 4];
		    c->comm_.buf[ipos2 * 5 - 1] = g[(i__ + ((k << 4) + 1) * 10) 
			    * 5 + 5];
		}
	    }
	    i__1 = c->cgcon_.nx * 10 * c->cgcon_.nz;
	    MPI_Send(c->comm_.buf, i__1, c->mpistuff_.dp_type__, 		    c->neigh_.west, 3, MPI_COMM_WORLD);
	}
	if (c->neigh_.east != -1) {
	    MPI_Wait(&mid, &status);
	    i__1 = c->cgcon_.nz;
	    for (k = 1; k <= i__1; ++k) {
		i__2 = c->cgcon_.nx;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    ipos1 = (k - 1) * c->cgcon_.nx + i__;
		    ipos2 = ipos1 + c->cgcon_.nx * c->cgcon_.nz;
		    g[(i__ + (c->cgcon_.ny + 2 + (k << 4)) * 10) * 5 + 1] = 
			    c->comm_.buf1[ipos1 * 5 - 5];
		    g[(i__ + (c->cgcon_.ny + 2 + (k << 4)) * 10) * 5 + 2] = 
			    c->comm_.buf1[ipos1 * 5 - 4];
		    g[(i__ + (c->cgcon_.ny + 2 + (k << 4)) * 10) * 5 + 3] = 
			    c->comm_.buf1[ipos1 * 5 - 3];
		    g[(i__ + (c->cgcon_.ny + 2 + (k << 4)) * 10) * 5 + 4] = 
			    c->comm_.buf1[ipos1 * 5 - 2];
		    g[(i__ + (c->cgcon_.ny + 2 + (k << 4)) * 10) * 5 + 5] = 
			    c->comm_.buf1[ipos1 * 5 - 1];
		    g[(i__ + (c->cgcon_.ny + 1 + (k << 4)) * 10) * 5 + 1] = 
			    c->comm_.buf1[ipos2 * 5 - 5];
		    g[(i__ + (c->cgcon_.ny + 1 + (k << 4)) * 10) * 5 + 2] = 
			    c->comm_.buf1[ipos2 * 5 - 4];
		    g[(i__ + (c->cgcon_.ny + 1 + (k << 4)) * 10) * 5 + 3] = 
			    c->comm_.buf1[ipos2 * 5 - 3];
		    g[(i__ + (c->cgcon_.ny + 1 + (k << 4)) * 10) * 5 + 4] = 
			    c->comm_.buf1[ipos2 * 5 - 2];
		    g[(i__ + (c->cgcon_.ny + 1 + (k << 4)) * 10) * 5 + 5] = 
			    c->comm_.buf1[ipos2 * 5 - 1];
		}
	    }
	}
    }
    return 0;
} /* exchange_3__ */

/* Subroutine */ int exact_(integer *i__, integer *j, integer *k, float *
	u000ijk, context *c)
{
    integer m;
    float xi, eta, zeta;

    --u000ijk;

    xi = (float) (*i__ - 1) / (c->cgcon_.nx0 - 1);
    eta = (float) (*j - 1) / (c->cgcon_.ny0 - 1);
    zeta = (float) (*k - 1) / (c->cgcon_.nz - 1);
    for (m = 1; m <= 5; ++m) {
	u000ijk[m] = c->cexact_.ce[m - 1] + c->cexact_.ce[m + 4] * xi + 
		c->cexact_.ce[m + 9] * eta + c->cexact_.ce[m + 14] * zeta + 
		c->cexact_.ce[m + 19] * xi * xi + c->cexact_.ce[m + 24] * eta * 
		eta + c->cexact_.ce[m + 29] * zeta * zeta + c->cexact_.ce[m + 34] 
		* xi * xi * xi + c->cexact_.ce[m + 39] * eta * eta * eta + 
		c->cexact_.ce[m + 44] * zeta * zeta * zeta + c->cexact_.ce[m + 49]
		 * xi * xi * xi * xi + c->cexact_.ce[m + 54] * eta * eta * eta *
		 eta + c->cexact_.ce[m + 59] * zeta * zeta * zeta * zeta;
    }
    return 0;
} /* exact_ */

/* Subroutine */ int error_(context *c)
{
    integer i__1, i__2, i__3;
    float d__1;

    integer i__, j, k, m;
    float tmp;
    integer iglob, jglob;
    float dummy[5];
    float u000ijk[5];

    for (m = 1; m <= 5; ++m) {
	c->ctscon_.errnm[m - 1] = 0.;
	dummy[m - 1] = 0.;
    }
    i__1 = c->cgcon_.nz - 1;
    for (k = 2; k <= i__1; ++k) {
	i__2 = c->cgcon_.jend;
	for (j = c->cgcon_.jst; j <= i__2; ++j) {
	    jglob = c->cgcon_.jpt + j;
	    i__3 = c->cgcon_.iend;
	    for (i__ = c->cgcon_.ist; i__ <= i__3; ++i__) {
		iglob = c->cgcon_.ipt + i__;
		exact_(&iglob, &jglob, &k, u000ijk, c);
		for (m = 1; m <= 5; ++m) {
		    tmp = u000ijk[m - 1] - c->cvar_.u[m + (i__ + (j + (k << 4)) 
			    * 10) * 5 - 746];
		    d__1 = tmp;
		    dummy[m - 1] += d__1 * d__1;
		}
	    }
	}
    }
    
    MPI_Allreduce(dummy, c->ctscon_.errnm, 5, c->mpistuff_.dp_type__, MPI_SUM, MPI_COMM_WORLD);
    for (m = 1; m <= 5; ++m) {
	c->ctscon_.errnm[m - 1] = sqrtf(c->ctscon_.errnm[m - 1] / ((c->cgcon_.nx0 - 
		2) * (c->cgcon_.ny0 - 2) * (c->cgcon_.nz0 - 2)));
    }
    return 0;
} /* error_ */

/* Subroutine */ int erhs_(context *c)
{
    integer i__1, i__2, i__3;
    float d__1, d__2, d__3, d__4, d__5, d__6, d__7, d__8;

    integer i__, j, k, m;
    float q;
    integer l1, l2;
    float u21, u31, u41, xi, eta, u21i, u31i, u41i, u51i;
    integer iex;
    float u21j, u31j, u41j, u51j, u21k, tmp, u31k, u41k, u51k;
    integer ist1, jst1;
    float zeta;
    integer iend1, jend1;
    float u21im1, u31im1, u41im1, u51im1, u21jm1, u31jm1, u41jm1, 
	    u51jm1, u21km1, u31km1, u41km1, u51km1;
    integer iglob, jglob;
    float dsspm;

    dsspm = c->disp_.dssp;
    i__1 = c->cgcon_.nz;
    for (k = 1; k <= i__1; ++k) {
	i__2 = c->cgcon_.ny;
	for (j = 1; j <= i__2; ++j) {
	    i__3 = c->cgcon_.nx;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		for (m = 1; m <= 5; ++m) {
		    c->cvar_.frct[m + (i__ + (j + (k << 4)) * 10) * 5 - 746] = 
			    0.;
		}
	    }
	}
    }
    i__1 = c->cgcon_.nz;
    for (k = 1; k <= i__1; ++k) {
	zeta = (float) (k - 1) / (c->cgcon_.nz - 1);
	i__2 = c->cgcon_.ny;
	for (j = 1; j <= i__2; ++j) {
	    jglob = c->cgcon_.jpt + j;
	    eta = (float) (jglob - 1) / (c->cgcon_.ny0 - 1);
	    i__3 = c->cgcon_.nx;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		iglob = c->cgcon_.ipt + i__;
		xi = (float) (iglob - 1) / (c->cgcon_.nx0 - 1);
		for (m = 1; m <= 5; ++m) {
		    c->cvar_.rsd[m + (i__ + (j + (k << 4)) * 10) * 5 - 746] = 
			    c->cexact_.ce[m - 1] + c->cexact_.ce[m + 4] * xi + 
			    c->cexact_.ce[m + 9] * eta + c->cexact_.ce[m + 14] * 
			    zeta + c->cexact_.ce[m + 19] * xi * xi + 
			    c->cexact_.ce[m + 24] * eta * eta + c->cexact_.ce[m + 
			    29] * zeta * zeta + c->cexact_.ce[m + 34] * xi * xi 
			    * xi + c->cexact_.ce[m + 39] * eta * eta * eta + 
			    c->cexact_.ce[m + 44] * zeta * zeta * zeta + 
			    c->cexact_.ce[m + 49] * xi * xi * xi * xi + 
			    c->cexact_.ce[m + 54] * eta * eta * eta * eta + 
			    c->cexact_.ce[m + 59] * zeta * zeta * zeta * zeta;
		}
	    }
	}
    }

    iex = 0;
    exchange_3__(c->cvar_.rsd, &iex, c);
    l1 = 0;
    if (c->neigh_.north == -1) {
	l1 = 1;
    }
    l2 = c->cgcon_.nx + 1;
    if (c->neigh_.south == -1) {
	l2 = c->cgcon_.nx;
    }
    ist1 = 1;
    iend1 = c->cgcon_.nx;
    if (c->neigh_.north == -1) {
	ist1 = 4;
    }
    if (c->neigh_.south == -1) {
	iend1 = c->cgcon_.nx - 3;
    }
    i__1 = c->cgcon_.nz - 1;
    for (k = 2; k <= i__1; ++k) {
	i__2 = c->cgcon_.jend;
	for (j = c->cgcon_.jst; j <= i__2; ++j) {
	    i__3 = l2;
	    for (i__ = l1; i__ <= i__3; ++i__) {
		c->cvar_.flux[(i__ + ((j + k * 14) << 3)) * 5 - 560] = c->cvar_.rsd[
			(i__ + (j + (k << 4)) * 10) * 5 - 744];
		u21 = c->cvar_.rsd[(i__ + (j + (k << 4)) * 10) * 5 - 744] / 
			c->cvar_.rsd[(i__ + (j + (k << 4)) * 10) * 5 - 745];
		q = (c->cvar_.rsd[(i__ + (j + (k << 4)) * 10) * 5 - 744] * 
			c->cvar_.rsd[(i__ + (j + (k << 4)) * 10) * 5 - 744] + 
			c->cvar_.rsd[(i__ + (j + (k << 4)) * 10) * 5 - 743] * 
			c->cvar_.rsd[(i__ + (j + (k << 4)) * 10) * 5 - 743] + 
			c->cvar_.rsd[(i__ + (j + (k << 4)) * 10) * 5 - 742] * 
			c->cvar_.rsd[(i__ + (j + (k << 4)) * 10) * 5 - 742]) * 
			.5 / c->cvar_.rsd[(i__ + (j + (k << 4)) * 10) * 5 - 745]
			;
		c->cvar_.flux[(i__ + ((j + k * 14) << 3)) * 5 - 559] = c->cvar_.rsd[
			(i__ + (j + (k << 4)) * 10) * 5 - 744] * u21 + (
			c->cvar_.rsd[(i__ + (j + (k << 4)) * 10) * 5 - 741] - q)
			 * .4;
		c->cvar_.flux[(i__ + ((j + k * 14) << 3)) * 5 - 558] = c->cvar_.rsd[
			(i__ + (j + (k << 4)) * 10) * 5 - 743] * u21;
		c->cvar_.flux[(i__ + ((j + k * 14) << 3)) * 5 - 557] = c->cvar_.rsd[
			(i__ + (j + (k << 4)) * 10) * 5 - 742] * u21;
		c->cvar_.flux[(i__ + ((j + k * 14) << 3)) * 5 - 556] = (
			c->cvar_.rsd[(i__ + (j + (k << 4)) * 10) * 5 - 741] * 
			1.4 - q * .4) * u21;
	    }
	}
    }
    i__1 = c->cgcon_.nz - 1;
    for (k = 2; k <= i__1; ++k) {
	i__2 = c->cgcon_.jend;
	for (j = c->cgcon_.jst; j <= i__2; ++j) {
	    i__3 = c->cgcon_.iend;
	    for (i__ = c->cgcon_.ist; i__ <= i__3; ++i__) {
		for (m = 1; m <= 5; ++m) {
		    c->cvar_.frct[m + (i__ + (j + (k << 4)) * 10) * 5 - 746] -= 
			    c->cgcon_.tx2 * (c->cvar_.flux[m + (i__ + 1 + ((j + k *
			     14) << 3)) * 5 - 561] - c->cvar_.flux[m + (i__ - 1 
			    + ((j + k * 14) << 3)) * 5 - 561]);
		}
	    }
	    i__3 = l2;
	    for (i__ = c->cgcon_.ist; i__ <= i__3; ++i__) {
		tmp = 1. / c->cvar_.rsd[(i__ + (j + (k << 4)) * 10) * 5 - 745];
		u21i = tmp * c->cvar_.rsd[(i__ + (j + (k << 4)) * 10) * 5 - 744]
			;
		u31i = tmp * c->cvar_.rsd[(i__ + (j + (k << 4)) * 10) * 5 - 743]
			;
		u41i = tmp * c->cvar_.rsd[(i__ + (j + (k << 4)) * 10) * 5 - 742]
			;
		u51i = tmp * c->cvar_.rsd[(i__ + (j + (k << 4)) * 10) * 5 - 741]
			;
		tmp = 1. / c->cvar_.rsd[(i__ - 1 + (j + (k << 4)) * 10) * 5 - 
			745];
		u21im1 = tmp * c->cvar_.rsd[(i__ - 1 + (j + (k << 4)) * 10) * 5 
			- 744];
		u31im1 = tmp * c->cvar_.rsd[(i__ - 1 + (j + (k << 4)) * 10) * 5 
			- 743];
		u41im1 = tmp * c->cvar_.rsd[(i__ - 1 + (j + (k << 4)) * 10) * 5 
			- 742];
		u51im1 = tmp * c->cvar_.rsd[(i__ - 1 + (j + (k << 4)) * 10) * 5 
			- 741];
		c->cvar_.flux[(i__ + ((j + k * 14) << 3)) * 5 - 559] = 
			c->cgcon_.tx3 * 1.3333333333333333 * (u21i - u21im1);
		c->cvar_.flux[(i__ + ((j + k * 14) << 3)) * 5 - 558] = 
			c->cgcon_.tx3 * (u31i - u31im1);
		c->cvar_.flux[(i__ + ((j + k * 14) << 3)) * 5 - 557] = 
			c->cgcon_.tx3 * (u41i - u41im1);
		d__1 = u21i;
		d__2 = u31i;
		d__3 = u41i;
		d__4 = u21im1;
		d__5 = u31im1;
		d__6 = u41im1;
		d__7 = u21i;
		d__8 = u21im1;
		c->cvar_.flux[(i__ + ((j + k * 14) << 3)) * 5 - 556] = 
			c->cgcon_.tx3 * -.47999999999999987 * (d__1 * d__1 + 
			d__2 * d__2 + d__3 * d__3 - (d__4 * d__4 + d__5 * 
			d__5 + d__6 * d__6)) + c->cgcon_.tx3 * 
			.16666666666666666 * (d__7 * d__7 - d__8 * d__8) + 
			c->cgcon_.tx3 * 1.9599999999999997 * (u51i - u51im1);
	    }
	    i__3 = c->cgcon_.iend;
	    for (i__ = c->cgcon_.ist; i__ <= i__3; ++i__) {
		c->cvar_.frct[(i__ + (j + (k << 4)) * 10) * 5 - 745] += 
			c->disp_.dx1 * c->cgcon_.tx1 * (c->cvar_.rsd[(i__ - 1 + (j 
			+ (k << 4)) * 10) * 5 - 745] - c->cvar_.rsd[(i__ + (j + 
			(k << 4)) * 10) * 5 - 745] * 2. + c->cvar_.rsd[(i__ + 1 
			+ (j + (k << 4)) * 10) * 5 - 745]);
		c->cvar_.frct[(i__ + (j + (k << 4)) * 10) * 5 - 744] = 
			c->cvar_.frct[(i__ + (j + (k << 4)) * 10) * 5 - 744] + 
			c->cgcon_.tx3 * .1 * 1. * (c->cvar_.flux[(i__ + 1 + ((j + 
			k * 14) << 3)) * 5 - 559] - c->cvar_.flux[(i__ + ((j + k *
			 14) << 3)) * 5 - 559]) + c->disp_.dx2 * c->cgcon_.tx1 * (
			c->cvar_.rsd[(i__ - 1 + (j + (k << 4)) * 10) * 5 - 744] 
			- c->cvar_.rsd[(i__ + (j + (k << 4)) * 10) * 5 - 744] * 
			2. + c->cvar_.rsd[(i__ + 1 + (j + (k << 4)) * 10) * 5 - 
			744]);
		c->cvar_.frct[(i__ + (j + (k << 4)) * 10) * 5 - 743] = 
			c->cvar_.frct[(i__ + (j + (k << 4)) * 10) * 5 - 743] + 
			c->cgcon_.tx3 * .1 * 1. * (c->cvar_.flux[(i__ + 1 + ((j + 
			k * 14) << 3)) * 5 - 558] - c->cvar_.flux[(i__ + ((j + k *
			 14) << 3)) * 5 - 558]) + c->disp_.dx3 * c->cgcon_.tx1 * (
			c->cvar_.rsd[(i__ - 1 + (j + (k << 4)) * 10) * 5 - 743] 
			- c->cvar_.rsd[(i__ + (j + (k << 4)) * 10) * 5 - 743] * 
			2. + c->cvar_.rsd[(i__ + 1 + (j + (k << 4)) * 10) * 5 - 
			743]);
		c->cvar_.frct[(i__ + (j + (k << 4)) * 10) * 5 - 742] = 
			c->cvar_.frct[(i__ + (j + (k << 4)) * 10) * 5 - 742] + 
			c->cgcon_.tx3 * .1 * 1. * (c->cvar_.flux[(i__ + 1 + ((j + 
			k * 14) << 3)) * 5 - 557] - c->cvar_.flux[(i__ + ((j + k *
			 14) << 3)) * 5 - 557]) + c->disp_.dx4 * c->cgcon_.tx1 * (
			c->cvar_.rsd[(i__ - 1 + (j + (k << 4)) * 10) * 5 - 742] 
			- c->cvar_.rsd[(i__ + (j + (k << 4)) * 10) * 5 - 742] * 
			2. + c->cvar_.rsd[(i__ + 1 + (j + (k << 4)) * 10) * 5 - 
			742]);
		c->cvar_.frct[(i__ + (j + (k << 4)) * 10) * 5 - 741] = 
			c->cvar_.frct[(i__ + (j + (k << 4)) * 10) * 5 - 741] + 
			c->cgcon_.tx3 * .1 * 1. * (c->cvar_.flux[(i__ + 1 + ((j + 
			k * 14) << 3)) * 5 - 556] - c->cvar_.flux[(i__ + ((j + k *
			 14) << 3)) * 5 - 556]) + c->disp_.dx5 * c->cgcon_.tx1 * (
			c->cvar_.rsd[(i__ - 1 + (j + (k << 4)) * 10) * 5 - 741] 
			- c->cvar_.rsd[(i__ + (j + (k << 4)) * 10) * 5 - 741] * 
			2. + c->cvar_.rsd[(i__ + 1 + (j + (k << 4)) * 10) * 5 - 
			741]);
	    }
	    if (c->neigh_.north == -1) {
		for (m = 1; m <= 5; ++m) {
		    c->cvar_.frct[m + ((j + (k << 4)) * 10 + 2) * 5 - 746] -= 
			    dsspm * (c->cvar_.rsd[m + ((j + (k << 4)) * 10 + 2) 
			    * 5 - 746] * 5. - c->cvar_.rsd[m + ((j + (k << 4)) *
			     10 + 3) * 5 - 746] * 4. + c->cvar_.rsd[m + ((j + (
			    k << 4)) * 10 + 4) * 5 - 746]);
		    c->cvar_.frct[m + ((j + (k << 4)) * 10 + 3) * 5 - 746] -= 
			    dsspm * (c->cvar_.rsd[m + ((j + (k << 4)) * 10 + 2) 
			    * 5 - 746] * -4. + c->cvar_.rsd[m + ((j + (k << 4)) 
			    * 10 + 3) * 5 - 746] * 6. - c->cvar_.rsd[m + ((j + (
			    k << 4)) * 10 + 4) * 5 - 746] * 4. + c->cvar_.rsd[m 
			    + ((j + (k << 4)) * 10 + 5) * 5 - 746]);
		}
	    }
	    i__3 = iend1;
	    for (i__ = ist1; i__ <= i__3; ++i__) {
		for (m = 1; m <= 5; ++m) {
		    c->cvar_.frct[m + (i__ + (j + (k << 4)) * 10) * 5 - 746] -= 
			    dsspm * (c->cvar_.rsd[m + (i__ - 2 + (j + (k << 4)) 
			    * 10) * 5 - 746] - c->cvar_.rsd[m + (i__ - 1 + (j + 
			    (k << 4)) * 10) * 5 - 746] * 4. + c->cvar_.rsd[m + (
			    i__ + (j + (k << 4)) * 10) * 5 - 746] * 6. - 
			    c->cvar_.rsd[m + (i__ + 1 + (j + (k << 4)) * 10) * 
			    5 - 746] * 4. + c->cvar_.rsd[m + (i__ + 2 + (j + (k 
			    << 4)) * 10) * 5 - 746]);
		}
	    }
	    if (c->neigh_.south == -1) {
		for (m = 1; m <= 5; ++m) {
		    c->cvar_.frct[m + (c->cgcon_.nx - 2 + (j + (k << 4)) * 10) * 
			    5 - 746] -= dsspm * (c->cvar_.rsd[m + (c->cgcon_.nx - 
			    4 + (j + (k << 4)) * 10) * 5 - 746] - c->cvar_.rsd[
			    m + (c->cgcon_.nx - 3 + (j + (k << 4)) * 10) * 5 - 
			    746] * 4. + c->cvar_.rsd[m + (c->cgcon_.nx - 2 + (j + 
			    (k << 4)) * 10) * 5 - 746] * 6. - c->cvar_.rsd[m + (
			    c->cgcon_.nx - 1 + (j + (k << 4)) * 10) * 5 - 746] *
			     4.);
		    c->cvar_.frct[m + (c->cgcon_.nx - 1 + (j + (k << 4)) * 10) * 
			    5 - 746] -= dsspm * (c->cvar_.rsd[m + (c->cgcon_.nx - 
			    3 + (j + (k << 4)) * 10) * 5 - 746] - c->cvar_.rsd[
			    m + (c->cgcon_.nx - 2 + (j + (k << 4)) * 10) * 5 - 
			    746] * 4. + c->cvar_.rsd[m + (c->cgcon_.nx - 1 + (j + 
			    (k << 4)) * 10) * 5 - 746] * 5.);
		}
	    }
	}
    }

    iex = 1;
    exchange_3__(c->cvar_.rsd, &iex, c);
    l1 = 0;
    if (c->neigh_.west == -1) {
	l1 = 1;
    }
    l2 = c->cgcon_.ny + 1;
    if (c->neigh_.east == -1) {
	l2 = c->cgcon_.ny;
    }
    jst1 = 1;
    jend1 = c->cgcon_.ny;
    if (c->neigh_.west == -1) {
	jst1 = 4;
    }
    if (c->neigh_.east == -1) {
	jend1 = c->cgcon_.ny - 3;
    }
    i__1 = c->cgcon_.nz - 1;
    for (k = 2; k <= i__1; ++k) {
	i__2 = l2;
	for (j = l1; j <= i__2; ++j) {
	    i__3 = c->cgcon_.iend;
	    for (i__ = c->cgcon_.ist; i__ <= i__3; ++i__) {
		c->cvar_.flux[(i__ + ((j + k * 14) << 3)) * 5 - 560] = c->cvar_.rsd[
			(i__ + (j + (k << 4)) * 10) * 5 - 743];
		u31 = c->cvar_.rsd[(i__ + (j + (k << 4)) * 10) * 5 - 743] / 
			c->cvar_.rsd[(i__ + (j + (k << 4)) * 10) * 5 - 745];
		q = (c->cvar_.rsd[(i__ + (j + (k << 4)) * 10) * 5 - 744] * 
			c->cvar_.rsd[(i__ + (j + (k << 4)) * 10) * 5 - 744] + 
			c->cvar_.rsd[(i__ + (j + (k << 4)) * 10) * 5 - 743] * 
			c->cvar_.rsd[(i__ + (j + (k << 4)) * 10) * 5 - 743] + 
			c->cvar_.rsd[(i__ + (j + (k << 4)) * 10) * 5 - 742] * 
			c->cvar_.rsd[(i__ + (j + (k << 4)) * 10) * 5 - 742]) * 
			.5 / c->cvar_.rsd[(i__ + (j + (k << 4)) * 10) * 5 - 745]
			;
		c->cvar_.flux[(i__ + ((j + k * 14) << 3)) * 5 - 559] = c->cvar_.rsd[
			(i__ + (j + (k << 4)) * 10) * 5 - 744] * u31;
		c->cvar_.flux[(i__ + ((j + k * 14) << 3)) * 5 - 558] = c->cvar_.rsd[
			(i__ + (j + (k << 4)) * 10) * 5 - 743] * u31 + (
			c->cvar_.rsd[(i__ + (j + (k << 4)) * 10) * 5 - 741] - q)
			 * .4;
		c->cvar_.flux[(i__ + ((j + k * 14) << 3)) * 5 - 557] = c->cvar_.rsd[
			(i__ + (j + (k << 4)) * 10) * 5 - 742] * u31;
		c->cvar_.flux[(i__ + ((j + k * 14) << 3)) * 5 - 556] = (
			c->cvar_.rsd[(i__ + (j + (k << 4)) * 10) * 5 - 741] * 
			1.4 - q * .4) * u31;
	    }
	}
    }
    i__1 = c->cgcon_.nz - 1;
    for (k = 2; k <= i__1; ++k) {
	i__2 = c->cgcon_.iend;
	for (i__ = c->cgcon_.ist; i__ <= i__2; ++i__) {
	    i__3 = c->cgcon_.jend;
	    for (j = c->cgcon_.jst; j <= i__3; ++j) {
		for (m = 1; m <= 5; ++m) {
		    c->cvar_.frct[m + (i__ + (j + (k << 4)) * 10) * 5 - 746] -= 
			    c->cgcon_.ty2 * (c->cvar_.flux[m + (i__ + ((j + 1 + k *
			     14) << 3)) * 5 - 561] - c->cvar_.flux[m + (i__ + ((j 
			    - 1 + k * 14) << 3)) * 5 - 561]);
		}
	    }
	}
	i__2 = l2;
	for (j = c->cgcon_.jst; j <= i__2; ++j) {
	    i__3 = c->cgcon_.iend;
	    for (i__ = c->cgcon_.ist; i__ <= i__3; ++i__) {
		tmp = 1. / c->cvar_.rsd[(i__ + (j + (k << 4)) * 10) * 5 - 745];
		u21j = tmp * c->cvar_.rsd[(i__ + (j + (k << 4)) * 10) * 5 - 744]
			;
		u31j = tmp * c->cvar_.rsd[(i__ + (j + (k << 4)) * 10) * 5 - 743]
			;
		u41j = tmp * c->cvar_.rsd[(i__ + (j + (k << 4)) * 10) * 5 - 742]
			;
		u51j = tmp * c->cvar_.rsd[(i__ + (j + (k << 4)) * 10) * 5 - 741]
			;
		tmp = 1. / c->cvar_.rsd[(i__ + (j - 1 + (k << 4)) * 10) * 5 - 
			745];
		u21jm1 = tmp * c->cvar_.rsd[(i__ + (j - 1 + (k << 4)) * 10) * 5 
			- 744];
		u31jm1 = tmp * c->cvar_.rsd[(i__ + (j - 1 + (k << 4)) * 10) * 5 
			- 743];
		u41jm1 = tmp * c->cvar_.rsd[(i__ + (j - 1 + (k << 4)) * 10) * 5 
			- 742];
		u51jm1 = tmp * c->cvar_.rsd[(i__ + (j - 1 + (k << 4)) * 10) * 5 
			- 741];
		c->cvar_.flux[(i__ + ((j + k * 14) << 3)) * 5 - 559] = 
			c->cgcon_.ty3 * (u21j - u21jm1);
		c->cvar_.flux[(i__ + ((j + k * 14) << 3)) * 5 - 558] = 
			c->cgcon_.ty3 * 1.3333333333333333 * (u31j - u31jm1);
		c->cvar_.flux[(i__ + ((j + k * 14) << 3)) * 5 - 557] = 
			c->cgcon_.ty3 * (u41j - u41jm1);
		d__1 = u21j;
		d__2 = u31j;
		d__3 = u41j;
		d__4 = u21jm1;
		d__5 = u31jm1;
		d__6 = u41jm1;
		d__7 = u31j;
		d__8 = u31jm1;
		c->cvar_.flux[(i__ + ((j + k * 14) << 3)) * 5 - 556] = 
			c->cgcon_.ty3 * -.47999999999999987 * (d__1 * d__1 + 
			d__2 * d__2 + d__3 * d__3 - (d__4 * d__4 + d__5 * 
			d__5 + d__6 * d__6)) + c->cgcon_.ty3 * 
			.16666666666666666 * (d__7 * d__7 - d__8 * d__8) + 
			c->cgcon_.ty3 * 1.9599999999999997 * (u51j - u51jm1);
	    }
	}
	i__2 = c->cgcon_.jend;
	for (j = c->cgcon_.jst; j <= i__2; ++j) {
	    i__3 = c->cgcon_.iend;
	    for (i__ = c->cgcon_.ist; i__ <= i__3; ++i__) {
		c->cvar_.frct[(i__ + (j + (k << 4)) * 10) * 5 - 745] += 
			c->disp_.dy1 * c->cgcon_.ty1 * (c->cvar_.rsd[(i__ + (j - 1 
			+ (k << 4)) * 10) * 5 - 745] - c->cvar_.rsd[(i__ + (j + 
			(k << 4)) * 10) * 5 - 745] * 2. + c->cvar_.rsd[(i__ + (
			j + 1 + (k << 4)) * 10) * 5 - 745]);
		c->cvar_.frct[(i__ + (j + (k << 4)) * 10) * 5 - 744] = 
			c->cvar_.frct[(i__ + (j + (k << 4)) * 10) * 5 - 744] + 
			c->cgcon_.ty3 * .1 * 1. * (c->cvar_.flux[(i__ + ((j + 1 + 
			k * 14) << 3)) * 5 - 559] - c->cvar_.flux[(i__ + ((j + k *
			 14) << 3)) * 5 - 559]) + c->disp_.dy2 * c->cgcon_.ty1 * (
			c->cvar_.rsd[(i__ + (j - 1 + (k << 4)) * 10) * 5 - 744] 
			- c->cvar_.rsd[(i__ + (j + (k << 4)) * 10) * 5 - 744] * 
			2. + c->cvar_.rsd[(i__ + (j + 1 + (k << 4)) * 10) * 5 - 
			744]);
		c->cvar_.frct[(i__ + (j + (k << 4)) * 10) * 5 - 743] = 
			c->cvar_.frct[(i__ + (j + (k << 4)) * 10) * 5 - 743] + 
			c->cgcon_.ty3 * .1 * 1. * (c->cvar_.flux[(i__ + ((j + 1 + 
			k * 14) << 3)) * 5 - 558] - c->cvar_.flux[(i__ + ((j + k *
			 14) << 3)) * 5 - 558]) + c->disp_.dy3 * c->cgcon_.ty1 * (
			c->cvar_.rsd[(i__ + (j - 1 + (k << 4)) * 10) * 5 - 743] 
			- c->cvar_.rsd[(i__ + (j + (k << 4)) * 10) * 5 - 743] * 
			2. + c->cvar_.rsd[(i__ + (j + 1 + (k << 4)) * 10) * 5 - 
			743]);
		c->cvar_.frct[(i__ + (j + (k << 4)) * 10) * 5 - 742] = 
			c->cvar_.frct[(i__ + (j + (k << 4)) * 10) * 5 - 742] + 
			c->cgcon_.ty3 * .1 * 1. * (c->cvar_.flux[(i__ + ((j + 1 + 
			k * 14) << 3)) * 5 - 557] - c->cvar_.flux[(i__ + ((j + k *
			 14) << 3)) * 5 - 557]) + c->disp_.dy4 * c->cgcon_.ty1 * (
			c->cvar_.rsd[(i__ + (j - 1 + (k << 4)) * 10) * 5 - 742] 
			- c->cvar_.rsd[(i__ + (j + (k << 4)) * 10) * 5 - 742] * 
			2. + c->cvar_.rsd[(i__ + (j + 1 + (k << 4)) * 10) * 5 - 
			742]);
		c->cvar_.frct[(i__ + (j + (k << 4)) * 10) * 5 - 741] = 
			c->cvar_.frct[(i__ + (j + (k << 4)) * 10) * 5 - 741] + 
			c->cgcon_.ty3 * .1 * 1. * (c->cvar_.flux[(i__ + ((j + 1 + 
			k * 14) << 3)) * 5 - 556] - c->cvar_.flux[(i__ + ((j + k *
			 14) << 3)) * 5 - 556]) + c->disp_.dy5 * c->cgcon_.ty1 * (
			c->cvar_.rsd[(i__ + (j - 1 + (k << 4)) * 10) * 5 - 741] 
			- c->cvar_.rsd[(i__ + (j + (k << 4)) * 10) * 5 - 741] * 
			2. + c->cvar_.rsd[(i__ + (j + 1 + (k << 4)) * 10) * 5 - 
			741]);
	    }
	}
	if (c->neigh_.west == -1) {
	    i__2 = c->cgcon_.iend;
	    for (i__ = c->cgcon_.ist; i__ <= i__2; ++i__) {
		for (m = 1; m <= 5; ++m) {
		    c->cvar_.frct[m + (i__ + ((k << 4) + 2) * 10) * 5 - 746] -= 
			    dsspm * (c->cvar_.rsd[m + (i__ + ((k << 4) + 2) * 
			    10) * 5 - 746] * 5. - c->cvar_.rsd[m + (i__ + ((k <<
			     4) + 3) * 10) * 5 - 746] * 4. + c->cvar_.rsd[m + (
			    i__ + ((k << 4) + 4) * 10) * 5 - 746]);
		    c->cvar_.frct[m + (i__ + ((k << 4) + 3) * 10) * 5 - 746] -= 
			    dsspm * (c->cvar_.rsd[m + (i__ + ((k << 4) + 2) * 
			    10) * 5 - 746] * -4. + c->cvar_.rsd[m + (i__ + ((k 
			    << 4) + 3) * 10) * 5 - 746] * 6. - c->cvar_.rsd[m + 
			    (i__ + ((k << 4) + 4) * 10) * 5 - 746] * 4. + 
			    c->cvar_.rsd[m + (i__ + ((k << 4) + 5) * 10) * 5 - 
			    746]);
		}
	    }
	}
	i__2 = jend1;
	for (j = jst1; j <= i__2; ++j) {
	    i__3 = c->cgcon_.iend;
	    for (i__ = c->cgcon_.ist; i__ <= i__3; ++i__) {
		for (m = 1; m <= 5; ++m) {
		    c->cvar_.frct[m + (i__ + (j + (k << 4)) * 10) * 5 - 746] -= 
			    dsspm * (c->cvar_.rsd[m + (i__ + (j - 2 + (k << 4)) 
			    * 10) * 5 - 746] - c->cvar_.rsd[m + (i__ + (j - 1 + 
			    (k << 4)) * 10) * 5 - 746] * 4. + c->cvar_.rsd[m + (
			    i__ + (j + (k << 4)) * 10) * 5 - 746] * 6. - 
			    c->cvar_.rsd[m + (i__ + (j + 1 + (k << 4)) * 10) * 
			    5 - 746] * 4. + c->cvar_.rsd[m + (i__ + (j + 2 + (k 
			    << 4)) * 10) * 5 - 746]);
		}
	    }
	}
	if (c->neigh_.east == -1) {
	    i__2 = c->cgcon_.iend;
	    for (i__ = c->cgcon_.ist; i__ <= i__2; ++i__) {
		for (m = 1; m <= 5; ++m) {
		    c->cvar_.frct[m + (i__ + (c->cgcon_.ny - 2 + (k << 4)) * 10) *
			     5 - 746] -= dsspm * (c->cvar_.rsd[m + (i__ + (
			    c->cgcon_.ny - 4 + (k << 4)) * 10) * 5 - 746] - 
			    c->cvar_.rsd[m + (i__ + (c->cgcon_.ny - 3 + (k << 4)) 
			    * 10) * 5 - 746] * 4. + c->cvar_.rsd[m + (i__ + (
			    c->cgcon_.ny - 2 + (k << 4)) * 10) * 5 - 746] * 6. 
			    - c->cvar_.rsd[m + (i__ + (c->cgcon_.ny - 1 + (k << 4)
			    ) * 10) * 5 - 746] * 4.);
		    c->cvar_.frct[m + (i__ + (c->cgcon_.ny - 1 + (k << 4)) * 10) *
			     5 - 746] -= dsspm * (c->cvar_.rsd[m + (i__ + (
			    c->cgcon_.ny - 3 + (k << 4)) * 10) * 5 - 746] - 
			    c->cvar_.rsd[m + (i__ + (c->cgcon_.ny - 2 + (k << 4)) 
			    * 10) * 5 - 746] * 4. + c->cvar_.rsd[m + (i__ + (
			    c->cgcon_.ny - 1 + (k << 4)) * 10) * 5 - 746] * 5.);
		}
	    }
	}
    }
    i__1 = c->cgcon_.nz;
    for (k = 1; k <= i__1; ++k) {
	i__2 = c->cgcon_.jend;
	for (j = c->cgcon_.jst; j <= i__2; ++j) {
	    i__3 = c->cgcon_.iend;
	    for (i__ = c->cgcon_.ist; i__ <= i__3; ++i__) {
		c->cvar_.flux[(i__ + ((j + k * 14) << 3)) * 5 - 560] = c->cvar_.rsd[
			(i__ + (j + (k << 4)) * 10) * 5 - 742];
		u41 = c->cvar_.rsd[(i__ + (j + (k << 4)) * 10) * 5 - 742] / 
			c->cvar_.rsd[(i__ + (j + (k << 4)) * 10) * 5 - 745];
		q = (c->cvar_.rsd[(i__ + (j + (k << 4)) * 10) * 5 - 744] * 
			c->cvar_.rsd[(i__ + (j + (k << 4)) * 10) * 5 - 744] + 
			c->cvar_.rsd[(i__ + (j + (k << 4)) * 10) * 5 - 743] * 
			c->cvar_.rsd[(i__ + (j + (k << 4)) * 10) * 5 - 743] + 
			c->cvar_.rsd[(i__ + (j + (k << 4)) * 10) * 5 - 742] * 
			c->cvar_.rsd[(i__ + (j + (k << 4)) * 10) * 5 - 742]) * 
			.5 / c->cvar_.rsd[(i__ + (j + (k << 4)) * 10) * 5 - 745]
			;
		c->cvar_.flux[(i__ + ((j + k * 14) << 3)) * 5 - 559] = c->cvar_.rsd[
			(i__ + (j + (k << 4)) * 10) * 5 - 744] * u41;
		c->cvar_.flux[(i__ + ((j + k * 14) << 3)) * 5 - 558] = c->cvar_.rsd[
			(i__ + (j + (k << 4)) * 10) * 5 - 743] * u41;
		c->cvar_.flux[(i__ + ((j + k * 14) << 3)) * 5 - 557] = c->cvar_.rsd[
			(i__ + (j + (k << 4)) * 10) * 5 - 742] * u41 + (
			c->cvar_.rsd[(i__ + (j + (k << 4)) * 10) * 5 - 741] - q)
			 * .4;
		c->cvar_.flux[(i__ + ((j + k * 14) << 3)) * 5 - 556] = (
			c->cvar_.rsd[(i__ + (j + (k << 4)) * 10) * 5 - 741] * 
			1.4 - q * .4) * u41;
	    }
	}
    }
    i__1 = c->cgcon_.nz - 1;
    for (k = 2; k <= i__1; ++k) {
	i__2 = c->cgcon_.jend;
	for (j = c->cgcon_.jst; j <= i__2; ++j) {
	    i__3 = c->cgcon_.iend;
	    for (i__ = c->cgcon_.ist; i__ <= i__3; ++i__) {
		for (m = 1; m <= 5; ++m) {
		    c->cvar_.frct[m + (i__ + (j + (k << 4)) * 10) * 5 - 746] -= 
			    c->cgcon_.tz2 * (c->cvar_.flux[m + (i__ + ((j + (k + 1)
			     * 14) << 3)) * 5 - 561] - c->cvar_.flux[m + (i__ + (
			    (j + (k - 1) * 14) << 3)) * 5 - 561]);
		}
	    }
	}
    }
    i__1 = c->cgcon_.nz;
    for (k = 2; k <= i__1; ++k) {
	i__2 = c->cgcon_.jend;
	for (j = c->cgcon_.jst; j <= i__2; ++j) {
	    i__3 = c->cgcon_.iend;
	    for (i__ = c->cgcon_.ist; i__ <= i__3; ++i__) {
		tmp = 1. / c->cvar_.rsd[(i__ + (j + (k << 4)) * 10) * 5 - 745];
		u21k = tmp * c->cvar_.rsd[(i__ + (j + (k << 4)) * 10) * 5 - 744]
			;
		u31k = tmp * c->cvar_.rsd[(i__ + (j + (k << 4)) * 10) * 5 - 743]
			;
		u41k = tmp * c->cvar_.rsd[(i__ + (j + (k << 4)) * 10) * 5 - 742]
			;
		u51k = tmp * c->cvar_.rsd[(i__ + (j + (k << 4)) * 10) * 5 - 741]
			;
		tmp = 1. / c->cvar_.rsd[(i__ + (j + ((k - 1) << 4)) * 10) * 5 - 
			745];
		u21km1 = tmp * c->cvar_.rsd[(i__ + (j + ((k - 1) << 4)) * 10) * 5 
			- 744];
		u31km1 = tmp * c->cvar_.rsd[(i__ + (j + ((k - 1) << 4)) * 10) * 5 
			- 743];
		u41km1 = tmp * c->cvar_.rsd[(i__ + (j + ((k - 1) << 4)) * 10) * 5 
			- 742];
		u51km1 = tmp * c->cvar_.rsd[(i__ + (j + ((k - 1) << 4)) * 10) * 5 
			- 741];
		c->cvar_.flux[(i__ + ((j + k * 14) << 3)) * 5 - 559] = 
			c->cgcon_.tz3 * (u21k - u21km1);
		c->cvar_.flux[(i__ + ((j + k * 14) << 3)) * 5 - 558] = 
			c->cgcon_.tz3 * (u31k - u31km1);
		c->cvar_.flux[(i__ + ((j + k * 14) << 3)) * 5 - 557] = 
			c->cgcon_.tz3 * 1.3333333333333333 * (u41k - u41km1);
		d__1 = u21k;
		d__2 = u31k;
		d__3 = u41k;
		d__4 = u21km1;
		d__5 = u31km1;
		d__6 = u41km1;
		d__7 = u41k;
		d__8 = u41km1;
		c->cvar_.flux[(i__ + ((j + k * 14) << 3)) * 5 - 556] = 
			c->cgcon_.tz3 * -.47999999999999987 * (d__1 * d__1 + 
			d__2 * d__2 + d__3 * d__3 - (d__4 * d__4 + d__5 * 
			d__5 + d__6 * d__6)) + c->cgcon_.tz3 * 
			.16666666666666666 * (d__7 * d__7 - d__8 * d__8) + 
			c->cgcon_.tz3 * 1.9599999999999997 * (u51k - u51km1);
	    }
	}
    }
    i__1 = c->cgcon_.nz - 1;
    for (k = 2; k <= i__1; ++k) {
	i__2 = c->cgcon_.jend;
	for (j = c->cgcon_.jst; j <= i__2; ++j) {
	    i__3 = c->cgcon_.iend;
	    for (i__ = c->cgcon_.ist; i__ <= i__3; ++i__) {
		c->cvar_.frct[(i__ + (j + (k << 4)) * 10) * 5 - 745] += 
			c->disp_.dz1 * c->cgcon_.tz1 * (c->cvar_.rsd[(i__ + (j + ((k 
			+ 1) << 4)) * 10) * 5 - 745] - c->cvar_.rsd[(i__ + (j + (
			k << 4)) * 10) * 5 - 745] * 2. + c->cvar_.rsd[(i__ + (j 
			+ ((k - 1) << 4)) * 10) * 5 - 745]);
		c->cvar_.frct[(i__ + (j + (k << 4)) * 10) * 5 - 744] = 
			c->cvar_.frct[(i__ + (j + (k << 4)) * 10) * 5 - 744] + 
			c->cgcon_.tz3 * .1 * 1. * (c->cvar_.flux[(i__ + 
                              ((j + (k + 1) * 14) << 3)) * 5 - 559] - c->cvar_.flux[(i__ + 
                              ((j + k * 14) << 3)) * 5 - 559]) + c->disp_.dz2 * c->cgcon_.tz1 * (
			c->cvar_.rsd[(i__ + (j + ((k + 1) << 4)) * 10) * 5 - 744] 
			- c->cvar_.rsd[(i__ + (j + (k << 4)) * 10) * 5 - 744] * 
			2. + c->cvar_.rsd[(i__ + (j + ((k - 1) << 4)) * 10) * 5 - 
			744]);
		c->cvar_.frct[(i__ + (j + (k << 4)) * 10) * 5 - 743] = 
			c->cvar_.frct[(i__ + (j + (k << 4)) * 10) * 5 - 743] + 
			c->cgcon_.tz3 * .1 * 1. * (c->cvar_.flux[(i__ + 
                              ((j + (k + 1) * 14) << 3)) * 5 - 558] - c->cvar_.flux[(i__ + 
                              ((j + k * 14) << 3)) * 5 - 558]) + c->disp_.dz3 * c->cgcon_.tz1 * (
			c->cvar_.rsd[(i__ + (j + ((k + 1) << 4)) * 10) * 5 - 743] 
			- c->cvar_.rsd[(i__ + (j + (k << 4)) * 10) * 5 - 743] * 
			2. + c->cvar_.rsd[(i__ + (j + ((k - 1) << 4)) * 10) * 5 - 
			743]);
		c->cvar_.frct[(i__ + (j + (k << 4)) * 10) * 5 - 742] = 
			c->cvar_.frct[(i__ + (j + (k << 4)) * 10) * 5 - 742] + 
			c->cgcon_.tz3 * .1 * 1. * (c->cvar_.flux[(i__ + 
                              ((j + (k + 1) * 14) << 3)) * 5 - 557] - c->cvar_.flux[(i__ + 
                              ((j + k * 14) << 3)) * 5 - 557]) + c->disp_.dz4 * c->cgcon_.tz1 * (
			c->cvar_.rsd[(i__ + (j + ((k + 1) << 4)) * 10) * 5 - 742] 
			- c->cvar_.rsd[(i__ + (j + (k << 4)) * 10) * 5 - 742] * 
			2. + c->cvar_.rsd[(i__ + (j + ((k - 1) << 4)) * 10) * 5 - 
			742]);
		c->cvar_.frct[(i__ + (j + (k << 4)) * 10) * 5 - 741] = 
			c->cvar_.frct[(i__ + (j + (k << 4)) * 10) * 5 - 741] + 
			c->cgcon_.tz3 * .1 * 1. * (c->cvar_.flux[(i__ + 
                              ((j + (k + 1) * 14) << 3)) * 5 - 556] - c->cvar_.flux[(i__ + 
                              ((j + k * 14) << 3)) * 5 - 556]) + c->disp_.dz5 * c->cgcon_.tz1 * (
			c->cvar_.rsd[(i__ + (j + ((k + 1) << 4)) * 10) * 5 - 741] 
			- c->cvar_.rsd[(i__ + (j + (k << 4)) * 10) * 5 - 741] * 
			2. + c->cvar_.rsd[(i__ + (j + ((k - 1) << 4)) * 10) * 5 - 
			741]);
	    }
	}
    }
    i__1 = c->cgcon_.jend;
    for (j = c->cgcon_.jst; j <= i__1; ++j) {
	i__2 = c->cgcon_.iend;
	for (i__ = c->cgcon_.ist; i__ <= i__2; ++i__) {
	    for (m = 1; m <= 5; ++m) {
		c->cvar_.frct[m + (i__ + (j + 32) * 10) * 5 - 746] -= dsspm * (
			c->cvar_.rsd[m + (i__ + (j + 32) * 10) * 5 - 746] * 5. 
			- c->cvar_.rsd[m + (i__ + (j + 48) * 10) * 5 - 746] * 
			4. + c->cvar_.rsd[m + (i__ + (j + 64) * 10) * 5 - 746]);
		c->cvar_.frct[m + (i__ + (j + 48) * 10) * 5 - 746] -= dsspm * (
			c->cvar_.rsd[m + (i__ + (j + 32) * 10) * 5 - 746] * -4. 
			+ c->cvar_.rsd[m + (i__ + (j + 48) * 10) * 5 - 746] * 
			6. - c->cvar_.rsd[m + (i__ + (j + 64) * 10) * 5 - 746] *
			 4. + c->cvar_.rsd[m + (i__ + (j + 80) * 10) * 5 - 746])
			;
	    }
	}
    }
    i__1 = c->cgcon_.nz - 3;
    for (k = 4; k <= i__1; ++k) {
	i__2 = c->cgcon_.jend;
	for (j = c->cgcon_.jst; j <= i__2; ++j) {
	    i__3 = c->cgcon_.iend;
	    for (i__ = c->cgcon_.ist; i__ <= i__3; ++i__) {
		for (m = 1; m <= 5; ++m) {
		    c->cvar_.frct[m + (i__ + (j + (k << 4)) * 10) * 5 - 746] -= 
			    dsspm * (c->cvar_.rsd[m + (i__ + (j + ((k - 2) << 4)) 
			    * 10) * 5 - 746] - c->cvar_.rsd[m + (i__ + (j + ((k 
			    - 1) << 4)) * 10) * 5 - 746] * 4. + c->cvar_.rsd[m + 
			    (i__ + (j + (k << 4)) * 10) * 5 - 746] * 6. - 
			    c->cvar_.rsd[m + (i__ + (j + ((k + 1) << 4)) * 10) * 
			    5 - 746] * 4. + c->cvar_.rsd[m + (i__ + (j + ((k + 2)
			    << 4)) * 10) * 5 - 746]);
		}
	    }
	}
    }
    i__1 = c->cgcon_.jend;
    for (j = c->cgcon_.jst; j <= i__1; ++j) {
	i__2 = c->cgcon_.iend;
	for (i__ = c->cgcon_.ist; i__ <= i__2; ++i__) {
	    for (m = 1; m <= 5; ++m) {
		c->cvar_.frct[m + (i__ + (j + ((c->cgcon_.nz - 2) << 4)) * 10) * 5 
			- 746] -= dsspm * (c->cvar_.rsd[m + (i__ + (j + (
			(c->cgcon_.nz - 4) << 4)) * 10) * 5 - 746] - c->cvar_.rsd[m 
			+ (i__ + (j + ((c->cgcon_.nz - 3) << 4)) * 10) * 5 - 746] 
			* 4. + c->cvar_.rsd[m + (i__ + (j + ((c->cgcon_.nz - 2) << 
			4)) * 10) * 5 - 746] * 6. - c->cvar_.rsd[m + (i__ + (j 
			+ ((c->cgcon_.nz - 1) << 4)) * 10) * 5 - 746] * 4.);
		c->cvar_.frct[m + (i__ + (j + ((c->cgcon_.nz - 1) << 4)) * 10) * 5 
			- 746] -= dsspm * (c->cvar_.rsd[m + (i__ + (j + (
			(c->cgcon_.nz - 3) << 4)) * 10) * 5 - 746] - c->cvar_.rsd[m 
			+ (i__ + (j + ((c->cgcon_.nz - 2) << 4)) * 10) * 5 - 746] 
			* 4. + c->cvar_.rsd[m + (i__ + (j + ((c->cgcon_.nz - 1) << 
			4)) * 10) * 5 - 746] * 5.);
	    }
	}
    }
    return 0;
} /* erhs_ */

/* Subroutine */ int buts_(integer *ldmx, integer *ldmy, integer *ldmz, 
	integer *nx, integer *ny, integer *nz, integer *k, float *omega, 
	float *v, float *tv, float *d__, float *udx, 
	float *udy, float *udz, integer *ist, integer *iend, 
	integer *jst, integer *jend, integer *nx0, integer *ny0, integer *ipt,
	 integer *jpt, context *c)
{
    integer v_dim2, v_dim3, v_offset, tv_dim2, tv_offset, d_dim3, d_offset, 
	    udx_dim3, udx_offset, udy_dim3, udy_offset, udz_dim3, udz_offset, 
	    i__1, i__2;

    integer i__, j, m;
    integer iex;
    float tmp, tmp1;
    float tmat[25]	/* was [5][5] */;

    udz_dim3 = *ldmx;
    udz_offset = 1 + 5 * (1 + 5 * (1 + udz_dim3));
    udz -= udz_offset;
    udy_dim3 = *ldmx;
    udy_offset = 1 + 5 * (1 + 5 * (1 + udy_dim3));
    udy -= udy_offset;
    udx_dim3 = *ldmx;
    udx_offset = 1 + 5 * (1 + 5 * (1 + udx_dim3));
    udx -= udx_offset;
    d_dim3 = *ldmx;
    d_offset = 1 + 5 * (1 + 5 * (1 + d_dim3));
    d__ -= d_offset;
    tv_dim2 = *ldmx;
    tv_offset = 1 + 5 * (1 + tv_dim2);
    tv -= tv_offset;
    v_dim2 = *ldmx + 2 + 1 + 1;
    v_dim3 = *ldmy + 2 + 1 + 1;
    v_offset = 1 + 5 * (-1 + v_dim2 * (-1 + v_dim3));
    v -= v_offset;

    if (c->timer_.timeron) {
	timer_start(9, c);
    }
    iex = 1;
    exchange_1__(&v[v_offset], k, &iex, c);
    if (c->timer_.timeron) {
	timer_stop(9, c);
    }
    if (c->timer_.timeron) {
	timer_start(4, c);
    }
    i__1 = *jst;
    for (j = *jend; j >= i__1; --j) {
	i__2 = *ist;
	for (i__ = *iend; i__ >= i__2; --i__) {
	    for (m = 1; m <= 5; ++m) {
		tv[m + (i__ + j * tv_dim2) * 5] = *omega * (udz[m + ((i__ + j 
			* udz_dim3) * 5 + 1) * 5] * v[(i__ + (j + (*k + 1) * 
			v_dim3) * v_dim2) * 5 + 1] + udz[m + ((i__ + j * 
			udz_dim3) * 5 + 2) * 5] * v[(i__ + (j + (*k + 1) * 
			v_dim3) * v_dim2) * 5 + 2] + udz[m + ((i__ + j * 
			udz_dim3) * 5 + 3) * 5] * v[(i__ + (j + (*k + 1) * 
			v_dim3) * v_dim2) * 5 + 3] + udz[m + ((i__ + j * 
			udz_dim3) * 5 + 4) * 5] * v[(i__ + (j + (*k + 1) * 
			v_dim3) * v_dim2) * 5 + 4] + udz[m + ((i__ + j * 
			udz_dim3) * 5 + 5) * 5] * v[(i__ + (j + (*k + 1) * 
			v_dim3) * v_dim2) * 5 + 5]);
	    }
	}
    }
    i__1 = *jst;
    for (j = *jend; j >= i__1; --j) {
	i__2 = *ist;
	for (i__ = *iend; i__ >= i__2; --i__) {
	    for (m = 1; m <= 5; ++m) {
		tv[m + (i__ + j * tv_dim2) * 5] += *omega * (udy[m + ((i__ + 
			j * udy_dim3) * 5 + 1) * 5] * v[(i__ + (j + 1 + *k * 
			v_dim3) * v_dim2) * 5 + 1] + udx[m + ((i__ + j * 
			udx_dim3) * 5 + 1) * 5] * v[(i__ + 1 + (j + *k * 
			v_dim3) * v_dim2) * 5 + 1] + udy[m + ((i__ + j * 
			udy_dim3) * 5 + 2) * 5] * v[(i__ + (j + 1 + *k * 
			v_dim3) * v_dim2) * 5 + 2] + udx[m + ((i__ + j * 
			udx_dim3) * 5 + 2) * 5] * v[(i__ + 1 + (j + *k * 
			v_dim3) * v_dim2) * 5 + 2] + udy[m + ((i__ + j * 
			udy_dim3) * 5 + 3) * 5] * v[(i__ + (j + 1 + *k * 
			v_dim3) * v_dim2) * 5 + 3] + udx[m + ((i__ + j * 
			udx_dim3) * 5 + 3) * 5] * v[(i__ + 1 + (j + *k * 
			v_dim3) * v_dim2) * 5 + 3] + udy[m + ((i__ + j * 
			udy_dim3) * 5 + 4) * 5] * v[(i__ + (j + 1 + *k * 
			v_dim3) * v_dim2) * 5 + 4] + udx[m + ((i__ + j * 
			udx_dim3) * 5 + 4) * 5] * v[(i__ + 1 + (j + *k * 
			v_dim3) * v_dim2) * 5 + 4] + udy[m + ((i__ + j * 
			udy_dim3) * 5 + 5) * 5] * v[(i__ + (j + 1 + *k * 
			v_dim3) * v_dim2) * 5 + 5] + udx[m + ((i__ + j * 
			udx_dim3) * 5 + 5) * 5] * v[(i__ + 1 + (j + *k * 
			v_dim3) * v_dim2) * 5 + 5]);
	    }
	    for (m = 1; m <= 5; ++m) {
		tmat[m - 1] = d__[m + ((i__ + j * d_dim3) * 5 + 1) * 5];
		tmat[m + 4] = d__[m + ((i__ + j * d_dim3) * 5 + 2) * 5];
		tmat[m + 9] = d__[m + ((i__ + j * d_dim3) * 5 + 3) * 5];
		tmat[m + 14] = d__[m + ((i__ + j * d_dim3) * 5 + 4) * 5];
		tmat[m + 19] = d__[m + ((i__ + j * d_dim3) * 5 + 5) * 5];
	    }
	    tmp1 = 1. / tmat[0];
	    tmp = tmp1 * tmat[1];
	    tmat[6] -= tmp * tmat[5];
	    tmat[11] -= tmp * tmat[10];
	    tmat[16] -= tmp * tmat[15];
	    tmat[21] -= tmp * tmat[20];
	    tv[(i__ + j * tv_dim2) * 5 + 2] -= tv[(i__ + j * tv_dim2) * 5 + 1]
		     * tmp;
	    tmp = tmp1 * tmat[2];
	    tmat[7] -= tmp * tmat[5];
	    tmat[12] -= tmp * tmat[10];
	    tmat[17] -= tmp * tmat[15];
	    tmat[22] -= tmp * tmat[20];
	    tv[(i__ + j * tv_dim2) * 5 + 3] -= tv[(i__ + j * tv_dim2) * 5 + 1]
		     * tmp;
	    tmp = tmp1 * tmat[3];
	    tmat[8] -= tmp * tmat[5];
	    tmat[13] -= tmp * tmat[10];
	    tmat[18] -= tmp * tmat[15];
	    tmat[23] -= tmp * tmat[20];
	    tv[(i__ + j * tv_dim2) * 5 + 4] -= tv[(i__ + j * tv_dim2) * 5 + 1]
		     * tmp;
	    tmp = tmp1 * tmat[4];
	    tmat[9] -= tmp * tmat[5];
	    tmat[14] -= tmp * tmat[10];
	    tmat[19] -= tmp * tmat[15];
	    tmat[24] -= tmp * tmat[20];
	    tv[(i__ + j * tv_dim2) * 5 + 5] -= tv[(i__ + j * tv_dim2) * 5 + 1]
		     * tmp;
	    tmp1 = 1. / tmat[6];
	    tmp = tmp1 * tmat[7];
	    tmat[12] -= tmp * tmat[11];
	    tmat[17] -= tmp * tmat[16];
	    tmat[22] -= tmp * tmat[21];
	    tv[(i__ + j * tv_dim2) * 5 + 3] -= tv[(i__ + j * tv_dim2) * 5 + 2]
		     * tmp;
	    tmp = tmp1 * tmat[8];
	    tmat[13] -= tmp * tmat[11];
	    tmat[18] -= tmp * tmat[16];
	    tmat[23] -= tmp * tmat[21];
	    tv[(i__ + j * tv_dim2) * 5 + 4] -= tv[(i__ + j * tv_dim2) * 5 + 2]
		     * tmp;
	    tmp = tmp1 * tmat[9];
	    tmat[14] -= tmp * tmat[11];
	    tmat[19] -= tmp * tmat[16];
	    tmat[24] -= tmp * tmat[21];
	    tv[(i__ + j * tv_dim2) * 5 + 5] -= tv[(i__ + j * tv_dim2) * 5 + 2]
		     * tmp;
	    tmp1 = 1. / tmat[12];
	    tmp = tmp1 * tmat[13];
	    tmat[18] -= tmp * tmat[17];
	    tmat[23] -= tmp * tmat[22];
	    tv[(i__ + j * tv_dim2) * 5 + 4] -= tv[(i__ + j * tv_dim2) * 5 + 3]
		     * tmp;
	    tmp = tmp1 * tmat[14];
	    tmat[19] -= tmp * tmat[17];
	    tmat[24] -= tmp * tmat[22];
	    tv[(i__ + j * tv_dim2) * 5 + 5] -= tv[(i__ + j * tv_dim2) * 5 + 3]
		     * tmp;
	    tmp1 = 1. / tmat[18];
	    tmp = tmp1 * tmat[19];
	    tmat[24] -= tmp * tmat[23];
	    tv[(i__ + j * tv_dim2) * 5 + 5] -= tv[(i__ + j * tv_dim2) * 5 + 4]
		     * tmp;
	    tv[(i__ + j * tv_dim2) * 5 + 5] /= tmat[24];
	    tv[(i__ + j * tv_dim2) * 5 + 4] -= tmat[23] * tv[(i__ + j * 
		    tv_dim2) * 5 + 5];
	    tv[(i__ + j * tv_dim2) * 5 + 4] /= tmat[18];
	    tv[(i__ + j * tv_dim2) * 5 + 3] = tv[(i__ + j * tv_dim2) * 5 + 3] 
		    - tmat[17] * tv[(i__ + j * tv_dim2) * 5 + 4] - tmat[22] * 
		    tv[(i__ + j * tv_dim2) * 5 + 5];
	    tv[(i__ + j * tv_dim2) * 5 + 3] /= tmat[12];
	    tv[(i__ + j * tv_dim2) * 5 + 2] = tv[(i__ + j * tv_dim2) * 5 + 2] 
		    - tmat[11] * tv[(i__ + j * tv_dim2) * 5 + 3] - tmat[16] * 
		    tv[(i__ + j * tv_dim2) * 5 + 4] - tmat[21] * tv[(i__ + j *
		     tv_dim2) * 5 + 5];
	    tv[(i__ + j * tv_dim2) * 5 + 2] /= tmat[6];
	    tv[(i__ + j * tv_dim2) * 5 + 1] = tv[(i__ + j * tv_dim2) * 5 + 1] 
		    - tmat[5] * tv[(i__ + j * tv_dim2) * 5 + 2] - tmat[10] * 
		    tv[(i__ + j * tv_dim2) * 5 + 3] - tmat[15] * tv[(i__ + j *
		     tv_dim2) * 5 + 4] - tmat[20] * tv[(i__ + j * tv_dim2) * 
		    5 + 5];
	    tv[(i__ + j * tv_dim2) * 5 + 1] /= tmat[0];
	    v[(i__ + (j + *k * v_dim3) * v_dim2) * 5 + 1] -= tv[(i__ + j * 
		    tv_dim2) * 5 + 1];
	    v[(i__ + (j + *k * v_dim3) * v_dim2) * 5 + 2] -= tv[(i__ + j * 
		    tv_dim2) * 5 + 2];
	    v[(i__ + (j + *k * v_dim3) * v_dim2) * 5 + 3] -= tv[(i__ + j * 
		    tv_dim2) * 5 + 3];
	    v[(i__ + (j + *k * v_dim3) * v_dim2) * 5 + 4] -= tv[(i__ + j * 
		    tv_dim2) * 5 + 4];
	    v[(i__ + (j + *k * v_dim3) * v_dim2) * 5 + 5] -= tv[(i__ + j * 
		    tv_dim2) * 5 + 5];
	}
    }
    if (c->timer_.timeron) {
	timer_stop(4, c);
    }
    if (c->timer_.timeron) {
	timer_start(9, c);
    }
    iex = 3;
    exchange_1__(&v[v_offset], k, &iex, c);
    if (c->timer_.timeron) {
	timer_stop(9, c);
    }
    return 0;
} /* buts_ */

/* Subroutine */ int blts_(integer *ldmx, integer *ldmy, integer *ldmz, 
	integer *nx, integer *ny, integer *nz, integer *k, float *omega, 
	float *v, float *ldz, float *ldy, float *ldx, 
	float *d__, integer *ist, integer *iend, integer *jst, integer *
	jend, integer *nx0, integer *ny0, integer *ipt, integer *jpt, context *c)
{
    integer v_dim2, v_dim3, v_offset, ldz_dim3, ldz_offset, ldy_dim3, 
	    ldy_offset, ldx_dim3, ldx_offset, d_dim3, d_offset, i__1, i__2;

    integer i__, j, m;
    integer iex;
    float tmp, tmp1;
    float tmat[25]	/* was [5][5] */;
#define timer_ timer_

    d_dim3 = *ldmx;
    d_offset = 1 + 5 * (1 + 5 * (1 + d_dim3));
    d__ -= d_offset;
    ldx_dim3 = *ldmx;
    ldx_offset = 1 + 5 * (1 + 5 * (1 + ldx_dim3));
    ldx -= ldx_offset;
    ldy_dim3 = *ldmx;
    ldy_offset = 1 + 5 * (1 + 5 * (1 + ldy_dim3));
    ldy -= ldy_offset;
    ldz_dim3 = *ldmx;
    ldz_offset = 1 + 5 * (1 + 5 * (1 + ldz_dim3));
    ldz -= ldz_offset;
    v_dim2 = *ldmx + 2 + 1 + 1;
    v_dim3 = *ldmy + 2 + 1 + 1;
    v_offset = 1 + 5 * (-1 + v_dim2 * (-1 + v_dim3));
    v -= v_offset;

    if (c->timer_.timeron) {
	timer_start(8, c);
    }
    iex = 0;
    exchange_1__(&v[v_offset], k, &iex, c);
    if (c->timer_.timeron) {
	timer_stop(8, c);
    }
    if (c->timer_.timeron) {
	timer_start(3, c);
    }
    i__1 = *jend;
    for (j = *jst; j <= i__1; ++j) {
	i__2 = *iend;
	for (i__ = *ist; i__ <= i__2; ++i__) {
	    for (m = 1; m <= 5; ++m) {
		v[m + (i__ + (j + *k * v_dim3) * v_dim2) * 5] -= *omega * (
			ldz[m + ((i__ + j * ldz_dim3) * 5 + 1) * 5] * v[(i__ 
			+ (j + (*k - 1) * v_dim3) * v_dim2) * 5 + 1] + ldz[m 
			+ ((i__ + j * ldz_dim3) * 5 + 2) * 5] * v[(i__ + (j + 
			(*k - 1) * v_dim3) * v_dim2) * 5 + 2] + ldz[m + ((i__ 
			+ j * ldz_dim3) * 5 + 3) * 5] * v[(i__ + (j + (*k - 1)
			 * v_dim3) * v_dim2) * 5 + 3] + ldz[m + ((i__ + j * 
			ldz_dim3) * 5 + 4) * 5] * v[(i__ + (j + (*k - 1) * 
			v_dim3) * v_dim2) * 5 + 4] + ldz[m + ((i__ + j * 
			ldz_dim3) * 5 + 5) * 5] * v[(i__ + (j + (*k - 1) * 
			v_dim3) * v_dim2) * 5 + 5]);
	    }
	}
    }
    i__1 = *jend;
    for (j = *jst; j <= i__1; ++j) {
	i__2 = *iend;
	for (i__ = *ist; i__ <= i__2; ++i__) {
	    for (m = 1; m <= 5; ++m) {
		v[m + (i__ + (j + *k * v_dim3) * v_dim2) * 5] -= *omega * (
			ldy[m + ((i__ + j * ldy_dim3) * 5 + 1) * 5] * v[(i__ 
			+ (j - 1 + *k * v_dim3) * v_dim2) * 5 + 1] + ldx[m + (
			(i__ + j * ldx_dim3) * 5 + 1) * 5] * v[(i__ - 1 + (j 
			+ *k * v_dim3) * v_dim2) * 5 + 1] + ldy[m + ((i__ + j 
			* ldy_dim3) * 5 + 2) * 5] * v[(i__ + (j - 1 + *k * 
			v_dim3) * v_dim2) * 5 + 2] + ldx[m + ((i__ + j * 
			ldx_dim3) * 5 + 2) * 5] * v[(i__ - 1 + (j + *k * 
			v_dim3) * v_dim2) * 5 + 2] + ldy[m + ((i__ + j * 
			ldy_dim3) * 5 + 3) * 5] * v[(i__ + (j - 1 + *k * 
			v_dim3) * v_dim2) * 5 + 3] + ldx[m + ((i__ + j * 
			ldx_dim3) * 5 + 3) * 5] * v[(i__ - 1 + (j + *k * 
			v_dim3) * v_dim2) * 5 + 3] + ldy[m + ((i__ + j * 
			ldy_dim3) * 5 + 4) * 5] * v[(i__ + (j - 1 + *k * 
			v_dim3) * v_dim2) * 5 + 4] + ldx[m + ((i__ + j * 
			ldx_dim3) * 5 + 4) * 5] * v[(i__ - 1 + (j + *k * 
			v_dim3) * v_dim2) * 5 + 4] + ldy[m + ((i__ + j * 
			ldy_dim3) * 5 + 5) * 5] * v[(i__ + (j - 1 + *k * 
			v_dim3) * v_dim2) * 5 + 5] + ldx[m + ((i__ + j * 
			ldx_dim3) * 5 + 5) * 5] * v[(i__ - 1 + (j + *k * 
			v_dim3) * v_dim2) * 5 + 5]);
	    }

	    for (m = 1; m <= 5; ++m) {
		tmat[m - 1] = d__[m + ((i__ + j * d_dim3) * 5 + 1) * 5];
		tmat[m + 4] = d__[m + ((i__ + j * d_dim3) * 5 + 2) * 5];
		tmat[m + 9] = d__[m + ((i__ + j * d_dim3) * 5 + 3) * 5];
		tmat[m + 14] = d__[m + ((i__ + j * d_dim3) * 5 + 4) * 5];
		tmat[m + 19] = d__[m + ((i__ + j * d_dim3) * 5 + 5) * 5];
	    }
	    tmp1 = 1. / tmat[0];
	    tmp = tmp1 * tmat[1];
	    tmat[6] -= tmp * tmat[5];
	    tmat[11] -= tmp * tmat[10];
	    tmat[16] -= tmp * tmat[15];
	    tmat[21] -= tmp * tmat[20];
	    v[(i__ + (j + *k * v_dim3) * v_dim2) * 5 + 2] -= v[(i__ + (j + *k 
		    * v_dim3) * v_dim2) * 5 + 1] * tmp;
	    tmp = tmp1 * tmat[2];
	    tmat[7] -= tmp * tmat[5];
	    tmat[12] -= tmp * tmat[10];
	    tmat[17] -= tmp * tmat[15];
	    tmat[22] -= tmp * tmat[20];
	    v[(i__ + (j + *k * v_dim3) * v_dim2) * 5 + 3] -= v[(i__ + (j + *k 
		    * v_dim3) * v_dim2) * 5 + 1] * tmp;
	    tmp = tmp1 * tmat[3];
	    tmat[8] -= tmp * tmat[5];
	    tmat[13] -= tmp * tmat[10];
	    tmat[18] -= tmp * tmat[15];
	    tmat[23] -= tmp * tmat[20];
	    v[(i__ + (j + *k * v_dim3) * v_dim2) * 5 + 4] -= v[(i__ + (j + *k 
		    * v_dim3) * v_dim2) * 5 + 1] * tmp;
	    tmp = tmp1 * tmat[4];
	    tmat[9] -= tmp * tmat[5];
	    tmat[14] -= tmp * tmat[10];
	    tmat[19] -= tmp * tmat[15];
	    tmat[24] -= tmp * tmat[20];
	    v[(i__ + (j + *k * v_dim3) * v_dim2) * 5 + 5] -= v[(i__ + (j + *k 
		    * v_dim3) * v_dim2) * 5 + 1] * tmp;
	    tmp1 = 1. / tmat[6];
	    tmp = tmp1 * tmat[7];
	    tmat[12] -= tmp * tmat[11];
	    tmat[17] -= tmp * tmat[16];
	    tmat[22] -= tmp * tmat[21];
	    v[(i__ + (j + *k * v_dim3) * v_dim2) * 5 + 3] -= v[(i__ + (j + *k 
		    * v_dim3) * v_dim2) * 5 + 2] * tmp;
	    tmp = tmp1 * tmat[8];
	    tmat[13] -= tmp * tmat[11];
	    tmat[18] -= tmp * tmat[16];
	    tmat[23] -= tmp * tmat[21];
	    v[(i__ + (j + *k * v_dim3) * v_dim2) * 5 + 4] -= v[(i__ + (j + *k 
		    * v_dim3) * v_dim2) * 5 + 2] * tmp;
	    tmp = tmp1 * tmat[9];
	    tmat[14] -= tmp * tmat[11];
	    tmat[19] -= tmp * tmat[16];
	    tmat[24] -= tmp * tmat[21];
	    v[(i__ + (j + *k * v_dim3) * v_dim2) * 5 + 5] -= v[(i__ + (j + *k 
		    * v_dim3) * v_dim2) * 5 + 2] * tmp;
	    tmp1 = 1. / tmat[12];
	    tmp = tmp1 * tmat[13];
	    tmat[18] -= tmp * tmat[17];
	    tmat[23] -= tmp * tmat[22];
	    v[(i__ + (j + *k * v_dim3) * v_dim2) * 5 + 4] -= v[(i__ + (j + *k 
		    * v_dim3) * v_dim2) * 5 + 3] * tmp;
	    tmp = tmp1 * tmat[14];
	    tmat[19] -= tmp * tmat[17];
	    tmat[24] -= tmp * tmat[22];
	    v[(i__ + (j + *k * v_dim3) * v_dim2) * 5 + 5] -= v[(i__ + (j + *k 
		    * v_dim3) * v_dim2) * 5 + 3] * tmp;
	    tmp1 = 1. / tmat[18];
	    tmp = tmp1 * tmat[19];
	    tmat[24] -= tmp * tmat[23];
	    v[(i__ + (j + *k * v_dim3) * v_dim2) * 5 + 5] -= v[(i__ + (j + *k 
		    * v_dim3) * v_dim2) * 5 + 4] * tmp;
	    v[(i__ + (j + *k * v_dim3) * v_dim2) * 5 + 5] /= tmat[24];
	    v[(i__ + (j + *k * v_dim3) * v_dim2) * 5 + 4] -= tmat[23] * v[(
		    i__ + (j + *k * v_dim3) * v_dim2) * 5 + 5];
	    v[(i__ + (j + *k * v_dim3) * v_dim2) * 5 + 4] /= tmat[18];
	    v[(i__ + (j + *k * v_dim3) * v_dim2) * 5 + 3] = v[(i__ + (j + *k *
		     v_dim3) * v_dim2) * 5 + 3] - tmat[17] * v[(i__ + (j + *k 
		    * v_dim3) * v_dim2) * 5 + 4] - tmat[22] * v[(i__ + (j + *
		    k * v_dim3) * v_dim2) * 5 + 5];
	    v[(i__ + (j + *k * v_dim3) * v_dim2) * 5 + 3] /= tmat[12];
	    v[(i__ + (j + *k * v_dim3) * v_dim2) * 5 + 2] = v[(i__ + (j + *k *
		     v_dim3) * v_dim2) * 5 + 2] - tmat[11] * v[(i__ + (j + *k 
		    * v_dim3) * v_dim2) * 5 + 3] - tmat[16] * v[(i__ + (j + *
		    k * v_dim3) * v_dim2) * 5 + 4] - tmat[21] * v[(i__ + (j + 
		    *k * v_dim3) * v_dim2) * 5 + 5];
	    v[(i__ + (j + *k * v_dim3) * v_dim2) * 5 + 2] /= tmat[6];
	    v[(i__ + (j + *k * v_dim3) * v_dim2) * 5 + 1] = v[(i__ + (j + *k *
		     v_dim3) * v_dim2) * 5 + 1] - tmat[5] * v[(i__ + (j + *k *
		     v_dim3) * v_dim2) * 5 + 2] - tmat[10] * v[(i__ + (j + *k 
		    * v_dim3) * v_dim2) * 5 + 3] - tmat[15] * v[(i__ + (j + *
		    k * v_dim3) * v_dim2) * 5 + 4] - tmat[20] * v[(i__ + (j + 
		    *k * v_dim3) * v_dim2) * 5 + 5];
	    v[(i__ + (j + *k * v_dim3) * v_dim2) * 5 + 1] /= tmat[0];
	}
    }
    if (c->timer_.timeron) {
	timer_stop(3, c);
    }
    if (c->timer_.timeron) {
	timer_start(8, c);
    }
    iex = 2;
    exchange_1__(&v[v_offset], k, &iex, c);

    if (c->timer_.timeron) {
	timer_stop(8, c);
    }
    return 0;
} /* blts_ */

/* Subroutine */ int exchange_1__(float *g, integer *k, integer *iex, context *c)
{
    integer i__1;
    int i;

    integer i__, j;
    float dum[5*(ISIZ1+ISIZ2)]	/* was [5][18] */, dum1[5*(ISIZ1+ISIZ2)]	/* was [5][18]
	     */;
    MPI_Status status[5];

    g -= 746;

    if (*iex == 0) {
	if (c->neigh_.north != -1) {
	    i__1 = (c->cgcon_.jend - c->cgcon_.jst + 1) * 5;
	    MPI_Recv(c->ftmp, i__1, c->mpistuff_.dp_type__, c->neigh_.north, 2, MPI_COMM_WORLD, status);
        for (i = 0; i < i__1; i++) dum1[c->cgcon_.jst * 5 - 5 + i] = c->ftmp[i];
//	    MPI_Recv(&dum1[c->cgcon_.jst * 5 - 5], i__1, c->mpistuff_.dp_type__, c->neigh_.north, 2, MPI_COMM_WORLD, status);
	    i__1 = c->cgcon_.jend;
	    for (j = c->cgcon_.jst; j <= i__1; ++j) {
		g[(j + (*k << 4)) * 50 + 1] = dum1[j * 5 - 5];
		g[(j + (*k << 4)) * 50 + 2] = dum1[j * 5 - 4];
		g[(j + (*k << 4)) * 50 + 3] = dum1[j * 5 - 3];
		g[(j + (*k << 4)) * 50 + 4] = dum1[j * 5 - 2];
		g[(j + (*k << 4)) * 50 + 5] = dum1[j * 5 - 1];
	    }
	}
	if (c->neigh_.west != -1) {
	    i__1 = (c->cgcon_.iend - c->cgcon_.ist + 1) * 5;
	    MPI_Recv(c->ftmp, i__1, c->mpistuff_.dp_type__, c->neigh_.west, 4, MPI_COMM_WORLD, status);
        for (i = 0; i < i__1; i++) dum1[c->cgcon_.ist * 5 - 5 + i] = c->ftmp[i];
//	    MPI_Recv(&dum1[c->cgcon_.ist * 5 - 5], i__1, c->mpistuff_.dp_type__, c->neigh_.west, 4, MPI_COMM_WORLD, status);
	    i__1 = c->cgcon_.iend;
	    for (i__ = c->cgcon_.ist; i__ <= i__1; ++i__) {
		g[(i__ + *k * 160) * 5 + 1] = dum1[i__ * 5 - 5];
		g[(i__ + *k * 160) * 5 + 2] = dum1[i__ * 5 - 4];
		g[(i__ + *k * 160) * 5 + 3] = dum1[i__ * 5 - 3];
		g[(i__ + *k * 160) * 5 + 4] = dum1[i__ * 5 - 2];
		g[(i__ + *k * 160) * 5 + 5] = dum1[i__ * 5 - 1];
	    }
	}
    } else if (*iex == 1) {
	if (c->neigh_.south != -1) {
	    i__1 = (c->cgcon_.jend - c->cgcon_.jst + 1) * 5;
	    MPI_Recv(c->ftmp, i__1, c->mpistuff_.dp_type__, c->neigh_.south, 1, MPI_COMM_WORLD, 		    status);
        for (i = 0; i < i__1; i++) dum1[c->cgcon_.jst * 5 - 5 + i] = c->ftmp[i];
	    i__1 = c->cgcon_.jend;
	    for (j = c->cgcon_.jst; j <= i__1; ++j) {
		g[(c->cgcon_.nx + 1 + (j + (*k << 4)) * 10) * 5 + 1] = dum1[j * 
			5 - 5];
		g[(c->cgcon_.nx + 1 + (j + (*k << 4)) * 10) * 5 + 2] = dum1[j * 
			5 - 4];
		g[(c->cgcon_.nx + 1 + (j + (*k << 4)) * 10) * 5 + 3] = dum1[j * 
			5 - 3];
		g[(c->cgcon_.nx + 1 + (j + (*k << 4)) * 10) * 5 + 4] = dum1[j * 
			5 - 2];
		g[(c->cgcon_.nx + 1 + (j + (*k << 4)) * 10) * 5 + 5] = dum1[j * 
			5 - 1];
	    }
	}
	if (c->neigh_.east != -1) {
	    i__1 = (c->cgcon_.iend - c->cgcon_.ist + 1) * 5;
	    MPI_Recv(c->ftmp, i__1, c->mpistuff_.dp_type__, c->neigh_.east, 3, MPI_COMM_WORLD, status);
        for (i = 0; i < i__1; i++) dum1[c->cgcon_.ist * 5 - 5 + i] = c->ftmp[i];
	    i__1 = c->cgcon_.iend;
	    for (i__ = c->cgcon_.ist; i__ <= i__1; ++i__) {
		g[(i__ + (c->cgcon_.ny + 1 + (*k << 4)) * 10) * 5 + 1] = dum1[
			i__ * 5 - 5];
		g[(i__ + (c->cgcon_.ny + 1 + (*k << 4)) * 10) * 5 + 2] = dum1[
			i__ * 5 - 4];
		g[(i__ + (c->cgcon_.ny + 1 + (*k << 4)) * 10) * 5 + 3] = dum1[
			i__ * 5 - 3];
		g[(i__ + (c->cgcon_.ny + 1 + (*k << 4)) * 10) * 5 + 4] = dum1[
			i__ * 5 - 2];
		g[(i__ + (c->cgcon_.ny + 1 + (*k << 4)) * 10) * 5 + 5] = dum1[
			i__ * 5 - 1];
	    }
	}
    } else if (*iex == 2) {
	if (c->neigh_.south != -1) {
	    i__1 = c->cgcon_.jend;
	    for (j = c->cgcon_.jst; j <= i__1; ++j) {
		dum[j * 5 - 5] = g[(c->cgcon_.nx + (j + (*k << 4)) * 10) * 5 + 
			1];
		dum[j * 5 - 4] = g[(c->cgcon_.nx + (j + (*k << 4)) * 10) * 5 + 
			2];
		dum[j * 5 - 3] = g[(c->cgcon_.nx + (j + (*k << 4)) * 10) * 5 + 
			3];
		dum[j * 5 - 2] = g[(c->cgcon_.nx + (j + (*k << 4)) * 10) * 5 + 
			4];
		dum[j * 5 - 1] = g[(c->cgcon_.nx + (j + (*k << 4)) * 10) * 5 + 
			5];
	    }
	    i__1 = (c->cgcon_.jend - c->cgcon_.jst + 1) * 5;
        for (i = 0; i < i__1; i++) c->ftmp[i] = dum[c->cgcon_.jst * 5 - 5 + i];
	    MPI_Send(c->ftmp, i__1, c->mpistuff_.dp_type__, c->neigh_.south, 2, MPI_COMM_WORLD);
	}
	if (c->neigh_.east != -1) {
	    i__1 = c->cgcon_.iend;
	    for (i__ = c->cgcon_.ist; i__ <= i__1; ++i__) {
		dum[i__ * 5 - 5] = g[(i__ + (c->cgcon_.ny + (*k << 4)) * 10) * 
			5 + 1];
		dum[i__ * 5 - 4] = g[(i__ + (c->cgcon_.ny + (*k << 4)) * 10) * 
			5 + 2];
		dum[i__ * 5 - 3] = g[(i__ + (c->cgcon_.ny + (*k << 4)) * 10) * 
			5 + 3];
		dum[i__ * 5 - 2] = g[(i__ + (c->cgcon_.ny + (*k << 4)) * 10) * 
			5 + 4];
		dum[i__ * 5 - 1] = g[(i__ + (c->cgcon_.ny + (*k << 4)) * 10) * 
			5 + 5];
	    }
	    i__1 = (c->cgcon_.iend - c->cgcon_.ist + 1) * 5;
        for (i = 0; i < i__1; i++) c->ftmp[i] = dum[c->cgcon_.ist * 5 - 5 + i];
	    MPI_Send(c->ftmp, i__1,	c->mpistuff_.dp_type__, c->neigh_.east, 4, MPI_COMM_WORLD);
	}
    } else {
	if (c->neigh_.north != -1) {
	    i__1 = c->cgcon_.jend;
	    for (j = c->cgcon_.jst; j <= i__1; ++j) {
		dum[j * 5 - 5] = g[((j + (*k << 4)) * 10 + 1) * 5 + 1];
		dum[j * 5 - 4] = g[((j + (*k << 4)) * 10 + 1) * 5 + 2];
		dum[j * 5 - 3] = g[((j + (*k << 4)) * 10 + 1) * 5 + 3];
		dum[j * 5 - 2] = g[((j + (*k << 4)) * 10 + 1) * 5 + 4];
		dum[j * 5 - 1] = g[((j + (*k << 4)) * 10 + 1) * 5 + 5];
	    }
	    i__1 = (c->cgcon_.jend - c->cgcon_.jst + 1) * 5;
        for (i = 0; i < i__1; i++) c->ftmp[i] = dum[c->cgcon_.jst * 5 - 5 + i];
	    MPI_Send(c->ftmp, i__1, c->mpistuff_.dp_type__, c->neigh_.north, 1, MPI_COMM_WORLD);
	}
	if (c->neigh_.west != -1) {
	    i__1 = c->cgcon_.iend;
	    for (i__ = c->cgcon_.ist; i__ <= i__1; ++i__) {
		dum[i__ * 5 - 5] = g[(i__ + ((*k << 4) + 1) * 10) * 5 + 1];
		dum[i__ * 5 - 4] = g[(i__ + ((*k << 4) + 1) * 10) * 5 + 2];
		dum[i__ * 5 - 3] = g[(i__ + ((*k << 4) + 1) * 10) * 5 + 3];
		dum[i__ * 5 - 2] = g[(i__ + ((*k << 4) + 1) * 10) * 5 + 4];
		dum[i__ * 5 - 1] = g[(i__ + ((*k << 4) + 1) * 10) * 5 + 5];
	    }
	    i__1 = (c->cgcon_.iend - c->cgcon_.ist + 1) * 5;
        for (i = 0; i < i__1; i++) c->ftmp[i] = dum[c->cgcon_.ist * 5 - 5 + i];
	    MPI_Send(c->ftmp, i__1, c->mpistuff_.dp_type__, c->neigh_.west, 3, MPI_COMM_WORLD);
	}
    }
    return 0;
} /* exchange_1__ */

/* Subroutine */ int bcast_inputs__(context *c)
{
    int *itmp;
    float *ftmp;

    itmp = malloc(sizeof(int));
    ftmp = malloc(sizeof(float));

    *itmp = c->cprcon_.ipr;
    MPI_Bcast(itmp, 1, MPI_INT, c->mpistuff_.root, MPI_COMM_WORLD);
    c->cprcon_.ipr = *itmp;
    *itmp = c->cprcon_.inorm;
    MPI_Bcast(itmp, 1, MPI_INT, c->mpistuff_.root, MPI_COMM_WORLD)	    ;
    c->cprcon_.inorm = *itmp;
    *itmp = c->ctscon_.itmax;
    MPI_Bcast(itmp, 1, MPI_INT, c->mpistuff_.root, MPI_COMM_WORLD)	    ;
    c->ctscon_.itmax = *itmp;
    *ftmp = c->ctscon_.dt;
    MPI_Bcast(ftmp, 1, c->mpistuff_.dp_type__, c->mpistuff_.root, 	    MPI_COMM_WORLD);
    c->ctscon_.dt = *ftmp;
    *ftmp = c->ctscon_.omega;
    MPI_Bcast(ftmp, 1, c->mpistuff_.dp_type__, 	    c->mpistuff_.root, MPI_COMM_WORLD);
    c->ctscon_.omega = *ftmp;
    MPI_Bcast(c->ctscon_.tolrsd, 5, c->mpistuff_.dp_type__, c->mpistuff_.root, MPI_COMM_WORLD);
    *itmp = c->cgcon_.nx0;
    MPI_Bcast(itmp, 1, MPI_INT, c->mpistuff_.root, MPI_COMM_WORLD);
    c->cgcon_.nx0 = *itmp;
    *itmp = c->cgcon_.ny0;
    MPI_Bcast(itmp, 1, MPI_INT, c->mpistuff_.root, MPI_COMM_WORLD);
    c->cgcon_.ny0 = *itmp;
    *itmp = c->cgcon_.nz0;
    MPI_Bcast(itmp, 1, MPI_INT, c->mpistuff_.root, MPI_COMM_WORLD);
    c->cgcon_.nz0 = *itmp;
    *itmp = c->timer_.timeron;
    MPI_Bcast(itmp, 1, MPI_INT, c->mpistuff_.root, MPI_COMM_WORLD);
    c->timer_.timeron = *itmp;

    free(itmp);
    free(ftmp);
    return 0;
} /* bcast_inputs__ */

/* Subroutine */ int read_input__(context *c)
{

    integer nnodes, fstatus;

    c->mpistuff_.root = 0;
    if (c->dim_.id == c->mpistuff_.root) {
	fstatus = -1;
	c->timer_.timeron = FALSE_;
	if (fstatus == 0) {
	} else {
	    c->cprcon_.ipr = 1;
	    c->cprcon_.inorm = INORM_DEFAULT;
	    c->ctscon_.itmax = ITMAX_DEFAULT;
	    c->ctscon_.dt = DT_DEFAULT;
	    c->ctscon_.omega = 1.2;
	    c->ctscon_.tolrsd[0] = 1e-8;
	    c->ctscon_.tolrsd[1] = 1e-8;
	    c->ctscon_.tolrsd[2] = 1e-8;
	    c->ctscon_.tolrsd[3] = 1e-8;
	    c->ctscon_.tolrsd[4] = 1e-8;
	    c->cgcon_.nx0 = ISIZ01;
	    c->cgcon_.ny0 = ISIZ02;
	    c->cgcon_.nz0 = ISIZ03;
	}
	MPI_Comm_size(MPI_COMM_WORLD, &nnodes);
	if (nnodes != NNODES_COMPILED) {
	}
	if (c->cgcon_.nx0 < 4 || c->cgcon_.ny0 < 4 || c->cgcon_.nz0 < 4) {
	    MPI_Abort(MPI_COMM_WORLD, 16);
	}
	if (c->cgcon_.nx0 > ISIZ01 || c->cgcon_.ny0 > ISIZ02 || c->cgcon_.nz0 > ISIZ03) {
	    MPI_Abort(MPI_COMM_WORLD, 16);
	}
        printf(" Size: %dx %dx %d\n", c->cgcon_.nx0, c->cgcon_.ny0, c->cgcon_.nz0);
        printf(" Iterations = %d\n", c->ctscon_.itmax);
        printf(" Number of processes: %d\n", nnodes);
    }
    bcast_inputs__(c);
    return 0;
} /* read_input__ */

#ifdef ARCH_MB
int nas_lu_mpi()
#else
int main(void)
#endif

{


    real r__1;

    logical verified;
    integer i__;
    float t1[12];
    float tsum[12];
    char class__[1];
    float tming[12], tmaxg[12];
    float mflops;
    context c;

    c.cvar_.u = malloc(5*(ISIZ1+4)*(ISIZ2+4)*ISIZ3*sizeof(float));
    c.cvar_.rsd = malloc(5*(ISIZ1+4)*(ISIZ2+4)*ISIZ3*sizeof(float));
    c.cvar_.frct = malloc(5*(ISIZ1+4)*(ISIZ2+4)*ISIZ3*sizeof(float));
    c.cvar_.flux = malloc(5*(ISIZ1+2)*(ISIZ2+2)*ISIZ3*sizeof(float));
    c.cjac_.a = malloc(5*5*ISIZ1*ISIZ2*sizeof(float));
    c.cjac_.b = malloc(5*5*ISIZ1*ISIZ2*sizeof(float));
    c.cjac_.c__ = malloc(5*5*ISIZ1*ISIZ2*sizeof(float));
    c.cjac_.d__ = malloc(5*5*ISIZ1*ISIZ2*sizeof(float));
    c.cexact_.ce = malloc(65*sizeof(float));
    c.comm_.buf = malloc(5*2*ISIZ2*ISIZ3*sizeof(float));
    c.comm_.buf1 = malloc(5*2*ISIZ2*ISIZ3*sizeof(float));
    c.ctscon_.tolrsd = malloc(5*sizeof(float));
    c.ctscon_.errnm = malloc(5*sizeof(float));
    c.ftmp = malloc(1024*sizeof(float));

    init_comm__(&c);
    read_input__(&c);
    for (i__ = 1; i__ <= 10; ++i__) {
	timer_clear(i__, &c);
    }
    proc_grid__(&c);
    neighbors_(&c);
    subdomain_(&c);
    setcoeff_(&c);
    setbv_(&c);
    setiv_(&c);
    erhs_(&c);
    integer c__1 = 1;
    ssor_(&c__1, &c);
    setbv_(&c);
    setiv_(&c);
    ssor_(&c.ctscon_.itmax, &c);
    error_(&c);
    pintgr_(&c);
    if (c.dim_.id == 0) {
	verify_(c.ctscon_.rsdnm, c.ctscon_.errnm, &c.ctscon_.frc, class__, &
		verified, (ftnlen)1, &c);
	r__1 = (real) (c.cgcon_.nx0 + c.cgcon_.ny0 + c.cgcon_.nz0) / 3.f;
	mflops = (real) c.ctscon_.itmax * ((real) c.cgcon_.nx0 * 1984.77f * (
		real) c.cgcon_.ny0 * (real) c.cgcon_.nz0 - r__1 * r__1 * 
		10923.3f + (real) (c.cgcon_.nx0 + c.cgcon_.ny0 + c.cgcon_.nz0) * 
		27770.9f / 3.f - 144010.f) / (c.timer_.maxtime * 1e6f);
        int c__2 = NNODES_COMPILED;
	print_results__("LU", class__, &c.cgcon_.nx0, &c.cgcon_.ny0, &
		c.cgcon_.nz0, &c.ctscon_.itmax, &c__2, &c.dim_.num, &
		c.timer_.maxtime, &mflops, "          floating point", &
		verified, "3.3.1", "25 May 2012", "f77 ", "f77", "-lmpi", 
		"(none)", "-O3", "(none)", "(none)", (ftnlen)2, (ftnlen)1, (
		ftnlen)24, (ftnlen)5, (ftnlen)11, (ftnlen)4, (ftnlen)3, (
		ftnlen)5, (ftnlen)6, (ftnlen)3, (ftnlen)6, (ftnlen)6);
    }
    if (! c.timer_.timeron) {
	goto L999;
    }
    for (i__ = 1; i__ <= 10; ++i__) {
	t1[i__ - 1] = timer_read(i__, &c);
    }
    t1[1] -= t1[6];
    t1[11] = t1[7] + t1[8] + t1[9] + t1[6];
    t1[10] = t1[0] - t1[11];
    MPI_Reduce(t1, &tsum, 12, c.mpistuff_.dp_type__, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(t1, &tming, 12, c.mpistuff_.dp_type__, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(t1, &tmaxg, 12, c.mpistuff_.dp_type__, MPI_MAX, 0, MPI_COMM_WORLD);
    if (c.dim_.id == 0) {
	for (i__ = 1; i__ <= 12; ++i__) {
	    tsum[i__ - 1] /= c.dim_.num;
	}
    }
L999:
while(1);
    MPI_Finalize();
    return 0;
} /* MAIN__ */

