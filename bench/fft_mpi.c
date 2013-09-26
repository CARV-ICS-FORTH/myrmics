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
// ==========================[ Static Information ]===========================
//
// Author        : Spyros Lyberis
// Abstract      : 2D-FFT kernel, MPI version. Adapted from "Introduction to
//                 Parallel Computing", Petersen & Arbenz.
//
// =============================[ CVS Variables ]=============================
// 
// File name     : $RCSfile: fft_mpi.c,v $
// CVS revision  : $Revision: 1.2 $
// Last modified : $Date: 2012/11/05 12:49:54 $
// Last author   : $Author: lyberis-spree $
// 
// ===========================================================================

#include <arch.h>
#include <kernel_toolset.h>
#include <fmpi.h>


// ===========================================================================
// ===========================================================================
void fft_mpi_step (int n, int mj, float *a, float *b, float *c, 
               float *d, float *w, float sign) {

  int j, k, jc, jw, lj, mj2;
  float rp, up;
  float wr[4], wu[4];
  float V0[4], V1[4], V2[4], V3[4], V4[4], V6[4], V7[4];

  mj2 = 2 * mj;
  lj = n / mj2;
  for (j = 0; j < lj; j++) {
    jw = j * mj;
    jc = j * mj2;
    rp = w[2 * jw + 0];
    up = w[2 * jw + 1];
    if (sign < 0.0)
      up = -up;
    if (mj < 2) {               /* special case mj=1 */
      d[2*jc+0] = rp * (a[2*jw+0] - b[2*jw+0]) - up * (a[2*jw+1] - b[2*jw+1]);
      d[2*jc+1] = up * (a[2*jw+0] - b[2*jw+0]) + rp * (a[2*jw+1] - b[2*jw+1]);
      c[2*jc+0] = a[2*jw+0] + b[2*jw+0];
      c[2*jc+1] = a[2*jw+1] + b[2*jw+1];
    }
    else {                      /* mj > 1 cases */
      wr[0] = rp;
      wr[1] = rp;
      wr[2] = rp;
      wr[3] = rp;

      wu[0] = -up;
      wu[1] = up;
      wu[2] = -up;
      wu[3] = up;

      //V6 = _mm_load_ps (wr);
      V6[0] = wr[0];
      V6[1] = wr[1];
      V6[2] = wr[2];
      V6[3] = wr[3];

      //V7 = _mm_load_ps (wu);
      V7[0] = wu[0];
      V7[1] = wu[1];
      V7[2] = wu[2];
      V7[3] = wu[3];

      for (k = 0; k < mj; k += 2) {
        //V0 = _mm_load_ps (&a[jw + k][0]);
        V0[0] = a[2 * jw + k + 0];
        V0[1] = a[2 * jw + k + 1];
        V0[2] = a[2 * jw + k + 2];
        V0[3] = a[2 * jw + k + 3];

        //V1 = _mm_load_ps (&b[jw + k][0]);
        V1[0] =b[2 * jw + k + 0];
        V1[1] =b[2 * jw + k + 1];
        V1[2] =b[2 * jw + k + 2];
        V1[3] =b[2 * jw + k + 3];

        //V2 = _mm_add_ps (V0, V1); /* a+b */
        V2[0] = V0[0] + V1[0];
        V2[1] = V0[1] + V1[1];
        V2[2] = V0[2] + V1[2];
        V2[3] = V0[3] + V1[3];

        //_mm_store_ps (&c[jc + k][0], V2); /* c to M */
        c[2 * jc + k + 0] = V2[0];
        c[2 * jc + k + 1] = V2[1];
        c[2 * jc + k + 2] = V2[2];
        c[2 * jc + k + 3] = V2[3];

        //V3 = _mm_sub_ps (V0, V1); /* a-b */
        V3[0] = V0[0] - V1[0];
        V3[1] = V0[1] - V1[1];
        V3[2] = V0[2] - V1[2];
        V3[3] = V0[3] - V1[3];

        //V4 = _mm_shuffle_ps (V3, V3, _MM_SHUFFLE (2, 3, 0, 1));
        V4[0] = V3[1];
        V4[1] = V3[0];
        V4[2] = V3[3];
        V4[3] = V3[2];

        //V0 = _mm_mul_ps (V6, V3);
        V0[0] = V6[0] * V3[0];
        V0[1] = V6[1] * V3[1];
        V0[2] = V6[2] * V3[2];
        V0[3] = V6[3] * V3[3];

        //V1 = _mm_mul_ps (V7, V4);
        V1[0] = V7[0] * V4[0];
        V1[1] = V7[1] * V4[1];
        V1[2] = V7[2] * V4[2];
        V1[3] = V7[3] * V4[3];

        //V2 = _mm_add_ps (V0, V1); /* w*(a-b) */
        V2[0] = V0[0] + V1[0];
        V2[1] = V0[1] + V1[1];
        V2[2] = V0[2] + V1[2];
        V2[3] = V0[3] + V1[3];

        //_mm_store_ps (&d[jc + k][0], V2); /* d to M */
        d[2 * jc + k + 0] = V2[0];
        d[2 * jc + k + 1] = V2[1];
        d[2 * jc + k + 2] = V2[2];
        d[2 * jc + k + 3] = V2[3];
      }
    }
  }
}

// ===========================================================================
/* x=in, y=out, w=exp(2*pi*i*k/n), k=0..n/2-1 */
// ===========================================================================
void fft_mpi_cfft2 (int n, float *x, float *y, float *w, float sign) {

  int m, j, mj, tgle, i;

  m = kt_int_log2(n); //(int) (log ((float) n) / log (1.99));
  mj = 1;
  tgle = 1;

  fft_mpi_step (n, mj, x, x + n, y, y + 2 * mj, w, sign);

  for (j = 0; j < m - 2; j++) {
    mj *= 2;
    if (tgle) {
      fft_mpi_step (n, mj, y, y + n, x, x + 2 * mj, w, sign);
      tgle = 0;
    }
    else {
      fft_mpi_step (n, mj, x, x + n, y, y + 2 * mj, w, sign);
      tgle = 1;
    }
  }
  if (tgle) {
    for (i = 0; i < n; i++) {
      y[i] = x[i];
    }
  }

  mj = n / 2;

  fft_mpi_step (n, mj, x, x + n, y, y + 2 * mj, w, sign);
}

// ===========================================================================
// ===========================================================================
void fft_mpi_Xpose (float *a, int n, int num_procs, int rank) {

  float t0, t1;
  float *buf_i, *buf_o;
  int i, ij, is, j, step, n2, nn, other;
  MPI_Status status;
  MPI_Request reqs[2];

  nn = n / num_procs;
  n2 = 2 * nn;

  /* allocate buffers */
  buf_i = (float *) kt_malloc (nn * n2 * sizeof (float));
  buf_o = (float *) kt_malloc (nn * n2 * sizeof (float));
  
  /* local transpose of first block (in-place) */
  for (j = 0; j < nn; j++) {
    for (i = 0; i < j; i++) {
      t0 = a[rank * n2 + i * 2 * n + j * 2];
      t1 = a[rank * n2 + i * 2 * n + j * 2 + 1];
      a[rank * n2 + i * 2 * n + j * 2] = a[rank * n2 + j * 2 * n + 2 * i];
      a[rank * n2 + i * 2 * n + j * 2 + 1] =
        a[rank * n2 + j * 2 * n + 2 * i + 1];
      a[rank * n2 + j * 2 * n + 2 * i] = t0;
      a[rank * n2 + j * 2 * n + 2 * i + 1] = t1;
    }
  }
  
  /* num_procs-1 communication steps */
  for (step = 1; step < num_procs; step++) {
    other = rank ^ step;    /* XOR trick */
    ij = 0;
    for (i = 0; i < nn; i++) { /* fill send buffer */
      is = other * n2 + i * 2 * n;
      for (j = 0; j < n2; j++) {
        buf_o[ij++] = a[is + j];
      }
    }
    
    /* exchange data */
    //MPI_Sendrecv_replace (buf_io, 2 * nn * nn, MPI_FLOAT, other, rank, 
    //                      other, other, MPI_COMM_WORLD, &stat); 
    MPI_Irecv(buf_i, 2 * nn * nn, MPI_FLOAT, other, other,
              MPI_COMM_WORLD, &(reqs[0]));
    MPI_Isend(buf_o, 2 * nn * nn, MPI_FLOAT, other, rank,
              MPI_COMM_WORLD, &(reqs[1]));
//kt_printf("%d: waiting for %d...\r\n", rank, other);
    MPI_Waitall(2, reqs, &status);
//kt_printf("%d: wait for %d ok\r\n", rank);
    
    /* write back recv buffer in transposed order */
    for (i = 0; i < nn; i++) {
      for (j = 0; j < nn; j++) {
        a[other * n2 + j * 2 * n + i * 2] = buf_i[i * n2 + j * 2];
        a[other * n2 + j * 2 * n + i * 2 + 1] =
          buf_i[i * n2 + j * 2 + 1];
      }
    }
  }

  /* free buffers */
  kt_free(buf_i);
  kt_free(buf_o);
}

// ===========================================================================
// n:  num columns
// ny: num my_lines
// ===========================================================================
void fft_FFT2D (float *a, float *w, float sign, int ny, int n, int num_procs, 
                int rank) {
  int i, off;
  float *pa;


  for (i = 0; i < ny; i++) {
    off = 2 * i * n;
    pa = a + off;
    fft_mpi_cfft2 (n, pa, pa, w, sign);

//kt_printf("%d: fft1 ok\r\n", rank);

  }

  fft_mpi_Xpose (a, n, num_procs, rank);
//kt_printf("%d: xpose1 ok\r\n", rank);

  for (i = 0; i < ny; i++) {
    off = 2 * i * n;
    pa = a + off;
    fft_mpi_cfft2 (n, pa, pa, w, sign);
//kt_printf("%d: fft2 ok\r\n", rank);
  }

  fft_mpi_Xpose (a, n, num_procs, rank);
//kt_printf("%d: xpose2 ok\r\n", rank);
}


// ===========================================================================
// function()                   FIXME comments
// ===========================================================================
// * INPUTS
//   unsigned char *arg1        Describe arg1
//   int arg2                   Describe arg2
//
// * OUTPUTS
//   int *arg3                  Describe arg3
//
// * RETURN VALUE
//   int                        0 for success
// ===========================================================================
int fft_mpi(int num_procs, int n) {

  int           num_cores;
  int           rank;
  float         *a = NULL;
  float         *w = NULL;
#ifdef ARCH_MB
  unsigned int  time_start = 0;
  unsigned int  time_stop;
  unsigned int  time;
#endif
  int           i;


  // Who are we?
  MPI_Comm_size(MPI_COMM_WORLD, &num_cores);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // Sanity checks
  if (num_cores < num_procs) {
    if (!rank) {
      kt_printf("Cannot run with %d cores, MPI setup has only %d cores\r\n",
                num_procs, num_cores);
    }
    return 1;
  }
  if (n & (n - 1)) {
    if (!rank) {
      kt_printf("N must be a power of 2\r\n");
    }
    return 1;
  }
  if (num_procs & (num_procs - 1)) {
    if (!rank) {
      kt_printf("Cores must be a power of 2\r\n");
    }
    return 1;
  }


  // Allocate buffers for our portion of the table
  if (rank < num_procs) {
    a = kt_malloc(2 * (n / num_procs) * n * sizeof(float));
    w = kt_malloc(2 * n * sizeof(float));

    // Initialize them with some values
    for (i = 0; i < 2 * n / num_procs * n; i++) {
      a[i] = i * 1.2;
    }
    for (i = 0; i < 2 * n; i++) {
      w[i] = i * 0.3;
    }
  }


  // Synchronize everyone and print infomercial
  MPI_Barrier(MPI_COMM_WORLD);
  if (!rank) {
    kt_printf("FFT 2D-block of %d x %d starting on %d core(s)\r\n",
              n, n, num_procs);
  }
  MPI_Barrier(MPI_COMM_WORLD);

  // Keep time
  if (!rank) {
#ifdef ARCH_MB
    time_start = ar_glob_timer_read();
#endif
  }

  // This kernel is reentrant for multiple MPI runs. If our core is not
  // part of the current num_procs setup, go to the next barrier directly.
  if (rank >= num_procs) {
    goto skip;
  }


  // Run FFT
  fft_FFT2D(a, w, 1.0, n / num_procs, n, num_procs, rank);


skip:

  MPI_Barrier(MPI_COMM_WORLD);

  // Keep time
  if (!rank) {
#ifdef ARCH_MB
    time_stop = ar_glob_timer_read();
    if (time_stop > time_start) {
      time = time_stop - time_start;
    }
    else {
      time = 0xFFFFFFFF - (time_start - time_stop);
    }
    kt_printf("Time: %10u cycles (%6u msec)\r\n", time, time / 10000);
#endif
  }

  // Free stuff
  if (rank < num_procs) {
    kt_free(a);
    kt_free(w);
  }

  return 0;
}
