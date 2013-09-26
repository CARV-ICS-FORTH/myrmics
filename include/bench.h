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
// Abstract      : Benchmarks
//
// =============================[ CVS Variables ]=============================
// 
// File name     : $RCSfile: bench.h,v $
// CVS revision  : $Revision: 1.20 $
// Last modified : $Date: 2013/03/27 09:55:03 $
// Last author   : $Author: jacob $
// 
// ===========================================================================

#ifndef BENCH_H
#define BENCH_H

// MPI test routines
extern int test_mpi();


// MPI small kernels
extern int matrix_mult_mpi(int num_procs, int tile_size, int blk_size,
                           int verify);

extern int jacobi_mpi(int num_procs, int num_rows, int num_cols, int num_iter);

extern int swat_mpi(int num_procs, int seq0_len, int seq1_len, int seq1_step,
                    int min_rows_per_stripe);

extern int bitonic_mpi(int num_procs, int num_elements, int verify);

extern int fft_mpi(int num_procs, int n);

extern int cray_mpi(int num_procs, int xres, int yres);

extern int kmeans_mpi(int num_procs, int num_clusters, int num_objects, 
                      int num_reps);

// Sandia MPI benchmark suite
extern int smb_overhead_mpi(int direction, int iterations, int data_size,
                            float threshold, float base_threshold, int nohdr,
                            int verbose);

extern int smb_msgrate_mpi(int npeers, int niters, int nmsgs, int nbytes,
                           int cache_size, int ppn, int machine_output);


// Bare-metal benchmarks
extern int jacobi_bare(int num_procs, int num_rows, int num_cols, int num_iter);

extern int stream_bare(int n);

extern int test_bare();


// NAS benchmarks
extern int nas_is_mpi();
extern int nas_cg_mpi();
extern int nas_dt_mpi();
extern int nas_ep_mpi();
extern int nas_lu_mpi();

// 3D Barnes-Hut application
extern int barnes_mpi(int num_procs, int num_particles, int reps);

// Myrmics benchmarks
extern void (**test_myrmics_task_table)();
extern void (**matrix_mult_myrmics_task_table)();
extern void (**jacobi_myrmics_task_table)();
extern void (**bitonic_myrmics_task_table)();
extern void (**fft_myrmics_task_table)();
extern void (**kmeans_myrmics_task_table)();
extern void (**cray_myrmics_task_table)();
extern void (**md5_myrmics_task_table)();
extern void (**multilevel_myrmics_task_table)();

// 3D Barnes-Hut application
extern void (**barnes_myrmics_task_table)();

#endif
