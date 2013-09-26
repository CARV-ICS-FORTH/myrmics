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
// This code is derived by the one by Michael Andersch, which was under
// the following copyright notice:
//
// Copyright (C) 2013 Michael Andersch <michael.andersch@mailbox.tu-berlin.de>
//
// This file is part of Starbench.
//
// Starbench is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Starbench is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Starbench.  If not, see <http://www.gnu.org/licenses/>.
//
// ==========================[ Static Information ]===========================
//
// Author        : Spyros Lyberis
// Abstract      : k-means clustering kernel, MPI version
//
// =============================[ CVS Variables ]=============================
// 
// File name     : $RCSfile: kmeans_mpi.c,v $
// CVS revision  : $Revision: 1.1 $
// Last modified : $Date: 2013/01/22 17:55:49 $
// Last author   : $Author: lyberis-spree $
// 
// ===========================================================================

#include <arch.h>
#include <kernel_toolset.h>
#include <fmpi.h>


#define COORDS  3   /* no. coordinates */



// ===========================================================================
// ===========================================================================
static inline float euclid_dist_2(float *coord1, float *coord2) {

  float ans;
  int   i;

  ans = 0.0F;

  for (i = 0; i < COORDS; i++)
    ans += (coord1[i] - coord2[i]) * (coord1[i] - coord2[i]);

  return (ans);
}


// ===========================================================================
// ===========================================================================
static inline int find_nearest_cluster(int num_clusters, float *object, 
                                       float *clusters) {
  int   index;
  int   i;
  float dist;
  float min_dist;

  // Find the cluster id that has min distance to object
  index = 0;
  min_dist = euclid_dist_2(object, clusters);

  for (i = 1; i < num_clusters; i++) {
    dist = euclid_dist_2(object, clusters + i * COORDS);

    // No need square root
    if (dist < min_dist) { // Find the min and its array index
      min_dist = dist;
      index = i;
    }
  }

  return index;
}


// ===========================================================================
// ===========================================================================
void kmeans_mpi_do_tile(float *objects, int *membership, 
                        float *partial_clusters, int *partial_sizes, 
                        float *clusters, int objects_per_tile, 
                        int num_clusters) {
  int   i;
  int   j;

//printf("%d: do_tile 0x%X 0x%X 0x%X 0x%X 0x%X\r\n", sys_get_worker_id(), objects, membership, partial_clusters, partial_sizes, clusters);

  // Clear partial arrays
  for (i = 0; i < num_clusters; i++) {
    partial_sizes[i] = 0;
    for (j = 0; j < COORDS; j++) {
      partial_clusters[i * COORDS + j] = 0.0F;
    }
  }

  // Cluster
  for (i = 0; i < objects_per_tile; i++) {
    membership[i] = find_nearest_cluster(num_clusters, objects + i * COORDS, 
                                         clusters);
    partial_sizes[membership[i]]++;
    for (j = 0; j < COORDS; j++) {
      partial_clusters[membership[i] * COORDS + j] += objects[i * COORDS + j];
    }
  }

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
int kmeans_mpi(int num_procs,          // MPI processors to use
               int num_clusters,       // number of output clusters
               int num_objects,        // number of total input objects
               int num_reps) {         // loop repetitions

  int           num_cores;
  int           rank;
  int           objects_per_core;
  float         *objects;
  int           *membership;
  float         *clusters;
  float         *partial_clusters;
  int           *partial_sizes;
  float         *reduce_clusters;
  int           *reduce_sizes;
  unsigned int  seed = 42;
  unsigned int  time_start = 0;
  unsigned int  time_stop;
  unsigned int  time;
  int           i;
  int           j;
  int           loop;


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
  if (num_objects % num_procs) {
    kt_printf("%d objects not divisible by %d cores\r\n", 
              num_objects, num_procs);
    return 1;
  }
  objects_per_core = num_objects / num_procs;


  // Synchronize everyone and print infomercial
  MPI_Barrier(MPI_COMM_WORLD);
  if (!rank) {
    kt_printf("k-means of %d -> %d starting on %d core(s)\r\n",
              num_objects, num_clusters, num_procs);
  }
  MPI_Barrier(MPI_COMM_WORLD);




  // Create random input objects per core
  objects = kt_malloc(objects_per_core * COORDS * sizeof(float));
  for (i = 0; i < objects_per_core; i++) {
    for (j = 0; j < COORDS; j++) {
      objects[i * COORDS + j] = (float) ((seed = kt_rand(seed)) % 1000) / 10.0F;
    }
  }

  // Create membership, no object belongs to any cluster yet
  membership = kt_malloc(objects_per_core * sizeof(int));
  for (i = 0; i < objects_per_core; i++) {
    membership[i] = -1;
  }

  // Create partial cluster arrays, to be used as temporary buffers, and their
  // size arrays. Init everything to a neutral value, so they can be used
  // for reductions even with cores that are not part of the setup.
  partial_clusters = kt_malloc(num_clusters * COORDS * sizeof(float));
  partial_sizes = kt_malloc(num_clusters * sizeof(int));
  for (i = 0; i < num_clusters; i++) {
    partial_sizes[i] = 0;
    for (j = 0; j < COORDS; j++) {
      partial_clusters[i * COORDS + j] = 0.0F;
    }
  }

  // Rank 0 allocates extra buffers for the reductions
  if (!rank) {
    reduce_clusters = kt_malloc(num_clusters * COORDS * sizeof(float));
    reduce_sizes = kt_malloc(num_clusters * sizeof(int));
  }
  else {
    reduce_clusters = NULL;
    reduce_sizes = NULL;
  }


  // Allocate stable clusters (will be changed across loop repetitions) for
  // everyone. Prepare initial clusters for core 0 with centers copied from
  // first objects.
  clusters = kt_malloc(num_clusters * COORDS * sizeof(float *));

  if (!rank) {
    ar_assert(objects_per_core >= num_clusters);

    for (i = 0; i < num_clusters; i++) {
      for (j = 0; j < COORDS; j++) {
        clusters[i * COORDS + j] = objects[i * COORDS + j];
      }
    }
  }


  // Keep time
  MPI_Barrier(MPI_COMM_WORLD);
  if (!rank) {
    time_start = ar_free_timer_get_ticks();
  }


  // For all repetitions
  for (loop = 0; loop < num_reps; loop++) {

    // Broadcast clusters for this rep from core 0 to everyone
    MPI_Bcast(clusters, num_clusters * COORDS, MPI_FLOAT, 0, MPI_COMM_WORLD);

    // Process our objects (only for cores belonging in the setup).
    if (rank < num_procs) {
      kmeans_mpi_do_tile(objects, membership, partial_clusters, partial_sizes, 
                         clusters, objects_per_core, num_clusters);
    }

    // Reduce all results to rank 0 (non-setup cores will add their 0 values)
    MPI_Reduce(partial_clusters, reduce_clusters, num_clusters * COORDS,
               MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(partial_sizes, reduce_sizes, num_clusters,
               MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    // Rank 0 prepares the next clusters, based on reduction values
    if (!rank) {
      kt_printf("done %d of %d\r\n", loop + 1, num_reps);
      for (i = 0; i < num_clusters; i++) {
        for (j = 0; j < COORDS; j++) {
          clusters[i * COORDS + j] = reduce_clusters[i * COORDS + j] /
                                     (float) reduce_sizes[i];
        }
      }
    }
  }


  MPI_Barrier(MPI_COMM_WORLD);

  // Compute elapsed time
  if (!rank) {
    time_stop = ar_free_timer_get_ticks();
    if (time_stop > time_start) {
      time = time_stop - time_start;
    }
    else {
      time = 0xFFFFFFFF - (time_start - time_stop);
    }
    kt_printf("Time: %10u cycles (%6u msec)\r\n", time, time / 10000);

    // Print clusters
    for (i = 0; i < num_clusters; i++) {
      kt_printf("Cluster %d: %f %f %f\r\n", 
                i, 
                clusters[i * COORDS], 
                clusters[i * COORDS + 1], 
                clusters[i * COORDS + 2]);
    }
  }

  // Free stuff
  kt_free(clusters);
  kt_free(partial_sizes);
  kt_free(partial_clusters);
  kt_free(reduce_sizes);
  kt_free(reduce_clusters);
  kt_free(objects);
  kt_free(membership);

  return 0;
}

