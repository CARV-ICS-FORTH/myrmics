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
// Abstract      : k-means clustering kernel, Myrmics version
//
// =============================[ CVS Variables ]=============================
// 
// File name     : $RCSfile: kmeans_myrmics.c,v $
// CVS revision  : $Revision: 1.7 $
// Last modified : $Date: 2013/01/14 17:17:09 $
// Last author   : $Author: lyberis-spree $
// 
// ===========================================================================

#include <myrmics.h>


#define COORDS  3   /* no. coordinates */

// object region structure
typedef struct {

  float         **objects;              // [tile][object x coord]
  int           **membership;           // [tile][object]

} per_obj_region_t;


// reduction region structure
typedef struct {

  float         *clusters;              // [cluster x coord]

  float         **partial_clusters;     // [tile][cluster x coord]
  int           **partial_sizes;        // [tile][cluster]

} per_red_region_t;


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
void kmeans_myrmics_do_tile(float *objects, int *membership, 
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
// ===========================================================================
void kmeans_myrmics_reduce_tiles(rid_t red_region, per_red_region_t *per_red,
                                 int tiles_per_region, int num_clusters) {

  int   i;
  int   j;
  int   k;

//printf("%d: reduce_tiles\r\n", sys_get_worker_id());

  // Merge tiles 1...N results to tile 0
  for (i = 1; i < tiles_per_region; i++) {
    for (j = 0; j < num_clusters; j++) {
      per_red->partial_sizes[0][j] += per_red->partial_sizes[i][j];
      for (k = 0; k < COORDS; k++) {
        per_red->partial_clusters[0][j * COORDS + k] += 
        per_red->partial_clusters[i][j * COORDS + k];
      }
    }
  }

}


// ===========================================================================
// ===========================================================================
void kmeans_myrmics_do_region(rid_t obj_region, per_obj_region_t *per_obj,
                              rid_t red_region, per_red_region_t *per_red,
                              int tiles_per_region, int objects_per_tile, 
                              int num_clusters) {
  int i;

//printf("%d: do_region\r\n", sys_get_worker_id());

  // Compute the tiles
  for (i = 0; i < tiles_per_region; i++) {

    float *scoop_objects          = per_obj->objects[i];
    int   *scoop_membership       = per_obj->membership[i];
    float *scoop_partial_clusters = per_red->partial_clusters[i];
    int   *scoop_partial_sizes    = per_red->partial_sizes[i];
    float *scoop_clusters         = per_red->clusters;

    #pragma myrmics task in(scoop_objects) inout(scoop_membership) \
                         inout(scoop_partial_clusters) \
                         inout(scoop_partial_sizes) in(scoop_clusters) \
                         in(objects_per_tile, num_clusters) \
                         safe(objects_per_tile, num_clusters)
    kmeans_myrmics_do_tile(scoop_objects, scoop_membership, 
                           scoop_partial_clusters, scoop_partial_sizes, 
                           scoop_clusters, objects_per_tile, num_clusters);
  }

  // Reduce their results
  #pragma myrmics task region inout(red_region) \
                       in(per_red) safe(per_red) \
                       in(tiles_per_region, num_clusters) \
                       safe(tiles_per_region, num_clusters)
  kmeans_myrmics_reduce_tiles(red_region, per_red, 
                              tiles_per_region, num_clusters);
}


// ===========================================================================
// ===========================================================================
void kmeans_myrmics_reduce(rid_t top_red_region, 
                           per_red_region_t **per_red_region,
                           int num_regions, int num_clusters) {
  int i;
  int j;
  int k;

printf("%d: reduce\r\n", sys_get_worker_id());

  // Merge all reduction results 1...N to 0
  for (i = 1; i < num_regions; i++) {
    for (j = 0; j < num_clusters; j++) {
      per_red_region[0]->partial_sizes[0][j] += 
      per_red_region[i]->partial_sizes[0][j];
      for (k = 0; k < COORDS; k++) {
        per_red_region[0]->partial_clusters[0][j * COORDS + k] += 
        per_red_region[i]->partial_clusters[0][j * COORDS + k];
      }
    }
  }

  // Replace clusters of region 0 with the average
  for (i = 0; i < num_clusters; i++) {
    for (j = 0; j < COORDS; j++) {
      per_red_region[0]->clusters[i * COORDS + j] = 
                per_red_region[0]->partial_clusters[0][i * COORDS + j] /
                (float) per_red_region[0]->partial_sizes[0][i];
    }
//printf("Cluster %d: %f %f %f\r\n", i, per_region[0]->clusters[i * COORDS], per_region[0]->clusters[i * COORDS + 1], per_region[0]->clusters[i * COORDS + 2]);
  }
//printf("\r\n");

  // Copy region 1...N clusters from region 0
  for (i = 1; i < num_regions; i++) {
    for (j = 0; j < num_clusters; j++) {
      for (k = 0; k < COORDS; k++) {
        per_red_region[i]->clusters[j * COORDS + k] =
        per_red_region[0]->clusters[j * COORDS + k];
      }
    }
  }

}


// ===========================================================================
// ===========================================================================
void kmeans_myrmics_print(rid_t top_red_region, 
                          per_red_region_t **per_red_region,
                          unsigned int time_start, int num_clusters) {
  
  unsigned int time_stop;
  unsigned int time;
  int i;
  int j;
  int k;

  // Compute elapsed time
  time_stop = sys_free_timer_get_ticks();
  if (time_stop > time_start) {
    time = time_stop - time_start;
  }
  else {
    time = 0xFFFFFFFF - (time_start - time_stop);
  }
  printf("Time: %10u cycles (%6u msec)\r\n", time, time / 10000);


  // Print clusters
  for (i = 0; i < num_clusters; i++) {
    printf("Cluster %d: %f %f %f\r\n", 
           i, 
           per_red_region[0]->clusters[i * COORDS], 
           per_red_region[0]->clusters[i * COORDS + 1], 
           per_red_region[0]->clusters[i * COORDS + 2]);
  }
}


// ===========================================================================
// total objects: num_regions x tiles_per_region x objects_per_tile
// ===========================================================================
void kmeans_myrmics(int num_clusters,       // number of output clusters
                    int num_regions,        // number of regions to use
                    int tiles_per_region,   // tiles per region
                    int objects_per_tile,   // input objects per tile
                    int num_reps) {         // loop repetitions

  rid_t                 top_obj_region;
  rid_t                 top_red_region;
  rid_t                 *obj_regions;       
  rid_t                 *red_regions;       
  per_obj_region_t      **per_obj_region;           // [region]
  per_obj_region_t      **per_obj_region_copy;
  per_red_region_t      **per_red_region;           // [region]
  per_red_region_t      **per_red_region_copy;
  unsigned int          seed = 42;
  unsigned int          time_start;
  int                   i;
  int                   j;
  int                   k;
  int                   l;
  int                   loop;


  // Infomercial
  printf(
    "k-means of %d -> %d starting in %d total tile(s) and %d region(s)\r\n",
    objects_per_tile * tiles_per_region * num_regions, 
    num_clusters,
    tiles_per_region * num_regions,
    num_regions);


  // Create top-level regions
  top_obj_region = sys_ralloc(0, 99); // highest level
  top_red_region = sys_ralloc(0, 99); // highest level

  // Create mid-level regions. Warning: to get decent performance,
  // obj_regions[x] and red_regions[x] must be located in the same scheduler
  // (i.e. make num_regions a multiple of the number of the low-level
  // schedulers)
  obj_regions = sys_alloc(num_regions * sizeof(rid_t), top_obj_region);
  red_regions = sys_alloc(num_regions * sizeof(rid_t), top_red_region);
  for (i = 0; i < num_regions; i++) {
    obj_regions[i] = sys_ralloc(top_obj_region, 0); // lowest level
  }
  for (i = 0; i < num_regions; i++) {
    red_regions[i] = sys_ralloc(top_red_region, 0); // lowest level
  }

  // Create the per_region structures and a copy for master task spawning
  per_obj_region = sys_alloc(num_regions * sizeof(per_obj_region_t *), 
                             top_obj_region);
  per_obj_region_copy = sys_alloc(num_regions * sizeof(per_obj_region_t *), 0);
  for (i = 0; i < num_regions; i++) {
    per_obj_region[i] = sys_alloc(sizeof(per_obj_region_t), obj_regions[i]);
    per_obj_region_copy[i] = per_obj_region[i];
  }

  per_red_region = sys_alloc(num_regions * sizeof(per_red_region_t *), 
                             top_red_region);
  per_red_region_copy = sys_alloc(num_regions * sizeof(per_red_region_t *), 0);
  for (i = 0; i < num_regions; i++) {
    per_red_region[i] = sys_alloc(sizeof(per_red_region_t), red_regions[i]);
    per_red_region_copy[i] = per_red_region[i];
  }


  // Fill the per_region structure
  for (i = 0; i < num_regions; i++) {

    // Create input random objects
    per_obj_region[i]->objects = sys_alloc(tiles_per_region * sizeof(float *), 
                                           obj_regions[i]);
    if (tiles_per_region > 1) {
      sys_balloc(objects_per_tile * COORDS * sizeof(float), obj_regions[i],
                 tiles_per_region, per_obj_region[i]->objects);
    }
    else {
      per_obj_region[i]->objects[0] = sys_alloc(
                objects_per_tile * COORDS * sizeof(float), obj_regions[i]);
    }
    for (j = 0; j < tiles_per_region; j++) {
      for (k = 0; k < objects_per_tile; k++) {
        for (l = 0; l < COORDS; l++) {
          per_obj_region[i]->objects[j][k * COORDS + l] = 
                                (float) ((seed = rand(seed)) % 1000) / 10.0F;
        }
      }
    }

    // Create membership, no object belongs to any cluster yet
    per_obj_region[i]->membership = sys_alloc(tiles_per_region * sizeof(int *), 
                                              obj_regions[i]);
    if (tiles_per_region > 1) {
      sys_balloc(objects_per_tile * sizeof(int), obj_regions[i],
                 tiles_per_region, per_obj_region[i]->membership);
    }
    else {
      per_obj_region[i]->membership[0] = sys_alloc(
                          objects_per_tile * sizeof(int), obj_regions[i]);
    }
    for (j = 0; j < tiles_per_region; j++) {
      for (k = 0; k < objects_per_tile; k++) {
        per_obj_region[i]->membership[j][k] = -1;
      }
    }

    // Create partial cluster arrays, to be used by tasks as temporary buffers,
    // and their size arrays. We leave these uninitialized.
    per_red_region[i]->partial_clusters = sys_alloc(
                          tiles_per_region * sizeof(float *), red_regions[i]);
    if (tiles_per_region > 1) {
      sys_balloc(num_clusters * COORDS * sizeof(float), red_regions[i],
                 tiles_per_region, per_red_region[i]->partial_clusters);
    }
    else {
      per_red_region[i]->partial_clusters[0] = sys_alloc(
                        num_clusters * COORDS * sizeof(float), red_regions[i]);
    }

    per_red_region[i]->partial_sizes = sys_alloc(
                          tiles_per_region * sizeof(int *), red_regions[i]);
    if (tiles_per_region > 1) {
      sys_balloc(num_clusters * sizeof(int), red_regions[i],
                 tiles_per_region, per_red_region[i]->partial_sizes);
    }
    else {
      per_red_region[i]->partial_sizes[0] = sys_alloc(
                              num_clusters * sizeof(int), red_regions[i]);
    }

    printf("init %d/%d done\r\n", i+1, num_regions);
  }


  // Create clusters for region 0, initially with centers copied from first
  // objects
  sys_assert(objects_per_tile >= num_clusters);
  per_red_region[0]->clusters = sys_alloc(
                      num_clusters * COORDS * sizeof(float *), red_regions[0]);
  for (i = 0; i < num_clusters; i++) {
    for (j = 0; j < COORDS; j++) {
      per_red_region[0]->clusters[i * COORDS + j] = 
        per_obj_region[0]->objects[0][i * COORDS + j];
    }
  }

  // Allocate and copy clusters to all other regions
  for (i = 1; i < num_regions; i++) {
    per_red_region[i]->clusters = sys_alloc(
                      num_clusters * COORDS * sizeof(float *), red_regions[i]);
    for (j = 0; j < num_clusters; j++) {
      for (k = 0; k < COORDS; k++) {
        per_red_region[i]->clusters[j * COORDS + k] = 
          per_red_region[0]->clusters[j * COORDS + k];
      }
    }
  }


  // Start time
  time_start = sys_free_timer_get_ticks();


  // For all repetitions
  for (loop = 0; loop < num_reps; loop++) {

    for (i = 0; i < num_regions; i++) {

      rid_t             scoop_obj_region      = obj_regions[i];
      per_obj_region_t  *scoop_per_obj_region = per_obj_region_copy[i];
      rid_t             scoop_red_region      = red_regions[i];
      per_red_region_t  *scoop_per_red_region = per_red_region_copy[i];

      #pragma myrmics task \
                      region inout(scoop_obj_region) \
                      in(scoop_per_obj_region) safe(scoop_per_obj_region) \
                      region inout(scoop_red_region) \
                      in(scoop_per_red_region) safe(scoop_per_red_region) \
                      in(tiles_per_region, objects_per_tile, num_clusters) \
                      safe(tiles_per_region, objects_per_tile, num_clusters)
      kmeans_myrmics_do_region(scoop_obj_region, scoop_per_obj_region, 
                               scoop_red_region, scoop_per_red_region, 
                               tiles_per_region, objects_per_tile, 
                               num_clusters);
    }

    #pragma myrmics task region inout(top_red_region) \
                         in(per_red_region) safe(per_red_region) \
                         in(num_regions, num_clusters) \
                         safe(num_regions, num_clusters)
    kmeans_myrmics_reduce(top_red_region, per_red_region, 
                          num_regions, num_clusters);
  }


  // Print time and results
  #pragma myrmics task region inout(top_red_region) \
                       in(per_red_region) safe(per_red_region) \
                       in(time_start, num_clusters) \
                       safe(time_start, num_clusters)
  kmeans_myrmics_print(top_red_region, per_red_region, 
                       time_start, num_clusters);

}

