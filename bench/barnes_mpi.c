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
// The code in this file is derived from the MMA version 1.0 project, which 
// was licensed under the following copyright:
//
// Copyright (c) 2011, Lawrence Livermore National Security, LLC.
// Produced at the Lawrence Livermore National Laboratory
// Written by Spyros Lyberis <lymperis1@llnl.gov>
// LLNL-CODE-637218, OCEC-13-184
// All rights reserved.
// 
// This file is part of MMA, version 1.0. 
// For details, please see http://myrmics.com/download.php
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
// - Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the disclaimer below.
// - Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the disclaimer (as noted below) in the
//   documentation and/or other materials provided with the distribution.
// - Neither the name of the LLNS/LLNL nor the names of its contributors may be
//   used to endorse or promote products derived from this software without
//   specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY, LLC,
// THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
// INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
// THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 
// Additional BSD Notice
// 
// 1. This notice is required to be provided under our contract with the U.S.
//    Department of Energy (DOE). This work was produced at Lawrence Livermore
//    National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
// 2. Neither the United States Government nor Lawrence Livermore National
//    Security, LLC nor any of their employees, makes any warranty, express or
//    implied, or assumes any liability or responsibility for the accuracy,
//    completeness, or usefulness of any information, apparatus, product, or
//    process disclosed, or represents that its use would not infringe
//    privately-owned rights.
// 3. Also, reference herein to any specific commercial products, process, or
//    services by trade name, trademark, manufacturer or otherwise does not
//    necessarily constitute or imply its endorsement, recommendation, or
//    favoring by the United States Government or Lawrence Livermore National
//    Security, LLC. The views and opinions of authors expressed herein do not
//    necessarily state or reflect those of the United States Government or
//    Lawrence Livermore National Security, LLC, and shall not be used for
//    advertising or product endorsement purposes.
//
// ============================================================================
// Code adapted for plain MPI by Iakovos Mavroidis
// ============================================================================

#include <arch.h>
#include <kernel_toolset.h>
#include <fmpi.h>

#define COST_EVAL 1            // 0: use plain, non-cost based load balancing
                               // 1: use cost-based load balancing.
                               // Cost-based load balancing distributes bodies
                               // to cores depending on the number of
                               // interactions the bodies had on the previous
                               // step. If they open less cells, they have
                               // less interactions and are thus "cheaper" to
                               // compute than bodies that open many cells.


                               // 1: use cost-based load balancing.
                               // Cost-based load balancing distributes bodies
                               // to cores depending on the number of
                               // interactions the bodies had on the previous
                               // step. If they open less cells, they have
                               // less interactions and are thus "cheaper" to
                               // compute than bodies that open many cells.

#define MAX_IMBALANCE 0.333    // During bisections, maximum imbalance factor.
                               // Select from 0.0 (no imbalance, equal number
                               // of bodies among two sides) up to 1.0 (full
                               // imbalance). Note that value 1.0 is clipped
                               // internally to mean that at least num_cores
                               // bodies will be left to a bisection, so no
                               // cores are left completely without bodies.
                               //
                               // 0.333 means that one side can get up to 
                               // double the other side's bodies.

#define ROUND_ERROR  0.000001F // Amount of floating-point rounding error
                               // we consider "normal" and adjust for when
                               // computing bounding boxes and processor
                               // tree levels exchanges based on them

#define ROUND_ERROR2 0.0001F   // Amount of floating-point rounding error

#define INIT_CRD_MAX 100
#define CRD_MAX 1000
#define VEL_MAX 50
#define MASS_MAX 20
#define MAX_PART_SIZE 800

#define OCT_TREE_INC_SIZE 20	// increment size of oct tree on every realloc

#define THETA 1.0F
#define TIME_STEP 0.125F

// How many Oct Tree levels we descend initially before changing region. Each
// time we descend, this number is halved. E.g. if set to 8, the first 8 Oct
// tree levels are in the first region, the next 4 are in the second region,
// the next 2 in the next etc. When it reaches 1, it stays there (we use one
// region for each additional oct tree level).
#define REG_HEIGHT 4

// Maximum number of regions for an Oct tree
#define MAX_REGIONS 256

// Stack size to operate recursively on the Oct tree. Increase this if the
// problem size is bigger and the assertion fails.
#define MAX_STACK_SIZE 1024
  
#define get_bit(num, pos)       ((num >> (pos)) & 1)
#define toggle_bit(num, pos)    (num ^ (1 << (pos)))


// particle info
typedef struct {
    float   pos[3];
    float   vel[3];
    float   mass;
    unsigned int cost;
} particle_t;

// Structure used for bisection-based load balancing among processors
typedef struct {

  float         dim;            // Particle chosen dimension (x or y or z)
  unsigned long cost;           // Number of force evaluations for particle

} dim_exchange_t;


// OctTree nodes
typedef struct OctTreeNodeStruct OctTreeNode;
struct OctTreeNodeStruct {
  
  // Pointers
//  int   parent;        // OctTree parent node
  int   children[8];   // Children pointers (see above for indexing)
  
  // Values
  float         pos[3];         // If leaf, particle location (x, y, z).
                                // If non-leaf, center of mass location.
  float         vel[3];         // If leaf, particle velocity vector (x, y, z).
                                // If non-leaf, center of mass velocity vector.
  float         mass;           // If leaf, particle mass.
                                // If non-leaf, sum of mass of all children.
  // Flags
  char          is_leaf;        // 1: node is leaf (refers to single particle,
                                //                  no children in array)
                                // 0: node is non-leaf (has some children,
                                //                      subtree has particles)
  char          visited;        // Used during center of mass updates
};


// Per-processor structure. Defines a bounding box and stores the root of the
// OctTree and the region array. All bodies in the local root belong to this
// processor.
typedef struct {

  float         origin[3];      // Bounding box origin location (x, y, z)
  float         size[3];        // Bounding box size (x, y, z)
  int			oct_root;       // OctTree root
  int			ltree_size;		// number of nodes in OctTree
  int           num_levels;     // Number of tree levels

} LocalRoot;



void qsort(dim_exchange_t *dims, int left, int right) {
    dim_exchange_t pivot;
    int tmp, l_hold, r_hold;
  
    l_hold = left;
    r_hold = right;
    pivot.dim = dims[left].dim;
    pivot.cost = dims[left].cost;
    while (left < right)
    {
      while ((dims[right].dim >= pivot.dim) && (left < right))
        right--;
      if (left != right)
      {
        dims[left].dim = dims[right].dim;
        dims[left].cost = dims[right].cost;
        left++;
      }
      while ((dims[left].dim <= pivot.dim) && (left < right))
        left++;
      if (left != right)
      {
        dims[right].dim = dims[left].dim;
        dims[right].cost = dims[left].cost;
        right--;
      }
    }
    dims[left].dim = pivot.dim;
    dims[left].cost = pivot.cost;
    tmp = left;
    left = l_hold;
    right = r_hold;
    if (left < tmp)
      qsort(dims, left, tmp-1);
    if (right > tmp)
      qsort(dims, tmp+1, right);
}

// initialize dims with selected dimension
void barnes_mpi_dims_init(int *dim_choice, particle_t *part, int part_size, dim_exchange_t *dims, int *dims_size) {
    int j;
  
    ar_assert(part_size <= MAX_PART_SIZE);
  
    for (j = 0; j < part_size; j++) {
      dims[j].dim = (*dim_choice == 0) ? part[j].pos[0] :
                       (*dim_choice == 1) ? part[j].pos[1] :
                                               part[j].pos[2];
      dims[j].cost = part[j].cost;
    }

    *dims_size = part_size;
}

// decide dimension to be bisected
void barnes_mpi_dim_choice(float *bbox, int *dim_choice) {
  
    // We'll cut along the dimension which is the longest one
    *dim_choice = 0;
    if (bbox[3] - bbox[2] >
        bbox[2 * *dim_choice + 1] - bbox[2 * *dim_choice]) {
        *dim_choice = 1;
    }
    if (bbox[5] - bbox[4] >
        bbox[2 * *dim_choice + 1] - bbox[2 * *dim_choice]) {
        *dim_choice = 2;
    }
}

// Merge the two bbox
void barnes_mpi_bbox_merge(float *bbox1, float *bbox2) {
    if (bbox1[0] < bbox2[0]) bbox2[0] = bbox1[0];
    if (bbox1[1] > bbox2[1]) bbox2[1] = bbox1[1];
    if (bbox1[2] < bbox2[2]) bbox2[2] = bbox1[2];
    if (bbox1[3] > bbox2[3]) bbox2[3] = bbox1[3];
    if (bbox1[4] < bbox2[4]) bbox2[4] = bbox1[4];
    if (bbox1[5] > bbox2[5]) bbox2[5] = bbox1[5];
}

// merge dims
void barnes_mpi_dims_merge(dim_exchange_t *dims1, int *dims1_size, dim_exchange_t *dims2, int *dims2_size) {
    int i;
  
    //ar_assert(*dims2_size + *dims1_size <= MAX_PART_SIZE);

    for (i = 0; i < *dims1_size; i++) {
        dims2[*dims2_size+i].dim = dims1[i].dim;
        dims2[*dims2_size+i].cost = dims1[i].cost;
    }

    *dims2_size += *dims1_size;
}


void barnes_mpi_bisect(dim_exchange_t *dims, int *dims_size, float *bisect) {
    int i;
    unsigned long total_cost, bisect_cost;
  
    // Sort dimensions 
    qsort(dims, 0, *dims_size-1);
  
    // find total cost and bisection cost
    total_cost = 0;
    for (i = 0; i < *dims_size; i++) {
        total_cost += dims[i].cost;
    }
    bisect_cost = total_cost >> 1;
    // find middle
    total_cost = 0;
    for (i = 0; i < *dims_size; i++) {
        total_cost += dims[i].cost;
        if (total_cost > bisect_cost) {
            *bisect = dims[i].dim;
            return;
        }
    }
}


// swap particles between nodes
void barnes_myrmics_part_swaps(particle_t *part1, int *part_size1, particle_t *part2, int *part_size2, int *dim_choice, float *bisect, int id, int depth) {
	int i;
	int part_pos1, part_pos2;

	ar_assert (*part_size1 + *part_size2 <= MAX_PART_SIZE);

	// first go through part1
	part_pos1 = 0;
	part_pos2 = *part_size2;
	for (i = 0; i < *part_size1; i++) {
		if (part1[i].pos[*dim_choice] < *bisect) {
			// move i to part_pos1 and advance part1_pos
			if (part_pos1 != i) {
				part1[part_pos1].pos[0] = part1[i].pos[0];
				part1[part_pos1].pos[1] = part1[i].pos[1];
				part1[part_pos1].pos[2] = part1[i].pos[2];
				part1[part_pos1].vel[0] = part1[i].vel[0];
				part1[part_pos1].vel[1] = part1[i].vel[1];
				part1[part_pos1].vel[2] = part1[i].vel[2];
				part1[part_pos1].mass = part1[i].mass;
				part1[part_pos1].cost = part1[i].cost;
			}
			part_pos1++;
		} else {
			// move i at the end of part2 and advance part_pos2
			part2[part_pos2].pos[0] = part1[i].pos[0];
			part2[part_pos2].pos[1] = part1[i].pos[1];
			part2[part_pos2].pos[2] = part1[i].pos[2];
			part2[part_pos2].vel[0] = part1[i].vel[0];
			part2[part_pos2].vel[1] = part1[i].vel[1];
			part2[part_pos2].vel[2] = part1[i].vel[2];
			part2[part_pos2].mass = part1[i].mass;
			part2[part_pos2].cost = part1[i].cost;
			part_pos2++;
		}
	}

	// fix sizes
	*part_size1 = part_pos1;
	*part_size2 = part_pos2;

	// now go through part2
	part_pos2 = 0;
	for (i = 0; i < *part_size2; i++) {
		if (part2[i].pos[*dim_choice] < *bisect) {
			// move i to part1_pos and advance part1_pos
			part1[part_pos1].pos[0] = part2[i].pos[0];
			part1[part_pos1].pos[1] = part2[i].pos[1];
			part1[part_pos1].pos[2] = part2[i].pos[2];
			part1[part_pos1].vel[0] = part2[i].vel[0];
			part1[part_pos1].vel[1] = part2[i].vel[1];
			part1[part_pos1].vel[2] = part2[i].vel[2];
			part1[part_pos1].mass = part2[i].mass;
			part1[part_pos1].cost = part2[i].cost;
			part_pos1++;
		} else {
			// move i part_pos2
			if (part_pos2 != i) {
				part2[part_pos2].pos[0] = part2[i].pos[0];
				part2[part_pos2].pos[1] = part2[i].pos[1];
				part2[part_pos2].pos[2] = part2[i].pos[2];
				part2[part_pos2].vel[0] = part2[i].vel[0];
				part2[part_pos2].vel[1] = part2[i].vel[1];
				part2[part_pos2].vel[2] = part2[i].vel[2];
				part2[part_pos2].mass = part2[i].mass;
				part2[part_pos2].cost = part2[i].cost;
			}
			part_pos2++;
		}
	}

	// fix sizes
	*part_size1 = part_pos1;
	*part_size2 = part_pos2;

	/*
	for (i = 0; i < *part_size1; i++) {
		kt_printf("id %d=%f\n\r", id, part1[i].pos[2]);
	}
	kt_printf("id %d---sz1 =%d az2=%d-------------\n\r", id, *part_size1, *part_size2);
	for (i = 0; i < *part_size2; i++) {
		kt_printf("id %d=%f\n\r", id, part2[i].pos[2]);
	}
	*/
}

void barnes_mpi_load_balance(int num_procs, particle_t *part, int *part_size, int *dim_choice, dim_exchange_t *dims, int *dims_size, float *bisect, int worker_id, float *my_minmax, float *foreign_minmax, dim_exchange_t *foreign_dims, int *foreign_dims_size, particle_t *foreign_part, int *foreign_part_size)  {
  int			i, j;
  int			tag_lb = 1;
  MPI_Status	status;

	int depth;
	for (depth = kt_int_log2(num_procs); depth > 0; depth--) {
		// Initialize our min/max array. If we temporarily have no bodies, the
		// dummy limits will be left, which are ok to be combined with other
		// cores (they will be ignored).
		my_minmax[0] = FLT_MAX;
		my_minmax[1] = -FLT_MAX;
		my_minmax[2] = FLT_MAX;
		my_minmax[3] = -FLT_MAX;
		my_minmax[4] = FLT_MAX;
		my_minmax[5] = -FLT_MAX;
		for (j = 0; j < *part_size; j++) {
		if (part[j].pos[0] < my_minmax[0])
			my_minmax[0] = part[j].pos[0];
		if (part[j].pos[0] > my_minmax[1])
	        my_minmax[1] = part[j].pos[0];
		
		if (part[j].pos[1] < my_minmax[2])
			my_minmax[2] = part[j].pos[1];
		if (part[j].pos[1] > my_minmax[3])
			my_minmax[3] = part[j].pos[1];

		if (part[j].pos[2] < my_minmax[4])
			my_minmax[4] = part[j].pos[2];
		if (part[j].pos[2] > my_minmax[5])
			my_minmax[5] = part[j].pos[2];
		}


		//merge bbox between regions up to level depth
		for (i = 0; i <= depth; i++) {
              int start = ((2<<i)>>1)-1;
              int step = (2<<i)>>1;
              for (j = start; j < num_procs; ) {
				if (i < depth) {
                    // merge when we are on intermidiate levels

					if (j == worker_id)  {
						// Send our array to j+step
						MPI_Send(my_minmax, 6, MPI_FLOAT, j+step, tag_lb, MPI_COMM_WORLD);
					}
					if (j + step == worker_id) {
				        // Receive from toggle_bit(worker_id, j)
						MPI_Recv(foreign_minmax, 6, MPI_FLOAT, j, tag_lb, MPI_COMM_WORLD, &status);
						barnes_mpi_bbox_merge(foreign_minmax, my_minmax);
					}
					j += 2*step;
				} else {
					if (j == worker_id) {
                      //nodes on last level decide which bisection will be bisected
                      barnes_mpi_dim_choice(my_minmax, dim_choice);
//                      kt_printf("dim_choice=%d\n\r", *dim_choice);
					}
					j += step;
                  }
              }
		  }

          // broadcast dim_choice
          int src = ((2<<depth)>>1)-1;
          int step = (2<<depth)>>1;
          for (i = 0; i < num_procs; i++) {
				if (i > src) src += step;
				if (i == worker_id) {
					if (worker_id == src) {
						for (j = src-step+1; j < src; j++) {
							MPI_Send(dim_choice, 1, MPI_INT, j, tag_lb, MPI_COMM_WORLD);
						}
					} else {
					MPI_Recv(dim_choice, 1, MPI_INT, src, tag_lb, MPI_COMM_WORLD, &status);
					}
				}
          }

		// Initialize dimensions
		barnes_mpi_dims_init(dim_choice, part, *part_size, dims, dims_size);

		//merge dimensions between regions up to level depth and find bisection
		for (i = 0; i <= depth; i++) {
			int start = ((2<<i)>>1)-1;
			int step = (2<<i)>>1;
			for (j = start; j < num_procs; ) {
				if (i < depth) {
					if (j == worker_id)  {
						// Send our array to j+step
						MPI_Send(dims_size, 1, MPI_INT, j+step, tag_lb, MPI_COMM_WORLD);
						if (*dims_size) MPI_Send(dims, *dims_size*sizeof(dim_exchange_t), MPI_CHAR, j+step, tag_lb, MPI_COMM_WORLD);
					}
					if (j + step == worker_id) {
				        // Receive from toggle_bit(worker_id, j)
						MPI_Recv(foreign_dims_size, 1, MPI_INT, j, tag_lb, MPI_COMM_WORLD, &status);
						if (*foreign_dims_size) {
							MPI_Recv(foreign_dims, *foreign_dims_size*sizeof(dim_exchange_t), MPI_CHAR, j, tag_lb, MPI_COMM_WORLD, &status);
							// merge when we are on intermidiate levels
	                        //printf("merge dims %d %d\n\r", j, j+step);
	                        barnes_mpi_dims_merge(foreign_dims, foreign_dims_size, dims, dims_size);
						}
					}
  
                    j += 2*step;
                  } else {
					if (j == worker_id) {
                      //nodes on last level decide which bisection will be bisected
                      barnes_mpi_bisect(dims, dims_size, bisect);
//                      kt_printf("bisect = %f\n\r", *bisect);
					}
					j += step;
                  }
              }
          }

          // broadcast bisect
          src = ((2<<depth)>>1)-1;
          step = (2<<depth)>>1;
          for (i = 0; i < num_procs; i++) {
                  if (i > src) src += step;
                  if (i == worker_id) {
                      if (worker_id == src) {
                          for (j = src-step+1; j < src; j++) {
                              MPI_Send(bisect, 1, MPI_FLOAT, j, tag_lb, MPI_COMM_WORLD);
                          }
                      } else {
                      MPI_Recv(bisect, 1, MPI_FLOAT, src, tag_lb, MPI_COMM_WORLD, &status);
                      }
                  }
			}


			// transfer particles between processors
			if (get_bit(worker_id, depth-1)) {	
				// send particles to peer
				MPI_Send(part_size, 1, MPI_INT, toggle_bit(worker_id, depth-1), tag_lb, MPI_COMM_WORLD);
				if (*part_size) MPI_Send(part, *part_size*sizeof(particle_t), MPI_CHAR, toggle_bit(worker_id, depth-1), tag_lb, MPI_COMM_WORLD);
				// receive updated particles from peer
				MPI_Recv(part_size, 1, MPI_INT, toggle_bit(worker_id, depth-1), tag_lb, MPI_COMM_WORLD, &status);
				if (*part_size) MPI_Recv(part, *part_size*sizeof(particle_t), MPI_CHAR, toggle_bit(worker_id, depth-1), tag_lb, MPI_COMM_WORLD, &status);
			}
			if (!get_bit(worker_id, depth-1)) {	
				// receive particles from peer
				MPI_Recv(foreign_part_size, 1, MPI_INT, toggle_bit(worker_id, depth-1), tag_lb, MPI_COMM_WORLD, &status);
				if (*foreign_part_size) MPI_Recv(foreign_part, *foreign_part_size*sizeof(particle_t), MPI_CHAR, toggle_bit(worker_id, depth-1), tag_lb, MPI_COMM_WORLD, &status);
				// fix particles
				barnes_myrmics_part_swaps(part, part_size, foreign_part, foreign_part_size, dim_choice, bisect, worker_id, depth);
				// send new particles to peer
				MPI_Send(foreign_part_size, 1, MPI_INT, toggle_bit(worker_id, depth-1), tag_lb, MPI_COMM_WORLD);
				if (*foreign_part_size) MPI_Send(foreign_part, *foreign_part_size*sizeof(particle_t), MPI_CHAR, toggle_bit(worker_id, depth-1), tag_lb, MPI_COMM_WORLD);
			}
				
	}

}

void barnes_mpi_part_init(unsigned int *seed, int part_size, particle_t *part) {
      int j;
  
      // set random coordinates
      for (j = 0; j < part_size; j++) {
          part[j].pos[0] = (float) ((*seed = kt_rand(*seed)) % (20*INIT_CRD_MAX)) / 10.0F - INIT_CRD_MAX;
          part[j].pos[1] = (float) ((*seed = kt_rand(*seed)) % (20*INIT_CRD_MAX)) / 10.0F - INIT_CRD_MAX;
          part[j].pos[2] = (float) ((*seed = kt_rand(*seed)) % (20*INIT_CRD_MAX)) / 10.0F - INIT_CRD_MAX;
          part[j].vel[0] = (float) ((*seed = kt_rand(*seed)) % (20*VEL_MAX)) / 10.0F - VEL_MAX;
          part[j].vel[1] = (float) ((*seed = kt_rand(*seed)) % (20*VEL_MAX)) / 10.0F - VEL_MAX;
          part[j].vel[2] = (float) ((*seed = kt_rand(*seed)) % (20*VEL_MAX)) / 10.0F - VEL_MAX;
          part[j].mass = (float) (1.0F + (*seed = kt_rand(*seed)) % (10*MASS_MAX)) / 10.0F;
          part[j].cost = 1;
      }
}

int create_new_child(LocalRoot *lroot, OctTreeNode **ltree_p, int *ltree_size, int *ltree_tot_size, int cur_oct_pos, int level,
                      int child, OctTreeNode **ret_new_oct, int *ret_level, int id) {

  OctTreeNode   *new_oct;
  int           new_level;
  int			i;
  OctTreeNode   *ltree;
  int r = 0;
 
  ltree = *ltree_p;

  // Change levels
  new_level = level + 1;

  // Create the new, empty Oct node. Oct tree nodes are always part of 
  // a region which is dependent on the oct tree level.
  if (*ltree_size == *ltree_tot_size) {
	*ltree_tot_size += OCT_TREE_INC_SIZE;
	ltree = kt_realloc(ltree, *ltree_tot_size * sizeof(OctTreeNode));
	*ltree_p = ltree;
	r = 1;
  }
  new_oct = &ltree[*ltree_size];
//  new_oct->parent = cur_oct;
  ltree[cur_oct_pos].children[child] = *ltree_size;
  for (i = 0; i < 8; i++) new_oct->children[i] = -1;
  new_oct->is_leaf = 0;
  new_oct->visited = 0;

  *ltree_size += 1;

  // Update max tree levels
  if (lroot->num_levels < new_level + 1) {
    lroot->num_levels = new_level + 1;
  }

  // Return stuff
  if (ret_new_oct) {
    *ret_new_oct = new_oct;
  }
  if (ret_level) {
    *ret_level = new_level;
  }

  return r;
}


// ===========================================================================
// ===========================================================================
void barnes_mpi_create_oct_tree(particle_t *part, int *part_size, LocalRoot *lroot, OctTreeNode **ltree_p, int *ltree_size, int *ltree_tot_size, int id) {
	OctTreeNode   *cur_oct;
	float         min[3];
	float         max[3];
	float         cur_centre[3];
	float         cur_size[3];
	float         new_centre[3];
	float         new_size[3];
	int           level;
	int           child;
	int           child_reloc;
	OctTreeNode   *new_oct;
	OctTreeNode   *new2_oct;
	int           new_level;
	OctTreeNode   *reloc;
	int           i;
	OctTreeNode   *ltree;
	int           cur_oct_pos;
	int           reloc_pos;
 
	ltree = *ltree_p;
 
  // Sanity checks
  ar_assert(*part_size > 0);

  // Find bounding box
  min[0] = min[1] = min[2] = CRD_MAX;
  max[0] = max[1] = max[2] = -CRD_MAX;
  for (i = 0; i < *part_size; i++) {
    if (part[i].pos[0] < min[0]) min[0] = part[i].pos[0];
    if (part[i].pos[1] < min[1]) min[1] = part[i].pos[1];
    if (part[i].pos[2] < min[2]) min[2] = part[i].pos[2];
    if (part[i].pos[0] > max[0]) max[0] = part[i].pos[0];
    if (part[i].pos[1] > max[1]) max[1] = part[i].pos[1];
    if (part[i].pos[2] > max[2]) max[2] = part[i].pos[2];
  }
  ar_assert((min[0] < CRD_MAX) && (min[0] > -CRD_MAX));
  ar_assert((min[1] < CRD_MAX) && (min[1] > -CRD_MAX));
  ar_assert((min[2] < CRD_MAX) && (min[2] > -CRD_MAX));
  ar_assert((max[0] < CRD_MAX) && (max[0] > -CRD_MAX));
  ar_assert((max[1] < CRD_MAX) && (max[1] > -CRD_MAX));
  ar_assert((max[2] < CRD_MAX) && (max[2] > -CRD_MAX));

  // Set bounding box according to min and max values
  lroot->origin[0] = min[0];
  lroot->origin[1] = min[1];
  lroot->origin[2] = min[2];
  lroot->size[0]   = max[0] - min[0];
  lroot->size[1]   = max[1] - min[1];
  lroot->size[2]   = max[2] - min[2];

  // Allocate oct tree root, using the first body
  cur_oct = &ltree[0];
  cur_oct_pos = 0;
  *ltree_size = 1;
  lroot->oct_root = 0;
//  cur_oct->parent = -1;
  for (i = 0; i < 8; i++) cur_oct->children[i] = -1;
  cur_oct->is_leaf = 1;
  cur_oct->visited = 0;
  cur_oct->pos[0]  = part[0].pos[0];
  cur_oct->pos[1]  = part[0].pos[1];
  cur_oct->pos[2]  = part[0].pos[2];
  cur_oct->vel[0]  = part[0].vel[0];
  cur_oct->vel[1]  = part[0].vel[1];
  cur_oct->vel[2]  = part[0].vel[2];
  cur_oct->mass    = part[0].mass;
  lroot->num_levels = 1;

  ar_assert(cur_oct->pos[0] >= lroot->origin[0]);
  ar_assert(cur_oct->pos[0] <= lroot->origin[0] +
                                lroot->size[0]+ROUND_ERROR2);
  ar_assert(cur_oct->pos[1] >= lroot->origin[1]);
  ar_assert(cur_oct->pos[1] <= lroot->origin[1] +
                                lroot->size[1]+ROUND_ERROR2);
  ar_assert(cur_oct->pos[2] >= lroot->origin[2]);
  ar_assert(cur_oct->pos[2] <= lroot->origin[2] +
                                lroot->size[2]+ROUND_ERROR2);

  // Deal with all the rest
  for (i = 1; i < *part_size; i++) {

    // Start from the top
    level = 0;
    cur_oct = &ltree[lroot->oct_root];
	cur_oct_pos = lroot->oct_root;
    cur_size[0] = lroot->size[0];
    cur_size[1] = lroot->size[1];
    cur_size[2] = lroot->size[2];
    cur_centre[0] = lroot->origin[0] + cur_size[0] / 2.0F;
    cur_centre[1] = lroot->origin[1] + cur_size[1] / 2.0F;
    cur_centre[2] = lroot->origin[2] + cur_size[2] / 2.0F;
    reloc = NULL;
	reloc_pos = -1;

    while (1) {

      // Find out into which bisection the body must go       
      child = 0;
      if (part[i].pos[0] >= cur_centre[0]) {
        child += 4;
      }
      if (part[i].pos[1] >= cur_centre[1]) {
        child += 2;
      }
      if (part[i].pos[2] >= cur_centre[2]) {
        child += 1;
      }

      // Compute the new centre of the bounding box for the target bisection
      new_size[0] = cur_size[0] / 2.0F;
      new_size[1] = cur_size[1] / 2.0F;
      new_size[2] = cur_size[2] / 2.0F;
      if (child & 4) {
        new_centre[0] = cur_centre[0] + new_size[0] / 2.0F;
      }
      else {
        new_centre[0] = cur_centre[0] - new_size[0] / 2.0F;
      }
      if (child & 2) {
        new_centre[1] = cur_centre[1] + new_size[1] / 2.0F;
      }
      else {
        new_centre[1] = cur_centre[1] - new_size[1] / 2.0F;
      }
      if (child & 1) {
        new_centre[2] = cur_centre[2] + new_size[2] / 2.0F;
      }
      else {
        new_centre[2] = cur_centre[2] - new_size[2] / 2.0F;
      }

      // If we can descend, do so
      if (cur_oct->children[child] != -1) {
		cur_oct_pos = cur_oct->children[child];
        cur_oct = &ltree[cur_oct->children[child]];

        // Change level. Region must exist for this one.
        level++;

        // Update current bounding box
        cur_size[0] = new_size[0];
        cur_size[1] = new_size[1];
        cur_size[2] = new_size[2];
        cur_centre[0] = new_centre[0];
        cur_centre[1] = new_centre[1];
        cur_centre[2] = new_centre[2];

        // Descend
        continue;
      }

      // We have to split this node. If it's a leaf, remember the oct node,
      // because we have to relocate the body in it.
      if (cur_oct->is_leaf) {
        ar_assert(!reloc);
        reloc = cur_oct;
		reloc_pos = cur_oct_pos;
        cur_oct->is_leaf = 0;
      }

      // Create a new child, possibly creating a new region as well
      if (create_new_child(lroot, ltree_p, ltree_size, ltree_tot_size, cur_oct_pos, level, child, &new_oct, &new_level, id)) {
			// if reallocated new space we have to fix the pointers
			ltree = *ltree_p;
			cur_oct = &ltree[cur_oct_pos];
      }

      // Do we have an old body to relocate?
      if (reloc) {
		reloc = &ltree[reloc_pos];

        // Compute on which bisection it needs to go.
        child_reloc = 0;
        if (reloc->pos[0] >= cur_centre[0]) {
          child_reloc += 4;
        }
        if (reloc->pos[1] >= cur_centre[1]) {
          child_reloc += 2;
        }
        if (reloc->pos[2] >= cur_centre[2]) {
          child_reloc += 1;
        }

        // If they go on a different bisections, leave the new node we created
        // as a leaf for the new body and create a new one for the old body.
        if (child_reloc != child) {

          // Create child for the old body ltree_r, 
          if (create_new_child(lroot, ltree_p, ltree_size, ltree_tot_size, cur_oct_pos, level, child_reloc, &new2_oct, NULL, id)) {
			  // if reallocated new space we have to fix the pointers
			  ltree = *ltree_p;
	          cur_oct = &ltree[cur_oct_pos];
	          reloc = &ltree[reloc_pos];
	          new_oct = &ltree[*ltree_size-2];
          }

          // Copy fields
          new2_oct->is_leaf = 1;
          new2_oct->pos[0] = reloc->pos[0];
          new2_oct->pos[1] = reloc->pos[1];
          new2_oct->pos[2] = reloc->pos[2];
          new2_oct->vel[0] = reloc->vel[0];
          new2_oct->vel[1] = reloc->vel[1];
          new2_oct->vel[2] = reloc->vel[2];
          new2_oct->mass   = reloc->mass;
          cur_oct->visited = reloc->visited;

          ar_assert(new2_oct->pos[0] >= lroot->origin[0]);
          ar_assert(new2_oct->pos[0] <= lroot->origin[0] +
                                        lroot->size[0]+ROUND_ERROR2);
          ar_assert(new2_oct->pos[1] >= lroot->origin[1]);
          ar_assert(new2_oct->pos[1] <= lroot->origin[1] +
                                        lroot->size[1]+ROUND_ERROR2);
          ar_assert(new2_oct->pos[2] >= lroot->origin[2]);
          ar_assert(new2_oct->pos[2] <= lroot->origin[2] +
                                        lroot->size[2]+ROUND_ERROR2);

          // Relocation finished
          reloc = NULL;
        }
      }

      // If we have no relocation pending, the body is put here and the new
      // node is now a leaf.
      if (!reloc) {
        new_oct->is_leaf = 1;
        new_oct->pos[0] = part[i].pos[0];
        new_oct->pos[1] = part[i].pos[1];
        new_oct->pos[2] = part[i].pos[2];
        new_oct->vel[0] = part[i].vel[0];
        new_oct->vel[1] = part[i].vel[1];
        new_oct->vel[2] = part[i].vel[2];
        new_oct->mass   = part[i].mass;
        new_oct->visited= 0;

        ar_assert(new_oct->pos[0] >= lroot->origin[0]);
        ar_assert(new_oct->pos[0] <= lroot->origin[0] +
                                      lroot->size[0]+ROUND_ERROR2);
        ar_assert(new_oct->pos[1] >= lroot->origin[1]);
        ar_assert(new_oct->pos[1] <= lroot->origin[1] +
                                      lroot->size[1]+ROUND_ERROR2);
        ar_assert(new_oct->pos[2] >= lroot->origin[2]);
        ar_assert(new_oct->pos[2] <= lroot->origin[2] +
                                      lroot->size[2]+ROUND_ERROR2);

        // We're done
        break;
      }

      // Otherwise, both the new and the old body need to go further down 
      // into the tree because they are too near to each other for this
      // granularity. Continue descending.

      // Leave the new node as non-leaf
      cur_oct = new_oct;
      cur_oct_pos = *ltree_size-1;
      level = new_level;


      // Update current bounding box
      cur_size[0] = new_size[0];
      cur_size[1] = new_size[1];
      cur_size[2] = new_size[2];
      cur_centre[0] = new_centre[0];
      cur_centre[1] = new_centre[1];
      cur_centre[2] = new_centre[2];
    }
  }


//  barnes_myrmics_update_mass_centers(id, ltree_r, lroot);

}



//////////////////////
////  MASS CENTER ////
//////////////////////

void barnes_mpi_update_mass_centers(LocalRoot *lroot, OctTreeNode *ltree, int *ltree_size, int id) {

  int	stack[MAX_STACK_SIZE];
  int	stack_size;
  int	i;
  OctTreeNode   *node;

  // Initialize the stack with the oct tree root
  stack[0] = lroot->oct_root;
  stack_size = 1;

  // Iterate on the stack
  while (stack_size) {

    // Peek next node
    ar_assert(stack_size > 0);
    node = &ltree[stack[stack_size - 1]];

    // Does it have children?
    if (!node->is_leaf) {

      // If we haven't visited this in the past...
      if (!node->visited) {

        // ... leave it on the stack and add all its children, to be updated
        // before this node.
        for (i = 0; i < 8; i++) {
          if (node->children[i] != -1) {
            stack[stack_size++] = node->children[i];
            ar_assert(stack_size < MAX_STACK_SIZE);
          }
        }

        // Mark that we've visited this one
        node->visited = 1;
      }

      // If we have visited it in the past...
      else {
        // ... its children are updated. Update this node.
        node->pos[0] = 0;
        node->pos[1] = 0;
        node->pos[2] = 0;
        node->mass   = 0;
        // We do it in two steps, because when using large masses and/or
        // positions the floats can overflow
        for (i = 0; i < 8; i++) {
          if (node->children[i] != -1) {
            node->mass   += ltree[node->children[i]].mass;
          }
        }

		ar_assert(node->mass);

        for (i = 0; i < 8; i++) {
          if (node->children[i] != -1) {
            node->pos[0] += ltree[node->children[i]].pos[0] *
                            (ltree[node->children[i]].mass / node->mass);
            node->pos[1] += ltree[node->children[i]].pos[1] *
                            (ltree[node->children[i]].mass / node->mass);
            node->pos[2] += ltree[node->children[i]].pos[2] *
                            (ltree[node->children[i]].mass / node->mass);
          }
        }

        // Sanity check 1: that's a dangerous spot for inf and nan
        ar_assert((node->pos[0] < CRD_MAX) && (node->pos[0] > -CRD_MAX));
        ar_assert((node->pos[1] < CRD_MAX) && (node->pos[1] > -CRD_MAX));
        ar_assert((node->pos[2] < CRD_MAX) && (node->pos[2] > -CRD_MAX));

        // Sanity check 2: it should be part of our bounding box
        ar_assert(node->pos[0] >= lroot->origin[0]);
        ar_assert(node->pos[0] <= lroot->origin[0] + lroot->size[0]+ROUND_ERROR);
        ar_assert(node->pos[1] >= lroot->origin[1]);
        ar_assert(node->pos[1] <= lroot->origin[1] + lroot->size[1]+ROUND_ERROR);
        ar_assert(node->pos[2] >= lroot->origin[2]);
        ar_assert(node->pos[2] <= lroot->origin[2] + lroot->size[2]+ROUND_ERROR);

        // Pop it
        stack_size--;
      }
    }
    else {
      // Make sure it's in our bounding box, and adjust for any rounding
      // errors. We can encounter these even in leaf nodes, due to the
      // original min/max which is converted to origin/size.
      for (i = 0; i < 3; i++) {
        ar_assert(node->pos[i] >= lroot->origin[i]);
        if (node->pos[i] > lroot->origin[i] + lroot->size[i]) {
          ar_assert(node->pos[i] <=
                 lroot->origin[i] + lroot->size[i] * (1.0F + ROUND_ERROR));
          lroot->size[i] *= (1.0F + ROUND_ERROR);
        }
      }

      // We don't care if it's a leaf. Just pop it.
      stack_size--;
    }
  }
}




void barnes_mpi_update_particles_common(LocalRoot **lroot_p, OctTreeNode **ltree_p, int id, int num_bodies, int *part_size, particle_t *part) {

  int			stack1[MAX_STACK_SIZE];
  int           stack1_level[MAX_STACK_SIZE];
  int           stack1_size;
  OctTreeNode   *node1;
  particle_t    *part_n;
  int           level1;
  int			stack2[MAX_STACK_SIZE];
  int           stack2_level[2 * MAX_STACK_SIZE];
  float         stack2_box_origin[3 * MAX_STACK_SIZE];
  float         stack2_box_size[3 * MAX_STACK_SIZE];
  int           stack2_size;
  OctTreeNode   *node2;
  int           level2;
  int           root_id2;
  float         box2_origin[3];
  float         box2_size[3];
  float         dist;
  float         dist_squared;
  float         ft;
  float         f[3];
  float         d[3];
  int           i;
  LocalRoot     *my_root;
  int           part_i;

	
  // Initialize stack1 with the local oct tree root. This will iterate on all
  // the local bodies that we're responsible to update forces and velocities.
  stack1[0] = lroot_p[id]->oct_root;
  stack1_level[0] = 0;
  stack1_size = 1;

  my_root = lroot_p[id];

  // index in particles
  part_i = 0;

  // Iterate on stack1
  while (stack1_size) {

    // Pop next node
    ar_assert(stack1_size > 0);
    node1 = &ltree_p[id][stack1[--stack1_size]];

    level1 = stack1_level[stack1_size];

    // Sanity check: it should be part of our bounding box
    ar_assert(node1->pos[0] >= my_root->origin[0]);
    ar_assert(node1->pos[0] <= my_root->origin[0] +
                            my_root->size[0]+ROUND_ERROR2);
    ar_assert(node1->pos[1] >= my_root->origin[1]);
    ar_assert(node1->pos[1] <= my_root->origin[1] +
                            my_root->size[1]+ROUND_ERROR2);
    ar_assert(node1->pos[2] >= my_root->origin[2]);
    ar_assert(node1->pos[2] <= my_root->origin[2] +
                            my_root->size[2]+ROUND_ERROR2);

    // Does it have children?
    if (!node1->is_leaf) {
      // Push them on the stack
      for (i = 0; i < 8; i++) {
        if (node1->children[i] != -1) {
          stack1_level[stack1_size] = level1 + 1;
          stack1[stack1_size++] = node1->children[i];
          ar_assert(stack1_size < MAX_STACK_SIZE);
          // Check that the children level is sane
          ar_assert(level1 + 1 < my_root->num_levels);
        }
      }
    }

    // Leaf node. This is one of our bodies, so we need to compute all the
    // interactions with all other bodies in the system.
    else {

      // copy particle in part array
      part_n = &part[part_i++];
      part_n->pos[0] = node1->pos[0];
      part_n->pos[1] = node1->pos[1];
      part_n->pos[2] = node1->pos[2];
      part_n->vel[0] = node1->vel[0];
      part_n->vel[1] = node1->vel[1];
      part_n->vel[2] = node1->vel[2];
      part_n->mass = node1->mass;
#if COST_EVAL
      part_n->cost = 0;
#else 
      part_n->cost = 1;
#endif

      // Initialize forces
      f[0] = 0;
      f[1] = 0;
      f[2] = 0;

      // Initialize stack2 with both local and remote oct tree roots. This
      // will iterate on all local and remote bodies that we must check for
      // interactions with node1.
      stack2_size = 0;

      for (i = 0; i < num_bodies; i++) {
        stack2[stack2_size]                  = lroot_p[i]->oct_root;
        stack2_level[2*stack2_size]          = 0; // level
        stack2_level[2*stack2_size + 1]      = i; // root id
        stack2_box_origin[3*stack2_size]     = lroot_p[i]->origin[0];
        stack2_box_origin[3*stack2_size + 1] = lroot_p[i]->origin[1];
        stack2_box_origin[3*stack2_size + 2] = lroot_p[i]->origin[2];
        stack2_box_size[3*stack2_size]       = lroot_p[i]->size[0];
        stack2_box_size[3*stack2_size + 1]   = lroot_p[i]->size[1];
        stack2_box_size[3*stack2_size + 2]   = lroot_p[i]->size[2];
        stack2_size++;
        ar_assert(stack2_size < MAX_STACK_SIZE);

        ar_assert(ltree_p[i][lroot_p[i]->oct_root].pos[0] >= lroot_p[i]->origin[0]);

      }

      // Iterate on stack2
      while (stack2_size) {

        // Pop next node
        ar_assert(stack2_size > 0);
		stack2_size--;
        node2 = &ltree_p[stack2_level[2*stack2_size + 1]][stack2[stack2_size]];
        level2         = stack2_level[2*stack2_size];
        root_id2       = stack2_level[2*stack2_size + 1];
        box2_origin[0] = stack2_box_origin[3*stack2_size];
        box2_origin[1] = stack2_box_origin[3*stack2_size + 1];
        box2_origin[2] = stack2_box_origin[3*stack2_size + 2];
        box2_size[0]   = stack2_box_size[3*stack2_size];
        box2_size[1]   = stack2_box_size[3*stack2_size + 1];
        box2_size[2]   = stack2_box_size[3*stack2_size + 2];

        // Sanity check: it must be inside root_id's bounding box
        ar_assert(node2->pos[0] >= lroot_p[root_id2]->origin[0]);
        ar_assert(node2->pos[0] <= lroot_p[root_id2]->origin[0] +
                                lroot_p[root_id2]->size[0]+ROUND_ERROR2);
        ar_assert(node2->pos[1] >= lroot_p[root_id2]->origin[1]);
        ar_assert(node2->pos[1] <= lroot_p[root_id2]->origin[1] +
                                lroot_p[root_id2]->size[1]+ROUND_ERROR2);
        ar_assert(node2->pos[2] >= lroot_p[root_id2]->origin[2]);
        ar_assert(node2->pos[2] <= lroot_p[root_id2]->origin[2] +
                                lroot_p[root_id2]->size[2]+ROUND_ERROR2);

        // We may bump upon ourselves. Don't self-interact, it's sick.
        if (node2 == node1) {
          continue;
        }

        // Compute node1 - node2 distance per dimension
        d[0] = node2->pos[0] - node1->pos[0];
        d[1] = node2->pos[1] - node1->pos[1];
        d[2] = node2->pos[2] - node1->pos[2];

        // Euclidean distance
        dist_squared = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
        dist = ar_float_sqrt(dist_squared);

        // Decide if the node can be used as-is, or it should be opened.
        if ((!node2->is_leaf) &&
            // Case 1: if node1 is contained in node2 bounding box, node2
            //         should be opened.
            (((node1->pos[0] >= box2_origin[0]) &&
              (node1->pos[0] <  box2_origin[0] + box2_size[0]) &&
              (node1->pos[1] >= box2_origin[1]) &&
              (node1->pos[1] <  box2_origin[1] + box2_size[1]) &&
              (node1->pos[2] >= box2_origin[2]) &&
              (node1->pos[2] <  box2_origin[2] + box2_size[2])) ||
              // Case 2: if node1 is close enough to node2, node2 should
              //         also be opened. The criterion dist <= length / theta
              //         is used here per dimension, because we have non-cubic
              //         bisections.
             //((dist <= (box2_size[0] + box2_size[1] + box2_size[2]) / 
             //          (3.0 * theta))))) 
             ((dist <= box2_size[0] / THETA) ||
              (dist <= box2_size[1] / THETA) ||
              (dist <= box2_size[2] / THETA)))) {

          // Open node2, i.e. push its children to the stack
          for (i = 0; i < 8; i++) {
            if (node2->children[i] != -1) {

              // Push the child
              stack2[stack2_size] = node2->children[i];
              stack2_level[2 * stack2_size] = level2 + 1; // next level
              stack2_level[2 * stack2_size + 1] = root_id2; // same root

              // Check that the child level is sane
              if (root_id2 == id) {
                // Check local level
                ar_assert(level2 + 1 < my_root->num_levels);
              }

              // Child bounding box origin
              if (i & 4) {
                stack2_box_origin[3 * stack2_size] =
                                        box2_origin[0] + box2_size[0] / 4.0F;
              }
              else {
                stack2_box_origin[3 * stack2_size] = box2_origin[0];
              }
              if (i & 2) {
                stack2_box_origin[3 * stack2_size + 1] =
                                        box2_origin[1] + box2_size[1] / 4.0F;
              }
              else {
                stack2_box_origin[3 * stack2_size + 1] = box2_origin[1];
              }
              if (i & 1) {
                stack2_box_origin[3 * stack2_size + 2] =
                                        box2_origin[2] + box2_size[2] / 4.0F;
              }
              else {
                stack2_box_origin[3 * stack2_size + 2] = box2_origin[2];
              }

              // Child bounding box size
              stack2_box_size[3 * stack2_size]     = box2_size[0] / 2.0;
              stack2_box_size[3 * stack2_size + 1] = box2_size[1] / 2.0;
              stack2_box_size[3 * stack2_size + 2] = box2_size[2] / 2.0;

              // Increase stack pointer
              stack2_size++;
              ar_assert(stack2_size < MAX_STACK_SIZE);
            }
          }
        }

        // If it's not opened, we use it as it is. It's already popped from
        // the stack.
        else {

          // Add to node1 (the outer loop body) cost that we did one more
          // force evaluation for its sake
#if COST_EVAL
          part_n->cost++;
#endif

          // Collision or rounding error? Disregard.
          if (dist == 0) {
            continue;
          }

          // Total force. It's difficult to do this with single-point accuracy
          // and still have realistic data models to run. Floats crap out at
          // e38. Typical earth/sun masses are at e24/e30, so the function
          // below is doomed to inf. We still try for this not to happen, by
          // interleaving (probably) small and big quantities. If your data
          // models fail here, try scaling them down.
          ft = ((6.67e-11F * node1->mass) / dist_squared) * node2->mass;

          // Project onto the three dimensions and accumulate to what other
          // bodies have given so far
          f[0] += ft * d[0] / dist;
          f[1] += ft * d[1] / dist;
          f[2] += ft * d[2] / dist;

          // Sanity checks for single-point overflows
          ar_assert((dist < FLT_MAX) && (dist > 0));
          ar_assert((ft   < FLT_MAX) && (ft   > -FLT_MAX));
          ar_assert((f[0] < FLT_MAX) && (f[0] > -FLT_MAX));
          ar_assert((f[1] < FLT_MAX) && (f[1] > -FLT_MAX));
          ar_assert((f[2] < FLT_MAX) && (f[2] > -FLT_MAX));
        }
      }

      // All other bodies accounted for. Adjust our velocity.
      part_n->vel[0] += (f[0] * TIME_STEP) / node1->mass;
      part_n->vel[1] += (f[1] * TIME_STEP) / node1->mass;
      part_n->vel[2] += (f[2] * TIME_STEP) / node1->mass;

      // Update its position
      part_n->pos[0] += part_n->vel[0] * TIME_STEP;
      part_n->pos[1] += part_n->vel[1] * TIME_STEP;
      part_n->pos[2] += part_n->vel[2] * TIME_STEP;

      // Sanity checks
      ar_assert((part_n->vel[0] < FLT_MAX) && (part_n->vel[0] > -FLT_MAX));
      ar_assert((part_n->vel[1] < FLT_MAX) && (part_n->vel[1] > -FLT_MAX));
      ar_assert((part_n->vel[2] < FLT_MAX) && (part_n->vel[2] > -FLT_MAX));
    }
  }

  ar_assert(part_i == *part_size);

}

void barnes_mpi_update_particles(LocalRoot **lroot_p, LocalRoot *lroot, OctTreeNode *ltree, int *ltree_size, particle_t *part, int *part_size, int num_procs, int worker_id) {

	int i;
	OctTreeNode **ltree_p;

	// fix local root in global array of lroots
	kt_memcpy(lroot_p[worker_id], lroot, sizeof(LocalRoot));
	lroot_p[worker_id]->ltree_size = *ltree_size;

	// broadcast local roots
	for (i = 0; i < num_procs; i++) 
		MPI_Bcast(lroot_p[i], sizeof(LocalRoot), MPI_CHAR, i, MPI_COMM_WORLD);

	ar_assert(lroot_p[worker_id]->ltree_size == *ltree_size);

	// broadcast trees 
	ltree_p = kt_malloc(num_procs * sizeof(OctTreeNode *));
	for (i = 0; i < num_procs; i++) 
		ltree_p[i] = kt_malloc(lroot_p[i]->ltree_size * sizeof(OctTreeNode));

	// copy local local tree
	kt_memcpy(ltree_p[worker_id], ltree, (*ltree_size) * sizeof(OctTreeNode));

	// broadcast oct trees
	for (i = 0; i < num_procs; i++) 
		MPI_Bcast(ltree_p[i], lroot_p[i]->ltree_size * sizeof(OctTreeNode), MPI_CHAR, i, MPI_COMM_WORLD);

	// call the real update particles rootine
	barnes_mpi_update_particles_common(lroot_p, ltree_p, worker_id, num_procs, part_size, part);

	// free buffers
	for (i = 0; i < num_procs; i++) {
		kt_free(ltree_p[i]);
	}
	kt_free(ltree_p);

}
  

int barnes_mpi(int num_procs, int num_particles, int reps) {
	int rank;
	unsigned int *seed;
	int rep;
	MPI_Status	status;

    particle_t  *part;      // particles in this region
    int			*part_size;
	int			num_cores;
	int			*dim_choice;
	dim_exchange_t *dims;  // merged dimensions used to find bisection
	int			*dims_size;
	float		*bisect;
	LocalRoot	*lroot;
	OctTreeNode	*ltree;
	int			*ltree_size;
	int			*ltree_tot_size;
	LocalRoot	**lroot_p;
	int			i;
	float         *my_minmax;  
	float         *foreign_minmax;  
	dim_exchange_t   *foreign_dims;
	int			*foreign_dims_size;
	particle_t	*foreign_part;
	int			*foreign_part_size;
#ifdef ARCH_MB
	unsigned int          time_start = 0;
	unsigned int          time_stop;
	unsigned int          time;
#endif


    // Initialize MPI
    //MPI_Init(NULL, NULL);

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

	// Number of workers must be a power of 2
	ar_assert((num_procs & (num_procs - 1)) == 0);

	// if we are not used exit
	if (rank >= num_procs) return 0;

	// initialize particles
	// allocate particles
	part = kt_malloc(MAX_PART_SIZE * sizeof(particle_t));
	part_size = kt_malloc(sizeof(int));
	*part_size = num_particles/num_procs;

	ar_assert(*part_size <= MAX_PART_SIZE);
	
	seed = kt_malloc(sizeof(int));
	if (!rank) {
		*seed = 42;
		barnes_mpi_part_init(seed, *part_size, part);
		if (num_procs > 1) MPI_Send(seed, 1, MPI_INT, rank+1, 0, MPI_COMM_WORLD);
	} else {
		MPI_Recv(seed, 1, MPI_INT, rank-1, 0, MPI_COMM_WORLD, &status);
		barnes_mpi_part_init(seed, *part_size, part);
		if (rank < num_procs-1) MPI_Send(seed, 1, MPI_INT, rank+1, 0, MPI_COMM_WORLD); 
	}

	// allocate dim_choice, dims, bisect, etc.
	dim_choice = kt_malloc(sizeof(int));
//	dims = kt_malloc(MAX_PART_SIZE * sizeof(dim_exchange_t));
	dims = kt_malloc(*part_size * num_procs * sizeof(dim_exchange_t));
	dims_size = kt_malloc(sizeof(int));
	bisect = kt_malloc(sizeof(float));
	
	lroot = kt_malloc(sizeof(LocalRoot));
	ltree = kt_malloc(OCT_TREE_INC_SIZE * sizeof(OctTreeNode));
	ltree_size = kt_malloc(sizeof(int));
	ltree_tot_size = kt_malloc(sizeof(int));
	*ltree_size = 0;
	*ltree_tot_size = OCT_TREE_INC_SIZE;
	
	// following variables are used in load_balance. Allocate once here outside the routine.
	my_minmax = kt_malloc(6 * sizeof(float));
	foreign_minmax = kt_malloc(6 * sizeof(float));
	foreign_dims = kt_malloc(MAX_PART_SIZE * sizeof(dim_exchange_t));
	foreign_dims_size = kt_malloc(sizeof(int));
	foreign_part = kt_malloc(MAX_PART_SIZE * sizeof(particle_t));
	foreign_part_size = kt_malloc(sizeof(int));

	// allocate global pointers to lroots
	// following variable is used in update_particles. Allocate once here outside the routine.
	lroot_p = kt_malloc(num_procs * sizeof(LocalRoot *));
	for (i = 0; i < num_procs; i++) 
		lroot_p[i] = kt_malloc(sizeof(LocalRoot));

	OctTreeNode	**ltree_p;
	ltree_p = &ltree;

	MPI_Barrier(MPI_COMM_WORLD);

	// Keep time
	if (!rank) {
#ifdef ARCH_MB
		time_start = ar_glob_timer_read();
#endif
	}

	for (rep = 0; rep < reps; rep++) {
		// LOAD BALANCING [Computation & Communication]
		//
		// Make cores agree on how to split the whole space into multiple
		// orthogonal recursive bisections; have them exchange particles that
		// comply with the newly agreed bisections.
		barnes_mpi_load_balance(num_procs, part, part_size, dim_choice, dims, dims_size, bisect, rank, my_minmax, foreign_minmax, foreign_dims, foreign_dims_size, foreign_part, foreign_part_size);

		// NEW OCT TREE [Computation]
		//
		// Build a new Oct tree from our local bodies
		barnes_mpi_create_oct_tree(part, part_size, lroot, ltree_p, ltree_size, ltree_tot_size, rank);

		// UPDATE MASS CENTERS [Computation]
		//
		// Walk the Oct tree and find the centers of mass of all non-leaf nodes
		barnes_mpi_update_mass_centers(lroot, ltree, ltree_size, rank);

		// CALCULATE NEW PARTICLES
		//  
		barnes_mpi_update_particles(lroot_p, lroot, ltree, ltree_size, part, part_size, num_procs, rank);

	}

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


	return 0;
}

