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
// Code adapted for Myrmics by Iakovos Mavroidis
// ============================================================================

#define FLAT 0                 // Execute flat code with a single scheduler or
                               // hierarchical code with multiple schedulers.

#if FLAT

#define IGNORE_FSTEP 1         // ignore first step since the particles are not
                               // balanced yet and the results do not scale well 
                               // for a small number of steps.

#define COST_EVAL 1            // 0: use plain, non-cost based load balancing
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

#define INIT_CRD_MAX 200
#define CRD_MAX 1000
#define VEL_MAX 5
#define MASS_MAX 2

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

#define get_bit(num, pos)       ((num >> pos) & 1)
#define toggle_bit(num, pos)    (num ^ (1 << pos))

#include <myrmics.h>

// particle info
typedef struct {
	float	pos[3];
	float	vel[3];
	float	mass;
	unsigned int cost;
} particle_t;

// OctTree nodes
typedef struct OctTreeNodeStruct OctTreeNode;
struct OctTreeNodeStruct {

  // Pointers
//  OctTreeNode   *parent;        // OctTree parent node
  OctTreeNode   *children[8];   // Children pointers (see above for indexing)

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


// Structure used for bisection-based load balancing among processors
typedef struct {
	float         dim;            // Particle chosen dimension (x or y or z)
	unsigned long cost;           // Number of force evaluations for particle
} dim_exchange_t;

// Per-processor structure. Defines a bounding box and stores the root of the
// OctTree and the region array. All bodies in the local root belong to this
// processor.
typedef struct {

  float         origin[3];      // Bounding box origin location (x, y, z)
  float         size[3];        // Bounding box size (x, y, z)
  OctTreeNode   *oct_root;      // OctTree root
  int           num_levels;     // Number of tree levels

} LocalRoot;


// Body of a region
typedef struct {
	particle_t	*part;		// particles in this region
	int		*part_size;

	float	*bbox;
	int		*dim_choice;	// select dimention to bisect

	dim_exchange_t	*dims;	// merged dimensions used to find bisection
	int		*dims_size;

	float	*bisect;		// bisection point using dims

	rid_t	ltree_r;		// includes local tree
	LocalRoot	*lroot;		// local oct tree
} body_t;


//////////////////////
//// LOAD BALANCE ////
//////////////////////
void barnes_myrmics_bbox_init(particle_t *part, int *part_size, float *bbox) {
	int j;

    bbox[0] = CRD_MAX;
    bbox[1] = -CRD_MAX;
    bbox[2] = CRD_MAX;
    bbox[3] = -CRD_MAX;
    bbox[4] = CRD_MAX;
    bbox[5] = -CRD_MAX;
    for (j = 0; j < *part_size; j++) {
      if (part[j].pos[0] < bbox[0])
        bbox[0] = part[j].pos[0];
      if (part[j].pos[0] > bbox[1])
        bbox[1] = part[j].pos[0];

      if (part[j].pos[1] < bbox[2])
        bbox[2] = part[j].pos[1];
      if (part[j].pos[1] > bbox[3])
        bbox[3] = part[j].pos[1];

      if (part[j].pos[2] < bbox[4])
        bbox[4] = part[j].pos[2];
      if (part[j].pos[2] > bbox[5])
        bbox[5] = part[j].pos[2];
    }


}

// Merge the two bbox
void barnes_myrmics_bbox_merge(float *bbox1, float *bbox2) {
	if (bbox1[0] < bbox2[0]) bbox2[0] = bbox1[0];
	if (bbox1[1] > bbox2[1]) bbox2[1] = bbox1[1];
	if (bbox1[2] < bbox2[2]) bbox2[2] = bbox1[2];
	if (bbox1[3] > bbox2[3]) bbox2[3] = bbox1[3];
	if (bbox1[4] < bbox2[4]) bbox2[4] = bbox1[4];
	if (bbox1[5] > bbox2[5]) bbox2[5] = bbox1[5];
}

// decide dimension to be bisected
void barnes_myrmics_dim_choice(float *bbox, int *dim_choice) {

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
	//printf("dim_choice=%d\n\r", *dim_choice);
}


void barnes_myrmics_part_init(unsigned int *seed, int part_size_init, particle_t *part, int *part_size) {
	int j;

	*part_size = part_size_init;

	// set random coordinates
	for (j = 0; j < part_size_init; j++) {
		part[j].pos[0] = (float) ((*seed = rand(*seed)) % (20*INIT_CRD_MAX)) / 10.0F - INIT_CRD_MAX;
		part[j].pos[1] = (float) ((*seed = rand(*seed)) % (20*INIT_CRD_MAX)) / 10.0F - INIT_CRD_MAX; 
		part[j].pos[2] = (float) ((*seed = rand(*seed)) % (20*INIT_CRD_MAX)) / 10.0F - INIT_CRD_MAX;
		part[j].vel[0] = (float) ((*seed = rand(*seed)) % (20*VEL_MAX)) / 10.0F - VEL_MAX;
		part[j].vel[1] = (float) ((*seed = rand(*seed)) % (20*VEL_MAX)) / 10.0F - VEL_MAX; 
		part[j].vel[2] = (float) ((*seed = rand(*seed)) % (20*VEL_MAX)) / 10.0F - VEL_MAX;
		part[j].mass = (float) (1.0F + (*seed = rand(*seed)) % (10*MASS_MAX)) / 10.0F;
		part[j].cost = 1;
	}
}

// initialize dims with selected dimension
void barnes_myrmics_dims_init(int *dim_choice, particle_t *part, int *part_size, dim_exchange_t *dims, int *dims_size) {
	int j;

    for (j = 0; j < *part_size; j++) {
      dims[j].dim = (*dim_choice == 0) ? part[j].pos[0] :
                       (*dim_choice == 1) ? part[j].pos[1] :
                                               part[j].pos[2];
      dims[j].cost = part[j].cost;
    }

	*dims_size = *part_size;
}

// merge dims
void barnes_myrmics_dims_merge(dim_exchange_t *dims1, int *dims1_size, dim_exchange_t *dims2, int *dims2_size) {
	int i;

	for (i = 0; i < *dims1_size; i++) {
		dims2[*dims2_size+i].dim = dims1[i].dim;
		dims2[*dims2_size+i].cost = dims1[i].cost;
	}

	*dims2_size += *dims1_size;
}


void qsort(dim_exchange_t *dims, int left, int right)
{
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

void barnes_myrmics_bisect(dim_exchange_t *dims, int *dims_size, float *bisect) {
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
//	printf("bisect = %f\n\t", *bisect);
			return;
		}
	}
}


void barnes_myrmics_part_swaps(particle_t *part1, int *part_size1, particle_t *part2, int *part_size2, int *dim_choice, float *bisect, int id, int depth, int max_part_size) {
	int i;
	int part_pos1, part_pos2;


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
			sys_assert (part_pos2+1 <= max_part_size);
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
			sys_assert (part_pos1+1 <= max_part_size);
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
			sys_assert (part_pos2+1 <= max_part_size);
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

void barnes_myrmics_load_balance(body_t **body, int num_bodies, int max_part_size) {
	int i, j;

	// For each Orthogonal Recursive Bisection phase
	int depth;
	for (depth = int_log2(num_bodies); depth > 0; depth--) {

		//initialize bbox according to particles' coordinates
		for (i = 0; i < num_bodies; i++) {
			particle_t	*scoop_part = body[i]->part;
			int			*scoop_part_size = body[i]->part_size;
			float		*scoop_bbox = body[i]->bbox;
			#pragma myrmics task in(scoop_part, scoop_part_size) out(scoop_bbox)
			barnes_myrmics_bbox_init(scoop_part, scoop_part_size, scoop_bbox);
		}

		//merge bbox between regions up to level depth
		for (i = 0; i <= depth; i++) {
			int start = ((2<<i)>>1)-1;
			int step = (2<<i)>>1;
			for (j = start; j < num_bodies; ) {
				float	*scoop_bbox1 = body[j]->bbox;
				if (i < depth) {
					float	*scoop_bbox2 = body[j+step]->bbox;
					// merge when we are on intermidiate levels
					//printf("merge bbox %d %d\n\r", j, j+step);
					#pragma myrmics task in(scoop_bbox1) inout(scoop_bbox2)
					barnes_myrmics_bbox_merge(scoop_bbox1, scoop_bbox2);
					j += 2*step;
				} else {
					//printf("decide %d\n\r", j, j+step);
					int	*scoop_dim_choice = body[j]->dim_choice;
					//nodes on last level decide which bisection will be bisected
					#pragma myrmics task in(scoop_bbox1) out(scoop_dim_choice)
					barnes_myrmics_dim_choice(scoop_bbox1, scoop_dim_choice);
					j += step;
				}
			}
		}

		// initialize dimensions
		int src = ((2<<depth)>>1)-1;
		int step = (2<<depth)>>1;
		for (i = 0; i < num_bodies; i++) {
			particle_t	*scoop_part = body[i]->part;
			int	*scoop_part_size = body[i]->part_size;
			if (i > src) src += step;
			int	*scoop_dim_choice = body[src]->dim_choice;
			dim_exchange_t	*scoop_dims = body[i]->dims;
			int	*scoop_dims_size = body[i]->dims_size;
			#pragma myrmics task in(scoop_dim_choice, scoop_part, scoop_part_size) out(scoop_dims, scoop_dims_size)
			barnes_myrmics_dims_init(scoop_dim_choice, scoop_part, scoop_part_size, scoop_dims, scoop_dims_size);
		}

		//merge dimensions between regions up to level depth and find bisection
		for (i = 0; i <= depth; i++) {
			int start = ((2<<i)>>1)-1;
			int step = (2<<i)>>1;
			for (j = start; j < num_bodies; ) {
				dim_exchange_t	*scoop_dims1 = body[j]->dims;
				int		*scoop_dims1_size = body[j]->dims_size;
				if (i < depth) {
					dim_exchange_t	*scoop_dims2 = body[j+step]->dims;
					int		*scoop_dims2_size = body[j+step]->dims_size;
					// merge when we are on intermidiate levels
					//printf("merge dims %d %d\n\r", j, j+step);
					#pragma myrmics task in(scoop_dims1, scoop_dims1_size) inout(scoop_dims2, scoop_dims2_size)
					barnes_myrmics_dims_merge(scoop_dims1, scoop_dims1_size, scoop_dims2, scoop_dims2_size);

					j += 2*step;
				} else {
					//printf("bisect %d\n\r", j, j+step);
					float	*scoop_bisect = body[j]->bisect;
					//nodes on last level decide which bisection will be bisected
					#pragma myrmics task in(scoop_dims1, scoop_dims1_size) out(scoop_bisect)
					barnes_myrmics_bisect(scoop_dims1, scoop_dims1_size, scoop_bisect);

					j += step;
				}
			}
		}

		// transfer particles between regions
		int src = ((2<<depth)>>1)-1;
		int step = (2<<depth)>>1;
		for (i = 0; i < num_bodies; i++) {
			if (!get_bit(i, depth-1)) {
				particle_t	*scoop_part1 = body[i]->part;
				int			*scoop_part_size1 = body[i]->part_size;
				particle_t	*scoop_part2 = body[toggle_bit(i, depth-1)]->part;
				int			*scoop_part_size2 = body[toggle_bit(i, depth-1)]->part_size;
				if (i > src) src += step;
				int			*scoop_dim_choice = body[src]->dim_choice;
				float		*scoop_bisect = body[src]->bisect;
				//printf("bisect %d %d\n\r", i, src);
				//printf("swap particles between %d to %d\n\r", i, toggle_bit(i, depth-1));
				#pragma myrmics task inout(scoop_part1, scoop_part_size1, scoop_part2, scoop_part_size2) in(scoop_dim_choice, scoop_bisect, i, depth, max_part_size)
				barnes_myrmics_part_swaps(scoop_part1, scoop_part_size1, scoop_part2, scoop_part_size2, scoop_dim_choice, scoop_bisect, i, depth, max_part_size);
			}
		}

	}
}

//////////////////////
////    OCT TREE  ////
//////////////////////

int level_to_region(int level) {

  int reg = 0;
  int reg_rem_levels = REG_HEIGHT - 1;
  int reg_limit = REG_HEIGHT;
  int i;

  for (i = 0; i < level; i++) {
    if (reg_rem_levels) {
      reg_rem_levels--;
    }
    else {
      if (reg_limit > 1) {
        reg_limit /= 2;
      }
      reg_rem_levels = reg_limit - 1;
      reg++;
    }
  }

  return reg;
}

// ===========================================================================
// ===========================================================================
void create_new_child(LocalRoot *lroot, rid_t ltree_r, OctTreeNode *cur_oct, int level,
                      int child, OctTreeNode **ret_new_oct, int *ret_level) {

  OctTreeNode   *new_oct;
  int           new_level;


  // Change levels
  new_level = level + 1;

  // Create the new, empty Oct node. Oct tree nodes are always part of 
  // a region which is dependent on the oct tree level.
  sys_assert(new_oct = sys_alloc(sizeof(OctTreeNode), ltree_r));
//  new_oct->parent = cur_oct;
  cur_oct->children[child] = new_oct;
  memset(new_oct->children, 0, 8 * sizeof(OctTreeNode *));
  new_oct->is_leaf = 0;
  new_oct->visited = 0;

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
}


// ===========================================================================
// ===========================================================================
void barnes_myrmics_create_oct_tree(int id, particle_t *part, int *part_size, rid_t ltree_r, LocalRoot *lroot) {
	rid_t         tmp;
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

 
  // Sanity checks
  sys_assert(*part_size > 0);

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
  sys_assert((min[0] < CRD_MAX) && (min[0] > -CRD_MAX));
  sys_assert((min[1] < CRD_MAX) && (min[1] > -CRD_MAX));
  sys_assert((min[2] < CRD_MAX) && (min[2] > -CRD_MAX));
  sys_assert((max[0] < CRD_MAX) && (max[0] > -CRD_MAX));
  sys_assert((max[1] < CRD_MAX) && (max[1] > -CRD_MAX));
  sys_assert((max[2] < CRD_MAX) && (max[2] > -CRD_MAX));

  // Set bounding box according to min and max values
  lroot->origin[0] = min[0];
  lroot->origin[1] = min[1];
  lroot->origin[2] = min[2];
  lroot->size[0]   = max[0] - min[0];
  lroot->size[1]   = max[1] - min[1];
  lroot->size[2]   = max[2] - min[2];

  // Allocate oct tree root, using the first body
  sys_assert(cur_oct = sys_alloc(sizeof(OctTreeNode), ltree_r));
  lroot->oct_root = cur_oct;
//  cur_oct->parent = NULL;
  memset(cur_oct->children, 0, 8 * sizeof(OctTreeNode *));
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

  sys_assert(cur_oct->pos[0] >= lroot->origin[0]);
  sys_assert(cur_oct->pos[0] <= lroot->origin[0] +
                                lroot->size[0]+ROUND_ERROR2);
  sys_assert(cur_oct->pos[1] >= lroot->origin[1]);
  sys_assert(cur_oct->pos[1] <= lroot->origin[1] +
                                lroot->size[1]+ROUND_ERROR2);
  sys_assert(cur_oct->pos[2] >= lroot->origin[2]);
  sys_assert(cur_oct->pos[2] <= lroot->origin[2] +
                                lroot->size[2]+ROUND_ERROR2);

  // Deal with all the rest
  for (i = 1; i < *part_size; i++) {

    // Start from the top
    level = 0;
    cur_oct = lroot->oct_root;
    cur_size[0] = lroot->size[0];
    cur_size[1] = lroot->size[1];
    cur_size[2] = lroot->size[2];
    cur_centre[0] = lroot->origin[0] + cur_size[0] / 2.0F;
    cur_centre[1] = lroot->origin[1] + cur_size[1] / 2.0F;
    cur_centre[2] = lroot->origin[2] + cur_size[2] / 2.0F;
    reloc = NULL;

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
      if (cur_oct->children[child]) {
        cur_oct = cur_oct->children[child];

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
        sys_assert(!reloc);
        reloc = cur_oct;
        cur_oct->is_leaf = 0;
      }

      // Create a new child, possibly creating a new region as well
      create_new_child(lroot, ltree_r, cur_oct, level, child, &new_oct, &new_level);

      // Do we have an old body to relocate?
      if (reloc) {

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
          create_new_child(lroot, ltree_r, cur_oct, level, child_reloc, &new2_oct, NULL);

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

          sys_assert(new2_oct->pos[0] >= lroot->origin[0]);
          sys_assert(new2_oct->pos[0] <= lroot->origin[0] +
                                        lroot->size[0]+ROUND_ERROR2);
          sys_assert(new2_oct->pos[1] >= lroot->origin[1]);
          sys_assert(new2_oct->pos[1] <= lroot->origin[1] +
                                        lroot->size[1]+ROUND_ERROR2);
          sys_assert(new2_oct->pos[2] >= lroot->origin[2]);
          sys_assert(new2_oct->pos[2] <= lroot->origin[2] +
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

        sys_assert(new_oct->pos[0] >= lroot->origin[0]);
        sys_assert(new_oct->pos[0] <= lroot->origin[0] +
                                      lroot->size[0]+ROUND_ERROR2);
        sys_assert(new_oct->pos[1] >= lroot->origin[1]);
        sys_assert(new_oct->pos[1] <= lroot->origin[1] +
                                      lroot->size[1]+ROUND_ERROR2);
        sys_assert(new_oct->pos[2] >= lroot->origin[2]);
        sys_assert(new_oct->pos[2] <= lroot->origin[2] +
                                      lroot->size[2]+ROUND_ERROR2);

        // We're done
        break;
      }

      // Otherwise, both the new and the old body need to go further down 
      // into the tree because they are too near to each other for this
      // granularity. Continue descending.

      // Leave the new node as non-leaf
      cur_oct = new_oct;
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

void barnes_myrmics_update_mass_centers(int id, rid_t ltree_r, LocalRoot *lroot) {

  OctTreeNode   *stack[MAX_STACK_SIZE];
  int           stack_size;
  int           i;
  OctTreeNode   *node;


  // Initialize the stack with the oct tree root
  stack[0] = lroot->oct_root;
  stack_size = 1;

  // Iterate on the stack
  while (stack_size) {

    // Peek next node
    sys_assert(stack_size > 0);
    node = stack[stack_size - 1];

    // Does it have children?
    if (!node->is_leaf) {

      // If we haven't visited this in the past...
      if (!node->visited) {

        // ... leave it on the stack and add all its children, to be updated
        // before this node.
        for (i = 0; i < 8; i++) {
          if (node->children[i]) {
            stack[stack_size++] = node->children[i];
            sys_assert(stack_size < MAX_STACK_SIZE);
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
          if (node->children[i]) {
            node->mass   += node->children[i]->mass;
          }
        }
        for (i = 0; i < 8; i++) {
          if (node->children[i]) {
            node->pos[0] += node->children[i]->pos[0] *
                            (node->children[i]->mass / node->mass);
            node->pos[1] += node->children[i]->pos[1] *
                            (node->children[i]->mass / node->mass);
            node->pos[2] += node->children[i]->pos[2] *
                            (node->children[i]->mass / node->mass);
          }
        }

        // Sanity check 1: that's a dangerous spot for inf and nan
        sys_assert((node->pos[0] < CRD_MAX) && (node->pos[0] > -CRD_MAX));
        sys_assert((node->pos[1] < CRD_MAX) && (node->pos[1] > -CRD_MAX));
        sys_assert((node->pos[2] < CRD_MAX) && (node->pos[2] > -CRD_MAX));

        // Sanity check 2: it should be part of our bounding box
        sys_assert(node->pos[0] >= lroot->origin[0]);
        sys_assert(node->pos[0] <= lroot->origin[0] + lroot->size[0]+ROUND_ERROR);
        sys_assert(node->pos[1] >= lroot->origin[1]);
        sys_assert(node->pos[1] <= lroot->origin[1] + lroot->size[1]+ROUND_ERROR);
        sys_assert(node->pos[2] >= lroot->origin[2]);
        sys_assert(node->pos[2] <= lroot->origin[2] + lroot->size[2]+ROUND_ERROR);

        // Pop it
        stack_size--;
      }
    }
    else {
      // Make sure it's in our bounding box, and adjust for any rounding
      // errors. We can encounter these even in leaf nodes, due to the
      // original min/max which is converted to origin/size.
      for (i = 0; i < 3; i++) {
        sys_assert(node->pos[i] >= lroot->origin[i]);
        if (node->pos[i] > lroot->origin[i] + lroot->size[i]) {
          sys_assert(node->pos[i] <=
                 lroot->origin[i] + lroot->size[i] * (1.0F + ROUND_ERROR));
          lroot->size[i] *= (1.0F + ROUND_ERROR);
        }
      }

      // We don't care if it's a leaf. Just pop it.
      stack_size--;
    }
  }
}



void barnes_myrmics_update_particles(rid_t roots_r, rid_t trees_r, LocalRoot **lroot_p, int id, int num_bodies, int *part_size, particle_t *part) {

  OctTreeNode   *stack1[MAX_STACK_SIZE];
  int           stack1_level[MAX_STACK_SIZE];
  int           stack1_size;
  OctTreeNode   *node1;
  particle_t    *part_n;
  int           level1;
  OctTreeNode   *stack2[MAX_STACK_SIZE];
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
    sys_assert(stack1_size > 0);
    node1 = stack1[--stack1_size];

    level1 = stack1_level[stack1_size];

    // Sanity check: it should be part of our bounding box
    sys_assert(node1->pos[0] >= my_root->origin[0]);
    sys_assert(node1->pos[0] <= my_root->origin[0] +
                            my_root->size[0]+ROUND_ERROR2);
    sys_assert(node1->pos[1] >= my_root->origin[1]);
    sys_assert(node1->pos[1] <= my_root->origin[1] +
                            my_root->size[1]+ROUND_ERROR2);
    sys_assert(node1->pos[2] >= my_root->origin[2]);
    sys_assert(node1->pos[2] <= my_root->origin[2] +
                            my_root->size[2]+ROUND_ERROR2);

    // Does it have children?
    if (!node1->is_leaf) {
      // Push them on the stack
      for (i = 0; i < 8; i++) {
        if (node1->children[i]) {
          stack1_level[stack1_size] = level1 + 1;
          stack1[stack1_size++] = node1->children[i];
          sys_assert(stack1_size < MAX_STACK_SIZE);
          // Check that the children level is sane
          sys_assert(level1 + 1 < my_root->num_levels);
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
        sys_assert(stack2_size < MAX_STACK_SIZE);

        sys_assert(lroot_p[i]->oct_root->pos[0] >= lroot_p[i]->origin[0]);

      }

      // Iterate on stack2
      while (stack2_size) {

        // Pop next node
        sys_assert(stack2_size > 0);
        node2 = stack2[--stack2_size];
        level2         = stack2_level[2*stack2_size];
        root_id2       = stack2_level[2*stack2_size + 1];
        box2_origin[0] = stack2_box_origin[3*stack2_size];
        box2_origin[1] = stack2_box_origin[3*stack2_size + 1];
        box2_origin[2] = stack2_box_origin[3*stack2_size + 2];
        box2_size[0]   = stack2_box_size[3*stack2_size];
        box2_size[1]   = stack2_box_size[3*stack2_size + 1];
        box2_size[2]   = stack2_box_size[3*stack2_size + 2];

        // Sanity check: it must be inside root_id's bounding box
        sys_assert(node2->pos[0] >= lroot_p[root_id2]->origin[0]);
        sys_assert(node2->pos[0] <= lroot_p[root_id2]->origin[0] +
                                lroot_p[root_id2]->size[0]+ROUND_ERROR2);
        sys_assert(node2->pos[1] >= lroot_p[root_id2]->origin[1]);
        sys_assert(node2->pos[1] <= lroot_p[root_id2]->origin[1] +
                                lroot_p[root_id2]->size[1]+ROUND_ERROR2);
        sys_assert(node2->pos[2] >= lroot_p[root_id2]->origin[2]);
        sys_assert(node2->pos[2] <= lroot_p[root_id2]->origin[2] +
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
        dist = sqrt(dist_squared);

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
            if (node2->children[i]) {

              // Push the child
              stack2[stack2_size] = node2->children[i];
              stack2_level[2 * stack2_size] = level2 + 1; // next level
              stack2_level[2 * stack2_size + 1] = root_id2; // same root

              // Check that the child level is sane
              if (root_id2 == id) {
                // Check local level
                sys_assert(level2 + 1 < my_root->num_levels);
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
              sys_assert(stack2_size < MAX_STACK_SIZE);
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
          sys_assert((dist < FLT_MAX) && (dist > 0));
          sys_assert((ft   < FLT_MAX) && (ft   > -FLT_MAX));
          sys_assert((f[0] < FLT_MAX) && (f[0] > -FLT_MAX));
          sys_assert((f[1] < FLT_MAX) && (f[1] > -FLT_MAX));
          sys_assert((f[2] < FLT_MAX) && (f[2] > -FLT_MAX));
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
      sys_assert((part_n->vel[0] < FLT_MAX) && (part_n->vel[0] > -FLT_MAX));
      sys_assert((part_n->vel[1] < FLT_MAX) && (part_n->vel[1] > -FLT_MAX));
      sys_assert((part_n->vel[2] < FLT_MAX) && (part_n->vel[2] > -FLT_MAX));
    }
  }

  sys_assert(part_i == *part_size);

}

void barnes_myrmics_free_trees(rid_t roots_r, rid_t trees_r) {
  sys_rfree(trees_r);
}

void barnes_myrmics_phase(rid_t r, rid_t roots_r, body_t **body, LocalRoot **lroot_p, unsigned int *time_start, int step, int num_bodies, int max_part_size) {
	int i;
	rid_t trees_r;


#if IGNORE_FSTEP
	if (step == 1) *time_start = sys_free_timer_get_ticks();
#endif

    // LOAD BALANCING [Computation & Communication]
    //
    // Make cores agree on how to split the whole space into multiple
    // orthogonal recursive bisections; have them exchange particles that
    // comply with the newly agreed bisections.
	//
	// Load balance is slow. Perform once every 4 steps. Effectiveness 
	// of load balance depends on the VEL_MAX, MASS_MAX: if they are small
	// then load balance does not have to be called so often.
	if (!(step & 3)) barnes_myrmics_load_balance(body, num_bodies, max_part_size);

	// allocate all trees
	trees_r = sys_ralloc(r, 0); // this region includes all ltree_r's
	for (i = 0; i < num_bodies; i++) {
		// allocate local tree region
		body[i]->ltree_r = sys_ralloc(trees_r, 0);
	}

	for (i = 0; i < num_bodies; i++) {
		rid_t		scoop_ltree_r = body[i]->ltree_r;
		particle_t	*scoop_part = body[i]->part;
		int			*scoop_part_size = body[i]->part_size;
		LocalRoot	*scoop_lroot = body[i]->lroot;

		// NEW OCT TREE [Computation]
	    //
	    // Build a new Oct tree from our local bodies
		#pragma myrmics task in(i, scoop_part, scoop_part_size) region inout(scoop_ltree_r) inout(scoop_lroot) 
		barnes_myrmics_create_oct_tree(i, scoop_part, scoop_part_size, scoop_ltree_r, scoop_lroot);

	    // UPDATE MASS CENTERS [Computation]
	    //
	    // Walk the Oct tree and find the centers of mass of all non-leaf nodes
		#pragma myrmics task in(i) region inout(scoop_ltree_r) in(scoop_lroot) 
		barnes_myrmics_update_mass_centers(i, scoop_ltree_r, scoop_lroot);
	}

    // CALCULATE NEW PARTICLES
    //
	for (i = 0; i < num_bodies;i++) {
		particle_t	*scoop_part = body[i]->part;
		int			*scoop_part_size = body[i]->part_size;
		// scoop_part_size used for sanity check
		#pragma myrmics task region in(roots_r, trees_r) in(lroot_p, i, num_bodies, scoop_part_size) out(scoop_part)
		barnes_myrmics_update_particles(roots_r, trees_r, lroot_p, i, num_bodies, scoop_part_size, scoop_part);
	}

	// free trees
	#pragma myrmics task region inout(roots_r) in(trees_r) safe(trees_r)
	barnes_myrmics_free_trees(roots_r, trees_r);

}

// ===========================================================================
// ===========================================================================
void barnes_myrmics_time_start(rid_t unused, unsigned int *time_start) {
	*time_start = sys_free_timer_get_ticks();
}

void barnes_myrmics_time_stop(rid_t unused, unsigned int *time_start) {
  unsigned int          time_stop;
  unsigned int          time;

  time_stop = sys_free_timer_get_ticks();
  if (time_stop > *time_start) {
    time = time_stop - *time_start;
  }
  else {
    time = 0xFFFFFFFF - (*time_start - time_stop);
  }

  printf("Time: %10u cycles (%6u msec)\r\n", time, time / 10000);
}


void barnes_myrmics_free_all(rid_t r, rid_t roots_r, body_t **body, LocalRoot **lroot_p, int num_bodies) {
	int i;

	// Free bodies
	for (i = 0; i < num_bodies; i++) {
		sys_free(body[i]->part);
		sys_free(body[i]->part_size);
		sys_free(body[i]->bbox);
		sys_free(body[i]->dim_choice);
		sys_free(body[i]->dims);
		sys_free(body[i]->dims_size);
		sys_free(body[i]->bisect);
		sys_free(body[i]->lroot);
		sys_free(body[i]);
	}
}

void barnes_myrmics(int num_bodies, int num_particles, int steps) {
	rid_t	r;
	rid_t	roots_r;
	LocalRoot **lroot_p;
	body_t	**body;
	int	particles_per_region;
	int i, j;
	unsigned int	*seed;
	unsigned int	*time_start;
	int				max_part_size = (num_particles / num_bodies)*3;

	// Create all-holding body
	r = sys_ralloc(0, 99); // highest level

	// initialize variables
	particles_per_region = num_particles/num_bodies;
	sys_assert(particles_per_region <= max_part_size);
	seed = sys_alloc(sizeof(unsigned int), r);
	*seed = 42;
	time_start = sys_alloc(sizeof(unsigned int), r);

	// Sanity checks
	if (num_bodies & (num_bodies - 1)) {
		printf("num_bodies must be a power of 2\r\n");
		return;
	}

	lroot_p = sys_alloc(num_bodies * sizeof(LocalRoot *), r); // points to lroots
	roots_r = sys_ralloc(r, 0); // this region includes all lroot's

	// Create regions with bodies
	body = sys_alloc(num_bodies * sizeof(body_t *), r);
	for (i = 0; i < num_bodies; i++) {
		// put body in region
		body[i] = sys_alloc(sizeof(body_t), r);
		// allocate particles 
		body[i]->part = sys_alloc(max_part_size * sizeof(particle_t), r);
		body[i]->part_size = sys_alloc(sizeof(int), r);
		particle_t	*scoop_part = body[i]->part;
		int	*scoop_part_size = body[i]->part_size;

		// initialize particles
		#pragma myrmics task inout(seed) in(particles_per_region) out(scoop_part, scoop_part_size)
		barnes_myrmics_part_init(seed, particles_per_region, scoop_part, scoop_part_size);
		
		// allocate bbox
		body[i]->bbox = sys_alloc(6 * sizeof(float), r);

		// allocate dim_choice, dims, bisect, etc.
		body[i]->dim_choice = sys_alloc(sizeof(int), r);
		body[i]->dims = sys_alloc(particles_per_region * num_bodies * sizeof(dim_exchange_t), r);
		body[i]->dims_size = sys_alloc(sizeof(int), r);
		body[i]->bisect = sys_alloc(sizeof(float), r);
		body[i]->lroot = sys_alloc(sizeof(LocalRoot), roots_r);
		lroot_p[i] = body[i]->lroot;
	}

	// Print we're starting
    printf("Barnes Hut of %d particles (%d per region) splitted into %d regions; running %d steps.\r\n",
           num_particles, particles_per_region, num_bodies, steps);

    // Start time
    #pragma myrmics task region inout(r) in(time_start) safe(time_start)
    barnes_myrmics_time_start(r, time_start);

	// Run main loop sequentially
	int step;
	for (step = 0; step < steps; step++) {
		#pragma myrmics task region inout(r) in(roots_r, body, lroot_p, time_start) safe(roots_r, body, lroot_p, time_start) in(step, num_bodies, max_part_size)
		barnes_myrmics_phase(r, roots_r, body, lroot_p, time_start, step, num_bodies, max_part_size);

	}

    // Stop time
    #pragma myrmics task region inout(r) in(time_start) safe(time_start)
    barnes_myrmics_time_stop(r, time_start);

	// Free everything
//	#pragma myrmics task region inout(r) in(roots_r, body, lroot_p) safe(roots_r, body, lroot_p) in(num_bodies)
//	barnes_myrmics_free_all(r, roots_r, body, lroot_p, num_bodies);

}

#else



#define IGNORE_FSTEP 1         // ignore first step since the particles are not
                               // balanced yet and the results do not scale well 
                               // for a small number of steps.

#define COST_EVAL 1            // 0: use plain, non-cost based load balancing
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

#define INIT_CRD_MAX 200
#define CRD_MAX 1000
#define VEL_MAX 5
#define MASS_MAX 2

#define THETA 1.0F
#define TIME_STEP 0.125F

// How many Oct Tree levels we descend initially before changing region. Each
// time we descend, this number is halved. E.g. if set to 8, the first 8 Oct
// tree levels are in the first region, the next 4 are in the second region,
// the next 2 in the next etc. When it reaches 1, it stays there (we use one
// region for each additional oct tree level).
#define REG_HEIGHT 4

// Maximum number of cpus/regions 
#define MAX_REGIONS 256

// Stack size to operate recursively on the Oct tree. Increase this if the
// problem size is bigger and the assertion fails.
#define MAX_STACK_SIZE 1024

#define get_bit(num, pos)       ((num >> pos) & 1)
#define toggle_bit(num, pos)    (num ^ (1 << pos))


#include <myrmics.h>

// particle info
typedef struct {
	float	pos[3];
	float	vel[3];
	float	mass;
	unsigned int cost;
} particle_t;

// OctTree nodes
typedef struct OctTreeNodeStruct OctTreeNode;
struct OctTreeNodeStruct {

  // Pointers
//  OctTreeNode   *parent;        // OctTree parent node
  OctTreeNode   *children[8];   // Children pointers (see above for indexing)

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


// Structure used for bisection-based load balancing among processors
typedef struct {
	float         dim;            // Particle chosen dimension (x or y or z)
	unsigned long cost;           // Number of force evaluations for particle
} dim_exchange_t;

// Per-processor structure. Defines a bounding box and stores the root of the
// OctTree and the region array. All bodies in the local root belong to this
// processor.
typedef struct {

  float         origin[3];      // Bounding box origin location (x, y, z)
  float         size[3];        // Bounding box size (x, y, z)
  OctTreeNode   *oct_root;      // OctTree root
  int           num_levels;     // Number of tree levels

} LocalRoot;


// Body of a region
typedef struct {
	particle_t	*part;		// particles in this region
	int		*part_size;

	float	*bbox;
	int		*dim_choice;	// select dimention to bisect

	dim_exchange_t	*dims;	// merged dimensions used to find bisection
	int		*dims_size;

	float	*bisect;		// bisection point using dims
} body0_t;

typedef struct {
	rid_t	ltree_r;		// includes local tree
	LocalRoot	*lroot;		// local oct tree
} body1_t;


//////////////////////
//// LOAD BALANCE ////
//////////////////////
//void barnes_myrmics_bbox_init(particle_t *part, int *part_size, float *bbox, int step, int id) {
void barnes_myrmics_bbox_init(particle_t *part, int *part_size, float *bbox) {
	int j;

    bbox[0] = CRD_MAX;
    bbox[1] = -CRD_MAX;
    bbox[2] = CRD_MAX;
    bbox[3] = -CRD_MAX;
    bbox[4] = CRD_MAX;
    bbox[5] = -CRD_MAX;
    for (j = 0; j < *part_size; j++) {
      if (part[j].pos[0] < bbox[0])
        bbox[0] = part[j].pos[0];
      if (part[j].pos[0] > bbox[1])
        bbox[1] = part[j].pos[0];

      if (part[j].pos[1] < bbox[2])
        bbox[2] = part[j].pos[1];
      if (part[j].pos[1] > bbox[3])
        bbox[3] = part[j].pos[1];

      if (part[j].pos[2] < bbox[4])
        bbox[4] = part[j].pos[2];
      if (part[j].pos[2] > bbox[5])
        bbox[5] = part[j].pos[2];
    }

//if (step && (id==0)) printf("pos0=%f\n\r", part[0].pos[0]);

}

// Merge the two bbox
void barnes_myrmics_bbox_merge(float *bbox1, float *bbox2) {
	if (bbox1[0] < bbox2[0]) bbox2[0] = bbox1[0];
	if (bbox1[1] > bbox2[1]) bbox2[1] = bbox1[1];
	if (bbox1[2] < bbox2[2]) bbox2[2] = bbox1[2];
	if (bbox1[3] > bbox2[3]) bbox2[3] = bbox1[3];
	if (bbox1[4] < bbox2[4]) bbox2[4] = bbox1[4];
	if (bbox1[5] > bbox2[5]) bbox2[5] = bbox1[5];
}

// decide dimension to be bisected
void barnes_myrmics_dim_choice(float *bbox, int *dim_choice, int depth) {

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
//	printf("dim_choice=%d, d=%d\n\r", *dim_choice, depth);
}


void barnes_myrmics_part_init(unsigned int *seed, int part_size_init, particle_t *part, int *part_size) {
	int j;

	*part_size = part_size_init;

	// set random coordinates
	for (j = 0; j < part_size_init; j++) {
		part[j].pos[0] = (float) ((*seed = rand(*seed)) % (20*INIT_CRD_MAX)) / 10.0F - INIT_CRD_MAX;
		part[j].pos[1] = (float) ((*seed = rand(*seed)) % (20*INIT_CRD_MAX)) / 10.0F - INIT_CRD_MAX; 
		part[j].pos[2] = (float) ((*seed = rand(*seed)) % (20*INIT_CRD_MAX)) / 10.0F - INIT_CRD_MAX;
		part[j].vel[0] = (float) ((*seed = rand(*seed)) % (20*VEL_MAX)) / 10.0F - VEL_MAX;
		part[j].vel[1] = (float) ((*seed = rand(*seed)) % (20*VEL_MAX)) / 10.0F - VEL_MAX; 
		part[j].vel[2] = (float) ((*seed = rand(*seed)) % (20*VEL_MAX)) / 10.0F - VEL_MAX;
		part[j].mass = (float) (1.0F + (*seed = rand(*seed)) % (10*MASS_MAX)) / 10.0F;
		part[j].cost = 1;
	}

}

// initialize dims with selected dimension
void barnes_myrmics_dims_init(int *dim_choice, particle_t *part, int *part_size, dim_exchange_t *dims, int *dims_size) {
	int j;

    for (j = 0; j < *part_size; j++) {
      dims[j].dim = (*dim_choice == 0) ? part[j].pos[0] :
                       (*dim_choice == 1) ? part[j].pos[1] :
                                               part[j].pos[2];
      dims[j].cost = part[j].cost;
    }

	*dims_size = *part_size;
}

// merge dims
void barnes_myrmics_dims_merge(dim_exchange_t *dims1, int *dims1_size, dim_exchange_t *dims2, int *dims2_size) {
	int i;

	for (i = 0; i < *dims1_size; i++) {
		dims2[*dims2_size+i].dim = dims1[i].dim;
		dims2[*dims2_size+i].cost = dims1[i].cost;
	}

	*dims2_size += *dims1_size;
}


void qsort(dim_exchange_t *dims, int left, int right)
{
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

void barnes_myrmics_bisect(dim_exchange_t *dims, int *dims_size, float *bisect, int depth) {
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
//printf("bisect=%f, d=%d\n\r", *bisect, depth);
			return;
		}
	}
}


void barnes_myrmics_part_swaps(particle_t *part1, int *part_size1, particle_t *part2, int *part_size2, int *dim_choice, float *bisect, int id, int depth, int max_part_size) {
	int i;
	int part_pos1, part_pos2;

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
			sys_assert (part_pos2+1 <= max_part_size);
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
			sys_assert (part_pos1+1 <= max_part_size);
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
			sys_assert (part_pos2+1 <= max_part_size);
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

void barnes_myrmics_load_balance_top(body0_t ***body, int num_regions, int region_size, int max_part_size, int step) {
	int i, j;

//	printf("Load Balance Top start\n\r");
	// For each Orthogonal Recursive Bisection phase
	int depth;
	for (depth = int_log2(num_regions*region_size); depth > int_log2(region_size); depth--) {

		//initialize bbox according to particles' coordinates
		for (i = 0; i < num_regions; i++) {
		 for (j = 0; j < region_size; j++) {
			particle_t	*scoop_part = body[i][j]->part;
			int			*scoop_part_size = body[i][j]->part_size;
			float		*scoop_bbox = body[i][j]->bbox;
			int			id = i*region_size+j;
			#pragma myrmics task in(scoop_part, scoop_part_size) out(scoop_bbox) 
			barnes_myrmics_bbox_init(scoop_part, scoop_part_size, scoop_bbox);
		 }
		}

		//merge bbox between regions up to level depth
		for (i = 0; i <= depth; i++) {
			int start = ((2<<i)>>1)-1;
			int step = (2<<i)>>1;
			for (j = start; j < region_size*num_regions; ) {
				float	*scoop_bbox1 = body[j/region_size][j%region_size]->bbox;
				if (i < depth) {
					float	*scoop_bbox2 = body[(j+step)/region_size][(j+step)%region_size]->bbox;
					// merge when we are on intermidiate levels
					//printf("merge bbox %d %d\n\r", j, j+step);
					#pragma myrmics task in(scoop_bbox1) inout(scoop_bbox2)
					barnes_myrmics_bbox_merge(scoop_bbox1, scoop_bbox2);
					j += 2*step;
				} else {
					//printf("decide %d\n\r", j, j+step);
					int	*scoop_dim_choice = body[j/region_size][j%region_size]->dim_choice;
					//nodes on last level decide which bisection will be bisected
					#pragma myrmics task in(scoop_bbox1) out(scoop_dim_choice) in(depth)
					barnes_myrmics_dim_choice(scoop_bbox1, scoop_dim_choice, depth);
					j += step;
				}
			}
		}

		// initialize dimensions
		int src = ((2<<depth)>>1)-1;
		int step = (2<<depth)>>1;
		for (i = 0; i < num_regions; i++) {
		 for (j = 0; j < region_size; j++) {
			particle_t	*scoop_part = body[i][j]->part;
			int	*scoop_part_size = body[i][j]->part_size;
			if (i > src) src += step;
			int	*scoop_dim_choice = body[src/region_size][src%region_size]->dim_choice;
			dim_exchange_t	*scoop_dims = body[i][j]->dims;
			int	*scoop_dims_size = body[i][j]->dims_size;
			#pragma myrmics task in(scoop_dim_choice, scoop_part, scoop_part_size) out(scoop_dims, scoop_dims_size)
			barnes_myrmics_dims_init(scoop_dim_choice, scoop_part, scoop_part_size, scoop_dims, scoop_dims_size);
		 }
		}

		//merge dimensions between regions up to level depth and find bisection
		for (i = 0; i <= depth; i++) {
			int start = ((2<<i)>>1)-1;
			int step = (2<<i)>>1;
			for (j = start; j < num_regions*region_size; ) {
				dim_exchange_t	*scoop_dims1 = body[j/region_size][j%region_size]->dims;
				int		*scoop_dims1_size = body[j/region_size][j%region_size]->dims_size;
				if (i < depth) {
					dim_exchange_t	*scoop_dims2 = body[(j+step)/region_size][(j+step)%region_size]->dims;
					int		*scoop_dims2_size = body[(j+step)/region_size][(j+step)%region_size]->dims_size;
					// merge when we are on intermidiate levels
					//printf("merge dims %d %d\n\r", j, j+step);
					#pragma myrmics task in(scoop_dims1, scoop_dims1_size) inout(scoop_dims2, scoop_dims2_size)
					barnes_myrmics_dims_merge(scoop_dims1, scoop_dims1_size, scoop_dims2, scoop_dims2_size);

					j += 2*step;
				} else {
					//printf("bisect %d\n\r", j, j+step);
					float	*scoop_bisect = body[j/region_size][j%region_size]->bisect;
					//nodes on last level decide which bisection will be bisected
					#pragma myrmics task in(scoop_dims1, scoop_dims1_size) out(scoop_bisect) in(depth)
					barnes_myrmics_bisect(scoop_dims1, scoop_dims1_size, scoop_bisect, depth);

					j += step;
				}
			}
		}

		// transfer particles between regions
		int src = ((2<<depth)>>1)-1;
		int step = (2<<depth)>>1;
		for (i = 0; i < num_regions*region_size; i++) {
			if (!get_bit(i, depth-1)) {
				particle_t	*scoop_part1 = body[i/region_size][i%region_size]->part;
				int			*scoop_part_size1 = body[i/region_size][i%region_size]->part_size;
				int			peer = toggle_bit(i, depth-1);
				particle_t	*scoop_part2 = body[peer/region_size][peer%region_size]->part;
				int			*scoop_part_size2 = body[peer/region_size][peer%region_size]->part_size;
				if (i > src) src += step;
				int			*scoop_dim_choice = body[src/region_size][src%region_size]->dim_choice;
				float		*scoop_bisect = body[src/region_size][src%region_size]->bisect;
				//printf("bisect %d %d\n\r", i, src);
				//printf("swap particles between %d to %d\n\r", i, toggle_bit(i, depth-1));
				#pragma myrmics task inout(scoop_part1, scoop_part_size1, scoop_part2, scoop_part_size2) in(scoop_dim_choice, scoop_bisect, i, depth, max_part_size)
				barnes_myrmics_part_swaps(scoop_part1, scoop_part_size1, scoop_part2, scoop_part_size2, scoop_dim_choice, scoop_bisect, i, depth, max_part_size);
			}
		}
	}
}


void barnes_myrmics_load_balance(rid_t body_region, body0_t **body, int region_size, int max_part_size) {
	int i, j;

//	printf("Load Balance start\n\r");
	// For each Orthogonal Recursive Bisection phase
	int depth;
	for (depth = int_log2(region_size); depth > 0; depth--) {

		//initialize bbox according to particles' coordinates
		for (i = 0; i < region_size; i++) {
			particle_t	*scoop_part = body[i]->part;
			int			*scoop_part_size = body[i]->part_size;
			float		*scoop_bbox = body[i]->bbox;
			#pragma myrmics task in(scoop_part, scoop_part_size) out(scoop_bbox)
			barnes_myrmics_bbox_init(scoop_part, scoop_part_size, scoop_bbox);
		}

		//merge bbox between regions up to level depth
		for (i = 0; i <= depth; i++) {
			int start = ((2<<i)>>1)-1;
			int step = (2<<i)>>1;
			for (j = start; j < region_size; ) {
				float	*scoop_bbox1 = body[j]->bbox;
				if (i < depth) {
					float	*scoop_bbox2 = body[j+step]->bbox;
					// merge when we are on intermidiate levels
					//printf("merge bbox %d %d\n\r", j, j+step);
					#pragma myrmics task in(scoop_bbox1) inout(scoop_bbox2)
					barnes_myrmics_bbox_merge(scoop_bbox1, scoop_bbox2);
					j += 2*step;
				} else {
					//printf("decide %d\n\r", j, j+step);
					int	*scoop_dim_choice = body[j]->dim_choice;
					//nodes on last level decide which bisection will be bisected
					#pragma myrmics task in(scoop_bbox1) out(scoop_dim_choice) in(depth)
					barnes_myrmics_dim_choice(scoop_bbox1, scoop_dim_choice, depth);
					j += step;
				}
			}
		}

		// initialize dimensions
		int src = ((2<<depth)>>1)-1;
		int step = (2<<depth)>>1;
		for (i = 0; i < region_size; i++) {
			particle_t	*scoop_part = body[i]->part;
			int	*scoop_part_size = body[i]->part_size;
			if (i > src) src += step;
			int	*scoop_dim_choice = body[src]->dim_choice;
			dim_exchange_t	*scoop_dims = body[i]->dims;
			int	*scoop_dims_size = body[i]->dims_size;
			#pragma myrmics task in(scoop_dim_choice, scoop_part, scoop_part_size) out(scoop_dims, scoop_dims_size)
			barnes_myrmics_dims_init(scoop_dim_choice, scoop_part, scoop_part_size, scoop_dims, scoop_dims_size);
		}

		//merge dimensions between regions up to level depth and find bisection
		for (i = 0; i <= depth; i++) {
			int start = ((2<<i)>>1)-1;
			int step = (2<<i)>>1;
			for (j = start; j < region_size; ) {
				dim_exchange_t	*scoop_dims1 = body[j]->dims;
				int		*scoop_dims1_size = body[j]->dims_size;
				if (i < depth) {
					dim_exchange_t	*scoop_dims2 = body[j+step]->dims;
					int		*scoop_dims2_size = body[j+step]->dims_size;
					// merge when we are on intermidiate levels
					//printf("merge dims %d %d\n\r", j, j+step);
					#pragma myrmics task in(scoop_dims1, scoop_dims1_size) inout(scoop_dims2, scoop_dims2_size)
					barnes_myrmics_dims_merge(scoop_dims1, scoop_dims1_size, scoop_dims2, scoop_dims2_size);

					j += 2*step;
				} else {
					//printf("bisect %d\n\r", j, j+step);
					float	*scoop_bisect = body[j]->bisect;
					//nodes on last level decide which bisection will be bisected
					#pragma myrmics task in(scoop_dims1, scoop_dims1_size) out(scoop_bisect) in(depth)
					barnes_myrmics_bisect(scoop_dims1, scoop_dims1_size, scoop_bisect, depth);

					j += step;
				}
			}
		}

		// transfer particles between regions
		int src = ((2<<depth)>>1)-1;
		int step = (2<<depth)>>1;
		for (i = 0; i < region_size; i++) {
			if (!get_bit(i, depth-1)) {
				particle_t	*scoop_part1 = body[i]->part;
				int			*scoop_part_size1 = body[i]->part_size;
				particle_t	*scoop_part2 = body[toggle_bit(i, depth-1)]->part;
				int			*scoop_part_size2 = body[toggle_bit(i, depth-1)]->part_size;
				if (i > src) src += step;
				int			*scoop_dim_choice = body[src]->dim_choice;
				float		*scoop_bisect = body[src]->bisect;
				//printf("bisect %d %d\n\r", i, src);
				//printf("swap particles between %d to %d\n\r", i, toggle_bit(i, depth-1));
				#pragma myrmics task inout(scoop_part1, scoop_part_size1, scoop_part2, scoop_part_size2) in(scoop_dim_choice, scoop_bisect, i, depth, max_part_size)
				barnes_myrmics_part_swaps(scoop_part1, scoop_part_size1, scoop_part2, scoop_part_size2, scoop_dim_choice, scoop_bisect, i, depth, max_part_size);
			}
		}
	}
}

//////////////////////
////    OCT TREE  ////
//////////////////////

int level_to_region(int level) {

  int reg = 0;
  int reg_rem_levels = REG_HEIGHT - 1;
  int reg_limit = REG_HEIGHT;
  int i;

  for (i = 0; i < level; i++) {
    if (reg_rem_levels) {
      reg_rem_levels--;
    }
    else {
      if (reg_limit > 1) {
        reg_limit /= 2;
      }
      reg_rem_levels = reg_limit - 1;
      reg++;
    }
  }

  return reg;
}

// ===========================================================================
// ===========================================================================
void create_new_child(LocalRoot *lroot, rid_t ltree_r, OctTreeNode *cur_oct, int level,
                      int child, OctTreeNode **ret_new_oct, int *ret_level) {

  OctTreeNode   *new_oct;
  int           new_level;


  // Change levels
  new_level = level + 1;

  // Create the new, empty Oct node. Oct tree nodes are always part of 
  // a region which is dependent on the oct tree level.
  sys_assert(new_oct = sys_alloc(sizeof(OctTreeNode), ltree_r));
//  new_oct->parent = cur_oct;
  cur_oct->children[child] = new_oct;
  memset(new_oct->children, 0, 8 * sizeof(OctTreeNode *));
  new_oct->is_leaf = 0;
  new_oct->visited = 0;

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
}


// ===========================================================================
// ===========================================================================
void barnes_myrmics_create_oct_tree(int id, particle_t *part, int *part_size, rid_t ltree_r, LocalRoot *lroot) {
	rid_t         tmp;
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

//printf("oct tree %d\n\r", id);

  // Sanity checks
  sys_assert(*part_size > 0);

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
  sys_assert((min[0] < CRD_MAX) && (min[0] > -CRD_MAX));
  sys_assert((min[1] < CRD_MAX) && (min[1] > -CRD_MAX));
  sys_assert((min[2] < CRD_MAX) && (min[2] > -CRD_MAX));
  sys_assert((max[0] < CRD_MAX) && (max[0] > -CRD_MAX));
  sys_assert((max[1] < CRD_MAX) && (max[1] > -CRD_MAX));
  sys_assert((max[2] < CRD_MAX) && (max[2] > -CRD_MAX));

  // Set bounding box according to min and max values
  lroot->origin[0] = min[0];
  lroot->origin[1] = min[1];
  lroot->origin[2] = min[2];
  lroot->size[0]   = max[0] - min[0];
  lroot->size[1]   = max[1] - min[1];
  lroot->size[2]   = max[2] - min[2];

  // Allocate oct tree root, using the first body
  sys_assert(cur_oct = sys_alloc(sizeof(OctTreeNode), ltree_r));
  lroot->oct_root = cur_oct;
//  cur_oct->parent = NULL;
  memset(cur_oct->children, 0, 8 * sizeof(OctTreeNode *));
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

  sys_assert(cur_oct->pos[0] >= lroot->origin[0]);
  sys_assert(cur_oct->pos[0] <= lroot->origin[0] +
                                lroot->size[0]+ROUND_ERROR2);
  sys_assert(cur_oct->pos[1] >= lroot->origin[1]);
  sys_assert(cur_oct->pos[1] <= lroot->origin[1] +
                                lroot->size[1]+ROUND_ERROR2);
  sys_assert(cur_oct->pos[2] >= lroot->origin[2]);
  sys_assert(cur_oct->pos[2] <= lroot->origin[2] +
                                lroot->size[2]+ROUND_ERROR2);

  // Deal with all the rest
  for (i = 1; i < *part_size; i++) {

    // Start from the top
    level = 0;
    cur_oct = lroot->oct_root;
    cur_size[0] = lroot->size[0];
    cur_size[1] = lroot->size[1];
    cur_size[2] = lroot->size[2];
    cur_centre[0] = lroot->origin[0] + cur_size[0] / 2.0F;
    cur_centre[1] = lroot->origin[1] + cur_size[1] / 2.0F;
    cur_centre[2] = lroot->origin[2] + cur_size[2] / 2.0F;
    reloc = NULL;

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
      if (cur_oct->children[child]) {
        cur_oct = cur_oct->children[child];

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
        sys_assert(!reloc);
        reloc = cur_oct;
        cur_oct->is_leaf = 0;
      }

      // Create a new child, possibly creating a new region as well
      create_new_child(lroot, ltree_r, cur_oct, level, child, &new_oct, &new_level);

      // Do we have an old body to relocate?
      if (reloc) {

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
          create_new_child(lroot, ltree_r, cur_oct, level, child_reloc, &new2_oct, NULL);

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

          sys_assert(new2_oct->pos[0] >= lroot->origin[0]);
          sys_assert(new2_oct->pos[0] <= lroot->origin[0] +
                                        lroot->size[0]+ROUND_ERROR2);
          sys_assert(new2_oct->pos[1] >= lroot->origin[1]);
          sys_assert(new2_oct->pos[1] <= lroot->origin[1] +
                                        lroot->size[1]+ROUND_ERROR2);
          sys_assert(new2_oct->pos[2] >= lroot->origin[2]);
          sys_assert(new2_oct->pos[2] <= lroot->origin[2] +
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

        sys_assert(new_oct->pos[0] >= lroot->origin[0]);
        sys_assert(new_oct->pos[0] <= lroot->origin[0] +
                                      lroot->size[0]+ROUND_ERROR2);
        sys_assert(new_oct->pos[1] >= lroot->origin[1]);
        sys_assert(new_oct->pos[1] <= lroot->origin[1] +
                                      lroot->size[1]+ROUND_ERROR2);
        sys_assert(new_oct->pos[2] >= lroot->origin[2]);
        sys_assert(new_oct->pos[2] <= lroot->origin[2] +
                                      lroot->size[2]+ROUND_ERROR2);

        // We're done
        break;
      }

      // Otherwise, both the new and the old body need to go further down 
      // into the tree because they are too near to each other for this
      // granularity. Continue descending.

      // Leave the new node as non-leaf
      cur_oct = new_oct;
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

  barnes_myrmics_update_mass_centers(id, ltree_r, lroot);
}



//////////////////////
////  MASS CENTER ////
//////////////////////

void barnes_myrmics_update_mass_centers(int id, rid_t ltree_r, LocalRoot *lroot) {

  OctTreeNode   *stack[MAX_STACK_SIZE];
  int           stack_size;
  int           i;
  OctTreeNode   *node;


  // Initialize the stack with the oct tree root
  stack[0] = lroot->oct_root;
  stack_size = 1;

  // Iterate on the stack
  while (stack_size) {

    // Peek next node
    sys_assert(stack_size > 0);
    node = stack[stack_size - 1];

    // Does it have children?
    if (!node->is_leaf) {

      // If we haven't visited this in the past...
      if (!node->visited) {

        // ... leave it on the stack and add all its children, to be updated
        // before this node.
        for (i = 0; i < 8; i++) {
          if (node->children[i]) {
            stack[stack_size++] = node->children[i];
            sys_assert(stack_size < MAX_STACK_SIZE);
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
          if (node->children[i]) {
            node->mass   += node->children[i]->mass;
          }
        }
        for (i = 0; i < 8; i++) {
          if (node->children[i]) {
            node->pos[0] += node->children[i]->pos[0] *
                            (node->children[i]->mass / node->mass);
            node->pos[1] += node->children[i]->pos[1] *
                            (node->children[i]->mass / node->mass);
            node->pos[2] += node->children[i]->pos[2] *
                            (node->children[i]->mass / node->mass);
          }
        }

        // Sanity check 1: that's a dangerous spot for inf and nan
        sys_assert((node->pos[0] < CRD_MAX) && (node->pos[0] > -CRD_MAX));
        sys_assert((node->pos[1] < CRD_MAX) && (node->pos[1] > -CRD_MAX));
        sys_assert((node->pos[2] < CRD_MAX) && (node->pos[2] > -CRD_MAX));

        // Sanity check 2: it should be part of our bounding box
        sys_assert(node->pos[0] >= lroot->origin[0]);
        sys_assert(node->pos[0] <= lroot->origin[0] + lroot->size[0]+ROUND_ERROR);
        sys_assert(node->pos[1] >= lroot->origin[1]);
        sys_assert(node->pos[1] <= lroot->origin[1] + lroot->size[1]+ROUND_ERROR);
        sys_assert(node->pos[2] >= lroot->origin[2]);
        sys_assert(node->pos[2] <= lroot->origin[2] + lroot->size[2]+ROUND_ERROR);

        // Pop it
        stack_size--;
      }
    }
    else {
      // Make sure it's in our bounding box, and adjust for any rounding
      // errors. We can encounter these even in leaf nodes, due to the
      // original min/max which is converted to origin/size.
      for (i = 0; i < 3; i++) {
        sys_assert(node->pos[i] >= lroot->origin[i]);
        if (node->pos[i] > lroot->origin[i] + lroot->size[i]) {
          sys_assert(node->pos[i] <=
                 lroot->origin[i] + lroot->size[i] * (1.0F + ROUND_ERROR));
          lroot->size[i] *= (1.0F + ROUND_ERROR);
        }
      }

      // We don't care if it's a leaf. Just pop it.
      stack_size--;
    }
  }
}


void barnes_myrmics_update_particles2(rid_t r, LocalRoot **lroot_p, int *part_size, particle_t *part, int id, int num_bodies) {

}

//void barnes_myrmics_update_particles(rid_t roots_r0, rid_t roots_r1, rid_t roots_r2, rid_t roots_r3, rid_t trees_r0, rid_t trees_r1, rid_t trees_r2, rid_t trees_r3, LocalRoot **lroot_p, int *part_size, particle_t *part, int id, int num_bodies) {
void barnes_myrmics_update_particles(rid_t r, LocalRoot **lroot_p, int *part_size, particle_t *part, int id, int num_bodies) {

  OctTreeNode   *stack1[MAX_STACK_SIZE];
  int           stack1_level[MAX_STACK_SIZE];
  int           stack1_size;
  OctTreeNode   *node1;
  particle_t    *part_n;
  int           level1;
  OctTreeNode   *stack2[MAX_STACK_SIZE];
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

	
//printf("update particles\n\r");
  // Initialize stack1 with the local oct tree root. This will iterate on all
  // the local bodies that we're responsible to update forces and velocities.
  stack1[0] = lroot_p[id]->oct_root;
  stack1_level[0] = 0;
  stack1_size = 1;

  my_root = lroot_p[id];

  // index in particles
  part_i = 0;

// printf("Update parts %d\n\r", id);
  // Iterate on stack1
  while (stack1_size) {

    // Pop next node
    sys_assert(stack1_size > 0);
    node1 = stack1[--stack1_size];

    level1 = stack1_level[stack1_size];

    // Sanity check: it should be part of our bounding box
    sys_assert(node1->pos[0] >= my_root->origin[0]);
    sys_assert(node1->pos[0] <= my_root->origin[0] +
                            my_root->size[0]+ROUND_ERROR2);
    sys_assert(node1->pos[1] >= my_root->origin[1]);
    sys_assert(node1->pos[1] <= my_root->origin[1] +
                            my_root->size[1]+ROUND_ERROR2);
    sys_assert(node1->pos[2] >= my_root->origin[2]);
    sys_assert(node1->pos[2] <= my_root->origin[2] +
                            my_root->size[2]+ROUND_ERROR2);

    // Does it have children?
    if (!node1->is_leaf) {
      // Push them on the stack
      for (i = 0; i < 8; i++) {
        if (node1->children[i]) {
          stack1_level[stack1_size] = level1 + 1;
          stack1[stack1_size++] = node1->children[i];
          sys_assert(stack1_size < MAX_STACK_SIZE);
          // Check that the children level is sane
          sys_assert(level1 + 1 < my_root->num_levels);
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
        sys_assert(stack2_size < MAX_STACK_SIZE);

        sys_assert(lroot_p[i]->oct_root->pos[0] >= lroot_p[i]->origin[0]);

      }

      // Iterate on stack2
      while (stack2_size) {

        // Pop next node
        sys_assert(stack2_size > 0);
        node2 = stack2[--stack2_size];
        level2         = stack2_level[2*stack2_size];
        root_id2       = stack2_level[2*stack2_size + 1];
        box2_origin[0] = stack2_box_origin[3*stack2_size];
        box2_origin[1] = stack2_box_origin[3*stack2_size + 1];
        box2_origin[2] = stack2_box_origin[3*stack2_size + 2];
        box2_size[0]   = stack2_box_size[3*stack2_size];
        box2_size[1]   = stack2_box_size[3*stack2_size + 1];
        box2_size[2]   = stack2_box_size[3*stack2_size + 2];

        // Sanity check: it must be inside root_id's bounding box
        sys_assert(node2->pos[0] >= lroot_p[root_id2]->origin[0]);
        sys_assert(node2->pos[0] <= lroot_p[root_id2]->origin[0] +
                                lroot_p[root_id2]->size[0]+ROUND_ERROR2);
        sys_assert(node2->pos[1] >= lroot_p[root_id2]->origin[1]);
        sys_assert(node2->pos[1] <= lroot_p[root_id2]->origin[1] +
                                lroot_p[root_id2]->size[1]+ROUND_ERROR2);
        sys_assert(node2->pos[2] >= lroot_p[root_id2]->origin[2]);
        sys_assert(node2->pos[2] <= lroot_p[root_id2]->origin[2] +
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
        dist = sqrt(dist_squared);

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
            if (node2->children[i]) {

              // Push the child
              stack2[stack2_size] = node2->children[i];
              stack2_level[2 * stack2_size] = level2 + 1; // next level
              stack2_level[2 * stack2_size + 1] = root_id2; // same root

              // Check that the child level is sane
              if (root_id2 == id) {
                // Check local level
                sys_assert(level2 + 1 < my_root->num_levels);
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
              sys_assert(stack2_size < MAX_STACK_SIZE);
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
          sys_assert((dist < FLT_MAX) && (dist > 0));
          sys_assert((ft   < FLT_MAX) && (ft   > -FLT_MAX));
          sys_assert((f[0] < FLT_MAX) && (f[0] > -FLT_MAX));
          sys_assert((f[1] < FLT_MAX) && (f[1] > -FLT_MAX));
          sys_assert((f[2] < FLT_MAX) && (f[2] > -FLT_MAX));
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
      sys_assert((part_n->vel[0] < FLT_MAX) && (part_n->vel[0] > -FLT_MAX));
      sys_assert((part_n->vel[1] < FLT_MAX) && (part_n->vel[1] > -FLT_MAX));
      sys_assert((part_n->vel[2] < FLT_MAX) && (part_n->vel[2] > -FLT_MAX));
    }
  }

  sys_assert(part_i == *part_size);
//if (id == 0) {
//	for (i = 0; i < part_i; i++) {
//		printf("pos0=%f\n\r", part[i].pos[0]);
//	}
//}

}


void barnes_myrmics_create_oct_tree_region(rid_t body0_region, rid_t body1_region, body0_t **body0, body1_t **body1, int region_size) {
	int i;

//    printf("create tree region\n\r");
	for (i = 0; i < region_size; i++) {
		rid_t		scoop_ltree_r = body1[i]->ltree_r;
		particle_t	*scoop_part = body0[i]->part;
		int			*scoop_part_size = body0[i]->part_size;
		LocalRoot	*scoop_lroot = body1[i]->lroot;

		// NEW OCT TREE [Computation]
	    //
	    // Build a new Oct tree from our local bodies
		#pragma myrmics task in(i, scoop_part, scoop_part_size) region out(scoop_ltree_r) out(scoop_lroot) 
		barnes_myrmics_create_oct_tree(i, scoop_part, scoop_part_size, scoop_ltree_r, scoop_lroot);

	    // UPDATE MASS CENTERS [Computation]
	    //
	    // Walk the Oct tree and find the centers of mass of all non-leaf nodes
//		#pragma myrmics task in(i) region inout(scoop_ltree_r) in(scoop_lroot) 
//		barnes_myrmics_update_mass_centers(i, scoop_ltree_r, scoop_lroot);
	}
//printf("Done\n\r");
}

void barnes_myrmics_phase(rid_t r, rid_t r0, rid_t r1, rid_t *body0_regions, rid_t *body1_regions, rid_t *roots_r, rid_t *trees_r, body0_t ***body0, body1_t ***body1, LocalRoot **lroot_p, unsigned int *time_start, int step, int region_size, int max_part_size, int num_regions) {
	int i, j;

	particle_t	*up_part[MAX_REGIONS];
	int			*up_part_size[MAX_REGIONS];

	printf("Step %d\n\r", step);
	if (num_regions * region_size > MAX_REGIONS) {
		printf("Error: Can't have more than %d regions\n\r", MAX_REGIONS);
		return;
	}

#if IGNORE_FSTEP
	if (step == 1) *time_start = sys_free_timer_get_ticks();
#endif

	// reallocate trees
	if (0) {
	 printf("Realloc trees\n\r");
	 for (i = 0; i < num_regions; i++) {
		sys_rfree(trees_r[i]);
		trees_r[i] = sys_ralloc(body1_regions[i], 0); 
		for (j = 0; j < region_size; j++) {
			body1[i][j]->ltree_r = sys_ralloc(trees_r[i], 0);
		}
	 }
	}

	// Remember pointers of part and part_size before you call load_balance in order to
	// use them in update_particles
	for (i = 0; i < num_regions; i++) {
	 for (j = 0; j < region_size; j++) {
		up_part[i*region_size+j] = body0[i][j]->part;
		up_part_size[i*region_size+j] = body0[i][j]->part_size;
	 }
	}

    // LOAD BALANCING [Computation & Communication]
    //
    // Make cores agree on how to split the whole space into multiple
    // orthogonal recursive bisections; have them exchange particles that
    // comply with the newly agreed bisections.
	//
	// Load balance is slow. Perform once every 4 steps. Effectiveness 
	// of load balance depends on the VEL_MAX, MASS_MAX: if they are small
	// then load balance does not have to be called so often.
	if (!(step & 3)) {
	barnes_myrmics_load_balance_top(body0, num_regions, region_size, max_part_size, step);


	for (i = 0; i < num_regions; i++) {
		rid_t	scoop_body_region = body0_regions[i];
		body0_t	**scoop_body = body0[i]; 
		#pragma myrmics task region inout(scoop_body_region) in(scoop_body) safe(scoop_body) in(region_size, max_part_size)
		barnes_myrmics_load_balance(scoop_body_region, scoop_body, region_size, max_part_size);
	}
	}


	for (i = 0; i < num_regions; i++) {
		rid_t	scoop_body0_region = body0_regions[i];
		rid_t	scoop_body1_region = body1_regions[i];
		body0_t	**scoop_body0 = body0[i]; 
		body1_t	**scoop_body1 = body1[i]; 
		#pragma myrmics task region in(scoop_body0_region) region out(scoop_body1_region) in(scoop_body0, scoop_body1) safe(scoop_body0, scoop_body1) in(region_size)
		barnes_myrmics_create_oct_tree_region(scoop_body0_region, scoop_body1_region, scoop_body0, scoop_body1, region_size);
	}

	// added delay in order to avoid the deadlock in Myrmics!
//	volatile unsigned int t;
//	for (t = 0; t < 300000; t++);

    // CALCULATE NEW PARTICLES
    //
	for (i = 0; i < num_regions; i++) {
	 for (j = 0; j < region_size; j++) {
		particle_t	*scoop_part = up_part[i*region_size+j];
		int			*scoop_part_size = up_part_size[i*region_size+j];
		int			num_bodies = num_regions * region_size;
		int			id = i * region_size + j;
		#pragma myrmics task region in(r1) in(lroot_p, scoop_part_size) out(scoop_part) in(id, num_bodies)
		barnes_myrmics_update_particles(r1, lroot_p, scoop_part_size, scoop_part, id, num_bodies);
	 }
	}


}

// ===========================================================================
// ===========================================================================
void barnes_myrmics_time_start(rid_t unused, unsigned int *time_start) {
	*time_start = sys_free_timer_get_ticks();
}

void barnes_myrmics_time_stop(rid_t unused, unsigned int *time_start) {
  unsigned int          time_stop;
  unsigned int          time;

  time_stop = sys_free_timer_get_ticks();
  if (time_stop > *time_start) {
    time = time_stop - *time_start;
  }
  else {
    time = 0xFFFFFFFF - (*time_start - time_stop);
  }

  printf("Time: %10u cycles (%6u msec)\r\n", time, time / 10000);
}



void barnes_myrmics(int num_bodies, int num_particles, int num_regions, int steps) {
	rid_t	r;
	rid_t	r0;
	rid_t	r1;
	rid_t	*roots_r;
	rid_t	*trees_r;
	LocalRoot **lroot_p;
	body0_t	***body0;
	body1_t	***body1;
	int	particles_per_region;
	int i, j;
	unsigned int	*seed;
	unsigned int	*time_start;
	int				max_part_size = (num_particles / num_bodies)*3;
	rid_t			*body0_regions;
	rid_t			*body1_regions;
	int				region_size;

	// Create all-holding body
	r = sys_ralloc(0, 99); // highest level
	r0 = sys_ralloc(r, 99); // highest level
	r1 = sys_ralloc(r, 99); // highest level

	// initialize variables
	particles_per_region = num_particles/num_bodies;
	sys_assert(particles_per_region <= max_part_size);
	seed = sys_alloc(sizeof(unsigned int), r);
	*seed = 42;
	time_start = sys_alloc(sizeof(unsigned int), r);

	// Sanity checks
	if (num_bodies & (num_bodies - 1)) {
		printf("num_bodies must be a power of 2\r\n");
		return;
	}
	if (num_regions & (num_regions - 1)) {
		printf("num_regions must be a power of 2\r\n");
		return;
	}

	lroot_p = sys_alloc(num_bodies * sizeof(LocalRoot *), r); // points to lroots

	// Create a region for each scheduler
	body0_regions = sys_alloc(num_regions * sizeof(rid_t), r0);
	body1_regions = sys_alloc(num_regions * sizeof(rid_t), r1);
	roots_r = sys_alloc(num_regions * sizeof(rid_t), r);
	trees_r = sys_alloc(num_regions * sizeof(rid_t), r);
	for (i = 0; i < num_regions; i++) {
		body0_regions[i] = sys_ralloc(r0, 0); // lowest level
		body1_regions[i] = sys_ralloc(r1, 0); // lowest level
		roots_r[i] = sys_ralloc(body1_regions[i], 0); // this region includes all lroot's
		trees_r[i] = sys_ralloc(body1_regions[i], 0); // this region includes all oct trees
	}
	
	// Create regions with bodies
	region_size = num_bodies / num_regions;
	body0 = sys_alloc(num_regions * sizeof(body0_t **), r0);
	body1 = sys_alloc(num_regions * sizeof(body1_t **), r1);
	for (i = 0; i < num_regions; i++) {
	 body0[i] = sys_alloc(region_size * sizeof(body0_t *), body0_regions[i]);
	 body1[i] = sys_alloc(region_size * sizeof(body0_t *), body1_regions[i]);
	 for (j = 0; j < region_size; j++) {
		// put body in region
		body0[i][j] = sys_alloc(sizeof(body0_t), body0_regions[i]);
		body1[i][j] = sys_alloc(sizeof(body1_t), body1_regions[i]);
		// allocate particles 
		body0[i][j]->part = sys_alloc(max_part_size * sizeof(particle_t), body0_regions[i]);
		body0[i][j]->part_size = sys_alloc(sizeof(int), body0_regions[i]);
		particle_t	*scoop_part = body0[i][j]->part;
		int	*scoop_part_size = body0[i][j]->part_size;

		// initialize particles
		#pragma myrmics task inout(seed) in(particles_per_region) out(scoop_part, scoop_part_size)
		barnes_myrmics_part_init(seed, particles_per_region, scoop_part, scoop_part_size);
		
		// allocate bbox
		body0[i][j]->bbox = sys_alloc(6 * sizeof(float), body0_regions[i]);

		// allocate dim_choice, dims, bisect, etc.
		body0[i][j]->dim_choice = sys_alloc(sizeof(int), body0_regions[i]);
		body0[i][j]->dims = sys_alloc(particles_per_region * num_bodies * sizeof(dim_exchange_t), body0_regions[i]);
		body0[i][j]->dims_size = sys_alloc(sizeof(int), body0_regions[i]);
		body0[i][j]->bisect = sys_alloc(sizeof(float), body0_regions[i]);
//		body1[i][j]->lroot = sys_alloc(sizeof(LocalRoot), roots_r[i]);
		body1[i][j]->lroot = sys_alloc(sizeof(LocalRoot), body1_regions[i]);
		lroot_p[i * region_size + j] = body1[i][j]->lroot;

//		body1[i][j]->ltree_r = sys_ralloc(trees_r[i], 0);
		body1[i][j]->ltree_r = sys_ralloc(body1_regions[i], 0);
	 }
	}

	// Print we're starting
    printf("Barnes Hut of %d particles (%d per region) splitted into %d regions; running %d steps.\r\n",
           num_particles, particles_per_region, num_bodies, steps);

    // Start time
    #pragma myrmics task region inout(r) in(time_start) safe(time_start)
    barnes_myrmics_time_start(r, time_start);

	// Run main loop sequentially
	int step;
//	for (step = 0; step < steps; step++) {
	for (step = 0; step < steps; step++) {
		#pragma myrmics task region inout(r) in(r0, r1, body0_regions, body1_regions, roots_r, trees_r, body0, body1, lroot_p, time_start) safe(r0, r1, body0_regions, body1_regions, roots_r, trees_r, body0, body1, lroot_p, time_start) in(step, region_size, max_part_size, num_regions)
		barnes_myrmics_phase(r, r0, r1,  body0_regions, body1_regions, roots_r, trees_r, body0, body1, lroot_p, time_start, step, region_size, max_part_size, num_regions);

	}

    // Stop time
    #pragma myrmics task region inout(r) in(time_start) safe(time_start)
    barnes_myrmics_time_stop(r, time_start);

	// Free everything
//	#pragma myrmics task region inout(r) in(roots_r, body, lroot_p) safe(roots_r, body, lroot_p) in(num_bodies)
//	barnes_myrmics_free_all(r, roots_r, body, lroot_p, num_bodies);

}

#endif
