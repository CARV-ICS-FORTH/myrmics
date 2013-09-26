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
// Abstract      : Raytracing filter, MPI version
// =============================[ CVS Variables ]=============================
// 
// File name     : $RCSfile: cray_mpi.c,v $
// CVS revision  : $Revision: 1.1 $
// Last modified : $Date: 2013/01/22 13:40:54 $
// Last author   : $Author: lyberis-spree $
// 
// ===========================================================================

#include <arch.h>
#include <kernel_toolset.h>
#include <fmpi.h>


struct vec3 {
  float x, y, z;
};

struct ray {
  struct vec3 orig, dir;
};

struct material {
  struct vec3 col;        /* color */
  float spow;             /* specular power */
  float refl;             /* reflection intensity */
};

struct sphere {
  struct vec3 pos;
  float rad;
  struct material mat;
};

struct spoint {
  struct vec3 pos, normal, vref;  /* position, normal and view reflection */
  float dist;        /* parametric distance of intersection along the ray */
};

struct camera {
  struct vec3 pos, targ;
  float fov;
};

struct vec3 cray_mpi_shade(struct sphere *obj, struct spoint *sp, int depth,
                           struct sphere *objects, int onum, 
                           struct vec3 *lights, int lnum);

#define MAX_OBJECTS     128             /* maximum number of objects */
#define MAX_LIGHTS      16              /* maximum number of lights */
#define RAY_MAG         1000.0F         /* trace rays of this magnitude */
#define MAX_RAY_DEPTH   5               /* raytrace recursion limit */
#define FOV             0.78539816F     /* field of view in rads (pi/4) */
#define HALF_FOV        (FOV * 0.5F)
#define ERR_MARGIN      1e-3F           /* an arbitrary error margin to avoid 
                                           surface acne */
#define RAYS_PER_PIXEL  1

/* bit-shift ammount for packing each color into a 32bit uint */
#define RSHIFT  16
#define BSHIFT  0
#define GSHIFT  8       /* this is the same in both byte orders */

/* some helpful macros... */
#define SQ(x)           ((x) * (x))
#define MAX(a, b)       ((a) > (b) ? (a) : (b))
#define MIN(a, b)       ((a) < (b) ? (a) : (b))
#define DOT(a, b)       ((a).x * (b).x + (a).y * (b).y + (a).z * (b).z)
#define NORMALIZE(a)  do {\
        float len = ar_float_sqrt(DOT(a, a));\
        (a).x /= len; (a).y /= len; (a).z /= len;\
} while(0);

#define NRAN    1024
#define MASK    (NRAN - 1)


typedef struct {
  struct sphere objects[MAX_OBJECTS];
  int           onum;
  struct vec3   lights[MAX_LIGHTS];
  int           lnum;
  struct camera cam;
  struct vec3   urand[NRAN];
  int           irand[NRAN];
} global_t;



// ===========================================================================
// Calculates reflection vector 
// ===========================================================================
struct vec3 cray_mpi_reflect(struct vec3 v, struct vec3 n) {
  struct vec3 res;
  float dot = v.x * n.x + v.y * n.y + v.z * n.z;
  res.x = -(2.0F * dot * n.x - v.x);
  res.y = -(2.0F * dot * n.y - v.y);
  res.z = -(2.0F * dot * n.z - v.z);
  return res;
}


// ===========================================================================
// ===========================================================================
struct vec3 cray_mpi_cross_product(struct vec3 v1, struct vec3 v2) {
  struct vec3 res;
  res.x = v1.y * v2.z - v1.z * v2.y;
  res.y = v1.z * v2.x - v1.x * v2.z;
  res.z = v1.x * v2.y - v1.y * v2.x;
  return res;
}


// ===========================================================================
// Jitter function taken from Graphics Gems I
// ===========================================================================
struct vec3 cray_mpi_jitter(int x, int y, int s, struct vec3 *urand, 
                            int *irand) {
  struct vec3 pt;
  pt.x = urand[(x + (y << 2) + irand[(x + s) & MASK]) & MASK].x;
  pt.y = urand[(y + (x << 2) + irand[(y + s) & MASK]) & MASK].y;
  return pt;
}


// ===========================================================================
// ===========================================================================
struct vec3 cray_mpi_get_sample_pos(int x, int y, int sample, int xres, 
                                    int yres, float aspect, 
                                    struct vec3 *urand, int *irand) {
  struct vec3 pt;
  float sf;

  sf = 1.5F / (float)xres;

  pt.x = ((float)x / (float)xres) - 0.5F;
  pt.y = -(((float)y / (float)yres) - 0.65F) / aspect;

  if(sample) {
    struct vec3 jt = cray_mpi_jitter(x, y, sample, urand, irand);
    pt.x += jt.x * sf;
    pt.y += jt.y * sf / aspect;
  }
  return pt;
}


// ===========================================================================
// Calculate ray-sphere intersection, and return {1, 0} to signify hit or no
// hit. Also the surface point parameters like position, normal, etc are
// returned through the sp pointer if it is not NULL.
// ===========================================================================
int cray_mpi_ray_sphere(const struct sphere *sph, struct ray ray, 
                        struct spoint *sp) {

  float a, b, c, d, sqrt_d, t1, t2;
  
  a = SQ(ray.dir.x) + SQ(ray.dir.y) + SQ(ray.dir.z);
  b = 2.0F * ray.dir.x * (ray.orig.x - sph->pos.x) +
      2.0F * ray.dir.y * (ray.orig.y - sph->pos.y) +
      2.0F * ray.dir.z * (ray.orig.z - sph->pos.z);
  c = SQ(sph->pos.x) + SQ(sph->pos.y) + SQ(sph->pos.z) +
      SQ(ray.orig.x) + SQ(ray.orig.y) + SQ(ray.orig.z) +
      2.0F * (-sph->pos.x * ray.orig.x - 
             sph->pos.y * ray.orig.y - 
             sph->pos.z * ray.orig.z) - 
      SQ(sph->rad);
  
  if((d = SQ(b) - 4.0F * a * c) < 0.0F) return 0;

  sqrt_d = ar_float_sqrt(d);
  t1 = (-b + sqrt_d) / (2.0F * a);
  t2 = (-b - sqrt_d) / (2.0F * a);

  if((t1 < ERR_MARGIN && t2 < ERR_MARGIN) || (t1 > 1.0F && t2 > 1.0F)) return 0;

  if(sp) {
    if(t1 < ERR_MARGIN) t1 = t2;
    if(t2 < ERR_MARGIN) t2 = t1;
    sp->dist = t1 < t2 ? t1 : t2;
    
    sp->pos.x = ray.orig.x + ray.dir.x * sp->dist;
    sp->pos.y = ray.orig.y + ray.dir.y * sp->dist;
    sp->pos.z = ray.orig.z + ray.dir.z * sp->dist;
    
    sp->normal.x = (sp->pos.x - sph->pos.x) / sph->rad;
    sp->normal.y = (sp->pos.y - sph->pos.y) / sph->rad;
    sp->normal.z = (sp->pos.z - sph->pos.z) / sph->rad;

    sp->vref = cray_mpi_reflect(ray.dir, sp->normal);
    NORMALIZE(sp->vref);
  }
  return 1;
}


// ===========================================================================
// Traces a ray throught the scene recursively (the recursion happens through
// shade() to calculate reflection rays if necessary).
// ===========================================================================
struct vec3 cray_mpi_trace(struct ray ray, int depth, struct sphere *objects, 
                           int onum, struct vec3 *lights, int lnum) {

  struct vec3 col;
  struct spoint sp, nearest_sp;
  struct sphere *nearest_obj = 0;
  int i;

  // if we reached the recursion limit, bail out
  if(depth >= MAX_RAY_DEPTH) {
    col.x = col.y = col.z = 0.0F;
    return col;
  }
  
  // find the nearest intersection ...
  for (i = 0; i < onum; i++) {
    if(cray_mpi_ray_sphere(objects + i, ray, &sp)) {
      if(!nearest_obj || sp.dist < nearest_sp.dist) {
        nearest_obj = objects + i;
        nearest_sp = sp;
      }
    }
  }

  // and perform shading calculations as needed by calling shade()
  if(nearest_obj) {
    col = cray_mpi_shade(nearest_obj, &nearest_sp, depth, objects, onum, 
                             lights, lnum);
  } else {
    col.x = col.y = col.z = 0.0F;
  }

  return col;
}


// ===========================================================================
// Calculates direct illumination with the phong reflectance model.
// Also handles reflections by calling trace again, if necessary.
// ===========================================================================
struct vec3 cray_mpi_shade(struct sphere *obj, struct spoint *sp, int depth,
                           struct sphere *objects, int onum, 
                           struct vec3 *lights, int lnum) {

  int i, j;
  struct vec3 col = {0, 0, 0};

  // for all lights ...
  for(i=0; i<lnum; i++) {
    float ispec, idiff;
    struct vec3 ldir;
    struct ray shadow_ray;
    int in_shadow = 0;

    ldir.x = lights[i].x - sp->pos.x;
    ldir.y = lights[i].y - sp->pos.y;
    ldir.z = lights[i].z - sp->pos.z;

    shadow_ray.orig = sp->pos;
    shadow_ray.dir = ldir;

    // shoot shadow rays to determine if we have a line of sight with the light
    for (j = 0; j < onum; j++) {
      if(cray_mpi_ray_sphere(objects + j, shadow_ray, 0)) {
        in_shadow = 1;
        break;
      }
    }

    // and if we're not in shadow, calculate direct illumination with the phong
    // model.
    if(!in_shadow) {
      NORMALIZE(ldir);

      idiff = MAX(DOT(sp->normal, ldir), 0.0F);
      //ispec = obj->mat.spow > 0.0F ? pow(MAX(DOT(sp->vref, ldir), 0.0F), obj->mat.spow) : 0.0F;
      // FIXME: pow implementation too slow, temporarily bypassing it
      ispec = 0.0F;

      col.x += idiff * obj->mat.col.x + ispec;
      col.y += idiff * obj->mat.col.y + ispec;
      col.z += idiff * obj->mat.col.z + ispec;
    }
  }

  // Also, if the object is reflective, spawn a reflection ray, and call trace()
  // to calculate the light arriving from the mirror direction.
  if(obj->mat.refl > 0.0F) {
    struct ray ray;
    struct vec3 rcol;

    ray.orig = sp->pos;
    ray.dir = sp->vref;
    ray.dir.x *= RAY_MAG;
    ray.dir.y *= RAY_MAG;
    ray.dir.z *= RAY_MAG;

    rcol = cray_mpi_trace(ray, depth + 1, objects, onum, lights, lnum);
    col.x += rcol.x * obj->mat.refl;
    col.y += rcol.y * obj->mat.refl;
    col.z += rcol.z * obj->mat.refl;
  }

  return col;
}


// ===========================================================================
// Determines the primary ray corresponding to the specified pixel (x, y)
// ===========================================================================
struct ray cray_mpi_get_primary_ray(int x, int y, int sample, int xres, 
                                    int yres, float aspect, struct vec3 *urand,
                                    int *irand, struct camera cam) {
  struct ray ray;
  float m[3][3];
  struct vec3 i, j = {0, 1, 0}, k, dir, orig, foo;

  k.x = cam.targ.x - cam.pos.x;
  k.y = cam.targ.y - cam.pos.y;
  k.z = cam.targ.z - cam.pos.z;
  NORMALIZE(k);

  i = cray_mpi_cross_product(j, k);
  j = cray_mpi_cross_product(k, i);
  m[0][0] = i.x; m[0][1] = j.x; m[0][2] = k.x;
  m[1][0] = i.y; m[1][1] = j.y; m[1][2] = k.y;
  m[2][0] = i.z; m[2][1] = j.z; m[2][2] = k.z;
  
  ray.orig.x = ray.orig.y = ray.orig.z = 0.0F;
  ray.dir = cray_mpi_get_sample_pos(x, y, sample, xres, yres, aspect, 
                                        urand, irand);
  ray.dir.z = 1.0F / HALF_FOV;
  ray.dir.x *= RAY_MAG;
  ray.dir.y *= RAY_MAG;
  ray.dir.z *= RAY_MAG;
  
  dir.x = ray.dir.x + ray.orig.x;
  dir.y = ray.dir.y + ray.orig.y;
  dir.z = ray.dir.z + ray.orig.z;
  foo.x = dir.x * m[0][0] + dir.y * m[0][1] + dir.z * m[0][2];
  foo.y = dir.x * m[1][0] + dir.y * m[1][1] + dir.z * m[1][2];
  foo.z = dir.x * m[2][0] + dir.y * m[2][1] + dir.z * m[2][2];

  orig.x = ray.orig.x * m[0][0] + ray.orig.y * m[0][1] + ray.orig.z * m[0][2] + cam.pos.x;
  orig.y = ray.orig.x * m[1][0] + ray.orig.y * m[1][1] + ray.orig.z * m[1][2] + cam.pos.y;
  orig.z = ray.orig.x * m[2][0] + ray.orig.y * m[2][1] + ray.orig.z * m[2][2] + cam.pos.z;

  ray.orig = orig;
  ray.dir.x = foo.x + orig.x;
  ray.dir.y = foo.y + orig.y;
  ray.dir.z = foo.z + orig.z;
  
  return ray;
}


// ===========================================================================
// ===========================================================================
void cray_mpi_render_scanline(int xsz, int ysz, int sl, unsigned int *fb, 
                              float aspect, struct vec3 *urand, int *irand,
                              struct sphere *objects, int onum, 
                              struct vec3 *lights, int lnum, 
                              struct camera cam) {

  int i, s;
  float rcp_samples = 1.0F / (float) RAYS_PER_PIXEL;
  float r, g, b;
  struct vec3 col;

  for(i=0; i<xsz; i++) {
    r = g = b = 0.0F;

    for(s=0; s<RAYS_PER_PIXEL; s++) {
      col = cray_mpi_trace(cray_mpi_get_primary_ray(i, sl, s, xsz, ysz, aspect,
                                                    urand, irand, cam), 
                           0, objects, onum, lights, lnum);
      r += col.x;
      g += col.y;
      b += col.z;
    }

    r = r * rcp_samples;
    g = g * rcp_samples;
    b = b * rcp_samples;
            
    fb[i] = ((unsigned int)(MIN(r, 1.0F) * 255.0F) & 0xff) << RSHIFT |
            ((unsigned int)(MIN(g, 1.0F) * 255.0F) & 0xff) << GSHIFT |
            ((unsigned int)(MIN(b, 1.0F) * 255.0F) & 0xff) << BSHIFT;
  }
}


// ===========================================================================
// Builds a scene same as in the "scene" benchmark input
// ===========================================================================
void cray_mpi_load_scene(struct sphere *obj_array, int *num_objs, 
                         struct vec3 *lights_array, int *num_lights, 
                         struct camera *cam) {

  *num_objs = 0;
  *num_lights = 0;

  obj_array[*num_objs].pos.x     = -1.5F; 
  obj_array[*num_objs].pos.y     = -0.3F;
  obj_array[*num_objs].pos.z     = -1.0F;
  obj_array[*num_objs].rad       = 0.7F;
  obj_array[*num_objs].mat.col.x = 1.0F;
  obj_array[*num_objs].mat.col.y = 0.2F;
  obj_array[*num_objs].mat.col.z = 0.05F;
  obj_array[*num_objs].mat.spow  = 50.0F;
  obj_array[*num_objs].mat.refl  = 0.3F;
  (*num_objs)++;

  obj_array[*num_objs].pos.x     = 1.5F; 
  obj_array[*num_objs].pos.y     = -0.4F;
  obj_array[*num_objs].pos.z     = 0.0F;
  obj_array[*num_objs].rad       = 0.6F;
  obj_array[*num_objs].mat.col.x = 0.1F;
  obj_array[*num_objs].mat.col.y = 0.85F;
  obj_array[*num_objs].mat.col.z = 1.0F;
  obj_array[*num_objs].mat.spow  = 50.0F;
  obj_array[*num_objs].mat.refl  = 0.4F;
  (*num_objs)++;

  obj_array[*num_objs].pos.x     = 0.0F; 
  obj_array[*num_objs].pos.y     = -1000.0F;
  obj_array[*num_objs].pos.z     = 2.0F;
  obj_array[*num_objs].rad       = 999.0F;
  obj_array[*num_objs].mat.col.x = 0.1F;
  obj_array[*num_objs].mat.col.y = 0.2F;
  obj_array[*num_objs].mat.col.z = 0.6F;
  obj_array[*num_objs].mat.spow  = 80.0F;
  obj_array[*num_objs].mat.refl  = 0.8F;
  (*num_objs)++;

  obj_array[*num_objs].pos.x     = 0.0F;
  obj_array[*num_objs].pos.y     = 0.0F;
  obj_array[*num_objs].pos.z     = 2.0F;
  obj_array[*num_objs].rad       = 1.0F;
  obj_array[*num_objs].mat.col.x = 1.0F;
  obj_array[*num_objs].mat.col.y = 0.5F;
  obj_array[*num_objs].mat.col.z = 0.1F;
  obj_array[*num_objs].mat.spow  = 60.0F;
  obj_array[*num_objs].mat.refl  = 0.7F;
  (*num_objs)++;

  lights_array[*num_lights].x    = -50.0F;
  lights_array[*num_lights].y    = 100.0F;
  lights_array[*num_lights].z    = -50.0F;
  (*num_lights)++;

  lights_array[*num_lights].x    = 40.0F;
  lights_array[*num_lights].y    = 40.0F;
  lights_array[*num_lights].z    = 150.0F;
  (*num_lights)++;

  (*cam).pos.x                   = 0.0F;
  (*cam).pos.y                   = 6.0F;
  (*cam).pos.z                   = -17.0F;
  (*cam).fov                     = 45.0F;
  (*cam).targ.x                  = 0.0F;
  (*cam).targ.y                  = -1.0F;
  (*cam).targ.z                  = 0.0F;
}


// ===========================================================================
// ===========================================================================
void cray_mpi_do_tile(unsigned int *tile_buffer, global_t *global,
                      int tile_y_offset, int lines_per_tile, int xres, 
                      int yres) {

  float aspect = (float) xres / (float) yres;
  int   i;

  for (i = 0; i < lines_per_tile; i++) {
    cray_mpi_render_scanline(xres, yres, tile_y_offset + i,
                             tile_buffer + i * xres, aspect, global->urand, 
                             global->irand, global->objects, global->onum, 
                             global->lights, global->lnum, global->cam);
  }
}




// ===========================================================================
// ===========================================================================
int cray_mpi(int num_procs,        // MPI processors to use
             int xres,             // x resolution (total picture cols)
             int yres              // y resolution (total picture rows)
            ) {

  unsigned int  seed;
  int           lines_per_core;
  global_t      *global;
  unsigned int  *pixels;
  unsigned int  time_start = 0;
  unsigned int  time_stop;
  unsigned int  time;
  int           num_cores;
  int           rank;
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
  if (yres % num_procs) {
    kt_printf("%d picture lines not divisible by %d cores\r\n",
              yres, num_procs);
    return 1;
  }
  lines_per_core = yres / num_procs;


  // Synchronize everyone and print infomercial
  MPI_Barrier(MPI_COMM_WORLD);
  if (!rank) {
    kt_printf("c-ray of %d x %d starting on %d core(s)\r\n",
              xres, yres, num_procs);
  }
  MPI_Barrier(MPI_COMM_WORLD);

  
  // Create per-core pixel buffers
  pixels = kt_malloc(xres * lines_per_core * sizeof(unsigned int));


  // Create per-core structure which will hold the scene objects and jitter
  // tables
  global = kt_malloc(sizeof(global_t));
    
  // Core 0 initializes global and broadcasts
  if (!rank) {
    // Build the scene
    cray_mpi_load_scene(global->objects, &(global->onum), 
                        global->lights, &(global->lnum), 
                        &(global->cam));

    // Initialize the random number tables for the jitter
    seed = 1;
    for (i = 0; i < NRAN; i++) {
      global->urand[i].x = 
                (float) ((seed = kt_rand(seed)) % 32768) / 32768.0F - 0.5F;
    }
    for(i = 0; i < NRAN; i++) {
      global->urand[i].y = 
                (float) ((seed = kt_rand(seed)) % 32768) / 32768.0F - 0.5F;
    }
    for(i = 0; i < NRAN; i++) {
      global->irand[i] = 
              (int) (NRAN * ((float)((seed = kt_rand(seed)) % 32768) / 32768));
    }
  }
  MPI_Bcast(global, sizeof(global_t), MPI_CHAR, 0, MPI_COMM_WORLD);


  // This kernel is reentrant for multiple MPI runs. If our core is not
  // part of the current num_procs setup, go to the next barrier directly.
  if (rank >= num_procs) {
    goto skip;
  }

  // Keep time
  if (!rank) {
    time_start = ar_free_timer_get_ticks();
  }

  // Run
  cray_mpi_do_tile(pixels, global, rank * lines_per_core, lines_per_core, 
                   xres, yres);
 
skip:

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
  }

  // Free stuff
  kt_free(pixels);
  kt_free(global);

  return 0;
}

