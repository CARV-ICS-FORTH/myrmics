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
// Abstract      : Raytracing filter, Myrmics version
// =============================[ CVS Variables ]=============================
// 
// File name     : $RCSfile: cray_myrmics.c,v $
// CVS revision  : $Revision: 1.3 $
// Last modified : $Date: 2013/01/14 17:17:09 $
// Last author   : $Author: lyberis-spree $
// 
// ===========================================================================

#include <myrmics.h>


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
        float len = sqrt(DOT(a, a));\
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
} per_region_t;



// ===========================================================================
// Calculates reflection vector 
// ===========================================================================
struct vec3 cray_myrmics_reflect(struct vec3 v, struct vec3 n) {
  struct vec3 res;
  float dot = v.x * n.x + v.y * n.y + v.z * n.z;
  res.x = -(2.0F * dot * n.x - v.x);
  res.y = -(2.0F * dot * n.y - v.y);
  res.z = -(2.0F * dot * n.z - v.z);
  return res;
}


// ===========================================================================
// ===========================================================================
struct vec3 cray_myrmics_cross_product(struct vec3 v1, struct vec3 v2) {
  struct vec3 res;
  res.x = v1.y * v2.z - v1.z * v2.y;
  res.y = v1.z * v2.x - v1.x * v2.z;
  res.z = v1.x * v2.y - v1.y * v2.x;
  return res;
}


// ===========================================================================
// Jitter function taken from Graphics Gems I
// ===========================================================================
struct vec3 cray_myrmics_jitter(int x, int y, int s, struct vec3 *urand, 
                                int *irand) {
  struct vec3 pt;
  pt.x = urand[(x + (y << 2) + irand[(x + s) & MASK]) & MASK].x;
  pt.y = urand[(y + (x << 2) + irand[(y + s) & MASK]) & MASK].y;
  return pt;
}


// ===========================================================================
// ===========================================================================
struct vec3 cray_myrmics_get_sample_pos(int x, int y, int sample, int xres, 
                                        int yres, float aspect, 
                                        struct vec3 *urand, int *irand) {
  struct vec3 pt;
  float sf;

  sf = 1.5F / (float)xres;

  pt.x = ((float)x / (float)xres) - 0.5F;
  pt.y = -(((float)y / (float)yres) - 0.65F) / aspect;

  if(sample) {
    struct vec3 jt = cray_myrmics_jitter(x, y, sample, urand, irand);
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
int cray_myrmics_ray_sphere(const struct sphere *sph, struct ray ray, 
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

  sqrt_d = sqrt(d);
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

    sp->vref = cray_myrmics_reflect(ray.dir, sp->normal);
    NORMALIZE(sp->vref);
  }
  return 1;
}


// ===========================================================================
// Traces a ray throught the scene recursively (the recursion happens through
// shade() to calculate reflection rays if necessary).
// ===========================================================================
struct vec3 cray_myrmics_trace(struct ray ray, int depth, 
                               struct sphere *objects, int onum, 
                               struct vec3 *lights, int lnum) {

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
    if(cray_myrmics_ray_sphere(objects + i, ray, &sp)) {
      if(!nearest_obj || sp.dist < nearest_sp.dist) {
        nearest_obj = objects + i;
        nearest_sp = sp;
      }
    }
  }

  // and perform shading calculations as needed by calling shade()
  if(nearest_obj) {
    col = cray_myrmics_shade(nearest_obj, &nearest_sp, depth, objects, onum, 
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
struct vec3 cray_myrmics_shade(struct sphere *obj, struct spoint *sp, int depth,
                  struct sphere *objects, int onum, struct vec3 *lights, 
                  int lnum) {

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
      if(cray_myrmics_ray_sphere(objects + j, shadow_ray, 0)) {
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

    rcol = cray_myrmics_trace(ray, depth + 1, objects, onum, lights, lnum);
    col.x += rcol.x * obj->mat.refl;
    col.y += rcol.y * obj->mat.refl;
    col.z += rcol.z * obj->mat.refl;
  }

  return col;
}


// ===========================================================================
// Determines the primary ray corresponding to the specified pixel (x, y)
// ===========================================================================
struct ray cray_myrmics_get_primary_ray(int x, int y, int sample, int xres, 
                                        int yres, float aspect, 
                                        struct vec3 *urand, int *irand,
                                        struct camera cam) {
  struct ray ray;
  float m[3][3];
  struct vec3 i, j = {0, 1, 0}, k, dir, orig, foo;

  k.x = cam.targ.x - cam.pos.x;
  k.y = cam.targ.y - cam.pos.y;
  k.z = cam.targ.z - cam.pos.z;
  NORMALIZE(k);

  i = cray_myrmics_cross_product(j, k);
  j = cray_myrmics_cross_product(k, i);
  m[0][0] = i.x; m[0][1] = j.x; m[0][2] = k.x;
  m[1][0] = i.y; m[1][1] = j.y; m[1][2] = k.y;
  m[2][0] = i.z; m[2][1] = j.z; m[2][2] = k.z;
  
  ray.orig.x = ray.orig.y = ray.orig.z = 0.0F;
  ray.dir = cray_myrmics_get_sample_pos(x, y, sample, xres, yres, aspect, 
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
void cray_myrmics_render_scanline(int xsz, int ysz, int sl, unsigned int *fb, 
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
      col = cray_myrmics_trace(cray_myrmics_get_primary_ray(i, sl, s, xsz, ysz,
                                                            aspect, urand, 
                                                            irand, cam), 
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
void cray_myrmics_load_scene(struct sphere *obj_array, int *num_objs, 
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
void cray_myrmics_do_tile(unsigned int *tile_buffer, per_region_t *per_region,
                          int tile_y_offset, int lines_per_tile, int xres, 
                          int yres) {

  float aspect = (float) xres / (float) yres;
  int   i;

  for (i = 0; i < lines_per_tile; i++) {
    cray_myrmics_render_scanline(xres, yres, tile_y_offset + i,
                                 tile_buffer + i * xres, aspect, 
                                 per_region->urand, per_region->irand,
                                 per_region->objects, per_region->onum, 
                                 per_region->lights, per_region->lnum, 
                                 per_region->cam);
  }
}


// ===========================================================================
// ===========================================================================
void cray_myrmics_do_region(rid_t region, unsigned int **region_buffers, 
                            int region_idx, int lines_per_region, 
                            int lines_per_tile, int tiles_per_region, int xres, 
                            int yres, per_region_t *per_region) {

  unsigned int  *tile_buffer;
  int           tile_y_offset;
  int           i;

  for (i = 0; i < tiles_per_region; i++) {

    tile_y_offset = region_idx * lines_per_region + i * lines_per_tile;
    tile_buffer   = region_buffers[i];

    #pragma myrmics task inout(tile_buffer) \
                         in(per_region) \
                         in(tile_y_offset, lines_per_tile, xres, yres) \
                         safe(tile_y_offset, lines_per_tile, xres, yres)
    cray_myrmics_do_tile(tile_buffer, per_region, tile_y_offset, 
                         lines_per_tile, xres, yres);
  }
}


// ===========================================================================
// ===========================================================================
void cray_myrmics_print(rid_t region, unsigned int time_start) {

  unsigned int time_stop;
  unsigned int time;


  // Compute elapsed time
  time_stop = sys_free_timer_get_ticks();
  if (time_stop > time_start) {
    time = time_stop - time_start;
  }
  else {
    time = 0xFFFFFFFF - (time_start - time_stop);
  }
  printf("Time: %10u cycles (%6u msec)\r\n", time, time / 10000);
}


// ===========================================================================
// ===========================================================================
void cray_myrmics(int xres,             // x resolution (total picture cols)
                  int yres,             // y resolution (total picture rows)
                  int num_regions,      // regions to split yres
                  int tiles_per_region  // tiles per region for further split
                 ) {

  unsigned int  seed;
  int           lines_per_region;
  int           lines_per_tile;
  rid_t         r;
  rid_t         *regions;
  per_region_t  **per_region; // [region]
  unsigned int  ***pixels; // [region][tile]
  unsigned int  time_start;
  int           i;
  int           j;


  // Sanity checks
  if (yres % num_regions) {
    printf("%d picture lines not divisible by %d regions\r\n",
           yres, num_regions);
    return;
  }
  lines_per_region = yres / num_regions;

  if (lines_per_region % tiles_per_region) {
    printf("%d lines per region not divisible by %d tiles\r\n",
           lines_per_region, tiles_per_region);
    return;
  }
  lines_per_tile = lines_per_region / tiles_per_region;


  // Infomercial
  printf("c-ray of %d x %d starting in %d total tile(s) and %d region(s)\r\n",
         xres, yres, tiles_per_region * num_regions, num_regions);


  // Create all-holding region
  r = sys_ralloc(0, 99); // highest level

  // Create regions
  regions = sys_alloc(num_regions * sizeof(rid_t), r);
  for (i = 0; i < num_regions; i++) {
    regions[i] = sys_ralloc(r, 0); // lowest level
  }

  // Create per-region pixel buffers
  pixels = sys_alloc(num_regions * sizeof(unsigned int **), r);
  for (i = 0; i < num_regions; i++) {
    pixels[i] = sys_alloc(tiles_per_region * sizeof(unsigned int *),
                          regions[i]);
    if (tiles_per_region > 1) {
      sys_balloc(xres * lines_per_tile * sizeof(unsigned int), regions[i],
                 tiles_per_region, pixels[i]);
    }
    else {
      pixels[i][0] = sys_alloc(xres * lines_per_tile * sizeof(unsigned int), 
                               regions[i]);
    }
  }

  // Create per-region structure which will hold the scene objects and jitter
  // tables, so they are private per leaf scheduler and we don't need to 
  // involve the top-level scheduler for leaf tasks
  per_region = sys_alloc(num_regions * sizeof(per_region_t *), r);
  for (i = 0; i < num_regions; i++) {
    per_region[i] = sys_alloc(sizeof(per_region_t), regions[i]);
    
    // Build the scene
    cray_myrmics_load_scene(per_region[i]->objects, &(per_region[i]->onum), 
                            per_region[i]->lights, &(per_region[i]->lnum), 
                            &(per_region[i]->cam));

    // Initialize the random number tables for the jitter
    seed = 1;
    for (j = 0; j < NRAN; j++) {
      per_region[i]->urand[j].x = 
                (float) ((seed = rand(seed)) % 32768) / 32768.0F - 0.5F;
    }
    for(j = 0; j < NRAN; j++) {
      per_region[i]->urand[j].y = 
                (float) ((seed = rand(seed)) % 32768) / 32768.0F - 0.5F;
    }
    for(j = 0; j < NRAN; j++) {
      per_region[i]->irand[j] = 
                (int) (NRAN * ((float)((seed = rand(seed)) % 32768) / 32768));
    }
  }


  // Start time
  time_start = sys_free_timer_get_ticks();

  // Spawn all region tasks
  for (i = 0; i < num_regions; i++) {

    rid_t          scoop_region      = regions[i];
    unsigned int   **scoop_pixels    = pixels[i];
    per_region_t   *scoop_per_region = per_region[i];

    #pragma myrmics task \
                    region inout(scoop_region) \
                    in(scoop_pixels, i, lines_per_region, lines_per_tile,\
                       tiles_per_region, xres, yres, scoop_per_region) \
                    safe(scoop_pixels, i, lines_per_region, lines_per_tile,\
                       tiles_per_region, xres, yres, scoop_per_region)
    cray_myrmics_do_region(scoop_region, scoop_pixels, i, lines_per_region, 
                           lines_per_tile, tiles_per_region, xres, yres, 
                           scoop_per_region);
  }

  // Spawn task to wait for everyone and print the elapsed time
  #pragma myrmics task region inout(r) in(time_start) safe(time_start)
  cray_myrmics_print(r, time_start);
 
}

