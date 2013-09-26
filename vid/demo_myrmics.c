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
// Abstract      : XUP-based video demo using Myrmics. This file contains the
//                 task-based code, run by the workers. Any low-level or
//                 arch/context-level functionality is separated and can
//                 be found in vid/demo_myrmics_vid.c.
//
// =============================[ CVS Variables ]=============================
// 
// File name     : $RCSfile: demo_myrmics.c,v $
// CVS revision  : $Revision: 1.16 $
// Last modified : $Date: 2013/03/22 12:25:41 $
// Last author   : $Author: lyberis-spree $
// 
// ===========================================================================

#include <myrmics.h>
#include <video.h>


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

struct vec3 demo_myrmics_shade(struct sphere *obj, struct spoint *sp, int depth,
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
#define RSHIFT  8
#define BSHIFT  24
#define GSHIFT  16       /* this is the same in both byte orders */

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
  unsigned int  vid_out_base;
  int           tiles_per_region;
  int           num_tiles;
} global_t;



// ===========================================================================
// Calculates reflection vector 
// ===========================================================================
struct vec3 demo_myrmics_reflect(struct vec3 v, struct vec3 n) {
  struct vec3 res;
  float dot = v.x * n.x + v.y * n.y + v.z * n.z;
  res.x = -(2.0F * dot * n.x - v.x);
  res.y = -(2.0F * dot * n.y - v.y);
  res.z = -(2.0F * dot * n.z - v.z);
  return res;
}


// ===========================================================================
// ===========================================================================
struct vec3 demo_myrmics_cross_product(struct vec3 v1, struct vec3 v2) {
  struct vec3 res;
  res.x = v1.y * v2.z - v1.z * v2.y;
  res.y = v1.z * v2.x - v1.x * v2.z;
  res.z = v1.x * v2.y - v1.y * v2.x;
  return res;
}


// ===========================================================================
// Jitter function taken from Graphics Gems I
// ===========================================================================
struct vec3 demo_myrmics_jitter(int x, int y, int s, struct vec3 *urand, 
                                int *irand) {
  struct vec3 pt;
  pt.x = urand[(x + (y << 2) + irand[(x + s) & MASK]) & MASK].x;
  pt.y = urand[(y + (x << 2) + irand[(y + s) & MASK]) & MASK].y;
  return pt;
}


// ===========================================================================
// ===========================================================================
struct vec3 demo_myrmics_get_sample_pos(int x, int y, int sample, int xres, 
                                        int yres, float aspect, 
                                        struct vec3 *urand, int *irand) {
  struct vec3 pt;
  float sf;

  sf = 1.5F / (float)xres;

  pt.x = ((float)x / (float)xres) - 0.5F;
  pt.y = -(((float)y / (float)yres) - 0.65F) / aspect;

  if(sample) {
    struct vec3 jt = demo_myrmics_jitter(x, y, sample, urand, irand);
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
int demo_myrmics_ray_sphere(const struct sphere *sph, struct ray ray, 
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

    sp->vref = demo_myrmics_reflect(ray.dir, sp->normal);
    NORMALIZE(sp->vref);
  }
  return 1;
}


// ===========================================================================
// Traces a ray throught the scene recursively (the recursion happens through
// shade() to calculate reflection rays if necessary).
// ===========================================================================
struct vec3 demo_myrmics_trace(struct ray ray, int depth, 
                               struct sphere *objects, 
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
    if(demo_myrmics_ray_sphere(objects + i, ray, &sp)) {
      if(!nearest_obj || sp.dist < nearest_sp.dist) {
        nearest_obj = objects + i;
        nearest_sp = sp;
      }
    }
  }

  // and perform shading calculations as needed by calling shade()
  if(nearest_obj) {
    col = demo_myrmics_shade(nearest_obj, &nearest_sp, depth, objects, onum, 
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
struct vec3 demo_myrmics_shade(struct sphere *obj, struct spoint *sp, int depth,
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
      if(demo_myrmics_ray_sphere(objects + j, shadow_ray, 0)) {
        in_shadow = 1;
        break;
      }
    }

    // and if we're not in shadow, calculate direct illumination with the phong
    // model.
    if(!in_shadow) {
      NORMALIZE(ldir);

      idiff = MAX(DOT(sp->normal, ldir), 0.0F);
      ispec = obj->mat.spow > 0.0F ? ar_float_pow(MAX(DOT(sp->vref, ldir), 0.0F),
                                                  obj->mat.spow) : 0.0F;

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

    rcol = demo_myrmics_trace(ray, depth + 1, objects, onum, lights, lnum);
    col.x += rcol.x * obj->mat.refl;
    col.y += rcol.y * obj->mat.refl;
    col.z += rcol.z * obj->mat.refl;
  }

  return col;
}


// ===========================================================================
// Determines the primary ray corresponding to the specified pixel (x, y)
// ===========================================================================
struct ray demo_myrmics_get_primary_ray(int x, int y, int sample, int xres, 
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

  i = demo_myrmics_cross_product(j, k);
  j = demo_myrmics_cross_product(k, i);
  m[0][0] = i.x; m[0][1] = j.x; m[0][2] = k.x;
  m[1][0] = i.y; m[1][1] = j.y; m[1][2] = k.y;
  m[2][0] = i.z; m[2][1] = j.z; m[2][2] = k.z;
  
  ray.orig.x = ray.orig.y = ray.orig.z = 0.0F;
  ray.dir = demo_myrmics_get_sample_pos(x, y, sample, xres, yres, aspect, 
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
void demo_myrmics_render_tile(int origin_x, int origin_y, int tile_width,
                              int tile_height, unsigned int *fb, float aspect, 
                              struct vec3 *urand, int *irand,
                              struct sphere *objects, int onum, 
                              struct vec3 *lights, int lnum, struct camera cam) {

  int x, y, s;
  float rcp_samples = 1.0F / (float) RAYS_PER_PIXEL;
  float r, g, b;
  struct vec3 col;

  for (y = origin_y; y < origin_y + tile_height; y++) {
    for (x = origin_x; x < origin_x + tile_width; x++) {
      r = g = b = 0.0F;

      for (s = 0; s < RAYS_PER_PIXEL; s++) {
        col = demo_myrmics_trace(demo_myrmics_get_primary_ray(x, y, s, 
                                                              VID_FULL_WIDTH, 
                                                              VID_STRIP_HEIGHT, 
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
              
      *fb++ = ((unsigned int)(MIN(r, 1.0F) * 255.0F) & 0xff) << RSHIFT |
              ((unsigned int)(MIN(g, 1.0F) * 255.0F) & 0xff) << GSHIFT |
              ((unsigned int)(MIN(b, 1.0F) * 255.0F) & 0xff) << BSHIFT;
    }
  }
}


// ===========================================================================
// Builds a scene same as in the "scene" benchmark input
// ===========================================================================
void demo_myrmics_load_scene(struct sphere *obj_array, int *num_objs, 
                             struct vec3 *lights_array, int *num_lights, 
                             struct camera *cam, int which_scene,
                             int red_present, int red_x,
                             int red_y, int green_present, int green_x,
                             int green_y) {
  float x;
  float z;
  int   i;
  float encore[] = {
    -2.405, 6.400, -1.162, 6.400,  1.069, 6.386,  2.088, 6.386, -0.072, 6.372,
     0.919, 6.359, -0.238, 6.345,  2.592, 6.345,  2.752, 6.345, -1.335, 6.331,
     1.918, 6.317,  1.706, 6.276, -2.586, 6.262, -2.286, 6.138, -1.689, 6.111,
     0.051, 6.097,  0.767, 6.069, -0.415, 6.028,  2.434, 5.973,  1.177, 5.959,
     2.865, 5.931, -1.491, 5.918, -1.077, 5.780, -2.761, 5.711,  1.795, 5.614,
    -2.176, 5.559,  0.625, 5.490,  0.151, 5.463, -0.553, 5.435,  2.963, 5.394,
     2.301, 5.366, -1.724, 5.339,  1.662, 5.297,  1.262, 5.270,  1.260, 5.256,
    -1.624, 5.242, -2.888, 4.980, -1.046, 4.925,  0.515, 4.746,  0.513, 4.732,
    -2.130, 4.677, -0.665, 4.608, -1.747, 4.553,  3.000, 4.553,  2.201, 4.511,
     1.296, 4.415,  1.631, 4.374, -2.973, 3.781, -2.821, 3.781, -2.688, 3.781,
    -2.530, 3.781, -2.363, 3.781, -2.228, 3.781, -2.109, 3.781, -0.717, 3.781,
     0.436, 3.781,  2.149, 3.781,  2.290, 3.781,  2.436, 3.781,  2.575, 3.781,
     2.717, 3.781,  2.867, 3.781,  2.992, 3.781, -1.067, 3.739, -1.793, 3.643,
     1.279, 3.422,  1.593, 3.395, -0.742, 2.871,  2.126, 2.871,  0.421, 2.843,
    -1.081, 2.816, -3.000, 2.788, -1.810, 2.609,  1.246, 2.609,  1.243, 2.595,
     1.572, 2.513,  2.161, 2.044, -2.213, 1.906, -1.098, 1.892,  0.451, 1.892,
    -0.704, 1.851, -2.950, 1.823,  1.156, 1.810,  2.886, 1.713, -1.835, 1.672,
     1.554, 1.575,  0.041, 1.561, -2.303, 1.217,  2.224, 1.217,  0.505, 1.203,
    -0.621, 1.134, -2.888, 1.024,  2.754, 1.024,  1.046, 0.996, -0.066, 0.941,
    -1.133, 0.817, -2.455, 0.803, -1.862, 0.776,  0.602, 0.762,  0.908, 0.734,
     2.317, 0.721,  1.518, 0.707, -0.523, 0.665,  2.596, 0.638, -2.763, 0.610,
    -0.232, 0.610,  2.444, 0.541,  0.729, 0.528, -2.613, 0.514, -0.386, 0.500
  };


  // Initialize
  *num_objs = 0;
  *num_lights = 0;

  // Balls: Z = -3.2f -> down, Z = 2.7f -> up
  //        X = -3.0f -> left, X = 3.0f -> right

  // Red ball
  if (red_present) {
    red_y = VID_STRIP_HEIGHT - 2 * DEMO_MYRMICS_BALL_SKETCH_RADIUS - red_y;
    obj_array[*num_objs].pos.x     = ((float) red_x * 5.5F) / 
                                      (float) (VID_FULL_WIDTH - 2 * 
                                               DEMO_MYRMICS_BALL_SKETCH_RADIUS)
                                      - 3.0f;
    obj_array[*num_objs].pos.y     = 0.0F;
    obj_array[*num_objs].pos.z     = ((float) red_y * 5.9F) / 
                                      (float) (VID_STRIP_HEIGHT - 2 * 
                                               DEMO_MYRMICS_BALL_SKETCH_RADIUS)
                                      - 3.2f;
    obj_array[*num_objs].rad       = 0.7F;
    obj_array[*num_objs].mat.col.x = 0.8F;
    obj_array[*num_objs].mat.col.y = 0.0F;
    obj_array[*num_objs].mat.col.z = 0.0F;
    obj_array[*num_objs].mat.spow  = 50.0F;
    obj_array[*num_objs].mat.refl  = 0.4F;
    (*num_objs)++;
  }

  // Green ball
  if (green_present) {
    green_y = VID_STRIP_HEIGHT - 2 * DEMO_MYRMICS_BALL_SKETCH_RADIUS - green_y;
    obj_array[*num_objs].pos.x     = ((float) green_x * 5.5F) / 
                                      (float) (VID_FULL_WIDTH - 2 * 
                                               DEMO_MYRMICS_BALL_SKETCH_RADIUS) 
                                      - 3.0f;
    obj_array[*num_objs].pos.y     = 0.0F;
    obj_array[*num_objs].pos.z     = ((float) green_y * 5.9F) / 
                                      (float) (VID_STRIP_HEIGHT - 2 * 
                                               DEMO_MYRMICS_BALL_SKETCH_RADIUS) 
                                      - 3.2f;
    obj_array[*num_objs].rad       = 0.7F;
    obj_array[*num_objs].mat.col.x = 0.0F;
    obj_array[*num_objs].mat.col.y = 0.6F;
    obj_array[*num_objs].mat.col.z = 0.0F;
    obj_array[*num_objs].mat.spow  = 30.0F;
    obj_array[*num_objs].mat.refl  = 0.5F;
    (*num_objs)++;
  }


  // Simple scene
  if (which_scene == 0) {

    // Put a black sphere in the center
    obj_array[*num_objs].pos.x     = 0.0F; 
    obj_array[*num_objs].pos.y     = 0.0F;
    obj_array[*num_objs].pos.z     = 0.0F;
    obj_array[*num_objs].rad       = 0.7F;
    obj_array[*num_objs].mat.col.x = 0.0F;
    obj_array[*num_objs].mat.col.y = 0.0F;
    obj_array[*num_objs].mat.col.z = 0.0F;
    obj_array[*num_objs].mat.spow  = 40.0F;
    obj_array[*num_objs].mat.refl  = 0.8F;
    (*num_objs)++;

    // Single large sphere as a floor
    obj_array[*num_objs].pos.x     = 0.0F; 
    obj_array[*num_objs].pos.y     = -1000.0F;
    obj_array[*num_objs].pos.z     = 2.0F;
    obj_array[*num_objs].rad       = 999.0F;
    obj_array[*num_objs].mat.col.x = 0.2F;
    obj_array[*num_objs].mat.col.y = 0.2F;
    obj_array[*num_objs].mat.col.z = 0.4F;
    obj_array[*num_objs].mat.spow  = 80.0F;
    obj_array[*num_objs].mat.refl  = 0.4F;
    (*num_objs)++;

    // 2 lights
    lights_array[*num_lights].x    = -50.0F;
    lights_array[*num_lights].y    = 50.0F;
    lights_array[*num_lights].z    = -50.0F;
    (*num_lights)++;

    lights_array[*num_lights].x    = 50.0F;
    lights_array[*num_lights].y    = 50.0F;
    lights_array[*num_lights].z    = 50.0F;
    (*num_lights)++;

  }

  // Medium-complexity scene
  else if (which_scene == 1) {

  // Small sphere floor (54 spheres)
    for (z = -1.0f; z <= 4.7f; z += 1.0f) {
      for (x = -2.75f; x <= 3.0f; x += 0.7f) {
        obj_array[*num_objs].pos.x     = x; 
        obj_array[*num_objs].pos.y     = -1.4F;
        obj_array[*num_objs].pos.z     = z;
        obj_array[*num_objs].rad       = 0.35F;
        obj_array[*num_objs].mat.col.x = 0.3F;
        obj_array[*num_objs].mat.col.y = 0.3F;
        obj_array[*num_objs].mat.col.z = 0.5F;
        obj_array[*num_objs].mat.spow  = 80.0F;
        obj_array[*num_objs].mat.refl  = 0.4F;
        (*num_objs)++;
        sys_assert(*num_objs <= MAX_OBJECTS);
      }
    }

    // 2 lights
    lights_array[*num_lights].x    = -50.0F;
    lights_array[*num_lights].y    = 50.0F;
    lights_array[*num_lights].z    = -50.0F;
    (*num_lights)++;

    lights_array[*num_lights].x    = 50.0F;
    lights_array[*num_lights].y    = 50.0F;
    lights_array[*num_lights].z    = 50.0F;
    (*num_lights)++;
  }


  // High-complexity scene
  else if (which_scene == 2) {

    // EnCORE logo (110 spheres)
    for (i = 0; i < 110; i++) {
      obj_array[*num_objs].pos.x     = encore[2*i] * 1.20f + 0.05f; 
      obj_array[*num_objs].pos.y     = -1.4F;
      obj_array[*num_objs].pos.z     = encore[2*i+1] + 1.50f;
      obj_array[*num_objs].rad       = 0.13F;
      if (encore[2*i] > -1.0f) {
        obj_array[*num_objs].mat.col.x = 0.65F;
        obj_array[*num_objs].mat.col.y = 0.0F;
        obj_array[*num_objs].mat.col.z = 0.0F;
      }
      else {
        obj_array[*num_objs].mat.col.x = 0.45F;
        obj_array[*num_objs].mat.col.y = 0.45F;
        obj_array[*num_objs].mat.col.z = 0.45F;
      }
      obj_array[*num_objs].mat.spow  = 90.0F;
      obj_array[*num_objs].mat.refl  = 0.4F;
      (*num_objs)++;
      sys_assert(*num_objs <= MAX_OBJECTS);
    }

    // Large sphere below logo
    obj_array[*num_objs].pos.x     = 0.0F; 
    obj_array[*num_objs].pos.y     = -21.7F;
    obj_array[*num_objs].pos.z     = 0.0F;
    obj_array[*num_objs].rad       = 19.65F;
    obj_array[*num_objs].mat.col.x = 0.8F;
    obj_array[*num_objs].mat.col.y = 0.8F;
    obj_array[*num_objs].mat.col.z = 1.0F;
    obj_array[*num_objs].mat.spow  = 0.0F;
    obj_array[*num_objs].mat.refl  = 0.8F;
    (*num_objs)++;

    // Very large sphere for horizon effect
    obj_array[*num_objs].pos.x     = 0.0F; 
    obj_array[*num_objs].pos.y     = -1045.0F;
    obj_array[*num_objs].pos.z     = 0.0F;
    obj_array[*num_objs].rad       = 999.65F;
    obj_array[*num_objs].mat.col.x = 0.6F;
    obj_array[*num_objs].mat.col.y = 0.6F;
    obj_array[*num_objs].mat.col.z = 1.0F;
    obj_array[*num_objs].mat.spow  = 0.0F;
    obj_array[*num_objs].mat.refl  = 0.0F;
    (*num_objs)++;

    // 2 lights
    lights_array[*num_lights].x    = -40.0F;
    lights_array[*num_lights].y    = 55.0F;
    lights_array[*num_lights].z    = -120.0F;
    (*num_lights)++;

    lights_array[*num_lights].x    = -40.0F;
    lights_array[*num_lights].y    = 55.0F;
    lights_array[*num_lights].z    = -120.0F;
    (*num_lights)++;
  }

  else {
    // Invalid scene
    sys_abort();
  }


  // Camera
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
void demo_myrmics_do_tile(unsigned int *tile_buffer, global_t *global,
                          int tile_id, int num_tiles) {

  float aspect;
  int   origin_x;
  int   origin_y;
  int   tile_width;
  int   tile_height;

  aspect = (float) VID_FULL_WIDTH / (float) VID_STRIP_HEIGHT;

  // Find out tile origin and size
  demo_myrmics_pix_buf_coords(tile_id, num_tiles, &origin_x,
                              &origin_y, &tile_width, &tile_height);

  // Render tile to local buffer
  demo_myrmics_render_tile(origin_x, origin_y, tile_width, tile_height,
                           tile_buffer, aspect, global->urand, 
                           global->irand, global->objects, global->onum, 
                           global->lights, global->lnum, global->cam);

  // Blit my boundaries
  vid_draw_line(tile_buffer, tile_width, tile_height, 
                0, 0, tile_width - 1, 0, 
                DEMO_MYRMICS_BOUNDARY_COLOR, DEMO_MYRMICS_BOUNDARY_TRANSP);
  vid_draw_line(tile_buffer, tile_width, tile_height, 
                0, 0, 0, tile_height - 1, 
                DEMO_MYRMICS_BOUNDARY_COLOR, DEMO_MYRMICS_BOUNDARY_TRANSP);
  if (origin_x + tile_width == VID_FULL_WIDTH) {
    vid_draw_line(tile_buffer, tile_width, tile_height, 
                  tile_width - 1, 0, tile_width - 1, tile_height - 1, 
                  DEMO_MYRMICS_BOUNDARY_COLOR, DEMO_MYRMICS_BOUNDARY_TRANSP);
  }
  if (origin_y + tile_height == VID_STRIP_HEIGHT) {
    vid_draw_line(tile_buffer, tile_width, tile_height, 
                  0, tile_height - 1, tile_width - 1, tile_height - 1, 
                  DEMO_MYRMICS_BOUNDARY_COLOR, DEMO_MYRMICS_BOUNDARY_TRANSP);
  }

  // Transfer buffer to XUP
  demo_myrmics_draw_output(tile_buffer, tile_width, tile_height,        
                           origin_x, origin_y, global->vid_out_base, tile_id,
                           num_tiles);
}


// ===========================================================================
// ===========================================================================
void demo_myrmics_do_region(rid_t r, unsigned int **pix_buf, int which_region,
                            global_t *global) {

  int tile_id;
  int num_tiles;
  int i;

  tile_id = global->tiles_per_region * which_region;
  num_tiles = global->num_tiles;
  for (i = 0; i < global->tiles_per_region; i++, tile_id++) {

    unsigned int *scoop_pix_buf = pix_buf[i];
    #pragma myrmics task inout(scoop_pix_buf) in(global) \
                          in(tile_id, num_tiles) safe(tile_id, num_tiles)
    demo_myrmics_do_tile(scoop_pix_buf, global, tile_id, num_tiles);
  }
}


// ===========================================================================
// ===========================================================================
void demo_myrmics_finish(rid_t r) {
  printf("all done, restarting\r\n");
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
void demo_myrmics(int num_regions) {

  rid_t                 r;
  int                   i;
  int                   j;
  unsigned int          seed;
  rid_t                 *regions;
  global_t              **global;
  unsigned int          ***pix_buf;
  int                   tile_width;
  int                   tile_height;
  int                   red_present;
  int                   red_x;
  int                   red_y;
  int                   green_present;
  int                   green_x;
  int                   green_y;
  int                   which_scene;


  // Create all-holding region
  r = sys_ralloc(0, 99); // highest level

  // Create regions
  regions = sys_alloc(num_regions * sizeof(rid_t), r);
  for (i = 0; i < num_regions; i++) {
    regions[i] = sys_ralloc(r, 0); // lowest level
  }

  // Allocate a global state structure per region
  global = sys_alloc(num_regions * sizeof(global_t *), r);
  for (i = 0; i < num_regions; i++) {
    global[i] = sys_alloc(sizeof(global_t), regions[i]);
  }

  // Communicate with video output and learn its worker buffer base address,
  // the red/green ball positions and which scene to render
  demo_myrmics_init(&global[0]->vid_out_base, &red_present, &red_x, &red_y,
                    &green_present, &green_x, &green_y, &which_scene, 
                    &global[0]->num_tiles);
  for (i = 1; i < num_regions; i++) {
    global[i]->vid_out_base = global[0]->vid_out_base;
    global[i]->num_tiles = global[0]->num_tiles;
  }

  // Compute tile size
  demo_myrmics_pix_buf_coords(0, global[0]->num_tiles, NULL, NULL, &tile_width, 
                              &tile_height);
  global[0]->tiles_per_region = global[0]->num_tiles / num_regions;
  sys_assert(num_regions * global[0]->tiles_per_region == global[0]->num_tiles);
  for (i = 1; i < num_regions; i++) {
    global[i]->tiles_per_region = global[0]->tiles_per_region;
  }

  // Create per-region pixel buffers
  pix_buf = sys_alloc(num_regions * sizeof(unsigned int **), r);
  for (i = 0; i < num_regions; i++) {
    pix_buf[i] = sys_alloc(global[0]->tiles_per_region * sizeof(unsigned int *),
                           regions[i]);
    if (global[0]->tiles_per_region > 1) {
      sys_balloc(tile_width * tile_height * sizeof(unsigned int), regions[i],
                 global[0]->tiles_per_region, pix_buf[i]);
    }
    else {
      pix_buf[i][0] = sys_alloc(tile_width * tile_height * sizeof(unsigned int),
                                regions[i]);
    }
  }

  for (i = 0; i < num_regions; i++) {

    // Build the scene
    demo_myrmics_load_scene(global[i]->objects, &(global[i]->onum), 
                            global[i]->lights, &(global[i]->lnum), 
                            &(global[i]->cam), which_scene, red_present, red_x, 
                            red_y, green_present, green_x, green_y);

    // Initialize the random number tables for the jitter
    seed = 1;
    for (j = 0; j < NRAN; j++) {
      global[i]->urand[j].x = 
                (float) ((seed = rand(seed)) % 32768) / 32768.0F - 0.5F;
    }
    for(j = 0; j < NRAN; j++) {
      global[i]->urand[j].y = 
                (float) ((seed = rand(seed)) % 32768) / 32768.0F - 0.5F;
    }
    for(j = 0; j < NRAN; j++) {
      global[i]->irand[j] = 
              (int) (NRAN * ((float)((seed = rand(seed)) % 32768) / 32768));
    }
  }

  // Run
  for (i = 0; i < num_regions; i++) {

    rid_t          scoop_region    = regions[i];
    unsigned int   **scoop_pix_buf = pix_buf[i];
    global_t       *scoop_global   = global[i];

    #pragma myrmics task \
                    region inout(scoop_region) \
                    in(scoop_pix_buf, i, scoop_global) \
                    safe(scoop_pix_buf, i, scoop_global)
    demo_myrmics_do_region(scoop_region, scoop_pix_buf, i, scoop_global);
  }

  #pragma myrmics task region inout(r)
  demo_myrmics_finish(r);
}


