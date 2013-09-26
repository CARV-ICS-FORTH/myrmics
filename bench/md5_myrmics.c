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
// The code in this file is derived by the public domain code by Alexander
// Peslyak. The code was under the following notice:
//
// This is an OpenSSL-compatible implementation of the RSA Data Security, Inc.
// MD5 Message-Digest Algorithm (RFC 1321).
//
// Homepage:
// http://openwall.info/wiki/people/solar/software/public-domain-source-code/md5
//
// Author:
// Alexander Peslyak, better known as Solar Designer <solar at openwall.com>
//
// This software was written by Alexander Peslyak in 2001.  No copyright is
// claimed, and the software is hereby placed in the public domain.
// In case this attempt to disclaim copyright and place the software in the
// public domain is deemed null and void, then the software is
// Copyright (c) 2001 Alexander Peslyak and it is hereby released to the
// general public under the following terms:
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted.
//
// There's ABSOLUTELY NO WARRANTY, express or implied.
//
// (This is a heavily cut-down "BSD license".)
//
// ==========================[ Static Information ]===========================
//
// Author        : Spyros Lyberis
// Abstract      : MD5 hash benchmark, Myrmics version
// =============================[ CVS Variables ]=============================
// 
// File name     : $RCSfile: md5_myrmics.c,v $
// CVS revision  : $Revision: 1.2 $
// Last modified : $Date: 2013/01/14 17:17:09 $
// Last author   : $Author: lyberis-spree $
// 
// ===========================================================================

#include <myrmics.h>

#define DIGEST_SIZE 16

typedef unsigned int MD5_u32plus;
typedef unsigned char uint8_t;

typedef struct {
  MD5_u32plus lo, hi;
  MD5_u32plus a, b, c, d;
  unsigned char buffer[64];
  MD5_u32plus block[16];
} MD5_CTX;


// The basic MD5 functions.
//
// F and G are optimized compared to their RFC 1321 definitions for
// architectures that lack an AND-NOT instruction, just like in Colin Plumb's
// implementation.
#define F(x, y, z)                      ((z) ^ ((x) & ((y) ^ (z))))
#define G(x, y, z)                      ((y) ^ ((z) & ((x) ^ (y))))
#define H(x, y, z)                      ((x) ^ (y) ^ (z))
#define I(x, y, z)                      ((y) ^ ((x) | ~(z)))

// The MD5 transformation for all four rounds.
#define STEP(f, a, b, c, d, x, t, s) \
        (a) += f((b), (c), (d)) + (x) + (t); \
        (a) = (((a) << (s)) | (((a) & 0xffffffff) >> (32 - (s)))); \
        (a) += (b);

// SET reads 4 input bytes in little-endian byte order and stores them
// in a properly aligned word in host byte order.
//
// The check for little-endian architectures that tolerate unaligned
// memory accesses is just an optimization.  Nothing will break if it
// doesn't work.
#define SET(n) \
        (*(MD5_u32plus *)&ptr[(n) * 4])
#define GET(n) \
        SET(n)


// ===========================================================================
// This processes one or more 64-byte data blocks, but does NOT update
// the bit counters.  There are no alignment requirements.
// ===========================================================================
void *md5_myrmics_body(MD5_CTX *ctx, void *data, unsigned int size) {

  unsigned char *ptr;
  MD5_u32plus a, b, c, d;
  MD5_u32plus saved_a, saved_b, saved_c, saved_d;

  ptr = data;

  a = ctx->a;
  b = ctx->b;
  c = ctx->c;
  d = ctx->d;

  do {
    saved_a = a;
    saved_b = b;
    saved_c = c;
    saved_d = d;

    /* Round 1 */
    STEP(F, a, b, c, d, SET(0), 0xd76aa478, 7)
    STEP(F, d, a, b, c, SET(1), 0xe8c7b756, 12)
    STEP(F, c, d, a, b, SET(2), 0x242070db, 17)
    STEP(F, b, c, d, a, SET(3), 0xc1bdceee, 22)
    STEP(F, a, b, c, d, SET(4), 0xf57c0faf, 7)
    STEP(F, d, a, b, c, SET(5), 0x4787c62a, 12)
    STEP(F, c, d, a, b, SET(6), 0xa8304613, 17)
    STEP(F, b, c, d, a, SET(7), 0xfd469501, 22)
    STEP(F, a, b, c, d, SET(8), 0x698098d8, 7)
    STEP(F, d, a, b, c, SET(9), 0x8b44f7af, 12)
    STEP(F, c, d, a, b, SET(10), 0xffff5bb1, 17)
    STEP(F, b, c, d, a, SET(11), 0x895cd7be, 22)
    STEP(F, a, b, c, d, SET(12), 0x6b901122, 7)
    STEP(F, d, a, b, c, SET(13), 0xfd987193, 12)
    STEP(F, c, d, a, b, SET(14), 0xa679438e, 17)
    STEP(F, b, c, d, a, SET(15), 0x49b40821, 22)

    /* Round 2 */
    STEP(G, a, b, c, d, GET(1), 0xf61e2562, 5)
    STEP(G, d, a, b, c, GET(6), 0xc040b340, 9)
    STEP(G, c, d, a, b, GET(11), 0x265e5a51, 14)
    STEP(G, b, c, d, a, GET(0), 0xe9b6c7aa, 20)
    STEP(G, a, b, c, d, GET(5), 0xd62f105d, 5)
    STEP(G, d, a, b, c, GET(10), 0x02441453, 9)
    STEP(G, c, d, a, b, GET(15), 0xd8a1e681, 14)
    STEP(G, b, c, d, a, GET(4), 0xe7d3fbc8, 20)
    STEP(G, a, b, c, d, GET(9), 0x21e1cde6, 5)
    STEP(G, d, a, b, c, GET(14), 0xc33707d6, 9)
    STEP(G, c, d, a, b, GET(3), 0xf4d50d87, 14)
    STEP(G, b, c, d, a, GET(8), 0x455a14ed, 20)
    STEP(G, a, b, c, d, GET(13), 0xa9e3e905, 5)
    STEP(G, d, a, b, c, GET(2), 0xfcefa3f8, 9)
    STEP(G, c, d, a, b, GET(7), 0x676f02d9, 14)
    STEP(G, b, c, d, a, GET(12), 0x8d2a4c8a, 20)

    /* Round 3 */
    STEP(H, a, b, c, d, GET(5), 0xfffa3942, 4)
    STEP(H, d, a, b, c, GET(8), 0x8771f681, 11)
    STEP(H, c, d, a, b, GET(11), 0x6d9d6122, 16)
    STEP(H, b, c, d, a, GET(14), 0xfde5380c, 23)
    STEP(H, a, b, c, d, GET(1), 0xa4beea44, 4)
    STEP(H, d, a, b, c, GET(4), 0x4bdecfa9, 11)
    STEP(H, c, d, a, b, GET(7), 0xf6bb4b60, 16)
    STEP(H, b, c, d, a, GET(10), 0xbebfbc70, 23)
    STEP(H, a, b, c, d, GET(13), 0x289b7ec6, 4)
    STEP(H, d, a, b, c, GET(0), 0xeaa127fa, 11)
    STEP(H, c, d, a, b, GET(3), 0xd4ef3085, 16)
    STEP(H, b, c, d, a, GET(6), 0x04881d05, 23)
    STEP(H, a, b, c, d, GET(9), 0xd9d4d039, 4)
    STEP(H, d, a, b, c, GET(12), 0xe6db99e5, 11)
    STEP(H, c, d, a, b, GET(15), 0x1fa27cf8, 16)
    STEP(H, b, c, d, a, GET(2), 0xc4ac5665, 23)

    /* Round 4 */
    STEP(I, a, b, c, d, GET(0), 0xf4292244, 6)
    STEP(I, d, a, b, c, GET(7), 0x432aff97, 10)
    STEP(I, c, d, a, b, GET(14), 0xab9423a7, 15)
    STEP(I, b, c, d, a, GET(5), 0xfc93a039, 21)
    STEP(I, a, b, c, d, GET(12), 0x655b59c3, 6)
    STEP(I, d, a, b, c, GET(3), 0x8f0ccc92, 10)
    STEP(I, c, d, a, b, GET(10), 0xffeff47d, 15)
    STEP(I, b, c, d, a, GET(1), 0x85845dd1, 21)
    STEP(I, a, b, c, d, GET(8), 0x6fa87e4f, 6)
    STEP(I, d, a, b, c, GET(15), 0xfe2ce6e0, 10)
    STEP(I, c, d, a, b, GET(6), 0xa3014314, 15)
    STEP(I, b, c, d, a, GET(13), 0x4e0811a1, 21)
    STEP(I, a, b, c, d, GET(4), 0xf7537e82, 6)
    STEP(I, d, a, b, c, GET(11), 0xbd3af235, 10)
    STEP(I, c, d, a, b, GET(2), 0x2ad7d2bb, 15)
    STEP(I, b, c, d, a, GET(9), 0xeb86d391, 21)

    a += saved_a;
    b += saved_b;
    c += saved_c;
    d += saved_d;

    ptr += 64;
  } while (size -= 64);

  ctx->a = a;
  ctx->b = b;
  ctx->c = c;
  ctx->d = d;

  return ptr;
}


// ===========================================================================
// ===========================================================================
void md5_myrmics_MD5_Init(MD5_CTX *ctx) {
  ctx->a = 0x67452301;
  ctx->b = 0xefcdab89;
  ctx->c = 0x98badcfe;
  ctx->d = 0x10325476;

  ctx->lo = 0;
  ctx->hi = 0;
}


// ===========================================================================
// ===========================================================================
void md5_myrmics_MD5_Update(MD5_CTX *ctx, void *data, unsigned int size) {
  MD5_u32plus saved_lo;
  unsigned int used, free;

  saved_lo = ctx->lo;
  if ((ctx->lo = (saved_lo + size) & 0x1fffffff) < saved_lo) {
    ctx->hi++;
  }
  ctx->hi += size >> 29;

  used = saved_lo & 0x3f;

  if (used) {
    free = 64 - used;

    if (size < free) {
      memcpy(&ctx->buffer[used], data, size);
      return;
    }

    memcpy(&ctx->buffer[used], data, free);
    data = (unsigned char *)data + free;
    size -= free;
    md5_myrmics_body(ctx, ctx->buffer, 64);
  }

  if (size >= 64) {
    data = md5_myrmics_body(ctx, data, size & ~(unsigned int)0x3f);
    size &= 0x3f;
  }

  memcpy(ctx->buffer, data, size);
}


// ===========================================================================
// ===========================================================================
void md5_myrmics_MD5_Final(unsigned char *result, MD5_CTX *ctx) {

  unsigned int used, free;

  used = ctx->lo & 0x3f;

  ctx->buffer[used++] = 0x80;

  free = 64 - used;

  if (free < 8) {
    memset(&ctx->buffer[used], 0, free);
    md5_myrmics_body(ctx, ctx->buffer, 64);
    used = 0;
    free = 64;
  }

  memset(&ctx->buffer[used], 0, free - 8);

  ctx->lo <<= 3;
  ctx->buffer[56] = ctx->lo;
  ctx->buffer[57] = ctx->lo >> 8;
  ctx->buffer[58] = ctx->lo >> 16;
  ctx->buffer[59] = ctx->lo >> 24;
  ctx->buffer[60] = ctx->hi;
  ctx->buffer[61] = ctx->hi >> 8;
  ctx->buffer[62] = ctx->hi >> 16;
  ctx->buffer[63] = ctx->hi >> 24;

  md5_myrmics_body(ctx, ctx->buffer, 64);

  result[0] = ctx->a;
  result[1] = ctx->a >> 8;
  result[2] = ctx->a >> 16;
  result[3] = ctx->a >> 24;
  result[4] = ctx->b;
  result[5] = ctx->b >> 8;
  result[6] = ctx->b >> 16;
  result[7] = ctx->b >> 24;
  result[8] = ctx->c;
  result[9] = ctx->c >> 8;
  result[10] = ctx->c >> 16;
  result[11] = ctx->c >> 24;
  result[12] = ctx->d;
  result[13] = ctx->d >> 8;
  result[14] = ctx->d >> 16;
  result[15] = ctx->d >> 24;

  memset(ctx, 0, sizeof(*ctx));
}


// ===========================================================================
// Processes one input buffer, delivering the digest into out.
// ===========================================================================
void md5_myrmics_process(uint8_t *in, uint8_t *out, int bufsize) {
  MD5_CTX context;
  uint8_t digest[DIGEST_SIZE];
  
  md5_myrmics_MD5_Init(&context);
  md5_myrmics_MD5_Update(&context, in, bufsize);
  md5_myrmics_MD5_Final(digest, &context);

  memcpy(out, digest, DIGEST_SIZE);
}


// ===========================================================================
// Processes all buffers in a region
// ===========================================================================
void myrmics_md5_do_region(rid_t region, uint8_t **in_buffers, 
                           uint8_t **out_buffers, int buffers_per_region, 
                           int bufsize) {

  int i;

  for (i = 0; i < buffers_per_region; i++) {
    uint8_t     *scoop_in_buffer  = in_buffers[i];
    uint8_t     *scoop_out_buffer = out_buffers[i];

    #pragma myrmics task in(scoop_in_buffer) out(scoop_out_buffer) \
                         in(bufsize) safe(bufsize)
    md5_myrmics_process(scoop_in_buffer, scoop_out_buffer, bufsize);
  }
}


// ===========================================================================
// ===========================================================================
void md5_myrmics_print(rid_t region, unsigned int time_start) {

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
void md5_myrmics(int num_buffers, int bufsize, int num_regions) {

  unsigned int  seed = 1;
  rid_t         r;
  rid_t         *regions;
  int           buffers_per_region;
  uint8_t       ***in_buffers;  // [region][buffer]
  uint8_t       ***out_buffers; // [region][buffer]
  unsigned int  time_start;
  uint8_t       tmp;
  int           tenth;
  int           i;
  int           j;
  int           k;


  // Sanity checks
  if (num_buffers % num_regions) {
    printf("%d buffers not divisible by %d regions\r\n", 
           num_buffers, num_regions);
    return;
  }
  buffers_per_region = num_buffers / num_regions;

  // Infomercial
  printf("MD5 of %d buffers starting split into %d region(s)\r\n",
         num_buffers, num_regions);

  // Create all-holding region
  r = sys_ralloc(0, 99); // highest level

  // Create regions
  regions = sys_alloc(num_regions * sizeof(rid_t), r);
  for (i = 0; i < num_regions; i++) {
    regions[i] = sys_ralloc(r, 0); // lowest level
  }

  // Create per-region input and output buffers
  in_buffers = sys_alloc(num_regions * sizeof(uint8_t **), r);
  out_buffers = sys_alloc(num_regions * sizeof(uint8_t **), r);
  for (i = 0; i < num_regions; i++) {
    in_buffers[i] = sys_alloc(buffers_per_region * sizeof(uint8_t *),
                          regions[i]);
    sys_balloc(bufsize * sizeof(uint8_t), regions[i],
               buffers_per_region, in_buffers[i]);

    out_buffers[i] = sys_alloc(buffers_per_region * sizeof(uint8_t *),
                          regions[i]);
    sys_balloc(DIGEST_SIZE * sizeof(uint8_t), regions[i],
               buffers_per_region, out_buffers[i]);
  }

  // Generate random input buffers
  tenth = 0;
  for (k = 0; k < bufsize; k++) {
    tmp = (seed = rand(seed)) & 0xFF;
    for (i = 0; i < num_regions; i++) {
      for (j = 0; j < buffers_per_region; j++) {
        in_buffers[i][j][k] = tmp;
      }
    }
    if (k * 10 / bufsize > tenth + 1) {
      tenth++;
      printf("serial init %d%% done\r\n", tenth * 10);
    }
  }

  // Start time
  printf("starting parallel phase...\r\n");
  time_start = sys_free_timer_get_ticks();

  // Spawn all region tasks
  for (i = 0; i < num_regions; i++) {

    rid_t       scoop_region        = regions[i];
    uint8_t     **scoop_in_buffers  = in_buffers[i];
    uint8_t     **scoop_out_buffers = out_buffers[i];

    #pragma myrmics task region inout(scoop_region) \
                         in(scoop_in_buffers, scoop_out_buffers, \
                            buffers_per_region, bufsize) \
                         safe(scoop_in_buffers, scoop_out_buffers, \
                            buffers_per_region, bufsize)
    myrmics_md5_do_region(scoop_region, scoop_in_buffers, scoop_out_buffers,
                          buffers_per_region, bufsize);
  }

  // Spawn task to wait for everyone and print the elapsed time
  #pragma myrmics task region inout(r) in(time_start) safe(time_start)
  md5_myrmics_print(r, time_start);

}
