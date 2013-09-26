/* -*- C -*-
 *
 * Copyright 2006 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, 
 * Boston, MA  02110-1301, USA.
 */

#include <kernel_toolset.h>
#include <arch.h>
#include <fmpi.h>

/* constants */
const int magic_tag = 1;


static void smb_msgrate_cache_invalidate(int *cache_buf, int cache_size) {
    int i;

    cache_buf[0] = 1;
    for (i = 1 ; i < cache_size ; ++i) {
        cache_buf[i] = cache_buf[i - 1];
    }
}


static void smb_msgrate_display_result(int rank, int machine_output, const char *test, const float result)
{
    if (0 == rank) {
        if (machine_output) {
            kt_printf("%f ", result);
        } else {
            kt_printf("%10s: %f messages per 1000cc\n", test, result);
        }
    }
}


void smb_msgrate_test_one_way(int world_size, int rank, int niters,
                              int *cache_buf, int cache_size, int nmsgs,
                              char *send_buf, int nbytes, MPI_Request *reqs,
                              MPI_Status *statuses, char *recv_buf, 
                              int machine_output) {

    int i, k, nreqs;
    unsigned int *tmp, *total;
    float rate;
    tmp   = kt_malloc(sizeof(unsigned int));
    total = kt_malloc(sizeof(unsigned int));

    *tmp = 0;
    *total = 0;


    MPI_Barrier(MPI_COMM_WORLD);

    ar_assert(world_size % 2 == 0);

    if (!(world_size % 2 == 1 && rank == (world_size - 1))) {
        if (rank < world_size / 2) {
            for (i = 0 ; i < niters ; ++i) {
                smb_msgrate_cache_invalidate(cache_buf, cache_size);

                MPI_Barrier(MPI_COMM_WORLD);

                ar_timer_reset();
                nreqs = 0;
                for (k = 0 ; k < nmsgs ; ++k) {
                    MPI_Isend(send_buf + (nbytes * k),
                              nbytes, MPI_CHAR, rank + (world_size / 2), 
                              magic_tag,
                              //magic_tag + k + i * nmsgs, 
                              MPI_COMM_WORLD, &reqs[nreqs++]);
                }
                MPI_Waitall(nreqs, reqs, statuses);
                *total += ar_timer_get_cycles();
            }
        } else {
            for (i = 0 ; i < niters ; ++i) {
                smb_msgrate_cache_invalidate(cache_buf, cache_size);

                MPI_Barrier(MPI_COMM_WORLD);

                ar_timer_reset();
                nreqs = 0;
                for (k = 0 ; k < nmsgs ; ++k) {
                    MPI_Irecv(recv_buf + (nbytes * k),
                              nbytes, MPI_CHAR, rank - (world_size / 2), 
                              magic_tag,
                              //magic_tag + k + i * nmsgs, 
                              //MPI_ANY_TAG,
                              MPI_COMM_WORLD, &reqs[nreqs++]);
                }
                MPI_Waitall(nreqs, reqs, statuses);
                *total += ar_timer_get_cycles();
            }
        }

        //MPI_Allreduce(&total, &tmp, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

        MPI_Reduce(total, tmp, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        if (0 == rank) {
          rate = (float) (niters * nmsgs) / (float) (*tmp / world_size / 1000);
          smb_msgrate_display_result(rank, machine_output, "one-way", rate);
        }


        kt_free(total);
        kt_free(tmp);

    }

    MPI_Barrier(MPI_COMM_WORLD);
}


void smb_msgrate_test_same_direction(int niters, int *cache_buf, int cache_size,
                                     int npeers, int nmsgs, char *recv_buf,
                                     int nbytes, int *recv_peers, 
                                     MPI_Request *reqs, MPI_Status *statuses,
                                     char *send_buf, int *send_peers, int rank,
                                     int machine_output, int world_size) {

    int i, j, k, nreqs;
    unsigned int *tmp, *total;
    float rate;
    tmp   = kt_malloc(sizeof(unsigned int));
    total = kt_malloc(sizeof(unsigned int));

    *tmp = 0;
    *total = 0;

    MPI_Barrier(MPI_COMM_WORLD);

    for (i = 0 ; i < niters ; ++i) {
        smb_msgrate_cache_invalidate(cache_buf, cache_size);

        MPI_Barrier(MPI_COMM_WORLD);

        ar_timer_reset();
        for (j = 0 ; j < npeers ; ++j) {
            nreqs = 0;
            for (k = 0 ; k < nmsgs ; ++k) {
                MPI_Irecv(recv_buf + (nbytes * (k + j * nmsgs)),
                          nbytes, MPI_CHAR, recv_peers[j], 
                          magic_tag, 
                          //magic_tag + k + i * nmsgs, 
                          //MPI_ANY_TAG,
                          MPI_COMM_WORLD, &reqs[nreqs++]);
            }
            for (k = 0 ; k < nmsgs ; ++k) {
                MPI_Isend(send_buf + (nbytes * (k + j * nmsgs)),
                          nbytes, MPI_CHAR, send_peers[npeers - j - 1], 
                          magic_tag, 
                          //magic_tag + k + i * nmsgs, 
                          MPI_COMM_WORLD, &reqs[nreqs++]);
            }
            MPI_Waitall(nreqs, reqs, statuses);
        }
        *total += ar_timer_get_cycles();
    }

    //MPI_Allreduce(&total, &tmp, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    //smb_msgrate_display_result(rank, machine_output, "pair-based", (niters * npeers * nmsgs * 2) / (tmp / world_size));
    
    MPI_Reduce(total, tmp, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    if (0 == rank) {
      rate = (float) (niters * npeers * nmsgs * 2) / (float) (*tmp / world_size / 1000);
      smb_msgrate_display_result(rank, machine_output, "pair-based", rate);
    }


    kt_free(total);
    kt_free(tmp);
}


void smb_msgrate_test_prepost(int npeers, int nmsgs, char *recv_buf, 
                              int nbytes, int *recv_peers, MPI_Request *reqs,
                              int niters, int *cache_buf, int cache_size,
                              char *send_buf, int *send_peers, 
                              MPI_Status *statuses, int rank, 
                              int machine_output, int world_size) {

    int i, j, k, nreqs = 0;
    unsigned int *tmp, *total;
    float rate;
    tmp   = kt_malloc(sizeof(unsigned int));
    total = kt_malloc(sizeof(unsigned int));

    *tmp = 0;
    *total = 0;

    MPI_Barrier(MPI_COMM_WORLD);

    ar_timer_reset();
    for (j = 0 ; j < npeers ; ++j) {
        for (k = 0 ; k < nmsgs ; ++k) {
            MPI_Irecv(recv_buf + (nbytes * (k + j * nmsgs)),
                      nbytes, MPI_CHAR, recv_peers[j], 
                      magic_tag, 
                      //magic_tag + k,
                      //MPI_ANY_TAG,
                      MPI_COMM_WORLD, &reqs[nreqs++]);
        }
    }
    *total += ar_timer_get_cycles();

    for (i = 0 ; i < niters - 1 ; ++i) {
        smb_msgrate_cache_invalidate(cache_buf, cache_size);

        MPI_Barrier(MPI_COMM_WORLD);

        ar_timer_reset();
        for (j = 0 ; j < npeers ; ++j) {
            for (k = 0 ; k < nmsgs ; ++k) {
                MPI_Isend(send_buf + (nbytes * (k + j * nmsgs)),
                          nbytes, MPI_CHAR, send_peers[npeers - j - 1], 
                          magic_tag, 
                          //magic_tag + k + i * nmsgs, 
                          MPI_COMM_WORLD, &reqs[nreqs++]);
            }
        }
        MPI_Waitall(nreqs, reqs, statuses);
        nreqs = 0;
        for (j = 0 ; j < npeers ; ++j) {
            for (k = 0 ; k < nmsgs ; ++k) {
                MPI_Irecv(recv_buf + (nbytes * (k + j * nmsgs)),
                          nbytes, MPI_CHAR, recv_peers[j], 
                          magic_tag, 
                          //magic_tag + k + (i + 1) * nmsgs, 
                          //MPI_ANY_TAG,
                          MPI_COMM_WORLD, &reqs[nreqs++]);
            }
        }
        *total += ar_timer_get_cycles();
    }
    ar_timer_reset();
    for (j = 0 ; j < npeers ; ++j) {
        for (k = 0 ; k < nmsgs ; ++k) {
            MPI_Isend(send_buf + (nbytes * (k + j * nmsgs)),
                      nbytes, MPI_CHAR, send_peers[npeers - j - 1], 
                      magic_tag, 
                      //magic_tag + k + i * nmsgs, 
                      MPI_COMM_WORLD, &reqs[nreqs++]);
        }
    }
    MPI_Waitall(nreqs, reqs, statuses);
    *total += ar_timer_get_cycles();

    //MPI_Allreduce(&total, &tmp, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    //smb_msgrate_display_result(rank, machine_output, "pre-post", (niters * npeers * nmsgs * 2) / (tmp / world_size));

    MPI_Reduce(total, tmp, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    if (0 == rank) {
      rate = (float) (niters * npeers * nmsgs * 2) / (float) (*tmp / world_size / 1000);
      smb_msgrate_display_result(rank, machine_output, "pre-post", rate);
    }


    kt_free(total);
    kt_free(tmp);
}


void smb_msgrate_test_allstart(int niters, int *cache_buf, int cache_size,
                               int npeers, int nmsgs, char *recv_buf,
                               int nbytes, int *recv_peers, MPI_Request *reqs,
                               char *send_buf, int *send_peers, 
                               MPI_Status *statuses, int rank, 
                               int machine_output, int world_size) {

    int i, j, k, nreqs = 0;
    unsigned int *tmp, *total;
    float rate;
    tmp   = kt_malloc(sizeof(unsigned int));
    total = kt_malloc(sizeof(unsigned int));

    *tmp = 0;
    *total = 0;

    MPI_Barrier(MPI_COMM_WORLD);

    for (i = 0 ; i < niters ; ++i) {
        smb_msgrate_cache_invalidate(cache_buf, cache_size);

        MPI_Barrier(MPI_COMM_WORLD);

        ar_timer_reset();
        nreqs = 0;
        for (j = 0 ; j < npeers ; ++j) {
            for (k = 0 ; k < nmsgs ; ++k) {
                MPI_Irecv(recv_buf + (nbytes * (k + j * nmsgs)),
                          nbytes, MPI_CHAR, recv_peers[j], 
                          magic_tag, 
                          //magic_tag + k + i * nmsgs, 
                          //MPI_ANY_TAG,
                          MPI_COMM_WORLD, &reqs[nreqs++]);
            }
            for (k = 0 ; k < nmsgs ; ++k) {
                MPI_Isend(send_buf + (nbytes * (k + j * nmsgs)),
                          nbytes, MPI_CHAR, send_peers[npeers - j - 1], 
                          magic_tag, 
                          //magic_tag + k + i * nmsgs, 
                          MPI_COMM_WORLD, &reqs[nreqs++]);
            }
        }
        MPI_Waitall(nreqs, reqs, statuses);
        *total += ar_timer_get_cycles();
    }

    //MPI_Allreduce(&total, &tmp, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    //smb_msgrate_display_result(rank, machine_output, "all-start", (niters * npeers * nmsgs * 2) / (tmp / world_size));

    MPI_Reduce(total, tmp, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    if (0 == rank) {
      rate = (float) (niters * npeers * nmsgs * 2) / (float) (*tmp / world_size / 1000);
      smb_msgrate_display_result(rank, machine_output, "all-start", rate);
    }


    kt_free(total);
    kt_free(tmp);

}


// Usage: msgrate -n <ppn> [OPTION]...
//   -h           Display this help message and exit
//   -p <num>     Number of peers used in communication
//   -i <num>     Number of iterations per test
//   -m <num>     Number of messages per peer per iteration
//   -s <size>    Number of bytes per message
//   -c <size>    Cache size in integers
//   -n <ppn>     Number of procs per node
//   -o           Format output to be machine readable
// \nReport bugs to <bwbarre@sandia.gov>

// Defaults:
// npeers = 6
// niters = 4096
// nmsgs = 128
// nbytes = 8
// cache_size = (8 * 1024 * 1024 / sizeof(int))
// ppn = -1
// machine_output = 0
int smb_msgrate_mpi(int npeers, int niters, int nmsgs, int nbytes,
                    int cache_size, int ppn, int machine_output) {

    int *send_peers;
    int *recv_peers;
    int *cache_buf;
    char *send_buf;
    char *recv_buf;
    MPI_Request *reqs;
    MPI_Status *statuses;

    int rank = -1;
    int world_size = -1;

    int i;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    /* sanity check */
    if (world_size < 3) {
        kt_printf( "Error: At least three processes are required\n");
        return 1;
    } else if (world_size <= npeers) {
        kt_printf( "Error: job size (%d) <= number of peers (%d)\n",
                world_size, npeers);
        return 1;
    } else if (ppn < 1) {
        kt_printf( "Error: must specify process per node (-n #)\n");
        return 1;
    } else if (world_size / ppn <= npeers) {
        kt_printf( "Error: node count <= number of peers\n");
        return 1;
    }


    if (0 == rank) {
        if (!machine_output) {
            kt_printf("job size:   %d\n", world_size);
            kt_printf("npeers:     %d\n", npeers);
            kt_printf("niters:     %d\n", niters);
            kt_printf("nmsgs:      %d\n", nmsgs);
            kt_printf("nbytes:     %d\n", nbytes);
            //kt_printf("cache size: %d\n", cache_size * (int)sizeof(int));
            //kt_printf("ppn:        %d\n", ppn);
        } else {
            kt_printf("%d %d %d %d %d %d %d ", 
                   world_size, npeers, niters, nmsgs, nbytes,
                   cache_size * (int)sizeof(int), ppn);
        }
    }

    /* allocate buffers */
    send_peers = kt_malloc(sizeof(int) * npeers);
    recv_peers = kt_malloc(sizeof(int) * npeers);
    cache_buf = kt_malloc(sizeof(int) * cache_size);
    send_buf = kt_malloc(npeers * nmsgs * nbytes);
    recv_buf = kt_malloc(npeers * nmsgs * nbytes);
    reqs = kt_malloc(sizeof(MPI_Request) * 2 * nmsgs * npeers);
    statuses = kt_malloc(sizeof(MPI_Status) * 2 * nmsgs * npeers);

    /* calculate peers */
    for (i = 0 ; i < npeers ; ++i) {
        if (i < npeers / 2) {
            send_peers[i] = (rank + world_size + ((i - npeers / 2) * ppn)) % world_size;
        } else {
            send_peers[i] = (rank + world_size + ((i - npeers / 2 + 1) * ppn)) % world_size;
        }
    }
    if (npeers % 2 == 0) {
        /* even */
        for (i = 0 ; i < npeers ; ++i) {
            if (i < npeers / 2) {
                recv_peers[i] = (rank + world_size + ((i - npeers / 2) *ppn)) % world_size;
            } else {
                recv_peers[i] = (rank + world_size + ((i - npeers / 2 + 1) * ppn)) % world_size;
            }
        } 
    } else {
        /* odd */
        for (i = 0 ; i < npeers ; ++i) {
            if (i < npeers / 2 + 1) {
                recv_peers[i] = (rank + world_size + ((i - npeers / 2 - 1) * ppn)) % world_size;
            } else {
                recv_peers[i] = (rank + world_size + ((i - npeers / 2) * ppn)) % world_size;
            }
        }
    }

    /* BWB: FIX ME: trash the free lists / malloc here */

    /* sync, although tests will do this on their own (in theory) */
    MPI_Barrier(MPI_COMM_WORLD);

    /* run tests */
#if 1
    smb_msgrate_test_one_way(world_size, rank, niters, cache_buf, cache_size, 
                             nmsgs, send_buf, nbytes, reqs, statuses, recv_buf, 
                             machine_output);
#endif

#if 1
    smb_msgrate_test_same_direction(niters, cache_buf, cache_size, npeers, 
                                    nmsgs, recv_buf, nbytes, recv_peers, reqs, 
                                    statuses, send_buf, send_peers, rank,
                                    machine_output, world_size);
#endif

#if 1
    smb_msgrate_test_prepost(npeers, nmsgs, recv_buf, nbytes, recv_peers, reqs,
                             niters, cache_buf, cache_size, send_buf, 
                             send_peers, statuses, rank, machine_output, 
                             world_size);
#endif

#if 1
    smb_msgrate_test_allstart(niters, cache_buf, cache_size, npeers, nmsgs, 
                              recv_buf, nbytes, recv_peers, reqs, send_buf, 
                              send_peers, statuses, rank, machine_output, 
                              world_size);
#endif

    if (rank == 0) kt_printf("\n");

    /* done */
    kt_free(send_peers);
    kt_free(recv_peers);
    kt_free(cache_buf);
    kt_free(send_buf);
    kt_free(recv_buf);
    kt_free(reqs);
    kt_free(statuses);
   
    MPI_Barrier(MPI_COMM_WORLD);

    return 0;
}
