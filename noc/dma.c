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
// Author        : Spyros LYBERIS
// Abstract      : Network-on-chip DMA transfer functions
//
// =============================[ CVS Variables ]=============================
// 
// File name     : $RCSfile: dma.c,v $
// CVS revision  : $Revision: 1.8 $
// Last modified : $Date: 2013/02/15 10:21:36 $
// Last author   : $Author: lyberis-spree $
// 
// ===========================================================================

#include <kernel_toolset.h>
#include <arch.h>
#include <memory_management.h>
#include <processing.h>


// ===========================================================================
// noc_dma_add_new()            Creates a new DMA structure, appends it to
//                              the list of active DMAs and tries to start
//                              it if there are available hw counters.
// ===========================================================================
// * INPUTS
//   int src_core_id            Unique core ID of the source core
//   void *src_adr              Source buffer address
//   int dst_core_id            Unique core ID of the destination core, or
//                              -1 to indicate the DRAM memory of the board
//                              where src_core_id belongs
//   void *dst_adr              Destination buffer address
//   int size                   Transfer size in bytes
//   int i_flag                 DMA engine I flag (ignore dirty on src)
//   int w_flag                 DMA engine W flag (write-through on dst)
//   int c_flag                 DMA engine C flag (clean on dst)
//   int *group_notif_var       If not NULL, the pointed variable will be
//                              decremented by 1 when the DMA has been
//                              successfully finished. This can be used to
//                              group multiple DMAs together, so that the
//                              completion of a DMA group can be monitored. The
//                              same group_notif_var must be given for each new
//                              DMA of a group and the initial value should be
//                              set by the caller to the number of DMAs that
//                              will be requested. Can be NULL, if this
//                              mechanism is not required.
//   void *group_notif_data     Generic data hook, not touched by the DMA
//                              functions. Will be returned to the caller
//                              when a DMA group finishes. Because the order
//                              of the DMA completions may be arbitrary, it 
//                              makes sense to give the same group_notif_data
//                              for all new DMAs of a group, so when the last
//                              of the group completes the data hook is valid.
//                              Can be NULL, if this mechanism is not required.
//
// * RETURN VALUE
//   int                        0: success and DMA was started
//                              1: success, but DMA was queued for later on
// ===========================================================================
int noc_dma_add_new(int src_core_id, void *src_adr, int dst_core_id, 
                    void *dst_adr, int size, int i_flag, int w_flag, 
                    int c_flag,int *group_notif_var, void *group_notif_data) {

  int           my_bid;
  int           my_cid;
  Context       *context;
  NocDma        *dma;
  ListNode      *cnt_node;


  // Get context and architectural board/core ID
  my_bid = ar_get_board_id();
  my_cid = ar_get_core_id();
  context = mm_get_context(my_cid);
  ar_assert(size);

  // Create a new NocDma structure and initialize it
  dma = kt_malloc(sizeof(NocDma));
  pr_core_arch_bid_cid(src_core_id, &(dma->src_bid), &(dma->src_cid));
  dma->src_adr = (unsigned int) src_adr;
  if (dst_core_id == -1) {
    dma->dst_bid = dma->src_bid;
    dma->dst_cid = 0xC; // DRAM destination
  }
  else {
    pr_core_arch_bid_cid(dst_core_id, &(dma->dst_bid), &(dma->dst_cid));
  }
  dma->dst_adr = (unsigned int) dst_adr;
  dma->size = size;
  dma->flags[0] = i_flag;
  dma->flags[1] = w_flag;
  dma->flags[2] = c_flag;
  dma->loc_notif_var = group_notif_var;
  dma->loc_notif_data = group_notif_data;

  // Add the DMA to the tail of the active DMAs list
  kt_list_insert(context->noc_active_dmas, dma, 
                 kt_list_tail(context->noc_active_dmas), 1);


  // Try to allocate a hardware counter
  cnt_node = kt_list_head(context->noc_cnt_free);
  if (cnt_node) {
    dma->local_cnt = (int) cnt_node->data;
    kt_list_delete(context->noc_cnt_free, cnt_node, NULL);
  }
  else {
    // Mark no counter was available and quit
    dma->local_cnt = -1;
    return 1;
  }


  // Wait until our DMA engine can support at least one more DMA
  while (!(ar_ni_status_get(my_cid) & 0xFF)) {
    ;
  }

  // Start the DMA
//ar_assert(dma->src_adr == dma->dst_adr); kt_printf("%X %X/%d->%X/%d %d%d%d\r\n", dma->src_adr, dma->src_bid, dma->src_cid, dma->dst_bid, dma->dst_cid, dma->flags[0], dma->flags[2], dma->flags[1]);
  ar_assert((dma->local_cnt >= 0) && (dma->local_cnt < NOC_MAX_COUNTERS));
  ar_cnt_set(my_cid, dma->local_cnt, -dma->size);
  ar_dma_with_ack(my_cid,
                  dma->src_bid, dma->src_cid, dma->src_adr,
                  dma->dst_bid, dma->dst_cid, dma->dst_adr,
                  my_bid,       my_cid,       dma->local_cnt,
                  dma->size, dma->flags[0], dma->flags[2], dma->flags[1]);

  // Success
  return 0;

}


// ===========================================================================
// noc_dma_check_progress()     Checks the active DMA list for any progress.
//                              The function tries to allocate hw counters
//                              for DMAs that don't have one. Failed DMAs
//                              that have got a Nack response are restarted.
//                              Completed DMAs are removed from the list and
//                              their group notification variables are
//                              decremented. 
// ===========================================================================
// * OUTPUTS
//   void ***ret_group_data     If not NULL, an array of all completed group
//                              data hooks is allocated and returned here.
//                              The array is allocated on *ret_group_data
//                              (type is void **) and must be freed by the
//                              caller.
//
// * RETURN VALUE
//   int                        Number of groups that all of their DMAs were 
//                              completed successfully (and size of the 
//                              *ret_group_data allocated array)
// ===========================================================================
int noc_dma_check_progress(void ***ret_group_data) {

  int           my_bid;
  int           my_cid;
  Context       *context;
  NocDma        *dma;
  ListNode      *cnt_node;
  ListNode      *n;
  ListNode      *next_n;
  int           collect;
  int           restart;
  int           groups_completed;


  // Get context and architectural board/core ID
  my_bid = ar_get_board_id();
  my_cid = ar_get_core_id();
  context = mm_get_context(my_cid);

  // Remember how many groups of DMAs we have completed; no array yet
  groups_completed = 0;
  if (ret_group_data) {
    *ret_group_data = NULL;
  }

  // For all active DMAs
  n = kt_list_head(context->noc_active_dmas);
  while (n) {

    // Discover the NocDma structure
    dma = n->data;
    ar_assert(dma);

    // No actions by default
    collect = 0;
    restart = 0;

    // Has the DMA been started in the past? If yes, check its current status.
    if (dma->local_cnt > -1) {
      switch (ar_cnt_get_triggered(my_cid, dma->local_cnt)) {
        case 0:
          // Still in progress, do nothing
          break;
        case 2:
          // Finished normally, collect this DMA
          collect = 1;
          ar_assert(!ar_cnt_get(my_cid, dma->local_cnt)); // sanity check
          break;
        case 3:
          // Got a Nack, restart the DMA
          restart = 1;
          break;
        default:
          // Unknown response
          ar_abort();
      }
    }

    // If the DMA doesn't have a hw counter, try to allocate one now
    else {
      cnt_node = kt_list_head(context->noc_cnt_free);
      if (cnt_node) {
        // Success, remember the counter and mark to start the DMA
        dma->local_cnt = (int) cnt_node->data;
        kt_list_delete(context->noc_cnt_free, cnt_node, NULL);
        restart = 1;
      }
//else{kt_printf("%d: OUT OF COUNTERS\r\n", context->pr_core_id);}
    }

    // Should we (re)start this DMA?
    if (restart) {

      // Wait until our DMA engine can support at least one more DMA
      while (!(ar_ni_status_get(my_cid) & 0xFF)) {
        ;
      }

      // Start the DMA
//ar_assert(dma->src_adr == dma->dst_adr); kt_printf("***RESTART %X %X/%d->%X/%d %d%d%d\r\n", dma->src_adr, dma->src_bid, dma->src_cid, dma->dst_bid, dma->dst_cid, dma->flags[0], dma->flags[2], dma->flags[1]);

      ar_assert((dma->local_cnt >= 0) && (dma->local_cnt < NOC_MAX_COUNTERS));
      ar_cnt_set(my_cid, dma->local_cnt, -dma->size);
      ar_dma_with_ack(my_cid,
                      dma->src_bid, dma->src_cid, dma->src_adr,
                      dma->dst_bid, dma->dst_cid, dma->dst_adr,
                      my_bid,       my_cid,       dma->local_cnt,
                      dma->size, dma->flags[0], dma->flags[2], dma->flags[1]);
    }

    // Remember next node in the DMA list, in case we delete the current one
    next_n = kt_list_next(n);

    // Should we delete this DMA?
    if (collect) {

      // Release the hw counter back onto the free list
      kt_list_insert(context->noc_cnt_free, (void *) dma->local_cnt,
                     kt_list_tail(context->noc_cnt_free), 1);
//kt_printf("%d: dma collect 0x%08X-- (now %d)\r\n", context->pr_core_id, dma->loc_notif_var, (dma->loc_notif_var) ? *dma->loc_notif_var - 1 : -1);
      // If there's a notification variable, notify it by decrementing by one
      if (dma->loc_notif_var) {
        ar_assert(*dma->loc_notif_var > 0);
        (*dma->loc_notif_var)--;

        // If the notification variable reached 0, one more DMA group is done.
        if (!(*dma->loc_notif_var)) {

          // Remember how many, and (re)allocate array and store group data
          // hook. We only do this if the user has supplied a data hook.
          if (ret_group_data) {
            *ret_group_data = kt_realloc(*ret_group_data,
                                  (groups_completed + 1) * sizeof (void *));
            (*ret_group_data)[groups_completed] = dma->loc_notif_data;
          }

          groups_completed++;
        }

      }

      // Remove the DMA from the list and free the data structure
      kt_list_delete(context->noc_active_dmas, n, kt_free);
    }

    // Go to the next DMA in the list
    n = next_n;
  }

  // Return how many DMA groups we've managed to complete
  return groups_completed;
}
