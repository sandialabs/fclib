/*
 * Copyright (2000) Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
 * certain rights in this software.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *
 *   * Redistributions in binary form must reproduce the above copyright
 *     notice, this list of conditions and the following disclaimer in the
 *     documentation and/or other materials provided with the
 *     distribution.
 *
 *   * Neither the name of Sandia nor the names of any contributors may
 *     be used to endorse or promote products derived from this software
 *     without specific prior written permission.
 *
 *   * Modified source versions must be plainly marked as such, and must
 *     not be misrepresented as being the original software.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHORS OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
 * DAMAGE.
 */

/** 
 * \file table.c
 * \brief Implementation for \ref PrivateTable module.
 *
 * $Source: /home/Repositories/fcdmf/fclib/modules/table.c,v $
 * $Revision: 1.124 $ 
 * $Date: 2006/11/09 22:34:35 $
 */

// C library includes
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// fc library dependencies
#include "base.h"
#include "storage.h"
#include "library.h"

// this module
#include "tableP.h"

/** 
 * \addtogroup PrivateTable
 * \brief Routines for manipulating the table data structures.
 *
 * \description
 *  
 *   The primary data structures for FCLib are managed as "tables" (arrays) of
 *   "slots" (structs), which are structures that hold metadata and data 
 *   corresponding to that data entity type (i.e. _FC_DsSlot, _FC_SeqSlot, 
 *   _FC_MeshSlot, _FC_SubSlot, _FC_VarSlot). For more
 *   information about what information each type of slot holds, see the
 *   documentation on the struct definitions.
 *
 *   All slot types have a common header (_FC_SlotHeader) as the
 *   their first member. This allows us to 
 *   make some general purpose utility routines which work on any slot type.
 *   The header also contains information about the state of a slot like
 *   whether it can be modified on disk (writable), and whether the disk 
 *   version matches the current version (committed).
 *
 *   Any table consists of an array of pointers to the appropriate slot
 *   type. The tables are expanded as necessary by doubling the number of
 *   entries at a time. Entries start with a NULL value and are allocated as
 *   needed. The slots of "Deleted" slots are freed and set back to NULL.
 *
 *   Each table has supporting data structures - xTableSize is the length
 *   of the table; xOpenSlots is a sorted array of available slots.
 *
 *   Generic versions of basic table operations are defined in this module
 *   (like get a new slot or delete a slot). These basic operations take
 *   a specific kind of table and the necessary supporting information.
 *   To simplify matters, each kind of data (dataset, mesh, etc) has
 *   a wrapper that does the operation on that type of table (e.g. 
 *   _fc_getNewDsSlot() does _fc_getNewSlot() with the dataset table).
 *
 *   The data tables are not exposed at the user level and are 
 *   instead manipulated by \ref Handles which contain just enough information
 *   for called routines to locate the appropriate slot in the tables
 *   (i.e. the index into the table and a unique sequence number to make sure 
 *   the slot still refers to the same entity).
 */

/**
 * \ingroup PrivateTable
 * \brief Master unique ID
 *
 * \description
 *
 *   Each slot has a unique ID which is checked against
 *   when dereferencing a handle. Each uID is only used
 *   once, so this master uID keeps track of the next
 *   available number. uIDs start at 1 and go up. (Start
 *   with 1 so anything initialized with 0 won't be accidentally o.k.
 */
static int next_uID = 1;   

/**
 * \ingroup  PrivateTable
 * \brief  Replaces the name for a slot header with a copy of the input
 *
 * \description
 *  
 *    Takes a name, allocates space for the name (realloc if
 *    a name is already set) and copies it (user is still
 *    responsible for input char string).
 *
 *    This replaces the current name (and frees that memory). 
 *    Make sure to call this only if header has been initialized
 *    at some point with header.name = NULL.
 *
 * \modifications
 *    - Nancy Collins Created.
 */
FC_ReturnCode _fc_setSlotHeaderName(
  _FC_SlotHeader *header,   /**< input/output - pointer to a slot header */
  char *name                /**< input - name to set */
) {
  // check input
  if (!header || !name) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // delete the current name
  if (header->name)
    free(header->name);
  
  // copy the name
  header->name = malloc(strlen(name) + 1);
  if (header->name == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  strcpy(header->name, name);
  
  return FC_SUCCESS;
}

/**
 * \ingroup  PrivateTable
 * \brief  Initialize a table header
 *
 * \description
 *  
 *    Sets all table header values to appropriate default values.
 *    Except the slot number and the unique ID which should only be set by 
 *    _fc_new_xxxSlot.
 *
 * \modifications  
 *   APR-29-2003  W Koegler  Created. Moved code from getnextslot
 *
 */
void _fc_initSlotHeader(
  _FC_SlotHeader *header  /**< input - the slotheader to initialize */
) {
  if (header != NULL) {
    // header->slot is NOT changed
    // header->uID is NOT changed 
    header->name = NULL;
    header->committed = 0;
  }   
}

/**
 * \ingroup  PrivateTable
 * \brief  Clear a slot header's
 *
 * \description
 *  
 *    Releases all dynamically allocated resources in a 
 *    table header, and reinitializes all members.
 *
 * \modifications  
 *   - APR-29-2003  W Koegler  Created. 
 *   - 2003-APR-31  W Koegler  Made more comprehensive
 */
void _fc_clearSlotHeader(
  _FC_SlotHeader *header  /**< input/output - the slotheader to clear */
) {
  if (header != NULL) {
    // free dynamically allocated arrays
    if (header->name)
      free(header->name);
    
    // reinitialize
    _fc_initSlotHeader(header);
    
    // make "empty"
    header->uID = -1;
    // do NOT change header->slotID: for as long as the slot exists
    // it's slot ID is set and should stay the same
  }   
}

/**
 * \ingroup  PrivateTable
 * \brief  Print the contents of an _FC_SlotHeader to stderr.
 *
 * \description
 *  
 *    This is a helper routine for the _fc_printXTable routines, used
 *    to print _FC_SlotHeader information (all table have an _FC_SlotHeader).
 *    It is meant to be called from these routines and wont' have pretty
 *    formatting if called on its own; but the info will be right.
 *
 * \modifications  
 *    - Nancy Collins Created.
 *    - 2003-NOV-13  WSK  Fixed up.
 *    - 12/03/03 WSK Fixed up some more.
 */
void _fc_printSlotHeader(
  _FC_SlotHeader header 
) {
  fprintf(stderr, "     name = '%s'\n", header.name);
  fprintf(stderr, "     slotID = %d, uID = %d, committed = %d\n", 
          header.slotID, header.uID, header.committed);
  fflush(NULL);
}

/**
 * \ingroup  PrivateTable
 * \brief   Get an empty slot from table.
 *
 * \description
 *  
 *    Code to get an empty slot from any type of table. It first checks
 *    for an unused spot (either the entry is NULL or has a unique
 *    ID less than 0). If there are no unused slots, it grows the table
 *    by factor of 2. If the unused slot is NULL, it allocates a
 *    new slot.
 *
 * \modifications 
 *    - 5/11/04 WSK  Moved common code from _fc_getNewDsSlot(), 
 *      _fc_getNewSeqSlot(), _fc_getNewMeshSlot() & _fc_getNewVarSlot() to
 *      here.
 */
_FC_SlotHeader* _fc_getNewSlot(
   int* tableSize,           /**< input/output: pointer to table size  */
   _FC_SlotHeader*** table,  /**< input/output: pointer to table */
   FC_SortedIntArray* openSlots,     /**< input/output:  sia holds
				indices of available slots */
   int slotSize              /**< input: size of a slot */
) {
  int i;
  int slotID, newSize, ret;
  void *temp;

  // get an unused slot
  if (openSlots->numVal > 0) {
    slotID = openSlots->vals[0];
    fc_popSortedIntArrayFront(openSlots);
  }
  else { // if an unused slot doesn't exist, make some more slots
    if (*tableSize == 0)
      newSize = 1;
    else
      newSize = 2*(*tableSize);
    temp = realloc(*table, newSize*sizeof(_FC_SlotHeader*));
    if (temp == NULL) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      return NULL;
    }
    *table = temp;
    for (i = (*tableSize); i < newSize; i++)
      (*table)[i] = NULL;
    // NOTE! the ID of the first new slot doesn't get added to openSlots
    slotID = (*tableSize);
    for (i = newSize-1; i > (*tableSize); i--) {
      ret = fc_addIntToSortedIntArray(openSlots, i);
      if (ret < 0) {
	fc_printfErrorMessage("Failed to add key to sia, memory error?");
	return NULL; // an error occured
      }
      else if (ret == 0) {
	fprintf(stderr, "***developer note: Algorithmic error! shouldn't reach here!");
	fflush(NULL);
      }
    }
    (*tableSize) = newSize;
  }

  // allocate slot if necessary
  if ((*table)[slotID] == NULL) {
    (*table)[slotID] = malloc(slotSize);
    if ((*table)[slotID] == NULL) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      return NULL;
    }
  }
        
  // setup the slot
  _fc_initSlotHeader((*table)[slotID]);
  (*table)[slotID]->slotID = slotID;
  (*table)[slotID]->uID = next_uID;
  next_uID++;

  return (*table)[slotID];
}

/**
 * \ingroup  PrivateTable
 * \brief   Check that a handle is valid.
 *
 * \description 
 * 
 *    Returns 1 if the handle is valid and 0 if it is not. Valid means that the
 *    object exists and that it has the same uID as the handle (this extra
 *    check is necessary because the slot may get reused, at which time the
 *    handle is invalid even though it refers to a real object).
 *
 * \modifications 
 *    - 05/25/04 WSK Created.
 */
int _fc_isHandleValid(
   int slotID,             /**< input: slot id from a handle  */
   int uID,                /**< input: unique id from a handle */
   int tableSize,          /**< input: the table size  */
   _FC_SlotHeader** table  /**< input: the table */
) {
  if (slotID < 0 || slotID >= tableSize || !table[slotID] ||
      table[slotID]->uID != uID)
    return 0;
  else
    return 1;
}

/**
 * \ingroup  PrivateTable
 * \brief   Get a slot by ID.
 *
 * \description  
 *
 *    Returns a slot if it exists (exists means has some thing in it and
 *    therefore the uID will be 1 or greater. Otherwise returns NULL.
 *
 * \modifications 
 *    - 05/25/04 WSK Created.
 */
_FC_SlotHeader* _fc_getSlot(
   int slotID,             /**< input: slot id  */
   int tableSize,          /**< input: the table size  */
   _FC_SlotHeader** table  /**< input: the table */
) {
  if (slotID < 0 || slotID >= tableSize || !table[slotID] ||
      table[slotID]->uID < 1)
    return NULL;
  else
    return table[slotID];
}

/**
 * \ingroup  PrivateTable
 * \brief   Delete the contents of a slot
 *
 * \description  
 *
 *    If the slot exists, deletes the memory associated with
 *    that slot and add the slotID to the openSlots lists.
 *
 *    Assumes that dynamic data in the slot has already been
 *    "cleared".
 *
 * \modifications 
 *    - 12/20/05 WSD Created.
 */
FC_ReturnCode _fc_deleteSlot(
   int slotID,             /**< input: slot id  */
   int tableSize,          /**< input: the table size  */
   _FC_SlotHeader** table, /**< input: the table */
   FC_SortedIntArray* openSlots  /**< input/output: the open slots */
) {
  int ret;

  if (slotID >= 0 && slotID < tableSize && table[slotID]) {
    free(table[slotID]);
    table[slotID] = NULL;
    ret = fc_addIntToSortedIntArray(openSlots, slotID);
    if (ret < 0) {
      fc_printfErrorMessage("Failed to add ID to sia");
      return ret; // an error occured
    }
    else if (ret == 0) {
      fprintf(stderr, "***developer note: Algorithmic error! shouldn't reach here!");
      fflush(NULL);
    }
  }
  
  return FC_SUCCESS;
}
