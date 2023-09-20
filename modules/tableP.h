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
 * \file tableP.h
 * \brief Private declarations for \ref PrivateTable module.
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/modules/tableP.h,v $
 * $Revision: 1.52 $ 
 * $Date: 2006/08/30 19:20:01 $
 */

#ifndef _FC_TABLE_P_H_
#define _FC_TABLE_P_H_

#ifdef __cplusplus
extern "C" {
#endif

/**
 * \ingroup  PrivateTable
 * \brief Indent string for _fc_printXTable() routines.
 */
#define INDENT_STR "    "

/**
 * \ingroup PrivateTable
 * \brief A string of hyphens to use as a table divider for 
          _fc_printXTable() routines 
*/
#define TABLE_DIVIDER_STRING "----------------------------------------------\n"

/** 
 * \ingroup  PrivateTable
 * \brief  Get a handle from a table slot
 *
 * \description
 * 
 *    Copies values from the table slot into the handle 
 *    Warning! No type checking
 */
#define _FC_GET_HANDLE(handle /* output - the handle to be made */, \
                       xSlot_pointer /* input - a pointer to a table slot */) \
         ( (handle).slotID = (xSlot_pointer)->header.slotID, \
           (handle).uID = (xSlot_pointer)->header.uID )

/**
 * \ingroup  PrivateTable
 * \brief Slot Header
 *
 * \description
 * 
 *   This header is used as the first entry for all types of slots in
 *   the data tables. It contains basic identifying information including
 *   the slot ID (index into a particular data table), the sequence
 *   ID (an identifier which is unique over all types of data slots),
 *   a name string and some state flags. The state flags are described below.
 *   Note that the writable flag only applies to the disk representation of
 *   the entity, all in-core entities can be modified. 
 */
typedef struct {
  int slotID;    /**< This slot's index into its owning table. */
  int uID;       /**< Unique identifier for this entity. Since slots can be 
                     reused, this is needed to make sure that handle is 
                     still referring to the same entity. uIDs are unique
                     over all tables and all time (unless uID overflows). */
  char *name;    /**< Name for entity described by this slot. */
  //int writable;  /**< Flag indicating whether entity is writable (e.g. can be 
  //                 modified on disk). 
  //                 1 = writable, 0 = read-only (0 is the default). */
  int committed; /**< Flag indicating whether current data on an entity is
                     the same as on disk. 1 = committed, 0 = data in file and 
                     data in memory differ (0 is the default). */
} _FC_SlotHeader;

// slot header routines
FC_ReturnCode _fc_setSlotHeaderName(_FC_SlotHeader *header, char *name);
void _fc_initSlotHeader(_FC_SlotHeader *header);
void _fc_clearSlotHeader(_FC_SlotHeader *header);
void _fc_printSlotHeader(_FC_SlotHeader header);

// table routines
_FC_SlotHeader* _fc_getNewSlot(int* tableSize, _FC_SlotHeader*** table,
			       FC_SortedIntArray* openSlots, int slotSize);
_FC_SlotHeader* _fc_getSlot(int slotID, int tableSize, _FC_SlotHeader** table);
FC_ReturnCode _fc_deleteSlot(int slotID, int tableSize, _FC_SlotHeader** table,
			     FC_SortedIntArray* openSlots); 
int _fc_isHandleValid(int slotID, int uID, int tableSize, 
		      _FC_SlotHeader** table);


#ifdef __cplusplus
}
#endif

#endif // _FC_TABLE_P_H_

