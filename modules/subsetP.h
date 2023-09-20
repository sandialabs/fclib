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
 * \file subsetP.h
 * \brief Private declarations for the \ref Subset module.
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/modules/subsetP.h,v $
 * $Revision: 1.41 $ 
 * $Date: 2006/08/30 19:20:01 $
 */

#ifndef _FC_SUBSET_P_H_
#define _FC_SUBSET_P_H_

#ifdef __cplusplus
extern "C" {
#endif

// include these in basic data structures to simplify includes for the
// the rest of the library which won't ever use them
#include "tableP.h"
#include "fileio.h"
#include "fileioP.h"

/**
 * \ingroup  PrivateSubset
 * \brief Struct describing subset elements
 *
 * \description 
 * 
 *    Structure to hold information about subset members. 
 *
 * \modifications 
 *    - 03/25/04 RM Created.
 *    - 05/11/04 WSK Modified to become like Dataset, Mesh, etc.
 *    - 08/10/06 WSD Changed to have only 1 type of storage = sorted int array.
 */
typedef struct {
  _FC_SlotHeader header; /**< Header, which contains general purpose info. */
  FC_Mesh mesh;          /**< Handle to owning mesh */
  //---File reference
  _FC_SubFileIOInfo fileInfo; /**< Holds info pertinent to file access. */
  //---About the subset: members
  FC_AssociationType assoc; /**< subset association (elements, vortices etc.)*/
  int numMember;         /**< number of members, have to have in addition
			   to sia.numVal because sia is big data and could be
			   empty */
  int maxNumMember;      /**< maximum value for a member, same as the 
			    number of elements in mask array */
  FC_SortedIntArray sia; /**< member IDs stored in sorted int array */
  // int numMember;      // this is now sia.numVal
  // int* members;       // this is sia.vals 
} _FC_SubSlot;

// slot routines
void _fc_initSubSlot(_FC_SubSlot *subSlot);
void _fc_releaseSubSlot(_FC_SubSlot *subSlot);
void _fc_clearSubSlot(_FC_SubSlot *subSlot);

// table access
int _fc_getSubTableSize(void);
_FC_SubSlot* _fc_getSubSlot(FC_Subset subset);
_FC_SubSlot* _fc_getSubSlotFromID(int subID);
_FC_SubSlot* _fc_getNewSubSlot(void);
FC_ReturnCode _fc_deleteSubSlot(FC_Subset subset);
void _fc_printSubTable(char* label);
void _fc_freeSubTable(void);

#ifdef __cplusplus
}
#endif

#endif // _FC_SUBSET_P_H_
