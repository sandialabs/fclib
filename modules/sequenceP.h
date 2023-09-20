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
 * \file sequenceP.h
 * \brief Private declarations for the \ref Sequence module.
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/modules/sequenceP.h,v $
 * $Revision: 1.14 $ 
 * $Date: 2006/08/30 19:20:01 $
 */

#ifndef _FC_SEQUENCE_P_H_
#define _FC_SEQUENCE_P_H_

#ifdef __cplusplus
extern "C" {
#endif

// include these in basic data structures to simplify includes for the
// the rest of the library which won't ever use them
#include "tableP.h"
#include "fileio.h"
#include "fileioP.h"

/**
 * \ingroup  PrivateSequence
 *
 * \brief Sequence Slot Structure for the Sequence Table
 *
 * \description 
 *
 *   This structure contains meta data and big data for a sequence. A sequence
 *   is an ordered series which is orthogonal to the spatial dimensions of the
 *   meshes. The most common type of sequence is a time series where the
 *   coordinates of the series are the time values. Besides the information
 *   stored in the header (see _FC_SlotHeader), this struct has a handle for
 *   the owning dataset and a structure containing file IO info. It also stores
 *   the coordinates of the sequence which are assumed to be 1D.
 */
typedef struct {
  _FC_SlotHeader header;  /**< Header, which contains general purpose info. */
  FC_Dataset ds;      /**< Handle to the owning dataset. */
  //---File reference
  _FC_SeqFileIOInfo fileInfo; /**< Holds info pertinent to file access. */
  //---About the sequence:  coordinates
  int numStep;        /**< The number of steps in this sequence (e.g. the 
			 number of time steps). */
  FC_DataType datatype;  /**< The data type of the sequence coords (e.g. int, 
                            float, etc.). */
  void *coords;       /**< Big Data: the coordinate values (assumed 1D). */
} _FC_SeqSlot;

// slot methods
void _fc_initSeqSlot(_FC_SeqSlot* seqSlot);
void _fc_releaseSeqSlot(_FC_SeqSlot *seqSlot);
void _fc_clearSeqSlot(_FC_SeqSlot *seqSlot);

// table access
int _fc_getSeqTableSize(void);
_FC_SeqSlot* _fc_getNewSeqSlot(void);
_FC_SeqSlot* _fc_getSeqSlot(FC_Sequence sequence);
_FC_SeqSlot* _fc_getSeqSlotFromID(int seqID);
FC_ReturnCode _fc_deleteSeqSlot(FC_Sequence sequence);
void _fc_printSeqTable(char *label);
void _fc_freeSeqTable(void);

#ifdef __cplusplus
}
#endif

#endif // _FC_SEQUENCE_P_H_

