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
 * \file datasetP.h
 * \brief Private declarations for the \ref Dataset module.
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/modules/datasetP.h,v $
 * $Revision: 1.19 $ 
 * $Date: 2006/09/22 02:16:52 $
 */

#ifndef _FC_DATASET_P_H_
#define _FC_DATASET_P_H_

#ifdef __cplusplus
extern "C" {
#endif

// include these in basic data structures to simplify includes for the
// the rest of the library which won't ever use them
#include "tableP.h"
#include "fileio.h"
#include "fileioP.h"

/**
 * \ingroup  PrivateDataset
 * \brief Dataset Slot Structure for the Dataset Table
 *
 * \description
 *
 *   This structure contains meta data for a dataset. Besides the information
 *   stored in the header (see _FC_SlotHeader), this struct has a structure
 *   containing file IO info. It also knows which sequences and meshes it owns.
 *
 * \todo ?someday make fileInfo dynamically allocated -- would only save
 *   a little space at the expensse of added complexity. Unless there
 *   is dynamically allocated stuff within fileInfo. Then should
 *   probably use special free routines.
 *
 * \modifications
 *    - 7/17/2005 WSD moved file IO stuff out into file io modules.
 */
typedef struct {
  _FC_SlotHeader header; /**<  Header, which contains general purpose info. */
  //---File reference
  FC_FileIOType fileType;  /**< Type of file format */
  char *baseFile;          /**< Path to primary file of format */
  _FC_DsFileIOInfo fileInfo; /**< Holds info pertinent to file access. */
  //---Flag for whether dataset is writable
  int writable;   /**< 1 = can write this dataset to disk. 0 = can NOT write 
		     dataset to disk (because it was loaded as read only) */
  //---Owned entities:
  //---sequences
  int numSeq;     /**< Number of sequences owned by this dataset. */
  int *seqIDs;    /**< Array of the slot numbers identifying the sequences
                      owned by this dataset. */
  //---meshes
  int numMesh;    /**< Number of meshes owned by this dataset. */
  int *meshIDs;   /**< Array of the slot numbers identifying the meshes
                      owned by this dataset. */
  //---global variables
  int numBasicVar; /**< Number of basic variables (i.e. non-sequence vars) */
  int* basicVarIDs; /** < Array of slot IDs for the basic variables */
  int numSeqVar;    /**< Number of sequence variables owned. */
  int *numStepPerSeqVar; /**< Number of steps in each seq var */
  int **seqVarIDs;   /**< Array of slot IDs identifying the vars that make 
                          up the sequence vars. This array has the dimensions
                          [numSeqVar][numStepPerSeq]. */
} _FC_DsSlot;

// basic slot 
void _fc_initDsSlot(_FC_DsSlot* dsSlot);
void _fc_clearDsSlot(_FC_DsSlot *dsSlot);

// table access
int _fc_getDsTableSize(void);
_FC_DsSlot* _fc_getNewDsSlot(void);
_FC_DsSlot* _fc_getDsSlot(FC_Dataset dataset);
_FC_DsSlot* _fc_getDsSlotFromID(int dsID);
FC_ReturnCode _fc_deleteDsSlot(FC_Dataset dataset);
void _fc_printDsTable(char *label);
void _fc_freeDsTable(void);

#ifdef __cplusplus
}
#endif

#endif // _FC_DATASET_P_H_

