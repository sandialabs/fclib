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
 * \file variableP.h
 * \brief Private declarations for \ref Variable module.
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/modules/variableP.h,v $
 * $Revision: 1.17 $ 
 * $Date: 2006/09/22 02:16:52 $
 */

#ifndef _FC_VARIABLE_P_H_
#define _FC_VARIABLE_P_H_

#ifdef __cplusplus
extern "C" {
#endif

// include these in basic data structures to simplify includes for the
// the rest of the library which won't ever use them
#include "tableP.h"
#include "fileio.h"
#include "fileioP.h"

/**
 * \ingroup  PrivateVariable
 * \brief Variable Slot Structure for the Variable Table
 *
 * \description
 * 
 *   This structure contains meta data and big data for a variable.
 *   First, it contains a header, a handle to the owning mesh,
 *   and a structure containing file IO info. In addition, it
 *   has the data for the variable.
 *
 * \todo? There is duplication between the mesh being NULL and
 *   and the assoc being FC_AT_DATASET. Is this fixable? 
 *   Maybe we don't need the FC_AT_DATASET version?
 */
typedef struct {
  _FC_SlotHeader header;  /**< Header, which contains general purpose info. */
  FC_Dataset ds;      /**< Handle to the owning dataset. */
  FC_Mesh mesh;       /**< Handle to the owning mesh. If this is a global
		       variable, this will be FC_NULL_MESH. */
  FC_Sequence sequence; /**< Handle to owning sequence (can be NULL) */
  int stepID;      /**< ID of step in a sequence, -1 if sequence != NULL */
  //---File reference
  _FC_VarFileIOInfo fileInfo; /**< Holds info pertinent to file access. */
   //---About the variable:  data
  int numDataPoint;   /**< The number of data points (should match the number 
			 of the mesh subentities of the data association type 
                         (e.g. should be the number vertices, or elements, 
                         etc.). */
  int numComponent;   /**< The number of components for each data point. E.g.
                          2 for a vector with x and y components or 9 for
                          a 3X3 tensor. */
  FC_AssociationType assoc;  /**< The subentities of the mesh that the data is 
                             associated with (e.g. vertices, elements, etc). */
  FC_MathType mathtype; /**< Mathematical ordering of the components of a 
                            a data point. E.g. scalar, vector, tensor, etc. */
  FC_DataType datatype;  /**< The data type of data (e.g. int, float, etc.) */
  void *data;   /**< Void pointer to the data which will be ordered 1st by 
		   the data point, then by the component (e.g. xyzxyz ...). */
} _FC_VarSlot;

// slot
void _fc_initVarSlot(_FC_VarSlot* varSlot);
void _fc_releaseVarSlot(_FC_VarSlot *varSlot);
void _fc_clearVarSlot(_FC_VarSlot *varSlot);

// table access
int _fc_getVarTableSize(void);
_FC_VarSlot* _fc_getVarSlot(FC_Variable variable);
_FC_VarSlot* _fc_getVarSlotFromID(int varID);
_FC_VarSlot* _fc_getNewVarSlot(void);
FC_ReturnCode _fc_deleteVarSlot(FC_Variable variable);
void _fc_printVarTable(char *label);
void  _fc_freeVarTable(void);

#ifdef __cplusplus
}
#endif

#endif // _FC_VARIABLE_P_H_

