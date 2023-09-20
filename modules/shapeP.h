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
 * \file shapeP.h
 * \brief Private declarations for \ref Shape module.
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/modules/shapeP.h,v $
 * $Revision: 1.14 $
 * $Date: 2006/08/30 19:20:01 $
 */

#ifndef _FC_SHAPE_P_H_
#define _FC_SHAPE_P_H_

#ifdef __cplusplus
extern "C" {
#endif

typedef struct _intpair{
    int ID1;
    int ID2;
}_FC_INTPAIR;

typedef struct _faceinfo{
  int faceID; /**< id of the face of a part of the skin. decide if want to
		 keep this since its also the key */
  int elemID; /**< id of the element giving rise to this face. there is only
		 1 since its a skin piece */
  FC_Vector normal; /**< oriented normal to the face */
}_FC_FACEINFO;

typedef struct _sfi_array{ //sorted array
  int numVals;
  int maxVals;               
  _FC_FACEINFO** vals;         /**< array of faceinfo ptrs*/ 
}_FC_SFI_ARRAY;


//Innards for building shapes
FC_ReturnCode _fc_getOutwardFaceNormal( FC_Mesh mesh, int faceID,
					int elemID, FC_Vector* normal);

FC_ReturnCode _fc_buildFaceInfoListFromSubset(FC_Subset subset,
					      _FC_SFI_ARRAY *sfi);

FC_ReturnCode _fc_buildFaceInfoListFromShape(FC_Shape *shape,
					     _FC_SFI_ARRAY *sfi);

FC_ReturnCode _fc_createShapePartsFromFaceInfoList( FC_Mesh mesh,
					    _FC_SFI_ARRAY *in_faceinfolist,
					    char **basename,            
					    double compangle,
					    int shared_dim,
					    int *numSides, 
					    FC_Subset** newFaceSubsets,
					    FC_Subset** newElemSubsets,
					    int ***adjmatrix);

//helpers for sorted face info array
FC_ReturnCode _fc_initSFI_Array(_FC_SFI_ARRAY *sfi);
int _fc_addToSFI_Array(_FC_SFI_ARRAY *sfi, _FC_FACEINFO* fi);
FC_ReturnCode _fc_lookupInSFI_Array(_FC_SFI_ARRAY *sfi, int i, _FC_FACEINFO** fi);
FC_ReturnCode _fc_removeFromSFI_Array(_FC_SFI_ARRAY *sfi, int i, _FC_FACEINFO** fi);
int _fc_reallocSFI_Array(_FC_SFI_ARRAY *sfi);
void _fc_printSFI_Array(_FC_SFI_ARRAY *sfi);
void _fc_freeSFI_Array(_FC_SFI_ARRAY *sfi);


#ifdef __cplusplus
}
#endif

#endif // _FC_SHAPE_P_H_
