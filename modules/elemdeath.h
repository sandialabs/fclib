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
 * \file elemdeath.h
 * \brief Declarations for \ref ElemDeath module.
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/modules/elemdeath.h,v $
 * $Revision: 1.30 $
 * $Date: 2006/08/30 19:20:01 $
 */

#ifndef _FC_ELEMDEATH_H_
#define _FC_ELEMDEATH_H_

#ifdef __cplusplus
extern "C" {
#endif

// Dead element regions
FC_ReturnCode fc_getExposedSkin(FC_Subset subset, FC_Subset *opt_meshSkin,
				FC_Subset* exposed);

FC_ReturnCode fc_calcTearLength(FC_Subset subset, FC_Subset *opt_meshSkin,
				FC_Variable* displ_coords, double* length);

//dead region segmenting mesh
FC_ReturnCode fc_subsetSegmentsMesh(FC_Subset subset, int shared_segdim,
   char* varname, int *numSubset, FC_Subset **newSubset);

FC_ReturnCode fc_subsetSegmentsSubset(FC_Subset subset_inner, 
				      FC_Subset subset_outer,
				      int shared_segdim,
				      char* varname, int *numSubset,
				      FC_Subset **newSubset);

FC_ReturnCode fc_subsetPlusNeighborGreaterSegmentationMesh(FC_Subset subset,
						   int shared_neighbordim,
						   int shared_segdim,
						 int* numNbr, int ** nBrIDs);

//decay of subsets
FC_ReturnCode fc_getDecayedSkin(FC_Subset subset, FC_Subset *opt_meshSkin,
				FC_Subset* decayed);
FC_ReturnCode fc_getSubsetDecayType(FC_Subset deadSubset, FC_Subset compSubset,
				 int* sidedecayflag);


 //decay of shapes
FC_ReturnCode fc_getDecayedShapeSides(FC_Subset subset, FC_Shape *shape,
				      FC_Subset **decayed);
FC_ReturnCode fc_getDecayedShapeSkin(FC_Subset subset, FC_Shape *shape,
				      FC_Subset *decayed);
FC_ReturnCode fc_getShapeSidesDecayType(FC_Subset deadSubset, FC_Shape *shape,
				 int** sidedecayflag);

#ifdef __cplusplus
}
#endif

#endif  // end of _FC_ELEMDEATH_H_
