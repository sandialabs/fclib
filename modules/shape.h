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
 * \file shape.h
 * \brief Declarations for \ref Shape module.
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/modules/shape.h,v $
 * $Revision: 1.52 $
 * $Date: 2006/09/28 21:05:27 $
 */

#ifndef _FC_SHAPE_H_
#define _FC_SHAPE_H_

#ifdef __cplusplus
extern "C" {
#endif

/**
 * \ingroup Shape
 * \brief Shape Struct
 *
 * \description
 *
 *  This struct holds information about shapes. Currently
 *  the info is about the sides of the shape only, which is
 *  skin information. If you want to keep the elements of the
 *  shape itself, then segment the
 *  mesh first to get that information and then call
 *  \ref fc_getSubsetShapes on each segment separately.
 *
 *  Note that the order of the struct arrays and matrix matters, that is
 *  The 0th element of the face array is the same side as the
 *  0th element of the elem array and is the same side as that
 *  referred to by the 0th row and col of the adjcency matrix.
 *
 *
 * \todo
 *  - may want to have a struct that does that the shape innards as well
 *    (e.g., I use the innards in projecting the dead element region
 *     onto a side)
 *  - may want a method to do the segment first as well so the
 *    user doesnt have to call that separately, or add a method
 *    that can generate the shape elements for you afterwards
 *    once you have the shape.
 *  - may want to add edge info into the struct
 *
 * \modifications
 *  - 6/XX/06 ACG created.
 *  - 7/16/06 ACG dropped "skin" in name
 */
typedef struct {
  int numSides;         /**< number of sides of the shape */
  FC_Subset *faces;     /**< array of faces making up the sides */
  FC_Subset *elems;     /**< array of elems making up the sides */
  int **adjmatrix;      /**< adjacency matrix or the sides */
} FC_Shape;


//***********elementary methods***********************/
//get sides

FC_ReturnCode fc_getSubsetShapes(FC_Subset subset, double angle,
			      int shared_dim, int *numShapes,
			      FC_Shape **shapes);
FC_ReturnCode fc_getMeshShapes(FC_Mesh mesh, double angle,
			      int shared_dim, int *numShapes,
			      FC_Shape **shapes);
FC_ReturnCode fc_reshapeByAngle(FC_Shape *inshape, double newangle, int shared_dim,
				  FC_Shape *newshape);
FC_ReturnCode fc_reduceShapeNumSides(FC_Shape *inshape, int endsides,
				    FC_Shape *reducedshape);

//adjacencies
FC_ReturnCode fc_createSideAdjacencyMatrix(int numSubsets,
					   FC_Subset *facesubsets,
					   int shared_dim,
					   int ***intmatrix);

FC_ReturnCode fc_segmentAdjacencyMatrix(int matrixdim, int **adjmatrix,
					int *numsegments, int**sidespersegment,
					int ***sideids);

FC_ReturnCode fc_calcDistanceMatrix(int nsides, int**adjmatrix,
				    int ***distmatrix);


//orders
FC_ReturnCode fc_createScrewShapeOrder(FC_Shape *shape, int **neworder);
FC_ReturnCode fc_createDescendingAreaOrder(FC_Shape *shape, int **neworder);
FC_ReturnCode fc_createLargeAndNonAdjacentSidesOrder(FC_Shape *shape,
						  int **neworder);
FC_ReturnCode fc_createLargeAndOpposingSidesOrder(FC_Shape *shape, 
					       double maxangle,
					       int **neworder);



//normal
FC_ReturnCode fc_getPlanarSideNormal(FC_Shape *shape, int sideid,
				       FC_Coords *normal);
FC_ReturnCode fc_getSideNormalMeanSdev(FC_Shape *shape,
				       int sideid,
				       FC_Vector *normal,
				       FC_Vector *std);

//misc 
FC_ReturnCode fc_getShapeEnds(FC_Shape *shape,
			      int *numEnds,
			      int **endarray);
FC_ReturnCode fc_getShapeSidesAreas(FC_Shape *shape, double **areas);


//helpers
FC_ReturnCode fc_copyShape(FC_Shape *inshape, FC_Shape *outshape);
FC_ReturnCode fc_freeShape(FC_Shape *fcss); 
FC_ReturnCode fc_clearShape(FC_Shape *fcss); 
FC_ReturnCode fc_initShape(FC_Shape *fcss); 
FC_ReturnCode fc_printShape(FC_Shape *fcss, char* label, int verb);
FC_ReturnCode fc_reorderShape(int *neworder, FC_Shape *shape);




#ifdef __cplusplus
}
#endif

#endif // _FC_SHAPE_H_

