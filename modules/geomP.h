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
 * \file geomP.h
 * \brief Private declarations for \ref GeometricRelations module.
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/modules/geomP.h,v $
 * $Revision: 1.36 $ 
 * $Date: 2007/07/06 04:58:34 $
 */

#ifndef _FC_GEOM_P_H_
#define _FC_GEOM_P_H_

#ifdef __cplusplus
extern "C" {
#endif


// Math helpers
double _fc_calcDeterminant2x2(double a11, double a12,
				double a21, double a22);
double _fc_calcDeterminant3x3(double a11, double a12, double a13,
			      double a21, double a22, double a23,
			      double a31, double a32, double a33);

// Diameter helpers
FC_ReturnCode _fc_getDiameterOutputArgsTest(double* diameter,  
                       int* numPair, int** pairIDs);
FC_ReturnCode _fc_getDiameterOutputArgsTestMeshes(double* diameter,  
		       int* numPair, int** pairIDs, FC_Mesh** pairMeshes); 
FC_ReturnCode _fc_getDiameterOutputArgsTestSubsets(double* diameter,  
		       int* numPair, int** pairIDs, FC_Subset** pairSubsets); 
FC_ReturnCode _fc_getSubsetDiameterCore(FC_Subset subset, FC_Variable* displ,
		       double* diameter, int* numPair, int** pairIDs);
FC_ReturnCode _fc_getSubsetsDiameterCore(int numSubset, FC_Subset* subsets, 
                       FC_Variable* displs, double* diameter, int* numPair, 
                       int** pairIDs, FC_Subset** pairSubsets);

// Normals helper
void _fc_calcSurfaceNormal(int numPoint, FC_Coords *points, FC_Vector *normal);

// Parameterization helpers
void _fc_calcQuadParams(FC_Coords *vertices, FC_Coords point, 
                       int* numParam, double** params); 
void _fc_calcQuadLocation(FC_Coords *vertices, int numParam, 
                       double* params, FC_Coords *point);

// Mesh entities' measurements helpers
void _fc_calcTriArea(int numDim, FC_Coords *vertices, double* area);
void _fc_calcQuadArea(int numDim, FC_Coords *vertices, double* area);
void _fc_calcTetVolume(FC_Coords *vertices, double* volume);
void _fc_calcPyramidVolume(FC_Coords *vertices, double* volume);
void _fc_calcPrismVolume(FC_Coords *vertices, double* volume);
void _fc_calcHexVolume(FC_Coords *vertices, double* volume);

// Proximity helpers
FC_ReturnCode _fc_getSubsetsProximityArgsTest(FC_Mesh* mesh1, 
	               FC_Subset* subset1, FC_Variable* displ1, 
                       FC_Mesh* mesh2, FC_Subset* subset2, 
		       FC_Variable* displ2, double min_dist, int* numPair,
                       int** mesh1VertIDs, int** mesh2VertIDs); 
FC_ReturnCode _fc_getSubsetsProximityCore(FC_Subset subset1, 
		       FC_Variable* displ1, FC_Subset subset2,
		       FC_Variable* displ2, double min_dist,
		       int* numPair, int** mesh1vertIDs, int** mesh2vertIDs);

// Range query helpers
FC_ReturnCode _fc_createVertexBinCore(int numVertex, int numDim, double* coords, 
		       FC_Coords lowers, FC_Coords uppers, 
                       FC_VertexBin** vertexBin);
FC_ReturnCode _fc_getVertexIDsWithinBoundingBoxUsingVertexBin(FC_Mesh mesh, FC_VertexBin* bin, 
		       FC_Coords lowpoint, FC_Coords highpoint, int* numVertex, 
		       int** vertexIDs);
FC_ReturnCode _fc_getVertexIDsWithinSphereUsingVertexBin(FC_Mesh mesh, 
                       FC_VertexBin* bin, FC_Coords point, double radius,
		       int* numVertex, int** vertexIDs);
FC_ReturnCode _fc_geomSmoothCore(int numMesh, FC_Mesh* meshes, int doDispl,
                       FC_Variable* coord_displs, FC_VertexBin** bins, 
                       FC_Variable* vars, double radius, int doLocalSmooth,
                       FC_Variable** smooth_vars);


// Containment function
FC_ReturnCode  _fc_doesElementContainPoint(
   FC_Mesh mesh, int element_id, double point[3], int dim,       
   int max_num_vert_per_face, int *num_vert_per_face,
   int *face_conns, double *vert_coords, int *hit_type);


#ifdef __cplusplus
}
#endif

#endif // _FC_GEOM_P_H_
