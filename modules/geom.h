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
 * \file geom.h
 * \brief Public declarations for \ref GeometricRelations module.
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/modules/geom.h,v $
 * $Revision: 1.86 $ 
 * $Date: 2007/07/06 05:54:35 $
 */

#ifndef _FC_GEOM_H_
#define _FC_GEOM_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "statistics.h"

// General geometry
FC_ReturnCode fc_calcEuclideanDistance(FC_Coords v1Coords, FC_Coords v2Coords,
                                       int dim, double* distance );
FC_ReturnCode fc_calcSquaredEuclideanDistance(FC_Coords v1Coords, 
                                       FC_Coords v2Coords,
                                       int dim, double* distance2 );
FC_ReturnCode fc_calcAngleBetweenVectors(FC_Vector vect1, FC_Vector vect2, 
					 double* angle);

// Displaced coordinates
int fc_isValidDisplacementVariable(FC_Mesh mesh, FC_Variable variable);
FC_ReturnCode fc_getDisplacedMeshCoords(FC_Mesh mesh, FC_Variable displs,
                       double** displaced_coords);

// Bounding box
int fc_isBoundingBoxValid(int dim, FC_Coords lowpoint, FC_Coords highpoint);
FC_ReturnCode fc_getBoundingBoxesOverlap(int dim, FC_Coords lowpoint1,
                       FC_Coords highpoint1, FC_Coords lowpoint2,
                       FC_Coords highpoint2, int* intersect_flag,
                       FC_Coords* intersect_lowpoint, 
                       FC_Coords* intersect_highpoint);
FC_ReturnCode fc_combineBoundingBoxes(int dim, FC_Coords lowpoint1,
                       FC_Coords highpoint1, FC_Coords lowpoint2, 
                       FC_Coords highpoint2, FC_Coords* combine_lowpoint,
                       FC_Coords* combine_highpoint);
FC_ReturnCode fc_getMeshBoundingBox(FC_Mesh mesh, int* dim, 
                       FC_Coords *lowpoint, FC_Coords *highpoint);
FC_ReturnCode fc_getDisplacedMeshBoundingBox(FC_Mesh mesh, 
                       FC_Variable coord_displ, int* dim, FC_Coords *lowpoint, 
                       FC_Coords* highpoint);
FC_ReturnCode fc_getSubsetBoundingBox(FC_Subset, int *dim, FC_Coords* lowpoint,
                       FC_Coords* highpoint);   
FC_ReturnCode fc_getDisplacedSubsetBoundingBox(FC_Subset, 
                       FC_Variable coord_displ, int *dim, FC_Coords* lowpoint,
                       FC_Coords* highpoint);
FC_ReturnCode fc_createBoundingBoxMesh(FC_Dataset dest_ds, char* mesh_name,
                       int dim, FC_Coords lowpoint, FC_Coords highpoint,
                       FC_Mesh* bb_mesh);
// ??
FC_ReturnCode fc_projectedElementBoundingBoxIntersection(
  FC_Mesh mesh, int elemID_1, int elemID_2,  FC_Coords proj_normal,
  int* intersectval); 

FC_ReturnCode fc_projectedSubsetElementBoundingBoxIntersection(
  FC_Subset compelemIDs, FC_Subset refelemIDs,
  FC_Coords proj_normal, FC_Subset* intersectingrefIDs);


// Diameter (max distance between any points in a point set)
FC_ReturnCode fc_getMeshDiameter(FC_Mesh mesh, double* diameter, int* numPair,
		       int** pairIDs);
FC_ReturnCode fc_getDisplacedMeshDiameter(FC_Mesh mesh, 
                       FC_Variable coord_displ, double* diameter, 
                       int* numPair, int** pairIDs);
FC_ReturnCode fc_getSubsetDiameter(FC_Subset subset, double* diameter, 
                       int* numPair, int** pairIDs);
FC_ReturnCode fc_getDisplacedSubsetDiameter(FC_Subset subset, 
                       FC_Variable coord_displ, double* diameter,
                       int* numPair, int** pairIDs);
FC_ReturnCode fc_getMeshesDiameter(int numMesh, FC_Mesh* mesh, 
		       double* diameter, int* numPair, int** pairIDs, 
                       FC_Mesh** pairMeshes);
FC_ReturnCode fc_getDisplacedMeshesDiameter(int numMesh, FC_Mesh* meshes, 
                       FC_Variable* coord_displ, double* diameter, 
                       int* numPair, int** pairIDs, FC_Mesh** pairMeshes);
FC_ReturnCode fc_getSubsetsDiameter(int numSubset, FC_Subset* subsets, 
		       double* diameter, int* numPair, int** pairIDs,
		       FC_Subset** pairSubsets);
FC_ReturnCode fc_getDisplacedSubsetsDiameter(int numSubset, FC_Subset* subsets,
                       FC_Variable* coord_displs, double* diameter,
                       int* numPair, int** pairIDs, FC_Subset** pairSubsets);

// Centroid
FC_ReturnCode fc_getMeshCentroid(FC_Mesh mesh, int* dim, FC_Coords *centroid);
FC_ReturnCode fc_getDisplacedMeshCentroid(FC_Mesh mesh, 
		       FC_Variable coord_displ, int* dim, FC_Coords *centroid);
FC_ReturnCode fc_getVariableWeightedMeshCentroid(FC_Mesh mesh, 
		       FC_Variable weights, int* dim, FC_Coords *centroid);
FC_ReturnCode fc_getVariableWeightedDisplacedMeshCentroid(FC_Mesh mesh,
		       FC_Variable coord_displ, FC_Variable weights, int* dim, 
                       FC_Coords *centroid);
FC_ReturnCode fc_getSubsetCentroid(FC_Subset subset, int* dim, 
                       FC_Coords *centroid);
FC_ReturnCode fc_getDisplacedSubsetCentroid(FC_Subset subset, 
		       FC_Variable coord_displ, int* dim, FC_Coords *centroid);
FC_ReturnCode fc_getVariableWeightedSubsetCentroid(FC_Subset subset, 
		       FC_Variable weights, int* dim, FC_Coords *centroid);
FC_ReturnCode fc_getVariableWeightedDisplacedSubsetCentroid(FC_Subset subset,
		       FC_Variable coord_displ, FC_Variable weights, int* dim, 
                       FC_Coords *centroid);

// Mesh entities' measurements
FC_ReturnCode fc_getEdgeLengths(FC_Mesh mesh, FC_Variable* edge_lengths);
FC_ReturnCode fc_getDisplacedEdgeLengths(FC_Mesh mesh, FC_Variable coord_displ,
                       FC_Variable* edge_lengths);
FC_ReturnCode fc_getFaceAreas(FC_Mesh mesh, FC_Variable* face_areas);
FC_ReturnCode fc_getDisplacedFaceAreas(FC_Mesh mesh, FC_Variable coord_displ, 
                       FC_Variable* face_areas);
FC_ReturnCode fc_getRegionVolumes(FC_Mesh mesh, FC_Variable* region_volumes);
FC_ReturnCode fc_getDisplacedRegionVolumes(FC_Mesh mesh,
                       FC_Variable coord_displ, FC_Variable* region_volumes);
FC_ReturnCode fc_getMeshVolume(FC_Mesh mesh,double*volume);
FC_ReturnCode fc_getDisplacedMeshVolume(FC_Mesh mesh, FC_Variable displs,
                       double* volume);
FC_ReturnCode fc_getSubsetVolume(FC_Subset, int doStrict, double* volume);
FC_ReturnCode fc_getDisplacedSubsetVolume(FC_Subset, int doStrict,
					  FC_Variable displs, double* volume);
FC_ReturnCode fc_getSubsetArea(FC_Subset, double* area);

// Return normals on all entries in a mesh as a data variable
FC_ReturnCode fc_getSurfaceNormals(FC_Mesh mesh, FC_Variable *new_var);
FC_ReturnCode fc_getDisplacedSurfaceNormals(FC_Mesh mesh, 
                       FC_Variable coord_displ, FC_Variable *new_var);

// Separate mesh deformation from mesh displacement
FC_ReturnCode fc_calcMeshDeformationUsingRefVertex(FC_Mesh mesh, int vertexID,
                       int numStep, FC_Variable* displ, FC_Variable** deform);
FC_ReturnCode fc_calcMeshDeformationUsingCentroid(FC_Mesh mesh, int numStep, 
		       FC_Variable* displ, FC_Variable** deform);
// future options: 1) using 3 vertices, 2) global fit

// Proximity
FC_ReturnCode fc_getMeshesProximity(FC_Mesh mesh1, FC_Mesh mesh2, 
		       double min_dist, int* numPair, int** mesh1VertIDs, 
		       int** mesh2VertIDs);
FC_ReturnCode fc_getDisplacedMeshesProximity(FC_Mesh mesh1, 
		       FC_Variable displ1, FC_Mesh mesh2, FC_Variable displ2,
		       double min_dist, int* numPair, int** mesh1VertIDs, 
		       int** mesh2VertIDs);
FC_ReturnCode fc_getSubsetsProximity(FC_Subset subset1, FC_Subset subset2,
		       double min_dist, int* numPair, int** mesh1VertIDs,
		       int** mesh2VertIDs);
FC_ReturnCode fc_getDisplacedSubsetsProximity(FC_Subset subset1, 
		       FC_Variable displ1, FC_Subset subset2, 
                       FC_Variable displ2, double min_dist, int* numPair, 
                       int** mesh1VertIDs, int** mesh2VertIDs);

/**
 * \ingroup GeometricRelations
 * \brief Spatially binned vertex IDs
 *
 * \description
 *
 *    This is a helper structure to look up vertices by spatial
 *    location. Because the bin structure is a regularly structured
 *    grid, it is easy to identify which bin a given set of 
 *    coordinates is in. Then search for a vertexID is reduced to
 *    searching only the vertexIDs in this single bin instead of
 *    the entire mesh.
 *
 * \todo ? someday move mesh bounding box into mesh
 * \modifications
 *    - 03/08/04 WSK, Created.
 */
typedef struct {
  FC_Mesh mesh;            /**< The owning mesh */
  FC_Variable displ;       /**< The displacement variable or FC_NULL_VARIABLE */
  int numDim;              /**< Dimensionality of coordinate system of mesh */
  // Had to change these to get them to swig
  //FC_CoordslowerCoords;   /**< Bounding box of owning mesh */
  //FC_Coords upperCoords;   /**< Bounding box of owning mesh */
  double lowerCoords[3];   /**< Bounding box of owning mesh */
  double upperCoords[3];   /**< Bounding box of owning mesh */
  double widths[3];        /**< Widths of bounding box of owning mesh */
  int numBins[3];          /**< Number of bins in each dimensions */
  double binWidths[3];     /**< Width of each dimension of the bins */
  int*** numVertex;        /**< Number of vertices per bin, eg num[i][j][k] */
  int**** vertexIDs;       /**< Vertex IDs per bin, eg IDs[i][j][k][n] */
} FC_VertexBin;

FC_ReturnCode fc_createMeshVertexBin(FC_Mesh mesh, FC_VertexBin** vertexBin);
FC_ReturnCode fc_createDisplacedMeshVertexBin(FC_Mesh mesh, 
		       FC_Variable coord_displ, FC_VertexBin** vertexBin);
void fc_freeVertexBin(FC_VertexBin* bin);

FC_ReturnCode fc_getEntitiesWithinBoundingBox(FC_Mesh mesh, FC_VertexBin* bin,
     	               FC_Coords lowpoint, FC_Coords highpoint,
  	               FC_AssociationType entity_type, int* numEntity, 
                       int** entityIDs);
FC_ReturnCode fc_getEntitiesWithinSphere(FC_Mesh mesh, FC_VertexBin* bin, 
                       FC_Coords center, double radius, 
                       FC_AssociationType entity_type, int *numEntity, 
                       int** entityIDs);

FC_ReturnCode fc_geomSmoothVariable(int numMesh, FC_Mesh* meshes, 
                       FC_VertexBin** bins, FC_Variable* vars, double radius, 
                       int doLocalSmooth, FC_Variable** smooth_vars);
FC_ReturnCode fc_displacedGeomSmoothVariable(int numMesh, FC_Mesh* meshes,
                       FC_Variable* coord_displs, FC_VertexBin** bins, 
                       FC_Variable* vars, double radius, int doLocalSmooth,
                       FC_Variable** smooth_vars);
FC_ReturnCode fc_geomSmoothSeqVariable(int numMesh, FC_Mesh* meshes, 
                       FC_VertexBin** bins, int numStep, 
                       FC_Variable** seqVars, double radius, 
                       int doLocalSmooth, FC_Variable*** smooth_seq_vars);
FC_ReturnCode fc_displacedGeomSmoothSeqVariable(int numMesh, FC_Mesh* meshes, 
		       int numStepDispl, FC_Variable** seqDispls,
                       FC_VertexBin*** bins, int numStep, 
                       FC_Variable** seqVars, double radius, 
                       int doLocalSmooth, FC_Variable*** smooth_seq_vars);


#ifdef __cplusplus
}
#endif

#endif  // _FC_GEOM_H_
