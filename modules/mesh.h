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
 * \file mesh.h
 * \brief Public declarations for the \ref Mesh Module.
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/modules/mesh.h,v $
 * $Revision: 1.20 $ 
 * $Date: 2006/09/22 04:40:45 $
 */

#ifndef _FC_MESH_H_
#define _FC_MESH_H_

#ifdef __cplusplus
extern "C" {
#endif

// create new mesh from scratch
//-------------------------------
FC_ReturnCode fc_createMesh(FC_Dataset dataset, char *meshname, 
		      FC_Mesh *mesh); 
FC_ReturnCode fc_setMeshCoords(FC_Mesh mesh, int dim, int numVert, 
                      double *coords);
FC_ReturnCode fc_setMeshElementConns(FC_Mesh mesh, FC_ElementType elemType, 
                      int numElem, int *conns);
FC_ReturnCode fc_setMeshCoordsPtr(FC_Mesh mesh, int dim, int numVert, 
                      double *coords_p);
FC_ReturnCode fc_setMeshElementConnsPtr(FC_Mesh mesh, FC_ElementType elemType, 
                      int numElem, int *conn_p);

// other ways to get new meshes
//-------------------------------
FC_ReturnCode fc_copyMesh(FC_Mesh src_mesh, FC_Dataset dest_ds, 
                      char* newMeshName, FC_Mesh* new_mesh);
FC_ReturnCode fc_createSubsetMesh(FC_Subset src_subset, FC_Dataset dest_ds, 
		      int doStrict, char* newName, FC_Mesh* subset_mesh);
FC_ReturnCode fc_createSimpleHexMesh(FC_Dataset dataset, char* newMeshName, 
                      int N, int M, int P, double* lower_coords, 
                      double* upper_coords, FC_Mesh* new_mesh);

// get meshes
//-------------------------------
FC_ReturnCode fc_getNumMesh(FC_Dataset dataset, int *numMesh);
FC_ReturnCode fc_getMeshes(FC_Dataset dataset, int *numMesh, FC_Mesh **meshes);
FC_ReturnCode fc_getMeshByName(FC_Dataset dataset, char *meshName,
			       int *numMeshes, FC_Mesh **mesh);  

// change the name of a mesh
//-------------------------------
FC_ReturnCode fc_changeMeshName(FC_Mesh mesh, char* newName);

// release/delete mesh
//-------------------------------
FC_ReturnCode fc_releaseMesh(FC_Mesh mesh);
FC_ReturnCode fc_deleteMesh(FC_Mesh mesh);

// get Mesh metadata
//-------------------------------
int fc_isMeshValid(FC_Mesh mesh);
FC_ReturnCode fc_getMeshName(FC_Mesh mesh, char** meshname);
FC_ReturnCode fc_getDatasetFromMesh(FC_Mesh mesh, FC_Dataset *dataset);
FC_ReturnCode fc_getMeshInfo(FC_Mesh mesh, int *topodim, int *dim, 
		      int *numVertex, int *numElement, 
		      FC_ElementType *elemType); 
FC_ReturnCode fc_getMeshTopodim(FC_Mesh mesh, int *topodim);
FC_ReturnCode fc_getMeshDim(FC_Mesh mesh, int *dim);             
FC_ReturnCode fc_getMeshNumVertex(FC_Mesh mesh, int *numVertex);           
FC_ReturnCode fc_getMeshNumElement(FC_Mesh mesh, int *numElement);
FC_ReturnCode fc_getMeshElementType(FC_Mesh mesh, FC_ElementType *elemType);
FC_ReturnCode fc_getMeshEdgeFaceInfo(FC_Mesh mesh, int* numEdge,
  	              int* numFace, FC_ElementType *facetype);
FC_ReturnCode fc_getMeshNumEdge(FC_Mesh mesh, int* numEdge);
FC_ReturnCode fc_getMeshNumFace(FC_Mesh mesh, int* numFace);
FC_ReturnCode fc_getMeshFaceType(FC_Mesh mesh, FC_ElementType *faceType);
FC_ReturnCode fc_getMeshNumEntity(FC_Mesh mesh, FC_AssociationType assoc, 
		      int *numEntity);
FC_ReturnCode fc_getMeshEntityElementType(FC_Mesh mesh, 
                      FC_AssociationType assoc, int ID, 
                      FC_ElementType* elemType);

// get big data
//-------------------------------
FC_ReturnCode fc_getMeshCoordsPtr(FC_Mesh mesh, double **coords_p);
FC_ReturnCode fc_getMeshCoordsAsVariable(FC_Mesh mesh, FC_Variable *coords);
FC_ReturnCode fc_getMeshEdgeConnsPtr(FC_Mesh mesh, int** conns_p);
FC_ReturnCode fc_getMeshFaceTypesPtr(FC_Mesh mesh, FC_ElementType** faceTypes);
FC_ReturnCode fc_getMeshFaceConnsPtr(FC_Mesh mesh, int** numVertexPerFace_p, 
		      int* stride, int** conns_p);
FC_ReturnCode fc_getMeshElementConnsPtr(FC_Mesh mesh, int **conns_p);
FC_ReturnCode fc_getMeshElementToEdgeConnsPtr(FC_Mesh mesh, int **elemToEdgeConns);
FC_ReturnCode fc_getMeshElementToFaceConnsPtr(FC_Mesh mesh, int **elemToFaceConns);
FC_ReturnCode fc_getMeshElementFaceOrientsPtr(FC_Mesh mesh, int **elemFaceOrients);

// Print
//-------------------------------
FC_ReturnCode fc_printMesh(FC_Mesh mesh, char* label, int print_edgeConns, 
		 int print_faceConns, int print_elemConns, int print_coords);

#ifdef __cplusplus
}
#endif

#endif // _FC_MESH_H_
