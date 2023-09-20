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
 * \file meshP.h
 * \brief Private declarations for the \ref Mesh module.
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/modules/meshP.h,v $
 * $Revision: 1.38 $ 
 * $Date: 2006/09/20 23:31:02 $
 */

#ifndef _FC_MESH_P_H_
#define _FC_MESH_P_H_

#ifdef __cplusplus
extern "C" {
#endif

// include these in basic data structures to simplify includes for the
// the rest of the library which won't ever use them
#include "tableP.h"
#include "fileio.h"
#include "fileioP.h"

/**
 * \ingroup  PrivateMesh
 *
 * \brief Mesh Slot Structure for the Mesh Table
 *
 * \description
 *
 *   This structure contains meta data and big data for a mesh.  First, it
 *   contains a header, a handle to the owning dataset, and a
 *   structure containing file IO info. In addition, it has information about
 *   the coordinates of it's vertices, it's topology (element to vertex
 *   relationships), and it knows what subsets and variables it owns (and what
 *   sequences it's associated with if any of the variables are sequence
 *   variables).
 *
 *   When a user asks for the coords as a variable, a reference to that
 *   variable is kept so that it can be re-fetched.  Also copy mesh, does not
 *   actually copy this variable but instead creates a coord variable on the
 *   new mesh.
 *
 *   A varSlot is associated with one step (e.g. time step) of a variable.
 *   However, access is through the more specific types of variables.  For
 *   instance, basic variables are those not associated with a sequence and
 *   have only varSlot each and can be accessed with 'numBasicVar' and
 *   'basicVarIDs'. The variable which provided the coordinates of the mesh is
 *   a basic variable but can also be accessed with member 'coordVarID' (there
 *   is only 1). (Not yet implemented: In addition, the non-coordinate basic
 *   vars can be accessed through 'numBasicNonCoordVar' and
 *   'basicNonCoordVarIDs'.)  Sequence variables are associated with a
 *   sequence. A sequence variable consists of an ordered group of
 *   varSlots. Some sequence information is repeated here from the seqSlots to
 *   ease access.
 *
 *   A mesh can only be made of one type of element.
 */
typedef struct {
  _FC_SlotHeader header;  /**< Header, which contains general purpose info. */
  FC_Dataset ds;   /**< Handle to the owning dataset. */
  //---File reference 
  _FC_MeshFileIOInfo fileInfo; /**< Holds info pertinent to file access. */
  //---About the mesh:  vertices and elements
  int dim;         /**< The dimensionality of the vertex coordinates. */
  int topodim;     /**< The topological dimensionality of the mesh (e.g. 2D 
		      for a quad mesh, even if the vertices have 3D coords). */
  int numVertex;   /**< The number of vertices in the mesh. */
  int numElement;  /**< The number of elements in the mesh. */
  FC_ElementType elemType; /**< The type of element mesh is made of. */
  //---About the mesh: edges and faces, this is built info
  int numEdge;     /**< The number of edges in the mesh. Default value of -1
		    means edge data has not been created. */
  int numFace;     /**< The number of faces in the mesh. Default value of -1
		    means face data has not been created. */
  //---Provided data about mesh coordinates and connectivity
  double *coords;  /**< Big Data: the coordinates of the vertices stored as an
                      array of coordinates values with order xyzxyzxyz... 
                      (or xyxyxy if only 2D). This array will be numVertex*dim
                      long. */
  int *elemToVertConns;  /**< Big Data: the connectivity of the mesh. This is 
		      stored as an array of vertex indices for each element 
		      (e.g., for a mesh of two triangles and four vertices, the
		      array could be { 0 1 2 1 3 2 } where the first three 
		      numbers are the IDs of the vertices which are  make up 
		      the first triangle and then next three numbers are the 
		      vertices of the second triangle). Note that PATRAN 
		      ordering of vertices of an element is assumed. This array
		      will be numElement*(number of vertices per element) 
		      long. */
  //---Built data -- for edges & faces (downward conns)
  int* edgeToVertConns;  /**< Big Data: the vertices that make up the edges 
			    listed as E1V1 E1V2 ... EnV1 EnV2 */
  FC_ElementType* faceTypes; /**< Big Data: For each face, the face type */
  int* numVertPerFace;   /**< Big Data: For each face, the number of verts */
  int* faceToVertConns;  /**< Big Data: the vertices that make up the faces.
			    listed as all verts for first face, then all verts
                            for next face, etc. NOTE: For a list of mixed type
			    faces, for each face there is space for the maximum
                            possible number of verts; and if a face
                            has fewer verts, the extra spots hold -1 */
  int* elemToEdgeConns;  /**< Big Data: element to edge connectivity. For 
			    each element, the IDs of it's edges */
  int* elemToFaceConns;  /**< Big Data: element to face connectivity. For 
			    each element, the IDs of it's faces */
  int* elemFaceOrients;  /**< Big Data: the relative orientation of the face
			    in the element compared to the global face's */
  //---Built data -- parents (upward "conns") for all
  int* numEdgePerVert;   /**< Big Data: for each vert, # of parent edges */
  int** edgeParentsPerVert; /**< Big Data: for each vert, array of parent edge
			     IDs */
  int* numFacePerVert;   /**< Big Data: for each vert, # of parent faces. */
  int** faceParentsPerVert; /**< Big Data: for each vert, array of parent face
			     IDs */
  int* numElemPerVert;   /**< Big Data: for each vert, # of parent elems. */
  int** elemParentsPerVert; /**< Big Data: for each vert, array of parent
			     element IDs */
  int* numElemPerEdge;   /**< Big Data: for each edge, # of parent elems. */
  int** elemParentsPerEdge; /**< Big Data: for each edge, array of parent
			    element IDs */
  int* numElemPerFace;   /**< Big Data: for each face, # of parent elems. */
  int** elemParentsPerFace; /**< Big Data: for each face, array of parent
			     element IDs */
  //---Built data -- neighbors
  int* numVertNeighsViaEdge; /**< Big Data: for each vert, # of neighbor
				verts connected by an edge */
  int** vertNeighsViaEdge;   /**< Big Data: for each vert, array of neighbor
				verts connected by an edge */ 
  int* numElemNeighsViaVert; /**< Big Data: for each elem, # of neighbor
				elems connected by an vert */
  int** elemNeighsViaVert;   /**< Big Data: for each elem, array of neighbor
				elems connected by an vert */ 
  int* numElemNeighsViaEdge; /**< Big Data: for each elem, # of neighbor
				elems connected by an edge */
  int** elemNeighsViaEdge;   /**< Big Data: for each elem, array of neighbor
				elems connected by an edge */ 
  int* numElemNeighsViaFace; /**< Big Data: for each elem, # of neighbor
				elems connected by an face */
  int** elemNeighsViaFace;   /**< Big Data: for each elem, array of neighbor
				elems connected by an face */ 
  //---Owned entities
  //---subsets
  int numSub;         /**< Number of subsets owned by this mesh */
  int *subIDs;        /**< Array of slot numbers identifying subsets */
  //---Variable breakdown: Access of variables by type = coord, basic, seq
  int coordVarID; /**< Slot ID of the coordinate array. This is -1 if
                       if there are no coords on the mesh, or the slotID 
                       of one of the basic variables. */
  int numBasicVar; /**< Number of basic variables (i.e. non-sequence vars) */
  int* basicVarIDs; /** < Array of slot IDs for the basic variables */
  int numSeqVar;    /**< Number of sequence variables owned by this mesh. */
  int *numStepPerSeqVar; /**< Number of steps in each seq var */
  int **seqVarIDs;   /**< Array of slot IDs identifying the vars that make 
                          up the sequence vars. This array has the dimensions
                          [numSeqVar][numStepPerSeq]. */
} _FC_MeshSlot;

// Private big data access - parents
FC_ReturnCode _fc_getMeshEdgeParentsOfVerticesPtr(FC_Mesh mesh, 
                      int** numEdgePerVert, int*** edgeParentsPerVert);
FC_ReturnCode _fc_getMeshFaceParentsOfVerticesPtr(FC_Mesh mesh, 
                      int** numFacePerVert, int*** edgeParentsPerVert);
FC_ReturnCode _fc_getMeshElementParentsOfVerticesPtr(FC_Mesh mesh, 
                      int** numElemPerVert, int*** elemParentsPerVert);
FC_ReturnCode _fc_getMeshElementParentsOfEdgesPtr(FC_Mesh mesh, 
                      int** numElemPerEdge, int*** elemParentsPerEdge);
FC_ReturnCode _fc_getMeshElementParentsOfFacesPtr(FC_Mesh mesh, 
                      int** numElemPerFace, int*** elemParentsPerFace);
// Private big data access - neighbors
FC_ReturnCode _fc_getMeshVertexNeighborsViaEdgePtr(FC_Mesh mesh, 
		      int** numVertNeighsViaEdge, int*** vertNeighsViaEdge);
FC_ReturnCode _fc_getMeshElementNeighborsViaVertexPtr(FC_Mesh mesh, 
		      int** numElemNeighsViaVert, int*** elemNeighsViaVert);
FC_ReturnCode _fc_getMeshElementNeighborsViaEdgePtr(FC_Mesh mesh, 
		      int** numElemNeighsViaEdge, int*** elemNeighsViaEdge);
FC_ReturnCode _fc_getMeshElementNeighborsViaFacePtr(FC_Mesh mesh, 
		      int** numElemNeighsViaFace, int*** elemNeighsViaFace);

// helper functions to build conns, parents & neighbors
FC_ReturnCode _fc_buildEdgeConns(FC_Mesh mesh);
FC_ReturnCode _fc_buildFaceConns(FC_Mesh mesh);
FC_ReturnCode _fc_buildParents(FC_Mesh mesh, FC_AssociationType childType,
		      FC_AssociationType parentType);
FC_ReturnCode _fc_buildMeshVertexNeighborsViaEdge(FC_Mesh mesh);
FC_ReturnCode _fc_buildMeshElementNeighborsViaEntity(FC_Mesh mesh, 
                      FC_AssociationType sharedType);

// slot stuff
void _fc_initMeshSlot(_FC_MeshSlot* meshSlot);
void _fc_releaseMeshSlot(_FC_MeshSlot *meshSlot);
void _fc_clearMeshSlot(_FC_MeshSlot *meshSlot);

// table access
int _fc_getMeshTableSize(void);
_FC_MeshSlot* _fc_getNewMeshSlot(void);
_FC_MeshSlot* _fc_getMeshSlot(FC_Mesh mesh);
_FC_MeshSlot* _fc_getMeshSlotFromID(int meshID);
FC_ReturnCode _fc_deleteMeshSlot(FC_Mesh mesh);
void _fc_printMeshTable(char *label);
void _fc_freeMeshTable(void);


#ifdef __cplusplus
}
#endif

#endif // _FC_MESH_P_H_

