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
 * \file topoP.h
 * \brief Private declarations for \ref TopologyRelations module.
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/modules/topoP.h,v $
 * $Revision: 1.40 $ 
 * $Date: 2006/08/30 19:20:01 $
 */

#ifndef _FC_TOPO_P_H_
#define _FC_TOPO_P_H_

#ifdef __cplusplus
extern "C" {
#endif

// membership relations by entity
FC_ReturnCode _fc_getVertexParents(FC_Mesh mesh, int vertexID, 
			FC_AssociationType parent_assoc, int* numParent,
			int** parentIDs);
FC_ReturnCode _fc_getEdgeParents(FC_Mesh mesh, int edgeID, 
			FC_AssociationType parent_assoc, int* numParent,
			int** parentIDs);
FC_ReturnCode _fc_getFaceParents(FC_Mesh mesh, int faceID, 
			FC_AssociationType parent_assoc, int* numParent,
			int** parentIDs);

// Neighbor relations by entity
FC_ReturnCode _fc_getVertexNeighbors(FC_Mesh mesh, int vertexID, int* numNeigh,
			int** neighborIDs);
FC_ReturnCode _fc_getEdgeNeighbors(FC_Mesh mesh, int edgeID, int* numNeigh,
			int** neighborIDs);
FC_ReturnCode _fc_getFaceNeighbors(FC_Mesh mesh, int faceID, 
			int shared_dim, int* numNeigh, int** neighborIDs);
FC_ReturnCode _fc_getElementNeighbors(FC_Mesh mesh, int elementID, 
			int shared_dim, int* numNeigh, int** neighborIDs);
FC_ReturnCode _fc_growNeighborList(FC_Mesh mesh, FC_SortedIntArray* newIDs, 
                        FC_AssociationType assoc, int level, int shared_dim, 
                        int* seen, FC_SortedIntArray* keepIDs);

#ifdef __cplusplus
}
#endif

#endif // _FC_TOPO_P_H_
