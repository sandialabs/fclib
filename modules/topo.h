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
 * \file topo.h
 * \brief Public declarations for \ref TopologyRelations module.
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/modules/topo.h,v $
 * $Revision: 1.46 $ 
 * $Date: 2006/08/30 19:20:01 $
 */

#ifndef _FC_TOPO_H_
#define _FC_TOPO_H_

#ifdef __cplusplus
extern "C" {
#endif

//??

// membership relations
FC_ReturnCode fc_getMeshEntityChildren(FC_Mesh mesh, 
                            FC_AssociationType parent_assoc,
			    int parentID, FC_AssociationType child_assoc,
			    int* numChild, int** childIDs);
FC_ReturnCode fc_getMeshEntityParents(FC_Mesh mesh, 
                            FC_AssociationType child_assoc, 
			    int childID, FC_AssociationType parent_assoc, 
			    int* numParent, int** parentIDs);
FC_ReturnCode fc_changeMeshEntityType(FC_Mesh mesh,
			    FC_AssociationType old_assoc, int numOldId,
			    int *oldIDs, FC_AssociationType new_assoc,
			    int doStrict, int* numNewID, int** newIDs);

// neighbors
FC_ReturnCode fc_getMeshEntityNeighbors(FC_Mesh mesh, int entityID, 
                            FC_AssociationType assoc, int shared_dim,
		            int* numNeighbor, int** neighborIDs);
FC_ReturnCode fc_getSubsetNeighbors(FC_Subset subset, int level, 
			    int shared_dim, FC_AssociationType *neigh_assoc, 
                            int* numNeighbor, int** neighborIDs);

// create connected components
FC_ReturnCode fc_segment(FC_Subset subset, int shared_dim, int *numSegment, 
			 FC_Subset** newSubsets);

// skin
FC_ReturnCode fc_getMeshSkin(FC_Mesh mesh, FC_Subset *skin);
FC_ReturnCode fc_getSubsetSkin(FC_Subset subset, FC_Subset *skin); 
  //FC_ReturnCode fc_getGenusOfSubsetSkin(FC_Subset subset, int* genus);

#ifdef __cplusplus
}
#endif

#endif // _FC_TOPO_H_
