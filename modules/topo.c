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
 * \file topo.c
 * \brief Implementation for \ref TopologyRelations module.
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/modules/topo.c,v $
 * $Revision: 1.151 $ 
 * $Date: 2006/10/19 03:14:51 $
 *
 * \modifications
 *    - 12/17/03 WSK Moved minimum spanning tree stuff from here
 *      to it's own file.
 */

// C library includes
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

// fc library dependencies
#include "base.h"
#include "storage.h"
#include "storageP.h"
#include "library.h"
#include "mesh.h"
#include "meshP.h"
#include "subset.h"
#include "subsetP.h"
#include "geom.h"
#include "geomP.h"

// this module
#include "topo.h"
#include "topoP.h"

/**
 * \addtogroup  TopologyRelations
 * \brief Relationships between elements and vertices of the mesh.
 * 
 * \description 
 *  
 *    The topology functions compute relationships between elements and
 *    vertices of the mesh itself, with no consideration of the actual
 *    physical coordinates of the vertices.  So topology functions include
 *    ones such as compute which elements share a vertex, 
 *    which vertices are part of an element, which vertices are 
 *    part of an element which shares a face with a given element, etc.
 */

/**
 * \ingroup TopologyRelations
 * \defgroup PrivateTopologyRelations (Private)
 */

/** \name Membership relations. */
//-------------------------------
//@{

/**
 * \ingroup  TopologyRelations
 * \brief Get the entities which make up the specified entity.
 *
 * \description
 *   
 *    Given an entity on a mesh (e.g. vertex, element, etc) find the specified
 *    type of children of that entity. The specified child has to have a
 *    topological dimensionality less than the parent entity. (E.g., vertices
 *    and edges can be children of a face, but elements cannot.)
 *
 *    Entities cannot be of association FC_AT_WHOLE_DATASET
 *
 * \todo Make sure this only calls the edge/face stuff when it's needed
 *       Clean up!
 *
 * \modifications
 *    - 07/16/04 WSK Created.
 *    - 11/01/2005 WSD. Fixed up so lazy building of edges & faces 
 *      doesn't get called where not needed.
 */
FC_ReturnCode fc_getMeshEntityChildren(
  FC_Mesh mesh,     /**< input - mesh */
  FC_AssociationType parent_assoc,  /**< input - the parent's entity type */
  int parentID,      /**< input - the parent's ID (e.g. element ID, etc) */
  FC_AssociationType child_assoc,  /**< input - the type of mesh entities to
                                      return (must have topo dimensionality < 
                                      parent_assoc. */
  int *numChild,    /**< output - the number of children */
  int **childIDs    /**< output - array of child IDs */
) {
  FC_ReturnCode rc;
  int i, j, k;
  int numParent;
  FC_ElementType elemType;
  int maxNumChildPerParent, *parentToChildConns;
  
  // default return values
  if (numChild)
    *numChild = -1;
  if (childIDs)
    *childIDs = NULL;

  // check input
  if (!fc_isMeshValid(mesh) || !fc_isAssociationTypeValid(parent_assoc) || 
      parent_assoc == FC_AT_UNKNOWN || parent_assoc == FC_AT_WHOLE_DATASET ||
      parentID < 0 || !fc_isAssociationTypeValid(child_assoc) || 
      child_assoc == FC_AT_UNKNOWN || parent_assoc <= child_assoc || 
      !numChild || !childIDs) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  rc = fc_getMeshNumEntity(mesh, parent_assoc, &numParent);
  if (rc != FC_SUCCESS || parentID >= numParent) {
        fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR)); 
    return FC_INPUT_ERROR;
  }

  // log message
  fc_printfLogMessage("Getting entities which make up the specified entity.");
  
  // special case 1: if whole, return list that is all of the child entities
  if (parent_assoc == FC_AT_WHOLE_MESH) {
    rc = fc_getMeshNumEntity(mesh, child_assoc, numChild);
    if (*numChild > 0) {
      *childIDs = (int*)malloc(*numChild*sizeof(int));
      if (*childIDs == NULL) {
        fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
        return FC_MEMORY_ERROR;
      }
      for (i = 0; i < *numChild; i++)
        (*childIDs)[i] = i;
    }
    return FC_SUCCESS;
  }

  // special case 2: face->edge conns don't exist on the mesh, so do
  // a local calc
  if (child_assoc == FC_AT_EDGE && parent_assoc == FC_AT_FACE) {
    int elemID, edgeID, found;
    int numEdgePerElem, *edgeConns, *elemToEdgeConns;
    int *numVertPerFace, maxNumVertPerFace, *faceConns;
    int *numElemPerFace, **elemParentsPerFace;
    FC_SortedIntArray sia;

    // setup
    rc = fc_getMeshElementType(mesh, &elemType);
    if (rc != FC_SUCCESS)
      return rc; 
    rc = fc_getMeshEdgeConnsPtr(mesh, &edgeConns);
    if (rc != FC_SUCCESS)
      return rc; 
    rc = fc_getMeshFaceConnsPtr(mesh, &numVertPerFace, &maxNumVertPerFace,
                                &faceConns);
    if (rc != FC_SUCCESS)
      return rc; 
    rc = fc_getMeshElementToEdgeConnsPtr(mesh, &elemToEdgeConns);
    if (rc != FC_SUCCESS)
      return rc; 
    rc = _fc_getMeshElementParentsOfFacesPtr(mesh, &numElemPerFace,
                                          &elemParentsPerFace);
    if (rc != FC_SUCCESS)
      return rc;
    numEdgePerElem = fc_getElementTypeNumEdge(elemType);

    // Get an element that owns the face
    elemID = elemParentsPerFace[parentID][0];
    // For each edge in element, if edge's verts are all in the face, add!
    fc_initSortedIntArray(&sia);
    for (i = 0; i < numEdgePerElem; i++) {
      edgeID = elemToEdgeConns[elemID*numEdgePerElem+i];
      found = 0;
      for (j = 0; j < 2; j++) 
        for (k = 0; k < numVertPerFace[parentID]; k++)
          if (edgeConns[edgeID*2+j] == 
              faceConns[parentID*maxNumVertPerFace+k]) 
            found++;
      if (found == 2) { // both verts matched!
        int ret = fc_addIntToSortedIntArray(&sia, edgeID);
        if (ret < 0) {
          fc_printfErrorMessage("Failed to add edge to sia");
          return ret;
        }
      }
    }

    // all done
    rc = fc_convertSortedIntArrayToIntArray(&sia, numChild, childIDs);
    if (rc != FC_SUCCESS)
      return rc;

    return FC_SUCCESS;
  }

  // Get numChild & ids
  rc = fc_getMeshEntityElementType(mesh, parent_assoc, parentID, &elemType);
  if (rc != FC_SUCCESS)
    return rc;
  *numChild = fc_getElementTypeNumEntity(elemType, child_assoc);
  if (*numChild > 0) {
    *childIDs = (int*)malloc(*numChild*sizeof(int));
    if (*childIDs == NULL) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      return FC_MEMORY_ERROR;
    }
    maxNumChildPerParent = *numChild;
    if (child_assoc == FC_AT_VERTEX && parent_assoc == FC_AT_EDGE)
      rc = fc_getMeshEdgeConnsPtr(mesh, &parentToChildConns);
    else if (child_assoc == FC_AT_VERTEX && parent_assoc == FC_AT_FACE)
      rc = fc_getMeshFaceConnsPtr(mesh, NULL, &maxNumChildPerParent,
                                  &parentToChildConns);
    else if (child_assoc == FC_AT_VERTEX && parent_assoc == FC_AT_ELEMENT) 
      rc = fc_getMeshElementConnsPtr(mesh, &parentToChildConns);
    else if (child_assoc == FC_AT_EDGE && parent_assoc == FC_AT_ELEMENT) 
      rc = fc_getMeshElementToEdgeConnsPtr(mesh, &parentToChildConns);
    else if (child_assoc == FC_AT_FACE && parent_assoc == FC_AT_ELEMENT) 
      rc = fc_getMeshElementToFaceConnsPtr(mesh, &parentToChildConns);
    memcpy(*childIDs, &parentToChildConns[parentID*maxNumChildPerParent], 
           *numChild*sizeof(int));
  }

  return FC_SUCCESS;
}     

/**
 * \ingroup  TopologyRelations
 * \brief Get the entities which own the specified entity.
 *
 * \description
 *   
 *    Given an entity on a mesh (e.g. vertex, element, etc) find the specified
 *    type of parents of that entity. The specified parent has to have a
 *    topological dimensionality greater than the child entity. (E.g, faces
 *    and elements can be parents of an edge, but vertices cannot.)
 *
 *    This type of information is sometimes called adjacencies, e.g. get
 *    the elements adjacent to the given vertex.
 *
 *    Entities cannot be of association type FC_AT_WHOLE_DATASET.
 *
 * \modifications
 *    - 07/09/04 WSK Created to replace fc_getElementIDsWithVertex with
 *      a more general purpose routine.
 */
FC_ReturnCode fc_getMeshEntityParents(
  FC_Mesh mesh,     /**< input - mesh */
  FC_AssociationType child_assoc,  /**< input - the child's entity type */
  int childID,      /**< input - the child's ID (e.g. vertex ID, etc) */
  FC_AssociationType parent_assoc,  /**< input - the type of mesh entities to
                                      return (must have topo dimensionality > 
                                      child_assoc. */
  int *numParent,    /**< output - the number of parents */
  int **parentIDs    /**< output - array of parent IDs */
) {
  FC_ReturnCode rc;
  int numEntity;
  
  // default return values
  if (numParent)
    *numParent = -1;
  if (parentIDs)
    *parentIDs = NULL;

  // check input
  if (!fc_isMeshValid(mesh) || !fc_isAssociationTypeValid(child_assoc) ||
      child_assoc == FC_AT_UNKNOWN || childID < 0 || 
      !fc_isAssociationTypeValid(parent_assoc) || parent_assoc == FC_AT_UNKNOWN
      || parent_assoc == FC_AT_WHOLE_DATASET || parent_assoc <= child_assoc ||
      !numParent || !parentIDs) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  rc = fc_getMeshNumEntity(mesh, child_assoc, &numEntity);
  if (rc != FC_SUCCESS || childID >= numEntity) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // log message
  fc_printfLogMessage("Getting the entities which own the specified entity.");

  // call a helper
  switch(child_assoc) {
  case FC_AT_VERTEX: 
    return _fc_getVertexParents(mesh, childID, parent_assoc, numParent, 
                                parentIDs);
  case FC_AT_EDGE:
    return _fc_getEdgeParents(mesh, childID, parent_assoc, numParent, 
                              parentIDs);
  case FC_AT_FACE:   
    return _fc_getFaceParents(mesh, childID, parent_assoc, numParent, 
                              parentIDs);
  case FC_AT_ELEMENT:
    *numParent = 1;
    *parentIDs = (int*)malloc(sizeof(int));
    if (*parentIDs == NULL) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      return FC_MEMORY_ERROR;
    }
    (*parentIDs)[0] = 0;
    return FC_SUCCESS;
  case FC_AT_WHOLE_MESH: // shouldn't reach here
  case FC_AT_WHOLE_DATASET: // shouldn't reach here
  case FC_AT_UNKNOWN: // shouldn't reach here
    fc_printfErrorMessage("Input Error: Wrong association type");
    return FC_INPUT_ERROR;
  }

  // should never get here
  return FC_ERROR;
}     

/**
 * \ingroup  TopologyRelations
 * \brief Change the type of mesh entities.
 *
 * \description
 *   
 *    Given an array of mesh entities, return an array of mesh entities
 *    of a different type that correspond to the same part of the mesh.
 *
 *    If the new association type is of higher dimensionality than the old
 *    (e.g. new subset will be on element and old was on vertices), the
 *    'doStrict' flag determines whether to keep new entities if not all of
 *    their child entities are in the old subset. 1 = be strict and only have
 *    the entity in the new subset if all of it's children are present in the
 *    old subset. 0 = be lenient and add all entities that have at least one
 *    child in the old subset. 
 *
 *    Entities cannot be of association type FC_AT_WHOLE_DATASET.
 *
 * \modifications
 *    - 08/24/04 WSK, created.
 */
FC_ReturnCode fc_changeMeshEntityType(
  FC_Mesh mesh,     /**< input - mesh */
  FC_AssociationType old_assoc,  /**< input - the input association type */
  int numOldID,     /**< input - the number of old entities */
  int* oldIDs,      /**< input - array of old entities' IDs */
  FC_AssociationType new_assoc,  /**< input - the output association type */
  int doStrict,     /**< input - flag for dropping incomplete entities */
  int *numNewID,    /**< output - the number of new entities */
  int **newIDs   /**< output - array of new entities' IDs */
) {
  FC_ReturnCode rc;
  int i, j;
  int numOldEntity, numTempID, *tempIDs;
  FC_SortedIntArray sia;
  int *mask;
  
  // default return values
  if (numNewID)
    *numNewID = -1;
  if (newIDs)
    *newIDs = NULL;

  // check input
  if (!fc_isMeshValid(mesh) || !fc_isAssociationTypeValid(old_assoc) || 
      old_assoc == FC_AT_UNKNOWN || old_assoc == FC_AT_WHOLE_DATASET ||
      numOldID < 0 || (numOldID > 0 && !oldIDs) || 
      !fc_isAssociationTypeValid(new_assoc) || 
      new_assoc == FC_AT_UNKNOWN || new_assoc == FC_AT_WHOLE_DATASET ||
      (doStrict != 0 && doStrict != 1) || !numNewID || !newIDs) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  rc = fc_getMeshNumEntity(mesh, old_assoc, &numOldEntity);
  for (i = 0; i < numOldID; i++) {
    if (oldIDs[i] < 0 || oldIDs[i] > numOldEntity) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
      return FC_INPUT_ERROR;
    }
  }

  // log message
  fc_printfLogMessage("Changing mesh entity type");

  // trivial case - empty old array
  if (numOldID == 0) {
    *numNewID = 0;
    return FC_SUCCESS;
  }
  
  // get list of all new association type (lenient)
  fc_initSortedIntArray(&sia);
  for (i = 0; i < numOldID; i++) {
    if (new_assoc < old_assoc)
      fc_getMeshEntityChildren(mesh, old_assoc, oldIDs[i], new_assoc,
                               &numTempID, &tempIDs);
    else if (new_assoc > old_assoc)
      fc_getMeshEntityParents(mesh, old_assoc, oldIDs[i], new_assoc,
                              &numTempID, &tempIDs);
    else {
      numTempID = 1;
      tempIDs = &oldIDs[i];
    }
    fc_addIntArrayToSortedIntArray(&sia, numTempID, tempIDs, 0);
    if (new_assoc != old_assoc)
      free(tempIDs);
  }
  
  // if strict, remove from members that don't have all of their children
  if (doStrict && new_assoc > old_assoc) {
    // make mask for easy lookup of old members
    mask = calloc(numOldEntity, sizeof(int));
    if (!mask)
      return FC_MEMORY_ERROR;
    for (i = 0; i < numOldID; i++)
      mask[oldIDs[i]] = 1;
    
    // walk sia and remove members
    j = 0;
    while (j < sia.numVal) {
      int notFound = 0;
      fc_getMeshEntityChildren(mesh, new_assoc, sia.vals[j], old_assoc,
                               &numTempID, &tempIDs);
      for (i = 0; i < numTempID; i++) {
        if (mask[tempIDs[i]] == 0) {
          notFound = 1;
          break;
        }
      }
      free(tempIDs);
      if (notFound) { // remove this node
	_fc_deleteEntryFromSortedIntArray(&sia, j);
      }
      else {
	j++;
      }
    }
    
    free(mask);
  }

  // convert sia to array
  return fc_convertSortedIntArrayToIntArray(&sia, numNewID, newIDs);
}     

//@}

/** \name Neighbors. */
//-------------------------------
//@{

/**
 * \ingroup  TopologyRelations
 * \brief Get the neighbors of a mesh subentity.
 *
 * \description
 *   
 *    Given an entity on a mesh (e.g. vertex, element, etc) find the neighbors
 *    of that entity. The returned entities will be of the same type as the
 *    initial entity. The order indicates the depth of the neighbors, i.e. an
 *    order of 1 will return the immediate neighbors of the entity, an order of
 *    2 will return the immediate neighbors and also the neighbors of the
 *    neighbors, etc.
 *
 *    The shared_dim parameter indicates the minimum dimensionality of the
 *    overlap (or shared parts) neighbors. In general, the shared_dim has to be
 *    less than the dimensionality of the entity and/or equal to 0 (vertices
 *    and point elements are special cases where shared_dim = dim = 0). Also,
 *    as the shared_dim increases, the condition for being neighbors becomes
 *    more stringent. For example, in a hex mesh, you can request element
 *    neighbors that share a face (share_dim = 2) or element neighbors that
 *    share an edge (share_dim = 1, this includes the face neighbors too) or
 *    element neighbors that share a vertex (share_dim = 0, this includes face
 *    and edge neighbors).
 *
 *    Entities cannot be of association type FC_AT_WHOLE_DATASET.
 *
 * \modifications
 *    - 07/09/04 WSK Created to replace all of the entity specific neighbor
 *      routines like fc_getElementVertexNeighbors() with a single general
 *      purpose routine.
 */
FC_ReturnCode fc_getMeshEntityNeighbors(
  FC_Mesh mesh,     /**< input - mesh */
  int entityID,     /**< input - mesh entity ID (e.g. vertex ID, etc) */
  FC_AssociationType assoc,  /**< input - the type of mesh entity */
  int shared_dim,   /**< input - the minimum dimensionality of shared part of
                       neighbors */
  int *numNbr,      /**< output - the number of neighbors */
  int **nbrIDs      /**< output - array of neighbor IDs */
) {
  FC_ReturnCode rc;
  int numEntity;
  
  // default return values
  if (numNbr)
    *numNbr = -1;
  if (nbrIDs)
    *nbrIDs = NULL;

  // check input
  if (!fc_isMeshValid(mesh) || entityID < 0 ||
      !fc_isAssociationTypeValid(assoc) || assoc == FC_AT_UNKNOWN || 
      assoc == FC_AT_WHOLE_DATASET ||
      shared_dim < 0 || shared_dim > 2 || !numNbr || !nbrIDs) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  rc = fc_getMeshNumEntity(mesh, assoc, &numEntity);
  if (rc != FC_SUCCESS || entityID >= numEntity) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // log message
  fc_printfLogMessage("Getting the neighbors of a mesh subentity.");

  // fill up list of neighbors by calling recursive functions
  switch (assoc) {
  case FC_AT_VERTEX:
    if (shared_dim != 0)
      return FC_INPUT_ERROR;
    else
      return _fc_getVertexNeighbors(mesh, entityID, numNbr, nbrIDs);
  case FC_AT_EDGE:
    if (shared_dim != 0)
      return FC_INPUT_ERROR;
    else
      return _fc_getEdgeNeighbors(mesh, entityID, numNbr, nbrIDs);
  case FC_AT_FACE:
    return _fc_getFaceNeighbors(mesh, entityID, shared_dim, numNbr, nbrIDs);
  case FC_AT_ELEMENT:
    return _fc_getElementNeighbors(mesh, entityID, shared_dim, numNbr, nbrIDs);
  case FC_AT_WHOLE_MESH:
    *numNbr = 0;
    return FC_SUCCESS;
  case FC_AT_WHOLE_DATASET: // shouldn't reach here
    return FC_ERROR;
  case FC_AT_UNKNOWN:
    fc_printfErrorMessage("Unknown association type");
    return FC_ERROR;
  }
  
  return FC_ERROR; // should never get here 
}     

/**
 * \ingroup  TopologyRelations
 * \brief Get the neighbors of a subset.
 *
 * \description
 *   
 *    This routine will return the neighbors of a subset.
 *    Neighbor type is determined by association type of the subset and
 *    the argument 'shared_dim' (see fc_getMeshEntityNeighbors for more on
 *    neighbors). Level indicates the depth of the neighbors where a
 *    level of 1 will return the the first layer of neighbors immediately
 *    surrounding the subset members, level 2 will return two layers of
 *    neighbors (neighbors of the neighbors), etc.
 *
 *    The shared_dim parameter indicates the minimum dimensionality of the
 *    overlap (or shared parts) neighbors. In general, the shared_dim has to be
 *    less than the dimensionality of the subset's entities and/or equal to 0
 *    (vertices and point elements are special cases where shared_dim = dim =
 *    0). Also, as the shared_dim increases, the condition for being neighbors
 *    becomes more stringent. For example, in a hex mesh, you can request
 *    element neighbors that share a face (share_dim = 2) or element neighbors
 *    that share an edge (share_dim = 1, this includes the face neighbors too)
 *    or element neighbors that share a vertex (share_dim = 0, this includes
 *    face and edge neighbors).
 *
 * \modifications
 *    - 07/023/04 WSK Created so that I can get rid of recursive getNeighbors
 *      calls.
 */
FC_ReturnCode fc_getSubsetNeighbors(
  FC_Subset subset, /**< input - subset */
  int level,        /**< input - depth of neighbors to return, 1 or greater */
  int shared_dim,   /**< input - the minimum dimensionality of shared part of
                       neighbors */
  FC_AssociationType* nbr_assoc,  /**< output - the type of mesh entity */
  int *numNbr,      /**< output - the number of neighbors */
  int **nbrIDs      /**< output - array of neighbor IDs */
) {
  FC_ReturnCode rc;
  int i;
  int maxNumMember, numMember, *memberIDs;
  int *checked;
  FC_Mesh mesh;
  FC_SortedIntArray members_list, nbrs_list;

  // default return values
  if (nbr_assoc)
    *nbr_assoc = FC_AT_UNKNOWN;
  if (numNbr)
    *numNbr = -1;
  if (nbrIDs)
    *nbrIDs = NULL;

  // check input
  if (!fc_isSubsetValid(subset) ||  level < 1 || shared_dim < 0 || 
      shared_dim > 2 || !nbr_assoc  || !numNbr || !nbrIDs) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // log message
  fc_printfLogMessage("Getting the neighbors of a subset.");

  // get info from subset
  rc = fc_getSubsetInfo(subset, &numMember, &maxNumMember, nbr_assoc);
  if (rc != FC_SUCCESS) {
    fc_printfErrorMessage("Error getting info from subset.");
    return rc;
  }
  rc = fc_getSubsetMembersAsArray(subset, &numMember, &memberIDs);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_getMeshFromSubset(subset, &mesh);
  if (rc != FC_SUCCESS)
    return rc;  
  fc_initSortedIntArray(&members_list);
  fc_addIntArrayToSortedIntArray(&members_list, numMember, memberIDs, 1);
  
  // Recursively expand the nbrs_list
  checked = calloc(maxNumMember, sizeof(int));
  if (checked == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  fc_initSortedIntArray(&nbrs_list);
  rc = _fc_growNeighborList(mesh, &members_list, *nbr_assoc, level, shared_dim,
                             checked, &nbrs_list);
  free(checked);
  fc_freeSortedIntArray(&members_list);
  if (rc != FC_SUCCESS) {
    *nbr_assoc = FC_AT_UNKNOWN;
    free(memberIDs);
    fc_printfErrorMessage("Error growing neighbor list.");
    return rc;
  }

  // make sure current members are not in list
  for (i = 0; i < numMember; i++)
    fc_deleteIntFromSortedIntArray(&nbrs_list, memberIDs[i]);
  free(memberIDs);

  // convert nbrs_list to array
  rc = fc_convertSortedIntArrayToIntArray(&nbrs_list, numNbr, nbrIDs);

  return rc;
}     

//@}

/**
 * \ingroup  PrivateTopologyRelations
 * \brief Get mesh entities which contain the given vertex.
 *
 * \description
 *  
 *    Specify a vertex and  mesh entity (e.g. element) and this routines
 *    returns IDs of that type of mesh entity which contain (are made from)
 *    the given vertex.
 *
 * \modifications 
 *    - 07/09/04 WSK created.
 */
FC_ReturnCode _fc_getVertexParents(
  FC_Mesh mesh,      /**< input - mesh */
  int vertexID,      /**< input - vertex number */
  FC_AssociationType parent_assoc, /**< input - association type of the parent */
  int *numParent,   /**< output - the number of parents */
  int **parentIDs   /**< output - array of parent IDs */
) {
  FC_ReturnCode rc;
  int i;
  int numVertex, *temp_IDs, whole_ids[1] = { 0 };
  int *numEntityPerVert, **entityParentsPerVert;

  // default return values
  if (numParent)
    *numParent = -1;
  if (parentIDs)
    *parentIDs = NULL;
  
  // check input
  if (!fc_isMeshValid(mesh) || vertexID < 0 || 
      !fc_isAssociationTypeValid(parent_assoc) || 
      parent_assoc == FC_AT_VERTEX || parent_assoc == FC_AT_UNKNOWN || 
      !numParent || !parentIDs) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // get mesh info
  rc = fc_getMeshNumVertex(mesh, &numVertex);
  if (rc != FC_SUCCESS || vertexID >= numVertex) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Getting mesh entities which contain the given vertex.");

  // get the number of parents & pointer to parent array to copy
  if (parent_assoc == FC_AT_WHOLE_MESH) {
    *numParent = 1;
    temp_IDs = whole_ids;
  }
  else {
    switch (parent_assoc) {
    case FC_AT_EDGE:
      rc = _fc_getMeshEdgeParentsOfVerticesPtr(mesh, &numEntityPerVert,
                                           &entityParentsPerVert);
      break;
    case FC_AT_FACE:
      rc = _fc_getMeshFaceParentsOfVerticesPtr(mesh, &numEntityPerVert,
                                           &entityParentsPerVert);
      break;
    case FC_AT_ELEMENT:
      rc = _fc_getMeshElementParentsOfVerticesPtr(mesh, &numEntityPerVert,
                                              &entityParentsPerVert);
      break;
    default:
      ; // do nothing, just here to quiet compiler
    }
    if (rc != FC_SUCCESS)
      return rc;
    // check to make sure entities of that type exist
    if (numEntityPerVert == NULL) { 
      *numParent = 0;
      *parentIDs = NULL;
      return FC_SUCCESS;
    }
    // finally, set values
    *numParent = numEntityPerVert[vertexID];
    temp_IDs = entityParentsPerVert[vertexID];
  }
  
  // make parent ID array
  if (*numParent > 0) {
    *parentIDs = (int*)malloc(*numParent*sizeof(int));
    if (*parentIDs == NULL) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      return FC_MEMORY_ERROR;
    }
    for (i = 0; i < *numParent; i++)
      (*parentIDs)[i] = temp_IDs[i];
  }
  
  return FC_SUCCESS;
}

/**
 * \ingroup  PrivateTopologyRelations
 * \brief Get mesh entities which contain the given edge.
 *
 * \description
 *  
 *    Specify an edge and  mesh entity (e.g. element) and this routines
 *    returns IDs of that type of mesh entity which contain (are made from)
 *    the given edge.
 *
 * \todo Make sure only calls edges/face stuff when it needs to. 
 *       Cleanup!
 *
 * \modifications 
 *    - 07/09/04 WSK created.
 */
FC_ReturnCode _fc_getEdgeParents(
  FC_Mesh mesh,      /**< input - mesh */
  int edgeID,        /**< input - edge id */
  FC_AssociationType parent_assoc, /**< input - association type of the parent */
  int *numParent,   /**< output - the number of parents */
  int **parentIDs   /**< output - array of parent IDs */
) {
  FC_ReturnCode rc;
  int i, j;
  int numEdge;
  FC_SortedIntArray sia;
  int* edgeConns;
  
  // default return values
  if (numParent)
    *numParent = -1;
  if (parentIDs)
    *parentIDs = NULL;
  
  // check input
  if (!fc_isMeshValid(mesh) || edgeID < 0 || 
      !fc_isAssociationTypeValid(parent_assoc) || 
      parent_assoc == FC_AT_VERTEX || parent_assoc == FC_AT_EDGE || 
      parent_assoc == FC_AT_UNKNOWN || !numParent || !parentIDs) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }  

  // get mesh info
  rc = fc_getMeshNumEdge(mesh, &numEdge);
  if (rc != FC_SUCCESS || edgeID >= numEdge) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  rc = fc_getMeshEdgeConnsPtr(mesh, &edgeConns);
  if (rc != FC_SUCCESS)
    return rc;

  // log message
  fc_printfLogMessage("Getting mesh entities which contain the given edge.");

  // type of parent
  if (parent_assoc == FC_AT_WHOLE_MESH) {
    *numParent = 1;
    *parentIDs = (int*)calloc(1, sizeof(int));
    if (*parentIDs == NULL) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      return FC_MEMORY_ERROR;
    }
    return FC_SUCCESS;
  }
  else if (parent_assoc == FC_AT_ELEMENT) {
    int *numElemPerEdge, **elemParentsPerEdge;
    rc = _fc_getMeshElementParentsOfEdgesPtr(mesh, &numElemPerEdge,
                                          &elemParentsPerEdge);
    if (rc != FC_SUCCESS)
      return rc;
    *numParent = numElemPerEdge[edgeID];
    *parentIDs = (int*)malloc(*numParent*sizeof(int));
    if (*parentIDs == NULL) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      return FC_MEMORY_ERROR;
    }
    for (i = 0; i < *numParent; i++)
      (*parentIDs)[i] = elemParentsPerEdge[edgeID][i];
    return FC_SUCCESS;
  }
  else if (parent_assoc == FC_AT_FACE) {
    int *numFacePerVert, **faceParentsPerVert;
    int vertID1, vertID2, numFace1, numFace2, *faceIDs1, *faceIDs2;
    fc_initSortedIntArray(&sia);
    rc = _fc_getMeshFaceParentsOfVerticesPtr(mesh, &numFacePerVert,
                                         &faceParentsPerVert);
    if (rc != FC_SUCCESS)
      return rc;
    if (numFacePerVert == 0) {
      *numParent = 0;
      *parentIDs = NULL;
      return FC_SUCCESS;
    }
    vertID1 = edgeConns[edgeID*2];
    vertID2 = edgeConns[edgeID*2+1];
    numFace1 = numFacePerVert[vertID1];
    numFace2 = numFacePerVert[vertID2];
    faceIDs1 = faceParentsPerVert[vertID1];
    faceIDs2 = faceParentsPerVert[vertID2];
    // WSK - thought about walking the two lists at the same time for
    // speedup, but that would require lists to be sorted, and since
    // these lists are usually small, speedup would probably be minimal. 
    for (i = 0; i < numFace1; i++) {
      for (j = 0; j < numFace2; j++) {
        if (faceIDs1[i] == faceIDs2[j]) {
          int ret = fc_addIntToSortedIntArray(&sia, faceIDs1[i]);
          if (ret < 0) {
            fc_printfErrorMessage("Failed to add int to sia");
            return ret;
          }
        }
      }
    }
    rc = fc_convertSortedIntArrayToIntArray(&sia, numParent, parentIDs);
    if (rc != FC_SUCCESS)
      return rc;
    return FC_SUCCESS;
  }
  else {
    fc_printfErrorMessage("shouldn't be able to get here."); 
    return FC_ERROR; // shouldn't be able to get here
  }
}

/**
 * \ingroup  PrivateTopologyRelations
 * \brief Get mesh entities which contain the given face.
 *
 * \description
 *  
 *    Specify a face and  mesh entity (e.g. element) and this routines
 *    returns IDs of that type of mesh entity which contain (are made from)
 *    the given face.
 *
 * \modifications 
 *    - 07/12/04 WSK created.
 */
FC_ReturnCode _fc_getFaceParents(
  FC_Mesh mesh,      /**< input - mesh */
  int faceID,        /**< input - face id */
  FC_AssociationType parent_assoc, /**< input - association type of the parent */
  int *numParent,   /**< output - the number of parents */
  int **parentIDs   /**< output - array of parent IDs */
) {
  FC_ReturnCode rc;
  int i;
  int numFace;
  
  // default return values
  if (numParent)
    *numParent = -1;
  if (parentIDs)
    *parentIDs = NULL;
  
  // check input
  if (!fc_isMeshValid(mesh) || faceID < 0 || 
      (parent_assoc != FC_AT_ELEMENT && parent_assoc != FC_AT_WHOLE_MESH) ||
      !numParent || !parentIDs) {
        fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR)); 
    return FC_INPUT_ERROR;
  }
  // get mesh info
  rc = fc_getMeshNumFace(mesh, &numFace);
  if (rc != FC_SUCCESS || faceID >= numFace) {
        fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // log message
  fc_printfLogMessage("Getting mesh entities which contain the given face.");

  // type of parent
  if (parent_assoc == FC_AT_WHOLE_MESH) {
    *numParent = 1;
    *parentIDs = (int*)calloc(1, sizeof(int));
    if (*parentIDs == NULL) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      return FC_MEMORY_ERROR;
    }
    return FC_SUCCESS;
  }
  else if (parent_assoc == FC_AT_ELEMENT) {
    int* numElemPerFace, **elemParentsPerFace;
    rc = _fc_getMeshElementParentsOfFacesPtr(mesh, &numElemPerFace,
                                          &elemParentsPerFace);
    if (rc != FC_SUCCESS)
      return rc;
    *numParent = numElemPerFace[faceID];
    if (*numParent > 0) {
      *parentIDs = (int*)malloc(*numParent*sizeof(int));
      if (*parentIDs == NULL) {
        fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
        return FC_MEMORY_ERROR;
      } 
      for (i = 0; i < *numParent; i++)
        (*parentIDs)[i] = elemParentsPerFace[faceID][i];
    }
    return FC_SUCCESS;
  }
  else {
        fc_printfErrorMessage("shouldn't be able to get here."); 
    return FC_ERROR; // shouldn't be able to get here
  }
}

/**
 * \ingroup  PrivateTopologyRelations
 * \brief Get vertices which are neighbors of this vertex.
 *
 * \description
 *  
 *    Returns the neighbor vertices of the given vertices. Vertices
 *    are considered to be neighbors if they are connected by
 *    an edge.
 *
 * \modifications 
 *    - Nancy Collins  Added this wrapper around steph's original routine to 
 *      use compute the neighbors and extract them in a subset.
 *    - 02/24/04 WSK, changed name and made more general purpose, like
 *      with _fc_getElementTopoNeighbors.
 *    - 07/15/04 WSK - Changed to not be recursive.
 */
FC_ReturnCode _fc_getVertexNeighbors(
  FC_Mesh mesh,      /**< input - the mesh */
  int vertexID,      /**< input - vertex ID */
  int* numNeigh,     /**< output - the number of neighbors */
  int** neighIDs    /**< output - the neighbor IDs */
) {
  FC_ReturnCode rc;
  int numVertex, *numVertNeighsViaEdge, **vertNeighsViaEdge;

  // default return
  if (numNeigh)
    *numNeigh = -1;
  if (neighIDs)
    *neighIDs = NULL;

  // check input
  if (!fc_isMeshValid(mesh) || vertexID < 0 || !numNeigh || !neighIDs) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  rc = fc_getMeshNumVertex(mesh, &numVertex);
  if (rc != FC_SUCCESS || vertexID >= numVertex) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // log message
  fc_printfLogMessage("Getting vertices which are neighbors of this vertex.");

  // copy array from Big Data struct
  rc = _fc_getMeshVertexNeighborsViaEdgePtr(mesh, &numVertNeighsViaEdge,
                                            &vertNeighsViaEdge);
  if (rc != FC_SUCCESS)
    return rc;
  // special case -- if edges don't exists
  if (numVertNeighsViaEdge == NULL) {
    *numNeigh = 0;
    neighIDs = NULL;
    return FC_SUCCESS;
  }
  *numNeigh = numVertNeighsViaEdge[vertexID];
  if (*numNeigh > 0) {
    *neighIDs = (int*)malloc((*numNeigh)*sizeof(int));
    if (*neighIDs == NULL) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      return FC_MEMORY_ERROR;
    }
    memcpy(*neighIDs, vertNeighsViaEdge[vertexID], *numNeigh*sizeof(int));
  }
    
  return FC_SUCCESS;
} 

/**
 * \ingroup  PrivateTopologyRelations
 * \brief Get edges which are neighbors of this edge.
 *
 * \description
 *  
 *    Returns the neighbor edges of the given edge. Edges
 *    are considered to be neighbors if they are connected by a vertex.
 *
 * \modifications 
 *    - 07/23/04 WSK, created.
 */
FC_ReturnCode _fc_getEdgeNeighbors(
  FC_Mesh mesh,     /**< input - the mesh */
  int edgeID,       /**< input - edge ID */
  int* numNeigh,    /**< output - the number of neighbors */
  int** neighIDs    /**< output - the neighbor IDs */
) {
  FC_ReturnCode rc;
  int i;
  int numEdge, *edgeConns, *numEdgePerVert, **edgeParentsPerVert;
  FC_SortedIntArray sia;

  // default return
  if (numNeigh)
    *numNeigh = -1;
  if (neighIDs)
    *neighIDs = NULL;

  // check input
  if (!fc_isMeshValid(mesh) || edgeID < 0 || !numNeigh || !neighIDs) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  rc = fc_getMeshNumEdge(mesh, &numEdge);
  if (rc != FC_SUCCESS || edgeID >= numEdge) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // log message
  fc_printfLogMessage("Getting edges which are neighbors of this edge.");

  // setup
  rc = fc_getMeshEdgeConnsPtr(mesh, &edgeConns);
  if (rc != FC_SUCCESS)
    return rc;
  rc = _fc_getMeshEdgeParentsOfVerticesPtr(mesh, &numEdgePerVert,
                                       &edgeParentsPerVert); 
  if (rc != FC_SUCCESS)
    return rc;  

  // Get all edges that own vertices from this edge
  fc_initSortedIntArray(&sia);
  for (i = 0; i < 2; i++) {
    int vertID = edgeConns[edgeID*2+i];
    int ret = fc_addIntArrayToSortedIntArray(&sia, numEdgePerVert[vertID],
					     edgeParentsPerVert[vertID], 0);
    if (ret < 0) {
      fc_printfErrorMessage("Failed to add int to sia");
      return ret;
    }
  }
  
  // remove this edge from list
  fc_deleteIntFromSortedIntArray(&sia, edgeID);

  // get array
  rc = fc_convertSortedIntArrayToIntArray(&sia, numNeigh, neighIDs);  
  
  return rc;
} 

/**
 * \ingroup  PrivateTopologyRelations
 * \brief Get faces which are neighbors of this faces.
 *
 * \description
 *  
 *    Returns the neighbor faces of the given face. The type of neighbor 
 *    depends on the input argument 'shared_dim' which is the minimum
 *    dimensionality of the overlap between the elements. For example,
 *    0 means that faces share at least a vertex while 1 means that
 *    faces must share an edge. The higher the shared_dim, the more
 *    more restrictive the neighbor condition.
 *
 * \modifications 
 *    - 07/23/04 WSK, created.
 */
FC_ReturnCode _fc_getFaceNeighbors(
  FC_Mesh mesh,     /**< input - the mesh */
  int faceID,       /**< input - face ID */
  int shared_dim,    /**< input - the minimum dimensionality of shared entity */
  int* numNeigh,    /**< output - the number of neighbors */
  int** neighIDs    /**< output - the neighbor IDs */
) {
  FC_ReturnCode rc;
  int i, j, k;
  int numFace, maxNumVertPerFace, *numVertPerFace, *faceConns;
  int *numFacePerVert, **faceParentsPerVert;
  //_FC_Vertex *vertices, *vertex_p;
  FC_SortedIntArray sia;

  // default return
  if (numNeigh)
    *numNeigh = -1;
  if (neighIDs)
    *neighIDs = NULL;

  // check input
  if (!fc_isMeshValid(mesh) || faceID < 0 || shared_dim < 0 || shared_dim > 1
      || !numNeigh || !neighIDs) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR)); 
    return FC_INPUT_ERROR;
  }
  rc = fc_getMeshNumFace(mesh, &numFace);
  if (rc != FC_SUCCESS || faceID >= numFace) {
        fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR)); 
    return FC_INPUT_ERROR;
  }

  // log message
  fc_printfLogMessage("Getting faces which are neighbors of this face.");

  // setup
  rc = fc_getMeshFaceConnsPtr(mesh, &numVertPerFace, &maxNumVertPerFace,
                              &faceConns);
  if (rc != FC_SUCCESS)
    return rc;
  rc = _fc_getMeshFaceParentsOfVerticesPtr(mesh, &numFacePerVert,
                                      &faceParentsPerVert);
  if (rc != FC_SUCCESS)
    return rc;

  // Vertex neighbors: Get all faces that own vertices from this face
  fc_initSortedIntArray(&sia);
  for (i = 0; i < numVertPerFace[faceID]; i++) {
    int vertID = faceConns[faceID*maxNumVertPerFace+i];
    int ret = fc_addIntArrayToSortedIntArray(&sia, numFacePerVert[vertID],
					     faceParentsPerVert[vertID], 0);
    if (ret < 0) {
      fc_printfErrorMessage("Failed to add ID to sia");
      return ret;
    }
  }
  
  // remove this face from list
  fc_deleteIntFromSortedIntArray(&sia, faceID);

  // Edge neighbors - remove faces that don't have 2 verts from this face
  if (shared_dim == 1) {
    k = 0;
    while (k < sia.numVal) {
      int currentID = sia.vals[k];
      int found = 0;
      for (i = 0; i < numVertPerFace[faceID]; i++) {
        for (j = 0; j < numVertPerFace[currentID]; j++) {
          if (faceConns[faceID*maxNumVertPerFace+i] == 
              faceConns[currentID*maxNumVertPerFace+j]) {
            found++;
            break;
          }
        }
        if (found == 2)
          break;
      }
      if (found < 2) {
        fc_deleteIntFromSortedIntArray(&sia, currentID);
      }
      else {
	k++;
      }
    }
  }
  
  // get array
  rc = fc_convertSortedIntArrayToIntArray(&sia, numNeigh, neighIDs);
    
  return rc;
} 

/**
 * \ingroup  PrivateTopologyRelations
 * \brief  Get the neighbors of an element
 *
 * \description
 *
 *    Returns the elements that are neighbors of the given element.  The type
 *    of neighbor depends on the input argument 'shared_dim' which is the
 *    minimum dimensionality of the overlap between the elements.  For example,
 *    0 means that elements that share a vertex, edge or face are neighbors,
 *    while 2 means that only elements that share a face are neighbors. The
 *    higher the shared_dim, the more restrictive the neighbor condition.
 *
 * \modifications 
 *    - 04/24/03 WSK, Create. A more general purpose function to replace 
 *      _fc_getEdgeNbrs and _fc_getVertexNbrs since they do almost
 *      exactly the same.
 */
FC_ReturnCode _fc_getElementNeighbors(
  FC_Mesh mesh,      /**< input - the mesh */
  int elemID,        /**< input - element number */
  int shared_dim,    /**< input - the minimum dimensionality of shared entity */
  int* numNeigh,     /**< output - the number of neighbors */
  int** neighIDs    /**< output - the neighbor IDs */
) {
  FC_ReturnCode rc;
  int *numElemNeighsViaEntity, **elemNeighsViaEntity;

  // default return
  if (numNeigh)
    *numNeigh = -1;
  if (neighIDs)
    *neighIDs = NULL;

  // check input
  if (!fc_isMeshValid(mesh) || elemID < 0 || shared_dim < 0 || shared_dim > 2
      || !numNeigh || !neighIDs) {
        fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // log message
  fc_printfLogMessage("Getting the neighbors of an element.");

  // Get entity specific info
  switch (shared_dim) {
  case 0: 
    rc = _fc_getMeshElementNeighborsViaVertexPtr(mesh, &numElemNeighsViaEntity, 
                                                 &elemNeighsViaEntity);
    break;
  case 1:
    rc = _fc_getMeshElementNeighborsViaEdgePtr(mesh, &numElemNeighsViaEntity, 
                                               &elemNeighsViaEntity);
    break;
  case 2:
    rc = _fc_getMeshElementNeighborsViaFacePtr(mesh, &numElemNeighsViaEntity, 
                                               &elemNeighsViaEntity);
    break;
  }
  if (rc != FC_SUCCESS)
    return rc;

  // Copy the ids
  if (numElemNeighsViaEntity == NULL) {
    *numNeigh = 0;
    *neighIDs = NULL;
    return FC_SUCCESS;
  }
  *numNeigh = numElemNeighsViaEntity[elemID];
  if (*numNeigh > 0) {
    *neighIDs = (int*)malloc(*numNeigh*sizeof(int));
    if (*neighIDs == NULL) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      return FC_MEMORY_ERROR;
    }
    memcpy(*neighIDs, elemNeighsViaEntity[elemID], *numNeigh*sizeof(int));
  }

  return FC_SUCCESS;
} 

/**
 * \ingroup  PrivateTopologyRelations
 * \brief Get the neighbors of a group of mesh entities
 *
 * \description
 *   
 *    This is a helper function for fc_getSubsetNeighbors().
 *    It takes a list of entities to get neighbors from and adds them
 *    to the list of neighbors. It is meant to be called recursively with
 *    the first call providing the seed IDs and a null list to add
 *    to. In addition, a mask for the association type is provided
 *    with all zero entries (can prevent adding specific IDs by setting
 *    their mask entry to one). 
 *
 *    Entities cannot be of association type FC_AT_WHOLE_DATASET.
 *
 * \modifications
 *    - 07/023/04 WSK Created.
 */
FC_ReturnCode _fc_growNeighborList(
  FC_Mesh mesh,       /**< input - the mesh the entities are in */
  FC_SortedIntArray* srcIDs,   /**< input - List of IDs to get neighbors from */
  FC_AssociationType assoc, /**< input - the type of entity */
  int level,          /**< input - depth of recursion (1 = no recursion) */
  int shared_dim,     /**< input - neighbor shared dimensionality */
  int* checked,       /**< input - mask allocated by caller */
  FC_SortedIntArray* keepIDs   /**< input/output - list of neighbor IDs */
) {
  FC_ReturnCode rc;
  int i, j;
  FC_SortedIntArray newIDs;
  int numNeigh, *neighIDs;
  
  // input no checking because it's a helper function
 
  // log message
  fc_printfLogMessage("Getting the neighbors of a group of mesh entities.");

  // go through srcIDs & get unexamined neighbors
  fc_initSortedIntArray(&newIDs);
  for (j = 0; j < srcIDs->numVal; j++) {
    int currentID = srcIDs->vals[j];

    if (checked[currentID])
      continue;

    // check it's neighbors
    rc = fc_getMeshEntityNeighbors(mesh, currentID, assoc, shared_dim,
				   &numNeigh, &neighIDs);
    if (rc != FC_SUCCESS)
      return rc;
    for (i = 0; i < numNeigh; i++) {
        if (!checked[neighIDs[i]]) {
          int ret = fc_addIntToSortedIntArray(&newIDs, neighIDs[i]);
          if (ret < 0) {
            fc_printfErrorMessage("Failed to add ID to list");
            return ret;
          }
        }
    }
    free(neighIDs);
    checked[currentID] = 1;
  }
  
  // add the new neighbors to our keep list
  rc = fc_addIntArrayToSortedIntArray(keepIDs, newIDs.numVal,
				      newIDs.vals, 1);
  if (rc < FC_SUCCESS)
    return rc;


  // call recursively
  if (level == 1 || newIDs.numVal == 0) {
    fc_freeSortedIntArray(&newIDs);
    return FC_SUCCESS;
  }
  else {
    rc = _fc_growNeighborList(mesh, &newIDs, assoc, level-1, 
                               shared_dim, checked, keepIDs);
    fc_freeSortedIntArray(&newIDs);
    return rc;
  }
}

/** \name Create connected components. */
//-------------------------------
//@{

/** 
 * \ingroup  TopologyRelations
 * \brief Created separate connected components
 *
 * \description 
 *
 *    This function takes a subset and separates into separate connected
 *    regions. It returns the number of regions and a new subset for each
 *    region. The new subsets will be named as the name of the original
 *    subset, postpended with "_Seg#" where # will be the index into
 *    the array of returned subsets.
 *
 *    Mesh subentities (e.g. elements, vertices, etc) are considered to
 *    be connected if they are neighbors. The type of neighborship relation
 *    is specified by the minimum shared dim (see fc_getMeshEntityNeighbors).
 *    The most common value of shared_dim is 0.
 *
 * \modifications 
 *   - 02/06/04 WSK, created.
 *   - 06/29/04 RM,  combined _fc_segment() into fc_segment(), changed 
 *       interface to take subset as input and return subsets
 */
FC_ReturnCode fc_segment(
  FC_Subset subset,       /**< input - subset */
  int shared_dim,         /**< input - the minimum dimensionality of shared 
                             part of neighbors */
  int *numSubset,         /**< output - number of new subsets */
  FC_Subset** newSubsets  /**< output - array of segmented regions */
) {
  int i;
  int ret;
  FC_ReturnCode rc;
  int numSegment = 0;
  int numMember, *memberIDs;
  FC_Subset *segments = NULL;
  FC_SortedIntArray orig_list;  // the starting list, will end empty
  FC_AssociationType assoc;
  _FC_SubSlot *subsetSlot;
  char* temp_name, *subsetName;
  FC_Mesh mesh;
  int key, numID, maxNumID, *IDs; // not using an sia for ids - don't need sort
  int numWorking, maxNumWorking, *workingIDs; // also don't need sort
  void* temp;

  //---default return values
  if (numSubset)
    *numSubset = -1;
  if (newSubsets)
    *newSubsets = NULL;

  // check input
  subsetSlot = _fc_getSubSlot(subset);
  if (!subsetSlot || shared_dim < 0 || shared_dim > 2 || !numSubset || 
      !newSubsets) { 
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  rc = fc_getSubsetInfo(subset, NULL, NULL, &assoc); 
  if (rc != FC_SUCCESS)
    return FC_ERROR;
  if ( ( (assoc == FC_AT_VERTEX || assoc == FC_AT_EDGE) && shared_dim > 0) ||
       (assoc == FC_AT_FACE && shared_dim > 1) ) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // log message
  fc_printfLogMessage("Segmenting subset '%s'", subsetSlot->header.name);

  // Get more info from subset
  rc = fc_getMeshFromSubset(subset, &mesh);
  if (rc != FC_SUCCESS)
    return FC_ERROR;
  rc = fc_getSubsetMembersAsArray(subset, &numMember, &memberIDs);
  if (rc != FC_SUCCESS)
    return FC_ERROR;
  subsetName = subsetSlot->header.name;

  // trivial case, if no members, not an error, return success
  if (numMember == 0) {
    *numSubset = 0;
    return FC_SUCCESS;
  }

  // convert memberIDs to orig list
  rc = fc_convertIntArrayToSortedIntArray(numMember, memberIDs, 1, &orig_list);
  if (rc < FC_SUCCESS)
    return rc;

  // setup space for names of segments
  temp_name = malloc((strlen(subsetName)+20)*sizeof(char));
  if (temp_name == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  
  //---segment the linked list
  numSegment = 0;
  while (orig_list.numVal > 0) { // each pass will create a new subset

    // FIX change all to back 'cause more efficient

    // initialize working list & segment list
    numWorking = 0;
    maxNumWorking = 0;
    workingIDs = NULL;
    numID = 0;
    maxNumID = 0;
    IDs = NULL;

    // remove starting ID (seed) from orig list ...
    fc_getSortedIntArrayFront(&orig_list, &key);
    fc_popSortedIntArrayFront(&orig_list);

    // ... & add to working & segment lists
    _fc_expandIntArray(&maxNumWorking, &workingIDs);
    _fc_expandIntArray(&maxNumID, &IDs);
    workingIDs[0] = key;
    numWorking++;
    IDs[0] = key;
    numID++;

    // "grow" connected region from seed by checking neighbors recursively
    while (numWorking > 0) {
      int numNeighbor, *neighborIDs;

      // can stop early if orig_list is empty
      if (orig_list.numVal == 0) 
        break;

      // remove 1st entity from working list
      key = workingIDs[numWorking-1];
      numWorking--;

      // Move neighbors of that entity that are in the orig
      // list to the segments list & also the working list
      rc = fc_getMeshEntityNeighbors(mesh, key, assoc, shared_dim,
                                     &numNeighbor, &neighborIDs);
      if (rc != FC_SUCCESS)
        return rc;
      for (i = 0; i < numNeighbor; i++) {
        ret = fc_deleteIntFromSortedIntArray(&orig_list, neighborIDs[i]);
        if (ret > 0) { // found one
	  if (numWorking == maxNumWorking) {
	    rc = _fc_expandIntArray(&maxNumWorking, &workingIDs);
	    if (rc != FC_SUCCESS)
	      return rc;
	  }
	  if (numID == maxNumID) {
	    rc = _fc_expandIntArray(&maxNumID, &IDs);
	    rc = _fc_expandIntArray(&maxNumWorking, &workingIDs);
	    if (rc != FC_SUCCESS)
	      return rc;
	  }
	  workingIDs[numWorking] = neighborIDs[i];
	  numWorking++;
	  IDs[numID] = neighborIDs[i];
	  numID++;
        }
      }
      free(neighborIDs);
    }

    // create subset from list of IDs
    temp = realloc(segments, (numSegment+1)*sizeof(FC_Subset));
    if (temp == NULL) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      return FC_MEMORY_ERROR;
    }
    segments = temp;
    sprintf(temp_name, "%s_Seg%d", subsetName, numSegment);
    rc = fc_createSubset(mesh, temp_name, assoc, &segments[numSegment]);
    if (rc != FC_SUCCESS) {
      fc_printfErrorMessage("cannot create newSubset[%d]", numSegment);
      return FC_ERROR; 
    }
    rc = fc_addArrayMembersToSubset(segments[numSegment], numID, IDs);
    if (rc != FC_SUCCESS) {
      fc_printfErrorMessage("cannot add members to newSubset[%d]", numSegment);
      return FC_ERROR; 
    }
    free(workingIDs);
    free(IDs);
    numSegment++;
  }

  // cleanup (orig_list is empty)
  free(temp_name);
  fc_freeSortedIntArray(&orig_list);

  // return values
  *numSubset = numSegment;
  *newSubsets = segments;

  return FC_SUCCESS;
}

//@}

/** \name Skins. */
//-------------------------------
//@{

/** 
 * \ingroup  TopologyRelations
 * \brief Create a skin of the subset.
 *
 * \description 
 *
 *    This function creates a subset which contains the "skin" of the original
 *    subset. The  topo dimensionality of the skin will always be 1 less than 
 *    the topo dimensionality of the original subset's members.
 *    (e.g. the skin of a subset containing 3D elements will be faces, the
 *    skin of an edge subset will be vertices).
 *
 *    This may change: for now, it is an error to try to create a skin
 *    of an empty subset.
 *
 *    The skin of a subset is found by iterating over the "skins" of the
 *    original mesh entities and determining which skin elements are only seen
 *    once. The list of unique skins of the original mesh entities will be the
 *    skin of the entire subset.
 *
 *    Be warned that behavior on subsets associated with entities that
 *    are not elements may not be what you expect. For example, if you
 *    have a subset which consists of the six faces of a hex element, the
 *    generated skin will be empty because the edges of the faces
 *    all belong to more than one face.
 *
 * \todo Decide what to do if no members in original subset. Error or not?
 *      If not, return null handle or empty subset? 
 *
 * \modifications 
 *   - 09/03/04 WSK, created.
 */
FC_ReturnCode fc_getSubsetSkin(
  FC_Subset subset,    /**< input - subset */
  FC_Subset* skin      /**< output - array of segmented regions */
) {
  FC_ReturnCode rc;
  int i, j;
  FC_AssociationType assoc_orig, assoc_skin;
  char* temp_name, *subset_name;
  int numMember, *memberIDs, maxNumSkin, *counts, numChild, *childIDs;
  int topodim_orig;
  FC_Mesh mesh;
  FC_ElementType elemType;

  //---default return values
  if (skin)
    *skin = FC_NULL_SUBSET;

  // check input
  if (!fc_isSubsetValid(subset) || !skin) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // log message
  fc_printfLogMessage("Creating new subset which is the skin of this subset");

  // Get info from subset and a little more checking
  rc = fc_getSubsetNumMember(subset, &numMember);
  if (rc != FC_SUCCESS)
    return FC_ERROR;
  rc = fc_getSubsetInfo(subset, NULL, NULL, &assoc_orig); 
  if (rc != FC_SUCCESS)
    return FC_ERROR;
  rc = fc_getMeshFromSubset(subset, &mesh);
  if (rc != FC_SUCCESS)
    return FC_ERROR;
  if (numMember < 1) {
    fc_printfErrorMessage("Can't create skin of empty subset");
    return FC_ERROR;
  }
  if (!fc_isAssociationTypeValid(assoc_orig) || assoc_orig == FC_AT_UNKNOWN ||
      assoc_orig == FC_AT_VERTEX) {
    fc_printfErrorMessage("Can't handle subset with association type %s",
                          fc_getAssociationTypeText(assoc_orig));
    return FC_ERROR;
  }

  // Get members -- If subset is whole, get all elements
  if (assoc_orig == FC_AT_WHOLE_MESH) {
    assoc_orig = FC_AT_ELEMENT;
    rc = fc_getMeshNumElement(mesh, &numMember);
    if (rc != FC_SUCCESS)
      return FC_ERROR;
    memberIDs = malloc(numMember*sizeof(int));
    if (!memberIDs) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      return FC_MEMORY_ERROR;
    }
    for (i = 0; i < numMember; i++)
      memberIDs[i] = i;
  }
  else {
    rc = fc_getSubsetMembersAsArray(subset, &numMember, &memberIDs);
    if (rc != FC_SUCCESS)
      return FC_ERROR;
  }

  // determine topodim of the members
  rc = fc_getMeshEntityElementType(mesh, assoc_orig, 0, &elemType);
  if (rc != FC_SUCCESS)
    return FC_ERROR;
  topodim_orig = fc_getElementTypeTopoDim(elemType);

  // count children of the members
  switch (topodim_orig) {
  case 3:  assoc_skin = FC_AT_FACE;    break;
  case 2:  assoc_skin = FC_AT_EDGE;    break;
  case 1:  assoc_skin = FC_AT_VERTEX;  break;
  default:
    fc_printfErrorMessage("Invalid topodim '%d'", topodim_orig);
    free(memberIDs);
    return FC_ERROR;
  }
  rc = fc_getMeshNumEntity(mesh, assoc_skin, &maxNumSkin);
  if (rc != FC_SUCCESS)
    return rc;
  counts = calloc(maxNumSkin, sizeof(int));
  if (!counts) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  for (i = 0; i < numMember; i++) {
    fc_getMeshEntityChildren(mesh, assoc_orig, memberIDs[i], assoc_skin,
                             &numChild, &childIDs);
    for (j = 0; j < numChild; j++)
      counts[childIDs[j]]++;
    free(childIDs);
  }
  free(memberIDs);

  // create new subset
  rc = fc_getSubsetName(subset, &subset_name);
  if (rc != FC_SUCCESS)
    return FC_ERROR;
  temp_name = malloc((strlen(subset_name)+10)*sizeof(char));
  sprintf(temp_name, "%s_skin", subset_name);
  free(subset_name);
  rc = fc_createSubset(mesh, temp_name, assoc_skin, skin);
  free(temp_name);
  if (rc != FC_SUCCESS) {
    fc_printfErrorMessage("Error creating subset.");
    return FC_ERROR;
  }

  // stuff in new members
  for (i = maxNumSkin-1; i >= 0; i--) {
    if (counts[i] == 1) {
      fc_addMemberToSubset(*skin, i);
    }
  }

  // cleanup
  free(counts);

  return FC_SUCCESS;
}


/** 
 * \ingroup  TopologyRelations
 * \brief Create skin of the mesh.
 *
 * \description 
 *
 *    This function creates a subset which contains the "skin" of the mesh.
 *    The  topo dimensionality of the skin will always be 1 less than 
 *    the topo dimensionality of the mesh's elements
 *    (e.g. the skin of a hex mesh will be faces, and the
 *    skin of an edge mesh will be verticies).
 *
 *    This routine is similar to fc_getSubsetSkin.
 *
 * \todo ?Let user name skin.
 *
 * \modifications 
 *   - 09/03/04 WSK, created.
 */
FC_ReturnCode fc_getMeshSkin(
  FC_Mesh mesh,        /**< input - mesh */
  FC_Subset* skin      /**< output - subset with skin of the mesh */
) {
  FC_ReturnCode rc;
  FC_Subset temp_subset;

  //---default return values
  if (skin)
    *skin = FC_NULL_SUBSET;

  // check input
  if (!fc_isMeshValid(mesh) || !skin) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // log message
  fc_printfLogMessage("Creating new subset which is the skin of the mesh");

  // create a temp subset which is the whole mesh & call fc_getSubsetSkin()
  rc = fc_createSubset(mesh, "fc_getMeshSkin temp subset", FC_AT_WHOLE_MESH,
                       &temp_subset);
  if (rc != FC_SUCCESS) {
    fc_deleteSubset(temp_subset);
    return rc;
  }
  fc_addMemberToSubset(temp_subset, 0);

  // do it
  rc = fc_getSubsetSkin(temp_subset, skin);
  fc_deleteSubset(temp_subset);
  return rc;
}

//@}

// /** 
// * \ingroup  TopologyRelations
// * \brief given a subset, skins it, and if the skin is a closed surface,
// *        returns its genus
// *
// * \description 
// * 
// *        given a subset, skins it, and if the skin is a closed surface,
// *        returns its genus. if the skin is not a closed surface(determined
// *        by segmenting the skin) it returns with an error
// *
// * NOTE: still workign on this. This is not guarenteed to be
// *       anything meanignful right now.
// *
// * \todo
// *   - untested
// *   - shoudl they jsut pass in the skin ?
// *   - check genus is an int. should we return the euler characteristic 
// *     instead ?
// *   - check what this means - getting some large genus values for
// *     some of the more deformed screw breaking regions.....
// *
// * \modifications 
// *   - 09/23/05 ACG, created.
// *   - 10/04/05 ACG, commenting out because way segmenting algorithm
// *                    works, it does not actaully guarentee that
// *                    the skin is a valid surface (can have contact
// *                    just at an edge for many elements)
// */
//
//FC_ReturnCode fc_getGenusOfSubsetSkin(
//  FC_Subset subset, /**< input - subset */
//  int* genus /** output - genus */
//){
/*
  FC_ReturnCode rc;
  FC_Subset skin, *segments,junk;
  int numSegments, numfaces, numedges, numverticies;
  int eulerchar;
  int i;

// default returns
  if (genus)
    *genus = -1;

  //using shared dim of 0 by default - check if want that.....
  // Test input 
  if (!fc_isSubsetValid(subset) ||  !genus){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }  

  // log message
  fc_printfLogMessage("genus of subset skin");

  rc = fc_getSubsetSkin(subset, &skin);
  if (rc != FC_SUCCESS){
    return rc;
  }

  rc = fc_segment(skin,0,&numSegments,&segments);
  if (rc != FC_SUCCESS){
    fc_deleteSubset(skin);
    return rc;
  }

  for (i = 0; i <numSegments; i++){
    fc_deleteSubset(segments[i]);
  }
  free(segments);

  if (numSegments != 1){
    fc_printfErrorMessage("Not a closed surface");
    fc_deleteSubset(skin);
    return FC_INPUT_ERROR;
  }

  //now get quantities
  rc = fc_getSubsetNumMember(skin,&numfaces);
  if (rc != FC_SUCCESS){
    fc_deleteSubset(skin);
    return rc;
  }

  rc = fc_copySubsetWithNewAssociation(skin,"junk",FC_AT_EDGE,0,&junk);
  if (rc != FC_SUCCESS){
    fc_deleteSubset(skin);
    return rc;
  }

  rc = fc_getSubsetNumMember(junk,&numedges);
  if (rc != FC_SUCCESS){
    fc_deleteSubset(skin);
    fc_deleteSubset(junk);
    return rc;
  }

  fc_deleteSubset(junk);


  rc = fc_copySubsetWithNewAssociation(skin,"junk",FC_AT_VERTEX,0,&junk);
  if (rc != FC_SUCCESS){
    fc_deleteSubset(skin);
    return rc;
  }

  rc = fc_getSubsetNumMember(junk,&numverticies);
  if (rc != FC_SUCCESS){
    fc_deleteSubset(skin);
    fc_deleteSubset(junk);
    return rc;
  }

  fc_deleteSubset(junk);
  fc_deleteSubset(skin);


  eulerchar = numverticies - numedges + numfaces;
  if (((eulerchar-2) % 2) != 0){
    fc_printfErrorMessage("Not integer genus");
    return FC_ERROR;
  }

  *genus = -((numverticies - numedges + numfaces)-2)/2; 

  return FC_SUCCESS;
}
*/
