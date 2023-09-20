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
 * \file mesh.c
 * \brief Implementation for \ref Mesh module.
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/modules/mesh.c,v $
 * $Revision: 1.89 $ 
 * $Date: 2006/10/19 03:14:50 $
 */

// C library includes
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// fc library dependencies
#include "base.h"
#include "storage.h"
#include "storageP.h"
#include "library.h"
#include "dataset.h"
#include "sequence.h"
#include "subset.h"
#include "variable.h"
#include "tableP.h"
#include "datasetP.h"
#include "sequenceP.h"
#include "subsetP.h"
#include "variableP.h"
#include "fileioP.h"
#include "topo.h"

// this module
#include "mesh.h"
#include "meshP.h"
#include "util.h"

/**
 * \addtogroup  Mesh
 * \brief  Operations on meshes.
 *
 * \description
 *
 *    (For an explanation of general data manipulations, see 
 *    \ref DataInterface.)
 */

/**
 * \ingroup  Mesh
 * \defgroup  PrivateMesh (Private)
 */

/** 
 * \ingroup PrivateMesh
 * \brief The number of allocated slots in the mesh table.
 */
static int meshTableSize = 0;
/** 
 * \ingroup PrivateMesh
 * \brief The mesh table.
 */
static _FC_MeshSlot **meshTable = NULL;
/**
 * \ingroup PrivateMesh
 * \brief The list of unused slots in the mesh table.
 */
static FC_SortedIntArray meshOpenSlots = { 0, 0, 0 };

/** \name Create new mesh from scratch. */
//-------------------------------
//@{

/**
 * \ingroup  Mesh
 * \brief  Create a new mesh
 *
 * \description
 * 
 *    Create a new, in-memory, mesh.
 *
 * \modifications  
 *    - Nancy Collins, Created.
 */
FC_ReturnCode fc_createMesh(
  FC_Dataset dataset,  /**< input - dataset handle to associate mesh with */
  char *meshname,      /**< input - new mesh name */
  FC_Mesh *mesh        /**< output - mesh handle */
) {
  int *tempslotIDs;
  _FC_DsSlot* dsSlot;
  _FC_MeshSlot* meshSlot;
  
  // default returns
  if (mesh)
    *mesh = FC_NULL_MESH;
  
  // check input
  dsSlot = _fc_getDsSlot(dataset);
  if (dsSlot == NULL || !meshname || !mesh) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Creating new mesh '%s'", meshname);
  
  // get an open slot
  meshSlot = _fc_getNewMeshSlot();
  if (meshSlot == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  
  // table header information
  _FC_GET_HANDLE(*mesh, meshSlot);
  _fc_setSlotHeaderName(&(meshSlot->header), meshname);
  
  // set back and forth references
  meshSlot->ds = dataset;
  
  // update owning dataset with info about this new mesh
  tempslotIDs = realloc(dsSlot->meshIDs, sizeof(int) * (dsSlot->numMesh+1));
  if (tempslotIDs == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  dsSlot->meshIDs = tempslotIDs;
  dsSlot->meshIDs[dsSlot->numMesh] = meshSlot->header.slotID;
  dsSlot->numMesh++;
  
  return FC_SUCCESS;
}

/**
 * \ingroup  Mesh
 * \brief  Set the mesh coordinates by copy.
 *
 * \description
 *  
 *    Set the coordinate values for the vertices in a mesh.  The dim parameter
 *    is the number of spatial dimensions a vertex exists in.  The coords 
 *    array should contain dim*numVertex floating point numbers. The order
 *    is interleaved, that is your give all coordinates values for each
 *    vertex (e.g. XYZXYZ ...).
 *
 *    This routine makes a copy of the vertices for the library,
 *    so the user remains responsible for freeing the original array.
 *    If you would like to give the coordinates array to the mesh instead
 *    of just copying them, use fc_setMeshCoordsPtr().
 *
 *    You cannot set coords if they already exist.
 *
 * \modifications  
 *    - Nancy Collins Created.
 *    - 12/09/03 WSK Changed so that it creates a coordinate variable
 *    - 09/13/04 WSK Changed so you cannot reset coords
 */
FC_ReturnCode fc_setMeshCoords(
  FC_Mesh mesh,       /**< input - mesh handle */
  int dim,            /**< input - number of dimensions in 1 vertex */
  int numVertex,      /**< input - vertex count */
  double *coords      /**< input - array of vertex data */
) {
  _FC_MeshSlot *meshSlot;
  
  // NOTE!: if you change stuff here, also change fc_setSetMeshCoordsPtr()

  // check input
  meshSlot = _fc_getMeshSlot(mesh);
  if (meshSlot == NULL || dim < 1 || dim > 3 || numVertex < 1 || 
      coords == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // a little more checking -- can't overwrite coords
  if (meshSlot->coords != NULL || meshSlot->header.committed == 1) {
    fc_printfErrorMessage("Coords already exist on mesh '%s'",
                          meshSlot->header.name);
    return FC_ERROR;
  }

  // log message
  fc_printfLogMessage("Setting vertex coords on mesh '%s'", 
                      meshSlot->header.name);

  // copy meta data
  meshSlot->numVertex = numVertex;
  meshSlot->dim = dim;
  
  // big data
  meshSlot->coords = malloc(sizeof(double) * numVertex * dim);
  if (meshSlot->coords == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  memcpy(meshSlot->coords, coords, sizeof(double)*numVertex*dim);
  
  return FC_SUCCESS;
}

/**
 * \ingroup  Mesh
 * \brief  Set the pointer to the mesh coordinates.
 *
 * \description
 *  
 *    Set the coordinate values for the vertices in a mesh.  The dim parameter
 *    is the number of spatial dimensions a vertex exists in.  The coords 
 *    array should contain dim*numVertex floating point numbers. The order
 *    is interleaved, that is your give all coordinates values for each
 *    vertex (e.g. XYZXYZ ...).
 *
 *    This routine passes the coords array directly to the library,
 *    so the user no longer has responsibility of the original buffer.
 *    If you would like to keep the coords array, you should use 
 *    fc_setMeshCoords() instead which makes a copy of the coords array.
 *
 *    You cannot set coords if they already exist.
 *
 * \modifications  
 *    - 10/01/04 WSK, created.
 */
FC_ReturnCode fc_setMeshCoordsPtr(
  FC_Mesh mesh,       /**< input - mesh handle */
  int dim,            /**< input - number of dimensions in 1 vertex */
  int numVertex,      /**< input - vertex count */
  double *coords      /**< input - array of vertex data */
) {
  _FC_MeshSlot *meshSlot;
  
  // NOTE!: if you change stuff here, also change fc_setSetMeshCoords()

  // check input
  meshSlot = _fc_getMeshSlot(mesh);
  if (meshSlot == NULL || dim < 1 || dim > 3 || numVertex < 1 || 
      coords == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // a little more checking -- can't overwrite coords
  if (meshSlot->coords != NULL || meshSlot->header.committed == 1) {
    fc_printfErrorMessage("Coords already exist on mesh '%s'",
                          meshSlot->header.name);
    return FC_ERROR;
  }

  // log message
  fc_printfLogMessage("Setting vertex coords on mesh '%s'", 
                      meshSlot->header.name);

  // copy meta data
  meshSlot->numVertex = numVertex;
  meshSlot->dim = dim;
  
  // big data
  meshSlot->coords = coords;
  
  return FC_SUCCESS;
}

/**
 * \ingroup  Mesh
 * \brief Set the mesh element to vertex connectivities by copy.
 *
 * \description
 *
 *    Set the element to vertex connectivities for the mesh. Meshes can only
 *    be made of one element type and the connectivity array should be
 *    numElem*(number of vertices per element) long. The order is interleave,
 *    that is you give the vertex IDs for each element (e.g. E1V1 E1V2 ...
 *    E1V6 E2V1 E2V2 ... E2V6 E3V1 ...).
 *
 *    This routine makes a copy of the conns for the library, so the
 *    user remains responsible for freeing the original array.
 *    If you would like to give the conns array to the mesh instead
 *    of just copying them, use fc_setMeshElementConnsPtr().
 *
 *    You cannot set conns if they already exist.
 *
 * \modifications  
 *    - Nancy Collins, Created.
 *    - 09/13/04 WSK Changed so you cannot reset conns
 */
FC_ReturnCode fc_setMeshElementConns (
  FC_Mesh mesh,            /**< input - mesh handle */
  FC_ElementType elemType, /**< input - element type */
  int numElem,             /**< input - count of elements */
  int *connectivities      /**< input - array of vertex IDs for each element */
) {
  //FC_ReturnCode rc;
  int num;
  _FC_MeshSlot *meshSlot;
  
  // NOTE!: if you change stuff here, also change fc_setSetMeshElementConnsPtr()

  // check input
  meshSlot = _fc_getMeshSlot(mesh);
  if (meshSlot == NULL || numElem < 1 || !fc_isElementTypeValid(elemType) ||
      elemType == FC_ET_UNKNOWN || elemType == FC_ET_MIXED ||
      connectivities == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));    
    return FC_INPUT_ERROR;
  }

  // a little more checking -- can't overwrite conns
  if (meshSlot->elemToVertConns != NULL || meshSlot->header.committed == 1) {
    fc_printfErrorMessage("Elem-to-vert conns already exist on mesh '%s'", 
                          meshSlot->header.name);
    return FC_ERROR;
  }

  // log message
  fc_printfLogMessage("Setting elements on mesh '%s'", meshSlot->header.name);

  // copy meta data
  meshSlot->numElement = numElem;
  meshSlot->elemType = elemType;
  meshSlot->topodim = fc_getElementTypeTopoDim(elemType);
  
  // copy element buffer
  num = numElem * fc_getElementTypeNumVertex(elemType);
  meshSlot->elemToVertConns = malloc(sizeof(int) * num);
  if (meshSlot->elemToVertConns == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  memcpy(meshSlot->elemToVertConns, connectivities, sizeof(int)*num);
  
  return FC_SUCCESS;
}

/**
 * \ingroup  Mesh
 * \brief Set the pointer to the mesh element connectivities.
 *
 * \description
 *
 *    Set the element to vertex connectivities for the mesh. Meshes can only
 *    be made of one element type and the connectivity array should be
 *    numElem*(number of vertices per element) long. The order is interleave,
 *    that is you give the vertex IDs for each element (e.g. E1V1 E1V2 ...
 *    E1V6 E2V1 E2V2 ... E2V6 E3V1 ...).
 *
 *    This routine passes the conns array directly to the library,
 *    so the user no longer has responsibility of the original buffer.
 *    If you would like to keep the conns array, you should use 
 *    fc_setMeshElementConns() instead which makes a copy of the conns array.
 *
 *    You cannot set conns if they already exist.
 *
 * \modifications  
 *    - 10/01/04 WSK, created.
 */
FC_ReturnCode fc_setMeshElementConnsPtr(
  FC_Mesh mesh,            /**< input - mesh handle */
  FC_ElementType elemType, /**< input - element type */
  int numElem,             /**< input - count of elements */
  int *connectivities      /**< input - array of vertex IDs for each element */
) {
  _FC_MeshSlot *meshSlot;
  
  // NOTE!: if you change stuff here, also change fc_setSetMeshElementConns()

  // check input
  meshSlot = _fc_getMeshSlot(mesh);
  if (meshSlot == NULL || numElem < 1 || !fc_isElementTypeValid(elemType) ||
      elemType == FC_ET_UNKNOWN || elemType == FC_ET_MIXED ||
      connectivities == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));    
    return FC_INPUT_ERROR;
  }

  // a little more checking -- can't overwrite conns
  if (meshSlot->elemToVertConns != NULL || meshSlot->header.committed == 1) {
    fc_printfErrorMessage("Elem-to-vert conns already exist on mesh '%s'", 
                          meshSlot->header.name);
    return FC_ERROR;
  }

  // log message
  fc_printfLogMessage("Setting elements on mesh '%s'", meshSlot->header.name);

  // copy meta data
  meshSlot->numElement = numElem;
  meshSlot->elemType = elemType;
  meshSlot->topodim = fc_getElementTypeTopoDim(elemType);
  
  // copy element buffer
  meshSlot->elemToVertConns = connectivities;
  
  return FC_SUCCESS;
}

//@}

/** \name Other ways to get new meshes. */
//-------------------------------
//@{

/**
 * \ingroup  Mesh
 * \brief  Copy a mesh
 *
 * \description
 *  
 *    Create a copy of a mesh, in the specified dataset, with the 
 *    given name. This also copies all variables on the mesh.
 *
 *    Note: if there are any sequence variables on the mesh, this
 *    routine may copy their related sequences if they do not
 *    exist in the destination dataset.
 *
 * \todo Should better define what constitutes an appropriate sequence
 *    when copying sequence variables. Right now very conservative:
 *    expects to find essentially the same sequence: same name &
 *    same data values.
 *
 * \modifications 
 *   - 2003-NOV-11 WSK Created
 *   - 10/4/04 WSK, changed internals so they call more routines and
 *       go to the structs less often.
 */
FC_ReturnCode fc_copyMesh(
  FC_Mesh src_mesh,   /**< Input - the source mesh to be copied */
  FC_Dataset dest_ds, /**< Input - the destination dataset to be copied into */
  char *newName,      /**< Input - the name of the new mesh, a NULL value
                            is a flag to use the name of the source mesh */
  FC_Mesh* new_mesh   /**< Output - handle to the new mesh */
) { 
  int i, j;
  FC_ReturnCode rc;
  int topodim, dim, numVert, numElem;
  int numSubset, numVariable, numSeqVar, *numStepPerSeq;
  FC_Subset *subsets, new_subset;
  FC_Variable *variables, new_var, **seqVars, *new_seqVar;
  FC_ElementType elemType;
  double *coords;
  int* conns;
  _FC_DsSlot* dest_dsSlot;
  _FC_MeshSlot* src_meshSlot;
  
  // default return value
  if (new_mesh)
    *new_mesh = FC_NULL_MESH;

  // check input
  dest_dsSlot = _fc_getDsSlot(dest_ds);
  src_meshSlot = _fc_getMeshSlot(src_mesh);
  if (dest_dsSlot == NULL || src_meshSlot == NULL || new_mesh == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // test for NULL name
  if (newName == NULL)
    newName = src_meshSlot->header.name;

  // log message
  fc_printfLogMessage("Copying mesh '%s' to '%s'", src_meshSlot->header.name, 
                      newName);

  // get src data
  rc = fc_getMeshInfo(src_mesh, &topodim, &dim, &numVert, &numElem,
                        &elemType);
  if (rc != FC_SUCCESS) {
    fc_printfErrorMessage("Failed to get info for mesh '%s'",
                          src_meshSlot->header.name);
    return rc;
  }
  rc = fc_getMeshCoordsPtr(src_mesh, &coords);
  if (rc != FC_SUCCESS) {
    fc_printfErrorMessage("Failed to get coords for mesh '%s'",
                          src_meshSlot->header.name);
    return rc;
  }
  rc = fc_getMeshElementConnsPtr(src_mesh, &conns);
  if (rc != FC_SUCCESS) {
    fc_printfErrorMessage("Failed to get connectivities for mesh '%s'",
                          src_meshSlot->header.name);
    return rc;
  }

  // create new mesh 
  rc = fc_createMesh(dest_ds, newName, new_mesh);
  if (rc != FC_SUCCESS) {
    fc_printfErrorMessage("Failed to get create new mesh '%s'", newName);
    return rc;
  }
  rc = fc_setMeshCoords(*new_mesh, dim, numVert, coords);
  if (rc != FC_SUCCESS) {
    fc_printfErrorMessage("Failed to set vertices on mesh '%s'", newName);
    return rc;
  }
  rc = fc_setMeshElementConns(*new_mesh, elemType, numElem, conns);
  if (rc != FC_SUCCESS) {
    fc_printfErrorMessage("Failed to set elements on mesh '%s'", newName);
    return rc;
  }

  // Note: if you change behavior here, might want to change in
  // fc_createSubsetMesh() too

  // copy the subsets
  rc = fc_getSubsets(src_mesh, &numSubset, &subsets);
  if (rc != FC_SUCCESS)
    return rc;
  for (i = 0; i < numSubset; i++) {
    rc = fc_copySubset(subsets[i], *new_mesh, NULL, &new_subset);
    if (rc != FC_SUCCESS) {
      fc_printfErrorMessage("Failed to copy member subset");
      return rc;
    }
  }
  free(subsets);

  // copy the basic variables
  rc = fc_getVariables(src_mesh, &numVariable, &variables);
  if (rc != FC_SUCCESS)
    return rc;
  for (i = 0; i < numVariable; i++) {
    // do something different for the mesh 'coords as variable'
    if (variables[i].slotID == src_meshSlot->coordVarID) {
      rc = fc_getMeshCoordsAsVariable(*new_mesh, &new_var);
      if (rc != FC_SUCCESS) {
        fc_printfErrorMessage("Failed to create copy of coords variable");
        return rc;
      }
    }
    // normal path for variables
    else {
      rc = fc_copyVariable(variables[i], *new_mesh, NULL, &new_var);
      if (rc != FC_SUCCESS) {
        fc_printfErrorMessage("Failed to copy member variable");
        return rc;
      }
    }
  }
  free(variables);

  // copy the sequence variables
  rc = fc_getSeqVariables(src_mesh, &numSeqVar, &numStepPerSeq, &seqVars);
  if (rc != FC_SUCCESS)
    return rc;
  for (i = 0; i < numSeqVar; i++) {
    FC_Sequence src_seq, dest_seq;

    // look for appropriate sequence on the destination dataset
    rc = fc_getSequenceFromSeqVariable(numStepPerSeq[i], seqVars[i],
                                       &src_seq);
    if (rc != FC_SUCCESS)
      return rc;
    dest_seq = FC_NULL_SEQUENCE;
    // if we are in the same dataset, reuse the same sequence
    if (FC_HANDLE_EQUIV(dest_ds, src_meshSlot->ds)) 
      dest_seq = src_seq;
    else { // look for an appropriate sequence in dest_ds
      int dest_numSequence, src_numStep, dest_numStep;
      FC_Sequence *dest_sequences;
      char* src_name, *dest_name;
      FC_DataType src_dataType, dest_dataType;
      void* src_coords_p, *dest_coords_p;

      // collect src sequence info      
      rc = fc_getSequenceName(src_seq, &src_name);
      if (rc != FC_SUCCESS)
        return rc;
      rc = fc_getSequenceInfo(src_seq, &src_numStep, &src_dataType); 
      if (rc != FC_SUCCESS)
        return rc;
      rc = fc_getSequenceCoordsPtr(src_seq, &src_coords_p);
      if (rc != FC_SUCCESS)
        return rc;
 
      // compare to dest sequences
      rc = fc_getSequences(dest_ds, &dest_numSequence, &dest_sequences);
      if (rc != FC_SUCCESS)
        return rc;
      for (j = 0; j < dest_numSequence; j++) {
        rc = fc_getSequenceName(dest_sequences[j], &dest_name);
        if (rc != FC_SUCCESS)
          return rc;
        rc = fc_getSequenceInfo(dest_sequences[j], &dest_numStep,
                                &dest_dataType); 
        if (rc != FC_SUCCESS)
          return rc;
        rc = fc_getSequenceCoordsPtr(dest_sequences[j], &dest_coords_p);
        if (rc != FC_SUCCESS)
          return rc;
        // compare (?checking actual coords values might be too stringent)
        if (!strcmp(dest_name, src_name) && dest_numStep == src_numStep && 
            dest_dataType == src_dataType && 
            !memcmp(dest_coords_p, src_coords_p, 
                    fc_sizeofDataType(src_dataType))) {

          free(dest_name);

          // passed all tests, this is the match
          dest_seq = dest_sequences[j];
          break;
        }
        free(dest_name);
      }
      free(dest_sequences);
      free(src_name);
    }
    // if failed to find one, just copy it
    if (FC_HANDLE_EQUIV(dest_seq, FC_NULL_SEQUENCE)) {
      rc = fc_copySequence(src_seq, dest_ds, NULL, &dest_seq);
      if (rc != FC_SUCCESS) {
        fc_printfErrorMessage("Failed to copy sequence needed for seqVar");
        return rc;
      }
    }

    // copy the sequence variables
    rc = fc_copySeqVariable(numStepPerSeq[i], seqVars[i], *new_mesh, dest_seq,
                            NULL, &new_seqVar);
    if (rc != FC_SUCCESS) {
      fc_printfErrorMessage("Failed to copy member sequence variable");
      return rc;
    }
    free(new_seqVar);
    free(seqVars[i]);
  }
  free(numStepPerSeq);
  free(seqVars);

  return FC_SUCCESS;
}

/**
 * \ingroup  Mesh
 * \brief Create a new mesh which is the subset of another mesh.
 * 
 * \description  
 *
 *    This function extracts a selected subset of a mesh as elements. 
 *    Point and element indices are renumbered to maintain the correct vertex
 *    and element relationships. (This does mean the correspondence 
 *    between vertices and elements extracted in the subset and the 
 *    original dataset are now lost.)
 *
 *    This routine automatically copies all subsets, variables and
 *    sequence variables.
 *
 *    If the association type of the subset is of lesser dimensionality than
 *    FC_AT_ELEMENT, the 'doStrict' flag determines whether to keep new
 *    entities if not all of their child entities are in the old subset. 1 = be
 *    strict and only have the entity in the new subset if all of it's children
 *    are present in the old subset. 0 = be lenient and add all entities that
 *    have at least one child in the old subset (see also
 *    fc_changeMeshEntityType).
 *
 *    NOTE: An FC_NULL_DATASET handle can be returned with no error if there
 *    is no appropriate dataset to make. (E.g. if the marked lists
 *    are all 0, or you want a 0 level mesh but not enough vertices
 *    are marked to get whole elements.)
 *
 *    NOTE: all subsets are copied even if they have no members in the
 *    new mesh (they'll just be empty subsets). The original subset is also 
 *    copied.
 *
 * \todo Create new variables that keep track of the mesh entity IDs
 *       from the originating mesh (trivial to implement).
 * 
 * \todo internals could probably be made more efficient spacewise?
 *
 * \modifications 
 *   - 2003-NOV-20  W Koegler  Replace all 6 of old subset routines
 *     with new multipurpose routine.
 *   - 08/20/04 WSK Changed to used new topo routine for figuring
 *     out which elements to keep.
 *   - 10/2/04 WSK fixing up recursive subsetting of member subsets
 *     variables & seq vars.
 */
FC_ReturnCode fc_createSubsetMesh(
  FC_Subset subset,     /**< input - subset */
  FC_Dataset dest_ds,   /**< input - dataset to put the new mesh on */
  int doStrict,         /**< input - flag for dropping incomplete entities */
  char *newName,        /**< Input - the name of the new mesh, a NULL value
                           is a flag to use the name of the source mesh */
  FC_Mesh *new_mesh  /**< output - new mesh handle */
) {
  FC_ReturnCode rc;
  int i, j, k;
  // original mesh
  int numVert, numElem, dim, numVertPerElem, numEdgePerElem, numFacePerElem;
  int subset_numVert, subset_numElem;
  // index by old ID, get new ID (-1 = doesn't exist)
  int *vertLUNewToOld, *elemLUNewToOld, *vertLUOldToNew, *elemLUOldToNew;
  int *edgeLUOldToNew = NULL, *faceLUOldToNew = NULL;
  int *edgeLUNewToOld = NULL, *faceLUNewToOld = NULL;
  int wholeLUOldToNew[1] = { 0 }, wholeLUNewToOld[1] = { 0 };
  int *LUOldToNew, *LUNewToOld;
  int* src_conns, *new_conns;
  double* src_coords, *new_coords;
  int numSubset, numVariable, numSeqVar, *numStepPerSeq;
  FC_Subset *subsets, new_subset;
  FC_Variable *variables, new_variable, **seqVars, *new_seqVar;
  _FC_DsSlot* dest_dsSlot;
  _FC_MeshSlot* src_meshSlot;
  _FC_SubSlot* subSlot;
  FC_Mesh src_mesh;
  FC_ElementType elemType;
  FC_AssociationType assoc;
  int numOrig, *origIDs;
  FC_SortedIntArray list;
  int usingEdge, usingFace;

  // default output
  if (new_mesh)
    *new_mesh = FC_NULL_MESH;
  
  // check input
  subSlot = _fc_getSubSlot(subset);
  dest_dsSlot = _fc_getDsSlot(dest_ds);
  if (subSlot == NULL || dest_dsSlot == NULL || new_mesh == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  src_meshSlot = _fc_getMeshSlot(subSlot->mesh);
  if (src_meshSlot == NULL) {
    fc_printfErrorMessage("Bad library state, subset o.k. but mesh bad, ??");
    return FC_ERROR;
  }
   
 // test for NULL name
  if (newName == NULL)
    newName = src_meshSlot->header.name;
  
  // log message
  fc_printfLogMessage("Creating new mesh '%s' from a subset '%s' of mesh "
                      "'%s'", newName, "subName", src_meshSlot->header.name);

  // get subset & mesh info
  rc = fc_getSubsetInfo(subset, &numOrig, NULL, &assoc);
  if (rc != FC_SUCCESS) {
    fc_printfErrorMessage("Failed to get subset info");
    return rc;
  }
  rc = fc_getMeshFromSubset(subset, &src_mesh);
  if (rc != FC_SUCCESS) {
    fc_printfErrorMessage("Failed to get subset's parent mesh");
    return rc;
  }
  // get some mesh info
  rc = fc_getMeshInfo(src_mesh, NULL, &dim, &numVert, &numElem, &elemType);
  if (rc != FC_SUCCESS) {
    fc_printfErrorMessage("Failed to get info for mesh '%s'",
                          src_meshSlot->header.name);
    return rc;
  }

  // trivial case: assoc = FC_AT_WHOLE_MESH
  if (assoc == FC_AT_WHOLE_MESH) {
    if (numOrig > 0) 
      return fc_copyMesh(src_mesh, dest_ds, newName, new_mesh);
    else 
      return FC_SUCCESS;
  }
  
  // get the elements that correspond with the subset
  rc = fc_getSubsetMembersAsArray(subset, &numOrig, &origIDs);
  if (rc != FC_SUCCESS) {
    fc_printfErrorMessage("Failed to get subset members");
    return rc;
  }
  rc = fc_changeMeshEntityType(src_mesh, assoc, numOrig, origIDs, 
                               FC_AT_ELEMENT, doStrict, &subset_numElem, 
                               &elemLUNewToOld);
  free(origIDs);
  if (rc != FC_SUCCESS) {
    return rc;
  }

  // more trivial cases: no elements or all elements
  if (subset_numElem == 0) {
    return FC_SUCCESS;
  }
  else if (subset_numElem == numElem) {
    free(elemLUNewToOld);
    return fc_copyMesh(src_mesh, dest_ds, newName, new_mesh);
  }

  // more mesh info
  numVertPerElem = fc_getElementTypeNumVertex(elemType);
  numEdgePerElem = fc_getElementTypeNumEdge(elemType);
  numFacePerElem = fc_getElementTypeNumFace(elemType);
  rc = fc_getMeshElementConnsPtr(src_mesh, &src_conns);
  if (rc != FC_SUCCESS)
    return rc;

  // make element lookup tables
  elemLUOldToNew = malloc(numElem*sizeof(int));
  if (elemLUOldToNew == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  for (i = 0; i < numElem; i++)
    elemLUOldToNew[i] = -1;
  for (i = 0; i < subset_numElem; i++) 
    elemLUOldToNew[elemLUNewToOld[i]] = i;

  // Figure out which vertices are going to be kept & make lookup tables
  fc_initSortedIntArray(&list);
  for (i = 0; i < subset_numElem; i++) {
    int ret = fc_addIntArrayToSortedIntArray(&list, numVertPerElem,
					     &src_conns[elemLUNewToOld[i]*numVertPerElem], 0);
    if (ret < 0) {
      fc_printfErrorMessage("Failed to add ids to sia");
      return ret;    
    }
  }
  rc = fc_convertSortedIntArrayToIntArray(&list, &subset_numVert, &vertLUNewToOld);
  if (rc != FC_SUCCESS)
    return rc;
  vertLUOldToNew = malloc(numVert*sizeof(int));
  if (vertLUOldToNew == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  for (i = 0; i < numVert; i++)
    vertLUOldToNew[i] = -1;
  for (i = 0; i < subset_numVert; i++)
    vertLUOldToNew[vertLUNewToOld[i]] = i;

  // make element->vertex conns
  new_conns = malloc(subset_numElem * numVertPerElem * sizeof(int));
  if (new_conns == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  for (i = 0; i < subset_numElem; i++)
    for (j = 0; j < numVertPerElem; j++) 
      new_conns[i*numVertPerElem + j] = 
        vertLUOldToNew[src_conns[elemLUNewToOld[i]*numVertPerElem+j]];

  // determine new coords
  new_coords = malloc(subset_numVert * dim * sizeof(double));
  if (new_coords == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  rc = fc_getMeshCoordsPtr(src_mesh, &src_coords);
  if (rc != FC_SUCCESS)
    return rc;
  for (i = 0; i < subset_numVert; i++)
    memcpy(&(new_coords[i*dim]), &(src_coords[vertLUNewToOld[i]*dim]), 
           dim*sizeof(double));
  
  // create new mesh
  rc = fc_createMesh(dest_ds, newName, new_mesh);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_setMeshCoordsPtr(*new_mesh, dim, subset_numVert, new_coords);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_setMeshElementConnsPtr(*new_mesh, elemType, subset_numElem, 
                                 new_conns);
  if (rc != FC_SUCCESS)
    return rc;

  // Get subsets & variables for use now & later
  rc = fc_getSubsets(src_mesh, &numSubset, &subsets);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_getVariables(src_mesh, &numVariable, &variables);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_getSeqVariables(src_mesh, &numSeqVar, &numStepPerSeq, &seqVars);
  if (rc != FC_SUCCESS)
    return rc;

  // Build edge & face lu tables only if there are edge or face 
  // subsets or variables
  usingEdge = 0;
  usingFace = 0;
  for (i = 0; i < numSubset; i++) {
    fc_getSubsetInfo(subsets[i], NULL, NULL, &assoc);
    if (assoc == FC_AT_EDGE)
      usingEdge = 1;
    else if (assoc == FC_AT_FACE)
      usingFace = 1;
  }
  for (i = 0; i < numVariable; i++) {
    fc_getVariableInfo(variables[i], NULL, NULL, &assoc, NULL, NULL);
    if (assoc == FC_AT_EDGE)
      usingEdge = 1;
    else if (assoc == FC_AT_FACE)
      usingFace = 1;
  }    
  for (i = 0; i < numSeqVar; i++) {
    fc_getVariableInfo(seqVars[i][0], NULL, NULL, &assoc, NULL, NULL);
    if (assoc == FC_AT_EDGE)
      usingEdge = 1;
    else if (assoc == FC_AT_FACE)
      usingFace = 1;
  }    
  if (usingEdge) {
    int numEdge, new_numElement, new_numEdge;
    int *src_elemToEdgeConns, *new_elemToEdgeConns;

    // initialize old to new lookup tables
    fc_getMeshNumEntity(src_mesh, FC_AT_EDGE, &numEdge);
    if (numEdge > 0) {
      edgeLUOldToNew = malloc(numEdge*sizeof(int));
      if (edgeLUOldToNew == NULL) {
        fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
        return FC_MEMORY_ERROR;
      }
      for (i = 0; i < numEdge; i++)
        edgeLUOldToNew[i] = -1;
    }

    // create new to old lu (& force building of edges & faces on new mesh)
    fc_getMeshNumEntity(*new_mesh, FC_AT_EDGE, &new_numEdge);
    if (new_numEdge > 0) {
      edgeLUNewToOld = malloc(new_numEdge*sizeof(int));
       if (edgeLUOldToNew == NULL) {
        fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
        return FC_MEMORY_ERROR;
      }
    }

    // create lookups
    fc_getMeshNumElement(*new_mesh, &new_numElement);
    rc = fc_getMeshElementToEdgeConnsPtr(src_mesh, &src_elemToEdgeConns);
    if (rc != FC_SUCCESS)
      return rc;
    rc = fc_getMeshElementToEdgeConnsPtr(*new_mesh, &new_elemToEdgeConns);
    if (rc != FC_SUCCESS)
      return rc;
    for (i = 0; i < new_numElement; i++) {
      int* new_edgeIDs = &new_elemToEdgeConns[i*numEdgePerElem];
      int* old_edgeIDs = &src_elemToEdgeConns[elemLUNewToOld[i]*numEdgePerElem];
      for (j = 0; j < numEdgePerElem; j++) {
        edgeLUNewToOld[new_edgeIDs[j]] = old_edgeIDs[j];
        edgeLUOldToNew[old_edgeIDs[j]] = new_edgeIDs[j];
      }
    }
  }
  if (usingFace) {
    int numFace, new_numElement, new_numFace;
    int *src_elemToFaceConns, *new_elemToFaceConns;

    // initialize old to new lookup tables
    fc_getMeshNumEntity(src_mesh, FC_AT_FACE, &numFace);
    if (numFace > 0) {
      faceLUOldToNew = malloc(numFace*sizeof(int));
      if (faceLUOldToNew == NULL) {
        fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
        return FC_MEMORY_ERROR;
      }
      for (i = 0; i < numFace; i++)
        faceLUOldToNew[i] = -1;
    }

    // create new to old lu (& force building of edges & faces on new mesh)
    fc_getMeshNumEntity(*new_mesh, FC_AT_FACE, &new_numFace);
    if (new_numFace > 0) {
      faceLUNewToOld = malloc(new_numFace*sizeof(int));
       if (faceLUOldToNew == NULL) {
        fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
        return FC_MEMORY_ERROR;
      }
    }     

    // create lookups
    fc_getMeshNumElement(*new_mesh, &new_numElement);
    rc = fc_getMeshElementToFaceConnsPtr(src_mesh, &src_elemToFaceConns);
    if (rc != FC_SUCCESS)
      return rc;
    rc = fc_getMeshElementToFaceConnsPtr(*new_mesh, &new_elemToFaceConns);
    if (rc != FC_SUCCESS)
      return rc;
    for (i = 0; i < new_numElement; i++) {
      int* new_faceIDs = &new_elemToFaceConns[i*numFacePerElem];
      int* old_faceIDs = &src_elemToFaceConns[elemLUNewToOld[i]*numFacePerElem];
      for (j = 0; j < numFacePerElem; j++) {
        faceLUNewToOld[new_faceIDs[j]] = old_faceIDs[j];
        faceLUOldToNew[old_faceIDs[j]] = new_faceIDs[j];
      }
    }
  }

  // copy stuff -- a lot like in fc_copyMesh

  // copy subsets
  for (i = 0; i < numSubset; i++) {
    int numMember, *memberIDs;
    char* subset_name;
    FC_AssociationType subset_assoc;

    // create the new subset
    rc = fc_getSubsetAssociationType(subsets[i], &subset_assoc);
    if (rc != FC_SUCCESS)
      return rc;
    rc = fc_getSubsetName(subsets[i], &subset_name);
    if (rc != FC_SUCCESS)
      return rc;
    rc = fc_createSubset(*new_mesh, subset_name, subset_assoc, &new_subset);
    free(subset_name);
    if (rc != FC_SUCCESS)
      return rc;
    
    // Add the members that are in the new subset (i.e. have a new ID)
    switch (subset_assoc) {
    case FC_AT_VERTEX:   LUOldToNew = vertLUOldToNew;   break;
    case FC_AT_EDGE:     LUOldToNew = edgeLUOldToNew;   break;
    case FC_AT_FACE:     LUOldToNew = faceLUOldToNew;   break;
    case FC_AT_ELEMENT:  LUOldToNew = elemLUOldToNew;   break;
    case FC_AT_WHOLE_MESH:    LUOldToNew = wholeLUOldToNew;  break;
    default:
      fc_printfErrorMessage("Should never reach this point!");
      return rc;
    }
    rc = fc_getSubsetMembersAsArray(subsets[i], &numMember, &memberIDs);
    if (rc != FC_SUCCESS)
      return rc;
    for (j = 0; j < numMember; j++) {
      if (LUOldToNew[memberIDs[j]] > -1) {
        rc = fc_addMemberToSubset(new_subset, LUOldToNew[memberIDs[j]]);
        if (rc != FC_SUCCESS)
          return rc;
      }
    }
    free(memberIDs);
  }
  free(subsets);

  // copy the basic variables
  for (i = 0; i < numVariable; i++) {
    // do something different for the mesh 'coords as variable'
    if (variables[i].slotID == src_meshSlot->coordVarID) {
      rc = fc_getMeshCoordsAsVariable(*new_mesh, &new_variable);
      if (rc != FC_SUCCESS) {
        fc_printfErrorMessage("Failed to create copy of coords variable");
        return rc;
      }
    }
    else { // normal path for variables
      int numDataPoint, numComp, new_numDataPoint, size;
      char* var_name;
      FC_AssociationType var_assoc;
      FC_MathType mathtype;
      FC_DataType datatype;
      void* src_data, *new_data;
      
      // make the variable
      rc = fc_getVariableName(variables[i], &var_name);
      if (rc != FC_SUCCESS)
        return rc;
      rc = fc_createVariable(*new_mesh, var_name, &new_variable);
      free(var_name);
      if (rc != FC_SUCCESS)
        return rc;

      // get the original data
      rc = fc_getVariableInfo(variables[i], &numDataPoint, &numComp, 
                              &var_assoc, &mathtype, &datatype);
      if (rc != FC_SUCCESS)
        return rc;
      fc_getVariableDataPtr(variables[i], &src_data);
      if (rc != FC_SUCCESS)
        return rc;

      // create subset of data
      rc = fc_getMeshNumEntity(*new_mesh, var_assoc, &new_numDataPoint);
      if (rc != FC_SUCCESS)
        return rc;
      size = fc_sizeofDataType(datatype)*numComp;
      new_data = malloc(new_numDataPoint*size);
      if (new_data == NULL) {
        fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
        return FC_MEMORY_ERROR;
      }
      switch (var_assoc) {
      case FC_AT_VERTEX:   LUNewToOld = vertLUNewToOld;   break;
      case FC_AT_EDGE:     LUNewToOld = edgeLUNewToOld;   break;
      case FC_AT_FACE:     LUNewToOld = faceLUNewToOld;   break;
      case FC_AT_ELEMENT:  LUNewToOld = elemLUNewToOld;   break;
      case FC_AT_WHOLE_MESH:    LUNewToOld = wholeLUNewToOld;  break;
      default:
        fc_printfErrorMessage("Should never reach this point!");
        return rc;
      }
      for (j = 0; j < new_numDataPoint; j++)
        memcpy((char*)new_data + j*size, 
               (char*)src_data + LUNewToOld[j]*size, size);
      
      // set the data
      rc = fc_setVariableDataPtr(new_variable, new_numDataPoint, numComp, 
                                 var_assoc, mathtype, datatype, new_data);
      if (rc != FC_SUCCESS)
        return rc;
    }
  }
  free(variables);

  // copy the sequence vars
  for (i = 0; i < numSeqVar; i++) {
    int numStep;
    FC_Sequence src_seq, dest_seq;
    char* var_name;

    // look for appropriate sequence on the destination dataset
    rc = fc_getSequenceFromSeqVariable(numStepPerSeq[i], seqVars[i],
                                       &src_seq);
    if (rc != FC_SUCCESS)
      return rc;
    dest_seq = FC_NULL_SEQUENCE;
    // if we are in the same dataset, reuse the same sequence
    if (FC_HANDLE_EQUIV(dest_ds, src_meshSlot->ds)) 
      dest_seq = src_seq;
    else { // look for an appropriate sequence in dest_ds
      int dest_numSequence, src_numStep, dest_numStep;
      FC_Sequence *dest_sequences;
      char* src_name, *dest_name;
      FC_DataType src_dataType, dest_dataType;
      void* src_coords_p, *dest_coords_p;

      // collect src sequence info      
      rc = fc_getSequenceName(src_seq, &src_name);
      if (rc != FC_SUCCESS)
        return rc;
      rc = fc_getSequenceInfo(src_seq, &src_numStep, &src_dataType); 
      if (rc != FC_SUCCESS)
        return rc;
      rc = fc_getSequenceCoordsPtr(src_seq, &src_coords_p);
      if (rc != FC_SUCCESS)
        return rc;
 
      // compare to dest sequences
      rc = fc_getSequences(dest_ds, &dest_numSequence, &dest_sequences);
      if (rc != FC_SUCCESS)
        return rc;
      for (j = 0; j < dest_numSequence; j++) {
        rc = fc_getSequenceName(dest_sequences[j], &dest_name);
        if (rc != FC_SUCCESS)
          return rc;
        rc = fc_getSequenceInfo(dest_sequences[j], &dest_numStep,
                                &dest_dataType); 
        if (rc != FC_SUCCESS)
          return rc;
        rc = fc_getSequenceCoordsPtr(dest_sequences[j], &dest_coords_p);
        if (rc != FC_SUCCESS)
          return rc;
        // compare (?checking actual coords values might be too stringent)
        if (!strcmp(dest_name, src_name) && dest_numStep == src_numStep && 
            dest_dataType == src_dataType && 
            !memcmp(dest_coords_p, src_coords_p, 
                    fc_sizeofDataType(src_dataType))) {
          free(dest_name);
          // passed all tests, this is the match
          dest_seq = dest_sequences[j];
          break;
        }
        free(dest_name);
      }
      free(dest_sequences);
      free(src_name);
    }
    // if failed to find one, just copy it
    if (FC_HANDLE_EQUIV(dest_seq, FC_NULL_SEQUENCE)) {
      rc = fc_copySequence(src_seq, dest_ds, NULL, &dest_seq);
      if (rc != FC_SUCCESS) {
        fc_printfErrorMessage("Failed to copy sequence needed for seqVar");
        return rc;
      }
    }

    // make the seq variable
    rc = fc_getVariableName(seqVars[i][0], &var_name);
    if (rc != FC_SUCCESS)
      return rc;
    rc = fc_createSeqVariable(*new_mesh, dest_seq, var_name, &numStep,
                              &new_seqVar);
    free(var_name);
    if (rc != FC_SUCCESS)
      return rc;

    // set data on steps
    for (j = 0; j < numStepPerSeq[i]; j++) {
      int numDataPoint, numComp, new_numDataPoint, size;
      FC_AssociationType var_assoc;
      FC_MathType mathtype;
      FC_DataType datatype;
      void* src_data, *new_data;

      // get the original data
      rc = fc_getVariableInfo(seqVars[i][j], &numDataPoint, &numComp, 
                              &var_assoc, &mathtype, &datatype);
      if (rc != FC_SUCCESS)
        return rc;
      fc_getVariableDataPtr(seqVars[i][j], &src_data);
      if (rc != FC_SUCCESS)
        return rc;
    
      // create subset of data
      rc = fc_getMeshNumEntity(*new_mesh, var_assoc, &new_numDataPoint);
      if (rc != FC_SUCCESS)
        return rc;
      size = fc_sizeofDataType(datatype)*numComp;
      new_data = malloc(new_numDataPoint*size);
      if (new_data == NULL) {
        fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
        return FC_MEMORY_ERROR;
      }
      switch (var_assoc) {
      case FC_AT_VERTEX:   LUNewToOld = vertLUNewToOld;   break;
      case FC_AT_EDGE:     LUNewToOld = edgeLUNewToOld;   break;
      case FC_AT_FACE:     LUNewToOld = faceLUNewToOld;   break;
      case FC_AT_ELEMENT:  LUNewToOld = elemLUNewToOld;   break;
      case FC_AT_WHOLE_MESH:    LUNewToOld = wholeLUNewToOld;  break;
      default:
        fc_printfErrorMessage("Should never reach this point!");
        return rc;
      }
      for (k = 0; k < new_numDataPoint; k++)
        memcpy((char*)new_data + k*size, 
               (char*)src_data + LUNewToOld[k]*size, size);

      // set the data
      rc = fc_setVariableDataPtr(new_seqVar[j], new_numDataPoint, numComp, 
                                 var_assoc, mathtype, datatype, new_data);
      if (rc != FC_SUCCESS)
        return rc;
    }
    free(new_seqVar);
    free(seqVars[i]);
  }
  free(numStepPerSeq);
  free(seqVars);

  // cleanup
  free(vertLUOldToNew);
  free(edgeLUOldToNew);
  free(faceLUOldToNew);
  free(elemLUOldToNew);
  free(vertLUNewToOld);
  free(edgeLUNewToOld);
  free(faceLUNewToOld);
  free(elemLUNewToOld);

  return FC_SUCCESS;
}

/**
 * \ingroup  Mesh
 * \brief Create a simple hex mesh.
 * 
 * \description  
 *
 *    This function creates a hex mesh which is a regular parallelpiped (if
 *    all the sides are the same length, a cube). The number of elements in
 *    each dimension are given by N, M, & K. The coordinates are specified by
 *    giving the coordinates of the lower and upper corners of the mesh.
 *
 * \modifications 
 *   - 11/09/2005  WSD  Created.
 */
FC_ReturnCode fc_createSimpleHexMesh(
  FC_Dataset dataset, /**< Input - dataset handle to associate mesh with */
  char *newName,     /**< Input - the name of the new mesh */
  int N,             /**< Input - number of elements in the x direction */
  int M,             /**< Input - number of elements in the y direction */
  int P,             /**< Input - number of elements in the z direction */  
  double* lower_coords, /**< Input - the coordinates of lower corner. */
  double* upper_coords, /**< Input - the coordinates of upper corner. */
  FC_Mesh *new_mesh  /**< output - new mesh handle */
) {
  FC_ReturnCode rc;
  int i, j, k;
  int numElemPerDim[3] = { N, M, P }, eindex, vindex;
  int numDim = 3, numVertPerElem = 8;
  int numVert, numElem;
  double side_lengths[3];
  int *conns;
  double* coords;

  // default output
  if (new_mesh)
    *new_mesh = FC_NULL_MESH;

  // check input
  if (!fc_isDatasetValid(dataset) || !newName || N < 1 || M < 1 || P < 1 || 
      !lower_coords || !upper_coords || !new_mesh) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  for (i = 0; i < numDim; i++) {
    if (lower_coords[i] >= upper_coords[i]) {
      fc_printfErrorMessage("Upper coords for dim%d must be greater than lower "
                            "coords", i);
      return FC_INPUT_ERROR;
    }
  }

  // log message
  fc_printfLogMessage("Creating simple hex mesh '%s'", newName);

  // Set up
  numVert = (M+1)*(N+1)*(P+1);
  numElem = M*N*P;
  conns = (int*)malloc(numElem*numVertPerElem*sizeof(int));
  coords = (double*)malloc(numVert*numDim*sizeof(double));
  if (!conns || !coords) {
    free(conns);
    free(coords);
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  for (i = 0; i < numDim; i++)
    side_lengths[i] = (upper_coords[i] - lower_coords[i])/numElemPerDim[i];

  // Make coords
  for (k = 0; k < P+1; k++) {
    for (j = 0; j < M+1; j++) {
      for (i = 0; i < N+1; i++) {
        vindex = k*(N+1)*(M+1) + j*(N+1) + i;
        coords[vindex*numDim] =   i*side_lengths[0] + lower_coords[0];
        coords[vindex*numDim+1] = j*side_lengths[1] + lower_coords[1];
        coords[vindex*numDim+2] = k*side_lengths[2] + lower_coords[2];
      }
    }
  }

  // Make conns
  for (k = 0; k < P; k++) {
    for (j = 0; j < M; j++) {
      for (i = 0; i < N; i++) {
        eindex = k*N*M + j*N + i;
        vindex = k*(N+1)*(M+1) + j*(N+1) + i; // index of lowest vert
        conns[eindex*numVertPerElem]   = vindex;
        conns[eindex*numVertPerElem+1] = vindex + 1;
        vindex = k*(N+1)*(M+1) + (j+1)*(N+1) + i;
        conns[eindex*numVertPerElem+2] = vindex + 1;
        conns[eindex*numVertPerElem+3] = vindex;
        vindex = (k+1)*(N+1)*(M+1) + j*(N+1) + i;
        conns[eindex*numVertPerElem+4] = vindex;
        conns[eindex*numVertPerElem+5] = vindex + 1;
        vindex = (k+1)*(N+1)*(M+1) + (j+1)*(N+1) + i;
        conns[eindex*numVertPerElem+6] = vindex + 1;
        conns[eindex*numVertPerElem+7] = vindex;
      }
    }
  }

  // Make the mesh
  rc = fc_createMesh(dataset, newName, new_mesh);
  if (rc != FC_SUCCESS) {
    free(coords);
    free(conns);
    return rc;
  }
  rc = fc_setMeshCoordsPtr(*new_mesh, numDim, numVert, coords);
  if (rc != FC_SUCCESS) {
    free(conns);
    fc_deleteMesh(*new_mesh);
    return rc;
  }
  rc = fc_setMeshElementConnsPtr(*new_mesh, FC_ET_HEX, numElem, conns);
  if (rc != FC_SUCCESS) {
    fc_deleteMesh(*new_mesh);
    return rc;
  }

  // All done!
  return FC_SUCCESS;
}

//@}

/** \name Get existing meshes. */
//-------------------------------
//@{

/**
 * \ingroup  Mesh
 * \brief  Get the meshes in a dataset.
 *
 * \description
 *  
 *   This call returns an array of meshes in the dataset. The caller is
 *   responsible for freeing the array of sequence handles.
 *
 * \modifications  
 *   - Nancy Collins  Created.
 *   - July 23, 2002  W Koegler  Changed so it no longer makes safcover calls
 *   - 11/21/03 RM, added default return nmeshes = -1, and meshnames = NULL 
 *                    if anything is wrong
 *   - 9/8/04 WSK, changed to return mesh handles instead of mesh
 *               names.
 */
FC_ReturnCode fc_getMeshes(
 FC_Dataset dataset,  /**< input - dataset */
 int *numMesh,        /**< output - number of meshes returned */
 FC_Mesh **meshes     /**< output - array of meshes */
) {
  int i;
  _FC_DsSlot* dsSlot;
  _FC_MeshSlot* meshSlot;
  
  // default return values
  if (numMesh != NULL) 
    *numMesh = -1;
  if (meshes != NULL)
    *meshes = NULL;  
  
  // check input
  dsSlot = _fc_getDsSlot(dataset);
  if (dsSlot == NULL || numMesh == NULL || meshes == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Getting meshes from dataset '%s'",
                      dsSlot->header.name);
  
  // get num meshes
  *numMesh = dsSlot->numMesh;
  
  // get meshes
  if (dsSlot->numMesh > 0) {
    *meshes = malloc(sizeof(FC_Mesh)*dsSlot->numMesh);
    if (*meshes == NULL) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      return FC_MEMORY_ERROR;
    }
    for (i = 0; i < dsSlot->numMesh; i++) {
      meshSlot = _fc_getMeshSlotFromID(dsSlot->meshIDs[i]);
      _FC_GET_HANDLE((*meshes)[i], meshSlot);
    }
  }
  
  return FC_SUCCESS;
}

/**
 * \ingroup  Mesh
 * \brief  Get the the number of meshes in a dataset.
 *
 * \modifications  
 *   - 9/28/04 WSK, created.
 */
FC_ReturnCode fc_getNumMesh(
 FC_Dataset dataset,  /**< input - dataset */
 int *numMesh         /**< output - number of meshes returned */
) {
  _FC_DsSlot* dsSlot;
  
  // default return values
  if (numMesh != NULL) 
    *numMesh = -1;
  
  // check input
  dsSlot = _fc_getDsSlot(dataset);
  if (dsSlot == NULL || numMesh == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Getting number of meshes from dataset '%s'",
                      dsSlot->header.name);
  
  // get num meshes
  *numMesh = dsSlot->numMesh;
  
  return FC_SUCCESS;
}

/**
 * \ingroup  Mesh
 * \brief  Get all meshes with a given name from the dataset.
 *
 * \description
 *  
 *    Returns all meshes on the dataset with the specified name.
 *
 * \modifications  
 *   - Nancy Collins  Created.
 *   - 2003-MAY-30  W Koegler  Moved getting mesh meta data from 
 *       _fc_setupDsSlot() to here. Moved getting variable meta data 
 *       to fc_openVariable().
 *   - 11/25/03 RM, added default return mesh = \ref FC_NULL_MESH
 *   - 11/26/03 RM, added checking for invalid mesh name, issue error if 
 *        fc_getLibraryVerbosity()  
 *   - 12/03/03 RM  added checking for parameter mesh == NULL 
 *   - 12/12/03 RM, eliminated _fc_lookup_meshSlot() and _fc_namematch()
 *       those input checking changed and added into this function  
 *   - 9/8/03 WSK, changed name and moved to dataset (to be next to
 *       fc_getMeshes) 
 *   - 9/21/06 ACG, now returns all meshes with matching names
 */
FC_ReturnCode fc_getMeshByName(
  FC_Dataset dataset, /**< input - dataset handle */
  char *meshname,     /**< input - name of mesh to get */
  int *numMeshes,       /**< output - number of matching meshes */
  FC_Mesh **meshes       /**< output - array of mesh handles */
) {
  int i;
  _FC_DsSlot* dsSlot;
  _FC_MeshSlot *meshSlot;
  FC_Mesh mesh;

  FC_Mesh *lmatch = NULL;
  int numlmatch, maxnumlmatch;
  
  // default values returned if anything goes wrong
  if (numMeshes)
    *numMeshes = -1;
  if (meshes)
    *meshes = NULL;
  
  // check input
  dsSlot = _fc_getDsSlot(dataset);
  if (dsSlot == NULL || meshname == NULL || !meshes || !numMeshes) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Getting mesh '%s'", meshname);

  numlmatch = 0;
  maxnumlmatch = 0;
  
  // look for existing slots with the matching name
  for (i = 0; i < dsSlot->numMesh; i++) {
    meshSlot = _fc_getMeshSlotFromID(dsSlot->meshIDs[i]);
    if (!strcmp(meshSlot->header.name, meshname)) {
      _FC_GET_HANDLE(mesh, meshSlot);

      if (lmatch == NULL){
	lmatch = (FC_Mesh*)malloc(sizeof(FC_Mesh));
	if (lmatch == NULL){
	  fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
	  return FC_MEMORY_ERROR;
	}
	numlmatch = 1;
	maxnumlmatch = 1;
	lmatch[0] = mesh;
      }else{
	if(numlmatch == maxnumlmatch){
	  FC_Mesh *temp;
	  temp = (FC_Mesh*)realloc(lmatch,(2*numlmatch)*sizeof(FC_Mesh));
	  if (temp == NULL){
	    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
	    return FC_MEMORY_ERROR;
	  }
	  maxnumlmatch*=2;
	  lmatch = temp;
	}
	lmatch[numlmatch] = mesh;
	numlmatch++;
      }
    }
  }

  *numMeshes = numlmatch;
  *meshes = lmatch;

  return FC_SUCCESS;
}

//@}

/** \name Change the name of a mesh. */
//-------------------------------
//@{

/**
 * \ingroup Mesh
 * \brief  Change the name of a mesh.
 *
 * \modifications
 *   - 8/24/2005 WSD. Created.
 */
FC_ReturnCode fc_changeMeshName(
  FC_Mesh mesh, /**< Input - the mesh. */
  char* newName       /**< Input - new name for the mesh. */
) {
  _FC_MeshSlot* meshSlot;

  // check input
  meshSlot = _fc_getMeshSlot(mesh);
  if (meshSlot == NULL || newName == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // log message
  fc_printfLogMessage("Changing name of mesh from '%s' to '%s'", 
                      meshSlot->header.name, newName);
  
  // Do it
  return _fc_setSlotHeaderName(&meshSlot->header, newName);
}

//@}

/** \name Release/delete meshes. */
//-------------------------------
//@{

/**
 * \ingroup Mesh
 * \brief Attempt to minimize mesh's memory usage.
 *
 * \description
 *
 *    Call to try to release un-needed memory in a mesh. If the mesh has been
 *    saved to disk, the large data is released. If the mesh has not be saved
 *    to disk, all large data is saved.
 *
 *    If internal mesh structures like elements and vertices exists, they
 *    may be freed as long as there are no subsets or variables associated with
 *    edges or faces.
 *
 *    This function will recursively call release on it's child subsets and
 *    variables
 *
 * \modifications
 *    - 9/27/04 WSK Created.
 */
FC_ReturnCode fc_releaseMesh(
 FC_Mesh mesh       /**< input - mesh handle */
) {
  int i;
  int *numSteps;
  FC_Subset *subsets;
  FC_Variable *variables, **seqVariables;
  _FC_MeshSlot* meshSlot;
  FC_AssociationType assoc;
  int usingEdge, usingFace, numSubset, numVar, numSeqVar;
  
  // special case: releasing a null handle is not an error
  if ( FC_HANDLE_EQUIV(mesh, FC_NULL_MESH)) {
    fc_printfWarningMessage("Releasing FC_NULL_MESH");
    return FC_SUCCESS;
  }

  // check input
  meshSlot = _fc_getMeshSlot(mesh);
  if (meshSlot == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Cleaning up mesh '%s'", meshSlot->header.name);

  // call release on children subsets & variables
  // also note whether any children use edges or faces
  usingEdge = 0;
  usingFace = 0;
  fc_getSubsets(mesh, &numSubset, &subsets);
  for (i = 0; i < numSubset; i++) {
    fc_releaseSubset(subsets[i]);
    fc_getSubsetInfo(subsets[i], NULL, NULL, &assoc);
    if (assoc == FC_AT_EDGE)
      usingEdge = 1;
    else if (assoc == FC_AT_FACE) 
      usingFace = 1;
  }
  free(subsets);
  fc_getVariables(mesh, &numVar, &variables);
  for (i = 0; i < numVar; i++) {
    fc_releaseVariable(variables[i]);
    fc_getVariableInfo(variables[i], NULL, NULL, &assoc, NULL, NULL); 
    if (assoc == FC_AT_EDGE)
      usingEdge = 1;
    else if (assoc == FC_AT_FACE) 
      usingFace = 1;
  }
  free(variables);
  fc_getSeqVariables(mesh, &numSeqVar, &numSteps, &seqVariables);
  for (i = 0; i < numSeqVar; i++) {
    fc_releaseSeqVariable(numSteps[i], seqVariables[i]);
    fc_getVariableInfo(seqVariables[i][0], NULL, NULL, &assoc, NULL, NULL); 
    if (assoc == FC_AT_EDGE)
      usingEdge = 1;
    else if (assoc == FC_AT_FACE) 
      usingFace = 1;
  }
  free(numSteps);
  for (i = 0; i < numSeqVar; i++)
    free(seqVariables[i]);
  free(seqVariables);
             
  // remove coords & conns if this slot is committed (they can be reloaded)
  if (meshSlot->header.committed == 1) {
    free(meshSlot->coords);
    meshSlot->coords = NULL;
    free(meshSlot->elemToVertConns);
    meshSlot->elemToVertConns = NULL;
  }
 
  //-- remove parent lists
  free(meshSlot->numEdgePerVert);
  meshSlot->numEdgePerVert = NULL;
  if (meshSlot->edgeParentsPerVert)
      for (i = 0; i < meshSlot->numVertex; i++)
        free(meshSlot->edgeParentsPerVert[i]);
  free(meshSlot->edgeParentsPerVert);
  meshSlot->edgeParentsPerVert = NULL;
  //--
  free(meshSlot->numFacePerVert);
  meshSlot->numFacePerVert = NULL;
  if (meshSlot->faceParentsPerVert)
    for (i = 0; i < meshSlot->numVertex; i++)
      free(meshSlot->faceParentsPerVert[i]);
  free(meshSlot->faceParentsPerVert);
  meshSlot->faceParentsPerVert = NULL;
  //--
  free(meshSlot->numElemPerVert);
  meshSlot->numElemPerVert = NULL;
  if (meshSlot->elemParentsPerVert)
    for (i = 0; i < meshSlot->numVertex; i++)
      free(meshSlot->elemParentsPerVert[i]);
  free(meshSlot->elemParentsPerVert);
  meshSlot->elemParentsPerVert = NULL;
  //--
  free(meshSlot->numElemPerEdge);
  meshSlot->numElemPerEdge = NULL;
  if (meshSlot->elemParentsPerEdge)
    for (i = 0; i < meshSlot->numEdge; i++)
      free(meshSlot->elemParentsPerEdge[i]);
  free(meshSlot->elemParentsPerEdge);
  meshSlot->elemParentsPerEdge = NULL;
  //--
  free(meshSlot->numElemPerFace);
  meshSlot->numElemPerFace = NULL;
  if (meshSlot->elemParentsPerFace)
    for (i = 0; i < meshSlot->numFace; i++)
      free(meshSlot->elemParentsPerFace[i]);
  free(meshSlot->elemParentsPerFace);
  meshSlot->elemParentsPerFace = NULL;

  //-- remove neighbor lists
  free(meshSlot->numVertNeighsViaEdge);
  meshSlot->numVertNeighsViaEdge = NULL;
  if (meshSlot->vertNeighsViaEdge)
    for (i = 0; i < meshSlot->numVertex; i++)
      free(meshSlot->vertNeighsViaEdge[i]);
  free(meshSlot->vertNeighsViaEdge);
  meshSlot->vertNeighsViaEdge = NULL;
  //--
  free(meshSlot->numElemNeighsViaVert);
  meshSlot->numElemNeighsViaVert = NULL;
  if (meshSlot->elemNeighsViaVert)
    for (i = 0; i < meshSlot->numElement; i++)
      free(meshSlot->elemNeighsViaVert[i]);
  free(meshSlot->elemNeighsViaVert);
  meshSlot->elemNeighsViaVert = NULL;
  //--
  free(meshSlot->numElemNeighsViaEdge);
  meshSlot->numElemNeighsViaEdge = NULL;
  if (meshSlot->elemNeighsViaEdge)
    for (i = 0; i < meshSlot->numElement; i++)
      free(meshSlot->elemNeighsViaEdge[i]);
  free(meshSlot->elemNeighsViaEdge);
  meshSlot->elemNeighsViaEdge = NULL;
  //--
  free(meshSlot->numElemNeighsViaFace);
  meshSlot->numElemNeighsViaFace = NULL;
  if (meshSlot->elemNeighsViaFace)
    for (i = 0; i < meshSlot->numElement; i++)
      free(meshSlot->elemNeighsViaFace[i]);
  free(meshSlot->elemNeighsViaFace);
  meshSlot->elemNeighsViaFace = NULL;

  // if no subsets or variables use edges or faces, conns & associated data
  // do this last because we might want numEdge or numFace
  if (!usingEdge) {
    free(meshSlot->edgeToVertConns);
    meshSlot->edgeToVertConns = NULL;
    free(meshSlot->elemToEdgeConns);
    meshSlot->elemToEdgeConns = NULL;
    //-- set edges back to "unknown"
    meshSlot->numEdge = -1;
  }
  if (!usingFace) {
    free(meshSlot->faceTypes);
    meshSlot->faceTypes = NULL;
    free(meshSlot->numVertPerFace);
    meshSlot->numVertPerFace = NULL;
    free(meshSlot->faceToVertConns);
    meshSlot->faceToVertConns = NULL;
    free(meshSlot->elemToFaceConns);
    meshSlot->elemToFaceConns = NULL;
    free(meshSlot->elemFaceOrients);
    meshSlot->elemFaceOrients = NULL;
    //-- set faces back to "unknown"
    meshSlot->numFace = -1;
  }

  return FC_SUCCESS;
}

/**
 * \ingroup  Mesh
 * \brief  Delete a mesh.
 *
 * \description
 *  
 *    Call when this mesh is no longer needed. This automatically deletes
 *    all of the mesh's subsets and variables.
 *
 * \modifications 
 *    - Nancy Collins, Created.
 *    - 12/01/03 WSK Changed so it removes an uncommitted mesh
 *      form the owning dataset's meshID list.
 *    - 12/03/03 RM, added more input checking, if mesh is FC_NULL_MESH,
 *         then do nothing and return success 
 *        (issue warning if fc_getLibraryVerbosity()). 
 *    - 9/9/04 WSK changed name to delete. Changed behavior so that
 *      delete always deletes no matter whether on disk or not.
  */
FC_ReturnCode fc_deleteMesh(
 FC_Mesh mesh       /**< input - mesh handle */
) {
  FC_ReturnCode rc;
  int i, j;
  _FC_MeshSlot* meshSlot;
  _FC_DsSlot* dsSlot;
  int numSubset, numVar, numSeqVar, *numStepPerVar;
  FC_Subset* subsets;
  FC_Variable *variables, **seqVariables;

  // special case: deleting a null handle is not an error
  if ( FC_HANDLE_EQUIV(mesh, FC_NULL_MESH)) {
    fc_printfWarningMessage("Deleting FC_NULL_MESH");
    return FC_SUCCESS;
  }

  // check input
  meshSlot = _fc_getMeshSlot(mesh);
  if (meshSlot == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
 
  // log message
  fc_printfLogMessage("Deleting mesh '%s'", meshSlot->header.name);

  // delete member subsets
  rc = fc_getSubsets(mesh, &numSubset, &subsets);
  if (rc != FC_SUCCESS)
    return rc;
  for (i = 0; i < numSubset; i++) {
    rc = fc_deleteSubset(subsets[i]);
    if (rc != FC_SUCCESS)
      return rc;
  }
  free(subsets);
  
  // delete member variables
  rc = fc_getVariables(mesh, &numVar, &variables);
  if (rc != FC_SUCCESS)
    return rc;
  for (i = 0; i < numVar; i++) { 
    rc = fc_deleteVariable(variables[i]);
    if (rc != FC_SUCCESS)
      return rc;
  }
  free(variables);
  
  // delete member seq vars
  rc = fc_getSeqVariables(mesh, &numSeqVar, &numStepPerVar, &seqVariables);
  if (rc != FC_SUCCESS)
    return rc;
  for (i = 0; i < numSeqVar; i++) {
    rc = fc_deleteSeqVariable(numStepPerVar[i], seqVariables[i]);
    if (rc != FC_SUCCESS)
      return rc;
  }
  free(numStepPerVar);
  for (i = 0; i < numSeqVar; i++)
    free(seqVariables[i]);
  free(seqVariables);
  
  // remove reference from dataset
  dsSlot = _fc_getDsSlot(meshSlot->ds);
  if (dsSlot->numMesh == 1) {
    free(dsSlot->meshIDs);
    dsSlot->meshIDs = NULL;
    dsSlot->numMesh = 0;
  }
  else {
    for (i = 0; i < dsSlot->numMesh; i++) {
      if (dsSlot->meshIDs[i] == meshSlot->header.slotID)
        break;
    }
    for (j = i+1; j < dsSlot->numMesh; j++) {
      dsSlot->meshIDs[j - 1] = dsSlot->meshIDs[j];
    }
    dsSlot->numMesh--;
  }

  // clear & delete slot
  return _fc_deleteMeshSlot(mesh);
}

//@}

/** \name Get mesh metadata. */
//-------------------------------
//@{

/**
 * \ingroup  Mesh
 * \brief Check that the handle refers to a valid mesh.
 *
 * \description
 *
 *    Returns 1 (true) is the handle is valid, and 0 (false) if it is not.
 *
 * \modifications  
 *    - 05/25/04 WSK Created.
 */
int fc_isMeshValid(
  FC_Mesh mesh  /**< input - mesh handle */
) {
  return _fc_isHandleValid(mesh.slotID, mesh.uID, meshTableSize,
                           (_FC_SlotHeader**)meshTable);
}

/**
 * \ingroup  Mesh
 * \brief  Return name of a mesh.
 *
 * \description
 *  
 *    Return the name associated with a mesh handle. Spaced
 *    is allocated for the string so it needs to be freed by the
 *    user when finished.
 *
 * \modifications  
 *     - Aug, 9, 2002.  W Koegler  Created.
 */
FC_ReturnCode fc_getMeshName(
  FC_Mesh mesh,         /**< input - mesh handle */
  char **meshname       /**< output - mesh name, space allocated here
                             and must be freed by user when finished */
) {
  char *cp;
  _FC_MeshSlot* meshSlot;
  
  // default
  if (meshname)
    *meshname = NULL;
  
  // check input
  meshSlot = _fc_getMeshSlot(mesh);
  if (meshSlot == NULL || meshname == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Getting name for mesh '%s'", meshSlot->header.name);

  // copy name
  cp = meshSlot->header.name;
  *meshname = malloc(strlen(cp) + 1);
  if (*meshname == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  strcpy(*meshname, cp);
  
  return FC_SUCCESS;
}

/**
 * \ingroup  Mesh
 * \brief  Get the parent dataset of the mesh
 *
 * \modifications  
 *    - Nancy Collins, Created.
 *    - 12/11/03 RM, added default return dataset = \ref FC_NULL_DATASET, 
 *            input check for dataset == NULL 
 */
FC_ReturnCode fc_getDatasetFromMesh(
  FC_Mesh mesh,          /**< input - mesh handle */
  FC_Dataset *dataset    /**< output - dataset handle */
) {
  _FC_MeshSlot* meshSlot;
  
  // set up default value
  if(dataset != NULL)
    *dataset = FC_NULL_DATASET;
  
  // check input
  meshSlot = _fc_getMeshSlot(mesh);
  if (meshSlot == NULL || dataset == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));    
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Getting dataset from mesh '%s'", meshSlot->header.name);
  
  // return requested value
  *dataset = meshSlot->ds;
  
  return FC_SUCCESS;
}

/**
 * \ingroup  Mesh
 * \brief  Get information about a mesh.
 *
 * \description
 *
 *    Returns the topological dimensionality of the mesh elements, the
 *    spatial dimensionality of the mesh coordinates, the number of
 *    vertices, the number of elements, and the element type.
 *
 * \modifications  
 *    - Nancy Collins Created.
 *    - 11/25/03 RM, added checking for NULL mesh handle, 
 *         default return: topodim = -1, dim = -1, numVertex = -1;
 *         numElement = -1, elemType = \ref FC_ET_UNKNOWN; issue an error if 
 *         fc_getLibraryVerbosity().
 */
FC_ReturnCode fc_getMeshInfo(
  FC_Mesh mesh,     /**< input - mesh handle */
  int *topodim,     /**< output - topological dimensionality of mesh (e.g. 0 
                            for a point mesh, 2 for a surface, 3 for volume) */
  int *dim,         /**< output - 2=2d coordinates, 3=3d coords */
  int *numVertex,   /**< output - number of vertices */
  int *numElement,  /**< output - number of elements */
  FC_ElementType *elemType    /**< output - element type (quad, tet) */
) {
  _FC_MeshSlot* meshSlot;
  
  // default values returned if anything goes wrong
  if (topodim != NULL) 
    *topodim = -1;   
  if (dim != NULL) 
    *dim = -1;
  if (numVertex != NULL) 
    *numVertex = -1;
  if (numElement != NULL)
    *numElement = -1;
  if (elemType != NULL) 
    *elemType = FC_ET_UNKNOWN; 
  
  // check input
  meshSlot = _fc_getMeshSlot(mesh);
  if (meshSlot == NULL || 
      (!topodim && !dim && !numVertex && !numElement && !elemType) ) {
     fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Getting info for mesh '%s'", meshSlot->header.name);

  // fill in requested values
  if (dim != NULL)
    *dim = meshSlot->dim;
  if (topodim != NULL)
    *topodim = meshSlot->topodim;
  if (numVertex != NULL)
    *numVertex = meshSlot->numVertex;
  if (numElement != NULL)
    *numElement = meshSlot->numElement;
  if (elemType != NULL)
    *elemType = meshSlot->elemType;
  
  return FC_SUCCESS;
}

/**
 * \ingroup  Mesh
 * \brief  Returns the topological dimensionality of the mesh.
 *
 * \description
 *
 *    The topological dimensionality of the mesh is not the same
 *    as the dimensionality of the space it is embedded in. Instead the
 *    topological dimensionality is dependent on the element type where 
 *    point meshes are 0D, lines meshes are 1D, surface meshes like 
 *    triangles and quads are 2D, and volume meshes like tets, pyramids, 
 *    prisms and hexs are 3D. 
 *
 * \modifications  
 *     - Nancy Collins, Created.
 *     - 12/11/03 RM, added default return dim=-1, and input checking   
 */
FC_ReturnCode fc_getMeshTopodim(
 FC_Mesh mesh,   /**< input - mesh  */
 int *topodim    /**< output - topological dimensionality of elements */
) {
  _FC_MeshSlot* meshSlot;
  
  // set up defaults 
  if (topodim != NULL)
    *topodim = -1; 
  
  // check input
  meshSlot = _fc_getMeshSlot(mesh);
  if (meshSlot == NULL || topodim == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Getting topo dim for mesh '%s'", 
                      meshSlot->header.name);

  // fill in requested value
  *topodim = meshSlot->topodim;
  
  return FC_SUCCESS;
}

/**
 * \ingroup  Mesh
 * \brief  Returns dimensionality of the mesh coordinates.
 *
 * \description
 *  
 *    I.e. 1D, 2D or 3D coordinates.
 *
 * \modifications  
 *    - Nancy Collins, Created.
 *    - 12/11/03 RM, added default return dim=-1, and check input for dim=NULL
 */
FC_ReturnCode fc_getMeshDim(
 FC_Mesh mesh,       /**< input - mesh  */
 int *dim            /**< output - count of dimensions */
) {
  _FC_MeshSlot* meshSlot;
  
  // set up defaults
  if (dim != NULL)
    *dim = -1;
  
  // check input
  meshSlot = _fc_getMeshSlot(mesh);
  if (meshSlot == NULL || dim == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Getting topo dimensionality for mesh '%s'", 
                      meshSlot->header.name);

  // fill in requested value
  *dim = meshSlot->dim;
  
  return FC_SUCCESS;
}

/**
 * \ingroup  Mesh
 * \brief  Return the number of vertices in a mesh
 *
 * \modifications  
 *     - Nancy Collins, Created.
 *     - 12/11/03 RM, added default return numVertex = -1, and input check 
 *               for numVertex = NULL
 */
FC_ReturnCode fc_getMeshNumVertex(
  FC_Mesh mesh,     /**< input - mesh */
  int *numVertex    /**< output - integer count of vertices */
) {
  _FC_MeshSlot* meshSlot;
  
  // set up defaults
  if (numVertex != NULL)
    *numVertex = -1;
  
  // check input
  meshSlot = _fc_getMeshSlot(mesh);
  if (meshSlot == NULL || numVertex == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Getting numVertex for mesh '%s'", 
                      meshSlot->header.name);

  // fill in requested value
  *numVertex = meshSlot->numVertex;
  
  return FC_SUCCESS;
}

/**
 * \ingroup  Mesh
 * \brief  Returns number of elements in a mesh.
 *
 * \modifications  
 *    - Nancy Collins, Created.
 *    - 12/11/03 RM, added default return numElem = -1, and input check 
 *               for numElem = NULL
 */
FC_ReturnCode fc_getMeshNumElement(
  FC_Mesh mesh,    /**< input - mesh */
  int *numElem     /**< output - number of elements in this set */
) {
  _FC_MeshSlot* meshSlot;
  
  // set up defaults
  if (numElem != NULL)
    *numElem = -1;
  
  // check input
  meshSlot = _fc_getMeshSlot(mesh);
  if (meshSlot == NULL || numElem == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // log message
  fc_printfLogMessage("Getting numElement for mesh '%s'", 
                      meshSlot->header.name);

  // fill in requested value
  *numElem = meshSlot->numElement;
  
  return FC_SUCCESS;
}

/**
 * \ingroup  Mesh
 * \brief  Returns type of elements in the mesh.
 *
 * \modifications  
 *    - Nancy Collins, Created.
 *    - 12/11/03 RM, added default return elemType = FC_ET_UNKNOWN, 
 *               and input check for elemType = NULL
 */
FC_ReturnCode fc_getMeshElementType(
  FC_Mesh mesh,              /**< input - mesh */
  FC_ElementType *elemType   /**< output - element type */
) {
  _FC_MeshSlot* meshSlot;
  
  // set up defaults
  if(elemType != NULL)
    *elemType = FC_ET_UNKNOWN;
  
  // check input
  meshSlot = _fc_getMeshSlot(mesh);
  if (meshSlot == NULL || elemType == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Getting dimensionality for mesh '%s'", 
                      meshSlot->header.name);

  // fill in requested value
  *elemType = meshSlot->elemType;
  
  return FC_SUCCESS;
}

/**
 * \ingroup  Mesh
 * \brief  Get information about edges and faces in a mesh.
 *
 * \description
 *
 *    If edges and faces do not yet exist, they are built.
 *
 * \modifications  
 *    - 09/21/04 WSK, created.
 */
FC_ReturnCode fc_getMeshEdgeFaceInfo(
  FC_Mesh mesh,   /**< input - mesh handle */
  int *numEdge,   /**< output - number of edges */
  int *numFace,   /**< output - number of faces */
  FC_ElementType *faceType    /**< output - face type (tri, quad or mixed) */
) {
  FC_ReturnCode rc;
  _FC_MeshSlot* meshSlot;
  
  // default values returned if anything goes wrong
  if (numEdge != NULL) 
    *numEdge = -1;
  if (numFace != NULL)
    *numFace = -1;
  if (faceType != NULL) 
    *faceType = FC_ET_UNKNOWN; 
  
  // check input
  meshSlot = _fc_getMeshSlot(mesh);
  if (meshSlot == NULL || (!numEdge && !numFace && !faceType) ) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Getting edge and face info for mesh '%s'", 
                      meshSlot->header.name);

  // build edges & faces if necessary
  if (numEdge != NULL && meshSlot->numEdge == -1) {
    rc = _fc_buildEdgeConns(mesh);
    if (rc != FC_SUCCESS)
      return rc;
  }
  if (numFace != NULL && meshSlot->numFace == -1) {
    rc = _fc_buildFaceConns(mesh);
    if (rc != FC_SUCCESS)
      return rc;
  }

  // fill in requested values
  if (numEdge != NULL)
    *numEdge = meshSlot->numEdge;
  if (numFace != NULL)
    *numFace = meshSlot->numFace;
  if (faceType != NULL)
    *faceType = fc_getElementTypeFaceType(meshSlot->elemType);
  
  return FC_SUCCESS;
}

/**
 * \ingroup  Mesh
 * \brief  Get number of edges in a mesh.
 *
 * \description
 *
 *    If edges and faces do not yet exist, they are built.
 *
 * \modifications  
 *    - 09/21/04 WSK, created.
 */
FC_ReturnCode fc_getMeshNumEdge(
  FC_Mesh mesh,   /**< input - mesh handle */
  int *numEdge    /**< output - number of edges */
) {
  FC_ReturnCode rc;
  _FC_MeshSlot* meshSlot;
  
  // default values returned if anything goes wrong
  if (numEdge != NULL) 
    *numEdge = -1;
  
  // check input
  meshSlot = _fc_getMeshSlot(mesh);
  if (meshSlot == NULL || !numEdge) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Getting num edge for mesh '%s'", 
                      meshSlot->header.name);

  // build edges & faces if necessary
  if (meshSlot->numEdge == -1) {
    rc = _fc_buildEdgeConns(mesh);
    if (rc != FC_SUCCESS)
      return rc;
  }

  // fill in requested values
  if (numEdge != NULL)
    *numEdge = meshSlot->numEdge;
  
  return FC_SUCCESS;
}

/**
 * \ingroup  Mesh
 * \brief  Get number of faces in a mesh.
 *
 * \description
 *
 *    If edges and faces do not yet exist, they are built.
 *
 * \modifications  
 *    - 09/21/04 WSK, created.
 */
FC_ReturnCode fc_getMeshNumFace(
  FC_Mesh mesh,   /**< input - mesh handle */
  int *numFace    /**< output - number of faces */
) {
  FC_ReturnCode rc;
  _FC_MeshSlot* meshSlot;
  
  // default values returned if anything goes wrong
  if (numFace != NULL)
    *numFace = -1;
  
  // check input
  meshSlot = _fc_getMeshSlot(mesh);
  if (meshSlot == NULL || !numFace) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Getting number of faces for mesh '%s'", 
                      meshSlot->header.name);

  // build faces if necessary
  if (meshSlot->numFace == -1) {
    rc = _fc_buildFaceConns(mesh);
    if (rc != FC_SUCCESS)
      return rc;
  }

  // fill in requested values
  if (numFace != NULL)
    *numFace = meshSlot->numFace;
  
  return FC_SUCCESS;
}

/**
 * \ingroup  Mesh
 * \brief  Get element type of the face in the mesh.
 *
 * \description
 *
 *    If edges and faces do not yet exist, they are built.
 *    Pyramids and prisms may return the type FC_ET_MIXED.
 *
 * \modifications  
 *    - 09/21/04 WSK, created.
 */
FC_ReturnCode fc_getMeshFaceType(
  FC_Mesh mesh,   /**< input - mesh handle */
  FC_ElementType *faceType    /**< output - face type (tri, quad or mixed) */
) {
  _FC_MeshSlot* meshSlot;
  
  // default values returned if anything goes wrong
  if (faceType != NULL) 
    *faceType = FC_ET_UNKNOWN; 
  
  // check input
  meshSlot = _fc_getMeshSlot(mesh);
  if (meshSlot == NULL || !faceType) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Getting edge and face info for mesh '%s'", 
                      meshSlot->header.name);

  // fill in requested values
  if (faceType != NULL)
    *faceType = fc_getElementTypeFaceType(meshSlot->elemType);
  
  return FC_SUCCESS;
}

/**
 * \ingroup  Mesh
 * \brief  Returns number of a type of mesh entity.
 *
 * \description
 *
 *    If you ask for edges or faces and they don't exist yet, they will
 *    be built.
 *
 *    It is an error to ask for the association type FC_AT_WHOLE_DATASET.
 *
 * \modifications  
 *    - 06/07/04 RM, Created.
 *    - 07/16/04 WSK, Changed name.
 */
FC_ReturnCode fc_getMeshNumEntity(
  FC_Mesh mesh,             /**< input - mesh */
  FC_AssociationType assoc, /**< input - association type of entity */
  int *numEntity            /**< output - number of entities in the mesh */
) {
  FC_ReturnCode rc;
  _FC_MeshSlot* meshSlot;
  // set up default return
  if (numEntity != NULL)
    *numEntity = -1;
    
  // check input
  meshSlot = _fc_getMeshSlot(mesh);
  if (meshSlot == NULL || !fc_isAssociationTypeValid(assoc) || 
      assoc == FC_AT_UNKNOWN || numEntity == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // log message
  fc_printfLogMessage("Getting number of entities of type %s for mesh '%s'", 
                      fc_getAssociationTypeText(assoc),
                      meshSlot->header.name);

  // get the number entities on the mesh for given assoc type
  // build edges & faces if necessary
  switch(assoc) {
  case FC_AT_VERTEX:   *numEntity = meshSlot->numVertex;   break;
  case FC_AT_EDGE:     
    if (meshSlot->numEdge == -1) {
      rc = _fc_buildEdgeConns(mesh);
      if (rc != FC_SUCCESS)
        return rc;
    }
    *numEntity = meshSlot->numEdge;     break;
  case FC_AT_FACE:     
    if (meshSlot->numFace == -1) {
      rc = _fc_buildFaceConns(mesh);
      if (rc != FC_SUCCESS)
        return rc;
    }
    *numEntity = meshSlot->numFace;     break; 
  case FC_AT_ELEMENT:  *numEntity = meshSlot->numElement;  break;
  case FC_AT_WHOLE_MESH:    *numEntity = 1;                     break;
  case FC_AT_UNKNOWN:  // fall through to default
  default:
    return FC_INPUT_ERROR;
  }
  
  // should never reach here
  //fc_printfErrorMessage("Should never reach here");
  //return FC_ERROR;
  return FC_SUCCESS;
}

/**
 * \ingroup Mesh
 * \brief Get the element type of the mesh entity
 *
 * \description
 *
 *    This is really only interesting for faces and elements. The association
 *    type of "whole" returns the type of the elements of the mesh.
 *
 *    If you ask for edges or faces and they don't exist yet, they will
 *    be built.
 *
 *    It is an error to ask for the association type FC_AT_WHOLE_DATASET.
 *
 * \modifications
 *    - 07/16/04 WSK Created.
 */
FC_ReturnCode fc_getMeshEntityElementType(
  FC_Mesh mesh,              /**< input - mesh */
  FC_AssociationType assoc,  /**< input - the entity's type */
  int ID,                    /**< input - the ID of the mesh entity */
  FC_ElementType* elemType   /**< output - the element type of the entity */
) {
  FC_ReturnCode rc;
  int num;
  FC_ElementType* faceTypes;

  // default return
  if (elemType)
    *elemType = FC_ET_UNKNOWN;

  // check input
  if (!fc_isMeshValid(mesh) || !fc_isAssociationTypeValid(assoc) || 
      assoc == FC_AT_UNKNOWN || ID < 0 || !elemType)
    return FC_INPUT_ERROR;
  rc = fc_getMeshNumEntity(mesh, assoc, &num);
  if (rc != FC_SUCCESS)
    return rc;
  if (ID >= num)
    return FC_INPUT_ERROR;

  switch (assoc) {
  case FC_AT_VERTEX: *elemType = FC_ET_POINT; break;
  case FC_AT_EDGE:   *elemType = FC_ET_LINE;  break;
  case FC_AT_FACE:
    fc_getMeshFaceTypesPtr(mesh, &faceTypes);
    *elemType = faceTypes[ID];
    break;
  case FC_AT_ELEMENT: // fall through to FC_AT_WHOLE_MESH
  case FC_AT_WHOLE_MESH: fc_getMeshElementType(mesh, elemType); break;
  case FC_AT_WHOLE_DATASET: // fall through to FC_AT_UNKNOWN
  case FC_AT_UNKNOWN:
    return FC_INPUT_ERROR;
  }
  
  return FC_SUCCESS;
}
    
//@}

/** \name Get mesh big data. */
//-------------------------------
//@{

/**
 * \ingroup  Mesh
 * \brief  Return the vertex coordinates of the mesh.
 *
 * \description
 *
 *    Return an array containing the values of each vertex location.  If the
 *    vertices are 2D, the array has X0, Y0, X1, Y1, etc - interleaved.  If the
 *    vertices are 3D, the array has X0, Y0, Z0, X1, Y1, Z1, ...  The user gets
 *    a pointer back to a shared array and should treat it as read only.  To
 *    modify any of the values the user must first make a copy.
 *
 *    Note, you can index to coordinate j of the ith vertex directly as
 *    coords[i*numDim + j];
 *
 * \modifications  
 *    - Nancy Collins, Created.
 */
FC_ReturnCode fc_getMeshCoordsPtr(
  FC_Mesh mesh,      /**< input - mesh */
  double **coords_p  /**< output - vertices. user should treat as read only */
) {
  FC_ReturnCode rc;
  _FC_MeshSlot* meshSlot;

  // default return
  if (coords_p)
    *coords_p = NULL;

  // check input
  meshSlot = _fc_getMeshSlot(mesh);
  if (meshSlot == NULL || coords_p == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Getting coordinates for mesh '%s'", 
                      meshSlot->header.name);


  // lazy fetching of coords
  if (meshSlot->coords == NULL && meshSlot->header.committed == 1) {
    rc = _fc_readMeshCoords(mesh);
    if (rc != FC_SUCCESS) {
      fc_printfErrorMessage("Lazy loading of coords for mesh '%s' failed",
                            meshSlot->header.name);
      return rc;
    }
  }

  // fill in requested value
  *coords_p = meshSlot->coords;
  
  return FC_SUCCESS;
}

/**
 * \ingroup  Mesh
 * \brief  Return the vertex coordinates of the mesh as a variable.
 *
 * \description
 * 
 *    Returns the variable which is the coords of the mesh, creating it on 
 *    the first call to this function (and returning the same handle on
 *    successive calls). The name of the variable will be 'coordinates'.
 *
 *    The official coords variable is treated a little differently from regular
 *    variables.  Calling this routine always provides a handle to the same
 *    variable.  However, making a copy gives back a normal variable with the
 *    same data (it's not an official coords variable).
 *
 *    This variable will not be written to disk when the dataset containing
 *    it is written to disk. If a copy of a mesh with a coords variable 
 *    is made, the copied mesh will have an existing coords variable.
 *
 * \modifications  
 *    - 09/15/04 WSK, Created.
 */
FC_ReturnCode fc_getMeshCoordsAsVariable(
  FC_Mesh mesh,      /**< input - mesh */
  FC_Variable *coords_var /**< output - variable containing coords */
) {
  FC_ReturnCode rc;
  FC_MathType mathType;
  _FC_MeshSlot* meshSlot;
  double* coords_p;

  // default return
  if (coords_var)
    *coords_var = FC_NULL_VARIABLE;

  // check input
  meshSlot = _fc_getMeshSlot(mesh);
  if (meshSlot == NULL || coords_var == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Creating coordinates variable for mesh '%s'", 
                      meshSlot->header.name);

  // If the coords variable already exists, return it
  if (meshSlot->coordVarID > -1) {
    _FC_GET_HANDLE(*coords_var, _fc_getVarSlotFromID(meshSlot->coordVarID));
    return FC_SUCCESS;
  }

  // create the variable
  rc = fc_getMeshCoordsPtr(mesh, &coords_p);
  if (rc != FC_SUCCESS)
    return rc;
  if (meshSlot->dim == 1)
    mathType = FC_MT_SCALAR;
  else
    mathType = FC_MT_VECTOR;
  rc = fc_createVariable(mesh, "coords", coords_var);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_setVariableData(*coords_var, meshSlot->numVertex, meshSlot->dim,
                          FC_AT_VERTEX, mathType, FC_DT_DOUBLE,
                          (void*)coords_p);
  if (rc != FC_SUCCESS)
    return FC_ERROR;

  // set the coord var
  meshSlot->coordVarID = coords_var->slotID;
  
  return FC_SUCCESS;
}

/**
 * \ingroup  Mesh
 * \brief  Return edge connectivities
 *
 * \description
 *  
 *    This routine returns a pointer to the array that
 *    has the IDs of the vertices that make up the edges. The array
 *    is { Edge0Vert0 Edge0Vert1 Edge1Vert0 Edge1Vert1 ... EdgeNVert0 EdgeNVert1 }.
 *    This is a big data array and should NOT be freed by the calling program.
 *
 *    You can index to the jth vertex of the ith edge directly as conns[i*2 + j];
 *
 *    Note that edges may not yet exist for the mesh, and if they don't,
 *    calling this function will build them (which may take some time).
 *
 * \modifications
 *    - 02/06/04 WSK, created.
 *    - 09/21/04 WSK, removed args for non-big data values which can
 *      now be gotten with fc_getMeshEdgeFaceInfo().
 */
FC_ReturnCode fc_getMeshEdgeConnsPtr(
  FC_Mesh mesh,    /**< input - mesh */
  int **conns_p    /**< output - edge to vertex connections, treat read only */
) {
  FC_ReturnCode rc;
  _FC_MeshSlot* meshSlot;
  
  // default return
  if (conns_p)
    *conns_p = NULL;

  // check input
  meshSlot = _fc_getMeshSlot(mesh);
  if (meshSlot == NULL || conns_p == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Getting edge connectivities on mesh '%s'", 
                      meshSlot->header.name);

  // if it doesn't exist, build vertex, element & neighbor stuff
  if (meshSlot->numEdge == -1) {
    rc = _fc_buildEdgeConns(mesh);
    if (rc != FC_SUCCESS)
      return rc;
  }

  // return pointer
  *conns_p = meshSlot->edgeToVertConns;
  
  return FC_SUCCESS;
}

/**
 * \ingroup  Mesh
 * \brief  Return face types
 *
 * \description
 *
 *    This routine returns a pointer to the array that
 *    has the element type of each face.
 *    This is a big data array and should NOT be freed by the calling program.
 *
 *    Note that faces may not yet exist for the mesh, and if they don't,
 *    calling this function will build them (which may take some time).
 *
 * \modifications
 *    - 10/25/05 WSD, created.
 */
FC_ReturnCode fc_getMeshFaceTypesPtr(
  FC_Mesh mesh,    /**< input - mesh */
  FC_ElementType **faceTypes_p    /**< output - face types, treat read only */
) {
  FC_ReturnCode rc;
  _FC_MeshSlot* meshSlot;
  
  // default return
  if (faceTypes_p)
    *faceTypes_p = NULL;

  // check input
  meshSlot = _fc_getMeshSlot(mesh);
  if (meshSlot == NULL || faceTypes_p == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Getting edge connectivities on mesh '%s'", 
                      meshSlot->header.name);

  // if it doesn't exist, build vertex, element & neighbor stuff
  if (meshSlot->numFace == -1) {
    rc = _fc_buildFaceConns(mesh);
    if (rc != FC_SUCCESS)
      return rc;
  }

  // return pointer
  *faceTypes_p = meshSlot->faceTypes;
  
  return FC_SUCCESS;
}

/**
 * \ingroup  Mesh
 * \brief  Return face connectivities
 *
 * \description  
 *
 *    This routine returns an array that has the IDs of the vertices that make
 *    up the faces. In addition there is an array with the number of vertices
 *    per face because some element types have multiple types of faces
 *    (i.e. pyramids and prisms).  For example, if the faces are triangles,
 *    every entry of the numVertPerFace array will be '3' and the connectivity
 *    array will look like { Face0Vert0 Face0Vert1 Face0Vert2 Face1Vert0
 *    ... FaceNVert0 FaceNVert1 FaceNVert2 }. 
 *
 *    If the face type is mixed, the connectivity array will be
 *    maxNumVertPerFace * numFace long and the unused entries (for when
 *    numVertPerFace < maxNumVertPerFace) will be equal to -1.  For example, a
 *    prism mesh's face connectivity might have numVertPerFace = { 3 4 3 ... }
 *    and conns = { Face0Vert0 Face0Vert1 Face0Vert2 -1 Face1Vert0 Face1Vert1
 *    Face1Vert2 Face1Vert3 Face2Vert0 Face2Vert1 Face2Vert2 -1 ... }.
 *    Note you can index directly to face n by conns[n*maxNumVertPerFace].
 *
 *    The number of vertices per face and the connectivity arrays are both
 *    considered to big data arrays and should NOT be freed by the calling
 *    program. Fetching the numVertPerFaceArray is optional.
 * 
 *    Note that faces may not yet exist for the mesh, and if they don't,
 *    calling this function will build them (which may take some time).
 *
 * \modifications
 *    - 02/06/04 WSK, created.
 *    - 09/21/04 WSK, removed args for non-big data values which can
 *      now be gotten with fc_getMeshEdgeFaceInfo().
 */
FC_ReturnCode fc_getMeshFaceConnsPtr(
  FC_Mesh mesh,          /**< input - mesh */
  int** numVertPerFace_p,  /**< output -  (optional) num vert per face, 
                              treat read only */
  int* maxNumVertPerFace_p,  /**< output - (optional) stride length of the
                              conns array */ 
  int **conns_p   /**< output - face to vertex connections, treat read only */
) {
  FC_ReturnCode rc;
  _FC_MeshSlot* meshSlot;
  
  // default return
  if (numVertPerFace_p)
    *numVertPerFace_p = NULL;
  if (maxNumVertPerFace_p)
    *maxNumVertPerFace_p = -1;
  if (conns_p)
    *conns_p = NULL;

  // check input
  meshSlot = _fc_getMeshSlot(mesh);
  if (meshSlot == NULL || !conns_p) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Getting face connectivities for mesh '%s'", 
                      meshSlot->header.name);

  // if it doesn't exist, build vertex, element & neighbor stuff
  if (meshSlot->numFace == -1) {
    rc = _fc_buildFaceConns(mesh);
    if (rc != FC_SUCCESS)
      return rc;
  }

  // return values
  if (numVertPerFace_p)
    *numVertPerFace_p = meshSlot->numVertPerFace;
  if (maxNumVertPerFace_p) {
    FC_ElementType faceType = fc_getElementTypeFaceType(meshSlot->elemType);
    if (faceType == FC_ET_MIXED)
      *maxNumVertPerFace_p = 4; // assume will be mix of tris and quads
    else if (faceType == FC_ET_UNKNOWN)
      *maxNumVertPerFace_p = 0;
    else
      *maxNumVertPerFace_p = meshSlot->numVertPerFace[0];
  } 
  *conns_p = meshSlot->faceToVertConns;
  
  return FC_SUCCESS;
}

/**
 * \ingroup  Mesh
 * \brief  Return element to vertex connectivities for the mesh
 *
 * \description
 *
 *    Return actual element indices.  The return value is a pointer to an array
 *    of integers.  The user must treat this array as read only.  To modify any
 *    values, the user must make a copy first.  If each element contains N
 *    vertices, the first N integers are for element 0, the next N integers are
 *    element 1, etc. To get the jth vertex ID of the ith element index
 *    as conns[i*numVertPrElem + j]. numVertPerElem can be find by calling
 *    fc_getElementTypeNumVertex() or fc_getElementTypeNumEntity().
 *
 * \modifications  
 *    - Nancy Collins, Created.
 */
FC_ReturnCode fc_getMeshElementConnsPtr(
  FC_Mesh mesh,   /**< input - mesh */
  int **conns_p   /**< output - element index buffer, treat read only */
) {
  FC_ReturnCode rc;
  _FC_MeshSlot* meshSlot;
  
  // default return
  if (conns_p)
    *conns_p = NULL;

  // check input
  meshSlot = _fc_getMeshSlot(mesh);
  if (meshSlot == NULL || conns_p == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Getting element connectivities for mesh '%s'", 
                      meshSlot->header.name);

  // lazy fetching of conns
  if (meshSlot->elemToVertConns == NULL && meshSlot->header.committed == 1) {
    rc = _fc_readMeshElemConns(mesh);
    if (rc != FC_SUCCESS) {
      fc_printfErrorMessage("Lazy loading of conns for mesh '%s' failed",
                            meshSlot->header.name);
      return rc;
    }
  }

  // fill in requested value
  *conns_p = meshSlot->elemToVertConns;
  
  return FC_SUCCESS;
}

/**
 * \ingroup  Mesh
 * \brief  Return element to edge connectivities.
 *
 * \description
 *
 *    This routine returns a pointer to the array that has the IDs of the edges
 *    that make up the elements. The array is { Elem0Edge0 Elem0Edge1
 *    ... Elem0EdgeN Elem1Edge0 Elem1Edge1 ... Elem1EdgeN ...... ElemMEdge0
 *    ElemMEdge1 ... ElemMEdgeN }.  This is a big data array and should NOT be
 *    freed by the calling program. To get the id of the jth edge in the ith
 *    element, index by conns[i*numEdgePerElem+j]. numEdgePerElem can be
 *    found by a call to fc_getElementTypeNumEdge() or
 *    fc_getElementTypeNumEntity().
 *
 *    Note that edges may not yet exist for the mesh, and if they don't,
 *    calling this function will build them (which may take some time).
 *
 * \modifications
 *    - 10/26/05 WSD, created.
 */
FC_ReturnCode fc_getMeshElementToEdgeConnsPtr(
  FC_Mesh mesh,    /**< input - mesh */
  int **elemToEdgeConns_p   /**< output - element to edge connections, 
                               treat read only */
) {
  FC_ReturnCode rc;
  _FC_MeshSlot* meshSlot;
  
  // default return
  if (elemToEdgeConns_p)
    *elemToEdgeConns_p = NULL;

  // check input
  meshSlot = _fc_getMeshSlot(mesh);
  if (meshSlot == NULL || elemToEdgeConns_p == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Getting element to edge connectivities on mesh '%s'", 
                      meshSlot->header.name);

  // if it doesn't exist, build vertex, element & neighbor stuff
  if (meshSlot->numEdge == -1) {
    rc = _fc_buildEdgeConns(mesh);
    if (rc != FC_SUCCESS)
      return rc;
  }

  // return pointer
  *elemToEdgeConns_p = meshSlot->elemToEdgeConns;
  
  return FC_SUCCESS;
}

/**
 * \ingroup  Mesh
 * \brief  Return element to face connectivities.
 *
 * \description
 *
 *    This routine returns a pointer to the array that has the IDs of the faces
 *    that make up the elements. The array is { Elem0Face0 Elem0Face1
 *    ... Elem0FaceN Elem1Face0 Elem1Face1 ... Elem1FaceN ...... ElemMFace0
 *    ElemMFace1 ... ElemMFaceN }.  This is a big data array and should NOT be
 *    freed by the calling program. To get the id of the jth face in the ith
 *    element, index by conns[i*numFacePerElem+j]. numFacePerElem can be
 *    found by a call to fc_getElementTypeNumFace() or
 *    fc_getElementTypeNumEntity().
 *
 *    See also fc_getMeshElementFaceOrientsPtr().
 *
 *    Note that faces may not yet exist for the mesh, and if they don't,
 *    calling this function will build them (which may take some time).
 *
 * \modifications
 *    - 10/26/05 WSD, created.
 */
FC_ReturnCode fc_getMeshElementToFaceConnsPtr(
  FC_Mesh mesh,    /**< input - mesh */
  int **elemToFaceConns_p   /**< output - element to face connections, 
                               treat read only */
) {
  FC_ReturnCode rc;
  _FC_MeshSlot* meshSlot;
  
  // default return
  if (elemToFaceConns_p)
    *elemToFaceConns_p = NULL;

  // check input
  meshSlot = _fc_getMeshSlot(mesh);
  if (meshSlot == NULL || elemToFaceConns_p == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Getting element to face connectivities on mesh '%s'", 
                      meshSlot->header.name);

  // if it doesn't exist, build vertex, element & neighbor stuff
  if (meshSlot->numFace == -1) {
    rc = _fc_buildFaceConns(mesh);
    if (rc != FC_SUCCESS)
      return rc;
  }

  // return pointer
  *elemToFaceConns_p = meshSlot->elemToFaceConns;
  
  return FC_SUCCESS;
}

/**
 * \ingroup  Mesh
 * \brief  Return element's face's orientations relative to global faces.
 *
 * \description
 *
 *    By definition, the vertices of the faces of an element are ordered
 *    (counter clockwise) so that the face's normal points out of the element.
 *    In the gobal face array for a mesh, however, the vertices are (somewhat)
 *    arbitrarily ordered and you cannot know for sure the direction of
 *    a face without reference to an element. For example, a face shared by
 *    two elements has a different orientation depending on which element
 *    you are talking about. Let's say that the face was defined so that
 *    it points out of the first element. Then the entry of the faceOrients
 *    array for that face of the first element is 1. The faceOrients value
 *    for that face in the second element is -1. The -1 tells you that you
 *    need to reverse the order of the vertices of that face in order to
 *    treat it as if it had come from the second element. 
 *
 *    This routine returns a pointer to the array that has orientation flags
 *    for the faces that make up the elements. The array is ordered {
 *    Elem0Face0 Elem0Face1 ... Elem0FaceN Elem1Face0 Elem1Face1 ... Elem1FaceN
 *    ...... ElemMFace0 ElemMFace1 ... ElemMFaceN }.  This is a big data array
 *    and should NOT be freed by the calling program. Possible values are 1:
 *    vertices of the face are in order expected for the element; or -1:
 *    vertices of the face are in the reverse order of that expected for the
 *    element. 
 *
 *    (I'm pretty sure that all unique faces will have orient = 1.)
 *
 *    Note that faces may not yet exist for the mesh, and if they don't,
 *    calling this function will build them (which may take some time).
 *
 * \modifications
 *    - 10/26/05 WSD, created.
 */
FC_ReturnCode fc_getMeshElementFaceOrientsPtr(
  FC_Mesh mesh,    /**< input - mesh */
  int **elemFaceOrients_p   /**< output - element's face's orientations, 
                               treat read only */
) {
  FC_ReturnCode rc;
  _FC_MeshSlot* meshSlot;
  
  // default return
  if (elemFaceOrients_p)
    *elemFaceOrients_p = NULL;

  // check input
  meshSlot = _fc_getMeshSlot(mesh);
  if (meshSlot == NULL || elemFaceOrients_p == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Getting element to face connectivities on mesh '%s'", 
                      meshSlot->header.name);

  // if it doesn't exist, build vertex, element & neighbor stuff
  if (meshSlot->numFace == -1) {
    rc = _fc_buildFaceConns(mesh);
    if (rc != FC_SUCCESS)
      return rc;
  }

  // return pointer
  *elemFaceOrients_p = meshSlot->elemFaceOrients;
  
  return FC_SUCCESS;
}

//@}

/** \name Print mesh. */
//-------------------------------
//@{

/**
 * \ingroup  Mesh
 * \brief Print mesh meta data to stdout, and optionally the connections
 *        and/or coordinates.
 *
 * \modifications 
 *    - 11/20/03 WSK Created.
 *    - 10/31/07 ACG optional precision
 */
FC_ReturnCode fc_printMesh(
  FC_Mesh mesh,         /**< input - mesh */
  char* label,          /**< input - label to prepend variable name with, 
                                     can be NULL */
  int print_edgeConns,  /**< input > 0  = print the edge connectivities. */
  int print_faceConns,  /**< input > 0 = print the face connectivities. */
  int print_elemConns,  /**< input > 0 = print the element connectivities. */
  int print_coords      /**< input - 1 = print the coordinates too
			   2 = print the coordinates in exponential notation w 5 digits precision */

) {
  FC_ReturnCode rc;
  int i, j;
  _FC_MeshSlot* meshSlot;
  int topodim, dim, numVertex, numElement, numSubset, numVar, numSeqVar;
  int numEdge, numFace;
  FC_ElementType elemType, faceType;
  int* ebuf, *numVertPerFace, maxNum;
  double* vbuf;

  // special case: printing a null handle is not an error
  if (FC_HANDLE_EQUIV(mesh, FC_NULL_MESH)) {
    printf("Mesh: FC_NULL_MESH\n");
    fflush(NULL);
    return FC_SUCCESS;
  }

  // check input
  meshSlot = _fc_getMeshSlot(mesh);
  if (meshSlot == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // log message
  fc_printfLogMessage("Printing mesh '%s'", meshSlot->header.name);

  // print label & var name
  if (label)
    printf("%s: ", label);
  else
    printf("Mesh: ");
  printf("'%s'\n", meshSlot->header.name);
  fflush(NULL);

  // print meta data
  rc = fc_getMeshInfo(mesh,  &topodim, &dim, &numVertex, &numElement,
                      &elemType);
  if (rc != FC_SUCCESS) {
    fc_printfErrorMessage("Failed to get info for mesh '%s'",
                          meshSlot->header.name);
    return rc;
  }
  rc = fc_getNumSubset(mesh, &numSubset);
  if (rc != FC_SUCCESS) {
    fc_printfErrorMessage("Failed to get subsets for mesh '%s'",
                          meshSlot->header.name);
    return rc;
  }
  rc = fc_getNumVariable(mesh, &numVar);
  if (rc != FC_SUCCESS) {
    fc_printfErrorMessage("Failed to get variables for mesh '%s'",
                          meshSlot->header.name);
    return rc;
  }
  rc = fc_getNumSeqVariable(mesh, &numSeqVar);
  if (rc != FC_SUCCESS) {
    fc_printfErrorMessage("Failed to get seq variables for mesh '%s'",
                          meshSlot->header.name);
    return rc;
  }
  printf("%snumVertex = %d, dim = %d\n", INDENT_STR, numVertex, dim);
  printf("%snumElement = %d, topodim = %d, elemType = %s\n", INDENT_STR,
         numElement, topodim, fc_getElementTypeText(elemType));
  printf("%snumSubset = %d\n", INDENT_STR, numSubset);
  printf("%snumVar = %d, numSeqVar = %d\n", INDENT_STR, numVar, numSeqVar);
  fflush(NULL);
  
  // print edge connectivities
  if (print_edgeConns > 0) {
    if (meshSlot->numEdge < 0) {
      printf("%sEdges have not been built for this mesh\n", INDENT_STR);
      fflush(NULL);
    }
    else {
      rc = fc_getMeshNumEdge(mesh, &numEdge);
      if (rc != FC_SUCCESS) {
        fc_printfErrorMessage("Failed to get numEdge for mesh '%s'",
                               meshSlot->header.name);
        return rc;
      }
      rc = fc_getMeshEdgeConnsPtr(mesh, &ebuf);
      if (rc != FC_SUCCESS) {
        fc_printfErrorMessage("Failed to get edge conns for mesh '%s'",
                              meshSlot->header.name);
        return rc;
      }
      
      // print the conns
      printf("%sEdge connectivities (2 vertices per edge):\n", INDENT_STR);
      for (i = 0; i < numEdge; i++) {
        printf("%s%s%d:  ", INDENT_STR, INDENT_STR, i);
        for (j = 0; j < 2; j++) {
          printf("%6d ", ebuf[i*2 + j]);
          fflush(NULL);
        }
        printf("\n");
        fflush(NULL);
      }
    }
  }
  
  // print face connectivities
  if (print_faceConns > 0) {
    if (meshSlot->numFace < 0) {
      printf("%sFaces have not been built for this mesh\n", INDENT_STR);
      fflush(NULL);
    }
    else {
      rc = fc_getMeshEdgeFaceInfo(mesh, NULL, &numFace, &faceType);
      if (rc != FC_SUCCESS) {
        fc_printfErrorMessage("Failed to get edge/face info for mesh '%s'",
                              meshSlot->header.name);
        return rc;
      }
      rc = fc_getMeshFaceConnsPtr(mesh, &numVertPerFace, &maxNum, &ebuf);
      if (rc != FC_SUCCESS) {
        fc_printfErrorMessage("Failed to get edge conns for mesh '%s'",
                              meshSlot->header.name);
        return rc;
      }
      
      // print the conns
      if (faceType == FC_ET_MIXED)
        printf("%sFace connectivities (faceType = %s, 3 or 4 vertices per "
               "element):\n", INDENT_STR, fc_getElementTypeText(faceType));
      else
        printf("%sFace connectivities (faceType = %s, %d vertices per "
               "element):\n", INDENT_STR, fc_getElementTypeText(faceType), 
               fc_getElementTypeNumVertex(faceType));
      fflush(NULL);
      for (i = 0; i < numFace; i++) {
        printf("%s%s%d:  ", INDENT_STR, INDENT_STR, i);
        for (j = 0; j < numVertPerFace[i]; j++) { 
          printf("%6d ", ebuf[i*maxNum+j]);
          fflush(NULL);
        }
        printf("\n");
        fflush(NULL);
      }
    }
  }
  
  // print element connectivities
  if (print_elemConns > 0) {
    int numVertPerElem = fc_getElementTypeNumVertex(elemType);
    rc = fc_getMeshElementConnsPtr(mesh, &ebuf);
    if (rc != FC_SUCCESS) {
      fc_printfErrorMessage("Failed to get element conns for mesh '%s'",
                            meshSlot->header.name);
      return FC_ERROR;
    }
    
    // print the conns
    printf("%sElement connectivities (%d vertices per element):\n", 
           INDENT_STR, numVertPerElem);
    fflush(NULL);
    for (i = 0; i < numElement; i++) {
      printf("%s%s%d:  ", INDENT_STR, INDENT_STR, i);
      for (j = 0; j < numVertPerElem; j++) {
        printf("%6d ", ebuf[i*numVertPerElem + j]);
        fflush(NULL);
      }
      printf("\n");
      fflush(NULL);
    }
  }
  
  // print coords
  if (print_coords > 0) {
    rc = fc_getMeshCoordsPtr(mesh, &vbuf);
    if (rc != FC_SUCCESS) {
      fc_printfErrorMessage("Failed to get coords for mesh '%s'",
                            meshSlot->header.name);
      return FC_ERROR;
    }
    
    // print the coordinates per vertex
    printf("%sVertex Coordinates:\n", INDENT_STR);
    fflush(NULL);
    for (i = 0; i < numVertex; i++) {
      printf("%s%s%d:  ", INDENT_STR, INDENT_STR, i);
      for (j = 0; j < dim; j++) {
        // see comment on precision in for DOUBLE case in fc_printVariable()
        //printf("%20.15g ", vbuf[i*dim + j]);
	switch (print_coords){
	case 2:
	  printf("%5le ", vbuf[i*dim + j]);
	  break;
	default:
	  printf("%16.11g ", vbuf[i*dim + j]);
	  break;
	}
        fflush(NULL);
      }
      printf("\n");
      fflush(NULL);
    }
  }
  
  return FC_SUCCESS;
}

//@}

/**
 * \ingroup  PrivateMesh
 * \brief  Return parent edges of the vertices.
 *
 * \description
 *
 *    This routine returns pointers to two arrays which hold the number of
 *    parent edges for each vertex and the IDs of the parent edges for
 *    each vertex. The array with the edge IDs is a 2D array. Each vertex's
 *    edge parent IDs are in sorted order. Get an edge parent ID by
 *
 *    edgeParentsPerVert[VertID][i];
 *
 *    where i is 0 ... numEdgePerVert[vertID] - 1. These are big data arrays
 *    and should NOT be freed by the calling program.
 *
 *    Note that edges may not yet exist for the mesh, and if they don't,
 *    calling this function will build them (which may take some time).
 *   
 *    Also note that if there cannot be any edges on the mesh (only happens
 *    for point meshes), NULLs will be returned instead of arrays.
 *
 * \modifications
 *    - 10/27/05 WSD, created.
 */
FC_ReturnCode _fc_getMeshEdgeParentsOfVerticesPtr(
  FC_Mesh mesh,             /**< input - mesh */
  int** numEdgePerVert,     /**< output - array of number of edge parents per
                               vertex */
  int*** edgeParentsPerVert /**< output - array of arrays of edge parent IDs, 
                               treat read only */
) {
  FC_ReturnCode rc;
  _FC_MeshSlot* meshSlot;
  
  // default return
  if (numEdgePerVert)
    *numEdgePerVert = NULL;
  if (edgeParentsPerVert)
    *edgeParentsPerVert = NULL;

  // check input
  meshSlot = _fc_getMeshSlot(mesh);
  if (meshSlot == NULL || numEdgePerVert == NULL ||
      edgeParentsPerVert == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Getting parent edges of the vertices on mesh '%s'", 
                      meshSlot->header.name);

  // if it doesn't exist, build stuff
  if (meshSlot->numEdge < 0 || 
      (meshSlot->numEdge > 0 && meshSlot->numEdgePerVert == NULL)) {
    rc = _fc_buildParents(mesh, FC_AT_VERTEX, FC_AT_EDGE);
    if (rc != FC_SUCCESS)
      return rc;
  }

  // return pointers
  *numEdgePerVert = meshSlot->numEdgePerVert;
  *edgeParentsPerVert = meshSlot->edgeParentsPerVert;
  
  return FC_SUCCESS;
}

/**
 * \ingroup  PrivateMesh
 * \brief  Return parent faces of the vertices.
 *
 * \description
 *
 *    This routine returns pointers to two arrays which hold the number of
 *    parent faces for each vertex and the IDs of the parent faces for
 *    each vertex. The array with the face IDs is a 2D array. Each vertex's
 *    face parent IDs are in sorted order. Get an face parent ID by
 *
 *    faceParentsPerVert[VertID][i];
 *
 *    where i is 0 ... numFacePerVert[vertID] - 1. These are big data arrays
 *    and should NOT be freed by the calling program.
 *
 *    Note that faces may not yet exist for the mesh, and if they don't,
 *    calling this function will build them (which may take some time).
 *
 *    Also note that if there cannot be any faces on the mesh (only happens
 *    for point and line meshes), NULLs will be returned instead of arrays.
 *
 * \modifications
 *    - 10/27/05 WSD, created.
 */
FC_ReturnCode _fc_getMeshFaceParentsOfVerticesPtr(
  FC_Mesh mesh,             /**< input - mesh */
  int** numFacePerVert,     /**< output - array of number of face parents per
                               vertex */
  int*** faceParentsPerVert /**< output - array of arrays of face parent IDs, 
                               treat read only */
) {
  FC_ReturnCode rc;
  _FC_MeshSlot* meshSlot;
  
  // default return
  if (numFacePerVert)
    *numFacePerVert = NULL;
  if (faceParentsPerVert)
    *faceParentsPerVert = NULL;

  // check input
  meshSlot = _fc_getMeshSlot(mesh);
  if (meshSlot == NULL || numFacePerVert == NULL ||
      faceParentsPerVert == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Getting parent faces of the vertices on mesh '%s'", 
                      meshSlot->header.name);

  // if it doesn't exist, build stuff
  if (meshSlot->numFace < 0 || 
      (meshSlot->numFace > 0 && meshSlot->numFacePerVert == NULL)) {
    rc = _fc_buildParents(mesh, FC_AT_VERTEX, FC_AT_FACE);
    if (rc != FC_SUCCESS)
      return rc;
  }

  // return pointers
  *numFacePerVert = meshSlot->numFacePerVert;
  *faceParentsPerVert = meshSlot->faceParentsPerVert;
  
  return FC_SUCCESS;
}

/**
 * \ingroup  PrivateMesh
 * \brief  Return parent elements of the vertices.
 *
 * \description
 *
 *    This routine returns pointers to two arrays which hold the number of
 *    parent elements for each vertex and the IDs of the parent elements for
 *    each vertex. The array with the element IDs is a 2D array. Each vertex's
 *    element parent IDs are in sorted order. Get an element parent ID by
 *
 *    elemParentsPerVert[VertID][i];
 *
 *    where i is 0 ... numElemPerVert[vertID] - 1. These are a big data arrays and
 *    should NOT be freed by the calling program.
 *
 * \modifications
 *    - 10/27/05 WSD, created.
 */
FC_ReturnCode _fc_getMeshElementParentsOfVerticesPtr(
  FC_Mesh mesh,             /**< input - mesh */
  int** numElemPerVert,     /**< output - array of number of elem parents per
                               vertex */
  int*** elemParentsPerVert /**< output - array of arrays of elem parent IDs, 
                               treat read only */
) {
  FC_ReturnCode rc;
  _FC_MeshSlot* meshSlot;
  
  // default return
  if (numElemPerVert)
    *numElemPerVert = NULL;
  if (elemParentsPerVert)
    *elemParentsPerVert = NULL;

  // check input
  meshSlot = _fc_getMeshSlot(mesh);
  if (meshSlot == NULL || numElemPerVert == NULL ||
      elemParentsPerVert == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Getting parent elems of the vertices on mesh '%s'", 
                      meshSlot->header.name);

  // if it doesn't exist, build stuff
  if (!meshSlot->numElemPerVert) {
    rc = _fc_buildParents(mesh, FC_AT_VERTEX, FC_AT_ELEMENT);
    if (rc != FC_SUCCESS)
      return rc;
  }

  // return pointers
  *numElemPerVert = meshSlot->numElemPerVert;
  *elemParentsPerVert = meshSlot->elemParentsPerVert;
  
  return FC_SUCCESS;
}

/**
 * \ingroup  PrivateMesh
 * \brief  Return parent elements of the edges.
 *
 * \description
 *
 *    This routine returns pointers to two arrays which hold the number of
 *    parent elements for each edge and the IDs of the parent elements for
 *    each edge. The array with the element IDs is a 2D array. Each edge's
 *    element parent IDs are in sorted order. Get an element parent ID by
 *
 *    elemParentsPerEdge[EdgeID][i];
 *
 *    where i is 0 ... numElemPerEdge[edgeID] - 1. These are big data arrays
 *    and should NOT be freed by the calling program.
 *
 *    Note that edges may not yet exist for the mesh, and if they don't,
 *    calling this function will build them (which may take some time).
 *
 * \modifications
 *    - 10/27/05 WSD, created.
 */
FC_ReturnCode _fc_getMeshElementParentsOfEdgesPtr(
  FC_Mesh mesh,             /**< input - mesh */
  int** numElemPerEdge,     /**< output - array of number of elem parents per
                               edge */
  int*** elemParentsPerEdge /**< output - array of arrays of elem parent IDs, 
                               treat read only */
) {
  FC_ReturnCode rc;
  _FC_MeshSlot* meshSlot;
  
  // default return
  if (numElemPerEdge)
    *numElemPerEdge = NULL;
  if (elemParentsPerEdge)
    *elemParentsPerEdge = NULL;

  // check input
  meshSlot = _fc_getMeshSlot(mesh);
  if (meshSlot == NULL || numElemPerEdge == NULL ||
      elemParentsPerEdge == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Getting parent elems of the edges on mesh '%s'", 
                      meshSlot->header.name);

  // if it doesn't exist, build stuff
  if (meshSlot->numEdge < 0 || 
      (meshSlot->numEdge > 0 && meshSlot->numElemPerEdge == NULL)) {
    rc = _fc_buildParents(mesh, FC_AT_EDGE, FC_AT_ELEMENT);
    if (rc != FC_SUCCESS)
      return rc;
  }

  // return pointers
  *numElemPerEdge = meshSlot->numElemPerEdge;
  *elemParentsPerEdge = meshSlot->elemParentsPerEdge;
  
  return FC_SUCCESS;
}

/**
 * \ingroup  PrivateMesh
 * \brief  Return parent elements of the faces.
 *
 * \description
 *
 *    This routine returns pointers to two arrays which hold the number of
 *    parent elements for each face and the IDs of the parent elements for
 *    each face. The array with the element IDs is a 2D array. Each face's
 *    element parent IDs are in sorted order. Get an element parent ID by
 *
 *    elemParentsPerFace[FaceID][i];
 *
 *    where i is 0 ... numElemPerFace[faceID] - 1. These are big data arrays
 *    and should NOT be freed by the calling program.
 *
 *    Note that faces may not yet exist for the mesh, and if they don't,
 *    calling this function will build them (which may take some time).
 *
 * \modifications
 *    - 10/27/05 WSD, created.
 */
FC_ReturnCode _fc_getMeshElementParentsOfFacesPtr(
  FC_Mesh mesh,             /**< input - mesh */
  int** numElemPerFace,     /**< output - array of number of elem parents per
                               face */
  int*** elemParentsPerFace /**< output - array of arrays of elem parent IDs, 
                               treat read only */
) {
  FC_ReturnCode rc;
  _FC_MeshSlot* meshSlot;
  
  // default return
  if (numElemPerFace)
    *numElemPerFace = NULL;
  if (elemParentsPerFace)
    *elemParentsPerFace = NULL;

  // check input
  meshSlot = _fc_getMeshSlot(mesh);
  if (meshSlot == NULL || numElemPerFace == NULL ||
      elemParentsPerFace == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Getting parent elems of the faces on mesh '%s'", 
                      meshSlot->header.name);

  // if it doesn't exist, build stuff
  if (meshSlot->numFace < 0 || 
      (meshSlot->numFace > 0 && meshSlot->numElemPerFace == NULL)) {
    rc = _fc_buildParents(mesh, FC_AT_FACE, FC_AT_ELEMENT);
    if (rc != FC_SUCCESS)
      return rc;
  }

  // return pointers
  *numElemPerFace = meshSlot->numElemPerFace;
  *elemParentsPerFace = meshSlot->elemParentsPerFace;
  
  return FC_SUCCESS;
}

/**
 * \ingroup  PrivateMesh
 * \brief  Return vertex neighbors that share an edge.
 *
 * \description
 *
 *    This routine returns pointers to two arrays which hold the number of
 *    vertex neighbors for each vertex and the IDs of the vertex neighbors for
 *    each vertex. The array with the vertex IDs is a 2D array. Each vert's
 *    neighbor IDs are in sorted order. Get an vetex neighbor ID by
 *
 *    vertNeighsViaEdge[vertID][i];
 *
 *    where i is 0 ... numVertNeighsViaEdge[vertID] - 1. These are big data
 *    arrays and should NOT be freed by the calling program.
 *
 *    Note that edges may not yet exist for the mesh, and if they don't,
 *    calling this function will build them (which may take some time).
 *
 *    Also note that if there cannot be any edges on the mesh (only happens
 *    for point meshes), NULLs will be returned instead of arrays.
 *
 * \modifications
 *    - 10/31/05 WSD, created.
 */
FC_ReturnCode _fc_getMeshVertexNeighborsViaEdgePtr(
  FC_Mesh mesh,             /**< input - mesh */
  int** numVertNeighsViaEdge,  /**< output - array of number of vert 
                                neighbors per vert (that share an edge)*/
  int*** vertNeighsViaEdge  /**< output - array of arrays of vert neighbor
                               IDs, treat read only */
) {
  FC_ReturnCode rc;
  _FC_MeshSlot* meshSlot;
  
  // default return
  if (numVertNeighsViaEdge)
    *numVertNeighsViaEdge = NULL;
  if (vertNeighsViaEdge)
    *vertNeighsViaEdge = NULL;

  // check input
  meshSlot = _fc_getMeshSlot(mesh);
  if (meshSlot == NULL || numVertNeighsViaEdge == NULL ||
      vertNeighsViaEdge == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Getting parent elems of the faces on mesh '%s'", 
                      meshSlot->header.name);

  // if it doesn't exist, build stuff
  if (meshSlot->numEdge < 0 || 
      (meshSlot->numEdge > 0 && meshSlot->numVertNeighsViaEdge == NULL)) {
    rc = _fc_buildMeshVertexNeighborsViaEdge(mesh);
    if (rc != FC_SUCCESS)
      return rc;
  }

  // return pointers
  *numVertNeighsViaEdge = meshSlot->numVertNeighsViaEdge;
  *vertNeighsViaEdge = meshSlot->vertNeighsViaEdge;
  
  return FC_SUCCESS;
}

/**
 * \ingroup  PrivateMesh
 * \brief  Return element neighbors that share an vertex.
 *
 * \description
 *
 *    This routine returns pointers to two arrays which hold the number of
 *    element neighbors for each element and the IDs of the element neighbors for
 *    each element. The array with the element IDs is a 2D array. Each elem's
 *    neighbor IDs are in sorted order. Get an vetex neighbor ID by
 *
 *    elemNeighsViaVert[elemID][i];
 *
 *    where i is 0 ... numElemNeighsViaVert[elemID] - 1. These are big data
 *    arrays and should NOT be freed by the calling program.
 *
 * \modifications
 *    - 10/31/05 WSD, created.
 */
FC_ReturnCode _fc_getMeshElementNeighborsViaVertexPtr(
  FC_Mesh mesh,             /**< input - mesh */
  int** numElemNeighsViaVert,  /**< output - array of number of elem 
                                neighbors per elem (that share an vert)*/
  int*** elemNeighsViaVert  /**< output - array of arrays of elem neighbor
                               IDs, treat read only */
) {
  FC_ReturnCode rc;
  _FC_MeshSlot* meshSlot;
  
  // default return
  if (numElemNeighsViaVert)
    *numElemNeighsViaVert = NULL;
  if (elemNeighsViaVert)
    *elemNeighsViaVert = NULL;

  // check input
  meshSlot = _fc_getMeshSlot(mesh);
  if (meshSlot == NULL || numElemNeighsViaVert == NULL ||
      elemNeighsViaVert == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Getting parent elems of the faces on mesh '%s'", 
                      meshSlot->header.name);

  // if it doesn't exist, build stuff
  if (meshSlot->numVertex > 0 && meshSlot->numElemNeighsViaVert == NULL) {
    rc = _fc_buildMeshElementNeighborsViaEntity(mesh, FC_AT_VERTEX);
    if (rc != FC_SUCCESS)
      return rc;
  }

  // return pointers
  *numElemNeighsViaVert = meshSlot->numElemNeighsViaVert;
  *elemNeighsViaVert = meshSlot->elemNeighsViaVert;
  
  return FC_SUCCESS;
}

/**
 * \ingroup  PrivateMesh
 * \brief  Return element neighbors that share an edge.
 *
 * \description
 *
 *    This routine returns pointers to two arrays which hold the number of
 *    element neighbors for each element and the IDs of the element neighbors for
 *    each element. The array with the element IDs is a 2D array. Each elem's
 *    neighbor IDs are in sorted order. Get an vetex neighbor ID by
 *
 *    elemNeighsViaEdge[elemID][i];
 *
 *    where i is 0 ... numElemNeighsViaEdge[elemID] - 1. These are big data
 *    arrays and should NOT be freed by the calling program.
 *
 *    Note that edges may not yet exist for the mesh, and if they don't,
 *    calling this function will build them (which may take some time).
 *
 *    Also note that if there cannot be any edges on the mesh (only happens
 *    for point meshes), NULLs will be returned instead of arrays.
 *
 * \modifications
 *    - 10/31/05 WSD, created.
 */
FC_ReturnCode _fc_getMeshElementNeighborsViaEdgePtr(
  FC_Mesh mesh,             /**< input - mesh */
  int** numElemNeighsViaEdge,  /**< output - array of number of elem 
                                neighbors per elem (that share an edge)*/
  int*** elemNeighsViaEdge  /**< output - array of arrays of elem neighbor
                               IDs, treat read only */
) {
  FC_ReturnCode rc;
  _FC_MeshSlot* meshSlot;
  
  // default return
  if (numElemNeighsViaEdge)
    *numElemNeighsViaEdge = NULL;
  if (elemNeighsViaEdge)
    *elemNeighsViaEdge = NULL;

  // check input
  meshSlot = _fc_getMeshSlot(mesh);
  if (meshSlot == NULL || numElemNeighsViaEdge == NULL ||
      elemNeighsViaEdge == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Getting parent elems of the faces on mesh '%s'", 
                      meshSlot->header.name);

  // if it doesn't exist, build stuff
  if (meshSlot->numEdge < 0 || 
      (meshSlot->numEdge > 0 && meshSlot->numElemNeighsViaEdge == NULL)) {
    rc = _fc_buildMeshElementNeighborsViaEntity(mesh, FC_AT_EDGE);
    if (rc != FC_SUCCESS)
      return rc;
  }

  // return pointers
  *numElemNeighsViaEdge = meshSlot->numElemNeighsViaEdge;
  *elemNeighsViaEdge = meshSlot->elemNeighsViaEdge;
  
  return FC_SUCCESS;
}

/**
 * \ingroup  PrivateMesh
 * \brief  Return element neighbors that share an face.
 *
 * \description
 *
 *    This routine returns pointers to two arrays which hold the number of
 *    element neighbors for each element and the IDs of the element neighbors for
 *    each element. The array with the element IDs is a 2D array. Each elem's
 *    neighbor IDs are in sorted order. Get an vetex neighbor ID by
 *
 *    elemNeighsViaFace[elemID][i];
 *
 *    where i is 0 ... numElemNeighsViaFace[elemID] - 1. These are big data
 *    arrays and should NOT be freed by the calling program.
 *
 *    Note that faces may not yet exist for the mesh, and if they don't,
 *    calling this function will build them (which may take some time).
 *
 *    Also note that if there cannot be any faces on the mesh (only happens
 *    for point and line meshes), NULLs will be returned instead of arrays.
 *
 * \modifications
 *    - 10/31/05 WSD, created.
 */
FC_ReturnCode _fc_getMeshElementNeighborsViaFacePtr(
  FC_Mesh mesh,             /**< input - mesh */
  int** numElemNeighsViaFace,  /**< output - array of number of elem 
                                neighbors per elem (that share an face)*/
  int*** elemNeighsViaFace  /**< output - array of arrays of elem neighbor
                               IDs, treat read only */
) {
  FC_ReturnCode rc;
  _FC_MeshSlot* meshSlot;
  
  // default return
  if (numElemNeighsViaFace)
    *numElemNeighsViaFace = NULL;
  if (elemNeighsViaFace)
    *elemNeighsViaFace = NULL;

  // check input
  meshSlot = _fc_getMeshSlot(mesh);
  if (meshSlot == NULL || numElemNeighsViaFace == NULL ||
      elemNeighsViaFace == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Getting parent elems of the faces on mesh '%s'", 
                      meshSlot->header.name);

  // if it doesn't exist, build stuff
  if (meshSlot->numFace < 0 || 
      (meshSlot->numFace > 0 && meshSlot->numElemNeighsViaFace == NULL)) {
    rc = _fc_buildMeshElementNeighborsViaEntity(mesh, FC_AT_FACE);
    if (rc != FC_SUCCESS)
      return rc;
  }

  // return pointers
  *numElemNeighsViaFace = meshSlot->numElemNeighsViaFace;
  *elemNeighsViaFace = meshSlot->elemNeighsViaFace;
  
  return FC_SUCCESS;
}

/**
 * \ingroup  PrivateMesh
 * \brief Build an arrays of parent IDs.
 *
 * \description
 *
 *     Build arrays that hold the parent IDs for the requested child/parent
 *     combination.
 *              
 * \modifications  
 *      - 11/01/2005 WSD. Creatd to replace building parents lists
 *      in _fc_buildElements(), _fc_buildEdges, & _fc_buildFaces.
 */
FC_ReturnCode _fc_buildParents(
  FC_Mesh mesh,                  /**< The mesh to build stuff on. */
  FC_AssociationType childType,  /**< Child entity type. */
  FC_AssociationType parentType  /**< Parent entity type. */
) {
  FC_ReturnCode rc;
  int i, j;
  int numChild, numParent, *parentToChildConns, *numVertPerFace = NULL;
  int numChildPerParent, maxNumChildPerParent;
  int *numParentPerChild, **parentsPerChild;
   _FC_MeshSlot *meshSlot;
   FC_ElementType elemType;
   FC_SortedIntArray* lists; // a list of parent IDs
  
  // check inputs
  if (!fc_isMeshValid(mesh) || !fc_isAssociationTypeValid(childType) ||
      childType == FC_AT_UNKNOWN || childType == FC_AT_WHOLE_MESH ||
      !fc_isAssociationTypeValid(parentType) || parentType == FC_AT_UNKNOWN || 
      parentType == FC_AT_WHOLE_MESH) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  if (childType >= parentType) {
    fc_printfErrorMessage("child type must be child of parent type");
    return FC_INPUT_ERROR;
  }
  
  // Setup I
  rc = fc_getMeshNumEntity(mesh, parentType, &numParent);
  if (rc != FC_SUCCESS)
    return rc;

  // Special case, if the parent type doesn't exist, then we're done
  if (numParent == 0) 
    return FC_SUCCESS;

  // log message
  meshSlot = _fc_getMeshSlot(mesh);
  fc_printfLogMessage("Building vertex neighbors via edge on mesh '%s'", 
                      meshSlot->header.name);

  // Setup II
  rc = fc_getMeshNumEntity(mesh, childType, &numChild);
  if (rc != FC_SUCCESS)
    return rc;

   rc = fc_getMeshEntityElementType(mesh, parentType, 0, &elemType);
  if (rc != FC_SUCCESS)
    return rc;
 // this # gets overwritten or vert-face case
  maxNumChildPerParent = fc_getElementTypeNumEntity(elemType, childType);
  // do entity specific stuff - conns
  if (childType == FC_AT_VERTEX && parentType == FC_AT_EDGE) 
    rc = fc_getMeshEdgeConnsPtr(mesh, &parentToChildConns);
  else if (childType == FC_AT_VERTEX && parentType == FC_AT_FACE) 
    rc = fc_getMeshFaceConnsPtr(mesh, &numVertPerFace, &maxNumChildPerParent,
                                &parentToChildConns);
  else if (childType == FC_AT_VERTEX && parentType == FC_AT_ELEMENT) 
    rc = fc_getMeshElementConnsPtr(mesh, &parentToChildConns);
  // child = edge, parent = face not supported
  else if (childType == FC_AT_EDGE && parentType == FC_AT_ELEMENT) 
    rc = fc_getMeshElementToEdgeConnsPtr(mesh, &parentToChildConns);
  else if (childType == FC_AT_FACE && parentType == FC_AT_ELEMENT) 
    rc = fc_getMeshElementToFaceConnsPtr(mesh, &parentToChildConns);
  else {
    fc_printfErrorMessage("The combination: child = %s, parent = %s is "
                          "not supported",
                          fc_getAssociationTypeText(childType), 
                          fc_getAssociationTypeText(parentType));
    return FC_INPUT_ERROR;
  }
  if (rc != FC_SUCCESS)
    return rc;
  lists = malloc(numChild*sizeof(FC_SortedIntArray));
  if (lists == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  for (i = 0; i < numChild; i++)
    fc_initSortedIntArray(&lists[i]);

  // for each Parent
  for (i = 0; i < numParent; i++) {
    if (numVertPerFace)
      numChildPerParent = numVertPerFace[i];
    else
      numChildPerParent = maxNumChildPerParent;
    // get parent's children & add parent ID to those children's lists
    for (j = 0; j < numChildPerParent; j++) {
      int ret = fc_addIntToSortedIntArray(&(lists[parentToChildConns[i*maxNumChildPerParent+j]]),
			    i);
      if (ret < 0) {
        fc_printfErrorMessage("Failed to add id to sia");
        return ret;
      }    
    }
  }

  // convert lists to arrays
  numParentPerChild = (int*)malloc(numChild*sizeof(int));
  parentsPerChild = (int**)calloc(numChild, sizeof(int*));
  if (numParentPerChild == NULL || parentsPerChild == NULL) {
    free(numParentPerChild);
    free(parentsPerChild);
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  for (i = 0; i < numChild; i++) {
    rc = fc_convertSortedIntArrayToIntArray(&(lists[i]), &(numParentPerChild[i]),
                                &(parentsPerChild[i]));
    if (rc != FC_SUCCESS)
      return rc;
  }

  // cleanup 
  free(lists);

  // set results on mesh
  if (childType == FC_AT_VERTEX && parentType == FC_AT_EDGE) {
    meshSlot->numEdgePerVert = numParentPerChild;
    meshSlot->edgeParentsPerVert = parentsPerChild;
  }
  else if (childType == FC_AT_VERTEX && parentType == FC_AT_FACE) {
    meshSlot->numFacePerVert = numParentPerChild;
    meshSlot->faceParentsPerVert = parentsPerChild;
  }
  else if (childType == FC_AT_VERTEX && parentType == FC_AT_ELEMENT) {
    meshSlot->numElemPerVert = numParentPerChild;
    meshSlot->elemParentsPerVert = parentsPerChild;
  }
  else if (childType == FC_AT_EDGE && parentType == FC_AT_ELEMENT) {
    meshSlot->numElemPerEdge = numParentPerChild;
    meshSlot->elemParentsPerEdge = parentsPerChild;
  }
  else if (childType == FC_AT_FACE && parentType == FC_AT_ELEMENT) {
    meshSlot->numElemPerFace = numParentPerChild;
    meshSlot->elemParentsPerFace = parentsPerChild;
  }

  return FC_SUCCESS;
}

/**
 * \ingroup PrivateMesh
 * \brief Special data blob for sorting edges
 *
 * \description
 *
 *    Used only in _fc_buildEdges().
 */
typedef struct {
  int endVert;
  int edgeID;
} _EdgeBlob;
 
/**
 * \ingroup  PrivateMesh
 * \brief Build the edge info for the mesh
 *
 * \description
 *  
 *    Finds the unique list of edges for the given set of elements.
 *    It also adds edgeIDs arrays to each element.
 *
 *    Assumes that all elements are the same type.
 *
 * \modifications   
 *    - Wendy Koegler 4/02  Created
 *    - 02/06/04 WSK, changed to return raw data arrays too.
 *    - 07/09/04 WSK, added building of upward and downward connections
 *      to vertices and elements (so could get rid of _fc_guildUpwardsConns).
 */
FC_ReturnCode _fc_buildEdgeConns(
  FC_Mesh mesh              /**< The mesh to build stuff on. */
) {
  FC_ReturnCode rc;
  int i, j, k;
  int numVertex, numElement, *elemConns;
  int numVertPerElem, numEdgePerElem, numEdge;
  int temp_ID;
  int *edgeConns, *elemToEdgeConns;
  _FC_MeshSlot* meshSlot;
  FC_ElementType elemType;
  typedef int VertexIDPair[2];
  VertexIDPair vertexIDs;
  VertexIDPair* map; // map element's vertices to it's edges vertices
  // for 2D or less, edges are traversed in similar order to vertices
  // for 3D, edges traversed in order of faces, then faces edges (of course,
  // ignoring edges already traversed) -- see fc_buildFaces for face order.
  VertexIDPair lineMap[1] =    { { 0, 1 } };
  VertexIDPair triMap[3] =     { { 0, 1 }, { 1, 2 }, { 2, 0 } };
  VertexIDPair quadMap[4] =    { { 0, 1 }, { 1, 2 }, { 2, 3 }, { 3, 0 } };
  // tetMap: 3 from 1st face, 2 from 2nd face, 1 from 3rd face, 0 from base
  VertexIDPair tetMap[6] =     { { 0, 1 }, { 1, 3 }, { 3, 0 }, 
                                 { 1, 2 }, { 2, 3 }, 
                                 { 2, 0 } };
  // pyramidMap: 3 from 1st face, 2 from 2nd &3rd faces, 1 from 4th, 0 base
  VertexIDPair pyramidMap[8] = { { 0, 1 }, { 1, 4 }, { 4, 0 }, 
                                 { 1, 2 }, { 2, 4 }, 
                                 { 2, 3 }, { 3, 4 }, 
                                 { 3, 0 } };
  // prismMap: 4 for 1st face, 3 for 2nd face, 2 for 4rd face, 0 for base & top
  VertexIDPair prismMap[9] =   { { 0, 1 }, { 1, 4 }, { 4, 3 }, { 3, 0 }, 
                                 { 1, 2 }, { 2, 5 }, { 5, 4 }, 
                                 { 2, 0 }, { 3, 5 } };
  // hexMap: 4 for 1st face, 3 for 2nd & 3rd face, 2 for 4th face, 0 for b&t 
  VertexIDPair hexMap[12] =    { { 0, 1 }, { 1, 5 }, { 5, 4 }, { 4, 0 },
                                 { 1, 2 }, { 2, 6 }, { 6, 5 }, 
                                 { 2, 3 }, { 3, 7 }, { 7, 6 }, 
                                 { 3, 0 }, { 4, 7 } };
  // store edges in set of blob arrays, one for each vertex where the vertex
  // it is stored in is the first vertex of the edge after "sorting"
  // A sorted edge has the smaller vertID first
  // some will have multiple edges, some will have none
  FC_SortedBlobArray* edgeLists; // numVertex long
  _EdgeBlob *edgeBlob;

  // log
  meshSlot = _fc_getMeshSlot(mesh);

  // setup I, assume all elements are of same type
  rc = fc_getMeshInfo(mesh, NULL, NULL, &numVertex, &numElement, &elemType);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_getMeshElementConnsPtr(mesh, &elemConns);
  if (rc != FC_SUCCESS)
    return rc;
  
  // handle the element type & pick the vertex map
  switch (elemType) {
  case FC_ET_POINT:
    // no edges for point elements, so no references to edges for verts
    meshSlot->numEdge = 0;    
    return FC_SUCCESS;  
  case FC_ET_LINE:     map = lineMap;     break;
  case FC_ET_TRI:      map = triMap;      break;
  case FC_ET_QUAD:     map = quadMap;     break;
  case FC_ET_TET:      map = tetMap;      break;
  case FC_ET_PYRAMID:  map = pyramidMap;  break;
  case FC_ET_PRISM:    map = prismMap;    break;
  case FC_ET_HEX:      map = hexMap;      break;
  case FC_ET_MIXED:    // fall through to FC_ET_UNKNOWN
  case FC_ET_UNKNOWN:
    fc_printfErrorMessage("can't do element type %s",
                          fc_getElementTypeText(elemType));
    return FC_ERROR;
  }

  // setup II
  numVertPerElem = fc_getElementTypeNumVertex(elemType);
  numEdgePerElem = fc_getElementTypeNumEdge(elemType);
  // assuming that calloc will get me same results as fc_initSortedBlobArray
  edgeLists = (FC_SortedBlobArray*)calloc(numVertex, sizeof(FC_SortedBlobArray));
  if (edgeLists == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }

  // Find unique edges, filling in each element's edgeIDs array as we go
  elemToEdgeConns = malloc(numElement*numEdgePerElem*sizeof(int));
  if (elemToEdgeConns == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  numEdge = 0;
  // for each element:
  for (i = 0; i < numElement; i++) {
    
    // for each edge
    for (j = 0; j < numEdgePerElem; j++) {
      int foundBlob, blobIdx;
      
      // get and sort edge's vertexIDs
      for (k = 0; k < 2; k++)
        vertexIDs[k] = elemConns[i*numVertPerElem+map[j][k]];
      if (vertexIDs[0] > vertexIDs[1]) {
        temp_ID = vertexIDs[0];
        vertexIDs[0] = vertexIDs[1];
        vertexIDs[1] = temp_ID;
      }
      
      // look in edgeList of first vertex for edge with 2nd vertex
      _fc_lookupBlobInSortedBlobArray(&edgeLists[vertexIDs[0]], &vertexIDs[1],
				      fc_intCompare, &foundBlob, &blobIdx);
      if (foundBlob)  {
	edgeBlob = edgeLists[vertexIDs[0]].blobs[blobIdx];
      }
      else { // make an edge
	edgeBlob = (_EdgeBlob*)malloc(sizeof(_EdgeBlob));
        if (edgeBlob == NULL) {
          fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
          return FC_MEMORY_ERROR;
        }  
        edgeBlob->endVert = vertexIDs[1];
        edgeBlob->edgeID = numEdge;
	_fc_addEntryToSortedBlobArray(&edgeLists[vertexIDs[0]], blobIdx, 
				      edgeBlob);
        numEdge++;
      }
      elemToEdgeConns[i*numEdgePerElem+j] = edgeBlob->edgeID;
    }
  }
  
  // now build the edgeConns & edges array by walking array of edge lists
  edgeConns = malloc(2*numEdge*sizeof(int));
  if (edgeConns == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  } 
  for (i = 0; i < numVertex; i++) {
    for (j = 0; j < edgeLists[i].numBlob; j++) {
      edgeBlob = (_EdgeBlob*)edgeLists[i].blobs[j];
      temp_ID = edgeBlob->edgeID;
      edgeConns[2*temp_ID] = i;
      edgeConns[2*temp_ID + 1] = edgeBlob->endVert;
      free(edgeBlob);
    }
    fc_freeSortedBlobArray(&edgeLists[i]);
  } 
  free(edgeLists);

  // set output
  meshSlot->numEdge = numEdge;
  meshSlot->edgeToVertConns = edgeConns;
  meshSlot->elemToEdgeConns = elemToEdgeConns;

  return FC_SUCCESS;
}

/**
 * \ingroup PrivateMesh
 * \brief Special data blob for sorting faces.
 *
 * \description
 *
 *    Used only by _fc_faceNode_compare() in _fc_buildFaces().
 */
typedef struct {
  FC_ElementType faceType;
  int vertexIDs[4];   // faces won't have more than 4 vertices
  int numVertex;
  int faceID;
} _FaceBlob;

/**
 * \ingroup PrivateMesh
 * \brief Compare two faces
 *
 * \description
 *
 *    To be used in this file (mesh.c) only for comparing faces.
 *    Expects faces to already be sorted so that smallest node ID
 *    is first and 2nd smallest vertex ID is next. (The order of the
 *    last 1 or 2 vertices is set by PATRAN ordering conventions.)
 *
 *    This method returns -1 if the first face's vertex IDs
 *    are smaller than the 2nd face's, 0 if the faces are the same;
 *    and 1 if the first face's vertex IDs are greater than the 2nd's.
 *
 *    Note, only the 1st 3 vertices are compared. In well behaved meshes, this
 *    should be suffient to compare faces. (I.e. we are assuming that quads
 *    will share 1, 2 or 4 vertices, but never just 3--if 3 are the same, all 4
 *    will be the same. We are also assuming that tris will not be "inside a
 *    quad" and therefore share 3 vertices.)
 *
 * \modifications
 *    - 10/10/2005 WSD. Promoted. Used to be inside _fc_buildFaces() but static
 *      functions not allowed inside funtions in gcc 4.0.
 */
static int _faceBlob_compare(const void* node1, const void* node2) {
  const _FaceBlob* faceNode1 = (const _FaceBlob*)node1;
  const _FaceBlob* faceNode2 = (const _FaceBlob*)node2;
  
  if (faceNode1->vertexIDs[0] < faceNode2->vertexIDs[0]) 
    return -1;
  else if (faceNode1->vertexIDs[0] > faceNode2->vertexIDs[0])
    return 1;
  else {
    if (faceNode1->vertexIDs[1] < faceNode2->vertexIDs[1])
      return -1;
    else if (faceNode1->vertexIDs[1] > faceNode2->vertexIDs[1])
      return 1;
    else {
      if (faceNode1->vertexIDs[2] < faceNode2->vertexIDs[2])
        return -1;
      else if (faceNode1->vertexIDs[2] > faceNode2->vertexIDs[2])
        return 1;
      else 
        return 0;
    }
  }
};

/**
 * \ingroup  PrivateMesh
 * \brief Build the face info for the mesh
 *
 * \description
 *  
 *    Finds the unique list of faces for the given set of elements.
 *    It also adds faceIDs arrays to each element.
 *    The vertices in each face are ordered so that the "smallest"
 *    edge is first. That is, the smallest ID is first, then the
 *    edge from that vertex that has the smaller 2nd ID is next.
 *
 *    Assumes that if comparing first 3 vertices of two face and
 *    they are the same, it is the same face (follows from assuming 
 *    planar faces or some kind of well behaved mesh).
 *
 *    When the face type is FC_ET_MIXED (i.e. mix of triangles and quads
 *    for pyramid and prism meshes), the returned connectivity array
 *    will be of the size numFace*maxNumVertPerFace and for faces
 *    that have fewer vertices than maxNumVertPerFace, the remaining
 *    entries for that face will have values of -1. Although this
 *    takes up more space, this organization makes it easier to
 *    index to a particular face.
 *
 * \todo Make it so that unique faces will have orient of 1 in their element.
 *
 * \modifications   
 *    - 02/03/04 WSK, Created
 *    - 02/06/04 WSK, changed to return raw data arrays too.
 *    - 07/09/04 WSK, added building of upward and downward connections
 *      to vertices and elements (so could get rid of _fc_guildUpwardsConns).
 */
FC_ReturnCode _fc_buildFaceConns(
  FC_Mesh mesh       /**< The mesh to build stuff on. */
) { 
  FC_ReturnCode rc;
  int i, j, k;
  int numVertex, numElement, *elemConns;
  int numVertPerElem, numFacePerElem, numFace, maxNumVertPerFace;
  int idx, temp_ID, numVertexID;
  int* numVertPerFace, *faceConns, *elemToFaceConns, *elemFaceOrients;
  _FC_MeshSlot* meshSlot;
  FC_ElementType elemType, faceType, *faceTypes;
  typedef int VertexIDSet[4]; // faces can have maximum of 4 vertices
  int temp_vertexIDs[7]; // 2*n -1
  FC_ElementType* faceTypeMap;
  FC_ElementType triFaceTypeMap[1] =     { FC_ET_TRI };
  FC_ElementType quadFaceTypeMap[1] =    { FC_ET_QUAD };
  FC_ElementType tetFaceTypeMap[4] =     { FC_ET_TRI, FC_ET_TRI, 
                                         FC_ET_TRI, FC_ET_TRI };
  FC_ElementType pyramidFaceTypeMap[5] = { FC_ET_TRI, FC_ET_TRI, 
                                         FC_ET_TRI, FC_ET_TRI, FC_ET_QUAD };
  FC_ElementType prismFaceTypeMap[5] =   { FC_ET_QUAD, FC_ET_QUAD, FC_ET_QUAD,
                                         FC_ET_TRI, FC_ET_TRI };
  FC_ElementType hexFaceTypeMap[6] =     { FC_ET_QUAD, FC_ET_QUAD, FC_ET_QUAD,
                                         FC_ET_QUAD, FC_ET_QUAD, FC_ET_QUAD,};
  VertexIDSet* map; // map element's vertex to it's faces vertices
  // Faces use exodus ordering 
  VertexIDSet triMap[1] =     { { 0, 1, 2, -1 } };
  VertexIDSet quadMap[1] =    { { 0, 1, 2, 3 } };
  // tetMap: 3 sides in same order as bottom verts, then base
  VertexIDSet tetMap[4] =     { { 0, 1, 3, -1 }, 
                                { 1, 2, 3, -1 }, 
                                { 0, 3, 2, -1 },
                                { 0, 2, 1, -1 } }; 
  // pyramidMap: 4 sides in same order as bottom verts, then base
  VertexIDSet pyramidMap[8] = { { 0, 1, 4, -1 }, 
                                { 1, 2, 4, -1 },
                                { 2, 3, 4, -1 },
                                { 0, 4, 3, -1 },
                                { 0, 3, 2, 1 } };
  // prismMap: 3 sides in same order as bottom verts, then base & top
  VertexIDSet prismMap[9] =   { { 0, 1, 4, 3 },
                                { 1, 2, 5, 4 }, 
                                { 0, 3, 5, 2 }, 
                                { 0, 2, 1, -1 }, 
                                { 3, 4, 5, -1 } }; 
  // hexMap: 4 sides in same order as bottom verts, then base & top
  VertexIDSet hexMap[6] =     { { 0, 1, 5, 4 }, 
                                { 1, 2, 6, 5 }, 
                                { 2, 3, 7, 6 }, 
                                { 0, 4, 7, 3 }, 
                                { 0, 3, 2, 1 },
                                { 4, 5, 6, 7 } };
  // store faces in set of blob arrays, one for each vertex where the vertex
  // it is stored in is the first vertex of the face after "sorting"
  // A sorted face as the smallest vertID first, and the second smallest
  // vertID 2nd, the traversal order may have to be "flipped" to do this.
  // some will have multiple faces, some will have none
  FC_SortedBlobArray* faceLists;
  _FaceBlob tempFaceBlob, *faceBlob;

  // log
  meshSlot = _fc_getMeshSlot(mesh);

  // setup I, assume all elements are of same type
  rc = fc_getMeshInfo(mesh, NULL, NULL, &numVertex, &numElement, &elemType);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_getMeshElementConnsPtr(mesh, &elemConns);
  if (rc != FC_SUCCESS)
    return rc;

  // handle the element type & pick the vertex map & faceType
  switch (elemType) {
  case FC_ET_POINT:    // fall through to FC_ET_LINE!   
  case FC_ET_LINE:
    // no faces
    meshSlot->numFace = 0;      
    return FC_SUCCESS; 
  case FC_ET_TRI:    map = triMap;    faceTypeMap = triFaceTypeMap;    break;
  case FC_ET_QUAD:   map = quadMap;   faceTypeMap = quadFaceTypeMap;   break;
  case FC_ET_TET:    map = tetMap;    faceTypeMap = tetFaceTypeMap;    break;
  case FC_ET_PYRAMID: map = pyramidMap; faceTypeMap = pyramidFaceTypeMap; break;
  case FC_ET_PRISM:  map = prismMap;  faceTypeMap = prismFaceTypeMap;  break;
  case FC_ET_HEX:    map = hexMap;    faceTypeMap = hexFaceTypeMap;    break;
  case FC_ET_MIXED:    // fall through to FC_ET_UNKNOWN
  case FC_ET_UNKNOWN:
    fc_printfErrorMessage("can't do element type %s",
                          fc_getElementTypeText(elemType));
    return FC_ERROR;
  }
  // if points or lines, don't need to make face list

  // setup II
  numVertPerElem = fc_getElementTypeNumVertex(elemType);
  numFacePerElem = fc_getElementTypeNumFace(elemType);
  // assuming that calloc will get me same results as fc_initSortedBlobArray
  faceLists = (FC_SortedBlobArray*)calloc(numVertex, sizeof(FC_SortedBlobArray));
  if (faceLists == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }

  // Find unique faces, filling in each element's faceIDs array as we go
  //printf("numVertex = %d, numElement = %d\n", numVertex, numElement);
  elemToFaceConns = malloc(numElement*numFacePerElem*sizeof(int));
  elemFaceOrients = malloc(numElement*numFacePerElem*sizeof(int));
  if (elemToFaceConns == NULL || elemFaceOrients == NULL) {
    free(elemToFaceConns);
    free(elemFaceOrients);
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  numFace = 0;
  // for each element:
  for (i = 0; i < numElement; i++) {

    // for each face
    for (j = 0; j < numFacePerElem; j++) {
      int foundBlob, blobIdx;
      numVertexID = fc_getElementTypeNumVertex(faceTypeMap[j]);
      
      // get face's vertexIDs
      for (k = 0; k < numVertexID; k++)
        temp_vertexIDs[k] = elemConns[i*numVertPerElem+map[j][k]];

      // rotate array so smallest ID is first (keep order)
      idx = 0;
      temp_ID = temp_vertexIDs[0];
      for (k = 1; k < numVertexID; k++) {
        if (temp_vertexIDs[k] < temp_ID) {
          idx = k;
          temp_ID = temp_vertexIDs[k];
        }
      }
      for (k = 0; k < idx; k++) 
        temp_vertexIDs[numVertexID + k] = temp_vertexIDs[k];

      // Fill in temporary face
      tempFaceBlob.faceType = faceTypeMap[j];
      tempFaceBlob.numVertex = numVertexID;

      // Make the vertex array, & flip orientation if needed
      if (temp_vertexIDs[idx + 1] < temp_vertexIDs[idx + numVertexID - 1]) {        
        elemFaceOrients[i*numFacePerElem+j] = 1;
        for (k = 0; k < numVertexID; k++)
          tempFaceBlob.vertexIDs[k] = temp_vertexIDs[k+idx];
      }
      else {
        elemFaceOrients[i*numFacePerElem+j] = -1;
        tempFaceBlob.vertexIDs[0] = temp_vertexIDs[idx];
        for (k = 1; k < numVertexID; k++)
          tempFaceBlob.vertexIDs[k] = temp_vertexIDs[idx + numVertexID - k];
      }

      // look for entry in faceLists based on first vertex
      _fc_lookupBlobInSortedBlobArray(&faceLists[tempFaceBlob.vertexIDs[0]],
				      &tempFaceBlob, _faceBlob_compare,
				      &foundBlob, &blobIdx);
      if (foundBlob) {
	faceBlob = faceLists[tempFaceBlob.vertexIDs[0]].blobs[blobIdx];
      }
      else { // make a new face
	faceBlob = (_FaceBlob*)malloc(sizeof(_FaceBlob));
	if (faceBlob == NULL) {
	  fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
          return FC_MEMORY_ERROR;
	}
	*faceBlob = tempFaceBlob; // check me!
        faceBlob->faceID = numFace;
	_fc_addEntryToSortedBlobArray(&faceLists[tempFaceBlob.vertexIDs[0]],
				      blobIdx, faceBlob);
        numFace++;
      }
      elemToFaceConns[i*numFacePerElem+j] = faceBlob->faceID;
    }
  }

  // now build the faceTypes, numVertPerFace, faceConns and faces arrays by walking array
  // of face lists
  faceType = fc_getElementTypeFaceType(elemType);
  maxNumVertPerFace = fc_getElementTypeNumVertex(faceType);
  if (maxNumVertPerFace < 0 && faceType == FC_ET_MIXED)
    maxNumVertPerFace = 4; // only pyramids & prisms have mixed face type
  faceTypes = malloc(numFace*sizeof(FC_ElementType));
  numVertPerFace = malloc(numFace*sizeof(int));
  faceConns = malloc(numFace*maxNumVertPerFace*sizeof(int));
  if (faceTypes == NULL || numVertPerFace == NULL || faceConns == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    free(faceTypes);
    free(numVertPerFace);
    free(faceConns);
    return FC_MEMORY_ERROR;
  }
  for (i = 0; i < numVertex; i++) {
    for (j = 0; j < faceLists[i].numBlob; j++) {
      faceBlob = (_FaceBlob*)faceLists[i].blobs[j];
      temp_ID = faceBlob->faceID;
      faceTypes[temp_ID] = faceBlob->faceType;
      numVertPerFace[temp_ID] = faceBlob->numVertex;
      memcpy(&faceConns[maxNumVertPerFace*temp_ID], faceBlob->vertexIDs,
             faceBlob->numVertex*sizeof(int));
      for (k = faceBlob->numVertex; k < maxNumVertPerFace; k++)
        faceConns[maxNumVertPerFace*temp_ID + k] = -1;
      free(faceBlob);
    }
    fc_freeSortedBlobArray(&faceLists[i]);
  } 
  free(faceLists);

  // set output
  meshSlot->numFace = numFace;
  meshSlot->faceTypes = faceTypes;
  meshSlot->numVertPerFace = numVertPerFace;
  meshSlot->faceToVertConns = faceConns;
  meshSlot->elemToFaceConns = elemToFaceConns;
  meshSlot->elemFaceOrients = elemFaceOrients;

  return FC_SUCCESS;
}

/**
 * \ingroup  PrivateMesh
 * \brief Build the neighbor information for the vertices array
 *
 * \description  
 *
 * \modifications   
 *    - WSK 02/24/04 WSK, created to build quick lookups for vertex neighbors
 *    - WSK 07/15/04 WSK, now only builds "edge" neighbors
 *    - WSK 10/31/05 WSK, takes mesh as argument instead of all the
 *      stuf it needs.
 */
FC_ReturnCode _fc_buildMeshVertexNeighborsViaEdge(
  FC_Mesh mesh   /**< The mesh to build stuff on. */
) {
  FC_ReturnCode rc;
  int i;
  int numVert, numEdge, *edgeConns;
  int* numVertNeighsViaEdge, **vertNeighsViaEdge;
  _FC_MeshSlot *meshSlot;
  FC_SortedIntArray *lists;  // a list of neighbor IDs
  
  // check inputs
  if (!fc_isMeshValid(mesh)) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  meshSlot = _fc_getMeshSlot(mesh);
  fc_printfLogMessage("Building vertex neighbors via edge on mesh '%s'", 
                      meshSlot->header.name);
  
  // setup
  rc = fc_getMeshNumVertex(mesh, &numVert);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_getMeshNumEdge(mesh, &numEdge);
  if (rc != FC_SUCCESS)
    return rc;

  // A little more checking
  if (numEdge < 1)
    return FC_SUCCESS; // Nothing to do

  // More setup
  rc = fc_getMeshEdgeConnsPtr(mesh, &edgeConns);
  if (rc != FC_SUCCESS)
    return rc;

  // for each edge, add each vertex to the other vertex's list
  lists = malloc(numVert*sizeof(FC_SortedIntArray));
  if (lists == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  for (i = 0; i < numVert; i++)
    fc_initSortedIntArray(&lists[i]);
  for (i = 0; i < numEdge; i++) {
    int ret = fc_addIntToSortedIntArray(&lists[edgeConns[i*2+1]], edgeConns[i*2]);
    if (ret < 0) {
      fc_printfErrorMessage("Failed to add id to sia");
      return ret;
    }
    ret = fc_addIntToSortedIntArray(&lists[edgeConns[i*2]], edgeConns[i*2+1]);
    if (ret < 0) {
      fc_printfErrorMessage("Failed to add id to sia");
      return ret;
    }
  }

  // convert neigh lists to arrays & store on mesh
  numVertNeighsViaEdge = (int*)malloc(numVert*sizeof(int));
  vertNeighsViaEdge = (int**)malloc(numVert*sizeof(int*));
  if (numVertNeighsViaEdge == NULL || vertNeighsViaEdge == NULL) {
    free(numVertNeighsViaEdge);
    free(vertNeighsViaEdge);
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  for (i = 0; i < numVert; i++) {
    rc = fc_convertSortedIntArrayToIntArray(&(lists[i]), &(numVertNeighsViaEdge[i]),
                                &(vertNeighsViaEdge[i]));
    if (rc != FC_SUCCESS)
      return rc;
  }
  free(lists);

  // set results on mesh
  meshSlot->numVertNeighsViaEdge = numVertNeighsViaEdge;
  meshSlot->vertNeighsViaEdge = vertNeighsViaEdge;

  return FC_SUCCESS;
}

/**
 * \ingroup PrivateMesh
 * \brief Build element neighbors that share the specified entity type.
 *
 * \description
 *
 *    Computes the big data arrays that store number of element neighbors
 *    per element and the ids of those neighbors.
 *
 *    This is a helper and should only be called internally. It doesn't
 *    have robust checking.
 *
 * \modifications
 *    - 10/31/2005 WSD. Created to replace _fc_buildElementNeighbors.
 *      Since the code is essentially the same for three possible
 *      share entity types, replace 3 specific versions with 1 general one.
 */
FC_ReturnCode _fc_buildMeshElementNeighborsViaEntity(
  FC_Mesh mesh,   /**< The mesh to build stuff on. */
  FC_AssociationType sharedType /**< The type of shared entity */ 
) {
  FC_ReturnCode rc, rc2;
  int i, j;
  int numElem, numEntity, numEntityPerElem, entityID;
  int *elemToEntityConns, *numElemPerEntity, **elemParentsPerEntity;
  int *numElemNeighsViaEntity, **elemNeighsViaEntity;
  _FC_MeshSlot *meshSlot;
  FC_ElementType elemType;
  FC_SortedIntArray list; // list for each element

  // check inputs
  if (!fc_isMeshValid(mesh) || !fc_isAssociationTypeValid(sharedType) || 
      sharedType == FC_AT_UNKNOWN || sharedType >= FC_AT_ELEMENT) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  // FIX?: Check that entity has topodim less than element

  // log message
  meshSlot = _fc_getMeshSlot(mesh);
  fc_printfLogMessage("Building element neighbors via '%s' on mesh '%s'", 
                      fc_getAssociationTypeText(sharedType),
                      meshSlot->header.name);

  // setup
  rc = fc_getMeshInfo(mesh, NULL, NULL, NULL, &numElem, &elemType);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_getMeshNumEntity(mesh, sharedType, &numEntity);
  if (rc != FC_SUCCESS)
    return rc;
  numEntityPerElem = fc_getElementTypeNumEntity(elemType, sharedType);

  // A little more checking
  if (numEntity < 1) 
    return FC_SUCCESS; // Nothing to do

  // More setup 
  // FIX: replace this if wrapper calls become available
  switch(sharedType) {
  case FC_AT_VERTEX:
    rc = fc_getMeshElementConnsPtr(mesh, &elemToEntityConns);
    rc2 = _fc_getMeshElementParentsOfVerticesPtr(mesh, &numElemPerEntity,
                                             &elemParentsPerEntity);
    break;
  case FC_AT_EDGE:
    rc = fc_getMeshElementToEdgeConnsPtr(mesh, &elemToEntityConns);
    rc2 = _fc_getMeshElementParentsOfEdgesPtr(mesh, &numElemPerEntity,
                                           &elemParentsPerEntity);
    break;
  case FC_AT_FACE:
    rc = fc_getMeshElementToFaceConnsPtr(mesh, &elemToEntityConns);
    rc2 = _fc_getMeshElementParentsOfFacesPtr(mesh, &numElemPerEntity,
                                           &elemParentsPerEntity);
    break;
  default :
    ; // nothing
  }
  if (rc != FC_SUCCESS)
    return rc;
  if (rc2 != FC_SUCCESS)
    return rc2;

  // FIX? different algorithm
  // This *might* be faster but will use more memory
  // For each entity, add it's parents to each other's lists
  // for (i = 0; i < numEntity; i++) {
  //   for (j = 0; j < numElemPerEntity[i]; j++) {
  //     for (k = 0; k < numElemPerEntity[i]; k++) {
  //       if (j != k) {
  //          fc_addIntToSortedIntArray(&list[k], elemParentsPerEntity[i][j]);
  //          fc_addIntToSortedIntArray(&list[j], elemParentsPerEntity[i][k]);
  //       }
  //     } 
  //   }
  // }

  // find neighbors
  numElemNeighsViaEntity = (int*)malloc(numElem*sizeof(int));
  elemNeighsViaEntity = (int**)malloc(numElem*sizeof(int*));
  if (numElemNeighsViaEntity == NULL || elemNeighsViaEntity == NULL) {
    free(numElemNeighsViaEntity);
    free(elemNeighsViaEntity);
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  for (i = 0; i < numElem; i++) {
    fc_initSortedIntArray(&list);
    for (j = 0; j < numEntityPerElem; j++) {
      entityID = elemToEntityConns[i*numEntityPerElem+j];
      int ret = fc_addIntArrayToSortedIntArray(&list,
					       numElemPerEntity[entityID],
					       elemParentsPerEntity[entityID],
					       0);
      if (ret < 0) {
	fc_printfErrorMessage("Failed to add entry to sia");
	return ret;    
      }
    }
    fc_deleteIntFromSortedIntArray(&list, i); // don't count this elem
    rc = fc_convertSortedIntArrayToIntArray(&list, &numElemNeighsViaEntity[i],
                                &elemNeighsViaEntity[i]);
    if (rc != FC_SUCCESS)
      return rc;
  }
 
  // Set results on mesh
  // FIX?: pass array pointer to function call?
  switch(sharedType) {
  case FC_AT_VERTEX:
    meshSlot->numElemNeighsViaVert = numElemNeighsViaEntity;
    meshSlot->elemNeighsViaVert = elemNeighsViaEntity;
    break;
  case FC_AT_EDGE:
    meshSlot->numElemNeighsViaEdge = numElemNeighsViaEntity;
    meshSlot->elemNeighsViaEdge = elemNeighsViaEntity;
    break;
  case FC_AT_FACE:
    meshSlot->numElemNeighsViaFace = numElemNeighsViaEntity;
    meshSlot->elemNeighsViaFace = elemNeighsViaEntity;
    break;
  default :
    ; // nothing
  }

  return FC_SUCCESS;
}

/**
 * \ingroup  PrivateMesh
 * \brief  Initialize a mesh slot
 *
 * \modifications  
 *   APR-30-2003  W Koegler  Created. 
 */
void _fc_initMeshSlot(
  _FC_MeshSlot* meshSlot    /**< input - the mtable slot to initialize */
) {
  // Note: tried to use the same order of members as the meshSlot declaration
  if (meshSlot != NULL) {
    _fc_initSlotHeader(&meshSlot->header);
    meshSlot->ds = FC_NULL_DATASET;
    //---
    meshSlot->fileInfo.id = -1;
    //---
    meshSlot->dim = 0;
    meshSlot->topodim = 0;
    meshSlot->numVertex = 0;
    meshSlot->numElement = 0;
    meshSlot->elemType = FC_ET_UNKNOWN;
    //---
    meshSlot->coords = NULL;
    meshSlot->elemToVertConns = NULL;
    //--
    meshSlot->numEdge = -1;
    meshSlot->numFace = -1;
    meshSlot->edgeToVertConns = NULL;
    meshSlot->faceTypes = NULL;
    meshSlot->numVertPerFace = NULL;
    meshSlot->faceToVertConns = NULL;
    meshSlot->elemToEdgeConns = NULL;
    meshSlot->elemToFaceConns = NULL;
    meshSlot->elemFaceOrients = NULL;
    //---
    meshSlot->numEdgePerVert = NULL;
    meshSlot->edgeParentsPerVert = NULL;
    meshSlot->numFacePerVert = NULL;
    meshSlot->faceParentsPerVert = NULL;
    meshSlot->numElemPerVert = NULL;
    meshSlot->elemParentsPerVert = NULL;
    meshSlot->numElemPerEdge = NULL;
    meshSlot->elemParentsPerEdge = NULL;
    meshSlot->numElemPerFace = NULL;
    meshSlot->elemParentsPerFace = NULL;
    //---
    meshSlot->numVertNeighsViaEdge = NULL;
    meshSlot->vertNeighsViaEdge = NULL;
    meshSlot->numElemNeighsViaVert = NULL;
    meshSlot->elemNeighsViaVert = NULL;
    meshSlot->numElemNeighsViaEdge = NULL;
    meshSlot->elemNeighsViaEdge = NULL;
    meshSlot->numElemNeighsViaFace = NULL;
    meshSlot->elemNeighsViaFace = NULL;
    //---
    meshSlot->numSub = 0;
    meshSlot->subIDs = NULL;
    //---
    meshSlot->coordVarID = -1;
    meshSlot->numBasicVar = 0;
    meshSlot->basicVarIDs = NULL;
    meshSlot->numSeqVar = 0;
    meshSlot->numStepPerSeqVar = NULL;
    meshSlot->seqVarIDs = NULL;
  }   
}

/**
 * \ingroup  PrivateMesh
 * \brief  Release non-necessary information in a mesh slot.
 *
 * \description
 *  
 *    Releases the coordinates, connectivities and topology
 *    structures.
 *
 * \todo Make this smarter ....
 *
 * \modifications  
 *   Aug 6, 2002  W Koegler  Created
 */
void _fc_releaseMeshSlot(
  _FC_MeshSlot* meshSlotp    /**< input - mtable pointer */
) {
  int i;

  // release big data
  //---
  if (meshSlotp->coords) {
    free(meshSlotp->coords);
    meshSlotp->coords = NULL;
  }
  if (meshSlotp->elemToVertConns) {
    free(meshSlotp->elemToVertConns);
    meshSlotp->elemToVertConns = NULL;
  }
  //---
  if (meshSlotp->edgeToVertConns) {
    free(meshSlotp->edgeToVertConns);
    meshSlotp->edgeToVertConns = NULL;
  }
  if (meshSlotp->faceTypes) {
    free(meshSlotp->faceTypes);
    meshSlotp->faceTypes = NULL;
  }
  if (meshSlotp->numVertPerFace) {
    free(meshSlotp->numVertPerFace);
    meshSlotp->numVertPerFace = NULL;
  }
  if (meshSlotp->faceToVertConns) {
    free(meshSlotp->faceToVertConns);
    meshSlotp->faceToVertConns = NULL;
  }
  if (meshSlotp->elemToEdgeConns) {
    free(meshSlotp->elemToEdgeConns);
    meshSlotp->elemToEdgeConns = NULL;
  }
  if (meshSlotp->elemToFaceConns) {
    free(meshSlotp->elemToFaceConns);
    meshSlotp->elemToFaceConns = NULL;
  }
  if (meshSlotp->elemFaceOrients) {
    free(meshSlotp->elemFaceOrients);
    meshSlotp->elemFaceOrients = NULL;
  }
  //---
  if (meshSlotp->numEdgePerVert) {
    free(meshSlotp->numEdgePerVert);
    meshSlotp->numEdgePerVert = NULL;
  }
  if (meshSlotp->edgeParentsPerVert) {
    for (i = 0; i < meshSlotp->numVertex; i++)
      free(meshSlotp->edgeParentsPerVert[i]);
    free(meshSlotp->edgeParentsPerVert);
    meshSlotp->edgeParentsPerVert = NULL;
  }
  if (meshSlotp->numFacePerVert) {
    free(meshSlotp->numFacePerVert);
    meshSlotp->numFacePerVert = NULL;
  }
  if (meshSlotp->faceParentsPerVert) {
    for (i = 0; i < meshSlotp->numVertex; i++)
      free(meshSlotp->faceParentsPerVert[i]);
    free(meshSlotp->faceParentsPerVert);
    meshSlotp->faceParentsPerVert = NULL;
  }
  if (meshSlotp->numElemPerVert) {
    free(meshSlotp->numElemPerVert);
    meshSlotp->numElemPerVert = NULL;
  }
  if (meshSlotp->elemParentsPerVert) {
    for (i = 0; i < meshSlotp->numVertex; i++)
      free(meshSlotp->elemParentsPerVert[i]);
    free(meshSlotp->elemParentsPerVert);
    meshSlotp->elemParentsPerVert = NULL;
  }
  if (meshSlotp->numElemPerEdge) {
    free(meshSlotp->numElemPerEdge);
    meshSlotp->numElemPerEdge = NULL;
  }
  if (meshSlotp->elemParentsPerEdge) {
    for (i = 0; i < meshSlotp->numEdge; i++)
      free(meshSlotp->elemParentsPerEdge[i]);
    free(meshSlotp->elemParentsPerEdge);
    meshSlotp->elemParentsPerEdge = NULL;
  }
   if (meshSlotp->numElemPerFace) {
    free(meshSlotp->numElemPerFace);
    meshSlotp->numElemPerFace = NULL;
  }
  if (meshSlotp->elemParentsPerFace) {
    for (i = 0; i < meshSlotp->numFace; i++)
      free(meshSlotp->elemParentsPerFace[i]);
    free(meshSlotp->elemParentsPerFace);
    meshSlotp->elemParentsPerFace = NULL;
  }
  //---
  if (meshSlotp->numVertNeighsViaEdge) {
    free(meshSlotp->numVertNeighsViaEdge);
    meshSlotp->numVertNeighsViaEdge = NULL;
  }
  if (meshSlotp->vertNeighsViaEdge) {
    for (i = 0; i < meshSlotp->numVertex; i++)
      free(meshSlotp->vertNeighsViaEdge[i]);
    free(meshSlotp->vertNeighsViaEdge);
    meshSlotp->vertNeighsViaEdge = NULL;
  }
  if (meshSlotp->numElemNeighsViaVert) {
    free(meshSlotp->numElemNeighsViaVert);
    meshSlotp->numElemNeighsViaVert = NULL;
  }
  if (meshSlotp->elemNeighsViaVert) {
    for (i = 0; i < meshSlotp->numElement; i++)
      free(meshSlotp->elemNeighsViaVert[i]);
    free(meshSlotp->elemNeighsViaVert);
    meshSlotp->elemNeighsViaVert = NULL;
  }
  if (meshSlotp->numElemNeighsViaEdge) {
    free(meshSlotp->numElemNeighsViaEdge);
    meshSlotp->numElemNeighsViaEdge = NULL;
  }
  if (meshSlotp->elemNeighsViaEdge) {
    for (i = 0; i < meshSlotp->numElement; i++)
      free(meshSlotp->elemNeighsViaEdge[i]);
    free(meshSlotp->elemNeighsViaEdge);
    meshSlotp->elemNeighsViaEdge = NULL;
  }
  if (meshSlotp->numElemNeighsViaFace) {
    free(meshSlotp->numElemNeighsViaFace);
    meshSlotp->numElemNeighsViaFace = NULL;
  }
  if (meshSlotp->elemNeighsViaFace) {
    for (i = 0; i < meshSlotp->numElement; i++)
      free(meshSlotp->elemNeighsViaFace[i]);
    free(meshSlotp->elemNeighsViaFace);
    meshSlotp->elemNeighsViaFace = NULL;
  }

  //---
  meshSlotp->numEdge = -1;
  meshSlotp->numFace = -1;
  
  return;
}

/**
 * \ingroup  PrivateMesh
 * \brief  Clear a mesh slot.
 *
 * \description
 *  
 *    Releases all dynamically allocated resources in a mesh,
 *    and reinitializes all members.
 *
 * \modifications  
 *   - Aug 6, 2002  W Koegler  Created
 *   - 2003-APR-31  W Koegler  Made more comprehensive
 */
void _fc_clearMeshSlot(
 _FC_MeshSlot* meshSlotp       /**< input - vtable pointer */
) {
  int i;

  if (meshSlotp != NULL) {
    // release big data
    _fc_releaseMeshSlot(meshSlotp);

    // free other dynamic arrays
    if (meshSlotp->subIDs)
      free(meshSlotp->subIDs);
    if (meshSlotp->basicVarIDs)
      free(meshSlotp->basicVarIDs);
    if (meshSlotp->numStepPerSeqVar)
      free(meshSlotp->numStepPerSeqVar);
    if (meshSlotp->seqVarIDs) {
      for (i = 0; i < meshSlotp->numSeqVar; i++)
	free(meshSlotp->seqVarIDs[i]);
      free(meshSlotp->seqVarIDs);
    }

    // clear table header
    _fc_clearSlotHeader(&meshSlotp->header);

    // reinit
    _fc_initMeshSlot(meshSlotp);
  }
}

/**
 * \ingroup  PrivateMesh
 * \brief  Get the size of the mesh table.
 *
 * \modifications  
 *   - 05/21/04 WSK Created.
 */
int _fc_getMeshTableSize() {
  return meshTableSize;
}

/**
 * \ingroup  PrivateMesh
 * \brief   Get an empty mesh slot from the mesh table.
 *
 * \description
 *  
 *    Will return a new slot or NULL. A NULL means a memory error occurred.
 *
 * \modifications  
 *    - 2003-MAY-01  W Koegler  Created
 *    - 12/08/03 WSK Changed because changed XTables from arrays of slots
 *      to arrays of handles to slots. The _fc_new_XSlot() routines do
 *      not call a common helper routine anymore.
 *    - 5/11/04 WSK  Moved common code from _fc_getNewDsSlot(), 
 *      _fc_getNewSeqSlot(), etc to _fc_getNewSlot().
 */
_FC_MeshSlot* _fc_getNewMeshSlot() {
  _FC_MeshSlot* meshSlot;

  meshSlot = (_FC_MeshSlot*) _fc_getNewSlot(&meshTableSize, 
					    (_FC_SlotHeader***)&meshTable,
				       &meshOpenSlots, sizeof(_FC_MeshSlot));
  _fc_initMeshSlot(meshSlot);

  return meshSlot;
}

/**
 * \ingroup  PrivateMesh
 *
 * \brief Get a mesh slot from a mesh handle.
 *
 * \description
 * 
 *    The slot id of the mesh handle is used to find the
 *    the slot in the mesh table. The identity is then
 *    verified by checking that the slot and the handle's
 *    uIDs are the same. On success, a pointer to the slot 
 *    is returned. Otherwise a NULL pointer is returned.
 *
 * \modifications 
 *   - 2003-OCT-19  WSK  Created to replace using _fc_checkhandle and GET_ID 
 *       for better type checking (using pointer's to slots is less bug
 *       prone than lot's of xxxTable[xxxID] references).
 *   - 05/25/04 WSK Changed to call new _fc_isHandleValid common to all tables.
 */
_FC_MeshSlot* _fc_getMeshSlot(
  FC_Mesh mesh    /**< input - mesh */
) {
  if (_fc_isHandleValid(mesh.slotID, mesh.uID, meshTableSize,
                        (_FC_SlotHeader**)meshTable))
    return meshTable[mesh.slotID];
  else
    return NULL;
}

/**
 * \ingroup  PrivateMesh
 *
 * \brief Get a mesh slot based on ID.
 *
 * \description
 * 
 *    The _FC_MeshSlot at the requested slotID is returned. This should
 *    only be used when you really know that it is the slot you want
 *    because there is no checking of the uID like with _fc_getMeshSlot.
 *    The only checking is that the slotID is exists in the table and
 *    that the uID is > 0 (i.e. it is not empty).
 *
 * \modifications 
 *    - 11/27/03 WSK Created so that external routines can index into
 *      tables without knowing about them.
 *    - 05/25/04 WSK Changed to call new _fc_getSlot() common to all tables.
 */
_FC_MeshSlot* _fc_getMeshSlotFromID(
  int meshID
) {
  return (_FC_MeshSlot*) _fc_getSlot(meshID, meshTableSize, 
                                     (_FC_SlotHeader**)meshTable);
}

/**
 * \ingroup  PrivateMesh
 * \brief Delete the mesh slot associated with the mesh handle.
 *
 * \description 
 *    
 *    This routine clears the dynamic data in the slot and then
 *    deletes the slot from the table.
 *
 * \modifications
 *   - 12/20/2005 WSD. Created.
 */
FC_ReturnCode _fc_deleteMeshSlot(
  FC_Mesh mesh    /**< input - mesh */
) {
  _FC_MeshSlot* meshSlot;
  
  meshSlot = _fc_getMeshSlot(mesh);
  if (meshSlot == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // do it
  _fc_clearMeshSlot(meshSlot);
  return _fc_deleteSlot(meshSlot->header.slotID, meshTableSize, 
			(_FC_SlotHeader**)meshTable, &meshOpenSlots);
}

/**
 * \ingroup  PrivateMesh
 * \brief  Print the contents of the Mesh Table (\ref meshTable) to stderr.
 *
 * \description
 *
 *    Prints the contents of the Mesh Table in a human readable form.
 *    It takes a string to use as a label for this print out of the table.
 *    Pass NULL if you don't want a label.
 *
 * \modifications 
 *    - Nancy Collins Created.
 *    - 2003-NOV-13  WSK  Fixed up.
 */
void _fc_printMeshTable(
  char *label   /**< Input - label for this version the table (without a
                           trailing \\n) pass NULL if no label is desired */
) {
  int i, j, k;
  _FC_MeshSlot* meshSlot;
  
  // print table heading
  fprintf(stderr, TABLE_DIVIDER_STRING);
  fprintf(stderr, "Mesh Table (meshTable):\n");
  if (label != NULL)
    fprintf(stderr, "%s\n", label);
  fprintf(stderr, TABLE_DIVIDER_STRING);
  fflush(NULL);

  // We're done if it's empty
  if (meshTableSize < 1) {
    fprintf(stderr, "(empty)\n");
    fprintf(stderr, "\n");
    fflush(NULL);
    return;
  }

  // print contents for each slot
  for (i = 0; i <  meshTableSize; i++) {
    meshSlot = meshTable[i];
    fprintf(stderr, "%3d: %s\n", i, (meshSlot) ? "exists" : "NULL");
    if (meshSlot == NULL)
      continue;
    _fc_printSlotHeader(meshSlot->header);
    fprintf(stderr, "     FC_Dataset = { %d, %d }\n",
            meshSlot->ds.slotID, meshSlot->ds.uID);
    fprintf(stderr, "     dim = %d, topodim = %d\n", 
            meshSlot->dim, meshSlot->topodim);
    fprintf(stderr, "     numVertex = %d, coords = %s\n",
            meshSlot->numVertex, (meshSlot->coords) ? "exists" : "NULL");
    fprintf(stderr, "     numElement = %d, elemType = %s, elemConns = %s\n",
            meshSlot->numElement, fc_getElementTypeText(meshSlot->elemType), 
            (meshSlot->elemToVertConns) ? "exists" : "NULL");
    fprintf(stderr, "     numEdge = %d, edgeConns = %s, elemToEdgeConns = %s\n",
            meshSlot->numEdge, (meshSlot->edgeToVertConns) ? "exists" : "NULL",
            (meshSlot->elemToEdgeConns) ? "exists" : "NULL");
    fprintf(stderr, "     numFace = %d, faceType = %s, faceTypes = %s\n", 
            meshSlot->numFace,
            fc_getElementTypeText(fc_getElementTypeFaceType(meshSlot->elemType)),
            (meshSlot->faceTypes) ? "exists" : "NULL");
    fprintf(stderr, "         numVertPerFace = %s, faceConns = %s\n",
            (meshSlot->numVertPerFace) ? "exists" : "NULL", 
            (meshSlot->faceToVertConns) ? "exists" : "NULL");
    fprintf(stderr, "         elemToFaceConns = %s, elemFaceOrients = %s\n",
            (meshSlot->elemToFaceConns) ? "exists" : "NULL",
            (meshSlot->elemFaceOrients) ? "exists" : "NULL");
    fprintf(stderr, "     numEdgePerVert = %s, edgeParentsPerVert = %s\n",
            (meshSlot->numEdgePerVert) ? "exists" : "NULL",
            (meshSlot->edgeParentsPerVert) ? "exists" : "NULL");
    fprintf(stderr, "     numFacePerVert = %s, faceParentsPerVert = %s\n",
            (meshSlot->numFacePerVert) ? "exists" : "NULL",
            (meshSlot->faceParentsPerVert) ? "exists" : "NULL");
    fprintf(stderr, "     numElemPerVert = %s, elemParentsPerVert = %s\n",
            (meshSlot->numElemPerVert) ? "exists" : "NULL",
            (meshSlot->elemParentsPerVert) ? "exists" : "NULL");
    fprintf(stderr, "     numElemPerEdge = %s, elemParentsPerEdge = %s\n",
            (meshSlot->numElemPerEdge) ? "exists" : "NULL",
            (meshSlot->elemParentsPerEdge) ? "exists" : "NULL");
    fprintf(stderr, "     numElemPerFace = %s, elemParentsPerFace = %s\n",
            (meshSlot->numElemPerFace) ? "exists" : "NULL",
            (meshSlot->elemParentsPerFace) ? "exists" : "NULL");
    fprintf(stderr, "     numVertNeighsViaEdge = %s, vertNeighsViaEdge = %s\n",
            (meshSlot->numVertNeighsViaEdge) ? "exists" : "NULL", 
            (meshSlot->vertNeighsViaEdge) ? "exists" : "NULL");
    fprintf(stderr, "     numElemNeighsViaVert = %s, elemNeighsViaVert = %s\n",
            (meshSlot->numElemNeighsViaVert) ? "exists" : "NULL", 
            (meshSlot->elemNeighsViaVert) ? "exists" : "NULL");
    fprintf(stderr, "     numElemNeighsViaEdge = %s, elemNeighsViaEdge = %s\n",
            (meshSlot->numElemNeighsViaEdge) ? "exists" : "NULL", 
            (meshSlot->elemNeighsViaEdge) ? "exists" : "NULL");
    fprintf(stderr, "     numElemNeighsViaFace = %s, elemNeighsViaFace = %s\n",
            (meshSlot->numElemNeighsViaFace) ? "exists" : "NULL", 
            (meshSlot->elemNeighsViaFace) ? "exists" : "NULL");
    fprintf(stderr, "     numSub = %d, subIDs = ", meshSlot->numSub);
    if (meshSlot->numSub == 0) 
      fprintf(stderr, "NULL\n");
    else {
      fprintf(stderr, "{ %d", meshSlot->subIDs[0]);
      for (j = 1; j < meshSlot->numSub; j++)
        fprintf(stderr, ", %d", meshSlot->subIDs[j]);
      fprintf(stderr, " }\n");
    }
    fprintf(stderr, "     coordVarID = %d\n", meshSlot->coordVarID); 
    fprintf(stderr, "     numBasicVar = %d, basicVarIDs = ", 
            meshSlot->numBasicVar);
    if (meshSlot->numBasicVar == 0) 
      fprintf(stderr, "NULL\n");
    else {
      fprintf(stderr, "{ %d", meshSlot->basicVarIDs[0]);
      for (j = 1; j < meshSlot->numBasicVar; j++)
        fprintf(stderr, ", %d", meshSlot->basicVarIDs[j]);
      fprintf(stderr, " }\n");
    }
    fprintf(stderr, "     numSeqVar = %d, numStepPerSeqVar = ", meshSlot->numSeqVar);
    if (!meshSlot->numStepPerSeqVar)
      fprintf(stderr, "NULL\n");
    else {
      fprintf(stderr, "{ %d", meshSlot->numStepPerSeqVar[0]);
      for (j = 1; j < meshSlot->numSeqVar; j++)
        fprintf(stderr, ", %d", meshSlot->numStepPerSeqVar[j]);
      fprintf(stderr, " }\n");
    }
    fprintf(stderr, "     seqVarIDs = ");
    if (!meshSlot->seqVarIDs) 
      fprintf(stderr, "NULL\n");
    else {
      fprintf(stderr, "{\n");
      for (j = 0; j < meshSlot->numSeqVar; j++) {
	fprintf(stderr, "         seqVarIDs[%d] = { %d", j, 
                  meshSlot->seqVarIDs[j][0]);
	for (k = 1; k < meshSlot->numStepPerSeqVar[j]; k++)
	  fprintf(stderr, ", %d", meshSlot->seqVarIDs[j][k]);
	fprintf(stderr, " }\n");
      }
      fprintf(stderr, "     }\n");
    }
    fflush(NULL);
  }
  fprintf(stderr, "\n");
  fflush(NULL);

  return;
}

/**
 * \ingroup  PrivateMesh
 * \brief  Free all entries in the mesh table
 *
 * \description
 *
 *    This should only be called after making sure that
 *    all slots have been cleared.
 *
 * \modifications  
 *   - 05/21/04 WSK Created.
 */
void _fc_freeMeshTable() {
  int i;
  for (i = 0; i < meshTableSize; i++) {
    _fc_clearMeshSlot(meshTable[i]);
    free(meshTable[i]);
  }
  free(meshTable);
  fc_freeSortedIntArray(&meshOpenSlots);
  meshTableSize = 0;
  meshTable = NULL;
}
