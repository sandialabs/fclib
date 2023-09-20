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
 * \file subset.c
 * \brief Implementation for the \ref Subset module.
 *
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/modules/subset.c,v $ 
 * $Revision: 1.154 $ 
 * $Date: 2006/10/19 03:14:50 $
 */

// C library includes
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

// fc library includes
#include "base.h"
#include "storage.h"
#include "library.h"
#include "mesh.h"
#include "variable.h"
#include "topo.h"
#include "tableP.h"
#include "datasetP.h"
#include "sequenceP.h"
#include "meshP.h"
#include "variableP.h"
#include "fileioP.h"
#include "util.h"

// this module
#include "subset.h"
#include "subsetP.h"
#include "sequence.h"

/**
 * \addtogroup  Subset
 * \brief Operations on subsets.
 *
 * \description
 *
 *   Sequence subsets are stored as an array of subsets. If a function 
 *   does not have a seqSubset version, in many cases you can call the
 *   subset version on the first varaible in the seqSubset array and 
 *   get the results you expect, e.g., fc_getSubsetInfo(seqSubset[0],...).
 *   In some cases you would need to loop over all subsets in the
 *   seqSubset array, .e.g, fc_getSubsetMembersAsArray(seqSubset[i], ...).
 *   (For an explanation of general data manipulations, see 
 *   \ref DataInterface.)
 */

/**
 * \ingroup Subset
 * \defgroup PrivateSubset (Private)
 */

// ----------- Global variables --------------

/** 
 * \ingroup PrivateSubset
 * \brief The number of allocated slots in the subset table.
 */
static int subTableSize = 0;
/** 
 * \ingroup PrivateSubset
 * \brief The subset table.
 */
static _FC_SubSlot **subTable = NULL;
/**
 * \ingroup PrivateSubset
 * \brief The list of unused slots in the subset table.
 */
static FC_SortedIntArray subOpenSlots = { 0, 0, 0 };

/** \name Create new subset from scratch. */
//-------------------------------
//@{

/**
 * \ingroup Subset
 * \brief  Create a new subset on a specified mesh.
 *
 * \description
 *
 *     Create a new (empty) subset to hold entities of the specified type on a
 *     specified mesh. Sets up subset data structure and initializes it to
 *     default values.  Mesh, subsetName and assoc must be given and valid.
 *     It is an error to try to create a subset on an entity type that does
 *     not exists (e.g. faces on a line mesh).  
 *
 *     It is not possible to create subsets with association
 *     FC_AT_WHOLE_DATASET. 
 *
 * \modifications  
 *    - 03/31/04 RM, Created.
 *    - 05/11/05 WSK, changed to reflect changes to FC_Subset.
 *    - 07/23/05 WSK, removed storage type from argument list.
 */
FC_ReturnCode fc_createSubset(
  FC_Mesh mesh,       /**< input - mesh handle */
  char *subsetName,   /**< input - name of a subset to be created  */
  FC_AssociationType assoc,  /**< input - association of subset entities */
  FC_Subset *subset   /**< output - subset handle */
) {
  FC_ReturnCode rc;
  int num;
  int *tempSlotIDs;
  _FC_MeshSlot* meshSlot;
  _FC_SubSlot* subSlot;

  // default return
  if (subset) 
    *subset = FC_NULL_SUBSET; 
   
  // check input
  meshSlot = _fc_getMeshSlot(mesh);
  if (meshSlot == NULL || !subsetName  || !assoc || 
      !fc_isAssociationTypeValid(assoc) || assoc == FC_AT_UNKNOWN || 
      assoc == FC_AT_WHOLE_DATASET || !subset) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // log message
  fc_printfLogMessage("Creating new subset '%s'", subsetName);

  // do as much checking as possible before creating new slot
  // ----- get number of mask elements depending on association
  rc = fc_getMeshNumEntity(mesh,assoc,&num); 
  if (rc != FC_SUCCESS) 
    return FC_INPUT_ERROR;
  if (num < 1) {
    fc_printfErrorMessage("Entity type %s does not exists on mesh '%s'",
                          fc_getAssociationTypeText(assoc),
                          meshSlot->header.name); 
    return FC_INPUT_ERROR;
  }

  // get an open slot
  subSlot = _fc_getNewSubSlot();
  if (!subSlot) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }

  // table header information
  _FC_GET_HANDLE(*subset, subSlot);
  _fc_setSlotHeaderName(&subSlot->header, subsetName);

  // set back and forth references to mesh
  subSlot->mesh = mesh;
  tempSlotIDs = realloc(meshSlot->subIDs, sizeof(int)*(meshSlot->numSub+1));
  if (tempSlotIDs == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  meshSlot->subIDs = tempSlotIDs;
  meshSlot->subIDs[meshSlot->numSub] = subSlot->header.slotID;
  meshSlot->numSub++;

  // Setup subset metadata
  subSlot->assoc = assoc;
  subSlot->maxNumMember = num;

  return FC_SUCCESS;
}


/**
 * \ingroup Subset
 * \brief  Create a new sequence subset on a specified mesh.
 *
 * \description
 *
 *     Create a new (empty) sequence subset as an array of subsets.
 *     A sequence must be provided in addition to the mesh that
 *     will be the owner of the subset. The sequence and mesh
 *     must be in the same dataset.
 *
 * \modifications  
 *    - 03-03-08 ACG, created
 */
FC_ReturnCode fc_createSeqSubset(
  FC_Mesh mesh,         /**< input - mesh handle */
  FC_Sequence sequence, /**< input - sequence handle */
  char *subsetName,     /**< input - sequence subset name */
  FC_AssociationType assoc,  /**< input - association of subset entities */
  int *numStep,         /**< output -  the number of steps */
  FC_Subset **seqSubs     /**< output - seq subset handles */
) {
  FC_ReturnCode rc;
  _FC_MeshSlot* meshSlot;
  _FC_SeqSlot* seqSlot;
  _FC_SubSlot* subSlot;
  int *tempNumPer, **tempSlotIDs;
  int num;
  int i;

  // default return
  if (numStep)
    *numStep = -1;
  if (seqSubs) 
    *seqSubs = NULL;
   
  // check input
  meshSlot = _fc_getMeshSlot(mesh);
  seqSlot = _fc_getSeqSlot(sequence);
  if (meshSlot == NULL || seqSlot == NULL || subsetName == NULL ||
      !assoc ||  !fc_isAssociationTypeValid(assoc) || assoc == FC_AT_UNKNOWN || 
      assoc == FC_AT_WHOLE_DATASET || numStep == NULL || seqSubs == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  if (!FC_HANDLE_EQUIV(meshSlot->ds, seqSlot->ds)) {
    fc_printfErrorMessage("%s",fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // log message
  fc_printfLogMessage("Creating new seq subset '%s'", subsetName);

  // do as much checking as possible before creating new slot
  // ----- get number of mask elements depending on association
  rc = fc_getMeshNumEntity(mesh,assoc,&num); 
  if (rc != FC_SUCCESS) 
    return FC_INPUT_ERROR;
  if (num < 1) {
    fc_printfErrorMessage("Entity type %s does not exists on mesh '%s'",
                          fc_getAssociationTypeText(assoc),
                          meshSlot->header.name); 
    return FC_INPUT_ERROR;
  }

  //make room for array & setup return values;
  *numStep = seqSlot->numStep;
  *seqSubs = malloc(sizeof(FC_Subset)*seqSlot->numStep);
  if (*seqSubs == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }

  for (i = 0; i < seqSlot->numStep; i++){
    // get an open slot
    subSlot = _fc_getNewSubSlot();
    if (!subSlot) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      return FC_MEMORY_ERROR;
    }

    // table header information
    _FC_GET_HANDLE((*seqSubs)[i], subSlot);
    _fc_setSlotHeaderName(&subSlot->header, subsetName);

    // Setup subset metadata
    subSlot->assoc = assoc;
    subSlot->maxNumMember = num;

    // set back and forth references
    subSlot->mesh = mesh;
    subSlot->sequence = sequence;
    subSlot->stepID = i;
  }

  //set more back and forth references
  tempNumPer = (int*)realloc(meshSlot->numStepPerSeqSub,
			     (meshSlot->numSeqSub+1)*sizeof(int));
  tempSlotIDs = (int**)realloc(meshSlot->seqSubIDs,
			       (meshSlot->numSeqSub+1)*sizeof(int*));
  if (!tempNumPer || ! tempSlotIDs) {
    free(tempNumPer);
    free(tempSlotIDs);
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  tempSlotIDs[meshSlot->numSeqSub] = (int*)malloc(seqSlot->numStep*sizeof(int));
  if (!tempSlotIDs[meshSlot->numSeqSub]) {
    free(tempNumPer);
    free(tempSlotIDs);
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  meshSlot->numStepPerSeqSub = tempNumPer;
  meshSlot->numStepPerSeqSub[meshSlot->numSeqSub] = seqSlot->numStep;
  meshSlot->seqSubIDs = tempSlotIDs;
  for (i = 0; i < seqSlot->numStep; i++) {
    meshSlot->seqSubIDs[meshSlot->numSeqSub][i] = (*seqSubs)[i].slotID;
  }
  meshSlot->numSeqSub++;

  return FC_SUCCESS;
}

//@}

/** \name Other ways to get new subsets. */
//-------------------------------
//@{

/**
 * \ingroup Subset
 * \brief  Copy subset to a new subset.
 *
 * \description
 *  
 *      Create new subset by copying the input subset onto a destination 
 *      mesh. If destination and source mesh is not close enough the 
 *      error will be returned, the new subset will not be created.
 *
 * \modifications
 *    - 05/24/04 RM, created. 
 *    - 06/02/04 RM reworked, added checking compatibility between the 
 *        destination mesh's number of association elements and the 
 *        maxNumMember number of subset 
 */
FC_ReturnCode fc_copySubset(
  FC_Subset subset,      /**< input - handle to subset to be copied        */
  FC_Mesh destMesh,      /**< input - the destination mesh to be copied to */
  char* newSubsetName,   /**< input - new subset's name                    */
  FC_Subset *newSubset   /**< output - handle to new subset                */
){
  FC_ReturnCode rc; 
  int numEntity, maxNumMember, numMember, *memberIDs;
  FC_AssociationType assoc;
  _FC_SubSlot *subSlot;

  // default returns
  if (newSubset)
    *newSubset = FC_NULL_SUBSET;        
 
  // check input - NULL name is o.k.
  if (!fc_isSubsetValid(subset) || !fc_isMeshValid(destMesh) || !newSubset) { 
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // check for null name
  subSlot = _fc_getSubSlot(subset);
  if (newSubsetName == NULL)
    newSubsetName = subSlot->header.name;

  // log message
  fc_printfLogMessage("Copying subset '%s' to '%s'",
                      subSlot->header.name, newSubsetName);

  // get other info from source subset & dest mesh
  rc = fc_getSubsetInfo(subset, NULL, &maxNumMember, &assoc);
  if (rc != FC_SUCCESS) {
    fc_printfErrorMessage("Failed to get info from subset.");
    return FC_ERROR;
  }
  rc = fc_getMeshNumEntity(destMesh, assoc, &numEntity);
  if (rc != FC_SUCCESS) 
    return FC_INPUT_ERROR;

  // more checking - destination mesh and subset assoc type are compatible
  if (numEntity != maxNumMember ) {
    fc_printfErrorMessage("%s: destination mesh not compatible with the subset", 
                          fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // --- do it
  rc = fc_getSubsetMembersAsArray(subset, &numMember, &memberIDs); 
  if (rc != FC_SUCCESS) {
    fc_printfErrorMessage("error in getting members from subset.");
    return FC_ERROR;
  }
  rc = fc_createSubset(destMesh, newSubsetName, assoc, &(*newSubset));
  if ( rc != FC_SUCCESS) {
    fc_printfErrorMessage("error creating subset.");
    return FC_ERROR;
  }
  rc = fc_addArrayMembersToSubset(*newSubset, numMember, memberIDs);
  if ( rc != FC_SUCCESS) {
    fc_printfErrorMessage("error adding members to a new subset.");
    return FC_ERROR;
  }  
  
  // cleanup
  free(memberIDs);

  return FC_SUCCESS;
}


/**
 * \ingroup Subset
 * \brief  Copy seq subset
 *
 * \description
 *  
 *      Create a copy of a seq subset into the specified mesh with the given name.
 *      This will return an erorr if the source and
 *      destination sequences and meshes are not similiar enough.
 *      This will copy a sequence if necessary. The destination mesh
 *      and sequence need to belong to the same dataset of an error
 *      will be returned.
 *
 * \modifications
 *    - 03/05/08 ACG, created. 
 */
FC_ReturnCode fc_copySeqSubset(
  int numStep, /**< input - the number of steps in the sequence variable */ 
  FC_Subset* src_seqSubset, /**< input - the source seq subsets to be copied */
  FC_Mesh dest_mesh,  /**< input - the destination mesh to be copied into */
  FC_Sequence dest_seq, /**< input - the destination seq to be copied into */
  char* newName, /**< input - the name of the new seqsubset, a NULL value is a flag 
		    to use the name of the source seqsubset */
  FC_Subset **new_seqSubset /**< output - array of handles to the new seq subset */
){
  FC_ReturnCode rc; 
  _FC_MeshSlot* dest_meshSlot;
  _FC_SeqSlot* dest_seqSlot;
  _FC_SubSlot* subSlot;
  FC_AssociationType assoc;
  int numEntity, maxNumMember, numMember, *memberIDs;
  int temp_numStep;
  int i; 

  // default returns
  if (new_seqSubset)
    *new_seqSubset =NULL;
 
  // check input - NULL name is o.k.
  dest_meshSlot = _fc_getMeshSlot(dest_mesh);
  dest_seqSlot = _fc_getSeqSlot(dest_seq);
  if (!fc_isSeqSubsetValid(numStep,src_seqSubset) || dest_meshSlot == NULL ||
      dest_seqSlot == NULL || ! FC_HANDLE_EQUIV(dest_meshSlot->ds, dest_seqSlot->ds) ||
      numStep != dest_seqSlot->numStep || new_seqSubset == NULL){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // check for null name
  subSlot = _fc_getSubSlot(src_seqSubset[0]);
  if (newName == NULL)
    newName = subSlot->header.name;

  // log message
  fc_printfLogMessage("Copying seq subset '%s' to '%s'",
                      subSlot->header.name, newName);

  // get other info from source subset & dest mesh
  rc = fc_getSubsetInfo(src_seqSubset[0], NULL, &maxNumMember, &assoc);
  if (rc != FC_SUCCESS) {
    fc_printfErrorMessage("Failed to get info from subset.");
    return FC_ERROR;
  }
  rc = fc_getMeshNumEntity(dest_mesh, assoc, &numEntity);
  if (rc != FC_SUCCESS) 
    return FC_INPUT_ERROR;

  // more checking - destination mesh and subset assoc type are compatible
  if (numEntity != maxNumMember ) {
    fc_printfErrorMessage("%s: destination mesh not compatible with the subset", 
                          fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // --- do it
  rc = fc_createSeqSubset(dest_mesh, dest_seq, newName, assoc, &temp_numStep,
			  new_seqSubset);
  if ( rc != FC_SUCCESS) {
    fc_printfErrorMessage("error creating new seq subset '%s'",newName);
    return rc;
  }

  //get and copy data for each step
  for (i = 0; i < numStep; i++){
    rc = fc_getSubsetMembersAsArray(src_seqSubset[i], &numMember, &memberIDs); 
    if (rc != FC_SUCCESS) {
      fc_printfErrorMessage("error in getting members from subset.");
      fc_deleteSeqSubset(numStep, *new_seqSubset);
      free(new_seqSubset);
      return rc;
    }
    rc = fc_addArrayMembersToSubset((*new_seqSubset)[i], numMember, memberIDs);
    if ( rc != FC_SUCCESS) {
      fc_printfErrorMessage("error adding members to a new subset.");
      fc_deleteSeqSubset(numStep, *new_seqSubset);
      free(new_seqSubset);
      if (memberIDs) free(memberIDs);
      return rc;
    }  
    
    // cleanup
    if (memberIDs) free(memberIDs); //could have been empty subset
  }

  return FC_SUCCESS;
}


/**
 * \ingroup  Subset
 * \brief  Convert an array of basic subsets into a sequence subset
 *
 * \description
 *
 *    This routine converts an array of basic subsets into a sequence subset.
 *    It is a LOT more efficient than copy the data of basic subsets into
 *    the steps of a newly created seq subset because no creation/destruction
 *    goes on -- it just changes how the mesh views the subset handles.
 *
 *    You should put the subsets in the order they will be on
 *    the sequence. There must be the same number of subsets as there
 *    are steps on the sequence. The subsets and the sequence have to be
 *    in the same mesh. All of the subsets must be on the
 *    same mesh and have the same metadata (except the name). The subsets
 *    cannot already be part of a sequence subset.
 *
 *    If no name is supplied, the seqSubset will have the name of the
 *    first subset.
 *
 *    (The array of subsets that is the seq subset should be exactly the
 *    array passed in. A new one is made in case the original array
 *    was automatically generated--fc_deleteSeqSubset would fail
 *    on an automatically allocated array because it tries to free
 *    the array of handles.)
 *
 *    \todo
 *     - Note that unlike the seqVariable version we do not check that
 *       they are on the same dataset. Is that a problem?
 *    - check what happens with the committed flag. since
 *      exodus doesnt store seqSubsets then a seq subset
 *      cannot be commited, but what happens when we convert?
 *      similiar issue with conversion of variables???
 *
 * \modifications  
 *   - 03/05/08 ACG created
 */
FC_ReturnCode fc_convertSubsetsToSeqSubset(
  int numStep,           /**< input - the number of steps */
  FC_Subset* subs,     /**< input - the subsets (array numStep long)
                            (this should be freed by user if malloc'd) */
  FC_Sequence sequence,  /**< input - the sequence (with numStep steps) */
  char* new_name,        /**< input - new name (optional) */
  FC_Subset** seqSubset   /**< output - pointer to seq's subsets handles */
) {
  int i, j;
  int temp_numStep;
  FC_Mesh mesh;
  _FC_MeshSlot* meshSlot;
  _FC_SubSlot *subSlots[numStep];
  int* subMask;
  int *tempNumPer, **tempSlotIDs;
  int *numBasicSub_p, **basicSubIDs_p;
  int *numSeqSub_p, **numStepPerSeqSub_p, ***seqSubIDs_p;

  // default return
  if (seqSubset)
    *seqSubset = NULL;

  // test input -- a lot of conditions!
  if (numStep < 0 || !subs || !fc_isSequenceValid(sequence) ||
      !seqSubset) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  fc_getSequenceNumStep(sequence, &temp_numStep);
  if (numStep != temp_numStep) {
    fc_printfErrorMessage("numStep (%d) does not agree with sequence's "
                          "numStep (%d)", numStep, temp_numStep);
    return FC_INPUT_ERROR;
  }
  for (i = 0; i < numStep; i++) {
    if (!fc_isSubsetValid(subs[i])) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
      return FC_INPUT_ERROR;
    }
    subSlots[i] = _fc_getSubSlot(subs[i]);
    if (subSlots[i]->stepID >= 0) {
      fc_printfErrorMessage("Subset is already part of a sequence subset");
      return FC_INPUT_ERROR;
    }
  }
  fc_getMeshFromSubset(subs[0], &mesh);
  for (i = 1; i < numStep; i++) {
    if (!FC_HANDLE_EQUIV(subSlots[i]->mesh, mesh)) {
      fc_printfErrorMessage("subsets are not on same mesh");
      return FC_INPUT_ERROR;
    }
    if (subSlots[i]->assoc        != subSlots[0]->assoc ||
	subSlots[i]->maxNumMember        != subSlots[0]->maxNumMember) {
      fc_printfErrorMessage("variables have differing metadata");
      return FC_INPUT_ERROR;
    }
  }

  // similar to fc_createSeqSubset()

  // Change names
  if (new_name) {
    for (i = 0; i < numStep; i++) 
      _fc_setSlotHeaderName(&subSlots[i]->header, new_name);
  }
  else {
    for (i = 1; i < numStep; i++)
      _fc_setSlotHeaderName(&subSlots[i]->header, subSlots[0]->header.name);
  }

  // change sub contents to reflect new state as seq subs
  for (i = 0; i < numStep; i++) {
    subSlots[i]->sequence = sequence;
    subSlots[i]->stepID = i;
  }

  // get pointers to stuff on mesh depending
  meshSlot = _fc_getMeshSlot(subSlots[0]->mesh);
  numBasicSub_p = &meshSlot->numSub;
  basicSubIDs_p = &meshSlot->subIDs;
  numSeqSub_p = &meshSlot->numSeqSub;
  numStepPerSeqSub_p = &meshSlot->numStepPerSeqSub;
  seqSubIDs_p = &meshSlot->seqSubIDs;

 // Add the handles to list of sequence vars
  tempNumPer = (int*)realloc(*numStepPerSeqSub_p,
			     ((*numSeqSub_p)+1)*sizeof(int)); 
  tempSlotIDs = (int**)realloc(*seqSubIDs_p,
			       (*numSeqSub_p+1)*sizeof(int*));
  if (!tempNumPer || !tempSlotIDs) {
    free(tempNumPer); free(tempSlotIDs);
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  tempSlotIDs[*numSeqSub_p] = (int*)malloc(numStep*sizeof(int)); 
  if (!tempSlotIDs[*numSeqSub_p]) {
    free(tempNumPer); free(tempSlotIDs);
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  tempNumPer[*numSeqSub_p] = numStep;
  for (i = 0; i < numStep; i++) 
    tempSlotIDs[*numSeqSub_p][i] = subs[i].slotID;
  *numStepPerSeqSub_p = tempNumPer;
  *seqSubIDs_p = tempSlotIDs;
  *numSeqSub_p = (*numSeqSub_p) + 1;

  // remove handles from list of basic subs
  subMask = calloc(*numBasicSub_p, sizeof(int));
  if (subMask == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  for (i = 0; i < numStep; i++) {
    for (j = 0; j < *numBasicSub_p; j++) {
      if (subs[i].slotID == (*basicSubIDs_p)[j]) {
        subMask[j] = 1;
        break;
      }
    }
    if (j == *numBasicSub_p) {
      printf("**** developer note: something went wrong!!\n");
      fflush(NULL);
      return FC_ERROR;
    }
  }
  j = 0;
  *numBasicSub_p = *numBasicSub_p - numStep;
  if (*numBasicSub_p == 0) {
    free(*basicSubIDs_p);
    *basicSubIDs_p = NULL;
  }
  else {
    for (i = 0; i < *numBasicSub_p + numStep; i++) {
      if (subMask[i] == 0) {
        (*basicSubIDs_p)[j] = (*basicSubIDs_p)[i];
        j++;
      }
    }
  }
  free(subMask);
  // FIX realloc to make smaller

  // make return values
  *seqSubset = (FC_Subset*)malloc(numStep*sizeof(FC_Subset));
  for (i = 0; i < numStep; i++)
    (*seqSubset)[i] = subs[i];
    
  return FC_SUCCESS;
}


/**
 * \ingroup  Subset
 * \brief  Convert a sequence subset into an array of normal subsets.
 *
 * \description
 *
 *    This routine converts a sequence subset into an array of basic subsets.
 *    It is a LOT more efficient than copy the data of a seq subset into
 *    the steps of newly created basic subsets because no creation/destruction
 *    goes on -- it just changes how the mesh views the subset handles.
 *
 *    The subsets will come out in the step order. Number of output subsets
 *    will be the number of steps in the original seqSubset. If a name is
 *    passed in the subsets will be named as the name postpended with "_#" where
 *    the number is the step number. Otherwise they will be named with the
 *    current name similiarly postpended. Note that if you are using this
 *    method for conversions to write to exodus, be sure that any name you
 *    use if short enough so that the name is under the exodus name character
 *    limit, or else the numbering will be lost.
 *
 * \todo
 *   - this has max length of digit postpending as 40 char - actaully
 *     calc the powers of 10 to set this correctly
 *
 * \modifications  
 *   - 03/13/08 ACG created
 */
FC_ReturnCode fc_convertSeqSubsetToSubsets(
  int numStep,             /**< input - the number of new steps in the seqSubset */
  FC_Subset* seqSubset,     /**< input - the seqSubset */
  char* new_name,        /**< input - new name (optional) */
  FC_Subset** subsets    /**< output - array of subsets */  
) {
  _FC_MeshSlot* meshSlot;
  _FC_SubSlot* subSlot;
  char* newSubsetName;
  int* tempSlotIDs;
  int i, j;

  // default return
  if (subsets)
    *subsets = NULL; 

  // check input  
  if (!subsets || !seqSubset || !fc_isSeqSubsetValid(numStep, seqSubset) || numStep == 0) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }


  subSlot = _fc_getSubSlot(seqSubset[0]);
  meshSlot = _fc_getMeshSlot(subSlot->mesh);

  // log message
  fc_printfLogMessage("Copying seq subset '%s' to basic subsets",
                      subSlot->header.name);

  //remove handles from the list of sequence subs
  if (meshSlot->numSeqSub == 1){
    free(meshSlot->numStepPerSeqSub);
    meshSlot->numStepPerSeqSub = NULL;
    free(meshSlot->seqSubIDs[0]);
    free(meshSlot->seqSubIDs);
    meshSlot->seqSubIDs = NULL;
    meshSlot->numSeqSub = 0;
  } else {
    for (i = 0; i < meshSlot->numSeqSub; i++) {
      if ((meshSlot->seqSubIDs)[i][0] == subSlot->header.slotID)
	break;
    }
    free((meshSlot->seqSubIDs)[i]);
    for (j = i+1; j < meshSlot->numSeqSub; j++){
      (meshSlot->numStepPerSeqSub)[j-1] = (meshSlot->numStepPerSeqSub)[j];
      (meshSlot->seqSubIDs)[j-1] = (meshSlot->seqSubIDs)[j];
    }
    meshSlot->numSeqSub--;
  }

  tempSlotIDs = realloc(meshSlot->subIDs, sizeof(int)*(meshSlot->numSub+numStep));
  if (tempSlotIDs == NULL){
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      return FC_MEMORY_ERROR;
  }
  meshSlot->subIDs = tempSlotIDs;

  if (new_name == NULL){
    newSubsetName = (char*)malloc((strlen(subSlot->header.name)+1)*
				  sizeof(char));
    strcpy(newSubsetName, subSlot->header.name);
  } else {
    newSubsetName = (char*)malloc((strlen(new_name)+1)*
				  sizeof(char));
    strcpy(newSubsetName, new_name);
  }

  //add handles to the list of basic subs
  for (i = 0; i < numStep; i++) {
    char* temp_name;
    //change sub contents to reflect new as state non seq subset
    subSlot = _fc_getSubSlot(seqSubset[i]);
    subSlot->sequence = FC_NULL_SEQUENCE;
    subSlot->stepID = -1;
    temp_name = malloc((strlen(newSubsetName)+40)*sizeof(char));
    if (temp_name == NULL) {
      free(newSubsetName);
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      return FC_MEMORY_ERROR;
    }
    snprintf(temp_name, strlen(newSubsetName)+40,
	       "%s_%d", newSubsetName, i);
    _fc_setSlotHeaderName(&subSlot->header, temp_name);
    free(temp_name);
    meshSlot->subIDs[meshSlot->numSub] = subSlot->header.slotID;
    meshSlot->numSub++;
  }
  free(newSubsetName);

  // make return values
  *subsets= (FC_Subset*)malloc(numStep*sizeof(FC_Subset));
  for (i = 0; i < numStep; i++)
    (*subsets)[i] = seqSubset[i];
    
  return FC_SUCCESS;
}


/**
 * \ingroup Subset
 * \brief  Create copy subset associated with a different mesh entity type.
 *
 * \description
 *
 *    This routine is essentially the same as
 *    fc_changeMeshEntityType() but for a subset instead of a list of entity
 *    IDs.
 *
 *    If the new association type is of higher dimensionality than the old
 *    (e.g. new subset will be on element and old was on vertices), the
 *    'doStrict' flag determines whether to keep new entities if not all of
 *    their child entities are in the old subset. 1 = be strict and only have
 *    the entity in the new subset if all of it's children are present in the
 *    old subset. 0 = be lenient and add all entities that have at least one
 *    child in the old subset. 
 *
 *    The new association type cannot be FC_AT_WHOLE_DATASET.
 *
 * \modifications
 *    - 11/29/04 WSK, created.
 */
FC_ReturnCode fc_copySubsetWithNewAssociation(
  FC_Subset subset,      /**< input - handle to subset to be copied */
  char* newSubsetName,   /**< input - new subset's name */
  FC_AssociationType new_assoc,  /**< input - the output association type  */
  int doStrict,          /**< input - flag for dropping incomplete entities */
  FC_Subset *newSubset   /**< output - handle to new subset */
){
  FC_ReturnCode rc; 
  int numMember, *memberIDs, newNumMember, *newMemberIDs;
  FC_Mesh mesh;
  FC_AssociationType assoc;
  _FC_SubSlot *subSlot;

  // default returns
  if (newSubset)
    *newSubset = FC_NULL_SUBSET;        
 
  // check input - NULL name is o.k.
  if (!fc_isSubsetValid(subset) || !fc_isAssociationTypeValid(new_assoc) ||
      new_assoc == FC_AT_UNKNOWN || new_assoc == FC_AT_WHOLE_DATASET ||
      (doStrict != 0 && doStrict != 1) || !newSubset) { 
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // check for null name
  subSlot = _fc_getSubSlot(subset);
  if (newSubsetName == NULL)
    newSubsetName = subSlot->header.name;

  // log message
  fc_printfLogMessage("Copying subset '%s' while changing association type",
                      subSlot->header.name);

  // get other info from source subset & dest mesh
  rc = fc_getSubsetAssociationType(subset, &assoc);
  if (rc != FC_SUCCESS) {
    fc_printfErrorMessage("Failed to get info from subset.");
    return FC_ERROR;
  }
  rc = fc_getMeshFromSubset(subset, &mesh);
  if (rc != FC_SUCCESS) {
    fc_printfErrorMessage("Failed to get mesh from subset.");
    return FC_ERROR;
  }

  // special case -- if assoc is the same, just call fc_copySubset()
  if (new_assoc == assoc)
    return fc_copySubset(subset, mesh, newSubsetName, newSubset);

  // --- do it
  rc = fc_createSubset(mesh, newSubsetName, new_assoc, newSubset);
  if ( rc != FC_SUCCESS) {
    fc_printfErrorMessage("error creating subset.");
    return FC_ERROR;
  }
  rc = fc_getSubsetMembersAsArray(subset, &numMember, &memberIDs); 
  if (rc != FC_SUCCESS) {
    fc_printfErrorMessage("error in getting members from subset.");
    return FC_ERROR;
  }
  if (numMember > 0) {
    rc = fc_changeMeshEntityType(mesh, assoc, numMember, memberIDs, new_assoc,
                                 doStrict, &newNumMember, &newMemberIDs);
    free(memberIDs);
    if (rc != FC_SUCCESS) {
      fc_printfErrorMessage("error changing member association");
      return FC_ERROR;
    }
    if (newNumMember > 0) {
      rc = fc_addArrayMembersToSubset(*newSubset, newNumMember, newMemberIDs);
      if (rc != FC_SUCCESS) {
        fc_printfErrorMessage("error adding members to a new subset.");
        return FC_ERROR;
      }
      free(newMemberIDs);
    }
  }

  return FC_SUCCESS;
}

/**
 * \ingroup Subset
 * \brief  Create subset containing complement of the subset (NOT operator).
 *
 * \description
 *  
 *      Create new complement subset by performing logical NOT operation 
 *      on a subset. New subset will be on the same mesh as original subset. 
 *
 * \modifications
 *    - 05/21/04 RM, created. 
 */
FC_ReturnCode fc_createSubsetComplement(
  FC_Subset subset,       /**< input - handle to original subset       */
  char* newSubsetName,    /**< input - new subset's name               */
  FC_Subset *newSubset    /**< output - handle to new subset           */
){
  FC_ReturnCode rc, rc1, rc2;
  int i;
  int maskLength, *mask;
  _FC_SubSlot *subSlot;
  FC_Mesh mesh;
  FC_AssociationType assoc;

  // default returns
  if (newSubset)
    *newSubset = FC_NULL_SUBSET;        
 
  // check input
  if (!fc_isSubsetValid(subset) || !newSubsetName || !newSubset) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  subSlot = _fc_getSubSlot(subset);
  fc_printfLogMessage("Creating complement to subset '%s' with name '%s'",
                      subSlot->header.name, newSubsetName);
  
  // get subset members
  rc = fc_getMeshFromSubset(subset, &mesh);
  rc1 = fc_getSubsetInfo(subset, NULL, NULL, &assoc);
  rc2 = fc_getSubsetMembersAsMask(subset, &maskLength, &mask); 
  if (rc != FC_SUCCESS || rc1 != FC_SUCCESS || rc2 != FC_SUCCESS) {
    fc_printfErrorMessage("error in accessing subset's fields");
    return FC_ERROR;
  } 
  
  // revert mask in place
  for (i = 0; i < maskLength; i++) 
    mask[i] = ( mask[i] == 0 ) ? 1 : 0;

  // create the new subset
  rc = fc_createSubset(mesh, newSubsetName, assoc, &(*newSubset));
  if (rc1 != FC_SUCCESS) {
    free(mask);
    fc_printfErrorMessage("error creating subset");
    return FC_ERROR;
  }
  rc = fc_addMaskMembersToSubset(*newSubset, maskLength, mask);
  free(mask);

  return rc;
}

/**
 * \ingroup Subset
 * \brief  Create subset based on intersection of two input subsets.
 *
 * \description
 *  
 *      Create new subset by performing logical operation (intersection)
 *      on two input subsets. The operation is given by a string and
 *      can have the following values: 
 *       - \b "AND", \b "and", \b "&&" or \b "&" ---  the logical AND operation
 *       - \b "OR", \b "or", \b "||" or \b "|"   ---  the logical OR op.
 *       - \b "XOR" \b "xor", or \b "^"        ---  the logical XOR  (OR-AND) op.
 *
 *      Subset's need to have the same association type.
 *
 *      FC_NULL_SUBSET is allowed for one of the subsets and behaves the
 *      same as if you had provided an empty subset. 
 *
 * \todo ?Right now can only intersect sets on same mesh.
 *
 * \modifications
 *    - 05/18/04 RM, created.
 *    - 09/21/04 WSK, changed input checking & cleanup.
 */
FC_ReturnCode fc_createSubsetIntersection(
  FC_Subset subset1,    /**< input - first subset's handle               */
  char* operation,      /**< input - string describing logical operation */
  FC_Subset subset2,    /**< input - second subset's handle              */
  char* newSubsetName,  /**< output - new subset's name                  */
  FC_Subset *newSubset  /**< output - new subset handle                  */
) {
  int ioper, i, j, k;
  FC_ReturnCode rc,rc1,rc2;
  FC_Subset tempSubset;
  FC_Mesh mesh1, mesh2;
  FC_AssociationType assoc1, assoc2;
  int numMember1, *members1, numMember2, *members2;
  int maxNumMember, *newMembers, newNumMember;

  // default returns
  if (newSubset)
    *newSubset = FC_NULL_SUBSET;

  // check input (1 NULL subset is o.k., but not both)
  if (!newSubset || !operation || !newSubsetName || !newSubset ||
      (!fc_isSubsetValid(subset1) && !FC_HANDLE_EQUIV(subset1,FC_NULL_SUBSET)) ||
      (!fc_isSubsetValid(subset2) && !FC_HANDLE_EQUIV(subset2,FC_NULL_SUBSET)) ||
      (!fc_isSubsetValid(subset1) && !fc_isSubsetValid(subset2))   ) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // check the ops
  if (!strcmp(operation,"AND") || 
      !strcmp(operation,"and") || 
      !strcmp(operation,"&&")  || 
      !strcmp(operation,"&")      ) 
    ioper = 0;
  else if (!strcmp(operation,"OR") ||
           !strcmp(operation,"or") ||  
           !strcmp(operation,"||") || 
           !strcmp(operation,"|")      ) 
    ioper = 1;
  else if (!strcmp(operation,"XOR") ||
           !strcmp(operation,"xor") ||  
           !strcmp(operation,"^")      )
    ioper = 2;
  else {
        fc_printfErrorMessage("unrecognized operation '%s'",operation);
    return FC_INPUT_ERROR;
  }

  // ----- special case: one subset is FC_NULL_SUBSET 

  if (FC_HANDLE_EQUIV(subset1, FC_NULL_SUBSET) || 
      FC_HANDLE_EQUIV(subset2, FC_NULL_SUBSET)    ) {
    // identify non-null subset & it's parent mesh
    if (!FC_HANDLE_EQUIV(subset1, FC_NULL_SUBSET))
      tempSubset = subset1;
    else
      tempSubset = subset2;
    rc = fc_getMeshFromSubset(tempSubset, &mesh1); 
    if (rc != FC_SUCCESS) {
      fc_printfErrorMessage("Can't get mesh from subset.");
      return  FC_ERROR;
    }

    // if operation is "AND" return an empty subset of appropriate assoc 
    if (ioper == 0 ) {
      rc = fc_getSubsetInfo(tempSubset, NULL, NULL, &assoc1);
      if (rc != FC_SUCCESS) 
        return FC_ERROR;
      rc = fc_createSubset(mesh1, newSubsetName, assoc1, newSubset);
      if (rc != FC_SUCCESS)
        return FC_ERROR;
    }
    // otherwise return a copy of the other subset  
    else {
      rc = fc_copySubset(tempSubset, mesh1, newSubsetName, &(*newSubset));
      if (rc != FC_SUCCESS)
        return  FC_ERROR;
    }     
    return FC_SUCCESS;  // return after special cases
  }
  
  // check if subsets are defined on the same mesh, if not return input error
  rc1 = fc_getMeshFromSubset(subset1, &mesh1);
  rc2 = fc_getMeshFromSubset(subset2, &mesh2);
  if (rc1 != FC_SUCCESS || rc2 != FC_SUCCESS || !FC_HANDLE_EQUIV(mesh1,mesh2)) {
    fc_printfErrorMessage("FC_INPUT_ERROR: Can't get mesh from subset, or "
                          "subsets defined on different meshes");
    return  FC_INPUT_ERROR;
  } 
  
  // get info about the subsets
  rc1 = fc_getSubsetInfo(subset1, NULL, NULL, &assoc1);
  rc2 = fc_getSubsetInfo(subset2, NULL, NULL, &assoc2);
  
  // check if both subsets have the same association
  if (rc1 != FC_SUCCESS || rc2 != FC_SUCCESS || assoc1 != assoc2 ) {
    fc_printfErrorMessage("FC_INPUT_ERROR: Can't get info from subset, or " 
                          "subsets defined with different association type");
    return  FC_INPUT_ERROR;
  }
  
  // --- get new members

  // setup
  rc = fc_getSubsetMaxNumMember(subset1, &maxNumMember);
  rc1 = fc_getSubsetMembersAsArray(subset1, &numMember1, &members1);
  rc2 = fc_getSubsetMembersAsArray(subset2, &numMember2, &members2);
  if (rc != FC_SUCCESS || rc1 != FC_SUCCESS || rc2 != FC_SUCCESS) {
    free(members1);
    free(members2);
    return FC_ERROR;
  }
  newMembers = (int*)malloc(sizeof(int)*maxNumMember);
  if (!newMembers) {
    free(members1);
    free(members2);
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }

  // Walk both arrays at same time, adding to newMembers when appropriate
  newNumMember = 0;
  i = 0;
  j = 0;
  while (i < numMember1 && j < numMember2) {
    if (members1[i] < members2[j]) { // not the same
      if (ioper == 1) { // "OR"
        newMembers[newNumMember] = members1[i];
        newNumMember++;
      }
      else if (ioper == 2) { // "XOR"
        newMembers[newNumMember] = members1[i];
        newNumMember++;
      }
      i++;
    }
    else if (members1[i] > members2[j]) { // not the same
      if (ioper == 1) { // "OR"
        newMembers[newNumMember] = members2[j];
        newNumMember++;
      }
      else if (ioper == 2) { // "XOR"
        newMembers[newNumMember] = members2[j];
        newNumMember++;
      }
      j++;
    }
    else { // the same
      if (ioper == 0) { // "AND"
        newMembers[newNumMember] = members1[i];
        newNumMember++;
      }
      else if (ioper == 1) { // "OR"
        newMembers[newNumMember] = members1[i];
        newNumMember++;
      }
      i++;
      j++;
    }
  }
  for (k = i; k < numMember1; k++) {
    if (ioper > 0) {
      newMembers[newNumMember] = members1[k];
      newNumMember++;
    }
  }
  for (k = j; k < numMember2; k++) {
    if (ioper > 0) {
      newMembers[newNumMember] = members2[k];
      newNumMember++;
    }
  }
  free(members1);
  free(members2);

  // --- make the new subset
  rc = fc_createSubset(mesh1, newSubsetName, assoc1, newSubset);
  if (rc != FC_SUCCESS) {
    free(newMembers);
    return FC_ERROR;
  }

  // --- done if there are no members
  if (newNumMember == 0) {
    free(newMembers);
    return FC_SUCCESS;
  }
  
  rc = fc_addArrayMembersToSubset(*newSubset, newNumMember, newMembers);
  free(newMembers); 
  return rc;
}

//@}

/** \name Get subset. */
//-------------------------------
//@{

/**
 * \ingroup  Subset
 * \brief  Get the subsets in a mesh.
 *
 * \description
 *  
 *    This returns an array of the subsets in the mesh. The caller is
 *    responsible for freeing the array of mesh handles.
 *
 * \modifications  
 *   - 05/13/04 WSK Created.  
 *   - 9/9/04 WSK, changed to return subset handles instead of subset
 *               names.
  */
FC_ReturnCode fc_getSubsets(
 FC_Mesh mesh,        /**< input - mesh */
 int *numSubset,      /**< output - number of subsets */
 FC_Subset **subsets  /**< output - array of subsets */
) {
  int i;
  _FC_MeshSlot* meshSlot;
  _FC_SubSlot* subSlot;
  
  // set up defaults
  if (numSubset != NULL)
    *numSubset = -1;
  if (subsets != NULL)
    *subsets = NULL;
  
  // check input
  meshSlot = _fc_getMeshSlot(mesh);
  if (meshSlot == NULL || numSubset == NULL || subsets == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Getting subsets from mesh '%s'", 
                      meshSlot->header.name);

  // get num subset
  *numSubset = meshSlot->numSub;
  
  // get subsets
  if (meshSlot->numSub > 0) {
    *subsets = malloc(sizeof(FC_Subset) * meshSlot->numSub);
    if (*subsets == NULL) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      return FC_MEMORY_ERROR;
    }
    for (i = 0; i < meshSlot->numSub; i++) {
      subSlot = _fc_getSubSlotFromID(meshSlot->subIDs[i]);
      _FC_GET_HANDLE((*subsets)[i], subSlot);
    }
  }
  
  return FC_SUCCESS;
}


/**
 * \ingroup  Subset
 * \brief  Get the seq subsets in a mesh.
 *
 * \description
 *  
 *    This returns an array of the seq subsets in the mesh (array of arrays of subsets)
 *    as well as an array of the number of steps per seq subset. The caller is
 *    responsible for freeing all returned arrays. (SeqSubs is a 2D array 
 *    so to free: for (i = 0; i < numSeqSub) { free(seqSub[i])}  )
 *
 * \modifications  
 *   - 03/03/08 ACG Created.  
 */
FC_ReturnCode fc_getSeqSubsets(
 FC_Mesh mesh,        /**< input - mesh */
 int *numSeqSubs,      /**< output - number of seq subsets */
 int **numStepPerSub,      /**< output - array with number of steps for
				each seq sub */
 FC_Subset ***seqSubs  /**< output - array of seq subsets (each
			   seq sub is an array of numStep subsets) */
) {
  _FC_MeshSlot* meshSlot;
  _FC_SubSlot* subSlot;
  int i,j;
  
  // set up defaults
  if (numSeqSubs != NULL)
    *numSeqSubs = -1;
  if (numStepPerSub != NULL)
    *numStepPerSub = NULL;
  if (seqSubs != NULL)
    *seqSubs = NULL;

  // check input
  meshSlot = _fc_getMeshSlot(mesh);
  if (meshSlot == NULL || numSeqSubs == NULL || numStepPerSub == NULL ||
      seqSubs == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Getting seq subsets from mesh '%s'", 
                      meshSlot->header.name);

  // get num seq subsets
  *numSeqSubs = meshSlot->numSeqSub;
  
  // get seq subsets
  if (*numSeqSubs > 0) {
    *numStepPerSub = malloc(sizeof(int) * (*numSeqSubs));
    if (*numStepPerSub == NULL) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      return FC_MEMORY_ERROR;
    }
    *seqSubs = malloc(sizeof(FC_Subset*) * (*numSeqSubs));
    if (*seqSubs == NULL) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      return FC_MEMORY_ERROR;
    }
    for (i = 0; i < *numSeqSubs; i++) {
      (*numStepPerSub)[i] = meshSlot->numStepPerSeqSub[i];
      (*seqSubs)[i] = malloc(meshSlot->numStepPerSeqSub[i]*sizeof(FC_Subset));
      if ( (*seqSubs)[i] == NULL) {
	fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
	return FC_MEMORY_ERROR;
      }
      for (j = 0; j < meshSlot->numStepPerSeqSub[i]; j++){
	subSlot = _fc_getSubSlotFromID(meshSlot->seqSubIDs[i][j]);
	_FC_GET_HANDLE((*seqSubs)[i][j],subSlot);
      }
    }
  }
  
  return FC_SUCCESS;
}


/**
 * \ingroup  Subset
 * \brief  Get the number of subsets in a mesh.
 *
 * \modifications  
 *    - 9/28/04 WSK created.
 */
FC_ReturnCode fc_getNumSubset(
 FC_Mesh mesh,        /**< input - mesh */
 int *numSubset       /**< output - number of subsets */
) {
  _FC_MeshSlot* meshSlot;
  
  // set up defaults
  if (numSubset != NULL)
    *numSubset = -1;
  
  // check input
  meshSlot = _fc_getMeshSlot(mesh);
  if (meshSlot == NULL || numSubset == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Getting number of subsets from mesh '%s'", 
                      meshSlot->header.name);

  // get num subset
  *numSubset = meshSlot->numSub;

  return FC_SUCCESS;
}


/**
 * \ingroup  Subset
 * \brief  Get the number of seq subsets in a mesh.
 *
 * \modifications  
 *    - 03/03/08 ACG created
 */
FC_ReturnCode fc_getNumSeqSubset(
 FC_Mesh mesh,        /**< input - mesh */
 int *numSeqSubset       /**< output - number of subsets */
) {
  _FC_MeshSlot* meshSlot;
  
  // set up defaults
  if (numSeqSubset != NULL)
    *numSeqSubset = -1;
  
  // check input
  meshSlot = _fc_getMeshSlot(mesh);
  if (meshSlot == NULL || numSeqSubset == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Getting number of seq subsets from mesh '%s'", 
                      meshSlot->header.name);

  // get num subset
  *numSeqSubset = meshSlot->numSeqSub;

  return FC_SUCCESS;
}

/**
 * \ingroup  Subset
 * \brief  Get subsets of a given name from a mesh.
 *
 * \description
 *  
 *    Returns subset on the mesh with the specified name.
 *
 * \modifications  
 *   - 05/27/04 WSK Created.
 *   - 9/9/03 WSK, changed name and moved to mesh (to be next to
 *       fc_getSubsets) 
 *   - 9/26/06 ACG returns array
 */
FC_ReturnCode fc_getSubsetByName(
  FC_Mesh mesh,           /**< input - mesh handle */
  char *subName,          /**< input - subset name */
  int *numSubsets,        /**< output - num matching */
  FC_Subset **subsets   /**< output - subset handles */
) {
  int i;
  _FC_MeshSlot* meshSlot;
  _FC_SubSlot *subSlot;
  FC_Subset subset;

  FC_Subset *lmatch = NULL;
  int numlmatch, maxnumlmatch;

  // default values returned if anything goes wrong
  if (numSubsets)
    *numSubsets = -1;
  if (subsets)
    *subsets = NULL;

  // check input
  meshSlot = _fc_getMeshSlot(mesh);
  if (meshSlot == NULL || subName == NULL || !numSubsets || !subsets){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Getting subset '%s'", subName);

  numlmatch = 0;
  maxnumlmatch = 0;

  
  // find the existing slot
  for (i = 0; i < meshSlot->numSub; i++) {
    subSlot = _fc_getSubSlotFromID(meshSlot->subIDs[i]); 
    if (!strcmp(subSlot->header.name, subName)) {
      _FC_GET_HANDLE(subset, subSlot);

      if (lmatch == NULL){
	lmatch = (FC_Subset*)malloc(sizeof(FC_Subset));
	if (lmatch == NULL){
	  fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
	  return FC_MEMORY_ERROR;
	}
	numlmatch = 1;
	maxnumlmatch = 1;
	lmatch[0] = subset;
      }else{
	if(numlmatch == maxnumlmatch){
	  FC_Subset *temp;
	  temp = (FC_Subset*)realloc(lmatch,(2*numlmatch)*sizeof(FC_Subset));
	  if (temp == NULL){
	    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
	    return FC_MEMORY_ERROR;
	  }
	  maxnumlmatch*=2;
	  lmatch = temp;
	}
	lmatch[numlmatch] = subset;
	numlmatch++;
      }
    }
  }

  *numSubsets = numlmatch;
  *subsets = lmatch;

  return FC_SUCCESS;

}


/**
 * \ingroup  Subset
 * \brief  Get sequence subsets of a given name from the mesh.
 *
 * \description
 *
 *    Returns an array of seq subset (array of arrays of subsets)
 *    as well as an array of the number of steps per sequence subset.
 *    The use is responsible for freeing all returned arrays.
 *    (SeqSubsets is a 2D array, so to free: for (i = 0; i < numSeq subs) {
 *    free(seqSubsets[i])); } )
 *
 *
 * \modifications  
 *   - 03/04/08 ACG created
 */
FC_ReturnCode fc_getSeqSubsetByName(
  FC_Mesh mesh,              /**< input - mesh handle */
  char *subName,             /**< input - seq subset name */
  int *numSeqSubset,       /**< output - number of sequence subsets */
  int** numStepPerSubset,  /**< output - array with number of steps for each
                          sequence subset */
  FC_Subset ***seqSubsets  /**< output - array of sequence subsets (each
				  sequence subset is array of numStep
				  subsets */
) {
  int i, j;
  _FC_MeshSlot *meshSlot;
  _FC_SubSlot *subSlot;

  int numStep;
  FC_Subset *seqsub;
  FC_Subset **lmatch = NULL;
  int *lstepmatch = NULL;
  int numlmatch, maxnumlmatch;
  
  // default return values
 if (numSeqSubset != NULL)
    *numSeqSubset = -1;
  if (numStepPerSubset != NULL)
    *numStepPerSubset = NULL;
  if (seqSubsets != NULL)
    *seqSubsets = NULL;

  // check input
  meshSlot = _fc_getMeshSlot(mesh);
  if (meshSlot == NULL || !subName || numSeqSubset == NULL || 
      numStepPerSubset == NULL || seqSubsets == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Getting sequence subset '%s'", subName);

  numlmatch = 0;
  maxnumlmatch = 0;

 // look for subName in this mesh's seq sub list
  for (i = 0; i < meshSlot->numSeqSub; i++) {
    subSlot = _fc_getSubSlotFromID(meshSlot->seqSubIDs[i][0]);
    if (!strcmp(subName, subSlot->header.name)) {
      numStep = meshSlot->numStepPerSeqSub[i];
      seqsub = malloc((numStep)*sizeof(FC_Subset));
      if (seqsub == NULL) {
	fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
	return FC_MEMORY_ERROR;
      }
      for (j = 0; j < numStep; j++) {
	subSlot = _fc_getSubSlotFromID(meshSlot->seqSubIDs[i][j]);
	_FC_GET_HANDLE(seqsub[j], subSlot);
      }

      if (lmatch == NULL){
	lmatch = (FC_Subset**)malloc(sizeof(FC_Subset*));
	lstepmatch = (int*)malloc(sizeof(int));
	if (lmatch == NULL || lstepmatch == NULL){
	  fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
	  return FC_MEMORY_ERROR;
	}
	numlmatch = 1;
	maxnumlmatch = 1;
	lmatch[0] = seqsub;
	lstepmatch[0] = numStep;
      }else{
	if(numlmatch == maxnumlmatch){
	  FC_Subset **temp;
	  int *tempstep;
	  temp = (FC_Subset**)realloc(lmatch,(2*numlmatch)*sizeof(FC_Subset*));
	  tempstep = (int*)realloc(lstepmatch,(2*numlmatch)*sizeof(FC_Subset*));
	  if (temp == NULL || tempstep == NULL){
	    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
	    return FC_MEMORY_ERROR;
	  }
	  maxnumlmatch*=2;
	  lmatch = temp;
	  lstepmatch = tempstep;
	}
	lmatch[numlmatch] = seqsub;
	lstepmatch[numlmatch] = numStep;
	numlmatch++;
      }
    }
  }

  if (numlmatch){
    *numSeqSubset = numlmatch;
    *numStepPerSubset = lstepmatch;
    *seqSubsets = lmatch;
    return FC_SUCCESS;
  }

  //no matches
  *numSeqSubset = 0;
  *numStepPerSubset = NULL;
  *seqSubsets = NULL;
  return FC_SUCCESS;

}


//@}

/** \name Subset comparisons. */
//-------------------------------
//@{

/**
 * \ingroup Subset
 * \brief Test if two subsets intersect.
 *
 * \description
 *  
 *      The indicator flag is set to true (1) if two subsets intersect and
 *      false (0) if they do not, and -1 in case of error.
 *
 *      Either input subset can be FC_NULL_SUBSET, but the routine will return
 *      early with an error if both are. Either subset can be empty, but the
 *      routine will return early with an error if both are.
 *
 *      Subsets must have the same associatation type and be defined on the
 *      same mesh.  These are the same constraints as in \ref
 *      fc_createSubsetIntersection().
 *
 * \modifications
 *    - 10/27/05 ACG created
 *    - 01/16/06 ACG making more efficient by copying guts of
 *            getSubsetIntersection but stopping when a match is found.
 */
FC_ReturnCode fc_doSubsetsIntersect(
  FC_Subset subset1,   /**< input - subset */
  FC_Subset subset2,   /**< input - subset */
  int* indicator        /** output - indicator flag, 1 = they intersect */
 ){
  FC_ReturnCode rc1,rc2;
  FC_Mesh mesh1, mesh2;
  FC_AssociationType assoc1, assoc2;
  int numMember1, *members1, numMember2, *members2;
  int i,j;

  if (indicator)
    *indicator = -1;

  // check input (1 NULL subset is o.k., but not both)
  if (!indicator || 
      (!fc_isSubsetValid(subset1) && !FC_HANDLE_EQUIV(subset1,FC_NULL_SUBSET)) ||
      (!fc_isSubsetValid(subset2) && !FC_HANDLE_EQUIV(subset2,FC_NULL_SUBSET)) ||
      (!fc_isSubsetValid(subset1) && !fc_isSubsetValid(subset2))   ) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  fc_printfLogMessage("doSubsetsIntersect ?");
  
  //if one is NULL then return 0
  if (FC_HANDLE_EQUIV(subset1,FC_NULL_SUBSET) ||
      FC_HANDLE_EQUIV(subset2,FC_NULL_SUBSET)){
    *indicator = 0;
    return FC_SUCCESS;
  }

  // check if subsets are defined on the same mesh, if not return input error
  rc1 = fc_getMeshFromSubset(subset1, &mesh1);
  rc2 = fc_getMeshFromSubset(subset2, &mesh2);
  if (rc1 != FC_SUCCESS || rc2 != FC_SUCCESS || !FC_HANDLE_EQUIV(mesh1,mesh2)) {
    fc_printfErrorMessage("FC_INPUT_ERROR: Can't get mesh from subset, or "
                          "subsets defined on different meshes");
    return  FC_INPUT_ERROR;
  } 
  
  // get info about the subsets
  rc1 = fc_getSubsetInfo(subset1, NULL, NULL, &assoc1);
  rc2 = fc_getSubsetInfo(subset2, NULL, NULL, &assoc2);
  
  // check if both subsets have the same association
  if (rc1 != FC_SUCCESS || rc2 != FC_SUCCESS || assoc1 != assoc2 ) {
    fc_printfErrorMessage("FC_INPUT_ERROR: Can't get info from subset, or " 
                          "subsets defined with different association type");
    return  FC_INPUT_ERROR;
  }

  // setup
  rc1 = fc_getSubsetMembersAsArray(subset1, &numMember1, &members1);
  rc2 = fc_getSubsetMembersAsArray(subset2, &numMember2, &members2);
  if (rc1 != FC_SUCCESS || rc2 != FC_SUCCESS) {
    free(members1);
    free(members2);
    return FC_ERROR;
  }

  //empty subsets
  if (numMember1 == 0 && numMember2 == 0){
    free(members1);
    free(members2);
    fc_printfErrorMessage("FC_INPUT_ERROR: Can't do two empty subsets");
    return  FC_INPUT_ERROR;
  }

  if (numMember1 == 0 || numMember2 == 0){
    free(members1);
    free(members2);
    *indicator = 0;
    return FC_SUCCESS;
  }


  // Walk both arrays at same time
  i = 0;
  j = 0;
  while (i < numMember1 && j < numMember2) {
    if (members1[i] < members2[j]) { // not the same
      i++;
    }
    else if (members1[i] > members2[j]) { // not the same
      j++;
    }
    else { // the same
      //then they intersect
      free(members1);
      free(members2);
      *indicator = 1;
      return FC_SUCCESS;
    }
  }

  //havent found anything in common
  *indicator = 0;
  free(members1);
  free(members2);

  return FC_SUCCESS;
}

/**
 * \ingroup Subset
 * \brief Test if first subset is a superset of the second subset.
 *
 * \description
 *  
 *      The indicator flag is set to true (1) if the first subset is a superset
 *      of the second and false (0) if, not and -1 in case of error.
 *
 *      Either input subset can be FC_NULL_SUBSET, but the routine will return
 *      early with an error if both are. Either subset can be empty, but the
 *      routine will return early with an error if both are.
 *
 *      Subsets must have the same associatation type and be defined on the
 *      same mesh.  These are the same constraints as in \ref
 *      fc_createSubsetIntersection().
 *
 * \modifications
 *    - 02/13/06 ACG created
 *    - 09/28/06 ACG bug fix in walking the array and end condition
 */
FC_ReturnCode fc_isSubsetSuperset(
  FC_Subset subset_super,   /**< input - subset */
  FC_Subset subset_sub,   /**< input - subset */
  int* indicator        /** output - indicator flag, 1 = is a superset */
 ){
  FC_ReturnCode rc1,rc2;
  FC_Mesh mesh1, mesh2;
  FC_AssociationType assoc1, assoc2;
  int numMember1, *members1, numMember2, *members2;
  int i,j;

  if (indicator)
    *indicator = -1;

  // check input (1 NULL subset is o.k., but not both)
  if (!indicator || 
      (!fc_isSubsetValid(subset_super) && !FC_HANDLE_EQUIV(subset_super,FC_NULL_SUBSET)) ||
      (!fc_isSubsetValid(subset_sub) && !FC_HANDLE_EQUIV(subset_sub,FC_NULL_SUBSET)) ||
      (!fc_isSubsetValid(subset_super) && !fc_isSubsetValid(subset_sub))   ) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  fc_printfLogMessage("is Subset Superset ?");
  
  //if one is NULL then return 0
  if (FC_HANDLE_EQUIV(subset_super,FC_NULL_SUBSET) ||
      FC_HANDLE_EQUIV(subset_sub,FC_NULL_SUBSET)){
    *indicator = 0;
    return FC_SUCCESS;
  }

  // check if subsets are defined on the same mesh, if not return input error
  rc1 = fc_getMeshFromSubset(subset_super, &mesh1);
  rc2 = fc_getMeshFromSubset(subset_sub, &mesh2);
  if (rc1 != FC_SUCCESS || rc2 != FC_SUCCESS || !FC_HANDLE_EQUIV(mesh1,mesh2)) {
    fc_printfErrorMessage("FC_INPUT_ERROR: Can't get mesh from subset, or "
                          "subsets defined on different meshes");
    return  FC_INPUT_ERROR;
  } 
  
  // get info about the subsets
  rc1 = fc_getSubsetInfo(subset_super, NULL, NULL, &assoc1);
  rc2 = fc_getSubsetInfo(subset_sub, NULL, NULL, &assoc2);
  
  // check if both subsets have the same association
  if (rc1 != FC_SUCCESS || rc2 != FC_SUCCESS || assoc1 != assoc2 ) {
    fc_printfErrorMessage("FC_INPUT_ERROR: Can't get info from subset, or " 
                          "subsets defined with different association type");
    return  FC_INPUT_ERROR;
  }

  // setup
  rc1 = fc_getSubsetMembersAsArray(subset_super, &numMember1, &members1);
  rc2 = fc_getSubsetMembersAsArray(subset_sub, &numMember2, &members2);
  if (rc1 != FC_SUCCESS || rc2 != FC_SUCCESS) {
    free(members1);
    free(members2);
    return FC_ERROR;
  }

  //empty
  if (numMember1 == 0 && numMember2 == 0){
    free(members1);
    free(members2);
    fc_printfErrorMessage("FC_INPUT_ERROR: Can't do two empty subsets");
    return  FC_INPUT_ERROR;
  }

  if (numMember1 == 0 || numMember2 == 0){
    free(members1);
    free(members2);
    *indicator = 0;
    return FC_SUCCESS;
  }

  if (numMember1 < numMember2){
    free(members1);
    free(members2);
    *indicator = 0;
    return FC_SUCCESS;
  }


  // Walk both arrays at same time
  i = 0;
  j = 0;
  while (i < numMember1 && j < numMember2){
    if (members1[i] > members2[j]){
      free(members1);
      free(members2);
      *indicator = 0;
      return FC_SUCCESS;
    }else if (members1[i] == members2[j]){
      i++;
      j++;
    }else{
      i++;
    }
  }

  free(members1);
  free(members2);
  //did we get the last one ?
  *indicator = (j == numMember2);
  return FC_SUCCESS;
}

//@}

/** \name Change the name of a subset. */
//-------------------------------
//@{

/**
 * \ingroup Subset
 * \brief  Change the name of a subset.
 *
 * \modifications
 *   - 8/24/2005 WSD. Created.
 */
FC_ReturnCode fc_changeSubsetName(
  FC_Subset subset, /**< Input - the subset. */
  char* newName       /**< Input - new name for the subset. */
) {
  _FC_SubSlot* subSlot;

  // check input
  subSlot = _fc_getSubSlot(subset);
  if (subSlot == NULL || newName == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // log message
  fc_printfLogMessage("Changing name of subset from '%s' to '%s'", 
                      subSlot->header.name, newName);
  
  // Do it
  return _fc_setSlotHeaderName(&subSlot->header, newName);
}

/**
 * \ingroup Subset
 * \brief  Change the name of a sequence subset
 *
 * \modifications
 *   - 03/04/08 ACG created
 */
FC_ReturnCode fc_changeSeqSubsetName(
  int numStep, /**< Input - Number of steps in the seq subset. */
  FC_Subset* subsets, /**< Input - the steps of the seq subset. */
  char* newName       /**< Input - new name for the seq subset. */
) {
  FC_ReturnCode rc, rc_keep = FC_SUCCESS;
  int i;
  _FC_SubSlot* subSlot;

  // check input
  if (!fc_isSeqSubsetValid(numStep, subsets) || newName == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // log message
  subSlot = _fc_getSubSlot(subsets[0]); 
  fc_printfLogMessage("Changing name of seq subset from '%s' to '%s'", 
                      subSlot->header.name, newName);
  
  // Do it
  for (i = 0; i < numStep; i++) {
    rc = _fc_setSlotHeaderName(&subSlot->header, newName);
    if (rc != FC_SUCCESS)
      rc_keep = rc;
  }
  return rc_keep;
}

//@}

/** \name Add members to subset */
//-------------------------------
//@{

/**
 * \ingroup Subset
 * \brief  Add member to a subset. 
 *
 * \description  
 *
 *    Add member to a given subset.  Returns with an error if the member is
 *    outside of the valid range.
 *
 *    Members cannot be added to a committed subset (i.e. exists on disk).
 *
 * \modifications  
 *    - 04/30/04 RM, created. 
 *    - 05/11/05 WSK, changed to reflect changes to FC_Subset.
 */
FC_ReturnCode fc_addMemberToSubset(
  FC_Subset subset,    /**< input - subset to add member to  */
  int memberID         /**< input - ID of member to add    */
) {
  _FC_SubSlot* subSlot;
  int ret;

  // check input, subset->storageType should never be FC_ST_UNKNOWN,
  // memberID has to be positive >= 0  
  subSlot = _fc_getSubSlot(subset);
  if(!subSlot || memberID < 0 || memberID > subSlot->maxNumMember-1 ) { 
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // can't modify a committed subset
  if (subSlot->header.committed) {
    fc_printfErrorMessage("Can't modify committed subset '%s'",
                          subSlot->header.name);
    return FC_INPUT_ERROR;
  }

  // no logging on purpose since this may get called over and over

  // add memberID to the subset
  ret = fc_addIntToSortedIntArray(&subSlot->sia, memberID);
  if (ret < 0) {
    fc_printfErrorMessage("Failed to add int to sorted int array");
    return ret; // an error occured
  }
  subSlot->numMember += ret;
 
  return FC_SUCCESS;
}

/**
 * \ingroup Subset
 * \brief  Add members from an array to a subset. 
 *
 * \description  
 *
 *    Add members given in an array to a subset. 
 *    The operation is done on native subset representation type. 
 *    \ref FC_INPUT_ERROR is returned if any of the input members are 
 *    outside the valid range.
 *
 *    An empty array is o.k and will not return an error.
 *
 *    Members cannot be added to a committed subset (i.e. written to disk).
 *
 * \modifications  
 *    - 04/30/04 RM, created. 
 *    - 05/11/05 WSK, changed to reflect changes to FC_Subset.
 *    - 06/02/06 WSD, added array specific imlementation instead of calling
 *                    fc_addMembersToSubset() which was not efficient enough.
 */
FC_ReturnCode fc_addArrayMembersToSubset(
  FC_Subset subset,/**< input/output - subset handle                       */
  int numMember,   /**< input - number of members in an input array       */
  int *array       /**< input - pointer to array of member IDs */
) {
  int i;
  int ret;
  _FC_SubSlot* subSlot;

  // check input
  subSlot = _fc_getSubSlot(subset);
  if( !subSlot || numMember < 0 || ( numMember > 0 && !array)) { 
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }     
  for (i = 0; i < numMember; i++) {
    if (array[i] < 0 || array[i] >= subSlot->maxNumMember) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
      return FC_INPUT_ERROR;
    }
  }

   // can't modify a committed subset
  if (subSlot->header.committed) {
    fc_printfErrorMessage("Can't modify committed subset '%s'",
                          subSlot->header.name);
    return FC_INPUT_ERROR;
  }

 // no logging on purpose since this may get called over and over

  // trivial case, no members in an array
  if (numMember == 0) 
    return FC_SUCCESS;

  // do it
  ret = fc_addIntArrayToSortedIntArray(&subSlot->sia, numMember, array, 0);
  if (ret < 0) {
    fc_printfErrorMessage("Failed to add int array to sorted int array");
    return ret;
  }
  subSlot->numMember += ret;

  return FC_SUCCESS;
}

/**
 * \ingroup Subset
 * \brief  Add members from a mask to a subset. 
 *
 * \description
 *  
 *    Add members given by a mask to a subset.  The operation is done on native
 *    subset representation type.  FC_INPUT_ERROR is returned if the given mask
 *    has different size that the mask in subset.
 *
 *    A mask of zero length is o.k. and will not return an error.
 *
 *    Members cannot be added to a committed subset (i.e. written to disk).
 *
 * \modifications  
 *    - 04/30/04 RM, created. 
 *    - 05/11/05 WSK, changed to reflect changes to FC_Subset.
 */
FC_ReturnCode fc_addMaskMembersToSubset(
  FC_Subset subset,  /**< input - subset handle */
  int maskLength,  /**< input - number of elements in an input mask */
  int *mask          /**< input - pointer to input mask with subset members */
) {
  int i;
  FC_ReturnCode rc;  
  _FC_SubSlot* subSlot;

  // check input
  subSlot = _fc_getSubSlot(subset);
  if ( !subSlot || maskLength < 0 || (maskLength > 0 && !mask) ||
       (maskLength > 0  && maskLength != subSlot->maxNumMember)) { 
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  } 

  // no logging on purpose since this may get called over and over

  // can't modify a committed subset
  if (subSlot->header.committed) {
    fc_printfErrorMessage("Can't modify committed subset '%s'",
                          subSlot->header.name);
    return FC_INPUT_ERROR;
  }

  // trivial case, no elements in a mask
  if (maskLength == 0) 
    return FC_SUCCESS;

  for (i = 0; i < maskLength; i++) {
    // add member where mask == 1
    if (mask[i] == 1) {
      rc = fc_addMemberToSubset(subset, i);
      if (rc != FC_SUCCESS ) {
        fc_printfErrorMessage("Failed to add member to subset '%s'",
                              subSlot->header.name);
        return FC_ERROR;
      }
    }
  }

  return FC_SUCCESS;
}

//@}

/** \name Delete members from a subset. */
//-------------------------------
//@{

/**
 * \ingroup Subset
 * \brief  Delete member from a subset. 
 *
 * \description 
 *  
 *    Delete member from a given subset, does not change subset type--operation
 *    is done on the native subset type.  \ref FC_INPUT_ERROR is returned if
 *    the member is outside of the valid range.
 *
 *    Members cannot be deleted from a committed subset (i.e. written to disk).
 *
 * \modifications  
 *    - 05/03/04 RM, created. 
 *    - 05/11/04 WSK, changed to reflect changes to FC_Subset.
 */
FC_ReturnCode fc_deleteMemberFromSubset(
  FC_Subset subset,  /**< input/output - subset handle   */
  int member         /**< input - member to delete       */
) {
  _FC_SubSlot* subSlot;
  int ret;

  // check input, subset->storageType should never be FC_ST_UNKNOWN,
  // member has to be positive >= 0  
  subSlot = _fc_getSubSlot(subset);
  if (!subSlot || member < 0 || member > subSlot->maxNumMember-1 ) { 
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // can't modify a committed subset
  if (subSlot->header.committed) {
    fc_printfErrorMessage("Can't modify committed subset '%s'",
                          subSlot->header.name);
    return FC_INPUT_ERROR;
  }

  // no logging on purpose since this may get called over and over

  // delete member from the subset
  ret = fc_deleteIntFromSortedIntArray(&subSlot->sia, member);
  if (ret < 0) {
    fc_printfErrorMessage("Failed to delete in from sorted int array");
    return ret;
  }
  subSlot->numMember -= ret;
  
  return FC_SUCCESS;
}

/**
 * \ingroup Subset
 * \brief  Delete array members from a subset. 
 *
 * \description
 *  
 *    Delete members given in an array from a subset. 
 *    The operation is done on native subset representation type. 
 *    \ref FC_INPUT_ERROR will be returned if any of the array members 
 *    are outside the valid range.
 *
 *    An empty array is o.k. and will not return an error.
 *
 *    Members cannot be deleted from a committed subset (i.e. written to disk).
 *
 * \modifications  
 *    - 05/03/04 RM, created. 
 *    - 05/11/05 WSK, changed to reflect changes to FC_Subset.
 */
FC_ReturnCode fc_deleteArrayMembersFromSubset(
  FC_Subset subset, /**< input - subset handle */
  int numElement,   /**< input - number of elements in an input array */
  int *array        /**< input - pointer to input array with subset members */
) {
  int i;
  FC_ReturnCode rc;
  _FC_SubSlot* subSlot;

  // check input
  subSlot = _fc_getSubSlot(subset);
  if ( !subSlot || numElement < 0 || (numElement > 0 && !array) ) { 
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }     
  for (i = 0; i < numElement; i++) {
    if (array[i] < 0 || array[i] >= subSlot->maxNumMember) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
      return FC_INPUT_ERROR;
    }
  }

  // can't modify a committed subset
  if (subSlot->header.committed) {
    fc_printfErrorMessage("Can't modify committed subset '%s'",
                          subSlot->header.name);
    return FC_INPUT_ERROR;
  }

  // trivial case, no members in an array
  if (numElement == 0) 
    return FC_SUCCESS;

  // do it
  for (i = 0; i < numElement; i++) {
    rc = fc_deleteMemberFromSubset(subset, array[i]);
    if(rc != FC_SUCCESS ) {
      fc_printfErrorMessage("Failed to delete member from subset '%s'",
                            subSlot->header.name);
      return FC_ERROR;
    }
  }

  return FC_SUCCESS;
}

/**
 * \ingroup Subset
 * \brief  Delete mask members from a subset. 
 *
 * \description
 *  
 *    Delete members given by a mask from a subset. 
 *    The operation is done on native subset representation type. 
 *    \ref FC_INPUT_ERROR is returned if the given mask
 *    has different size that the mask in subset.
 *
 *    Members cannot be deleted from a committed subset (i.e. written to disk).
 *
 * \modifications  
 *    - 05/03/04 RM, created. 
 *    - 05/11/05 WSK, changed to reflect changes to FC_Subset.
 */
FC_ReturnCode fc_deleteMaskMembersFromSubset(
  FC_Subset subset,  /**< input - subset handle */
  int maskLength,  /**< input - number of elements in an input mask */
  int *mask          /**< input - pointer to input mask with subset members */
) {
  int i;
  FC_ReturnCode rc;  
  _FC_SubSlot* subSlot;

  // check input
  subSlot = _fc_getSubSlot(subset);
  if ( !subSlot || maskLength < 0 || (maskLength > 0 && !mask) || 
       (maskLength > 0 && maskLength != subSlot->maxNumMember) ) { 
   fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }  
   
  // can't modify a committed subset
  if (subSlot->header.committed) {
    fc_printfErrorMessage("Can't modify committed subset '%s'",
                          subSlot->header.name);
    return FC_INPUT_ERROR;
  }

  // no logging on purpose since this may get called over and over

  // trivial case, no elements in a mask
  if (maskLength == 0) 
    return FC_SUCCESS;

  for (i = 0; i < maskLength; i++)
    // delete member where mask == 1
    if (mask[i] == 1) {
      rc = fc_deleteMemberFromSubset(subset, i);
      if (rc != FC_SUCCESS ) {
        fc_printfErrorMessage("Failed to delete member to subset '%s'",
                              subSlot->header.name);
        return FC_ERROR;
      }
    }

  return FC_SUCCESS;
}

//@}

/** \name Release/delete subsets. */
//-------------------------------
//@{

/**
 * \ingroup  Subset
 * \brief  Attempt to minimize subset's memory usage.
 *
 * \description
 *  
 *    Call to try and release un-needed memory in a  subset. 
 *    If the subset has been saved to disk, the large data is 
 *    released. If the subset has not been saved to disk, this 
 *    will do nothing.
 *
 * \modifications  
 *    - 04/19/04 WSK Created.
 */
FC_ReturnCode fc_releaseSubset(
 FC_Subset subset       /**< input - subset handle */
) {
  _FC_SubSlot* subSlot;
  
  // special case: releasing a null handle is not an error
  if (FC_HANDLE_EQUIV(subset, FC_NULL_SUBSET)) {
    fc_printfWarningMessage("Releasing FC_NULL_SUBSET");
    return FC_SUCCESS;
  }

  // check input
  subSlot = _fc_getSubSlot(subset);
  if (subSlot == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Cleaning up subset '%s'", subSlot->header.name);

  // remove big data if this slot is committed
  if (subSlot->header.committed == 1) {
    fc_freeSortedIntArray(&subSlot->sia);
  }

  return FC_SUCCESS;
}


/**
 * \ingroup  Subset
 * \brief  Attempt to minimize seq subset's memory usage.
 *
 * \description
 *  
 *    Call to try and release un-needed memory in a seq subset. 
 *    If the subset has been saved to disk, the large data is 
 *    released. If the subset has not been saved to disk, this 
 *    will do nothing.
 *
 *   fc_releaseSubset() can be called on specific steps of the
 *   seq subset.
 *
 *   \todo
 *    - test this. note that there is no way to have
 *      commmitted seqSubsets in exodus
 *
 * \modifications  
 *    - 03/03/08 ACG Created.
 */
FC_ReturnCode fc_releaseSeqSubset(
  int numStep, /**< input - number of steps in the seq subset */
  FC_Subset* subsets       /**< input - subset handles */
) {
  FC_ReturnCode rc, rc_keep = FC_SUCCESS;
  int i;
  _FC_SubSlot* subSlot;

  // special case: releasing a null handle is not an error
  if ( numStep < 1 || subsets == NULL) {
    fc_printfWarningMessage("Releasing NULL seq subset");
    return FC_SUCCESS;
  }

  // check input
  if (!fc_isSeqSubsetValid(numStep, subsets)) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  subSlot = _fc_getSubSlot(subsets[0]);
  if (subSlot == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Cleaning up subset '%s'", subSlot->header.name);

  //do it
  rc_keep = FC_SUCCESS;
  for (i = 0; i < numStep; i++) {
    rc = fc_releaseSubset(subsets[i]);
    if (rc != FC_SUCCESS)
      rc_keep = rc;
  }

  return rc_keep;

}

/**
 * \ingroup Subset
 * \brief  Delete a subset.
 *
 * \description
 * 
 *    Call when this subset is no longer needed.
 *
 * \modifications  
 *    - 04/20/04 RM, created.
 *    - 05/11/05 WSK, changed to reflect changes to FC_Subset.
 *    - 9/9/04 WSK changed name to delete. Changed behavior so that
 *      delete always deletes no matter whether on disk or not.
 */
FC_ReturnCode fc_deleteSubset(
  FC_Subset subset         /**< input - subset handle */
) {
  int i, j;
  _FC_SubSlot* subSlot;
  _FC_MeshSlot* meshSlot;
         
  // special case: deleting a null handle is not an error
  if (FC_HANDLE_EQUIV(subset, FC_NULL_SUBSET)) {
    fc_printfWarningMessage("Deleting FC_NULL_SUBSET");
    return FC_SUCCESS;
  }

  // check input
  subSlot = _fc_getSubSlot(subset);
  if (subSlot == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
        
  // log message
  fc_printfLogMessage("Deleting subset '%s'", subSlot->header.name);

  // setup
  meshSlot = _fc_getMeshSlot(subSlot->mesh);

  // remove reference from mesh
  if (meshSlot->numSub == 1) {
    free(meshSlot->subIDs);
    meshSlot->subIDs = NULL;
  }
  else {
    for (i = 0; i < meshSlot->numSub; i++) {
      if (meshSlot->subIDs[i] == subSlot->header.slotID)
        break;
    }
    for (j = i+1; j < meshSlot->numSub; j++) {
      meshSlot->subIDs[j-1] = meshSlot->subIDs[j];
    }
  }
  meshSlot->numSub--;

  // clear & delete slot
  return _fc_deleteSubSlot(subset);
}


/**
 * \ingroup Subset
 * \brief  Delete a seq subset.
 *
 * \description
 * 
 *    Call when this seq subset is no longer needed.
 *
 *    NOTE: this does not delete the array of subset handles that makes
 *    up the seq subset, the user is still responsible for this.
 *
 * \modifications  
 *    - 03/03/08 ACG created
 */
FC_ReturnCode fc_deleteSeqSubset(
  int numStep,		/**< input - the number of subset handles*/		 
  FC_Subset* subsets    /**< input - subset handles */
) {
  FC_ReturnCode rc, rc_keep;
  _FC_SubSlot* first_subSlot;
  _FC_MeshSlot* meshSlot;
  int i, j;
         
  // special case: deleting a null list is not an error
  if (numStep <= 0 && subsets == NULL) {
    fc_printfWarningMessage("Deleting NULL seq subset");
    return FC_SUCCESS;
  }

  // check input
  if (!fc_isSeqSubsetValid(numStep, subsets)) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  first_subSlot = _fc_getSubSlot(subsets[0]);
  if (first_subSlot == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  // log message
  fc_printfLogMessage("Deleting subset '%s'", first_subSlot->header.name);

  //remove references from mesh
  meshSlot = _fc_getMeshSlot(first_subSlot->mesh);
  if (meshSlot->numSeqSub == 1 ) {
    free(meshSlot->numStepPerSeqSub);
    meshSlot->numStepPerSeqSub = NULL;
    free(meshSlot->seqSubIDs[0]);
    free(meshSlot->seqSubIDs);
    meshSlot->seqSubIDs = NULL;
    meshSlot->numSeqSub = 0;
  } else {
    for (i = 0; i < meshSlot->numSeqSub; i++){
      if ((meshSlot->seqSubIDs)[i][0] == first_subSlot->header.slotID)
        break;
    }
    free((meshSlot->seqSubIDs)[i]);
    for (j = i+1; j < meshSlot->numSeqSub; j++) {
      (meshSlot->numStepPerSeqSub)[j-1] = (meshSlot->numStepPerSeqSub)[j];
      (meshSlot->seqSubIDs)[j-1] = (meshSlot->seqSubIDs)[j];
    }
    meshSlot->numSeqSub--;
  }

  // clear & delete slots
  rc_keep = FC_SUCCESS;
  for (i = 0; i < numStep; i++){
    rc = _fc_deleteSubSlot(subsets[i]);
    if (rc != FC_SUCCESS)
      rc_keep = rc;
  }

  return rc_keep;
}

//@}

/** \name Get subset metadata. */
//-------------------------------
//@{

/**
 * \ingroup  Subset
 * \brief Check that the handle refers to a valid subset.
 *
 * \description
 *
 *    Returns 1 (true) if the handle is valid, and 0 (false) if it is not.
 *
 * \modifications  
 *    - 05/25/04 WSK Created.
 */
int fc_isSubsetValid(
  FC_Subset subset  /**< input - subset handle */
) {
  return _fc_isHandleValid(subset.slotID, subset.uID, subTableSize,
                           (_FC_SlotHeader**)subTable);
}


/**
 * \ingroup  Subset
 * \brief Check that the handle refers to a valid subset.
 *
 * \description
 *
 *    Returns 1 (true) if the seq subset handle is valid, and 0 (false)
 *    if it is not. Valid means numStep is correct and all subset handles
 *    are valid and in correct order)
 *
 * \modifications  
 *    - 03/03/08 ACG created
 */
int fc_isSeqSubsetValid(
  int numStep,	    /**< input - num of steps in seq sub */
  FC_Subset *seqSub  /**< input - seq subset (array of subsets) */
) {
  int i;
  _FC_SubSlot* subSlot;
  _FC_SeqSlot* seqSlot;
  _FC_MeshSlot* meshSlot;
  int sub_id;  //id's into mesh's slotID arrays

  //quick check
  if (numStep < 0 || seqSub == NULL || !fc_isSubsetValid(seqSub[0]))
    return 0;
  // check that numStep agrees with sequence
  subSlot = _fc_getSubSlot(seqSub[0]);
  seqSlot = _fc_getSeqSlot(subSlot->sequence);
  if (numStep != seqSlot->numStep)
    return 0;
  //check that all steps are valid and fromt he same subset
  //instead of checking metadata, we check against the parent's list
  meshSlot = _fc_getMeshSlot(subSlot->mesh);
  sub_id = -1;
  for (i = 0; i < meshSlot->numSeqSub; i++){
    if (subSlot->header.slotID == (meshSlot->seqSubIDs)[i][0]){
      sub_id = i;
      break;
    }
  }
  if (sub_id == -1){
    fc_printfErrorMessage("bad code? shouldn't ever happen?");
    return FC_ERROR;
  }
  for (i = 1; i < numStep; i++){
    if (!fc_isSubsetValid(seqSub[i]) ||
	seqSub[i].slotID != (meshSlot->seqSubIDs)[sub_id][i])
      return 0;
  }

  //passed all tests
  return 1;
}



/**
 * \ingroup Subset
 * \brief  Get subset name.
 *
 * \description
 *  
 *    Return copy of the name from the subset. User should free
 *    this memory later.
 *
 * \modifications  
 *    - 04/20/04 RM, created. 
 *    - 05/11/05 WSK, changed to reflect changes to FC_Subset.
 */
FC_ReturnCode fc_getSubsetName(
  FC_Subset subset,    /**< input - subset handle */
  char **subsetName    /**< output - subset name  */
) {
  char *name;
  _FC_SubSlot* subSlot;

  // default return
  if(subsetName)
    *subsetName = NULL;
        
  // check input
  subSlot = _fc_getSubSlot(subset);
  if(subSlot == NULL || subsetName == NULL) { 
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Getting name for subset '%s'", subSlot->header.name);

  // copy name
  name = subSlot->header.name;
  *subsetName = malloc(strlen(name) + 1);
  if (*subsetName == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  strcpy(*subsetName, name); 

  return FC_SUCCESS;
} 

/**
 * \ingroup Subset
 * \brief  Get the parent mesh of the subset.
 *
 * \modifications  
 *    - 04/20/04 RM, created. 
 *    - 05/11/05 WSK, changed to reflect changes to FC_Subset.
 */
FC_ReturnCode fc_getMeshFromSubset(
  FC_Subset subset,   /**< input - subset handle   */
  FC_Mesh *mesh       /**< output - mesh handle    */ 
) {
  _FC_SubSlot* subSlot;
  
  // default return
  if (mesh) 
    *mesh = FC_NULL_MESH;

  // check input   
  subSlot = _fc_getSubSlot(subset);
  if(!subSlot || !mesh) { 
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }  

  // log message
  fc_printfLogMessage("Getting mesh from subset '%s'", subSlot->header.name);
  
  // return requested value
  *mesh = subSlot->mesh; 

  return FC_SUCCESS;
} 

/**
 * \ingroup Subset
 * \brief  Get information about a subset
 *
 * \description
 *  
 *    Get number of elements in a subset, subset type and its association. 
 *    Set an output argument to NULL if you don't need it to be returned.
 *
 * \modifications  
 *    - 04/20/04 RM, created. 
 *    - 05/11/04 WSK, changed to reflect changes to FC_Subset.
 *    - 05/17/04 RM, added maxNumMember into return parameters
 *    - 09/21/04 WSK, removed storage type from arguments.
 */
FC_ReturnCode fc_getSubsetInfo(
  FC_Subset subset,          /**< input - subset handle                   */
  int *numMember,    /**< output - optional, number of members in subset */
  int *maxNumMember,         /**< output - optional, maximum value Member */
  FC_AssociationType *assoc  /**< output - optional, association type     */
) {
  _FC_SubSlot* subSlot;

  // default return values
  if (numMember)
    *numMember = -1;
  if (maxNumMember)
    *maxNumMember = -1;
  if (assoc)
    *assoc = FC_AT_UNKNOWN;

  // check input   
  subSlot = _fc_getSubSlot(subset);
  if (subSlot == NULL || (!numMember && !maxNumMember && !assoc)) {
     fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
 
  // log message
  fc_printfLogMessage("Getting info for subset '%s'", subSlot->header.name);

  // fill in requested values
  if (numMember) 
    *numMember = subSlot->numMember;
  if (maxNumMember)
    *maxNumMember = subSlot->maxNumMember;
  if (assoc)
    *assoc = subSlot->assoc;

  return FC_SUCCESS;
        
}       
        
/**
 * \ingroup Subset
 * \brief  Get number of members in a subset.
 *
 * \modifications  
 *    - 04/20/04 RM, created. 
 *    - 05/11/05 WSK, changed to reflect changes to FC_Subset.
 */
FC_ReturnCode fc_getSubsetNumMember(
  FC_Subset subset,   /**< input - subset handle                     */
  int *numMember      /**< output - number of members in the subset  */ 
) {
  _FC_SubSlot* subSlot;

  // default return value
  if (numMember)
    *numMember = -1;

  // check input
  subSlot = _fc_getSubSlot(subset);
  if(subSlot == NULL || !numMember) { 
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
     return FC_INPUT_ERROR;
  }  

  // log message
  fc_printfLogMessage("Getting number of members in the subset '%s'", 
                       subSlot->header.name);

  *numMember = subSlot->numMember; 

  return FC_SUCCESS;
} 

/**
 * \ingroup Subset
 * \brief  Get the maximum possible number of members a subset could have. 
 *
 * \description 
 * 
 *    This number should be equal to the number of mesh subentities
 *    (e.g. vertices, edges, elements, etc). This is the maximum 
 *    possible number of members and is also the length of the subset
 *    mask.
 *
 * \modifications  
 *    - 05/17/04 RM, created.
 */
FC_ReturnCode fc_getSubsetMaxNumMember(
  FC_Subset subset,  /**< input  - subset handle                             */
  int *maxNumMember  /**< output - maximum number the subset member can have */ 
) {
  _FC_SubSlot* subSlot;

  // default return value
  if (maxNumMember)
    *maxNumMember = -1;

  // check input   
  subSlot = _fc_getSubSlot(subset);
  if(subSlot == NULL || !maxNumMember) { 
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }  

  // log message
  fc_printfLogMessage("Getting maximum possible number of members for subset "
                      "'%s'", subSlot->header.name);

  // fill in requested value
  *maxNumMember = subSlot->maxNumMember; 
    
  return FC_SUCCESS;
} 

/**
 * \ingroup Subset
 * \brief  Get subset's association type.
 *
 * \modifications  
 *    - 09/20/04 WSK, created.
 */
FC_ReturnCode fc_getSubsetAssociationType(
  FC_Subset subset,   /**< input - subset handle */
  FC_AssociationType *assoc  /**< output - association type */ 
) {
  _FC_SubSlot *subSlot;

  // default return value
  if (assoc)
    *assoc = FC_AT_UNKNOWN;

  // check input   
  subSlot = _fc_getSubSlot(subset);
  if (subSlot == NULL || !assoc) { 
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }  

  // log message
  fc_printfLogMessage("Getting association type for subset '%s'", 
                      subSlot->header.name);

  // fill in requested value
  *assoc = subSlot->assoc;

  return FC_SUCCESS;
} 

//@}

/** \name Get/query subset big data. */
//-------------------------------
//@{

/**
 * \ingroup Subset
 * \brief  Check to see if a member is in a subset.
 *
 * \description
 *  
 *    Check to see if a member is in a subset. Returns 1 if the member is
 *    present and 0 if not. If an error occurs, return an error code (always
 *    less than 1).
 *
 * \modifications  
 *    - 2006/8/17 WSK, created.
 */
int fc_isMemberInSubset(
  FC_Subset subset, /**< input - subset handle */
  int memberID      /**< input - member ID */
) {
  FC_ReturnCode rc;
  _FC_SubSlot* subSlot;

  // check args
  subSlot = _fc_getSubSlot(subset);
  if (!subSlot) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // easy case 1: out of range
  if (memberID < 0 || memberID >= subSlot->maxNumMember) {
    fc_printfWarningMessage("Asked for member ID (%d) out of range (0-%d)",
			    memberID, subSlot->maxNumMember-1);
    return 0;
  }

  // easy case 2: empty
  if (subSlot->numMember < 1) 
    return 0;

  // lazy member fetching done here
  if (subSlot->header.committed == 1 && subSlot->sia.numVal == 0) {
    rc = _fc_readSubsetMembers(subset);
    if (rc != FC_SUCCESS) {
      fc_printfErrorMessage("Lazy loading of members of subset '%s' failed",
                            subSlot->header.name);
      return rc;
    }
  }

  // general case
  return fc_isIntInSortedIntArray(&subSlot->sia, memberID);
}

/**
 * \ingroup Subset
 * \brief  Get a copy of subset members as array. 
 *
 * \description
 *  
 *    Get a copy of a subset's members as array.  Note 
 *    It's user's responsibility to deallocate array's memory.
 *    The operation is done on native subset representation type.
 *
 * \modifications  
 *    - 04/30/04 RM, created. 
 *    - 05/11/05 WSK, changed to reflect changes to FC_Subset.
 */
FC_ReturnCode fc_getSubsetMembersAsArray (
  FC_Subset subset, /**< input - subset handle */
  int *numMember,  /**< output - number of members in output array */
  int **array     /**< output - pointer to output array with subset members */
) {
  FC_ReturnCode rc;
  int i;
  _FC_SubSlot* subSlot;

  // default returns
  if(numMember)
    *numMember = -1;
  if(array)
    *array = NULL;

  // check input
  subSlot = _fc_getSubSlot(subset);
  if (!subSlot || !numMember || !array) { 
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // FIX? log?

  // trivial case, no members in subset 
  if (subSlot->numMember == 0) {
    *numMember = 0;
    return FC_SUCCESS;
  }

  // lazy member fetching done here
  if (subSlot->header.committed == 1 && subSlot->sia.numVal == 0) {
    rc = _fc_readSubsetMembers(subset);
    if (rc != FC_SUCCESS) {
      fc_printfErrorMessage("Lazy loading of members of subset '%s' failed",
                            subSlot->header.name);
      return rc;
    }
  }

  // allocate memory for resulting array, copy subset members into it
  *array = malloc(sizeof(int)*subSlot->numMember);  
  if(!(*array)) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  //FIX: why doesn't this work?
  //memcpy(*array, subSlot->sia.vals, subSlot->numMember*sizeof(int));
  for (i = 0; i < subSlot->numMember; i++)
    (*array)[i] = subSlot->sia.vals[i];
  *numMember = subSlot->numMember;
 
  return FC_SUCCESS;
}

/**
 * \ingroup Subset
 * \brief  Get a copy of subset members as mask. 
 *
 * \description
 *  
 *    Get a copy of a subset's members as mask array.  Note 
 *    It's user's responsibility to deallocate mask's memory.  
 *
 * \modifications  
 *    - 04/30/04 RM, created. 
 *    - 05/11/05 WSK, changed to reflect changes to FC_Subset.
 */
FC_ReturnCode fc_getSubsetMembersAsMask(
  FC_Subset subset,   /**< input - subset handle */
  int *maskLength,  /**< output - length of mask array */
  int **mask          /**< output - pointer to output mask array */
) {
  FC_ReturnCode rc;
  int i;  
  _FC_SubSlot* subSlot;

  // default returns
  if(maskLength)
    *maskLength = -1;
  if(mask)
    *mask = NULL;

  // check input
  subSlot = _fc_getSubSlot(subset);
  if ( !subSlot || !maskLength || !mask) { 
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }     

  // FIX? log?

  // lazy member fetching done here
  if (subSlot->header.committed == 1 && subSlot->sia.numVal == 0) {
    rc = _fc_readSubsetMembers(subset);
    if (rc != FC_SUCCESS) {
      fc_printfErrorMessage("Lazy loading of members of subset '%s' failed",
                            subSlot->header.name);
      return rc;
    }
  }

  // allocate memory for resulting mask array, initialize to zeros
  *mask = (int*)calloc(subSlot->maxNumMember, sizeof(int));
  if (!(*mask)) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }

  // and length
  *maskLength = subSlot->maxNumMember;

  // trivial case, no members in subset
  if (subSlot->numMember == 0) 
    return FC_SUCCESS;

  for (i = 0; i < subSlot->numMember; i++) {
    (*mask)[subSlot->sia.vals[i]] = 1; 
  }

  return FC_SUCCESS;
}

/**
 * \ingroup Subset
 * \brief Get the parent sequence of a seq subset
 *
 * \modifications
 * - 03/03/08 ACG created
 */
FC_ReturnCode fc_getSequenceFromSeqSubset(
  int numStep,    /**< input - number of steps in sequence variable */
  FC_Subset *seqSubset, /**< input - sequence subset (array of subsets) */
  FC_Sequence *sequence /**< output - sequence handle */
  ){
  _FC_SubSlot* subSlot;

  //set default value
  if (sequence != NULL)
    *sequence = FC_NULL_SEQUENCE;

  //check input
  if (!fc_isSeqSubsetValid(numStep, seqSubset) || sequence == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  //log message
  subSlot = _fc_getSubSlot(seqSubset[0]);
  fc_printfLogMessage("Getting sequence from seq subset '%s'",
		      subSlot->header.name);
  
  //return requested value
  *sequence = subSlot->sequence;

  return FC_SUCCESS;
}


//@}

/** \name Print subset. */
//-------------------------------
//@{

/**
 * \ingroup  Subset
 * \brief Print subset meta data to stdout, and optionally the member IDs.
 *
 * \modifications 
 *    - 05/24/04 WSK, created.
 */
FC_ReturnCode fc_printSubset(
  FC_Subset subset,   /**< input - the subset */
  char* label,        /**< input - label to prepend subset name with, 
                                     can be NULL */
  int print_members   /**< input - 1 = print the member IDs. */
) {
  FC_ReturnCode rc;
  int i;
  _FC_SubSlot* subSlot;
  int maxNumMember, numMember;
  FC_AssociationType assoc;
  int* members;

  // special case: printing a null handle is not an error
  if (FC_HANDLE_EQUIV(subset, FC_NULL_SUBSET)) {
    printf("Subset: FC_NULL_SUBSET\n");
    fflush(NULL);
    return FC_SUCCESS;
  }

  // check input
  subSlot = _fc_getSubSlot(subset);
  if (subSlot == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // log message
  fc_printfLogMessage("Printing subset '%s'", subSlot->header.name);

  // print label & var name
  if (label)
    printf("%s: ", label);
  else
    printf("Subset: ");
  printf("'%s'\n", subSlot->header.name);
  fflush(NULL);

  // print meta data
  rc = fc_getSubsetInfo(subset, &numMember, &maxNumMember, &assoc);
  if (rc != FC_SUCCESS) {
    fc_printfErrorMessage("Failed to get info for subset '%s'",
                          subSlot->header.name);
    return rc;
  }
  printf("%snumMember = %d, maxNumMember = %d\n", INDENT_STR, 
         numMember, maxNumMember);
  printf("%sassoc = %s\n", INDENT_STR, fc_getAssociationTypeText(assoc));
  fflush(NULL);

  // print values
  if (print_members && numMember > 0) {
    rc = fc_getSubsetMembersAsArray(subset, &numMember, &members);
    if (rc != FC_SUCCESS) {
      fc_printfErrorMessage("Failed to get members for subset '%s'",
                            subSlot->header.name);
      return FC_ERROR;
    }
    
    // print the members
    printf("%sMember IDs:\n", INDENT_STR);
    for (i = 0; i < numMember; i++) 
      printf("%s%s%d:  %6d\n", INDENT_STR, INDENT_STR, i, members[i]);
    fflush(NULL);

    free(members);
  }

  return FC_SUCCESS;
}

/** \name Print Sequence Subset Range. */
//-------------------------------
//@{

/**
 * \ingroup  Subset
 * \brief Print sequence subset meta data to stdout, and optionally the member IDs
 *        (with start and stop index values)
 *
 * \modifications 
 *    - 04/29/08 CDU, created.
 */
FC_ReturnCode fc_printSeqSubsetRange(
  int numStep,          /**< input - number of steps in sequence subset */
  FC_Subset *seqSubset, /**< input - the sequence subset */
  int range_start,      /**< input - offset of first variable to print  
		                     or -1 for first, -2 for last */
  int range_stop,       /**< input - offset of last variable to print 
		                     or -1 for first, -2 for last */
  char *label,          /**< input - label to prepend variable name with, 
		                     can be NULL */
  int print_members     /**< input - print the data values.*/
		                     
){

  FC_ReturnCode rc;
  int i,k;
  _FC_SubSlot* subSlots[numStep];
  int start, stop, mod;
  int maxNumMember, numMember;
  FC_AssociationType assoc;
  int* members;

  // special case: printing a null handle is not an error
  if(numStep <= 0 || seqSubset == NULL) {
    fc_printfWarningMessage("Printing a NULL sequence subset");
    printf("Sequence Subset: FC_NULL_VARIABLE\n");
    return FC_SUCCESS;
  }

  //Check input
  for(i=0; i<numStep; i++){
    subSlots[i] = _fc_getSubSlot(seqSubset[i]);
    if(subSlots[i] == NULL){
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
      return FC_INPUT_ERROR;
    }
  }

  //log message
  fc_printfLogMessage("Printing SeqSubset '%s'", subSlots[0]->header.name); 

  // print label & var name
  if (label)
    printf("%s: ", label);
  else
    printf("Sequence Subset: ");
  printf("'%s'\n", subSlots[0]->header.name);
  fflush(NULL);

  // print meta data
  rc = fc_getSubsetInfo(seqSubset[0], &numMember, &maxNumMember, &assoc);
  if (rc != FC_SUCCESS) {
    fc_printfErrorMessage("Failed to get info for subset '%s'",
                          subSlots[0]->header.name);
    return rc;
  }

  printf("%snumMember = %d, maxNumMember = %d\n", INDENT_STR, 
         numMember, maxNumMember);
  printf("%sassoc = %s\n", INDENT_STR, fc_getAssociationTypeText(assoc));
  fflush(NULL);  

  if( print_members && numStep){

    //Handle special case where start/stop are the first/last steps
    if(range_start<0){
      start = (range_start==-2) ? numStep-1 : 0;
    } else {
      start = (range_start>numStep-1) ? numStep-1 : range_start;
    }
    if(range_stop<0){
      stop = (range_stop==-2) ? numStep-1 : 0;
    } else {
      stop = (range_stop>numStep-1) ? numStep-1 : range_stop;
    }

    //Calculate the step size's direction
    mod   = (start<=stop)   ? 1 : -1;

    for(k=start; (mod>0) ? (k<=stop) : (k>=stop) ; k+= mod){
      rc = fc_getSubsetMembersAsArray(seqSubset[k], &numMember, &members);
      if (rc != FC_SUCCESS) {
	fc_printfErrorMessage("Failed to get members for sequence subset '%s' step %d",
			      subSlots[k]->header.name, k);
	return FC_ERROR;
      }

      // print the members
      printf("%sMember IDs (Step %d of seqSubset '%s'):\n", 
	     INDENT_STR, k, subSlots[0]->header.name);
      for (i = 0; i < numMember; i++) 
	printf("%s%s%d:  %6d\n", INDENT_STR, INDENT_STR, i, members[i]);
      fflush(NULL);
      
      free(members);
    }
  }

  return FC_SUCCESS;  
}
//}@

/** \name Print Sequence Subset*/
//-------------------------------
//@{

/**
 * \ingroup  Subset
 * \brief Print sequence subset meta data to stdout, and optionally the member IDs
 *
 * \modifications 
 *    - 04/29/08 CDU, created.
 */
FC_ReturnCode fc_printSeqSubset(
  int numStep,          /**< input - number of steps in sequence subset */
  FC_Subset *seqSubset, /**< input - the sequence subset */
  char *label,          /**< input - label to prepend variable name with, 
		                     can be NULL */
  int print_members     /**< input - print the data values. */
){
  return fc_printSeqSubsetRange(numStep, seqSubset, -1, -2, label, print_members);
}
//}@



/**
 * \ingroup  Subset
 * \brief Print the subset of a mesh.
 *
 * \description
 * 
 *    If the association is the vertices, this will print the coordinates.
 *    If the association is not the vertices, this will print the
 *    appropriate connectivities.
 *
 *    This routine can be used with a subset on a different mesh,
 *    as long as the meshes are similar enough (i.e. same number of
 *    mesh subentities for the association type).
 *
 * \modifications 
 *    - 11/23/03 WSK Created to replace lots and lots of too specific
 *      printBLAH routines.
 *    - 05/24/04 WSK Changed to printing a subset of a mesh.
 *    - 12/15/04 WSK Fixed up - can handle edges & faces.
 */
FC_ReturnCode fc_printSubsetOfMesh(
  FC_Subset subset,   /**< input - subset */
  FC_Mesh mesh        /**< input - mesh */
) {
  FC_ReturnCode rc;
  int i, j;
  int dim, numMember, maxNumMember, *members, numEntity;
  char* meshName, *subsetName;
  FC_AssociationType assoc;
  FC_ElementType elemType;
  
  // check input
  if (!fc_isSubsetValid(subset) || !fc_isMeshValid(mesh)) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }    

  // get meta data
  rc = fc_getSubsetInfo(subset, &numMember, &maxNumMember, &assoc);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_getMeshNumEntity(mesh, assoc, &numEntity);
  if (rc != FC_SUCCESS)
    return rc;
  
  // a little more checking - have to have same total number of entity
  if (maxNumMember != numEntity) {
    fc_printfErrorMessage("Subset and mesh do not have similar enough "
                          "associations");
    return FC_INPUT_ERROR;
  }

  // get some more info
  rc = fc_getMeshName(mesh, &meshName);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_getSubsetName(subset, &subsetName);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_getSubsetMembersAsArray(subset, &numMember, &members);
  if (rc != FC_SUCCESS)
    return rc;

  // log message
  fc_printfLogMessage("Printing subset '%s' of mesh '%s'", subsetName,
                      meshName);
  
  // print meta data
  printf("Subset '%s' of Mesh '%s'\n", subsetName, meshName);
  printf("%snumMember = %d (of %d), assoc = %s, ", INDENT_STR, numMember, 
         maxNumMember, fc_getAssociationTypeText(assoc));
  fflush(NULL);
  switch(assoc) {
  case FC_AT_EDGE: // fall through to FC_AT_WHOLE_MESH
  case FC_AT_WHOLE_MESH:
    printf("\n");
    break;
  case FC_AT_VERTEX: 
    fc_getMeshDim(mesh, &dim);
    printf("dim = %d\n", dim);
    break;
  case FC_AT_FACE:
    fc_getMeshFaceType(mesh, &elemType);
    printf("faceType = %s\n", fc_getElementTypeText(elemType));
    break;
  case FC_AT_ELEMENT:
    fc_getMeshElementType(mesh, &elemType);
    printf("elemType = %s\n", fc_getElementTypeText(elemType));
    break;
  default:
    return FC_ERROR;
  }
  fflush(NULL);

  // print the data
  if (numMember > 0) {

    if (assoc == FC_AT_WHOLE_MESH) { 
      ; // do nothing
    }

    else if (assoc == FC_AT_VERTEX) { // print coords of members
      double* coords;      
      rc = fc_getMeshCoordsPtr(mesh, &coords);
      if (rc != FC_SUCCESS)
        return rc;
      printf("%sVertex Coordinates:\n", INDENT_STR);
      for (i = 0; i < numMember; i++) {
        printf("%s%s%d:  ", INDENT_STR, INDENT_STR, members[i]);
        for (j = 0; j < dim; j++)
          printf("%11g ", coords[members[i]*dim +j]);
        printf("\n");
      }
      fflush(NULL);
    }

    else { // print connectivities of members
      int* conns;
      int numVertPerEntity = 0, *numVertPerFace = NULL, numVert;

      if (assoc == FC_AT_EDGE) {
        numVertPerEntity = 2;
        rc = fc_getMeshEdgeConnsPtr(mesh, &conns);
        if (rc != FC_SUCCESS)
          return rc;
        printf("%sEdge connectivities:\n", INDENT_STR);
      }
      else if (assoc == FC_AT_FACE && elemType == FC_ET_MIXED) {
        rc = fc_getMeshFaceConnsPtr(mesh, &numVertPerFace, &numVertPerEntity, 
                                    &conns);
        if (rc != FC_SUCCESS)
          return rc;
        printf("%sFace connectivities (3 or 4 vertices per face):\n",
               INDENT_STR);
      }
      else if (assoc == FC_AT_FACE) {
        rc = fc_getMeshFaceConnsPtr(mesh, NULL, &numVertPerEntity, &conns);
        if (rc != FC_SUCCESS)
          return rc;
        printf("%sFace connectivities (%d vertices per face):\n",
               INDENT_STR, numVertPerEntity);
      }
      else {
        numVertPerEntity = fc_getElementTypeNumVertex(elemType);
        rc = fc_getMeshElementConnsPtr(mesh, &conns);
        if (rc != FC_SUCCESS)
          return rc;
        printf("%sElement connectivities (%d vertices per element):\n", 
               INDENT_STR, numVertPerEntity);
      }
      fflush(NULL);

      for (i = 0; i < numMember; i++) {
        printf("%s%s%d:  ", INDENT_STR, INDENT_STR, members[i]);
        if (assoc == FC_AT_FACE && elemType == FC_ET_MIXED) 
          numVert = numVertPerFace[members[i]];
        else 
          numVert = numVertPerEntity;
        for (j = 0; j < numVert; j++)
          printf("%6d ", conns[members[i]*numVertPerEntity + j]);
        printf("\n");
      }
      fflush(NULL);
    }
  } 

  // cleanup
  free(members);
  free(meshName);
  free(subsetName);

  return FC_SUCCESS;
}

/**
 * \ingroup  Subset
 * \brief Print the subset of a variable.
 *
 * \description
 * 
 *    Print the data values on the subset of the mesh. This only works
 *    if the subset and variable share the same association type, and
 *    their parent meshes are similar enough (i.e. same number of
 *    mesh subentities for association type).
 *
 * \modifications 
 *    - 11/26/03 WSK Created.
 *    - 05/24/04 WSK Changed to printing a subset of a variable.
 */
FC_ReturnCode fc_printSubsetOfVariable(
  FC_Subset subset,       /**< input - subset */
  FC_Variable variable    /**< input - variable */
) {
  FC_ReturnCode rc;
  int i, j;
  int numDataPoint, numComp, numMember, maxNumMember, *members;
  char* varName, *subsetName;
  FC_AssociationType assoc, temp_assoc;
  FC_DataType dataType;
  FC_MathType mathType;
  
  // check input
  if (!fc_isSubsetValid(subset) || !fc_isVariableValid(variable)) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // get meta data
  rc = fc_getSubsetInfo(subset, &numMember, &maxNumMember, &temp_assoc);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_getVariableInfo(variable, &numDataPoint, &numComp, &assoc, 
                          &mathType, &dataType);
  if (rc != FC_SUCCESS)
    return rc;

  // a little more checking - have to have same total number of entity
  if (temp_assoc != assoc || maxNumMember != numDataPoint) {
    fc_printfErrorMessage("Subset and mesh do not have similar enough "
                          "associations");
    return FC_INPUT_ERROR;
  }

  rc = fc_getVariableName(variable, &varName);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_getSubsetName(subset, &subsetName);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_getSubsetMembersAsArray(subset, &numMember, &members);
  if (rc != FC_SUCCESS)
    return rc;

  // log message
  fc_printfLogMessage("Printing subset '%s' of variable '%s'", subsetName,
                      varName);

  // print meta data
  printf("Subset '%s' of Variable '%s'\n", subsetName, varName);
  printf("%snumMember = %d (of %d), numComponent = %d\n", INDENT_STR,
         numMember, maxNumMember, numComp);
  printf("%sassoc = %s, mathType = %s, dataType = %s\n", INDENT_STR, 
         fc_getAssociationTypeText(assoc), fc_getMathTypeText(mathType), 
         fc_getDataTypeText(dataType));
  fflush(NULL);
  
  if (numMember > 0) {

    // print the values
    void* dbuf;
    rc = fc_getVariableDataPtr(variable, &dbuf);
    if (rc != FC_SUCCESS)
      return rc;
    printf("%sData:\n", INDENT_STR);
    fflush(NULL);
    switch (dataType) {
    case FC_DT_CHAR: 
      {
        char* cast_buf = (char *)dbuf;
        for (i = 0; i < numMember; i++) {
          printf("%s%s%d:  ", INDENT_STR, INDENT_STR, members[i]);
          for (j = 0; j < numComp; j++)
            printf("%c ", cast_buf[members[i]*numComp + j]);
          printf("\n");
        }
      }
      break;
    case FC_DT_INT:
      {
        int* cast_buf = (int *)dbuf;
        for (i = 0; i < numMember; i++) {
          printf("%s%s%d:  ", INDENT_STR, INDENT_STR, members[i]);
          for (j = 0; j < numComp; j++)
            printf("%6d ", cast_buf[members[i]*numComp + j]);
          printf("\n");
        }
      }
      break;
    case FC_DT_FLOAT:
      {
        float* cast_buf = (float *)dbuf;
        for (i = 0; i < numMember; i++) {
          printf("%s%s%d:  ", INDENT_STR, INDENT_STR, members[i]);
          for (j = 0; j < numComp; j++)
            printf("%11.6g ", cast_buf[members[i]*numComp + j]);
          printf("\n");
        }
      }
      break;
    case FC_DT_DOUBLE:
      {
        double* cast_buf = (double *)dbuf;
        for (i = 0; i < numMember; i++) {
          printf("%s%s%d:  ", INDENT_STR, INDENT_STR, members[i]);
          for (j = 0; j < numComp; j++)
            // print the same way as in fc_printVariable()
            printf("%16.11g ", cast_buf[members[i]*numComp + j]);
          printf("\n");
        }
      }
      break;
    case FC_DT_UNKNOWN:
      printf("Cannot print data type '%s'\n", fc_getDataTypeText(dataType));
    }
    fflush(NULL);
  }

  // cleanup
  free(members);
  free(varName);
  free(subsetName);

  return FC_SUCCESS;
}

//@}

/**
 * \ingroup  PrivateSubset
 * \brief  Initialize a subset slot
 *
 * \modifications  
 *    - 05/11/04 WSK. Created.
 *    - 03/03/08 ACG support sequence subset
 */
void _fc_initSubSlot(
  _FC_SubSlot* subSlot    /**< input - the subset slot to initialize */
) {
  if (subSlot != NULL) {
    _fc_initSlotHeader(&subSlot->header);
    subSlot->mesh = FC_NULL_MESH;
    subSlot->sequence = FC_NULL_SEQUENCE;
    subSlot->stepID = -1;
    subSlot->assoc = FC_AT_UNKNOWN;
    subSlot->numMember = 0;
    subSlot->maxNumMember = 0;
    fc_initSortedIntArray(&subSlot->sia);
  }
}

/**
 * \ingroup  PrivateSubset
 * \brief Release non-necessary resources in a subset slot.
 *
 * \description
 * 
 *    Releases the members
 *
 * \modifications
 *    - 05/11/04 WSK Created.
 */
void _fc_releaseSubSlot(
  _FC_SubSlot* subSlotp       /**< input - subset slot */
) {

  if (subSlotp) {  
    // release big data
    fc_freeSortedIntArray(&subSlotp->sia);
  }
  return;
}

/**
 * \ingroup  PrivateSubset
 * \brief   Clear a subset slot.
 *
 * \description
 *  
 *    Releases all dynamically allocated resources in a variable,
 *    and reinitializes all members.
 *
 * \modifications  
 *    - 05/11/04 WSK Created.
 */
void _fc_clearSubSlot(
 _FC_SubSlot* subSlotp       /**< input - subset slot */
) {
  if (subSlotp != NULL) {
    // free big data
    _fc_releaseSubSlot(subSlotp);

    // free other dynamic arrays, if any

    // clear table header
    _fc_clearSlotHeader(&subSlotp->header);

    // reinit
    _fc_initSubSlot(subSlotp);
  }
}


/**
 * \ingroup  PrivateSubset
 * \brief  Get the size of the subset table.
 *
 * \modifications  
 *   - 05/21/04 WSK Created.
 */
int _fc_getSubTableSize() {
  return subTableSize;
}

/**
 * \ingroup  PrivateSubset
 * \brief Get a subset slot from a subset handle.
 *
 * \description
 * 
 *    The slot id of the subset handle is used to find the
 *    the slot in the subset table. The identity is then
 *    verified by checking that the slot and the handle's
 *    uIDs are the same. On success, a pointer to the slot 
 *    is returned. Otherwise a NULL pointer is returned.
 *
 * \modifications 
 *   - 05/11/04 WSK Created
 *   - 05/25/04 WSK Changed to call new _fc_isHandleValid common to all tables.
 */
_FC_SubSlot* _fc_getSubSlot(
  FC_Subset subset    /**< input - subset */
) {
  if (_fc_isHandleValid(subset.slotID, subset.uID, subTableSize,
                        (_FC_SlotHeader**)subTable))
    return subTable[subset.slotID];
  else
    return NULL;
}

/**
 * \ingroup  PrivateSubset
 * \brief Get a subset slot based on ID.
 *
 * \description
 * 
 *    The _FC_SubSlot at the requested slotID is returned. This should
 *    only be used when you really know that it is the slot you want
 *    because there is no checking of the uID like with _fc_getSubSlot.
 *    The only checking is that the slotID is exists in the table and
 *    that the uID is > 0 (i.e. it is not empty).
 *
 * \modifications 
 *    - 05/11/04 WSK Created.
 *    - 05/25/04 WSK Changed to call new _fc_getSlot() common to all tables.
 */
_FC_SubSlot* _fc_getSubSlotFromID(
  int subID
) {
  return (_FC_SubSlot*) _fc_getSlot(subID, subTableSize, 
                                   (_FC_SlotHeader**)subTable);
}

/**
 * \ingroup  PrivateSubset
 * \brief Delete the subset slot associated with the subset handle.
 *
 * \description 
 *    
 *    This routine clears the dynamic data in the slot and then
 *    deletes the slot from the table.
 *
 * \modifications
 *   - 12/20/2005 WSD. Created.
 */
FC_ReturnCode _fc_deleteSubSlot(
  FC_Subset subset    /**< input - subset */
) {
  _FC_SubSlot* subSlot;
  
  subSlot = _fc_getSubSlot(subset);
  if (subSlot == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // do it
  _fc_clearSubSlot(subSlot);
  return _fc_deleteSlot(subSlot->header.slotID, subTableSize, 
			(_FC_SlotHeader**)subTable, &subOpenSlots);
}

/**
 * \ingroup  PrivateSubset
 * \brief   Get an empty subset slot from the subset table.
 *
 * \description
 *  
 *    Will return a new slot or NULL. A NULL means a memory error occurred.
 *
 * \modifications  
 *    - 5/11/04 WSK  Created.
 */
_FC_SubSlot* _fc_getNewSubSlot() {
  _FC_SubSlot* subSlot;

  subSlot = (_FC_SubSlot*) _fc_getNewSlot(&subTableSize, 
					  (_FC_SlotHeader***)&subTable, 
					  &subOpenSlots, sizeof(_FC_SubSlot));
  _fc_initSubSlot(subSlot);

  return subSlot;
}

/**
 * \ingroup  PrivateSubset
 * \brief  Print the contents of the Subset Table (\ref subTable) to stderr.
 *
 * \description
 *
 *    Prints the contents of the Subset Table in a human readable form.
 *    It takes a string to use as a label for this print out of the table.
 *    Pass NULL if you don't want a label.
 *
 * \modifications 
 *    - 12/8/04 WSK created.
 */
void _fc_printSubTable(
  char *label  /**< Input - label for this version the table (without a
                           trailing \\n) pass NULL if no label is desired */
) {
  int i;
  _FC_SubSlot* subSlot;
  
  // print table heading
  fprintf(stderr, TABLE_DIVIDER_STRING);
  fprintf(stderr, "Subset Table (subTable):\n");
  if (label != NULL)
    fprintf(stderr, "%s\n", label);
  fprintf(stderr, TABLE_DIVIDER_STRING);
  fflush(NULL);

  // We're done if it's empty
  if (subTableSize < 1) {
    fprintf(stderr, "(empty)\n");
    fprintf(stderr, "\n");
    fflush(NULL);
    return;
  }

  // print contents for each slot
  for (i = 0; i < subTableSize; i++) {
    subSlot = subTable[i];
    fprintf(stderr, "%3d: %s\n", i, (subSlot) ? "exists" : "NULL");
    if (subSlot == NULL)
      continue;
    _fc_printSlotHeader(subSlot->header);
    fprintf(stderr, "     FC_Mesh = { %d, %d }\n", 
            subSlot->mesh.slotID, subSlot->mesh.uID);
    fprintf(stderr, "     assoc = %s, numMember = %d, maxNumMember = %d\n",
            fc_getAssociationTypeText(subSlot->assoc), subSlot->numMember,
            subSlot->maxNumMember);
    fprintf(stderr, "     sia.numVal = %d, sia.maxNumVal = %d, sia.vals = %s\n",
	    subSlot->sia.numVal, subSlot->sia.maxNumVal,
	    (subSlot->sia.vals) ? "exists" :"NULL"); 
  }
  fprintf(stderr, "\n");
  fflush(NULL);

  return;
}

/**
 * \ingroup  PrivateSubset
 * \brief  Free all entries in the subset table
 *
 * \description
 *
 *    This should only be called after making sure that
 *    all slots have been cleared.
 *
 * \modifications  
 *   - 05/21/04 WSK Created.
 */
void _fc_freeSubTable() {
  int i;
  for (i = 0; i < subTableSize; i++) {
    _fc_clearSubSlot(subTable[i]);
    free(subTable[i]);
  }
  free(subTable);
  fc_freeSortedIntArray(&subOpenSlots);
  subTableSize = 0;
  subTable = NULL;
}
