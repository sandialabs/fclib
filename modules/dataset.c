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
 * \file dataset.c
 * \brief Implementation for \ref Dataset module.
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/modules/dataset.c,v $
 * $Revision: 1.57 $ 
 * $Date: 2006/10/11 00:19:50 $
 */

// C library includes
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// fc library dependencies
#include "base.h"
#include "storage.h"
#include "library.h"
#include "sequence.h"
#include "mesh.h"
#include "variable.h"
#include "tableP.h"
#include "libraryP.h"
#include "sequenceP.h"
#include "meshP.h"
#include "subsetP.h"
#include "variableP.h"

// this module
#include "dataset.h"
#include "datasetP.h"

/**
 * \addtogroup  Dataset
 * \brief  Operations on datasets
 *
 * \description
 *
 *    (For an explanation of general data manipulations, see 
 *    \ref DataInterface.)
 */

/**
 * \ingroup Dataset
 * \defgroup PrivateDataset (Private)
 */

/** 
 * \ingroup PrivateDataset
 * \brief The number of allocated slots in the dataset table.
 */
static int dsTableSize = 0;
/** 
 * \ingroup PrivateDataset
 * \brief The datasest table.
 */
static _FC_DsSlot **dsTable = NULL;
/**
 * \ingroup PrivateDataset
 * \brief The list of unused slots in the dataset table.
 */
static FC_SortedIntArray dsOpenSlots = { 0, 0, 0 };

/** \name Create new dataset from scratch. */
//-------------------------------
//@{


/**
 * \ingroup  Dataset
 * \brief  Create a new dataset.
 *
 * \description
 *
 *    Create a new, in memory, dataset. 
 *
 * \modifications  
 *    - Nancy Collins, Created.
 */
FC_ReturnCode fc_createDataset(
  char *name,          /**< input - name for new dataset */
  FC_Dataset *dataset   /**< output - new dataset  */
) {
  _FC_DsSlot* dsSlot;
  
  // default return
  if (dataset != NULL)
    *dataset = FC_NULL_DATASET;
  
  // check input
  if (name == NULL || dataset == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // make sure library is ready, return error if it's not
  if (!_fc_isLibraryInit()) {
    fc_printfErrorMessage("Cannot create a dataset until library is initialized.");
    return FC_ERROR;
  }
    
  // log message
  fc_printfLogMessage("Creating new dataset '%s'", name);

  // get an open slot
  dsSlot = _fc_getNewDsSlot();
  if (dsSlot == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  
  // table header information
  _FC_GET_HANDLE(*dataset, dsSlot);
  _fc_setSlotHeaderName(&(dsSlot->header), name);
  
  return FC_SUCCESS;
}

/** \name Other ways to get new datasets. */
//-------------------------------
//@{

/**
 * \ingroup  Dataset
 * \brief  Copy a dataset
 *
 * \description
 *  
 *    Create a copy of a dataset with the given name. This routine will
 *    copy all subentities (sequences & meshes).
 *
 * \modifications 
 *   - 2003-NOV-11 WSK Created
 *   - 9/21/2006 WSD Added copying of global vars.
 */
FC_ReturnCode fc_copyDataset(
  FC_Dataset src_ds,  /**< Input - the source dataset to be copied */
  char *newName,      /**< Input - the name of the new dataset, a NULL value
		       is a flag to use the name of the source dataset */
  FC_Dataset* new_ds  /**< Output - handle to the new dataset */
) { 
  int i, j;
  FC_ReturnCode rc;
  _FC_DsSlot* src_dsSlot;

  // default return
  if (new_ds)
    *new_ds = FC_NULL_DATASET;

  // check input
  src_dsSlot = _fc_getDsSlot(src_ds);
  if (src_dsSlot == NULL || new_ds == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // test for NULL name
  if (newName == NULL)
    newName = src_dsSlot->header.name;

  // log message
  fc_printfLogMessage("Copying dataset '%s' to '%s'", src_dsSlot->header.name,
                      newName);

  // create new dataset with new name
  rc = fc_createDataset(newName, new_ds);
  if (rc != FC_SUCCESS) {
    fc_printfErrorMessage("Failed to create new dataset '%s'",
                          newName);
    return FC_ERROR;
  }

  // copy the sequences
  for (i = 0; i < src_dsSlot->numSeq; i++) {
    FC_Sequence src_seq, new_seq;
    _FC_SeqSlot* src_seqSlot = _fc_getSeqSlotFromID(src_dsSlot->seqIDs[i]);
    _FC_GET_HANDLE(src_seq, src_seqSlot);
    rc = fc_copySequence(src_seq, *new_ds, NULL, &new_seq);
    if (rc != FC_SUCCESS) {
      fc_deleteDataset(*new_ds);
      *new_ds = FC_NULL_DATASET;
      fc_printfErrorMessage("Failed to copy sequence '%s'",
                            src_seqSlot->header.name);
      return rc;
    }
  }

  // copy the global vars
  for (i = 0; i < src_dsSlot->numBasicVar; i++) {
    FC_Variable src_var, new_var;
    _FC_VarSlot* src_varSlot = _fc_getVarSlotFromID(src_dsSlot->basicVarIDs[i]);
    _FC_GET_HANDLE(src_var, src_varSlot);
    rc = fc_copyGlobalVariable(src_var, *new_ds, NULL, &new_var);
    if (rc != FC_SUCCESS) {
      fc_deleteDataset(*new_ds);
      *new_ds = FC_NULL_DATASET;
      fc_printfErrorMessage("Failed to copy global var '%s'",
			    src_varSlot->header.name);
      return rc;
    }
  }
  
  // copy the global seq vars
  // (similar to code in fc_copyMesh())
  for (i = 0; i < src_dsSlot->numSeqVar; i++) {
    int numStep = src_dsSlot->numStepPerSeqVar[i];
    FC_Variable src_seqVar[numStep], *new_seqVar;
    int temp_numSeq;
    FC_Sequence dest_seq, *temp_seqs;
    _FC_VarSlot* src_varSlot = _fc_getVarSlotFromID(src_dsSlot->seqVarIDs[i][0]);
    _FC_SeqSlot* src_seqSlot = _fc_getSeqSlot(src_varSlot->sequence);

    // find the appropriate new sequence from the dest sequences
    rc = fc_getSequences(*new_ds, &temp_numSeq, &temp_seqs);
    for (j = 0; j < temp_numSeq; j++) {
      char* temp_name;
      int temp_numStep;
      FC_DataType temp_dataType;
      void* temp_coords;
      rc = fc_getSequenceName(temp_seqs[j], &temp_name);
      if (rc != FC_SUCCESS)
          return rc;
      rc = fc_getSequenceInfo(temp_seqs[j], &temp_numStep, &temp_dataType);
      if (rc != FC_SUCCESS)
	return rc;
      rc = fc_getSequenceCoordsPtr(temp_seqs[j], &temp_coords);
      if (rc != FC_SUCCESS)
          return rc;
      if (!strcmp(temp_name, src_seqSlot->header.name) &&
	  temp_numStep == src_seqSlot->numStep &&
	  temp_dataType == src_seqSlot->datatype &&
	  !memcmp(temp_coords, src_seqSlot->coords, src_seqSlot->numStep * 
		  fc_sizeofDataType(src_seqSlot->datatype))) {
	  free(temp_name);
	  dest_seq = temp_seqs[j];
	  break;
      }
      free(temp_name);
    }
    free(temp_seqs);

    // assemble handles
    for (j = 0; j < numStep; j++) {
      src_varSlot = _fc_getVarSlotFromID(src_dsSlot->seqVarIDs[i][j]);
      _FC_GET_HANDLE(src_seqVar[j], src_varSlot);
    }


    rc = fc_copyGlobalSeqVariable(numStep, src_seqVar, *new_ds, dest_seq,
				  NULL, &new_seqVar);
    if (rc != FC_SUCCESS) {
      fc_deleteDataset(*new_ds);
        *new_ds = FC_NULL_DATASET;
      fc_printfErrorMessage("Failed to copy global var '%s'",
			    src_varSlot->header.name);
      return rc;
    }
    free(new_seqVar);
  }

  // copy the meshes
  for (i = 0; i < src_dsSlot->numMesh; i++) {
    FC_Mesh src_mesh, new_mesh;
    _FC_MeshSlot* src_meshSlot = _fc_getMeshSlotFromID(src_dsSlot->meshIDs[i]);
    _FC_GET_HANDLE(src_mesh, src_meshSlot);
    rc = fc_copyMesh(src_mesh, *new_ds, NULL, 1, 1, 1, 1, &new_mesh);
    if (rc != FC_SUCCESS) {
      fc_deleteDataset(*new_ds);
      *new_ds = FC_NULL_DATASET;
      fc_printfErrorMessage("Failed to copy mesh '%s'",
                            src_meshSlot->header.name);
      return rc;
    }
  }

  return FC_SUCCESS;
}

//@}

/** \name Get existing datasets. */
//-------------------------------
//@{

/**
 * \ingroup  Dataset
 * \brief  Return a list of datasets in the library.
 *
 * \description
 *  
 *    Returns an array of handles for the existing datasets. The user is 
 *    responsible for freeing the array of handles.
 *
 * \modifications  
 *   - 9/8/04 WSK Created.
 */
FC_ReturnCode fc_getDatasets(
 int *numDataset,       /**< output - number of existing datasets */
 FC_Dataset **datasets  /**< output - array of dataset (optional) */
) {
  int i;
  int num;
  _FC_DsSlot* dsSlot;
  
  // default return values
  if (numDataset != NULL) 
    *numDataset = -1;
  if (datasets != NULL)
    *datasets = NULL;  
  
  // call on uninitialized library is an error
  if (!_fc_isLibraryInit()) {
    fc_printfErrorMessage("Library has not been initialized");
    return FC_ERROR;
  }

  // check input
  if (numDataset == NULL || datasets == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Getting datasets");
  
  // count datasets
  num = 0;
  for (i = 0; i < _fc_getDsTableSize(); i++) {
    dsSlot = _fc_getDsSlotFromID(i);
    if (dsSlot) 
      num++;
  }
  *numDataset = num;

  // get handles
  if (num > 0) {
    *datasets = (FC_Dataset*)malloc(num*sizeof(FC_Dataset));
    if (*datasets == NULL) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      return FC_MEMORY_ERROR;
    }
    num = 0;
    for (i = 0; i < _fc_getDsTableSize(); i++) {
      dsSlot = _fc_getDsSlotFromID(i);
      if (dsSlot) {
        _FC_GET_HANDLE((*datasets)[num], dsSlot);
        num++;
      }
    }
  }
  
  return FC_SUCCESS;
}

/**
 * \ingroup  Dataset
 * \brief  Return number of datasets in the library.
 *
 * \modifications  
 *   - 9/28/04 WSK Created.
 */
FC_ReturnCode fc_getNumDataset(
 int *numDataset      /**< output - number of existing datasets */
) {
  int i, num;
  _FC_DsSlot* dsSlot;
  
  // default return values
  if (numDataset != NULL) 
    *numDataset = -1;
  
  // call on uninitialized library is an error
  if (!_fc_isLibraryInit()) {
    fc_printfErrorMessage("Library has not been initialized");
    return FC_ERROR;
  }

  // check input
  if (numDataset == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Getting number of datasets");
  
  // count datasets
  num = 0;
  for (i = 0; i < _fc_getDsTableSize(); i++) {
    dsSlot = _fc_getDsSlotFromID(i);
    if (dsSlot) 
      num++;
  }
  *numDataset = num;
  
  return FC_SUCCESS;
}

/**
 * \ingroup  Dataset
 * \brief  Gets all datasets with a given name
 *
 * \description
 *
 *    Return the handles of all datasets with the given name.  
 *
 * \modifications  
 *   - 9/8/04 WSK, created.
 *   - 9/21/06 ACG now returns multiple matching items.
 */
FC_ReturnCode fc_getDatasetByName(
  char *datasetName,    /**< input - name of dataset to find */
  int *numDatasets,     /**< output - num matching datasets */
  FC_Dataset **datasets   /**< output - dataset handle */
) {
  int i;
  _FC_DsSlot* dsSlot;
  FC_Dataset dataset;
  FC_Dataset *lmatch = NULL;
  int numlmatch, maxnumlmatch;
  
  // default values returned if anything goes wrong
  if (numDatasets)
    *numDatasets = -1;
  if (datasets)
    *datasets = NULL;
  
  // call on uninitialized library is an error
  if (!_fc_isLibraryInit()) {
    fc_printfErrorMessage("Library has not been initialized");
    return FC_ERROR;
  }

  // check input
  if (datasetName == NULL || !datasets || !numDatasets) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // log message
  fc_printfLogMessage("Getting dataset '%s'", datasetName);

  numlmatch = 0;
  maxnumlmatch = 0;
  
  // look for existing slot with the matching name
  for (i = 0; i < _fc_getDsTableSize(); i++) {
    dsSlot = _fc_getDsSlotFromID(i);
    if (dsSlot) {
      if (!strcmp(datasetName, dsSlot->header.name)) {
        _FC_GET_HANDLE(dataset, dsSlot);

	if (lmatch == NULL){
	  lmatch = (FC_Dataset*)malloc(sizeof(FC_Dataset));
	  if (lmatch == NULL){
	    fc_printfErrorMessage("%s", 
				  fc_getReturnCodeText(FC_MEMORY_ERROR));
	    return FC_MEMORY_ERROR;
	  }
	  numlmatch = 1;
	  maxnumlmatch = 1;
	  lmatch[0] = dataset;
	}else{
	  if(numlmatch == maxnumlmatch){
	    FC_Dataset *temp;
	    temp = (FC_Dataset*)realloc(lmatch,
					(2*numlmatch)*sizeof(FC_Dataset));
	    if (temp == NULL){
	      fc_printfErrorMessage("%s",
				    fc_getReturnCodeText(FC_MEMORY_ERROR));
	      return FC_MEMORY_ERROR;
	    }
	    maxnumlmatch*=2;
	    lmatch = temp;
	  }
	  lmatch[numlmatch] = dataset;
	  numlmatch++;
	}
      }
    }
  }

  *numDatasets = numlmatch;
  *datasets = lmatch;

  return FC_SUCCESS;
}

//@}

/** \name Change the name of a dataset. */
//-------------------------------
//@{

/**
 * \ingroup Dataset
 * \brief  Change the name of a dataset.
 *
 * \description
 *
 *    This does not affect the name of the file the dataset is associated with.
 *
 * \modifications
 *   - 8/24/2005 WSD. Created.
 */
FC_ReturnCode fc_changeDatasetName(
  FC_Dataset dataset, /**< Input - the dataset. */
  char* newName       /**< Input - new name for the dataset. */
) {
  _FC_DsSlot* dsSlot;

  // check input
  dsSlot = _fc_getDsSlot(dataset);
  if (dsSlot == NULL || newName == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // log message
  fc_printfLogMessage("Changing name of dataset from '%s' to '%s'", 
                      dsSlot->header.name, newName);
  
  // Do it
  return _fc_setSlotHeaderName(&dsSlot->header, newName);
}

//@}

/** \name Release/delete datasets. */
//-------------------------------
//@{

/**
 * \ingroup Dataset
 * \brief Attempt to minimize dataset's memory usage.
 *
 * \description
 *
 *    Call to try to release un-needed memory in a dataset. It the dataset has
 *    been saved to disk, the large data is released. If the dataset has not be
 *    saved to disk, all large data is saved.
 *
 *    This function will recursively call release on it's child sequences and
 *    meshes (and they will call their children).
 *
 * \modifications
 *    - 9/27/04 WSK Created.
 */
FC_ReturnCode fc_releaseDataset(
 FC_Dataset dataset       /**< input - dataset handle */
) {
  int i;
  FC_Sequence *sequences;
  FC_Mesh *meshes;
  FC_Variable *variables, **seqVariables;
  _FC_DsSlot* dsSlot;
  int numSequence, numMesh, numVar, numSeqVar, *numSteps;
  
  // special case: releasing a null handle is not an error
  if (FC_HANDLE_EQUIV(dataset, FC_NULL_DATASET)) {
    fc_printfWarningMessage("Releasing FC_NULL_DATASET");
    return FC_SUCCESS;
  }

  // check input
  dsSlot = _fc_getDsSlot(dataset);
  if (dsSlot == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Cleaning up dataset '%s'", dsSlot->header.name);

  // call release on children sequences & meshes & global vars
  fc_getSequences(dataset, &numSequence, &sequences);
  for (i = 0; i < numSequence; i++) 
    fc_releaseSequence(sequences[i]); 
  free(sequences);
  fc_getMeshes(dataset, &numMesh, &meshes);
  for (i = 0; i < numMesh; i++) 
    fc_releaseMesh(meshes[i]);
  free(meshes);
  fc_getGlobalVariables(dataset, &numVar, &variables);
  for (i = 0; i < numVar; i++)
    fc_releaseVariable(variables[i]);
  free(variables);
  fc_getGlobalSeqVariables(dataset, &numSeqVar, &numSteps, &seqVariables);
  for (i = 0; i < numSeqVar; i++)
    fc_releaseSeqVariable(numSteps[i], seqVariables[i]);
  free(numSteps);
  for (i = 0; i < numSeqVar; i++)
    free(seqVariables[i]);
  free(seqVariables);
  
  // no big data
  //if (dsSlot->header.committed == 1) {
  //}

  return FC_SUCCESS;
}

/**
 * \ingroup  Dataset
 * \brief  Delete a dataset.
 *
 * \description
 *  
 *    Call when this dataset is no longer needed.  It deletes
 *    the dataset and all of it's children including
 *    it's meshes, variables, etc.
 *
 * \modifications  
 *   - Nancy Collins  Created.
 *   - Aug 5, 2002  W Koegler  Now calls _fc_killDsSlot so that 
 *     fc_closeDataset and fc_openDataset can use the same code to 
 *     release the dsTable slot.
 *   - 11/20/03 RM, added more input checking, if dataset closed with default 
 *     handle \ref FC_NULL_DATASET then do nothing and return success (issue 
 *     warning if fc_getLibraryVerbosity()).
 *   - 9/21/2006 WSD added deleting of member global vars
 */
FC_ReturnCode fc_deleteDataset(
 FC_Dataset dataset  /**< input - dataset to be deleted */
) {
  FC_ReturnCode rc;
  int i;
  _FC_DsSlot* dsSlot;
  int numSequence, numMesh, numVar, numSeqVar, *numStepPerVar;
  FC_Sequence *sequences;
  FC_Mesh *meshes;
  FC_Variable *variables, **seqVariables;
  
  // special case: deleting a null handle is not an error
  if ( FC_HANDLE_EQUIV(dataset, FC_NULL_DATASET)) {
    fc_printfWarningMessage("Deleting FC_NULL_DATASET.");
    return FC_SUCCESS;
  }
  
  // check input
  dsSlot = _fc_getDsSlot(dataset);
  if (dsSlot == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Deleting dataset '%s'", dsSlot->header.name);

  // free up any file IO related stuff
  if (dsSlot->header.committed) {
    rc = _fc_closeDataset(dataset);
    if (rc != FC_SUCCESS)
      return rc;
  }

  // delete member global vars
  rc = fc_getGlobalVariables(dataset, &numVar, &variables);
  if (rc != FC_SUCCESS)
    return rc;
  for (i = 0; i < numVar; i++) {
    rc = fc_deleteVariable(variables[i]);
    if (rc != FC_SUCCESS)
      return rc;
  }
  free(variables);

  // delete member global seq vars
  rc = fc_getGlobalSeqVariables(dataset, &numSeqVar, &numStepPerVar, 
				&seqVariables);
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
  
  // delete member meshes
  rc = fc_getMeshes(dataset, &numMesh, &meshes);
  if (rc != FC_SUCCESS)
    return rc;
  for (i = 0; i < numMesh; i++) {
    rc = fc_deleteMesh(meshes[i]);
    if (rc != FC_SUCCESS)
      return rc;
  }
  free(meshes);

  // delete member sequences
  rc = fc_getSequences(dataset, &numSequence, &sequences);
  if (rc != FC_SUCCESS)
    return rc;
  for (i = 0; i < numSequence; i++) {
    rc = fc_deleteSequence(sequences[i]);
    if (rc != FC_SUCCESS)
      return rc;
  }
  free(sequences);
  
  // clear & delete slot
  return _fc_deleteDsSlot(dataset);
}

//@}

/** \name Get dataset metadata. */
//-------------------------------
//@{

/**
 * \ingroup  Dataset
 * \brief Check that the handle refers to a valid dataset.
 *
 * \description
 *
 *    Returns 1 (true) if the handle is valid, and 0 (false) if it is not.
 *
 * \modifications  
 *    - 05/25/04 WSK Created.
 */
int fc_isDatasetValid(
  FC_Dataset dataset  /**< input - dataset handle */
) {
  return _fc_isHandleValid(dataset.slotID, dataset.uID, dsTableSize,
                           (_FC_SlotHeader**)dsTable);
}

/**
 * \ingroup  Dataset
 * \brief Return name of a dataset
 *
 * \description
 *  
 *    Return the name associated with a dataset handle. Spaced
 *    is allocated for the string so it needs to be freed by the
 *    user when finished.
 *
 * \modifications  
 *    - Aug, 9, 2002.  W Koegler  Created.
 *    - 11/21/03 RM, added checking for NULL dataset handle
 *        default return filename as NULL if anything is wrong
 */
FC_ReturnCode fc_getDatasetName(
  FC_Dataset dataset,  /**< input - dataset handle */
  char **dsName        /**< output - dataset name, space allocated here
                             and must be freed by user when finished */
) {
  char *cp;
  _FC_DsSlot* dsSlot;
  
  // default return value
  if (dsName)
    *dsName = NULL;
  
  // check input
  dsSlot = _fc_getDsSlot(dataset);
  if (dsSlot == NULL || dsName == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Getting name for dataset '%s'", dsSlot->header.name);

  // copy name
  cp = dsSlot->header.name;
  *dsName = malloc(strlen(cp) + 1);
  if (*dsName == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  strcpy(*dsName, cp);
  
  return FC_SUCCESS;
}

//@}

/** \name Print dataset. */
//-------------------------------
//@{

/**
 * \ingroup  Dataset
 * \brief Print dataset meta data to stdout.
 *
 * \modifications 
 *    - 11/20/03 WSK Created.
 */
FC_ReturnCode fc_printDataset(
  FC_Dataset dataset,   /**< input - dataset */
  char* label           /**< input - label to prepend dataset name with, 
                                     can be NULL */
) {
  FC_ReturnCode rc;
  int numMesh, numSeq;
  int numGlbVar, numGlbSeqVar;
  _FC_DsSlot* dsSlot;

  // special case: printing a null handle is not an error
  if (FC_HANDLE_EQUIV(dataset, FC_NULL_DATASET)) {
    printf("Dataset: FC_NULL_DATASET\n");
    fflush(NULL);
    return FC_SUCCESS;
  }

  // check input
  dsSlot = _fc_getDsSlot(dataset);
  if (dsSlot == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // log message
  fc_printfLogMessage("Printing dataset '%s'", dsSlot->header.name);

  // print label & var name
  if (label)
    printf("%s: ", label);
  else
    printf("Dataset: ");
  printf("'%s'\n", dsSlot->header.name);
  fflush(NULL);

  // print meta data
  rc = fc_getNumSequence(dataset, &numSeq);
  if (rc != FC_SUCCESS) {
    fc_printfErrorMessage("Failed to get sequences for dataset '%s'",
                          dsSlot->header.name);
    return rc;
  }
  rc = fc_getNumGlobalVariable(dataset, &numGlbVar);
  if (rc != FC_SUCCESS) {
    fc_printfErrorMessage("Failed to get global vars for dataset '%s'",
                          dsSlot->header.name);
    return rc;
  }
  rc = fc_getNumGlobalSeqVariable(dataset, &numGlbSeqVar);
  if (rc != FC_SUCCESS) {
    fc_printfErrorMessage("Failed to get global seq vars for dataset '%s'",
                          dsSlot->header.name);
    return rc;
  }
  rc = fc_getNumMesh(dataset, &numMesh);
  if (rc != FC_SUCCESS) {
    fc_printfErrorMessage("Failed to get meshes for dataset '%s'",
                          dsSlot->header.name);
    return rc;
  }
  printf("%snumSeq = %d, numGlbVar = %d, numGlbSeqVar = %d, numMesh = %d\n",
         INDENT_STR, numSeq, numGlbVar, numGlbSeqVar, numMesh);
  fflush(NULL);

  return FC_SUCCESS;
}

//@}

/**
 * \ingroup  PrivateDataset
 * \brief  Initialize a dataset slot
 *
 * \modifications  
 *   APR-30-2003  W Koegler  Created. 
 */
void _fc_initDsSlot(
  _FC_DsSlot* dsSlot    /**< input - the dataset slot to initialize */
) {
  if (dsSlot != NULL) {
    _fc_initSlotHeader(&dsSlot->header);
    dsSlot->fileType = FC_FT_NONE;
    dsSlot->baseFile = NULL;
    dsSlot->writable = 1;
    dsSlot->numSeq = 0;
    dsSlot->seqIDs = NULL;
    dsSlot->numMesh = 0;
    dsSlot->meshIDs = NULL;
    dsSlot->numBasicVar = 0;
    dsSlot->basicVarIDs = NULL;
    dsSlot->numSeqVar = 0;
    dsSlot->numStepPerSeqVar = 0;
    dsSlot->seqVarIDs = 0;
  }   
}

/**
 * \ingroup  PrivateDataset
 * \brief  Clear a dataset slot
 *
 * \description
 * 
 *    Releases all dynamically allocated resources in a dataset,
 *    and reinitializes all members.
 *
 * \modifications  
 *   - 2003-APR-17  W Koegler  Created
 *   - 2003-APR-31  W Koegler  Made more comprehensive
 */
void _fc_clearDsSlot(
 _FC_DsSlot* dsSlotp       /**< input - dtable pointer */
) {
  int i;

  if (dsSlotp != NULL) {
    // free strings
    if (dsSlotp->baseFile)
      free(dsSlotp->baseFile);

    // free dynamic arrays
    if (dsSlotp->meshIDs)
      free(dsSlotp->meshIDs);
    if (dsSlotp->seqIDs)
      free(dsSlotp->seqIDs);
    if (dsSlotp->basicVarIDs)
      free(dsSlotp->basicVarIDs);
    if (dsSlotp->numStepPerSeqVar)
      free(dsSlotp->numStepPerSeqVar);
    if (dsSlotp->seqVarIDs) {
      for (i = 0; i < dsSlotp->numSeqVar; i++)
	free(dsSlotp->seqVarIDs[i]);
      free(dsSlotp->seqVarIDs);
    }

    // clear table header
    _fc_clearSlotHeader(&dsSlotp->header);

    // reinit
    _fc_initDsSlot(dsSlotp);
  }
}

/**
 * \ingroup  PrivateDataset
 * \brief  Get the size of the dataset table.
 *
 * \modifications  
 *   - 05/21/04 WSK Created.
 */
int _fc_getDsTableSize() {
  return dsTableSize;
}

/**
 * \ingroup  PrivateDataset
 * \brief   Get an empty dataset slot from the dataset table.
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
_FC_DsSlot* _fc_getNewDsSlot() {
  _FC_DsSlot* dsSlot;

  dsSlot = (_FC_DsSlot*) _fc_getNewSlot(&dsTableSize, 
					(_FC_SlotHeader***)&dsTable, 
					&dsOpenSlots, sizeof(_FC_DsSlot));
  _fc_initDsSlot(dsSlot);

  return dsSlot;
}

/**
 * \ingroup  PrivateDataset
 * \brief Get a dataset slot from a dataset handle.
 *
 * \description
 * 
 *    The slot id of the dataset handle is used to find the
 *    the slot in the dataset table. The identity is then
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
_FC_DsSlot* _fc_getDsSlot(
  FC_Dataset dataset    /**< input - dataset */
) {
  if (_fc_isHandleValid(dataset.slotID, dataset.uID, dsTableSize,
                        (_FC_SlotHeader**)dsTable))
    return dsTable[dataset.slotID];
  else
    return NULL;
}

/**
 * \ingroup  PrivateDataset
 * \brief Get a dataset slot based on ID.
 *
 * \description
 * 
 *    The _FC_DsSlot at the requested slotID is returned. This should
 *    only be used when you really know that it is the slot you want
 *    because there is no checking of the uID like with _fc_getDsSlot.
 *    The only checking is that the slotID is exists in the table and
 *    that the uID is > 0 (i.e. it is not empty).
 *
 * \modifications 
 *    - 11/27/03 WSK Created so that external routines can index into
 *      tables without knowing about them.
 *    - 05/25/04 WSK Changed to call new _fc_getSlot() common to all tables.
 */
_FC_DsSlot* _fc_getDsSlotFromID(
  int dsID
) {
  return (_FC_DsSlot*) _fc_getSlot(dsID, dsTableSize, 
                                   (_FC_SlotHeader**)dsTable);
}

/**
 * \ingroup  PrivateDataset
 * \brief Delete the dataset slot associated with the dataset handle.
 *
 * \description 
 *    
 *    This routine clears the dynamic data in the slot and then
 *    deletes the slot from the table.
 *
 * \modifications
 *   - 12/20/2005 WSD. Created.
 */
FC_ReturnCode _fc_deleteDsSlot(
  FC_Dataset dataset    /**< input - dataset */
) {
  _FC_DsSlot* dsSlot;
  
  dsSlot = _fc_getDsSlot(dataset);
  if (dsSlot == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // do it
  _fc_clearDsSlot(dsSlot);
  return _fc_deleteSlot(dsSlot->header.slotID, dsTableSize, 
			(_FC_SlotHeader**)dsTable, &dsOpenSlots);
}

/**
 * \ingroup  PrivateDataset
 * \brief  Print the contents of the Dataset Table (\ref dsTable) to stderr.
 *
 * \description
 *
 *    Prints the contents of the Datatset Table in a human readable form.
 *    It takes a string to use as a label for this print out of the table.
 *    Pass NULL if you don't want a label.
 *
 * \modifications 
 *    - Nancy Collins Created.
 *    - 2003-NOV-13  WSK  Fixed up.
 */
void _fc_printDsTable(
  char *label /**< Input - label for this version the table (without a
                           trailing \\n) pass NULL if no label is desired */
) {
  int i, j, k;
  _FC_DsSlot* dsSlot;

  // print table heading
  fprintf(stderr, TABLE_DIVIDER_STRING);
  fprintf(stderr, "Dataset Table (dsTable):\n");
  if (label != NULL)
    fprintf(stderr, "%s\n", label);
  fprintf(stderr, TABLE_DIVIDER_STRING);
  fflush(NULL);
  
  // We're done if it's empty
  if (dsTableSize < 1) {
    fprintf(stderr, "(empty)\n");
    fprintf(stderr, "\n");
    fflush(NULL);
    return;
  }

  // print contents for each slot 
  for (i = 0; i < dsTableSize; i++) {
    dsSlot = dsTable[i];
    fprintf(stderr, "%3d: %s\n", i, (dsSlot) ? "exists" : "NULL");
    if (dsSlot == NULL)
      continue;
    _fc_printSlotHeader(dsSlot->header);
    fprintf(stderr, "     fileType = %s\n", 
            fc_getFileIOTypeText(dsSlot->fileType));
    fprintf(stderr, "     baseFile = %s\n", (dsSlot->baseFile) ?
            dsSlot->baseFile : "NULL");
    fprintf(stderr, "     writable = %d\n", dsSlot->writable);
    fprintf(stderr, "     numSeq = %d, seqIDs = ", dsSlot->numSeq);
    if (dsSlot->numSeq == 0) 
      fprintf(stderr, "NULL\n");
    else {
      fprintf(stderr, "{ %d", dsSlot->seqIDs[0]);
      for (j = 1; j < dsSlot->numSeq; j++)
        fprintf(stderr, ", %d", dsSlot->seqIDs[j]);
      fprintf(stderr, " }\n");
    }
    fprintf(stderr, "     numMesh = %d, meshIDs = ", dsSlot->numMesh);
    if (dsSlot->numMesh == 0) 
      fprintf(stderr, "NULL\n");
    else {
      fprintf(stderr, " { %d", dsSlot->meshIDs[0]);
      for (j = 1; j < dsSlot->numMesh; j++)
        fprintf(stderr, ", %d", dsSlot->meshIDs[j]);
      fprintf(stderr, " }\n");
    }
    fprintf(stderr, "     numBasicVar = %d, basicVarIDs = ",
	    dsSlot->numBasicVar);
    if (dsSlot->numBasicVar == 0) 
      fprintf(stderr, "NULL\n");
    else {
      fprintf(stderr, "{ %d", dsSlot->basicVarIDs[0]);
      for (j = 1; j < dsSlot->numBasicVar; j++)
        fprintf(stderr, ", %d", dsSlot->basicVarIDs[j]);
      fprintf(stderr, " }\n");
    }
    fprintf(stderr, "     numSeqVar = %d, numStepPerSeqVar = ", dsSlot->numSeqVar);
    if (!dsSlot->numStepPerSeqVar)
      fprintf(stderr, "NULL\n");
    else {
      fprintf(stderr, "{ %d", dsSlot->numStepPerSeqVar[0]);
      for (j = 1; j < dsSlot->numSeqVar; j++)
        fprintf(stderr, ", %d", dsSlot->numStepPerSeqVar[j]);
      fprintf(stderr, " }\n");
    }
    fprintf(stderr, "     seqVarIDs = ");
    if (!dsSlot->seqVarIDs) 
      fprintf(stderr, "NULL\n");
    else {
      fprintf(stderr, "{\n");
      for (j = 0; j < dsSlot->numSeqVar; j++) {
	fprintf(stderr, "         seqVarIDs[%d] = { %d", j, 
                  dsSlot->seqVarIDs[j][0]);
	for (k = 1; k < dsSlot->numStepPerSeqVar[j]; k++)
	  fprintf(stderr, ", %d", dsSlot->seqVarIDs[j][k]);
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
 * \ingroup  PrivateDataset
 * \brief  Free all entries in the dataset table
 *
 * \description
 *
 *    This should only be called after making sure that
 *    all slots have been cleared.
 *
 * \modifications  
 *   - 05/21/04 WSK Created.
 */
void _fc_freeDsTable() {
  int i;
  for (i = 0; i < dsTableSize; i++) {
    _fc_clearDsSlot(dsTable[i]);
    free(dsTable[i]);
  }
  free(dsTable);
  fc_freeSortedIntArray(&dsOpenSlots);
  dsTableSize = 0;
  dsTable = NULL;
}
