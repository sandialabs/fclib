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
 * \file variable.c
 * \brief Implementation for \ref Variable module.
 *
 * $Source: /home/Repositories/fcdmf/fclib/modules/variable.c,v $
 * $Revision: 1.76 $ 
 * $Date: 2007/02/28 20:18:28 $
 */

// C library includes
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

// fc library dependencies
#include "base.h"
#include "storage.h"
#include "library.h"
#include "dataset.h"
#include "sequence.h"
#include "mesh.h"
#include "tableP.h"
#include "datasetP.h"
#include "sequenceP.h"
#include "meshP.h"
#include "fileioP.h"
#include "topo.h"

// this module
#include "variable.h"
#include "variableP.h"

/**
 * \addtogroup  Variable
 * \brief  Operations on variables.
 *
 * \description
 *
 *    Sequence variables are stored as an array of variables. If a function
 *    does not have a seqVariable version, in many cases you can call the
 *    variable version on the first variable in the seqVariable array and
 *    get the results you expect, e.g. fc_getVariableInfo(seqVar[0], ...).
 *    In some cases you would need to loop over all variables in the
 *    seqVariable array, e.g. fc_setVariableData(seqVar[i], ...).
 *    (For an explanation of general data manipulations, see 
 *    \ref DataInterface.)
 */

/**
 * \ingroup  Variable
 * \defgroup   PrivateVariable (Private)
 */

/** 
 * \ingroup PrivateVariable
 * \brief The number of allocated slots in the variable table.
 */
static int varTableSize = 0;
/** 
 * \ingroup PrivateVariable
 * \brief The variable table.
 */
static _FC_VarSlot **varTable = NULL;
/**
 * \ingroup PrivateVariable
 * \brief The list of unused slots in the variable table.
 */
static FC_SortedIntArray varOpenSlots = { 0, 0, 0 };

/** \name Create new variable from scratch. */
//-------------------------------
//@{

/**
 * \ingroup Variable
 * \brief  Create a new global variable.
 *
 * \description
 *
 *    Create a new, in-memory, variable associated with the entire dataset.
 *
 * \modifications  
 *    - 9/19/2006 WSD Created. Essentially a copy of fc_createVariable,
 *        except parent is the dataset.
 */
FC_ReturnCode fc_createGlobalVariable(
  FC_Dataset dataset,    /**< input - dataset handle */
  char *varName,         /**< input - variable name */
  FC_Variable *variable  /**< output - variable handle */
) {
  int *tempslotIDs;
  _FC_DsSlot* dsSlot;
  _FC_VarSlot* varSlot;
  
  // default return
  if (variable)
    *variable = FC_NULL_VARIABLE;
  
  // check input
  dsSlot = _fc_getDsSlot(dataset);
  if (dsSlot == NULL || !varName || !variable) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Creating new global variable '%s'", varName);
  
  // get an open slot
  varSlot = _fc_getNewVarSlot();
  if (varSlot == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  
  // table header information
  _FC_GET_HANDLE(*variable, varSlot);
  _fc_setSlotHeaderName(&varSlot->header, varName);
  
  // set back and forth references
  varSlot->ds = dataset;
  tempslotIDs = realloc(dsSlot->basicVarIDs, 
                        sizeof(int) * (dsSlot->numBasicVar+1));
  if (tempslotIDs == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  dsSlot->basicVarIDs = tempslotIDs;
  dsSlot->basicVarIDs[dsSlot->numBasicVar] = varSlot->header.slotID;
  dsSlot->numBasicVar++;
  
  return FC_SUCCESS;
}

/**
 * \ingroup  Variable 
 * \brief  Create a new global sequence variable.
 *
 * \description
 *  
 *    Creates a new, in-memory, global sequence variable as an array
 *    of variables. A sequence must be provided in addition
 *    to the dataset that will be the owner of the variable.
 *    The sequence must be in the provide dataset.
 *
 * \modifications  
 *    - 9/21/2006 Created.
 */
FC_ReturnCode fc_createGlobalSeqVariable(
  FC_Dataset dataset,    /**< input - dataset handle */
  FC_Sequence sequence,  /**< input - sequence handle */
  char *varName,         /**< input - sequence variable name */
  int* numStep,          /**< output - the number of steps */
  FC_Variable **seqVars  /**< output - seq variable handles */
) {
  int i;
  _FC_DsSlot* dsSlot;
  _FC_SeqSlot* seqSlot;
  _FC_VarSlot* varSlot;
  int *tempNumPer, **tempSlotIDs;
  
  // default output
  if (numStep)
    *numStep = -1;
  if (seqVars)
    *seqVars = NULL;
  
  // check input
  dsSlot = _fc_getDsSlot(dataset);
  seqSlot = _fc_getSeqSlot(sequence);
  if (dsSlot == NULL || seqSlot == NULL || varName == NULL || 
      numStep == NULL || seqVars == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  if (!FC_HANDLE_EQUIV(dataset, seqSlot->ds)) {
    fc_printfErrorMessage("Mesh and sequence need to be in same dataset");
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Creating new global seq variable '%s'", varName);
  
  // make room for array & setup return values
  *numStep = seqSlot->numStep;
  *seqVars = malloc(sizeof(FC_Variable)*seqSlot->numStep);
  if (*seqVars == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  
  for (i = 0; i < seqSlot->numStep; i++) {
    // get an open slot
    varSlot = _fc_getNewVarSlot();
    if (varSlot == NULL) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));      
      return FC_MEMORY_ERROR;
    }
    
    // table header information
    _FC_GET_HANDLE((*seqVars)[i], varSlot);
    _fc_setSlotHeaderName(&varSlot->header, varName);
    
    // set back and forth references
    varSlot->ds = dataset;
    varSlot->sequence = sequence;
    varSlot->stepID = i;
  }
  
  // set more back and forth references
  tempNumPer = (int*)realloc(dsSlot->numStepPerSeqVar, 
			     (dsSlot->numSeqVar+1)*sizeof(int));
  tempSlotIDs = (int**)realloc(dsSlot->seqVarIDs,
			       (dsSlot->numSeqVar+1)*sizeof(int*));
  if (!tempNumPer || !tempSlotIDs) {
    free(tempNumPer); free(tempSlotIDs);
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  tempSlotIDs[dsSlot->numSeqVar] = (int*)malloc(seqSlot->numStep*sizeof(int)); 
  if (!tempSlotIDs[dsSlot->numSeqVar]) {
    free(tempNumPer); free(tempSlotIDs);
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  dsSlot->numStepPerSeqVar = tempNumPer;
  dsSlot->numStepPerSeqVar[dsSlot->numSeqVar] = seqSlot->numStep;
  dsSlot->seqVarIDs = tempSlotIDs;
  for (i = 0; i < seqSlot->numStep; i++) {
    dsSlot->seqVarIDs[dsSlot->numSeqVar][i] = (*seqVars)[i].slotID;
  }
  dsSlot->numSeqVar++;

  return FC_SUCCESS;
}

/**
 * \ingroup Variable
 * \brief  Create a new variable.
 *
 * \description
 *
 *    Create a new, in-memory, variable.
 *
 * \modifications  
 *    - Nancy Collins, Created.
 */
FC_ReturnCode fc_createVariable(
  FC_Mesh mesh,          /**< input - mesh handle */
  char *varName,         /**< input - variable name */
  FC_Variable *variable  /**< output - variable handle */
) {
  int *tempslotIDs;
  _FC_MeshSlot* meshSlot;
  _FC_VarSlot* varSlot;
  
  // default return
  if (variable)
    *variable = FC_NULL_VARIABLE;
  
  // check input
  meshSlot = _fc_getMeshSlot(mesh);
  if (meshSlot == NULL || !varName || !variable) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Creating new variable '%s'", varName);
  
  // get an open slot
  varSlot = _fc_getNewVarSlot();
  if (varSlot == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  
  // table header information
  _FC_GET_HANDLE(*variable, varSlot);
  _fc_setSlotHeaderName(&varSlot->header, varName);
  
  // set back and forth references
  varSlot->ds = meshSlot->ds;
  varSlot->mesh = mesh;
  tempslotIDs = realloc(meshSlot->basicVarIDs, 
                        sizeof(int) * (meshSlot->numBasicVar+1));
  if (tempslotIDs == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  meshSlot->basicVarIDs = tempslotIDs;
  meshSlot->basicVarIDs[meshSlot->numBasicVar] = varSlot->header.slotID;
  meshSlot->numBasicVar++;
  
  return FC_SUCCESS;
}

/**
 * \ingroup  Variable 
 * \brief  Create a new sequence variable.
 *
 * \description
 *  
 *    Creates a new, in-memory, sequence variable as an array
 *    of variables. A sequence must be provided in addition
 *    to the mesh that will be the owner of the variable.
 *    The sequence and mesh must be in the same dataset.
 *
 * \modifications  
 *    - 2003-NOV-18 WSK Created.
 */
FC_ReturnCode fc_createSeqVariable(
  FC_Mesh mesh,          /**< input - mesh handle */
  FC_Sequence sequence,  /**< input - sequence handle */
  char *varName,         /**< input - sequence variable name */
  int* numStep,          /**< output - the number of steps */
  FC_Variable **seqVars  /**< output - seq variable handles */
) {
  int i;
  _FC_MeshSlot* meshSlot;
  _FC_SeqSlot* seqSlot;
  _FC_VarSlot* varSlot;
  int *tempNumPer, **tempSlotIDs;
  
  // default output
  if (numStep)
    *numStep = -1;
  if (seqVars)
    *seqVars = NULL;
  
  // check input
  meshSlot = _fc_getMeshSlot(mesh);
  seqSlot = _fc_getSeqSlot(sequence);
  if (meshSlot == NULL || seqSlot == NULL || varName == NULL || 
      numStep == NULL || seqVars == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  if (!FC_HANDLE_EQUIV(meshSlot->ds, seqSlot->ds)) {
    fc_printfErrorMessage("Mesh and sequence need to be in same dataset");
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Creating new seq variable '%s'", varName);
  
  // make room for array & setup return values
  *numStep = seqSlot->numStep;
  *seqVars = malloc(sizeof(FC_Variable)*seqSlot->numStep);
  if (*seqVars == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  
  for (i = 0; i < seqSlot->numStep; i++) {
    // get an open slot
    varSlot = _fc_getNewVarSlot();
    if (varSlot == NULL) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));      
      return FC_MEMORY_ERROR;
    }
    
    // table header information
    _FC_GET_HANDLE((*seqVars)[i], varSlot);
    _fc_setSlotHeaderName(&varSlot->header, varName);
    
    // set back and forth references
    varSlot->ds = meshSlot->ds;
    varSlot->mesh = mesh;
    varSlot->sequence = sequence;
    varSlot->stepID = i;
  }
  
  // set more back and forth references
  tempNumPer = (int*)realloc(meshSlot->numStepPerSeqVar, 
			     (meshSlot->numSeqVar+1)*sizeof(int));
  tempSlotIDs = (int**)realloc(meshSlot->seqVarIDs,
			       (meshSlot->numSeqVar+1)*sizeof(int*));
  if (!tempNumPer || !tempSlotIDs) {
    free(tempNumPer); free(tempSlotIDs);
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  tempSlotIDs[meshSlot->numSeqVar] = (int*)malloc(seqSlot->numStep*sizeof(int)); 
  if (!tempSlotIDs[meshSlot->numSeqVar]) {
    free(tempNumPer); free(tempSlotIDs);
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  meshSlot->numStepPerSeqVar = tempNumPer;
  meshSlot->numStepPerSeqVar[meshSlot->numSeqVar] = seqSlot->numStep;
  meshSlot->seqVarIDs = tempSlotIDs;
  for (i = 0; i < seqSlot->numStep; i++) {
    meshSlot->seqVarIDs[meshSlot->numSeqVar][i] = (*seqVars)[i].slotID;
  }
  meshSlot->numSeqVar++;
  
  return FC_SUCCESS;
}

/**
 * \ingroup  Variable
 * \brief  Set the variable data by copy.
 *
 * \description 
 * 
 *    Sets the data values associated with this new
 *    variable. Must provide the number of data points, the number of
 *    components per data point, the math type of the data, and the
 *    data type of the array. (For example, if you have a 3D vector field
 *    associated with the elements of a 10 element mesh, and it'd data
 *    type is double, the call would be: 
 *
 *    fc_setVariableData(var, 10, 3, FC_AT_ELEMENT, FC_MT_VECTOR,
 *    FC_DT_DOUBLE, data);
 *
 *    The data array should be numDataPoint*numComponent long. The order is
 *    interleave, that is you give all the components for each point
 *    (e.g. P1C1 P1C2 ... P1Cn P2C1 P2C2 ... P1Cn ... PmC1 ... PmCn).
 *
 *    This routine makes a copy of the data for the library,
 *    so the user still has responsibility of the original buffer.
 *    If you would like to give the data array to the variable instead
 *    of just copying them, use fc_setVariableDataPtr().
 *
 * \modifications  
 *    - Nancy Collins, Created.
 *    - 02/03/04 WSK, Now handles edge & face associated data.
 *    - 01/30/06 CDU, Changed logging to support global variables
 */
FC_ReturnCode fc_setVariableData(
  FC_Variable variable,      /**< input - variable handle */
  int numDataPoint,          /**< input - the number of data points */
  int numComponent,          /**< input - the number of components */
  FC_AssociationType assoc,  /**< input - data association (e.g. vertices, 
                                            elements, etc) */
  FC_MathType mathtype,      /**< input - the organization of components (e.g.
                                            scalar, vector, etc.) */
  FC_DataType datatype,      /**< input - data type (e.g. int, float, etc.) */
  void *data                 /**< input - data buffer */
) {
  FC_ReturnCode rc;
  size_t size;
  _FC_VarSlot* varSlot;
  int temp_numDataPoint;
  
  // NOTE!: if you change stuff here, also change fc_setVariableDataPtr()

  // check input
  varSlot = _fc_getVarSlot(variable);
  if (varSlot == NULL || numDataPoint < 1 || numComponent < 1 || 
      !fc_isAssociationTypeValid(assoc) || assoc == FC_AT_UNKNOWN || 
      !fc_isMathTypeValid(mathtype) || mathtype == FC_MT_UNKNOWN ||
      !fc_isDataTypeValid(datatype) || datatype == FC_DT_UNKNOWN || 
      data == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // a little more checking --
  // a global variable can only have an assoc of FC_AT_WHOLE_DATASET
  if (fc_isVariableGlobal(variable) && assoc != FC_AT_WHOLE_DATASET) {
    fc_printfErrorMessage("Trying to write non-global data to global var '%s'",
			  varSlot->header.name);
    return FC_INPUT_ERROR;
  }
  else if (!fc_isVariableGlobal(variable) && assoc == FC_AT_WHOLE_DATASET) {
    fc_printfErrorMessage("Trying to write global data to non-global var '%s'",
			  varSlot->header.name);
  }
  // can't overwrite data
  if (varSlot->data != NULL || varSlot->header.committed == 1) {
    fc_printfErrorMessage("Data already exists on variable '%s'",
                          varSlot->header.name);
    return FC_ERROR;
  }
  // make sure mathtype & numComp are consistent
  if ( (numComponent == 1 && mathtype != FC_MT_SCALAR) ||
       (numComponent > 1  && mathtype == FC_MT_SCALAR) ) {
    fc_printfErrorMessage("Inconsistent number of components and math type "
                          "for variable '%s'", varSlot->header.name);
    return FC_INPUT_ERROR;
  }
  // Check that numDataPoint is consistent
  // (this will also force build edges & faces if necessary)
  if (fc_isVariableGlobal(variable)) {
    if (numDataPoint != 1) {
      fc_printfErrorMessage("Global var '%s' can only have 1 datapoint",
			    varSlot->header.name);
      return FC_INPUT_ERROR;
    }
  }
  else {
    rc = fc_getMeshNumEntity(varSlot->mesh, assoc, &temp_numDataPoint);
    if (rc != FC_SUCCESS)
      return rc;
    if (temp_numDataPoint != numDataPoint) {
      fc_printfErrorMessage("Inconsistent number of data points and association "
                          "type for variable '%s'", varSlot->header.name);
      return FC_INPUT_ERROR;
    }
  }
  
  // log message
  if(!fc_isVariableGlobal(variable)){
    //For Non-global variables, get the mesh's name
    _FC_MeshSlot* meshSlot;
    meshSlot = _fc_getMeshSlot(varSlot->mesh);
    fc_printfLogMessage("Setting vertices on mesh '%s'", meshSlot->header.name);
  } else {
    //For global variables, just get the dataset's name
    _FC_DsSlot *datasetSlot;
    datasetSlot = _fc_getDsSlot(varSlot->ds);
    fc_printfLogMessage("Setting vertices for dataset '%s'",datasetSlot->header.name);
  }

  // copy metadata
  varSlot->numDataPoint = numDataPoint;
  varSlot->numComponent = numComponent;
  varSlot->assoc = assoc;
  varSlot->mathtype = mathtype;
  varSlot->datatype = datatype;

  // copy data buffer
  size = numDataPoint * numComponent * fc_sizeofDataType(datatype);
  varSlot->data = malloc(size);
  if (varSlot->data == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  memcpy(varSlot->data, data, size);
  
  return FC_SUCCESS;
}

/**
 * \ingroup  Variable
 * \brief  Set the pointer to the variable data.
 *
 * \description
 *  
 *    Sets the data values associated with this new
 *    variable. Must provide the number of data points, the number of
 *    components per data point, the math type of the data, and the
 *    data type of the array. (For example, if you have a 3D vector field
 *    associated with the elements of a 10 element mesh, and it'd data
 *    type is double, the call would be: 
 *
 *    fc_setVariableData(var, 10, 3, FC_AT_ELEMENT, FC_MT_VECTOR,
 *    FC_DT_DOUBLE, data);
 *
 *    The data array should be numDataPoint*numComponent long. The order is
 *    interleave, that is you give all the components for each point
 *    (e.g. P1C1 P1C2 ... P1Cn P2C1 P2C2 ... P1Cn ... PmC1 ... PmCn).
 *
 *    This routine passes the data array directly to the library,
 *    so the user no longer has responsibility of the original buffer.
 *    If you would like to keep the data array, you should use 
 *    fc_setVariableData() instead which makes a copy of the data array.
 *
 *
 * \modifications  
 *    - 10/01/04 WSK, created.
 *    - 01/30/06 CDU, changed logging to support global variables
 */
FC_ReturnCode fc_setVariableDataPtr(
  FC_Variable variable,      /**< input - variable handle */
  int numDataPoint,          /**< input - the number of data points */
  int numComponent,          /**< input - the number of components */
  FC_AssociationType assoc,  /**< input - data association (e.g. vertices, 
                                            elements, etc) */
  FC_MathType mathtype,      /**< input - the organization of components (e.g.
                                            scalar, vector, etc.) */
  FC_DataType datatype,      /**< input - data type (e.g. int, float, etc.) */
  void *data                 /**< input - data buffer */
) {
  FC_ReturnCode rc;
  _FC_VarSlot* varSlot;
  int temp_numDataPoint;
  
  // NOTE!: if you change stuff here, also change fc_setVariableData()

  // check input
  varSlot = _fc_getVarSlot(variable);
  if (varSlot == NULL || numDataPoint < 1 || numComponent < 1 || 
      !fc_isAssociationTypeValid(assoc) || assoc == FC_AT_UNKNOWN || 
      !fc_isMathTypeValid(mathtype) || mathtype == FC_MT_UNKNOWN ||
      !fc_isDataTypeValid(datatype) || datatype == FC_DT_UNKNOWN || 
      data == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // a little more checking -- can't overwrite data
  // a global variable can only have an assoc of FC_AT_WHOLE_DATASET
  if (fc_isVariableGlobal(variable) && assoc != FC_AT_WHOLE_DATASET) {
    fc_printfErrorMessage("Trying to write non-global data to global var '%s'",
			  varSlot->header.name);
    return FC_INPUT_ERROR;
  }
  else if (!fc_isVariableGlobal(variable) && assoc == FC_AT_WHOLE_DATASET) {
    fc_printfErrorMessage("Trying to write global data to non-global var '%s'",
			  varSlot->header.name);
  }
  // can't overwrite data
  if (varSlot->data != NULL || varSlot->header.committed == 1) {
    fc_printfErrorMessage("Data already exists on variable '%s'",
                          varSlot->header.name);
    return FC_ERROR;
  }
  // make sure mathtype & numComp are consistent
  if ( (numComponent == 1 && mathtype != FC_MT_SCALAR) ||
       (numComponent > 1  && mathtype == FC_MT_SCALAR) ) {
    fc_printfErrorMessage("Inconsistent number of components and math type "
                          "for variable '%s'", varSlot->header.name);
    return FC_INPUT_ERROR;
  }
  // Check that numDataPoint is consistent
  // (this will also force build edges & faces if necessary)
  if (fc_isVariableGlobal(variable)) {
    if (numDataPoint != 1) {
      fc_printfErrorMessage("Global var '%s' can only have 1 datapoint",
			    varSlot->header.name);
      return FC_INPUT_ERROR;
    }
  }
  else {
    rc = fc_getMeshNumEntity(varSlot->mesh, assoc, &temp_numDataPoint);
    if (rc != FC_SUCCESS)
      return rc;
    if (temp_numDataPoint != numDataPoint) {
      fc_printfErrorMessage("Inconsistent number of data points and association "
			    "type for variable '%s'", varSlot->header.name);
      return FC_INPUT_ERROR;
    }
  }

  // log message
  if(!fc_isVariableGlobal(variable)){
    //Non-global variable
    _FC_MeshSlot* meshSlot;
    meshSlot = _fc_getMeshSlot(varSlot->mesh);
    fc_printfLogMessage("Setting vertices on mesh '%s'", meshSlot->header.name);
  } else {
    //Global variables don't have a mesh. Use the dataset name
    _FC_DsSlot *datasetSlot;
    datasetSlot = _fc_getDsSlot(varSlot->ds);
    fc_printfLogMessage("Setting vertices for dataset '%s'",datasetSlot->header.name);
  }

  // copy metadata
  varSlot->numDataPoint = numDataPoint;
  varSlot->numComponent = numComponent;
  varSlot->assoc = assoc;
  varSlot->mathtype = mathtype;
  varSlot->datatype = datatype;

  // copy data buffer
  varSlot->data = data;
  
  return FC_SUCCESS;
}

//@}

/** \name Other ways to get new variables. */
//-------------------------------
//@{

/**
 * \ingroup  Variable
 * \brief  Copy a global variable
 *
 * \description
 *  
 *    Create a copy of a global variable into the specified dataset using 
 *    the  given name. 
 *
 *    The input variable can be a non-global associated with a whole
 *    mesh (FC_AT_WHOLE_MESH).
 *
 * \modifications 
 *   - 2006/09/19 WSD Created.
 *   - 11/08/2006 WSD added ability to whole mesh var onto a dataset.
 */
FC_ReturnCode fc_copyGlobalVariable(
  FC_Variable src_var, /**< Input - the source variable to be copied */
  FC_Dataset dest_dataset, /**< Input - the destination dataset to be copied into */
  char *newName,       /**< Input - the name of the new variable, a NULL value
                            is a flag to use the name of the source variable */
  FC_Variable* new_var /**< Output - handle to the new variable */
  ) { 
  FC_ReturnCode rc;
  _FC_VarSlot* varSlot;
  int numDataPoint, numComp;
  FC_AssociationType assoc;
  FC_MathType mathtype;
  FC_DataType datatype;
  void* data;
  
  // default return value
  if (new_var)
    *new_var = FC_NULL_VARIABLE;

  // check input - NULL name is o.k.
  varSlot = _fc_getVarSlot(src_var);
  if (varSlot == NULL || !fc_isDatasetValid(dest_dataset) || new_var == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // test for NULL name
  if (newName == NULL)
    newName = varSlot->header.name;

  // log message
  fc_printfLogMessage("Copying global variable '%s' to '%s'", 
                      varSlot->header.name, newName);

  // get other info from source variable &
  rc = fc_getVariableInfo(src_var, &numDataPoint, &numComp, &assoc,
                          &mathtype, &datatype);
  if (rc != FC_SUCCESS) {
    fc_printfErrorMessage("Failed to get info for variable '%s'",
                          varSlot->header.name);
    return rc;
  }

  // more checking
  if (assoc == FC_AT_WHOLE_MESH || assoc == FC_AT_WHOLE_DATASET) {
    assoc = FC_AT_WHOLE_DATASET;
  }
  else {
    fc_printfErrorMessage("var of assoc type %s cannot be copied to a dataset",
                          fc_getAssociationTypeText(assoc));
    return FC_INPUT_ERROR;
  }

  // --- do it
  rc = fc_getVariableDataPtr(src_var, &data);
  if (rc != FC_SUCCESS) {
    fc_printfErrorMessage("Failed to get data for variable '%s'",
                          varSlot->header.name);
    return rc;
  }
  rc = fc_createGlobalVariable(dest_dataset, newName, new_var);
  if (rc != FC_SUCCESS) {
    fc_printfErrorMessage("Failed to create new global variable '%s'", newName);
    return rc;
  }
  rc = fc_setVariableData(*new_var, numDataPoint, numComp, assoc, mathtype, 
                   datatype, data);
  if (rc != FC_SUCCESS) {
    fc_printfErrorMessage("Failed to set data for variable '%s'", newName);
    return rc;
  }

  return rc;
}

/**
 * \ingroup  Variable
 * \brief  Copy a global sequence variable
 *
 * \description
 *  
 *    Create a copy of a global seq variable into the specified dataset using 
 *    the  given name. This will return an error if the source and
 *    destination sequences are not similar enough. This will copy a sequence
 *    if necessary. The destination sequence's dataset needs to be
 *    the same as the destintation dataset or an error will be returned.
 *
 *    The input variable can be a non-global associated with a whole
 *    mesh (FC_AT_WHOLE_MESH).
 *
 * \modifications 
 *   - 9/21/2006 WSD Created.
 *   - 11/08/2006 WSD added ability to whole mesh var onto a dataset.
 */
FC_ReturnCode fc_copyGlobalSeqVariable(
  int numStep, /**< Input - the number of steps in the sequence variable */
  FC_Variable* src_vars, /**< Input - the source seq variables to be copied */
  FC_Dataset dest_ds,   /**< Input - the destination dataset to be copied into */
  FC_Sequence dest_seq, /**< Input - the destination seq to be copied into */
  char *newName,       /**< Input - the name of the new variable, a NULL value
                            is a flag to use the name of the source variable */
  FC_Variable** new_vars /**< Output - array of handles to the new variable */
) {
  int i; 
  FC_ReturnCode rc;
  int temp_numStep;
  _FC_DsSlot* dest_dsSlot;
  _FC_SeqSlot* dest_seqSlot;
  _FC_VarSlot* varSlot;
  int numDataPoint, numComp;
  FC_AssociationType assoc;
  FC_MathType mathtype;
  FC_DataType datatype;
  void* data;
  
  // default return value
  if (new_vars)
    *new_vars = NULL;
  
  // check input
  dest_dsSlot = _fc_getDsSlot(dest_ds);
  dest_seqSlot = _fc_getSeqSlot(dest_seq);
  if (!fc_isSeqVariableValid(numStep, src_vars) || dest_dsSlot == NULL || 
      dest_seqSlot == NULL || !FC_HANDLE_EQUIV(dest_ds, dest_seqSlot->ds) ||
      numStep != dest_seqSlot->numStep  || new_vars == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  varSlot = _fc_getVarSlot(src_vars[0]);

  // test for NULL name
  if (newName == NULL)
    newName = varSlot->header.name;

  // log message
  fc_printfLogMessage("Copying global sequence variable '%s' to '%s'", 
                      varSlot->header.name, newName);

  // get src meta data
  rc = fc_getVariableInfo(src_vars[0], &numDataPoint, &numComp, &assoc,
			  &mathtype, &datatype);
  if (rc != FC_SUCCESS) {
    fc_printfErrorMessage("Failed to get info for seq variable '%s'",
                          varSlot->header.name);
    return rc;
  }

  // more checking
  if (assoc == FC_AT_WHOLE_MESH || assoc == FC_AT_WHOLE_DATASET) {
    assoc = FC_AT_WHOLE_DATASET;
  }
  else {
    fc_printfErrorMessage("var of assoc type %s cannot be copied to a dataset",
                          fc_getAssociationTypeText(assoc));
    return FC_INPUT_ERROR;
  }

  // create new variable
  rc = fc_createGlobalSeqVariable(dest_ds, dest_seq, newName, &temp_numStep,
				  new_vars);
  if (rc != FC_SUCCESS) {
    fc_printfErrorMessage("Failed to create new seq variable '%s'", newName);
    return rc;
  }
  
  // get and copy data for each step
  for (i = 0; i < numStep; i++) {
    rc = fc_getVariableDataPtr(src_vars[i], &data);
    if (rc != FC_SUCCESS) {
      fc_printfErrorMessage("Failed to get data for step %d of seq variable "
                            "'%s'", i, varSlot->header.name);
      return rc;
    }
    
    rc = fc_setVariableData((*new_vars)[i], numDataPoint, numComp, assoc, 
                            mathtype, datatype, data);
    if (rc != FC_SUCCESS) {
      fc_printfErrorMessage("Failed to set data for step %d of seq variable "
                            "'%s'", i, newName);
      return rc;
    }
  }
  
  return rc;
}

/**
 * \ingroup  Variable
 * \brief  Copy a variable
 *
 * \description
 *  
 *    Create a copy of a variable into the specified mesh using with the 
 *    given name. This will return an error if the source and
 *    destination meshes are not similar enough. 
 *
 *    The input variable can be a global variable; the output will be
 *    a variable associated with the mesh (FC_AT_WHOLE_MESH).
 *
 * \modifications 
 *   - 2003-NOV-11 WSK Created
 *   - 06/07/04 RM, added checking if destination mesh and a source variable
 *          are compatible to perform copy before creating the variable.
 *   - 11/08/2006 WSD added ability to copy global var onto a mesh.
 */
FC_ReturnCode fc_copyVariable(
  FC_Variable src_var, /**< Input - the source variable to be copied */
  FC_Mesh dest_mesh,   /**< Input - the destination mesh to be copied into */
  char *newName,       /**< Input - the name of the new variable, a NULL value
                            is a flag to use the name of the source variable */
  FC_Variable* new_var /**< Output - handle to the new variable */
  ) { 
  FC_ReturnCode rc;
  _FC_VarSlot* varSlot;
  int numDataPoint, numComp, numEntity;
  FC_AssociationType assoc;
  FC_MathType mathtype;
  FC_DataType datatype;
  void* data;
  
  // default return value
  if (new_var)
    *new_var = FC_NULL_VARIABLE;

  // check input - NULL name is o.k.
  varSlot = _fc_getVarSlot(src_var);
  if (varSlot == NULL || !fc_isMeshValid(dest_mesh) || new_var == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // test for NULL name
  if (newName == NULL)
    newName = varSlot->header.name;

  // log message
  fc_printfLogMessage("Copying variable '%s' to '%s'", 
                      varSlot->header.name, newName);

  // Do nothing special for a coord variable -WSK 9/15/04

  // get other info from source variable & dest mesh
  rc = fc_getVariableInfo(src_var, &numDataPoint, &numComp, &assoc,
                          &mathtype, &datatype);
  if (rc != FC_SUCCESS) {
    fc_printfErrorMessage("Failed to get info for variable '%s'",
                          varSlot->header.name);
    return rc;
  }
  // more checking - destination mesh and variable assoc type are compatible
  if (assoc == FC_AT_WHOLE_DATASET || assoc == FC_AT_WHOLE_MESH) {
    assoc = FC_AT_WHOLE_MESH;
  }
  else {
    rc = fc_getMeshNumEntity(dest_mesh, assoc, &numEntity); 
    if (rc != FC_SUCCESS) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
      return FC_INPUT_ERROR;
    }
    if (numEntity != numDataPoint) {
      fc_printfErrorMessage("%s: destination mesh not compatible with src variable", 
                            fc_getReturnCodeText(FC_INPUT_ERROR));
      return FC_INPUT_ERROR;
    }
  }

  // --- do it
  rc = fc_getVariableDataPtr(src_var, &data);
  if (rc != FC_SUCCESS) {
    fc_printfErrorMessage("Failed to get data for variable '%s'",
                          varSlot->header.name);
    return rc;
  }
  rc = fc_createVariable(dest_mesh, newName, new_var);
  if (rc != FC_SUCCESS) {
    fc_printfErrorMessage("Failed to create new variable '%s'", newName);
    return rc;
  }
  rc = fc_setVariableData(*new_var, numDataPoint, numComp, assoc, mathtype, 
                   datatype, data);
  if (rc != FC_SUCCESS) {
    fc_printfErrorMessage("Failed to set data for variable '%s'", newName);
    return rc;
  }


  return rc;
}


/**
 * \ingroup  Variable
 * \brief  Copy a variable defined on a smaller mesh onto a bigger one
 *
 * \description
 *  
 *    Copy a variable defined on a smaller mesh onto a bigger one.
 *    Meshes can be on different datasets.
 *    This is intended for cases where you may promote a subset to
 *    a mesh (becuase, for instance, it corresponds to a ROI) 
 *    and later want to map that area back onto the original mesh.
 *
 *    newName is the name of the var on the dest mesh.
 *    If the variable does not exist on the dest mesh, a new
 *    var is created there and values for the entities on the
 *    dest mesh that do not exist on the original mesh are
 *    filled in with the fillval. If the variable does exist on
 *    the dest mesh, then the vals from the srv_var are copied over
 *    according to the mapping. This latter action allows one to
 *    call this function in a loop, adding in the same var
 *    from several regions sequentially. (Note this is unlike
 *    copyVarToRegionMesh where the dest var cannot preexist.)
 *
 *    The entity mapping takes the form of another variable of the
 *    same association that contains the ID numbers of the corresponding
 *    entity type in the larger mesh. It is assumed that this is created
 *    at the time of promoting the subset to the mesh. Note that this
 *    is different that the lookup table ids since these go from
 *    currmesh to prev mesh, as opposed to the original on disk 
 *    numbering to the new fclib numbering.
 *
 *    Note this doest not check the coords or the conns to see if they
 *    are the same on the 2 meshes, it only goes by the mapping numbers.
 *
 *    Note that since we can only write out variables as atrributes in 
 *    exodus and therefore as doubles, the mapping var may not be explictly
 *    FC_DT_INT. We may want to explictly convert this on read in. For now
 *    we will just cast/round the val to an int if it is a float or double.
 *
 *    Note: if there is an invalid mapping data pt, it just wont copy that
 *    one over. should change this so the call fails. 
 *
 * \todo
 *   - fix the rounding. couldnt get round or rint to work.
 *   - prob should put the conversion for the region variable
 *     when it is read in from exodus.
 *   - have it fail if there is an invalid mapping data pt
 *
 * \modifications 
 *   - 04/24/08 ACG Created
 */
FC_ReturnCode fc_copyVariableFromRegionMesh(
  FC_Variable src_var, /**< Input - the source variable to be copied */
  FC_Mesh dest_mesh,   /**< Input - the destination mesh to be copied into */
  FC_Variable mapping, /**< Input - int var that maps entities from this mesh to the dest mesh */
  void* fillval,       /**< Input - optional variable fill value for entities not in the map */
  char *newName,       /**< Input - the name of the new variable, a NULL value
                            is a flag to use the name of the source variable */
  FC_Variable* new_var /**< Output - handle to the new variable */
) { 
  FC_ReturnCode rc;
  _FC_VarSlot* varSlot;
  FC_Variable* existVars;
  int numExist;
  int numDataPoint,currnumDataPoint, numComponent, currnumComponent,numEntity;
  int mappingnumdp,mappingnumc;
  FC_AssociationType assoc, currassoc, mappingassoc;
  FC_MathType mathtype, currmathtype, mappingmathtype;
  FC_DataType datatype, currdatatype, mappingdatatype;
  void *data, *new_data, *currdata, *mappingdata;
  int i, j;

  
  // default return value
  if (new_var)
    *new_var = FC_NULL_VARIABLE;

  // check input - NULL name is o.k.
  varSlot = _fc_getVarSlot(src_var);
  if (varSlot == NULL || !fc_isVariableValid(mapping) || !fc_isMeshValid(dest_mesh)
      || new_var == NULL){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // test for NULL name
  if (newName == NULL)
    newName = varSlot->header.name;

  // log message
  fc_printfLogMessage("Copying region variable '%s' to '%s'", 
                      varSlot->header.name, newName);

  //more checking
  rc = fc_getVariableInfo(mapping, &mappingnumdp, &mappingnumc,
			  &mappingassoc,&mappingmathtype,&mappingdatatype);
  if (rc != FC_SUCCESS) {
    fc_printfErrorMessage("Failed to get info for mapping variable");
    return rc;
  }
  rc = fc_getVariableInfo(src_var, &numDataPoint, &numComponent,
			  &assoc,&mathtype,&datatype);
  if (rc != FC_SUCCESS) {
    fc_printfErrorMessage("Failed to get info for src variable");
    return rc;
  }

  //NOTE: from writeout - this may be a double
  //  if (mappingdatatype != FC_DT_INT){
  //    fc_printfErrorMessage("mapping variable must be FC_DT_INT");
  //    return FC_INPUT_ERROR;
  //  }
  if (mappingdatatype != FC_DT_INT && mappingdatatype != FC_DT_DOUBLE &&
      mappingdatatype != FC_DT_FLOAT){
    fc_printfErrorMessage("mappingdatatype must be numeric");
    return FC_INPUT_ERROR;
  }
  if (mappingnumc != 1){
    fc_printfErrorMessage("mapping variable must be single component");
    return FC_INPUT_ERROR;
  }
  if (mappingassoc != assoc){
    fc_printfErrorMessage("mapping variable and src variable must have same assoc");
    return FC_INPUT_ERROR;
  }
  if (datatype != FC_DT_INT && datatype != FC_DT_DOUBLE && datatype != FC_DT_FLOAT &&
      datatype != FC_DT_CHAR){
    fc_printfErrorMessage("src variable must have known type");
    return FC_INPUT_ERROR;
  }
    
  //temporary
  if (assoc == FC_AT_WHOLE_DATASET || assoc == FC_AT_WHOLE_MESH) {
    fc_printfErrorMessage("not handling globals at this time");
    return FC_INPUT_ERROR;
  }
  if (assoc != FC_AT_VERTEX && assoc != FC_AT_ELEMENT){
    fc_printfErrorMessage("only handling vertex and element vars at this time");
    return FC_INPUT_ERROR;
  }
  //check appropriateness of meshes - are there enough of this entity in
  //the dest mesh (should we check dimension etc?)
  rc = fc_getMeshNumEntity(dest_mesh, assoc, &numEntity); 
  if (rc != FC_SUCCESS) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  if (numEntity < numDataPoint) {
    fc_printfErrorMessage("%s: destination mesh not compatible with src variable", 
			  fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  //does the var already exist?
  rc = fc_getVariableByName(dest_mesh, newName, &numExist, &existVars);
  if (rc != FC_SUCCESS){
    fc_printfErrorMessage("Failed to get info for curr variable on dest mesh");
    return rc;
  }    
  if (numExist > 1){
    fc_printfErrorMessage("Cannot copy region var onto new mesh since there is more than 1 instance there");
    free(existVars);
    return FC_ERROR;
  }
  if (numExist == 1){
    //check to see if it is compatable
    rc = fc_getVariableInfo(existVars[0], &currnumDataPoint, &currnumComponent,
			  &currassoc,&currmathtype,&currdatatype);
    
    if (rc != FC_SUCCESS) {
      fc_printfErrorMessage("Failed to get info for curr variable on dest mesh");
      free(existVars);
      return rc;
    }
    if (currnumComponent != numComponent){
      fc_printfErrorMessage("Cannot copy region var onto mesh because of existing numComponent incompatability");
      free(existVars);
      return FC_INPUT_ERROR;
    }
    if (currassoc != assoc){
      fc_printfErrorMessage("Cannot copy region var onto mesh because of existing assoc incompatability");
      free(existVars);
      return FC_INPUT_ERROR;
    }
    if (currmathtype != mathtype){
      fc_printfErrorMessage("Cannot copy region var onto mesh because of existing mathtype incompatability");
      free(existVars);
      return FC_INPUT_ERROR;
    }
    if (currdatatype != datatype){
      fc_printfErrorMessage("Cannot copy region var onto mesh because of existing datatype incompatability");
      free(existVars);
      return FC_INPUT_ERROR;
    }
  }

  // --- do it
  rc = fc_getVariableDataPtr(src_var, &data);
  if (rc != FC_SUCCESS) {
    fc_printfErrorMessage("Failed to get data for variable '%s'",
                          varSlot->header.name);
    return rc;
  }
  rc = fc_getVariableDataPtr(mapping, &mappingdata);
  if (rc != FC_SUCCESS) {
    fc_printfErrorMessage("Failed to get data for mapping variable");
    free(data);
    return rc;
  }

  if (numExist == 0){
    // create the new data
    switch(datatype){
    case FC_DT_INT:
      new_data = calloc(numEntity*numComponent,sizeof(int));
      if (new_data == NULL){
	fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
	return FC_MEMORY_ERROR;
      }
      if (fillval != NULL){
	for (i = 0; i < numEntity*numComponent; i++){
	  ((int*)new_data)[i] = *((int*)fillval);
	}
      }
      break;
    case FC_DT_DOUBLE:
      new_data = calloc(numEntity*numComponent,sizeof(double));
      if (new_data == NULL){
	fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
	return FC_MEMORY_ERROR;
      }
      if (fillval != NULL){
	for (i = 0; i < numEntity*numComponent; i++){
	  ((double*)new_data)[i] = *((double*)fillval);
	}
      }
      break;
    case FC_DT_FLOAT:
      new_data = calloc(numEntity*numComponent,sizeof(float));
      if (new_data == NULL){
	fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
	return FC_MEMORY_ERROR;
      }
      if (fillval != NULL){
	for (i = 0; i < numEntity*numComponent; i++){
	  ((float*)new_data)[i] = *((float*)fillval);
	}
      }
      break;
    case FC_DT_CHAR:
      new_data = calloc(numEntity*numComponent,sizeof(char));
      if (new_data == NULL){
	fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
	return FC_MEMORY_ERROR;
      }
      if (fillval != NULL){
	for (i = 0; i < numEntity*numComponent; i++){
	  ((char*)new_data)[i] = *((char*)fillval);
	}
      }
      break;
    default:
      //shouldnt happen
      fc_printfErrorMessage("Cannot handle datatype %d", datatype);
      return FC_ERROR;
    }
  } else {
    //will write it into the current data
    rc = fc_getVariableDataPtr(existVars[0], &currdata);
    if (rc != FC_SUCCESS) {
      fc_printfErrorMessage("Failed to get data for variable '%s'",
			    varSlot->header.name);
      free(existVars);
      return rc;
    }
    new_data = currdata;
  }

  for (i = 0; i < numDataPoint; i++){
    int meshpt;
    switch(mappingdatatype){
    case FC_DT_DOUBLE:
      meshpt = (int)floor((((double*) mappingdata)[i])+0.5);
      break;
    case FC_DT_FLOAT:
      meshpt = (int)floor((((float*) mappingdata)[i])+0.5);
      break;
    case FC_DT_INT:
      meshpt = ((int*) mappingdata)[i];
      break;
    default:
      //wont happen
      break;
    }
    //is it a valid pt?
    if (meshpt >= numEntity){
      fc_printfWarningMessage("Invalid mapping val.not copying this point");
    } else {
      switch (datatype){
      case FC_DT_INT:
	for (j = 0; j < numComponent; j++){
	  ((int*)new_data)[meshpt*numComponent+j] = ((int*)data)[i*numComponent+j];
	}
	break;
      case FC_DT_DOUBLE:
	for (j = 0; j < numComponent; j++){
	  //	printf("replacing big mesh pt %d with curr val %g with small mesh pt %d with val %g\n",
	  //	       (meshpt*numComponent+j), ((double*)new_data)[meshpt*numComponent+j],
	  //	       (i*numComponent+j), ((double*)data)[i*numComponent+j]);
	  ((double*)new_data)[meshpt*numComponent+j] = ((double*)data)[i*numComponent+j];
	}
	break;
      case FC_DT_FLOAT:
	for (j = 0; j < numComponent; j++){
	  ((float*)new_data)[meshpt*numComponent+j] = ((float*)data)[i*numComponent+j];
	}
	break;
      case FC_DT_CHAR:
	for (j = 0; j < numComponent; j++){
	  ((char*)new_data)[meshpt*numComponent+j] = ((char*)data)[i*numComponent+j];
	}
	break;
      default:
	//wont happen
	break;
      }
    }
  }

  if (numExist == 0){
    rc = fc_createVariable(dest_mesh, newName, new_var);
    if (rc != FC_SUCCESS) {
      fc_printfErrorMessage("Failed to create new variable '%s'", newName);
      return rc;
    }
    
    rc = fc_setVariableDataPtr(*new_var, numEntity, numComponent, assoc, mathtype, 
			    datatype, new_data);
    if (rc != FC_SUCCESS) {
      fc_printfErrorMessage("Failed to set data for variable '%s'", newName);
      return rc;
    }
  } else {
    //this was newly allocated memory
    *new_var = existVars[0];
    free(existVars);
  }

  return FC_SUCCESS;
}


/**
 * \ingroup  Variable
 * \brief  Copy a variable defined on a smaller mesh into a var on a bigger one
 *
 * \description
 *   This is like fc_copyVariableFromRegionMesh, but where you
 *   are copying into a specific step of a seq var.
 *
 *   Note; if there is an invlaid mapping meshpt, it just wont copy that val.
 *
 *   \todo
 *   - change so fails if invalid mapping
 *
 *  \modifications
 *   - 07/20/08 ACG - revamp from prev version that passed in the handle to the var
 */
FC_ReturnCode fc_copySeqVariableStepFromRegionMesh( FC_Variable src_var, /**< Input - the source variable to be copied */
						    FC_Mesh dest_mesh, /**< Input - the dest mesh */
						    FC_Sequence dest_seq, /**< Input - the dest seq */
						    int targetStep, /**< Input - which step in the seq to put the data */
						    FC_Variable mapping, /**< Input - mapping var */
						    void* fillval, /**< Input - optional default fill var only applicable if creating the var */
						    char* dest_seqvarName, /**< Input - name of new seq var on dest mesh. If NULL the name of the src_var is used */
						    FC_Variable** dest_seqvar /**< Output - the resulting seq var */
) { 
  FC_ReturnCode rc;
  _FC_VarSlot *varSlot_s;
  FC_Variable **seqVars, *dest_sv;
  char* meshname;
  char* dest_svName = NULL;
  int numSeqVar, *numStepPerSeqVar, numStep;
  int numDataPoint,currnumDataPoint, numComponent, currnumComponent;
  int mappingnumdp,mappingnumc;
  FC_AssociationType assoc, currassoc, mappingassoc;
  FC_MathType mathtype, currmathtype, mappingmathtype;
  FC_DataType datatype, currdatatype, mappingdatatype;
  void *data, *new_data, *mappingdata;
  int numEntity;
  int i, j;


  // default return
  if (dest_seqvar)
    *dest_seqvar = NULL;


  //test input
  varSlot_s = _fc_getVarSlot(src_var);
  if (varSlot_s == NULL || !fc_isVariableValid(mapping) || !fc_isMeshValid(dest_mesh) ||
      !fc_isSequenceValid(dest_seq) || targetStep < 0 || !fc_isVariableValid(mapping)){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  rc = fc_getSequenceNumStep(dest_seq, &numStep);
  if (rc != FC_SUCCESS){
    fc_printfErrorMessage("Cannot get seq num step");
    return rc;
  }

  if (numStep < targetStep){
    fc_printfErrorMessage("Target step is greater than num step in sequence");
    return FC_INPUT_ERROR;
  }

  // log message
  rc = fc_getMeshName(dest_mesh, &meshname);
  if (rc != FC_SUCCESS){
    fc_printfErrorMessage("Cannot get mesh name");
    return FC_ERROR;
  }
  if (dest_seqvarName == NULL){
    rc = fc_getVariableName(src_var, &dest_svName);
    if (rc != FC_SUCCESS){
      fc_printfErrorMessage("Cannot get var name");
      free(meshname);
      return FC_ERROR;
    }
  } else {
    dest_svName = dest_seqvarName;
  }

  
  fc_printfLogMessage("Copying region variable '%s' to '%s' of mesh '%s'", 
                      varSlot_s->header.name, meshname, dest_svName);
  free(meshname);

  //check consistency of src_var and mapping var
  rc = fc_getVariableInfo(mapping, &mappingnumdp, &mappingnumc,
			  &mappingassoc,&mappingmathtype,&mappingdatatype);
  if (rc != FC_SUCCESS) {
    fc_printfErrorMessage("Failed to get info for mapping variable");
    if (dest_seqvarName == NULL) free (dest_svName);
    return rc;
  }
  rc = fc_getVariableInfo(src_var, &numDataPoint, &numComponent,
			  &assoc,&mathtype,&datatype);
  if (rc != FC_SUCCESS) {
    fc_printfErrorMessage("Failed to get info for src variable");
    if (dest_seqvarName == NULL) free (dest_svName);
    return rc;
  }

    //NOTE: from writeout - this may be a double
  //  if (mappingdatatype != FC_DT_INT){
  //    fc_printfErrorMessage("mapping variable must be FC_DT_INT");
  //    return FC_INPUT_ERROR;
  //  }
  if (mappingdatatype != FC_DT_INT && mappingdatatype != FC_DT_DOUBLE &&
      mappingdatatype != FC_DT_FLOAT){
    fc_printfErrorMessage("mappingdatatype must be numeric");
    if (dest_seqvarName == NULL) free (dest_svName);
    return FC_INPUT_ERROR;
  }
  if (mappingnumc != 1){
    fc_printfErrorMessage("mapping variable must be single component");
    if (dest_seqvarName == NULL) free (dest_svName);
    return FC_INPUT_ERROR;
  }
  if (mappingassoc != assoc){
    fc_printfErrorMessage("mapping variable and src variable must have same assoc");
    if (dest_seqvarName == NULL) free (dest_svName);
    return FC_INPUT_ERROR;
  }
  if (datatype != FC_DT_INT && datatype != FC_DT_DOUBLE && datatype != FC_DT_FLOAT &&
      datatype != FC_DT_CHAR){
    fc_printfErrorMessage("src variable must have known type");
    if (dest_seqvarName == NULL) free (dest_svName);
    return FC_INPUT_ERROR;
  }
    

  switch (assoc){
  case FC_AT_ELEMENT:
    rc = fc_getMeshInfo(dest_mesh, NULL, NULL, NULL, &numEntity, NULL);
    break;
  case FC_AT_VERTEX:
    rc = fc_getMeshInfo(dest_mesh, NULL, NULL, &numEntity, NULL, NULL);
    break;
  default:
    //temporary
    if (assoc == FC_AT_WHOLE_DATASET || assoc == FC_AT_WHOLE_MESH) {
      fc_printfErrorMessage("not handling globals at this time");
    } else {
      fc_printfErrorMessage("only handling vertex and element vars at this time");
    }
    if (dest_seqvarName == NULL) free (dest_svName);
    return FC_ERROR; //shouldnt happen at this point
  }


  //check consistency of dest var
  //do we have a dest var?
  //only get it if it exists, dont do the component coalescing
  rc = fc_getSeqVariableByName(dest_mesh,dest_seqvarName, &numSeqVar, &numStepPerSeqVar,
				    &seqVars);
  if (rc != FC_SUCCESS){
    fc_printfErrorMessage("Cannot get seq var by name");
    if (dest_seqvarName == NULL) free (dest_svName);
    return rc;
  }
  if (numSeqVar > 1){
    fc_printfErrorMessage("Cannot find unique seq var to copy to");
    for (i = 0; i < numSeqVar; i++)
      free(seqVars[i]);
    free(seqVars);
    free(numStepPerSeqVar);
    if (dest_seqvarName == NULL) free (dest_svName);
    return FC_INPUT_ERROR;
  }
  if (numSeqVar == 0){
    int num;

    //create the seq var
    // printf("Creating main var <%s>\n",dest_svName);
    
    rc = fc_createSeqVariable(dest_mesh, dest_seq,
			     dest_svName, &num, &dest_sv);
    if (rc != FC_SUCCESS){
      fc_printfErrorMessage("Cannot make seq var");
      if (dest_seqvarName == NULL) free (dest_svName);
      return rc;
    }

    switch (datatype){
    case FC_DT_INT:
      new_data = calloc(numEntity*numComponent,sizeof(int));
     if (new_data == NULL){
       fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
       return FC_MEMORY_ERROR;
     }
     if (fillval != NULL){
       for (i = 0; i < numEntity*numComponent; i++){
	 ((int*)new_data)[i] = *((int*)fillval);
       }
     }
     break;
   case FC_DT_DOUBLE:
     new_data = calloc(numEntity*numComponent,sizeof(double));
     if (new_data == NULL){
       fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
       return FC_MEMORY_ERROR;
     }
     if (fillval != NULL){
       for (i = 0; i < numEntity*numComponent; i++){
	 ((double*)new_data)[i] = *((double*)fillval);
       }
     }
     break;
   case FC_DT_FLOAT:
     new_data = calloc(numEntity*numComponent,sizeof(float));
     if (new_data == NULL){
       fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
       return FC_MEMORY_ERROR;
     }
     if (fillval != NULL){
       for (i = 0; i < numEntity*numComponent; i++){
	 ((float*)new_data)[i] = *((float*) fillval);
       }
     }
     break;
   case FC_DT_CHAR:
     new_data = calloc(numEntity*numComponent,sizeof(char));
     if (new_data == NULL){
       fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
       return FC_MEMORY_ERROR;
     }
     if (fillval != NULL){
       for (i = 0; i < numEntity*numComponent; i++){
	 ((char*)new_data)[i] = *((char*) fillval);
       }
     }
     break;
   default:
     //shouldnt happen
     fc_printfErrorMessage("Cannot handle datatype %d", datatype);
     free(dest_sv); 
     free(new_data);
     return FC_ERROR;
    }
  
    for (i = 0; i < numStep; i++){
      rc = fc_setVariableData(dest_sv[i], numEntity, numComponent,
			      assoc, mathtype, datatype, new_data);
      if (rc != FC_SUCCESS){
	fc_printfErrorMessage("cannot set variable data");
	free(dest_sv);
	free(new_data);
	return rc;
      }
    }
    free(new_data);
    if (dest_seqvarName == NULL) free (dest_svName);
  } else {
    //we have the seq var
    FC_Sequence tempseq;

    rc = fc_getSequenceFromSeqVariable(numStep,seqVars[0], &tempseq);
    if (rc != FC_SUCCESS){
      fc_printfErrorMessage("Cannot get sequence from seq var");
      for (i = 0; i < numSeqVar; i++)
	free(seqVars[i]);
      free(seqVars);
      free(numStepPerSeqVar);
      if (dest_seqvarName == NULL) free (dest_svName);
      return rc;
    } 
    if (!FC_HANDLE_EQUIV(tempseq, dest_seq)){
      fc_printfErrorMessage("Already have a seq var by this name by on another seq");
      for (i = 0; i < numSeqVar; i++)
	free(seqVars[i]);
      free(seqVars);
      free(numStepPerSeqVar);
      if (dest_seqvarName == NULL) free (dest_svName);
      return rc;
    }
    dest_sv = (FC_Variable*)malloc(numStep*sizeof(FC_Variable));
    for (i = 0; i < numStep; i++){
      dest_sv[i] = seqVars[0][i];
    }
    free(numStepPerSeqVar);
    free(seqVars[0]);
    free(seqVars);
    if (dest_seqvarName == NULL) free (dest_svName);

    //now check this var (will check at all time steps)
    for (i = 0; i < numStep; i++){
      rc = fc_getVariableInfo(dest_sv[i], &currnumDataPoint, &currnumComponent,
			      &currassoc,&currmathtype,&currdatatype);
      if (rc != FC_SUCCESS) {
	fc_printfErrorMessage("Failed to get info for curr variable on dest mesh");
	free(dest_sv);
	return rc;
      }

      if (currnumDataPoint < numDataPoint) {
	printf("dest = %d src = %d\n",currnumDataPoint, numDataPoint);
	fc_printfErrorMessage("%s: destination variable not compatible with src variable", 
			      fc_getReturnCodeText(FC_INPUT_ERROR));
	free(dest_sv);
	return FC_INPUT_ERROR;
      }

      if (currnumComponent != numComponent){
	fc_printfErrorMessage("Cannot copy region var onto mesh because of existing numComponent incompatability");
	free(dest_sv);
	return FC_INPUT_ERROR;
      }

      if (currassoc != assoc){
	fc_printfErrorMessage("Cannot copy region var onto mesh because of existing assoc incompatability");
	free(dest_sv);

	return FC_INPUT_ERROR;
      }

      if (currmathtype != mathtype){
	fc_printfErrorMessage("Cannot copy region var onto mesh because of existing mathtype incompatability");
	free(dest_sv);
	return FC_INPUT_ERROR;
      }

      if (currdatatype != datatype){
	fc_printfErrorMessage("Cannot copy region var onto mesh because of existing datatype incompatability");
	free(dest_sv);
	return FC_INPUT_ERROR;
      }
    }

  } //create or get the dest var

  // --- now do it
  rc = fc_getVariableDataPtr(src_var, &data);
  if (rc != FC_SUCCESS) {
    fc_printfErrorMessage("Failed to get data for variable '%s'",
                          varSlot_s->header.name);
    free(dest_sv);
    return rc;
  }
  rc = fc_getVariableDataPtr(mapping, &mappingdata);
  if (rc != FC_SUCCESS) {
    fc_printfErrorMessage("Failed to get data for mapping variable");
    free(data);
    free(dest_sv);
    return rc;
  }

  //will write it into the current data
  rc = fc_getVariableDataPtr(dest_sv[targetStep], &new_data);
  if (rc != FC_SUCCESS) {
    fc_printfErrorMessage("Failed to get data for dest variable");
    free(dest_sv);
    return rc;
  }
  for (i = 0; i < numDataPoint; i++){
    int meshpt;

    switch(mappingdatatype){
    case FC_DT_DOUBLE:
      meshpt = (int)floor((((double*) mappingdata)[i])+0.5);
      break;
    case FC_DT_FLOAT:
      meshpt = (int)floor((((float*) mappingdata)[i])+0.5);
      break;
    case FC_DT_INT:
      meshpt = ((int*) mappingdata)[i];
      break;
    default:
      fc_printfErrorMessage("DEVELOPER ERROR: should not happen\n");
      return FC_ERROR;
      break;
    }
    //is it a valid pt?
    if (meshpt >= numEntity){
      fc_printfWarningMessage("Invalid mapping val.not copying this point");
    } else {
      switch (datatype){
      case FC_DT_INT:
	for (j = 0; j < numComponent; j++){
	  //	  printf("replacing big mesh pt %d with small mesh pt %d curr val %d to val %d\n",
	  //		 (meshpt*numComponent+j), (i*numComponent+j), 
	  //		 ((int*)new_data)[meshpt*numComponent+j],((int*)data)[i*numComponent+j]);
	  ((int*)new_data)[meshpt*numComponent+j] = ((int*)data)[i*numComponent+j];
	}
	break;
      case FC_DT_DOUBLE:
	for (j = 0; j < numComponent; j++){
	  ((double*)new_data)[meshpt*numComponent+j] = ((double*)data)[i*numComponent+j];
	}
	break;
      case FC_DT_FLOAT:
	for (j = 0; j < numComponent; j++){
	  ((float*)new_data)[meshpt*numComponent+j] = ((float*)data)[i*numComponent+j];
	}
	break;
      case FC_DT_CHAR:
	for (j = 0; j < numComponent; j++){
	  ((char*)new_data)[meshpt*numComponent+j] = ((char*)data)[i*numComponent+j];
	}
	break;
      default:
	fc_printfErrorMessage("DEVELOPER ERROR: should not happen\n");
	return FC_ERROR;
	break;
      }
    }
  }

  *dest_seqvar = dest_sv;
  return FC_SUCCESS;
}


/**
 * \ingroup  Variable
 * \brief  Copy a variable defined on a larger mesh onto a smaller one
 *
 * \description
 *  
 *    Copy a variable defined on a larger mesh onto a smaller one.
 *    Meshes can be on different datasets.
 *    This is intended for cases where you may promote a subset to
 *    a mesh (becuase, for instance, it corresponds to a ROI) 
 *    and want to copy vars from the orginal mesh onto the new one
 *    using the mapping var that is on the smaller one.
 *    THis is particularly useful when the smaller one is created
 *    from a time step of the larger mesh and therefore you
 *    want only the data of a seq var at that time in the
 *    larger mesh to become a var on the smaller mesh.
 *
 *    newName is the name of the var on the dest mesh.
 *    it must not preexist there. 
 *    The vals from the srv_var are copied over
 *    according to the mapping. 
 *
 *    The entity mapping takes the form of another variable of the
 *    same association that contains the ID numbers of the corresponding
 *    entity type in the larger mesh. It is assumed that this is created
 *    at the time of promoting the subset to the mesh. Note that this
 *    is different that the lookup table ids since these go from
 *    currmesh to prev mesh, as opposed to the original on disk 
 *    numbering to the new fclib numbering.
 *
 *    Note this does not check the coords or the conns to see if they
 *    are the same on the 2 meshes, it only goes by the mapping numbers.
 *
 *    Note that since we can only write out variables as attributes in 
 *    exodus and therefore as doubles, the mapping var may not be explictly
 *    FC_DT_INT. We may want to explictly convert this on read in. For now
 *    we will just cast/round the val to an int if it is a float or double.
 *
 *    Note if mapping var in orig (not region) mesh is invalid that
 *    point is not copied.
 *
 *    Fails if the variable already exists on the region mesh.
 *
 * \todo
 *   - fix the rounding. couldnt get round or rint to work.
 *   - prob should put the conversion for the region variable
 *     when it is read in from exodus.
 *   - change so function fails in invalid mapping data.
 *
 * \modifications 
 *   - 06/20/08 ACG Created
 */
FC_ReturnCode fc_copyVariableToRegionMesh(
  FC_Variable src_var, /**< Input - the source variable to be copied */
  FC_Mesh dest_mesh,   /**< Input - the destination mesh to be copied into */
  FC_Variable mapping, /**< Input - int var that contains the orig numbers corresponding to the region mesh */
  void* fillval,       /**< Input - optional variable fill value for entities not in the map */
  char *newName,       /**< Input - the name of the new variable, a NULL value
                            is a flag to use the name of the source variable */
  FC_Variable* new_var /**< Output - handle to the new variable */
) { 
  FC_ReturnCode rc;
  _FC_VarSlot* varSlot;
  FC_Variable* existVars;
  int numExist;
  int numDataPoint, numComponent, numEntity;
  int mappingnumdp,mappingnumc;
  FC_AssociationType assoc, mappingassoc;
  FC_MathType mathtype, mappingmathtype;
  FC_DataType datatype, mappingdatatype;
  void *data, *new_data, *mappingdata;
  int i, j;
  
  // default return value
  if (new_var)
    *new_var = FC_NULL_VARIABLE;

  // check input - NULL name is o.k.
  varSlot = _fc_getVarSlot(src_var);
  if (varSlot == NULL || !fc_isVariableValid(mapping) || !fc_isMeshValid(dest_mesh)
      || new_var == NULL){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // test for NULL name
  if (newName == NULL)
    newName = varSlot->header.name;

  // log message
  fc_printfLogMessage("Copying region variable '%s' to '%s'", 
                      varSlot->header.name, newName);

  //more checking
  rc = fc_getVariableInfo(mapping, &mappingnumdp, &mappingnumc,
			  &mappingassoc,&mappingmathtype,&mappingdatatype);
  if (rc != FC_SUCCESS) {
    fc_printfErrorMessage("Failed to get info for mapping variable");
    return rc;
  }
  rc = fc_getVariableInfo(src_var, &numDataPoint, &numComponent,
			  &assoc,&mathtype,&datatype);
  if (rc != FC_SUCCESS) {
    fc_printfErrorMessage("Failed to get info for src variable");
    return rc;
  }

  //NOTE: from writeout - this may be a double
  //  if (mappingdatatype != FC_DT_INT){
  //    fc_printfErrorMessage("mapping variable must be FC_DT_INT");
  //    return FC_INPUT_ERROR;
  //  }
  if (mappingdatatype != FC_DT_INT && mappingdatatype != FC_DT_DOUBLE &&
      mappingdatatype != FC_DT_FLOAT){
    fc_printfErrorMessage("mappingdatatype must be numeric");
    return FC_INPUT_ERROR;
  }
  if (mappingnumc != 1){
    fc_printfErrorMessage("mapping variable must be single component");
    return FC_INPUT_ERROR;
  }
  if (mappingassoc != assoc){
    fc_printfErrorMessage("mapping variable and src variable must have same assoc");
    return FC_INPUT_ERROR;
  }
  if (datatype != FC_DT_INT && datatype != FC_DT_DOUBLE && datatype != FC_DT_FLOAT &&
      datatype != FC_DT_CHAR){
    fc_printfErrorMessage("src variable must have known type");
    return FC_INPUT_ERROR;
  }
    
  //temporary
  if (assoc == FC_AT_WHOLE_DATASET || assoc == FC_AT_WHOLE_MESH) {
    fc_printfErrorMessage("not handling globals at this time");
    return FC_INPUT_ERROR;
  }
  if (assoc != FC_AT_VERTEX && assoc != FC_AT_ELEMENT){
    fc_printfErrorMessage("only handling vertex and elem vars at this time");
    return FC_INPUT_ERROR;
  }
  //check appropriateness of meshes - are there enough of this entity in
  //the dest mesh (should we check dimension etc?)
  rc = fc_getMeshNumEntity(dest_mesh, assoc, &numEntity); 
  if (rc != FC_SUCCESS) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  if (numEntity > numDataPoint) {
    fc_printfErrorMessage("%s: destination mesh not compatible with src variable", 
			  fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  //does the var already exist?
  rc = fc_getVariableByName(dest_mesh, newName, &numExist, &existVars);
  if (rc != FC_SUCCESS){
    fc_printfErrorMessage("Failed to get info for curr variable on dest mesh");
    return rc;
  }    
  if (numExist > 0){
    fc_printfErrorMessage("Cannot copy var onto new (region) mesh since there is more than 0 instance there");
    free(existVars);
    return FC_ERROR;
  }


  // --- do it
  rc = fc_getVariableDataPtr(src_var, &data);
  if (rc != FC_SUCCESS) {
    fc_printfErrorMessage("Failed to get data for variable '%s'",
                          varSlot->header.name);
    return rc;
  }
  rc = fc_getVariableDataPtr(mapping, &mappingdata);
  if (rc != FC_SUCCESS) {
    fc_printfErrorMessage("Failed to get data for mapping variable");
    free(data);
    return rc;
  }

  // create the new data
  switch(datatype){
  case FC_DT_INT:
    new_data = calloc(numEntity*numComponent,sizeof(int));
    if (new_data == NULL){
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      return FC_MEMORY_ERROR;
    }
    if (fillval != NULL){
      for (i = 0; i < numEntity*numComponent; i++){
	((int*)new_data)[i] = *((int*)fillval);
      }
  }
  break;
 case FC_DT_DOUBLE:
   new_data = calloc(numEntity*numComponent,sizeof(double));
   if (new_data == NULL){
     fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
     return FC_MEMORY_ERROR;
   }
   if (fillval != NULL){
     for (i = 0; i < numEntity*numComponent; i++){
       ((double*)new_data)[i] = *((double*)fillval);
     }
   }
   break;
 case FC_DT_FLOAT:
   new_data = calloc(numEntity*numComponent,sizeof(float));
   if (new_data == NULL){
     fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
     return FC_MEMORY_ERROR;
   }
   if (fillval != NULL){
     for (i = 0; i < numEntity*numComponent; i++){
       ((float*)new_data)[i] = *((float*)fillval);
     }
   }
   break;
 case FC_DT_CHAR:
   new_data = calloc(numEntity*numComponent,sizeof(char));
   if (new_data == NULL){
     fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
     return FC_MEMORY_ERROR;
   }
   if (fillval != NULL){
     for (i = 0; i < numEntity*numComponent; i++){
       ((char*)new_data)[i] = *((char*)fillval);
     }
   }
   break;
 default:
   //shouldnt happen
   fc_printfErrorMessage("Cannot handle datatype %d", datatype);
   return FC_ERROR;
  }

  //look thru the map to find the orig and dest mapping for this var
  for (i = 0; i < numEntity; i++){
    int meshpt;
    switch(mappingdatatype){
    case FC_DT_DOUBLE:
      meshpt = (int)floor((((double*) mappingdata)[i])+0.5);
      break;
    case FC_DT_FLOAT:
      meshpt = (int)floor((((float*) mappingdata)[i])+0.5);
      break;
    case FC_DT_INT:
      meshpt = ((int*) mappingdata)[i];
      break;
    default:
      //wont happen
      break;
    }
    if (meshpt >= numDataPoint){
      fc_printfWarningMessage("Invalid mapping val.not copying this point");
    } else {
      switch (datatype){
      case FC_DT_INT:
	for (j = 0; j < numComponent; j++){
	  ((int*)new_data)[i*numComponent+j] = ((int*)data)[meshpt*numComponent+j];
	}
	break;
      case FC_DT_DOUBLE:
	for (j = 0; j < numComponent; j++){
	  ((double*)new_data)[i*numComponent+j] = ((double*)data)[meshpt*numComponent+j];
	}
	break;
      case FC_DT_FLOAT:
	for (j = 0; j < numComponent; j++){
	  ((float*)new_data)[i*numComponent+j] = ((float*)data)[meshpt*numComponent+j];
	}
	break;
      case FC_DT_CHAR:
	for (j = 0; j < numComponent; j++){
	  ((char*)new_data)[i*numComponent+j] = ((char*)data)[meshpt*numComponent+j];
	}
	break;
      default:
	//wont happen
	break;
      }
    }
  }

  rc = fc_createVariable(dest_mesh, newName, new_var);
  if (rc != FC_SUCCESS) {
    fc_printfErrorMessage("Failed to create new variable '%s'", newName);
    return rc;
  }
    
  rc = fc_setVariableDataPtr(*new_var, numEntity, numComponent, assoc, mathtype, 
			     datatype, new_data);
  if (rc != FC_SUCCESS) {
    fc_printfErrorMessage("Failed to set data for variable '%s'", newName);
    return rc;
  }

  return FC_SUCCESS;
}


/**
 * \ingroup  Variable
 * \brief  Copies a Variable and (optionally) converts to a new association type
 *
 * \description
 *
 *   Create a copy of a variable and (if necessary) convert its data to a new
 *   association type. This conversion simply averages all the values of source entities
 *   that touch a destination entity (eg, in Element-to-Vertex conversion, the an element
 *   value is the average of the values of all vertices that make up the element).
 *
 *   While this function supports conversions to/from any of the association types for
 *   a mesh, it will return an input error if the conversion cannot take place due to
 *   due to dimension problems (eg, trying to convert a variable to a Face type when
 *   the elements are lines).
 *
 *   Additionally, this function only works when the source and destination variables
 *   utilize the same mesh. Global variables are not supported.
 *
 *   Results are always of data type FC_DT_DOUBLE.
 *
 * \modifications  
 *   - 01/30/07 CDU Created
 *   - 02/14/07 CDU Updated to support broader range of associations
 *   - 02/26/07 CDU removed mesh from api and cleaned up internals
 */
FC_ReturnCode fc_copyVariableWithNewAssociation(
    FC_Variable        src_var,     /**< input - the source variable  */
    FC_AssociationType newAssoc,    /**< input - the new association Type  */
    char *             newName,     /**< input - new name for the variable (NULL uses source's name)*/
    FC_Variable*       new_var      /**< output - handle to the new variable. Result is always double precision. */
){

  int  numDataPoint, numComponent;
  double *dnew = NULL;
  int numObj, *objIDs;
  void *var_dptr;
  double r;
  int spot;
  int i,j,k;
  int needs_parent; //For figuring out whether we need parent or child call

  FC_Mesh            mesh; 
  FC_MathType        oldMathType;
  FC_AssociationType oldAssoc;

  FC_DataType dataType;
  FC_ReturnCode rc;
  _FC_VarSlot*  varSlot;


  //Default return value
  if (new_var)
    *new_var = FC_NULL_VARIABLE;

  //Check Inputs: first get basic info and bail if obious problems
  varSlot = _fc_getVarSlot(src_var);
  if((varSlot == NULL)||(new_var==NULL)) {
    fc_printfErrorMessage("%s: problem with %s", 
			  fc_getReturnCodeText(FC_INPUT_ERROR),
			  (varSlot==NULL) ? "input variable" : 
			  "new variable pointer");
    return FC_INPUT_ERROR;
  }

  //Get mesh. We now force the user to use same src/dst mesh
  rc = fc_getMeshFromVariable(src_var, &mesh);
  if(rc!=FC_SUCCESS)
    return rc;
  
  //Get info on original variable
  rc = fc_getVariableInfo(src_var, NULL, &numComponent, 
			  &oldAssoc, &oldMathType, &dataType);
  if(rc!=FC_SUCCESS)
    return rc;
  

  //If no conversion necessary, let copyVariable do all the work
  if (newAssoc == oldAssoc) 
    return fc_copyVariable(src_var, mesh, newName, new_var);
  
  
  //Figure out default name if none given
  if(newName == NULL)
    newName = varSlot->header.name;


  //Look at the destination assoc and figure out how many data points
  //it has. If it's a face/edge, make sure we've created the appropriate conns.
  switch(newAssoc){
  case FC_AT_ELEMENT : 
    rc = fc_getMeshNumElement(mesh, &numDataPoint);
    if(rc!=FC_SUCCESS)
      return rc;
    break;

  case FC_AT_FACE:
    //Build faces if necessary
    rc = fc_getMeshNumFace(mesh, &numDataPoint);
    if(rc!=FC_SUCCESS)
      return rc;
    break;

  case FC_AT_EDGE: 
    //Build edges if necessary
    rc = fc_getMeshNumEdge(mesh, &numDataPoint);
    if(rc != FC_SUCCESS)
	return rc;
    break;

  case FC_AT_VERTEX:
    rc = fc_getMeshNumVertex(mesh, &numDataPoint);
    if(rc!=FC_SUCCESS)
      return rc;
    break;

  default:
    fc_printfErrorMessage("Unexpected conversion target type");
    return FC_INPUT_ERROR;
  }

  //We need at least one data point to continue
  if(numDataPoint<1){
    fc_printfErrorMessage("%s Conversion would not give at least one data point",
			  fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }


  //Figure out what direction we're going. If the new association
  //is less than the old one (newAssoc < oldAssoc), we will need
  //to do a parent-of-newAssoc's ids to find out which oldAssoc ids
  //need to be averaged.
  needs_parent = 0;
  if(newAssoc < oldAssoc)
    needs_parent = 1;


  //Get at the source variable's data
  rc = fc_getVariableDataPtr(src_var, &var_dptr);               
  if(rc!=FC_SUCCESS){
    fc_printfErrorMessage("Failed to get info for variable '%s'",
                          varSlot->header.name);
    return rc;
  } 

  //Make room for the new data
  dnew = (double *)malloc(numDataPoint * numComponent * sizeof(double));
  if (dnew == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }

  //Go through all the new items and calculate their value
  for(i=0; i<numDataPoint; i++){

    //Get the source ids that contain the destination id
    if(needs_parent){
      rc = fc_getMeshEntityParents( mesh, newAssoc,  i, oldAssoc,
				    &numObj, &objIDs);
    } else {
      rc = fc_getMeshEntityChildren(mesh, newAssoc, i, oldAssoc,
				    &numObj, &objIDs);
    }
    if(rc!=FC_SUCCESS){
      fc_printfErrorMessage("Failed to get %s for variable '%s'",
			    (needs_parent)?"Parents":"Children", varSlot->header.name);
      free(dnew);
      return rc;
    }

    //Work on each component individually
    for(j=0; j<numComponent; j++){
      r = 0.0;
      for(k=0; k<numObj; k++){
	spot = (objIDs[k] * numComponent) + j;
	switch(dataType){
	case FC_DT_CHAR:   r += (double) ((char *)  var_dptr)[spot]; break;
	case FC_DT_INT :   r += (double) ((int *)   var_dptr)[spot]; break;
	case FC_DT_FLOAT:  r += (double) ((float *) var_dptr)[spot]; break;
	case FC_DT_DOUBLE: r += (double) ((double *)var_dptr)[spot]; break;
	default:
	  fc_printfErrorMessage("Requested processing of variable with an unknown data type");
	  free(dnew);
	  return FC_INPUT_ERROR;
	}
      }
      if(numObj)
	r = r/(double)numObj;

      //Store result for this component
      dnew[ i*numComponent + j ] = r;

    }
    //Cleanup
    free(objIDs);
  }

  //Make a new variable to stuff this in
  rc = fc_createVariable(mesh, newName, new_var);
  if(rc!=FC_SUCCESS){
    fc_printfErrorMessage("Could not create variable '%s'", newName);
    free(dnew);
    return rc;
  }

  //printf("About to set : numpt %d   numcompon %d  mathtype=%s newassoc=%s\n",
  //	 numDataPoint, numComponent, fc_getMathTypeText(varSlot->mathtype), fc_getAssociationTypeText(newAssoc));
 
  //Fill in the data
  rc = fc_setVariableData(*new_var, numDataPoint, numComponent,
			  newAssoc, oldMathType, FC_DT_DOUBLE,
			  (void *)dnew);

  free(dnew);

  if(rc!=FC_SUCCESS){ 
    fc_printfErrorMessage("Failed to set data for variable '%s'",
			  newName);
    *new_var = FC_NULL_VARIABLE;
    return rc;
  }
  return FC_SUCCESS; 
}



/**
 * \ingroup  Variable
 * \brief  Copy a sequence variable
 *
 * \description
 *  
 *    Create a copy of a sequence variable into the specified mesh using with 
 *    the  given name. This will return an error if the source and
 *    destination sequences and meshes are not similar enough. This will copy a 
 *    sequence if necessary. The destination mesh and sequence need to belong to
 *    the same dataset or an error will be returned.
 *
 *    The input variable can be a global sequence variable; the output will be
 *    a sequence variable associated with the mesh (FC_AT_WHOLE_MESH).
 *
 * \modifications 
 *   - 2003-NOV-11 WSK Created
 *   - 06/07/04 RM, added checking if destination mesh and a source variable
 *          are compatible to perform copy before creating the variable.
 *   - 11/08/2006 WSD added ability to copy global seq var onto a mesh.
 */
FC_ReturnCode fc_copySeqVariable(
  int numStep, /**< Input - the number of steps in the sequence variable */
  FC_Variable* src_seqVar, /**< Input - the source seq variables to be copied */
  FC_Mesh dest_mesh,   /**< Input - the destination mesh to be copied into */
  FC_Sequence dest_seq, /**< Input - the destination seq to be copied into */
  char *newName,       /**< Input - the name of the new variable, a NULL value
                            is a flag to use the name of the source variable */
  FC_Variable** new_seqVar /**< Output - array of handles to the new variable */
) {
  int i; 
  FC_ReturnCode rc;
  int temp_numStep;
  _FC_MeshSlot* dest_meshSlot;
  _FC_SeqSlot* dest_seqSlot;
  _FC_VarSlot* varSlot;
  int numDataPoint, numComp, numEntity;
  FC_AssociationType assoc;
  FC_MathType mathtype;
  FC_DataType datatype;
  void* data;
  
  // default return value
  if (new_seqVar)
    *new_seqVar = NULL;

  // check input
  dest_meshSlot = _fc_getMeshSlot(dest_mesh);
  dest_seqSlot  = _fc_getSeqSlot(dest_seq);
  if (!fc_isSeqVariableValid(numStep, src_seqVar) || dest_meshSlot == NULL || 
      dest_seqSlot == NULL || !FC_HANDLE_EQUIV(dest_meshSlot->ds, dest_seqSlot->ds) ||
      numStep != dest_seqSlot->numStep  || new_seqVar == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // test for NULL name
  varSlot = _fc_getVarSlot(src_seqVar[0]);
  if (newName == NULL)
    newName = varSlot->header.name;

  // log message
  fc_printfLogMessage("Copying sequence variable '%s' to '%s'", 
                      varSlot->header.name, newName);

  // get src meta data
  rc = fc_getVariableInfo(src_seqVar[0], &numDataPoint, &numComp, &assoc,
                            &mathtype, &datatype);
  if (rc != FC_SUCCESS) {
    fc_printfErrorMessage("Failed to get info for seq variable '%s'",
                          varSlot->header.name);
    return rc;
  }
  // more checking - destination mesh and variable assoc type are compatible
  if (assoc == FC_AT_WHOLE_DATASET || assoc == FC_AT_WHOLE_MESH) {
    assoc = FC_AT_WHOLE_MESH;
  }
  else {
    rc = fc_getMeshNumEntity(dest_mesh, assoc, &numEntity); 
    if (rc != FC_SUCCESS) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
      return FC_INPUT_ERROR;
    }
    if (numEntity != numDataPoint) {
      fc_printfErrorMessage("%s: destination mesh not compatible with src variable", 
                            fc_getReturnCodeText(FC_INPUT_ERROR));
      return FC_INPUT_ERROR;
    }
  }

  // --- do it

  // create new variable
  rc = fc_createSeqVariable(dest_mesh, dest_seq, newName, &temp_numStep,
                            new_seqVar);
  if (rc != FC_SUCCESS) {
    fc_printfErrorMessage("Failed to create new seq variable '%s'", newName);
    return rc;
  }
  
  // get and copy data for each step
  for (i = 0; i < numStep; i++) {
    rc = fc_getVariableDataPtr(src_seqVar[i], &data);
    if (rc != FC_SUCCESS) {
      fc_printfErrorMessage("Failed to get data for step %d of seq variable "
                            "'%s'", i, varSlot->header.name);
      return rc;
    }
    
    rc = fc_setVariableData((*new_seqVar)[i], numDataPoint, numComp, assoc, 
                            mathtype, datatype, data);
    if (rc != FC_SUCCESS) {
      fc_printfErrorMessage("Failed to set data for step %d of seq variable "
                            "'%s'", i, newName);
      return rc;
    }
  }

  return rc;
}


/**
 * \ingroup  Variable
 * \brief  Copies a seqVariable and (optionally) converts to a new association type 
 *
 * \description
 *
 *   Create a copy of a sequence variable and (if necessary) convert its data to a new
 *   association type. This conversion simply averages all the values of source entities
 *   that touch a destination entity (eg, in Element-to-Vertex conversion, the an element
 *   value is the average of the values of all vertices that make up the element).
 *
 *   While this function supports conversions to/from any of the association types for
 *   a mesh, it will return an input error if the conversion cannot take place due to
 *   due to dimension problems (eg, trying to convert a variable to a Face type when
 *   the elements are lines).
 *
 *   Additionally, this function only works when the source and destination variables
 *   utilize the same mesh. Global variables are not supported.
 *
 * \modifications  
 *   - 01/30/07 CDU Created
 *   - 02/14/07 CDU Updated to support broader range of associations
 *   - 02/28/07 CDU removed mesh from api and cleaned up internals 
 */

FC_ReturnCode fc_copySeqVariableWithNewAssociation(

    int                src_numStep, /**< input - the number of steps for the sequence variable*/
    FC_Variable*       src_seqVar,  /**< input - the source variable  */
    FC_AssociationType newAssoc,    /**< input - the new association Type */
    char *             newName,     /**< input - new name for the variable (NULL uses source's name)*/
    FC_Variable**      new_seqVar   /**< output - handle to the new sequence variable */
){


  int  numDataPoint, numComponent;
  double *dnew = NULL;
  int numObj, *objIDs;
  void **var_dptr_array;
  void *var_dptr;
  int numStep_out;
  double r;
  int spot;
  int i,j,k,s;
  int needs_parent; //For figuring out whether we need parent or child call

  FC_AssociationType oldAssoc;
  FC_MathType        mathType;
  FC_DataType dataType;
  FC_Sequence seq;
  FC_Variable* seqVar_out;
  FC_Mesh            mesh;
  FC_ReturnCode rc;

  _FC_VarSlot*  varSlot;
 


  //Default return value
  if (new_seqVar)
    *new_seqVar = NULL;

  //Check Inputs: first get basic info and bail if obvious problems
  // note: fc_issv also does a variableValid check on seqvar[0], so we don't need to.
  if((!fc_isSeqVariableValid(src_numStep, src_seqVar)) || 
     (new_seqVar == NULL) ){
    fc_printfErrorMessage("%s: problem with %s",
			  fc_getReturnCodeText(FC_INPUT_ERROR),
			  (new_seqVar == NULL) ? "new seqVariable pointer" :
			  "input sequence variable");
    return FC_INPUT_ERROR;
  }

  //Get var slot. We only use this for accessing the name
  varSlot = _fc_getVarSlot(src_seqVar[0]);
  if(varSlot==NULL){
    fc_printfErrorMessage("Problem with sequence variable\n");
    return FC_INPUT_ERROR;
  }

  //Get mesh. We now force the user to use same src/dst mesh
  rc = fc_getMeshFromVariable(src_seqVar[0], &mesh);
  if(rc!=FC_SUCCESS)
    return rc;


  //Get basic info from variable
  rc = fc_getVariableInfo(src_seqVar[0], NULL, &numComponent, &oldAssoc, &mathType, &dataType);
  if(rc!=FC_SUCCESS)
    return rc;

  //Get Sequence info. Same as old
  rc = fc_getSequenceFromSeqVariable(src_numStep, src_seqVar, &seq);
  if(rc!=FC_SUCCESS)
    return rc;


  //If no conversion necessary, let copySeqVariable do all the work
  if (newAssoc == oldAssoc) {
    return fc_copySeqVariable(src_numStep, src_seqVar, mesh, seq, newName, new_seqVar);
  }

  //Figure out default name if none given
  if(newName == NULL)
    newName = varSlot->header.name;


  //Look at the destination assoc and figure out how many data points
  //it has. If it's a face/edge, make sure we've created the appropriate conns.
  switch(newAssoc){
  case FC_AT_ELEMENT :
    rc = fc_getMeshNumElement(mesh, &numDataPoint);
    if(rc != FC_SUCCESS)
      return rc;
    break;

  case FC_AT_FACE:
    //Builds the faces if necessary
     rc = fc_getMeshNumFace(mesh, &numDataPoint);
    if(rc!=FC_SUCCESS)
      return rc;
    break;

  case FC_AT_EDGE: 
    //Builds the edges if necessary
    rc = fc_getMeshNumEdge(mesh, &numDataPoint);
    if(rc != FC_SUCCESS)
      return rc;
    break;

  case FC_AT_VERTEX:
    rc = fc_getMeshNumVertex(mesh, &numDataPoint);
    if(rc!=FC_SUCCESS)
      return rc;
    break;

  default:
    fc_printfErrorMessage("Unexpected conversion target type");
    return FC_INPUT_ERROR;
  }

  //We need at least one data point to continue
  if(numDataPoint<1){
    fc_printfErrorMessage("%s Conversion would not give at least one data point",
			  fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }


  //Figure out what direction we're going. If the new association
  //is less than the old one (newAssoc < oldAssoc), we will need
  //to do a parent-of-newAssoc's ids to find out which oldAssoc ids
  //need to be averaged.
  needs_parent = 0;
  if(newAssoc < oldAssoc)
    needs_parent = 1;



  // We get all of the varaible data pointers in advance to reduce overhead and
  // mid-loop disaster
  var_dptr_array = (void *) malloc( src_numStep * sizeof(void *));
  if(var_dptr_array == NULL){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  for(s=0; s<src_numStep; s++){
    rc = fc_getVariableDataPtr(src_seqVar[s], &var_dptr_array[s]);               
    if(rc!=FC_SUCCESS){
      fc_printfErrorMessage("Failed to get info for SeqVariable '%s' step %d",
			    varSlot->header.name, s);
      free(var_dptr_array);
      return rc;
    } 
  }


  //Make room for the new data --> dnew[step][vertex][component]
  dnew = (double *)malloc(src_numStep * numDataPoint * numComponent * sizeof(double));
  if(rc!=FC_SUCCESS){
    fc_printfErrorMessage("Failed to get info for variable '%s'",
                          varSlot->header.name);
    free(var_dptr_array);
    return rc;
  } 
  
 
  // In order to reduce the amount of vertex-to-element conversion we have live
  // at one time, we do all sequence steps for a particular vertex. 
  for(i=0; i<numDataPoint; i++){   
    
    if(needs_parent){
      rc = fc_getMeshEntityParents(  mesh, newAssoc, i, oldAssoc,
				     &numObj, &objIDs);
    } else {
      rc = fc_getMeshEntityChildren( mesh, newAssoc, i, oldAssoc,
				     &numObj, &objIDs);
    }
    if(rc!=FC_SUCCESS){
      fc_printfErrorMessage("Failed to get %s for variable '%s'",
			    (needs_parent)?"Parents":"Children", varSlot->header.name);
      free(dnew);
      free(var_dptr_array);
      return rc;
    }

    //Work on all steps for this vertex (or element)
    for(s=0; s<src_numStep; s++){      
      var_dptr = var_dptr_array[s]; //Get at the source variable's data for this step

      //Work on each component individually
      for(j=0; j<numComponent; j++){
	r = 0.0;
	for(k=0; k<numObj; k++){
	  spot = (objIDs[k] * numComponent) + j;
	  switch(dataType){
	  case FC_DT_CHAR:   r += (double) ((char *)  var_dptr)[spot]; break;
	  case FC_DT_INT :   r += (double) ((int *)   var_dptr)[spot]; break;
	  case FC_DT_FLOAT:  r += (double) ((float *) var_dptr)[spot]; break;
	  case FC_DT_DOUBLE: r += (double) ((double *)var_dptr)[spot]; break;
	  default:
	    fc_printfErrorMessage("Requested processing of variable with an unknown data type");
	    free(objIDs);
	    free(dnew);
	    free(var_dptr_array);
	    return FC_INPUT_ERROR;
	  }
	}
	if(numObj)
	  r = r/(double)numObj;
	
	//Store result for this component
	dnew[ (s*numDataPoint*numComponent) + (i*numComponent) + j ] = r;
      }
    }
    //Cleanup this output spot's memory
    free(objIDs);
  }  
  free(var_dptr_array);

  // Create an empty seqVariable to store this
  rc = fc_createSeqVariable(mesh, seq, newName, &numStep_out, &seqVar_out); 
  if(rc!=FC_SUCCESS){
    fc_printfErrorMessage("Could not create seqVariable");
    free(dnew);
    return rc;
  }

  //Fill in all variables
  for(s=0; s<src_numStep; s++){

    //Fill in the data
    rc = fc_setVariableData(seqVar_out[s], numDataPoint, numComponent,
			    newAssoc, mathType, FC_DT_DOUBLE,
			    (void *)&dnew[s*numDataPoint*numComponent]);

    if(rc!=FC_SUCCESS){

      fc_printfErrorMessage("Failed to set data for variable '%s' step %d",
			    newName, s);
      fc_deleteSeqVariable(src_numStep,seqVar_out);
      free(seqVar_out);
      free(dnew);
      return rc;
    }
  }

  free(dnew);
  
  *new_seqVar = seqVar_out;
  return FC_SUCCESS;

}


/**
 * \ingroup  Variable
 * \brief  Convert an array of basic vars into a sequence var
 *
 * \description
 *
 *    This routine converts an array of basic vars into a sequence var.
 *    It is a LOT more efficient than copy the data of basic vars into
 *    the steps of a newly created seq var because no creation/destruction
 *    goes on -- it just changes how the mesh views the variable handles.
 *
 *    You should put the variables in the order they will be on
 *    the sequence. There must be the same number of variables as there
 *    are steps on the sequence. The variables and the sequence have to be
 *    in the same dataset. All of the variables must be on the
 *    same mesh and have the same metadata (except the name). The variables
 *    cannot already be part of a sequence variable.
 *
 *    If no name is supplied, the seqVar will have the name of the
 *    first variable.
 *
 *    (The array of vars that is the seq var should be exactly the
 *    array passed in. A new one is made in case the original array
 *    was automatically generated--fc_deleteSeqVariable would fail
 *    on an automatically allocated array because it tries to free
 *    the array of handles.)
 *
 * \modifications  
 *   - 10/21/04 WSK created.
 *   - 9/21/2006 WSD modified to handle global vars
 */
FC_ReturnCode fc_convertVariablesToSeqVariable(
  int numStep,           /**< input - the number of steps */
  FC_Variable* vars,     /**< input - the variables (array numStep long)
                            (this should be freed by user if malloc'd) */
  FC_Sequence sequence,  /**< input - the sequence (with numStep steps) */
  char* new_name,        /**< input - new name (optional) */
  FC_Variable** seqVar   /**< output - pointer to seq's var handles */
) {
  int i, j;
  int temp_numStep;
  FC_Dataset dataset1, dataset2;
  FC_Mesh mesh;
  _FC_VarSlot *varSlots[numStep];
  int* varMask;
  int *tempNumPer, **tempSlotIDs;
  int *numBasicVar_p, **basicVarIDs_p;
  int *numSeqVar_p, **numStepPerSeqVar_p, ***seqVarIDs_p;

  // default return
  if (seqVar)
    *seqVar = NULL;

  // test input -- a lot of conditions!
  if (numStep < 0 || !vars || !fc_isSequenceValid(sequence) ||
      !seqVar) {
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
    if (!fc_isVariableValid(vars[i])) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
      return FC_INPUT_ERROR;
    }
    varSlots[i] = _fc_getVarSlot(vars[i]);
    if (varSlots[i]->stepID >= 0) {
      fc_printfErrorMessage("Variable is already part of a sequence variable");
      return FC_INPUT_ERROR;
    }
  }
  fc_getDatasetFromSequence(sequence, &dataset1);
  fc_getDatasetFromVariable(vars[0], &dataset2);
  fc_getMeshFromVariable(vars[0], &mesh);
  if (!FC_HANDLE_EQUIV(dataset1, dataset2)) {
    fc_printfErrorMessage("variable and sequence are not in same dataset");
    return FC_INPUT_ERROR;
  }
  for (i = 1; i < numStep; i++) {
    if (!FC_HANDLE_EQUIV(varSlots[i]->ds, dataset1)) {
      fc_printfErrorMessage("variables are not on same dataset");
      return FC_INPUT_ERROR;
    }
    if (!FC_HANDLE_EQUIV(varSlots[i]->mesh, mesh)) {
      fc_printfErrorMessage("variables are not on same mesh");
      return FC_INPUT_ERROR;
    }
    if (varSlots[i]->numDataPoint != varSlots[0]->numDataPoint ||
        varSlots[i]->numComponent != varSlots[0]->numComponent ||
        varSlots[i]->assoc        != varSlots[0]->assoc ||
        varSlots[i]->mathtype     != varSlots[0]->mathtype ||
        varSlots[i]->datatype     != varSlots[0]->datatype) {
      fc_printfErrorMessage("variables have differing metadata");
      return FC_INPUT_ERROR;
    }
  }

  // similar to fc_createSeqVariable()

  // Change names
  if (new_name) {
    for (i = 0; i < numStep; i++) 
      _fc_setSlotHeaderName(&varSlots[i]->header, new_name);
  }
  else {
    for (i = 1; i < numStep; i++)
      _fc_setSlotHeaderName(&varSlots[i]->header, varSlots[0]->header.name);
  }

  // change var contents to reflect new state as seq vars
  for (i = 0; i < numStep; i++) {
    varSlots[i]->sequence = sequence;
    varSlots[i]->stepID = i;
  }

  // get pointers to stuff on mesh or dataset, depending
  if (FC_HANDLE_EQUIV(mesh, FC_NULL_MESH)) {
    _FC_DsSlot *dsSlot = _fc_getDsSlot(varSlots[0]->ds);
    numBasicVar_p = &dsSlot->numBasicVar;
    basicVarIDs_p = &dsSlot->basicVarIDs;
    numSeqVar_p = &dsSlot->numSeqVar;
    numStepPerSeqVar_p = &dsSlot->numStepPerSeqVar;
    seqVarIDs_p = &dsSlot->seqVarIDs;
  }
  else {
    _FC_MeshSlot *meshSlot = _fc_getMeshSlot(varSlots[0]->mesh);
    numBasicVar_p = &meshSlot->numBasicVar;
    basicVarIDs_p = &meshSlot->basicVarIDs;
    numSeqVar_p = &meshSlot->numSeqVar;
    numStepPerSeqVar_p = &meshSlot->numStepPerSeqVar;
    seqVarIDs_p = &meshSlot->seqVarIDs;
  }

  // Add the handles to list of sequence vars
  tempNumPer = (int*)realloc(*numStepPerSeqVar_p,
			     ((*numSeqVar_p)+1)*sizeof(int)); 
  tempSlotIDs = (int**)realloc(*seqVarIDs_p,
			       (*numSeqVar_p+1)*sizeof(int*));
  if (!tempNumPer || !tempSlotIDs) {
    free(tempNumPer); free(tempSlotIDs);
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  tempSlotIDs[*numSeqVar_p] = (int*)malloc(numStep*sizeof(int)); 
  if (!tempSlotIDs[*numSeqVar_p]) {
    free(tempNumPer); free(tempSlotIDs);
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  tempNumPer[*numSeqVar_p] = numStep;
  for (i = 0; i < numStep; i++) 
    tempSlotIDs[*numSeqVar_p][i] = vars[i].slotID;
  *numStepPerSeqVar_p = tempNumPer;
  *seqVarIDs_p = tempSlotIDs;
  *numSeqVar_p = (*numSeqVar_p) + 1;

  // remove handles from list of basic vars
  varMask = calloc(*numBasicVar_p, sizeof(int));
  if (varMask == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  for (i = 0; i < numStep; i++) {
    for (j = 0; j < *numBasicVar_p; j++) {
      if (vars[i].slotID == (*basicVarIDs_p)[j]) {
        varMask[j] = 1;
        break;
      }
    }
    if (j == *numBasicVar_p) {
      printf("**** developer note: something went wrong!!\n");
      fflush(NULL);
      return FC_ERROR;
    }
  }
  j = 0;
  *numBasicVar_p = *numBasicVar_p - numStep;
  if (*numBasicVar_p == 0) {
    free(*basicVarIDs_p);
    *basicVarIDs_p = NULL;
  }
  else {
    for (i = 0; i < *numBasicVar_p + numStep; i++) {
      if (varMask[i] == 0) {
        (*basicVarIDs_p)[j] = (*basicVarIDs_p)[i];
        j++;
      }
    }
  }
  free(varMask);
  // FIX realloc to make smaller

  // make return values
  *seqVar = (FC_Variable*)malloc(numStep*sizeof(FC_Variable));
  for (i = 0; i < numStep; i++)
    (*seqVar)[i] = vars[i];
    
  return FC_SUCCESS;
}

/**
 * \ingroup Variable
 * \brief Creates new variables which are the components of a given variable.
 *    
 * \description 
 * 
 *    Given a vector variable, this routine will create a new variable for
 *    each of its components. For example, if given a vector field like
 *    coordinates, this will create 3 new variables of just the X coordinate,
 *    the Y coordinate and the Z coordinate.
 *
 *    Naming of the new variables is automatic: the name of the original
 *    variable with "_cX" where X is the component ID indexing from 0. 
 *
 * \modifications 
 *    - 05/03/04 WSK, Created.
 *    - 9/22/2006 WSD modified to handle global vars
 *    - 10/12/2007 ACG changed var name from "_componentX" to "_cX"
 */
FC_ReturnCode fc_createComponentVariables(
  FC_Variable var,      /**< input - a variable with non-scalar data */
  int* numComponent,    /**< output - the number of components */
  FC_Variable** components  /**< output - the components as vars */
) {
        
  FC_ReturnCode rc;
  int i, j;
  int numData, dim;
  FC_AssociationType assoc;
  FC_DataType datatype;
  FC_MathType mathtype;
  void *data, **new_data;
  size_t size;
  FC_Dataset dataset;
  FC_Mesh mesh;
  char* var_name, *new_var_name;
  FC_Variable new_var;
  int isGlobal;

  // set default return values
  if (numComponent)
    *numComponent = -1;
  if (components)
    *components = NULL;

  // check input
  if (!fc_isVariableValid(var) || !numComponent || !components){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));   
    return FC_INPUT_ERROR;
  }

  // get variable info
  rc = fc_getVariableInfo(var, &numData, &dim, &assoc, &mathtype, &datatype);
  if (rc != FC_SUCCESS)
    return rc;

  if (dim < 2) {
    fc_printfErrorMessage("Input Error: won't do if only 1 component");
    return FC_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Creating new component variables.");

  // get data
  rc = fc_getVariableDataPtr(var, &data);
  if (rc != FC_SUCCESS)
    return rc; 

  // create the new data
  size = fc_sizeofDataType(datatype);
  new_data = malloc(sizeof(void*)*dim);
  if (new_data == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  for (i = 0; i < dim; i++) {
    new_data[i] = malloc(size*numData);
    if (new_data[i] == NULL) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      return FC_MEMORY_ERROR;
    }
  }
  for (i = 0; i < numData; i++) 
    for (j = 0; j < dim; j++) 
      memcpy((char*)new_data[j]+ size*i, (char*)data + size*(i*dim+j), size); 

  // create the new variables
  isGlobal = fc_isVariableGlobal(var);
  rc = fc_getVariableName(var, &var_name);
  if (rc != FC_SUCCESS)
    return rc;
  *numComponent = dim;
  rc = fc_getDatasetFromVariable(var, &dataset);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_getMeshFromVariable(var, &mesh);
  if (rc != FC_SUCCESS)
    return rc;
  *components = malloc(sizeof(FC_Variable)*dim);
  new_var_name = malloc(sizeof(char)*(strlen(var_name)+20));
  if (*components == NULL || new_var_name == NULL){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  for (i = 0; i < dim; i++) {
    sprintf(new_var_name, "%s_c%d", var_name, i);
    if (isGlobal)
      rc = fc_createGlobalVariable(dataset, new_var_name, &new_var);
    else
      rc = fc_createVariable(mesh, new_var_name, &new_var);
    if (rc != FC_SUCCESS)
      return rc;
    rc = fc_setVariableDataPtr(new_var, numData, 1, assoc, FC_MT_SCALAR, 
                               datatype, (void*)(new_data[i]));
    if (rc != FC_SUCCESS)
      return rc;
    (*components)[i] = new_var;
    // do not free free new_data[i]
  }
  free(new_data);
  free(var_name);
  free(new_var_name);

  return FC_SUCCESS;    
}       

/**
 * \ingroup Variable
 * \brief Creates new sequence variables which are the components of a given
 *     sequence variable.
 *    
 * \description
 *
 *    Given a vector sequence variable, this routine will create a new sequence
 *    variable for each of its components. For example, if given a vector field
 *    like velocity, this will create 3 new sequence variables of just the X
 *    velocity, the Y velocity and the Z velocity.
 *
 *    Naming of the new variables is automatic: the name of the original
 *    variable with "_cX" where X is the component ID indexing from 0. 
 *
 * \modifications 
 *    - 01/19/05 WSK, Created.
 *    - 9/22/2006 WSD modified to handle global vars
 *    - 10/12/2007 ACG changed var name from "_componentX" to "_cX"
 */
FC_ReturnCode fc_createComponentSeqVariables(
  int numStep,    /**< input - number of steps in the sequence variable */
  FC_Variable* seqVar,  /**< input - a seq variable with non-scalar data */
  int* numComponent,    /**< output - the number of components */
  FC_Variable*** compSeqVars  /**< output - the components as seq vars
                               (array of numComp seqVar arrays) */
) {
  FC_ReturnCode rc;
  int i, j;
  int numComp, temp_numComp;
  FC_Sequence sequence;
  char* var_name;
  FC_Variable **comp_vars, *temp_comps;

  // set default return values
  if (numComponent)
    *numComponent = -1;
  if (compSeqVars)
    *compSeqVars = NULL;

  // check input
  if (!fc_isSeqVariableValid(numStep, seqVar) || !numComponent || 
      !compSeqVars) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));   
    return FC_INPUT_ERROR;
  }

  // get variable info
  rc = fc_getVariableNumComponent(seqVar[0], &numComp);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_getSequenceFromSeqVariable(numStep, seqVar, &sequence);
  if (rc != FC_SUCCESS)
    return rc;

  if (numComp < 2) {
    fc_printfErrorMessage("Input Error: won't do if only 1 component");
    return FC_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Creating new component seq variables.");

  // Basic procedure: use fc_createComponentVariables() to break each
  // step of sequence var into components and then convert components
  // back into sequence vars.

  // Make room to save component vars
  comp_vars = (FC_Variable**)malloc(numComp*sizeof(FC_Variable*));
  if (!comp_vars) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  for (i = 0; i < numComp; i++) {
    comp_vars[i] = (FC_Variable*)malloc(numStep*sizeof(FC_Variable));
    if (!comp_vars[i]) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      return FC_MEMORY_ERROR;
    }
  }

  // make component vars
  for (i = 0; i < numStep; i++) {
    rc = fc_createComponentVariables(seqVar[i], &temp_numComp, &temp_comps);
    if (rc != FC_SUCCESS) 
      return FC_ERROR;
    for (j = 0; j < numComp; j++)
      comp_vars[j][i] = temp_comps[j];
    free(temp_comps);
  }

  // make component vars into component seqVars
  *compSeqVars = (FC_Variable**)malloc(numComp*sizeof(FC_Variable*));
  for (i = 0; i < numComp; i++) {
    rc = fc_getVariableName(comp_vars[i][0], &var_name);
    if (rc != FC_SUCCESS)
      return rc;
    rc = fc_convertVariablesToSeqVariable(numStep, comp_vars[i], sequence,
                                          var_name, &((*compSeqVars)[i]));
    free(var_name);
    if (rc != FC_SUCCESS)
      return rc;
    free(comp_vars[i]);
  }
  free(comp_vars);
  *numComponent = numComp;

  return FC_SUCCESS;    
}       

/**
 * \ingroup  Variable
 * \brief Creates a new multicomponent variable from a group of scalars.
 *    
 * \description 
 * 
 *    This routines creates a new multicomponent variable from a
 *    group of scalar variables. The components must all be on the
 *    same mesh and have the same association type and data type.
 *
 * \modifications 
 *    - 05/03/04 WSK, Created.
 *    - 9/22/2006 WSD modified to handle global vars
 */
FC_ReturnCode fc_mergeComponentVariables(
  int numComponent,        /**< input - the number of components */
  FC_Variable* components, /**< input - the components vars */
  char* name,              /**< input - the name of the new merged var */
  FC_MathType mathtype,    /**< input - the mathtype of the new merged var */
  FC_Variable* new_var     /**< output - the new merged var */
) {
        
  FC_ReturnCode rc;
  int i, j;
  int numData, temp_numData, temp_numComp;
  FC_Dataset dataset, temp_dataset;
  FC_Mesh mesh, temp_mesh;
  FC_AssociationType assoc, temp_assoc;
  FC_DataType datatype, temp_datatype;
  void **data, *new_data;
  size_t size;

  // set default return values
  if (new_var)
    *new_var = FC_NULL_VARIABLE;

  // check input
  if (numComponent < 2 || components == NULL || name == NULL || 
      !fc_isMathTypeValid(mathtype) || mathtype == FC_MT_UNKNOWN || 
      mathtype == FC_MT_SCALAR || new_var == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  for (i = 0; i < numComponent; i++) {
    if (!fc_isVariableValid(components[i])) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
      return FC_INPUT_ERROR;
    }
  }
  
  // log message
  fc_printfLogMessage("Merging component variables into new variable.");

  // get & test component variables info
  for (i = 0; i < numComponent; i++) {
    rc = fc_getDatasetFromVariable(components[i], &temp_dataset);
    if (rc != FC_SUCCESS)
      return rc;
    rc = fc_getMeshFromVariable(components[i], &temp_mesh);
    if (rc != FC_SUCCESS)
      return rc;
    if (i == 0) {
      dataset = temp_dataset;
      mesh = temp_mesh;
    }
    else {
      if (!FC_HANDLE_EQUIV(mesh, temp_mesh) || 
	  !FC_HANDLE_EQUIV(dataset, temp_dataset)) {
        fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
        return FC_INPUT_ERROR;
      }  
    }
    rc = fc_getVariableInfo(components[i], &temp_numData, &temp_numComp, 
                              &temp_assoc, NULL, &temp_datatype);
    if (rc != FC_SUCCESS)
      return rc;
    if (temp_numComp != 1) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
      return FC_INPUT_ERROR;
    }
    if (i == 0) {
      numData = temp_numData;
      assoc = temp_assoc;
      datatype = temp_datatype;
    }
    else {
      if (temp_numData != numData || temp_assoc != assoc || 
          temp_datatype != datatype) {
        fc_printfErrorMessage("Components did not have all the same metadata");
        return FC_INPUT_ERROR;
      }
    }
  }

  // get data
  data = malloc(sizeof(void*)*numComponent);
  if (data == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  for (i = 0; i < numComponent; i++) {
    rc = fc_getVariableDataPtr(components[i], &data[i]);
    if (rc != FC_SUCCESS)
      return rc; 
  }

  // create the new data
  size = fc_sizeofDataType(datatype);
  new_data = malloc(size*numComponent*numData);
  if (new_data == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  for (i = 0; i < numComponent; i++) 
    for (j = 0; j < numData; j++) 
      memcpy((char*)new_data + size*(j*numComponent+i), 
             (char*)data[i] + size*j, size); 

  // create the new variable
  if (fc_isVariableGlobal(components[0])) 
    rc = fc_createGlobalVariable(dataset, name, new_var);
  else
    rc = fc_createVariable(mesh, name, new_var);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_setVariableDataPtr(*new_var, numData, numComponent, assoc, mathtype, 
                             datatype, (void*)new_data);
  if (rc != FC_SUCCESS)
    return rc;

  // cleanup
  // don't free new_data; 
  // don't free data's components
  free(data); 

  return FC_SUCCESS;    
}       

/**
 * \ingroup  Variable
 * \brief Creates a new multicomponent seq variable from a group of scalars.
 *    
 * \description
 *
 *    This routines creates a new multicomponent sequence variable from a
 *    group of scalar sequence variables. The components must all be on the
 *    same sequence and on the same mesh, they must all have the same 
 *    association type and data type.
 *
 * \modifications 
 *    - 01/19/05 WSK, Created.
 *    - 9/22/2006 WSD modified to handle global vars
 */
FC_ReturnCode fc_mergeComponentSeqVariables(
  int numComponent,        /**< input - the number of components */
  int numStep,      /**< input - number of steps in the sequence variables */
  FC_Variable** compSeqVars, /**< input - the array of component seqVars */
  char* name,              /**< input - the name of the new merged var */
  FC_MathType mathtype,    /**< input - the mathtype of the new merged var */
  FC_Variable** newSeqVar  /**< output - the new merged seq var */
) {
  FC_ReturnCode rc;
  int i, j;
  int numData, temp_numData, temp_numComp;
  FC_Sequence sequence, temp_sequence;
  FC_Mesh mesh, temp_mesh;
  FC_Variable *merged_vars, *temp_vars;
  FC_AssociationType assoc, temp_assoc;
  FC_DataType datatype, temp_datatype;

  // set default return values
  if (newSeqVar)
    *newSeqVar = NULL;

  // check input
  if (numComponent < 2 || !compSeqVars || !name || 
      !fc_isMathTypeValid(mathtype) || mathtype == FC_MT_UNKNOWN || 
      mathtype == FC_MT_SCALAR || !newSeqVar) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  for (i = 0; i < numComponent; i++) {
    if (!fc_isSeqVariableValid(numStep, compSeqVars[i])) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
      return FC_INPUT_ERROR;
    }
  }
  
  // log message
  fc_printfLogMessage("Merging component seq variables into new seq var.");

  // get & test component seq variable's info
  for (i = 0; i < numComponent; i++) {
    rc = fc_getMeshFromVariable(compSeqVars[i][0], &temp_mesh);
    if (rc != FC_SUCCESS)
      return rc;
    if (i == 0)
      mesh = temp_mesh;
    else {
      if (!FC_HANDLE_EQUIV(mesh, temp_mesh)) {
        fc_printfErrorMessage("Components were not all on the same mesh");
        return FC_INPUT_ERROR;
      }  
    }
    rc = fc_getSequenceFromSeqVariable(numStep, compSeqVars[i],
                                       &temp_sequence);
    if (rc != FC_SUCCESS)
      return rc;
    if (i == 0)
      sequence = temp_sequence;
    else {
      if (!FC_HANDLE_EQUIV(sequence, temp_sequence)) {
        fc_printfErrorMessage("Components were not all on the same sequence");
        return FC_INPUT_ERROR;
      }  
    }
    rc = fc_getVariableInfo(compSeqVars[i][0], &temp_numData, &temp_numComp, 
                            &temp_assoc, NULL, &temp_datatype);
    if (rc != FC_SUCCESS)
      return rc;
    if (temp_numComp != 1) {
      fc_printfErrorMessage("Component %d was not a scalar", i);
      return FC_INPUT_ERROR;
    }
    if (i == 0) {
      numData = temp_numData;
      assoc = temp_assoc;
      datatype = temp_datatype;
    }
    else {
      if (temp_numData != numData || temp_assoc != assoc || 
          temp_datatype != datatype) {
        fc_printfErrorMessage("Components did not have all the same metadata");
        return FC_INPUT_ERROR;
      }
    }
  }

  // Basic procedure: use fc_mergeComponentVariables() to create a merged 
  // variable for each step of the seq vars, then convert merged variables
  // back into a sequence var.

  // Make room to save merged vars
  merged_vars = (FC_Variable*)malloc(numStep*sizeof(FC_Variable));
  if (!merged_vars) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }

  // make merged vars
  temp_vars = (FC_Variable*)malloc(numComponent*sizeof(FC_Variable));
  for (i = 0; i < numStep; i++) {
    for (j = 0; j < numComponent; j++)
      temp_vars[j] = compSeqVars[j][i];
    rc = fc_mergeComponentVariables(numComponent, temp_vars, name, mathtype,
                                    &merged_vars[i]);
    if (rc != FC_SUCCESS)
      return rc;
  }
  free(temp_vars);

  // make merged vars into merged seqVar;
  rc = fc_convertVariablesToSeqVariable(numStep, merged_vars, sequence,
                                        name, newSeqVar);
  free(merged_vars);
  return rc;    
}       

//@}

/** \name Get existing variables. */
//-------------------------------
//@{

/**
 * \ingroup  Variable
 * \brief  Get the global variables in a dataset
 *
 * \description
 *  
 *    Returns an array of variables in the dataset. The caller is
 *    responsible for freeing the array of variable handles.
 *
 * \modifications  
 *   - 9/19/06 WSD created.
 */
FC_ReturnCode fc_getGlobalVariables(
 FC_Dataset dataset,     /**< input - dataset */
 int *numVar,      /**< output - number of variables */
 FC_Variable **variables  /**< output - array of variable handles */
) {
  int i;
  _FC_DsSlot* dsSlot;
  _FC_VarSlot* varSlot;
  
  // set up defaults
  if (numVar != NULL)
    *numVar = -1;
  if (variables != NULL)
    *variables = NULL;
  
  // check input
  dsSlot = _fc_getDsSlot(dataset);
  if (dsSlot == NULL || numVar == NULL || variables == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Getting global variables for dataset '%s'",
                      dsSlot->header.name);

  // get num variable
  *numVar = dsSlot->numBasicVar;
  
  // get variables
  if (dsSlot->numBasicVar > 0) {
    *variables = malloc(sizeof(FC_Variable) * dsSlot->numBasicVar);
    if (*variables == NULL) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      return FC_MEMORY_ERROR;
    }
    for (i = 0; i < dsSlot->numBasicVar; i++) {
      varSlot = _fc_getVarSlotFromID(dsSlot->basicVarIDs[i]);
      _FC_GET_HANDLE((*variables)[i], varSlot);
    }
  }
  
  return FC_SUCCESS;
}

/**
 * \ingroup  Variable
 * \brief  Get the global sequence variables in a dataset.
 *
 * \description
 *  
 *    Returns an array of global seq variables (array of arrays of variables)
 *    as well as an array of the number of steps per sequence variable.
 *    The use is responsible for freeing all returned arrays.
 *    (SeqVars is a 2D array, so to free: for (i = 0; i < numSeq var) {
 *    free(seqVars[i])); } )
 *
 * \modifications  
 *   - 9/21/2006 WSD Created.
 */
FC_ReturnCode fc_getGlobalSeqVariables(
 FC_Dataset dataset,   /**< input - dataset */
 int *numSeqVar,       /**< output - number of sequence variables */
 int** numStepPerVar,  /**< output - array with number of steps for each
                          sequence variable */
 FC_Variable ***seqVariables  /**< output - array of sequence variables (each
                          sequence variable  is array of numStep variables */
) {
  int i, j;
  _FC_DsSlot* dsSlot;
  _FC_VarSlot* varSlot;
  
  // set up defaults
  if (numSeqVar != NULL)
    *numSeqVar = -1;
  if (numStepPerVar != NULL)
    *numStepPerVar = NULL;
  if (seqVariables != NULL)
    *seqVariables = NULL;
  
  // check input
  dsSlot = _fc_getDsSlot(dataset);
  if (dsSlot == NULL || numSeqVar == NULL || numStepPerVar == NULL ||
      seqVariables == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Getting global sequence variables for dataset '%s'",
                      dsSlot->header.name);

  // get num seq variable
  *numSeqVar = dsSlot->numSeqVar;
  
  // get sequence variables
  if (*numSeqVar > 0) {
    *numStepPerVar = malloc(sizeof(int) * (*numSeqVar));
    if (*numStepPerVar == NULL) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      return FC_MEMORY_ERROR;
    }
    *seqVariables = malloc(sizeof(FC_Variable*) * (*numSeqVar));
    if (*seqVariables == NULL) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      return FC_MEMORY_ERROR;
    }
    for (i = 0; i < *numSeqVar; i++) {
      (*numStepPerVar)[i] = dsSlot->numStepPerSeqVar[i];
      (*seqVariables)[i] = malloc(dsSlot->numStepPerSeqVar[i]*sizeof(FC_Variable));
      if ( (*seqVariables)[i] == NULL) {
	fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
	return FC_MEMORY_ERROR;
      }
      for (j = 0; j < dsSlot->numStepPerSeqVar[i]; j++) {
	varSlot = _fc_getVarSlotFromID(dsSlot->seqVarIDs[i][j]);
	_FC_GET_HANDLE((*seqVariables)[i][j], varSlot);
      }
    }
  }
  
  return FC_SUCCESS;
}

/**
 * \ingroup  Variable
 * \brief  Get the variables in a mesh.
 *
 * \description
 *  
 *    Returns an array of variables in the mesh. The caller is
 *    responsible for freeing the array of variable handles.
 *
 * \modifications  
 *   - Nancy Collins Created.
 *   - 2003-NOV-14  WSK  Changed to only names of basic variables, not all.
 *   - 12/11/03 RM, added default return values, numVar = -1, varNames = NULL
 *            check input/output for numVar = NULL      
 *   - 9/8/04 WSK, changed to return variable handles instead of variable
 *               names.
 */
FC_ReturnCode fc_getVariables(
 FC_Mesh mesh,     /**< input - mesh */
 int *numVar,      /**< output - number of variables */
 FC_Variable **variables  /**< output - array of variable handles */
) {
  int i;
  _FC_MeshSlot* meshSlot;
  _FC_VarSlot* varSlot;
  
  // set up defaults
  if (numVar != NULL)
    *numVar = -1;
  if (variables != NULL)
    *variables = NULL;
  
  // check input
  meshSlot = _fc_getMeshSlot(mesh);
  if (meshSlot == NULL || numVar == NULL || variables == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Getting variables for mesh '%s'",
                      meshSlot->header.name);

  // get num variable
  *numVar = meshSlot->numBasicVar;
  
  // get variables
  if (meshSlot->numBasicVar > 0) {
    *variables = malloc(sizeof(FC_Variable) * meshSlot->numBasicVar);
    if (*variables == NULL) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      return FC_MEMORY_ERROR;
    }
    for (i = 0; i < meshSlot->numBasicVar; i++) {
      varSlot = _fc_getVarSlotFromID(meshSlot->basicVarIDs[i]);
      _FC_GET_HANDLE((*variables)[i], varSlot);
    }
  }
  
  return FC_SUCCESS;
}

/**
 * \ingroup  Variable
 * \brief  Get the sequence variables in a mesh.
 *
 * \description
 *  
 *    Returns an array of sequence variables (array of arrays of variables)
 *    as well as an array of the number of steps per sequence variable.
 *    The use is responsible for freeing all returned arrays.
 *    (SeqVars is a 2D array, so to free: for (i = 0; i < numSeq var) {
 *    free(seqVars[i])); } )
 *
 * \modifications  
 *   - 2003-APR-18  W Koegler  Created
 *   - 12/11/03 RM, added default return values, numSeqVar = -1, 
 *               seqVarNames = NULL, check input/output for numSeqVar = NULL.
 *   - 9/9/04 WSK, changed to return seq variable handles instead of names
 *   - 10/15/04 WSK, changed to not return sequence handles
 */
FC_ReturnCode fc_getSeqVariables(
 FC_Mesh mesh,        /**< input - mesh */
 int *numSeqVar,      /**< output - number of sequence variables */
 int** numStepPerVar,  /**< output - array with number of steps for each
                          sequence variable */
 FC_Variable ***seqVariables  /**< output - array of sequence variables (each
                          sequence variable  is array of numStep variables */
) {
  int i, j;
  _FC_MeshSlot* meshSlot;
  _FC_VarSlot* varSlot;
  
  // set up defaults
  if (numSeqVar != NULL)
    *numSeqVar = -1;
  if (numStepPerVar != NULL)
    *numStepPerVar = NULL;
  if (seqVariables != NULL)
    *seqVariables = NULL;
  
  // check input
  meshSlot = _fc_getMeshSlot(mesh);
  if (meshSlot == NULL || numSeqVar == NULL || numStepPerVar == NULL ||
      seqVariables == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Getting sequence variables for mesh '%s'",
                      meshSlot->header.name);

  // get num seq variable
  *numSeqVar = meshSlot->numSeqVar;
  
  // get sequence variables
  if (*numSeqVar > 0) {
    *numStepPerVar = malloc(sizeof(int) * (*numSeqVar));
    if (*numStepPerVar == NULL) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      return FC_MEMORY_ERROR;
    }
    *seqVariables = malloc(sizeof(FC_Variable*) * (*numSeqVar));
    if (*seqVariables == NULL) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      return FC_MEMORY_ERROR;
    }
    for (i = 0; i < *numSeqVar; i++) {
      (*numStepPerVar)[i] = meshSlot->numStepPerSeqVar[i];
      (*seqVariables)[i] = malloc(meshSlot->numStepPerSeqVar[i]*sizeof(FC_Variable));
      if ( (*seqVariables)[i] == NULL) {
	fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
	return FC_MEMORY_ERROR;
      }
      for (j = 0; j < meshSlot->numStepPerSeqVar[i]; j++) {
	varSlot = _fc_getVarSlotFromID(meshSlot->seqVarIDs[i][j]);
	_FC_GET_HANDLE((*seqVariables)[i][j], varSlot);
      }
    }
  }
  
  return FC_SUCCESS;
}

/**
 * \ingroup  Variable
 * \brief  Get the number of global variables in a dataset.
 *
 * \modifications  
 *   - 9/19/06 WSD created.
 */
FC_ReturnCode fc_getNumGlobalVariable(
 FC_Dataset dataset, /**< input - dataset */
 int *numVar         /**< output - number of variables */
) {
  _FC_DsSlot* dsSlot;
  
  // set up defaults
  if (numVar != NULL)
    *numVar = -1;
  
  // check input
  dsSlot = _fc_getDsSlot(dataset);
  if (dsSlot == NULL || numVar == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Getting number of global variables for dataset '%s'",
                      dsSlot->header.name);

  // get num variable
  *numVar = dsSlot->numBasicVar;
  
  return FC_SUCCESS;
}

/**
 * \ingroup  Variable
 * \brief  Get the number of global sequence variables in a dataset.
 *
 * \modifications  
 *   - 9/21/2006 WSD created.
 */
FC_ReturnCode fc_getNumGlobalSeqVariable(
 FC_Dataset dataset,  /**< input - dataset */
 int *numSeqVar       /**< output - number of sequence variables */
) {
  _FC_DsSlot* dsSlot;
  
  // set up defaults
  if (numSeqVar != NULL)
    *numSeqVar = -1;
  
  // check input
  dsSlot = _fc_getDsSlot(dataset);
  if (dsSlot == NULL || numSeqVar == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Getting sequence variables for dataset '%s'",
                      dsSlot->header.name);

  // get num seq variable
  *numSeqVar = dsSlot->numSeqVar;
    
  return FC_SUCCESS;
}

/**
 * \ingroup  Variable
 * \brief  Get the number of variables in a mesh.
 *
 * \modifications  
 *   - 9/28/04 WSK created.
 */
FC_ReturnCode fc_getNumVariable(
 FC_Mesh mesh,     /**< input - mesh */
 int *numVar       /**< output - number of variables */
) {
  _FC_MeshSlot* meshSlot;
  
  // set up defaults
  if (numVar != NULL)
    *numVar = -1;
  
  // check input
  meshSlot = _fc_getMeshSlot(mesh);
  if (meshSlot == NULL || numVar == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Getting number of variables for mesh '%s'",
                      meshSlot->header.name);

  // get num variable
  *numVar = meshSlot->numBasicVar;
  
  return FC_SUCCESS;
}

/**
 * \ingroup  Variable
 * \brief  Get the number of sequence variables in a mesh.
 *
 * \modifications  
 *   - 9/28/04 WSK, created.
 */
FC_ReturnCode fc_getNumSeqVariable(
 FC_Mesh mesh,        /**< input - mesh */
 int *numSeqVar       /**< output - number of sequence variables */
) {
  _FC_MeshSlot* meshSlot;
  
  // set up defaults
  if (numSeqVar != NULL)
    *numSeqVar = -1;
  
  // check input
  meshSlot = _fc_getMeshSlot(mesh);
  if (meshSlot == NULL || numSeqVar == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Getting sequence variables for mesh '%s'",
                      meshSlot->header.name);

  // get num seq variable
  *numSeqVar = meshSlot->numSeqVar;
    
  return FC_SUCCESS;
}

/**
 * \ingroup  Variable
 * \brief  Get global variables from the dataset.
 *
 * \description
 *  
 *    Returns array of global variables on the dataset with the specified name.
 *
 * \modifications  
 *   - 9/19/06 WSD created.
 *   = 9/22/06 ACG returns array
 */
FC_ReturnCode fc_getGlobalVariableByName(
  FC_Dataset dataset,     /**< input - dataset handle */
  char *varname,          /**< input - variable name */
  int *numVariables,       /**< output - num variables */
  FC_Variable **variables   /**< output - variable handle */
) {
  _FC_DsSlot* dsSlot;
  _FC_VarSlot *varSlot;
  FC_Variable variable;
  int i;
  
  FC_Variable *lmatch = NULL;
  int numlmatch, maxnumlmatch;

  // default values returned if anything goes wrong
  if (numVariables)
    *numVariables = -1;
  if (variables)
    *variables = NULL;
  
  // check input
  dsSlot = _fc_getDsSlot(dataset);
  if (dsSlot == NULL || varname == NULL || !variables || ! numVariables){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Getting variable '%s'", varname);

  numlmatch = 0;
  maxnumlmatch = 0;
  
  // find the existing slot
  for (i = 0; i < dsSlot->numBasicVar; i++) {
    varSlot = _fc_getVarSlotFromID(dsSlot->basicVarIDs[i]); 
    if (!strcmp(varSlot->header.name, varname)) {
      _FC_GET_HANDLE(variable, varSlot);

      if (lmatch == NULL){
	lmatch = (FC_Variable*)malloc(sizeof(FC_Variable));
	if (lmatch == NULL){
	  fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
	  return FC_MEMORY_ERROR;
	}
	numlmatch = 1;
	maxnumlmatch = 1;
	lmatch[0] = variable;
      }else{
	if(numlmatch == maxnumlmatch){
	  FC_Variable *temp;
	  temp = (FC_Variable*)realloc(lmatch,(2*numlmatch)*sizeof(FC_Variable));
	  if (temp == NULL){
	    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
	    return FC_MEMORY_ERROR;
	  }
	  maxnumlmatch*=2;
	  lmatch = temp;
	}
	lmatch[numlmatch] = variable;
	numlmatch++;
      }
    }
  }

  *numVariables = numlmatch;
  *variables = lmatch;

  return FC_SUCCESS;
}

/**
 * \ingroup  Variable
 * \brief  Get global sequence variables with a given name from the dataset.
 *
 * \description
 *  
 *    Returns an array of global seq variables (array of arrays of variables)
 *    as well as an array of the number of steps per sequence variable.
 *    The use is responsible for freeing all returned arrays.
 *    (SeqVars is a 2D array, so to free: for (i = 0; i < numSeq var) {
 *    free(seqVars[i])); } )
 *
 * \modifications  
 *   - 9/21/2006 WSD created
 *   - 9/25/2006 ACG returns array
 */
FC_ReturnCode fc_getGlobalSeqVariableByName(
  FC_Dataset dataset,        /**< input - dataset handle */
  char *varName,             /**< input - seq variable name */
  int *numSeqVar,       /**< output - number of sequence variables */
  int** numStepPerVar,  /**< output - array with number of steps for each
                          sequence variable */
  FC_Variable ***seqVariables  /**< output - array of sequence variables (each
                          sequence variable  is array of numStep variables */
) {
  int i, j;
  _FC_DsSlot *dsSlot;
  _FC_VarSlot *varSlot;

  int numStep;
  FC_Variable *seqvar;
  FC_Variable **lmatch = NULL;
  int *lstepmatch = NULL;
  int numlmatch, maxnumlmatch;
  
  // default return values
 if (numSeqVar != NULL)
    *numSeqVar = -1;
  if (numStepPerVar != NULL)
    *numStepPerVar = NULL;
  if (seqVariables != NULL)
    *seqVariables = NULL;
  
  // check input
  dsSlot = _fc_getDsSlot(dataset);
  if (dsSlot == NULL || !varName ||  numSeqVar == NULL || numStepPerVar == NULL ||
      seqVariables == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Getting sequence variable '%s'", varName);

  numlmatch = 0;
  maxnumlmatch = 0;
  
  // look for varName in this dataset's seq var list
  for (i = 0; i < dsSlot->numSeqVar; i++) {
    varSlot = _fc_getVarSlotFromID(dsSlot->seqVarIDs[i][0]);
    if (!strcmp(varName, varSlot->header.name)) {
      numStep = dsSlot->numStepPerSeqVar[i];
      seqvar = malloc((numStep)*sizeof(FC_Variable));
      if (seqvar == NULL) {
	fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
	return FC_MEMORY_ERROR;
      }
      for (j = 0; j < numStep; j++) {
	varSlot = _fc_getVarSlotFromID(dsSlot->seqVarIDs[i][j]);
	_FC_GET_HANDLE(seqvar[j], varSlot);
      }

      if (lmatch == NULL){
	lmatch = (FC_Variable**)malloc(sizeof(FC_Variable*));
	lstepmatch = (int*)malloc(sizeof(int));
	if (lmatch == NULL || lstepmatch == NULL){
	  fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
	  return FC_MEMORY_ERROR;
	}
	numlmatch = 1;
	maxnumlmatch = 1;
	lmatch[0] = seqvar;
	lstepmatch[0] = numStep;
      }else{
	if(numlmatch == maxnumlmatch){
	  FC_Variable **temp;
	  int *tempstep;
	  temp = (FC_Variable**)realloc(lmatch,(2*numlmatch)*sizeof(FC_Variable*));
	  tempstep = (int*)realloc(lstepmatch,(2*numlmatch)*sizeof(FC_Variable*));
	  if (temp == NULL || tempstep == NULL){
	    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
	    return FC_MEMORY_ERROR;
	  }
	  maxnumlmatch*=2;
	  lmatch = temp;
	  lstepmatch = tempstep;
	}
	lmatch[numlmatch] = seqvar;
	lstepmatch[numlmatch] = numStep;
	numlmatch++;
      }
    }
  }

  *numSeqVar = numlmatch;
  *numStepPerVar = lstepmatch;
  *seqVariables = lmatch;

  return FC_SUCCESS;
}

/**
 * \ingroup  Variable
 * \brief  Get variables from the mesh.
 *
 * \description
 *  
 *    Returns variables on the mesh with the specified name.
 *    Note that this doesnt merge components like the seqvar
 *    version doesx
 *
 * \modifications  
 *   - Nancy Collins  Created.
 *   - 2003-MAY-30  W Koegler  Moved getting variable meta data from 
 *        fc_getMeshByName() to here. Actually moved to _fc_readVariable() 
 *        which is called here.
 *   - 12/12/03 RM, added default return mesh = FC_NULL_VARIABLE  
 *        and checking for invalid variable name. Eliminated 
 *        _fc_lookup_varSlot() and _fc_namematch(), those input checking 
 *        changed and added into the body of this function  
 *   - 09/22/06 ACG returns array
 */
FC_ReturnCode fc_getVariableByName(
  FC_Mesh mesh,           /**< input - mesh handle */
  char *varname,          /**< input - variable name */
  int *numVariables,      /**< output - num matching variables */
  FC_Variable **variables   /**< output - variable handles */
) {
  _FC_MeshSlot* meshSlot;
  _FC_VarSlot *varSlot;
  FC_Variable variable;
  int i;

  FC_Variable *lmatch = NULL;
  int numlmatch, maxnumlmatch;
  
  // default values returned if anything goes wrong
  if (numVariables)
    *numVariables = -1;
  if (variables)
    *variables = NULL;
  
  // check input
  meshSlot = _fc_getMeshSlot(mesh);
  if (meshSlot == NULL || varname == NULL || !variables || !numVariables){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Getting variable '%s'", varname);
  
  numlmatch = 0;
  maxnumlmatch = 0;
  
  // find the existing slot
  for (i = 0; i < meshSlot->numBasicVar; i++) {
    varSlot = _fc_getVarSlotFromID(meshSlot->basicVarIDs[i]); 
    if (!strcmp(varSlot->header.name, varname)) {
      _FC_GET_HANDLE(variable, varSlot);

      if (lmatch == NULL){
	lmatch = (FC_Variable*)malloc(sizeof(FC_Variable));
	if (lmatch == NULL){
	  fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
	  return FC_MEMORY_ERROR;
	}
	numlmatch = 1;
	maxnumlmatch = 1;
	lmatch[0] = variable;
      }else{
	if(numlmatch == maxnumlmatch){
	  FC_Variable *temp;
	  temp = (FC_Variable*)realloc(lmatch,(2*numlmatch)*sizeof(FC_Variable));
	  if (temp == NULL){
	    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
	    return FC_MEMORY_ERROR;
	  }
	  maxnumlmatch*=2;
	  lmatch = temp;
	}
	lmatch[numlmatch] = variable;
	numlmatch++;
      }
    }
  }

  *numVariables = numlmatch;
  *variables = lmatch;

  return FC_SUCCESS;
}


/**
 * \ingroup  Variable
 * \brief  Get variables that can be components from the mesh.
 *
 * \description
 *  
 *    Returns variables on the mesh with the specified base 
 *    name and _x,_y,_z and similiar extensions.
 *    This is meant to be used to build an array of components
 *    split with the exodus component naming convention.
 *
 *
 * \modifications  
 *   - 10/11/2007 ACG created
 */
FC_ReturnCode fc_getVariableComponentsByName(
  FC_Mesh mesh,           /**< input - mesh handle */
  char *varname,          /**< input - variable name */
  FC_MathType *type,      /**< output - type */
  int *numVariables,      /**< output - num matching variables */
  FC_Variable **variables   /**< output - variable handles */
) {
  _FC_MeshSlot* meshSlot;
  FC_Variable* returnVars;
  char* compName;


  int numComp;
  unsigned int lastCompInGuess;
  unsigned int lastCompInSearch;
  unsigned int maxExtLen;
  unsigned int maxComp;
  unsigned int varNameLen;
  
  unsigned int extState;


  FC_ReturnCode rc;
  
  // default values returned if anything goes wrong
  if (numVariables)
    *numVariables = -1;
  if (type)
    *type = FC_MT_UNKNOWN; 
  if (variables)
    *variables = NULL;
  
  // check input
  meshSlot = _fc_getMeshSlot(mesh);
  if (meshSlot == NULL || varname == NULL || !variables || !numVariables
      || !type){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Getting variable components '%s'", varname);


  //First, query component function to see what the guesses look like
  rc = _fc_getNextComponentExtension(NULL, &maxExtLen, &maxComp, NULL);
  if(rc!=FC_SUCCESS){
    return rc;
  }

  //Allocate working space for temporary data
  compName   = (char *)       malloc(strlen(varname)+maxExtLen);
  returnVars = (FC_Variable *)malloc(maxComp*sizeof(FC_Variable));
  if (!compName || !returnVars){
    free(compName); free(returnVars);
    return FC_MEMORY_ERROR;
  }
  varNameLen = strlen(varname); 
  strcpy(compName, varname); //search name always starts with varname

  for(numComp=0, extState=0, lastCompInGuess=0, lastCompInSearch=0; !lastCompInSearch; ){
    int numVar;
    FC_Variable *tempVars;

    rc = _fc_getNextComponentExtension(&extState,  
				       &lastCompInGuess, &lastCompInSearch, 
				       &compName[varNameLen]);
    if(rc!=FC_SUCCESS){
      free(compName);
      return rc;
    }
    //printf("Checking '%s'\n",compName);

    //Look for this variable
    rc = fc_getVariableByName(mesh, compName,&numVar, &tempVars);
    if (rc != FC_SUCCESS){
      free(compName);
      return rc;
    }
   
    if(numVar>1){
      //Multiple occurrences of this variable
      free(compName); free(returnVars); free(tempVars); 
      return FC_INPUT_ERROR;

    } else if (numVar == 1){
      //Found another proper match in this guess
      returnVars[numComp] = tempVars[0];
      free(tempVars);
      numComp++;

      if(lastCompInGuess){
	//We hit all components, stop searching. Resize the data
	//printf("Found complete set\n");
	*variables = (FC_Variable *)realloc(returnVars, numComp*sizeof(FC_Variable));
	free(compName);
	if(!(*variables)){
	  free(returnVars); 
	  return FC_MEMORY_ERROR;
	}
	*type = FC_MT_VECTOR;
	*numVariables = numComp;
	return FC_SUCCESS;
      }

    } else {
      //found nothing, forget this guess
      //Nothing to free, just tell guessing engine to skip to next
      while(!lastCompInGuess){
	rc = _fc_getNextComponentExtension(&extState, &lastCompInGuess, &lastCompInSearch, NULL);
	if(rc!=FC_SUCCESS){
	  free(compName); free(returnVars); 
	  return rc;
	}
      }
      numComp=0;
    }
  }
  //printf("Not found\n");
  //Ran through all guesses, but no matches
  *numVariables = 0;
  free(compName); free(returnVars);
  return FC_SUCCESS;
}




/**
 * \ingroup Variable
 * \brief Gets a unique var if it exists, otherwise generate it from components if possible
 *
 * \description 
 *
 *  Gets a unique var if it exists, otherwise generate it
 *  from components if possible. 
 *
 *  Warning: If the unique var is generated
 *  from components, the original components are deleted. 
 *
 * \todo 
 *  - should send back a flag if it was generated or not
 *    (that may make a difference if you go to release it)
 *
 * \modifications
 * - 10/15/2007 ACG created
 */
FC_ReturnCode fc_getOrGenerateUniqueVariableByName(
  FC_Mesh mesh,           /**< input - mesh handle */
  char *varname,          /**< input - variable name */
  FC_Variable *variable   /**< output - variable handle */
) {
  _FC_MeshSlot* meshSlot;
  FC_Variable* returnVars;
  int numReturnVars;
  FC_MathType mtype;
  FC_Variable generatedVar;

  int i;
  FC_ReturnCode rc;
  
  // default values returned if anything goes wrong
  if (variable)
    *variable = FC_NULL_VARIABLE;
  
  // check input
  meshSlot = _fc_getMeshSlot(mesh);
  if (meshSlot == NULL || varname == NULL || !variable){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Getting unique variable by name '%s'", varname);
  rc = fc_getVariableByName(mesh, varname, &numReturnVars,
			    &returnVars);
  if (rc != FC_SUCCESS) {
    fc_printfErrorMessage("can't get variable by name");
    return rc;
  }    

  if (numReturnVars > 1){
    fc_printfWarningMessage("More than 1 variable named '%s'",varname);
    free(returnVars);
    return FC_SUCCESS;
  }

  if (numReturnVars == 1){
    *variable = returnVars[0];
    free(returnVars);
    return FC_SUCCESS;
  }

  //try to make it from components
  rc = fc_getVariableComponentsByName(mesh, varname, &mtype,
				    &numReturnVars, &returnVars);
  if (rc != FC_SUCCESS){
    fc_printfErrorMessage("can't get variable components by name '%s'",varname);
    return rc;
  }

  if (numReturnVars == 0){
    //no components
    return FC_SUCCESS;
  }


  rc = fc_mergeComponentVariables(numReturnVars,
				  returnVars, varname,
				  mtype, &generatedVar);
  if (rc != FC_SUCCESS){
    fc_printfErrorMessage("can't merge component variables");
    free(returnVars);
    return rc;
  }

  //now kill off the component vars
  for (i = 0; i < numReturnVars; i++){
    rc = fc_deleteVariable(returnVars[i]);
    if (rc != FC_SUCCESS){
      fc_printfErrorMessage("cannot delete component variable");
      return rc;
    }
  }
  free(returnVars);


  *variable = generatedVar;
  return FC_SUCCESS;
}


/**
 * \ingroup  Variable
 * \brief  Get sequence variables of a given name from the mesh.
 *
 * \description
 *
 *    Returns an array of seq variables (array of arrays of variables)
 *    as well as an array of the number of steps per sequence variable.
 *    The use is responsible for freeing all returned arrays.
 *    (SeqVars is a 2D array, so to free: for (i = 0; i < numSeq var) {
 *    free(seqVars[i])); } )
 *
 *
 * \modifications  
 *   - 2003-APR-21  W Koegler  Created.
 *   - 12/15/03 RM, added checking if seq variable with the specified 
 *        name was found
 *   - 10/15/04 WSK, changed to not return sequence handles
 *   - 09/25/06 ACG returns array
 *   - 10/14/07 ACG removed hack 
 */
FC_ReturnCode fc_getSeqVariableByName(
  FC_Mesh mesh,              /**< input - mesh handle */
  char *varName,             /**< input - seq variable name */
  int *numSeqVar,       /**< output - number of sequence variables */
  int** numStepPerVar,  /**< output - array with number of steps for each
                          sequence variable */
  FC_Variable ***seqVariables  /**< output - array of sequence variables (each
				  sequence variable is array of numStep
				  variables */
) {
  int i, j;
  _FC_MeshSlot *meshSlot;
  _FC_VarSlot *varSlot;

  int numStep;
  FC_Variable *seqvar;
  FC_Variable **lmatch = NULL;
  int *lstepmatch = NULL;
  int numlmatch, maxnumlmatch;
  
  // default return values
 if (numSeqVar != NULL)
    *numSeqVar = -1;
  if (numStepPerVar != NULL)
    *numStepPerVar = NULL;
  if (seqVariables != NULL)
    *seqVariables = NULL;

    // check input
  meshSlot = _fc_getMeshSlot(mesh);
  if (meshSlot == NULL || !varName || numSeqVar == NULL || 
      numStepPerVar == NULL || seqVariables == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Getting sequence variable '%s'", varName);

  numlmatch = 0;
  maxnumlmatch = 0;

  
  // look for varName in this mesh's seq var list
  for (i = 0; i < meshSlot->numSeqVar; i++) {
    varSlot = _fc_getVarSlotFromID(meshSlot->seqVarIDs[i][0]);
    if (!strcmp(varName, varSlot->header.name)) {
      numStep = meshSlot->numStepPerSeqVar[i];
      seqvar = malloc((numStep)*sizeof(FC_Variable));
      if (seqvar == NULL) {
	fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
	return FC_MEMORY_ERROR;
      }
      for (j = 0; j < numStep; j++) {
	varSlot = _fc_getVarSlotFromID(meshSlot->seqVarIDs[i][j]);
	_FC_GET_HANDLE(seqvar[j], varSlot);
      }

      if (lmatch == NULL){
	lmatch = (FC_Variable**)malloc(sizeof(FC_Variable*));
	lstepmatch = (int*)malloc(sizeof(int));
	if (lmatch == NULL || lstepmatch == NULL){
	  fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
	  return FC_MEMORY_ERROR;
	}
	numlmatch = 1;
	maxnumlmatch = 1;
	lmatch[0] = seqvar;
	lstepmatch[0] = numStep;
      }else{
	if(numlmatch == maxnumlmatch){
	  FC_Variable **temp;
	  int *tempstep;
	  temp = (FC_Variable**)realloc(lmatch,(2*numlmatch)*sizeof(FC_Variable*));
	  tempstep = (int*)realloc(lstepmatch,(2*numlmatch)*sizeof(FC_Variable*));
	  if (temp == NULL || tempstep == NULL){
	    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
	    return FC_MEMORY_ERROR;
	  }
	  maxnumlmatch*=2;
	  lmatch = temp;
	  lstepmatch = tempstep;
	}
	lmatch[numlmatch] = seqvar;
	lstepmatch[numlmatch] = numStep;
	numlmatch++;
      }
    }
  }

  if (numlmatch){
    *numSeqVar = numlmatch;
    *numStepPerVar = lstepmatch;
    *seqVariables = lmatch;
    return FC_SUCCESS;
  }

  //no matches
  *numSeqVar = 0;
  *numStepPerVar = NULL;
  *seqVariables = NULL;
  return FC_SUCCESS;

}


/**
 * \ingroup  Variable
 * \brief  Get seq variables that can be components from the mesh.
 *
 * \description
 *  
 *    Returns seq variables on the mesh with the specified base 
 *    name and _x,_y,_z and similiar extensions.
 *    This is meant to be used to build an array of components
 *    split with the exodus component naming convention.
 *    This function supports a few common naming extensions
 *    (eg, xyz, xy, c0c1c2, tuvw) and tries a few spacing
 *    options (eg, ending in _x, __x, or just x). 
 *
 * \modifications  
 *   - 10/12/2007 ACG created
 *   - 12/11/2007 CDU added basic guessing logic that tries different names
 *   - 12/12/2007 CDU Moved guessing logic to external function
 */


FC_ReturnCode fc_getSeqVariableComponentsByName(
  FC_Mesh mesh,           /**< input - mesh handle */
  char *varname,          /**< input - variable name */
  FC_MathType *type,      /**< output - math type */ 
  int *numVariables,      /**< output - num matching variables */
  int **numStepsPerVariable,  /**< output - numSteps per variable */     
  FC_Variable ***variables   /**< output - variable handles */
) {  
  _FC_MeshSlot* meshSlot;
  FC_Variable** returnVars;

  int *nSteps;
  char* compName;


  int  j;
  int numComp;
  unsigned int lastCompInGuess;
  unsigned int lastCompInSearch;
  unsigned int maxExtLen;
  unsigned int maxComp;
  unsigned int varNameLen;
  
  unsigned int extState;
  FC_ReturnCode rc;
  
  // default values returned if anything goes wrong
  if (type)
    *type = FC_MT_UNKNOWN;
  if (numVariables)
    *numVariables = -1;
  if (variables)
    *variables = NULL;
  if (numStepsPerVariable)
    *numStepsPerVariable = NULL;
  
  // check input
  meshSlot = _fc_getMeshSlot(mesh);
  if (meshSlot == NULL || varname == NULL || !variables || !numVariables
      || ! numStepsPerVariable || !type){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Getting variable components '%s'", varname);

  //First, query component function to see what the guesses look like
  rc = _fc_getNextComponentExtension(NULL, &maxExtLen, &maxComp, NULL);
  if(rc!=FC_SUCCESS){
    return rc;
  }


  //Allocate working space for temporary data
  compName   = (char *)        malloc(strlen(varname)+maxExtLen);
  returnVars = (FC_Variable **)malloc(maxComp*sizeof(FC_Variable*));
  nSteps     = (int *)         malloc(maxComp*sizeof(int));
  if (!compName || !returnVars || ! nSteps){
    free(compName); free(returnVars); free(nSteps);
    return FC_MEMORY_ERROR;
  }
  varNameLen = strlen(varname);
  strcpy(compName, varname); //search name always starts with varname
  

  //Keep asking for extensions until we run out of choices
  for(numComp=0, extState=0, lastCompInGuess=0, lastCompInSearch=0; !lastCompInSearch; ){
    int numVar;
    int* numStep;
    FC_Variable** tempVar;     

    rc = _fc_getNextComponentExtension(&extState,  
				       &lastCompInGuess, &lastCompInSearch, 
				       &compName[varNameLen]);
    if(rc!=FC_SUCCESS){
      free(compName); free(returnVars); free(nSteps);
      return rc;
    }
    //printf("Checking '%s'\n",compName);
    
    //Search for this component
    rc = fc_getSeqVariableByName(mesh, compName, &numVar, &numStep, &tempVar);
    if(rc != FC_SUCCESS){
      for (j = 0; j < numComp; j++)
	free(returnVars[j]);	  
      free(compName); free(nSteps); free(returnVars);
      return rc;
    }

    if (numVar > 1){
      //printf("found too many %d\n",i);
      for (j = 0; j < numVar; j++)
	free(tempVar[j]);	
      free(tempVar);  free(numStep);
      for (j = 0; j < numComp; j++)
	free(returnVars[j]);      
      free(compName); free(returnVars); free(nSteps);
      return FC_INPUT_ERROR;

    } else if (numVar == 1 ){
      //printf("found %d\n",numComp);
      returnVars[numComp] = tempVar[0];
      nSteps[numComp] = numStep[0];
      free(tempVar);
      free(numStep);
      numComp++;

      if(lastCompInGuess){
	//We hit on all components, stop searching, trim data down to size
	*variables =  (FC_Variable **)realloc(returnVars, numComp*sizeof(FC_Variable **));
	*numStepsPerVariable = (int *)realloc(nSteps,     numComp*sizeof(int));
	free(compName);
	if(!(*variables)||!(*numStepsPerVariable)){
	  if(!(*variables))           free(returnVars);
	  else                        free(*variables);
	  if(!(*numStepsPerVariable)) free(nSteps);
	  else                        free(*numStepsPerVariable);
	  return FC_MEMORY_ERROR;
	}

	*type = FC_MT_VECTOR;
	*numVariables = numComp;
	return FC_SUCCESS;

      }

    } else {
      //Didn't find, bail out
      for(j=0; j<numComp; j++)
	free(returnVars[j]);
      
      //Skip through all of this row's remaining components
      while(!lastCompInGuess){
	rc = _fc_getNextComponentExtension(&extState, &lastCompInGuess, &lastCompInSearch, NULL);
	if(rc!=FC_SUCCESS){
	  free(compName); free(returnVars); free(nSteps);
	  return rc;
	}
      }
      numComp=0;
    }


  }

  //Ran through all guesses, but no matches
  *numVariables = 0;
  free(compName); free(returnVars); free(nSteps);
  return FC_SUCCESS;

}





/**
 * \ingroup Variable
 * \brief Gets a unique seq var if it exists, otherwise generate it from components if possible
 *
 * \description 
 *
 *  Gets a unique seq var if it exists, otherwise generate it
 *  from components if possible. 
 *
 *  Warning: If the unique seq var is generated
 *  from components, the original components are deleted. 
 *
 * \todo 
 *  - should send back a flag if it was generated or not
 *    (that may make a difference if you go to release it)
 *
 * \modifications
 * - 10/15/2007 ACG created
 */
FC_ReturnCode fc_getOrGenerateUniqueSeqVariableByName(
  FC_Mesh mesh,           /**< input - mesh handle */
  char *varname,          /**< input - variable name */
  int *numStep,           /**< out - num step */
  FC_Variable **variable   /**< output - variable handle */
) {
  _FC_MeshSlot* meshSlot;
  FC_Variable** returnVars;
  int numReturnVars, *numReturnSteps;
  FC_MathType mtype;
  FC_Variable *generatedVar;

  int i;
  FC_ReturnCode rc;
  
  // default values returned if anything goes wrong
  if (variable)
    *variable = NULL;
  if (numStep)
    *numStep = -1;
  
  // check input
  meshSlot = _fc_getMeshSlot(mesh);
  if (meshSlot == NULL || varname == NULL || !variable || !numStep){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Getting unique seq variable by name '%s'", varname);
  rc = fc_getSeqVariableByName(mesh, varname, &numReturnVars,
			       &numReturnSteps, &returnVars);
  if (rc != FC_SUCCESS) {
    fc_printfErrorMessage("can't get seq variable by name");
    return rc;
  }    

  if (numReturnVars > 1){
    fc_printfWarningMessage("More than 1 seq variable named '%s'",varname);
    free(returnVars);
    free(numReturnSteps);
    return FC_SUCCESS;
  }

  if (numReturnVars == 1){
    *variable = returnVars[0];
    *numStep = numReturnSteps[0];
    free(numReturnSteps);
    free(returnVars);
    return FC_SUCCESS;
  }

  //try to make it from components
  rc = fc_getSeqVariableComponentsByName(mesh, varname, &mtype,
					 &numReturnVars, &numReturnSteps,
					 &returnVars);
  if (rc != FC_SUCCESS){
    fc_printfErrorMessage("can't get seq variable components by name '%s'",varname);
    return rc;
  }

  if (numReturnVars == 0){
    //no components
    return FC_SUCCESS;
  }

  rc = fc_mergeComponentSeqVariables(numReturnVars,
				     numReturnSteps[0],
				     returnVars, varname,
				     mtype, &generatedVar);
  if (rc != FC_SUCCESS){
    fc_printfErrorMessage("can't merge component variables");
    free(numReturnSteps);
    free(returnVars);
    return rc;
  }

  //now kill off the component vars
  for (i = 0; i < numReturnVars; i++){
    rc = fc_deleteSeqVariable(numReturnSteps[i],returnVars[i]);
    if (rc != FC_SUCCESS){
      fc_printfErrorMessage("cannot delete component variable");
      return rc;
    }
    free(returnVars[i]);
  }
  free(returnVars);

  *variable = generatedVar;
  *numStep = numReturnSteps[0];
  free(numReturnSteps);
  return FC_SUCCESS;
}

//@}

/** \name Change the name of a mesh. */
//-------------------------------
//@{

/**
 * \ingroup Variable
 * \brief  Change the name of a variable.
 *
 * \modifications
 *   - 8/24/2005 WSD. Created.
 */
FC_ReturnCode fc_changeVariableName(
  FC_Variable variable, /**< Input - the variable. */
  char* newName       /**< Input - new name for the variable. */
) {
  _FC_VarSlot* varSlot;

  // check input
  varSlot = _fc_getVarSlot(variable);
  if (varSlot == NULL || newName == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // log message
  fc_printfLogMessage("Changing name of variable from '%s' to '%s'", 
                      varSlot->header.name, newName);
  
  // Do it
  return _fc_setSlotHeaderName(&varSlot->header, newName);
}

/**
 * \ingroup Variable
 * \brief  Change the name of a sequence variable.
 *
 * \modifications
 *   - 8/24/2005 WSD. Created.
 */
FC_ReturnCode fc_changeSeqVariableName(
  int numStep, /**< Input - Number of steps in the seq variable. */
  FC_Variable* variables, /**< Input - the steps of the seq variable. */
  char* newName       /**< Input - new name for the seq variable. */
) {
  FC_ReturnCode rc, rc_keep = FC_SUCCESS;
  int i;
  _FC_VarSlot* varSlot;

  // check input
  if (!fc_isSeqVariableValid(numStep, variables) || newName == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // log message
  varSlot = _fc_getVarSlot(variables[0]); 
  fc_printfLogMessage("Changing name of seq variable from '%s' to '%s'", 
                      varSlot->header.name, newName);
  
  // Do it
  for (i = 0; i < numStep; i++) {
    rc = _fc_setSlotHeaderName(&varSlot->header, newName);
    if (rc != FC_SUCCESS)
      rc_keep = rc;
  }
  return rc_keep;
}

//@}

/** \name Release/delete meshes. */
//-------------------------------
//@{

/**
 * \ingroup  Variable
 * \brief  Attempt to minimize variable's memory usage.
 *
 * \description
 *  
 *    Call to try and release un-needed memory in a  variable. 
 *    If the variable has been saved to disk, the large data is 
 *    released. If the variable has not been saved to disk, this 
 *    will do nothing.
 *
 * \modifications  
 *    - 04/19/04 WSK Created.
 */
FC_ReturnCode fc_releaseVariable(
 FC_Variable variable       /**< input - variable handle */
) {
  _FC_VarSlot* varSlot;
  
  // special case: releasing a null handle is not an error
  if (FC_HANDLE_EQUIV(variable, FC_NULL_VARIABLE)) {
    fc_printfWarningMessage("Releasing FC_NULL_VARIABLE");
    return FC_SUCCESS;
  }

  // check input
  varSlot = _fc_getVarSlot(variable);
  if (varSlot == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Cleaning up variable '%s'", varSlot->header.name);

  // remove big data if this slot is committed
  if (varSlot->header.committed == 1) {
    free(varSlot->data);
    varSlot->data = NULL;
  }

  return FC_SUCCESS;
}

/**
 * \ingroup  Variable
 * \brief  Attempt to minimize variable's memory usage.
 *
 * \description
 *  
 *    Call to try and release un-needed memory in a sequence variable. 
 *    If the variable has been saved to disk, the large data is 
 *    released. If the variable has not been saved to disk, this 
 *    will do nothing.
 *
 *    fc_releaseVariable() can be called on specific steps of the 
 *    sequence variable.
 *
 * \modifications  
 *    - 04/19/04 WSK Created.
 */
FC_ReturnCode fc_releaseSeqVariable(
  int numStep, /**< Input - Number of steps in the seq variable. */
  FC_Variable* variables /**< Input - the seq var. */
) {
  FC_ReturnCode rc, rc_keep = FC_SUCCESS;
  int i;  
  _FC_VarSlot* varSlot;
  
  // special case: releasing a null handle is not an error
  if (numStep < 1 || variables == NULL) {
    fc_printfWarningMessage("Releasing NULL seq variable");
    return FC_SUCCESS;
  }

  // check input
  if (!fc_isSeqVariableValid(numStep, variables)) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // log message
  varSlot = _fc_getVarSlot(variables[0]);
  fc_printfLogMessage("Cleaning up variable '%s'", varSlot->header.name);

  // Do it
  rc_keep = FC_SUCCESS;
  for (i = 0; i < numStep; i++) {
    rc = fc_releaseVariable(variables[i]);
    if (rc != FC_SUCCESS)
      rc_keep = rc;
  }

  return rc_keep;
}

/**
 * \ingroup  Variable
 * \brief  Delete a variable.
 *
 * \description
 *  
 *    Call when this variable is no longer needed.
 *
 * \modifications  
 *    - Nancy Collins Created.
 *    - 12/01/03 WSK Changed so it removes uncommitted variables
 *      form the owning mesh's varID lists.
 *    - 11/25/03 RM, if variable is FC_NULL_VARIABLE, then do nothing and 
 *         return success (issue warning if fc_getLibraryVerbosity()). 
 *    - 9/9/04 WSK changed name to delete. Changed behavior so that
 *      delete always deletes no matter whether on disk or not.
 *    - 9/19/2006 WSD Added behavior so handles variables owned by
 *        datasets as well.
 */
FC_ReturnCode fc_deleteVariable(
 FC_Variable variable       /**< input - variable handle */
) {
  int i, j;
  _FC_VarSlot* varSlot;
  int* num_p;     // pointer to num location
  int** slots_p;  // pointer to slotIDs location
  char *extra = "";
 
  // special case: deleting a null handle is not an error
  if ( FC_HANDLE_EQUIV(variable, FC_NULL_VARIABLE)) {
    fc_printfWarningMessage("Deleting FC_NULL_VARIABLE");
    return FC_SUCCESS;
  }
 
  // check input
  varSlot = _fc_getVarSlot(variable);
  if (varSlot == NULL || varSlot->stepID > -1) { // can't be seq var
    fc_printfErrorMessage("%s, invalid variable", 
             fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  if (FC_HANDLE_EQUIV(varSlot->mesh, FC_NULL_MESH))
    extra = " global";
  fc_printfLogMessage("Deleting%s variable '%s'", extra, 
		      varSlot->header.name);

  // remove reference from dataset or mesh
  if (FC_HANDLE_EQUIV(varSlot->mesh, FC_NULL_MESH)) {
    _FC_DsSlot *dsSlot = _fc_getDsSlot(varSlot->ds);
    num_p = &dsSlot->numBasicVar;
    slots_p = &dsSlot->basicVarIDs;
  }
  else {
    _FC_MeshSlot* meshSlot = _fc_getMeshSlot(varSlot->mesh);
    num_p = &meshSlot->numBasicVar;
    slots_p = &meshSlot->basicVarIDs;
    // special case: the coord field
    if (meshSlot->coordVarID == varSlot->header.slotID)
      meshSlot->coordVarID = -1;
  }
  if (*num_p == 1) {
    free(*slots_p);
    *slots_p = NULL;
    *num_p = 0;
  }
  else {
    for (i = 0; i < *num_p; i++) {
      if ((*slots_p)[i] == varSlot->header.slotID)
	break;
    }
    for (j = i+1; j < *num_p; j++)
      (*slots_p)[j-1] = (*slots_p)[j];
    *num_p = *num_p - 1;
  }

  // clear & delete slot
  return _fc_deleteVarSlot(variable);
}

/**
 * \ingroup  Variable
 * \brief  Delete a sequence variable.
 *
 * \description
 *  
 *    Call when this sequence variable is no longer needed.
 *
 *    NOTE: this does not delete the array of variable handles that makes
 *    up the sequence variable, the user is still responsible for this.
 *
 * \modifications  
 *    - 2003-APR-21  W Koegler  Created.
 *    - 12/01/03 WSK Changed so it removes uncommitted variables
 *      form the owning mesh's varID lists.
 *    - 9/9/04 WSK changed name to delete. Changed behavior so that
 *      delete always deletes no matter whether on disk or not.
 *    - 9/19/2006 WSD Added behavior so handles variables owned by
 *        datasets as well.
 */
FC_ReturnCode fc_deleteSeqVariable(
 int numStep,               /**< input - the number of variable handles */
 FC_Variable* variables     /**< input - variable handles */
) {
  FC_ReturnCode rc, rc_keep;
  int i, j;
  _FC_VarSlot* first_varSlot;
  int* num_p;     // pointer to num
  int** numPer_p; // pointer to numStep per seqVar
  int*** slots_p;  // pointer to slotIDs
  char *extra = "";

  // special case: deleting a null list is not an error
  if (numStep <= 0 && variables == NULL) {
    fc_printfWarningMessage("Deleting a NULL sequence variable");
    return FC_SUCCESS;
  }

  // check input
  if (!fc_isSeqVariableValid(numStep, variables)) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  first_varSlot = _fc_getVarSlot(variables[0]);
  
  // log message
  if (FC_HANDLE_EQUIV(first_varSlot->mesh, FC_NULL_MESH))
    extra = " global";
  fc_printfLogMessage("Deleting%s seq variable '%s'", extra, 
                      first_varSlot->header.name);

  // remove references from dataset or mesh
  if (FC_HANDLE_EQUIV(first_varSlot->mesh, FC_NULL_MESH)) {
    _FC_DsSlot *dsSlot = _fc_getDsSlot(first_varSlot->ds);
    num_p = &dsSlot->numSeqVar;
    numPer_p = &dsSlot->numStepPerSeqVar;
    slots_p = &dsSlot->seqVarIDs;
  }
  else {
    _FC_MeshSlot *meshSlot = _fc_getMeshSlot(first_varSlot->mesh);
    num_p = &meshSlot->numSeqVar;
    numPer_p = &meshSlot->numStepPerSeqVar;
    slots_p = &meshSlot->seqVarIDs;
  }
  if (*num_p == 1) {
    free(*numPer_p);
    *numPer_p = NULL;
    free((*slots_p)[0]);
    free(*slots_p);
    *slots_p = NULL;
    *num_p = 0;
  }
  else {
    for (i = 0; i < *num_p; i++) {
      if ((*slots_p)[i][0] == first_varSlot->header.slotID)
	break;
    }
    free((*slots_p)[i]);
    for (j = i+1; j < *num_p; j++) {
      (*numPer_p)[j-1] = (*numPer_p)[j];
      (*slots_p)[j-1] = (*slots_p)[j];
    }
    *num_p = *num_p - 1;
  }
  
  // clear & delete slots
  rc_keep = FC_SUCCESS;
  for (i = 0; i < numStep; i++) {
    rc = _fc_deleteVarSlot(variables[i]);
    if (rc != FC_SUCCESS)
      rc_keep = rc;
  }
  
  return rc_keep;
}

//@}

/** \name Get variable metadata. */
//-------------------------------
//@{

/**
 * \ingroup  Variable
 * \brief Check that the handle refers to a valid variable.
 *
 * \description
 *
 *    Returns 1 (true) if the handle is valid, and 0 (false) if it is not.
 *
 * \modifications  
 *    - 05/25/04 WSK Created.
 */
int fc_isVariableValid(
  FC_Variable variable  /**< input - variable handle */
) {
  return _fc_isHandleValid(variable.slotID, variable.uID, varTableSize,
                           (_FC_SlotHeader**)varTable);
}

/**
 * \ingroup  Variable
 * \brief Check that the handle refers to a valid variable.
 *
 * \description
 *
 *    Returns 1 (true) if the seq variable is valid, and 0 (false) if it is 
 *    not. (Valid means numStep is correct and all variable handles are
 *    valid and in correct order).
 *
 * \modifications  
 *    - 09/23/04 WSK Created.
 *    - 9/19/2006 WSD Added behavior so handles variables owned by
 *        datasets as well.
 */
int fc_isSeqVariableValid(
  int numStep,         /**< input - number of steps in seq var */
  FC_Variable *seqVar  /**< input - seq var (array of vars) */
) {
  int i;
  _FC_VarSlot* varSlot;
  _FC_SeqSlot* seqSlot;
  int var_id;  // id's into mesh's slotID arrays
  int num, **slots;

  // quick check
  if (numStep < 0 || seqVar == NULL || !fc_isVariableValid(seqVar[0])) 
    return 0;
  // check that numStep agrees with sequence
  varSlot = _fc_getVarSlot(seqVar[0]);
  seqSlot = _fc_getSeqSlot(varSlot->sequence);
  if (numStep != seqSlot->numStep)
    return 0;
  // check that all steps are valid and are from the same variable
  // instead of checking metadata, we check against the parent's list
  if (FC_HANDLE_EQUIV(varSlot->mesh, FC_NULL_MESH)) {
    _FC_DsSlot *dsSlot =_fc_getDsSlot(varSlot->ds);
    num = dsSlot->numSeqVar;
    slots = dsSlot->seqVarIDs;
  }
  else {
    _FC_MeshSlot *meshSlot = _fc_getMeshSlot(varSlot->mesh);
    num = meshSlot->numSeqVar;
    slots = meshSlot->seqVarIDs;
  }
  var_id = -1;
  for (i = 0; i < num; i++) {
    if (varSlot->header.slotID == slots[i][0]) {
      var_id = i;
      break;
    }
  }
  if (var_id == -1) {
    fc_printfErrorMessage("bad code? shouldn't ever happen?");
    return FC_ERROR;
  }
  for (i = 1; i < numStep; i++) {
    if (!fc_isVariableValid(seqVar[i]) ||
	seqVar[i].slotID != slots[var_id][i])
      return 0;
  }

  // passed all tests
  return 1;
}

/**
 * \ingroup  Variable
 * \brief Check that whether a variable is global.
 *
 * \description
 *
 *    Returns 1 (true) if the variable is a global variable, 0 (false) if it is
 *    not, and an error code (< 0) if the handle is not valid.
 *
 *    NOTE: The test for a global variable is that the owning mesh is NULL.
 *
 * \modifications  
 *    - 09/21/2006
 */
int fc_isVariableGlobal(
  FC_Variable variable  /**< input - variable handle */
) {
  _FC_VarSlot *varSlot;

  // check input
  varSlot = _fc_getVarSlot(variable);
  if (varSlot == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  if (FC_HANDLE_EQUIV(varSlot->mesh, FC_NULL_MESH))
    return 1;
  else
    return 0;
}

/**
 * \ingroup  Variable
 * \brief  Return name of a variable.
 *
 * \description
 *  
 *    Return the name associated with a variable handle. Space
 *    is allocated for the string so it needs to be freed by the
 *    user when finished.
 *
 * \modifications  
 *   - Nancy Collins, Created.  
 *   - Aug, 9, 2002. W Koegler  Used to be _fc_get_vid_name. Now it's public.
 *   - 12/16/03 RM, added default NULL return value
 */
FC_ReturnCode fc_getVariableName(
  FC_Variable variable,   /**< input - variable handle */
  char **varName          /**< output - variable name, space allocated here
                             and must be freed by user when finished */
) {
  char *cp;
  _FC_VarSlot* varSlot;
  
  // set default return value
  if (varName)
    *varName = NULL;
  
  // check input
  varSlot = _fc_getVarSlot(variable);
  if (varSlot == NULL || varName == NULL) { 
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Getting name for variable '%s'", varSlot->header.name);

  // copy name
  cp = varSlot->header.name;
  *varName = malloc(strlen(cp) + 1);
  if (*varName == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  strcpy(*varName, cp);

  return FC_SUCCESS;
}

/**
 * \ingroup  Variable
 * \brief  Get the parent sequence of a sequence variable.
 *
 * \modifications
 *   - 09/23/04 WSK, created.
 */
FC_ReturnCode fc_getSequenceFromSeqVariable(
  int numStep,           /**< input - number of steps in sequence variable */
  FC_Variable* seqVar,   /**< input - sequence variable (array of vars) */
  FC_Sequence *sequence  /**< output - sequence handle */
) {
  _FC_VarSlot* varSlot;
  
  // set default value
  if (sequence != NULL)
    *sequence = FC_NULL_SEQUENCE;
  
  // check input
  if (!fc_isSeqVariableValid(numStep, seqVar) || sequence == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  varSlot = _fc_getVarSlot(seqVar[0]);
  fc_printfLogMessage("Getting sequence from seq variable '%s'", 
                      varSlot->header.name);
  
  // return requested value
  *sequence = varSlot->sequence;
  
  return FC_SUCCESS;
}

/**
 * \ingroup  Variable
 * \brief  Get the parent dataset of the variable.
 *
 * \modifications  
 *   - 9/19/2006 WSD Created.
 */
FC_ReturnCode fc_getDatasetFromVariable(
  FC_Variable variable,   /**< input - variable handle */
  FC_Dataset *dataset     /**< output - dataset handle */
) {
  _FC_VarSlot* varSlot;
  
  // set default value
  if (dataset != NULL)
    *dataset = FC_NULL_DATASET;
  
  // check input
  varSlot = _fc_getVarSlot(variable);
  if (varSlot == NULL || dataset == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Getting dataset from variable '%s'", varSlot->header.name);
  
  // return requested value
  *dataset = varSlot->ds;
  
  return FC_SUCCESS;
}

/**
 * \ingroup  Variable
 * \brief  Get the parent mesh of the variable.
 *
 * \description
 *
 *    If the variable is a global variable, this will return FC_NULL_MESH
 *    but will not error.
 *
 * \modifications  
 *   - Nancy Collins, Created.
 *   - 12/16/03 RM, added setting default value for mesh = FC_NULL_MESH
 *       and checking for NULL input/output
 */
FC_ReturnCode fc_getMeshFromVariable(
  FC_Variable variable,   /**< input - variable handle */
  FC_Mesh *mesh           /**< output - mesh handle */
) {
  _FC_VarSlot* varSlot;
  
  // set default value
  if (mesh != NULL)
    *mesh = FC_NULL_MESH;
  
  // check input
  varSlot = _fc_getVarSlot(variable);
  if (varSlot == NULL || mesh == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Getting mesh from variable '%s'", varSlot->header.name);
  
  // return requested value
  *mesh = varSlot->mesh;
  
  return FC_SUCCESS;
}

/**
 * \ingroup  Variable
 * \brief  Get information about a variable.
 *
 * \description
 *
 *    Returns the number of data points and number of components in the
 *    variable (the total number of values is numDataPoint*numComponent).
 *    Also returns the association type, math type, and data type.
 *
 * \modifications  
 *   - Nancy Collins, Created.
 *   - 12/16/03 RM, added setting default values to FC_**_UNKNOWN for
 *       data types and -1 for number of data points and number of components
 */
FC_ReturnCode fc_getVariableInfo(
  FC_Variable variable,     /**< input - variable handle */
  int *numDataPoint,        /**< output - the number of data points */
  int *numComponent,        /**< output - the number of components */
  FC_AssociationType *assoc,  /**< output - data association (e.g. vertex,
                                          element, etc.) */
  FC_MathType *mathtype,    /**< output - the organization of components 
                                 (e.g. scalar, vector, etc.) */
  FC_DataType *datatype     /**< output - data type (e.g. int, float, etc.) */
) {
  _FC_VarSlot* varSlot;
  
  // set default values
  if (datatype != NULL)
    *datatype = FC_DT_UNKNOWN;
  if (assoc != NULL)
    *assoc = FC_AT_UNKNOWN;
  if (numDataPoint != NULL)
    *numDataPoint = -1;
  if (numComponent != NULL)
    *numComponent = -1;
  if (mathtype != NULL)
    *mathtype = FC_MT_UNKNOWN;
  
  // check input
  varSlot = _fc_getVarSlot(variable);
  if (varSlot == NULL || 
      (!numDataPoint && !numComponent && !assoc && !mathtype && !datatype)) { 
     fc_printfErrorMessage("%s, invalid variable", 
              fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Getting info for variable '%s'", varSlot->header.name);

  // fill in requested values
  if (datatype != NULL)
    *datatype = varSlot->datatype;
  if (assoc != NULL)
    *assoc = varSlot->assoc;
  if (numDataPoint != NULL)
    *numDataPoint = varSlot->numDataPoint;
  if (numComponent != NULL)
    *numComponent = varSlot->numComponent;
  if (mathtype != NULL)
    *mathtype = varSlot->mathtype;
  
  return FC_SUCCESS;
}

/**
 * \ingroup  Variable
 * \brief  Return the data type of a variable
 *
 * \description
 *
 *    Returns the type that the variable data is stored as (e.g.
 *    float, int, double, etc>
 *
 * \modifications  
 *     - Nancy Collins, Created.
 *     - 12/11/03 RM, added default return datatype = FC_DT_UNKNOWN, and 
 *               input check for datatype = NULL
 */
FC_ReturnCode fc_getVariableDataType(
  FC_Variable variable,        /**< input - variable  */
  FC_DataType *datatype        /**< output - type float, int, char, etc */
) {
  _FC_VarSlot* varSlot;
  
  // set up defaults
  if (datatype != NULL) 
    *datatype = FC_DT_UNKNOWN;
  
  // check input
  varSlot = _fc_getVarSlot(variable);
  if (varSlot == NULL || datatype == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Getting data type for variable '%s'", 
                      varSlot->header.name);

  // fill in requested value
  *datatype = varSlot->datatype;
  
  return FC_SUCCESS;
}

/**
 * \ingroup  Variable
 * \brief  Return data association of a variable.
 *
 * \description
 *  
 *    Return what subentity of the mesh (e.g. vertex, element, 
 *    face, edge) the data values are associated with.
 *
 * \modifications   
 *     - Nancy Collins, Created.
 *     - 12/16/03 RM, added default return datatype = FC_AT_UNKNOWN, and 
 *               input check for datatype = NULL
 */
FC_ReturnCode fc_getVariableAssociationType(
  FC_Variable variable,     /**< input - variable  */
  FC_AssociationType *assoc /**< output - data is per vertex, per element, 
                               etc */
) {
  _FC_VarSlot* varSlot;
  
  // set up default value
  if (assoc != NULL)
      *assoc = FC_AT_UNKNOWN; 
  
  // check input
  varSlot = _fc_getVarSlot(variable);
  if (varSlot == NULL || assoc == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
    
  // log message
  fc_printfLogMessage("Getting association type for variable '%s'", 
                      varSlot->header.name);

  // fill in requested value
  *assoc = varSlot->assoc;
  
  return FC_SUCCESS;
}

/**
 * \ingroup  Variable
 * \brief  Return number of components of a variable.
 *
 * \description 
 * 
 *    Return the number of components that make up a data point.
 *    This will always be 1 for scalar data, 2+ for vectors, etc.
 *    The total number of data values is numDataPoint * numComponent.
 *
 * \modifications  
 *     - Nancy Collins, Created.
 *     - 12/11/03 RM, added default return numComponent = -1, and input check 
 *               for numComponent = NULL
 */
FC_ReturnCode fc_getVariableNumComponent(
  FC_Variable variable,  /**< input - variable */
  int *numComponent      /**< output - number of components */
) {
  _FC_VarSlot* varSlot;
  
  // set up defaults
  if(numComponent != NULL)
    *numComponent = -1 ; 
  
  // check input
  varSlot = _fc_getVarSlot(variable);
  if (varSlot == NULL || numComponent == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Getting numComponent for variable '%s'", 
                      varSlot->header.name);

  // fill in requested value
  *numComponent = varSlot->numComponent;
  
  return FC_SUCCESS;
}

/**
 * \ingroup  Variable
 * \brief  Return the math type of a variable
 *
 * \description
 *  
 *    The math type of a data point describes the organization
 *    of the data point's components. E.g. FC_MT_VECTOR describes
 *    a vector organization of the components like X, Y & Z.
 *
 * \modifications  
 *     - Nancy Collins, Created.
 *     - 12/11/03 RM, added default return mathtype = FC_MT_UNKNOWN, and 
 *               input check for mathtype = NULL
 */
FC_ReturnCode fc_getVariableMathType(
  FC_Variable variable,  /**< input - variable */
  FC_MathType* mathtype  /**< output - the math type of the the data points */
) {
  _FC_VarSlot* varSlot;
  
  // set default values
  if (mathtype != NULL)
    *mathtype = FC_MT_UNKNOWN; 
  
  // check input
  varSlot = _fc_getVarSlot(variable);
  if (varSlot == NULL || mathtype == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Getting math type for variable '%s'", 
                      varSlot->header.name);

  // fill in requested value
  *mathtype = varSlot->mathtype;
  
  return FC_SUCCESS;
}

/**
 * \ingroup  Variable
 * \brief  Return number of data points in a variable.
 *
 * \description
 *
 *    The number of data points will be equal to the number of
 *    mesh subentities that the variable is associated with.
 *    NumComponent is the number of values per point. The total
 *    number of data values is numDataPoint * numComponent;
 *
 * \modifications  
 *     - Nancy Collins, Created.
 *     - 12/11/03 RM, added default return numDataPoint = -1, and input check 
 *               for numDataPoint = NULL
 */
FC_ReturnCode fc_getVariableNumDataPoint(
  FC_Variable variable,    /**< input - variable handle */
  int *numDataPoint        /**< output - count of data items */
) { 
  _FC_VarSlot* varSlot;
  
  // set default values
  if(numDataPoint != NULL)
    *numDataPoint = -1;
  
  // check input
  varSlot = _fc_getVarSlot(variable);
  if (varSlot == NULL || numDataPoint == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Getting numDataPoint for variable '%s'", 
                      varSlot->header.name);

  // fill in requested value
  *numDataPoint = varSlot->numDataPoint;
  
  return FC_SUCCESS;
}

//@}

/** \name Get variable big data. */
//-------------------------------
//@{

/**
 * \ingroup  Variable
 * \brief  Get the data in a variable.
 *
 * \description
 *  
 *    Return a pointer to the data array.  This must be treated as read-only.  
 *    If the user wants to change the data they must make a copy first.
 *
 *    Can also be used to get the data in a step of a sequence variable
 *    by calling with the variable handle for the desired step.
 *
 * \modifications  
 *    - Nancy Collins, Created.
 */
FC_ReturnCode fc_getVariableDataPtr(
 FC_Variable variable,  /**< input - variable  */
 void **data_p          /**< output - pointer to data array, treat as read only */
) {
  FC_ReturnCode rc;
  _FC_VarSlot* varSlot;
  
  // default
  if (data_p)
    *data_p = NULL;
  
  // check input
  varSlot = _fc_getVarSlot(variable);
  if (varSlot == NULL || data_p == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Getting data for variable '%s'", varSlot->header.name);

  // lazy fetching of data done here
  if (varSlot->data == NULL && varSlot->header.committed == 1) {
    rc = _fc_readVariableData(variable);
    if (rc != FC_SUCCESS) {
      fc_printfErrorMessage("Lazy loading of data for variable '%s' failed",
                            varSlot->header.name);
      return rc;
    }
  }

  // fill in requested value
  *data_p = varSlot->data;
  
  return FC_SUCCESS;
}

/**
 * \ingroup  Variable
 * \brief  Get a copy of a variable's data as the requested data type.
 *
 * \description
 *
 *    Returns the variable's data as the type requested by the user.
 *    The user is responsible for freeing the array. This routine casts
 *    to the requested user type so there is a potential for lost information
 *    (e.g. if you cast a floating point number to int it is truncated, not
 *    rounded).
 *
 *    If you are worried about memory or casting, you should use
 *    fc_getVariableDataPtr() to get a pointer to the raw data array.
 *
 *    Can also be used to get the data in a step of a sequence variable
 *    by calling with the variable handle for the desired step.
 *
 * \modifications  
 *    - 02/02/05 WSK created.
 */
FC_ReturnCode fc_getVariableDataAsDataType(
 FC_Variable variable,  /**< input - variable  */
 FC_DataType dataType,  /**< input - the data type for the returned data */
 void **data            /**< output - returned data */
) {
  FC_ReturnCode rc;
  int i;
  int numDataPoint, numComp, numTotal;
  FC_DataType dataType_orig;
  void *data_orig;
  
  // default
  if (data)
    *data = NULL;
  
  // check input
  if (!fc_isVariableValid(variable) || !fc_isDataTypeValid(dataType) || 
      dataType == FC_DT_UNKNOWN || data == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Getting data for variable '%s' as type %s", 
                      _fc_getVarSlot(variable)->header.name,
                      fc_getDataTypeText(dataType));

  // get variable info and data
  rc = fc_getVariableInfo(variable, &numDataPoint, &numComp, NULL, 
                          NULL, &dataType_orig);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_getVariableDataPtr(variable, &data_orig);
  if (rc != FC_SUCCESS)
    return rc;
  numTotal = numDataPoint*numComp;
  
  // allocate return array
  *data = malloc(numTotal*fc_sizeofDataType(dataType));
  if (!*data) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }

  // special case - the same data type
  if (dataType == dataType_orig) {
    memcpy(*data, data_orig, numTotal*fc_sizeofDataType(dataType));
    return FC_SUCCESS;
  }

  // Fun!
  if (dataType == FC_DT_CHAR) {
    if (dataType_orig == FC_DT_INT) 
      for (i = 0; i < numTotal; i++)
        ((char*)(*data))[i] = ((int*)data_orig)[i];
    else if (dataType_orig == FC_DT_FLOAT) 
      for (i = 0; i < numTotal; i++)
        ((char*)(*data))[i] = ((float*)data_orig)[i];
    else if (dataType_orig == FC_DT_DOUBLE) 
      for (i = 0; i < numTotal; i++)
        ((char*)(*data))[i] = ((double*)data_orig)[i];
  }
  else if (dataType == FC_DT_INT) {
    if (dataType_orig == FC_DT_CHAR) 
      for (i = 0; i < numTotal; i++)
        ((int*)(*data))[i] = ((char*)data_orig)[i];
    else if (dataType_orig == FC_DT_FLOAT) 
      for (i = 0; i < numTotal; i++)
        ((int*)(*data))[i] = ((float*)data_orig)[i];
    else if (dataType_orig == FC_DT_DOUBLE) 
      for (i = 0; i < numTotal; i++)
        ((int*)(*data))[i] = ((double*)data_orig)[i];
  }
  else if (dataType == FC_DT_FLOAT) {
    if (dataType_orig == FC_DT_CHAR) 
      for (i = 0; i < numTotal; i++)
        ((float*)(*data))[i] = ((char*)data_orig)[i];
    else if (dataType_orig == FC_DT_INT) 
      for (i = 0; i < numTotal; i++)
        ((float*)(*data))[i] = ((int*)data_orig)[i];
    else if (dataType_orig == FC_DT_DOUBLE) 
      for (i = 0; i < numTotal; i++)
        ((float*)(*data))[i] = ((double*)data_orig)[i];
  }
  else if (dataType == FC_DT_DOUBLE) {
    if (dataType_orig == FC_DT_CHAR) 
      for (i = 0; i < numTotal; i++)
        ((double*)(*data))[i] = ((char*)data_orig)[i];
    else if (dataType_orig == FC_DT_INT) 
      for (i = 0; i < numTotal; i++)
        ((double*)(*data))[i] = ((int*)data_orig)[i];
    else if (dataType_orig == FC_DT_FLOAT) 
      for (i = 0; i < numTotal; i++)
        ((double*)(*data))[i] = ((float*)data_orig)[i];
  }
  
  return FC_SUCCESS;
}


/**\ingroup Variable
 * \brief  In a sequence variable, get data for steps at a single data point
 *         returned as a sequence variable with association FC_AT_WHOLE_MESH
 *
 * \description
 *
 *    Given a sequence variable and the id of a particular data point
 *    (corresponds to a particular entity id for the variable), this will
 *    return a seqvar with data values for that point only with association
 *    FC_AT_WHOLE_MESH. This method calls fc_getVariableDataPointSlice and then
 *    packs the return array into a sequence Variable.
 *
 * \modifications  
 *    - 2/10/05 ACG, Created.
 */
FC_ReturnCode fc_getSeqVariableDataPointSliceAsSeqVariable(
 int numStep,         /**< input - number of steps in the sequence variable */ 
 FC_Variable* seqVar, /**< input - seq variable handles  */
 int dataPointID,     /**< input - ID of the data point/entity of interest */
 char* slicevarname,    /**< input - name of output seqvar */
 FC_Variable** sliceVar /**< output - output seqvar */
) {
  FC_ReturnCode rc;
  FC_Mesh mesh;
  FC_Sequence seq;

  int numComponent;
  FC_MathType mathtype;
  FC_DataType datatype;
  FC_AssociationType assoc;

  void* history;
  void* dcounter;
  int i;

  //only need to check slicevar, fc_getSeqVariableDataPointSlice checks the
  //rest

  if (sliceVar)
    *sliceVar = NULL;

  if (sliceVar == NULL || slicevarname == NULL){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  };

  //get data & its does input checks
  rc = fc_getSeqVariableDataPointSlice(numStep,seqVar, dataPointID, &history);
  if (rc != FC_SUCCESS){
    return rc;
  }


  // Get info about seq variable
  rc = fc_getVariableInfo(seqVar[0], NULL, &numComponent, &assoc, &mathtype,
                          &datatype);
  if (rc != FC_SUCCESS){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(rc));
    free(history);
    return rc;
  }

  //make slice var
  rc = fc_getMeshFromVariable(seqVar[0],&mesh);
  rc = fc_getSequenceFromSeqVariable(numStep,seqVar,&seq);
  if (rc != FC_SUCCESS){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(rc));
    free(history);
    return rc;
  }

  rc = fc_createSeqVariable(mesh, seq, slicevarname, &numStep,sliceVar);
  if (rc != FC_SUCCESS){
      fc_deleteSeqVariable(numStep,*sliceVar);
      free(*sliceVar);
      *sliceVar = NULL;
      free(history);
      return rc;
  }

  for (i = 0; i < numStep; i++){
    switch (datatype){
    case FC_DT_INT:
      dcounter = &(((int*)history)[i*numComponent]);
      break;
    case FC_DT_FLOAT:
      dcounter = &(((float*)history)[i*numComponent]);
      break;
    case FC_DT_DOUBLE:
      dcounter = &(((double*)history)[i*numComponent]);
      break;
    case FC_DT_CHAR:
      dcounter = &(((char*)history)[i*numComponent]);
      break;
    default:
      printf("** Developer Alert! Should never reach this point!**\n");
      fflush(NULL);
      fc_printfErrorMessage("%s", "Invalid data type");
      return FC_ERROR;
      break;
    }

    rc = fc_setVariableData((*sliceVar)[i],1,numComponent,
                            FC_AT_WHOLE_MESH,mathtype,
                            datatype,dcounter);
    if (rc != FC_SUCCESS){
      fc_deleteSeqVariable(numStep,*sliceVar);
      free(*sliceVar);
      *sliceVar = NULL;
      free(history);
      return rc;
    }
  }

  dcounter = NULL;
  free(history);
  return FC_SUCCESS;
}


/**
 * \ingroup  Variable
 * \brief  In a sequence variable, get data for steps at a single data.
 *
 * \description
 *
 *    Given a sequence variable and the id of a particular data point
 *    (corresponds to a particular entity id for the variable), this will
 *    return an array of data values at that point for each step in the
 *    sequence. This type of information for time series is often called a
 *    point history. If the data has more than one component, the array will be
 *    interleaved; that is, for each step all components will be listed before
 *    the next step (e.g.  Step1Comp1 Step2Comp2 .. Step1CompN ... StepMComp1
 *    StepMCompN).
 *
 *    This routine returns an allocated data array.  The user is responsible 
 *    for freeing it.
 *
 * \modifications  
 *    - 10/21/04 WSK, Created.
 */
FC_ReturnCode fc_getSeqVariableDataPointSlice(
 int numStep,         /**< input - number of steps in the sequence variable */ 
 FC_Variable* seqVar, /**< input - seq variable handles  */
 int dataPointID,     /**< input - ID of the data point/entity of interest */
 void **history       /**< output - data array */
) {
  FC_ReturnCode rc;
  int i;
  _FC_VarSlot* varSlot;
  size_t size;
  int numDataPoint, numComp;
  FC_DataType dataType;
  void* temp_data, *data_p;

  // default
  if (history)
    *history = NULL;
  
  // check input
  if (!fc_isSeqVariableValid(numStep, seqVar) || dataPointID < 0 ||
      !history) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // Get info about seq variable
  rc = fc_getVariableInfo(seqVar[0], &numDataPoint, &numComp, NULL, NULL,
                          &dataType);
  if (rc != FC_SUCCESS)
    return rc;

  // a little more checking
  if (dataPointID >= numDataPoint) {
    fc_printfErrorMessage("Data point ID is outside range");
    return FC_INPUT_ERROR;
  }

  // log message
  varSlot = _fc_getVarSlot(seqVar[0]);
  fc_printfLogMessage("Getting history for seq variable '%s'", 
                      varSlot->header.name);

  // Collect history data
  size = numComp * fc_sizeofDataType(dataType);
  temp_data = malloc(numStep*size);
  if (temp_data == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  for (i = 0; i < numStep; i++) {
    rc = fc_getVariableDataPtr(seqVar[i], &data_p);
    if (rc != FC_SUCCESS) {
      free(temp_data);
      return rc;
    }
    memcpy((char*)temp_data + i*size, (char*)data_p + dataPointID*size, size);
  }

  // set return   
  *history = temp_data;

  return FC_SUCCESS;
}

//@}

/** \name Print variable. */
//-------------------------------
//@{

/**
 * \ingroup  Variable
 * \brief Print variable meta data to stdout, and optionally the data values.
 *
 * \modifications 
 *    - Nancy Collins  Created.
 *    - 11/20/2003 WSK Fixed up. Want to replace a bunch of
 *                     specific print routines (like printSurfaceNormals).
 *    - 10/31/2007 CDU Added flag to allow printing of reduced number of digits
 */
FC_ReturnCode fc_printVariable(
  FC_Variable variable, /**< input - data variable */
  char* label,          /**< input - label to prepend variable name with, 
                                     can be NULL */
  int print_data        /**< input - 1 = print the data values. 
                                     2 = print the data values in exponent form with 5 digits of precision */
) {
  FC_ReturnCode rc;
  int i, j;
  _FC_VarSlot* varSlot;
  int numDataPoint, numComp;
  FC_AssociationType assoc;
  FC_MathType mathType;
  FC_DataType dataType;
  void* dbuf;

  // special case: printing a null handle is not an error
  if (FC_HANDLE_EQUIV(variable, FC_NULL_VARIABLE)) {
    printf("Variable: FC_NULL_VARIABLE\n");
    fflush(NULL);
    return FC_SUCCESS;
  }

  // check input
  varSlot = _fc_getVarSlot(variable);
  if (varSlot == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // log message
  fc_printfLogMessage("Printing variable '%s'", varSlot->header.name);

  // print label & var name
  if (label)
    printf("%s: ", label);
  else
    printf("Variable: ");
  printf("'%s'\n", varSlot->header.name);

  // print meta data
  rc = fc_getVariableInfo(variable, &numDataPoint, &numComp, &assoc, 
                            &mathType, &dataType);
  if (rc != FC_SUCCESS) {
    fc_printfErrorMessage("Failed to get info for variable '%s'",
                          varSlot->header.name);
    return rc;
  }
  printf("%snumDataPoint = %d, numComponent = %d\n", INDENT_STR, 
         numDataPoint, numComp);
  printf("%sassoc = %s, mathType = %s, dataType = %s\n", INDENT_STR, 
         fc_getAssociationTypeText(assoc), fc_getMathTypeText(mathType), 
         fc_getDataTypeText(dataType));

  // print values
  if (print_data) {
    rc = fc_getVariableDataPtr(variable, &dbuf);
    if (rc != FC_SUCCESS) {
      fc_printfErrorMessage("Failed to get data for variable '%s'",
                            varSlot->header.name);
      return FC_ERROR;
    }
    
    // print the data 
    printf("%sData:\n", INDENT_STR);
    switch (dataType) {
    case FC_DT_CHAR: 
      {
        char* cast_buf = (char *)dbuf;
        for (i = 0; i < numDataPoint; i++) {
          printf("%s%s%d:  ", INDENT_STR, INDENT_STR, i);
          for (j = 0; j < numComp; j++)
            printf("%c ", cast_buf[i*numComp + j]);
          printf("\n");
        }
      }
      break;
    case FC_DT_INT:
      {
        int* cast_buf = (int *)dbuf;
        for (i = 0; i < numDataPoint; i++) {
          printf("%s%s%d:  ", INDENT_STR, INDENT_STR, i);
          for (j = 0; j < numComp; j++)
            printf("%6d ", cast_buf[i*numComp + j]);
          printf("\n");
        }
      }
      break;
    case FC_DT_FLOAT:
      {
        float* cast_buf = (float *)dbuf;
        for (i = 0; i < numDataPoint; i++) {
          printf("%s%s%d:  ", INDENT_STR, INDENT_STR, i);
          for (j = 0; j < numComp; j++)
            printf("%11.6g ", cast_buf[i*numComp + j]);
          printf("\n");
        }
      }
      break;
    case FC_DT_DOUBLE:
      {
        double* cast_buf = (double *)dbuf;
        for (i = 0; i < numDataPoint; i++) {
          printf("%s%s%d:  ", INDENT_STR, INDENT_STR, i);
          for (j = 0; j < numComp; j++)
            //Should be 15 significant digits but there were
            //differences between Mac and Linux (sun same as Mac)
            //these may be due to rounding differences in printf()
            //printf("%20.15g ", cast_buf[i*numComp + j]);
	    switch(print_data){
	    case 2:  printf("%.5le ",   cast_buf[i*numComp + j]); break;
	    default: printf("%16.11g ", cast_buf[i*numComp + j]); break;
	    }
          printf("\n");
        }
      }
      break;
    case FC_DT_UNKNOWN:
      printf("Cannot print data type '%s'\n", fc_getDataTypeText(dataType));
    }
  }

  return FC_SUCCESS;
}

/**
 * \ingroup  Variable
 * \brief Print sequence variable meta data to stdout, and optionally the 
 *        data values.
 *
 * \modifications 
 *    - 11/20/03 WSK Created.
 *    - 07/23/2007 CDU Turned into a stub that calls fc_printSeqVariableRange
 */
FC_ReturnCode fc_printSeqVariable(
  int numStep,        /**< input - number of steps in the sequence variable */
  FC_Variable* seqVars, /**< input - the steps of the sequence variable */
  char* label,          /**< input - label to prepend variable name with, 
                                     can be NULL */
  int print_data        /**< input - 1 = print the data values. */
) {
  return fc_printSeqVariableRange(numStep, seqVars, -1, -2, label, print_data);
}
/**
 * \ingroup  Variable
 * \brief Print sequence variable meta data to stdout, and optionally the 
 *        data values (with start and stop index values)
 *
 * \modifications 
 *    - 11/20/03   WSK Created.
 *    - 07/23/2007 CDU Moved guts of printSeqVariable to here and added ranges
 *    - 10/31/2007 CDU Added flag to allow printing of reduced number of digits
 */
FC_ReturnCode fc_printSeqVariableRange(
  int numStep,        /**< input - number of steps in the sequence variable */
  FC_Variable* seqVars, /**< input - the steps of the sequence variable */
  int range_start,      /**< input - offset of first variable to print  
                                      or -1 for first, -2 for last */
  int range_stop,       /**< input - offset of last variable to print 
                                      or -1 for first, -2 for last */
  char* label,          /**< input - label to prepend variable name with, 
                                     can be NULL */
  int print_data        /**< input - 1 = print the data values.
                                     2 = print the data values in exponent form with 5 digits of precision */
) {
  FC_ReturnCode rc;
  int i, j, k;
  _FC_VarSlot* varSlots[numStep];
  int numDataPoint, numComp;
  FC_AssociationType assoc;
  FC_MathType mathType;
  FC_DataType dataType;
  void* dbuf;
  int start, stop, mod;

  // special case: printing a null list is not an error
  if (numStep <= 0 && seqVars == NULL) {
    fc_printfWarningMessage("Printing a NULL sequence variable");
    printf("Sequence Variable: FC_NULL_VARIABLE\n");
    return FC_SUCCESS;
  }

  // check input
  for (i = 0; i < numStep; i++) {
    varSlots[i] = _fc_getVarSlot(seqVars[i]);
    if (varSlots[i] == NULL) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
      return FC_INPUT_ERROR;
    }
  }

  // log message
  fc_printfLogMessage("Printing variable '%s'", varSlots[0]->header.name);

  // print label & var name
  if (label)
    printf("%s: ", label);
  else
    printf("Sequence Variable: ");
  printf("'%s'\n", varSlots[0]->header.name);

  // print meta data
  rc = fc_getVariableInfo(seqVars[0], &numDataPoint, &numComp, &assoc, 
                            &mathType, &dataType);
  if (rc != FC_SUCCESS) {
    fc_printfErrorMessage("Failed to get info for seq variable '%s'",
                          varSlots[0]->header.name);
    return rc;
  }
  printf("%snumDataPoint = %d, numComponent = %d, numStep = %d\n", INDENT_STR, 
         numDataPoint, numComp, numStep);
  printf("%sassoc = %s, mathType = %s, dataType = %s\n", INDENT_STR, 
         fc_getAssociationTypeText(assoc), fc_getMathTypeText(mathType), 
         fc_getDataTypeText(dataType));


  // print values
  if ((print_data) && (numStep)) {

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

      rc = fc_getVariableDataPtr(seqVars[k], &dbuf);
      if (rc != FC_SUCCESS) {
        fc_printfErrorMessage("Failed to get data for step %d of seq variable"
                              " '%s'", k, varSlots[0]->header.name);
        return FC_ERROR;
      }
      
      // print the data 
      printf("%sData (Step %d of seqVar '%s'):\n", INDENT_STR, 
             k, varSlots[0]->header.name);
      switch (dataType) {
      case FC_DT_CHAR: 
        {
          char* cast_buf = (char *)dbuf;
          for (i = 0; i < numDataPoint; i++) {
            printf("%s%s%d:  ", INDENT_STR, INDENT_STR, i);
            for (j = 0; j < numComp; j++)
              printf("%c ", cast_buf[i*numComp + j]);
            printf("\n");
          }
        }
        break;
      case FC_DT_INT:
        {
          int* cast_buf = (int *)dbuf;
          for (i = 0; i < numDataPoint; i++) {
            printf("%s%s%d:  ", INDENT_STR, INDENT_STR, i);
            for (j = 0; j < numComp; j++)
              printf("%6d ", cast_buf[i*numComp + j]);
            printf("\n");
          }
        }
        break;
      case FC_DT_FLOAT:
        {
          float* cast_buf = (float *)dbuf;
          for (i = 0; i < numDataPoint; i++) {
            printf("%s%s%d:  ", INDENT_STR, INDENT_STR, i);
            for (j = 0; j < numComp; j++)
              printf("%11.6g ", cast_buf[i*numComp + j]);
            printf("\n");
          }
        }
        break;
      case FC_DT_DOUBLE:
        {
          double* cast_buf = (double *)dbuf;
          for (i = 0; i < numDataPoint; i++) {
            printf("%s%s%d:  ", INDENT_STR, INDENT_STR, i);
            for (j = 0; j < numComp; j++)
              // FIX? maybe %20.15g so always align
              //printf("%11.15g ", cast_buf[i*numComp + j]);
	      switch(print_data){
	      case 2:  printf("%.5le ",   cast_buf[i*numComp + j]); break;
	      default: printf("%11.15g ", cast_buf[i*numComp + j]); break;
	      }
            printf("\n");
          }
        }
        break;
      case FC_DT_UNKNOWN:
        printf("Cannot print data type '%s'\n", fc_getDataTypeText(dataType));
      }
    }
  }

  return FC_SUCCESS;
}

//@}

/**
 * \ingroup  PrivateVariable
 * \brief  Initialize a variable slot.
 *
 * \modifications  
 *   APR-30-2003  W Koegler  Created. 
 */
void _fc_initVarSlot(
  _FC_VarSlot* varSlot    /**< input - the vtable slot to initialize */
) {
  if (varSlot != NULL) {
    _fc_initSlotHeader(&varSlot->header);
    varSlot->ds = FC_NULL_DATASET;
    varSlot->mesh = FC_NULL_MESH;
    varSlot->sequence = FC_NULL_SEQUENCE;
    varSlot->stepID = -1;
    varSlot->numDataPoint = 0;
    varSlot->numComponent = 0;
    varSlot->mathtype = FC_MT_UNKNOWN;
    varSlot->datatype = FC_DT_UNKNOWN;
    varSlot->assoc = FC_AT_UNKNOWN;
    varSlot->data = NULL;
  }   
}

/**
 * \ingroup  PrivateVariable
 * \brief  Release non-necessary resources in a variable slot.
 *
 * \description
 *  
 *    Releases the data values.
 *
 * \modifications  
 *   Aug 6, 2002  W Koegler  Created
 */
void _fc_releaseVarSlot(
  _FC_VarSlot* varSlotp       /**< input - vtable pointer */
) {
    // release big data
    if (varSlotp->data) {
      free(varSlotp->data);
      varSlotp->data = NULL;
    }

    return;
}

/**
 * \ingroup  PrivateVariable
 * \brief   Clear a variable slot.
 *
 * \description
 *
 *    Releases all dynamically allocated resources in a variable,
 *    and reinitializes all members.
 *
 * \modifications  
 *   - Aug 6, 2002  W Koegler  Created
 *   - 2003-APR-31  W Koegler  Made more comprehensive
 */
void _fc_clearVarSlot(
 _FC_VarSlot* varSlotp       /**< input - vtable pointer */
) {
  if (varSlotp != NULL) {
    // free big data
    _fc_releaseVarSlot(varSlotp);

    // free other dynamic arrays, if any

    // clear table header
    _fc_clearSlotHeader(&varSlotp->header);

    // reinit
    _fc_initVarSlot(varSlotp);
  }
}

/**
 * \ingroup  PrivateVariable
 * \brief  Get the size of the variable table.
 *
 * \modifications  
 *   - 05/21/04 WSK Created.
 */
int _fc_getVarTableSize() {
  return varTableSize;
}

/**
 * \ingroup  PrivateVariable
 * \brief   Get an empty variable slot from the variable table.
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
_FC_VarSlot* _fc_getNewVarSlot() {
  _FC_VarSlot* varSlot;

  varSlot = (_FC_VarSlot*) _fc_getNewSlot(&varTableSize, 
					  (_FC_SlotHeader***)&varTable,
					  &varOpenSlots, sizeof(_FC_VarSlot));
  _fc_initVarSlot(varSlot);

  return varSlot;
}

/**
 * \ingroup  PrivateVariable
 *
 * \brief Get a variable slot from a variable handle.
 *
 * \description
 * 
 *    The slot id of the variable handle is used to find the
 *    the slot in the variable table. The identity is then
 *    verified by checking that the slot and the handle's
 *    uIDs are the same. On success, a pointer to the slot 
 *    is returned. Otherwise a NULL pointer is returned.
 *
 * \modifications 
 *   - 2003-OCT-19  WSK  Created to replace using _fc_checkhandle and GET_ID 
 *      for better type checking (using pointer's to slots is less bug
 *      prone than lot's of xxxTable[xxxID] references).
 *   - 05/25/04 WSK Changed to call new _fc_isHandleValid common to all tables.
 */
_FC_VarSlot* _fc_getVarSlot(
  FC_Variable variable    /**< input - variable */
) {
  if (_fc_isHandleValid(variable.slotID, variable.uID, varTableSize,
                        (_FC_SlotHeader**)varTable))
    return varTable[variable.slotID];
  else
    return NULL;
}


/**
 * \ingroup  PrivateVariable
 *
 * \brief Get a variable slot based on ID.
 *
 * \description
 * 
 *    The _FC_VarSlot at the requested slotID is returned. This should
 *    only be used when you really know that it is the slot you want
 *    because there is no checking of the uID like with _fc_getVarSlot.
 *    The only checking is that the slotID is exists in the table and
 *    that the uID is > 0 (i.e. it is not empty).
 *
 * \modifications 
 *    - 11/27/03 WSK Created so that external routines can index into
 *      tables without knowing about them.
 *    - 05/25/04 WSK Changed to call new _fc_getSlot() common to all tables.
 */
_FC_VarSlot* _fc_getVarSlotFromID(
  int varID
) {
  return (_FC_VarSlot*) _fc_getSlot(varID, varTableSize, 
                                   (_FC_SlotHeader**)varTable);
}

/**
 * \ingroup  PrivateVariable
 * \brief Delete the variable slot associated with the variable handle.
 *
 * \description 
 *    
 *    This routine clears the dynamic data in the slot and then
 *    deletes the slot from the table.
 *
 * \modifications
 *   - 12/20/2005 WSD. Created.
 */
FC_ReturnCode _fc_deleteVarSlot(
  FC_Variable variable    /**< input - variable */
) {
  _FC_VarSlot* varSlot;
  
  varSlot = _fc_getVarSlot(variable);
  if (varSlot == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // do it
  _fc_clearVarSlot(varSlot);
  return _fc_deleteSlot(varSlot->header.slotID, varTableSize, 
			(_FC_SlotHeader**)varTable, &varOpenSlots);
}

/**
 * \ingroup  PrivateVariable
 * \brief  Print the contents of the Variable Table (varTable) to stderr.
 *
 * \description
 *
 *    Prints the contents of the Variable Table in a human readable form.
 *    It takes a string to use as a label for this print out of the table.
 *    Pass NULL if you don't want a label.
 *
 * \modifications 
 *    - Nancy Collins Created.
 *    - 2003-NOV-13  WSK  Fixed up.
 */
void _fc_printVarTable(
  char *label     /**< Input - label for this version the table (without a
                           trailing \\n) pass NULL if no label is desired */
) {
  int i;
  _FC_VarSlot* varSlot;
  
  // print table heading
  fprintf(stderr, TABLE_DIVIDER_STRING);
  fprintf(stderr, "Variable Table (varTable):\n");
  if (label != NULL)
    fprintf(stderr, "%s\n", label);
  fprintf(stderr, TABLE_DIVIDER_STRING);
  fflush(NULL);
  
  // We're done if it's empty
  if (varTableSize < 1) {
    fprintf(stderr, "(empty)\n");
    fprintf(stderr, "\n");
    fflush(NULL);
    return;
  }

  // print contents for each slot
  for (i = 0; i < varTableSize; i++) {
    varSlot = varTable[i];
    fprintf(stderr, "%3d: %s\n", i, (varSlot) ? "exists" : "NULL");
    if (varSlot == NULL)
      continue; 
    _fc_printSlotHeader(varSlot->header);
    fprintf(stderr, "     FC_Dataset = { %d, %d }, FC_Mesh = { %d, %d }\n",
            varSlot->ds.slotID, varSlot->ds.uID,
	    varSlot->mesh.slotID, varSlot->mesh.uID);
    fprintf(stderr, "     FC_Sequence = { %d, %d }, stepID = %d\n",
            varSlot->sequence.slotID, varSlot->sequence.uID, varSlot->stepID);
    fprintf(stderr, "     numDataPoint = %d, numComponent = %d\n", 
            varSlot->numDataPoint, varSlot->numComponent);
    fprintf(stderr, "     association = %s, mathtype = %s\n",
            fc_getAssociationTypeText(varSlot->assoc),
            fc_getMathTypeText(varSlot->mathtype));
    fprintf(stderr, "     datatype = %s, data = %s\n", 
            fc_getDataTypeText(varSlot->datatype), 
            (varSlot->data) ? "exists" : "NULL");
  }
  fprintf(stderr, "\n");
  fflush(NULL);
  
  return;
}

/**
 * \ingroup  PrivateVariable
 * \brief  Free all entries in the variable table
 *
 * \description
 *
 *    This should only be called after making sure that
 *    all slots have been cleared.
 *
 * \modifications  
 *   - 05/21/04 WSK Created.
 */
void _fc_freeVarTable() {
  int i;
  for (i = 0; i < varTableSize; i++) {
    _fc_clearVarSlot(varTable[i]);
    free(varTable[i]);
  }
  free(varTable);
  fc_freeSortedIntArray(&varOpenSlots);
  varTableSize = 0;
  varTable = NULL;
}


/**
 * \ingroup  PrivateVariable
 * \brief  Get a common component extension for combining vars
 *
 * \description
 *
 *    This function is used when searching for single-component
 *    vars/seqVars that should be combined into a multicomponent
 *    var/seqVar. This function holds a table of common extensions
 *    and is designed to generate several permutations of 
 *    the extension format to make guessing the component names
 *    easier. A state variable is passed to the function in order
 *    to keep track of where the search is. The idea is that a user
 *    keeps calling this function until a good match is found or
 *    the function runs out of guesses.
 *
 *    This function produces one component extension at a time,
 *    sequentially (ie, for the first extension guess, it gives
 *    _x, _y, and then _z). The lastCompInGuess variable indicates
 *    that the result is the last component of a particular
 *    extension. The lastCompInSearch indicates that the result
 *    is the last component of the last possible guess.
 *
 *    The user can find out info about the guesses by passing 
 *    NULL into state. In this case, lastCompInGuess returns the
 *    maximum number of characters in an individual extension,
 *    while lastCompInSearch returns the largest number of 
 *    components a guess can have.
 *
 * \modifications  
 *   - 12/12/2007 CDU Created.
 */
FC_ReturnCode _fc_getNextComponentExtension(
  	               unsigned int *state,             /**< input - state information retained between calls */
		       unsigned int *lastCompInGuess,   /**< output - flag indicating this is the last component of a guess */
		       unsigned int *lastCompInSearch,  /**< output - flag indicating this is the last component of last guess */
		       char         *extension          /**< input - pre-allocated array where extension is written */ 
){



  //Extension Table: If you want to add more extension labels to the below
  //table, make sure you do the following:
  // - Update num_rows to reflect number of different extension patterns
  // - Update num_cols to reflect max number of components in a pattern
  // - Update max_ext_len if your extension symbols become longer
  // - Insert longer vectors first: Make sure tuvw comes before
  //    tuv so longer number of components are found first
  const unsigned int num_rows    = 6; //How many rows in table below
  const unsigned int num_cols    = 4; //How many columns in table below
  const unsigned int num_rounds  = 6; //All underline/capital letter cases
  const unsigned int max_ext_len = 2+2+1; //underscores + extension + \0
  
  char *vector_endings[][4] = { {"x",  "y",  "z",  NULL },
				{"x",  "y",  NULL, NULL },
				{"c0", "c1", "c2", "c3" },
				{"c0", "c1", "c2", NULL },
				{"t",  "u",  "v",  "w"  },
				{"t",  "u",  "v",  NULL } };
  

  char *underscore, *s;
  int upper;
  unsigned int round, row, col;


  //Always must have at least these two return pointers
  if((!lastCompInGuess)||(!lastCompInSearch)){
    fc_printfErrorMessage("Null variable in input");
    return FC_INPUT_ERROR;
  }
  
  //When user passes NULL in for state, they're actually querying to
  //see what the dimensions of the table are. Just pass these back
  if(!state){
    *lastCompInGuess  = max_ext_len;
    *lastCompInSearch = num_cols;
    return FC_SUCCESS;
  }

  //Unpack the state
  col   = (*state>> 0) & 0x0FF;
  row   = (*state>> 8) & 0x0FF;
  round = (*state>>16) & 0x0FF;
  
  //printf("Working on round/row/col : %d/%d/%d\n", round, row, col);
  
  //Check bounds
  if( (round >= num_rounds) || (row >= num_rows) || (col >= num_cols) || 
      (!vector_endings[row][col])){
    fc_printfErrorMessage("Input boundary check failure");
    return FC_INPUT_ERROR;
  }


  //The user may not want to actually retrieve the extension. They may
  //just be looping through the components in a row because one of the
  //components wasn't found.
  if(extension){

    switch(round){
    case 0: underscore="_";  upper=0; break;
    case 1: underscore="_";  upper=1; break;
    case 2: underscore="__"; upper=0; break;
    case 3: underscore="__"; upper=1; break;
    case 4: underscore="";   upper=1; break;
    case 5: underscore="";   upper=0; break;
    default:
      fc_printfErrorMessage("Bad case statement?");
      return FC_ERROR;
    }
    
    snprintf(extension, max_ext_len, "%s%s",underscore,vector_endings[row][col]);
    if(upper){ //Convert to uppercase
      s=extension;
      while(*s!='\0'){
	*s = toupper(*s);
	s++;
      }
    }
  }

  //Done with this entry: Update the state
  col++;
  *lastCompInGuess=0;
  *lastCompInSearch=0;
  if((col==num_cols) || (!vector_endings[row][col])){
    //End of column
    *lastCompInGuess=1;
    col=0;
    row++;
    if(row==num_rows){
      //End of rows
      row=0;
      round++;
      if(round==num_rounds){
	//End of search
	*lastCompInSearch=1;
      }
    }
  }

  //Pack the state again
  *state = (round<<16) | (row<<8) | (col);

  return FC_SUCCESS;
}
