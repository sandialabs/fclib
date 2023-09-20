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
 * \file threshold.c
 * \brief Implementation for \ref Threshold module.
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/modules/threshold.c,v $
 * $Revision: 1.69 $ 
 * $Date: 2006/09/23 02:20:15 $
 */

// C library includes
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// fc library dependencies
#include "base.h"
#include "storage.h"
#include "library.h"
#include "variable.h"
#include "variableP.h"
#include "subset.h"
#include "util.h"

// this module
#include "threshold.h"

/**
 * \addtogroup Threshold
 * \brief Apply thresholds to variables.
 */

/**
 * \ingroup  Threshold
 * \brief Create thresholded region.
 *
 * \description
 * 
 *     Takes a variable and a reference value, and creates a subset of entities
 *     whose data values pass the threshold test given by the comparison
 *     string. Possible
 *     values of this string: <b> ">","<",">=","<=","=","==","=>","=<" </b> .
 *
 * \modifications 
 *    - Nancy Collins  Created.
 *    - 11/26/03 WSK added return of association type & cleaned up
 *    - 04/12/04 RM, changed input int compare, to string (one of "<",">"
 *          "<=", ">=","=<","=>","=","=="). Changed comparison section to do 
 *          all comparisons on doubles.    
 *    - 06/15/04 RM changed output to return subset instead of marked list
 *    = 11/8/04 WSK simplified internals. Added '!=' case.
 */
FC_ReturnCode fc_createThresholdSubset(
  FC_Variable variable,  /**< input - variable */
  char* compare_str,     /**< input - compare string (eg. "<","<=", etc.) */
  double ref_value,      /**< input - reference data value */
  char* subsetName,      /**< input - new subset's name */
  FC_Subset *subset      /**< output - subset with thresholded elements */
) {  
  FC_ReturnCode rc;
  int i;
  int numDataPoint, numComp;
  void *datap;
  double epsilon, minimum; 
  FC_DataType datatype;
  FC_AssociationType assoc;
  FC_Mesh mesh;
  // comparisons
  int (*compare)(double, double, double, double);

  // default return
  if (subset) 
    *subset = FC_NULL_SUBSET; 
   
  // check input (can't be a global var)
  if (!fc_isVariableValid(variable) || fc_isVariableGlobal(variable) ||
      !compare_str || !subsetName || !subset) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // determine comparison
  if (!strcmp(compare_str, "=") || !strcmp(compare_str, "=="))
    compare = fc_eq;
  else if (!strcmp(compare_str, "!=") || !strcmp(compare_str, "=!"))
    compare = fc_neq;
  else if (!strcmp(compare_str, "<"))
    compare = fc_lt;
  else if (!strcmp(compare_str, "<=") || !strcmp(compare_str, "=<"))
    compare = fc_lteq;
  else if (!strcmp(compare_str, ">"))
    compare = fc_gt;
  else if (!strcmp(compare_str, ">=") || !strcmp(compare_str, "=>"))
    compare = fc_gteq;
  else {
    fc_printfErrorMessage("Unknown comparison string: '%s'", compare_str);
    return FC_INPUT_ERROR;
  }

  // log message
  fc_printfLogMessage("Creating new threshold subset '%s'", subsetName);

  // get variable info
  rc = fc_getVariableInfo(variable, &numDataPoint, &numComp, &assoc, NULL, 
			  &datatype);
  if (rc != FC_SUCCESS)
    return rc;

  // more checking
  if (numComp > 1) {
	fc_printfErrorMessage("can only do variables with 1 component.");
    return FC_ERROR;
  }
  if (datatype == FC_DT_INT || datatype == FC_DT_DOUBLE) {
    epsilon = DBL_EPSILON;
    minimum = DBL_MIN;
  }
  else if (datatype == FC_DT_FLOAT) {
    epsilon = FLT_EPSILON;
    minimum = FLT_MIN;
  }
  else {
    fc_printfErrorMessage("Can't do data type '%s'.",
			  fc_getDataTypeText(datatype));
    return FC_ERROR;
  }
	 
  rc = fc_getVariableDataPtr(variable, &datap);
  if (rc != FC_SUCCESS) {
    fc_printfErrorMessage("error getting data.");
    return rc;
  }  
  
  // get mesh that the variable belongs to
  rc = fc_getMeshFromVariable(variable, &mesh);
  if (rc != FC_SUCCESS)
    return rc;
  // create new subset
  rc = fc_createSubset(mesh, subsetName, assoc, &(*subset));
  if (rc != FC_SUCCESS)
    return rc;

  // determine subset members
  for (i = 0; i < numDataPoint; i++) {
    double value;
    if (datatype == FC_DT_INT)
      value = ((int*)datap)[i];
    else if (datatype == FC_DT_FLOAT)
      value = ((float*)datap)[i];
    else if (datatype == FC_DT_DOUBLE)
      value = ((double*)datap)[i];
    if (compare(value, ref_value, epsilon, minimum))
      fc_addMemberToSubset(*subset, i);
  }

  return FC_SUCCESS;
}


/**
 * \ingroup  Threshold
 * \brief Create subset according to range conditions.
 *
 * \description
 * 
 *     Creates a subset of the entities that pass the test specified
 *     by the given range condition; that is the reference variable
 *     satisfies two comparisons connected by
 *     a logical operator. Comparisons are defined by a string,
 *     with possible values: <b> compare = ">",">=",
 *     "==","<","<=", "!=" </b>. The operation string connecting
 *     the two comparisons can be: <b> operation = "OR"="or"="||",
 *     "AND"="and"="&&"="&","XOR"="xor"="^" (exclusive or) </b>.
 *
 *     Example:  to get the range ((var >= a) && (var < b)) 
 *       - fc_markThresholdedRange(var, ">=", a, "&&", "<", "b", &subsetName, 
 *         &subset);
 *
 * \modifications 
 *    - Nancy Collins  Created.
 *    - 11/26/03 WSK added return of association type & cleaned up
 *    - 04/12/04 RM, rewritten, added comparison strings to input argument 
 *         list and an operation string to link both comparisons.
 *         f.e. we can now check condition:
 *          variable "<" reference1  &&  variable ">=" reference2.
 *    - 06/25/04 RM changed to operate on subsets and return subset 
 *               instead of marked list  
 */
FC_ReturnCode fc_createThresholdRangeSubset(
  FC_Variable variable,/**< input - variable */
  char* compare1,    /**< input - comparison (string eg.">") with reference1*/
  double ref_value1, /**< input - reference data value                      */
  char* operation,   /**< input - operation that links the two comparisons
			eg. "OR","AND","XOR", "||","&&" etc.       */
  char* compare2,    /**< input - comparison (eg."<=") with reference2      */
  double ref_value2, /**< input - reference data value                      */
  char* subsetName,  /**< input - new subset's name                       */
  FC_Subset *subset  /**< output - subset with thresholded elements       */
) {
  FC_ReturnCode rc ;
  FC_Subset subset1;
  FC_Subset subset2;
  char *subset1Name = "temp_subset", *subset2Name ="temp_subset";

  // default return
  if (subset) 
    *subset = FC_NULL_SUBSET; 
   
  // check input
  if (!fc_isVariableValid(variable) || fc_isVariableGlobal(variable) ||
      !compare1 || !operation || !compare2 || !subsetName || !subset) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // log message
  fc_printfLogMessage("Creating new threshold range subset '%s'",
		      subsetName);

  // for each comparison call fc_getThreshSubset() function, creating new
  // subsets, then perform operation on the two subsets  
  rc = fc_createThresholdSubset(variable, compare1, ref_value1, subset1Name,
			     &subset1);
  if (rc != FC_SUCCESS) 
    return rc;
  
  rc = fc_createThresholdSubset(variable, compare2, ref_value2, subset2Name, 
			     &subset2);
  if (rc != FC_SUCCESS) {
    fc_deleteSubset(subset1);
    return rc;
  }

  // perform desired operation on the two subsets
  rc = fc_createSubsetIntersection(subset1, operation, subset2, subsetName,
       subset);

  // delete temporary subsets
  fc_deleteSubset(subset1);
  fc_deleteSubset(subset2);

  // report error of intersection operation
  if (rc != FC_SUCCESS) {
    fc_printfErrorMessage("error creating subset intersection.");
    return rc;
  }  
 
  return FC_SUCCESS;
}

