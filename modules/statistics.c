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
 * \file statistics.c
 * \brief Implementation for \ref Statistics module.
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/modules/statistics.c,v $
 * $Revision: 1.92 $ 
 * $Date: 2006/09/19 00:57:58 $
 */

// C library includes
#include <stdio.h> 
#include <stdlib.h>
#include <string.h>
#include <math.h>

// fc library dependencies
#include "base.h"
#include "storage.h"
#include "library.h"
#include "sequenceP.h"
#include "sequence.h"
#include "subset.h"
#include "variable.h"
#include "variableP.h"
#include "util.h"

// this module
#include "statistics.h"

/**
 * \addtogroup Statistics Statistics Routines
 * \brief  These routines compute simple statistical values of the 
 *         data variables.
 */

/** \name For a variable (spatial domain). */
//-------------------------------
//@{

/**
 * \ingroup  Statistics
 * \brief  Find min & max data values of a variable.
 *
 * \description
 *
 *    Returns the min and max values of a variable field.  Also returns the id
 *    of the entity which had the min and max values (the first one encountered
 *    if more than one node or element have the same value).  This routine only
 *    operates on fields with one component and error value is returned for
 *    vector and tensor fields.
 *
 * \modifications
 *   - Andre Smirnov  Created.
 *   - Nancy Collins  revamped.
 *   - Wendy Koegler  Added checking datatype. Also changed initial moving_max
 *   from MINFLOAT to -MAXFLOAT.
 *   - APR-08-2003  Wendy Koegler.  Added magnitude calculation for vectors.
 *   - 2003-JUL-08  W Koegler  Removed dependence on MINFLOAT & MAXFLOAT
 *   - 5/3/04 WSK, changed to only work on one component fields.
 */
FC_ReturnCode fc_getVariableMinMax(
  FC_Variable var,   /**< input - variable handle */
  double *min,       /**< output - min data value */
  int *min_index,    /**< output - index at which first min value was found */
  double *max,       /**< output - max data value */
  int *max_index     /**< output - index as which first max value was found */
) {
  FC_ReturnCode rc;
  int i;
  int numDataPoint, dim;
  void *datap;
  FC_DataType datatype;
  double value, moving_min, moving_max;
  int mini, maxi;
  
  // default output
  if (min)
    *min = -1;
  if (min_index)
    *min_index = -1;
  if (max)
    *max = -1;
  if (max_index)
    *max_index = -1;

  // check input - they don't have to take 
  if (!fc_isVariableValid(var) || (!min && !min_index && !max && !max_index)) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
      return FC_INPUT_ERROR;
    }

  // log message
  fc_printfLogMessage("Finding min and max values.");

  // get a count & pointer to the data array. 
  rc = fc_getVariableInfo(var, &numDataPoint, &dim, NULL, NULL, &datatype);
  if (rc != FC_SUCCESS)
    return rc;

  // can only do 1 component fields
  if (dim != 1 || numDataPoint < 1)  {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
      return FC_INPUT_ERROR;
    }
  
  // get data
  rc = fc_getVariableDataPtr(var, &datap);
  if (rc != FC_SUCCESS)
    return rc;
  
  // initialize with first point
  switch (datatype) {
  case FC_DT_INT :    value = ((int*)datap)[0];     break;
  case FC_DT_FLOAT:   value = ((float*)datap)[0];   break;
  case FC_DT_DOUBLE:  value = ((double*)datap)[0];  break;
  default: 
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
      return FC_INPUT_ERROR;
  }
  moving_min = moving_max = value;
  mini = maxi = 0;

  // look at the remaining points
  for (i = 1; i < numDataPoint; i++) {
    switch(datatype) {
    case FC_DT_INT :    value = ((int*)datap)[i];     break;
    case FC_DT_FLOAT:   value = ((float*)datap)[i];   break;
    case FC_DT_DOUBLE:  value = ((double*)datap)[i];  break;
    default: ;// nothing
    }
    if (value < moving_min) {
      moving_min = value;
      mini = i;
    }
    else if (value > moving_max) {
      moving_max = value;
      maxi = i;
    }
  }
  
  // only set the ones the user requested
  // ignore any which are null
  
  if (min != NULL)
    *min = moving_min;
  if (min_index != NULL)
    *min_index = mini;
  if (max != NULL)
    *max = moving_max;
  if (max_index != NULL)
    *max_index = maxi;
  
  return FC_SUCCESS;
}

/**
 * \ingroup  Statistics
 * \brief  Compute mean and standard deviation of a variable.
 *
 * \description
 *  
 *    Returns the mean and standard deviation of the values of a variable
 *    field. This routine only operates on fields with one component and
 *    and an error value is returned for vector and tensor fields.
 *
 *    In the case of variable with a single data point, the standard dev
 *    is technically undefined, but this routine will return 0.
 *
 * \modifications 
 *    - Andre Smirnov Created.
 *    - Nancy Collins  revamped.
 *    - Wendy Koegler
 *       Added checking datatype. Also changed initial moving_max
 *    - APR-08-2003  Wendy Koegler. Added magnitude calculation for vectors.
 */
FC_ReturnCode fc_getVariableMeanSdev(
  FC_Variable var,           /**< input - variable handle */
  double *mean,              /**< output - average */
  double *sdev               /**< output - standard deviation (optional) */
) {
  FC_ReturnCode rc;
  int i, numDataPoint, dim;
  FC_DataType datatype;
  void *datap;
  double temp, value;
  double moving_sum = 0.0;
  double moving_sumsq = 0.0;
  
  // default return
  if (mean)
    *mean = -1;
  if (sdev)
    *sdev = -1;

  // check input
  if (!fc_isVariableValid(var) || !mean) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
      return FC_INPUT_ERROR;
    }

  // log message
  fc_printfLogMessage("Computing mean and standard deviation.");

  // get a count & pointer to the data array
  rc = fc_getVariableInfo(var, &numDataPoint, &dim, NULL, NULL, &datatype);
  if (rc != FC_SUCCESS)
    return rc;

  // can only do 1 component fields
  if (dim != 1 || numDataPoint < 1)  {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
      return FC_INPUT_ERROR;
    }

  // get data
  rc = fc_getVariableDataPtr(var, &datap);
  if (rc != FC_SUCCESS)
    return rc;
  
  // Perform the calculation
  for (i = 0; i < numDataPoint; i++) {
    switch(datatype) {
    case FC_DT_INT:     value = ((int*)datap)[i];     break;
    case FC_DT_FLOAT:   value = ((float*)datap)[i];   break;
    case FC_DT_DOUBLE:  value = ((double*)datap)[i];  break;
    default: 
        fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
        return FC_INPUT_ERROR;
    }
    moving_sum += value;
    moving_sumsq += value*value;
  }
  
  // if user asked for mean, compute it here
  if (mean != NULL) 
    *mean = moving_sum / numDataPoint;
  
  // if user asked for std dev, compute it here
  // This is probably the "best" way to computer the
  // stdev as it minimizes the number of '-' operations that might
  // cause roundoff error.
  if (sdev != NULL) {
    if (numDataPoint > 1) { 
      temp = (moving_sumsq - (moving_sum * moving_sum / numDataPoint)) 
        / (numDataPoint - 1);
      *sdev = sqrt(fabs(temp)); // might be negative if roundoff error
    }
    else
      *sdev = 0;  //FIX? report NAN?
  }
  
  return FC_SUCCESS;
}

/**
 * \ingroup  Statistics
 * \brief  Compute sum of all data values of a variable.
 *
 * \description
 *  
 *    This routine only operates on fields with one component and
 *    and an error value is returned for vector and tensor fields. 
 *
 * \modifications
 *   - 02/19/04 WSK, Created
 */
FC_ReturnCode fc_getVariableSum(
  FC_Variable var,   /**< input - variable handle */
  double *sum_p      /**< output - sum of the data values */
) {
  FC_ReturnCode rc;
  int i, numDataPoint, dim;
  void *datap;
  FC_DataType datatype;
  double sum;

  // default output
  if (sum_p)
    *sum_p = -1;

  // check input
  if (!fc_isVariableValid(var) || !sum_p) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
      return FC_INPUT_ERROR;
  }

  // get a count & pointer to the data array. 
  rc = fc_getVariableInfo(var, &numDataPoint, &dim, NULL, NULL, &datatype);
  if (rc != FC_SUCCESS)
    return rc;

  // can only do 1 component fields
  if (dim != 1 || numDataPoint < 1)  {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // log message
  fc_printfLogMessage("Computing sum of all data values.");

  // get data
  rc = fc_getVariableDataPtr(var, &datap);
  if (rc != FC_SUCCESS)
    return rc;
  
  // Perform the calculation
  sum = 0.;
  for (i = 0; i < numDataPoint; i++) {
    switch(datatype) {
    case FC_DT_INT:     sum += ((int*)datap)[i];     break;
    case FC_DT_FLOAT:   sum += ((float*)datap)[i];   break;
    case FC_DT_DOUBLE:  sum += ((double*)datap)[i];  break;
    default: {
        fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
        return FC_INPUT_ERROR;
      }
    }
  }
  
  *sum_p = sum;
  return FC_SUCCESS;
}

//@}

/** \name For a subset of a variable (spatial domain). */
//-------------------------------
//@{

/**
 * \ingroup  Statistics
 * \brief  Find min & max data values in a subset of a variable.
 *
 * \description
 *
 *    Returns the min and max values of the subset of a variable field.  Also
 *    returns the id of the entity which had the min and max values (the first
 *    one encountered if more than one node or element have the same value).
 *    This routine only operates on fields with one component and error value
 *    is returned for vector and tensor fields.
 *
 *    The variable and the subset must be defined on the same type of
 *    entity (e.g. both on vertices or both on elements). They don't
 *    necessarily have to be on the same mesh as long as the total number
 *    of entities is the same. 
 *
 *    The function returns an error if the subset is empty.
 *
 * \modifications
 *   - Andre Smirnov  Created.
 *   - Nancy Collins  revamped.
 *   - Wendy Koegler  Added checking datatype. Also changed initial moving_max
 *   from MINFLOAT to -MAXFLOAT.
 *   - APR-08-2003  Wendy Koegler.  Added magnitude calculation for vectors.
 *   - 2003-JUL-08  W Koegler  Removed dependence on MINFLOAT & MAXFLOAT
 *   - 5/3/04 WSK, changed to only work on one component fields.
 */
FC_ReturnCode fc_getVariableSubsetMinMax(
  FC_Variable var,   /**< input - variable handle */
  FC_Subset subset,  /**< input - subset handle */
  double *min,       /**< output - min data value (optional) */
  int *min_index,    /**< output - index at which first min value
                        was found (optional) */
  double *max,       /**< output - max data value (optional) */
  int *max_index     /**< output - index as which first max value
                        was found  (optional) */
) {
  FC_ReturnCode rc;
  int i;
  int numDataPoint, dim, maxNumMember, numMember, *memberIDs;
  void *datap;
  FC_AssociationType assoc_var, assoc_subset;
  FC_DataType datatype;
  double value, moving_min, moving_max;
  int mini, maxi;
  
  // default output
  if (min)
    *min = -1;
  if (min_index)
    *min_index = -1;
  if (max)
    *max = -1;
  if (max_index)
    *max_index = -1;

  // check input - they don't have to take all of the returned values
  if (!fc_isVariableValid(var) || !fc_isSubsetValid(subset) ||
      (!min && !min_index && !max && !max_index) ) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
      return FC_INPUT_ERROR;
    }

  // log message
  fc_printfLogMessage("Computing min and max of subset of a variable");

  // get info from var & subset
  rc = fc_getVariableInfo(var, &numDataPoint, &dim, &assoc_var, NULL, 
                          &datatype);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_getSubsetInfo(subset, &numMember, &maxNumMember, &assoc_subset);
  if (rc != FC_SUCCESS)
    return rc;

  // A little more checking
  if (dim != 1 || numDataPoint < 1)  {
    fc_printfErrorMessage("Variable can only have 1 component not %d", dim);
    return FC_INPUT_ERROR;
  }
  if (numDataPoint < 1)  {
    fc_printfErrorMessage("Variable has no data");
    return FC_INPUT_ERROR;
  }
  if (assoc_var != assoc_subset || numDataPoint != maxNumMember) {
    fc_printfErrorMessage("Mismatch of entities");
    return FC_INPUT_ERROR;
  }
  if (numMember < 1) {
    fc_printfErrorMessage("Subset cannot be empty");
    return FC_INPUT_ERROR;
  }
  
  // get variable data && subset members
  rc = fc_getVariableDataPtr(var, &datap);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_getSubsetMembersAsArray(subset, &numMember, &memberIDs);
  if (rc != FC_SUCCESS)
    return rc;
  
  // initialize with first point
  switch (datatype) {
  case FC_DT_INT :    value = ((int*)datap)[memberIDs[0]];     break;
  case FC_DT_FLOAT:   value = ((float*)datap)[memberIDs[0]];   break;
  case FC_DT_DOUBLE:  value = ((double*)datap)[memberIDs[0]];  break;
  default: 
    free(memberIDs);
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  moving_min = moving_max = value;
  mini = maxi = memberIDs[0];

  // look at the remaining points
  for (i = 1; i < numMember; i++) {
    switch(datatype) {
    case FC_DT_INT :    value = ((int*)datap)[memberIDs[i]];     break;
    case FC_DT_FLOAT:   value = ((float*)datap)[memberIDs[i]];   break;
    case FC_DT_DOUBLE:  value = ((double*)datap)[memberIDs[i]];  break;
    default: ;// nothing
    }
    if (value < moving_min) {
      moving_min = value;
      mini = memberIDs[i];
    }
    else if (value > moving_max) {
      moving_max = value;
      maxi = memberIDs[i];
    }
  }
  free(memberIDs);

  // only set the ones the user requested
  // ignore any which are null
  
  if (min != NULL)
    *min = moving_min;
  if (min_index != NULL)
    *min_index = mini;
  if (max != NULL)
    *max = moving_max;
  if (max_index != NULL)
    *max_index = maxi;
  
  return FC_SUCCESS;
}

/**
 * \ingroup  Statistics
 * \brief  Compute mean and standard deviation of a subset of a variable.
 *
 * \description
 *  
 *    Returns the mean and standard deviation of the values of a subset of a
 *    variable field. This routine only operates on fields with one component 
 *    and and an error value is returned for vector and tensor fields.
 *
 *    In the case of subset with a single data point, the standard dev
 *    is technically undefined, but this routine will return 0.
 *
 *    The variable and the subset must be defined on the same type of
 *    entity (e.g. both on vertices or both on elements). They don't
 *    necessarily have to be on the same mesh as long as the total number
 *    of entities is the same. 
 *
 *    The function returns an error if the subset is empty.
 *
 * \modifications 
 *    - Andre Smirnov Created.
 *    - Nancy Collins  revamped.
 *    - Wendy Koegler
 *       Added checking datatype. Also changed initial moving_max
 *    - APR-08-2003  Wendy Koegler. Added magnitude calculation for vectors.
 */
FC_ReturnCode fc_getVariableSubsetMeanSdev(
  FC_Variable var,           /**< input - variable handle */
  FC_Subset subset,          /**< input - subset handle */
  double *mean,              /**< output - average */
  double *sdev               /**< output - standard deviation (optional) */
) {
  FC_ReturnCode rc;
  int i, numDataPoint, dim, maxNumMember, numMember, *memberIDs;
  FC_AssociationType assoc_var, assoc_subset;
  FC_DataType datatype;
  void *datap;
  double temp, value;
  double moving_sum = 0.0;
  double moving_sumsq = 0.0;
  
  // default return
  if (mean)
    *mean = -1;
  if (sdev)
    *sdev = -1;

  // check input
  if (!fc_isVariableValid(var) || !fc_isSubsetValid(subset) || !mean) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // log message
  fc_printfLogMessage("Computing mean and standard deviation of a subset "
                      " of a variable");

  // get info from var & subset
  rc = fc_getVariableInfo(var, &numDataPoint, &dim, &assoc_var, NULL, 
                          &datatype);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_getSubsetInfo(subset, &numMember, &maxNumMember, &assoc_subset);
  if (rc != FC_SUCCESS)
    return rc;

  // A little more checking
  if (dim != 1 || numDataPoint < 1)  {
    fc_printfErrorMessage("Variable can only have 1 component not %d", dim);
    return FC_INPUT_ERROR;
  }
  if (numDataPoint < 1)  {
    fc_printfErrorMessage("Variable has no data");
    return FC_INPUT_ERROR;
  }
  if (assoc_var != assoc_subset || numDataPoint != maxNumMember) {
    fc_printfErrorMessage("Mismatch of entities");
    return FC_INPUT_ERROR;
  }
  if (numMember < 1) {
    fc_printfErrorMessage("Subset cannot be empty");
    return FC_INPUT_ERROR;
  }
  
  // get variable data && subset members
  rc = fc_getVariableDataPtr(var, &datap);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_getSubsetMembersAsArray(subset, &numMember, &memberIDs);
  if (rc != FC_SUCCESS)
    return rc;
    
  // Perform the calculation
  for (i = 0; i < numMember; i++) {
    switch(datatype) {
    case FC_DT_INT:     value = ((int*)datap)[memberIDs[i]];     break;
    case FC_DT_FLOAT:   value = ((float*)datap)[memberIDs[i]];   break;
    case FC_DT_DOUBLE:  value = ((double*)datap)[memberIDs[i]];  break;
    default: 
      free(memberIDs);
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
      return FC_INPUT_ERROR;
    }
    moving_sum += value;
    moving_sumsq += value*value;
  }
  free(memberIDs);
  
  // if user asked for mean, compute it here
  if (mean != NULL) 
    *mean = moving_sum / numMember;
  
  // if user asked for std dev, compute it here
  // This is probably the "best" way to computer the
  // stdev as it minimizes the number of '-' operations that might
  // cause roundoff error.
  if (sdev != NULL) {
    if (numMember > 1) { 
      temp = (moving_sumsq - (moving_sum * moving_sum / numMember)) 
        / (numMember - 1);
      *sdev = sqrt(fabs(temp)); // might be negative if roundoff error
    }
    else
      *sdev = 0;  //FIX? report NAN?
  }
  
  return FC_SUCCESS;
}

/**
 * \ingroup  Statistics
 * \brief  Compute sum of all data values of a subset of a variable.
 *
 * \description
 *  
 *    This routine only operates on fields with one component and
 *    and an error value is returned for vector and tensor fields. 
 *
 *    The variable and the subset must be defined on the same type of
 *    entity (e.g. both on vertices or both on elements). They don't
 *    necessarily have to be on the same mesh as long as the total number
 *    of entities is the same. 
 *
 *    The function returns an error if the subset is empty.
 *
 * \modifications
 *   - 02/19/04 WSK, Created
 */
FC_ReturnCode fc_getVariableSubsetSum(
  FC_Variable var,   /**< input - variable handle */
  FC_Subset subset,  /**< input - subset handle */
  double *sum_p      /**< output - sum of the data values */
) {
  FC_ReturnCode rc;
  int i, numDataPoint, dim, maxNumMember, numMember, *memberIDs;
  void *datap;
  FC_AssociationType assoc_var, assoc_subset;
  FC_DataType datatype;
  double sum;
  
  // default output
  if (sum_p)
    *sum_p = -1;

  // check input
  if (!fc_isVariableValid(var) || !fc_isSubsetValid(subset) || !sum_p) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // get info from var & subset
  rc = fc_getVariableInfo(var, &numDataPoint, &dim, &assoc_var, NULL, 
                          &datatype);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_getSubsetInfo(subset, &numMember, &maxNumMember, &assoc_subset);
  if (rc != FC_SUCCESS)
    return rc;

  // A little more checking
  if (dim != 1 || numDataPoint < 1)  {
    fc_printfErrorMessage("Variable can only have 1 component not %d", dim);
    return FC_INPUT_ERROR;
  }
  if (numDataPoint < 1)  {
    fc_printfErrorMessage("Variable has no data");
    return FC_INPUT_ERROR;
  }
  if (assoc_var != assoc_subset || numDataPoint != maxNumMember) {
    fc_printfErrorMessage("Mismatch of entities");
    return FC_INPUT_ERROR;
  }
  if (numMember < 1) {
    fc_printfErrorMessage("Subset cannot be empty");
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Computing sum of all data values of subset of variable.");

  // get variable data && subset members
  rc = fc_getVariableDataPtr(var, &datap);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_getSubsetMembersAsArray(subset, &numMember, &memberIDs);
  if (rc != FC_SUCCESS)
    return rc;
  
  // Perform the calculation
  sum = 0.;
  for (i = 0; i < numMember; i++) {
    switch(datatype) {
    case FC_DT_INT:     sum += ((int*)datap)[memberIDs[i]];     break;
    case FC_DT_FLOAT:   sum += ((float*)datap)[memberIDs[i]];   break;
    case FC_DT_DOUBLE:  sum += ((double*)datap)[memberIDs[i]];  break;
    default: 
      free(memberIDs);
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
      return FC_INPUT_ERROR;
    }
  }
  free(memberIDs);
  
  *sum_p = sum;
  return FC_SUCCESS;
}

//@}

/** \name For a sequence variable (spatial & sequence domain). */
//-------------------------------
//@{

/**
 * \ingroup  Statistics
 * \brief  Find min & max data values of a sequence variable.
 *
 * \description
 *
 *    Returns the min and max values of a sequence variable. Also returns
 *    sequence step id and the id of the entity which had the min and max
 *    values (the first one encountered if more than one node or element have
 *    the same value).  This routine only operates on fields with one component
 *    and error value is returned for vector and tensor fields.
 *
 * \modifications
 *    - 2006-07-24 WSD Created.
 */
FC_ReturnCode fc_getSeqVariableMinMax(
  int numStep,         /**< input - the number of steps in the seq variable */
  FC_Variable* seqVar, /**< input - the seq variable's handles */
  double *min,         /**< output - the min data value */
  int *min_seq_id,     /**< output - the step ID at which the first min value
			  was found*/
  int *min_entity_id,  /**< output - the entity ID at which first min value was 
			 found */
  double *max,         /**< output - max data value */
  int *max_seq_id,     /**< output - the step ID at which the first max value
			  was found */
  int *max_entity_id   /**< output - the entity ID at which first max value was 
			  found */
) {
  FC_ReturnCode rc;
  int i, j;
  int numDataPoint, dim;
  void *datap;
  FC_DataType datatype;
  double value, moving_min, moving_max;
  int minsi, minei, maxsi, maxei;
  
  // default output
  if (min)
    *min = -1;
  if (min_seq_id)
    *min_seq_id = -1;
  if (min_entity_id)
    *min_entity_id = -1;
  if (max)
    *max = -1;
  if (max_seq_id)
    *max_seq_id = -1;
  if (max_entity_id)
    *max_entity_id = -1;

  // check input - they don't have to take 
  if (!fc_isSeqVariableValid(numStep, seqVar) || 
      (!min && !min_seq_id && !min_entity_id && !max && !max_seq_id && !max_entity_id)) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // log message
  fc_printfLogMessage("Finding min and max values.");

  // get a count & pointer to the data array. 
  rc = fc_getVariableInfo(seqVar[0], &numDataPoint, &dim, NULL, NULL, &datatype);
  if (rc != FC_SUCCESS)
    return rc;

  // can only do 1 component fields
  if (dim != 1 || numDataPoint < 1)  {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // initialize the first point
  rc = fc_getVariableDataPtr(seqVar[0], &datap);
  if (rc != FC_SUCCESS)
    return rc;
  switch (datatype) {
  case FC_DT_INT :    value = ((int*)datap)[0];     break;
  case FC_DT_FLOAT:   value = ((float*)datap)[0];   break;
  case FC_DT_DOUBLE:  value = ((double*)datap)[0];  break;
  default: 
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  moving_min = moving_max = value;
  minsi = minei = maxsi = maxei = 0;

  // Look at all points (will repeat first one because less coding work)		
  for (i = 0; i < numStep; i++) {
    rc = fc_getVariableDataPtr(seqVar[i], &datap);
    if (rc != FC_SUCCESS)
      return rc;
    for (j = 0; j < numDataPoint; j++) {
      switch(datatype) {
      case FC_DT_INT :    value = ((int*)datap)[j];     break;
      case FC_DT_FLOAT:   value = ((float*)datap)[j];   break;
      case FC_DT_DOUBLE:  value = ((double*)datap)[j];  break;
      default: ;// nothing
      }
      if (value < moving_min) {
	moving_min = value;
	minsi = i;
	minei = j;
      }
      else if (value > moving_max) {
	moving_max = value;
	maxsi = i;
	maxei = j;
      }
    }
  }
  
  // only set the ones the user requested  
  if (min != NULL)
    *min = moving_min;
  if (min_seq_id != NULL)
    *min_seq_id = minsi;
  if (min_entity_id != NULL)
    *min_entity_id = minei;
  if (max != NULL)
    *max = moving_max;
  if (max_seq_id != NULL)
    *max_seq_id = maxsi;
  if (max_entity_id != NULL)
    *max_entity_id = maxei;
  
  return FC_SUCCESS;
}

/**
 * \ingroup  Statistics
 * \brief  Compute mean and standard deviation of a sequence variable.
 *
 * \description
 *
 *    Returns the mean and standard deviation of the values of all steps of a
 *    sequence variable. This routine only operates on fields with one
 *    component and error value is returned for vector and tensor fields.
 *
 * \modifications
 *    - 2006-07-24 WSD Created.
 */
FC_ReturnCode fc_getSeqVariableMeanSdev(
  int numStep,         /**< input - the number of steps in the seq variable */
  FC_Variable* seqVar, /**< input - the seq variable's handles */
  double *mean,        /**< output - average */
  double *sdev         /**< output - standard deviation(optional) */
) {
  FC_ReturnCode rc;
  int i, j;
  int numDataPoint, numTotal, dim;
  void *datap;
  FC_DataType datatype;
  double temp, value;
  double moving_sum = 0.0;
  double moving_sumsq = 0.0;

  // default output
  if (mean)
    *mean = -1;
  if (sdev)
    *sdev = -1;

  // check input - they don't have to take 
  if (!fc_isSeqVariableValid(numStep, seqVar) || !mean) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  fc_printfLogMessage("Computing mean and standard deviation.");

  // get a count & pointer to the data array. 
  rc = fc_getVariableInfo(seqVar[0], &numDataPoint, &dim, NULL, NULL, &datatype);
  if (rc != FC_SUCCESS)
    return rc;

  // can only do 1 component fields
  if (dim != 1 || numDataPoint < 1)  {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  

  // collect sums
  numTotal = numStep*numDataPoint;
  for (i = 0; i < numStep; i++) {
    // get data
    rc = fc_getVariableDataPtr(seqVar[i], &datap);
    if (rc != FC_SUCCESS)
      return rc;
  
    // Perform the calculation
    for (j = 0; j < numDataPoint; j++) {
      switch(datatype) {
      case FC_DT_INT:     value = ((int*)datap)[j];     break;
      case FC_DT_FLOAT:   value = ((float*)datap)[j];   break;
      case FC_DT_DOUBLE:  value = ((double*)datap)[j];  break;
      default: 
        fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
        return FC_INPUT_ERROR;
      }
      moving_sum += value;
      moving_sumsq += value*value;
    }
  }
  
  // if user asked for mean, compute it here
  if (mean != NULL) 
    *mean = moving_sum / numTotal;
  
  // if user asked for std dev, compute it here
  // This is probably the "best" way to computer the
  // stdev as it minimizes the number of '-' operations that might
  // cause roundoff error.
  if (sdev != NULL) {
    if (numTotal > 1) { 
      temp = (moving_sumsq - (moving_sum * moving_sum / numTotal)) 
        / (numTotal - 1);
      *sdev = sqrt(fabs(temp)); // might be negative if roundoff error
    }
    else
      *sdev = 0;  //FIX? report NAN?
  }
  
  return FC_SUCCESS;
}

//@}

/** \name For each data point of a sequence variable (sequence domain). */
//-------------------------------
//@{

/**
 * \ingroup  Statistics
 * \brief  Find min & max data values of a sequence variable across
 *         its sequence (series) on a per data point, per component basis.
 *
 * \description
 *
 *    Returns the min and max values of a variable field (and
 *    sequence (e.g., time) indicies of them - see below). 
 *    This routine treats each data point and each component of
 *    the sequence variable separately. Thus, the min and max are calc
 *    across a series (e.g., time) rather than over the data
 *    points as it is in the fc_getVariableXXX functions. 
 *    This routine also returns the sequence index of the min and max
 *    (the first one encountered if the var has teh same value at more
 *    than one point in "time").
 *
 *    NOTE:
 *    - Outputs are  FC_Variables with each value at the corresponding
 *    data array location, rather than raw data - unlike the 
 *    fc_getVariableMinMax function. The output variable names are
 *    the name of the original subset, postpended with "_SeriesMin",
 *    "_SeriesMax", "_SeriesMinIndex" and "_SeriesMaxIndex".   
 *    - Return type of variable data are FC_DT_DOUBLE as in fc_getVariableMinMax,
 *    rather than native type.
 *    - Unlike fc_getVariableMinMax, none of the output variables
 *        is optional
 *
 *  \todo  see about making the output options optional.
 *
 *
 * \modifications
 *   - 1/18/05 ACG Created.
 *   - 2/3/05 ACG changed name to have "Time" in title. While the sequence
 *            is not necessarily a time sequence, this makes it clear that
 *            the function is over the sequence and not over the data
 *            points, as in the fc_getVariableXXX functions in this module.
 *   - 2/4/05 ACG changing interface to return FC_Variables at request of WSK. 
 *   - 2/9/05 ACG removed output name options as the signature was too
 *     unwieldy.
 *   - 6/1/05 ACG renamed to "series" rather than "time" to reflect change
 *     of timeseries module to series module
 *   - 1/09/06 ACG removing cleanup in case of developer error at suggest of WSK
 *     becuase its confusing hte results of the coverage tests. if things reach
 *     that point, it should really be killed so cleanup doesnt matter.
 */
FC_ReturnCode fc_getSeqVariableSeriesMinMax(
  int numStep,     /**<input - numStep in sequenceVariable */
  FC_Variable *seqvar, /**< input - seq variable handle */
  FC_Variable *minvar, /**< output - min var variable */
  FC_Variable *maxvar, /**< output - max var variable */
  FC_Variable *minindexvar, /**< output - min series index var variable */
  FC_Variable *maxindexvar /**< output - max series index var variable */
 ){
  FC_ReturnCode rc;
  FC_Mesh mesh;
  _FC_VarSlot* varSlot;
  int numDataPoint;
  int numComponent;
  FC_AssociationType assoc;
  FC_MathType mathtype;
  FC_DataType datatype;
  int var_dim;
  int i,j;
  double *mindata, *maxdata;
  int *minindexdata, *maxindexdata;
  void *data;

  char* varname;
  char* name;
  void* mdata;
  FC_Variable* mvar;
  FC_DataType mdt;


  typedef enum {
    MMIN = 0,
    MMAX = 1,
    MMINI = 2,
    MMAXI = 3
  } junktype;


  if (minvar)
    *minvar = FC_NULL_VARIABLE;
  if (maxvar)
    *maxvar = FC_NULL_VARIABLE;
  if (minindexvar)
    *minindexvar = FC_NULL_VARIABLE;
  if (maxindexvar)
    *maxindexvar = FC_NULL_VARIABLE;



  //check input
  if (minvar == NULL || 
      maxvar == NULL || 
      minindexvar == NULL || 
      maxindexvar == NULL || 
      (!fc_isSeqVariableValid(numStep,seqvar))){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }


  varSlot = _fc_getVarSlot(seqvar[0]);
  fc_printfLogMessage("Getting min max for seq variable '%s'",
                      varSlot->header.name);

  rc = fc_getVariableInfo(seqvar[0], &numDataPoint, &numComponent, &assoc,
                           &mathtype, &datatype);
  if (rc != FC_SUCCESS){
    return rc;
  }

  if (!(datatype == FC_DT_FLOAT ||
        datatype == FC_DT_DOUBLE ||
        datatype == FC_DT_INT )
      ){
    return FC_INPUT_ERROR;
  } 

  //finally ok

  var_dim = numDataPoint*numComponent;
  

  mindata = (double*)malloc(var_dim*sizeof(double));
  maxdata = (double*)malloc(var_dim*sizeof(double));
  minindexdata = (int*)malloc(var_dim*sizeof(int));
  maxindexdata = (int*)malloc(var_dim*sizeof(int));
  if (mindata == NULL || maxdata == NULL || minindexdata == NULL ||
      maxindexdata == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    free(mindata);
    free(maxdata);
    free(minindexdata);
    return FC_MEMORY_ERROR;
  }

  
  //initialize with first point
  rc = fc_getVariableDataPtr(seqvar[0], &data);
  if (rc != FC_SUCCESS){
    data = NULL;
    free(mindata);
    free(maxdata);
    free(minindexdata);
    free(maxindexdata);
    return rc;
  }
  for (j = 0; j < var_dim; j++){
    double quan = 0.0;
    switch(datatype){
    case FC_DT_DOUBLE:
      quan= ((double*)data)[j];
      break;
    case FC_DT_INT:
      quan= (double)(((int*)data)[j]);
      break;
    case FC_DT_FLOAT:
      quan= (double)(((float*)data)[j]);
      break;
    default:
      printf("** Developer Alert! Should never reach this point!**\n");
      fflush(NULL);
      fc_printfErrorMessage("%s", "Invalid id type");
      return FC_INPUT_ERROR;
      break;
    }
    mindata[j] = quan; 
    maxdata[j] = quan; 
    minindexdata[j]= 0;
    maxindexdata[j]= 0;
    }

  
  for (i = 1; i < numStep; i++){
    rc = fc_getVariableDataPtr(seqvar[i], &data);
    if (rc != FC_SUCCESS){
      data = NULL;
      free(mindata);
      free(maxdata);
      free(minindexdata);
      free(maxindexdata);
      return rc;
    }
    for (j = 0; j < var_dim; j++){
      double quan = 0.0;
      switch(datatype){
      case FC_DT_DOUBLE:
        quan = ((double*)data)[j];
        break;
      case FC_DT_INT:
        quan = (double)(((int*)data)[j]);
        break;
      case FC_DT_FLOAT:
        quan = (double)(((float*)data)[j]);
        break;
      default:
        printf("** Developer Alert! Should never reach this point!**\n");
        fflush(NULL);
        fc_printfErrorMessage("%s", "Invalid id type");
        return FC_INPUT_ERROR;
        break;
      }
      if (quan < mindata[j]){
        mindata[j] = quan;
        minindexdata[j] = i;
      }
      if (quan > maxdata[j]){
        maxdata[j] = quan;
        maxindexdata[j] = i;
      }
    }
  }

  rc = fc_getMeshFromVariable(seqvar[0],&mesh);
  if (rc != FC_SUCCESS){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(rc));
    free(mindata);
    free(maxdata);
    free(minindexdata);
    free(maxindexdata);
    return rc;
  }


  varname = varSlot->header.name;
  name = malloc((strlen(varname)+20)*sizeof(char));
  if (name == NULL){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  for (i = 0; i < 4; i++){
    switch (i){
    case MMIN:
      sprintf(name, "%s_SeriesMin", varname);
      mdata = mindata;
      mvar = minvar;
      mdt = FC_DT_DOUBLE;
      break;
    case MMAX:
      sprintf(name, "%s_SeriesMax", varname);
      mdata = maxdata;
      mvar = maxvar;
      mdt = FC_DT_DOUBLE;
      break;
    case MMINI:
      sprintf(name, "%s_SeriesMinIndex", varname);
      mdata = minindexdata;
      mvar = minindexvar;
      mdt = FC_DT_INT;
      break;
    case MMAXI:
      sprintf(name, "%s_SeriesMaxIndex", varname);
      mdata = maxindexdata;
      mvar = maxindexvar;
      mdt = FC_DT_INT;
      break;
    default:
      printf("** Developer Alert! Should never reach this point!**\n");
      fc_printfErrorMessage("%s", "Invalid id type");
      return FC_ERROR;
      break;
    }
    rc = fc_createVariable(mesh,name,mvar);
    if (rc != FC_SUCCESS){
      free(name);
      fc_deleteVariable(*minvar);
      fc_deleteVariable(*maxvar);
      fc_deleteVariable(*minindexvar);
      fc_deleteVariable(*maxindexvar);
      if (mindata) free(mindata);
      if (maxdata) free(maxdata);
      if (minindexdata) free(minindexdata);
      if (maxindexdata) free(maxindexdata);
      return rc;
    }

    rc = fc_setVariableDataPtr(*mvar,numDataPoint, numComponent,
                             assoc,mathtype,mdt,mdata);
    if (rc != FC_SUCCESS){
      free(name);
      fc_deleteVariable(*minvar);
      fc_deleteVariable(*maxvar);
      fc_deleteVariable(*minindexvar);
      fc_deleteVariable(*maxindexvar);
      if (mindata) free(mindata);
      if (maxdata) free(maxdata);
      if (minindexdata) free(minindexdata);
      if (maxindexdata) free(maxindexdata);
      return rc;
    }
  }

  free(name);

  return FC_SUCCESS;
}

/**
 * \ingroup  Statistics
 * \brief  Compute mean and standard deviation of a seq variable across
 *         its sequence (series) on a per data point, per component basis.
 *
 * \description
 *  
 *    Returns the mean and standard deviation of the values of a variable
 *    field. This routine treats each data point and each component of the
 *    sequence variable separately. Thus, the computation is over the
 *    sequence quantity (e.g., time) rather than over the data
 *    points as it is in the fc_getVariableXXX functions. 
 *
 *    NOTE:
 *    - Output is an FC_Variable with each sum at the corresponding
 *    data array location, rather than raw data - unlike the 
 *    fc_getVariableMeanSdevfunction. The output variable names are
 *    the name of the original subset, postpended with "_SeriesMean" or
 *    "_SeriesSdev". 
 *    - Return type of data array is FC_DT_DOUBLE as in fc_getVariableMeanSdev.
 *    - Like other methods in this module, this is the square root
 *    of the bias corrected variance, i.e.,
 *    sqrt(((1/(N-1)) * sum((xi - xave)**2))
 *    so sample variance =sqrt((N-1)/N) * stddev_from_here where N = 
 *    number of steps in the series.
 *    - Stdev is an optional arg.
 *    - In the case of variable with a single data point, the standard dev
 *    is technically undefined, but this routine will return 0.
 *
 *
 * \modifications 
 *    - 1/19/05 ACG Created.
 *    - 2/3/05 ACG changed name to have "Time" in title. While the sequence
 *            is not necessarily a time sequence, this makes it clear that
 *            the function is over the sequence and not over the data
 *            points, as in the fc_getVariableXXX functions in this module.
 *    - 2/3/05 ACG changed output type to FC_Variable at WSK request.
 *    - 2/9/05 ACG removed output name option to be consistent with
 *            other sequence variable functions in the module
 *    - 6/2/05 ACG changing name from "TimeMean" to "SeriesMean"
 */
FC_ReturnCode fc_getSeqVariableSeriesMeanSdev(
  int numStep,     /**< input - numStep in sequenceVariable */
  FC_Variable *seqvar, /**< input - seq variable handle */
  FC_Variable *meanvar,    /**< output - mean data value array*/
  FC_Variable *sdevvar    /**< output - sdev data value array*/
  ){
  FC_ReturnCode rc;
  _FC_VarSlot* varSlot;
  FC_Mesh mesh;
  int numDataPoint;
  int numComponent;
  FC_AssociationType assoc;
  FC_MathType mathtype;
  FC_DataType datatype;
  int sum_dim;
  int i,j;
  void *data;
  double *sum, *sumsq, *temp, *mean, *sdev;

  char *meanname, *sdevname, *varname;


  if (meanvar)
    *meanvar = FC_NULL_VARIABLE;


  if (sdevvar)
    *sdevvar = FC_NULL_VARIABLE;

  //check input
  if (meanvar == NULL ||
      !fc_isSeqVariableValid(numStep,seqvar)){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }


  varSlot = _fc_getVarSlot(seqvar[0]);
  fc_printfLogMessage("Getting mean and sdev for seq variable '%s'",
                      varSlot->header.name);


  rc = fc_getVariableInfo(seqvar[0], &numDataPoint, &numComponent, &assoc,
                           &mathtype, &datatype);
  if (rc != FC_SUCCESS){
    return rc;
  }

  if (!(datatype == FC_DT_FLOAT ||
        datatype == FC_DT_DOUBLE ||
        datatype == FC_DT_INT )
      ){
    return FC_INPUT_ERROR;
  }

  sum_dim = numDataPoint*numComponent;
  sum = (double*)malloc(sum_dim*sizeof(double));
  if (sum == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  sumsq = (double*)malloc(sum_dim*sizeof(double));
  if (sumsq == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }

  for (i =0; i < sum_dim; i++){
    sum[i] = 0.0;
    sumsq[i] = 0.0;
  }

  for (i = 0; i < numStep; i++){
    rc = fc_getVariableDataPtr(seqvar[i], &data);
    if (rc != FC_SUCCESS){
      data = NULL;
      free(sum);
      free(sumsq);
      return rc;
    }
    for (j = 0; j < sum_dim; j++){
      double quan = 0.0;
      switch(datatype){
      case FC_DT_DOUBLE:
        quan = ((double*)data)[j];
        break;
      case FC_DT_INT:
        quan = (double)(((int*)data)[j]);
        break;
      case FC_DT_FLOAT:
        quan = (double)(((float*)data)[j]);
        break;
      default:
        printf("** Developer Alert! Should never reach this point!**\n");
        fflush(NULL);
        fc_printfErrorMessage("%s", "Invalid data type");
        return FC_ERROR;
        break;
      }
      sum[j]+= quan;
      sumsq[j]+= quan*quan;
    }
  }


  temp = (double*)malloc(sum_dim*sizeof(double));
  if (temp == NULL) {
    free(sum);
    free(sumsq);
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }

  mean = (double*)malloc(sum_dim*sizeof(double));
  if (mean == NULL) {
    free(temp);
    free(sum);
    free(sumsq);
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }

  if (sdevvar != NULL){
    sdev = (double*)malloc(sum_dim*sizeof(double));
    if (sdev == NULL) {
      free(temp);
      free(mean);
      free(sum);
      free(sumsq);
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      return FC_MEMORY_ERROR;
    }
  }

  for (j = 0; j < sum_dim; j++){
    mean[j] = sum[j]/numStep;
    
    if (sdevvar != NULL){
      if (numStep > 1) { 
        temp[j] = (sumsq[j] - (sum[j] * sum[j] / numStep)) / (numStep - 1);
        sdev[j] = sqrt(fabs(temp[j]));
        // might be negative if roundoff error
      }
      else{
        sdev[j] = 0.0;  //FIX? report NAN?
      }
    }
  }

  free(temp);
  free(sum);
  free(sumsq);

  rc = fc_getMeshFromVariable(seqvar[0],&mesh);
  if (rc != FC_SUCCESS){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(rc));
    free(mean);
    if (sdev) free(sdev);
    return rc;
  }


  varname = varSlot->header.name;
  meanname = malloc((strlen(varname)+20)*sizeof(char));
  if (meanname == NULL) {
    free(mean);
    if (sdev) free(sdev);
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  sprintf(meanname, "%s_SeriesMean", varname);



  rc = fc_createVariable(mesh,meanname,meanvar);
  if (rc != FC_SUCCESS){
    free(mean);
    if (meanname) free(meanname);
    if (sdev) free(sdev);
    return rc;
  }

  free(meanname);

  rc = fc_setVariableDataPtr(*meanvar,numDataPoint, numComponent,
                             assoc,mathtype,FC_DT_DOUBLE,mean);
  if (rc != FC_SUCCESS){
    fc_deleteVariable(*meanvar);
    if (mean) free(mean);
    if (sdev) free(sdev);
    return rc;
  }


  if (sdevvar != NULL){
    sdevname = malloc((strlen(varname)+20)*sizeof(char));
    if (sdevname == NULL) {
      free(sdev);
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      return FC_MEMORY_ERROR;
    }
    sprintf(sdevname, "%s_SeriesSdev", varname);
    
    rc = fc_createVariable(mesh,sdevname,sdevvar);
    if (rc != FC_SUCCESS){
      fc_deleteVariable(*meanvar);
      free(sdevname);
      free(sdev);
      return rc;
    }

    free(sdevname);

    rc = fc_setVariableDataPtr(*sdevvar,numDataPoint, numComponent,
                               assoc,mathtype,FC_DT_DOUBLE,sdev);
    if (rc != FC_SUCCESS){
      fc_deleteVariable(*meanvar);
      fc_deleteVariable(*sdevvar);
      if (sdev) free(sdev);
      return rc;
    }
  }

  return FC_SUCCESS;
}

/**
 * \ingroup  Statistics
 * \brief  Compute sum of all data values of a sequence variable 
 *         across its sequence (series) on a per data point,
 *         per component basis.
 *
 * \description
 *
 *    Returns the sum of a sequence variable field. This
 *    routine treats each data point and each component of
 *    the sequence variable separately. Thus, the sum is over the
 *    sequence quantity (e.g., time) rather than over the data
 *    points as it is in the fc_getVariableXXX functions. 
 *
 *    NOTE:
 *    - Output is an FC_Variable with each sum at the corresponding
 *    data array location, rather than raw data - unlike the 
 *    fc_getVariableSum function. The output variable name is
 *    the name of the original subset, postpended with "_SeriesSum" 
 *    - Return type of variable data is FC_DT_DOUBLE as in fc_getVariableSum,
 *    rather than native type.
 *
 *
 * \modifications
 *   - 1/18/05 ACG Created.
 *   - 2/3/05 ACG changed name to have "Time" in title. While the sequence
 *            is not necessarily a time sequence, this makes it clear that
 *            the function is over the sequence and not over the data
 *            points, as in the fc_getVariableXXX functions in this module.
 *   - 2/3/05 ACG changed output type to FC_Variable at WSK request.
 *   - 2/9/05 ACG removed output name option to be consistent with
 *            other sequence variable functions in the module
 *   - 6/2/05 ACG changed title from "TimeSum" to "SeriesSum"
 */
FC_ReturnCode fc_getSeqVariableSeriesSum(
  int numStep,     /**< input - numStep in sequenceVariable */
  FC_Variable *seqvar, /**< input - seq variable handle */
  FC_Variable *sumvar    /**< output - sum data value car*/
 ){
  FC_ReturnCode rc;
  _FC_VarSlot* varSlot;
  FC_Mesh mesh;
  int numDataPoint;
  int numComponent;
  FC_AssociationType assoc;
  FC_MathType mathtype;
  FC_DataType datatype;
  int sum_dim;
  int i,j;
  void *data;
  double* sum;
  char *varname, *returnname;


  if (sumvar)
    *sumvar = FC_NULL_VARIABLE;

  //check input
  if (sumvar == NULL ||
      (!fc_isSeqVariableValid(numStep,seqvar))){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }


  varSlot = _fc_getVarSlot(seqvar[0]);
  fc_printfLogMessage("Getting sum for seq variable '%s'",
                      varSlot->header.name);


  rc = fc_getVariableInfo(seqvar[0], &numDataPoint, &numComponent, &assoc,
                           &mathtype, &datatype);
  if (rc != FC_SUCCESS){
    return rc;
  }

  if (!(datatype == FC_DT_FLOAT ||
        datatype == FC_DT_DOUBLE ||
        datatype == FC_DT_INT )
      ){
    return FC_INPUT_ERROR;
  } 

  sum_dim = numDataPoint*numComponent;
  sum = (double*)malloc(sum_dim*sizeof(double));
  if (sum == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }

  for (i =0; i < sum_dim; i++){
    sum[i] = 0;
  }

  for (i = 0; i < numStep; i++){
    rc = fc_getVariableDataPtr(seqvar[i], &data);
    if (rc != FC_SUCCESS){
      data = NULL;
      free(sum);
      return rc;
    }
    for (j = 0; j < sum_dim; j++){
      switch(datatype){
      case FC_DT_DOUBLE:
        sum[j]+= ((double*)data)[j];
        break;
      case FC_DT_INT:
        sum[j]+= (double)(((int*)data)[j]);
        break;
      case FC_DT_FLOAT:
        sum[j]+= (double)(((float*)data)[j]);
        break;
      default:
        printf("** Developer Alert! Should never reach this point!**\n");
        fflush(NULL);
        fc_printfErrorMessage("%s", "Invalid data type");
        return FC_ERROR;
        break;
      }
    }
  }

  rc = fc_getMeshFromVariable(seqvar[0],&mesh);
  if (rc != FC_SUCCESS){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(rc));
    free(sum);
    return rc;
  }


  varname = varSlot->header.name;
  returnname = malloc((strlen(varname)+20)*sizeof(char));
  if (returnname == NULL) {
    free(sum);
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  sprintf(returnname, "%s_SeriesSum", varname);

  rc = fc_createVariable(mesh,returnname,sumvar);
  if (rc != FC_SUCCESS){
    free(sum);
    free(returnname);
    return rc;
  }

  free(returnname);
  rc = fc_setVariableDataPtr(*sumvar,numDataPoint, numComponent,
                             assoc,mathtype,FC_DT_DOUBLE,sum);
  if (rc != FC_SUCCESS){
    fc_deleteVariable(*sumvar);
    return rc;
  }

  return FC_SUCCESS;
}

//@}

/** \name For a sequence (sequence domain). */
//-------------------------------
//@{

/**
 * \ingroup Statistics
 * \brief Returns min and max of a sequence and the indicies where they occur 
 *        and returns monotonicity indicator.
 *
 * \description 
 *
 *        Returns min and max of a sequence and indicies of where
 *        they occur. in case of repeated values returns the index
 *        of the first occurrence. 
 *
 *        Return vals for mono are: 
 *        - montonically increasing = 1,
 *        - monotonically decreasing = -1,
 *        - single point sequence can be either so = 2,
 *        - and 0 otherwise.
 *
 *        Nonnumerical sequences return with an error. In case of error,
 *        return vals are -1, except for mono which returns 0.
 *
 *        Any return arguments can be passed in NULL. Return vals
 *        of min and max are doubles, despite the type of the
 *        original sequence. 
 * 
 *        If you are doing any comparisons afterward you should do them
 *        relative to the native tye of the sequence (e.g.< float) even
 *        though the return type is double becuase of casting errors.
 *
 * \todo
 *    - is there any reason we would want to do this in its native type
 *      throughout ? 
 *
 *
 * \modifications
 *    - 5/12/05 ACG Created
 *    - 5/16/05 ACG changed so single point seq returns 2 for monotonicity
 *      (this gives way for it to match either increasing or decreasing seqs)
 *    - 7/01/05 ACG added regularity check
 *    - 7/07/05 ACG removed regularity check
 */
FC_ReturnCode fc_getSequenceMinMaxMono(
  FC_Sequence sequence, /**< input - sequence */
  double *min,       /**< output - min data value */
  int *min_index,    /**< output - index at which first min value was found */
  double *max,       /**< output - max data value */
  int *max_index,     /**< output - index as which first max value was found */
  int* mono          /**< output - integer indicating monotonicity */
  ){
  FC_ReturnCode rc;
  _FC_SeqSlot* seqSlot;
  FC_DataType datatype;
  int numStep;
  void *coordsp;
  double prev,next;
  double movingmin, movingmax;
  int movingminindex, movingmaxindex;
  int onlymono;
  int i, movingmono;

  //default
  if (mono)
    *mono = 0;
  if (min)
    *min = -1;
  if (min_index)
    *min_index = -1;
  if (max)
    *max = -1;
  if (max_index)
    *max_index = -1;

  seqSlot = _fc_getSeqSlot(sequence);
  // check input - they don't have to take 
  if (seqSlot == NULL || !fc_isSequenceValid(sequence)  ||
      (!min && !min_index && !max && !max_index && !mono)){
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
      return FC_INPUT_ERROR;
  }

  fc_printfLogMessage("getting minmax and montonicity of sequence '%s'", 
                      seqSlot->header.name);

  rc = fc_getSequenceInfo (sequence, &numStep, &datatype);
  if (rc!= FC_SUCCESS){
    return rc;
  }

  if (datatype != FC_DT_DOUBLE && datatype !=FC_DT_FLOAT &&
      datatype != FC_DT_INT){
    //non-numerical input sequence
    return FC_INPUT_ERROR;
  }

  if (numStep == 0){
    //zero step seq
    return FC_INPUT_ERROR;
  }

  rc = fc_getSequenceCoordsPtr(sequence,&coordsp);
  if (rc!= FC_SUCCESS){
    //return 0
    return rc;
  }
  switch(datatype){
  case FC_DT_INT:
    prev = (double)(((int*)coordsp)[0]);
    break;
  case FC_DT_FLOAT:
    prev = (double)(((float*)coordsp)[0]);
    break;
  default: //know its a double
    prev = ((double*)coordsp)[0];
    break;
  }
  if (numStep == 1){
    if (min)
      *min = prev;
    if (min_index)
      *min_index = 0;
    if (max)
      *max = prev;
    if (max_index)
      *max_index = 0;
    if (mono)
      *mono = 2;
    return FC_SUCCESS;
  } 
  movingmin = prev;
  movingmax = prev;
  movingminindex = 0;
  movingmaxindex = 0;
  movingmono = 2;

  //at least two steps
  onlymono = (mono && (!min && !min_index && !max && !max_index));

  next = prev;
  for (i = 1; i < numStep; i++){
    prev = next;
    switch(datatype){
    case FC_DT_INT:
      next = (double)(((int*)coordsp)[i]);
      break;
    case FC_DT_FLOAT:
      next = (double)(((float*)coordsp)[i]);
      break;
    default: //know its a double
      next = ((double*)coordsp)[i];
      break;
    }

    //this is inefficient
    if (i == 1){
      //2nd pt
      if (FC_DBL_EQUIV(prev,next)){
        movingmono = 0;
      } else {
        movingmono = (next > prev ? 1: -1);
      }
    } else {
      if (movingmono != 0 && mono){ //not already bad & we actaully care about it
        if (FC_DBL_EQUIV(prev,next)){
          movingmono = 0;
        } else {
          if (next > prev){
            movingmono = (movingmono == 1 ? 1: 0);
          } else {
            movingmono = (movingmono == -1 ? -1: 0);
          }
        }
      }
    }

    if (next <movingmin){
      movingmin = next;
      movingminindex = i;
    }
    if (next > movingmax){
      movingmax = next;
      movingmaxindex = i;
    }  //if equal keep the original one

    //returns if we know we are done

    if (onlymono && (movingmono == 0)){
      //exit if already know its zero
      *mono = 0;
      return FC_SUCCESS;
    }
  }

  //finished

  if (mono)
    *mono = movingmono;
  if (min)
    *min = movingmin;
  if (max)
    *max = movingmax;
  if (min_index)
    *min_index = movingminindex;
  if (max_index)
    *max_index = movingmaxindex;
  return FC_SUCCESS;
}


/**
 * \ingroup Statistics
 * \brief Returns min, max, mean, and std of spacing of a sequence. Returns
 *        indicies of the min and max spacings.
 *
 * \description
 * 
 *        Returns min, max, mean, and std of spacing of a sequence. Returns
 *        indicies of the min and max spacings. Indicies are those of the beginning
 *        of the step (ie. return val of index 1 means the step between 1&2).
 *        If two or more steps have the same spacing and it is the min or
 *        max spacing on the sequence then the index of the first occurence is
 *        the one returned.
 *
 *        Return vals of stats are doubles, regardless of the native datatype
 *        of the sequence (although the calculation is done in the native type
 *        becuase of round off errors). If you want to check afterward for
 *        the regularity of your sequence, you will want to know that min & max differ
 *        by no more than the epsilon of the native type - i.e., if your original
 *        sequence was FC_DT_FLOAT then use FC_FLT_EQUIV even though you are
 *        returned double values.
 *
 *        Nonnumerical sequences return with an error. Sequence must have at least 2 vals.
 *        All invalid return vals are -1. Must have at least one non-NULL output parameter.
 *
 * \todo
 *    - see if want optional args
 *
 * \modifications
 *    - 7/07/05 ACG Created
 */
FC_ReturnCode fc_getSequenceSpacingMinMaxMeanSdev(
  FC_Sequence sequence, /**< input - sequence */
  double *min,       /**< output - min spacing value */
  int *min_index,    /**< output - index at which first min spacing was found */
  double *max,       /**< output - max spacing value */
  int *max_index,     /**< output - index as which first max spacing was found */
  double *mean,       /**< output - mean spacing value */
  double *std       /**< output -  std of the spacing */
  ){
  FC_ReturnCode rc;
  _FC_SeqSlot* seqSlot;
  FC_DataType datatype;
  int numStep;
  void *coordsp;
  int i;

  //default
  if (min)
    *min = -1;
  if (min_index)
    *min_index = -1;
  if (max)
    *max = -1;
  if (max_index)
    *max_index = -1;
  if (mean)
    *mean = -1;
  if (std)
    *std = -1;

  seqSlot = _fc_getSeqSlot(sequence);
  if (seqSlot == NULL || !fc_isSequenceValid(sequence)  ||
      (!min && !min_index && !max && !max_index && !mean && ! std)){
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
      return FC_INPUT_ERROR;
  }

  fc_printfLogMessage("getting spacing minmaxmeanstd of sequence '%s'", 
                      seqSlot->header.name);

  rc = fc_getSequenceInfo (sequence, &numStep, &datatype);
  if (rc!= FC_SUCCESS){
    return rc;
  }

  if (datatype != FC_DT_DOUBLE && datatype !=FC_DT_FLOAT &&
      datatype != FC_DT_INT){
    //non-numerical input sequence
    fc_printfErrorMessage("%s","invalid seq datatype");
    return FC_INPUT_ERROR;
  }

  if (numStep < 2){
    fc_printfErrorMessage("%s","too few steps in seq");
    return FC_ERROR;
  }

  rc = fc_getSequenceCoordsPtr(sequence,&coordsp);
  if (rc!= FC_SUCCESS){
    return rc;
  }

  switch(datatype){
  case FC_DT_INT:
    {
      int *spacings;
      int movmin,movminindex,movmax,movmaxindex, movsum, movsum2;
      
      spacings = (int*)malloc((numStep-1)*sizeof(int));
      
      //first get spacings in native type
      for (i = 0; i < numStep-1; i++){
        spacings[i] = ((int*)coordsp)[i+1] - ((int*)coordsp)[i];
      }
      
      //then do comparisons in native type
      
      movsum = spacings[0];
      movsum2 = spacings[0]*spacings[0];
      movmin = spacings[0];
      movminindex = 0;
      movmax = spacings[0];
      movmaxindex = 0;
      for (i = 1; i < numStep-1; i++){
        if (spacings[i] < movmin){
          movmin = spacings[i];
          movminindex = i;
        }
        if (spacings[i] > movmax){
          movmax = spacings[i];
          movmaxindex = i;
        }
        movsum+=spacings[i];
        movsum2+=spacings[i]*spacings[i];
      }
      if (min)  *min = (double)(movmin);
      if (min_index) *min_index = movminindex;
      if (max) *max = (double)(movmax);
      if (max_index) *max_index = movmaxindex;
      if (mean) *mean = (double)((movsum)/(double)(numStep-1));
      if (std){
        if (numStep-2 > 0) {
          if (movsum2*(numStep-1) == movsum*movsum){
            *std = 0.0;
          }else {
            int temp = (numStep-1)*movsum2 - (movsum * movsum);
            double temp2 = (double)(temp)/(double)((numStep-1)*(numStep-2));
            *std = sqrt(fabs(temp2)); // might be negative if roundoff error
          }
        } else {
          *std = 0.;
        }
      }

      free (spacings);
      return FC_SUCCESS;
    }
    break;
  case FC_DT_FLOAT:
    {
      float *spacings;
      float movmin,movmax, movsum, movsum2;
      int movminindex, movmaxindex;
    
      spacings = (float*)malloc((numStep-1)*sizeof(float));

      //first get spacings in native type
      for (i = 0; i < numStep-1; i++){
        spacings[i] = ((float*)coordsp)[i+1] - ((float*)coordsp)[i];
      }
    
      //then do comparisons in native type
      movsum = spacings[0];
      movsum2 = spacings[0]*spacings[0];
      movmin = spacings[0];
      movminindex = 0;
      movmax = spacings[0];
      movmaxindex = 0;
      for (i = 1; i < numStep-1; i++){
        if (spacings[i] < movmin){
          movmin = spacings[i];
          movminindex = i;
        }
        if (spacings[i] > movmax){
          movmax = spacings[i];
          movmaxindex = i;
        }
        movsum+=spacings[i];
        movsum2+=spacings[i]*spacings[i];
      }
      if (min) *min = (double)(movmin);
      if (min_index) *min_index = movminindex;
      if (max) *max = (double)(movmax);
      if (max_index) *max_index = movmaxindex;
      if (mean) *mean = (double)(movsum)/(float)(numStep-1);
      if (std){
        if (numStep-2 > 0) {
          if (FC_FLT_EQUIV((numStep-1)*movsum2,(movsum*movsum))){
            //      printf("eval std to zero : %10.8g %10.8g\n",(numStep-1)*movsum2,(movsum*movsum));
            *std = 0.0;
          }else {
            float temp = (numStep-1)*movsum2 - (movsum * movsum);
            temp/=(float)((numStep-2)*(numStep-1));
            *std = (double)sqrt(fabs(temp)); // might be negative if roundoff error
          }
        } else {
          *std = 0.;
        }
      }
      free (spacings);
      return FC_SUCCESS;
    }
    break;
  case FC_DT_DOUBLE:
    {
      double *spacings;
      double movmin,movmax, movsum, movsum2;
      int movminindex, movmaxindex;
      
      spacings = (double*)malloc((numStep-1)*sizeof(double));
      
      //first get spacings in native type
      for (i = 0; i < numStep-1; i++){
        spacings[i] = ((double*)coordsp)[i+1] - ((double*)coordsp)[i];
      }
      
      //then do comparisons in native type
      
      movsum = spacings[0];
      movsum2 = spacings[0]*spacings[0];
      movmin = spacings[0];
      movminindex = 0;
      movmax = spacings[0];
      movmaxindex = 0;
      for (i = 1; i < numStep-1; i++){
        if (spacings[i] < movmin){
          movmin = spacings[i];
          movminindex = i;
        }
        if (spacings[i] > movmax){
          movmax = spacings[i];
          movmaxindex = i;
        }
        movsum+=spacings[i];
        movsum2+=spacings[i]*spacings[i];
      }
      if (min) *min = movmin;
      if (min_index) *min_index = movminindex;
      if (max) *max = movmax;
      if (max_index)  *max_index = movmaxindex;
      if (mean) *mean = movsum/(double)(numStep-1);
      if (std){
        if (numStep-2 > 0) {
          if (FC_DBL_EQUIV(movsum2,(movsum*movsum/(float)(numStep-1)))){
            *std = 0.0;
          }else {
            double temp = (numStep-1)*movsum2 - (movsum * movsum);
            temp/=(double)((numStep-2)*(numStep-1));
            *std = sqrt(fabs(temp)); // might be negative if roundoff error
          }
        } else {
          *std = 0.;
        }
      }
      free (spacings);
      return FC_SUCCESS;
    }
    break;
  default:
    fc_printfErrorMessage("%s","WARNING: devloper error ! should never reach this point!");
    return FC_ERROR;
    break;
  }

  //should never reach this point
  return FC_ERROR;
}


/**
 * \ingroup Statistics
 * \brief Given a value and a criteria for comparison (<, <=, =, !=, >=, >) returns the
 *        closest value meeting the conditions of the comparer a sequence has to that
 *        value and its index.
 *
 * \description 
 *
 *       Given a value and a criteria for comparison (<,=,>) returns the
 *       closest value meeting the conditions of the comparer a sequence has to
 *       that value and its index. e.g., use to get the closest value in a
 *       sequence greater than or equal to 3.
 *
 *        Sequence does not have to be ordered. In case of
 *        a tie, returns the first occurence of the value. In case there is
 *        no value matching the criteria or in case of error,
 *        return vals are -1. Return value is of type double regardless
 *        of the type of the orignal sequence.
 *
 * \todo
 *    - this does comparison as doubles. should we want something that treats INTs differently ? 
 *    - = criteria is based on dbl comparison, should there be an epsilon you can pass in ? 
 *      prob not as you can do the <= and >= and get the info you need. if so though can
 *      change this to use fc_leq etc macros.
 *
 * \modifications
 *    - 5/17/05 ACG Created
 *    - 6/27/05 ACG changed to support comparisons other than equals
 */
FC_ReturnCode fc_getClosestSequenceValue(
  FC_Sequence sequence, /**< input - sequence */
  char *comparer, /**<input - character string indicating the comparison criteria */
  double compval,      /**< input - val comparing to */
  double *seqval,       /**< output - seq dataval closest to the compval */
  int *stepIndex           /**< output - index of the seqval */
  ){
  FC_ReturnCode rc;
  _FC_SeqSlot* seqSlot;
  FC_DataType datatype;
  int numStep;
  void *coordsp;
  int holdindex;
  double holdval, currval;
  double holddiff, currdiff;
  int i, icomparer;

  //default
  if (stepIndex)
    *stepIndex = -1;
  if (seqval)
    *seqval = -1;

  fc_printfLogMessage("getting closest seq val for sequence '%s'", 
                      seqSlot->header.name);

  seqSlot = _fc_getSeqSlot(sequence);
  if (seqSlot == NULL){
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
      return FC_INPUT_ERROR;
  }

  if (seqSlot == NULL || !fc_isSequenceValid(sequence)  ||
      !stepIndex || !seqval || !comparer){
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
      return FC_INPUT_ERROR;
  }

  icomparer = -1;
  if (strcmp(comparer,"<") == 0) icomparer = 0;
  if (strcmp(comparer,"<=") == 0) icomparer = 1;
  if (strcmp(comparer,"=") == 0) icomparer = 2;
  if (strcmp(comparer,">=") == 0) icomparer = 3;
  if (strcmp(comparer,">") == 0) icomparer = 4;
  if (strcmp(comparer,"!=") == 0) icomparer = 5;
  if (icomparer ==- 1){
    fc_printfErrorMessage("%s", "invalid comparison criteria");
    return FC_INPUT_ERROR;
  }
 
  rc = fc_getSequenceInfo (sequence, &numStep, &datatype);
  if (rc!= FC_SUCCESS){
    return rc;
  }

  if (datatype != FC_DT_DOUBLE && datatype !=FC_DT_FLOAT &&
      datatype != FC_DT_INT){
    //non-numerical input sequence
    return FC_INPUT_ERROR;
  }

  //ready to go

  rc = fc_getSequenceCoordsPtr(sequence,&coordsp);
  if (rc!= FC_SUCCESS){
    return rc;
  }

  holdindex = -1;
  for (i = 0; i <numStep; i++){
    switch(datatype){
    case FC_DT_INT:
      currval = (double)(((int*)coordsp)[i]);
      break;
    case FC_DT_FLOAT:
      currval = (double)(((float*)coordsp)[i]);
      break;
    default: //know its a double
      currval = ((double*)coordsp)[i];
      break;
    }

    switch (icomparer){ //closest value strictly less than
    case 0:
      //avoid round off error
      if (!FC_DBL_EQUIV(currval,compval) && currval < compval){
        currdiff = fabs(currval-compval);
        if (holdindex == -1 || (holdindex != -1 && currdiff < holddiff)){
          holdindex = i;
          holdval = currval;
          holddiff = currdiff;
        }
      }
    break;
    case 1: //closest value lt or = to
      if (FC_DBL_EQUIV(currval,compval)){
        //call them equal and return
        *stepIndex = i;
        *seqval = currval;
        return FC_SUCCESS;
      }
      if (currval < compval){
        currdiff = fabs(currval-compval);
        if (holdindex == -1 || (holdindex != -1 && currdiff < holddiff)){
          holdindex = i;
          holdval = currval;
          holddiff = currdiff;
        }
      }
      break;
    case 2: //closest value, including possibly = to
      if (FC_DBL_EQUIV(currval,compval)){
        //call them equal and return
        *stepIndex = i;
        *seqval = compval;
        return FC_SUCCESS;
      } else {
        currdiff = fabs(currval-compval);
        if (holdindex == -1 || (holdindex != -1 && currdiff < holddiff)){
          holdindex = i;
          holdval = currval;
          holddiff = currdiff;
        }
      }
      break;
    case 3:  //closest value gt or = to
      if (FC_DBL_EQUIV(currval,compval)){
        //call them equal and return
        *stepIndex = i;
        *seqval = currval;
        return FC_SUCCESS;
      }
      if (currval > compval){
        currdiff = fabs(currval-compval);
        if (holdindex == -1 || (holdindex != -1 && currdiff < holddiff)){
          holdindex = i;
          holdval = currval;
          holddiff = currdiff;
        }
      }
      break;
    case 4: //closest value strictly greater than
      //avoid round off error
      if (!FC_DBL_EQUIV(currval,compval) && currval > compval){
        currdiff = fabs(currval-compval);
        if (holdindex == -1 || (holdindex != -1 && currdiff < holddiff)){
          holdindex = i;
          holdval = currval;
          holddiff = currdiff;
        }
      }
      break;
    case 5: //closest value either < or > but not equal to
      if (!FC_DBL_EQUIV(currval,compval)){
        currdiff = fabs(currval-compval);
        if (holdindex == -1 || (holdindex != -1 && currdiff < holddiff)){
          holdindex = i;
          holdval = currval;
          holddiff = currdiff;
        }
      }
      break;
    default:
      fc_printfErrorMessage("%s","WARNING: devloper error ! should never reach this point!");
      return FC_ERROR;
    }
  }

  if (holdindex != -1){
    *stepIndex = holdindex;
    *seqval = holdval;
  }
  return FC_SUCCESS;
}

//@}
