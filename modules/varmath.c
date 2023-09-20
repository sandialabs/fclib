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
 * \file varmath.c
 * \brief Implementation for \ref VariableMath module.
 *
 * $Source: /home/Repositories/fcdmf/fclib/modules/varmath.c,v $
 * $Revision: 1.84 $ 
 * $Date: 2007/02/28 23:50:29 $
 *
 * \todo Add const versions of functions for seq var and mix of seq var & var. 
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
#include "variable.h"
#include "variableP.h"
#include "util.h"

// this module
#include "varmath.h"
#include "varmathP.h"

/**
 * \addtogroup VariableMath Variable Math
 * \brief Mathematic operations on and between variables.
 *
 * \description
 *
 *    This module provides three sets of functions for doing mathematic
 *    operations on and between variables. In general, these functions always
 *    create new variables which contain the desired result. There are
 *    versions for normal variables and analogous versions for sequence
 *    variables.
 *
 *    <b>1) Built-in operators</b>
 *
 *    The first group of functions, e.g. fc_operatorVar() and
 *    fc_varOperatorVar(), can be used to do built-in operations like '+' and
 *    '*'. Like built-in operations, they preserve or promote type as necessary
 *    with the major expectation that the resultants are only ever int
 *    or double. 
 *
 *    The FCLib function name indicates the type of variables to
 *    which the operator is applied (e.g., Var, SeqVar) and their order
 *    relative to the operator: i.e., fc_seqVarOperatorVar() can be used
 *    to take a seqvar and add a regular var to each of its values at
 *    each of its sequence steps. Thie mimics an operator's usage:
 *    arg1 op arg2. Currently all functions in this group
 *    that take a seqvar as at least of the arguments return a seqvar; all
 *    others return a var.
 *
 *    <b>2) User supplied functions</b>
 *
 *    The second group of functions, e.g. fc_varUnaryFunction() and
 *    fc_varBinaryFunctionVarVar(), allow the user to pass in a function that
 *    does the desired operation. The following example produces a new
 *    variables whose values are the square root of the values of the original
 *    variable ('sqrt' is a standard C math function):
 *
 *    \code 
 * fc_varUnaryFunction(var, sqrt, "sqrt", &newVar);
 *    \endcode
 *
 *    The FCLib function name indicates the type of variables
 *    (e.g., Var, SeqVar) to which the user's function is applied, their
 *     order as arguments to the function, and the return var type. 
 *    i.e., fc_seqVarBinaryFunctionSeqVarVar() can be used to take a 
 *    seqvar and, with a user supplied function to do an add, add a
 *    regular var to it each of its values at each of its
 *    sequence steps and the return value is a seqvar. This mimics
 *    a function's declaration: returntype functionname(arg1,arg2).
 *
 *
 *    <b>(1 & 2 Aside: Binary operations with a "constant")</b>
 *
 *    There are additional versions of the binary versions of these first type
 *    sets of functions that allow operations between a variable and a
 *    'constant'.  That is, a single (but can be multi-component) data point
 *    can be applied to all data points of the variable.  For example, the
 *    following scales a 3D vector by a factor of two:
 *    
 *    \code
 * int scale[3] = { 2, 2, 2 };
 * fc_varOperatorConst(vectorVar, '*', 3, scale, FC_DT_INT, "scaled vector", &newVar);
 *    \endcode
 *
 *    <b>3) Other</b>
 *
 *    The final set of functions provide operations that don't fit into
 *    the previous twos. Mostly these are functions that treat the different
 *    components of multicomponent variables differently. 
 *
 *
 * \modifications
 *    - 2003-JUN-23  W. Koegler  Created.
 *    - 02/02/04 RM, Redone, new implementation, added unary operators and 
 *                   general functions for unary and binary operations.
 *    - 02/23/04 RM, Added fc_varOperatorConst() 
 *    - 02/26/04 RM, Added user operator functions and examples for 
 *                   them, and functions for more advanced users 
 *    - 03/04/04 RM, Added examples of unary and binary user functions.
 *    - 11/10/04 WSK, reorg & rename. Also removed 'more sophisticated' unary &
 *                    binary functions as they were too general purpose.
 *                    (Might as well just write an entire new function.)
 *    - 6/03/05 ACG, adding analgous functions but for seqvar
 */

/**
 * \ingroup  VariableMath
 * \defgroup  PrivateVariableMath (Private)
 *
 * \description
 *
 *    <b>Helper functions to perform built-in operations on data arrays.</b>
 *
 *       - _fc_negateInts()
 *       - _fc_negateDoubles()
 *       - _fc_addInts()
 *       - _fc_subtractInts()
 *       - _fc_multiplyInts()
 *       - _fc_divideInts()
 *       - _fc_addDoubles()
 *       - _fc_subtractDoubles()
 *       - _fc_multiplyDoubles()
 *       - _fc_divideDoubles()
 *   
 *    The motivation for using these functions (as opposed to using
 *    fc_varUnaryFunction() and similar functions) are: 1) it is more
 *    efficient to do built-in operation over entire array than to embed the
 *    operation in a function that is called repeatedly; and 2) we are also
 *    using them help to retain type when appropriate (i.e. if you add two
 *    ints, you get an int).
 *
 *    All of the built-in helper functions assume that the calling program has
 *    allocated room for the results.
 *
 *    The binary built-in helper functions operate on interleaved data
 *    arrays--that is if you have numData data points and numComp components,
 *    the jth component of the ith data point is found by data[i*numComp+j].
 *    All data arrays are assume to have the same number of components.  If
 *    they also have the same number of data points, the operation is performed
 *    for each data value. A special case is allowed for when one of the data
 *    arrays has just one data point.  In that case, that single data point is
 *    used for all data points of the other array.
 */

/** \name Built in operators on variables. */
//-------------------------------
//@{

/**
 * \ingroup  VariableMath
 * \brief Performs built-in unary operations on a variable.
 *
 * \description
 *
 *    Given a string that represents the desired unary operation, return
 *    a new variable with the result of that operation performed on the
 *    given variable.
 *
 *    Supported operations are:
 *      - <tt> "-"  =>  -var </tt>
 *
 *    The new variable will be on the same mesh as the original and will
 *    have the same association and math type as the original. The operation
 *    is performed on each component of each data point. 
 *  
 *    If the original variable had data type FC_DT_INT, then the result
 *    will have data type FC_DT_INT. If the original variable had data type
 *    FC_DT_FLOAT or FC_DT_DOUBLE, the result will have data type FC_DT_DOUBLE.
 *    The function returns with an error if the data type is FC_DT_CHAR.
 *
 * \modifications 
 *    - 02/02/04  R.McCoy  Created.
 *    - 11/09/03  WSK, restricted to built-in types.
 */
FC_ReturnCode fc_operatorVar(
  FC_Variable var,      /**< input - variable */
  char* operation,      /**< input - operation to be performed on variable */
  char* new_var_name,   /**< input - name for the new variable */
  FC_Variable* new_var  /**< output - the new variable */
) {
  FC_ReturnCode rc;
  int i;
  int opID;
  int numData, numComp;
  FC_AssociationType assoc, new_datatype;
  FC_DataType datatype;
  FC_MathType mathtype;
  void *data, *resultPtr;
  double *dblDataPtr;
  // create array of pointers to int & double array functions
  FC_ReturnCode (*IntUnaryOper[])(int, int*, int*) = { 
                         _fc_negateInts
                        };
  FC_ReturnCode (*DblUnaryOper[])(int, double*, double*) = {
                          _fc_negateDoubles 
                        }; 

  // set default return value
  if (new_var)
    *new_var = FC_NULL_VARIABLE;
  
  // check input
  if (!(fc_isVariableValid(var)) || !operation || !new_var || !new_var_name) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // pick operation
  if (!strcmp(operation, "-")) 
    opID = 0;
  else {
    fc_printfErrorMessage("operation '%s' not recognized.", operation); 
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Performing unary operation '%s' on variable '%s'",
            operation, _fc_getVarSlot(var)->header.name);

  // get variable info
  rc = fc_getVariableInfo(var, &numData, &numComp, &assoc, 
			  &mathtype, &datatype);
  if (rc != FC_SUCCESS)
    return rc;

  // check that the datatype can be handled here
  if (datatype != FC_DT_INT && datatype != FC_DT_FLOAT && 
      datatype != FC_DT_DOUBLE) {
    fc_printfErrorMessage("data type '%s' cannot be handled.", 
			  fc_getDataTypeText(datatype));
    return FC_ERROR;
  } 

  // get data
  rc = fc_getVariableDataPtr(var, &data);
  if (rc != FC_SUCCESS)
    return rc;

  if (datatype == FC_DT_INT) {
    new_datatype = FC_DT_INT;
    resultPtr = (int*)malloc(sizeof(int)*numComp*numData);
    if (resultPtr == NULL) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));  
      return FC_MEMORY_ERROR;
    }
    rc = (*IntUnaryOper[opID])(numComp*numData, (int*)data, (int*)resultPtr);
    if (rc != FC_SUCCESS)
      return rc;
  }
  else {
    new_datatype = FC_DT_DOUBLE;

     // get double representation of data
    if (datatype == FC_DT_DOUBLE)
      dblDataPtr = (double*)data;
    else if (datatype == FC_DT_FLOAT) {
      float* cast_data = (float*)data;
      dblDataPtr = malloc(sizeof(double)*numComp*numData);
      if (dblDataPtr == NULL) {
	fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));  
	return FC_MEMORY_ERROR;
      }
      for (i = 0; i < numComp*numData; i++)
	dblDataPtr[i] = cast_data[i];
    }

    // do it
    resultPtr = malloc(sizeof(double)*numComp*numData);
    if (resultPtr == NULL) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));  
      return FC_MEMORY_ERROR;
    }
    rc = (*DblUnaryOper[opID])(numComp*numData, (double*)dblDataPtr, 
			       (double*)resultPtr);
    if (rc != FC_SUCCESS)
      return rc;
    
    // deallocate memory if it was allocated
    if (datatype != FC_DT_DOUBLE)
      free(dblDataPtr);
  }
 
  // create the new variable
  if (fc_isVariableGlobal(var)) {
    FC_Dataset dataset;
    rc = fc_getDatasetFromVariable(var, &dataset);
    if (rc != FC_SUCCESS)
      return rc;
    rc = fc_createGlobalVariable(dataset, new_var_name, new_var);
    if (rc != FC_SUCCESS)
      return rc; 
  }
  else {
    FC_Mesh mesh;
    rc = fc_getMeshFromVariable(var, &mesh);
    if (rc != FC_SUCCESS)
      return rc;
    rc = fc_createVariable(mesh, new_var_name, new_var);
    if (rc != FC_SUCCESS)
      return rc; 
  }
  rc = fc_setVariableDataPtr(*new_var, numData, numComp, assoc, mathtype, 
			     new_datatype, resultPtr);
  return rc;    

} 

/**
 * \ingroup  VariableMath
 * \brief Performs built-in binary operation on two variables.
 *
 * \description 
 *
 *   Given a string that represents the desired binary operation, return a new
 *   variable with the result of that operation performed on the the given
 *   variables.
 *
 *    Supported operations are:
 *     - <tt> "+"  =>  var1 + var2 </tt>
 *     - <tt> "-"  =>  var1 - var2 </tt>
 *     - <tt> "*"  =>  var1 * var2 </tt>
 *     - <tt> "/"  =>  var1 / var2 </tt>
 * 
 *    The two variables must have the same number of data points and the same
 *    number of components.
 *
 *    The new variable will be on the same mesh and have the same association
 *    and mathtype as the first variable. The operation is performed per
 *    component and data point.
 *
 *    If both original variables have the data type FC_DT_INT, then the 
 *    operation will be an integer operation and the result
 *    will have data type FC_DT_INT. Otherwise, all variable values are
 *    promoted to double and the operation will be a double operation and
 *    the result will have data type FC_DT_DOUBLE.
 *
 * \modifications 
 *    - 01/29/04  R.McCoy  Created.
 *    - 11/09/03  WSK, restricted to built-in types.
 */
FC_ReturnCode fc_varOperatorVar(
  FC_Variable var1,     /**< input - variable 1  */
  char* operation,      /**< input - operation to be performed on variables */
  FC_Variable var2,     /**< input - variable 2  */
  char* new_var_name,   /**< input - name for the new variable */
  FC_Variable* new_var  /**< output - the new variable */
) {        
  FC_ReturnCode rc;
  int numData1, numData2, numComp1, numComp2;
  FC_DataType datatype1, datatype2;
  void *data1, *data2;
   
  // set default return value
  if (new_var)
    *new_var = FC_NULL_VARIABLE;

  // check input
  if (!(fc_isVariableValid(var1)) || !operation || 
      !(fc_isVariableValid(var2)) || !new_var_name || !new_var) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));   
    return FC_INPUT_ERROR;
  }
 
  // log message
  fc_printfLogMessage("Performing operation (var1 %s var2), with var1 = '%s', "
           "var2 = '%s'" , operation, _fc_getVarSlot(var1)->header.name,
           _fc_getVarSlot(var2)->header.name);  

  // get variable info
  rc = fc_getVariableInfo(var1, &numData1, &numComp1, NULL, NULL, &datatype1);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_getVariableInfo(var2, &numData2, &numComp2, NULL, NULL, &datatype2);
  if (rc != FC_SUCCESS)
    return rc;

  // check that variables are similar enough
  if (numData1 != numData2 || numComp1 != numComp2) {
    fc_printfErrorMessage("Vars need to have same number of data points "
			  "and same number of components");
    return FC_INPUT_ERROR;
  }
  if ((fc_isVariableGlobal(var1) && !fc_isVariableGlobal(var2)) ||
      (!fc_isVariableGlobal(var1) && fc_isVariableGlobal(var2))) {
    fc_printfErrorMessage("Cannot process a mix of global and non-global vars");
    return FC_INPUT_ERROR;
  }

  // get data
  rc = fc_getVariableDataPtr(var1, &data1);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_getVariableDataPtr(var2, &data2);
  if (rc != FC_SUCCESS)
    return rc;

  // do it
  return _fc_binaryOperatorCore(numComp1, numData1, data1, datatype1,
				operation, numData2, data2, datatype2,
				var1, new_var_name, new_var);
}

/**
 * \ingroup  VariableMath
 * \brief Performs built-in binary operation on a variable and a constant. 
 *
 * \description  
 *
 *   Given a string that represents the desired binary operation, return a new
 *   variable with the result of that operation performed with the constant
 *   on each data point of the variable. Use fc_constOperatorVar() to apply
 *   operation in reverse order.
 *
 *    Supported operations are:
 *     - <tt> "+"  =>  var + const </tt>
 *     - <tt> "-"  =>  var - const </tt>
 *     - <tt> "*"  =>  var * const </tt>
 *     - <tt> "/"  =>  var / const </tt>
 * 
 *    The variables and the constant must have the same number of components.
 *
 *    The new variable will be on the same mesh and have the same association
 *    and mathtype as the variable. The operation is performed per component
 *    and data point.
 *
 *    If both the variable and the constant have the data type FC_DT_INT, then
 *    the operation will be an integer operation and the result will have data
 *    type FC_DT_INT. Otherwise, all values are promoted to double and the
 *    operation will be a double operation and the result will have data type
 *    FC_DT_DOUBLE.
 *
 * \modifications 
 *    - 02/19/04  R.McCoy  Created.
 *    - 11/09/03  WSK, restricted to built-in types.
 */
FC_ReturnCode fc_varOperatorConst(
  FC_Variable var,       /**< input - FCLib variable */
  char* operation,       /**< input - operation to be performed on variables */
  int constNumComp,      /**< input - number of components in constVector */
  void* constVector,     /**< input - constant vector  */
  FC_DataType constType, /**< input - data type of constVector*/
  char* new_var_name,    /**< input - name for the new variable */
  FC_Variable* new_var   /**< output - the new variable */
) {
  FC_ReturnCode rc;
  int numComp, numData;
  FC_DataType datatype; 
  void *data; 
   
  // set default return value
  if (new_var)
    *new_var = FC_NULL_VARIABLE;
    
  // check input
  if (!(fc_isVariableValid(var)) || !operation || constNumComp < 1 ||
      !constVector || !fc_isDataTypeValid(constType) || 
      constType == FC_DT_UNKNOWN || !new_var_name  || !new_var) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));   
    return FC_INPUT_ERROR;
  }

  // log message
  fc_printfLogMessage("Performing operation (var %s const), with var '%s'",
		      operation, _fc_getVarSlot(var)->header.name);

  // get variable info
  rc = fc_getVariableInfo(var, &numData, &numComp, NULL, NULL, &datatype);
  if (rc != FC_SUCCESS)
    return rc;

  // check that const is similar enough
  if (numComp != constNumComp) {
    fc_printfErrorMessage("constant need sot have same number of components "
			  "as var");
    return FC_INPUT_ERROR;
  }

  // get data
  rc = fc_getVariableDataPtr(var, &data);
  if (rc != FC_SUCCESS)
    return rc;

  // do it
  return _fc_binaryOperatorCore(numComp, numData, data, datatype,
				operation, 1, constVector, constType,
				var, new_var_name, new_var);
}

/**
 * \ingroup  VariableMath
 * \brief Performs built-in binary operation on a constant and a variable. 
 *
 * \description
 * 
 *   Given a string that represents the desired binary operation, return a new
 *   variable with the result of that operation performed with the constant
 *   on each data point of the variable. Use fc_varOperatorConst() to apply
 *   operation in reverse order.
 *
 *    Supported operations are:
 *     - <tt> "+"  =>  const + var </tt>
 *     - <tt> "-"  =>  const - var </tt>
 *     - <tt> "*"  =>  const * var </tt>
 *     - <tt> "/"  =>  const / var </tt>
 * 
 *    The variables and the constant must have the same number of components.
 *
 *    The new variable will be on the same mesh and have the same association
 *    and mathtype as the variable. The operation is performed per component
 *    and data point.
 *
 *    If both the variable and the constant have the data type FC_DT_INT, then
 *    the operation will be an integer operation and the result will have data
 *    type FC_DT_INT. Otherwise, all values are promoted to double and the
 *    operation will be a double operation and the result will have data type
 *    FC_DT_DOUBLE.
 *
 * \modifications 
 *    - 11/09/03  WSK, created to complement fc_varOperatorConst.
 */
FC_ReturnCode fc_constOperatorVar(
  int constNumComp,      /**< input - number of components in constVector */
  void* constVector,     /**< input - constant vector  */
  FC_DataType constType, /**< input - data type of constVector*/
  char* operation,       /**< input - operation to be performed on variables */
  FC_Variable var,       /**< input - FCLib variable */
  char* new_var_name,    /**< input - name for the new variable */
  FC_Variable* new_var   /**< output - the new variable */
) {
  FC_ReturnCode rc;
  int numComp, numData;
  FC_DataType datatype;
  void *data; 
   
  // set default return value
  if (new_var)
    *new_var = FC_NULL_VARIABLE;
    
  // check input
  if (!(fc_isVariableValid(var)) || !operation || constNumComp < 1 ||
      !constVector || !fc_isDataTypeValid(constType) || 
      constType == FC_DT_UNKNOWN || !new_var_name  || !new_var) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));   
    return FC_INPUT_ERROR;
  }

  // log message
  fc_printfLogMessage("Performing operation (var %s const), with var '%s'",
		      operation, _fc_getVarSlot(var)->header.name);

  // get variable info
  rc = fc_getVariableInfo(var, &numData, &numComp, NULL, NULL, &datatype);
  if (rc != FC_SUCCESS)
    return rc;

  // check that const is similar enough
  if (numComp != constNumComp) {
    fc_printfErrorMessage("constant need sot have same number of components "
			  "as var");
    return FC_INPUT_ERROR;
  }

  // get data
  rc = fc_getVariableDataPtr(var, &data);
  if (rc != FC_SUCCESS)
    return rc;
 
  return _fc_binaryOperatorCore(numComp, 1, constVector, constType, operation,
				numData, data, datatype,
				var, new_var_name, new_var);
}

//@}

/** \name Built-in operators on sequence variables. */
//-------------------------------
//@{

/**
 * \ingroup  VariableMath
 * \brief Performs built-in unary operations on a sequence variable.
 *
 * \description
 *
 *    Given a string that represents the desired unary operation, return
 *    a new sequence variable with the result of that operation
 *    performed on each step of the given sequence variable. This is the
 *    sequence variable analogy to fc_operatorVar();
 *
 *    Supported operations are:
 *      - <tt> "-"  =>  -var </tt>
 *
 *    The new sequence variable will be on the same mesh and sequence as the 
 *    original and will have the same association and math type as
 *    the original. The operation is performed on each component of each 
 *    data point at each sequence step. 
 *  
 *    If the original sequence variable had data type FC_DT_INT, 
 *    then the result
 *    will have data type FC_DT_INT. If the original sequence variable had
 *    data type FC_DT_FLOAT or FC_DT_DOUBLE, the result will have data
 *    type FC_DT_DOUBLE. The function returns with an error if the data
 *    type is FC_DT_CHAR.
 *
 * \modifications 
 *    - 06/03/04  ACG  Created.
 *    - 06/07/05  ACG removing datatype check and leaving that to internal call
 *    - 01/24/06  WSK Rewrote innards to mimimize data creation.
 */
FC_ReturnCode fc_operatorSeqVar(
  int numStep, 	/**< input - numstep in the input seq var */			
  FC_Variable *seq_var,      /**< input - seq variable */
  char* operation,      /**< input - operation to be performed on variable */
  char* new_seq_var_name,   /**< input - name for the new seq variable */
  FC_Variable** new_seq_var  /**< output - the new seq variable */
) {
  FC_ReturnCode rc;
  FC_Sequence seq;
  FC_Variable* temp_vars;
  int i, j;

  // default return
  if (new_seq_var)
    *new_seq_var = NULL;

  // check input
  if ((!fc_isSeqVariableValid(numStep,seq_var)) || (new_seq_var == NULL)
      || new_seq_var_name == NULL){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  //other checks are done internally in operatorvar

  // log message
  fc_printfLogMessage("Performing operatorSeqVar '%s' on variable '%s'",
            operation, _fc_getVarSlot(seq_var[0])->header.name);

  // Get info
  rc = fc_getSequenceFromSeqVariable(numStep,seq_var,&seq);
  if (rc != FC_SUCCESS)
    return rc;

  // Make space to save temp_vars
  temp_vars = malloc(numStep*sizeof(FC_Variable));
  if (!temp_vars) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  
  // Do calculation on each step of the seq variables
  for (i = 0; i < numStep; i++) {
    rc = fc_operatorVar(seq_var[i], operation, "junk_var", &temp_vars[i]);
    if (rc != FC_SUCCESS) {
      for (j = 0; j < i-1; j++)
	fc_deleteVariable(temp_vars[j]);
      fc_printfErrorMessage("Operation failed on step %d", i);
      free(temp_vars);
      return rc;
    }
  }

  // stitch steps into the new seq var
  rc = fc_convertVariablesToSeqVariable(numStep, temp_vars, seq,
				     new_seq_var_name, new_seq_var);
  free(temp_vars);
  return rc;
}

/**
 * \ingroup  VariableMath
 * \brief Performs built-in binary operation on two seq variables.
 *
 * \description 
 *
 *   Given a string that represents the desired binary operation, return a new
 *   seq variable with the result of that operation performed on the the given
 *   seq variables. This is the seq var version of fc_varOperatorVar.
 *
 *    Supported operations are:
 *     - <tt> "+"  =>  var1 + var2 </tt>
 *     - <tt> "-"  =>  var1 - var2 </tt>
 *     - <tt> "*"  =>  var1 * var2 </tt>
 *     - <tt> "/"  =>  var1 / var2 </tt>
 * 
 *    The two seq variables must have the same number of data points and the same
 *    number of components and be on the same sequence.
 *
 *    The seq new variable will be on the same mesh and have the same association
 *    and mathtype as the first variable. The operation is performed per
 *    component and data point at each sequence step.
 *
 *    If both original variables have the data type FC_DT_INT, then the 
 *    operation will be an integer operation and the result
 *    will have data type FC_DT_INT. Otherwise, all variable values are
 *    promoted to double and the operation will be a double operation and
 *    the result will have data type FC_DT_DOUBLE.
 *
 * \modifications 
 *    - 06/07/05  ACG  Created.
 *    - 01/24/06  WSK Rewrote innards to mimimize data creation.
 */
FC_ReturnCode fc_seqVarOperatorSeqVar(
  int numStep, /**< input - num step in input seq var */
  FC_Variable *seq_var1, /**< input seq var1 */
  char* operation,      /**< input - operation to be performed on variables */
  FC_Variable *seq_var2, /**< input seq var2 */
  char* new_seq_var_name,   /**< input - name for the new variable */ 
  FC_Variable** new_seq_var  /**< output - the new seq variable */
) {   
  FC_ReturnCode rc;
  FC_Sequence seq1, seq2;
  FC_Variable* temp_vars;
  int i, j;

  // default return
  if (new_seq_var)
    *new_seq_var = NULL;

  // check input (this indirectly checks that numStep is the same on both seqs
  if ((!fc_isSeqVariableValid(numStep,seq_var1)) ||
      (!fc_isSeqVariableValid(numStep,seq_var2)) ||
      (new_seq_var == NULL) || new_seq_var_name == NULL){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  //other checks are done internally in varoperatorvar

  // log message
  fc_printfLogMessage("Performing SeqVarOperatorSeqVar '%s' on"
		      "seq variables '%s' and '%s'", operation,
		      _fc_getVarSlot(seq_var1[0])->header.name,
		      _fc_getVarSlot(seq_var2[0])->header.name);

  // get info
  rc = fc_getSequenceFromSeqVariable(numStep,seq_var1,&seq1);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_getSequenceFromSeqVariable(numStep,seq_var2,&seq2);
  if (rc != FC_SUCCESS)
    return rc;
  // FIX? remove this?
  if (!(FC_HANDLE_EQUIV(seq1,seq2))){
    fc_printfErrorMessage("%s", "Variables are not on the same sequence");
    return FC_INPUT_ERROR;
  }

  // Make space to save temp_vars
  temp_vars = malloc(numStep*sizeof(FC_Variable));
  if (!temp_vars) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }

  // Do calculation on each step of the seq variables
  for (i = 0; i < numStep; i++) {
    rc = fc_varOperatorVar(seq_var1[i], operation, seq_var2[i], "junk_var",
			   &temp_vars[i]); 
    if (rc != FC_SUCCESS) { 
      for (j = 0; j < i-1; j++)
	fc_deleteVariable(temp_vars[j]);
      free(temp_vars);
      fc_printfErrorMessage("Operation failed on step %d", i);
      return rc;
    }
  }

  // stitch steps into the new seq var
  rc = fc_convertVariablesToSeqVariable(numStep, temp_vars, seq1,
				     new_seq_var_name, new_seq_var);
  free(temp_vars);
  return rc;
}

/**
 * \ingroup  VariableMath
 * \brief Performs built-in binary operation on a seq variable a seq of constants.
 *
 * \description 
 *
 *   Given a string that represents the desired binary operation, return a new
 *   seq variable with the result of that operation performed on the the given
 *   seq of constants. Use \ref fc_seqConstOperatorSeqVar() to apply the
 *   operation in the reverse order.
 *
 *   This is the seq var version of \ref fc_varOperatorConst() where a
 *   different constant vector is applied to each step of the sequence
 *   variable. The constant values should be an array with all components for
 *   the first step first, then all components for the next step, etc.
 *   If you would like to apply the same constant to all steps, use
 *   fc_seqVarOperatorConst(). 
 *
 *    Supported operations are:
 *     - <tt> "+"  =>  var + const </tt>
 *     - <tt> "-"  =>  var - const </tt>
 *     - <tt> "*"  =>  var * const </tt>
 *     - <tt> "/"  =>  var / const </tt>
 * 
 *    The seq variable and the constants must have the same number components.
 *
 *    The seq new variable will be on the same mesh and have the same association
 *    and mathtype as the seq variable. The operation is performed per
 *    component and data point at each sequence step.
 *
 *    If both the original variable and the constants have the data type
 *    FC_DT_INT, then the operation will be an integer operation and the result
 *    will have data type FC_DT_INT. Otherwise, all variable values are
 *    promoted to double and the operation will be a double operation and the
 *    result will have data type FC_DT_DOUBLE.
 *
 * \modifications 
 *    - 11/20/2006 WSD Created.
 */
FC_ReturnCode fc_seqVarOperatorSeqConst(
  int numStep,           /**< input - num step in input seq var and const */
  FC_Variable *seq_var1, /**< input seq var1 */
  char* operation,       /**< input - operation to be performed on variables */
  int constNumComp,      /**< input - number of components in constVector */
  void* constVectors,    /**< input - constant vectors, stored as vector for
                          first step, then vector for next step, etc. */
  FC_DataType constType, /**< input - data type of constVector */
  char* new_seq_var_name,   /**< input - name for the new variable */ 
  FC_Variable** new_seq_var  /**< output - the new seq variable */
) {   
  FC_ReturnCode rc;
  FC_Sequence seq;
  FC_Variable* temp_vars;
  int i, j;

  // default return
  if (new_seq_var)
    *new_seq_var = NULL;

  // check input (this indirectly checks that numStep is the same on both seqs
  if (!fc_isSeqVariableValid(numStep, seq_var1) || !operation || 
      constNumComp < 1 || !constVectors || !fc_isDataTypeValid(constType) ||
      constType == FC_DT_UNKNOWN || new_seq_var_name == NULL ||
      new_seq_var == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  //other checks are done internally in varoperatorconst

  // log message
  fc_printfLogMessage("Performing SeqVarOperatorSeqConst '%s' with"
		      "seq variable '%s'", operation,
		      _fc_getVarSlot(seq_var1[0])->header.name);

  // get info
  rc = fc_getSequenceFromSeqVariable(numStep, seq_var1, &seq);
  if (rc != FC_SUCCESS)
    return rc;

  // Make space to save temp_vars
  temp_vars = malloc(numStep*sizeof(FC_Variable));
  if (!temp_vars) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }

  // Do calculation on each step of the seq variables
  for (i = 0; i < numStep; i++) {
    void* temp_constVector;
    switch(constType) {
    case FC_DT_CHAR:
      temp_constVector = &((char*)constVectors)[i*constNumComp]; break;
    case  FC_DT_INT: 
      temp_constVector = &((int*)constVectors)[i*constNumComp]; break;
    case FC_DT_FLOAT:
      temp_constVector = &((float*)constVectors)[i*constNumComp]; break;
    case FC_DT_DOUBLE:
      temp_constVector = &((double*)constVectors)[i*constNumComp]; break;
    case FC_DT_UNKNOWN:
      fc_printfErrorMessage("should never get here");
      return FC_ERROR;
    }
    rc = fc_varOperatorConst(seq_var1[i], operation, constNumComp, 
                            temp_constVector, constType,
                             "junk_var", &temp_vars[i]); 
    if (rc != FC_SUCCESS) { 
      for (j = 0; j < i-1; j++)
	fc_deleteVariable(temp_vars[j]);
      free(temp_vars);
      fc_printfErrorMessage("Operation failed on step %d", i);
      return rc;
    }
  }

  // stitch steps into the new seq var
  rc = fc_convertVariablesToSeqVariable(numStep, temp_vars, seq,
				     new_seq_var_name, new_seq_var);
  free(temp_vars);
  return rc;
}

/**
 * \ingroup  VariableMath
 * \brief Performs built-in binary operation on a seq variable a seq of constants.
 *
 * \description 
 *
 *   Given a string that represents the desired binary operation, return a new
 *   seq variable with the result of that operation performed on the the given
 *   seq of constants. Use \ref fc_seqVarOperatorSeqConst() to apply the
 *   operation in the reverse order.
 *
 *   This is the seq var version of \ref fc_constOperatorVar() where a
 *   different constant vector is applied to each step of the sequence
 *   variable. The constant values should be an array with all components for
 *   the first step first, then all components for the next step, etc.
 *   If you would like to apply the same constant to all steps, use
 *   fc_constOperatorSeqVar(). 
 *
 *    Supported operations are:
 *     - <tt> "+"  =>  const + var </tt>
 *     - <tt> "-"  =>  const - var </tt>
 *     - <tt> "*"  =>  const * var </tt>
 *     - <tt> "/"  =>  const / var </tt>
 * 
 *    The seq variable and the constants must have the same number components.
 *
 *    The seq new variable will be on the same mesh and have the same association
 *    and mathtype as the seq variable. The operation is performed per
 *    component and data point at each sequence step.
 *
 *    If both the original variable and the constants have the data type
 *    FC_DT_INT, then the operation will be an integer operation and the result
 *    will have data type FC_DT_INT. Otherwise, all variable values are
 *    promoted to double and the operation will be a double operation and the
 *    result will have data type FC_DT_DOUBLE.
 *
 * \modifications 
 *    - 11/20/2006 WSD Created.
 */
FC_ReturnCode fc_seqConstOperatorSeqVar(
  int numStep,           /**< input - num step in input seq var and const */
  int constNumComp,      /**< input - number of components in constVector */
  void* constVectors,    /**< input - constant vectors, stored as vector for
                          first step, then vector for next step, etc. */
  FC_DataType constType, /**< input - data type of constVector */
  char* operation,       /**< input - operation to be performed on variables */
  FC_Variable *seq_var2, /**< input seq var2 */
  char* new_seq_var_name,   /**< input - name for the new variable */ 
  FC_Variable** new_seq_var  /**< output - the new seq variable */
) {   
  FC_ReturnCode rc;
  FC_Sequence seq;
  FC_Variable* temp_vars;
  int i, j;

  // default return
  if (new_seq_var)
    *new_seq_var = NULL;

  // check input (this indirectly checks that numStep is the same on both seqs
  if (!fc_isSeqVariableValid(numStep, seq_var2) || !operation || 
      constNumComp < 1 || !constVectors || !fc_isDataTypeValid(constType) ||
      constType == FC_DT_UNKNOWN || new_seq_var_name == NULL ||
      new_seq_var == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  //other checks are done internally in varoperatorconst

  // log message
  fc_printfLogMessage("Performing SeqVarOperatorSeqConst '%s' with"
		      "seq variable '%s'", operation,
		      _fc_getVarSlot(seq_var2[0])->header.name);

  // get info
  rc = fc_getSequenceFromSeqVariable(numStep, seq_var2, &seq);
  if (rc != FC_SUCCESS)
    return rc;

  // Make space to save temp_vars
  temp_vars = malloc(numStep*sizeof(FC_Variable));
  if (!temp_vars) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }

  // Do calculation on each step of the seq variables
  for (i = 0; i < numStep; i++) {
    void* temp_constVector;
    switch(constType) {
    case FC_DT_CHAR:
      temp_constVector = &((char*)constVectors)[i*constNumComp]; break;
    case  FC_DT_INT: 
      temp_constVector = &((int*)constVectors)[i*constNumComp]; break;
    case FC_DT_FLOAT:
      temp_constVector = &((float*)constVectors)[i*constNumComp]; break;
    case FC_DT_DOUBLE:
      temp_constVector = &((double*)constVectors)[i*constNumComp]; break;
    case FC_DT_UNKNOWN:
      fc_printfErrorMessage("should never get here");
      return FC_ERROR;
    }
    rc = fc_constOperatorVar(constNumComp, temp_constVector, constType,
                             operation, seq_var2[i], 
                             "junk_var", &temp_vars[i]); 
    if (rc != FC_SUCCESS) { 
      for (j = 0; j < i-1; j++)
	fc_deleteVariable(temp_vars[j]);
      free(temp_vars);
      fc_printfErrorMessage("Operation failed on step %d", i);
      return rc;
    }
  }

  // stitch steps into the new seq var
  rc = fc_convertVariablesToSeqVariable(numStep, temp_vars, seq,
				     new_seq_var_name, new_seq_var);
  free(temp_vars);
  return rc;
}

/** \name Built-in operators between variables and sequence variables. */
//-------------------------------
//@{

/**
 * \ingroup  VariableMath
 * \brief Performs built-in binary operation on a seq variable and a normal var.
 *
 * \description 
 *
 *   Given a string that represents the desired binary operation, return a new
 *   seq variable with the result of that operation between a given
 *   seq variable and a normal var applied at each sequence step.
 *   For example, one can use this to scale a seqvar's values by their mean over
 *   the sequence (e.g., time).
 *
 *    Supported operations are:
 *     - <tt> "+"  =>  var1 + var2 </tt>
 *     - <tt> "-"  =>  var1 - var2 </tt>
 *     - <tt> "*"  =>  var1 * var2 </tt>
 *     - <tt> "/"  =>  var1 / var2 </tt>
 * 
 *    The two seq variable and the normal variable must have the same number of data
 *    points and the same number of components.
 *
 *    The seq new variable will be on the same mesh and sequence and will have the
 *    same association and mathtype as the seq variable. The operation is 
 *    performed per component and data point at each sequence step.
 *
 *    If both original variables have the data type FC_DT_INT, then the 
 *    operation will be an integer operation and the result
 *    will have data type FC_DT_INT. Otherwise, all variable values are
 *    promoted to double and the operation will be a double operation and
 *    the result will have data type FC_DT_DOUBLE.
 *    - 01/24/06  WSK Rewrote innards to mimimize data creation.
 *
 * \modifications 
 *    - 06/07/05  ACG  Created.
 *    - 01/24/06  WSK Rewrote innards to mimimize data creation.
 */
FC_ReturnCode fc_seqVarOperatorVar(
  int numStep, /**< input - num step in input seq var */
  FC_Variable *seq_var1, /**< input seq var1 */
  char* operation,      /**< input - operation to be performed on variables */
  FC_Variable var2, /**< input var */
  char* new_seq_var_name,   /**< input - name for the new variable */ 
  FC_Variable** new_seq_var  /**< output - the new seq variable */
) {   
  FC_ReturnCode rc;
  FC_Sequence seq;
  FC_Variable* temp_vars;
  int i, j;

  // default return
  if (new_seq_var)
    *new_seq_var = NULL;

  // check input
  if ((!fc_isSeqVariableValid(numStep, seq_var1)) ||
      (!fc_isVariableValid(var2)) ||
      (new_seq_var == NULL) || new_seq_var_name == NULL){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  //other checks are done internally in varoperatorvar

  // log message
  fc_printfLogMessage("Performing operation '%s' on seq variable '%s'"
		      " and variable '%s'", operation,
		      _fc_getVarSlot(seq_var1[0])->header.name,
		      _fc_getVarSlot(var2)->header.name);

  // Get info
  rc = fc_getSequenceFromSeqVariable(numStep, seq_var1, &seq);
  if (rc != FC_SUCCESS)
    return rc;

  // Make space to save temp_vars
  temp_vars = malloc(numStep*sizeof(FC_Variable));
  if (!temp_vars) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }

  // Do calculation on each step of the seq variables
  for (i = 0; i < numStep; i++) {
    rc = fc_varOperatorVar(seq_var1[i], operation, var2, "junk_var",
			   &temp_vars[i]); 
    if (rc != FC_SUCCESS) { 
      for (j = 0; j < i-1; j++)
	fc_deleteVariable(temp_vars[j]);
      free(temp_vars);
      fc_printfErrorMessage("Operation failed on step %d", i);
      return rc;
    }
  }

  // stitch steps into the new seq var
  rc = fc_convertVariablesToSeqVariable(numStep, temp_vars, seq,
				     new_seq_var_name, new_seq_var);
  free(temp_vars);
  return rc;
}

/**
 * \ingroup  VariableMath
 * \brief Performs built-in binary operation on a seq variable and a normal var.
 *
 * \description 
 *
 *   Given a string that represents the desired binary operation, return a new
 *   seq variable with the result of that operation between a given
 *   seq variable and a normal var applied at each sequence step.
 *   For example, one can use this to scale a seqvar's values by their mean over
 *   the sequence (e.g., time).
 *
 *    Supported operations are:
 *     - <tt> "+"  =>  var1 + var2 </tt>
 *     - <tt> "-"  =>  var1 - var2 </tt>
 *     - <tt> "*"  =>  var1 * var2 </tt>
 *     - <tt> "/"  =>  var1 / var2 </tt>
 * 
 *    The two seq variable and the normal variable must have the same number of data
 *    points and the same number of components.
 *
 *    The seq new variable will be on the same mesh and sequence and will have the
 *    same association and mathtype as the seq variable. The operation is 
 *    performed per component and data point at each sequence step.
 *
 *    If both original variables have the data type FC_DT_INT, then the 
 *    operation will be an integer operation and the result
 *    will have data type FC_DT_INT. Otherwise, all variable values are
 *    promoted to double and the operation will be a double operation and
 *    the result will have data type FC_DT_DOUBLE.
 *    - 01/24/06  WSK Rewrote innards to mimimize data creation.
 *
 * \modifications 
 *    - 06/07/05  ACG  Created.
 *    - 01/24/06  WSK Rewrote innards to mimimize data creation.
 */
FC_ReturnCode fc_varOperatorSeqVar(
  FC_Variable var1, /**< input var1 */
  char* operation,      /**< input - operation to be performed on variables */
  int numStep, /**< input - num step in input seq var */
  FC_Variable *seq_var2, /**< input seq var */
  char* new_seq_var_name,   /**< input - name for the new variable */ 
  FC_Variable** new_seq_var  /**< output - the new seq variable */
) {   
  FC_ReturnCode rc;
  FC_Sequence seq;
  FC_Variable* temp_vars;
  int i, j;

  // default return
  if (new_seq_var)
    *new_seq_var = NULL;

  // check input
  if ((!fc_isVariableValid(var1)) ||
      (!fc_isSeqVariableValid(numStep, seq_var2)) ||
      (new_seq_var == NULL) || new_seq_var_name == NULL){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  //other checks are done internally in varoperatorvar

  // log message
  fc_printfLogMessage("Performing operation '%s' on variable '%s'"
		      " and seq variable '%s'", operation,
		      _fc_getVarSlot(var1)->header.name,
		      _fc_getVarSlot(seq_var2[0])->header.name);

  // Get info
  rc = fc_getSequenceFromSeqVariable(numStep, seq_var2, &seq);
  if (rc != FC_SUCCESS)
    return rc;

  // Make space to save temp_vars
  temp_vars = malloc(numStep*sizeof(FC_Variable));
  if (!temp_vars) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }

  // Do calculation on each step of the seq variables
  for (i = 0; i < numStep; i++) {
    rc = fc_varOperatorVar(var1, operation, seq_var2[i],"junk_var",
			   &temp_vars[i]); 
    if (rc != FC_SUCCESS) { 
      for (j = 0; j < i-1; j++)
	fc_deleteVariable(temp_vars[j]);
      free(temp_vars);
      fc_printfErrorMessage("Operation failed on step %d", i);
      return rc;
    }
  }

  // stitch steps into the new seq var
  rc = fc_convertVariablesToSeqVariable(numStep, temp_vars, seq,
				     new_seq_var_name, new_seq_var);
  free(temp_vars);
  return rc;
}

/**
 * \ingroup  VariableMath
 * \brief Performs built-in binary operation on a seq variable and constant.
 *
 * \description 
 *
 *   Given a string that represents the desired binary operation, return a new
 *   seq variable with the result of that operation between a given
 *   seq variable and a constant applied at each sequence step.
 *
 *    Supported operations are:
 *     - <tt> "+"  =>  var + const </tt>
 *     - <tt> "-"  =>  var - const </tt>
 *     - <tt> "*"  =>  var * const </tt>
 *     - <tt> "/"  =>  var / const </tt>
 * 
 *    The seq variable and the constant must have the same number of components.
 *
 *    The seq new variable will be on the same mesh and sequence and will have
 *    the same association and mathtype as the seq variable. The operation is
 *    performed per component at each sequence step.
 *
 *    If both original variables have the data type FC_DT_INT, then the
 *    operation will be an integer operation and the result will have data type
 *    FC_DT_INT. Otherwise, all variable values are promoted to double and the
 *    operation will be a double operation and the result will have data type
 *    FC_DT_DOUBLE.
 *
 * \modifications 
 *    - 11/20/2006 WSD Created.
 */
FC_ReturnCode fc_seqVarOperatorConst(
  int numStep, /**< input - num step in input seq var */
  FC_Variable *seq_var1, /**< input seq var1 */
  char* operation,       /**< input - operation to be performed on variables */
  int constNumComp,      /**< input - number of components in constVector */
  void* constVector,     /**< input - constant vector */
  FC_DataType constType, /**< input - data type of constVector */
  char* new_seq_var_name,   /**< input - name for the new variable */ 
  FC_Variable** new_seq_var  /**< output - the new seq variable */
) {   
  FC_ReturnCode rc;
  FC_Sequence seq;
  FC_Variable* temp_vars;
  int i, j;

  // default return
  if (new_seq_var)
    *new_seq_var = NULL;

  // check input
  if ((!fc_isSeqVariableValid(numStep, seq_var1)) || !operation ||
      constNumComp < 1 || !constVector || !fc_isDataTypeValid(constType) ||
      constType == FC_DT_UNKNOWN ||
      (new_seq_var == NULL) || new_seq_var_name == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  //other checks are done internally in varoperatorvar

  // log message
  fc_printfLogMessage("Performing operation '%s' on seq variable '%s'",
                      operation, 
		      _fc_getVarSlot(seq_var1[0])->header.name);

  // Get info
  rc = fc_getSequenceFromSeqVariable(numStep, seq_var1, &seq);
  if (rc != FC_SUCCESS)
    return rc;

  // Make space to save temp_vars
  temp_vars = malloc(numStep*sizeof(FC_Variable));
  if (!temp_vars) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }

  // Do calculation on each step of the seq variables
  for (i = 0; i < numStep; i++) {
    rc = fc_varOperatorConst(seq_var1[i], operation, constNumComp, 
                             constVector, constType,
                             "junk_var", &temp_vars[i]); 
    if (rc != FC_SUCCESS) { 
      for (j = 0; j < i-1; j++)
	fc_deleteVariable(temp_vars[j]);
      free(temp_vars);
      fc_printfErrorMessage("Operation failed on step %d", i);
      return rc;
    }
  }

  // stitch steps into the new seq var
  rc = fc_convertVariablesToSeqVariable(numStep, temp_vars, seq,
				     new_seq_var_name, new_seq_var);
  free(temp_vars);
  return rc;
}

/**
 * \ingroup  VariableMath
 * \brief Performs built-in binary operation on a seq variable and constant.
 *
 * \description 
 *
 *   Given a string that represents the desired binary operation, return a new
 *   seq variable with the result of that operation between a given
 *   seq variable and a constant applied at each sequence step.
 *
 *    Supported operations are:
 *     - <tt> "+"  =>  const + var </tt>
 *     - <tt> "-"  =>  const - var </tt>
 *     - <tt> "*"  =>  const * var </tt>
 *     - <tt> "/"  =>  const / var </tt>
 * 
 *    The seq variable and the constant must have the same number of components.
 *
 *    The seq new variable will be on the same mesh and sequence and will have
 *    the same association and mathtype as the seq variable. The operation is
 *    performed per component at each sequence step.
 *
 *    If both original variables have the data type FC_DT_INT, then the
 *    operation will be an integer operation and the result will have data type
 *    FC_DT_INT. Otherwise, all variable values are promoted to double and the
 *    operation will be a double operation and the result will have data type
 *    FC_DT_DOUBLE.
 *
 * \modifications 
 *    - 11/20/2006 WSD Created.
 */
FC_ReturnCode fc_constOperatorSeqVar(
  int constNumComp,      /**< input - number of components in constVector */
  void* constVector,     /**< input - constant vector */
  FC_DataType constType, /**< input - data type of constVector */
  char* operation,       /**< input - operation to be performed on variables */
  int numStep, /**< input - num step in input seq var */
  FC_Variable *seq_var2, /**< input seq var1 */
  char* new_seq_var_name,   /**< input - name for the new variable */ 
  FC_Variable** new_seq_var  /**< output - the new seq variable */
) {   
  FC_ReturnCode rc;
  FC_Sequence seq;
  FC_Variable* temp_vars;
  int i, j;

  // default return
  if (new_seq_var)
    *new_seq_var = NULL;

  // check input
  if ((!fc_isSeqVariableValid(numStep, seq_var2)) || !operation ||
      constNumComp < 1 || !constVector || !fc_isDataTypeValid(constType) ||
      constType == FC_DT_UNKNOWN ||
      (new_seq_var == NULL) || new_seq_var_name == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  //other checks are done internally in varoperatorvar

  // log message
  fc_printfLogMessage("Performing operation '%s' on seq variable '%s'",
                      operation, 
		      _fc_getVarSlot(seq_var2[0])->header.name);

  // Get info
  rc = fc_getSequenceFromSeqVariable(numStep, seq_var2, &seq);
  if (rc != FC_SUCCESS)
    return rc;

  // Make space to save temp_vars
  temp_vars = malloc(numStep*sizeof(FC_Variable));
  if (!temp_vars) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }

  // Do calculation on each step of the seq variables
  for (i = 0; i < numStep; i++) {
    rc = fc_constOperatorVar(constNumComp, constVector, constType, operation, 
                             seq_var2[i], "junk_var", &temp_vars[i]); 
    if (rc != FC_SUCCESS) { 
      for (j = 0; j < i-1; j++)
	fc_deleteVariable(temp_vars[j]);
      free(temp_vars);
      fc_printfErrorMessage("Operation failed on step %d", i);
      return rc;
    }
  }

  // stitch steps into the new seq var
  rc = fc_convertVariablesToSeqVariable(numStep, temp_vars, seq,
                                        new_seq_var_name, new_seq_var);
  free(temp_vars);
  return rc;
}

//@}

/** \name User supplied functions on variables. */
//-------------------------------
//@{

/**
 * \ingroup  VariableMath
 * \brief Performs unary function on a variable.
 *
 * \description  
 *
 *    Given a unary function that takes a double and returns a double, this
 *    function creates a new variable whose values are the results of the
 *    function applied to the original variable's values.
 *
 *    The new variable will be on the same mesh as the original and will have
 *    the same association and math type as the original. The function is
 *    performed on each component of each data point. The new variable will
 *    always have data type FC_DT_DOUBLE. This function returns with an error
 *    if the data is FC_DT_CHAR.
 *
 * \modifications 
 *    - 02/26/04  R.McCoy  Created.
 */	
FC_ReturnCode fc_varUnaryFunction(
  FC_Variable var,     /**< input - variable  */ 
  double (*unaryFunction)(double),
                        /**< input - function that operates on a double
			 and returns a double */     
  char* new_var_name,   /**< input - name for the new variable */ 
  FC_Variable* new_var  /**< output - the new variable */
){
  FC_ReturnCode rc;
  int i;
  int numData;
  int numComp;
  int numTotalValues;
  FC_AssociationType assoc;
  FC_DataType datatype;
  FC_MathType mathtype;
  void *data;
  double *dblResultPtr;

  // set default return value
  if (new_var)
    *new_var = FC_NULL_VARIABLE;

  // check input
  if (!(fc_isVariableValid(var)) || !new_var_name || !new_var) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  if (!unaryFunction) {
    fc_printfErrorMessage("Function was not provided");
    return FC_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Performing unary user operation on variable '%s'",
            _fc_getVarSlot(var)->header.name);

  // get variable info
  rc = fc_getVariableInfo(var, &numData, &numComp, &assoc, &mathtype, 
			  &datatype);
  if (rc != FC_SUCCESS)
    return rc;

  // check the datatype that can be handled here
  if (datatype != FC_DT_INT && datatype != FC_DT_FLOAT && 
      datatype != FC_DT_DOUBLE) {
    fc_printfErrorMessage("the datatype = %s type cannot be handled",
			  fc_getDataTypeText(datatype) );
    return FC_INPUT_ERROR;
  }  

  numTotalValues = numData * numComp;

  // get data
  rc = fc_getVariableDataPtr(var, &data);
  if (rc != FC_SUCCESS)
    return rc; 

  // calc values
  dblResultPtr = malloc(sizeof(double)*numTotalValues);
  if (dblResultPtr == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;   
  }
  if (datatype == FC_DT_INT) {
    int* cast_data = (int*)data;
    for (i = 0; i < numTotalValues; i++)
      dblResultPtr[i] = (*unaryFunction)(cast_data[i]);
  }
  else if (datatype == FC_DT_FLOAT) {
    float* cast_data = (float*)data;
    for (i = 0; i < numTotalValues; i++)
      dblResultPtr[i] = (*unaryFunction)(cast_data[i]);
  }
  if (datatype == FC_DT_DOUBLE) {
    double* cast_data = (double*)data;
    for (i = 0; i < numTotalValues; i++)
      dblResultPtr[i] = (*unaryFunction)(cast_data[i]);
  }

  // create the new variable
  if (fc_isVariableGlobal(var)) {
    FC_Dataset dataset;
    rc = fc_getDatasetFromVariable(var, &dataset);
    if (rc != FC_SUCCESS)
      return rc;
    rc = fc_createGlobalVariable(dataset, new_var_name, new_var);
    if (rc != FC_SUCCESS)
      return rc; 
  }
  else {
    FC_Mesh mesh;
    rc = fc_getMeshFromVariable(var, &mesh);
    if (rc != FC_SUCCESS)
      return rc;
    rc = fc_createVariable(mesh, new_var_name, new_var);
    if (rc != FC_SUCCESS)
      return rc; 
  }
  rc = fc_setVariableDataPtr(*new_var, numData, numComp, assoc, mathtype, 
			     FC_DT_DOUBLE, dblResultPtr);
  return rc;	 
}

/**
 * \ingroup  VariableMath
 * \brief Performs binary function on two variables.
 *
 * \description
 *
 *    Given a binary function that takes two doubles and returns a double, this
 *    function creates a new variable whose values are the results of the
 *    function applied to the original variables' values.
 *
 *    The two variables must have the same number of data points and the same
 *    number of components.
 *
 *    The new variable will be on the same mesh and will have the same
 *    association and math type as the first variable. The function is applied
 *    per component and data point. The new variable will always have data type
 *    FC_DT_DOUBLE. This function returns with an error if the data is
 *    FC_DT_CHAR.
 *
 * \modifications 
 *    - 02/26/04  R.McCoy  Created.
 */	
FC_ReturnCode fc_varBinaryFunctionVarVar(
  FC_Variable var1,     /**< input - variable 1  */
  FC_Variable var2,     /**< input - variable 2  */
  double (*binaryFunction)(double ddata1, double ddata2),
                        /**< input - function that operates on two doubles
			   and returns a double */
  char* new_var_name,   /**< input - name for the new variable */
  FC_Variable* new_var  /**< output - the new variable */
){
  FC_ReturnCode rc;
  int i;
  int numData1, numData2, numComp1, numComp2, numTotalValues;
  FC_AssociationType assoc1, assoc2;
  FC_DataType datatype1, datatype2;
  FC_MathType mathtype1, mathtype2; 
  void *data1, *data2; 
  double *dblResultPtr;

  // set default return value
  if (new_var)
    *new_var = FC_NULL_VARIABLE;

  // check input
  if (!(fc_isVariableValid(var1)) || !(fc_isVariableValid(var2)) ||
      !new_var_name || !new_var) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));   
    return FC_INPUT_ERROR;
  }
  if (!binaryFunction) {
    fc_printfErrorMessage("Function was not provided");	 
    return FC_ERROR;
  }

  // log message
  fc_printfLogMessage("Performing binary user operation on var1 = '%s' and "
           "var2 = '%s'" , _fc_getVarSlot(var1)->header.name,
           _fc_getVarSlot(var2)->header.name  );  

  // get variable infos
  rc = fc_getVariableInfo(var1, &numData1, &numComp1, &assoc1, 
                            &mathtype1, &datatype1);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_getVariableInfo(var2, &numData2, &numComp2, &assoc2, 
                            &mathtype2, &datatype2);
  if (rc != FC_SUCCESS)
    return rc;

  // check that variables are similar enough
  if (numData1 != numData2 || numComp1 != numComp2) {
    fc_printfErrorMessage("Variables are not similar enough.");
    return FC_INPUT_ERROR;
  }

  // check the data types that can be handled here
  if (datatype1 != FC_DT_INT && datatype1 != FC_DT_FLOAT && 
      datatype1 != FC_DT_DOUBLE) {
    fc_printfErrorMessage("datatype for variable 1 '%s' cannot be handled",
			  fc_getDataTypeText(datatype1));
    return FC_ERROR;
  } 
  if (datatype2 != FC_DT_INT && datatype2 != FC_DT_FLOAT && 
      datatype2 != FC_DT_DOUBLE) {
    fc_printfErrorMessage("datatype2 for variable 2 '%s' cannot be handled",
			  fc_getDataTypeText(datatype2));
    return FC_ERROR;
  }

  numTotalValues = numData1 * numComp1;

  // get data
  rc = fc_getVariableDataPtr(var1, &data1);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_getVariableDataPtr(var2, &data2);
  if (rc != FC_SUCCESS)
    return rc;

  // calc values -- 9 combinations of type
  dblResultPtr = malloc(sizeof(double)*numTotalValues);
  if (dblResultPtr == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;   
  }
  if (datatype1 == FC_DT_INT && datatype2 == FC_DT_INT) {
    int* cast_data1 = (int*)data1;
    int* cast_data2 = (int*)data2;
    for (i = 0; i < numTotalValues; i++)
      dblResultPtr[i] = (*binaryFunction)(cast_data1[i], cast_data2[i]);
  }
  else if (datatype1 == FC_DT_INT && datatype2 == FC_DT_FLOAT) {
    int* cast_data1 = (int*)data1;
    float* cast_data2 = (float*)data2;
    for (i = 0; i < numTotalValues; i++)
      dblResultPtr[i] = (*binaryFunction)(cast_data1[i], cast_data2[i]);
  }
  else if (datatype1 == FC_DT_INT && datatype2 == FC_DT_DOUBLE) {
    int* cast_data1 = (int*)data1;
    double* cast_data2 = (double*)data2;
    for (i = 0; i < numTotalValues; i++)
      dblResultPtr[i] = (*binaryFunction)(cast_data1[i], cast_data2[i]);
  }
  else if (datatype1 == FC_DT_FLOAT && datatype2 == FC_DT_INT) {
    float* cast_data1 = (float*)data1;
    int* cast_data2 = (int*)data2;
    for (i = 0; i < numTotalValues; i++)
      dblResultPtr[i] = (*binaryFunction)(cast_data1[i], cast_data2[i]);
  }
  else if (datatype1 == FC_DT_FLOAT && datatype2 == FC_DT_FLOAT) {
    float* cast_data1 = (float*)data1;
    float* cast_data2 = (float*)data2;
    for (i = 0; i < numTotalValues; i++)
      dblResultPtr[i] = (*binaryFunction)(cast_data1[i], cast_data2[i]);
  }
  else if (datatype1 == FC_DT_FLOAT && datatype2 == FC_DT_DOUBLE) {
    float* cast_data1 = (float*)data1;
    double* cast_data2 = (double*)data2;
    for (i = 0; i < numTotalValues; i++)
      dblResultPtr[i] = (*binaryFunction)(cast_data1[i], cast_data2[i]);
  }
  else if (datatype1 == FC_DT_DOUBLE && datatype2 == FC_DT_INT) {
    double* cast_data1 = (double*)data1;
    int* cast_data2 = (int*)data2;
    for (i = 0; i < numTotalValues; i++)
      dblResultPtr[i] = (*binaryFunction)(cast_data1[i], cast_data2[i]);
  }
  else if (datatype1 == FC_DT_DOUBLE && datatype2 == FC_DT_FLOAT) {
    double* cast_data1 = (double*)data1;
    float* cast_data2 = (float*)data2;
    for (i = 0; i < numTotalValues; i++)
      dblResultPtr[i] = (*binaryFunction)(cast_data1[i], cast_data2[i]);
  }
  else if (datatype1 == FC_DT_DOUBLE && datatype2 == FC_DT_DOUBLE) {
    double* cast_data1 = (double*)data1;
    double* cast_data2 = (double*)data2;
    for (i = 0; i < numTotalValues; i++)
      dblResultPtr[i] = (*binaryFunction)(cast_data1[i], cast_data2[i]);
  }
  else {
    fc_printfErrorMessage("Should never reach here, error in code!\n");
    return FC_ERROR;
  }

  // create the new variable -- put on same mesh as first variable
  if (fc_isVariableGlobal(var1)) {
    FC_Dataset dataset;
    rc = fc_getDatasetFromVariable(var1, &dataset);
    if (rc != FC_SUCCESS)
      return rc;
    rc = fc_createGlobalVariable(dataset, new_var_name, new_var);
    if (rc != FC_SUCCESS)
      return rc; 
  }
  else {
    FC_Mesh mesh;
    rc = fc_getMeshFromVariable(var1, &mesh);
    if (rc != FC_SUCCESS)
      return rc;
    rc = fc_createVariable(mesh, new_var_name, new_var);
    if (rc != FC_SUCCESS)
      return rc; 
  }
  rc = fc_setVariableDataPtr(*new_var, numData1, numComp1, assoc1, mathtype1, 
			     FC_DT_DOUBLE, dblResultPtr);
  return rc;
}

/**
 * \ingroup  VariableMath
 * \brief Performs binary function on a variable and a constant.
 *
 * \description
 *
 *    Given a binary function that takes two doubles and returns a double, this
 *    function creates a new variable whose values are the results of the
 *    function applied with the constant to each data point of the variable.
 *    Use fc_varBinaryFunctionConstVar() to apply the function in reverse
 *    order.
 *
 *    The variable and the constant must have the same number of components.
 *
 *    The new variable will be on the same mesh and will have the same
 *    association and math type as the variable. The function is applied
 *    per component and data point. The new variable will always have data type
 *    FC_DT_DOUBLE. This function returns with an error if the data is
 *    FC_DT_CHAR.
 *
 * \modifications 
 *    - 11/09/04 WSK, created.
 */	
FC_ReturnCode fc_varBinaryFunctionVarConst(
  FC_Variable var,      /**< input - variable */
  int constNumComp,     /**< input - number of components in constant */
  double* constValues,  /**< input - the values of the constant */
  double (*binaryFunction)(double, double),
                        /**< input - function that operates on two doubles
			   and returns a double */
  char* new_var_name,   /**< input - name for the new variable */
  FC_Variable* new_var  /**< output - the new variable */
){
  FC_ReturnCode rc;
  int i, j;
  int numData, numComp;
  FC_AssociationType assoc;
  FC_DataType datatype;
  FC_MathType mathtype; 
  void *data; 
  double *dblResultPtr;

  // set default return value
  if (new_var)
    *new_var = FC_NULL_VARIABLE;

  // check input
  if (!(fc_isVariableValid(var)) || constNumComp < 1 || !constValues ||
      !new_var_name || !new_var) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));   
    return FC_INPUT_ERROR;
  }
  if (!binaryFunction) {
    fc_printfErrorMessage("Function was not provided");	 
    return FC_ERROR;
  }

  // log message
  fc_printfLogMessage("Performing binary user operation on var = '%s' and "
           "a constant" , _fc_getVarSlot(var)->header.name);  

  // get variable info
  rc = fc_getVariableInfo(var, &numData, &numComp, &assoc, 
			  &mathtype, &datatype);
  if (rc != FC_SUCCESS)
    return rc;

  // check that variables and const are similar enough
  if (numComp != constNumComp) {
    fc_printfErrorMessage("Variables and const do not have the same numComp");
    return FC_INPUT_ERROR;
  }

  // check that the datatype can be handled here
  if (datatype != FC_DT_INT && datatype != FC_DT_FLOAT && 
      datatype != FC_DT_DOUBLE) {
    fc_printfErrorMessage("datatype '%s' cannot be handled",
			  fc_getDataTypeText(datatype));
    return FC_ERROR;
  }

  // get data
  rc = fc_getVariableDataPtr(var, &data);
  if (rc != FC_SUCCESS)
    return rc;

  // calc values
  dblResultPtr = malloc(sizeof(double)*numData*numComp);
  if (dblResultPtr == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;   
  }
  if (datatype == FC_DT_INT) {
    int* cast_data = (int*)data;
    for (i = 0; i < numData; i++)
      for (j = 0; j < numComp; j++)
	dblResultPtr[i*numComp+j] = (*binaryFunction)(cast_data[i*numComp+j], 
						    constValues[j]);
  }
  else if (datatype == FC_DT_FLOAT) {
    float* cast_data = (float*)data;
    for (i = 0; i < numData; i++)
      for (j = 0; j < numComp; j++)
	dblResultPtr[i*numComp+j] = (*binaryFunction)(cast_data[i*numComp+j], 
						    constValues[j]);
  }
  else if (datatype == FC_DT_DOUBLE) {
    double* cast_data = (double*)data;
    for (i = 0; i < numData; i++)
      for (j = 0; j < numComp; j++)
	dblResultPtr[i*numComp+j] = (*binaryFunction)(cast_data[i*numComp+j], 
						    constValues[j]);
  }
  else {
    fc_printfErrorMessage("Should never reach here, error in code!\n");
    return FC_ERROR;
  }

  // create the new variable
  if (fc_isVariableGlobal(var)) {
    FC_Dataset dataset;
    rc = fc_getDatasetFromVariable(var, &dataset);
    if (rc != FC_SUCCESS)
      return rc;
    rc = fc_createGlobalVariable(dataset, new_var_name, new_var);
    if (rc != FC_SUCCESS)
      return rc; 
  }
  else {
    FC_Mesh mesh;
    rc = fc_getMeshFromVariable(var, &mesh);
    if (rc != FC_SUCCESS)
      return rc;
    rc = fc_createVariable(mesh, new_var_name, new_var);
    if (rc != FC_SUCCESS)
      return rc; 
  }
  rc = fc_setVariableDataPtr(*new_var, numData, numComp, assoc, mathtype, 
			     FC_DT_DOUBLE, dblResultPtr);
  return rc;
}

/**
 * \ingroup  VariableMath
 * \brief Performs binary function on a constant and a variable.
 *
 * \description
 *
 *    Given a binary function that takes two doubles and returns a double, this
 *    function creates a new variable whose values are the results of the
 *    function applied with the constant to each data point of the variable.
 *    Use fc_varBinaryFunctionConstVar() to apply the function in reverse
 *    order.
 *
 *    The variable and the constant must have the same number of components.
 *
 *    The new variable will be on the same mesh and will have the same
 *    association and math type as the variable. The function is applied
 *    per component and data point. The new variable will always have data type
 *    FC_DT_DOUBLE. This function returns with an error if the data is
 *    FC_DT_CHAR.
 *
 * \modifications 
 *    - 11/09/04 WSK, created.
 */
FC_ReturnCode fc_varBinaryFunctionConstVar(
  int constNumComp,     /**< input - number of components in constant */
  double* constValues,  /**< input - the values of the constant */
  FC_Variable var,      /**< input - variable */
  double (*binaryFunction)(double, double),
                        /**< input - function that operates on two doubles
			   and returns a double */
  char* new_var_name,   /**< input - name for the new variable */
  FC_Variable* new_var  /**< output - the new variable */
){
  FC_ReturnCode rc;
  int i, j;
  int numData, numComp;
  FC_AssociationType assoc;
  FC_DataType datatype;
  FC_MathType mathtype; 
  void *data; 
  double *dblResultPtr;

  // set default return value
  if (new_var)
    *new_var = FC_NULL_VARIABLE;

  // check input
  if (!(fc_isVariableValid(var)) || constNumComp < 1 || !constValues ||
      !new_var_name || !new_var) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));   
    return FC_INPUT_ERROR;
  }
  if (!binaryFunction) {
    fc_printfErrorMessage("Function was not provided");	 
    return FC_ERROR;
  }

  // log message
  fc_printfLogMessage("Performing binary user operation on var = '%s' and "
           "a constant" , _fc_getVarSlot(var)->header.name);  

  // get variable info
  rc = fc_getVariableInfo(var, &numData, &numComp, &assoc, 
			  &mathtype, &datatype);
  if (rc != FC_SUCCESS)
    return rc;

  // check that variables and const are similar enough
  if (numComp != constNumComp) {
    fc_printfErrorMessage("Variables and const do not have the same numComp");
    return FC_INPUT_ERROR;
  }

  // check that the datatype can be handled here
  if (datatype != FC_DT_INT && datatype != FC_DT_FLOAT && 
      datatype != FC_DT_DOUBLE) {
    fc_printfErrorMessage("datatype '%s' cannot be handled",
			  fc_getDataTypeText(datatype));
    return FC_ERROR;
  }

  // get data
  rc = fc_getVariableDataPtr(var, &data);
  if (rc != FC_SUCCESS)
    return rc;

  // calc values
  dblResultPtr = malloc(sizeof(double)*numData*numComp);
  if (dblResultPtr == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;   
  }
  if (datatype == FC_DT_INT) {
    int* cast_data = (int*)data;
    for (i = 0; i < numData; i++)
      for (j = 0; j < numComp; j++)
	dblResultPtr[i*numComp+j] = (*binaryFunction)(constValues[j], 
						    cast_data[i*numComp+j]);
  }
  else if (datatype == FC_DT_FLOAT) {
    float* cast_data = (float*)data;
    for (i = 0; i < numData; i++)
      for (j = 0; j < numComp; j++)
	dblResultPtr[i*numComp+j] = (*binaryFunction)(constValues[j], 
						    cast_data[i*numComp+j]);
  }
  else if (datatype == FC_DT_DOUBLE) {
    double* cast_data = (double*)data;
    for (i = 0; i < numData; i++)
      for (j = 0; j < numComp; j++)
	dblResultPtr[i*numComp+j] = (*binaryFunction)(constValues[j], 
						    cast_data[i*numComp+j]);
  }
  else {
    fc_printfErrorMessage("Should never reach here, error in code!\n");
    return FC_ERROR;
  }

  // create the new variable
  if (fc_isVariableGlobal(var)) {
    FC_Dataset dataset;
    rc = fc_getDatasetFromVariable(var, &dataset);
    if (rc != FC_SUCCESS)
      return rc;
    rc = fc_createGlobalVariable(dataset, new_var_name, new_var);
    if (rc != FC_SUCCESS)
      return rc; 
  }
  else {
    FC_Mesh mesh;
    rc = fc_getMeshFromVariable(var, &mesh);
    if (rc != FC_SUCCESS)
      return rc;
    rc = fc_createVariable(mesh, new_var_name, new_var);
    if (rc != FC_SUCCESS)
      return rc; 
  }
  rc = fc_setVariableDataPtr(*new_var, numData, numComp, assoc, mathtype, 
			     FC_DT_DOUBLE, dblResultPtr);
  return rc;
}

//@}

/** \name User functions on sequence variables. */
//-------------------------------
//@{

/**
 * \ingroup  VariableMath
 * \brief Performs unary function on a seq variable.
 *
 * \description  
 *
 *    Given a unary function that takes a double and returns a double, this
 *    function creates a new seq variable whose values are the results of the
 *    function applied to the original seq variable's values. This is the
 *    seqvar analogy of fc_varUnaryFunction
 *
 *    The new seq variable will be on the same mesh and sequence as the original
 *    and will have the same association and math type as the original. 
 *    The function is performed on each component of each data point at each sequence
 *    step. The new variable will
 *    always have data type FC_DT_DOUBLE. This function returns with an error
 *    if the data is FC_DT_CHAR.
 *
 *   NOTE: this uses createSeqVariable so user will have to deleteSeqVariable 
 *   *and* free it unless deleteseqvar is changed.
 *
 * \modifications 
 *    - 06/07/05  ACG  Created.
 *    - 06/07/05  ACG removing datatype check and leaving that to internal call
 *    - 01/24/06  WSK Rewrote innards to mimimize data creation.
 */	
FC_ReturnCode fc_seqVarUnaryFunction(
  int numStep, /**< input - num step in input seq var */
  FC_Variable *seq_var, /**< input seq var */
  double (*unaryFunction)(double),
                        /**< input - function that operates on a double
			 and returns a double */     
  char* new_seq_var_name,   /**< input - name for the new variable */ 
  FC_Variable** new_seq_var  /**< output - the new seq variable */
) {
  FC_ReturnCode rc;
  FC_Sequence seq;
  FC_Variable* temp_vars;
  int i, j;
  
  // check input
  if (new_seq_var)
    *new_seq_var = NULL;

  // check input
  if ((!fc_isSeqVariableValid(numStep,seq_var)) || (new_seq_var == NULL)
      || new_seq_var_name == NULL){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  //other checks are done internally in varUnaryFunction (incl numerical type),

  // log message
  fc_printfLogMessage("Performing UnaryFunctionSeqVar on seq variable '%s'",
		      _fc_getVarSlot(seq_var[0])->header.name);

  // Get info
  rc = fc_getSequenceFromSeqVariable(numStep,seq_var,&seq);
  if (rc != FC_SUCCESS)
    return rc;

  // Make space to save temp_vars
  temp_vars = malloc(numStep*sizeof(FC_Variable));
  if (!temp_vars) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  
  // Do calculation on each step of the seq variables
  for (i = 0; i < numStep; i++) {
    rc = fc_varUnaryFunction(seq_var[i],unaryFunction, "junk_var", &temp_vars[i]);
    if (rc != FC_SUCCESS) {
      for (j = 0; j < i-1; j++)
	fc_deleteVariable(temp_vars[j]);
      free(temp_vars);
      fc_printfErrorMessage("Operation failed on step %d", i);
      return rc;
    }
  }

  // stitch steps into the new seq var
  rc = fc_convertVariablesToSeqVariable(numStep, temp_vars, seq,
				     new_seq_var_name, new_seq_var);
  free(temp_vars);
  return rc;
}

/**
 * \ingroup  VariableMath
 * \brief Performs binary function on two seq variables.
 *
 * \description
 *
 *    Given a binary function that takes two doubles and returns a double, this
 *    function creates a newseq  variable whose values are the results of the
 *    function applied to the original seq variables' values. This is the
 *    seq var version of varBinaryFunctionVarVar.
 *
 *    The two variables must have the same number of data points and the same
 *    number of components and be on the same sequence.
 *
 *    The new variable will be on the same mesh and sequence and will have the same
 *    association and math type as the first variable. The function is applied
 *    per component and data point. The new variable will always have data type
 *    FC_DT_DOUBLE. This function returns with an error if the data is
 *    FC_DT_CHAR.
 *
 * \modifications 
 *    - 06/07/05  ACG  Created.
 *    - 01/24/06  WSK Rewrote innards to mimimize data creation.
 */	
FC_ReturnCode fc_seqVarBinaryFunctionSeqVarSeqVar(
  int numStep, /**< input - num step in input seq var */
  FC_Variable *seq_var1, /**< input seq var1 */
  FC_Variable *seq_var2, /**< input seq var2 */
  double (*binaryFunction)(double ddata1, double ddata2),
                        /**< input - function that operates on two doubles
			   and returns a double */
  char* new_seq_var_name,   /**< input - name for the new variable */ 
  FC_Variable** new_seq_var  /**< output - the new seq variable */
){
  FC_ReturnCode rc;
  FC_Sequence seq1, seq2;
  FC_Variable* temp_vars;
  int i, j;

  // default return
  if (new_seq_var)
    *new_seq_var = NULL;

  // check input (this indirectly checks that numStep is the same on both seqs
  if ((!fc_isSeqVariableValid(numStep,seq_var1)) ||
      (!fc_isSeqVariableValid(numStep,seq_var2)) ||
      (new_seq_var == NULL) || new_seq_var_name == NULL){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  //other checks are done internally in varoperatorvar

  // log message
  fc_printfLogMessage("Performing SeqVarBinaryFunctionSeqVarSeqVar on"
		      "seq variables '%s' and  '%s'",
		      _fc_getVarSlot(seq_var1[0])->header.name,
		      _fc_getVarSlot(seq_var2[0])->header.name);



  // get info
  rc = fc_getSequenceFromSeqVariable(numStep,seq_var1,&seq1);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_getSequenceFromSeqVariable(numStep,seq_var2,&seq2);
  if (rc != FC_SUCCESS)
    return rc;
  // FIX? remove this?
  if (!(FC_HANDLE_EQUIV(seq1,seq2))){
    fc_printfErrorMessage("%s", "Variables are not on the same sequence");
    return FC_INPUT_ERROR;
  }
 
  // Make space to save temp_vars
  temp_vars = malloc(numStep*sizeof(FC_Variable));
  if (!temp_vars) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }

  // Do calculation on each step of the seq variables
  for (i = 0; i < numStep; i++) {
    rc = fc_varBinaryFunctionVarVar(seq_var1[i], seq_var2[i], binaryFunction,
				    "junk_var", &temp_vars[i]); 
    if (rc != FC_SUCCESS) { 
      for (j = 0; j < i-1; j++)
	fc_deleteVariable(temp_vars[j]);
      free(temp_vars);
      fc_printfErrorMessage("Operation failed on step %d", i);
      return rc;
    }
  }

  // stitch steps into the new seq var
  rc = fc_convertVariablesToSeqVariable(numStep, temp_vars, seq1,
				     new_seq_var_name, new_seq_var);
  free(temp_vars);
  return rc;
}

//@}


/** \name User functions on sequence variables. */
//-------------------------------
/**
 * \ingroup  VariableMath
 * \brief Performs binary function on a seqVar and seqConstant.
 *
 * \description
 *
 *    Given a binary function that takes two doubles and returns a double, this
 *    function creates a new sequence variable whose values are the results of the
 *    function applied to each value and component of the input sequence variable
 *    and the input sequence constant. This is the sequence variable version of 
 *    the varBinaryFunctionVarConst.
 *
 *    The seqVariable and seqConstant must have the same number of components, as
 *    well as the same number of steps.
 *
 *    The new variable will be on the same mesh and sequence and will have the same
 *    association and math type as the sequence variable. The function is applied
 *    per component and data point. The new variable will always have data type
 *    FC_DT_DOUBLE. This function returns with an error if the data is
 *    FC_DT_CHAR. 
 *
 * \modifications 
 *    - 03/01/2007 CDU Created
 */	
FC_ReturnCode fc_seqVarBinaryFunctionSeqVarSeqConst(
  int          numStep,       /**< input  - num step in input seq var */
  FC_Variable *seq_var,       /**< input  - seq var1 */
  int          constNumComp,  /**< input  - number of components in constant*/
  void*        constVectors,  /**< input  - constant vectors, stored as vector for
                                            first step, then vector for next step, etc. */
  FC_DataType  constType,     /**< input  - data type of constVector */

  double (*binaryFunction)(double ddata1, double ddata2),
                              /**< input  - function that operates on two doubles
			                    and returns a double */
  char* new_seq_var_name,     /**< input  - name for the new variable */ 
  FC_Variable** new_seq_var   /**< output - the new seq variable */
){
  FC_ReturnCode rc;
  FC_Sequence seq;
  FC_Variable* temp_vars;
  double *temp_const;
  int i, j;

  // default return
  if (new_seq_var)
    *new_seq_var = NULL;

  // check input (this indirectly checks that numStep is the same on both seqs
  if ((!fc_isSeqVariableValid(numStep,seq_var)) ||
      (constNumComp < 1) || (constVectors==NULL)  || 
      (!fc_isDataTypeValid(constType)) || (constType == FC_DT_UNKNOWN) ||
      (new_seq_var == NULL) || new_seq_var_name == NULL){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  //other checks are done internally in varoperatorvar

  // log message
  fc_printfLogMessage("Performing SeqVarBinaryFunctionSeqVarSeqConst on"
		      "seq variable '%s' and a constant",
		      _fc_getVarSlot(seq_var[0])->header.name);

  // get info
  rc = fc_getSequenceFromSeqVariable(numStep,seq_var,&seq);
  if (rc != FC_SUCCESS)
    return rc;
 
  // Make space to save temp_vars
  temp_vars = malloc(numStep*sizeof(FC_Variable));
  if (!temp_vars) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }

  // Convert the constant to double
  temp_const = malloc(constNumComp * sizeof(double));  
  if(!temp_const){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    free(temp_vars);
    return FC_MEMORY_ERROR;
  }

  // Do calculation on each step of the seq variables
  for (i = 0; i < numStep; i++) {

    //Convert this step of the constVector
    for(j=0;j<constNumComp;j++){
      switch(constType){
      case FC_DT_CHAR:
	temp_const[j] = (double) ((char  *)constVectors)[i*constNumComp+j]; break;
      case  FC_DT_INT: 
	temp_const[j] = (double) ((int   *)constVectors)[i*constNumComp+j]; break;
      case FC_DT_FLOAT:
	temp_const[j] = (double) ((float *)constVectors)[i*constNumComp+j]; break;
      case FC_DT_DOUBLE:
	temp_const[j] = (double) ((double*)constVectors)[i*constNumComp+j]; break;
      case FC_DT_UNKNOWN:
	fc_printfErrorMessage("should never get here");
	return FC_ERROR;
      }
    }

    //Make the actual call
    rc = fc_varBinaryFunctionVarConst(seq_var[i], 
				      constNumComp, temp_const,
				      binaryFunction,
				      "junk_var", &temp_vars[i]); 
    if (rc != FC_SUCCESS) {
      for (j = 0; j < i-1; j++)
	fc_deleteVariable(temp_vars[j]);
      free(temp_const);
      free(temp_vars);
      fc_printfErrorMessage("Operation failed on step %d", i);
      return rc;
    }
  }
    
  free(temp_const);

  // stitch steps into the new seq var
  rc = fc_convertVariablesToSeqVariable(numStep, temp_vars, seq,
					new_seq_var_name, new_seq_var);
  free(temp_vars);
  return rc;
}
//@}

/** \name User functions on sequence variables. */
//-------------------------------
/**
 * \ingroup  VariableMath
 * \brief Performs binary function on a seqConst and a seqVar.
 *
 * \description
 *
 *    Given a binary function that takes two doubles and returns a double, this
 *    function creates a new sequence variable whose values are the results of the
 *    function applied to each value and component of the input sequence constant
 *    and the input sequence variable. This is the sequence variable version of 
 *    the varBinaryFunctionConstVar.
 *
 *    The seqConstant and seqVariable must have the same number of components, as
 *    well as the same number of steps.
 *
 *    The new variable will be on the same mesh and sequence and will have the same
 *    association and math type as the sequence variable. The function is applied
 *    per component and data point. The new variable will always have data type
 *    FC_DT_DOUBLE. This function returns with an error if the data is
 *    FC_DT_CHAR.
 *
 * \modifications 
 *    - 03/01/2007 CDU Created
 */	
FC_ReturnCode fc_seqVarBinaryFunctionSeqConstSeqVar(
  int          numStep,       /**< input - num step in input seq var */
  int          constNumComp,  /**< input - number of components in constant*/
  void*        constVectors,  /**< input - constant vectors, stored as vector for
                                   first step, then vector for next step, etc. */
  FC_DataType  constType,     /**< input - data type of constVector */
  FC_Variable *seq_var,      /**< input seq var1 */

  double (*binaryFunction)(double ddata1, double ddata2),
                              /**< input - function that operates on two doubles
			          and returns a double */
  char* new_seq_var_name,     /**< input - name for the new variable */ 
  FC_Variable** new_seq_var   /**< output - the new seq variable */
){
  FC_ReturnCode rc;
  FC_Sequence seq;
  FC_Variable* temp_vars;
  double *temp_const;
  int i, j;

  // default return
  if (new_seq_var)
    *new_seq_var = NULL;

  // check input (this indirectly checks that numStep is the same on both seqs
  if ((!fc_isSeqVariableValid(numStep,seq_var)) ||
      (constNumComp < 1) || (constVectors==NULL)  || 
      (!fc_isDataTypeValid(constType)) || (constType == FC_DT_UNKNOWN) ||
      (new_seq_var == NULL) || new_seq_var_name == NULL){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  //other checks are done internally in varoperatorvar

  // log message
  fc_printfLogMessage("Performing SeqVarBinaryFunctionSeqVarSeqConst on"
		      "seq variable '%s' and a constant",
		      _fc_getVarSlot(seq_var[0])->header.name);

  // get info
  rc = fc_getSequenceFromSeqVariable(numStep,seq_var,&seq);
  if (rc != FC_SUCCESS)
    return rc;
 
  // Make space to save temp_vars
  temp_vars = malloc(numStep*sizeof(FC_Variable));
  if (!temp_vars) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }

  // Convert the constant to double
  temp_const = malloc(constNumComp * sizeof(double));  
  if(!temp_const){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    free(temp_vars);
    return FC_MEMORY_ERROR;
  }


  // Do calculation on each step of the seq variables
  for (i = 0; i < numStep; i++) {

    //Convert this step of the constVector
    for(j=0;j<constNumComp;j++){
      switch(constType){
      case FC_DT_CHAR:
	temp_const[j] = (double) ((char  *)constVectors)[i*constNumComp+j]; break;
      case  FC_DT_INT: 
	temp_const[j] = (double) ((int   *)constVectors)[i*constNumComp+j]; break;
      case FC_DT_FLOAT:
	temp_const[j] = (double) ((float *)constVectors)[i*constNumComp+j]; break;
      case FC_DT_DOUBLE:
	temp_const[j] = (double) ((double*)constVectors)[i*constNumComp+j]; break;
      case FC_DT_UNKNOWN:
	fc_printfErrorMessage("should never get here");
	return FC_ERROR;
      }
    }

    //Make the actual call
    rc = fc_varBinaryFunctionConstVar(constNumComp, temp_const,
				      seq_var[i], 
				      binaryFunction,
				      "junk_var", &temp_vars[i]); 
    if (rc != FC_SUCCESS) {
      for (j = 0; j < i-1; j++)
	fc_deleteVariable(temp_vars[j]);
      free(temp_const);
      free(temp_vars);
      fc_printfErrorMessage("Operation failed on step %d", i);
      return rc;
    }
  }
    
  free(temp_const);

  // stitch steps into the new seq var
  rc = fc_convertVariablesToSeqVariable(numStep, temp_vars, seq,
					new_seq_var_name, new_seq_var);
  free(temp_vars);
  return rc;
}
//@}



/** \name User functions on sequence variables. */
//-------------------------------
/**
 * \ingroup  VariableMath
 * \brief Performs binary function on a seqVar and a constant.
 *
 * \description
 *
 *    Given a binary function that takes two doubles and returns a double, this
 *    function creates a new sequence variable whose values are the results of the
 *    function applied to each value and component of the sequence variable and
 *    the constant. This is the sequence variable version of the 
 *    varBinaryFunctionVarConst.
 *
 *    The sequence variable and constant must have the same number of components.
 *
 *    The new variable will be on the same mesh and sequence and will have the same
 *    association and math type as the sequence variable. The function is applied
 *    per component and data point. The new variable will always have data type
 *    FC_DT_DOUBLE. This function returns with an error if the data is
 *    FC_DT_CHAR.
 *
 * \modifications 
 *    - 03/01/2007 CDU Created
 */	
FC_ReturnCode fc_seqVarBinaryFunctionSeqVarConst(
  int numStep,              /**< input - num step in input seq var */
  FC_Variable *seq_var,     /**< input - sequence variable */
  int constNumComp,         /**< input - number of components in constant*/
  double* constValues,      /**< input - values of constant */
  double (*binaryFunction)(double ddata1, double ddata2),
                            /**< input - function that operates on two doubles
			         and returns a double */
  char* new_seq_var_name,   /**< input - name for the new variable */ 
  FC_Variable** new_seq_var /**< output - the new seq variable */
){
  FC_ReturnCode rc;
  FC_Sequence seq;
  FC_Variable* temp_vars;
  int i, j;

  // default return
  if (new_seq_var)
    *new_seq_var = NULL;

  // check input (this indirectly checks that numStep is the same on both seqs
  if ((!fc_isSeqVariableValid(numStep,seq_var)) ||
      (constNumComp < 1) || (constValues==NULL)  ||
      (new_seq_var == NULL) || new_seq_var_name == NULL){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  //other checks are done internally in varoperatorvar

  // log message
  fc_printfLogMessage("Performing SeqVarBinaryFunctionSeqVarConst on"
		      "seq variable '%s' and a constant",
		      _fc_getVarSlot(seq_var[0])->header.name);

  // get info
  rc = fc_getSequenceFromSeqVariable(numStep,seq_var,&seq);
  if (rc != FC_SUCCESS)
    return rc;
 
  // Make space to save temp_vars
  temp_vars = malloc(numStep*sizeof(FC_Variable));
  if (!temp_vars) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }

  // Do calculation on each step of the seq variables
  for (i = 0; i < numStep; i++) {
    rc = fc_varBinaryFunctionVarConst(seq_var[i], 
				      constNumComp, constValues,
				      binaryFunction,
				      "junk_var", &temp_vars[i]); 
    if (rc != FC_SUCCESS) {
      for (j = 0; j < i-1; j++)
	fc_deleteVariable(temp_vars[j]);
      free(temp_vars);
      fc_printfErrorMessage("Operation failed on step %d", i);
      return rc;
    }
  }

  // stitch steps into the new seq var
  rc = fc_convertVariablesToSeqVariable(numStep, temp_vars, seq,
				     new_seq_var_name, new_seq_var);
  free(temp_vars);
  return rc;
}
//@}

/** \name User functions on sequence variables. */
//-------------------------------
/**
 * \ingroup  VariableMath
 * \brief Performs binary function on a constant and a seqVar.
 *
 * \description
 *
 *    Given a binary function that takes two doubles and returns a double, this
 *    function creates a new sequence variable whose values are the results of the
 *    function applied to each value and component of the sequence variable and
 *    the constant. This is the sequence variable version of the
 *    varBinaryFunctionConstVar function.
 *
 *    The sequence variable and constant must have the same number of components.
 *
 *    The new variable will be on the same mesh and sequence and will have the same
 *    association and math type as the sequence variable. The function is applied
 *    per component and data point. The new variable will always have data type
 *    FC_DT_DOUBLE. This function returns with an error if the data is
 *    FC_DT_CHAR.
 *
 * \modifications 
 *    - 03/01/2007 CDU Created
 */	
FC_ReturnCode fc_seqVarBinaryFunctionConstSeqVar(
  int constNumComp,          /**< input - number of components in constant*/
  double* constValues,       /**< input - values of constant */
  int numStep,               /**< input - num step in input seq var */
  FC_Variable *seq_var,      /**< input - sequence variable */
  double (*binaryFunction)(double ddata1, double ddata2),
                             /**< input - function that operates on two doubles
			          and returns a double */
  char* new_seq_var_name,    /**< input - name for the new variable */ 
  FC_Variable** new_seq_var  /**< output - the new seq variable */
){
  FC_ReturnCode rc;
  FC_Sequence seq;
  FC_Variable* temp_vars;
  int i, j;

  // default return
  if (new_seq_var)
    *new_seq_var = NULL;

  // check input (this indirectly checks that numStep is the same on both seqs
  if ((!fc_isSeqVariableValid(numStep,seq_var))  ||
      (constNumComp < 1) || (constValues==NULL)  ||
      (new_seq_var == NULL) || new_seq_var_name == NULL){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  //other checks are done internally in varoperatorvar

  // log message
  fc_printfLogMessage("Performing SeqVarBinaryFunctionConstSeqVar on"
		      "seq variable '%s' and a constant",
		      _fc_getVarSlot(seq_var[0])->header.name);

  // get info
  rc = fc_getSequenceFromSeqVariable(numStep,seq_var,&seq);
  if (rc != FC_SUCCESS)
    return rc;
 
  // Make space to save temp_vars
  temp_vars = malloc(numStep*sizeof(FC_Variable));
  if (!temp_vars) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }

  // Do calculation on each step of the seq variables
  for (i = 0; i < numStep; i++) {
    rc = fc_varBinaryFunctionConstVar(constNumComp, constValues,
				      seq_var[i], 
				      binaryFunction,
				      "junk_var", &temp_vars[i]); 
    if (rc != FC_SUCCESS) {
      for (j = 0; j < i-1; j++)
	fc_deleteVariable(temp_vars[j]);
      free(temp_vars);
      fc_printfErrorMessage("Operation failed on step %d", i);
      return rc;
    }
  }

  // stitch steps into the new seq var
  rc = fc_convertVariablesToSeqVariable(numStep, temp_vars, seq,
					new_seq_var_name, new_seq_var);
  free(temp_vars);
  return rc;
}
//@}


/** \name User functions between variables and sequence variables. */
//-------------------------------
//@{

/**
 * \ingroup  VariableMath
 * \brief Performs binary function on a seq variable and a var.
 *
 * \description
 *
 *    Given a binary function that takes two doubles and returns a double, this
 *    function creates a newseq variable whose values are the results of the
 *    function applied to the original seq variable and a normal variable's
 *    values. 
 *
 *    The two variables must have the same number of data points and the same
 *    number of components.
 *
 *    The new variable will be on the same mesh and sequence and
 *    will have the same
 *    association and math type as the first variable. The function is applied
 *    per component and data point. The new variable will always have data type
 *    FC_DT_DOUBLE. This function returns with an error if the data is
 *    FC_DT_CHAR.
 *
 * \modifications 
 *    - 06/08/05  ACG  Created.
 *    - 06/09/05  ACG renamed for seqvarBFvarseqvar to seqvarBFseqvarvar
 *                in order to agre with naming convention that for functions
 *                the first identifier is the return type.
 *    - 01/24/06  WSK Rewrote innards to mimimize data creation.
 */	
FC_ReturnCode fc_seqVarBinaryFunctionSeqVarVar(
  int numStep, /**< input - num step in input seq var */
  FC_Variable *seq_var1, /**< input seq var1 */
  FC_Variable var2, /**< input var2 */
  double (*binaryFunction)(double ddata1, double ddata2),
                        /**< input - function that operates on two doubles
			   and returns a double */
  char* new_seq_var_name,   /**< input - name for the new variable */ 
  FC_Variable** new_seq_var  /**< output - the new seq variable */
){
  FC_ReturnCode rc;
  FC_Sequence seq;
  FC_Variable* temp_vars;
  int i, j;

  // default return
  if (new_seq_var)
    *new_seq_var = NULL;

  // check input
  if ((!fc_isSeqVariableValid(numStep,seq_var1)) ||
      (!fc_isVariableValid(var2)) ||
      (new_seq_var == NULL) || new_seq_var_name == NULL){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  //other checks are done internally in varoperatorvar

  // log message
  fc_printfLogMessage("Performing SeqVarBinaryFunctionSeqVarVar on"
		      "seq variable '%s' and  variable '%s'",
		      _fc_getVarSlot(seq_var1[0])->header.name,
		      _fc_getVarSlot(var2)->header.name);

  //get info 
  rc = fc_getSequenceFromSeqVariable(numStep, seq_var1, &seq);
  if (rc != FC_SUCCESS)
    return rc;

  temp_vars = malloc(numStep*sizeof(FC_Variable));
  if (!temp_vars) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }

  // Do calculation on each step of the seq variables
  for (i = 0; i < numStep; i++) {
    rc = fc_varBinaryFunctionVarVar(seq_var1[i], var2, binaryFunction,
				    "junk_var", &temp_vars[i]); 
    if (rc != FC_SUCCESS) { 
      for (j = 0; j < i-1; j++)
	fc_deleteVariable(temp_vars[j]);
      free(temp_vars);
      fc_printfErrorMessage("Operation failed on step %d", i);
      return rc;
    }
  }

  // stitch steps into the new seq var
  rc = fc_convertVariablesToSeqVariable(numStep, temp_vars, seq,
				     new_seq_var_name, new_seq_var);
  free(temp_vars);
  return rc;
}

/**
 * \ingroup  VariableMath
 * \brief Performs binary function on a seq variable and a var.
 *
 * \description
 *
 *    Given a binary function that takes two doubles and returns a double, this
 *    function creates a newseq variable whose values are the results of the
 *    function applied to the original seq variable and a normal variable's
 *    values. 
 *
 *    The two variables must have the same number of data points and the same
 *    number of components.
 *
 *    The new variable will be on the same mesh and sequence and
 *    will have the same
 *    association and math type as the first variable. The function is applied
 *    per component and data point. The new variable will always have data type
 *    FC_DT_DOUBLE. This function returns with an error if the data is
 *    FC_DT_CHAR.
 *
 * \modifications 
 *    - 06/08/05  ACG  Created.
 *    - 06/09/05  ACG renamed for seqvarBFvarseqvar to seqvarBFseqvarvar
 *                in order to agre with naming convention that for functions
 *                the first identifier is the return type.
 *    - 01/24/06  WSK Rewrote innards to mimimize data creation.
 */	
FC_ReturnCode fc_seqVarBinaryFunctionVarSeqVar(
  FC_Variable var1, /**< input var1 */
  int numStep, /**< input - num step in input seq var */
  FC_Variable *seq_var2, /**< input seq var2 */
  double (*binaryFunction)(double ddata1, double ddata2),
                        /**< input - function that operates on two doubles
			   and returns a double */
  char* new_seq_var_name,   /**< input - name for the new variable */ 
  FC_Variable** new_seq_var  /**< output - the new seq variable */
){
  FC_ReturnCode rc;
  FC_Sequence seq;
  FC_Variable* temp_vars;
  int i, j;

  // default return
  if (new_seq_var)
    *new_seq_var = NULL;

  // check input
  if ((!fc_isVariableValid(var1)) ||
      (!fc_isSeqVariableValid(numStep, seq_var2)) ||
      (new_seq_var == NULL) || new_seq_var_name == NULL){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  //other checks are done internally in varoperatorvar

  // log message
  fc_printfLogMessage("Performing SeqVarBinaryFunctionSeqVarVar on"
		      "seq variable '%s' and  variable '%s'",
		      _fc_getVarSlot(var1)->header.name,
		      _fc_getVarSlot(seq_var2[0])->header.name);

  //get info 
  rc = fc_getSequenceFromSeqVariable(numStep, seq_var2, &seq);
  if (rc != FC_SUCCESS)
    return rc;

  temp_vars = malloc(numStep*sizeof(FC_Variable));
  if (!temp_vars) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }

  // Do calculation on each step of the seq variables
  for (i = 0; i < numStep; i++) {
    rc = fc_varBinaryFunctionVarVar(var1, seq_var2[i], binaryFunction,
				    "junk_var", &temp_vars[i]); 
    if (rc != FC_SUCCESS) { 
      for (j = 0; j < i-1; j++)
	fc_deleteVariable(temp_vars[j]);
      free(temp_vars);
      fc_printfErrorMessage("Operation failed on step %d", i);
      return rc;
    }
  }

  // stitch steps into the new seq var
  rc = fc_convertVariablesToSeqVariable(numStep, temp_vars, seq,
				     new_seq_var_name, new_seq_var);
  free(temp_vars);
  return rc;
}

//@}

/**
 * \ingroup  PrivateVariableMath
 * \brief Helper function to do negation of int array
 *
 * \description
 *
 *    See \ref PrivateVariableMath.
 *
 * \modifications 
 *    - 11/09/04  WSK, created.
 */
FC_ReturnCode _fc_negateInts(
  int num,         /**< input - length of array */
  int* data,       /**< input - input array */
  int* resultData  /**< input/output - array for results */
) {
  int i;
 
  if (!data || !resultData)
    return FC_INPUT_ERROR;

  for (i = 0; i < num; i++)
      resultData[i] = - data[i];

  return FC_SUCCESS;
}

/**
 * \ingroup  PrivateVariableMath
 * \brief Helper function to do negation of double array
 *
 * \description
 *
 *    See \ref PrivateVariableMath.
 *
 * \modifications 
 *    - 11/09/04  WSK, created.
 */
FC_ReturnCode _fc_negateDoubles(
  int num,         /**< input - length of array */
  double* data,       /**< input - input array */
  double* resultData  /**< input/output - array for results */
) {
  int i;
 
  if (!data || !resultData)
    return FC_INPUT_ERROR;

  for (i = 0; i < num; i++)
      resultData[i] = - data[i];

  return FC_SUCCESS;
}

/**
 * \ingroup  PrivateVariableMath
 * \brief Helper function to do addition of int data arrays.
 *
 * \description
 *
 *    See \ref PrivateVariableMath.
 *
 * \modifications 
 *    - 01/29/04  R.McCoy  Created.
 *    - 02/20/04  changed to enable operations with "constant" arrays.
 *    - 11/09/04  WSK, changed to allow operations with constants in either
 *        position, and no longer allocates memory.
 */
FC_ReturnCode _fc_addInts(
  int numComp,    /**< input - number of components in all data arrays */
  int numDataPoint1, /**< input - number of data points in data1 */
  int* data1,        /**< input - data1 array */
  int numDataPoint2, /**< input - number of data points in data2 */
  int* data2,        /**< input - data2 array */
  int* resultData    /**< input/output - data array for results */
) {
  int i, j;
 
  if (!data1 || !data2 || !resultData)
    return FC_INPUT_ERROR;

  if (numDataPoint1 == numDataPoint2) {
    for (i = 0; i < numDataPoint1*numComp; i++)
      resultData[i] = data1[i] + data2[i];
  }
  else if (numDataPoint1 == 1) { // 1st entry is a constant
    for (i = 0; i < numDataPoint2; i++)
      for (j = 0; j < numComp; j++)
	resultData[i*numComp+j] = data1[j] + data2[i*numComp+j];
  }
  else if (numDataPoint2 == 1) { // 2nd entry is a constant
    for (i = 0; i < numDataPoint1; i++)
      for (j = 0; j < numComp; j++)
	resultData[i*numComp+j] = data1[i*numComp+j] + data2[j]; 
  }
  else 
    return FC_ERROR;

  return FC_SUCCESS;
}

/**
 * \ingroup  PrivateVariableMath
 * \brief Helper function to do subtraction of int data arrays.
 *
 * \description
 *
 *    See \ref PrivateVariableMath.
 *
 * \modifications 
 *    - 01/29/04  R.McCoy  Created.
 *    - 02/20/04  changed to enable operations with "constant" arrays.
 *    - 11/09/04  WSK, changed to allow operations with constants in either
 *        position, and no longer allocates memory.
 */
FC_ReturnCode _fc_subtractInts(
  int numComp,    /**< input - number of components in all data arrays */
  int numDataPoint1, /**< input - number of data points in data1 */
  int* data1,        /**< input - data1 array */
  int numDataPoint2, /**< input - number of data points in data2 */
  int* data2,        /**< input - data2 array */
  int* resultData    /**< input/output - data array for results */
) {
  int i, j;

  if (!data1 || !data2 || !resultData)
    return FC_INPUT_ERROR;

  if (numDataPoint1 == numDataPoint2) {
    for (i = 0; i < numDataPoint1*numComp; i++)
      resultData[i] = data1[i] - data2[i];
  }
  else if (numDataPoint1 == 1) { // 1st entry is a constant
    for (i = 0; i < numDataPoint2; i++)
      for (j = 0; j < numComp; j++)
	resultData[i*numComp+j] = data1[j] - data2[i*numComp+j];
  }
  else if (numDataPoint2 == 1) { // 2nd entry is a constant
    for (i = 0; i < numDataPoint1; i++)
      for (j = 0; j < numComp; j++)
	resultData[i*numComp+j] = data1[i*numComp+j] - data2[j]; 
  }
  else 
    return FC_ERROR;

  return FC_SUCCESS;
} 

/**
 * \ingroup  PrivateVariableMath
 * \brief Helper function to do multiplication of int data arrays.
 *
 * \description
 *
 *    See \ref PrivateVariableMath.
 *
 * \modifications 
 *    - 01/29/04  R.McCoy  Created.
 *    - 02/20/04  changed to enable operations with "constant" arrays.
 *    - 11/09/04  WSK, changed to allow operations with constants in either
 *        position, and no longer allocates memory.
 */
FC_ReturnCode _fc_multiplyInts(
  int numComp,    /**< input - number of components in all data arrays */
  int numDataPoint1, /**< input - number of data points in data1 */
  int* data1,        /**< input - data1 array */
  int numDataPoint2, /**< input - number of data points in data2 */
  int* data2,        /**< input - data2 array */
  int* resultData    /**< input/output - data array for results */
) {
  int i, j;
 
  if (!data1 || !data2 || !resultData)
    return FC_INPUT_ERROR;

  if (numDataPoint1 == numDataPoint2) {
    for (i = 0; i < numDataPoint1*numComp; i++)
      resultData[i] = data1[i] * data2[i];
  }
  else if (numDataPoint1 == 1) { // 1st entry is a constant
    for (i = 0; i < numDataPoint2; i++)
      for (j = 0; j < numComp; j++)
	resultData[i*numComp+j] = data1[j] * data2[i*numComp+j];
  }
  else if (numDataPoint2 == 1) { // 2nd entry is a constant
    for (i = 0; i < numDataPoint1; i++)
      for (j = 0; j < numComp; j++)
	resultData[i*numComp+j] = data1[i*numComp+j] * data2[j]; 
  }
  else 
    return FC_ERROR;

  return FC_SUCCESS;
} 

/**
 * \ingroup  PrivateVariableMath
 * \brief Helper function to do division of int data arrays.
 *
 * \description
 *
 *    See \ref PrivateVariableMath.
 *
 * \modifications 
 *    - 01/29/04  R.McCoy  Created.
 *    - 02/20/04  changed to enable operations with "constant" arrays.
 *    - 11/09/04  WSK, changed to allow operations with constants in either
 *        position, and no longer allocates memory.
 */
FC_ReturnCode _fc_divideInts(
  int numComp,    /**< input - number of components in all data arrays */
  int numDataPoint1, /**< input - number of data points in data1 */
  int* data1,        /**< input - data1 array */
  int numDataPoint2, /**< input - number of data points in data2 */
  int* data2,        /**< input - data2 array */
  int* resultData    /**< input/output - data array for results */
) {
  int i, j;
 
  if (!data1 || !data2 || !resultData)
    return FC_INPUT_ERROR;

  if (numDataPoint1 == numDataPoint2) {
    for (i = 0; i < numDataPoint1*numComp; i++)
      resultData[i] = data1[i] / data2[i];
  }
  else if (numDataPoint1 == 1) { // 1st entry is a constant
    for (i = 0; i < numDataPoint2; i++)
      for (j = 0; j < numComp; j++)
	resultData[i*numComp+j] = data1[j] / data2[i*numComp+j];
  }
  else if (numDataPoint2 == 1) { // 2nd entry is a constant
    for (i = 0; i < numDataPoint1; i++)
      for (j = 0; j < numComp; j++)
	resultData[i*numComp+j] = data1[i*numComp+j] / data2[j]; 
  }
  else 
    return FC_ERROR;

  return FC_SUCCESS;
} 

/**
 * \ingroup  PrivateVariableMath
 * \brief Helper function to do addition of double data arrays.
 *
 * \description
 *
 *    See \ref PrivateVariableMath.
 *
 * \modifications 
 *    - 01/29/04  R.McCoy  Created.
 *    - 02/20/04  changed to enable operations with "constant" arrays.
 *    - 11/09/04  WSK, changed to allow operations with constants in either
 *        position, and no longer allocates memory.
 */
FC_ReturnCode _fc_addDoubles(
  int numComp,    /**< input - number of components in all data arrays */
  int numDataPoint1, /**< input - number of data points in data1 */
  double* data1,     /**< input - data1 array */
  int numDataPoint2, /**< input - number of data points in data2 */
  double* data2,     /**< input - data2 array */
  double* resultData /**< input/output - data array for results */
) {
  int i, j;
 
  if (!data1 || !data2 || !resultData)
    return FC_INPUT_ERROR;

  if (numDataPoint1 == numDataPoint2) {
    for (i = 0; i < numDataPoint1*numComp; i++)
      resultData[i] = data1[i] + data2[i];
  }
  else if (numDataPoint1 == 1) { // 1st entry is a constant
    for (i = 0; i < numDataPoint2; i++)
      for (j = 0; j < numComp; j++)
	resultData[i*numComp+j] = data1[j] + data2[i*numComp+j];
   }
  else if (numDataPoint2 == 1) { // 2nd entry is a constant
    for (i = 0; i < numDataPoint1; i++)
      for (j = 0; j < numComp; j++)
	resultData[i*numComp+j] = data1[i*numComp+j] + data2[j]; 
  }
  else 
    return FC_ERROR;

  return FC_SUCCESS;
}

/**
 * \ingroup  PrivateVariableMath
 * \brief Helper function to do subtraction of double data arrays.
 *
 * \description
 *
 *    See \ref PrivateVariableMath.
 *
 * \modifications 
 *    - 01/29/04  R.McCoy  Created.
 *    - 02/20/04  changed to enable operations with "constant" arrays.
 *    - 11/09/04  WSK, changed to allow operations with constants in either
 *        position, and no longer allocates memory.
 */
FC_ReturnCode _fc_subtractDoubles(
  int numComp,    /**< input - number of components in all data arrays */
  int numDataPoint1, /**< input - number of data points in data1 */
  double* data1,     /**< input - data1 array */
  int numDataPoint2, /**< input - number of data points in data2 */
  double* data2,     /**< input - data2 array */
  double* resultData /**< input/output - data array for results */
) {
  int i, j;
 
  if (!data1 || !data2 || !resultData)
    return FC_INPUT_ERROR;

  if (numDataPoint1 == numDataPoint2) {
    for (i = 0; i < numDataPoint1*numComp; i++)
      resultData[i] = data1[i] - data2[i];
  }
  else if (numDataPoint1 == 1) { // 1st entry is a constant
    for (i = 0; i < numDataPoint2; i++)
      for (j = 0; j < numComp; j++)
	resultData[i*numComp+j] = data1[j] - data2[i*numComp+j];
  }
  else if (numDataPoint2 == 1) { // 2nd entry is a constant
    for (i = 0; i < numDataPoint1; i++)
      for (j = 0; j < numComp; j++)
	resultData[i*numComp+j] = data1[i*numComp+j] - data2[j]; 
  }
  else 
    return FC_ERROR;

  return FC_SUCCESS;
}

/**
 * \ingroup  PrivateVariableMath
 * \brief Helper function to do multiplication of double data arrays.
 *
 * \description
 *
 *    See \ref PrivateVariableMath.
 *
 * \modifications 
 *    - 01/29/04  R.McCoy  Created.
 *    - 02/20/04  changed to enable operations with "constant" arrays.
 *    - 11/09/04  WSK, changed to allow operations with constants in either
 *        position, and no longer allocates memory.
 */
FC_ReturnCode _fc_multiplyDoubles(
  int numComp,    /**< input - number of components in all data arrays */
  int numDataPoint1, /**< input - number of data points in data1 */
  double* data1,     /**< input - data1 array */
  int numDataPoint2, /**< input - number of data points in data2 */
  double* data2,     /**< input - data2 array */
  double* resultData /**< input/output - data array for results */
) {
  int i, j;
 
  if (!data1 || !data2 || !resultData)
    return FC_INPUT_ERROR;

  if (numDataPoint1 == numDataPoint2) {
    for (i = 0; i < numDataPoint1*numComp; i++)
      resultData[i] = data1[i] * data2[i];
  }
  else if (numDataPoint1 == 1) { // 1st entry is a constant
    for (i = 0; i < numDataPoint2; i++)
      for (j = 0; j < numComp; j++)
	resultData[i*numComp+j] = data1[j] * data2[i*numComp+j];
  }
  else if (numDataPoint2 == 1) { // 2nd entry is a constant
    for (i = 0; i < numDataPoint1; i++)
      for (j = 0; j < numComp; j++)
	resultData[i*numComp+j] = data1[i*numComp+j] * data2[j]; 
  }
  else
    return FC_ERROR;

  return FC_SUCCESS;
}

/**
 * \ingroup  PrivateVariableMath
 * \brief Helper function to do addition of double data arrays.
 *
 * \description
 *
 *    See \ref PrivateVariableMath.
 *
 * \modifications 
 *    - 01/29/04  R.McCoy  Created.
 *    - 02/20/04  changed to enable operations with "constant" arrays.
 *    - 11/09/04  WSK, changed to allow operations with constants in either
 *        position, and no longer allocates memory.
 */
FC_ReturnCode _fc_divideDoubles(
  int numComp,    /**< input - number of components in all data arrays */
  int numDataPoint1, /**< input - number of data points in data1 */
  double* data1,     /**< input - data1 array */
  int numDataPoint2, /**< input - number of data points in data2 */
  double* data2,     /**< input - data2 array */
  double* resultData /**< input/output - data array for results */
) {
  int i, j;
 
  if (!data1 || !data2 || !resultData)
    return FC_INPUT_ERROR;

  if (numDataPoint1 == numDataPoint2) {
    for (i = 0; i < numDataPoint1*numComp; i++)
      resultData[i] = data1[i] / data2[i];
  }
  else if (numDataPoint1 == 1) { // 1st entry is a constant
    for (i = 0; i < numDataPoint2; i++)
      for (j = 0; j < numComp; j++)
	resultData[i*numComp+j] = data1[j] / data2[i*numComp+j];
  }
  else if (numDataPoint2 == 1) { // 2nd entry is a constant
    for (i = 0; i < numDataPoint1; i++)
      for (j = 0; j < numComp; j++)
	resultData[i*numComp+j] = data1[i*numComp+j] / data2[j]; 
  }
  else 
    return FC_ERROR;

  return FC_SUCCESS;
}

/**
 * \ingroup  PrivateVariableMath
 * \brief Performs built-in binary operation on two data arrays.
 *
 * \description 
 *
 *   This is the core of fc_varOperatorVar(), fc_varOperatorConst(),
 *   and fc_constOperatorVar().
 *
 * \modifications 
 *    - 11/10/03  WSK, created.
 */
FC_ReturnCode _fc_binaryOperatorCore(
  int numComp,        /**< input - number of components in both data arrays */
  int numData1,       /**< input - number of data points in data 1 */
  void* data1,        /**< input - data array 1 */ 
  FC_DataType datatype1, /**< input - type of data 1 */
  char* operation,    /**< input - operation string */
  int numData2,       /**< input - number of data points in data 2 */
  void* data2,        /**< input - data array 2 */ 
  FC_DataType datatype2, /**< input - type of data 2 */
  FC_Variable origVar,   /**< input - reference var */
  char* new_var_name, /**< input - name for the new variable */
  FC_Variable* new_var   /**< output - the new variable */
) {        
  FC_ReturnCode rc;
  int i;
  int opID, new_numData;
  FC_DataType new_datatype;
  FC_MathType mathtype;
  FC_AssociationType assoc;
  void* resultPtr;
  double *dblData1Ptr, *dblData2Ptr;
  // create array of pointers to int binary functions
  FC_ReturnCode (*IntBinaryOper[])(int, int, int*, int, int*, int*) = {       
                          _fc_addInts, 
                          _fc_subtractInts, 
                          _fc_multiplyInts, 
                          _fc_divideInts
                        };
  // corresponding double operation functions: 
  FC_ReturnCode (*DblBinaryOper[])(int, int, double*, int, double*, double*) = {
                          _fc_addDoubles, 
                          _fc_subtractDoubles, 
                          _fc_multiplyDoubles, 
                          _fc_divideDoubles
                        };

  // pick operation
  if (!strcmp(operation, "+")) 
    opID = 0;
  else if (!strcmp(operation, "-"))
    opID = 1;
  else if (!strcmp(operation, "*"))
    opID = 2;
  else if (!strcmp(operation, "/"))
    opID = 3;
  //else if (!strcmp(operation, "%"))
  //opID = 4;
  else {
    fc_printfErrorMessage("Input Error:operation %s not recognized.",
			  operation); 
    return FC_INPUT_ERROR;
  }   

  // no log message

  // check that the data types can be handled here
  if (datatype1 != FC_DT_INT && datatype1 != FC_DT_FLOAT && 
      datatype1 != FC_DT_DOUBLE) {
    fc_printfErrorMessage("data type 1 '%s' cannot be handled.", 
			  fc_getDataTypeText(datatype1));
    return FC_ERROR;
  } 
  if (datatype2 != FC_DT_INT && datatype2 != FC_DT_FLOAT && 
      datatype2 != FC_DT_DOUBLE){
    fc_printfErrorMessage("data type 2 '%s' cannot be handled",
			  fc_getDataTypeText(datatype2));
    return FC_ERROR;
  }

  // get variable info
  rc = fc_getVariableInfo(origVar, &new_numData, &numComp, &assoc, 
			  &mathtype, NULL);
  
  // Do the operation over arrays - results could be int or double    
  if (datatype1 == FC_DT_INT && datatype2 == FC_DT_INT) {
    new_datatype = FC_DT_INT;
    resultPtr = (int*)malloc(sizeof(int)*numComp*new_numData);
    if (resultPtr == NULL) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));  
      return FC_MEMORY_ERROR;
    }
    rc = (*IntBinaryOper[opID])(numComp, numData1, (int*)data1, numData2,
				(int*)data2, (int*)resultPtr);
    if (rc != FC_SUCCESS)
      return rc;
  }
  else {
    new_datatype = FC_DT_DOUBLE;
    // get double representation of data1
    if (datatype1 != FC_DT_DOUBLE) {
      dblData1Ptr = malloc(sizeof(double)*numComp*numData1);
      if (dblData1Ptr == NULL) {
        fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));  
        return FC_MEMORY_ERROR;
      } 
    } 
    if (datatype1 == FC_DT_DOUBLE) 
      dblData1Ptr = (double*) data1;
    else if (datatype1 == FC_DT_INT) { 
      int* cast_data = (int*)data1;  
      for(i = 0; i < numComp*numData1; i++)
	dblData1Ptr[i] = (double)((cast_data)[i]);
    }
    else if (datatype1 == FC_DT_FLOAT) {   // data1 is float
      float* cast_data = (float*)data1;
      for(i = 0; i < numComp*numData1; i++)    
	dblData1Ptr[i] = (double)((cast_data)[i]);
    }
    // get double representation of data2
    if (datatype2 != FC_DT_DOUBLE) {
      dblData2Ptr = malloc(sizeof(double)*numComp*numData2);
      if (dblData2Ptr == NULL) {
        fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
        return FC_MEMORY_ERROR;
      }
    }
    if (datatype2 == FC_DT_DOUBLE)
      dblData2Ptr = (double*)data2;
    else if (datatype2 == FC_DT_INT) {
      int* cast_data = (int*)data2;
      for(i = 0; i < numComp*numData2; i++)    
	dblData2Ptr[i] = (double)((cast_data)[i]);
    }
    else if (datatype2 == FC_DT_FLOAT) {
      float* cast_data = (float*)data2;  
      for(i = 0; i < numComp*numData2; i++)    
	dblData2Ptr[i] = (double)((cast_data)[i]);
    }
    
    // do it
    resultPtr = (int*)malloc(sizeof(double)*numComp*new_numData);
    if (resultPtr == NULL) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));  
      return FC_MEMORY_ERROR;
    }
    rc = (*DblBinaryOper[opID])(numComp, numData1, dblData1Ptr, numData2,
				dblData2Ptr, (double*)resultPtr);
    if (rc != FC_SUCCESS)
      return rc;

    // deallocate memory if it was allocated
    if (datatype1 != FC_DT_DOUBLE)
      free(dblData1Ptr);
    if (datatype2 != FC_DT_DOUBLE)
      free(dblData2Ptr);  
  }
  
  // create the new variable
  if (fc_isVariableGlobal(origVar)) {
    FC_Dataset dataset;
    rc = fc_getDatasetFromVariable(origVar, &dataset);
    if (rc != FC_SUCCESS)
      return rc;
    rc = fc_createGlobalVariable(dataset, new_var_name, new_var);
    if (rc != FC_SUCCESS)
      return rc; 
  }
  else {
    FC_Mesh mesh;
    rc = fc_getMeshFromVariable(origVar, &mesh);
    if (rc != FC_SUCCESS)
      return rc;
    rc = fc_createVariable(mesh, new_var_name, new_var);
    if (rc != FC_SUCCESS)
      return rc; 
  }
  rc = fc_setVariableDataPtr(*new_var, new_numData, numComp, assoc, mathtype, 
			     new_datatype, resultPtr);
  return rc;    
}

/** \name Component aware manipulations. */
//-------------------------------
//@{

/**
 * \ingroup  VariableMath
 * \brief Creates a new variable which is the magnitude of the current variable
 *    
 * \description
 *  
 *    Given a vector variable, this routine will create a new variable whose
 *    values are the magnitudes of the vector variable. Magnitudes are
 *    calculated as the distance from the origin to the point. This routine
 *    only works for vector variables with more than 1 components.
 *
 * \modifications 
 *    - 05/03/04 WSK, Created.
 */
FC_ReturnCode fc_createMagnitudeVariable(
  FC_Variable var,      /**< input - a variable with vector data*/
  char* new_var_name,   /**< input - name for the new variable */
  FC_Variable* new_var  /**< output - the new magnitude of the variable */
) {
        
  FC_ReturnCode rc;
  int i, j;
  int numData, dim;
  FC_AssociationType assoc;
  FC_DataType datatype;
  FC_MathType mathtype;
  void *data;
  double *new_data;

  // set default return value
  if (new_var)
    *new_var = FC_NULL_VARIABLE;

  // check input
  if (!fc_isVariableValid(var) || !new_var_name || !new_var) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));   
    return FC_INPUT_ERROR;
  }

  // get variable info
  rc = fc_getVariableInfo(var, &numData, &dim, &assoc, &mathtype, &datatype);
  if (rc != FC_SUCCESS)
    return rc;
  // check the datatype that can be handled here
  // the datatype can be only FC_DT_INT or FC_DT_FLOAT or FC_DT_DOUBLE
  if (dim < 2 || mathtype != FC_MT_VECTOR || 
      !(datatype == FC_DT_INT || datatype == FC_DT_FLOAT || datatype == FC_DT_DOUBLE)) {
    fc_printfErrorMessage("Input Error: not a vector or bad datatype");
    return FC_ERROR;
  }  

  // log message
  fc_printfLogMessage("Creating new magnitude variable.");

  // get data
  rc = fc_getVariableDataPtr(var, &data);
  if (rc != FC_SUCCESS)
    return rc; 

  // create the new data
  new_data = malloc(sizeof(double)*numData);
  if (new_data == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  for (i = 0; i < numData; i++) {
    new_data[i] = 0;
    for (j = 0; j < dim; j++) {
      switch(datatype) {
      case FC_DT_INT: 
	new_data[i] += ((int*)data)[i*dim+j] * ((int*)data)[i*dim+j];
	break;
      case FC_DT_FLOAT:
	new_data[i] += ((float*)data)[i*dim+j] * ((float*)data)[i*dim+j];
	break;
      case FC_DT_DOUBLE:
	new_data[i] += ((double*)data)[i*dim+j] * ((double*)data)[i*dim+j];
	break;
      default:
	; // nothing
      }
    }
    new_data[i] = sqrt(new_data[i]);
  }

  // create the new variable
  if (fc_isVariableGlobal(var)) {
    FC_Dataset dataset;
    rc = fc_getDatasetFromVariable(var, &dataset);
    if (rc != FC_SUCCESS)
      return rc;
    rc = fc_createGlobalVariable(dataset, new_var_name, new_var);
    if (rc != FC_SUCCESS)
      return rc; 
  }
  else {
    FC_Mesh mesh;
    rc = fc_getMeshFromVariable(var, &mesh);
    if (rc != FC_SUCCESS)
      return rc;
    rc = fc_createVariable(mesh, new_var_name, new_var);
    if (rc != FC_SUCCESS)
      return rc; 
  }
 
  rc = fc_setVariableDataPtr(*new_var, numData, 1, assoc, FC_MT_SCALAR, 
			     FC_DT_DOUBLE, (void*)new_data);
  // do not free new_data 
       
  return rc;    
}       

/**
 * \ingroup  VariableMath
 * \brief Decompose a vector into normal and tangent components respective
 *          to a reference vector.
 *
 * \description
 *  
 *    Most common usage will be to decompose a force on a surface into a
 *    component in the direction of the surface normal, and a component
 *    tangential to the surface.
 *
 *    Note, an infinite number of vectors are tangent to the reference
 *    vector. The calculated tangential component will be in the plane defined
 *    by the input vector and the reference vector. (Perhaps we should also
 *    calculate and return this?)
 *
 *    Note that the dimensions of the input and reference vector do not have to
 *    be equal. The extra dimensions are assumed to have component values of
 *    zero. (In fact, one of the vectors could be a scalar. This would
 *    correspond to a vector aligned with the X axis).
 *
 * \modifications  
 *   - 2003-MAY-15  W Koegler  Created.
 *   - 2003-SEP-05  W Koegler  Moved core of decompose algorithm to its own
 *   function _fc_decomposeVector to better share w/ others functions.
 *   - 11/01/04 WSK Moved from geom to varmath
 */
FC_ReturnCode fc_createNormalTangentVariables(
  FC_Variable inputVector, /**< input - vector variable to decompose */
  int dim_ref,    /**< input - the number of components in reference vector */
  double* referenceVector, /**< input - the components of the reference vector */
  FC_Variable* normalComponent, /**< output - new variable with normal components */
  FC_Variable* tangentComponent /**< output - new variable with tangent components */
) {
  FC_ReturnCode rc;
  int i, j;
  int numData; // number of data points
  int dim_in; // number of components per point
  FC_AssociationType assoc;
  FC_DataType datatype;
  FC_MathType mathtype;
  double magnitude;
  double *normRefVect = NULL;
  double *normals = NULL, *tangents = NULL;
  void *datap = NULL;
  char *vname = NULL;
  char *tempstring1 = NULL, *tempstring2 = NULL;
  
  // set default return values
  if (normalComponent)
    *normalComponent = FC_NULL_VARIABLE;
  if (tangentComponent)
    *tangentComponent = FC_NULL_VARIABLE;
  
  // check input
  if (!fc_isVariableValid(inputVector) || dim_ref < 1 || !referenceVector ||
      !normalComponent || !tangentComponent)  {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Decomposing a vector into normal and tangent "
                      "components.");   
  
  // get info about inputVector
  rc = fc_getVariableInfo(inputVector, &numData, &dim_in, &assoc,  
                          &mathtype, &datatype);
  if (rc != FC_SUCCESS)
    return rc;
  if (mathtype != FC_MT_SCALAR && mathtype != FC_MT_VECTOR)
    return FC_ERROR; // can't handle tensor
  
  // get input vector data
  rc = fc_getVariableDataPtr(inputVector, &datap);
  if (rc != FC_SUCCESS)
    return rc;
  
  // calculate magnitude of reference vector
  magnitude = 0;
  for (i = 0; i < dim_ref; i++)
    magnitude += referenceVector[i]*referenceVector[i];
  magnitude = sqrt(magnitude);
  // FIX!
  if (magnitude <= 0)
    return FC_ERROR; 
  
  // make normalized reference vector
  normRefVect = malloc(sizeof(double)*dim_ref);
  if (normRefVect == NULL ) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  } 
  if (magnitude > 0) {
    for (i = 0; i < dim_ref; i++)
      normRefVect[i] = referenceVector[i]/magnitude;
  }
  
  // make space for temporary normal and tangent component arrays
  normals = malloc(sizeof(double)*numData);
  tangents = malloc(sizeof(double)*numData);
  if (normals == NULL || tangents == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  } 
  
  switch (datatype) {
  case FC_DT_INT: 
    {
      int *castp = (int*)datap;
      double temp[dim_in];
      for (i = 0; i < numData; i++) {
        for (j = 0; j < dim_in; j++)
          temp[j] = castp[i*dim_in + j];
        _fc_decomposeVector(dim_in, temp, dim_ref, normRefVect, &normals[i],
                            &tangents[i]);
      }
    } 
    break;
  case FC_DT_FLOAT: 
    {
      float *castp = (float*)datap;
      double temp[dim_in];
      for (i = 0; i < numData; i++) {
        for (j = 0; j < dim_in; j++)
          temp[j] = castp[i*dim_in + j];
        _fc_decomposeVector(dim_in, temp, dim_ref, normRefVect, &normals[i],
                            &tangents[i]);
      }
    } 
    break;
  case FC_DT_DOUBLE:
    {
      double *castp = (double*)datap;
      for (i = 0; i < numData; i++) 
        _fc_decomposeVector(dim_in, &castp[i*dim_in], dim_ref,
                            normRefVect, &normals[i], &tangents[i]);
    } 
    break;
  default:
    return FC_ERROR;
  }
  free(normRefVect);
  
  // add new variables to current mesh/dataset
  rc = fc_getVariableName(inputVector, &vname);
  if (rc != FC_SUCCESS)
    return rc;
  tempstring1 = malloc(sizeof(char)*(strlen(vname) + 55));
  tempstring2 = malloc(sizeof(char)*(strlen(vname) + 55));
  if (tempstring1 == NULL || tempstring2 == NULL) {
    free(tempstring1);
    free(tempstring2);
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  } 
  sprintf(tempstring1, "%s Normal Component decomposed on a reference vector",
          vname);
  sprintf(tempstring2, "%s Tangent Component decomposed on a reference vector",
          vname);
  free(vname);
  if (fc_isVariableGlobal(inputVector)) {
    FC_Dataset dataset;
    rc = fc_getDatasetFromVariable(inputVector, &dataset);
    if (rc != FC_SUCCESS)
      return rc;
    rc = fc_createGlobalVariable(dataset, tempstring1, normalComponent);
    if (rc != FC_SUCCESS)
      return rc;
    rc = fc_createGlobalVariable(dataset, tempstring2, tangentComponent);
    if (rc != FC_SUCCESS)
      return rc;
  }
  else {
    FC_Mesh mesh;
    rc = fc_getMeshFromVariable(inputVector, &mesh);
    if (rc != FC_SUCCESS)
      return rc;
    rc = fc_createVariable(mesh, tempstring1, normalComponent);
    if (rc != FC_SUCCESS)
      return rc;
    rc = fc_createVariable(mesh, tempstring2, tangentComponent);
    if (rc != FC_SUCCESS)
      return rc;
  }
  free(tempstring1);
  free(tempstring2);
  rc = fc_setVariableDataPtr(*normalComponent, numData, 1, assoc, FC_MT_SCALAR, 
                             FC_DT_DOUBLE, normals);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_setVariableDataPtr(*tangentComponent, numData, 1, assoc, FC_MT_SCALAR, 
                             FC_DT_DOUBLE, tangents);
  if (rc != FC_SUCCESS)
    return rc;
  
  return FC_SUCCESS;
}

/**
 * \ingroup  VariableMath
 * \brief Decompose a vector into normal and tangent components respective
 *          to a reference vector.
 *
 * \description
 *  
 *    See fc_createNormalTangentVariables. This function is similar except
 *    that it takes a field of reference vectors instead of a
 *    single reference vector.
 *
 * \modifications  
 *   - 2003-MAY-15  W Koegler  Created.
 *   - 2003-SEP-05  W Koegler  Moved core of decompose algorithm to its own
 *   function _fc_decomposeVector to better share w/ others functions.
 *   - 11/01/04 WSK Moved from geom to varmath
 */
FC_ReturnCode fc_createNormalTangentVariables2(
  FC_Variable inputVector,      /**< input - vector variable to decompose */
  FC_Variable referenceVector,  /**< input - variable of reference vectors */
  FC_Variable* normalComponent, /**< output - new variable with normal components */
  FC_Variable* tangentComponent /**< output - new variable with tangent components */
) {
  FC_ReturnCode rc;
  int i, j;
  int dim_in;  // number of components in input vector
  int dim_ref; // number of components in reference vector
  int numData_in, numData_ref;
  FC_AssociationType assoc_in, assoc_ref;
  FC_DataType datatype_in, datatype_ref;
  //FC_DataType datatype_min;
  FC_MathType mathtype_in, mathtype_ref;
  double *normRefVect = NULL;
  double *normals = NULL, *tangents = NULL;
  void *data_in = NULL, *data_ref;
  char *name_in = NULL, *name_ref = NULL;
  char *tempstring1 = NULL, *tempstring2 = NULL;
  
  // set default return values
  if (normalComponent)
    *normalComponent = FC_NULL_VARIABLE;
  if (tangentComponent)
    *tangentComponent = FC_NULL_VARIABLE;
  
  // check input
  if (!fc_isVariableValid(inputVector) || 
      !fc_isVariableValid(referenceVector) || !normalComponent || 
      !tangentComponent) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Decomposing vector into a normal and tangent "
                      "components with respect to a reference vector.");

  // get info about inputVector
  rc = fc_getVariableInfo(inputVector, &numData_in, &dim_in, &assoc_in, 
                     &mathtype_in, &datatype_in);
  if (rc != FC_SUCCESS)
    return rc;
  if (mathtype_in != FC_MT_SCALAR && mathtype_in != FC_MT_VECTOR) 
    return FC_ERROR; // can't handle tensors
  
  // get info about reference Vector
  rc = fc_getVariableInfo(referenceVector, &numData_ref, &dim_ref, &assoc_ref, 
                     &mathtype_ref, &datatype_ref);
  if (rc != FC_SUCCESS)
    return rc;
  if (numData_in != numData_ref)
    return FC_ERROR; // Must have same number of data points
  if (mathtype_ref != FC_MT_SCALAR && mathtype_ref != FC_MT_VECTOR) 
    return FC_ERROR; // can't handle tensors
  
  // get vector data
  rc = fc_getVariableDataPtr(inputVector, &data_in);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_getVariableDataPtr(referenceVector, &data_ref);
  if (rc != FC_SUCCESS)
    return rc;
  
  // normalize the reference vector
  // also converts to double
  normRefVect = malloc(sizeof(double)*numData_ref*dim_ref);
  if (normRefVect == NULL ) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  } 
  switch (datatype_ref) {
  case FC_DT_INT:
    {
      int* castp = (int*)data_ref;
      for (i = 0; i < numData_in; i++) {
        for (j = 0; j < dim_ref; j++)
          normRefVect[i*dim_ref + j] = castp[i*dim_ref + j];
      }
    }
    break;
  case FC_DT_FLOAT:
    {
      float* castp = (float*)data_ref;
      for (i = 0; i < numData_in; i++) {
        for (j = 0; j < dim_ref; j++)
          normRefVect[i*dim_ref + j] = castp[i*dim_ref + j];
      }
    }
    break;
  case FC_DT_DOUBLE:
    {
      double* castp = (double*)data_ref;
      for (i = 0; i < numData_in; i++) {
        for (j = 0; j < dim_ref; j++)
          normRefVect[i*dim_ref + j] = castp[i*dim_ref + j];
      }
    }
    break;
  default:
    return FC_ERROR;
  }
  
  // make space for temporary normal and tangent component arrays
  normals = malloc(sizeof(double)*numData_in);
  tangents = malloc(sizeof(double)*numData_in);
  if (normals == NULL || tangents == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  } 
  
  switch (datatype_in) {
  case FC_DT_INT: 
    {
      int *castp = (int*)data_in;
      double temp[dim_in];
      for (i = 0; i < numData_in; i++) {
        for (j = 0; j < dim_in; j++)
          temp[j] = castp[i*dim_in + j];
        _fc_decomposeVector(dim_in, temp, dim_ref, &normRefVect[i*dim_ref], 
                            &normals[i], &tangents[i]);
      }
    } 
    break;
  case FC_DT_FLOAT: 
    {
      float *castp = (float*)data_in;
      double temp[dim_in];
      for (i = 0; i < numData_in; i++) {
        for (j = 0; j < dim_in; j++)
          temp[j] = castp[i*dim_in + j];
        _fc_decomposeVector(dim_in, temp, dim_ref, &normRefVect[i*dim_ref], 
                            &normals[i], &tangents[i]);
      }
    } 
    break;
  case FC_DT_DOUBLE:
    {
      double *castp = (double*)data_in;
      for (i = 0; i < numData_in; i++) 
        _fc_decomposeVector(dim_in, &castp[i*dim_in], dim_ref,
                            &normRefVect[i*dim_ref], &normals[i], &tangents[i]);
    } 
    break;
  default:
    return FC_ERROR;
  }
  free(normRefVect);
  
  // add new variables to current mesh
  // FIX check for name collisions, better if done in create_variable
  rc = fc_getVariableName(inputVector, &name_in);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_getVariableName(referenceVector, &name_ref);
  if (rc != FC_SUCCESS)
    return rc;
  tempstring1 = malloc(sizeof(char)*(strlen(name_in) + strlen(name_ref) +40));
  tempstring2 = malloc(sizeof(char)*(strlen(name_in) + strlen(name_ref) +40));
  if (tempstring1 == NULL || tempstring2 == NULL) {
    free(tempstring1);
    free(tempstring2);
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  sprintf(tempstring1, "%s Normal Component with respect to %s", name_in, name_ref);
  sprintf(tempstring2, "%s Tangent Component with respect to %s", name_in, name_ref);
  free(name_in);
  free(name_ref);
  if (fc_isVariableGlobal(inputVector)) {
    FC_Dataset dataset;
    rc = fc_getDatasetFromVariable(inputVector, &dataset);
    if (rc != FC_SUCCESS)
      return rc;
    rc = fc_createGlobalVariable(dataset, tempstring1, normalComponent);
    if (rc != FC_SUCCESS)
      return rc;
    rc = fc_createGlobalVariable(dataset, tempstring2, tangentComponent);
    if (rc != FC_SUCCESS)
      return rc;
  }
  else {
    FC_Mesh mesh;
    rc = fc_getMeshFromVariable(inputVector, &mesh);
    if (rc != FC_SUCCESS)
      return rc;
    rc = fc_createVariable(mesh, tempstring1, normalComponent);
    if (rc != FC_SUCCESS)
      return rc;
    rc = fc_createVariable(mesh, tempstring2, tangentComponent);
    if (rc != FC_SUCCESS)
      return rc;
  }
  free(tempstring1);
  free(tempstring2);
  rc = fc_setVariableDataPtr(*normalComponent, numData_in, 1, assoc_in, 
                        FC_MT_SCALAR, FC_DT_DOUBLE, normals);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_setVariableDataPtr(*tangentComponent, numData_in, 1, assoc_in, 
                        FC_MT_SCALAR, FC_DT_DOUBLE, tangents);
  if (rc != FC_SUCCESS)
    return rc;
  
  return FC_SUCCESS;
}

/**
 * \ingroup VariableMath
 * \brief Compute Von Mises stress and pressure from the stress matrix.
 *
 * \description
 *
 *    Given the six variables that make up the symmetric stress tensor,
 *    this function computes the Von Mises stress and the pressure.
 *    The expected order of the stress tensor components (passed in as
 *    an array) are { xx, yy, zz, xy, yz, zx };
 *
 *    The calculations are based on equations from "Intermediate Mechanics
 *    of Materials", J. R. Barber, McGraw-Hill, 2001, pgs 33-48. Pressure
 *    is defined such that it is negative for compression and positive
 *    for tension.
 *
 * \todo Is this the best home for this function?
 *
 * \todo this doesn't handle global vars yet
 *
 * \modifications
 *    - 8/08/2006 WSD Created.
 */
FC_ReturnCode fc_createStressPressureVariables(
  FC_Variable* sigmaVars, /**< input - array of 6 stress tensor components
			   in order { xx, yy, zz, xy, yz, zx } */
  FC_Variable* stressVar,   /**< output - the von misses stress */
  FC_Variable* pressureVar  /**< output - the pressure */
) 
{
  FC_ReturnCode rc;
  int i;
  // sigmas: 0 - xx, 1 - yy, 2 - zz, 3 - xy , 4 - yz, 5 - zx
  double* sigmas[6];
  double* stress;
  double* pressure;
  int numDataPoint;
  FC_Mesh mesh;
  FC_AssociationType assoc;

  // default return
  if (stressVar)
    *stressVar = FC_NULL_VARIABLE;
  if (pressureVar)
    *pressureVar = FC_NULL_VARIABLE;

  // check input 
  if (!sigmaVars || !stressVar || !pressureVar) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  for (i = 0; i < 6; i++) {
    int numComp;
    FC_Mesh temp_mesh;
    if (!fc_isVariableValid(sigmaVars[i])) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
      return FC_INPUT_ERROR;
    }
    fc_getVariableNumComponent(sigmaVars[i], &numComp);
    if (numComp != 1) {
      fc_printfErrorMessage("input stress components must be scalars");
      return FC_INPUT_ERROR;
    }
    fc_getMeshFromVariable(sigmaVars[i], &temp_mesh);
    if (i == 0)
      mesh = temp_mesh;
    else if (!FC_HANDLE_EQUIV(mesh, temp_mesh)) {
      fc_printfErrorMessage("all stress components must be on same mesh");
      return FC_INPUT_ERROR;
    }
  }

  // log messaeg
  fc_printfLogMessage("Creating Von Misses stress and pressure variables");

  // get the starting data
  rc = fc_getMeshFromVariable(sigmaVars[0], &mesh);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_getVariableAssociationType(sigmaVars[0], &assoc);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_getVariableNumDataPoint(sigmaVars[0], &numDataPoint);
  if (rc != FC_SUCCESS)
    return rc;
  for (i = 0; i < 6; i++) {
    rc = fc_getVariableDataAsDataType(sigmaVars[i], FC_DT_DOUBLE, 
				      (void*)&sigmas[i]);
    if (rc != FC_SUCCESS)
      return rc;
  }

  // make room for new data
  stress = (double*)malloc(numDataPoint*sizeof(double));
  pressure = (double*)malloc(numDataPoint*sizeof(double));
  if (!stress || !pressure) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }

  // calculate data
  for (i = 0; i < numDataPoint; i++) {
    double I1, I2, I3;
    double phi;
    
    I1 = sigmas[0][i] + sigmas[1][i] + sigmas[2][i];
    I2 =   sigmas[0][i]*sigmas[1][i]
         + sigmas[1][i]*sigmas[2][i]
         + sigmas[2][i]*sigmas[0][i] 
         - sigmas[3][i]*sigmas[3][i] 
         - sigmas[4][i]*sigmas[4][i]
         - sigmas[5][i]*sigmas[5][i];
    I3 =   sigmas[0][i]*sigmas[1][i]*sigmas[2][i]
         - sigmas[0][i]*sigmas[4][i]*sigmas[4][i]
         - sigmas[1][i]*sigmas[5][i]*sigmas[5][i]
         - sigmas[2][i]*sigmas[3][i]*sigmas[3][i]
         + 2*sigmas[3][i]*sigmas[4][i]*sigmas[5][i];

    stress[i] = sqrt(I1*I1 - 3*I2);

    if (FC_DBL_EQUIV(stress[i], 0)) {
      pressure[i] = 1/3. * I1;
    }
    else {
      phi = 1/3. * acos((2*I1*I1*I1 - 9*I1*I2 + 27*I3)/2./
			(stress[i]*stress[i]*stress[i]));      
      pressure[i] = 1/3. * (I1 + 2.*stress[i]*(cos(phi) +
					       cos(phi + 2./3.*FC_PI) + 
					       cos(phi + 4./3.*FC_PI)));
    }
  }
  for (i = 0; i < 6; i++)
    free(sigmas[i]);

  // create new vars
  rc = fc_createVariable(mesh, "von mises stress", stressVar);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_setVariableDataPtr(*stressVar, numDataPoint, 1, assoc,
			     FC_MT_SCALAR, FC_DT_DOUBLE, stress);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_createVariable(mesh, "pressure", pressureVar);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_setVariableDataPtr(*pressureVar, numDataPoint, 1, assoc,
			     FC_MT_SCALAR, FC_DT_DOUBLE, pressure);
  if (rc != FC_SUCCESS)
    return rc;

  return FC_SUCCESS;
}

//@}

/**
 * \ingroup  PrivateVariableMath
 * \brief Decompose a vector into normal and tangent components respective
 *          to a reference vector.
 *
 * \description
 *  
 *    Most common usage will be to decompose a force on a surface into a
 *    component in the direction of the surface normal, and a component
 *    tangential to the surface.
 *
 *    Note, an infinite number of vectors are tangent to the reference
 *    vector. The calculated tangential component will be in the plane defined
 *    by the input vector and the reference vector. (Perhaps we should also
 *    calculate and return this?)
 *
 *    Note, the dimensions of the input and reference vector do not have to be
 *    equal. The extra dimensions are assumed to have component values of
 *    zero. (In fact, one of the vectors could be a scalar. This would
 *    correspond to a vector aligned with the X axis).
 *
 * \modifications  
 *    - 2003-SEP-05  W Koegler  Created by moving code from fc_decomposeVector
 *    - 11/01/04 WSK Moved from geom to varmath
 */
void _fc_decomposeVector(
  int dim_in,    /**< input - the number of components in the input vector */
  double* inputVector, /**< input - the vector to be decomposed */
  int dim_ref,   /**< input - the number of components in the reference vector */
  double* referenceVector, /**< input - the components of the reference vector */
  double* normalComponent, /**< output - the normal component */
  double* tangentComponent /**< output - the tangent component */
) {
  int i;
  int dim_min;
  double inMagnitude_sq;            // squared magnitude of input vector (^2)
  double refMagnitude;              // magnitude of reference vector
  double normalizedRefVec[dim_ref]; // normalized reference vector
  double* refVector_p;              // pointer to the reference vector
  double normal, tangent;
  
  // default output values  
  *normalComponent = 0.;
  *tangentComponent = 0.;
  
  // log message
  fc_printfLogMessage("Decomposing vector.");  
  
  // calculate magnitudes of the vectors & return if they are zero
  inMagnitude_sq = 0;
  for (i = 0; i < dim_in; i++)
    inMagnitude_sq += inputVector[i]*inputVector[i];
  if (sqrt(inMagnitude_sq) < DBL_MIN)
    return;
  refMagnitude = 0;
  for (i = 0; i < dim_ref; i++)
    refMagnitude += referenceVector[i]*referenceVector[i];
  refMagnitude = sqrt(refMagnitude);
  if (refMagnitude < DBL_MIN)
    return;
  
  // if necessary, normalize the reference vector
  if (!FC_DBL_EQUIV(refMagnitude, 1.)) {
    for (i = 0; i < dim_ref; i++)
      normalizedRefVec[i] = referenceVector[i]/refMagnitude;
    refVector_p = normalizedRefVec;
  }
  else
    refVector_p = referenceVector;
  
  // calculate Normal Component, use smaller of dims
  // Fn = F dot N = SUM(fi*ni)   
  dim_min = dim_in < dim_ref ? dim_in : dim_ref;
  normal = 0;
  for (i = 0; i < dim_min; i++)
    normal += inputVector[i]*refVector_p[i];
  
  // calculate tangent component, use
  // Ft = SQRT(|F|^2 - Fn^2) = SQRT( SUM(fi^2) - Fn^2 )
  tangent = inMagnitude_sq - normal*normal;
  if (FC_VALUE_EQUIV(inMagnitude_sq, normal*normal, 4*DBL_EPSILON, DBL_MIN))
    tangent = 0.;
  else
    tangent = sqrt(tangent);
  
  *normalComponent = normal;
  *tangentComponent = tangent;
}

