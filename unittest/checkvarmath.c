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
 * \file checkvarmath.c
 * \brief Unit tests for \ref VariableMath module.
 *
 * $Source: /home/Repositories/fcdmf/fclib/unittest/checkvarmath.c,v $
 * $Revision: 1.74 $ 
 * $Date: 2006/11/21 18:40:09 $
 *
 * \description
 *
 *     To test the functions that allow mathematical operations on
 *     fclib variables. 
 *
 * \todo 
 *     - test fc_exampleUnaryUserFn() &  fc_exampleBinaryUserFn().
 *       They are now being tested by being used in fc_varUnaryUserFunction()
 *       tests and fc_varBinaryUserFunction() tests.
 *
 * \modifications
 *    - 03/05/04 RM, created.
 *    - 05/04/04 WSK, added test for components split/merge.
 *    - 06/03/05 ACG adding analogous things for seqvars
 *    - 06/08/05 ACG revised checks for seqvar fctns. using trimesh
 *      since the data is smaller.these only have to check that the
 *      errors are being passed up thru the wrapper properly and that
 *      things are being called properly and placed int eh return seq val
 *      in the proper order, so im not checking from the start data, but
 *      rather by comparing it on a step by step basis with the non-seqvar
 *      versions
 */

#include <stdlib.h>
#include <string.h>
#include "check.h"
#include "fc.h"
#include "seriesP.h"
#include "varmathP.h"
#include "checkall.h"
    
// *********************************************
// ***** fc_VariableAddVariables tests
// *********************************************

// **** global variables
static FC_Dataset glb_dataset;
static FC_Sequence glb_sequence;
static FC_Mesh glb_mesh;
static int glb_numDataType = 4;
static FC_DataType glb_dataTypes[4] = { FC_DT_CHAR, FC_DT_INT,
					FC_DT_FLOAT, FC_DT_DOUBLE };
static int glb_numStep = 3;
static FC_Variable glb_dataTypeVars[4], glb_vectorVar, glb_elemVar;
static FC_Variable *glb_dataTypeSeqVars[4], *glb_vectorSeqVar, *glb_elemSeqVar;
static FC_Variable glb_glbVar, *glb_glbSeqVar;

// **** helper funtions

// operations
static int negInt(int x) { return -x; }
static double negDbl(double x) { return -x; }
static int addInts(int x, int y) { return x + y; }
static int subInts(int x, int y) { return x - y; }
static int multInts(int x, int y) { return x * y; }
static int divInts(int x, int y) { return x / y; }
static double addDbls(double x, double y) { return x + y; }
static double subDbls(double x, double y) { return x - y; }
static double multDbls(double x, double y) { return x * y; }
static double divDbls(double x, double y) { return x / y; }
// functions (couldn't think of good binary math functions)
static double subtraction(double x, double y) { return x - y; }
static double magnitude(double x, double y) { return sqrt(x*x + y*y); }

// **** test fixtures

static void varmath_setup(void) {
  FC_ReturnCode rc;
  if (isForking) {
    fc_setLibraryVerbosity(fc_messages);
    rc = fc_initLibrary();
    fail_unless(rc == FC_SUCCESS, "test aborted: failed to init library");
  }
}

static void varmath_teardown(void) {
  FC_ReturnCode rc;
  if (isForking) {
    rc = fc_finalLibrary();
    fail_unless(rc == FC_SUCCESS, "failed to final library at end of test");
  }
}

static void varmath_setup_dataset(void) {
  FC_ReturnCode rc;
  int i;
  int numElemPerSide = 3, numVert = 64, numElem = 27, numDim = 3;
  int maxNumData =  192; // 4*4*4 entities * 3 components
  FC_Coords lowers = { 0, 0, 0 }, uppers = { 1, 1, 1 };
  char charData[192];
  int intData[192];
  float fltData[192];
  double dblData[192];
  void* datas[4] = { charData, intData, fltData, dblData };
  char* names[6] = { "vert_scalar_char", "vert_scalar_int", "vert_scalar_flt",
		     "vert_scalar_dbl", "vert_vector_dbl", "elem_scalar_dbl" };
  
  if (isForking) {
    fc_setLibraryVerbosity(fc_messages);
    rc = fc_initLibrary();
    fail_unless(rc == FC_SUCCESS, "test aborted: failed to init library");
  }
  
  // fabricate data
  // no zeros as they mess up some tests (have to check for nan)
  for (i = 0; i < maxNumData; i++) {
    charData[i] = 'a' + i%26;
    intData[i] = 1+i;
    fltData[i] = 1+i/10.;
    dblData[i] = 1+i/100.;
  }

  // create dataset
  rc = fc_createDataset("dataset", &glb_dataset);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create dataset");
  // create global var
  rc = fc_createGlobalVariable(glb_dataset, "global var", &glb_glbVar);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create global var");
  rc = fc_setVariableData(glb_glbVar, 1, 1, FC_AT_WHOLE_DATASET, 
                          FC_MT_SCALAR, FC_DT_DOUBLE, dblData);
  fail_unless(rc == FC_SUCCESS, "abort: failed to set global var data");
  // create mesh
  rc = fc_createSimpleHexMesh(glb_dataset, "hex mesh", numElemPerSide,
			      numElemPerSide, numElemPerSide, lowers,
			      uppers, &glb_mesh);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create hex mesh");
  // creat vertex scalar variables of diff datatypes
  for (i = 0; i < glb_numDataType; i++) {
    rc = fc_createVariable(glb_mesh, names[i], &glb_dataTypeVars[i]);
    fail_unless(rc == FC_SUCCESS, "abort: failed to create var");
    rc = fc_setVariableData(glb_dataTypeVars[i], numVert, 1, FC_AT_VERTEX, 
			    FC_MT_SCALAR, glb_dataTypes[i], datas[i]);
    fail_unless(rc == FC_SUCCESS, "abort: failed to set var data");
  }
  // create vertex vector double variable
  rc = fc_createVariable(glb_mesh, names[4], &glb_vectorVar);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create var");
  rc = fc_setVariableData(glb_vectorVar, numVert, numDim, FC_AT_VERTEX, 
			  FC_MT_VECTOR, FC_DT_DOUBLE, dblData);
  fail_unless(rc == FC_SUCCESS, "abort: failed to set var data");
  // create element scalar double variable
  rc = fc_createVariable(glb_mesh, names[5], &glb_elemVar);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create var");
  rc = fc_setVariableData(glb_elemVar, numElem, 1, FC_AT_ELEMENT, 
			  FC_MT_SCALAR, FC_DT_DOUBLE, dblData);
  fail_unless(rc == FC_SUCCESS, "abort: failed to set var data");
  //fc_printVariable(glb_glbVar, "", 1);
  //for (i = 0; i < glb_numDataType; i++)
    //fc_printVariable(glb_dataTypeVars[i], "", 1);
  //fc_printVariable(glb_vectorVar, "", 1);
  //fc_printVariable(glb_elemVar, "", 1);
}

static void varmath_teardown_dataset(void) {
  FC_ReturnCode rc;

  // cleanup
  rc = fc_deleteDataset(glb_dataset);
  fail_unless(rc == FC_SUCCESS, "abort: failed to delete dataset");

  if (isForking) {
    rc = fc_finalLibrary();
    fail_unless(rc == FC_SUCCESS, "failed to final library at end of test");
  }
}

static void varmath_setup_seqdataset(void) {
  FC_ReturnCode rc;
  int i, j;
  int numElemPerSide = 3, numVert = 64, numElem = 27;
  int temp_numStep, numDim = 3;
  int maxNumData =  195; // (4*4*4 entities * 3 components) + 3 steps
  FC_Coords lowers = { 0, 0, 0 }, uppers = { 1, 1, 1 };
  char charData[195];
  int intData[195];
  float fltData[195];
  double dblData[195];
  double seqCoords[3] = { 0, .1, .25 };
  void* datas[4] = { charData, intData, fltData, dblData };
  char* names[6] = { "vert_scalar_char", "vert_scalar_int", "vert_scalar_flt",
		     "vert_scalar_dbl", "vert_vector_dbl", "elem_scalar_dbl" };
  
  if (isForking) {
    fc_setLibraryVerbosity(fc_messages);
    rc = fc_initLibrary();
    fail_unless(rc == FC_SUCCESS, "test aborted: failed to init library");
  }
  
  // fabricate data
  // no zeros as they mess up some tests (have to check for nan)
  for (i = 0; i < maxNumData; i++) {
    charData[i] = 'a' + i%26;
    intData[i] = 1+i;
    fltData[i] = 1+i/10.;
    dblData[i] = 1+i/100.;
  }

  // create dataset
  rc = fc_createDataset("dataset", &glb_dataset);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create dataset");
  // create sequence
  rc = fc_createSequence(glb_dataset, "sequence", &glb_sequence);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create seuqence");
  rc = fc_setSequenceCoords(glb_sequence, glb_numStep, FC_DT_DOUBLE, seqCoords);
  fail_unless(rc == FC_SUCCESS, "abort: failed to set seq coords");
  // create global seq var
  rc = fc_createGlobalSeqVariable(glb_dataset, glb_sequence, "global seqVar",
                                  &temp_numStep, &glb_glbSeqVar);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create global seqVar");
  for (i = 0; i < glb_numStep; i++) {
    rc = fc_setVariableData(glb_glbSeqVar[i], 1, 1, FC_AT_WHOLE_DATASET,
                            FC_MT_SCALAR, FC_DT_DOUBLE, &dblData[i]);
    fail_unless(rc == FC_SUCCESS, "abort: failed to set glb seq var data");
  }
  // create mesh
  rc = fc_createSimpleHexMesh(glb_dataset, "hex mesh", numElemPerSide,
			      numElemPerSide, numElemPerSide, lowers,
			      uppers, &glb_mesh);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create hex mesh");
  // creat vertex scalar variables of diff datatypes
  for (i = 0; i < glb_numDataType; i++) {
    rc = fc_createSeqVariable(glb_mesh, glb_sequence, names[i], 
			      &temp_numStep, &glb_dataTypeSeqVars[i]);
    fail_unless(rc == FC_SUCCESS, "abort: failed to createSeq var");
    for (j = 0; j < glb_numStep; j++) {
      rc = fc_setVariableData(glb_dataTypeSeqVars[i][j], numVert, 1, 
			      FC_AT_VERTEX, FC_MT_SCALAR, glb_dataTypes[i], 
			      (char*)datas[i]+j*fc_sizeofDataType(glb_dataTypes[i]));
      fail_unless(rc == FC_SUCCESS, "abort: failed to set var data");
    }
  }
  // create vertex vector double variable
  rc = fc_createSeqVariable(glb_mesh, glb_sequence, names[4], 
			    &temp_numStep, &glb_vectorSeqVar);
  fail_unless(rc == FC_SUCCESS, "abort: failed to createSeq var");
  for (j = 0; j < glb_numStep; j++) {
    rc = fc_setVariableData(glb_vectorSeqVar[j], numVert, numDim, FC_AT_VERTEX, 
			    FC_MT_VECTOR, FC_DT_DOUBLE, &dblData[j]);
    fail_unless(rc == FC_SUCCESS, "abort: failed to set var data");
  }
  // create element scalar double variable
  rc = fc_createSeqVariable(glb_mesh, glb_sequence, names[5], 
			    &temp_numStep, &glb_elemSeqVar);
  fail_unless(rc == FC_SUCCESS, "abort: failed to createSeq var");
  for (j = 0; j < glb_numStep; j++) {
    rc = fc_setVariableData(glb_elemSeqVar[j], numElem, 1, FC_AT_ELEMENT, 
			    FC_MT_SCALAR, FC_DT_DOUBLE, &dblData[j]);
    fail_unless(rc == FC_SUCCESS, "abort: failed to set var data");
  }
  //fc_printSeqVariable(glb_numStep, glb_glbSeqVar, "", 1);
  //for (i = 0; i < glb_numDataType; i++)
  //fc_printSeqVariable(glb_numStep, glb_dataTypeSeqVars[i], "", 1);
  //fc_printSeqVariable(glb_numStep, glb_vectorSeqVar, "", 1);
  //fc_printSeqVariable(glb_numStep, glb_elemSeqVar, "", 1);
}

static void varmath_teardown_seqdataset(void) {
  FC_ReturnCode rc;
  int i;

  // cleanup
  rc = fc_deleteDataset(glb_dataset);
  fail_unless(rc == FC_SUCCESS, "abort: failed to delete dataset");
  free(glb_glbSeqVar);
  for (i = 0; i < glb_numDataType; i++)
    free(glb_dataTypeSeqVars[i]);
  free(glb_vectorSeqVar);
  free(glb_elemSeqVar);

  if (isForking) {
    rc = fc_finalLibrary();
    fail_unless(rc == FC_SUCCESS, "failed to final library at end of test");
  }
}

// test unary operations on arrays
START_TEST(unary_array_ops)
{
  FC_ReturnCode rc;
  int i;
  int num = 10;
  int data[10] = { 5, 0, 10, 20, -3, -4000, -9999, 300, 10000, 1} ;
  double ddata[10];
  int results[10];
  double dresults[10];

  // setup doubles just like ints
  for (i = 0; i < num; i++)
    ddata[i] = data[i];

  // --- int operations

  // negatation
  rc = _fc_negateInts(num, data, results);
  fail_unless(rc == FC_SUCCESS, "shouldn't fail to negate ints");
  for (i = 0; i < num; i++) 
    fail_unless(results[i] == -data[i], "_fc_negateInts() wrong result");

  // --- double operations

  // negation
  rc = _fc_negateDoubles(num, ddata, dresults);
  fail_unless(rc == FC_SUCCESS, "shouldn't fail to negate ints");
  for (i = 0; i < num; i++) 
    fail_unless(dresults[i] == -ddata[i], "_fc_negateInts() wrong result");
}
END_TEST

// test binary operations on arrays
START_TEST(binary_array_ops)
{
  FC_ReturnCode rc;
  int i;
  int num = 10;
  int data1[10] = { 5, 0, 10, 20, -3, -4000, -9999, 300, 10000, 1};
  int data2[10] = { 10, 50, 1, -3, -4, 10, -99, 33, -55, 2};
  // (no 0's in data2 because some machines throw floating point exceptions)
  double ddata1[10], ddata2[10];
  int goodResult, result[10];
  double dgoodResult, dresult[10];

  // setup doubles just like ints
  for (i = 0; i < num; i++) {
    ddata1[i] = data1[i];
    ddata2[i] = data2[i];
  }

  // Add ints
  rc = _fc_addInts(1, num, data1, num, data2, result);
  for (i = 0; i < num; i++) {
    goodResult = data1[i] + data2[i];
    fail_unless(result[i] == goodResult, "_fc_addInts() wrong result");
  }

  // subtract ints
  rc = _fc_subtractInts(1, num, data1, num, data2, result);
  for (i = 0; i < num; i++) {
    goodResult = data1[i] - data2[i];
    fail_unless(result[i] == goodResult, "_fc_subtractInts() wrong result");
  }

  // multiply ints
  rc = _fc_multiplyInts(1, num, data1, num, data2, result);
  for (i = 0; i < num; i++) {
    goodResult = data1[i] * data2[i];
    fail_unless(result[i] == goodResult, "_fc_multiplyInts() wrong result");
  }

  // divide ints
  rc = _fc_divideInts(1, num, data1, num, data2, result);
  for (i = 0; i < num; i++) {
    goodResult = data1[i] / data2[i];
    fail_unless(result[i] == goodResult, "_fc_divideInts() wrong result");
  }

  // --- now do the doubles

  // Add doubles
  rc = _fc_addDoubles(1, num, ddata1, num, ddata2, dresult);
  for (i = 0; i < num; i++) {
    dgoodResult = ddata1[i] + ddata2[i];
    fail_unless(dresult[i] == dgoodResult, "_fc_addDoubles() wrong result");
  }

  // subtract doubles
  rc = _fc_subtractDoubles(1, num, ddata1, num, ddata2, dresult);
  for (i = 0; i < num; i++) {
    dgoodResult = ddata1[i] - ddata2[i];
    fail_unless(dresult[i] == dgoodResult, 
		"_fc_subtractDoubles() wrong result");
  }

  // multiply doubles
  rc = _fc_multiplyDoubles(1, num, ddata1, num, ddata2, dresult);
  for (i = 0; i < num; i++) {
    dgoodResult = ddata1[i] * ddata2[i];
    fail_unless(dresult[i] == dgoodResult, 
		"_fc_multiplyDoubles() wrong result");
  }

  // divide doubles
  rc = _fc_divideDoubles(1, num, ddata1, num, ddata2, dresult);
  for (i = 0; i < num; i++) {
    dgoodResult = ddata1[i] / ddata2[i];
    fail_unless(dresult[i] == dgoodResult, "_fc_divideDoubles() wrong result");
  }
}
END_TEST

// test fc_operatorVar()
START_TEST(opVar)
{
  FC_ReturnCode rc;
  int i, j, k;
  FC_Variable var, newVar, badVar = { 999, 999 };
  int numOp = 1;
  char opStrings[1][10] = { "-" };
  int (*intOps[1])(int) = { negInt };
  double (*dblOps[1])(double) = { negDbl };
  int numDataPoint, numComp, temp_numDataPoint, temp_numComp;
  char newVarName[20] = "new var";
  FC_AssociationType assoc, temp_assoc;
  FC_MathType mathType, temp_mathType;
  FC_DataType dataType, temp_dataType;
  void* startData, *resultData;

  // --- use global variables created in the fixture

  // --- test all ops on all datatypes & different # of components
  for (i = 0; i < glb_numDataType + 1; i++) {
    if (i < glb_numDataType)
      var = glb_dataTypeVars[i];
    else
      var = glb_vectorVar;

    // get the start data
    rc = fc_getVariableInfo(var, &numDataPoint, &numComp, &assoc, &mathType, 
			    &dataType);
    fail_unless(rc == FC_SUCCESS, "abort: failed to get var info");
    rc = fc_getVariableDataPtr(var, &startData);
    fail_unless(rc == FC_SUCCESS, "abort: failed to get data pointer");

    for (j = 0; j < numOp; j++) {

      // do it
      rc = fc_operatorVar(var, opStrings[j], newVarName, &newVar);
      if (dataType == FC_DT_CHAR) {
	fail_unless(rc != FC_SUCCESS, "should fail to do char data");
	fail_unless(FC_HANDLE_EQUIV(newVar, FC_NULL_VARIABLE),
		    "fail should return NULL");
	continue;
      }
      else
	fail_unless(rc == FC_SUCCESS, "failed to do unary operation");
      
      // test result meta data
      rc = fc_getVariableInfo(newVar, &temp_numDataPoint, &temp_numComp, 
			      &temp_assoc, &temp_mathType, &temp_dataType);
      fail_unless(rc == FC_SUCCESS, "failed to get new  variable info");
      fail_unless(temp_numDataPoint == numDataPoint, 
		  "mismatch of numDataPoint");
      fail_unless(temp_numComp == numComp, "mismatch of numComp");
      fail_unless(temp_assoc == assoc, "mismatch of assoc");
      fail_unless(temp_mathType == mathType, "mismatch of mathType");
      if (dataType == FC_DT_INT)
	fail_unless(temp_dataType == FC_DT_INT, "mismatch of datatype");
      else
	fail_unless(temp_dataType == FC_DT_DOUBLE, "mismatch of dataType");

      // test result big data
      rc = fc_getVariableDataPtr(newVar, &resultData);
      fail_unless(rc == FC_SUCCESS, "failed to get data ptr from new var");
      for (k = 0; k < numDataPoint*numComp; k++) {
	if (dataType == FC_DT_INT)
	  fail_unless(((int*)resultData)[k] == intOps[j](((int*)startData)[k]),
		      "mismatch of reusult data");
	else if (dataType == FC_DT_FLOAT)
	  fail_unless(((double*)resultData)[k] == dblOps[j](((float*)startData)[k]),
		      "mismatch of reusult data");
	else if (dataType == FC_DT_DOUBLE)
	  fail_unless(((double*)resultData)[k] == dblOps[j](((double*)startData)[k]),
		      "mismatch of reusult data");
      }
      rc = fc_deleteVariable(newVar);
      fail_unless(rc == FC_SUCCESS, "failed to delete new var");
    }
  }

  // --- things should still work with global vars
  rc = fc_operatorVar(glb_glbVar, opStrings[0], newVarName, &newVar);
  fail_unless(rc == FC_SUCCESS, "should work with global var");
  fail_unless(fc_isVariableGlobal(newVar), "should be global");
  rc = fc_getVariableDataPtr(glb_glbVar, &startData);
  fail_unless(rc == FC_SUCCESS, "failed to get data ptr from old var");
  rc = fc_getVariableDataPtr(newVar, &resultData);
  fail_unless(rc == FC_SUCCESS, "failed to get data ptr from new var");
  // just checking first value
  fail_unless(((double*)resultData)[0] == dblOps[0](((double*)startData)[0]),
              "mismatch of reusult data");
  rc = fc_deleteVariable(newVar);
  fail_unless(rc == FC_SUCCESS, "failed to delete new var");

  // --- test bad inputs
  var = glb_dataTypeVars[3];
  rc = fc_operatorVar(badVar, opStrings[0], newVarName, &newVar);
  fail_unless(rc != FC_SUCCESS, "should fail with bad var");
  fail_unless(FC_HANDLE_EQUIV(newVar, FC_NULL_VARIABLE),
	      "fail should return NULL");
  rc = fc_operatorVar(var, "badOpString", newVarName, &newVar);
  fail_unless(rc != FC_SUCCESS, "should fail with bad op string");
  fail_unless(FC_HANDLE_EQUIV(newVar, FC_NULL_VARIABLE),
	      "fail should return NULL");
  rc = fc_operatorVar(var, opStrings[0], NULL, &newVar);
  fail_unless(rc != FC_SUCCESS, "should fail with NULL name");
  fail_unless(FC_HANDLE_EQUIV(newVar, FC_NULL_VARIABLE),
	      "fail should return NULL");
  rc = fc_operatorVar(var, opStrings[0], newVarName, NULL);
  fail_unless(rc != FC_SUCCESS, "should fail if NULL new var");
}
END_TEST

// test fc_varOperatorVar()
START_TEST(varOpVar)
{
  FC_ReturnCode rc;
  int i, j, k, m;
  FC_Variable var1, var2, newVar, badVar = { 999, 999 };
  int numDataPoint1, numComp1, numDataPoint2, numComp2;
  int temp_numDataPoint, temp_numComp;
  char newVarName[20] = "new var";
  FC_AssociationType assoc1, assoc2, temp_assoc;
  FC_MathType mathType1, mathType2, temp_mathType;
  FC_DataType dataType1, dataType2, temp_dataType;
  void* startData1, *startData2, *resultData;
  int numOp = 4;
  char opStrings[16][10] = { "+", "-", "*", "/" };
  int (*intOps[4])(int, int) = { addInts, subInts, multInts, divInts };
  double (*dblOps[4])(double, double) = { addDbls, subDbls, multDbls, divDbls };

  // --- use global variables created in the fixture

  // --- test all operations on all combinations of all var types
  for (i = 0; i < glb_numDataType + 2; i++) {

    if (i < glb_numDataType)
      var1 = glb_dataTypeVars[i];
    else if (i == glb_numDataType)
      var1 = glb_vectorVar;
    else
      var1 = glb_elemVar;

    for (j = 0; j < glb_numDataType + 2; j++) {
      if (j < glb_numDataType)
	var2 = glb_dataTypeVars[j];
      else if (j == glb_numDataType)
	var2 = glb_vectorVar;
      else
	var2 = glb_elemVar;

      // get the start data
      rc = fc_getVariableInfo(var1, &numDataPoint1, &numComp1, &assoc1, &mathType1, 
			      &dataType1);
      fail_unless(rc == FC_SUCCESS, "abort: failed to get var info");
      rc = fc_getVariableDataPtr(var1, &startData1);
      fail_unless(rc == FC_SUCCESS, "abort: failed to get data pointer");
      rc = fc_getVariableInfo(var2, &numDataPoint2, &numComp2, &assoc2, &mathType2, 
			 &dataType2);
      fail_unless(rc == FC_SUCCESS, "abort: failed to get var info");
      rc = fc_getVariableDataPtr(var2, &startData2);
      fail_unless(rc == FC_SUCCESS, "abort: failed to get data pointer");

      // loop over ops
      for (k = 0; k < numOp; k++) {

	// do it
	rc = fc_varOperatorVar(var1, opStrings[k], var2, newVarName, &newVar);
	if (dataType1 == FC_DT_CHAR || dataType2 == FC_DT_CHAR ||
	    numDataPoint1 != numDataPoint2 || numComp1 != numComp2) {
	  fail_unless(rc != FC_SUCCESS, "should fail if char data or if "
		      "vars don't have same numDataPoint or don't have same "
		      "numComp");
	  fail_unless(FC_HANDLE_EQUIV(newVar, FC_NULL_VARIABLE),
		      "fail should return NULL");
	  continue;
	}
	else
	  fail_unless(rc == FC_SUCCESS, "failed to do binary operation");

	// test result meta data
	rc = fc_getVariableInfo(newVar, &temp_numDataPoint, &temp_numComp, 
				&temp_assoc, &temp_mathType, &temp_dataType);
	fail_unless(rc == FC_SUCCESS, "failed to get new  variable info");
	fail_unless(temp_numDataPoint == numDataPoint1, 
		    "mismatch of numDataPoint");
	fail_unless(temp_numComp == numComp1, "mismatch of numComp");
	fail_unless(temp_assoc == assoc1, "mismatch of assoc");
	fail_unless(temp_mathType == mathType1, "mismatch of mathType");
	if (dataType1 == FC_DT_INT && dataType2 == FC_DT_INT)
	  fail_unless(temp_dataType == FC_DT_INT, "mismatch of dataType");
	else
	  fail_unless(temp_dataType == FC_DT_DOUBLE, "mismatch of dataType");

	// test result big data
	rc = fc_getVariableDataPtr(newVar, &resultData);
	fail_unless(rc == FC_SUCCESS, "failed to get data ptr from new var");
	for (m = 0; m < numDataPoint1*numComp1; m++) {
	  if (temp_dataType == FC_DT_INT) {
	    fail_unless(((int*)resultData)[m] ==
		   intOps[k](((int*)startData1)[m], ((int*)startData2)[m]),
			"mimstch of result data");
	  }
	  else {
	    double start1, start2;
	    switch (dataType1) {
	    case FC_DT_INT:     start1 =    ((int*)startData1)[m];  break;
	    case FC_DT_FLOAT:   start1 =  ((float*)startData1)[m];  break;
	    case FC_DT_DOUBLE:  start1 = ((double*)startData1)[m];  break;
	    default: ;
	    }
	    switch (dataType2) {
	    case FC_DT_INT:     start2 =    ((int*)startData2)[m];  break;
	    case FC_DT_FLOAT:   start2 =  ((float*)startData2)[m];  break;
	    case FC_DT_DOUBLE:  start2 = ((double*)startData2)[m];  break;
	    default: ;
	    }
	    fail_unless(FC_DBL_EQUIV(((double*)resultData)[m],
				     dblOps[k](start1, start2)),
			"mismatch of result data");
	  }
	}
	rc = fc_deleteVariable(newVar);
	fail_unless(rc == FC_SUCCESS, "failed to delete new var");
      }
    }
  }

  // --- things should still work with global vars
  rc = fc_varOperatorVar(glb_glbVar, opStrings[0], glb_glbVar, newVarName, &newVar);
  fail_unless(rc == FC_SUCCESS, "should work with global var");
  fail_unless(fc_isVariableGlobal(newVar), "should be global");
  rc = fc_getVariableDataPtr(glb_glbVar, &startData1);
  fail_unless(rc == FC_SUCCESS, "failed to get data ptr from old var");
  rc = fc_getVariableDataPtr(newVar, &resultData);
  fail_unless(rc == FC_SUCCESS, "failed to get data ptr from new var");
  // just checking first value
  fail_unless(((double*)resultData)[0] == dblOps[0](((double*)startData1)[0],
                                                    ((double*)startData1)[0]),
              "mismatch of reusult data");
  rc = fc_deleteVariable(newVar);
  fail_unless(rc == FC_SUCCESS, "failed to delete new var");

  // --- things should fail with mixed global non-global vars
  var1 = glb_dataTypeVars[3];
  rc = fc_varOperatorVar(glb_glbVar, opStrings[0], var1, newVarName, &newVar);
  fail_unless(rc != FC_SUCCESS, "should fail with mix of global/nonglobal");
  fail_unless(FC_HANDLE_EQUIV(newVar, FC_NULL_VARIABLE),
	      "fail should return NULL");
  rc = fc_varOperatorVar(var1, opStrings[0], glb_glbVar, newVarName, &newVar);
  fail_unless(rc != FC_SUCCESS, "should fail with mix of global/nonglobal");
  fail_unless(FC_HANDLE_EQUIV(newVar, FC_NULL_VARIABLE),
	      "fail should return NULL");
   
  // --- test bad inputs
  var1 = var2 = glb_dataTypeVars[3];
  rc = fc_varOperatorVar(badVar, opStrings[0], var2, newVarName, &newVar);
  fail_unless(rc != FC_SUCCESS, "should fail with bad var 1");
  fail_unless(FC_HANDLE_EQUIV(newVar, FC_NULL_VARIABLE),
	      "fail should return NULL");
  rc = fc_varOperatorVar(var1, "bad op", var2, newVarName, &newVar);
  fail_unless(rc != FC_SUCCESS, "should fail with bad op string");
  fail_unless(FC_HANDLE_EQUIV(newVar, FC_NULL_VARIABLE),
	      "fail should return NULL");
  rc = fc_varOperatorVar(var1, opStrings[0], badVar, newVarName, &newVar);
  fail_unless(rc != FC_SUCCESS, "should fail with bad var 2");
  fail_unless(FC_HANDLE_EQUIV(newVar, FC_NULL_VARIABLE),
	      "fail should return NULL");
  rc = fc_varOperatorVar(var1, opStrings[0], var2, NULL, &newVar);
  fail_unless(rc != FC_SUCCESS, "should fail with null new name");
  fail_unless(FC_HANDLE_EQUIV(newVar, FC_NULL_VARIABLE),
	      "fail should return NULL");
  rc = fc_varOperatorVar(var1, opStrings[0], var2, newVarName, NULL);
  fail_unless(rc != FC_SUCCESS, "should fail with null newVar");
}
END_TEST

// test fc_varOperatorConst()  
START_TEST(varOpConst)
{
  FC_ReturnCode rc;
  int i, j, k, m, n;
  FC_Variable var, newVar, badVar = { 999, 999 };
  int numDataPoint, numComp, temp_numDataPoint, temp_numComp;
  char newVarName[20] = "new var";
  FC_AssociationType assoc, temp_assoc;
  FC_MathType mathType, temp_mathType;
  FC_DataType dataType, temp_dataType;
  void* startData, *resultData;
  int numOp = 4;
  char opStrings[16][10] = { "+", "-", "*", "/" };
  int (*intOps[4])(int, int) = { addInts, subInts, multInts, divInts };
  double (*dblOps[4])(double, double) = { addDbls, subDbls, multDbls, divDbls };
  int numConst = 4*2; // four data types * whether scalar or vector
  char charConst[3] = { 'a', 'z', '$' };
  int intConst[3] = { 20, -2, 6 };
  float fltConst[3] = { 3.5, -99.99, 1.001 };
  double dblConst[3] = { 1.2, .00001, -29.9 };
  int constNumComps[4*2] = { 1, 3, 1, 3, 1, 3, 1, 3 };
  void *constValues[4*2] = { &charConst, &charConst, &intConst, &intConst,
			     &fltConst, &fltConst, &dblConst, &dblConst };
  FC_DataType constTypes[4*2] = { FC_DT_CHAR, FC_DT_CHAR, FC_DT_INT,
				   FC_DT_INT, FC_DT_FLOAT, FC_DT_FLOAT,
				   FC_DT_DOUBLE, FC_DT_DOUBLE };

  // --- use global variables created in the fixture

  // --- test alls operation on all combinations of vars & consts
  for (i = 0; i < glb_numDataType + 1; i++) {
    if (i < glb_numDataType)
      var = glb_dataTypeVars[i];
    else
      var = glb_vectorVar;
    
    for (j = 0; j < numConst; j++) {
   
      // get the start data
      rc = fc_getVariableInfo(var, &numDataPoint, &numComp, &assoc, &mathType, 
			      &dataType);
      fail_unless(rc == FC_SUCCESS, "abort: failed to get var info");
      rc = fc_getVariableDataPtr(var, &startData);
      fail_unless(rc == FC_SUCCESS, "abort: failed to get data pointer");

      // loop over ops
      for (k = 0; k < numOp; k++) {

	// do it
	rc = fc_varOperatorConst(var, opStrings[k], constNumComps[j],
				 constValues[j], constTypes[j], newVarName, 
				 &newVar);
	if (dataType == FC_DT_CHAR || constTypes[j] == FC_DT_CHAR ||
	    numComp != constNumComps[j]) {
	  fail_unless(rc != FC_SUCCESS, "should fail if char data or if "
		      "don't have same numComp");
	  fail_unless(FC_HANDLE_EQUIV(newVar, FC_NULL_VARIABLE),
		      "fail should return NULL");
	  continue;
	}
	else
	  fail_unless(rc == FC_SUCCESS, "failed to do binary operation");

	// test result meta data
	rc = fc_getVariableInfo(newVar, &temp_numDataPoint, &temp_numComp, 
				&temp_assoc, &temp_mathType, &temp_dataType);
	fail_unless(rc == FC_SUCCESS, "failed to get new  variable info");
	fail_unless(temp_numDataPoint == numDataPoint, 
		    "mismatch of numDataPoint");
	fail_unless(temp_numComp == numComp, "mismatch of numComp");
	fail_unless(temp_assoc == assoc, "mismatch of assoc");
	fail_unless(temp_mathType == mathType, "mismatch of mathType");
	if (dataType == FC_DT_INT && constTypes[j] == FC_DT_INT)
	  fail_unless(temp_dataType == FC_DT_INT, "mismatch of dataType");
	else
	  fail_unless(temp_dataType == FC_DT_DOUBLE, "mismatch of dataType");

	// test result big data
	rc = fc_getVariableDataPtr(newVar, &resultData);
	fail_unless(rc == FC_SUCCESS, "failed to get data ptr from new var");
	for (m = 0; m < numDataPoint; m++) {
	  for (n = 0; n < numComp; n++) {
	    int idx = m*numComp + n;
	    if (temp_dataType == FC_DT_INT) {
	      fail_unless(((int*)resultData)[idx] ==
			  intOps[k](((int*)startData)[idx], intConst[n]),
			  "mismatch of result data");
	    }
	    else {
	      double varStart, constStart;
	      switch (dataType) {
	      case FC_DT_INT:    varStart =    ((int*)startData)[idx]; break;
	      case FC_DT_FLOAT:  varStart =  ((float*)startData)[idx]; break;
	      case FC_DT_DOUBLE: varStart = ((double*)startData)[idx]; break;
	      default: ;
	      }
	      switch (constTypes[j]) {
	      case FC_DT_INT:    constStart =    intConst[n];  break;
	      case FC_DT_FLOAT:  constStart =    fltConst[n];  break;
	      case FC_DT_DOUBLE: constStart =    dblConst[n];  break;
	      default: ;
	      }
	      fail_unless(FC_DBL_EQUIV(((double*)resultData)[idx], 
				       dblOps[k](varStart, constStart)),
			  "mismatch of result data");
	    }
	  }
	}
	rc = fc_deleteVariable(newVar);
	fail_unless(rc == FC_SUCCESS, "failed to delete new var");
      }
    }
  }
 
  // --- things should still work with global vars
  rc = fc_varOperatorConst(glb_glbVar, opStrings[0], 1, dblConst, FC_DT_DOUBLE,
                           newVarName, &newVar);
  fail_unless(rc == FC_SUCCESS, "should work with global var");
  fail_unless(fc_isVariableGlobal(newVar), "should be global");
  rc = fc_getVariableDataPtr(glb_glbVar, &startData);
  fail_unless(rc == FC_SUCCESS, "failed to get data ptr from old var");
  rc = fc_getVariableDataPtr(newVar, &resultData);
  fail_unless(rc == FC_SUCCESS, "failed to get data ptr from new var");
  // just checking first value
  fail_unless(FC_DBL_EQUIV(((double*)resultData)[0],
			   dblOps[0](((double*)startData)[0], dblConst[0])),
              "mismatch of reusult data");
  rc = fc_deleteVariable(newVar);
  fail_unless(rc == FC_SUCCESS, "failed to delete new var");

  // --- test bad inputs
  var = glb_dataTypeVars[3];
  rc = fc_varOperatorConst(badVar, opStrings[0], 1, dblConst, FC_DT_DOUBLE, 
			   newVarName, &newVar);
  fail_unless(rc != FC_SUCCESS, "should fail with bad var");
  fail_unless(FC_HANDLE_EQUIV(newVar, FC_NULL_VARIABLE),
	      "fail should return NULL");
  rc = fc_varOperatorConst(var, "bad op", 1, dblConst, FC_DT_DOUBLE,
			 newVarName, &newVar);
  fail_unless(rc != FC_SUCCESS, "should fail with bad op string");
  fail_unless(FC_HANDLE_EQUIV(newVar, FC_NULL_VARIABLE),
	      "fail should return NULL");
  rc = fc_varOperatorConst(var, opStrings[0], 0, dblConst, FC_DT_DOUBLE,
			   newVarName, &newVar);
  fail_unless(rc != FC_SUCCESS, "should fail if constLen < 1");
  fail_unless(FC_HANDLE_EQUIV(newVar, FC_NULL_VARIABLE),
	      "fail should return NULL");
  rc = fc_varOperatorConst(var, opStrings[0], 1, NULL, FC_DT_DOUBLE,
			   newVarName, &newVar);
  fail_unless(rc != FC_SUCCESS, "should fail if constValues = NULL");
  fail_unless(FC_HANDLE_EQUIV(newVar, FC_NULL_VARIABLE),
	      "fail should return NULL");
  rc = fc_varOperatorConst(var, opStrings[0], 1, dblConst, FC_DT_UNKNOWN,
			   newVarName, &newVar);
  fail_unless(rc != FC_SUCCESS, "should fail if data type unknown");
  fail_unless(FC_HANDLE_EQUIV(newVar, FC_NULL_VARIABLE),
	      "fail should return NULL");
  rc = fc_varOperatorConst(var, opStrings[0], 1, dblConst, -99,
			   newVarName, &newVar);
  fail_unless(rc != FC_SUCCESS, "should fail if bad datatype");
  fail_unless(FC_HANDLE_EQUIV(newVar, FC_NULL_VARIABLE),
	      "fail should return NULL");
  rc = fc_varOperatorConst(var, opStrings[0], 1, dblConst, FC_DT_DOUBLE,
			   NULL, &newVar);
  fail_unless(rc != FC_SUCCESS, "should fail with null new name");
  fail_unless(FC_HANDLE_EQUIV(newVar, FC_NULL_VARIABLE),
	      "fail should return NULL");
  rc = fc_varOperatorConst(var, opStrings[0], 1, dblConst, FC_DT_DOUBLE,
			   newVarName, NULL);
  fail_unless(rc != FC_SUCCESS, "should fail with null newVar");
}
END_TEST

// test fc_constOperatorVar() 
START_TEST(constOpVar)
{
  FC_ReturnCode rc;
  int i, j, k, m, n;
  FC_Variable var, newVar, badVar = { 999, 999 };
  int numDataPoint, numComp, temp_numDataPoint, temp_numComp;
  char newVarName[20] = "new var";
  FC_AssociationType assoc, temp_assoc;
  FC_MathType mathType, temp_mathType;
  FC_DataType dataType, temp_dataType;
  void* startData, *resultData;
  int numOp = 4;
  char opStrings[16][10] = { "+", "-", "*", "/" };
  int (*intOps[4])(int, int) = { addInts, subInts, multInts, divInts };
  double (*dblOps[4])(double, double) = { addDbls, subDbls, multDbls, divDbls };
  int numConst = 4*2;
  char charConst[3] = { 'a', 'z', '$' };
  int intConst[3] = { 20, -2, 6 };
  float fltConst[3] = { 3.5, -99.99, 1.001 };
  double dblConst[3] = { 1.2, .00001, -29.9 };
  int constNumComps[4*2] = { 1, 3, 1, 3, 1, 3, 1, 3 };
  void *constValues[4*2] = { &charConst, &charConst, &intConst, &intConst,
			     &fltConst, &fltConst, &dblConst, &dblConst };
  FC_DataType constTypes[4*2] = { FC_DT_CHAR, FC_DT_CHAR, FC_DT_INT,
				   FC_DT_INT, FC_DT_FLOAT, FC_DT_FLOAT,
				   FC_DT_DOUBLE, FC_DT_DOUBLE };

  // --- use global variables created in the fixture

  // --- test alls operation on all combinations of vars & consts
  for (i = 0; i < glb_numDataType + 1; i++) {
    if (i < glb_numDataType)
      var = glb_dataTypeVars[i];
    else
      var = glb_vectorVar;

    for (j = 0; j < numConst; j++) {

      // get the start data
      fc_getVariableInfo(var, &numDataPoint, &numComp, &assoc, &mathType, 
			 &dataType);
      rc = fc_getVariableDataPtr(var, &startData);
      fail_unless(rc == FC_SUCCESS, "abort: failed to get data pointer");

      // loop over ops
      for (k = 0; k < numOp; k++) {

      // do it
	rc = fc_constOperatorVar(constNumComps[j], constValues[j], constTypes[j], 
				 opStrings[k], var, newVarName, &newVar);
	if (dataType == FC_DT_CHAR || constTypes[j] == FC_DT_CHAR ||
	    numComp != constNumComps[j]) {
	  fail_unless(rc != FC_SUCCESS, "should fail if char data or if "
		      "don't have same numComp");
	  fail_unless(FC_HANDLE_EQUIV(newVar, FC_NULL_VARIABLE),
		      "fail should return NULL");
	  continue;
	}
	else
	  fail_unless(rc == FC_SUCCESS, "failed to do binary operation");

	// test result meta data
	rc = fc_getVariableInfo(newVar, &temp_numDataPoint, &temp_numComp, 
				&temp_assoc, &temp_mathType, &temp_dataType);
	fail_unless(rc == FC_SUCCESS, "failed to get new  variable info");
	fail_unless(temp_numDataPoint == numDataPoint, 
		    "mismatch of numDataPoint");
	fail_unless(temp_numComp == numComp, "mismatch of numComp");
	fail_unless(temp_assoc == assoc, "mismatch of assoc");
	fail_unless(temp_mathType == mathType, "mismatch of mathType");
	if (dataType == FC_DT_INT && constTypes[j] == FC_DT_INT)
	  fail_unless(temp_dataType == FC_DT_INT, "mismatch of dataType");
	else
	  fail_unless(temp_dataType == FC_DT_DOUBLE, "mismatch of dataType");

	// test result big data
	rc = fc_getVariableDataPtr(newVar, &resultData);
	fail_unless(rc == FC_SUCCESS, "failed to get data ptr from new var");
	for (m = 0; m < numDataPoint; m++) {
	  for (n = 0; n < numComp; n++) {
	    int idx = m*numComp + n;
	    if (temp_dataType == FC_DT_INT) {
	      fail_unless(((int*)resultData)[idx] ==
			  intOps[k](intConst[n], ((int*)startData)[idx]),
			  "mismatch of result data");
	    }
	    else {
	      double varStart, constStart;
	      switch (dataType) {
	      case FC_DT_INT:    varStart =    ((int*)startData)[idx]; break;
	      case FC_DT_FLOAT:  varStart =  ((float*)startData)[idx]; break;
	      case FC_DT_DOUBLE: varStart = ((double*)startData)[idx]; break;
	      default: ;
	      }
	      switch (constTypes[j]) {
	      case FC_DT_INT:    constStart =    intConst[n];  break;
	      case FC_DT_FLOAT:  constStart =    fltConst[n];  break;
	      case FC_DT_DOUBLE: constStart =    dblConst[n];  break;
	      default: ;
	      }
	      fail_unless(FC_DBL_EQUIV(((double*)resultData)[idx], 
				       dblOps[k](constStart, varStart)),
			  "mismatch of result data");
	    }
	  }
	}
	rc = fc_deleteVariable(newVar);
	fail_unless(rc == FC_SUCCESS, "failed to delete new var");
      }
    }
  }

  // --- things should still work with global vars
  rc = fc_constOperatorVar(1, dblConst, FC_DT_DOUBLE, opStrings[0], glb_glbVar,
                           newVarName, &newVar);
  fail_unless(rc == FC_SUCCESS, "should work with global var");
  fail_unless(fc_isVariableGlobal(newVar), "should be global");
  rc = fc_getVariableDataPtr(glb_glbVar, &startData);
  fail_unless(rc == FC_SUCCESS, "failed to get data ptr from old var");
  rc = fc_getVariableDataPtr(newVar, &resultData);
  fail_unless(rc == FC_SUCCESS, "failed to get data ptr from new var");
  // just checking first value
  fail_unless(FC_DBL_EQUIV(((double*)resultData)[0],
			   dblOps[0](dblConst[0],((double*)startData)[0])),
              "mismatch of reusult data");
  rc = fc_deleteVariable(newVar);
  fail_unless(rc == FC_SUCCESS, "failed to delete new var");

  // --- test bad inputs
  var = glb_dataTypeVars[3];
  rc = fc_constOperatorVar(0, dblConst, FC_DT_DOUBLE, opStrings[0], var,
			   newVarName, &newVar);
  fail_unless(rc != FC_SUCCESS, "should fail if constLen < 1");
  fail_unless(FC_HANDLE_EQUIV(newVar, FC_NULL_VARIABLE),
	      "fail should return NULL");
  rc = fc_constOperatorVar(1, NULL, FC_DT_DOUBLE, opStrings[0], var,
			   newVarName, &newVar);
  fail_unless(rc != FC_SUCCESS, "should fail if constValues = NULL");
  fail_unless(FC_HANDLE_EQUIV(newVar, FC_NULL_VARIABLE),
	      "fail should return NULL");
  rc = fc_constOperatorVar(1, dblConst, FC_DT_UNKNOWN, opStrings[0], var,
			   newVarName, &newVar);
  fail_unless(rc != FC_SUCCESS, "should fail if data type unknown");
  fail_unless(FC_HANDLE_EQUIV(newVar, FC_NULL_VARIABLE),
	      "fail should return NULL");
  rc = fc_constOperatorVar(1, dblConst, -99, opStrings[0], var,
			   newVarName, &newVar);
  fail_unless(rc != FC_SUCCESS, "should fail if bad datatype");
  fail_unless(FC_HANDLE_EQUIV(newVar, FC_NULL_VARIABLE),
	      "fail should return NULL");
  rc = fc_constOperatorVar(1, dblConst, FC_DT_DOUBLE, "bad op", var,
			 newVarName, &newVar);
  fail_unless(rc != FC_SUCCESS, "should fail with bad op string");
  fail_unless(FC_HANDLE_EQUIV(newVar, FC_NULL_VARIABLE),
	      "fail should return NULL");
  rc = fc_constOperatorVar(1, dblConst, FC_DT_DOUBLE, opStrings[0], badVar,
			   newVarName, &newVar);
  fail_unless(rc != FC_SUCCESS, "should fail with bad var");
  fail_unless(FC_HANDLE_EQUIV(newVar, FC_NULL_VARIABLE),
	      "fail should return NULL");
  rc = fc_constOperatorVar(1, dblConst, FC_DT_DOUBLE, opStrings[0], var,
			   NULL, &newVar);
  fail_unless(rc != FC_SUCCESS, "should fail with null new name");
  fail_unless(FC_HANDLE_EQUIV(newVar, FC_NULL_VARIABLE),
	      "fail should return NULL");
  rc = fc_constOperatorVar(1, dblConst, FC_DT_DOUBLE, opStrings[0],var,  
			   newVarName, NULL);
  fail_unless(rc != FC_SUCCESS, "should fail with null newVar");
}
END_TEST

// test fc_operatorSeqVar()
START_TEST(opSeqVar)
{
  FC_ReturnCode rc;
  int i, j, k, m;
  FC_Sequence temp_seq;
  FC_Variable *seqVar, *newSeqVar, badSeqVar[glb_numStep], badVar = { 999, 999 };
  int numOp = 1;
  char opStrings[1][10] = { "-" };
  int numDataPoint, numComp, temp_numDataPoint, temp_numComp;
  char newSeqVarName[20] = "new seq var";
  FC_AssociationType assoc, temp_assoc;
  FC_MathType mathType, temp_mathType;
  FC_DataType dataType, temp_dataType;
  void* startData, *resultData;
  int (*intOps[1])(int) = { negInt };
  double (*dblOps[1])(double) = { negDbl };

  // --- use global variables created in the fixture

  // --- test all ops on all datatypes & different # of components
  for (i = 0; i < glb_numDataType + 1; i++) {
    if (i < glb_numDataType)
      seqVar = glb_dataTypeSeqVars[i];
    else
      seqVar = glb_vectorSeqVar;

    // get the start metadata
    rc = fc_getVariableInfo(seqVar[0], &numDataPoint, &numComp, &assoc,
			    &mathType, &dataType);
    fail_unless(rc == FC_SUCCESS, "cant get var info");

    for (j = 0; j < numOp; j++) {

      // do it
      rc = fc_operatorSeqVar(glb_numStep, seqVar, opStrings[j],
			     newSeqVarName, &newSeqVar);
      if (dataType == FC_DT_CHAR) {
	fail_unless(rc != FC_SUCCESS, "should fail to do char data");
	fail_unless(newSeqVar == NULL, "fail should return NULL");
	continue;
      }
      else
	fail_unless(rc == FC_SUCCESS, "failed to do unary operation");
      
      // test result meta data
      rc = fc_getVariableInfo(newSeqVar[0], &temp_numDataPoint, &temp_numComp, 
			      &temp_assoc, &temp_mathType, &temp_dataType);
      fail_unless(rc == FC_SUCCESS, "failed to get new variable info");
      fail_unless(temp_numDataPoint == numDataPoint, 
		  "mismatch of numDataPoint");
      fail_unless(temp_numComp == numComp, "mismatch of numComp");
      fail_unless(temp_assoc == assoc, "mismatch of assoc");
      fail_unless(temp_mathType == mathType, "mismatch of mathType");
      if (dataType == FC_DT_INT)
	fail_unless(temp_dataType == FC_DT_INT, "mismatch of datatype");
      else
	fail_unless(temp_dataType == FC_DT_DOUBLE, "mismatch of dataType");
      rc = fc_getSequenceFromSeqVariable(glb_numStep, newSeqVar, &temp_seq);
      fail_unless(rc == FC_SUCCESS, "abort: failed to get new sequence");
      fail_unless(FC_HANDLE_EQUIV(temp_seq, glb_sequence), "mismatch of sequence");

      // test result big data
      for (k = 0; k < glb_numStep; k++) {
	rc = fc_getVariableDataPtr(seqVar[k], &startData);
	fail_unless(rc == FC_SUCCESS, "abort: failed to get data pointer");
	rc = fc_getVariableDataPtr(newSeqVar[k], &resultData);
	fail_unless(rc == FC_SUCCESS, "failed to get data ptr from new var");
	for (m = 0; m < numDataPoint*numComp; m++) {
	  if (dataType == FC_DT_INT)
	    fail_unless(((int*)resultData)[m] == intOps[j](((int*)startData)[m]),
			"mismatch of reusult data");
	  else if (dataType == FC_DT_FLOAT)
	    fail_unless(((double*)resultData)[m] == dblOps[j](((float*)startData)[m]),
			"mismatch of reusult data");
	  else if (dataType == FC_DT_DOUBLE)
	    fail_unless(((double*)resultData)[m] == dblOps[j](((double*)startData)[m]),
			"mismatch of reusult data");
	}
      }
      rc = fc_deleteSeqVariable(glb_numStep, newSeqVar);
      fail_unless(rc == FC_SUCCESS, "failed to delete new seqvar");
      free(newSeqVar);
    }
  }
 
  // --- things should still work with global vars
  rc = fc_operatorSeqVar(glb_numStep, glb_glbSeqVar, opStrings[0], 
                         newSeqVarName, &newSeqVar);
  fail_unless(rc == FC_SUCCESS, "should work with global var");
  fail_unless(fc_isVariableGlobal(newSeqVar[0]), "should be global");
  for (i = 0; i < glb_numStep; i++) {
    rc = fc_getVariableDataPtr(glb_glbSeqVar[i], &startData);
    fail_unless(rc == FC_SUCCESS, "failed to get data ptr from old var");
    rc = fc_getVariableDataPtr(newSeqVar[i], &resultData);
    fail_unless(rc == FC_SUCCESS, "failed to get data ptr from new var");
    // just checking first value
    fail_unless(((double*)resultData)[0] == dblOps[0](((double*)startData)[0]),
                "mismatch of result data");
  }
  rc = fc_deleteSeqVariable(glb_numStep, newSeqVar);
  fail_unless(rc == FC_SUCCESS, "failed to delete new var");
  free(newSeqVar);

  // --- test bad inputs
  seqVar = glb_dataTypeSeqVars[3];
  for (i = 0; i < glb_numStep; i++)
    badSeqVar[i] = glb_dataTypeSeqVars[3][i];
  badSeqVar[1] = badVar;
  rc = fc_operatorSeqVar(glb_numStep, NULL, opStrings[0], newSeqVarName,
			 &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail with bad input var");
  fail_unless(newSeqVar == NULL,"fail should return NULL");
   rc = fc_operatorSeqVar(glb_numStep, badSeqVar, opStrings[0], newSeqVarName,
			 &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail with bad input var");
  fail_unless(newSeqVar == NULL,"fail should return NULL");
  rc = fc_operatorSeqVar(glb_numStep, seqVar, "badOpString", newSeqVarName,
			 &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail with bad op string");
  fail_unless(newSeqVar == NULL,"fail should return NULL");
  rc = fc_operatorSeqVar(glb_numStep, seqVar, opStrings[0], NULL, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail with NULL name");
  fail_unless(newSeqVar == NULL,"fail should return NULL");
  rc = fc_operatorSeqVar(glb_numStep, seqVar, opStrings[0], newSeqVarName, NULL);
  fail_unless(rc != FC_SUCCESS, "should fail if NULL new var");
}
END_TEST

// test fc_seqVarOperatorSeqVar()
START_TEST(seqVarOpSeqVar)
{
  FC_ReturnCode rc;
  int i, j, k, m, n;
  FC_Sequence temp_seq, seq2;
  FC_Variable *seqVar1, *seqVar2, *newSeqVar, badSeqVar[3], badVar = { 999, 999 };
  FC_Variable *otherSeqVar;
  int numDataPoint1, numComp1, numDataPoint2, numComp2;
  int temp_numDataPoint, temp_numComp;
  char newSeqVarName[20] = "new seq var";
  FC_AssociationType assoc1, assoc2, temp_assoc;
  FC_MathType mathType1, mathType2, temp_mathType;
  FC_DataType dataType1, dataType2, temp_dataType;
  void* startData1, *startData2, *resultData;
  int numOp = 4;
  char opStrings[16][10] = { "+", "-", "*", "/" };
  int (*intOps[4])(int, int) = { addInts, subInts, multInts, divInts };
  double (*dblOps[4])(double, double) = { addDbls, subDbls, multDbls, divDbls };

  // --- use global variables created in the fixture

  // --- test all ops on all combinations of all var types
  for (i = 0; i < glb_numDataType + 2; i++) {

    if (i < glb_numDataType)
      seqVar1 = glb_dataTypeSeqVars[i];
    else if (i == glb_numDataType)
      seqVar1 = glb_vectorSeqVar;
    else 
      seqVar1 = glb_vectorSeqVar;

    for (j = 0; j < glb_numDataType + 2; j++) {
      if (j < glb_numDataType)
	seqVar2 = glb_dataTypeSeqVars[j];
      else if (j == glb_numDataType)
	seqVar2 = glb_vectorSeqVar;
      else
	seqVar2 = glb_elemSeqVar;

      // get the start meta data
      rc = fc_getVariableInfo(seqVar1[0], &numDataPoint1, &numComp1, &assoc1, 
			      &mathType1,  &dataType1);
      fail_unless(rc == FC_SUCCESS, "abort: failed to get var info");
      rc = fc_getVariableInfo(seqVar2[0], &numDataPoint2, &numComp2, &assoc2, 
			      &mathType2, &dataType2);
      fail_unless(rc == FC_SUCCESS, "abort: failed to get var info");

      // loop over ops
      for (k = 0; k < numOp; k++) {

	// do it
	rc = fc_seqVarOperatorSeqVar(glb_numStep, seqVar1, opStrings[k],
				     seqVar2, newSeqVarName, &newSeqVar);
	if (dataType1 == FC_DT_CHAR || dataType2 == FC_DT_CHAR ||
	    numDataPoint1 != numDataPoint2 || numComp1 != numComp2) {
	  fail_unless(rc != FC_SUCCESS, "should fail to do char data or if "
		      "vars don't have some numDataPoint or don't have same "
		      "numComp");
	  fail_unless(newSeqVar == NULL, "fail should return NULL");
	  continue;
	}
	else
	  fail_unless(rc == FC_SUCCESS, "failed to do binary operation");
       
	// test result meta data
	rc = fc_getVariableInfo(newSeqVar[0], &temp_numDataPoint, &temp_numComp, 
			      &temp_assoc, &temp_mathType, &temp_dataType);
	fail_unless(rc == FC_SUCCESS, "failed to get new variable info");
	fail_unless(temp_numDataPoint == numDataPoint1, 
		    "mismatch of numDataPoint");
	fail_unless(temp_numComp == numComp1, "mismatch of numComp");
	fail_unless(temp_assoc == assoc1, "mismatch of assoc");
	fail_unless(temp_mathType == mathType1, "mismatch of mathType");
	if (dataType1 == FC_DT_INT && dataType2 == FC_DT_INT)
	  fail_unless(temp_dataType == FC_DT_INT, "mismatch of datatype");
	else
	  fail_unless(temp_dataType == FC_DT_DOUBLE, "mismatch of dataType");
	rc = fc_getSequenceFromSeqVariable(glb_numStep, newSeqVar, &temp_seq);
	fail_unless(rc == FC_SUCCESS, "abort: failed to get new sequence");
	fail_unless(FC_HANDLE_EQUIV(temp_seq, glb_sequence), "mismatch of sequence");

	// test result big data
	for (m = 0; m < glb_numStep; m++) {	
	  rc = fc_getVariableDataPtr(seqVar1[m], &startData1);
	  fail_unless(rc == FC_SUCCESS, "abort: failed to get data pointer");
	  rc = fc_getVariableDataPtr(seqVar2[m], &startData2);
	  fail_unless(rc == FC_SUCCESS, "abort: failed to get data pointer");
	  rc = fc_getVariableDataPtr(newSeqVar[m], &resultData);
	  fail_unless(rc == FC_SUCCESS, "failed to get data ptr from new var");
	  for (n = 0; n < numDataPoint1*numComp1; n++) {
	    if (temp_dataType == FC_DT_INT) {
	      fail_unless(((int*)resultData)[n] == 
			  intOps[k](((int*)startData1)[n], ((int*)startData2)[n]),
			  "mismatch of reusult data");
	    }
	    else {
	      double start1, start2;
	      switch (dataType1) {
	      case FC_DT_INT:     start1 =    ((int*)startData1)[n];  break;
	      case FC_DT_FLOAT:   start1 =  ((float*)startData1)[n];  break;
	      case FC_DT_DOUBLE:  start1 = ((double*)startData1)[n];  break;
	      default: ;
	      }
	      switch (dataType2) {
	      case FC_DT_INT:     start2 =    ((int*)startData2)[n];  break;
	      case FC_DT_FLOAT:   start2 =  ((float*)startData2)[n];  break;
	      case FC_DT_DOUBLE:  start2 = ((double*)startData2)[n];  break;
	      default: ;
	      }
	      fail_unless(FC_DBL_EQUIV(((double*)resultData)[n],
				       dblOps[k](start1, start2)),
			  "mismatch of result data");
	    }
	  }
	}
	rc = fc_deleteSeqVariable(glb_numStep, newSeqVar);
	fail_unless(rc == FC_SUCCESS, "failed to delete new seqvar");
	free(newSeqVar);
      }
    }
  }

  // --- things should still work with global vars
  rc = fc_seqVarOperatorSeqVar(glb_numStep, glb_glbSeqVar, opStrings[0], 
                               glb_glbSeqVar, newSeqVarName, &newSeqVar);
  fail_unless(rc == FC_SUCCESS, "should work with global var");
  fail_unless(fc_isVariableGlobal(newSeqVar[0]), "should be global");
  for (i = 0; i < glb_numStep; i++) {
    rc = fc_getVariableDataPtr(glb_glbSeqVar[i], &startData1);
    fail_unless(rc == FC_SUCCESS, "failed to get data ptr from old var");
    rc = fc_getVariableDataPtr(newSeqVar[i], &resultData);
    fail_unless(rc == FC_SUCCESS, "failed to get data ptr from new var");
    // just checking first value
    fail_unless(((double*)resultData)[0] == dblOps[0](((double*)startData1)[0],
                                                      ((double*)startData1)[0]),
                "mismatch of result data");
  }
  rc = fc_deleteSeqVariable(glb_numStep, newSeqVar);
  fail_unless(rc == FC_SUCCESS, "failed to delete new var");
  free(newSeqVar);

  // --- things should fail with mixed global non-global vars
  seqVar1 = glb_dataTypeSeqVars[3];
  rc = fc_seqVarOperatorSeqVar(glb_numStep, glb_glbSeqVar, opStrings[0], seqVar1,
                               newSeqVarName, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail with mix of global/nonglobal");
  fail_unless(newSeqVar == NULL, "fail should return NULL");
  rc = fc_seqVarOperatorSeqVar(glb_numStep, seqVar1, opStrings[0], glb_glbSeqVar, 
                               newSeqVarName, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail with mix of global/nonglobal");
  fail_unless(newSeqVar == NULL, "fail should return NULL");

  // --- It is an error to have seq variables on different sequences
  seqVar1 = glb_dataTypeSeqVars[3];
  rc = fc_copySequence(glb_sequence, glb_dataset, "seq2", &seq2);
  fail_unless(rc == FC_SUCCESS, "failed to copy sequence");
  rc = fc_copySeqVariable(glb_numStep, seqVar1, glb_mesh, seq2, 
			  "var on another seq", &otherSeqVar);
  fail_unless(rc == FC_SUCCESS, "failed to create seq var on another seq");
  rc = fc_seqVarOperatorSeqVar(glb_numStep, seqVar1, opStrings[0],
			       otherSeqVar, newSeqVarName, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail with bad input var");
  fail_unless(newSeqVar == NULL,"fail should return NULL");
  free(otherSeqVar);
 
  // --- test bad inputs
  seqVar1 = seqVar2 = glb_dataTypeSeqVars[3];
  rc = fc_seqVarOperatorSeqVar(glb_numStep, NULL, opStrings[0], 
			       seqVar2, newSeqVarName, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail with bad input var");
  fail_unless(newSeqVar == NULL,"fail should return NULL");
  for (i = 0; i < glb_numStep; i++)
    badSeqVar[i] = glb_dataTypeSeqVars[3][i];
  badSeqVar[1] = badVar;
  rc = fc_seqVarOperatorSeqVar(glb_numStep, badSeqVar, opStrings[0], 
			       seqVar2, newSeqVarName, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail with bad input var");
  fail_unless(newSeqVar == NULL,"fail should return NULL");
  rc = fc_seqVarOperatorSeqVar(glb_numStep, seqVar1, "badOpString", 
			       seqVar2, newSeqVarName, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail with bad op string");
  fail_unless(newSeqVar == NULL,"fail should return NULL");
  rc = fc_seqVarOperatorSeqVar(glb_numStep, seqVar1, opStrings[0], 
			       NULL, newSeqVarName, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail with bad input var");
  fail_unless(newSeqVar == NULL,"fail should return NULL");
  rc = fc_seqVarOperatorSeqVar(glb_numStep, seqVar1, opStrings[0], 
			       badSeqVar, newSeqVarName, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail with bad input var");
  fail_unless(newSeqVar == NULL,"fail should return NULL");
  rc = fc_seqVarOperatorSeqVar(glb_numStep, seqVar1, opStrings[0], 
			       seqVar2, NULL, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail with NULL name");
  fail_unless(newSeqVar == NULL,"fail should return NULL");
  rc = fc_seqVarOperatorSeqVar(glb_numStep, seqVar1, opStrings[0], 
			       seqVar2, newSeqVarName, NULL);
  fail_unless(rc != FC_SUCCESS, "should fail if NULL new var");
}
END_TEST

// test fc_seqVarOperatorSeqConst()
START_TEST(seqVarOpSeqConst)
{
  FC_ReturnCode rc;
  int i, j, k, m, n, p;
  FC_Sequence temp_seq;
  FC_Variable *seqVar, *newSeqVar, badSeqVar[3], badVar = { 999, 999 };
  int numDataPoint, numComp, temp_numDataPoint, temp_numComp;
  char newSeqVarName[20] = "new seq var";
  FC_AssociationType assoc, temp_assoc;
  FC_MathType mathType, temp_mathType;
  FC_DataType dataType, temp_dataType;
  void* startData, *resultData;
  int numOp = 4;
  char opStrings[16][10] = { "+", "-", "*", "/" };
  int (*intOps[4])(int, int) = { addInts, subInts, multInts, divInts };
  double (*dblOps[4])(double, double) = { addDbls, subDbls, multDbls, divDbls };
  int numConst = 4*2; // four data types * whether scalar or vector
  // index by step, then by component
  char charConst[3*3] = { 'a', 'z', '$',
                          'b', 'y', '@', 
                          'c', 'x', '@' };
  int intConst[3*3] = { 20, -2,   6, 
                        21, -3,  16, 
                        22, -4, 160 };
  float fltConst[3*3] = { 3.5,  -99.99,   1.001, 
                          4.5, -100.99,  10.01, 
                          5.5, -101.99, 100.1 };
  double dblConst[3*3] = { 1.2,  0.00001,   -29.9,
                           2.2, -1.00001,  -299,
                           3.2, -2.00001, -2990 };
  int constNumComps[4*2] = { 1, 3, 1, 3, 1, 3, 1, 3 };
  void *constValues[4*2] = { &charConst, &charConst, &intConst, &intConst,
			     &fltConst, &fltConst, &dblConst, &dblConst };  
  FC_DataType constTypes[4*2] = { FC_DT_CHAR, FC_DT_CHAR, FC_DT_INT,
                                  FC_DT_INT, FC_DT_FLOAT, FC_DT_FLOAT,
                                  FC_DT_DOUBLE, FC_DT_DOUBLE };

  // --- use global variables created in the fixture

  // --- test all ops on all combinations of vars & consts
  for (i = 0; i < glb_numDataType + 1; i++) {

    if (i < glb_numDataType)
      seqVar = glb_dataTypeSeqVars[i];
    else 
      seqVar = glb_vectorSeqVar;

    for (j = 0; j < numConst; j++) {

      // get the start meta data
      rc = fc_getVariableInfo(seqVar[0], &numDataPoint, &numComp, &assoc, 
			      &mathType,  &dataType);
      fail_unless(rc == FC_SUCCESS, "abort: failed to get var info");

      // loop over ops
      for (k = 0; k < numOp; k++) {

	// do it
	rc = fc_seqVarOperatorSeqConst(glb_numStep, seqVar, opStrings[k],
                                       constNumComps[j], constValues[j],
                                       constTypes[j], newSeqVarName, &newSeqVar);
	if (dataType == FC_DT_CHAR || constTypes[j] == FC_DT_CHAR ||
	    numComp != constNumComps[j]) {
	  fail_unless(rc != FC_SUCCESS, "should fail to do char data or if "
		      "don't have same numComp");
	  fail_unless(newSeqVar == NULL, "fail should return NULL");
	  continue;
	}
	else
	  fail_unless(rc == FC_SUCCESS, "failed to do binary operation");
       
	// test result meta data
	rc = fc_getVariableInfo(newSeqVar[0], &temp_numDataPoint, &temp_numComp, 
			      &temp_assoc, &temp_mathType, &temp_dataType);
	fail_unless(rc == FC_SUCCESS, "failed to get new variable info");
	fail_unless(temp_numDataPoint == numDataPoint, 
		    "mismatch of numDataPoint");
	fail_unless(temp_numComp == numComp, "mismatch of numComp");
	fail_unless(temp_assoc == assoc, "mismatch of assoc");
	fail_unless(temp_mathType == mathType, "mismatch of mathType");
	if (dataType == FC_DT_INT && constTypes[j] == FC_DT_INT)
	  fail_unless(temp_dataType == FC_DT_INT, "mismatch of datatype");
	else
	  fail_unless(temp_dataType == FC_DT_DOUBLE, "mismatch of dataType");
	rc = fc_getSequenceFromSeqVariable(glb_numStep, newSeqVar, &temp_seq);
	fail_unless(rc == FC_SUCCESS, "abort: failed to get new sequence");
	fail_unless(FC_HANDLE_EQUIV(temp_seq, glb_sequence), "mismatch of sequence");
  
	// test result big data
	for (m = 0; m < glb_numStep; m++) {	
	  rc = fc_getVariableDataPtr(seqVar[m], &startData);
	  fail_unless(rc == FC_SUCCESS, "abort: failed to get data pointer");
	  rc = fc_getVariableDataPtr(newSeqVar[m], &resultData);
	  fail_unless(rc == FC_SUCCESS, "failed to get data ptr from new var");
	  for (n = 0; n < numDataPoint; n++) {
            for (p = 0; p < numComp; p++) {
              int var_idx = n*numComp + p;
              int const_idx = m*numComp + p;
              if (temp_dataType == FC_DT_INT) {
                fail_unless(((int*)resultData)[var_idx] == 
                            intOps[k](((int*)startData)[var_idx], 
                                      intConst[const_idx]),
                            "mismatch of reusult data");
              }
              else {
                double varStart, constStart;
                switch (dataType) {
                case FC_DT_INT:     varStart =    ((int*)startData)[var_idx];  break;
                case FC_DT_FLOAT:   varStart =  ((float*)startData)[var_idx];  break;
                case FC_DT_DOUBLE:  varStart = ((double*)startData)[var_idx];  break;
                default: ;
                }
                switch (constTypes[j]) {
                case FC_DT_INT:     constStart = intConst[const_idx];  break;
                case FC_DT_FLOAT:   constStart = fltConst[const_idx];  break;
                case FC_DT_DOUBLE:  constStart = dblConst[const_idx];  break;
                default: ;
                }
                fail_unless(FC_DBL_EQUIV(((double*)resultData)[var_idx],
                                         dblOps[k](varStart, constStart)),
                            "mismatch of result data");
              }
            }
	  }
	}
       
        
	rc = fc_deleteSeqVariable(glb_numStep, newSeqVar);
	fail_unless(rc == FC_SUCCESS, "failed to delete new seqvar");
	free(newSeqVar);
      }  
    }
  }

  // --- things should still work with global vars
  rc = fc_seqVarOperatorSeqConst(glb_numStep, glb_glbSeqVar, opStrings[0],
                                 1, dblConst, FC_DT_DOUBLE, 
                                 newSeqVarName, &newSeqVar);
  fail_unless(rc == FC_SUCCESS, "should work with global var");
  fail_unless(fc_isVariableGlobal(newSeqVar[0]), "should be global");
  for (i = 0; i < glb_numStep; i++) {
    rc = fc_getVariableDataPtr(glb_glbSeqVar[i], &startData);
    fail_unless(rc == FC_SUCCESS, "failed to get data ptr from old var");
    rc = fc_getVariableDataPtr(newSeqVar[i], &resultData);
    fail_unless(rc == FC_SUCCESS, "failed to get data ptr from new var");
    // just checking first value
    fail_unless(FC_DBL_EQUIV(((double*)resultData)[0],
                             dblOps[0](((double*)startData)[0], dblConst[i])),
                "mismatch of result data");
  }
  rc = fc_deleteSeqVariable(glb_numStep, newSeqVar);
  fail_unless(rc == FC_SUCCESS, "failed to delete new var");
  free(newSeqVar);
 
  // --- test bad inputs
  seqVar = glb_dataTypeSeqVars[3];
  rc = fc_seqVarOperatorSeqConst(glb_numStep, NULL, opStrings[0], 
                                 1, dblConst, FC_DT_DOUBLE, 
                                 newSeqVarName, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail with bad input var");
  fail_unless(newSeqVar == NULL,"fail should return NULL");
  for (i = 0; i < glb_numStep; i++)
    badSeqVar[i] = glb_dataTypeSeqVars[3][i];
  badSeqVar[1] = badVar;
  rc = fc_seqVarOperatorSeqConst(glb_numStep, badSeqVar, opStrings[0], 
                                 1, dblConst, FC_DT_DOUBLE, 
                                 newSeqVarName, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail with bad input var");
  fail_unless(newSeqVar == NULL,"fail should return NULL");
  rc = fc_seqVarOperatorSeqConst(glb_numStep, seqVar, "badOpString", 
                                 1, dblConst, FC_DT_DOUBLE,
                                 newSeqVarName, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail with bad op string");
  fail_unless(newSeqVar == NULL,"fail should return NULL");
  rc = fc_seqVarOperatorSeqConst(glb_numStep, seqVar, opStrings[0], 
                                 0, dblConst, FC_DT_DOUBLE, 
                                 newSeqVarName, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail if numConstComp < 1");
  fail_unless(newSeqVar == NULL,"fail should return NULL");
  rc = fc_seqVarOperatorSeqConst(glb_numStep, seqVar, opStrings[0], 
                                 1, NULL, FC_DT_DOUBLE,
                                 newSeqVarName, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail if NULL constValues");
  fail_unless(newSeqVar == NULL,"fail should return NULL");
  rc = fc_seqVarOperatorSeqConst(glb_numStep, seqVar, opStrings[0], 
                                 1, dblConst, FC_DT_UNKNOWN,
                                 newSeqVarName, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail if unknown data type");
  fail_unless(newSeqVar == NULL,"fail should return NULL");
  rc = fc_seqVarOperatorSeqConst(glb_numStep, seqVar, opStrings[0], 
                                 1, dblConst, FC_DT_DOUBLE, NULL, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail with NULL name");
  fail_unless(newSeqVar == NULL,"fail should return NULL");
  rc = fc_seqVarOperatorSeqConst(glb_numStep, seqVar, opStrings[0], 
                                 1, dblConst, FC_DT_DOUBLE, newSeqVarName, NULL);
  fail_unless(rc != FC_SUCCESS, "should fail if NULL new var");
}
END_TEST

// test fc_seqConstOperatorSeqVar()
START_TEST(seqConstOpSeqVar)
{
  FC_ReturnCode rc;
  int i, j, k, m, n, p;
  FC_Sequence temp_seq;
  FC_Variable *seqVar, *newSeqVar, badSeqVar[3], badVar = { 999, 999 };
  int numDataPoint, numComp, temp_numDataPoint, temp_numComp;
  char newSeqVarName[20] = "new seq var";
  FC_AssociationType assoc, temp_assoc;
  FC_MathType mathType, temp_mathType;
  FC_DataType dataType, temp_dataType;
  void* startData, *resultData;
  int numOp = 4;
  char opStrings[16][10] = { "+", "-", "*", "/" };
  int (*intOps[4])(int, int) = { addInts, subInts, multInts, divInts };
  double (*dblOps[4])(double, double) = { addDbls, subDbls, multDbls, divDbls };
  int numConst = 4*2; // four data types * whether scalar or vector
  // index by step, then by component
  char charConst[3*3] = { 'a', 'z', '$',
                          'b', 'y', '@', 
                          'c', 'x', '@' };
  int intConst[3*3] = { 20, -2,   6, 
                        21, -3,  16, 
                        22, -4, 160 };
  float fltConst[3*3] = { 3.5,  -99.99,   1.001, 
                          4.5, -100.99,  10.01, 
                          5.5, -101.99, 100.1 };
  double dblConst[3*3] = { 1.2,  0.00001,   -29.9,
                           2.2, -1.00001,  -299,
                           3.2, -2.00001, -2990 };
  int constNumComps[4*2] = { 1, 3, 1, 3, 1, 3, 1, 3 };
  void *constValues[4*2] = { &charConst, &charConst, &intConst, &intConst,
			     &fltConst, &fltConst, &dblConst, &dblConst };  
  FC_DataType constTypes[4*2] = { FC_DT_CHAR, FC_DT_CHAR, FC_DT_INT,
                                  FC_DT_INT, FC_DT_FLOAT, FC_DT_FLOAT,
                                  FC_DT_DOUBLE, FC_DT_DOUBLE };
  
  // --- use global variables created in the fixture

  // --- test all ops on all combinations of vars & consts
  for (i = 0; i < glb_numDataType + 1; i++) {

    if (i < glb_numDataType)
      seqVar = glb_dataTypeSeqVars[i];
    else 
      seqVar = glb_vectorSeqVar;

    for (j = 0; j < numConst; j++) {

      // get the start meta data
      rc = fc_getVariableInfo(seqVar[0], &numDataPoint, &numComp, &assoc, 
			      &mathType,  &dataType);
      fail_unless(rc == FC_SUCCESS, "abort: failed to get var info");

      // loop over ops
      for (k = 0; k < numOp; k++) {

	// do it
	rc = fc_seqConstOperatorSeqVar(glb_numStep, constNumComps[j], 
                                       constValues[j], constTypes[j], 
                                       opStrings[k], seqVar, 
                                       newSeqVarName, &newSeqVar);
	if (dataType == FC_DT_CHAR || constTypes[j] == FC_DT_CHAR ||
	    numComp != constNumComps[j]) {
	  fail_unless(rc != FC_SUCCESS, "should fail to do char data or if "
		      "don't have same numComp");
	  fail_unless(newSeqVar == NULL, "fail should return NULL");
	  continue;
	}
	else
	  fail_unless(rc == FC_SUCCESS, "failed to do binary operation");
       
	// test result meta data
	rc = fc_getVariableInfo(newSeqVar[0], &temp_numDataPoint, &temp_numComp, 
			      &temp_assoc, &temp_mathType, &temp_dataType);
	fail_unless(rc == FC_SUCCESS, "failed to get new variable info");
	fail_unless(temp_numDataPoint == numDataPoint, 
		    "mismatch of numDataPoint");
	fail_unless(temp_numComp == numComp, "mismatch of numComp");
	fail_unless(temp_assoc == assoc, "mismatch of assoc");
	fail_unless(temp_mathType == mathType, "mismatch of mathType");
	if (dataType == FC_DT_INT && constTypes[j] == FC_DT_INT)
	  fail_unless(temp_dataType == FC_DT_INT, "mismatch of datatype");
	else
	  fail_unless(temp_dataType == FC_DT_DOUBLE, "mismatch of dataType");
	rc = fc_getSequenceFromSeqVariable(glb_numStep, newSeqVar, &temp_seq);
	fail_unless(rc == FC_SUCCESS, "abort: failed to get new sequence");
	fail_unless(FC_HANDLE_EQUIV(temp_seq, glb_sequence), "mismatch of sequence");
  
	// test result big data
	for (m = 0; m < glb_numStep; m++) {	
	  rc = fc_getVariableDataPtr(seqVar[m], &startData);
	  fail_unless(rc == FC_SUCCESS, "abort: failed to get data pointer");
	  rc = fc_getVariableDataPtr(newSeqVar[m], &resultData);
	  fail_unless(rc == FC_SUCCESS, "failed to get data ptr from new var");
	  for (n = 0; n < numDataPoint; n++) {
            for (p = 0; p < numComp; p++) {
              int var_idx = n*numComp + p;
              int const_idx = m*numComp + p;
              if (temp_dataType == FC_DT_INT) {
                fail_unless(((int*)resultData)[var_idx] == 
                            intOps[k](intConst[const_idx],
                                      ((int*)startData)[var_idx]),
                            "mismatch of reusult data");
              }
              else {
                double varStart, constStart;
                switch (dataType) {
                case FC_DT_INT:     varStart =    ((int*)startData)[var_idx];  break;
                case FC_DT_FLOAT:   varStart =  ((float*)startData)[var_idx];  break;
                case FC_DT_DOUBLE:  varStart = ((double*)startData)[var_idx];  break;
                default: ;
                }
                switch (constTypes[j]) {
                case FC_DT_INT:     constStart = intConst[const_idx];  break;
                case FC_DT_FLOAT:   constStart = fltConst[const_idx];  break;
                case FC_DT_DOUBLE:  constStart = dblConst[const_idx];  break;
                default: ;
                }
                fail_unless(FC_DBL_EQUIV(((double*)resultData)[var_idx],
                                         dblOps[k](constStart, varStart)),
                            "mismatch of result data");
              }
            }
	  }
	}
       
        
	rc = fc_deleteSeqVariable(glb_numStep, newSeqVar);
	fail_unless(rc == FC_SUCCESS, "failed to delete new seqvar");
	free(newSeqVar);
      }  
    }
  }

  // --- things should still work with global vars
  rc = fc_seqConstOperatorSeqVar(glb_numStep, 1, dblConst, FC_DT_DOUBLE, 
                                 opStrings[0], glb_glbSeqVar, 
                                 newSeqVarName, &newSeqVar);
  fail_unless(rc == FC_SUCCESS, "should work with global var");
  fail_unless(fc_isVariableGlobal(newSeqVar[0]), "should be global");
  for (i = 0; i < glb_numStep; i++) {
    rc = fc_getVariableDataPtr(glb_glbSeqVar[i], &startData);
    fail_unless(rc == FC_SUCCESS, "failed to get data ptr from old var");
    rc = fc_getVariableDataPtr(newSeqVar[i], &resultData);
    fail_unless(rc == FC_SUCCESS, "failed to get data ptr from new var");
    // just checking first value
    fail_unless(FC_DBL_EQUIV(((double*)resultData)[0],
                             dblOps[0](dblConst[i], ((double*)startData)[0])),
                "mismatch of result data");
  }
  rc = fc_deleteSeqVariable(glb_numStep, newSeqVar);
  fail_unless(rc == FC_SUCCESS, "failed to delete new var");
  free(newSeqVar);

  // --- test bad inputs
  seqVar = glb_dataTypeSeqVars[3];
  rc = fc_seqConstOperatorSeqVar(glb_numStep, 0, dblConst, FC_DT_DOUBLE, 
                                 opStrings[0], seqVar, 
                                 newSeqVarName, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail if numConstComp < 1");
  fail_unless(newSeqVar == NULL,"fail should return NULL");
  rc = fc_seqConstOperatorSeqVar(glb_numStep, 1, NULL, FC_DT_DOUBLE, 
                                 opStrings[0], seqVar, 
                                 newSeqVarName, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail if NULL constValues");
  fail_unless(newSeqVar == NULL,"fail should return NULL");
  rc = fc_seqConstOperatorSeqVar(glb_numStep, 1, dblConst, FC_DT_UNKNOWN,
                                 opStrings[0], seqVar, 
                                 newSeqVarName, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail if unknown data type");
  fail_unless(newSeqVar == NULL,"fail should return NULL");
  rc = fc_seqConstOperatorSeqVar(glb_numStep, 1, dblConst, FC_DT_DOUBLE,
                                 "badOpString", seqVar, 
                                 newSeqVarName, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail with bad op string");
  fail_unless(newSeqVar == NULL,"fail should return NULL");
  rc = fc_seqConstOperatorSeqVar(glb_numStep, 1, dblConst, FC_DT_DOUBLE, 
                                 opStrings[0], NULL, 
                                 newSeqVarName, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail with bad input var");
  fail_unless(newSeqVar == NULL,"fail should return NULL");
  for (i = 0; i < glb_numStep; i++)
    badSeqVar[i] = glb_dataTypeSeqVars[3][i];
  badSeqVar[1] = badVar;
  rc = fc_seqConstOperatorSeqVar(glb_numStep, 1, dblConst, FC_DT_DOUBLE, 
                                 opStrings[0], badSeqVar, 
                                 newSeqVarName, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail with bad input var");
  fail_unless(newSeqVar == NULL,"fail should return NULL");
  rc = fc_seqConstOperatorSeqVar(glb_numStep,  1, dblConst, FC_DT_DOUBLE,
                                 opStrings[0], seqVar, NULL, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail with NULL name");
  fail_unless(newSeqVar == NULL,"fail should return NULL");
  rc = fc_seqConstOperatorSeqVar(glb_numStep, 1, dblConst, FC_DT_DOUBLE, 
                                 opStrings[0], seqVar, newSeqVarName, NULL);
  fail_unless(rc != FC_SUCCESS, "should fail if NULL new var");
}
END_TEST

// test fc_seqVarOperatorVar()
START_TEST(seqVarOpVar)
{
  FC_ReturnCode rc;
  int i, j, k, m, n;
  FC_Sequence newSeq;
  FC_Variable *seqVar1, var2, *newSeqVar, badSeqVar[3], badVar = { 999, 999 };
  int numDataPoint1, numComp1, numDataPoint2, numComp2;
  int temp_numDataPoint, temp_numComp;
  char newSeqVarName[20] = "new seq var";
  FC_AssociationType assoc1, assoc2, temp_assoc;
  FC_MathType mathType1, mathType2, temp_mathType;
  FC_DataType dataType1, dataType2, temp_dataType;
  void* startData1, *startData2, *resultData;
  int numOp = 4;
  char opStrings[16][10] = { "+", "-", "*", "/" };
  int (*intOps[4])(int, int) = { addInts, subInts, multInts, divInts };
  double (*dblOps[4])(double, double) = { addDbls, subDbls, multDbls, divDbls };

  // --- use global variables created in the fixture

  // --- test all ops on all combinations of all var types
  for (i = 0; i < glb_numDataType + 2; i++) {

    if (i < glb_numDataType)
      seqVar1 = glb_dataTypeSeqVars[i];
    else if (i == glb_numDataType)
      seqVar1 = glb_vectorSeqVar;
    else 
      seqVar1 = glb_vectorSeqVar;

    for (j = 0; j < glb_numDataType + 2; j++) {
      if (j < glb_numDataType)
	var2 = glb_dataTypeSeqVars[j][0];
      else if (j == glb_numDataType)
	var2 = glb_vectorSeqVar[0];
      else
	var2 = glb_elemSeqVar[0];

      // get the start meta data
      rc = fc_getVariableInfo(seqVar1[0], &numDataPoint1, &numComp1, &assoc1, 
			      &mathType1,  &dataType1);
      fail_unless(rc == FC_SUCCESS, "abort: failed to get var info");
      rc = fc_getVariableInfo(var2, &numDataPoint2, &numComp2, &assoc2, 
			      &mathType2, &dataType2);
      fail_unless(rc == FC_SUCCESS, "abort: failed to get var info");

      // loop over ops
      for (k = 0; k < numOp; k++) {

	// do it
	rc = fc_seqVarOperatorVar(glb_numStep, seqVar1, opStrings[k],
				  var2, newSeqVarName, &newSeqVar);
	if (dataType1 == FC_DT_CHAR || dataType2 == FC_DT_CHAR ||
	    numDataPoint1 != numDataPoint2 || numComp1 != numComp2) {
	  fail_unless(rc != FC_SUCCESS, "should fail to do char data or if "
		      "vars don't have some numDataPoint or don't have same "
		      "numComp");
	  fail_unless(newSeqVar == NULL, "fail should return NULL");
	  continue;
	}
	else
	  fail_unless(rc == FC_SUCCESS, "failed to do binary operation");
       
	// test result meta data
	rc = fc_getVariableInfo(newSeqVar[0], &temp_numDataPoint, &temp_numComp, 
			      &temp_assoc, &temp_mathType, &temp_dataType);
	fail_unless(rc == FC_SUCCESS, "failed to get new variable info");
	fail_unless(temp_numDataPoint == numDataPoint1, 
		    "mismatch of numDataPoint");
	fail_unless(temp_numComp == numComp1, "mismatch of numComp");
	fail_unless(temp_assoc == assoc1, "mismatch of assoc");
	fail_unless(temp_mathType == mathType1, "mismatch of mathType");
	if (dataType1 == FC_DT_INT && dataType2 == FC_DT_INT)
	  fail_unless(temp_dataType == FC_DT_INT, "mismatch of datatype");
	else
	  fail_unless(temp_dataType == FC_DT_DOUBLE, "mismatch of dataType");
	rc = fc_getSequenceFromSeqVariable(glb_numStep, newSeqVar, &newSeq);
	fail_unless(rc == FC_SUCCESS, "abort: failed to get new sequence");
	fail_unless(FC_HANDLE_EQUIV(newSeq, glb_sequence), "mismatch of sequence");

	// test result big data
	for (m = 0; m < glb_numStep; m++) {	
	  rc = fc_getVariableDataPtr(seqVar1[m], &startData1);
	  fail_unless(rc == FC_SUCCESS, "abort: failed to get data pointer");
	  rc = fc_getVariableDataPtr(var2, &startData2);
	  fail_unless(rc == FC_SUCCESS, "abort: failed to get data pointer");
	  rc = fc_getVariableDataPtr(newSeqVar[m], &resultData);
	  fail_unless(rc == FC_SUCCESS, "failed to get data ptr from new var");
	  for (n = 0; n < numDataPoint1*numComp1; n++) {
	    if (temp_dataType == FC_DT_INT) {
	      fail_unless(((int*)resultData)[n] == 
			  intOps[k](((int*)startData1)[n], ((int*)startData2)[n]),
			  "mismatch of reusult data");
	    }
	    else {
	      double start1, start2;
	      switch (dataType1) {
	      case FC_DT_INT:     start1 =    ((int*)startData1)[n];  break;
	      case FC_DT_FLOAT:   start1 =  ((float*)startData1)[n];  break;
	      case FC_DT_DOUBLE:  start1 = ((double*)startData1)[n];  break;
	      default: ;
	      }
	      switch (dataType2) {
	      case FC_DT_INT:     start2 =    ((int*)startData2)[n];  break;
	      case FC_DT_FLOAT:   start2 =  ((float*)startData2)[n];  break;
	      case FC_DT_DOUBLE:  start2 = ((double*)startData2)[n];  break;
	      default: ;
	      }
	      fail_unless(FC_DBL_EQUIV(((double*)resultData)[n],
				       dblOps[k](start1, start2)),
			  "mismatch of result data");
	    }
	  }
	}
	rc = fc_deleteSeqVariable(glb_numStep, newSeqVar);
	fail_unless(rc == FC_SUCCESS, "failed to delete new seqvar");
	free(newSeqVar);
      }
    }
  }

  // --- things should still work with global vars
  rc = fc_seqVarOperatorVar(glb_numStep, glb_glbSeqVar, opStrings[0], 
                            glb_glbSeqVar[0], newSeqVarName, &newSeqVar);
  fail_unless(rc == FC_SUCCESS, "should work with global var");
  fail_unless(fc_isVariableGlobal(newSeqVar[0]), "should be global");
  for (i = 0; i < glb_numStep; i++) {
    rc = fc_getVariableDataPtr(glb_glbSeqVar[i], &startData1);
    fail_unless(rc == FC_SUCCESS, "failed to get data ptr from old var");
    rc = fc_getVariableDataPtr(glb_glbSeqVar[0], &startData2);
    fail_unless(rc == FC_SUCCESS, "failed to get data ptr from old var");
    rc = fc_getVariableDataPtr(newSeqVar[i], &resultData);
    fail_unless(rc == FC_SUCCESS, "failed to get data ptr from new var");
    // just checking first value
    fail_unless(FC_DBL_EQUIV(((double*)resultData)[0],
			     dblOps[0](((double*)startData1)[0],
				       ((double*)startData2)[0])),
                "mismatch of result data");
  }
  rc = fc_deleteSeqVariable(glb_numStep, newSeqVar);
  fail_unless(rc == FC_SUCCESS, "failed to delete new var");
  free(newSeqVar);

  // --- things should fail with mixed global non-global vars
  seqVar1 = glb_dataTypeSeqVars[3];
  var2 = glb_dataTypeSeqVars[3][0];
  rc = fc_seqVarOperatorVar(glb_numStep, glb_glbSeqVar, opStrings[0], var2,
                            newSeqVarName, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail with mix of global/nonglobal");
  fail_unless(newSeqVar == NULL, "fail should return NULL");
  rc = fc_seqVarOperatorVar(glb_numStep, seqVar1, opStrings[0], glb_glbSeqVar[0], 
                            newSeqVarName, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail with mix of global/nonglobal");
  fail_unless(newSeqVar == NULL, "fail should return NULL");

  // --- test bad inputs
  seqVar1 = glb_dataTypeSeqVars[3];
  var2 = glb_dataTypeSeqVars[3][0];
  for (i = 0; i < glb_numStep; i++)
    badSeqVar[i] = glb_dataTypeSeqVars[3][i];
  badSeqVar[1] = badVar;
  rc = fc_seqVarOperatorVar(glb_numStep, NULL, opStrings[0], 
			    var2, newSeqVarName, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail with bad input var");
  fail_unless(newSeqVar == NULL,"fail should return NULL");
  rc = fc_seqVarOperatorVar(glb_numStep, badSeqVar, opStrings[0], 
			    var2, newSeqVarName, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail with bad input var");
  fail_unless(newSeqVar == NULL,"fail should return NULL");
  rc = fc_seqVarOperatorVar(glb_numStep, seqVar1, "badOpString", 
			    var2, newSeqVarName, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail with bad op string");
  fail_unless(newSeqVar == NULL,"fail should return NULL");
  rc = fc_seqVarOperatorVar(glb_numStep, seqVar1, opStrings[0], 
			    badVar, newSeqVarName, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail with bad input var");
  fail_unless(newSeqVar == NULL,"fail should return NULL");
  rc = fc_seqVarOperatorVar(glb_numStep, seqVar1, opStrings[0], 
			    var2, NULL, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail with NULL name");
  fail_unless(newSeqVar == NULL,"fail should return NULL");
  rc = fc_seqVarOperatorVar(glb_numStep, seqVar1, opStrings[0], 
			    var2, newSeqVarName, NULL);
  fail_unless(rc != FC_SUCCESS, "should fail if NULL new var");
}
END_TEST

// test fc_varOperatorSeqVar()
START_TEST(varOpSeqVar)
{
  FC_ReturnCode rc;
  int i, j, k, m, n;
  FC_Sequence newSeq;
  FC_Variable var1, *seqVar2, *newSeqVar, badSeqVar[3], badVar = { 999, 999 };
  int numDataPoint1, numComp1, numDataPoint2, numComp2;
  int temp_numDataPoint, temp_numComp;
  char newSeqVarName[20] = "new seq var";
  FC_AssociationType assoc1, assoc2, temp_assoc;
  FC_MathType mathType1, mathType2, temp_mathType;
  FC_DataType dataType1, dataType2, temp_dataType;
  void* startData1, *startData2, *resultData;
  int numOp = 4;
  char opStrings[16][10] = { "+", "-", "*", "/" };
  int (*intOps[4])(int, int) = { addInts, subInts, multInts, divInts };
  double (*dblOps[4])(double, double) = { addDbls, subDbls, multDbls, divDbls };

  // --- use global variables created in the fixture

  // --- test all ops on all combinations of all var types
  for (i = 0; i < glb_numDataType + 2; i++) {

    if (i < glb_numDataType)
      var1 = glb_dataTypeSeqVars[i][0];
    else if (i == glb_numDataType)
      var1 = glb_vectorSeqVar[0];
    else 
      var1 = glb_vectorSeqVar[0];

    for (j = 0; j < glb_numDataType + 2; j++) {
      if (j < glb_numDataType)
	seqVar2 = glb_dataTypeSeqVars[j];
      else if (j == glb_numDataType)
	seqVar2 = glb_vectorSeqVar;
      else
	seqVar2 = glb_elemSeqVar;

      // get the start meta data
      rc = fc_getVariableInfo(var1, &numDataPoint1, &numComp1, &assoc1, 
			      &mathType1,  &dataType1);
      fail_unless(rc == FC_SUCCESS, "abort: failed to get var info");
      rc = fc_getVariableInfo(seqVar2[0], &numDataPoint2, &numComp2, &assoc2, 
			      &mathType2, &dataType2);
      fail_unless(rc == FC_SUCCESS, "abort: failed to get var info");

      // loop over ops
      for (k = 0; k < numOp; k++) {

	// do it
	rc = fc_varOperatorSeqVar(var1, opStrings[k], glb_numStep, 
				  seqVar2, newSeqVarName, &newSeqVar);
	if (dataType1 == FC_DT_CHAR || dataType2 == FC_DT_CHAR ||
	    numDataPoint1 != numDataPoint2 || numComp1 != numComp2) {
	  fail_unless(rc != FC_SUCCESS, "should fail to do char data or if "
		      "vars don't have some numDataPoint or don't have same "
		      "numComp");
	  fail_unless(newSeqVar == NULL, "fail should return NULL");
	  continue;
	}
	else
	  fail_unless(rc == FC_SUCCESS, "failed to do binary operation");
       
	// test result meta data
	rc = fc_getVariableInfo(newSeqVar[0], &temp_numDataPoint, &temp_numComp, 
			      &temp_assoc, &temp_mathType, &temp_dataType);
	fail_unless(rc == FC_SUCCESS, "failed to get new variable info");
	fail_unless(temp_numDataPoint == numDataPoint1, 
		    "mismatch of numDataPoint");
	fail_unless(temp_numComp == numComp1, "mismatch of numComp");
	fail_unless(temp_assoc == assoc1, "mismatch of assoc");
	fail_unless(temp_mathType == mathType1, "mismatch of mathType");
	if (dataType1 == FC_DT_INT && dataType2 == FC_DT_INT)
	  fail_unless(temp_dataType == FC_DT_INT, "mismatch of datatype");
	else
	  fail_unless(temp_dataType == FC_DT_DOUBLE, "mismatch of dataType");
	rc = fc_getSequenceFromSeqVariable(glb_numStep, newSeqVar, &newSeq);
	fail_unless(rc == FC_SUCCESS, "abort: failed to get new sequence");
	fail_unless(FC_HANDLE_EQUIV(newSeq, glb_sequence), "mismatch of sequence");

	// test result big data
	for (m = 0; m < glb_numStep; m++) {	
	  rc = fc_getVariableDataPtr(var1, &startData1);
	  fail_unless(rc == FC_SUCCESS, "abort: failed to get data pointer");
	  rc = fc_getVariableDataPtr(seqVar2[m], &startData2);
	  fail_unless(rc == FC_SUCCESS, "abort: failed to get data pointer");
	  rc = fc_getVariableDataPtr(newSeqVar[m], &resultData);
	  fail_unless(rc == FC_SUCCESS, "failed to get data ptr from new var");
	  for (n = 0; n < numDataPoint1*numComp1; n++) {
	    if (temp_dataType == FC_DT_INT) {
	      fail_unless(((int*)resultData)[n] == 
			  intOps[k](((int*)startData1)[n], ((int*)startData2)[n]),
			  "mismatch of reusult data");
	    }
	    else {
	      double start1, start2;
	      switch (dataType1) {
	      case FC_DT_INT:     start1 =    ((int*)startData1)[n];  break;
	      case FC_DT_FLOAT:   start1 =  ((float*)startData1)[n];  break;
	      case FC_DT_DOUBLE:  start1 = ((double*)startData1)[n];  break;
	      default: ;
	      }
	      switch (dataType2) {
	      case FC_DT_INT:     start2 =    ((int*)startData2)[n];  break;
	      case FC_DT_FLOAT:   start2 =  ((float*)startData2)[n];  break;
	      case FC_DT_DOUBLE:  start2 = ((double*)startData2)[n];  break;
	      default: ;
	      }
	      fail_unless(FC_DBL_EQUIV(((double*)resultData)[n],
				       dblOps[k](start1, start2)),
			  "mismatch of result data");
	    }
	  }
	}
	rc = fc_deleteSeqVariable(glb_numStep, newSeqVar);
	fail_unless(rc == FC_SUCCESS, "failed to delete new seqvar");
	free(newSeqVar);
      }
    }
  }

  // --- things should still work with global vars
  rc = fc_varOperatorSeqVar(glb_glbSeqVar[0], opStrings[0], glb_numStep, 
                            glb_glbSeqVar, newSeqVarName, &newSeqVar);
  fail_unless(rc == FC_SUCCESS, "should work with global var");
  fail_unless(fc_isVariableGlobal(newSeqVar[0]), "should be global");
  for (i = 0; i < glb_numStep; i++) {
    rc = fc_getVariableDataPtr(glb_glbSeqVar[0], &startData1);
    fail_unless(rc == FC_SUCCESS, "failed to get data ptr from old var");
    rc = fc_getVariableDataPtr(glb_glbSeqVar[i], &startData2);
    fail_unless(rc == FC_SUCCESS, "failed to get data ptr from old var");
    rc = fc_getVariableDataPtr(newSeqVar[i], &resultData);
    fail_unless(rc == FC_SUCCESS, "failed to get data ptr from new var");
    // just checking first value
    fail_unless(FC_DBL_EQUIV(((double*)resultData)[0],
			     dblOps[0](((double*)startData1)[0],
				       ((double*)startData2)[0])),
                "mismatch of result data");
  }
  rc = fc_deleteSeqVariable(glb_numStep, newSeqVar);
  fail_unless(rc == FC_SUCCESS, "failed to delete new var");
  free(newSeqVar);

  // --- things should fail with mixed global non-global vars
  var1 = glb_dataTypeSeqVars[3][0];
  seqVar2 = glb_dataTypeSeqVars[3];
  rc = fc_varOperatorSeqVar(var1, opStrings[0], glb_numStep, glb_glbSeqVar,
                            newSeqVarName, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail with mix of global/nonglobal");
  fail_unless(newSeqVar == NULL, "fail should return NULL");
  rc = fc_varOperatorSeqVar(glb_glbSeqVar[0], opStrings[0], glb_numStep, seqVar2, 
                            newSeqVarName, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail with mix of global/nonglobal");
  fail_unless(newSeqVar == NULL, "fail should return NULL");

  // --- test bad inputs
  var1 = glb_dataTypeSeqVars[3][0];
  seqVar2 = glb_dataTypeSeqVars[3];
  for (i = 0; i < glb_numStep; i++)
    badSeqVar[i] = glb_dataTypeSeqVars[3][i];
  badSeqVar[1] = badVar;
  rc = fc_varOperatorSeqVar(badVar, opStrings[0], glb_numStep, 
			    seqVar2, newSeqVarName, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail with bad input var");
  fail_unless(newSeqVar == NULL,"fail should return NULL");
  rc = fc_varOperatorSeqVar(var1, "badOpString", glb_numStep,  
			    seqVar2, newSeqVarName, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail with bad op string");
  fail_unless(newSeqVar == NULL,"fail should return NULL");
  rc = fc_varOperatorSeqVar(var1, opStrings[0], glb_numStep, 
			    NULL, newSeqVarName, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail with bad input var");
  fail_unless(newSeqVar == NULL,"fail should return NULL");
  rc = fc_varOperatorSeqVar(var1, opStrings[0], glb_numStep, 
			    badSeqVar, newSeqVarName, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail with bad input var");
  fail_unless(newSeqVar == NULL,"fail should return NULL");
  rc = fc_varOperatorSeqVar(var1, opStrings[0], glb_numStep, 
			    seqVar2, NULL, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail with NULL name");
  fail_unless(newSeqVar == NULL,"fail should return NULL");
  rc = fc_varOperatorSeqVar(var1, opStrings[0], glb_numStep, 
			    seqVar2, newSeqVarName, NULL);
  fail_unless(rc != FC_SUCCESS, "should fail if NULL new var");
}
END_TEST

// test fc_seqVarOperatorConst()
START_TEST(seqVarOpConst)
{
  FC_ReturnCode rc;
  int i, j, k, m, n, p;
  FC_Sequence temp_seq;
  FC_Variable *seqVar, *newSeqVar, badSeqVar[3], badVar = { 999, 999 };
  int numDataPoint, numComp, temp_numDataPoint, temp_numComp;
  char newSeqVarName[20] = "new seq var";
  FC_AssociationType assoc, temp_assoc;
  FC_MathType mathType, temp_mathType;
  FC_DataType dataType, temp_dataType;
  void* startData, *resultData;
  int numOp = 4;
  char opStrings[16][10] = { "+", "-", "*", "/" };
  int (*intOps[4])(int, int) = { addInts, subInts, multInts, divInts };
  double (*dblOps[4])(double, double) = { addDbls, subDbls, multDbls, divDbls };
  int numConst = 4*2; // four data types * whether scalar or vector
  // index by step, then by component
  char charConst[3] = { 'a', 'z', '$' };
  int intConst[3] = { 20, -2, 6 };
  float fltConst[3] = { 3.5, -99.99, 1.001 };
  double dblConst[3] = { 1.2, 0.00001, -29.9 };
  int constNumComps[4*2] = { 1, 3, 1, 3, 1, 3, 1, 3 };
  void *constValues[4*2] = { &charConst, &charConst, &intConst, &intConst,
			     &fltConst, &fltConst, &dblConst, &dblConst };  
  FC_DataType constTypes[4*2] = { FC_DT_CHAR, FC_DT_CHAR, FC_DT_INT,
                                  FC_DT_INT, FC_DT_FLOAT, FC_DT_FLOAT,
                                  FC_DT_DOUBLE, FC_DT_DOUBLE };

  // --- use global variables created in the fixture

  // --- test all ops on all combinations of vars & consts
  for (i = 0; i < glb_numDataType + 1; i++) {

    if (i < glb_numDataType)
      seqVar = glb_dataTypeSeqVars[i];
    else 
      seqVar = glb_vectorSeqVar;

    for (j = 0; j < numConst; j++) {

      // get the start meta data
      rc = fc_getVariableInfo(seqVar[0], &numDataPoint, &numComp, &assoc, 
			      &mathType,  &dataType);
      fail_unless(rc == FC_SUCCESS, "abort: failed to get var info");

      // loop over ops
      for (k = 0; k < numOp; k++) {

	// do it
	rc = fc_seqVarOperatorConst(glb_numStep, seqVar, opStrings[k],
                                    constNumComps[j], constValues[j],
                                    constTypes[j], newSeqVarName, &newSeqVar);
	if (dataType == FC_DT_CHAR || constTypes[j] == FC_DT_CHAR ||
	    numComp != constNumComps[j]) {
	  fail_unless(rc != FC_SUCCESS, "should fail to do char data or if "
		      "don't have same numComp");
	  fail_unless(newSeqVar == NULL, "fail should return NULL");
	  continue;
	}
	else
	  fail_unless(rc == FC_SUCCESS, "failed to do binary operation");
       
	// test result meta data
	rc = fc_getVariableInfo(newSeqVar[0], &temp_numDataPoint, &temp_numComp, 
                                &temp_assoc, &temp_mathType, &temp_dataType);
	fail_unless(rc == FC_SUCCESS, "failed to get new variable info");
	fail_unless(temp_numDataPoint == numDataPoint, 
		    "mismatch of numDataPoint");
	fail_unless(temp_numComp == numComp, "mismatch of numComp");
	fail_unless(temp_assoc == assoc, "mismatch of assoc");
	fail_unless(temp_mathType == mathType, "mismatch of mathType");
	if (dataType == FC_DT_INT && constTypes[j] == FC_DT_INT)
	  fail_unless(temp_dataType == FC_DT_INT, "mismatch of datatype");
	else
	  fail_unless(temp_dataType == FC_DT_DOUBLE, "mismatch of dataType");
	rc = fc_getSequenceFromSeqVariable(glb_numStep, newSeqVar, &temp_seq);
	fail_unless(rc == FC_SUCCESS, "abort: failed to get new sequence");
	fail_unless(FC_HANDLE_EQUIV(temp_seq, glb_sequence), "mismatch of sequence");
  
	// test result big data
	for (m = 0; m < glb_numStep; m++) {	
	  rc = fc_getVariableDataPtr(seqVar[m], &startData);
	  fail_unless(rc == FC_SUCCESS, "abort: failed to get data pointer");
	  rc = fc_getVariableDataPtr(newSeqVar[m], &resultData);
	  fail_unless(rc == FC_SUCCESS, "failed to get data ptr from new var");
	  for (n = 0; n < numDataPoint; n++) {
            for (p = 0; p < numComp; p++) {
              int var_idx = n*numComp + p;
              if (temp_dataType == FC_DT_INT) {
                fail_unless(((int*)resultData)[var_idx] == 
                            intOps[k](((int*)startData)[var_idx], intConst[p]),
                            "mismatch of reusult data");
              }
              else {
                double varStart, constStart;
                switch (dataType) {
                case FC_DT_INT:     varStart =    ((int*)startData)[var_idx];  break;
                case FC_DT_FLOAT:   varStart =  ((float*)startData)[var_idx];  break;
                case FC_DT_DOUBLE:  varStart = ((double*)startData)[var_idx];  break;
                default: ;
                }
                switch (constTypes[j]) {
                case FC_DT_INT:     constStart = intConst[p];  break;
                case FC_DT_FLOAT:   constStart = fltConst[p];  break;
                case FC_DT_DOUBLE:  constStart = dblConst[p];  break;
                default: ;
                }
                fail_unless(FC_DBL_EQUIV(((double*)resultData)[var_idx],
                                         dblOps[k](varStart, constStart)),
                            "mismatch of result data");
              }
            }
	  }
	}
       
        
	rc = fc_deleteSeqVariable(glb_numStep, newSeqVar);
	fail_unless(rc == FC_SUCCESS, "failed to delete new seqvar");
	free(newSeqVar);
      }  
    }
  }
 
  // --- things should still work with global vars
  rc = fc_seqVarOperatorConst(glb_numStep, glb_glbSeqVar, opStrings[0],
                                 1, dblConst, FC_DT_DOUBLE, 
                                 newSeqVarName, &newSeqVar);
  fail_unless(rc == FC_SUCCESS, "should work with global var");
  fail_unless(fc_isVariableGlobal(newSeqVar[0]), "should be global");
  for (i = 0; i < glb_numStep; i++) {
    rc = fc_getVariableDataPtr(glb_glbSeqVar[i], &startData);
    fail_unless(rc == FC_SUCCESS, "failed to get data ptr from old var");
    rc = fc_getVariableDataPtr(newSeqVar[i], &resultData);
    fail_unless(rc == FC_SUCCESS, "failed to get data ptr from new var");
    // just checking first value
    fail_unless(FC_DBL_EQUIV(((double*)resultData)[0],
                             dblOps[0](((double*)startData)[0], dblConst[0])),
                "mismatch of result data");
  }
  rc = fc_deleteSeqVariable(glb_numStep, newSeqVar);
  fail_unless(rc == FC_SUCCESS, "failed to delete new var");
  free(newSeqVar);
 
  // --- test bad inputs
  seqVar = glb_dataTypeSeqVars[3];
  rc = fc_seqVarOperatorConst(glb_numStep, NULL, opStrings[0], 
                              1, dblConst, FC_DT_DOUBLE, 
                              newSeqVarName, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail with bad input var");
  fail_unless(newSeqVar == NULL,"fail should return NULL");
  for (i = 0; i < glb_numStep; i++)
    badSeqVar[i] = glb_dataTypeSeqVars[3][i];
  badSeqVar[1] = badVar;
  rc = fc_seqVarOperatorConst(glb_numStep, badSeqVar, opStrings[0], 
                              1, dblConst, FC_DT_DOUBLE, 
                              newSeqVarName, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail with bad input var");
  fail_unless(newSeqVar == NULL,"fail should return NULL");
  rc = fc_seqVarOperatorConst(glb_numStep, seqVar, "badOpString", 
                              1, dblConst, FC_DT_DOUBLE,
                              newSeqVarName, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail with bad op string");
  fail_unless(newSeqVar == NULL,"fail should return NULL");
  rc = fc_seqVarOperatorConst(glb_numStep, seqVar, opStrings[0], 
                              0, dblConst, FC_DT_DOUBLE, 
                              newSeqVarName, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail if numConstComp < 1");
  fail_unless(newSeqVar == NULL,"fail should return NULL");
  rc = fc_seqVarOperatorConst(glb_numStep, seqVar, opStrings[0], 
                              1, NULL, FC_DT_DOUBLE,
                              newSeqVarName, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail if NULL constValues");
  fail_unless(newSeqVar == NULL,"fail should return NULL");
  rc = fc_seqVarOperatorConst(glb_numStep, seqVar, opStrings[0], 
                              1, dblConst, FC_DT_UNKNOWN,
                              newSeqVarName, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail if unknown data type");
  fail_unless(newSeqVar == NULL,"fail should return NULL");
  rc = fc_seqVarOperatorConst(glb_numStep, seqVar, opStrings[0], 
                              1, dblConst, FC_DT_DOUBLE, NULL, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail with NULL name");
  fail_unless(newSeqVar == NULL,"fail should return NULL");
  rc = fc_seqVarOperatorConst(glb_numStep, seqVar, opStrings[0], 
                              1, dblConst, FC_DT_DOUBLE, newSeqVarName, NULL);
  fail_unless(rc != FC_SUCCESS, "should fail if NULL new var");
}
END_TEST

// test fc_constOperatorSeqVar()
START_TEST(constOpSeqVar)
{
  FC_ReturnCode rc;
  int i, j, k, m, n, p;
  FC_Sequence temp_seq;
  FC_Variable *seqVar, *newSeqVar, badSeqVar[3], badVar = { 999, 999 };
  int numDataPoint, numComp, temp_numDataPoint, temp_numComp;
  char newSeqVarName[20] = "new seq var";
  FC_AssociationType assoc, temp_assoc;
  FC_MathType mathType, temp_mathType;
  FC_DataType dataType, temp_dataType;
  void* startData, *resultData;
  int numOp = 4;
  char opStrings[16][10] = { "+", "-", "*", "/" };
  int (*intOps[4])(int, int) = { addInts, subInts, multInts, divInts };
  double (*dblOps[4])(double, double) = { addDbls, subDbls, multDbls, divDbls };
  int numConst = 4*2; // four data types * whether scalar or vector
  // index by step, then by component
  char charConst[3] = { 'a', 'z', '$' };
  int intConst[3] = { 20, -2, 6 };
  float fltConst[3] = { 3.5, -99.99, 1.001 };
  double dblConst[3] = { 1.2, 0.00001, -29.9 };
  int constNumComps[4*2] = { 1, 3, 1, 3, 1, 3, 1, 3 };
  void *constValues[4*2] = { &charConst, &charConst, &intConst, &intConst,
			     &fltConst, &fltConst, &dblConst, &dblConst };  
  FC_DataType constTypes[4*2] = { FC_DT_CHAR, FC_DT_CHAR, FC_DT_INT,
                                  FC_DT_INT, FC_DT_FLOAT, FC_DT_FLOAT,
                                  FC_DT_DOUBLE, FC_DT_DOUBLE };

  // --- use global variables created in the fixture

  // --- test all ops on all combinations of vars & consts
  for (i = 0; i < glb_numDataType + 1; i++) {

    if (i < glb_numDataType)
      seqVar = glb_dataTypeSeqVars[i];
    else 
      seqVar = glb_vectorSeqVar;

    for (j = 0; j < numConst; j++) {

      // get the start meta data
      rc = fc_getVariableInfo(seqVar[0], &numDataPoint, &numComp, &assoc, 
			      &mathType,  &dataType);
      fail_unless(rc == FC_SUCCESS, "abort: failed to get var info");

      // loop over ops
      for (k = 0; k < numOp; k++) {

	// do it
	rc = fc_constOperatorSeqVar(constNumComps[j], constValues[j],
                                    constTypes[j], opStrings[k], 
                                    glb_numStep, seqVar, 
                                    newSeqVarName, &newSeqVar);
	if (dataType == FC_DT_CHAR || constTypes[j] == FC_DT_CHAR ||
	    numComp != constNumComps[j]) {
	  fail_unless(rc != FC_SUCCESS, "should fail to do char data or if "
		      "don't have same numComp");
	  fail_unless(newSeqVar == NULL, "fail should return NULL");
	  continue;
	}
	else
	  fail_unless(rc == FC_SUCCESS, "failed to do binary operation");
       
	// test result meta data
	rc = fc_getVariableInfo(newSeqVar[0], &temp_numDataPoint, &temp_numComp, 
			      &temp_assoc, &temp_mathType, &temp_dataType);
	fail_unless(rc == FC_SUCCESS, "failed to get new variable info");
	fail_unless(temp_numDataPoint == numDataPoint, 
		    "mismatch of numDataPoint");
	fail_unless(temp_numComp == numComp, "mismatch of numComp");
	fail_unless(temp_assoc == assoc, "mismatch of assoc");
	fail_unless(temp_mathType == mathType, "mismatch of mathType");
	if (dataType == FC_DT_INT && constTypes[j] == FC_DT_INT)
	  fail_unless(temp_dataType == FC_DT_INT, "mismatch of datatype");
	else
	  fail_unless(temp_dataType == FC_DT_DOUBLE, "mismatch of dataType");
	rc = fc_getSequenceFromSeqVariable(glb_numStep, newSeqVar, &temp_seq);
	fail_unless(rc == FC_SUCCESS, "abort: failed to get new sequence");
	fail_unless(FC_HANDLE_EQUIV(temp_seq, glb_sequence), "mismatch of sequence");
  
	// test result big data
	for (m = 0; m < glb_numStep; m++) {	
	  rc = fc_getVariableDataPtr(seqVar[m], &startData);
	  fail_unless(rc == FC_SUCCESS, "abort: failed to get data pointer");
	  rc = fc_getVariableDataPtr(newSeqVar[m], &resultData);
	  fail_unless(rc == FC_SUCCESS, "failed to get data ptr from new var");
	  for (n = 0; n < numDataPoint; n++) {
            for (p = 0; p < numComp; p++) {
              int var_idx = n*numComp + p;
              if (temp_dataType == FC_DT_INT) {
                fail_unless(((int*)resultData)[var_idx] == 
                            intOps[k](intConst[p], ((int*)startData)[var_idx]),
                            "mismatch of reusult data");
              }
              else {
                double varStart, constStart;
                switch (dataType) {
                case FC_DT_INT:     varStart =    ((int*)startData)[var_idx];  break;
                case FC_DT_FLOAT:   varStart =  ((float*)startData)[var_idx];  break;
                case FC_DT_DOUBLE:  varStart = ((double*)startData)[var_idx];  break;
                default: ;
                }
                switch (constTypes[j]) {
                case FC_DT_INT:     constStart = intConst[p];  break;
                case FC_DT_FLOAT:   constStart = fltConst[p];  break;
                case FC_DT_DOUBLE:  constStart = dblConst[p];  break;
                default: ;
                }
                fail_unless(FC_DBL_EQUIV(((double*)resultData)[var_idx],
                                         dblOps[k](constStart, varStart)),
                            "mismatch of result data");
              }
            }
	  }
	}
       
        
	rc = fc_deleteSeqVariable(glb_numStep, newSeqVar);
	fail_unless(rc == FC_SUCCESS, "failed to delete new seqvar");
	free(newSeqVar);
      }  
    }
  }
 
  // --- things should still work with global vars
  rc = fc_constOperatorSeqVar(1, dblConst, FC_DT_DOUBLE, opStrings[0],
                              glb_numStep, glb_glbSeqVar,    
                              newSeqVarName, &newSeqVar);
  fail_unless(rc == FC_SUCCESS, "should work with global var");
  fail_unless(fc_isVariableGlobal(newSeqVar[0]), "should be global");
  for (i = 0; i < glb_numStep; i++) {
    rc = fc_getVariableDataPtr(glb_glbSeqVar[i], &startData);
    fail_unless(rc == FC_SUCCESS, "failed to get data ptr from old var");
    rc = fc_getVariableDataPtr(newSeqVar[i], &resultData);
    fail_unless(rc == FC_SUCCESS, "failed to get data ptr from new var");
    // just checking first value
    fail_unless(FC_DBL_EQUIV(((double*)resultData)[0],
                             dblOps[0](((double*)startData)[0], dblConst[0])),
                "mismatch of result data");
  }
  rc = fc_deleteSeqVariable(glb_numStep, newSeqVar);
  fail_unless(rc == FC_SUCCESS, "failed to delete new var");
  free(newSeqVar);

  // --- test bad inputs
  seqVar = glb_dataTypeSeqVars[3];
  rc = fc_constOperatorSeqVar(1, dblConst, FC_DT_DOUBLE, opStrings[0], 
                              glb_numStep, NULL,    
                              newSeqVarName, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail with bad input var");
  fail_unless(newSeqVar == NULL,"fail should return NULL");
  for (i = 0; i < glb_numStep; i++)
    badSeqVar[i] = glb_dataTypeSeqVars[3][i];
  badSeqVar[1] = badVar;
  rc = fc_constOperatorSeqVar(1, dblConst, FC_DT_DOUBLE,  opStrings[0], 
                              glb_numStep, badSeqVar,   
                              newSeqVarName, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail with bad input var");
  fail_unless(newSeqVar == NULL,"fail should return NULL");
  rc = fc_constOperatorSeqVar(1, dblConst, FC_DT_DOUBLE, "badOpString", 
                              glb_numStep, seqVar,  
                              newSeqVarName, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail with bad op string");
  fail_unless(newSeqVar == NULL,"fail should return NULL");
  rc = fc_constOperatorSeqVar( 0, dblConst, FC_DT_DOUBLE,  opStrings[0],  
                               glb_numStep, seqVar,
                               newSeqVarName, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail if numConstComp < 1");
  fail_unless(newSeqVar == NULL,"fail should return NULL");
  rc = fc_constOperatorSeqVar(1, NULL, FC_DT_DOUBLE, opStrings[0],
                              glb_numStep, seqVar,
                              newSeqVarName, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail if NULL constValues");
  fail_unless(newSeqVar == NULL,"fail should return NULL");
  rc = fc_constOperatorSeqVar(1, dblConst, FC_DT_UNKNOWN, opStrings[0],
                              glb_numStep, seqVar,
                              newSeqVarName, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail if unknown data type");
  fail_unless(newSeqVar == NULL,"fail should return NULL");
  rc = fc_constOperatorSeqVar(1, dblConst, FC_DT_DOUBLE, opStrings[0],
                              glb_numStep, seqVar, NULL, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail with NULL name");
  fail_unless(newSeqVar == NULL,"fail should return NULL");
  rc = fc_constOperatorSeqVar(1, dblConst, FC_DT_DOUBLE, opStrings[0], 
                              glb_numStep, seqVar, newSeqVarName, NULL);
  fail_unless(rc != FC_SUCCESS, "should fail if NULL new var");
}
END_TEST

// test fc_varUnaryFunction()
START_TEST(unaryVar)
{
  FC_ReturnCode rc;
  int i, j, k;
  FC_Variable var, newVar, badVar = { 999, 999 };
  int numFunc = 3;
  double (*dblFuncs[3])(double) = { sin, sqrt, log };
  int numDataPoint, numComp, temp_numDataPoint, temp_numComp;
  char newVarName[20] = "new var";
  FC_AssociationType assoc, temp_assoc;
  FC_MathType mathType, temp_mathType;
  FC_DataType dataType, temp_dataType;
  void* startData, *resultData;

  // --- use global variables created in the fixture

  // --- test all ops on all datatypes & different # of components
  for (i = 0; i < glb_numDataType + 1; i++) {
    if (i < glb_numDataType)
      var = glb_dataTypeVars[i];
    else
      var = glb_vectorVar;

    // get the start data
    rc = fc_getVariableInfo(var, &numDataPoint, &numComp, &assoc, &mathType, 
			    &dataType);
    fail_unless(rc == FC_SUCCESS, "abort: failed to get var info");
    rc = fc_getVariableDataPtr(var, &startData);
    fail_unless(rc == FC_SUCCESS, "abort: failed to get data pointer");

    for (j = 0; j < numFunc; j++) {

      // do it
      rc = fc_varUnaryFunction(var, dblFuncs[j], newVarName, &newVar);
      if (dataType == FC_DT_CHAR) {
	fail_unless(rc != FC_SUCCESS, "should fail to do char data");
	fail_unless(FC_HANDLE_EQUIV(newVar, FC_NULL_VARIABLE),
		    "fail should return NULL");
	continue;
      }
      else
	fail_unless(rc == FC_SUCCESS, "failed to do unary operation");
      
      // test result meta data
      rc = fc_getVariableInfo(newVar, &temp_numDataPoint, &temp_numComp, 
			      &temp_assoc, &temp_mathType, &temp_dataType);
      fail_unless(rc == FC_SUCCESS, "failed to get new  variable info");
      fail_unless(temp_numDataPoint == numDataPoint, 
		  "mismatch of numDataPoint");
      fail_unless(temp_numComp == numComp, "mismatch of numComp");
      fail_unless(temp_assoc == assoc, "mismatch of assoc");
      fail_unless(temp_mathType == mathType, "mismatch of mathType");
      fail_unless(temp_dataType == FC_DT_DOUBLE, "mismatch of dataType");

      // test result big data
      rc = fc_getVariableDataPtr(newVar, &resultData);
      fail_unless(rc == FC_SUCCESS, "failed to get data ptr from new var");
      for (k = 0; k < numDataPoint*numComp; k++) {
	if (dataType == FC_DT_INT)
	  fail_unless(FC_DBL_EQUIV(((double*)resultData)[k],
				   dblFuncs[j](((int*)startData)[k])),
		      "mismatch of reusult data");
	else if (dataType == FC_DT_FLOAT)
	  fail_unless(FC_DBL_EQUIV(((double*)resultData)[k],
				   dblFuncs[j](((float*)startData)[k])),
		      "mismatch of reusult data");
	else if (dataType == FC_DT_DOUBLE)
	  fail_unless(FC_DBL_EQUIV(((double*)resultData)[k],
				   dblFuncs[j](((double*)startData)[k])),
		      "mismatch of reusult data");
      }
      rc = fc_deleteVariable(newVar);
      fail_unless(rc == FC_SUCCESS, "failed to delete new var");
    }
  }

  // --- things should still work with global vars
  rc = fc_varUnaryFunction(glb_glbVar, dblFuncs[0], newVarName, &newVar);
  fail_unless(rc == FC_SUCCESS, "should work with global var");
  fail_unless(fc_isVariableGlobal(newVar), "should be global");
  rc = fc_getVariableDataPtr(glb_glbVar, &startData);
  fail_unless(rc == FC_SUCCESS, "failed to get data ptr from old var");
  rc = fc_getVariableDataPtr(newVar, &resultData);
  fail_unless(rc == FC_SUCCESS, "failed to get data ptr from new var");
  // just checking first value
  fail_unless(FC_DBL_EQUIV(((double*)resultData)[0],
			   dblFuncs[0](((double*)startData)[0])),
              "mismatch of result data");
  rc = fc_deleteVariable(newVar);
  fail_unless(rc == FC_SUCCESS, "failed to delete new var");

  // --- test bad inputs
  var = glb_dataTypeVars[3];
  rc = fc_varUnaryFunction(badVar, dblFuncs[0], newVarName, &newVar);
  fail_unless(rc != FC_SUCCESS, "should fail with bad var");
  fail_unless(FC_HANDLE_EQUIV(newVar, FC_NULL_VARIABLE),
	      "fail should return NULL");
  rc = fc_varUnaryFunction(var, NULL, newVarName, &newVar);
  fail_unless(rc != FC_SUCCESS, "should fail with bad op string");
  fail_unless(FC_HANDLE_EQUIV(newVar, FC_NULL_VARIABLE),
	      "fail should return NULL");
  rc = fc_varUnaryFunction(var, dblFuncs[0], NULL, &newVar);
  fail_unless(rc != FC_SUCCESS, "should fail with NULL name");
  fail_unless(FC_HANDLE_EQUIV(newVar, FC_NULL_VARIABLE),
	      "fail should return NULL");
  rc = fc_varUnaryFunction(var, dblFuncs[0], newVarName, NULL);
  fail_unless(rc != FC_SUCCESS, "should fail if NULL new var");
}
END_TEST

// test fc_varBinaryFunctionVarVar()
START_TEST(binaryVarVar)
{
  FC_ReturnCode rc;
  int i, j, k, m;
  FC_Variable var1, var2, newVar, badVar = { 999, 999 };
  int numDataPoint1, numComp1, numDataPoint2, numComp2;
  int temp_numDataPoint, temp_numComp;
  char newVarName[20] = "new var";
  FC_AssociationType assoc1, assoc2, temp_assoc;
  FC_MathType mathType1, mathType2, temp_mathType;
  FC_DataType dataType1, dataType2, temp_dataType;
  void* startData1, *startData2, *resultData;
  int numFunc = 3; // arbitrary number, pick stuff that order matters
  double (*dblFuncs[3])(double, double) = { pow, subtraction, magnitude };

  // --- use global variables created in the fixture

  // --- test all operations on all combinations of all var types
  for (i = 0; i < glb_numDataType + 2; i++) {

    if (i < glb_numDataType)
      var1 = glb_dataTypeVars[i];
    else if (i == glb_numDataType)
      var1 = glb_vectorVar;
    else
      var1 = glb_elemVar;

    for (j = 0; j < glb_numDataType + 2; j++) {
      if (j < glb_numDataType)
	var2 = glb_dataTypeVars[j];
      else if (j == glb_numDataType)
	var2 = glb_vectorVar;
      else
	var2 = glb_elemVar;

      // get the start data
      rc = fc_getVariableInfo(var1, &numDataPoint1, &numComp1, &assoc1, &mathType1, 
			      &dataType1);
      fail_unless(rc == FC_SUCCESS, "abort: failed to get var info");
      rc = fc_getVariableDataPtr(var1, &startData1);
      fail_unless(rc == FC_SUCCESS, "abort: failed to get data pointer");
      rc = fc_getVariableInfo(var2, &numDataPoint2, &numComp2, &assoc2, &mathType2, 
			 &dataType2);
      fail_unless(rc == FC_SUCCESS, "abort: failed to get var info");
      rc = fc_getVariableDataPtr(var2, &startData2);
      fail_unless(rc == FC_SUCCESS, "abort: failed to get data pointer");

      // loop over funcs
      for (k = 0; k < numFunc; k++) {
	
	// do it
	rc = fc_varBinaryFunctionVarVar(var1, var2, dblFuncs[k], newVarName, 
					&newVar);
	if (dataType1 == FC_DT_CHAR || dataType2 == FC_DT_CHAR ||
	    numDataPoint1 != numDataPoint2 || numComp1 != numComp2) {
	  fail_unless(rc != FC_SUCCESS, "should fail if char data or if "
		      "vars don't have same numDataPoint or don't have same "
		      "numComp");
	  fail_unless(FC_HANDLE_EQUIV(newVar, FC_NULL_VARIABLE),
		      "fail should return NULL");
	  continue;
	}
	else
	  fail_unless(rc == FC_SUCCESS, "failed to do binary operation");

	// test result meta data
	rc = fc_getVariableInfo(newVar, &temp_numDataPoint, &temp_numComp, 
				&temp_assoc, &temp_mathType, &temp_dataType);
	fail_unless(rc == FC_SUCCESS, "failed to get new variable info");
	fail_unless(temp_numDataPoint == numDataPoint1,
		    "mismatch of numDataPoint");
	fail_unless(temp_numComp == numComp1, "mismatch of numComp");
	fail_unless(temp_assoc == assoc1, "mismatch of assoc");
	fail_unless(temp_mathType == mathType1, "mismatch of mathType");
	fail_unless(temp_dataType == FC_DT_DOUBLE, "mismatch of dataType");

	// test result big data
	rc = fc_getVariableDataPtr(newVar, &resultData);
	fail_unless(rc == FC_SUCCESS, "failed to get data ptr from new var");
	for (m = 0; m < numDataPoint1*numComp1; m++) {
	  double start1, start2;
	  switch (dataType1) {
	  case FC_DT_INT:     start1 =    ((int*)startData1)[m];  break;
	  case FC_DT_FLOAT:   start1 =  ((float*)startData1)[m];  break;
	  case FC_DT_DOUBLE:  start1 = ((double*)startData1)[m];  break;
	  default: ;
	  }
	  switch (dataType2) {
	  case FC_DT_INT:     start2 =    ((int*)startData2)[m];  break;
	  case FC_DT_FLOAT:   start2 =  ((float*)startData2)[m];  break;
	  case FC_DT_DOUBLE:  start2 = ((double*)startData2)[m];  break;
	  default: ;
	  }
	  fail_unless(FC_DBL_EQUIV(((double*)resultData)[m],
				   dblFuncs[k](start1, start2)),
		      "mismatch of result data");
	}
	rc = fc_deleteVariable(newVar);
	fail_unless(rc == FC_SUCCESS, "failed to delete new var");
      }
    }   
  }
 
  // --- things should still work with global vars
  rc = fc_varBinaryFunctionVarVar(glb_glbVar, glb_glbVar, dblFuncs[0], 
                                  newVarName, &newVar);
  fail_unless(rc == FC_SUCCESS, "should work with global var");
  fail_unless(fc_isVariableGlobal(newVar), "should be global");
  rc = fc_getVariableDataPtr(glb_glbVar, &startData1);
  fail_unless(rc == FC_SUCCESS, "failed to get data ptr from old var");
  rc = fc_getVariableDataPtr(newVar, &resultData);
  fail_unless(rc == FC_SUCCESS, "failed to get data ptr from new var");
  // just checking first value
  fail_unless(((double*)resultData)[0] == dblFuncs[0](((double*)startData1)[0],
                                                    ((double*)startData1)[0]),
              "mismatch of reusult data");
  rc = fc_deleteVariable(newVar);
  fail_unless(rc == FC_SUCCESS, "failed to delete new var");

  // --- things should fail with mixed global non-global vars
  var1 = glb_dataTypeVars[3];
  rc = fc_varBinaryFunctionVarVar(glb_glbVar, var1, dblFuncs[0], 
                                  newVarName, &newVar);
  fail_unless(rc != FC_SUCCESS, "should fail with mix of global/nonglobal");
  fail_unless(FC_HANDLE_EQUIV(newVar, FC_NULL_VARIABLE),
	      "fail should return NULL");
  rc = fc_varBinaryFunctionVarVar(var1, glb_glbVar, dblFuncs[0],
                                  newVarName, &newVar);
  fail_unless(rc != FC_SUCCESS, "should fail with mix of global/nonglobal");
  fail_unless(FC_HANDLE_EQUIV(newVar, FC_NULL_VARIABLE),
	      "fail should return NULL");
   
  // --- test bad inputs
  var1 = var2 = glb_dataTypeVars[3];
  rc = fc_varBinaryFunctionVarVar(badVar, var2, dblFuncs[0], newVarName, &newVar);
  fail_unless(rc != FC_SUCCESS, "should fail with bad var 1");
  fail_unless(FC_HANDLE_EQUIV(newVar, FC_NULL_VARIABLE),
	      "fail should return NULL");
  rc = fc_varBinaryFunctionVarVar(var1, badVar, dblFuncs[0], newVarName, &newVar);
  fail_unless(rc != FC_SUCCESS, "should fail with bad var 1");
  fail_unless(FC_HANDLE_EQUIV(newVar, FC_NULL_VARIABLE),
	      "fail should return NULL");
  rc = fc_varBinaryFunctionVarVar(var1, var2, NULL, newVarName, &newVar);
  fail_unless(rc != FC_SUCCESS, "should fail with NULL op");
  fail_unless(FC_HANDLE_EQUIV(newVar, FC_NULL_VARIABLE),
	      "fail should return NULL");
  rc = fc_varBinaryFunctionVarVar(var1, var2, dblFuncs[0], NULL, &newVar);
  fail_unless(rc != FC_SUCCESS, "should fail with null new name");
  fail_unless(FC_HANDLE_EQUIV(newVar, FC_NULL_VARIABLE),
	      "fail should return NULL");
  rc = fc_varBinaryFunctionVarVar(var1, var2, dblFuncs[0], newVarName, NULL);
  fail_unless(rc != FC_SUCCESS, "should fail with null newVar");
}
END_TEST

// test fc_varBinaryOperatorVarConst()  
START_TEST(binaryVarConst)
{
  FC_ReturnCode rc;
  int i, j, k, m, n;
  FC_Variable var, newVar, badVar = { 999, 999 };
  int numDataPoint, numComp, temp_numDataPoint, temp_numComp;
  char newVarName[20] = "new var";
  FC_AssociationType assoc, temp_assoc;
  FC_MathType mathType, temp_mathType;
  FC_DataType dataType, temp_dataType;
  void* startData, *resultData;
  int numFunc = 3;
  double (*dblFuncs[3])(double, double) = { pow, subtraction, magnitude };
  int numConst = 2; // whether scalar or vector
  double dblConst[3] = { 1.2, .00001, 29.9 };
  int constNumComps[2] = { 1, 3 };

  // --- use global variables created in the fixture

  // --- test alls operation on all combinations of vars & consts
  for (i = 0; i < glb_numDataType + 1; i++) {
    if (i < glb_numDataType)
      var = glb_dataTypeVars[i];
    else
      var = glb_vectorVar;
    
    for (j = 0; j < numConst; j++) {
   
      // get the start data
      rc = fc_getVariableInfo(var, &numDataPoint, &numComp, &assoc, &mathType, 
			      &dataType);
      fail_unless(rc == FC_SUCCESS, "abort: failed to get var info");
      rc = fc_getVariableDataPtr(var, &startData);
      fail_unless(rc == FC_SUCCESS, "abort: failed to get data pointer");

      // loop over ops
      for (k = 0; k < numFunc; k++) {

	// do it
	rc = fc_varBinaryFunctionVarConst(var, constNumComps[j], dblConst, 
					  dblFuncs[k], newVarName, &newVar);
	if (dataType == FC_DT_CHAR || numComp != constNumComps[j]) {
	  fail_unless(rc != FC_SUCCESS, "should fail if char data or if "
		      "don't have same numComp");
	  fail_unless(FC_HANDLE_EQUIV(newVar, FC_NULL_VARIABLE),
		      "fail should return NULL");
	  continue;
	}
	else
	  fail_unless(rc == FC_SUCCESS, "failed to do binary operation");

	// test result meta data
	rc = fc_getVariableInfo(newVar, &temp_numDataPoint, &temp_numComp, 
				&temp_assoc, &temp_mathType, &temp_dataType);
	fail_unless(rc == FC_SUCCESS, "failed to get new  variable info");
	fail_unless(temp_numDataPoint == numDataPoint, 
		    "mismatch of numDataPoint");
	fail_unless(temp_numComp == numComp, "mismatch of numComp");
	fail_unless(temp_assoc == assoc, "mismatch of assoc");
	fail_unless(temp_mathType == mathType, "mismatch of mathType");
	fail_unless(temp_dataType == FC_DT_DOUBLE, "mismatch of dataType");

	// test result big data
	rc = fc_getVariableDataPtr(newVar, &resultData);
	fail_unless(rc == FC_SUCCESS, "failed to get data ptr from new var");
	for (m = 0; m < numDataPoint; m++) {
	  for (n = 0; n < numComp; n++) {
	    int idx = m*numComp + n;
	    if (dataType == FC_DT_INT) {
	      fail_unless(FC_DBL_EQUIV(((double*)resultData)[idx],
			  dblFuncs[k](((int*)startData)[idx], dblConst[n])),
			  "mismatch of result data");
	    }
	    else if (dataType == FC_DT_FLOAT) {
	      fail_unless(FC_DBL_EQUIV(((double*)resultData)[idx],
			  dblFuncs[k](((float*)startData)[idx], dblConst[n])),
			  "mismatch of result data");
	    }
	    else if (dataType == FC_DT_DOUBLE) {
	      fail_unless(FC_DBL_EQUIV(((double*)resultData)[idx],
			  dblFuncs[k](((double*)startData)[idx], dblConst[n])),
			  "mismatch of result data");
	    }
	  }
	}
	rc = fc_deleteVariable(newVar);
	fail_unless(rc == FC_SUCCESS, "failed to delete new var");
      }
    }
  }
 
  // --- things should still work with global vars
  rc = fc_varBinaryFunctionVarConst(glb_glbVar, 1, dblConst, dblFuncs[0],
                           newVarName, &newVar);
  fail_unless(rc == FC_SUCCESS, "should work with global var");
  fail_unless(fc_isVariableGlobal(newVar), "should be global");
  rc = fc_getVariableDataPtr(glb_glbVar, &startData);
  fail_unless(rc == FC_SUCCESS, "failed to get data ptr from old var");
  rc = fc_getVariableDataPtr(newVar, &resultData);
  fail_unless(rc == FC_SUCCESS, "failed to get data ptr from new var");
  // just checking first value
  fail_unless(((double*)resultData)[0] == dblFuncs[0](((double*)startData)[0],
                                                    dblConst[0]),
              "mismatch of reusult data");
  rc = fc_deleteVariable(newVar);
  fail_unless(rc == FC_SUCCESS, "failed to delete new var");

  // --- test bad inputs
  var = glb_dataTypeVars[3];
  rc = fc_varBinaryFunctionVarConst(badVar, 1, dblConst, dblFuncs[0], 
				    newVarName, &newVar);
  fail_unless(rc != FC_SUCCESS, "should fail with bad var");
  fail_unless(FC_HANDLE_EQUIV(newVar, FC_NULL_VARIABLE),
	      "fail should return NULL");
  rc = fc_varBinaryFunctionVarConst(var, 1, dblConst, NULL,
				    newVarName, &newVar);
  fail_unless(rc != FC_SUCCESS, "should fail with bad op string");
  fail_unless(FC_HANDLE_EQUIV(newVar, FC_NULL_VARIABLE),
	      "fail should return NULL");
  rc = fc_varBinaryFunctionVarConst(var, 0, dblConst, dblFuncs[0], 
				    newVarName, &newVar);
  fail_unless(rc != FC_SUCCESS, "should fail if constLen < 1");
  fail_unless(FC_HANDLE_EQUIV(newVar, FC_NULL_VARIABLE),
	      "fail should return NULL");
  rc = fc_varBinaryFunctionVarConst(var, 1, NULL, dblFuncs[0], 
				    newVarName, &newVar);
  fail_unless(rc != FC_SUCCESS, "should fail if constValues = NULL");
  fail_unless(FC_HANDLE_EQUIV(newVar, FC_NULL_VARIABLE),
	      "fail should return NULL");
  rc = fc_varBinaryFunctionVarConst(var, 1, dblConst, dblFuncs[0],
			   NULL, &newVar);
  fail_unless(rc != FC_SUCCESS, "should fail with null new name");
  fail_unless(FC_HANDLE_EQUIV(newVar, FC_NULL_VARIABLE),
	      "fail should return NULL");
  rc = fc_varBinaryFunctionVarConst(var, 1, dblConst, dblFuncs[0],
			   newVarName, NULL);
  fail_unless(rc != FC_SUCCESS, "should fail with null newVar");
}
END_TEST

// test fc_varBinaryOperatorConstVar()  
START_TEST(binaryConstVar)
{
  FC_ReturnCode rc;
  int i, j, k, m, n;
  FC_Variable var, newVar, badVar = { 999, 999 };
  int numDataPoint, numComp, temp_numDataPoint, temp_numComp;
  char newVarName[20] = "new var";
  FC_AssociationType assoc, temp_assoc;
  FC_MathType mathType, temp_mathType;
  FC_DataType dataType, temp_dataType;
  void* startData, *resultData;
  int numFunc = 3;
  double (*dblFuncs[3])(double, double) = { pow, subtraction, magnitude };
  int numConst = 2; // whether scalar or vector
  double dblConst[3] = { 1.2, .00001, 29.9 };
  int constNumComps[2] = { 1, 3 };

  // --- use global variables created in the fixture

  // --- test alls operation on all combinations of vars & consts
  for (i = 0; i < glb_numDataType + 1; i++) {
    if (i < glb_numDataType)
      var = glb_dataTypeVars[i];
    else
      var = glb_vectorVar;
    
    for (j = 0; j < numConst; j++) {
   
      // get the start data
      rc = fc_getVariableInfo(var, &numDataPoint, &numComp, &assoc, &mathType, 
			      &dataType);
      fail_unless(rc == FC_SUCCESS, "abort: failed to get var info");
      rc = fc_getVariableDataPtr(var, &startData);
      fail_unless(rc == FC_SUCCESS, "abort: failed to get data pointer");

      // loop over ops
      for (k = 0; k < numFunc; k++) {

	// do it
	rc = fc_varBinaryFunctionConstVar(constNumComps[j], dblConst, var,
					  dblFuncs[k], newVarName, &newVar);
	if (dataType == FC_DT_CHAR || numComp != constNumComps[j]) {
	  fail_unless(rc != FC_SUCCESS, "should fail if char data or if "
		      "don't have same numComp");
	  fail_unless(FC_HANDLE_EQUIV(newVar, FC_NULL_VARIABLE),
		      "fail should return NULL");
	  continue;
	}
	else
	  fail_unless(rc == FC_SUCCESS, "failed to do binary operation");

	// test result meta data
	rc = fc_getVariableInfo(newVar, &temp_numDataPoint, &temp_numComp, 
				&temp_assoc, &temp_mathType, &temp_dataType);
	fail_unless(rc == FC_SUCCESS, "failed to get new  variable info");
	fail_unless(temp_numDataPoint == numDataPoint, 
		    "mismatch of numDataPoint");
	fail_unless(temp_numComp == numComp, "mismatch of numComp");
	fail_unless(temp_assoc == assoc, "mismatch of assoc");
	fail_unless(temp_mathType == mathType, "mismatch of mathType");
	fail_unless(temp_dataType == FC_DT_DOUBLE, "mismatch of dataType");

	// test result big data
	rc = fc_getVariableDataPtr(newVar, &resultData);
	fail_unless(rc == FC_SUCCESS, "failed to get data ptr from new var");
	for (m = 0; m < numDataPoint; m++) {
	  for (n = 0; n < numComp; n++) {
	    int idx = m*numComp + n;
	    if (dataType == FC_DT_INT) {
	      fail_unless(FC_DBL_EQUIV(((double*)resultData)[idx],
			  dblFuncs[k](dblConst[n], ((int*)startData)[idx])),
			  "mismatch of result data");
	    }
	    else if (dataType == FC_DT_FLOAT) {
	      fail_unless(FC_DBL_EQUIV(((double*)resultData)[idx],
			  dblFuncs[k](dblConst[n], ((float*)startData)[idx])),
			  "mismatch of result data");
	    }
	    else if (dataType == FC_DT_DOUBLE) {
	      fail_unless(FC_DBL_EQUIV(((double*)resultData)[idx],
			  dblFuncs[k](dblConst[n], ((double*)startData)[idx])),
			  "mismatch of result data");
	    }
	  }
	}
	rc = fc_deleteVariable(newVar);
	fail_unless(rc == FC_SUCCESS, "failed to delete new var");
      }
    }
  }
 
  // --- things should still work with global vars
  rc = fc_varBinaryFunctionConstVar(1, dblConst, glb_glbVar, dblFuncs[0],
                           newVarName, &newVar);
  fail_unless(rc == FC_SUCCESS, "should work with global var");
  fail_unless(fc_isVariableGlobal(newVar), "should be global");
  rc = fc_getVariableDataPtr(glb_glbVar, &startData);
  fail_unless(rc == FC_SUCCESS, "failed to get data ptr from old var");
  rc = fc_getVariableDataPtr(newVar, &resultData);
  fail_unless(rc == FC_SUCCESS, "failed to get data ptr from new var");
  // just checking first value
  fail_unless(((double*)resultData)[0] == dblFuncs[0](dblConst[0],
                                                    ((double*)startData)[0]),
              "mismatch of reusult data");
  rc = fc_deleteVariable(newVar);
  fail_unless(rc == FC_SUCCESS, "failed to delete new var");

  // --- test bad inputs
  var = glb_dataTypeVars[3];
  rc = fc_varBinaryFunctionConstVar(1, dblConst, badVar, dblFuncs[0], 
				    newVarName, &newVar);
  fail_unless(rc != FC_SUCCESS, "should fail with bad var");
  fail_unless(FC_HANDLE_EQUIV(newVar, FC_NULL_VARIABLE),
	      "fail should return NULL");
  rc = fc_varBinaryFunctionConstVar(1, dblConst, var, NULL,
				    newVarName, &newVar);
  fail_unless(rc != FC_SUCCESS, "should fail with bad op string");
  fail_unless(FC_HANDLE_EQUIV(newVar, FC_NULL_VARIABLE),
	      "fail should return NULL");
  rc = fc_varBinaryFunctionConstVar(0, dblConst, var, dblFuncs[0], 
				    newVarName, &newVar);
  fail_unless(rc != FC_SUCCESS, "should fail if constLen < 1");
  fail_unless(FC_HANDLE_EQUIV(newVar, FC_NULL_VARIABLE),
	      "fail should return NULL");
  rc = fc_varBinaryFunctionConstVar(1, NULL, var, dblFuncs[0],
				    newVarName, &newVar);
  fail_unless(rc != FC_SUCCESS, "should fail if constValues = NULL");
  fail_unless(FC_HANDLE_EQUIV(newVar, FC_NULL_VARIABLE),
	      "fail should return NULL");
  rc = fc_varBinaryFunctionConstVar(1, dblConst, var, dblFuncs[0],
			   NULL, &newVar);
  fail_unless(rc != FC_SUCCESS, "should fail with null new name");
  fail_unless(FC_HANDLE_EQUIV(newVar, FC_NULL_VARIABLE),
	      "fail should return NULL");
  rc = fc_varBinaryFunctionConstVar(1, dblConst, var, dblFuncs[0],
			   newVarName, NULL);
  fail_unless(rc != FC_SUCCESS, "should fail with null newVar");
}
END_TEST

// test fc_varUnaryFunction()
START_TEST(unarySeqVar)
{
  FC_ReturnCode rc;
  int i, j, k, m;
  FC_Sequence temp_seq;
  FC_Variable *seqVar, *newSeqVar, badSeqVar[glb_numStep], badVar = { 999, 999 };
  int numFunc = 3; // just pick four C math functions
  double (*dblFuncs[3])(double) = { sin, sqrt, log };
  int numDataPoint, numComp, temp_numDataPoint, temp_numComp;
  char newSeqVarName[20] = "new seq var";
  FC_AssociationType assoc, temp_assoc;
  FC_MathType mathType, temp_mathType;
  FC_DataType dataType, temp_dataType;
  void* startData, *resultData;

  // --- use global variables created in the fixture

  // --- test all ops on all datatypes & different # of components
  for (i = 0; i < glb_numDataType + 1; i++) {
   if (i < glb_numDataType)
      seqVar = glb_dataTypeSeqVars[i];
    else
      seqVar = glb_vectorSeqVar;

    // get the start metadata
    rc = fc_getVariableInfo(seqVar[0], &numDataPoint, &numComp, &assoc,
		       &mathType, &dataType);
    fail_unless(rc == FC_SUCCESS, "cant get var info");

    for (j = 0; j < numFunc; j++) {

      // do it
      rc = fc_seqVarUnaryFunction(glb_numStep, seqVar, dblFuncs[j],
				  newSeqVarName, &newSeqVar);
      if (dataType == FC_DT_CHAR) {
	fail_unless(rc != FC_SUCCESS, "should fail to do char data");
	fail_unless(newSeqVar == NULL, "fail should return NULL");
	continue;
      } 
      else {
	fail_unless(rc == FC_SUCCESS, "failed to do unary operation");
      }
      
      // test result meta data
      rc = fc_getVariableInfo(newSeqVar[0], &temp_numDataPoint, &temp_numComp, 
			      &temp_assoc, &temp_mathType, &temp_dataType);
      fail_unless(rc == FC_SUCCESS, "failed to get new variable info");
      fail_unless(temp_numDataPoint == numDataPoint, 
		  "mismatch of numDataPoint");
      fail_unless(temp_numComp == numComp, "mismatch of numComp");
      fail_unless(temp_assoc == assoc, "mismatch of assoc");
      fail_unless(temp_mathType == mathType, "mismatch of mathType");
      fail_unless(temp_dataType == FC_DT_DOUBLE, "mismatch of datatype");
      rc = fc_getSequenceFromSeqVariable(glb_numStep, newSeqVar, &temp_seq);
      fail_unless(rc == FC_SUCCESS, "abort: failed to get new sequence");
      fail_unless(FC_HANDLE_EQUIV(temp_seq, glb_sequence), "mismatch of sequence");

      // test result big data
      for (k = 0; k < glb_numStep; k++) {
	rc = fc_getVariableDataPtr(seqVar[k], &startData);
	fail_unless(rc == FC_SUCCESS, "abort: failed to get data pointer");
	rc = fc_getVariableDataPtr(newSeqVar[k], &resultData);
	fail_unless(rc == FC_SUCCESS, "failed to get data ptr from new var");
	for (m = 0; m < numDataPoint*numComp; m++) {
	  if (dataType == FC_DT_INT)
	    fail_unless(FC_DBL_EQUIV(((double*)resultData)[m],
				     dblFuncs[j](((int*)startData)[m])),
			"mismatch of reusult data");
	  else if (dataType == FC_DT_FLOAT)
	    fail_unless(FC_DBL_EQUIV(((double*)resultData)[m],
				     dblFuncs[j](((float*)startData)[m])),
			"mismatch of reusult data");
	  else if (dataType == FC_DT_DOUBLE)
	    fail_unless(FC_DBL_EQUIV(((double*)resultData)[m],
				     dblFuncs[j](((double*)startData)[m])),
			"mismatch of reusult data");
	}
      }
      rc = fc_deleteSeqVariable(glb_numStep, newSeqVar);
      fail_unless(rc == FC_SUCCESS, "failed to delete new seqvar");
      free(newSeqVar);
    }
  }

  // --- things should still work with global vars
  rc = fc_seqVarUnaryFunction(glb_numStep, glb_glbSeqVar, dblFuncs[0], 
                         newSeqVarName, &newSeqVar);
  fail_unless(rc == FC_SUCCESS, "should work with global var");
  fail_unless(fc_isVariableGlobal(newSeqVar[0]), "should be global");
  for (i = 0; i < glb_numStep; i++) {
    rc = fc_getVariableDataPtr(glb_glbSeqVar[i], &startData);
    fail_unless(rc == FC_SUCCESS, "failed to get data ptr from old var");
    rc = fc_getVariableDataPtr(newSeqVar[i], &resultData);
    fail_unless(rc == FC_SUCCESS, "failed to get data ptr from new var");
    // just checking first value
    fail_unless(FC_DBL_EQUIV(((double*)resultData)[0],
			     dblFuncs[0](((double*)startData)[0])),
                "mismatch of result data");
  }
  rc = fc_deleteSeqVariable(glb_numStep, newSeqVar);
  fail_unless(rc == FC_SUCCESS, "failed to delete new var");
  free(newSeqVar);

  // --- test bad inputs
  seqVar = glb_dataTypeSeqVars[3];
  for (i = 0; i < glb_numStep; i++)
    badSeqVar[i] = glb_dataTypeSeqVars[3][i];
  badSeqVar[1] = badVar;
  rc = fc_seqVarUnaryFunction(glb_numStep, NULL, dblFuncs[0],
			      newSeqVarName, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail with bad input var");
  fail_unless(newSeqVar == NULL,"fail should return NULL");
  rc = fc_seqVarUnaryFunction(glb_numStep, badSeqVar, dblFuncs[0],
			      newSeqVarName, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail with bad input var");
  fail_unless(newSeqVar == NULL,"fail should return NULL");
  rc = fc_seqVarUnaryFunction(glb_numStep, seqVar, NULL,
			      newSeqVarName, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail with NULL op");
  fail_unless(newSeqVar == NULL,"fail should return NULL");
  rc = fc_seqVarUnaryFunction(glb_numStep, seqVar, dblFuncs[0],
			      NULL, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail with NULL name");
  fail_unless(newSeqVar == NULL,"fail should return NULL");
  rc = fc_seqVarUnaryFunction(glb_numStep, seqVar, dblFuncs[0],
			      newSeqVarName, NULL);
  fail_unless(rc != FC_SUCCESS, "should fail with NULL new var");
}
END_TEST

// test fc_seqVarBinaryFunctionSeqVarSeqVar()
START_TEST(binarySeqVarSeqVar)
{
  FC_ReturnCode rc;
  int i, j, k, m, n;
  FC_Sequence temp_seq, seq2;
  FC_Variable *seqVar1, *seqVar2, *newSeqVar, badSeqVar[3], badVar = { 999, 999 };
  FC_Variable *otherSeqVar;
  int numDataPoint1, numComp1, numDataPoint2, numComp2;
  int temp_numDataPoint, temp_numComp;
  char newSeqVarName[20] = "new seq var";
  FC_AssociationType assoc1, assoc2, temp_assoc;
  FC_MathType mathType1, mathType2, temp_mathType;
  FC_DataType dataType1, dataType2, temp_dataType;
  void* startData1, *startData2, *resultData;
  int numFunc = 3;
  double (*dblFuncs[3])(double, double) = { pow, subtraction, magnitude };

  // --- use global variables created in the fixture

  // --- test all ops on all combinations of all var types
  for (i = 0; i < glb_numDataType + 2; i++) {

    if (i < glb_numDataType)
      seqVar1 = glb_dataTypeSeqVars[i];
    else if (i == glb_numDataType)
      seqVar1 = glb_vectorSeqVar;
    else 
      seqVar1 = glb_vectorSeqVar;

    for (j = 0; j < glb_numDataType + 2; j++) {
      if (j < glb_numDataType)
	seqVar2 = glb_dataTypeSeqVars[j];
      else if (j == glb_numDataType)
	seqVar2 = glb_vectorSeqVar;
      else
	seqVar2 = glb_elemSeqVar;

      // get the start meta data
      rc = fc_getVariableInfo(seqVar1[0], &numDataPoint1, &numComp1, &assoc1, 
			      &mathType1,  &dataType1);
      fail_unless(rc == FC_SUCCESS, "abort: failed to get var info");
      rc = fc_getVariableInfo(seqVar2[0], &numDataPoint2, &numComp2, &assoc2, 
			      &mathType2, &dataType2);
      fail_unless(rc == FC_SUCCESS, "abort: failed to get var info");

      // loop over ops
      for (k = 0; k < numFunc; k++) {

	// do it
	rc = fc_seqVarBinaryFunctionSeqVarSeqVar(glb_numStep, seqVar1, seqVar2,
					dblFuncs[k], newSeqVarName, &newSeqVar);
	if (dataType1 == FC_DT_CHAR || dataType2 == FC_DT_CHAR ||
	    numDataPoint1 != numDataPoint2 || numComp1 != numComp2) {
	  fail_unless(rc != FC_SUCCESS, "should fail to do char data or if "
		      "vars don't have some numDataPoint or don't have same "
		      "numComp");
	  fail_unless(newSeqVar == NULL, "fail should return NULL");
	  continue;
	}
	else
	  fail_unless(rc == FC_SUCCESS, "failed to do binary operation");

	// test result meta data
	rc = fc_getVariableInfo(newSeqVar[0], &temp_numDataPoint, &temp_numComp, 
			      &temp_assoc, &temp_mathType, &temp_dataType);
	fail_unless(rc == FC_SUCCESS, "failed to get new variable info");
	fail_unless(temp_numDataPoint == numDataPoint1, 
		    "mismatch of numDataPoint");
	fail_unless(temp_numComp == numComp1, "mismatch of numComp");
	fail_unless(temp_assoc == assoc1, "mismatch of assoc");
	fail_unless(temp_mathType == mathType1, "mismatch of mathType");
	fail_unless(temp_dataType == FC_DT_DOUBLE, "mismatch of dataType");
	rc = fc_getSequenceFromSeqVariable(glb_numStep, newSeqVar, &temp_seq);
	fail_unless(rc == FC_SUCCESS, "abort: failed to get new sequence");
	fail_unless(FC_HANDLE_EQUIV(temp_seq, glb_sequence), "mismatch of sequence");

	// test result big data
	for (m = 0; m < glb_numStep; m++) {	
	  rc = fc_getVariableDataPtr(seqVar1[m], &startData1);
	  fail_unless(rc == FC_SUCCESS, "abort: failed to get data pointer");
	  rc = fc_getVariableDataPtr(seqVar2[m], &startData2);
	  fail_unless(rc == FC_SUCCESS, "abort: failed to get data pointer");
	  rc = fc_getVariableDataPtr(newSeqVar[m], &resultData);
	  fail_unless(rc == FC_SUCCESS, "failed to get data ptr from new var");
	  for (n = 0; n < numDataPoint1*numComp1; n++) {
	    double start1, start2;
	    switch (dataType1) {
	    case FC_DT_INT:     start1 =    ((int*)startData1)[n];  break;
	    case FC_DT_FLOAT:   start1 =  ((float*)startData1)[n];  break;
	    case FC_DT_DOUBLE:  start1 = ((double*)startData1)[n];  break;
	    default: ;
	    }
	    switch (dataType2) {
	    case FC_DT_INT:     start2 =    ((int*)startData2)[n];  break;
	    case FC_DT_FLOAT:   start2 =  ((float*)startData2)[n];  break;
	    case FC_DT_DOUBLE:  start2 = ((double*)startData2)[n];  break;
	    default: ;
	    }
	    fail_unless(FC_DBL_EQUIV(((double*)resultData)[n],
				     dblFuncs[k](start1, start2)),
			"mismatch of result data");
	  }
	}
	rc = fc_deleteSeqVariable(glb_numStep, newSeqVar);
	fail_unless(rc == FC_SUCCESS, "failed to delete new seqvar");
	free(newSeqVar);
      }
    }
  }

  // --- things should still work with global vars
  rc = fc_seqVarBinaryFunctionSeqVarSeqVar(glb_numStep, glb_glbSeqVar, glb_glbSeqVar, 
                                           dblFuncs[0], newSeqVarName, &newSeqVar);
  fail_unless(rc == FC_SUCCESS, "should work with global var");
  fail_unless(fc_isVariableGlobal(newSeqVar[0]), "should be global");
  for (i = 0; i < glb_numStep; i++) {
    rc = fc_getVariableDataPtr(glb_glbSeqVar[i], &startData1);
    fail_unless(rc == FC_SUCCESS, "failed to get data ptr from old var");
    rc = fc_getVariableDataPtr(newSeqVar[i], &resultData);
    fail_unless(rc == FC_SUCCESS, "failed to get data ptr from new var");
    // just checking first value
    fail_unless(((double*)resultData)[0] == dblFuncs[0](((double*)startData1)[0],
                                                      ((double*)startData1)[0]),
                "mismatch of result data");
  }
  rc = fc_deleteSeqVariable(glb_numStep, newSeqVar);
  fail_unless(rc == FC_SUCCESS, "failed to delete new var");
  free(newSeqVar);

  // --- things should fail with mixed global non-global vars
  seqVar1 = glb_dataTypeSeqVars[3];
  rc = fc_seqVarBinaryFunctionSeqVarSeqVar(glb_numStep, glb_glbSeqVar, seqVar1,
                                           dblFuncs[0], newSeqVarName, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail with mix of global/nonglobal");
  fail_unless(newSeqVar == NULL, "fail should return NULL");
  rc = fc_seqVarBinaryFunctionSeqVarSeqVar(glb_numStep, seqVar1, glb_glbSeqVar, 
                                           dblFuncs[0], newSeqVarName, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail with mix of global/nonglobal");
  fail_unless(newSeqVar == NULL, "fail should return NULL");

  // --- It is an error to have seq variables on different sequences
  seqVar1 = glb_dataTypeSeqVars[3];
  rc = fc_copySequence(glb_sequence, glb_dataset, "seq2", &seq2);
  fail_unless(rc == FC_SUCCESS, "failed to copy sequence");
  rc = fc_copySeqVariable(glb_numStep, seqVar1, glb_mesh, seq2, 
			  "var on another seq", &otherSeqVar);
  fail_unless(rc == FC_SUCCESS, "failed to create seq var on another seq");
  rc = fc_seqVarBinaryFunctionSeqVarSeqVar(glb_numStep, seqVar1, otherSeqVar,
				   dblFuncs[0], newSeqVarName, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail with bad input var");
  fail_unless(newSeqVar == NULL,"fail should return NULL");
  free(otherSeqVar);

  // --- test bad inputs
  seqVar1 = seqVar2 = glb_dataTypeSeqVars[3];
  rc = fc_seqVarBinaryFunctionSeqVarSeqVar(glb_numStep, NULL, seqVar2, 
                                      dblFuncs[0], newSeqVarName, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail with bad input var");
  fail_unless(newSeqVar == NULL,"fail should return NULL");
  for (i = 0; i < glb_numStep; i++)
    badSeqVar[i] = glb_dataTypeSeqVars[3][i];
  badSeqVar[1] = badVar;
  rc = fc_seqVarBinaryFunctionSeqVarSeqVar(glb_numStep, badSeqVar, seqVar2,
				   dblFuncs[0], newSeqVarName, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail with bad input var");
  fail_unless(newSeqVar == NULL,"fail should return NULL");
  rc = fc_seqVarBinaryFunctionSeqVarSeqVar(glb_numStep, seqVar1, seqVar2, 
					   NULL, newSeqVarName, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail with bad op string");
  fail_unless(newSeqVar == NULL,"fail should return NULL");
  rc = fc_seqVarBinaryFunctionSeqVarSeqVar(glb_numStep, seqVar1, NULL,
				   dblFuncs[0], newSeqVarName, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail with bad input var");
  fail_unless(newSeqVar == NULL,"fail should return NULL");
  rc = fc_seqVarBinaryFunctionSeqVarSeqVar(glb_numStep, seqVar1, badSeqVar,
				       dblFuncs[0], newSeqVarName, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail with bad input var");
  fail_unless(newSeqVar == NULL,"fail should return NULL");
  rc = fc_seqVarBinaryFunctionSeqVarSeqVar(glb_numStep, seqVar1, seqVar2,
					   dblFuncs[0], NULL, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail with NULL name");
  fail_unless(newSeqVar == NULL,"fail should return NULL");
  rc = fc_seqVarBinaryFunctionSeqVarSeqVar(glb_numStep, seqVar1, seqVar2,
					   dblFuncs[0], newSeqVarName, NULL);
  fail_unless(rc != FC_SUCCESS, "should fail if NULL new var");
}
END_TEST


//CDU
// test fc_seqVarBinaryFunctionSeqVarConst()
START_TEST(binarySeqVarConst)
{
  FC_ReturnCode rc;
  int i, j, k, m, n;
  FC_Sequence temp_seq;
  FC_Variable *seqVar,  *newSeqVar, badSeqVar[3], badVar = { 999, 999 };
  int numDataPoint, numComp;
  int temp_numDataPoint, temp_numComp;
  char newSeqVarName[20] = "new seq var";
  FC_AssociationType assoc, temp_assoc;
  FC_MathType mathType, temp_mathType;
  FC_DataType dataType, temp_dataType;
  void* startData, *resultData;
  int numFunc = 3;
  double (*dblFuncs[3])(double, double) = { pow, subtraction, magnitude };
  int numConst = 2;
  double dblConst[3] = { 1.2, .00001, 29.9 };
  int constNumComps[2] = { 1, 3 };


  // --- use global variables created in the fixture

  // --- test all ops on all combinations of all var types
  for (i = 0; i < glb_numDataType + 2; i++) {

    if      (i <  glb_numDataType)  seqVar = glb_dataTypeSeqVars[i];
    else if (i == glb_numDataType)  seqVar = glb_vectorSeqVar;
    else                            seqVar = glb_vectorSeqVar;

    for (j = 0; j < numConst; j++) {

      // get the start meta data
      rc = fc_getVariableInfo(seqVar[0], &numDataPoint, &numComp, &assoc, 
			      &mathType,  &dataType);

      fail_unless(rc == FC_SUCCESS, "abort: failed to get var info");

      // loop over ops
      for (k = 0; k < numFunc; k++) {

	// do it
	rc = fc_seqVarBinaryFunctionSeqVarConst(glb_numStep, seqVar, 
						constNumComps[j], dblConst,
						dblFuncs[k], newSeqVarName, &newSeqVar);

	if (dataType == FC_DT_CHAR || numComp != constNumComps[j]) {
	  fail_unless(rc != FC_SUCCESS, "should fail to do char data or if "
		      "vars don't have same numComponent");
	  fail_unless(newSeqVar == NULL, "fail should return NULL");
	  continue;
	}
	else
	  fail_unless(rc == FC_SUCCESS, "failed to do binary operation");

	// test result meta data
	rc = fc_getVariableInfo(newSeqVar[0], &temp_numDataPoint, &temp_numComp, 
			      &temp_assoc, &temp_mathType, &temp_dataType);
	fail_unless(rc == FC_SUCCESS, "failed to get new variable info");
	fail_unless(temp_numDataPoint == numDataPoint, "mismatch of numDataPoint");
	fail_unless(temp_numComp      == numComp,      "mismatch of numComp");
	fail_unless(temp_assoc        == assoc,        "mismatch of assoc");
	fail_unless(temp_mathType     == mathType,     "mismatch of mathType");
	fail_unless(temp_dataType     == FC_DT_DOUBLE, "mismatch of dataType");
	rc = fc_getSequenceFromSeqVariable(glb_numStep, newSeqVar, &temp_seq);
	fail_unless(rc == FC_SUCCESS, "abort: failed to get new sequence");
	fail_unless(FC_HANDLE_EQUIV(temp_seq, glb_sequence), "mismatch of sequence");

	// test result big data
	for (m = 0; m < glb_numStep; m++) {	
	  rc = fc_getVariableDataPtr(seqVar[m], &startData);
	  fail_unless(rc == FC_SUCCESS, "abort: failed to get data pointer");
	  rc = fc_getVariableDataPtr(newSeqVar[m], &resultData);
	  fail_unless(rc == FC_SUCCESS, "failed to get data ptr from new var");
	  for (n = 0; n < numDataPoint*numComp; n++) {
	    double start1, start2;
	    switch (dataType) {
	    case FC_DT_INT:     start1 =    ((int*)startData)[n];  break;
	    case FC_DT_FLOAT:   start1 =  ((float*)startData)[n];  break;
	    case FC_DT_DOUBLE:  start1 = ((double*)startData)[n];  break;
	    default: ;
	    }
	    start2 = dblConst[ n%numComp ];
	    fail_unless(FC_DBL_EQUIV(((double*)resultData)[n],
				     dblFuncs[k](start1, start2)),
			"mismatch of result data");
	  }
	}
	rc = fc_deleteSeqVariable(glb_numStep, newSeqVar);
	fail_unless(rc == FC_SUCCESS, "failed to delete new seqvar");
	free(newSeqVar);
      }
    }
  }


  // --- things should still work with global vars
  rc = fc_seqVarBinaryFunctionSeqVarConst(glb_numStep, glb_glbSeqVar, //gsv component is 1
					  constNumComps[0], dblConst, 
                                          dblFuncs[0], newSeqVarName, &newSeqVar);
  fail_unless(rc == FC_SUCCESS, "should work with global var");
  fail_unless(fc_isVariableGlobal(newSeqVar[0]), "should be global");
  for (i = 0; i < glb_numStep; i++) {
    rc = fc_getVariableDataPtr(glb_glbSeqVar[i], &startData);
    fail_unless(rc == FC_SUCCESS, "failed to get data ptr from old var");
    
    rc = fc_getVariableDataPtr(newSeqVar[i], &resultData);
    fail_unless(rc == FC_SUCCESS, "failed to get data ptr from new var");
    // just checking first value
    fail_unless(((double*)resultData)[0] == dblFuncs[0](((double*)startData)[0],
                                                        ((double*)dblConst)[0]),
                "mismatch of result data");
  }
  rc = fc_deleteSeqVariable(glb_numStep, newSeqVar);
  fail_unless(rc == FC_SUCCESS, "failed to delete new var");
  free(newSeqVar);


  // --- test bad inputs
  seqVar = glb_dataTypeSeqVars[3];
  for (i = 0; i < glb_numStep; i++)
    badSeqVar[i] = glb_dataTypeSeqVars[3][i];
  badSeqVar[1] = badVar;
  rc = fc_seqVarBinaryFunctionSeqVarConst(glb_numStep, NULL, constNumComps[0], dblConst, 
					  dblFuncs[0], newSeqVarName, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail with bad input var");
  fail_unless(newSeqVar == NULL,"fail should return NULL");
  rc = fc_seqVarBinaryFunctionSeqVarConst(glb_numStep, badSeqVar, constNumComps[0], dblConst, 
				   dblFuncs[0], newSeqVarName, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail with bad input var");
  fail_unless(newSeqVar == NULL,"fail should return NULL");
  rc = fc_seqVarBinaryFunctionSeqVarConst(glb_numStep, seqVar, constNumComps[0], dblConst, 
					   NULL, newSeqVarName, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail with bad op string");
  fail_unless(newSeqVar == NULL,"fail should return NULL");
  rc = fc_seqVarBinaryFunctionSeqVarConst(glb_numStep, seqVar, constNumComps[0], NULL, 
				       dblFuncs[0], newSeqVarName, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail with bad input var");
  fail_unless(newSeqVar == NULL,"fail should return NULL");
  rc = fc_seqVarBinaryFunctionSeqVarConst(glb_numStep, seqVar, 0, dblConst, 
				       dblFuncs[0], newSeqVarName, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail with bad input var");
  fail_unless(newSeqVar == NULL,"fail should return NULL");

  rc = fc_seqVarBinaryFunctionSeqVarConst(glb_numStep, seqVar, constNumComps[0], dblConst, 
					   dblFuncs[0], NULL, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail with NULL name");
  fail_unless(newSeqVar == NULL,"fail should return NULL");
  rc = fc_seqVarBinaryFunctionSeqVarConst(glb_numStep, seqVar, constNumComps[0], dblConst, 
					   dblFuncs[0], newSeqVarName, NULL);
  fail_unless(rc != FC_SUCCESS, "should fail if NULL new var");


  //Sanity check to make sure the above calls COULD have worked with right params
  rc = fc_seqVarBinaryFunctionSeqVarConst(glb_numStep, seqVar, constNumComps[0], dblConst, 
					   dblFuncs[0], newSeqVarName, &newSeqVar);
  fail_unless(rc == FC_SUCCESS, "copy should have worked?");
  rc = fc_deleteSeqVariable(glb_numStep, newSeqVar);
  fail_unless(rc == FC_SUCCESS, "failed to delete new var");
  free(newSeqVar);

}
END_TEST


//CDU
// test fc_seqVarBinaryFunctionConstSeqVar()
START_TEST(binaryConstSeqVar)
{
  FC_ReturnCode rc;
  int i, j, k, m, n;
  FC_Sequence temp_seq;
  FC_Variable *seqVar,  *newSeqVar, badSeqVar[3], badVar = { 999, 999 };
  int numDataPoint, numComp;
  int temp_numDataPoint, temp_numComp;
  char newSeqVarName[20] = "new seq var";
  FC_AssociationType assoc, temp_assoc;
  FC_MathType mathType, temp_mathType;
  FC_DataType dataType, temp_dataType;
  void* startData, *resultData;
  int numFunc = 3;
  double (*dblFuncs[3])(double, double) = { pow, subtraction, magnitude };
  int numConst = 2;
  double dblConst[3] = { 1.2, .00001, 29.9 };
  int constNumComps[2] = { 1, 3 };


  // --- use global variables created in the fixture

  // --- test all ops on all combinations of all var types
  for (i = 0; i < glb_numDataType + 2; i++) {

    if      (i <  glb_numDataType)  seqVar = glb_dataTypeSeqVars[i];
    else if (i == glb_numDataType)  seqVar = glb_vectorSeqVar;
    else                            seqVar = glb_vectorSeqVar;

    for (j = 0; j < numConst; j++) {

      // get the start meta data
      rc = fc_getVariableInfo(seqVar[0], &numDataPoint, &numComp, &assoc, 
			      &mathType,  &dataType);

      fail_unless(rc == FC_SUCCESS, "abort: failed to get var info");

      // loop over ops
      for (k = 0; k < numFunc; k++) {

	// do it
	rc = fc_seqVarBinaryFunctionConstSeqVar(constNumComps[j], dblConst,
						glb_numStep, seqVar, 
						dblFuncs[k], newSeqVarName, &newSeqVar);

	if (dataType == FC_DT_CHAR || numComp != constNumComps[j]) {
	  fail_unless(rc != FC_SUCCESS, "should fail to do char data or if "
		      "vars don't have same numComponent");
	  fail_unless(newSeqVar == NULL, "fail should return NULL");
	  continue;
	}
	else
	  fail_unless(rc == FC_SUCCESS, "failed to do binary operation");

	// test result meta data
	rc = fc_getVariableInfo(newSeqVar[0], &temp_numDataPoint, &temp_numComp, 
			      &temp_assoc, &temp_mathType, &temp_dataType);
	fail_unless(rc == FC_SUCCESS, "failed to get new variable info");
	fail_unless(temp_numDataPoint == numDataPoint, "mismatch of numDataPoint");
	fail_unless(temp_numComp      == numComp,      "mismatch of numComp");
	fail_unless(temp_assoc        == assoc,        "mismatch of assoc");
	fail_unless(temp_mathType     == mathType,     "mismatch of mathType");
	fail_unless(temp_dataType     == FC_DT_DOUBLE, "mismatch of dataType");
	rc = fc_getSequenceFromSeqVariable(glb_numStep, newSeqVar, &temp_seq);
	fail_unless(rc == FC_SUCCESS, "abort: failed to get new sequence");
	fail_unless(FC_HANDLE_EQUIV(temp_seq, glb_sequence), "mismatch of sequence");

	// test result big data
	for (m = 0; m < glb_numStep; m++) {	
	  rc = fc_getVariableDataPtr(seqVar[m], &startData);
	  fail_unless(rc == FC_SUCCESS, "abort: failed to get data pointer");
	  rc = fc_getVariableDataPtr(newSeqVar[m], &resultData);
	  fail_unless(rc == FC_SUCCESS, "failed to get data ptr from new var");
	  for (n = 0; n < numDataPoint*numComp; n++) {
	    double start1, start2;
	    start1 = dblConst[ n%numComp ];
	    switch (dataType) {
	    case FC_DT_INT:     start2 =    ((int*)startData)[n];  break;
	    case FC_DT_FLOAT:   start2 =  ((float*)startData)[n];  break;
	    case FC_DT_DOUBLE:  start2 = ((double*)startData)[n];  break;
	    default: ;
	    }
	    
	    fail_unless(FC_DBL_EQUIV(((double*)resultData)[n],
				     dblFuncs[k](start1, start2)),
			"mismatch of result data");
	  }
	}
	rc = fc_deleteSeqVariable(glb_numStep, newSeqVar);
	fail_unless(rc == FC_SUCCESS, "failed to delete new seqvar");
	free(newSeqVar);
      }
    }
  }


  // --- things should still work with global vars
  rc = fc_seqVarBinaryFunctionConstSeqVar(constNumComps[0], dblConst, 
					  glb_numStep, glb_glbSeqVar, //gsv component is 1
                                          dblFuncs[0], newSeqVarName, &newSeqVar);
  fail_unless(rc == FC_SUCCESS, "should work with global var");
  fail_unless(fc_isVariableGlobal(newSeqVar[0]), "should be global");
  for (i = 0; i < glb_numStep; i++) {
    rc = fc_getVariableDataPtr(glb_glbSeqVar[i], &startData);
    fail_unless(rc == FC_SUCCESS, "failed to get data ptr from old var");
    
    rc = fc_getVariableDataPtr(newSeqVar[i], &resultData);
    fail_unless(rc == FC_SUCCESS, "failed to get data ptr from new var");
    // just checking first value
    fail_unless(((double*)resultData)[0] == dblFuncs[0](((double*)dblConst)[0],
                                                        ((double*)startData)[0]),
                "mismatch of result data");
  }
  rc = fc_deleteSeqVariable(glb_numStep, newSeqVar);
  fail_unless(rc == FC_SUCCESS, "failed to delete new var");
  free(newSeqVar);


  // --- test bad inputs
  seqVar = glb_dataTypeSeqVars[3];
  for (i = 0; i < glb_numStep; i++)
    badSeqVar[i] = glb_dataTypeSeqVars[3][i];
  badSeqVar[1] = badVar;
  rc = fc_seqVarBinaryFunctionConstSeqVar(constNumComps[0], dblConst, glb_numStep, NULL,  
					  dblFuncs[0], newSeqVarName, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail with bad input var");
  fail_unless(newSeqVar == NULL,"fail should return NULL");
  rc = fc_seqVarBinaryFunctionConstSeqVar(constNumComps[0], dblConst, glb_numStep, badSeqVar,  
					  dblFuncs[0], newSeqVarName, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail with bad input var");
  fail_unless(newSeqVar == NULL,"fail should return NULL");
  rc = fc_seqVarBinaryFunctionConstSeqVar(constNumComps[0], dblConst, glb_numStep, seqVar,  
					  NULL, newSeqVarName, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail with bad op string");
  fail_unless(newSeqVar == NULL,"fail should return NULL");
  rc = fc_seqVarBinaryFunctionConstSeqVar(constNumComps[0], NULL, glb_numStep, seqVar,  
					  dblFuncs[0], newSeqVarName, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail with bad input var");
  fail_unless(newSeqVar == NULL,"fail should return NULL");
  rc = fc_seqVarBinaryFunctionConstSeqVar(0, dblConst, glb_numStep, seqVar, 
					  dblFuncs[0], newSeqVarName, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail with bad input var");
  fail_unless(newSeqVar == NULL,"fail should return NULL");

  rc = fc_seqVarBinaryFunctionConstSeqVar(constNumComps[0], dblConst, glb_numStep, seqVar, 
					  dblFuncs[0], NULL, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail with NULL name");
  fail_unless(newSeqVar == NULL,"fail should return NULL");
  rc = fc_seqVarBinaryFunctionConstSeqVar(constNumComps[0], dblConst, glb_numStep, seqVar,  
					  dblFuncs[0], newSeqVarName, NULL);
  fail_unless(rc != FC_SUCCESS, "should fail if NULL new var");


  //Sanity check to make sure the above calls COULD have worked with right params
  rc = fc_seqVarBinaryFunctionConstSeqVar(constNumComps[0], dblConst, glb_numStep, seqVar,  
					   dblFuncs[0], newSeqVarName, &newSeqVar);
  fail_unless(rc == FC_SUCCESS, "copy should have worked?");
  rc = fc_deleteSeqVariable(glb_numStep, newSeqVar);
  fail_unless(rc == FC_SUCCESS, "failed to delete new var");
  free(newSeqVar);

}
END_TEST




//CDU
// test fc_seqVarBinaryFunctionSeqVarSeqConst()
START_TEST(binarySeqVarSeqConst)
{
  FC_ReturnCode rc;
  int i, j, k, m, n, p;
  FC_Sequence temp_seq;
  FC_Variable *seqVar,  *newSeqVar, badSeqVar[3], badVar = { 999, 999 };
  int numDataPoint, numComp;
  int temp_numDataPoint, temp_numComp;
  char newSeqVarName[20] = "new seq var";
  FC_AssociationType assoc, temp_assoc;
  FC_MathType mathType, temp_mathType;
  FC_DataType dataType, temp_dataType;
  void* startData, *resultData;
  int numFunc = 3;
  double (*dblFuncs[3])(double, double) = { pow, subtraction, magnitude };
  int numConst = 4*2; // four data types * whether scalar or vector
  // index by step, then by component
  char charConst[3*3] = { 'a', 'z', '$',
                          'b', 'y', '@', 
                          'c', 'x', '@' };
  int intConst[3*3] = { 20, -2,   6, 
                        21, -3,  16, 
                        22, -4, 160 };
  float fltConst[3*3] = { 3.5,  -99.99,   1.001, 
                          4.5, -100.99,  10.01, 
                          5.5, -101.99, 100.1 };
  double dblConst[3*3] = { 1.2,  0.00001,   -29.9,
                           2.2, -1.00001,  -299,
                           3.2, -2.00001, -2990 };
  int constNumComps[4*2] = { 1, 3, 1, 3, 1, 3, 1, 3 };
  void *constValues[4*2] = { &charConst, &charConst, &intConst, &intConst,
			     &fltConst, &fltConst, &dblConst, &dblConst };  
  FC_DataType constTypes[4*2] = { FC_DT_CHAR, FC_DT_CHAR, FC_DT_INT,
                                  FC_DT_INT, FC_DT_FLOAT, FC_DT_FLOAT,
                                  FC_DT_DOUBLE, FC_DT_DOUBLE };




  // --- use global variables created in the fixture

  // --- test all ops on all combinations of all var types
  for (i = 0; i < glb_numDataType + 1; i++) {

    if      (i <  glb_numDataType)  seqVar = glb_dataTypeSeqVars[i];
    else                            seqVar = glb_vectorSeqVar;

    for (j = 0; j < numConst; j++) {

      // get the start meta data
      rc = fc_getVariableInfo(seqVar[0], &numDataPoint, &numComp, &assoc, 
			      &mathType,  &dataType);

      fail_unless(rc == FC_SUCCESS, "abort: failed to get var info");

      // loop over ops
      for (k = 0; k < numFunc; k++) {

	// do it
	rc = fc_seqVarBinaryFunctionSeqVarSeqConst(glb_numStep, seqVar, 
						   constNumComps[j], constValues[j], constTypes[j],
						   dblFuncs[k], newSeqVarName, &newSeqVar);

	if (dataType == FC_DT_CHAR || numComp != constNumComps[j]) {
	  fail_unless(rc != FC_SUCCESS, "should fail to do char data or if "
		      "vars don't have same numComponent");
	  fail_unless(newSeqVar == NULL, "fail should return NULL");
	  continue;
	}
	else
	  fail_unless(rc == FC_SUCCESS, "failed to do binary operation");

	// test result meta data
	rc = fc_getVariableInfo(newSeqVar[0], &temp_numDataPoint, &temp_numComp, 
			      &temp_assoc, &temp_mathType, &temp_dataType);
	fail_unless(rc == FC_SUCCESS, "failed to get new variable info");
	fail_unless(temp_numDataPoint == numDataPoint, "mismatch of numDataPoint");
	fail_unless(temp_numComp      == numComp,      "mismatch of numComp");
	fail_unless(temp_assoc        == assoc,        "mismatch of assoc");
	fail_unless(temp_mathType     == mathType,     "mismatch of mathType");
	fail_unless(temp_dataType     == FC_DT_DOUBLE, "mismatch of dataType");
	rc = fc_getSequenceFromSeqVariable(glb_numStep, newSeqVar, &temp_seq);
	fail_unless(rc == FC_SUCCESS, "abort: failed to get new sequence");
	fail_unless(FC_HANDLE_EQUIV(temp_seq, glb_sequence), "mismatch of sequence");
	
	// test result big data
	for (m = 0; m < glb_numStep; m++) {	
	  rc = fc_getVariableDataPtr(seqVar[m], &startData);
	  fail_unless(rc == FC_SUCCESS, "abort: failed to get data pointer");
	  rc = fc_getVariableDataPtr(newSeqVar[m], &resultData);
	  fail_unless(rc == FC_SUCCESS, "failed to get data ptr from new var");

	  for(n=0; n<numDataPoint; n++)
	    for(p=0; p<numComp; p++){

	      int vidx, cidx;
	      double start1, start2;
	      vidx = n*numComp + p;
	      cidx = m*numComp + p;

	      switch (dataType) {
	      case FC_DT_CHAR:    start1 = (double) ((char *) startData)[vidx];  break;
	      case FC_DT_INT:     start1 = (double) ((int*)   startData)[vidx];  break;
	      case FC_DT_FLOAT:   start1 = (double) ((float*) startData)[vidx];  break;
	      case FC_DT_DOUBLE:  start1 = (double) ((double*)startData)[vidx];  break;
	      default: fail_unless(0, "Unknown test case?");
	      }
	      switch (constTypes[j]){
	      case FC_DT_CHAR:    start2 = (double) ((char *)  constValues[j])[cidx]; break;
	      case FC_DT_INT:     start2 = (double) ((int *)   constValues[j])[cidx]; break;
	      case FC_DT_FLOAT:   start2 = (double) ((float *) constValues[j])[cidx]; break;
	      case FC_DT_DOUBLE:  start2 = (double) ((double *)constValues[j])[cidx]; break;
	      default: fail_unless(0, "Unknown test case?");
	      }
	    
	      fail_unless(FC_DBL_EQUIV(((double*)resultData)[vidx],
				       dblFuncs[k](start1, start2)),
			  "mismatch of result data");	    
	  }
	}
	rc = fc_deleteSeqVariable(glb_numStep, newSeqVar);
	fail_unless(rc == FC_SUCCESS, "failed to delete new seqvar");
	free(newSeqVar);
      }
    }
  }


  // --- things should still work with global vars
  rc = fc_seqVarBinaryFunctionSeqVarSeqConst(glb_numStep, glb_glbSeqVar, //gsv component is 1
					     1, dblConst, FC_DT_DOUBLE, 
					     dblFuncs[0], newSeqVarName, &newSeqVar);
  fail_unless(rc == FC_SUCCESS, "should work with global var");
  fail_unless(fc_isVariableGlobal(newSeqVar[0]), "should be global");
  for (i = 0; i < glb_numStep; i++) {
    rc = fc_getVariableDataPtr(glb_glbSeqVar[i], &startData);
    fail_unless(rc == FC_SUCCESS, "failed to get data ptr from old var");
    
    rc = fc_getVariableDataPtr(newSeqVar[i], &resultData);
    fail_unless(rc == FC_SUCCESS, "failed to get data ptr from new var");
    // just checking first value
    fail_unless(((double*)resultData)[0] == dblFuncs[0](((double*)startData)[0],
                                                        ((double*)dblConst)[i*1]),
                "mismatch of result data");
  }
  rc = fc_deleteSeqVariable(glb_numStep, newSeqVar);
  fail_unless(rc == FC_SUCCESS, "failed to delete new var");
  free(newSeqVar);


  // --- test bad inputs
  seqVar = glb_dataTypeSeqVars[3];
  for (i = 0; i < glb_numStep; i++)
    badSeqVar[i] = glb_dataTypeSeqVars[3][i];
  badSeqVar[1] = badVar;
  rc = fc_seqVarBinaryFunctionSeqVarSeqConst(glb_numStep, NULL, 1, dblConst, FC_DT_DOUBLE, 
					  dblFuncs[0], newSeqVarName, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail with bad input var");
  fail_unless(newSeqVar == NULL,"fail should return NULL");
  rc = fc_seqVarBinaryFunctionSeqVarSeqConst(glb_numStep, badSeqVar, 1, dblConst, FC_DT_DOUBLE,
				   dblFuncs[0], newSeqVarName, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail with bad input var");
  fail_unless(newSeqVar == NULL,"fail should return NULL");
  rc = fc_seqVarBinaryFunctionSeqVarSeqConst(glb_numStep, seqVar, 1, dblConst, FC_DT_DOUBLE,
					   NULL, newSeqVarName, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail with bad op string");
  fail_unless(newSeqVar == NULL,"fail should return NULL");
  rc = fc_seqVarBinaryFunctionSeqVarSeqConst(glb_numStep, seqVar, 1, NULL, FC_DT_DOUBLE,
				       dblFuncs[0], newSeqVarName, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail with bad input var");
  fail_unless(newSeqVar == NULL,"fail should return NULL");
  rc = fc_seqVarBinaryFunctionSeqVarSeqConst(glb_numStep, seqVar, 0, dblConst, FC_DT_DOUBLE,
				       dblFuncs[0], newSeqVarName, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail with bad input var");
  fail_unless(newSeqVar == NULL,"fail should return NULL");
  rc = fc_seqVarBinaryFunctionSeqVarSeqConst(glb_numStep, seqVar, 1, dblConst, FC_DT_UNKNOWN,
				       dblFuncs[0], newSeqVarName, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail with bad input var");
  fail_unless(newSeqVar == NULL,"fail should return NULL");

  rc = fc_seqVarBinaryFunctionSeqVarSeqConst(glb_numStep, seqVar, 1, dblConst, FC_DT_DOUBLE,
				       dblFuncs[0], NULL, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail with NULL name");
  fail_unless(newSeqVar == NULL,"fail should return NULL");

  rc = fc_seqVarBinaryFunctionSeqVarSeqConst(glb_numStep, seqVar, 0, dblConst, FC_DT_DOUBLE,
				       dblFuncs[0], newSeqVarName, NULL);
  fail_unless(rc != FC_SUCCESS, "should fail if NULL new var");



  //Sanity check to make sure the above calls COULD have worked with right params
  rc = fc_seqVarBinaryFunctionSeqVarSeqConst(glb_numStep, seqVar, 1, dblConst, FC_DT_DOUBLE,
				       dblFuncs[0], newSeqVarName, &newSeqVar);
  fail_unless(rc == FC_SUCCESS, "copy should have worked?");
  rc = fc_deleteSeqVariable(glb_numStep, newSeqVar);
  fail_unless(rc == FC_SUCCESS, "failed to delete new var");
  free(newSeqVar);
}
END_TEST


//CDU
// test fc_seqVarBinaryFunctionSeqConstSeqVar()
START_TEST(binarySeqConstSeqVar)
{
  FC_ReturnCode rc;
  int i, j, k, m, n, p;
  double r1,r2;
  FC_Sequence temp_seq;
  FC_Variable *seqVar,  *newSeqVar, badSeqVar[3], badVar = { 999, 999 };
  int numDataPoint, numComp;
  int temp_numDataPoint, temp_numComp;
  char newSeqVarName[20] = "new seq var";
  FC_AssociationType assoc, temp_assoc;
  FC_MathType mathType, temp_mathType;
  FC_DataType dataType, temp_dataType;
  void* startData, *resultData;
  int numFunc = 3;
  double (*dblFuncs[3])(double, double) = { pow, subtraction, magnitude };
  int numConst = 4*2; // four data types * whether scalar or vector
  // index by step, then by component
  char charConst[3*3] = { 'a', 'z', '$',
                          'b', 'y', '@', 
                          'c', 'x', '@' };
  int intConst[3*3] = { 20, -2,   6, 
                        21, -3,  16, 
                        22, -4, 160 };
  float fltConst[3*3] = { 3.5,  -99.99,   1.001, 
                          4.5, -100.99,  10.01, 
                          5.5, -101.99, 100.1 };
  double dblConst[3*3] = { 1.2,  0.00001,   -29.9,
                           2.2, -1.00001,  -299,
                           3.2, -2.00001, -2990 };
  int constNumComps[4*2] = { 1, 3, 1, 3, 1, 3, 1, 3 };
  void *constValues[4*2] = { &charConst, &charConst, &intConst, &intConst,
			     &fltConst, &fltConst, &dblConst, &dblConst };  
  FC_DataType constTypes[4*2] = { FC_DT_CHAR, FC_DT_CHAR, FC_DT_INT,
                                  FC_DT_INT, FC_DT_FLOAT, FC_DT_FLOAT,
                                  FC_DT_DOUBLE, FC_DT_DOUBLE };




  // --- use global variables created in the fixture

  // --- test all ops on all combinations of all var types
  for (i = 0; i < glb_numDataType + 1; i++) {

    if      (i <  glb_numDataType)  seqVar = glb_dataTypeSeqVars[i];
    else                            seqVar = glb_vectorSeqVar;

    for (j = 0; j < numConst; j++) {

      // get the start meta data
      rc = fc_getVariableInfo(seqVar[0], &numDataPoint, &numComp, &assoc, 
			      &mathType,  &dataType);

      fail_unless(rc == FC_SUCCESS, "abort: failed to get var info");

      // loop over ops
      for (k = 0; k < numFunc; k++) {


	// do it
	rc = fc_seqVarBinaryFunctionSeqConstSeqVar(glb_numStep, 
						   constNumComps[j], constValues[j], constTypes[j],
						   seqVar,
						   dblFuncs[k], newSeqVarName, &newSeqVar);

	if (dataType == FC_DT_CHAR || numComp != constNumComps[j]) {
	  fail_unless(rc != FC_SUCCESS, "should fail to do char data or if "
		      "vars don't have same numComponent");
	  fail_unless(newSeqVar == NULL, "fail should return NULL");
	  continue;
	}
	else
	  fail_unless(rc == FC_SUCCESS, "failed to do binary operation");

	// test result meta data
	rc = fc_getVariableInfo(newSeqVar[0], &temp_numDataPoint, &temp_numComp, 
			      &temp_assoc, &temp_mathType, &temp_dataType);
	fail_unless(rc == FC_SUCCESS, "failed to get new variable info");
	fail_unless(temp_numDataPoint == numDataPoint, "mismatch of numDataPoint");
	fail_unless(temp_numComp      == numComp,      "mismatch of numComp");
	fail_unless(temp_assoc        == assoc,        "mismatch of assoc");
	fail_unless(temp_mathType     == mathType,     "mismatch of mathType");
	fail_unless(temp_dataType     == FC_DT_DOUBLE, "mismatch of dataType");
	rc = fc_getSequenceFromSeqVariable(glb_numStep, newSeqVar, &temp_seq);
	fail_unless(rc == FC_SUCCESS, "abort: failed to get new sequence");
	fail_unless(FC_HANDLE_EQUIV(temp_seq, glb_sequence), "mismatch of sequence");
	
	// test result big data
	for (m = 0; m < glb_numStep; m++) {	
	  rc = fc_getVariableDataPtr(seqVar[m], &startData);
	  fail_unless(rc == FC_SUCCESS, "abort: failed to get data pointer");
	  rc = fc_getVariableDataPtr(newSeqVar[m], &resultData);
	  fail_unless(rc == FC_SUCCESS, "failed to get data ptr from new var");

	  for(n=0; n<numDataPoint; n++)
	    for(p=0; p<numComp; p++){

	      int vidx, cidx;
	      double start1, start2;
	      vidx = n*numComp + p;
	      cidx = m*numComp + p;

	      switch (constTypes[j]){
	      case FC_DT_CHAR:    start1 = (double) ((char *)  constValues[j])[cidx]; break;
	      case FC_DT_INT:     start1 = (double) ((int *)   constValues[j])[cidx]; break;
	      case FC_DT_FLOAT:   start1 = (double) ((float *) constValues[j])[cidx]; break;
	      case FC_DT_DOUBLE:  start1 = (double) ((double *)constValues[j])[cidx]; break;
	      default: fail_unless(0, "Unknown test case?");
	      }
	      switch (dataType) {
	      case FC_DT_CHAR:    start2 = (double) ((char *) startData)[vidx];  break;
	      case FC_DT_INT:     start2 = (double) ((int*)   startData)[vidx];  break;
	      case FC_DT_FLOAT:   start2 = (double) ((float*) startData)[vidx];  break;
	      case FC_DT_DOUBLE:  start2 = (double) ((double*)startData)[vidx];  break;
	      default: fail_unless(0, "Unknown test case?");
	      }
	    
	      //There are some NaNs in here, have to do a special check
	      r1 = dblFuncs[k](start1, start2);
	      r2 = ((double*)resultData)[vidx];
	      fail_unless( *(int*)&r1 == *(int *)&r2,"mismatch of result data");
	      
	  }
	}
	rc = fc_deleteSeqVariable(glb_numStep, newSeqVar);
	fail_unless(rc == FC_SUCCESS, "failed to delete new seqvar");
	free(newSeqVar);
      }
    }
  }


  // --- things should still work with global vars
  rc = fc_seqVarBinaryFunctionSeqConstSeqVar(glb_numStep, 
					     1, dblConst, FC_DT_DOUBLE,
					     glb_glbSeqVar, //gsv component is 1 
					     dblFuncs[0], newSeqVarName, &newSeqVar);
  fail_unless(rc == FC_SUCCESS, "should work with global var");
  fail_unless(fc_isVariableGlobal(newSeqVar[0]), "should be global");
  for (i = 0; i < glb_numStep; i++) {
    rc = fc_getVariableDataPtr(glb_glbSeqVar[i], &startData);
    fail_unless(rc == FC_SUCCESS, "failed to get data ptr from old var");
    
    rc = fc_getVariableDataPtr(newSeqVar[i], &resultData);
    fail_unless(rc == FC_SUCCESS, "failed to get data ptr from new var");

    r1 = ((double*)resultData)[0];
    r2 = dblFuncs[0](((double*)dblConst)[i*1], ((double*)startData)[0]);

    // just checking first value : note, last op makes a nan so look at ints in comparison
    fail_unless( *(int*)&r1==*(int*)&r2,"mismatch of result data");
  }
  rc = fc_deleteSeqVariable(glb_numStep, newSeqVar);
  fail_unless(rc == FC_SUCCESS, "failed to delete new var");
  free(newSeqVar);


  // --- test bad inputs
  seqVar = glb_dataTypeSeqVars[3];
  for (i = 0; i < glb_numStep; i++)
    badSeqVar[i] = glb_dataTypeSeqVars[3][i];
  badSeqVar[1] = badVar;
  rc = fc_seqVarBinaryFunctionSeqConstSeqVar(glb_numStep, 1, dblConst, FC_DT_DOUBLE, NULL,  
					  dblFuncs[0], newSeqVarName, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail with bad input var");
  fail_unless(newSeqVar == NULL,"fail should return NULL");
  rc = fc_seqVarBinaryFunctionSeqConstSeqVar(glb_numStep, 1, dblConst, FC_DT_DOUBLE, badSeqVar,
				   dblFuncs[0], newSeqVarName, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail with bad input var");
  fail_unless(newSeqVar == NULL,"fail should return NULL");
  rc = fc_seqVarBinaryFunctionSeqConstSeqVar(glb_numStep, 1, dblConst, FC_DT_DOUBLE, seqVar,
					   NULL, newSeqVarName, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail with bad op string");
  fail_unless(newSeqVar == NULL,"fail should return NULL");
  rc = fc_seqVarBinaryFunctionSeqConstSeqVar(glb_numStep, 1, NULL, FC_DT_DOUBLE, seqVar,
				       dblFuncs[0], newSeqVarName, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail with bad input var");
  fail_unless(newSeqVar == NULL,"fail should return NULL");
  rc = fc_seqVarBinaryFunctionSeqConstSeqVar(glb_numStep, 0, dblConst, FC_DT_DOUBLE, seqVar,
				       dblFuncs[0], newSeqVarName, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail with bad input var");
  fail_unless(newSeqVar == NULL,"fail should return NULL");
  rc = fc_seqVarBinaryFunctionSeqConstSeqVar(glb_numStep, 1, dblConst, FC_DT_UNKNOWN, seqVar, 
				       dblFuncs[0], newSeqVarName, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail with bad input var");
  fail_unless(newSeqVar == NULL,"fail should return NULL");

  rc = fc_seqVarBinaryFunctionSeqConstSeqVar(glb_numStep, 1, dblConst, FC_DT_DOUBLE, seqVar,
				       dblFuncs[0], NULL, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail with NULL name");
  fail_unless(newSeqVar == NULL,"fail should return NULL");

  rc = fc_seqVarBinaryFunctionSeqConstSeqVar(glb_numStep, 0, dblConst, FC_DT_DOUBLE, seqVar, 
				       dblFuncs[0], newSeqVarName, NULL);
  fail_unless(rc != FC_SUCCESS, "should fail if NULL new var");



  //Sanity check to make sure the above calls COULD have worked with right params
  rc = fc_seqVarBinaryFunctionSeqConstSeqVar(glb_numStep, 1, dblConst, FC_DT_DOUBLE, seqVar, 
				       dblFuncs[0], newSeqVarName, &newSeqVar);
  fail_unless(rc == FC_SUCCESS, "copy should have worked?");
  rc = fc_deleteSeqVariable(glb_numStep, newSeqVar);
  fail_unless(rc == FC_SUCCESS, "failed to delete new var");
  free(newSeqVar);



}
END_TEST



// test fc_seqVarBinaryFunctionSeqVarVar()
START_TEST(binarySeqVarVar)
{
  FC_ReturnCode rc;
  int i, j, k, m, n;
  FC_Sequence temp_seq;
  FC_Variable *seqVar1, var2, *newSeqVar, badSeqVar[3], badVar = { 999, 999 };
  int numDataPoint1, numComp1, numDataPoint2, numComp2;
  int temp_numDataPoint, temp_numComp;
  char newSeqVarName[20] = "new seq var";
  FC_AssociationType assoc1, assoc2, temp_assoc;
  FC_MathType mathType1, mathType2, temp_mathType;
  FC_DataType dataType1, dataType2, temp_dataType;
  void* startData1, *startData2, *resultData;
  int numFunc = 3;
  double (*dblFuncs[3])(double, double) = { pow, subtraction, magnitude };

  // --- use global variables created in the fixture

  // --- test all ops on all combinations of all var types
  for (i = 0; i < glb_numDataType + 2; i++) {

    if (i < glb_numDataType)
      seqVar1 = glb_dataTypeSeqVars[i];
    else if (i == glb_numDataType)
      seqVar1 = glb_vectorSeqVar;
    else 
      seqVar1 = glb_vectorSeqVar;

    for (j = 0; j < glb_numDataType + 2; j++) {
      if (j < glb_numDataType)
	var2 = glb_dataTypeSeqVars[j][0];
      else if (j == glb_numDataType)
	var2 = glb_vectorSeqVar[0];
      else
	var2 = glb_elemSeqVar[0];

      // get the start meta data
      rc = fc_getVariableInfo(seqVar1[0], &numDataPoint1, &numComp1, &assoc1, 
			      &mathType1,  &dataType1);
      fail_unless(rc == FC_SUCCESS, "abort: failed to get var info");
      rc = fc_getVariableInfo(var2, &numDataPoint2, &numComp2, &assoc2, 
			      &mathType2, &dataType2);
      fail_unless(rc == FC_SUCCESS, "abort: failed to get var info");

      // loop over ops
      for (k = 0; k < numFunc; k++) {

	// do it
	rc = fc_seqVarBinaryFunctionSeqVarVar(glb_numStep, seqVar1, var2,
					dblFuncs[k], newSeqVarName, &newSeqVar);
	if (dataType1 == FC_DT_CHAR || dataType2 == FC_DT_CHAR ||
	    numDataPoint1 != numDataPoint2 || numComp1 != numComp2) {
	  fail_unless(rc != FC_SUCCESS, "should fail to do char data or if "
		      "vars don't have some numDataPoint or don't have same "
		      "numComp");
	  fail_unless(newSeqVar == NULL, "fail should return NULL");
	  continue;
	}
	else
	  fail_unless(rc == FC_SUCCESS, "failed to do binary operation");

	// test result meta data
	rc = fc_getVariableInfo(newSeqVar[0], &temp_numDataPoint, &temp_numComp, 
			      &temp_assoc, &temp_mathType, &temp_dataType);
	fail_unless(rc == FC_SUCCESS, "failed to get new variable info");
	fail_unless(temp_numDataPoint == numDataPoint1, 
		    "mismatch of numDataPoint");
	fail_unless(temp_numComp == numComp1, "mismatch of numComp");
	fail_unless(temp_assoc == assoc1, "mismatch of assoc");
	fail_unless(temp_mathType == mathType1, "mismatch of mathType");
	fail_unless(temp_dataType == FC_DT_DOUBLE, "mismatch of dataType");
	rc = fc_getSequenceFromSeqVariable(glb_numStep, newSeqVar, &temp_seq);
	fail_unless(rc == FC_SUCCESS, "abort: failed to get new sequence");
	fail_unless(FC_HANDLE_EQUIV(temp_seq, glb_sequence), "mismatch of sequence");

	// test result big data
	for (m = 0; m < glb_numStep; m++) {	
	  rc = fc_getVariableDataPtr(seqVar1[m], &startData1);
	  fail_unless(rc == FC_SUCCESS, "abort: failed to get data pointer");
	  rc = fc_getVariableDataPtr(var2, &startData2);
	  fail_unless(rc == FC_SUCCESS, "abort: failed to get data pointer");
	  rc = fc_getVariableDataPtr(newSeqVar[m], &resultData);
	  fail_unless(rc == FC_SUCCESS, "failed to get data ptr from new var");
	  for (n = 0; n < numDataPoint1*numComp1; n++) {
	    double start1, start2;
	    switch (dataType1) {
	    case FC_DT_INT:     start1 =    ((int*)startData1)[n];  break;
	    case FC_DT_FLOAT:   start1 =  ((float*)startData1)[n];  break;
	    case FC_DT_DOUBLE:  start1 = ((double*)startData1)[n];  break;
	    default: ;
	    }
	    switch (dataType2) {
	    case FC_DT_INT:     start2 =    ((int*)startData2)[n];  break;
	    case FC_DT_FLOAT:   start2 =  ((float*)startData2)[n];  break;
	    case FC_DT_DOUBLE:  start2 = ((double*)startData2)[n];  break;
	    default: ;
	    }
	    fail_unless(FC_DBL_EQUIV(((double*)resultData)[n],
				     dblFuncs[k](start1, start2)),
			"mismatch of result data");
	  }
	}
	rc = fc_deleteSeqVariable(glb_numStep, newSeqVar);
	fail_unless(rc == FC_SUCCESS, "failed to delete new seqvar");
	free(newSeqVar);
      }
    }
  }
 
  // --- things should still work with global vars
  rc = fc_seqVarBinaryFunctionSeqVarVar(glb_numStep, glb_glbSeqVar, glb_glbSeqVar[0],
                                        dblFuncs[0], newSeqVarName, &newSeqVar);
  fail_unless(rc == FC_SUCCESS, "should work with global var");
  fail_unless(fc_isVariableGlobal(newSeqVar[0]), "should be global");
  for (i = 0; i < glb_numStep; i++) {
    rc = fc_getVariableDataPtr(glb_glbSeqVar[i], &startData1);
    fail_unless(rc == FC_SUCCESS, "failed to get data ptr from old var");
    rc = fc_getVariableDataPtr(glb_glbSeqVar[0], &startData2);
    fail_unless(rc == FC_SUCCESS, "failed to get data ptr from old var");
    rc = fc_getVariableDataPtr(newSeqVar[i], &resultData);
    fail_unless(rc == FC_SUCCESS, "failed to get data ptr from new var");
    // just checking first value
    fail_unless(((double*)resultData)[0] == dblFuncs[0](((double*)startData1)[0],
                                                      ((double*)startData2)[0]),
                "mismatch of result data");
  }
  rc = fc_deleteSeqVariable(glb_numStep, newSeqVar);
  fail_unless(rc == FC_SUCCESS, "failed to delete new var");
  free(newSeqVar);

  // --- things should fail with mixed global non-global vars
  seqVar1 = glb_dataTypeSeqVars[3];
  var2 = glb_dataTypeSeqVars[3][0];
  rc = fc_seqVarBinaryFunctionSeqVarVar(glb_numStep, glb_glbSeqVar, var2,
                                        dblFuncs[0], newSeqVarName, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail with mix of global/nonglobal");
  fail_unless(newSeqVar == NULL, "fail should return NULL");
  rc = fc_seqVarBinaryFunctionSeqVarVar(glb_numStep, seqVar1, glb_glbSeqVar[0], 
                                        dblFuncs[0], newSeqVarName, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail with mix of global/nonglobal");
  fail_unless(newSeqVar == NULL, "fail should return NULL");

  // --- test bad inputs
  seqVar1 = glb_dataTypeSeqVars[3];
  var2 = glb_dataTypeSeqVars[3][0];
  for (i = 0; i < glb_numStep; i++)
    badSeqVar[i] = glb_dataTypeSeqVars[3][i];
  badSeqVar[1] = badVar;
  rc = fc_seqVarBinaryFunctionSeqVarVar(glb_numStep, NULL, var2, 
                                      dblFuncs[0], newSeqVarName, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail with bad input var");
  fail_unless(newSeqVar == NULL,"fail should return NULL");
  rc = fc_seqVarBinaryFunctionSeqVarVar(glb_numStep, badSeqVar, var2,
				   dblFuncs[0], newSeqVarName, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail with bad input var");
  fail_unless(newSeqVar == NULL,"fail should return NULL");
  rc = fc_seqVarBinaryFunctionSeqVarVar(glb_numStep, seqVar1, var2, 
					   NULL, newSeqVarName, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail with bad op string");
  fail_unless(newSeqVar == NULL,"fail should return NULL");
  rc = fc_seqVarBinaryFunctionSeqVarVar(glb_numStep, seqVar1, badVar,
				       dblFuncs[0], newSeqVarName, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail with bad input var");
  fail_unless(newSeqVar == NULL,"fail should return NULL");
  rc = fc_seqVarBinaryFunctionSeqVarVar(glb_numStep, seqVar1, var2,
					   dblFuncs[0], NULL, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail with NULL name");
  fail_unless(newSeqVar == NULL,"fail should return NULL");
  rc = fc_seqVarBinaryFunctionSeqVarVar(glb_numStep, seqVar1, var2,
					   dblFuncs[0], newSeqVarName, NULL);
  fail_unless(rc != FC_SUCCESS, "should fail if NULL new var");
}
END_TEST

// test fc_seqVarBinaryFunctionVarSeqVar()
START_TEST(binaryVarSeqVar)
{
  FC_ReturnCode rc;
  int i, j, k, m, n;
  FC_Sequence temp_seq;
  FC_Variable var1, *seqVar2, *newSeqVar, badSeqVar[3], badVar = { 999, 999 };
  int numDataPoint1, numComp1, numDataPoint2, numComp2;
  int temp_numDataPoint, temp_numComp;
  char newSeqVarName[20] = "new seq var";
  FC_AssociationType assoc1, assoc2, temp_assoc;
  FC_MathType mathType1, mathType2, temp_mathType;
  FC_DataType dataType1, dataType2, temp_dataType;
  void* startData1, *startData2, *resultData;
  int numFunc = 3;
  double (*dblFuncs[3])(double, double) = { pow, subtraction, magnitude };

  // --- use global variables created in the fixture

  // --- test all ops on all combinations of all var types
  for (i = 0; i < glb_numDataType + 2; i++) {

    if (i < glb_numDataType)
      var1 = glb_dataTypeSeqVars[i][0];
    else if (i == glb_numDataType)
      var1 = glb_vectorSeqVar[0];
    else 
      var1 = glb_vectorSeqVar[0];

    for (j = 0; j < glb_numDataType + 2; j++) {
      if (j < glb_numDataType)
	seqVar2 = glb_dataTypeSeqVars[j];
      else if (j == glb_numDataType)
	seqVar2 = glb_vectorSeqVar;
      else
	seqVar2 = glb_elemSeqVar;

      // get the start meta data
      rc = fc_getVariableInfo(var1, &numDataPoint1, &numComp1, &assoc1, 
			      &mathType1,  &dataType1);
      fail_unless(rc == FC_SUCCESS, "abort: failed to get var info");
      rc = fc_getVariableInfo(seqVar2[0], &numDataPoint2, &numComp2, &assoc2, 
			      &mathType2, &dataType2);
      fail_unless(rc == FC_SUCCESS, "abort: failed to get var info");

      // loop over ops
      for (k = 0; k < numFunc; k++) {

	// do it
	rc = fc_seqVarBinaryFunctionVarSeqVar(var1, glb_numStep, seqVar2,
					dblFuncs[k], newSeqVarName, &newSeqVar);
	if (dataType1 == FC_DT_CHAR || dataType2 == FC_DT_CHAR ||
	    numDataPoint1 != numDataPoint2 || numComp1 != numComp2) {
	  fail_unless(rc != FC_SUCCESS, "should fail to do char data or if "
		      "vars don't have some numDataPoint or don't have same "
		      "numComp");
	  fail_unless(newSeqVar == NULL, "fail should return NULL");
	  continue;
	}
	else
	  fail_unless(rc == FC_SUCCESS, "failed to do binary operation");

	// test result meta data
	rc = fc_getVariableInfo(newSeqVar[0], &temp_numDataPoint, &temp_numComp, 
			      &temp_assoc, &temp_mathType, &temp_dataType);
	fail_unless(rc == FC_SUCCESS, "failed to get new variable info");
	fail_unless(temp_numDataPoint == numDataPoint1, 
		    "mismatch of numDataPoint");
	fail_unless(temp_numComp == numComp1, "mismatch of numComp");
	fail_unless(temp_assoc == assoc1, "mismatch of assoc");
	fail_unless(temp_mathType == mathType1, "mismatch of mathType");
	fail_unless(temp_dataType == FC_DT_DOUBLE, "mismatch of dataType");
	rc = fc_getSequenceFromSeqVariable(glb_numStep, newSeqVar, &temp_seq);
	fail_unless(rc == FC_SUCCESS, "abort: failed to get new sequence");
	fail_unless(FC_HANDLE_EQUIV(temp_seq, glb_sequence), "mismatch of sequence");

	// test result big data
	for (m = 0; m < glb_numStep; m++) {	
	  rc = fc_getVariableDataPtr(var1, &startData1);
	  fail_unless(rc == FC_SUCCESS, "abort: failed to get data pointer");
	  rc = fc_getVariableDataPtr(seqVar2[m], &startData2);
	  fail_unless(rc == FC_SUCCESS, "abort: failed to get data pointer");
	  rc = fc_getVariableDataPtr(newSeqVar[m], &resultData);
	  fail_unless(rc == FC_SUCCESS, "failed to get data ptr from new var");
	  for (n = 0; n < numDataPoint1*numComp1; n++) {
	    double start1, start2;
	    switch (dataType1) {
	    case FC_DT_INT:     start1 =    ((int*)startData1)[n];  break;
	    case FC_DT_FLOAT:   start1 =  ((float*)startData1)[n];  break;
	    case FC_DT_DOUBLE:  start1 = ((double*)startData1)[n];  break;
	    default: ;
	    }
	    switch (dataType2) {
	    case FC_DT_INT:     start2 =    ((int*)startData2)[n];  break;
	    case FC_DT_FLOAT:   start2 =  ((float*)startData2)[n];  break;
	    case FC_DT_DOUBLE:  start2 = ((double*)startData2)[n];  break;
	    default: ;
	    }
	    fail_unless(FC_DBL_EQUIV(((double*)resultData)[n],
				     dblFuncs[k](start1, start2)),
			"mismatch of result data");
	  }
	}
	rc = fc_deleteSeqVariable(glb_numStep, newSeqVar);
	fail_unless(rc == FC_SUCCESS, "failed to delete new seqvar");
	free(newSeqVar);
      }
    }
  }
 
  // --- things should still work with global vars
  rc = fc_seqVarBinaryFunctionVarSeqVar(glb_glbSeqVar[0], glb_numStep, glb_glbSeqVar, 
                                        dblFuncs[0], newSeqVarName, &newSeqVar);
  fail_unless(rc == FC_SUCCESS, "should work with global var");
  fail_unless(fc_isVariableGlobal(newSeqVar[0]), "should be global");
  for (i = 0; i < glb_numStep; i++) {
    rc = fc_getVariableDataPtr(glb_glbSeqVar[0], &startData1);
    fail_unless(rc == FC_SUCCESS, "failed to get data ptr from old var");
    rc = fc_getVariableDataPtr(glb_glbSeqVar[i], &startData2);
    fail_unless(rc == FC_SUCCESS, "failed to get data ptr from old var");
    rc = fc_getVariableDataPtr(newSeqVar[i], &resultData);
    fail_unless(rc == FC_SUCCESS, "failed to get data ptr from new var");
    // just checking first value
    fail_unless(((double*)resultData)[0] == dblFuncs[0](((double*)startData1)[0],
                                                      ((double*)startData2)[0]),
                "mismatch of result data");
  }
  rc = fc_deleteSeqVariable(glb_numStep, newSeqVar);
  fail_unless(rc == FC_SUCCESS, "failed to delete new var");
  free(newSeqVar);

  // --- things should fail with mixed global non-global vars
  var1 = glb_dataTypeSeqVars[3][0];
  seqVar2 = glb_dataTypeSeqVars[3];
  rc = fc_seqVarBinaryFunctionVarSeqVar(var1, glb_numStep, glb_glbSeqVar,
                            dblFuncs[0], newSeqVarName, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail with mix of global/nonglobal");
  fail_unless(newSeqVar == NULL, "fail should return NULL");
  rc = fc_seqVarBinaryFunctionVarSeqVar(glb_glbSeqVar[0], glb_numStep, seqVar2, 
                            dblFuncs[0], newSeqVarName, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail with mix of global/nonglobal");
  fail_unless(newSeqVar == NULL, "fail should return NULL");

  // --- test bad inputs
  var1 = glb_dataTypeSeqVars[3][0];
  seqVar2 = glb_dataTypeSeqVars[3];
  for (i = 0; i < glb_numStep; i++)
    badSeqVar[i] = glb_dataTypeSeqVars[3][i];
  badSeqVar[1] = badVar;
  rc = fc_seqVarBinaryFunctionVarSeqVar(badVar, glb_numStep, seqVar2,
				   dblFuncs[0], newSeqVarName, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail with bad input var");
  fail_unless(newSeqVar == NULL,"fail should return NULL");
  rc = fc_seqVarBinaryFunctionVarSeqVar(var1, glb_numStep, seqVar2,  
					   NULL, newSeqVarName, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail with bad op string");
  fail_unless(newSeqVar == NULL,"fail should return NULL");
  rc = fc_seqVarBinaryFunctionVarSeqVar(var1, glb_numStep, NULL,
				   dblFuncs[0], newSeqVarName, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail with bad input var");
  fail_unless(newSeqVar == NULL,"fail should return NULL");
  rc = fc_seqVarBinaryFunctionVarSeqVar(var1, glb_numStep, badSeqVar,
				       dblFuncs[0], newSeqVarName, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail with bad input var");
  fail_unless(newSeqVar == NULL,"fail should return NULL");
  rc = fc_seqVarBinaryFunctionVarSeqVar(var1, glb_numStep, seqVar2,
					   dblFuncs[0], NULL, &newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail with NULL name");
  fail_unless(newSeqVar == NULL,"fail should return NULL");
  rc = fc_seqVarBinaryFunctionVarSeqVar(var1, glb_numStep, seqVar2,
					   dblFuncs[0], newSeqVarName, NULL);
  fail_unless(rc != FC_SUCCESS, "should fail if NULL new var");
}
END_TEST

// test magnitude calculation
START_TEST(other_magnitude)
{
  FC_ReturnCode rc;
  int i, j;
  FC_Variable vector_var, scalar_var, mag_var, badVar = { 999, 999 };
  FC_Variable glb_vector_var;
  char mag_name[] = "magnitude", *temp_name;
  int numDataPoint, numComp, temp_numDataPoint, temp_numComp;
  FC_AssociationType assoc, temp_assoc;
  FC_MathType mathtype, temp_mathtype;
  FC_DataType datatype, temp_datatype;
  double* data, *temp_data, *magnitudes;

  // --- use global mesh and variables created in the fixture

  rc = fc_getMeshCoordsAsVariable(glb_mesh, &vector_var);
  fail_unless(rc == FC_SUCCESS, "abort: failed to get mesh coord as var");
  rc = fc_getVariableInfo(vector_var, &numDataPoint, &numComp, &assoc, 
			  &mathtype, &datatype);
  fail_unless(rc == FC_SUCCESS, "abort: failed to get vector var info");
  rc = fc_getVariableDataPtr(vector_var, (void**)&data);
  fail_unless(rc == FC_SUCCESS, "abort: failed to get vector var data");
  scalar_var = glb_dataTypeVars[3];
  rc = fc_createGlobalVariable(glb_dataset, "global vector var",
                               &glb_vector_var);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create global vector var");
  rc = fc_setVariableData(glb_vector_var, 1, numComp, FC_AT_WHOLE_DATASET, 
			  FC_MT_VECTOR, FC_DT_DOUBLE, data);
  fail_unless(rc == FC_SUCCESS, "abort: failed to set glb vector var data");

  // calculate magnitude data
  magnitudes = malloc(sizeof(double)*numDataPoint);
  for (i = 0; i < numDataPoint; i++) {
    magnitudes[i] = 0.;
    for (j = 0; j < numComp; j++)
      magnitudes[i] += data[i*numComp+j]*data[i*numComp+j];
    magnitudes[i] = sqrt(magnitudes[i]);
  }

  // use fclib to get magnitude data & compare
  rc = fc_createMagnitudeVariable(vector_var, mag_name, &mag_var);
  fail_unless(rc == FC_SUCCESS, "failed to get magnitude variable");
  rc = fc_getVariableName(mag_var, &temp_name);
  fail_unless(rc == FC_SUCCESS, "failed to get mag var name");
  fail_unless(!strcmp(mag_name, temp_name), "mismatch of name");
  free(temp_name);
  rc = fc_getVariableInfo(mag_var, &temp_numDataPoint, &temp_numComp,
			    &temp_assoc, &temp_mathtype, &temp_datatype);
  fail_unless(rc == FC_SUCCESS, "failed to describe mag var");
  fail_unless(temp_numDataPoint == numDataPoint, "numDataPoint mismatch");
  fail_unless(temp_numComp == 1, "numComp should = 1");
  fail_unless(temp_assoc == assoc, "assoc mismatch");
  fail_unless(temp_mathtype == FC_MT_SCALAR, "mathtype should be scalar");
  fail_unless(temp_datatype == FC_DT_DOUBLE, "mathtype should be double");
  rc = fc_getVariableDataPtr(mag_var, (void**)&temp_data);
  fail_unless(rc == FC_SUCCESS, "failed to get data from mag var");
  for (i = 0; i < numDataPoint; i++)
    fail_unless(magnitudes[i] == temp_data[i], "mismatch in data");
  fc_deleteVariable(mag_var);

  // --- should work on global var
  rc = fc_createMagnitudeVariable(glb_vector_var, mag_name, &mag_var);
  fail_unless(rc == FC_SUCCESS, "should work with global var");
  fail_unless(fc_isVariableGlobal(mag_var), "should be global");
  rc = fc_getVariableDataPtr(mag_var, (void**)&temp_data);
  fail_unless(rc == FC_SUCCESS, "failed to get data from mag var");
  // just checking first value
  fail_unless(FC_DBL_EQUIV(magnitudes[0], temp_data[0]), "mismatch in data");

  // --- test errors
  // Should not work on scalar var
  rc = fc_createMagnitudeVariable(scalar_var, mag_name, &mag_var);
  fail_unless(rc != FC_SUCCESS, "should fail if scalar var");
  fail_unless(FC_HANDLE_EQUIV(mag_var, FC_NULL_VARIABLE),
	      "fail should return NULL");
  
  // fc_createMagnitudeVariable() --- bad args
  rc = fc_createMagnitudeVariable(badVar, mag_name, &mag_var);
  fail_unless(rc != FC_SUCCESS, "should fail with bad var");
  fail_unless(FC_HANDLE_EQUIV(mag_var, FC_NULL_VARIABLE),
	      "fail should return NULL");
  rc = fc_createMagnitudeVariable(vector_var, NULL, &mag_var);
  fail_unless(rc != FC_SUCCESS, "should fail if NULL name");
  fail_unless(FC_HANDLE_EQUIV(mag_var, FC_NULL_VARIABLE),
	      "fail should return NULL");
  rc = fc_createMagnitudeVariable(vector_var, mag_name, NULL);
  fail_unless(rc != FC_SUCCESS, "should fail if NULL name");

  // cleanup
  free(magnitudes);
}
END_TEST	

// test decomposing vector into normal & tangent
/// \FIX? might want to test on scalar var (should still work)
START_TEST(decompose_vector)
{
  FC_ReturnCode rc;
  int i, j, k, i_type;
  int numDataType = 3;
  FC_DataType dataTypes[3] = { FC_DT_DOUBLE, FC_DT_INT, FC_DT_FLOAT };
  FC_Variable vector_vars[3], glb_vector_var, normals, tangents;
  FC_Variable scalar_var, elem_var, badVar = { 999, 999 };
  int numDataPoint, numComp, temp_numDataPoint, temp_numComp;
  FC_AssociationType assoc, temp_assoc;
  FC_MathType mathtype, temp_mathtype;
  FC_DataType datatype, temp_datatype;
  double *temp_data, *coords_data, *scalar_data, *normal_data, *tangent_data;
  double axes[3][3] = { { 1, }, { 0, 1 }, { 0, 0, 1 } };

  // get a double vector var & some scalar vars
  rc = fc_getMeshCoordsAsVariable(glb_mesh, &vector_vars[0]);
  fail_unless(rc == FC_SUCCESS, "abort: failed to get mesh coord as var");
  rc = fc_getVariableInfo(vector_vars[0], &numDataPoint, &numComp, &assoc, 
			  &mathtype, &datatype);
  fail_unless(datatype == FC_DT_DOUBLE, "abort: var should be double");
  rc = fc_getVariableDataPtr(vector_vars[0], (void**)&coords_data);
  fail_unless(rc == FC_SUCCESS, "abort: failed to get coords var data");
  scalar_var = glb_dataTypeVars[3]; 
  rc = fc_getVariableDataPtr(scalar_var, (void**)&scalar_data);
  fail_unless(rc == FC_SUCCESS, "abort: failed to get scalar var data");
  elem_var = glb_elemVar;
  rc = fc_createGlobalVariable(glb_dataset, "global vector var",
                               &glb_vector_var);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create global vector var");
  rc = fc_setVariableData(glb_vector_var, 1, numComp, FC_AT_WHOLE_DATASET, 
			  FC_MT_VECTOR, FC_DT_DOUBLE, coords_data);
  fail_unless(rc == FC_SUCCESS, "abort: failed to set glb vector var data");

  // make the int & float vector_vars
  rc = fc_getVariableDataAsDataType(vector_vars[0], FC_DT_INT, (void**)&temp_data);
  fail_unless(rc == FC_SUCCESS, "abort: failed to get data as ints");
  rc = fc_createVariable(glb_mesh, "ints", &vector_vars[1]);
  rc = fc_setVariableDataPtr(vector_vars[1], numDataPoint, numComp, assoc,
			     mathtype, dataTypes[1], temp_data);
  fail_unless(rc == FC_SUCCESS, "abort: failed to set data on ints var");
  rc = fc_getVariableDataAsDataType(vector_vars[0], FC_DT_FLOAT, (void**)&temp_data);
  fail_unless(rc == FC_SUCCESS, "abort: failed to get data as flts");
  rc = fc_createVariable(glb_mesh, "flts", &vector_vars[2]);
  rc = fc_setVariableDataPtr(vector_vars[2], numDataPoint, numComp, assoc,
			     mathtype, dataTypes[2], temp_data);
  fail_unless(rc == FC_SUCCESS, "abort: failed to set data on ints var");

  // --- test fc_createNormalTangentVariables() with good data

  // get tangents and normals of coords to reference vectors (the axes)
  // also tests against different dimensionality of reference vectors
  for (i = 0; i < 3; i++) {
    for (i_type = 0; i_type < numDataType; i_type++) {
      // do it
      rc = fc_createNormalTangentVariables(vector_vars[i_type], i+1, axes[i], &normals, 
					   &tangents);
      fail_unless(rc == FC_SUCCESS, "failed to decompose by axis vector");
      
      // check metadata
      //rc = fc_getVariableName(normals, temp_name);
      rc = fc_getVariableInfo(normals, &temp_numDataPoint, &temp_numComp,
			      &temp_assoc, &temp_mathtype, &temp_datatype);
      fail_unless(rc == FC_SUCCESS, "failed to describe normal var");
      fail_unless(temp_numDataPoint == numDataPoint, "numDataPoint mismatch");
      fail_unless(temp_numComp == 1, "numComp should = 1");
      fail_unless(temp_assoc == assoc, "assoc mismatch");
      fail_unless(temp_mathtype == FC_MT_SCALAR, "mathtype should be scalar");
      fail_unless(temp_datatype == FC_DT_DOUBLE, "mathtype should be double");
      rc = fc_getVariableInfo(tangents, &temp_numDataPoint, &temp_numComp,
			      &temp_assoc, &temp_mathtype, &temp_datatype);
      fail_unless(rc == FC_SUCCESS, "failed to describe tangent var");
      fail_unless(temp_numDataPoint == numDataPoint, "numDataPoint mismatch");
      fail_unless(temp_numComp == 1, "numComp should = 1");
      fail_unless(temp_assoc == assoc, "assoc mismatch");
      fail_unless(temp_mathtype == FC_MT_SCALAR, "mathtype should be scalar");
      fail_unless(temp_datatype == FC_DT_DOUBLE, "mathtype should be double");
      
      // check data -- normal should be distance to point along the axes,
      // tangent should be distance to point in plan perpendicular to axis
      rc = fc_getVariableDataPtr(normals, (void**)&normal_data);
      fail_unless(rc == FC_SUCCESS, "failed to describe mag var");
      rc = fc_getVariableDataPtr(tangents, (void**)&tangent_data);
      fail_unless(rc == FC_SUCCESS, "failed to describe mag var");
      for (j = 0; j < numDataPoint; j++) {
	double temp_sq_sum = 0.0;
	if (i_type == 0) {
	  fail_unless(normal_data[j] == coords_data[j*numComp+i], 
		      "mismatch of normal data");
	  for (k = 0; k < 3; k++) 
	    if (k != i)
	      temp_sq_sum += coords_data[j*numComp+k]*coords_data[j*numComp+k];
	  fail_unless(FC_VALUE_EQUIV(tangent_data[j], sqrt(temp_sq_sum),
				     40*DBL_EPSILON, DBL_MIN), "mismatch of tangent data");
	}
	else if (i_type == 1) {
	  fail_unless(normal_data[j] == (int)coords_data[j*numComp+i],
		      "msimatch of normal data");
	  for (k =0; k < 3; k++) 
	    if (k != i)
	      temp_sq_sum += ((int)coords_data[j*numComp+k])*((int)coords_data[j*numComp+k]);
	  fail_unless(FC_DBL_EQUIV(tangent_data[j], sqrt(temp_sq_sum)),
		      "mismatch of tangent data");
	}
	else if (i_type == 2) {
	  fail_unless(FC_FLT_EQUIV(normal_data[j], coords_data[j*numComp+i]),
		      "mismatch of normal data");
	  for (k = 0; k < 3; k++)
	    if (k != i)
	      temp_sq_sum += coords_data[j*numComp+k]*coords_data[j*numComp+k];
	  fail_unless(FC_VALUE_EQUIV(tangent_data[j], sqrt(temp_sq_sum),
				     FLT_EPSILON, FLT_MIN), "mismatch of tangent data");
	}
      }
    }
  }

  // get tangents and normals of coords to a scalar input vector
  // also tests against different dimensionality of reference vectors
  for (i = 0; i < 3; i++) {
    // do it
    rc = fc_createNormalTangentVariables(scalar_var, i+1, axes[i], &normals, 
					 &tangents);
    fail_unless(rc == FC_SUCCESS, "failed to decompose by axis vector");

    // not repeating metadata check
    rc = fc_getVariableDataPtr(normals, (void**)&normal_data);
    fail_unless(rc == FC_SUCCESS, "failed to describe mag var");
    rc = fc_getVariableDataPtr(tangents, (void**)&tangent_data);
    fail_unless(rc == FC_SUCCESS, "failed to describe mag var");

    // check data -- normal should be same as 
    // tangent should be distanc to point in plan perpendicular to axis
    for (j = 0; j < numDataPoint; j++) {
      if (i == 0) {
	fail_unless(normal_data[j] == scalar_data[j], "mismatch of normal data");
	fail_unless(tangent_data[j] == 0., "msimatch of tangent data");
      }
      else {
	fail_unless(normal_data[j] == 0, "mismatch of normal data");
	fail_unless(tangent_data[j] == scalar_data[j], 
		    "mismatch of tangent data");
      }
    }
  }

  // --- test fc_createNormalTangentVariables2() with good data

  for (i_type = 0; i_type < numDataType; i_type++) {
    // get component-wise normals & tangents between two vectors - to itself
    rc = fc_createNormalTangentVariables2(vector_vars[i_type], vector_vars[i_type], &normals,
					  &tangents);
    fail_unless(rc == FC_SUCCESS, "failed to decompose by axis vector");
    
    // check metadata
    //rc = fc_getVariableName()
    rc = fc_getVariableInfo(normals, &temp_numDataPoint, &temp_numComp,
			    &temp_assoc, &temp_mathtype, &temp_datatype);
    fail_unless(rc == FC_SUCCESS, "failed to describe normal var");
    fail_unless(temp_numDataPoint == numDataPoint, "numDataPoint mismatch");
    fail_unless(temp_numComp == 1, "numComp should = 1");
    fail_unless(temp_assoc == assoc, "assoc mismatch");
    fail_unless(temp_mathtype == FC_MT_SCALAR, "mathtype should be scalar");
    fail_unless(temp_datatype == FC_DT_DOUBLE, "mathtype should be double");
    rc = fc_getVariableInfo(tangents, &temp_numDataPoint, &temp_numComp,
			    &temp_assoc, &temp_mathtype, &temp_datatype);
    fail_unless(rc == FC_SUCCESS, "failed to describe tangent var");
    fail_unless(temp_numDataPoint == numDataPoint, "numDataPoint mismatch");
    fail_unless(temp_numComp == 1, "numComp should = 1");
    fail_unless(temp_assoc == assoc, "assoc mismatch");
    fail_unless(temp_mathtype == FC_MT_SCALAR, "mathtype should be scalar");
    fail_unless(temp_datatype == FC_DT_DOUBLE, "mathtype should be double");
   
    // check data -- normal should be distance to the point, tangent = 0
    for (j = 0; j < numDataPoint; j++) {
      double temp_sq_sum = 0.0;
      rc = fc_getVariableDataPtr(normals, (void**)&normal_data);
      fail_unless(rc == FC_SUCCESS, "failed to describe mag var");
      rc = fc_getVariableDataPtr(tangents, (void**)&tangent_data);
      fail_unless(rc == FC_SUCCESS, "failed to describe mag var");
      
      if (i_type == 0) {
	for (k = 0; k < 3; k++) 
	  temp_sq_sum += coords_data[j*numComp+k]*coords_data[j*numComp+k];
	fail_unless(FC_VALUE_EQUIV(normal_data[j], sqrt(temp_sq_sum),
				   2*DBL_EPSILON, DBL_MIN), "mismatch of normal data");

      }
      else if (i_type == 1) {
	for (k = 0; k < 3; k++)
	temp_sq_sum += ((int)coords_data[j*numComp+k])*((int)coords_data[j*numComp+k]);
	fail_unless(FC_VALUE_EQUIV(normal_data[j], sqrt(temp_sq_sum),
				   2*DBL_EPSILON, DBL_MIN),
		    "mismatch of normal data");
      }
      else if (i_type == 2) {
	for (k = 0; k < 3; k++)
	  temp_sq_sum += coords_data[j*numComp+k]*coords_data[j*numComp+k];
	fail_unless(FC_FLT_EQUIV(normal_data[j], sqrt(temp_sq_sum)),
		    "mismatch of normal data");
      }
      fail_unless(tangent_data[j] == 0., "msimatch of tangent data");
    }
  }

  // --- should work on global vars
  rc = fc_createNormalTangentVariables(glb_vector_var, 3, axes[0], &normals,
                                     &tangents);
  fail_unless(rc == FC_SUCCESS, "should work for global vars");
  fail_unless(fc_isVariableGlobal(normals), "normals should be global");
  fail_unless(fc_isVariableGlobal(tangents), "tangents should be global");
  rc = fc_getVariableDataPtr(normals, (void**)&normal_data);
  fail_unless(rc == FC_SUCCESS, "failed to get normal data");
  rc = fc_getVariableDataPtr(tangents, (void**)&tangent_data);
  fail_unless(rc == FC_SUCCESS, "failed to get tangent data");
  fail_unless(normal_data[0] == coords_data[0], "mismatch of normal data");
  {
    double temp_sq_sum = 0.0;
    for (i = 2; i < 3; i++)
      temp_sq_sum += coords_data[i]*coords_data[i];
    fail_unless(FC_DBL_EQUIV(tangent_data[0], sqrt(temp_sq_sum)),
                "mismatch of tangent data");
  }
  rc = fc_createNormalTangentVariables2(glb_vector_var, glb_vector_var, &normals,
                                       &tangents);
  fail_unless(rc == FC_SUCCESS, "should work for global vars");
  fail_unless(fc_isVariableGlobal(normals), "normals should be global");
  fail_unless(fc_isVariableGlobal(tangents), "tangents should be global");
  rc = fc_getVariableDataPtr(normals, (void**)&normal_data);
  fail_unless(rc == FC_SUCCESS, "failed to get normal data");
  rc = fc_getVariableDataPtr(tangents, (void**)&tangent_data);
  fail_unless(rc == FC_SUCCESS, "failed to get tangent data");
  {
    double temp_sq_sum = 0.0;
    for (i = 0; i < 3; i++)
      temp_sq_sum += coords_data[i]*coords_data[i];
    fail_unless(FC_DBL_EQUIV(normal_data[0], sqrt(temp_sq_sum)), 
                "mismatch of normal data");
  }
  fail_unless(FC_DBL_EQUIV(tangent_data[0], 0), "mismatch of tangent data");

  // --- test errors

  // fc_createNormalTangentVariables2() - should fail if assoc mismatch
  rc = fc_createNormalTangentVariables2(vector_vars[0], elem_var, &normals, &tangents);
  fail_unless(rc != FC_SUCCESS, "should fail if vector_vars don't have same assoc");
  fail_unless(FC_HANDLE_EQUIV(normals, FC_NULL_VARIABLE) &&
	      FC_HANDLE_EQUIV(tangents, FC_NULL_VARIABLE),
	      "fail should return nulls");

  // fc_createNormalTangentVariables()--- bad args
  rc = fc_createNormalTangentVariables(badVar, 3, axes[2], &normals, 
				       &tangents);
  fail_unless(rc != FC_SUCCESS, "should fail with bar var");
  fail_unless(FC_HANDLE_EQUIV(normals, FC_NULL_VARIABLE) &&
	      FC_HANDLE_EQUIV(tangents, FC_NULL_VARIABLE),
	      "fail should return nulls");
  rc = fc_createNormalTangentVariables(vector_vars[0], 0, axes[2], &normals, &tangents);
  fail_unless(rc != FC_SUCCESS, "should fail if dim < 1");
  fail_unless(FC_HANDLE_EQUIV(normals, FC_NULL_VARIABLE) &&
	      FC_HANDLE_EQUIV(tangents, FC_NULL_VARIABLE),
	      "fail should return nulls");
  rc = fc_createNormalTangentVariables(vector_vars[0], 3, NULL, &normals, &tangents);
  fail_unless(rc != FC_SUCCESS, "should fail if null ref_vector");
  fail_unless(FC_HANDLE_EQUIV(normals, FC_NULL_VARIABLE) &&
	      FC_HANDLE_EQUIV(tangents, FC_NULL_VARIABLE),
	      "fail should return nulls");
  rc = fc_createNormalTangentVariables(vector_vars[0], 3, axes[2], NULL, &tangents);
  fail_unless(rc != FC_SUCCESS, "should fail if null normal var");
  fail_unless(FC_HANDLE_EQUIV(tangents, FC_NULL_VARIABLE),
	      "fail should return nulls");
  rc = fc_createNormalTangentVariables(vector_vars[0], 3, axes[2], &normals, NULL);
  fail_unless(rc != FC_SUCCESS, "should fail if null tangent var");
  fail_unless(FC_HANDLE_EQUIV(normals, FC_NULL_VARIABLE),
	      "fail should return nulls");

  // fc_createNormalTangentVariables2()--- bad args
  rc = fc_createNormalTangentVariables2(badVar, vector_vars[0], &normals, &tangents);
  fail_unless(rc != FC_SUCCESS, "should fail with bad var 1");
  fail_unless(FC_HANDLE_EQUIV(normals, FC_NULL_VARIABLE) &&
	      FC_HANDLE_EQUIV(tangents, FC_NULL_VARIABLE),
	      "fail should return nulls");
  rc = fc_createNormalTangentVariables2(vector_vars[0], badVar, &normals, &tangents);
  fail_unless(rc != FC_SUCCESS, "should fail with bad var 2");
  fail_unless(FC_HANDLE_EQUIV(normals, FC_NULL_VARIABLE) &&
	      FC_HANDLE_EQUIV(tangents, FC_NULL_VARIABLE),
	      "fail should return nulls");
  rc = fc_createNormalTangentVariables2(vector_vars[0], vector_vars[0], NULL, &tangents);
  fail_unless(rc != FC_SUCCESS, "should fail if null normal var");
  fail_unless(FC_HANDLE_EQUIV(tangents, FC_NULL_VARIABLE),
	      "fail should return nulls");
  rc = fc_createNormalTangentVariables2(vector_vars[0], vector_vars[0], &normals, NULL);
  fail_unless(rc != FC_SUCCESS, "should fail if null normal var");
  fail_unless(FC_HANDLE_EQUIV(normals, FC_NULL_VARIABLE),
	      "fail should return nulls");
}
END_TEST	

// test calc of von mises stress and pressure
// Assuming will work for all assocation types and all data types
// (based on simplicity of implementation)
START_TEST(stress_pressure)
{
  FC_ReturnCode rc;
  int i, j;
  FC_Dataset dataset;
  FC_Mesh mesh;
  FC_Variable sigmaVars[6], stressVar, pressureVar, coordsVar;
  FC_Variable badVar = { 99, 99 }, badVars[6];
  int numDim = 3;
  int numDataPoint = 16;
  //{ xx, yy, zz, xy, yz, zx };
  double sigmas[16][6] = {
    // - all zeros
    { 0, 0, 0, 0, 0, 0 },
    // - all equal
    { 1, 1, 1, 1, 1, 1 },
    // - just the diagonals
    { 1, 1, 1, 0, 0, 0 },
    // - just the off diagonals
    { 0, 0, 0, 1, 1, 1 },
    // - just 1 diagonal component
    { 1, 0, 0, 0, 0, 0 },
    { 0, 1, 0, 0, 0, 0 },
    { 0, 0, 1, 0, 0, 0 },
    // - just 1 off diagonal component
    { 0, 0, 0, 1, 0, 0 },
    { 0, 0, 0, 0, 1, 0 },
    { 0, 0, 0, 0, 0, 1 },
    // - compression examples
    { 60.9214057922363, -10652.4248046875, 61.492244720459,
      1.60643982887268, -1.57024359703064, 0.724126398563385 },
    { 222.0478515625, -11164.4267578125, 221.66569519043,
      -0.871581554412842, -0.743675410747528, -0.813497066497803 },
    { -1229.19311523438, -13006.4228515625, -1226.36450195312,
      44.4218330383301,-40.9876098632812, -6.86527395248413  },
    // - tension examples
    { -62.4811553955078, 10664.8291015625, -63.1964683532715,
      -2.1255030632019, 2.08135604858398, 0.225947350263596 },
    { -221.840698242188, 11178.1005859375, -221.220443725586,
      0.876150906085968, 0.854998588562012, 1.31932067871094 },
    { 1234.462890625, 13016.591796875, 1231.99865722656,
      -45.3536911010742, 41.5615577697754, 7.44467353820801 } };
  char* stress_name = "von mises stress";
  char* pressure_name = "pressure";
  double stress[16] =   { 
    0, 3, 0, 3,
    1, 1, 1,    
    sqrt(3), sqrt(3), sqrt(3),
    10713.63242129134, 11386.28379611244, 11779.11553292312,
    10727.66917588874, 11399.63159581044, 11783.84999828698
  };
  double pressure[16] = { 
    0, 1, 1, 0,
    1./3., 1./3., 1./3.,
    0, 0, 0,
    -3510.00371805827, -3573.57107035319, -5153.99348958334,
    3513.05049260457, 3578.34648132324, 5161.017781575514  };
  int* dummy_conns;
  double* dummy_coords;
  char* temp_name;
  int temp_numDataPoint, temp_numComp;
  FC_AssociationType temp_assoc;
  FC_MathType temp_mathType;
  FC_DataType temp_dataType;
  double *temp_data;

  // Create dataset
  fc_createDataset("dummy", &dataset);

  // Create mesh
  dummy_conns = calloc(numDataPoint, sizeof(int));
  dummy_coords = calloc(numDim*numDataPoint, sizeof(double));
  fc_createMesh(dataset, "dummy", &mesh);
  rc = fc_setMeshCoordsPtr(mesh, numDim, numDataPoint, dummy_coords);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create mesh");
  rc = fc_setMeshElementConnsPtr(mesh, FC_ET_POINT, numDataPoint, dummy_conns);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create mesh");
  rc = fc_getMeshCoordsAsVariable(mesh, &coordsVar);
  fail_unless(rc == FC_SUCCESS, "abort: failed to get coords as var");

  // create input tensor
  for (i = 0; i < 6; i++) {
    temp_data = malloc(numDataPoint*sizeof(double));
    for (j = 0; j < numDataPoint; j++)
      temp_data[j] = sigmas[j][i];
    fc_createVariable(mesh, "stress", &sigmaVars[i]);
    rc = fc_setVariableDataPtr(sigmaVars[i], numDataPoint, 1,
			       FC_AT_ELEMENT, FC_MT_SCALAR,
			       FC_DT_DOUBLE, temp_data);
    fail_unless(rc == FC_SUCCESS, "abort: failed to set stress data");
  }
  
  // test
  rc = fc_createStressPressureVariables(sigmaVars, &stressVar, &pressureVar);
  fail_unless(rc == FC_SUCCESS, "failed to calc new vars");
  for (i = 0; i < 2; i++) {
    FC_Variable thisVar;
    char* thisName;
    double* thisData;
    if (i == 0) {
      thisVar = stressVar;
      thisName = stress_name;
      thisData = stress;
    }
    else {
      thisVar = pressureVar;
      thisName = pressure_name;
      thisData = pressure;
    }
    rc = fc_getVariableName(thisVar, &temp_name);
    fail_unless(!strcmp(temp_name, thisName), "mismatch of name");
    free(temp_name);
    rc = fc_getVariableInfo(thisVar, &temp_numDataPoint, &temp_numComp,
			    &temp_assoc, &temp_mathType, &temp_dataType);
    fail_unless(temp_numDataPoint == numDataPoint, "mismatch of datapoint");
    fail_unless(temp_numComp == 1, "mismatch of numComp");
    fail_unless(temp_assoc = FC_AT_ELEMENT, "mismatch of assoc");
    fail_unless(temp_mathType == FC_MT_SCALAR, "mismatch of mathtype");
    fail_unless(temp_dataType == FC_DT_DOUBLE, "mismatch of datatype");
    rc = fc_getVariableDataPtr(thisVar, (void**)&temp_data);
    fail_unless(rc == FC_SUCCESS, "failed to get data");
    for (j = 0; j < numDataPoint; j++) {
      //printf("i = %d, j = %d   temp = %.20g, good = %.20g\n", 
      //     i, j, temp_data[j], thisData[j]);
      if (thisData[j] == 0)
	fail_unless(fabs(temp_data[j]) < 3*DBL_EPSILON,
		    "mismatch of data");
      else
	fail_unless(FC_VALUE_EQUIV(temp_data[j], thisData[j], 8*DBL_EPSILON,
				   DBL_MIN), "mismatch of data");
    }
  }
  
  // ---test bad args
  for (i = 0; i < 6; i++)
    badVars[i] = sigmaVars[i];
  badVars[2] = badVar;
  rc = fc_createStressPressureVariables(badVars, &stressVar, &pressureVar);
  fail_unless(rc != FC_SUCCESS, "should fail if bad input");
  fail_unless(FC_HANDLE_EQUIV(stressVar, FC_NULL_VARIABLE),
	      "fail should return NULL");
  fail_unless(FC_HANDLE_EQUIV(pressureVar, FC_NULL_VARIABLE),
	      "fail should return NULL");
  badVars[2] = coordsVar;
  rc = fc_createStressPressureVariables(badVars, &stressVar, &pressureVar);
  fail_unless(rc != FC_SUCCESS, "should fail if multicomponent input");
  fail_unless(FC_HANDLE_EQUIV(stressVar, FC_NULL_VARIABLE),
	      "fail should return NULL");
  fail_unless(FC_HANDLE_EQUIV(pressureVar, FC_NULL_VARIABLE),
	      "fail should return NULL");
  rc = fc_createStressPressureVariables(NULL, &stressVar, &pressureVar);
  fail_unless(rc != FC_SUCCESS, "should fail if null args");
  fail_unless(FC_HANDLE_EQUIV(stressVar, FC_NULL_VARIABLE),
	      "fail should return NULL");
  fail_unless(FC_HANDLE_EQUIV(pressureVar, FC_NULL_VARIABLE),
	      "fail should return NULL");
  rc = fc_createStressPressureVariables(sigmaVars, NULL, &pressureVar);
  fail_unless(rc != FC_SUCCESS, "should fail if bad input");
  fail_unless(FC_HANDLE_EQUIV(pressureVar, FC_NULL_VARIABLE),
	      "fail should return NULL");
  rc = fc_createStressPressureVariables(sigmaVars, &stressVar, NULL);
  fail_unless(rc != FC_SUCCESS, "should fail if bad input");
  fail_unless(FC_HANDLE_EQUIV(stressVar, FC_NULL_VARIABLE),
	      "fail should return NULL");

  // all done
  rc = fc_deleteDataset(dataset);
  fail_unless(rc == FC_SUCCESS, "failed to delete dataset at end of tests");
}
END_TEST

// *********************************************
// ***** Populate the Suite with the tests
// *********************************************

Suite *varmath_suite(void)
{
  Suite *suite = suite_create("VarMath");

  TCase *tc_helpers = tcase_create(" - Helpers ");
  TCase *tc_operators = tcase_create(" - Built-in Operators (Basic)");
  TCase *tc_seqoperators = tcase_create(" - Built-in Operators (Seq) ");
  TCase *tc_mixoperators = tcase_create(" - Built-in Operators (Mixed) ");
  TCase *tc_functions = tcase_create(" - Simple Functions (Basic)");
  TCase *tc_seqfunctions = tcase_create(" - Simple Functions (Seq) ");
  TCase *tc_mixfunctions = tcase_create(" - Simple Functions (Mixed) ");
  TCase *tc_varOther = tcase_create(" - Specialized Functions ");
  TCase *tc_varOther2 = tcase_create(" - Specialized Functions2 ");

  // test helpers
  suite_add_tcase(suite, tc_helpers);
  tcase_add_test(tc_helpers, unary_array_ops);
  tcase_add_test(tc_helpers, binary_array_ops);

  // test built-in operators on vars
  suite_add_tcase(suite, tc_operators);
  tcase_add_checked_fixture(tc_operators, varmath_setup_dataset,
			    varmath_teardown_dataset);
  tcase_add_test(tc_operators, opVar);
  tcase_add_test(tc_operators, varOpVar);
  tcase_add_test(tc_operators, varOpConst);
  tcase_add_test(tc_operators, constOpVar);

  // test built-in operators on seq vars
  suite_add_tcase(suite, tc_seqoperators);
  tcase_add_checked_fixture(tc_seqoperators, varmath_setup_seqdataset, 
			    varmath_teardown_seqdataset);
  tcase_add_test(tc_seqoperators, opSeqVar);
  tcase_add_test(tc_seqoperators, seqVarOpSeqVar);
  tcase_add_test(tc_seqoperators, seqVarOpSeqConst);
  tcase_add_test(tc_seqoperators, seqConstOpSeqVar);

  // test built-in operators on mix of vars and seq vars
  suite_add_tcase(suite, tc_mixoperators);
  tcase_add_checked_fixture(tc_mixoperators, varmath_setup_seqdataset,
                            varmath_teardown_seqdataset);
  tcase_add_test(tc_mixoperators, seqVarOpVar);
  tcase_add_test(tc_mixoperators, varOpSeqVar);
  tcase_add_test(tc_mixoperators, seqVarOpConst);
  tcase_add_test(tc_mixoperators, constOpSeqVar);
 
  // test using supplied functions
  suite_add_tcase(suite, tc_functions);
  tcase_add_checked_fixture(tc_functions, varmath_setup_dataset, 
			    varmath_teardown_dataset);
  tcase_add_test(tc_functions, unaryVar);
  tcase_add_test(tc_functions, binaryVarVar);
  tcase_add_test(tc_functions, binaryVarConst);
  tcase_add_test(tc_functions, binaryConstVar);


  // test using supplied functions - seq vars
  suite_add_tcase(suite, tc_seqfunctions);
  tcase_add_checked_fixture(tc_seqfunctions, varmath_setup_seqdataset, 
  		    varmath_teardown_seqdataset);
  tcase_add_test(tc_seqfunctions, unarySeqVar);
  tcase_add_test(tc_seqfunctions, binarySeqVarSeqVar);
  
  tcase_add_test(tc_seqfunctions, binarySeqVarConst); //cdu
  tcase_add_test(tc_seqfunctions, binaryConstSeqVar); //cdu
  tcase_add_test(tc_seqfunctions, binarySeqVarSeqConst); //cdu
  tcase_add_test(tc_seqfunctions, binarySeqConstSeqVar); //cdu
 

  // test using supplied functions - mix of var and seq vars
  suite_add_tcase(suite, tc_mixfunctions);
  tcase_add_checked_fixture(tc_mixfunctions, varmath_setup_seqdataset,
                            varmath_teardown_seqdataset);
  tcase_add_test(tc_mixfunctions, binarySeqVarVar);
  tcase_add_test(tc_mixfunctions, binaryVarSeqVar);


  

  // other, specialized functions
  suite_add_tcase(suite, tc_varOther);
  tcase_add_checked_fixture(tc_varOther, varmath_setup_dataset, 
			    varmath_teardown_dataset);
  tcase_add_test(tc_varOther, other_magnitude);
  tcase_add_test(tc_varOther, decompose_vector);

  // other, specialized functions - not using global dataset
  suite_add_tcase(suite, tc_varOther2);
  tcase_add_checked_fixture(tc_varOther2, varmath_setup, varmath_teardown);
  tcase_add_test(tc_varOther2, stress_pressure);


  return suite;
}

