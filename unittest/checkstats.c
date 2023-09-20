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
 * \file checkstats.c
 * \brief Unit testing of \ref Statistics module.
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/unittest/checkstats.c,v $
 * $Revision: 1.43 $ 
 * $Date: 2006/10/19 03:14:53 $
 */

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <check.h>
#include "fc.h"
#include "seriesP.h"
#include "series.h"
#include "checkall.h"

// **** global variables

static int global_min_id = 0;
static int global_max_id = 3;
static double global_min = -13;
static double global_max = 99;
static int global_numMember = 4;
static int global_numAssoc = 5, global_numDataType = 3;
static FC_AssociationType global_assocs[5] = { FC_AT_VERTEX, FC_AT_EDGE,
					       FC_AT_FACE, FC_AT_ELEMENT,
					       FC_AT_WHOLE_MESH };
static FC_DataType global_dataTypes[3] = { FC_DT_INT, FC_DT_FLOAT, 
					   FC_DT_DOUBLE };
static FC_Dataset global_dataset;
static FC_Mesh global_mesh;
static FC_Variable global_vars[5][3];
static FC_Variable global_char_var;
static FC_Variable global_vector_var;
static FC_Subset global_subsets[5];
static FC_Subset global_empty_subset;

//helper fucntion specific for testing
//make from single data point seqvar, make it an N data point seq variable
//same vals over all the data points except first one will be all 1's
static FC_ReturnCode ag_multiDataPointSeqVariable(int orignumStep,
						FC_Variable* orig_seq_var,
						FC_Sequence seq,
						int numDataPoint, 
						FC_AssociationType assoc,
						char* seq_var_name,
						int* numStep,
						FC_Variable** seqvar){

  int orignumDataPoint, orignumComponent;
  FC_AssociationType origassoc;
  FC_MathType origmathtype;
  FC_DataType origdatatype;
  int i,j,k, var_dim;
  void* data;

  FC_ReturnCode rc;

  rc = fc_getVariableInfo(orig_seq_var[0],&orignumDataPoint,
			  &orignumComponent,
			  &origassoc, &origmathtype, &origdatatype);

  rc = fc_createSeqVariable(global_mesh,seq,seq_var_name,
			    numStep, seqvar);

  if (rc != FC_SUCCESS){
    return rc;
  }

  var_dim = numDataPoint*orignumComponent;
  for (i=0; i < orignumStep; i++){
    rc = fc_getVariableDataPtr(orig_seq_var[i],&data);
    if (rc != FC_SUCCESS){
      return rc;
    }
    //neanderthal
    switch(origdatatype){
    case FC_DT_INT:
      {
	int *tempdata = (int*) malloc(var_dim*sizeof(int));
	for (j = 0; j < numDataPoint; j++){
	  for (k = 0; k < orignumComponent; k++){
	    if (j == 0){
	      tempdata[orignumComponent*j+k] = 1;
	    }else {
	      tempdata[orignumComponent*j+k] = ((int*)data)[k];
	    }
	  }
	}
	rc = fc_setVariableData((*seqvar)[i],numDataPoint,orignumComponent,
				assoc,origmathtype,
				origdatatype,tempdata);
	free(tempdata);
      }
    break;
  case FC_DT_FLOAT:
    {
      float *tempdata = (float*) malloc(var_dim*sizeof(float));
      for (j = 0; j < numDataPoint; j++){
	for (k = 0; k < orignumComponent; k++){
	  if (j == 0){
	    tempdata[orignumComponent*j+k] = 1.0;
	  }else {
	    tempdata[orignumComponent*j+k] = ((float*)data)[k];
	  }
	}
      }
      rc = fc_setVariableData((*seqvar)[i],numDataPoint,orignumComponent,
			      assoc,origmathtype,
			      origdatatype,tempdata);
      free(tempdata);
    }
    break;
  case FC_DT_DOUBLE:
    {
      double *tempdata = (double*) malloc(var_dim*sizeof(double));
      for (j = 0; j < numDataPoint; j++){
	for (k = 0; k < orignumComponent; k++){
	  if (j == 0){
	    tempdata[orignumComponent*j+k] = 1.0;
	  }else {
	    tempdata[orignumComponent*j+k] = ((double*)data)[k];
	  }
	}
      }
      rc = fc_setVariableData((*seqvar)[i],numDataPoint,orignumComponent,
			      assoc,origmathtype,
			      origdatatype,tempdata);
      free(tempdata);
    }
    break;
  case FC_DT_CHAR:
    {
      char *tempdata = (char*) malloc(var_dim*sizeof(char));
      for (j = 0; j < numDataPoint; j++){
	for (k = 0; k < orignumComponent; k++){
	  tempdata[orignumComponent*j+k] = ((char*)data)[k];
	}
      }
      rc = fc_setVariableData((*seqvar)[i],numDataPoint,orignumComponent,
			      assoc,origmathtype,
			      origdatatype,tempdata);
      free(tempdata);
    }
    break;
  default:
    //invalid type
    fc_deleteSeqVariable(*numStep,*seqvar);
    free(*seqvar);
    return FC_INPUT_ERROR;
    break;
    }
    data = NULL;
  }

  return FC_SUCCESS;
}



// **** test fixtures - create some variables to test

static void stats_setup(void) {
  FC_ReturnCode rc;
  int i, j;
  char var_name[200], subset_name[200];
  double lowerCoords[3] = { 0., 0., 0. };
  double upperCoords[3] = { 1., 1., 1. };
  int numMember, memberIDs[4] = { global_min_id, 1, global_max_id, 10 };
  int numEntity, maxNumEntity = 144;
  int int_data[144];
  float float_data[144];
  double double_data[144];
  char char_data[144];
  void* data_p[3] = { int_data, float_data, double_data };

  if (isForking) {
    fc_setLibraryVerbosity(fc_messages);
    rc = fc_initLibrary();
    fail_unless(rc == FC_SUCCESS, "test aborted: failed to init library");
  }

  // setup data  
  for (i = 0; i < maxNumEntity; i++) {
    int_data[i] = 0;
    float_data[i] = 0.;
    double_data[i] = 0.;
    char_data[i] = 'a';
  }
  int_data[global_min_id] = global_min;
  int_data[global_max_id] = global_max;
  float_data[global_min_id] = global_min;
  float_data[global_max_id] = global_max;
  double_data[global_min_id] = global_min;
  double_data[global_max_id] = global_max;

  // create test dataset
  fc_createDataset("temp.xxx", &global_dataset);

  // create 3x3 hex mesh
  fc_createSimpleHexMesh(global_dataset, "simple hex mesh", 3, 3, 3,
			 lowerCoords, upperCoords, &global_mesh);

  // create variables
  for (i = 0; i < global_numAssoc; i++) {
    for (j = 0; j < global_numDataType; j++) {
      fc_getMeshNumEntity(global_mesh, global_assocs[i], &numEntity);
      sprintf(var_name, "%s-%s", fc_getAssociationTypeText(global_assocs[i]),
              fc_getDataTypeText(global_dataTypes[j]));
      fc_createVariable(global_mesh, var_name, &global_vars[i][j]);
      rc = fc_setVariableData(global_vars[i][j], numEntity, 1,
                              global_assocs[i], FC_MT_SCALAR, 
                              global_dataTypes[j], data_p[j]);
      fail_unless(rc == FC_SUCCESS, "abort: failed to create variable");
    }
  }
  fc_getMeshNumEntity(global_mesh, FC_AT_VERTEX, &numEntity);
  fc_createVariable(global_mesh, "char var", &global_char_var);
  rc = fc_setVariableData(global_char_var, numEntity, 1, FC_AT_VERTEX,
                          FC_MT_SCALAR, FC_DT_CHAR, char_data);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create char variable");
  fc_createVariable(global_mesh, "vector var", &global_vector_var);
  rc = fc_setVariableData(global_vector_var, numEntity, 2, FC_AT_VERTEX,
                          FC_MT_VECTOR, FC_DT_DOUBLE, double_data);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create vector variable");

  // create subsets
  for (i = 0; i < global_numAssoc; i++) {
    if (global_assocs[i] == FC_AT_WHOLE_MESH) {
      numMember = 1;
      sprintf(subset_name, "the %s",
              fc_getAssociationTypeText(global_assocs[i]));
    }
    else {
      numMember = global_numMember;
      sprintf(subset_name, "first two %s",
              fc_getAssociationTypeText(global_assocs[i]));
    }
    fc_createSubset(global_mesh, subset_name, global_assocs[i], 
                    &global_subsets[i]);
    rc = fc_addArrayMembersToSubset(global_subsets[i], numMember, memberIDs);
    fail_unless(rc == FC_SUCCESS, "abort: failed to create subset");
  }
  rc = fc_createSubset(global_mesh, "empty subset", FC_AT_VERTEX, 
                       &global_empty_subset);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create empty subset");
}

static void stats_teardown(void) {
  FC_ReturnCode rc;

  rc = fc_deleteDataset(global_dataset);
  fail_unless(rc == FC_SUCCESS, "abort: failed to close dataset");

  if (isForking) {
    rc = fc_finalLibrary();
    fail_unless(rc == FC_SUCCESS, "failed to final library at end of test");
  }
}

START_TEST(var_stats)
{
  FC_ReturnCode rc;
  FC_Variable badVar = { 999, 999 };
  int i, j;
  int numEntity;
  double min, max, temp_min, temp_max;
  int min_id, max_id, temp_min_id, temp_max_id;
  double mean, stdev, temp_mean, temp_stdev;
  double sum, temp_sum;

  for (i = 0; i < global_numAssoc; i++) {
    fc_getMeshNumEntity(global_mesh, global_assocs[i], &numEntity);
    if (global_assocs[i] == FC_AT_WHOLE_MESH) {
      min = max = global_min;
      min_id = max_id = global_min_id;
      sum = global_min;
    }
    else {
      min = global_min;
      max = global_max;
      min_id = global_min_id;
      max_id = global_max_id;
      sum = global_min + global_max;
    }
    mean = sum / numEntity;
    if (numEntity == 1)
      stdev = 0.;
    else 
      stdev = sqrt( ( (global_min - mean)*(global_min - mean) + 
                      (global_max - mean)*(global_max - mean) + 
                      (numEntity-2)*mean*mean ) / (numEntity - 1) );
  
    for (j = 0; j < global_numDataType; j++) {
      //printf("i = %d, j = %d\n", i, j);

      // test MinMax
      rc = fc_getVariableMinMax(global_vars[i][j], &temp_min, &temp_min_id, 
                                &temp_max, &temp_max_id);
      fail_unless(rc == FC_SUCCESS, "failed to get min/max");
      fail_unless(temp_min == min, "mismatch for min");
      fail_unless(temp_min_id == min_id, "mismatch for min id");
      fail_unless(temp_max == max, "mismatch for max");
      fail_unless(temp_max_id == max_id, "mismatch for max id");

      // test Sum
      rc = fc_getVariableSum(global_vars[i][j], &temp_sum);
      fail_unless(rc == FC_SUCCESS, "failed to get sum");
      fail_unless(temp_sum == sum, "mismatch for sum");

      // test MeanSdev
      rc = fc_getVariableMeanSdev(global_vars[i][j], &temp_mean,
                                   &temp_stdev);
      fail_unless(rc == FC_SUCCESS, "failed to get mean/stdev");
      fail_unless(temp_mean == mean, "mismatch for mean");
      fail_unless(FC_DBL_EQUIV(temp_stdev, stdev), "mismatch for stdev");

      // --- test special cases

      // fc_getVariableMinMax can take any of the return args as optional
      temp_min = temp_min_id = temp_max = temp_max_id = -1;
      rc = fc_getVariableMinMax(global_vars[i][j], &temp_min, NULL, NULL, 
                                NULL);
      fail_unless(rc == FC_SUCCESS, "failed to get min/max");
      fail_unless(temp_min == min, "mismatch for min");
      rc = fc_getVariableMinMax(global_vars[i][j], NULL, &temp_min_id, 
                                NULL, NULL);
      fail_unless(rc == FC_SUCCESS, "failed to get min/max");
      fail_unless(temp_min_id == min_id, "mismatch for min id");
      rc = fc_getVariableMinMax(global_vars[i][j], NULL, NULL, 
                                &temp_max, NULL);
      fail_unless(rc == FC_SUCCESS, "failed to get min/max");
      fail_unless(temp_max == max, "mismatch for max");
      rc = fc_getVariableMinMax(global_vars[i][j], NULL, NULL, NULL, 
                                &temp_max_id);
      fail_unless(rc == FC_SUCCESS, "failed to get min/max");
      fail_unless(temp_max_id == max_id, "mismatch for max id");

      // fc_getVariableMeanSdev - stdev is optional argument
      temp_mean = -1;
      rc = fc_getVariableMeanSdev(global_vars[i][j], &temp_mean, NULL);
      fail_unless(rc == FC_SUCCESS, "failed to get mean/stdev");
      fail_unless(temp_mean == mean, "mismatch for mean");

      // --- test errors

      // bad args to fc_getVariableMinMax()
      rc = fc_getVariableMinMax(badVar, &temp_min, &temp_min_id, 
                                &temp_max, &temp_max_id);
      fail_unless(rc != FC_SUCCESS, "should fail with bad var");
      fail_unless(temp_min == -1 && temp_min_id == -1 && temp_max == -1 &&
                  temp_max_id == -1, "fail should return nulls");
      rc = fc_getVariableMinMax(global_vars[i][j], NULL, NULL, NULL, NULL);
      fail_unless(rc != FC_SUCCESS, "should fail if all return args are null");

      // bad args to fc_getVariableSum()
      rc = fc_getVariableSum(badVar, &temp_sum);
      fail_unless(rc != FC_SUCCESS, "should fail with bad var");
      fail_unless(temp_sum == -1, "fail should return null");
      rc = fc_getVariableSum(global_vars[i][j], NULL);
      fail_unless(rc != FC_SUCCESS, "should if no");

      // bad args to fc_getVariableMeanSdev()
      rc = fc_getVariableMeanSdev(badVar, &temp_mean, &temp_stdev);
      fail_unless(rc != FC_SUCCESS, "should fail if bad var");
      fail_unless(temp_mean == -1 && temp_stdev == -1, 
                  "fail should return null");
      rc = fc_getVariableMeanSdev(global_vars[i][j], NULL, &temp_stdev);
      fail_unless(rc != FC_SUCCESS, "should fail if no mean");
      fail_unless(temp_stdev == -1, "fail should return null");
    }
  }

  // --- more test of bad input

  // stats can't do char data
  rc = fc_getVariableMinMax(global_char_var, &temp_min, &temp_min_id, 
                            &temp_max, &temp_max_id);
  fail_unless(rc != FC_SUCCESS, "should fail to get min/max for char var");
  fail_unless(temp_min == -1 && temp_min_id == -1 && temp_max == -1 &&
              temp_max_id == -1, "fail should return nulls");
  rc = fc_getVariableSum(global_char_var, &temp_sum);
  fail_unless(rc != FC_SUCCESS, "should fail to get sum for char var");
  fail_unless(temp_sum == -1, "fail should return null");
  rc = fc_getVariableMeanSdev(global_char_var, &temp_mean, &temp_stdev);
  fail_unless(rc != FC_SUCCESS, "should fail to get mean/stdev for char var");
  fail_unless(temp_mean == -1 && temp_stdev == -1, "fail should return null");
  
  // stats can't do multi component data
  rc = fc_getVariableMinMax(global_vector_var, &temp_min, &temp_min_id, 
                            &temp_max, &temp_max_id);
  fail_unless(rc != FC_SUCCESS, "should fail to get min/max for vec var");
  fail_unless(temp_min == -1, "mismatch for min");
  fail_unless(temp_min_id == -1, "mismatch for min id");
  fail_unless(temp_max == -1, "mismatch for max");
  fail_unless(temp_max_id == -1, "mismatch for max id");
  rc = fc_getVariableSum(global_vector_var, &temp_sum);
  fail_unless(rc != FC_SUCCESS, "should fail to get sum for vec var");
  fail_unless(temp_sum == -1, "mismatch for sum");
  rc = fc_getVariableMeanSdev(global_vector_var, &temp_mean, &temp_stdev);
  fail_unless(rc != FC_SUCCESS, "should fail to get mean/stdev for vec var");
  fail_unless(temp_mean == -1, "mismatch for mean");
  fail_unless(temp_stdev == -1, "mismatch for stdev");
}
END_TEST
  
START_TEST(varSubset_stats)
{
  FC_ReturnCode rc;
  FC_Subset badSubset = { 999, 999 };
  FC_Variable badVar = { 999, 999 };
  int i, j, k;
  int numMember;
  double min, max, temp_min, temp_max;
  int min_id, max_id, temp_min_id, temp_max_id;
  double mean, stdev, temp_mean, temp_stdev;
  double sum, temp_sum;

  for (i = 0; i < global_numAssoc; i++) {
    if (global_assocs[i] == FC_AT_WHOLE_MESH) {
      numMember = 1;
      min = max = global_min;
      min_id = max_id = global_min_id;
      sum = global_min;
    }
    else {
      numMember = global_numMember;
      min = global_min;
      max = global_max;
      min_id = global_min_id;
      max_id = global_max_id;
      sum = global_min + global_max;
    }
    mean = sum / numMember;
    if (numMember == 1)
      stdev = 0.;
    else 
      stdev = sqrt( ( (global_min - mean)*(global_min - mean) + 
                      (global_max - mean)*(global_max - mean) + 
                      (numMember-2)*mean*mean ) / (numMember - 1) );

    for (j = 0; j < global_numDataType; j++) {
      //printf("i = %d, j = %d\n", i, j);

      // test MinMax
      rc = fc_getVariableSubsetMinMax(global_vars[i][j], global_subsets[i],
                                      &temp_min, &temp_min_id, 
                                      &temp_max, &temp_max_id);
      fail_unless(rc == FC_SUCCESS, "failed to get min/max");
      fail_unless(temp_min == min, "mismatch for min");
      fail_unless(temp_min_id == min_id, "mismatch for min id");
      fail_unless(temp_max == max, "mismatch for max");
      fail_unless(temp_max_id == max_id, "mismatch for max id");

      // test Sum
      rc = fc_getVariableSubsetSum(global_vars[i][j], global_subsets[i],
                                   &temp_sum);
      fail_unless(rc == FC_SUCCESS, "failed to get sum");
      fail_unless(temp_sum == sum, "mismatch for sum");

      // test MeanSdev
      rc = fc_getVariableSubsetMeanSdev(global_vars[i][j], global_subsets[i],
                                        &temp_mean, &temp_stdev);
      fail_unless(rc == FC_SUCCESS, "failed to get mean/stdev");
      fail_unless(temp_mean == mean, "mismatch for mean");
      fail_unless(FC_DBL_EQUIV(temp_stdev, stdev), "mismatch for stdev");

      // --- test special cases

      // fc_getVariableSubsetMinMax can take any of the return args as optional
      temp_min = temp_min_id = temp_max = temp_max_id = -1;
      rc = fc_getVariableSubsetMinMax(global_vars[i][j], global_subsets[i],
                                      &temp_min, NULL, NULL, NULL);
      fail_unless(rc == FC_SUCCESS, "failed to get min/max");
      fail_unless(temp_min == min, "mismatch for min");
      rc = fc_getVariableSubsetMinMax(global_vars[i][j], global_subsets[i],
                                      NULL, &temp_min_id, NULL, NULL);
      fail_unless(rc == FC_SUCCESS, "failed to get min/max");
      fail_unless(temp_min_id == min_id, "mismatch for min id");
      rc = fc_getVariableSubsetMinMax(global_vars[i][j], global_subsets[i],
                                      NULL, NULL, &temp_max, NULL);
      fail_unless(rc == FC_SUCCESS, "failed to get min/max");
      fail_unless(temp_max == max, "mismatch for max");
      rc = fc_getVariableSubsetMinMax(global_vars[i][j], global_subsets[i],
                                      NULL, NULL, NULL, &temp_max_id);
      fail_unless(rc == FC_SUCCESS, "failed to get min/max");
      fail_unless(temp_max_id == max_id, "mismatch for max id");

      // fc_getVariableSubsetMeanSdev - stdev is optional argument
      temp_mean = -1;
      rc = fc_getVariableSubsetMeanSdev(global_vars[i][j], global_subsets[i],
                                        &temp_mean, NULL);
      fail_unless(rc == FC_SUCCESS, "failed to get mean/stdev");
      fail_unless(temp_mean == mean, "mismatch for mean");

      // --- test errors

      // bad args to fc_getVariableSubsetMinMax()
      rc = fc_getVariableSubsetMinMax(badVar, global_subsets[i], &temp_min, 
                                      &temp_min_id, &temp_max, &temp_max_id);
      fail_unless(rc != FC_SUCCESS, "should fail with bad var");
      fail_unless(temp_min == -1 && temp_min_id == -1 && temp_max == -1 &&
                  temp_max_id == -1, "fail should return nulls");
      rc = fc_getVariableSubsetMinMax(global_vars[i][j], badSubset, &temp_min, 
                                      &temp_min_id, &temp_max, &temp_max_id);
      fail_unless(rc != FC_SUCCESS, "should fail with bad subset");
      fail_unless(temp_min == -1 && temp_min_id == -1 && temp_max == -1 &&
                  temp_max_id == -1, "fail should return nulls");
      rc = fc_getVariableSubsetMinMax(global_vars[i][j], global_subsets[i], 
                                      NULL, NULL, NULL, NULL);
      fail_unless(rc != FC_SUCCESS, "should fail if all return args are null");

      // bad args to fc_getVariableSubsetSum()
      rc = fc_getVariableSubsetSum(badVar, global_subsets[i], &temp_sum);
      fail_unless(rc != FC_SUCCESS, "should fail with bad var");
      fail_unless(temp_sum == -1, "fail should return null"); 
      rc = fc_getVariableSubsetSum(global_vars[i][j], badSubset, &temp_sum);
      fail_unless(rc != FC_SUCCESS, "should fail with bad subset");
      fail_unless(temp_sum == -1, "fail should return null");
      rc = fc_getVariableSubsetSum(global_vars[i][j], global_subsets[i], NULL);
      fail_unless(rc != FC_SUCCESS, "should if no");

      // bad args to fc_getVariableSubsetMeanSdev()
      rc = fc_getVariableSubsetMeanSdev(badVar, global_subsets[i], &temp_mean,
                                        &temp_stdev);
      fail_unless(rc != FC_SUCCESS, "should fail if bad var");
      fail_unless(temp_mean == -1 && temp_stdev == -1, 
                  "fail should return null");
      rc = fc_getVariableSubsetMeanSdev(global_vars[i][j], badSubset, 
                                        &temp_mean, &temp_stdev);
      fail_unless(rc != FC_SUCCESS, "should fail if bad subset");
      fail_unless(temp_mean == -1 && temp_stdev == -1, 
                  "fail should return null");
      rc = fc_getVariableSubsetMeanSdev(global_vars[i][j], global_subsets[i],
                                        NULL, &temp_stdev);
      fail_unless(rc != FC_SUCCESS, "should fail if no mean");
      fail_unless(temp_stdev == -1, "fail should return null");

      // assoc type has to match
      for (k = 0; k < global_numAssoc; k++) {
        if (k == i)
          continue;
        rc = fc_getVariableSubsetMinMax(global_vars[k][j], global_subsets[i], 
                             &temp_min, &temp_min_id, &temp_max, &temp_max_id);
        fail_unless(rc != FC_SUCCESS, "should fail with mixed assoc");
        fail_unless(temp_min == -1 && temp_min_id == -1 && temp_max == -1 &&
                    temp_max_id == -1, "fail should return nulls");
        rc = fc_getVariableSubsetSum(global_vars[k][j], global_subsets[i], 
                                     &temp_sum);
        fail_unless(rc != FC_SUCCESS, "should fail with mixed assoc");
        fail_unless(temp_sum == -1, "fail should return null");
        rc = fc_getVariableSubsetMeanSdev(global_vars[k][j], global_subsets[i], 
                                         &temp_mean, &temp_stdev);
        fail_unless(rc != FC_SUCCESS, "should fail if mixed assoc");
        fail_unless(temp_mean == -1 && temp_stdev == -1, 
                    "fail should return null");
      }
    }
  }

  // --- more test of bad input

  // stats can't do char data
  rc = fc_getVariableSubsetMinMax(global_char_var, global_subsets[0], &temp_min,
                                  &temp_min_id, &temp_max, &temp_max_id);
  fail_unless(rc != FC_SUCCESS, "should fail to get min/max for char var");
  fail_unless(temp_min == -1 && temp_min_id == -1 && temp_max == -1 &&
              temp_max_id == -1, "fail should return nulls");
  rc = fc_getVariableSubsetSum(global_char_var, global_subsets[0], &temp_sum);
  fail_unless(rc != FC_SUCCESS, "should fail to get sum for char var");
  fail_unless(temp_sum == -1, "fail should return null");
  rc = fc_getVariableSubsetMeanSdev(global_char_var, global_subsets[0], 
                                    &temp_mean, &temp_stdev);
  fail_unless(rc != FC_SUCCESS, "should fail to get mean/stdev for char var");
  fail_unless(temp_mean == -1 && temp_stdev == -1, "fail should return null");
  
  // stats can't do multi component data
  rc = fc_getVariableSubsetMinMax(global_vector_var, global_subsets[0], 
                          &temp_min, &temp_min_id, &temp_max, &temp_max_id);
  fail_unless(rc != FC_SUCCESS, "should fail to get min/max for vec var");
  fail_unless(temp_min == -1, "mismatch for min");
  fail_unless(temp_min_id == -1, "mismatch for min id");
  fail_unless(temp_max == -1, "mismatch for max");
  fail_unless(temp_max_id == -1, "mismatch for max id");
  rc = fc_getVariableSubsetSum(global_vector_var, global_subsets[0], &temp_sum);
  fail_unless(rc != FC_SUCCESS, "should fail to get sum for vec var");
  fail_unless(temp_sum == -1, "mismatch for sum");
  rc = fc_getVariableSubsetMeanSdev(global_vector_var, global_subsets[0], 
                                    &temp_mean, &temp_stdev);
  fail_unless(rc != FC_SUCCESS, "should fail to get mean/stdev for vec var");
  fail_unless(temp_mean == -1, "mismatch for mean");
  fail_unless(temp_stdev == -1, "mismatch for stdev");
}
END_TEST

START_TEST(seqVar_stats)
{
  FC_ReturnCode rc;
  FC_Variable badVar = { 999, 999 };
  int i, j, k, m;
  int numEntity;
  double min, max, temp_min, temp_max;
  int min_id, temp_min_step_id, temp_min_entity_id;
  int max_id, temp_max_step_id, temp_max_entity_id;
  double mean, stdev, temp_mean, temp_stdev;
  double sum/*, *temp_sum*/;
  FC_Sequence sequence;
  FC_Variable badSeqVar[5];
  int numStep = 5, seq_coords[5] = { 0, 1, 2, 3, 4 };
  int min_step_id = 1, max_step_id = 3;

  // make a sequence
  rc = fc_createSequence(global_dataset, "time", &sequence);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create sequence");
  rc = fc_setSequenceCoords(sequence, numStep, FC_DT_INT, seq_coords);
  fail_unless(rc == FC_SUCCESS, "abort: failed to set seq coords"); 

  // The seqvar we end up making will be
  // step ID 0 - will be our original
  // step ID 1 - subtract 1, will have smallest value
  // step ID 2 - will be our original
  // step ID 3 - add 1, will have largest value
  // step ID 4 - will be our original

  for (i = 0; i < global_numAssoc; i++) {
    fc_getMeshNumEntity(global_mesh, global_assocs[i], &numEntity);
    if (global_assocs[i] == FC_AT_WHOLE_MESH) {
      min = global_min - 1;
      min_id = global_min_id;
      max = global_min + 1;
      max_id = global_min_id;
      sum = global_min*numStep;
    }
    else {
      min = global_min - 1;
      max = global_max + 1;
      min_id = global_min_id;
      max_id = global_max_id;
      sum = (global_min + global_max)*numStep;
    }
    mean = sum / (numEntity*numStep);
    if (numEntity == 1)
      stdev = sqrt( ( (global_min - mean)*(global_min - mean)*(numStep - 2) +
		      (global_min - 1 - mean)*(global_min - 1 - mean) +
		      (global_min + 1 - mean)*(global_min + 1 - mean) ) / 
		    (numStep - 1));
    else 
      // purposefully different way of computing compared to implementation
      // stdev = sqrt ( (sum (xi - mean)^2)/(num - 1) )
      stdev = sqrt( ( (global_min - mean)*(global_min - mean)*(numStep - 2) + 
                      (global_max - mean)*(global_max - mean)*(numStep - 2) + 
		      (global_min - 1 - mean)*(global_min - 1 - mean) +
		      (global_min + 1 - mean)*(global_min + 1 - mean) +
		      (global_max - 1 - mean)*(global_max - 1 - mean) +
		      (global_max + 1 - mean)*(global_max + 1 - mean) +
		      (numEntity - 2)*mean*mean*(numStep - 2) +
		      (numEntity - 2)*(1 - mean)*(1 - mean) +
		      (numEntity - 2)*(1 + mean)*(1 + mean) ) /
		    (numEntity*numStep - 1) );
  
    for (j = 0; j < global_numDataType; j++) {
      FC_Variable *seqVar, vars[5];
      void* data, *data_copy;
      rc = fc_getVariableDataPtr(global_vars[i][j], &data);
      for (k = 0; k < numStep; k++) {
	data_copy = malloc(numEntity*fc_sizeofDataType(global_dataTypes[j]));
	fail_unless(data_copy != NULL, "abort: malloc error");
	memcpy(data_copy, data, numEntity*fc_sizeofDataType(global_dataTypes[j]));
	if (k == min_step_id || k == max_step_id) {
	  int adj = -1;
	  if (k == max_step_id)
	    adj = 1;
	  switch(global_dataTypes[j]) {
	  case FC_DT_INT: 
	    for (m = 0; m < numEntity; m++)
	      ((int*)data_copy)[m] += adj;
	    break;
	  case FC_DT_FLOAT:
	    for (m = 0; m < numEntity; m++)
	      ((float*)data_copy)[m] += adj;
	    break;
	  case FC_DT_DOUBLE:
	    for (m = 0; m < numEntity; m++)
	      ((double*)data_copy)[m] += adj;
	    break;
	  default:
	    ; // nothing
	  }	
	}
	rc = fc_createVariable(global_mesh, "temp", &vars[k]);
	fail_unless(rc == FC_SUCCESS, "abort: failed to create variable");
	rc = fc_setVariableDataPtr(vars[k], numEntity, 1, global_assocs[i],
			       FC_MT_SCALAR, global_dataTypes[j], data_copy);
	fail_unless(rc == FC_SUCCESS, "abort: faile to set variable data");
      }
      rc = fc_convertVariablesToSeqVariable(numStep, vars, sequence, "temp",
					    &seqVar);
      fail_unless(rc == FC_SUCCESS, 
		  "abort: failed to convert vars to seq var");

      //printf("i = %d, j = %d\n", i, j);

      // test MinMax
      rc = fc_getSeqVariableMinMax(numStep, seqVar, &temp_min,
				   &temp_min_step_id, &temp_min_entity_id, 
				   &temp_max, &temp_max_step_id, 
				   &temp_max_entity_id);
      fail_unless(rc == FC_SUCCESS, "failed to get min/max");
      fail_unless(temp_min == min, "mismatch for min");
      fail_unless(temp_min_step_id == min_step_id, "mismatch of min step id");
      fail_unless(temp_min_entity_id == min_id, 
		  "mismatch for min entity id");
      fail_unless(temp_max == max, "mismatch for max");
      fail_unless(temp_max_step_id == max_step_id, "mismatch of max step id");
      fail_unless(temp_max_entity_id == max_id, 
		  "mismatch for max entity id");
      /*
      // test Sum
      rc = fc_getVariableSum(global_vars[i][j], &temp_sum);
      fail_unless(rc == FC_SUCCESS, "failed to get sum");
      fail_unless(temp_sum == sum, "mismatch for sum");
      */

      // test MeanSdev
      rc = fc_getSeqVariableMeanSdev(numStep, seqVar, &temp_mean, &temp_stdev);
      fail_unless(rc == FC_SUCCESS, "failed to get mean/stdev");
      fail_unless(temp_mean == mean, "mismatch for mean");
      fail_unless(FC_DBL_EQUIV(temp_stdev, stdev), "mismatch for stdev");
      
      // --- test special cases
      
      // fc_getSeqVariableMinMax can take any of the return args as optional
      temp_min = temp_max = -1;
      temp_min_step_id = temp_max_step_id = -1;
      temp_min_entity_id = temp_max_entity_id = -1;
      rc = fc_getSeqVariableMinMax(numStep, seqVar, &temp_min, NULL, NULL, 
				   NULL, NULL, NULL);
      fail_unless(rc == FC_SUCCESS, "failed to get min/max");
      fail_unless(temp_min == min, "mismatch for min");
      rc = fc_getSeqVariableMinMax(numStep, seqVar, NULL, &temp_min_step_id,
				   NULL, NULL, NULL, NULL);
      fail_unless(rc == FC_SUCCESS, "failed to get min/max");
      fail_unless(temp_min_step_id == min_step_id, "mismatch for min id");
      rc = fc_getSeqVariableMinMax(numStep, seqVar, NULL, NULL,
				   &temp_min_entity_id, NULL, NULL, NULL);
      fail_unless(rc == FC_SUCCESS, "failed to get min/max");
      fail_unless(temp_min_entity_id == min_id, "mismatch for min id");
      rc = fc_getSeqVariableMinMax(numStep, seqVar, NULL, NULL, NULL, 
				   &temp_max, NULL, NULL);
      fail_unless(rc == FC_SUCCESS, "failed to get min/max");
      fail_unless(temp_max == max, "mismatch for max");
      rc = fc_getSeqVariableMinMax(numStep, seqVar, NULL, NULL, NULL,
				   NULL, &temp_max_step_id, NULL);
      fail_unless(rc == FC_SUCCESS, "failed to get min/max");
      fail_unless(temp_max_step_id == max_step_id, "mismatch for max id");
      rc = fc_getSeqVariableMinMax(numStep, seqVar, NULL, NULL, NULL, 
				   NULL, NULL, &temp_max_entity_id);
      fail_unless(rc == FC_SUCCESS, "failed to get min/max");
      fail_unless(temp_max_entity_id == max_id, "mismatch for max id");

      // fc_getVariableMeanSdev - stdev is optional argument
      temp_mean = -1;
      rc = fc_getSeqVariableMeanSdev(numStep, seqVar, &temp_mean, NULL);
      fail_unless(rc == FC_SUCCESS, "failed to get mean/stdev");
      fail_unless(temp_mean == mean, "mismatch for mean");
      
      // --- test errors

      // setup bad seqVar
      for (k = 0; k < numStep; k++)
	badSeqVar[k] = seqVar[k];
      badSeqVar[numStep/2] = badVar;

      // bad args to fc_getVariableMinMax()
      rc = fc_getSeqVariableMinMax(numStep, badSeqVar, &temp_min, 
				   &temp_min_step_id, &temp_min_entity_id, 
				   &temp_max, 
				   &temp_max_step_id, &temp_max_entity_id);
      fail_unless(rc != FC_SUCCESS, "should fail with bad var");
      fail_unless(temp_min == -1 && temp_min_step_id == -1 &&
		  temp_min_entity_id == -1 && temp_max == -1 &&
                  temp_max_step_id == -1 && temp_max_entity_id, 
		  "fail should return nulls");
      rc = fc_getSeqVariableMinMax(numStep, seqVar, NULL, NULL, NULL, NULL,
				  NULL, NULL);
      fail_unless(rc != FC_SUCCESS, "should fail if all return args are null");
      /*
      // bad args to fc_getVariableSum()
      rc = fc_getVariableSum(badVar, &temp_sum);
      fail_unless(rc != FC_SUCCESS, "should fail with bad var");
      fail_unless(temp_sum == -1, "fail should return null");
      rc = fc_getVariableSum(global_vars[i][j], NULL);
      fail_unless(rc != FC_SUCCESS, "should if no");
      */
      // bad args to fc_getVariableMeanSdev()
      rc = fc_getSeqVariableMeanSdev(numStep, badSeqVar, &temp_mean, &temp_stdev);
      fail_unless(rc != FC_SUCCESS, "should fail if bad var");
      fail_unless(temp_mean == -1 && temp_stdev == -1, 
                  "fail should return null");
      rc = fc_getSeqVariableMeanSdev(numStep, seqVar, NULL, &temp_stdev);
      fail_unless(rc != FC_SUCCESS, "should fail if no mean");
      fail_unless(temp_stdev == -1, "fail should return null");

      // cleanup
      fc_deleteSeqVariable(numStep, seqVar);
      free(seqVar);
    }
  }

  // --- more test of bad input
  /*
  // stats can't do char data
  rc = fc_getVariableMinMax(global_char_var, &temp_min, &temp_min_id, 
                            &temp_max, &temp_max_id);
  fail_unless(rc != FC_SUCCESS, "should fail to get min/max for char var");
  fail_unless(temp_min == -1 && temp_min_id == -1 && temp_max == -1 &&
              temp_max_id == -1, "fail should return nulls");
  rc = fc_getVariableSum(global_char_var, &temp_sum);
  fail_unless(rc != FC_SUCCESS, "should fail to get sum for char var");
  fail_unless(temp_sum == -1, "fail should return null");
  rc = fc_getVariableMeanSdev(global_char_var, &temp_mean, &temp_stdev);
  fail_unless(rc != FC_SUCCESS, "should fail to get mean/stdev for char var");
  fail_unless(temp_mean == -1 && temp_stdev == -1, "fail should return null");
  
  // stats can't do multi component data
  rc = fc_getVariableMinMax(global_vector_var, &temp_min, &temp_min_id, 
                            &temp_max, &temp_max_id);
  fail_unless(rc != FC_SUCCESS, "should fail to get min/max for vec var");
  fail_unless(temp_min == -1, "mismatch for min");
  fail_unless(temp_min_id == -1, "mismatch for min id");
  fail_unless(temp_max == -1, "mismatch for max");
  fail_unless(temp_max_id == -1, "mismatch for max id");
  rc = fc_getVariableSum(global_vector_var, &temp_sum);
  fail_unless(rc != FC_SUCCESS, "should fail to get sum for vec var");
  fail_unless(temp_sum == -1, "mismatch for sum");
  rc = fc_getVariableMeanSdev(global_vector_var, &temp_mean, &temp_stdev);
  fail_unless(rc != FC_SUCCESS, "should fail to get mean/stdev for vec var");
  fail_unless(temp_mean == -1, "mismatch for mean");
  fail_unless(temp_stdev == -1, "mismatch for stdev");
*/
}
END_TEST

START_TEST(seqVar_series_stats)
{
  FC_ReturnCode rc;
  FC_Variable badVar = { 999, 999 };
  int i,j,k,l;
  int numEntity;
  double min, max, *temp_min, *temp_max;
  int min_id, max_id, *temp_minindex, *temp_maxindex;
  double mean, sdev, *temp_mean, *temp_sdev;
  double sum, *temp_sum;
  FC_MathType newmathtype;
  FC_DataType newdatatype;
  FC_AssociationType newassoc;
  int newnumdatapoint, newnumcomponent;

  FC_Sequence seq;
  FC_Variable *seqVar, *singleseqVar;
  FC_Variable componentVar;
  FC_Variable sumVar, *returnVars;
  FC_Variable meanVar, sdevVar;
  FC_Variable minVar,maxVar,minindexVar,maxindexVar;
  int numReturnVars;

  int numStep, numDataPoint;

  for (i = 0; i < global_numAssoc; i++) {
    fc_getMeshNumEntity(global_mesh, global_assocs[i], &numEntity);
    if (global_assocs[i] == FC_AT_WHOLE_MESH) {
      min = max = global_min;
      min_id = max_id = global_min_id;
      sum = global_min;
    }
    else {
      min = global_min;
      max = global_max;
      min_id = global_min_id;
      max_id = global_max_id;
      sum = global_min + global_max;
    }
    mean = sum / numEntity;
    if (numEntity == 1)
      sdev = 0.;
    else 
      sdev = sqrt( ( (global_min - mean)*(global_min - mean) + 
                      (global_max - mean)*(global_max - mean) + 
                      (numEntity-2)*mean*mean ) / (numEntity - 1) );
  
    for (j = 0; j < global_numDataType; j++) {
      //need multiple component var

      int numComponent = 2;
      FC_Variable* tempvars = (FC_Variable*)malloc(numComponent*
						   sizeof(FC_Variable));
      for (k=0; k < numComponent; k++){
	tempvars[k] = global_vars[i][j];      
      }

      //make multiple component var
      rc = fc_mergeComponentVariables(numComponent,tempvars,
				      "junk_component_var",
				      FC_MT_VECTOR,&componentVar);
      fail_unless(rc == FC_SUCCESS, "couldnt merge component variables");

      //change it into a seq var
      rc = fc_getVariableNumDataPoint(componentVar,&numDataPoint);
      fail_unless(rc == FC_SUCCESS, "couldnt get numdatapoint");
      rc = fc_createRegularSequence(global_dataset,numDataPoint,0,1,
				    "junk_seq",&seq);
      fail_unless(rc == FC_SUCCESS, "couldnt create regular sequnce");
      rc = _fc_createSeqVariableFromVariable(numDataPoint,seq,componentVar,
					      "junk_seq_var",
					      &singleseqVar);
      fail_unless(rc == FC_SUCCESS, "couldnt make seq var from regular var");

      //clean up
      free(tempvars);
      rc = fc_deleteVariable(componentVar);
      fail_unless(rc == FC_SUCCESS, "couldnt delete component variable");

      //change it into a multidata point seq var
      rc = ag_multiDataPointSeqVariable(numDataPoint, singleseqVar,
					     seq,numEntity,
					     global_assocs[i],
					     "seqv", &numStep,
					     &seqVar);
      fail_unless(rc == FC_SUCCESS, 
		  "couldnt make multidata point component seq variable");

      //clean up
      rc = fc_deleteSeqVariable(numStep,singleseqVar);
      fail_unless(rc == FC_SUCCESS, "couldnt delete single seq variable");
      free(singleseqVar);

      // test MinMax
      //      fc_printVariable(global_vars[i][j],"",0);
      rc = fc_getSeqVariableSeriesMinMax(numStep,seqVar,
					    &minVar,
					    &maxVar,
					    &minindexVar,
					    &maxindexVar);
      fail_unless(rc == FC_SUCCESS, "failed to get min/max");

      rc = fc_getVariableByName(global_mesh,"seqv_SeriesMin",
				&numReturnVars,&returnVars);
      fail_unless(rc == FC_SUCCESS, "failed to get min variable by name");
      fail_unless(numReturnVars == 1, "wont number of matching vars");
      fail_unless(FC_HANDLE_EQUIV(returnVars[0],minVar),
		  "mismatch of variable handles");
      free(returnVars);

      rc = fc_getVariableInfo(minVar,&newnumdatapoint, &newnumcomponent,
			      &newassoc,&newmathtype,&newdatatype);
      fail_unless(rc == FC_SUCCESS, "failed to get min variable info");
      fail_unless(newnumdatapoint == numDataPoint,"mismatch of data point");
      fail_unless(newnumcomponent == numComponent,"mismatch of numcomponent");
      fail_unless(newassoc == global_assocs[i],"mismatch of association");
      fail_unless(newmathtype == FC_MT_VECTOR,"mismatch of mathtype");
      fail_unless(newdatatype == FC_DT_DOUBLE,"mismatch of datatype");
      rc = fc_getVariableDataPtr(minVar,(void**)&temp_min);
      fail_unless(rc == FC_SUCCESS, "failed to get min variable data");

      rc = fc_getVariableByName(global_mesh,"seqv_SeriesMax",&numReturnVars,
				&returnVars);
      fail_unless(rc == FC_SUCCESS, "failed to get max variable by name");
      fail_unless(numReturnVars == 1, "wont number of matching vars");
      fail_unless(FC_HANDLE_EQUIV(returnVars[0],maxVar),
		  "mismatch of variable handles");
      free(returnVars);

      rc = fc_getVariableInfo(maxVar,&newnumdatapoint, &newnumcomponent,
			      &newassoc,&newmathtype,&newdatatype);
      fail_unless(rc == FC_SUCCESS, "failed to get max variable info");
      fail_unless(newnumdatapoint == numDataPoint,"mismatch of data point");
      fail_unless(newnumcomponent == numComponent,"mismatch of numcomponent");
      fail_unless(newassoc == global_assocs[i],"mismatch of association");
      fail_unless(newmathtype == FC_MT_VECTOR,"mismatch of mathtype");
      fail_unless(newdatatype == FC_DT_DOUBLE,"mismatch of datatype");
      rc = fc_getVariableDataPtr(maxVar,(void**)&temp_max);
      fail_unless(rc == FC_SUCCESS, "failed to get max variable data");


      rc = fc_getVariableByName(global_mesh,"seqv_SeriesMinIndex",
				&numReturnVars,&returnVars);
      fail_unless(rc == FC_SUCCESS, 
		  "failed to get min index variable by name");
      fail_unless(numReturnVars == 1, "wrong number of return vars");
      fail_unless(FC_HANDLE_EQUIV(returnVars[0],minindexVar),
		  "mismatch of variable handles");
      free(returnVars);

      rc = fc_getVariableInfo(minindexVar,&newnumdatapoint, &newnumcomponent,
			      &newassoc,&newmathtype,&newdatatype);
      fail_unless(rc == FC_SUCCESS, "failed to get min index variable info");
      fail_unless(newnumdatapoint == numDataPoint,"mismatch of data point");
      fail_unless(newnumcomponent == numComponent,"mismatch of numcomponent");
      fail_unless(newassoc == global_assocs[i],"mismatch of association");
      fail_unless(newmathtype == FC_MT_VECTOR,"mismatch of mathtype");
      fail_unless(newdatatype == FC_DT_INT,"mismatch of datatype");
      rc = fc_getVariableDataPtr(minindexVar,(void**)&temp_minindex);
      fail_unless(rc == FC_SUCCESS, "failed to get min indexvariable data");


      rc = fc_getVariableByName(global_mesh,"seqv_SeriesMaxIndex",
				&numReturnVars,&returnVars);
      fail_unless(rc == FC_SUCCESS, 
		  "failed to get max index variable by name");
      fail_unless(numReturnVars == 1, "wrong number of return vars");
      fail_unless(FC_HANDLE_EQUIV(returnVars[0],maxindexVar),
		  "mismatch of variable handles");
      free(returnVars);
      rc = fc_getVariableInfo(maxindexVar,&newnumdatapoint, &newnumcomponent,
			      &newassoc,&newmathtype,&newdatatype);
      fail_unless(rc == FC_SUCCESS, "failed to get max index variable info");
      fail_unless(newnumdatapoint == numDataPoint,"mismatch of data point");
      fail_unless(newnumcomponent == numComponent,"mismatch of numcomponent");
      fail_unless(newassoc == global_assocs[i],"mismatch of association");
      fail_unless(newmathtype == FC_MT_VECTOR,"mismatch of mathtype");
      fail_unless(newdatatype == FC_DT_INT,"mismatch of datatype");
      rc = fc_getVariableDataPtr(maxindexVar,(void**)&temp_maxindex);
      fail_unless(rc == FC_SUCCESS, "failed to get max indexvariable data");

      for (l = 0; l < numEntity; l++){
	for (k = 0; k < numComponent; k++){
	  //	printf("DP=%d COMP = %d minmax %f %d %f %d - %f %d %f %d\n", 
	  //	       l,k,min,min_id,max,max_id,temp_min[numComponent*l+k],
	  //	       temp_min_id[numComponent*l+k],
	  //	       temp_max[numComponent*l+k],
	  //	       temp_max_id[numComponent*l+k]);
	  if (l == 0){
	    fail_unless(temp_min[numComponent*l+k] == 1.0, "mismatch for min");
	    fail_unless(temp_minindex[numComponent*l+k] == 0,
			"mismatch for min id");
	    fail_unless(temp_max[numComponent*l+k] == 1.0, "mismatch for max");
	    fail_unless(temp_maxindex[numComponent*l+k] == 0,
			"mismatch for max id");
	  }else {
	    fail_unless(temp_min[numComponent*l+k] == min, "mismatch for min");
	    fail_unless(temp_minindex[numComponent*l+k] == min_id,
			"mismatch for min id");
	    fail_unless(temp_max[numComponent*l+k] == max, "mismatch for max");
	    fail_unless(temp_maxindex[numComponent*l+k] == max_id,
			"mismatch for max id");
	  }
	}
      }

      rc = fc_deleteVariable(minVar);
      fail_unless(rc == FC_SUCCESS, "failed to get delete min var");
      rc = fc_deleteVariable(maxVar);
      fail_unless(rc == FC_SUCCESS, "failed to get delete max var");
      rc = fc_deleteVariable(minindexVar);
      fail_unless(rc == FC_SUCCESS, "failed to get delete minindex var");
      rc = fc_deleteVariable(maxindexVar);
      fail_unless(rc == FC_SUCCESS, "failed to get delete maxindex var");


      // test Sum
      rc = fc_getSeqVariableSeriesSum(numStep,seqVar, 
					 &sumVar);
      fail_unless(rc == FC_SUCCESS, "failed to get sum");
      rc = fc_getVariableByName(global_mesh,"seqv_SeriesSum",&numReturnVars,
				&returnVars);
      fail_unless(rc == FC_SUCCESS, "failed to get sum variable by name");
      fail_unless(numReturnVars == 1, "failed to get matching var");
      fail_unless(FC_HANDLE_EQUIV(returnVars[0],sumVar),
		  "mismatch of variable handles");
      free(returnVars);
      rc = fc_getVariableInfo(sumVar,&newnumdatapoint, &newnumcomponent,
			      &newassoc,&newmathtype,&newdatatype);
      fail_unless(rc == FC_SUCCESS, "failed to get sum variable info");
      fail_unless(newnumdatapoint == numDataPoint,"mismatch of data point");
      fail_unless(newnumcomponent == numComponent,"mismatch of numcomponent");
      fail_unless(newassoc == global_assocs[i],"mismatch of association");
      fail_unless(newmathtype == FC_MT_VECTOR,"mismatch of mathtype");
      fail_unless(newdatatype == FC_DT_DOUBLE,"mismatch of datatype");

      rc = fc_getVariableDataPtr(sumVar,(void**)&temp_sum);
      fail_unless(rc == FC_SUCCESS, "failed to get sum variable data");
      for (l = 0; l < numEntity; l++){
	for (k = 0; k < numComponent; k++){
	  if (l == 0){
	    fail_unless(FC_DBL_EQUIV(temp_sum[numComponent*l+k],
				     numEntity), "mismatch for sum");
	  }else {
	    fail_unless(FC_DBL_EQUIV(temp_sum[numComponent*l+k],
				     sum), "mismatch for sum");
	  }
	}
      }

      rc = fc_deleteVariable(sumVar);
      fail_unless(rc == FC_SUCCESS,"failed to delete sumvar");

      // test MeanSdev
      rc = fc_getSeqVariableSeriesMeanSdev(numStep,seqVar,
					      &meanVar, &sdevVar);
      fail_unless(rc == FC_SUCCESS, "failed to get mean/stdev");
      rc = fc_getVariableByName(global_mesh,"seqv_SeriesMean",
				&numReturnVars,&returnVars);
      fail_unless(rc == FC_SUCCESS, "failed to get mean variable by name");
      fail_unless(numReturnVars == 1, "wrong number of matching vars");
      fail_unless(FC_HANDLE_EQUIV(returnVars[0],meanVar),
		  "mismatch of variable handles");
      free(returnVars);
      rc = fc_getVariableInfo(meanVar,&newnumdatapoint, &newnumcomponent,
			      &newassoc,&newmathtype,&newdatatype);
      fail_unless(rc == FC_SUCCESS, "failed to get mean variable info");
      fail_unless(newnumdatapoint == numDataPoint,"mismatch of data point");
      fail_unless(newnumcomponent == numComponent,"mismatch of numcomponent");
      fail_unless(newassoc == global_assocs[i],"mismatch of association");
      fail_unless(newmathtype == FC_MT_VECTOR,"mismatch of mathtype");
      fail_unless(newdatatype == FC_DT_DOUBLE,"mismatch of datatype");
      rc = fc_getVariableDataPtr(meanVar,(void**)&temp_mean);
      fail_unless(rc == FC_SUCCESS,"failed to get mean data");

      rc = fc_getVariableByName(global_mesh,"seqv_SeriesSdev",
				&numReturnVars,&returnVars);
      fail_unless(rc == FC_SUCCESS, "failed to get mean variable by name");
      fail_unless(numReturnVars == 1, "wrong number of matching vars");
      fail_unless(FC_HANDLE_EQUIV(returnVars[0],sdevVar),
		  "mismatch of variable handles");
      free(returnVars);
      rc = fc_getVariableInfo(sdevVar,&newnumdatapoint, &newnumcomponent,
			      &newassoc,&newmathtype,&newdatatype);
      fail_unless(rc == FC_SUCCESS, "failed to get mean variable info");
      fail_unless(newnumdatapoint == numDataPoint,"mismatch of data point");
      fail_unless(newnumcomponent == numComponent,"mismatch of numcomponent");
      fail_unless(newassoc == global_assocs[i],"mismatch of association");
      fail_unless(newmathtype == FC_MT_VECTOR,"mismatch of mathtype");
      fail_unless(newdatatype == FC_DT_DOUBLE,"mismatch of datatype");
      rc = fc_getVariableDataPtr(sdevVar,(void**)&temp_sdev);
      fail_unless(rc == FC_SUCCESS,"failed to get sdev data");

      for (l = 0; l < numEntity; l++){
	for (k = 0; k < numComponent; k++){
	  //	  printf("mean sdev %f %f - %f %f\n",
	  //		 mean,stdev,temp_mean[numComponent*l+k],
	  //		 temp_stdev[numComponent*l+k]);
	  if (l == 0){
	    fail_unless(FC_DBL_EQUIV(temp_mean[numComponent*l+k],1.0),
			"mismatch for mean");
	    fail_unless(FC_DBL_EQUIV(temp_sdev[numComponent*l+k], 0.0),
			"mismatch for stdev");
	  } else{
	    fail_unless(temp_mean[numComponent*l+k] == mean,
			"mismatch for mean");
	    fail_unless(FC_DBL_EQUIV(temp_sdev[numComponent*l+k], sdev),
		    "mismatch for stdev");
	  }
	}
      }

      rc = fc_deleteVariable(meanVar);
      fail_unless(rc == FC_SUCCESS, "failed to get delete mean var");
      rc = fc_deleteVariable(sdevVar);
      fail_unless(rc == FC_SUCCESS, "failed to get delete sdev var");


      // --- test special cases
      //unlike fc_getVariableMinMax, there are no optional arguments

      //      fc_getSeqVariableSeriesMeanSdev - stdev is optional argument
      rc = fc_getSeqVariableSeriesMeanSdev(numStep,seqVar,
					      &meanVar, NULL);
      fail_unless(rc == FC_SUCCESS, "failed to get mean/stdev");
      rc = fc_getVariableDataPtr(meanVar,(void**)&temp_mean);
      fail_unless(rc == FC_SUCCESS,"failed to get mean data");
      for (l = 0; l < numEntity; l++){
	for (k = 0; k < numComponent; k++){
	  if (l == 0){
	    fail_unless(FC_DBL_EQUIV(temp_mean[numComponent*l+k],1.0),
			"mismatch for mean");
	  }else {
	    fail_unless(FC_DBL_EQUIV(temp_mean[numComponent*l+k],mean),
			"mismatch for mean");
	  }
	}
      }
      rc = fc_deleteVariable(meanVar);
      fail_unless(rc == FC_SUCCESS, "failed to get delete mean var");

      rc = fc_deleteSeqVariable(numStep,seqVar);
      fail_unless(rc == FC_SUCCESS,
		  "failed to delete multi-data point seq var");
      free(seqVar);
      rc = fc_deleteSequence(seq);
      fail_unless(rc == FC_SUCCESS,"failed to delete sequence");
    }
  }


  // --- test of bad input

  //multiple component,single data pt var
  rc = fc_getVariableNumDataPoint(global_vars[0][0],&numDataPoint);
  fail_unless(rc == FC_SUCCESS, "couldnt get numdatapoint");

  numStep = numDataPoint;

  rc = fc_createRegularSequence(global_dataset,numStep,0,1,
				"junk_seq",&seq);
  fail_unless(rc == FC_SUCCESS, "couldnt create regular sequnce");

  rc = _fc_createSeqVariableFromVariable(numStep,seq,global_vars[0][0],
					"junk_seq_var",
					&seqVar);
  fail_unless(rc == FC_SUCCESS, "couldnt make seq var from regular var");


  // bad args to fc_getSeqVariableSeriesMinMax()
  rc = fc_getSeqVariableSeriesMinMax(numStep,&badVar, 
					&minVar,
					&maxVar,
					&minindexVar,
					&maxindexVar);
  fail_unless(rc != FC_SUCCESS, "should fail with bad var");
  fail_unless(FC_HANDLE_EQUIV(minVar,FC_NULL_VARIABLE),
	      "var should be null when fails");
  fail_unless(FC_HANDLE_EQUIV(maxVar,FC_NULL_VARIABLE),
	      "var should be null when fails");
  fail_unless(FC_HANDLE_EQUIV(minindexVar,FC_NULL_VARIABLE),
	      "var should be null when fails");
  fail_unless(FC_HANDLE_EQUIV(maxindexVar,FC_NULL_VARIABLE),
	      "var should be null when fails");



  rc = fc_getSeqVariableSeriesMinMax(numStep,seqVar, 
					NULL,
					&maxVar,
					&minindexVar,
					&maxindexVar);
  fail_unless(rc != FC_SUCCESS, "should fail with NULL var");
  fail_unless(FC_HANDLE_EQUIV(minVar,FC_NULL_VARIABLE),
	      "var should be null when fails");
  fail_unless(FC_HANDLE_EQUIV(maxVar,FC_NULL_VARIABLE),
	      "var should be null when fails");
  fail_unless(FC_HANDLE_EQUIV(minindexVar,FC_NULL_VARIABLE),
	      "var should be null when fails");
  fail_unless(FC_HANDLE_EQUIV(maxindexVar,FC_NULL_VARIABLE),
	      "var should be null when fails");



  rc = fc_getSeqVariableSeriesMinMax(numStep,seqVar, 
					&minVar,
					NULL,
					&minindexVar,
					&maxindexVar);
  fail_unless(rc != FC_SUCCESS, "should fail with NULL var");
  fail_unless(FC_HANDLE_EQUIV(minVar,FC_NULL_VARIABLE),
	      "var should be null when fails");
  fail_unless(FC_HANDLE_EQUIV(maxVar,FC_NULL_VARIABLE),
	      "var should be null when fails");
  fail_unless(FC_HANDLE_EQUIV(minindexVar,FC_NULL_VARIABLE),
	      "var should be null when fails");
  fail_unless(FC_HANDLE_EQUIV(maxindexVar,FC_NULL_VARIABLE),
	      "var should be null when fails");



  rc = fc_getSeqVariableSeriesMinMax(numStep,seqVar, 
					&minVar,
					&maxVar,
					NULL,	
					&maxindexVar);
  fail_unless(rc != FC_SUCCESS, "should fail with NULL var");
  fail_unless(FC_HANDLE_EQUIV(minVar,FC_NULL_VARIABLE),
	      "var should be null when fails");
  fail_unless(FC_HANDLE_EQUIV(maxVar,FC_NULL_VARIABLE),
	      "var should be null when fails");
  fail_unless(FC_HANDLE_EQUIV(minindexVar,FC_NULL_VARIABLE),
	      "var should be null when fails");
  fail_unless(FC_HANDLE_EQUIV(maxindexVar,FC_NULL_VARIABLE),
	      "var should be null when fails");


  rc = fc_getSeqVariableSeriesMinMax(numStep,seqVar, 
					&minVar,
					&maxVar,
					&minindexVar,
					NULL);
  fail_unless(rc != FC_SUCCESS, "should fail with NULL var");
  fail_unless(FC_HANDLE_EQUIV(minVar,FC_NULL_VARIABLE),
	      "var should be null when fails");
  fail_unless(FC_HANDLE_EQUIV(maxVar,FC_NULL_VARIABLE),
	      "var should be null when fails");
  fail_unless(FC_HANDLE_EQUIV(minindexVar,FC_NULL_VARIABLE),
	      "var should be null when fails");
  fail_unless(FC_HANDLE_EQUIV(maxindexVar,FC_NULL_VARIABLE),
	      "var should be null when fails");


  // bad args to fc_getVariableSeriesMeanSdev()
  rc = fc_getSeqVariableSeriesMeanSdev(numStep,&badVar,
					  &meanVar, 
					  &sdevVar);
  fail_unless(rc != FC_SUCCESS, "should fail if bad var");
  fail_unless(FC_HANDLE_EQUIV(sdevVar,FC_NULL_VARIABLE),
	      "var should be null when fails");
  fail_unless(FC_HANDLE_EQUIV(meanVar,FC_NULL_VARIABLE),
	      "var should be null when fails");


  rc = fc_getSeqVariableSeriesMeanSdev(numStep,seqVar,
					  NULL, &sdevVar);
  fail_unless(rc != FC_SUCCESS, "should fail if NULL mean var");
  fail_unless(FC_HANDLE_EQUIV(sdevVar,FC_NULL_VARIABLE),
	      "var should be null when fails");
  fail_unless(FC_HANDLE_EQUIV(meanVar,FC_NULL_VARIABLE),
	      "var should be null when fails");
  //checked sdev in special cases above



  // bad args to fc_getSeqVariableSeriesSum()
  rc = fc_getSeqVariableSeriesSum(numStep,&badVar, &sumVar);
  fail_unless(rc != FC_SUCCESS, "should fail with bad var");
  fail_unless(FC_HANDLE_EQUIV(sumVar,FC_NULL_VARIABLE),
	      "should return NULL var when fail");
  rc = fc_getSeqVariableSeriesSum(numStep,seqVar, NULL);
  fail_unless(rc != FC_SUCCESS, "should fail with NULL sumvar");



  rc = fc_deleteSeqVariable(numStep,seqVar);
  fail_unless(rc == FC_SUCCESS,"failed to delete seq var");
  free(seqVar);
  rc = fc_deleteSequence(seq);
  fail_unless(rc == FC_SUCCESS,"failed to delete sequence");

  //cant do char
  rc = fc_getVariableNumDataPoint(global_char_var,&numDataPoint);
  fail_unless(rc == FC_SUCCESS, "couldnt get numdatapoint");

  numStep = numDataPoint;

  rc = fc_createRegularSequence(global_dataset,numStep,0,1,
				"junk_seq",&seq);
  fail_unless(rc == FC_SUCCESS, "couldnt create regular sequnce");

  rc = _fc_createSeqVariableFromVariable(numStep,seq,global_char_var,
					"junk_seq_var",&seqVar);
  fail_unless(rc == FC_SUCCESS, "couldnt make seq var from regular var");


  rc = fc_getSeqVariableSeriesMinMax(numStep,seqVar, 
					&minVar,&maxVar,
					&minindexVar,
					&maxindexVar);
  fail_unless(rc != FC_SUCCESS, "should fail to get min/max for char var");
  fail_unless(FC_HANDLE_EQUIV(minVar,FC_NULL_VARIABLE),
	      "var should be null when fails");
  fail_unless(FC_HANDLE_EQUIV(maxVar,FC_NULL_VARIABLE),
	      "var should be null when fails");
  fail_unless(FC_HANDLE_EQUIV(minindexVar,FC_NULL_VARIABLE),
	      "var should be null when fails");
  fail_unless(FC_HANDLE_EQUIV(maxindexVar,FC_NULL_VARIABLE),
	      "var should be null when fails");



  rc = fc_getSeqVariableSeriesSum(numStep,seqVar,&sumVar);
  fail_unless(rc != FC_SUCCESS, "should fail to get sum for char var");
  fail_unless(FC_HANDLE_EQUIV(sumVar,FC_NULL_VARIABLE),
	      "should return NULL var when fail");

  rc = fc_getSeqVariableSeriesMeanSdev(numStep,seqVar, 
					  &meanVar,&sdevVar);
  fail_unless(rc != FC_SUCCESS, "should fail to get mean/stdev for char var");
  fail_unless(FC_HANDLE_EQUIV(sdevVar,FC_NULL_VARIABLE),
	      "var should be null when fails");
  fail_unless(FC_HANDLE_EQUIV(meanVar,FC_NULL_VARIABLE),
	      "var should be null when fails");


  rc = fc_deleteSeqVariable(numStep,seqVar);
  fail_unless(rc == FC_SUCCESS,"failed to delete seq var");
  free(seqVar);
  rc = fc_deleteSequence(seq);
  fail_unless(rc == FC_SUCCESS,"failed to delete sequence");

}
END_TEST


START_TEST(seq_stats)
{
  FC_ReturnCode rc;
  int i,j;
  double min, max, compval;
  int minindex, maxindex, mono;
  double ssmin, ssmax, ssmean,sssdev;
  int ssminindex, ssmaxindex;
  int compindex;
  FC_Sequence seq;
  int numStep = 4;


  char* comp[6] = {"<", "<=", "=", ">=", ">", "!="};
  double coords[4] = {0.0,1.0,1.5,2.0};
  float* fcoords;
  int* icoords;
  char* ccoords;
  
  //check validity

  //increasing seq
  rc= fc_createSequence(global_dataset,"seq",&seq);
  fail_unless(rc == FC_SUCCESS, "couldnt create sequence");
  rc = fc_setSequenceCoords(seq,numStep,FC_DT_DOUBLE,(void*)coords);
  fail_unless(rc == FC_SUCCESS, "couldnt set sequence coords");

  rc = fc_getSequenceMinMaxMono(seq,&min,&minindex,&max,
				&maxindex,&mono);
  fail_unless(rc== FC_SUCCESS, "failed to get sequence minmaxmono");
  fail_unless(FC_DBL_EQUIV(min,0.0), "invalid return value");
  fail_unless(FC_DBL_EQUIV(max,2.0), "invalid return value");
  fail_unless(minindex ==0, "invalid return value");
  fail_unless(maxindex ==3, "invalid return value");
  fail_unless(mono == 1, "invalid return value");

  rc = fc_getSequenceSpacingMinMaxMeanSdev(seq,&ssmin,&ssminindex,&ssmax,&ssmaxindex,
					&ssmean,&sssdev);
  fail_unless(rc== FC_SUCCESS, "failed to get sequence spacing minmax");
  fail_unless(FC_DBL_EQUIV(ssmin,0.5), "invalid return value");
  fail_unless(FC_DBL_EQUIV(ssmax,1.0), "invalid return value");
  fail_unless(ssminindex ==1, "invalid return value");
  fail_unless(ssmaxindex ==0, "invalid return value");
  fail_unless(FC_DBL_EQUIV(ssmean,2.0/3.0), "invalid return value");
  fail_unless(FC_DBL_EQUIV(sssdev,sqrt((1.5-4.0/3.0)/2.0)), "invalid return value");

  rc = fc_deleteSequence(seq);
  
  //decreasing
  for (i = 0; i < numStep; i++){
    coords[i] *= -1;
  }
  rc= fc_createSequence(global_dataset,"seq",&seq);
  fail_unless(rc == FC_SUCCESS, "couldnt create sequence");
  rc = fc_setSequenceCoords(seq,numStep,FC_DT_DOUBLE,(void*)coords);
  fail_unless(rc == FC_SUCCESS, "couldnt set sequence coords");
  rc = fc_getSequenceMinMaxMono(seq,&min,&minindex,&max,
				&maxindex,&mono);
  fail_unless(rc== FC_SUCCESS, "failed to get sequence minmaxmono");
  fail_unless(FC_DBL_EQUIV(max,0.0), "invalid return value");
  fail_unless(FC_DBL_EQUIV(min,-2.0), "invalid return value");
  fail_unless(maxindex ==0, "invalid return value");
  fail_unless(minindex ==3, "invalid return value");
  fail_unless(mono == -1, "invalid return value");

  rc = fc_getSequenceSpacingMinMaxMeanSdev(seq,&ssmin,&ssminindex,&ssmax,&ssmaxindex,
					&ssmean,&sssdev);
  fail_unless(rc== FC_SUCCESS, "failed to get sequence spacing minmax");
  fail_unless(FC_DBL_EQUIV(ssmin,-1.0), "invalid return value");
  fail_unless(FC_DBL_EQUIV(ssmax,-0.5), "invalid return value");
  fail_unless(ssminindex ==0, "invalid return value");
  fail_unless(ssmaxindex ==1, "invalid return value");
  fail_unless(FC_DBL_EQUIV(ssmean,-2.0/3.0), "invalid return value");
  fail_unless(FC_DBL_EQUIV(sssdev,sqrt((1.5-4.0/3.0)/2.0)), "invalid return value");

  rc = fc_deleteSequence(seq);

  //neither & repeating vals
  coords[0] = 1.0;
  coords[1] = 0.5;
  coords[2] = 0.5;
  coords[3] = 2.0;

  rc= fc_createSequence(global_dataset,"seq",&seq);
  fail_unless(rc == FC_SUCCESS, "couldnt create sequence");
  rc = fc_setSequenceCoords(seq,numStep,FC_DT_DOUBLE,(void*)coords);
  fail_unless(rc == FC_SUCCESS, "couldnt set sequence coords");
  rc = fc_getSequenceMinMaxMono(seq,&min,&minindex,&max,
				&maxindex,&mono);
  fail_unless(rc== FC_SUCCESS, "failed to get sequence minmaxmono");
  fail_unless(FC_DBL_EQUIV(max,2.0), "invalid return value");
  fail_unless(FC_DBL_EQUIV(min,0.5), "invalid return value");
  fail_unless(maxindex ==3, "invalid return value");
  fail_unless(minindex ==1, "invalid return value");
  fail_unless(mono == 0, "invalid return value");

  rc = fc_getSequenceSpacingMinMaxMeanSdev(seq,&ssmin,&ssminindex,&ssmax,&ssmaxindex,
					&ssmean,&sssdev);
  fail_unless(rc== FC_SUCCESS, "failed to get sequence spacing minmax");
  fail_unless(FC_DBL_EQUIV(ssmin,-0.5), "invalid return value");
  fail_unless(FC_DBL_EQUIV(ssmax,1.5), "invalid return value");
  fail_unless(ssminindex ==0, "invalid return value");
  fail_unless(ssmaxindex ==2, "invalid return value");
  fail_unless(FC_DBL_EQUIV(ssmean,1.0/3.0), "invalid return value");
  fail_unless(FC_DBL_EQUIV(sssdev,sqrt((2.5-1.0/3.0)/2.0)), "invalid return value");

  //comparison tests use known values
  for (i = 0; i < 6; i++){
    char* crit = comp[i];
    rc = fc_getClosestSequenceValue(seq,crit,min,&compval,&compindex);
    fail_unless(rc== FC_SUCCESS, "failed to get sequence closest value");
    switch (i) {
    case 0:
      fail_unless(compindex == -1, "invalid return value");
      fail_unless(FC_DBL_EQUIV(compval,-1), "invalid return value");
      break;
    case 1:
      fail_unless(compindex == minindex, "invalid return value");
      fail_unless(FC_DBL_EQUIV(compval,min), "invalid return value");
      break;
    case 2:
      fail_unless(compindex == minindex, "invalid return value");
      fail_unless(FC_DBL_EQUIV(compval,min), "invalid return value");
      break;
    case 3:
      fail_unless(compindex == minindex, "invalid return value");
      fail_unless(FC_DBL_EQUIV(compval,min), "invalid return value");
      break;
    case 4:
      fail_unless(compindex == 0, "invalid return value");
      fail_unless(FC_DBL_EQUIV(compval,1.0), "invalid return value");
      break;
    case 5:
      fail_unless(compindex == 0, "invalid return value");
      fail_unless(FC_DBL_EQUIV(compval,1.0), "invalid return value");
      break;
    default:
      break;
    }

    rc = fc_getClosestSequenceValue(seq,crit,0.9,&compval,&compindex);
    fail_unless(rc== FC_SUCCESS, "failed to get sequence closest value");
    switch (i) {
    case 0:
      fail_unless(compindex == 1, "invalid return value");
      fail_unless(FC_DBL_EQUIV(compval,0.5), "invalid return value");
      break;
    case 1:
      fail_unless(compindex == 1, "invalid return value");
      fail_unless(FC_DBL_EQUIV(compval,0.5), "invalid return value");
      break;
    case 2:
      fail_unless(compindex == 0, "invalid return value");
      fail_unless(FC_DBL_EQUIV(compval,1.0), "invalid return value");
      break;
    case 3:
      fail_unless(compindex == 0, "invalid return value");
      fail_unless(FC_DBL_EQUIV(compval,1.0), "invalid return value");
      break;
    case 4:
      fail_unless(compindex == 0, "invalid return value");
      fail_unless(FC_DBL_EQUIV(compval,1.0), "invalid return value");
      break;
    case 5:
      fail_unless(compindex == 0, "invalid return value");
      fail_unless(FC_DBL_EQUIV(compval,1.0), "invalid return value");
      break;
    default:
      break;
    }

    rc = fc_getClosestSequenceValue(seq,crit,2.2,&compval,&compindex);
    fail_unless(rc== FC_SUCCESS, "failed to get sequence closest value");
    switch (i) {
    case 0:
      fail_unless(compindex == 3, "invalid return value");
      fail_unless(FC_DBL_EQUIV(compval,2.0), "invalid return value");
      break;
    case 1:
      fail_unless(compindex == 3, "invalid return value");
      fail_unless(FC_DBL_EQUIV(compval,2.0), "invalid return value");
      break;
    case 2:
      fail_unless(compindex == 3, "invalid return value");
      fail_unless(FC_DBL_EQUIV(compval,2.0), "invalid return value");
      break;
    case 3:
      fail_unless(compindex == -1, "invalid return value");
      fail_unless(FC_DBL_EQUIV(compval,-1.0), "invalid return value");
      break;
    case 4:
      fail_unless(compindex == -1, "invalid return value");
      fail_unless(FC_DBL_EQUIV(compval,-1.0), "invalid return value");
      break;
    case 5:
      fail_unless(compindex == 3, "invalid return value");
      fail_unless(FC_DBL_EQUIV(compval,2.0), "invalid return value");
      break;
    default:
      break;
    }
  }


  //make sure NULL args work
  rc = fc_getSequenceMinMaxMono(seq,&min,NULL,NULL,NULL,NULL);
  fail_unless(rc== FC_SUCCESS, "failed to get sequence minmaxmono");
  fail_unless(FC_DBL_EQUIV(min,0.5), "invalid return value");
  rc = fc_getSequenceMinMaxMono(seq,NULL,&minindex,NULL,NULL,NULL);
  fail_unless(rc== FC_SUCCESS, "failed to get sequence minmaxmono");
  fail_unless(minindex ==1, "invalid return value");
  rc = fc_getSequenceMinMaxMono(seq,NULL,NULL,&max,NULL,NULL);
  fail_unless(rc== FC_SUCCESS, "failed to get sequence minmaxmono");
  fail_unless(FC_DBL_EQUIV(max,2.0), "invalid return value");
  rc = fc_getSequenceMinMaxMono(seq,NULL,NULL,NULL,&maxindex,NULL);
  fail_unless(rc== FC_SUCCESS, "failed to get sequence minmaxmono");
  fail_unless(maxindex ==3, "invalid return value");
  rc = fc_getSequenceMinMaxMono(seq,NULL,NULL,NULL,NULL,&mono);
  fail_unless(rc== FC_SUCCESS, "failed to get sequence minmaxmonoreg");
  fail_unless(mono == 0, "invalid return value");

  rc = fc_getSequenceSpacingMinMaxMeanSdev(seq,&ssmin,NULL,NULL,NULL,NULL,NULL);
  fail_unless(rc== FC_SUCCESS, "failed to get sequence spacing minmax");
  fail_unless(FC_DBL_EQUIV(ssmin,-0.5), "invalid return value");
  rc = fc_getSequenceSpacingMinMaxMeanSdev(seq,NULL,&ssminindex,NULL,NULL,NULL,NULL);
  fail_unless(rc== FC_SUCCESS, "failed to get sequence spacing minmax");
  fail_unless(ssminindex ==0, "invalid return value");
  rc = fc_getSequenceSpacingMinMaxMeanSdev(seq,NULL,NULL,&ssmax,NULL,NULL,NULL);
  fail_unless(rc== FC_SUCCESS, "failed to get sequence spacing minmax");
  fail_unless(FC_DBL_EQUIV(ssmax,1.5), "invalid return value");
  rc = fc_getSequenceSpacingMinMaxMeanSdev(seq,NULL,NULL,NULL,&ssmaxindex,NULL,NULL);
  fail_unless(rc== FC_SUCCESS, "failed to get sequence spacing minmax");
  fail_unless(ssmaxindex ==2, "invalid return value");
  rc = fc_getSequenceSpacingMinMaxMeanSdev(seq,NULL,NULL,NULL,NULL,&ssmean,NULL);
  fail_unless(rc== FC_SUCCESS, "failed to get sequence spacing minmax");
  fail_unless(FC_DBL_EQUIV(ssmean,1.0/3.0), "invalid return value");
  rc = fc_getSequenceSpacingMinMaxMeanSdev(seq,NULL,NULL,NULL,NULL,NULL,&sssdev);
  fail_unless(rc== FC_SUCCESS, "failed to get sequence spacing minmax");
  fail_unless(FC_DBL_EQUIV(sssdev,sqrt((2.5-1.0/3.0)/2.0)), "invalid return value");
  rc = fc_deleteSequence(seq);


  //other data types
  for (j = 0; j < 5; j++){
    float spacingx;
    if (j == 0){
      spacingx = 10.;
    }else {
      spacingx/=10;
    }

    fcoords = (float*)malloc(numStep*sizeof(float));
    for (i = 0; i < numStep; i++){
      fcoords[i] = spacingx*(float)(i+1);
    }
    rc= fc_createSequence(global_dataset,"seq",&seq);
    fail_unless(rc == FC_SUCCESS, "couldnt create sequence");
    rc = fc_setSequenceCoords(seq,numStep,FC_DT_FLOAT,(void*)fcoords);
    fail_unless(rc == FC_SUCCESS, "couldnt set sequence coords");
    rc = fc_getSequenceMinMaxMono(seq,&min,&minindex,&max,
				  &maxindex,&mono);
    fail_unless(rc== FC_SUCCESS, "failed to get sequence minmaxmono");
    fail_unless(FC_FLT_EQUIV(max,fcoords[3]), "invalid return value");
    fail_unless(FC_FLT_EQUIV(min,fcoords[0]), "invalid return value");
    fail_unless(maxindex ==3, "invalid return value");
    fail_unless(minindex ==0, "invalid return value");
    fail_unless(mono == 1, "invalid return value");
    
    rc = fc_getClosestSequenceValue(seq,"=",min,&compval,&compindex);
    fail_unless(rc== FC_SUCCESS, "failed to get sequence closest value");
    fail_unless(FC_DBL_EQUIV(min,compval), "invalid return value");
    fail_unless(compindex == minindex, "invalid return index");
    
    rc = fc_getSequenceSpacingMinMaxMeanSdev(seq,&ssmin,&ssminindex,&ssmax,&ssmaxindex,&ssmean,&sssdev);
    fail_unless(rc== FC_SUCCESS, "failed to get sequence spacing minmax");
    //    printf("\nmin %10g minindex %d max %10g maxindex %d mean %10g sdev %10g max-min %10g\n",
    //	   ssmin,ssminindex,ssmax,ssmaxindex,ssmean,sssdev, ssmax-ssmin);
    fail_unless(FC_FLT_EQUIV(ssmin,spacingx), "invalid return value");
    fail_unless(FC_FLT_EQUIV(ssmax,spacingx), "invalid return value");
    //    fail_unless(FC_FLT_EQUIV(ssmax,ssmmin), "invalid return value");  wont work
    fail_unless(FC_VALUE_EQUIV(ssmin,ssmax,10*FLT_EPSILON,FLT_MIN), "invalid return vals");
    //these cant be distinguished to the degree given
    //    fail_unless(ssminindex ==0, "invalid return value");
    //    fail_unless(ssmaxindex ==0, "invalid return value");
    fail_unless(FC_FLT_EQUIV(ssmean,spacingx), "invalid return value");
    fail_unless(FC_FLT_EQUIV(ssmax,ssmean), "invalid return value"); 
    fail_unless(FC_FLT_EQUIV(ssmin,ssmean), "invalid return value"); 
    fail_unless(FC_FLT_EQUIV(sssdev,0.0), "invalid return val" );


    free(fcoords);
    fc_deleteSequence(seq);
  }

  icoords = (int*)malloc(numStep*sizeof(int));
  icoords[0] = 1;
  icoords[1] = 1;
  icoords[2] = 2;
  icoords[3] = 3;

  rc= fc_createSequence(global_dataset,"seq",&seq);
  fail_unless(rc == FC_SUCCESS, "couldnt create sequence");
  rc = fc_setSequenceCoords(seq,numStep,FC_DT_INT,(void*)icoords);
  fail_unless(rc == FC_SUCCESS, "couldnt set sequence coords");
  rc = fc_getSequenceMinMaxMono(seq,&min,&minindex,&max,
				&maxindex,&mono);
  fail_unless(rc== FC_SUCCESS, "failed to get sequence minmaxmono");
  fail_unless((int)max == 3, "invalid return value");
  fail_unless((int)min == 1, "invalid return value");
  fail_unless(maxindex ==3, "invalid return value");
  fail_unless(minindex ==0, "invalid return value");
  fail_unless(mono == 0, "invalid return value");

  rc = fc_getClosestSequenceValue(seq,"=",0.0,&compval,&compindex);
  fail_unless(rc== FC_SUCCESS, "failed to get sequence closest value");
  fail_unless(FC_DBL_EQUIV(min,compval), "invalid return value");
  fail_unless(compindex == minindex, "invalid return index");

  free(icoords);
  fc_deleteSequence(seq);

  //--------------special cases
  // 1 pt seq
  icoords = (int*)malloc(1*sizeof(int));
  icoords[0] = 1.0;
  rc= fc_createSequence(global_dataset,"seq",&seq);
  fail_unless(rc == FC_SUCCESS, "couldnt create sequence");
  rc = fc_setSequenceCoords(seq,1,FC_DT_INT,(void*)icoords);
  fail_unless(rc == FC_SUCCESS, "couldnt set sequence coords");
  rc = fc_getSequenceMinMaxMono(seq,&min,&minindex,&max,
				&maxindex,&mono);
  fail_unless(rc== FC_SUCCESS, "failed to get sequence minmaxmono");
  fail_unless(FC_DBL_EQUIV(max,1.0), "invalid return value");
  fail_unless(FC_DBL_EQUIV(min,1.0), "invalid return value");
  fail_unless(maxindex ==0, "invalid return value");
  fail_unless(minindex ==0, "invalid return value");
  //special return val for 1 pt seq
  fail_unless(mono == 2, "invalid return value");

  rc = fc_getClosestSequenceValue(seq,"=",100,&compval,&compindex);
  fail_unless(rc== FC_SUCCESS, "failed to get sequence closest value");
  fail_unless(FC_DBL_EQUIV(1.0,compval), "invalid return value");
  fail_unless(compindex == 0, "invalid return index");


  rc = fc_getSequenceSpacingMinMaxMeanSdev(seq,&ssmin,&ssminindex,
					   &ssmax,&ssmaxindex,&ssmean,&sssdev);
  fail_unless(rc != FC_SUCCESS, "should fail for single pt seq");
  fail_unless(FC_DBL_EQUIV(ssmin,-1.), "invalid return value");
  fail_unless(FC_DBL_EQUIV(ssmax,-1.), "invalid return value");
  fail_unless(ssminindex == -1, "invalid return value");
  fail_unless(ssmaxindex == -1, "invalid return value");
  fail_unless(FC_DBL_EQUIV(ssmean,-1.), "invalid return value");
  fail_unless(FC_DBL_EQUIV(sssdev,-1.0), "invalid return val" );

  free(icoords);
  fc_deleteSequence(seq);

  //--------------- 2 pt seq
  icoords = (int*)malloc(2*sizeof(int));
  icoords[0] = 1.0;
  icoords[1] = 2.0;
  rc= fc_createSequence(global_dataset,"seq",&seq);
  fail_unless(rc == FC_SUCCESS, "couldnt create sequence");
  rc = fc_setSequenceCoords(seq,2,FC_DT_INT,(void*)icoords);
  fail_unless(rc == FC_SUCCESS, "couldnt set sequence coords");

  rc = fc_getSequenceSpacingMinMaxMeanSdev(seq,&ssmin,&ssminindex,
					   &ssmax,&ssmaxindex,&ssmean,&sssdev);
  fail_unless(rc == FC_SUCCESS, "2pt seq should be valid");
  //  printf("\nmin %10g minindex %d max %10g maxindex %d mean %10g sdev %10g max-min %10g\n",
  //	   ssmin,ssminindex,ssmax,ssmaxindex,ssmean,sssdev, ssmax-ssmin);
  fail_unless(FC_DBL_EQUIV(ssmin,1.), "invalid return value");
  fail_unless(FC_DBL_EQUIV(ssmax,1.), "invalid return value");
  fail_unless(ssminindex == 0, "invalid return value");
  fail_unless(ssmaxindex == 0, "invalid return value");
  fail_unless(FC_DBL_EQUIV(ssmean,1.), "invalid return value");
  fail_unless(FC_DBL_EQUIV(sssdev,0.0), "invalid return val" );

  free(icoords);
  //keep seq for a minute

  //bad args
  rc = fc_getSequenceMinMaxMono(seq,NULL,NULL,NULL,NULL,NULL);
  fail_unless(rc != FC_SUCCESS, "should fail for all NULL return vals");

  rc = fc_getSequenceSpacingMinMaxMeanSdev(seq,NULL,NULL,NULL,NULL,NULL,NULL);
  fail_unless(rc != FC_SUCCESS, "should fail for all NULL return vals");

  rc = fc_getClosestSequenceValue(seq,"=",100,NULL,&compindex);
  fail_unless(rc!= FC_SUCCESS, "should fail for NULL return val");
  fail_unless(compindex == -1, "should return -1 when fail");
  rc = fc_getClosestSequenceValue(seq,"=",100,&compval,NULL);
  fail_unless(rc!= FC_SUCCESS, "should fail for NULL return index");
  fail_unless(FC_DBL_EQUIV(compval,-1), "should return -1 when fail");
  rc = fc_getClosestSequenceValue(seq,NULL,100,&compval,&compindex);
  fail_unless(rc!= FC_SUCCESS, "should fail for NULL comparison criteria");
  fail_unless(FC_DBL_EQUIV(compval,-1), "should return -1 when fail");
  fail_unless(compindex == -1, "should return -1 when fail");
  rc = fc_getClosestSequenceValue(seq,"junk",100,&compval,&compindex);
  fail_unless(rc!= FC_SUCCESS, "should fail for bad comparison criteria");
  fail_unless(FC_DBL_EQUIV(compval,-1), "should return -1 when fail");
  fail_unless(compindex == -1, "should return -1 when fail");

  rc = fc_getSequenceMinMaxMono(FC_NULL_SEQUENCE,&min,&minindex,&max,
				&maxindex,&mono);
  fail_unless(rc!= FC_SUCCESS, "should fail for NULL seq");
  fail_unless(FC_DBL_EQUIV(min,-1.0), "should return -1 when fail");
  fail_unless(FC_DBL_EQUIV(max,-1.0), "should return -1 when fail");
  fail_unless(maxindex ==-1,"should return -1 when fail");
  fail_unless(minindex ==-1,"should return -1 when fail");
  fail_unless(mono == 0, "should return 0 when fail");

  rc = fc_getSequenceSpacingMinMaxMeanSdev(FC_NULL_SEQUENCE,&ssmin,&ssminindex,
					   &ssmax,&ssmaxindex,&ssmean,&sssdev);
  fail_unless(rc != FC_SUCCESS, "should fail for NULL seq");
  fail_unless(FC_DBL_EQUIV(ssmin,-1.), "invalid return value");
  fail_unless(FC_DBL_EQUIV(ssmax,-1.), "invalid return value");
  fail_unless(ssminindex == -1, "invalid return value");
  fail_unless(ssmaxindex == -1, "invalid return value");
  fail_unless(FC_DBL_EQUIV(ssmean,-1.), "invalid return value");
  fail_unless(FC_DBL_EQUIV(sssdev,-1.0), "invalid return val" );

  rc = fc_getClosestSequenceValue(FC_NULL_SEQUENCE,"=",100,&compval,&compindex);
  fail_unless(rc!= FC_SUCCESS, "should fail for NULL seq");
  fail_unless(compval == -1, "should return -1 when fail");
  fail_unless(compindex == -1, "should return -1 when fail");

  fc_deleteSequence(seq);
  rc = fc_getSequenceMinMaxMono(seq,&min,&minindex,&max,
				&maxindex,&mono);
  fail_unless(rc!= FC_SUCCESS, "should fail for invalid seq");
  fail_unless(FC_DBL_EQUIV(min,-1.0), "should return -1 when fail");
  fail_unless(FC_DBL_EQUIV(max,-1.0), "should return -1 when fail");
  fail_unless(maxindex ==-1,"should return -1 when fail");
  fail_unless(minindex ==-1,"should return -1 when fail");
  fail_unless(mono == 0, "should return 0 when fail");

  rc = fc_getSequenceSpacingMinMaxMeanSdev(seq,&ssmin,&ssminindex,
					   &ssmax,&ssmaxindex,&ssmean,&sssdev);
  fail_unless(rc != FC_SUCCESS, "should fail for NULL seq");
  fail_unless(FC_DBL_EQUIV(ssmin,-1.), "invalid return value");
  fail_unless(FC_DBL_EQUIV(ssmax,-1.), "invalid return value");
  fail_unless(ssminindex == -1, "invalid return value");
  fail_unless(ssmaxindex == -1, "invalid return value");
  fail_unless(FC_DBL_EQUIV(ssmean,-1.), "invalid return value");
  fail_unless(FC_DBL_EQUIV(sssdev,-1.0), "invalid return val" );

  rc = fc_getClosestSequenceValue(seq,"=",100,&compval,&compindex);
  fail_unless(rc!= FC_SUCCESS, "should fail for invalid seq");
  fail_unless(compval == -1, "should return -1 when fail");
  fail_unless(compindex == -1, "should return -1 when fail");

  ccoords = (char*)malloc(numStep*sizeof(char));
  for (i = 0; i < numStep; i++){
    ccoords[i] = 'a';
  }
  rc= fc_createSequence(global_dataset,"seq",&seq);
  fail_unless(rc == FC_SUCCESS, "couldnt create sequence");
  rc = fc_setSequenceCoords(seq,numStep,FC_DT_CHAR,(void*)ccoords);
  fail_unless(rc == FC_SUCCESS, "couldnt set sequence coords");
  rc = fc_getSequenceMinMaxMono(seq,&min,&minindex,&max,
				&maxindex,&mono);
  fail_unless(rc!= FC_SUCCESS, "should fail to get char sequence minmaxmono");
  fail_unless(FC_DBL_EQUIV(min,-1.0), "should return -1 when fail");
  fail_unless(FC_DBL_EQUIV(max,-1.0), "should return -1 when fail");
  fail_unless(maxindex ==-1,"should return -1 when fail");
  fail_unless(minindex ==-1,"should return -1 when fail");
  fail_unless(mono == 0, "should return 0 when fail");

  rc = fc_getSequenceSpacingMinMaxMeanSdev(seq,&ssmin,&ssminindex,
					   &ssmax,&ssmaxindex,&ssmean,&sssdev);
  fail_unless(rc != FC_SUCCESS, "should fail for char seq");
  fail_unless(FC_DBL_EQUIV(ssmin,-1.), "invalid return value");
  fail_unless(FC_DBL_EQUIV(ssmax,-1.), "invalid return value");
  fail_unless(ssminindex == -1, "invalid return value");
  fail_unless(ssmaxindex == -1, "invalid return value");
  fail_unless(FC_DBL_EQUIV(ssmean,-1.), "invalid return value");
  fail_unless(FC_DBL_EQUIV(sssdev,-1.0), "invalid return val" );

  rc = fc_getClosestSequenceValue(seq,"=",100,&compval,&compindex);
  fail_unless(rc!= FC_SUCCESS, "should fail for char seq");
  fail_unless(compval == -1, "should return -1 when fail");
  fail_unless(compindex == -1, "should return -1 when fail");

  free(ccoords);
  fc_deleteSequence(seq);

  //other clean up - none

}
END_TEST

  
// Populate the Suite with the tests

Suite *stats_suite(void)
{
  Suite *suite = suite_create("Statistics");

  TCase *tc_variable = tcase_create(" - Variable ");
  TCase *tc_varSubset = tcase_create(" - Variable Subset ");
  TCase *tc_seqVar = tcase_create(" - SeqVariable ");
  TCase *tc_seqVar_series = tcase_create(" - SeqVariable Series");
  TCase *tc_seq = tcase_create(" - Sequence ");

  suite_add_tcase(suite, tc_variable);
  tcase_add_checked_fixture(tc_variable, stats_setup, stats_teardown);
  tcase_add_test(tc_variable, var_stats);

  suite_add_tcase(suite, tc_varSubset);
  tcase_add_checked_fixture(tc_varSubset, stats_setup, stats_teardown);
  tcase_add_test(tc_varSubset, varSubset_stats);

  suite_add_tcase(suite, tc_seqVar);
  tcase_add_checked_fixture(tc_seqVar, stats_setup, stats_teardown);
  tcase_add_test(tc_seqVar, seqVar_stats);

  suite_add_tcase(suite, tc_seqVar_series);
  tcase_add_checked_fixture(tc_seqVar_series, stats_setup, stats_teardown);
  tcase_add_test(tc_seqVar_series, seqVar_series_stats);

  suite_add_tcase(suite, tc_seq);
  tcase_add_checked_fixture(tc_seq, stats_setup, stats_teardown);
  tcase_add_test(tc_seq, seq_stats);

  return suite;
}
