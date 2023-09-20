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
 * \file checkvariable.c
 * \brief Unit tests for data \ref Variable module.
 *
 * $Source: /home/Repositories/fcdmf/fclib/unittest/checkvariable.c,v $
 * $Revision: 1.54 $ 
 * $Date: 2007/02/26 22:13:59 $
 *
 * \modifications
 *    9/10/04 WSK, split off of checklibrary
 */

#include <stdlib.h>
#include <string.h>
#include <check.h>
#include "fc.h"
#include "fcP.h"
#include "checkall.h"

// name of test database created in test fixtures
static char dataset_name[40] = "fixture created dataset";

// **** test fixtures
// setup: init library & create a dataset with a mesh of each element type
// teardown: final library & delete the dataset

static void variable_setup(void) {
  FC_ReturnCode rc;
  int i;
  FC_Dataset dataset;
  FC_Mesh mesh;
  FC_Sequence sequence;
  int numElemType = 8, numElem = 2;
  int numVertPerType[8] = { 2, 3, 4, 6, 5, 6, 8, 12 }; 
  char meshNames[8][40] = { "point mesh", "line mesh", "tri mesh", 
			    "quad mesh", "tet mesh", "pyramid mesh", 
			    "prism mesh", "hex mesh" };
  int numSequence = 2, numSteps[2] = { 5, 10 };
  char seqNames[2][40] = { "time", "time2" };
  FC_ElementType elemTypes[8] = { FC_ET_POINT,  FC_ET_LINE,  FC_ET_TRI,
		      	          FC_ET_QUAD,   FC_ET_TET,   FC_ET_PYRAMID,
			          FC_ET_PRISM,  FC_ET_HEX };
  int conns[8][16] = { { 0,    1 },
		       { 0, 1,    1, 2 },
		       { 0, 1, 2,    1, 3, 2 },
		       { 0, 1, 2, 3,    1, 4, 5, 2 },
		       { 0, 1, 2, 3,    0, 2, 1, 4 },
		       { 0, 1, 2, 3, 4,    0, 3, 2, 1, 5 },
		       { 0, 1, 2, 3, 4, 5,    1, 6, 2, 4, 7, 5 },
		       { 0, 1, 2, 3, 4, 5, 6, 7,    1, 8, 9, 2, 5, 10, 11, 6 } };
  double coords[3*12];
  

  if (isForking) {
    fc_setLibraryVerbosity(fc_messages);
    rc = fc_initLibrary();
    fail_unless(rc == FC_SUCCESS, "test aborted: failed to init library");
  }
  // make a dataset to play in
  rc = fc_createDataset(dataset_name, &dataset);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create dataset");

  // setup coords (these make no physical sense)
  for (i = 0; i < 3*numVertPerType[numElemType-1]; i++)
    coords[i] = i + i/10;

  // make meshes, 1 for each element type
  for (i = 0; i < numElemType; i++) {
    rc = fc_createMesh(dataset, meshNames[i], &mesh);
    fail_unless(rc == FC_SUCCESS, "abort: failed to create mesh");
    rc = fc_setMeshCoords(mesh, 3, numVertPerType[i], coords);
    fail_unless(rc == FC_SUCCESS, "failed to set vertex coords");
    rc = fc_setMeshElementConns(mesh, elemTypes[i], numElem, conns[i]);
    fail_unless(rc == FC_SUCCESS, "failed to set element conns");
  }

  // sequences
  for (i = 0; i < numSequence; i++) {
    rc = fc_createSequence(dataset, seqNames[i], &sequence);
    fail_unless(rc == FC_SUCCESS, "abort: failed to create sequence");
    rc = fc_setSequenceCoords(sequence, numSteps[i], FC_DT_DOUBLE,
			      coords);
    fail_unless(rc == FC_SUCCESS, "abort: failed to set sequence coords");
  }
}

static void variable_teardown(void) {
  FC_ReturnCode rc;
  FC_Dataset *dataset;
  int numDataset;

  rc = fc_getDatasetByName(dataset_name, &numDataset,&dataset);
  fail_unless(rc == FC_SUCCESS && numDataset == 1,
	      "expected to find fixture dataset");
  rc = fc_deleteDataset(dataset[0]);
  free(dataset);
  fail_unless(rc == FC_SUCCESS, 
	      "should not fail to close fixture dataset");

  if (isForking) {
    rc = fc_finalLibrary();
    fail_unless(rc == FC_SUCCESS, "failed to final library at end of test");
  }
}

// --- private variable

// Test that creating and deleting slots modifies table as expected
// Deleting should open up slots that get reused.
START_TEST(slot_new_delete)
{
  FC_ReturnCode rc;
  int i;
  int num = 7;
  // uID_incs = diff between current uID and startID
  int uID_start, uID_incs[7] = { 0, 7, 2, 8, 4, 9, 6 };
  FC_Dataset dataset, *returnDatasets;
  FC_Mesh mesh, *returnMeshes;
  int numReturnDatasets,numReturnMeshes;
  FC_Variable variables[7];
  
  // get dataset and mesh
  rc = fc_getDatasetByName(dataset_name, &numReturnDatasets,
			   &returnDatasets);
  fail_unless(rc == FC_SUCCESS, "aborted - cant get dataset by name");
  fail_unless(numReturnDatasets == 1,
	      "failed to find unique dataset by name");
  dataset = returnDatasets[0];
  free(returnDatasets);

  rc = fc_getMeshByName(dataset, "point mesh", &numReturnMeshes,
			&returnMeshes);
  fail_unless(rc == FC_SUCCESS, "aborted - cant get mesh by name");
  fail_unless(numReturnMeshes == 1,
	      "failed to find unique mesh by name");
  mesh = returnMeshes[0];
  free(returnMeshes);

  // create parent datset & mesh
  //rc = fc_createDataset("fake", &dataset);
  //fail_unless(rc == FC_SUCCESS, "failed to create dataset");
  //rc = fc_createMesh(dataset, "fake", &mesh);
  //fail_unless(rc == FC_SUCCESS, "failed to create mesh");

  // Create num
  for (i = 0; i < num; i++) {
    rc = fc_createVariable(mesh, "fake", &variables[i]);
    fail_unless(rc == FC_SUCCESS, "failed to create variable");
  }
  uID_start = variables[0].uID;

  // Delete "odd" slots
  for (i = 1; i < num; i+=2) {
    rc = fc_deleteVariable(variables[i]);
    fail_unless(rc == FC_SUCCESS, "failed to delete variable");
  }

  // Create some more (fills in the cracks)
  for (i = 1; i < num; i+=2) {
    rc = fc_createVariable(mesh, "fake", &variables[i]);
    fail_unless(rc == FC_SUCCESS, "failed to create variable");
  }

  // Check slotID and uID
  for (i = 0; i < num; i++) {
    fail_unless(variables[i].slotID == i, "mismatch of slot id");
    fail_unless(variables[i].uID == uID_start + uID_incs[i],
		"mismatch of uID");
  }

  // cleanup is done in the fixture
  // cleanup
  //rc = fc_deleteDataset(dataset);
  //fail_unless(rc == FC_SUCCESS, "failed to delete dataset at end of test");
}
END_TEST

// Test that creating and deleting slots modifies table as expected
// Deleting should open up slots that get reused.
START_TEST(slot_new_delete_seq)
{
  FC_ReturnCode rc;
  int i, j;
  int num = 7;
  // uID_incs = diff between current uID and startID
  int uID_start, uID_incs[7] = { 0, 7, 2, 8, 4, 9, 6 };
  FC_Dataset dataset, *returnDatasets;
  FC_Sequence sequence, *returnSequences;
  FC_Mesh mesh, *returnMeshes;
  int numReturnDatasets, numReturnMeshes, numReturnSequences;
  int numStep;
  FC_Variable* seqVars[7];
  
  // get dataset and mesh
  rc = fc_getDatasetByName(dataset_name, &numReturnDatasets,
			   &returnDatasets);
  fail_unless(rc == FC_SUCCESS, "aborted - cant get dataset by name");
  fail_unless(numReturnDatasets == 1,
	      "failed to find unique dataset by name");
  dataset = returnDatasets[0];
  free(returnDatasets);
  rc = fc_getSequenceByName(dataset,"time", &numReturnSequences,
			   &returnSequences);
  fail_unless(rc == FC_SUCCESS, "aborted - cant get sequence by name");
  fail_unless(numReturnSequences == 1,
	      "failed to find unique sequences by name");
  sequence = returnSequences[0];
  free(returnSequences);
  rc = fc_getMeshByName(dataset, "point mesh", &numReturnMeshes,
			&returnMeshes);
  fail_unless(rc == FC_SUCCESS, "aborted - cant get mesh by name");
  fail_unless(numReturnMeshes == 1, "failed to find unique mesh by name");
  mesh = returnMeshes[0];
  free(returnMeshes);
  // create parent datset & mesh
  //rc = fc_createDataset("fake", &dataset);
  //fail_unless(rc == FC_SUCCESS, "failed to create dataset");
  //rc = fc_createMesh(dataset, "fake", &mesh);
  //fail_unless(rc == FC_SUCCESS, "failed to create mesh");

  // Create num
  for (i = 0; i < num; i++) {
    rc = fc_createSeqVariable(mesh, sequence, "fake", &numStep, &seqVars[i]);
    fail_unless(rc == FC_SUCCESS, "failed to create variable");
  }
  uID_start = seqVars[0][0].uID;

  // Delete "odd" slots
  for (i = 1; i < num; i+=2) {
    rc = fc_deleteSeqVariable(numStep, seqVars[i]);
    fail_unless(rc == FC_SUCCESS, "failed to delete variable");
    free(seqVars[i]);
  }
 
  // Create some more (fills in the cracks)
  for (i = 1; i < num; i+=2) {
    rc = fc_createSeqVariable(mesh, sequence, "fake", &numStep,
			      &seqVars[i]);
    fail_unless(rc == FC_SUCCESS, "failed to create variable");
  }
 
  // Check slotID and uID
  for (i = 0; i < num; i++) {
    for (j = 0; j < numStep; j++) {
      fail_unless(seqVars[i][j].slotID == i*numStep+j, "mismatch of slot id");
      fail_unless(seqVars[i][j].uID == uID_start + uID_incs[i]*numStep+j,
		  "mismatch of uID");
    }
  }

  // cleanup
  for (i = 0; i < num; i++)
    free(seqVars[i]);
  // cleanup is done in the fixture
  //rc = fc_deleteDataset(dataset);
  //fail_unless(rc == FC_SUCCESS, "failed to delete dataset at end of test");
}
END_TEST

// -- variable interface

// test creating and deleting variables, also test querying mesh for vars
// for variables before, inbetween and after.
START_TEST(create_get_delete)
{
  FC_ReturnCode rc;
  int i, j;
  FC_Dataset dataset,*returnDatasets;
  FC_Mesh *meshes, mesh;
  int numReturnDatasets,numMesh;
  int numVariable = 10, temp_numVariable;
  FC_Variable variables[10], *temp_variables, temp_variable;
  FC_Variable badVariable = { 999, 999 };
  char names[10][100] = { "one", "two", "three", "four", "five", "six",
			  "seven", "eight", "nine", "ten" };
  char newName[100] = { "sparkly new" };

  // get dataset and meshes
  rc = fc_getDatasetByName(dataset_name, &numReturnDatasets,
			   &returnDatasets);
  fail_unless(rc == FC_SUCCESS, "aborted - cant get dataset by name");
  fail_unless(numReturnDatasets == 1,
	      "failed to find unique dataset by name");
  dataset = returnDatasets[0];
  free(returnDatasets);
  rc = fc_getMeshes(dataset, &numMesh, &meshes);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create mesh");
  fail_unless(numMesh == 8, "abort: expected 8 meshes");

  // loop over all mesh types & association types
  for (i = 0; i < numMesh; i++) { 
    mesh = meshes[i];
    
    // create some variables
    for (j = 0; j < numVariable; j++) {
      rc = fc_createVariable(mesh, names[j], &variables[j]);
      fail_unless(rc == FC_SUCCESS, "failed to create variable");
      fail_unless(!FC_HANDLE_EQUIV(variables[j], FC_NULL_VARIABLE),
		  "created variable should not = FC_NULL_VARIABLE");
    }
    
    // get variables (should be numVariable);
    rc = fc_getVariables(mesh, &temp_numVariable, &temp_variables);
    fail_unless(rc == FC_SUCCESS, "failed to get variables");
    fail_unless(temp_numVariable == numVariable, "mismatch of numVariable");
    for (j = 0; j < numVariable; j++) 
      fail_unless(FC_HANDLE_EQUIV(temp_variables[j], variables[j]), 
		  "mismatch of variable handles");
    free(temp_variables);
    rc = fc_getNumVariable(mesh, &temp_numVariable);
    fail_unless(rc == FC_SUCCESS, "failed to get variables");
    fail_unless(temp_numVariable == numVariable, "mismatch of numVariable");
    for (j = 0; j < numVariable; j++) {
      rc = fc_getVariableByName(mesh, names[j], &temp_numVariable,&temp_variables);
      fail_unless(rc == FC_SUCCESS, "failed to get variable by name");
      fail_unless(temp_numVariable == 1, "failed to get matching variable");
      fail_unless(FC_HANDLE_EQUIV(temp_variables[0], variables[j]),
		  "mismatch of variable handles");
      free(temp_variables);
    }

    //create temporary variable to check multiple variable return
    rc = fc_createVariable(mesh, names[0], &temp_variable);
    fail_unless(rc == FC_SUCCESS, "failed to create variable");
    fail_unless(!FC_HANDLE_EQUIV(temp_variable, FC_NULL_VARIABLE),
		"created variable should not = FC_NULL_VARIABLE");
    rc = fc_getVariableByName(mesh, names[0], &temp_numVariable,&temp_variables);
    fail_unless(rc == FC_SUCCESS, "failed to get variable by name");
    fail_unless(temp_numVariable == 2, "failed to get matching variable");
    fail_unless((FC_HANDLE_EQUIV(temp_variables[0], variables[0]) &&
		 FC_HANDLE_EQUIV(temp_variables[1], temp_variable)) ||
		(FC_HANDLE_EQUIV(temp_variables[1], variables[0]) &&
		 FC_HANDLE_EQUIV(temp_variables[0], temp_variable)),
		"mismatch of variable handles");
    fc_deleteVariable(temp_variable);
    free(temp_variables);
    
    
    // change the name of the first variable
    // also test negative example of fc_getVariableByName
    rc = fc_changeVariableName(variables[0], newName);
    fail_unless(rc == FC_SUCCESS, "failed to change variable name");
    rc = fc_getVariableByName(mesh, names[0], &temp_numVariable,&temp_variables);
    fail_unless(rc == FC_SUCCESS, "old name should succeed");
    fail_unless(temp_numVariable == 0, "old name shouldn't find match");
    fail_unless(temp_variables == NULL, "old name should return NULL");
    rc = fc_getVariableByName(mesh, newName, &temp_numVariable,&temp_variables);
    fail_unless(rc == FC_SUCCESS, "new name should succeed");
    fail_unless(temp_numVariable == 1, "new name should find match");
    fail_unless(FC_HANDLE_EQUIV(temp_variables[0], variables[0]),
		"mismatch of variable handle");
    free(temp_variables);

    // delete half of the variables (alternate)
    for (j = 0; j < numVariable; j+=2) {
      rc = fc_deleteVariable(variables[j]);
      fail_unless(rc == FC_SUCCESS, "failed to delete variable");
    }
    
    // get variables (should be numVariable/2);
    fc_getVariables(mesh, &temp_numVariable, &temp_variables);
    fail_unless(rc == FC_SUCCESS, "failed to get variables");
    fail_unless(temp_numVariable == numVariable/2, "mismatch of numVariable");
    for (j = 0; j < numVariable/2; j++) 
      fail_unless(FC_HANDLE_EQUIV(temp_variables[j], variables[j*2+1]), 
		  "mismatch of variable handles");
    free(temp_variables);
    for (j = 0; j < numVariable/2; j++) {
      rc = fc_getVariableByName(mesh, names[j*2+1], &temp_numVariable,
				&temp_variables);
      fail_unless(rc == FC_SUCCESS, "failed to get variable by name");
      fail_unless(temp_numVariable == 1, "failed to find match");
      fail_unless(FC_HANDLE_EQUIV(temp_variables[0], variables[j*2+1]),
		  "mismatch of variable handles");
      free(temp_variables);
    }
    
    // delete remaining variables
    for (j = 1; j < numVariable; j+=2) { 
      rc = fc_deleteVariable(variables[j]);
      fail_unless(rc == FC_SUCCESS, "failed to delete variable");
    }
    
    // get variables (should be none)
    rc = fc_getVariables(mesh, &temp_numVariable, &temp_variables);
    fail_unless(rc == FC_SUCCESS, 
		"failed to get variables from empty library");
    fail_unless(temp_numVariable == 0 && temp_variables == NULL,
		"should return 0 if all variables deleted");
    
    // make one last variable for further testing
    rc = fc_createVariable(mesh, names[0], &variables[0]);
    fail_unless(rc == FC_SUCCESS, 
		"aborted: failed to create variable");

    // ---- test special cases
    
    // not an error to delete FC_NULL_VARIABLE
    rc = fc_deleteVariable(FC_NULL_VARIABLE);
    fail_unless(rc == FC_SUCCESS, "should not error to delete NULL variable");
    
    // ---- test error conditions
    
    // bad args to fc_createVariable()
    temp_variable = badVariable;
    rc = fc_createVariable(FC_NULL_MESH, names[0], &temp_variable);
    fail_unless(rc != FC_SUCCESS, 
		"should have failed to create variable with null database");
    fail_unless(FC_HANDLE_EQUIV(temp_variable, FC_NULL_VARIABLE),
		"fail should return NULL variable");
    temp_variable = badVariable;
    rc = fc_createVariable(mesh, NULL, &temp_variable);
    fail_unless(rc != FC_SUCCESS, 
		"should have failed to create variable with null name");
    fail_unless(FC_HANDLE_EQUIV(temp_variable, FC_NULL_VARIABLE),
		"fail should return NULL variable");
    rc = fc_createVariable(mesh, names[0], NULL);
    fail_unless(rc != FC_SUCCESS, 
		"should have failed to create variable with null handle");
    
    // bad args to fc_changeVariableName()
    rc = fc_changeVariableName(variables[0], NULL);
    fail_unless(rc != FC_SUCCESS, "should fail if new name is NULL");
    rc = fc_getVariableByName(mesh, names[0], &temp_numVariable,
				&temp_variables);
    fail_unless(rc == FC_SUCCESS && temp_numVariable == 1 &&
		FC_HANDLE_EQUIV(variables[0], temp_variables[0]),
		"should not change name of variable");
    free(temp_variables);

     // bad args to fc_getVariables()
    temp_numVariable = 99;
    temp_variables = (FC_Variable*)1;
    rc = fc_getVariables(FC_NULL_MESH, &temp_numVariable, &temp_variables);
    fail_unless(rc != FC_SUCCESS, 
		"should have failed to create variable with null mesh");
    fail_unless(temp_numVariable == -1,
		"fail should return -1 for numVariable"); 
    fail_unless(temp_variables == NULL, "fail should return NULL");
    temp_variables = (FC_Variable*)1;
    rc = fc_getVariables(mesh, NULL, &temp_variables);
    fail_unless(rc != FC_SUCCESS, 
		"should have failed to get variable with NULL numVariable");
    fail_unless(temp_variables == NULL, "fail should return null");
    temp_numVariable = 999;
    rc = fc_getVariables(mesh, &temp_numVariable, NULL);
    fail_unless(rc != FC_SUCCESS, "should be error to get just numVariable");
    fail_unless(temp_numVariable == -1, "mismatch of numVariable");    
    
    // bad args to fc_getNumVariable()
    temp_numVariable = 99;
    rc = fc_getNumVariable(FC_NULL_MESH, &temp_numVariable);
    fail_unless(rc != FC_SUCCESS, 
		"should have failed to get numVariable with null mesh");
    fail_unless(temp_numVariable == -1,
		"fail should return -1 for numVariable"); 
    rc = fc_getNumVariable(mesh, NULL);
    fail_unless(rc != FC_SUCCESS, 
		"should have failed to get numVariable with NULL numVariable");
    
    // bad args to fc_getVariableByName()
    rc = fc_getVariableByName(FC_NULL_MESH, names[0], &temp_numVariable,
			      &temp_variables);
    fail_unless(rc != FC_SUCCESS, 
		"should have failed to get variable with null mesh");
    fail_unless(temp_numVariable == -1, "should return -1 when fail");
    fail_unless(temp_variables == NULL, "should return NULL when fail");
    rc = fc_getVariableByName(mesh, NULL, &temp_numVariable,
			      &temp_variables);
    fail_unless(rc != FC_SUCCESS, 
		"should have failed to get variable with null name");
    fail_unless(temp_numVariable == -1, "should return -1 when fail");
    fail_unless(temp_variables == NULL, "should return NULL when fail");
    rc = fc_getVariableByName(mesh, names[0], NULL, &temp_variables);
    fail_unless(rc != FC_SUCCESS, 
		"should have failed to get variable with null arg");
    fail_unless(temp_variables == NULL, "should return NULL when fail");
    rc = fc_getVariableByName(mesh, names[0], &temp_numVariable, NULL);
    fail_unless(rc != FC_SUCCESS, 
		"should have failed to get variable with null arg");
    fail_unless(temp_numVariable == -1, "should return -1 when fail");
    
    // bad args to fc_deleteVariable()
    rc = fc_deleteVariable(badVariable);
    fail_unless(rc != FC_SUCCESS, 
		"should error to delete nonexistent variable");

    // --- done

    // delete last variable
    rc = fc_deleteVariable(variables[0]);
    fail_unless(rc == FC_SUCCESS, "failed to delete last variable");
  }

  free(meshes);
}
END_TEST

// query meta data
START_TEST(metadata_query)
{
  FC_ReturnCode rc;
  int i;
  FC_Dataset dataset, temp_dataset, badDataset = { 999, 999 };
  FC_Dataset *returnDatasets;
  int numReturnDatasets,numMesh;
  FC_Mesh *meshes, mesh, temp_mesh, badMesh = {999, 999 };
  char name[20] = "blue berry", *temp_name;
  FC_Variable variable, badVariable = { 999, 999 };

  // get dataset and meshes
  rc = fc_getDatasetByName(dataset_name, &numReturnDatasets,
			   &returnDatasets);
  fail_unless(rc == FC_SUCCESS, "aborted - cant get dataset by name");
  fail_unless(numReturnDatasets == 1,
	      "failed to find unique dataset by name");
  dataset = returnDatasets[0];
  free(returnDatasets);
  rc = fc_getMeshes(dataset, &numMesh, &meshes);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create mesh");
  fail_unless(numMesh == 8, "abort: expected 8 meshes");

  // loop over all mesh types & association types
  for (i = 0; i < numMesh; i++) { 
    mesh = meshes[i];

    // create a variable to play with
    rc = fc_createVariable(mesh, name, &variable);
    fail_unless(rc == FC_SUCCESS, "failed to create variable");
    
    // test fc_isVariableValid()
    fail_unless(fc_isVariableValid(variable),
		"failed to validate a valid variable");
    
    // test fc_isVariableGlobal()
    fail_unless(!fc_isVariableGlobal(variable),
		"variable should be non-global");

    // test fc_getVariableName();
    temp_name = NULL;
    rc = fc_getVariableName(variable, &temp_name);
    fail_unless(rc == FC_SUCCESS, "failed to get variable name");
    fail_unless(!strcmp(name, temp_name), "mismatch of name");
    // should have returned copy, mess with temp_name & try again to make sure
    temp_name[0] = 'q';
    free(temp_name);
    rc = fc_getVariableName(variable, &temp_name);
    fail_unless(rc == FC_SUCCESS, "failed to get variable name");
    fail_unless(!strcmp(name, temp_name), "mismatch of name");
    free(temp_name);
    
    // test fc_getDatasetFromVariable()
    temp_dataset = badDataset;
    rc = fc_getDatasetFromVariable(variable, &temp_dataset);
    fail_unless(rc == FC_SUCCESS, "failed to get parent dataset");
    fail_unless(FC_HANDLE_EQUIV(temp_dataset, dataset), 
		"mismatch of parent dataset");

    // test fc_getMeshFromVariable()
    temp_mesh = badMesh;
    rc = fc_getMeshFromVariable(variable, &temp_mesh);
    fail_unless(rc == FC_SUCCESS, "failed to get parent mesh");
    fail_unless(FC_HANDLE_EQUIV(temp_mesh, mesh), "mismatch of parent mesh");
    
    // --- check with bad args
    
    // fc_isVariableValid()
    fail_unless(!fc_isVariableValid(badVariable), 
		"badVariable should NOT be valid");
    fail_unless(!fc_isVariableValid(FC_NULL_VARIABLE), 
		"FC_NULL_VARIABLE should not be valid");

    // fc_isVariableGlobal()
    fail_unless(fc_isVariableGlobal(badVariable) < FC_SUCCESS,
		"badVariable returns an error");
    
    // fc_getVariableName()
    temp_name = (char*)1;
    rc = fc_getVariableName(badVariable, &temp_name);
    fail_unless(rc != FC_SUCCESS, "badVariable should NOT return a name");
    fail_unless(temp_name == NULL, "fail should return NULL name");
    rc = fc_getVariableName(variable, NULL);
    fail_unless(rc != FC_SUCCESS, "should fail for NULL name");
    
    // fc_getDatasetFromVariable()
    temp_dataset = badDataset;
    rc = fc_getDatasetFromVariable(badVariable, &temp_dataset);
    fail_unless(rc != FC_SUCCESS, "badVariable should NOT return a dataset");
    fail_unless(FC_HANDLE_EQUIV(temp_dataset, FC_NULL_DATASET),
		"failure should return NULL dataset");
    rc = fc_getDatasetFromVariable(variable, NULL);
    fail_unless(rc != FC_SUCCESS, "should fail for NULL dataset");

    // fc_getMeshFromVariable()
    temp_mesh = badMesh;
    rc = fc_getMeshFromVariable(badVariable, &temp_mesh);
    fail_unless(rc != FC_SUCCESS, "badVariable should NOT return a mesh");
    fail_unless(FC_HANDLE_EQUIV(temp_mesh, FC_NULL_MESH),
		"failure should return NULL mesh");
    rc = fc_getMeshFromVariable(variable, NULL);
    fail_unless(rc != FC_SUCCESS, "should fail for NULL mesh");
    
    // --- all done
    
    // one last test
    fc_deleteVariable(variable);
    fail_unless(!fc_isVariableValid(variable), 
		"handle should not be valid after a delete");
  }

  free(meshes);
}
END_TEST

// query the variable data - all assocs, mathtypes & datatypes
START_TEST(data_query)
{
  FC_ReturnCode rc;
  int i, j, k, m;
  char name[20] = "blue berry", name2[20] = "banana nut";
  FC_Dataset dataset, *returnDatasets;
  int numReturnDatasets,numMesh;
  FC_Mesh *meshes, mesh;
  int maxNumEntity = 20; // for all the meshes and all assoc types
  FC_Variable variable, variable2, empty_variable, badVariable = { 999, 999 };
  FC_Variable global_var;
  _FC_VarSlot* varSlot;
  int numAssoc = 5, numDataPoint;  // numDataPoint determined later
  FC_AssociationType assocs[5] = { FC_AT_VERTEX, FC_AT_EDGE, FC_AT_FACE,
				   FC_AT_ELEMENT, FC_AT_WHOLE_MESH };
  int numMathType = 4, numComps[4] = { 1, 2, 3, 4 }, maxNumComp = 4;
  FC_MathType mathTypes[4] = { FC_MT_SCALAR, FC_MT_VECTOR, 
			       FC_MT_SYMTENSOR, FC_MT_TENSOR };
  int numDataType = 4;
  FC_DataType dataTypes[4] = { FC_DT_CHAR, FC_DT_INT, 
			       FC_DT_FLOAT, FC_DT_DOUBLE };
  char charData[20*4] = { 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j',
			'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't' };
  int intData[20*4];
  float floatData[20*4];
  double doubleData[20*4];
  int sizes[4] = {sizeof(char), sizeof(int), sizeof(float), sizeof(double) };
  void *data[4] = { charData, intData, floatData, doubleData };
  int temp_numDataPoint, temp_numComp;
  FC_AssociationType temp_assoc;
  FC_MathType temp_mathType;
  FC_DataType temp_dataType;
  void *temp_data, *data_cpy;
  
  // get dataset and meshes & global var
  rc = fc_getDatasetByName(dataset_name, &numReturnDatasets,
			   &returnDatasets);
  fail_unless(rc == FC_SUCCESS, "aborted - cant get dataset by name");
  fail_unless(numReturnDatasets == 1,
	      "failed to find unique dataset by name");
  dataset = returnDatasets[0];
  free(returnDatasets);
  rc = fc_getMeshes(dataset, &numMesh, &meshes);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create mesh");
  fail_unless(numMesh == 8, "abort: expected 8 meshes");
  rc = fc_createGlobalVariable(dataset, "global", &global_var);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create global var");

  // setup some data
  for (i = 0; i < maxNumEntity*maxNumComp; i++) {
    if (i > 0)
      charData[i] = charData[i%maxNumEntity];
    intData[i] = i;
    floatData[i] = intData[i] + 0.1;
    doubleData[i] = intData[i] + 0.00000000001;
  }

  // loop over all mesh types, association types, math types & data types
  for (i = 0; i < numMesh; i++) { 

    mesh = meshes[i];

    for (j = 0; j < numAssoc; j++) {

      // setup numDataPoint
      fc_getMeshNumEntity(mesh, assocs[j], &numDataPoint);

      for (k = 0; k < numMathType; k++) {
	for (m = 0; m < numDataType; m++) {

	  //printf("mesh %d: assoc %s, mathtype %s, datatype %s\n", i,
	  //     fc_getAssociationTypeText(assocs[j]),
	  //     fc_getMathTypeText(mathType[k]),
	  //     fc_getDataTypeText(dataTypes[m]));

	  // create three variables
	  rc = fc_createVariable(mesh, name, &variable);
	  fail_unless(rc == FC_SUCCESS, "abort: failed to create variable");
	  rc = fc_createVariable(mesh, name2, &variable2);
	  fail_unless(rc == FC_SUCCESS, "abort: failed to create variable");
	  rc = fc_createVariable(mesh, "empty variable", &empty_variable);
	  fail_unless(rc == FC_SUCCESS, "abort: failed to create variable");
	  
	  // check default values
	  rc = fc_getVariableInfo(variable, &temp_numDataPoint, &temp_numComp, 
				  &temp_assoc, &temp_mathType, &temp_dataType);
	  fail_unless(rc == FC_SUCCESS, "failed to get variable info");
	  fail_unless(temp_numDataPoint == 0 && temp_numComp == 0 && 
		      temp_assoc == FC_AT_UNKNOWN && temp_mathType ==
		      FC_MT_UNKNOWN && temp_dataType == FC_DT_UNKNOWN,
		      "empty variable should have empty values");
	  
	  // setting global vars should fail
	  if (assocs[j] != FC_AT_WHOLE_MESH) {
	    rc = fc_setVariableData(global_var, numDataPoint, numComps[k],
				    assocs[j], mathTypes[k], dataTypes[m],
				    data[m]);
	    fail_unless(rc < FC_SUCCESS, 
			"should fail to set global var w/ nonglobal assoc");
	  }

	  // variable1: set/get data by copy and check
	  temp_numDataPoint = 999;
	  temp_dataType = FC_DT_UNKNOWN;
	  temp_data = NULL;
	  rc = fc_setVariableData(variable, numDataPoint, numComps[k],
				  assocs[j],  mathTypes[k], dataTypes[m],
				  data[m]);
	  if ( (i == 0 && (j == 1 || j == 2)) // point mesh - no edges,faces
	       || (i == 1 && j == 2) ) { // line mesh - no faces
	    fail_unless(rc != FC_SUCCESS, "should fail to create variable");
	    fc_deleteVariable(variable);
	    continue; // early end to this loop
	  }
	  else
	    fail_unless(rc == FC_SUCCESS, "failed to set data");
	  rc = fc_getVariableInfo(variable, &temp_numDataPoint, &temp_numComp, 
				  &temp_assoc, &temp_mathType, &temp_dataType);
	  fail_unless(rc == FC_SUCCESS, "failed to get variable info");
	  fail_unless(temp_numDataPoint == numDataPoint, 
		      "mismatch of numDataPoint");
	  fail_unless(temp_numComp == numComps[k], "mismatch of numComp");
	  fail_unless(temp_assoc == assocs[j], "mismatch of assoc");
	  fail_unless(temp_mathType == mathTypes[k], "mismatch of mathType");
	  fail_unless(temp_dataType == dataTypes[m], "mismatch of datatype");
	  rc = fc_getVariableDataPtr(variable, &temp_data);
	  fail_unless(rc == FC_SUCCESS, "failed to get data");
	  fail_unless(!memcmp(temp_data, data[m], 
			      numDataPoint*numComps[k]*sizes[m]),
		      "mismatch of data");

	  // variable2: set/get data pointer and check
	  data_cpy = malloc(numDataPoint*numComps[k]*sizes[m]);
	  memcpy(data_cpy, data[m], numDataPoint*numComps[k]*sizes[m]);
	  temp_numDataPoint = 999;
	  temp_dataType = FC_DT_UNKNOWN;
	  temp_data = NULL;
	  rc = fc_setVariableDataPtr(variable2, numDataPoint, numComps[k],
				  assocs[j],  mathTypes[k], dataTypes[m],
				  data_cpy);
	  if ( (i == 0 && (j == 1 || j == 2)) // point mesh - no edges,faces
	       || (i == 1 && j == 2) ) { // line mesh - no faces
	    fail_unless(rc != FC_SUCCESS, "should fail to create variable");
	    fc_deleteVariable(variable);
	    continue; // early end to this loop
	  }
	  else
	    fail_unless(rc == FC_SUCCESS, "failed to set data");
	  rc = fc_getVariableInfo(variable2, &temp_numDataPoint, &temp_numComp, 
				  &temp_assoc, &temp_mathType, &temp_dataType);
	  fail_unless(rc == FC_SUCCESS, "failed to get variable info");
	  fail_unless(temp_numDataPoint == numDataPoint, 
		      "mismatch of numDataPoint");
	  fail_unless(temp_numComp == numComps[k], "mismatch of numComp");
	  fail_unless(temp_assoc == assocs[j], "mismatch of assoc");
	  fail_unless(temp_mathType == mathTypes[k], "mismatch of mathType");
	  fail_unless(temp_dataType == dataTypes[m], "mismatch of datatype");
	  rc = fc_getVariableDataPtr(variable2, &temp_data);
	  fail_unless(rc == FC_SUCCESS, "failed to get data");
	  fail_unless(!memcmp(temp_data, data[m], 
			      numDataPoint*numComps[k]*sizes[m]),
		      "mismatch of data");

	  // test set by copy vs. set the pointer
	  varSlot = _fc_getVarSlot(variable);
	  fail_unless(varSlot->data != data[m], 
		      "should have copied the coords");
	  varSlot = _fc_getVarSlot(variable2);
	  fail_unless(varSlot->data == data_cpy, "should have ptr to coords");

	  // check other query simpler forms of fc_getVariableInfo
	  rc = fc_getVariableNumDataPoint(variable, &temp_numDataPoint);
	  fail_unless(rc == FC_SUCCESS, "failed to get numDataPoint");
	  fail_unless(temp_numDataPoint == numDataPoint, 
		      "mismatch of numDataPoint");
	  rc = fc_getVariableNumComponent(variable, &temp_numComp);
	  fail_unless(rc == FC_SUCCESS && temp_numComp == numComps[k],
		      "failed to get numComponent");
	  rc = fc_getVariableAssociationType(variable, &temp_assoc);
	  fail_unless(rc == FC_SUCCESS && temp_assoc == assocs[j],
		      "failed to get assoc");
	  rc = fc_getVariableMathType(variable, &temp_mathType);
	  fail_unless(rc == FC_SUCCESS && temp_mathType == mathTypes[k],
		      "failed to get mathType");
	  rc = fc_getVariableDataType(variable, &temp_dataType);
	  fail_unless(rc == FC_SUCCESS, "failed to get data type");
	  fail_unless(temp_dataType == dataTypes[m], "mismatch of datatype");
	  
	  // ---- check special cases
	  
	  // fc_getVariableInfo -- most arguments are optional
	  rc = fc_getVariableInfo(variable, &temp_numDataPoint, NULL, NULL, 
				  NULL, NULL);
	  fail_unless(rc == FC_SUCCESS, "failed to get variable info");
	  fail_unless(temp_numDataPoint == numDataPoint, 
		      "mismatch of numDataPoint");
	  rc = fc_getVariableInfo(variable, NULL, &temp_numComp, NULL, NULL,
				  NULL); 
	  fail_unless(rc == FC_SUCCESS, "failed to get variable info");
	  fail_unless(temp_numComp == numComps[k], "mismatch of numComp");
	  rc = fc_getVariableInfo(variable, NULL, NULL, &temp_assoc, NULL, 
				  NULL);
	  fail_unless(rc == FC_SUCCESS, "failed to get variable info");
	  fail_unless(temp_assoc == assocs[j], "mismatch of assoc");
	  rc = fc_getVariableInfo(variable, NULL, NULL, NULL, &temp_mathType,
				  NULL); 
	  fail_unless(rc == FC_SUCCESS, "failed to get variable info");
	  fail_unless(temp_mathType == mathTypes[k], "mismatch of mathType");
	  rc = fc_getVariableInfo(variable, NULL, NULL, NULL, NULL,
				  &temp_dataType); 
	  fail_unless(rc == FC_SUCCESS, "failed to get variable info");
	  fail_unless(temp_dataType == dataTypes[m], "mismatch of datatype");
	  
	  // ---- check errors
	  
	  // fc_getVariableInfo() -- bad args
	  rc = fc_getVariableInfo(badVariable, &temp_numDataPoint, 
				  &temp_numComp, &temp_assoc, &temp_mathType,
				  &temp_dataType);
	  fail_unless(rc != FC_SUCCESS, "should fail with bad variable");
	  fail_unless(temp_numDataPoint == -1 && temp_numComp == -1 && 
		      temp_assoc == FC_AT_UNKNOWN && 
		      temp_mathType == FC_MT_UNKNOWN && 
		      temp_dataType == FC_DT_UNKNOWN,
		      "fail should return NULLs");
	  rc = fc_getVariableInfo(variable, NULL, NULL, NULL, NULL, NULL);
	  fail_unless(rc != FC_SUCCESS, "can't have all args be NULL");

	  // functions similar to fc_getVariableInfo() -- bad args
	  rc = fc_getVariableNumDataPoint(badVariable, &temp_numDataPoint);
	  fail_unless(rc != FC_SUCCESS && temp_numDataPoint == -1, 
		      "should fail with bad variable");
	  rc = fc_getVariableNumDataPoint(variable, NULL);
	  fail_unless(rc != FC_SUCCESS, "should faile with null numDataPoint");
	  rc = fc_getVariableNumComponent(badVariable, &temp_numComp);
	  fail_unless(rc != FC_SUCCESS && temp_numComp == -1, 
		      "should fail with bad variable");
	  rc = fc_getVariableNumComponent(variable, NULL);
	  fail_unless(rc != FC_SUCCESS, "should faile with null numComponent");
	  rc = fc_getVariableAssociationType(badVariable, &temp_assoc);
	  fail_unless(rc != FC_SUCCESS && temp_assoc == FC_AT_UNKNOWN,
		      "should fail with bad variable");
	  rc = fc_getVariableAssociationType(variable, NULL);
	  fail_unless(rc != FC_SUCCESS, "should faile with null assoc");
	  rc = fc_getVariableMathType(badVariable, &temp_mathType);
	  fail_unless(rc != FC_SUCCESS && temp_mathType == FC_MT_UNKNOWN,
		      "should fail with bad variable");
	  rc = fc_getVariableMathType(variable, NULL);
	  fail_unless(rc != FC_SUCCESS, "should faile with null mathType");
	  rc = fc_getVariableDataType(badVariable, &temp_dataType);
	  fail_unless(rc != FC_SUCCESS && temp_dataType == FC_DT_UNKNOWN,
		      "should fail with bad variable");
	  rc = fc_getVariableDataType(variable, NULL);
	  fail_unless(rc != FC_SUCCESS, "should faile with null dataType");
	  
	  // fc_getVariableDataPtr() -- bad args
	  rc = fc_getVariableDataPtr(badVariable, &temp_data);
	  fail_unless(rc != FC_SUCCESS, "should fail with bad variable");
	  fail_unless(temp_data == NULL, "fail should return nulls");
	  rc = fc_getVariableDataPtr(variable, NULL);
	  fail_unless(rc != FC_SUCCESS, "should fail if data = NULL");
	  
	  // fc_setVariableData() -- error to try to reset
	  rc = fc_setVariableData(variable, numDataPoint, numComps[k],
				  assocs[j], mathTypes[k], dataTypes[m], 
				  data[m]);
	  fail_unless(rc != FC_SUCCESS, "shouldn't be able to reset data");
	  
	  // fc_setVariableData() -- bad args
	  rc = fc_setVariableData(badVariable, numDataPoint, numComps[k], 
				  assocs[j], mathTypes[k], dataTypes[m], 
				  data[m]);
	  fail_unless(rc != FC_SUCCESS, "should fail with bad variable");
	  rc = fc_setVariableData(empty_variable, 0, numComps[k], assocs[j], 
				  mathTypes[k], dataTypes[m], data[m]);
	  fail_unless(rc != FC_SUCCESS, "should fail if numDataPoint < 1");
	  rc = fc_setVariableData(empty_variable, numDataPoint + 1, 
				  numComps[k], assocs[j], mathTypes[k], 
				  dataTypes[m], data[m]);
	  fail_unless(rc != FC_SUCCESS,
		      "should fail if numDataPoint & assoc don't agree");
	  rc = fc_setVariableData(empty_variable, numDataPoint, 0, assocs[j],
				  mathTypes[k], dataTypes[m], data[m]);
	  fail_unless(rc != FC_SUCCESS, "should fail if numComponent < 1");
	  if (mathTypes[k] == FC_MT_SCALAR) {
	    rc = fc_setVariableData(empty_variable, numDataPoint, 2, assocs[j],
				    mathTypes[k], dataTypes[m], data[m]);
	    fail_unless(rc != FC_SUCCESS, 
			"should fail if numComponent & mathtype don't agree");
	  }
	  else {
	    rc = fc_setVariableData(empty_variable, numDataPoint, 1, assocs[j],
				    mathTypes[k], dataTypes[m], data[m]);
	    fail_unless(rc != FC_SUCCESS, 
			"should fail if numComponent & mathtype don't agree");
	  }
	  rc = fc_setVariableData(empty_variable, numDataPoint, numComps[k], 
				  FC_AT_UNKNOWN, mathTypes[k], dataTypes[m], 
				  data[m]);
	  fail_unless(rc != FC_SUCCESS, "should fail with unknown assoc");
	  rc = fc_setVariableData(empty_variable, numDataPoint, numComps[k], 
				  -999, mathTypes[k], dataTypes[m], data[m]);
	  fail_unless(rc != FC_SUCCESS, "should fail with invalid assoc");
	  rc = fc_setVariableData(empty_variable, numDataPoint, numComps[k], 
				  assocs[j], FC_MT_UNKNOWN, dataTypes[m], 
				  data[m]);
	  fail_unless(rc != FC_SUCCESS, "should fail with unknown mathType");
	  rc = fc_setVariableData(empty_variable, numDataPoint, numComps[k], 
				  assocs[j], -999, dataTypes[m], data[m]);
	  fail_unless(rc != FC_SUCCESS, "should fail with invalid mathType");
	  rc = fc_setVariableData(empty_variable, numDataPoint, numComps[k], 
				  assocs[j], mathTypes[k], FC_DT_UNKNOWN, 
				  data[m]);
	  fail_unless(rc != FC_SUCCESS, "should fail with unknown dateType");
	  rc = fc_setVariableData(empty_variable, numDataPoint, numComps[k], 
				  assocs[j], mathTypes[k], -999, data[m]);
	  fail_unless(rc != FC_SUCCESS, "should fail with invalid dateType");
	  rc = fc_setVariableData(empty_variable, numDataPoint, numComps[k], 
				  assocs[j], mathTypes[k], dataTypes[m], NULL);
	  fail_unless(rc != FC_SUCCESS, "should fail with no data");

	  // fc_setVariableDataPtr() -- bad args
	  rc = fc_setVariableDataPtr(badVariable, numDataPoint, numComps[k], 
				  assocs[j], mathTypes[k], dataTypes[m], 
				  data[m]);
	  fail_unless(rc != FC_SUCCESS, "should fail with bad variable");
	  rc = fc_setVariableDataPtr(empty_variable, 0, numComps[k], assocs[j],
				  mathTypes[k], dataTypes[m], data[m]);
	  fail_unless(rc != FC_SUCCESS, "should fail if numDataPoint < 1");
	  rc = fc_setVariableDataPtr(empty_variable, numDataPoint + 1, 
				     numComps[k], assocs[j], mathTypes[k], 
				     dataTypes[m], data[m]);
	  fail_unless(rc != FC_SUCCESS,
		      "should fail if numDataPoint & assoc don't agree");
	  rc = fc_setVariableDataPtr(empty_variable, numDataPoint, 0, assocs[j],
				  mathTypes[k], dataTypes[m], data[m]);
	  fail_unless(rc != FC_SUCCESS, "should fail if numComponent < 1");
	  if (mathTypes[k] == FC_MT_SCALAR) {
	    rc = fc_setVariableDataPtr(empty_variable, numDataPoint, 2, 
				       assocs[j], mathTypes[k], dataTypes[m], 
				       data[m]);
	    fail_unless(rc != FC_SUCCESS, 
			"should fail if numComponent & mathtype don't agree");
	  }
	  else {
	    rc = fc_setVariableDataPtr(empty_variable, numDataPoint, 1, 
				       assocs[j], mathTypes[k], dataTypes[m], 
				       data[m]);
	    fail_unless(rc != FC_SUCCESS, 
			"should fail if numComponent & mathtype don't agree");
	  }
	  rc = fc_setVariableDataPtr(empty_variable, numDataPoint, numComps[k],
				  FC_AT_UNKNOWN, mathTypes[k], dataTypes[m], 
				  data[m]);
	  fail_unless(rc != FC_SUCCESS, "should fail with unknown assoc");
	  rc = fc_setVariableDataPtr(empty_variable, numDataPoint, numComps[k],
				     -999, mathTypes[k], dataTypes[m], 
				     data[m]);
	  fail_unless(rc != FC_SUCCESS, "should fail with invalid assoc");
	  rc = fc_setVariableDataPtr(empty_variable, numDataPoint, numComps[k],
				  assocs[j], FC_MT_UNKNOWN, dataTypes[m], 
				  data[m]);
	  fail_unless(rc != FC_SUCCESS, "should fail with unknown mathType");
	  rc = fc_setVariableDataPtr(empty_variable, numDataPoint, numComps[k],
				  assocs[j], -999, dataTypes[m], data[m]);
	  fail_unless(rc != FC_SUCCESS, "should fail with invalid mathType");
	  rc = fc_setVariableDataPtr(empty_variable, numDataPoint, numComps[k],
				  assocs[j], mathTypes[k], FC_DT_UNKNOWN, 
				  data[m]);
	  fail_unless(rc != FC_SUCCESS, "should fail with unknown dateType");
	  rc = fc_setVariableDataPtr(empty_variable, numDataPoint, numComps[k],
				  assocs[j], mathTypes[k], -999, data[m]);
	  fail_unless(rc != FC_SUCCESS, "should fail with invalid dateType");
	  rc = fc_setVariableDataPtr(empty_variable, numDataPoint, numComps[k],
				  assocs[j], mathTypes[k], dataTypes[m], NULL);
	  fail_unless(rc != FC_SUCCESS, "should fail with no data");

	  // empty variable should still be empty (all set's failed)
	  rc = fc_getVariableInfo(empty_variable, &temp_numDataPoint, 
				  &temp_numComp, &temp_assoc, &temp_mathType, 
				  &temp_dataType);
	  fail_unless(rc == FC_SUCCESS, "failed to get variable info");
	  fail_unless(temp_numDataPoint == 0 && temp_numComp == 0 && 
		      temp_assoc == FC_AT_UNKNOWN && temp_mathType ==
		      FC_MT_UNKNOWN && temp_dataType == FC_DT_UNKNOWN,
		      "empty variable should have empty values");

	  // --- all done
	  // cleanup
	  fc_deleteVariable(variable);
	}
      }
    }
  }
  free(meshes);
}
END_TEST

// Request variable data be returned as a specific type
START_TEST(get_data_as)
{
  FC_ReturnCode rc;
  int i, j, k;
  char name[20] = "blue berry";
  FC_Dataset dataset, *returnDatasets;
  FC_Mesh mesh, *returnMeshes;
  int numReturnDatasets,numReturnMeshes;
  FC_AssociationType assoc = FC_AT_VERTEX;
  int numDataPoint; // numDataPoint determined later
  FC_Variable variables[4], badVariable = { 999, 999 };
  int numMathType = 4, numComps[4] = { 1, 2, 3, 4 }, maxNumComp = 4;
  FC_MathType mathTypes[4] = { FC_MT_SCALAR, FC_MT_VECTOR, 
			       FC_MT_SYMTENSOR, FC_MT_TENSOR };
  int numDataType = 4;
  FC_DataType dataTypes[4] = { FC_DT_CHAR, FC_DT_INT, 
			       FC_DT_FLOAT, FC_DT_DOUBLE };
  char charData[20*4] = { 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j',
			'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't' };
  int intData[20*4];
  float floatData[20*4];
  double doubleData[20*4];
  void *data[4] = { charData, intData, floatData, doubleData };
  void *temp_data;
  
  // get dataset and mesh
  rc = fc_getDatasetByName(dataset_name, &numReturnDatasets,
			   &returnDatasets);
  fail_unless(rc == FC_SUCCESS, "aborted - cant get dataset by name");
  fail_unless(numReturnDatasets == 1,
	      "failed to find unique dataset by name");
  dataset = returnDatasets[0];
  free(returnDatasets);
  rc = fc_getMeshByName(dataset, "hex mesh", &numReturnMeshes,
			&returnMeshes);
  fail_unless(rc == FC_SUCCESS, "aborted - cant get mesh by name");
  fail_unless(numReturnMeshes == 1, "failed to find unique mesh by name");
  mesh = returnMeshes[0];
  free(returnMeshes);

  // setup numDataPoint
  fc_getMeshNumEntity(mesh, assoc, &numDataPoint);

  // setup some data
  for (i = 0; i < numDataPoint*maxNumComp; i++) {
    if (i == 0)
      charData[i] = charData[i%numDataPoint];
    intData[i] = i;
    floatData[i] = intData[i] + 0.1;
    doubleData[i] = intData[i] + 0.00000000001;
  }

    
  for (i = 0; i < numMathType; i++) {

    // create variables
    for (j = 0; j < numDataType; j++) {
      rc = fc_createVariable(mesh, name, &variables[j]);
      fail_unless(rc == FC_SUCCESS, "abort: failed to create variable");
      rc = fc_setVariableData(variables[j], numDataPoint, numComps[i],
			      assoc,  mathTypes[i], dataTypes[j],
			      data[j]);
    }
      
    // (Couldn't figure out good way to make this loop because of the
    // casting)

    // check getting char data as different types
    rc = fc_getVariableDataAsDataType(variables[0], FC_DT_CHAR, &temp_data);
    fail_unless(rc == FC_SUCCESS, "failed to get char data as char");
    for (k = 0; k < numDataPoint*numComps[i]; k++)
      fail_unless(((char*)(temp_data))[k] == charData[k],
		  "mismatch of char to char");
    free(temp_data);
    rc = fc_getVariableDataAsDataType(variables[0], FC_DT_INT, &temp_data);
    fail_unless(rc == FC_SUCCESS, "failed to get char data as int");
    for (k = 0; k < numDataPoint*numComps[i]; k++)
      fail_unless(((int*)(temp_data))[k] == (int)charData[k],
		  "mismatch of char to int");
    free(temp_data);
    rc = fc_getVariableDataAsDataType(variables[0], FC_DT_FLOAT, &temp_data);
    fail_unless(rc == FC_SUCCESS, "failed to get char data as float");
    for (k = 0; k < numDataPoint*numComps[i]; k++)
      fail_unless(((float*)(temp_data))[k] == (float)charData[k],
		  "mismatch of char to float");
    free(temp_data);
    rc = fc_getVariableDataAsDataType(variables[0], FC_DT_DOUBLE, &temp_data);
    fail_unless(rc == FC_SUCCESS, "failed to get char data as double");
    for (k = 0; k < numDataPoint*numComps[i]; k++)
      fail_unless(((double*)(temp_data))[k] == (double)charData[k],
		  "mismatch of char to double");
    free(temp_data);
   
    // check getting int data as different types
    rc = fc_getVariableDataAsDataType(variables[1], FC_DT_CHAR, &temp_data);
    fail_unless(rc == FC_SUCCESS, "failed to get int data as char");
    for (k = 0; k < numDataPoint*numComps[i]; k++)
      fail_unless(((char*)(temp_data))[k] == (char)intData[k],
		  "mismatch of int to char");
    free(temp_data);
    rc = fc_getVariableDataAsDataType(variables[1], FC_DT_INT, &temp_data);
    fail_unless(rc == FC_SUCCESS, "failed to get int data as int");
    for (k = 0; k < numDataPoint*numComps[i]; k++)
      fail_unless(((int*)(temp_data))[k] == intData[k],
		  "mismatch of int to int");
    free(temp_data);
    rc = fc_getVariableDataAsDataType(variables[1], FC_DT_FLOAT, &temp_data);
    fail_unless(rc == FC_SUCCESS, "failed to get int data as float");
    for (k = 0; k < numDataPoint*numComps[i]; k++)
      fail_unless(((float*)(temp_data))[k] == (float)intData[k],
		  "mismatch of int to float");
    free(temp_data);
    rc = fc_getVariableDataAsDataType(variables[1], FC_DT_DOUBLE, &temp_data);
    fail_unless(rc == FC_SUCCESS, "failed to get int data as double");
    for (k = 0; k < numDataPoint*numComps[i]; k++)
      fail_unless(((double*)(temp_data))[k] == (double)intData[k],
		  "mismatch of int to double");
    free(temp_data);
 
    // check getting float data as different types
    rc = fc_getVariableDataAsDataType(variables[2], FC_DT_CHAR, &temp_data);
    fail_unless(rc == FC_SUCCESS, "failed to get float data as char");
    for (k = 0; k < numDataPoint*numComps[i]; k++)
      fail_unless(((char*)(temp_data))[k] == (char)floatData[k],
		  "mismatch of float to char");
    free(temp_data);
    rc = fc_getVariableDataAsDataType(variables[2], FC_DT_INT, &temp_data);
    fail_unless(rc == FC_SUCCESS, "failed to get float data as int");
    for (k = 0; k < numDataPoint*numComps[i]; k++)
      fail_unless(((int*)(temp_data))[k] == (int)floatData[k],
		  "mismatch of float to int");
    free(temp_data);
    rc = fc_getVariableDataAsDataType(variables[2], FC_DT_FLOAT, &temp_data);
    fail_unless(rc == FC_SUCCESS, "failed to get float data as float");
    for (k = 0; k < numDataPoint*numComps[i]; k++)
      fail_unless(((float*)(temp_data))[k] == floatData[k],
		  "mismatch of float to float");
    free(temp_data);
    rc = fc_getVariableDataAsDataType(variables[2], FC_DT_DOUBLE, &temp_data);
    fail_unless(rc == FC_SUCCESS, "failed to get float data as double");
    for (k = 0; k < numDataPoint*numComps[i]; k++)
      fail_unless(((double*)(temp_data))[k] == (double)floatData[k],
		  "mismatch of float to double");
    free(temp_data);
 
    // check getting double data as different types
    rc = fc_getVariableDataAsDataType(variables[3], FC_DT_CHAR, &temp_data);
    fail_unless(rc == FC_SUCCESS, "failed to get double data as char");
    for (k = 0; k < numDataPoint*numComps[i]; k++)
      fail_unless(((char*)(temp_data))[k] == (char)doubleData[k],
		  "mismatch of double to char");
    free(temp_data);
    rc = fc_getVariableDataAsDataType(variables[3], FC_DT_INT, &temp_data);
    fail_unless(rc == FC_SUCCESS, "failed to get double data as int");
    for (k = 0; k < numDataPoint*numComps[i]; k++)
      fail_unless(((int*)(temp_data))[k] == (int)doubleData[k],
		  "mismatch of double to int");
    free(temp_data);
    rc = fc_getVariableDataAsDataType(variables[3], FC_DT_FLOAT, &temp_data);
    fail_unless(rc == FC_SUCCESS, "failed to get double data as float");
    for (k = 0; k < numDataPoint*numComps[i]; k++)
      fail_unless(((float*)(temp_data))[k] == (float)doubleData[k],
		  "mismatch of double to float");
    free(temp_data);
    rc = fc_getVariableDataAsDataType(variables[3], FC_DT_DOUBLE, &temp_data);
    fail_unless(rc == FC_SUCCESS, "failed to get double data as double");
    for (k = 0; k < numDataPoint*numComps[i]; k++)
      fail_unless(((double*)(temp_data))[k] == doubleData[k],
		  "mismatch of double to double");
    free(temp_data);
    
    // ---- check errors
    
    // fc_getVariableDataAsDataType() -- bad args
    rc = fc_getVariableDataAsDataType(badVariable, FC_DT_CHAR, &temp_data);
    fail_unless(rc != FC_SUCCESS, "should fail with bad variable");
    fail_unless(temp_data == NULL, "fail should return NULL");
    rc = fc_getVariableDataAsDataType(variables[0], FC_DT_UNKNOWN, &temp_data);
    fail_unless(rc != FC_SUCCESS, "should fail with FC_DT_UNKNOWN");
    fail_unless(temp_data == NULL, "fail should return NULL");
    rc = fc_getVariableDataAsDataType(variables[0], -99, &temp_data);
    fail_unless(rc != FC_SUCCESS, "should fail with bad datatype");
    fail_unless(temp_data == NULL, "fail should return NULL");
  }
}
END_TEST

// copy
START_TEST(copy_test)
{
  FC_ReturnCode rc;
  int i, j, k, m;
  char name[20] = "blue berry", copy_name[20] = "banana nut", *temp_name;
  FC_Dataset dataset, *returnDatasets;
  FC_Mesh mesh1, mesh2, weirdMesh, *returnMeshes;
  int numReturnDatasets, numReturnMeshes;
  int maxNumEntity = 20; // for all the meshes and all assoc types
  int numVariable1, numVariable2, temp_numVariable;
  FC_Variable variable, copy_variable, badVariable = { 999, 999 };
  FC_Variable glbVar;
  int numAssoc = 5, numDataPoint;  // numDataPoint determined later
  FC_AssociationType assocs[5] = { FC_AT_VERTEX, FC_AT_EDGE, FC_AT_FACE,
				   FC_AT_ELEMENT, FC_AT_WHOLE_MESH };
  int numMathType = 4, numComps[4] = { 1, 2, 3, 4 }, maxNumComp = 4;
  FC_MathType mathTypes[4] = { FC_MT_SCALAR, FC_MT_VECTOR, 
			       FC_MT_SYMTENSOR, FC_MT_TENSOR };
  int numDataType = 4;
  FC_DataType dataTypes[4] = { FC_DT_CHAR, FC_DT_INT, 
			       FC_DT_FLOAT, FC_DT_DOUBLE };
  char charData[20*4] = { 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j',
			'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't' };
  int intData[20*4];
  float floatData[20*4];
  double doubleData[20*4];
  int sizes[4] = {sizeof(char), sizeof(int), sizeof(float), sizeof(double) };
  void *data[4] = { charData, intData, floatData, doubleData };
  int temp_numDataPoint, temp_numComp;
  FC_AssociationType temp_assoc;
  FC_MathType temp_mathType;
  FC_DataType temp_dataType;
  void *temp_data;
  
  // get dataset and mesh
  rc = fc_getDatasetByName(dataset_name, &numReturnDatasets,
			   &returnDatasets);
  fail_unless(rc == FC_SUCCESS, "aborted - cant get dataset by name");
  fail_unless(numReturnDatasets == 1,
	      "failed to find unique dataset by name");
  dataset = returnDatasets[0];
  free(returnDatasets);
  rc = fc_getMeshByName(dataset, "hex mesh", &numReturnMeshes,
			&returnMeshes);
  fail_unless(rc == FC_SUCCESS, "aborted - cant get mesh by name");
  fail_unless(numReturnMeshes == 1, "failed to find unique mesh by name");
  mesh1 = returnMeshes[0];
  free(returnMeshes);

  // make a second mesh to copy to
  rc = fc_copyMesh(mesh1, dataset, "copy mesh", &mesh2);
  fail_unless(rc == FC_SUCCESS, "abort: failed to copy mesh");

  // make a third mesh we shouldn't be able to copy to - just 1 element
  {
    int topodim, dim, numVertex, numElement, *conns;
    double* coords;
    FC_ElementType elemType;
    fc_getMeshInfo(mesh1, &topodim, &dim, &numVertex, &numElement, &elemType);
    fc_getMeshCoordsPtr(mesh1, &coords);
    fc_getMeshElementConnsPtr(mesh1, &conns);
    rc = fc_createMesh(dataset, "different mesh", &weirdMesh);
    fail_unless(rc == FC_SUCCESS, "abort: failed to create mesh");
    rc = fc_setMeshCoords(weirdMesh, dim, fc_getElementTypeNumVertex(elemType),
			  coords);
    fail_unless(rc == FC_SUCCESS, "abort: failed to set coords");
    rc = fc_setMeshElementConns(weirdMesh, elemType, 1, conns);
    fail_unless(rc == FC_SUCCESS, "abort: failed to set conns");
  }

  // setup some data
  for (i = 0; i < maxNumEntity*maxNumComp; i++) {
    if (i == 0)
      charData[i] = charData[i%maxNumEntity];
    intData[i] = i;
    floatData[i] = intData[i] + 0.1;
    doubleData[i] = intData[i] + 0.00000000001;
  }

  // Make a global variable (double)
  rc = fc_createGlobalVariable(dataset, name, &glbVar);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create glb variable");
  rc = fc_setVariableData(glbVar, 1, 1, FC_AT_WHOLE_DATASET,  FC_MT_SCALAR,
                          FC_DT_DOUBLE, doubleData);
  fail_unless(rc == FC_SUCCESS, "abort: failed to set glb var data");

  numVariable1 = 0;
  numVariable2 = 0;
  // association types, math types & data types
  for (j = 0; j < numAssoc; j++) {
    
    // setup numDataPoint
    fc_getMeshNumEntity(mesh1, assocs[j], &numDataPoint);
    
    for (k = 0; k < numMathType; k++) {
      for (m = 0; m < numDataType; m++) {
	
	//printf("assoc %s, mathtype %s, datatype %s\n",
	//     fc_getAssociationTypeText(assocs[j]),
	//     fc_getMathTypeText(mathTypes[k]),
	//     fc_getDataTypeText(dataTypes[m]));
	
	// create a variable
	rc = fc_createVariable(mesh1, name, &variable);
	fail_unless(rc == FC_SUCCESS, "abort: failed to create variable");
	rc = fc_setVariableData(variable, numDataPoint, numComps[k],
				assocs[j],  mathTypes[k], dataTypes[m],
				data[m]);
	numVariable1++;

	// copy to same mesh & check it
	rc = fc_copyVariable(variable, mesh1, copy_name, &copy_variable);
	fail_unless(rc == FC_SUCCESS, "failed to copy to same mesh");
	numVariable1++;
	rc = fc_getVariableName(copy_variable, &temp_name);
	fail_unless(rc == FC_SUCCESS, "failed to get name of copy");
	fail_unless(!strcmp(temp_name, copy_name), "mismatch of name");
	free(temp_name);
	rc = fc_getVariableInfo(copy_variable, &temp_numDataPoint, 
				&temp_numComp, &temp_assoc, &temp_mathType, 
				&temp_dataType);
	fail_unless(rc == FC_SUCCESS, "failed to get variable info");
	fail_unless(temp_numDataPoint == numDataPoint, 
		    "mismatch of numDataPoint");
	fail_unless(temp_numComp == numComps[k], "mismatch of numComp");
	fail_unless(temp_assoc == assocs[j], "mismatch of assoc");
	fail_unless(temp_mathType == mathTypes[k], "mismatch of mathType");
	fail_unless(temp_dataType == dataTypes[m], "mismatch of datatype");
	rc = fc_getVariableDataPtr(copy_variable, &temp_data);
	fail_unless(rc == FC_SUCCESS, "failed to get data");
	fail_unless(!memcmp(temp_data, data[m], 
			    numDataPoint*numComps[k]*sizes[m]),
		    "mismatch of data");

 	// copy to different mesh & check it
	rc = fc_copyVariable(variable, mesh2, copy_name, &copy_variable);
	fail_unless(rc == FC_SUCCESS, "failed to copy to different mesh");
	numVariable2++;
	rc = fc_getVariableName(copy_variable, &temp_name);
	fail_unless(rc == FC_SUCCESS, "failed to get name of copy");
	fail_unless(!strcmp(temp_name, copy_name), "mismatch of name");
	free(temp_name);
	rc = fc_getVariableInfo(copy_variable, &temp_numDataPoint, 
				&temp_numComp, &temp_assoc, &temp_mathType, 
				&temp_dataType);
	fail_unless(rc == FC_SUCCESS, "failed to get variable info");
	fail_unless(temp_numDataPoint == numDataPoint, 
		    "mismatch of numDataPoint");
	fail_unless(temp_numComp == numComps[k], "mismatch of numComp");
	fail_unless(temp_assoc == assocs[j], "mismatch of assoc");
	fail_unless(temp_mathType == mathTypes[k], "mismatch of mathType");
	fail_unless(temp_dataType == dataTypes[m], "mismatch of datatype");
	rc = fc_getVariableDataPtr(copy_variable, &temp_data);
	fail_unless(rc == FC_SUCCESS, "failed to get data");
	fail_unless(!memcmp(temp_data, data[m], 
			    numDataPoint*numComps[k]*sizes[m]),
		    "mismatch of data");

	// fc_copyVariable() to a very different mesh should fail for non WHOLE
	if (assocs[j] != FC_AT_WHOLE_MESH) {
	  rc = fc_copyVariable(variable, weirdMesh, copy_name, &copy_variable);
	  fail_unless(rc != FC_SUCCESS, 
		      "should fail to copy to very different mesh");
	  fail_unless(FC_HANDLE_EQUIV(copy_variable, FC_NULL_VARIABLE),
		      "fail should return NULL handle");
	}

	// --- special cases
	
	// copy with NULL name will use original's name
	rc = fc_copyVariable(variable, mesh1, NULL, &copy_variable);
	fail_unless(rc == FC_SUCCESS, "should not fail");
	rc = fc_getVariableName(copy_variable, &temp_name);
	fail_unless(rc == FC_SUCCESS, "failed to get name of copy");
	fail_unless(!strcmp(temp_name, name), "mismatch of name");
	free(temp_name);
	rc = fc_deleteVariable(copy_variable);
	fail_unless(rc == FC_SUCCESS, "failed to delete copied variable");
    
	// ---- check errors
	
	// fc_copyVariable() -- bad args
	rc = fc_copyVariable(badVariable, mesh1, copy_name, &copy_variable);
	fail_unless(rc != FC_SUCCESS, "should fail to copy bad variable");
	fail_unless(FC_HANDLE_EQUIV(copy_variable, FC_NULL_VARIABLE),
		    "fail should return NULL");
	rc = fc_copyVariable(variable, FC_NULL_MESH, copy_name, 
			     &copy_variable);
	fail_unless(rc != FC_SUCCESS, "should fail to copy to bad mesh");
	fail_unless(FC_HANDLE_EQUIV(copy_variable, FC_NULL_VARIABLE),
		    "fail should return NULL");
	rc = fc_copyVariable(variable, mesh1, copy_name, NULL);
	fail_unless(rc != FC_SUCCESS, "should fail to copy if NULL handle");
      }
    }
  }
  
  // It's o.k. to copy a global var onto a mesh
  rc = fc_copyVariable(glbVar, mesh1, copy_name, &copy_variable);
  fail_unless(rc == FC_SUCCESS, "failed to copy glb var to a mesh");
  numVariable1++;
  rc = fc_getVariableName(copy_variable, &temp_name);
  fail_unless(rc == FC_SUCCESS, "failed to get name of copy");
  fail_unless(!strcmp(temp_name, copy_name), "mismatch of name");
  free(temp_name);
  rc = fc_getVariableInfo(copy_variable, &temp_numDataPoint, 
                          &temp_numComp, &temp_assoc, &temp_mathType, 
                          &temp_dataType);
  fail_unless(rc == FC_SUCCESS, "failed to get variable info");
  fail_unless(temp_numDataPoint == 1, "mismatch of numDataPoint");
  fail_unless(temp_numComp == 1, "mismatch of numComp");
  fail_unless(temp_assoc == FC_AT_WHOLE_MESH, "mismatch of assoc");
  fail_unless(temp_mathType == FC_MT_SCALAR, "mismatch of mathType");
  fail_unless(temp_dataType == FC_DT_DOUBLE, "mismatch of datatype");
  rc = fc_getVariableDataPtr(copy_variable, &temp_data);
  fail_unless(rc == FC_SUCCESS, "failed to get data");
  fail_unless(!memcmp(temp_data, doubleData, sizeof(double)),
              "mismatch of data");

  // a little more checking
  rc = fc_deleteMesh(weirdMesh);
  fail_unless(rc == FC_SUCCESS, "failed to delete weird mesh");  
  rc = fc_getNumVariable(mesh2, &temp_numVariable);
  fail_unless(rc == FC_SUCCESS, "shouldn't fail to get num variables");
  fail_unless(temp_numVariable == numVariable2, "mismatch of numVariable");
  rc = fc_deleteMesh(mesh2);
  fail_unless(rc == FC_SUCCESS, "failed to delete second mesh");
  rc = fc_getNumVariable(mesh1, &temp_numVariable);
  fail_unless(rc == FC_SUCCESS, "shouldn't fail to get num variable");
  fail_unless(temp_numVariable == numVariable1, "mismatch of numVariable");
}
END_TEST

// release vars
START_TEST(release_var)
{
  FC_ReturnCode rc;
  FC_Dataset dataset;
  FC_Mesh mesh, *returnMeshes;
  FC_Variable *returnVariables;
  int numReturnMeshes, numReturnVariables;
  FC_Variable var, var2, badVar = { 999, 999 };
  _FC_VarSlot* varSlot, *varSlot2;
  void *data;

  // setup
  rc = fc_loadDataset("../data/gen_hex.ex2", &dataset);
  fail_unless(rc == FC_SUCCESS, "abort: failed to load dataset");
  rc = fc_getMeshByName(dataset, "hex mesh", &numReturnMeshes,&returnMeshes);
  fail_unless(rc == FC_SUCCESS, "aborted - cant get mesh by name");
  fail_unless(numReturnMeshes == 1, "failed to find unique mesh by name");
  mesh = returnMeshes[0];
  free(returnMeshes);

  rc = fc_getVariableByName(mesh, "temperature", &numReturnVariables,&returnVariables);
  fail_unless(rc == FC_SUCCESS, "aborted - cant get variables by name");
  fail_unless(numReturnVariables == 1, "failed to find unique variable by name");
  var = returnVariables[0];
  free(returnVariables);

  varSlot = _fc_getVarSlot(var);
  fail_unless(varSlot != NULL, "abort: failed to get varslot");

  // test that lazily loaded data is not there yet
  fail_unless(varSlot->data == NULL, "abort: should start with no data");

  // make copy (uncommitted) variable
  rc = fc_copyVariable(var, mesh, "new_var", &var2);
  fail_unless(rc == FC_SUCCESS, "failed to copy var");
  varSlot2 = _fc_getVarSlot(var2);
  fail_unless(varSlot2 != NULL, "abort: failed to get varslot2");

  // --- data

  // load & unload variable data
  rc = fc_getVariableDataPtr(var, &data); // force loading
  fail_unless(rc == FC_SUCCESS, "failed to get data");
  fail_unless(varSlot->data != NULL, "abort: data should have loaded");
  fail_unless(data == varSlot->data, "mismatch of data pointers");
  rc = fc_releaseVariable(var);
  fail_unless(rc == FC_SUCCESS, "failed to release data");
  fail_unless(varSlot->data == NULL, "data should now be null");

  // uncommitted var's shouldn't get deleted
  fail_unless(varSlot2->data != NULL, "should start with data");
  fc_releaseVariable(var); // make sure that copy doesn't just copy pointers
  rc = fc_releaseVariable(var2);
  fail_unless(rc == FC_SUCCESS, "failed to release data");
  fail_unless(varSlot2->data != NULL, "should not delete temp var data");

  // --- special cases
  
  // releasing a null handle does not cause an error
  rc = fc_releaseVariable(FC_NULL_VARIABLE);
  fail_unless(rc == FC_SUCCESS, "NULL handle should not fail");
  
  // --- errors
  
  // fc_releaseVariable() --- bad args
  rc = fc_releaseVariable(badVar);
  fail_unless(rc != FC_SUCCESS, "bad variable should fail");
 
  // cleanup
  fc_deleteDataset(dataset);
}
END_TEST 

// test error input for fc_printVariable()
START_TEST(print)
{
  FC_ReturnCode rc;
  FC_Variable badVariable = { 999, 999 };

  // fc_printVariable()
  rc = fc_printVariable(badVariable, "Bad Variable!", 1);
  fail_unless(rc != FC_SUCCESS, "should fail if bad variable");
}
END_TEST 

// --- seq variable interface

// test creating and deleting seq variables, also test querying mesh for vars
// for variables before, inbetween and after.
START_TEST(seq_create_get_delete)
{
  FC_ReturnCode rc;
  int i, j, k, m;
  FC_Dataset dataset, *returnDatasets;
  FC_Sequence *sequences;
  FC_Mesh *meshes, mesh;
  int numReturnDatasets,numMesh, numSequence;
  int numSeqVar = 10, temp_numSeqVar, numSteps[2], temp_numStep, *temp_numSteps;
  FC_Variable *seqVars[2*10], **temp_seqVars, *temp_seqVar;
  FC_Variable *badSeqVar, badVariable = { 999, 999 };
  char names[2][10][100] = { { "one", "two", "three", "four", "five", 
			       "six", "seven", "eight", "nine", "ten" },
			     { "un", "deux", "trois", "quatre", "cinq",
			       "six (french)", "sept", "huit", "neuf", "dix"}};
  char newName[100] = { "sparkly new" };

  _FC_VarSlot *varSlot, *varSlot2;

  // get dataset, sequences and meshes
  rc = fc_getDatasetByName(dataset_name, &numReturnDatasets,
			   &returnDatasets);
  fail_unless(rc == FC_SUCCESS, "aborted - cant get dataset by name");
  fail_unless(numReturnDatasets == 1,
	      "failed to find unique dataset by name");
  dataset = returnDatasets[0];
  free(returnDatasets);
  rc = fc_getSequences(dataset, &numSequence, &sequences);
  fail_unless(rc == FC_SUCCESS, "abort: failed to get sequences");
  fail_unless(numSequence == 2, "abort: expected 2 sequences");
  rc = fc_getMeshes(dataset, &numMesh, &meshes);
  fail_unless(rc == FC_SUCCESS, "abort: failed to get meshes");
  fail_unless(numMesh == 8, "abort: expected 8 meshes");
  for (i = 0; i < numSequence; i++) 
    fc_getSequenceNumStep(sequences[i], &numSteps[i]);

  // loop over all mesh types 
  for (i = 0; i < numMesh; i++) { 
    mesh = meshes[i];

    // create some seq variables for each sequence
    for (j = 0; j < numSequence; j++) {
      for (k = 0; k < numSeqVar; k++) {
	int idx = j*numSeqVar + k;
	rc = fc_createSeqVariable(mesh, sequences[j], names[j][k], 
				  &temp_numStep, &seqVars[idx]);
	fail_unless(rc == FC_SUCCESS, "failed to create seqVar");
	fail_unless(temp_numStep == numSteps[j], "mismatch in numStep");
	for (m = 0; m < numSteps[j]; m++)
	  fail_unless(!FC_HANDLE_EQUIV(seqVars[idx][m], FC_NULL_VARIABLE),
		      "created variable should not = FC_NULL_VARIABLE");
      }
    }

    // get seq variables (should be numSeqVar*numSequence);
    // (will come out in same order they were written)
    rc = fc_getSeqVariables(mesh, &temp_numSeqVar, &temp_numSteps, 
			    &temp_seqVars);
    fail_unless(rc == FC_SUCCESS, "failed to get seq variables");
    fail_unless(temp_numSeqVar == numSeqVar*numSequence, 
		"mismatch of numSeqVariable");
    for (j = 0; j < numSequence; j++) {
      for (k = 0; k < numSeqVar; k++) {
	int idx = j*numSeqVar + k;
	varSlot = _fc_getVarSlot(temp_seqVars[idx][0]);
	fail_unless(FC_HANDLE_EQUIV(varSlot->sequence, sequences[j]),
		    "mismatch of sequence handle");
	fail_unless(temp_numSteps[idx] == numSteps[j],
		    "mismatch of numStep");
	for (m = 0; m < temp_numSteps[idx]; m++)
	  fail_unless(FC_HANDLE_EQUIV(temp_seqVars[idx][m], 
				      seqVars[idx][m]), 
		      "mismatch of variable handles");
	free(temp_seqVars[idx]);
      }
    }
    free(temp_numSteps);
    free(temp_seqVars);
    rc = fc_getNumSeqVariable(mesh, &temp_numSeqVar);
    fail_unless(rc == FC_SUCCESS, "failed to get num seq variables");
    fail_unless(temp_numSeqVar == numSeqVar*numSequence, 
		"mismatch of numSeqVariable");

    // test fc_getSeqVariableByName()
    for (j = 0; j < numSequence; j++) {
      for (k = 0; k < numSeqVar; k++) {
	int idx = j*numSeqVar + k;
	rc = fc_getSeqVariableByName(mesh, names[j][k], &temp_numSeqVar,
				     &temp_numSteps, 
				     &temp_seqVars);
	fail_unless(rc == FC_SUCCESS, "failed to get variable by name");
	fail_unless(temp_numSeqVar == 1, "wrong number of matching vars");
	varSlot = _fc_getVarSlot(temp_seqVars[0][0]);
	fail_unless(FC_HANDLE_EQUIV(varSlot->sequence, sequences[j]),
		    "mismatch of sequence handle");
	fail_unless(temp_numSteps[0] == numSteps[j], "mismatch of numStep");
	for (m = 0; m < temp_numSteps[0]; m++)
	  fail_unless(FC_HANDLE_EQUIV(temp_seqVars[0][m], seqVars[idx][m]),
		      "mismatch of variable handles");
	free(temp_seqVars[0]);
	free(temp_seqVars);
	free(temp_numSteps);
      }
    }

    //temporarily create a new seq var to check array return
    //make one the same as the 0th one, but on the first sequence
    rc = fc_createSeqVariable(mesh, sequences[1], names[0][0], 
			      &temp_numStep, &temp_seqVar);
    fail_unless(rc == FC_SUCCESS, "failed to create seqVar");
    fail_unless(temp_numStep == numSteps[1], "mismatch in numStep");
    for (m = 0; m < numSteps[1]; m++)
      fail_unless(!FC_HANDLE_EQUIV(temp_seqVar[m], FC_NULL_VARIABLE),
		  "created variable should not = FC_NULL_VARIABLE");

    rc = fc_getSeqVariableByName(mesh, names[0][0],
				 &temp_numSeqVar,
				 &temp_numSteps,
				 &temp_seqVars);
    fail_unless(rc == FC_SUCCESS, "failed to get variable by name");
    fail_unless(temp_numSeqVar == 2, "wrong number of matching seq vars");
    varSlot = _fc_getVarSlot(temp_seqVars[0][0]);
    varSlot2 = _fc_getVarSlot(temp_seqVars[1][0]);
    fail_unless(((FC_HANDLE_EQUIV(varSlot->sequence, sequences[0]) &&
		  FC_HANDLE_EQUIV(varSlot2->sequence, sequences[1])) ||
		 (FC_HANDLE_EQUIV(varSlot->sequence, sequences[1]) &&
		  FC_HANDLE_EQUIV(varSlot2->sequence, sequences[0]))),
		"mismatch of sequence handle");
    j = !(FC_HANDLE_EQUIV(varSlot->sequence, sequences[0]));
    fail_unless(temp_numSteps[0] == numSteps[j] &&
		temp_numSteps[1] == numSteps[!j], "mismatch of numStep");
    for (m = 0; m < temp_numSteps[j]; m++){
      fail_unless(FC_HANDLE_EQUIV(temp_seqVars[j][m], seqVars[0][m]),
		  "mismatch of variable handles");
    }
    for (m = 0; m < temp_numSteps[!j]; m++){
      fail_unless(FC_HANDLE_EQUIV(temp_seqVars[!j][m], temp_seqVar[m]),
		  "mismatch of variable handles");
    }

    fc_deleteSeqVariable(temp_numStep,temp_seqVar);
    free(temp_seqVar);
    
    for (m = 0; m < temp_numSeqVar; m++){
      free(temp_seqVars[m]);
    }
    free(temp_seqVars);
    free(temp_numSteps);

    // change the name of the first sequence
    // also test negative example of fc_getVariableByName
    rc = fc_changeSeqVariableName(numSteps[0], seqVars[0], newName);
    fail_unless(rc == FC_SUCCESS, "failed to change seq variable name");
    rc = fc_getSeqVariableByName(mesh, names[0][0], &temp_numSeqVar,
				 &temp_numSteps,
				 &temp_seqVars);
    fail_unless(rc == FC_SUCCESS, "old name should work");
    fail_unless(temp_numSeqVar == 0 && temp_numSteps == NULL &&
		temp_seqVars == NULL, "old name should find no match");
    rc = fc_getSeqVariableByName(mesh, newName, &temp_numSeqVar,
				 &temp_numSteps, &temp_seqVars);
    fail_unless(rc == FC_SUCCESS, "new name should work for find");
    fail_unless(temp_numSeqVar ==1 ,"wrong number of matching vars");
    fail_unless(temp_numSteps[0] == numSteps[0], "mismatch of numStep");
    for (j = 0; j < numSteps[0]; j++)
      fail_unless(FC_HANDLE_EQUIV(temp_seqVars[0][j], seqVars[0][j]),
		"mismatch of variable handle");
    free(temp_seqVars[0]);
    free(temp_seqVars);
    free(temp_numSteps);


    // delete half of the seq variables (alternate)
    for (j = 0; j < numSequence; j++) {
      for (k = 0; k < numSeqVar; k+=2) {
	int idx = j*numSeqVar + k;
	rc = fc_deleteSeqVariable(numSteps[j], seqVars[idx]);
	fail_unless(rc == FC_SUCCESS, "failed to delete seq variable");
	free(seqVars[idx]);
      }
    }

    // get seq variables (should be numSeqVar*numSequence/2);
    fc_getSeqVariables(mesh, &temp_numSeqVar, &temp_numSteps, &temp_seqVars);
    fail_unless(rc == FC_SUCCESS, "failed to get seq variables");
    fail_unless(temp_numSeqVar == numSeqVar*numSequence/2, 
		"mismatch of numSeqVariable");
    for (j = 0; j < numSequence; j++) {
      for (k = 0; k < numSeqVar/2; k++) {
	int idx = j*numSeqVar/2 + k; // into temp array
	int idx_orig = j*numSeqVar + 2*k + 1; // into original array
	varSlot = _fc_getVarSlot(temp_seqVars[idx][0]);
	fail_unless(FC_HANDLE_EQUIV(varSlot->sequence, sequences[j]),
		    "mismatch of sequence handle");
	fail_unless(temp_numSteps[idx] == numSteps[j],
		    "mismatch of numStep");
	for (m = 0; m < temp_numSteps[idx]; m++)
	  fail_unless(FC_HANDLE_EQUIV(temp_seqVars[idx][m],
				      seqVars[idx_orig][m]), 
		      "mismatch of variable handles");
	free(temp_seqVars[idx]);
      }
    }
    free(temp_numSteps);
    free(temp_seqVars);
    
    // delete remaining seq variables on first sequence
    for (k = 1; k < numSeqVar; k+=2) {
      rc = fc_deleteSeqVariable(numSteps[0], seqVars[k]);
      fail_unless(rc == FC_SUCCESS, "failed to delete variable");
      free(seqVars[k]);
    }
    
    // get seq variables (should be numSeqVar*(numSequence-1)/2);
    fc_getSeqVariables(mesh, &temp_numSeqVar, &temp_numSteps, &temp_seqVars);
    fail_unless(rc == FC_SUCCESS, "failed to get seq variables");
    fail_unless(temp_numSeqVar == numSeqVar*(numSequence-1)/2, 
		"mismatch of numSeqVariable");
    for (j = 1; j < numSequence; j++) {
      for (k = 0; k < numSeqVar/2; k++) {
	int idx = (j-1)*numSeqVar/2 + k; // into temp array
	int idx_orig = j*numSeqVar + 2*k + 1; // into original array
	varSlot = _fc_getVarSlot(temp_seqVars[idx][0]);
	fail_unless(FC_HANDLE_EQUIV(varSlot->sequence, sequences[j]),
		    "mismatch of sequence handle");
	fail_unless(temp_numSteps[idx] == numSteps[j],
		    "mismatch of numStep");
	for (m = 0; m < temp_numSteps[idx]; m++)
	  fail_unless(FC_HANDLE_EQUIV(temp_seqVars[idx][m],
				      seqVars[idx_orig][m]), 
		      "mismatch of variable handles");
	free(temp_seqVars[idx]);
      }
    }
    free(temp_numSteps);
    free(temp_seqVars);
   
    // delete remaining variables
    for (j = 1; j < numSequence; j++) {
      for (k = 1; k < numSeqVar; k+=2) { 
	rc = fc_deleteSeqVariable(numSteps[j], seqVars[j*numSeqVar+k]);
	fail_unless(rc == FC_SUCCESS, "failed to delete variable");
	free(seqVars[j*numSeqVar+k]);
      }
    }

    // get seq variables (should be none);
    rc = fc_getSeqVariables(mesh, &temp_numSeqVar, &temp_numSteps, 
			    &temp_seqVars);
    fail_unless(rc == FC_SUCCESS, "failed to get seq variables");
    fail_unless(temp_numSeqVar == 0 && temp_numSteps == NULL && 
		temp_seqVars == NULL, 
		"should return 0 & nulls if all seqVars delete");
    
    // create 1 seq var for further testing
    rc = fc_createSeqVariable(mesh, sequences[0], names[0][0],
			      &temp_numStep, &seqVars[0]);
    fail_unless(rc == FC_SUCCESS, "abort: failed to create seq variable");
 
    // ---- test special cases

    // rc_deleteSeqVariable() - not an error to delete (0, NULL)
    rc = fc_deleteSeqVariable(0, NULL);
    fail_unless(rc == FC_SUCCESS, "should not error to delete NULL seqVar");
    
    // ---- test error conditions
    
    // construct a bad seqVar (make one of the vars bad)
    badSeqVar = (FC_Variable*)malloc(numSteps[0]*sizeof(FC_Variable));
    for (j = 0; j < numSteps[0]; j++)
      badSeqVar[j] = seqVars[0][j];
    badSeqVar[numSteps[0]/2] = badVariable;

    // bad args to fc_createSeqVariable()
    rc = fc_createSeqVariable(FC_NULL_MESH, sequences[0], names[0][0], 
			      &temp_numStep, &temp_seqVar);
    fail_unless(rc != FC_SUCCESS, "should fail with bad mesh");
    fail_unless(temp_numStep == -1 && temp_seqVar == NULL, 
		"fail should return NULL seqVar");
    rc = fc_createSeqVariable(meshes[0], FC_NULL_SEQUENCE, names[0][0], 
			      &temp_numStep, &temp_seqVar);
    fail_unless(rc != FC_SUCCESS, "should fail with bad sequence");
    fail_unless(temp_numStep == -1 && temp_seqVar == NULL, 
		"fail should return NULL seqVar");
    rc = fc_createSeqVariable(meshes[0], sequences[0], NULL, 
			      &temp_numStep, &temp_seqVar);
    fail_unless(rc != FC_SUCCESS, "should fail with null name");
    fail_unless(temp_numStep == -1 && temp_seqVar == NULL, 
		"fail should return NULL seqVar");
    rc = fc_createSeqVariable(meshes[0], sequences[0], names[0][0], 
			      NULL, &temp_seqVar);
    fail_unless(rc != FC_SUCCESS, "should fail with null numStep");
    fail_unless(temp_seqVar == NULL, "fail should return NULL seqVar");
    rc = fc_createSeqVariable(meshes[0], sequences[0], names[0][0], 
			      &temp_numStep, NULL);
    fail_unless(rc != FC_SUCCESS, "should fail with null seqVar");
    fail_unless(temp_numStep == -1, "fail should return NULL");
    
    // bad args to fc_changeSeqVariableName()
    rc = fc_changeSeqVariableName(numSteps[0], seqVars[0], NULL);
    fail_unless(rc != FC_SUCCESS, "should fail if new name is NULL");
    rc = fc_getSeqVariableByName(mesh, names[0][0], &temp_numSeqVar,
				 &temp_numSteps, 
				 &temp_seqVars);
    fail_unless(rc == FC_SUCCESS && temp_numSeqVar == 1 &&
		temp_numSteps[0] == numSteps[0],
		"should not change name of sequence");
    for (j = 0; j < numSteps[0]; j++)
      fail_unless(FC_HANDLE_EQUIV(seqVars[0][j], temp_seqVars[0][j]),
		  "should not change name of sequence");
    free(temp_seqVars[0]);
    free(temp_seqVars);
    free(temp_numSteps);


    // bad args to fc_getSeqVariables()
    temp_numSteps = (int*)1;
    temp_seqVars = (FC_Variable**)1;
    rc = fc_getSeqVariables(FC_NULL_MESH, &temp_numSeqVar, &temp_numSteps, 
			    &temp_seqVars);
    fail_unless(rc != FC_SUCCESS, 
		"should have failed to get variables with null mesh");
    fail_unless(temp_numSeqVar == -1 && temp_numSteps == NULL && 
		temp_seqVars == NULL, "fail should return nulls"); 
    temp_numSteps = (int*)1;
    temp_seqVars = (FC_Variable**)1;
    rc = fc_getSeqVariables(mesh, NULL, &temp_numSteps, &temp_seqVars);
    fail_unless(rc != FC_SUCCESS, 
		"should have failed with NULL numSeqVar");
    fail_unless(temp_numSteps == NULL && temp_seqVars == NULL, 
		"fail should return nulls"); 
    temp_seqVars = (FC_Variable**)1;
    rc = fc_getSeqVariables(mesh, &temp_numSeqVar, NULL, &temp_seqVars);
    fail_unless(rc != FC_SUCCESS, 
		"should have failed with NULL numStepPerVar");
    fail_unless(temp_numSeqVar == -1 && temp_seqVars == NULL, 
		"fail should return nulls"); 
    temp_numSteps = (int*)1;
    rc = fc_getSeqVariables(mesh, &temp_numSeqVar, &temp_numSteps, NULL);
    fail_unless(rc != FC_SUCCESS, 
		"should have failed with NULL seqVars");
    fail_unless(temp_numSeqVar == -1 && temp_numSteps == NULL, 
		"fail should return nulls"); 

    // bad args to fc_getNumSeqVariable()
    rc = fc_getNumSeqVariable(FC_NULL_MESH, &temp_numSeqVar);
    fail_unless(rc != FC_SUCCESS, 
		"should have failed to get num seq variables with null mesh");
    fail_unless(temp_numSeqVar == -1, "fail should return nulls"); 
    rc = fc_getNumSeqVariable(mesh, NULL);
    fail_unless(rc != FC_SUCCESS, 
		"should have failed with NULL numSeqVar");

    // bad args to fc_getSeqVariableByName()
  // bad args to fc_getGlobalSeqVarByName()
    rc = fc_getSeqVariableByName(FC_NULL_MESH, names[0][0],
				       &temp_numSeqVar,
				       &temp_numSteps,
				       &temp_seqVars);
    fail_unless(rc != FC_SUCCESS, "should fail with NULL MESH");
    fail_unless(temp_numSeqVar == -1 && temp_numSteps == NULL &&
		temp_seqVars == NULL, "should return -1 and NULL when fail");

    rc = fc_getSeqVariableByName(mesh, NULL,
				       &temp_numSeqVar,
				       &temp_numSteps,
				       &temp_seqVars);
    fail_unless(rc != FC_SUCCESS, "should fail with NULL name");
    fail_unless(temp_numSeqVar == -1 && temp_numSteps == NULL &&
		temp_seqVars == NULL, "should return -1 and NULL when fail");

    rc = fc_getSeqVariableByName(mesh, names[0][0],
				       NULL,
				     &temp_numSteps,
				     &temp_seqVars);
    fail_unless(rc != FC_SUCCESS, "should fail with NULL arg");
    fail_unless(temp_numSteps == NULL &&
		temp_seqVars == NULL, "should return -1 and NULL when fail");

    rc = fc_getSeqVariableByName(mesh, names[0][0],
				 &temp_numSeqVar,
				 NULL,
				 &temp_seqVars);
    fail_unless(rc != FC_SUCCESS, "should fail with NULL arg");
    fail_unless(temp_numSeqVar == -1 && temp_seqVars == NULL,
		"should return -1 and NULL when fail");

    rc = fc_getSeqVariableByName(mesh, names[0][0],
				 &temp_numSeqVar,
				 &temp_numSteps,
				 NULL);
    fail_unless(rc != FC_SUCCESS, "should fail with NULL name");
    fail_unless(temp_numSeqVar == -1 && temp_numSteps == NULL,
		"should return -1 and NULL when fail");

    // bad args to fc_deleteVariable()
    rc = fc_deleteSeqVariable(numSteps[0]-1, seqVars[0]);
    fail_unless(rc != FC_SUCCESS, 
		"should error if numStep doesn't correspond to seqVar");
    rc = fc_deleteSeqVariable(numSteps[0], NULL);
    fail_unless(rc != FC_SUCCESS, "should error if seqVar = NULL");
    rc = fc_deleteSeqVariable(numSteps[0], badSeqVar);
    fail_unless(rc != FC_SUCCESS, "should error if badSeqVar");

    // --- done

    // delete last sequence variable
    rc = fc_deleteSeqVariable(numSteps[0], seqVars[0]);
    fail_unless(rc == FC_SUCCESS, "failed to delete last seq var");

    // cleanup
    free(seqVars[0]);
    free(badSeqVar);
  }
  free(sequences);
  free(meshes);
}
END_TEST

// test creating and deleting seq variables, also test querying mesh for vars
// for variables before, inbetween and after.
START_TEST(seq_metadata_query)
{
  FC_ReturnCode rc;
  int i, j, k;
  FC_Dataset dataset, temp_dataset, badDataset = { 999, 999 };
  FC_Dataset *returnDatasets;
  int numReturnDatasets;
  int numMesh, numSequence;
  FC_Sequence *sequences, temp_sequence, badSequence = { 999, 999 };
  FC_Mesh *meshes, mesh, temp_mesh, badMesh = { 999, 999 };
  int numStep;
  char name[20] = "blue berry", *temp_name;
  FC_Variable *seqVar, *badSeqVar, badVariable = { 999, 999 };

  // get dataset, sequences and meshes
  rc = fc_getDatasetByName(dataset_name, &numReturnDatasets,
			   &returnDatasets);
  fail_unless(rc == FC_SUCCESS, "aborted - cant get dataset by name");
  fail_unless(numReturnDatasets == 1,
	      "failed to find unique dataset by name");
  dataset = returnDatasets[0];
  free(returnDatasets);
  rc = fc_getSequences(dataset, &numSequence, &sequences);
  fail_unless(rc == FC_SUCCESS, "abort: failed to get sequences");
  fail_unless(numSequence == 2, "abort: expected 2 sequences");
  rc = fc_getMeshes(dataset, &numMesh, &meshes);
  fail_unless(rc == FC_SUCCESS, "abort: failed to get meshes");
  fail_unless(numMesh == 8, "abort: expected 8 meshes");

  // loop over all mesh types & different sequences
  for (i = 0; i < numMesh; i++) {
    mesh = meshes[i];
    for (j = 0; j < numSequence; j++) {
      
      // create a seq variable to play with
      rc = fc_createSeqVariable(mesh, sequences[j], name, &numStep, 
				&seqVar);
      fail_unless(rc == FC_SUCCESS, "failed to create seqVar");

      // test fc_isSeqVariableValid()
      fail_unless(fc_isSeqVariableValid(numStep, seqVar), 
		  "failed to validate a valid seq variable");

      // test fc_isVariableGlobal()
      fail_unless(fc_isVariableGlobal(seqVar[0]) == 0, "should be non-global"); 

      // test fc_getVariableName()
      for (k = 0; k < numStep; k++) {
	temp_name = NULL;
	rc = fc_getVariableName(seqVar[k], &temp_name);
	fail_unless(rc == FC_SUCCESS, "failed to get variable name");
	fail_unless(!strcmp(name, temp_name), "mismatch of name");
	// should have returned copy, mess with temp_name & make sure
	temp_name[0] = 'q';
	free(temp_name);
	rc = fc_getVariableName(seqVar[k], &temp_name);
	fail_unless(rc == FC_SUCCESS, "failed to get variable name");
	fail_unless(!strcmp(name, temp_name), "mismatch of name");
	free(temp_name);
      }

       // test fc_getDatasetFromVariable()
      for (k = 0; k < numStep; k++) {
	temp_dataset = badDataset;
	rc = fc_getDatasetFromVariable(seqVar[k], &temp_dataset);
	fail_unless(rc == FC_SUCCESS, "failed to get parent dataset");
	fail_unless(FC_HANDLE_EQUIV(temp_dataset, dataset), 
		    "mismatch of parent dataset");
      }

      // test fc_getMeshFromVariable()
      for (k = 0; k < numStep; k++) {
	temp_mesh = badMesh;
	rc = fc_getMeshFromVariable(seqVar[k], &temp_mesh);
	fail_unless(rc == FC_SUCCESS, "failed to get parent mesh");
	fail_unless(FC_HANDLE_EQUIV(temp_mesh, mesh), 
		    "mismatch of parent mesh");
      }

      // test fc_getSequenceFromSeqVariable()
      temp_sequence = badSequence;
      rc = fc_getSequenceFromSeqVariable(numStep, seqVar,
					      &temp_sequence);
      fail_unless(rc == FC_SUCCESS, "failed to get parent sequence");
      fail_unless(FC_HANDLE_EQUIV(temp_sequence, sequences[j]), 
		    "mismatch of parent sequence");

      // --- check with bad args

      // construct a bad seqVar (make one of the vars bad)
      badSeqVar = (FC_Variable*)malloc(numStep*sizeof(FC_Variable));
      for (k = 0; k < numStep; k++)
	badSeqVar[k] = seqVar[k];
      badSeqVar[numStep/2] = badVariable;

      // fc_isSeqVariableValid()
      fail_unless(!fc_isSeqVariableValid(numStep-1, seqVar),
		  "wrong numStep should NOT be valid");
      fail_unless(!fc_isSeqVariableValid(numStep, NULL),
		  "NULL should NOT be valid");
      fail_unless(!fc_isSeqVariableValid(numStep, badSeqVar),
		  "badSeqVar should NOT be valid");
      fail_unless(!fc_isSeqVariableValid(0, NULL), 
		  "Null seq var should not be valid");

      // already tested fc_getVariableName() & fc_getMeshFromVariable() in
      // variable interface 

      // fc_getSequenceFromSeqVariable()
      temp_sequence = badSequence;
      rc = fc_getSequenceFromSeqVariable(numStep-1, seqVar, &temp_sequence);
      fail_unless(rc != FC_SUCCESS, "wrong numStep should fail");
      fail_unless(FC_HANDLE_EQUIV(temp_sequence, FC_NULL_SEQUENCE),
		  "failure should return NULL sequence");
      rc = fc_getSequenceFromSeqVariable(numStep, NULL, &temp_sequence);
      fail_unless(rc != FC_SUCCESS, "NULL seqVar should fail");
      fail_unless(FC_HANDLE_EQUIV(temp_sequence, FC_NULL_SEQUENCE),
		  "failure should return NULL sequence");
      rc = fc_getSequenceFromSeqVariable(numStep, badSeqVar, &temp_sequence);
      fail_unless(rc != FC_SUCCESS, "bad seqVar should fail");
      fail_unless(FC_HANDLE_EQUIV(temp_sequence, FC_NULL_SEQUENCE),
		  "failure should return NULL sequence");
      rc = fc_getSequenceFromSeqVariable(numStep, seqVar, NULL);
      fail_unless(rc != FC_SUCCESS, "should fail with NULL sequence");

      // --- all done
      
      // one last test
      fc_deleteSeqVariable(numStep, seqVar);
      fail_unless(!fc_isSeqVariableValid(numStep, seqVar),
		  "seq var should not be valid after a delete");

      // cleanup
      free(badSeqVar);
      free(seqVar);
    }
  }

  free(sequences);
  free(meshes);
}
END_TEST

// query the seq variable data by step slice & point slice
// don't bother testing set data/get data interface except where it
// being a sequence variable matters (since fully tested above)
START_TEST(seq_data_query)
{
  FC_ReturnCode rc;
  int i, j, k, m, n, p, q;
  char name[20] = "blue berry", *temp_name;
  FC_Dataset dataset, *returnDatasets;
  FC_Sequence sequence, *returnSequences;
  FC_Mesh mesh, *returnMeshes;
  int numReturnDatasets, numReturnMeshes, numReturnSequences;
  int maxNumEntity = 20; // for all the meshes and all assoc types
  int numStep;
  FC_Variable *seqVar, *badSeqVar;
  FC_Variable badVariable = { 999, 999 };
  int numAssoc = 5, numDataPoint;  // numDataPoint determined later
  FC_AssociationType assocs[5] = { FC_AT_VERTEX, FC_AT_EDGE, FC_AT_FACE,
				   FC_AT_ELEMENT, FC_AT_WHOLE_MESH };
  int numMathType = 4, numComps[4] = { 1, 2, 3, 4 }, maxNumComp = 4;
  FC_MathType mathTypes[4] = { FC_MT_SCALAR, FC_MT_VECTOR, 
			       FC_MT_SYMTENSOR, FC_MT_TENSOR };
  int numDataType = 4;
  FC_DataType dataTypes[4] = { FC_DT_CHAR, FC_DT_INT, 
			       FC_DT_FLOAT, FC_DT_DOUBLE };
  char charData[10][20*4] = { {'a', 'b', 'c', 'd', 'e', 
			       'f', 'g', 'h', 'i', 'j',
			       'k', 'l', 'm', 'n', 'o', 
			       'p', 'q', 'r', 's', 't' } };
  int intData[10][20*4];
  float floatData[10][20*4];
  double doubleData[10][20*4];
  int sizes[4] = {sizeof(char), sizeof(int), sizeof(float), sizeof(double) };
  void *data[4][10];
  int temp_numDataPoint, temp_numComp;
  FC_AssociationType temp_assoc;
  FC_MathType temp_mathType;
  FC_DataType temp_dataType;
  void *temp_data, *temp_history;

  // get dataset and mesh
  rc = fc_getDatasetByName(dataset_name, &numReturnDatasets,
			   &returnDatasets);
  fail_unless(rc == FC_SUCCESS, "aborted - cant get dataset by name");
  fail_unless(numReturnDatasets == 1,
	      "failed to find unique dataset by name");
  dataset = returnDatasets[0];
  free(returnDatasets);
  rc = fc_getSequenceByName(dataset,"time", &numReturnSequences,
			   &returnSequences);
  fail_unless(rc == FC_SUCCESS, "aborted - cant get sequence by name");
  fail_unless(numReturnSequences == 1,
	      "failed to find unique sequences by name");
  sequence = returnSequences[0];
  free(returnSequences);
  rc = fc_getMeshByName(dataset, "hex mesh", &numReturnMeshes,
			&returnMeshes);
  fail_unless(rc == FC_SUCCESS, "aborted - cant get mesh by name");
  fail_unless(numReturnMeshes == 1, "failed to find unique mesh by name");
  mesh = returnMeshes[0];
  free(returnMeshes);

  // setup some arbitrary data (hopefully different for every step, except char)
  fc_getSequenceNumStep(sequence, &numStep);
  for (i = 0; i < numStep; i++) {
    data[0][i] = &charData[i];
    data[1][i] = &intData[i];
    data[2][i] = &floatData[i];
    data[3][i] = &doubleData[i];
    for (j = 0; j < maxNumEntity*maxNumComp; j++) {
      if (i > 0)
	charData[i][j] = charData[0][j%maxNumEntity];
      intData[i][j] = i*j;
      floatData[i][j] = intData[i][j] + 0.1;
      doubleData[i][j] = intData[i][j] + 0.00000000001;
    }
  }

  // association types, math types & data types
  for (j = 0; j < numAssoc; j++) {
    
    // setup numDataPoint
    fc_getMeshNumEntity(mesh, assocs[j], &numDataPoint);
    
    for (k = 0; k < numMathType; k++) {
      for (m = 0; m < numDataType; m++) {
	
	//printf("assoc %s, mathtype %s, datatype %s\n",
	//   fc_getAssociationTypeText(assocs[j]),
	//   fc_getMathTypeText(mathTypes[k]),
	//   fc_getDataTypeText(dataTypes[m]));
	
	// create a sequence variable
	rc = fc_createSeqVariable(mesh, sequence, name, &numStep, &seqVar);
	fail_unless(rc == FC_SUCCESS, "abort: failed to create seqVar");

	// test setting data
	for (n = 0; n < numStep; n++) {
	  rc = fc_setVariableData(seqVar[n], numDataPoint, numComps[k],
				assocs[j],  mathTypes[k], dataTypes[m],
				data[m][n]);
	  fail_unless(rc == FC_SUCCESS, "failed to set variable data");
	}

	// check meta data & get data by step
	for (n = 0; n < numStep; n++) {
	  rc = fc_getVariableName(seqVar[n], &temp_name);
	  fail_unless(rc == FC_SUCCESS, "failed to get name");
	  fail_unless(!strcmp(temp_name, name), "mismatch of name");
	  free(temp_name);
	  rc = fc_getVariableInfo(seqVar[n], &temp_numDataPoint, &temp_numComp,
				  &temp_assoc, &temp_mathType, &temp_dataType);
	  fail_unless(rc == FC_SUCCESS, "failed to get seqVar info");
	  fail_unless(temp_numDataPoint == numDataPoint, 
		      "mismatch of numDataPoint");
	  fail_unless(temp_numComp == numComps[k], "mismatch of numComp");
	  fail_unless(temp_assoc == assocs[j], "mismatch of assoc");
	  fail_unless(temp_mathType == mathTypes[k], "mismatch of mathType");
	  fail_unless(temp_dataType == dataTypes[m], "mismatch of datatype");
	  rc = fc_getVariableDataPtr(seqVar[n], &temp_data);
	  fail_unless(rc == FC_SUCCESS, "failed to get data on a step");
	  fail_unless(!memcmp(temp_data, data[m][n], 
			      numDataPoint*numComps[k]*sizes[m]),
		      "mismatch of data");
	}

	// check get data by data point
	for (n = 0; n < numDataPoint; n++) {
	  rc = fc_getSeqVariableDataPointSlice(numStep, seqVar, n, 
					       &temp_history);
	  fail_unless(rc == FC_SUCCESS, "failed to get history on a point");
	  for (p = 0; p < numStep; p++) {
	    for (q = 0; q < numComps[k]; q++) {
	      if (dataTypes[m] == FC_DT_CHAR) {
		fail_unless(((char*)temp_history)[p*numComps[k]+q] ==
			    charData[p][n*numComps[k]+q],
			    "mismatch of data");	      
	      }
	      else if (dataTypes[m] == FC_DT_INT) {
		fail_unless(((int*)temp_history)[p*numComps[k]+q] ==
			    intData[p][n*numComps[k]+q],
			    "mismatch of data");	      
	      }
	      else if (dataTypes[m] == FC_DT_FLOAT) {
		fail_unless(((float*)temp_history)[p*numComps[k]+q] ==
			    floatData[p][n*numComps[k]+q],
			    "mismatch of data");	      
	      }
	      else if (dataTypes[m] == FC_DT_DOUBLE) {
		fail_unless(((double*)temp_history)[p*numComps[k]+q] ==
			    doubleData[p][n*numComps[k]+q],
			    "mismatch of data");	      
	      }
	      else
		fail_unless(0, "abort: missed a datatype");
	    }
	  }
	  free(temp_history);
	}

	// same as above but for slice as seqvar
	for (n = 0; n < numDataPoint; n++) {
	  FC_Variable *slicevar, **checkvars;
	  void* slicedata;
	  FC_DataType slicedt;
	  FC_MathType slicemt;
	  FC_AssociationType sliceat;
	  FC_Sequence sliceseq;
	  int slicenc, slicedp;
	  int numcheckvars, *checksteps;

	  rc = fc_getSeqVariableDataPointSliceAsSeqVariable(numStep,
				    seqVar, n, "slicevarname", &slicevar);
	  fail_unless(rc == FC_SUCCESS, "failed to get slicevar");
	  fail_unless(fc_isSeqVariableValid(numStep,slicevar),
		"didn't create valid seqvar");
	  fc_getVariableInfo(slicevar[0],&slicedp,&slicenc,
		       &sliceat,&slicemt,&slicedt);
	  fail_unless(slicedp == 1,"mismatch of numDataPoint");
	  fail_unless(slicenc == numComps[k],"mismatch of numComponent");
	  fail_unless(sliceat == FC_AT_WHOLE_MESH,"mismatch of assoc");
	  fail_unless(slicemt == mathTypes[k],"mismatch of mathtype");
	  fail_unless(slicedt == dataTypes[m],"mismatch of datatype");

	  fc_getSequenceFromSeqVariable(numStep,slicevar,&sliceseq);
	  fail_unless(FC_HANDLE_EQUIV(sliceseq,sequence),
		      "mismatch of sequence handles");
	  fc_getSeqVariableByName(mesh,"slicevarname",&numcheckvars,
				  &checksteps,&checkvars);
	  fail_unless(numcheckvars ==1, "wrong number of matching vars");
	  for (i = 0; i < numStep; i++){
	    fail_unless(FC_HANDLE_EQUIV(slicevar[i],checkvars[0][i]),
			"mismatch of variable handles");
	  }
	  free(checkvars[0]);
	  free(checkvars);
	  free(checksteps);

	  for (p = 0; p < numStep; p++) {
	    rc = fc_getVariableDataPtr(slicevar[p],&slicedata);
	    fail_unless(rc == FC_SUCCESS, "failed to get slicevar data");

	    for (q = 0; q < numComps[k]; q++) {
	      if (dataTypes[m] == FC_DT_CHAR) {
		fail_unless(((char*)slicedata)[q] ==
			    charData[p][n*numComps[k]+q],
			    "mismatch of data");	      
	      }
	      else if (dataTypes[m] == FC_DT_INT) {
		fail_unless(((int*)slicedata)[q] ==
			    intData[p][n*numComps[k]+q],
			    "mismatch of data");	      
	      }
	      else if (dataTypes[m] == FC_DT_FLOAT) {
		fail_unless(((float*)slicedata)[q] ==
			    floatData[p][n*numComps[k]+q],
			    "mismatch of data");	      
	      }
	      else if (dataTypes[m] == FC_DT_DOUBLE) {
		fail_unless(((double*)slicedata)[q] ==
			    doubleData[p][n*numComps[k]+q],
			    "mismatch of data");	      
	      }
	      else
		fail_unless(0, "abort: missed a datatype");
	    }
	  }
	  fc_deleteSeqVariable(numStep,slicevar);
	  free(slicevar);
	}

	 
	// ----- check errors

	// construct a bad seqVar (make one of the vars bad)
	badSeqVar = (FC_Variable*)malloc(numStep*sizeof(FC_Variable));
	for (n = 0; n < numStep; n++)
	  badSeqVar[n] = seqVar[n];
	badSeqVar[numStep/2] = badVariable;

	// fc_getSeqVariableDataPointSlice() -- bad args
	rc = fc_getSeqVariableDataPointSlice(numStep-1, seqVar, 0, 
					     &temp_history);
	fail_unless(rc != FC_SUCCESS, "should fail with wrong numStep");
	fail_unless(temp_history == NULL, "fail should return null");
	rc = fc_getSeqVariableDataPointSlice(numStep, NULL, 0, &temp_history);
	fail_unless(rc != FC_SUCCESS, "should fail with null seqVar");
	fail_unless(temp_history == NULL, "fail should return null");
	rc = fc_getSeqVariableDataPointSlice(numStep, badSeqVar, 0, 
					     &temp_history);
	fail_unless(rc != FC_SUCCESS, "should fail with bad seqVar");
	fail_unless(temp_history == NULL, "fail should return null");
	rc = fc_getSeqVariableDataPointSlice(numStep, seqVar, -1, 
					     &temp_history);
	fail_unless(rc != FC_SUCCESS, "should failed if ID < 0");
	fail_unless(temp_history == NULL, "fail should return null");
	rc = fc_getSeqVariableDataPointSlice(numStep, seqVar, numDataPoint, 
					     &temp_history);
	fail_unless(rc != FC_SUCCESS, "should failed if ID too big");
	fail_unless(temp_history == NULL, "fail should return null");
	rc = fc_getSeqVariableDataPointSlice(numStep, seqVar, 0, NULL); 
	fail_unless(rc != FC_SUCCESS, "should failed if no history");

	{
	  // fc_getSeqVariableDataPointSliceAsSeqVariable() -- bad args
	  FC_Variable* slicevar;
	  char* slicename = "slicename";

	  rc = fc_getSeqVariableDataPointSliceAsSeqVariable(numStep-1,
					    seqVar, 0, slicename,&slicevar);
	  fail_unless(rc != FC_SUCCESS, "should fail with wrong numStep");
	  fail_unless(slicevar == NULL, "fail should return null");
	  rc = fc_getSeqVariableDataPointSliceAsSeqVariable(numStep,
                                            NULL, 0, slicename,&slicevar);
	  fail_unless(rc != FC_SUCCESS, "should fail with null seqVar");
	  fail_unless(slicevar == NULL, "fail should return null");
	  rc = fc_getSeqVariableDataPointSliceAsSeqVariable(numStep,
                                            badSeqVar, 0, slicename,&slicevar);
	  fail_unless(rc != FC_SUCCESS, "should fail with bad seqVar");
	  fail_unless(slicevar == NULL, "fail should return null");
	  rc = fc_getSeqVariableDataPointSliceAsSeqVariable(numStep, seqVar,
                                            -1, slicename,&slicevar);
	  fail_unless(rc != FC_SUCCESS, "should failed if ID < 0");
	  fail_unless(slicevar == NULL, "fail should return null");
	  rc = fc_getSeqVariableDataPointSliceAsSeqVariable(numStep, seqVar,
                                            numDataPoint, 
                                            slicename,&slicevar);
	  fail_unless(rc != FC_SUCCESS, "should failed if ID too big");
	  fail_unless(slicevar == NULL, "fail should return null");
	  rc = fc_getSeqVariableDataPointSliceAsSeqVariable(numStep, seqVar,
							  0, NULL, &slicevar);
	  fail_unless(rc != FC_SUCCESS, "should failed if NULL var name");
	  fail_unless(slicevar == NULL, "fail should return null");
	  rc = fc_getSeqVariableDataPointSliceAsSeqVariable(numStep, seqVar,
							  0, slicename,
							  NULL);
	  fail_unless(rc != FC_SUCCESS, "should failed if NULL var handle");
	  fail_unless(slicevar == NULL, "fail should return null");
	}

	// cleanup
	rc = fc_deleteSeqVariable(numStep, seqVar);
	fail_unless(rc == FC_SUCCESS, "failed to delete seqVar");
	free(badSeqVar);
	free(seqVar);
      }
    }
  }
}
END_TEST

// copy
START_TEST(seq_copy_test)
{
  FC_ReturnCode rc;
  int i, j, k, m, n;
  char name[20] = "blue berry", copy_name[20] = "banana nut", *temp_name;
  FC_Dataset dataset, *returnDatasets;
  FC_Sequence sequence1, sequence2, weirdSequence, badSequence = { 999, 999 };
  FC_Sequence *returnSequences;
  FC_Mesh mesh1, mesh2, weirdMesh, badMesh = { 999, 999 };
  FC_Mesh *returnMeshes;
  int numReturnDatasets, numReturnMeshes, numReturnSequences;
  int maxNumEntity = 20; // for all the meshes and all assoc types
  int numSeqVar1, numSeqVar2, temp_numSeqVar, numStep;
  FC_Variable *seqVar, *copy_seqVar, *badSeqVar;
  FC_Variable badVariable = { 999, 999 };
  FC_Variable *glbSeqVar;
  int numAssoc = 5, numDataPoint;  // numDataPoint determined later
  FC_AssociationType assocs[5] = { FC_AT_VERTEX, FC_AT_EDGE, FC_AT_FACE,
				   FC_AT_ELEMENT, FC_AT_WHOLE_MESH };
  int numMathType = 4, numComps[4] = { 1, 2, 3, 4 }, maxNumComp = 4;
  FC_MathType mathTypes[4] = { FC_MT_SCALAR, FC_MT_VECTOR, 
			       FC_MT_SYMTENSOR, FC_MT_TENSOR };
  int numDataType = 4;
  FC_DataType dataTypes[4] = { FC_DT_CHAR, FC_DT_INT, 
			       FC_DT_FLOAT, FC_DT_DOUBLE };
  char charData[10][20*4] = { {'a', 'b', 'c', 'd', 'e', 
			       'f', 'g', 'h', 'i', 'j',
			       'k', 'l', 'm', 'n', 'o', 
			       'p', 'q', 'r', 's', 't' } };
  int intData[10][20*4];
  float floatData[10][20*4];
  double doubleData[10][20*4];
  int sizes[4] = {sizeof(char), sizeof(int), sizeof(float), sizeof(double) };
  void *data[4][10];
  int temp_numDataPoint, temp_numComp;
  FC_AssociationType temp_assoc;
  FC_MathType temp_mathType;
  FC_DataType temp_dataType;
  void *temp_data;


  // get dataset and mesh
  rc = fc_getDatasetByName(dataset_name, &numReturnDatasets,
			   &returnDatasets);
  fail_unless(rc == FC_SUCCESS, "aborted - cant get dataset by name");
  fail_unless(numReturnDatasets == 1,
	      "failed to find unique dataset by name");
  dataset = returnDatasets[0];
  free(returnDatasets);
  rc = fc_getSequenceByName(dataset,"time", &numReturnSequences,
			   &returnSequences);
  fail_unless(rc == FC_SUCCESS, "aborted - cant get sequence by name");
  fail_unless(numReturnSequences == 1,
	      "failed to find unique sequences by name");
  sequence1 = returnSequences[0];
  free(returnSequences);
  rc = fc_getMeshByName(dataset, "hex mesh", &numReturnMeshes,
			&returnMeshes);
  fail_unless(rc == FC_SUCCESS, "aborted - cant get mesh by name");
  fail_unless(numReturnMeshes == 1, "failed to find unique mesh by name");
  mesh1 = returnMeshes[0];
  free(returnMeshes);

  // make a second sequence & second mesh to copy to
  rc = fc_copySequence(sequence1, dataset, "copy sequence", &sequence2);
  fail_unless(rc == FC_SUCCESS, "abort: failed to copy sequence");
  rc = fc_copyMesh(mesh1, dataset, "copy mesh", &mesh2);
  fail_unless(rc == FC_SUCCESS, "abort: failed to copy mesh");

  // get a third sequence we shouldn't be able to copy to - different numStep
  // (don't delete when we're done)
  rc = fc_getSequenceByName(dataset,"time2", &numReturnSequences,
			   &returnSequences);
  fail_unless(rc == FC_SUCCESS, "aborted - cant get sequence by name");
  fail_unless(numReturnSequences == 1,
	      "failed to find unique sequences by name");
  weirdSequence = returnSequences[0];
  free(returnSequences);

  // make a third mesh we shouldn't be able to copy to - just 1 element
  {
    int topodim, dim, numVertex, numElement, *conns;
    double* coords;
    FC_ElementType elemType;
    fc_getMeshInfo(mesh1, &topodim, &dim, &numVertex, &numElement, &elemType);
    fc_getMeshCoordsPtr(mesh1, &coords);
    fc_getMeshElementConnsPtr(mesh1, &conns);
    rc = fc_createMesh(dataset, "different mesh", &weirdMesh);
    fail_unless(rc == FC_SUCCESS, "abort: failed to create mesh");
    rc = fc_setMeshCoords(weirdMesh, dim, fc_getElementTypeNumVertex(elemType),
			  coords);
    fail_unless(rc == FC_SUCCESS, "abort: failed to set coords");
    rc = fc_setMeshElementConns(weirdMesh, elemType, 1, conns);
    fail_unless(rc == FC_SUCCESS, "abort: failed to set conns");
  }

  // setup some arbitrary data (hopefully different for every step, except char)
  fc_getSequenceNumStep(sequence1, &numStep);
  for (i = 0; i < numStep; i++) {
    data[0][i] = &charData[i];
    data[1][i] = &intData[i];
    data[2][i] = &floatData[i];
    data[3][i] = &doubleData[i];
    for (j = 0; j < maxNumEntity*maxNumComp; j++) {
      if (j > 0)
	charData[i][j] = charData[0][j%maxNumEntity];
      intData[i][j] = i*j;
      floatData[i][j] = intData[i][j] + 0.1;
      doubleData[i][j] = intData[i][j] + 0.00000000001;
    }
  }

  // create a global sequence variable (double)
  rc = fc_createGlobalSeqVariable(dataset, sequence1, name, &numStep, &glbSeqVar);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create glb seqVar");
  for (n = 0; n < numStep; n++) {
    rc = fc_setVariableData(glbSeqVar[n], 1, 1, FC_AT_WHOLE_DATASET,
                            FC_MT_SCALAR, FC_DT_DOUBLE, doubleData[n]);
    fail_unless(rc == FC_SUCCESS, "abort: failed to set data");
  }

  numSeqVar1 = 0;
  numSeqVar2 = 0;
  // association types, math types & data types
  for (j = 0; j < numAssoc; j++) {
    
    // setup numDataPoint
    fc_getMeshNumEntity(mesh1, assocs[j], &numDataPoint);
    
    for (k = 0; k < numMathType; k++) {
      for (m = 0; m < numDataType; m++) {
	
	//printf("assoc %s, mathtype %s, datatype %s\n",
	//     fc_getAssociationTypeText(assocs[j]),
	//     fc_getMathTypeText(mathTypes[k]),
	//     fc_getDataTypeText(dataTypes[m]));
	
	// create a sequence variable
	rc = fc_createSeqVariable(mesh1, sequence1, name, &numStep, &seqVar);
	fail_unless(rc == FC_SUCCESS, "abort: failed to create seqVar");
	for (n = 0; n < numStep; n++) {
	  rc = fc_setVariableData(seqVar[n], numDataPoint, numComps[k],
				assocs[j],  mathTypes[k], dataTypes[m],
				data[m][n]);
	  fail_unless(rc == FC_SUCCESS, "abort: failed to set data");
	}
	numSeqVar1++;

	// copy to same mesh & sequence & check it
	rc = fc_copySeqVariable(numStep, seqVar, mesh1, sequence1, 
				copy_name, &copy_seqVar);
	fail_unless(rc == FC_SUCCESS, "failed to copy to same mesh");
	numSeqVar1++;
	for (n = 0; n < numStep; n++) {
	  rc = fc_getVariableName(copy_seqVar[n], &temp_name);
	  fail_unless(rc == FC_SUCCESS, "failed to get name of copy");
	  fail_unless(!strcmp(temp_name, copy_name), "mismatch of name");
	  free(temp_name);
	  rc = fc_getVariableInfo(copy_seqVar[n], &temp_numDataPoint, 
				  &temp_numComp, &temp_assoc, &temp_mathType, 
				  &temp_dataType);
	  fail_unless(rc == FC_SUCCESS, "failed to get seqVar info");
	  fail_unless(temp_numDataPoint == numDataPoint, 
		      "mismatch of numDataPoint");
	  fail_unless(temp_numComp == numComps[k], "mismatch of numComp");
	  fail_unless(temp_assoc == assocs[j], "mismatch of assoc");
	  fail_unless(temp_mathType == mathTypes[k], "mismatch of mathType");
	  fail_unless(temp_dataType == dataTypes[m], "mismatch of datatype");
	  rc = fc_getVariableDataPtr(copy_seqVar[n], &temp_data);
	  fail_unless(rc == FC_SUCCESS, "failed to get data");
	  fail_unless(!memcmp(temp_data, data[m][n], 
			      numDataPoint*numComps[k]*sizes[m]),
		      "mismatch of data");
	}
	free(copy_seqVar);

	// copy to different mesh & different sequence & check it
	rc = fc_copySeqVariable(numStep, seqVar, mesh2, sequence2, 
				copy_name, &copy_seqVar);
	fail_unless(rc == FC_SUCCESS, "failed to copy to same mesh");
	numSeqVar2++;
	for (n = 0; n < numStep; n++) {
	  rc = fc_getVariableName(copy_seqVar[n], &temp_name);
	  fail_unless(rc == FC_SUCCESS, "failed to get name of copy");
	  fail_unless(!strcmp(temp_name, copy_name), "mismatch of name");
	  free(temp_name);
	  rc = fc_getVariableInfo(copy_seqVar[n], &temp_numDataPoint, 
				  &temp_numComp, &temp_assoc, &temp_mathType, 
				  &temp_dataType);
	  fail_unless(rc == FC_SUCCESS, "failed to get seqVar info");
	  fail_unless(temp_numDataPoint == numDataPoint, 
		      "mismatch of numDataPoint");
	  fail_unless(temp_numComp == numComps[k], "mismatch of numComp");
	  fail_unless(temp_assoc == assocs[j], "mismatch of assoc");
	  fail_unless(temp_mathType == mathTypes[k], "mismatch of mathType");
	  fail_unless(temp_dataType == dataTypes[m], "mismatch of datatype");
	  rc = fc_getVariableDataPtr(copy_seqVar[n], &temp_data);
	  fail_unless(rc == FC_SUCCESS, "failed to get data");
	  fail_unless(!memcmp(temp_data, data[m][n], 
			      numDataPoint*numComps[k]*sizes[m]),
		      "mismatch of data");
	}
	free(copy_seqVar);

	// fc_copySeqVariable() to very differnt mesh should fail for non WHOLE
	if (assocs[j] != FC_AT_WHOLE_MESH) {
	  rc = fc_copySeqVariable(numStep, seqVar, weirdMesh, sequence1, 
				  copy_name, &copy_seqVar);
	  fail_unless(rc != FC_SUCCESS,
		      "should fail to copy to a very different mesh");
	  fail_unless(copy_seqVar == NULL, "should return null seqvar");
	}

	// fc_copySeqVariable() to a very different sequence should fail
	rc = fc_copySeqVariable(numStep, seqVar, mesh1, weirdSequence, 
				copy_name, &copy_seqVar);
	fail_unless(rc != FC_SUCCESS,
		    "should fail to copy to a very different mesh");
	fail_unless(copy_seqVar == NULL, "should return null seqvar");

	// --- special cases
	
	// copy with NULL name will use original's name
	rc = fc_copySeqVariable(numStep, seqVar, mesh1, sequence1, 
				NULL, &copy_seqVar);
	fail_unless(rc == FC_SUCCESS, "failed to copy to same mesh");
	for (n = 0; n < numStep; n++) {
	  rc = fc_getVariableName(copy_seqVar[n], &temp_name);
	  fail_unless(rc == FC_SUCCESS, "failed to get name of copy");
	  fail_unless(!strcmp(temp_name, name), "mismatch of name");
	  free(temp_name);
	}
	rc = fc_deleteSeqVariable(numStep, copy_seqVar);
	fail_unless(rc == FC_SUCCESS, "failed to delete copied seqVar");
	free(copy_seqVar);

	// ---- check errors
	
	// construct a bad seqVar (make one of the vars bad)
	badSeqVar = (FC_Variable*)malloc(numStep*sizeof(FC_Variable));
	for (n = 0; n < numStep; n++)
	  badSeqVar[n] = seqVar[n];
	badSeqVar[numStep/2] = badVariable;

	// fc_copyVariable() -- bad args
	rc = fc_copySeqVariable(numStep-1, seqVar, mesh1, sequence1, 
				copy_name, &copy_seqVar);
	fail_unless(rc != FC_SUCCESS, "should fail with wrong numStep");
	fail_unless(copy_seqVar == NULL, "fail should return null seqVar");
	rc = fc_copySeqVariable(numStep, NULL, mesh1, sequence1, 
				copy_name, &copy_seqVar);
	fail_unless(rc != FC_SUCCESS, "should fail with null seqVar");
	fail_unless(copy_seqVar == NULL, "fail should return null seqVar");
	rc = fc_copySeqVariable(numStep, badSeqVar, mesh1, sequence1, 
				copy_name, &copy_seqVar);
	fail_unless(rc != FC_SUCCESS, "should fail with badSeqVar");
	fail_unless(copy_seqVar == NULL, "fail should return null seqVar");
	rc = fc_copySeqVariable(numStep, seqVar, badMesh, sequence1, 
				copy_name, &copy_seqVar);
	fail_unless(rc != FC_SUCCESS, "should fail with bad mesh");
	fail_unless(copy_seqVar == NULL, "fail should return null seqVar");
	rc = fc_copySeqVariable(numStep, seqVar, mesh1, badSequence, 
				copy_name, &copy_seqVar);
	fail_unless(rc != FC_SUCCESS, "should fail with bad sequence");
	fail_unless(copy_seqVar == NULL, "fail should return null seqVar");
	rc = fc_copySeqVariable(numStep, seqVar, mesh1, sequence1, 
				copy_name, NULL);
	fail_unless(rc != FC_SUCCESS, "should fail with null seqvar");

	// --- done
	
	// cleanup
	free(seqVar);
	free(badSeqVar);
      }
    }
  }

  // it's o.k. to copy a global seq var onto a mesh
  rc = fc_copySeqVariable(numStep, glbSeqVar, mesh1, sequence1, 
                          copy_name, &copy_seqVar);
  fail_unless(rc == FC_SUCCESS, "failed to copy glb seq var to a mesh");
  numSeqVar1++;
  for (n = 0; n < numStep; n++) {
    rc = fc_getVariableName(copy_seqVar[n], &temp_name);
    fail_unless(rc == FC_SUCCESS, "failed to get name of copy");
    fail_unless(!strcmp(temp_name, copy_name), "mismatch of name");
    free(temp_name);
    rc = fc_getVariableInfo(copy_seqVar[n], &temp_numDataPoint, 
                            &temp_numComp, &temp_assoc, &temp_mathType, 
                            &temp_dataType);
    fail_unless(rc == FC_SUCCESS, "failed to get seqVar info");
    fail_unless(temp_numDataPoint == 1, "mismatch of numDataPoint");
    fail_unless(temp_numComp == 1, "mismatch of numComp");
    fail_unless(temp_assoc == FC_AT_WHOLE_MESH, "mismatch of assoc");
    fail_unless(temp_mathType == FC_MT_SCALAR, "mismatch of mathType");
    fail_unless(temp_dataType == FC_DT_DOUBLE, "mismatch of datatype");
    rc = fc_getVariableDataPtr(copy_seqVar[n], &temp_data);
    fail_unless(rc == FC_SUCCESS, "failed to get data");
    fail_unless(!memcmp(temp_data, doubleData[n], sizeof(double)), 
                "mismatch of data");
  }
  free(copy_seqVar);

  // a little more checking & cleanup
  rc = fc_deleteMesh(weirdMesh);
  fail_unless(rc == FC_SUCCESS, "failed to delete weird mesh");  
  rc = fc_getNumSeqVariable(mesh2, &temp_numSeqVar);
  fail_unless(rc == FC_SUCCESS, "shouldn't fail to get numSeqVars");
  fail_unless(temp_numSeqVar == numSeqVar2, "mismatch of numSeqVar");
  rc = fc_deleteMesh(mesh2);
  fail_unless(rc == FC_SUCCESS, "failed to delete second mesh");
  rc = fc_getNumSeqVariable(mesh1, &temp_numSeqVar);
  fail_unless(rc == FC_SUCCESS, "shouldn't fail to get numSeqVar");
  fail_unless(temp_numSeqVar == numSeqVar1, "mismatch of numSeqVar");

  // cleanup
  free(glbSeqVar);
}
END_TEST

// release vars
START_TEST(seq_release_var)
{
  int i;
  FC_ReturnCode rc;
  FC_Dataset dataset;
  FC_Sequence sequence;
  FC_Mesh mesh, *returnMeshes;
  int numReturnMeshes;
  int numStep;
  FC_Variable *seqVar, *seqVar2, *badSeqVar, badVar = { 999, 999 };
  FC_Variable **returnSeqVars;
  int numReturnSeqVars, *returnSeqVarsSteps;
  _FC_VarSlot* varSlot;
  void *data;

  // setup
  rc = fc_loadDataset("../data/gen_hex_seq.ex2", &dataset);
  fail_unless(rc == FC_SUCCESS, "abort: failed to load dataset");
  rc = fc_getMeshByName(dataset, "hex mesh", &numReturnMeshes,
			&returnMeshes);
  fail_unless(rc == FC_SUCCESS, "aborted - cant get mesh by name");
  fail_unless(numReturnMeshes == 1, "failed to find unique mesh by name");
  mesh = returnMeshes[0];
  free(returnMeshes);

  rc = fc_getSeqVariableByName(mesh, "temperature per vertex", &numReturnSeqVars,
			       &returnSeqVarsSteps, &returnSeqVars);
  fail_unless(rc == FC_SUCCESS, "abort: failed to get var");
  fail_unless(numReturnSeqVars == 1, "wrong number of matching vars");
  numStep = returnSeqVarsSteps[0];
  seqVar = returnSeqVars[0];
  free(returnSeqVarsSteps);
  free(returnSeqVars);

  rc = fc_getSequenceFromSeqVariable(numStep, seqVar, &sequence);
  fail_unless(rc == FC_SUCCESS, "abort: failed to get sequence");

  // test that lazily loaded data is not there yet
  for (i = 0; i < numStep; i++) {
    varSlot = _fc_getVarSlot(seqVar[i]);
    fail_unless(varSlot->data == NULL, "abort: should start with no data");
  }
  
  // make copy (uncommitted) variable
  rc = fc_copySeqVariable(numStep, seqVar, mesh, sequence, "new_var", &seqVar2);
  fail_unless(rc == FC_SUCCESS, "failed to copy var");

  // test that created data is there
  // --- data

  // load & unload variable data
  for (i = 0; i < numStep; i++) {
    varSlot = _fc_getVarSlot(seqVar[i]);
    rc = fc_getVariableDataPtr(seqVar[i], &data); // force loading
    fail_unless(rc == FC_SUCCESS, "failed to get data");
    fail_unless(varSlot->data != NULL, "abort: data should have loaded");
    fail_unless(data == varSlot->data, "mismatch of data pointers");
  }
  rc = fc_releaseSeqVariable(numStep, seqVar);
  fail_unless(rc == FC_SUCCESS, "failed to release data");
  for (i = 0; i < numStep; i++) {
    varSlot = _fc_getVarSlot(seqVar[i]);
    fail_unless(varSlot->data == NULL, "data should now be null");
  }

  // uncommitted var's shouldn't get deleted
  for (i = 0; i < numStep; i++) { 
    varSlot = _fc_getVarSlot(seqVar2[i]);
    fail_unless(varSlot->data != NULL, "abort: start with data");
  }
  fc_releaseSeqVariable(numStep, seqVar); // make sure that copy doesn't just copy pointers
  rc = fc_releaseSeqVariable(numStep, seqVar2);
  fail_unless(rc == FC_SUCCESS, "failed to release data");
  for (i = 0; i < numStep; i++) { 
    varSlot = _fc_getVarSlot(seqVar2[i]);
    fail_unless(varSlot->data != NULL, "should not delete temp var data");
  }

  // --- special cases
  
  // releasing a null handle does not cause an error
  rc = fc_releaseSeqVariable(0, NULL);
  fail_unless(rc == FC_SUCCESS, "NULL handle should not fail");
  
  // --- errors

  // construct a bad seqVar (make one of the vars bad)
  badSeqVar = (FC_Variable*)malloc(numStep*sizeof(FC_Variable));
  for (i = 0; i < numStep; i++)
    badSeqVar[i] = seqVar[i];
  badSeqVar[numStep/2] = badVar;
  
  // fc_releaseSeqVariable() --- bad args
  rc = fc_releaseSeqVariable(numStep, badSeqVar);
  fail_unless(rc != FC_SUCCESS, "bad variable should fail");

  // cleanup
  free(seqVar);
  free(seqVar2);
  free(badSeqVar);
  fc_deleteDataset(dataset);
}
END_TEST 

// --- global variable interface 

// only stuff specific to global variables is tested here since the
// basic variable stuff is test in the variable interface tests above

// test creating and deleting global variables, also test querying dataset
// for lobal vars before, inbetween and after.
START_TEST(glb_create_get_delete)
{
  FC_ReturnCode rc;
  int j;
  FC_Dataset dataset, *returnDatasets;
  int numReturnDatasets;
  int numVariable = 10, temp_numVariable;
  FC_Variable variables[10], *temp_variables, temp_variable;
  FC_Variable badVariable = { 999, 999 };
  char names[10][100] = { "one", "two", "three", "four", "five", "six",
			  "seven", "eight", "nine", "ten" };
  char newName[100] = { "sparkly new" };

  // get dataset
  rc = fc_getDatasetByName(dataset_name, &numReturnDatasets,
			   &returnDatasets);
  fail_unless(rc == FC_SUCCESS, "aborted - cant get dataset by name");
  fail_unless(numReturnDatasets == 1,
	      "failed to find unique dataset by name");
  dataset = returnDatasets[0];
  free(returnDatasets);

  // create some global variables
  for (j = 0; j < numVariable; j++) {
    rc = fc_createGlobalVariable(dataset, names[j], &variables[j]);
    fail_unless(rc == FC_SUCCESS, "failed to create variable");
    fail_unless(!FC_HANDLE_EQUIV(variables[j], FC_NULL_VARIABLE),
		"created variable should not = FC_NULL_VARIABLE");
  }
 
  // get global variables (should be numVariable);
  rc = fc_getGlobalVariables(dataset, &temp_numVariable, &temp_variables);
  fail_unless(rc == FC_SUCCESS, "failed to get variables");
  fail_unless(temp_numVariable == numVariable, "mismatch of numVariable");
  for (j = 0; j < numVariable; j++) 
    fail_unless(FC_HANDLE_EQUIV(temp_variables[j], variables[j]), 
		"mismatch of variable handles");
  free(temp_variables);
  rc = fc_getNumGlobalVariable(dataset, &temp_numVariable);
  fail_unless(rc == FC_SUCCESS, "failed to get variables");
  fail_unless(temp_numVariable == numVariable, "mismatch of numVariable");
  for (j = 0; j < numVariable; j++) {
    rc = fc_getGlobalVariableByName(dataset, names[j], &temp_numVariable, &temp_variables);
    fail_unless(rc == FC_SUCCESS, "failed to get variable by name");
    fail_unless(temp_numVariable == 1, "failed to find unique variable by name");
    fail_unless(FC_HANDLE_EQUIV(temp_variables[0], variables[j]),
		"mismatch of variable handles");
    free(temp_variables);
  }

  // test negative example of fc_getGlobalVariableByName
  rc = fc_getGlobalVariableByName(dataset, newName, &temp_numVariable, &temp_variables);
  fail_unless(rc == FC_SUCCESS, "failed to get variable by name");
  fail_unless(temp_numVariable == 0, "should fail to find unique variable by name");
  fail_unless(temp_variables == NULL, "should return NULL");

  //make a temporary variable to check multivariable return
  rc = fc_createGlobalVariable(dataset, names[0], &temp_variable);
  fail_unless(rc == FC_SUCCESS, "failed to create variable");
  fail_unless(!FC_HANDLE_EQUIV(temp_variable, FC_NULL_VARIABLE),
	      "created variable should not = FC_NULL_VARIABLE");
  rc = fc_getGlobalVariableByName(dataset, names[0], &temp_numVariable, &temp_variables);
  fail_unless(rc == FC_SUCCESS, "failed to get variable by name");
  fail_unless(temp_numVariable == 2, "should fail to find unique variable by name");
  fail_unless((FC_HANDLE_EQUIV(temp_variables[0], variables[0]) &&
	       FC_HANDLE_EQUIV(temp_variables[1], temp_variable)) ||
	      (FC_HANDLE_EQUIV(temp_variables[1], variables[0]) &&
	       FC_HANDLE_EQUIV(temp_variables[0], temp_variable)),
	      "mismatch of variable handles");
  fc_deleteVariable(temp_variable);
  free(temp_variables);

  // delete half of the variables (alternate)
  for (j = 0; j < numVariable; j+=2) {
    rc = fc_deleteVariable(variables[j]);
    fail_unless(rc == FC_SUCCESS, "failed to delete variable");
  }
 
  // get global variables (should be numVariable/2);
  fc_getGlobalVariables(dataset, &temp_numVariable, &temp_variables);
  fail_unless(rc == FC_SUCCESS, "failed to get variables");
  fail_unless(temp_numVariable == numVariable/2, "mismatch of numVariable");
  for (j = 0; j < numVariable/2; j++) 
    fail_unless(FC_HANDLE_EQUIV(temp_variables[j], variables[j*2+1]), 
		"mismatch of variable handles");
  free(temp_variables);
  for (j = 0; j < numVariable/2; j++) {
    rc = fc_getGlobalVariableByName(dataset, names[j*2+1], 
				    &temp_numVariable, &temp_variables);
    fail_unless(rc == FC_SUCCESS, "failed to get variable by name");
    fail_unless(temp_numVariable == 1, "failed to find unique variable by name");
    fail_unless(FC_HANDLE_EQUIV(temp_variables[0], variables[j*2+1]),
		"mismatch of variable handles");
    free(temp_variables);
  }
 
  // delete remaining variables
  for (j = 1; j < numVariable; j+=2) { 
    rc = fc_deleteVariable(variables[j]);
    fail_unless(rc == FC_SUCCESS, "failed to delete variable");
  }
  
  // get variables (should be none)
  rc = fc_getGlobalVariables(dataset, &temp_numVariable, &temp_variables);
  fail_unless(rc == FC_SUCCESS, 
	      "failed to get variables from empty library");
  fail_unless(temp_numVariable == 0 && temp_variables == NULL,
	      "should return 0 if all variables deleted");
  
  // make one last variable for further testing
  rc = fc_createGlobalVariable(dataset, names[0], &variables[0]);
  fail_unless(rc == FC_SUCCESS, 
	      "aborted: failed to create variable");
 
  // ---- test error conditions
  
  // bad args to fc_createGlobalVariable()
  temp_variable = badVariable;
  rc = fc_createGlobalVariable(FC_NULL_DATASET, names[0], &temp_variable);
  fail_unless(rc != FC_SUCCESS, 
	      "should have failed to create variable with null database");
  fail_unless(FC_HANDLE_EQUIV(temp_variable, FC_NULL_VARIABLE),
	      "fail should return NULL variable");
  temp_variable = badVariable;
  rc = fc_createGlobalVariable(dataset, NULL, &temp_variable);
  fail_unless(rc != FC_SUCCESS, 
	      "should have failed to create variable with null name");
  fail_unless(FC_HANDLE_EQUIV(temp_variable, FC_NULL_VARIABLE),
	      "fail should return NULL variable");
  rc = fc_createGlobalVariable(dataset, names[0], NULL);
  fail_unless(rc != FC_SUCCESS, 
	      "should have failed to create variable with null handle");

  // bad args to fc_getGlobalVariables()
  temp_numVariable = 99;
  temp_variables = (FC_Variable*)1;
  rc = fc_getGlobalVariables(FC_NULL_DATASET, &temp_numVariable, &temp_variables);
  fail_unless(rc != FC_SUCCESS, 
	      "should have failed to create variable with null database");
  fail_unless(temp_numVariable == -1,
	      "fail should return -1 for numVariable"); 
  fail_unless(temp_variables == NULL, "fail should return NULL");
  temp_variables = (FC_Variable*)1;
  rc = fc_getGlobalVariables(dataset, NULL, &temp_variables);
  fail_unless(rc != FC_SUCCESS, 
	      "should have failed to get variable with NULL numVariable");
  fail_unless(temp_variables == NULL, "fail should return null");
  temp_numVariable = 999;
  rc = fc_getGlobalVariables(dataset, &temp_numVariable, NULL);
  fail_unless(rc != FC_SUCCESS, "should be error to get just numVariable");
  fail_unless(temp_numVariable == -1, "mismatch of numVariable");    
    
  // bad args to fc_getNumGlobalVariable()
  temp_numVariable = 99;
  rc = fc_getNumGlobalVariable(FC_NULL_DATASET, &temp_numVariable);
  fail_unless(rc != FC_SUCCESS, 
	      "should have failed to get numVariable with null dataset");
  fail_unless(temp_numVariable == -1,
	      "fail should return -1 for numVariable"); 
  rc = fc_getNumGlobalVariable(dataset, NULL);
  fail_unless(rc != FC_SUCCESS, 
	      "should have failed to get numVariable with NULL numVariable");
  
  // --- done
  
  // delete last variable
  rc = fc_deleteVariable(variables[0]);
  fail_unless(rc == FC_SUCCESS, "failed to delete last variable");
}
END_TEST

// query meta data
START_TEST(glb_metadata_query)
{
  FC_ReturnCode rc;
  FC_Dataset dataset, temp_dataset, *returnDatasets;
  int numReturnDatasets;
  FC_Mesh temp_mesh, badMesh = { 999, 999 };
  char name[20] = "blue berry";
  FC_Variable variable;

  // get dataset
  rc = fc_getDatasetByName(dataset_name, &numReturnDatasets,
			   &returnDatasets);
  fail_unless(rc == FC_SUCCESS, "aborted - cant get dataset by name");
  fail_unless(numReturnDatasets == 1,
	      "failed to find unique dataset by name");
  dataset = returnDatasets[0];
  free(returnDatasets);

  // create a variable to play with
  rc = fc_createGlobalVariable(dataset, name, &variable);
  fail_unless(rc == FC_SUCCESS, "failed to create variable");
  
  // test fc_isVariableValid()
  fail_unless(fc_isVariableValid(variable),
	      "failed to validate a valid variable");

  // test fc_isVariableGlobal()
  fail_unless(fc_isVariableGlobal(variable), "should be global");

  // test fc_getDatasetFromVariable()
  temp_dataset = FC_NULL_DATASET;
  rc = fc_getDatasetFromVariable(variable, &temp_dataset);
  fail_unless(rc == FC_SUCCESS, "failed to get parent dataset");
  fail_unless(FC_HANDLE_EQUIV(temp_dataset, dataset), "mismatch of parent dataset"); 
 
  // test fc_getMeshFromVariable()
  temp_mesh = badMesh;
  rc = fc_getMeshFromVariable(variable, &temp_mesh);
  fail_unless(rc == FC_SUCCESS, "failed to get parent mesh");
  fail_unless(FC_HANDLE_EQUIV(temp_mesh, FC_NULL_MESH), "mismatch of parent mesh");
  
  // --- all done
  
  // one last test
  fc_deleteVariable(variable);
  fail_unless(!fc_isVariableValid(variable), 
	      "handle should not be valid after a delete");
}
END_TEST

// query the variable data - all assocs, mathtypes & datatypes
START_TEST(glb_data_query)
{
  FC_ReturnCode rc;
  int k, m;
  char name[20] = "blue berry", name2[20] = "banana nut";
  FC_Dataset dataset, *returnDatasets;
  int numReturnDatasets;
  FC_Variable variable, variable2;
  int numDataPoint = 1; 
  int numMathType = 4, numComps[4] = { 1, 2, 3, 4 }, maxNumComp = 4;
  FC_AssociationType assoc = FC_AT_WHOLE_DATASET;
  FC_MathType mathTypes[4] = { FC_MT_SCALAR, FC_MT_VECTOR, 
			       FC_MT_SYMTENSOR, FC_MT_TENSOR };
  int numDataType = 4;
  FC_DataType dataTypes[4] = { FC_DT_CHAR, FC_DT_INT, 
			       FC_DT_FLOAT, FC_DT_DOUBLE };
  char charData[4] = { 'a', 'b', 'c', 'd' };
  int intData[4];
  float floatData[4];
  double doubleData[4];
  int sizes[4] = {sizeof(char), sizeof(int), sizeof(float), sizeof(double) };
  void *data[4] = { charData, intData, floatData, doubleData };
  int temp_numDataPoint, temp_numComp;
  FC_AssociationType temp_assoc;
  FC_MathType temp_mathType;
  FC_DataType temp_dataType;
  void *temp_data;
  
  // get dataset
  rc = fc_getDatasetByName(dataset_name, &numReturnDatasets,
			   &returnDatasets);
  fail_unless(rc == FC_SUCCESS, "aborted - cant get dataset by name");
  fail_unless(numReturnDatasets == 1,
	      "failed to find unique dataset by name");
  dataset = returnDatasets[0];
  free(returnDatasets);

  // setup some data
  for (k = 0; k < maxNumComp; k++) {
    intData[k] = k;
    floatData[k] = intData[k] + 0.1;
    doubleData[k] = intData[k] + 0.00000000001;
  }

  // loop over all math types & data types
  for (k = 0; k < numMathType; k++) {
    for (m = 0; m < numDataType; m++) {
      
      // create two variables
      rc = fc_createGlobalVariable(dataset, name, &variable);
      fail_unless(rc == FC_SUCCESS, "abort: failed to create variable");
      rc = fc_createGlobalVariable(dataset, name2, &variable2);
      fail_unless(rc == FC_SUCCESS, "abort: failed to create variable");
	
      // check default values
      rc = fc_getVariableInfo(variable, &temp_numDataPoint, &temp_numComp, 
			      &temp_assoc, &temp_mathType, &temp_dataType);
      fail_unless(rc == FC_SUCCESS, "failed to get variable info");
      fail_unless(temp_numDataPoint == 0 && temp_numComp == 0 && 
		  temp_assoc == FC_AT_UNKNOWN && temp_mathType ==
		  FC_MT_UNKNOWN && temp_dataType == FC_DT_UNKNOWN,
		  "empty variable should have empty values");
      
      // variable1: set/get data by copy and check
      temp_numDataPoint = 999;
      temp_dataType = FC_DT_UNKNOWN;
      temp_data = NULL;
      rc = fc_setVariableData(variable, numDataPoint, numComps[k],assoc,
			      mathTypes[k], dataTypes[m], data[m]);
      fail_unless(rc == FC_SUCCESS, "failed to set data");
      rc = fc_getVariableInfo(variable, &temp_numDataPoint, &temp_numComp, 
			      &temp_assoc, &temp_mathType, &temp_dataType);
      fail_unless(rc == FC_SUCCESS, "failed to get variable info");
      fail_unless(temp_numDataPoint == numDataPoint, "mismatch of numDataPoint");
      fail_unless(temp_numComp == numComps[k], "mismatch of numComp");
      fail_unless(temp_assoc == assoc, "mismatch of assoc");
      fail_unless(temp_mathType == mathTypes[k], "mismatch of mathType");
      fail_unless(temp_dataType == dataTypes[m], "mismatch of datatype");
      rc = fc_getVariableDataPtr(variable, &temp_data);
      fail_unless(rc == FC_SUCCESS, "failed to get data");
      fail_unless(!memcmp(temp_data, data[m], 
			  numDataPoint*numComps[k]*sizes[m]),
		  "mismatch of data");
      
      
      // ---- check errors
      
      // set with assoc = vertex should fail
      rc = fc_setVariableData(variable2, numDataPoint, numComps[k],
			      FC_AT_VERTEX, mathTypes[k], dataTypes[m],
			      data[m]);
      fail_unless(rc < FC_SUCCESS, "should fail to set data");
      
      // --- all done

      // cleanup
      fc_deleteVariable(variable);
    }
  }
}
END_TEST

// copy
START_TEST(glb_copy_test)
{
  FC_ReturnCode rc;
  int k, m;
  char name[20] = "blue berry", copy_name[20] = "banana nut", *temp_name;
  FC_Dataset dataset1, dataset2, *returnDatasets;
  FC_Mesh mesh, *returnMeshes;
  int numReturnDatasets, numReturnMeshes, numElem;
  int numVariable1, numVariable2, temp_numVariable;
  FC_Variable variable, copy_variable, badVariable = { 999, 999 };
  FC_Variable wholeMeshVar, elemVar;
  int numDataPoint = 1; 
  FC_AssociationType assoc = FC_AT_WHOLE_DATASET;
  int numMathType = 4, numComps[4] = { 1, 2, 3, 4 }, maxNumComp = 4;
  FC_MathType mathTypes[4] = { FC_MT_SCALAR, FC_MT_VECTOR, 
			       FC_MT_SYMTENSOR, FC_MT_TENSOR };
  int numDataType = 4;
  FC_DataType dataTypes[4] = { FC_DT_CHAR, FC_DT_INT, 
			       FC_DT_FLOAT, FC_DT_DOUBLE };
  char charData[4] = { 'a', 'b', 'c', 'd' };
  int intData[4];
  float floatData[4];
  double doubleData[4];
  int sizes[4] = {sizeof(char), sizeof(int), sizeof(float), sizeof(double) };
  void *data[4] = { charData, intData, floatData, doubleData };
  int temp_numDataPoint, temp_numComp;
  FC_AssociationType temp_assoc;
  FC_MathType temp_mathType;
  FC_DataType temp_dataType;
  void *temp_data;
  
  // get dataset and mesh
  rc = fc_getDatasetByName(dataset_name, &numReturnDatasets,
			   &returnDatasets);
  fail_unless(rc == FC_SUCCESS, "aborted - cant get dataset by name");
  fail_unless(numReturnDatasets == 1,
	      "failed to find unique dataset by name");
  dataset1 = returnDatasets[0];
  free(returnDatasets);
  rc = fc_getMeshByName(dataset1, "hex mesh", &numReturnMeshes,
			&returnMeshes);
  fail_unless(rc == FC_SUCCESS, "aborted - cant get mesh by name");
  fail_unless(numReturnMeshes == 1, "failed to find unique mesh by name");
  mesh = returnMeshes[0];
  free(returnMeshes);

  // make a second dataset to copy to
  rc = fc_copyDataset(dataset1, "copy dataset", &dataset2);
  fail_unless(rc == FC_SUCCESS, "abort: failed to copy dataset");

  // setup some data
  for (k = 0; k < maxNumComp; k++) {
    intData[k] = k;
    floatData[k] = intData[k] + 0.1;
    doubleData[k] = intData[k] + 0.00000000001;
  }

  // Make a whole mesh variable (double)
  rc = fc_createVariable(mesh, name, &wholeMeshVar);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create whole mesh var");
  rc = fc_setVariableData(wholeMeshVar, 1, 1, FC_AT_WHOLE_MESH, FC_MT_SCALAR,
                          FC_DT_DOUBLE, doubleData);
  fail_unless(rc == FC_SUCCESS, "abort: failed to set whole mesh var data");
  // Make an element variable (double) 
  rc = fc_createVariable(mesh, name, &elemVar);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create vert var");
  fc_getMeshNumElement(mesh, &numElem);
  fail_unless(numElem < maxNumComp, "abort: developer assumption changed");
  rc = fc_setVariableData(elemVar, numElem, 1, FC_AT_ELEMENT, FC_MT_SCALAR,
                          FC_DT_DOUBLE, doubleData);
  fail_unless(rc == FC_SUCCESS, "abort: failed to set vert var data");

  numVariable1 = 0;
  numVariable2 = 0;
  // math types & data types
  for (k = 0; k < numMathType; k++) {
    for (m = 0; m < numDataType; m++) {
      
      // create a variable
      rc = fc_createGlobalVariable(dataset1, name, &variable);
      fail_unless(rc == FC_SUCCESS, "abort: failed to create variable");
      rc = fc_setVariableData(variable, numDataPoint, numComps[k], assoc,
			      mathTypes[k], dataTypes[m], data[m]);
      numVariable1++;

      // copy to same dataset & check it
      rc = fc_copyGlobalVariable(variable, dataset1, copy_name, &copy_variable);
      fail_unless(rc == FC_SUCCESS, "failed to copy to same dataset");
      numVariable1++;
      rc = fc_getVariableName(copy_variable, &temp_name);
      fail_unless(rc == FC_SUCCESS, "failed to get name of copy");
      fail_unless(!strcmp(temp_name, copy_name), "mismatch of name");
      free(temp_name);
      rc = fc_getVariableInfo(copy_variable, &temp_numDataPoint, 
			      &temp_numComp, &temp_assoc, &temp_mathType, 
			      &temp_dataType);
      fail_unless(rc == FC_SUCCESS, "failed to get variable info");
      fail_unless(temp_numDataPoint == numDataPoint, 
		  "mismatch of numDataPoint");
      fail_unless(temp_numComp == numComps[k], "mismatch of numComp");
      fail_unless(temp_assoc == assoc, "mismatch of assoc");
      fail_unless(temp_mathType == mathTypes[k], "mismatch of mathType");
      fail_unless(temp_dataType == dataTypes[m], "mismatch of datatype");
      rc = fc_getVariableDataPtr(copy_variable, &temp_data);
      fail_unless(rc == FC_SUCCESS, "failed to get data");
      fail_unless(!memcmp(temp_data, data[m], 
			  numDataPoint*numComps[k]*sizes[m]),
		    "mismatch of data");

      // copy to different dataset & check it
      rc = fc_copyGlobalVariable(variable, dataset2, copy_name, &copy_variable);
      fail_unless(rc == FC_SUCCESS, "failed to copy to different dataset");
      numVariable2++;
      rc = fc_getVariableName(copy_variable, &temp_name);
      fail_unless(rc == FC_SUCCESS, "failed to get name of copy");
      fail_unless(!strcmp(temp_name, copy_name), "mismatch of name");
      free(temp_name);
      rc = fc_getVariableInfo(copy_variable, &temp_numDataPoint, 
			      &temp_numComp, &temp_assoc, &temp_mathType, 
			      &temp_dataType);
      fail_unless(rc == FC_SUCCESS, "failed to get variable info");
      fail_unless(temp_numDataPoint == numDataPoint, 
		  "mismatch of numDataPoint");
      fail_unless(temp_numComp == numComps[k], "mismatch of numComp");
      fail_unless(temp_assoc == assoc, "mismatch of assoc");
      fail_unless(temp_mathType == mathTypes[k], "mismatch of mathType");
      fail_unless(temp_dataType == dataTypes[m], "mismatch of datatype");
      rc = fc_getVariableDataPtr(copy_variable, &temp_data);
      fail_unless(rc == FC_SUCCESS, "failed to get data");
      fail_unless(!memcmp(temp_data, data[m], 
			    numDataPoint*numComps[k]*sizes[m]),
		  "mismatch of data");
      
      // --- special cases
      
      // copy with NULL name will use original's name
      rc = fc_copyGlobalVariable(variable, dataset1, NULL, &copy_variable);
      fail_unless(rc == FC_SUCCESS, "should not fail");
      rc = fc_getVariableName(copy_variable, &temp_name);
      fail_unless(rc == FC_SUCCESS, "failed to get name of copy");
      fail_unless(!strcmp(temp_name, name), "mismatch of name");
      free(temp_name);
      rc = fc_deleteVariable(copy_variable);
      fail_unless(rc == FC_SUCCESS, "failed to delete copied variable");
      
      // ---- check errors
      
      // fc_copyVariable() -- bad args
      rc = fc_copyGlobalVariable(badVariable, dataset1, copy_name, &copy_variable);
      fail_unless(rc != FC_SUCCESS, "should fail to copy bad variable");
      fail_unless(FC_HANDLE_EQUIV(copy_variable, FC_NULL_VARIABLE),
		  "fail should return NULL");
      rc = fc_copyGlobalVariable(variable, FC_NULL_DATASET, copy_name, 
			   &copy_variable);
      fail_unless(rc != FC_SUCCESS, "should fail to copy to bad dataset");
      fail_unless(FC_HANDLE_EQUIV(copy_variable, FC_NULL_VARIABLE),
		  "fail should return NULL");
      rc = fc_copyGlobalVariable(variable, dataset1, copy_name, NULL);
	fail_unless(rc != FC_SUCCESS, "should fail to copy if NULL handle");
    }
  }
  
  // It's o.k. to copy a whole mesh var onto a dataset
  rc = fc_copyGlobalVariable(wholeMeshVar, dataset1, copy_name, &copy_variable);
  fail_unless(rc == FC_SUCCESS, "failed to copy whole mesh var to a dataset");
  numVariable1++;
  rc = fc_getVariableName(copy_variable, &temp_name);
  fail_unless(rc == FC_SUCCESS, "failed to get name of copy");
  fail_unless(!strcmp(temp_name, copy_name), "mismatch of name");
  free(temp_name);
  rc = fc_getVariableInfo(copy_variable, &temp_numDataPoint, 
                          &temp_numComp, &temp_assoc, &temp_mathType, 
                          &temp_dataType);
  fail_unless(rc == FC_SUCCESS, "failed to get variable info");
  fail_unless(temp_numDataPoint == 1, "mismatch of numDataPoint");
  fail_unless(temp_numComp == 1, "mismatch of numComp");
  fail_unless(temp_assoc == FC_AT_WHOLE_DATASET, "mismatch of assoc");
  fail_unless(temp_mathType == FC_MT_SCALAR, "mismatch of mathType");
  fail_unless(temp_dataType == FC_DT_DOUBLE, "mismatch of datatype");
  rc = fc_getVariableDataPtr(copy_variable, &temp_data);
  fail_unless(rc == FC_SUCCESS, "failed to get data");
  fail_unless(!memcmp(temp_data, doubleData, sizeof(double)), 
              "mismatch of data");

  // It's NOT o.k. to copy an elem var onto a dataset
  rc = fc_copyGlobalVariable(elemVar, dataset1, copy_name, &copy_variable);
  fail_unless(rc != FC_SUCCESS, "should fail to copy elem var to dataset");
  fail_unless(FC_HANDLE_EQUIV(copy_variable, FC_NULL_VARIABLE),
              "fail should return NULL");

  // a little more checking
  rc = fc_getNumGlobalVariable(dataset2, &temp_numVariable);
  fail_unless(rc == FC_SUCCESS, "shouldn't fail to get num variables");
  fail_unless(temp_numVariable == numVariable2, "mismatch of numVariable");
  rc = fc_deleteDataset(dataset2);
  fail_unless(rc == FC_SUCCESS, "failed to delete second dataset");
  rc = fc_getNumGlobalVariable(dataset1, &temp_numVariable);
  fail_unless(rc == FC_SUCCESS, "shouldn't fail to get num variable");
  fail_unless(temp_numVariable == numVariable1, "mismatch of numVariable");
}
END_TEST

// --- global seq variable interface

// only stuff specific to global seq variables is tested here since
// the basic variable & seq variable stuff is tested in variable interace

// test creating and deleting global seq variables, also test querying 
// for global seq vars before, inbetween and after.
START_TEST(glb_seq_create_get_del)
{
  FC_ReturnCode rc;
  int j, k, m;
  FC_Dataset dataset, *returnDatasets;
  int numReturnDatasets;
  FC_Sequence *sequences;
  int numSequence;
  int numSeqVar = 10, temp_numSeqVar, numSteps[2];
  int temp_numStep, *temp_numSteps;
  FC_Variable *seqVars[2*10], **temp_seqVars, *temp_seqVar;
  FC_Variable *badSeqVar, badVariable = { 999, 999 };
  char names[2][10][100] = { { "one", "two", "three", "four", "five", 
			       "six", "seven", "eight", "nine", "ten" },
			     { "un", "deux", "trois", "quatre", "cinq",
			       "six (fr.)", "sept", "huit", "neuf", "dix"}};
  char newName[100] = { "sparkly new" };

  _FC_VarSlot *varSlot, *varSlot2;

  // get dataset, sequences 
  rc = fc_getDatasetByName(dataset_name, &numReturnDatasets,
			   &returnDatasets);
  fail_unless(rc == FC_SUCCESS, "aborted - cant get dataset by name");
  fail_unless(numReturnDatasets == 1,
	      "failed to find unique dataset by name");
  dataset = returnDatasets[0];
  free(returnDatasets);
  rc = fc_getSequences(dataset, &numSequence, &sequences);
  fail_unless(rc == FC_SUCCESS, "abort: failed to get sequences");
  fail_unless(numSequence == 2, "abort: expected 2 sequences");
  for (j = 0; j < numSequence; j++) 
    fc_getSequenceNumStep(sequences[j], &numSteps[j]);

  // create some seq variables for each sequence
  for (j = 0; j < numSequence; j++) {
    for (k = 0; k < numSeqVar; k++) {
      int idx = j*numSeqVar + k;
      rc = fc_createGlobalSeqVariable(dataset, sequences[j], names[j][k], 
				      &temp_numStep, &seqVars[idx]);
      fail_unless(rc == FC_SUCCESS, "failed to create seqVar");
      fail_unless(temp_numStep == numSteps[j], "mismatch in numStep");
      for (m = 0; m < numSteps[j]; m++)
	fail_unless(!FC_HANDLE_EQUIV(seqVars[idx][m], FC_NULL_VARIABLE),
		      "created variable should not = FC_NULL_VARIABLE");
    }
  }

  // get seq variables (should be numSeqVar*numSequence);
  // (will come out in same order they were written)
  rc = fc_getGlobalSeqVariables(dataset, &temp_numSeqVar, &temp_numSteps, 
				&temp_seqVars);
  fail_unless(rc == FC_SUCCESS, "failed to get seq variables");
  fail_unless(temp_numSeqVar == numSeqVar*numSequence, 
	      "mismatch of numSeqVariable");
  for (j = 0; j < numSequence; j++) {
    for (k = 0; k < numSeqVar; k++) {
      int idx = j*numSeqVar + k;
      varSlot = _fc_getVarSlot(temp_seqVars[idx][0]);
      fail_unless(FC_HANDLE_EQUIV(varSlot->sequence, sequences[j]),
		  "mismatch of sequence handle");
      fail_unless(temp_numSteps[idx] == numSteps[j],
		  "mismatch of numStep");
      for (m = 0; m < temp_numSteps[idx]; m++)
	fail_unless(FC_HANDLE_EQUIV(temp_seqVars[idx][m], 
				    seqVars[idx][m]), 
		    "mismatch of variable handles");
      free(temp_seqVars[idx]);
    }
  }
  free(temp_numSteps);
  free(temp_seqVars);

  rc = fc_getNumGlobalSeqVariable(dataset, &temp_numSeqVar);
  fail_unless(rc == FC_SUCCESS, "failed to get num seq variables");
  fail_unless(temp_numSeqVar == numSeqVar*numSequence, 
	      "mismatch of numSeqVariable");

  // test fc_getSeqVariableByName()
  for (j = 0; j < numSequence; j++) {
    for (k = 0; k < numSeqVar; k++) {
      int idx = j*numSeqVar + k;
      rc = fc_getGlobalSeqVariableByName(dataset, names[j][k],
					 &temp_numSeqVar,
					 &temp_numSteps,
					 &temp_seqVars);
      fail_unless(rc == FC_SUCCESS, "failed to get variable by name");
      fail_unless(temp_numSeqVar == 1, "wrong number of matching seq vars");
      varSlot = _fc_getVarSlot(temp_seqVars[0][0]);
      fail_unless(FC_HANDLE_EQUIV(varSlot->sequence, sequences[j]),
		  "mismatch of sequence handle");
      fail_unless(temp_numSteps[0] == numSteps[j], "mismatch of numStep");
      for (m = 0; m < temp_numSteps[0]; m++)
	fail_unless(FC_HANDLE_EQUIV(temp_seqVars[0][m], seqVars[idx][m]),
		    "mismatch of variable handles");
      free(temp_numSteps);
      free(temp_seqVars[0]);
      free(temp_seqVars);
    }
  }

  // test negative example of fc_getGlobalSeqVariableByName
  rc = fc_getGlobalSeqVariableByName(dataset, newName, 
				     &temp_numSeqVar,
				     &temp_numSteps,
				     &temp_seqVars);
  fail_unless(rc == FC_SUCCESS, "new name should succeed");
  fail_unless(temp_numSeqVar == 0 && temp_numSteps == NULL &&
	      temp_seqVars == NULL, 
	      "new name shouldn't find anything and return NULL");


  //temporarily create a new seq var to check array return
  //make one the same as the 0th one, but on the first sequence
  rc = fc_createGlobalSeqVariable(dataset, sequences[1], names[0][0], 
				      &temp_numStep, &temp_seqVar);
  fail_unless(rc == FC_SUCCESS, "failed to create seqVar");
  fail_unless(temp_numStep == numSteps[1], "mismatch in numStep");
  for (m = 0; m < numSteps[1]; m++)
    fail_unless(!FC_HANDLE_EQUIV(temp_seqVar[m], FC_NULL_VARIABLE),
		"created variable should not = FC_NULL_VARIABLE");

  rc = fc_getGlobalSeqVariableByName(dataset, names[0][0],
				     &temp_numSeqVar,
				     &temp_numSteps,
				     &temp_seqVars);
  fail_unless(rc == FC_SUCCESS, "failed to get variable by name");
  fail_unless(temp_numSeqVar == 2, "wrong number of matching seq vars");
  varSlot = _fc_getVarSlot(temp_seqVars[0][0]);
  varSlot2 = _fc_getVarSlot(temp_seqVars[1][0]);
  fail_unless(((FC_HANDLE_EQUIV(varSlot->sequence, sequences[0]) &&
		FC_HANDLE_EQUIV(varSlot2->sequence, sequences[1])) ||
	       (FC_HANDLE_EQUIV(varSlot->sequence, sequences[1]) &&
		FC_HANDLE_EQUIV(varSlot2->sequence, sequences[0]))),
	      "mismatch of sequence handle");
  j = !(FC_HANDLE_EQUIV(varSlot->sequence, sequences[0]));
  fail_unless(temp_numSteps[0] == numSteps[j] &&
	      temp_numSteps[1] == numSteps[!j], "mismatch of numStep");
  for (m = 0; m < temp_numSteps[j]; m++){
    fail_unless(FC_HANDLE_EQUIV(temp_seqVars[j][m], seqVars[0][m]),
		"mismatch of variable handles");
  }
  for (m = 0; m < temp_numSteps[!j]; m++){
    fail_unless(FC_HANDLE_EQUIV(temp_seqVars[!j][m], temp_seqVar[m]),
		"mismatch of variable handles");
  }

  fc_deleteSeqVariable(temp_numStep,temp_seqVar);
  free(temp_seqVar);

  for (m = 0; m < temp_numSeqVar; m++){
      free(temp_seqVars[m]);
  }
  free(temp_seqVars);
  free(temp_numSteps);

  
  // delete half of the seq variables (alternate)
  for (j = 0; j < numSequence; j++) {
    for (k = 0; k < numSeqVar; k+=2) {
      int idx = j*numSeqVar + k;
      rc = fc_deleteSeqVariable(numSteps[j], seqVars[idx]);
      fail_unless(rc == FC_SUCCESS, "failed to delete seq variable");
      free(seqVars[idx]);
    }
  }

  // get global seq variables (should be numSeqVar*numSequence/2);
  fc_getGlobalSeqVariables(dataset, &temp_numSeqVar, &temp_numSteps, 
			   &temp_seqVars);
  fail_unless(rc == FC_SUCCESS, "failed to get seq variables");
  fail_unless(temp_numSeqVar == numSeqVar*numSequence/2, 
	      "mismatch of numSeqVariable");
  for (j = 0; j < numSequence; j++) {
    for (k = 0; k < numSeqVar/2; k++) {
      int idx = j*numSeqVar/2 + k; // into temp array
      int idx_orig = j*numSeqVar + 2*k + 1; // into original array
      varSlot = _fc_getVarSlot(temp_seqVars[idx][0]);
      fail_unless(FC_HANDLE_EQUIV(varSlot->sequence, sequences[j]),
		  "mismatch of sequence handle");
      fail_unless(temp_numSteps[idx] == numSteps[j],
		  "mismatch of numStep");
      for (m = 0; m < temp_numSteps[idx]; m++)
	fail_unless(FC_HANDLE_EQUIV(temp_seqVars[idx][m],
				    seqVars[idx_orig][m]), 
		    "mismatch of variable handles");
      free(temp_seqVars[idx]);
    }
  }
  free(temp_numSteps);
  free(temp_seqVars);

  // delete remaining seq variables on first sequence
  for (k = 1; k < numSeqVar; k+=2) {
    rc = fc_deleteSeqVariable(numSteps[0], seqVars[k]);
    fail_unless(rc == FC_SUCCESS, "failed to delete variable");
    free(seqVars[k]);
  }

  // get seq variables (should be numSeqVar*(numSequence-1)/2);
  fc_getGlobalSeqVariables(dataset, &temp_numSeqVar, &temp_numSteps, 
			   &temp_seqVars);
  fail_unless(rc == FC_SUCCESS, "failed to get seq variables");
  fail_unless(temp_numSeqVar == numSeqVar*(numSequence-1)/2, 
	      "mismatch of numSeqVariable");
  for (j = 1; j < numSequence; j++) {
    for (k = 0; k < numSeqVar/2; k++) {
      int idx = (j-1)*numSeqVar/2 + k; // into temp array
      int idx_orig = j*numSeqVar + 2*k + 1; // into original array
      varSlot = _fc_getVarSlot(temp_seqVars[idx][0]);
      fail_unless(FC_HANDLE_EQUIV(varSlot->sequence, sequences[j]),
		  "mismatch of sequence handle");
      fail_unless(temp_numSteps[idx] == numSteps[j],
		    "mismatch of numStep");
      for (m = 0; m < temp_numSteps[idx]; m++)
	fail_unless(FC_HANDLE_EQUIV(temp_seqVars[idx][m],
				    seqVars[idx_orig][m]), 
		      "mismatch of variable handles");
      free(temp_seqVars[idx]);
    }
  }
  free(temp_numSteps);
  free(temp_seqVars);
  
  // delete remaining variables
  for (j = 1; j < numSequence; j++) {
    for (k = 1; k < numSeqVar; k+=2) { 
      rc = fc_deleteSeqVariable(numSteps[j], seqVars[j*numSeqVar+k]);
      fail_unless(rc == FC_SUCCESS, "failed to delete variable");
      free(seqVars[j*numSeqVar+k]);
    }
  }
  
  // get seq variables (should be none);
  rc = fc_getGlobalSeqVariables(dataset, &temp_numSeqVar, &temp_numSteps, 
				&temp_seqVars);
  fail_unless(rc == FC_SUCCESS, "failed to get seq variables");
  fail_unless(temp_numSeqVar == 0 && temp_numSteps == NULL && 
	      temp_seqVars == NULL, 
		"should return 0 & nulls if all seqVars delete");
  
  // create 1 seq var for further testing
  rc = fc_createGlobalSeqVariable(dataset, sequences[0], names[0][0],
				  &temp_numStep, &seqVars[0]);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create seq variable");

  // ---- test error conditions
    
  // construct a bad seqVar (make one of the vars bad)
  badSeqVar = (FC_Variable*)malloc(numSteps[0]*sizeof(FC_Variable));
  for (j = 0; j < numSteps[0]; j++)
    badSeqVar[j] = seqVars[0][j];
  badSeqVar[numSteps[0]/2] = badVariable;
  
  // bad args to fc_createGlobalSeqVariable()
  rc = fc_createGlobalSeqVariable(FC_NULL_DATASET, sequences[0], names[0][0],
				  &temp_numStep, &temp_seqVar);
  fail_unless(rc != FC_SUCCESS, "should fail with bad dataset");
  fail_unless(temp_numStep == -1 && temp_seqVar == NULL, 
	      "fail should return NULL seqVar");
  rc = fc_createGlobalSeqVariable(dataset, FC_NULL_SEQUENCE, names[0][0], 
			    &temp_numStep, &temp_seqVar);
  fail_unless(rc != FC_SUCCESS, "should fail with bad sequence");
  fail_unless(temp_numStep == -1 && temp_seqVar == NULL, 
	      "fail should return NULL seqVar");
  rc = fc_createGlobalSeqVariable(dataset, sequences[0], NULL, 
			    &temp_numStep, &temp_seqVar);
  fail_unless(rc != FC_SUCCESS, "should fail with null name");
  fail_unless(temp_numStep == -1 && temp_seqVar == NULL, 
	      "fail should return NULL seqVar");
  rc = fc_createGlobalSeqVariable(dataset, sequences[0], names[0][0], 
			    NULL, &temp_seqVar);
  fail_unless(rc != FC_SUCCESS, "should fail with null numStep");
  fail_unless(temp_seqVar == NULL, "fail should return NULL seqVar");
  rc = fc_createGlobalSeqVariable(dataset, sequences[0], names[0][0], 
			    &temp_numStep, NULL);
  fail_unless(rc != FC_SUCCESS, "should fail with null seqVar");
  fail_unless(temp_numStep == -1, "fail should return NULL");
    
  // bad args to fc_getGlobalSeqVariables()
  temp_numSteps = (int*)1;
  temp_seqVars = (FC_Variable**)1;
  rc = fc_getGlobalSeqVariables(FC_NULL_DATASET, &temp_numSeqVar,
				&temp_numSteps, &temp_seqVars);
  fail_unless(rc != FC_SUCCESS, 
		"should have failed to get variables with null dataset");
  fail_unless(temp_numSeqVar == -1 && temp_numSteps == NULL && 
	      temp_seqVars == NULL, "fail should return nulls"); 
  temp_numSteps = (int*)1;
  temp_seqVars = (FC_Variable**)1;
  rc = fc_getGlobalSeqVariables(dataset, NULL, &temp_numSteps, &temp_seqVars);
  fail_unless(rc != FC_SUCCESS, 
	      "should have failed with NULL numSeqVar");
  fail_unless(temp_numSteps == NULL && temp_seqVars == NULL, 
	      "fail should return nulls"); 
  temp_seqVars = (FC_Variable**)1;
  rc = fc_getGlobalSeqVariables(dataset, &temp_numSeqVar, NULL, &temp_seqVars);
  fail_unless(rc != FC_SUCCESS, 
	      "should have failed with NULL numStepPerVar");
  fail_unless(temp_numSeqVar == -1 && temp_seqVars == NULL, 
	      "fail should return nulls"); 
  temp_numSteps = (int*)1;
  rc = fc_getGlobalSeqVariables(dataset, &temp_numSeqVar, &temp_numSteps,
				NULL);
  fail_unless(rc != FC_SUCCESS, 
	      "should have failed with NULL seqVars");
  fail_unless(temp_numSeqVar == -1 && temp_numSteps == NULL, 
	      "fail should return nulls"); 
 
  // bad args to fc_getNumGlobalSeqVariable()
  rc = fc_getNumGlobalSeqVariable(FC_NULL_DATASET, &temp_numSeqVar);
  fail_unless(rc != FC_SUCCESS, 
	      "should have failed to get num seq variables with null dataset");
  fail_unless(temp_numSeqVar == -1, "fail should return nulls"); 
  rc = fc_getNumGlobalSeqVariable(dataset, NULL);
  fail_unless(rc != FC_SUCCESS, 
	      "should have failed with NULL numSeqVar");

  // bad args to fc_getGlobalSeqVarByName()
  rc = fc_getGlobalSeqVariableByName(FC_NULL_DATASET, names[0][0],
				     &temp_numSeqVar,
				     &temp_numSteps,
				     &temp_seqVars);
  fail_unless(rc != FC_SUCCESS, "should fail with NULL dataset");
  fail_unless(temp_numSeqVar == -1 && temp_numSteps == NULL &&
	      temp_seqVars == NULL, "should return -1 and NULL when fail");

  rc = fc_getGlobalSeqVariableByName(dataset, NULL,
				     &temp_numSeqVar,
				     &temp_numSteps,
				     &temp_seqVars);
  fail_unless(rc != FC_SUCCESS, "should fail with NULL name");
  fail_unless(temp_numSeqVar == -1 && temp_numSteps == NULL &&
	      temp_seqVars == NULL, "should return -1 and NULL when fail");

  rc = fc_getGlobalSeqVariableByName(dataset, names[0][0],
				     NULL,
				     &temp_numSteps,
				     &temp_seqVars);
  fail_unless(rc != FC_SUCCESS, "should fail with NULL arg");
  fail_unless(temp_numSteps == NULL &&
	      temp_seqVars == NULL, "should return -1 and NULL when fail");

  rc = fc_getGlobalSeqVariableByName(dataset, names[0][0],
				     &temp_numSeqVar,
				     NULL,
				     &temp_seqVars);
  fail_unless(rc != FC_SUCCESS, "should fail with NULL arg");
  fail_unless(temp_numSeqVar == -1 && temp_seqVars == NULL,
	      "should return -1 and NULL when fail");

  rc = fc_getGlobalSeqVariableByName(dataset, names[0][0],
				     &temp_numSeqVar,
				     &temp_numSteps,
				     NULL);
  fail_unless(rc != FC_SUCCESS, "should fail with NULL name");
  fail_unless(temp_numSeqVar == -1 && temp_numSteps == NULL,
	      "should return -1 and NULL when fail");


  // --- done
  
  // delete last sequence variable
  rc = fc_deleteSeqVariable(numSteps[0], seqVars[0]);
  fail_unless(rc == FC_SUCCESS, "failed to delete last seq var");
  
  // cleanup
  free(seqVars[0]);
  free(badSeqVar);
  free(sequences);
}
END_TEST

// query meta data
START_TEST(glb_seq_metadata_query)
{
  FC_ReturnCode rc;
  int j, k;
  FC_Dataset dataset, temp_dataset, badDataset = { 999, 999 };
  FC_Dataset *returnDatasets;
  int numReturnDatasets;
  int numSequence;
  FC_Sequence *sequences, temp_sequence, badSequence = { 999, 999 };
  FC_Mesh temp_mesh, badMesh = { 999, 999 };
  int numStep;
  char name[20] = "blue berry";
  FC_Variable *seqVar;

  // get dataset, sequences
  rc = fc_getDatasetByName(dataset_name, &numReturnDatasets,
			   &returnDatasets);
  fail_unless(rc == FC_SUCCESS, "aborted - cant get dataset by name");
  fail_unless(numReturnDatasets == 1,
	      "failed to find unique dataset by name");
  dataset = returnDatasets[0];
  free(returnDatasets);
  rc = fc_getSequences(dataset, &numSequence, &sequences);
  fail_unless(rc == FC_SUCCESS, "abort: failed to get sequences");
  fail_unless(numSequence == 2, "abort: expected 2 sequences");

  // loop over  sequences
  for (j = 0; j < numSequence; j++) {
      
    // create a seq variable to play with
    rc = fc_createGlobalSeqVariable(dataset, sequences[j], name, &numStep, 
				    &seqVar);
    fail_unless(rc == FC_SUCCESS, "failed to create seqVar");
    
    // test fc_isSeqVariableValid()
    fail_unless(fc_isSeqVariableValid(numStep, seqVar), 
		"failed to validate a valid seq variable");

    // test fc_isVariableGlobal()
    fail_unless(fc_isVariableGlobal(seqVar[0]), "should be global");
    
    // test fc_getDatasetFromVariable()
    for (k = 0; k < numStep; k++) {
      temp_dataset = badDataset;
      rc = fc_getDatasetFromVariable(seqVar[k], &temp_dataset);
      fail_unless(rc == FC_SUCCESS, "failed to get parent dataset");
	fail_unless(FC_HANDLE_EQUIV(temp_dataset, dataset), 
		    "mismatch of parent dataset");
    }
    
    // test fc_getMeshFromVariable()
    for (k = 0; k < numStep; k++) {
      temp_mesh = badMesh;
      rc = fc_getMeshFromVariable(seqVar[k], &temp_mesh);
      fail_unless(rc == FC_SUCCESS, "failed to get parent mesh");
      fail_unless(FC_HANDLE_EQUIV(temp_mesh, FC_NULL_MESH), 
		  "mismatch of parent mesh");
    }

    // test fc_getSequenceFromSeqVariable()
    temp_sequence = badSequence;
    rc = fc_getSequenceFromSeqVariable(numStep, seqVar,
				       &temp_sequence);
    fail_unless(rc == FC_SUCCESS, "failed to get parent sequence");
    fail_unless(FC_HANDLE_EQUIV(temp_sequence, sequences[j]), 
		"mismatch of parent sequence");

    // --- all done
    
    // one last test
    fc_deleteSeqVariable(numStep, seqVar);
    fail_unless(!fc_isSeqVariableValid(numStep, seqVar),
		"seq var should not be valid after a delete");
    
    // cleanup
    free(seqVar);
  }

  free(sequences);
}
END_TEST

// query the seq variable data by step slice & point slice
// don't bother testing set data/get data interface except where it
// being a sequence variable matters (since fully tested above)
START_TEST(glb_seq_data_query)
{
  FC_ReturnCode rc;
  int k, m, n;
  char name[20] = "blue berry", name2[20] = "banana nut", *temp_name;
  FC_Dataset dataset, *returnDatasets;
  int numReturnDatasets;
  FC_Sequence sequence, *returnSequences;
  int numReturnSequences;
  int numStep;
  FC_Variable *seqVar, *seqVar2;
  int numDataPoint = 1; 
  FC_AssociationType assoc = FC_AT_WHOLE_DATASET;
  int numMathType = 4, numComps[4] = { 1, 2, 3, 4 }, maxNumComp = 4;
  FC_MathType mathTypes[4] = { FC_MT_SCALAR, FC_MT_VECTOR, 
			       FC_MT_SYMTENSOR, FC_MT_TENSOR };
  int numDataType = 4;
  FC_DataType dataTypes[4] = { FC_DT_CHAR, FC_DT_INT, 
			       FC_DT_FLOAT, FC_DT_DOUBLE };
  char charData[10][4] = { {'a', 'b', 'c', 'd' } };
  int intData[10][4];
  float floatData[10][4];
  double doubleData[10][4];
  int sizes[4] = {sizeof(char), sizeof(int), sizeof(float), sizeof(double) };
  void *data[4][10];
  int temp_numDataPoint, temp_numComp;
  FC_AssociationType temp_assoc;
  FC_MathType temp_mathType;
  FC_DataType temp_dataType;
  void *temp_data;

  // get dataset & squence
  rc = fc_getDatasetByName(dataset_name, &numReturnDatasets,
			   &returnDatasets);
  fail_unless(rc == FC_SUCCESS, "aborted - cant get dataset by name");
  fail_unless(numReturnDatasets == 1,
	      "failed to find unique dataset by name");
  dataset = returnDatasets[0];
  free(returnDatasets);
  rc = fc_getSequenceByName(dataset,"time", &numReturnSequences,
			   &returnSequences);
  fail_unless(rc == FC_SUCCESS, "aborted - cant get sequence by name");
  fail_unless(numReturnSequences == 1,
	      "failed to find unique sequences by name");
  sequence = returnSequences[0];
  free(returnSequences);

  // setup some arbitrary data (hopefully different for every step, except char)
  fc_getSequenceNumStep(sequence, &numStep);
  for (k = 0; k < numStep; k++) {
    data[0][k] = &charData[k];
    data[1][k] = &intData[k];
    data[2][k] = &floatData[k];
    data[3][k] = &doubleData[k];
    for (m = 0; m < maxNumComp; m++) {
      if (k > 0)
	charData[k][m] = charData[0][m];
      intData[k][m] = k*m;
      floatData[k][m] = intData[k][m] + 0.1;
      doubleData[k][m] = intData[k][m] + 0.00000000001;
    }
  }
  
  // math types & data types    
  for (k = 0; k < numMathType; k++) {
    for (m = 0; m < numDataType; m++) {
	
      // create two sequence variables
      rc = fc_createGlobalSeqVariable(dataset, sequence, name, &numStep,
				      &seqVar);
      fail_unless(rc == FC_SUCCESS, "abort: failed to create seqVar");
      rc = fc_createGlobalSeqVariable(dataset, sequence, name2, &numStep,
				      &seqVar2);
      fail_unless(rc == FC_SUCCESS, "abort: failed to create seqVar");


      // variable 1: set data
      for (n = 0; n < numStep; n++) {
	rc = fc_setVariableData(seqVar[n], numDataPoint, numComps[k],
				assoc,  mathTypes[k], dataTypes[m],
				data[m][n]);
	fail_unless(rc == FC_SUCCESS, "failed to set data");
      }

      // variable 1: check getting by step
      for (n = 0; n < numStep; n++) {
	rc = fc_getVariableName(seqVar[n], &temp_name);
	fail_unless(rc == FC_SUCCESS, "failed to get name");
	fail_unless(!strcmp(temp_name, name), "mismatch of name");
	free(temp_name);
	rc = fc_getVariableInfo(seqVar[n], &temp_numDataPoint, &temp_numComp,
				&temp_assoc, &temp_mathType, &temp_dataType);
	fail_unless(rc == FC_SUCCESS, "failed to get seqVar info");
	fail_unless(temp_numDataPoint == numDataPoint, 
		    "mismatch of numDataPoint");
	fail_unless(temp_numComp == numComps[k], "mismatch of numComp");
	fail_unless(temp_assoc == assoc, "mismatch of assoc");
	fail_unless(temp_mathType == mathTypes[k], "mismatch of mathType");
	fail_unless(temp_dataType == dataTypes[m], "mismatch of datatype");
	rc = fc_getVariableDataPtr(seqVar[n], &temp_data);
	fail_unless(rc == FC_SUCCESS, "failed to get data on a step");
	fail_unless(!memcmp(temp_data, data[m][n], 
			    numDataPoint*numComps[k]*sizes[m]),
		    "mismatch of data");
      }


      /*FIX
      // same as above but for slice as seqvar
      for (n = 0; n < numDataPoint; n++) {
	FC_Variable *slicevar, *checkvar;
	void* slicedata;
	FC_DataType slicedt;
	FC_MathType slicemt;
	FC_AssociationType sliceat;
	FC_Sequence sliceseq;
	int slicenc, slicedp;
	int checkstep;
	
	rc = fc_getSeqVariableDataPointSliceAsSeqVariable(numStep,
							  seqVar, n, "slicevarname", &slicevar);
	fail_unless(rc == FC_SUCCESS, "failed to get slicevar");
	fail_unless(fc_isSeqVariableValid(numStep,slicevar),
		    "didn't create valid seqvar");
	fc_getVariableInfo(slicevar[0],&slicedp,&slicenc,
			   &sliceat,&slicemt,&slicedt);
	fail_unless(slicedp == 1,"mismatch of numDataPoint");
	fail_unless(slicenc == numComps[k],"mismatch of numComponent");
	fail_unless(sliceat == FC_AT_WHOLE_MESH,"mismatch of assoc");
	fail_unless(slicemt == mathTypes[k],"mismatch of mathtype");
	fail_unless(slicedt == dataTypes[m],"mismatch of datatype");
	
	fc_getSequenceFromSeqVariable(numStep,slicevar,&sliceseq);
	fail_unless(FC_HANDLE_EQUIV(sliceseq,sequence),
		    "mismatch of sequence handles");
	fc_getSeqVariableByName(dataset,"slicevarname",
				&checkstep,&checkvar);
	for (i = 0; i < numStep; i++){
	  fail_unless(FC_HANDLE_EQUIV(slicevar[i],checkvar[i]),
		      "mismatch of variable handles");
	}
	free(checkvar);
	
	for (p = 0; p < numStep; p++) {
	  rc = fc_getVariableDataPtr(slicevar[p],&slicedata);
	  fail_unless(rc == FC_SUCCESS, "failed to get slicevar data");
	  
	  for (q = 0; q < numComps[k]; q++) {
	    if (dataTypes[m] == FC_DT_CHAR) {
	      fail_unless(((char*)slicedata)[q] ==
			  charData[p][n*numComps[k]+q],
			  "mismatch of data");	      
	    }
	    else if (dataTypes[m] == FC_DT_INT) {
	      fail_unless(((int*)slicedata)[q] ==
			  intData[p][n*numComps[k]+q],
			  "mismatch of data");	      
	    }
	    else if (dataTypes[m] == FC_DT_FLOAT) {
	      fail_unless(((float*)slicedata)[q] ==
			  floatData[p][n*numComps[k]+q],
			  "mismatch of data");	      
	    }
	    else if (dataTypes[m] == FC_DT_DOUBLE) {
	      fail_unless(((double*)slicedata)[q] ==
			  doubleData[p][n*numComps[k]+q],
			  "mismatch of data");	      
	    }
	    else
	      fail_unless(0, "abort: missed a datatype");
	  }
	}
	fc_deleteSeqVariable(numStep,slicevar);
	free(slicevar);
      }
      */
	 
      // ----- check errors
      
      // set with assoc = vertes should fail
      rc = fc_setVariableData(seqVar2[0], numDataPoint, numComps[k],
			      FC_AT_VERTEX, mathTypes[k], dataTypes[m],
			      data[m]);
      fail_unless(rc < FC_SUCCESS, "should fail to set data");

      // --- all done

      
      // cleanup
      rc = fc_deleteSeqVariable(numStep, seqVar);
      fail_unless(rc == FC_SUCCESS, "failed to delete seqVar");
      free(seqVar);
      rc = fc_deleteSeqVariable(numStep, seqVar2);
      fail_unless(rc == FC_SUCCESS, "failed to delete seqVar");
      free(seqVar2);
    }
  }
}
END_TEST

// copy
START_TEST(glb_seq_copy_test)
{
  FC_ReturnCode rc;
  int k, m, n;
  char name[20] = "blue berry", copy_name[20] = "banana nut", *temp_name;
  FC_Dataset dataset1, dataset2, badDataset = { 999, 999 };
  FC_Dataset *returnDatasets;
  FC_Mesh mesh, *returnMeshes;
  int numReturnDatasets, numReturnMeshes, numElem;
  FC_Sequence sequence1, sequence2, weirdSequence, badSequence = { 999, 999 };
  FC_Sequence *returnSequences;
  int numReturnSequences;
  int numSeqVar1, numSeqVar2, temp_numSeqVar, numStep;
  FC_Variable *seqVar, *copy_seqVar, *badSeqVar;
  FC_Variable badVariable = { 999, 999 };
  FC_Variable *wholeMeshSeqVar, *elemSeqVar;
  int numDataPoint = 1; 
  FC_AssociationType assoc = FC_AT_WHOLE_DATASET;
  int numMathType = 4, numComps[4] = { 1, 2, 3, 4 }, maxNumComp = 4;
  FC_MathType mathTypes[4] = { FC_MT_SCALAR, FC_MT_VECTOR, 
			       FC_MT_SYMTENSOR, FC_MT_TENSOR };
  int numDataType = 4;
  FC_DataType dataTypes[4] = { FC_DT_CHAR, FC_DT_INT, 
			       FC_DT_FLOAT, FC_DT_DOUBLE };
  char charData[10][20*4] = { {'a', 'b', 'c', 'd', 'e', 
			       'f', 'g', 'h', 'i', 'j',
			       'k', 'l', 'm', 'n', 'o', 
			       'p', 'q', 'r', 's', 't' } };
  int intData[10][20*4];
  float floatData[10][20*4];
  double doubleData[10][20*4];
  int sizes[4] = {sizeof(char), sizeof(int), sizeof(float), sizeof(double) };
  void *data[4][10];
  int temp_numDataPoint, temp_numComp;
  FC_AssociationType temp_assoc;
  FC_MathType temp_mathType;
  FC_DataType temp_dataType;
  void *temp_data;

  // get dataset and mesh
  rc = fc_getDatasetByName(dataset_name, &numReturnDatasets,
			   &returnDatasets);
  fail_unless(rc == FC_SUCCESS, "aborted - cant get dataset by name");
  fail_unless(numReturnDatasets == 1,
	      "failed to find unique dataset by name");
  dataset1 = returnDatasets[0];
  free(returnDatasets);
  rc = fc_getSequenceByName(dataset1,"time", &numReturnSequences,
			   &returnSequences);
  fail_unless(rc == FC_SUCCESS, "aborted - cant get sequence by name");
  fail_unless(numReturnSequences == 1,
	      "failed to find unique sequences by name");
  sequence1 = returnSequences[0];
  free(returnSequences);
  rc = fc_getMeshByName(dataset1, "hex mesh", &numReturnMeshes,
			&returnMeshes);
  fail_unless(rc == FC_SUCCESS, "aborted - cant get mesh by name");
  fail_unless(numReturnMeshes == 1, "failed to find unique mesh by name");
  mesh = returnMeshes[0];
  free(returnMeshes);

  // make a second dataset & sequence to copy to
  rc = fc_copyDataset(dataset1, "copy dataset", &dataset2);
  fail_unless(rc == FC_SUCCESS, "abort: failed to copy dataset");
  rc = fc_copySequence(sequence1, dataset2, "copy sequence", &sequence2);
  fail_unless(rc == FC_SUCCESS, "abort: failed to copy sequence");

  // get a third sequence we shouldn't be able to copy to - different numStep
  // (don't delete when we're done)
  rc = fc_getSequenceByName(dataset1,"time2", &numReturnSequences,
			   &returnSequences);
  fail_unless(rc == FC_SUCCESS, "aborted - cant get sequence by name");
  fail_unless(numReturnSequences == 1,
	      "failed to find unique sequences by name");
  weirdSequence = returnSequences[0];
  free(returnSequences);

  // setup some arbitrary data (hopefully different for every step, except char)
  fc_getSequenceNumStep(sequence1, &numStep);
  for (k = 0; k < numStep; k++) {
    data[0][k] = &charData[k];
    data[1][k] = &intData[k];
    data[2][k] = &floatData[k];
    data[3][k] = &doubleData[k];
    for (m = 0; m < maxNumComp; m++) {
      if (k > 0)
	charData[k][m] = charData[0][m];
      intData[k][m] = k*m;
      floatData[k][m] = intData[k][m] + 0.1;
      doubleData[k][m] = intData[k][m] + 0.00000000001;
    }
  }

  // Make a whole mesh seq variable (double)
  rc = fc_createSeqVariable(mesh, sequence1, name, &numStep, &wholeMeshSeqVar);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create whole mesh seqvar");
  for (n = 0; n < numStep; n++) {
    rc = fc_setVariableData(wholeMeshSeqVar[n], 1, 1, FC_AT_WHOLE_MESH, 
                            FC_MT_SCALAR, FC_DT_DOUBLE, doubleData[n]);
    fail_unless(rc == FC_SUCCESS, "abort: failed to set whole mesh var data");
  }
  // Make an element seq variable (double) 
  rc = fc_createSeqVariable(mesh, sequence1, name, &numStep, &elemSeqVar);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create vert var");
  fc_getMeshNumElement(mesh, &numElem);
  fail_unless(numElem < maxNumComp, "abort: developer assumption changed");
  for (n = 0; n < numStep; n++) {
    rc = fc_setVariableData(elemSeqVar[n], numElem, 1, FC_AT_ELEMENT, 
                            FC_MT_SCALAR, FC_DT_DOUBLE, doubleData[n]);
    fail_unless(rc == FC_SUCCESS, "abort: failed to set vert var data");
  }

  numSeqVar1 = 0;
  numSeqVar2 = 0;
  // math types & data types
  for (k = 0; k < numMathType; k++) {
    for (m = 0; m < numDataType; m++) {
      
      // create a sequence variable
      rc = fc_createGlobalSeqVariable(dataset1, sequence1, name, &numStep,
				      &seqVar);
      fail_unless(rc == FC_SUCCESS, "abort: failed to create seqVar");
      for (n = 0; n < numStep; n++) {
	rc = fc_setVariableData(seqVar[n], numDataPoint, numComps[k],
				assoc,  mathTypes[k], dataTypes[m],
				data[m][n]);
	fail_unless(rc == FC_SUCCESS, "abort: failed to set data");
      }
      numSeqVar1++;
      
      // copy to same dataset & sequence & check it
      rc = fc_copyGlobalSeqVariable(numStep, seqVar, dataset1, sequence1, 
				    copy_name, &copy_seqVar);
      fail_unless(rc == FC_SUCCESS, "failed to copy to same dataset");
      numSeqVar1++;
      for (n = 0; n < numStep; n++) {
	rc = fc_getVariableName(copy_seqVar[n], &temp_name);
	fail_unless(rc == FC_SUCCESS, "failed to get name of copy");
	fail_unless(!strcmp(temp_name, copy_name), "mismatch of name");
	free(temp_name);
	rc = fc_getVariableInfo(copy_seqVar[n], &temp_numDataPoint, 
				&temp_numComp, &temp_assoc, &temp_mathType, 
				&temp_dataType);
	fail_unless(rc == FC_SUCCESS, "failed to get seqVar info");
	fail_unless(temp_numDataPoint == numDataPoint, 
		    "mismatch of numDataPoint");
	fail_unless(temp_numComp == numComps[k], "mismatch of numComp");
	fail_unless(temp_assoc == assoc, "mismatch of assoc");
	fail_unless(temp_mathType == mathTypes[k], "mismatch of mathType");
	fail_unless(temp_dataType == dataTypes[m], "mismatch of datatype");
	rc = fc_getVariableDataPtr(copy_seqVar[n], &temp_data);
	fail_unless(rc == FC_SUCCESS, "failed to get data");
	fail_unless(!memcmp(temp_data, data[m][n], 
			    numDataPoint*numComps[k]*sizes[m]),
		    "mismatch of data");
      }
      free(copy_seqVar);
      
      // copy to different dataset & different sequence & check it
      rc = fc_copyGlobalSeqVariable(numStep, seqVar, dataset2, sequence2, 
				    copy_name, &copy_seqVar);
      fail_unless(rc == FC_SUCCESS, "failed to copy to same dataset");
      numSeqVar2++;
      for (n = 0; n < numStep; n++) {
	rc = fc_getVariableName(copy_seqVar[n], &temp_name);
	fail_unless(rc == FC_SUCCESS, "failed to get name of copy");
	fail_unless(!strcmp(temp_name, copy_name), "mismatch of name");
	free(temp_name);
	rc = fc_getVariableInfo(copy_seqVar[n], &temp_numDataPoint, 
				&temp_numComp, &temp_assoc, &temp_mathType, 
				&temp_dataType);
	fail_unless(rc == FC_SUCCESS, "failed to get seqVar info");
	fail_unless(temp_numDataPoint == numDataPoint, 
		    "mismatch of numDataPoint");
	fail_unless(temp_numComp == numComps[k], "mismatch of numComp");
	fail_unless(temp_assoc == assoc, "mismatch of assoc");
	fail_unless(temp_mathType == mathTypes[k], "mismatch of mathType");
	fail_unless(temp_dataType == dataTypes[m], "mismatch of datatype");
	rc = fc_getVariableDataPtr(copy_seqVar[n], &temp_data);
	fail_unless(rc == FC_SUCCESS, "failed to get data");
	fail_unless(!memcmp(temp_data, data[m][n], 
			    numDataPoint*numComps[k]*sizes[m]),
		    "mismatch of data");
      }
      free(copy_seqVar);
      
      // fc_copyGlobalSeqVariable() to a very different sequence should fail
      rc = fc_copyGlobalSeqVariable(numStep, seqVar, dataset1, weirdSequence, 
			      copy_name, &copy_seqVar);
      fail_unless(rc != FC_SUCCESS,
		  "should fail to copy to a very different dataset");
      fail_unless(copy_seqVar == NULL, "should return null seqvar");
      
      // --- special cases
      
      // copy with NULL name will use original's name
      rc = fc_copyGlobalSeqVariable(numStep, seqVar, dataset1, sequence1, 
			      NULL, &copy_seqVar);
      fail_unless(rc == FC_SUCCESS, "failed to copy to same dataset");
      for (n = 0; n < numStep; n++) {
	rc = fc_getVariableName(copy_seqVar[n], &temp_name);
	fail_unless(rc == FC_SUCCESS, "failed to get name of copy");
	fail_unless(!strcmp(temp_name, name), "mismatch of name");
	free(temp_name);
      }
      rc = fc_deleteSeqVariable(numStep, copy_seqVar);
      fail_unless(rc == FC_SUCCESS, "failed to delete copied seqVar");
      free(copy_seqVar);
      
      // ---- check errors
      
      // construct a bad seqVar (make one of the vars bad)
      badSeqVar = (FC_Variable*)malloc(numStep*sizeof(FC_Variable));
      for (n = 0; n < numStep; n++)
	badSeqVar[n] = seqVar[n];
      badSeqVar[numStep/2] = badVariable;
      
      // fc_copyVariable() -- bad args
      rc = fc_copyGlobalSeqVariable(numStep-1, seqVar, dataset1, sequence1, 
			      copy_name, &copy_seqVar);
      fail_unless(rc != FC_SUCCESS, "should fail with wrong numStep");
      fail_unless(copy_seqVar == NULL, "fail should return null seqVar");
      rc = fc_copyGlobalSeqVariable(numStep, NULL, dataset1, sequence1, 
			      copy_name, &copy_seqVar);
      fail_unless(rc != FC_SUCCESS, "should fail with null seqVar");
      fail_unless(copy_seqVar == NULL, "fail should return null seqVar");
      rc = fc_copyGlobalSeqVariable(numStep, badSeqVar, dataset1, sequence1, 
			      copy_name, &copy_seqVar);
      fail_unless(rc != FC_SUCCESS, "should fail with badSeqVar");
      fail_unless(copy_seqVar == NULL, "fail should return null seqVar");
      rc = fc_copyGlobalSeqVariable(numStep, seqVar, badDataset, sequence1, 
			      copy_name, &copy_seqVar);
      fail_unless(rc != FC_SUCCESS, "should fail with bad dataset");
      fail_unless(copy_seqVar == NULL, "fail should return null seqVar");
      rc = fc_copyGlobalSeqVariable(numStep, seqVar, dataset1, badSequence, 
			      copy_name, &copy_seqVar);
      fail_unless(rc != FC_SUCCESS, "should fail with bad sequence");
      fail_unless(copy_seqVar == NULL, "fail should return null seqVar");
      rc = fc_copyGlobalSeqVariable(numStep, seqVar, dataset1, sequence1, 
			      copy_name, NULL);
      fail_unless(rc != FC_SUCCESS, "should fail with null seqvar");
      
      // --- done
      
      // cleanup
      free(seqVar);
      free(badSeqVar);
    }
  }

  // It's o.k. to copy a whole mesh var onto a dataset
  rc = fc_copyGlobalSeqVariable(numStep, wholeMeshSeqVar, dataset1, sequence1,
                                copy_name, &copy_seqVar);
  fail_unless(rc == FC_SUCCESS, "failed to copy whole mesh seq var to a dataset");
  numSeqVar1++;
  for (n = 0; n < numStep; n++) {
    rc = fc_getVariableName(copy_seqVar[n], &temp_name);
    fail_unless(rc == FC_SUCCESS, "failed to get name of copy");
    fail_unless(!strcmp(temp_name, copy_name), "mismatch of name");
    free(temp_name);
    rc = fc_getVariableInfo(copy_seqVar[n], &temp_numDataPoint, 
                            &temp_numComp, &temp_assoc, &temp_mathType, 
                            &temp_dataType);
    fail_unless(rc == FC_SUCCESS, "failed to get variable info");
    fail_unless(temp_numDataPoint == 1, "mismatch of numDataPoint");
    fail_unless(temp_numComp == 1, "mismatch of numComp");
    fail_unless(temp_assoc == FC_AT_WHOLE_DATASET, "mismatch of assoc");
    fail_unless(temp_mathType == FC_MT_SCALAR, "mismatch of mathType");
    fail_unless(temp_dataType == FC_DT_DOUBLE, "mismatch of datatype");
    rc = fc_getVariableDataPtr(copy_seqVar[n], &temp_data);
    fail_unless(rc == FC_SUCCESS, "failed to get data");
    fail_unless(!memcmp(temp_data, doubleData[n], sizeof(double)), 
                "mismatch of data");
  }
  free(copy_seqVar);

  // It's NOT o.k. to copy an elem seq var onto a dataset
  rc = fc_copyGlobalSeqVariable(numStep, elemSeqVar, dataset1, sequence1,
                                copy_name, &copy_seqVar);
  fail_unless(rc != FC_SUCCESS, "should fail to copy elem seq var to dataset");
  fail_unless(copy_seqVar == NULL, "fail should return NULL");

  // a little more checking & cleanup
  rc = fc_getNumGlobalSeqVariable(dataset2, &temp_numSeqVar);
  fail_unless(rc == FC_SUCCESS, "shouldn't fail to get numSeqVars");
  fail_unless(temp_numSeqVar == numSeqVar2, "mismatch of numSeqVar");
  rc = fc_deleteDataset(dataset2);
  fail_unless(rc == FC_SUCCESS, "failed to delete second dataset");
  rc = fc_getNumGlobalSeqVariable(dataset1, &temp_numSeqVar);
  fail_unless(rc == FC_SUCCESS, "shouldn't fail to get numSeqVar");
  fail_unless(temp_numSeqVar == numSeqVar1, "mismatch of numSeqVar");

  // cleanup
  free(wholeMeshSeqVar);
  free(elemSeqVar);
}
END_TEST

// --- creating new vars from other vars

// convert from basic var to seq var test
// test representative case
START_TEST(vars_to_seq_var)
{
  FC_ReturnCode rc;
  int i, j, k;
  char varName[20] = "peanut", seqVarName[20] = "butter", *temp_name;
  char copyName[20] = "fudge"; // o.k., I may be hungry
  FC_Dataset dataset, *returnDatasets;
  int numReturnDatasets;
  FC_Sequence *sequences, temp_sequence, badSequence = { 999, 999 };
  FC_Mesh *meshes, mesh;
  FC_Variable vars[2][10]; // index by seq then step
  FC_Variable copyVars[10], badVars[10], badVar = { 999, 999 };
  int numMesh, numSequence, numSteps[3], numElem = 2;
  int temp_numVar, *temp_numSteps, temp_numSeqVar;
  FC_Variable *new_seqVars[3], **temp_seqVars, *temp_seqVar;
  double datas[10][2] = { { 0.0, 12.1 } }, *temp_data;
  _FC_MeshSlot *meshSlot;
  _FC_VarSlot* varSlot;

  // get dataset, sequences and meshes
  rc = fc_getDatasetByName(dataset_name, &numReturnDatasets,
			   &returnDatasets);
  fail_unless(rc == FC_SUCCESS, "aborted - cant get dataset by name");
  fail_unless(numReturnDatasets == 1,
	      "failed to find unique dataset by name");
  dataset = returnDatasets[0];
  free(returnDatasets);
  rc = fc_getSequences(dataset, &numSequence, &sequences);
  fail_unless(rc == FC_SUCCESS, "abort: failed to get sequences");
  fail_unless(numSequence == 2, "abort: expected 2 sequences");
  rc = fc_getMeshes(dataset, &numMesh, &meshes);
  fail_unless(rc == FC_SUCCESS, "abort: failed to get meshes");
  fail_unless(numMesh == 8, "abort: expected 8 meshes"); 
  for (i = 0; i < numSequence; i++)
    fc_getSequenceNumStep(sequences[i], &numSteps[i]);
  fail_unless(numSteps[0] != numSteps[1], 
	      "abort: expect sequences to have different number of steps");
  mesh = meshes[7];
  meshSlot = _fc_getMeshSlot(mesh);

  // create variables
  for (i = 0; i < numSequence; i++) {
    for (j = 0; j < numSteps[i]; j++) {
      if (j > 0) 
	for (k = 0; k < numElem; k++)
	  datas[j][k] = datas[0][k] + j;
      fc_createVariable(mesh, varName, &vars[i][j]);
      rc = fc_setVariableData(vars[i][j], numElem, 1, FC_AT_ELEMENT, 
			      FC_MT_SCALAR, FC_DT_DOUBLE, datas[j]);
      fail_unless(rc == FC_SUCCESS, "failed to create basic variable");
    }
  }

  // ---- First converted seq var - give new name
  rc = fc_convertVariablesToSeqVariable(numSteps[0], vars[0], sequences[0], 
					 seqVarName, &new_seqVars[0]);
  fail_unless(rc == FC_SUCCESS, "failed to convert to seq var");
  fail_unless(fc_isSeqVariableValid(numSteps[0], new_seqVars[0]),
	      "seq var failed valid test");
  for (i = 0; i < numSteps[0]; i++) {
    fail_unless(FC_HANDLE_EQUIV(new_seqVars[0][i], vars[0][i]),
		"mismatch of handle");
    fc_getVariableName(new_seqVars[0][i], &temp_name);
    fail_unless(!strcmp(temp_name, seqVarName), "mismatch of name");
    free(temp_name);
    fc_getVariableDataPtr(new_seqVars[0][i], (void**)&temp_data);
    fail_unless(!memcmp(temp_data, datas[i], numElem*sizeof(double)),
		"mismatch of data");
  }
  // test state of mesh
  rc = fc_getNumVariable(mesh, &temp_numVar);
  fail_unless(rc == FC_SUCCESS, "failed to get number of variables");
  fail_unless(temp_numVar == numSteps[1], "should have vars for 2nd seq");
  for (i = 0; i < numSteps[1]; i++)
    fail_unless(meshSlot->basicVarIDs[i] == vars[1][i].slotID,
		"bad varIDs on meshSlot");
  rc = fc_getSeqVariables(mesh, &temp_numSeqVar, &temp_numSteps, 
			  &temp_seqVars);
  fail_unless(rc == FC_SUCCESS, "failed to get seq vars");
  fail_unless(temp_numSeqVar == 1, "should be one seq var");
  fail_unless(temp_numSteps[0] == numSteps[0], "mismatch of step");
  for (i = 0; i < numSteps[0]; i++) 
    fail_unless(FC_HANDLE_EQUIV(temp_seqVars[0][i], new_seqVars[0][i]),
		"mismatch of handle");
  free(temp_numSteps);
  for (i = 0; i < temp_numSeqVar; i++)
    free(temp_seqVars[i]);
  free(temp_seqVars);

  // ---- Second converted seq var - use old name
  rc = fc_convertVariablesToSeqVariable(numSteps[1], vars[1], sequences[1], 
					 NULL, &new_seqVars[1]);
  fail_unless(rc == FC_SUCCESS, "failed to convert to seq var");
  fail_unless(fc_isSeqVariableValid(numSteps[1], new_seqVars[1]),
	      "seq var failed valid test");
  for (i = 0; i < numSteps[1]; i++) {
    fail_unless(FC_HANDLE_EQUIV(new_seqVars[1][i], vars[1][i]),
		"mismatch of handle");
    fc_getVariableName(new_seqVars[1][i], &temp_name);
    fail_unless(!strcmp(temp_name, varName), "mismatch of name");
    free(temp_name);
    fc_getVariableDataPtr(new_seqVars[1][i], (void**)&temp_data);
    fail_unless(!memcmp(temp_data, datas[i], numElem*sizeof(double)),
		"mismatch of data");
  }
  // test state of mesh
  rc = fc_getNumVariable(mesh, &temp_numVar);
  fail_unless(rc == FC_SUCCESS, "failed to get number of variables");
  fail_unless(temp_numVar == 0, "should not have any vars left");
  fail_unless(meshSlot->basicVarIDs == NULL, "bad varIDs on meshSlot");
  rc = fc_getSeqVariables(mesh, &temp_numSeqVar, &temp_numSteps, 
			  &temp_seqVars);
  fail_unless(rc == FC_SUCCESS, "failed to get seq vars");
  fail_unless(temp_numSeqVar == 2, "should be two seq vars");
  for (i = 0; i < temp_numSeqVar; i++) {
    fail_unless(temp_numSteps[i] == numSteps[i], "mismatch of step");
    fc_getSequenceFromSeqVariable(temp_numSteps[i], temp_seqVars[i],
				  &temp_sequence);
    fail_unless(FC_HANDLE_EQUIV(temp_sequence, sequences[i]),
		"mismatch of sequence");
    for (j = 0; j < temp_numSteps[i]; j++)
      fail_unless(FC_HANDLE_EQUIV(temp_seqVars[i][j], new_seqVars[i][j]),
		  "mismatch of handle");
  }
  free(temp_numSteps);
  for (i = 0; i < temp_numSeqVar; i++)
    free(temp_seqVars[i]);
  free(temp_seqVars);

  // ---- Third converted seq var - put one more on first sequence
  // also test copying of steps of sequence var -> should get basic var
  for (i = 0; i < numSteps[1]; i++) {
    rc = fc_copyVariable(new_seqVars[1][i], mesh, copyName, &copyVars[i]);
    fail_unless(rc == FC_SUCCESS, "failed to copy var");
    varSlot = _fc_getVarSlot(copyVars[i]);
    fail_unless(FC_HANDLE_EQUIV(varSlot->sequence, FC_NULL_SEQUENCE) && 
		varSlot->stepID == -1, "copy of step should not be seqVar");
  }
  numSteps[2] = numSteps[1];
  rc = fc_convertVariablesToSeqVariable(numSteps[1], copyVars, sequences[1], 
					 NULL, &new_seqVars[2]);
  fail_unless(rc == FC_SUCCESS, "failed to convert to seq var");
  fail_unless(fc_isSeqVariableValid(numSteps[2], new_seqVars[2]),
	      "seq var failed valid test");
  for (i = 0; i < numSteps[2]; i++) {
    fail_unless(FC_HANDLE_EQUIV(new_seqVars[2][i], copyVars[i]),
		"mismatch of handle");
    fc_getVariableName(new_seqVars[2][i], &temp_name);
    fail_unless(!strcmp(temp_name, copyName), "mismatch of name");
    free(temp_name);
    fc_getVariableDataPtr(new_seqVars[2][i], (void**)&temp_data);
    fail_unless(!memcmp(temp_data, datas[i], numElem*sizeof(double)),
		"mismatch of data");
  }
  // test state of mesh
  rc = fc_getNumVariable(mesh, &temp_numVar);
  fail_unless(rc == FC_SUCCESS, "failed to get number of variables");
  fail_unless(temp_numVar == 0, "should not have any vars left");
  fail_unless(meshSlot->basicVarIDs == NULL, "bad varIDs on meshSlot");
  rc = fc_getSeqVariables(mesh, &temp_numSeqVar, &temp_numSteps, 
			  &temp_seqVars);
  fail_unless(rc == FC_SUCCESS, "failed to get seq vars");
  fail_unless(temp_numSeqVar == 3, "should be three seq var");
  for (i = 0; i < 3; i++) {
    if (i < 2)
      k = i;
    else 
      k = 1;
    fail_unless(temp_numSteps[i] == numSteps[k], "mismatch of step");
    fc_getSequenceFromSeqVariable(temp_numSteps[i], temp_seqVars[i],
				  &temp_sequence);
    fail_unless(FC_HANDLE_EQUIV(temp_sequence, sequences[k]),
		"mismatch of sequence");
    for (j = 0; j < temp_numSteps[i]; j++)
	fail_unless(FC_HANDLE_EQUIV(temp_seqVars[i][j], new_seqVars[i][j]),
		    "mismatch of handle");
  }
  free(temp_numSteps);
  for (i = 0; i < temp_numSeqVar; i++)
    free(temp_seqVars[i]);
  free(temp_seqVars);

  // --- test errors

  // make a set of basic vars for testing purposes
  for (i = 0; i < numSteps[1]; i++) {
    rc = fc_copyVariable(new_seqVars[1][i], mesh, copyName, &copyVars[i]);
    fail_unless(rc == FC_SUCCESS, "abort: failed to copy var");
    varSlot = _fc_getVarSlot(copyVars[i]);
    fail_unless(FC_HANDLE_EQUIV(varSlot->sequence, FC_NULL_SEQUENCE) && 
		varSlot->stepID == -1, "copy of step should not be seqVar");
    badVars[i] = copyVars[i];
  }
  // bad vars has one var on another mesh
  rc = fc_copyVariable(new_seqVars[1][5], meshes[1], copyName, &badVars[5]);
  fail_unless(rc == FC_SUCCESS, "abort: failed to copy to diff mesh");

  // bad input -- bad vars
  rc = fc_convertVariablesToSeqVariable(numSteps[0], copyVars, sequences[1], 
					seqVarName, &temp_seqVar);
  fail_unless(rc != FC_SUCCESS, "should fail if numStep & seq do not agree");
  fail_unless(temp_seqVar == NULL, "fail should return nulls");
  rc = fc_convertVariablesToSeqVariable(numSteps[1], badVars, sequences[1], 
					 seqVarName, &temp_seqVar);
  fail_unless(rc != FC_SUCCESS, "should fail if vars not all on same mesh");
  fail_unless(temp_seqVar == NULL, "fail should return nulls");
  badVars[5] = badVar;
  rc = fc_convertVariablesToSeqVariable(numSteps[1], badVars, sequences[1], 
					 seqVarName, &temp_seqVar);
  fail_unless(rc != FC_SUCCESS, "should fail if one bad var in list");
  fail_unless(temp_seqVar == NULL, "fail should return nulls");

  // bad input -- bad sequence
  rc = fc_convertVariablesToSeqVariable(numSteps[1], copyVars, badSequence, 
					seqVarName, &temp_seqVar);
  fail_unless(rc != FC_SUCCESS, "should fail if bad sequence");
  fail_unless(temp_seqVar == NULL, "fail should return nulls");

  // bad input -- missing return variables
  rc = fc_convertVariablesToSeqVariable(numSteps[1], copyVars, sequences[1], 
					 seqVarName, NULL);
  fail_unless(rc != FC_SUCCESS, "should fail if no seqVar");

  free(meshes);
  free(sequences);
  for(i = 0; i < 3; i++)
    free(new_seqVars[i]);
 }
END_TEST

// convert from global var to global seq var test
// test representative case
// this is essentially a copy of vars_to_seq_var
START_TEST(glb_vars_to_seq_var)
{
  FC_ReturnCode rc;
  int i, j, k;
  char varName[20] = "peanut", seqVarName[20] = "butter", *temp_name;
  char copyName[20] = "fudge"; // o.k., I may be hungry
  FC_Dataset dataset, *returnDatasets, dataset2;
  int numReturnDatasets;
  FC_Sequence *sequences, temp_sequence, badSequence = { 999, 999 };
  FC_Variable vars[2][10]; // index by seq then step
  FC_Variable copyVars[10], badVars[10], badVar = { 999, 999 };
  int numSequence, numSteps[3];
  int temp_numVar, *temp_numSteps, temp_numSeqVar;
  FC_Variable *new_seqVars[3], **temp_seqVars, *temp_seqVar;
  double datas[10], *temp_data;
  _FC_DsSlot *dsSlot;
  _FC_VarSlot* varSlot;

  // get datasets & sequences
  rc = fc_getDatasetByName(dataset_name, &numReturnDatasets,
			   &returnDatasets);
  fail_unless(numReturnDatasets == 1,
	      "failed to find unique dataset by name");
  dataset = returnDatasets[0];
  free(returnDatasets);
  fail_unless(rc == FC_SUCCESS, "abort: failed to get dataset by name");
  rc = fc_getSequences(dataset, &numSequence, &sequences);
  fail_unless(rc == FC_SUCCESS, "abort: failed to get sequences");
  fail_unless(numSequence == 2, "abort: expected 2 sequences");
  for (i = 0; i < numSequence; i++)
    fc_getSequenceNumStep(sequences[i], &numSteps[i]);
  fail_unless(numSteps[0] != numSteps[1], 
	      "abort: expect sequences to have different number of steps");
  dsSlot = _fc_getDsSlot(dataset);
  rc = fc_createDataset("new dataset", &dataset2);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create dataset");

  // create variables
  for (i = 0; i < numSequence; i++) {
    for (j = 0; j < numSteps[i]; j++) {
      datas[j] = j;
      fc_createGlobalVariable(dataset, varName, &vars[i][j]);
      rc = fc_setVariableData(vars[i][j], 1, 1, FC_AT_WHOLE_DATASET, 
			      FC_MT_SCALAR, FC_DT_DOUBLE, &datas[j]);
      fail_unless(rc == FC_SUCCESS, "failed to create global variable");
    }
  }
 
  // ---- First converted seq var - give new name
  rc = fc_convertVariablesToSeqVariable(numSteps[0], vars[0], sequences[0], 
					 seqVarName, &new_seqVars[0]);
  fail_unless(rc == FC_SUCCESS, "failed to convert to seq var");
  fail_unless(fc_isSeqVariableValid(numSteps[0], new_seqVars[0]),
	      "seq var failed valid test");
  for (i = 0; i < numSteps[0]; i++) {
    fail_unless(FC_HANDLE_EQUIV(new_seqVars[0][i], vars[0][i]),
		"mismatch of handle");
    fc_getVariableName(new_seqVars[0][i], &temp_name);
    fail_unless(!strcmp(temp_name, seqVarName), "mismatch of name");
    free(temp_name);
    fc_getVariableDataPtr(new_seqVars[0][i], (void**)&temp_data);
    fail_unless(!memcmp(temp_data, &datas[i], sizeof(double)),
		"mismatch of data");
  }
  // test state of dataset
  rc = fc_getNumGlobalVariable(dataset, &temp_numVar);
  fail_unless(rc == FC_SUCCESS, "failed to get number of variables");
  fail_unless(temp_numVar == numSteps[1], "should have vars for 2nd seq");
  for (i = 0; i < numSteps[1]; i++)
    fail_unless(dsSlot->basicVarIDs[i] == vars[1][i].slotID,
		"bad varIDs on dsSlot");
  rc = fc_getGlobalSeqVariables(dataset, &temp_numSeqVar, &temp_numSteps, 
				&temp_seqVars);
  fail_unless(rc == FC_SUCCESS, "failed to get seq vars");
  fail_unless(temp_numSeqVar == 1, "should be one seq var");
  fail_unless(temp_numSteps[0] == numSteps[0], "mismatch of step");
  for (i = 0; i < numSteps[0]; i++) 
    fail_unless(FC_HANDLE_EQUIV(temp_seqVars[0][i], new_seqVars[0][i]),
		"mismatch of handle");
  free(temp_numSteps);
  for (i = 0; i < temp_numSeqVar; i++)
    free(temp_seqVars[i]);
  free(temp_seqVars);

  // ---- Second converted seq var - use old name
  rc = fc_convertVariablesToSeqVariable(numSteps[1], vars[1], sequences[1], 
					 NULL, &new_seqVars[1]);
  fail_unless(rc == FC_SUCCESS, "failed to convert to seq var");
  fail_unless(fc_isSeqVariableValid(numSteps[1], new_seqVars[1]),
	      "seq var failed valid test");
  for (i = 0; i < numSteps[1]; i++) {
    fail_unless(FC_HANDLE_EQUIV(new_seqVars[1][i], vars[1][i]),
		"mismatch of handle");
    fc_getVariableName(new_seqVars[1][i], &temp_name);
    fail_unless(!strcmp(temp_name, varName), "mismatch of name");
    free(temp_name);
    fc_getVariableDataPtr(new_seqVars[1][i], (void**)&temp_data);
    fail_unless(!memcmp(temp_data, &datas[i], sizeof(double)),
		"mismatch of data");
  }
  // test state of dataset
  rc = fc_getNumGlobalVariable(dataset, &temp_numVar);
  fail_unless(rc == FC_SUCCESS, "failed to get number of variables");
  fail_unless(temp_numVar == 0, "should not have any vars left");
  fail_unless(dsSlot->basicVarIDs == NULL, "bad varIDs on dsSlot");
  rc = fc_getGlobalSeqVariables(dataset, &temp_numSeqVar, &temp_numSteps, 
			  &temp_seqVars);
  fail_unless(rc == FC_SUCCESS, "failed to get seq vars");
  fail_unless(temp_numSeqVar == 2, "should be two seq vars");
  for (i = 0; i < temp_numSeqVar; i++) {
    fail_unless(temp_numSteps[i] == numSteps[i], "mismatch of step");
    fc_getSequenceFromSeqVariable(temp_numSteps[i], temp_seqVars[i],
				  &temp_sequence);
    fail_unless(FC_HANDLE_EQUIV(temp_sequence, sequences[i]),
		"mismatch of sequence");
    for (j = 0; j < temp_numSteps[i]; j++)
      fail_unless(FC_HANDLE_EQUIV(temp_seqVars[i][j], new_seqVars[i][j]),
		  "mismatch of handle");
  }
  free(temp_numSteps);
  for (i = 0; i < temp_numSeqVar; i++)
    free(temp_seqVars[i]);
  free(temp_seqVars);

  // ---- Third converted seq var - put one more on first sequence
  // also test copying of steps of sequence var -> should get basic var
  for (i = 0; i < numSteps[1]; i++) {
    rc = fc_copyGlobalVariable(new_seqVars[1][i], dataset, copyName, &copyVars[i]);
    fail_unless(rc == FC_SUCCESS, "failed to copy var");
    varSlot = _fc_getVarSlot(copyVars[i]);
    fail_unless(FC_HANDLE_EQUIV(varSlot->sequence, FC_NULL_SEQUENCE) && 
		varSlot->stepID == -1, "copy of step should not be seqVar");
  }
  numSteps[2] = numSteps[1];
  rc = fc_convertVariablesToSeqVariable(numSteps[1], copyVars, sequences[1], 
					 NULL, &new_seqVars[2]);
  fail_unless(rc == FC_SUCCESS, "failed to convert to seq var");
  fail_unless(fc_isSeqVariableValid(numSteps[2], new_seqVars[2]),
	      "seq var failed valid test");
  for (i = 0; i < numSteps[2]; i++) {
    fail_unless(FC_HANDLE_EQUIV(new_seqVars[2][i], copyVars[i]),
		"mismatch of handle");
    fc_getVariableName(new_seqVars[2][i], &temp_name);
    fail_unless(!strcmp(temp_name, copyName), "mismatch of name");
    free(temp_name);
    fc_getVariableDataPtr(new_seqVars[2][i], (void**)&temp_data);
    fail_unless(!memcmp(temp_data, &datas[i], sizeof(double)),
		"mismatch of data");
  }
  // test state of dataset
  rc = fc_getNumGlobalVariable(dataset, &temp_numVar);
  fail_unless(rc == FC_SUCCESS, "failed to get number of variables");
  fail_unless(temp_numVar == 0, "should not have any vars left");
  fail_unless(dsSlot->basicVarIDs == NULL, "bad varIDs on dsSlot");
  rc = fc_getGlobalSeqVariables(dataset, &temp_numSeqVar, &temp_numSteps, 
				&temp_seqVars);
  fail_unless(rc == FC_SUCCESS, "failed to get seq vars");
  fail_unless(temp_numSeqVar == 3, "should be three seq var");
  for (i = 0; i < 3; i++) {
    if (i < 2)
      k = i;
    else 
      k = 1;
    fail_unless(temp_numSteps[i] == numSteps[k], "mismatch of step");
    fc_getSequenceFromSeqVariable(temp_numSteps[i], temp_seqVars[i],
				  &temp_sequence);
    fail_unless(FC_HANDLE_EQUIV(temp_sequence, sequences[k]),
		"mismatch of sequence");
    for (j = 0; j < temp_numSteps[i]; j++)
	fail_unless(FC_HANDLE_EQUIV(temp_seqVars[i][j], new_seqVars[i][j]),
		    "mismatch of handle");
  }
  free(temp_numSteps);
  for (i = 0; i < temp_numSeqVar; i++)
    free(temp_seqVars[i]);
  free(temp_seqVars);
 
  // --- test errors

  // make a set of basic vars for testing purposes
  for (i = 0; i < numSteps[1]; i++) {
    rc = fc_copyGlobalVariable(new_seqVars[1][i], dataset, copyName, &copyVars[i]);
    fail_unless(rc == FC_SUCCESS, "abort: failed to copy var");
    varSlot = _fc_getVarSlot(copyVars[i]);
    fail_unless(FC_HANDLE_EQUIV(varSlot->sequence, FC_NULL_SEQUENCE) && 
		varSlot->stepID == -1, "copy of step should not be seqVar");
    badVars[i] = copyVars[i];
  }
  // bad vars has one var on another dataset
  rc = fc_copyGlobalVariable(new_seqVars[1][5], dataset2, copyName, &badVars[5]);
  fail_unless(rc == FC_SUCCESS, "abort: failed to copy to diff dataset");

  // bad input -- bad vars
  rc = fc_convertVariablesToSeqVariable(numSteps[0], copyVars, sequences[1], 
					seqVarName, &temp_seqVar);
  fail_unless(rc != FC_SUCCESS, "should fail if numStep & seq do not agree");
  fail_unless(temp_seqVar == NULL, "fail should return nulls");
  rc = fc_convertVariablesToSeqVariable(numSteps[1], badVars, sequences[1], 
					 seqVarName, &temp_seqVar);
  fail_unless(rc != FC_SUCCESS, "should fail if vars not all on same dataset");
  fail_unless(temp_seqVar == NULL, "fail should return nulls");
  badVars[5] = badVar;
  rc = fc_convertVariablesToSeqVariable(numSteps[1], badVars, sequences[1], 
					 seqVarName, &temp_seqVar);
  fail_unless(rc != FC_SUCCESS, "should fail if one bad var in list");
  fail_unless(temp_seqVar == NULL, "fail should return nulls");

  // bad input -- bad sequence
  rc = fc_convertVariablesToSeqVariable(numSteps[1], copyVars, badSequence, 
					seqVarName, &temp_seqVar);
  fail_unless(rc != FC_SUCCESS, "should fail if bad sequence");
  fail_unless(temp_seqVar == NULL, "fail should return nulls");

  // bad input -- missing return variables
  rc = fc_convertVariablesToSeqVariable(numSteps[1], copyVars, sequences[1], 
					 seqVarName, NULL);
  fail_unless(rc != FC_SUCCESS, "should fail if no seqVar");

  free(sequences);
  for(i = 0; i < 3; i++)
    free(new_seqVars[i]);
 }
END_TEST

// split variable into components
START_TEST(split_merge_comps)
{
  FC_ReturnCode rc;
  int i, j;
  FC_Dataset dataset, *returnDatasets;
  int numReturnDatasets;
  FC_Mesh mesh1, mesh2;
  FC_Variable vector_var;
  FC_Variable scalar_var1, scalar_var2, scalar_var3, scalar_var4;
  FC_Mesh *returnMeshes;
  int numReturnMeshes;
  FC_Variable *comps, temp_comps[2], new_var;
  int numElem, numVert;
  int numDataPoint, numComp, temp_numDataPoint, temp_numComp;
  FC_AssociationType assoc, temp_assoc;
  FC_MathType mathtype, temp_mathtype, new_mathtype = FC_MT_VECTOR;
  FC_DataType datatype, temp_datatype;
  double* data, *temp_data;
  char new_name[] = "hopeful", *old_name, *temp_name;
  char compNames[3][20] = { "_c0", "_c1", "_c2" }; 
  int name_len;
  double dbl_data[50] = { 0 }; // this data isn't really used
  float  flt_data[50] = { 0 }; // this data isn't really used

  // get dataset and meshes
  rc = fc_getDatasetByName(dataset_name, &numReturnDatasets,
			   &returnDatasets);
  fail_unless(rc == FC_SUCCESS, "aborted - cant get dataset by name");
  fail_unless(numReturnDatasets == 1,
	      "failed to find unique dataset by name");
  dataset = returnDatasets[0];
  free(returnDatasets);
  rc = fc_getMeshByName(dataset, "hex mesh", &numReturnMeshes,
			&returnMeshes);
  fail_unless(rc == FC_SUCCESS, "aborted - cant get mesh by name");
  fail_unless(numReturnMeshes == 1, "failed to find unique mesh by name");
  mesh1 = returnMeshes[0];
  free(returnMeshes);
  rc = fc_getMeshInfo(mesh1, NULL, NULL, &numVert, &numElem, NULL);
  fail_unless(rc == FC_SUCCESS, "abort: failed to get mesh info");
  rc = fc_getMeshByName(dataset, "tri mesh", &numReturnMeshes,
			&returnMeshes);
  fail_unless(rc == FC_SUCCESS, "aborted - cant get mesh by name");
  fail_unless(numReturnMeshes == 1, "failed to find unique mesh by name");
  mesh2 = returnMeshes[0];
  free(returnMeshes);
  // create some vars
  // vector var is on verts, of type double
  rc = fc_getMeshCoordsAsVariable(mesh1, &vector_var);
  fail_unless(rc == FC_SUCCESS, "abort: failed to get coords as a var");
  fc_getVariableName(vector_var, &old_name);
  name_len = strlen(old_name);
  fc_getVariableInfo(vector_var, &numDataPoint, &numComp, &assoc, &mathtype,
		     &datatype);
  rc = fc_getVariableDataPtr(vector_var, (void**)&data);
  fail_unless(rc == FC_SUCCESS, "abort: failed to query vector variable");
  // scalar var1 is on verts, of type double
  rc = fc_createVariable(mesh1, "scalar vert var", &scalar_var1);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create var");
  rc = fc_setVariableData(scalar_var1, numVert, 1, FC_AT_VERTEX,
			  FC_MT_SCALAR, FC_DT_DOUBLE, dbl_data);
  fail_unless(rc == FC_SUCCESS, "abort: failed to set var data");
  // scalar var2 is on verts, of type float
  rc = fc_createVariable(mesh1, "scalar vert var float", &scalar_var2);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create var");
  rc = fc_setVariableData(scalar_var2, numVert, 1, FC_AT_VERTEX,
			  FC_MT_SCALAR, FC_DT_FLOAT, flt_data);
  fail_unless(rc == FC_SUCCESS, "abort: failed to set var data");
  // scalar var3 is on elems, of type double
  rc = fc_createVariable(mesh1, "scalar elem var", &scalar_var3);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create var");
  rc = fc_setVariableData(scalar_var3, numElem, 1, FC_AT_ELEMENT,
			  FC_MT_SCALAR, FC_DT_DOUBLE, dbl_data);
  fail_unless(rc == FC_SUCCESS, "abort: failed to set var data");
  // scalar var4 is on elems of mesh2, of type double
  // (is on elems because mesh1 and mesh 2 have same # of elems)
  rc = fc_createVariable(mesh2, "scalar elem var on other mesh", &scalar_var4);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create var");
  rc = fc_setVariableData(scalar_var4, numElem, 1, FC_AT_ELEMENT,
			  FC_MT_SCALAR, FC_DT_DOUBLE, dbl_data);
  fail_unless(rc == FC_SUCCESS, "abort: failed to set var data");
  
  // --- test fc_createComponentVariables()
  // should fail to work on scalar variable
  rc = fc_createComponentVariables(scalar_var1, &temp_numComp, &comps);
  fail_unless(rc != FC_SUCCESS, "should not be able to get comps of scalar var");
  fail_unless(temp_numComp == -1, "failure should return -1 numComp");
  fail_unless(comps == NULL, "failure should return comps = NULL");
  // should will work on vector variable
  rc = fc_createComponentVariables(vector_var, &temp_numComp, &comps);
  fail_unless(rc == FC_SUCCESS, "fc_fetComponentVariables failed");
  fail_unless(temp_numComp == numComp, "mismatch in number of components");
  // check components
  for (i = 0; i < numComp; i++) {
    rc = fc_getVariableName(comps[i], &temp_name);
    fail_unless(rc == FC_SUCCESS, "failed to get name of component var");
    fail_unless(!strncmp(temp_name, old_name, name_len), 
		"mismatch of root name");
    fail_unless(!strcmp(&temp_name[name_len], compNames[i]), 
		"mismatch of name postfix");
    free(temp_name);
    rc = fc_getVariableInfo(comps[i], &temp_numDataPoint, &temp_numComp,
			      &temp_assoc, &temp_mathtype, &temp_datatype);
    fail_unless(rc == FC_SUCCESS, "failed to describe a component");
    fail_unless(temp_numDataPoint == numDataPoint, "mismatch of numDataPoint");
    fail_unless(temp_numComp == 1, "numComp should be 1 for component");
    fail_unless(temp_assoc == assoc, "mismatch of assoc");
    fail_unless(temp_mathtype == FC_MT_SCALAR, "mathtype should be SCALAR");
    fail_unless(temp_datatype == datatype, "mismatch of datatype");
    rc = fc_getVariableDataPtr(comps[i], (void**)&temp_data);
    fail_unless(rc == FC_SUCCESS, "failed to get data of a component");
    for (j = 0; j < numDataPoint; j++)
      fail_unless(data[j*numComp+i] == temp_data[j], "mismatch of data");
  }

  // --- test fc_mergeComponentVariables()
  // should fail if they aren't all scalars
  temp_comps[0] = scalar_var1;
  temp_comps[1] = vector_var;
  rc = fc_mergeComponentVariables(2, temp_comps, new_name, new_mathtype,
				  &new_var);
  fail_unless(rc != FC_SUCCESS, "should fail if not scalar");
  fail_unless(FC_HANDLE_EQUIV(new_var, FC_NULL_VARIABLE), 
	      "fail should return FC_NULL_VARIABLE");
  // should fail if components don't match - different datatypes
  temp_comps[0] = scalar_var1;
  temp_comps[1] = scalar_var2;
  rc = fc_mergeComponentVariables(2, temp_comps, new_name, new_mathtype,
				  &new_var);
  fail_unless(rc != FC_SUCCESS, "should fail to merge diff assoc");
  fail_unless(FC_HANDLE_EQUIV(new_var, FC_NULL_VARIABLE), 
	      "fail should return FC_NULL_VARIABLE");
  // should fail if components don't match - different assocs
  temp_comps[0] = scalar_var1;
  temp_comps[1] = scalar_var3;
  rc = fc_mergeComponentVariables(2, temp_comps, new_name, new_mathtype,
				  &new_var);
  fail_unless(rc != FC_SUCCESS, "should fail to merge diff assoc");
  fail_unless(FC_HANDLE_EQUIV(new_var, FC_NULL_VARIABLE), 
	      "fail should return FC_NULL_VARIABLE");
  // should fail if they aren't on same mesh
  temp_comps[0] = scalar_var3;
  temp_comps[1] = scalar_var4;
  rc = fc_mergeComponentVariables(2, temp_comps, new_name, new_mathtype,
				  &new_var);
  fail_unless(rc != FC_SUCCESS, "should fail if not on same mesh");
  fail_unless(FC_HANDLE_EQUIV(new_var, FC_NULL_VARIABLE), 
	      "fail should return FC_NULL_VARIABLE");
  // merge components
  rc = fc_mergeComponentVariables(numComp, comps, new_name, new_mathtype,
				  &new_var);
  fail_unless(rc == FC_SUCCESS, "failed to merge components");
  // check merged variable
  rc = fc_getVariableName(new_var, &temp_name);
  fail_unless(rc == FC_SUCCESS, "failed to get name of merged var");
  fail_unless(!strcmp(new_name, temp_name), "mismatch of name");
  free(temp_name);
  rc = fc_getVariableInfo(new_var, &temp_numDataPoint, &temp_numComp,
			    &temp_assoc, &temp_mathtype, &temp_datatype);
  fail_unless(rc == FC_SUCCESS, "failed to describe merged variable");
  fail_unless(temp_numDataPoint == numDataPoint, "mismatch of numDataPoint");
  fail_unless(temp_numComp == numComp, "mismatch of numComp");
  fail_unless(temp_assoc == assoc, "mismatch of assoc");
  fail_unless(temp_mathtype == new_mathtype, "mismatch of mathtype");
  fail_unless(temp_datatype == datatype, "mismatch of datatype");
  rc = fc_getVariableDataPtr(new_var, (void**)&temp_data);
  fail_unless(rc == FC_SUCCESS, "failed to get data from merged var");
  fail_unless(!memcmp(data, temp_data, numDataPoint*numComp*fc_sizeofDataType(datatype)), 
	      "mismatch in data");

  // all done
  free(old_name);
  free(comps);
}
END_TEST	

// split variable into components
// essentially a copy of split_merge_comps
START_TEST(glb_split_merge_comps)
{
  FC_ReturnCode rc;
  int i, j;
  FC_Dataset dataset1, dataset2, *returnDatasets;
  int numReturnDatasets;
  FC_Variable vector_var;
  FC_Variable scalar_var1, scalar_var2, scalar_var3;
  FC_Variable *comps, temp_comps[2], new_var;
  int numDataPoint = 1, numComp = 3, temp_numDataPoint, temp_numComp;
  FC_AssociationType assoc = FC_AT_WHOLE_DATASET, temp_assoc;
  FC_MathType mathtype = FC_MT_SYMTENSOR, temp_mathtype, new_mathtype = FC_MT_VECTOR;
  FC_DataType datatype = FC_DT_DOUBLE, temp_datatype;
  double data[3] = { 1, 2, 3 }, *temp_data;
  char new_name[] = "hopeful", *old_name, *temp_name;
  char compNames[3][20] = { "_c0", "_c1", "_c2" }; 
  int name_len;
  double dbl_data[1] = { 0 }; // this data isn't really used
  float  flt_data[1] = { 0 }; // this data isn't really used

  // get dataset
  rc = fc_getDatasetByName(dataset_name, &numReturnDatasets,
			   &returnDatasets);
  fail_unless(rc == FC_SUCCESS, "aborted - cant get dataset by name");
  fail_unless(numReturnDatasets == 1,
	      "failed to find unique dataset by name");
  dataset1 = returnDatasets[0];
  free(returnDatasets);
  rc = fc_createDataset("another dataset", &dataset2);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create dataset");

  // create some vars
  // vector var is of type double
  rc = fc_createGlobalVariable(dataset1, "vector var", &vector_var);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create vector var");
  rc = fc_setVariableData(vector_var, 1, numComp, assoc, mathtype, 
			  datatype, data);
  fail_unless(rc == FC_SUCCESS, "abort: failed to set vector var data");
  fc_getVariableName(vector_var, &old_name);
  name_len = strlen(old_name);
  // scalar var1 is of type double
  rc = fc_createGlobalVariable(dataset1, "scalar var", &scalar_var1);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create var");
  rc = fc_setVariableData(scalar_var1, 1, 1, assoc, FC_MT_SCALAR, 
			  FC_DT_DOUBLE, dbl_data);
  fail_unless(rc == FC_SUCCESS, "abort: failed to set var data");
  // scalar var2 is of type float
  rc = fc_createGlobalVariable(dataset1, "scalar var float",
			       &scalar_var2); 
  fail_unless(rc == FC_SUCCESS, "abort: failed to create var");
  rc = fc_setVariableData(scalar_var2, 1, 1, assoc, FC_MT_SCALAR, FC_DT_FLOAT, 
			  flt_data);
  fail_unless(rc == FC_SUCCESS, "abort: failed to set var data");
  // scalar var3 is on dataset2
  rc = fc_createGlobalVariable(dataset2, "scalar var on other dataset",
			       &scalar_var3); 
  fail_unless(rc == FC_SUCCESS, "abort: failed to create var");
  rc = fc_setVariableData(scalar_var3, 1, 1, assoc, FC_MT_SCALAR, FC_DT_DOUBLE,
			  dbl_data);
  fail_unless(rc == FC_SUCCESS, "abort: failed to set var data");

  // --- test fc_createComponentVariables()
  // should fail to work on scalar variable
  rc = fc_createComponentVariables(scalar_var1, &temp_numComp, &comps);
  fail_unless(rc != FC_SUCCESS, "should not be able to get comps of scalar var");
  fail_unless(temp_numComp == -1, "failure should return -1 numComp");
  fail_unless(comps == NULL, "failure should return comps = NULL");
  // should will work on vector variable
  rc = fc_createComponentVariables(vector_var, &temp_numComp, &comps);
  fail_unless(rc == FC_SUCCESS, "fc_createComponentVariables failed");
  fail_unless(temp_numComp == numComp, "mismatch in number of components");
  // check components
  for (i = 0; i < numComp; i++) {
    fail_unless(fc_isVariableGlobal(comps[i]), "should be a global var");
    rc = fc_getVariableName(comps[i], &temp_name);
    fail_unless(rc == FC_SUCCESS, "failed to get name of component var");
    fail_unless(!strncmp(temp_name, old_name, name_len), 
		"mismatch of root name");
    fail_unless(!strcmp(&temp_name[name_len], compNames[i]), 
		"mismatch of name postfix");
    free(temp_name);
    rc = fc_getVariableInfo(comps[i], &temp_numDataPoint, &temp_numComp,
			      &temp_assoc, &temp_mathtype, &temp_datatype);
    fail_unless(rc == FC_SUCCESS, "failed to describe a component");
    fail_unless(temp_numDataPoint == numDataPoint, "mismatch of numDataPoint");
    fail_unless(temp_numComp == 1, "numComp should be 1 for component");
    fail_unless(temp_assoc == assoc, "mismatch of assoc");
    fail_unless(temp_mathtype == FC_MT_SCALAR, "mathtype should be SCALAR");
    fail_unless(temp_datatype == datatype, "mismatch of datatype");
    rc = fc_getVariableDataPtr(comps[i], (void**)&temp_data);
    fail_unless(rc == FC_SUCCESS, "failed to get data of a component");
    for (j = 0; j < numDataPoint; j++)
      fail_unless(data[j*numComp+i] == temp_data[j], "mismatch of data");
  }

  // --- test fc_mergeComponentVariables()
  // should fail if they aren't all scalars
  temp_comps[0] = scalar_var1;
  temp_comps[1] = vector_var;
  rc = fc_mergeComponentVariables(2, temp_comps, new_name, new_mathtype,
				  &new_var);
  fail_unless(rc != FC_SUCCESS, "should fail if not scalar");
  fail_unless(FC_HANDLE_EQUIV(new_var, FC_NULL_VARIABLE), 
	      "fail should return FC_NULL_VARIABLE");
  // should fail if components don't match - different datatypes
  temp_comps[0] = scalar_var1;
  temp_comps[1] = scalar_var2;
  rc = fc_mergeComponentVariables(2, temp_comps, new_name, new_mathtype,
				  &new_var);
  fail_unless(rc != FC_SUCCESS, "should fail to merge diff assoc");
  fail_unless(FC_HANDLE_EQUIV(new_var, FC_NULL_VARIABLE), 
	      "fail should return FC_NULL_VARIABLE");
  // should fail if they aren't on same dataset
  temp_comps[0] = scalar_var1;
  temp_comps[1] = scalar_var3;
  rc = fc_mergeComponentVariables(2, temp_comps, new_name, new_mathtype,
				  &new_var);
  fail_unless(rc != FC_SUCCESS, "should fail if not on same dataset");
  fail_unless(FC_HANDLE_EQUIV(new_var, FC_NULL_VARIABLE), 
	      "fail should return FC_NULL_VARIABLE");
  // merge components
  rc = fc_mergeComponentVariables(numComp, comps, new_name, new_mathtype,
				  &new_var);
  fail_unless(rc == FC_SUCCESS, "failed to merge components");
  // check merged variable
  fail_unless(fc_isVariableGlobal(new_var), "should be a global variable");
  rc = fc_getVariableName(new_var, &temp_name);
  fail_unless(rc == FC_SUCCESS, "failed to get name of merged var");
  fail_unless(!strcmp(new_name, temp_name), "mismatch of name");
  free(temp_name);
  rc = fc_getVariableInfo(new_var, &temp_numDataPoint, &temp_numComp,
			    &temp_assoc, &temp_mathtype, &temp_datatype);
  fail_unless(rc == FC_SUCCESS, "failed to describe merged variable");
  fail_unless(temp_numDataPoint == numDataPoint, "mismatch of numDataPoint");
  fail_unless(temp_numComp == numComp, "mismatch of numComp");
  fail_unless(temp_assoc == assoc, "mismatch of assoc");
  fail_unless(temp_mathtype == new_mathtype, "mismatch of mathtype");
  fail_unless(temp_datatype == datatype, "mismatch of datatype");
  rc = fc_getVariableDataPtr(new_var, (void**)&temp_data);
  fail_unless(rc == FC_SUCCESS, "failed to get data from merged var");
  fail_unless(!memcmp(data, temp_data, numDataPoint*numComp*fc_sizeofDataType(datatype)), 
	      "mismatch in data");

  // all done
  free(old_name);
  free(comps);
}
END_TEST	

// split seq variable into components
START_TEST(seq_split_merge_comps)
{
  FC_ReturnCode rc;
  int i, j, k;
  FC_Dataset dataset, *returnDatasets;
  int numReturnDatasets;
  FC_Sequence sequence1, sequence2, *returnSequences;
  int numStep = 5, temp_numStep, numReturnSequences;
  FC_Mesh mesh1, mesh2, *returnMeshes;
  int numReturnMeshes;
  FC_Variable coords_var, vector_vars[5], *vector_seqVar, *scalar_seqVar1;
  FC_Variable *scalar_seqVar2, *scalar_seqVar3, *scalar_seqVar4, *scalar_seqVar5;
  FC_Variable **compSeqVars, *temp_comps[2], *new_seqVar;
  int numDataPoint, numComp, temp_numDataPoint, temp_numComp;
  FC_AssociationType assoc, temp_assoc;
  FC_MathType mathtype, temp_mathtype, new_mathtype = FC_MT_VECTOR;
  FC_DataType datatype, temp_datatype;
  double* data, *temp_data;
  char new_name1[] = "hopeful", new_name2[] = "pessimistic", *temp_name;
  char compNames[3][20] = { "_c0", "_c1", "_c2" }; 
  int name_len;
  double dbl_data[50] = { 0 }; // this data isn't used
  float flt_data[50] = { 0 }; // this data isn't used
  int numVert, numElem;

  // get dataset sequences and meshes
  rc = fc_getDatasetByName(dataset_name, &numReturnDatasets,
			   &returnDatasets);
  fail_unless(rc == FC_SUCCESS, "aborted - cant get dataset by name");
  fail_unless(numReturnDatasets == 1,
	      "failed to find unique dataset by name");
  dataset = returnDatasets[0];
  free(returnDatasets);
  rc = fc_getSequenceByName(dataset,"time", &numReturnSequences,
			   &returnSequences);
  fail_unless(rc == FC_SUCCESS, "aborted - cant get sequence by name");
  fail_unless(numReturnSequences == 1,
	      "failed to find unique sequences by name");
  sequence1 = returnSequences[0];
  free(returnSequences);
  fc_getSequenceNumStep(sequence1, &temp_numStep);
  fail_unless(temp_numStep == numStep, "abort: mismatch of numStep");
  rc = fc_getSequenceByName(dataset,"time2", &numReturnSequences,
			   &returnSequences);
  fail_unless(rc == FC_SUCCESS, "aborted - cant get sequence by name");
  fail_unless(numReturnSequences == 1,
	      "failed to find unique sequences by name");
  sequence2 = returnSequences[0];
  free(returnSequences);
  rc = fc_getMeshByName(dataset, "hex mesh", &numReturnMeshes,
			&returnMeshes);
  fail_unless(rc == FC_SUCCESS, "aborted - cant get mesh by name");
  fail_unless(numReturnMeshes == 1, "failed to find unique mesh by name");
  mesh1 = returnMeshes[0];
  free(returnMeshes);
  rc = fc_getMeshInfo(mesh1, NULL, NULL, &numVert, &numElem, NULL);
  fail_unless(rc == FC_SUCCESS, "aborted: failed to get mesh info");
  rc = fc_getMeshByName(dataset, "tri mesh", &numReturnMeshes,
			&returnMeshes);
  fail_unless(rc == FC_SUCCESS, "aborted - cant get mesh by name");
  fail_unless(numReturnMeshes == 1, "failed to find unique mesh by name");
  mesh2 = returnMeshes[0];
  free(returnMeshes);
  // create a sequence vector var, data = stepID * coords;
  fc_getMeshCoordsAsVariable(mesh1, &coords_var);
  fc_getVariableInfo(coords_var, &numDataPoint, &numComp, &assoc, &mathtype,
		     &datatype);
  rc = fc_getVariableDataPtr(coords_var, (void**)&data);
  fail_unless(rc == FC_SUCCESS, "abort: failed to query coords variable");
  for (i = 0; i < numStep; i++) {
    int constant[numComp];
    for (j = 0; j < numComp; j++)
      constant[j] = i;
    fc_constOperatorVar(numComp, constant, FC_DT_INT, "*", coords_var, "temp", 
			&vector_vars[i]);
  } 
  name_len = strlen(new_name1);
  rc = fc_convertVariablesToSeqVariable(numStep, vector_vars, sequence1, 
					new_name1, &vector_seqVar);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create vector seq var");
  // scalar seqVar1 is on verts, of type double
  rc = fc_createSeqVariable(mesh1, sequence1, "scalar vert seqVar", &temp_numStep,
			    &scalar_seqVar1);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create scalar seqVar");
  for (i = 0; i < numStep; i++) {
    rc = fc_setVariableData(scalar_seqVar1[i], numVert, 1, FC_AT_VERTEX,
			    FC_MT_SCALAR, FC_DT_DOUBLE, dbl_data);
    fail_unless(rc == FC_SUCCESS, "abort: failed to set seqVar data");
  }
  // scalar seqVar2 is on verts, of type float
  rc = fc_createSeqVariable(mesh1, sequence1, "scalar vert seqVar float", 
			    &temp_numStep, &scalar_seqVar2);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create scalar seqVar");
  for (i = 0; i < numStep; i++) {
    rc = fc_setVariableData(scalar_seqVar2[i], numVert, 1, FC_AT_VERTEX,
			    FC_MT_SCALAR, FC_DT_FLOAT, flt_data);
    fail_unless(rc == FC_SUCCESS, "abort: failed to set seqVar data");
  }
  // scalar seqVar3 is on elems, of type double
  rc = fc_createSeqVariable(mesh1, sequence1, "scalar elem seqVar", 
			       &temp_numStep, &scalar_seqVar3);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create scalar seqVar");
  for (i = 0; i < numStep; i++) {
    rc = fc_setVariableData(scalar_seqVar3[i], numElem, 1, FC_AT_ELEMENT,
			    FC_MT_SCALAR, FC_DT_DOUBLE, dbl_data);
    fail_unless(rc == FC_SUCCESS, "abort: failed to set seqVar data");
  }
  // scalar seqVar4 is on elems on mesh2, of type double 
  // (is on elems because mesh1 and mesh 2 have same # of elems)
  rc = fc_createSeqVariable(mesh2, sequence1, "scalar vert seqVar on other mesh",
			    &temp_numStep, &scalar_seqVar4);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create scalar seqVar");
  for (i = 0; i < numStep; i++) {
    rc = fc_setVariableData(scalar_seqVar4[i], numElem, 1, FC_AT_ELEMENT,
			    FC_MT_SCALAR, FC_DT_DOUBLE, dbl_data);
    fail_unless(rc == FC_SUCCESS, "abort: failed to set seqVar data");
  }
  // scalar seqVar5 is on verts on sequence2, of type double
  rc = fc_createSeqVariable(mesh1, sequence2, "scalar vert seqVar on other seq",
			       &temp_numStep, &scalar_seqVar5);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create scalar seqVar");
  for (i = 0; i < numStep; i++) {
    rc = fc_setVariableData(scalar_seqVar5[i], numVert, 1, FC_AT_VERTEX,
			    FC_MT_SCALAR, FC_DT_DOUBLE, dbl_data);
    fail_unless(rc == FC_SUCCESS, "abort: failed to set seqVar data");
  }

  // --- test fc_createComponentSeqVariables()
  // should fail to work on a scalar variable
  rc = fc_createComponentSeqVariables(numStep, scalar_seqVar1, &temp_numComp, 
				      &compSeqVars);
  fail_unless(rc != FC_SUCCESS, "should not be able to get comps of scalar");
  fail_unless(temp_numComp == -1, "failure should return -1 numComp");
  fail_unless(compSeqVars == NULL, "failure should return compSeqVars = NULL");
  // should work on vector variable
  rc = fc_createComponentSeqVariables(numStep, vector_seqVar, &temp_numComp, 
				      &compSeqVars);
  fail_unless(rc == FC_SUCCESS, "fc_fetComponentVariables failed");
  fail_unless(temp_numComp == numComp, "mismatch in number of components");
  fail_unless(compSeqVars != NULL, "shouldn't be null");
  for (i = 0; i < numComp; i++)
    fail_unless(compSeqVars[i] != NULL, "shoudln't be null");
  // check components
  for (i = 0; i < numComp; i++) {
    for (j = 0; j < numStep; j++) {
      rc = fc_getVariableName(compSeqVars[i][j], &temp_name);
      fail_unless(rc == FC_SUCCESS, "failed to get name of component var");
      fail_unless(!strncmp(temp_name, new_name1, name_len), 
		  "mismatch of root name");
      fail_unless(!strcmp(&temp_name[name_len], compNames[i]), 
		  "mismatch of name postfix");
      free(temp_name);
      rc = fc_getVariableInfo(compSeqVars[i][j], &temp_numDataPoint, 
			      &temp_numComp, &temp_assoc, &temp_mathtype, 
			      &temp_datatype);
      fail_unless(rc == FC_SUCCESS, "failed to describe a component");
      fail_unless(temp_numDataPoint == numDataPoint, 
		  "mismatch of numDataPoint");
      fail_unless(temp_numComp == 1, "numComp should be 1 for component");
      fail_unless(temp_assoc == assoc, "mismatch of assoc");
      fail_unless(temp_mathtype == FC_MT_SCALAR, "mathtype should be SCALAR");
      fail_unless(temp_datatype == datatype, "mismatch of datatype");
      rc = fc_getVariableDataPtr(compSeqVars[i][j], (void**)&temp_data);
      fail_unless(rc == FC_SUCCESS, "failed to get data of a component");
      for (k = 0; k < numDataPoint; k++)
	fail_unless(FC_DBL_EQUIV(j*data[k*numComp+i], temp_data[k]), 
		    "mismatch of data");
    }
  }
  
  // --- test fc_mergeComponentSeqVariables()
  // should fail if they aren't all scalars
  temp_comps[0] = scalar_seqVar1;
  temp_comps[1] = vector_seqVar;
  rc = fc_mergeComponentSeqVariables(2, numStep, temp_comps, new_name2, 
				     new_mathtype, &new_seqVar);
  fail_unless(rc != FC_SUCCESS, "should fail to merge if not scalar");
  fail_unless(new_seqVar == NULL, "fail should return NULL");
  // should fail if components don't math - different assocs
  temp_comps[0] = scalar_seqVar1;
  temp_comps[1] = scalar_seqVar3;
  rc = fc_mergeComponentSeqVariables(2, numStep, temp_comps, new_name2, 
				     new_mathtype, &new_seqVar);
  fail_unless(rc != FC_SUCCESS, "should fail to merge diff assocs");
  fail_unless(new_seqVar == NULL, "fail should return NULL");
  // should fail if they aren't on same mesh
  temp_comps[0] = scalar_seqVar3;
  temp_comps[1] = scalar_seqVar4;
  rc = fc_mergeComponentSeqVariables(2, numStep, temp_comps, new_name2, 
				     new_mathtype, &new_seqVar);
  fail_unless(rc != FC_SUCCESS, "should fail to merge if not scalar");
  fail_unless(new_seqVar == NULL, "fail should return NULL");
  // should fail if they aren't on same sequence
  temp_comps[0] = scalar_seqVar1;
  temp_comps[1] = scalar_seqVar5;
  rc = fc_mergeComponentSeqVariables(2, numStep, temp_comps, new_name2, 
				     new_mathtype, &new_seqVar);
  fail_unless(rc != FC_SUCCESS, "should fail to merge if not scalar");
  fail_unless(new_seqVar == NULL, "fail should return NULL");
  
  // merge components
  rc = fc_mergeComponentSeqVariables(numComp, numStep, compSeqVars, new_name2, 
				     new_mathtype, &new_seqVar);
  fail_unless(rc == FC_SUCCESS, "failed to merge components");
  // check merged variable
  for (i = 0; i < numStep; i++) {
    rc = fc_getVariableName(new_seqVar[i], &temp_name);
    fail_unless(rc == FC_SUCCESS, "failed to get name of merged var");
    fail_unless(!strcmp(new_name2, temp_name), "mismatch of name");
    free(temp_name);
    rc = fc_getVariableInfo(new_seqVar[i], &temp_numDataPoint, &temp_numComp,
			    &temp_assoc, &temp_mathtype, &temp_datatype);
    fail_unless(rc == FC_SUCCESS, "failed to describe merged variable");
    fail_unless(temp_numDataPoint == numDataPoint, "mismatch of numDataPoint");
    fail_unless(temp_numComp == numComp, "mismatch of numComp");
    fail_unless(temp_assoc == assoc, "mismatch of assoc");
    fail_unless(temp_mathtype == new_mathtype, "mismatch of mathtype");
    fail_unless(temp_datatype == datatype, "mismatch of datatype");
    rc = fc_getVariableDataPtr(new_seqVar[i], (void**)&temp_data);
    fail_unless(rc == FC_SUCCESS, "failed to get data from merged var");
    for (j = 0; j < numDataPoint*numComp; j++)
      fail_unless(FC_DBL_EQUIV(i*data[j], temp_data[j]), "mismatch in data");
  }

  // all done
  free(vector_seqVar);
  free(scalar_seqVar1);
  free(scalar_seqVar2);
  free(scalar_seqVar3);
  free(scalar_seqVar4);
  free(scalar_seqVar5);
  for (i = 0; i < numComp; i++)
    free(compSeqVars[i]);
  free(compSeqVars);
  free(new_seqVar); 
}
END_TEST


// this must go after the merge test
START_TEST(comps_by_name)
{
  FC_ReturnCode rc;
  FC_Dataset dataset[2];
  FC_Mesh *returnMeshes, mesh[2], badMesh = {999,999};
  FC_MathType type;
  FC_Variable* vars = NULL;
  FC_Variable genvar;
  FC_Variable** seqvars = NULL;
  FC_Variable* genseqvar = NULL;
  int numReturnMeshes;
  char* names[4]= {"displacement", "velocity", "vertex ids", "stress tensor"};
  char vector_endings[3][2] = { "x", "y", "z" };
  int numRetVals[4] = {3,3,0,0};
  int numSeqRetVals[4] = {3,0,0,0};
  int numNames = 4;
  int numgenSteps;
  int *numSteps, numVars;
  int i, j;

  rc = fc_loadDataset("../data/gen_small_tri.ex2", &dataset[0]);
  fail_unless(rc == FC_SUCCESS, "abort: failed to load dataset");
  rc = fc_getMeshByName(dataset[0], "small_tri", &numReturnMeshes,&returnMeshes);
  fail_unless(rc == FC_SUCCESS, "aborted - cant get mesh by name");
  fail_unless(numReturnMeshes == 1, "failed to find unique mesh by name");
  mesh[0] = returnMeshes[0];
  free(returnMeshes);
  rc = fc_loadDataset("../data/gen_gaussians.ex2", &dataset[1]);
  fail_unless(rc == FC_SUCCESS, "abort: failed to load dataset");
  rc = fc_getMeshByName(dataset[1], "grid", &numReturnMeshes,&returnMeshes);
  fail_unless(rc == FC_SUCCESS, "aborted - cant get mesh by name");
  fail_unless(numReturnMeshes == 1, "failed to find unique mesh by name");
  mesh[1] = returnMeshes[0];
  free(returnMeshes);

  for (i = 0; i < numNames; i++){
    rc = fc_getVariableComponentsByName(mesh[0], names[i], &type, &numVars,
					&vars);
    fail_unless(rc == FC_SUCCESS, "failed to get var components");
    fail_unless(numVars == numRetVals[i], "wrong number of return vars");
    if (numVars > 0){
      fail_unless(type == FC_MT_VECTOR, "wrong math type");
    } else {
      fail_unless(type == FC_MT_UNKNOWN, "wrong math type");
    }

    if (numRetVals[i] > 0){
      fail_unless( vars != NULL, "should be return vars");
      char* compname = (char*)malloc((strlen(names[i])+3)*sizeof(char));
      for (j = 0; j < numVars; j++){
	char* retname;
	sprintf(compname, "%s_%s", names[i], vector_endings[j]);
	rc = fc_getVariableName(vars[j], &retname);
	fail_unless(strcmp(compname, retname)== 0 , "wrong var name");
	free(retname);
      }
      free(vars);

      // we can assume the vals for this works since its just a wrapper function
      rc = fc_getOrGenerateUniqueVariableByName(mesh[0], names[i], &genvar);
      fail_unless(rc == FC_SUCCESS, "failed to get var");
      fail_unless(!FC_HANDLE_EQUIV(genvar,FC_NULL_VARIABLE), "wrong return val");
      //make sure the prev vals are wiped out
      for (j = 0; j < numVars; j++){
	sprintf(compname, "%s_%s", names[i], vector_endings[j]);
	rc = fc_getVariableByName(mesh[0], compname, &numVars, &vars);
	fail_unless(numVars == 0, "should have wiped out vars");
      }
      free(compname);
    }
  }


  for (i = 0; i < numNames; i++){
    rc = fc_getSeqVariableComponentsByName(mesh[1], names[i], &type,
					   &numVars, &numSteps, &seqvars);
    fail_unless(rc == FC_SUCCESS, "failed to get var components");
    fail_unless(numVars == numSeqRetVals[i], "wrong number of return vars");
    if (numVars > 0){
      fail_unless(type == FC_MT_VECTOR, "wrong math type");
    } else {
      fail_unless(type == FC_MT_UNKNOWN, "wrong math type");
    }

    if (numVars > 0){
      fail_unless( seqvars != NULL, "should be return vars");
      char* compname = (char*)malloc((strlen(names[i])+3)*sizeof(char));
      for (j = 0; j < numVars; j++){
	char* retname;
	sprintf(compname, "%s_%s", names[i], vector_endings[j]);
	rc = fc_getVariableName(seqvars[j][0], &retname);
	fail_unless(strcmp(compname,retname) == 0, "wrong var name");
	if (seqvars[j]) free(seqvars[j]);
	free(retname);
      }
      free(seqvars);   
      free(numSteps);

      rc = fc_getOrGenerateUniqueSeqVariableByName(mesh[1], names[i],
						   &numgenSteps, &genseqvar);
      fail_unless(rc == FC_SUCCESS, "failed to get var");
      fail_unless( genseqvar != NULL, "wrong return val");
      //make sure the prev vals are wiped out
      for (j = 0; j < numVars; j++){
	sprintf(compname, "%s_%s", names[i], vector_endings[j]);
	rc = fc_getSeqVariableByName(mesh[1], compname, &numVars, 
				     &numSteps, &seqvars);
	fail_unless(numVars == 0, "should have wiped out vars");
      }
      free(compname);
      free(genseqvar);
    }
  }

  //non component tests - temp is nonseqvar in one file, seqvar in other
  rc = fc_getOrGenerateUniqueVariableByName(mesh[0], "temperature", &genvar);
  fail_unless(rc == FC_SUCCESS, "should work for noncomponents");
  rc = fc_getVariableByName(mesh[0], "temperature", &numVars, &vars);
  fail_unless(numVars == 1, "wrong number return vars");
  fail_unless(FC_HANDLE_EQUIV(vars[0], genvar), "wrong return var");
  free(vars);  
  rc = fc_getOrGenerateUniqueVariableByName(mesh[1], "temperature", &genvar);
  fail_unless(rc == FC_SUCCESS, "should work bad name");
  fail_unless(FC_HANDLE_EQUIV(genvar, FC_NULL_VARIABLE), "should return null var");

  rc = fc_getOrGenerateUniqueSeqVariableByName(mesh[0], "temperature", &numgenSteps,
					       &genseqvar);
  fail_unless(rc == FC_SUCCESS, "should work for bad name");
  fail_unless(genseqvar == NULL, "should return null var");
  fail_unless(numgenSteps == -1, "should return -1 steps");
  rc = fc_getOrGenerateUniqueSeqVariableByName(mesh[1], "temperature", &numgenSteps,
					       &genseqvar);
  fail_unless(rc == FC_SUCCESS, "should work for non components");
  rc = fc_getSeqVariableByName(mesh[1],"temperature", &numVars, &numSteps,
				    &seqvars);
  fail_unless(numVars == 1, "wrong number return vars");
  fail_unless(FC_HANDLE_EQUIV(seqvars[0][0], genseqvar[0]), "wrong return var");
  for (i = 0; i < numVars; i++){
    free(seqvars[i]);
  }
  free(numSteps);
  free(seqvars);
  free(genseqvar); 


  //bad args
  rc = fc_getVariableComponentsByName(FC_NULL_MESH, names[0], &type, &numVars, &vars);
  fail_unless(rc != FC_SUCCESS, "should fail for null mesh");
  rc = fc_getVariableComponentsByName(badMesh, names[0], &type, &numVars, &vars);
  fail_unless(rc != FC_SUCCESS, "should fail for bad mesh");
  rc = fc_getVariableComponentsByName(mesh[0], NULL, &type, &numVars, &vars);
  fail_unless(rc != FC_SUCCESS, "should fail for NULL name");
  rc = fc_getVariableComponentsByName(mesh[0], names[0], NULL, &numVars, &vars);
  fail_unless(rc != FC_SUCCESS, "should fail for NULL args");
  rc = fc_getVariableComponentsByName(mesh[0], names[0], &type, NULL, &vars);
  fail_unless(rc != FC_SUCCESS, "should fail for NULL args");
  rc = fc_getVariableComponentsByName(mesh[0], names[0], &type, &numVars, NULL);
  fail_unless(rc != FC_SUCCESS, "should fail for NULL args");

  rc = fc_getSeqVariableComponentsByName(FC_NULL_MESH, names[0], &type, &numVars,
					 &numSteps, &seqvars);
  fail_unless(rc != FC_SUCCESS, "should fail for null mesh");
  rc = fc_getSeqVariableComponentsByName(badMesh, names[0], &type, &numVars,
					 &numSteps, &seqvars);
  fail_unless(rc != FC_SUCCESS, "should fail for bad mesh");
  rc = fc_getSeqVariableComponentsByName(mesh[0], NULL, &type, &numVars, 
					 &numSteps, &seqvars);
  fail_unless(rc != FC_SUCCESS, "should fail for NULL name");
  rc = fc_getSeqVariableComponentsByName(mesh[0], names[0], NULL, &numVars,
					 &numSteps, &seqvars);
  fail_unless(rc != FC_SUCCESS, "should fail for NULL args");
  rc = fc_getSeqVariableComponentsByName(mesh[0], names[0], &type, NULL, 
					 &numSteps, &seqvars);
  fail_unless(rc != FC_SUCCESS, "should fail for NULL args");
  rc = fc_getSeqVariableComponentsByName(mesh[0], names[0], &type, &numVars,
					 NULL, &seqvars);
  fail_unless(rc != FC_SUCCESS, "should fail for NULL args");
  rc = fc_getSeqVariableComponentsByName(mesh[0], names[0], &type, &numVars,
					 &numSteps, NULL);
  fail_unless(rc != FC_SUCCESS, "should fail for NULL args");

  rc = fc_getOrGenerateUniqueVariableByName(FC_NULL_MESH, names[0], &genvar);
  fail_unless(rc != FC_SUCCESS, "should fail for null mesh");
  rc = fc_getOrGenerateUniqueVariableByName(badMesh, names[0], &genvar);
  fail_unless(rc != FC_SUCCESS, "should fail for bad mesh");
  rc = fc_getOrGenerateUniqueVariableByName(mesh[0], NULL, &genvar);
  fail_unless(rc != FC_SUCCESS, "should fail for NULL args");
  rc = fc_getOrGenerateUniqueVariableByName(mesh[0], names[0], NULL);
  fail_unless(rc != FC_SUCCESS, "should fail for NULL args");

  rc = fc_getOrGenerateUniqueSeqVariableByName(FC_NULL_MESH, names[0],
					       &numgenSteps, &genseqvar);
  fail_unless(rc != FC_SUCCESS, "should fail for null mesh");
  rc = fc_getOrGenerateUniqueSeqVariableByName(badMesh, names[0],
					       &numgenSteps, &genseqvar);
  fail_unless(rc != FC_SUCCESS, "should fail for bad mesh");
  rc = fc_getOrGenerateUniqueSeqVariableByName(mesh[0], NULL,
					       &numgenSteps, &genseqvar);
  fail_unless(rc != FC_SUCCESS, "should fail for NULL args");
  rc = fc_getOrGenerateUniqueSeqVariableByName(mesh[0], names[0],
					       NULL, &genseqvar);
  fail_unless(rc != FC_SUCCESS, "should fail for NULL args");
  rc = fc_getOrGenerateUniqueSeqVariableByName(mesh[0], names[0],
					       &numgenSteps, NULL);
  fail_unless(rc != FC_SUCCESS, "should fail for NULL args");
}
END_TEST
	
// split seq variable into components
// essentially a copy of seq_split_merge_comps
START_TEST(glb_seq_split_merge_comps)
{
  FC_ReturnCode rc;
  int i, j, k;
  FC_Dataset dataset, *returnDatasets;
  int numReturnDatasets;
  FC_Sequence sequence1, sequence2, *returnSequences;
  int numReturnSequences;
  int numStep = 5, temp_numStep;
  FC_Variable *vector_seqVar, *scalar_seqVar1, *scalar_seqVar2, *scalar_seqVar3;
  FC_Variable **compSeqVars, *temp_comps[2], *new_seqVar;
  int numDataPoint = 1, numComp = 3, temp_numDataPoint, temp_numComp;
  FC_AssociationType assoc = FC_AT_WHOLE_DATASET, temp_assoc;
  FC_MathType mathtype = FC_MT_SYMTENSOR, temp_mathtype, new_mathtype = FC_MT_VECTOR;
  FC_DataType datatype = FC_DT_DOUBLE, temp_datatype;
  double data[5][3] = { { 0, 1, 2 }, { 1, 2, 3 }, { 2, 3, 4 }, { 3, 4, 5 },
			 { 4, 5, 6 } };
  double *temp_data;
  char new_name1[] = "hopeful", new_name2[] = "pessimistic", *temp_name;
  char compNames[3][20] = { "_c0", "_c1", "_c2" }; 
  int name_len;
  double dbl_data[1] = { 0 }; // this data isn't used
  float flt_data[1] = { 0 }; // this data isn't used

  // get datasets and sequences
  rc = fc_getDatasetByName(dataset_name, &numReturnDatasets,
			   &returnDatasets);
  fail_unless(rc == FC_SUCCESS, "aborted - cant get dataset by name");
  fail_unless(numReturnDatasets == 1,
	      "failed to find unique dataset by name");
  dataset = returnDatasets[0];
  free(returnDatasets);
  rc = fc_getSequenceByName(dataset,"time", &numReturnSequences,
			   &returnSequences);
  fail_unless(rc == FC_SUCCESS, "aborted - cant get sequence by name");
  fail_unless(numReturnSequences == 1,
	      "failed to find unique sequences by name");
  sequence1 = returnSequences[0];
  free(returnSequences);
  fc_getSequenceNumStep(sequence1, &temp_numStep);
  fail_unless(temp_numStep == numStep, "abort: mismatch of numStep");
  rc = fc_getSequenceByName(dataset,"time2", &numReturnSequences,
			   &returnSequences);
  fail_unless(rc == FC_SUCCESS, "aborted - cant get sequence by name");
  fail_unless(numReturnSequences == 1,
	      "failed to find unique sequences by name");
  sequence2 = returnSequences[0];
  free(returnSequences);

  // create some vars
  // vector var is of type double
  rc = fc_createGlobalSeqVariable(dataset, sequence1, new_name1, 
				  &temp_numStep, &vector_seqVar);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create vector var");
  for (i = 0; i < numStep; i++) {
    rc = fc_setVariableData(vector_seqVar[i], 1, numComp, assoc, mathtype, 
			    datatype, data[i]);
    fail_unless(rc == FC_SUCCESS, "abort: failed to set vector var data");
  }
  name_len = strlen(new_name1);
  // scalar seqVar1 of type double
  rc = fc_createGlobalSeqVariable(dataset, sequence1, "scalar seqVar", 
				  &temp_numStep, &scalar_seqVar1);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create scalar seqVar");
  for (i = 0; i < numStep; i++) {
    rc = fc_setVariableData(scalar_seqVar1[i], 1, 1, assoc,
			    FC_MT_SCALAR, FC_DT_DOUBLE, dbl_data);
    fail_unless(rc == FC_SUCCESS, "abort: failed to set seqVar data");
  }
  // scalar seqVar2 is of type float
  rc = fc_createGlobalSeqVariable(dataset, sequence1, "scalar seqVar float", 
			    &temp_numStep, &scalar_seqVar2);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create scalar seqVar");
  for (i = 0; i < numStep; i++) {
    rc = fc_setVariableData(scalar_seqVar2[i], 1, 1, assoc,
			    FC_MT_SCALAR, FC_DT_FLOAT, flt_data);
    fail_unless(rc == FC_SUCCESS, "abort: failed to set seqVar data");
  }
  // scalar seqVar3 is on sequence2, of type double
  rc = fc_createGlobalSeqVariable(dataset, sequence2, 
				  "scalar seqVar on other seq",
			    &temp_numStep, &scalar_seqVar3);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create scalar seqVar");
  for (i = 0; i < numStep; i++) {
    rc = fc_setVariableData(scalar_seqVar3[i], 1, 1, assoc,
			    FC_MT_SCALAR, FC_DT_DOUBLE, dbl_data);
    fail_unless(rc == FC_SUCCESS, "abort: failed to set seqVar data");
  }

  // --- test fc_createComponentSeqVariables()
  // should fail to work on a scalar variable
  rc = fc_createComponentSeqVariables(numStep, scalar_seqVar1, &temp_numComp, 
				      &compSeqVars);
  fail_unless(rc != FC_SUCCESS, "should not be able to get comps of scalar");
  fail_unless(temp_numComp == -1, "failure should return -1 numComp");
  fail_unless(compSeqVars == NULL, "failure should return compSeqVars = NULL");
  // should work on vector variable
  rc = fc_createComponentSeqVariables(numStep, vector_seqVar, &temp_numComp, 
				      &compSeqVars);
  fail_unless(rc == FC_SUCCESS, "fc_fetComponentVariables failed");
  fail_unless(temp_numComp == numComp, "mismatch in number of components");
  fail_unless(compSeqVars != NULL, "shouldn't be null");
  for (i = 0; i < numComp; i++)
    fail_unless(compSeqVars[i] != NULL, "shoudln't be null");
  // check components
  for (i = 0; i < numComp; i++) {
    for (j = 0; j < numStep; j++) {
      rc = fc_getVariableName(compSeqVars[i][j], &temp_name);
      fail_unless(rc == FC_SUCCESS, "failed to get name of component var");
      fail_unless(!strncmp(temp_name, new_name1, name_len), 
		  "mismatch of root name");
      fail_unless(!strcmp(&temp_name[name_len], compNames[i]), 
		  "mismatch of name postfix");
      free(temp_name);
      rc = fc_getVariableInfo(compSeqVars[i][j], &temp_numDataPoint, 
			      &temp_numComp, &temp_assoc, &temp_mathtype, 
			      &temp_datatype);
      fail_unless(rc == FC_SUCCESS, "failed to describe a component");
      fail_unless(temp_numDataPoint == numDataPoint, 
		  "mismatch of numDataPoint");
      fail_unless(temp_numComp == 1, "numComp should be 1 for component");
      fail_unless(temp_assoc == assoc, "mismatch of assoc");
      fail_unless(temp_mathtype == FC_MT_SCALAR, "mathtype should be SCALAR");
      fail_unless(temp_datatype == datatype, "mismatch of datatype");
      rc = fc_getVariableDataPtr(compSeqVars[i][j], (void**)&temp_data);
      fail_unless(rc == FC_SUCCESS, "failed to get data of a component");
      for (k = 0; k < numDataPoint; k++)
	fail_unless(FC_DBL_EQUIV(data[j][k*numComp+i], temp_data[k]), 
		    "mismatch of data");
    }
  }

  // --- test fc_mergeComponentSeqVariables()
  // should fail if they aren't all scalars
  temp_comps[0] = scalar_seqVar1;
  temp_comps[1] = vector_seqVar;
  rc = fc_mergeComponentSeqVariables(2, numStep, temp_comps, new_name2, 
				     new_mathtype, &new_seqVar);
  fail_unless(rc != FC_SUCCESS, "should fail to merge if not scalar");
  fail_unless(new_seqVar == NULL, "fail should return NULL");
  // should fail if components don't math - different assocs
  temp_comps[0] = scalar_seqVar1;
  temp_comps[1] = scalar_seqVar2;
  rc = fc_mergeComponentSeqVariables(2, numStep, temp_comps, new_name2, 
				     new_mathtype, &new_seqVar);
  fail_unless(rc != FC_SUCCESS, "should fail to merge diff assocs");
  fail_unless(new_seqVar == NULL, "fail should return NULL");
  // should fail if they aren't on same sequence
  temp_comps[0] = scalar_seqVar1;
  temp_comps[1] = scalar_seqVar3;
  rc = fc_mergeComponentSeqVariables(2, numStep, temp_comps, new_name2, 
				     new_mathtype, &new_seqVar);
  fail_unless(rc != FC_SUCCESS, "should fail to merge if not scalar");
  fail_unless(new_seqVar == NULL, "fail should return NULL");
  
  // merge components
  rc = fc_mergeComponentSeqVariables(numComp, numStep, compSeqVars, new_name2, 
				     new_mathtype, &new_seqVar);
  fail_unless(rc == FC_SUCCESS, "failed to merge components");
  // check merged variable
  for (i = 0; i < numStep; i++) {
    rc = fc_getVariableName(new_seqVar[i], &temp_name);
    fail_unless(rc == FC_SUCCESS, "failed to get name of merged var");
    fail_unless(!strcmp(new_name2, temp_name), "mismatch of name");
    free(temp_name);
    rc = fc_getVariableInfo(new_seqVar[i], &temp_numDataPoint, &temp_numComp,
			    &temp_assoc, &temp_mathtype, &temp_datatype);
    fail_unless(rc == FC_SUCCESS, "failed to describe merged variable");
    fail_unless(temp_numDataPoint == numDataPoint, "mismatch of numDataPoint");
    fail_unless(temp_numComp == numComp, "mismatch of numComp");
    fail_unless(temp_assoc == assoc, "mismatch of assoc");
    fail_unless(temp_mathtype == new_mathtype, "mismatch of mathtype");
    fail_unless(temp_datatype == datatype, "mismatch of datatype");
    rc = fc_getVariableDataPtr(new_seqVar[i], (void**)&temp_data);
    fail_unless(rc == FC_SUCCESS, "failed to get data from merged var");
    for (j = 0; j < numDataPoint*numComp; j++)
      fail_unless(FC_DBL_EQUIV(data[i][j], temp_data[j]), "mismatch in data");
  }

  // all done
  free(vector_seqVar);
  free(scalar_seqVar1);
  free(scalar_seqVar2);
  free(scalar_seqVar3);
  for (i = 0; i < numComp; i++)
    free(compSeqVars[i]);
  free(compSeqVars);
  free(new_seqVar); 
}
END_TEST	



//Association type conversion with variables
START_TEST(vars_copy_with_association)
{
  int    numComponent = 3;
  int    MAX_DP       = 128;

  FC_ReturnCode rc, rc2;
  int i, j, k, l;
  int numReturnDatasets;
  int numMesh;

  double *vvals;
  FC_Variable *src_vars, *dst_vars;
  FC_Mesh *meshes;
  FC_Dataset dataset, *returnDatasets;

  double *src_data, *dst_data;
  double *res;
  int spot;


  FC_AssociationType src_assoc, dst_assoc, tst_type;
  int numDataPoint;
  int tst_num;
  char nameBuf[1024];
  char *tst_name;
  int  numObj, *objs;
  int numVertex, numElement;
  FC_Variable *var_list;


  // get dataset and meshes
  rc = fc_getDatasetByName(dataset_name, &numReturnDatasets,
			   &returnDatasets);
  fail_unless(rc == FC_SUCCESS, "aborted - cant get dataset by name");
  dataset = returnDatasets[0];
  free(returnDatasets);
  rc = fc_getMeshes(dataset, &numMesh, &meshes);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create mesh");

  //build some fake values to assign to vertices and elements
  vvals = (double *)malloc(MAX_DP*numComponent*sizeof(double));
  for(i=0; i<MAX_DP; i++)
    for(j=0; j<numComponent; j++)
      vvals[i*numComponent + j] = 1000.0 + (double)(100*i) + (double)j;


  //Allocate some space to hold on to a src and dst variable for each mesh
  src_vars = (FC_Variable *)malloc(numMesh      * sizeof(FC_Variable)); //orig vars we fill in
  dst_vars = (FC_Variable *)malloc(numMesh      * sizeof(FC_Variable)); //copy we create
  res      = (double *)     malloc(numComponent * sizeof(double));      //holds computed results
 

  //First test: convert every combination of associations on every input mesh
  //Outter loop: Try all source types
  for(src_assoc=FC_AT_VERTEX; src_assoc<=FC_AT_ELEMENT; src_assoc++){

    //printf("Starting with new source => %s\n", fc_getAssociationTypeText(src_assoc)); 

    //Create a source var for each mesh. Zero out the vars where the
    //format doesn't make sense.
    for(i=0; i<numMesh; i++){
      
      sprintf(nameBuf,"orig_%s", fc_getAssociationTypeText(src_assoc));
      
      switch(src_assoc){
      case FC_AT_VERTEX:
	rc = fc_getMeshNumVertex(meshes[i], &numDataPoint);
	fail_unless(rc == FC_SUCCESS, "Abort: failed to get mesh's edges");
	break;	  
      case FC_AT_EDGE:    
	rc = fc_getMeshNumEdge(meshes[i], &numDataPoint);
	fail_unless(rc == FC_SUCCESS, "Abort: failed to get mesh's edges");
	break;
      case FC_AT_FACE:   
	rc = fc_getMeshNumFace(meshes[i], &numDataPoint);
	fail_unless(rc == FC_SUCCESS, "Abort: failed to get mesh's faces");
	break;
      case FC_AT_ELEMENT:
	rc = fc_getMeshNumElement(meshes[i], &numDataPoint);
	fail_unless(rc == FC_SUCCESS, "Abort: failed to get mesh's edges");
	break;
      default:
	fail_unless(0, "abort: Unknown src_assoc in test..?");
      }
      
      if(numDataPoint==0){
	//Source format didn't make sense for this mesh.
	//printf("CREATE: SKIP  mesh %d due to src dims\n",i);
	src_vars[i] = FC_NULL_VARIABLE; //null out for later use
	continue; //Don't create the var
      } else {
	//printf("CREATE: doing mesh %d\n", i);
      }
      
      //Src Var ok. Try creating it
      rc = fc_createVariable(meshes[i], nameBuf, &src_vars[i]);
      fail_unless(rc == FC_SUCCESS, "abort: failed to create variable");
      rc = fc_setVariableData(src_vars[i], numDataPoint, numComponent, src_assoc,
			      FC_MT_VECTOR, FC_DT_DOUBLE,
			      (void *)vvals);
      fail_unless(rc == FC_SUCCESS, "abort: failed to set variable");
      
    }
    
    //Inner Loop: Now test against every type of conversion
    for(dst_assoc=FC_AT_VERTEX; dst_assoc<=FC_AT_ELEMENT; dst_assoc++){
      
      //printf("Working on type = %s ==> %s\n", 
      //       fc_getAssociationTypeText(src_assoc),
      //       fc_getAssociationTypeText(dst_assoc));
      
      //Do the copy operation
      for(i=0; i<numMesh; i++){

	//Skip the cases where the mesh dims prevented us from creating source
	if(FC_HANDLE_EQUIV(src_vars[i],FC_NULL_VARIABLE)) {
	  dst_vars[i] = FC_NULL_VARIABLE;
	  continue; 
	}
	sprintf(nameBuf,"copy_%s", fc_getAssociationTypeText(src_assoc));
	
	
	rc = fc_copyVariableWithNewAssociation(src_vars[i],
					       dst_assoc, nameBuf , &dst_vars[i]);
	//note: rc check done in the following switch block
	
	//In some meshes, output cannot be created because of dims
	switch(dst_assoc){
	case FC_AT_VERTEX:   
	  fail_unless(rc == FC_SUCCESS, "abort: unexpected failure for copyVarWithNewAssoc");
	  rc = fc_getMeshNumVertex(meshes[i], &numDataPoint);
	  fail_unless(rc == FC_SUCCESS, "Abort: failed to get mesh's edges");
	  break;
	  
	case FC_AT_EDGE:
	  rc2 = fc_getMeshNumEdge(meshes[i], &numDataPoint);
	  fail_unless(rc2 == FC_SUCCESS, "abort: couldn't get number of edges");
	  if(numDataPoint>0){
	    fail_unless(rc==FC_SUCCESS, "abort: copy failed when it should have passed");
	  } else {
	    fail_unless(rc!=FC_SUCCESS, "abort: copy passed when it should have failed");
	    fail_unless(FC_HANDLE_EQUIV(dst_vars[i],FC_NULL_VARIABLE), "operation failed ok, but didn't null the var");
	  }
	  break;
	  
	case FC_AT_FACE:
	  rc2 = fc_getMeshNumFace(meshes[i], &numDataPoint);
	  fail_unless(rc2 == FC_SUCCESS, "abort: couldn't get number of faces");
	  if(numDataPoint>0){
	    fail_unless(rc==FC_SUCCESS, "abort: copy failed when it should have passed");
	  } else {
	    fail_unless(rc!=FC_SUCCESS, "abort: copy passed when it should have failed");
	    fail_unless(FC_HANDLE_EQUIV(dst_vars[i],FC_NULL_VARIABLE), "operation failed ok, but didn't null the var");
	  }
	  break;
	  
	case FC_AT_ELEMENT:   
	  fail_unless(rc == FC_SUCCESS, "abort: unexpected failure for copyVarWithNewAssoc");	    
	  rc = fc_getMeshNumElement(meshes[i], &numDataPoint);
	  fail_unless(rc == FC_SUCCESS, "Abort: failed to get mesh's edges");	   
	  break;
	  
	default:
	  fail_unless(0, "whaa? unexpected dest type?");
	  break;
	}
	
	//Skip out if we didn't expect to get a result (and correctly got an error)
	if(numDataPoint<=0) continue;	  
	
	//Compare the var's num with expected num
	rc = fc_getVariableNumDataPoint(dst_vars[i], &tst_num);
	fail_unless(rc==FC_SUCCESS, "abort: couldn't get the var's number of points");
	fail_unless(numDataPoint == tst_num, "abort: var's number of points different than expected");
	
	//Make sure destination is what it's supposed to be
	rc = fc_getVariableAssociationType(dst_vars[i], &tst_type);
	fail_unless(rc==FC_SUCCESS, "couldn't get association type");
	fail_unless(tst_type == dst_assoc, "got wrong association type");

        //Check the name
	rc = fc_getVariableName(dst_vars[i], &tst_name);
	fail_unless(rc==FC_SUCCESS, "couldn't get variable name");
	fail_unless(!strcmp(tst_name,nameBuf), "new variable didn't have the expected name");
	free(tst_name);
	
	//Get a pointer to the created data
	rc  = fc_getVariableDataPtr(src_vars[i], (void **) &src_data);
	rc2 = fc_getVariableDataPtr(dst_vars[i], (void **) &dst_data);
	fail_unless(rc ==FC_SUCCESS, "couldn't get source data pointer");
	fail_unless(rc2==FC_SUCCESS, "couldn't get dest data pointer");
	
	//printf("==> Mesh %d  (%s -> %s)\n",i,
	//	 fc_getAssociationTypeText(src_assoc),
	//	 fc_getAssociationTypeText(dst_assoc));
	//fc_printMesh(meshes[i],"Mesh", 
	//	       ((src_assoc==FC_AT_EDGE)||(dst_assoc==FC_AT_EDGE)),
	//	       ((src_assoc==FC_AT_FACE)||(dst_assoc==FC_AT_FACE)),
	//	       ((src_assoc==FC_AT_ELEMENT)||(dst_assoc==FC_AT_ELEMENT)),
	//	       0); //Coords
	//fc_printVariable(src_vars[i], "src_var",1);
	//fc_printVariable(dst_vars[i], "dst_var",1);
	
	//Check equal src/dst types
	if(src_assoc == dst_assoc){
	  for(j=0;j<numDataPoint*numComponent;j++)
	    fail_unless(src_data[j] == dst_data[j], "data different in vars that should have been copied");
	  //printf("Checked same ok...\n");
	  continue;
	}
	
	//A conversion took place, so we need to check all the data
	for(j=0;j<numDataPoint; j++){
	  //find all source that contain this particular destination
	  if(src_assoc < dst_assoc)
	    rc = fc_getMeshEntityChildren(meshes[i], dst_assoc, j, src_assoc, &numObj, &objs);
	  else
	    rc = fc_getMeshEntityParents( meshes[i], dst_assoc, j, src_assoc, &numObj, &objs);
	  fail_unless(rc == FC_SUCCESS, "Could not get parent/child");
	  
	  //Add them up
	  memset(res, 0, numComponent*sizeof(double));
	  for(k=0; k<numComponent; k++){
	    for(l=0; l<numObj; l++){
	      spot = (objs[l] * numComponent) + k;
	      //res[k] += src_data[spot]; //if you compare to src_var's values
	      res[k] += vvals[spot];      //if you compare to original data
	    }
	  }
	  if(numObj)
	    for(k=0;k<numComponent;k++)
	      res[k] /= (double)numObj;
	  
	  for(k=0;k<numComponent; k++)
	    fail_unless(res[k] == dst_data[j*numComponent + k], "result data doesn't look right");
	  
	  free(objs);
	}

	//No longer need this var
	fc_deleteVariable(dst_vars[i]);

      }
    } //All dst's 
    //Cleanup the sources
    for(i=0;i<numMesh;i++)
      fc_deleteVariable(src_vars[i]);

  }//All src's

  //Basic input checking

  //Null input
  src_vars[0] = FC_NULL_VARIABLE;
  rc = fc_copyVariableWithNewAssociation(src_vars[0], FC_AT_VERTEX, "Booya", &dst_vars[0]);
  fail_unless(rc != FC_SUCCESS, "failed null input variable test");
  fail_unless(FC_HANDLE_EQUIV(dst_vars[0], FC_NULL_VARIABLE), "failed to set null variable when input is null");

  //Next couple of tests work on mesh[3] and need a source var
  rc = fc_getMeshInfo(meshes[3], NULL, NULL, &numVertex, &numElement, NULL);
  fail_unless(rc == FC_SUCCESS, "couldn't get mesh info");
  rc = fc_createVariable(meshes[3], "shammalangadingdong", &src_vars[3]);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create variable");
  rc = fc_setVariableData(src_vars[3], numVertex, numComponent, FC_AT_VERTEX,
			  FC_MT_VECTOR, FC_DT_DOUBLE,
			  (void *)vvals);
  fail_unless(rc == FC_SUCCESS, "abort: failed to set variable");


  //Null output test  
  rc = fc_copyVariableWithNewAssociation(src_vars[3], FC_AT_ELEMENT, "Booya", NULL);
  fail_unless(rc != FC_SUCCESS, "failed null output variable test");


  //Null name: create and query new variable
  rc = fc_copyVariableWithNewAssociation(src_vars[3], FC_AT_ELEMENT, NULL, &dst_vars[3]);
  fail_unless(rc == FC_SUCCESS, "failed null variable name test");
  rc = fc_getVariableName(dst_vars[3], &tst_name);
  fail_unless(rc==FC_SUCCESS, "couldn't get variable name");
  fail_unless(!strcmp(tst_name, "shammalangadingdong"), "new variable didn't have the expected name");
  free(tst_name);

  //Null name: look at mesh and make sure ok
  rc = fc_getVariableByName(meshes[3], "shammalangadingdong", &tst_num, &var_list);
  fail_unless(rc==FC_SUCCESS, "couldn't get variable by name");
  fail_unless(tst_num==2, "Was expecting only two versions of a named variable");
  free(var_list);
  fc_deleteVariable(dst_vars[3]); //get rid of copy


  //Copy with a new name
  rc = fc_copyVariableWithNewAssociation(src_vars[3], FC_AT_ELEMENT, "raskolnikov", &dst_vars[3]);
  fail_unless(rc == FC_SUCCESS, "failed to copy variable");
  rc = fc_getVariableByName(meshes[3], "raskolnikov", &tst_num, &var_list);
  fail_unless(tst_num==1, "Was expecting only one versions of a named variable");
  free(var_list);
  fc_deleteVariable(dst_vars[3]); //get rid of copy
  


  //Global Variable should be weeded out
  rc = fc_createGlobalVariable(dataset, "globbo", &src_vars[0]);
  fail_unless(rc==FC_SUCCESS,"Couldn't create a global variable");
  rc = fc_setVariableData(src_vars[0], 1, numComponent, FC_AT_WHOLE_DATASET,
			  FC_MT_VECTOR, FC_DT_DOUBLE,
			  (void *)vvals);
  fail_unless(rc==FC_SUCCESS,"Couldn't set a global variable");
  rc = fc_copyVariableWithNewAssociation(src_vars[0], FC_AT_ELEMENT, "Booya", &dst_vars[1]);
  fail_unless(rc != FC_SUCCESS, "Global var copy worked when it should have failed");

  fc_deleteVariable(src_vars[0]); //get rid of global


  //Cleanup everything
  free(res);
  free(dst_vars);
  free(src_vars);
  free(vvals);
  free(meshes);



  //Ok finally, lets test the calls with a simple dataset that's been
  //processed by hand offline.
  {
    double *co;
    FC_Mesh m;
    FC_Variable var_result, var_verts;
    double *tst_data;

    //We're building a simple 3x3x2 structure here with 4 cube elements
    //
    //  Vertex are 00-17  and elements A-D
    //
    // z=0:          z=1:
    //     00--01--02       09--10--11         
    //      | A | B |        | A | B |
    //     03--04--05       12--13--14 
    //      | C | D |        | C | D |
    //     06--07--08       15--16--17
    


    int newconns[] = { 0,1,4,3,9,10,13,12,  
		       1,2,5,4,10,11,14,13,
		       3,4,7,6,12,13,16,15,
		       4,5,8,7,13,14,17,16 };

    //First test is vertex->element
    double exp_elem[] = { (double)(0+1+4+3+9+10+13+12)/8.0,
			  (double)(1+2+5+4+10+11+14+13)/8.0,
			  (double)(3+4+7+6+12+13+16+15)/8.0,
			  (double)(4+5+8+7+13+14+17+16)/8.0 };
    int exp_elem_num = 4;

    //Second test is vertex->face.
    // Note: we build the face list by extracting all faces from each element
    //       and then removing redundant faces.
    //
    // This list is built cube-by-cube, and then redundant faces are commented out
    double exp_face[] = { //First cube
                          (double)(0+1+9+10)/4.0,
			  (double)(1+4+10+13)/4.0,
			  (double)(3+4+12+13)/4.0,
			  (double)(0+3+9+12)/4.0,
			  (double)(0+1+4+3)/4.0,
			  (double)(9+10+13+12)/4.0,
			  //Second cube
			  (double)(1+2+10+11)/4.0,
			  (double)(2+5+11+14)/4.0,
			  (double)(4+5+13+14)/4.0,
			  //(double)(4+1+13+10)/4.0, //Remove: same as second in first cube
			  (double)(1+2+4+5)/4.0,
			  (double)(10+11+13+14)/4.0,
			  //Third cube
			  //(double)(3+4+12+13)/4.0, //Remove: same as third in first cube
			  (double)(4+7+13+16)/4.0,
			  (double)(6+7+15+16)/4.0,
			  (double)(3+6+12+15)/4.0,
			  (double)(3+4+7+6)/4.0,
			  (double)(12+13+16+15)/4.0,
			  //Fourth cube
			  //(double)(4+5+13+14)/4.0, //Remove: same as third in second cube
			  (double)(5+8+14+17)/4.0,
			  (double)(8+7+16+17)/4.0,
			  //(double)(7+4+16+13)/4.0, //Remove: same as first in third cube
			  (double)(4+5+8+7)/4.0,
			  (double)(13+14+17+16)/4.0};


    int exp_face_num = 4*6 - 4;


    double exp_edge[] =  { //First cube
                           //Face 0 (0,1,10,9): (do all)
                           (double) (0+1)/2.0,
			   (double) (1+10)/2.0,
			   (double) (9+10)/2.0,
			   (double) (0+9)/2.0,
			   //Face 1 (1,4,13,10): Already did 10-1
			   (double) (1+4)/2.0,
			   (double) (4+13)/2.0,
			   (double) (13+10)/2.0,	   
			   //Face 2 (4,3,12,13): Already did 13-4
			   (double) (4+3)/2.0,
			   (double) (3+12)/2.0,
			   (double) (12+13)/2.0,
			   //Face 3 (3,0,9,12): already did 0-9,12-3 
			   (double) (3+0)/2.0,			  
			   (double) (9+12)/2.0,
			   //Face 4,5: Top/bottom faces already done

			   //Second cube
			   //Face 0 (1,2,11,10): already did 1-10
			   (double) (1+2)/2.0,
			   (double) (2+11)/2.0,
			   (double) (11+10)/2.0,
			   //Face 1 (2,5,14,11): already did 11-2
			   (double) (2+5)/2.0,
			   (double) (5+14)/2.0,
			   (double) (14+11)/2.0,
			   //Face 2 (5,4,13,14): already did 4-13,14-5
			   (double) (5+4)/2.0,
			   (double) (13+14)/2.0,
			   //Face 3 (4,1,10,13) same as cube 0, face 1
			   //Face 4,5: Top/bottom faces already done

			   //Third cube
			   //Face 0 (3,4,13,12): same as cube 0 face 2
			   //Face 1 (4,7,16,13): already did 4-7
			   (double) (4+7)/2.0,
			   (double) (7+16)/2.0,
			   (double) (16+13)/2.0,
			   //Face 2 (7,6,15,16): already did 16-7
			   (double) (7+6)/2.0,
			   (double) (6+15)/2.0,
			   (double) (15+16)/2.0,
			   //Face 3 (6,3,12,15): already did 3-12,15-6
			   (double) (6+3)/2.0,
			   (double) (12+15)/2.0,
			   //Face 4,5: top/bottom faces already done

			   //Fourth cube
			   //Face 0 (4,5,14,13): same as cube 1 face 2
			   //Face 1 (5,8,17,14): already did 5-14
			   (double) (5+8)/2.0,
			   (double) (8+17)/2.0,
			   (double) (17+14)/2.0,
			   //Face 2 (8,7,16,17): already did 8-17, 7-16
			   (double) (8+7)/2.0,
			   (double) (16+17)/2.0
			   //Face 3 (7,4,13,16) same as cube 2 face 1
			   //Face 4,5: Top/bottom faces already done
                         };

    int exp_edge_num = 12*2+9; /* 12 edges per z-plane, 9 edges to connect planes */

    typedef struct {
      FC_AssociationType assocType;
      int                num;
      double            *vals;
      char              *test_name;
    } exp_t;


    //Stuff all the tests into a simple struct so it's easier to manage
    exp_t expected[] = { {FC_AT_ELEMENT, exp_elem_num, exp_elem, "vertex2element"},
                         {FC_AT_FACE,    exp_face_num, exp_face, "vertex2face"},
			 {FC_AT_EDGE,    exp_edge_num, exp_edge, "vertex2edge"} };
    int expected_num = 3;



    //Coords don't really matter
    co = (double *) malloc(3*18*sizeof(double));
    fail_unless(co!=NULL, "no malloc");
    for(i=0;i<18;i++)
      for(j=0;j<3;j++)
	co[i*3 + j] = (double) i;

    rc = fc_createMesh(dataset, "hexer mesher", &m);
    fail_unless(rc==FC_SUCCESS, "failed to create mesh");
    rc = fc_setMeshCoords(m, 3, 18, co);
    fail_unless(rc==FC_SUCCESS, "failed to set mesh");
    rc = fc_setMeshElementConns(m, FC_ET_HEX, 4, newconns);
    fail_unless(rc == FC_SUCCESS, "failed to set element conns");

    //Set each vertex to have it's id value 
    rc = fc_createVariable(m, "bosephus", &var_verts);
    fail_unless(rc == FC_SUCCESS, "abort: failed to create variable");
    rc = fc_setVariableData(var_verts, 18, 3, FC_AT_VERTEX,
			    FC_MT_VECTOR, FC_DT_DOUBLE,
			    (void *)co);
    fail_unless(rc == FC_SUCCESS, "abort: failed to set variable");

    //fc_printMesh(m, "mesh", 0,0,1,1);
    //fc_printVariable(var_verts, "vert data",1);


    //Try all the tests
    for(i=0;i<expected_num; i++){
      //printf("Trying test %s\n", expected[i].test_name);

      rc = fc_copyVariableWithNewAssociation(var_verts, expected[i].assocType, expected[i].test_name, &var_result);
      fail_unless(rc == FC_SUCCESS, "Couldn't copy variable to a particular association type");

      //fc_printVariable(var_elems, expected[i].test_name, 1);
      
      rc = fc_getVariableNumDataPoint(var_result, &tst_num);
      fail_unless(rc == FC_SUCCESS, "Failed to get data points");
      fail_unless(tst_num == expected[i].num, "Didn't get the expected number of values");
      
      rc = fc_getVariableDataPtr(var_result, (void **) &tst_data);
      fail_unless(rc == FC_SUCCESS, "Failed to get data pointer");
      
      for(j=0; j<expected[i].num; j++) {
	//printf("%s edge %d:: %.3lf %.3lf\n", expected[i].test_name, j, tst_data[j*numComponent + 0], expected[i].vals[j]);
	for(k=0; k<numComponent;k++){  
	  fail_unless( tst_data[j*numComponent + k]==expected[i].vals[j], "Converted data did not match expected value");
	}
      }
      
      //Done with particular test. Cleanup this round
      fc_deleteVariable(var_result);
    }

    //Free up stuff
    fc_deleteVariable(var_verts);
    fc_deleteMesh(m);
    free(co);

  } //End of canned test




}
END_TEST
  



//Association type conversion with variables
START_TEST(seqvars_copy_with_association)
{
  FC_ReturnCode rc, rc2;
  int i, j, k, l, s;
  int numReturnDatasets;
  int numMesh;

  double *vvals;
  int    numComponent = 3;
  int    MAX_DP = 128;
  
  FC_Mesh *meshes;
  FC_Dataset dataset, *returnDatasets;


  double *src_data, *dst_data;
  double *res;
  int spot;


  FC_AssociationType src_assoc, dst_assoc, tst_type;
  int numDataPoint;
  int tst_num;
  char nameBuf[1024];
  int  numObj, *objs;
  char *tst_name;
   

  //seqversion
  int numReturnSequences;
  FC_Sequence *returnSequences, seq;
  int numStep;

  FC_Variable **src_vars, **dst_vars;
  int          test_numStep;
  FC_Variable *test_gsv;



  // get dataset and meshes
  rc = fc_getDatasetByName(dataset_name, &numReturnDatasets,
			   &returnDatasets);
  fail_unless(rc == FC_SUCCESS, "aborted - cant get dataset by name");
  dataset = returnDatasets[0];
  free(returnDatasets);
  rc = fc_getMeshes(dataset, &numMesh, &meshes);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create mesh");
  
  rc = fc_getSequenceByName(dataset,"time", 
			    &numReturnSequences, &returnSequences);
  fail_unless(rc == FC_SUCCESS, "aborted - cant get sequence by name");
  fail_unless(numReturnSequences == 1,
	      "failed to find unique sequences by name");
  seq = returnSequences[0];
  free(returnSequences);
  rc = fc_getSequenceNumStep(seq, &numStep);
  fail_unless(rc == FC_SUCCESS, "Couldn't get numSteps");
  fail_unless(numStep>0, "number of steps was negative");

  //build some fake values to assign to vertices and elements
  vvals = (double *)malloc(numStep*MAX_DP*numComponent*sizeof(double));
  for(i=0; i<numStep; i++)
    for(j=0; j<MAX_DP; j++)
      for(k=0; k<numComponent; k++)
	vvals[i*MAX_DP*numComponent + j*numComponent + k] = 10000.0 + (double) ((1000*i) + (100*j) + k);


  //Allocate some space to hold on to a src and dst variable for each mesh
  src_vars = (FC_Variable **)malloc(numMesh      * sizeof(FC_Variable *)); //orig vars we fill in
  dst_vars = (FC_Variable **)malloc(numMesh      * sizeof(FC_Variable *)); //copy we create
  res      = (double *)      malloc(numComponent * sizeof(double));      //holds computed results
 

  //First test: convert every combination of associations on every input mesh
  //Outter loop: Try all source types
  for(src_assoc=FC_AT_VERTEX; src_assoc<=FC_AT_ELEMENT; src_assoc++){

    //printf("Starting with new source => %s\n", fc_getAssociationTypeText(src_assoc)); 

    //Create a source var for each mesh. Zero out the vars where the
    //format doesn't make sense.
    for(i=0; i<numMesh; i++){
      
      sprintf(nameBuf,"orig_%s", fc_getAssociationTypeText(src_assoc));
      
      switch(src_assoc){
      case FC_AT_VERTEX:
	rc = fc_getMeshNumVertex(meshes[i], &numDataPoint);
	fail_unless(rc == FC_SUCCESS, "Abort: failed to get mesh's edges");
	break;	  
      case FC_AT_EDGE:    
	rc = fc_getMeshNumEdge(meshes[i], &numDataPoint);
	fail_unless(rc == FC_SUCCESS, "Abort: failed to get mesh's edges");
	break;
      case FC_AT_FACE:   
	rc = fc_getMeshNumFace(meshes[i], &numDataPoint);
	fail_unless(rc == FC_SUCCESS, "Abort: failed to get mesh's faces");
	break;
      case FC_AT_ELEMENT:
	rc = fc_getMeshNumElement(meshes[i], &numDataPoint);
	fail_unless(rc == FC_SUCCESS, "Abort: failed to get mesh's edges");
	break;
      default:
	fail_unless(0, "abort: Unknown src_assoc in test..?");
      }
      
      if(numDataPoint==0){
	//Source format didn't make sense for this mesh.
	//printf("CREATE: SKIP  mesh %d due to src dims\n",i);
	src_vars[i] = NULL; //FC_NULL_VARIABLE; //null out for later use
	continue; //Don't create the var
      } else {
	//printf("CREATE: doing mesh %d\n", i);
      }
      
      //Src Var ok. Try creating it
      rc = fc_createSeqVariable(meshes[i], seq,  nameBuf, &test_numStep,  &src_vars[i]);
      fail_unless(rc == FC_SUCCESS, "abort: failed to create variable");
      fail_unless(test_numStep == numStep, "seqvar create got wrong number of steps");
      //Fill in all the steps
      for(j=0; j<numStep; j++){

	rc = fc_setVariableData((src_vars[i])[j], numDataPoint, numComponent, src_assoc,
				FC_MT_VECTOR, FC_DT_DOUBLE,
				(void *) &vvals[j*MAX_DP*numComponent]);
	fail_unless(rc == FC_SUCCESS, "abort: failed to set variable");

      }
      
    }
    
    //Inner Loop: Now test against every type of conversion
    for(dst_assoc=FC_AT_VERTEX; dst_assoc<=FC_AT_ELEMENT; dst_assoc++){
      
      //printf("Working on type = %s ==> %s\n", 
      //       fc_getAssociationTypeText(src_assoc),
      //       fc_getAssociationTypeText(dst_assoc));
      
      //Do the copy operation
      for(i=0; i<numMesh; i++){

	//Skip the cases where the mesh dims prevented us from creating source
	if(src_vars[i]==NULL) {
	  dst_vars[i] = NULL;
	  continue; 
	}
	sprintf(nameBuf,"copy_%s", fc_getAssociationTypeText(src_assoc));
	
	rc = fc_copySeqVariableWithNewAssociation(numStep, src_vars[i], 
						  dst_assoc, nameBuf ,
						  &dst_vars[i]);
	//note: rc check done in the following switch block
	
	//In some meshes, output cannot be created because of dims
	switch(dst_assoc){
	case FC_AT_VERTEX:   
	  fail_unless(rc == FC_SUCCESS, "abort: unexpected failure for copyVarWithNewAssoc");
	  rc = fc_getMeshNumVertex(meshes[i], &numDataPoint);
	  fail_unless(rc == FC_SUCCESS, "Abort: failed to get mesh's edges");
	  break;
	  
	case FC_AT_EDGE:
	  rc2 = fc_getMeshNumEdge(meshes[i], &numDataPoint);
	  fail_unless(rc2 == FC_SUCCESS, "abort: couldn't get number of edges");
	  if(numDataPoint>0){
	    fail_unless(rc==FC_SUCCESS, "abort: copy failed when it should have passed");
	  } else {
	    fail_unless(rc!=FC_SUCCESS, "abort: copy passed when it should have failed");
	    fail_unless(dst_vars[i]==NULL, "operation failed ok, but didn't null the var");
	  }
	  break;
	  
	case FC_AT_FACE:
	  rc2 = fc_getMeshNumFace(meshes[i], &numDataPoint);
	  fail_unless(rc2 == FC_SUCCESS, "abort: couldn't get number of faces");
	  if(numDataPoint>0){
	    fail_unless(rc==FC_SUCCESS, "abort: copy failed when it should have passed");
	  } else {
	    fail_unless(rc!=FC_SUCCESS, "abort: copy passed when it should have failed");
	    fail_unless(dst_vars[i]==NULL, "operation failed ok, but didn't null the var");
	  }
	  break;
	  
	case FC_AT_ELEMENT:   
	  fail_unless(rc == FC_SUCCESS, "abort: unexpected failure for copyVarWithNewAssoc");	    
	  rc = fc_getMeshNumElement(meshes[i], &numDataPoint);
	  fail_unless(rc == FC_SUCCESS, "Abort: failed to get mesh's edges");	   
	  break;
	  
	default:
	  fail_unless(0, "whaa? unexpected dest type?");
	  break;
	}
	
	//Skip out if we didn't expect to get a result
	if(numDataPoint<=0) 	  
	  continue;	  
	
	//Compare the var's num with expected num
	rc = fc_getVariableNumDataPoint(dst_vars[i][0], &tst_num);
	fail_unless(rc==FC_SUCCESS, "abort: couldn't get the var's number of points");
	fail_unless(numDataPoint == tst_num, "abort: var's number of points different than expected");
	
	//Make sure the name is what we expect
	rc = fc_getVariableName(dst_vars[i][0], &tst_name);
	fail_unless(rc==FC_SUCCESS, "couldn't get new variable name");
	fail_unless(!strcmp(tst_name, nameBuf), "new variable name didn't match");
	free(tst_name);

	//Make sure destination is what it's supposed to be
	rc = fc_getVariableAssociationType(dst_vars[i][0], &tst_type);
	fail_unless(rc==FC_SUCCESS, "couldn't get association type");
	fail_unless(tst_type == dst_assoc, "got wrong association type");
	
	//Check all steps
	for(s=0; s<numStep; s++){	  

	  //Get a pointer to the created data
	  rc  = fc_getVariableDataPtr(src_vars[i][s], (void **) &src_data);
	  rc2 = fc_getVariableDataPtr(dst_vars[i][s], (void **) &dst_data);
	  fail_unless(rc ==FC_SUCCESS, "couldn't get source data pointer");
	  fail_unless(rc2==FC_SUCCESS, "couldn't get dest data pointer");
	

	  //Check equal src/dst types
	  if(src_assoc == dst_assoc){
	    for(j=0;j<numDataPoint*numComponent;j++)
	      fail_unless(src_data[j] == dst_data[j], "data different in vars that should have been copied");
	    continue;
	  }
	
	  //A conversion took place, so we need to check all the data
	  for(j=0;j<numDataPoint; j++){
	    //find all source that contain this particular destination
	    if(src_assoc < dst_assoc)
	      rc = fc_getMeshEntityChildren(meshes[i], dst_assoc, j, src_assoc, &numObj, &objs);
	    else
	      rc = fc_getMeshEntityParents( meshes[i], dst_assoc, j, src_assoc, &numObj, &objs);
	    fail_unless(rc == FC_SUCCESS, "Could not get parent/child");
	  
	    //Add them up
	    memset(res, 0, numComponent*sizeof(double));
	    for(k=0; k<numComponent; k++){
	      for(l=0; l<numObj; l++){
		spot = (objs[l] * numComponent) + k;
		//res[k] += src_data[spot];                         //if you compare to src_var's values
		res[k] += vvals[(s*MAX_DP*numComponent) + spot];  //if you compare to original data
	      }
	    }
	    if(numObj)
	      for(k=0;k<numComponent;k++)
		res[k] /= (double)numObj;
	    
	    for(k=0;k<numComponent; k++){
	      //printf("Step %d: %.3f vs %.3f\n",s, res[k], dst_data[j*numComponent + k]);
	      fail_unless(res[k] == dst_data[j*numComponent + k], "result data doesn't look right");
	    }
	    free(objs);
	  }
	}//steps

	//Cleanup destination vars
	rc = fc_deleteSeqVariable(numStep, dst_vars[i]);
	fail_unless(rc==FC_SUCCESS, "failed to delete dest seqVar");
	free(dst_vars[i]);
	

      }//all meshes
    }//all dst's
    //Cleanup all sources
    for(i=0;i<numMesh;i++){
      if(src_vars[i]!=NULL){
	rc = fc_deleteSeqVariable(numStep, src_vars[i]);
	fail_unless(rc==FC_SUCCESS, "failed to delete src seqVar");
	free(src_vars[i]);
      }
    }
  }//all src's 
  //Done checking all combinations of source/destination associations


  //Basic input checking

  //Null input
  rc = fc_copySeqVariableWithNewAssociation(numStep, NULL, FC_AT_VERTEX, "Booya", &dst_vars[0]);
  fail_unless(rc != FC_SUCCESS, "failed null input variable test");
  fail_unless((dst_vars[0]==NULL), "failed to set null variable when input is null");

  //Setup a fake variable for next two tests
  rc = fc_createSeqVariable(meshes[0], seq,  "clover", &test_numStep,  &src_vars[0]);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create variable");
  for(j=0; j<test_numStep; j++){
    rc = fc_setVariableData(src_vars[0][j], numDataPoint, numComponent, FC_AT_VERTEX,
			    FC_MT_VECTOR, FC_DT_DOUBLE,
			    (void *) &vvals[j*MAX_DP*numComponent]);
    fail_unless(rc == FC_SUCCESS, "abort: failed to set variable");
  }

  //Null output variable
  rc = fc_copySeqVariableWithNewAssociation(numStep,src_vars[0], FC_AT_ELEMENT, "Booya", NULL);
  fail_unless(rc != FC_SUCCESS, "failed null output variable test");

  //Null name: should produce the same name as source
  rc = fc_copySeqVariableWithNewAssociation(numStep,src_vars[0], FC_AT_ELEMENT, NULL, &dst_vars[0] );
  fail_unless(rc == FC_SUCCESS, "failed null output variable test");
  rc = fc_getVariableName(dst_vars[0][0], &tst_name);
  fail_unless(rc==FC_SUCCESS, "couldn't get new variable name");
  fail_unless(!strcmp(tst_name, "clover"), "new variable name didn't match expected");
  free(tst_name);
  fc_deleteSeqVariable(numStep, dst_vars[0]);
  free(dst_vars[0]);

  //Cleanup seqvar
  fc_deleteSeqVariable(numStep, src_vars[0]);
  free(src_vars[0]);


  //Global Variable should be weeded out
  rc = fc_createGlobalSeqVariable(dataset, seq, "globboramma", &test_numStep, &test_gsv);
  fail_unless(rc==FC_SUCCESS,"Couldn't create a global variable");
  for(i=0; i<test_numStep; i++){
    rc = fc_setVariableData(test_gsv[i], 1, numComponent, FC_AT_WHOLE_DATASET,
			    FC_MT_VECTOR, FC_DT_DOUBLE,
			    (void *)vvals);\
    fail_unless(rc==FC_SUCCESS,"Couldn't set a global variable");
  }
  
  rc = fc_copySeqVariableWithNewAssociation(test_numStep,test_gsv, FC_AT_ELEMENT, "Booya", &dst_vars[1]);
  fail_unless(rc != FC_SUCCESS, "Global var copy worked when it should have failed");

  fc_deleteSeqVariable(test_numStep, test_gsv);
  free(test_gsv);

  //Cleanup
  free(res);
  free(dst_vars);
  free(src_vars);
  free(vvals);
  free(meshes);






  //Ok finally, lets test the calls with a simple dataset that's been
  //processed by hand offline.
  {
    double *co;
    FC_Mesh m;
    FC_Variable *var_results, *var_verts;
    double *tst_data;

    //We're building a simple 3x3x2 structure here with 4 cube elements
    //
    //  Vertex are 00-17  and elements A-D
    //
    // z=0:          z=1:
    //     00--01--02       09--10--11         
    //      | A | B |        | A | B |
    //     03--04--05       12--13--14 
    //      | C | D |        | C | D |
    //     06--07--08       15--16--17
    


    int newconns[] = { 0,1,4,3,9,10,13,12,  
		       1,2,5,4,10,11,14,13,
		       3,4,7,6,12,13,16,15,
		       4,5,8,7,13,14,17,16 };

    //First test is vertex->element
    double exp_elem[] = { (double)(0+1+4+3+9+10+13+12)/8.0,
			  (double)(1+2+5+4+10+11+14+13)/8.0,
			  (double)(3+4+7+6+12+13+16+15)/8.0,
			  (double)(4+5+8+7+13+14+17+16)/8.0 };
    int exp_elem_num = 4;

    //Second test is vertex->face.
    // Note: we build the face list by extracting all faces from each element
    //       and then removing redundant faces.
    //
    // This list is built cube-by-cube, and then redundant faces are commented out
    double exp_face[] = { //First cube
                          (double)(0+1+9+10)/4.0,
			  (double)(1+4+10+13)/4.0,
			  (double)(3+4+12+13)/4.0,
			  (double)(0+3+9+12)/4.0,
			  (double)(0+1+4+3)/4.0,
			  (double)(9+10+13+12)/4.0,
			  //Second cube
			  (double)(1+2+10+11)/4.0,
			  (double)(2+5+11+14)/4.0,
			  (double)(4+5+13+14)/4.0,
			  //(double)(4+1+13+10)/4.0, //Remove: same as second in first cube
			  (double)(1+2+4+5)/4.0,
			  (double)(10+11+13+14)/4.0,
			  //Third cube
			  //(double)(3+4+12+13)/4.0, //Remove: same as third in first cube
			  (double)(4+7+13+16)/4.0,
			  (double)(6+7+15+16)/4.0,
			  (double)(3+6+12+15)/4.0,
			  (double)(3+4+7+6)/4.0,
			  (double)(12+13+16+15)/4.0,
			  //Fourth cube
			  //(double)(4+5+13+14)/4.0, //Remove: same as third in second cube
			  (double)(5+8+14+17)/4.0,
			  (double)(8+7+16+17)/4.0,
			  //(double)(7+4+16+13)/4.0, //Remove: same as first in third cube
			  (double)(4+5+8+7)/4.0,
			  (double)(13+14+17+16)/4.0};


    int exp_face_num = 4*6 - 4;


    double exp_edge[] =  { //First cube
                           //Face 0 (0,1,10,9): (do all)
                           (double) (0+1)/2.0,
			   (double) (1+10)/2.0,
			   (double) (9+10)/2.0,
			   (double) (0+9)/2.0,
			   //Face 1 (1,4,13,10): Already did 10-1
			   (double) (1+4)/2.0,
			   (double) (4+13)/2.0,
			   (double) (13+10)/2.0,	   
			   //Face 2 (4,3,12,13): Already did 13-4
			   (double) (4+3)/2.0,
			   (double) (3+12)/2.0,
			   (double) (12+13)/2.0,
			   //Face 3 (3,0,9,12): already did 0-9,12-3 
			   (double) (3+0)/2.0,			  
			   (double) (9+12)/2.0,
			   //Face 4,5: Top/bottom faces already done

			   //Second cube
			   //Face 0 (1,2,11,10): already did 1-10
			   (double) (1+2)/2.0,
			   (double) (2+11)/2.0,
			   (double) (11+10)/2.0,
			   //Face 1 (2,5,14,11): already did 11-2
			   (double) (2+5)/2.0,
			   (double) (5+14)/2.0,
			   (double) (14+11)/2.0,
			   //Face 2 (5,4,13,14): already did 4-13,14-5
			   (double) (5+4)/2.0,
			   (double) (13+14)/2.0,
			   //Face 3 (4,1,10,13) same as cube 0, face 1
			   //Face 4,5: Top/bottom faces already done

			   //Third cube
			   //Face 0 (3,4,13,12): same as cube 0 face 2
			   //Face 1 (4,7,16,13): already did 4-7
			   (double) (4+7)/2.0,
			   (double) (7+16)/2.0,
			   (double) (16+13)/2.0,
			   //Face 2 (7,6,15,16): already did 16-7
			   (double) (7+6)/2.0,
			   (double) (6+15)/2.0,
			   (double) (15+16)/2.0,
			   //Face 3 (6,3,12,15): already did 3-12,15-6
			   (double) (6+3)/2.0,
			   (double) (12+15)/2.0,
			   //Face 4,5: top/bottom faces already done

			   //Fourth cube
			   //Face 0 (4,5,14,13): same as cube 1 face 2
			   //Face 1 (5,8,17,14): already did 5-14
			   (double) (5+8)/2.0,
			   (double) (8+17)/2.0,
			   (double) (17+14)/2.0,
			   //Face 2 (8,7,16,17): already did 8-17, 7-16
			   (double) (8+7)/2.0,
			   (double) (16+17)/2.0
			   //Face 3 (7,4,13,16) same as cube 2 face 1
			   //Face 4,5: Top/bottom faces already done
                         };

    int exp_edge_num = 12*2+9; /* 12 edges per z-plane, 9 edges to connect planes */

    typedef struct {
      FC_AssociationType assocType;
      int                num;
      double            *vals;
      char              *test_name;
    } exp_t;


    //Stuff all the tests into a simple struct so it's easier to manage
    exp_t expected[] = { {FC_AT_ELEMENT, exp_elem_num, exp_elem, "vertex2element"},
                         {FC_AT_FACE,    exp_face_num, exp_face, "vertex2face"},
			 {FC_AT_EDGE,    exp_edge_num, exp_edge, "vertex2edge"} };
    int expected_num = 3;



    //Coords don't really matter
    co = (double *) malloc(3*18*sizeof(double));
    fail_unless(co!=NULL, "no malloc");
    for(i=0;i<18;i++)
      for(j=0;j<3;j++)
	co[i*3 + j] = (double) i;

    rc = fc_createMesh(dataset, "hexer mesher", &m);
    fail_unless(rc==FC_SUCCESS, "failed to create mesh");
    rc = fc_setMeshCoords(m, 3, 18, co);
    fail_unless(rc==FC_SUCCESS, "failed to set mesh");
    rc = fc_setMeshElementConns(m, FC_ET_HEX, 4, newconns);
    fail_unless(rc == FC_SUCCESS, "failed to set element conns");

    //Get an existing sequence
    rc = fc_getSequenceByName(dataset,"time", &numReturnSequences, &returnSequences);
    fail_unless(rc == FC_SUCCESS, "failed to get sequences");
    seq = returnSequences[0];
    free(returnSequences);

    //Create seqVar and fill in the data
    rc = fc_createSeqVariable(m, seq, "bosequeous", &numStep, &var_verts);
    fail_unless(rc == FC_SUCCESS, "failed to create seqVar");
    for(s=0; s<numStep; s++){
      //Scale data for each step
      for(i=0; i<18; i++)
	for(j=0; j<3;j++)
	  co[i*3 + j] = (double) (1000*s + i);

      rc = fc_setVariableData(var_verts[s], 18, 3, FC_AT_VERTEX,
			      FC_MT_VECTOR, FC_DT_DOUBLE,
			      (void *)co);
      fail_unless(rc == FC_SUCCESS, "abort: failed to set variable");
    }




    //fc_printMesh(m, "mesh", 0,0,1,1);
    //fc_printVariable(var_verts, "vert data",1);


    //Try all the tests
    for(i=0;i<expected_num; i++){
      //printf("Trying test %s\n", expected[i].test_name);

      rc = fc_copySeqVariableWithNewAssociation(numStep, var_verts, expected[i].assocType, expected[i].test_name, &var_results);
      fail_unless(rc == FC_SUCCESS, "Couldn't copy seqVariable to a particular association type");

      //fc_printVariable(var_elems[0], expected[i].test_name, 1);
      
      rc = fc_getVariableNumDataPoint(var_results[0], &tst_num);
      fail_unless(rc == FC_SUCCESS, "Failed to get data points");
      fail_unless(tst_num == expected[i].num, "Didn't get the expected number of values");
      
      for(s=0; s<numStep; s++){
	rc = fc_getVariableDataPtr(var_results[s], (void **) &tst_data);
	fail_unless(rc == FC_SUCCESS, "Failed to get data pointer");
      
	for(j=0; j<expected[i].num; j++) {
	  //printf("%s edge %d:: %.3lf %.3lf\n", expected[i].test_name, j, tst_data[j*numComponent + 0], expected[i].vals[j]);
	  for(k=0; k<numComponent;k++){  
	    fail_unless( tst_data[j*numComponent + k]==(s*1000.0 + expected[i].vals[j]), "Converted data did not match expected value");
	  }
	}
      }
      
      //Done with particular test. Cleanup this round
      fc_deleteSeqVariable(numStep, var_results);
      free(var_results);
    }

    //Free up stuff
    fc_deleteSeqVariable(numStep, var_verts);
    free(var_verts);
    fc_deleteMesh(m);
    fc_deleteSequence(seq);
    free(co);

  } //End of canned test





}
END_TEST





// *********************************************
// ***** Populate the Suite with the tests
// *********************************************

Suite *variable_suite(void)
{
  Suite *suite = suite_create("Variable");

  TCase *tc_private_var = tcase_create(" - Private Variable ");
  TCase *tc_var = tcase_create(" - Variable Interface ");
  TCase *tc_seq_var = tcase_create(" - SeqVar Interface ");
  TCase *tc_glb_var = tcase_create(" - Global Var Interface ");
  TCase *tc_glb_seq_var = tcase_create(" - Global SeqVar Intf. ");
  TCase *tc_convert = tcase_create(" - Convert vars ");

  // private variable
  suite_add_tcase(suite, tc_private_var);
  tcase_add_checked_fixture(tc_private_var, variable_setup, 
			    variable_teardown);

  tcase_add_test(tc_private_var, slot_new_delete);
  tcase_add_test(tc_private_var, slot_new_delete_seq);


  // variable interface  
  suite_add_tcase(suite, tc_var);
  tcase_add_checked_fixture(tc_var, variable_setup, variable_teardown);
  tcase_add_test(tc_var, create_get_delete);
  tcase_add_test(tc_var, metadata_query);
  tcase_add_test(tc_var, data_query);
  tcase_add_test(tc_var, get_data_as);
  tcase_add_test(tc_var, copy_test);
  tcase_add_test(tc_var, release_var);
  tcase_add_test(tc_var, print);

  // seq variable interface  
  suite_add_tcase(suite, tc_seq_var);
  tcase_add_checked_fixture(tc_seq_var, variable_setup, variable_teardown);
  tcase_add_test(tc_seq_var, seq_create_get_delete);
  tcase_add_test(tc_seq_var, seq_metadata_query);
  tcase_add_test(tc_seq_var, seq_data_query);
  tcase_add_test(tc_seq_var, seq_copy_test);
  tcase_add_test(tc_seq_var, seq_release_var);

  // global variable interface
  suite_add_tcase(suite, tc_glb_var);
  tcase_add_checked_fixture(tc_glb_var, variable_setup, variable_teardown);
  tcase_add_test(tc_glb_var, glb_create_get_delete);
  tcase_add_test(tc_glb_var, glb_metadata_query);
  tcase_add_test(tc_glb_var, glb_data_query);
  tcase_add_test(tc_glb_var, glb_copy_test);

  // gloval seq variable interface
  suite_add_tcase(suite, tc_glb_seq_var);
  tcase_add_checked_fixture(tc_glb_seq_var, variable_setup, variable_teardown);
  tcase_add_test(tc_glb_seq_var, glb_seq_create_get_del);
  tcase_add_test(tc_glb_seq_var, glb_seq_metadata_query);
  tcase_add_test(tc_glb_seq_var, glb_seq_data_query);
  tcase_add_test(tc_glb_seq_var, glb_seq_copy_test);

  // test creating new vars from other vars
  suite_add_tcase(suite, tc_convert);
  tcase_add_checked_fixture(tc_convert, variable_setup, variable_teardown);
  tcase_add_test(tc_convert, vars_to_seq_var);
  tcase_add_test(tc_convert, glb_vars_to_seq_var);
  tcase_add_test(tc_convert, split_merge_comps);
  tcase_add_test(tc_convert, comps_by_name);
  tcase_add_test(tc_convert, glb_split_merge_comps);
  tcase_add_test(tc_convert, seq_split_merge_comps);
  tcase_add_test(tc_convert, glb_seq_split_merge_comps);
  tcase_add_test(tc_convert, vars_copy_with_association);
  tcase_add_test(tc_convert, seqvars_copy_with_association);


  return suite;
}
