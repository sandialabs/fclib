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
 * \file checkdataset.c
 * \brief Unit tests for \ref Dataset Module.
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/unittest/checkdataset.c,v $
 * $Revision: 1.31 $ 
 * $Date: 2006/10/19 03:14:52 $
 *
 * \modifications
 *    9/8/04 WSK, split off of checklibrary
 */

#include <stdlib.h>
#include <string.h>
#include <check.h>
#include "fc.h"
#include "fcP.h"
#include "checkall.h"

// **** Fixtures
static void dataset_setup(void) {
  FC_ReturnCode rc;

  if (isForking) {
    fc_setLibraryVerbosity(fc_messages);
    rc = fc_initLibrary();
    fail_unless(rc == FC_SUCCESS, "test aborted: failed to init library");
  }
}

static void dataset_teardown(void) {
  FC_ReturnCode rc;

  if (isForking) {
    rc = fc_finalLibrary();
    fail_unless(rc == FC_SUCCESS, "failed to final library at end of test");
  }
}

// Test that creating and deleting slots modifies table as expected
// Deleting should open up slots that get reused.
START_TEST(slot_new_delete)
{
  FC_ReturnCode rc;
  int i;
  int num = 7;
  // uID_incs = diff between current uID and startID
  int uID_start, uID_incs[7] = { 0, 7, 2, 8, 4, 9, 6 };
  FC_Dataset datasets[7];
  
  // Create num
  for (i = 0; i < num; i++) {
    rc = fc_createDataset("fake", &datasets[i]);
    fail_unless(rc == FC_SUCCESS, "failed to create dataset");
  }
  uID_start = datasets[0].uID;

  // Delete "odd" slots
  for (i = 1; i < num; i+=2) {
    rc = fc_deleteDataset(datasets[i]);
    fail_unless(rc == FC_SUCCESS, "failed to delete dataset");
  }

  // Create some more (fills in the cracks)
  for (i = 1; i < num; i+=2) {
    rc = fc_createDataset("fake", &datasets[i]);
    fail_unless(rc == FC_SUCCESS, "failed to create dataset");
  }

  // Check slotID and uID
  for (i = 0; i < num; i++) {
    fail_unless(datasets[i].slotID == i, "mismatch of slot id");
    fail_unless(datasets[i].uID == uID_start + uID_incs[i],
		"mismatch of uID");
  }

  // cleanup
  for (i = 0; i < num; i++) {
    rc = fc_deleteDataset(datasets[i]);
    fail_unless(rc == FC_SUCCESS, "failed to delete dataset at end of test");
  }
}
END_TEST

// test creating and deleting datasets, also test querying library
// for datasets before, inbetween and after.
START_TEST(create_get_delete)
{
  FC_ReturnCode rc;
  int i;
  int numDataset = 10, temp_numDataset;
  FC_Dataset datasets[10], *temp_datasets, temp_dataset;
  FC_Dataset badDataset = { 999, 999 };
  char names[10][100] = { "one", "../two", "../three", "four", "five", "six",
			  "seven", "eight", "/a/b/nine", "a b ten" };
  char newName[100] = { "sparkly new" };
  
  if (isForking) {
    // create without initializing library is not ok
    rc = fc_createDataset(names[0], &temp_dataset);
    fail_unless(rc == FC_ERROR, 
		"fc_createDataset should fail if library was not inited");
    fail_unless(FC_HANDLE_EQUIV(temp_dataset, FC_NULL_DATASET),
		"fail should return FC_NULL_DATASET");
    
    // get without initializing library is not ok
    temp_numDataset = 999; temp_datasets = (FC_Dataset*)1;
    rc = fc_getDatasets(&temp_numDataset, &temp_datasets);
    fail_unless(rc == FC_ERROR,
		"fc_getDatasets should fail if library was not inited");
    fail_unless(temp_numDataset == -1 && temp_datasets == NULL,
		"fail should return -1 & NULL");
    temp_dataset = badDataset;
    rc = fc_getNumDataset(&temp_numDataset);
    fail_unless(rc == FC_ERROR,
		"fc_getNumDataset should fail if library was not inited");
    fail_unless(temp_numDataset == -1 && temp_datasets == NULL,
		"fail should return -1 & NULL");
    rc = fc_getDatasetByName(names[0], &temp_numDataset,&temp_datasets);
    fail_unless(rc == FC_ERROR,
		"fc_getDatasetByName should fail if library was not inited");
    fail_unless((temp_numDataset == -1 && temp_datasets == NULL),
		"fail should return numDatasets = -1 && NULL dataset handle");
    
    // init library
    fc_setLibraryVerbosity(fc_messages);
    fc_initLibrary();
  }
  else 
    printf("**CK_NOFORK: Can't test create/get on un-inited library\n");
 
  // get datasets (should be none)
  temp_numDataset = 999; temp_datasets = (FC_Dataset*)1;
  rc = fc_getDatasets(&temp_numDataset, &temp_datasets);
  fail_unless(rc == FC_SUCCESS, "failed to get datasets from empty library");
  fail_unless(temp_numDataset == 0 && temp_datasets == NULL,
	      "should return 0 if no datasets created yet");
  rc = fc_getNumDataset(&temp_numDataset);
  fail_unless(rc == FC_SUCCESS, "failed to get numDataset from empty library");
  fail_unless(temp_numDataset == 0,
	      "should return 0 if no datasets created yet");
  
  // create some datasets
  for (i = 0; i < numDataset; i++) {
    rc = fc_createDataset(names[i], &datasets[i]);
    fail_unless(rc == FC_SUCCESS, "failed to create dataset");
    fail_unless(!FC_HANDLE_EQUIV(datasets[i], FC_NULL_DATASET),
		"created dataset should not = FC_NULL_DATASET");
  }

  // get datasets (should be numDataset);
  rc = fc_getDatasets(&temp_numDataset, &temp_datasets);
  fail_unless(rc == FC_SUCCESS, "failed to get datasets");
  fail_unless(temp_numDataset == numDataset, "mismatch of numDataset");
  for (i = 0; i < numDataset; i++) 
    fail_unless(FC_HANDLE_EQUIV(temp_datasets[i], datasets[i]), 
		"mismatch of dataset handles");
  free(temp_datasets);
  rc = fc_getNumDataset(&temp_numDataset);
  fail_unless(rc == FC_SUCCESS, "failed to get numDataset");
  fail_unless(temp_numDataset == numDataset, "mismatch of numDataset");
  for (i = 0; i < numDataset; i++) {
    rc = fc_getDatasetByName(names[i],&temp_numDataset,&temp_datasets);
    fail_unless(rc == FC_SUCCESS, "failed to get dataset by name");
    fail_unless(temp_numDataset == 1, "failed to find matching dataset");
    fail_unless(FC_HANDLE_EQUIV(temp_datasets[0], datasets[i]),
		"mismatch of dataset handles");
    free(temp_datasets);
  }

  //create temporary dataset to check for multiple return
  rc = fc_createDataset(names[0], &temp_dataset);
  fail_unless(rc == FC_SUCCESS, "failed to create dataset");
  fail_unless(!FC_HANDLE_EQUIV(temp_dataset, FC_NULL_DATASET),
	      "created dataset should not = FC_NULL_DATASET");
  rc = fc_getDatasetByName(names[0],&temp_numDataset,&temp_datasets);
  fail_unless(rc == FC_SUCCESS, "failed to get dataset by name");
  fail_unless(temp_numDataset == 2, "failed to find matching dataset");
  fail_unless(((FC_HANDLE_EQUIV(temp_datasets[0], datasets[0]) &&
		FC_HANDLE_EQUIV(temp_datasets[1], temp_dataset)) ||
	       (FC_HANDLE_EQUIV(temp_datasets[1], datasets[0]) &&
		FC_HANDLE_EQUIV(temp_datasets[0], temp_dataset))),
	      "mismatch of dataset handles");
  fc_deleteDataset(temp_dataset);
  free(temp_datasets);
  

  // change the name of the first dataset
  rc = fc_changeDatasetName(datasets[0], newName);
  fail_unless(rc == FC_SUCCESS, "failed to change dataset name");
  rc = fc_getDatasetByName(names[0], &temp_numDataset,&temp_datasets);
  fail_unless(rc == FC_SUCCESS, "old name should work");
  fail_unless(temp_numDataset == 0, "old name shouldnt find anything");
  fail_unless(temp_datasets == NULL, "should be NULL array");
  rc = fc_getDatasetByName(newName, &temp_numDataset,&temp_datasets);
  fail_unless(rc == FC_SUCCESS, "new name should work");
  fail_unless(temp_numDataset == 1, "new name should work for find");
  fail_unless(FC_HANDLE_EQUIV(temp_datasets[0], datasets[0]),
	      "mismatch of dataset handle");
  free(temp_datasets);

  // delete half of the datasets (alternate)
  for (i = 0; i < numDataset; i+=2) {
    rc = fc_deleteDataset(datasets[i]);
    fail_unless(rc == FC_SUCCESS, "failed to delete dataset");
  }

  // get datasets (should be numDataset/2);
  fc_getDatasets(&temp_numDataset, &temp_datasets);
  fail_unless(rc == FC_SUCCESS, "failed to get datasets");
  fail_unless(temp_numDataset == numDataset/2, "mismatch of numDataset");
  for (i = 0; i < numDataset/2; i++) 
    fail_unless(FC_HANDLE_EQUIV(temp_datasets[i], datasets[i*2+1]), 
		"mismatch of dataset handles");
  free(temp_datasets);
  for (i = 0; i < numDataset/2; i++) {
    rc = fc_getDatasetByName(names[i*2+1], &temp_numDataset,&temp_datasets);
    fail_unless(rc == FC_SUCCESS, "failed to get dataset by name");
    fail_unless(temp_numDataset == 1, "should find matching dataset");
    fail_unless(FC_HANDLE_EQUIV(temp_datasets[0], datasets[i*2+1]),
		"mismatch of dataset handles");
    free(temp_datasets);
  }

  // delete remaining datasets
  for (i = 1; i < numDataset; i+=2) { 
    rc = fc_deleteDataset(datasets[i]);
    fail_unless(rc == FC_SUCCESS, "failed to delete dataset");
  }
  
  // get datasets (should be none)
  rc = fc_getDatasets(&temp_numDataset, &temp_datasets);
  fail_unless(rc == FC_SUCCESS, "failed to get datasets from empty library");
  fail_unless(temp_numDataset == 0 && temp_datasets == NULL,
	      "should return 0 if all datasets deleted");

  // create one more dataset for further testing
  rc = fc_createDataset(names[0], &datasets[0]);
  fail_unless(rc == FC_SUCCESS, "aborted: failed to create dataset dataset");

  // ---- test special cases

  // not an error to delete FC_NULL_DATASET
  rc = fc_deleteDataset(FC_NULL_DATASET);
  fail_unless(rc == FC_SUCCESS, "should not error to delete NULL dataset");

  // fc_getDatasetByName() - returns none if it doesnt exist
  rc = fc_getDatasetByName("this name doesn't exist", &temp_numDataset,
			   &temp_datasets);
  fail_unless(rc == FC_SUCCESS, 
	      "should work for a non existent dataset");
  fail_unless(temp_numDataset == 0 , "should find no matching datasets");
  fail_unless(temp_datasets == NULL, "should return NULL array");
   // ---- test error conditions

  // bad args to fc_createDataset()
  temp_dataset = badDataset;
  rc = fc_createDataset(NULL, &temp_dataset);
  fail_unless(rc != FC_SUCCESS, 
	      "should have failed to create dataset with null name");
  fail_unless(FC_HANDLE_EQUIV(temp_dataset, FC_NULL_DATASET),
	      "fail should return NULL dataset");
  rc = fc_createDataset(names[0], NULL);
  fail_unless(rc != FC_SUCCESS, 
	      "should have failed to create dataset with null handle");

  // bad args to fc_changeDatasetName()
  rc = fc_changeDatasetName(datasets[0], NULL);
  fail_unless(rc != FC_SUCCESS, "should fail if new name is NULL");
  rc = fc_getDatasetByName(names[0], &temp_numDataset,&temp_datasets);
  fail_unless(rc == FC_SUCCESS && temp_numDataset == 1 &&
	      FC_HANDLE_EQUIV(temp_datasets[0], datasets[0]),
	      "should not change name of dataset");
  free(temp_datasets);

  // bad args to fc_getDatasets()
  temp_datasets = (FC_Dataset*)1;
  rc = fc_getDatasets(NULL, &temp_datasets);
  fail_unless(rc != FC_SUCCESS, 
	      "should have failed to get dataset with NULL numDataset");
  fail_unless(temp_datasets == NULL, "mismatch of datasets");
  
  // bad args to fc_getNumDataset()
  rc = fc_getNumDataset(NULL);
  fail_unless(rc != FC_SUCCESS, 
	      "should have failed to get numDataset with NULL numDataset");
  temp_numDataset = 999;
  rc = fc_getDatasets(&temp_numDataset, NULL);
  fail_unless(rc != FC_SUCCESS, "should error id datasets = NULL");
  fail_unless(temp_numDataset == -1, "fail should return NULL");
  
  // bad args to fc_getDatasetByName()
  temp_dataset = badDataset;
  rc = fc_getDatasetByName(NULL, &temp_numDataset,&temp_datasets);
  fail_unless(rc != FC_SUCCESS, 
	      "should have failed to get dataset with null name");
  fail_unless(temp_numDataset == -1, "fail should return -1");
  fail_unless(temp_datasets == NULL, "fail should return NULL array");

  rc = fc_getDatasetByName(names[0], NULL,&temp_datasets);
  fail_unless(rc != FC_SUCCESS, 
	      "should have failed to get dataset with null number arg");
  fail_unless(temp_datasets == NULL, "fail should return NULL array");

  rc = fc_getDatasetByName(names[0], &temp_numDataset,NULL);
  fail_unless(rc != FC_SUCCESS, 
	      "should have failed to get dataset with null return arg");
  fail_unless(temp_numDataset == -1, "fail should return -1");

  // bad args to fc_deleteDataset()
  rc = fc_deleteDataset(badDataset);
  fail_unless(rc != FC_SUCCESS, 
	      "should error to delete nonexistent dataset");

  // --- done

  // cleanup
  rc = fc_deleteDataset(datasets[0]);
  fail_unless(rc == FC_SUCCESS, "final delete dataset failed");

  if (isForking) {
    fc_finalLibrary();

    // create on finaled library is not ok
    rc = fc_createDataset(names[0], &temp_dataset);
    fail_unless(rc == FC_ERROR, 
		"fc_createDataset should fail if library was finaled");
    fail_unless(FC_HANDLE_EQUIV(temp_dataset, FC_NULL_DATASET),
		"fail should return FC_NULL_DATASET");
    
    // get without initializing library is not ok
    temp_numDataset = 999; temp_datasets = (FC_Dataset*)1;
    rc = fc_getDatasets(&temp_numDataset, &temp_datasets);
    fail_unless(rc == FC_ERROR,
		"fc_getDatasets should fail if library was finaled");
    fail_unless(temp_numDataset == -1 && temp_datasets == NULL,
		"fail should return -1 & NULL");
    temp_dataset = badDataset;
    rc = fc_getDatasetByName(names[0], &temp_numDataset,
			     &temp_datasets);
    fail_unless(rc == FC_ERROR && temp_numDataset == -1 && temp_datasets == NULL,
		"fc_getDatasets should fail if library was finaled");
  }
  else 
    printf("**CK_NOFORK: Can't test create/get on finaled library\n");
}
END_TEST

// query meta data
START_TEST(metadata_query)
{
  FC_ReturnCode rc;
  char name[20] = "blue berry", *temp_name;
  FC_Dataset dataset, badDataset = { 999, 999 };

  if (isForking) {
    fc_setLibraryVerbosity(fc_messages);
    rc = fc_initLibrary();
    fail_unless(rc == FC_SUCCESS, "test aborted: failed to init library");
  }

  // create a dataset to play with
  rc = fc_createDataset(name, &dataset);
  fail_unless(rc == FC_SUCCESS, "failed to created dataset");
 
  // test fc_isDatasetValid()
  fail_unless(fc_isDatasetValid(dataset), 
	      "failed to validate a valid dataset");

  // test fc_getDatasetName()
  temp_name = NULL;
  rc = fc_getDatasetName(dataset, &temp_name);
  fail_unless(rc == FC_SUCCESS, "failed to get dataset name");
  fail_unless(!strcmp(name, temp_name), "mismatch of name");
  // should have returned a copy, mess with temp_name & try again to make sure
  temp_name[0] = 'q';
  free(temp_name);
  rc = fc_getDatasetName(dataset, &temp_name);
  fail_unless(rc == FC_SUCCESS, "failed to get dataset name");
  fail_unless(!strcmp(name, temp_name), "mismatch of name");
  free(temp_name);

  // --- check with bad args

  // fc_isDatasetValid()
  fail_unless(!fc_isDatasetValid(badDataset), 
	      "badDataset should NOT be valid");
  fail_unless(!fc_isDatasetValid(FC_NULL_DATASET), 
	      "FC_NULL_DATASET should not be valid");

  // fc_getDatasetName()
  temp_name = (char*)1;
  rc = fc_getDatasetName(badDataset, &temp_name);
  fail_unless(rc != FC_SUCCESS, "badDataset should NOT return a name");
  fail_unless(temp_name == NULL, "fail should return NULL name");
  rc = fc_getDatasetName(dataset, NULL);
  fail_unless(rc != FC_SUCCESS, "should fail for NULL name");

  // --- all done

  // cleanup & one last test
  fc_deleteDataset(dataset);
  fail_unless(!fc_isDatasetValid(dataset),
	      "handle should not be valid after a delete");

  if (isForking) {
    rc = fc_finalLibrary();
    fail_unless(rc == FC_SUCCESS, "abort: failed to final library");
  }
}
END_TEST

// test error input for fc_printDataset()
START_TEST(print)
{
  FC_ReturnCode rc;
  FC_Dataset badDataset = { 999, 999 };

  // fc_printDataset()
  rc = fc_printDataset(badDataset, "Bad Dataset!");
  fail_unless(rc != FC_SUCCESS, "should fail if bad dataset");
}
END_TEST 

// ************ The following functions should probably be run after
// **********   testing the interfaces for seq, mesh, subset, & var

// test that children get copied
START_TEST(copy_test)
{
  FC_ReturnCode rc;
  int i, j, k;
  char name[20] = "blue berry", copy_name[20] = "banana nut";
  char* temp_name;
  FC_Dataset dataset, copy_dataset, badDataset = { 999, 999 };
  int numSequence = 3, temp_numSequence;
  FC_Sequence sequences[3], *temp_sequences;
  char* seq_names[3] = { "seq1", "seq2", "seq3" };
  int numSteps[3] = { 2, 3, 4 }, temp_numStep;
  FC_DataType dataType = FC_DT_DOUBLE, temp_dataType;
  double seqCoords[3][4] = { { 0, 1 }, { 2, 3, 4 }, { 5, 6, 7, 8 } };
  void* temp_data;
  int numMesh = 3, temp_numMesh;
  FC_Mesh meshes[3], *temp_meshes;
  char* mesh_names[3] = { "mesh1", "mesh2", "mesh3" };
  int numDim = 3, temp_numDim;
  int topodims[3] = { 0, 2, 3 }, temp_topodim;
  FC_ElementType elemTypes[3] = { FC_ET_POINT, FC_ET_QUAD, FC_ET_HEX };
  FC_ElementType temp_elemType;
  int numVertPerElems[3] = { 1, 4, 8 };
  int numVerts[3] = { 3, 6, 8 }, temp_numVert;
  int numElems[3] = { 3, 2, 1 }, temp_numElem;
  int conns[3][3*8] = { {  0, 1, 2 }, 
			{  0, 1, 3, 2, 1, 2, 5, 4 },
			{  0, 1, 3, 2, 4, 5, 7, 6 } };
  double coords[3][3*8] = { { 0, 0, 0, 1, 1, 1, 2, 2, 2 },
			    { 0, 0, 0, 1, 0, 0, 2, 0, 0, 0, 1, 0, 1, 1, 0, 2, 1, 0 },
			    { 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0,
			      0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1 } };
  int* temp_conns;
  double* temp_coords;
  int numSubset = 3, temp_numSubset;
  FC_Subset subset, *temp_subsets;
  char* subset_names[3] = { "sub1", "sub2", "sub3" };
  FC_AssociationType assocs[3] = { FC_AT_VERTEX, FC_AT_EDGE, FC_AT_ELEMENT };
  FC_AssociationType temp_assoc;
  int numMembers[3] = { 2, 3, 1 }, temp_numMember;
  int maxNumMembers[3] = { 8, 12, 1 }, temp_maxNumMember;
  int memberIDs[3][3] = { { 0, 2 }, { 1, 3, 7 }, { 0 } }, *temp_memberIDs;
  int numVar = 3, temp_numVar;
  FC_Variable var, *temp_vars;
  char* var_names[3] = { "var1", "var2", "var3" };
  // resuse assocs from subsets
  int numDataPoints[3] = { 8, 12, 1 }, temp_numDataPoint; // same as maxNumMember
  int numComps[3] = { 1, 3, 4 }, temp_numComp;
  FC_MathType mathTypes[3] = { FC_MT_SCALAR, FC_MT_VECTOR, FC_MT_SYMTENSOR };
  FC_MathType temp_mathType;
  double datas[3][12] = { { .1, .2, .3, .4, .5, .6, .7, .8 }, 
			  { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 }, 
			  { 999 } };
  int numSeqVar = 3, temp_numSeqVar, *temp_numStepPerSeqVar;
  FC_Variable *seqVar, **temp_seqVars;
  char* seqVar_names[3] = { "seqVar1", "seqVar2", "seqVar3" };
  int numGlbVar = 3, temp_numGlbVar;
  FC_Variable glbVar, *temp_glbVars;
  char* glbVar_names[3] = { "glbVar1", "glbVar2", "glbVar3" };
  int numGlbSeqVar = 3, temp_numGlbSeqVar;
  FC_Variable *glbSeqVar, **temp_glbSeqVars;
  char* glbSeqVar_names[3] = { "glbSeqVar1", "glbSeqVar2", "glbSeqVar3" };
  
  // setup
  rc = fc_createDataset(name, &dataset);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create dataset");

  // special case - copy with NULL name will use original's name
  // (testing with empty dataset so will be fast)
  rc = fc_copyDataset(dataset, NULL, &copy_dataset);
  fail_unless(rc == FC_SUCCESS, "should NOT fail if no name");
  rc = fc_getDatasetName(copy_dataset, &temp_name);
  fail_unless(rc == FC_SUCCESS, "failed to get name of copy");
  fail_unless(!strcmp(temp_name, name), "mismatch of name");
  free(temp_name);
  rc = fc_deleteDataset(copy_dataset);
  fail_unless(rc == FC_SUCCESS, "failed to delete copied dataset");

  // fc_copyDataset() -- bad args
  rc = fc_copyDataset(badDataset, copy_name, &copy_dataset);
  fail_unless(rc != FC_SUCCESS, "should fail if bad dataset");
  fail_unless(FC_HANDLE_EQUIV(copy_dataset, FC_NULL_DATASET), 
	      "fail should return NULL");
  rc = fc_copyDataset(dataset, copy_name, NULL); 
  fail_unless(rc != FC_SUCCESS, "should fail if no dataset handle");

  // -- do it

  // put some of everything on dataset
  // sequences
  for (i = 0; i < numSequence; i++) {
    rc = fc_createSequence(dataset, seq_names[i], &sequences[i]);
    fail_unless(rc == FC_SUCCESS, "abort: failed to create sequence");
    rc = fc_setSequenceCoords(sequences[i], numSteps[i], dataType, seqCoords[i]);
    fail_unless(rc == FC_SUCCESS, "abort: failed to set sequence coords");
  }
  // global vars
  for (i = 0; i < numGlbVar; i++) {
    rc = fc_createGlobalVariable(dataset, glbVar_names[i], &glbVar);
    fail_unless(rc == FC_SUCCESS, "abort: failed to create global variable");
    rc = fc_setVariableData(glbVar, 1, numComps[i], FC_AT_WHOLE_DATASET,
			    mathTypes[i], dataType, datas[i]);
    fail_unless(rc == FC_SUCCESS, "abort: failed to set global var data");
  }
  // global seq vars - just on last sequence
  for (i = 0; i < numGlbSeqVar; i++) {
    rc = fc_createGlobalSeqVariable(dataset, sequences[numSequence-1],
				    glbSeqVar_names[i], &temp_numStep, 
				    &glbSeqVar);
    fail_unless(rc == FC_SUCCESS, "abort: failed to create global seqVar");
    for (j = 0; j < temp_numStep; j++) {
      rc = fc_setVariableData(glbSeqVar[j], 1, numComps[i], FC_AT_WHOLE_DATASET,
			      mathTypes[i], dataType, datas[i]);
      fail_unless(rc == FC_SUCCESS, "abort: failed to set global seqVar data");
    }
    free(glbSeqVar);
  }
  // meshes
  for (i = 0; i < numMesh; i++) {
    rc = fc_createMesh(dataset, mesh_names[i], &meshes[i]);
    fail_unless(rc == FC_SUCCESS, "abort: failed to create mesh");
    rc = fc_setMeshCoords(meshes[i], numDim, numVerts[i], coords[i]);
    fail_unless(rc == FC_SUCCESS, "abort: failed to set mesh coords");
    rc = fc_setMeshElementConns(meshes[i], elemTypes[i], numElems[i],
				conns[i]);
    fail_unless(rc == FC_SUCCESS, "abort: failed to set mesh conns");
  }
  // subsets - just on last mesh
  for (i = 0; i < numSubset; i++) {
    rc = fc_createSubset(meshes[numMesh-1], subset_names[i], assocs[i],
			 &subset);
    fail_unless(rc == FC_SUCCESS, "abort: failed to create subset");
    for (j = 0; j < numMembers[i]; j++) {
      rc = fc_addMemberToSubset(subset, memberIDs[i][j]);
      fail_unless(rc == FC_SUCCESS, "abort: failed to add members");
    }
  }
  // vars - just on last mesh
  for (i = 0; i < numVar; i++) {
    rc = fc_createVariable(meshes[numMesh-1], var_names[i], &var);
    fail_unless(rc == FC_SUCCESS, "abort: failed to create variable");
    rc = fc_setVariableData(var, numDataPoints[i], numComps[i], assocs[i],
			    mathTypes[i], dataType, datas[i]);
    fail_unless(rc == FC_SUCCESS, "abort: failed to set var data");
  }
  // seqVars - just on last mesh and last sequence
  for (i = 0; i < numSeqVar; i++) {
    rc = fc_createSeqVariable(meshes[numMesh-1], sequences[numSequence-1],
			      seqVar_names[i], &temp_numStep, &seqVar);
    fail_unless(rc == FC_SUCCESS, "abort: failed to create seqVariable");
    for (j = 0; j < temp_numStep; j++) {
      rc = fc_setVariableData(seqVar[j], numDataPoints[i], numComps[i], 
			      assocs[i], mathTypes[i], dataType, datas[i]);
      fail_unless(rc == FC_SUCCESS, "abort: failed to set seqVar data");
    }
    free(seqVar);
  }

  // copy (and delete old so we don't accidentally test it)
  rc = fc_copyDataset(dataset, copy_name, &copy_dataset);
  fail_unless(rc == FC_SUCCESS, "should NOT fail");
  rc = fc_deleteDataset(dataset);
  fail_unless(rc == FC_SUCCESS, "failed to delete copied dataset");

  // test name
  rc = fc_getDatasetName(copy_dataset, &temp_name);
  fail_unless(rc == FC_SUCCESS, "failed to get name of copy");
  fail_unless(!strcmp(temp_name, copy_name), "mismatch of name");
  free(temp_name);
  // test sequences
  rc = fc_getSequences(copy_dataset, &temp_numSequence, &temp_sequences);
  fail_unless(rc == FC_SUCCESS, "should not fail to get sequences");
  fail_unless(temp_numSequence == numSequence, "mismatch of numSequence");
  for (i = 0; i < numSequence; i++) {
    rc = fc_getSequenceName(temp_sequences[i], &temp_name);
    fail_unless(rc == FC_SUCCESS, "failed to get sequence name");
    fail_unless(!strcmp(temp_name, seq_names[i]), "seq name mismatch");
    free(temp_name);
    rc = fc_getSequenceInfo(temp_sequences[i], &temp_numStep, &temp_dataType);
    fail_unless(rc == FC_SUCCESS, "failed to get sequence info");
    fail_unless(temp_numStep == numSteps[i], "mismatch of numStep");
    fail_unless(temp_dataType == dataType, "mismatch of dataType");
    rc = fc_getSequenceCoordsPtr(temp_sequences[i], &temp_data);
    fail_unless(rc == FC_SUCCESS, "failed to get sequence coords");
    fail_unless(!memcmp(temp_data, seqCoords[i], numSteps[i]*sizeof(double)),
		"mismatch of data");
  }
  free(temp_sequences);
  // test global vars
  rc = fc_getGlobalVariables(copy_dataset, &temp_numGlbVar, &temp_glbVars);
  fail_unless(rc == FC_SUCCESS, "should not fail to get global vars");
  fail_unless(temp_numGlbVar == numGlbVar, "mismatch of numGlbVar");
  for (i = 0; i < numGlbVar; i++) {
    rc = fc_getVariableName(temp_glbVars[i], &temp_name);
    fail_unless(rc == FC_SUCCESS, "failed to get var name");
    fail_unless(!strcmp(temp_name, glbVar_names[i]), "var name mismatch");
    free(temp_name);
    rc = fc_getVariableInfo(temp_glbVars[i], &temp_numDataPoint,
			    &temp_numComp, &temp_assoc, &temp_mathType, 
			    &temp_dataType);
    fail_unless(rc == FC_SUCCESS, "failed to get var info");
    fail_unless(temp_numDataPoint == 1, "mismatch of numDataPoints");
    fail_unless(temp_numComp == numComps[i], "mismatch of numComp");
    fail_unless(temp_assoc == FC_AT_WHOLE_DATASET, "mismatch of assoc");
    fail_unless(temp_mathType == mathTypes[i], "mismatch of mathType");
    fail_unless(temp_dataType == dataType, "mismatch of dataType");
    rc = fc_getVariableDataPtr(temp_glbVars[i], &temp_data);
    fail_unless(rc == FC_SUCCESS, "failed to get var data");
    fail_unless(!memcmp(temp_data, datas[i], sizeof(double)),
		"mismatch of var data");
  }
  free(temp_glbVars);
  // test global seq vars - are only on last sequence
  rc = fc_getGlobalSeqVariables(copy_dataset, &temp_numGlbSeqVar,
				&temp_numStepPerSeqVar, &temp_glbSeqVars);
  fail_unless(rc == FC_SUCCESS, "should not fail to get glbSeqVars");
  fail_unless(temp_numGlbSeqVar == numGlbSeqVar, "mismatch of numGlbSeqVar");
  for (i = 0; i < numGlbSeqVar; i++) {
    fail_unless(temp_numStepPerSeqVar[i] == numSteps[numSequence-1],
		"mismatch of numStep & probably the sequence");
    rc = fc_getVariableName(temp_glbSeqVars[i][0], &temp_name);
    fail_unless(rc == FC_SUCCESS, "failed to get glbSeqVar name");
    fail_unless(!strcmp(temp_name, glbSeqVar_names[i]), 
		"glbSeqVar name mismatch");
    free(temp_name);
    rc = fc_getVariableInfo(temp_glbSeqVars[i][0], &temp_numDataPoint,
			    &temp_numComp, &temp_assoc, &temp_mathType, 
			    &temp_dataType);
    fail_unless(rc == FC_SUCCESS, "failed to get glbSeqVar info");
    fail_unless(temp_numDataPoint == 1, "mismatch of numDataPoints");
    fail_unless(temp_numComp == numComps[i], "mismatch of numComp");
    fail_unless(temp_assoc == FC_AT_WHOLE_DATASET, "mismatch of assoc");
    fail_unless(temp_mathType == mathTypes[i], "mismatch of mathType");
    fail_unless(temp_dataType == dataType, "mismatch of dataType");
    for (j = 0; j < temp_numStepPerSeqVar[i]; j++) {
      rc = fc_getVariableDataPtr(temp_glbSeqVars[i][j], &temp_data);
      fail_unless(rc == FC_SUCCESS, "failed to get glbSeqVar data");
      fail_unless(!memcmp(temp_data, datas[i], sizeof(double)),
		  "mismatch of glbSeqVar data");
    }
    free(temp_glbSeqVars[i]);
  }
  free(temp_numStepPerSeqVar);
  free(temp_glbSeqVars);
  // test meshes
  rc = fc_getMeshes(copy_dataset, &temp_numMesh, &temp_meshes);
  fail_unless(rc == FC_SUCCESS, "should not fail to get meshes");
  fail_unless(temp_numMesh == numMesh, "mismatch of numMesh");
  for (i = 0; i < numMesh; i++) {
    rc = fc_getMeshName(temp_meshes[i], &temp_name);
    fail_unless(rc == FC_SUCCESS, "failed to get mesh name");
    fail_unless(!strcmp(temp_name, mesh_names[i]), "mesh name mismatch");
    free(temp_name);
    rc = fc_getMeshInfo(temp_meshes[i], &temp_topodim, &temp_numDim,
			&temp_numVert, &temp_numElem, &temp_elemType);
    fail_unless(rc == FC_SUCCESS, "failed to get mesh info");
    fail_unless(temp_topodim == topodims[i], "mismatch of topodim");
    fail_unless(temp_numDim == numDim, "mismatch of numdim");
    fail_unless(temp_numVert == numVerts[i], "mismatch of numVert");
    fail_unless(temp_numElem == numElems[i], "mismatch of numElem");
    fail_unless(temp_elemType == elemTypes[i], "mismatch of elemType");
    rc = fc_getMeshCoordsPtr(temp_meshes[i], &temp_coords);
    fail_unless(rc == FC_SUCCESS, "failed to get mesh coords");
    fail_unless(!memcmp(temp_coords, coords[i], numDim*numVerts[i]*sizeof(double)),
		"mismatch of coords");
    rc = fc_getMeshElementConnsPtr(temp_meshes[i], &temp_conns);
    fail_unless(rc == FC_SUCCESS, "failed to get mesh coords");
    fail_unless(!memcmp(temp_conns, conns[i], 
			numVertPerElems[i]*numElems[i]*sizeof(int)),
		"mismatch of conns");
    // test subsets - are only on last mesh
    rc = fc_getSubsets(temp_meshes[i], &temp_numSubset, &temp_subsets);
    fail_unless(rc == FC_SUCCESS, "should not fail to get subsets");
    if (i < numMesh-1) {
      fail_unless(temp_numSubset == 0, "mismatch of numsubset");
    }
    else {
      fail_unless(temp_numSubset == numSubset, "mismatch of numSubset");
      for (j = 0; j < numSubset; j++) {
	rc = fc_getSubsetName(temp_subsets[j], &temp_name);
	fail_unless(rc == FC_SUCCESS, "failed to get subset name");
	fail_unless(!strcmp(temp_name, subset_names[j]), "subset name mismatch");
	free(temp_name);
	rc = fc_getSubsetInfo(temp_subsets[j], &temp_numMember,
			      &temp_maxNumMember, &temp_assoc);
	fail_unless(rc == FC_SUCCESS, "failed to get subset info");
	fail_unless(temp_numMember == numMembers[j], "mismatch of numMember");
	fail_unless(temp_maxNumMember == maxNumMembers[j], 
		    "mismatch of maxNumMember");
	fail_unless(temp_assoc == assocs[j], "mismatch of assoc");
	rc = fc_getSubsetMembersAsArray(temp_subsets[j], &temp_numMember, 
					&temp_memberIDs);
	fail_unless(rc == FC_SUCCESS, "failed to get member ids");
	fail_unless(!memcmp(temp_memberIDs, memberIDs[j], numMembers[j]*sizeof(int)),
		    "mismatch of memberIDs");
	free(temp_memberIDs);
      }
    }
    free(temp_subsets);
    // test vars - are only on last mesh
    rc = fc_getVariables(temp_meshes[i], &temp_numVar, &temp_vars);
    fail_unless(rc == FC_SUCCESS, "should not fail to get vars");
    if (i < numMesh-1) {
      fail_unless(temp_numVar == 0, "mismatch of numvar");
    }
    else {
      fail_unless(temp_numVar == numVar, "mismatch of numVar");
      for (j = 0; j < numVar; j++) {
	rc = fc_getVariableName(temp_vars[j], &temp_name);
	fail_unless(rc == FC_SUCCESS, "failed to get var name");
	fail_unless(!strcmp(temp_name, var_names[j]), "var name mismatch");
	free(temp_name);
	rc = fc_getVariableInfo(temp_vars[j], &temp_numDataPoint,
				&temp_numComp, &temp_assoc, &temp_mathType, 
				&temp_dataType);
	fail_unless(rc == FC_SUCCESS, "failed to get var info");
	fail_unless(temp_numDataPoint == numDataPoints[j], 
		    "mismatch of numDataPoints");
	fail_unless(temp_numComp == numComps[j], "mismatch of numComp");
	fail_unless(temp_assoc == assocs[j], "mismatch of assoc");
	fail_unless(temp_mathType == mathTypes[j], "mismatch of mathType");
	fail_unless(temp_dataType == dataType, "mismatch of dataType");
	rc = fc_getVariableDataPtr(temp_vars[j], &temp_data);
	fail_unless(rc == FC_SUCCESS, "failed to get var data");
	fail_unless(!memcmp(temp_data, datas[j], numDataPoints[j]*sizeof(double)),
		    "mismatch of var data");      }
    }
    free(temp_vars);
    // test seq vars - are only on last mesh & last sequence
    rc = fc_getSeqVariables(temp_meshes[i], &temp_numSeqVar,
			    &temp_numStepPerSeqVar, &temp_seqVars);
    fail_unless(rc == FC_SUCCESS, "should not fail to get seqVars");
    if (i < numMesh-1) {
      fail_unless(temp_numSeqVar == 0, "mismatch of numseqVar");
    }
    else {
      fail_unless(temp_numSeqVar == numSeqVar, "mismatch of numSeqVar");
      for (j = 0; j < numSeqVar; j++) {
	fail_unless(temp_numStepPerSeqVar[j] == numSteps[numSequence-1],
		    "mismatch of numStep & probably the sequence");
	rc = fc_getVariableName(temp_seqVars[j][0], &temp_name);
	fail_unless(rc == FC_SUCCESS, "failed to get seqVar name");
	fail_unless(!strcmp(temp_name, seqVar_names[j]), "seqVar name mismatch");
	free(temp_name);
	rc = fc_getVariableInfo(temp_seqVars[j][0], &temp_numDataPoint,
				&temp_numComp, &temp_assoc, &temp_mathType, 
				&temp_dataType);
	fail_unless(rc == FC_SUCCESS, "failed to get seqVar info");
	fail_unless(temp_numDataPoint == numDataPoints[j], 
		    "mismatch of numDataPoints");
	fail_unless(temp_numComp == numComps[j], "mismatch of numComp");
	fail_unless(temp_assoc == assocs[j], "mismatch of assoc");
	fail_unless(temp_mathType == mathTypes[j], "mismatch of mathType");
	fail_unless(temp_dataType == dataType, "mismatch of dataType");
	for (k = 0; k < temp_numStepPerSeqVar[j]; k++) {
	  rc = fc_getVariableDataPtr(temp_seqVars[j][k], &temp_data);
	  fail_unless(rc == FC_SUCCESS, "failed to get seqVar data");
	  fail_unless(!memcmp(temp_data, datas[j], numDataPoints[j]*sizeof(double)),
		      "mismatch of seqVar data");
	}
	free(temp_seqVars[j]);
      }
    }
    free(temp_numStepPerSeqVar);
    free(temp_seqVars);
  }
  free(temp_meshes);
  

 
  // cleanup
  fc_deleteDataset(copy_dataset);
}
END_TEST 

// This only tests error & boundary conditions
// The main load & write tests are done in the library regression tests
START_TEST(load_write)
{
  FC_ReturnCode rc;
  char *good_file_name = "../data/gen_multivar_seq.ex2";
  char *bad_file_name = "doesnt/exist.ex2";
  char new_name[20] = "blue berry";
  FC_Dataset dataset, newDataset, badDataset = { 999, 999 };
  
  // fc_loadDataset() -- bad args
  rc = fc_loadDataset(NULL, &dataset);
  fail_unless(rc != FC_SUCCESS, "should fail if no filename");
  fail_unless(FC_HANDLE_EQUIV(dataset, FC_NULL_DATASET), 
	      "fail should return NULL");
  rc = fc_loadDataset(good_file_name, NULL);
  fail_unless(rc != FC_SUCCESS, "should fail if no dataset to return");

  // fc_loadDataset() -- should fail if file doesn't exits
  rc = fc_loadDataset(bad_file_name, &dataset);
  fail_unless(rc != FC_SUCCESS, "should fail if dataset doesn't exist");
  fail_unless(FC_HANDLE_EQUIV(dataset, FC_NULL_DATASET), 
	      "fail should return NULL");

  // setup for write test
  rc = fc_createDataset(good_file_name, &dataset);
  fail_unless(rc == FC_SUCCESS, "abort: failed to load dataset");
  rc = fc_createDataset(new_name, &newDataset);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create dataset");

  // fc_writeDataset() -- bad args
  rc = fc_writeDataset(badDataset, good_file_name, FC_FT_EXODUS);
  fail_unless(rc != FC_SUCCESS, "should fail if bad dataset");
  rc = fc_writeDataset(dataset, good_file_name, FC_FT_EXODUS);
  fail_unless(rc != FC_SUCCESS, "should fail if read only dataset");

  // fc_writeDataset() -- should fail if an empty dataset
  rc = fc_writeDataset(newDataset, good_file_name, FC_FT_EXODUS); 
  fail_unless(rc != FC_SUCCESS, "should fail if dataset is empty");

  // cleanup
  fc_deleteDataset(dataset);
}
END_TEST 

// release dataset
// since dataset has no big data, this essentially test recursiveness
START_TEST(release_dataset)
{
  FC_ReturnCode rc;
  int i, j, k;
  FC_Dataset dataset, badDataset = { 999, 999 };
  int numSequence, numMesh, numSubset, numVar, numSeqVar, *numSteps;
  FC_Sequence *sequences;
  FC_Mesh *meshes;
  FC_Subset *subsets;
  FC_Variable *variables, **seqVars;

  // setup
  rc = fc_loadDataset("../data/gen_multivar_seq.ex2", &dataset);
  fail_unless(rc == FC_SUCCESS, "abort: failed to load dataset");


  // --- test that releaseDataset propagates to children

  // sequences
  rc = fc_getSequences(dataset, &numSequence, &sequences);
  fail_unless(rc == FC_SUCCESS && numSequence > 0, 
	      "abort: expected sequences");
  for (i = 0; i < numSequence; i++) {
    void* coords;
    _FC_SeqSlot* seqSlot = _fc_getSeqSlot(sequences[i]);
    fail_unless(seqSlot->coords == NULL, "should start with NULL");
    // force loading
    rc = fc_getSequenceCoordsPtr(sequences[i], &coords); 
    fail_unless(rc == FC_SUCCESS, "abort: failed to get seq coords");
    fail_unless(seqSlot->coords != NULL, "big data should NOT be NULL");
    // test
    rc = fc_releaseDataset(dataset);
    fail_unless(rc == FC_SUCCESS, "failed to release dataset");
    fail_unless(seqSlot->coords == NULL, 
		"after release sequence big data should be NULL");
  }
  free(sequences);

  // meshes - A comprehensive testing of releasing mesh is in checkmesh,
  // this just checks that release gets passed on to mesh
  rc = fc_getMeshes(dataset, &numMesh, &meshes);
  fail_unless(rc == FC_SUCCESS && numMesh > 0, "abort: expected meshes");
  for (i = 0; i < numMesh; i++) {
    double* coords;
    _FC_MeshSlot* meshSlot = _fc_getMeshSlot(meshes[i]);

    fail_unless(meshSlot->coords == NULL && 
		meshSlot->elemToVertConns == NULL &&
		meshSlot->numEdge == -1 && meshSlot->numFace == -1 && 
		meshSlot->edgeToVertConns == NULL &&
		meshSlot->faceTypes == NULL &&
		meshSlot->numVertPerFace == NULL &&
		meshSlot->faceToVertConns == NULL &&
		meshSlot->elemToEdgeConns == NULL &&
		meshSlot->elemToFaceConns == NULL &&
		meshSlot->elemFaceOrients == NULL,
		"should start with NULL");    
    // force loading & building
    rc = fc_getMeshCoordsPtr(meshes[i], &coords);
    fail_unless(rc == FC_SUCCESS, "abort: failed to get mesh coords");
    rc = _fc_buildEdgeConns(meshes[i]);
    fail_unless(rc == FC_SUCCESS, "abort: failed to build edges");
    rc = _fc_buildFaceConns(meshes[i]);
    fail_unless(rc == FC_SUCCESS, "abort: failed to build faces");
    fail_unless(meshSlot->coords != NULL && 
		meshSlot->elemToVertConns != NULL &&
		meshSlot->numEdge != -1 && meshSlot->numFace != -1 && 
		meshSlot->edgeToVertConns != NULL &&
		meshSlot->faceTypes != NULL &&
		meshSlot->numVertPerFace != NULL &&
		meshSlot->faceToVertConns != NULL &&
		meshSlot->elemToEdgeConns != NULL &&
		meshSlot->elemToFaceConns != NULL &&
		meshSlot->elemFaceOrients != NULL,
		"abort: big data should NOT be NULL");

    // test
    rc = fc_releaseDataset(dataset);
    fail_unless(rc == FC_SUCCESS, "failed to release dataset");
    fail_unless(meshSlot->coords == NULL && 
		meshSlot->elemToVertConns == NULL &&
		meshSlot->numEdge == -1 && meshSlot->numFace == -1 && 
		meshSlot->edgeToVertConns == NULL &&
		meshSlot->faceTypes == NULL &&
		meshSlot->numVertPerFace == NULL &&
		meshSlot->faceToVertConns == NULL &&
		meshSlot->elemToEdgeConns == NULL &&
		meshSlot->elemToFaceConns == NULL &&
		meshSlot->elemFaceOrients == NULL,
		"should start with NULL");

    // subsets  
    rc = fc_getSubsets(meshes[i], &numSubset, &subsets);
    fail_unless(rc == FC_SUCCESS && numSubset > 0, 
		"abort: expected some somesets");
    for (j = 0; j < numSubset; j++) {
      _FC_SubSlot* subSlot = _fc_getSubSlot(subsets[j]);
      fail_unless(subSlot->sia.vals != NULL, "exisiting subsets are not lazily loaded");
    }
    //lets cheat and copy a subset and make it committed
    {
      FC_Subset tempSubset;
      _FC_SubSlot* subSlot;
      rc = fc_copySubset(subsets[i], meshes[i], "copysub", &tempSubset);
      fail_unless(rc == FC_SUCCESS, "failed to create subset");
      subSlot = _fc_getSubSlot(tempSubset);
      fail_unless(subSlot->sia.vals != NULL, "exisiting subsets are not lazily loaded");
      subSlot->header.committed = 1; //cheating

      rc = fc_releaseDataset(dataset);
      fail_unless(rc == FC_SUCCESS, "failed to release dataset");
      //all but the committed one should still not be null
      subSlot = _fc_getSubSlot(tempSubset);
      fail_unless(subSlot->sia.vals == NULL, "created subset should be released");
      fc_deleteSubset(tempSubset);
      
      for (j = 0; j < numSubset; j++){
	subSlot = _fc_getSubSlot(subsets[i]);
	fail_unless(subSlot->sia.vals != NULL, "after release these subsets are still not null");
      }
    }
    free(subsets);


    // variables
    rc = fc_getVariables(meshes[i], &numVar, &variables);
    fail_unless(rc == FC_SUCCESS && numVar > 0, 
		"abort: expected some variables"); 
    for (j = 0; j < numVar; j++) {
      void* data;
      char* tempname;
      _FC_VarSlot* varSlot = _fc_getVarSlot(variables[j]);
      rc = fc_getVariableName(variables[j], &tempname);
      fail_unless(rc == FC_SUCCESS, "cant get variable name");
      if (!strcmp(tempname,"ID")){
	//ID var is a whole mesh var and is not lazily loaded
	fail_unless(varSlot->data != NULL, "should not start with NULL");
      } else {
	fail_unless(varSlot->data == NULL, "should start with NULL");
	// force loading
	rc = fc_getVariableDataPtr(variables[j], &data);
	fail_unless(rc == FC_SUCCESS, "abort: failed to get var data");
	fail_unless(varSlot->data != NULL, "big data should NOT be NULL");
      }
      // test
      rc = fc_releaseDataset(dataset);
      fail_unless(rc == FC_SUCCESS, "failed to release dataset");
      if (!strcmp(tempname,"ID")){
	//ID var is a whole mesh var and is not lazily loaded
	fail_unless(varSlot->data != NULL, "after release ID variable big data should not be NULL");
      } else {
	fail_unless(varSlot->data == NULL,
		    "after release variable big data should be NULL");
      }
      free(tempname);
    }
    free(variables);

    // sequence variables
    rc = fc_getSeqVariables(meshes[i], &numSeqVar, &numSteps, &seqVars);
    fail_unless(rc == FC_SUCCESS && numVar > 0, 
		"abort: expected some variables"); 
    for (j = 0; j < numSeqVar; j++) {
      for (k = 0; k < numSteps[j]; k++) {
	void* data;
	_FC_VarSlot* varSlot = _fc_getVarSlot(seqVars[j][k]);
	fail_unless(varSlot->data == NULL, "should start with NULL");
	// force loading
	rc = fc_getVariableDataPtr(seqVars[j][k], &data);
	fail_unless(rc == FC_SUCCESS, "abort: failed to get var data");
	fail_unless(varSlot->data != NULL, "big data should NOT be NULL");
	// test
	rc = fc_releaseDataset(dataset);
	fail_unless(rc == FC_SUCCESS, "failed to release dataset");
	fail_unless(varSlot->data == NULL,
		    "after release seqVar big data should be NULL");
      }
    }
    free(numSteps);
    for (j = 0; j < numSeqVar; j++)
      free(seqVars[j]);
    free(seqVars);
  }
  free(meshes);
		
  // --- special cases
  
  // releasing a null handle does not cause an error
  rc = fc_releaseDataset(FC_NULL_DATASET);
  fail_unless(rc == FC_SUCCESS, "NULL handle should not fail");
  
  // --- errors
  
  // fc_releaseDataset() --- bad args
  rc = fc_releaseDataset(badDataset);
  fail_unless(rc != FC_SUCCESS, "bad dataset should fail");
 
  // cleanup
  fc_deleteDataset(dataset);
}
END_TEST 


// *********************************************
// ***** Populate the Suite with the tests
// *********************************************

Suite *dataset_suite(void)
{
  Suite *suite = suite_create("Dataset");

  TCase *tc_private_ds = tcase_create(" - Private Dataset ");
  TCase *tc_dataset = tcase_create(" - Dataset Interface ");
  TCase *tc_recursive = tcase_create(" - Recursive Dataset Intf.");

  // private dataset
  suite_add_tcase(suite, tc_private_ds);
  tcase_add_checked_fixture(tc_private_ds, dataset_setup, dataset_teardown);
  tcase_add_test(tc_private_ds, slot_new_delete);
  
  // dataset interface
  suite_add_tcase(suite, tc_dataset);
  tcase_add_test(tc_dataset, create_get_delete);
  tcase_add_test(tc_dataset, metadata_query);
  tcase_add_test(tc_dataset, print);

  // dataset calls that recurse through children : seq, mesh, var etc
  suite_add_tcase(suite, tc_recursive);
  tcase_add_checked_fixture(tc_recursive, dataset_setup, dataset_teardown);
  tcase_add_test(tc_recursive, copy_test);
  tcase_add_test(tc_recursive, load_write);
  tcase_add_test(tc_recursive, release_dataset);

  return suite;
}
