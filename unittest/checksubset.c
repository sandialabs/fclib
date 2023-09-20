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
 * \file checksubset.c
 * \brief Unit tests for \ref Subset module.
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/unittest/checksubset.c,v $
 * $Revision: 1.61 $ 
 * $Date: 2006/10/19 03:14:53 $
 *
 * \description
 *
 *     To test the creation and operations on subsets of entities living on
 *     a specified mesh.
 *
 * \modifications
 *    - 03/31/04 RM, created.
 */

#include <stdlib.h>
#include <string.h>
#include "check.h"
#include "fc.h"
#include "fcP.h"
#include "checkall.h"

// name of test database created in test fixtures
static char dataset_name[40] = "fixture created dataset";

// **** test fixtures
// setup: init library & create a dataset with a mesh of each element type
// teardown: final library & delete the dataset

static void subset_setup(void) {
  FC_ReturnCode rc;
  int i;
  FC_Dataset dataset;
  FC_Mesh mesh;
  int numElemType = 8, numElem = 2;
  int numVertPerType[8] = { 2, 3, 4, 6, 5, 6, 8, 12 }; 
  char names[8][40] = { "point mesh", "line mesh", "tri mesh", "quad mesh",
			"tet mesh", "pyramid mesh", "prism mesh", "hex mesh" };
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
    rc = fc_createMesh(dataset, names[i], &mesh);
    fail_unless(rc == FC_SUCCESS, "abort: failed to create mesh");
    rc = fc_setMeshCoords(mesh, 3, numVertPerType[i], coords);
    fail_unless(rc == FC_SUCCESS, "failed to set vertex coords");
    rc = fc_setMeshElementConns(mesh, elemTypes[i], numElem, conns[i]);
    fail_unless(rc == FC_SUCCESS, "failed to set element conns");
  }
}

static void subset_teardown(void) {
  FC_ReturnCode rc;
  FC_Dataset *datasets;
  int numDatasets;

  rc = fc_getDatasetByName(dataset_name, &numDatasets,&datasets);
  fail_unless(rc == FC_SUCCESS && numDatasets == 1,
	      "expected to find fixture dataset");
  rc = fc_deleteDataset(datasets[0]);
  free(datasets);
  fail_unless(rc == FC_SUCCESS,
	      "should not fail to close fixture dataset");

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
  FC_Dataset dataset, *returnDatasets;
  FC_Mesh mesh, *returnMeshes;
  int numReturnMeshes, numReturnDatasets;
  FC_Subset subsets[7];
  
  // get dataset and mesh
  rc = fc_getDatasetByName(dataset_name,&numReturnDatasets,
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
    rc = fc_createSubset(mesh, "fake", FC_AT_VERTEX, &subsets[i]);
    fail_unless(rc == FC_SUCCESS, "failed to create subset");
  }
  uID_start = subsets[0].uID;

  // Delete "odd" slots
  for (i = 1; i < num; i+=2) {
    rc = fc_deleteSubset(subsets[i]);
    fail_unless(rc == FC_SUCCESS, "failed to delete subset");
  }

  // Create some more (fills in the cracks)
  for (i = 1; i < num; i+=2) {
    rc = fc_createSubset(mesh, "fake", FC_AT_VERTEX, &subsets[i]);
    fail_unless(rc == FC_SUCCESS, "failed to create subset");
  }

  // Check slotID and uID
  for (i = 0; i < num; i++) {
    fail_unless(subsets[i].slotID == i, "mismatch of slot id");
    fail_unless(subsets[i].uID == uID_start + uID_incs[i],
		"mismatch of uID");
  }

  // cleanup is done in the fixture
  // cleanup
  //rc = fc_deleteDataset(dataset);
  //fail_unless(rc == FC_SUCCESS, "failed to delete dataset at end of test");
}
END_TEST

// test creating and deleting subsets, also test querying mesh
// for subsets before, inbetween and after.
START_TEST(create_get_delete)
{
  FC_ReturnCode rc;
  int i, j, k;
  FC_Dataset dataset,*returnDatasets;
  FC_Mesh *meshes, mesh;
  int numReturnDatasets, numMesh, numAssoc = 5;
  FC_AssociationType assocs[5] = { FC_AT_VERTEX, FC_AT_EDGE, FC_AT_FACE,
				   FC_AT_ELEMENT, FC_AT_WHOLE_MESH };
  int numSubset = 10, temp_numSubset;
  FC_Subset subsets[10], *temp_subsets, temp_subset;
  FC_Subset badSubset = { 999, 999 };
  char names[10][100] = { "one", "two", "three", "four", "five", "six",
			  "seven", "eight", "nine", "ten" };
  char newName[100] = { "sparkly new" };

  // get dataset and meshes
  rc = fc_getDatasetByName(dataset_name,&numReturnDatasets,
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

    // Expect error to create subset of assoc FC_AT_WHOLE_SUBSET
    rc = fc_createSubset(mesh, names[0], FC_AT_WHOLE_DATASET, &subsets[0]);
    fail_unless(rc != FC_SUCCESS, "should fail to create WHOLE_DATASET sub");
    fail_unless(FC_HANDLE_EQUIV(subsets[0], FC_NULL_SUBSET), 
		"fail should return NULL");

    for (j = 0; j < numAssoc; j++) {
      
      // create first subset (continue if supposed to fail)
      rc = fc_createSubset(mesh, names[0], assocs[j], &subsets[0]);
      if ( (i == 0 && (j == 1 || j == 2)) /* point mesh - no edges,faces */
	   || (i == 1 && j == 2) ) { /* line mesh - no faces */
	fail_unless(rc != FC_SUCCESS, "should fail to create subset");
	continue;
      }
      else {
	fail_unless(rc == FC_SUCCESS, "failed to create subset");
	fail_unless(!FC_HANDLE_EQUIV(subsets[0], FC_NULL_SUBSET),
		    "created subset should not = FC_NULL_SUBSET");
      }

      // create some more meshes
      for (k = 1; k < numSubset; k++) {
	rc = fc_createSubset(mesh, names[k], assocs[j], &subsets[k]);
	fail_unless(rc == FC_SUCCESS, "failed to create subset");
	fail_unless(!FC_HANDLE_EQUIV(subsets[k], FC_NULL_SUBSET),
		    "created subset should not = FC_NULL_SUBSET");
      }

      // get subsets (should be numSubset);
      rc = fc_getSubsets(mesh, &temp_numSubset, &temp_subsets);
      fail_unless(rc == FC_SUCCESS, "failed to get subsets");
      fail_unless(temp_numSubset == numSubset, "mismatch of numSubset");
      for (k = 0; k < numSubset; k++) 
	fail_unless(FC_HANDLE_EQUIV(temp_subsets[k], subsets[k]), 
		    "mismatch of subset handles");
      free(temp_subsets);
      rc = fc_getNumSubset(mesh, &temp_numSubset);
      fail_unless(rc == FC_SUCCESS, "failed to get numSubset");
      fail_unless(temp_numSubset == numSubset, "mismatch of numSubset");
      for (k = 0; k < numSubset; k++) {
	rc = fc_getSubsetByName(mesh, names[k], &temp_numSubset,&temp_subsets);
	fail_unless(rc == FC_SUCCESS, "failed to get subset by name");
	fail_unless(temp_numSubset == 1, "wrong number of matching subsets");
	fail_unless(FC_HANDLE_EQUIV(temp_subsets[0], subsets[k]),
		    "mismatch of subset handles");
	free(temp_subsets);
      }

      //temporarily create a new subset to check array return
      rc = fc_createSubset(mesh, names[0], assocs[j], &temp_subset);
      fail_unless(rc == FC_SUCCESS, "failed to create subset");
      fail_unless(!FC_HANDLE_EQUIV(temp_subset, FC_NULL_SUBSET),
		  "created subset should not = FC_NULL_SUBSET");
      rc = fc_getSubsetByName(mesh, names[0], &temp_numSubset,&temp_subsets);
      fail_unless(rc == FC_SUCCESS, "failed to get subset by name");
      fail_unless(temp_numSubset == 2, "wrong number of matching subsets");
      fail_unless(((FC_HANDLE_EQUIV(temp_subsets[0], subsets[0]) &&
		    FC_HANDLE_EQUIV(temp_subsets[1], temp_subset)) ||
		   (FC_HANDLE_EQUIV(temp_subsets[1], subsets[0]) &&
		    FC_HANDLE_EQUIV(temp_subsets[0], temp_subset))),
		  "mismatch of subset handles");
      free(temp_subsets);
      fc_deleteSubset(temp_subset);

      // change the name of the first subset
      rc = fc_changeSubsetName(subsets[0], newName);
      fail_unless(rc == FC_SUCCESS, "failed to change subset name");
      rc = fc_getSubsetByName(mesh, names[0], &temp_numSubset,&temp_subsets);
      fail_unless(rc == FC_SUCCESS, "old name should work");
      fail_unless(temp_numSubset == 0, "old name shouldn't find anything");
      fail_unless(temp_subsets == NULL, "old name should return NULL");
      rc = fc_getSubsetByName(mesh, newName, &temp_numSubset,&temp_subsets);
      fail_unless(rc == FC_SUCCESS, "new name should work for find");
      fail_unless(temp_numSubset == 1, "wrong number of matching subsets");
      fail_unless(FC_HANDLE_EQUIV(temp_subsets[0], subsets[0]),
		  "mismatch of subset handle");
      free(temp_subsets);


      // delete half of the subsets (alternate)
      for (k = 0; k < numSubset; k+=2) {
	rc = fc_deleteSubset(subsets[k]);
	fail_unless(rc == FC_SUCCESS, "failed to delete subset");
      }
      
      // get subsets (should be numSubset/2);
      fc_getSubsets(mesh, &temp_numSubset, &temp_subsets);
      fail_unless(rc == FC_SUCCESS, "failed to get subsets");
      fail_unless(temp_numSubset == numSubset/2, "mismatch of numSubset");
      for (k = 0; k < numSubset/2; k++) 
	fail_unless(FC_HANDLE_EQUIV(temp_subsets[k], subsets[k*2+1]), 
		    "mismatch of subset handles");
      free(temp_subsets);
      for (k = 0; k < numSubset/2; k++) {
	rc = fc_getSubsetByName(mesh, names[k*2+1], &temp_numSubset,&temp_subsets);
	fail_unless(rc == FC_SUCCESS, "failed to get subset by name");
	fail_unless(temp_numSubset == 1, "wrong number of matching subsets");	
	fail_unless(FC_HANDLE_EQUIV(temp_subsets[0], subsets[k*2+1]),
		    "mismatch of subset handles");
	free(temp_subsets);
      }
      
      // delete remaining subsets
      for (k = 1; k < numSubset; k+=2) { 
	rc = fc_deleteSubset(subsets[k]);
	fail_unless(rc == FC_SUCCESS, "failed to delete subset");
      }
      
      // get subsets (should be none)
      rc = fc_getSubsets(mesh, &temp_numSubset, &temp_subsets);
      fail_unless(rc == FC_SUCCESS, "failed to get subsets from empty library");
      fail_unless(temp_numSubset == 0 && temp_subsets == NULL,
		  "should return 0 if all subsets deleted");

      // make one more subset for further testing
      rc = fc_createSubset(mesh, names[0], assocs[j], &subsets[0]);
      fail_unless(rc == FC_SUCCESS, "aborted: failed to create subset subset");

      // ---- test special cases
      
      // not an error to delete FC_NULL_SUBSET
      rc = fc_deleteSubset(FC_NULL_SUBSET);
      fail_unless(rc == FC_SUCCESS, "should not error to delete NULL subset");
      
      // ---- test error conditions
      
      // bad args to fc_createSubset()
      temp_subset = badSubset;
      rc = fc_createSubset(FC_NULL_MESH, names[0], assocs[j], &temp_subset);
      fail_unless(rc != FC_SUCCESS, 
		  "should have failed to create subset with null database");
      fail_unless(FC_HANDLE_EQUIV(temp_subset, FC_NULL_SUBSET),
		  "fail should return NULL subset");
      temp_subset = badSubset;
      rc = fc_createSubset(mesh, NULL, assocs[j], &temp_subset);
      fail_unless(rc != FC_SUCCESS, 
		  "should have failed to create subset with null name");
      fail_unless(FC_HANDLE_EQUIV(temp_subset, FC_NULL_SUBSET),
		  "fail should return NULL subset");
      temp_subset = badSubset;
      rc = fc_createSubset(mesh, NULL, FC_AT_UNKNOWN, &temp_subset);
      fail_unless(rc != FC_SUCCESS, 
		  "should have failed with unknown assoc");
      fail_unless(FC_HANDLE_EQUIV(temp_subset, FC_NULL_SUBSET),
		  "fail should return NULL subset");
      temp_subset = badSubset;
      rc = fc_createSubset(mesh, NULL, -999, &temp_subset);
      fail_unless(rc != FC_SUCCESS, 
		  "should have failed with invalid assoc");
      fail_unless(FC_HANDLE_EQUIV(temp_subset, FC_NULL_SUBSET),
		  "fail should return NULL subset");
      rc = fc_createSubset(mesh, names[0], assocs[j], NULL);
      fail_unless(rc != FC_SUCCESS, 
		  "should have failed to create subset with null handle");

      // bad args to fc_changeSubsetName()
      rc = fc_changeSubsetName(subsets[0], NULL);
      fail_unless(rc != FC_SUCCESS, "should fail if new name is NULL");
      rc = fc_getSubsetByName(mesh, names[0], &temp_numSubset, &temp_subsets);
      fail_unless(rc == FC_SUCCESS && temp_numSubset == 1 &&
		  FC_HANDLE_EQUIV(subsets[0], temp_subsets[0]),
		  "should not change name of subset");
      free(temp_subsets);

      // bad args to fc_getSubsets()
      temp_numSubset = 99;
      temp_subsets = (FC_Subset*)1;
      rc = fc_getSubsets(FC_NULL_MESH, &temp_numSubset, &temp_subsets);
      fail_unless(rc != FC_SUCCESS, 
		  "should have failed to create subset with null mesh");
      fail_unless(temp_numSubset == -1,
		  "fail should return -1 for numSubset"); 
      fail_unless(temp_subsets == NULL, "fail should return null");
      temp_subsets = (FC_Subset*)1;
      rc = fc_getSubsets(mesh, NULL, &temp_subsets);
      fail_unless(rc != FC_SUCCESS, 
		  "should have failed to get subset with NULL numSubset");
      fail_unless(temp_subsets == NULL, "fail should return null");
      temp_numSubset = 999;
      rc = fc_getSubsets(mesh, &temp_numSubset, NULL);
      fail_unless(rc != FC_SUCCESS, "should error to get just numSubset");
      fail_unless(temp_numSubset == -1, "mismatch of numSubset");
      
      // bad args to fc_getNumSubset()
      temp_numSubset = 99;
      rc = fc_getNumSubset(FC_NULL_MESH, &temp_numSubset);
      fail_unless(rc != FC_SUCCESS, 
		  "should have failed to get numSubset with null mesh");
      fail_unless(temp_numSubset == -1,
		  "fail should return -1 for numSubset"); 
      rc = fc_getNumSubset(mesh, NULL);
      fail_unless(rc != FC_SUCCESS, 
		  "should have failed to get numSubset with NULL numSubset");
      fail_unless(temp_subsets == NULL, "mismatch of subsets");
      
      // bad args to fc_getSubsetByName()
      rc = fc_getSubsetByName(FC_NULL_MESH, names[0], &temp_numSubset,
			      &temp_subsets);
      fail_unless(rc != FC_SUCCESS, 
		  "should have failed to get subset with null mesh");
      fail_unless(temp_numSubset == -1, "fail should return -1");
      fail_unless(temp_subsets == NULL, "fail should return NULL");
      rc = fc_getSubsetByName(mesh, NULL, &temp_numSubset,
			      &temp_subsets);
      fail_unless(rc != FC_SUCCESS, 
		  "should have failed to get subset with null name");
      fail_unless(temp_numSubset == -1, "fail should return -1");
      fail_unless(temp_subsets == NULL, "fail should return NULL");
      rc = fc_getSubsetByName(mesh, names[0], NULL, &temp_subsets);
      fail_unless(rc != FC_SUCCESS, 
		  "should have failed to get subset with null arg");
      fail_unless(temp_subsets == NULL, "fail should return NULL");
      rc = fc_getSubsetByName(mesh, names[0], &temp_numSubset, NULL);
      fail_unless(rc != FC_SUCCESS, 
		  "should have failed to get subset with null name");
      fail_unless(temp_numSubset == -1, "fail should return -1");

      
      // bad args to fc_deleteSubset()
      rc = fc_deleteSubset(badSubset);
      fail_unless(rc != FC_SUCCESS, 
		  "should error to delete nonexistent subset");

      // --- done
      
      // delete last subset
      rc = fc_deleteSubset(subsets[0]);
      fail_unless(rc == FC_SUCCESS, "failed to delete last subset");
    }
  }
  free(meshes);
}
END_TEST

// query meta data
START_TEST(metadata_query)
{
  FC_ReturnCode rc;
  int i, j;
  FC_Dataset dataset, *returnDatasets;
  FC_Mesh *meshes, mesh, temp_mesh, badMesh = {999, 999 };
  int numReturnDatasets, numMesh, numAssoc = 5;
  FC_AssociationType assocs[5] = { FC_AT_VERTEX, FC_AT_EDGE, FC_AT_FACE,
				   FC_AT_ELEMENT, FC_AT_WHOLE_MESH };
  char name[20] = "blue berry", *temp_name;
  FC_Subset subset, badSubset = { 999, 999 };

  // get dataset and meshes
  rc = fc_getDatasetByName(dataset_name,&numReturnDatasets,
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
    for (j = 0; j < numAssoc; j++) {
      
      //printf("mesh = %d, assoc = %s\n",
      //     i,fc_getAssociationTypeText(assocs[j]));

      // skip uncreateable subsets
      if ( (i == 0 && (j == 1 || j == 2)) /* point mesh - no edges,faces */
	   || (i == 1 && j == 2) ) { /* line mesh - no faces */
	continue;
      }

      // create a subset to play with
      rc = fc_createSubset(mesh, name, assocs[j], &subset);
      fail_unless(rc == FC_SUCCESS, "failed to create subset");

      // test fc_isSubsetValid()
      fail_unless(fc_isSubsetValid(subset),
		  "failed to validate a valid subset");

      // test fc_getSubsetName();
      temp_name = NULL;
      rc = fc_getSubsetName(subset, &temp_name);
      fail_unless(rc == FC_SUCCESS, "failed to get subset name");
      fail_unless(!strcmp(name, temp_name), "mismatch of name");
      // should have returned a copy, mess with temp_name & try again to make sure
      temp_name[0] = 'q';
      free(temp_name);
      rc = fc_getSubsetName(subset, &temp_name);
      fail_unless(rc == FC_SUCCESS, "failed to get subset name");
      fail_unless(!strcmp(name, temp_name), "mismatch of name");
      free(temp_name);
      
      // test fc_getMeshFromSubset()
      temp_mesh = badMesh;
      rc = fc_getMeshFromSubset(subset, &temp_mesh);
      fail_unless(rc == FC_SUCCESS, "failed to get parent mesh");
      fail_unless(FC_HANDLE_EQUIV(temp_mesh, mesh), "mismatch of parent mesh");

      // --- check with bad args
      
      // fc_isMeshValid()
      fail_unless(!fc_isSubsetValid(badSubset), 
		  "badSubset should NOT be valid");
      fail_unless(!fc_isSubsetValid(FC_NULL_SUBSET), 
		  "FC_NULL_SUBSET should not be valid");
      
      // fc_getSubsetName()
      temp_name = (char*)1;
      rc = fc_getSubsetName(badSubset, &temp_name);
      fail_unless(rc != FC_SUCCESS, "badSubset should NOT return a name");
      fail_unless(temp_name == NULL, "fail should return NULL name");
      rc = fc_getSubsetName(subset, NULL);
      fail_unless(rc != FC_SUCCESS, "should fail for NULL name");
      
      // fc_getMeshFromSubset()
      temp_mesh = badMesh;
      rc = fc_getMeshFromSubset(badSubset, &temp_mesh);
      fail_unless(rc != FC_SUCCESS, "badSubset should NOT return a mesh");
      fail_unless(FC_HANDLE_EQUIV(temp_mesh, FC_NULL_MESH),
		  "failure should return NULL mesh");
      rc = fc_getMeshFromSubset(subset, NULL);
      fail_unless(rc != FC_SUCCESS, "should fail for NULL mesh");
      
      // --- all done
      
      // cleanup & one last test
      fc_deleteSubset(subset);
      fail_unless(!fc_isSubsetValid(subset), 
		  "handle should not be valid after a delete");
    }
  }
  free(meshes);
}
END_TEST

// add/delete single members and get members as array for
//   all association types
START_TEST(one_member_add_del)
{
  FC_ReturnCode rc;
  int i, k;
  FC_Dataset dataset, *returnDatasets;
  FC_Mesh mesh, *returnMeshes;
  int numReturnDatasets,numReturnMeshes;
  int numAssoc = 5;
  FC_AssociationType assocs[5] = { FC_AT_VERTEX, FC_AT_EDGE, FC_AT_FACE,
				   FC_AT_ELEMENT, FC_AT_WHOLE_MESH };
  char name[20] = "blue berry";
  FC_Subset subset, badSubset = { 999, 999 };
  int maxNumMembers[5] = { 12, 20, 11, 2, 1 };
  int temp_numMember, temp_maxNumMember, *temp_members;
  FC_AssociationType temp_assoc;

  // get dataset and mesh
  rc = fc_getDatasetByName(dataset_name,&numReturnDatasets,
			   &returnDatasets);
  fail_unless(rc == FC_SUCCESS, "aborted - cant get dataset by name");
  fail_unless(numReturnDatasets == 1,
	      "failed to find unique dataset by name");
  dataset = returnDatasets[0];
  free(returnDatasets);
  rc = fc_getMeshByName(dataset, "hex mesh", &numReturnMeshes,
			&returnMeshes);
  fail_unless(rc == FC_SUCCESS, "aborted - cant get mesh by name");
  fail_unless(numReturnMeshes == 1,
	      "failed to find unique mesh by name");
  mesh = returnMeshes[0];
  free(returnMeshes);

  // loop over association types
  for (i = 0; i < numAssoc; i++) {

    // create a subset
    rc = fc_createSubset(mesh, name, assocs[i], &subset);
    fail_unless(rc == FC_SUCCESS, "failed to create subset");
    
    // check default values
    rc = fc_getSubsetInfo(subset, &temp_numMember, &temp_maxNumMember,
			  &temp_assoc);
    fail_unless(rc == FC_SUCCESS, "failed to get subset info");
    fail_unless(temp_numMember == 0 && temp_maxNumMember == maxNumMembers[i] 
		&& temp_assoc == assocs[i], 
		"did not have default subset values");

    // add half of the members (evens) (in reverse order for fun) and check
    for (k = maxNumMembers[i] - 1 - (maxNumMembers[i]+1)%2; k >=  0; k-=2) {
      rc = fc_addMemberToSubset(subset, k);
      fail_unless(rc == FC_SUCCESS, "failed to add member");
    }
    rc = fc_getSubsetInfo(subset, &temp_numMember, &temp_maxNumMember,
			  &temp_assoc);
    fail_unless(rc == FC_SUCCESS, "failed to get subset info");
    fail_unless(temp_numMember == (maxNumMembers[i]+1)/2, 
		"mismatch of numMember");
    fail_unless(temp_maxNumMember == maxNumMembers[i] && 
		temp_assoc == assocs[i],
		"this part of info should not change");
    rc = fc_getSubsetMembersAsArray(subset, &temp_numMember, &temp_members);
    fail_unless(rc == FC_SUCCESS, "failed to get subset info");
    fail_unless(temp_numMember == (maxNumMembers[i]+1)/2, 
		"mismatch of numMember");
    for (k = 0; k < temp_numMember; k++) 
      fail_unless(temp_members[k] == k*2, "mismatch of member");
    free(temp_members);
    
    // check other query functions that are simple forms of fc_getMeshInfo
    rc = fc_getSubsetNumMember(subset, &temp_numMember);
    fail_unless(rc == FC_SUCCESS, "failed to get numMember");
    fail_unless(temp_numMember == (maxNumMembers[i]+1)/2, 
		"mismatch of numMember");
    rc = fc_getSubsetMaxNumMember(subset, &temp_maxNumMember);
    fail_unless(rc == FC_SUCCESS, "failed to get maxNumMember");
    fail_unless(temp_maxNumMember == maxNumMembers[i],
		"mismatch of maxNumMember");
    rc = fc_getSubsetAssociationType(subset, &temp_assoc);
    fail_unless(rc == FC_SUCCESS, "failed to get assoc");
    fail_unless(temp_assoc == assocs[i], "mismatch of assoc");
    
    // add all members (with repeats) and check
    for (k = 0; k < maxNumMembers[i]; k++) {
      rc = fc_addMemberToSubset(subset, k);
      fail_unless(rc == FC_SUCCESS, "failed to add member");
    }
    rc = fc_getSubsetInfo(subset, &temp_numMember, &temp_maxNumMember,
			  &temp_assoc);
    fail_unless(rc == FC_SUCCESS, "failed to get subset info");
    fail_unless(temp_numMember == maxNumMembers[i], "mismatch of numMember");
    fail_unless(temp_maxNumMember == maxNumMembers[i] && 
		temp_assoc == assocs[i],
		"this part of info should not change");
    rc = fc_getSubsetMembersAsArray(subset, &temp_numMember, &temp_members);
    fail_unless(rc == FC_SUCCESS, "failed to get subset info");
    fail_unless(temp_numMember == maxNumMembers[i], "mismatch of numMember");
    for (k = 0; k < temp_numMember; k++) 
      fail_unless(temp_members[k] == k, "mismatch of member");
    free(temp_members);
    
    
    // delete half of the members (evens) and check
    for (k = 0; k < maxNumMembers[i]; k+=2) {
      rc = fc_deleteMemberFromSubset(subset, k);
      fail_unless(rc == FC_SUCCESS, "failed to delete member");
    }
    rc = fc_getSubsetInfo(subset, &temp_numMember, &temp_maxNumMember,
			  &temp_assoc);
    fail_unless(rc == FC_SUCCESS, "failed to get subset info");
    fail_unless(temp_numMember == maxNumMembers[i]/2, 
		"mismatch of numMember");
    fail_unless(temp_maxNumMember == maxNumMembers[i] && 
		temp_assoc == assocs[i],
		"this part of info should not change");
    rc = fc_getSubsetMembersAsArray(subset, &temp_numMember, &temp_members);
    fail_unless(rc == FC_SUCCESS, "failed to get subset info");
    fail_unless(temp_numMember == maxNumMembers[i]/2, 
		"mismatch of numMember");
    if (temp_numMember > 0) {
      for (k = 0; k < temp_numMember; k++) 
	fail_unless(temp_members[k] == k*2+1, "mismatch of member");
      free(temp_members);
    }
    else
      fail_unless(temp_members == NULL, "empty should be null");
    
    // delete all members (with repeats) and check
    for (k = 0; k < maxNumMembers[i]; k++) {
      rc = fc_deleteMemberFromSubset(subset, k);
      fail_unless(rc == FC_SUCCESS, "failed to delete member");
    }
    rc = fc_getSubsetInfo(subset, &temp_numMember, &temp_maxNumMember,
			  &temp_assoc);
    fail_unless(rc == FC_SUCCESS, "failed to get subset info");
    fail_unless(temp_numMember == 0, "mismatch of numMember");
    fail_unless(temp_maxNumMember == maxNumMembers[i] && 
		temp_assoc == assocs[i],
		"this part of info should not change");
    rc = fc_getSubsetMembersAsArray(subset, &temp_numMember, &temp_members);
    fail_unless(rc == FC_SUCCESS, "failed to get subset info");
    fail_unless(temp_numMember == 0, "mismatch of numMember");
    fail_unless(temp_members == NULL, "empty should be null");
    
    // ---- check special cases
    
    // fc_getSubsetInfo() -- most arguments are optional
    temp_numMember = temp_maxNumMember = -999;
    temp_assoc = FC_AT_UNKNOWN;
    rc = fc_getSubsetInfo(subset, &temp_numMember, NULL, NULL);
    fail_unless(rc == FC_SUCCESS, "failed to get info");
    fail_unless(temp_numMember == 0, "mismatch of numMember");
    rc = fc_getSubsetInfo(subset, NULL, &temp_maxNumMember, NULL);
    fail_unless(rc == FC_SUCCESS, "failed to get info");
    fail_unless(temp_maxNumMember == maxNumMembers[i], 
		"mismatch of maxNumMember");
    rc = fc_getSubsetInfo(subset, NULL, NULL, &temp_assoc);
    fail_unless(rc == FC_SUCCESS, "failed to get info");
    fail_unless(temp_assoc == assocs[i], "mismatch of assoc");
    
    // ---- check errors
    
    // fc_getSubsetInfo() -- bad args
    rc = fc_getSubsetInfo(badSubset, &temp_numMember, &temp_maxNumMember,
			  &temp_assoc);
    fail_unless(rc != FC_SUCCESS, "should fail with bad subset");
    fail_unless(temp_numMember == -1 && temp_maxNumMember == -1 &&
		temp_assoc == FC_AT_UNKNOWN,
		"failure should return nulls");
    rc = fc_getSubsetInfo(subset, NULL, NULL, NULL);
    fail_unless(rc != FC_SUCCESS, "can't have all args be NULL");
    
    // functions similar to fc_getSubsetInof() -- bad args
    rc = fc_getSubsetNumMember(badSubset, &temp_numMember);
    fail_unless(rc != FC_SUCCESS, "should fail with bad subset");
    fail_unless(temp_numMember == -1, "failure should return -1");
    rc = fc_getSubsetNumMember(subset, NULL);
    fail_unless(rc != FC_SUCCESS, "should fail if numMember = NULL");
    rc = fc_getSubsetMaxNumMember(badSubset, &temp_maxNumMember);
    fail_unless(rc != FC_SUCCESS, "should fail with bad subset");
    fail_unless(temp_maxNumMember == -1, "failure should return -1");
    rc = fc_getSubsetMaxNumMember(subset, NULL);
    fail_unless(rc != FC_SUCCESS, "should fail if maxNumMember = NULL");
    rc = fc_getSubsetAssociationType(badSubset, &temp_assoc);
    fail_unless(rc != FC_SUCCESS, "should fail with bad subset");
    fail_unless(temp_assoc == FC_AT_UNKNOWN, "fail should return null");
    rc = fc_getSubsetAssociationType(subset, NULL);
    fail_unless(rc != FC_SUCCESS, "should fail with assoc = NULL");
    
    // fc_addMemberToSubset() -- bad args
    rc = fc_addMemberToSubset(badSubset, 0);
    fail_unless(rc != FC_SUCCESS, "should fail with bad subset");
    rc = fc_addMemberToSubset(subset, -1);
    fail_unless(rc != FC_SUCCESS, "should fail to add member out of range");
    rc = fc_addMemberToSubset(subset, maxNumMembers[i]+1);
    fail_unless(rc != FC_SUCCESS, "should fail to add member out of range");
    
    // fc_getSubsetMembersAsArray() -- bad args
    rc = fc_getSubsetMembersAsArray(badSubset, &temp_numMember, &temp_members);
    fail_unless(rc != FC_SUCCESS, "should fail with bad subset");
    fail_unless(temp_numMember == -1 && temp_members == NULL,
		"failure should return nulls");
    rc = fc_getSubsetMembersAsArray(subset, NULL, &temp_members);
    fail_unless(rc != FC_SUCCESS, "should fail if numMember = NULL");
    fail_unless(temp_members == NULL, "failure should return nulls");
    rc = fc_getSubsetMembersAsArray(subset, &temp_numMember, NULL);
    fail_unless(rc != FC_SUCCESS, "should fail if members = NULL");
    fail_unless(temp_numMember == -1, "failure should return nulls");
    
    // fc_deleteMemberFromSubset() -- bad args
    rc = fc_deleteMemberFromSubset(badSubset, 0);
    fail_unless(rc != FC_SUCCESS, "should fail with bad subset");
    rc = fc_deleteMemberFromSubset(subset, -1);
    fail_unless(rc != FC_SUCCESS, "should fail to delete member out of range");
    rc = fc_deleteMemberFromSubset(subset, maxNumMembers[i]+1);
    fail_unless(rc != FC_SUCCESS, "should fail to delete member out of range");
    
    // --- all done
    // cleanup
    fc_deleteSubset(subset);
  }
}
END_TEST

// query getting Subset Members as array and mask
START_TEST(get_members)
{
  FC_ReturnCode rc;
  int i;
  FC_Dataset dataset, *returnDatasets;
  FC_Mesh mesh, *returnMeshes;
  int numReturnDatasets,numReturnMeshes;
  char subsetName[] = "New subset";
  FC_Subset subset, badSubset = { 999, 999 };
  FC_AssociationType assoc = FC_AT_VERTEX;
  int numMember = 5, maxNumMember;
  int members[5] = { 1, 0, 2, 10, 7 };
  int orderedMembers[5] = { 0, 1, 2, 7, 10 };
  int count, *array, *mask, sum;

  // get dataset and mesh
  rc = fc_getDatasetByName(dataset_name,&numReturnDatasets,
			   &returnDatasets);
  fail_unless(rc == FC_SUCCESS, "aborted - cant get dataset by name");
  fail_unless(numReturnDatasets == 1,
	      "failed to find unique dataset by name");
  dataset = returnDatasets[0];
  free(returnDatasets);
  rc = fc_getMeshByName(dataset, "hex mesh", &numReturnMeshes,
			&returnMeshes);
  fail_unless(rc == FC_SUCCESS, "aborted - cant get mesh by name");
  fail_unless(numReturnMeshes == 1,
	      "failed to find unique mesh by name");
  mesh = returnMeshes[0];
  free(returnMeshes);

  fc_getMeshNumVertex(mesh, &maxNumMember);
  fail_unless(maxNumMember == 12, "abort: dataset is not as expected");
        
  // ----------- test all gets on empty subset
  
  // create empty subset
  rc = fc_createSubset(mesh, subsetName, assoc, &subset);
  fail_unless(rc == FC_SUCCESS, "failed to create subset");
  
  // get array
  rc = fc_getSubsetMembersAsArray (subset, &count, &array);
  fail_unless(rc == FC_SUCCESS, "failed to get array from empty subset");
  fail_unless(count == 0 && array == NULL, "empty array should return NULL");
  // get mask
  rc = fc_getSubsetMembersAsMask(subset, &count, &mask);
  fail_unless(rc == FC_SUCCESS, "failed to get mask from empty subset");
  fail_unless(count == maxNumMember, "wrong maxNumMember");
  for(i = 0; i < maxNumMember; i++)
    fail_unless(mask[i] == 0, "empty mask should have all 0's"); 
  free(mask);
  
  // ----- test all gets on subset with some members
  
  // add some members to subset
  for (i = 0; i < 5; i++) { 
    rc = fc_addMemberToSubset(subset,members[i]);
    fail_unless(rc == FC_SUCCESS, 
		"getSubsetMembers: add to subset should not fail");
  }
  
  // get array
  rc = fc_getSubsetMembersAsArray (subset, &count, &array);
  fail_unless(rc == FC_SUCCESS, "should not fail to get members as array");
  fail_unless(count == numMember, "mismatch of numMember"); 
  for (i = 0; i < numMember; i++) { 
    fail_unless( array[i] == orderedMembers[i], "mismatch of member");
  }   
  free(array);
  // get mask
  rc = fc_getSubsetMembersAsMask(subset, &count, &mask);
  fail_unless(rc == FC_SUCCESS, "should not fail to get members as mask");
  fail_unless(count == maxNumMember, "wrong maxNumMember");
  sum = 0;
  for (i = 0; i < maxNumMember; i++)
    if (mask[i] == 1) { 
      fail_unless(i  == orderedMembers[sum], "mismatch of member"); 
	sum++; 
    }
  fail_unless(sum == numMember, "wrong number of members in mask");
  free(mask);
  
  // ---- check errors
  
    // fc_getSubsetMembersAsArray() -- bad args
  rc = fc_getSubsetMembersAsArray(badSubset, &count, &array);
  fail_unless(rc != FC_SUCCESS, "should fail with bad subset");
  fail_unless(count == -1 && array == NULL, "default values");
  rc = fc_getSubsetMembersAsArray(subset, &count, NULL);
  fail_unless(rc != FC_SUCCESS, "get subset members with NULL args bad");
  fail_unless(count == -1, "should return nulls");
  rc = fc_getSubsetMembersAsArray(subset, NULL, &array);
  fail_unless(rc != FC_SUCCESS, "get subset members with NULL args bad");
  fail_unless(array == NULL, "should return nulls");
  
  // fc_getSubsetMembersAsMask() -- bad args
  rc = fc_getSubsetMembersAsMask(badSubset, &count, &mask);
  fail_unless(rc != FC_SUCCESS, "should fail with bad subset");
  fail_unless(count == -1 && mask == NULL, "default values");
  rc = fc_getSubsetMembersAsMask(subset, &count, NULL);
  fail_unless(rc != FC_SUCCESS, "get subset members with NULL args bad");
  fail_unless(count == -1, "should return nulls");
  rc = fc_getSubsetMembersAsMask(subset, NULL, &mask);
  fail_unless(rc != FC_SUCCESS, "get subset members with NULL args bad");
  fail_unless(mask == NULL, "should return nulls");
  
  // ---- cleanup
  rc = fc_deleteSubset(subset);
  fail_unless(rc == FC_SUCCESS, "failed to delete subset");
}
END_TEST


// query subset about members
START_TEST(query_members)
{
  FC_ReturnCode rc;
  int i;
  FC_Dataset dataset, *returnDatasets;
  FC_Mesh mesh, *returnMeshes;
  int numReturnDatasets,numReturnMeshes;
  char subsetName[] = "New subset";
  FC_Subset subset, badSubset = { 999, 999 };
  FC_AssociationType assoc = FC_AT_VERTEX;
  int numMember = 5, maxNumMember;
  int members[5] = { 1, 0, 2, 10, 7 };
  int orderedMembers[5] = { 0, 1, 2, 7, 10 };
  int ret;

  // get dataset and mesh
  rc = fc_getDatasetByName(dataset_name,&numReturnDatasets,
			   &returnDatasets);
  fail_unless(rc == FC_SUCCESS, "aborted - cant get dataset by name");
  fail_unless(numReturnDatasets == 1,
	      "failed to find unique dataset by name");
  dataset = returnDatasets[0];
  free(returnDatasets);
  rc = fc_getMeshByName(dataset, "hex mesh", &numReturnMeshes,
			&returnMeshes);
  fail_unless(rc == FC_SUCCESS, "aborted - cant get mesh by name");
  fail_unless(numReturnMeshes == 1,
	      "failed to find unique mesh by name");
  mesh = returnMeshes[0];
  free(returnMeshes);

  fc_getMeshNumVertex(mesh, &maxNumMember);
  fail_unless(maxNumMember == 12, "abort: dataset is not as expected");
        
  // ----------- test on empty subset
  
  // create empty subset
  rc = fc_createSubset(mesh, subsetName, assoc, &subset);
  fail_unless(rc == FC_SUCCESS, "failed to create subset");
  
  // look for member in empty array
  ret = fc_isMemberInSubset(subset, members[0]);
  fail_unless(ret == 0, "should not find one");
  
  // ----- n subset with some members
  
  // add some members to subset
  for (i = 0; i < 5; i++) { 
    rc = fc_addMemberToSubset(subset, members[i]);
    fail_unless(rc == FC_SUCCESS, 
		"getSubsetMembers: add to subset should not fail");
  }
  
  // look for members that exist in subset
  for (i = 0; i < 5; i++) {
    ret = fc_isMemberInSubset(subset, members[0]);
    fail_unless(ret == 1, "should find one");
  }
  // look for members that are out of range of subset
  ret = fc_isMemberInSubset(subset, -1);
  fail_unless(ret == 0, "should not find out of range");
  ret = fc_isMemberInSubset(subset, maxNumMember);
  fail_unless(ret == 0, "should not find out of range");
  // look for possible members that are not present
  ret = fc_isMemberInSubset(subset, orderedMembers[numMember-1]+1);
  fail_unless(ret == 0, "should not find not present member");
  
  // ---- check errors
  
  // fc_isMemberInSubset -- bad args
  ret = fc_isMemberInSubset(badSubset, 1);
  fail_unless(ret < 0, "should fail");
  
  // ---- cleanup
  rc = fc_deleteSubset(subset);
  fail_unless(rc == FC_SUCCESS, "failed to delete subset");
}
END_TEST


// query adding Subset Members from array and mask
START_TEST(multi_members_add)
{
  FC_ReturnCode rc;
  int i;   
  char subsetName[] = "New subset";
  FC_Subset subset;
  FC_Mesh mesh, *returnMeshes;
  int numReturnDatasets,numReturnMeshes;
  FC_Dataset dataset, *returnDatasets;
  FC_AssociationType assoc = FC_AT_VERTEX;
  int members[6] = {1,0,5,10,7,5};
  int orderedMembers[5] = {0,1,5,7,10};
  int numMembers = 6, numOrderedMembers = 5;
  int members2[4] = {10,2,7,6};
  int ordered2[7] = {0,1,2,5,6,7,10}; // combined with 1
  int numMembers2 = 4, numOrderedMembers2 = 7;
  int *mask1, *mask2, maxNumMember;
  int num, *array;
  int numBad = 2;
  int bad1[2] = { 0, -1 };  // one bad entry makes them all bad
  int bad2[2] = { 0, 10000000 };

  // get dataset and mesh
  rc = fc_getDatasetByName(dataset_name,&numReturnDatasets,
			   &returnDatasets);
  fail_unless(rc == FC_SUCCESS, "aborted - cant get dataset by name");
  fail_unless(numReturnDatasets == 1,
	      "failed to find unique dataset by name");
  dataset = returnDatasets[0];
  free(returnDatasets);
  rc = fc_getMeshByName(dataset, "hex mesh", &numReturnMeshes,
			&returnMeshes);
  fail_unless(rc == FC_SUCCESS, "aborted - cant get mesh by name");
  fail_unless(numReturnMeshes == 1, 
	      "failed to find unique mesh by name");
  mesh = returnMeshes[0];
  free(returnMeshes);

  fc_getMeshNumVertex(mesh, &maxNumMember);
  fail_unless(maxNumMember == 12, "abort: dataset is not as expected");

  // make masks that correspond to the arrays
  // this is out here so we don't do it in every loop
  mask1 = calloc(maxNumMember, sizeof(int));
  mask2 = calloc(maxNumMember, sizeof(int));
  for (i = 0; i < numMembers; i++)
    mask1[members[i]] = 1;
  for (i = 0; i < numMembers2; i++)
    mask2[members2[i]] = 1;

  // ------------ add FC_ST_ARRAY members ---------------
  
  // create new subset
  rc = fc_createSubset(mesh, subsetName, assoc, &subset);
  fail_unless(rc == FC_SUCCESS, "test aborted, failed to create subset");

  // add members to empty subset
  rc = fc_addArrayMembersToSubset(subset, numMembers, members);
  fail_unless(rc == FC_SUCCESS, "add array should not fail");
  rc = fc_getSubsetMembersAsArray(subset, &num, &array);
  fail_unless(rc == FC_SUCCESS, "failed to get members as array");
  fail_unless(num == numOrderedMembers, "wrong numMembers");
  for (i = 0; i < numOrderedMembers; i++)
    fail_unless(array[i] == orderedMembers[i], 
		"wrong value in the subset array");
  free(array);
  
  // add more members to subset
  rc = fc_addArrayMembersToSubset(subset, numMembers2, members2);
  fail_unless(rc == FC_SUCCESS, "adding another array should not fail");
  rc = fc_getSubsetMembersAsArray(subset, &num, &array);
  fail_unless(rc == FC_SUCCESS, "failed to get members as array");
  fail_unless(num == numOrderedMembers2, "wrong numMembers");
  for (i = 0; i < numOrderedMembers2; i++)
    fail_unless(array[i] == ordered2[i], 
		"addSubsetMembers: wrong value in the subset array");
  free(array);
  
  // bad input
  rc = fc_addArrayMembersToSubset(FC_NULL_SUBSET, numMembers, members);
  fail_unless(rc != FC_SUCCESS, "should fail with null subset handle");
  rc = fc_addArrayMembersToSubset(subset, numMembers, NULL);
  fail_unless(rc != FC_SUCCESS, "should fail if don't pass members");
  rc = fc_addArrayMembersToSubset(subset, -1, members);
  fail_unless(rc != FC_SUCCESS, "shouldn't be able to add -1 members");
  rc = fc_addArrayMembersToSubset(subset, numBad, bad1);
  fail_unless(rc != FC_SUCCESS, "should fail if bad members");
  rc = fc_getSubsetMembersAsArray(subset, &num, &array);
  fail_unless(num = numOrderedMembers2, "should not change");
  for (i = 0; i < numOrderedMembers2; i++)
    fail_unless(array[i] == ordered2[i], "mismatch of member");
  free(array);
  rc = fc_addArrayMembersToSubset(subset, numBad, bad2);
  fail_unless(rc != FC_SUCCESS, "should fail if bad members");
  rc = fc_getSubsetMembersAsArray(subset, &num, &array);
  fail_unless(num = numOrderedMembers2, "should not change");    
  for (i = 0; i < numOrderedMembers2; i++)
    fail_unless(array[i] == ordered2[i], "mismatch of members");
  free(array);
  
  // o.k. input - zero members
  rc = fc_addArrayMembersToSubset(subset, 0, NULL);
  fail_unless(rc == FC_SUCCESS, "can add 0 members");
  rc = fc_addArrayMembersToSubset(subset, 0, members);
  fail_unless(rc == FC_SUCCESS, "can add 0 members");
  rc = fc_getSubsetMembersAsArray(subset, &num, &array);
  fail_unless(num = numOrderedMembers2, "should not change");    
  for (i = 0; i < numOrderedMembers2; i++)
    fail_unless(array[i] == ordered2[i], "mismatch of members");
  free(array);
  
  // delete subset
  fc_deleteSubset(subset);
  
  // ---------------- add FC_ST_MASK members --------------------

  // create new subset
  rc = fc_createSubset(mesh, subsetName, assoc, &subset);
  fail_unless(rc == FC_SUCCESS, "failed to create another subset");

  // add members to empty subset
  rc = fc_addMaskMembersToSubset(subset, maxNumMember, mask1);
  fail_unless(rc == FC_SUCCESS, "failed to add mask members to subset");
  rc = fc_getSubsetMembersAsArray(subset, &num, &array);
  fail_unless(rc == FC_SUCCESS, "failed to get subset members as array");
  fail_unless(num == numOrderedMembers, "num members wrong");
  for (i = 0; i < numOrderedMembers; i++)
    fail_unless(array[i] == orderedMembers[i], "wrong members");
  free(array);
  
  // add more members to subet, make sure we are testing adding to same type
  rc = fc_addMaskMembersToSubset(subset, maxNumMember, mask2);
  fail_unless(rc == FC_SUCCESS, "failed to add mask members to subset");
  rc = fc_getSubsetMembersAsArray(subset, &num, &array);
  fail_unless(rc == FC_SUCCESS, "failed to get subset members as array");
  fail_unless(num == numOrderedMembers2, "num members wrong");
  for (i = 0; i < numOrderedMembers2; i++)
    fail_unless(array[i] == ordered2[i], "wrong members");
  free(array);
  
  // bad input
  rc = fc_addMaskMembersToSubset(FC_NULL_SUBSET, maxNumMember, mask1);
  fail_unless(rc != FC_SUCCESS, "should fail with null subset handle");
  rc = fc_addMaskMembersToSubset(subset, -1, mask1);
  fail_unless(rc != FC_SUCCESS, "can't add negative");
  rc = fc_addMaskMembersToSubset(subset, 1, mask1);
  fail_unless(rc != FC_SUCCESS, "can't have mismatch in maxNumMember");
  
  // o.k. input - zero members
  rc = fc_addArrayMembersToSubset(subset, 0, NULL);
  fail_unless(rc == FC_SUCCESS, "can add 0 members");
  rc = fc_addArrayMembersToSubset(subset, 0, mask1);
  fail_unless(rc == FC_SUCCESS, "can add 0 members");
  rc = fc_getSubsetMembersAsArray(subset, &num, &array);
  fail_unless(num = numOrderedMembers2, "should not change");    
  for (i = 0; i < numOrderedMembers2; i++)
    fail_unless(array[i] == ordered2[i], 
		"addSubsetMembers: wrong value in the subset array");
  free(array); 
  
  // delete subset
  rc = fc_deleteSubset(subset);
  fail_unless(rc == FC_SUCCESS, "failed to added deleted subset");

  // --------------- done ---------------------- 

  // cleanup
  free(mask1);
  free(mask2);
}
END_TEST


// test delete groups of members from a subset
START_TEST(multi_members_del)
{
  FC_ReturnCode rc;
  int i;
  FC_Dataset dataset, *returnDatasets;
  FC_Mesh mesh, *returnMeshes;
  int numReturnDatasets, numReturnMeshes;
  char subsetName[] = "New subset";
  FC_Subset subset;
  FC_AssociationType assoc = FC_AT_VERTEX;
  int numMembers = 5, numMembers2 = 2;
  int members[5] = { 0, 1, 5, 7, 10 };
  int members2[5] = { 1, 7 };
  // delte from front, then back, then a non member,then from middle
  int delMembers[4] = { 0, 10, 6, 5 }, numDelMembers = 4;
  int bad1[3] = {1, 0, -1}, numBad = 3; // one bad member makes all bad
  int bad2[3] = {1, 0, 100000 };
  int *del_mask1, *del_mask2, maxNumMember;
  int num, *array;

  // get dataset and mesh
  rc = fc_getDatasetByName(dataset_name,&numReturnDatasets,
			   &returnDatasets);
  fail_unless(rc == FC_SUCCESS, "aborted - cant get dataset by name");
  fail_unless(numReturnDatasets == 1,
	      "failed to find unique dataset by name");
  dataset = returnDatasets[0];
  free(returnDatasets);
  rc = fc_getMeshByName(dataset, "hex mesh", &numReturnMeshes,
			&returnMeshes);
  fail_unless(rc == FC_SUCCESS, "aborted - cant get mesh by name");
  fail_unless(numReturnMeshes == 1,
	      "failed to find unique mesh by name");
  mesh = returnMeshes[0];
  free(returnMeshes);

  fc_getMeshNumVertex(mesh, &maxNumMember);
  fail_unless(maxNumMember == 12, "abort: dataset is not as expected");

  // make masks that correspond to the arrays
  // this is out here so we don't do it in every loop
  del_mask1 = calloc(maxNumMember, sizeof(int));
  del_mask2 = calloc(maxNumMember, sizeof(int));
  for (i = 0; i < numDelMembers; i++)
    del_mask1[delMembers[i]] = 1;
  for (i = 0; i < numMembers; i++)
    del_mask2[members[i]] = 1;

  // setup a new subset to play with
  rc = fc_createSubset(mesh, subsetName, assoc, &subset);
  fail_unless(rc == FC_SUCCESS, "failed to create subset");
  rc = fc_addArrayMembersToSubset(subset, numMembers, members);
  fail_unless(rc == FC_SUCCESS, "failed to add members to subset");
  rc = fc_getSubsetMembersAsArray(subset, &num, &array);
  fail_unless(rc == FC_SUCCESS, "failed to get members as array");
  fail_unless(num == numMembers, "wrong numMembers");
  for (i = 0; i < numMembers; i++)
    fail_unless(array[i] == members[i], 
		"addSubsetMembers: wrong value in the subset array");
  free(array);
  
  // --------------- delete FC_ST_ARRAY members
  
  // delete some members
  rc = fc_deleteArrayMembersFromSubset(subset,numDelMembers, delMembers);
  fail_unless(rc == FC_SUCCESS, 
	      "delFromSubset: deleting from subset should not fail 1");
  rc = fc_getSubsetMembersAsArray(subset, &num, &array);
  fail_unless(rc == FC_SUCCESS, "failed to get members as array");
  fail_unless(num == numMembers2, "wrong numMembers");
  for (i = 0; i < numMembers2; i++)
    fail_unless(array[i] == members2[i],
		"addSubsetMembers: wrong value in the subset array");
  free(array);
  
  // delete all members
  rc = fc_deleteArrayMembersFromSubset(subset, numMembers, members);
  fail_unless(rc == FC_SUCCESS, 
	      "delFromSubset: deleting from subset should not fail 1");
  rc = fc_getSubsetMembersAsArray(subset, &num, &array);
  fail_unless(rc == FC_SUCCESS, "failed to get members as array");
  fail_unless(num == 0, "wrong numMembers");
  fail_unless(array == NULL, 
	      "addSubsetMembers: wrong value in the subset array"); 
  
  // add members to this subset, again
  rc = fc_addArrayMembersToSubset(subset, numMembers, members);
  fail_unless(rc == FC_SUCCESS, "failed to add members to subset");
  
  // bad input: deleting members with negative #
  rc = fc_deleteArrayMembersFromSubset(subset, numBad, bad1);
  fail_unless(rc == FC_INPUT_ERROR, 
	      "delFromSubset: deleting from subset should fail 1");
  rc = fc_getSubsetMembersAsArray(subset, &num, &array);
  fail_unless(rc == FC_SUCCESS, "failed to get members as array");
  fail_unless(num == numMembers, "wrong numMembers");
  for (i = 0; i < numMembers; i++)
    fail_unless(array[i] == members[i], 
		"addSubsetMembers: wrong value in the subset array");
  free(array);
  
  // fail: deleting members with too big #
  rc = fc_deleteArrayMembersFromSubset(subset, numBad, bad2);
  fail_unless(rc == FC_INPUT_ERROR,
	      "delFromSubset: deleting from subset should fail 1");
  rc = fc_getSubsetMembersAsArray(subset, &num, &array);
  fail_unless(rc == FC_SUCCESS, "failed to get members as array");
  fail_unless(num == numMembers, "wrong numMembers");
  for (i = 0; i < numMembers; i++)
    fail_unless(array[i] == members[i], 
		"addSubsetMembers: wrong value in the subset array");
  free(array);
  
  // bad input
  rc = fc_deleteArrayMembersFromSubset(FC_NULL_SUBSET, numMembers, members);
  fail_unless(rc != FC_SUCCESS, "should fail with null subset handle");
  rc = fc_deleteArrayMembersFromSubset(subset, numMembers, NULL);
  fail_unless(rc != FC_SUCCESS, "should fail if don't pass members");
  rc = fc_deleteArrayMembersFromSubset(subset, -1, members);
  fail_unless(rc != FC_SUCCESS, "shouldn't be able to add -1 members");
  
  // o.k. input - zero members
  rc = fc_deleteArrayMembersFromSubset(subset, 0, NULL);
  fail_unless(rc == FC_SUCCESS, "can delete 0 members");
  rc = fc_deleteArrayMembersFromSubset(subset, 0, members);
  fail_unless(rc == FC_SUCCESS, "can delete 0 members");
  rc = fc_getSubsetMembersAsArray(subset, &num, &array);
  fail_unless(num = numMembers, "should not change");    
  for (i = 0; i < numMembers; i++)
    fail_unless(array[i] == members[i], 
		"addSubsetMembers: wrong value in the subset array");
  free(array);
  
  // --------------- delete FC_ST_MASK members
  
  // delete some members
  rc = fc_deleteMaskMembersFromSubset(subset, maxNumMember, del_mask1);
  fail_unless(rc == FC_SUCCESS, 
	      "delFromSubset: deleting from subset should not fail 1");
  rc = fc_getSubsetMembersAsArray(subset, &num, &array);
  fail_unless(rc == FC_SUCCESS, "failed to get members as array");
  fail_unless(num == numMembers2, "wrong numMembers");
  for (i = 0; i < numMembers2; i++)
    fail_unless(array[i] == members2[i], 
		"addSubsetMembers: wrong value in the subset array");
  free(array);
  
  // delete all members
  rc = fc_deleteMaskMembersFromSubset(subset, maxNumMember, del_mask2);
  fail_unless(rc == FC_SUCCESS, 
	      "delFromSubset: deleting from subset should not fail 1");
  rc = fc_getSubsetMembersAsArray(subset, &num, &array);
  fail_unless(rc == FC_SUCCESS, "failed to get members as array");
  fail_unless(num == 0, "wrong numMembers");
  fail_unless(array == NULL, 
	      "addSubsetMembers: wrong value in the subset array"); 
  
  // add members to this subset, again
  rc = fc_addArrayMembersToSubset(subset, numMembers, members);
  fail_unless(rc == FC_SUCCESS, "failed to add members to subset");
  
  // can't use masks to add negative or too big numbers
  
  // bad input
  rc = fc_deleteMaskMembersFromSubset(FC_NULL_SUBSET, maxNumMember, del_mask1);
  fail_unless(rc != FC_SUCCESS, "should fail with null subset handle");
  rc = fc_deleteMaskMembersFromSubset(subset, -1, del_mask1);
  fail_unless(rc != FC_SUCCESS, "can't delete negative");
  rc = fc_deleteMaskMembersFromSubset(subset, 1, del_mask1);
  fail_unless(rc != FC_SUCCESS, "can't have mismatch in maxNumMember");
  
  // o.k. input - zero members
  rc = fc_deleteArrayMembersFromSubset(subset, 0, NULL);
  fail_unless(rc == FC_SUCCESS, "can delete 0 members");
  rc = fc_deleteArrayMembersFromSubset(subset, 0, del_mask1);
  fail_unless(rc == FC_SUCCESS, "can delete 0 members");    
  rc = fc_getSubsetMembersAsArray(subset, &num, &array);
  fail_unless(num = numMembers, "should not change");    
  for (i = 0; i < numMembers; i++)
    fail_unless(array[i] == members[i], 
		"addSubsetMembers: wrong value in the subset array");
  free(array);
    
  // -------------- done ---------------------------
  
  // delete subset
  rc = fc_deleteSubset(subset);
  fail_unless(rc == FC_SUCCESS, "failed to delete deleted subset");

  free(del_mask1);
  free(del_mask2);
}
END_TEST

// copy
START_TEST(copy_test)
{
  FC_ReturnCode rc;
  char name[20] = "blue berry", copy_name[20] = "banana nut", *temp_name;
  FC_Dataset dataset, *returnDatasets;
  FC_Mesh mesh1, mesh2, weirdMesh, *returnMeshes;
  int numReturnDatasets,numReturnMeshes;
  int numSubset1, numSubset2, temp_numSubset;
  FC_Subset subset, copy_subset, badSubset = { 999, 999 };
  int numMember = 5, maxNumMember, members[5] = { 0, 1, 5, 7, 10 };
  int temp_numMember, temp_maxNumMember, *temp_members;
  FC_AssociationType assoc = FC_AT_VERTEX, temp_assoc;

  // get dataset and mesh
  rc = fc_getDatasetByName(dataset_name,&numReturnDatasets,
			   &returnDatasets);
  fail_unless(rc == FC_SUCCESS, "aborted - cant get dataset by name");
  fail_unless(numReturnDatasets == 1,
	      "failed to find unique dataset by name");
  dataset = returnDatasets[0];
  free(returnDatasets);
  rc = fc_getMeshByName(dataset, "hex mesh", &numReturnMeshes,
			&returnMeshes);
  fail_unless(rc == FC_SUCCESS, "aborted - cant get mesh by name");
  fail_unless(numReturnMeshes == 1,
	      "failed to find unique mesh by name");
  mesh1 = returnMeshes[0];
  free(returnMeshes);

  fc_getMeshNumVertex(mesh1, &maxNumMember);
  fail_unless(maxNumMember == 12, "abort: dataset is not as expected");
  
  // make second mesh to copy to
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
  
  numSubset1 = 0;
  numSubset2 = 0;
    
  // create a subset & add members
  rc = fc_createSubset(mesh1, name, assoc, &subset);
  fail_unless(rc == FC_SUCCESS, "failed to create subset");
  rc = fc_addArrayMembersToSubset(subset, numMember, members);
  fail_unless(rc == FC_SUCCESS, "abort: failed to creat start subset");
  numSubset1++;
  
  // fc_copySubset() to same mesh & check it
  rc = fc_copySubset(subset, mesh1, copy_name, &copy_subset);
  fail_unless(rc == FC_SUCCESS, "failed to copy to same mesh");
  numSubset1++;
  rc = fc_getSubsetName(copy_subset, &temp_name);
  fail_unless(rc == FC_SUCCESS, "failed to get name of copy");
  fail_unless(!strcmp(temp_name, copy_name), "mismatch of name");
  free(temp_name);
  rc = fc_getSubsetInfo(copy_subset, &temp_numMember, &temp_maxNumMember,
			&temp_assoc);
  fail_unless(rc == FC_SUCCESS, "failed to get subset info");
  fail_unless(temp_numMember == numMember, "mismatch of numMember");
  fail_unless(temp_maxNumMember == maxNumMember, 
	      "mismatch of maxNumMember"); 
  fail_unless(temp_assoc == assoc, "mismatch of assoc");
  rc = fc_getSubsetMembersAsArray(copy_subset, &temp_numMember, 
				  &temp_members);
  fail_unless(rc == FC_SUCCESS, "failed to get subset info");
  fail_unless(temp_numMember == numMember, "mismatch of numMember");
  fail_unless(!memcmp(temp_members, members, sizeof(int)*numMember),
	      "mismatch of members");
  free(temp_members);
  
  // fc_copySubset() to different mesh & check it
  rc = fc_copySubset(subset, mesh2, copy_name, &copy_subset);
  fail_unless(rc == FC_SUCCESS, "failed to copy to different mesh");
  numSubset2++;
  rc = fc_getSubsetName(copy_subset, &temp_name);
  fail_unless(rc == FC_SUCCESS, "failed to get name of copy");
  fail_unless(!strcmp(temp_name, copy_name), "mismatch of name");
  free(temp_name);
  rc = fc_getSubsetInfo(copy_subset, &temp_numMember, &temp_maxNumMember,
			&temp_assoc);
  fail_unless(rc == FC_SUCCESS, "failed to get subset info");
  fail_unless(temp_numMember == numMember, "mismatch of numMember");
  fail_unless(temp_maxNumMember == maxNumMember, 
	      "mismatch of maxNumMember"); 
  fail_unless(temp_assoc == assoc, "mismatch of assoc");
  rc = fc_getSubsetMembersAsArray(copy_subset, &temp_numMember, &temp_members);
  fail_unless(rc == FC_SUCCESS, "failed to get subset info");
  fail_unless(temp_numMember == numMember, "mismatch of numMember");
  fail_unless(!memcmp(temp_members, members, sizeof(int)*numMember),
	      "mismatch of members");
  free(temp_members);
  
  // fc_copySubset() to a very different mesh should fail
  rc = fc_copySubset(subset, weirdMesh, copy_name, &copy_subset);
  fail_unless(rc != FC_SUCCESS, 
	      "should fail to copy to very different mesh");
  fail_unless(FC_HANDLE_EQUIV(copy_subset, FC_NULL_SUBSET),
	      "fail should return NULL handle");
  
  // --- special cases
  
  // copy with NULL name will use original's name
  rc = fc_copySubset(subset, mesh1, NULL, &copy_subset);
  fail_unless(rc == FC_SUCCESS, "should not fail");
  rc = fc_getSubsetName(copy_subset, &temp_name);
  fail_unless(rc == FC_SUCCESS, "failed to get name of copy");
  fail_unless(!strcmp(temp_name, name), "mismatch of name");
  free(temp_name);
  rc = fc_deleteSubset(copy_subset);
  fail_unless(rc == FC_SUCCESS, "failed to delete copied subset");
  
  // ---- check errors
  
  // fc_copySubset() -- bad args
  rc = fc_copySubset(badSubset, mesh1, copy_name, &copy_subset);
  fail_unless(rc != FC_SUCCESS, "should fail to copy bad subset");
  fail_unless(FC_HANDLE_EQUIV(copy_subset, FC_NULL_SUBSET),
	      "fail should return NULL");
  rc = fc_copySubset(subset, FC_NULL_MESH, copy_name, &copy_subset);
  fail_unless(rc != FC_SUCCESS, "should fail to copy to bad mesh");
  fail_unless(FC_HANDLE_EQUIV(copy_subset, FC_NULL_SUBSET),
	      "fail should return NULL");
  rc = fc_copySubset(subset, mesh1, copy_name, NULL);
  fail_unless(rc != FC_SUCCESS, "should fail to copy if NULL handle");

  // a little more checking & cleanup
  rc = fc_deleteMesh(weirdMesh);
  fail_unless(rc == FC_SUCCESS, "failed to delete weird mesh");
  rc = fc_getNumSubset(mesh2, &temp_numSubset);
  fail_unless(rc == FC_SUCCESS, "shouldn't fail to get subsets");
  fail_unless(temp_numSubset == numSubset2, "mismatch of numSubset");
  rc = fc_deleteMesh(mesh2);
  fail_unless(rc == FC_SUCCESS, "failed to delete copy mesh");
  rc = fc_getNumSubset(mesh1, &temp_numSubset);
  fail_unless(rc == FC_SUCCESS, "shouldn't fail to get subset");
  fail_unless(temp_numSubset == numSubset1, "mismatch of numSubset");

  // delete dataset is called in the fixture
}
END_TEST
  
// copy while changing association type
// This assumes fc_changeMeshEntityType() has been tested and is correct
START_TEST(change_assoc)
{
  FC_ReturnCode rc;
  int i, j, k;
  char name[20] = "blue berry", copy_name[20] = "banana nut", *temp_name;
  FC_Dataset dataset, *returnDatasets;
  FC_Mesh mesh, *returnMeshes;
  int numReturnDatasets, numReturnMeshes;
  FC_Subset subset, copy_subset, badSubset = { 999, 999 };
  int numMember, totalNumMember = 5, members[5] = { 0, 1, 5, 7, 10 };
  int numChangedMember, *changedMembers;
  int numAssocType = 5;
  FC_AssociationType assocs[5] = { FC_AT_VERTEX, FC_AT_EDGE, FC_AT_FACE,
				   FC_AT_ELEMENT, FC_AT_WHOLE_MESH };
  int maxNumMember, temp_maxNumMember, temp_numMember, *temp_members;
  FC_AssociationType temp_assoc;

  // get dataset and mesh
  rc = fc_getDatasetByName(dataset_name,&numReturnDatasets,
			   &returnDatasets);
  fail_unless(rc == FC_SUCCESS, "aborted - cant get dataset by name");
  fail_unless(numReturnDatasets == 1,
	      "failed to find unique dataset by name");
  dataset = returnDatasets[0];
  free(returnDatasets);

  rc = fc_getMeshByName(dataset, "hex mesh", &numReturnMeshes,
			&returnMeshes);
  fail_unless(rc == FC_SUCCESS, "aborted - cant get mesh by name");
  fail_unless(numReturnMeshes == 1,
	      "failed to find unique mesh by name");
  mesh = returnMeshes[0];
  free(returnMeshes);

  fc_getMeshNumVertex(mesh, &maxNumMember);
  fail_unless(maxNumMember == 12, "abort: dataset is not as expected");
  
  // for each combination of association types, for each strictness flag
  for (i = 0; i < numAssocType; i++) {

    // create a subset & add members
    rc = fc_createSubset(mesh, name, assocs[i], &subset);
    fail_unless(rc == FC_SUCCESS, "failed to create subset");
    if (assocs[i] == FC_AT_ELEMENT || assocs[i] == FC_AT_WHOLE_MESH) 
      numMember = 1;
    else
      numMember = totalNumMember;
    rc = fc_addArrayMembersToSubset(subset, numMember, members);
    fail_unless(rc == FC_SUCCESS, "abort: failed to create start subset");
    
    for (j = 0; j < numAssocType; j++) {
      for (k = 0; k < 2; k++) { // this is doStrict flag

	// do it
	rc = fc_copySubsetWithNewAssociation(subset, copy_name, assocs[j], k,
					     &copy_subset);
	fail_unless(rc == FC_SUCCESS, "failed to copy");
	rc = fc_getSubsetName(copy_subset, &temp_name);
	fail_unless(rc == FC_SUCCESS, "failed to get name of copy");
	fail_unless(!strcmp(temp_name, copy_name), "mismatch of name");
	free(temp_name);
    
	// test answers
	fc_getMeshNumEntity(mesh, assocs[j], &maxNumMember);
	rc = fc_changeMeshEntityType(mesh, assocs[i], numMember, members,
				     assocs[j], k, &numChangedMember,
				     &changedMembers);
	rc = fc_getSubsetInfo(copy_subset, &temp_numMember, &temp_maxNumMember,
			      &temp_assoc);
	fail_unless(rc == FC_SUCCESS, "failed to get subset info");
	fail_unless(temp_numMember == numChangedMember, 
		    "mismatch of numMember");
	fail_unless(temp_maxNumMember == maxNumMember, 
		    "mismatch of maxNumMember");
	fail_unless(temp_assoc == assocs[j], "mismatch of assoc");
	rc = fc_getSubsetMembersAsArray(copy_subset, &temp_numMember, 
					&temp_members);
	fail_unless(rc == FC_SUCCESS, "failed to get subset info");
	fail_unless(temp_numMember == numChangedMember, 
		    "mismatch of numMember");
	fail_unless(!memcmp(temp_members, changedMembers, 
		    sizeof(int)*numChangedMember), "mismatch of members");
	free(temp_members);
      
	// cleanup
	rc = fc_deleteSubset(copy_subset);
	fail_unless(rc == FC_SUCCESS, "failed to delete copied subset");
	free(changedMembers);

	// --- special cases
	
	// copy with NULL name will use original's name
	rc = fc_copySubsetWithNewAssociation(subset, NULL, assocs[j], k,
					     &copy_subset);
	fail_unless(rc == FC_SUCCESS, "should not fail");
	rc = fc_getSubsetName(copy_subset, &temp_name);
	fail_unless(rc == FC_SUCCESS, "failed to get name of copy");
	fail_unless(!strcmp(temp_name, name), "mismatch of name");
	free(temp_name);
	fc_deleteSubset(copy_subset);

	// ---- check errors
	
	// fc_copySubsetWithNewAssociation() -- bad args
	rc = fc_copySubsetWithNewAssociation(badSubset, copy_name, assocs[j],
					     k, &copy_subset);
	fail_unless(rc != FC_SUCCESS, "should fail to copy bad subset");
	fail_unless(FC_HANDLE_EQUIV(copy_subset, FC_NULL_SUBSET),
		    "fail should return NULL");
	rc = fc_copySubsetWithNewAssociation(subset, copy_name, FC_AT_UNKNOWN,
					     k, &copy_subset);
	fail_unless(rc != FC_SUCCESS, "should fail if unknown assoc");
	fail_unless(FC_HANDLE_EQUIV(copy_subset, FC_NULL_SUBSET),
		    "fail should return NULL");
	rc = fc_copySubsetWithNewAssociation(subset, copy_name, -1, k,
					     &copy_subset);
	fail_unless(rc != FC_SUCCESS, "should fail to copy bad assoc");
	fail_unless(FC_HANDLE_EQUIV(copy_subset, FC_NULL_SUBSET),
		    "fail should return NULL");
	rc = fc_copySubsetWithNewAssociation(subset, copy_name, assocs[j], -1,
					     &copy_subset);
	fail_unless(rc != FC_SUCCESS, "should fail to copy if bad doStrict");
	fail_unless(FC_HANDLE_EQUIV(copy_subset, FC_NULL_SUBSET),
		    "fail should return NULL");
	rc = fc_copySubsetWithNewAssociation(subset, copy_name, assocs[j], 2,
					     &copy_subset);
	fail_unless(rc != FC_SUCCESS, "should fail to copy if bad doStrict");
	fail_unless(FC_HANDLE_EQUIV(copy_subset, FC_NULL_SUBSET),
		    "fail should return NULL");
	rc = fc_copySubsetWithNewAssociation(subset, copy_name, assocs[j], k,
					     NULL);
	fail_unless(rc != FC_SUCCESS, "should fail to copy if NULL handle");
      }
    }
  }

  // delete dataset is called in the fixture
}
END_TEST

START_TEST(release_subset)
{
  FC_ReturnCode rc;
  FC_Dataset dataset, newDataset;
  FC_Mesh mesh, *meshes;
  int numMeshes, numReturnSubsets;
  _FC_SubSlot* subSlot;
  int numMember, *members, wholemembers, compMember;
  int numAssoc = 5;
  FC_AssociationType assocs[5] = { FC_AT_VERTEX, FC_AT_EDGE, FC_AT_FACE,
				   FC_AT_ELEMENT, FC_AT_WHOLE_MESH };
  FC_Subset subsets[5], *returnSubsets, tempSubset, badSubset = { 999, 999 };
  char* subsetNames[5] = { "every 5 verts", "every 5 edges", 
			   "every 5 faces", "every 5 elements", 
			   "the whole thing" };
  FC_Coords lowers = { 0., 0., 0.};
  FC_Coords uppers = { 1., 1., 1.};
  char* filename = "junk.ex2";
  int i, j;

  // setup
  rc = fc_createDataset("dataset", &dataset);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create dataset");
  rc = fc_createSimpleHexMesh(dataset, "mesh", 10, 10 ,10, lowers, uppers, &mesh);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create mesh");
  for (i = 0; i < numAssoc; i++) {
    rc = fc_createSubset(mesh, subsetNames[i], assocs[i], &subsets[i]);
    fail_unless( rc == FC_SUCCESS, "failed to create subset");
    if (assocs[i] == FC_AT_VERTEX){
      rc = fc_getMeshNumVertex(mesh, &numMember);
      fail_unless(rc == FC_SUCCESS, "failed to get num member");
    } else if (assocs[i] == FC_AT_EDGE){
      rc = fc_getMeshNumEdge(mesh, &numMember);
      fail_unless(rc == FC_SUCCESS, "failed to get num member");
    } else if (assocs[i] == FC_AT_FACE){
      rc = fc_getMeshNumFace(mesh, &numMember);
      fail_unless(rc == FC_SUCCESS, "failed to get num member");
    } else if (assocs[i] == FC_AT_ELEMENT){
      rc = fc_getMeshNumElement(mesh, &numMember);
      fail_unless(rc == FC_SUCCESS, "failed to get num member");
    } else {
      numMember = 1;
    }

    for (j = 0; j < numMember; j+=5){
      fc_addMemberToSubset(subsets[i], j);
    }
  }

  //now write to file and reload
  rc = fc_writeDataset(dataset, filename, FC_FT_EXODUS);
  fail_unless(rc == FC_SUCCESS, "cant write dataset");
  rc = fc_loadDataset(filename, &newDataset);
  fail_unless(rc == FC_SUCCESS, "cant load dataset");
  
  //make sure all the lazy loaded things are not there and
  //that the others are there
  rc = fc_getMeshes(newDataset, &numMeshes, &meshes);
  fail_unless(rc == FC_SUCCESS, "cant get meshes");
  fail_unless(numMeshes == 1, "wrong num meshes");
  rc = fc_getSubsets(meshes[0], &numReturnSubsets, &returnSubsets);
  fail_unless(rc == FC_SUCCESS, "cant get subsets");
  //no at_whole, only face subset on 3D mesh
  fail_unless(numReturnSubsets == numAssoc-2, "wrong num subsets"); 

  for (i = 0; i < numReturnSubsets; i++){
    FC_AssociationType retassoc;
    subSlot = _fc_getSubSlot(returnSubsets[i]);
    fail_unless(subSlot != NULL, "abort: failed to get subslot");
    retassoc = subSlot->assoc;

    if (retassoc == FC_AT_FACE || retassoc == FC_AT_EDGE){
      fail_unless(subSlot->sia.vals == NULL, "members should now be null");

      // load & unload subset members
      rc = fc_getSubsetMembersAsArray(returnSubsets[i], &numMember, &members); // force loading
      fail_unless(rc == FC_SUCCESS, "failed to get members");
      fail_unless(subSlot->sia.vals != NULL, "abort: members should have loaded");
      fail_unless(!memcmp(members, subSlot->sia.vals, numMember*sizeof(int)),
		  "mismatch of members");
      if (retassoc == FC_AT_EDGE){
	rc = fc_getMeshNumEdge(mesh, &compMember);
      } else {
	rc = fc_getMeshNumFace(mesh, &compMember);
      }
      fail_unless(rc == FC_SUCCESS, "failed to get num member");
      for (j = 0; j < numMember; j++){
	fail_unless(members[j] == j*5, "mismatch of members");
      }
      wholemembers = (int)(compMember/5);
      fail_unless(wholemembers == numMember, "wrong number of members");
      free(members);

      if (retassoc == FC_AT_FACE){
	// make copy (uncommitted) subset
	rc = fc_copySubset(returnSubsets[i], meshes[0], "copy_sub", &tempSubset);
	fail_unless(rc == FC_SUCCESS, "failed to copy sub");	
	fail_unless(subSlot->sia.vals != NULL, "abort: members should have loaded");
      }

      rc = fc_releaseSubset(returnSubsets[i]);
      fail_unless(rc == FC_SUCCESS, "failed to release members");
      fail_unless(subSlot->sia.vals == NULL, "members should now be null");

      if (retassoc == FC_AT_FACE){
	//check that the uncommitted subset does nor release
	subSlot = _fc_getSubSlot(tempSubset);
	fail_unless(subSlot != NULL, "abort: failed to get subslot");
	fail_unless(subSlot->sia.vals != NULL, "abort: members should have loaded");

	rc = fc_releaseSubset(tempSubset);
	fail_unless(rc == FC_SUCCESS, "failed to release members");
	fail_unless(subSlot->sia.vals != NULL, "members should not be null");
	fc_deleteSubset(tempSubset);
      }
    } else if (retassoc == FC_AT_VERTEX || retassoc == FC_AT_ELEMENT){
      //not lazily loaded
      fail_unless(subSlot->sia.vals != NULL, "members should not be null");      
      rc = fc_releaseSubset(returnSubsets[i]);
      fail_unless(rc == FC_SUCCESS, "failed to release members");
      fail_unless(subSlot->sia.vals != NULL, "members should not be null");
    } else {
      fail_unless( 0 == 1, "should not be another other subset");
    }
  }
  free(returnSubsets);
  free(meshes);

  // --- special cases
  
  // releasing a null handle does not cause an error
  rc = fc_releaseSubset(FC_NULL_SUBSET);
  fail_unless(rc == FC_SUCCESS, "NULL handle should not fail");
  
  // --- errors
  
  // fc_releaseSubset() --- bad args
  rc = fc_releaseSubset(badSubset);
  fail_unless(rc != FC_SUCCESS, "bad subset should fail");
 
  // cleanup
  fc_deleteDataset(dataset);
  fc_deleteDataset(newDataset);
}
END_TEST 

// test error input for fc_printSubset()
START_TEST(print)
{
  FC_ReturnCode rc;
  int i, j;
  int numAssoc = 5, numEntity;
  FC_AssociationType assocs[5] = { FC_AT_VERTEX, FC_AT_EDGE, FC_AT_FACE,
				   FC_AT_ELEMENT, FC_AT_WHOLE_MESH };
  int numMember = 5, memberIDs[5] = { 1, 3, 5, 7, 9 };
  FC_Dataset dataset;
  FC_Subset subsets[5], badSubset = { 999, 999 };
  FC_Mesh *returnMeshes;
  int numReturnMeshes;
  FC_Mesh mesh, wrongMesh, badMesh = { 999, 999 };
  FC_Variable variables[5], badVariable = { 999, 999 };
  char* subsetNames[5] = { "first 5 odd verts", "first 5 odd edges", 
			   "first 5 odd faces", "first 5 odd elements", 
			   "the whole thing" };
  char* varNames[5] = { "var verts", "var edges", "var faces", "var elements", 
			"var whole" };
  void* data;

  // setup
  fc_loadDataset("../data/gen_multimesh.ex2", &dataset);
  rc = fc_getMeshByName(dataset, "hex mesh", &numReturnMeshes,
			&returnMeshes);
  fail_unless(rc == FC_SUCCESS, "aborted - cant get mesh by name");
  fail_unless(numReturnMeshes == 1,
	      "failed to find unique mesh by name");
  mesh = returnMeshes[0];
  free(returnMeshes);

  fc_getMeshByName(dataset, "pyramid mesh", &numReturnMeshes,
		   &returnMeshes);
  fail_unless(rc == FC_SUCCESS, "aborted - cant get mesh by name");
  fail_unless(numReturnMeshes == 1,
	      "failed to find unique mesh by name");
  wrongMesh = returnMeshes[0];
  free(returnMeshes);

  for (i = 0; i < numAssoc; i++) 
    fc_createSubset(mesh, subsetNames[i], assocs[i], &subsets[i]);
  for (i = 0; i < numAssoc-1; i++) 
    fc_addArrayMembersToSubset(subsets[i], numMember, memberIDs);
  fc_addMemberToSubset(subsets[numAssoc-2], 0);
  for (i = 0; i < numAssoc; i++) {
    fc_createVariable(mesh, varNames[i], &variables[i]);
    fc_getMeshNumEntity(mesh, assocs[i], &numEntity);
    data = calloc(numEntity, sizeof(int)); // all zeros
    fc_setVariableDataPtr(variables[i], numEntity, 1, assocs[i], FC_MT_SCALAR,
			  FC_DT_INT, data);
  }

  // fc_printSubset() -- good args
  //for (i = 0; i < numAssoc; i++) {
  //rc = fc_printSubset(subsets[i], "", 1);
  //fail_unless(rc == FC_SUCCESS, "shouldn't fail");
  //}
  
   // fc_printSubset() -- bad args
  rc = fc_printSubset(badSubset, "Bad Subset!", 1);
  fail_unless(rc != FC_SUCCESS, "should fail if bad subset");
  
  // fc_printSubsetOfMesh() -- good args
  //for (i = 0; i < numAssoc; i++) {
  //rc = fc_printSubsetOfMesh(subsets[i], mesh);
  //fail_unless(rc == FC_SUCCESS, "shouldn't fail");
  //}

  // fc_printSubsetOfMesh() -- bad args
  rc = fc_printSubsetOfMesh(badSubset, mesh);
  fail_unless(rc != FC_SUCCESS, "should fail if bad subset");
  for (i = 0; i < numAssoc; i++) {
    rc = fc_printSubsetOfMesh(subsets[i], badMesh);
    fail_unless(rc != FC_SUCCESS, "should fail if bad mesh");
    if (i == numAssoc - 1) // whole assoc won't fail
      continue;
    rc = fc_printSubsetOfMesh(subsets[i], wrongMesh);
    fail_unless(rc != FC_SUCCESS, "should fail if subset/mesh mismatched");
  }

  // fc_printSubsetOfVariable() -- good args
  //for (i = 0; i < numAssoc; i++) {
  //rc = fc_printSubsetOfVariable(subsets[i], variables[i]);
  //fail_unless(rc != FC_SUCCESS, "shouldn't fail");
  //}

  // fc_printSubsetOfVariable() -- bad args
  for (i = 0; i < numAssoc; i++) {
    rc = fc_printSubsetOfVariable(badSubset, variables[i]);
    fail_unless(rc != FC_SUCCESS, "should fail if bad subset");
    rc = fc_printSubsetOfVariable(subsets[i], badVariable);
    fail_unless(rc != FC_SUCCESS, "should fail if bad variable");
    for (j = 0; j < numAssoc; j++) {
      if (i == j)
	continue;
      rc = fc_printSubsetOfVariable(subsets[i], variables[j]);
      fail_unless(rc != FC_SUCCESS, "should fail if subset/var mismatched");
    }
  }
  
  // cleanup
  fc_deleteDataset(dataset);
}
END_TEST 

// complement
START_TEST(subset_complement)
{
  FC_ReturnCode rc;
  char name[20] = "blue berry", compl_name[20] = "banana nut", *temp_name;
  FC_Dataset dataset, *returnDatasets;
  FC_Mesh mesh, *returnMeshes;
  int numReturnDatasets,numReturnMeshes;
  FC_Subset subset, compl_subset, badSubset = { 999, 999 };
  int numMember = 5, maxNumMember, members[5] = { 0, 1, 5, 7, 10 };
  int compl_numMember = 7, compl_members[7] = { 2, 3, 4, 6, 8, 9, 11 };
  int temp_numMember, temp_maxNumMember, *temp_members;
  FC_AssociationType assoc = FC_AT_VERTEX, temp_assoc;

  // get dataset and mesh
  rc = fc_getDatasetByName(dataset_name,&numReturnDatasets,
			   &returnDatasets);
  fail_unless(rc == FC_SUCCESS, "aborted - cant get dataset by name");
  fail_unless(numReturnDatasets == 1,
	      "failed to find unique dataset by name");
  dataset = returnDatasets[0];
  free(returnDatasets);
  rc = fc_getMeshByName(dataset, "hex mesh", &numReturnMeshes,
			&returnMeshes);
  fail_unless(rc == FC_SUCCESS, "aborted - cant get mesh by name");
  fail_unless(numReturnMeshes == 1,
	      "failed to find unique mesh by name");
  mesh = returnMeshes[0];
  free(returnMeshes);

  fc_getMeshNumVertex(mesh, &maxNumMember);
  fail_unless(maxNumMember == 12, "abort: dataset is not as expected");
    
  // create a subset & add members
  rc = fc_createSubset(mesh, name, assoc, &subset);
  fail_unless(rc == FC_SUCCESS, "failed to create subset");
  rc = fc_addArrayMembersToSubset(subset, numMember, members);
  fail_unless(rc == FC_SUCCESS, "abort: failed to creat start subset");
  
  // create subset complement & check it
  rc = fc_createSubsetComplement(subset, compl_name, &compl_subset);
  fail_unless(rc == FC_SUCCESS, "failed to complement to same mesh");
  rc = fc_getSubsetName(compl_subset, &temp_name);
  fail_unless(rc == FC_SUCCESS, "failed to get name of compl");
  fail_unless(!strcmp(temp_name, compl_name), "mismatch of name");
  free(temp_name);
  rc = fc_getSubsetInfo(compl_subset, &temp_numMember, &temp_maxNumMember,
			&temp_assoc);
  fail_unless(rc == FC_SUCCESS, "failed to get subset info");
  fail_unless(temp_numMember == compl_numMember, "mismatch of numMember");
  fail_unless(temp_maxNumMember == maxNumMember, 
	      "mismatch of maxNumMember"); 
  fail_unless(temp_assoc == assoc, "mismatch of assoc");
  rc = fc_getSubsetMembersAsArray(compl_subset, &temp_numMember, 
				  &temp_members);
  fail_unless(rc == FC_SUCCESS, "failed to get subset info");
  fail_unless(temp_numMember == compl_numMember, "mismatch of numMember");
  fail_unless(!memcmp(temp_members, compl_members, 
		      sizeof(int)*compl_numMember), "mismatch of members");
  free(temp_members);
  
  // ---- check errors
  
  // fc_createSubsetComplement() -- bad args
  rc = fc_createSubsetComplement(badSubset, compl_name, &compl_subset);
  fail_unless(rc != FC_SUCCESS, "should fail with bad subset");
  fail_unless(FC_HANDLE_EQUIV(compl_subset, FC_NULL_SUBSET),
	      "fail should return NULL");
  rc = fc_createSubsetComplement(subset, NULL, &compl_subset);
  fail_unless(rc != FC_SUCCESS, "should fail with NULL name");
  fail_unless(FC_HANDLE_EQUIV(compl_subset, FC_NULL_SUBSET),
	      "fail should return NULL");
  rc = fc_createSubsetComplement(subset, compl_name, NULL);
  fail_unless(rc != FC_SUCCESS, "should fail to copy if NULL handle");
}
END_TEST

// test creating subset intersection
START_TEST(subset_intersection)
{
  FC_ReturnCode rc;
  int j, k, ioper;
  FC_Dataset dataset, *returnDatasets;
  FC_Mesh mesh, *returnMeshes;
  int numReturnDatasets,numReturnMeshes;
  int maxNumMember;   
  char *temp_name;
  char  subsetAName[] = "subset A";
  char  subsetBName[] = "subset B"; 
  char  emptySubsetName[] = "empty subset";
  char  newSubsetName[] = "intersection subset";
  FC_Subset subsetA, subsetB, newSubset, emptySubset, badSubset = { 999, 999 };
  FC_AssociationType temp_assoc, assoc = FC_AT_VERTEX;
  int numSetA = 6, numSetB = 6;
  int membersSetA[] =  { 0,1,2,5,8,10 };
  int membersSetB[] =  { 0,3,5,7,8, 9 };
  int numOps = 3;
  int numVarPerOp[3] = { 4, 4, 3 }; // number of variations of the op
  char ops[3][4][10] = { {"AND", "and", "&&", "&"}, 
			 {"OR", "or", "||", "|"}, 
			 {"XOR", "xor", "^"} };
  int numInters[] = { 3, 9, 6};
  int membersInters[3][12]={ { 0,5,8,-1,-1,-1,-1,-1,-1 },
			     { 0,1,2,3,5,7,8,9,10 },
			     { 1,2,3,7,9,10,-1,-1,-1 } };
  int temp_numMember, temp_maxNumMember, *temp_members; 
  int intersectflag;

  // get dataset and mesh
  rc = fc_getDatasetByName(dataset_name,&numReturnDatasets,
			   &returnDatasets);
  fail_unless(rc == FC_SUCCESS, "aborted - cant get dataset by name");
  fail_unless(numReturnDatasets == 1,
	      "failed to find unique dataset by name");
  dataset = returnDatasets[0];
  free(returnDatasets);
  rc = fc_getMeshByName(dataset, "hex mesh", &numReturnMeshes,
			&returnMeshes);
  fail_unless(rc == FC_SUCCESS, "aborted - cant get mesh by name");
  fail_unless(numReturnMeshes == 1, 
	      "failed to find unique mesh by name");
  mesh = returnMeshes[0];
  free(returnMeshes);

  fc_getMeshNumVertex(mesh, &maxNumMember);
  fail_unless(maxNumMember == 12, "abort: dataset is not as expected");
 
  // create 3 new subsets: A, B, and empty
  rc = fc_createSubset(mesh, subsetAName, assoc, &subsetA);
  fail_unless(rc == FC_SUCCESS, "test aborted, failed to create subset A");
  rc = fc_addArrayMembersToSubset(subsetA, numSetA, membersSetA);
  fail_unless(rc == FC_SUCCESS,"add members to subset A should not fail");
  rc = fc_createSubset(mesh, subsetBName, assoc, &subsetB);
  fail_unless(rc == FC_SUCCESS, "test aborted, failed to create subset B");
  rc = fc_addArrayMembersToSubset(subsetB, numSetB, membersSetB);
  fail_unless(rc == FC_SUCCESS,"add members to subset B should not fail");
  rc = fc_createSubset(mesh, emptySubsetName, assoc, &emptySubset);
  fail_unless(rc == FC_SUCCESS, "abort:, failed to create empty subset");

  for (ioper = 0; ioper < numOps; ioper++) {
    
    // ---- create intersection of A & B
    for (j = 0; j < 2; j++) { // switch order of operands
      for (k = 0; k < numVarPerOp[ioper]; k++) { // try variations of string
	if (j == 0) {
	  if (ioper == 0){ //AND
	    rc = fc_doSubsetsIntersect(subsetA,subsetB,&intersectflag);
	    fail_unless(rc == FC_SUCCESS, "should not fail");		
	  } //AND
	  rc = fc_createSubsetIntersection(subsetA, ops[ioper][k], subsetB,
					   newSubsetName, &newSubset);
	  
	}else {
	  if (ioper == 0){
	    rc = fc_doSubsetsIntersect(subsetB,subsetA,&intersectflag);
	    fail_unless(rc == FC_SUCCESS, "should not fail");		
	    
	  }
	  rc = fc_createSubsetIntersection(subsetB, ops[ioper][k], subsetA,
					   newSubsetName, &newSubset);
	}
	fail_unless(rc == FC_SUCCESS, "should not fail");		
	
	
	// check new subset
	rc = fc_getSubsetName(newSubset, &temp_name);
	fail_unless(rc == FC_SUCCESS,"get name should not fail");  
	fail_unless(!strcmp(temp_name, newSubsetName), "mismatch of name");
	free(temp_name);
	rc = fc_getSubsetInfo(newSubset, &temp_numMember,  
			      &temp_maxNumMember, &temp_assoc);
	fail_unless(rc == FC_SUCCESS,"get info should not fail");
	fail_unless(temp_numMember == numInters[ioper] && 
		    temp_maxNumMember == maxNumMember && 
		    temp_assoc == assoc, "info mismatch");
	rc = fc_getSubsetMembersAsArray(newSubset, &temp_numMember, 
					&temp_members);
	fail_unless(rc == FC_SUCCESS, "failed to get members as array");
	fail_unless(temp_numMember == numInters[ioper], 
		    "mismatch of numMember");
	fail_unless(!memcmp(temp_members, membersInters[ioper],
			    sizeof(int)*numInters[ioper]), "mismatch of members");
	if (ioper == 0){
	  fail_unless(intersectflag == (temp_numMember > 0 ? 1: 0), 
		      "wrong value for intersect flag");
	}
	free(temp_members);
	// delete new subset
	fc_deleteSubset(newSubset);
	fail_unless(rc == FC_SUCCESS, "failed delete complement subset");    
      }
    }
    
    //---- intersect with empty subset (AND = empty, OR,XOR = original set)
    // (don't need to test ops variations again)
    for (j = 0; j < 2; j++) { // switch order of operands
      if (j == 0) {
	rc = fc_createSubsetIntersection(subsetA, ops[ioper][0],
					 emptySubset, newSubsetName, &newSubset);
	
      } else{
	rc = fc_createSubsetIntersection(emptySubset, ops[ioper][0],
					 subsetA, newSubsetName, &newSubset);
      }
      fail_unless(rc == FC_SUCCESS, "should not fail with NULL subset");
      
      
      // check new subset
      rc = fc_getSubsetName(newSubset, &temp_name);
      fail_unless(rc == FC_SUCCESS,"get name should not fail");  
      fail_unless(!strcmp(temp_name, newSubsetName), "mismatch of name");
      free(temp_name);
      rc = fc_getSubsetInfo(newSubset, &temp_numMember,  &temp_maxNumMember, 
			    &temp_assoc);
      fail_unless(rc == FC_SUCCESS,"geting info should not fail");
      if (ioper == 0)  // AND = empty set
	fail_unless(temp_numMember == 0, "mismatch of numMember");
      else  // OR,XOR = original set
	fail_unless(temp_numMember == numSetA, "mismatch of numMember");
      fail_unless(temp_assoc == assoc && temp_maxNumMember == maxNumMember,
		  "info mismatch");
      rc = fc_getSubsetMembersAsArray(newSubset, &temp_numMember, 
				      &temp_members);
      fail_unless(rc == FC_SUCCESS, "failed to get members as array");
      if (ioper == 0) {
	fail_unless(temp_numMember == 0 && temp_members == NULL,
		    "AND should return an empty subset");
      }
      else {
	fail_unless(temp_numMember == numSetA, "mismatch of numMember");
	fail_unless(!memcmp(temp_members, membersSetA, sizeof(int)*numSetA), 
		    "msimatch of members");
	free(temp_members);
      }
      // delete new subset
      fc_deleteSubset(newSubset);
      fail_unless(rc == FC_SUCCESS, "failed delete complement subset");
    }
    
    // ---- intersect with FC_NULL_SUBSET should behave like empty subset
    for (j = 0; j < 2; j++) { // switch order of operands
      if (j == 0) {
	rc = fc_createSubsetIntersection(subsetA, ops[ioper][0], 
					 FC_NULL_SUBSET, newSubsetName, &newSubset);
      } else{
	rc = fc_createSubsetIntersection(FC_NULL_SUBSET, ops[ioper][0], 
					 subsetA, newSubsetName, &newSubset);
      }
      fail_unless(rc == FC_SUCCESS, "should not fail with NULL subset");
      
      
      // check new subset
      rc = fc_getSubsetName(newSubset, &temp_name);
      fail_unless(rc == FC_SUCCESS,"get name should not fail");  
      fail_unless(!strcmp(temp_name, newSubsetName), "mismatch of name");
      free(temp_name);
      rc = fc_getSubsetInfo(newSubset, &temp_numMember,  &temp_maxNumMember, 
			    &temp_assoc);
      fail_unless(rc == FC_SUCCESS,"geting info should not fail");
      if (ioper == 0)  // AND = empty set
	fail_unless(temp_numMember == 0, "mismatch of numMember");
      else  // OR,XOR = original set
	fail_unless(temp_numMember == numSetA, "mismatch of numMember");
      fail_unless(temp_assoc == assoc && temp_maxNumMember == maxNumMember,
		  "info mismatch");
      rc = fc_getSubsetMembersAsArray(newSubset, &temp_numMember, 
				      &temp_members);
      fail_unless(rc == FC_SUCCESS, "failed to get members as array");
      if (ioper == 0) { // AND = empty set
	fail_unless(temp_numMember == 0 && temp_members == NULL,
		    "AND should return an empty subset");
      }
      else { // OR,XOR = original set
	fail_unless(temp_numMember == numSetA, "mismatch of numMember");
	fail_unless(!memcmp(temp_members, membersSetA, sizeof(int)*numSetA), 
		    "msimatch of members");
	free(temp_members);
      }
      // delete new subset
      fc_deleteSubset(newSubset);
      fail_unless(rc == FC_SUCCESS, "failed delete complement subset");
    }
    
    // --- intersect with itself is o.k.
    rc = fc_createSubsetIntersection(subsetA, ops[ioper][0], subsetA, 
				     newSubsetName, &newSubset);
    fail_unless(rc == FC_SUCCESS, "should not fail with same subset");
    
    // check new subset
    rc = fc_getSubsetName(newSubset, &temp_name);
    fail_unless(rc == FC_SUCCESS,"get name should not fail");  
    fail_unless(!strcmp(temp_name, newSubsetName), "mismatch of name");
    free(temp_name);
    rc = fc_getSubsetInfo(newSubset, &temp_numMember,  &temp_maxNumMember, 
			  &temp_assoc);
    fail_unless(rc == FC_SUCCESS,"geting info should not fail");
    if (ioper < 2)  // AND,OR = get back same set
      fail_unless(temp_numMember == numSetA, "mismatch of numMember");
    else  // XOR = empty set
      fail_unless(temp_numMember == 0, "mismatch of numMember");
    fail_unless(temp_assoc == assoc && temp_maxNumMember == maxNumMember,
		"info mismatch");
    rc = fc_getSubsetMembersAsArray(newSubset, &temp_numMember, &temp_members);
    fail_unless(rc == FC_SUCCESS, "failed to get members as array");
    if (ioper < 2) { // AND,OR = get back same set
      fail_unless(temp_numMember == numSetA, "mismatch of numMember");
      fail_unless(!memcmp(temp_members, membersSetA, sizeof(int)*numSetA), 
		  "msimatch of members");
      free(temp_members);
    }
    else { // XOR = empty set
      fail_unless(temp_numMember == 0 && temp_members == NULL,
		  "AND should return an empty subset");
    }
    // delete new subset
    fc_deleteSubset(newSubset);
    fail_unless(rc == FC_SUCCESS, "failed delete complement subset");

    // --- check errors 
    
    // fc_getSubsetMembersAsArray() - can't have two nulls
    rc = fc_createSubsetIntersection(FC_NULL_SUBSET, ops[0][0], 
				     FC_NULL_SUBSET, newSubsetName, &newSubset);
    fail_unless(rc != FC_SUCCESS, "shouldn't do two NULL subsets");
    fail_unless(FC_HANDLE_EQUIV(newSubset, FC_NULL_SUBSET), 
		"fail should return NULL");

    rc = fc_createSubsetIntersection(badSubset, ops[0][0], subsetB,
				     newSubsetName, &newSubset);
    fail_unless(rc != FC_SUCCESS, "shouldn't work with bad first subset");
    fail_unless(FC_HANDLE_EQUIV(newSubset, FC_NULL_SUBSET), 
		"fail should return NULL");

    rc = fc_createSubsetIntersection(subsetA, ops[0][0], badSubset,
				     newSubsetName, &newSubset);
    fail_unless(rc != FC_SUCCESS, "shouldn't work with bad second subset");
    fail_unless(FC_HANDLE_EQUIV(newSubset, FC_NULL_SUBSET), 
		"fail should return NULL");

    rc = fc_createSubsetIntersection(subsetA, NULL, subsetB,
				     newSubsetName, &newSubset);
    fail_unless(rc != FC_SUCCESS, "shouldn't work with NULL operation");
    fail_unless(FC_HANDLE_EQUIV(newSubset, FC_NULL_SUBSET), 
		"fail should return NULL");

    rc = fc_createSubsetIntersection(subsetA, "unknown", subsetB,
				     newSubsetName, &newSubset);
    fail_unless(rc != FC_SUCCESS, "shouldn't work with bad operation");
    fail_unless(FC_HANDLE_EQUIV(newSubset, FC_NULL_SUBSET), 
		"fail should return NULL");

    rc = fc_createSubsetIntersection(subsetA, ops[0][0], subsetB,
				     NULL, &newSubset);
    fail_unless(rc != FC_SUCCESS, "shouldn't work with NULL name");
    fail_unless(FC_HANDLE_EQUIV(newSubset, FC_NULL_SUBSET), 
		"fail should return NULL");

    rc = fc_createSubsetIntersection(subsetA, ops[0][0], subsetB,
				       newSubsetName, NULL);
    fail_unless(rc != FC_SUCCESS, "shouldn't work with NULL subset");
  } //ioper

  //testing for isSubsetSuperset and doSubsetsIntersect
  {
    FC_Subset subsetC, subsetD, subsetE;
    char  subsetCName[] = "subset C";
    char  subsetDName[] = "subset D";
    char  subsetEName[] = "subset E";
    int numSetC = 3, numSetD = 4, numSetE = 4;
    int membersSetC[] =  { 8,9,10};
    int membersSetD[] =  { 0,1,2,3};
    int membersSetE[] =  { 7,8,9,10};
    int flag;
    
    rc = fc_createSubset(mesh, subsetCName, assoc, &subsetC);
    fail_unless(rc == FC_SUCCESS, "test aborted, failed to create subset C");
    rc = fc_addArrayMembersToSubset(subsetC, numSetC, membersSetC);
    fail_unless(rc == FC_SUCCESS,"add members to subset C should not fail");
    rc = fc_createSubset(mesh, subsetDName, assoc, &subsetD);
    fail_unless(rc == FC_SUCCESS, "test aborted, failed to create subset D");
    rc = fc_addArrayMembersToSubset(subsetD, numSetD, membersSetD);
    fail_unless(rc == FC_SUCCESS,"add members to subset D should not fail");
    rc = fc_createSubset(mesh, subsetEName, assoc, &subsetE);
    fail_unless(rc == FC_SUCCESS, "test aborted, failed to create subset E");
    rc = fc_addArrayMembersToSubset(subsetE, numSetE, membersSetE);
    fail_unless(rc == FC_SUCCESS,"add members to subset E should not fail");
    
    rc = fc_isSubsetSuperset(subsetE,subsetC,&flag);
    fail_unless(rc == FC_SUCCESS, "can't determine superset");
    fail_unless(flag == 1, "wrong return val for superset");
    rc = fc_isSubsetSuperset(subsetC,subsetE,&flag);
    fail_unless(rc == FC_SUCCESS, "can't determine superset");
    fail_unless(flag == 0, "wrong return val for superset");
    rc = fc_doSubsetsIntersect(subsetC,subsetE,&intersectflag);
    fail_unless(rc == FC_SUCCESS, "should not fail");		
    fail_unless(intersectflag == 1, "subsets should intersect");
    rc = fc_doSubsetsIntersect(subsetE,subsetC,&intersectflag);
    fail_unless(rc == FC_SUCCESS, "should not fail");		
    fail_unless(intersectflag == 1, "subsets should intersect");

    rc = fc_isSubsetSuperset(subsetC,subsetD,&flag);
    fail_unless(rc == FC_SUCCESS, "can't determine superset");
    fail_unless(flag == 0, "wrong return val for superset");
    rc = fc_isSubsetSuperset(subsetD,subsetC,&flag);
    fail_unless(rc == FC_SUCCESS, "can't determine superset");
    fail_unless(flag == 0, "wrong return val for superset");
    rc = fc_doSubsetsIntersect(subsetD,subsetC,&intersectflag);
    fail_unless(rc == FC_SUCCESS, "should not fail");		
    fail_unless(intersectflag == 0, "subsets should not intersect");
    rc = fc_doSubsetsIntersect(subsetC,subsetD,&intersectflag);
    fail_unless(rc == FC_SUCCESS, "should not fail");		
    fail_unless(intersectflag == 0, "subsets should not intersect");
    
    //with itself
    rc = fc_doSubsetsIntersect(subsetC,subsetC,&intersectflag);
    fail_unless(rc == FC_SUCCESS, "should not fail");		
    fail_unless(intersectflag == 1, "postive intersection with itself");
    
    rc = fc_isSubsetSuperset(subsetC,subsetC,&intersectflag);
    fail_unless(rc == FC_SUCCESS, "should not fail");		
    fail_unless(intersectflag == 1, "postive superset with itself");
      
    // one null
    rc = fc_doSubsetsIntersect(FC_NULL_SUBSET,subsetC,&intersectflag);
    fail_unless(rc == FC_SUCCESS, "should not fail");		
    fail_unless(intersectflag == 0, "zero intersection with empty subset");
      
    rc = fc_isSubsetSuperset(FC_NULL_SUBSET,subsetC,&intersectflag);
    fail_unless(rc == FC_SUCCESS, "should not fail");		
    fail_unless(intersectflag == 0, "no superset with empty subset");
      
    rc = fc_doSubsetsIntersect(subsetC,FC_NULL_SUBSET,&intersectflag);
    fail_unless(rc == FC_SUCCESS, "should not fail");		
    fail_unless(intersectflag == 0, "zero intersection with empty subset");
      
    rc = fc_isSubsetSuperset(subsetC,FC_NULL_SUBSET,&intersectflag);
    fail_unless(rc == FC_SUCCESS, "should not fail");		
    fail_unless(intersectflag == 0, "no superset with empty subset");
      
    //one empty
    rc = fc_doSubsetsIntersect(emptySubset,subsetC,&intersectflag);
    fail_unless(rc == FC_SUCCESS, "should not fail");		
    fail_unless(intersectflag == 0, "zero intersection with empty subset");
    
    rc = fc_isSubsetSuperset(emptySubset,subsetC,&intersectflag);
    fail_unless(rc == FC_SUCCESS, "should not fail");		
    fail_unless(intersectflag == 0, "no superset with empty subset");
    
    rc = fc_doSubsetsIntersect(subsetC,emptySubset,&intersectflag);
    fail_unless(rc == FC_SUCCESS, "should not fail");		
    fail_unless(intersectflag == 0, "zero intersection with empty subset");
    
    rc = fc_isSubsetSuperset(subsetC,emptySubset,&intersectflag);
    fail_unless(rc == FC_SUCCESS, "should not fail");		
    fail_unless(intersectflag == 0, "no superset with empty subset");

    // -- bad args -- //
    rc = fc_doSubsetsIntersect(FC_NULL_SUBSET, FC_NULL_SUBSET, &intersectflag);
    fail_unless(rc != FC_SUCCESS, "shouldn't do two NULL subsets");
    fail_unless(intersectflag == -1, "fail should return -1");
    rc = fc_isSubsetSuperset(FC_NULL_SUBSET, FC_NULL_SUBSET, &intersectflag);
    fail_unless(rc != FC_SUCCESS, "shouldn't do two NULL subsets");
    fail_unless(intersectflag == -1, "fail should return -1");

    rc = fc_doSubsetsIntersect(badSubset, subsetC, &intersectflag);
    fail_unless(rc != FC_SUCCESS, "shouldn't work with bad first subset");
    fail_unless(intersectflag == -1, "fail should return -1");
    rc = fc_isSubsetSuperset(badSubset, subsetC, &intersectflag);
    fail_unless(rc != FC_SUCCESS, "shouldn't work with bad first subset");
    fail_unless(intersectflag == -1, "fail should return -1");

    rc = fc_doSubsetsIntersect(subsetC, badSubset, &intersectflag);
    fail_unless(rc != FC_SUCCESS, "shouldn't work with bad second subset");
    fail_unless(intersectflag == -1, "fail should return -1");
    rc = fc_isSubsetSuperset(subsetC, badSubset, &intersectflag);
    fail_unless(rc != FC_SUCCESS, "shouldn't work with bad second subset");
    fail_unless(intersectflag == -1, "fail should return -1");

    rc = fc_doSubsetsIntersect(subsetC, subsetC, NULL);
    fail_unless(rc != FC_SUCCESS, "shouldn't work with NULL return val");
    rc = fc_isSubsetSuperset(subsetC, subsetD, NULL);
    fail_unless(rc != FC_SUCCESS, "shouldn't work with NULL return val");

    rc = fc_deleteSubset(subsetC);
    fail_unless(rc == FC_SUCCESS, "failed to delete subset C");
    rc = fc_deleteSubset(subsetD);
    fail_unless(rc == FC_SUCCESS, "failed to delete subset D");
    rc = fc_deleteSubset(subsetE);
    fail_unless(rc == FC_SUCCESS, "failed to delete subset E");
  }
    
  
  // ---- cleanup
  rc = fc_deleteSubset(subsetA);
  fail_unless(rc == FC_SUCCESS, "failed to delete subset A");
  rc = fc_deleteSubset(subsetB);
  fail_unless(rc == FC_SUCCESS, "failed to delete subset B");
  rc = fc_deleteSubset(emptySubset);
  fail_unless(rc == FC_SUCCESS, "failed to delete empty subset");
}
END_TEST


// *********************************************
// ***** Populate the Suite with the tests
// *********************************************
Suite *subset_suite(void)
{
  Suite *suite = suite_create("Subset");

  TCase *tc_private_sub = tcase_create(" - Private Subset ");
  TCase *tc_subset = tcase_create(" - Subset Interface ");
  TCase *tc_subsetOps = tcase_create(" - Subset Operations ");

  // private subset
  suite_add_tcase(suite, tc_private_sub);
  tcase_add_checked_fixture(tc_private_sub, subset_setup, subset_teardown);
  tcase_add_test(tc_private_sub, slot_new_delete);

  // subset interface
  suite_add_tcase(suite, tc_subset);
  tcase_add_checked_fixture(tc_subset, subset_setup, subset_teardown);
  tcase_add_test(tc_subset, create_get_delete);
  tcase_add_test(tc_subset, metadata_query);
  tcase_add_test(tc_subset, one_member_add_del);
  tcase_add_test(tc_subset, get_members);
  tcase_add_test(tc_subset, query_members);
  tcase_add_test(tc_subset, multi_members_add);
  tcase_add_test(tc_subset, multi_members_del);
  tcase_add_test(tc_subset, copy_test);
  tcase_add_test(tc_subset, change_assoc);
  tcase_add_test(tc_subset, release_subset);
  tcase_add_test(tc_subset, print);

  // subset operations
  suite_add_tcase(suite, tc_subsetOps);
  tcase_add_checked_fixture(tc_subsetOps, subset_setup, subset_teardown);
  tcase_add_test(tc_subsetOps, subset_complement);
  tcase_add_test(tc_subsetOps, subset_intersection);

  return suite;
}
