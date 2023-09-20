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
 * \file checkelemdeath.c
 * \brief Unit testing of \ref ElemDeath module
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/unittest/checkelemdeath.c,v $
 * $Revision: 1.44 $
 * $Date: 2006/10/19 03:14:52 $
 *
 */

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <check.h>
#include "fc.h"
#include "fcP.h"
#include "checkall.h"

// **** test fixtures

static void elemdeath_setup(void) {
  FC_ReturnCode rc;
  if (isForking) {
    fc_setLibraryVerbosity(fc_messages);
    rc = fc_initLibrary();
    fail_unless(rc == FC_SUCCESS, "test aborted: failed to init library");
  }
}

static void elemdeath_teardown(void) {
  FC_ReturnCode rc;
  if (isForking) {
    rc = fc_finalLibrary();
    fail_unless(rc == FC_SUCCESS, "failed to final library at end of test");
  }
}

// basic layout of vertices is on 2 X 2 X 2 block of hexes
// for 1D use only lower edge, for 2D use bottom quads
START_TEST(exposed_skin)
{
  FC_ReturnCode rc;
  int i, j, k;
  char names[4][20] = { "point mesh", "line mesh", "quad mesh", "hex mesh" };
  FC_Dataset dataset;
  FC_Mesh meshes[4];
  FC_Subset badSubset = { 999, 999 };
  int numDim = 3;
  int numElemType = 4;
  FC_ElementType elemTypes[4] = { FC_ET_POINT,  FC_ET_LINE,
                                  FC_ET_QUAD,   FC_ET_HEX };
  int numVertPerType[4] = { 3, 3, 9, 27}, numElemPerType[4] = { 3, 2, 4, 8 };
  int conns[4][27*8] = { { 0,    1,    2 },
			 { 0, 1,    1, 2 },
			 { 0, 1, 4, 3,    1, 2, 5, 4,
			   3, 4, 7, 6,    4, 5, 8, 7 },
			 { } };
  double* coords;
  FC_Coords lowers = { 0., 0., 0. };
  FC_Coords uppers = { 1., 1., 1. };
  FC_AssociationType assocs[4] = { FC_AT_UNKNOWN, FC_AT_VERTEX, FC_AT_EDGE,
				   FC_AT_FACE };
  // The answer if remove 2nd element
  int numMembers[4] = { 0, 1, 2, 3 };
  int numMemberVerts[4] = { 0, 1, 3, 7 };
  int memberVerts[4][9] = { { }, { 1 }, { 1, 4, 5 }, 
			    { 1, 4, 5, 10, 11, 13, 14 } }; 

  // setup - create test dataset
  rc = fc_createDataset("temp.xxx", &dataset);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create dataset");

  // Make the hex mesh, then reuse coords for other meshes  
  rc = fc_createSimpleHexMesh(dataset, names[3], 2, 2, 2, lowers, uppers,
			      &meshes[3]);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create simple hex mesh");
  fc_getMeshCoordsPtr(meshes[3], &coords);
  for (i = 0; i < numElemType-1; i++) {
    rc = fc_createMesh(dataset, names[i], &meshes[i]);
    fail_unless(rc == FC_SUCCESS, "abort: failed to created mesh");
    rc = fc_setMeshCoords(meshes[i], numDim, numVertPerType[i], coords);
    fail_unless(rc == FC_SUCCESS, "failed to set vertex coords");
    rc = fc_setMeshElementConns(meshes[i], elemTypes[i], numElemPerType[i], 
				conns[i]);
    fail_unless(rc == FC_SUCCESS, "failed to set element conns");
  }

  // test meshes of different types
  for (i = 0; i < numElemType; i++) { // loop over elem types
    FC_Mesh mesh;
    FC_Subset meshSkin, dead, exposed, emptySubset;
    FC_AssociationType temp_assoc;
    int temp_numMember, *temp_memberIDs, temp_numVertex, *temp_vertexIDs;
    
    mesh = meshes[i];

    // --- create data structures for testing: skin, displVar, dead elements
    
    // create skin
    rc = fc_getMeshSkin(mesh, &meshSkin);
    if (i > 0)
      fail_unless(rc == FC_SUCCESS, "failed to get mesh Skin");
    // Create dead element subset -> the 2nd element
    rc = fc_createSubset(mesh, "dead", FC_AT_ELEMENT, &dead);
    fail_unless(rc == FC_SUCCESS, "abort: failed to create dead elem subset"); 
    rc = fc_addMemberToSubset(dead, 1);
    fail_unless(rc == FC_SUCCESS, "abort: failed to add element to subset");
    // Create empty subset
    rc = fc_createSubset(mesh, "empty", FC_AT_ELEMENT, &emptySubset);
    fail_unless(rc == FC_SUCCESS, "abort: failed to create empty subset"); 

    // --- test w/o & w/ pasing mesh skin
    for (j = 0; j < 2; j++) {
      if (j == 0)
	rc = fc_getExposedSkin(dead, NULL, &exposed);
      else
	rc = fc_getExposedSkin(dead, &meshSkin, &exposed);
      if (i == 0) {
	fail_unless(rc != FC_SUCCESS, "should fail for point mesh - no skin");
	fail_unless(FC_HANDLE_EQUIV(exposed, FC_NULL_SUBSET),
		    "fail should return null");
	continue;
      }
      fail_unless(rc == FC_SUCCESS, "failed to get exposed skin");
      fc_getSubsetAssociationType(exposed, &temp_assoc);
      fail_unless(temp_assoc == assocs[i], "mismatch of assoc");
      fc_getSubsetMembersAsArray(exposed, &temp_numMember, &temp_memberIDs);
      fail_unless(temp_numMember = numMembers[i], "mismatch of numMember");
      // since we don't know how edges & faces numbered, check the verts
      fc_changeMeshEntityType(mesh, temp_assoc, temp_numMember, temp_memberIDs,
			      FC_AT_VERTEX, 1, &temp_numVertex, &temp_vertexIDs);
      free(temp_memberIDs);
      fail_unless(temp_numVertex == numMemberVerts[i], "mismatch of numVertex");
      for (k = 0; k < temp_numVertex; k++)
	fail_unless(temp_vertexIDs[k] == memberVerts[i][k], 
		    "mismatch of vert IDs");
      free(temp_vertexIDs);
    }

    // --- Test special cases

    // Error if empty subset
    rc = fc_getExposedSkin(emptySubset, NULL, &exposed);
    fail_unless(rc != FC_SUCCESS, "Should fail for empty subset");
    fail_unless(FC_HANDLE_EQUIV(exposed, FC_NULL_SUBSET),
		"fail should return null");

    // --- Test bad args
    rc = fc_getExposedSkin(badSubset, NULL, &exposed);
    fail_unless(rc != FC_SUCCESS, "should fail if bad subset");
    fail_unless(FC_HANDLE_EQUIV(exposed, FC_NULL_SUBSET),
		"fail should return null");
    rc = fc_getExposedSkin(dead, NULL, NULL);
    fail_unless(rc != FC_SUCCESS, "should fail if null return arg");


    if (elemTypes[i] == FC_ET_HEX){
      FC_Subset wrongassoc;
      FC_Subset wrongmeshskin;
      FC_Subset whole, allelems;
      int numID,*IDs;

      rc = fc_createSubset(mesh,"wrongassoc",FC_AT_FACE,&wrongassoc);
      fail_unless(rc == FC_SUCCESS, "cant make new assoc to test");
      rc = fc_addMemberToSubset(wrongassoc,1);
      fail_unless(rc == FC_SUCCESS, "cant add member to subset");
      rc = fc_getExposedSkin(wrongassoc,&meshSkin, &exposed);
      fail_unless(rc != FC_SUCCESS, "should fail if subset has wrong assoc");
      fail_unless(FC_HANDLE_EQUIV(exposed, FC_NULL_SUBSET),
		  "fail should return null");

      rc = fc_getMeshSkin(meshes[i-1], &wrongmeshskin);
      fail_unless(rc == FC_SUCCESS, "cant get mesh skin for testing");
      rc = fc_getExposedSkin(dead,&wrongmeshskin,&exposed);
      fail_unless(rc != FC_SUCCESS, "should fail if subset/skin is on different meshes");
      fail_unless(FC_HANDLE_EQUIV(exposed, FC_NULL_SUBSET),
		  "fail should return null");

      //special case - all elements are dead
      fc_createSubset(mesh, "whole", FC_AT_WHOLE_MESH, &whole);
      fc_addMemberToSubset(whole, 0);
      fc_createSubset(mesh, "all elements", FC_AT_ELEMENT, &allelems);
      fc_getMeshNumElement(mesh,&numID);
      IDs = malloc(numID*sizeof(int));
      for (j = 0; j < numID; j++)
	IDs[j] = j;
      fc_addArrayMembersToSubset(allelems, numID, IDs);
      free(IDs);

      rc = fc_getExposedSkin(whole,NULL,&exposed);
      fail_unless(rc == FC_SUCCESS, "should work for full dead subset");
      fail_unless(FC_HANDLE_EQUIV(exposed, FC_NULL_SUBSET),
		  "full dead subset should return NULL");
      rc = fc_getExposedSkin(allelems,NULL,&exposed);
      fail_unless(rc == FC_SUCCESS, "should work for full dead subset");
      fail_unless(FC_HANDLE_EQUIV(exposed, FC_NULL_SUBSET),
		  "full dead subset should return NULL");

      fc_deleteSubset(wrongmeshskin);
      fc_deleteSubset(whole);
      fc_deleteSubset(allelems);
    }
  }

  // --- cleanup
  rc = fc_deleteDataset(dataset); 
  fail_unless(rc == FC_SUCCESS, "failed to delete dataset at end of tests");
}
END_TEST

// basic layout of vertices is on 2 X 2 X 2 block of hexes
// for 1D use only lower edge, for 2D use bottom quads
// Testing a wrapper so only testing that it returns same thing as
// calling underyling functions, no need to test bad args or special cases
// since wrapper just passes args through
START_TEST(tear_length)
{
  FC_ReturnCode rc, rc2, rc3;
  int i, j;
  char names[4][20] = { "point mesh", "line mesh", "quad mesh", "hex mesh" };
  FC_Dataset dataset;
  FC_Mesh meshes[4];
  FC_Variable variables[4];
  int numDim = 3;
  int numElemType = 4;
  FC_ElementType elemTypes[4] = { FC_ET_POINT,  FC_ET_LINE,
                                  FC_ET_QUAD,   FC_ET_HEX };
  int numVertPerType[4] = { 3, 3, 9, 27}, numElemPerType[4] = { 3, 2, 4, 8 };
  int conns[4][27*8] = { { 0,    1,    2 },
			 { 0, 1,    1, 2 },
			 { 0, 1, 4, 3,    1, 2, 5, 4,
			   3, 4, 7, 6,    4, 5, 8, 7 },
			 { } };
  double* coords;
  FC_Coords lowers = { 0., 0., 0. };
  FC_Coords uppers = { 1., 1., 1. };
  double temp_length, length;

  // setup - create test dataset
  rc = fc_createDataset("temp.xxx", &dataset);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create dataset");

  // Make the hex mesh, then reuse coords for other meshes  
  rc = fc_createSimpleHexMesh(dataset, names[3], 2, 2, 2, lowers, uppers,
			      &meshes[3]);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create simple hex mesh");
  fc_getMeshCoordsPtr(meshes[3], &coords);
  for (i = 0; i < numElemType-1; i++) {
    rc = fc_createMesh(dataset, names[i], &meshes[i]);
    fail_unless(rc == FC_SUCCESS, "abort: failed to created mesh");
    rc = fc_setMeshCoords(meshes[i], numDim, numVertPerType[i], coords);
    fail_unless(rc == FC_SUCCESS, "failed to set vertex coords");
    rc = fc_setMeshElementConns(meshes[i], elemTypes[i], numElemPerType[i], 
				conns[i]);
    fail_unless(rc == FC_SUCCESS, "failed to set element conns");
  }
  for (i = 0; i < numElemType; i++) {
    rc = fc_createVariable(meshes[i], "displacement", &variables[i]);
    rc = fc_setVariableData(variables[i], numVertPerType[i], numDim, 
			    FC_AT_VERTEX, FC_MT_VECTOR, FC_DT_DOUBLE, coords);
    fail_unless(rc == FC_SUCCESS, "failed to set displ var data");
  }

  // test meshes of different types
  for (i = 0; i < numElemType; i++) { // loop over elem types
    FC_Mesh mesh;
    FC_Subset meshSkin, dead, exposed, emptySubset;
    FC_Variable displ;
    
    mesh = meshes[i];
    displ = variables[i];

    // --- create data structures for testing: skin, displVar, dead elements
    
    // create skin
    rc = fc_getMeshSkin(mesh, &meshSkin);
    if (i > 0)
      fail_unless(rc == FC_SUCCESS, "failed to get mesh Skin");
    // Create dead element subset -> the 2nd element
    rc = fc_createSubset(mesh, "dead", FC_AT_ELEMENT, &dead);
    fail_unless(rc == FC_SUCCESS, "abort: failed to create dead elem subset"); 
    rc = fc_addMemberToSubset(dead, 1);
    fail_unless(rc == FC_SUCCESS, "abort: failed to add element to subset");
    // Create empty subset
    rc = fc_createSubset(mesh, "empty", FC_AT_ELEMENT, &emptySubset);
    fail_unless(rc == FC_SUCCESS, "abort: failed to create empty subset"); 

    // --- test w/o & w/ pasing mesh skin, AND w & w/o displ
    for (j = 0; j < 4; j++) {
      if (j == 0) {
	rc = fc_calcTearLength(dead, NULL, NULL, &temp_length);
	rc2 = fc_getExposedSkin(dead, NULL, &exposed);
	rc3 = fc_getSubsetDiameter(exposed, &length, NULL, NULL);
      }
      else if (j == 1) {
	rc = fc_calcTearLength(dead, &meshSkin, NULL, &temp_length);
	rc2 = fc_getExposedSkin(dead, &meshSkin, &exposed);
	rc3 = fc_getSubsetDiameter(exposed, &length, NULL, NULL);
      }
      else if (j == 2) {
	rc = fc_calcTearLength(dead, NULL, &displ, &temp_length);
	rc2 = fc_getExposedSkin(dead, NULL, &exposed);
	rc3 = fc_getDisplacedSubsetDiameter(exposed, displ, &length, NULL,
					   NULL);
      }
      else {
	rc = fc_calcTearLength(dead, &meshSkin, &displ, &temp_length);
	rc2 = fc_getExposedSkin(dead, &meshSkin, &exposed);
	rc3 = fc_getDisplacedSubsetDiameter(exposed, displ, &length, NULL,
					    NULL);
      }
      if (i == 0) {
	fail_unless(rc != FC_SUCCESS, "should fail for point mesh - no skin");
	fail_unless(rc2 != FC_SUCCESS && rc3 != FC_SUCCESS, "abort: shoudl fail");
	fail_unless(temp_length == length, "fail should return null");
	continue;
      }
      fail_unless(rc == FC_SUCCESS, "failed to get tear length");
      fail_unless(rc2 == FC_SUCCESS && rc3 == FC_SUCCESS, "abort: alt method failed");
      fail_unless(temp_length == length, "should get same length either way");
    }

    // --- Test special cases
    // No need to test special cases of the underlying calls (fc_getExposedSkin
    // & fc_getSubsetDiameter) since this is just a wrapper

    // --- Test bad args
    // No need to test bad args since this is just a wrapper (args are
    // tested by wrapped functions)
  }

  // --- cleanup
  rc = fc_deleteDataset(dataset); 
  fail_unless(rc == FC_SUCCESS, "failed to delete dataset at end of tests");
}
END_TEST


//methods related to breaking a mesh
//this piece tests fc_subsetSegmentsMesh &
// fc_subsetPlusneighborGreaterSegementationMesh
// (removed fc_subsetSegmentsMesh_Time as a faunction, but some of the
// setup here is left from when that existed)
// NOTE still have to finish neighbors for hex - need a case where
// it really does result in greater segmentation....
START_TEST(breaks)
{
  FC_ReturnCode rc;
  int i,j,jj,k;
  char names[4][20] = { "point mesh", "line mesh", "quad mesh", "hex mesh" };
  FC_Dataset dataset;
  FC_Mesh meshes[4];

  int numDim = 3;
  int numElemType = 4;
  FC_ElementType elemTypes[4] = { FC_ET_POINT,  FC_ET_LINE,
				  FC_ET_QUAD,   FC_ET_HEX };

  int numVertPerType[4] = { 4, 4, 16, 64}, numElemPerType[4] = { 4, 3, 9, 27 };
  int conns[4][64*27] = { { 0,       1,       2,      3},
			  { 0, 1,    1, 2,    2, 3},
			  { 0, 1, 5, 4,    1, 2, 6, 5,    2, 3, 7, 6,
			    4, 5, 9, 8,    5, 6, 10, 9,   6, 7, 11, 10,
			    8, 9, 13, 12,  9, 10, 14,13,  10, 11,15, 14
			  },
			  { } };
  double* coords;
  FC_Coords lowers = { 0., 0., 0. };
  FC_Coords uppers = { 1., 1., 1. };

  int shared_dim = 0;  //not checking shared_dim (segmentation)
  //since it is just passed into fc_segment and that is checked there
  int shared_neighbordim = 0; //decide if want to iterate thru this...

  // setup - create test dataset
  rc = fc_createDataset("temp.xxx", &dataset);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create dataset");

  // Make the hex mesh, then reuse coords for other meshes  
  rc = fc_createSimpleHexMesh(dataset, names[3], 3, 3, 3, lowers, uppers,
			      &meshes[3]);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create simple hex mesh");
  fc_getMeshCoordsPtr(meshes[3], &coords);
  for (i = 0; i < numElemType-1; i++) {
    rc = fc_createMesh(dataset, names[i], &meshes[i]);
    fail_unless(rc == FC_SUCCESS, "abort: failed to created mesh");
    rc = fc_setMeshCoords(meshes[i], numDim, numVertPerType[i], coords);
    fail_unless(rc == FC_SUCCESS, "failed to set vertex coords");
    rc = fc_setMeshElementConns(meshes[i], elemTypes[i], numElemPerType[i], 
				conns[i]);
    fail_unless(rc == FC_SUCCESS, "failed to set element conns");
  }


  // test meshes of different types
  for (i = 0; i < numElemType; i++) { // loop over elem types
    //for the seqvar
    FC_Mesh mesh;
    FC_Subset whole;
    int numMem;

    int numStep = 3;
    int numDataPoint;
    int times[3] = {0,1,2};
    FC_Variable* seqvar;
    int checkSeqStep;
    FC_Sequence seq;

    mesh = meshes[i];
    //make a subset of the whole mesh
    rc = fc_getMeshNumElement(mesh, &numMem);
    fail_unless(rc == FC_SUCCESS, "cant get mesh num elem");
    rc = fc_createSubset(mesh, "whole", FC_AT_ELEMENT, &whole);
    fail_unless(rc == FC_SUCCESS, "cant make mesh subset");
    for (k = 0; k < numMem; k++){
      rc = fc_addMemberToSubset(whole,k);
      fc_exitIfErrorPrintf(rc, "cant add member to subset");
    }

    //now make seqvar with dead
    rc = fc_createSequence(dataset,"new_seq",&seq);
    fail_unless(rc == FC_SUCCESS, "failed to create sequence");
    rc = fc_setSequenceCoords(seq,numStep,FC_DT_INT,times);
    fail_unless(rc == FC_SUCCESS, "failed to set sequence coords");
    rc = fc_createSeqVariable(mesh,seq,"deadvar",&checkSeqStep,
		       &seqvar);
    fail_unless(rc == FC_SUCCESS, "failed to create seqVar");

    //now set dead data
    numDataPoint = numElemPerType[i];
    for (k = 0; k < numStep; k++){
      int *data = (int*)malloc(numDataPoint*sizeof(int));
      switch(k){
      case 0: //empty set
	for (j = 0; j < numDataPoint; j++)
	  data[j] = 0;
	break;
      case 1: //complete set
	for (j = 0; j < numDataPoint; j++)
	  data[j] = 1;
	break;
      case 2: //specified
	switch (elemTypes[i]){
	case FC_ET_QUAD:
	  for (j = 0; j < numDataPoint; j++)
	    data[j] = 0;
	  data[1] = 1;
	  data[4] = 1;
	  data[5] = 1;
	  break;
	case FC_ET_HEX:
	  for (j = 0; j < numDataPoint; j++)
	    data[j] = 0;
	  for (j = 9; j <18; j++)
	    data[j] = 1;
	  break;
	default: 
	  for (j = 0; j < numDataPoint; j++)
	    data[j] = 0;
	  data[numDataPoint/2] = 1;
	  break;
	}
	break;
      default: //empty - shouldnt occur
	for (j = 0; j < numDataPoint; j++)
	  data[j] = 0;
	break;
      }

      rc = fc_setVariableDataPtr(seqvar[k],numDataPoint,1,FC_AT_ELEMENT,
				 FC_MT_SCALAR, FC_DT_INT,data);
      fail_unless(rc == FC_SUCCESS, "failed to set var data");
    }

    //test return values here - these are known quantitities
    for (j = 0; j < numStep; j++){
      switch(elemTypes[i]){
      case FC_ET_POINT:
	switch (j){
	case 0: //empty
	  {
	    FC_Subset threshsubs; 
	    int numSubsets, numSubsets2, numNbr, *nbrIDs;

	    rc = fc_createThresholdSubset(seqvar[j],"=",1,"temp_thresh",
					  &threshsubs);
	    fail_unless(rc == FC_SUCCESS, "failed to create threshold subset");
	    rc = fc_subsetSegmentsMesh(threshsubs,shared_dim,NULL,&numSubsets,
				       NULL);
	    fail_unless(rc == FC_SUCCESS, 
			"should work for empty input subset and optional args");
	    fail_unless(numSubsets == 4,
			"empty set should return same seg as full mesh");
	    rc = fc_subsetSegmentsSubset(threshsubs,whole,shared_dim,NULL,
					 &numSubsets2,NULL);
	    fail_unless(rc == FC_SUCCESS, 
			"should work for empty input subset and optional args");
	    fail_unless(numSubsets2 == numSubsets,
			"empty set should return same seg as full mesh");
	    //could check to see that they are the same subsets here

	    rc = fc_subsetPlusNeighborGreaterSegmentationMesh(threshsubs,
							      shared_neighbordim,
							      shared_dim,&numNbr,
							      &nbrIDs);
	    fail_unless(rc != FC_SUCCESS, "should not work for empty input subset");

	    fc_deleteSubset(threshsubs);
	  }
	  break;
	case 1: //full
	  {
	    FC_Subset threshsubs; 
	    FC_Subset *newSubsets,*newSubsets2;
	    int numSubsets, numSubsets2,numNbr, *nbrIDs;
	    rc = fc_createThresholdSubset(seqvar[j],"=",1,"temp_thresh",
					  &threshsubs);
	    fail_unless(rc == FC_SUCCESS, "failed to create threshold subset");
	    rc = fc_subsetSegmentsMesh(threshsubs,shared_dim,"live",&numSubsets,
				       &newSubsets);
	    fail_unless(rc == FC_SUCCESS, "should work for full input subset");
	    fail_unless(numSubsets == 0, "should return 0 for full set");
	    rc = fc_subsetSegmentsSubset(threshsubs,whole,shared_dim,"live",
					 &numSubsets2,&newSubsets2);
	    fail_unless(rc == FC_SUCCESS, "should work for full input subset");

	    fail_unless(numSubsets2 == numSubsets, "should return 0 for full set");
	    //clean up
	    for (jj = 0; jj < numSubsets; jj++){
	      fc_deleteSubset(newSubsets[jj]);
	      fc_deleteSubset(newSubsets2[jj]);
	    }
	    free(newSubsets);
	    free(newSubsets2);

	    rc = fc_subsetPlusNeighborGreaterSegmentationMesh(threshsubs,
							      shared_neighbordim,
						      shared_dim,&numNbr,&nbrIDs);
	    fail_unless(rc == FC_SUCCESS, "should work for full set");
	    fail_unless(numNbr == 0, "should return no neighbors");
	    fail_unless(nbrIDs == NULL, "should return null neighbor set");

	    fc_deleteSubset(threshsubs);
	  }
	  break;
	case 2: //specified
	  {
	    int ans[3] = {0,1,3};
	    FC_Subset threshsubs; 
	    FC_Subset *newSubsets, *newSubsets2;
	    int numSubsets, numSubsets2, numNbr, *nbrIDs;

	    rc = fc_createThresholdSubset(seqvar[j],"=",1,"temp_thresh",
					  &threshsubs);
	    fail_unless(rc == FC_SUCCESS, "failed to create threshold subset");
	    rc = fc_subsetSegmentsMesh(threshsubs,shared_dim,"live",
				       &numSubsets,&newSubsets);
	    fail_unless(rc == FC_SUCCESS,"failed to subsetSegmentsMesh");
	    fail_unless(numSubsets == 3,
			"failed to get correct number of subsets");
	    rc = fc_subsetSegmentsSubset(threshsubs,whole,shared_dim,"live",
					 &numSubsets2,&newSubsets2);
	    fail_unless(rc == FC_SUCCESS,"failed to subsetSegmentsMesh");
	    fail_unless(numSubsets2 == numSubsets,
			"failed to get correct number of subsets");
	    for (j = 0; j < numSubsets; j++){
	      int temp_numMember;
	      int *temp_memberIDs;

	      fc_getSubsetMembersAsArray(newSubsets[j], &temp_numMember,
					 &temp_memberIDs);
	      fail_unless(rc ==  FC_SUCCESS,
			  "unable to get subset members as array");
	      fail_unless(temp_numMember ==1 ,
			  "wrong number of members in the subset");
	      fail_unless(temp_memberIDs[0] == ans[j], "wrong id in this subset");
	      free (temp_memberIDs);

	      fc_getSubsetMembersAsArray(newSubsets2[j], &temp_numMember,
					 &temp_memberIDs);
	      fail_unless(rc ==  FC_SUCCESS,
			  "unable to get subset members as array");
	      fail_unless(temp_numMember ==1 ,
			  "wrong number of members in the subset");
	      fail_unless(temp_memberIDs[0] == ans[j], "wrong id in this subset");
	      free (temp_memberIDs);
	    }
	    //clean up
	    for (jj = 0; jj < numSubsets; jj++){
	      fc_deleteSubset(newSubsets[jj]);
	      fc_deleteSubset(newSubsets2[jj]);
	    }
	    free(newSubsets);
	    free(newSubsets2);

	    rc = fc_subsetPlusNeighborGreaterSegmentationMesh(threshsubs,
							      shared_neighbordim,
							      shared_dim,
							      &numNbr,&nbrIDs);
	    fail_unless(rc == FC_SUCCESS, "should work");
	    fail_unless(numNbr == 0, "should return no neighbors");
	    fail_unless(nbrIDs == NULL, "should return null neighbor set");

	    fc_deleteSubset(threshsubs);
	  }
	  break;
	default:
	  break;
	}
	break;
      case FC_ET_LINE:
	switch (j){
	case 0:
	  {
	    FC_Subset threshsubs; 
	    FC_Subset *newSubsets, *newSubsets2;
	    int numSubsets, numSubsets2, numNbr, *nbrIDs;
	    rc = fc_createThresholdSubset(seqvar[j],"=",1,"temp_thresh", &threshsubs);
	    fail_unless(rc == FC_SUCCESS, "failed to create threshold subset");
	    rc = fc_subsetSegmentsMesh(threshsubs,shared_dim,"live",&numSubsets,&newSubsets);
	    fail_unless(rc == FC_SUCCESS, "should work for empty input subset");
	    fail_unless(numSubsets == 1, "empty set should return same seg as full mesh");
	    rc = fc_subsetSegmentsSubset(threshsubs,whole,shared_dim,"live",
					 &numSubsets2,&newSubsets2);
	    fail_unless(rc == FC_SUCCESS, "should work for empty input subset");
	    fail_unless(numSubsets2 == numSubsets,
			"empty set should return same seg as full mesh");
	    //could check to see that they are the same subsets here
	    //clean up
	    for (jj = 0; jj < numSubsets; jj++){
	      fc_deleteSubset(newSubsets[jj]);
	      fc_deleteSubset(newSubsets2[jj]);
	    }
	    free(newSubsets);
	    free(newSubsets2);
	    rc = fc_subsetPlusNeighborGreaterSegmentationMesh(threshsubs,shared_neighbordim,
							      shared_dim,&numNbr,&nbrIDs);
	    fail_unless(rc != FC_SUCCESS, "should not work for empty input subset");

	    fc_deleteSubset(threshsubs);
	  }
	  break;
	case 1:
	  {
	    FC_Subset threshsubs; 
	    int numSubsets, numSubsets2, numNbr, *nbrIDs;
	    rc = fc_createThresholdSubset(seqvar[j],"=",1,"temp_thresh", &threshsubs);
	    fail_unless(rc == FC_SUCCESS, "failed to create threshold subset");
	    rc = fc_subsetSegmentsMesh(threshsubs,shared_dim,NULL,&numSubsets,NULL);
	    fail_unless(rc == FC_SUCCESS,
			"should work for full input subset and optional args");
	    fail_unless(numSubsets == 0, "should return 0 for full set");
	    rc = fc_subsetSegmentsSubset(threshsubs,whole,shared_dim,NULL,
					 &numSubsets2,NULL);
	    fail_unless(rc == FC_SUCCESS,
			"should work for full input subset and optional args");
	    fail_unless(numSubsets2 == numSubsets, "should return 0 for full set");
	    //clean up

	    rc = fc_subsetPlusNeighborGreaterSegmentationMesh(threshsubs,shared_neighbordim,
							      shared_dim,&numNbr,&nbrIDs);
	    fail_unless(rc == FC_SUCCESS, "should work for full set");
	    fail_unless(numNbr == 0, "should return no neighbors");
	    fail_unless(nbrIDs == NULL, "should return null neighbor set");
	    fc_deleteSubset(threshsubs);
	  }
	  break;
	case 2:
	  {
	    FC_Subset threshsubs; 
	    FC_Subset *newSubsets, *newSubsets2;
	    int numSubsets, numSubsets2,numNbr, *nbrIDs;
	    int ans[2] = {0,2};
	    rc = fc_createThresholdSubset(seqvar[j],"=",1,"temp_thresh", &threshsubs);
	    fail_unless(rc == FC_SUCCESS, "failed to create threshold subset");

	    rc = fc_subsetSegmentsMesh(threshsubs,shared_dim,"live",&numSubsets,&newSubsets);
	    fail_unless(rc == FC_SUCCESS,"failed to subsetSegmentsMesh");
	    fail_unless(numSubsets == 2,"failed to get correct number of subsets");
	    rc = fc_subsetSegmentsSubset(threshsubs,whole,shared_dim,"live",
					 &numSubsets2,&newSubsets2);
	    fail_unless(rc == FC_SUCCESS,"failed to subsetSegmentsMesh");
	    fail_unless(numSubsets2 == numSubsets,
			"failed to get correct number of subsets");
	    for (jj = 0; jj < numSubsets; jj++){
	      int temp_numMember;
	      int *temp_memberIDs;
	      fc_getSubsetMembersAsArray(newSubsets[jj], &temp_numMember, &temp_memberIDs);
	      fail_unless(rc ==  FC_SUCCESS, "unable to get subset members as array");
	      fail_unless(temp_numMember ==1 , "wrong number of members in the subset");
	      fail_unless(temp_memberIDs[0] == ans[jj], "wrong id in this subset");
	      free (temp_memberIDs);

	      fc_getSubsetMembersAsArray(newSubsets2[jj], &temp_numMember, &temp_memberIDs);
	      fail_unless(rc ==  FC_SUCCESS, "unable to get subset members as array");
	      fail_unless(temp_numMember ==1 , "wrong number of members in the subset");
	      fail_unless(temp_memberIDs[0] == ans[jj], "wrong id in this subset");
	      free (temp_memberIDs);
	    }
	    //clean up
	    for (jj = 0; jj < numSubsets; jj++){
	      fc_deleteSubset(newSubsets[jj]);
	      fc_deleteSubset(newSubsets2[jj]);
	    }
	    free(newSubsets);
	    free(newSubsets2);

	    rc = fc_subsetPlusNeighborGreaterSegmentationMesh(threshsubs,shared_neighbordim,
							      shared_dim,&numNbr,&nbrIDs);
	    fail_unless(rc == FC_SUCCESS, "should work");
	    fail_unless(numNbr == 0, "should return no neighbors");
	    fail_unless(nbrIDs == NULL, "should return null neighbor set");

	    fc_deleteSubset(threshsubs);
	  }
	  break;
	default:
	  break;
	}
	break;
      case FC_ET_QUAD:
	switch (j){
	case 0:
	  {
	    FC_Subset threshsubs; 
	    FC_Subset *newSubsets, *newSubsets2;
	    int numSubsets, numSubsets2, numNbr, *nbrIDs;
	    rc = fc_createThresholdSubset(seqvar[j],"=",1,"temp_thresh", &threshsubs);
	    fail_unless(rc == FC_SUCCESS, "failed to create threshold subset");
	    rc = fc_subsetSegmentsMesh(threshsubs,shared_dim,"live",&numSubsets,&newSubsets);
	    fail_unless(rc == FC_SUCCESS, "should work for empty input subset");
	    fail_unless(numSubsets == 1, "empty set should return same seg as full mesh");
	    rc = fc_subsetSegmentsSubset(threshsubs,whole,shared_dim,"live",
					 &numSubsets2,&newSubsets2);
	    fail_unless(rc == FC_SUCCESS, "should work for empty input subset");
	    fail_unless(numSubsets2 == numSubsets, "empty set should return same seg as full mesh");
	    //could check to see that they are the same subsets here
	    //clean up
	    for (jj = 0; jj < numSubsets; jj++){
	      fc_deleteSubset(newSubsets[jj]);
	      fc_deleteSubset(newSubsets2[jj]);
	    }
	    free(newSubsets);
	    free(newSubsets2);

	    rc = fc_subsetPlusNeighborGreaterSegmentationMesh(threshsubs,shared_neighbordim,
							      shared_dim,&numNbr,&nbrIDs);
	    fail_unless(rc != FC_SUCCESS, "should not work for empty input subset");
	    fc_deleteSubset(threshsubs);
	  }
	  break;
	case 1:
	  {
	    FC_Subset threshsubs; 
	    FC_Subset *newSubsets, *newSubsets2;
	    int numSubsets, numSubsets2,numNbr, *nbrIDs;
	    rc = fc_createThresholdSubset(seqvar[j],"=",1,"temp_thresh", &threshsubs);
	    fail_unless(rc == FC_SUCCESS, "failed to create threshold subset");
	    rc = fc_subsetSegmentsMesh(threshsubs,shared_dim,"live",&numSubsets,&newSubsets);
	    fail_unless(rc == FC_SUCCESS, "should work for full input subset");
	    fail_unless(numSubsets == 0, "should return 0 for full set");
	    rc = fc_subsetSegmentsSubset(threshsubs,whole,shared_dim,"live",&numSubsets2,&newSubsets2);
	    fail_unless(rc == FC_SUCCESS, "should work for full input subset");
	    fail_unless(numSubsets2 == numSubsets, "should return 0 for full set");
	    //clean up
	    for (jj = 0; jj < numSubsets; jj++){
	      fc_deleteSubset(newSubsets[jj]);
	      fc_deleteSubset(newSubsets2[jj]);
	    }
	    free(newSubsets);
	    rc = fc_subsetPlusNeighborGreaterSegmentationMesh(threshsubs,shared_neighbordim,
							      shared_dim,&numNbr,&nbrIDs);
	    fail_unless(rc == FC_SUCCESS, "should work for full set");
	    fail_unless(numNbr == 0, "should return no neighbors");
	    fail_unless(nbrIDs == NULL, "should return null neighbor set");
	    fc_deleteSubset(threshsubs);
	  }
	  break;
	case 2:
	  {
	    int ans[5] = {0,3,6,7,8};
	    FC_Subset threshsubs; 
	    FC_Subset *newSubsets, *newSubsets2;
	    int numSubsets, numSubsets2, temp_numMember;
	    int numNbr, *nbrIDs;
	    int jk, found = 0;
	    int *temp_memberIDs;
	    rc = fc_createThresholdSubset(seqvar[j],"=",1,"temp_thresh", &threshsubs);
	    fail_unless(rc == FC_SUCCESS, "failed to create threshold subset");
	    rc = fc_subsetSegmentsMesh(threshsubs,shared_dim,"live",&numSubsets,&newSubsets);
	    fail_unless(rc == FC_SUCCESS, "failed to segment subset");
	    fail_unless(numSubsets == 2 ,
			"failed to get correct number of subsets");
	    rc = fc_subsetSegmentsSubset(threshsubs,whole,shared_dim,"live",
					 &numSubsets2,&newSubsets2);
	    fail_unless(rc == FC_SUCCESS, "failed to segment subset");
	    fail_unless(numSubsets2 == numSubsets ,
			"failed to get correct number of subsets");
	    //first subset
	    rc = fc_getSubsetMembersAsArray(newSubsets[0],
					    &temp_numMember, &temp_memberIDs);
	    fail_unless(rc ==  FC_SUCCESS, "unable to get subset members as array");
	    fail_unless(temp_numMember == 5 , "wrong number of members in the subset");
	    for (jk = 0; jk < 5; jk++){
	      for (jj = 0; jj < 5; jj++){
		if(ans[jk] == temp_memberIDs[jj]){
		  found = 1;
		  break;
		}
	      }
	      fail_unless(found, "missing id in this subset");
	      found = 0;
	    }
	    free (temp_memberIDs);
	    rc = fc_getSubsetMembersAsArray(newSubsets2[0],
					    &temp_numMember, &temp_memberIDs);
	    fail_unless(rc ==  FC_SUCCESS, "unable to get subset members as array");
	    fail_unless(temp_numMember == 5 , "wrong number of members in the subset");
	    for (jk = 0; jk < 5; jk++){
	      for (jj = 0; jj < 5; jj++){
		if(ans[jk] == temp_memberIDs[jj]){
		  found = 1;
		  break;
		}
	      }
	      fail_unless(found, "missing id in this subset");
	      found = 0;
	    }
	    free (temp_memberIDs);
	    //second subset
	    rc = fc_getSubsetMembersAsArray(newSubsets[1],
					    &temp_numMember, &temp_memberIDs);
	    fail_unless(rc ==  FC_SUCCESS, "unable to get subset members as array");
	    fail_unless(temp_numMember == 1 , "wrong number of members in the subset");
	    fail_unless(temp_memberIDs[0] ==  2, "wrong id in this subset");
	    free (temp_memberIDs);
	    rc = fc_getSubsetMembersAsArray(newSubsets2[1],
					    &temp_numMember, &temp_memberIDs);
	    fail_unless(rc ==  FC_SUCCESS, "unable to get subset members as array");
	    fail_unless(temp_numMember == 1 , "wrong number of members in the subset");
	    fail_unless(temp_memberIDs[0] ==  2, "wrong id in this subset");
	    free (temp_memberIDs);
	    
	    //clean up
	    for (jj = 0; jj < numSubsets; jj++){
	      fc_deleteSubset(newSubsets[jj]);
	      fc_deleteSubset(newSubsets2[jj]);
	    }
	    free(newSubsets);
	    free(newSubsets2);
	    
	    //this set is a case where the assumptions of the recursive case fail
	    rc = fc_subsetPlusNeighborGreaterSegmentationMesh(threshsubs,shared_neighbordim,
							      shared_dim,&numNbr,&nbrIDs);
	    fail_unless(rc == FC_SUCCESS, "should work");
	    fail_unless(numNbr == 0, "should return no neighbors");
	    fail_unless(nbrIDs == NULL, "should return null neighbor set");


	    //try diff case
	    fc_addMemberToSubset(threshsubs,7);
	    fc_deleteMemberFromSubset(threshsubs,4);
	    fc_deleteMemberFromSubset(threshsubs,5);
	    rc = fc_subsetPlusNeighborGreaterSegmentationMesh(threshsubs,1,
							      1,&numNbr,&nbrIDs);
	    
	    fail_unless(rc==FC_SUCCESS,"neighbor seg should work");
	    fail_unless(numNbr==1,"wrong number of neighbors");
	    fail_unless(nbrIDs[0]==4,"wrong number of neighbors");
	    free(nbrIDs);
	    
	    fc_deleteSubset(threshsubs);
	  }
	  break;
	default:
	  break;
	}
	break;
      case FC_ET_HEX:
	switch (j){
	case 0:
	  {
	    FC_Subset threshsubs; 
	    FC_Subset *newSubsets, *newSubsets2;
	    int numSubsets, numSubsets2,numNbr, *nbrIDs;
	    rc = fc_createThresholdSubset(seqvar[j],"=",1,"temp_thresh", &threshsubs);
	    fail_unless(rc == FC_SUCCESS, "failed to create threshold subset");
	    rc = fc_subsetSegmentsMesh(threshsubs,shared_dim,"live",&numSubsets,&newSubsets);
	    fail_unless(rc == FC_SUCCESS, "should work for empty input subset");
	    fail_unless(numSubsets == 1, "empty set should return same seg as full mesh");
	    rc = fc_subsetSegmentsSubset(threshsubs,whole,shared_dim,"live",
					 &numSubsets2,&newSubsets2);
	    fail_unless(rc == FC_SUCCESS, "should work for empty input subset");
	    fail_unless(numSubsets2 == numSubsets,
			"empty set should return same seg as full mesh");
	    //could check to see that they are the same subsets here
	    //clean up
	    for (jj = 0; jj < numSubsets; jj++){
	      fc_deleteSubset(newSubsets[jj]);
	      fc_deleteSubset(newSubsets2[jj]);
	    }
	    free(newSubsets);
	    free(newSubsets2);
	    rc = fc_subsetPlusNeighborGreaterSegmentationMesh(threshsubs,shared_neighbordim,
							      shared_dim,&numNbr,&nbrIDs);
	    fail_unless(rc != FC_SUCCESS, "should not work for empty input subset");
	    fc_deleteSubset(threshsubs);
	  }
	  break;
	case 1:
	  {
	    FC_Subset threshsubs; 
	    FC_Subset *newSubsets, *newSubsets2;
	    int numSubsets, numSubsets2,numNbr, *nbrIDs;
	    rc = fc_createThresholdSubset(seqvar[j],"=",1,"temp_thresh", &threshsubs);
	    fail_unless(rc == FC_SUCCESS, "failed to create threshold subset");
	    rc = fc_subsetSegmentsMesh(threshsubs,shared_dim,"live",
				       &numSubsets,&newSubsets);
	    fail_unless(rc == FC_SUCCESS, "should work for full input subset");
	    fail_unless(numSubsets == 0, "should return 0 for full set");
	    rc = fc_subsetSegmentsSubset(threshsubs,whole,shared_dim,"live",
					 &numSubsets2,&newSubsets2);
	    fail_unless(rc == FC_SUCCESS, "should work for full input subset");
	    fail_unless(numSubsets2 == numSubsets, "should return 0 for full set");
	    //could check to see that they are the same subsets here
	    //clean up
	    for (jj = 0; jj < numSubsets; jj++){
	      fc_deleteSubset(newSubsets[jj]);
	      fc_deleteSubset(newSubsets2[jj]);
	    }
	    free(newSubsets);
	    free(newSubsets2);
	    rc = fc_subsetPlusNeighborGreaterSegmentationMesh(threshsubs,shared_neighbordim,
							      shared_dim,&numNbr,&nbrIDs);
	    fail_unless(rc == FC_SUCCESS, "should work for full set");
	    fail_unless(numNbr == 0, "should return no neighbors");
	    fail_unless(nbrIDs == NULL, "should return null neighbor set");
	    fc_deleteSubset(threshsubs);
	  }
	  break;
	case 2:
	  {
	    int ans[9] = {0,1,2,3,4,5,6,7,8};
	    int ansnbr[3] = {10,13,16};
	    FC_Subset threshsubs; 
	    FC_Subset *newSubsets,*newSubsets2;
	    int numSubsets, numSubsets2,temp_numMember;
	    int jk,subs, found = 0;
	    int *temp_memberIDs;
	    int numNbr, *nbrIDs;
	    rc = fc_createThresholdSubset(seqvar[j],"=",1,"temp_thresh", &threshsubs);
	    fail_unless(rc == FC_SUCCESS, "failed to create threshold subset");
	    rc = fc_subsetSegmentsMesh(threshsubs,shared_dim,"live",&numSubsets,&newSubsets);
	    fail_unless(rc == FC_SUCCESS, "failed to segment subset");
	    fail_unless(numSubsets == 2,
			"failed to get correct number of subsets");
	    rc = fc_subsetSegmentsSubset(threshsubs,whole,shared_dim,"live",
					 &numSubsets2,&newSubsets2);
	    fail_unless(rc == FC_SUCCESS, "failed to segment subset");
	    fail_unless(numSubsets == 2,
			"failed to get correct number of subsets");

	    for (subs = 0; subs < 2; subs++){
	      //		fc_printSubset(newSubsets[subs],"",1);
	      fc_getSubsetMembersAsArray(newSubsets[subs], &temp_numMember,
					 &temp_memberIDs);
	      fail_unless(rc ==  FC_SUCCESS, "unable to get subset members as array");
	      fail_unless(temp_numMember == 9 , "wrong number of members in the subset");
	      for (jk = 0; jk < 9; jk++){
		for (jj = 0; jj < 9; jj++){
		  if( (subs == 0 ? ans[jk] : ans[jk]+18)  == temp_memberIDs[jj]){
		    found = 1;
		    break;
		  }
		}
		fail_unless(found, "missing id in this subset");
		found = 0;
	      }
	      free (temp_memberIDs);
	      fc_getSubsetMembersAsArray(newSubsets2[subs], &temp_numMember,
					 &temp_memberIDs);
	      fail_unless(rc ==  FC_SUCCESS, "unable to get subset members as array");
	      fail_unless(temp_numMember == 9 , "wrong number of members in the subset");
	      for (jk = 0; jk < 9; jk++){
		for (jj = 0; jj < 9; jj++){
		  if( (subs == 0 ? ans[jk] : ans[jk]+18)  == temp_memberIDs[jj]){
		    found = 1;
		    break;
		  }
		}
		fail_unless(found, "missing id in this subset");
		found = 0;
	      }
	      free (temp_memberIDs);
	    }

	    //check names here - only doing this this one time
	    for (jj = 0; jj < numSubsets; jj++){
	      char testname[20];
	      char *newname;
	      rc = fc_getSubsetName(newSubsets[jj],&newname);
	      fail_unless(rc == FC_SUCCESS, "failed to get new subsetname");
	      sprintf(testname, "%s%d","live_Seg",jj);
	      fail_unless(strcmp(testname,newname) == 0,
			  "bad name for new subset");
	      free(newname);
	      rc = fc_getSubsetName(newSubsets2[jj],&newname);
	      fail_unless(rc == FC_SUCCESS, "failed to get new subsetname");
	      sprintf(testname, "%s%d","live_Seg",jj);
	      fail_unless(strcmp(testname,newname) == 0,
			  "bad name for new subset");
	      free(newname);
	    }

	    //clean up
	    for (jj = 0; jj < numSubsets; jj++){
	      fc_deleteSubset(newSubsets[jj]);
	      fc_deleteSubset(newSubsets2[jj]);
	    }
	    free(newSubsets);
	    free(newSubsets2);

	    
	    // this naturally a special case where the recursive assumption fails.
	    //change the case. note this is still a prob if dont do the rotation
	    //becuase of the way the numbers fall.
	    fc_addMemberToSubset(threshsubs,3);
	    fc_addMemberToSubset(threshsubs,4);
	    fc_addMemberToSubset(threshsubs,5);
	    fc_addMemberToSubset(threshsubs,21);
	    fc_addMemberToSubset(threshsubs,22);
	    fc_addMemberToSubset(threshsubs,23);
	    for (jj =9; jj < 18; jj++){
	      fc_deleteMemberFromSubset(threshsubs,jj);
	    }
	    fc_addMemberToSubset(threshsubs,12);
	    fc_addMemberToSubset(threshsubs,14);
	    fc_deleteMemberFromSubset(threshsubs,13);
	    
	    rc = fc_subsetPlusNeighborGreaterSegmentationMesh(threshsubs,1,
							      2,&numNbr,&nbrIDs);
	    fail_unless(rc==FC_SUCCESS,"neighbor seg should work");
	    fail_unless(numNbr==3,"wrong number of neighbors");
	    found = 0;
	    for (jk = 0; jk < 3; jk++){
	      for (jj = 0; jj < 3; jj++){
		if(nbrIDs[jk] == ansnbr[jj]){
		  found = 1;
		  break;
		}
	      }
	      fail_unless(found, "missing id in return vals");
	      found = 0;
	    }
	    free(nbrIDs);

	    fc_deleteSubset(threshsubs);

	    //sticking this in here - non intersecting subsets
	    {
	      FC_Subset a,b;

	      rc = fc_createSubset(mesh,"a",FC_AT_ELEMENT,&a);
	      fail_unless(rc == FC_SUCCESS, "cant make new subset");
	      rc = fc_addMemberToSubset(a,1);
	      fail_unless(rc == FC_SUCCESS, "cant add member to subset");
	      rc = fc_createSubset(mesh,"b",FC_AT_ELEMENT,&b);
	      fail_unless(rc == FC_SUCCESS, "cant make new subset");
	      rc = fc_addMemberToSubset(b,2);
	      fail_unless(rc == FC_SUCCESS, "cant add member to subset");

	      rc = fc_subsetSegmentsSubset(a,b,shared_dim,"live",
					 &numSubsets,&newSubsets);
	      fail_unless(rc == FC_SUCCESS, "failed to segement subset");
	      fail_unless(numSubsets == 1,
			  "failed to get correct number of segments");
	      fc_deleteSubset(newSubsets[0]);
	      free(newSubsets);
	      fc_deleteSubset(a);
	      fc_deleteSubset(b);
	    }
	  }
	  break;
	default:
	  break;
	}
	break;
      default:
	//shouldnt occur
	break;
      }
    } //numstep



    // --- Test bad args 
    {
      FC_Subset threshsubs, *newSubsets;
      FC_Subset emptySubset;
      FC_Subset badSubset= {999,999};
      int numSubsets, numNbr, *nbrIDs;
      
      rc = fc_createThresholdSubset(seqvar[0],"=",1,"temp_thresh", &threshsubs);
      fail_unless(rc == FC_SUCCESS, "failed to create threshold subset");
      //note some empty cases were checkes above, but im explictly testing
      //the bad ones here
      rc = fc_createSubset(mesh,"empty",FC_AT_ELEMENT,&emptySubset);
      fail_unless(rc == FC_SUCCESS, "cant make empty subset to test");

      if (i == 3){//hex mesh - need case with faces
	FC_Subset wrongassoc;
	FC_Subset wrongmesh;

	rc = fc_createSubset(mesh,"wrongassoc",FC_AT_FACE,&wrongassoc);
	fail_unless(rc == FC_SUCCESS, "cant make new assoc to test");
	rc = fc_addMemberToSubset(wrongassoc,1);
	fail_unless(rc == FC_SUCCESS, "cant add member to subset");

	rc = fc_createSubset(meshes[1],"wrongmesh",FC_AT_ELEMENT,&wrongmesh);
	fail_unless(rc == FC_SUCCESS, "cant make new subset to test");
	rc = fc_addMemberToSubset(wrongmesh,1);
	fail_unless(rc == FC_SUCCESS, "cant add member to subset");

	/*
	rc = fc_subsetSegmentsMesh(wrongassoc,shared_dim,"live",&numSubsets,
				   &newSubsets);
	fail_unless(rc != FC_SUCCESS, "Should fail for wrong assoc subset");
	fail_unless(newSubsets == NULL,
		    "fail should return null array");
	fail_unless(numSubsets == -1, "fail should return -1");
	
	rc = fc_subsetSegmentsSubset(whole,wrongassoc,shared_dim,"live",
				     &numSubsets, &newSubsets);
	fail_unless(rc != FC_SUCCESS, "Should fail for wrong assoc subset");
	fail_unless(newSubsets == NULL,
		    "fail should return null array");
	fail_unless(numSubsets == -1, "fail should return -1");

	rc = fc_subsetSegmentsSubset(wrongassoc,whole,shared_dim,"live",
				     &numSubsets, &newSubsets);
	fail_unless(rc != FC_SUCCESS, "Should fail for wrong assoc subset");
	fail_unless(newSubsets == NULL,
		    "fail should return null array");
	fail_unless(numSubsets == -1, "fail should return -1");
	*/

	rc = fc_subsetPlusNeighborGreaterSegmentationMesh(wrongassoc,
							  shared_neighbordim,
							  shared_dim,&numNbr,
							  &nbrIDs);
	fail_unless(rc != FC_SUCCESS, "Should fail for wrong assoc subset");
	fail_unless(nbrIDs == NULL,
		    "fail should return null array");
	fail_unless(numNbr == -1, "fail should return -1");

	rc = fc_subsetSegmentsSubset(whole,wrongmesh,shared_dim,"live",
				     &numSubsets, &newSubsets);
	fail_unless(rc != FC_SUCCESS, "Should fail for mismatched mesh subset");
	fail_unless(newSubsets == NULL,
		    "fail should return null array");
	fail_unless(numSubsets == -1, "fail should return -1");

	//clean up
	rc = fc_deleteSubset(wrongassoc);
      }

      rc = fc_subsetSegmentsSubset(threshsubs,emptySubset,shared_dim,
				   "live", &numSubsets,&newSubsets);
      fail_unless(rc != FC_SUCCESS, "Should fail for empty outer subset");
      fail_unless(newSubsets == NULL,
		  "fail should return null array");
      fail_unless(numSubsets == -1, "fail should return -1");

      rc = fc_subsetSegmentsSubset(threshsubs,FC_NULL_SUBSET,shared_dim,
				   "live", &numSubsets,&newSubsets);
      fail_unless(rc != FC_SUCCESS, "Should fail for NULL outer subset");
      fail_unless(newSubsets == NULL,
		  "fail should return null array");
      fail_unless(numSubsets == -1, "fail should return -1");

      rc = fc_subsetSegmentsMesh(badSubset,shared_dim,"live",
				 &numSubsets,&newSubsets);
      fail_unless(rc != FC_SUCCESS, "Should fail for bad subset");
      fail_unless(newSubsets == NULL,
		  "fail should return null array");
      fail_unless(numSubsets == -1, "fail should return -1");


      rc = fc_subsetSegmentsSubset(badSubset,whole, shared_dim,"live",
				   &numSubsets,&newSubsets);
      fail_unless(rc != FC_SUCCESS, "Should fail for bad subset");
      fail_unless(newSubsets == NULL,
		  "fail should return null array");
      fail_unless(numSubsets == -1, "fail should return -1");

      rc = fc_subsetSegmentsSubset(whole,badSubset,shared_dim,"live",
				   &numSubsets,&newSubsets);
      fail_unless(rc != FC_SUCCESS, "Should fail for bad subset");
      fail_unless(newSubsets == NULL,
		  "fail should return null array");
      fail_unless(numSubsets == -1, "fail should return -1");


      rc = fc_subsetPlusNeighborGreaterSegmentationMesh(badSubset,
							shared_neighbordim,
							shared_dim,&numNbr,
							&nbrIDs);
      fail_unless(rc != FC_SUCCESS, "Should fail for bad subset");
      fail_unless(nbrIDs == NULL,
		  "fail should return null array");
      fail_unless(numNbr == -1, "fail should return -1");


      rc = fc_subsetPlusNeighborGreaterSegmentationMesh(threshsubs,-1,
							shared_dim,
							&numNbr,&nbrIDs);
      fail_unless(rc != FC_SUCCESS, "Should fail for bad neighbordim");
      fail_unless(nbrIDs == NULL,
		  "fail should return null array");
      fail_unless(numNbr == -1, "fail should return -1");

      rc = fc_subsetPlusNeighborGreaterSegmentationMesh(threshsubs,3,
							shared_dim,
							&numNbr,&nbrIDs);
      fail_unless(rc != FC_SUCCESS, "Should fail for bad neighbordim");
      fail_unless(nbrIDs == NULL,
		  "fail should return null array");
      fail_unless(numNbr == -1, "fail should return -1");

      //should check shared segdim here. also if dims of both make sense given the
      //subset. in the function they are just passed into other methods that check
      //them so maybe its unnecessary...

      rc = fc_subsetSegmentsMesh(threshsubs,shared_dim,"live",
				 &numSubsets,NULL);
      fail_unless(rc != FC_SUCCESS,
		  "Should fail for NULL subset but non-NULL name");
      fail_unless(newSubsets == NULL,
		  "fail should return null array");
      fail_unless(numSubsets == -1, "fail should return -1");

      rc = fc_subsetSegmentsSubset(threshsubs,whole,shared_dim,"live",
				   &numSubsets,NULL);
      fail_unless(rc != FC_SUCCESS,
		  "Should fail for NULL subset but non-NULL name");
      fail_unless(newSubsets == NULL,
		  "fail should return null array");
      fail_unless(numSubsets == -1, "fail should return -1");

      rc = fc_subsetSegmentsMesh(threshsubs,shared_dim,NULL,&numSubsets,
				 &newSubsets);
      fail_unless(rc != FC_SUCCESS,
		  "Should fail for non-NULL subset but NULL name");
      fail_unless(newSubsets == NULL,
		  "fail should return null array");
      fail_unless(numSubsets == -1, "fail should return -1");

      rc = fc_subsetSegmentsSubset(threshsubs,whole,shared_dim,NULL,&numSubsets,
				   &newSubsets);
      fail_unless(rc != FC_SUCCESS,
		  "Should fail for non-NULL subset but NULL name");
      fail_unless(newSubsets == NULL,
		  "fail should return null array");
      fail_unless(numSubsets == -1, "fail should return -1");

      // --- cleanup
      rc = fc_deleteSubset(threshsubs);
      rc = fc_deleteSubset(emptySubset);
    } //bad args

    // --- cleanup
    rc = fc_deleteSeqVariable(numStep,seqvar);
    free(seqvar);
    fc_deleteSubset(whole);

  }
  for (i = 0; i < numElemType; i++){
    fc_deleteMesh(meshes[i]); // also cleans up subset
  }
  fc_deleteDataset(dataset);
}
END_TEST


//assumes that shape works, not other shape methods
START_TEST(decay){
  FC_ReturnCode rc;
  FC_Dataset dataset;
  FC_Mesh mesh;

  //coord info for the mesh
  FC_Coords lowers = { 0., 0., 0. };
  FC_Coords uppers = { 1., 1., 1. };
  int numElemPerType = 27;
  FC_Subset deadregion, decayedSkin, decayedSkin2,meshskin;
  FC_Subset *decayedSides;
  FC_Subset emptyf, emptye, whole, keep;
  FC_Shape *shapes, badshape;
  int numShapes;
  FC_SortedIntArray sia;

  FC_Subset badSubset={999,999};
  int *arr, *arr2, *arr3;
  int decayflag, *decayflagarray;
  int i,j,k;

  
  
  // setup - create test mesh
  rc = fc_createDataset("dset", &dataset);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create dataset");

  rc = fc_createSimpleHexMesh(dataset,"hex_mesh",3,3,3,lowers,uppers,&mesh);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create mesh");

  rc = fc_getMeshSkin(mesh,&meshskin);
  fail_unless(rc == FC_SUCCESS,"abort: failed to get mesh skin");

  rc = fc_getMeshShapes(mesh,30,0,&numShapes,&shapes);
  fail_unless(rc == FC_SUCCESS, "failed to surface segment whole");
  fail_unless(numShapes == 1, "wrong number of shapes");


  // --- empty subset -----/
  rc = fc_createSubset(mesh,"emptye",FC_AT_ELEMENT,&emptye);
  fc_exitIfErrorPrintf(rc, "failed to make subset");

  rc = fc_createSubset(mesh,"emptyf",FC_AT_FACE,&emptyf);
  fc_exitIfErrorPrintf(rc, "failed to make subset");

  rc = fc_getDecayedSkin(emptye,&meshskin,&decayedSkin);
  fail_unless(rc == FC_SUCCESS, "failed to get decayedSkin");
  fail_unless(FC_HANDLE_EQUIV(decayedSkin, FC_NULL_SUBSET),
	      "bad return val decayedSkin");

  rc = fc_getDecayedShapeSides(emptye,&shapes[0], &decayedSides);
  fail_unless(rc == FC_SUCCESS, "failed to get decayedShapeSides");
  for (i = 0; i < shapes[0].numSides; i++){
    fail_unless(FC_HANDLE_EQUIV(decayedSides[i], FC_NULL_SUBSET),
		"bad return val decayedSkin");
  }
  free(decayedSides);
    
  rc = fc_getDecayedShapeSkin(emptye,&shapes[0], &decayedSkin);
  fail_unless(rc == FC_SUCCESS, "failed to get decayedSkin");
  fail_unless(FC_HANDLE_EQUIV(decayedSkin, FC_NULL_SUBSET),
	      "bad return val decayedSkin");

  // --- complete mesh -----/
  rc = fc_createSubset(mesh, "whole", FC_AT_ELEMENT, &whole);
  fail_unless(rc == FC_SUCCESS,"cant create whole subset"); 
  for (j = 0; j < numElemPerType; j++){
    fc_addMemberToSubset(whole,j);
  }

  rc = fc_getDecayedSkin(whole,&meshskin,&decayedSkin);
  fail_unless(rc == FC_SUCCESS, "failed to get decayedSkin");
  rc = fc_getDecayedShapeSkin(whole,&shapes[0], &decayedSkin2);
  fail_unless(rc == FC_SUCCESS, "failed to get decayedSkin");
  rc = fc_getSubsetMembersAsArray(decayedSkin,&i,&arr);
  fail_unless(rc == FC_SUCCESS, "cant get decayedskin members");
  rc = fc_getSubsetMembersAsArray(decayedSkin2,&j,&arr2);
  fail_unless(rc == FC_SUCCESS, "cant get decayedskin2 members");
  rc = fc_getSubsetMembersAsArray(meshskin,&k,&arr3);
  fail_unless(rc == FC_SUCCESS, "cant get mesh skin members");
  fail_unless(i == k, "wrong num members in the decayedSkin");
  fail_unless(j == k, "wrong num members in the decayedSkin");
  for (i = 0; i <j ; i++){
    //subset array are in ascending order
    fail_unless(arr[i] == arr3[i], "wrong id in decayed side");
    fail_unless(arr2[i] == arr3[i], "wrong id in decayed side");
  }
  free(arr);
  free(arr2);
  free(arr3);
  rc = fc_deleteSubset(decayedSkin);
  rc = fc_deleteSubset(decayedSkin2);


  rc = fc_getDecayedShapeSides(whole,&shapes[0], &decayedSides);
  fail_unless(rc == FC_SUCCESS, "failed to get decayedShapeSides");
  //each whole side should be decayed
  for (i = 0; i < shapes[0].numSides; i++){
    rc = fc_getSubsetMembersAsArray(shapes[0].faces[i],&j,&arr);
    fail_unless(rc == FC_SUCCESS, "cant get members of the side");
    rc = fc_getSubsetMembersAsArray(decayedSides[i],&k,&arr2);
    fail_unless(rc == FC_SUCCESS, "cant get members of the decayed shape side");
    fail_unless(j == k, "wrong number of members in the decayed side");
    for (j = 0; j <k; j++){
      //subset array are in ascending order
      fail_unless(arr[j] == arr2[j], "wrong id in decayed side");
    }
    free(arr);
    free(arr2);
    fc_deleteSubset(decayedSides[i]);
  }
  free(decayedSides);

  // --- internal element --- /
  //should be no intersections

  rc = fc_createSubset(mesh, "middle", FC_AT_ELEMENT, &deadregion);
  rc = fc_addMemberToSubset(deadregion,13);
  fail_unless(rc == FC_SUCCESS, "failed to add member to subset");

  for (i= 0 ; i < shapes[0].numSides; i++){
    rc = fc_getDecayedSkin(deadregion,&shapes[0].faces[i],&decayedSkin);
    fail_unless(rc == FC_SUCCESS, "failed to get decayedSkin");
    rc = fc_getSubsetMembersAsArray(decayedSkin,&j,&arr);
    fail_unless(j == 0 , "wrong number of decays in this side");
    free(arr);
    fc_deleteSubset(decayedSkin);
  }

  rc = fc_getDecayedShapeSides(deadregion,&shapes[0], &decayedSides);
  fail_unless(rc == FC_SUCCESS, "failed to get decayedShapeSides");
  for (i = 0 ; i < shapes[0].numSides; i++){
    rc = fc_getSubsetMembersAsArray(decayedSides[i],&j,&arr);
    fail_unless(j == 0, "wrong number of decays this side");
    free(arr);
    fc_deleteSubset(decayedSides[i]);
  }
  free(decayedSides);

   
  rc = fc_getDecayedShapeSkin(deadregion,&shapes[0], &decayedSkin);
  fail_unless(rc == FC_SUCCESS, "failed to get decayedSkin");
  rc = fc_getSubsetMembersAsArray(decayedSkin,&i,&arr);
  fail_unless(i == 0, "wrong number of members in decayedskin");
  fc_deleteSubset(decayedSkin);
  free(arr);

   
  rc = fc_getDecayedShapeSkin(deadregion,&shapes[0], &decayedSkin);
  fail_unless(rc == FC_SUCCESS, "failed to get decayedSkin");
  rc = fc_getSubsetMembersAsArray(decayedSkin,&i,&arr);
  fail_unless(i == 0, "wrong number of members in decayedskin");
  fc_deleteSubset(decayedSkin);
  free(arr);

  // --- case with creation of mesh skin -----/
  rc = fc_getDecayedSkin(deadregion,NULL, &decayedSkin);
  fail_unless(rc == FC_SUCCESS, "failed to get decayedSkin");
  rc = fc_getSubsetMembersAsArray(decayedSkin,&i,&arr);
  fail_unless(i == 0, "wrong number of members in decayedskin");
  fc_deleteSubset(decayedSkin);
  free(arr);

  fc_deleteSubset(deadregion);


  // --- single element -----/
  //since these calls are essentially wrappers around subsetIntersections
  //i wont check the vals
  rc = fc_createSubset(mesh, "middle", FC_AT_ELEMENT, &deadregion);
  rc = fc_addMemberToSubset(deadregion,4);
  fail_unless(rc == FC_SUCCESS, "failed to add member to subset");


  rc = fc_getDecayedShapeSides(deadregion,&shapes[0], &decayedSides);
  fail_unless(rc == FC_SUCCESS, "failed to get decayedShapeSides");
  for (i = 0 ; i < shapes[0].numSides; i++){
    rc = fc_getDecayedSkin(deadregion,&shapes[0].faces[i],&decayedSkin);
    fail_unless(rc == FC_SUCCESS, "failed to get decayedSkin");
    rc = fc_getSubsetNumMember(decayedSkin,&j);
    fail_unless(rc == FC_SUCCESS, "can't get subset num memebrs");
    rc = fc_doSubsetsIntersect(decayedSkin,shapes[0].faces[i],&k);
    fail_unless(rc == FC_SUCCESS, "cant determine if subsets intersect");
    fail_unless(j == k,"wrong return for decayed skin on this side");
    rc = fc_getSubsetMembersAsArray(decayedSides[i],&k,&arr);
    fail_unless(j == k,"wrong return for decayed skin on this side");
    free(arr);
    fc_deleteSubset(decayedSkin);  
    fc_deleteSubset(decayedSides[i]);
  }
  free(decayedSides);
   

  //not chekcing hte val since its essentially a wrapper
  rc = fc_getDecayedShapeSkin(deadregion,&shapes[0], &decayedSkin);
  fail_unless(rc == FC_SUCCESS, "failed to get decayedSkin");
  rc = fc_getSubsetMembersAsArray(decayedSkin,&i,&arr);
  fail_unless(i == 1, "wrong number of members in decayedskin");
  fc_deleteSubset(decayedSkin);
  free(arr);

  // --- case with creation of mesh skin -----/
  rc = fc_getDecayedSkin(deadregion,NULL, &decayedSkin);
  fail_unless(rc == FC_SUCCESS, "failed to get decayedSkin");
  rc = fc_getSubsetMembersAsArray(decayedSkin,&i,&arr);
  fail_unless(i == 1, "wrong number of members in decayedskin");
  fc_deleteSubset(decayedSkin);
  free(arr);


  // --- thru the middle -----/
  //since these calls are essentially wrappers around subsetIntersections
  //i wont check the vals
  rc = fc_createSubset(mesh, "middle", FC_AT_ELEMENT, &deadregion);
  rc = fc_addMemberToSubset(deadregion,4);
  fail_unless(rc == FC_SUCCESS, "failed to add member to subset");
  rc = fc_addMemberToSubset(deadregion,13);
  fail_unless(rc == FC_SUCCESS, "failed to add member to subset");
  rc = fc_addMemberToSubset(deadregion,22);
  fail_unless(rc == FC_SUCCESS, "failed to add member to subset");


  //compare the dead region to all sides, only 2 of those should
  //have intersections. known numbers for the intersections.
  rc = fc_initSortedIntArray(&sia);
  fail_unless(rc == FC_SUCCESS, "can't init sorted int array");
  for (i= 0 ; i < shapes[0].numSides; i++){
    rc = fc_getDecayedSkin(deadregion,&shapes[0].faces[i],&decayedSkin);
    fail_unless(rc == FC_SUCCESS, "failed to get decayedSkin");
    rc = fc_getSubsetMembersAsArray(decayedSkin,&j,&arr);
    fail_unless((j == 0 || j == 1), "wrong number of decays in this side");
    //since this call is a wrapper around getting intersection
    //i wont check the results - for now, just check that they are non
    //adj sides
    if (j ==1){
      fc_addIntToSortedIntArray(&sia,i);
      fail_unless(rc == FC_SUCCESS, "failed to add int to sorted int array"); 
    }
    free(arr);
    fc_deleteSubset(decayedSkin);
  }
  fail_unless(sia.numVal == 2, "should only be 2 intersections");
  fail_unless(shapes[0].adjmatrix[sia.vals[0]][sia.vals[1]] == 0,
	      "should intersect non adj sides");

  //now whole shape call
  rc = fc_getDecayedShapeSides(deadregion,&shapes[0], &decayedSides);
  fail_unless(rc == FC_SUCCESS, "failed to get decayedShapeSides");
  //compare the dead region to all sides, only 2 of those should
  //have intersections. known numbers for the intersections.
  for (i = 0 ; i < shapes[0].numSides; i++){
    rc = fc_getSubsetMembersAsArray(decayedSides[i],&j,&arr);
    if (fc_isIntInSortedIntArray(&sia,i)){
      fail_unless(j == 1, "wrong number of decays this side");
      //do we want to check the vals ?
    }else{
      fail_unless(j == 0, "wrong number of decays this side");
    }
    free(arr);
    fc_deleteSubset(decayedSides[i]);
  }
  free(decayedSides);

   
  rc = fc_getDecayedShapeSkin(deadregion,&shapes[0], &decayedSkin);
  fail_unless(rc == FC_SUCCESS, "failed to get decayedSkin");
  rc = fc_getSubsetMembersAsArray(decayedSkin,&i,&arr);
  fail_unless(i == 2, "wrong number of members in decayedskin");
  //do we want to check the vals ?
  fc_deleteSubset(decayedSkin);
  free(arr);

  // --- case with creation of mesh skin -----/
  rc = fc_getDecayedSkin(deadregion,NULL, &decayedSkin);
  fail_unless(rc == FC_SUCCESS, "failed to get decayedSkin");
  rc = fc_getSubsetMembersAsArray(decayedSkin,&i,&arr);
  fail_unless(i == 2, "wrong number of members in decayedskin");
  //do we want to check the vals ?
  fc_deleteSubset(decayedSkin);
  free(arr);

  // --- decay type tests ---/
  // this conceptually out of order but it uses the dead region
  rc = fc_getShapeSidesDecayType(deadregion,&shapes[0],&decayflagarray);
  fail_unless(rc == FC_SUCCESS, "failed to get decaytype");
  for (i = 0 ; i < shapes[0].numSides; i++){
    if (fc_isIntInSortedIntArray(&sia,i)){
      fail_unless(decayflagarray[i] == 1, "wrong decay type");
    }else{
      fail_unless(decayflagarray[i] == 0, "wrong decay type");
    }
  }
  free(decayflagarray);
  fc_freeSortedIntArray(&sia);

  rc = fc_getSubsetDecayType(emptye,shapes[0].elems[0],&decayflag);
  fail_unless(rc == FC_SUCCESS, "failed to get decaytype");
  fail_unless(decayflag == 0, "wrong val for decay type");

  rc = fc_getSubsetDecayType(emptyf,shapes[0].faces[0],&decayflag);
  fail_unless(rc == FC_SUCCESS, "failed to get decaytype");
  fail_unless(decayflag == 0, "wrong val for decay type");

  rc = fc_getShapeSidesDecayType(emptye,&shapes[0],&decayflagarray);
  fail_unless(rc == FC_SUCCESS, "failed to get decaytype");
  for (i = 0; i < shapes[0].numSides; i++){
    fail_unless(decayflagarray[i] == 0, "wrong val for decay type");
  }
  free(decayflagarray);

  rc = fc_getShapeSidesDecayType(emptyf,&shapes[0],&decayflagarray);
  fail_unless(rc == FC_SUCCESS, "failed to get decaytype");
  for (i = 0; i < shapes[0].numSides; i++){
    fail_unless(decayflagarray[i] == 0, "wrong val for decay type");
  }
  free(decayflagarray);

  rc = fc_getSubsetDecayType(whole,shapes[0].elems[0],&decayflag);
  fail_unless(rc == FC_SUCCESS, "failed to get decaytype");
  fail_unless(decayflag == 2, "wrong val for decay type");

  rc = fc_getSubsetDecayType(meshskin,shapes[0].faces[0],&decayflag);
  fail_unless(rc == FC_SUCCESS, "failed to get decaytype");
  fail_unless(decayflag == 2, "wrong val for decay type");

  rc = fc_getShapeSidesDecayType(whole,&shapes[0],&decayflagarray);
  fail_unless(rc == FC_SUCCESS, "failed to get decaytype");
  for (i = 0; i < shapes[0].numSides; i++){
    fail_unless(decayflagarray[i] == 2, "wrong val for decay type");
  }
  free(decayflagarray);

  rc = fc_getShapeSidesDecayType(meshskin,&shapes[0],&decayflagarray);
  fail_unless(rc == FC_SUCCESS, "failed to get decaytype");
  for (i = 0; i < shapes[0].numSides; i++){
    fail_unless(decayflagarray[i] == 2, "wrong val for decay type");
  }
  free(decayflagarray);

  rc = fc_getSubsetDecayType(shapes[0].elems[0],whole,&decayflag);
  fail_unless(rc == FC_SUCCESS, "failed to get decaytype");
  fail_unless(decayflag == 1, "wrong val for decay type");

  rc = fc_getSubsetDecayType(shapes[0].faces[0],meshskin,&decayflag);
  fail_unless(rc == FC_SUCCESS, "failed to get decaytype");
  fail_unless(decayflag == 1, "wrong val for decay type");

  rc = fc_getSubsetDecayType(shapes[0].elems[0],shapes[0].elems[0],&decayflag);
  fail_unless(rc == FC_SUCCESS, "failed to get decaytype");
  fail_unless(decayflag == 2, "wrong val for decay type");

  rc = fc_getSubsetDecayType(shapes[0].faces[0],shapes[0].faces[0],&decayflag);
  fail_unless(rc == FC_SUCCESS, "failed to get decaytype");
  fail_unless(decayflag == 2, "wrong val for decay type");

  rc = fc_getShapeSidesDecayType(shapes[0].elems[0],&shapes[0],
				 &decayflagarray);
  fail_unless(rc == FC_SUCCESS, "failed to get decaytype");
  fail_unless(decayflagarray[0] == 2, "wrong val for decay type");
  for (i = 1; i < shapes[0].numSides; i++){
    if (shapes[0].adjmatrix[0][i] == 1){
      //edge elements overlap adj sides of course, and we know its only
      //partial
      fail_unless(decayflagarray[i] == 1, "wrong val for decay type");
    }else{
      fail_unless(decayflagarray[i] == 0, "wrong val for decay type");
    }
  }
  free(decayflagarray);

  //bad args (other bad args are checked within the inner calls)
  rc = fc_getDecayedSkin(badSubset,&meshskin, &decayedSkin);
  fail_unless(rc != FC_SUCCESS, "should fail for bad subset");
  fail_unless(FC_HANDLE_EQUIV(decayedSkin,FC_NULL_SUBSET), "bad return val");

  rc = fc_getDecayedShapeSides(badSubset,&shapes[0], &decayedSides);
  fail_unless(rc != FC_SUCCESS, "should fail for bad subset");
  fail_unless(decayedSides == NULL, "bad return val");

  rc = fc_getDecayedShapeSkin(badSubset,&shapes[0], &decayedSkin);
  fail_unless(rc != FC_SUCCESS, "should fail for bad subset");
  fail_unless(decayedSides == NULL, "bad return val");

  rc = fc_getSubsetDecayType(badSubset,shapes[0].elems[0],&decayflag);
  fail_unless(rc != FC_SUCCESS, "should fail for bad subset");
  fail_unless(decayflag == -1, "wrong val for decay type");

  rc = fc_getSubsetDecayType(badSubset,shapes[0].faces[0],&decayflag);
  fail_unless(rc != FC_SUCCESS, "should fail for bad subset");
  fail_unless(decayflag == -1, "wrong val for decay type");

  rc = fc_getShapeSidesDecayType(badSubset,&shapes[0],
				 &decayflagarray);
  fail_unless(rc != FC_SUCCESS, "should fail for bad subset");
  fail_unless(decayflagarray == NULL, "bad return val");

  rc = fc_getDecayedSkin(FC_NULL_SUBSET,&meshskin, &decayedSkin);
  fail_unless(rc != FC_SUCCESS, "should fail for NULL subset");
  fail_unless(FC_HANDLE_EQUIV(decayedSkin,FC_NULL_SUBSET), "bad return val");

  rc = fc_getDecayedShapeSides(FC_NULL_SUBSET,&shapes[0], &decayedSides);
  fail_unless(rc != FC_SUCCESS, "should fail for NULL subset");
  fail_unless(decayedSides == NULL, "bad return val");

  rc = fc_getDecayedShapeSkin(FC_NULL_SUBSET,&shapes[0], &decayedSkin);
  fail_unless(rc != FC_SUCCESS, "should fail for NULL subset");
  fail_unless(decayedSides == NULL, "bad return val");

  rc = fc_getSubsetDecayType(FC_NULL_SUBSET,shapes[0].elems[0],&decayflag);
  fail_unless(rc != FC_SUCCESS, "should fail for bad subset");
  fail_unless(decayflag == -1, "wrong val for decay type");

  rc = fc_getSubsetDecayType(FC_NULL_SUBSET,shapes[0].faces[0],&decayflag);
  fail_unless(rc != FC_SUCCESS, "should fail for bad subset");
  fail_unless(decayflag == -1, "wrong val for decay type");

  rc = fc_getShapeSidesDecayType(FC_NULL_SUBSET,&shapes[0],
				 &decayflagarray);
  fail_unless(rc != FC_SUCCESS, "should fail for bad subset");
  fail_unless(decayflagarray == NULL, "bad return val");

  //mutate the shape - empty subset
  rc = fc_copyShape(&shapes[0],&badshape);
  fail_unless(rc == FC_SUCCESS, "cant copy shape");
  keep = badshape.faces[1];
  badshape.faces[1] = emptyf;

  rc = fc_getDecayedSkin(deadregion,&emptyf, &decayedSkin);
  fail_unless(rc != FC_SUCCESS, "should fail for empty skin");
  fail_unless(FC_HANDLE_EQUIV(decayedSkin,FC_NULL_SUBSET), "bad return val");

  rc = fc_getDecayedShapeSides(deadregion,&badshape, &decayedSides);
  fail_unless(rc != FC_SUCCESS, "should fail for empty subset in shape");
  fail_unless(decayedSides == NULL, "bad return val");

  rc = fc_getDecayedShapeSkin(deadregion,&badshape, &decayedSkin);
  fail_unless(rc != FC_SUCCESS, "should fail for empty subset in shape");
  fail_unless(decayedSides == NULL, "bad return val");

  rc = fc_getSubsetDecayType(shapes[0].elems[0],emptye,&decayflag);
  fail_unless(rc != FC_SUCCESS, "should fail for empty subset");
  fail_unless(decayflag == -1, "wrong val for decay type");

  rc = fc_getSubsetDecayType(shapes[0].faces[0],emptyf,&decayflag);
  fail_unless(rc != FC_SUCCESS, "should fail for empty subset");
  fail_unless(decayflag == -1, "wrong val for decay type");

  //test with the bad face set
  rc = fc_getShapeSidesDecayType(meshskin,&badshape,
				 &decayflagarray);
  fail_unless(rc != FC_SUCCESS, "should fail for empty subset in shape");
  fail_unless(decayflagarray == NULL, "bad return val");

  //get a bad elem set
  badshape.faces[1] = keep;
  keep = badshape.elems[1];
  badshape.elems[1] = emptye;

  rc = fc_getShapeSidesDecayType(deadregion,&badshape,
				 &decayflagarray);
  fail_unless(rc != FC_SUCCESS, "should fail for empty subset in shape");
  fail_unless(decayflagarray == NULL, "bad return val");

  //now put it back
  badshape.elems[1] = keep;
  keep = badshape.faces[1];
  badshape.faces[1] = emptye;


  //NULL subset
  rc = fc_getSubsetDecayType(shapes[0].elems[0],FC_NULL_SUBSET,&decayflag);
  fail_unless(rc != FC_SUCCESS, "should fail for NULL subset");
  fail_unless(decayflag == -1, "wrong val for decay type");

  rc = fc_getSubsetDecayType(shapes[0].faces[0],FC_NULL_SUBSET,&decayflag);
  fail_unless(rc != FC_SUCCESS, "should fail for NULL subset");
  fail_unless(decayflag == -1, "wrong val for decay type");

  rc = fc_getDecayedShapeSides(deadregion,NULL, &decayedSides);
  fail_unless(rc != FC_SUCCESS, "should fail for NULL shape");
  fail_unless(decayedSides == NULL, "bad return val");

  rc = fc_getDecayedShapeSkin(deadregion,NULL, &decayedSkin);
  fail_unless(rc != FC_SUCCESS, "should fail for NULL shape");
  fail_unless(decayedSides == NULL, "bad return val");

  rc = fc_getShapeSidesDecayType(deadregion,NULL,
				 &decayflagarray);
  fail_unless(rc != FC_SUCCESS, "should fail for NULL shape");
  fail_unless(decayflagarray == NULL, "bad return val");

  //mutate the shape - bad subset
  badshape.faces[1] = badSubset;

  rc = fc_getDecayedSkin(deadregion,&badSubset, &decayedSkin);
  fail_unless(rc != FC_SUCCESS, "should fail for badSubset");
  fail_unless(FC_HANDLE_EQUIV(decayedSkin,FC_NULL_SUBSET), "bad return val");

  rc = fc_getDecayedShapeSides(deadregion,&badshape, &decayedSides);
  fail_unless(rc != FC_SUCCESS, "should fail for bad subset in shape");
  fail_unless(decayedSides == NULL, "bad return val");

  rc = fc_getDecayedShapeSkin(deadregion,&badshape, &decayedSkin);
  fail_unless(rc != FC_SUCCESS, "should fail for bad subset in shape");
  fail_unless(decayedSides == NULL, "bad return val");

  rc = fc_getSubsetDecayType(shapes[0].elems[0],badSubset,&decayflag);
  fail_unless(rc != FC_SUCCESS, "should fail for bad subset");
  fail_unless(decayflag == -1, "wrong val for decay type");

  rc = fc_getSubsetDecayType(shapes[0].faces[0],badSubset,&decayflag);
  fail_unless(rc != FC_SUCCESS, "should fail for bad subset");
  fail_unless(decayflag == -1, "wrong val for decay type");

  //bad faces
  rc = fc_getShapeSidesDecayType(meshskin,&badshape,
				 &decayflagarray);
  fail_unless(rc != FC_SUCCESS, "should fail for bad subset in shape");
  fail_unless(decayflagarray == NULL, "bad return val");

  //get a bad elem set
  badshape.faces[1] = keep;
  keep = badshape.elems[1];
  badshape.elems[1] = badSubset;

  rc = fc_getShapeSidesDecayType(deadregion,&badshape,
				 &decayflagarray);
  fail_unless(rc != FC_SUCCESS, "should fail for bad subset in shape");
  fail_unless(decayflagarray == NULL, "bad return val");

  //now put it back
  badshape.elems[1] = keep;
  keep = badshape.faces[1];
  badshape.faces[1] = badSubset;

  //wrong, mismatches assoc
  //could have a few more of these with difff assoc
  rc = fc_getDecayedSkin(deadregion,&deadregion, &decayedSkin);
  fail_unless(rc != FC_SUCCESS, "should fail for wrong assoc");
  fail_unless(FC_HANDLE_EQUIV(decayedSkin,FC_NULL_SUBSET), "bad return val");

  rc = fc_getDecayedSkin(shapes[0].faces[0],&(shapes[0].faces[0]),
			 &decayedSkin);
  fail_unless(rc != FC_SUCCESS, "should fail for wrong assoc");
  fail_unless(FC_HANDLE_EQUIV(decayedSkin,FC_NULL_SUBSET), "bad return val");

  rc = fc_getDecayedShapeSides(shapes[0].faces[0],&shapes[0], &decayedSides);
  fail_unless(rc != FC_SUCCESS, "should fail for wrong assoc");
  fail_unless(decayedSides == NULL, "bad return val");

  rc = fc_getDecayedShapeSkin(shapes[0].faces[0],&badshape, &decayedSkin);
  fail_unless(rc != FC_SUCCESS, "should fail for wrong assoc");
  fail_unless(decayedSides == NULL, "bad return val");

  rc = fc_getSubsetDecayType(shapes[0].elems[0],shapes[0].faces[0],&decayflag);
  fail_unless(rc != FC_SUCCESS, "should fail for mismatched assoc");
  fail_unless(decayflag == -1, "wrong val for decay type");

  rc = fc_getSubsetDecayType(shapes[0].faces[0],shapes[0].elems[0],&decayflag);
  fail_unless(rc != FC_SUCCESS, "should fail for mismatched assoc");
  fail_unless(decayflag == -1, "wrong val for decay type");

  //NULL  return arg
  rc = fc_getDecayedSkin(deadregion,&(shapes[0].faces[0]), NULL);
  fail_unless(rc != FC_SUCCESS, "should fail for NULL return arg");

  rc = fc_getDecayedShapeSides(deadregion,&shapes[0], NULL);
  fail_unless(rc != FC_SUCCESS, "should fail for NULL return arg");

  rc = fc_getDecayedShapeSkin(deadregion,&shapes[0], NULL);
  fail_unless(rc != FC_SUCCESS, "should fail for NULL return arg");

  rc = fc_getSubsetDecayType(shapes[0].elems[0],shapes[0].elems[0],NULL);
  fail_unless(rc != FC_SUCCESS, "should fail for NULL return arg");

  rc = fc_getSubsetDecayType(shapes[0].faces[0],shapes[0].faces[0],NULL);
  fail_unless(rc != FC_SUCCESS, "should fail for NULL return arg");

  rc = fc_getShapeSidesDecayType(deadregion,&shapes[0],
				 NULL);
  fail_unless(rc != FC_SUCCESS, "should fail for NULL return arg");

  //clean up
  rc = fc_deleteSubset(emptye); 
  rc = fc_deleteSubset(emptyf); 
  rc = fc_deleteSubset(deadregion);
  rc = fc_deleteSubset(meshskin);
  rc = fc_deleteSubset(whole);

  for (i = 0; i < numShapes; i++){
    rc = fc_freeShape(&shapes[i]);
  }
  free(shapes);
  badshape.faces[1] = keep;
  rc = fc_freeShape(&badshape);
  if (rc != FC_SUCCESS){
    fc_printfErrorMessage("cant free bad shape");
  }

  fc_deleteMesh(mesh);
  fc_deleteDataset(dataset);
}
END_TEST
  
// Populate the Suite with the tests

Suite *elemdeath_suite(void)
{
  Suite *suite = suite_create("ElemDeath");

  TCase *tc_exposed = tcase_create(" - Exposed ");
  TCase *tc_tears = tcase_create(" - Tears ");
  TCase *tc_breaks = tcase_create(" - Breaks ");
  TCase *tc_decay = tcase_create(" - Decay ");

  suite_add_tcase(suite, tc_exposed);
  tcase_add_checked_fixture(tc_exposed, elemdeath_setup, elemdeath_teardown);
  tcase_add_test(tc_exposed, exposed_skin);

  suite_add_tcase(suite, tc_tears);
  tcase_add_checked_fixture(tc_tears, elemdeath_setup, elemdeath_teardown);
  tcase_add_test(tc_tears, tear_length);

  suite_add_tcase(suite, tc_breaks);
  tcase_add_checked_fixture(tc_breaks, elemdeath_setup, elemdeath_teardown);
  tcase_add_test(tc_breaks, breaks);

  suite_add_tcase(suite, tc_decay);
  tcase_add_checked_fixture(tc_decay, elemdeath_setup, elemdeath_teardown);
  tcase_add_test(tc_decay, decay);

  return suite;
}
