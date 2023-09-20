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
 * \file checktopo.c
 * \brief Unit testing of \ref TopologyRelations module
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/unittest/checktopo.c,v $
 * $Revision: 1.37 $ 
 * $Date: 2006/10/19 03:14:53 $
 *
 * \description
 *
 *    Tests the topo routines. It assumes that the internal mesh structs
 *    (_FC_Vertex, _FC_Edge, _FC_Face, _FC_Element) are correct (should
 *    have been tested in checkmesh). 
 *
 * \note 
 *    Tests are not currently being performed on the tet and pyramid
 *    meshes in gen_multimesh to save on time.
 *
 * \modifications
 *    - 07/11/04 WSK
 */

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <check.h>
#include "fc.h"
#include "fcP.h"
#include "checkall.h"

// **** global variables

static int numAssoc = 5;
static FC_AssociationType assocs[5] = { FC_AT_VERTEX, 
					FC_AT_EDGE, 
					FC_AT_FACE,
					FC_AT_ELEMENT, 
					FC_AT_WHOLE_MESH };

// **** test fixtures

static void topo_setup(void) {
  FC_ReturnCode rc;
  if (isForking) {
    fc_setLibraryVerbosity(fc_messages);
    rc = fc_initLibrary();
    fail_unless(rc == FC_SUCCESS, "test aborted: failed to init library");
  }
}

static void topo_teardown(void) {
  FC_ReturnCode rc;
  if (isForking) {
    rc = fc_finalLibrary();
    fail_unless(rc == FC_SUCCESS, "failed to final library at end of test");
  }
}

// FIX: This is not that rigorous of a test, mostly checks internal
// consistency and not absolute rightness.
START_TEST(get_children)
{
  FC_ReturnCode rc;
  int i, j, k, m, p;
  FC_Dataset dataset;
  FC_Mesh mesh, *meshes;
  int numMesh, i_mesh;
  int nums[5] = { 0, 0, 0, 0, 1 }, numChild, *childIDs;
  int numChild_temp;
  int numOwnerVertex, numChildVertex;
  int *owner_verts, *child_verts;
  FC_SortedIntArray sia;
  FC_ElementType elemType;
  
  // setup
  fc_loadDataset("../data/gen_multimesh.ex2", &dataset);
  fc_getMeshes(dataset, &numMesh, &meshes);
  
  // test all meshes
  for (i_mesh = 0; i_mesh < numMesh; i_mesh++) {
    int numVertPerElem, maxNumVertPerFace, *numVertPerFace;
    int *edgeConns, *faceConns, *elemConns;

    mesh = meshes[i_mesh];
    fc_getMeshElementType(mesh, &elemType);
    numVertPerElem = fc_getElementTypeNumVertex(elemType);
    fc_getMeshEdgeConnsPtr(mesh, &edgeConns);
    fc_getMeshFaceConnsPtr(mesh, &numVertPerFace, &maxNumVertPerFace,
			   &faceConns);
    fc_getMeshElementConnsPtr(mesh, &elemConns);

    // skip the tet mesh & the pyramid hex to save time
    if (i_mesh == 4 || i_mesh == 6)
      continue;

    //printf("i_mesh = %d, %s\n", i_mesh, meshnames[i_mesh]);

    for (i = 0; i < numAssoc; i++)
      fc_getMeshNumEntity(mesh, assocs[i], &nums[i]);
    //printf("num vertex = %d, numEdge = %d, numFace = %d, numElement = %d\n",
    // nums[0], nums[1], nums[2], nums[3]);
    
    // Expected error - can't call with FC_AT_WHOLE_DATASET
    rc = fc_getMeshEntityChildren(mesh, FC_AT_WHOLE_DATASET, 0, FC_AT_VERTEX,
				  &numChild, &childIDs);
    fail_unless(rc != FC_SUCCESS, "should fail for parent WHOLE_DATASET");
    fail_unless(numChild == -1 && childIDs == NULL, "fail should return NULL");

    // test
    for (i = 0; i < numAssoc; i++) {
      for (j = 0; j < nums[i]; j++) {

	// make owner verts
	fc_initSortedIntArray(&sia);
	switch (assocs[i]) {
	case FC_AT_EDGE:
	  numOwnerVertex = 2;
	  fc_addIntToSortedIntArray(&sia, edgeConns[j*2]);
	  fc_addIntToSortedIntArray(&sia, edgeConns[j*2+1]);
	  break;
	case FC_AT_FACE:
	  numOwnerVertex = numVertPerFace[j];
	  for (m = 0; m < numVertPerFace[j]; m++)
	    fc_addIntToSortedIntArray(&sia, faceConns[j*maxNumVertPerFace+m]);
	  break;
	case FC_AT_ELEMENT:
	  numOwnerVertex = numVertPerElem;
	  for (m = 0; m < numVertPerElem; m++)
	    fc_addIntToSortedIntArray(&sia, elemConns[j*numVertPerElem+m]);
	  break;
	default:
	  fail_unless(1, "should never reach here");
	}
	fc_convertSortedIntArrayToIntArray(&sia, &numOwnerVertex, &owner_verts);

	for (k = 0; k < numAssoc; k++) {
	  // for speed up, only test 1st & last 250 & middle 100 of whatever
	  //if ( (j > 250 && j < nums[i]/2 - 50) ||
	  //   (j > nums[i]/2+50 && j < nums[i]-250) )
	  //continue;

	  // get childs
	  rc = fc_getMeshEntityChildren(mesh, assocs[i], j, assocs[k], 
					&numChild, &childIDs);
	  if (i <= k) {
	    fail_unless(rc != FC_SUCCESS,
			"should fail if child assoc <= owner assoc");
	    fail_unless(numChild == -1, "fail should return null");
	    fail_unless(childIDs == NULL, "fail should return null");
	  }
	  else if (assocs[i] == FC_AT_WHOLE_MESH) {
	    fail_unless(rc == FC_SUCCESS, "shouldn't fail");
	    fail_unless(numChild == nums[k], "mismatch of numChild");
	    if (numChild > 0) {
	      for (m = 0; m < nums[k]; m++)
		fail_unless(childIDs[m] == m, "mismatch of childIDs");
	    }
	    else 
	      fail_unless(childIDs == NULL, "mismatch of childIDs");
	  }
	  else {
	    fail_unless(rc == FC_SUCCESS, "shouldn't fail");

	    // check that numChild is correct
	    fc_getMeshEntityElementType(mesh, assocs[i], j, &elemType);
	    switch(assocs[k]) {
	    case FC_AT_VERTEX: 
	      numChild_temp = fc_getElementTypeNumVertex(elemType);
	      break;
	    case FC_AT_EDGE:
	      numChild_temp = fc_getElementTypeNumEdge(elemType);
	      break;
	    case FC_AT_FACE:
	      numChild_temp = fc_getElementTypeNumFace(elemType);
	      break;
	    default:
	      fail_unless(1, "should never reach here");
	    }
	    fail_unless(numChild == numChild_temp, "mismatch of numchild");

	    // special case if numChild is zero
	    if (numChild == 0) {
	      fail_unless(childIDs == NULL, "mismatch of child IDs");
	      continue; // don't do following
	    }

	    // make composite list of child verts
	    fc_initSortedIntArray(&sia);
	    switch(assocs[k]) {
	    case FC_AT_VERTEX: 
	      numChild_temp = fc_getElementTypeNumVertex(elemType);
	      for (m = 0; m < numChild; m++)
		fc_addIntToSortedIntArray(&sia, childIDs[m]);
	      break;
	    case FC_AT_EDGE:
	      numChild_temp = fc_getElementTypeNumEdge(elemType);
	      for (m = 0; m < numChild; m++) {
		fc_addIntToSortedIntArray(&sia, edgeConns[childIDs[m]*2]);
		fc_addIntToSortedIntArray(&sia, edgeConns[childIDs[m]*2+1]);
	      }
	      break;
	    case FC_AT_FACE:
	      numChild_temp = fc_getElementTypeNumFace(elemType);
	      for (m = 0; m < numChild; m++) {
		for (p = 0; p < numVertPerFace[childIDs[m]]; p++) {
		  fc_addIntToSortedIntArray(&sia,
				 faceConns[childIDs[m]*maxNumVertPerFace+p]);
		}
	      }
	      break;
	    default:
	      fail_unless(1, "should never reach here");
	    }
	    fc_convertSortedIntArrayToIntArray(&sia, &numChildVertex, 
					       &child_verts);
	    fail_unless(numChildVertex == numOwnerVertex, 
			"mismatch in number of vertices");
	    fail_unless(!memcmp(child_verts, owner_verts,
	       		numChildVertex*sizeof(int)), "mismatch of vertIDs");
	    free(child_verts);
	  }
	  free(childIDs);
	}// k loop over assoc
	free(owner_verts);
      } // j loop over num
    } // i loop over assoc
    fail_unless(i == numAssoc, "test error, didn't finish looping");

    fc_deleteMesh(mesh);
  } // loop over i_mesh

  free(meshes);
  fc_deleteDataset(dataset);
}
END_TEST

START_TEST(get_parents)
{
  FC_ReturnCode rc;
  int i, j, k, m;
  FC_Dataset dataset;
  FC_Mesh mesh, *meshes;
  int numMesh, i_mesh;
  int nums[5] = { 0, 0, 0, 0, 1 }, numParent, *parentIDs;
  int numParent_temp, *parentIDs_temp;
  FC_SortedIntArray* lists;
  int numChild, *childIDs;

  // setup
  fc_loadDataset("../data/gen_multimesh.ex2", &dataset);
  fc_getMeshes(dataset, &numMesh, &meshes);

  // test all meshes
  for (i_mesh = 0; i_mesh < numMesh; i_mesh++) {
    mesh = meshes[i_mesh];

    // skip the tet mesh & the pyramid hex to save time
    if (i_mesh == 4 || i_mesh == 6)
      continue;

    //printf("i_mesh = %d, %s\n", i_mesh, meshnames[i_mesh]);

    // setup
    for (i = 0; i < numAssoc; i++)
      fc_getMeshNumEntity(mesh, assocs[i], &nums[i]);
    //printf("num vertex = %d, numEdge = %d, numFace = %d, numElement = %d\n",
    //   nums[0], nums[1], nums[2], nums[3]);

    // Expected error - can't call with FC_AT_WHOLE_DATASET
    rc = fc_getMeshEntityParents(mesh, FC_AT_VERTEX, 0, FC_AT_WHOLE_DATASET,
				  &numParent, &parentIDs);
    fail_unless(rc != FC_SUCCESS, "should fail for parent WHOLE_DATASET");
    fail_unless(numParent == -1 && parentIDs == NULL, "fail should return NULL");

    // test
    for (i = 0; i < numAssoc; i++) {
      for (j = 0; j < numAssoc; j++) {

	// make adjacencies
	lists = NULL;
	if (i < 3 && j > i && j < 4) {
	  lists = (FC_SortedIntArray*)calloc(nums[i], 
					     sizeof(FC_SortedIntArray));
	  for (k = 0; k < nums[j]; k++) {
	    fc_getMeshEntityChildren(mesh, assocs[j], k, assocs[i], &numChild,
				     &childIDs);
	    for (m = 0; m < numChild; m++)
	      fc_addIntToSortedIntArray(&lists[childIDs[m]], k);
	    free(childIDs);
	  }
	}
 
	for (k = 0; k < nums[i]; k++) {
	  // for speed up, only test 1st & last 250 & middle 100 of whatever
	  //if ( (k > 250 && k < nums[i]/2 - 50) ||
	  // (k > nums[i]/2+50 && k < nums[i]-250) )
	  //continue;

	  // get parents
	  rc = fc_getMeshEntityParents(mesh, assocs[i], k, assocs[j], 
				       &numParent, &parentIDs);
	  if (j <= i) {
	    fail_unless(rc != FC_SUCCESS, 
			"should fail if parent assoc <= owned assoc");
	    fail_unless(numParent == -1, "fail should return null");
	    fail_unless(parentIDs == NULL, "fail should return null");
	  }
	  else if (assocs[j] == FC_AT_WHOLE_MESH) {
	    fail_unless(rc == FC_SUCCESS, "shouldn't fail");
	    fail_unless(numParent == 1, "mismatch of numParent");
	    fail_unless(parentIDs[0] == 0, "mismatch of parentIDs");
	  }
	  else {
	    fail_unless(rc == FC_SUCCESS, "shouldn't fail");
	    fc_convertSortedIntArrayToIntArray(&lists[k], &numParent_temp, &parentIDs_temp);
	    
	    fail_unless(numParent == numParent_temp, "mismatch of numparent");
	    fail_unless(!memcmp(parentIDs, parentIDs_temp, 
			     numParent*sizeof(int)), "mismatch of parentIDs");
	    free(parentIDs_temp);
	  }
	  free(parentIDs);
	  
	}// k loop over num

	free(lists);

      } // j loop over assoc
    } // i loop over assoc
    fc_deleteMesh(mesh);
  } // loop over i_mesh

  free(meshes);
  fc_deleteDataset(dataset);
}
END_TEST

// Only testing a representative case
// Start with the verts, edges, faces or element that make up a hex
// For lenient testing, all hexes that one of the starting entites will
// be in the new list. (They will come out of the 27 hex block centered
// on the hex.
// For strict testing, only the entities within the hex will be in new list. 
START_TEST(change_entity_type)
{
  FC_ReturnCode rc;
  int i, j, k, m;
  FC_Dataset dataset;
  FC_Mesh mesh, *returnMeshes;
  int numReturnMeshes;
  int numEntity;
  int numTempID, *tempIDs;
  int numOldID, *oldIDs, numNewID, *newIDs, numNewID_good, *newIDs_good;
  FC_SortedIntArray sia;
  int hexID = 145;  // this hex must lie in interior
  // number of expected new type for lenient case, indexed by old then new
  // the diagonal is number of that entity in a hex
  int nums[4][4] = { { 8, 36, 54, 27 }, { 8, 12, 30, 19 },
		     { 8, 12,  6, 7 }, { 8, 12,  6,  1 } };
  // note: 36 = 12 + 3*8, 54 = 3*27 - 9*3, 30 = ?, 19 = 27-8, 7 = 27-8-12

  // setup
  fc_loadDataset("../data/gen_multimesh.ex2", &dataset);
  rc = fc_getMeshByName(dataset, "hex mesh", &numReturnMeshes,&returnMeshes);
  fail_unless(rc == FC_SUCCESS, "aborted - cant get mesh by name");
  fail_unless(numReturnMeshes == 1, "failed to find unique mesh by name");
  mesh = returnMeshes[0];
  free(returnMeshes);

  // Loop over all combinations of association (except neglect FC_AT_WHOLE_MESH)
  for (i = 0; i < numAssoc-1; i++) { // old
    for (j = 0; j < numAssoc-1; j++) { // new

      // some setup
      if (assocs[i] < FC_AT_ELEMENT) {
	fc_getMeshEntityChildren(mesh, FC_AT_ELEMENT, hexID, assocs[i],
				 &numOldID, &oldIDs);
	fc_sortIntArray(numOldID, oldIDs);
      }
      else {
	numOldID = 1;
	oldIDs = malloc(sizeof(int));
	oldIDs[0] = hexID;
      }
 
      // create expected answers
      if (assocs[j] == assocs[i]) {
	numNewID_good = numOldID;
	newIDs_good = malloc(numOldID*sizeof(int));
	memcpy(newIDs_good, oldIDs, numOldID*sizeof(int));
      }
      else { 
	fc_initSortedIntArray(&sia);
	for (k = 0; k < numOldID; k++) {
	  if (assocs[j] < assocs[i]) 
	    fc_getMeshEntityChildren(mesh, assocs[i], oldIDs[k], assocs[j],
				     &numTempID, &tempIDs);
	  else 
	    fc_getMeshEntityParents(mesh, assocs[i], oldIDs[k], assocs[j],
				    &numTempID, &tempIDs);
	  for (m = 0; m < numTempID; m++)
	    fc_addIntToSortedIntArray(&sia, tempIDs[m]);
	  free(tempIDs);
	}
	fc_convertSortedIntArrayToIntArray(&sia, &numNewID_good, &newIDs_good);
      }
      fail_unless(numNewID_good == nums[i][j], "test aborted: error in calc"); 

      // do it - "lenient"
      rc = fc_changeMeshEntityType(mesh, assocs[i], numOldID, oldIDs,
				   assocs[j], 0, &numNewID, &newIDs);
      fail_unless(rc == FC_SUCCESS, "failed to get different types");
      fail_unless(numNewID == numNewID_good, "mismatch of numNewID");
      fail_unless(!memcmp(newIDs, newIDs_good, numNewID*sizeof(int)),
		  "mismatch of IDs");
      free(newIDs);
      
      // do it - "strict"
      rc = fc_changeMeshEntityType(mesh, assocs[i], numOldID, oldIDs,
				   assocs[j], 1, &numNewID, &newIDs);
      fail_unless(rc == FC_SUCCESS, "failed to get different types");
      //printf("strict: numNewID = %d, numTempID = %d\n", numNewID, numTempID);
      if (j > i) {
	fail_unless(numNewID == nums[j][j], "mismatch of numNewID");
	if (j == 3)
	  fail_unless(numNewID == 1 && newIDs[0] == hexID, "mismatch of IDs");
	else {
	  fc_getMeshEntityChildren(mesh, FC_AT_ELEMENT, hexID, assocs[j],
				  &numTempID, &tempIDs);
	  fc_sortIntArray(numTempID, tempIDs);
	  fail_unless(numTempID == numNewID && 
		      !memcmp(newIDs, tempIDs, numNewID*sizeof(int)),
		      "mismatch of IDs");
	  free(tempIDs);
	}
      }
      else {
	fail_unless(numNewID == numNewID_good, "mismatch of numNewID");
	fail_unless(!memcmp(newIDs, newIDs_good, numNewID*sizeof(int)),
		    "mismatch of IDs");      
      }
      free(newIDs);

      // special case - empty old list is o.k.
      for (k = 0; k < 2; k++) {
	rc = fc_changeMeshEntityType(mesh, assocs[i], 0, NULL, assocs[j], k,
				     &numNewID, &newIDs);
	fail_unless(rc == FC_SUCCESS, "should not fail for 0 oldIDs");
	fail_unless(numNewID == 0 && newIDs == NULL, 
		    "mismatch of numNewID and newIDs"); 
      }

      // test bad arguments
      rc = fc_changeMeshEntityType(FC_NULL_MESH, assocs[i], numOldID, oldIDs,
				   assocs[j], 0, &numNewID, &newIDs);
      fail_unless(rc != FC_SUCCESS, "should fail with bad mesh handle");
      fail_unless(numNewID == -1 && newIDs == NULL, 
		  "fail should return nulls"); 
      rc = fc_changeMeshEntityType(mesh, -2, numOldID, oldIDs,
				   assocs[j], 0, &numNewID, &newIDs);
      fail_unless(rc != FC_SUCCESS, "should fail with bad old assoc type");
      fail_unless(numNewID == -1 && newIDs == NULL, 
		  "fail should return nulls"); 
      rc = fc_changeMeshEntityType(mesh, assocs[i], -1, oldIDs,
				   assocs[j], 0, &numNewID, &newIDs);
      fail_unless(rc != FC_SUCCESS, "should fail if zero entities");
      fail_unless(numNewID == -1 && newIDs == NULL, 
		  "fail should return nulls"); 
      rc = fc_changeMeshEntityType(mesh, assocs[i], numOldID, NULL,
				   assocs[j], 0, &numNewID, &newIDs);
      fail_unless(rc != FC_SUCCESS, "should fail if no entites given");
      fail_unless(numNewID == -1 && newIDs == NULL, 
		  "fail should return nulls"); 
      rc = fc_changeMeshEntityType(mesh, assocs[i], numOldID, oldIDs,
				   -2, 0, &numNewID, &newIDs);
      fail_unless(rc != FC_SUCCESS, "should fail if bad new assoc type");
      fail_unless(numNewID == -1 && newIDs == NULL, 
		  "fail should return nulls"); 
      rc = fc_changeMeshEntityType(mesh, assocs[i], numOldID, oldIDs,
				   assocs[j], -1, &numNewID, &newIDs);
      fail_unless(rc != FC_SUCCESS, "should fail if bad doStrict flag");
      fail_unless(numNewID == -1 && newIDs == NULL, 
		  "fail should return nulls"); 
      rc = fc_changeMeshEntityType(mesh, assocs[i], numOldID, oldIDs,
				   assocs[j], 2, &numNewID, &newIDs);
      fail_unless(rc != FC_SUCCESS, "should fail if bad doStrict flag");
      fail_unless(numNewID == -1 && newIDs == NULL, 
		  "fail should return nulls"); 
      rc = fc_changeMeshEntityType(mesh, assocs[i], numOldID, oldIDs,
				   assocs[j], 0, NULL, &newIDs);
      fail_unless(rc != FC_SUCCESS, "should fail if null numNewID");
      fail_unless(numNewID == -1 && newIDs == NULL, 
		  "fail should return nulls"); 
      rc = fc_changeMeshEntityType(mesh, assocs[i], numOldID, oldIDs,
				   assocs[j], 0, &numNewID, NULL);
      fail_unless(rc != FC_SUCCESS, "should fail if null newIDs");
      fail_unless(numNewID == -1 && newIDs == NULL, 
		  "fail should return nulls"); 

      // cleanup
      free(oldIDs);
      free(newIDs_good);
    }
  }

  // Special tests for FC_AT_WHOLE_MESH
  for (i = 0; i < numAssoc - 1; i++) { 
    // setup
    fc_getMeshNumEntity(mesh, assocs[i], &numEntity);

    // Start with FC_AT_WHOLE_MESH = get all entities of that type
    numOldID = 1;
    oldIDs = malloc(sizeof(int));
    oldIDs[0] = 0;
    for (j = 0; j < 2; j++) { // both lenient and strict have same answer
      rc = fc_changeMeshEntityType(mesh, FC_AT_WHOLE_MESH, numOldID, oldIDs,
				   assocs[i], j, &numNewID, &newIDs);
      fail_unless(rc == FC_SUCCESS, "failed to get different types");
      fail_unless(numNewID == numEntity, "mismatch of numNewID");
      for (k = 0; k < numNewID; k++)
	fail_unless(newIDs[k] == k, "mismatch of IDs");
      free(newIDs);
    }
    free(oldIDs);

    // End with FC_AT_WHOLE_MESH - if lenient, get whole, if strict, get nothing
    if (i == 3) {
      numOldID = 1;
      oldIDs = malloc(sizeof(int));
      oldIDs[0] = hexID;
    }
    else 
      fc_getMeshEntityChildren(mesh, FC_AT_ELEMENT, hexID, assocs[i],
			       &numOldID, &oldIDs);
    rc = fc_changeMeshEntityType(mesh, assocs[i], numOldID, oldIDs,
				   FC_AT_WHOLE_MESH, 0, &numNewID, &newIDs);
    fail_unless(rc == FC_SUCCESS, "failed to get different types");
    fail_unless(numNewID == 1, "mismatch of numNewID");
    fail_unless(newIDs[0] == 0, "mismatch of IDs");
    free(newIDs);
    rc = fc_changeMeshEntityType(mesh, assocs[i], numOldID, oldIDs,
				   FC_AT_WHOLE_MESH, 1, &numNewID, &newIDs);
    fail_unless(rc == FC_SUCCESS, "failed to get different types");
    fail_unless(numNewID == 0, "mismatch of numNewID");
    fail_unless(newIDs == NULL, "mismatch of IDs");
    free(oldIDs);
  }

  // Expect failures for FC_AT_WHOLE_DATASET
  numOldID = 1;
  oldIDs = malloc(sizeof(int));
  oldIDs[0] = 0;
  rc = fc_changeMeshEntityType(mesh, FC_AT_WHOLE_DATASET, numOldID, oldIDs,
			       FC_AT_VERTEX, 0, &numNewID, &newIDs);
  fail_unless(rc != FC_SUCCESS, "should fail for FC_AT_WHOLE_DATASET");
  fail_unless(numNewID == -1 && newIDs == NULL, "fail should return NULL");
  rc = fc_changeMeshEntityType(mesh, FC_AT_VERTEX, numOldID, oldIDs,
			       FC_AT_WHOLE_DATASET, 0, &numNewID, &newIDs);
  fail_unless(rc != FC_SUCCESS, "should fail for FC_AT_WHOLE_DATASET");
  fail_unless(numNewID == -1 && newIDs == NULL, "fail should return NULL");
  free(oldIDs);

  // cleanup
  fc_deleteMesh(mesh);
  fc_deleteDataset(dataset);
}
END_TEST
  
START_TEST(entity_neighbors)
{
  FC_ReturnCode rc;
  int i, j, k, m, dim_i;
  FC_Dataset dataset;
  FC_Mesh mesh, *meshes;
  int numMesh, i_mesh;
  int nums[5] = { 0, 0, 0, 0, 1 }, numNeigh, *neighIDs;
  int numNeigh_temp, *neighIDs_temp;

  // setup
  fc_loadDataset("../data/gen_multimesh.ex2", &dataset);
  fc_getMeshes(dataset, &numMesh, &meshes);

  // test all meshes
  for (i_mesh = 0; i_mesh < numMesh; i_mesh++) {
    int *edgeConns;
    mesh = meshes[i_mesh];
    
    // skip the tet mesh & the pyramid hex to save time
    if (i_mesh == 4 || i_mesh == 6)
      continue;
    
    //printf("i_mesh = %d, %s\n", i_mesh, meshnames[i_mesh]);

    // setup
    for (i = 0; i < numAssoc; i++)
      fc_getMeshNumEntity(mesh, assocs[i], &nums[i]);
    fc_getMeshEdgeConnsPtr(mesh, &edgeConns);
    //printf("num vertex = %d, numEdge = %d, numFace = %d, numElement = %d\n",
    //   nums[0], nums[1], nums[2], nums[3]);

    // Expected error - can't call with FC_AT_WHOLE_DATASET
    rc = fc_getMeshEntityNeighbors(mesh, 0, FC_AT_WHOLE_DATASET, 0,
				  &numNeigh, &neighIDs);
    fail_unless(rc != FC_SUCCESS, "should fail for  WHOLE_DATASET");
    fail_unless(numNeigh == -1 && neighIDs == NULL, "fail should return NULL");

    // test
    for (i = 0; i < numAssoc; i++) {
      FC_SortedIntArray* neighLists = calloc(nums[i], sizeof(FC_SortedIntArray));
      
      for (dim_i = 0; dim_i < 3; dim_i++) {

	// create neighbor lists for of the assoc
	if (dim_i == 0 && i == 0) {
	  for (j = 0; j < nums[1]; j++) {
	    fc_addIntToSortedIntArray(&neighLists[edgeConns[j*2+1]],
				      edgeConns[j*2]);
	    fc_addIntToSortedIntArray(&neighLists[edgeConns[j*2]],
				      edgeConns[j*2+1]);
	  }
	}
	else if (dim_i < i) {
	  // WARNING: the following depends on enum values not changing (7/04)
	  FC_AssociationType childAssoc = dim_i + 1;
	  int numParent, *parentIDs;
	  for (j = 0; j < nums[dim_i]; j++) {
	    rc = fc_getMeshEntityParents(mesh, childAssoc, j, assocs[i],
					 &numParent, &parentIDs);
	    fail_unless(rc == FC_SUCCESS, "test aborted: bug in test code");
	    for (k = 0; k < numParent; k++) {
	      for (m = 0; m < numParent; m++) {
		if (m == k)
		  continue; // don't add to self
		fc_addIntToSortedIntArray(&neighLists[parentIDs[m]],
					  parentIDs[k]);
	      }
	    }
	    free(parentIDs);
	  }
	}

	// test each entity
	for (j = 0; j < nums[i]; j++) {
	  rc = fc_getMeshEntityNeighbors(mesh, j, assocs[i], dim_i,
					 &numNeigh, &neighIDs);
	  // test return values
	  if ( (i == 0 && dim_i > 0) || (i > 0 && dim_i >= i) ) {
	    fail_unless(rc != FC_SUCCESS, "should have failed");
	    fail_unless(numNeigh == -1, "fail should return -1 for num");
	    fail_unless(neighIDs == NULL, "fail should return NULL IDs");
	    continue;
	  }
	  else {
	    fail_unless(rc == FC_SUCCESS, "shouldn't fail");
	    fc_convertSortedIntArrayToIntArray(&neighLists[j], &numNeigh_temp,
				   &neighIDs_temp);
	    fail_unless(numNeigh == numNeigh_temp, "mismatch of numNeighbor");
	    fail_unless(!memcmp(neighIDs, neighIDs_temp, numNeigh*sizeof(int)),
			  "mismatch of neighborIDs");
	    free(neighIDs_temp);
	  }
	  free(neighIDs);

	} // looper over j (nums[i]) 
      } // loop over i_dim

      free(neighLists);

    } // i loop over assoc

    fc_deleteMesh(mesh);

  } // loop over i_mesh

  free(meshes);
  fc_deleteDataset(dataset);
}
END_TEST


// assuming previous tests passed, a few represntative cases are probably
// enough
START_TEST(subset_neighbors)
{
  FC_ReturnCode rc;
  int i, j, k, dim_i;
  FC_Dataset dataset;
  FC_Mesh mesh, *meshes;
  FC_Subset subset;
  FC_AssociationType neigh_assoc;
  int numMesh, i_mesh;
  int nums[5] = { 0, 0, 0, 0, 1 }, numNeigh, *neighIDs;

  // setup
  fc_loadDataset("../data/gen_multimesh.ex2", &dataset);
  fc_getMeshes(dataset, &numMesh, &meshes);

  // test all meshes
  for (i_mesh = 0; i_mesh < numMesh; i_mesh++) {
    mesh = meshes[i_mesh];

    // skip the tet mesh & the pyramid hex to save time
    //if (i_mesh == 4 || i_mesh == 6)
    //continue;

    //printf("i_mesh = %d, %s\n", i_mesh, meshnames[i_mesh]);

    // setup
    for (i = 0; i < numAssoc; i++)
      fc_getMeshNumEntity(mesh, assocs[i], &nums[i]);
    //printf("num vertex = %d, numEdge = %d, numFace = %d, numElement = %d\n",
    //   nums[0], nums[1], nums[2], nums[3]);

    // test all association types
    for (i = 0; i < numAssoc; i++) {
      int seedID;

      // skip non existent entity types
      if (nums[i] == 0)
	continue;

      // make a subset with just one entity, one in the middle
      seedID = nums[i]/2;
      fc_createSubset(mesh, "temp subset", assocs[i], &subset);
      fc_addMemberToSubset(subset, seedID);

      // test different shared_dims
      for (dim_i = 0; dim_i < 3; dim_i++) {
	FC_SortedIntArray *list, *newlist, *sia_p;

	list = malloc(sizeof(FC_SortedIntArray));
	newlist = malloc(sizeof(FC_SortedIntArray));
	fc_initSortedIntArray(list);
	fc_initSortedIntArray(newlist);

	// start up a list of neighbors (inclusive of seed)
	fc_addIntToSortedIntArray(list, seedID);

	// test different levels
	for (j = 1; j < 4; j++) {
	  rc = fc_getSubsetNeighbors(subset, j, dim_i, &neigh_assoc,
					 &numNeigh, &neighIDs);
	  // test return values
	  if ( j == 0 || (i == 0 && dim_i > 0) || (i > 0 && dim_i >= i) ) {
	    fail_unless(rc != FC_SUCCESS, "should have failed");
	    fail_unless(neigh_assoc == FC_AT_UNKNOWN, 
			"fail should return NULL assoc");
	    fail_unless(numNeigh == -1, "fail should return -1 for num");
	    fail_unless(neighIDs == NULL, "fail should return NULL IDs");
	  }
	  else {
	    fail_unless(rc == FC_SUCCESS, "shouldn't fail");
	    fail_unless(neigh_assoc == assocs[i], "mismatch of assoc");
	    // each time through j loop, add more neighbors
	    for (k = 0; k < list->numVal; k++) {
	      int key = list->vals[k];
	      int num, *IDs;
	      fc_getMeshEntityNeighbors(mesh, key, assocs[i], dim_i,
					&num, &IDs);
	      fc_addIntArrayToSortedIntArray(newlist, num, IDs, 0);
	      fc_addIntToSortedIntArray(newlist, key);
	      free(IDs);
	    }
	    // remove seed & convert to array
	    fc_deleteIntFromSortedIntArray(newlist, seedID);
	    fail_unless(numNeigh == newlist->numVal, "mismatch of numNeighbor");
	    fail_unless(!memcmp(neighIDs, newlist->vals, numNeigh*sizeof(int)),
			"mismatch of neighborIDs");
	    sia_p = list;
	    list = newlist;
	    newlist = sia_p;
	    fc_freeSortedIntArray(newlist);
	  }
	  free(neighIDs);
	} // loop over levels
	fc_freeSortedIntArray(list);
	free(list);
	free(newlist);
      } // loop over i_dim
      fc_deleteSubset(subset);
    } // i loop over assoc
    fc_deleteMesh(mesh);
  } // loop over i_mesh

  free(meshes);
  fc_deleteDataset(dataset);
}
END_TEST

// Since underlying routines tested above, will just test in hex mesh
// * Choose 2 regions of elements, 1 will always be 1 region but 2nd may be 
//     differently connected with different "minimum shared dims".
// * Tested different association types by just using the children of
//     elements.
// * Have to test FC_AT_WHOLE_MESH separately because it doesn't behave the same.
START_TEST(segmentation)
{
  FC_ReturnCode rc;
  int i, j, k, m, n;
  FC_Dataset dataset;
  FC_Mesh mesh, *returnMeshes;
  int numReturnMeshes;
  FC_Subset subset, *segments;
  int numSegment;
  int shared_dims[5] = { 0, 0, 1, 2, 2 };
  int numBlock1_hex = 12;
  // indexed by shared_dim, then the segment#
  int numBlock2_hex[3][3] = { { 3, 0, 0 }, { 2, 1, 0 }, { 1, 1, 1 } };
  int block1_hex[12] = { 0, 1, 2, 11, 12, 13, 132, 133, 134, 143, 144, 145 };
  int block2_hex[3] = { 17, 29, 159 };
  int *block2_hex_ptr;
  int numSegment_good[3] = { 2, 3, 4 };
  int numTemp, *temps;
  int numBlock, *block;
  FC_SortedIntArray block_list;

  // setup
  fc_loadDataset("../data/gen_multimesh.ex2", &dataset);
  rc = fc_getMeshByName(dataset, "hex mesh", &numReturnMeshes,&returnMeshes);
  fail_unless(rc == FC_SUCCESS, "aborted - cant get mesh by name");
  fail_unless(numReturnMeshes == 1, "failed to find unique mesh by name");
  mesh = returnMeshes[0];
  free(returnMeshes);
 
  // For each association type except FC_AT_WHOLE_MESH
  for (i = 0; i < numAssoc -1 ; i++) {

    // test empty subset
    rc = fc_createSubset(mesh, "temp subset", assocs[i], &subset);
    fail_unless(rc == FC_SUCCESS, "failed to create empty subset");
    for (j = 0; j <= shared_dims[i]; j++) {
      rc = fc_segment(subset, j, &numSegment, &segments);
      fail_unless(rc == FC_SUCCESS, "failed to segment empty subset");
      fail_unless(numSegment == 0, "empty subset should have no segments");
      fail_unless(segments == NULL, "empty subset should have null segements");
    }

    // populate subset with specific hex elements or their children
    if (assocs[i] == FC_AT_ELEMENT) {
      fc_addArrayMembersToSubset(subset, numBlock1_hex, block1_hex);
      fc_addArrayMembersToSubset(subset, numBlock2_hex[0][0], block2_hex);
    }
    else {
      for (j = 0; j < numBlock1_hex; j++) {
	fc_getMeshEntityChildren(mesh, FC_AT_ELEMENT, block1_hex[j],
				 assocs[i], &numTemp, &temps);
	fc_addArrayMembersToSubset(subset, numTemp, temps);
	free(temps);
      }
      for (j = 0; j < numBlock2_hex[0][0]; j++) {
	fc_getMeshEntityChildren(mesh, FC_AT_ELEMENT, block2_hex[j],
				 assocs[i], &numTemp, &temps);
	fc_addArrayMembersToSubset(subset, numTemp, temps);
	free(temps);
      }
    }

    // test for each possible shared_dim
    for (j = 0; j <= shared_dims[i]; j++) {
      // do the segmentation
      rc = fc_segment(subset, j, &numSegment, &segments);
      fail_unless(rc == FC_SUCCESS, "failed to segment subset");
      fail_unless(numSegment == numSegment_good[j], "mismatch of numSegment");

      // check block1 (will always be first segment)
      if (assocs[i] == FC_AT_ELEMENT) {
	numBlock = numBlock1_hex;
	block = block1_hex;
      }
      else {
	fc_initSortedIntArray(&block_list);
	for (k = 0; k < numBlock1_hex; k++) {
	  fc_getMeshEntityChildren(mesh, FC_AT_ELEMENT, block1_hex[k],
				   assocs[i], &numTemp, &temps);
	  for (m = 0; m < numTemp; m++)
	    fc_addIntToSortedIntArray(&block_list, temps[m]);
	  free(temps);
	}
	fc_convertSortedIntArrayToIntArray(&block_list, &numBlock, &block);
      }
      fc_getSubsetMembersAsArray(segments[0], &numTemp, &temps);
      fail_unless(numTemp == numBlock, "mismatch of numMember");
      fail_unless(!memcmp(temps, block, numBlock*sizeof(int)),
		  "mismatch of members");
      free(temps);
      if (assocs[i] != FC_AT_ELEMENT)
	free(block);

      // check block2 (number of segments depends on shared_dim)
      block2_hex_ptr = block2_hex;
      for (k = 0; k < j+1; k++) {
	if (assocs[i] == FC_AT_ELEMENT) {
	  numBlock = numBlock2_hex[j][k];
	  block = malloc(numBlock*sizeof(int));
	  memcpy(block, block2_hex_ptr, numBlock*sizeof(int));
	}
	else {
	  fc_initSortedIntArray(&block_list);
	  for (m = 0; m < numBlock2_hex[j][k]; m++) {
	    fc_getMeshEntityChildren(mesh, FC_AT_ELEMENT, block2_hex_ptr[m],
				     assocs[i], &numTemp, &temps);
	    for (n = 0; n < numTemp; n++)
	      fc_addIntToSortedIntArray(&block_list, temps[n]);
	    free(temps);
	  }
	  fc_convertSortedIntArrayToIntArray(&block_list, &numBlock, &block);
	}
	fc_getSubsetMembersAsArray(segments[1+k], &numTemp, &temps);
	fail_unless(numTemp == numBlock, "mismatch of numMember");
	fail_unless(!memcmp(temps, block, numBlock*sizeof(int)),
		    "mismatch of members");
	free(block);
	free(temps);
	block2_hex_ptr += numBlock2_hex[j][k];
      }

      // cleanup
      for (k = 0; k < numSegment; k++)
        fc_deleteSubset(segments[k]);
      free(segments);
    }

    // test bad arguments
    rc = fc_segment(FC_NULL_SUBSET, 0, &numSegment, &segments);
    fail_unless(rc != FC_SUCCESS, "should fail if bad subset");
    fail_unless(numSegment = -1 && segments == NULL, 
		"fail should return nulls"); 
    rc = fc_segment(subset, -1, &numSegment, &segments);
    fail_unless(rc != FC_SUCCESS, "should fail if bad shared_dim");
    fail_unless(numSegment = -1 && segments == NULL, 
		"fail should return nulls");
    rc = fc_segment(subset, shared_dims[i]+1, &numSegment, &segments);
    fail_unless(rc != FC_SUCCESS, "should fail if bad shared_dim");
    fail_unless(numSegment = -1 && segments == NULL, 
		"fail should return nulls"); 
    rc = fc_segment(subset, 0, NULL, &segments);
    fail_unless(rc != FC_SUCCESS, "should fail if null numSegment");
    fail_unless(segments == NULL, "fail should return nulls"); 
    rc = fc_segment(subset, shared_dims[i]+1, &numSegment, NULL);
    fail_unless(rc != FC_SUCCESS, "should fail if null segmetns");
    fail_unless(numSegment = -1, "fail should return nulls"); 

    // cleanup
    fc_deleteSubset(subset);
  }

  // test FC_AT_WHOLE_MESH 
  // test empty subset
  rc = fc_createSubset(mesh, "temp subset", FC_AT_WHOLE_MESH, &subset);
  fail_unless(rc == FC_SUCCESS, "failed to create empty subset");
  for (j = 0; j <= shared_dims[i]; j++) {
    rc = fc_segment(subset, j, &numSegment, &segments);
    fail_unless(rc == FC_SUCCESS, "failed to segment empty subset");
    fail_unless(numSegment == 0, "empty subset should have no segments");
    fail_unless(segments == NULL, "empty subset should have null segements");
  }
  // test full subuset
  fc_addMemberToSubset(subset, 0);
  for (j = 0; j <= shared_dims[i]; j++) {
    rc = fc_segment(subset, j, &numSegment, &segments);
    fail_unless(rc == FC_SUCCESS, "failed to segment empty subset");
    fail_unless(numSegment == 1, "empty subset should have no segments");
    fc_getSubsetMembersAsArray(segments[0], &numTemp, &temps);
    fail_unless(numTemp == 1, "mismatch of numMember");
    fail_unless(temps[0] == 0, "mismatch of members");
    free(temps);
    free(segments);
  }
  fc_deleteSubset(subset);

  // cleanup
  fc_deleteMesh(mesh);
  fc_deleteDataset(dataset);
}
END_TEST
  
// Since underlying routines tested above, will just test on elements & whole
START_TEST(get_skin)
{
  FC_ReturnCode rc;
  int i, j;
  int M = 11, N = 12, P = 6;    // dimensions from data/mesh_generator.c
  FC_Dataset dataset;
  int numMesh = 4;
  char meshNames[4][1028] = { "point mesh", "line mesh", 
			      "quad mesh", "hex mesh" };
  FC_Mesh mesh, *returnMeshes;
  int numReturnMeshes;
  int numSubset = 2;
  FC_Subset subsets[2], skin;
  int numID, *IDs, numID_good, *IDs_good, numID_temp, *IDs_temp;
  FC_SortedIntArray sia;

  // setup
  fc_loadDataset("../data/gen_multimesh.ex2", &dataset);

  for (i = 0; i < numMesh; i++) {
    //printf("mesh: %s\n", meshNames[i]);
    rc = fc_getMeshByName(dataset, meshNames[i],  &numReturnMeshes,&returnMeshes);
    fail_unless(rc == FC_SUCCESS, "aborted - cant get mesh by name");
    fail_unless(numReturnMeshes == 1, "failed to find unique mesh by name");
    mesh = returnMeshes[0];
    free(returnMeshes);

    fc_getMeshNumEntity(mesh, FC_AT_ELEMENT, &numID);

    // create subsets for testing
    fc_createSubset(mesh, "whole", FC_AT_WHOLE_MESH, &subsets[0]);
    fc_addMemberToSubset(subsets[0], 0);
    fc_createSubset(mesh, "all elements", FC_AT_ELEMENT, &subsets[1]);
    IDs = malloc(numID*sizeof(int));
    for (j = 0; j < numID; j++)
      IDs[j] = j;
    fc_addArrayMembersToSubset(subsets[1], numID, IDs);
    free(IDs);

    // calculate good answers
    if (i == 1) {
      numID_good = 3;
      IDs_good = malloc(numID_good*sizeof(int));
      IDs_good[0] = 11;
      IDs_good[1] = 144;
      IDs_good[2] = 155;
    }
    else if (i == 2) {
      fc_initSortedIntArray(&sia);
      for (j = 0; j < M*N; j++) {
	fc_getMeshEntityChildren(mesh, FC_AT_ELEMENT, j, FC_AT_EDGE,
				 &numID_temp, &IDs_temp);
	if (j < M) // bottom edges
	  fc_addIntToSortedIntArray(&sia, IDs_temp[0]);
	if (j >= M*(N-1)) // top edges
	  fc_addIntToSortedIntArray(&sia, IDs_temp[2]);
	if (j%M == 0) // left edges
	  fc_addIntToSortedIntArray(&sia, IDs_temp[3]);
	if (j%M == M-1) // right edges
	  fc_addIntToSortedIntArray(&sia, IDs_temp[1]);
	free(IDs_temp);
      }
      fc_convertSortedIntArrayToIntArray(&sia, &numID_good, &IDs_good);
      fail_unless(numID_good == 2*(M+N), "test aborted: failure in calc");
    }
    else if (i == 3) {
      fc_initSortedIntArray(&sia);
      for (j = 0; j < M*N*P; j++) {
	fc_getMeshEntityChildren(mesh, FC_AT_ELEMENT, j, FC_AT_FACE,
				 &numID_temp, &IDs_temp);
	if (j < M*N) // bottom faces
	  fc_addIntToSortedIntArray(&sia, IDs_temp[4]);
	if (j >= M*N*(P-1)) // top faces
	  fc_addIntToSortedIntArray(&sia, IDs_temp[5]);
	if (j%M == 0) // left faces
	  fc_addIntToSortedIntArray(&sia, IDs_temp[3]);
	if (j%M == M-1) // right faces
	  fc_addIntToSortedIntArray(&sia, IDs_temp[1]);
	if (j%(M*N) < M) // front faces
	  fc_addIntToSortedIntArray(&sia, IDs_temp[0]);
	if (j%(M*N) >= M*(N-1)) // back faces
	  fc_addIntToSortedIntArray(&sia, IDs_temp[2]);
	free(IDs_temp);
      }
      fc_convertSortedIntArrayToIntArray(&sia, &numID_good, &IDs_good);
      fail_unless(numID_good == 2*(M*N + M*P + N*P),
		  "test aborted: failure in calc");
    }

    // get skin of mesh
    rc = fc_getMeshSkin(mesh, &skin);
    if (i == 0) {
      fail_unless(rc != FC_SUCCESS, "should fail on vertex mesh");
      fail_unless(FC_HANDLE_EQUIV(skin, FC_NULL_SUBSET), 
		  "fail should return null"); 
    }
    else {
      fail_unless(rc == FC_SUCCESS, "failed to get mesh skin");
      fc_getSubsetMembersAsArray(skin, &numID_temp, &IDs_temp);
      fail_unless(numID_temp == numID_good, "mismatch of numID");
      for (j = 0; j < numID_good; j++) 
	fail_unless(IDs_temp[j] == IDs_good[j], "mismatch of IDs");
      free(IDs_temp);
      fc_deleteSubset(skin);
    }

    // get skins of all subsets
    for (j = 0; j < numSubset; j++) {
      rc = fc_getSubsetSkin(subsets[j], &skin);
      if (i == 0) {
	fail_unless(rc != FC_SUCCESS, "should fail on vertex mesh");
	fail_unless(FC_HANDLE_EQUIV(skin, FC_NULL_SUBSET), 
		    "fail should return null"); 
      }
      else {
	fail_unless(rc == FC_SUCCESS, "failed to get subset skin");
	fc_getSubsetMembersAsArray(skin, &numID_temp, &IDs_temp);
	fail_unless(numID_temp == numID_good, "mismatch of numID");
	for (j = 0; j < numID_good; j++) 
	  fail_unless(IDs_temp[j] == IDs_good[j], "mismatch of IDs");
	free(IDs_temp);
	fc_deleteSubset(skin);
      }
    }

    // get skin should fail on empty subset
    fc_deleteSubset(subsets[1]);
    fc_createSubset(mesh, "empty subset", FC_AT_ELEMENT, &subsets[1]);
    rc = fc_getSubsetSkin(subsets[1], &skin);
    fail_unless(rc != FC_SUCCESS, "should fail on empty subset");
    fail_unless(FC_HANDLE_EQUIV(skin, FC_NULL_SUBSET), 
		"fail should return null");

    // skins of subsets which are lower dim than elements may be zero
    for (j = 1; j < i; j++) {
      FC_AssociationType assoc_temp;
      switch(j) {
      case 1:  assoc_temp = FC_AT_EDGE;    break;
      case 2:  assoc_temp = FC_AT_FACE;    break;
      }
      fc_getMeshEntityChildren(mesh, FC_AT_ELEMENT, j, assoc_temp,
				 &numID_temp, &IDs_temp);
      fc_deleteSubset(subsets[1]);
      fc_createSubset(mesh, "test subset", assoc_temp, &subsets[1]);
      fc_addArrayMembersToSubset(subsets[1], numID_temp, IDs_temp);
      free(IDs_temp);
      rc = fc_getSubsetSkin(subsets[1], &skin);
      fail_unless(rc == FC_SUCCESS, "should not fail");
      fc_getSubsetMembersAsArray(skin, &numID_temp, &IDs_temp);
      fail_unless(numID_temp == 0, "numID should be zero");
      fail_unless(IDs_temp == NULL, "IDs should be null");
      fc_deleteSubset(skin);
    }

    // test bad arguments for fc_getMeshSkin()
    rc = fc_getMeshSkin(FC_NULL_MESH, &skin);
    fail_unless(rc != FC_SUCCESS, "should fail to get NULL mesh");
    fail_unless(FC_HANDLE_EQUIV(skin, FC_NULL_SUBSET), 
		"fail should return null");
    rc = fc_getMeshSkin(mesh, NULL);
    fail_unless(rc != FC_SUCCESS, "should fail if skin == NULL");
   
    // test bad arguments for fc_getSubsetSkin()
    rc = fc_getSubsetSkin(FC_NULL_SUBSET, &skin);
    fail_unless(rc != FC_SUCCESS, "should fail to get NULL subset");
    fail_unless(FC_HANDLE_EQUIV(skin, FC_NULL_SUBSET), 
		"fail should return null");
    rc = fc_getSubsetSkin(subsets[0], NULL);
    fail_unless(rc != FC_SUCCESS, "should fail if skin == NULL");

    // cleanup
    if (i > 0)
      free(IDs_good);
    fc_deleteMesh(mesh);
  }
 
  // cleanup
  fc_deleteDataset(dataset);
}
END_TEST


// Populate the Suite with the tests

Suite *topo_suite(void)
{
  Suite *suite = suite_create("Topo");

  TCase *tc_membership = tcase_create(" - Membership Relations ");
  TCase *tc_neighbors = tcase_create(" - Get Neighbors ");
  TCase *tc_segment = tcase_create(" - Segmentation ");
  TCase *tc_skin = tcase_create(" - Skin ");

  suite_add_tcase(suite, tc_membership);
  tcase_add_checked_fixture(tc_membership, topo_setup, topo_teardown);
  tcase_add_test(tc_membership, get_children);
  tcase_add_test(tc_membership, get_parents);
  tcase_add_test(tc_membership, change_entity_type);

  suite_add_tcase(suite, tc_neighbors);
  tcase_add_checked_fixture(tc_neighbors, topo_setup, topo_teardown);
  tcase_add_test(tc_neighbors, entity_neighbors);
  tcase_add_test(tc_neighbors, subset_neighbors);

  suite_add_tcase(suite, tc_segment);
  tcase_add_checked_fixture(tc_segment, topo_setup, topo_teardown);
  tcase_add_test(tc_segment, segmentation);

  suite_add_tcase(suite, tc_skin);
  tcase_add_checked_fixture(tc_skin, topo_setup, topo_teardown);
  tcase_add_test(tc_skin, get_skin);

  return suite;
}
