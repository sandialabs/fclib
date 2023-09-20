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
 * \file checkmesh.c
 * \brief Unit testing for \ref Mesh module.
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/unittest/checkmesh.c,v $
 * $Revision: 1.49 $ 
 * $Date: 2006/10/19 03:14:53 $
 *
 * \description
 *
 *    To test the building of the mesh structs topo info (
 *    like edges and faces, and members and neighbors) we take
 *    advantage of the fact that some of the generated data
 *    are really structured meshes and we calculate the same
 *    info in a different way for comparison. So far we're
 *    only going to compare one representative mesh type from
 *    each topo dimension and use non-simplex (i.e. not tris
 *    and tets): POINT, LINE, QUAD & HEX.
 *
 * \modifications
 *    - 05/07/04 WSK created. Moved checking of mesh structures here from 
 *      checkbase.c
 *    - 09/13/04 WSK added checking of mesh interface (used to be in
 *      checklibrary.c).
 */

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <check.h>
#include "fc.h"
#include "fcP.h"
#include "checkall.h"

// *** Helpers for build conns/neighbors & fc_createSimpleHex mesh tests

// ICoords are i, j, k coordinates into a block structure
typedef int ICoords[3];

// copy ICoords
static void ICoordscpy(ICoords src, ICoords* dest) {
  int i;
  for (i = 0; i < 3; i++)
    (*dest)[i] = src[i];
}

// Convert from index to icoords (if number elements starting with icoords = 0,
// and varying x fastest)
static void getICoords(int idx, int numDim, int* dims, ICoords* icoords) {
  int i, j;
  int temp = idx;
  for (i = numDim -1; i > 0; i--) {
    int stride = 1;
    for (j = 0; j < i; j++)
      stride *= dims[j];
    (*icoords)[i] = temp/stride;
    temp = temp%stride;
  }
  (*icoords)[0] = temp;
}

// Convert from icoords to index, return index, or -1 if out of bounds
static int getIndex(int numDim, int* dims, ICoords icoords) {
  int i, j;
  int idx = 0;
  for (i = 0; i < numDim; i++) {
    if (icoords[i] < 0 || icoords[i] >= dims[i])
      return -1;
  }
  for (i = 0; i < numDim; i++) {
    int stride = 1;
    for (j = 0; j < i; j++)
      stride *= dims[j];
    idx += icoords[i]*stride;
  }
  return idx;
}

// **** Fixtures
static void mesh_setup(void) {
  FC_ReturnCode rc;
  if (isForking) {
    fc_setLibraryVerbosity(fc_messages);
    rc = fc_initLibrary();
    fail_unless(rc == FC_SUCCESS, "test aborted: failed to init library");
  }
}

static void mesh_teardown(void) {
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
  FC_Dataset dataset;
  FC_Mesh meshes[7];
  
  // create parent datset
  rc = fc_createDataset("fake", &dataset);
  fail_unless(rc == FC_SUCCESS, "failed to create dataset");

  // Create num
  for (i = 0; i < num; i++) {
    rc = fc_createMesh(dataset, "fake", &meshes[i]);
    fail_unless(rc == FC_SUCCESS, "failed to create mesh");
  }
  uID_start = meshes[0].uID;

  // Delete "odd" slots
  for (i = 1; i < num; i+=2) {
    rc = fc_deleteMesh(meshes[i]);
    fail_unless(rc == FC_SUCCESS, "failed to delete mesh");
  }

  // Create some more (fills in the cracks)
  for (i = 1; i < num; i+=2) {
    rc = fc_createMesh(dataset, "fake", &meshes[i]);
    fail_unless(rc == FC_SUCCESS, "failed to create mesh");
  }

  // Check slotID and uID
  for (i = 0; i < num; i++) {
    fail_unless(meshes[i].slotID == i, "mismatch of slot id");
    fail_unless(meshes[i].uID == uID_start + uID_incs[i],
		"mismatch of uID");
  }

  // cleanup
  rc = fc_deleteDataset(dataset);
  fail_unless(rc == FC_SUCCESS, "failed to delete dataset at end of test");
}
END_TEST

// test creating and deleting meshes, also test querying dataset
// for meshes before, inbetween and after.
START_TEST(create_get_delete)
{
  FC_ReturnCode rc;
  int i;
  int numMesh = 10;
  int temp_numMesh = 0;
  FC_Dataset dataset;
  FC_Mesh meshes[10], *temp_meshes, temp_mesh;
  FC_Mesh badMesh = { 999, 999 };
  char names[10][100] = { "one", "two", "three", "four", "five", "six",
			  "seven", "eight", "nine", "ten" };
  char newName[100] = { "sparkly new" };

  // make a dataset to play in
  rc = fc_createDataset("temp dataset", &dataset);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create dataset");

  // create some meshes
  for (i = 0; i < numMesh; i++) {
    rc = fc_createMesh(dataset, names[i], &meshes[i]);
    fail_unless(rc == FC_SUCCESS, "failed to create mesh");
    fail_unless(!FC_HANDLE_EQUIV(meshes[i], FC_NULL_MESH),
		"created mesh should not = FC_NULL_MESH");
  }

  // get meshes (should be numMesh);
  rc = fc_getMeshes(dataset, &temp_numMesh, &temp_meshes);
  fail_unless(rc == FC_SUCCESS, "failed to get meshes");
  fail_unless(temp_numMesh == numMesh, "mismatch of numMesh");
  for (i = 0; i < numMesh; i++) 
    fail_unless(FC_HANDLE_EQUIV(temp_meshes[i], meshes[i]), 
		"mismatch of mesh handles");
  free(temp_meshes);
  rc = fc_getNumMesh(dataset, &temp_numMesh);
  fail_unless(rc == FC_SUCCESS, "failed to get numMesh");
  fail_unless(temp_numMesh == numMesh, "mismatch of numMesh");
  for (i = 0; i < numMesh; i++) {
    rc = fc_getMeshByName(dataset, names[i], &temp_numMesh,&temp_meshes);
    fail_unless(rc == FC_SUCCESS, "failed to get mesh by name");
    fail_unless(temp_numMesh == 1, "wrong number of matching meshes by name");
    fail_unless(FC_HANDLE_EQUIV(temp_meshes[0], meshes[i]),
		"mismatch of mesh handles");
    free(temp_meshes);
  }

  // change the name of the first mesh
  rc = fc_changeMeshName(meshes[0], newName);
  fail_unless(rc == FC_SUCCESS, "failed to change mesh name");
  temp_numMesh = -5;
  rc = fc_getMeshByName(dataset, names[0],&temp_numMesh,
			&temp_meshes);
  fail_unless(rc == FC_SUCCESS, "should work for old name");
  fail_unless(temp_numMesh == 0, "should find no matches for old name");
  fail_unless(temp_meshes == NULL, "should return emtpy array for old name");

  rc = fc_getMeshByName(dataset, newName, &temp_numMesh,&temp_meshes);
  fail_unless(rc == FC_SUCCESS, "should work for new name");
  fail_unless(temp_numMesh == 1, "should find match for new name");
  fail_unless(FC_HANDLE_EQUIV(temp_meshes[0], meshes[0]),
	      "mismatch of mesh handle");
  free(temp_meshes);

  //make another mesh temporarily with the same name
  rc = fc_createMesh(dataset, newName, &temp_mesh);
  fail_unless(rc == FC_SUCCESS, "aborted: failed to create mesh mesh");
  rc = fc_getMeshByName(dataset, newName, &temp_numMesh,&temp_meshes);
  fail_unless(rc == FC_SUCCESS, "should work for new name");
  fail_unless(temp_numMesh == 2, "should find match for new name");
  fail_unless(((FC_HANDLE_EQUIV(temp_meshes[0], meshes[0]) &&
		(FC_HANDLE_EQUIV(temp_meshes[1], temp_mesh))) ||
	       (FC_HANDLE_EQUIV(temp_meshes[1], meshes[0]) &&
		(FC_HANDLE_EQUIV(temp_meshes[0], temp_mesh)))),
	      "mismatch of mesh handles");
  free(temp_meshes);
  fc_deleteMesh(temp_mesh);

  // delete half of the meshes (alternate)
  for (i = 0; i < numMesh; i+=2) {
    rc = fc_deleteMesh(meshes[i]);
    fail_unless(rc == FC_SUCCESS, "failed to delete mesh");
  }

  // get meshes (should be numMesh/2);
  fc_getMeshes(dataset, &temp_numMesh, &temp_meshes);
  fail_unless(rc == FC_SUCCESS, "failed to get meshes");
  fail_unless(temp_numMesh == numMesh/2, "mismatch of numMesh");
  for (i = 0; i < numMesh/2; i++) 
    fail_unless(FC_HANDLE_EQUIV(temp_meshes[i], meshes[i*2+1]), 
		"mismatch of mesh handles");
  free(temp_meshes);
  for (i = 0; i < numMesh/2; i++) {
    rc = fc_getMeshByName(dataset, names[i*2+1], &temp_numMesh,&temp_meshes);
    fail_unless(rc == FC_SUCCESS, "failed to get mesh by name");
    fail_unless(temp_numMesh == 1, "should find match for name");    
    fail_unless(FC_HANDLE_EQUIV(temp_meshes[0], meshes[i*2+1]),
		"mismatch of mesh handles");
    free(temp_meshes);
  }

  // delete remaining meshes
  for (i = 1; i < numMesh; i+=2) { 
    rc = fc_deleteMesh(meshes[i]);
    fail_unless(rc == FC_SUCCESS, "failed to delete mesh");
  }
  
  // get meshes (should be none)
  rc = fc_getMeshes(dataset, &temp_numMesh, &temp_meshes);
  fail_unless(rc == FC_SUCCESS, "failed to get meshes from empty library");
  fail_unless(temp_numMesh == 0 && temp_meshes == NULL,
	      "should return 0 if all meshes deleted");

  // make one more mesh for further testing
  rc = fc_createMesh(dataset, names[0], &meshes[0]);
  fail_unless(rc == FC_SUCCESS, "aborted: failed to create mesh mesh");

  // ---- test special cases

  // not an error to delete FC_NULL_MESH
  rc = fc_deleteMesh(FC_NULL_MESH);
  fail_unless(rc == FC_SUCCESS, "should not error to delete NULL mesh");

  // ---- test error conditions

  // bad args to fc_createMesh()
  temp_mesh = badMesh;
  rc = fc_createMesh(FC_NULL_DATASET, names[0], &temp_mesh);
  fail_unless(rc != FC_SUCCESS, 
	      "should have failed to create mesh with null database");
  fail_unless(FC_HANDLE_EQUIV(temp_mesh, FC_NULL_MESH),
	      "fail should return NULL mesh");
  temp_mesh = badMesh;
  rc = fc_createMesh(dataset, NULL, &temp_mesh);
  fail_unless(rc != FC_SUCCESS, 
	      "should have failed to create mesh with null name");
  fail_unless(FC_HANDLE_EQUIV(temp_mesh, FC_NULL_MESH),
	      "fail should return NULL mesh");
  rc = fc_createMesh(dataset, names[0], NULL);
  fail_unless(rc != FC_SUCCESS, 
	      "should have failed to create mesh with null handle");

  // bad args to fc_changeMeshName()
  rc = fc_changeMeshName(meshes[0], NULL);
  fail_unless(rc != FC_SUCCESS, "should fail if new name is NULL");
  rc = fc_getMeshByName(dataset, names[0], &temp_numMesh,&temp_meshes);
  fail_unless(rc == FC_SUCCESS && temp_numMesh == 1 &&
	      FC_HANDLE_EQUIV(meshes[0], temp_meshes[0]),
	      "should not change name of mesh");
  free(temp_meshes);

  // bad args to fc_getMeshes()
  temp_numMesh = 99;
  temp_meshes = (FC_Mesh*)1;
  rc = fc_getMeshes(FC_NULL_DATASET, &temp_numMesh, &temp_meshes);
  fail_unless(rc != FC_SUCCESS, 
	      "should have failed to create mesh with null dataset");
  fail_unless(temp_numMesh == -1, "fail should return -1 for numMesh"); 
  fail_unless(temp_meshes == NULL, "fail should return null");
  temp_meshes = (FC_Mesh*)1;
  rc = fc_getMeshes(dataset, NULL, &temp_meshes);
  fail_unless(rc != FC_SUCCESS, 
	      "should have failed to get mesh with NULL numMesh");
  fail_unless(temp_meshes == NULL, "fail should return null");
  temp_numMesh = 999;
  rc = fc_getMeshes(dataset, &temp_numMesh, NULL);
  fail_unless(rc != FC_SUCCESS, "should error to get just numMesh");
  fail_unless(temp_numMesh == -1, "mismatch of numMesh");
  
  // bad args to fc_getNumMesh()
  temp_numMesh = 99;
  rc = fc_getNumMesh(FC_NULL_DATASET, &temp_numMesh);
  fail_unless(rc != FC_SUCCESS, 
	      "should have failed to get numMesh with null dataset");
  fail_unless(temp_numMesh == -1, "fail should return -1 for numMesh"); 
  rc = fc_getNumMesh(dataset, NULL);
  fail_unless(rc != FC_SUCCESS, 
	      "should have failed to get numMesh with NULL numMesh");
  
  // bad args to fc_getMeshByName()
  rc = fc_getMeshByName(FC_NULL_DATASET, names[0], &temp_numMesh,&temp_meshes);
  fail_unless(rc != FC_SUCCESS, 
	      "should have failed to get mesh with null dataset");
  fail_unless(temp_numMesh == -1, "should return -1 when fail");
  fail_unless(temp_meshes == NULL, "fail should return NULL mesh array");

  rc = fc_getMeshByName(dataset, NULL, &temp_numMesh,&temp_meshes);
  fail_unless(rc != FC_SUCCESS, 
	      "should have failed to get mesh with null name");
  fail_unless(temp_numMesh == -1, "should return -1 when fail");
  fail_unless(temp_meshes == NULL, "fail should return NULL mesh array");

  rc = fc_getMeshByName(dataset, names[0], NULL,&temp_meshes);
  fail_unless(rc != FC_SUCCESS, 
	      "should have failed to get mesh with null arg");
  fail_unless(temp_meshes == NULL, "fail should return NULL mesh array");

  rc = fc_getMeshByName(dataset, names[0],&temp_numMesh,NULL);
  fail_unless(rc != FC_SUCCESS, 
	      "should have failed to get mesh with null arg");
  fail_unless(temp_numMesh == -1, "fail should return -1");

  // bad args to fc_deleteMesh()
  rc = fc_deleteMesh(badMesh);
  fail_unless(rc != FC_SUCCESS, 
	      "should error to delete nonexistent mesh");

  // --- done

  // delete last mesh
  rc = fc_deleteMesh(meshes[0]);
  fail_unless(rc == FC_SUCCESS, "failed to delete last mesh");

  // delete the dataset
  rc = fc_deleteDataset(dataset);
  fail_unless(rc == FC_SUCCESS, "failed to delete dataset at end of test");
}
END_TEST

// query meta data
START_TEST(metadata_query)
{
  FC_ReturnCode rc;
  char name[20] = "blue berry", *temp_name;
  FC_Dataset dataset, temp_dataset, badDataset = { 999, 999 };
  FC_Mesh mesh, badMesh = { 999, 999 };

  // make a dataset to play in
  rc = fc_createDataset("temp dataset", &dataset);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create dataset");

  // create a mesh to play with
  rc = fc_createMesh(dataset, name, &mesh);
  fail_unless(rc == FC_SUCCESS, "failed to created mesh");
 
  // test fc_isMeshValid()
  fail_unless(fc_isMeshValid(mesh), 
	      "failed to validate a valid mesh");

  // test fc_getMeshName()
  temp_name = NULL;
  rc = fc_getMeshName(mesh, &temp_name);
  fail_unless(rc == FC_SUCCESS, "failed to get mesh name");
  fail_unless(!strcmp(name, temp_name), "mismatch of name");
  // should have returned a copy, mess with temp_name & try again to make sure
  temp_name[0] = 'q';
  free(temp_name);
  rc = fc_getMeshName(mesh, &temp_name);
  fail_unless(rc == FC_SUCCESS, "failed to get mesh name");
  fail_unless(!strcmp(name, temp_name), "mismatch of name");
  free(temp_name);

  // test fc_getDatasetFromMesh()
  temp_dataset = badDataset;
  rc = fc_getDatasetFromMesh(mesh, &temp_dataset);
  fail_unless(rc == FC_SUCCESS, "failed to get parent dataset");
  fail_unless(FC_HANDLE_EQUIV(temp_dataset, dataset),
	      "mismatch of parent dataset");

  // --- check with bad args

  // fc_isMeshValid()
  fail_unless(!fc_isMeshValid(badMesh), 
	      "badMesh should NOT be valid");
  fail_unless(!fc_isMeshValid(FC_NULL_MESH), 
	      "FC_NULL_MESH should not be valid");

  // fc_getMeshName()
  temp_name = (char*)1;
  rc = fc_getMeshName(badMesh, &temp_name);
  fail_unless(rc != FC_SUCCESS, "badMesh should NOT return a name");
  fail_unless(temp_name == NULL, "fail should return NULL name");
  rc = fc_getMeshName(mesh, NULL);
  fail_unless(rc != FC_SUCCESS, "should fail for NULL name");

  // fc_getDatasetFromMesh()
  temp_dataset = badDataset;
  rc = fc_getDatasetFromMesh(badMesh, &temp_dataset);
  fail_unless(rc != FC_SUCCESS, "badMesh should NOT return a dataset");
  fail_unless(FC_HANDLE_EQUIV(temp_dataset, FC_NULL_DATASET),
	      "failure should return NULL dataset");
  rc = fc_getDatasetFromMesh(mesh, NULL);
  fail_unless(rc != FC_SUCCESS, "should fail for NULL dataset");

  // --- all done

  // cleanup & one last test
  fc_deleteMesh(mesh);
  fail_unless(!fc_isMeshValid(mesh), 
	      "handle should not be valid after a delete"); 

  // delete the dataset
  rc = fc_deleteDataset(dataset);
  fail_unless(rc == FC_SUCCESS, "failed to delete dataset at end of test");
}
END_TEST

// query the mesh coords and element conns
START_TEST(coords_conns_query)
{
  FC_ReturnCode rc;
  int i, j, k, m;
  char name[20] = "blue berry", name2[20] = "banana nut";
  FC_Dataset dataset;
  FC_Mesh mesh, mesh2, empty_mesh, badMesh = { 999, 999 };
  _FC_MeshSlot* meshSlot;
  int numElemType = 8;
  FC_ElementType elemTypes[8] = { FC_ET_POINT,  FC_ET_LINE,  FC_ET_TRI,
		      	          FC_ET_QUAD,   FC_ET_TET,   FC_ET_PYRAMID,
			          FC_ET_PRISM,  FC_ET_HEX };
  int topoDimPerType[8] = { 0, 1, 2, 2, 3, 3, 3, 3 };
  int numVertPerElem[8] = { 1, 2, 3, 4, 4, 5, 6, 8 };
  int numElem = 2; // two elements sharing edge/face
  int numVertPerType[8] = { 2, 3, 4, 6, 5, 6, 8, 12 };
  int conns[8][16] = { { 0,    1 },
		       { 0, 1,    1, 2 },
		       { 0, 1, 2,    1, 3, 2 },
		       { 0, 1, 2, 3,    1, 4, 5, 2 },
		       { 0, 1, 2, 3,    0, 2, 1, 4 },
		       { 0, 1, 2, 3, 4,    0, 3, 2, 1, 5 },
		       { 0, 1, 2, 3, 4, 5,    1, 6, 2, 4, 7, 5 },
		       { 0, 1, 2, 3, 4, 5, 6, 7,    1, 8, 9, 2, 5, 10, 11, 6 } };
  int numAssoc = 5;
  FC_AssociationType assocs[5] = { FC_AT_VERTEX,  FC_AT_EDGE, FC_AT_FACE,
				   FC_AT_ELEMENT, FC_AT_WHOLE_MESH };
  int numEntities[8][5] = { {  2,  0,  0, 2, 1 },
			    {  3,  2,  0, 2, 1 },
			    {  4,  5,  2, 2, 1 },
			    {  6,  7,  2, 2, 1 },
			    {  5,  9,  7, 2, 1 },
			    {  6, 12,  9, 2, 1 },
			    {  8, 14,  9, 2, 1 },
			    { 12, 20, 11, 2, 1 } };
  FC_ElementType entityFaceTypes[8] = { FC_ET_UNKNOWN, FC_ET_UNKNOWN,
					FC_ET_TRI, FC_ET_QUAD,
					FC_ET_TRI, FC_ET_MIXED, 
					FC_ET_MIXED, FC_ET_QUAD };
  FC_ElementType pyramidFaceTypes[9] = { FC_ET_TRI, FC_ET_TRI, FC_ET_TRI,
					 FC_ET_TRI, FC_ET_QUAD, FC_ET_TRI,
					 FC_ET_TRI, FC_ET_TRI, FC_ET_TRI };
  FC_ElementType prismFaceTypes[9] = { FC_ET_QUAD, FC_ET_QUAD, FC_ET_QUAD,
				       FC_ET_TRI, FC_ET_TRI, FC_ET_QUAD,
				       FC_ET_QUAD, FC_ET_TRI, FC_ET_TRI };  
  double coords[3*12];
  int temp_topoDim, temp_dim, temp_numVert, temp_numElem, temp_numEntity;
  int temp_numEdge, temp_numFace;
  FC_ElementType temp_elemType, temp_faceType;
  double *temp_coords, *coords_cpy;
  int *temp_conns, *temp_conns_helper, *conns_cpy;
  
  // make a dataset to play in
  rc = fc_createDataset("temp dataset", &dataset);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create dataset");

  // setup coords (these make no physical sense)
  for (i = 0; i < 3*numVertPerType[numElemType-1]; i++)
    coords[i] = i + i/10.;

  for (i = 0; i < numElemType; i++) { // loop over elem types
    for (j = 1; j <= 3; j++) { // loop over space dimensions
      // skip if it doesn't make sense to make this elem type at this dim
      if (topoDimPerType[i] > j)
	continue;

      //printf("type = %s, dim = %d\n", fc_getElementTypeText(elemTypes[i]),j); 

      // create three meshes
      rc = fc_createMesh(dataset, name, &mesh);
      fail_unless(rc == FC_SUCCESS, "abort: failed to created mesh");
      rc = fc_createMesh(dataset, name2, &mesh2);
      fail_unless(rc == FC_SUCCESS, "abort: failed to created mesh");
      rc = fc_createMesh(dataset, "empty mesh", &empty_mesh);
      fail_unless(rc == FC_SUCCESS, "abort: failed to created mesh");
      
      // check default values
      rc = fc_getMeshInfo(mesh, &temp_topoDim, &temp_dim, &temp_numVert,
			  &temp_numElem, &temp_elemType);
      fail_unless(rc == FC_SUCCESS, "failed to get mesh info");
      fail_unless(temp_topoDim == 0 && temp_dim == 0 && temp_numVert == 0 &&
		  temp_numElem == 0 && temp_elemType == FC_ET_UNKNOWN,
		  "empty mesh should have empty values");

      // mesh1: set/get coords & conns by copy and check
      rc = fc_setMeshCoords(mesh, j, numVertPerType[i], coords);
      fail_unless(rc == FC_SUCCESS, "failed to set vertex coords");
      rc = fc_setMeshElementConns(mesh, elemTypes[i], numElem, conns[i]);
      fail_unless(rc == FC_SUCCESS, "failed to set element conns");
      rc = fc_getMeshInfo(mesh, &temp_topoDim, &temp_dim, &temp_numVert,
			  &temp_numElem, &temp_elemType);
      fail_unless(rc == FC_SUCCESS, "failed to get info");
      fail_unless(temp_topoDim == topoDimPerType[i], "mismatch of topoDim");
      fail_unless(temp_dim == j, "mismatch of dim");
      fail_unless(temp_numVert == numVertPerType[i], "mistmatch of numVert");
      fail_unless(temp_numElem == numElem, "mismatch of numElem");
      fail_unless(temp_elemType == elemTypes[i], "mismatch of elemType");
      temp_coords = NULL;
      rc = fc_getMeshCoordsPtr(mesh, &temp_coords);
      fail_unless(rc == FC_SUCCESS, "failed to get coords");
      fail_unless(!memcmp(temp_coords, coords,
		  j*numVertPerType[i]*sizeof(double)), "mismatch of coords");
      temp_conns = NULL;
      rc = fc_getMeshElementConnsPtr(mesh, &temp_conns);
      fail_unless(rc == FC_SUCCESS, "failed to get conns");
      fail_unless(!memcmp(temp_conns, conns[i],
		  numElem*numVertPerElem[i]*sizeof(int)), "mismatch of conns");

      // mesh2: set/get coords & conns pointer and check
      coords_cpy = malloc(j*numVertPerType[i]*sizeof(double));
      memcpy(coords_cpy, coords, j*numVertPerType[i]*sizeof(double));
      conns_cpy = malloc(numElem*numVertPerElem[i]*sizeof(int));
      memcpy(conns_cpy, conns[i], numElem*numVertPerElem[i]*sizeof(int));
      rc = fc_setMeshCoordsPtr(mesh2, j, numVertPerType[i], coords_cpy);
      fail_unless(rc == FC_SUCCESS, "failed to set vertex coords");
      rc = fc_setMeshElementConnsPtr(mesh2, elemTypes[i], numElem, conns_cpy);
      fail_unless(rc == FC_SUCCESS, "failed to set element conns");
      rc = fc_getMeshInfo(mesh2, &temp_topoDim, &temp_dim, &temp_numVert,
			  &temp_numElem, &temp_elemType);
      fail_unless(rc == FC_SUCCESS, "failed to get info");
      fail_unless(temp_topoDim == topoDimPerType[i], "mismatch of topoDim");
      fail_unless(temp_dim == j, "mismatch of dim");
      fail_unless(temp_numVert == numVertPerType[i], "mistmatch of numVert");
      fail_unless(temp_numElem == numElem, "mismatch of numElem");
      fail_unless(temp_elemType == elemTypes[i], "mismatch of elemType");
      temp_coords = NULL;
      rc = fc_getMeshCoordsPtr(mesh2, &temp_coords);
      fail_unless(rc == FC_SUCCESS, "failed to get coords");
      fail_unless(!memcmp(temp_coords, coords,
		  j*numVertPerType[i]*sizeof(double)), "mismatch of coords");
      temp_conns = NULL;
      rc = fc_getMeshElementConnsPtr(mesh2, &temp_conns);
      fail_unless(rc == FC_SUCCESS, "failed to get conns");
      fail_unless(!memcmp(temp_conns, conns[i],
		  numElem*numVertPerElem[i]*sizeof(int)), "mismatch of conns");

      // test set by copy vs. set the pointer
      meshSlot = _fc_getMeshSlot(mesh);
      fail_unless(meshSlot->coords != coords, "should have copied the coords");
      meshSlot = _fc_getMeshSlot(mesh2);
      fail_unless(meshSlot->coords == coords_cpy, "should have ptr to coords");

      // check other query functions that are simpler forms of fc_getMeshInfo
      rc = fc_getMeshTopodim(mesh, &temp_topoDim);
      fail_unless(rc == FC_SUCCESS, "failed to get topodim");
      fail_unless(temp_topoDim == topoDimPerType[i], "mismatch of topoDim");
      rc = fc_getMeshDim(mesh, &temp_dim);
      fail_unless(rc == FC_SUCCESS, "failed to get dim");
      fail_unless(temp_dim == j, "mismatch of dim");
      rc = fc_getMeshNumVertex(mesh, &temp_numVert);
      fail_unless(rc == FC_SUCCESS, "failed to get numVertex");
      fail_unless(temp_numVert == numVertPerType[i], "mistmatch of numVert");
      rc = fc_getMeshNumElement(mesh, &temp_numElem);
      fail_unless(rc == FC_SUCCESS, "failed to get numElement");
      fail_unless(temp_numElem == numElem, "mismatch of numElem");
      rc = fc_getMeshElementType(mesh, &temp_elemType);
      fail_unless(rc == FC_SUCCESS, "failed to get elementType");
      fail_unless(temp_elemType == elemTypes[i], "mismatch of elemType");

      // check getting edges and faces
      rc = fc_getMeshEdgeFaceInfo(mesh, &temp_numEdge, &temp_numFace, 
				  &temp_faceType);
      fail_unless(rc == FC_SUCCESS, "failed to get edge/face info");
      fail_unless(temp_numEdge == numEntities[i][1], "mismatch of numEdge");
      fail_unless(temp_numFace == numEntities[i][2], "mismatch of numFace");
      fail_unless(temp_faceType == entityFaceTypes[i], "mismatch of faceType");
      // (just test that conns returns something, will test actual values
      // later in this file on specific meshes (hard to calc))
      rc = fc_getMeshEdgeConnsPtr(mesh, &temp_conns);
      fail_unless(rc == FC_SUCCESS, "should not fail");
      if (numEntities[i][1] > 0) 
	fail_unless(temp_conns != NULL, "conns should not be null");
      else
	fail_unless(temp_conns == NULL, "conns should be null");
      rc = fc_getMeshFaceConnsPtr(mesh, &temp_conns_helper, &temp_numVert,
				  &temp_conns);
      fail_unless(rc == FC_SUCCESS, "should not fail");
      if (numEntities[i][2] > 0) {
	int maxNumVert = fc_getElementTypeNumVertex(temp_faceType);
	if (temp_faceType == FC_ET_MIXED)
	  maxNumVert = 4;
	fail_unless(temp_numVert == maxNumVert, "mismatch of maxNumVertPerFace");
	for (k = 0; k < temp_numFace; k++) {
	  for (m = temp_conns_helper[k]; m < maxNumVert; m++)
	    fail_unless(temp_conns[maxNumVert*k+m] == -1, "mismatch of conns");
	}
      }
      else 
	fail_unless(temp_conns_helper == NULL && temp_numVert == 0 &&
		    temp_conns == NULL, "conns should be NULL");

      // check query functions that are simpler forms of fc_getEdgeFaceMeshInfo
      rc = fc_getMeshNumEdge(mesh, &temp_numEdge);
      fail_unless(rc == FC_SUCCESS && temp_numEdge == numEntities[i][1],
		  "mismatch of numEdge");
      rc = fc_getMeshNumFace(mesh, &temp_numFace);
      fail_unless(rc == FC_SUCCESS && temp_numFace == numEntities[i][2],
		  "mismatch of numFace");
      rc = fc_getMeshFaceType(mesh, &temp_faceType);
      fail_unless(rc == FC_SUCCESS && temp_faceType == entityFaceTypes[i],
		  "mismatch of faceType");

      // check mesh entity queries
      for (k = 0; k < numAssoc; k++) {
	rc = fc_getMeshNumEntity(mesh, assocs[k], &temp_numEntity);
	fail_unless(rc == FC_SUCCESS, "failed to get numEntity");
	fail_unless(temp_numEntity == numEntities[i][k],
		    "mismatch of numEntity");
	for (m = 0; m < numEntities[i][k]; m++) {
	  FC_ElementType goodType;
	  switch (assocs[k]) {
	  case FC_AT_VERTEX:  goodType = FC_ET_POINT;  break;
	  case FC_AT_EDGE:    goodType = FC_ET_LINE;   break;
	  case FC_AT_FACE: {
	    if (entityFaceTypes[i] == FC_ET_MIXED) {
	      if (elemTypes[i] == FC_ET_PYRAMID)
		goodType = pyramidFaceTypes[m];
	      else if (elemTypes[i] == FC_ET_PRISM)
		goodType = prismFaceTypes[m];
	      else
	      fail_unless(1, "abort: error in testing code");
	    }
	    else
	      goodType = entityFaceTypes[i];
	  } break;
	  case FC_AT_ELEMENT: // fall through to FC_AT_WHOLE_MESH
	  case FC_AT_WHOLE_MESH: goodType = elemTypes[i]; break;
	  default: fail_unless(1, "abort: error in testing code");
	  }
	  rc = fc_getMeshEntityElementType(mesh, assocs[k], m, &temp_elemType);
	  fail_unless(rc == FC_SUCCESS, "failed to get entity type");
	  fail_unless(temp_elemType == goodType, "mismatch of entity type");
	}
      }

      // ---- check special cases
      
      // fc_getMeshInfo -- most arguments are optional
      rc = fc_getMeshInfo(mesh, &temp_topoDim, NULL, NULL, NULL, NULL);
      fail_unless(rc == FC_SUCCESS, "failed to get info");
      fail_unless(temp_topoDim == topoDimPerType[i], "mismatch of topoDim");
      rc = fc_getMeshInfo(mesh, NULL, &temp_dim, NULL, NULL, NULL);
      fail_unless(rc == FC_SUCCESS, "failed to get info");
      fail_unless(temp_dim == j, "mismatch of dim");
      rc = fc_getMeshInfo(mesh, NULL, NULL, &temp_numVert, NULL, NULL);
      fail_unless(rc == FC_SUCCESS, "failed to get info");
      fail_unless(temp_numVert == numVertPerType[i], "mistmatch of numVert");
      rc = fc_getMeshInfo(mesh, NULL, NULL, NULL, &temp_numElem, NULL);
      fail_unless(rc == FC_SUCCESS, "failed to get info");
      fail_unless(temp_numElem == numElem, "mismatch of numElem");
      rc = fc_getMeshInfo(mesh, NULL, NULL, NULL, NULL, &temp_elemType);
      fail_unless(rc == FC_SUCCESS, "failed to get info");
      fail_unless(temp_elemType == elemTypes[i], "mismatch of elemType");

      // fc_getMeshEdgeFaceInfo -- most arguments are optional
      rc = fc_getMeshEdgeFaceInfo(mesh, &temp_numEdge, NULL, NULL);
      fail_unless(rc == FC_SUCCESS, "failed to get edge/face info");
      fail_unless(temp_numEdge == numEntities[i][1], "mismatch of numEdge");
      rc = fc_getMeshEdgeFaceInfo(mesh, NULL, &temp_numFace, NULL);
      fail_unless(rc == FC_SUCCESS, "failed to get edge/face info");
      fail_unless(temp_numFace == numEntities[i][2], "mismatch of numFace");
      rc = fc_getMeshEdgeFaceInfo(mesh, NULL, NULL, &temp_faceType);
      fail_unless(rc == FC_SUCCESS, "failed to get edge/face info");
      fail_unless(temp_faceType == entityFaceTypes[i], "mismatch of faceType");

      // fc_getMeshFaceConnsPtr -- numVertPerFace is optional
      rc = fc_getMeshFaceConnsPtr(mesh, NULL, NULL, &temp_conns);
      fail_unless(rc == FC_SUCCESS, "should not fail");
      if (numEntities[i][2] > 0)
	fail_unless(temp_conns != NULL, "should return some conns");
      else 
	fail_unless(temp_conns == NULL, "conns should be NULL");

      // ---- check errors
      
      // fc_getMeshInfo() -- bad args
      rc = fc_getMeshInfo(badMesh, &temp_topoDim, &temp_dim, &temp_numVert,
			  &temp_numElem, &temp_elemType);
      fail_unless(rc != FC_SUCCESS, "should fail with bad mesh");
      fail_unless(temp_topoDim == -1 && temp_dim == -1 && temp_numVert == -1 &&
		  temp_numElem == -1 && temp_elemType == FC_ET_UNKNOWN,
		  "fail should return NULLs");
      rc = fc_getMeshInfo(mesh, NULL, NULL, NULL, NULL, NULL);
      fail_unless(rc != FC_SUCCESS, "can't have all args be NULL");

      // functions similar to fc_getMeshInfo() -- bad args
      rc = fc_getMeshTopodim(badMesh, &temp_topoDim);
      fail_unless(rc != FC_SUCCESS && temp_topoDim == -1, 
		  "should fail to get topodim from bad mesh");
      rc = fc_getMeshTopodim(mesh, NULL);
      fail_unless(rc != FC_SUCCESS, "should fail if topoDim = NULL");
      rc = fc_getMeshDim(badMesh, &temp_dim);
      fail_unless(rc != FC_SUCCESS && temp_dim == -1, 
		  "should failed to get dim from bad mesh");
      rc = fc_getMeshDim(mesh, NULL);
      fail_unless(rc != FC_SUCCESS, "should fail if dim = NULL");
      rc = fc_getMeshNumVertex(badMesh, &temp_numVert);
      fail_unless(rc != FC_SUCCESS && temp_numVert == -1, 
		  "should fail to get numVertex from bad mesh");
      rc = fc_getMeshNumVertex(mesh, NULL);
      fail_unless(rc != FC_SUCCESS, "should fail if numVert = NULL");
      rc = fc_getMeshNumElement(badMesh, &temp_numElem);
      fail_unless(rc != FC_SUCCESS && temp_numElem == -1, 
		  "should fail to get numElement from bad mesh");
      rc = fc_getMeshNumElement(mesh, NULL);
      fail_unless(rc != FC_SUCCESS, "should fail if numElement = NULL");
      rc = fc_getMeshElementType(badMesh, &temp_elemType);
      fail_unless(rc != FC_SUCCESS && temp_elemType == FC_ET_UNKNOWN, 
		  "should fail to get elementType from bad mesh"); 
      rc = fc_getMeshElementType(mesh, NULL);
      fail_unless(rc != FC_SUCCESS, "should failed if elementType = unknown");

      // fc_getMeshEdgeFaceInfo() -- bad args
      rc = fc_getMeshEdgeFaceInfo(badMesh, &temp_numEdge, &temp_numFace, 
				  &temp_faceType);
      fail_unless(rc != FC_SUCCESS, "should fail with bad mesh");
      fail_unless(temp_numEdge == -1 && temp_numFace == -1 && 
		  temp_faceType == FC_ET_UNKNOWN, "should return NULLs");
      rc = fc_getMeshEdgeFaceInfo(mesh, NULL, NULL, NULL);
      fail_unless(rc != FC_SUCCESS, "can't have all args be NULL");

      // functions similar to fc_getMeshEdgeFaceInfo() -- bad args
      rc = fc_getMeshNumEdge(badMesh, &temp_numEdge);
      fail_unless(rc != FC_SUCCESS && temp_numEdge == -1,
		  "should fail with bad mesh");
      rc = fc_getMeshNumEdge(mesh, NULL);
      fail_unless(rc != FC_SUCCESS, "should fail with NULL numEdge");
      rc = fc_getMeshNumFace(badMesh, &temp_numFace);
      fail_unless(rc != FC_SUCCESS && temp_numFace == -1,
		  "should fail with bad mesh");
      rc = fc_getMeshNumFace(mesh, NULL);
      fail_unless(rc != FC_SUCCESS, "should fail with NULL numFace");
      rc = fc_getMeshFaceType(badMesh, &temp_faceType);
      fail_unless(rc != FC_SUCCESS && temp_faceType == FC_ET_UNKNOWN,
		  "should fail with bad mesh");
      rc = fc_getMeshFaceType(mesh, NULL);
      fail_unless(rc != FC_SUCCESS, "should fail with NULL faceType");

      // fc_getMeshNumEntity() -- bad args
      rc = fc_getMeshNumEntity(badMesh, assocs[0], &temp_numEntity);
      fail_unless(rc != FC_SUCCESS && temp_numEntity == -1,
		  "should fail to get numEntity from bad mesh");
      rc = fc_getMeshNumEntity(mesh, FC_AT_UNKNOWN, &temp_numEntity);
      fail_unless(rc != FC_SUCCESS && temp_numEntity == -1,
		  "should fail to get numEntity of unknown association");
      rc = fc_getMeshNumEntity(mesh, -99, &temp_numEntity);
      fail_unless(rc != FC_SUCCESS && temp_numEntity == -1,
		  "should fail to get numEntity of impossible association");
      rc = fc_getMeshNumEntity(mesh, assocs[0], NULL);
      fail_unless(rc != FC_SUCCESS, "should fail if numEntity = NULL");

      // fc_getMeshEntityElemenType() -- bad args
      rc = fc_getMeshEntityElementType(badMesh, assocs[0], 0, &temp_elemType);
      fail_unless(rc != FC_SUCCESS && temp_elemType == FC_ET_UNKNOWN,
		  "should fail to get entityType from bad mesh");
      rc = fc_getMeshEntityElementType(mesh, FC_AT_UNKNOWN, 0, &temp_elemType);
      fail_unless(rc != FC_SUCCESS && temp_elemType == FC_ET_UNKNOWN,
		  "should fail to get entityType of unknown association");
      rc = fc_getMeshEntityElementType(mesh, -1, 0, &temp_elemType);
      fail_unless(rc != FC_SUCCESS && temp_elemType == FC_ET_UNKNOWN,
		  "should fail to get numEntity of impossible association");
      rc = fc_getMeshEntityElementType(mesh, assocs[0], -1, &temp_elemType);
      fail_unless(rc != FC_SUCCESS && temp_elemType == FC_ET_UNKNOWN,
		  "should fail to get entity type if ID out of range");
      rc = fc_getMeshEntityElementType(mesh, assocs[0], 999, &temp_elemType);
      fail_unless(rc != FC_SUCCESS && temp_elemType == FC_ET_UNKNOWN,
		  "should fail to get entity type if ID out of range");
      rc = fc_getMeshEntityElementType(mesh, assocs[0], 0, NULL);
      fail_unless(rc != FC_SUCCESS, "should fail if entity type = NULL");

      // fc_getMeshCoordsPtr() & fc_getMeshElementConnsPtr() -- bad args
      rc = fc_getMeshCoordsPtr(badMesh, &temp_coords);
      fail_unless(rc != FC_SUCCESS && temp_coords == NULL, 
		  "should fail with bad mesh");
      rc = fc_getMeshCoordsPtr(mesh, NULL);
      fail_unless(rc != FC_SUCCESS, "should fail if coords = NULL");
      rc = fc_getMeshElementConnsPtr(badMesh, &temp_conns);
      fail_unless(rc != FC_SUCCESS && temp_conns == NULL, 
		  "should fail with bad mesh");
      rc = fc_getMeshElementConnsPtr(mesh, NULL);
      fail_unless(rc != FC_SUCCESS, "should fail if conns = NULL");

      // fc_getMeshEdgeConnsPtr() & fc_getMeshFaceConnsPtr() -- bad args
      rc = fc_getMeshEdgeConnsPtr(badMesh, &temp_conns);
      fail_unless(rc != FC_SUCCESS && temp_conns == NULL,
		  "should fail with bad mesh");
      rc = fc_getMeshEdgeConnsPtr(mesh, NULL);
      fail_unless(rc != FC_SUCCESS, "should fail if conns = NULL");
      rc = fc_getMeshFaceConnsPtr(badMesh, &temp_conns_helper, &temp_numVert,
				  &temp_conns);
      fail_unless(rc != FC_SUCCESS && temp_conns_helper == NULL &&
		  temp_numVert == -1 && temp_conns == NULL, 
		  "should fail with bad mesh");
      rc = fc_getMeshFaceConnsPtr(mesh, &temp_conns_helper, &temp_numVert,
				  NULL);
      fail_unless(rc != FC_SUCCESS && temp_conns_helper == NULL &&
		  temp_numVert == -1,
		  "should fail with NULL conns");
      
      // fc_setMeshCoords() & fc_setMeshElementConns()-- error to reset
      rc = fc_setMeshCoords(mesh, j, numVertPerType[i], coords);
      fail_unless(rc != FC_SUCCESS, "shouldn't be able to reset coords");
      rc = fc_setMeshElementConns(mesh, elemTypes[i], numElem, conns[i]);
      fail_unless(rc != FC_SUCCESS, "shouldn't be able to reset conns");

      // fc_setMeshCoords() & fc_setMeshElementConns()--bad args 
      rc = fc_setMeshCoords(badMesh, j, numVertPerType[i], coords);
      fail_unless(rc != FC_SUCCESS, "should fail with bad mesh");
      rc = fc_setMeshElementConns(badMesh, elemTypes[i], numElem, conns[i]);
      fail_unless(rc != FC_SUCCESS, "should fail with bad mesh");
      rc = fc_setMeshCoords(empty_mesh, 0, numVertPerType[i], coords);
      fail_unless(rc != FC_SUCCESS, "should fail if dim out of range");
      rc = fc_setMeshCoords(empty_mesh, 4, numVertPerType[i], coords);
      fail_unless(rc != FC_SUCCESS, "should fail if dim out of range");
      rc = fc_setMeshCoords(empty_mesh, j, 0, coords);
      fail_unless(rc != FC_SUCCESS, "should fail if numVertex < 1");
      rc = fc_setMeshCoords(empty_mesh, j, numVertPerType[i], NULL);
      fail_unless(rc != FC_SUCCESS, "should fail if no coords");
      rc = fc_setMeshElementConns(empty_mesh, FC_ET_UNKNOWN, numElem, conns[i]);
      fail_unless(rc != FC_SUCCESS, "should fail if element type UNKNOWN");
      rc = fc_setMeshElementConns(empty_mesh, FC_ET_MIXED, numElem, conns[i]);
      fail_unless(rc != FC_SUCCESS, "should fail if element type is MIXED");
      rc = fc_setMeshElementConns(empty_mesh, -999, numElem, conns[i]);
      fail_unless(rc != FC_SUCCESS, "should fail if impossible elem type"); 
      rc = fc_setMeshElementConns(empty_mesh, elemTypes[i], 0, conns[i]);
      fail_unless(rc != FC_SUCCESS, "should fail if numElem < 1");
      rc = fc_setMeshElementConns(empty_mesh, elemTypes[i], numElem, NULL);
      fail_unless(rc != FC_SUCCESS, "should fail if no conns");

      // fc_setMeshCoordsPtr() & fc_setMeshElementConnsPtr)-- error to reset
      rc = fc_setMeshCoordsPtr(mesh, j, numVertPerType[i], coords);
      fail_unless(rc != FC_SUCCESS, "shouldn't be able to reset coords");
      rc = fc_setMeshElementConnsPtr(mesh, elemTypes[i], numElem, conns[i]);
      fail_unless(rc != FC_SUCCESS, "shouldn't be able to reset conns");

      // fc_setMeshCoordsPtr() & fc_setMeshElementConnsPtr)--bad args 
      rc = fc_setMeshCoordsPtr(badMesh, j, numVertPerType[i], coords);
      fail_unless(rc != FC_SUCCESS, "should fail with bad mesh");
      rc = fc_setMeshElementConnsPtr(badMesh, elemTypes[i], numElem, conns[i]);
      fail_unless(rc != FC_SUCCESS, "should fail with bad mesh");
      rc = fc_setMeshCoordsPtr(empty_mesh, 0, numVertPerType[i], coords);
      fail_unless(rc != FC_SUCCESS, "should fail if dim out of range");
      rc = fc_setMeshCoordsPtr(empty_mesh, 4, numVertPerType[i], coords);
      fail_unless(rc != FC_SUCCESS, "should fail if dim out of range");
      rc = fc_setMeshCoordsPtr(empty_mesh, j, 0, coords);
      fail_unless(rc != FC_SUCCESS, "should fail if numVertex < 1");
      rc = fc_setMeshCoordsPtr(empty_mesh, j, numVertPerType[i], NULL);
      fail_unless(rc != FC_SUCCESS, "should fail if no coords");
      rc = fc_setMeshElementConnsPtr(empty_mesh, FC_ET_UNKNOWN, numElem, conns[i]);
      fail_unless(rc != FC_SUCCESS, "should fail if element type UNKNOWN");
      rc = fc_setMeshElementConnsPtr(empty_mesh, FC_ET_MIXED, numElem, conns[i]);
      fail_unless(rc != FC_SUCCESS, "should fail if element type is MIXED");
      rc = fc_setMeshElementConnsPtr(empty_mesh, -999, numElem, conns[i]);
      fail_unless(rc != FC_SUCCESS, "should fail if impossible elem type"); 
      rc = fc_setMeshElementConnsPtr(empty_mesh, elemTypes[i], 0, conns[i]);
      fail_unless(rc != FC_SUCCESS, "should fail if numElem < 1");
      rc = fc_setMeshElementConnsPtr(empty_mesh, elemTypes[i], numElem, NULL);
      fail_unless(rc != FC_SUCCESS, "should fail if no conns");

      // empty mesh should still be empty (all set's failed)
      rc = fc_getMeshInfo(empty_mesh, &temp_topoDim, &temp_dim, &temp_numVert,
			  &temp_numElem, &temp_elemType);
      fail_unless(rc == FC_SUCCESS, "failed to get mesh info");
      fail_unless(temp_topoDim == 0 && temp_dim == 0 && temp_numVert == 0 &&
		  temp_numElem == 0 && temp_elemType == FC_ET_UNKNOWN,
		  "empty mesh should have empty values");
      
      // --- all done
      // cleanup
      rc = fc_deleteMesh(mesh);
      fail_unless(rc == FC_SUCCESS, 
		  "failed to delete mesh at end of tests");
      rc = fc_deleteMesh(mesh2);
      fail_unless(rc == FC_SUCCESS, 
		  "failed to delete mesh at end of tests");
      rc = fc_deleteMesh(empty_mesh);
      fail_unless(rc == FC_SUCCESS, 
		  "failed to delete mesh at end of tests");
    }
  }

  // delete the dataset
  rc = fc_deleteDataset(dataset);
  fail_unless(rc == FC_SUCCESS, "failed to delete dataset at end of test");
}
END_TEST

// copy
START_TEST(copy_test)
{
  FC_ReturnCode rc;
  int i, j;
  char name[20] = "blue berry", copy_name[20] = "banana nut", *temp_name;
  FC_Dataset dataset1, dataset2;
  int numMesh1, numMesh2, temp_numMesh;
  FC_Mesh mesh, copy_mesh, badMesh = { 999, 999 };
  int numElemType = 8;
  FC_ElementType elemTypes[8] = { FC_ET_POINT,  FC_ET_LINE,  FC_ET_TRI,
		      	          FC_ET_QUAD,   FC_ET_TET,   FC_ET_PYRAMID,
			          FC_ET_PRISM,  FC_ET_HEX };
  int topoDimPerType[8] = { 0, 1, 2, 2, 3, 3, 3, 3 };
  int numVertPerElem[8] = { 1, 2, 3, 4, 4, 5, 6, 8 };
  int numElem = 2; // two elements sharing edge/face
  int numVertPerType[8] = { 2, 3, 4, 6, 5, 6, 8, 12 };
  int conns[8][16] = { { 0,    1 },
		       { 0, 1,    1, 2 },
		       { 0, 1, 2,    1, 3, 2 },
		       { 0, 1, 2, 3,    1, 4, 5, 2 },
		       { 0, 1, 2, 3,    0, 2, 1, 4 },
		       { 0, 1, 2, 3, 4,    0, 3, 2, 1, 5 },
		       { 0, 1, 2, 3, 4, 5,    1, 6, 2, 4, 7, 5 },
		       { 0, 1, 2, 3, 4, 5, 6, 7,    1, 8, 9, 2, 5, 10, 11, 6 } };
  double coords[3*12];
  int temp_topoDim, temp_dim, temp_numVert, temp_numElem;
  FC_ElementType temp_elemType;
  double *temp_coords;
  int *temp_conns;

  FC_Coords lowers = { 0., 0., 0. };
  FC_Coords uppers = { 1., 1., 1. };
  FC_DataType returndatatype;
  FC_Sequence seq, retseq1, retseq2;
  FC_Variable var, *seqvar1, *seqvar2, globalvar, *globalseqvar, **returnseqvars; 
  FC_Subset sub, *seqsub, **returnseqsubs;
  int *data;
  int numStep, numMembers, numVar;
  int numreturn, *numstepperseq;
  
  // make two datasets to play in
  rc = fc_createDataset("temp dataset1", &dataset1);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create dataset");
  rc = fc_createDataset("temp dataset2", &dataset2);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create dataset");

  // setup coords (these make no physical sense)
  for (i = 0; i < 3*numVertPerType[numElemType-1]; i++)
    coords[i] = i + i/10.;

  numMesh1 = 0;
  numMesh2 = 0;
  for (i = 0; i < numElemType; i++) { // loop over elem types
    for (j = 1; j <= 3; j++) { // loop over space dimensions
      // skip if it doesn't make sense to make this elem type at this dim
      if (topoDimPerType[i] > j)
	continue;

      // create a mesh & set coords & conns
      rc = fc_createMesh(dataset1, name, &mesh);
      fail_unless(rc == FC_SUCCESS, "abort: failed to created mesh");
      rc = fc_setMeshCoords(mesh, j, numVertPerType[i], coords);
      fail_unless(rc == FC_SUCCESS, "failed to set vertex coords");
      rc = fc_setMeshElementConns(mesh, elemTypes[i], numElem, conns[i]);
      fail_unless(rc == FC_SUCCESS, "failed to set element conns");
      numMesh1++;

      // fc_copyMesh() to same dataset & check it
      rc = fc_copyMesh(mesh, dataset1, copy_name, 1, 1, 1, 1, &copy_mesh);
      fail_unless(rc == FC_SUCCESS, "failed to copy to same dataset");
      numMesh1++;
      rc = fc_getMeshName(copy_mesh, &temp_name);
      fail_unless(rc == FC_SUCCESS, "failed to get name of copy");
      fail_unless(!strcmp(temp_name, copy_name), "mismatch of name");
      free(temp_name);
      rc = fc_getMeshInfo(copy_mesh, &temp_topoDim, &temp_dim, &temp_numVert,
			  &temp_numElem, &temp_elemType);
      fail_unless(rc == FC_SUCCESS, "failed to get info");
      fail_unless(temp_topoDim == topoDimPerType[i], "mismatch of topoDim");
      fail_unless(temp_dim == j, "mismatch of dim");
      fail_unless(temp_numVert == numVertPerType[i], "mistmatch of numVert");
      fail_unless(temp_numElem == numElem, "mismatch of numElem");
      fail_unless(temp_elemType == elemTypes[i], "mismatch of elemType");
      temp_coords = NULL;
      rc = fc_getMeshCoordsPtr(copy_mesh, &temp_coords);
      fail_unless(rc == FC_SUCCESS, "failed to get coords");
      fail_unless(!memcmp(temp_coords, coords,
		  numVertPerType[i]*sizeof(double)), "mismatch of coords");
      temp_conns = NULL;
      rc = fc_getMeshElementConnsPtr(copy_mesh, &temp_conns);
      fail_unless(rc == FC_SUCCESS, "failed to get conns");
      fail_unless(!memcmp(temp_conns, conns[i],
		  numElem*numVertPerElem[i]*sizeof(int)), "mismatch of conns");

      // fc_copyMesh() to same different dataset & check it
      rc = fc_copyMesh(mesh, dataset2, copy_name, 1, 1, 1, 1, &copy_mesh);
      fail_unless(rc == FC_SUCCESS, "failed to copy to same dataset");
      numMesh2++;
      rc = fc_getMeshName(copy_mesh, &temp_name);
      fail_unless(rc == FC_SUCCESS, "failed to get name of copy");
      fail_unless(!strcmp(temp_name, copy_name), "mismatch of name");
      free(temp_name);
      rc = fc_getMeshInfo(copy_mesh, &temp_topoDim, &temp_dim, &temp_numVert,
			  &temp_numElem, &temp_elemType);
      fail_unless(rc == FC_SUCCESS, "failed to get info");
      fail_unless(temp_topoDim == topoDimPerType[i], "mismatch of topoDim");
      fail_unless(temp_dim == j, "mismatch of dim");
      fail_unless(temp_numVert == numVertPerType[i], "mistmatch of numVert");
      fail_unless(temp_numElem == numElem, "mismatch of numElem");
      fail_unless(temp_elemType == elemTypes[i], "mismatch of elemType");
      temp_coords = NULL;
      rc = fc_getMeshCoordsPtr(copy_mesh, &temp_coords);
      fail_unless(rc == FC_SUCCESS, "failed to get coords");
      fail_unless(!memcmp(temp_coords, coords,
		  numVertPerType[i]*sizeof(double)), "mismatch of coords");
      temp_conns = NULL;
      rc = fc_getMeshElementConnsPtr(copy_mesh, &temp_conns);
      fail_unless(rc == FC_SUCCESS, "failed to get conns");
      fail_unless(!memcmp(temp_conns, conns[i],
		  numElem*numVertPerElem[i]*sizeof(int)), "mismatch of conns");
 
      // --- special cases
      
      // copy with NULL name will use original's name
      rc = fc_copyMesh(mesh, dataset1, NULL, 1, 1, 1, 1, &copy_mesh);
      fail_unless(rc == FC_SUCCESS, "should fail to copy bad mesh");
      rc = fc_getMeshName(copy_mesh, &temp_name);
      fail_unless(rc == FC_SUCCESS, "failed to get name of copy");
      fail_unless(!strcmp(temp_name, name), "mismatch of name");
      free(temp_name);
      rc = fc_deleteMesh(copy_mesh);
      fail_unless(rc == FC_SUCCESS, "failed to delete copied mesh");
  
      // ---- check errors
      
      // fc_copyMesh() -- bad args
      rc = fc_copyMesh(badMesh, dataset1, copy_name, 1, 1, 1, 1,
		       &copy_mesh);
      fail_unless(rc != FC_SUCCESS, "should fail to copy bad mesh");
      fail_unless(FC_HANDLE_EQUIV(copy_mesh, FC_NULL_MESH),
		  "fail should return NULL");
      rc = fc_copyMesh(mesh, FC_NULL_DATASET, copy_name, 1, 1, 1, 1, 
		       &copy_mesh);
      fail_unless(rc != FC_SUCCESS, "should fail to copy to bad dataset");
      fail_unless(FC_HANDLE_EQUIV(copy_mesh, FC_NULL_MESH),
		  "fail should return NULL");
      rc = fc_copyMesh(mesh, dataset1, copy_name, 1, 1, 1, 1, NULL);
      fail_unless(rc != FC_SUCCESS, "should fail to copy if NULL handle");
    }
  }

  // a little more checking & cleanup
  rc = fc_getNumMesh(dataset2, &temp_numMesh);
  fail_unless(rc == FC_SUCCESS, "shouldn't fail to get mesh");
  fail_unless(temp_numMesh == numMesh2, "mismatch of numMesh");
  rc = fc_deleteDataset(dataset2);
  fail_unless(rc == FC_SUCCESS, "failed to delete local dataset");
  rc = fc_getNumMesh(dataset1, &temp_numMesh);
  fail_unless(rc == FC_SUCCESS, "shouldn't fail to get mesh");
  fail_unless(temp_numMesh == numMesh1, "mismatch of numMesh");

  // delete the dataset
  rc = fc_deleteDataset(dataset1);
  fail_unless(rc == FC_SUCCESS, "failed to delete dataset");

  //----- now check the flags but only on a hex mesh-----
  //since these are simple copy calls for each type, make this as simple as possible
  rc = fc_createDataset("dataset1", &dataset1);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create dataset");
  rc = fc_createDataset("dataset2", &dataset2);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create dataset");

  rc = fc_createSimpleHexMesh(dataset1, "hex mesh", 3, 3, 3, lowers, uppers, &mesh);
  fail_unless(rc == FC_SUCCESS, "failed to create mesh");
  rc = fc_getMeshInfo(mesh, &temp_topoDim, &temp_dim, &temp_numVert, &temp_numElem, &temp_elemType);
  fail_unless(rc == FC_SUCCESS, "failed to get mesh info");

  //give it some var and subsets
  rc = fc_createVariable(mesh, "var", &var);
  fail_unless(rc == FC_SUCCESS, "failed to create var");
  rc = fc_createGlobalVariable(dataset1, "global var", &globalvar);
  fail_unless(rc == FC_SUCCESS, "failed to create global var");
  data = (int*)malloc(temp_numElem*sizeof(int));
  for (i = 0; i < temp_numElem; i++)
    data[i] = i;
  rc = fc_setVariableData(var, temp_numElem, 1, FC_AT_ELEMENT, FC_MT_SCALAR,
			  FC_DT_INT, (void*)data);
  fail_unless(rc == FC_SUCCESS, "failed to set variable data");
  rc = fc_setVariableData(globalvar, 1, 1, FC_AT_WHOLE_DATASET, FC_MT_SCALAR,
			  FC_DT_INT, (void*)data);
  fail_unless(rc == FC_SUCCESS, "failed to set glb variable data");
  free(data);

  numStep = 3;
  rc = fc_createSequence(dataset1, "seq", &seq);
  fail_unless(rc == FC_SUCCESS, "failed to create seq");
  data = (int*)malloc(numStep*sizeof(int));
  for (i = 0; i < numStep; i++)
    data[i] = i;
  rc = fc_setSequenceCoords(seq, numStep, FC_DT_INT, (void*)(data));
  fail_unless(rc == FC_SUCCESS, "failed to set seq coords");
  free(data);
  rc = fc_createSeqVariable(mesh, seq, "seqvar1", &numStep, &seqvar1);
  fail_unless(rc == FC_SUCCESS, "failed to create seqvar");
  rc = fc_createSeqVariable(mesh, seq, "seqvar2", &numStep, &seqvar2);
  fail_unless(rc == FC_SUCCESS, "failed to create seqvar");
  rc = fc_createGlobalSeqVariable(dataset1, seq, "global seqvar", &numStep, &globalseqvar);
  fail_unless(rc == FC_SUCCESS, "failed to create seqvar");
  data = (int*)malloc(temp_numVert*sizeof(int));
  for (i = 0; i < temp_numVert; i++)
    data[i] = i;
  for (i = 0; i < numStep; i++){
    rc = fc_setVariableData(seqvar1[i], temp_numVert, 1, FC_AT_VERTEX,
			    FC_MT_SCALAR, FC_DT_INT, (void*)data);
    fail_unless(rc == FC_SUCCESS, "failed to set var data");
    rc = fc_setVariableData(seqvar2[i], temp_numVert, 1, FC_AT_VERTEX,
			    FC_MT_SCALAR, FC_DT_INT, (void*)data);
    fail_unless(rc == FC_SUCCESS, "failed to set var data");
    rc = fc_setVariableData(globalseqvar[i], 1, 1, FC_AT_WHOLE_DATASET,
			    FC_MT_SCALAR, FC_DT_INT, (void*)data);
    fail_unless(rc == FC_SUCCESS, "failed to set var data");
  }
  free(data);

  numMembers = 2;
  rc = fc_createSubset(mesh, "sub", FC_AT_ELEMENT, &sub);
  fail_unless(rc == FC_SUCCESS, "failed to create Subset");
  for (i = 0; i < numMembers; i++){
    rc = fc_addMemberToSubset(sub,i);
    fail_unless(rc == FC_SUCCESS, "failed to add member to subset");
  }
  rc = fc_createSeqSubset(mesh, seq, "seqsub", FC_AT_VERTEX, &numStep, &seqsub);
  fail_unless(rc == FC_SUCCESS, "failed to create seq Subset");
  for (i = 0; i < numStep; i++){
    for (j = 0; j < numMembers; j++){
    rc = fc_addMemberToSubset(seqsub[i],j);
    fail_unless(rc == FC_SUCCESS, "failed to add member to subset");
    }
  }

  // try these in 3 separate calls (geom only and 2 pair)
  // geom copy tested above
  rc = fc_copyMesh(mesh, dataset2, "new mesh", 0, 0, 0, 0, &copy_mesh);
  fail_unless(rc == FC_SUCCESS, "cannot copy mesh");
  rc = fc_getMeshName(copy_mesh, &temp_name);
  fail_unless(rc == FC_SUCCESS, "failed to get name of copy");
  fail_unless(!strcmp(temp_name, "new mesh"), "mismatch of name");
  free(temp_name);

  rc = fc_getNumGlobalVariable(dataset2, &numVar);
  fail_unless(rc == FC_SUCCESS, "failed to get num global var");
  fail_unless(numVar == 0, "wrong num global vars");
  rc = fc_getNumGlobalSeqVariable(dataset2, &numVar);
  fail_unless(rc == FC_SUCCESS, "failed to get num global seqvar");
  fail_unless(numVar == 0, "wrong num global seq vars");
  rc = fc_getNumVariable(copy_mesh, &numVar);
  fail_unless(rc == FC_SUCCESS, "failed to get num var");
  fail_unless(numVar == 0, "wrong num vars");
  rc = fc_getNumSeqVariable(copy_mesh, &numVar);
  fail_unless(rc == FC_SUCCESS, "failed to get num seqvar");
  fail_unless(numVar == 0, "wrong num seq vars");
  rc = fc_getNumSubset(copy_mesh, &numVar);
  fail_unless(rc == FC_SUCCESS, "failed to get num subset");
  fail_unless(numVar == 0, "wrong num subsets");
  rc = fc_getNumSeqSubset(copy_mesh, &numVar);
  fail_unless(rc == FC_SUCCESS, "failed to get num seq subs");
  fail_unless(numVar == 0, "wrong num seq subs");
  fc_deleteMesh(copy_mesh);

  //just check numbers not vals, since for simple cases, it is a simple copy
  rc = fc_copyMesh(mesh, dataset2, "new mesh", 1, 0, 1, 0, &copy_mesh);
  fail_unless(rc == FC_SUCCESS, "cannot copy mesh");
  rc = fc_getMeshName(copy_mesh, &temp_name);
  fail_unless(rc == FC_SUCCESS, "failed to get name of copy");
  fail_unless(!strcmp(temp_name, "new mesh"), "mismatch of name");
  free(temp_name);

  rc = fc_getNumGlobalVariable(dataset2, &numVar);
  fail_unless(rc == FC_SUCCESS, "failed to get num global var");
  fail_unless(numVar == 0, "wrong num global vars");
  rc = fc_getNumGlobalSeqVariable(dataset2, &numVar);
  fail_unless(rc == FC_SUCCESS, "failed to get num global seqvar");
  fail_unless(numVar == 0, "wrong num global seq vars");
  rc = fc_getNumVariable(copy_mesh, &numVar);
  fail_unless(rc == FC_SUCCESS, "failed to get num var");
  fail_unless(numVar == 1, "wrong num vars");
  rc = fc_getNumSeqVariable(copy_mesh, &numVar);
  fail_unless(rc == FC_SUCCESS, "failed to get num seqvar");
  fail_unless(numVar == 0, "wrong num seq vars");
  rc = fc_getNumSubset(copy_mesh, &numVar);
  fail_unless(rc == FC_SUCCESS, "failed to get num subset");
  fail_unless(numVar == 1, "wrong num subsets");
  rc = fc_getNumSeqSubset(copy_mesh, &numVar);
  fail_unless(rc == FC_SUCCESS, "failed to get num seq subs");
  fail_unless(numVar == 0, "wrong num seq subs");
  fc_deleteMesh(copy_mesh);

  //for sequence cases test the sequences and that the vars and subsets
  //are on a matching sequence
  rc = fc_copyMesh(mesh, dataset2, "new mesh", 0, 1, 0, 1, &copy_mesh);
  fail_unless(rc == FC_SUCCESS, "cannot copy mesh");
  rc = fc_getMeshName(copy_mesh, &temp_name);
  fail_unless(rc == FC_SUCCESS, "failed to get name of copy");
  fail_unless(!strcmp(temp_name, "new mesh"), "mismatch of name");
  free(temp_name);

  rc = fc_getNumGlobalVariable(dataset2, &numVar);
  fail_unless(rc == FC_SUCCESS, "failed to get num global var");
  fail_unless(numVar == 0, "wrong num global vars");
  rc = fc_getNumGlobalSeqVariable(dataset2, &numVar);
  fail_unless(rc == FC_SUCCESS, "failed to get num global seqvar");
  fail_unless(numVar == 0, "wrong num global seq vars");
  rc = fc_getNumVariable(copy_mesh, &numVar);
  fail_unless(rc == FC_SUCCESS, "failed to get num var");
  fail_unless(numVar == 0, "wrong num vars");
  rc = fc_getNumSeqVariable(copy_mesh, &numVar);
  fail_unless(rc == FC_SUCCESS, "failed to get num seqvar");
  fail_unless(numVar == 2, "wrong num seq vars");
  rc = fc_getSeqVariableByName(copy_mesh, "seqvar1", &numreturn,
				    &numstepperseq, &returnseqvars);
  fail_unless(rc == FC_SUCCESS, "can't get seq var by name");
  fail_unless(numreturn == 1, "wrong num matching var");
  rc = fc_getSequenceFromSeqVariable(numStep, returnseqvars[0], &retseq1);
  fail_unless(rc == FC_SUCCESS, "cant get seq");
  free(numstepperseq);
  free(returnseqvars[0]);
  free(returnseqvars);
  rc = fc_getSeqVariableByName(copy_mesh, "seqvar2", &numreturn,
				    &numstepperseq, &returnseqvars);
  fail_unless(rc == FC_SUCCESS, "can't get seq var by name");
  fail_unless(numreturn == 1, "wrong num matching var");
  rc = fc_getSequenceFromSeqVariable(numStep, returnseqvars[0], &retseq2);
  fail_unless(rc == FC_SUCCESS, "cant get seq");
  free(numstepperseq);
  free(returnseqvars[0]);
  free(returnseqvars);
  fail_unless(FC_HANDLE_EQUIV(retseq1, retseq2), "wrong ret seqs");
  rc = fc_getSequenceInfo(retseq1,&numreturn, &returndatatype);
  fail_unless(rc == FC_SUCCESS, "cant get seq info");
  fail_unless(returndatatype == FC_DT_INT, "wrong data type");
  fail_unless(numreturn == numStep, "wrong num step");
  //check known seq vals
  rc = fc_getSequenceCoordsPtr(retseq1, (void*) &data);
  fail_unless(rc == FC_SUCCESS, "can get seq coords");
  for (i = 0; i < numStep; i++){
    fail_unless(((int*)(data))[i] == i , "wrong data val");
  }
  rc = fc_getNumSubset(copy_mesh, &numVar);
  fail_unless(rc == FC_SUCCESS, "failed to get num subset");
  fail_unless(numVar == 0, "wrong num subsets");
  rc = fc_getNumSeqSubset(copy_mesh, &numVar);
  fail_unless(rc == FC_SUCCESS, "failed to get num seq subs");
  fail_unless(numVar == 1, "wrong num seq subs");
  rc = fc_getSeqSubsetByName(copy_mesh, "seqsub", &numreturn,
				    &numstepperseq, &returnseqsubs);
  fail_unless(rc == FC_SUCCESS, "can't get seq sub by name");
  fail_unless(numreturn == 1, "wrong num matching seq sub");
  rc = fc_getSequenceFromSeqSubset(numStep, returnseqsubs[0], &retseq1);
  fail_unless(rc == FC_SUCCESS, "cant get seq");
  free(numstepperseq);
  free(returnseqsubs[0]);
  free(returnseqsubs);
  rc = fc_getSequenceInfo(retseq1,&numreturn, &returndatatype);
  fail_unless(rc == FC_SUCCESS, "cant get seq info");
  fail_unless(returndatatype == FC_DT_INT, "wrong data type");
  fail_unless(numreturn == numStep, "wrong num step");
  //check known seq vals
  rc = fc_getSequenceCoordsPtr(retseq1, (void*) &data);
  fail_unless(rc == FC_SUCCESS, "can get seq coords");
  for (i = 0; i < numStep; i++){
    fail_unless(((int*)(data))[i] == i , "wrong data val");
  }
  fc_deleteMesh(copy_mesh);

  //one last test to same dataset - just make sure dont recopy seq
  rc = fc_copyMesh(mesh, dataset1, "copy mesh", 1, 1, 1, 1, &copy_mesh);
  fail_unless(rc == FC_SUCCESS, "cannot copy mesh");
  rc = fc_getMeshName(copy_mesh, &temp_name);
  fail_unless(rc == FC_SUCCESS, "failed to get name of copy");
  fail_unless(!strcmp(temp_name, "copy mesh"), "mismatch of name");
  free(temp_name);

  rc = fc_getNumGlobalVariable(dataset1, &numVar);
  fail_unless(rc == FC_SUCCESS, "failed to get num global var");
  fail_unless(numVar == 1, "wrong num global vars");
  rc = fc_getNumGlobalSeqVariable(dataset1, &numVar);
  fail_unless(rc == FC_SUCCESS, "failed to get num global seqvar");
  fail_unless(numVar == 1, "wrong num global seq vars");
  rc = fc_getNumVariable(copy_mesh, &numVar);
  fail_unless(rc == FC_SUCCESS, "failed to get num var");
  fail_unless(numVar == 1, "wrong num vars");
  rc = fc_getNumSeqVariable(copy_mesh, &numVar);
  fail_unless(rc == FC_SUCCESS, "failed to get num seqvar");
  fail_unless(numVar == 2, "wrong num seq vars");
  rc = fc_getNumSubset(copy_mesh, &numVar);
  fail_unless(rc == FC_SUCCESS, "failed to get num subset");
  fail_unless(numVar == 1, "wrong num subsets");
  rc = fc_getNumSeqSubset(copy_mesh, &numVar);
  fail_unless(rc == FC_SUCCESS, "failed to get num seq subs");
  fail_unless(numVar == 1, "wrong num seq subs");
  rc = fc_getNumSequence(dataset1, &numVar);
  fail_unless(rc == FC_SUCCESS, "can't get num seq"); 
  fail_unless(numVar == 1, "wrong num seq");
  fc_deleteMesh(copy_mesh);


  //bad args - only the flags
  rc = fc_copyMesh(mesh, dataset1, "new_mesh", -1, 1, 1, 1, &copy_mesh);
  fail_unless(rc != FC_SUCCESS, "Should fail for bad flag");
  fail_unless(FC_HANDLE_EQUIV(copy_mesh, FC_NULL_MESH),"bad arg should return null mesh");
  rc = fc_copyMesh(mesh, dataset1, "new_mesh", 2, 1, 1, 1, &copy_mesh);
  fail_unless(rc != FC_SUCCESS, "Should fail for bad flag");
  fail_unless(FC_HANDLE_EQUIV(copy_mesh, FC_NULL_MESH),"bad arg should return null mesh");
  rc = fc_copyMesh(mesh, dataset1, "new_mesh", 0, -1, 0, 0, &copy_mesh);
  fail_unless(rc != FC_SUCCESS, "Should fail for bad flag");
  fail_unless(FC_HANDLE_EQUIV(copy_mesh, FC_NULL_MESH),"bad arg should return null mesh");
  rc = fc_copyMesh(mesh, dataset1, "new_mesh", 0, 2, 0, 0, &copy_mesh);
  fail_unless(rc != FC_SUCCESS, "Should fail for bad flag");
  fail_unless(FC_HANDLE_EQUIV(copy_mesh, FC_NULL_MESH),"bad arg should return null mesh");
  rc = fc_copyMesh(mesh, dataset1, "new_mesh", 0, 0, -1, 0, &copy_mesh);
  fail_unless(rc != FC_SUCCESS, "Should fail for bad flag");
  fail_unless(FC_HANDLE_EQUIV(copy_mesh, FC_NULL_MESH),"bad arg should return null mesh");
  rc = fc_copyMesh(mesh, dataset1, "new_mesh", 0, 0, 2, 0, &copy_mesh);
  fail_unless(rc != FC_SUCCESS, "Should fail for bad flag");
  fail_unless(FC_HANDLE_EQUIV(copy_mesh, FC_NULL_MESH),"bad arg should return null mesh");
  rc = fc_copyMesh(mesh, dataset1, "new_mesh", 0, 0, 0, -1, &copy_mesh);
  fail_unless(rc != FC_SUCCESS, "Should fail for bad flag");
  fail_unless(FC_HANDLE_EQUIV(copy_mesh, FC_NULL_MESH),"bad arg should return null mesh");
  rc = fc_copyMesh(mesh, dataset1, "new_mesh", 0, 0, 0, 2, &copy_mesh);
  fail_unless(rc != FC_SUCCESS, "Should fail for bad flag");
  fail_unless(FC_HANDLE_EQUIV(copy_mesh, FC_NULL_MESH),"bad arg should return null mesh");

  // clean up
  free(seqsub);
  free(seqvar1);
  free(seqvar2);
  free(globalseqvar);
  rc = fc_deleteDataset(dataset1);
  fail_unless(rc == FC_SUCCESS, "failed to delete dataset at end of test");
  rc = fc_deleteDataset(dataset2);
  fail_unless(rc == FC_SUCCESS, "failed to delete dataset at end of test");
}
END_TEST

// test error input for fc_printMesh()
START_TEST(print)
{
  FC_ReturnCode rc;
  FC_Mesh badMesh = { 999, 999 };

  // fc_printMesh()
  rc = fc_printMesh(badMesh, "Bad Mesh!", 1, 1, 1, 1);
  fail_unless(rc != FC_SUCCESS, "should fail if bad mesh");
}
END_TEST 

// release mesh
// don't bother testing that release propagate to child variables
// because we test that in dataset interface
START_TEST(release_mesh)
{
  FC_ReturnCode rc;
  int i;
  FC_Dataset dataset;
  FC_Sequence sequence;
  FC_Mesh mesh1, mesh2, badMesh = { 999, 999 };
  FC_Mesh *returnMeshes; 
  _FC_MeshSlot* meshSlot, *meshSlot2;
  int numStep = 1, time_coords[1] = { 0 };
  double *coords_p;
  int* conns_p;
  int numReturnMeshes;

  // setup
  fc_loadDataset("../data/gen_hex.ex2", &dataset);
  rc = fc_getMeshByName(dataset, "hex mesh", &numReturnMeshes,&returnMeshes);
  fail_unless(rc == FC_SUCCESS, "aborted - cant get mesh by name");
  fail_unless(numReturnMeshes == 1, "failed to find unique mesh by name");
  mesh1 = returnMeshes[0];
  free(returnMeshes);
  meshSlot = _fc_getMeshSlot(mesh1);
  fail_unless(meshSlot != NULL, "abort: failed to get meshslot");
  fc_createSequence(dataset, "temp_sequence", &sequence);
  rc = fc_setSequenceCoords(sequence, numStep, FC_DT_INT, time_coords);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create sequence");

  // --- special cases
  
  // releasing a null handle does not cause an error
  rc = fc_releaseMesh(FC_NULL_MESH);
  fail_unless(rc == FC_SUCCESS, "NULL handle should not fail");
  
  // --- errors
  
  // fc_releaseMesh() --- bad args
  rc = fc_releaseMesh(badMesh);
  fail_unless(rc != FC_SUCCESS, "bad mesh should fail");
  
  // --- coords & elem conns
  
  // test that lazily loaded data is not there yet
  fail_unless(meshSlot->coords == NULL, "abort: should start with no coords");
  fail_unless(meshSlot->elemToVertConns == NULL, 
	      "abort: should start with no conns");
  
  // load & unload mesh coords & conns
  rc = fc_getMeshCoordsPtr(mesh1, &coords_p); // force loading
  fail_unless(rc == FC_SUCCESS, "failed to get coords");
  rc = fc_getMeshElementConnsPtr(mesh1, &conns_p);
  fail_unless(coords_p == meshSlot->coords, "mismatch of coords pointers");
  fail_unless(conns_p == meshSlot->elemToVertConns, 
	      "mismatch of conns pointers");
  rc = fc_releaseMesh(mesh1);
  fail_unless(rc == FC_SUCCESS, "failed to release data");
  fail_unless(meshSlot->coords == NULL && meshSlot->elemToVertConns == NULL, 
	      "coords & conns should now be null"); 

  // --- make copy (uncommitted) mesh
  rc = fc_copyMesh(mesh1, dataset, "new_mesh", 1, 1, 1, 1, &mesh2);
  fail_unless(rc == FC_SUCCESS, "failed to copy mesh");
  meshSlot2 = _fc_getMeshSlot(mesh2);

  // uncommitted mesh's coords & conns shouldn't get deleted
  fail_unless(meshSlot2 != NULL, "abort: failed to get meshslot2");
  fail_unless(meshSlot2->coords != NULL && meshSlot2->elemToVertConns != NULL,
	      "should start with data");
  rc = fc_releaseMesh(mesh2);
  fail_unless(rc == FC_SUCCESS, "failed to release data");
  fail_unless(meshSlot2->coords != NULL && meshSlot2->elemToVertConns != NULL, 
	      "should not delete temp mesh data");

  // --- edges, faces & other big data

  // Behavior should be the same on committed or uncommitted mesh
  for (i = 0; i < 2; i++) { // committed or uncommitted
    FC_Mesh mesh;
    FC_Subset edgeSubset, faceSubset, *edgeSeqSub, *faceSeqSub;
    FC_Variable edgeVar, *edgeSeqVar, faceVar, *faceSeqVar;
    int numEdge, numFace, *edgeData, *faceData;
    int temp_int, *temp_int_p, **temp_int_p_p;
    
    // choose the mesh
    if (i == 0) 
      mesh = mesh1;
    else
      mesh = mesh2;
    meshSlot = _fc_getMeshSlot(mesh);
    
    // --- Check release w/o any subsets or variables
    
    // create conns, parents & neighbors & check
    fc_getMeshEdgeConnsPtr(mesh, &temp_int_p); // force build of edges
    fc_getMeshFaceConnsPtr(mesh, NULL, &temp_int, &temp_int_p); // force build of faces
    fail_unless(meshSlot->numEdge > -1, "numEdge should be > -1");
    fail_unless(meshSlot->edgeToVertConns != NULL, "should have edge conns");
    fail_unless(meshSlot->elemToEdgeConns != NULL, "should have elemToEdgeConns");
    fail_unless(meshSlot->numFace > -1, "numFace should be > -1");
    fail_unless(meshSlot->faceTypes != NULL, "should have facetypes");
    fail_unless(meshSlot->numVertPerFace != NULL, "should have numVertPerFace");
    fail_unless(meshSlot->faceToVertConns != NULL, "should have face conns");
    fail_unless(meshSlot->elemToFaceConns != NULL, "should have elemToFaceConns");
    fail_unless(meshSlot->elemFaceOrients != NULL, "should have elemFaceOrients");
    
    // release -- edges & faces should go away
    rc = fc_releaseMesh(mesh);
    fail_unless(rc == FC_SUCCESS, "failed to release data");
    fail_unless(meshSlot->numEdge == -1, "should have numEdge = -1");
    fail_unless(meshSlot->edgeToVertConns == NULL, "should have no edge conns");
    fail_unless(meshSlot->elemToEdgeConns == NULL, "should have no elemToEdgeConns");
    fail_unless(meshSlot->numFace == -1, "should have numFace = -1");
    fail_unless(meshSlot->faceTypes == NULL, "should have no facetypes");
    fail_unless(meshSlot->numVertPerFace == NULL, "should have no numVertPerFace");
    fail_unless(meshSlot->faceToVertConns == NULL, "should have no face conns");
    fail_unless(meshSlot->elemToFaceConns == NULL, "should have no elemToFaceConns");
    fail_unless(meshSlot->elemFaceOrients == NULL, "should have no elemFaceOrients");
    
    // --- check release w/ various combinations of subsets
    
    // add edge & face subsets
    rc = fc_createSubset(mesh, "temp subset", FC_AT_EDGE, &edgeSubset);
    fail_unless(rc == FC_SUCCESS, "abort: failed to create edge subset");
    rc = fc_createSubset(mesh, "temp subset", FC_AT_FACE, &faceSubset);
    fail_unless(rc == FC_SUCCESS, "abort: failed to create face subset");
    
    // Have both subsets, everything should exist
    rc = fc_releaseMesh(mesh);
    fail_unless(rc == FC_SUCCESS, "failed to release data");
    fail_unless(meshSlot->numEdge > -1, "numEdge should be > -1");
    fail_unless(meshSlot->edgeToVertConns != NULL, "should have edge conns");
    fail_unless(meshSlot->elemToEdgeConns != NULL, "should have elemToEdgeConns");
    fail_unless(meshSlot->numFace > -1, "numFace should be > -1");
    fail_unless(meshSlot->faceTypes != NULL, "should have facetypes");
    fail_unless(meshSlot->numVertPerFace != NULL, "should have numVertPerFace");
    fail_unless(meshSlot->faceToVertConns != NULL, "should have face conns");
    fail_unless(meshSlot->elemToFaceConns != NULL, "should have elemToFaceConns");
    fail_unless(meshSlot->elemFaceOrients != NULL, "should have elemFaceOrients");
    
    // delete just the edge subset & release
    // egde stuff should go away but face stuff stays
    fc_deleteSubset(edgeSubset);
    rc = fc_releaseMesh(mesh);
    fail_unless(rc == FC_SUCCESS, "failed to release data");
    fail_unless(meshSlot->numEdge == -1, "should have numEdge = -1");
    fail_unless(meshSlot->edgeToVertConns == NULL, "should have no edge conns");
    fail_unless(meshSlot->elemToEdgeConns == NULL, "should have no elemToEdgeConns");
    fail_unless(meshSlot->numFace > -1, "numFace should be > -1");
    fail_unless(meshSlot->faceTypes != NULL, "should have facetypes");
    fail_unless(meshSlot->numVertPerFace != NULL, "should have numVertPerFace");
    fail_unless(meshSlot->faceToVertConns != NULL, "should have face conns");
    fail_unless(meshSlot->elemToFaceConns != NULL, "should have elemToFaceConns");
    fail_unless(meshSlot->elemFaceOrients != NULL, "should have elemFaceOrients");
    
    // recreate edge subset,delete face subset & release
    // face stuff should go away but edge stuff stays
    fc_createSubset(mesh, "temp subset", FC_AT_EDGE, &edgeSubset);
    fc_deleteSubset(faceSubset);
    rc = fc_releaseMesh(mesh);
    fail_unless(rc == FC_SUCCESS, "failed to release data");
    fail_unless(meshSlot->numEdge > -1, "numEdge should be > -1");
    fail_unless(meshSlot->edgeToVertConns != NULL, "should have edge conns");
    fail_unless(meshSlot->elemToEdgeConns != NULL, "should have elemToEdgeConns");
    fail_unless(meshSlot->numFace == -1, "should have numFace = -1");
    fail_unless(meshSlot->faceTypes == NULL, "should have no facetypes");
    fail_unless(meshSlot->numVertPerFace == NULL, "should have no numVertPerFace");
    fail_unless(meshSlot->faceToVertConns == NULL, "should have no face conns");
    fail_unless(meshSlot->elemToFaceConns == NULL, "should have no elemToFaceConns");
    fail_unless(meshSlot->elemFaceOrients == NULL, "should have no elemFaceOrients");
    
    // delete all remaining subsets before going on
    fc_deleteSubset(edgeSubset);

    // --- check release w/ various combinations of seq subsets
    
    // create edge & face seq subsets
    fc_createSeqSubset(mesh, sequence, "temp subset", FC_AT_EDGE, &numStep, &edgeSeqSub);
    fail_unless(rc == FC_SUCCESS, "abort: failed to create edge seq subset");
    fc_createSeqSubset(mesh, sequence, "temp subset", FC_AT_FACE, &numStep, &faceSeqSub);
    fail_unless(rc == FC_SUCCESS, "abort: failed to create face seq variable");

    // Have both subsets, everything should exist
    rc = fc_releaseMesh(mesh);
    fail_unless(rc == FC_SUCCESS, "failed to release data");
    fail_unless(meshSlot->numEdge > -1, "numEdge should be > -1");
    fail_unless(meshSlot->edgeToVertConns != NULL, "should have edge conns");
    fail_unless(meshSlot->elemToEdgeConns != NULL, "should have elemToEdgeConns");
    fail_unless(meshSlot->numFace > -1, "numFace should be > -1");
    fail_unless(meshSlot->faceTypes != NULL, "should have facetypes");
    fail_unless(meshSlot->numVertPerFace != NULL, "should have numVertPerFace");
    fail_unless(meshSlot->faceToVertConns != NULL, "should have face conns");
    fail_unless(meshSlot->elemToFaceConns != NULL, "should have elemToFaceConns");
    fail_unless(meshSlot->elemFaceOrients != NULL, "should have elemFaceOrients");
    
    // delete just the edge subset & release
    // egde stuff should go away but face stuff stays
    fc_deleteSeqSubset(numStep,edgeSeqSub);
    free(edgeSeqSub);
    rc = fc_releaseMesh(mesh);
    fail_unless(rc == FC_SUCCESS, "failed to release data");
    fail_unless(meshSlot->numEdge == -1, "should have numEdge = -1");
    fail_unless(meshSlot->edgeToVertConns == NULL, "should have no edge conns");
    fail_unless(meshSlot->elemToEdgeConns == NULL, "should have no elemToEdgeConns");
    fail_unless(meshSlot->numFace > -1, "numFace should be > -1");
    fail_unless(meshSlot->faceTypes != NULL, "should have facetypes");
    fail_unless(meshSlot->numVertPerFace != NULL, "should have numVertPerFace");
    fail_unless(meshSlot->faceToVertConns != NULL, "should have face conns");
    fail_unless(meshSlot->elemToFaceConns != NULL, "should have elemToFaceConns");
    fail_unless(meshSlot->elemFaceOrients != NULL, "should have elemFaceOrients");
    
    // recreate edge subset,delete face subset & release
    // face stuff should go away but edge stuff stays
    fc_createSeqSubset(mesh, sequence, "temp subset", FC_AT_EDGE, &numStep, &edgeSeqSub);
    fc_deleteSeqSubset(numStep,faceSeqSub);
    free(faceSeqSub);
    rc = fc_releaseMesh(mesh);
    fail_unless(rc == FC_SUCCESS, "failed to release data");
    fail_unless(meshSlot->numEdge > -1, "numEdge should be > -1");
    fail_unless(meshSlot->edgeToVertConns != NULL, "should have edge conns");
    fail_unless(meshSlot->elemToEdgeConns != NULL, "should have elemToEdgeConns");
    fail_unless(meshSlot->numFace == -1, "should have numFace = -1");
    fail_unless(meshSlot->faceTypes == NULL, "should have no facetypes");
    fail_unless(meshSlot->numVertPerFace == NULL, "should have no numVertPerFace");
    fail_unless(meshSlot->faceToVertConns == NULL, "should have no face conns");
    fail_unless(meshSlot->elemToFaceConns == NULL, "should have no elemToFaceConns");
    fail_unless(meshSlot->elemFaceOrients == NULL, "should have no elemFaceOrients");
    
    // delete all remaining subsets before going on
    fc_deleteSeqSubset(numStep,edgeSeqSub);
    free(edgeSeqSub);
    
    // --- setup some data for vars

    fc_getMeshNumEdge(mesh, &numEdge);
    fc_getMeshNumFace(mesh, &numFace);
    edgeData = calloc(numEdge, sizeof(int));
    faceData = calloc(numFace, sizeof(int));

    // --- check release w/ various combinations of variables
    
    // create edge & face variables
    fc_createVariable(mesh, "temp variable", &edgeVar);
    rc = fc_setVariableData(edgeVar, numEdge, 1, FC_AT_EDGE,
			    FC_MT_SCALAR, FC_DT_INT, edgeData);
    fail_unless(rc == FC_SUCCESS, "abort: failed to create edge variable");
    fc_createVariable(mesh, "temp variable", &faceVar);
    rc = fc_setVariableData(faceVar, numFace, 1, FC_AT_FACE,
			    FC_MT_SCALAR, FC_DT_INT, faceData);
    fail_unless(rc == FC_SUCCESS, "abort: failed to create face variable");
    
    // Have both vars, everything should exist
    rc = fc_releaseMesh(mesh);
    fail_unless(rc == FC_SUCCESS, "failed to release data");
    fail_unless(meshSlot->numEdge > -1, "numEdge should be > -1");
    fail_unless(meshSlot->edgeToVertConns != NULL, "should have edge conns");
    fail_unless(meshSlot->elemToEdgeConns != NULL, "should have elemToEdgeConns");
    fail_unless(meshSlot->numFace > -1, "numFace should be > -1");
    fail_unless(meshSlot->faceTypes != NULL, "should have facetypes");
    fail_unless(meshSlot->numVertPerFace != NULL, "should have numVertPerFace");
    fail_unless(meshSlot->faceToVertConns != NULL, "should have face conns");
    fail_unless(meshSlot->elemToFaceConns != NULL, "should have elemToFaceConns");
    fail_unless(meshSlot->elemFaceOrients != NULL, "should have elemFaceOrients");
    
    // delete just the edge variable & release
    // egde stuff should go away but face stuff stays
    fc_deleteVariable(edgeVar);
    rc = fc_releaseMesh(mesh);
    fail_unless(rc == FC_SUCCESS, "failed to release data");
    fail_unless(meshSlot->numEdge == -1, "should have numEdge = -1");
    fail_unless(meshSlot->edgeToVertConns == NULL, "should have no edge conns");
    fail_unless(meshSlot->elemToEdgeConns == NULL, "should have no elemToEdgeConns");
    fail_unless(meshSlot->numFace > -1, "numFace should be > -1");
    fail_unless(meshSlot->faceTypes != NULL, "should have facetypes");
    fail_unless(meshSlot->numVertPerFace != NULL, "should have numVertPerFace");
    fail_unless(meshSlot->faceToVertConns != NULL, "should have face conns");
    fail_unless(meshSlot->elemToFaceConns != NULL, "should have elemToFaceConns");
    fail_unless(meshSlot->elemFaceOrients != NULL, "should have elemFaceOrients");
    
    // recreate edge var,delete face var & release
    // face stuff should go away but edge stuff stays
    fc_createVariable(mesh, "temp variable", &edgeVar);
    fc_setVariableData(edgeVar, numEdge, 1, FC_AT_EDGE,
		       FC_MT_SCALAR, FC_DT_INT, edgeData);
    fc_deleteVariable(faceVar);
    rc = fc_releaseMesh(mesh);
    fail_unless(rc == FC_SUCCESS, "failed to release data");
    fail_unless(meshSlot->numEdge > -1, "numEdge should be > -1");
    fail_unless(meshSlot->edgeToVertConns != NULL, "should have edge conns");
    fail_unless(meshSlot->elemToEdgeConns != NULL, "should have elemToEdgeConns");
    fail_unless(meshSlot->numFace == -1, "should have numFace = -1");
    fail_unless(meshSlot->faceTypes == NULL, "should have no facetypes");
    fail_unless(meshSlot->numVertPerFace == NULL, "should have no numVertPerFace");
    fail_unless(meshSlot->faceToVertConns == NULL, "should have no face conns");
    fail_unless(meshSlot->elemToFaceConns == NULL, "should have no elemToFaceConns");
    fail_unless(meshSlot->elemFaceOrients == NULL, "should have no elemFaceOrients");
    
    // delete all remaining variables before going on
    fc_deleteVariable(edgeVar);
    
    // --- check release w/ various combinations of seq variables
    
    // create edge & face seq variables
    fc_createSeqVariable(mesh, sequence, "temp variable", &numStep, &edgeSeqVar);
    rc = fc_setVariableData(edgeSeqVar[0], numEdge, 1, FC_AT_EDGE,
			    FC_MT_SCALAR, FC_DT_INT, edgeData);
    fail_unless(rc == FC_SUCCESS, "abort: failed to create edge seq variable");
    fc_createSeqVariable(mesh, sequence, "temp variable", &numStep, &faceSeqVar);
    rc = fc_setVariableData(faceSeqVar[0], numFace, 1, FC_AT_FACE,
			    FC_MT_SCALAR, FC_DT_INT, faceData);
    fail_unless(rc == FC_SUCCESS, "abort: failed to create face seq variable");
    
    // Have both seq vars, everything should exist
    rc = fc_releaseMesh(mesh);
    fail_unless(rc == FC_SUCCESS, "failed to release data");
    fail_unless(meshSlot->numEdge > -1, "numEdge should be > -1");
    fail_unless(meshSlot->edgeToVertConns != NULL, "should have edge conns");
    fail_unless(meshSlot->elemToEdgeConns != NULL, "should have elemToEdgeConns");
    fail_unless(meshSlot->numFace > -1, "numFace should be > -1");
    fail_unless(meshSlot->faceTypes != NULL, "should have facetypes");
    fail_unless(meshSlot->numVertPerFace != NULL, "should have numVertPerFace");
    fail_unless(meshSlot->faceToVertConns != NULL, "should have face conns");
    fail_unless(meshSlot->elemToFaceConns != NULL, "should have elemToFaceConns");
    fail_unless(meshSlot->elemFaceOrients != NULL, "should have elemFaceOrients");
    
    // delete just the edge seq variable & release
    // egde stuff should go away but face stuff stays
    fc_deleteSeqVariable(numStep, edgeSeqVar);
    free(edgeSeqVar);
    rc = fc_releaseMesh(mesh);
    fail_unless(rc == FC_SUCCESS, "failed to release data");
    fail_unless(meshSlot->numEdge == -1, "should have numEdge = -1");
    fail_unless(meshSlot->edgeToVertConns == NULL, "should have no edge conns");
    fail_unless(meshSlot->elemToEdgeConns == NULL, "should have no elemToEdgeConns");
    fail_unless(meshSlot->numFace > -1, "numFace should be > -1");
    fail_unless(meshSlot->faceTypes != NULL, "should have facetypes");
    fail_unless(meshSlot->numVertPerFace != NULL, "should have numVertPerFace");
    fail_unless(meshSlot->faceToVertConns != NULL, "should have face conns");
    fail_unless(meshSlot->elemToFaceConns != NULL, "should have elemToFaceConns");
    fail_unless(meshSlot->elemFaceOrients != NULL, "should have elemFaceOrients");
    
    // recreate edge seq var,delete face seq var & release
    // face stuff should go away but edge stuff stays
    fc_createSeqVariable(mesh, sequence, "temp variable", &numStep, &edgeSeqVar);
    fc_setVariableData(edgeSeqVar[0], numEdge, 1, FC_AT_EDGE,
		       FC_MT_SCALAR, FC_DT_INT, edgeData);
    fc_deleteSeqVariable(numStep, faceSeqVar);
    free(faceSeqVar);
    rc = fc_releaseMesh(mesh);
    fail_unless(rc == FC_SUCCESS, "failed to release data");
    fail_unless(meshSlot->numEdge > -1, "numEdge should be > -1");
    fail_unless(meshSlot->edgeToVertConns != NULL, "should have edge conns");
    fail_unless(meshSlot->elemToEdgeConns != NULL, "should have elemToEdgeConns");
    fail_unless(meshSlot->numFace == -1, "should have numFace = -1");
    fail_unless(meshSlot->faceTypes == NULL, "should have no facetypes");
    fail_unless(meshSlot->numVertPerFace == NULL, "should have no numVertPerFace");
    fail_unless(meshSlot->faceToVertConns == NULL, "should have no face conns");
    fail_unless(meshSlot->elemToFaceConns == NULL, "should have no elemToFaceConns");
    fail_unless(meshSlot->elemFaceOrients == NULL, "should have no elemFaceOrients");
    
    // delete all remaining variables before going on
    fc_deleteSeqVariable(numStep, edgeSeqVar);
    free(edgeSeqVar);

    // --- cleanup from vars

    free(edgeData);
    free(faceData);

    // --- check release of parents
    
    // create parents & check
    _fc_getMeshEdgeParentsOfVerticesPtr(mesh, &temp_int_p, &temp_int_p_p);
    _fc_getMeshFaceParentsOfVerticesPtr(mesh, &temp_int_p, &temp_int_p_p);
    _fc_getMeshElementParentsOfVerticesPtr(mesh, &temp_int_p, &temp_int_p_p);
    _fc_getMeshElementParentsOfEdgesPtr(mesh, &temp_int_p, &temp_int_p_p);
    _fc_getMeshElementParentsOfFacesPtr(mesh, &temp_int_p, &temp_int_p_p);
    fail_unless(meshSlot->numEdgePerVert != NULL &&
		meshSlot->edgeParentsPerVert != NULL &&
		meshSlot->numFacePerVert != NULL &&
		meshSlot->faceParentsPerVert != NULL &&
		meshSlot->numElemPerVert != NULL &&
		meshSlot->elemParentsPerVert != NULL &&
		meshSlot->numElemPerEdge != NULL &&
		meshSlot->elemParentsPerEdge != NULL &&
		meshSlot->numElemPerFace != NULL &&
		meshSlot->elemParentsPerFace != NULL,
		"should have all parents");

    // releasee & check
    rc = fc_releaseMesh(mesh);
    fail_unless(rc == FC_SUCCESS, "failed to release data");
    fail_unless(meshSlot->numEdgePerVert == NULL &&
		meshSlot->edgeParentsPerVert == NULL &&
		meshSlot->numFacePerVert == NULL &&
		meshSlot->faceParentsPerVert == NULL &&
		meshSlot->numElemPerVert == NULL &&
		meshSlot->elemParentsPerVert == NULL &&
		meshSlot->numElemPerEdge == NULL &&
		meshSlot->elemParentsPerEdge == NULL &&
		meshSlot->numElemPerFace == NULL &&
		meshSlot->elemParentsPerFace == NULL,
		"should have no parents");

    // --- check release of neighbors
    
    // create neighbors & check
    _fc_getMeshVertexNeighborsViaEdgePtr(mesh, &temp_int_p, &temp_int_p_p);
    _fc_getMeshElementNeighborsViaVertexPtr(mesh, &temp_int_p, &temp_int_p_p);
    _fc_getMeshElementNeighborsViaEdgePtr(mesh, &temp_int_p, &temp_int_p_p);
    _fc_getMeshElementNeighborsViaFacePtr(mesh, &temp_int_p, &temp_int_p_p);
    fail_unless(meshSlot->numVertNeighsViaEdge != NULL &&
		meshSlot->vertNeighsViaEdge != NULL &&
		meshSlot->numElemNeighsViaVert != NULL &&
		meshSlot->elemNeighsViaVert != NULL &&
		meshSlot->numElemNeighsViaEdge != NULL &&
		meshSlot->elemNeighsViaEdge != NULL &&
		meshSlot->numElemNeighsViaFace != NULL &&
		meshSlot->elemNeighsViaFace != NULL,
		"should have all neighbors");
    rc = fc_releaseMesh(mesh);
    fail_unless(rc == FC_SUCCESS, "failed to release data");
    fail_unless(meshSlot->numVertNeighsViaEdge == NULL &&
		meshSlot->vertNeighsViaEdge == NULL &&
		meshSlot->numElemNeighsViaVert == NULL &&
		meshSlot->elemNeighsViaVert == NULL &&
		meshSlot->numElemNeighsViaEdge == NULL &&
		meshSlot->elemNeighsViaEdge == NULL &&
		meshSlot->numElemNeighsViaFace == NULL &&
		meshSlot->elemNeighsViaFace == NULL,
		"should have no neighbors");
   }

  // cleanup
  fc_deleteDataset(dataset);
}
END_TEST

// test getting the coords as a variable
// The coords variable is treated specially - only 1 official coords variable
// exists and the same one is always gotten by fc_getMeshCoordsAsVariable.
// However, a copy of the Coords variable gives you a normal variable.
// If you make a copy of a mesh with a coords variable, the copy will have
// it's coord variable present.
START_TEST(coords_as_var)
{
  FC_ReturnCode rc;
  int i, j;
  char name[20] = "blue berry";
  FC_Dataset dataset;
  FC_Mesh mesh, mesh2, badMesh = { 999, 999 };
  int numVar;
  FC_Variable coords_var, coords_var2, *variables;
  int numElemType = 8;
  FC_ElementType elemTypes[8] = { FC_ET_POINT,  FC_ET_LINE,  FC_ET_TRI,
		      	          FC_ET_QUAD,   FC_ET_TET,   FC_ET_PYRAMID,
			          FC_ET_PRISM,  FC_ET_HEX };
  int topoDimPerType[8] = { 0, 1, 2, 2, 3, 3, 3, 3 };
  int numElem = 2; // two elements sharing edge/face
  int numVertPerType[8] = { 2, 3, 4, 6, 5, 6, 8, 12 };
  int conns[8][16] = { { 0,    1 },
		       { 0, 1,    1, 2 },
		       { 0, 1, 2,    1, 3, 2 },
		       { 0, 1, 2, 3,    1, 4, 5, 2 },
		       { 0, 1, 2, 3,    0, 2, 1, 4 },
		       { 0, 1, 2, 3, 4,    0, 3, 2, 1, 5 },
		       { 0, 1, 2, 3, 4, 5,    1, 6, 2, 4, 7, 5 },
		       { 0, 1, 2, 3, 4, 5, 6, 7,    1, 8, 9, 2, 5, 10, 11, 6 } };
  double coords[3*12];
  int temp_numDataPoint, temp_numComp;
  double *temp_coords;
  FC_AssociationType temp_assoc;
  FC_MathType temp_mathType;
  FC_DataType temp_dataType;
  
  // make a dataset to play in
  rc = fc_createDataset("temp dataset", &dataset);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create dataset");

  // setup coords (these make no physical sense)
  for (i = 0; i < 3*numVertPerType[numElemType-1]; i++)
    coords[i] = i + i/10.;

  for (i = 0; i < numElemType; i++) { // loop over elem types
    for (j = 1; j <= 3; j++) { // loop over space dimensions
      // skip if it doesn't make sense to make this elem type at this dim
      if (topoDimPerType[i] > j)
	continue;

      // create a mesh & set coords & conns
      rc = fc_createMesh(dataset, name, &mesh);
      fail_unless(rc == FC_SUCCESS, "abort: failed to created mesh");
      rc = fc_setMeshCoords(mesh, j, numVertPerType[i], coords);
      fail_unless(rc == FC_SUCCESS, "failed to set vertex coords");
      rc = fc_setMeshElementConns(mesh, elemTypes[i], numElem, conns[i]);
      fail_unless(rc == FC_SUCCESS, "failed to set element conns");

      // get coords as variable & test
      rc = fc_getMeshCoordsAsVariable(mesh, &coords_var);
      fail_unless(rc == FC_SUCCESS, "failed to get coords as variable");
      rc = fc_getVariableInfo(coords_var, &temp_numDataPoint, &temp_numComp,
			      &temp_assoc, &temp_mathType, &temp_dataType);
      fail_unless(rc == FC_SUCCESS, "failed to get coords_var info");
      fail_unless(temp_numDataPoint == numVertPerType[i], 
		  "mismatch of numDataPoint");
      fail_unless(temp_numComp == j, "mismatch of numComp");
      fail_unless(temp_assoc == FC_AT_VERTEX, "mismatch of assoc");
      if (j == 1)
	fail_unless(temp_mathType == FC_MT_SCALAR, "mismatch of mathtype"); 
      else
	fail_unless(temp_mathType == FC_MT_VECTOR, "mismatch of mathtype"); 
      fail_unless(temp_dataType == FC_DT_DOUBLE, "mismatch of datatype");
      rc = fc_getVariableDataPtr(coords_var, (void**)&temp_coords);
      fail_unless(rc == FC_SUCCESS, "failed to get coords_var data");
      fail_unless(!memcmp(temp_coords, coords,
		  numVertPerType[i]*sizeof(double)), "mismatch of coords");

      // Calling get again should return same variable (and not create new one)
      rc = fc_getMeshCoordsAsVariable(mesh, &coords_var2);
      fail_unless(rc == FC_SUCCESS, "failed to get coords variable twice");
      fail_unless(FC_HANDLE_EQUIV(coords_var, coords_var2),
		  "did not return same variable handle");
      fc_getVariables(mesh, &numVar, &variables);
      fail_unless(numVar == 1, "should only be 1 variable at this point");
      free(variables);

      // But calling copy does create a new var
      rc = fc_copyVariable(coords_var, mesh, "new var", &coords_var2);
      fail_unless(rc == FC_SUCCESS, "failed to copy coords var");
      fail_unless(!FC_HANDLE_EQUIV(coords_var, coords_var2), 
		  "copy should NOT return same handle");
      fc_getVariables(mesh, &numVar, &variables);
      fail_unless(numVar == 2, "should be two variables at this point");
      free(variables);
      fc_deleteVariable(coords_var2);

      // delete & reget
      rc = fc_deleteVariable(coords_var);
      fail_unless(rc == FC_SUCCESS, "failed to delete coords variable");
      fail_unless(!fc_isVariableValid(coords_var), "handle should be invalid");
      fc_getVariables(mesh, &numVar, &variables);
      fail_unless(numVar == 0, "should be no vars");
      rc = fc_getMeshCoordsAsVariable(mesh, &coords_var);
      fail_unless(rc == FC_SUCCESS, 
		  "failed to get coords var after it was deleted");
      rc = fc_getVariableDataPtr(coords_var, (void**)&temp_coords);
      fail_unless(rc == FC_SUCCESS, "failed to get coords_var data");
      fail_unless(!memcmp(temp_coords, coords,
		  numVertPerType[i]*sizeof(double)), "mismatch of coords");

      // copy the entire mesh w/ the coords variable
      rc = fc_copyMesh(mesh, dataset, "new mesh", 1, 1, 1, 1, &mesh2);
      fail_unless(rc == FC_SUCCESS, "failed to copy mesh");
      fc_getVariables(mesh2, &numVar, &variables);
      fail_unless(numVar == 1, "should have the coords var");
      rc = fc_getMeshCoordsAsVariable(mesh2, &coords_var);
      fail_unless(rc == FC_SUCCESS, "should get coords var");
      fail_unless(FC_HANDLE_EQUIV(coords_var, variables[0]),
		  "coords var should be the variable that was there");
      free(variables);
      rc = fc_getMeshCoordsAsVariable(mesh2, &coords_var2);
      fail_unless(FC_HANDLE_EQUIV(coords_var, coords_var2),
		  "should still get same handle");
      fc_getVariables(mesh2, &numVar, &variables);
      fail_unless(numVar == 1, "should only have 1 var");
      free(variables);
      fc_deleteMesh(mesh2);

      // ------ check errors

      // fc_getMeshCoordsAsVariable() -- bad args
      rc = fc_getMeshCoordsAsVariable(badMesh, &coords_var2);
      fail_unless(rc != FC_SUCCESS, "should fail for bad mesh");
      fail_unless(FC_HANDLE_EQUIV(coords_var2, FC_NULL_VARIABLE),
		  "failure should return null variable");
      rc = fc_getMeshCoordsAsVariable(mesh, NULL);
      fail_unless(rc != FC_SUCCESS, "should fail if var is NULL");
    }
  }

  // delete the dataset
  rc = fc_deleteDataset(dataset);
  fail_unless(rc == FC_SUCCESS, "failed to delete dataset at end of test");
}
END_TEST

// test creating a mesh which is a subset of the original mesh
START_TEST(create_subset_mesh)
{
  FC_ReturnCode rc;
  int i, j, k, m, n, p, q;
  char name[20] = "blue berry", subset_mesh_name[20] = "banana nut";
  FC_Dataset datasets[2];
  FC_Mesh mesh, subset_mesh;
  FC_Subset subsets[7], *temp_subsets;
  FC_Variable variables[6], *temp_variables;
  int numSubset, temp_numSubset, numVariable, temp_numVariable;
  int numElemType = 8;
  FC_ElementType elemTypes[8] = { FC_ET_POINT,  FC_ET_LINE,  FC_ET_TRI,
		      	          FC_ET_QUAD,   FC_ET_TET,   FC_ET_PYRAMID,
			          FC_ET_PRISM,  FC_ET_HEX };
  int topoDimPerType[8] = { 0, 1, 2, 2, 3, 3, 3, 3 };
  int numVertPerElem[8] = { 1, 2, 3, 4, 4, 5, 6, 8 };
  int numElem = 2; // two elements sharing edge/face
  int numVertPerType[8] = { 2, 3, 4, 6, 5, 6, 8, 12 };
  int conns[8][16] = { { 0,    1 },
		       { 0, 1,    1, 2 },
		       { 0, 1, 2,    1, 3, 2 },
		       { 0, 1, 2, 3,    1, 4, 5, 2 },
		       { 0, 1, 2, 3,    0, 2, 1, 4 },
		       { 0, 1, 2, 3, 4,    0, 3, 2, 1, 5 },
		       { 0, 1, 2, 3, 4, 5,    1, 6, 2, 4, 7, 5 },
		       { 0, 1, 2, 3, 4, 5, 6, 7,    1, 8, 9, 2, 5, 10, 11, 6 }
		       };
  // two verts that are only in one or the other element
  int extremeVerts[8][2] = { { 0, 1 },
			     { 0, 2 },
			     { 0, 3 },
			     { 0, 5 },
			     { 3, 4 },
			     { 4, 5 },
			     { 0, 7 },
			     { 0, 11 } };
  double coords[3*12], *temp_coords;
  FC_Sequence sequence, *temp_sequences;
  int temp_numSequence, numStep = 7, temp_numStep, temp_numSeqVar, temp_numSeqSub;
  int *temp_numStepPerVar, *temp_numStepPerSub;
  FC_Variable *seqVars[6], **temp_seqVars;
  FC_Subset *seqSubs[7], **temp_seqSubs;
  FC_DataType seq_dataType = FC_DT_FLOAT;
  float seq_coords[7] = { .11, .22, .33, .44, .55, .66, .77 };
  // index by type, then new vertID - gets old vertID
  int vertLookup[8][8] = { { 1 }, 
			   { 1, 2 },
			   { 1, 2, 3 },
			   { 1, 2, 4, 5 },
			   { 0, 1, 2, 4 },
			   { 0, 1, 2, 3, 5 },
			   { 1, 2, 4, 5, 6, 7 },
			   { 1, 2, 5, 6, 8, 9, 10, 11 } };
  int temp_dim, temp_numVertex, temp_numElement, *temp_conns;
  FC_ElementType temp_elemType;
 
  // make two datasets to play in
  rc = fc_createDataset("same dataset", &datasets[0]);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create dataset");
  rc = fc_createDataset("different dataset", &datasets[1]);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create dataset");

  // create seq for sequence variables
  fc_createSequence(datasets[0], "only sequence", &sequence);
  rc = fc_setSequenceCoords(sequence, numStep, seq_dataType, seq_coords);
  fail_unless(rc == FC_SUCCESS, "abort: failed to set seq coords");

  // setup coords (these make no physical sense)
  for (i = 0; i < 3*numVertPerType[numElemType-1]; i++)
    coords[i] = i + i/10.;

  // do it
  for (i = 0; i < numElemType; i++) { // loop over elem types
    for (j = 1; j <= 3; j++) { // loop over space dimensions
     // skip if it doesn't make sense to make this elem type at this dim
      if (topoDimPerType[i] > j)
	continue;

      // create a mesh & set coords & conns
      rc = fc_createMesh(datasets[0], name, &mesh);
      fail_unless(rc == FC_SUCCESS, "abort: failed to created mesh");
      rc = fc_setMeshCoords(mesh, j, numVertPerType[i], coords);
      fail_unless(rc == FC_SUCCESS, "failed to set vertex coords");
      rc = fc_setMeshElementConns(mesh, elemTypes[i], numElem, conns[i]);
      fail_unless(rc == FC_SUCCESS, "failed to set element conns");

      // create variety of subsets
      // - (0) empty vertex subset - to test empty
      // - (1) just first and last vertex - test permissive subset & subsetting
      //     of a subset
      // - (2) whole, (3) vert & (4) elem subsets always get made
      // - (5) edge + (6) face subsets get made if edges or faces exist for 
      //    elem type 
      numSubset = 0; 
      rc = fc_createSubset(mesh, "empty subset", FC_AT_ELEMENT, 
			   &subsets[numSubset]);
      fail_unless(rc == FC_SUCCESS, "abort: failed to create subset");
      numSubset++;
      rc = fc_createSubset(mesh, "1st and last vertices", FC_AT_VERTEX,
			   &subsets[numSubset]);
      fail_unless(rc == FC_SUCCESS, "abort: failed to create subset");
      rc = fc_addMemberToSubset(subsets[numSubset], extremeVerts[i][0]);
      rc = fc_addMemberToSubset(subsets[numSubset], extremeVerts[i][1]);
      numSubset++;
      rc = fc_createSubset(mesh, "whole", FC_AT_WHOLE_MESH, &subsets[numSubset]);
      fail_unless(rc == FC_SUCCESS, "abort: failed to create subset");
      rc = fc_addMemberToSubset(subsets[numSubset], 0);
      numSubset++;
      rc = fc_createSubset(mesh, "verts in 2nd element", FC_AT_VERTEX,
			   &subsets[numSubset]);
      fail_unless(rc == FC_SUCCESS, "abort: failed to create subset");
      rc = fc_addArrayMembersToSubset(subsets[numSubset], numVertPerElem[i],
				      &conns[i][numVertPerElem[i]]);
      numSubset++;
      if (topoDimPerType[i] > 1) {
	int numChild, *childIDs;
	rc = fc_createSubset(mesh, "edges in 2nd element", FC_AT_EDGE,
			     &subsets[numSubset]);
	fail_unless(rc == FC_SUCCESS, "abort: failed to create subset");
 	rc = fc_getMeshEntityChildren(mesh, FC_AT_ELEMENT, 1,
				      FC_AT_EDGE, &numChild, &childIDs);
	rc = fc_addArrayMembersToSubset(subsets[numSubset], numChild, 
					childIDs);
	free(childIDs);
	numSubset++;
      }
      if (topoDimPerType[i] > 2) {
	int numChild, *childIDs;
	rc = fc_createSubset(mesh, "faces in 2nd element", FC_AT_FACE,
			     &subsets[numSubset]);
	fail_unless(rc == FC_SUCCESS, "abort: failed to create subset");
 	rc = fc_getMeshEntityChildren(mesh, FC_AT_ELEMENT, 1,
				      FC_AT_FACE, &numChild, &childIDs);
	rc = fc_addArrayMembersToSubset(subsets[numSubset], numChild,
					childIDs);
	free(childIDs);
	numSubset++;
      }
      rc = fc_createSubset(mesh, "2nd element", FC_AT_ELEMENT,
			   &subsets[numSubset]);
      fail_unless(rc == FC_SUCCESS, "abort: failed to create subset");
      rc = fc_addMemberToSubset(subsets[numSubset], 1);
      numSubset++;
      if (topoDimPerType[i] <= 1)
	fail_unless(numSubset == 5, 
		 "abort: didn't create expected number of subsets");
      else
	fail_unless(numSubset == 5+topoDimPerType[i]-1, 
		    "abort: didn't create expected number of subsets");

      // create seq subsets that are essentially copies of subsets
      for (k = 0; k < numSubset; k++) {
	char *sub_name;
	FC_AssociationType sub_assoc;
	int numMember,*memberIDs;

	// get info & data from basic subset
	rc = fc_getSubsetAssociationType(subsets[k], &sub_assoc);
	fail_unless(rc == FC_SUCCESS, "abort: failed to get subset assoc type");
	rc = fc_getSubsetName(subsets[k], &sub_name);
	fail_unless(rc == FC_SUCCESS, "abort: failed to get subset name");
	rc = fc_createSeqSubset(mesh, sequence, sub_name, sub_assoc,
				&temp_numStep,&seqSubs[k]);
	fail_unless(rc == FC_SUCCESS, "abort: failed to create seq subset");
	fail_unless(temp_numStep == numStep, "abort: mimsatch in numStep"); 
	free(sub_name);
	rc = fc_getSubsetMembersAsArray(subsets[k], &numMember, &memberIDs);
	fail_unless(temp_numStep == numStep, "abort: failed to get subset members"); 
	// write each step
	for (m = 0; m < numStep; m++) {
	  rc = fc_addArrayMembersToSubset(seqSubs[k][m], numMember, memberIDs);
	  fail_unless(rc == FC_SUCCESS, "abort: failed to add members to subset");
	}
	if (memberIDs) free (memberIDs);
	free(seqSubs[k]);
      }

      // create some variables
      // 0) coords as variable
      // 1) copy of coords
      // 2) whole variable (tensor)
      // 3) edge lengths
      // 4) face areas
      // 5) region volumes
      numVariable = 0;
      fc_getMeshCoordsAsVariable(mesh, &variables[numVariable]);
      numVariable++;
      fc_copyVariable(variables[numVariable-1], mesh, "Coords copy", 
		      &variables[numVariable]);
      numVariable++;
      // use arbitrary data
      fc_createVariable(mesh, "whole_data", &variables[numVariable]);
      fc_setVariableData(variables[numVariable], 1, 4, FC_AT_WHOLE_MESH, 
			 FC_MT_TENSOR, FC_DT_INT, numVertPerType);
      numVariable++;
      rc = fc_getEdgeLengths(mesh, &variables[numVariable]);
      if (rc == FC_SUCCESS)
	numVariable++;
      rc = fc_getFaceAreas(mesh, &variables[numVariable]);
      if (rc == FC_SUCCESS)
	numVariable++;
      rc = fc_getRegionVolumes(mesh, &variables[numVariable]);
      if (rc == FC_SUCCESS)
	numVariable++;
      fail_unless(numVariable == 3 + topoDimPerType[i], 
		  "abort: mismatch of numVariable");

      // create seq variables that are essentially copies of variables
      for (k = 0; k < numVariable; k++) {
	char *var_name;
	int numDataPoint, numComp;
	FC_AssociationType var_assoc;
	FC_MathType mathType;
	FC_DataType datatype;
	void *data;

	// get info & data from basic variable
	fc_getVariableName(variables[k], &var_name);
	fc_createSeqVariable(mesh, sequence, var_name, &temp_numStep,
			     &seqVars[k]);
	free(var_name);
	fail_unless(temp_numStep == numStep, "abort: mimsatch in numStep"); 
	fc_getVariableInfo(variables[k], &numDataPoint, &numComp, &var_assoc,
			   &mathType, &datatype);
	fc_getVariableDataPtr(variables[k], &data);

	// write each step
	for (m = 0; m < numStep; m++) {
	  rc = fc_setVariableData(seqVars[k][m], numDataPoint, numComp,
				  var_assoc, mathType, datatype, data);
	  fail_unless(rc == FC_SUCCESS, "abort: failed to set step data");
	}
	free(seqVars[k]);
      }


      // --- finally, the testing!

      // create subset meshes  
      for (k = 0; k < numSubset; k++) { // loop over subsets	
	for (m = 0; m < 2; m++) { // loop over strictness options
	  for (n = 0; n < 2; n++) { // loop over datasets (same or different)
	    //printf("type = %d, j=%d, subset = %d, strictness = %d, dataset=%d\n", i,j,k,m,n);
	    FC_Variable elemmap, vertmap;
	    rc = fc_createSubsetMesh(subsets[k], datasets[n], m, 
				     1, 1, 1, 1,
				     subset_mesh_name, &subset_mesh,
				     &vertmap, &elemmap);
	    fail_unless(rc == FC_SUCCESS, "failed to create subset mesh");
	    
	    // we set things up so only 3 possible outcomes
	    // 1) subsetting returns empty mesh
	    // 2) subsetting returns the entire mesh
	    // 3) subsetting returns just the 2nd element

	    // testing involves querying the created mesh, it's subset
	    //   and it's variables for the different cases
	    // because our seq variables were copies of the basic
	    //   vars, we can general test seq vars on all cases at same time
	    //   (by comparing to basic vars which were already tested)
	    //  -WSK 10/5/04

	    // empty mesh
	    if (k == 0 || (i != 0 && k == 1 && m == 1)) { // empty mesh
	      fail_unless(FC_HANDLE_EQUIV(subset_mesh, FC_NULL_MESH),
			  "empty subset should return a NULL handle");
	      fail_unless(FC_HANDLE_EQUIV(vertmap, FC_NULL_VARIABLE),
			  "empty subset should return a NULL vertmap variable");
	      fail_unless(FC_HANDLE_EQUIV(elemmap, FC_NULL_VARIABLE),
			  "empty subset should return a NULL elemmap variable");
	      fc_deleteVariable(vertmap);
	      fc_deleteVariable(elemmap);
	      continue;
	    }
	    
	    // cases that return the whole mesh
	    if (k == 1 || k == 2 || 
		     (i != 0 && k < numSubset-1 && m == 0) ) {
	      // query the mesh
	      rc = fc_getMeshInfo(subset_mesh, NULL, &temp_dim, &temp_numVertex,
				  &temp_numElement, &temp_elemType);
	      fail_unless(rc == FC_SUCCESS, "failed to get subset mesh info");
	      fail_unless(temp_dim == j && temp_numVertex == numVertPerType[i]
			  && temp_numElement == numElem && 
			  temp_elemType == elemTypes[i], 
			  "mismatch of mesh info");
	      rc = fc_getMeshCoordsPtr(subset_mesh, &temp_coords);
	      fail_unless(rc == FC_SUCCESS, "failed to get subset mesh coords");
	      fail_unless(!memcmp(temp_coords, coords,
				  sizeof(double)*numVertPerType[i]*j), 
			  "mismatch of coords");
	      rc = fc_getMeshElementConnsPtr(subset_mesh, &temp_conns);
	      fail_unless(rc == FC_SUCCESS, "failed to get subset mesh conns");
	      fail_unless(!memcmp(temp_conns, conns[i],
				  sizeof(int)*numElem*numVertPerElem[i]),
			  "mismatch of conns");

	      {
		char* name_orig;
		int numDP_orig, numComp_orig;
		FC_AssociationType assoc_orig;
		FC_MathType mathType_orig;
		FC_DataType dataType_orig;
		int* data_orig;

		fail_unless(!FC_HANDLE_EQUIV(elemmap,FC_NULL_VARIABLE),
			    "elemmap should not be null");
		fc_getVariableName(elemmap, &name_orig);
		fail_unless(strcmp(name_orig,"FC_ELEMENTMAP") == 0,
			    "mismatch of name");
		free(name_orig);
		fc_getVariableInfo(elemmap, &numDP_orig, &numComp_orig,
				   &assoc_orig, &mathType_orig,
				   &dataType_orig);
		fail_unless(dataType_orig == FC_DT_INT, "elemmap should be int");
		fail_unless(mathType_orig == FC_MT_SCALAR, "elemmap should be scalar");
		fail_unless(assoc_orig == FC_AT_ELEMENT, "elemmap should be at elem");
		fail_unless(numComp_orig == 1, "elemmap should have 1 component");
		fail_unless(numDP_orig == temp_numElement, "mismatch of num datapoint");
		fc_getVariableDataPtr(elemmap, (void*)&data_orig);
		for (p = 0; p < temp_numElement; p++){
		  fail_unless(data_orig[p] == p, "mismatch of elems");
		}
		fc_deleteVariable(elemmap);

		fail_unless(!FC_HANDLE_EQUIV(vertmap,FC_NULL_VARIABLE),
			    "vertmmap should not be null");
		fc_getVariableName(vertmap, &name_orig);
		fail_unless(strcmp(name_orig,"FC_VERTEXMAP") == 0,
			    "mismatch of name");
		free(name_orig);
		fc_getVariableInfo(vertmap, &numDP_orig, &numComp_orig,
				   &assoc_orig, &mathType_orig,
				   &dataType_orig);
		fail_unless(dataType_orig == FC_DT_INT, "vertmap should be int");
		fail_unless(mathType_orig == FC_MT_SCALAR, "vertmap should be scalar");
		fail_unless(assoc_orig == FC_AT_VERTEX, "vertmap should be at vertex");
		fail_unless(numComp_orig == 1, "vertmap should have 1 component");
		fail_unless(numDP_orig == temp_numVertex, "mismatch of num datapoint");
		fc_getVariableDataPtr(vertmap, (void*)&data_orig);
		for (p = 0; p < temp_numVertex; p++){
		  fail_unless(data_orig[p] == p, "mismatch of verts");
		}
		fc_deleteVariable(vertmap);
	      }

	      // query the subsets
	      rc = fc_getSubsets(subset_mesh, &temp_numSubset, &temp_subsets);
	      fail_unless(rc == FC_SUCCESS, "failed to get subsets");
	      fail_unless(temp_numSubset == numSubset, "mismatch of numSubset");
	      for (p = 0; p < numSubset; p++) {
		char* name_orig, *name_new;
		int numMember_orig, numMember_new, *members_orig, *members_new;
		int maxNumMember_orig, maxNumMember_new; 
		FC_AssociationType assoc_orig, assoc_new;
		fc_getSubsetName(subsets[p], &name_orig);
		fc_getSubsetName(temp_subsets[p], &name_new);
		fail_unless(!strcmp(name_orig, name_new), "mismatch of name");
		free(name_orig);
		free(name_new);
		fc_getSubsetInfo(subsets[p], &numMember_orig, 
				 &maxNumMember_orig, &assoc_orig);
		fc_getSubsetInfo(temp_subsets[p], &numMember_new,
				 &maxNumMember_new, &assoc_new);
		fail_unless(numMember_orig == numMember_new, 
			    "mismatch of numMember");
		fail_unless(maxNumMember_orig == maxNumMember_new, 
			    "mismatch of maxNumMember");
		fail_unless(assoc_orig == assoc_new, "mismatch of assoc");
		fc_getSubsetMembersAsArray(subsets[p], &numMember_orig,
					   &members_orig);
		fc_getSubsetMembersAsArray(temp_subsets[p], &numMember_new,
					   &members_new);
		fail_unless(!memcmp(members_orig, members_new,
				    numMember_orig*sizeof(int)),
			    "mismatch of members");
		free(members_orig);
		free(members_new);
	      }

	      // query the variables
	      rc = fc_getVariables(subset_mesh, &temp_numVariable,
				   &temp_variables); 
	      fail_unless(rc == FC_SUCCESS, "failed to get variables");
	      fail_unless(temp_numVariable == numVariable, 
			  "mismatch of numVariable");	      
	      for (p = 0; p < numVariable; p++) {
		char* name_orig, *name_new;
		int numDP_orig, numDP_new, numComp_orig, numComp_new;
		FC_AssociationType assoc_orig, assoc_new;
		FC_MathType mathType_orig, mathType_new;
		FC_DataType dataType_orig, dataType_new;
		void* data_orig, *data_new;
		fc_getVariableName(variables[p], &name_orig);
		fc_getVariableName(temp_variables[p], &name_new);
		fail_unless(!strcmp(name_orig, name_new), "mismatch of name");
		free(name_orig);
		free(name_new);
		fc_getVariableInfo(variables[p], &numDP_orig, &numComp_orig,
				   &assoc_orig, &mathType_orig,
				   &dataType_orig);
		fc_getVariableInfo(temp_variables[p], &numDP_new, &numComp_new,
				   &assoc_new, &mathType_new, &dataType_new);
		fail_unless(numDP_orig == numDP_new, 
			    "mismatch of numDataPoint");
		fail_unless(numComp_orig == numComp_new, 
			    "mismatch of numComp");
		fail_unless(assoc_orig == assoc_new, "mismatch of assoc");
		fail_unless(mathType_orig == mathType_new, 
			    "mismatch of mathType");
		fail_unless(dataType_orig == dataType_new,
			    "mismatch of dataType");
		fc_getVariableDataPtr(variables[p], &data_orig);
		fc_getVariableDataPtr(temp_variables[p], &data_new);
		fail_unless(!memcmp(data_orig, data_new, numDP_new * 
       			          numComp_new*fc_sizeofDataType(dataType_new)),
			    "mismatch of data");
	      }
	    } // end of case where whole mesh is copied
	    
	    // remaining cases return the 2nd element
	    else { // the subset
	      {
		char* name_orig;
		int numDP_orig, numComp_orig;
		FC_AssociationType assoc_orig;
		FC_MathType mathType_orig;
		FC_DataType dataType_orig;
		int* data_orig;

		fail_unless(!FC_HANDLE_EQUIV(elemmap,FC_NULL_VARIABLE),
			    "elemmap should not be null");
		fc_getVariableName(elemmap, &name_orig);
		fail_unless(strcmp(name_orig,"FC_ELEMENTMAP") == 0,
			    "mismatch of name");
		free(name_orig);
		fc_getVariableInfo(elemmap, &numDP_orig, &numComp_orig,
				   &assoc_orig, &mathType_orig,
				   &dataType_orig);
		fail_unless(dataType_orig == FC_DT_INT, "elemmap should be int");
		fail_unless(mathType_orig == FC_MT_SCALAR, "elemmap should be scalar");
		fail_unless(assoc_orig == FC_AT_ELEMENT, "elemmap should be at elem");
		fail_unless(numComp_orig == 1, "elemmap should have 1 component");
		fail_unless(numDP_orig == 1, "elemmap should have 1 datapoint");
		fc_getVariableDataPtr(elemmap, (void*)&data_orig);
		fail_unless(data_orig[0] == 1,
			    "elemmap should have the original second element");

		fail_unless(!FC_HANDLE_EQUIV(vertmap,FC_NULL_VARIABLE),
			    "vertmap should not be null");
		fc_getVariableName(vertmap, &name_orig);
		fail_unless(strcmp(name_orig,"FC_VERTEXMAP") == 0,
			    "mismatch of name");
		free(name_orig);
		fc_getVariableInfo(vertmap, &numDP_orig, &numComp_orig,
				   &assoc_orig, &mathType_orig,
				   &dataType_orig);
		fail_unless(dataType_orig == FC_DT_INT, "vertmap should be int");
		fail_unless(mathType_orig == FC_MT_SCALAR, "vertmap should be scalar");
		fail_unless(assoc_orig == FC_AT_VERTEX, "vertmap should be at vertex");
		fail_unless(numComp_orig == 1, "vertmap should have 1 component");
		fail_unless(numDP_orig == numVertPerElem[i],
			    "vertmap should have num datapoint in one elem");
		fc_getVariableDataPtr(vertmap, (void*)&data_orig);
		for (p = 0; p < numVertPerElem[i]; p++){
		  fail_unless(data_orig[p] == vertLookup[i][p],
			      "mismatch in vertmap");
		}
		fc_deleteVariable(elemmap);
		fc_deleteVariable(vertmap);
	      }

	      // query the mesh
	      rc = fc_getMeshInfo(subset_mesh, NULL, &temp_dim, &temp_numVertex,
				  &temp_numElement, &temp_elemType);
	      fail_unless(rc == FC_SUCCESS, "failed to get subset mesh info");
	      fail_unless(temp_dim == j && temp_numVertex == numVertPerElem[i]
			  && temp_numElement == 1 && 
			  temp_elemType == elemTypes[i], 
			  "mismatch of mesh info");
	      rc = fc_getMeshCoordsPtr(subset_mesh, &temp_coords);
	      fail_unless(rc == FC_SUCCESS, "failed to get subset mesh coords");
	      for (p = 0; p < numVertPerElem[i]; p++) 
		fail_unless(!memcmp(&temp_coords[p*j], 
				    &coords[j*vertLookup[i][p]], 
				    j*sizeof(double)), 
			    "mismatch of coords");
	      rc = fc_getMeshElementConnsPtr(subset_mesh, &temp_conns);
	      fail_unless(rc == FC_SUCCESS, "failed to get subset mesh conns");
	      for (p = 0; p < numVertPerElem[i]; p++) 
		fail_unless(vertLookup[i][temp_conns[p]] == 
			    conns[i][numVertPerElem[i]+p],
			    "mismatch of conns");
	      
	      // query the subsets
	      rc = fc_getSubsets(subset_mesh, &temp_numSubset, &temp_subsets);
	      fail_unless(rc == FC_SUCCESS, "failed to get subsets");
	      fail_unless(temp_numSubset == numSubset, "mismatch of numSubset");
	      for (p = 0; p < numSubset; p++) {
		char* name_orig, *name_new;
		int numEntity;
		int numMember_new, maxNumMember_new, *members_new; 
		FC_AssociationType assoc_orig, assoc_new;
		fc_getSubsetName(subsets[p], &name_orig);
		fc_getSubsetName(temp_subsets[p], &name_new);
		fail_unless(!strcmp(name_orig, name_new), "mismatch of name");
		free(name_orig);
		free(name_new);
		fc_getSubsetAssociationType(subsets[p], &assoc_orig);
		fc_getSubsetInfo(temp_subsets[p], NULL, &maxNumMember_new, 
				 &assoc_new);
		fc_getMeshNumEntity(subset_mesh, assoc_orig, &numEntity);
		fail_unless(maxNumMember_new == numEntity, 
			    "mismatch of maxNumMember");
		fail_unless(assoc_orig == assoc_new, "mismatch of assoc");
		fc_getSubsetMembersAsArray(temp_subsets[p], &numMember_new,
					   &members_new);
		if (p == 0) { // first subset is empty
		  fail_unless(numMember_new == 0, "mismatch of numMember");
		  fail_unless(members_new == NULL, "mismatch of members");
		}
		else if (p == 1) {// orig subset had 2 verts but only 1 here
		  fail_unless(numMember_new == 1, "mismatch of numMember");
		  fail_unless(members_new[0] == numVertPerElem[i]-1,
			      "mismatch of members");
		}
		else { // remaining subsets are all entities in only element
		  fail_unless(numMember_new == numEntity, 
			      "mismatch of numMember");
		  for (q = 0; q < numMember_new; q++)
		    fail_unless(members_new[q] == q, "mismatch of members");
		}
		free(members_new);
	      }

	      // query the variables
	      rc = fc_getVariables(subset_mesh, &temp_numVariable,
				   &temp_variables); 
	      fail_unless(rc == FC_SUCCESS, "failed to get variables");
	      fail_unless(temp_numVariable == numVariable, 
			  "mismatch of numVariable");
	      // in loop just check metadata
	      for (p = 0; p < numVariable; p++) {
		char* name_orig, *name_new;
		int numEntity;
		int numDP_new, numComp_orig, numComp_new;
		FC_AssociationType assoc_orig, assoc_new;
		FC_MathType mathType_orig, mathType_new;
		FC_DataType dataType_orig, dataType_new;
		void *data_good, *data_new;
		fc_getVariableName(variables[p], &name_orig);
		fc_getVariableName(temp_variables[p], &name_new);
		fail_unless(!strcmp(name_orig, name_new), "mismatch of name");
		free(name_orig);
		free(name_new);
		fc_getVariableInfo(variables[p], NULL, &numComp_orig,
				   &assoc_orig, &mathType_orig,
				   &dataType_orig);
		fc_getVariableInfo(temp_variables[p], &numDP_new, &numComp_new,
				   &assoc_new, &mathType_new, &dataType_new);
		fc_getMeshNumEntity(subset_mesh, assoc_orig, &numEntity);
		fail_unless(numDP_new == numEntity, "mismatch of numDataPoint");
		fail_unless(numComp_orig == numComp_new, 
			    "mismatch of numComp");
		fail_unless(assoc_orig == assoc_new, "mismatch of assoc");
		fail_unless(mathType_orig == mathType_new, 
			    "mismatch of mathType");
		fail_unless(dataType_orig == dataType_new,
			    "mismatch of dataType");
		// check data
		if (p == 0) { // coords as var
		  FC_Variable coords_var;
		  fc_getMeshCoordsAsVariable(subset_mesh, &coords_var);
		  fail_unless(FC_HANDLE_EQUIV(coords_var, temp_variables[0]),
					    "first var should be coords var");
		  fc_getMeshCoordsPtr(subset_mesh, (double**)&data_good);
		}
		else if (p == 1) { // copy of coords var
		  fc_getMeshCoordsPtr(subset_mesh, (double**)&data_good);
		}
		else if (p == 2) { // whole variable
		  fc_getVariableDataPtr(variables[p], &data_good);
		}
		else if (p == 3) { // edge lengths
		  FC_Variable edge_lengths;
		  fc_getEdgeLengths(subset_mesh, &edge_lengths);
		  fc_getVariableDataPtr(edge_lengths, &data_good);
		}
		else if (p == 4) { // face areas
		  FC_Variable face_areas;
		  fc_getFaceAreas(subset_mesh, &face_areas);
		  fc_getVariableDataPtr(face_areas, &data_good);
		} 
		else if (p == 5) { // region volumes
		  FC_Variable region_volumes;
		  fc_getRegionVolumes(subset_mesh, &region_volumes);
		  fc_getVariableDataPtr(region_volumes, &data_good);
		}
		fc_getVariableDataPtr(temp_variables[p], &data_new);
		fail_unless(!memcmp(data_good, data_new, numEntity * 
			       numComp_new * fc_sizeofDataType(dataType_new)),
			    "mismatch of data");
	      }
	    } // end of case where only part of mesh was copied

	    // For both cases:

	    // query sequences - should always be 1
	    fc_getSequences(datasets[n], &temp_numSequence, &temp_sequences);
	    fail_unless(temp_numSequence == 1, "should always be 1");
	    
	    // query the sequence variables
	    // check against the basic variables which we already checked
	    rc = fc_getSeqVariables(subset_mesh, &temp_numSeqVar,
				    &temp_numStepPerVar, &temp_seqVars);
	    fail_unless(rc == FC_SUCCESS, "failed to get seqVars");
	    fail_unless(temp_numSeqVar == numVariable,
			"mismatch of numSeqVar");
	    for (p = 0; p < numVariable; p++) {
	      FC_Sequence temp_seq;
	      char* name_orig, *name_new;
	      int numDP_orig, numDP_new, numComp_orig, numComp_new;
	      FC_AssociationType assoc_orig, assoc_new;
	      FC_MathType mathType_orig, mathType_new;
	      FC_DataType dataType_orig, dataType_new;
	      void* data_orig, *data_new;
	      
	      fc_getSequenceFromSeqVariable(temp_numStepPerVar[p],
					 temp_seqVars[p], &temp_seq);
	      fail_unless(FC_HANDLE_EQUIV(temp_seq, temp_sequences[0]), 
			  "mismatch of sequence");
	      fail_unless(temp_numStepPerVar[p] == numStep, 
			  "mismatch of numStep");

	      // get "correct" answer from tested basic variables
	      fc_getVariableName(temp_variables[p], &name_orig);
	      fc_getVariableInfo(temp_variables[p], &numDP_orig, &numComp_orig,
				 &assoc_orig, &mathType_orig,
				 &dataType_orig);
	      fc_getVariableDataPtr(temp_variables[p], &data_orig);

	      for (q = 0; q < numStep; q++) {
		fc_getVariableName(temp_seqVars[p][q], &name_new);
		fail_unless(!strcmp(name_orig, name_new), "mismatch of name");
		free(name_new);
		fc_getVariableInfo(temp_seqVars[p][q], &numDP_new, &numComp_new,
				   &assoc_new, &mathType_new, &dataType_new);
		fail_unless(numDP_orig == numDP_new, 
			    "mismatch of numDataPoint");
		fail_unless(numComp_orig == numComp_new, 
			    "mismatch of numComp");
		fail_unless(assoc_orig == assoc_new, "mismatch of assoc");
		fail_unless(mathType_orig == mathType_new, 
			    "mismatch of mathType");
		fail_unless(dataType_orig == dataType_new,
			    "mismatch of dataType");
		fc_getVariableDataPtr(temp_seqVars[p][q], &data_new);
		fail_unless(!memcmp(data_orig, data_new, numDP_new * 
				numComp_new*fc_sizeofDataType(dataType_new)),
			    "mismatch of data");
	      }
	      free(name_orig);
	    }

	    // query the sequence subsets
	    // check against the basic subsets which we already checked
	    rc = fc_getSeqSubsets(subset_mesh, &temp_numSeqSub,
				    &temp_numStepPerSub, &temp_seqSubs);
	    fail_unless(rc == FC_SUCCESS, "failed to get seqSubs");
	    fail_unless(temp_numSeqSub == numSubset,
			"mismatch of numSeqSub");
	    for (p = 0; p < numSubset; p++) {
	      FC_Sequence temp_seq;
	      char* name_orig, *name_new;
	      FC_AssociationType assoc_orig, assoc_new;
	      int numMember_orig, numMember_new, maxNumMember_orig, maxNumMember_new;
	      int *members_orig, *members_new;

	      fc_getSequenceFromSeqSubset(temp_numStepPerSub[p],
					 temp_seqSubs[p], &temp_seq);
	      fail_unless(FC_HANDLE_EQUIV(temp_seq, temp_sequences[0]), 
			  "mismatch of sequence");
	      fail_unless(temp_numStepPerSub[p] = numStep, 
			  "mismatch of numStep");

	      // get "correct" answer from tested basic subsets
	      fc_getSubsetName(temp_subsets[p], &name_orig);
	      fc_getSubsetInfo(temp_subsets[p], &numMember_orig, &maxNumMember_orig,
			       &assoc_orig);
	      fc_getSubsetMembersAsArray(temp_subsets[p], &numMember_orig, &members_orig);
	      for (q = 0; q < numStep; q++) {
		fc_getSubsetName(temp_seqSubs[p][q], &name_new);
		fail_unless(!strcmp(name_orig, name_new), "mismatch of name");
		free(name_new);
		fc_getSubsetInfo(temp_seqSubs[p][q], &numMember_new, &maxNumMember_new,
				 &assoc_new);
		fail_unless(numMember_orig == numMember_new,
			    "mismatch of numMember");
		fail_unless(maxNumMember_orig == maxNumMember_new,
			    "mismatch of maxNumMember");
		fail_unless(assoc_orig == assoc_new, "mismatch of assoc");
		fc_getSubsetMembersAsArray(temp_seqSubs[p][q], &numMember_new, &members_new);
		fail_unless(!memcmp(members_orig, members_new, numMember_new*sizeof(int)),
			    "mismatch of members");
		free(members_new);
	      }
	      free(name_orig);
	      free(members_orig);
	    }

	    free(temp_sequences);
	    free(temp_numStepPerVar);
	    free(temp_numStepPerSub);
	    for (p = 0; p < temp_numSeqVar; p++)
	      free(temp_seqVars[p]);
	    free(temp_seqVars);
	    for (p = 0; p < temp_numSeqSub; p++)
	      free(temp_seqSubs[p]);
	    free(temp_seqSubs);

	    // this came from different places depending on k case
	    free(temp_variables);
	    free(temp_subsets);

	    // done with subset mesh
	    rc = fc_deleteMesh(subset_mesh);
	    fail_unless(rc == FC_SUCCESS, "failed to delete subset mesh");

	    //test no copy options (try 3 cases)
	    rc = fc_createSubsetMesh(subsets[k], datasets[n], m, 
				     0,0,0,0,
				     subset_mesh_name, &subset_mesh,
				     &vertmap, &elemmap);
	    fail_unless(rc == FC_SUCCESS, "failed to create subset mesh");
	    //make sure nothing there, except maps
	    rc = fc_getVariables(subset_mesh, &temp_numVariable,
				 &temp_variables); 
	    fail_unless(temp_numVariable == 2, "wrong num vars");
	    for (p = 0; p < temp_numVariable; p++){
	      char *tempname;
	      rc  = fc_getVariableName(temp_variables[p],&tempname);
	      fail_unless(rc == FC_SUCCESS, "can't get variable name");
	      fail_unless(strcmp(tempname,"FC_ELEMENTMAP") == 0 ||
			  strcmp(tempname,"FC_VERTEXMAP") == 0, 
			  "wrong variable name");
	      free(tempname);
	    }
	    free(temp_variables);
	    rc = fc_getNumSeqVariable(subset_mesh,&temp_numVariable);
	    fail_unless(rc == FC_SUCCESS, "failed to get num seqvars");
	    fail_unless(temp_numVariable == 0, "wrong num of variables");
	    rc = fc_getNumSubset(subset_mesh,&temp_numVariable);
	    fail_unless(rc == FC_SUCCESS, "failed to get num subsets");
	    fail_unless(temp_numVariable == 0, "wrong num of subsets");
	    rc = fc_getNumSeqSubset(subset_mesh,&temp_numSubset);
	    fail_unless(rc == FC_SUCCESS, "failed to get num subsets");
	    fail_unless(temp_numVariable == 0, "wrong num of seq subsets");
	    rc = fc_deleteMesh(subset_mesh);
	    fail_unless(rc == FC_SUCCESS, "failed to delete subset mesh");


	    rc = fc_createSubsetMesh(subsets[k], datasets[n], m, 
				     1,0,0,1,
				     subset_mesh_name, &subset_mesh,
				     &vertmap, &elemmap);
	    fail_unless(rc == FC_SUCCESS, "failed to create subset mesh");
	    rc = fc_getNumVariable(subset_mesh,&temp_numVariable);
	    fail_unless(rc == FC_SUCCESS, "failed to get num vars");
	    fail_unless(temp_numVariable == 3 + topoDimPerType[i]+2, "wrong num of vars"); 
	    rc = fc_getNumSeqVariable(subset_mesh,&temp_numVariable);
	    fail_unless(rc == FC_SUCCESS, "failed to get num seqvars");
	    fail_unless(temp_numVariable == 0, "wrong num of seq vars");
	    rc = fc_getNumSubset(subset_mesh,&temp_numVariable);
	    fail_unless(rc == FC_SUCCESS, "failed to get num subsets");
	    fail_unless(temp_numVariable == 0, "wrong num of subsets");
	    rc = fc_getNumSeqSubset(subset_mesh,&temp_numVariable);
	    fail_unless(rc == FC_SUCCESS, "failed to get num seq subsets");
	    if (topoDimPerType[i] <= 1)
	      fail_unless(temp_numVariable == 5, "wrong number of seq subsets");
	    else
	      fail_unless(temp_numVariable == 5+topoDimPerType[i]-1, 
			  "wrong number of seq subsets");

	    rc = fc_createSubsetMesh(subsets[k], datasets[n], m, 
				     0,1,1,0,
				     subset_mesh_name, &subset_mesh,
				     &vertmap, &elemmap);
	    fail_unless(rc == FC_SUCCESS, "failed to create subset mesh");
	    rc = fc_getNumVariable(subset_mesh,&temp_numVariable);
	    fail_unless(rc == FC_SUCCESS, "failed to get num vars");
	    fail_unless(temp_numVariable == 2, "wrong num of vars");
	    rc = fc_getNumSeqVariable(subset_mesh,&temp_numVariable);
	    fail_unless(rc == FC_SUCCESS, "failed to get num seqvars");
	    fail_unless(temp_numVariable == 3 + topoDimPerType[i], "wrong num of seq vars"); 
	    rc = fc_getNumSubset(subset_mesh,&temp_numVariable);
	    fail_unless(rc == FC_SUCCESS, "failed to get num subsets");
	    if (topoDimPerType[i] <= 1)
	      fail_unless(temp_numVariable == 5, "wrong number of subsets");
	    else
	      fail_unless(temp_numVariable == 5+topoDimPerType[i]-1, 
			  "wrong number of subsets");
	    rc = fc_getNumSeqSubset(subset_mesh,&temp_numVariable);
	    fail_unless(rc == FC_SUCCESS, "failed to get num subsets");
	    fail_unless(temp_numVariable == 0, "wrong num of seq subsets");

	  } // n - dataset
	} // m - strickness
      } // k - subset used to subset

      // done with the mesh
      rc = fc_deleteMesh(mesh);
      fail_unless(rc == FC_SUCCESS, "failed to delete orig mesh");

    } // j - dim
  } // i - mesh type

  // delete the datasets
  rc = fc_deleteDataset(datasets[0]);
  fail_unless(rc == FC_SUCCESS, "failed to delete dataset 0 at end of test");
  rc = fc_deleteDataset(datasets[1]);
  fail_unless(rc == FC_SUCCESS, "failed to delete dataset 1 at end of test");
}
END_TEST

// test fc_createSimpleHexMesh()
START_TEST(simple_hex_mesh)
{
  FC_ReturnCode rc;
  int i, j;
  char name[20] = "blue berry", *temp_name;
  FC_Dataset dataset, bad_dataset = { -99, -99 };
  FC_Mesh mesh;
  int numElemPerDim[3] = { 2, 4, 3 }; // all different to fully test
  int numVertPerDim[3];
  double lowers[3] = { 1, 1.1, 1.01 }, uppers[3] = { 1.99, 1.9, 2.0 };
  double lens[3], temp_coords[3];
  int numDim = 3, numVert, numElem, numVertPerElem = 8, vertID;
  int temp_topoDim, temp_numDim, temp_numVert, temp_numElem;
  FC_ElementType temp_elemType;
  double* coords;
  int* conns;
  int firstElemConns[8] = {  0,  1,  4,  3, 15, 16, 19, 18 };
  int lastElemConns[8] =  { 40, 41, 44, 43, 55, 56, 59, 58 };
  ICoords icoords;
  // adjust icoords (compared to lower vertex) to get icoords of other verts in hex
  int vertAdjust[8][3] = { { 0, 0, 0 },
			   { 1, 0, 0 },
			   { 1, 1, 0 },
			   { 0, 1, 0 },
			   { 0, 0, 1 },
			   { 1, 0, 1 },
			   { 1, 1, 1 },
			   { 0, 1, 1 } };

  // setup
  numVert = 1;
  numElem = 1;
  for (i = 0; i < numDim; i++) {
    numVertPerDim[i] = numElemPerDim[i] + 1;
    numVert *= numVertPerDim[i];
    numElem *= numElemPerDim[i];
    lens[i] = (uppers[i] - lowers[i])/numElemPerDim[i];
  }

  // *** make a dataset to play in
  rc = fc_createDataset("same dataset", &dataset);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create dataset");

  // *** Do it right
  rc = fc_createSimpleHexMesh(dataset, name, numElemPerDim[0], numElemPerDim[1],
			      numElemPerDim[2], lowers, uppers, &mesh);
  fail_unless(rc == FC_SUCCESS, "failed to create simple hex mesh");
  // check meta data
  rc = fc_getMeshName(mesh, &temp_name);
  fail_unless(rc == FC_SUCCESS, "failed to get mesh name");
  fail_unless(!strcmp(name, temp_name), "mismatch of name");
  free(temp_name);
  rc = fc_getMeshInfo(mesh, &temp_topoDim, &temp_numDim, &temp_numVert, 
		      &temp_numElem, &temp_elemType);
  fail_unless(rc == FC_SUCCESS, "failed to get mesh info");
  fail_unless(temp_topoDim == numDim && temp_numDim == numDim,
	      "mismatch of topodim or dim");
  fail_unless(temp_numVert == numVert, "mismatch of numVert");
  fail_unless(temp_numElem == numElem, "mismatch of numElem");
  fail_unless(temp_elemType == FC_ET_HEX, "mismatch of elemType");
  // check coords - 1st and last, then all
  rc = fc_getMeshCoordsPtr(mesh, &coords);
  fail_unless(rc == FC_SUCCESS, "failed to get coords");
  for (i = 0; i < numDim; i++) { // 1st and last are easy
    fail_unless(coords[i] == lowers[i], "mismatch of 1st vert coords");
    fail_unless(coords[(numVert-1)*numDim+i] == uppers[i],
		"mismatch of last vert coords");
  }
  for (i = 0; i < numVert; i++) {
    getICoords(i, numDim, numVertPerDim, &icoords);
    for (j = 0; j < numDim; j++)
      fail_unless(FC_DBL_EQUIV(coords[i*numDim+j],
			       lowers[j] + icoords[j]*lens[j]),
		  "mismatch of coords");
  }
  // check conns
  rc = fc_getMeshElementConnsPtr(mesh, &conns);
  fail_unless(rc == FC_SUCCESS, "failed to get conns");
  for (i = 0; i < numVertPerElem; i++) { // 1st and last known
    fail_unless(conns[i] == firstElemConns[i], "mismatch of 1st elem conns"); 
    fail_unless(conns[(numElem-1)*numVertPerElem + i] == lastElemConns[i], 
		"mismatch of last elem conns");
  }
  for (i = 0; i < numElem; i++) {
    getICoords(i, numDim, numElemPerDim, &icoords);
    for (j = 0; j < numVertPerElem; j++) {
      vertID = (icoords[2]+vertAdjust[j][2])*numVertPerDim[1]*numVertPerDim[0] +
	       (icoords[1]+vertAdjust[j][1])*numVertPerDim[0] + 
	       (icoords[0]+vertAdjust[j][0]);
      fail_unless(conns[i*numVertPerElem+j] == vertID, "mismatch of conns");
    }
  }

  // *** test doing it wrong

  // bad dataset
  rc = fc_createSimpleHexMesh(bad_dataset, name, numElemPerDim[0], numElemPerDim[1], 
			      numElemPerDim[2], lowers, uppers, &mesh);
  fail_unless(rc != FC_SUCCESS, "should fail with bad dataset");
  fail_unless(FC_HANDLE_EQUIV(mesh, FC_NULL_MESH), "fail should return null");
  // bad name
  rc = fc_createSimpleHexMesh(dataset, NULL, numElemPerDim[0], numElemPerDim[1], 
			      numElemPerDim[2], lowers, uppers, &mesh);
  fail_unless(rc != FC_SUCCESS, "should fail with bad name");
  fail_unless(FC_HANDLE_EQUIV(mesh, FC_NULL_MESH), "fail should return null");
  // bad N, M, or P
  rc = fc_createSimpleHexMesh(dataset, name, 0, numElemPerDim[1], 
			      numElemPerDim[2], lowers, uppers, &mesh);
  fail_unless(rc != FC_SUCCESS, "should fail with N=0");
  fail_unless(FC_HANDLE_EQUIV(mesh, FC_NULL_MESH), "fail should return null");
  rc = fc_createSimpleHexMesh(dataset, name, numElemPerDim[0], 0, 
			      numElemPerDim[2], lowers, uppers, &mesh);
  fail_unless(rc != FC_SUCCESS, "should fail with M=0");
  fail_unless(FC_HANDLE_EQUIV(mesh, FC_NULL_MESH), "fail should return null");
  rc = fc_createSimpleHexMesh(dataset, name, numElemPerDim[0], numElemPerDim[1], 
			      0, lowers, uppers, &mesh);
  fail_unless(rc != FC_SUCCESS, "should fail with P=0");
  fail_unless(FC_HANDLE_EQUIV(mesh, FC_NULL_MESH), "fail should return null");
  rc = fc_createSimpleHexMesh(bad_dataset, name, numElemPerDim[0], numElemPerDim[1], 
			      numElemPerDim[2], lowers, uppers, &mesh);
  fail_unless(rc != FC_SUCCESS, "should fail with bad dataset");
  fail_unless(FC_HANDLE_EQUIV(mesh, FC_NULL_MESH), "fail should return null");
  // bad uppers or lowers
  rc = fc_createSimpleHexMesh(dataset, name, numElemPerDim[0], numElemPerDim[1], 
			      numElemPerDim[2], NULL, uppers, &mesh);
  fail_unless(rc != FC_SUCCESS, "should fail with bad lowers");
  fail_unless(FC_HANDLE_EQUIV(mesh, FC_NULL_MESH), "fail should return null");
  rc = fc_createSimpleHexMesh(dataset, name, numElemPerDim[0], numElemPerDim[1], 
			      numElemPerDim[2], lowers, NULL, &mesh);
  fail_unless(rc != FC_SUCCESS, "should fail with bad uppers");
  fail_unless(FC_HANDLE_EQUIV(mesh, FC_NULL_MESH), "fail should return null");
  // uppers have to be more than lowers
  for (i = 0; i < numDim; i++) {
    for (j = 0; j < numDim; j++)
      temp_coords[j] = lowers[j];
    temp_coords[i] = uppers[i];
    rc = fc_createSimpleHexMesh(dataset, name, numElemPerDim[0], numElemPerDim[1], 
				numElemPerDim[2], temp_coords, uppers, &mesh);
    fail_unless(rc != FC_SUCCESS, "should fail if lowers = uppers");
    fail_unless(FC_HANDLE_EQUIV(mesh, FC_NULL_MESH), "fail should return null");
  }
  // bad mesh pointer
  rc = fc_createSimpleHexMesh(bad_dataset, name, numElemPerDim[0], numElemPerDim[1], 
			      numElemPerDim[2], lowers, uppers, NULL);
  fail_unless(rc != FC_SUCCESS, "should fail with if mesh is NULL");

  // *** all done
  rc = fc_deleteDataset(dataset);
  fail_unless(rc == FC_SUCCESS, "failed to delete dataset at end of testing");
}
END_TEST


// test fc_createSphereMesh()
START_TEST(sphere_mesh)
{
  FC_ReturnCode rc;
  int i, j;
  char name[] = "booya";
  FC_Dataset dataset, bad_dataset = { -99, -99 };
  FC_Coords center = { 1.0, 2.0, 3.0 };
  FC_Mesh m1,m2;
  double *points; //ro
  const double init_radius = 6.0; 
  double radius;

  double rad_vals[] = {4.0, 5.0, 6.0, 7.0};
  double centers[]  = {  1.0,  7.0,  6.0,
			 8.0,  2.0, -3.0,
			-2.0,  8.0,  5.0,
			-5.0, -3.0,  9.0 };
  const int num_step=4;

  int sv_step;
  FC_Variable *sv;
  FC_Sequence seq;
  FC_Variable  var_bad = {-100,-100};
  double *disp_coords;
  double  avg[3];

  int topodim, dim, numVertex, numElement;
  FC_ElementType elemType;



  // *** make a dataset to play in
  rc = fc_createDataset("ball dataset", &dataset);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create dataset");

  // *** do it right
  rc = fc_createSphereMesh(dataset, name, center, init_radius, 102, &m1);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create dataset");
 
  rc = fc_getMeshInfo(m1, &topodim, &dim, &numVertex, &numElement, &elemType);
  fail_unless(rc == FC_SUCCESS, "abort: failed to get mesh info");
  //Note: ball is made of several planes (rings) and a top and bottom node
  //thus,
  // num_planes = sqrt(x-2)
  // nodes = sqrt(x-2)^2 + 2 = num_planes^2+2 
  // triangles = (num_planes-1)*(2*num_planes) + 2*num_planes
  // for 102 nodes-> num_planes = 10
  // exoect num tri to be 200

  fail_unless(topodim==2, "abort: topodim not right");
  fail_unless(dim==3, "abort: dim not right");
  fail_unless(numVertex==102, "abort: topodim not right");
  fail_unless(numElement==200, "abort: topodim not right");
  fail_unless(elemType==FC_ET_TRI, "abort: wrong element type");

  rc = fc_getMeshCoordsPtr(m1, &points);
  fail_unless(rc == FC_SUCCESS, "abort: couldn't get coords");
  for(i=0; i<numVertex; i++){
    //All vertices should be the same as the radius distance away
    rc = fc_calcEuclideanDistance(center, &points[i*3], 3, &radius);
    fail_unless(rc == FC_SUCCESS, "abort: failed to get node distance");
    fail_unless(fc_eqd(radius, init_radius), "abort: radius didn't work out");
  }
  

  //Now try the displacement 
  rc = fc_createSequence(dataset, "booya seq", &seq);
  fail_unless(rc == FC_SUCCESS, "abort: failed to get node distance");
  rc = fc_setSequenceCoords(seq,num_step, FC_DT_DOUBLE, centers);
  fail_unless(rc == FC_SUCCESS, "abort: failed to set sequence");

  rc = fc_createSeqVariable(m1, seq, "Booya SeqVar", &sv_step, &sv);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create seqvar");
  fail_unless(sv_step==num_step, "abort: bad step count");

  //Set the values
  for(i=0; i<sv_step; i++){
    rc = fc_setSphereDisplacementVariable( sv[i], &centers[i*3], rad_vals[i]);
    fail_unless(rc == FC_SUCCESS, "abort: problem setting displacement variable");
  }

  //Check
  for(i=0; i<sv_step; i++){
    rc = fc_getDisplacedMeshCoords(m1, sv[i], &disp_coords); 
    fail_unless(rc == FC_SUCCESS, "abort: problem setting displacement variable");


    //Note: error accumulates here because we average lots of points. Thus
    //      we need to crank epsilon up to help tolerance. 1e-15 still worked
    //      on my machine, but I'm setting the tolerance to 1e-10 for older
    //      machines. In reality, we're dealing with integers, so this should be fine.

    //Find distance from expected center to all values
    for(j=0; j<numVertex; j++){
      rc = fc_calcEuclideanDistance(&centers[i*3], &disp_coords[j*3], 3, &radius);
      //printf("[%d] %.16lf vs %.16lf\n",i, rad_vals[i], radius);
      fail_unless(rc == FC_SUCCESS, "abort: failed to get node distance");
      fail_unless(fc_eq(radius, rad_vals[i], 1e-10, DBL_MIN), "abort: radius didn't work out");
    }
    
    //Averaging all of the points together should give the center. This is similar
    //to the above, BUT, the center should only work out if the points are equally
    //distributed along the sphere (eg, an avg wouldn't give you the center if
    //all the points were on one side).
    memset(avg, 0, 3*sizeof(double));
    for(j=0; j<numVertex*3; j++){
      avg[ j%3 ] += disp_coords[j];
    }

    //printf("\n[%d] ",i);
    for(j=0;j<3;j++){
      avg[j] /= numVertex;
      //printf(" (%.16lf|%.16lf) ",avg[j], centers[i*3+j]);
      fail_unless( fc_eq(avg[j], centers[i*3+j], 1e-10, DBL_MIN), "abort: sphere centers didn't match up");
    }
    //printf("\n");
    free(disp_coords);
  }

  //Should be ok if we give 0 for number of nodes
  rc = fc_createSphereMesh(dataset, "shiznet", center, init_radius, 0, &m2);
  fail_unless(rc == FC_SUCCESS, "abort: didn't use default number of nodes");
  rc = fc_getMeshInfo(m2, NULL, NULL, &numVertex, NULL, NULL);
  fail_unless(rc == FC_SUCCESS, "abort: failed to get mesh info");
  fail_unless(numVertex > 1, "abort: Didn't generate a mesh with points");  



  // *** test doing it wrong
  // bad dataset
  rc = fc_createSphereMesh(bad_dataset, name, center, init_radius, 102, &m2);
  fail_unless(rc != FC_SUCCESS, "should fail on bad dataset");
  fail_unless(FC_HANDLE_EQUIV(m2, FC_NULL_MESH), "fail should return null");

  rc = fc_createSphereMesh(dataset, NULL, center, init_radius, 102, &m2);
  fail_unless(rc != FC_SUCCESS, "should faile on bad name");
  fail_unless(FC_HANDLE_EQUIV(m2, FC_NULL_MESH), "fail should return null");

  rc = fc_createSphereMesh(dataset, name, center, init_radius, 102, NULL);
  fail_unless(rc != FC_SUCCESS, "should fail on bad output variable");

  rc = fc_createSphereMesh(dataset, name, center, -2.0, 102, &m2);
  fail_unless(rc != FC_SUCCESS, "should fail on negative radius");
  fail_unless(FC_HANDLE_EQUIV(m2, FC_NULL_MESH), "fail should return null");


  //Fail the displacement
  rc = fc_setSphereDisplacementVariable( var_bad, &centers[0], rad_vals[0]);
  fail_unless(rc != FC_SUCCESS, "should fail on bad variable");

  free(sv);

  //setup - make a new sv that we can fill in
  rc = fc_createSeqVariable(m1, seq, "Booya SeqVar2", &sv_step, &sv);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create seqvar");
  fail_unless(sv_step==num_step, "abort: bad step count");

  rc = fc_setSphereDisplacementVariable( sv[0], NULL, rad_vals[0]);
  fail_unless(rc != FC_SUCCESS, "should fail on bad variable");
  rc = fc_setSphereDisplacementVariable( sv[0], &centers[0], 0.0);
  fail_unless(rc != FC_SUCCESS, "should fail on bad variable");

  //Try setting a variable twice.. should fail
  rc = fc_setSphereDisplacementVariable( sv[0], &centers[0], rad_vals[0]);
  fail_unless(rc == FC_SUCCESS, "abort: failed to set var");
  rc = fc_setSphereDisplacementVariable( sv[0], &centers[0], rad_vals[0]);
  fail_unless(rc != FC_SUCCESS, "should fail when setting a displacement twice?");
 

  // *** all done
  free(sv);

  rc = fc_deleteDataset(dataset);
  fail_unless(rc == FC_SUCCESS, "failed to delete dataset at end of testing"); 

}
END_TEST

// ******  Globals for the build structs tests

// -- these should be same as in mesh_generator.c
//static double A = 10., B = 8., C = 5.;
static int M = 11, N = 12, P = 6;

// test building edges, faces & neighbors
START_TEST(build_point_mesh)
{
  FC_ReturnCode rc;
  int i;
  int topoDim, numDim, numVertex, numEdge, numFace, numElement;
  int numVertex_good, numEdge_good, numFace_good, numElement_good;
  FC_Dataset dataset;
  FC_Mesh mesh, *returnMeshes;
  int numReturnMeshes;
  FC_ElementType elemType, *faceTypes;
  int* edgeConns, maxNumVertPerFace, *numVertPerFace, *faceConns;
  int *elemToEdgeConns, *elemToFaceConns, *elemFaceOrients; 
  int* numEdgePerVert, **edgeParentsPerVert;
  int* numFacePerVert, **faceParentsPerVert;
  int* numElemPerVert, **elemParentsPerVert;
  int* numElemPerEdge, **elemParentsPerEdge;
  int* numElemPerFace, **elemParentsPerFace;
  int* numVertNeighsViaEdge, **vertNeighsViaEdge;
  int* numElemNeighsViaVert, **elemNeighsViaVert;
  int* numElemNeighsViaEdge, **elemNeighsViaEdge;
  int* numElemNeighsViaFace, **elemNeighsViaFace;

  // pre calc "good" values
  numVertex_good = (M+1)*(N+1);
  numEdge_good = 0;
  numFace_good = 0;
  numElement_good = (M+1)*(N+1);

  // setup
  fc_loadDataset("../data/gen_point.ex2", &dataset);
  rc = fc_getMeshByName(dataset, "point mesh", &numReturnMeshes,&returnMeshes);
  fail_unless(rc == FC_SUCCESS, "Test aborted: mesh not found");
  fail_unless(numReturnMeshes == 1, "failed to find unique mesh by name");
  mesh = returnMeshes[0];
  free(returnMeshes);

  rc = fc_getMeshInfo(mesh, &topoDim, &numDim, &numVertex, &numElement, 
			&elemType);
  fail_unless(rc == FC_SUCCESS, "Test aborted: mesh could not be described");
  fail_unless(topoDim == 0, "Test aborted: expected mesh with topodim of 0");
  fail_unless(numVertex == numVertex_good, "mismatch number of vertices");
  fail_unless(elemType == FC_ET_POINT, "Test aborted: expected point elements");
  fail_unless(rc == FC_SUCCESS, "Test aborted: failed to build mesh strutures");

  // vertices
  fc_getMeshNumVertex(mesh, &numVertex);
  fail_unless(numVertex == numVertex_good, "mismatch number of vertices");
  rc = _fc_getMeshEdgeParentsOfVerticesPtr(mesh, &numEdgePerVert,
				       &edgeParentsPerVert); 
  fail_unless(rc == FC_SUCCESS, "failed to get vertices' edge parents");
  fail_unless(numEdgePerVert == NULL, "mismatch of numEdgePerVert");
  fail_unless(edgeParentsPerVert == NULL, "mismatch of edgeParentsPerVert");
  rc = _fc_getMeshFaceParentsOfVerticesPtr(mesh, &numFacePerVert,
				       &faceParentsPerVert); 
  fail_unless(rc == FC_SUCCESS, "failed to get vertices' face parents");
  fail_unless(numFacePerVert == NULL, "mismatch of numFacePerVert");
  fail_unless(faceParentsPerVert == NULL, "mismatch of faceParentsPerVert");
  rc = _fc_getMeshElementParentsOfVerticesPtr(mesh, &numElemPerVert,
					  &elemParentsPerVert); 
  fail_unless(rc == FC_SUCCESS, "failed to get vertices' elem parents");
  rc = _fc_getMeshVertexNeighborsViaEdgePtr(mesh, &numVertNeighsViaEdge,
					    &vertNeighsViaEdge);
  fail_unless(rc == FC_SUCCESS, "failed to get vert neighs via edge");
  fail_unless(numVertNeighsViaEdge == NULL &&
	      vertNeighsViaEdge == NULL, 
	      "mismatch of vert neighs via edge");
  for (i = 0; i < numVertex; i++) {
    fail_unless(numElemPerVert[i] == 1, "mismatch in numElemPerVert");
    fail_unless(elemParentsPerVert[i][0] == i, 
		"mismatch of elemParentsPerVert");
  }

  // edges
  rc = fc_getMeshNumEdge(mesh, &numEdge);
  fail_unless(numEdge == numEdge_good, "mismatch in number of edges");
  rc = fc_getMeshEdgeConnsPtr(mesh, &edgeConns);
  fail_unless(rc == FC_SUCCESS, "failed to get edge conns");
  fail_unless(edgeConns == NULL, "edge conns should be NULL");
  rc = _fc_getMeshElementParentsOfEdgesPtr(mesh, &numElemPerEdge,
					&elemParentsPerEdge);
  fail_unless(rc == FC_SUCCESS, "failed to get edge's elem parents");
  fail_unless(numElemPerEdge == NULL, "mismatch of numElemPerEdge");
  fail_unless(elemParentsPerEdge == NULL, "mismatch of numElemPerEdge");

  // faces
  rc = fc_getMeshNumFace(mesh, &numFace);
  fail_unless(numFace == numFace_good, "mismatch in number of faces");
  rc = fc_getMeshFaceTypesPtr(mesh, &faceTypes);
  fail_unless(rc == FC_SUCCESS, "failed to get face types");
  fail_unless(faceTypes == NULL, "face types should be NULL");
  rc = fc_getMeshFaceConnsPtr(mesh, &numVertPerFace, &maxNumVertPerFace,
			      &faceConns);
  fail_unless(rc == FC_SUCCESS, "failed to get face conns");
  fail_unless(numVertPerFace == NULL && maxNumVertPerFace == 0 &&
	      faceConns == NULL, "face conns should be NULL");
  rc = _fc_getMeshElementParentsOfFacesPtr(mesh, &numElemPerFace,
					&elemParentsPerFace);
  fail_unless(rc == FC_SUCCESS, "failed to get face's elem parents");
  fail_unless(numElemPerFace == NULL, "mismatch of numElemPerFace");
  fail_unless(elemParentsPerFace == NULL, "mismatch of numElemPerFace");

  // elements
  rc = fc_getMeshNumElement(mesh, &numElement);
  fail_unless(numElement == numElement_good, "mismatch in number of elements");
  rc = fc_getMeshElementToEdgeConnsPtr(mesh, &elemToEdgeConns);
  fail_unless(rc == FC_SUCCESS, "failed to get elem to edge conns");
  fail_unless(elemToEdgeConns == NULL, "elem to edge conns should be NULL");
  rc = fc_getMeshElementToFaceConnsPtr(mesh, &elemToFaceConns);
  fail_unless(rc == FC_SUCCESS, "failed to get elem to edge conns");
  fail_unless(elemToFaceConns == NULL, "elem to edge conns should be NULL");
  rc = fc_getMeshElementFaceOrientsPtr(mesh, &elemFaceOrients);
  fail_unless(rc == FC_SUCCESS, "failed to get elem to edge conns");
  fail_unless(elemFaceOrients == NULL, "elem to edge conns should be NULL");
  rc = _fc_getMeshElementNeighborsViaVertexPtr(mesh, &numElemNeighsViaVert,
					       &elemNeighsViaVert);
  fail_unless(rc == FC_SUCCESS, "failed to get elemNeighsViaVert");
  rc = _fc_getMeshElementNeighborsViaEdgePtr(mesh, &numElemNeighsViaEdge,
					     &elemNeighsViaEdge);
  fail_unless(rc == FC_SUCCESS, "failed to get elemNeighsViaEdge");
  fail_unless(numElemNeighsViaEdge == NULL &&
	      elemNeighsViaEdge == NULL, "mismatch of elemNeighsViaEdge");
  rc = _fc_getMeshElementNeighborsViaFacePtr(mesh, &numElemNeighsViaFace,
					     &elemNeighsViaFace);
  fail_unless(rc == FC_SUCCESS, "failed to get elemNeighsViaFace");
  fail_unless(numElemNeighsViaFace == NULL &&
	      elemNeighsViaFace == NULL, "mismatch of elemNeighsViaEdge");
  for (i = 0; i < numElement; i++) {
    fail_unless(numElemNeighsViaVert[i] == 0, 
		"mismatch of numElemNeighsViaVert");
    fail_unless(elemNeighsViaVert[i] == NULL,
		"mismatch of elemNeighsViaVert");
  }

  fc_deleteDataset(dataset);
}
END_TEST
  
// test building edges, faces & neighbors
START_TEST(build_line_mesh)
{
  FC_ReturnCode rc;
  int i, j, k;
  int topoDim, numDim, numVertex, numEdge, numFace, numElement;
  int numVertex_good, numEdge_good, numFace_good, numElement_good;
  FC_Dataset dataset;
  FC_Mesh mesh, *returnMeshes;
  int numReturnMeshes;
  FC_ElementType elemType, *faceTypes;
  int* edgeConns, maxNumVertPerFace, *numVertPerFace, *faceConns;
  int* elemToEdgeConns, *elemToFaceConns, *elemFaceOrients;
  int* numEdgePerVert, **edgeParentsPerVert;
  int* numFacePerVert, **faceParentsPerVert;
  int* numElemPerVert, **elemParentsPerVert;
  int* numElemPerEdge, **elemParentsPerEdge;
  int* numElemPerFace, **elemParentsPerFace;
  int* numVertNeighsViaEdge, **vertNeighsViaEdge;
  int* numElemNeighsViaVert, **elemNeighsViaVert;
  int* numElemNeighsViaEdge, **elemNeighsViaEdge;
  int* numElemNeighsViaFace, **elemNeighsViaFace;
  int vert_dims[2] = { M + 1, N + 1 };
  int quad_dims[2] = { M, N };


  // NOTE: the line dataset was a quad mesh with every quad element 
  //    replaced by 3 line elements radiating from the "lowest" vertex
  //    In a quad, with counter clockwise numbering of the vertices,
  //    the edge to vertex conns were { 0, 1 } { 0, 2 } { 3, 0 }

  // pre calc "good" values
  numVertex_good = (M+1)*(N+1);
  numEdge_good = 3*M*N;
  numFace_good = 0;
  numElement_good = numEdge_good;

  // setup
  fc_loadDataset("../data/gen_line.ex2", &dataset);
  rc = fc_getMeshByName(dataset, "line mesh",  &numReturnMeshes,&returnMeshes);
  fail_unless(rc == FC_SUCCESS, "Test aborted: mesh not found");
  fail_unless(numReturnMeshes == 1, "failed to find unique mesh by name");
  mesh = returnMeshes[0];
  free(returnMeshes);

  rc = fc_getMeshInfo(mesh, &topoDim, &numDim, &numVertex, &numElement,
			&elemType);
  fail_unless(rc == FC_SUCCESS, "Test aborted: mesh could not be described");
  fail_unless(topoDim == 1, "Text aborted: expected mesh with topodim of 1");
  fail_unless(elemType == FC_ET_LINE, "Test aborted: expected line elements");

  // vertices
  fc_getMeshNumVertex(mesh, &numVertex);
  fail_unless(numVertex == numVertex_good, "mismatch number of vertices");
  rc = _fc_getMeshEdgeParentsOfVerticesPtr(mesh, &numEdgePerVert,
				       &edgeParentsPerVert); 
  fail_unless(rc == FC_SUCCESS, "failed to get vertices' edge parents");
  rc = _fc_getMeshFaceParentsOfVerticesPtr(mesh, &numFacePerVert,
				       &faceParentsPerVert); 
  fail_unless(rc == FC_SUCCESS, "failed to get vertices' face parents");
  fail_unless(numFacePerVert == NULL, "mismatch of numFacePerVert");
  fail_unless(faceParentsPerVert == NULL, "mismatch of faceParentsPerVert");
  rc = _fc_getMeshElementParentsOfVerticesPtr(mesh, &numElemPerVert,
				       &elemParentsPerVert); 
  fail_unless(rc == FC_SUCCESS, "failed to get vertices' elem parents");
  rc = _fc_getMeshVertexNeighborsViaEdgePtr(mesh, &numVertNeighsViaEdge,
					    &vertNeighsViaEdge);
  fail_unless(rc == FC_SUCCESS, "failed to get vert neighs via edge");
  for (i = 0; i < numVertex; i++) {
    // calculate stuff
    // element neighbors will be same as edge neighbors since 1D
    // face neighbors do not exists
    // setup 
    ICoords vert_icoords, temp_icoords;
    int quad_index;
    int maxNumElem = 6;
    int actualNumElem;
    int elemIDs[6] = { -1, -1, -1, -1, -1, -1 };
    int maxNumEdgeOwner = 6;
    int maxNumEdgeNeigh = 6;
    //int maxNumFaceOwner = 0;
    int actualNumEdgeOwner;
    int actualNumEdgeNeigh;
    //int actualNumFaceOwner = maxNumFaceOwner;
    int edgeOwnerIDs[6] = { -1, -1, -1, -1, -1, -1 };
    int edgeNeighIDs[6] = { -1, -1, -1, -1, -1, -1 };
    // get index from current quad
    int edgeOwnerStencil[6][2] = { { -1, -1 }, { 0, -1 }, { -1, 0 },
				   {  0,  0 }, { 0,  0 }, {  0, 0 } };
    //get relative edge
    int edgeLU[6] = { 1, 2, 0, 0, 1, 2 };
    // this stencil adjusts vertex ijk coords
    int edgeNeighStencil[6][2] = { { -1, -1 }, { 0, -1 }, { -1, 0 }, 
				   {  1,  0 }, { 0,  1 }, {  1, 1 } };
    // setup II -shared calcs
    getICoords(i, 2, vert_dims, &vert_icoords);
    // calc edge owners
    for (j = 0; j < maxNumEdgeOwner; j++) {
      for (k = 0; k < 2; k++)
	temp_icoords[k] = vert_icoords[k] + edgeOwnerStencil[j][k];
      edgeOwnerIDs[j] = getIndex(2, quad_dims, temp_icoords);
      if (edgeOwnerIDs[j] > -1)
	edgeOwnerIDs[j] = edgeOwnerIDs[j]*3 + edgeLU[j];
    }
    actualNumEdgeOwner = 0;
    for (j = 0; j < maxNumEdgeOwner; j++)
      if (edgeOwnerIDs[j] > -1)
	actualNumEdgeOwner++;
    // calc edge neighbors
    for (j = 0; j < maxNumEdgeNeigh; j++) {
      for (k = 0; k < 2; k++)
	temp_icoords[k] = vert_icoords[k] + edgeNeighStencil[j][k];
      edgeNeighIDs[j] = getIndex(2, vert_dims, temp_icoords);
    }
    // (connectivity is different on edge of mesh, some edge elems don't exist)
    if (vert_icoords[0] == vert_dims[0] - 1) {
      edgeNeighIDs[1] = -1;
      edgeNeighIDs[4] = -1;
    }
    if (vert_icoords[1] == vert_dims[1] - 1) {
      edgeNeighIDs[2] = -1;
      edgeNeighIDs[3] = -1;
    }
    actualNumEdgeNeigh = 0;
    for (j = 0; j < maxNumEdgeNeigh; j++)
      if (edgeNeighIDs[j] > -1)
	actualNumEdgeNeigh++;
    // calc owning elements
    if (edgeNeighIDs[0] > -1) {
      getICoords(edgeNeighIDs[0], 2, vert_dims, &temp_icoords);
      quad_index = getIndex(2, quad_dims, temp_icoords);
      elemIDs[0] = 3*quad_index + 1;
    }
    if (edgeNeighIDs[1] > -1) {
      getICoords(edgeNeighIDs[1], 2, vert_dims, &temp_icoords);
      quad_index = getIndex(2, quad_dims, temp_icoords);
      elemIDs[1] = 3*quad_index + 2;
    }
    if (edgeNeighIDs[2] > -1) {
      getICoords(edgeNeighIDs[2], 2, vert_dims, &temp_icoords);
      quad_index = getIndex(2, quad_dims, temp_icoords);
      elemIDs[2] = 3*quad_index;
    }
    quad_index = getIndex(2, quad_dims, vert_icoords);
    if (edgeNeighIDs[3] > -1)
      elemIDs[3] = 3*quad_index;
    if (edgeNeighIDs[4] > -1)
      elemIDs[4] = 3*quad_index + 1;
    if (edgeNeighIDs[5] > -1)
      elemIDs[5] = 3*quad_index + 2;
    actualNumElem = 0;
    for (j = 0; j < maxNumElem; j++)
      if (elemIDs[j] > -1)
	actualNumElem++;

    // tests
    fail_unless(numEdgePerVert[i] == actualNumEdgeOwner,
		"mismatch in numEdgePerVert");
    k = 0; 
    for (j = 0; j < numEdgePerVert[i]; j++) {
      while (edgeOwnerIDs[k] < 0)
	k++;
      fail_unless(edgeParentsPerVert[i][j] == edgeOwnerIDs[k],
		  "vertex: mismatch in edgeIDs");
      k++;
    }
    fail_unless(numElemPerVert[i] == actualNumElem,
		"mismatch of numElemPerVert");
    k = 0;
    for (j = 0; j < numElemPerVert[i]; j++) {
      while (elemIDs[k] < 0)
	k++;
      fail_unless(elemParentsPerVert[i][j] == elemIDs[k],
    	"vertex: mismatch in elementIDs");
      k++;
    }
    fail_unless(numVertNeighsViaEdge[i] == actualNumEdgeNeigh, 
		"mismatch numVertNeighsViaEdge");
    k = 0;
    for (j = 0; j < numVertNeighsViaEdge[i]; j++) {
      while (edgeNeighIDs[k] < 0)
	k++;
      fail_unless(vertNeighsViaEdge[i][j] == edgeNeighIDs[k],
		  "mismatch in vertNeighsViaEdge");
      k++;
    }
  }

  // faces
  fc_getMeshNumFace(mesh, &numFace);
  fail_unless(numFace == numFace_good, "mismatch in number faces");
  rc = fc_getMeshFaceTypesPtr(mesh, &faceTypes);
  fail_unless(rc == FC_SUCCESS, "failed to get face types");
  fail_unless(faceTypes == NULL, "face types should be NULL");
  rc = fc_getMeshFaceConnsPtr(mesh, &numVertPerFace, &maxNumVertPerFace,
			      &faceConns);
  fail_unless(rc == FC_SUCCESS, "failed to get face conns");
  fail_unless(numVertPerFace == NULL && maxNumVertPerFace == 0 &&
	      faceConns == NULL, "face conns should be NULL");
  rc = _fc_getMeshElementParentsOfFacesPtr(mesh, &numElemPerFace,
					&elemParentsPerFace);
  fail_unless(rc == FC_SUCCESS, "failed to get face's elem parents");
  fail_unless(numElemPerFace == NULL, "mismatch of numElemPerFace");
  fail_unless(elemParentsPerFace == NULL, "mismatch of numElemPerFace");

  // edges & elements
  fc_getMeshNumEdge(mesh, &numEdge);
  fail_unless(numEdge == numEdge_good, "mismatch in number edges");
  fc_getMeshNumElement(mesh, &numElement);
  fail_unless(numElement == numElement_good, "mismatch in number of elements");
  fail_unless(numElement == numEdge, "abort test: numEdge should eq numElement");
  rc = fc_getMeshEdgeConnsPtr(mesh, &edgeConns);
  fail_unless(rc == FC_SUCCESS, "failed to get edge conns");
  rc = _fc_getMeshElementParentsOfEdgesPtr(mesh, &numElemPerEdge,
					&elemParentsPerEdge);
  fail_unless(rc == FC_SUCCESS, "failed to get edge's elem parents");
  rc = fc_getMeshElementToEdgeConnsPtr(mesh, &elemToEdgeConns);
  fail_unless(rc == FC_SUCCESS, "failed to get elem to edge conns");
  rc = fc_getMeshElementToFaceConnsPtr(mesh, &elemToFaceConns);
  fail_unless(rc == FC_SUCCESS, "failed to get elem to edge conns");
  fail_unless(elemToFaceConns == NULL, "elem to edge conns should be NULL");
  rc = fc_getMeshElementFaceOrientsPtr(mesh, &elemFaceOrients);
  fail_unless(rc == FC_SUCCESS, "failed to get elem to edge conns");
  fail_unless(elemFaceOrients == NULL, "elem to edge conns should be NULL");
  rc = _fc_getMeshElementNeighborsViaVertexPtr(mesh, &numElemNeighsViaVert,
					       &elemNeighsViaVert);
  fail_unless(rc == FC_SUCCESS, "failed to get elemNeighsViaVert");
  rc = _fc_getMeshElementNeighborsViaEdgePtr(mesh, &numElemNeighsViaEdge,
					     &elemNeighsViaEdge);
  fail_unless(rc == FC_SUCCESS, "failed to get elemNeighsViaEdge");
  rc = _fc_getMeshElementNeighborsViaFacePtr(mesh, &numElemNeighsViaFace,
					     &elemNeighsViaFace);
  fail_unless(rc == FC_SUCCESS, "failed to get elemNeighsViaFace");
  fail_unless(numElemNeighsViaFace == NULL &&
	      elemNeighsViaFace == NULL, "mismatch of elemNeighsViaEdge");
  for (i = 0; i < numElement; i++) {
    // calculate stuff
    // setup
    int quad_index, temp_quad_index;
    int relative_edgeID;
    ICoords quad_icoords, temp_icoords;
    int vertIDs[2];
    //int maxNumVertexNeigh = 10;
    int actualNumVertexNeigh;
    int vertexNeighIDs[10] =  { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 }; 
    // this stencil adjusts ijk coords for quads
    int quadNeighStencil[9][2] = { { -1, -1 }, { 0, -1 }, {  1, -1 },  
				   { -1,  0 }, { 0,  0 }, {  1,  0 },
				   { -1,  1 }, { 0,  1 }, {  1,  1 } };
    // relative edge ID in the neighbor quad
    // first index is the relative edge ID in current quad, second index
    // is the quad ID as found with quad NeighStencil above, third index
    // is the relative edgeID in the neighbor quad. 1 means edge is 
    // a neighbor of the originating edge.
    int edgeLU[3][9][3] = {{ { 0, 1, 0 }, { 0, 1, 1 }, { 0, 0, 1 },
			     { 1, 0, 0 }, { 0, 1, 1 }, { 1, 1, 1 },
			     { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 } }, 
			   { { 0, 1, 0 }, { 0, 0, 1 }, { 0, 0, 0 },
			     { 1, 0, 0 }, { 1, 0, 1 }, { 0, 0, 1 },
			     { 0, 0, 0 }, { 1, 0, 0 }, { 1, 1, 1 } }, 
			   { { 0, 1, 0 }, { 0, 0, 1 }, { 0, 0, 0 },
			     { 1, 1, 0 }, { 1, 1, 0 }, { 0, 0, 0 },
			     { 1, 0, 0 }, { 1, 1, 1 }, { 0, 0, 0 } } };
    // setup II - shared calcs      
    quad_index = i/3;
    relative_edgeID = i%3;
    getICoords(quad_index, 2, quad_dims, &quad_icoords);
    // owned vertices
    temp_icoords[0] = quad_icoords[0];
    temp_icoords[1] = quad_icoords[1];
    vertIDs[0] = getIndex(2, vert_dims, temp_icoords);
    switch (relative_edgeID) {
    case 0:  temp_icoords[0]++;                      break;
    case 1:  temp_icoords[0]++;  temp_icoords[1]++;  break;
    case 2:                      temp_icoords[1]++;  break;
    }
    vertIDs[1] = getIndex(2, vert_dims, temp_icoords);
    // vertex neighbors
    actualNumVertexNeigh = 0;
    for (j = 0; j < 9; j++) {
      for (k = 0; k < 2; k++)
	temp_icoords[k] = quad_icoords[k] + quadNeighStencil[j][k];
      temp_quad_index = getIndex(2, quad_dims, temp_icoords);
      if (temp_quad_index > -1) {
	for (k = 0; k < 3; k++) {
	  if (edgeLU[relative_edgeID][j][k] > 0) {
	    vertexNeighIDs[actualNumVertexNeigh] = temp_quad_index*3 + k;
	    actualNumVertexNeigh++;
	  }
	}
      }
    }
    // edge and face neighbors are not possible on edge mesh

    // edge tests
    fail_unless(numElemPerEdge[i] == 1, "mismatch of numElemPerEdge");
    fail_unless(elemParentsPerEdge[i][0] == i,
		"mismatch of elemParentsPerEdge");
    fail_unless(edgeConns[i*2] == vertIDs[0], "mismatch of edgeConns");
    fail_unless(edgeConns[i*2+1] == vertIDs[1], "mismatch of edgeConns");

    // element tests
    fail_unless(elemToEdgeConns[i] == i, "mismatch of elemToEdgeConns");
    fail_unless(numElemNeighsViaVert[i] == actualNumVertexNeigh,
		"mismatch of numElemNeighsViaVert");
    fail_unless(!memcmp(elemNeighsViaVert[i], vertexNeighIDs,
		actualNumVertexNeigh*sizeof(int)), 
		"mismatch in elemNeighsViaVert");
    fail_unless(numElemNeighsViaEdge[i] == 0,
		"mismatch of numElemNeighsViaEdge");
    fail_unless(elemNeighsViaEdge[i] == NULL,
		"mismatch of elemNeighsViaEdge");
  }

  fc_deleteDataset(dataset);
}
END_TEST
  
// test building edges, faces & neighbors
START_TEST(build_quad_mesh)
{
  FC_ReturnCode rc;
  int i, j, k;
  int topoDim, numDim, numVertex, numEdge, numFace, numElement;
  int numVertex_good, numEdge_good, numFace_good, numElement_good;
  FC_Dataset dataset;
  FC_Mesh mesh, *returnMeshes;
  int numReturnMeshes;
  FC_ElementType elemType, *faceTypes;
  int* edgeConns, *faceConns, maxNumVertPerFace, *numVertPerFace;
  int* elemToEdgeConns, *elemToFaceConns, *elemFaceOrients;
  int* numEdgePerVert, **edgeParentsPerVert;
  int* numFacePerVert, **faceParentsPerVert;
  int* numElemPerVert, **elemParentsPerVert;
  int* numElemPerEdge, **elemParentsPerEdge;
  int* numElemPerFace, **elemParentsPerFace;
  int* numVertNeighsViaEdge, **vertNeighsViaEdge;
  int* numElemNeighsViaVert, **elemNeighsViaVert;
  int* numElemNeighsViaEdge, **elemNeighsViaEdge;
  int* numElemNeighsViaFace, **elemNeighsViaFace;
  int vert_dims[2] = { M + 1, N + 1 };
  int quad_dims[2] = { M, N };
  int edge_vertIDs[2*M*N + M + N][2]; // number of edges
  int edge_numElem[2*M*N + M + N];    // number of edges
  int edge_elemIDs[2*M*N + M + N][2]; // number of edges
  int elem_edgeIDs[M*N][4]; // number of elements


  // pre calc "good" values
  numVertex_good = (M+1)*(N+1);
  numEdge_good = 2*M*N + M + N;
  numFace_good = M*N;
  numElement_good = numFace_good;

  // setup
  fc_loadDataset("../data/gen_quad.ex2", &dataset);
  rc = fc_getMeshByName(dataset, "quad mesh",  &numReturnMeshes,&returnMeshes);
  fail_unless(rc == FC_SUCCESS, "Test aborted: mesh not found\n");
  fail_unless(numReturnMeshes == 1, "failed to find unique mesh by name");
  mesh = returnMeshes[0];
  free(returnMeshes);

  rc = fc_getMeshInfo(mesh, &topoDim, &numDim, &numVertex, &numElement,
			&elemType);
  fail_unless(rc == FC_SUCCESS, "Test aborted: mesh could not be described\n");
  fail_unless(topoDim == 2, "Text aborted: expected mesh with topodim of 2");
  fail_unless(elemType == FC_ET_QUAD, "Test aborted: expected quad elements");

  // build edges info - vertex members, and back and forth 
  // conns to elems
  numEdge = 0; 
  for (i = 0; i < numElement; i++) {
    // setup
    int temp_ID;
    ICoords quad_icoords, temp_icoords;
    getICoords(i, topoDim, quad_dims, &quad_icoords);

    // Edge#0 - sometimes
    if (quad_icoords[1]%quad_dims[1] == 0) {
      edge_elemIDs[numEdge][0] = i;
      edge_numElem[numEdge] = 1;
      elem_edgeIDs[i][0] = numEdge;
      ICoordscpy(quad_icoords, &temp_icoords);
      edge_vertIDs[numEdge][0] = getIndex(topoDim, vert_dims, temp_icoords);
      temp_icoords[0]++;
      edge_vertIDs[numEdge][1] = getIndex(topoDim, vert_dims, temp_icoords);
      numEdge++;
    }
    // Edge#1 (same as edge#3 in quad neighbor x++)
    {
      edge_elemIDs[numEdge][0] = i;
      edge_numElem[numEdge] = 1;
      elem_edgeIDs[i][1] = numEdge;
      ICoordscpy(quad_icoords, &temp_icoords);
      temp_icoords[0]++;
      temp_ID = getIndex(topoDim, quad_dims, temp_icoords);
      if (temp_ID > -1) {
	edge_elemIDs[numEdge][1] = temp_ID;
	edge_numElem[numEdge]++;
	elem_edgeIDs[temp_ID][3] = numEdge;
      }
      ICoordscpy(quad_icoords, &temp_icoords);
      temp_icoords[0]++;
      edge_vertIDs[numEdge][0] = getIndex(topoDim, vert_dims, temp_icoords);
      temp_icoords[1]++;
      edge_vertIDs[numEdge][1] = getIndex(topoDim, vert_dims, temp_icoords);
      numEdge++;
    }
    // Edge#2 (same as edge#0 in quad neighbor y++)
    {
      edge_elemIDs[numEdge][0] = i;
      edge_numElem[numEdge] = 1;
      elem_edgeIDs[i][2] = numEdge;
      ICoordscpy(quad_icoords, &temp_icoords);
      temp_icoords[1]++;
      temp_ID = getIndex(topoDim, quad_dims, temp_icoords);
      if (temp_ID > -1) {
	edge_elemIDs[numEdge][1] = temp_ID;
	edge_numElem[numEdge]++;
	elem_edgeIDs[temp_ID][0] = numEdge;
      }
      ICoordscpy(quad_icoords, &temp_icoords);
      temp_icoords[1]++;
      edge_vertIDs[numEdge][0] = getIndex(topoDim, vert_dims, temp_icoords);
      temp_icoords[0]++;
      edge_vertIDs[numEdge][1] = getIndex(topoDim, vert_dims, temp_icoords);
      numEdge++;
    }
    // Edge#3
    if (quad_icoords[0]%quad_dims[0] == 0) {
      edge_elemIDs[numEdge][0] = i;
      edge_numElem[numEdge] = 1;
      elem_edgeIDs[i][3] = numEdge;
      ICoordscpy(quad_icoords, &temp_icoords);
      edge_vertIDs[numEdge][0] = getIndex(topoDim, vert_dims, temp_icoords);
      temp_icoords[1]++;
      edge_vertIDs[numEdge][1] = getIndex(topoDim, vert_dims, temp_icoords);
      numEdge++;	
    }
  }
  fail_unless(numEdge == numEdge_good, 
	      "failed to build edges for testing, should never happen");

  // vertices
  fc_getMeshNumVertex(mesh, &numVertex);
  fail_unless(numVertex == numVertex_good, "mismatch number of vertices");
  rc = _fc_getMeshEdgeParentsOfVerticesPtr(mesh, &numEdgePerVert,
				       &edgeParentsPerVert); 
  fail_unless(rc == FC_SUCCESS, "failed to get vertices' edge parents");
  rc = _fc_getMeshFaceParentsOfVerticesPtr(mesh, &numFacePerVert,
				       &faceParentsPerVert); 
  fail_unless(rc == FC_SUCCESS, "failed to get vertices' face parents");
  rc = _fc_getMeshElementParentsOfVerticesPtr(mesh, &numElemPerVert,
				       &elemParentsPerVert); 
  fail_unless(rc == FC_SUCCESS, "failed to get vertices' elem parents");
  rc = _fc_getMeshVertexNeighborsViaEdgePtr(mesh, &numVertNeighsViaEdge,
					    &vertNeighsViaEdge);
  fail_unless(rc == FC_SUCCESS, "failed to get vert neighs via edge");
  for (i = 0; i < numVertex; i++) {
    // calculate stuff
    // element neighbors will be same as face neighbors since 2D
    // setup 
    ICoords vert_icoords, temp_icoords;
    int temp_ID;
    int maxNumElem = 4;
    int actualNumElem;
    int actualNumEdgeOwner;
    int* elemIDs;   // faceOwner is same as elements
    int* edgeOwnerIDs;
    FC_SortedIntArray elemOwnerList;
    FC_SortedIntArray edgeOwnerList;
    int maxNumEdgeNeigh = 4;
    int actualNumEdgeNeigh;
    int edgeNeighIDs[4];
    // this stencil adjusts owning quad ijk coords (from the vertex coords)
    int elemOwnerStencil[4][2] = {  { -1, -1 }, {  0, -1 }, 
				    { -1,  0 }, {  0,  0 } };
    // lookup owning edges in the quads
    int edgeLU[4][2] = { { 1, 2 }, { 2, 3 }, { 0, 1 }, { 0, 3 } };
    // this stencil adjusts vertex ijk coords
    int edgeNeighStencil[4][2] = { {  0, -1 }, { -1,  0 }, 
				   {  1,  0 }, {  0,  1 } };
    // setup II - shared calcs
    getICoords(i, topoDim, vert_dims, &vert_icoords);
    // calc owners
    fc_initSortedIntArray(&elemOwnerList);
    fc_initSortedIntArray(&edgeOwnerList);
    for (j = 0; j < maxNumElem; j++) {
      for (k = 0; k < 2; k++)
	temp_icoords[k] = vert_icoords[k] + elemOwnerStencil[j][k];
      temp_ID = getIndex(topoDim, quad_dims, temp_icoords);
      if (temp_ID > -1) {
	fc_addIntToSortedIntArray(&elemOwnerList, temp_ID);
	for (k = 0; k < 2; k++)
	  fc_addIntToSortedIntArray(&edgeOwnerList,
				    elem_edgeIDs[temp_ID][edgeLU[j][k]]);
      }
    }
    fc_convertSortedIntArrayToIntArray(&elemOwnerList, &actualNumElem, 
				       &elemIDs);
    fc_convertSortedIntArrayToIntArray(&edgeOwnerList, &actualNumEdgeOwner, 
				    &edgeOwnerIDs);
    // calc edge neighbors
    actualNumEdgeNeigh = 0;
    for (j = 0; j < maxNumEdgeNeigh; j++) {
      for (k = 0; k < 2; k++)
	temp_icoords[k] = vert_icoords[k] + edgeNeighStencil[j][k];
      temp_ID = getIndex(topoDim, vert_dims, temp_icoords);
      if (temp_ID > -1) {
	edgeNeighIDs[actualNumEdgeNeigh] = temp_ID;
	actualNumEdgeNeigh++;
      }
    }

    // tests
    fail_unless(numEdgePerVert[i] == actualNumEdgeOwner,
		"mismatch in numEdgePerVert");
    fail_unless(!memcmp(edgeParentsPerVert[i], edgeOwnerIDs,
			numEdgePerVert[i]*sizeof(int)), 
		"mismatch in edgeParentsPerVert");
    free(edgeOwnerIDs);
    fail_unless(numFacePerVert[i] == actualNumElem, 
		"mismatch in numFacePerVert");
    fail_unless(numElemPerVert[i] == actualNumElem, 
		"vertex: mismatch in number of elements");
    fail_unless(!memcmp(faceParentsPerVert[i], elemIDs,
			numElemPerVert[i]*sizeof(int)),
		"vertex: mismatch in elementIDs");
    fail_unless(!memcmp(elemParentsPerVert[i], elemIDs,
			numElemPerVert[i]*sizeof(int)), 
		"vertex: mismatch in elementIDs");
    free(elemIDs);
    fail_unless(numVertNeighsViaEdge[i] == actualNumEdgeNeigh,
		"mismatch of numVertNeighsViaEdge");
    fail_unless(!memcmp(vertNeighsViaEdge[i], edgeNeighIDs,
			numVertNeighsViaEdge[i]*sizeof(int)), 
		"mismatch of vertNeighsViaEdge");
  }
  
  // edges
  fc_getMeshNumEdge(mesh, &numEdge);
  fail_unless(rc == FC_SUCCESS, "failed to get edges");
  rc = fc_getMeshEdgeConnsPtr(mesh, &edgeConns);
  fail_unless(rc == FC_SUCCESS, "failed to get edge conns");
  rc = _fc_getMeshElementParentsOfEdgesPtr(mesh, &numElemPerEdge,
					&elemParentsPerEdge);
  fail_unless(rc == FC_SUCCESS, "failed to get edge's elem parents");
  for (i = 0; i < numEdge; i++) {
    // tests
    fail_unless(!memcmp(&edgeConns[i*2], edge_vertIDs[i], 2*sizeof(int)), 
		  "edge: mismatch of vertIDs");
    fail_unless(numElemPerEdge[i] == edge_numElem[i], 
		"mismatch in numElemPerEdge");
    fail_unless(!memcmp(elemParentsPerEdge[i], edge_elemIDs[i],
			numElemPerEdge[i]*sizeof(int)), 
		"edge: mismatch of elemIDs");
  }

  // faces & elements
  fc_getMeshNumFace(mesh, &numFace);
  fail_unless(numFace == numFace_good, "mismatch in number faces");
  rc = fc_getMeshFaceTypesPtr(mesh, &faceTypes);
  fail_unless(rc == FC_SUCCESS, "failed to get face types");
  rc = fc_getMeshFaceConnsPtr(mesh, &numVertPerFace, &maxNumVertPerFace,
			      &faceConns); 
  fail_unless(rc == FC_SUCCESS, "failed to get face conns");
  rc = _fc_getMeshElementParentsOfFacesPtr(mesh, &numElemPerFace,
					&elemParentsPerFace);
  fail_unless(rc == FC_SUCCESS, "failed to get face's elem parents");
  fc_getMeshNumElement(mesh, &numElement);
  fail_unless(numElement == numElement_good, "mismatch in number of elements");
  fail_unless(numElement == numFace, 
	      "abort test: numFace should equal numElement");
  rc = fc_getMeshElementToEdgeConnsPtr(mesh, &elemToEdgeConns);
  fail_unless(rc == FC_SUCCESS, "failed to get elem to edge conns");
  rc = fc_getMeshElementToFaceConnsPtr(mesh, &elemToFaceConns);
  fail_unless(rc == FC_SUCCESS, "failed to get elem to edge conns");
  rc = fc_getMeshElementFaceOrientsPtr(mesh, &elemFaceOrients);
  fail_unless(rc == FC_SUCCESS, "failed to get elem to edge conns");
  rc = _fc_getMeshElementNeighborsViaVertexPtr(mesh, &numElemNeighsViaVert,
					       &elemNeighsViaVert);
  fail_unless(rc == FC_SUCCESS, "failed to get elemNeighsViaVert");
  rc = _fc_getMeshElementNeighborsViaEdgePtr(mesh, &numElemNeighsViaEdge,
					     &elemNeighsViaEdge);
  fail_unless(rc == FC_SUCCESS, "failed to get elemNeighsViaEdge");
  rc = _fc_getMeshElementNeighborsViaFacePtr(mesh, &numElemNeighsViaFace,
					     &elemNeighsViaFace);
  fail_unless(rc == FC_SUCCESS, "failed to get elemNeighsViaFace");
  for (i = 0; i < numElement; i++) {
    // calculate stuff
    // setup
    int quad_index;
    ICoords quad_icoords, temp_icoords;
    int numVert = 4;
    int vertIDs[4];
    int maxNumVertNeigh = 8; // 4 + 4
    int maxNumEdgeNeigh = 4;
    int actualNumVertNeigh;
    int actualNumEdgeNeigh;
    int vertNeighIDs[8]; 
    int edgeNeighIDs[4];
    // this stencil adjusts ijk coords for quads
    int vertNeighStencil[8][2] = { { -1, -1 }, {  0, -1 }, {  1, -1 },  
				   { -1,  0 }, {  1,  0 },
				   { -1,  1 }, {  0,  1 }, {  1,  1 } };
    int edgeNeighStencil[4][2] = { {  0, -1 }, { -1,  0 }, 
				   {  1,  0 }, {  0,  1 } };
    // setup II - shared calcs      
    getICoords(i, topoDim, quad_dims, &quad_icoords);
    // owned vertices
    for (j = 0; j < numVert; j++) { // walk around quad
      switch (j) {
      case 0:  
	temp_icoords[0] = quad_icoords[0]; 
	temp_icoords[1] = quad_icoords[1]; 
	break;
      case 1:  temp_icoords[0]++;  break;
      case 2:  temp_icoords[1]++;  break;
      case 3:  temp_icoords[0]--;  break;
      }
      vertIDs[j] = getIndex(topoDim, vert_dims, temp_icoords);
    }
    // vertex neighbors
    actualNumVertNeigh = 0;
    for (j = 0; j < maxNumVertNeigh; j++) {
      for (k = 0; k < 2; k++)
	temp_icoords[k] = quad_icoords[k] + vertNeighStencil[j][k];
      quad_index = getIndex(topoDim, quad_dims, temp_icoords);
      if (quad_index > -1) {
	vertNeighIDs[actualNumVertNeigh] = quad_index;
	actualNumVertNeigh++;
      }
    }
    // edge neighbors
    actualNumEdgeNeigh = 0;
    for (j = 0; j < maxNumEdgeNeigh; j++) {
      for (k = 0; k < 2; k++)
	temp_icoords[k] = quad_icoords[k] + edgeNeighStencil[j][k];
      quad_index = getIndex(topoDim, quad_dims, temp_icoords);
      if (quad_index > -1) {
	edgeNeighIDs[actualNumEdgeNeigh] = quad_index;
	actualNumEdgeNeigh++;
      }
    }
    // face neighbors are not possible on a 2D mesh

    // face tests
    fail_unless(numElemPerFace[i] == 1, "mismatch in numElemPerFace");
    fail_unless(elemParentsPerFace[i][0] == i, 
		"mismatch in elemParentsPerFace");

    // joint tests
    fail_unless(faceTypes[i] == FC_ET_QUAD, "face: mismatch in face type");
    fail_unless(numVertPerFace[i] = numVert, 
		"face: mismatch in number of vertices");
    for (j = 0; j < numVert; j++) {
      fail_unless(faceConns[i*maxNumVertPerFace+j] == vertIDs[j],
		  "face: mismatch in vertex IDs");
    }
    fail_unless(elemToFaceConns[i] == i, "mismatch of elem to face conns");
    fail_unless(elemFaceOrients[i] == 1, "mismatch of face orients");

    // element tests
    fail_unless(!memcmp(&elemToEdgeConns[i*4], elem_edgeIDs[i], 4*sizeof(int)),
		"mismatch of elem to edge conns");
    fail_unless(numElemNeighsViaVert[i] == actualNumVertNeigh,
		"mismatch of numElemNeighsViaVert");
    fail_unless(!memcmp(elemNeighsViaVert[i], vertNeighIDs,
			actualNumVertNeigh*sizeof(int)), 
		"element: mismatch in vertex neighbors");
    fail_unless(numElemNeighsViaEdge[i] == actualNumEdgeNeigh,
		"mismatch of numElemNeighsViaEdge");
    fail_unless(!memcmp(elemNeighsViaEdge[i], edgeNeighIDs,
			actualNumEdgeNeigh*sizeof(int)), 
		"mismatch of elemNeighsViaEdge");
    fail_unless(numElemNeighsViaFace[i] == 0, 
		"mismatch of numElemNeighsViaFace");
    fail_unless(elemNeighsViaFace[i] == NULL, 
		"mismatch of elemNeighsViaFace");
  }
  fc_deleteDataset(dataset);
}
END_TEST
  
// test building edges, faces & neighbors
START_TEST(build_hex_mesh)
{
  FC_ReturnCode rc;
  int i, j, k;
  int topoDim,numDim, numVertex, numEdge, numFace, numElement;
  int numVertex_good, numEdge_good, numFace_good, numElement_good;
  int numEdgePerHex = 12;
  int numFacePerHex = 6;
  FC_Dataset dataset;
  FC_Mesh mesh, *returnMeshes;
  int numReturnMeshes;
  FC_ElementType elemType, *faceTypes;
  int* edgeConns, *faceConns, maxNumVertPerFace, *numVertPerFace;
  int* elemToEdgeConns, *elemToFaceConns, *elemFaceOrients;
  int* numEdgePerVert, **edgeParentsPerVert;
  int* numFacePerVert, **faceParentsPerVert;
  int* numElemPerVert, **elemParentsPerVert;
  int* numElemPerEdge, **elemParentsPerEdge;
  int* numElemPerFace, **elemParentsPerFace;
  int* numVertNeighsViaEdge, **vertNeighsViaEdge;
  int* numElemNeighsViaVert, **elemNeighsViaVert;
  int* numElemNeighsViaEdge, **elemNeighsViaEdge;
  int* numElemNeighsViaFace, **elemNeighsViaFace;
  int vert_dims[3] = { M + 1, N + 1, P + 1 };
  int hex_dims[3] = { M, N, P };
  int edge_vertIDs[3*M*N*P + 2*M*N + 2*N*P + 2*M*P + M + N + P][2]; 
  int edge_numElem[3*M*N*P + 2*M*N + 2*N*P + 2*M*P + M + N + P];  
  int edge_elemIDs[3*M*N*P + 2*M*N + 2*N*P + 2*M*P + M + N + P][4];
  int face_vertIDs[3*M*N*P + M*N + N*P + M*P][4];
  int face_numElem[3*M*N*P + M*N + N*P + M*P];
  int face_elemIDs[3*M*N*P + M*N + N*P + M*P][2];
  int elem_edgeIDs[M*N*P][12]; // number of elements
  int elem_faceIDs[M*N*P][6];
  int elem_faceOrients[M*N*P][6];

  // pre calc "good" values
  numVertex_good = (M+1)*(N+1)*(P+1);
  numEdge_good = 3*M*N*P + 2*M*N + 2*N*P + 2*M*P + M + N + P;
  numFace_good = 3*M*N*P + M*N + N*P + M*P;
  numElement_good = M*N*P;

  // setup
  fc_loadDataset("../data/gen_hex.ex2", &dataset);
  rc = fc_getMeshByName(dataset, "hex mesh", &numReturnMeshes,&returnMeshes);
  fail_unless(rc == FC_SUCCESS, "Test aborted: mesh not found");
  fail_unless(numReturnMeshes == 1, "failed to find unique mesh by name");
  mesh = returnMeshes[0];
  free(returnMeshes);

  rc = fc_getMeshInfo(mesh, &topoDim, &numDim, &numVertex, &numElement,
			&elemType);
  fail_unless(rc == FC_SUCCESS, "Test aborted: mesh could not be described\n");
  fail_unless(topoDim == 3, "Text aborted: expected mesh with topodim of 3");
  fail_unless(elemType == FC_ET_HEX, "Test aborted: expected hex elements");

  // build edges & faces info - vertex members, and back and forth
  // connections to elems
  //
  // Basic procedure: walk every element, then every face of that element,
  // then every edge, making new faces and edges as we go.  Because edges and
  // faces are shared, we don't make every face and edge of every element. The
  // idea is that "downstream" faces and edges of a hex (like the top of the
  // hex) are always made and then pushed onto downstream hexes. The only time
  // "upstream" edges  and faces are made is if they are on the boundary of
  // the mesh. Procedurally, roughly
  //    for each element { 
  //      do the faces {
  //        skip upstream faces if element not on the boundary
  //        make face on this element
  //        add face to downstream elements
  //        do the edges {
  //          skip upstream edges if element not on the boundary
  //          make edge
  //          add edge to downstream elements
  //        }
  //      }
  //    }
  // It's probably helpful to examine mesh.c's fc_buildFaces and fc_buildEdges
  // to get an idea of the order used.
  numEdge = 0; 
  numFace = 0;
  for (i = 0; i < numEdge_good; i++) {
    edge_numElem[i] = 0;
  }
  for (i = 0; i < numElement; i++) {
    // setup
    int temp_ID;
    ICoords hex_icoords, temp_icoords;
    getICoords(i, topoDim, hex_dims, &hex_icoords);

    // Face#0 - sometimes
    if (hex_icoords[1]%hex_dims[1] == 0) {
      face_elemIDs[numFace][0] = i;
      face_numElem[numFace] = 1;
      elem_faceIDs[i][0] = numFace;
      elem_faceOrients[i][0] = 1;
      // Vertices
      ICoordscpy(hex_icoords, &temp_icoords);
      temp_ID = getIndex(topoDim, vert_dims, temp_icoords);
      face_vertIDs[numFace][0] = temp_ID;
      temp_icoords[0]++;
      temp_ID = getIndex(topoDim, vert_dims, temp_icoords);
      face_vertIDs[numFace][1] = temp_ID;
      temp_icoords[2]++;
      temp_ID = getIndex(topoDim, vert_dims, temp_icoords);
      face_vertIDs[numFace][2] = temp_ID;
      temp_icoords[0]--;
      temp_ID = getIndex(topoDim, vert_dims, temp_icoords);
      face_vertIDs[numFace][3] = temp_ID;
      // edges - Do 0-3
      // Edge#0 (#0 local to face)
      if (hex_icoords[1]%hex_dims[1] == 0 && hex_icoords[2]%hex_dims[2] == 0) {
	edge_elemIDs[numEdge][edge_numElem[numEdge]] = i;
	edge_numElem[numEdge]++;
	elem_edgeIDs[i][0] = numEdge;
	edge_vertIDs[numEdge][0] = face_vertIDs[numFace][0]; 
	edge_vertIDs[numEdge][1] = face_vertIDs[numFace][1];
	numEdge++;
      }
      // Edge#1 (#1 local to face) (same as Edge#3 in hex neighbor x++) 
      if (hex_icoords[1]%hex_dims[1] == 0) {
	edge_elemIDs[numEdge][edge_numElem[numEdge]] = i;
	edge_numElem[numEdge]++;
	elem_edgeIDs[i][1] = numEdge;
	ICoordscpy(hex_icoords, &temp_icoords);
	temp_icoords[0]++;
	temp_ID = getIndex(topoDim, hex_dims, temp_icoords);
	if (temp_ID > -1) {
	  edge_elemIDs[numEdge][edge_numElem[numEdge]] = temp_ID;
	  edge_numElem[numEdge]++;
	  elem_edgeIDs[temp_ID][3] = numEdge;
	}
	edge_vertIDs[numEdge][0] = face_vertIDs[numFace][1];
	edge_vertIDs[numEdge][1] = face_vertIDs[numFace][2];
	numEdge++;
      }
      // Edge#2 (#2 local to face) (same as Edge#0 in hex neighbor z++)
      if (hex_icoords[1]%hex_dims[1] == 0) {
	edge_elemIDs[numEdge][edge_numElem[numEdge]] = i;
	edge_numElem[numEdge]++;
	elem_edgeIDs[i][2] = numEdge;
	ICoordscpy(hex_icoords, &temp_icoords);
	temp_icoords[2]++;
	temp_ID = getIndex(topoDim, hex_dims, temp_icoords);
	if (temp_ID > -1) {
	  edge_elemIDs[numEdge][edge_numElem[numEdge]] = temp_ID;
	  edge_numElem[numEdge]++;
	  elem_edgeIDs[temp_ID][0] = numEdge;
	}
	edge_vertIDs[numEdge][0] = face_vertIDs[numFace][2];
	edge_vertIDs[numEdge][1] = face_vertIDs[numFace][3];
	numEdge++;
      }
      // Edge#3 (#3 local to face)
      if (hex_icoords[0]%hex_dims[0] == 0 && hex_icoords[1]%hex_dims[1] == 0) {
	edge_elemIDs[numEdge][edge_numElem[numEdge]] = i;
	edge_numElem[numEdge]++;
	elem_edgeIDs[i][3] = numEdge;	
	edge_vertIDs[numEdge][0] = face_vertIDs[numFace][3];
	edge_vertIDs[numEdge][1] = face_vertIDs[numFace][0];
	numEdge++;
      }
      numFace++;
    }
    // Face#1 (same as Face#3 in Hex neighbor with coords x++)
    { 
      face_elemIDs[numFace][0] = i;
      face_numElem[numFace] = 1;
      elem_faceIDs[i][1] = numFace;
      elem_faceOrients[i][1] = 1;
      ICoordscpy(hex_icoords, &temp_icoords);
      temp_icoords[0]++;
      temp_ID = getIndex(topoDim, hex_dims, temp_icoords);
      if (temp_ID > -1) {
	face_elemIDs[numFace][1] = temp_ID;
	face_numElem[numFace]++;
	elem_faceIDs[temp_ID][3] = numFace;
	elem_faceOrients[temp_ID][3] = -1;
      }
      // Vertices
      ICoordscpy(hex_icoords, &temp_icoords);
      temp_icoords[0]++;
      temp_ID = getIndex(topoDim, vert_dims, temp_icoords);
      face_vertIDs[numFace][0] = temp_ID;
      temp_icoords[1]++;
      temp_ID = getIndex(topoDim, vert_dims, temp_icoords);
      face_vertIDs[numFace][1] = temp_ID;
      temp_icoords[2]++;
      temp_ID = getIndex(topoDim, vert_dims, temp_icoords);
      face_vertIDs[numFace][2] = temp_ID;
      temp_icoords[1]--;
      temp_ID = getIndex(topoDim, vert_dims, temp_icoords);
      face_vertIDs[numFace][3] = temp_ID;
      // edges - Do 4-6
      // Edge#4 (#0 local to face) (same as edge#10 in hex neighbor x++) 
      if (hex_icoords[2]%hex_dims[2] == 0) {
	edge_elemIDs[numEdge][edge_numElem[numEdge]] = i;
	edge_numElem[numEdge]++;
	elem_edgeIDs[i][4] = numEdge;
	ICoordscpy(hex_icoords, &temp_icoords);
	temp_icoords[0]++;
	temp_ID = getIndex(topoDim, hex_dims, temp_icoords);
	if (temp_ID > -1) {
	  edge_elemIDs[numEdge][edge_numElem[numEdge]] = temp_ID;
	  edge_numElem[numEdge]++;
	  elem_edgeIDs[temp_ID][10] = numEdge;
	}
	edge_vertIDs[numEdge][0] = face_vertIDs[numFace][0];
	edge_vertIDs[numEdge][1] = face_vertIDs[numFace][1];
	numEdge++;
      }
      // Edge#5 (#1 local to face) (same as Edge#8 in Hex neighbor x++,
      //                                    Edge#1 in Hex neighbor y++,
      //                                    Edge#3 in Hex neighbor x++, y++) 
      {
	edge_elemIDs[numEdge][edge_numElem[numEdge]] = i;
	edge_numElem[numEdge]++;
	elem_edgeIDs[i][5] = numEdge;
	ICoordscpy(hex_icoords, &temp_icoords);
	temp_icoords[0]++;
	temp_ID = getIndex(topoDim, hex_dims, temp_icoords);
	if (temp_ID > -1) {
	  edge_elemIDs[numEdge][edge_numElem[numEdge]] = temp_ID;
	  edge_numElem[numEdge]++;
	  elem_edgeIDs[temp_ID][8] = numEdge;
	}
	ICoordscpy(hex_icoords, &temp_icoords);
	temp_icoords[1]++;
	temp_ID = getIndex(topoDim, hex_dims, temp_icoords);
	if (temp_ID > -1) {
	  edge_elemIDs[numEdge][edge_numElem[numEdge]] = temp_ID;
	  edge_numElem[numEdge]++;
	  elem_edgeIDs[temp_ID][1] = numEdge;
	}
	ICoordscpy(hex_icoords, &temp_icoords);
	temp_icoords[0]++;
	temp_icoords[1]++;
	temp_ID = getIndex(topoDim, hex_dims, temp_icoords);
	if (temp_ID > -1) {
	  edge_elemIDs[numEdge][edge_numElem[numEdge]] = temp_ID;
	  edge_numElem[numEdge]++;
	  elem_edgeIDs[temp_ID][3] = numEdge;
	}
	edge_vertIDs[numEdge][0] = face_vertIDs[numFace][1];
	edge_vertIDs[numEdge][1] = face_vertIDs[numFace][2];
	numEdge++;
      }
      // Edge#6 (#2 local to face) (same as Edge#11 in Hex neighbor x++,
      //                                    Edge#4 in Hex neighbor z++,
      //                                    Edge#10 in Hex neighbor x++, z++) 
      {
	edge_elemIDs[numEdge][edge_numElem[numEdge]] = i;
	edge_numElem[numEdge]++;
	elem_edgeIDs[i][6] = numEdge;
	ICoordscpy(hex_icoords, &temp_icoords);
	temp_icoords[0]++;
	temp_ID = getIndex(topoDim, hex_dims, temp_icoords);
	if (temp_ID > -1) {
	  edge_elemIDs[numEdge][edge_numElem[numEdge]] = temp_ID;
	  edge_numElem[numEdge]++;
	  elem_edgeIDs[temp_ID][11] = numEdge;
	}
	ICoordscpy(hex_icoords, &temp_icoords);
	temp_icoords[2]++;
	temp_ID = getIndex(topoDim, hex_dims, temp_icoords);
	if (temp_ID > -1) {
	  edge_elemIDs[numEdge][edge_numElem[numEdge]] = temp_ID;
	  edge_numElem[numEdge]++;
	  elem_edgeIDs[temp_ID][4] = numEdge;
	}
	ICoordscpy(hex_icoords, &temp_icoords);
	temp_icoords[0]++;
	temp_icoords[2]++;
	temp_ID = getIndex(topoDim, hex_dims, temp_icoords);
	if (temp_ID > -1) {
	  edge_elemIDs[numEdge][edge_numElem[numEdge]] = temp_ID;
	  edge_numElem[numEdge]++;
	  elem_edgeIDs[temp_ID][10] = numEdge;
	}
	edge_vertIDs[numEdge][0] = face_vertIDs[numFace][2];
	edge_vertIDs[numEdge][1] = face_vertIDs[numFace][3];
	numEdge++;
      }
      numFace++;
    }
    // Face#2 (same as Face#0 in Hex neighbor with coords y++)
    {
      face_elemIDs[numFace][0] = i;
      face_numElem[numFace] = 1;
      elem_faceIDs[i][2] = numFace;
      elem_faceOrients[i][2] = 1;
      ICoordscpy(hex_icoords, &temp_icoords);
      temp_icoords[1]++;
      temp_ID = getIndex(topoDim, hex_dims, temp_icoords);
      if (temp_ID > -1) {
	face_elemIDs[numFace][1] = temp_ID;
	face_numElem[numFace]++;
	elem_faceIDs[temp_ID][0] = numFace;
	elem_faceOrients[temp_ID][0] = -1;
      }
      // Vertices
      ICoordscpy(hex_icoords, &temp_icoords);
      temp_icoords[0]++;
      temp_icoords[1]++;
      temp_ID = getIndex(topoDim, vert_dims, temp_icoords);
      face_vertIDs[numFace][0] = temp_ID;
      temp_icoords[0]--;
      temp_ID = getIndex(topoDim, vert_dims, temp_icoords);
      face_vertIDs[numFace][1] = temp_ID;
      temp_icoords[2]++;
      temp_ID = getIndex(topoDim, vert_dims, temp_icoords);
      face_vertIDs[numFace][2] = temp_ID;
      temp_icoords[0]++;
      temp_ID = getIndex(topoDim, vert_dims, temp_icoords);
      face_vertIDs[numFace][3] = temp_ID;
      // edges - Do 7-9
      // Edge#7 (#0 local to face) (same as edge#0 in hex neighbor y++) 
      if (hex_icoords[2]%hex_dims[2] == 0) {
	edge_elemIDs[numEdge][edge_numElem[numEdge]] = i;
	edge_numElem[numEdge]++;
	elem_edgeIDs[i][7] = numEdge;
	ICoordscpy(hex_icoords, &temp_icoords);
	temp_icoords[1]++;
	temp_ID = getIndex(topoDim, hex_dims, temp_icoords);
	if (temp_ID > -1) {
	  edge_elemIDs[numEdge][edge_numElem[numEdge]] = temp_ID;
	  edge_numElem[numEdge]++;
	  elem_edgeIDs[temp_ID][0] = numEdge;
	}
	edge_vertIDs[numEdge][0] = face_vertIDs[numFace][0];
	edge_vertIDs[numEdge][1] = face_vertIDs[numFace][1];
	numEdge++;
      }
      // Edge#8 (#1 local to face) (same as Edge#3 in Hex neighbor y++) 
      if (hex_icoords[0]%hex_dims[0] == 0) {
	edge_elemIDs[numEdge][edge_numElem[numEdge]] = i;
	edge_numElem[numEdge]++;
	elem_edgeIDs[i][8] = numEdge;
	ICoordscpy(hex_icoords, &temp_icoords);
	temp_icoords[1]++;
	temp_ID = getIndex(topoDim, hex_dims, temp_icoords);
	if (temp_ID > -1) {
	  edge_elemIDs[numEdge][edge_numElem[numEdge]] = temp_ID;
	  edge_numElem[numEdge]++;
	  elem_edgeIDs[temp_ID][3] = numEdge;
	}
	edge_vertIDs[numEdge][0] = face_vertIDs[numFace][1];
	edge_vertIDs[numEdge][1] = face_vertIDs[numFace][2];
	numEdge++;
      }
      // Edge#9 (#2 local to face) (same as Edge#2 in Hex neighbor y++,
      //                                    Edge#7 in Hex neighbor z++,
      //                                    Edge#0 in Hex neighbor y++, z++) 
      {
	edge_elemIDs[numEdge][edge_numElem[numEdge]] = i;
	edge_numElem[numEdge]++;
	elem_edgeIDs[i][9] = numEdge;
	ICoordscpy(hex_icoords, &temp_icoords);
	temp_icoords[1]++;
	temp_ID = getIndex(topoDim, hex_dims, temp_icoords);
	if (temp_ID > -1) {
	  edge_elemIDs[numEdge][edge_numElem[numEdge]] = temp_ID;
	  edge_numElem[numEdge]++;
	  elem_edgeIDs[temp_ID][2] = numEdge;
	}
	ICoordscpy(hex_icoords, &temp_icoords);
	temp_icoords[2]++;
	temp_ID = getIndex(topoDim, hex_dims, temp_icoords);
	if (temp_ID > -1) {
	  edge_elemIDs[numEdge][edge_numElem[numEdge]] = temp_ID;
	  edge_numElem[numEdge]++;
	  elem_edgeIDs[temp_ID][7] = numEdge;
	}
	ICoordscpy(hex_icoords, &temp_icoords);
	temp_icoords[1]++;
	temp_icoords[2]++;
	temp_ID = getIndex(topoDim, hex_dims, temp_icoords);
	if (temp_ID > -1) {
	  edge_elemIDs[numEdge][edge_numElem[numEdge]] = temp_ID;
	  edge_numElem[numEdge]++;
	  elem_edgeIDs[temp_ID][0] = numEdge;
	}
	edge_vertIDs[numEdge][0] = face_vertIDs[numFace][2];
	edge_vertIDs[numEdge][1] = face_vertIDs[numFace][3];
	numEdge++;
      }
      numFace++;
    }
    // Face#3 - sometimes
    if (hex_icoords[0]%hex_dims[0] == 0) {
      face_elemIDs[numFace][0] = i;
      face_numElem[numFace] = 1;
      elem_faceIDs[i][3] = numFace;
      elem_faceOrients[i][3] = 1;
      // Vertices
      ICoordscpy(hex_icoords, &temp_icoords);
      temp_icoords[1]++;
      temp_ID = getIndex(topoDim, vert_dims, temp_icoords);
      face_vertIDs[numFace][0] = temp_ID;
      temp_icoords[1]--;
      temp_ID = getIndex(topoDim, vert_dims, temp_icoords);
      face_vertIDs[numFace][1] = temp_ID;
      temp_icoords[2]++;
      temp_ID = getIndex(topoDim, vert_dims, temp_icoords);
      face_vertIDs[numFace][2] = temp_ID;
      temp_icoords[1]++;
      temp_ID = getIndex(topoDim, vert_dims, temp_icoords);
      face_vertIDs[numFace][3] = temp_ID;
      // edges - Do 10&11
      // Edge#10 (#0 local to face)
      if (hex_icoords[0]%hex_dims[0] == 0 && hex_icoords[2]%hex_dims[2] == 0) {
	edge_elemIDs[numEdge][edge_numElem[numEdge]] = i;
	edge_numElem[numEdge]++;
	elem_edgeIDs[i][10] = numEdge;
	edge_vertIDs[numEdge][0] = face_vertIDs[numFace][0];
	edge_vertIDs[numEdge][1] = face_vertIDs[numFace][1];
	numEdge++;
      }
      // Edge#11 (#2 local to face) (same as Edge#10 in Hex neighbor z++)
      if (hex_icoords[0]%hex_dims[0] == 0) {
	edge_elemIDs[numEdge][edge_numElem[numEdge]] = i;
	edge_numElem[numEdge]++;
	elem_edgeIDs[i][11] = numEdge;
	ICoordscpy(hex_icoords, &temp_icoords);
	temp_icoords[2]++;
	temp_ID = getIndex(topoDim, hex_dims, temp_icoords);
	if (temp_ID > -1) {
	  edge_elemIDs[numEdge][edge_numElem[numEdge]] = temp_ID;
	  edge_numElem[numEdge]++;
	  elem_edgeIDs[temp_ID][10] = numEdge;
	}
	edge_vertIDs[numEdge][0] = face_vertIDs[numFace][2];
	edge_vertIDs[numEdge][1] = face_vertIDs[numFace][3];
	numEdge++;
      }
      numFace++;
    }
    // Face#4 - sometimes
    if (hex_icoords[2]%hex_dims[2] == 0) {
      face_elemIDs[numFace][0] = i;
      face_numElem[numFace] = 1;
      elem_faceIDs[i][4] = numFace;
      elem_faceOrients[i][4] = 1;
      // Vertices
      ICoordscpy(hex_icoords, &temp_icoords);
      temp_ID = getIndex(topoDim, vert_dims, temp_icoords);
      face_vertIDs[numFace][0] = temp_ID;
      temp_icoords[1]++;
      temp_ID = getIndex(topoDim, vert_dims, temp_icoords);
      face_vertIDs[numFace][1] = temp_ID;
      temp_icoords[0]++;
      temp_ID = getIndex(topoDim, vert_dims, temp_icoords);
      face_vertIDs[numFace][2] = temp_ID;
      temp_icoords[1]--;
      temp_ID = getIndex(topoDim, vert_dims, temp_icoords);
      face_vertIDs[numFace][3] = temp_ID;
      // no edges
      numFace++;
    }
    // Face#5 (same as Face#4 in Hex neighbor with coords z++)
    {
      face_elemIDs[numFace][0] = i;
      face_numElem[numFace] = 1;
      elem_faceIDs[i][5] = numFace;
      elem_faceOrients[i][5] = 1;
      ICoordscpy(hex_icoords, &temp_icoords);
      temp_icoords[2]++;
      temp_ID = getIndex(topoDim, hex_dims, temp_icoords);
      if (temp_ID > -1) {
	face_elemIDs[numFace][1] = temp_ID;
	face_numElem[numFace]++;
	elem_faceIDs[temp_ID][4] = numFace;
	elem_faceOrients[temp_ID][4] = -1;
      }
      // Vertices
      ICoordscpy(hex_icoords, &temp_icoords);
      temp_icoords[2]++;
      temp_ID = getIndex(topoDim, vert_dims, temp_icoords);
      face_vertIDs[numFace][0] = temp_ID;
      temp_icoords[0]++;
      temp_ID = getIndex(topoDim, vert_dims, temp_icoords);
      face_vertIDs[numFace][1] = temp_ID;
      temp_icoords[1]++;
      temp_ID = getIndex(topoDim, vert_dims, temp_icoords);
      face_vertIDs[numFace][2] = temp_ID;
      temp_icoords[0]--;
      temp_ID = getIndex(topoDim, vert_dims, temp_icoords);
      face_vertIDs[numFace][3] = temp_ID;
      // no edges
      numFace++;
    }
  }

  // quick check that build makes sense
  fail_unless(numFace == numFace_good,
	      "failed to build faces for testing, bug in test code");
  fail_unless(numEdge == numEdge_good, 
	      "failed to build edges for testing, bug in test code");
  
  // sort edges' vertices
  for (i = 0; i < numEdge_good; i++) {
    if (edge_vertIDs[i][1] < edge_vertIDs[i][0]) {
      int temp_int = edge_vertIDs[i][0];
      edge_vertIDs[i][0] = edge_vertIDs[i][1];
      edge_vertIDs[i][1] = temp_int;
    }
  }
  // sort & reorder faces' vertices 
  for (i = 0; i < numFace_good; i++) {
    int temp_vertIDs[4];
    int idx = 0;
    int min = face_vertIDs[i][0];
    for (j = 1; j < 4; j++) {
      if (face_vertIDs[i][j] < min) {
	min = face_vertIDs[i][j];
	idx = j;
      }
    }
    // "rotate" to get smallest vertID first
    for (j = idx; j < 4; j++) 
      temp_vertIDs[j-idx] = face_vertIDs[i][j];
    for (j = 0; j < idx; j++)
      temp_vertIDs[j+(4-idx)] = face_vertIDs[i][j];
    // swap order to get second smallest vertID next
    if (temp_vertIDs[3] < temp_vertIDs[1]) {
      face_vertIDs[i][0] = temp_vertIDs[0];
      for (j = 1; j < 4; j++) // loop over vertices
	face_vertIDs[i][j] = temp_vertIDs[4-j];
      for (j = 0; j < face_numElem[i]; j++) {  
	int elemID = face_elemIDs[i][j];
	int relative_faceID = -1;
	for (k = 0; k < numFacePerHex; k++) 
	  if (elem_faceIDs[elemID][k] == i)
	    relative_faceID = k;
	fail_unless(relative_faceID > -1, "didn't find faceID, bug in test code");
	elem_faceOrients[elemID][relative_faceID] *= -1;
      }
    }
    else 
      for (j = 0; j < 4; j++)
	face_vertIDs[i][j] = temp_vertIDs[j];
  }

  // vertices
  fc_getMeshNumVertex(mesh, &numVertex);
  fail_unless(numVertex == numVertex_good, "mismatch number of vertices");
  rc = _fc_getMeshEdgeParentsOfVerticesPtr(mesh, &numEdgePerVert,
				       &edgeParentsPerVert); 
  fail_unless(rc == FC_SUCCESS, "failed to get vertices' edge parents");
  rc = _fc_getMeshFaceParentsOfVerticesPtr(mesh, &numFacePerVert,
				       &faceParentsPerVert); 
  fail_unless(rc == FC_SUCCESS, "failed to get vertices' face parents");
  rc = _fc_getMeshElementParentsOfVerticesPtr(mesh, &numElemPerVert,
				       &elemParentsPerVert); 
  fail_unless(rc == FC_SUCCESS, "failed to get vertices' elem parents");
  rc = _fc_getMeshVertexNeighborsViaEdgePtr(mesh, &numVertNeighsViaEdge,
					    &vertNeighsViaEdge);
  fail_unless(rc == FC_SUCCESS, "failed to get vert neighs via edge");
  for (i = 0; i < numVertex; i++) {
    // calculate stuff
    // setup 
    ICoords vert_icoords, temp_icoords;
    int temp_ID;
    int maxNumElem = 8;
    int actualNumElem;
    int actualNumFace;
    int actualNumEdge;
    int* elemIDs;
    int* faceOwnerIDs;
    int* edgeOwnerIDs;
    FC_SortedIntArray elemOwnerList;
    FC_SortedIntArray faceOwnerList;
    FC_SortedIntArray edgeOwnerList;
    int maxNumEdgeNeigh = 6;
    int actualNumEdgeNeigh;
    int edgeNeighIDs[6];
    // this stencil adjusts owning hex coords (from the vertex coords);
    int elemOwnerStencil[8][3] = 
      { { -1, -1, -1 }, { 0, -1, -1 }, { -1, 0, -1 }, { 0, 0, -1 },
	{ -1, -1,  0 }, { 0, -1,  0 }, { -1, 0,  0 }, { 0, 0,  0 } };
    int faceLU[8][3] = 
      { { 1, 2, 5 }, { 2, 3, 5 }, { 0, 1, 5 }, { 0, 3, 5 },
	{ 1, 2, 4 }, { 2, 3, 4 }, { 0, 1, 4 }, { 0, 3, 4 } } ;
    int edgeLU[8][3] = 
      { { 5, 6, 9 }, { 8, 9, 11 }, { 1, 2, 6 }, { 2, 3, 11 },
	{ 4, 5, 7 }, { 7, 8, 10 }, { 0, 1, 4 }, { 0, 3, 10 } };
    // this stencil adjusts vertex ijk coords
    int edgeNeighStencil[6][3] = 
      { {  0,  0, -1 }, 
	{  0, -1,  0 }, { -1, 0, 0 }, {  1, 0, 0 }, { 0, 1, 0 },
	{  0,  0,  1 } };
    // setup II - shared calcs
    getICoords(i, topoDim, vert_dims, &vert_icoords);
    // calc owners
    fc_initSortedIntArray(&elemOwnerList);
    fc_initSortedIntArray(&faceOwnerList);
    fc_initSortedIntArray(&edgeOwnerList);
    for (j = 0; j < maxNumElem; j++) {
      for (k = 0; k < 3; k++)
	temp_icoords[k] = vert_icoords[k] + elemOwnerStencil[j][k];
      temp_ID = getIndex(topoDim, hex_dims, temp_icoords);
      if (temp_ID > -1) {
	fc_addIntToSortedIntArray(&elemOwnerList, temp_ID);
	for (k = 0; k < 3; k++) {
	  fc_addIntToSortedIntArray(&faceOwnerList, 
				    elem_faceIDs[temp_ID][faceLU[j][k]]);
	  fc_addIntToSortedIntArray(&edgeOwnerList,
				    elem_edgeIDs[temp_ID][edgeLU[j][k]]);	
	}
      }
    }
    fc_convertSortedIntArrayToIntArray(&elemOwnerList, &actualNumElem, &elemIDs);
    fc_convertSortedIntArrayToIntArray(&faceOwnerList, &actualNumFace, &faceOwnerIDs);
    fc_convertSortedIntArrayToIntArray(&edgeOwnerList, &actualNumEdge, &edgeOwnerIDs);
    // calc edge neighbors
    actualNumEdgeNeigh = 0;
    for (j = 0; j < maxNumEdgeNeigh; j++) {
      for (k = 0; k < 3; k++)
	temp_icoords[k] = vert_icoords[k] + edgeNeighStencil[j][k];
      temp_ID = getIndex(topoDim, vert_dims, temp_icoords);
      if (temp_ID > -1) {
	edgeNeighIDs[actualNumEdgeNeigh] = temp_ID;
	actualNumEdgeNeigh++;
      }
    }

    // tests
    fail_unless(numEdgePerVert[i] == actualNumEdge,
		"mismatch in numEdgePerVert");
    fail_unless(!memcmp(edgeParentsPerVert[i], edgeOwnerIDs,
			numEdgePerVert[i]*sizeof(int)), 
		"mismatch in edgeParentsPerVert");
    free(edgeOwnerIDs);
    fail_unless(numFacePerVert[i] == actualNumFace, 
		"mismatch in numFacePerVert");
    fail_unless(!memcmp(faceParentsPerVert[i], faceOwnerIDs,
			numFacePerVert[i]*sizeof(int)), 
		"mismatch in faceParentsPerVert");
    free(faceOwnerIDs);
    fail_unless(numElemPerVert[i] == actualNumElem, 
		"vertex: mismatch in number of elements");
    fail_unless(!memcmp(elemParentsPerVert[i], elemIDs, 
			numElemPerVert[i]*sizeof(int)),
		"vertex: mismatch in elementIDs");
    free(elemIDs);
    fail_unless(numVertNeighsViaEdge[i] == actualNumEdgeNeigh,
		"mismatch of numVertNeighsViaEdge");
    fail_unless(!memcmp(vertNeighsViaEdge[i], edgeNeighIDs,
			numVertNeighsViaEdge[i]*sizeof(int)), 
		"mismatch of vertNeighsViaEdge");
  }

  // edges
  fc_getMeshNumEdge(mesh, &numEdge);
  fail_unless(numEdge == numEdge_good, "mismatch in number edges");
  rc = fc_getMeshEdgeConnsPtr(mesh, &edgeConns);
  fail_unless(rc == FC_SUCCESS, "failed to get edge conns");
  rc = _fc_getMeshElementParentsOfEdgesPtr(mesh, &numElemPerEdge,
					&elemParentsPerEdge);
  fail_unless(rc == FC_SUCCESS, "failed to get edge's elem parents");
  for (i = 0; i < numEdge; i++) {
    // tests
    fail_unless(!memcmp(&edgeConns[i*2], edge_vertIDs[i], 2*sizeof(int)),
		"edge: mismatch of vertIDs");
    fail_unless(numElemPerEdge[i] == edge_numElem[i], 
		"mismatch in numElemPerEdge");
    fail_unless(!memcmp(elemParentsPerEdge[i], edge_elemIDs[i],
			numElemPerEdge[i]*sizeof(int)), 
		"mismatch of elemParentsPerEdge");
  }

  // faces
  fc_getMeshNumFace(mesh, &numFace);
  fail_unless(numFace == numFace_good, "mismatch in number faces");
  rc = fc_getMeshFaceTypesPtr(mesh, &faceTypes);
  fail_unless(rc == FC_SUCCESS, "failed to get face types");
  rc = fc_getMeshFaceConnsPtr(mesh, &numVertPerFace, &maxNumVertPerFace,
			      &faceConns);
  fail_unless(rc == FC_SUCCESS, "failed to get face conns");
  rc = _fc_getMeshElementParentsOfFacesPtr(mesh, &numElemPerFace,
					&elemParentsPerFace);
  fail_unless(rc == FC_SUCCESS, "failed to get face's elem parents");
  for (i = 0; i < numFace; i++) {
    // tests
    fail_unless(faceTypes[i] == FC_ET_QUAD, "face: mismatch in face type");
    fail_unless(numVertPerFace[i] == 4, 
		"face: mismatch in number of vertices");
    fail_unless(!memcmp(&faceConns[i*maxNumVertPerFace], face_vertIDs[i],
			 4*sizeof(int)), "face: mismatch in vertex IDs");
    fail_unless(numElemPerFace[i] == face_numElem[i], 
		"mismatch in numElemPerFace");
    fail_unless(!memcmp(elemParentsPerFace[i], face_elemIDs[i],
			numElemPerFace[i]*sizeof(int)), 
		"mismatch in elemParentsPerFace");
  }

  // elements
  fc_getMeshNumElement(mesh, &numElement);
  fail_unless(numElement == numElement_good, "mismatch in number of elements");
  rc = fc_getMeshElementToEdgeConnsPtr(mesh, &elemToEdgeConns);
  fail_unless(rc == FC_SUCCESS, "failed to get elem to edge conns");
  rc = fc_getMeshElementToFaceConnsPtr(mesh, &elemToFaceConns);
  fail_unless(rc == FC_SUCCESS, "failed to get elem to edge conns");
  rc = fc_getMeshElementFaceOrientsPtr(mesh, &elemFaceOrients);
  fail_unless(rc == FC_SUCCESS, "failed to get elem to edge conns");
  rc = _fc_getMeshElementNeighborsViaVertexPtr(mesh, &numElemNeighsViaVert,
					       &elemNeighsViaVert);
  fail_unless(rc == FC_SUCCESS, "failed to get elemNeighsViaVert");
  rc = _fc_getMeshElementNeighborsViaEdgePtr(mesh, &numElemNeighsViaEdge,
					     &elemNeighsViaEdge);
  fail_unless(rc == FC_SUCCESS, "failed to get elemNeighsViaEdge");
  rc = _fc_getMeshElementNeighborsViaFacePtr(mesh, &numElemNeighsViaFace,
					     &elemNeighsViaFace);
  fail_unless(rc == FC_SUCCESS, "failed to get elemNeighsViaFace");
  for (i = 0; i < numElement; i++) {
    // calculate stuff
    // setup
    int hex_index;
    ICoords hex_icoords, temp_icoords;
    int numVert = 8;
    int vertIDs[8];
    int maxNumVertNeigh = 26;  // 6 + 12 + 8
    int maxNumEdgeNeigh = 18;  // 6 + 12
    int maxNumFaceNeigh = 6;
    int actualNumVertNeigh;
    int actualNumEdgeNeigh;
    int actualNumFaceNeigh;
    int vertNeighIDs[26]; 
    int edgeNeighIDs[18];
    int faceNeighIDs[6];
    // this stencil adjusts ijk coords for quads
    int vertNeighStencil[26][3] = {
      // bottom 9
      { -1, -1, -1 }, {  0, -1, -1 }, {  1, -1, -1 },
      { -1,  0, -1 }, {  0,  0, -1 }, {  1,  0, -1 },
      { -1,  1, -1 }, {  0,  1, -1 }, {  1,  1, -1 },
      // middle 8
      { -1, -1,  0 }, {  0, -1,  0 }, {  1, -1,  0 },
      { -1,  0,  0 },                 {  1,  0,  0 },
      { -1,  1 , 0 }, {  0,  1,  0 }, {  1,  1,  0 },
      // top 9
      { -1, -1,  1 }, {  0, -1,  1 }, {  1, -1,  1 },
      { -1,  0,  1 }, {  0,  0,  1 }, {  1,  0,  1 },
      { -1,  1,  1 }, {  0,  1,  1 }, {  1,  1,  1 } };
    int edgeNeighStencil[18][3] = { 
      // bottom 5
      {  0, -1, -1 },
      { -1,  0, -1 }, {  0,  0, -1 }, {  1,  0, -1 },
      {  0,  1, -1 },
      // middle 8
      { -1, -1,  0 }, {  0, -1,  0 }, {  1, -1,  0 },
      { -1,  0,  0 },                 {  1,  0,  0 },
      { -1,  1 , 0 }, {  0,  1,  0 }, {  1,  1,  0 },
      // top 5
      {  0, -1,  1 },
      { -1,  0,  1 }, {  0,  0,  1 }, {  1,  0,  1 },
      {  0,  1,  1 } };
    int faceNeighStencil[6][3] = {
      // bottom 1 
      {  0,  0, -1 },
      // middle 4 
      {  0, -1,  0 }, { -1,  0,  0 },
      {  1,  0,  0 }, {  0,  1,  0 }, 
      // top 1
      {  0,  0,  1 } };
    // setup II - shared calcs      
    getICoords(i, topoDim, hex_dims, &hex_icoords);
    // owned vertices
    for (j = 0; j < numVert; j++) { // walk around hex
      switch (j) {
      case 0: // start of bottom face
	ICoordscpy(hex_icoords, &temp_icoords);
	break;
      case 4: // start of top face
	ICoordscpy(hex_icoords, &temp_icoords);
	temp_icoords[2]++;
	break;
      case 1:  
      case 5: 
	temp_icoords[0]++;  break;
      case 2:
      case 6:
	temp_icoords[1]++;  break;
      case 3:  
      case 7:
	temp_icoords[0]--;  break;
      }
      vertIDs[j] = getIndex(topoDim, vert_dims, temp_icoords);
    }
    // vertex neighbors
    actualNumVertNeigh = 0;
    for (j = 0; j < maxNumVertNeigh; j++) {
      for (k = 0; k < 3; k++)
      temp_icoords[k] = hex_icoords[k] + vertNeighStencil[j][k];
      hex_index = getIndex(topoDim, hex_dims, temp_icoords);
      if (hex_index > -1) {
	vertNeighIDs[actualNumVertNeigh] = hex_index;
	actualNumVertNeigh++;
      }
    }
    // edge neighbors
    actualNumEdgeNeigh = 0;
    for (j = 0; j < maxNumEdgeNeigh; j++) {
      for (k = 0; k < 3; k++)
	temp_icoords[k] = hex_icoords[k] + edgeNeighStencil[j][k];
      hex_index = getIndex(topoDim, hex_dims, temp_icoords);
      if (hex_index > -1) {
	edgeNeighIDs[actualNumEdgeNeigh] = hex_index;
	actualNumEdgeNeigh++;
      }
    }
    // face neighbors
    actualNumFaceNeigh = 0;
    for (j = 0; j < maxNumFaceNeigh; j++) {
      for (k = 0; k < 3; k++)
	temp_icoords[k] = hex_icoords[k] + faceNeighStencil[j][k];
      hex_index = getIndex(topoDim, hex_dims, temp_icoords);
      if (hex_index > -1) {
	faceNeighIDs[actualNumFaceNeigh] = hex_index;
	actualNumFaceNeigh++;
      }
    }

    // element tests
    fail_unless(!memcmp(&elemToEdgeConns[i*numEdgePerHex], elem_edgeIDs[i],
			numEdgePerHex*sizeof(int)), 
		"element: mismatch in edgeIDs");
    fail_unless(!memcmp(&elemToFaceConns[i*numFacePerHex], elem_faceIDs[i],
			numFacePerHex*sizeof(int)),
		"element: mismatch in face IDs");
    fail_unless(!memcmp(&elemFaceOrients[i*numFacePerHex], elem_faceOrients[i],
			numFacePerHex*sizeof(int)),
		"element: mismatch in face orients");
    fail_unless(numElemNeighsViaVert[i] == actualNumVertNeigh,
		"mismatch of numElemNeighsViaVert");
    fail_unless(!memcmp(elemNeighsViaVert[i], vertNeighIDs,
			actualNumVertNeigh*sizeof(int)), 
		"element: mismatch in vertex neighbors");
    fail_unless(numElemNeighsViaEdge[i] == actualNumEdgeNeigh,
		"mismatch of numElemNeighsViaEdge");
    fail_unless(!memcmp(elemNeighsViaEdge[i], edgeNeighIDs,
			actualNumEdgeNeigh*sizeof(int)),
		"mismatch of elemNeighsViaEdge");
    fail_unless(numElemNeighsViaFace[i] == actualNumFaceNeigh, 
		"mismatch of numElemNeighsViaFace");
    fail_unless(!memcmp(elemNeighsViaFace[i], faceNeighIDs,
			actualNumFaceNeigh*sizeof(int)),
		"mismatch of elemNeighsViaFace");
  }

  fc_deleteDataset(dataset);
}
END_TEST
  
// test building edges, faces & neighbors
START_TEST(build_prism_mesh)
{
  FC_ReturnCode rc;
  int i, j, k;
  int topoDim,numDim, numVertex, numEdge, numFace, numElement;
  int numVertex_good, numEdge_good, numFace_good, numElement_good;
  int numEdgePerPrism = 9;
  int numFacePerPrism = 5;
  FC_Dataset dataset;
  FC_Mesh mesh, *returnMeshes;
  int numReturnMeshes;
  FC_ElementType elemType, *faceTypes;
  int* edgeConns, *faceConns, maxNumVertPerFace, *numVertPerFace;
  int* elemToEdgeConns, *elemToFaceConns, *elemFaceOrients;
  int* numEdgePerVert, **edgeParentsPerVert;
  int* numFacePerVert, **faceParentsPerVert;
  int* numElemPerVert, **elemParentsPerVert;
  int* numElemPerEdge, **elemParentsPerEdge;
  int* numElemPerFace, **elemParentsPerFace;
  int* numVertNeighsViaEdge, **vertNeighsViaEdge;
  int* numElemNeighsViaVert, **elemNeighsViaVert;
  int* numElemNeighsViaEdge, **elemNeighsViaEdge;
  int* numElemNeighsViaFace, **elemNeighsViaFace;
  int vert_dims[3] = { M + 1, N + 1, P + 1 };
  int hex_dims[3] = { M, N, P };
  int edge_vertIDs[4*M*N*P + 3*M*N + 2*N*P + 2*M*P + M + N + P][2]; 
  int edge_numElem[4*M*N*P + 3*M*N + 2*N*P + 2*M*P + M + N + P];  
  int edge_elemIDs[4*M*N*P + 3*M*N + 2*N*P + 2*M*P + M + N + P][6];
  int face_numVert[5*M*N*P + 2*M*N + N*P + M*P];
  int face_vertIDs[5*M*N*P + 2*M*N + N*P + M*P][4];
  int face_numElem[5*M*N*P + 2*M*N + N*P + M*P];
  int face_elemIDs[5*M*N*P + 2*M*N + N*P + M*P][2];
  int elem_edgeIDs[2*M*N*P][9]; // number of elements
  int elem_faceIDs[2*M*N*P][5];
  int elem_faceOrients[2*M*N*P][5];

  // pre calc "good" values = hex values plus some for extra face in center
  numVertex_good = (M+1)*(N+1)*(P+1);
  numEdge_good = 3*M*N*P + 2*M*N + 2*N*P + 2*M*P + M + N + P + (M*N*P + M*N);
  numFace_good = 3*M*N*P + M*N + N*P + M*P + (2*M*N*P + M*N);
  numElement_good = 2*M*N*P;

  // setup
  fc_loadDataset("../data/gen_prism.ex2", &dataset);
  rc = fc_getMeshByName(dataset, "prism mesh", &numReturnMeshes,&returnMeshes);
  fail_unless(rc == FC_SUCCESS, "Test aborted: mesh not found");
  fail_unless(numReturnMeshes == 1, "failed to find unique mesh by name");
  mesh = returnMeshes[0];
  free(returnMeshes);

  rc = fc_getMeshInfo(mesh, &topoDim, &numDim, &numVertex, &numElement,
			&elemType);
  fail_unless(rc == FC_SUCCESS, "Test aborted: mesh could not be described\n");
  fail_unless(topoDim == 3, "Text aborted: expected mesh with topodim of 3");
  fail_unless(elemType == FC_ET_PRISM, "Test aborted: expected hex elements");
  fail_unless(numElement == numElement_good, "bad");
  
  // build edges & faces info - vertex members, and back and forth
  // connections to elems
  //
  // Basic procedure: walk every element, then every face of that element,
  // then every edge, making new faces and edges as we go.  Because edges and
  // faces are shared, we don't make every face and edge of every element. The
  // idea is that "downstream" faces and edges of a hex (like the top of the
  // hex) are always made and then pushed onto downstream hexes. The only time
  // "upstream" edges  and faces are made is if they are on the boundary of
  // the mesh. Procedurally, roughly
  //    for each element { 
  //      do the faces {
  //        skip upstream faces if element not on the boundary
  //        make face on this element
  //        add face to downstream elements
  //        do the edges {
  //          skip upstream edges if element not on the boundary
  //          make edge
  //          add edge to downstream elements
  //        }
  //      }
  //    }
  // It's probably helpful to examine mesh.c's fc_buildFaces and fc_buildEdges
  // to get an idea of the order used.
  numEdge = 0; 
  numFace = 0;
  for (i = 0; i < numEdge_good; i++) {
    edge_numElem[i] = 0;
  }
  for (i = 0; i < numElement_good; i++)
    for (j = 0; j < numEdgePerPrism; j++)
      elem_edgeIDs[i][j] = -1;
  // loop over underlying hex's
  for (i = 0; i < numElement_good/2; i++) { // is over hex coords
    // setup
    int temp_ID;
    ICoords hex_icoords, temp_icoords;
    getICoords(i, topoDim, hex_dims, &hex_icoords);

    // Prism#0Face#0 - sometimes
    if (hex_icoords[1]%hex_dims[1] == 0) {
      face_elemIDs[numFace][0] = i*2;
      face_numElem[numFace] = 1;
      elem_faceIDs[i*2][0] = numFace;
      elem_faceOrients[i*2][0] = 1;
      // Vertices
      face_numVert[numFace] = 4;
      ICoordscpy(hex_icoords, &temp_icoords);
      temp_ID = getIndex(topoDim, vert_dims, temp_icoords);
      face_vertIDs[numFace][0] = temp_ID;
      temp_icoords[0]++;
      temp_ID = getIndex(topoDim, vert_dims, temp_icoords);
      face_vertIDs[numFace][1] = temp_ID;
      temp_icoords[2]++;
      temp_ID = getIndex(topoDim, vert_dims, temp_icoords);
      face_vertIDs[numFace][2] = temp_ID;
      temp_icoords[0]--;
      temp_ID = getIndex(topoDim, vert_dims, temp_icoords);
      face_vertIDs[numFace][3] = temp_ID;
      // edges - Do 0-3
      // Edge#0 (#0 local to face)
      if (hex_icoords[1]%hex_dims[1] == 0 && hex_icoords[2]%hex_dims[2] == 0) {
	edge_elemIDs[numEdge][edge_numElem[numEdge]] = i*2;
	edge_numElem[numEdge]++;
	elem_edgeIDs[i*2][0] = numEdge;
	edge_vertIDs[numEdge][0] = face_vertIDs[numFace][0]; 
	edge_vertIDs[numEdge][1] = face_vertIDs[numFace][1];
	numEdge++;
      }
      // Edge#1 (#1 local to face) (same as 3 in prism2,
      //                            same as Edge#3 in neighbor0 x++) 
      if (hex_icoords[1]%hex_dims[1] == 0) {
	edge_elemIDs[numEdge][edge_numElem[numEdge]] = i*2;
	edge_numElem[numEdge]++;
	elem_edgeIDs[i*2][1] = numEdge;
	edge_elemIDs[numEdge][edge_numElem[numEdge]] = i*2+1;
	edge_numElem[numEdge]++;
	elem_edgeIDs[i*2+1][3] = numEdge;
	ICoordscpy(hex_icoords, &temp_icoords);
	temp_icoords[0]++;
	temp_ID = getIndex(topoDim, hex_dims, temp_icoords);
	if (temp_ID > -1) {
	  edge_elemIDs[numEdge][edge_numElem[numEdge]] = temp_ID*2;
	  edge_numElem[numEdge]++;
	  elem_edgeIDs[temp_ID*2][3] = numEdge;
	}
	edge_vertIDs[numEdge][0] = face_vertIDs[numFace][1];
	edge_vertIDs[numEdge][1] = face_vertIDs[numFace][2];
	numEdge++;
      }
      // Edge#2 (#2 local to face) (same as Edge#0 in neighbor0 z++)
      if (hex_icoords[1]%hex_dims[1] == 0) {
	edge_elemIDs[numEdge][edge_numElem[numEdge]] = i*2;
	edge_numElem[numEdge]++;
	elem_edgeIDs[i*2][2] = numEdge;
	ICoordscpy(hex_icoords, &temp_icoords);
	temp_icoords[2]++;
	temp_ID = getIndex(topoDim, hex_dims, temp_icoords);
	if (temp_ID > -1) {
	  edge_elemIDs[numEdge][edge_numElem[numEdge]] = temp_ID*2;
	  edge_numElem[numEdge]++;
	  elem_edgeIDs[temp_ID*2][0] = numEdge;
	}
	edge_vertIDs[numEdge][0] = face_vertIDs[numFace][2];
	edge_vertIDs[numEdge][1] = face_vertIDs[numFace][3];
	numEdge++;
      }
      // Edge#3 (#3 local to face)
      if (hex_icoords[0]%hex_dims[0] == 0 && hex_icoords[1]%hex_dims[1] == 0) {
	edge_elemIDs[numEdge][edge_numElem[numEdge]] = i*2;
	edge_numElem[numEdge]++;
	elem_edgeIDs[i*2][3] = numEdge;	
	edge_vertIDs[numEdge][0] = face_vertIDs[numFace][3];
	edge_vertIDs[numEdge][1] = face_vertIDs[numFace][0];
	numEdge++;
      }
      numFace++;
    }
    // Prism#0Face#1 - same as face 2 in prism2
    {
      face_elemIDs[numFace][0] = i*2;
      face_numElem[numFace] = 1;
      elem_faceIDs[i*2][1] = numFace;
      elem_faceOrients[i*2][1] = 1;
      face_elemIDs[numFace][1] = i*2+1;
      face_numElem[numFace]++;
      elem_faceIDs[i*2+1][2] = numFace;
      elem_faceOrients[i*2+1][2] = -1;
      // vertices
      face_numVert[numFace] = 4;
      ICoordscpy(hex_icoords, &temp_icoords);
      temp_icoords[0]++;
      temp_ID = getIndex(topoDim, vert_dims, temp_icoords);
      face_vertIDs[numFace][0] = temp_ID;
      temp_icoords[0]--;
      temp_icoords[1]++;
      temp_ID = getIndex(topoDim, vert_dims, temp_icoords);
      face_vertIDs[numFace][1] = temp_ID;
      temp_icoords[2]++;
      temp_ID = getIndex(topoDim, vert_dims, temp_icoords);
      face_vertIDs[numFace][2] = temp_ID;
      temp_icoords[0]++;
      temp_icoords[1]--;
      temp_ID = getIndex(topoDim, vert_dims, temp_icoords);
      face_vertIDs[numFace][3] = temp_ID;
      // edges 4-6
      // Edge#4 (0 local to face) (same as edge 7 in prism2)
      if (hex_icoords[2]%hex_dims[2] == 0) {
	edge_elemIDs[numEdge][edge_numElem[numEdge]] = i*2;
	edge_numElem[numEdge]++;
	elem_edgeIDs[i*2][4] = numEdge;
	edge_elemIDs[numEdge][edge_numElem[numEdge]] = i*2+1;
	edge_numElem[numEdge]++;
	elem_edgeIDs[i*2+1][7] = numEdge;
	edge_vertIDs[numEdge][0] = face_vertIDs[numFace][0];
	edge_vertIDs[numEdge][1] = face_vertIDs[numFace][1];
	numEdge++;
      }
      // Edge #5 (1 local to face) (same as edge 5 in prism2,
      //                            same as edge 3 in neighbor0 y++,)
      if (hex_icoords[0]%hex_dims[0] == 0) {
	edge_elemIDs[numEdge][edge_numElem[numEdge]] = i*2;
	edge_numElem[numEdge]++;
	elem_edgeIDs[i*2][5] = numEdge;
	edge_elemIDs[numEdge][edge_numElem[numEdge]] = i*2+1;
	edge_numElem[numEdge]++;
	elem_edgeIDs[i*2+1][5] = numEdge;
	ICoordscpy(hex_icoords, &temp_icoords);
	temp_icoords[1]++;
	temp_ID = getIndex(topoDim, hex_dims, temp_icoords);
	if (temp_ID > -1) {
	  edge_elemIDs[numEdge][edge_numElem[numEdge]] = temp_ID*2;
	  edge_numElem[numEdge]++;
	  elem_edgeIDs[temp_ID*2][3] = numEdge;
	}
	edge_vertIDs[numEdge][0] = face_vertIDs[numFace][1];
	edge_vertIDs[numEdge][1] = face_vertIDs[numFace][2];
	numEdge++;
      }
      // Edge#6 (2 local to face) (same as edge 8 in prism2,
      //                           same as edge 4 in neighbor0 z++,
      //                           same as edge 7 in neighbor1 z++)
      {
	edge_elemIDs[numEdge][edge_numElem[numEdge]] = i*2;
	edge_numElem[numEdge]++;
	elem_edgeIDs[i*2][6] = numEdge;
	edge_elemIDs[numEdge][edge_numElem[numEdge]] = i*2+1;
	edge_numElem[numEdge]++;
	elem_edgeIDs[i*2+1][8] = numEdge;
	ICoordscpy(hex_icoords, &temp_icoords);
	temp_icoords[2]++;
	temp_ID = getIndex(topoDim, hex_dims, temp_icoords);
	if (temp_ID > -1) {
	  edge_elemIDs[numEdge][edge_numElem[numEdge]] = temp_ID*2;
	  edge_numElem[numEdge]++;
	  elem_edgeIDs[temp_ID*2][4] = numEdge;
	  edge_elemIDs[numEdge][edge_numElem[numEdge]] = temp_ID*2+1;
	  edge_numElem[numEdge]++;
	  elem_edgeIDs[temp_ID*2+1][7] = numEdge;
	}
	edge_vertIDs[numEdge][0] = face_vertIDs[numFace][2];
	edge_vertIDs[numEdge][1] = face_vertIDs[numFace][3];
	numEdge++;
      }
      numFace++;
    }
    // Prism#0Face#2 - sometimes
    if (hex_icoords[0]%hex_dims[0] == 0) {
      face_elemIDs[numFace][0] = i*2;
      face_numElem[numFace] = 1;
      elem_faceIDs[i*2][2] = numFace;
      elem_faceOrients[i*2][2] = 1;
      // vertices
      face_numVert[numFace] = 4;
      ICoordscpy(hex_icoords, &temp_icoords);
      temp_icoords[1]++;
      temp_ID = getIndex(topoDim, vert_dims, temp_icoords);
      face_vertIDs[numFace][0] = temp_ID;
      temp_icoords[1]--;
      temp_ID = getIndex(topoDim, vert_dims, temp_icoords);
      face_vertIDs[numFace][1] = temp_ID;
      temp_icoords[2]++;
      temp_ID = getIndex(topoDim, vert_dims, temp_icoords);
      face_vertIDs[numFace][2] = temp_ID;
      temp_icoords[1]++;
      temp_ID = getIndex(topoDim, vert_dims, temp_icoords);
      face_vertIDs[numFace][3] = temp_ID;
      // edges - do 7&8
      // Edge#7 (#0 local to face)
      if (hex_icoords[2]%hex_dims[2] == 0) {
	edge_elemIDs[numEdge][edge_numElem[numEdge]] = i*2;
	edge_numElem[numEdge]++;
	elem_edgeIDs[i*2][7] = numEdge;
	edge_vertIDs[numEdge][0] = face_vertIDs[numFace][0];
	edge_vertIDs[numEdge][1] = face_vertIDs[numFace][1];
	numEdge++;
      }
      // Edge#8 (#2 local to face) (same as 7 in neighbor0 z++)
      {
	edge_elemIDs[numEdge][edge_numElem[numEdge]] = i*2;
	edge_numElem[numEdge]++;
	elem_edgeIDs[i*2][8] = numEdge;
	ICoordscpy(hex_icoords, &temp_icoords);
	temp_icoords[2]++;
	temp_ID = getIndex(topoDim, hex_dims, temp_icoords);
	if (temp_ID > -1) {
	  edge_elemIDs[numEdge][edge_numElem[numEdge]] = temp_ID*2;
	  edge_numElem[numEdge]++;
	  elem_edgeIDs[temp_ID*2][7] = numEdge;
	}
	edge_vertIDs[numEdge][0] = face_vertIDs[numFace][2];
	edge_vertIDs[numEdge][1] = face_vertIDs[numFace][3];
	numEdge++;
      }
      numFace++;
    }
    // Prism#0Face#3 - sometimes
    if (hex_icoords[2]%hex_dims[2] == 0) {
      face_elemIDs[numFace][0] = i*2;
      face_numElem[numFace] = 1;
      elem_faceIDs[i*2][3] = numFace;
      elem_faceOrients[i*2][3] = -1;
      //Vertices
      face_numVert[numFace] = 3;
      ICoordscpy(hex_icoords, &temp_icoords);
      temp_ID = getIndex(topoDim, vert_dims, temp_icoords);
      face_vertIDs[numFace][0] = temp_ID;
      temp_icoords[0]++;
      temp_ID = getIndex(topoDim, vert_dims, temp_icoords);
      face_vertIDs[numFace][1] = temp_ID;
      temp_icoords[0]--;
      temp_icoords[1]++;
      temp_ID = getIndex(topoDim, vert_dims, temp_icoords);
      face_vertIDs[numFace][2] = temp_ID;
      // edges - no edges
      numFace++;
    }
    // Prism#0Face#4 (same as 3 in neighbor0 z++)
    {
      face_elemIDs[numFace][0] = i*2;
      face_numElem[numFace] = 1;
      elem_faceIDs[i*2][4] = numFace;
      elem_faceOrients[i*2][4] = 1;
      ICoordscpy(hex_icoords, &temp_icoords);
      temp_icoords[2]++;
      temp_ID = getIndex(topoDim, hex_dims, temp_icoords);
      if (temp_ID > -1) {
	face_elemIDs[numFace][1] = temp_ID*2;
	face_numElem[numFace]++;
	elem_faceIDs[temp_ID*2][3] = numFace;
	elem_faceOrients[temp_ID*2][3] = -1;
      }
      //Vertices
      face_numVert[numFace] = 3;
      ICoordscpy(hex_icoords, &temp_icoords);
      temp_icoords[2]++;
      temp_ID = getIndex(topoDim, vert_dims, temp_icoords);
      face_vertIDs[numFace][0] = temp_ID;
      temp_icoords[0]++;
      temp_ID = getIndex(topoDim, vert_dims, temp_icoords);
      face_vertIDs[numFace][1] = temp_ID;
      temp_icoords[0]--;
      temp_icoords[1]++;
      temp_ID = getIndex(topoDim, vert_dims, temp_icoords);
      face_vertIDs[numFace][2] = temp_ID;
      // edges - no edges
      numFace++;
    }

    // --- 2nd prism

    // Prism#1Face#0 (same as 2 in neighbor0 x++)
    {
      face_elemIDs[numFace][0] = i*2+1;
      face_numElem[numFace] = 1;
      elem_faceIDs[i*2+1][0] = numFace;
      elem_faceOrients[i*2+1][0] = 1;
      ICoordscpy(hex_icoords, &temp_icoords);
      temp_icoords[0]++;
      temp_ID = getIndex(topoDim, hex_dims, temp_icoords);
      if (temp_ID > -1 ) {
	face_elemIDs[numFace][1] = temp_ID*2;
	face_numElem[numFace]++;
	elem_faceIDs[temp_ID*2][2] = numFace;
	elem_faceOrients[temp_ID*2][2] = -1;
      }     
      // Vertices
      face_numVert[numFace] = 4;
      ICoordscpy(hex_icoords, &temp_icoords);
      temp_icoords[0]++;
      temp_ID = getIndex(topoDim, vert_dims, temp_icoords);
      face_vertIDs[numFace][0] = temp_ID;
      temp_icoords[1]++;
      temp_ID = getIndex(topoDim, vert_dims, temp_icoords);
      face_vertIDs[numFace][1] = temp_ID;
      temp_icoords[2]++;
      temp_ID = getIndex(topoDim, vert_dims, temp_icoords);
      face_vertIDs[numFace][2] = temp_ID;
      temp_icoords[1]--;
      temp_ID = getIndex(topoDim, vert_dims, temp_icoords);
      face_vertIDs[numFace][3] = temp_ID;
      // edges do 0-2 (skip 3) 
      // Edge#0 (0 local to face) (same as 7 neighbor0 x++)
      if (hex_icoords[2]%hex_dims[2] == 0) {
	edge_elemIDs[numEdge][edge_numElem[numEdge]] = i*2+1;
	edge_numElem[numEdge]++;
	elem_edgeIDs[i*2+1][0] = numEdge;
	ICoordscpy(hex_icoords, &temp_icoords);
	temp_icoords[0]++;
	temp_ID = getIndex(topoDim, hex_dims, temp_icoords);
	if (temp_ID > -1) {
	  edge_elemIDs[numEdge][edge_numElem[numEdge]] = temp_ID*2;
	  edge_numElem[numEdge]++;
	  elem_edgeIDs[temp_ID*2][7] = numEdge;
	}
	edge_vertIDs[numEdge][0] = face_vertIDs[numFace][0];
	edge_vertIDs[numEdge][1] = face_vertIDs[numFace][1];
	numEdge++;
      }
      // Edge#1 (1 local to face) (same as 5 neighbor0 & 5 neighbor1 x++,
      //                           same as 1 in neighbor0 & 3 in neighbor1 y++,
      //                           same as 3 neighbor0 x++, y++)
      {
	edge_elemIDs[numEdge][edge_numElem[numEdge]] = i*2+1;
	edge_numElem[numEdge]++;
	elem_edgeIDs[i*2+1][1] = numEdge;
	ICoordscpy(hex_icoords, &temp_icoords);
	temp_icoords[0]++;
	temp_ID = getIndex(topoDim, hex_dims, temp_icoords);
	if (temp_ID > -1) {
	  edge_elemIDs[numEdge][edge_numElem[numEdge]] = temp_ID*2;
	  edge_numElem[numEdge]++;
	  elem_edgeIDs[temp_ID*2][5] = numEdge;
	  edge_elemIDs[numEdge][edge_numElem[numEdge]] = temp_ID*2+1;
	  edge_numElem[numEdge]++;
	  elem_edgeIDs[temp_ID*2+1][5] = numEdge;
	}
	temp_icoords[0]--;
	temp_icoords[1]++;
	temp_ID = getIndex(topoDim, hex_dims, temp_icoords);
	if (temp_ID > -1) {
	  edge_elemIDs[numEdge][edge_numElem[numEdge]] = temp_ID*2;
	  edge_numElem[numEdge]++;
	  elem_edgeIDs[temp_ID*2][1] = numEdge;
	  edge_elemIDs[numEdge][edge_numElem[numEdge]] = temp_ID*2+1;
	  edge_numElem[numEdge]++;
	  elem_edgeIDs[temp_ID*2+1][3] = numEdge;
	}
	temp_icoords[0]++;
	temp_ID = getIndex(topoDim, hex_dims, temp_icoords);
	if (temp_ID > -1) {
	  edge_elemIDs[numEdge][edge_numElem[numEdge]] = temp_ID*2;
	  edge_numElem[numEdge]++;
	  elem_edgeIDs[temp_ID*2][3] = numEdge;
	}
	edge_vertIDs[numEdge][0] = face_vertIDs[numFace][1];
	edge_vertIDs[numEdge][1] = face_vertIDs[numFace][2];
	numEdge++;
      }
      // Edge#2 (2 local to face) (same as 8 neighbor0 x++,
      //                           same as 0 neighbor1 z++
      //                           same as 7 neighbor0 x++, z++)
      {
	edge_elemIDs[numEdge][edge_numElem[numEdge]] = i*2+1;
	edge_numElem[numEdge]++;
	elem_edgeIDs[i*2+1][2] = numEdge;
	ICoordscpy(hex_icoords, &temp_icoords);
	temp_icoords[0]++;
	temp_ID = getIndex(topoDim, hex_dims, temp_icoords);
	if (temp_ID > -1) {
	  edge_elemIDs[numEdge][edge_numElem[numEdge]] = temp_ID*2;
	  edge_numElem[numEdge]++;
	  elem_edgeIDs[temp_ID*2][8] = numEdge;
	}
	temp_icoords[0]--;
	temp_icoords[2]++;
	temp_ID = getIndex(topoDim, hex_dims, temp_icoords);
	if (temp_ID > -1) {
	  edge_elemIDs[numEdge][edge_numElem[numEdge]] = temp_ID*2+1;
	  edge_numElem[numEdge]++;
	  elem_edgeIDs[temp_ID*2+1][0] = numEdge;
	}
	temp_icoords[0]++;
	temp_ID = getIndex(topoDim, hex_dims, temp_icoords);
	if (temp_ID > -1) {
	  edge_elemIDs[numEdge][edge_numElem[numEdge]] = temp_ID*2;
	  edge_numElem[numEdge]++;
	  elem_edgeIDs[temp_ID*2][7] = numEdge;
	}
 	edge_vertIDs[numEdge][0] = face_vertIDs[numFace][2];
	edge_vertIDs[numEdge][1] = face_vertIDs[numFace][3];
	numEdge++;
     }
      numFace++;
    }
    // Prims#1Face#1 - (same as face 0 in neighbor0 y++)
    {
      face_elemIDs[numFace][0] = i*2+1;
      face_numElem[numFace] = 1;
      elem_faceIDs[i*2+1][1] = numFace;
      elem_faceOrients[i*2+1][1] = 1;
      ICoordscpy(hex_icoords, &temp_icoords);
      temp_icoords[1]++;
      temp_ID = getIndex(topoDim, hex_dims, temp_icoords);
      if (temp_ID > -1) {
	face_elemIDs[numFace][1] = temp_ID*2;
	face_numElem[numFace]++;
	elem_faceIDs[temp_ID*2][0] = numFace;
	elem_faceOrients[temp_ID*2][0] = -1;
      }
      // Vertices
      face_numVert[numFace] = 4;
      ICoordscpy(hex_icoords, &temp_icoords);
      temp_icoords[0]++;
      temp_icoords[1]++;
      temp_ID = getIndex(topoDim, vert_dims, temp_icoords);
      face_vertIDs[numFace][0] = temp_ID;
      temp_icoords[0]--;
      temp_ID = getIndex(topoDim, vert_dims, temp_icoords);
      face_vertIDs[numFace][1] = temp_ID;
      temp_icoords[2]++;
      temp_ID = getIndex(topoDim, vert_dims, temp_icoords);
      face_vertIDs[numFace][2] = temp_ID;
      temp_icoords[0]++;
      temp_ID = getIndex(topoDim, vert_dims, temp_icoords);
      face_vertIDs[numFace][3] = temp_ID;
      // edges - do 4 & 6 (skip 5)
      // Edge#4 (0 local to face) (same as 0 in neighbor0 y++)
      if (hex_icoords[2]%hex_dims[2] == 0) {
	edge_elemIDs[numEdge][edge_numElem[numEdge]] = i*2+1;
	edge_numElem[numEdge]++;
	elem_edgeIDs[i*2+1][4] = numEdge;
	ICoordscpy(hex_icoords, &temp_icoords);
	temp_icoords[1]++;
	temp_ID = getIndex(topoDim, hex_dims, temp_icoords);
	if (temp_ID > -1) {
	  edge_elemIDs[numEdge][edge_numElem[numEdge]] = temp_ID*2;
	  edge_numElem[numEdge]++;
	  elem_edgeIDs[temp_ID*2][0] = numEdge;
	}
	edge_vertIDs[numEdge][0] = face_vertIDs[numFace][0];
	edge_vertIDs[numEdge][1] = face_vertIDs[numFace][1];
	numEdge++;
      }
      // Edge#6 (2 local to face) (same as 2 neighbor0 y++,
      //                           same as 4 neighbor1 z++,
      //                           same as 0 neighbor0 y++, z++)
      {
	edge_elemIDs[numEdge][edge_numElem[numEdge]] = i*2+1;
	edge_numElem[numEdge]++;
	elem_edgeIDs[i*2+1][6] = numEdge;
	ICoordscpy(hex_icoords, &temp_icoords);
	temp_icoords[1]++;
	temp_ID = getIndex(topoDim, hex_dims, temp_icoords);
	if (temp_ID > -1) {
	  edge_elemIDs[numEdge][edge_numElem[numEdge]] = temp_ID*2;
	  edge_numElem[numEdge]++;
	  elem_edgeIDs[temp_ID*2][2] = numEdge;
	}
	temp_icoords[1]--;
	temp_icoords[2]++;
	temp_ID = getIndex(topoDim, hex_dims, temp_icoords);
	if (temp_ID > -1) {
	  edge_elemIDs[numEdge][edge_numElem[numEdge]] = temp_ID*2+1;
	  edge_numElem[numEdge]++;
	  elem_edgeIDs[temp_ID*2+1][4] = numEdge;
	}	
	temp_icoords[1]++;
	temp_ID = getIndex(topoDim, hex_dims, temp_icoords);
	if (temp_ID > -1) {
	  edge_elemIDs[numEdge][edge_numElem[numEdge]] = temp_ID*2;
	  edge_numElem[numEdge]++;
	  elem_edgeIDs[temp_ID*2][0] = numEdge;
	}
	edge_vertIDs[numEdge][0] = face_vertIDs[numFace][2];
	edge_vertIDs[numEdge][1] = face_vertIDs[numFace][3];
	numEdge++;
      }
      numFace++;
    }
    // Prism#1Face#2 - same as Prism#0Face#1
    {}
    // Prism#1Face#3 - somtimes
    if (hex_icoords[2]%hex_dims[2] == 0) {
      face_elemIDs[numFace][0] = i*2+1;
      face_numElem[numFace] = 1;
      elem_faceIDs[i*2+1][3] = numFace;
      elem_faceOrients[i*2+1][3] = -1;
      //Vertices
      face_numVert[numFace] = 3;
      ICoordscpy(hex_icoords, &temp_icoords);
      temp_icoords[0]++;
      temp_ID = getIndex(topoDim, vert_dims, temp_icoords);
      face_vertIDs[numFace][0] = temp_ID;
      temp_icoords[1]++;
      temp_ID = getIndex(topoDim, vert_dims, temp_icoords);
      face_vertIDs[numFace][1] = temp_ID;
      temp_icoords[0]--;
      temp_ID = getIndex(topoDim, vert_dims, temp_icoords);
      face_vertIDs[numFace][2] = temp_ID;
      // edges - no edges
      numFace++;
    }
    // Prism#1Face#4 (same as 3 in neighbor1 z++)
    {
      face_elemIDs[numFace][0] = i*2+1;
      face_numElem[numFace] = 1;
      elem_faceIDs[i*2+1][4] = numFace;
      elem_faceOrients[i*2+1][4] = 1;
      ICoordscpy(hex_icoords, &temp_icoords);
      temp_icoords[2]++;
      temp_ID = getIndex(topoDim, hex_dims, temp_icoords);
      if (temp_ID > -1) {
	face_elemIDs[numFace][1] = temp_ID*2+1;
	face_numElem[numFace]++;
	elem_faceIDs[temp_ID*2+1][3] = numFace;
	elem_faceOrients[temp_ID*2+1][3] = -1;
      }
      //Vertices
      face_numVert[numFace] = 3;
      ICoordscpy(hex_icoords, &temp_icoords);
      temp_icoords[0]++;
      temp_icoords[2]++;
      temp_ID = getIndex(topoDim, vert_dims, temp_icoords);
      face_vertIDs[numFace][0] = temp_ID;
      temp_icoords[1]++;
      temp_ID = getIndex(topoDim, vert_dims, temp_icoords);
      face_vertIDs[numFace][1] = temp_ID;
      temp_icoords[0]--;
      temp_ID = getIndex(topoDim, vert_dims, temp_icoords);
      face_vertIDs[numFace][2] = temp_ID;
      // edges - no edges
      numFace++;
    }
  }

  // quick checks that build makes sense
  fail_unless(numFace == numFace_good,
	      "failed to build faces for testing, bug in test code");
  fail_unless(numEdge == numEdge_good, 
	      "failed to build edges for testing, bug in test code");
  for (i = 0; i < numElement_good; i++) {
    for (j = 0; j < numFacePerPrism; j++) 
      fail_unless(elem_faceIDs[i][j] > -1, "abort: bug in test code");
    for (j = 0; j < numEdgePerPrism; j++)
      fail_unless(elem_edgeIDs[i][j] > -1, "abort: bug in test code");
  }
  
  // sort edges' vertices
  for (i = 0; i < numEdge_good; i++) {
    if (edge_vertIDs[i][1] < edge_vertIDs[i][0]) {
      int temp_int = edge_vertIDs[i][0];
      edge_vertIDs[i][0] = edge_vertIDs[i][1];
      edge_vertIDs[i][1] = temp_int;
    }
  }
  // sort & reorder faces' vertices 
  for (i = 0; i < numFace_good; i++) {
    int temp_vertIDs[4];
    // find smallest vertID
    int idx = 0;
    int min = face_vertIDs[i][0];
    for (j = 1; j < face_numVert[i]; j++) {
      if (face_vertIDs[i][j] < min) {
	min = face_vertIDs[i][j];
	idx = j;
      }
    }
    // "rotate" to get smallest vertID first
    for (j = idx; j < face_numVert[i]; j++) 
      temp_vertIDs[j-idx] = face_vertIDs[i][j];
    for (j = 0; j < idx; j++)
      temp_vertIDs[j+(face_numVert[i]-idx)] = face_vertIDs[i][j];
    // swap order to get second smallest vertID next
    if ( (face_numVert[i] == 3 && temp_vertIDs[2] < temp_vertIDs[1]) ||
	 (face_numVert[i] == 4 && temp_vertIDs[3] < temp_vertIDs[1]) ) {
      face_vertIDs[i][0] = temp_vertIDs[0];
      for (j = 1; j < face_numVert[i]; j++) // loop over vertices
	face_vertIDs[i][j] = temp_vertIDs[face_numVert[i]-j];
      for (j = 0; j < face_numElem[i]; j++) {  
	int elemID = face_elemIDs[i][j];
	int relative_faceID = -1;
	for (k = 0; k < numFacePerPrism; k++) 
	  if (elem_faceIDs[elemID][k] == i)
	    relative_faceID = k;
	fail_unless(relative_faceID > -1, "didn't find faceID, bug in test code");
	elem_faceOrients[elemID][relative_faceID] *= -1;
      }
    }
    else 
      for (j = 0; j < face_numVert[i]; j++)
	face_vertIDs[i][j] = temp_vertIDs[j];
  }

  // vertices
  fc_getMeshNumVertex(mesh, &numVertex);
  fail_unless(numVertex == numVertex_good, "mismatch number of vertices");
  rc = _fc_getMeshEdgeParentsOfVerticesPtr(mesh, &numEdgePerVert,
				       &edgeParentsPerVert); 
  fail_unless(rc == FC_SUCCESS, "failed to get vertices' edge parents");
  rc = _fc_getMeshFaceParentsOfVerticesPtr(mesh, &numFacePerVert,
				       &faceParentsPerVert); 
  fail_unless(rc == FC_SUCCESS, "failed to get vertices' face parents");
  rc = _fc_getMeshEdgeParentsOfVerticesPtr(mesh, &numEdgePerVert,
				       &edgeParentsPerVert); 
  fail_unless(rc == FC_SUCCESS, "failed to get vertices' edge parents");
  rc = _fc_getMeshElementParentsOfVerticesPtr(mesh, &numElemPerVert,
				       &elemParentsPerVert); 
  fail_unless(rc == FC_SUCCESS, "failed to get vertices' elem parents");
  rc = _fc_getMeshVertexNeighborsViaEdgePtr(mesh, &numVertNeighsViaEdge,
					    &vertNeighsViaEdge);
  fail_unless(rc == FC_SUCCESS, "failed to get vert neighs via edge");
  for (i = 0; i < numVertex; i++) {
    // calculate stuff
    // setup 
    ICoords vert_icoords, temp_icoords;
    int temp_ID;
    int numHex = 8;
    int actualNumElem;
    int actualNumFace;
    int actualNumEdge;
    int* elemIDs;
    int* faceOwnerIDs;
    int* edgeOwnerIDs;
    FC_SortedIntArray elemOwnerList;
    FC_SortedIntArray faceOwnerList;
    FC_SortedIntArray edgeOwnerList;
    int maxNumEdgeNeigh = 8;
    int actualNumEdgeNeigh;
    int edgeNeighIDs[8];
    // this stencil adjusts owning hex coords (from the vertex coords);
    int hexOwnerStencil[8][3] = 
      { { -1, -1, -1 }, { 0, -1, -1 }, { -1, 0, -1 }, { 0, 0, -1 },
	{ -1, -1,  0 }, { 0, -1,  0 }, { -1, 0,  0 }, { 0, 0,  0 } };
    // whether prism j from quad i has the vertex
    int prismHasVert1[8] = { 0, 1, 1, 1, 0, 1, 1, 1 };
    int prismHasVert2[8] = { 1, 1, 1, 0, 1, 1, 1, 0 };
    int faceLU1[8][3] = 
      { { -1, -1, -1 }, { 1, 2, 4 }, { 0, 1, 4 }, { 0, 2, 4 },
	{ -1, -1, -1 }, { 1, 2, 3 }, { 0, 1, 3 }, { 0, 2, 3 } } ;
    int faceLU2[8][3] = 
      { { 0, 1, 4 }, { 1, 2, 4 }, { 0, 2, 4 }, { -1, -1, -1 },
	{ 0, 1, 3 }, { 1, 2, 3 }, { 0, 2, 3 }, { -1, -1, -1 } } ;
    int edgeLU1[8][3] = 
      { { -1, -1, -1 }, { 5, 6, 8 }, { 1, 2, 6 }, { 2, 3, 8 },
	{ -1, -1, -1 }, { 4, 5, 7 }, { 0, 1, 4 }, { 0, 3, 7 } };
    int edgeLU2[8][3] = 
      { { 1, 2, 6 }, { 5, 6, 8 }, { 2, 3, 8 }, { -1, -1, -1 },
	{ 0, 1, 4 }, { 4, 5, 7 }, { 0, 3, 7 }, { -1, -1, -1 } };
    // this stencil adjusts vertex ijk coords
    int edgeNeighStencil[8][3] = 
      { { 0, 0, -1 },
	{ 0, -1, 0 }, { 1, -1, 0 }, { -1, 0, 0 }, 
	{ 1, 0, 0 }, { -1, 1, 0 }, { 0, 1, 0 },
	{ 0, 0, 1 } };
 
    // setup II - shared calcs
    getICoords(i, topoDim, vert_dims, &vert_icoords);
    // calc owners
    fc_initSortedIntArray(&elemOwnerList);
    fc_initSortedIntArray(&faceOwnerList);
    fc_initSortedIntArray(&edgeOwnerList);
    for (j = 0; j < numHex; j++) {
      for (k = 0; k < 3; k++)
	temp_icoords[k] = vert_icoords[k] + hexOwnerStencil[j][k];
      temp_ID = getIndex(topoDim, hex_dims, temp_icoords);
      if (temp_ID > -1) {
	if (prismHasVert1[j])
	  fc_addIntToSortedIntArray(&elemOwnerList, temp_ID*2);
	if (prismHasVert2[j])
	  fc_addIntToSortedIntArray(&elemOwnerList, temp_ID*2+1);
	for (k = 0; k < 3; k++) {
	  if (faceLU1[j][k] > -1) 
	    fc_addIntToSortedIntArray(&faceOwnerList,
				      elem_faceIDs[temp_ID*2][faceLU1[j][k]]);
	  if (faceLU2[j][k] > -1)
	    fc_addIntToSortedIntArray(&faceOwnerList,
				      elem_faceIDs[temp_ID*2+1][faceLU2[j][k]]);
	  if (edgeLU1[j][k] > -1)
	    fc_addIntToSortedIntArray(&edgeOwnerList,
				      elem_edgeIDs[temp_ID*2][edgeLU1[j][k]]);
	  if (edgeLU2[j][k] > -1)
	    fc_addIntToSortedIntArray(&edgeOwnerList,
				      elem_edgeIDs[temp_ID*2+1][edgeLU2[j][k]]);
	}
      }
    }
    fc_convertSortedIntArrayToIntArray(&elemOwnerList, &actualNumElem, &elemIDs);
    fc_convertSortedIntArrayToIntArray(&faceOwnerList, &actualNumFace, &faceOwnerIDs);
    fc_convertSortedIntArrayToIntArray(&edgeOwnerList, &actualNumEdge, &edgeOwnerIDs);
    // calc edge neighbors
    actualNumEdgeNeigh = 0;
    for (j = 0; j < maxNumEdgeNeigh; j++) {
      for (k = 0; k < 3; k++)
	temp_icoords[k] = vert_icoords[k] + edgeNeighStencil[j][k];
      temp_ID = getIndex(topoDim, vert_dims, temp_icoords);
      if (temp_ID > -1) {
	edgeNeighIDs[actualNumEdgeNeigh] = temp_ID;
	actualNumEdgeNeigh++;
      }
    }

    // tests
    fail_unless(numEdgePerVert[i] == actualNumEdge,
		"mismatch in numEdgePerVert");
    fail_unless(!memcmp(edgeParentsPerVert[i], edgeOwnerIDs,
			numEdgePerVert[i]*sizeof(int)), 
		"mismatch in edgeParentsPerVert");
    free(edgeOwnerIDs);
    fail_unless(numFacePerVert[i] == actualNumFace, 
		"mismatch in numFacePerVert");
    fail_unless(!memcmp(faceParentsPerVert[i], faceOwnerIDs,
			numFacePerVert[i]*sizeof(int)), 
		"mismatch in faceParentsPerVert");
    free(faceOwnerIDs);
    fail_unless(numElemPerVert[i] == actualNumElem, 
		"vertex: mismatch in number of elements");
    fail_unless(!memcmp(elemParentsPerVert[i], elemIDs, 
			numElemPerVert[i]*sizeof(int)),
		"vertex: mismatch in elementIDs");
    free(elemIDs);
    fail_unless(numVertNeighsViaEdge[i] == actualNumEdgeNeigh,
		"mismatch of numVertNeighsViaEdge");
    fail_unless(!memcmp(vertNeighsViaEdge[i], edgeNeighIDs,
			numVertNeighsViaEdge[i]*sizeof(int)),
		"mismatch of vertNeighsViaEdge");
  }

  // edges
  fc_getMeshNumEdge(mesh, &numEdge);
  fail_unless(numEdge == numEdge_good, "mismatch in number edges");
  rc = fc_getMeshEdgeConnsPtr(mesh, &edgeConns);
  fail_unless(rc == FC_SUCCESS, "failed to get edge conns");
  rc = _fc_getMeshElementParentsOfEdgesPtr(mesh, &numElemPerEdge,
					&elemParentsPerEdge);
  fail_unless(rc == FC_SUCCESS, "failed to get edge's elem parents");
  for (i = 0; i < numEdge; i++) {
    // tests
    fail_unless(!memcmp(&edgeConns[i*2], edge_vertIDs[i], 2*sizeof(int)),
		"edge: mismatch of vertIDs");
    fail_unless(numElemPerEdge[i] == edge_numElem[i], 
		"mismatch in numElemPerEdge");
    fail_unless(!memcmp(elemParentsPerEdge[i], edge_elemIDs[i],
			numElemPerEdge[i]*sizeof(int)), 
		"mismatch of elemParentsPerEdge");
  }

  // faces
  fc_getMeshNumFace(mesh, &numFace);
  fail_unless(numFace == numFace_good, "mismatch in number faces");
  rc = fc_getMeshFaceTypesPtr(mesh, &faceTypes);
  fail_unless(rc == FC_SUCCESS, "failed to get face types");
  rc = fc_getMeshFaceConnsPtr(mesh, &numVertPerFace, &maxNumVertPerFace,
			      &faceConns);
  fail_unless(rc == FC_SUCCESS, "failed to get face conns");
  rc = _fc_getMeshElementParentsOfFacesPtr(mesh, &numElemPerFace,
					&elemParentsPerFace);
  fail_unless(rc == FC_SUCCESS, "failed to get face's elem parents");
  for (i = 0; i < numFace; i++) {
    FC_ElementType temp_type = FC_ET_QUAD;
    if (face_numVert[i] == 3)
      temp_type = FC_ET_TRI;
    // tests
    fail_unless(faceTypes[i] == temp_type, "face: mismatch in face type");
    fail_unless(numVertPerFace[i] == face_numVert[i], 
		"face: mismatch in number of vertices");
    fail_unless(!memcmp(&faceConns[i*maxNumVertPerFace], face_vertIDs[i],
			face_numVert[i]*sizeof(int)),
		  "face: mismatch in vertex IDs");
    for (j = face_numVert[i]; j < maxNumVertPerFace; j++)
      fail_unless(faceConns[i*maxNumVertPerFace+j] == -1,
		  "face: mismatch in vertexIDs");
    fail_unless(numElemPerFace[i] == face_numElem[i], 
		"mismatch in numElemPerFace");
    fail_unless(!memcmp(elemParentsPerFace[i], face_elemIDs[i],
			numElemPerFace[i]*sizeof(int)),
		"mismatch in elemParentsPerFace");
  }

  // elements
  fc_getMeshNumElement(mesh, &numElement);
  fail_unless(numElement == numElement_good, "mismatch in number of elements");
  rc = fc_getMeshElementToEdgeConnsPtr(mesh, &elemToEdgeConns);
  fail_unless(rc == FC_SUCCESS, "failed to get elem to edge conns");
  rc = fc_getMeshElementToFaceConnsPtr(mesh, &elemToFaceConns);
  fail_unless(rc == FC_SUCCESS, "failed to get elem to edge conns");
  rc = fc_getMeshElementFaceOrientsPtr(mesh, &elemFaceOrients);
  fail_unless(rc == FC_SUCCESS, "failed to get elem to edge conns");
  rc = _fc_getMeshElementNeighborsViaVertexPtr(mesh, &numElemNeighsViaVert,
					       &elemNeighsViaVert);
  fail_unless(rc == FC_SUCCESS, "failed to get elemNeighsViaVert");
  rc = _fc_getMeshElementNeighborsViaEdgePtr(mesh, &numElemNeighsViaEdge,
					     &elemNeighsViaEdge);
  fail_unless(rc == FC_SUCCESS, "failed to get elemNeighsViaEdge");
  rc = _fc_getMeshElementNeighborsViaFacePtr(mesh, &numElemNeighsViaFace,
					     &elemNeighsViaFace);
  fail_unless(rc == FC_SUCCESS, "failed to get elemNeighsViaFace");
  for (i = 0; i < numElement; i++) {
    // calculate stuff
    // setup
    int hex_index;
    ICoords hex_icoords, temp_icoords;
    int numVert = 6;
    int vertIDs[6];
    // hex neighbors
    int maxNumVertNeigh = 26;  // 6 + 12 + 8
    int maxNumEdgeNeigh = 18;  // 6 + 12
    int maxNumFaceNeigh = 6;
    int actualNumVertNeigh;
    int actualNumEdgeNeigh;
    int actualNumFaceNeigh;
    int vertNeighIDs[38];    // 
    int edgeNeighIDs[20];    // 
    int faceNeighIDs[5];     // prism has 5
    // this stencil adjusts ijk coords for hex
    int vertNeighStencil[26][3] = {
      // bottom 9
      { -1, -1, -1 }, {  0, -1, -1 }, {  1, -1, -1 },
      { -1,  0, -1 }, {  0,  0, -1 }, {  1,  0, -1 },
      { -1,  1, -1 }, {  0,  1, -1 }, {  1,  1, -1 },
      // middle 8
      { -1, -1,  0 }, {  0, -1,  0 }, {  1, -1,  0 },
      { -1,  0,  0 },                 {  1,  0,  0 },
      { -1,  1 , 0 }, {  0,  1,  0 }, {  1,  1,  0 },
      // top 9
      { -1, -1,  1 }, {  0, -1,  1 }, {  1, -1,  1 },
      { -1,  0,  1 }, {  0,  0,  1 }, {  1,  0,  1 },
      { -1,  1,  1 }, {  0,  1,  1 }, {  1,  1,  1 } };
    int vertLU1[26] = {
      // bottom 9
      1, 2, 2, 2, 2, 0, 2, 0, -1,
      // middle 8
      1, 2, 2, 2,    0, 2, 0, -1,
      // top 9
      1, 2, 2, 2, 2, 0, 2, 0, -1 };
    int vertLU2[26] = {
      // bottom 9
      -1, 1, 2, 1, 2, 2, 2, 2, 0,
      // middle 8
      -1, 1, 2, 1,    2, 2, 2, 0,
      // top 9
      -1, 1, 2, 1, 2, 2, 2, 2, 0 };
    int edgeNeighStencil[18][3] = { 
      // bottom 5
      {  0, -1, -1 },
      { -1,  0, -1 }, {  0,  0, -1 }, {  1,  0, -1 },
      {  0,  1, -1 },
      // middle 8
      { -1, -1,  0 }, {  0, -1,  0 }, {  1, -1,  0 },
      { -1,  0,  0 },                 {  1,  0,  0 },
      { -1,  1 , 0 }, {  0,  1,  0 }, {  1,  1,  0 },
      // top 5
      {  0, -1,  1 },
      { -1,  0,  1 }, {  0,  0,  1 }, {  1,  0,  1 },
      {  0,  1,  1 } };
    int edgeLU1[18] = { 
      // bottom 5
      1, 1, 2, -1, -1,
      // middle 8
      1, 2, 2, 2, 0, 2, 0, -1, 
      // top 5
      1, 1, 2, -1, -1 };
    int edgeLU2[18] = { 
      // bottom 5
      -1, -1, 2, 0, 0,
      // middle 8
      -1, 1, 2, 1, 2, 2, 2, 0, 
      // top 5
      -1, -1, 2, 0, 0 };
    int faceNeighStencil[6][3] = {
      // bottom 1 
      {  0,  0, -1 },
      // middle 4 
      {  0, -1,  0 }, { -1,  0,  0 },
      {  1,  0,  0 }, {  0,  1,  0 }, 
      // top 1
      {  0,  0,  1 } };
    int faceLU1[6] = { 0,   1, 1, -1, -1,    0 };
    int faceLU2[6] = { 1,   -1, -1, 0, 0,    1 };
    // setup II - shared calcs      
    getICoords(i/2, topoDim, hex_dims, &hex_icoords);
    // owned vertices
    for (j = 0; j < numVert; j++) { // walk around prism
      // prism0
      if (i%2 == 0) {
	switch (j) {
	case 0: // start of bottom face
	  ICoordscpy(hex_icoords, &temp_icoords);
          break;
	case 3: // start of top face
	  ICoordscpy(hex_icoords, &temp_icoords);
	  temp_icoords[2]++;
	  break;
	case 1: case 4: 
	  temp_icoords[0]++;  
	  break;
	case 2: case 5:
	  temp_icoords[0]--;
	  temp_icoords[1]++;  
	  break;
	}
      }
      // prism1
      else {
      	switch (j) {
	case 0: // start of bottom face
	  ICoordscpy(hex_icoords, &temp_icoords);
	  temp_icoords[0]++;
          break;
	case 3: // start of top face
	  ICoordscpy(hex_icoords, &temp_icoords);
	  temp_icoords[0]++;
	  temp_icoords[2]++;
	  break;
	case 1: case 4: 
	  temp_icoords[1]++;  
	  break;
	case 2: case 5:
	  temp_icoords[0]--;  
	  break;
	}
      }
      vertIDs[j] = getIndex(topoDim, vert_dims, temp_icoords);
    }
    // vertex neighbors
    actualNumVertNeigh = 0;
    for (j = 0; j < maxNumVertNeigh; j++) {
      if ( j == 13) { // insert the neighbor sitting in the same quad
	if (i%2 == 0) {
	  vertNeighIDs[actualNumVertNeigh] = i+1;
	  actualNumVertNeigh++;
	}
	else {
	  vertNeighIDs[actualNumVertNeigh] = i-1;
	  actualNumVertNeigh++;
	}
      }
      for (k = 0; k < 3; k++)
	temp_icoords[k] = hex_icoords[k] + vertNeighStencil[j][k];
      hex_index = getIndex(topoDim, hex_dims, temp_icoords);
      if (hex_index > -1) {
	if (i%2 == 0 && vertLU1[j] > -1) {
	  if (vertLU1[j] == 0 || vertLU1[j] == 2) {
	    vertNeighIDs[actualNumVertNeigh] = hex_index*2;
	    actualNumVertNeigh++;
	  }
	  if (vertLU1[j] == 1 || vertLU1[j] == 2) {
	    vertNeighIDs[actualNumVertNeigh] = hex_index*2+1;
	    actualNumVertNeigh++;
	  }
	}
	else if (i%2 == 1 && vertLU2[j] > -1) {
	  if (vertLU2[j] == 0 || vertLU2[j] == 2) {
	    vertNeighIDs[actualNumVertNeigh] = hex_index*2;
	    actualNumVertNeigh++;
	  }
	  if (vertLU2[j] == 1 || vertLU2[j] == 2) {
	    vertNeighIDs[actualNumVertNeigh] = hex_index*2+1;
	    actualNumVertNeigh++;
	  }
	}
      }
    }
    // edge neighbors
    actualNumEdgeNeigh = 0;
    for (j = 0; j < maxNumEdgeNeigh; j++) {
      if (j == 9) { // insert the neighbor sitting in the same quad
	if (i%2 == 0) {
	  edgeNeighIDs[actualNumEdgeNeigh] = i+1;
	  actualNumEdgeNeigh++;
	}
	else {
	  edgeNeighIDs[actualNumEdgeNeigh] = i-1;
	  actualNumEdgeNeigh++;
	}
      }
      for (k = 0; k < 3; k++)
	temp_icoords[k] = hex_icoords[k] + edgeNeighStencil[j][k];
      hex_index = getIndex(topoDim, hex_dims, temp_icoords);
      if (hex_index > -1) {
	if (i%2 == 0 && edgeLU1[j] > -1) {
	  if (edgeLU1[j] == 0 || edgeLU1[j] == 2) {
	    edgeNeighIDs[actualNumEdgeNeigh] = hex_index*2;
	    actualNumEdgeNeigh++;
	  }
	  if (edgeLU1[j] == 1 || edgeLU1[j] == 2) {
	    edgeNeighIDs[actualNumEdgeNeigh] = hex_index*2+1;
	    actualNumEdgeNeigh++;
	  }
	}
	else if (i%2 == 1 && edgeLU2[j] > -1) {
	  if (edgeLU2[j] == 0 || edgeLU2[j] == 2) {
	    edgeNeighIDs[actualNumEdgeNeigh] = hex_index*2;
	    actualNumEdgeNeigh++;
	  }
	  if (edgeLU2[j] == 1 || edgeLU2[j] == 2) {
	    edgeNeighIDs[actualNumEdgeNeigh] = hex_index*2+1;
	    actualNumEdgeNeigh++;
	  }
	}
      }
    }
    // face neighbors
    actualNumFaceNeigh = 0;
    for (j = 0; j < maxNumFaceNeigh; j++) {
      if (j == 3) { // insert the neighbor sitting in same quad
	if (i%2 == 0) {
	  faceNeighIDs[actualNumFaceNeigh] = i+1;
	  actualNumFaceNeigh++;
	}
	else {
	  faceNeighIDs[actualNumFaceNeigh] = i-1;
	  actualNumFaceNeigh++;
	}
      }
      for (k = 0; k < 3; k++)
	temp_icoords[k] = hex_icoords[k] + faceNeighStencil[j][k];
      hex_index = getIndex(topoDim, hex_dims, temp_icoords);
      if (hex_index > -1) {
	if (i%2 == 0 && faceLU1[j] > -1) {
	  faceNeighIDs[actualNumFaceNeigh] = hex_index*2+faceLU1[j];
	  actualNumFaceNeigh++;
	}
	else if (i%2 == 1 && faceLU2[j] > -1) {
	  faceNeighIDs[actualNumFaceNeigh] = hex_index*2+faceLU2[j];
	  actualNumFaceNeigh++;
	}
      }
    }

    // element tests
    fail_unless(!memcmp(&elemToEdgeConns[i*numEdgePerPrism], elem_edgeIDs[i],
			numEdgePerPrism*sizeof(int)),
		"element: mismatch in edgeIDs");
    fail_unless(!memcmp(&elemToFaceConns[i*numFacePerPrism], elem_faceIDs[i],
			numFacePerPrism*sizeof(int)),
		"element: mismatch in face IDs");
    fail_unless(!memcmp(&elemFaceOrients[i*numFacePerPrism],
			elem_faceOrients[i], numFacePerPrism*sizeof(int)),
		"element: mismatch in face orients");
    fail_unless(numElemNeighsViaVert[i] == actualNumVertNeigh,
		"mismatch of numElemNeighsViaVert");
    fail_unless(!memcmp(elemNeighsViaVert[i], vertNeighIDs,
			actualNumVertNeigh*sizeof(int)), 
		"element: mismatch in vertex neighbors");
    fail_unless(numElemNeighsViaEdge[i] == actualNumEdgeNeigh,
		"mismatch of numElemNeighsViaEdge");
    fail_unless(!memcmp(elemNeighsViaEdge[i], edgeNeighIDs,
			actualNumEdgeNeigh*sizeof(int)),  
		  "mismatch of elemNeighsViaEdge");
    fail_unless(numElemNeighsViaFace[i] == actualNumFaceNeigh, 
		"mismatch of numElemNeighsViaFace");
    fail_unless(!memcmp(elemNeighsViaFace[i], faceNeighIDs,
			actualNumFaceNeigh*sizeof(int)), 
		  "mismatch of elemNeighsViaFace");
  }
  fail_unless(i == numElement, "test error: didn't finish loop!");
  
  fc_deleteDataset(dataset);
}
END_TEST
  

// Populate the Suite with the tests

Suite *mesh_suite(void)
{
  Suite *suite = suite_create("Mesh");

  TCase *tc_private_mesh = tcase_create(" - Private Mesh ");
  TCase *tc_mesh = tcase_create(" - Mesh Interface ");
  TCase *tc_build_structs = tcase_create(" - Build Conns/Neighs");

  // private sequence
  suite_add_tcase(suite, tc_private_mesh);
  tcase_add_checked_fixture(tc_private_mesh, mesh_setup, mesh_teardown);
  tcase_add_test(tc_private_mesh, slot_new_delete);

  // mesh interface
  suite_add_tcase(suite, tc_mesh);
  tcase_add_checked_fixture(tc_mesh, mesh_setup, mesh_teardown);
  tcase_add_test(tc_mesh, create_get_delete);
  tcase_add_test(tc_mesh, metadata_query);
  tcase_add_test(tc_mesh, coords_conns_query);
  tcase_add_test(tc_mesh, copy_test);
  tcase_add_test(tc_mesh, release_mesh);
  tcase_add_test(tc_mesh, coords_as_var);
  tcase_add_test(tc_mesh, create_subset_mesh);
  tcase_add_test(tc_mesh, simple_hex_mesh);
  tcase_add_test(tc_mesh, sphere_mesh);
  tcase_add_test(tc_mesh, print);

  // connectivities and neighbors
  suite_add_tcase(suite, tc_build_structs);
  tcase_add_checked_fixture(tc_build_structs, mesh_setup, mesh_teardown);
  tcase_add_test(tc_build_structs, build_point_mesh);
  tcase_add_test(tc_build_structs, build_line_mesh);
  tcase_add_test(tc_build_structs, build_quad_mesh);
  tcase_add_test(tc_build_structs, build_prism_mesh);
  tcase_add_test(tc_build_structs, build_hex_mesh);

  return suite;
}
