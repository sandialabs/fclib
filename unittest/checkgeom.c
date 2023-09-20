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
 * \file checkgeom.c
 * \brief Unit testing of \ref GeometricRelations module.
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/unittest/checkgeom.c,v $
 * $Revision: 1.115 $ 
 * $Date: 2007/07/17 19:25:39 $
 */

#include <stdlib.h>
#include <string.h>
#include <check.h>
#include <math.h>
#include "fc.h"
#include "fcP.h"
#include "checkall.h"

// **** global data

static char exo_file[25] = "../data/gen_tri.ex2";
static char mesh_name[25] = "tri mesh";

// *** helper routines

// Helper routine to calculate arbirary rotation of coordinates
static void rotate2D(double* x_p, double* y_p) {
  double x = *x_p, y = *y_p;
  double new_x, new_y;
  
  // FIX calculate these from degrees using PI = d*pi/180
  // arbitrary rotations about each axis:
  double THETA = 0.366519143; // 21 degrees

  new_x = x*cos(THETA) - y*sin(THETA);
  new_y = x*sin(THETA) + y*cos(THETA);

  *x_p = new_x;
  *y_p = new_y;
}

// Helper routine to calculate arbirary rotation of coordinates
static void rotate3D(double* x_p, double* y_p, double* z_p) {
  double x = *x_p, y = *y_p, z = *z_p;
  double new_x, new_y, new_z;
  
  // FIX calculate these from degrees using PI = d*pi/180
  // arbitrary rotations about each axis:
  double THETA_X = 0.366519143; // 21 degrees
  double THETA_Y = 0.296705973; // 17 degrees
  double THETA_Z = 0.191986218; // 11 degrees
  double sinx = sin(THETA_X), siny = sin(THETA_Y), sinz = sin(THETA_Z);
  double cosx = cos(THETA_X), cosy = cos(THETA_Y), cosz = cos(THETA_Z);

  new_x = x*cosy*cosz - y*cosy*sinz + z*siny;
  new_y = x*(sinx*siny*cosz + cosx*sinz) + y*(-sinx*siny*sinz + cosx*cosz) -
    z*sinx*cosy;
  new_z = x*(-cosx*siny*cosz + sinx*sinz) + y*(cosx*siny*sinz + sinx*cosz) +
    z*cosx*cosy;

  *x_p = new_x;
  *y_p = new_y;
  *z_p = new_z;
}

// **** test fixtures

static void geom_setup(void) {
  FC_ReturnCode rc;
  if (isForking) {
    fc_setLibraryVerbosity(fc_messages);
    rc = fc_initLibrary();
    fail_unless(rc == FC_SUCCESS, "test aborted: failed to init library");
  }
}

static void geom_teardown(void) {
  FC_ReturnCode rc;
  if (isForking) {
    rc = fc_finalLibrary();
    fail_unless(rc == FC_SUCCESS, "failed to final library at end of test");
  }
}

// test calculation of euclidean distance
START_TEST(euclid_distances)
{
  FC_ReturnCode rc;
  int i;
  FC_Coords coords1 = { 1., 2., 3. };
  FC_Coords coords2 = { 2.25, 7.101, 99.0001 };
  double distances[3], temp_distance;
  double dist2s[3], temp_dist2;

  // calculate answers
  distances[0] = coords2[0] - coords1[0];
  dist2s[0] = distances[0]*distances[0];
  for (i = 1; i < 3; i++) {
    dist2s[i] = distances[i-1]*distances[i-1] + 
                     (coords2[i]-coords1[i])*(coords2[i]-coords1[i]);
    distances[i] = sqrt(dist2s[i]);
  }

  // test good cases
  for (i = 0; i < 3; i++) {
    rc = fc_calcEuclideanDistance(coords1, coords2, i+1, &temp_distance);
    fail_unless(rc == FC_SUCCESS, "failed to calc");
    fail_unless(temp_distance == distances[i], "mismatch of distance");
    rc = fc_calcEuclideanDistance(coords2, coords1, i+1, &temp_distance);
    fail_unless(rc == FC_SUCCESS, "failed to calc");
    fail_unless(temp_distance == distances[i], "mismatch of distance");
    rc = fc_calcSquaredEuclideanDistance(coords1, coords2, i+1, &temp_dist2);
    fail_unless(rc == FC_SUCCESS, "failed to calc");
    fail_unless(temp_dist2 == dist2s[i], "mismatch of squared distance");
    rc = fc_calcSquaredEuclideanDistance(coords2, coords1, i+1, &temp_dist2);
    fail_unless(rc == FC_SUCCESS, "failed to calc");
    fail_unless(temp_dist2 == dist2s[i], "mismatch of squared distance");
  }
  
  // test bad inputs to fc_calcEuclideanDistance
  rc = fc_calcEuclideanDistance(NULL, coords2, 3, &temp_distance);
  fail_unless(rc != FC_SUCCESS, "should fail with bad coord1");
  fail_unless(temp_distance == -1, "fail should return -1");
  rc = fc_calcEuclideanDistance(coords1, NULL, 3, &temp_distance);
  fail_unless(rc != FC_SUCCESS, "should fail with bad coord2");
  fail_unless(temp_distance == -1, "fail should return -1");
  rc = fc_calcEuclideanDistance(coords1, coords2, 0, &temp_distance);
  fail_unless(rc != FC_SUCCESS, "should fail if dim < 1");
  fail_unless(temp_distance == -1, "fail should return -1");
  rc = fc_calcEuclideanDistance(coords1, coords2, 4, &temp_distance);
  fail_unless(rc != FC_SUCCESS, "should fail with dim > 3");
  fail_unless(temp_distance == -1, "fail should return -1");
  rc = fc_calcEuclideanDistance(coords1, coords2, 3, NULL);
  fail_unless(rc != FC_SUCCESS, "should fail if no distance");

  // test bad inputs to fc_calcEuclideanDistance
  rc = fc_calcSquaredEuclideanDistance(NULL, coords2, 3, &temp_dist2);
  fail_unless(rc != FC_SUCCESS, "should fail with bad coord1");
  fail_unless(temp_dist2 == -1, "fail should return -1");
  rc = fc_calcSquaredEuclideanDistance(coords1, NULL, 3, &temp_dist2);
  fail_unless(rc != FC_SUCCESS, "should fail with bad coord2");
  fail_unless(temp_dist2 == -1, "fail should return -1");
  rc = fc_calcSquaredEuclideanDistance(coords1, coords2, 0, &temp_dist2);
  fail_unless(rc != FC_SUCCESS, "should fail if dim < 1");
  fail_unless(temp_dist2 == -1, "fail should return -1");
  rc = fc_calcSquaredEuclideanDistance(coords1, coords2, 4, &temp_dist2);
  fail_unless(rc != FC_SUCCESS, "should fail with dim > 3");
  fail_unless(temp_dist2 == -1, "fail should return -1");
  rc = fc_calcSquaredEuclideanDistance(coords1, coords2, 3, NULL);
  fail_unless(rc != FC_SUCCESS, "should fail if no distance");
}
END_TEST


//test calculation of angle between vectors
START_TEST(angle_between)
{
  FC_ReturnCode rc;

  FC_Vector v1;
  FC_Vector v2;
  double angle;

  int i;

  //test cases
  for (i = 0; i < 3; i++){
    v1[i] = 0;
    v2[i] = 0;
  }
  v1[2] = 1;
  v2[0] = 1;
  rc = fc_calcAngleBetweenVectors(v1,v2, &angle);
  fail_unless(rc == FC_SUCCESS, "Should work for 90 degree angle");
  fail_unless(FC_DBL_EQUIV(angle,90), "Should work for 90 degree angle");

  for (i = 0; i < 3; i++){
    v1[i] = 0;
    v2[i] = 0;
  }
  v1[2] = 1;
  v2[0] = -1;
  rc = fc_calcAngleBetweenVectors(v1,v2, &angle);
  fail_unless(rc == FC_SUCCESS, "Should work for 90 degree angle");
  fail_unless(FC_DBL_EQUIV(angle,90), "Wrong return val");

  for (i = 0; i < 3; i++){
    v1[i] = 0;
    v2[i] = 0;
  }
  v1[2] = -1;
  v2[0] = -1;
  rc = fc_calcAngleBetweenVectors(v1,v2, &angle);
  fail_unless(rc == FC_SUCCESS, "Should work for 90 degree angle");
  fail_unless(FC_DBL_EQUIV(angle,90), "Wrong return val");

  for (i = 0; i < 3; i++){
    v1[i] = 0;
    v2[i] = 0;
  }
  v1[2] = 1;
  v2[2] = 1;
  rc = fc_calcAngleBetweenVectors(v1,v2, &angle);
  fail_unless(rc == FC_SUCCESS, "Should work for zero degree angle");
  fail_unless(FC_DBL_EQUIV(angle,0), "Wrong return val"); 

  for (i = 0; i < 3; i++){
    v1[i] = 0;
    v2[i] = 0;
  }
  v1[2] = -1;
  v2[2] = -1;
  rc = fc_calcAngleBetweenVectors(v1,v2, &angle);
  fail_unless(rc == FC_SUCCESS, "Should work for zero degree angle");
  fail_unless(FC_DBL_EQUIV(angle,0), "Wrong return val"); 

  for (i = 0; i < 3; i++){
    v1[i] = 0;
    v2[i] = 0;
  }
  v1[2] = -1;
  v2[2] = 1;
  rc = fc_calcAngleBetweenVectors(v1,v2, &angle);
  fail_unless(rc == FC_SUCCESS, "Should work for 180 degree angle");
  fail_unless(FC_DBL_EQUIV(angle,180), "Wrong return val"); 

  v1[0] = 1;
  v1[1] = 0;
  v1[2] = 0;
  v2[0] = 1;
  v2[1] = 1;
  v2[2] = 0;
  rc = fc_calcAngleBetweenVectors(v1,v2, &angle);
  fail_unless(rc == FC_SUCCESS, "Should work off axis, non normalized vector");
  fail_unless(FC_DBL_EQUIV(angle,45), "Wrong return val"); 

  v1[0] = 1;
  v1[1] = 0;
  v1[2] = 0;
  v2[0] = -1;
  v2[1] = -1;
  v2[2] = 0;
  rc = fc_calcAngleBetweenVectors(v1,v2, &angle);
  fail_unless(rc == FC_SUCCESS, "Should work off axis, non normalized vector");
  fail_unless(FC_DBL_EQUIV(angle,135), "Wrong return val"); 

  v1[0] = 1;
  v1[1] = 0;
  v1[2] = 0;
  v2[0] = 1;
  v2[1] = -1;
  v2[2] = 0;
  rc = fc_calcAngleBetweenVectors(v1,v2, &angle);
  fail_unless(rc == FC_SUCCESS, "Should work off axis, non normalized vector");
  fail_unless(FC_DBL_EQUIV(angle,45), "Wrong return val"); 

  v1[0] = 1;
  v1[1] = 0;
  v1[2] = 0;
  v2[0] = -1;
  v2[1] = 1;
  v2[2] = 0;
  rc = fc_calcAngleBetweenVectors(v1,v2, &angle);
  fail_unless(rc == FC_SUCCESS, "Should work off axis, non normalized vector");
  fail_unless(FC_DBL_EQUIV(angle,135), "Wrong return val"); 


  //bad cases
  rc = fc_calcAngleBetweenVectors(NULL,v2, &angle);
  fail_unless(rc != FC_SUCCESS, "Should fail for NULL vector");
  fail_unless(angle == -1, "should return -1 angle when fail");

  rc = fc_calcAngleBetweenVectors(v1,NULL, &angle);
  fail_unless(rc != FC_SUCCESS, "Should fail for NULL vector");
  fail_unless(angle == -1, "should return -1 angle when fail");

  rc = fc_calcAngleBetweenVectors(v1,v2,NULL);
  fail_unless(rc != FC_SUCCESS, "Should fail for NULL return param");


  v1[0] = 0;
  v1[1] = 0;
  v1[2] = 0;
  v2[0] = -1;
  v2[1] = 1;
  v2[2] = 0;
  rc = fc_calcAngleBetweenVectors(v1,v2, &angle);
  fail_unless(rc != FC_SUCCESS, "Should fail for zero length vector");
  fail_unless(angle == -1, "should return -1 angle when fail");

  rc = fc_calcAngleBetweenVectors(v2,v1, &angle);
  fail_unless(rc != FC_SUCCESS, "Should fail for zero length vector");
  fail_unless(angle == -1, "should return -1 angle when fail");
}END_TEST

// test the valid displ var tester
START_TEST(valid_displ)
{
  FC_ReturnCode rc;
  int i;
  FC_Dataset dataset;
  FC_Mesh mesh, otherMesh, badMesh = { 999, 999 };
  FC_Variable var, otherVar, badVar = {999, 999 };
  int numVert = 5;
  int conns[5] = { 0, 1, 2, 3, 4 };
  double coords[4*5] = { }; // I think forces to initialize to zeros
  
  // create a dataset
  rc = fc_createDataset("dataset", &dataset);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create dataset");
  
  // loop over dimensions
  for (i = 1; i <= 3; i++) {
    int otherNumDim = i+1;
    FC_MathType mathType = FC_MT_VECTOR, otherMathType = FC_MT_TENSOR;
    if (i == 1) 
      mathType = FC_MT_SCALAR;
    if (i == 3)
      otherNumDim = 2;

    // create mesh & "perfect" displacement var
    fc_createMesh(dataset, "mesh", &mesh);
    rc = fc_setMeshCoords(mesh, i, numVert, coords);
    fail_unless(rc == FC_SUCCESS, "abort: failed to set coords");
    rc = fc_setMeshElementConns(mesh, FC_ET_POINT, numVert, conns);
    fail_unless(rc == FC_SUCCESS, "abort: failed to set conns");
    fc_createVariable(mesh, "var", &var);
    rc = fc_setVariableData(var, numVert, i, FC_AT_VERTEX, mathType,
			    FC_DT_DOUBLE, coords);
    fail_unless(rc == FC_SUCCESS, "abort: failed to set variable data");
 
    // test "perfect" displacement var
    fail_unless(fc_isValidDisplacementVariable(mesh, var) == 1,
		"perfect displacement var should pass");

    // create otherMesh that is the same as current mesh
    fc_createMesh(dataset, "otherMesh", &otherMesh);
    rc = fc_setMeshCoords(otherMesh, i, numVert, coords);
    fail_unless(rc == FC_SUCCESS, "abort: failed to set coords");
    rc = fc_setMeshElementConns(otherMesh, FC_ET_POINT, numVert, conns);
    fail_unless(rc == FC_SUCCESS, "abort: failed to set conns");
   
    // A different mesh is o.k.
    fc_createVariable(otherMesh, "otherVar", &otherVar);
    rc = fc_setVariableData(otherVar, numVert, i, FC_AT_VERTEX, mathType,
			    FC_DT_DOUBLE, coords);
    fail_unless(rc == FC_SUCCESS, "abort: failed to set variable data");
    fail_unless(fc_isValidDisplacementVariable(mesh, otherVar) == 1,
		"perfect displacment var on another mesh should pass");

    // A different assoc is o.k.
    fc_createVariable(otherMesh, "otherVar", &otherVar);
    rc = fc_setVariableData(otherVar, numVert, i, FC_AT_ELEMENT, mathType,
			    FC_DT_DOUBLE, coords);
    fail_unless(rc == FC_SUCCESS, "abort: failed to set variable data");
    fail_unless(fc_isValidDisplacementVariable(mesh, otherVar) == 1,
		"displacement var w/ diff assoc can still be valid");

    // A different mathtype is o.k.
    if (i != 1) { // can't really test for 1D
      fc_createVariable(otherMesh, "otherVar", &otherVar);
      rc = fc_setVariableData(otherVar, numVert, i, FC_AT_VERTEX,
			      otherMathType, FC_DT_DOUBLE, coords);
      fail_unless(rc == FC_SUCCESS, "abort: failed to set variable data");
      fail_unless(fc_isValidDisplacementVariable(mesh, otherVar) == 1,
		  "displacement var w/ diff math type can still be valid");
    }

    // A different datatype is o.k.
    fc_createVariable(otherMesh, "otherVar", &otherVar);
    rc = fc_setVariableData(otherVar, numVert, i, FC_AT_VERTEX, mathType,
			    FC_DT_DOUBLE, coords);
    fail_unless(rc == FC_SUCCESS, "abort: failed to set variable data");
    fail_unless(fc_isValidDisplacementVariable(mesh, otherVar) == 1,
		"displacment var w/ diff datatype can still be valid");

    // done with otherMesh that is just like this current mesh
    fc_deleteMesh(otherMesh);

    // if numVert != numDataPoint, should fail
    fc_createMesh(dataset, "otherMesh", &otherMesh);
    rc = fc_setMeshCoords(otherMesh, i, numVert-1, coords);
    fail_unless(rc == FC_SUCCESS, "abort: failed to set coords");
    rc = fc_setMeshElementConns(otherMesh, FC_ET_POINT, numVert-1, conns);
    fail_unless(rc == FC_SUCCESS, "abort: failed to set conns");
    fc_createVariable(otherMesh, "otherVar", &otherVar);
    rc = fc_setVariableData(otherVar, numVert-1, i, FC_AT_VERTEX, mathType,
			    FC_DT_DOUBLE, coords);
    fail_unless(rc == FC_SUCCESS, "abort: failed to set variable data");
    fail_unless(fc_isValidDisplacementVariable(mesh, otherVar) == 0,
		"if numVert != numDataPoint, should fail");
    fc_deleteMesh(otherMesh);

    // if numDim != numComp, should fail
    fc_createMesh(dataset, "otherMesh", &otherMesh);
    rc = fc_setMeshCoords(otherMesh, otherNumDim, numVert, coords);
    fail_unless(rc == FC_SUCCESS, "abort: failed to set coords");
    rc = fc_setMeshElementConns(otherMesh, FC_ET_POINT, numVert, conns);
    fail_unless(rc == FC_SUCCESS, "abort: failed to set conns");
    fc_createVariable(otherMesh, "otherVar", &otherVar);
    rc = fc_setVariableData(otherVar, numVert, otherNumDim, FC_AT_VERTEX, 
			    FC_MT_VECTOR, FC_DT_DOUBLE, coords);
    fail_unless(rc == FC_SUCCESS, "abort: failed to set variable data");
    fail_unless(fc_isValidDisplacementVariable(mesh, otherVar) == 0,
		"if numDim != numComp, should fail");
    fc_deleteMesh(otherMesh);

    // error conditions! bad mesh or bad variable
    fail_unless(fc_isValidDisplacementVariable(mesh, var) == 1,
		"abort: was making sure that mesh and var still o.k.");
    fail_unless(fc_isValidDisplacementVariable(badMesh, var) == 0,
		"variable cannot be valid for a bad mesh");
    fail_unless(fc_isValidDisplacementVariable(mesh, badVar) == 0,
		"an invalid variable cannot be valid");
  }

  // all done
  rc = fc_deleteDataset(dataset);
  fail_unless(rc == FC_SUCCESS, "failed to delete dataset at end of tests");
}
END_TEST
 
// test calculation of displaced coords
START_TEST(coords_displ)
{
  FC_ReturnCode rc;
  int i;
  int numDataPoint, numComp;
  FC_Dataset dataset;
  FC_Mesh *returnMeshes;
  int numReturnMeshes;
  FC_Mesh mesh, badMesh = { 999, 999 };
  FC_Variable coordsVar, edgeLengths, badVar = { 999, 999 };
  FC_Variable intCoordsVar, floatCoordsVar, charCoordsVar;
  FC_AssociationType assoc;
  FC_MathType mathtype;
  double* coords_p, *displ_coords;
  char* charCoords;
  int* intCoords;
  float* floatCoords;
  
  // --- setup

  fc_loadDataset("../data/gen_hex.ex2", &dataset);
  rc = fc_getMeshByName(dataset, "hex mesh",&numReturnMeshes,&returnMeshes);
  fail_unless(rc == FC_SUCCESS, "aborted - cant get mesh by name");
  fail_unless(numReturnMeshes == 1, "failed to find unique mesh by name");
  mesh = returnMeshes[0];
  free(returnMeshes);

  rc = fc_getMeshCoordsPtr(mesh, &coords_p);
  // get mesh coords (double)
  rc = fc_getMeshCoordsAsVariable(mesh, &coordsVar);
  fail_unless(rc == FC_SUCCESS, "abort: failed to get coords as var");
  // get a variable with different association
  rc = fc_getEdgeLengths(mesh, &edgeLengths);
  fail_unless(rc == FC_SUCCESS, "abort: failed to get edge lengths");
  // create some variables of different data types
  fc_getVariableInfo(coordsVar, &numDataPoint, &numComp, &assoc, &mathtype,
		     NULL);
  charCoords = malloc(numDataPoint*numComp*sizeof(char));
  intCoords = malloc(numDataPoint*numComp*sizeof(int));
  floatCoords = malloc(numDataPoint*numComp*sizeof(float));
  for (i = 0; i < numDataPoint*numComp; i++) {
    charCoords[i] = 'a';
    intCoords[i] = coords_p[i];
    floatCoords[i] = coords_p[i];
  }
  fc_createVariable(mesh, "char coords", &charCoordsVar);
  fc_createVariable(mesh, "int coords", &intCoordsVar);
  fc_createVariable(mesh, "float coords", &floatCoordsVar);
  rc = fc_setVariableDataPtr(charCoordsVar, numDataPoint, numComp, assoc,
			     mathtype, FC_DT_CHAR, (void*)charCoords);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create char var");
  rc = fc_setVariableDataPtr(intCoordsVar, numDataPoint, numComp, assoc,
			     mathtype, FC_DT_INT, (void*)intCoords);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create int var");
  rc = fc_setVariableDataPtr(floatCoordsVar, numDataPoint, numComp, assoc,
			     mathtype, FC_DT_FLOAT, (void*)floatCoords);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create float var");

  // --- test good args

  // if we displace with coords, should get 2*'s value
  rc = fc_getDisplacedMeshCoords(mesh, coordsVar, &displ_coords);
  fail_unless(rc == FC_SUCCESS, "shouldn't fail");
  for (i = 0; i < numDataPoint*numComp; i++)
    fail_unless(displ_coords[i] == 2.*coords_p[i],
                "mismatch of displaced coords");
  free(displ_coords);

  // should still work if displ_var is int or float
  rc = fc_getDisplacedMeshCoords(mesh, intCoordsVar, &displ_coords);
  fail_unless(rc == FC_SUCCESS, "shouldn't fail");
  for (i = 0; i < numDataPoint*numComp; i++)
    fail_unless(FC_DBL_EQUIV(displ_coords[i], (int)coords_p[i] + coords_p[i]),
    "mismatch of displaced coords");
  free(displ_coords);
  rc = fc_getDisplacedMeshCoords(mesh, floatCoordsVar, &displ_coords);
  fail_unless(rc == FC_SUCCESS, "shouldn't fail");
  for (i = 0; i < numDataPoint*numComp; i++)
    fail_unless(FC_FLT_EQUIV(displ_coords[i], 2.*coords_p[i]),
                "mismatch of displaced coords");
  free(displ_coords);

  // --- test bad args

  // shouldn't work if displ_var is a char
  rc = fc_getDisplacedMeshCoords(mesh, charCoordsVar, &displ_coords);
  fail_unless(rc != FC_SUCCESS, "should fail with char var");
  fail_unless(displ_coords == NULL, "fail should return NULL");

  // bad mesh
  rc = fc_getDisplacedMeshCoords(badMesh, coordsVar, &displ_coords);
  fail_unless(rc != FC_SUCCESS, "should fail with bad mesh");
  fail_unless(displ_coords == NULL, "fail should return NULL");

  // bad var
  rc = fc_getDisplacedMeshCoords(mesh, badVar, &displ_coords);
  fail_unless(rc != FC_SUCCESS, "should fail with bad var");
  fail_unless(displ_coords == NULL, "fail should return NULL");
  rc = fc_getDisplacedMeshCoords(mesh, edgeLengths, &displ_coords);
  fail_unless(rc != FC_SUCCESS, "should fail if var not on vertices");
  fail_unless(displ_coords == NULL, "fail should return NULL");

  // no displ_coords
  rc = fc_getDisplacedMeshCoords(mesh, edgeLengths, NULL);
  fail_unless(rc != FC_SUCCESS, "should fail if no coords");

  // cleanup
  rc = fc_deleteDataset(dataset);
  fail_unless(rc == FC_SUCCESS, "failed to delete dataset at end of test");
}
END_TEST


// params

// quad param & locate
// each case is for a parameter location
// also test "above" version of each case where the point is not in the plane
START_TEST(quad_param) {
  int i, j, k;
  int numCase = 11;
  double params[11][2] = { { 0, 0 }, { 1, 0 }, { 1, 1 }, { 0, 1 }, // corners
			   { 0.5, 0 }, { 1, 0.5 }, { 0.5, 1 }, { 0, 0.5 }, // edges
			   { 0.5, 0.5 }, // center
			   { 0.0001, 0.0002 },  // arbitrary middle near a
						 // problem point (0,0)
			   { 0.9999, 0.9998 } };
  FC_Coords rectQuad[4], rotRect[4], funkyQuad[4], rotFunky[4];
  FC_Coords point, apoint;
  FC_Coords start = { 1.1, 1.2, 0 };
  double xwidth = 1.2;
  double ywidth = 3.3;
  double funk = 1.1;
  FC_Coords temp_point, *point_p;
  int temp_numParam;
  double *temp_params;

  //---setup quads

  // rectangle in x,y plane
  for (i = 0; i < 4; i++)
    for (j = 0; j < 3; j++)
      rectQuad[i][j] = start[j];
  rectQuad[1][0] += xwidth;
  rectQuad[2][0] += xwidth;
  rectQuad[2][1] += ywidth;
  rectQuad[3][1] += ywidth;

  // funky quad in x,y plane
  for (i = 0; i < 4; i++)
    for (j = 0; j < 3; j++)
      funkyQuad[i][j] = start[j];
  funkyQuad[1][0] += xwidth*funk;
  funkyQuad[1][1] += funk;
  funkyQuad[2][0] += xwidth;
  funkyQuad[2][1] += ywidth;
  funkyQuad[3][0] += funk;
  funkyQuad[3][1] += ywidth*funk;

  // make arbitrary rotated versions
  for (i = 0; i < 4; i++) {
    for (j = 0; j < 3; j++) {
      rotRect[i][j] = rectQuad[i][j];
      rotFunky[i][j] = funkyQuad[i][j];
    }
    rotate3D(&rotRect[i][0], &rotRect[i][1], &rotRect[i][2]);
    rotate3D(&rotFunky[i][0], &rotFunky[i][1], &rotFunky[i][2]);
  }
  
  //---rectangle quad in x,y plane
  for (i = 0; i < numCase; i++) {
    for (j = 0; j < 3; j++)
      point[j] = rectQuad[0][j];
    point[0] += xwidth*params[i][0];
    point[1] += ywidth*params[i][1];
    for (j = 0; j < 2; j++)
      apoint[j] = point[j];
    apoint[2] = .1;
    
    _fc_calcQuadLocation(rectQuad, 2, params[i], &temp_point);
    //printf("%d: point = %g %g %g\n", i, point[0], point[1], point[2]);
    //printf("%d: temp  = %g %g %g\n", i, temp_point[0], temp_point[1], 
    //temp_point[2]);
    for (j = 0; j < 3; j++) {
      fail_unless(FC_VALUE_EQUIV(temp_point[j], point[j], 2*DBL_EPSILON, 
				 DBL_MIN), "mismatch of located point");
    }

    // 1st case, in plane, 2nd case, above plane    
    for (k = 0; k < 2; k++) { 
      if (k == 0)
	point_p = &point;
      else
	point_p = &apoint;

      _fc_calcQuadParams(rectQuad, *point_p, &temp_numParam, &temp_params);
      fail_unless(temp_params != NULL, "failed to converge");
      fail_unless(temp_numParam == 2, "should find two params");
      //printf("%d: params = %.20g %.20g\n", i, params[i][0], params[i][1]);
      //printf("%d: temp   = %.20g %.20g\n", i, temp_params[0],temp_params[1]);
      for (j = 0; j < 2; j++) {
	if (params[i][j] > .1)
	  fail_unless(FC_VALUE_EQUIV(temp_params[j], params[i][j], 
				     4*DBL_EPSILON, DBL_MIN), 
		      "mismatch of params");
	else 
	  fail_unless(fabs(temp_params[j] - params[i][j]) < DBL_EPSILON, 
		      "mismatch of params near zero");
      }
      free(temp_params);
    }
  }

  //---points outside quad in x,y plane
  for (i = 0; i < 2; i++)
    point[i] = rectQuad[0][i] - 10*DBL_EPSILON;
  point[2] = 0;
  for (i = 0; i < 2; i++)
    apoint[i] = point[i];
  apoint[2] = .1;
  for (k = 0; k < 2; k++) {
    if (k == 0)
      point_p = &point;
    else
      point_p = &apoint;
    _fc_calcQuadParams(rectQuad, *point_p, &temp_numParam, &temp_params);
    fail_unless(temp_params != NULL, "?? expect to converge??");
    for (i = 0; i < 2; i++)
      fail_unless(fc_ltd(temp_params[i], 0), "expect to be less than zero");
    free(temp_params);
  }
  for (i = 0; i < 2; i++)
    point[i] = rectQuad[2][i] + 10*DBL_EPSILON;
  point[2] = 0;
  for (i = 0; i < 2; i++)
    apoint[i] = point[i];
  apoint[2] = .1;
  for (k = 0; k < 2; k++) {
    if (k == 0)
      point_p = &point;
    else
      point_p = &apoint;
    _fc_calcQuadParams(rectQuad, *point_p, &temp_numParam, &temp_params);
    fail_unless(temp_params != NULL, "?? expect to converge??");
    for (i = 0; i < 2; i++)
      fail_unless(fc_gtd(temp_params[i], 1), "expect to be greater than 1");
    free(temp_params);
  }

  //--rotated rectangle quad
  for (i = 0; i < numCase; i++) {
    for (j = 0; j < 3; j++) 
      point[j] = rectQuad[0][j];
    point[0] += xwidth*params[i][0];
    point[1] += ywidth*params[i][1];
    for (j = 0; j < 2; j++)
      apoint[j] = point[j];
    apoint[2] = .1;
    rotate3D(&point[0], &point[1], &point[2]);
    rotate3D(&apoint[0], &apoint[1], &apoint[2]);

    _fc_calcQuadLocation(rotRect, 2, params[i], &temp_point);
    //printf("%d: point = %g %g %g\n", i, point[0], point[1], point[2]);
    //printf("%d: temp  = %g %g %g\n", i, temp_point[0], temp_point[1], 
    //temp_point[2]);
    for (j = 0; j < 3; j++) {
      fail_unless(FC_VALUE_EQUIV(temp_point[j], point[j], 4*DBL_EPSILON, 
				 DBL_MIN), "mismatch of located point");
    }

    for (k = 0; k < 2; k++) {
      if (k == 0)
	point_p = &point;
      else
	point_p = &apoint;

      _fc_calcQuadParams(rotRect, *point_p, &temp_numParam, &temp_params);
      fail_unless(temp_params != NULL, "failed to converge");
      //printf("%d: params = %.20g %.20g\n", i, params[i][0], params[i][1]);
      //printf("%d: temp   = %.20g %.20g\n", i, temp_params[0], temp_params[1]);
      fail_unless(temp_numParam == 2, "should find two params");
      for (j = 0; j < 2; j++) {
	if (params[i][j] > .1)
	  fail_unless(FC_VALUE_EQUIV(temp_params[j], params[i][j], 
				     4*DBL_EPSILON, DBL_MIN), 
		      "mismatch of params");
	else 
	  fail_unless(fabs(temp_params[j] - params[i][j]) < 2*DBL_EPSILON, 
		      "mismatch of params near zero");
      }
      free(temp_params);
    }
  }
    
  //---points outside rotated quad 
  for (i = 0; i < 2; i++)
    point[i] = rectQuad[0][i] - 10*DBL_EPSILON;
  point[2] = 0;
  for (i = 0; i < 2; i++)
    apoint[i] = point[i];
  apoint[2] = .1;
  rotate3D(&point[0], &point[1], &point[2]);
  rotate3D(&apoint[0], &apoint[1], &apoint[2]);
  for (k = 0; k < 2; k++) {
    if (k == 0)
      point_p = &point;
    else
      point_p = &apoint;
    _fc_calcQuadParams(rotRect, *point_p, &temp_numParam, &temp_params);
    fail_unless(temp_params != NULL, "?? expect to converge??");
    for (i = 0; i < 2; i++)
      fail_unless(fc_ltd(temp_params[i], 0), "expect to be less than zero");
    free(temp_params);
  }
  for (i = 0; i < 2; i++)
    point[i] = rectQuad[2][i] + 10*DBL_EPSILON;
  point[2] = 0;
  for (i = 0; i < 2; i++)
    apoint[i] = point[i];
  apoint[2] = .1;
  rotate3D(&point[0], &point[1], &point[2]);
  rotate3D(&apoint[0], &apoint[1], &apoint[2]);
  for (k = 0; k < 2; k++) {
    if (k == 0)
      point_p = &point;
    else
      point_p = &apoint;
    _fc_calcQuadParams(rotRect, *point_p, &temp_numParam, &temp_params);
    fail_unless(temp_params != NULL, "?? expect to converge??");
    for (i = 0; i < 2; i++)
      fail_unless(fc_gtd(temp_params[i], 1), "expect to be greater than 1");
    free(temp_params);
  }

  //---funky quad in x,y plane
  for (i = 0; i < 9; i++) {
    if (i < 4) {
      for (j = 0; j < 2; j++)
	point[j] = funkyQuad[i][j];
    }
    else if (i < 8) {
      FC_Coords *point1_p, *point2_p;
      if (i == 4) {
	point1_p = &funkyQuad[0];
	point2_p = &funkyQuad[1];
      }
      else if (i == 5) {
	point1_p = &funkyQuad[1];
	point2_p = &funkyQuad[2];
      }
      else if (i == 6) {
	point1_p = &funkyQuad[3];
	point2_p = &funkyQuad[2];
      }
      else if (i == 7) {
	point1_p = &funkyQuad[0];
	point2_p = &funkyQuad[3];
      }
      for (j = 0; j < 2; j++) {
	double diff = (*point2_p)[j] - (*point1_p)[j];
	point[j] = (*point1_p)[j] + 0.5*diff;
      }
    }
    else if (i == 8) {
      for (j = 0; j < 2; j++)
	point[j] = 0;
      for (j = 0; j < 4; j++)
	for (k = 0; k < 2; k++)
	  point[k] += funkyQuad[j][k];
      for (j = 0; j < 2; j++)
	point[j] /= 4.;
    }
    // FIX last cases?
    point[2] = 0;
    for (j = 0; j < 2; j++)
      apoint[j] = point[j];
    apoint[2] = .1;

    _fc_calcQuadLocation(funkyQuad, 2, params[i], &temp_point);
    //printf("%d: point = %g %g %g\n", i, point[0], point[1], point[2]);
    //printf("%d: temp  = %g %g %g\n", i, temp_point[0], temp_point[1], temp_point[2]);
    if (i < 8) {
      for (j = 0; j < 3; j++) {
	fail_unless(FC_VALUE_EQUIV(temp_point[j], point[j], 2*DBL_EPSILON, DBL_MIN), 
		    "mismatch of located point");
      }
    }

    // 1st case, in plane, 2nd case, above plane    
    for (k = 0; k < 2; k++) {
      if (k == 0)
	point_p = &point;
      else
	point_p = &apoint;

      _fc_calcQuadParams(funkyQuad, *point_p, &temp_numParam, &temp_params);
      fail_unless(temp_params != NULL, "failed to converge");
      //printf("%d: params = %g %g\n", i, params[i][0], params[i][1]);
      //printf("%d: temp   = %g %g\n", i, temp_params[0], temp_params[1]);
      for (j = 0; j < 2; j++) {
	if (params[i][j] > .1)
	  fail_unless(FC_VALUE_EQUIV(temp_params[j], params[i][j], 
				     20*DBL_EPSILON, DBL_MIN), 
		      "mismatch of params");
	else 
	  fail_unless(fabs(temp_params[j] - params[i][j]) < 2*DBL_EPSILON, 
		      "mismatch of params near zero");
      }
      free(temp_params);
    }
  }
 
  //---points outside funky quad in x,y plane
  for (i = 0; i < 2; i++)
    point[i] = funkyQuad[0][i] - 100*DBL_EPSILON;
  point[2] = 0;
  for (i = 0; i < 2; i++)
    apoint[i] = point[i];
  apoint[2] = .1;
  for (k = 0; k < 2; k++) {
    if (k == 0)
      point_p = &point;
    else
      point_p = &apoint;
    _fc_calcQuadParams(funkyQuad, *point_p, &temp_numParam, &temp_params);
    fail_unless(temp_params != NULL, "?? expect to converge??");
    for (i = 0; i < 2; i++)
      fail_unless(fc_ltd(temp_params[i], 0), "expect to be less than zero");
    free(temp_params);
  }
  for (i = 0; i < 2; i++)
    point[i] = funkyQuad[2][i] + 100*DBL_EPSILON;
  point[2] = 0;
  for (i = 0; i < 2; i++)
    apoint[i] = point[i];
  apoint[2] = .1;
  for (k = 0; k < 2; k++) {
    if (k == 0)
      point_p = &point;
    else
      point_p = &apoint;
    _fc_calcQuadParams(funkyQuad, *point_p, &temp_numParam, &temp_params);
    fail_unless(temp_params != NULL, "?? expect to converge??");
    for (i = 0; i < 2; i++)
      fail_unless(fc_gtd(temp_params[i], 1), "expect to be greater than 1");
    free(temp_params);
  }

  //---rotated funky quad 
  for (i = 0; i < 9; i++) {
    if (i < 4)
      for (j = 0; j < 2; j++)
	point[j] = funkyQuad[i][j];
    else if (i < 8) {
      FC_Coords *point1_p, *point2_p;
      if (i == 4) {
	point1_p = &funkyQuad[0];
	point2_p = &funkyQuad[1];
      }
      else if (i == 5) {
	point1_p = &funkyQuad[1];
	point2_p = &funkyQuad[2];
      }
      else if (i == 6) {
	point1_p = &funkyQuad[3];
	point2_p = &funkyQuad[2];
      }
      else if (i == 7) {
	point1_p = &funkyQuad[0];
	point2_p = &funkyQuad[3];
      }
      for (j = 0; j < 2; j++) {
	double diff = (*point2_p)[j] - (*point1_p)[j];
	point[j] = (*point1_p)[j] + 0.5*diff;
      }
    }
    else if (i == 8) {
      for (j = 0; j < 2; j++)
	point[j] = 0;
      for (j = 0; j < 4; j++)
	for (k = 0; k < 2; k++)
	  point[k] += funkyQuad[j][k];
      for (j = 0;j < 2; j++)
	point[j] /= 4.;
    }
    // FIX last case?
    point[2] = 0;
    for (j = 0; j < 2; j++)
      apoint[j] = point[j];
    apoint[2] = .1;
    rotate3D(&point[0], &point[1], &point[2]);
    rotate3D(&apoint[0], &apoint[1], &apoint[2]);
    
    _fc_calcQuadLocation(rotFunky, 2, params[i], &temp_point);
    //printf("%d: point = %g %g %g\n", i, point[0], point[1], point[2]);
    //printf("%d: temp  = %g %g %g\n", i, temp_point[0], temp_point[1], temp_point[2]);
    if (i < 8) {
      for (j = 0; j < 3; j++) {
	fail_unless(FC_VALUE_EQUIV(temp_point[j], point[j], 2*DBL_EPSILON, DBL_MIN), 
		    "mismatch of located point");
      }
    }

    // 1st case, in plane, 2nd case, above plane    
    for (k = 0; k < 2; k++) { 
      if (k == 0)
	point_p = &point;
      else
	point_p = &apoint;

      if (i == 3 && k == 1)
	continue; // FIX!
      _fc_calcQuadParams(rotFunky, *point_p, &temp_numParam, &temp_params);
      fail_unless(temp_params != NULL, "failed to converge");
      //printf("%d: params = %.20g %.20g\n", i, params[i][0], params[i][1]);
      //printf("%d: temp   = %.20g %.20g\n", i, temp_params[0], temp_params[1]);
      fail_unless(temp_numParam == 2, "should find two params");
      for (j = 0; j < 2; j++) {
	if (params[i][j] > .1)
	  fail_unless(FC_VALUE_EQUIV(temp_params[j], params[i][j], 
				     25*DBL_EPSILON, DBL_MIN), 
		      "mismatch of params");
	else 
	  fail_unless(fabs(temp_params[j] - params[i][j]) < 2*DBL_EPSILON, 
		      "mismatch of params near zero");
      }
      free(temp_params);
    }
  }

  //---points outside rotated funky quad 
  for (i = 0; i < 2; i++)
    point[i] = funkyQuad[0][i] - 100*DBL_EPSILON;
  point[2] = 0;
  for (i = 0; i < 2; i++)
    apoint[i] = point[i];
  apoint[2] = .1;
  rotate3D(&point[0], &point[1], &point[2]);
  rotate3D(&apoint[0], &apoint[1], &apoint[2]);
  for (k = 0; k < 2; k++) {
    if (k == 0)
      point_p = &point;
    else
      point_p = &apoint;
    _fc_calcQuadParams(rotFunky, *point_p, &temp_numParam, &temp_params);
    fail_unless(temp_params != NULL, "?? expect to converge??");
    for (i = 0; i < 2; i++)
    fail_unless(fc_ltd(temp_params[i], 0), "expect to be less than zero");
    free(temp_params);
  }
  for (i = 0; i < 2; i++)
    point[i] = funkyQuad[2][i] + 10*DBL_EPSILON;
  point[2] = 0;
  for (i = 0; i < 2; i++)
    apoint[i] = point[i];
  apoint[2] = .1;
  rotate3D(&point[0], &point[1], &point[2]);
  rotate3D(&apoint[0], &apoint[1], &apoint[2]);
  for (k = 0; k < 2; k++) {
    if (k == 0)
      point_p = &point;
    else
      point_p = &apoint;
    _fc_calcQuadParams(rotFunky, *point_p, &temp_numParam, &temp_params);
    fail_unless(temp_params != NULL, "?? expect to converge??");
    for (i = 0; i < 2; i++)
      fail_unless(fc_gtd(temp_params[i], 1), "expect to be greater than 1");
    free(temp_params);
  }
}
END_TEST


// test bounding box operations
START_TEST(bb_ops)
{
  FC_ReturnCode rc;
  int i, j, k;
  FC_Coords low = { 1.0, 1.1, 1.2 };
  FC_Coords low2 = { 0.9999999999999999, 1.1, 1.2 }; // ~equal to low1
  FC_Coords mid1 = { 3.1, 3.11, 3.111 };
  FC_Coords mid2 = { 6.2, 6.22, 6.222 };
  FC_Coords high = { 9.9, 9.7, 9.5 };
  int numCase = 8;
  // Cases: 0) the same, 1) 1 inside 2, on boundary, 2) 1 inside 2, not on
  // boundary, 3) 2 inside 1, on boundary, 4) 2 inside1, not on boundary, 
  // 5) no intersection, 6) partial overlap, 7) partial overlap - swap args
  int overlap_flags[8] = { 3, 1, 1, 2, 2, 0, 4, 4 }; 
  FC_Coords* lows1[8] =  { &low,  &low,  &mid1, &low,  &low,  &low,  &low, &mid1 };
  FC_Coords* highs1[8] = { &high, &mid2, &mid2, &high, &high, &mid1, &mid2, &high };
  FC_Coords* lows2[8] =  { &low,  &low,  &low,  &mid1, &mid1, &mid2, &mid1, &low };
  FC_Coords* highs2[8] = { &high, &high, &high, &high, &mid2, &high, &high, &mid2 };
  FC_Coords* overlap_lows[8] =  { &low, &low, &mid1, &mid1, &mid1, NULL, &mid1, &mid1 };
  FC_Coords* overlap_highs[8] = { &high, &mid2, &mid2, &high, &mid2, NULL, &mid2, &mid2 };
  // combine is always low/high
  int temp_flag;
  FC_Coords temp_low, temp_high;
  
  // test fc_isBoundingBoxValid()
  temp_flag = fc_isBoundingBoxValid(3, low, high);
  fail_unless(temp_flag == 1, "bb should be valid");
  temp_flag = fc_isBoundingBoxValid(3, low, low);
  fail_unless(temp_flag == 1, "bb with lowpoint = highpoint should be valid");
  temp_flag = fc_isBoundingBoxValid(3, high, low);
  fail_unless(temp_flag == 0, "reversed bb should not be valid");
  temp_flag = fc_isBoundingBoxValid(0, low, high); 
  fail_unless(temp_flag == 0, "bb should not be valid of numDim = 0");
  temp_flag = fc_isBoundingBoxValid(4, low, high); 
  fail_unless(temp_flag == 0, "bb should not be valid of numDim = 4");
  temp_flag = fc_isBoundingBoxValid(3, NULL, high); 
  fail_unless(temp_flag == 0, "bb should not be valid if lowpoint is null");
  temp_flag = fc_isBoundingBoxValid(3, low, NULL); 
  fail_unless(temp_flag == 0, "bb should not be valid if highpoint is null");
  temp_flag = fc_isBoundingBoxValid(3, low, low2);
  fail_unless(temp_flag == 1, "bb should be valid if doing floats right");

  // loop over cases to test overlap & combine bounding boxes
  for (i = 0; i < numCase; i++) {
    // loop over possible dimensions
    for (j = 1; j <= 3; j++) {
      
      // test fc_getBoundingBoxesOverlap()
      rc = fc_getBoundingBoxesOverlap(j, *lows1[i], *highs1[i], *lows2[i], 
                                      *highs2[i], &temp_flag, &temp_low, 
                                      &temp_high); 
      fail_unless(rc == FC_SUCCESS, "failed to get overlap");
      fail_unless(temp_flag == overlap_flags[i], "mismatch of flag");
      if (temp_flag > 0) {
        for (k = 0; k < j; k++) {
          fail_unless(temp_low[k] == (*overlap_lows[i])[k], 
                      "mismatch of low point");
          fail_unless(temp_high[k] == (*overlap_highs[i])[k], 
                      "mismatch of high point");
        }
        for (k = j; k < 3; k++) {
          fail_unless(temp_low[k] == 0., "unused should be 0");
          fail_unless(temp_high[k] == 0., "unused should be 0");
        }
      }
      else { // no overlap
        for (k = 0; k < 3; k++) {
          fail_unless(temp_low[k] == -1., "not overlap should be -1");
          fail_unless(temp_high[k] == -2., "not overlap should be -1");
        }
      }

      // test fc_combineBoundingBoxes()
      rc = fc_combineBoundingBoxes(j, *lows1[i], *highs1[i], *lows2[i], 
                                   *highs2[i], &temp_low, &temp_high); 
      fail_unless(rc == FC_SUCCESS, "failed to combine");
      // combine is always low/high
      for (k = 0; k < j; k++) {
        fail_unless(temp_low[k] == low[k], "mismatch of low point");
        fail_unless(temp_high[k] == high[k], "mismatch of high point");
      }
      for (k = j; k < 3; k++) {
        fail_unless(temp_low[k] == 0., "unused should be 0");
        fail_unless(temp_high[k] == 0., "unused should be 0");
      }
    }
  }

  // --- test special cases

  // fc_getBoundingBoxOverlap() -- o.k. not to ask for overlap bb
  rc = fc_getBoundingBoxesOverlap(3, low, high, low, high, &temp_flag, 
                                  NULL, &temp_high); 
  fail_unless(rc == FC_SUCCESS, "shoudl ont fail with null overlap_low");
  fail_unless(temp_flag == 3, "mismatch of flag");
  for (i = 0; i < 3; i++)
    fail_unless(temp_high[i] == high[i], "mismatch of high point");
  rc = fc_getBoundingBoxesOverlap(3, low, high, low, high, &temp_flag, 
                                  &temp_low, NULL); 
  fail_unless(rc == FC_SUCCESS, "shoudl ont fail with null overlap_low");
  fail_unless(temp_flag == 3, "mismatch of flag");
  for (i = 0; i < 3; i++)
    fail_unless(temp_low[i] == low[i], "mismatch of low point");

  // --- test errors

  // fc_getBoundingBoxesOverlap() -- bad args
  rc = fc_getBoundingBoxesOverlap(0, low, high, low, high, &temp_flag,
                                  &temp_low, &temp_high);
  fail_unless(rc != FC_SUCCESS, "should fail if numDim < 1");
  fail_unless(temp_flag == -1, "fail should return null");
  for (i = 0; i < 3; i++) 
    fail_unless(temp_low[i] == -1. && temp_high[i] == -2., 
                "fail should return null");
  rc = fc_getBoundingBoxesOverlap(4, low, high, low, high, &temp_flag,
                                  &temp_low, &temp_high);
  fail_unless(rc != FC_SUCCESS, "should fail if numDim > 3");
  fail_unless(temp_flag == -1, "fail should return null");
  for (i = 0; i < 3; i++)
    fail_unless(temp_low[i] == -1. && temp_high[i] == -2., 
                "fail should return null");
  rc = fc_getBoundingBoxesOverlap(3, NULL, high, low, high, &temp_flag,
                                  &temp_low, &temp_high);
  fail_unless(rc != FC_SUCCESS, "should fail if no low1");
  fail_unless(temp_flag == -1, "fail should return null");
  for (i = 0; i < 3; i++) 
    fail_unless(temp_low[i] == -1. && temp_high[i] == -2., 
                "fail should return null");
  rc = fc_getBoundingBoxesOverlap(3, low, NULL, low, high, &temp_flag,
                                  &temp_low, &temp_high);
  fail_unless(rc != FC_SUCCESS, "should fail if no high1");
  fail_unless(temp_flag == -1, "fail should return null");
  for (i = 0; i < 3; i++) 
    fail_unless(temp_low[i] == -1. && temp_high[i] == -2., 
                "fail should return null");
  rc = fc_getBoundingBoxesOverlap(3, low, high, NULL, high, &temp_flag,
                                  &temp_low, &temp_high);
  fail_unless(rc != FC_SUCCESS, "should fail if no low2");
  fail_unless(temp_flag == -1, "fail should return null");
  for (i = 0; i < 3; i++) 
    fail_unless(temp_low[i] == -1. && temp_high[i] == -2., 
                "fail should return null");
  rc = fc_getBoundingBoxesOverlap(3, low, high, low, NULL, &temp_flag,
                                  &temp_low, &temp_high);
  fail_unless(rc != FC_SUCCESS, "should fail if no high2");
  fail_unless(temp_flag == -1, "fail should return null");
  for (i = 0; i < 3; i++) 
     fail_unless(temp_low[i] == -1. && temp_high[i] == -2., 
                "fail should return null");
  rc = fc_getBoundingBoxesOverlap(3, low, high, low, high, NULL, 
                                  &temp_low,  &temp_high);
  fail_unless(rc != FC_SUCCESS, "should fail if no overlap_flag");
  fail_unless(temp_flag == -1, "fail should return null");
  for (i = 0; i < 3; i++) 
     fail_unless(temp_low[i] == -1. && temp_high[i] == -2., 
                "fail should return null");

  // fc_combineBoundingBoxes() -- bad args
  rc = fc_combineBoundingBoxes(0, low, high, low, high, &temp_low, &temp_high);
  fail_unless(rc != FC_SUCCESS, "should fail if numDim < 1");
  fail_unless(temp_flag == -1, "fail should return null");
  for (i = 0; i < 3; i++) 
    fail_unless(temp_low[i] == -1. && temp_high[i] == -2., 
                "fail should return null");
  rc = fc_combineBoundingBoxes(4, low, high, low, high, &temp_low, &temp_high);
  fail_unless(rc != FC_SUCCESS, "should fail if numDim > 3");
  fail_unless(temp_flag == -1, "fail should return null");
  for (i = 0; i < 3; i++)
    fail_unless(temp_low[i] == -1. && temp_high[i] == -2., 
                "fail should return null");
  rc = fc_combineBoundingBoxes(3, NULL, high, low, high, &temp_low, 
			       &temp_high);
  fail_unless(rc != FC_SUCCESS, "should fail if no low1");
  fail_unless(temp_flag == -1, "fail should return null");
  for (i = 0; i < 3; i++) 
    fail_unless(temp_low[i] == -1. && temp_high[i] == -2., 
                "fail should return null");
  rc = fc_combineBoundingBoxes(3, low, NULL, low, high, &temp_low, &temp_high);
  fail_unless(rc != FC_SUCCESS, "should fail if no high1");
  fail_unless(temp_flag == -1, "fail should return null");
  for (i = 0; i < 3; i++) 
    fail_unless(temp_low[i] == -1. && temp_high[i] == -2., 
                "fail should return null");
  rc = fc_combineBoundingBoxes(3, low, high, NULL, high, &temp_low, 
			       &temp_high);
  fail_unless(rc != FC_SUCCESS, "should fail if no low2");
  fail_unless(temp_flag == -1, "fail should return null");
  for (i = 0; i < 3; i++) 
    fail_unless(temp_low[i] == -1. && temp_high[i] == -2., 
                "fail should return null");
  rc = fc_combineBoundingBoxes(3, low, high, low, NULL, &temp_low, &temp_high);
  fail_unless(rc != FC_SUCCESS, "should fail if no high2");
  fail_unless(temp_flag == -1, "fail should return null");
  for (i = 0; i < 3; i++) 
     fail_unless(temp_low[i] == -1. && temp_high[i] == -2., 
                "fail should return null");
  rc = fc_combineBoundingBoxes(3, low, high, low, high, NULL, &temp_high);
  fail_unless(rc != FC_SUCCESS, "should fail if no overlap_low");
  fail_unless(temp_flag == -1, "fail should return null");
  for (i = 0; i < 3; i++) 
     fail_unless(temp_high[i] == -2., "fail should return null");
  rc = fc_combineBoundingBoxes(3, low, high, low, high, &temp_low, NULL);
  fail_unless(rc != FC_SUCCESS, "should fail if no overlap_high");
  fail_unless(temp_flag == -1, "fail should return null");
  for (i = 0; i < 3; i++) 
     fail_unless(temp_low[i] == -1., "fail should return null");
}
END_TEST

// test get mesh bounding box & get subset bounding box w/ &w/o displ
// basic layout of vertices is on two side-by-side hexes
// for 1D use only lower edge, for 2D use bottom faces
START_TEST(get_bb)
{
  FC_ReturnCode rc;
  int i, j, k, m;
  char name[20] = "blue berry", iName[30] = "inverted blue berry";
  FC_Dataset dataset;
  FC_Mesh mesh, iMesh, badMesh = { 999, 999 };
  FC_Subset subset, iSubset, badSubset = { 999, 999 };
  FC_Variable displVar, iDisplVar, badVar = { 999, 999 };
  int numElemType = 4;
  FC_ElementType elemTypes[4] = { FC_ET_POINT,  FC_ET_LINE,
                                  FC_ET_QUAD,   FC_ET_HEX };
  int topoDimPerType[4] = { 0, 1, 2, 3 };
  int numVertPerType[4] = { 3, 3, 6, 12}, numElemPerType[4] = { 3, 2, 2, 2 };
  int conns[4][16] = { { 0,    1,    2 },
                       { 0, 1,    1, 2 },
                       { 0, 1, 4, 3,    1, 2, 5, 4 },
                       { 0, 1, 4, 3, 6, 7, 10, 9,    1, 2, 5, 4, 7, 8, 11, 10
                       } };
  int dims[3] = { 3, 2, 2 };
  double comp_coords[3][3] = { { 1.1, 2.4, 3.9 },
                               { 0.3, 1.5 }, 
                               { 99., 99.9 } };
  double comp_displs[3][3] = { { 0.0, 0.01, 0.02 },
                               { 0.03, 0.04 },
                               { 0.05, 0.06 } };
  double *xcoords = comp_coords[0], *xdispls = comp_displs[0];
  double *ycoords = comp_coords[1], *ydispls = comp_displs[1];
  double *zcoords = comp_coords[2], *zdispls = comp_displs[2];
  double coords[3][3*12], displs[3][3*12]; // index by space dimensionality 
  // inverted coords = just multiply coords & displ by -1;
  // Had to add this so we could see vertex coords that decrease w/ incrs ID
  double iCoords[3*12], iDispls[3*12];
  int temp_dim;
  FC_Coords lowpoint, highpoint, iLowpoint, iHighpoint;
  FC_Coords nullPoint1 = { -1., -1., -1 }, nullPoint2 = { -2., -2., -2. };
  int numAssoc = 4, numEntity;
  FC_AssociationType assocs[4] = { FC_AT_VERTEX, FC_AT_EDGE, 
                                   FC_AT_FACE, FC_AT_ELEMENT };  
  FC_MathType mathType;

  // setup coords & displs
  // 1D coords
  for (i = 0; i < numVertPerType[1]; i++) {
    coords[0][i] = xcoords[i];
    displs[0][i] = xdispls[i];
  }
  // 2D coords
  for (i = 0; i < numVertPerType[2]; i++) {
    coords[1][i*2] = xcoords[i%dims[0]];
    coords[1][i*2+1] = ycoords[i/dims[0]];
    displs[1][i*2] = xdispls[i%dims[0]];
    displs[1][i*2+1] = ydispls[i/dims[0]];
  }
  // 3D coords
  for (i = 0; i < numVertPerType[3]; i++) {
    coords[2][i*3] = xcoords[i%dims[0]];
    coords[2][i*3+1] = ycoords[(i/dims[0])%dims[1]];
    coords[2][i*3+2] = zcoords[i/dims[0]/dims[1]];
    displs[2][i*3] = xdispls[i%dims[0]];
    displs[2][i*3+1] = ydispls[(i/dims[0])%dims[1]];
    displs[2][i*3+2] = zdispls[i/dims[0]/dims[1]];
  }

  // make dataset
  rc = fc_createDataset("temp dataset", &dataset);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create dataset");

  // test different types
  for (i = 0; i < numElemType; i++) { // loop over elem types
    for (j = 1; j <= 3; j++) { // loop over space dimensions
      // skip if it doesn't make sense to make this elem type at this dim
      if (topoDimPerType[i] > j)
        continue;

      // --- create data structures for testing: mesh & displVar

      // create a mesh & set coords & conns
      rc = fc_createMesh(dataset, name, &mesh);
      fail_unless(rc == FC_SUCCESS, "abort: failed to created mesh");
      rc = fc_setMeshCoords(mesh, j, numVertPerType[i], coords[j-1]);
      fail_unless(rc == FC_SUCCESS, "failed to set vertex coords");
      rc = fc_setMeshElementConns(mesh, elemTypes[i], numElemPerType[i], 
                                  conns[i]);
      fail_unless(rc == FC_SUCCESS, "failed to set element conns");
      // create displacement variable
      rc = fc_createVariable(mesh, "displ", &displVar);
      fail_unless(rc == FC_SUCCESS, "abort: failed to create variable");
      if (j == 1)
        mathType = FC_MT_SCALAR;
      else
        mathType = FC_MT_VECTOR;
      rc = fc_setVariableData(displVar, numVertPerType[i], j, FC_AT_VERTEX, 
                              mathType, FC_DT_DOUBLE, (void**)displs[j-1]);
      fail_unless(rc == FC_SUCCESS, "failed to set displVar data");
      // setup inverted mesh
      for (k = 0; k < numVertPerType[i]; k++) {
	for (m = 0; m < j; m++) {
	  iCoords[k*j + m] = -1.*coords[j-1][k*j + m];
	  iDispls[k*j + m] = -1.*displs[j-1][k*j + m];
	}
      }
      // create an inverted mesh & set coords & conns
      rc = fc_createMesh(dataset, iName, &iMesh);
      fail_unless(rc == FC_SUCCESS, "abort: failed to created mesh");
      rc = fc_setMeshCoords(iMesh, j, numVertPerType[i], iCoords);
      fail_unless(rc == FC_SUCCESS, "failed to set vertex coords");
      rc = fc_setMeshElementConns(iMesh, elemTypes[i], numElemPerType[i], 
                                  conns[i]);
      fail_unless(rc == FC_SUCCESS, "failed to set element conns");
      // create ied displacement variable
      rc = fc_createVariable(iMesh, "inverted displ", &iDisplVar);
      fail_unless(rc == FC_SUCCESS, "abort: failed to create variable");
      if (j == 1)
        mathType = FC_MT_SCALAR;
      else
        mathType = FC_MT_VECTOR;
      rc = fc_setVariableData(iDisplVar, numVertPerType[i], j, FC_AT_VERTEX, 
                              mathType, FC_DT_DOUBLE, (void**)iDispls);
      fail_unless(rc == FC_SUCCESS, "failed to set displVar data");

      // --- test get bounding box on whole mesh (& subset on FC_AT_WHOLE_MESH)

      // test fc_getMeshBoundingBox()
      rc = fc_getMeshBoundingBox(mesh, &temp_dim, &lowpoint, &highpoint);
      fail_unless(rc == FC_SUCCESS, "failed to get mesh bb");
      fail_unless(temp_dim == j, "mismatch of dim");
      for (k = 0; k < temp_dim; k++) {
        fail_unless(lowpoint[k] == comp_coords[k][0], 
                    "mismatch of lowpoint");
	//fail_unless(
        if (k > i || (k == i && i > 0) ) // cases where bb has no depth 
          fail_unless(highpoint[k] == lowpoint[k], "mismatch of highpoint"); 
        else
          fail_unless(highpoint[k] == comp_coords[k][dims[k]-1],
                      "mismatch of highpoint");
      }
      rc = fc_getMeshBoundingBox(iMesh, &temp_dim, &iLowpoint, &iHighpoint);
      fail_unless(rc == FC_SUCCESS, "failed to get imesh bb");
      fail_unless(temp_dim == j, "mismatch of idim");
      for (k = 0; k < temp_dim; k++) {
	fail_unless(iLowpoint[k] == -1.*highpoint[k], "mismatch of iLowpoint");
	fail_unless(iHighpoint[k] == -1.*lowpoint[k], "mismatch of iHighpoint");
      }

      // test fc_getDisplacedMeshBoundingBox()
      rc = fc_getDisplacedMeshBoundingBox(mesh, displVar, &temp_dim, &lowpoint,
                                          &highpoint); 
      fail_unless(rc == FC_SUCCESS, "failed to get displaced mesh bb");
      fail_unless(temp_dim == j, "mismatch of dim");
      for (k = 0; k < temp_dim; k++) {
        fail_unless(FC_DBL_EQUIV(lowpoint[k], comp_coords[k][0] +
                                 comp_displs[k][0]), "msimatch of lowpoint");
        if (k > i || (k == i && i > 0) ) // cases where bb has no depth 
          fail_unless(highpoint[k] == lowpoint[k], "mismatch of highpoint"); 
        else
          fail_unless(FC_DBL_EQUIV(highpoint[k],comp_coords[k][dims[k]-1] + 
                      comp_displs[k][dims[k]-1]), "mismatch of highpoint");
      }      
      rc = fc_getDisplacedMeshBoundingBox(iMesh, iDisplVar, &temp_dim,
					  &iLowpoint, &iHighpoint);
      fail_unless(rc == FC_SUCCESS, "failed to get displaced imesh bb");
      fail_unless(temp_dim == j, "mismatch of idim");
      for (k = 0; k < temp_dim; k++) {
	fail_unless(iLowpoint[k] == -1.*highpoint[k], "mismatch of iLowpoint");
	fail_unless(iHighpoint[k] == -1.*lowpoint[k], "mismatch of iHighpoint");
      }

      // setup for FC_AT_WHOLE_MESH subset tests
      rc = fc_createSubset(mesh, "temp subset", FC_AT_WHOLE_MESH, &subset);
      fail_unless(rc == FC_SUCCESS, "abort: failed to create subset");
      fc_addMemberToSubset(subset, 0);
      rc = fc_createSubset(iMesh, "inverted temp subset", FC_AT_WHOLE_MESH, &iSubset);
      fail_unless(rc == FC_SUCCESS, "abort: failed to create subset");
      fc_addMemberToSubset(iSubset, 0);

      // test fc_getSubsetBoundingBox() if FC_AT_WHOLE_MESH
      rc = fc_getSubsetBoundingBox(subset, &temp_dim, &lowpoint, &highpoint);
      fail_unless(rc == FC_SUCCESS, "failed to get subset bb");
      fail_unless(temp_dim == j, "mismatch of dim");
      for (k = 0; k < temp_dim; k++) {
        fail_unless(lowpoint[k] == comp_coords[k][0], "mismatch of lowpoint");
        if (k > i || (k == i && i > 0) ) // cases where bb has no depth 
          fail_unless(highpoint[k] == lowpoint[k], "mismatch of highpoint"); 
        else
          fail_unless(highpoint[k] == comp_coords[k][dims[k]-1],
                      "mismatch of highpoint");
      }
      rc = fc_getSubsetBoundingBox(iSubset, &temp_dim, &iLowpoint, &iHighpoint);
      fail_unless(rc == FC_SUCCESS, "failed to get isubset bb");
      fail_unless(temp_dim == j, "mismatch of idim");
      for (k = 0; k < temp_dim; k++) {
	fail_unless(iLowpoint[k] == -1.*highpoint[k], "mismatch of iLowpoint");
	fail_unless(iHighpoint[k] == -1.*lowpoint[k], "mismatch of iHighpoint");
      }

      // test fc_getDisplacedSubsetBoundingBox() if FC_AT_WHOLE_MESH
      rc = fc_getDisplacedSubsetBoundingBox(subset, displVar, &temp_dim,
                                            &lowpoint, &highpoint);
      fail_unless(rc == FC_SUCCESS, "failed to get displaced subset bb");
      fail_unless(temp_dim == j, "mismatch of dim");
      for (k = 0; k < temp_dim; k++) {
        //printf("    %d: %f %f\n", k, highpoint[k],
        //     comp_coords[k][dims[k]-1]+comp_displs[k][dims[k]-1]);
        fail_unless(FC_DBL_EQUIV(lowpoint[k], comp_coords[k][0] +
                                 comp_displs[k][0]), "mismatch of lowpoint");
        if (k > i || (k == i && i > 0) ) // cases where bb has no depth 
          fail_unless(highpoint[k] == lowpoint[k], "mismatch of highpoint"); 
        else
          fail_unless(FC_DBL_EQUIV(highpoint[k], comp_coords[k][dims[k]-1] + 
                      comp_displs[k][dims[k]-1]), "mismatch of highpoint");
      }
      rc = fc_getDisplacedSubsetBoundingBox(iSubset, iDisplVar, &temp_dim,
					    &iLowpoint, &iHighpoint);
      fail_unless(rc == FC_SUCCESS, "failed to get displaced isubset bb");
      fail_unless(temp_dim == j, "mismatch of idim");
      for (k = 0; k < temp_dim; k++) {
	fail_unless(iLowpoint[k] == -1.*highpoint[k], "mismatch of iLowpoint");
	fail_unless(iHighpoint[k] == -1.*lowpoint[k], "mismatch of iHighpoint");
      }

      // cleanup after FC_AT_WHOLE_MESH subset tests
      fc_deleteSubset(subset);
      fc_deleteSubset(iSubset);

      // --- now test different subset types

      // loop over assocs except FC_AT_WHOLE_MESH (already tested above)
      for (k = 0; k < numAssoc; k++) {
        int numMember, *memberIDs, elemID = numElemPerType[i] - 1;
        // skip this association type if it doesn't exist on the mesh
        fc_getMeshNumEntity(mesh, assocs[k], &numEntity);
        if (numEntity < 1)
          continue;

        // create subset - entities that make up last element
        rc = fc_createSubset(mesh, "temp subset", assocs[k], &subset);
        fail_unless(rc == FC_SUCCESS, "abort: failed to create subset");
        fc_changeMeshEntityType(mesh, FC_AT_ELEMENT, 1, &elemID,
                                assocs[k], 1, &numMember, &memberIDs);
        rc = fc_addArrayMembersToSubset(subset, numMember, memberIDs);
        fail_unless(rc == FC_SUCCESS, "abort: failed to add subset members");
	// inverted subset - use the same ids
	rc = fc_createSubset(iMesh, "inverted temp subset", assocs[k],
			     &iSubset);
	fail_unless(rc == FC_SUCCESS, "abort: failed to create subset");
        rc = fc_addArrayMembersToSubset(iSubset, numMember, memberIDs);
        fail_unless(rc == FC_SUCCESS, "abort: failed to add subset members");
	free(memberIDs);

        // test fc_getSubsetBoundingBox()
        rc = fc_getSubsetBoundingBox(subset, &temp_dim, &lowpoint, &highpoint);
        fail_unless(rc == FC_SUCCESS, "failed to get subset bb");
        fail_unless(temp_dim == j, "mismatch of dim");
        for (m = 0; m < temp_dim; m++) {
          if (i == 0 && m == 0) 
            fail_unless(lowpoint[m] == comp_coords[m][2],
                        "mismatch of lowpoint");
          else if (m == 0)
            fail_unless(lowpoint[m] == comp_coords[m][1],
                        "mismatch of lowpoint");
          else
            fail_unless(lowpoint[m] == comp_coords[m][0], 
                        "mismatch of lowpoint");
          if (m > i || (m == i && i > 0) ) // cases where bb has no depth 
            fail_unless(highpoint[m] == lowpoint[m], "mismatch of highpoint"); 
          else
            fail_unless(highpoint[m] == comp_coords[m][dims[m]-1],
                        "mismatch of highpoint");
        }
        rc = fc_getSubsetBoundingBox(iSubset, &temp_dim, &iLowpoint, &iHighpoint);
        fail_unless(rc == FC_SUCCESS, "failed to get isubset bb");
        fail_unless(temp_dim == j, "mismatch of idim");
	for (m = 0; m < temp_dim; m++) {
	  fail_unless(iLowpoint[m] == -1.*highpoint[m], "mismatch of iLowpoint");
	  fail_unless(iHighpoint[m] == -1.*lowpoint[m], "mismatch of iHighpoint");
	}

        // test fc_getDisplacedSubsetBoundingBox()
        rc = fc_getDisplacedSubsetBoundingBox(subset, displVar, &temp_dim,
                                              &lowpoint,  &highpoint);
        fail_unless(rc == FC_SUCCESS, "failed to get displaced subset bb");
        fail_unless(temp_dim == j, "mismatch of dim");
        for (m = 0; m < temp_dim; m++) {
          if (i == 0 && m == 0) 
            fail_unless(FC_DBL_EQUIV(lowpoint[m], comp_coords[m][2] +
                        comp_displs[m][2]), "mismatch of lowpoint");
          else if (m == 0)
            fail_unless(FC_DBL_EQUIV(lowpoint[m], comp_coords[m][1] + 
                        comp_displs[m][1]), "mismatch of lowpoint");
          else
            fail_unless(FC_DBL_EQUIV(lowpoint[m], comp_coords[m][0] + 
                        comp_displs[m][0]), "mismatch of lowpoint");
          if (m > i || (m == i && i > 0) ) // cases where bb has no depth 
            fail_unless(highpoint[m] == lowpoint[m], "mismatch of highpoint"); 
          else
            fail_unless(FC_DBL_EQUIV(highpoint[m], comp_coords[m][dims[m]-1] +
                        comp_displs[m][dims[m]-1]), "mismatch of highpoint");
        }
        rc = fc_getDisplacedSubsetBoundingBox(iSubset, iDisplVar, &temp_dim, 
					      &iLowpoint, &iHighpoint);
        fail_unless(rc == FC_SUCCESS, "failed to get displaced isubset bb");
        fail_unless(temp_dim == j, "mismatch of idim");
	for (m = 0; m < temp_dim; m++) {
	  fail_unless(iLowpoint[m] == -1.*highpoint[m], "mismatch of iLowpoint");
	  fail_unless(iHighpoint[m] == -1.*lowpoint[m], "mismatch of iHighpoint");
	}

        // --- test error conditions

        // fc_getSubsetBoundingBox() -- bad args
        rc = fc_getSubsetBoundingBox(badSubset, &temp_dim, &lowpoint, 
                                     &highpoint);
        fail_unless(rc != FC_SUCCESS, "should fail with bad subset");
        fail_unless(temp_dim == -1, "failure should return null");
        fail_unless(!memcmp(lowpoint, nullPoint1, 3*sizeof(double)),
                    "failure should return null");
        fail_unless(!memcmp(highpoint, nullPoint2, 3*sizeof(double)),
                    "failure should return null");
        rc = fc_getSubsetBoundingBox(subset, NULL, &lowpoint, &highpoint);
        fail_unless(rc != FC_SUCCESS, "should fail with no dim");
        fail_unless(!memcmp(lowpoint, nullPoint1, 3*sizeof(double)),
                    "failure should return null");
        fail_unless(!memcmp(highpoint, nullPoint2, 3*sizeof(double)),
                    "failure should return null");
        rc = fc_getSubsetBoundingBox(subset, &temp_dim, NULL, &highpoint);
        fail_unless(rc != FC_SUCCESS, "should fail with no lowpoint");
        fail_unless(temp_dim == -1, "failure should return null");
        fail_unless(!memcmp(highpoint, nullPoint2, 3*sizeof(double)),
                    "failure should return null");
        rc = fc_getSubsetBoundingBox(subset, &temp_dim, &lowpoint, NULL);
        fail_unless(rc != FC_SUCCESS, "should fail with no highpoint");
        fail_unless(temp_dim == -1, "failure should return null");
        fail_unless(!memcmp(lowpoint, nullPoint1, 3*sizeof(double)),
                    "failure should return null");

        // fc_getDisplacedSubsetBoundingBox() -- bad args
        rc = fc_getDisplacedSubsetBoundingBox(badSubset, displVar, &temp_dim,
                                              &lowpoint,  &highpoint);
        fail_unless(rc != FC_SUCCESS, "should fail with bad subset");
        fail_unless(temp_dim == -1, "failure should return null");
        fail_unless(!memcmp(lowpoint, nullPoint1, 3*sizeof(double)),
                    "failure should return null");
        fail_unless(!memcmp(highpoint, nullPoint2, 3*sizeof(double)),
                    "failure should return null");
        rc = fc_getDisplacedSubsetBoundingBox(subset, badVar, &temp_dim,
                                              &lowpoint,  &highpoint);
        fail_unless(rc != FC_SUCCESS, "should fail with bad var");
        fail_unless(temp_dim == -1, "failure should return null");
        fail_unless(!memcmp(lowpoint, nullPoint1, 3*sizeof(double)),
                    "failure should return null");
        fail_unless(!memcmp(highpoint, nullPoint2, 3*sizeof(double)),
                    "failure should return null");
        rc = fc_getDisplacedSubsetBoundingBox(subset, displVar, NULL,
                                              &lowpoint, &highpoint);
        fail_unless(rc != FC_SUCCESS, "should fail with no dim");
        fail_unless(!memcmp(lowpoint, nullPoint1, 3*sizeof(double)),
                    "failure should return null");
        fail_unless(!memcmp(highpoint, nullPoint2, 3*sizeof(double)),
                    "failure should return null");
        rc = fc_getDisplacedSubsetBoundingBox(subset, displVar, &temp_dim,
                                              NULL, &highpoint);
        fail_unless(rc != FC_SUCCESS, "should fail with no lowpoint");
        fail_unless(temp_dim == -1, "failure should return null");
        fail_unless(!memcmp(highpoint, nullPoint2, 3*sizeof(double)),
                    "failure should return null");
        rc = fc_getDisplacedSubsetBoundingBox(subset, displVar, &temp_dim,
                                              &lowpoint, NULL);
        fail_unless(rc != FC_SUCCESS, "should fail with no highpoint");
        fail_unless(temp_dim == -1, "failure should return null");
        fail_unless(!memcmp(lowpoint, nullPoint1, 3*sizeof(double)),
                    "failure should return null");

        // replace subset with empty subset for more testing
        fc_deleteSubset(subset);
        rc = fc_createSubset(mesh, "temp subset", assocs[k], &subset);
        fail_unless(rc == FC_SUCCESS, "abort: failed to create subset");

	// both getSubsetBoundingBox routines should fail if subset is empty
        rc = fc_getSubsetBoundingBox(subset, &temp_dim, &lowpoint, 
                                     &highpoint);
        fail_unless(rc != FC_SUCCESS, "should fail with empty subset");
        fail_unless(temp_dim == -1, "failure should return null");
        fail_unless(!memcmp(lowpoint, nullPoint1, 3*sizeof(double)),
                    "failure should return null");
        fail_unless(!memcmp(highpoint, nullPoint2, 3*sizeof(double)),
                    "failure should return null");
        rc = fc_getDisplacedSubsetBoundingBox(subset, displVar, &temp_dim,
                                              &lowpoint,  &highpoint);
        fail_unless(rc != FC_SUCCESS, "should fail with empty subset");
        fail_unless(temp_dim == -1, "failure should return null");
        fail_unless(!memcmp(lowpoint, nullPoint1, 3*sizeof(double)),
                    "failure should return null");
        fail_unless(!memcmp(highpoint, nullPoint2, 3*sizeof(double)),
                    "failure should return null");

        // cleanup
        fc_deleteSubset(subset);

      }

      // --- test error conditions 

      // fc_getMeshBoundingBox() -- bad args
      rc = fc_getMeshBoundingBox(badMesh, &temp_dim, &lowpoint, &highpoint);
      fail_unless(rc != FC_SUCCESS, "should fail with bad mesh");
      fail_unless(temp_dim == -1, "failure should return null");
      fail_unless(!memcmp(lowpoint, nullPoint1, 3*sizeof(double)),
                  "failure should return null");
      fail_unless(!memcmp(highpoint, nullPoint2, 3*sizeof(double)),
                  "failure should return null");
      rc = fc_getMeshBoundingBox(mesh, NULL, &lowpoint, &highpoint);
      fail_unless(rc != FC_SUCCESS, "should fail with no dim");
      fail_unless(!memcmp(lowpoint, nullPoint1, 3*sizeof(double)),
                  "failure should return null");
      fail_unless(!memcmp(highpoint, nullPoint2, 3*sizeof(double)),
                  "failure should return null");
      rc = fc_getMeshBoundingBox(mesh, &temp_dim, NULL, &highpoint);
      fail_unless(rc != FC_SUCCESS, "should fail with no lowpoint");
      fail_unless(temp_dim == -1, "failure should return null");
      fail_unless(!memcmp(highpoint, nullPoint2, 3*sizeof(double)),
                  "failure should return null");
      rc = fc_getMeshBoundingBox(mesh, &temp_dim, &lowpoint, NULL);
      fail_unless(rc != FC_SUCCESS, "should fail with no highpoint");
      fail_unless(temp_dim == -1, "failure should return null");
      fail_unless(!memcmp(lowpoint, nullPoint1, 3*sizeof(double)),
                  "failure should return null");

      // fc_getDisplacedMeshBoundingBox() -- bad args
      rc = fc_getDisplacedMeshBoundingBox(badMesh, displVar, &temp_dim, 
                                          &lowpoint, &highpoint);
      fail_unless(rc != FC_SUCCESS, "should fail with bad mesh");
      fail_unless(temp_dim == -1, "failure should return null");
      fail_unless(!memcmp(lowpoint, nullPoint1, 3*sizeof(double)),
                  "failure should return null");
      fail_unless(!memcmp(highpoint, nullPoint2, 3*sizeof(double)),
                  "failure should return null");
      rc = fc_getDisplacedMeshBoundingBox(mesh, badVar, &temp_dim, 
                                          &lowpoint, &highpoint);
      fail_unless(rc != FC_SUCCESS, "should fail with bad var");
      fail_unless(temp_dim == -1, "failure should return null");
      fail_unless(!memcmp(lowpoint, nullPoint1, 3*sizeof(double)),
                  "failure should return null");
      fail_unless(!memcmp(highpoint, nullPoint2, 3*sizeof(double)),
                  "failure should return null");
      rc = fc_getDisplacedMeshBoundingBox(mesh, displVar, NULL, &lowpoint, 
					  &highpoint);
      fail_unless(rc != FC_SUCCESS, "should fail with no dim");
      fail_unless(!memcmp(lowpoint, nullPoint1, 3*sizeof(double)),
                  "failure should return null");
      fail_unless(!memcmp(highpoint, nullPoint2, 3*sizeof(double)),
                  "failure should return null");
      rc = fc_getDisplacedMeshBoundingBox(mesh, displVar, &temp_dim, NULL, 
					  &highpoint);
      fail_unless(rc != FC_SUCCESS, "should fail with no lowpoint");
      fail_unless(temp_dim == -1, "failure should return null");
      fail_unless(!memcmp(highpoint, nullPoint2, 3*sizeof(double)),
                  "failure should return null");
      rc = fc_getDisplacedMeshBoundingBox(mesh, displVar, &temp_dim, &lowpoint,
					  NULL);
      fail_unless(rc != FC_SUCCESS, "should fail with no highpoint");
      fail_unless(temp_dim == -1, "failure should return null");
      fail_unless(!memcmp(lowpoint, nullPoint1, 3*sizeof(double)),
                  "failure should return null");  
    }

    fc_deleteMesh(mesh);
  }

  // cleanup
  rc = fc_deleteDataset(dataset);
  fail_unless(rc == FC_SUCCESS, "failed to delete subset at end of test");
}
END_TEST

// making a mesh from a bounding box
START_TEST(bb_mesh)
{
  FC_ReturnCode rc;
  int i;
  char name[20] = "blue berry", *temp_name;
  FC_Dataset dataset, badDataset = { 999, 999 };
  FC_Mesh mesh;
  FC_Coords low = { 1.0, 1.1, 1.2 };
  FC_Coords high = { 9.9, 9.7, 9.5 };
  int numVerts[3] = { 2, 4, 8 };
  int numElems[3] = { 1, 4, 12 };
  int temp_topodim, temp_dim, temp_numVert, temp_numElem;
  FC_ElementType temp_elemtype;
  double *coords_p;
  int *conns_p;

  // make dataset
  rc = fc_createDataset("temp dataset", &dataset);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create dataset");

  // successfully make bounding box mesh
  for (i = 1; i <= 3; i++) {
    // do it
    rc = fc_createBoundingBoxMesh(dataset, name, i, low, high, &mesh);
    fail_unless(rc == FC_SUCCESS, "failed to create mesh");
    fail_unless(fc_isMeshValid(mesh), "did not get a valid mesh");

    // check meta data
    fc_getMeshName(mesh, &temp_name);
    fail_unless(!strcmp(temp_name, name), "mismatch of name");
    free(temp_name);
    fc_getMeshInfo(mesh, &temp_topodim, &temp_dim, &temp_numVert,
		   &temp_numElem, &temp_elemtype);
    fail_unless(temp_topodim == 1, "mismatch of topodim");
    fail_unless(temp_dim == i, "mismatch of dim");
    fail_unless(temp_numVert = numVerts[i-1], "mismatch of numVertex");
    fail_unless(temp_numElem = numElems[i-1], "mismatch of numElement");
    fail_unless(temp_elemtype == FC_ET_LINE, "mismatch of element type");

    // check big data
    fc_getMeshCoordsPtr(mesh, &coords_p);
    fc_getMeshElementConnsPtr(mesh, &conns_p);

    if (i == 1) {
      fail_unless(conns_p[0] == 0 && conns_p[1] == 1, "mismatch of 1D conns");
      fail_unless(coords_p[0] == low[0], "mismatch of 1D coords");
      fail_unless(coords_p[1] == high[0], "mismatch of 1D coords");
    }
    else if (i == 2) {
      // vertices XY XY XY XY counter-clockwise around quad
      fail_unless(coords_p[0] == low[0], "mismatch of 2D coords");
      fail_unless(coords_p[1] == low[1], "mismatch of 2D coords");
      fail_unless(coords_p[2] == high[0], "mismatch of 2D coords");
      fail_unless(coords_p[3] == low[1], "mismatch of 2D coords");
      fail_unless(coords_p[4] == high[0], "mismatch of 2D coords");
      fail_unless(coords_p[5] == high[1], "mismatch of 2D coords");
      fail_unless(coords_p[6] == low[0], "mismatch of 2D coords");
      fail_unless(coords_p[7] == high[1], "mismatch of 2D coords");
      // Edges - counter clockwise around quad
      fail_unless(conns_p[0] == 0 && conns_p[1] == 1, "mismatch of 2D conns");
      fail_unless(conns_p[2] == 1 && conns_p[3] == 2, "mismatch of 2D conns");
      fail_unless(conns_p[4] == 2 && conns_p[5] == 3, "mismatch of 2D conns");
      fail_unless(conns_p[6] == 3 && conns_p[7] == 0, "mismatch of 2D conns");
    }
    else if (i == 3) {
      // vertices - XYZ XYZ XYZ XYZ counter-clockwise in XY plan at low[2] (Z)
      fail_unless(coords_p[0] == low[0], "mismatch of 3D coords");
      fail_unless(coords_p[1] == low[1], "mismatch of 3D coords");
      fail_unless(coords_p[2] == low[2], "mismatch of 3D coords");
      fail_unless(coords_p[3] == high[0], "mismatch of 3D coords");
      fail_unless(coords_p[4] == low[1], "mismatch of 3D coords");
      fail_unless(coords_p[5] == low[2], "mismatch of 3D coords");
      fail_unless(coords_p[6] == high[0], "mismatch of 3D coords");
      fail_unless(coords_p[7] == high[1], "mismatch of 3D coords");
      fail_unless(coords_p[8] == low[2], "mismatch of 3D coords");
      fail_unless(coords_p[9] == low[0], "mismatch of 3D coords");
      fail_unless(coords_p[10] == high[1], "mismatch of 3D coords");
      fail_unless(coords_p[11] == low[2], "mismatch of 3D coords");
      // vertices - XYZ XYZ XYZ XYZ counter-clockwise in XY plan at high[2] (Z)
      fail_unless(coords_p[12] == low[0], "mismatch of 3D coords");
      fail_unless(coords_p[13] == low[1], "mismatch of 3D coords");
      fail_unless(coords_p[14] == high[2], "mismatch of 3D coords");
      fail_unless(coords_p[15] == high[0], "mismatch of 3D coords");
      fail_unless(coords_p[16] == low[1], "mismatch of 3D coords");
      fail_unless(coords_p[17] == high[2], "mismatch of 3D coords");
      fail_unless(coords_p[18] == high[0], "mismatch of 3D coords");
      fail_unless(coords_p[19] == high[1], "mismatch of 3D coords");
      fail_unless(coords_p[20] == high[2], "mismatch of 3D coords");
      fail_unless(coords_p[21] == low[0], "mismatch of 3D coords");
      fail_unless(coords_p[22] == high[1], "mismatch of 3D coords");
      fail_unless(coords_p[23] == high[2], "mismatch of 3D coords");
      // edges - walk around bottom quad
      fail_unless(conns_p[0] == 0 && conns_p[1] == 1, "mismatch of 3D conns");
      fail_unless(conns_p[2] == 1 && conns_p[3] == 2, "mismatch of 3D conns");
      fail_unless(conns_p[4] == 2 && conns_p[5] == 3, "mismatch of 3D conns");
      fail_unless(conns_p[6] == 3 && conns_p[7] == 0, "mismatch of 3D conns");
      // edges - walk around top quad
      fail_unless(conns_p[8] == 4 && conns_p[9] == 5, "mismatch of 3D conns");
      fail_unless(conns_p[10] == 5 && conns_p[11] == 6, "mismatch of 3D conns");
      fail_unless(conns_p[12] == 6 && conns_p[13] == 7, "mismatch of 3D conns");
      fail_unless(conns_p[14] == 7 && conns_p[15] == 4, "mismatch of 3D conns");
      // edges - connecting bottom to top quad
      fail_unless(conns_p[16] == 0 && conns_p[17] == 4, "mismatch of 3D conns");
      fail_unless(conns_p[18] == 1 && conns_p[19] == 5, "mismatch of 3D conns");
      fail_unless(conns_p[20] == 2 && conns_p[21] == 6, "mismatch of 3D conns");
      fail_unless(conns_p[22] == 3 && conns_p[23] == 7, "mismatch of 3D conns");
     }
  }
   
  // bad args to function
  rc = fc_createBoundingBoxMesh(badDataset, name, 3, low, high, &mesh);
  fail_unless(rc != FC_SUCCESS, "should fail if bad dataset");
  fail_unless(FC_HANDLE_EQUIV(mesh, FC_NULL_MESH), "fail should return null");
  rc = fc_createBoundingBoxMesh(dataset, NULL, 3, low, high, &mesh);
  fail_unless(rc != FC_SUCCESS, "should fail if null name");
  fail_unless(FC_HANDLE_EQUIV(mesh, FC_NULL_MESH), "fail should return null");
  rc = fc_createBoundingBoxMesh(dataset, name, -1, low, high, &mesh);
  fail_unless(rc != FC_SUCCESS, "should fail if bad dim");
  fail_unless(FC_HANDLE_EQUIV(mesh, FC_NULL_MESH), "fail should return null");
  rc = fc_createBoundingBoxMesh(dataset, name, 4, low, high, &mesh);
  fail_unless(rc != FC_SUCCESS, "should fail if bad dim");
  fail_unless(FC_HANDLE_EQUIV(mesh, FC_NULL_MESH), "fail should return null");
  rc = fc_createBoundingBoxMesh(dataset, name, 3, NULL, high, &mesh);
  fail_unless(rc != FC_SUCCESS, "should fail if null lowpoint");
  fail_unless(FC_HANDLE_EQUIV(mesh, FC_NULL_MESH), "fail should return null");
  rc = fc_createBoundingBoxMesh(dataset, name, 3, low, NULL, &mesh);
  fail_unless(rc != FC_SUCCESS, "should fail if null highpoint");
  fail_unless(FC_HANDLE_EQUIV(mesh, FC_NULL_MESH), "fail should return null");
  rc = fc_createBoundingBoxMesh(dataset, name, 3, low, high, NULL);
  fail_unless(rc != FC_SUCCESS, "should fail if null mesh");
  fail_unless(FC_HANDLE_EQUIV(mesh, FC_NULL_MESH), "fail should return null");

  // cleanup
  rc = fc_deleteDataset(dataset);
  fail_unless(rc == FC_SUCCESS, "failed to delete subset at end of test");
}
END_TEST

START_TEST(bb_projection) {
  FC_ReturnCode rc;
  FC_Dataset dataset;
  FC_Mesh mesh, otherMesh, badMesh = {999,999};
  FC_Subset compSubset, refSubset1, refSubset2, faceSubset, intersectSubset,
    otherSubset, emptySubset, badSubset = {999,999};
  FC_Coords lowers = { 0., 0., 0. };
  FC_Coords uppers = { 1., 1., 1. };
  FC_Coords proj_normal;
  int elem_ID1, elem_ID2, elem_ID3;
  int intersectval,i,j,k;

  //known ids
  elem_ID1 = 0;
  elem_ID2 = 9;
  elem_ID3 = 10;

  rc = fc_createDataset("temp.xxx", &dataset);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create dataset");

  rc = fc_createSimpleHexMesh(dataset, "junk", 3, 3, 3, lowers, uppers,
                              &mesh);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create simple hex mesh");

  rc = fc_createSubset(mesh,"empty",FC_AT_ELEMENT,&emptySubset);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create subset");

  rc = fc_createSubset(mesh,"comp",FC_AT_ELEMENT,&compSubset);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create subset");

  rc = fc_addMemberToSubset(compSubset,elem_ID1);
  fail_unless(rc == FC_SUCCESS, "abort: failed to add member to subset");

  rc = fc_createSubset(mesh,"ref",FC_AT_ELEMENT,&refSubset1);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create subset");

  rc = fc_addMemberToSubset(refSubset1,elem_ID2);
  fail_unless(rc == FC_SUCCESS, "abort: failed to add member to subset");

  rc = fc_createSubset(mesh,"ref2",FC_AT_ELEMENT,&refSubset2);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create subset");

  rc = fc_addMemberToSubset(refSubset2,elem_ID3);
  fail_unless(rc == FC_SUCCESS, "abort: failed to add member to subset");

  rc = fc_createSimpleHexMesh(dataset, "junk", 3, 3, 3, lowers, uppers,
                              &otherMesh);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create simple hex mesh");

  rc = fc_createSubset(otherMesh,"comp",FC_AT_ELEMENT,&otherSubset);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create subset");

  rc = fc_addMemberToSubset(otherSubset,elem_ID1);
  fail_unless(rc == FC_SUCCESS, "abort: failed to add member to subset");

  for (i = 0; i < 3; i++){
    for (j = 0; j < 2; j++){
      int *elems;
      proj_normal[0] = 0;
      proj_normal[1] = 0;
      proj_normal[2] = 0;
      proj_normal[i] = 1;
      if (j == 1 ){
	proj_normal[i]*=-1;
      }

      //projected to itself
      rc = fc_projectedElementBoundingBoxIntersection(mesh, elem_ID1,
						      elem_ID1,
						      proj_normal,
						      &intersectval);
      fail_unless(rc == FC_SUCCESS, "should work for self");
      fail_unless(intersectval == 1, "wrong val for intersect");
     
       rc = fc_projectedSubsetElementBoundingBoxIntersection(compSubset, 
						     compSubset,
						     proj_normal,
						     &intersectSubset);
      fail_unless(rc == FC_SUCCESS, "should work for self");
      rc = fc_getSubsetMembersAsArray(intersectSubset,&k,&elems);
      fail_unless(rc == FC_SUCCESS, "can't get subset members");
      fail_unless(k == 1, "should be one member");
      fail_unless(elems[0] == elem_ID1 , "mismatch of subset member");
      free(elems);
      fc_deleteSubset(intersectSubset);

      //projected known stacked elements
      rc = fc_projectedElementBoundingBoxIntersection(mesh, elem_ID1,
						      elem_ID2,
						      proj_normal,
						      &intersectval);
      fail_unless(rc == FC_SUCCESS, "should work for stacked elems in z");
      if (proj_normal[2] != 0){
	fail_unless(intersectval == 1, "should intersect");
      }else{
	fail_unless(intersectval == 0, "should not intersect");
      }

      rc = fc_projectedSubsetElementBoundingBoxIntersection(compSubset, 
							    refSubset1,
							    proj_normal,
							    &intersectSubset);
      fail_unless(rc == FC_SUCCESS, "should work for stacked elems in z");
      rc = fc_getSubsetMembersAsArray(intersectSubset,&k,&elems);
      fail_unless(rc == FC_SUCCESS, "can't get subset members");
      if (proj_normal[2] != 0){
	fail_unless(k == 1, "should be one member");
	fail_unless(elems[0] == elem_ID2 , "mismatch of subset member");
	free(elems);
      }else{
	fail_unless(k == 0, "should be no members");
      }
      fc_deleteSubset(intersectSubset);

      // elems 1 off should never work along an axis(line intersection)
      rc = fc_projectedElementBoundingBoxIntersection(mesh, elem_ID1,
						      elem_ID3,
						      proj_normal,
						      &intersectval);
      fail_unless(rc == FC_SUCCESS, "should work for 1 off elems");
      fail_unless(intersectval == 0, "should not intersect");

      rc = fc_projectedSubsetElementBoundingBoxIntersection(compSubset, 
							    refSubset2,
							    proj_normal,
							    &intersectSubset);
      fail_unless(rc == FC_SUCCESS, "should work for 1 off elems");
      rc = fc_getSubsetMembersAsArray(intersectSubset,&k,&elems);
      fail_unless(rc == FC_SUCCESS, "can't get subset members");
      fail_unless(k == 0, "should be no members");
      fc_deleteSubset(intersectSubset);

      //and now make a diagonal
      proj_normal[0] = 1/sqrt(2.0);
      proj_normal[1] = 0;
      proj_normal[2] = 1/sqrt(2.0);

      if (j){
	proj_normal[0]*=-1;
	proj_normal[2]*=-1;
      }

      rc = fc_projectedElementBoundingBoxIntersection(mesh, elem_ID1,
						      elem_ID3,
						      proj_normal,
						      &intersectval);
      fail_unless(rc == FC_SUCCESS, "should work for diag elems in z");
      fail_unless(intersectval == 1, "should intersect");      


      rc = fc_projectedSubsetElementBoundingBoxIntersection(compSubset, 
							    refSubset2,
							    proj_normal,
							    &intersectSubset);
      fail_unless(rc == FC_SUCCESS, "should work for 1 off elems");
      rc = fc_getSubsetMembersAsArray(intersectSubset,&k,&elems);
      fail_unless(rc == FC_SUCCESS, "can't get subset members");
      fail_unless(k == 1, "should be one member");
      fail_unless(elems[0] == elem_ID3 , "mismatch of subset member");
      free(elems);
      fc_deleteSubset(intersectSubset);
    }
  }


  //empty compsubset
  rc = fc_projectedSubsetElementBoundingBoxIntersection(emptySubset,
							refSubset1,
							proj_normal,
							&intersectSubset);
  fail_unless(rc == FC_SUCCESS, "should work for empty subset");
  fail_unless(fc_isSubsetValid(intersectSubset), 
	      "should return a valid subset");
  rc = fc_getSubsetNumMember(intersectSubset,&i);
  fail_unless(rc == FC_SUCCESS, "can't get subset nummembers");
  fail_unless(i == 0, "should return empty subset");
  fc_deleteSubset(intersectSubset);

  //wil want something htat project down is ok, and proj at an angle is not

  //bad args
  rc = fc_createSubset(mesh,"face",FC_AT_FACE,&faceSubset);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create subset");

  rc = fc_addMemberToSubset(faceSubset,1);
  fail_unless(rc == FC_SUCCESS, "abort: failed to add member to subset");

  proj_normal[0] = 1;
  proj_normal[1] = 0;
  proj_normal[2] = 0;
  rc = fc_projectedElementBoundingBoxIntersection(FC_NULL_MESH, 1, 2,
						  proj_normal,
						  &intersectval);
  fail_unless(rc != FC_SUCCESS, "should fail for null mesh");
  fail_unless(intersectval == -1, "failure should return -1");

  rc = fc_projectedElementBoundingBoxIntersection(badMesh, 1, 2,
						  proj_normal,
						  &intersectval);
  fail_unless(rc != FC_SUCCESS, "should fail for bad mesh");
  fail_unless(intersectval == -1, "failure should return -1");


  rc = fc_projectedElementBoundingBoxIntersection(mesh, -1, 2,
						  proj_normal,
						  &intersectval);
  fail_unless(rc != FC_SUCCESS, "should fail for bad elemid");
  fail_unless(intersectval == -1, "failure should return -1");

  rc = fc_projectedElementBoundingBoxIntersection(mesh, 100, 2,
						  proj_normal,
						  &intersectval);
  fail_unless(rc != FC_SUCCESS, "should fail for bad elemid");
  fail_unless(intersectval == -1, "failure should return -1");

  rc = fc_projectedElementBoundingBoxIntersection(mesh, 1, -1,
						  proj_normal,
						  &intersectval);
  fail_unless(rc != FC_SUCCESS, "should fail for bad elemid");
  fail_unless(intersectval == -1, "failure should return -1");

  rc = fc_projectedElementBoundingBoxIntersection(mesh, 1, 100,
						  proj_normal,
						  &intersectval);
  fail_unless(rc != FC_SUCCESS, "should fail for bad elemid");
  fail_unless(intersectval == -1, "failure should return -1");

  proj_normal[2] = 2;
  rc = fc_projectedElementBoundingBoxIntersection(mesh, 1, 2,
						  proj_normal,
						  &intersectval);
  fail_unless(rc != FC_SUCCESS, "should fail for normal bad magnitude");
  fail_unless(intersectval == -1, "failure should return -1");
  proj_normal[2] = 0;   //fix it back

  rc = fc_projectedElementBoundingBoxIntersection(mesh, 1, 100,
						  proj_normal,
						  NULL);
  fail_unless(rc != FC_SUCCESS, "should fail for NULL return arg");


  rc = fc_projectedSubsetElementBoundingBoxIntersection(FC_NULL_SUBSET,
							refSubset1,
							proj_normal,
							&intersectSubset);
  fail_unless(rc != FC_SUCCESS, "should fail for NULL subset");
  fail_unless(FC_HANDLE_EQUIV(intersectSubset,FC_NULL_SUBSET),
	      "failure should return NULL");

  rc = fc_projectedSubsetElementBoundingBoxIntersection(badSubset,
							refSubset1,
							proj_normal,
							&intersectSubset);
  fail_unless(rc != FC_SUCCESS, "should fail for bad subset");
  fail_unless(FC_HANDLE_EQUIV(intersectSubset,FC_NULL_SUBSET),
	      "failure should return NULL");

  rc = fc_projectedSubsetElementBoundingBoxIntersection(compSubset,
							badSubset,
							proj_normal,
							&intersectSubset);
  fail_unless(rc != FC_SUCCESS, "should fail for bad subset");
  fail_unless(FC_HANDLE_EQUIV(intersectSubset,FC_NULL_SUBSET),
	      "failure should return NULL");

  rc = fc_projectedSubsetElementBoundingBoxIntersection(compSubset,
							FC_NULL_SUBSET,
							proj_normal,
							&intersectSubset);
  fail_unless(rc != FC_SUCCESS, "should fail for NULL subset");
  fail_unless(FC_HANDLE_EQUIV(intersectSubset,FC_NULL_SUBSET),
	      "failure should return NULL");

  rc = fc_projectedSubsetElementBoundingBoxIntersection(compSubset,
							emptySubset,
							proj_normal,
							&intersectSubset);
  fail_unless(rc != FC_SUCCESS, "should fail for empty subset");
  fail_unless(FC_HANDLE_EQUIV(intersectSubset,FC_NULL_SUBSET),
	      "failure should return NULL");

  rc = fc_projectedSubsetElementBoundingBoxIntersection(faceSubset,
							refSubset1,
							proj_normal,
							&intersectSubset);
  fail_unless(rc != FC_SUCCESS, "should fail for wrong subset assoc");
  fail_unless(FC_HANDLE_EQUIV(intersectSubset,FC_NULL_SUBSET),
	      "failure should return NULL");

  rc = fc_projectedSubsetElementBoundingBoxIntersection(compSubset,
							otherSubset,
							proj_normal,
							&intersectSubset);
  fail_unless(rc != FC_SUCCESS, "should fail for subsets on different meshes");
  fail_unless(FC_HANDLE_EQUIV(intersectSubset,FC_NULL_SUBSET),
	      "failure should return NULL");

  proj_normal[2] = 2;
  rc = fc_projectedSubsetElementBoundingBoxIntersection(compSubset,
							refSubset1,
							proj_normal,
							&intersectSubset);
  fail_unless(rc != FC_SUCCESS, "should fail for normal bad magnitude");
  fail_unless(FC_HANDLE_EQUIV(intersectSubset,FC_NULL_SUBSET),
	      "failure should return NULL");
  proj_normal[2] = 0;   //fix it back

  rc = fc_projectedSubsetElementBoundingBoxIntersection(compSubset,
							refSubset1,
							proj_normal,
							NULL);
  fail_unless(rc != FC_SUCCESS, "should fail for NULL arg");

  fc_deleteSubset(refSubset1); 
  fc_deleteSubset(refSubset2);
  fc_deleteSubset(compSubset);
  fc_deleteSubset(faceSubset);
  fc_deleteSubset(otherSubset);
  fc_deleteMesh(mesh);
  fc_deleteMesh(otherMesh);
  fc_deleteDataset(dataset);
}
END_TEST

// test getting diameter
// basic layout of vertices is on two side-by-side hexes
// for 1D use only lower edge, for 2D use bottom faces
// test getting diameter of multiple meshes, subsets
// duplicate the mesh, but move 3*side_length in x direction
START_TEST(get_diameter) {
  FC_ReturnCode rc;
  int i, j, k, m;
  double side_length = 2.4;
  char name[20] = "blue berry", name2[20] = "banana frost";
  FC_Dataset dataset;
  FC_Mesh mesh, badMesh = { 999, 999 }, mesh2, meshes[2];
  FC_Subset subset, emptySubset, badSubset = { 999, 999 };
  FC_Subset subset2, subset3, subsets[2], subs_same_m[2];
  FC_Variable displVar, badVar = { 999, 999 }, displVar2;
  FC_Variable displVars[2], displs_same_m[2];
  FC_MathType mathType;
  int numElemType = 4;
  FC_ElementType elemTypes[4] = { FC_ET_POINT,  FC_ET_LINE,
                                  FC_ET_QUAD,   FC_ET_HEX };
  int dims[3] = { 3, 2, 2 };
  int minDimPerType[4] = { 1, 1, 2, 3 };
  int numVertPerType[4] = { 3, 3, 6, 12}, numElemPerType[4] = { 3, 2, 2, 2 };
  int conns[4][16] = { { 0,    1,    2 },
                       { 0, 1,    1, 2 },
                       { 0, 1, 4, 3,    1, 2, 5, 4 },
                       { 0, 1, 4, 3, 6, 7, 10, 9,    1, 2, 5, 4, 7, 8, 11, 10
                       } };
  double coords[3][3*12], coords2[3][3*12];
  double temp_diameter;
  int temp_numPair, *temp_pairIDs;
  double diameters[4] = { 2*side_length, 2*side_length, sqrt(5)*side_length, 
			  sqrt(6)*side_length }; 
  int numPairs[4] = { 1, 1, 2, 4 };
  int pairIDs[4][8] = { { 0, 2 }, { 0, 2 }, { 0, 5, 2, 3 },
			{ 0, 11, 2, 9, 3, 8, 5, 6 } };
  double subset_diameters[4] = { 0., side_length, sqrt(2)*side_length,
				  sqrt(3)*side_length } ;
  int subset_numPairs[4] = { 0, 1, 2, 4 };
  int subset_pairIDs[4][8] = { { }, { 1, 2 }, { 1, 5, 2, 4 }, 
			       { 1, 11, 2, 10, 4, 8, 5, 7 } };
  int numAssoc = 4, numEntity;
  FC_AssociationType assocs[4] = { FC_AT_VERTEX, FC_AT_EDGE, 
                                   FC_AT_FACE, FC_AT_ELEMENT };  
  FC_Mesh* temp_pairMeshes;
  FC_Subset* temp_pairSubsets;
  double meshes_diameters[4] = { 5*side_length, 5*side_length, 
				 sqrt(26)*side_length, sqrt(27)*side_length };
  // NOTE: meshes_pairIDs is the same as pairIDs, just slightly diff order
  int meshes_pairIDs[4][8] = { { 0, 2 }, { 0, 2 }, { 0, 5, 3, 2 },
			       { 0, 11, 3, 8, 6, 5, 9, 2 } };
  double subsets_diameters[4] = { 3*side_length, 4*side_length,
				  sqrt(17)*side_length, 
				  sqrt(18)*side_length };
  int subsets_numPairs[4] = { 1, 1, 2, 4 };
  int subsets_pairIDs[4][8] = { { 2, 2 }, { 1, 2 }, { 1, 5, 4, 2}, 
				{ 1, 11, 4, 8, 7, 5, 10, 2 } };
  // setup coords
  // 1D coords
  for (i = 0; i < dims[0]; i++) 
    coords[0][i] = i*side_length;
  // 2D coords
  for (j = 0; j < dims[1]; j++) {
    for (i = 0; i < dims[0]; i++) {
      int ID = j*dims[0] + i;
      coords[1][2*ID] = i*side_length;
      coords[1][2*ID + 1] = j*side_length;
    }
  }
  // 3D coords
  for (k = 0; k < dims[2]; k++) {
    for (j = 0; j < dims[1]; j++) {
      for (i = 0; i < dims[0]; i++) {
	int ID = k*dims[0]*dims[1] + j*dims[0] + i;
	coords[2][3*ID] = i*side_length;
	coords[2][3*ID+1] = j*side_length;
	coords[2][3*ID+2] = k*side_length;
      }
    }
  }
  // the 2nd mesh is translated 3*side_length in the x direction
  // (separation of side-length between meshes
  for (i = 0; i < 3; i++)
    for (j = 0; j < 3*12; j++)
      coords2[i][j] = coords[i][j];
  for (i = 0; i < dims[0]; i++)
    coords2[0][i] += 3*side_length;
  for (i = 0; i < dims[1]*dims[0]; i++)
    coords2[1][2*i] += 3*side_length;
  for (i = 0; i < dims[2]*dims[1]*dims[0]; i++)
    coords2[2][3*i] += 3*side_length;

  // setup - create test mesh
  rc = fc_createDataset("temp.xxx", &dataset);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create dataset");

  // test different types
  for (i = 0; i < numElemType; i++) { // loop over elem types
    for (j = 1; j <= 3; j++) { // loop over space dimensions
      // skip if it doesn't make sense to make this elem type at this dim
      if (minDimPerType[i] > j)
        continue;

      // --- create data structures for testing: mesh & displVar

      // create meshes
      rc = fc_createMesh(dataset, name, &mesh);
      fail_unless(rc == FC_SUCCESS, "abort: failed to create mesh");
      rc = fc_setMeshCoords(mesh, j, numVertPerType[i], coords[j-1]);
      fail_unless(rc == FC_SUCCESS, "failed to set vertex coords");
      rc = fc_setMeshElementConns(mesh, elemTypes[i], numElemPerType[i], 
                                  conns[i]);
      fail_unless(rc == FC_SUCCESS, "failed to set element conns");
      rc = fc_createMesh(dataset, name2, &mesh2);
      fail_unless(rc == FC_SUCCESS, "abort: failed to create 2nd mesh");
      rc = fc_setMeshCoords(mesh2, j, numVertPerType[i], coords2[j-1]);
      fail_unless(rc == FC_SUCCESS, "failed to set vertex coords");
      rc = fc_setMeshElementConns(mesh2, elemTypes[i], numElemPerType[i], 
                                  conns[i]);
      fail_unless(rc == FC_SUCCESS, "failed to set element conns");
      meshes[0] = mesh;
      meshes[1] = mesh2;
      // create displacement variables
      if (j == 1)
	mathType = FC_MT_SCALAR;
      else
	mathType = FC_MT_VECTOR;
      rc = fc_createVariable(mesh, "displ", &displVar);
      fail_unless(rc == FC_SUCCESS, "abort: failed to create displ variable");
      rc = fc_setVariableData(displVar, numVertPerType[i], j, FC_AT_VERTEX, 
			      mathType, FC_DT_DOUBLE, (void*)coords[j-1]);
      fail_unless(rc == FC_SUCCESS, "failed to set displVar data");
      rc = fc_createVariable(mesh2, "displ", &displVar2);
      fail_unless(rc == FC_SUCCESS, "abort: failed to create displ variable");
      rc = fc_setVariableData(displVar2, numVertPerType[i], j, FC_AT_VERTEX, 
			      mathType, FC_DT_DOUBLE, (void*)coords2[j-1]);
      fail_unless(rc == FC_SUCCESS, "failed to set displVar data");
      displVars[0] = displVar;
      displVars[1] = displVar2;
      // create empty subset
      rc = fc_createSubset(mesh, "empty subset", FC_AT_VERTEX, &emptySubset);
      fail_unless(rc == FC_SUCCESS, "abort: faile to create empty subset");
 
      // -- test get diameter on whole mesh(es) (& subset(s) on FC_AT_WHOLE_MESH)

      // fc_getMeshDiameter()
      rc = fc_getMeshDiameter(mesh, &temp_diameter, &temp_numPair, 
			      &temp_pairIDs);
      fail_unless(rc == FC_SUCCESS, "failed to get diameter of mesh");
      fail_unless(FC_DBL_EQUIV(temp_diameter, diameters[i]),
		  "mismatch of diameter");
      fail_unless(temp_numPair == numPairs[i], "mismatch of numPair");
      for (k = 0; k < 2*numPairs[i]; k++)
	fail_unless(temp_pairIDs[k] == pairIDs[i][k], "mismatch of pair ids");
      free(temp_pairIDs);
      // fc_getMeshesDiameter() on single mesh
      rc = fc_getMeshesDiameter(1, &mesh, &temp_diameter, &temp_numPair,
				&temp_pairIDs, &temp_pairMeshes);
      fail_unless(rc == FC_SUCCESS, "failed to get diameter of mesh");
      fail_unless(FC_DBL_EQUIV(temp_diameter, diameters[i]),
		  "mismatch of diameter");
      fail_unless(temp_numPair == numPairs[i], "mismatch of numPair");
      for (k = 0; k < 2*numPairs[i]; k++) {
	fail_unless(temp_pairIDs[k] == pairIDs[i][k], "mismatch of pair ids");
	fail_unless(FC_HANDLE_EQUIV(temp_pairMeshes[k], mesh),
		    "mismatch of pair meshes");
      }
      free(temp_pairIDs);
      free(temp_pairMeshes);
      // fc_getMeshesDiameter() on two meshes
      rc = fc_getMeshesDiameter(2, meshes, &temp_diameter, &temp_numPair,
				&temp_pairIDs, &temp_pairMeshes);
      fail_unless(rc == FC_SUCCESS, "failed to get diameter of mesh");
      fail_unless(FC_DBL_EQUIV(temp_diameter, meshes_diameters[i]),
		  "mismatch of diameter");
      fail_unless(temp_numPair == numPairs[i], "mismatch of numPair");
      for (k = 0; k < 2*numPairs[i]; k++) 
	fail_unless(temp_pairIDs[k] == meshes_pairIDs[i][k], 
		    "mismatch of pair ids");
      for (k = 0; k < numPairs[i]; k++) {
	fail_unless(FC_HANDLE_EQUIV(temp_pairMeshes[2*k], meshes[0]),
		    "mismatch of pair meshes");
	fail_unless(FC_HANDLE_EQUIV(temp_pairMeshes[2*k+1], meshes[1]),
		    "mismatch of pair meshes");
      }
      free(temp_pairIDs);
      free(temp_pairMeshes);

      // fc_getDisplacedMeshDiameter()
      rc = fc_getDisplacedMeshDiameter(mesh, displVar, &temp_diameter, 
				       &temp_numPair, &temp_pairIDs);
      fail_unless(rc == FC_SUCCESS, "failed to get diameter of displ mesh");
      fail_unless(FC_DBL_EQUIV(temp_diameter, 2*diameters[i]),
		  "mismatch of displ diameter");
      fail_unless(temp_numPair == numPairs[i], "mismatch of displ numPair");
      for (k = 0; k < 2*numPairs[i]; k++)
	fail_unless(temp_pairIDs[k] == pairIDs[i][k], 
		    "mismatch of displ pair ids");
      free(temp_pairIDs);
      // fc_getDisplacedMeshesDiameter() on single mesh
      rc = fc_getDisplacedMeshesDiameter(1, &mesh, &displVar, &temp_diameter, 
					 &temp_numPair, &temp_pairIDs,
					 &temp_pairMeshes);
      fail_unless(rc == FC_SUCCESS, "failed to get diameter of displ mesh");
      fail_unless(FC_DBL_EQUIV(temp_diameter, 2*diameters[i]),
		  "mismatch of displ diameter");
      fail_unless(temp_numPair == numPairs[i], "mismatch of displ numPair");
      for (k = 0; k < 2*numPairs[i]; k++) {
	fail_unless(temp_pairIDs[k] == pairIDs[i][k], 
		    "mismatch of displ pair ids");
	fail_unless(FC_HANDLE_EQUIV(temp_pairMeshes[k], mesh),
		    "mismatch of displ pair meshes");
      }
      free(temp_pairIDs);
      free(temp_pairMeshes);
      // fc_getDisplacedMeshesDiameter() on two meshes
      rc = fc_getDisplacedMeshesDiameter(2, meshes, displVars, &temp_diameter,
					 &temp_numPair,
					 &temp_pairIDs, &temp_pairMeshes);
      fail_unless(rc == FC_SUCCESS, "failed to get diameter of mesh");
      fail_unless(FC_DBL_EQUIV(temp_diameter, 2*meshes_diameters[i]),
		  "mismatch of diameter");
      fail_unless(temp_numPair == numPairs[i], "mismatch of numPair");
      for (k = 0; k < 2*numPairs[i]; k++) 
	fail_unless(temp_pairIDs[k] == meshes_pairIDs[i][k], 
		    "mismatch of pair ids");
      for (k = 0; k < numPairs[i]; k++) {
	fail_unless(FC_HANDLE_EQUIV(temp_pairMeshes[2*k], meshes[0]),
		    "mismatch of pair meshes");
	fail_unless(FC_HANDLE_EQUIV(temp_pairMeshes[2*k+1], meshes[1]),
		    "mismatch of pair meshes");
      }
      free(temp_pairIDs);
      free(temp_pairMeshes);

      // -- special cases : don't have to return pairs

      // fc_getMeshDiameter() w/o pairs is o.k.
      rc = fc_getMeshDiameter(mesh, &temp_diameter, NULL, NULL);
      fail_unless(rc == FC_SUCCESS, "failed to get diameter with NULL pairs");
      fail_unless(FC_DBL_EQUIV(temp_diameter, diameters[i]),
		 "mismatch of diameter");
      // fc_getMeshesDiameter() w/o pairs is o.k.
      rc = fc_getMeshesDiameter(1, &mesh, &temp_diameter, NULL, NULL, NULL);
      fail_unless(rc == FC_SUCCESS, "failed to get diameter with NULL pairs");
      fail_unless(FC_DBL_EQUIV(temp_diameter, diameters[i]),
		 "mismatch of diameter");
      // fc_getDisplacedMeshDiameter() w/o pairs is o.k.
      rc = fc_getDisplacedMeshDiameter(mesh, displVar, &temp_diameter, NULL, 
				       NULL);
      fail_unless(rc == FC_SUCCESS, 
		  "failed to get displ diameter with NULL pairs");
      fail_unless(FC_DBL_EQUIV(temp_diameter, 2*diameters[i]),
		  "mismatch of displ diameter");
      // fc_getDisplacedMeshesDiameter() w/o pairs is o.k.
      rc = fc_getDisplacedMeshesDiameter(1, &mesh, &displVar, &temp_diameter, 
					 NULL, NULL, NULL);
      fail_unless(rc == FC_SUCCESS, 
		  "failed to get displ diameter with NULL pairs");
      fail_unless(FC_DBL_EQUIV(temp_diameter, 2*diameters[i]),
		  "mismatch of displ diameter");

      // --- now test different subset types

      // -- Deal with FC_AT_WHOLE_MESH subset test differently than over subsets
      // -- tests below -- results are the same as for mesh stuff
 
      // setup for FC_AT_WHOLE_MESH subset tests
      rc = fc_createSubset(mesh, "temp subset", FC_AT_WHOLE_MESH, &subset);
      fail_unless(rc == FC_SUCCESS, "abort: failed to create subset");
      fc_addMemberToSubset(subset, 0);
      rc = fc_createSubset(mesh2, "temp subset", FC_AT_WHOLE_MESH, &subset2);
      fail_unless(rc == FC_SUCCESS, "abort: failed to create subset");
      fc_addMemberToSubset(subset2, 0);
      subsets[0] = subset;
      subsets[1] = subset2;

      // fc_getSubsetDiameter() on whole
      rc = fc_getSubsetDiameter(subset, &temp_diameter, &temp_numPair,
				&temp_pairIDs);
      fail_unless(rc == FC_SUCCESS, "failed to get diameter of subset");
      fail_unless(FC_DBL_EQUIV(temp_diameter, diameters[i]),
		  "mismatch of diameter");
      fail_unless(temp_numPair == numPairs[i], "mismatch of numPair");
      for (k = 0; k < 2*numPairs[i]; k++)
	fail_unless(temp_pairIDs[k] == pairIDs[i][k], "mismatch of pair ids");
      free(temp_pairIDs);
      // fc_getSubsetsDiameter() on whole single mesh
      rc = fc_getSubsetsDiameter(1, &subset, &temp_diameter, &temp_numPair,
				 &temp_pairIDs, &temp_pairSubsets);
      fail_unless(rc == FC_SUCCESS, "failed to get diameter of subset");
      fail_unless(FC_DBL_EQUIV(temp_diameter, diameters[i]),
		  "mismatch of diameter");
      fail_unless(temp_numPair == numPairs[i], "mismatch of numPair");
      for (k = 0; k < 2*numPairs[i]; k++) {
	fail_unless(temp_pairIDs[k] == pairIDs[i][k], "mismatch of pair ids");
	fail_unless(FC_HANDLE_EQUIV(temp_pairSubsets[k], subset), 
	  "mismatch of pair subset");
      }
      free(temp_pairIDs);
      free(temp_pairSubsets);
      // fc_getSubsetsDiameter() on two whole meshes
      rc = fc_getSubsetsDiameter(2, subsets, &temp_diameter, &temp_numPair,
				 &temp_pairIDs, &temp_pairSubsets);
      fail_unless(rc == FC_SUCCESS, "failed to get diameter of subset");
      fail_unless(FC_DBL_EQUIV(temp_diameter, meshes_diameters[i]),
		  "mismatch of diameter");
      fail_unless(temp_numPair == numPairs[i], "mismatch of numPair");
      for (k = 0; k < 2*numPairs[i]; k++) 
	fail_unless(temp_pairIDs[k] == meshes_pairIDs[i][k], 
		    "mismatch of pair ids");
      for (k = 0; k < numPairs[i]; k++) {
	fail_unless(FC_HANDLE_EQUIV(temp_pairSubsets[2*k], subsets[0]), 
		    "mismatch of pair subset");
	fail_unless(FC_HANDLE_EQUIV(temp_pairSubsets[2*k+1], subsets[1]), 
		    "mismatch of pair subset");
      }
      free(temp_pairIDs);
      free(temp_pairSubsets);

      // fc_getDisplacedSubsetDiameter() 
      rc = fc_getDisplacedSubsetDiameter(subset, displVar, &temp_diameter, 
					 &temp_numPair, &temp_pairIDs);
      fail_unless(rc == FC_SUCCESS, "failed to get diameter of subset");
      fail_unless(FC_DBL_EQUIV(temp_diameter, 2*diameters[i]),
		  "mismatch of diameter");
      fail_unless(temp_numPair == numPairs[i], "mismatch of numPair");
      for (k = 0; k < 2*numPairs[i]; k++)
	fail_unless(temp_pairIDs[k] == pairIDs[i][k], "mismatch of pair ids");
      free(temp_pairIDs);
      // fc_getDisplacedSubsetsDiameter() on single whole mesh
      rc = fc_getDisplacedSubsetsDiameter(1, &subset, &displVar,
					  &temp_diameter, &temp_numPair, 
					  &temp_pairIDs, &temp_pairSubsets);
      fail_unless(rc == FC_SUCCESS, "failed to get diameter of subset");
      fail_unless(FC_DBL_EQUIV(temp_diameter, 2*diameters[i]),
		  "mismatch of diameter");
      fail_unless(temp_numPair == numPairs[i], "mismatch of numPair");
      for (k = 0; k < 2*numPairs[i]; k++) {
	fail_unless(temp_pairIDs[k] == pairIDs[i][k], "mismatch of pair ids");
	fail_unless(FC_HANDLE_EQUIV(temp_pairSubsets[k], subset), 
		    "mismatch of pair subset");
      }
      free(temp_pairIDs);
      free(temp_pairSubsets);
      // fc_getDisplacedSubsetsDiameter() on two whole meshes
      rc = fc_getDisplacedSubsetsDiameter(2, subsets, displVars, 
					  &temp_diameter, &temp_numPair,
					  &temp_pairIDs, &temp_pairSubsets);
      fail_unless(rc == FC_SUCCESS, "failed to get diameter of subset");
      fail_unless(FC_DBL_EQUIV(temp_diameter, 2*meshes_diameters[i]),
		  "mismatch of diameter");
      fail_unless(temp_numPair == numPairs[i], "mismatch of numPair");
      for (k = 0; k < 2*numPairs[i]; k++) 
	fail_unless(temp_pairIDs[k] == meshes_pairIDs[i][k], 
		    "mismatch of pair ids");
      for (k = 0; k < numPairs[i]; k++) {
	fail_unless(FC_HANDLE_EQUIV(temp_pairSubsets[2*k], subsets[0]), 
		    "mismatch of pair subset");
	fail_unless(FC_HANDLE_EQUIV(temp_pairSubsets[2*k+1], subsets[1]), 
		    "mismatch of pair subset");
      }
      free(temp_pairIDs);
      free(temp_pairSubsets);

      // cleanup after FC_AT_WHOLE_MESH subset tests
      fc_deleteSubset(subset);
      fc_deleteSubset(subset2);

      // --- now test different subset types besides WHOLE_MESH

      // loop over assocs except FC_AT_WHOLE_MESH (already tested above)
      for (k = 0; k < numAssoc; k++) {
        int numMember, *memberIDs, elemID;
        // skip this association type if it doesn't exist on the mesh
        fc_getMeshNumEntity(mesh, assocs[k], &numEntity);
        if (numEntity < 1)
          continue;

        // create subset & subset2 - entities that make up last element
	elemID = numElemPerType[i] - 1;
        fc_changeMeshEntityType(mesh, FC_AT_ELEMENT, 1, &elemID,
                                assocs[k], 1, &numMember, &memberIDs);
        rc = fc_createSubset(mesh, "temp subset", assocs[k], &subset);
        fail_unless(rc == FC_SUCCESS, "abort: failed to create subset");
        rc = fc_addArrayMembersToSubset(subset, numMember, memberIDs);
        fail_unless(rc == FC_SUCCESS, 
		    "abort: failed to add subset members");
        rc = fc_createSubset(mesh2, "temp subset", assocs[k], &subset2);
        fail_unless(rc == FC_SUCCESS, "abort: failed to create subset");
        rc = fc_addArrayMembersToSubset(subset2, numMember, memberIDs);
        fail_unless(rc == FC_SUCCESS, 
		    "abort: failed to add subset members");
	free(memberIDs);
	// create subset 3 - entities that make up first element
	elemID = 0;
        fc_changeMeshEntityType(mesh, FC_AT_ELEMENT, 1, &elemID,
                                assocs[k], 1, &numMember, &memberIDs);
        rc = fc_createSubset(mesh, "temp subset", assocs[k], &subset3);
        fail_unless(rc == FC_SUCCESS, "abort: failed to create subset");
        rc = fc_addArrayMembersToSubset(subset3, numMember, memberIDs);
        fail_unless(rc == FC_SUCCESS, 
		    "abort: failed to add subset members");
 	free(memberIDs);
	// subsets from different meshes
	subsets[0] = subset;
	subsets[1] = subset2;
	// subset from same mesh
	subs_same_m[0] = subset3;
	subs_same_m[1] = subset;
	displs_same_m[0] = displVar;
	displs_same_m[1] = displVar;

	// fc_getSubsetDiameter()
	rc = fc_getSubsetDiameter(subset, &temp_diameter, &temp_numPair, 
				  &temp_pairIDs);
	fail_unless(rc == FC_SUCCESS, "failed to get diameter of subset");
	fail_unless(FC_DBL_EQUIV(temp_diameter, subset_diameters[i]),
		    "mismatch of diameter");
	fail_unless(temp_numPair == subset_numPairs[i], "mismatch of numPair");
	for (m = 0; m < 2*temp_numPair; m++)
	  fail_unless(temp_pairIDs[m] == subset_pairIDs[i][m],
		      "mismatch of pair ids");
	free(temp_pairIDs);
	// fc_getSubsetsDiameter() of single subset
	rc = fc_getSubsetsDiameter(1, &subset, &temp_diameter, &temp_numPair, 
				   &temp_pairIDs, &temp_pairSubsets);
	fail_unless(rc == FC_SUCCESS, "failed to get diameter of subset");
	fail_unless(FC_DBL_EQUIV(temp_diameter, subset_diameters[i]),
		    "mismatch of diameter");
	fail_unless(temp_numPair == subset_numPairs[i], "mismatch of numPair");
	for (m = 0; m < 2*temp_numPair; m++) {
	  fail_unless(temp_pairIDs[m] == subset_pairIDs[i][m],
		      "mismatch of pair ids");
	  fail_unless(FC_HANDLE_EQUIV(temp_pairSubsets[m], subset),
		      "mismatch of pair subsets");
	}
	free(temp_pairIDs);
	free(temp_pairSubsets);
	// fc_getSubsetsDiameter() of subsets from two diff meshes
	rc = fc_getSubsetsDiameter(2, subsets, &temp_diameter, &temp_numPair, 
				   &temp_pairIDs, &temp_pairSubsets);
	fail_unless(rc == FC_SUCCESS, "failed to get diameter of subset");
	fail_unless(FC_DBL_EQUIV(temp_diameter, subsets_diameters[i]),
		    "mismatch of diameter");
	fail_unless(temp_numPair == subsets_numPairs[i], "mismatch of numPair");
	for (m = 0; m < 2*temp_numPair; m++) 
	  fail_unless(temp_pairIDs[m] == subsets_pairIDs[i][m],
		      "mismatch of pair ids");
	for (m = 0; m < temp_numPair; m++) {
	  fail_unless(FC_HANDLE_EQUIV(temp_pairSubsets[2*m], subsets[0]), 
		      "mismatch of pair subset");
	  fail_unless(FC_HANDLE_EQUIV(temp_pairSubsets[2*m+1], subsets[1]), 
		      "mismatch of pair subset");
	}
	free(temp_pairIDs);
	free(temp_pairSubsets);
	// fc_getSubsetsDiameter() of subsets from same mesh
	rc = fc_getSubsetsDiameter(2, subs_same_m, &temp_diameter, &temp_numPair, 
				   &temp_pairIDs, &temp_pairSubsets);
	fail_unless(rc == FC_SUCCESS, "failed to get diameter of subset");
	fail_unless(FC_DBL_EQUIV(temp_diameter, diameters[i]),
		    "mismatch of diameter");
	fail_unless(temp_numPair == numPairs[i], "mismatch of numPair");
	for (m = 0; m < 2*temp_numPair; m++) 
	  fail_unless(temp_pairIDs[m] == meshes_pairIDs[i][m],
		      "mismatch of pair ids");
	for (m = 0; m < temp_numPair; m++) {
	  fail_unless(FC_HANDLE_EQUIV(temp_pairSubsets[2*m], subs_same_m[0]), 
		      "mismatch of pair subset");
	  fail_unless(FC_HANDLE_EQUIV(temp_pairSubsets[2*m+1], subs_same_m[1]), 
		      "mismatch of pair subset");
	}
	free(temp_pairIDs);
	free(temp_pairSubsets);

	// fc_getDisplacedSubsetDiameter()
	rc = fc_getDisplacedSubsetDiameter(subset, displVar, &temp_diameter, 
					   &temp_numPair, &temp_pairIDs);
	fail_unless(rc == FC_SUCCESS, "failed to get diameter of displ subset");
	fail_unless(FC_DBL_EQUIV(temp_diameter, 2*subset_diameters[i]),
		    "mismatch of diameter");
	fail_unless(temp_numPair == subset_numPairs[i], "mismatch of numPair");
	for (m = 0; m < 2*temp_numPair; m++)
	  fail_unless(temp_pairIDs[m] == subset_pairIDs[i][m],
		      "mismatch of pair ids");
	free(temp_pairIDs);
	// fc_getDisplacedSubsetsDiameter() on single subset
	rc = fc_getDisplacedSubsetsDiameter(1, &subset, &displVar, 
					    &temp_diameter, &temp_numPair, 
					    &temp_pairIDs, &temp_pairSubsets);
	fail_unless(rc == FC_SUCCESS, "failed to get diameter of displ subset");
	fail_unless(FC_DBL_EQUIV(temp_diameter, 2*subset_diameters[i]),
		    "mismatch of diameter");
	fail_unless(temp_numPair == subset_numPairs[i], "mismatch of numPair");
	for (m = 0; m < 2*temp_numPair; m++) {
	  fail_unless(temp_pairIDs[m] == subset_pairIDs[i][m],
		      "mismatch of pair ids");
	  fail_unless(FC_HANDLE_EQUIV(temp_pairSubsets[m], subset),
		      "mismatch of pair subsets");
	}
	free(temp_pairIDs);
	free(temp_pairSubsets);
	// fc_getDisplacedSubsetsDiameter() of subsets from two meshes
	rc = fc_getDisplacedSubsetsDiameter(2, subsets, displVars, 
					   &temp_diameter, &temp_numPair, 
					   &temp_pairIDs, &temp_pairSubsets);
	fail_unless(rc == FC_SUCCESS, "failed to get diameter of subset");
	fail_unless(FC_DBL_EQUIV(temp_diameter, 2*subsets_diameters[i]),
		    "mismatch of diameter");
	fail_unless(temp_numPair == subsets_numPairs[i], "mismatch of numPair");
	for (m = 0; m < 2*temp_numPair; m++) 
	  fail_unless(temp_pairIDs[m] == subsets_pairIDs[i][m],
		      "mismatch of pair ids");
	for (m = 0; m < temp_numPair; m++) {
	  fail_unless(FC_HANDLE_EQUIV(temp_pairSubsets[2*m], subsets[0]), 
		      "mismatch of pair subset");
	  fail_unless(FC_HANDLE_EQUIV(temp_pairSubsets[2*m+1], subsets[1]), 
		      "mismatch of pair subset");
	}
	free(temp_pairIDs);
	free(temp_pairSubsets);
	// fc_getDisplacedSubsetsDiameter() of subsets from same mesh
	rc = fc_getDisplacedSubsetsDiameter(2, subs_same_m, displs_same_m,
				   &temp_diameter, &temp_numPair, 
				   &temp_pairIDs, &temp_pairSubsets);
	fail_unless(rc == FC_SUCCESS, "failed to get diameter of subset");
	fail_unless(FC_DBL_EQUIV(temp_diameter, 2*diameters[i]),
		    "mismatch of diameter");
	fail_unless(temp_numPair == numPairs[i], "mismatch of numPair");
	for (m = 0; m < 2*temp_numPair; m++) 
	  fail_unless(temp_pairIDs[m] == meshes_pairIDs[i][m],
		      "mismatch of pair ids");
	for (m = 0; m < temp_numPair; m++) {
	  fail_unless(FC_HANDLE_EQUIV(temp_pairSubsets[2*m], subs_same_m[0]), 
		      "mismatch of pair subset");
	  fail_unless(FC_HANDLE_EQUIV(temp_pairSubsets[2*m+1], subs_same_m[1]), 
		      "mismatch of pair subset");
	}
	free(temp_pairIDs);
	free(temp_pairSubsets);

	// --- special cases : don't have to return pairs
	
	// fc_getSubsetDiameter() w/o pairs is o.k.
	rc = fc_getSubsetDiameter(subset, &temp_diameter, NULL, NULL);
	fail_unless(rc == FC_SUCCESS, "failed to get diameter with null pairs");
	fail_unless(FC_DBL_EQUIV(temp_diameter, subset_diameters[i]),
		    "mismatch of diameter");
	// fc_getSubsetsDiameter() w/o pairs is o.k.
	rc = fc_getSubsetsDiameter(1, &subset, &temp_diameter, NULL, NULL, NULL);
	fail_unless(rc == FC_SUCCESS, "failed to get diameter with null pairs");
	fail_unless(FC_DBL_EQUIV(temp_diameter, subset_diameters[i]),
		    "mismatch of diameter");
 	// fc_getDisplacedSubsetDiameter() w/o pairs is o.k.
	rc = fc_getDisplacedSubsetDiameter(subset, displVar, &temp_diameter, 
					   NULL, NULL);
	fail_unless(rc == FC_SUCCESS, "failed to get diameter with null pairs");
	fail_unless(FC_DBL_EQUIV(temp_diameter, 2*subset_diameters[i]),
		    "mismatch of diameter");
	// fc_getDisplacedSubsetsDiameter() w/o pairs is o.k.
	rc = fc_getDisplacedSubsetsDiameter(1, &subset, &displVar, &temp_diameter,
					    NULL, NULL, NULL);
	fail_unless(rc == FC_SUCCESS, "failed to get diameter with null pairs");
	fail_unless(FC_DBL_EQUIV(temp_diameter, 2*subset_diameters[i]),
		    "mismatch of diameter");
	
        // --- special case2 : empty subsets may cause errors

	// fc_getSubsetDiameter() w/ empty subset is error
	rc = fc_getSubsetDiameter(emptySubset, &temp_diameter, &temp_numPair, 
				  &temp_pairIDs);
	fail_unless(rc != FC_SUCCESS, "should fail for empty subset");
	fail_unless(temp_diameter == -1 && temp_numPair == -1 &&
		    temp_pairIDs == NULL, "failure should return nulls");
	// fc_getSubsetsDiameter() w/ all empty subsets is error
	rc = fc_getSubsetsDiameter(1, &emptySubset, &temp_diameter,
			  &temp_numPair, &temp_pairIDs, &temp_pairSubsets);
	fail_unless(rc != FC_SUCCESS, "should fail for empty subset");
	fail_unless(temp_diameter == -1 && temp_numPair == -1 &&
		    temp_pairIDs == NULL && temp_pairSubsets == NULL, 
		    "failure should return nulls");
	// fc_getSubsetsDiameter() w/ 1 empty and 1 nonempty is ok
	subsets[0] = emptySubset;
	subsets[1] = subset2;
	rc = fc_getSubsetsDiameter(2, subsets, &temp_diameter,
			  &temp_numPair, &temp_pairIDs, &temp_pairSubsets);
	fail_unless(rc == FC_SUCCESS, "some empty subsets is ok");
	fail_unless(FC_DBL_EQUIV(temp_diameter, subset_diameters[i]),
		    "mismatch of diameter");
	fail_unless(temp_numPair == subset_numPairs[i], "mismatch of numPair");
	for (m = 0; m < 2*temp_numPair; m++) {
	  fail_unless(temp_pairIDs[m] == subset_pairIDs[i][m],
		      "mismatch of pair ids");
	  fail_unless(FC_HANDLE_EQUIV(temp_pairSubsets[m], subset2),
		      "mismatch of pair subsets");
	}
	free(temp_pairIDs);
	free(temp_pairSubsets);
	// fc_getDisplacedSubsetDiameter() w/ empty subset is error
	rc = fc_getDisplacedSubsetDiameter(emptySubset, displVar, &temp_diameter, 
					   &temp_numPair, &temp_pairIDs);
	fail_unless(rc != FC_SUCCESS, "should fail for empty subset");
	fail_unless(temp_diameter == -1 && temp_numPair == -1 &&
		    temp_pairIDs == NULL, "failure should return nulls");
	// fc_getDisplacedSubsetsDiameter() w/ all empty subsets is error
	rc = fc_getDisplacedSubsetsDiameter(1, &emptySubset, &displVar,
					    &temp_diameter,
			  &temp_numPair, &temp_pairIDs, &temp_pairSubsets);
	fail_unless(rc != FC_SUCCESS, "should fail for empty subset");
	fail_unless(temp_diameter == -1 && temp_numPair == -1 &&
		    temp_pairIDs == NULL && temp_pairSubsets == NULL, 
		    "failure should return nulls");
	// fc_getSubsetsDiameter() w/ 1 empty and 1 nonempty is ok
	subsets[0] = emptySubset;
	subsets[1] = subset2;
	rc = fc_getDisplacedSubsetsDiameter(2, subsets, displVars, &temp_diameter,
			  &temp_numPair, &temp_pairIDs, &temp_pairSubsets);
	fail_unless(rc == FC_SUCCESS, "some empty subsets is ok");
	fail_unless(FC_DBL_EQUIV(temp_diameter, 2*subset_diameters[i]),
		    "mismatch of diameter");
	fail_unless(temp_numPair == subset_numPairs[i], "mismatch of numPair");
	for (m = 0; m < 2*temp_numPair; m++) {
	  fail_unless(temp_pairIDs[m] == subset_pairIDs[i][m],
		      "mismatch of pair ids");
	  fail_unless(FC_HANDLE_EQUIV(temp_pairSubsets[m], subset2),
		      "mismatch of pair subsets");
	}
	free(temp_pairIDs);
	free(temp_pairSubsets);

	// --- test error conditions

	// fc_getSubsetDiameter() -- bad args
	rc = fc_getSubsetDiameter(badSubset, &temp_diameter, &temp_numPair, 
				  &temp_pairIDs);
	fail_unless(rc != FC_SUCCESS, "should fail for bad subset");
	fail_unless(temp_diameter == -1 && temp_numPair == -1 &&
		    temp_pairIDs == NULL, "failure should return nulls");
	rc = fc_getSubsetDiameter(subset, NULL, &temp_numPair, &temp_pairIDs);
	fail_unless(rc != FC_SUCCESS, "should fail if diameter = NULL");
	fail_unless(temp_numPair == -1 && temp_pairIDs == NULL, 
		    "failure should return nulls");
	rc = fc_getSubsetDiameter(subset, &temp_diameter, NULL, &temp_pairIDs);
	fail_unless(rc != FC_SUCCESS, "should fail if 1 pair arg is NULL");
	fail_unless(temp_diameter == -1 && temp_pairIDs == NULL, 
		    "failure should return nulls");
	rc = fc_getSubsetDiameter(subset, &temp_diameter, &temp_numPair, NULL);
	fail_unless(rc != FC_SUCCESS, "should fail if 1 pair arg is NULL");
	fail_unless(temp_diameter == -1 && temp_numPair == -1,
		    "failure should return nulls");

	// fc_getSubsetDiameters() -- bad args
	rc = fc_getSubsetsDiameter(0, &subset, &temp_diameter, &temp_numPair, 
				   &temp_pairIDs, &temp_pairSubsets);
	fail_unless(rc != FC_SUCCESS, "should fail if no subsets");
	fail_unless(temp_diameter == -1 && temp_numPair == -1 &&
		    temp_pairIDs == NULL && temp_pairSubsets == NULL, 
		    "failure should return nulls");
	rc = fc_getSubsetsDiameter(1, NULL, &temp_diameter, &temp_numPair, 
				   &temp_pairIDs, &temp_pairSubsets);
	fail_unless(rc != FC_SUCCESS, "should fail if no subsets");
	fail_unless(temp_diameter == -1 && temp_numPair == -1 &&
		    temp_pairIDs == NULL && temp_pairSubsets == NULL, 
		    "failure should return nulls");
	rc = fc_getSubsetsDiameter(1, &badSubset, &temp_diameter, &temp_numPair, 
				   &temp_pairIDs, &temp_pairSubsets);
	fail_unless(rc != FC_SUCCESS, "should fail for bad subset");
	fail_unless(temp_diameter == -1 && temp_numPair == -1 &&
		    temp_pairIDs == NULL && temp_pairSubsets == NULL, 
		    "failure should return nulls");
	rc = fc_getSubsetsDiameter(1, &subset, NULL, &temp_numPair, 
				   &temp_pairIDs, &temp_pairSubsets);
	fail_unless(rc != FC_SUCCESS, "should fail if diameter = NULL");
	fail_unless(temp_numPair == -1 && temp_pairIDs == NULL &&
		    temp_pairSubsets == NULL, "failure should return nulls");
	rc = fc_getSubsetsDiameter(1, &subset, &temp_diameter, NULL, 
				   &temp_pairIDs, &temp_pairSubsets);
	fail_unless(rc != FC_SUCCESS, "should fail if 1 pair arg is NULL");
	fail_unless(temp_diameter == -1 && temp_pairIDs == NULL &&
		    temp_pairSubsets == NULL, 
		    "failure should return nulls");
	rc = fc_getSubsetsDiameter(1, &subset, &temp_diameter, &temp_numPair, 
				   NULL, &temp_pairSubsets);
	fail_unless(rc != FC_SUCCESS, "should fail if 1 pair arg is NULL");
	fail_unless(temp_diameter == -1 && temp_numPair == -1 &&
		    temp_pairSubsets == NULL, "failure should return nulls");
	rc = fc_getSubsetsDiameter(1, &subset, &temp_diameter, &temp_numPair, 
				   &temp_pairIDs, NULL);
	fail_unless(rc != FC_SUCCESS, "should fail if 1 pair arg is NULL");
	fail_unless(temp_diameter == -1 && temp_numPair == -1 &&
		    temp_pairIDs == NULL, "failure should return nulls");

	// fc_getDisplacedSubsetDiameter() -- bad args
	rc = fc_getDisplacedSubsetDiameter(badSubset, displVar, &temp_diameter, 
					   &temp_numPair, &temp_pairIDs);
	fail_unless(rc != FC_SUCCESS, "should fail for bad subset");
	fail_unless(temp_diameter == -1 && temp_numPair == -1 &&
		    temp_pairIDs == NULL, "failure should return nulls");
	rc = fc_getDisplacedSubsetDiameter(subset, badVar, &temp_diameter, 
					   &temp_numPair, &temp_pairIDs);
	fail_unless(rc != FC_SUCCESS, "should fail for bad var");
	fail_unless(temp_diameter == -1 && temp_numPair == -1 &&
		    temp_pairIDs == NULL, "failure should return nulls");
	rc = fc_getDisplacedSubsetDiameter(subset, displVar, NULL,
					   &temp_numPair, &temp_pairIDs);
	fail_unless(rc != FC_SUCCESS, "should fail if diameter = NULL");
	fail_unless(temp_numPair == -1 && temp_pairIDs == NULL, 
		    "failure should return nulls");
	rc = fc_getDisplacedSubsetDiameter(subset, displVar, &temp_diameter,
					   NULL, &temp_pairIDs);
	fail_unless(rc != FC_SUCCESS, "should fail if 1 pair arg is NULL");
	fail_unless(temp_diameter == -1 && temp_pairIDs == NULL, 
		    "failure should return nulls");
	rc = fc_getDisplacedSubsetDiameter(subset, displVar, &temp_diameter, 
					   &temp_numPair, NULL);
	fail_unless(rc != FC_SUCCESS, "should fail if 1 pair arg is NULL");
	fail_unless(temp_diameter == -1 && temp_numPair == -1,
		    "failure should return nulls");

	// fc_getDisplacedSubsetDiameters() -- bad args
	rc = fc_getDisplacedSubsetsDiameter(0, &subset, &displVar, 
				   &temp_diameter, &temp_numPair, 
				   &temp_pairIDs, &temp_pairSubsets);
	fail_unless(rc != FC_SUCCESS, "should fail if no subsets");
	fail_unless(temp_diameter == -1 && temp_numPair == -1 &&
		    temp_pairIDs == NULL && temp_pairSubsets == NULL, 
		    "failure should return nulls");
	rc = fc_getDisplacedSubsetsDiameter(1, NULL, &displVar, &temp_diameter, 
				   &temp_numPair, &temp_pairIDs, &temp_pairSubsets);
	fail_unless(rc != FC_SUCCESS, "should fail if no subsets");
	fail_unless(temp_diameter == -1 && temp_numPair == -1 &&
		    temp_pairIDs == NULL && temp_pairSubsets == NULL, 
		    "failure should return nulls");
	rc = fc_getDisplacedSubsetsDiameter(1, &badSubset, &displVar,
				   &temp_diameter, &temp_numPair, 
				   &temp_pairIDs, &temp_pairSubsets);
	fail_unless(rc != FC_SUCCESS, "should fail for bad subset");
	fail_unless(temp_diameter == -1 && temp_numPair == -1 &&
		    temp_pairIDs == NULL && temp_pairSubsets == NULL, 
		    "failure should return nulls");
	rc = fc_getDisplacedSubsetsDiameter(1, &subset, NULL,
				   &temp_diameter, &temp_numPair, 
				   &temp_pairIDs, &temp_pairSubsets);
	fail_unless(rc != FC_SUCCESS, "should fail for null var");
	fail_unless(temp_diameter == -1 && temp_numPair == -1 &&
		    temp_pairIDs == NULL && temp_pairSubsets == NULL, 
		    "failure should return nulls");
	rc = fc_getDisplacedSubsetsDiameter(1, &subset, &badVar,
				   &temp_diameter, &temp_numPair, 
				   &temp_pairIDs, &temp_pairSubsets);
	fail_unless(rc != FC_SUCCESS, "should fail for bad var");
	fail_unless(temp_diameter == -1 && temp_numPair == -1 &&
		    temp_pairIDs == NULL && temp_pairSubsets == NULL, 
		    "failure should return nulls");
	rc = fc_getDisplacedSubsetsDiameter(1, &subset, &displVar, NULL, 
				   &temp_numPair, &temp_pairIDs, &temp_pairSubsets);
	fail_unless(rc != FC_SUCCESS, "should fail if diameter = NULL");
	fail_unless(temp_numPair == -1 && temp_pairIDs == NULL &&
		    temp_pairSubsets == NULL, "failure should return nulls");
	rc = fc_getDisplacedSubsetsDiameter(1, &subset, &displVar, 
				   &temp_diameter, NULL, 
				   &temp_pairIDs, &temp_pairSubsets);
	fail_unless(rc != FC_SUCCESS, "should fail if 1 pair arg is NULL");
	fail_unless(temp_diameter == -1 && temp_pairIDs == NULL &&
		    temp_pairSubsets == NULL, 
		    "failure should return nulls");
	rc = fc_getDisplacedSubsetsDiameter(1, &subset, &displVar, 
				   &temp_diameter, &temp_numPair, 
				   NULL, &temp_pairSubsets);
	fail_unless(rc != FC_SUCCESS, "should fail if 1 pair arg is NULL");
	fail_unless(temp_diameter == -1 && temp_numPair == -1 &&
		    temp_pairSubsets == NULL, "failure should return nulls");
	rc = fc_getDisplacedSubsetsDiameter(1, &subset, &displVar,
				   &temp_diameter, &temp_numPair, 
				   &temp_pairIDs, NULL);
	fail_unless(rc != FC_SUCCESS, "should fail if 1 pair arg is NULL");
	fail_unless(temp_diameter == -1 && temp_numPair == -1 &&
		    temp_pairIDs == NULL, "failure should return nulls");

        // cleanup
        fc_deleteSubset(subset);
	fc_deleteSubset(subset2);
      }

      // --- test error conditions
      
      // fc_getMeshDiameter() -- bad args
      rc = fc_getMeshDiameter(badMesh, &temp_diameter, &temp_numPair, 
			      &temp_pairIDs);
      fail_unless(rc != FC_SUCCESS, "should fail for bad mesh");
      fail_unless(temp_diameter == -1 && temp_numPair == -1 &&
		  temp_pairIDs == NULL, "failure should return nulls");
      rc = fc_getMeshDiameter(mesh, NULL, &temp_numPair, &temp_pairIDs);
      fail_unless(rc != FC_SUCCESS, "should fail if diameter = NULL");
      fail_unless(temp_diameter == -1 && temp_numPair == -1 &&
		  temp_pairIDs == NULL, "failure should return nulls");
      rc = fc_getMeshDiameter(mesh, &temp_diameter, NULL, &temp_pairIDs);
      fail_unless(rc != FC_SUCCESS, "should fail if 1 pair arg is NULL");
      fail_unless(temp_diameter == -1 && temp_pairIDs == NULL, 
		  "failure should return nulls");
      rc = fc_getMeshDiameter(mesh, &temp_diameter, &temp_numPair, NULL);
      fail_unless(rc != FC_SUCCESS, "should fail if 1 pair arg is NULL");
      fail_unless(temp_diameter == -1 && temp_numPair == -1,
		  "failure should return nulls");

      // fc_getMeshesDiameter() -- bad args
      rc = fc_getMeshesDiameter(0, &mesh, &temp_diameter, &temp_numPair, 
			      &temp_pairIDs, &temp_pairMeshes);
      fail_unless(rc != FC_SUCCESS, "should fail for 0 meshes");
      fail_unless(temp_diameter == -1 && temp_numPair == -1 &&
		  temp_pairIDs == NULL && temp_pairMeshes == NULL,
		  "failure should return nulls");
      rc = fc_getMeshesDiameter(1, NULL, &temp_diameter, &temp_numPair, 
			      &temp_pairIDs, &temp_pairMeshes);
      fail_unless(rc != FC_SUCCESS, "should fail for NULL meshes");
      fail_unless(temp_diameter == -1 && temp_numPair == -1 &&
		  temp_pairIDs == NULL && temp_pairMeshes == NULL,
		  "failure should return nulls");
      rc = fc_getMeshesDiameter(1, &badMesh, &temp_diameter, &temp_numPair,
				&temp_pairIDs, &temp_pairMeshes);
      fail_unless(rc != FC_SUCCESS, "should fail bad meshes");
      fail_unless(temp_diameter == -1 && temp_numPair == -1 &&
		  temp_pairIDs == NULL && temp_pairMeshes == NULL,
		  "failure should return nulls");
      rc = fc_getMeshesDiameter(1, &mesh, NULL, &temp_numPair, &temp_pairIDs,
				&temp_pairMeshes);
      fail_unless(rc != FC_SUCCESS, "should fail if diameter = NULL");
      fail_unless(temp_diameter == -1 && temp_numPair == -1 &&
		  temp_pairIDs == NULL, "failure should return nulls");
      rc = fc_getMeshesDiameter(1, &mesh, &temp_diameter, NULL, &temp_pairIDs,
				&temp_pairMeshes);
      fail_unless(rc != FC_SUCCESS, "should fail if 1 pair arg is NULL");
      fail_unless(temp_diameter == -1 && temp_pairIDs == NULL &&
		  temp_pairMeshes == NULL, "failure should return nulls");
      rc = fc_getMeshesDiameter(1, &mesh, &temp_diameter, &temp_numPair, 
				NULL, &temp_pairMeshes);
      fail_unless(rc != FC_SUCCESS, "should fail if 1 pair arg is NULL");
      fail_unless(temp_diameter == -1 && temp_numPair == -1 &&
		  temp_pairMeshes == NULL, "failure should return nulls");
      rc = fc_getMeshesDiameter(1, &mesh, &temp_diameter, &temp_numPair, 
				&temp_pairIDs, NULL);
      fail_unless(rc != FC_SUCCESS, "should fail if 1 pair arg is NULL");
      fail_unless(temp_diameter == -1 && temp_numPair == -1 &&
		  temp_pairMeshes == NULL, "failure should return nulls");

      // fc_getDisplacedMeshDiameter() -- bad args
      rc = fc_getDisplacedMeshDiameter(badMesh, displVar, &temp_diameter, 
				       &temp_numPair, &temp_pairIDs);
      fail_unless(rc != FC_SUCCESS, "should fail for bad mesh");
      fail_unless(temp_diameter == -1 && temp_numPair == -1 &&
		  temp_pairIDs == NULL, "failure should return nulls");
      rc = fc_getDisplacedMeshDiameter(mesh, badVar, &temp_diameter, 
				       &temp_numPair, &temp_pairIDs);
      fail_unless(rc != FC_SUCCESS, "should fail for bad mesh");
      fail_unless(temp_diameter == -1 && temp_numPair == -1 &&
		  temp_pairIDs == NULL, "failure should return nulls");
      rc = fc_getDisplacedMeshDiameter(mesh, displVar, NULL, &temp_numPair, 
				       &temp_pairIDs);
      fail_unless(rc != FC_SUCCESS, "should fail if diameter = NULL");
      fail_unless(temp_numPair == -1 && temp_pairIDs == NULL, 
		  "failure should return nulls");
      rc = fc_getDisplacedMeshDiameter(mesh, displVar, &temp_diameter, NULL, 
				       &temp_pairIDs);
      fail_unless(rc != FC_SUCCESS, "should fail if 1 pair arg is NULL");
      fail_unless(temp_diameter == -1 && temp_pairIDs == NULL, 
		  "failure should return nulls");
      rc = fc_getDisplacedMeshDiameter(mesh, displVar, &temp_diameter, 
				       &temp_numPair, NULL);
      fail_unless(rc != FC_SUCCESS, "should fail if 1 pair arg is NULL");
      fail_unless(temp_diameter == -1 && temp_numPair == -1,
		  "failure should return nulls");
 
      // fc_getDisplacedMeshesDiameter() -- bad args
      rc = fc_getDisplacedMeshesDiameter(0, &mesh, &displVar, &temp_diameter,
			      &temp_numPair, &temp_pairIDs, &temp_pairMeshes);
      fail_unless(rc != FC_SUCCESS, "should fail for 0 meshes");
      fail_unless(temp_diameter == -1 && temp_numPair == -1 &&
		  temp_pairIDs == NULL && temp_pairMeshes == NULL,
		  "failure should return nulls");
      rc = fc_getDisplacedMeshesDiameter(1, NULL, &displVar, &temp_diameter, 
			      &temp_numPair, &temp_pairIDs, &temp_pairMeshes);
      fail_unless(rc != FC_SUCCESS, "should fail for NULL meshes");
      fail_unless(temp_diameter == -1 && temp_numPair == -1 &&
		  temp_pairIDs == NULL && temp_pairMeshes == NULL,
		  "failure should return nulls");
      rc = fc_getDisplacedMeshesDiameter(1, &badMesh, &displVar,
				&temp_diameter, &temp_numPair,
				&temp_pairIDs, &temp_pairMeshes);
      fail_unless(rc != FC_SUCCESS, "should fail bad meshes");
      fail_unless(temp_diameter == -1 && temp_numPair == -1 &&
		  temp_pairIDs == NULL && temp_pairMeshes == NULL,
		  "failure should return nulls");
      rc = fc_getDisplacedMeshesDiameter(1, &mesh, NULL, &temp_diameter, 
                                &temp_numPair, &temp_pairIDs, &temp_pairMeshes);
      fail_unless(rc != FC_SUCCESS, "should fail null displ var");
      fail_unless(temp_diameter == -1 && temp_numPair == -1 &&
		  temp_pairIDs == NULL && temp_pairMeshes == NULL,
		  "failure should return nulls");
      rc = fc_getDisplacedMeshesDiameter(1, &mesh, &badVar, &temp_diameter, 
                                &temp_numPair, &temp_pairIDs, &temp_pairMeshes);
      fail_unless(rc != FC_SUCCESS, "should fail if bad var");
      fail_unless(temp_diameter == -1 && temp_numPair == -1 &&
		  temp_pairIDs == NULL && temp_pairMeshes == NULL,
		  "failure should return nulls");
      rc = fc_getDisplacedMeshesDiameter(1, &mesh, &displVar, NULL, 
				&temp_numPair, &temp_pairIDs, &temp_pairMeshes);
      fail_unless(rc != FC_SUCCESS, "should fail if diameter = NULL");
      fail_unless(temp_diameter == -1 && temp_numPair == -1 &&
		  temp_pairIDs == NULL, "failure should return nulls");
      rc = fc_getDisplacedMeshesDiameter(1, &mesh, &displVar, &temp_diameter,
				NULL, &temp_pairIDs, &temp_pairMeshes);
      fail_unless(rc != FC_SUCCESS, "should fail if 1 pair arg is NULL");
      fail_unless(temp_diameter == -1 && temp_pairIDs == NULL &&
		  temp_pairMeshes == NULL, "failure should return nulls");
      rc = fc_getDisplacedMeshesDiameter(1, &mesh, &displVar, &temp_diameter,
				&temp_numPair, NULL, &temp_pairMeshes);
      fail_unless(rc != FC_SUCCESS, "should fail if 1 pair arg is NULL");
      fail_unless(temp_diameter == -1 && temp_numPair == -1 &&
		  temp_pairMeshes == NULL, "failure should return nulls");
      rc = fc_getDisplacedMeshesDiameter(1, &mesh, &displVar, &temp_diameter,
				&temp_numPair, &temp_pairIDs, NULL);
      fail_unless(rc != FC_SUCCESS, "should fail if 1 pair arg is NULL");
      fail_unless(temp_diameter == -1 && temp_numPair == -1 &&
		  temp_pairMeshes == NULL, "failure should return nulls");

      // --- cleanup
      
      fc_deleteMesh(mesh);
    }
  }

  rc = fc_deleteDataset(dataset);
  fail_unless(rc == FC_SUCCESS, "failed to close dataset at end of test");
}
END_TEST

// test getting centroid
// basic layout of vertices is on two side-by-side hexes
// for 1D use only lower edge, for 2D use bottom faces
START_TEST(get_centroid) {
  FC_ReturnCode rc;
  int i, j, k, m, n, p;
  double side_length = 2.3;
  char name[20] = "blue berry";
  FC_Dataset dataset;
  FC_Mesh mesh, badMesh = { 999, 999 };
  FC_Subset subset, zeroSubset, badSubset = { 999, 999 };
  FC_Variable displVar, weightsVars[3], badVar = { 999, 999 };
  FC_MathType mathType;
  int numElemType = 4;
  FC_ElementType elemTypes[4] = { FC_ET_POINT,  FC_ET_LINE,
                                  FC_ET_QUAD,   FC_ET_HEX };
  int dims[3] = { 3, 2, 2 };
  int minDimPerType[4] = { 1, 1, 2, 3 };
  int numVertPerType[4] = { 3, 3, 6, 12}, numElemPerType[4] = { 3, 2, 2, 2 };
  int conns[4][16] = { { 0,    1,    2 },
                       { 0, 1,    1, 2 },
                       { 0, 1, 4, 3,    1, 2, 5, 4 },
                       { 0, 1, 4, 3, 6, 7, 10, 9,    1, 2, 5, 4, 7, 8, 11, 10
                       } };
  double coords[3][3*12];
  int intWeights[12] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 };
  float floatWeights[12] = { 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12. };
  double doubleWeights[12] = { 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12. };
  void *weights[3] = { intWeights, floatWeights, doubleWeights };
  double weight_sum;
  int temp_dim;
  FC_Coords mesh_centroid = { (dims[0]-1)*side_length/2., 
			      (dims[1]-1)*side_length/2.,
			      (dims[2]-1)*side_length/2. };
  FC_Coords subset_centroid = { (side_length + (dims[0]-1)*side_length)/2.,
				mesh_centroid[1],
				mesh_centroid[2] };
  FC_Coords weighted_centroid, displ_weighted_centroid, temp_centroid;
  FC_Coords nullPoint = { -1., -1., -1 };
  int numAssoc = 4, numEntity;
  FC_AssociationType assocs[4] = { FC_AT_VERTEX, FC_AT_EDGE, 
                                   FC_AT_FACE, FC_AT_ELEMENT };
  int numDataType = 3;
  FC_DataType dataTypes[3] = { FC_DT_INT, FC_DT_FLOAT, FC_DT_DOUBLE };
  
  // setup coords
  // 1D coords
  for (i = 0; i < dims[0]; i++) 
    coords[0][i] = i*side_length;
  // 2D coords
  for (j = 0; j < dims[1]; j++) {
    for (i = 0; i < dims[0]; i++) {
      int ID = j*dims[0] + i;
      coords[1][2*ID] = i*side_length;
      coords[1][2*ID + 1] = j*side_length;
    }
  }
  // 3D coords
  for (k = 0; k < dims[2]; k++) {
    for (j = 0; j < dims[1]; j++) {
      for (i = 0; i < dims[0]; i++) {
	int ID = k*dims[0]*dims[1] + j*dims[0] + i;
	coords[2][3*ID] = i*side_length;
	coords[2][3*ID+1] = j*side_length;
	coords[2][3*ID+2] = k*side_length;
      }
    }
  }

  // setup - create test mesh
  rc = fc_createDataset("temp.xxx", &dataset);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create dataset");

  // test different types
  for (i = 0; i < numElemType; i++) { // loop over elem types
    for (j = 1; j <= 3; j++) { // loop over space dimensions
      // skip if it doesn't make sense to make this elem type at this dim
      if (minDimPerType[i] > j)
        continue;

      // --- create data structures for testing: mesh & displVar

      // create a mesh 
      rc = fc_createMesh(dataset, name, &mesh);
      fail_unless(rc == FC_SUCCESS, "abort: failed to created mesh");
      rc = fc_setMeshCoords(mesh, j, numVertPerType[i], coords[j-1]);
      fail_unless(rc == FC_SUCCESS, "failed to set vertex coords");
      rc = fc_setMeshElementConns(mesh, elemTypes[i], numElemPerType[i], 
                                  conns[i]);
      fail_unless(rc == FC_SUCCESS, "failed to set element conns");
      // create displacement variable
      rc = fc_createVariable(mesh, "displ", &displVar);
      fail_unless(rc == FC_SUCCESS, "abort: failed to create displ variable");
      if (j == 1)
	mathType = FC_MT_SCALAR;
      else
	mathType = FC_MT_VECTOR;
      rc = fc_setVariableData(displVar, numVertPerType[i], j, FC_AT_VERTEX, 
			      mathType, FC_DT_DOUBLE, (void*)coords[j-1]);
      fail_unless(rc == FC_SUCCESS, "failed to set displVar data");
      // create weights variables
      for (k = 0; k < numDataType; k++) {
	rc = fc_createVariable(mesh, "weights", &weightsVars[k]);
	fail_unless(rc == FC_SUCCESS, "abort: failed to create weights var");
	rc = fc_setVariableData(weightsVars[k], numVertPerType[i], 1,
				FC_AT_VERTEX, FC_MT_SCALAR, dataTypes[k],
				weights[k]);
	fail_unless(rc == FC_SUCCESS, "failed to set weightsVar data");
      }

      // setup weighted centroids
      weight_sum = 0.;
      for (k = 0; k < 3; k++) {
	weighted_centroid[k] = 0.;
	displ_weighted_centroid[k] = 0.;
      }
      for (k = 0; k < numVertPerType[i]; k++) {
	weight_sum += doubleWeights[k];
	for (m = 0; m < j; m++) {
	  weighted_centroid[m] += coords[j-1][j*k+m]*doubleWeights[k];
	  displ_weighted_centroid[m] += 2*coords[j-1][j*k+m]*doubleWeights[k];
	}
      }
      for (k = 0; k < j; k++) {
	weighted_centroid[k] /= weight_sum;
	displ_weighted_centroid[k] /= weight_sum;
      }

      // -- test get centroid on whole mesh (& subset on FC_AT_WHOLE_MESH)

      // for each weightVar dataType
      for (p = 0; p < numDataType; p++) {
	// fc_getMeshCentroid()
	rc = fc_getMeshCentroid(mesh, &temp_dim, &temp_centroid);
	fail_unless(rc == FC_SUCCESS, "failed to get centroid of mesh");
	fail_unless(temp_dim == j, "mismatch of temp_dim");
	for (k = 0; k < minDimPerType[i]; k++) 
	  fail_unless(FC_DBL_EQUIV(temp_centroid[k], mesh_centroid[k]), 
		      "mismatch of centroid");
	for (k = minDimPerType[i]; k < 3; k++)
	  fail_unless(temp_centroid[k] == 0., "mismatch of centroid");
	
	// fc_getDisplacedMeshCentroid()
	rc = fc_getDisplacedMeshCentroid(mesh, displVar, &temp_dim,
					 &temp_centroid);
	fail_unless(rc == FC_SUCCESS, 
		    "failed to get centroid of displaced mesh");
	fail_unless(temp_dim == j, "mismatch of temp_dim");
	for (k = 0; k < minDimPerType[i]; k++) 
	  fail_unless(FC_DBL_EQUIV(temp_centroid[k], 2.*mesh_centroid[k]), 
		      "mismatch of centroid");
	for (k = minDimPerType[i]; k < 3; k++)
	  fail_unless(temp_centroid[k] == 0., "mismatch of centroid");
        
	// fc_getVariableWeightedMeshCentroid()
	rc = fc_getVariableWeightedMeshCentroid(mesh, weightsVars[p], &temp_dim,
						&temp_centroid);
	fail_unless(rc == FC_SUCCESS, 
		    "failed to get weighted centroid of mesh");
	fail_unless(temp_dim == j, "mismatch of temp_dim");
	for (k = 0; k < 3; k++) 
	  fail_unless(FC_DBL_EQUIV(temp_centroid[k], weighted_centroid[k]), 
		      "mismatch of centroid");
	
	// fc_getVariableWeightedDisplacedMeshCentroid()
	rc = fc_getVariableWeightedDisplacedMeshCentroid(mesh, displVar, 
							 weightsVars[p], &temp_dim, &temp_centroid);
	fail_unless(rc == FC_SUCCESS, 
		    "failed to get weighted centroid of displaced mesh");
	fail_unless(temp_dim == j, "mismatch of temp_dim");
	for (k = 0; k < 3; k++) 
	  fail_unless(FC_DBL_EQUIV(temp_centroid[k], 
				   displ_weighted_centroid[k]), "mismatch of centroid");
	
	// setup for FC_AT_WHOLE_MESH subset tests
	rc = fc_createSubset(mesh, "temp subset", FC_AT_WHOLE_MESH, &subset);
	fail_unless(rc == FC_SUCCESS, "abort: failed to create subset");
	fc_addMemberToSubset(subset, 0);
	
	// fc_getSubsetCentroid()
	rc = fc_getSubsetCentroid(subset, &temp_dim, &temp_centroid);
	fail_unless(rc == FC_SUCCESS, "failed to get centroid of subset");
	fail_unless(temp_dim == j, "mismatch of temp_dim");
	for (k = 0; k < minDimPerType[i]; k++) 
	  fail_unless(FC_DBL_EQUIV(temp_centroid[k], mesh_centroid[k]), 
		      "mismatch of centroid");
	for (k = minDimPerType[i]; k < 3; k++)
	  fail_unless(temp_centroid[k] == 0., "mismatch of centroid");
	
	// fc_getDisplacedSubsetCentroid()
	rc = fc_getDisplacedSubsetCentroid(subset, displVar, &temp_dim,
					   &temp_centroid);
	fail_unless(rc == FC_SUCCESS, 
		  "failed to get centroid of displaced subset");
	fail_unless(temp_dim == j, "mismatch of temp_dim");
	for (k = 0; k < minDimPerType[i]; k++) 
	  fail_unless(FC_DBL_EQUIV(temp_centroid[k], 2.*mesh_centroid[k]), 
		      "mismatch of centroid");
	for (k = minDimPerType[i]; k < 3; k++)
	  fail_unless(temp_centroid[k] == 0., "mismatch of centroid");
	
	// fc_getVariableWeightedSubsetCentroid()
	rc = fc_getVariableWeightedSubsetCentroid(subset, weightsVars[p], &temp_dim,
						  &temp_centroid);
	fail_unless(rc == FC_SUCCESS, 
		    "failed to get weighted centroid of subset");
	fail_unless(temp_dim == j, "mismatch of temp_dim");
	for (k = 0; k < 3; k++) 
	  fail_unless(FC_DBL_EQUIV(temp_centroid[k], weighted_centroid[k]), 
		      "mismatch of centroid");
	
	// fc_getVariableWeightedDisplacedSubsetCentroid()
	rc = fc_getVariableWeightedDisplacedSubsetCentroid(subset, displVar, 
							   weightsVars[p], &temp_dim, &temp_centroid);
	fail_unless(rc == FC_SUCCESS, 
		    "failed to get weighted centroid of displaced subset");
	fail_unless(temp_dim == j, "mismatch of temp_dim");
	for (k = 0; k < 3; k++) 
	  fail_unless(FC_DBL_EQUIV(temp_centroid[k], displ_weighted_centroid[k]),
		      "mismatch of centroid");
	
	// cleanup after FC_AT_WHOLE_MESH subset tests
	fc_deleteSubset(subset);
      } // end loop over p numDataType

      // --- now test different subset types
      
      // loop over assocs except FC_AT_WHOLE_MESH (already tested above)
      for (k = 0; k < numAssoc; k++) {
	int numMember, *memberIDs, elemID = numElemPerType[i] - 1;
	int numVert, *vertIDs;
	// skip this association type if it doesn't exist on the mesh
	fc_getMeshNumEntity(mesh, assocs[k], &numEntity);
	if (numEntity < 1)
	  continue;
	
	// create subset - entities that make up last element
	rc = fc_createSubset(mesh, "temp subset", assocs[k], &subset);
	fail_unless(rc == FC_SUCCESS, "abort: failed to create subset");
	fc_changeMeshEntityType(mesh, FC_AT_ELEMENT, 1, &elemID,
				assocs[k], 1, &numMember, &memberIDs);
	rc = fc_addArrayMembersToSubset(subset, numMember, memberIDs);
	fail_unless(rc == FC_SUCCESS, 
		    "abort: failed to add subset members");
	free(memberIDs);
	
	// setup weighted centroids
	fc_changeMeshEntityType(mesh, FC_AT_ELEMENT, 1, &elemID,
				FC_AT_VERTEX, 1, &numVert, &vertIDs);
	weight_sum = 0.;
	for (m = 0; m < 3; m++) {
	  weighted_centroid[m] = 0.;
	  displ_weighted_centroid[m] = 0.;
	}
	for (m = 0; m < numVert; m++) {
	  double weight = doubleWeights[vertIDs[m]];
	  weight_sum += weight;
	  for (n = 0; n < j; n++) {
	    weighted_centroid[n] += coords[j-1][j*vertIDs[m]+n]*weight;
	    displ_weighted_centroid[n] += 2*coords[j-1][j*vertIDs[m]+n]*weight;
	  }
	}
	for (m = 0; m < j; m++) {
	  weighted_centroid[m] /= weight_sum;
	  displ_weighted_centroid[m] /= weight_sum;
	}
	free(vertIDs);
	
	// for each weightVar dataType
	for (p = 0; p < numDataType; p++) {
	  
	  // fc_getSubsetCentroid()
	  rc = fc_getSubsetCentroid(subset, &temp_dim, &temp_centroid);
	  fail_unless(rc == FC_SUCCESS, "failed to get centroid of subset");
	  fail_unless(temp_dim == j, "mismatch of temp_dim");
	  if (i == 0) 
	    fail_unless(temp_centroid[0] == coords[0][dims[0]-1],
			"mismatch of cenroid");
	  else
	    for (m = 0; m < minDimPerType[i]; m++) 
	      fail_unless(FC_DBL_EQUIV(temp_centroid[m], subset_centroid[m]),
			  "mismatch of centroid");
	  for (m = minDimPerType[i]; m < 3; m++)
	    fail_unless(temp_centroid[m] == 0., "mismatch of centroid");
	  
	  // fc_getDisplacedSubsetCentroid()
	  rc = fc_getDisplacedSubsetCentroid(subset, displVar, &temp_dim, 
					     &temp_centroid);
	  fail_unless(rc == FC_SUCCESS, 
		      "failed to get centroid of displaced subset"); 
	  fail_unless(temp_dim == j, "mismatch of temp_dim");
	  if (i == 0) 
	    fail_unless(temp_centroid[0] == 2.*coords[0][dims[0]-1],
			"mismatch of cenroid");
	  else
	    for (m = 0; m < minDimPerType[i]; m++) 
	      fail_unless(FC_DBL_EQUIV(temp_centroid[m], 2.*subset_centroid[m]),
			  "mismatch of centroid");
	  for (m = minDimPerType[i]; m < 3; m++)
	    fail_unless(temp_centroid[m] == 0., "mismatch of centroid");
	  
	  // fc_getVariableWeightedSubsetCentroid()
	  rc = fc_getVariableWeightedSubsetCentroid(subset, weightsVars[p], 
						    &temp_dim, &temp_centroid);
	  fail_unless(rc == FC_SUCCESS, 
		      "failed to get weighted centroid of subset");
	  fail_unless(temp_dim == j, "mismatch of temp_dim");
	  for (m = 0; m < 3; m++) 
	    fail_unless(FC_DBL_EQUIV(temp_centroid[m], weighted_centroid[m]), 
			"mismatch of centroid");
	  
	  // fc_getVariableWeightedDisplacedSubsetCentroid()
	  rc = fc_getVariableWeightedDisplacedSubsetCentroid(subset, displVar,
				      weightsVars[p], &temp_dim, &temp_centroid);
	  fail_unless(rc == FC_SUCCESS, 
		      "failed to get weighted centroid of subset");
	  fail_unless(temp_dim == j, "mismatch of temp_dim");
	  for (m = 0; m < 3; m++) 
	    fail_unless(FC_DBL_EQUIV(temp_centroid[m], 
			 displ_weighted_centroid[m]), "mismatch of centroid");
	}// end loop over p numDataType

	fc_deleteSubset(subset);
      
      } // loop over k numAssoc
      
      
      // --- test error conditions for mesh versions
      
      // fc_getMeshCentroid() -- bad args
      rc = fc_getMeshCentroid(badMesh, &temp_dim, &temp_centroid);
      fail_unless(rc != FC_SUCCESS, "should fail for bad mesh");
      fail_unless(!memcmp(temp_centroid, nullPoint, 3*sizeof(double)) &&
                  temp_dim == -1, "failure should return null");
      rc = fc_getMeshCentroid(mesh, NULL, &temp_centroid);
      fail_unless(rc != FC_SUCCESS, "should fail if dim = NULL");
      fail_unless(!memcmp(temp_centroid, nullPoint, 3*sizeof(double)),
		  "failure should return null");
      rc = fc_getMeshCentroid(mesh, &temp_dim, NULL);
      fail_unless(rc != FC_SUCCESS, "should fail if centroid = null");
      fail_unless(temp_dim == -1, "failure should return null");

      // fc_getDisplacedMeshCentroid() -- bad args
      rc = fc_getDisplacedMeshCentroid(badMesh, displVar, &temp_dim, 
				       &temp_centroid);
      fail_unless(rc != FC_SUCCESS, "should fail for bad mesh");
      fail_unless(!memcmp(temp_centroid, nullPoint, 3*sizeof(double)) &&
                  temp_dim == -1, "failure should return null");
      rc = fc_getDisplacedMeshCentroid(mesh, badVar, &temp_dim, 
				       &temp_centroid);
      fail_unless(rc != FC_SUCCESS, "should fail if bad displ var");
      fail_unless(!memcmp(temp_centroid, nullPoint, 3*sizeof(double)) &&
		  temp_dim == -1, "failure should return null");
      if (j > 1) {
	rc = fc_getDisplacedMeshCentroid(mesh, weightsVars[0], &temp_dim, 
					 &temp_centroid);
	fail_unless(rc != FC_SUCCESS, "should fail if improper displ var");
	fail_unless(!memcmp(temp_centroid, nullPoint, 3*sizeof(double)) &&
		    temp_dim == -1, "failure should return null");
      }
      rc = fc_getDisplacedMeshCentroid(mesh, displVar, NULL, &temp_centroid);
      fail_unless(rc != FC_SUCCESS, "should fail if dim = NULL");
      fail_unless(!memcmp(temp_centroid, nullPoint, 3*sizeof(double)),
		  "failure should return null");
      rc = fc_getDisplacedMeshCentroid(mesh, displVar, &temp_dim, NULL);
      fail_unless(rc != FC_SUCCESS, "should fail if centroid = null");
      fail_unless(temp_dim == -1, "failure should return null");

      // fc_getVariableWeightedMeshCentroid() -- bad args
      rc = fc_getVariableWeightedMeshCentroid(badMesh, weightsVars[0], &temp_dim, 
				       &temp_centroid);
      fail_unless(rc != FC_SUCCESS, "should fail for bad mesh");
      fail_unless(!memcmp(temp_centroid, nullPoint, 3*sizeof(double)) &&
                  temp_dim == -1, "failure should return null");
      rc = fc_getVariableWeightedMeshCentroid(mesh, badVar, &temp_dim, 
				       &temp_centroid);
      fail_unless(rc != FC_SUCCESS, "should fail if bad displ var");
      fail_unless(!memcmp(temp_centroid, nullPoint, 3*sizeof(double)) &&
		  temp_dim == -1, "failure should return null");
      if (j > 1) {
	rc = fc_getVariableWeightedMeshCentroid(mesh, displVar, &temp_dim, 
					 &temp_centroid);
	fail_unless(rc != FC_SUCCESS, "should fail if improper weights var");
	fail_unless(!memcmp(temp_centroid, nullPoint, 3*sizeof(double)) &&
		    temp_dim == -1, "failure should return null");
      }
      rc = fc_getVariableWeightedMeshCentroid(mesh, displVar, NULL, 
					      &temp_centroid);
      fail_unless(rc != FC_SUCCESS, "should fail if dim = NULL");
      fail_unless(!memcmp(temp_centroid, nullPoint, 3*sizeof(double)),
		  "failure should return null");
      rc = fc_getVariableWeightedMeshCentroid(mesh, displVar, &temp_dim, NULL);
      fail_unless(rc != FC_SUCCESS, "should fail if centroid = null");
      fail_unless(temp_dim == -1, "failure should return null");

      // fc_getVariableWeightedDisplacedMeshCentroid() -- bad args
      rc = fc_getVariableWeightedDisplacedMeshCentroid(badMesh, displVar,
				       weightsVars[0], &temp_dim, &temp_centroid);
      fail_unless(rc != FC_SUCCESS, "should fail for bad mesh");
      fail_unless(!memcmp(temp_centroid, nullPoint, 3*sizeof(double)) &&
                  temp_dim == -1, "failure should return null");
      rc = fc_getVariableWeightedDisplacedMeshCentroid(mesh, badVar, 
				     weightsVars[0], &temp_dim, &temp_centroid);
      fail_unless(rc != FC_SUCCESS, "should fail if bad displ var");
      fail_unless(!memcmp(temp_centroid, nullPoint, 3*sizeof(double)) &&
		  temp_dim == -1, "failure should return null");
      rc = fc_getVariableWeightedDisplacedMeshCentroid(mesh, displVar, badVar,
                                                  &temp_dim, &temp_centroid);
      fail_unless(rc != FC_SUCCESS, "should fail if bad weights var");
      fail_unless(!memcmp(temp_centroid, nullPoint, 3*sizeof(double)) &&
		  temp_dim == -1, "failure should return null");
      if (j > 1) {
	rc = fc_getVariableWeightedDisplacedMeshCentroid(mesh, weightsVars[0], 
                                       weightsVars[0], &temp_dim, &temp_centroid);
	fail_unless(rc != FC_SUCCESS, "should fail if improper displ var");
	fail_unless(!memcmp(temp_centroid, nullPoint, 3*sizeof(double)) &&
		    temp_dim == -1, "failure should return null");
	rc = fc_getVariableWeightedDisplacedMeshCentroid(mesh, displVar, 
                                          displVar, &temp_dim, &temp_centroid);
	fail_unless(rc != FC_SUCCESS, "should fail if improper weights var");
 	fail_unless(!memcmp(temp_centroid, nullPoint, 3*sizeof(double)) &&
		    temp_dim == -1, "failure should return null");
      }
      rc = fc_getVariableWeightedDisplacedMeshCentroid(mesh, displVar,
                                             weightsVars[0], NULL, &temp_centroid);
      fail_unless(rc != FC_SUCCESS, "should fail if dim = NULL");
      fail_unless(!memcmp(temp_centroid, nullPoint, 3*sizeof(double)),
		  "failure should return null");
      rc = fc_getVariableWeightedDisplacedMeshCentroid(mesh, displVar, 
					     weightsVars[0], &temp_dim, NULL);
      fail_unless(rc != FC_SUCCESS, "should fail if centroid = null");
      fail_unless(temp_dim == -1, "failure should return null");

      // --- test error conditions for subset versions

      // make a subsets
      rc = fc_createSubset(mesh, "temp subset", FC_AT_VERTEX, &subset);
      fail_unless(rc == FC_SUCCESS, "abort: failed to create subset");
      fc_addMemberToSubset(subset, 0);
      rc = fc_createSubset(mesh, "zero subset", FC_AT_VERTEX, &zeroSubset);
      fail_unless(rc == FC_SUCCESS, "abort: failed to create zero subset");

      // all -- zero members in subset is error
      rc = fc_getSubsetCentroid(zeroSubset, &temp_dim, &temp_centroid);
      fail_unless(rc != FC_SUCCESS, "should fail for zero subset");
      fail_unless(!memcmp(temp_centroid, nullPoint, 3*sizeof(double)) &&
		  temp_dim == -1, "failure should return null");
      rc = fc_getDisplacedSubsetCentroid(zeroSubset, displVar, &temp_dim, 
					 &temp_centroid);
      fail_unless(rc != FC_SUCCESS, "should fail for zero subset");
      fail_unless(!memcmp(temp_centroid, nullPoint, 3*sizeof(double)) &&
		  temp_dim == -1, "failure should return null");
      rc = fc_getVariableWeightedDisplacedSubsetCentroid(zeroSubset, displVar,
				 weightsVars[0], &temp_dim, &temp_centroid);
      fail_unless(rc != FC_SUCCESS, "should fail for zero subset");
      fail_unless(!memcmp(temp_centroid, nullPoint, 3*sizeof(double)) &&
		  temp_dim == -1, "failure should return null");

      // fc_getSubsetCentroid() -- bad args
      rc = fc_getSubsetCentroid(badSubset, &temp_dim, &temp_centroid);
      fail_unless(rc != FC_SUCCESS, "should fail for bad subset");
      fail_unless(!memcmp(temp_centroid, nullPoint, 3*sizeof(double)) &&
		  temp_dim == -1, "failure should return null");
      rc = fc_getSubsetCentroid(subset, NULL, &temp_centroid);
      fail_unless(rc != FC_SUCCESS, "should fail if dim = NULL");
      fail_unless(!memcmp(temp_centroid, nullPoint, 3*sizeof(double)),
		  "failure should return null");
      rc = fc_getSubsetCentroid(subset, &temp_dim, NULL);
      fail_unless(rc != FC_SUCCESS, "should fail if centroid = null");
      fail_unless(temp_dim == -1, "failure should return null");
      rc = fc_getVariableWeightedSubsetCentroid(zeroSubset, weightsVars[0], 
						&temp_dim, &temp_centroid);
      fail_unless(rc != FC_SUCCESS, "should fail for zero subset");
      fail_unless(!memcmp(temp_centroid, nullPoint, 3*sizeof(double)) &&
		  temp_dim == -1, "failure should return null");
      
      // fc_getDisplacedSubsetCentroid() -- bad args
      rc = fc_getDisplacedSubsetCentroid(badSubset, displVar, &temp_dim, 
					 &temp_centroid);
      fail_unless(rc != FC_SUCCESS, "should fail for bad subset");
      fail_unless(!memcmp(temp_centroid, nullPoint, 3*sizeof(double)) &&
		  temp_dim == -1, "failure should return null");
      rc = fc_getDisplacedSubsetCentroid(subset, badVar, &temp_dim, 
					 &temp_centroid);
      fail_unless(rc != FC_SUCCESS, "should fail if bad displ var");
      fail_unless(!memcmp(temp_centroid, nullPoint, 3*sizeof(double)) &&
		  temp_dim == -1, "failure should return null");
      if (j > 1) {
	rc = fc_getDisplacedSubsetCentroid(subset, weightsVars[0], &temp_dim, 
					   &temp_centroid);
	fail_unless(rc != FC_SUCCESS, "should fail if improper displ var");
	fail_unless(!memcmp(temp_centroid, nullPoint, 3*sizeof(double)) &&
		    temp_dim == -1, "failure should return null");
      }
      rc = fc_getDisplacedSubsetCentroid(subset, displVar, NULL, 
					 &temp_centroid);
      fail_unless(rc != FC_SUCCESS, "should fail if dim = NULL");
      fail_unless(!memcmp(temp_centroid, nullPoint, 3*sizeof(double)),
		  "failure should return null");
      rc = fc_getDisplacedSubsetCentroid(subset, displVar, &temp_dim, NULL);
      fail_unless(rc != FC_SUCCESS, "should fail if centroid = null");
      fail_unless(temp_dim == -1, "failure should return null");
      
      // fc_getVariableWeightedSubsetCentroid() -- bad args
      rc = fc_getVariableWeightedSubsetCentroid(badSubset, weightsVars[0], 
						&temp_dim, &temp_centroid);
      fail_unless(rc != FC_SUCCESS, "should fail for bad subset");
      fail_unless(!memcmp(temp_centroid, nullPoint, 3*sizeof(double)) &&
		  temp_dim == -1, "failure should return null");
      rc = fc_getVariableWeightedSubsetCentroid(subset, badVar, &temp_dim, 
						&temp_centroid);
      fail_unless(rc != FC_SUCCESS, "should fail if bad displ var");
      fail_unless(!memcmp(temp_centroid, nullPoint, 3*sizeof(double)) &&
		  temp_dim == -1, "failure should return null");
      if (j > 1) {
	rc = fc_getVariableWeightedSubsetCentroid(subset, displVar, 
						  &temp_dim, &temp_centroid);
	fail_unless(rc != FC_SUCCESS, "should fail if improper weights var");
	fail_unless(!memcmp(temp_centroid, nullPoint, 3*sizeof(double)) &&
		    temp_dim == -1, "failure should return null");
      }
      rc = fc_getVariableWeightedSubsetCentroid(subset, displVar, NULL, 
						&temp_centroid);
      fail_unless(rc != FC_SUCCESS, "should fail if dim = NULL");
      fail_unless(!memcmp(temp_centroid, nullPoint, 3*sizeof(double)),
		  "failure should return null");
      rc = fc_getVariableWeightedSubsetCentroid(subset, displVar, &temp_dim,
						  NULL);
      fail_unless(rc != FC_SUCCESS, "should fail if centroid = null");
      fail_unless(temp_dim == -1, "failure should return null");
      
      // fc_getVariableWeightedDisplacedSubsetCentroid() -- bad args
      rc = fc_getVariableWeightedDisplacedSubsetCentroid(badSubset, displVar,
				 weightsVars[0], &temp_dim, &temp_centroid);
      fail_unless(rc != FC_SUCCESS, "should fail for bad subset");
      fail_unless(!memcmp(temp_centroid, nullPoint, 3*sizeof(double)) &&
		  temp_dim == -1, "failure should return null");
      rc = fc_getVariableWeightedDisplacedSubsetCentroid(subset, badVar, 
		  weightsVars[0], &temp_dim, &temp_centroid);
      fail_unless(rc != FC_SUCCESS, "should fail if bad displ var");
      fail_unless(!memcmp(temp_centroid, nullPoint, 3*sizeof(double)) &&
		  temp_dim == -1, "failure should return null");
      rc = fc_getVariableWeightedDisplacedSubsetCentroid(subset, displVar, 
					badVar, &temp_dim, &temp_centroid);
      fail_unless(rc != FC_SUCCESS, "should fail if bad weights var");
      fail_unless(!memcmp(temp_centroid, nullPoint, 3*sizeof(double)) &&
		  temp_dim == -1, "failure should return null");
      if (j > 1) {
	rc = fc_getVariableWeightedDisplacedSubsetCentroid(subset, 
		 weightsVars[0], weightsVars[0], &temp_dim, &temp_centroid);
	fail_unless(rc != FC_SUCCESS, "should fail if improper displ var");
	fail_unless(!memcmp(temp_centroid, nullPoint, 3*sizeof(double)) &&
		    temp_dim == -1, "failure should return null");
	rc = fc_getVariableWeightedDisplacedSubsetCentroid(subset, displVar, 
				     displVar, &temp_dim, &temp_centroid);
	fail_unless(rc != FC_SUCCESS, "should fail if improper weights var");
	fail_unless(!memcmp(temp_centroid, nullPoint, 3*sizeof(double)) &&
		    temp_dim == -1, "failure should return null");
      }
      rc = fc_getVariableWeightedDisplacedSubsetCentroid(subset, displVar,
				   weightsVars[0], NULL, &temp_centroid);
      fail_unless(rc != FC_SUCCESS, "should fail if dim = NULL");
      fail_unless(!memcmp(temp_centroid, nullPoint, 3*sizeof(double)),
		  "failure should return null");
      rc = fc_getVariableWeightedDisplacedSubsetCentroid(subset, displVar, 
				        weightsVars[0], &temp_dim, NULL);
      fail_unless(rc != FC_SUCCESS, "should fail if centroid = null");
      fail_unless(temp_dim == -1, "failure should return null");

      // --- cleanup
      fc_deleteMesh(mesh);
    }
  }

  rc = fc_deleteDataset(dataset);
  fail_unless(rc == FC_SUCCESS, "failed to close dataset at end of test");
}
END_TEST

// Not testing extreme shapes
START_TEST(measure_cores)
{
  int numDim;
  int i, j;
  double lens[3] =  { 1.1, 1.3, 1.7 }; // side lengths
  double A = 1.43;  // 1.1 * 1.3
  double V = 2.431; // 1.1 * 1.3 * 1.7
  FC_Coords boxPts[8] = { { 0, 0, 0 }, { 1, 0, 0 }, { 1, 1, 0 }, { 0, 1, 0 },
			  { 0, 0, 1 }, { 1, 0, 1 }, { 1, 1, 1 }, { 0, 1, 1 } };
  FC_Coords trans2DBoxPts[8], trans3DBoxPts[8], *transBoxPts;
  int numExtra = 11;
  FC_Coords extraPts[11] = { { 0.5, 0.5, 0 }, { 0.5, 0.5, 0.5 },
			     { 0.25, 0.25, 0 }, { 0.75, 0.25, 0 },
			     { 0.75, 0.75, 0 }, { 0.25, 0.75, 0 },
			     { 0.5, 0, 0 }, { 0.5, 1, 0 }, 
			     { 0.5, 0, 1 }, { 0.5, 1, 1 },  { 0, 0.5, 1 }, };
  FC_Coords trans2DExtraPts[11], trans3DExtraPts[11], *transExtraPts;
  FC_Coords testPts[8];
  double result;

  // FIX? not happy that FC_DBL_EQUIV is failing

  // --- setup, scale, translate & rotate points
  // scale
  for (j = 0; j < 3; j++) {
    for (i = 0; i < 8; i++) {
      trans2DBoxPts[i][j] = boxPts[i][j] * lens[j];
      trans3DBoxPts[i][j] = boxPts[i][j] * lens[j];
    }
    for (i = 0; i < numExtra; i++) {
      trans2DExtraPts[i][j] = extraPts[i][j] * lens[j];
      trans3DExtraPts[i][j] = extraPts[i][j] * lens[j];
    }
  }
  // translate
  for (j = 0; j < 3; j++) {
    for (i = 0; i < 8; i++) {
      trans2DBoxPts[i][j] += (j+1.)/10.;
      trans3DBoxPts[i][j] += (j+1.)/10.;
    }
    for (i = 0; i < numExtra; i++) {
      trans2DExtraPts[i][j] += (j+1.)/10.;
      trans3DExtraPts[i][j] += (j+1.)/10.;
    }
  }
  // rotate
  for (i = 0; i < 8; i++) {
    rotate2D(&trans2DBoxPts[i][0], &trans2DBoxPts[i][1]);
    rotate3D(&trans3DBoxPts[i][0], &trans3DBoxPts[i][1], &trans3DBoxPts[i][2]);
  }
  for (i = 0; i < numExtra; i++) {
    rotate2D(&trans2DExtraPts[i][0], &trans2DExtraPts[i][1]);
    rotate3D(&trans3DExtraPts[i][0], &trans3DExtraPts[i][1], &trans3DExtraPts[i][2]);
  }
  
  // --- tri calcs
  // 1/2 * base * height
  for (numDim = 2; numDim <= 3; numDim++) {
    if (numDim == 2) {
      transBoxPts = trans2DBoxPts;
      transExtraPts = trans2DExtraPts;
    }
    else {
      transBoxPts = trans3DBoxPts;
      transExtraPts = trans3DExtraPts;
    }
    // tri on unit box
    _fc_calcTriArea(numDim, boxPts, &result);
    fail_unless(result == 0.5, "mismatch of tri area");
    // tri on transformed box
    _fc_calcTriArea(numDim, transBoxPts, &result);
    fail_unless(FC_VALUE_EQUIV(result, 0.5*A, 4*DBL_EPSILON, DBL_MIN),
		"mismatch of trans tri area");
    // another tri w/o 90 degree angle
    for (i = 0; i < 2; i++) 
      for (j = 0; j < numDim; j++) 
	testPts[i][j] = transBoxPts[i][j];  
    for (j = 0; j < numDim; j++)
      testPts[2][j] = transExtraPts[0][j];
    _fc_calcTriArea(numDim, testPts, &result);
    fail_unless(FC_VALUE_EQUIV(result, 0.25*A, 4*DBL_EPSILON, DBL_MIN), 
		"mismatch of non-right tri area");
  }

  // --- quad calcs
  for (numDim = 2; numDim <= 3; numDim++) {
    if (numDim == 2) {
      transBoxPts = trans2DBoxPts;
      transExtraPts = trans2DExtraPts;
    }
    else {
      transBoxPts = trans3DBoxPts;
      transExtraPts = trans3DExtraPts;
    }
    // quad on unit box
    _fc_calcQuadArea(numDim, boxPts, &result);
    fail_unless(result == 1., "mismatch of quad area");
    // quad on transformed box
    _fc_calcQuadArea(numDim, transBoxPts, &result);
    fail_unless(FC_VALUE_EQUIV(result, 1.*A, 4*DBL_EPSILON, DBL_MIN), 
		"mismatch of trans quad area");
    // non convex quad on transformed box
    for (i = 0; i < 4; i++)
      for (j = 0; j < numDim; j++) 
	testPts[i][j] = transBoxPts[i][j];
    for (j = 0; j < numDim; j++)
      testPts[3][j] = transExtraPts[3][j];
    _fc_calcQuadArea(numDim, testPts, &result);
    fail_unless(FC_VALUE_EQUIV(result, 0.25*A, 6*DBL_EPSILON, DBL_MIN),
		"mismatch of concave quad area");
    // non convex quad II on transformed box
    for (i = 0; i < 4; i++)
      for (j = 0; j < numDim; j++)
	testPts[i][j] = transBoxPts[i][j];
    for (j = 0; j < numDim; j++)
      testPts[2][j] = transExtraPts[2][j];
    _fc_calcQuadArea(numDim, testPts, &result);
    fail_unless(FC_VALUE_EQUIV(result, 0.25*A, 4*DBL_EPSILON, DBL_MIN), 
		"mismatch of concave quad area II");
    // another quad w/o 90 degree angle
    for (i = 0; i < 2; i++)
      for (j = 0; j < numDim; j++)
	testPts[i][j] = transExtraPts[i+2][j];
    for (i = 2; i < 4; i++)
      for (j = 0; j < numDim; j++)
      testPts[i][j] = transBoxPts[i][j];
    _fc_calcQuadArea(numDim, testPts, &result);
    fail_unless(FC_DBL_EQUIV(result, 9./16*A), "mismatch of non-right quad");
  }

  // --- tet calcs
  // Vol = 1/3 x BaseArea x Height
  // tet on unit box
  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
      testPts[i][j] = boxPts[i][j];
  for (j = 0; j < 3; j++)
    testPts[3][j] = boxPts[5][j];
  _fc_calcTetVolume(testPts, &result);
  fail_unless(FC_DBL_EQUIV(result, 1./6.), "mismatch of tet volume");
  // tet on transformed box
  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
      testPts[i][j] = trans3DBoxPts[i][j];
  for (j = 0; j < 3; j++)
    testPts[3][j] = trans3DBoxPts[5][j];
  _fc_calcTetVolume(testPts, &result);
  fail_unless(FC_VALUE_EQUIV(result, 1./6.*V, 4*DBL_EPSILON, DBL_MIN), 
	      "mismatch of trans3D tet volume");
  // another tet - w/o 90 degree angle
  for (i = 0; i < 2; i++)
    for (j = 0; j < 3; j++)
      testPts[i][j] = trans3DBoxPts[i][j];
  for (j = 0; j < 3; j++) {
    testPts[2][j] = trans3DExtraPts[7][j];
    testPts[3][j] = trans3DExtraPts[1][j];
  }
  _fc_calcTetVolume(testPts, &result);
  fail_unless(FC_VALUE_EQUIV(result, 1./12.*V, 2*DBL_EPSILON, DBL_MIN), 
	      "mismatch of non-right tet volume");

  // --- pyramid calcs
  // Vol = 1/3 x BaseArea x Height
  // pyramid on unit box
  _fc_calcPyramidVolume(boxPts, &result);
  fail_unless(result == 1./3, "mismatch of pyramid volume");
  // pyramid on transformed box
  _fc_calcPyramidVolume(trans3DBoxPts, &result);
  fail_unless(FC_VALUE_EQUIV(result, 1./3.*V, 4*DBL_EPSILON, DBL_MIN), 
	      "mismatch of trans pyramid vol");
  // non convex pyramid
  for (i = 0; i < 5; i++)
    for (j = 0; j < 3; j++)
      testPts[i][j] = trans3DBoxPts[i][j];
  for (j = 0; j < 3; j++)
    testPts[3][j] = trans3DExtraPts[3][j];
  _fc_calcPyramidVolume(testPts, &result);
  fail_unless(FC_VALUE_EQUIV(result, 1./12.*V, 8*DBL_EPSILON, DBL_MIN), 
	      "mismatch of concave pyramid vol");
  // non convex pyramid II
  for (i = 0; i < 5; i++)
    for (j = 0; j < 3; j++)
      testPts[i][j] = trans3DBoxPts[i][j];
  for (j = 0; j < 3; j++)
    testPts[2][j] = trans3DExtraPts[2][j];
  _fc_calcPyramidVolume(testPts, &result);
  fail_unless(FC_VALUE_EQUIV(result, 1./12.*V, 2*DBL_EPSILON, DBL_MIN), 
	      "mismatch of concave pyramid vol II");
  // another pyramid w/o 90 degree angle (use non-right quad for base)
  for (i = 0; i < 2; i++)
    for (j = 0; j < 3; j++)
      testPts[i][j] = trans3DExtraPts[i+2][j];
  for (i = 2; i < 4; i++)
    for (j = 0; j < 3; j++)
      testPts[i][j] = trans3DBoxPts[i][j];
  for (j = 0; j < 3; j++)
    testPts[4][j] = trans3DExtraPts[1][j];
  _fc_calcPyramidVolume(testPts, &result);
  fail_unless(FC_DBL_EQUIV(result, 3./32.*V), "mismatch of pyramid volume");

  // --- prism calcs
  // prism on unit box
  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
      testPts[i][j] = boxPts[i][j];
  for (i = 3; i < 6; i++)
    for (j = 0; j < 3; j++)
      testPts[i][j] = boxPts[i+1][j];
  _fc_calcPrismVolume(testPts, &result);
  fail_unless(result == 1./2, "mismatch of prism volume");
  // prism on transformed box
  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
      testPts[i][j] = trans3DBoxPts[i][j];
  for (i = 3; i < 6; i++)
    for (j = 0; j < 3; j++)
      testPts[i][j] = trans3DBoxPts[i+1][j];
  _fc_calcPrismVolume(testPts, &result);
  fail_unless(FC_VALUE_EQUIV(result, .5*V, 4*DBL_EPSILON, DBL_MIN),
	      "mismatch of prism volume");
  // another prism with differing sized top & bottom
  // volume is vol of big tet projecting through prism - vol of
  // little tet on top = 1/3 - 1/(3*8)
  for (j = 0; j < 3 ; j++) {
    testPts[4][j] = trans3DExtraPts[8][j];
    testPts[5][j] = trans3DExtraPts[10][j];
  }
  _fc_calcPrismVolume(testPts, &result);
  fail_unless(FC_VALUE_EQUIV(result, 7./24.*V, 2*DBL_EPSILON, DBL_MIN), 
	      "mismatch of prism volume");
  // another prism, skewed (will have same volume as non skewed)
  for (j = 0; j < 3; j++) {
    testPts[0][j] = trans3DBoxPts[0][j];
    testPts[1][j] = trans3DExtraPts[6][j];
    testPts[2][j] = trans3DBoxPts[3][j];
    testPts[3][j] = trans3DExtraPts[8][j];
    testPts[4][j] = trans3DBoxPts[5][j];
    testPts[5][j] = trans3DExtraPts[9][j];
  }
  _fc_calcPrismVolume(testPts, &result);
  fail_unless(FC_VALUE_EQUIV(result, .25*V, 4*DBL_EPSILON, DBL_MIN), 
	      "mismatch of screwd prism volume");
  // another prism w/o 90 degree angle? hard to get so faces are planar ...

  // --- hex calcs
  // hex on unit box
  _fc_calcHexVolume(boxPts, &result);
  fail_unless(result == 1., "mismatch of hex volume");
  // hex on transformed box
  _fc_calcHexVolume(trans3DBoxPts, &result);
  fail_unless(FC_VALUE_EQUIV(result, V, 4*DBL_EPSILON, DBL_MIN),
	      "mismatch of hex volume");

  // another hex with differing sized top & bottom
  for (i = 0; i < 4; i++) 
    for (j = 0; j < 3; j++) 
      testPts[i][j] = trans3DExtraPts[i+2][j];
  for (i = 4; i < 8; i++)
    for (j = 0; j < 3 ; j++)
      testPts[i][j] = trans3DBoxPts[i][j];
  _fc_calcHexVolume(testPts, &result);
  fail_unless(FC_VALUE_EQUIV(result, 7./12.*V, 2*DBL_EPSILON, DBL_MIN), 
	      "mismatch of hex volume");
  // another hex, skewed (will have same volume as non skewed)
  for (j = 0; j < 3; j++) {
    testPts[0][j] = trans3DBoxPts[0][j];
    testPts[1][j] = trans3DExtraPts[6][j];
    testPts[2][j] = trans3DExtraPts[7][j];
    testPts[3][j] = trans3DBoxPts[3][j];
    testPts[4][j] = trans3DExtraPts[8][j];
    testPts[5][j] = trans3DBoxPts[5][j];
    testPts[6][j] = trans3DBoxPts[6][j];
    testPts[7][j] = trans3DExtraPts[9][j];
  }
  _fc_calcHexVolume(testPts, &result);
  //printf("** result = %g (%g)\n", result, 0.5*V);
  fail_unless(FC_VALUE_EQUIV(result, 0.5*V, 4*DBL_EPSILON, DBL_MIN), 
	      "mismatch of hex volume");
  // another hex w/o 90 degree angle? hard to get so faces are planar ...
}
END_TEST

// test calculation of edge lengths, face areas & region volumes,
// mesh and mesh subset volumes, also test displaced versions 
START_TEST(length_area_vol)
{
  FC_ReturnCode rc, rc_displ;
  int i, j, k, n;
  FC_Dataset dataset;
  int numMesh = 8;
  char meshNames[8][1024] = { "point mesh", "line mesh", "tri mesh", 
                              "quad mesh", "tet mesh", "pyramid mesh",
                              "prism mesh", "hex mesh" };
  FC_Mesh mesh, badMesh = { 999, 999 };
  FC_Mesh *returnMeshes;
  int numReturnMeshes;
  FC_Subset subset,badSubset = { 999, 999 };
  // from data/gen_tet_prism_hex.vert
  double A = 10., B = 8., C = 5.;  // from data/mesh_generator.c
  int M = 11, N = 12, P = 6;       // from data/mesh_generator.c
  int numKnown_lens, numKnown_areas, numKnown_vols;
  double known_lens[6];  // buffer 
  double known_areas[6]; // buffer
  double known_vols[2];  // buffer
  double known_mesh_vol;
  int numEdge, numFace, numElement;
  FC_Variable coords_var, nonDisplVar, badVar = { 999, 999 };
  FC_Variable edge_lengths, face_areas, element_volumes;
  FC_Variable edge_lengths_displ, face_areas_displ, element_volumes_displ;
  int numData, numComp;
  FC_DataType datatype;
  FC_AssociationType assoc;

  FC_MathType mathtype;
  double *edge_data, *face_data, *element_data;
  double *edge_data_displ, *face_data_displ, *element_data_displ;
  double temp;
  double mesh_volume, mesh_volume_displ, subset_volume, subset_volume_displ;


  // setup
  fc_loadDataset("../data/gen_multimesh.ex2", &dataset);
  rc = fc_getNumMesh(dataset, &numMesh);
  fail_unless(rc == FC_SUCCESS, "Test aborted: couldn't get mesh names");
  fail_unless(numMesh >= 8, "test aborted, expected at least 8 meshes");

  for (i = 0; i < numMesh; i++) {
    // get mesh and make sure it is what we expect
    rc = fc_getMeshByName(dataset, meshNames[i],&numReturnMeshes,&returnMeshes);
    fail_unless(rc == FC_SUCCESS, "aborted - cant get mesh by name");
    fail_unless(numReturnMeshes == 1, "failed to find unique mesh by name");
    mesh = returnMeshes[0];
    free(returnMeshes);

    fail_unless(rc == FC_SUCCESS, "Test aborted: mesh not found");
    rc = fc_getMeshCoordsAsVariable(mesh, &coords_var);
    fail_unless(rc == FC_SUCCESS, "test aborted: could not get coords as var");

    // --- calc answers for 3D meshes

    // setup
    fc_getMeshNumEntity(mesh, FC_AT_EDGE, &numEdge);
    fc_getMeshNumEntity(mesh, FC_AT_FACE, &numFace);
    fc_getMeshNumElement(mesh, &numElement);
    // base case is the hex mesh
    numKnown_lens = 3;
    known_lens[0] = A/M;
    known_lens[1] = B/N;
    known_lens[2] = C/P;
    numKnown_areas = 3;
    known_areas[0] = known_lens[0]*known_lens[1];
    known_areas[1] = known_lens[1]*known_lens[2];
    known_areas[2] = known_lens[2]*known_lens[0];
    numKnown_vols = 1;
    known_vols[0] = known_lens[0]*known_lens[1]*known_lens[2];
    known_mesh_vol = A*B*C;
    switch(i) {
    case 4: // tet mesh - pop off 4 corners and get 5th center tet
      numKnown_lens = 6;
      // add diagonals on every side
      known_lens[3] = sqrt(known_lens[0]*known_lens[0] +
                           known_lens[1]*known_lens[1]);
      known_lens[4] = sqrt(known_lens[1]*known_lens[1] +
                           known_lens[2]*known_lens[2]);
      known_lens[5] = sqrt(known_lens[0]*known_lens[0] +
                           known_lens[2]*known_lens[2]);
      numKnown_areas = 4;
      // divide quad sides in half
      known_areas[0] = known_areas[0]/2.0;
      known_areas[1] = known_areas[1]/2.0;
      known_areas[2] = known_areas[2]/2.0;
      // the new tri on each of the 4 corner tets, and all sides of the center
      // are the same - Use Heron's formula to calc tri volume
      temp = (known_lens[3] + known_lens[4] + known_lens[5]) / 2.0;
      known_areas[3] = sqrt( temp * (temp - known_lens[3]) *
                             (temp - known_lens[4]) * (temp - known_lens[5]));
      // don't know why, but this happens to work
      numKnown_vols = 2;
      temp = known_vols[0];
      known_vols[0] = temp/6.0;
      known_vols[1] = temp/3.0;
      break;
    case 5: // pyramid mesh
      numKnown_lens = 4;
      known_lens[3] = sqrt(known_lens[0]*known_lens[0] +
                           known_lens[1]*known_lens[1] +
                           known_lens[2]*known_lens[2]) / 2.0;
      numKnown_areas += 3;
      known_areas[3] = known_lens[0] * 
                               sqrt(known_lens[1]*known_lens[1]/4.0 +
                                    known_lens[2]*known_lens[2]/4.0) / 2.0;
      known_areas[4] = known_lens[1] * 
                               sqrt(known_lens[0]*known_lens[0]/4.0 +
                                    known_lens[2]*known_lens[2]/4.0) / 2.0;
      known_areas[5] = known_lens[2] * 
                               sqrt(known_lens[0]*known_lens[0]/4.0 +
                                    known_lens[1]*known_lens[1]/4.0) / 2.0;
      known_vols[0] = known_vols[0]/6.;
      break;
    case 6: // prism mesh
      numKnown_lens = 4;
      known_lens[3] = sqrt(known_lens[0]*known_lens[0] +
                           known_lens[1]*known_lens[1]);
      numKnown_areas = 4;
      known_areas[0] = known_areas[0]/2.;
      known_areas[3] = known_lens[2]*known_lens[3];
      known_vols[0] = known_vols[0]/2.;
    }

    // --- test comands (test actual values later)

    // test getting edge lengths
    rc = fc_getEdgeLengths(mesh, &edge_lengths);
    rc_displ = fc_getDisplacedEdgeLengths(mesh, coords_var,
					  &edge_lengths_displ); 
    if (i == 0) { // point mesh
      // regular
      fail_unless(rc != FC_SUCCESS, "should fail to get edges on point mesh");
      fail_unless(FC_HANDLE_EQUIV(edge_lengths, FC_NULL_VARIABLE),
                  "failure should return null handle");
      // displaced
      fail_unless(rc_displ != FC_SUCCESS, 
		  "should fail to get edges on displaced point mesh");
      fail_unless(FC_HANDLE_EQUIV(edge_lengths_displ, FC_NULL_VARIABLE),
                  "failure should return null handle");
    }
    else {
      // regular
      fail_unless(rc == FC_SUCCESS, "failed to get edges on line mesh");
      rc = fc_getVariableInfo(edge_lengths, &numData, &numComp, &assoc, 
                              &mathtype, &datatype);
      fail_unless(numData == numEdge, "mismatch of numData");
      fail_unless(numComp == 1 && mathtype == FC_MT_SCALAR && 
                  datatype == FC_DT_DOUBLE, "should return scalar double");
      if (i == 1) // line mesh
        fail_unless(assoc == FC_AT_ELEMENT, "mismatch of association");
      else
        fail_unless(assoc == FC_AT_EDGE, "mismatch of association");
      // displaced
      fail_unless(rc_displ == FC_SUCCESS, 
		  "failed to get edges on displaced line mesh");
      rc = fc_getVariableInfo(edge_lengths_displ, &numData, &numComp, &assoc, 
                              &mathtype, &datatype);
      fail_unless(numData == numEdge, "mismatch of numData");
      fail_unless(numComp == 1 && mathtype == FC_MT_SCALAR && 
                  datatype == FC_DT_DOUBLE, "should return scalar double");
      if (i == 1) // line mesh
        fail_unless(assoc == FC_AT_ELEMENT, "mismatch of association");
      else
        fail_unless(assoc == FC_AT_EDGE, "mismatch of association");
    }

    // test getting face areas
    rc = fc_getFaceAreas(mesh, &face_areas);
    rc_displ = fc_getDisplacedFaceAreas(mesh, coords_var, &face_areas_displ);
    if (i < 2) { // point mesh or line mesh
      // regular
      fail_unless(rc != FC_SUCCESS, 
                  "should fail to get faces on point or line mesh");
      fail_unless(FC_HANDLE_EQUIV(face_areas, FC_NULL_VARIABLE),
                  "failure should return null handle");
      // displaced
      fail_unless(rc_displ != FC_SUCCESS, 
                  "should fail to get faces on displaced point or line mesh");
      fail_unless(FC_HANDLE_EQUIV(face_areas_displ, FC_NULL_VARIABLE),
                  "failure should return null handle");
    }
    else {
      // regular
      fail_unless(rc == FC_SUCCESS, "failed get faces on quad mesh");
      rc = fc_getVariableInfo(face_areas, &numData, &numComp, &assoc, 
                              &mathtype, &datatype);
      fail_unless(numData == numFace, "mismatch of numData");
      fail_unless(numComp == 1 && mathtype == FC_MT_SCALAR && 
                  datatype == FC_DT_DOUBLE, "should return scalar double");
      if (i == 2 || i == 3) // tri mesh or quad mesh
        fail_unless(assoc == FC_AT_ELEMENT, "mismatch of association");
      else
        fail_unless(assoc == FC_AT_FACE, "mismatch of association");
      // displaced
      fail_unless(rc_displ == FC_SUCCESS, 
		  "failed get faces on displaced quad mesh");
      rc = fc_getVariableInfo(face_areas_displ, &numData, &numComp, &assoc, 
                              &mathtype, &datatype);
      fail_unless(numData == numFace, "mismatch of numData");
      fail_unless(numComp == 1 && mathtype == FC_MT_SCALAR && 
                  datatype == FC_DT_DOUBLE, "should return scalar double");
      if (i == 2 || i == 3) // tri mesh or quad mesh
        fail_unless(assoc == FC_AT_ELEMENT, "mismatch of association");
      else
        fail_unless(assoc == FC_AT_FACE, "mismatch of association");
    }

    //printf("i = %d (%s)\n", i, meshNames[i]);

    // test getting region volumes
    rc = fc_getRegionVolumes(mesh, &element_volumes);
    rc_displ = fc_getDisplacedRegionVolumes(mesh, coords_var,
					    &element_volumes_displ);
    if (i < 4) { // point mesh, line mesh, tri mesh, or quad mesh
      // regular
      fail_unless(rc != FC_SUCCESS, 
                  "should fail to get regions on 2D or less mesh");
      fail_unless(FC_HANDLE_EQUIV(element_volumes, FC_NULL_VARIABLE),
                  "failure should return null handle");
      // displaced
      fail_unless(rc != FC_SUCCESS, 
                  "should fail to get regions on 2D or less displaced mesh");
      fail_unless(FC_HANDLE_EQUIV(element_volumes_displ, FC_NULL_VARIABLE),
                  "failure should return null handle");
    }
    else {
      // regular
      fail_unless(rc == FC_SUCCESS, "failed to get regions on hex mesh");
      rc = fc_getVariableInfo(element_volumes, &numData, &numComp, &assoc, 
                              &mathtype, &datatype);
      fail_unless(numData == numElement, "mismatch of numData");
      fail_unless(numComp == 1 && mathtype == FC_MT_SCALAR && 
                  datatype == FC_DT_DOUBLE, "should return scalar double");
      fail_unless(assoc == FC_AT_ELEMENT, "mismatch of association");
      // displaced
      fail_unless(rc_displ == FC_SUCCESS, 
		  "failed to get regions on displaced hex mesh");
      rc = fc_getVariableInfo(element_volumes_displ, &numData, &numComp, 
			      &assoc, &mathtype, &datatype);
      fail_unless(numData == numElement, "mismatch of numData");
      fail_unless(numComp == 1 && mathtype == FC_MT_SCALAR && 
                  datatype == FC_DT_DOUBLE, "should return scalar double");
      fail_unless(assoc == FC_AT_ELEMENT, "mismatch of association");
    }

    // test getting mesh volume
    rc = fc_getMeshVolume(mesh, &mesh_volume);
    rc_displ = fc_getDisplacedMeshVolume(mesh, coords_var,
					 &mesh_volume_displ);
    if (i < 4) { // point mesh, line mesh, tri mesh, or quad mesh
      // regular
      fail_unless(rc != FC_SUCCESS, 
                  "should fail to get volumes on 2D or less mesh");
      fail_unless(mesh_volume == -1,
		  "should return mesh volume of -1 by side effect when fails");
      fail_unless(rc_displ != FC_SUCCESS, 
                  "should fail to get displaced volumes on 2D or less mesh");
      fail_unless(mesh_volume_displ == -1,
		  "should return displaced mesh volume of -1 by side effect when fails");
    } else {
      fail_unless(rc == FC_SUCCESS, "failed to get mesh volume on 3D mesh");
      fail_unless(rc_displ == FC_SUCCESS, "failed to get displaced mesh volume on 3D mesh");
    }

    //test getting subset and displaced subset volume
    //setup for FC_AT_WHOLE_MESH subset tests
    rc = fc_createSubset(mesh, "temp subset", FC_AT_WHOLE_MESH, &subset);
    fail_unless(rc == FC_SUCCESS, "abort: failed to create subset");
    fc_addMemberToSubset(subset, 0);

    rc = fc_getSubsetVolume(subset,1, &subset_volume);
    rc_displ = fc_getDisplacedSubsetVolume(subset,1, coords_var,&subset_volume_displ);
    if (i < 4) { // point mesh, line mesh, tri mesh, or quad mesh
      // regular
      fail_unless(rc != FC_SUCCESS, 
                  "should fail to get subset volumes on less than 3D mesh");
      fail_unless(subset_volume == -1,
		  "should return subset volume of -1 by side effect when fails");
      fail_unless(rc_displ != FC_SUCCESS, 
                  "should fail to get displaced subset volumes on less than 3D mesh");
      fail_unless(subset_volume_displ == -1,
		  "should return displaced subset volume of -1 by side effect when fails");
    } else {
      fail_unless(rc == FC_SUCCESS, "failed to get subset volume on 3D mesh");
      fail_unless(rc_displ == FC_SUCCESS,
		  "failed to get displaced subset volume on 3D mesh");
    }

    // test return answers for 3D meshes

    // test values for regular & displaced mesh
    if (i > 3) { // tet mesh, pyramid mesh, prism mesh or hex mesh
      //printf("\nlengths (%d): %.20g", numKnown_lens, known_lens[0]);
      //for (j = 1; j < numKnown_lens; j++)
      //printf(", %.20g", known_lens[j]);
      //printf("\n");
      //printf("areas (%d): %.20g", numKnown_areas, known_areas[0]);
      //for (j = 1; j < numKnown_areas; j++)
      //printf(", %.20g", known_areas[j]);
      //printf("\n");
      //printf("volumes (%d): %.20g", numKnown_vols, known_vols[0]);
      //for (j = 1; j < numKnown_vols; j++)
      //printf(", %.20g", known_vols[j]);
      //printf("\n");
      fc_getVariableDataPtr(edge_lengths, (void*)&edge_data);
      fc_getVariableDataPtr(edge_lengths_displ, (void*)&edge_data_displ);
      for (j = 0; j < numEdge; j++) {
        int foundMatch = 0;
        for (k = 0; k < numKnown_lens; k++) {
          if (FC_VALUE_EQUIV(edge_data[j], known_lens[k], 
                             100*DBL_EPSILON, DBL_MIN)) {
            foundMatch = 1;
            break;
          }
        }
        //printf("%d (%d): %.20g\n", j, numEdge, edge_data[j]);
        fail_unless(foundMatch == 1, "edge length wasn't of known length");
	fail_unless(FC_VALUE_EQUIV(edge_data_displ[j], 2.*known_lens[k],
				   100*DBL_EPSILON, DBL_MIN),
		    "displaced length mismatch");
      }
      fc_getVariableDataPtr(face_areas, (void*)&face_data);
      fc_getVariableDataPtr(face_areas_displ, (void*)&face_data_displ);
      for (j = 0; j < numFace; j++) {
        int foundMatch = 0;
        for (k = 0; k < numKnown_areas; k++) {
          if (FC_VALUE_EQUIV(face_data[j], known_areas[k], 
                             200*DBL_EPSILON, DBL_MIN)) {
            foundMatch = 1;
            break;
          }
        }
        //printf("%d (%d): %.20g\n", j, numFace, face_data[j]);
        fail_unless(foundMatch == 1, "face area wasn't of known area");
        fail_unless(FC_VALUE_EQUIV(face_data_displ[j], 4.*known_areas[k], 
				   200*DBL_EPSILON, DBL_MIN),
		    "displaced area mismatch");
      }
      fc_getVariableDataPtr(element_volumes, (void*)&element_data);
      fc_getVariableDataPtr(element_volumes_displ, (void*)&element_data_displ);
      for (j = 0; j < numElement; j++) {
        int foundMatch = 0;
        for (k = 0; k < numKnown_vols; k++) {
          if (FC_VALUE_EQUIV(element_data[j], known_vols[k], 
                             200*DBL_EPSILON, DBL_MIN)) {
            foundMatch = 1;
            break;
          }
        }
        //printf("%d (%d): %.20g\n", j, numElement, element_data[j]);
        fail_unless(foundMatch == 1, "region volume wasn't of known volume");
	fail_unless(FC_VALUE_EQUIV(element_data_displ[j], 8.*known_vols[k], 
				   200*DBL_EPSILON, DBL_MIN),
		    "displaced volume mismatch");
      }

      //checking result for fc_getMeshVolume
      fail_unless(FC_VALUE_EQUIV(mesh_volume,known_mesh_vol,500*DBL_EPSILON,
				 DBL_MIN),"mesh volume mismatch");

      //checking result for fc_getDisplacedMeshVolume
      fail_unless(FC_VALUE_EQUIV(mesh_volume_displ,known_mesh_vol*8,
				 500*DBL_EPSILON,DBL_MIN),
		  "displaced mesh volume mismatch");

      //checking result for fc_getSubsetVolume for FC_AT_WHOLE_MESH
      fail_unless(FC_VALUE_EQUIV(subset_volume,known_mesh_vol,
				 500*DBL_EPSILON,DBL_MIN),
		  "subset volume mismatch");

      //checking result for fc_getDisplacedSubsetVolume for FC_AT_WHOLE_MESH
      //      printf("\n%d %.20g %.20g\n",i,subset_volume_displ,known_mesh_vol);
      fail_unless(FC_VALUE_EQUIV(subset_volume_displ,8*known_mesh_vol,
				 500*DBL_EPSILON,DBL_MIN),
		  "displaced subset volume mismatch");

      //for fc_getSubsetVolume and fc_getDisplacedSubsetVolume for other associations
      //if it is the well known hex mesh do more rigorous testing
      if (i == 7){
	//145 is in the middle - see check_entity_type in checktopo.c
	int hexID = 145; 

	FC_AssociationType subset_assocs[4] = { FC_AT_VERTEX, FC_AT_EDGE, 
						FC_AT_FACE, FC_AT_ELEMENT };  

	for (k = 0; k < 4; k++){
	  FC_Subset junkSubset, newJunkSubset;	

	  //regular subset
	  rc = fc_createSubset(mesh,"junkSubset",FC_AT_ELEMENT,&junkSubset);
	  fail_unless(rc == FC_SUCCESS, "abort: failed to create subset");
	  rc = fc_addMemberToSubset(junkSubset,hexID);
	  fail_unless(rc == FC_SUCCESS,
		      "abort: failed to add members to subset");
	  rc = fc_copySubsetWithNewAssociation
	    (junkSubset,"newJunkSubset",subset_assocs[k],0,&newJunkSubset);
	  fail_unless(rc == FC_SUCCESS,
		      "abort: failed to change association of subset");

	  //test with newJunkSubset that has the assoc we are testing
	  for (n = 0; n < 2; n++){
	    double vol = -100.0;
	    double vol_displ = -100.0;
	    rc = fc_getSubsetVolume(newJunkSubset,n,&vol);
	    fail_unless(rc == FC_SUCCESS,
		      "abort: failed to get subset volume to check");
	    rc_displ = fc_getDisplacedSubsetVolume(newJunkSubset,n,coords_var,
						   &vol_displ);
	    fail_unless(rc_displ == FC_SUCCESS,
		      "abort: failed to get displaced subset volume to check");

	    if (n == 1){
	      //strict
	      fail_unless(FC_VALUE_EQUIV(vol,known_vols[0],
					 DBL_EPSILON,DBL_MIN),
			  "subset volume mismatch");
	      fail_unless(FC_VALUE_EQUIV(vol_displ,8*known_vols[0],
					 DBL_EPSILON,DBL_MIN),
			  "displaced subset volume mismatch");
	    } else {
	      //lenient
	      switch (k){
	      case 0:
		//its the element + the 26 elements around it
		fail_unless(FC_VALUE_EQUIV(vol,27*known_vols[0],
					 10*DBL_EPSILON,DBL_MIN),
			  "subset volume mismatch");
		fail_unless(FC_VALUE_EQUIV(vol_displ,8*27*known_vols[0],
					 10*DBL_EPSILON,DBL_MIN),
			  "displaced subset volume mismatch");
		break;
	      case 1:
		// its the 27 above minus the 8 very corner ones
		fail_unless(FC_VALUE_EQUIV(vol,19*known_vols[0],
					 10*DBL_EPSILON,DBL_MIN),
			  "subset volume mismatch");
		fail_unless(FC_VALUE_EQUIV(vol_displ,8*19*known_vols[0],
					 10*DBL_EPSILON,DBL_MIN),
			  "displaced subset volume mismatch");
		break;
	      case 2:
		//its the element and the 6 that it shares faces with
		fail_unless(FC_VALUE_EQUIV(vol,7*known_vols[0],
					 10*DBL_EPSILON,DBL_MIN),
			  "subset volume mismatch");
		fail_unless(FC_VALUE_EQUIV(vol_displ,8*7*known_vols[0],
					 10*DBL_EPSILON,DBL_MIN),
			  "displaced subset volume mismatch");
		break;
	      case 3:
		//its the element itself
		fail_unless(FC_VALUE_EQUIV(vol,known_vols[0],
					   DBL_EPSILON,DBL_MIN),
			    "subset volume mismatch");

		fail_unless(FC_VALUE_EQUIV(vol_displ,8*known_vols[0],
					   DBL_EPSILON,DBL_MIN),
			    "displaced subset volume mismatch");
		break;
	      }
	    }
	  }

	  fc_deleteSubset(junkSubset);
	  fc_deleteSubset(newJunkSubset);
	}
      }
    }

    // Save var that can't be interpreted as displ for tests below
    if (i > 0) {
      rc = fc_copyVariable(edge_lengths, mesh, "non-coordsVar", &nonDisplVar);
      fail_unless(rc == FC_SUCCESS, "abort: copy failed");
    }
    else
      nonDisplVar = badVar; 

    // cleanup -- except subset -- using in tests below
    fc_deleteVariable(edge_lengths);
    fc_deleteVariable(face_areas);
    fc_deleteVariable(element_volumes);
    fc_deleteVariable(edge_lengths_displ);
    fc_deleteVariable(face_areas_displ);
    fc_deleteVariable(element_volumes_displ);

    // --- test bad arguments

    // fc_getEdgeLengths()
    rc = fc_getEdgeLengths(badMesh, &edge_lengths);
    fail_unless(rc != FC_SUCCESS, "should fail if bad mesh");
    fail_unless(FC_HANDLE_EQUIV(edge_lengths, FC_NULL_VARIABLE),
                "failure should return null handle");
    rc = fc_getEdgeLengths(mesh, NULL);
    fail_unless(rc != FC_SUCCESS, "should fail if bad mesh");

    // fc_getFaceAreas()
    rc = fc_getFaceAreas(badMesh, &face_areas);
    fail_unless(rc != FC_SUCCESS, "should fail if bad mesh");
    fail_unless(FC_HANDLE_EQUIV(face_areas, FC_NULL_VARIABLE),
                "failure should return null handle");
    rc = fc_getFaceAreas(mesh, NULL);
    fail_unless(rc != FC_SUCCESS, "should fail if bad mesh");

    // fc_getRegionVolumes()
    rc = fc_getRegionVolumes(badMesh, &element_volumes);
    fail_unless(rc != FC_SUCCESS, "should fail if bad mesh");
    fail_unless(FC_HANDLE_EQUIV(element_volumes, FC_NULL_VARIABLE),
                "failure should return null handle");
    rc = fc_getRegionVolumes(mesh, NULL);
    fail_unless(rc != FC_SUCCESS, "should fail if bad mesh");    

    // fc_getDisplacedEdgeLengths()
    rc = fc_getDisplacedEdgeLengths(badMesh, coords_var, &edge_lengths_displ);
    fail_unless(rc != FC_SUCCESS, "should fail if bad mesh");
    fail_unless(FC_HANDLE_EQUIV(edge_lengths, FC_NULL_VARIABLE),
                "failure should return null handle");
    rc = fc_getDisplacedEdgeLengths(mesh, badVar, &edge_lengths);
    fail_unless(rc != FC_SUCCESS, "should fail if bad displ var");
    fail_unless(FC_HANDLE_EQUIV(edge_lengths, FC_NULL_VARIABLE),
                "failure should return null handle");
    rc = fc_getDisplacedEdgeLengths(mesh, nonDisplVar, &edge_lengths_displ);
    fail_unless(rc != FC_SUCCESS, "should fail if inappropriate displ var");
    fail_unless(FC_HANDLE_EQUIV(edge_lengths, FC_NULL_VARIABLE),
                "failure should return null handle");
    rc = fc_getDisplacedEdgeLengths(mesh, coords_var, NULL);
    fail_unless(rc != FC_SUCCESS, "should fail if bad mesh");

    // fc_getDisplacedFaceAreas()
    rc = fc_getDisplacedFaceAreas(badMesh, coords_var, &face_areas);
    fail_unless(rc != FC_SUCCESS, "should fail if bad mesh");
    fail_unless(FC_HANDLE_EQUIV(face_areas, FC_NULL_VARIABLE),
                "failure should return null handle");
     rc = fc_getDisplacedFaceAreas(mesh, badVar, &face_areas);
    fail_unless(rc != FC_SUCCESS, "should fail if bad displ var");
    fail_unless(FC_HANDLE_EQUIV(face_areas, FC_NULL_VARIABLE),
                "failure should return null handle");
    rc = fc_getDisplacedFaceAreas(mesh, nonDisplVar, &face_areas);
    fail_unless(rc != FC_SUCCESS, "should fail if inappropriate displ var");
    fail_unless(FC_HANDLE_EQUIV(face_areas, FC_NULL_VARIABLE),
                "failure should return null handle");
    rc = fc_getDisplacedFaceAreas(mesh, coords_var, NULL);
    fail_unless(rc != FC_SUCCESS, "should fail if bad mesh");

    // fc_getDisplacedRegionVolumes()
    rc = fc_getDisplacedRegionVolumes(badMesh, coords_var, &element_volumes);
    fail_unless(rc != FC_SUCCESS, "should fail if bad mesh");
    fail_unless(FC_HANDLE_EQUIV(element_volumes, FC_NULL_VARIABLE),
                "failure should return null handle");
    rc = fc_getDisplacedRegionVolumes(mesh, badVar, &element_volumes);
    fail_unless(rc != FC_SUCCESS, "should fail if bad var");
    fail_unless(FC_HANDLE_EQUIV(element_volumes, FC_NULL_VARIABLE),
                "failure should return null handle");
    rc = fc_getDisplacedRegionVolumes(mesh, nonDisplVar, &element_volumes);
    fail_unless(rc != FC_SUCCESS, "should fail if inappropriate displ var");
    fail_unless(FC_HANDLE_EQUIV(element_volumes, FC_NULL_VARIABLE),
                "failure should return null handle");
    rc = fc_getDisplacedRegionVolumes(mesh, coords_var, NULL);
    fail_unless(rc != FC_SUCCESS, "should fail if bad mesh");    

    //fc_getMeshVolume()
    rc = fc_getMeshVolume(badMesh, &mesh_volume);
    fail_unless(rc != FC_SUCCESS, "should fail if bad mesh");
    fail_unless(FC_DBL_EQUIV(mesh_volume,-1.0),
                "failure should return volume == -1 ");
    rc = fc_getMeshVolume(mesh, NULL);
    fail_unless(rc != FC_SUCCESS, "should fail if passed NULL volume");    
    fail_unless(FC_DBL_EQUIV(mesh_volume,-1.0),
                "failure should return volume == -1 ");

    //  fc_getDisplacedMeshVolume()
    rc = fc_getDisplacedMeshVolume(badMesh, coords_var,&mesh_volume);
    fail_unless(rc != FC_SUCCESS, "should fail if bad mesh");
    fail_unless(FC_DBL_EQUIV(mesh_volume,-1.0),
		"failure should return volume == -1 ");
    rc = fc_getDisplacedMeshVolume(mesh, coords_var,NULL);
    fail_unless(rc != FC_SUCCESS, "should fail if passed NULL volume");    
    fail_unless(FC_DBL_EQUIV(mesh_volume,-1.0),
                "failure should return volume == -1 ");
    rc = fc_getDisplacedMeshVolume(mesh, nonDisplVar, &mesh_volume);
    fail_unless(rc != FC_SUCCESS, "should fail if inappropriate displ var");
    fail_unless(FC_DBL_EQUIV(mesh_volume,-1.0),
		"failure should return volume == -1 ");

    //fc_getSubsetVolume()
    rc = fc_getSubsetVolume(badSubset, 1, &subset_volume);
    fail_unless(rc != FC_SUCCESS, "should fail if bad subset");
    fail_unless(FC_DBL_EQUIV(subset_volume,-1.0),
                "failure should return volume == -1 ");
    //just checking the most likely bad doStrict flags here
    rc = fc_getSubsetVolume(subset, -1, &subset_volume);
    fail_unless(rc != FC_SUCCESS, "should fail if passed invalid dostrict flag");
    fail_unless(FC_DBL_EQUIV(subset_volume,-1.0),
                "failure should return volume == -1 ");
    rc = fc_getSubsetVolume(subset, 2, &subset_volume);
    fail_unless(rc != FC_SUCCESS, "should fail if passed invalid dostrict flag");
    fail_unless(FC_DBL_EQUIV(subset_volume,-1.0),
                "failure should return volume == -1 ");
    rc = fc_getSubsetVolume(subset, 1, NULL);
    fail_unless(rc != FC_SUCCESS, "should fail if passed NULL volume");    

    //fc_getDisplacedSubsetVolume()
    rc = fc_getDisplacedSubsetVolume(badSubset, 1, coords_var, &subset_volume);
    fail_unless(rc != FC_SUCCESS, "should fail if bad subset");
    fail_unless(FC_DBL_EQUIV(subset_volume,-1.0),
                "failure should return volume == -1 ");
    //just checking the most likely bad doStrict flags here
    rc = fc_getDisplacedSubsetVolume(subset, -1, coords_var, &subset_volume);
    fail_unless(rc != FC_SUCCESS, "should fail if passed invalid dostrict flag");
    fail_unless(FC_DBL_EQUIV(subset_volume,-1.0),
                "failure should return volume == -1 ");
    rc = fc_getDisplacedSubsetVolume(subset, 2, coords_var, &subset_volume);
    fail_unless(rc != FC_SUCCESS, "should fail if passed invalid dostrict flag");
    fail_unless(FC_DBL_EQUIV(subset_volume,-1.0),
                "failure should return volume == -1 ");
    rc = fc_getDisplacedSubsetVolume(subset, 1, nonDisplVar, &subset_volume);
    fail_unless(rc != FC_SUCCESS, "should fail if inappropriate displ var");
    fail_unless(FC_DBL_EQUIV(subset_volume,-1.0),
                "failure should return volume == -1 ");
    rc = fc_getDisplacedSubsetVolume(subset, 1, coords_var, NULL);
    fail_unless(rc != FC_SUCCESS, "should fail if passed NULL volume");    

    // cleanup
    fc_deleteMesh(mesh);
  }

  rc = fc_deleteDataset(dataset);
  fail_unless(rc == FC_SUCCESS, "failed to close dataset at end of test");
}
END_TEST


START_TEST(subset_area){
  //testing subsetarea here - I can't find a good place to put it in above
  FC_ReturnCode rc;
  FC_Dataset dataset;
  FC_Mesh mesh, *returnMeshes;
  int numReturnMeshes;
  FC_Subset emptyf, wrongassoc,skin;
  FC_Subset badSubset={999,999};
  double area, tot;
  double A = 10., B = 8., C = 5.;  // from data/mesh_generator.c
  
  fc_loadDataset("../data/gen_multimesh.ex2", &dataset);
  rc = fc_getMeshByName(dataset, "hex mesh",&numReturnMeshes,&returnMeshes);
  fail_unless(rc == FC_SUCCESS, "aborted - cant get mesh by name");
  fail_unless(numReturnMeshes == 1, "failed to find unique mesh by name");
  mesh = returnMeshes[0];
  free(returnMeshes);
  
  rc = fc_createSubset(mesh,"emptyf",FC_AT_FACE,&emptyf);
  fc_exitIfErrorPrintf(rc, "failed to make subset");
  rc = fc_createSubset(mesh,"wrongassoc",FC_AT_ELEMENT,&wrongassoc);
  fc_exitIfErrorPrintf(rc, "failed to make subset");
  rc = fc_addMemberToSubset(wrongassoc,1);
  fc_exitIfErrorPrintf(rc, "failed to addmemeber to subset");

  rc = fc_getSubsetArea(emptyf, &area);
  fail_unless(rc == FC_SUCCESS, "cant get subset area");
  fail_unless(FC_DBL_EQUIV(0.0,area), "bad result for area");
  fc_deleteSubset(emptyf);
  
  rc = fc_getMeshSkin(mesh,&skin);
  fail_unless(rc == FC_SUCCESS, "cant get subset skin");    
  rc = fc_getSubsetArea(skin, &area);
  fail_unless(rc == FC_SUCCESS, "cant get subset area");
  
  //calc the area from the known_areas;
  tot = 2*A*C+2*B*C+2*A*B;
  fail_unless(FC_VALUE_EQUIV(tot,area, 2*DBL_EPSILON, DBL_MIN), 
	      "bad result for area");

  //bad args
  rc = fc_getSubsetArea(badSubset, &area);
  fail_unless(rc != FC_SUCCESS, "should fail if bad subset");
  fail_unless(FC_DBL_EQUIV(area,-1.0),
	      "failure should return area == -1 ");
  rc = fc_getSubsetArea(FC_NULL_SUBSET, &area);
  fail_unless(rc != FC_SUCCESS, "should fail if bad subset");
  fail_unless(FC_DBL_EQUIV(area,-1.0),
	      "failure should return area == -1 ");
  rc = fc_getSubsetArea(wrongassoc, &area);
  fail_unless(rc != FC_SUCCESS, "should fail for bad assoc");    
  fail_unless(FC_DBL_EQUIV(area,-1.0),
	      "failure should return area == -1 ");
  rc = fc_getSubsetArea(skin, NULL);
  fail_unless(rc != FC_SUCCESS, "should fail if passed NULL area");    

  //clean up
  fc_deleteSubset(skin);
  fc_deleteSubset(emptyf);
  fc_deleteSubset(wrongassoc);
  fc_deleteDataset(dataset);
}
END_TEST



START_TEST(polygon_normal)
{
  int i, j;
  // planar quads
  FC_Coords points1[4] = { { 0., 0., 0. }, { 1., 0., 0. }, 
			   { 1., 1., 0. }, { 0., 1., 0. } };
  FC_Coords points2[4] = { { 0., 0., 0. }, { 0., 1., 0. }, 
			   { 0., 1., 1. }, { 0., 0., 1. } };
  FC_Coords points3[4] = { { 0., 0., 0. }, { 0., 0., 1. }, 
			   { 1., 0., 1. }, { 1., 0., 0. } };
  // tri with extra colinear point
  FC_Coords points4[4] = { { 1., 0., 0. }, { 0., 1., 0. }, 
			   { 0., 0., 1. }, { 0.5, 0., 0.5 } };
  // non planar quad (wrapping the first 2 points to make testing easier)
  FC_Coords points5[6] = { { 1., 0., 0. }, { 0., 1., 0. }, 
			   { 0., 0., 1. }, { 0.501, 0.001, 0.501 } };
  FC_Coords normal;
  
  for (i = 3; i <= 4; i++) {
    _fc_calcSurfaceNormal(i, points1, &normal);
    fail_unless(normal[0] == 0. && normal[1] == 0, "mismatch in normal");
    fail_unless(normal[2] == 1., "mismatch in normal");
    _fc_calcSurfaceNormal(i, points2, &normal);
    fail_unless(normal[1] == 0. && normal[2] == 0, "mismatch in normal");
    fail_unless(normal[0] == 1., "mismatch in normal");
    _fc_calcSurfaceNormal(i, points3, &normal);
    fail_unless(normal[0] == 0. && normal[2] == 0, "mismatch in normal");
    fail_unless(normal[1] == 1., "mismatch in normal");
    _fc_calcSurfaceNormal(i, points4, &normal);
    for (j = 0; j < 3; j++)
      fail_unless(FC_DBL_EQUIV(normal[j], 1./sqrt(3)), "mismatch in normal");
  }

   // compute answer as average of triangles
  _fc_calcSurfaceNormal(4, points5, &normal);
  for (i = 0; i < 3; i++) {
    fail_unless(!FC_DBL_EQUIV(normal[i], 1./sqrt(3)),
		"should not be same as planar normal");
    fail_unless(FC_VALUE_EQUIV(normal[i], 1./sqrt(3), 1e-2, DBL_MIN),
		"should be close to normal of first 3 points");
  }
}
END_TEST

// test calculation of a displaced_surface
START_TEST(surface_normals)
{
  FC_ReturnCode rc;
  int i, j, k;
  FC_Dataset dataset;
  FC_Mesh mesh, *returnMeshes;
  int numReturnMeshes;
  FC_Variable coords;
  int numVertex, numElement, dim, temp_numDataPoint, temp_numComp;
  int numVertPerElem;
  FC_ElementType elemType;
  FC_AssociationType temp_assoc;
  FC_MathType temp_mathtype;
  FC_DataType temp_datatype;
  double *coord_data;
  double *norm_data, *displaced_norm_data;
  FC_Variable norms, displacedNorms;
  FC_Coords temp_points[4], normal;
  int* conns;

  // setup
  fc_loadDataset(exo_file, &dataset);
  rc = fc_getMeshByName(dataset, mesh_name,&numReturnMeshes,&returnMeshes);
  fail_unless(rc == FC_SUCCESS, "aborted - cant get mesh by name");
  fail_unless(numReturnMeshes == 1, "failed to find unique mesh by name");
  mesh = returnMeshes[0];
  free(returnMeshes);

  fc_getMeshInfo(mesh, NULL, NULL, NULL, &numElement, &elemType);
  rc = fc_getMeshCoordsAsVariable(mesh, &coords);
  fail_unless(rc == FC_SUCCESS, "Test aborted: couldn't find coord field");
  fc_getVariableInfo(coords, &numVertex, &dim, NULL, NULL, NULL);
  fail_unless(dim > 1, "Test aborted: expected vector coords");
  rc = fc_getVariableDataPtr(coords, (void**)&coord_data);
  fail_unless(rc == FC_SUCCESS, "Test aborted: could not get coord data");
  rc = fc_getMeshElementConnsPtr(mesh, &conns);
  fail_unless(rc == FC_SUCCESS, "Aborted: could nnot get elem conns");
  numVertPerElem = fc_getElementTypeNumVertex(elemType);

  // --- test fc_getSurfaceNormals()

  // check doing it
  rc = fc_getSurfaceNormals(mesh, &norms);
  fail_unless(rc == FC_SUCCESS, "failed to get surface normals");

  // check metadata
  rc = fc_getVariableInfo(norms, &temp_numDataPoint, &temp_numComp,
			  &temp_assoc, &temp_mathtype, &temp_datatype);
  fail_unless(rc == FC_SUCCESS, "failed to describe normal var");
  fail_unless(temp_numDataPoint == numElement, "numDataPoint mismatch");
  fail_unless(temp_numComp == 3, "numComp should = 3");
  fail_unless(temp_assoc == FC_AT_ELEMENT, "assoc mismatch");
  fail_unless(temp_mathtype == FC_MT_VECTOR, "mathtype should be vector");
  fail_unless(temp_datatype == FC_DT_DOUBLE, "mathtype should be double");

  // check values
  rc = fc_getVariableDataPtr(norms, (void**)&norm_data);
  fail_unless(rc == FC_SUCCESS, "failed to get norm data");
  for (i = 0; i < numElement; i++) {
    for (j = 0; j < numVertPerElem; j++)
      for (k = 0; k < dim; k++)
	temp_points[j][k] = coord_data[conns[i*numVertPerElem+j]*dim + k];
    _fc_calcSurfaceNormal(numVertPerElem, temp_points, &normal);
    for (j = 0; j < 3; j++)
      fail_unless(norm_data[i*3+j] == normal[j], "mismatch of normal value");
  }

  // --- test fc_getDisplacedSurfaceNormals()

  // check doing it
  rc = fc_getDisplacedSurfaceNormals(mesh, coords, &displacedNorms);
  fail_unless(rc == FC_SUCCESS, "failed to get displaced surface normals");

  // check metadata
  rc = fc_getVariableInfo(displacedNorms, &temp_numDataPoint, &temp_numComp,
			  &temp_assoc, &temp_mathtype, &temp_datatype);
  fail_unless(rc == FC_SUCCESS, "failed to describe normal var");
  fail_unless(temp_numDataPoint == numElement, "numDataPoint mismatch");
  fail_unless(temp_numComp == 3, "numComp should = 3");
  fail_unless(temp_assoc == FC_AT_ELEMENT, "assoc mismatch");
  fail_unless(temp_mathtype == FC_MT_VECTOR, "mathtype should be vector");
  fail_unless(temp_datatype == FC_DT_DOUBLE, "mathtype should be double");

  // check values -- displacement should not have affected normals
  rc = fc_getVariableDataPtr(displacedNorms, (void**)&displaced_norm_data);
  fail_unless(rc == FC_SUCCESS, "failed to get displaced norm data");
  for (i = 0; i < numElement*3; i++) 
    fail_unless(displaced_norm_data[i] == norm_data[i],
		"mismatch of displaced normal value");

  fc_deleteDataset(dataset);
}
END_TEST

// test calculation of extract deform via a reference vertex
START_TEST(deform_vertex)
{
  FC_ReturnCode rc;
  int i, j;
  FC_Dataset dataset;
  FC_Sequence sequence, *returnSequences;
  FC_Mesh mesh, *returnMeshes;
  int numReturnMeshes, numReturnSequences;
  FC_Variable* displSeqVar, *temp_seqVar;
  int numStep, numVertex, numDim, temp_numDataPoint, temp_numComp;
  FC_AssociationType temp_assoc;
  FC_MathType temp_mathtype;
  FC_DataType temp_datatype;
  double *trans, **deforms, *displs, *temp_data;

  // setup
  fc_loadDataset("../data/gen_hex_seq.ex2", &dataset);
  rc = fc_getSequenceByName(dataset, "time", &numReturnSequences, 
			    &returnSequences);
  fail_unless(rc == FC_SUCCESS, "Test aborted: can't get sequence by name");
  fail_unless(numReturnSequences == 1, "Test aborted: no unique sequence");
  sequence = returnSequences[0];
  free(returnSequences);

  fc_getSequenceNumStep(sequence, &numStep);
  rc = fc_getMeshByName(dataset, "hex mesh", &numReturnMeshes,&returnMeshes);
  fail_unless(rc == FC_SUCCESS, "aborted - cant get mesh by name");
  fail_unless(numReturnMeshes == 1, "failed to find unique mesh by name");
  mesh = returnMeshes[0];
  free(returnMeshes);

  fail_unless(rc == FC_SUCCESS, "Test aborted: mesh not found");
  fc_getMeshNumVertex(mesh, &numVertex);
  fc_getMeshDim(mesh, &numDim);

  // make space for calc arrays
  trans = (double*)malloc(numDim*numVertex*sizeof(double));
  deforms = (double**)malloc(numStep*sizeof(double*));
  displs = (double*)malloc(numDim*numVertex*sizeof(double));
  fail_unless(trans && deforms && displs, "malloc error");
  for (i = 0; i < numStep; i++) {
    deforms[i] = (double*)malloc(numDim*numVertex*sizeof(double));
    fail_unless(deforms[i] != NULL, "malloc error");
  }

  // make arbitrary trans data (the rate of change)
  for (i = 0; i < numVertex; i++) {
    trans[i*numDim] = 0.1;
    trans[i*numDim + 1] = 0.2;
    trans[i*numDim + 2] = 0.3;
  }
  // make arbitrary deform data - none at t0, then move 1st 10 vertices
  for (i = 0; i < numVertex*numDim; i++)
    deforms[0][i] = 0.;
  for (i = 1; i < numStep; i++) {
    for (j = 0; j < 10; j++) {
      deforms[i][j*numDim] =     0.01*i;
      deforms[i][j*numDim + 1] = 0.03*i;
      deforms[i][j*numDim + 2] = 0.05*i;
    }
    for (j = 10*numDim; j < numVertex*numDim; j++)
      deforms[i][j] = 0.;
  }

  // Make displacement seq variable = trans + deform 
  rc = fc_createSeqVariable(mesh, sequence, "diplacement", &numStep, 
			    &displSeqVar);
  fail_unless(rc == FC_SUCCESS, "failed to create seq variable");
  for (i = 0; i < numStep; i++) {
    for (j = 0; j < numDim*numVertex; j++)
      displs[j] = i*trans[j] + deforms[i][j];
    rc = fc_setVariableData(displSeqVar[i], numVertex, numDim, FC_AT_VERTEX, 
			    FC_MT_VECTOR, FC_DT_DOUBLE, (void*)displs);
    fail_unless(rc == FC_SUCCESS, "failed to set seq var data");
  }
  free(trans);
  free(displs);

  // Do it & check values

  // Centroid
  rc = fc_calcMeshDeformationUsingRefVertex(mesh, 11, numStep, displSeqVar,
					 &temp_seqVar);
  fail_unless(rc == FC_SUCCESS, "Failed to get deforms via centroid");
  for (i = 0; i < numStep; i++) {
    rc = fc_getVariableInfo(temp_seqVar[i], &temp_numDataPoint, &temp_numComp, 
			    &temp_assoc, &temp_mathtype, &temp_datatype);
    fail_unless(rc == FC_SUCCESS, "failed to describe  var");
    fail_unless(temp_numDataPoint == numVertex, "numDataPoint mismatch");
    fail_unless(temp_numComp == numDim, "numComp should = 3");
    fail_unless(temp_assoc == FC_AT_VERTEX, "assoc mismatch");
    fail_unless(temp_mathtype == FC_MT_VECTOR, "mathtype should be vector");
    fail_unless(temp_datatype == FC_DT_DOUBLE, "mathtype should be double");
    rc = fc_getVariableDataPtr(temp_seqVar[i], (void**)&temp_data);
    fail_unless(rc == FC_SUCCESS, "failed to get data");
    for (j = 0; j < numDim*numVertex; j++) {
      //printf("%d:%d  %.16f %.16f \n", i, j, temp_data[j], deforms[i][j]);
      fail_unless(FC_VALUE_EQUIV(temp_data[j], deforms[i][j],
		  500*DBL_EPSILON, DBL_MIN), "did not get expected data");
    }
  }
  fc_deleteSeqVariable(numStep, temp_seqVar);
  free(temp_seqVar);
  fc_deleteSeqVariable(numStep, displSeqVar);
  free(displSeqVar);

  for (i = 0; i < numStep; i++)
    free(deforms[i]);
  free(deforms);

  rc = fc_deleteDataset(dataset);
  fail_unless(rc == FC_SUCCESS, "failed to close dataset at end of test");
}
END_TEST

// test calculation of extract deform via centroid
START_TEST(deform_centroid)
{
  FC_ReturnCode rc;
  int i, j;
  FC_Dataset dataset;
  FC_Sequence sequence, *returnSequences;
  FC_Mesh mesh, *returnMeshes;
  int numReturnMeshes, numReturnSequences;
  FC_Variable* displSeqVar, *temp_seqVar;
  int numStep, numVertex, numDim, temp_numDataPoint, temp_numComp;
  FC_AssociationType temp_assoc;
  FC_MathType temp_mathtype;
  FC_DataType temp_datatype;
  double *trans, **deforms, *displs, *temp_data;

  // setup
  fc_loadDataset("../data/gen_hex_seq.ex2", &dataset);
  rc = fc_getSequenceByName(dataset, "time", &numReturnSequences, 
			    &returnSequences);
  fail_unless(rc == FC_SUCCESS, "Test aborted: can't get sequence by name");
  fail_unless(numReturnSequences == 1, "Test aborted: no unique sequence");
  sequence = returnSequences[0];
  free(returnSequences);

  fc_getSequenceNumStep(sequence, &numStep);
  rc = fc_getMeshByName(dataset, "hex mesh", &numReturnMeshes,&returnMeshes);
  fail_unless(rc == FC_SUCCESS, "aborted - cant get mesh by name");
  fail_unless(numReturnMeshes == 1, "failed to find unique mesh by name");
  mesh = returnMeshes[0];
  free(returnMeshes);

  fail_unless(rc == FC_SUCCESS, "Test aborted: mesh not found");
  fc_getMeshNumVertex(mesh, &numVertex);
  fc_getMeshDim(mesh, &numDim);

  // make space for calc arrays
  trans = (double*)malloc(numDim*numVertex*sizeof(double));
  deforms = (double**)malloc(numStep*sizeof(double*));
  displs = (double*)malloc(numDim*numVertex*sizeof(double));
  fail_unless(trans && deforms && displs, "malloc error");
  for (i = 0; i < numStep; i++) {
    deforms[i] = (double*)malloc(numDim*numVertex*sizeof(double));
    fail_unless(deforms[i] != NULL, "malloc error");
  }

  // make arbitrary trans data (the rate of change)
  for (i = 0; i < numVertex; i++) {
    trans[i*numDim] = 0.1;
    trans[i*numDim + 1] = 0.2;
    trans[i*numDim + 2] = 0.3;
  }
  // make arbitrary deform data - move opposite corners (to keep centroid the 
  // same)
  for (i = 0; i < numStep; i++)
    for (j = 0; j < numVertex*numDim; j++)
      deforms[i][j] = 0.;
  for (i = 1; i < numStep; i++) {
    deforms[i][0] = 0.1*i;
    deforms[i][1] = 0.1*i;
    deforms[i][2] = 0.1*i;
    deforms[i][3273] = -0.1*i;
    deforms[i][3274] = -0.1*i;
    deforms[i][3275] = -0.1*i;
  }

  // Make displacement seq variable = trans + deform 
  rc = fc_createSeqVariable(mesh, sequence, "diplacement", &numStep, 
			    &displSeqVar);
  fail_unless(rc == FC_SUCCESS, "failed to create seq variable");
  for (i = 0; i < numStep; i++) {
    for (j = 0; j < numDim*numVertex; j++)
      displs[j] = i*trans[j] + deforms[i][j];
    rc = fc_setVariableData(displSeqVar[i], numVertex, numDim, FC_AT_VERTEX, 
			    FC_MT_VECTOR, FC_DT_DOUBLE, (void*)displs);
    fail_unless(rc == FC_SUCCESS, "failed to set seq var data");
  }
  free(trans);
  free(displs);

  // Do it & check values

  // Centroid
  rc = fc_calcMeshDeformationUsingCentroid(mesh, numStep, displSeqVar,
					 &temp_seqVar);
  fail_unless(rc == FC_SUCCESS, "Failed to get deforms via centroid");
  for (i = 0; i < numStep; i++) {
    rc = fc_getVariableInfo(temp_seqVar[i], &temp_numDataPoint, &temp_numComp, 
			    &temp_assoc, &temp_mathtype, &temp_datatype);
    fail_unless(rc == FC_SUCCESS, "failed to describe  var");
    fail_unless(temp_numDataPoint == numVertex, "numDataPoint mismatch");
    fail_unless(temp_numComp == numDim, "numComp should = 3");
    fail_unless(temp_assoc == FC_AT_VERTEX, "assoc mismatch");
    fail_unless(temp_mathtype == FC_MT_VECTOR, "mathtype should be vector");
    fail_unless(temp_datatype == FC_DT_DOUBLE, "mathtype should be double");
    rc = fc_getVariableDataPtr(temp_seqVar[i], (void**)&temp_data);
    fail_unless(rc == FC_SUCCESS, "failed to get data");
    for (j = 0; j < numDim*numVertex; j++) {
      //printf("%d:%d  %.16f %.16f \n", i, j, temp_data[j], deforms[i][j]);
      // IF this test fails, try loosening the FC_VALUE_EQUIV() test within
      // fc_calcMeshDeformationUsingCentroid()
      // FIX?? don't know why had to make epsilon this big ??
      fail_unless(FC_VALUE_EQUIV(temp_data[j], deforms[i][j],
		  2000*DBL_EPSILON, DBL_MIN), "did not get expected data");
     }
  }
  fc_deleteSeqVariable(numStep, temp_seqVar);
  free(temp_seqVar);
  fc_deleteSeqVariable(numStep, displSeqVar);
  free(displSeqVar);

  for (i = 0; i < numStep; i++)
    free(deforms[i]);
  free(deforms);

  rc = fc_deleteDataset(dataset);
  fail_unless(rc == FC_SUCCESS, "failed to close dataset at end of test");
}
END_TEST

// Assuming that if works in 3D will work in all Ds
// FIX - the o.k. test at the end - should really test that
// the results are correct, not just that that don't fail
START_TEST(mesh_proximity)
{
  FC_ReturnCode rc;
  int i, j, k, m, n;
  int numElemPerSide[3] = { 3, 2, 2 };
  int numVerts[3] = { 4*4*4, 3*3*3, 3*3*3 };
  int numDim = 3;
  FC_Dataset dataset;
  int numMesh = 3;
  FC_Mesh meshes[3], badMesh = { 999, 999 };
  FC_Variable displs[3], badVar = { 999, 999 };
  FC_Subset subsets[3][5], badSubset = { 999, 999 };
  // the first mesh is 3x3 while the 2 & 3rd mesh are 2x2
  // the second mesh is right next to the first (x direction)
  // third mesh had gap of 0.5 between it and the first mesh (x direction)
  FC_Coords lowers[3] = { { 0.,  0.,  0. }, 
			  { 3.0, 0.,  0.  }, 
			  { 3.5, 0.,  0.  } };
  FC_Coords uppers[3] = { { 3.,  3.,  3.  }, 
			  { 5.0, 2.0, 2.0 }, 
			  { 5.5, 2.0, 2.0 } };
  char meshNames[3][10] = { "mesh0", "mesh1", "mesh2" };
  // the subsets are the elements of the meshes that are closest together
  // choosen so results will be about the same.
  int numElemPerSubset[3] = { 3*3, 2*2, 2*2 };
  int elemIDsPerSubset[3][9] = { { 2, 5, 8, 11, 14, 17, 20, 23, 26 },
				 { 0, 2, 4, 6 },
				 { 0, 2, 4, 6 } };
  int numAssocType = 5;
  FC_AssociationType assocs[5] = { FC_AT_WHOLE_MESH, FC_AT_ELEMENT, FC_AT_FACE,
				   FC_AT_EDGE, FC_AT_VERTEX };
  char subsetNames[5][20] = { "whole subset", "elem subset", "face subset",
			      "edge subset", "vert subset" };
  int numCase = 6;
  int meshPairs[6][2] = { { 0, 1 },
			  { 0, 1 },
			  { 0, 2 },
			  { 0, 2 },
			  { 0, 2 },
			  { 0, 2 } };
  double min_distances[6] = { 10.,     // all pairs
			      0.,      // find overlapping verts
			      0.4999,  // should find none
			      0.5,     // just barely get closest verts
			      1.118,   // still just get closest verts
                              1.1181 };  // get some verts in adjacent faces
                                         // sqrt(1^2 + 0.5^2) = 1.1180339
  int numPairs[6] = { 64*27, 9, 0, 9, 9, 39 };
  int numSubsetPairs[6] = { 32*18, 9, 0, 9, 9, 39 }; // the same execpt first one
  int vertIDs1[6][39] = { { }, //too many, so done in situ below
			    { 3, 7, 11, 19, 23, 27, 35, 39, 43 },
			    { }, // none
			    { 3, 7, 11, 19, 23, 27, 35, 39, 43 },
			    { 3, 7, 11, 19, 23, 27, 35, 39, 43 },
			    { 3, 3, 3, 7, 7, 7, 7, 11, 11, 11, 15, 
			      19, 19, 19, 19, 23, 23, 23, 23, 23,  
			      27, 27, 27, 27, 31,
			      35, 35, 35, 39, 39, 39, 39, 43, 43, 43, 47,
			      51, 55, 59 }
  };
  int vertIDs2[6][39] = { { }, // too many, so done in situ below
			    { 0, 3, 6, 9, 12, 15, 18, 21, 24 },
			    { }, // none
			    { 0, 3, 6, 9, 12, 15, 18, 21, 24 },
			    { 0, 3, 6, 9, 12, 15, 18, 21, 24 },
			    { 0, 3, 9, 0, 3, 6, 12, 3, 6, 15, 6,
			      0, 9, 12, 18, 3, 9, 12, 15, 21,
			      6, 12, 15, 24, 15,
			      9, 18, 21, 12, 18, 21, 24, 15, 21, 24, 24,
			      18, 21, 24 }
  };
  int temp_numPair, *temp_vertIDs1, *temp_vertIDs2;
  double data_buf[4*4*4*3];
  int numMembers[3];
  int *memberIDs[3];

  // Make displacement data - the same for all meshes, all verts
  for (i = 0; i < numVerts[0]*numDim; i++)
    data_buf[i] = 1.0;

  // Make the data 
  rc = fc_createDataset("dataset", &dataset);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create dataset");
  for (i = 0; i < numMesh; i++) {
    rc = fc_createSimpleHexMesh(dataset, meshNames[i], numElemPerSide[i], 
				numElemPerSide[i], numElemPerSide[i], 
				lowers[i], uppers[i], &meshes[i]);
    fail_unless(rc == FC_SUCCESS, "abort: failed to create simple mesh");
    fc_createVariable(meshes[i], "displ", &displs[i]);
    rc = fc_setVariableData(displs[i], numVerts[i], numDim, FC_AT_VERTEX,
			    FC_MT_VECTOR, FC_DT_DOUBLE, data_buf);
    fail_unless(rc == FC_SUCCESS, "abort: failed to create displ var");
    // create whole subset
    fc_createSubset(meshes[i], subsetNames[0], FC_AT_WHOLE_MESH, &subsets[i][0]);
    rc = fc_addMemberToSubset(subsets[i][0], 0);
    fail_unless(rc == FC_SUCCESS, "abort: failed to add member to subset");
    // create element subset
    fc_createSubset(meshes[i], subsetNames[1], assocs[1], &subsets[i][1]);
    rc = fc_addArrayMembersToSubset(subsets[i][1], numElemPerSubset[i],
				    elemIDsPerSubset[i]);
    fail_unless(rc == FC_SUCCESS, "abort: failed to add elems to subset");
    // create the other subsets as different assocs of element subset
    for (j = 2; j < numAssocType; j++) {
      rc = fc_copySubsetWithNewAssociation(subsets[i][1], subsetNames[j],
					   assocs[j], 0, &subsets[i][j]);
      fail_unless(rc == FC_SUCCESS, "abort: failed to copy subset");
    }
    // save members from vert subset for later testing
    fc_getSubsetMembersAsArray(subsets[i][numAssocType-1], &numMembers[i],
			       &memberIDs[i]);
  }

  // test each case
  for (i = 0; i < numCase; i++) {
    int meshID0 = meshPairs[i][0];
    int meshID1 = meshPairs[i][1];

    // test mesh proximity & WHOLE_MESH subset proximity
    // & whether displaced or not for each = 4 cases
    for (m = 0; m < 4; m++) {
      // test full usage
      if (m == 0)
	rc = fc_getMeshesProximity(meshes[meshID0], meshes[meshID1],
				   min_distances[i], &temp_numPair,
				   &temp_vertIDs1, &temp_vertIDs2);
      else if (m == 1)
	rc = fc_getDisplacedMeshesProximity(meshes[meshID0], displs[meshID0],
				   meshes[meshID1], displs[meshID1],
				   min_distances[i], &temp_numPair,
				   &temp_vertIDs1, &temp_vertIDs2);
      else if (m == 2)
	rc = fc_getSubsetsProximity(subsets[meshID0][0], subsets[meshID1][0],
				    min_distances[i], &temp_numPair,
				    &temp_vertIDs1, &temp_vertIDs2);
      else 
	rc = fc_getDisplacedSubsetsProximity(subsets[meshID0][0], displs[meshID0],
					     subsets[meshID1][0], displs[meshID1],
					     min_distances[i], &temp_numPair,
					     &temp_vertIDs1, &temp_vertIDs2);
      fail_unless(rc == FC_SUCCESS, "failed to get meshes proximity");
      fail_unless(temp_numPair == numPairs[i], "mismatch of numPair");
      if (i == 0) { // this is the "all" case
	int count = 0;
	for (j = 0; j < numVerts[0]; j++) {
	  for (k = 0; k < numVerts[1]; k++) {
	    fail_unless(temp_vertIDs1[count] == j, "mismatch of vertIDs1");
	    fail_unless(temp_vertIDs2[count] == k, "mismatch of vertIDs2");
	    count++;
	  }
	} 
      }
      else if (numPairs[i] == 0) { 
	fail_unless(temp_vertIDs1 == NULL, "mismatch of vertIDs1");
	fail_unless(temp_vertIDs2 == NULL, "mismatch of vertIDs2");
      }
      else {
	for (j = 0; j < temp_numPair; j++) {
	  fail_unless(temp_vertIDs1[j] == vertIDs1[i][j], "mismatch of vertIDs1");
	  fail_unless(temp_vertIDs2[j] == vertIDs2[i][j], "mismatch of vertIDs2");
	}
      }
      free(temp_vertIDs1);
      free(temp_vertIDs2);

      // special case -- test early return (both vertIDs = NULL)
      if (m == 0)
	rc = fc_getMeshesProximity(meshes[meshID0], meshes[meshID1],
				   min_distances[i], &temp_numPair,
				   NULL, NULL);
      else if (m == 1)
	rc = fc_getDisplacedMeshesProximity(meshes[meshID0], displs[meshID0],
				   meshes[meshID1], displs[meshID1],
				   min_distances[i], &temp_numPair,
				   NULL, NULL);
      else if (m == 2)
	rc = fc_getSubsetsProximity(subsets[meshID0][0], subsets[meshID1][0],
				    min_distances[i], &temp_numPair,
				    NULL, NULL);
      else 
	rc = fc_getDisplacedSubsetsProximity(subsets[meshID0][0], displs[meshID0],
					     subsets[meshID1][0], displs[meshID1],
					     min_distances[i], &temp_numPair,
					     NULL, NULL);
      fail_unless(rc == FC_SUCCESS, "failed to get meshes proximity");
      if (numPairs[i] > 0)
	fail_unless(temp_numPair == 1, "mismatch of proximal flag");
      else 
	fail_unless(temp_numPair == 0, "mismatch of proximal flag");
    }
   
    // test subset proximity on element derived subsets
    for (n = 1; n < numAssocType; n++) { 
      for (m = 0; m < 2; m++) { // whether displaced or not
	// test full usage
	if (m == 0)
	  rc = fc_getSubsetsProximity(subsets[meshID0][n], subsets[meshID1][n],
				     min_distances[i], &temp_numPair,
				     &temp_vertIDs1, &temp_vertIDs2);
	else 
	  rc = fc_getDisplacedSubsetsProximity(subsets[meshID0][n], displs[meshID0],
				      subsets[meshID1][n], displs[meshID1],
				      min_distances[i], &temp_numPair,
				      &temp_vertIDs1, &temp_vertIDs2);	
	fail_unless(rc == FC_SUCCESS, "failed to get meshes proximity");
	fail_unless(temp_numPair == numSubsetPairs[i], "mismatch of numPair");
	if (i == 0) { // this is the "all" case
	  int count = 0;
	  for (j = 0; j < numMembers[meshID0]; j++) {
	    for (k = 0; k < numMembers[meshID1]; k++) {
	      fail_unless(temp_vertIDs1[count] == memberIDs[meshID0][j],
			  "mismatch of vertIDs1");
	      fail_unless(temp_vertIDs2[count] == memberIDs[meshID1][k],
			  "mismatch of vertIDs2");
	      count++;
	    }
	  } 
	}
	else if (numPairs[i] == 0) { 
	  fail_unless(temp_vertIDs1 == NULL, "mismatch of vertIDs1");
	  fail_unless(temp_vertIDs2 == NULL, "mismatch of vertIDs2");
	}
	else {
	  for (j = 0; j < temp_numPair; j++) {
	    fail_unless(temp_vertIDs1[j] == vertIDs1[i][j], "mismatch of vertIDs1");
	    fail_unless(temp_vertIDs2[j] == vertIDs2[i][j], "mismatch of vertIDs2");
	  }
	}
	free(temp_vertIDs1);
	free(temp_vertIDs2);

	// special case -- test early return (both vertIDs = NULL)
	if (m == 0)
	  rc = fc_getSubsetsProximity(subsets[meshID0][n], subsets[meshID1][n],
				     min_distances[i], &temp_numPair,
				     NULL, NULL);
	else 
	  rc = fc_getDisplacedSubsetsProximity(subsets[meshID0][n], displs[meshID0],
				      subsets[meshID1][n], displs[meshID1],
				      min_distances[i], &temp_numPair,
				      NULL, NULL);	
	fail_unless(rc == FC_SUCCESS, "failed to get meshes proximity");
	if (numSubsetPairs[i] > 0)
	  fail_unless(temp_numPair == 1, "mismatch of proximal flag");
	else 
	  fail_unless(temp_numPair == 0, "mismatch of proximal flag");
      }
    }
  }

  // --- test boundary cases and errors

  // fc_getMeshesProximity -- an error to call on same mesh
  rc = fc_getMeshesProximity(meshes[0], meshes[0], min_distances[0],
			     &temp_numPair, &temp_vertIDs1, &temp_vertIDs2);
  fail_unless(rc != FC_SUCCESS, "should fail if mesh1 == mesh2");
  fail_unless(temp_numPair == -1 && temp_vertIDs1 == NULL && 
	      temp_vertIDs2 == NULL,  "fail should return nulls");

  // fc_getMeshesProximity -- test bad args
  rc = fc_getMeshesProximity(badMesh, meshes[1], min_distances[0],
			     &temp_numPair, &temp_vertIDs1, &temp_vertIDs2);
  fail_unless(rc != FC_SUCCESS, "should fail if bad mesh1");
  fail_unless(temp_numPair == -1 && temp_vertIDs1 == NULL && 
	      temp_vertIDs2 == NULL,  "fail should return nulls");
  rc = fc_getMeshesProximity(meshes[0], badMesh, min_distances[0],
			     &temp_numPair, &temp_vertIDs1, &temp_vertIDs2);
  fail_unless(rc != FC_SUCCESS, "should fail if bad mesh2");
  fail_unless(temp_numPair == -1 && temp_vertIDs1 == NULL && 
	      temp_vertIDs2 == NULL,  "fail should return nulls");
  rc = fc_getMeshesProximity(meshes[0], meshes[1], -1.,
			     &temp_numPair, &temp_vertIDs1, &temp_vertIDs2);
  fail_unless(rc != FC_SUCCESS, "should fail if negative min distance");
  fail_unless(temp_numPair == -1 && temp_vertIDs1 == NULL && 
	      temp_vertIDs2 == NULL,  "fail should return nulls");
  rc = fc_getMeshesProximity(meshes[0], meshes[1], min_distances[0],
			     NULL, &temp_vertIDs1, &temp_vertIDs2);
  fail_unless(rc != FC_SUCCESS, "should fail if no numPair");
  fail_unless(temp_vertIDs1 == NULL && 
	      temp_vertIDs2 == NULL,  "fail should return nulls");
  rc = fc_getMeshesProximity(meshes[0], meshes[1], min_distances[0],
			     &temp_numPair, NULL, &temp_vertIDs2);
  fail_unless(rc != FC_SUCCESS, "should fail if no vertIDs1");
  fail_unless(temp_numPair == -1 && 
	      temp_vertIDs2 == NULL,  "fail should return nulls");
  rc = fc_getMeshesProximity(meshes[0], meshes[1], min_distances[0],
			     &temp_numPair, &temp_vertIDs1, NULL);
  fail_unless(rc != FC_SUCCESS, "should fail if no vertIDs2");
  fail_unless(temp_numPair == -1 && temp_vertIDs1 == NULL,
	      "fail should return nulls");

  // fc_getDisplacedMeshesProximity -- an error to call on same mesh AND var,
  // but o.k. to have same mesh or same var (as long as not both)
  rc = fc_getDisplacedMeshesProximity(meshes[0], displs[0], meshes[0], 
                             displs[0], min_distances[0],
			     &temp_numPair, &temp_vertIDs1, &temp_vertIDs2);
  fail_unless(rc != FC_SUCCESS, "should fail if mesh1,displ1 == mesh2,displ2");
  fail_unless(temp_numPair == -1 && temp_vertIDs1 == NULL && 
	      temp_vertIDs2 == NULL,  "fail should return nulls");
  rc = fc_getDisplacedMeshesProximity(meshes[1], displs[1], meshes[1], 
                             displs[2], min_distances[0],
			     &temp_numPair, &temp_vertIDs1, &temp_vertIDs2);
  fail_unless(rc == FC_SUCCESS, "o.k. if mesh1 = mesh2 if displ1 != displ2");
  free(temp_vertIDs1);
  free(temp_vertIDs2);
  rc = fc_getDisplacedMeshesProximity(meshes[1], displs[1], meshes[2], 
                             displs[2], min_distances[0],
			     &temp_numPair, &temp_vertIDs1, &temp_vertIDs2);
  fail_unless(rc == FC_SUCCESS, "o.k. if displ1 = displ2 if mesh1 != mesh2");
  free(temp_vertIDs1);
  free(temp_vertIDs2);

  // fc_getDisplacedMeshesProximity -- test bad args
  rc = fc_getDisplacedMeshesProximity(badMesh, displs[0], meshes[1],
			    displs[1], min_distances[0],
			    &temp_numPair, &temp_vertIDs1, &temp_vertIDs2);
  fail_unless(rc != FC_SUCCESS, "should fail if bad mesh1");
  fail_unless(temp_numPair == -1 && temp_vertIDs1 == NULL && 
	      temp_vertIDs2 == NULL,  "fail should return nulls");
  rc = fc_getDisplacedMeshesProximity(meshes[0], badVar, meshes[1],
			    displs[1], min_distances[0],
			    &temp_numPair, &temp_vertIDs1, &temp_vertIDs2);
  fail_unless(rc != FC_SUCCESS, "should fail if bad displ1");
  fail_unless(temp_numPair == -1 && temp_vertIDs1 == NULL && 
	      temp_vertIDs2 == NULL,  "fail should return nulls");
  rc = fc_getDisplacedMeshesProximity(meshes[0], displs[1], meshes[1],
			    displs[1], min_distances[0],
			    &temp_numPair, &temp_vertIDs1, &temp_vertIDs2);
  fail_unless(rc != FC_SUCCESS, "should fail if displ1 not valid displ for mesh1");
  fail_unless(temp_numPair == -1 && temp_vertIDs1 == NULL && 
	      temp_vertIDs2 == NULL,  "fail should return nulls");
  rc = fc_getDisplacedMeshesProximity(meshes[0], displs[0], badMesh, 
			     displs[1], min_distances[0],
			     &temp_numPair, &temp_vertIDs1, &temp_vertIDs2);
  fail_unless(rc != FC_SUCCESS, "should fail if bad mesh2");
  fail_unless(temp_numPair == -1 && temp_vertIDs1 == NULL && 
	      temp_vertIDs2 == NULL,  "fail should return nulls");
  rc = fc_getDisplacedMeshesProximity(meshes[0], displs[0], meshes[1],
			    badVar, min_distances[0],
			    &temp_numPair, &temp_vertIDs1, &temp_vertIDs2);
  fail_unless(rc != FC_SUCCESS, "should fail if bad displ2");
  fail_unless(temp_numPair == -1 && temp_vertIDs1 == NULL && 
	      temp_vertIDs2 == NULL,  "fail should return nulls");
  rc = fc_getDisplacedMeshesProximity(meshes[0], displs[0], meshes[1],
			    displs[0], min_distances[0],
			    &temp_numPair, &temp_vertIDs1, &temp_vertIDs2);
  fail_unless(rc != FC_SUCCESS, "should fail if displ2 not valid displ for mesh2");
  fail_unless(temp_numPair == -1 && temp_vertIDs1 == NULL && 
	      temp_vertIDs2 == NULL,  "fail should return nulls");
  rc = fc_getDisplacedMeshesProximity(meshes[0], displs[0], meshes[1],
			     displs[1], -1.,
			     &temp_numPair, &temp_vertIDs1, &temp_vertIDs2);
  fail_unless(rc != FC_SUCCESS, "should fail if negative min distance");
  fail_unless(temp_numPair == -1 && temp_vertIDs1 == NULL && 
	      temp_vertIDs2 == NULL,  "fail should return nulls");
  rc = fc_getDisplacedMeshesProximity(meshes[0], displs[0], meshes[1],
			     displs[1], min_distances[0],
			     NULL, &temp_vertIDs1, &temp_vertIDs2);
  fail_unless(rc != FC_SUCCESS, "should fail if no numPair");
  fail_unless(temp_vertIDs1 == NULL && 
	      temp_vertIDs2 == NULL,  "fail should return nulls");
  rc = fc_getDisplacedMeshesProximity(meshes[0], displs[0], meshes[1], 
                             displs[1], min_distances[0],
			     &temp_numPair, NULL, &temp_vertIDs2);
  fail_unless(rc != FC_SUCCESS, "should fail if no vertIDs1");
  fail_unless(temp_numPair == -1 && 
	      temp_vertIDs2 == NULL,  "fail should return nulls");
  rc = fc_getDisplacedMeshesProximity(meshes[0], displs[0], meshes[1], 
			     displs[1], min_distances[0],
			     &temp_numPair, &temp_vertIDs1, NULL);
  fail_unless(rc != FC_SUCCESS, "should fail if no vertIDs2");
  fail_unless(temp_numPair == -1 && temp_vertIDs1 == NULL,
	      "fail should return nulls");

  // fc_getSubsetsProximity -- an error to call on same subset
  rc = fc_getSubsetsProximity(subsets[0][0], subsets[0][0], min_distances[0],
			     &temp_numPair, &temp_vertIDs1, &temp_vertIDs2);
  fail_unless(rc != FC_SUCCESS, "should fail if subset1 == subset2");
  fail_unless(temp_numPair == -1 && temp_vertIDs1 == NULL && 
	      temp_vertIDs2 == NULL,  "fail should return nulls");

  // fc_getSubsetProximity -- o.k. to call on diff subsets on same mesh
  rc = fc_getSubsetsProximity(subsets[0][0], subsets[0][1], 0,
			      &temp_numPair, &temp_vertIDs1, &temp_vertIDs2);
  fail_unless(rc == FC_SUCCESS, "should be o.k if subsets are on same mesh");
  free(temp_vertIDs1);
  free(temp_vertIDs2);

  // fc_getSubsetsProximity -- test bad args
  rc = fc_getSubsetsProximity(badSubset, subsets[1][0], min_distances[0],
			     &temp_numPair, &temp_vertIDs1, &temp_vertIDs2);
  fail_unless(rc != FC_SUCCESS, "should fail if bad subset1");
  fail_unless(temp_numPair == -1 && temp_vertIDs1 == NULL && 
	      temp_vertIDs2 == NULL,  "fail should return nulls");
  rc = fc_getSubsetsProximity(subsets[0][0], badSubset, min_distances[0],
			     &temp_numPair, &temp_vertIDs1, &temp_vertIDs2);
  fail_unless(rc != FC_SUCCESS, "should fail if bad subset2");
  fail_unless(temp_numPair == -1 && temp_vertIDs1 == NULL && 
	      temp_vertIDs2 == NULL,  "fail should return nulls");
  rc = fc_getSubsetsProximity(subsets[0][0], subsets[1][0], -1.,
			     &temp_numPair, &temp_vertIDs1, &temp_vertIDs2);
  fail_unless(rc != FC_SUCCESS, "should fail if negative min distance");
  fail_unless(temp_numPair == -1 && temp_vertIDs1 == NULL && 
	      temp_vertIDs2 == NULL,  "fail should return nulls");
  rc = fc_getSubsetsProximity(subsets[0][0], subsets[1][0], min_distances[0],
			     NULL, &temp_vertIDs1, &temp_vertIDs2);
  fail_unless(rc != FC_SUCCESS, "should fail if no numPair");
  fail_unless(temp_vertIDs1 == NULL && 
	      temp_vertIDs2 == NULL,  "fail should return nulls");
  rc = fc_getSubsetsProximity(subsets[0][0], subsets[1][0], min_distances[0],
			     &temp_numPair, NULL, &temp_vertIDs2);
  fail_unless(rc != FC_SUCCESS, "should fail if no vertIDs1");
  fail_unless(temp_numPair == -1 && 
	      temp_vertIDs2 == NULL,  "fail should return nulls");
  rc = fc_getSubsetsProximity(subsets[0][0], subsets[1][0], min_distances[0],
			     &temp_numPair, &temp_vertIDs1, NULL);
  fail_unless(rc != FC_SUCCESS, "should fail if no vertIDs2");
  fail_unless(temp_numPair == -1 && temp_vertIDs1 == NULL,
	      "fail should return nulls");

  // fc_getDisplacedSubsetsProximity -- an error to call on same subset AND var,
  // but o.k. to have same subset or same var (as long as not both)
  rc = fc_getDisplacedSubsetsProximity(subsets[0][0], displs[0], subsets[0][0], 
                             displs[0], min_distances[0],
			     &temp_numPair, &temp_vertIDs1, &temp_vertIDs2);
  fail_unless(rc != FC_SUCCESS, "should fail if subset1,displ1 == subset2,displ2");
  fail_unless(temp_numPair == -1 && temp_vertIDs1 == NULL && 
	      temp_vertIDs2 == NULL,  "fail should return nulls");
  rc = fc_getDisplacedSubsetsProximity(subsets[1][0], displs[1], subsets[1][0], 
                             displs[2], min_distances[0],
			     &temp_numPair, &temp_vertIDs1, &temp_vertIDs2);
  fail_unless(rc == FC_SUCCESS, "o.k. if subset1 = subset2 if displ1 != displ2");
  free(temp_vertIDs1);
  free(temp_vertIDs2);
  rc = fc_getDisplacedSubsetsProximity(subsets[1][0], displs[1], subsets[2][0], 
                             displs[2], min_distances[0],
			     &temp_numPair, &temp_vertIDs1, &temp_vertIDs2);
  fail_unless(rc == FC_SUCCESS, "o.k. if displ1 = displ2 if subset1 != subset2");
  free(temp_vertIDs1);
  free(temp_vertIDs2);

  // fc_getDisplacedSubsetProximity -- o.k. to call on diff subsets on same mesh
  rc = fc_getDisplacedSubsetsProximity(subsets[0][0], displs[0], subsets[0][1],
				       displs[0], min_distances[0],
			     &temp_numPair, &temp_vertIDs1, &temp_vertIDs2);
  fail_unless(rc == FC_SUCCESS, "o.k. if subset1 = subset2 if displ1 != displ2");
  free(temp_vertIDs1);
  free(temp_vertIDs2);
 
  // fc_getDisplacedSubsetsProximity -- test bad args
  rc = fc_getDisplacedSubsetsProximity(badSubset, displs[0], subsets[1][0],
			    displs[1], min_distances[0],
			    &temp_numPair, &temp_vertIDs1, &temp_vertIDs2);
  fail_unless(rc != FC_SUCCESS, "should fail if bad subset1");
  fail_unless(temp_numPair == -1 && temp_vertIDs1 == NULL && 
	      temp_vertIDs2 == NULL,  "fail should return nulls");
  rc = fc_getDisplacedSubsetsProximity(subsets[0][0], badVar, subsets[1][0],
			    displs[1], min_distances[0],
			    &temp_numPair, &temp_vertIDs1, &temp_vertIDs2);
  fail_unless(rc != FC_SUCCESS, "should fail if bad displ1");
  fail_unless(temp_numPair == -1 && temp_vertIDs1 == NULL && 
	      temp_vertIDs2 == NULL,  "fail should return nulls");
  rc = fc_getDisplacedSubsetsProximity(subsets[0][0], displs[1], subsets[1][0],
			    displs[1], min_distances[0],
			    &temp_numPair, &temp_vertIDs1, &temp_vertIDs2);
  fail_unless(rc != FC_SUCCESS, "should fail if displ1 not valid displ for subset1");
  fail_unless(temp_numPair == -1 && temp_vertIDs1 == NULL && 
	      temp_vertIDs2 == NULL,  "fail should return nulls");
  rc = fc_getDisplacedSubsetsProximity(subsets[0][0], displs[0], badSubset, 
			     displs[1], min_distances[0],
			     &temp_numPair, &temp_vertIDs1, &temp_vertIDs2);
  fail_unless(rc != FC_SUCCESS, "should fail if bad subset2");
  fail_unless(temp_numPair == -1 && temp_vertIDs1 == NULL && 
	      temp_vertIDs2 == NULL,  "fail should return nulls");
  rc = fc_getDisplacedSubsetsProximity(subsets[0][0], displs[0], subsets[1][0],
			    badVar, min_distances[0],
			    &temp_numPair, &temp_vertIDs1, &temp_vertIDs2);
  fail_unless(rc != FC_SUCCESS, "should fail if bad displ2");
  fail_unless(temp_numPair == -1 && temp_vertIDs1 == NULL && 
	      temp_vertIDs2 == NULL,  "fail should return nulls");
  rc = fc_getDisplacedSubsetsProximity(subsets[0][0], displs[0], subsets[1][0],
			    displs[0], min_distances[0],
			    &temp_numPair, &temp_vertIDs1, &temp_vertIDs2);
  fail_unless(rc != FC_SUCCESS, "should fail if displ2 not valid displ for subset2");
  fail_unless(temp_numPair == -1 && temp_vertIDs1 == NULL && 
	      temp_vertIDs2 == NULL,  "fail should return nulls");
  rc = fc_getDisplacedSubsetsProximity(subsets[0][0], displs[0], subsets[1][0],
			     displs[1], -1.,
			     &temp_numPair, &temp_vertIDs1, &temp_vertIDs2);
  fail_unless(rc != FC_SUCCESS, "should fail if negative min distance");
  fail_unless(temp_numPair == -1 && temp_vertIDs1 == NULL && 
	      temp_vertIDs2 == NULL,  "fail should return nulls");
  rc = fc_getDisplacedSubsetsProximity(subsets[0][0], displs[0], subsets[1][0],
			     displs[1], min_distances[0],
			     NULL, &temp_vertIDs1, &temp_vertIDs2);
  fail_unless(rc != FC_SUCCESS, "should fail if no numPair");
  fail_unless(temp_vertIDs1 == NULL && 
	      temp_vertIDs2 == NULL,  "fail should return nulls");
  rc = fc_getDisplacedSubsetsProximity(subsets[0][0], displs[0], subsets[1][0], 
                             displs[1], min_distances[0],
			     &temp_numPair, NULL, &temp_vertIDs2);
  fail_unless(rc != FC_SUCCESS, "should fail if no vertIDs1");
  fail_unless(temp_numPair == -1 && 
	      temp_vertIDs2 == NULL,  "fail should return nulls");
  rc = fc_getDisplacedSubsetsProximity(subsets[0][0], displs[0], subsets[1][0], 
			     displs[1], min_distances[0],
			     &temp_numPair, &temp_vertIDs1, NULL);
  fail_unless(rc != FC_SUCCESS, "should fail if no vertIDs2");
  fail_unless(temp_numPair == -1 && temp_vertIDs1 == NULL,
	      "fail should return nulls");
  // cleanup
  for (i = 0; i < numMesh; i++)
    free(memberIDs[i]);
  rc = fc_deleteDataset(dataset);
  fail_unless(rc == FC_SUCCESS, "failed to delete dataset at end of test");
}
END_TEST

// test getting vertices in bounding box & compare to known answer
START_TEST(verts_in_bb)
{
  FC_ReturnCode rc;
  int i, j;
  double side_length = 2.3;
  FC_Dataset dataset;
  FC_Mesh mesh, badMesh = { 999, 999 };
  FC_Variable displ;
  int numVert = 27; // numElem = 8; // 2 x 2 x 2 hex mesh
  FC_VertexBin* bin, *displBin, badBin;
  FC_Coords low =    { -0.1, -0.1, -0.1 };
  FC_Coords bottom = { 0., 0., 0. };
  FC_Coords hey =    { 0.5*side_length, 0.5*side_length, 0.5*side_length };
  FC_Coords mid1 =   { 0.9*side_length, 0.9*side_length, 0.9*side_length };
  FC_Coords center = {     side_length,     side_length,     side_length };
  FC_Coords mid2 =   { 1.1*side_length, 1.1*side_length, 1.1*side_length };
  FC_Coords top =    { 2.0*side_length, 2.0*side_length, 2.0*side_length };
  FC_Coords high =   { 2.1*side_length, 2.1*side_length, 2.1*side_length };
  FC_Coords high2 =  { 5.0*side_length, 5.0*side_length, 5.0*side_length };
  int numCase = 7;
  // Cases: 0) center point - on boundary of bb, 1) center - well inside bb,
  // 2) the lower octant - on boundary of bb, 3) lower octant - well inside bb,
  // 4) nothing (but within mesh bb) 5) nothing (outside of mesh bb) 
  // 6) the whole thing
  FC_Coords* lows[7] =  { &center, &mid1, &bottom, &low,  &hey,  &high,  &low };
  FC_Coords* highs[7] = { &center, &mid2, &center, &mid2, &mid1, &high2, &high };
  int numWithin[7] = { 1, 1, 8, 8, 0, 0, 27 };
  int withinIDs[7][27] = { { 13 },
			   { 13 },
                           { 0, 1, 3, 4, 9, 10, 12, 13 },
                           { 0, 1, 3, 4, 9, 10, 12, 13 } };
  int* ids, numID;
  FC_Coords temp_lows, temp_highs;
  
  // setup - finish filling in withinIDs
  for (j = 0; j < numVert; j++)
    withinIDs[numCase-1][j] = j;

  // setup - create test mesh
  fc_createDataset("temp.xxx", &dataset);
  rc = fc_createSimpleHexMesh(dataset, "simple hex mesh", 2, 2, 2, bottom, top,
			      &mesh);
  fail_unless(rc == FC_SUCCESS, "test aborted: failed to create simple mesh");
  rc = fc_getMeshCoordsAsVariable(mesh, &displ);
  fail_unless(rc == FC_SUCCESS, "abort: failed to get coords as var");

  // create bins
  badBin.mesh = badMesh;
  rc = fc_createMeshVertexBin(mesh, &bin);
  fail_unless(rc == FC_SUCCESS, "test aborted: failed to create vertex bin");
  rc = fc_createDisplacedMeshVertexBin(mesh, displ, &displBin);
  fail_unless(rc == FC_SUCCESS, "test aborted: failed to create displ vertex bin");
  // can't do this test - bins do not turn out the same on calico (?)
  // but just roundoff error - we still get the right answers
  // expect the bins to come out exactly the same
  //for (i = 0; i < 3; i++)
  //  for (j = 0; j < 3; j++)
  //    for (k = 0; k < 3; k++) 
  //      fail_unless(bin->vertexIDs[i][j][k][0] == displBin->vertexIDs[i][j][k][0],
  //	    "abort: bin and displBin should be similar");

  // check the different cases
  for (i = 0; i < numCase; i++) {
    rc = _fc_getVertexIDsWithinBoundingBoxUsingVertexBin(mesh, bin, *lows[i], 
                                                   *highs[i], &numID, &ids);
    fail_unless(rc == FC_SUCCESS, "failed to get vertex IDs");
    fail_unless(numID == numWithin[i], "mismatch in numID");
    for (j = 0; j < numID; j++) 
      fail_unless(withinIDs[i][j] == ids[j], "mismach in neighbor ids");
    free(ids);

    // displaced case (multiply by 2)
    for (j = 0; j < 3; j++) {
      temp_lows[j] = 2.*(*lows[i])[j];
      temp_highs[j] = 2.*(*highs[i])[j];
    }
    rc = _fc_getVertexIDsWithinBoundingBoxUsingVertexBin(mesh, displBin, 
                                   temp_lows, temp_highs, &numID, &ids);
    fail_unless(rc == FC_SUCCESS, "failed to get displ vertex IDs");
    //printf("%d: numID = %d, numID_correct = %d\n", i, numID, numWithin[i]);
    fail_unless(numID == numWithin[i], "mismatch in displ numID");
    for (j = 0; j < numID; j++) 
    fail_unless(withinIDs[i][j] == ids[j], "mismach in displ neighbor ids");
    free(ids);
  }

  // test error conditions
  rc = _fc_getVertexIDsWithinBoundingBoxUsingVertexBin(badMesh, bin, *lows[0], 
						    *highs[0], &numID, &ids);
  fail_unless(rc != FC_SUCCESS, "should fail with bad mesh");
  fail_unless(numID == -1 && ids == NULL, "fail should return nulls");
  rc = _fc_getVertexIDsWithinBoundingBoxUsingVertexBin(mesh, &badBin, *lows[0],
						     *highs[0], &numID, &ids);
  fail_unless(rc != FC_SUCCESS, "should fail with bad bin");
  fail_unless(numID == -1 && ids == NULL, "fail should return nulls");
  rc = _fc_getVertexIDsWithinBoundingBoxUsingVertexBin(mesh, NULL, *lows[0], 
						     *highs[0], &numID, &ids);
  fail_unless(rc != FC_SUCCESS, "should fail with NULL bin");
  fail_unless(numID == -1 && ids == NULL, "fail should return nulls");
  rc = _fc_getVertexIDsWithinBoundingBoxUsingVertexBin(mesh, bin, NULL, 
						     *highs[0], &numID, &ids);
  fail_unless(rc != FC_SUCCESS, "should fail with null lowpoint");
  fail_unless(numID == -1 && ids == NULL, "fail should return nulls");
  rc = _fc_getVertexIDsWithinBoundingBoxUsingVertexBin(mesh, bin, *lows[0], 
						       NULL, &numID, &ids);
  fail_unless(rc != FC_SUCCESS, "should fail with null high point");
  fail_unless(numID == -1 && ids == NULL, "fail should return nulls");
  rc = _fc_getVertexIDsWithinBoundingBoxUsingVertexBin(mesh, bin, *lows[0], 
						       *highs[0], NULL, &ids);
  fail_unless(rc != FC_SUCCESS, "should fail with null numID");
  fail_unless(ids == NULL, "fail should return nulls");
  rc = _fc_getVertexIDsWithinBoundingBoxUsingVertexBin(mesh, bin, *lows[0], 
						     *highs[0], &numID, NULL);
  fail_unless(rc != FC_SUCCESS, "should fail with null ids");
  fail_unless(numID == -1, "fail should return nulls");

  fc_freeVertexBin(bin);
  fc_freeVertexBin(displBin);
  fc_deleteDataset(dataset);
}
END_TEST

// test getting entities in bounding box & compare to brute force method
START_TEST(entities_in_bb)
{
  FC_ReturnCode rc;
  int i, j, k, m, n;
  double side_length = 2.3;
  FC_Dataset dataset;
  FC_Mesh mesh, badMesh = { 999, 999 };
  //int numVert = 27, numElem = 8; // 2 x 2 x 2 hex mesh
  FC_VertexBin *bin, badBin;
  double* coords;
  FC_Coords low =    { -0.1, -0.1, -0.1 };
  FC_Coords bottom = { 0., 0., 0. };
  FC_Coords hey =    { 0.5*side_length, 0.5*side_length, 0.5*side_length };
  FC_Coords mid1 =   { 0.9*side_length, 0.9*side_length, 0.9*side_length };
  FC_Coords center = {     side_length,     side_length,     side_length };
  FC_Coords mid2 =   { 1.1*side_length, 1.1*side_length, 1.1*side_length };
  FC_Coords top =    { 2.0*side_length, 2.0*side_length, 2.0*side_length };
  FC_Coords high =   { 2.1*side_length, 2.1*side_length, 2.1*side_length };
  int numCase = 6;
  // Cases: 0) center point - on boundary of bb, 1) center - well inside bb,
  // 2) the lower octant - on boundary of bb, 3) lower octant - well inside bb,
  // 4) nothing 5) the whole thing
  // Cases: 0) the center, 1) the lower octant, 2) the whole thing, 3) nothing
  FC_Coords* lows[6] =  { &center, &mid1, &bottom, &low,  &hey, &low };
  FC_Coords* highs[6] = { &center, &mid2, &center, &mid2, &mid1, &high };
  int* ids, numID, numVertPerEnt, *vertIDs, numID_calc, ids_calc[8*12];
  int numEntityType = 5, numEntity;
  FC_AssociationType assocs[5] = { FC_AT_VERTEX, FC_AT_EDGE, FC_AT_FACE, 
                                   FC_AT_ELEMENT, FC_AT_WHOLE_MESH };

  // setup - create test mesh
  fc_createDataset("temp.xxx", &dataset);
  rc = fc_createSimpleHexMesh(dataset, "simple hex mesh", 2, 2, 2, bottom, top,
			      &mesh);
  fail_unless(rc == FC_SUCCESS, "test aborted: failed to create simple mesh");
  fc_getMeshCoordsPtr(mesh, &coords);

  // create bins
  rc = fc_createMeshVertexBin(mesh, &bin);
  fail_unless(rc == FC_SUCCESS, "test aborted: failed to create vertex bin");
  badBin.mesh = badMesh;

  // test all of the cases for all assocs on hex mesh
  for (i = 0; i < numEntityType; i++) {
    fc_getMeshNumEntity(mesh, assocs[i], &numEntity);
    
    for (j = 0; j < numCase; j++) {

      // calc brute force answer
      numID_calc = 0;
      for (k = 0; k < numEntity; k++) {
	int isInside = 1;
	if (assocs[i] == FC_AT_VERTEX) {
	  numVertPerEnt = 1;
	  vertIDs = malloc(sizeof(int));
	  vertIDs[0] = k;
	}
	else {
	  fc_getMeshEntityChildren(mesh, assocs[i], k, FC_AT_VERTEX, 
				   &numVertPerEnt, &vertIDs);
	}
	for (m = 0; m < numVertPerEnt; m++) {
	  for (n = 0; n < 3; n++) {
	    if (coords[vertIDs[m]*3 + n] < (*lows[j])[n] ||
		coords[vertIDs[m]*3 + n] > (*highs[j])[n]) {
	      isInside = 0;
	      break;
	    }
	  }
	  if (isInside == 0)
	    break;
	}
	free(vertIDs);
	if (isInside) {
	  ids_calc[numID_calc] = k;
	  numID_calc++;
	}
      }
 
      // get binned answer
      rc = fc_getEntitiesWithinBoundingBox(mesh, bin, *lows[j], *highs[j],
					   assocs[i], &numID, &ids);
      fail_unless(rc == FC_SUCCESS, "failed to get entities");
      
      // compare
      fail_unless(numID == numID_calc, "mismatch in numID");
      if (numID == 0)
	fail_unless(ids == NULL, "mismatch of ids");
      else
	fail_unless(!memcmp(ids, ids_calc, numID*sizeof(int)), 
		    "mismatch of ids");
      free(ids);
    }

    // --- test bad args

    // should fail if bad mesh
    rc = fc_getEntitiesWithinBoundingBox(badMesh, bin, *lows[0], *highs[0], 
                                    assocs[i], &numID, &ids);
    fail_unless(rc != FC_SUCCESS, "should fail if bad mesh");
    fail_unless(numID == -1 && ids == NULL, "failure should return nulls");

    // should fail if bad bin
    rc = fc_getEntitiesWithinBoundingBox(mesh, NULL, *lows[0], *highs[0], 
					 assocs[i], &numID, &ids);
    fail_unless(rc != FC_SUCCESS, "should fail if null bin");
    fail_unless(numID == -1 && ids == NULL, "failure should return nulls");
    rc = fc_getEntitiesWithinBoundingBox(mesh, &badBin, *lows[0], *highs[0],
                                    assocs[i], &numID, &ids);
    fail_unless(rc != FC_SUCCESS, "should fail if wrong bin");
    fail_unless(numID == -1 && ids == NULL, "failure should return nulls");

    // should fail if bad bounding box
    rc = fc_getEntitiesWithinBoundingBox(mesh, bin, NULL, *highs[0], assocs[i],
                                    &numID, &ids);
    fail_unless(rc != FC_SUCCESS, "should fail if null lowpoint");
    fail_unless(numID == -1 && ids == NULL, "failure should return nulls");    
    rc = fc_getEntitiesWithinBoundingBox(mesh, bin, *lows[0], NULL, assocs[i], 
                                    &numID, &ids);
    fail_unless(rc != FC_SUCCESS, "should fail if null highpoint");
    fail_unless(numID == -1 && ids == NULL, "failure should return nulls");
 
    // should fail if bad assocs
    rc = fc_getEntitiesWithinBoundingBox(mesh, bin, *lows[0], *highs[0], -1, 
                                    &numID, &ids);
    fail_unless(rc != FC_SUCCESS, "should fail if bad assoc");
    fail_unless(numID == -1 && ids == NULL, "failure should return nulls");  
    rc = fc_getEntitiesWithinBoundingBox(mesh, bin, *lows[0], *highs[0], 
                                    FC_AT_UNKNOWN, &numID, &ids);
    fail_unless(rc != FC_SUCCESS, "should fail if assoc = FC_AT_UNKNOWN");
    fail_unless(numID == -1 && ids == NULL, "failure should return nulls");  

    // should fail if not numID or not ids
    rc = fc_getEntitiesWithinBoundingBox(mesh, bin, *lows[0], *highs[0], 
					 assocs[i], NULL, &ids);
    fail_unless(rc != FC_SUCCESS, "should faild if no numID");
    fail_unless(ids == NULL, "faild should return nulls");
    rc = fc_getEntitiesWithinBoundingBox(mesh, bin, *lows[0], *highs[0],
					 assocs[i], &numID, NULL);
    fail_unless(rc != FC_SUCCESS, "should faild if no ids");
    fail_unless(numID == -1, "faild should return nulls");
  }

  fc_freeVertexBin(bin);
  fc_deleteDataset(dataset);
}
END_TEST

// test getting vertices within radius of point & compare to known answer
START_TEST(verts_in_sphere)
{
  FC_ReturnCode rc;
  int i, j;
  double side_length = 2.3;
  FC_Dataset dataset;
  FC_Mesh mesh, badMesh = { 999, 999 };
  FC_Variable displ;
  int numVert = 27; // numElem = 8; // 2 x 2 x 2 hex mesh
  FC_VertexBin* bin, *displBin, badBin;
  FC_Coords bottom =    { 0., 0., 0. };
  FC_Coords offcenter = { 0.5*side_length, 0.5*side_length, 0.5*side_length };
  FC_Coords center =    {     side_length,     side_length,     side_length };
  FC_Coords top =       { 2.0*side_length, 2.0*side_length, 2.0*side_length };
  FC_Coords faraway =   { 5.0*side_length, 5.0*side_length, 5.0*side_length };
  int numCase = 9;
  FC_Coords* points[9] = { &offcenter, &faraway, &center, &center, &center, 
			   &center, &center, &center, &center }; 
  double radius[9] = { side_length*0.1,  // case 0: get nothing, overlaps mesh
		       side_length*0.1,  // case 1: get nothing, no overlap
		       0.,               // case 2: get center point
                       side_length*0.99, // case 3:             r < s
		       side_length*1.0,  // case 4:             r = s
                       side_length*1.2,  // case 5:         s < r < s*sqrt(2)
                       side_length*1.5,  // case 6: s*sqrt(2) < r < s*sqrt(3)
                       side_length*1.8,  // case 7: 2*sqrt(3) < r < 2*s
                       side_length*5. }; // case 8:             r >> 2*s
  // just look at results using center vertex as center
  int numWithin[9] = { 0, 0, 1, 1, 7, 7, 19, 27, 27 };
  int withinIDs[9][27] = { { },
			   { },
			   { 13 },
			   { 13 },
                           { 4, 10, 12, 13, 14, 16, 22 },
                           { 4, 10, 12, 13, 14, 16, 22 },
                           { 1, 3, 4, 5, 7, 9, 10, 11, 12, 13, 
                             14, 15, 16, 17, 19, 21, 22, 23, 25 } };
  int* ids, numID;
  FC_Coords temp_point;
  
  // setup - finish filling in withinIDs
  for (i = numCase-2; i < numCase; i++)
    for (j = 0; j < numVert; j++)
      withinIDs[i][j] = j;

  // setup - create test mesh
  fc_createDataset("temp.xxx", &dataset);
  rc = fc_createSimpleHexMesh(dataset, "simple hex mesh", 2, 2, 2, bottom, top,
			      &mesh);
  fail_unless(rc == FC_SUCCESS, "test aborted: failed to create simple mesh");
  rc = fc_getMeshCoordsAsVariable(mesh, &displ);
  fail_unless(rc == FC_SUCCESS, "abort: failed to get coords as var");

  // create bins
  badBin.mesh = badMesh;
  rc = fc_createMeshVertexBin(mesh, &bin);
  fail_unless(rc == FC_SUCCESS, "test aborted: failed to create vertex bin");
  rc = fc_createDisplacedMeshVertexBin(mesh, displ, &displBin);
  fail_unless(rc == FC_SUCCESS, "test aborted: failed to create displ vertex bin");
  // can't do this test - bins do not turn out the same on calico (?)
  // but just roundoff error - we still get the right answers
  // expect the bins to come out exactly the same
  //for (i = 0; i < 3; i++)
  //  for (j = 0; j < 3; j++)
  //    for (k = 0; k < 3; k++) 
  //      fail_unless(bin->vertexIDs[i][j][k][0] == displBin->vertexIDs[i][j][k][0],
  //		    "abort: bin and displBin should be similar");

  // check the different cases
  for (i = 0; i < numCase; i++) {
    rc = _fc_getVertexIDsWithinSphereUsingVertexBin(mesh, bin, *points[i], 
                                                   radius[i], &numID, &ids);
    fail_unless(rc == FC_SUCCESS, "failed to get vertex IDs");
    fail_unless(numID == numWithin[i], "mismatch in numID");
    for (j = 0; j < numID; j++) 
      fail_unless(withinIDs[i][j] == ids[j], "mismach in neighbor ids");
    free(ids);

    // displaced case (multiply by 2)
    for (j = 0; j < 3; j++) 
      temp_point[j] = 2.*(*points[i])[j];
    rc = _fc_getVertexIDsWithinSphereUsingVertexBin(mesh, displBin, temp_point, 
                                                   2.*radius[i], &numID, &ids);
    fail_unless(rc == FC_SUCCESS, "failed to get vertex IDs");
    fail_unless(numID == numWithin[i], "mismatch in numID");
    for (j = 0; j < numID; j++) 
      fail_unless(withinIDs[i][j] == ids[j], "mismach in neighbor ids");
    free(ids);
  }

  // test error conditions
  rc = _fc_getVertexIDsWithinSphereUsingVertexBin(badMesh, bin, center, 
						  radius[0], &numID, &ids);
  fail_unless(rc != FC_SUCCESS, "should fail with bad mesh");
  fail_unless(numID == -1 && ids == NULL, "fail should return nulls");
  rc = _fc_getVertexIDsWithinSphereUsingVertexBin(mesh, &badBin, center, 
						  radius[0], &numID, &ids);
  fail_unless(rc != FC_SUCCESS, "should fail with bad bin");
  fail_unless(numID == -1 && ids == NULL, "fail should return nulls");
  rc = _fc_getVertexIDsWithinSphereUsingVertexBin(mesh, NULL, center, 
						  radius[0], &numID, &ids);
  fail_unless(rc != FC_SUCCESS, "should fail with NULL bin");
  fail_unless(numID == -1 && ids == NULL, "fail should return nulls");
  rc = _fc_getVertexIDsWithinSphereUsingVertexBin(mesh, bin, NULL, 
						  radius[0], &numID, &ids);
  fail_unless(rc != FC_SUCCESS, "should fail with null center");
  fail_unless(numID == -1 && ids == NULL, "fail should return nulls");
  rc = _fc_getVertexIDsWithinSphereUsingVertexBin(mesh, bin, center, 
						  -1., &numID, &ids);
  fail_unless(rc != FC_SUCCESS, "should fail with negative radius");
  fail_unless(numID == -1 && ids == NULL, "fail should return nulls");
  rc = _fc_getVertexIDsWithinSphereUsingVertexBin(mesh, bin, center, 
						  radius[0], NULL, &ids);
  fail_unless(rc != FC_SUCCESS, "should fail with null numID");
  fail_unless(ids == NULL, "fail should return nulls");
  rc = _fc_getVertexIDsWithinSphereUsingVertexBin(mesh, bin, center, 
						  radius[0], &numID, NULL);
  fail_unless(rc != FC_SUCCESS, "should fail with null ids");
  fail_unless(numID == -1, "fail should return nulls");

  fc_freeVertexBin(bin);
  fc_freeVertexBin(displBin);
  fc_deleteDataset(dataset);
}
END_TEST

// test getting elements within radius of point and compare to known answer
// WSK 10/04 - probably don't need this test anymore, but I worked hard on it!
START_TEST(elems_in_sphere)
{
  FC_ReturnCode rc;
  int i, j;
  double side_length = 2.3;
  FC_Dataset dataset;
  FC_Mesh mesh;
  FC_Coords bottom = { 0., 0., 0. };
  FC_Coords center = { 1.5*side_length, 1.5*side_length, 1.5*side_length };
  FC_Coords top =    { 3.0*side_length, 3.0*side_length, 3.0*side_length };
  int numElem = 27; // numVert = 64 // 3 x 3 x 3 hex mesh
  FC_VertexBin* bin;
  int numCase = 5;
  // pick radius a little less than upper boundary, except last case
  double radius[5] = { side_length*0.8,  // case 0: r < s*sqrt(3)/2
                       side_length*1.6,  // case 1: " < r < s*sqrt(11)/2
                       side_length*2.1,  // case 2: " < r < s*sqrt(19)/2
                       side_length*2.5,  // case 3: " < r < s*sqrt(3)*3/2
                       side_length*2.6 }; // case 4: r > s*sqrt(3)*3/2
  // just look at results using center vertex as center
  int numWithin[5] = { 0, 1, 7, 19, 27 };
  int withinIDs[5][27] = { {},
                           { 13 },
                           { 4, 10, 12, 13, 14, 16, 22 },
                           { 1, 3, 4, 5, 7, 9, 10, 11, 12, 13, 
                             14, 15, 16, 17, 19, 21, 22, 23, 25 } };
  int* ids, numID;
 
  // setup - finish filling in withinIDs
  for (i = 0; i < numElem; i++)
    withinIDs[4][i] = i;
 
  // setup - create test mesh
  fc_createDataset("temp.xxx", &dataset);
  rc = fc_createSimpleHexMesh(dataset, "simple hex mesh", 3, 3, 3, bottom, top,
			      &mesh);
  fail_unless(rc == FC_SUCCESS, "test aborted: failed to create simple mesh");
  rc = fc_createMeshVertexBin(mesh, &bin);
  fail_unless(rc == FC_SUCCESS, "test aborted: failed to create vertex bin");

  // test all of the cases
  for (i = 0; i < numCase; i++) {
    rc = fc_getEntitiesWithinSphere(mesh, bin, center, 
                                    radius[i], FC_AT_ELEMENT, &numID, &ids);
    fail_unless(rc == FC_SUCCESS, "failed to get vertices");
    fail_unless(numID == numWithin[i], "mismatch in numID");
    for (j = 0; j < numID; j++) 
      fail_unless(withinIDs[i][j] == ids[j], "mismach in ids");
    free(ids);
  }

  fc_freeVertexBin(bin);
  fc_deleteDataset(dataset);
}
END_TEST

// test getting entities within radius of point and compare to brute force
// method
START_TEST(entities_in_sphere)
{
  FC_ReturnCode rc;
  int i, j, k, m;
  double side_length = 2.3;
  FC_Dataset dataset;
  FC_Mesh mesh, badMesh = { 999, 999 };
  //int numVert = 64, numElem = 27; // 3 x 3 x 3 hex mesh
  double* coords;
  FC_VertexBin* bin, badBin;
  FC_Coords bottom =    { 0., 0., 0. };
  FC_Coords offcenter = { 0.5*side_length, 0.5*side_length, 0.5*side_length };
  FC_Coords center =    { 1.5*side_length, 1.5*side_length, 1.5*side_length };
  FC_Coords top =       { 3.0*side_length, 3.0*side_length, 3.0*side_length };
  int numCase = 8;
  // pick radius a little less than upper boundary, except last case
  FC_Coords* points[8] = { &offcenter, &center, &center, &center, &center,
			   &center, &center, &center }; 
  double radius[8] = { side_length*0.1,  // case 0: get nothing
		       0.,               // case 1: get center point
		       side_length*0.8,  // case 1: r < s*sqrt(3)/2
                       side_length*1.6,  // case 2: " < r < s*sqrt(11)/2
                       side_length*2.1,  // case 3: " < r < s*sqrt(19)/2
                       side_length*2.5,  // case 4: " < r < s*sqrt(3)*3/2
                       side_length*2.6,  // case 5: r > s*sqrt(3)*3/2
                       side_length*10 }; // case 6: the whole thing
  int numEntityType = 5, numEntity;
  FC_AssociationType assocs[5] = { FC_AT_VERTEX, FC_AT_EDGE, FC_AT_FACE, 
                                   FC_AT_ELEMENT, FC_AT_WHOLE_MESH };
  int numID, *ids, numID_calc, ids_calc[27*12], numVertPerEnt, *vertIDs;
  double distance;
  double* point;
 
  // setup - create test mesh
  fc_createDataset("temp.xxx", &dataset);
  rc = fc_createSimpleHexMesh(dataset, "simple hex mesh", 3, 3, 3, bottom, top,
			      &mesh);
  fail_unless(rc == FC_SUCCESS, "test aborted: failed to create simple mesh");
  fc_getMeshCoordsPtr(mesh, &coords);

  // create bins
  rc = fc_createMeshVertexBin(mesh, &bin);
  fail_unless(rc == FC_SUCCESS, "test aborted: failed to create vertex bin");
  badBin.mesh = badMesh;

  // test all of the radius' for all assocs on hex mesh
  for (i = 0; i < numEntityType; i++) {
    fc_getMeshNumEntity(mesh, assocs[i], &numEntity);

    for (j = 0; j < numCase; j++) {

      // calc brute force answer
      numID_calc = 0;
      for (k = 0; k < numEntity; k++) {
	int isInside = 1;
	if (assocs[i] == FC_AT_VERTEX) {
	  numVertPerEnt = 1;
	  vertIDs = malloc(sizeof(int));
	  vertIDs[0] = k;
	}
	else {
	  fc_getMeshEntityChildren(mesh, assocs[i], k, FC_AT_VERTEX, 
				   &numVertPerEnt, &vertIDs);
	}
	for (m = 0; m < numVertPerEnt; m++) {
	  point = &coords[vertIDs[m]*3];
	  rc = fc_calcEuclideanDistance(*points[j], point, 3, &distance);
	  fail_unless(rc == FC_SUCCESS, "abort: failed to calc distance");
	  if (distance > radius[j]) {
	    isInside = 0;
	    break;
	  }
	}
	free(vertIDs);
	if (isInside) {
	  ids_calc[numID_calc] = k;
	  numID_calc++;
	}
      }

      // get binned answer
      rc = fc_getEntitiesWithinSphere(mesh, bin, *points[j], radius[j], assocs[i],
                                      &numID, &ids);
      fail_unless(rc == FC_SUCCESS, "failed to get entities");
      
      // compare
      fail_unless(numID == numID_calc, "mismatch in numID");
      if (numID == 0)
        fail_unless(ids == NULL, "mismatch of ids");
      else
        fail_unless(!memcmp(ids, ids_calc, numID*sizeof(int)), 
                    "mismatch of ids");
      free(ids);
    }
    
    // --- test bad args

    // should fail if bad mesh
    rc = fc_getEntitiesWithinSphere(badMesh, bin, center, radius[0], 
                                    assocs[i], &numID, &ids);
    fail_unless(rc != FC_SUCCESS, "should fail if bad mesh");
    fail_unless(numID == -1 && ids == NULL, "failure should return nulls");

    // should fail if bad bin
    rc = fc_getEntitiesWithinSphere(mesh, NULL, center, radius[0], assocs[i], 
                                    &numID, &ids);
    fail_unless(rc != FC_SUCCESS, "should fail if null bin");
    fail_unless(numID == -1 && ids == NULL, "failure should return nulls");
    rc = fc_getEntitiesWithinSphere(mesh, &badBin, center, radius[0],
                                    assocs[i], &numID, &ids);
    fail_unless(rc != FC_SUCCESS, "should fail if wrong bin");
    fail_unless(numID == -1 && ids == NULL, "failure should return nulls");

    // should fail if bad center or radius
    rc = fc_getEntitiesWithinSphere(mesh, bin, NULL, radius[0], assocs[i], 
                                    &numID, &ids);
    fail_unless(rc != FC_SUCCESS, "should fail if null center");
    fail_unless(numID == -1 && ids == NULL, "failure should return nulls");
    rc = fc_getEntitiesWithinSphere(mesh, bin, center, -1., assocs[i], 
                                    &numID, &ids);
    fail_unless(rc != FC_SUCCESS, "should fail if radius is negative");
    fail_unless(numID == -1 && ids == NULL, "failure should return nulls");

    // should fail if bad assocs
    rc = fc_getEntitiesWithinSphere(mesh, bin, center, radius[0], -1, 
                                    &numID, &ids);
    fail_unless(rc != FC_SUCCESS, "should fail if bad assoc");
    fail_unless(numID == -1 && ids == NULL, "failure should return nulls");  
    rc = fc_getEntitiesWithinSphere(mesh, bin, center, radius[0], 
                                    FC_AT_UNKNOWN, &numID, &ids);
    fail_unless(rc != FC_SUCCESS, "should fail if assoc = FC_AT_UNKNOWN");
    fail_unless(numID == -1 && ids == NULL, "failure should return nulls");  

    // should fail if not numID or not ids
    rc = fc_getEntitiesWithinSphere(mesh, bin, center, radius[0], assocs[i], 
                                    NULL, &ids);
    fail_unless(rc != FC_SUCCESS, "should faild if no numID");
    fail_unless(ids == NULL, "faild should return nulls");
    rc = fc_getEntitiesWithinSphere(mesh, bin, center, radius[0], assocs[i],
                                    &numID, NULL);
    fail_unless(rc != FC_SUCCESS, "should faild if no ids");
    fail_unless(numID == -1, "faild should return nulls");
  }

  fc_freeVertexBin(bin);
  fc_deleteDataset(dataset);
}
END_TEST

// test smoothing of vertex associated variable
// Given a 2 x 2 x 2 hex mesh (perfect cubes), set the variable on the
// center vertex to a large value, and 0 everywhere else. Then smoothing
// results should be easy to predict by choosing radius relative to
// the cube side lengths. 
START_TEST(vertex_smooth)
{
  FC_ReturnCode rc;
  int i, j;
  char var_name[20] = "arbitrary var";
  double side_length = 2.3;
  double big_value = 10.0;
  FC_Dataset dataset;
  FC_Mesh mesh;
  FC_Variable var, *smooth_vars;
  int numVert = 27; // numElem = 8; // 2 x 2 x 2 hex mesh
  FC_Coords bottom = { 0., 0., 0. };
  FC_Coords top =    { 2.0*side_length, 2.0*side_length, 2.0*side_length };
  double data[27];
  double* temp_data;
  FC_VertexBin* bin;
  int numCase = 3;
  double radius[3] = { side_length*1.2,   // case 0:         s < r < s*sqrt(2)
                       side_length*1.5,   // case 1: s*sqrt(2) < r < s*sqrt(3)
                       side_length*1.8 }; // case 2: 2*sqrt(3) < r < 2*s
  // each value is the big_value divided by the number of vertices that
  // fall within the radius.
  //                           case 0        case 1        case 2
  double center_values[3] = { big_value/7, big_value/19, big_value/27 };
  double face_values[3] = { big_value/6, big_value/14, big_value/18 };
  double edge_values[3] = { 0., big_value/10, big_value/12 };
  double corner_values[3] = { 0., 0., big_value/8 };
  //int numCenter = 1;
  int numFace = 6, numEdge = 12, numCorner = 8;
  int center_ids[1] = { 13 };
  int face_ids[6] = { 4, 10, 12, 14, 16, 22 };
  int edge_ids[12] = { 1, 3, 5, 7, 9, 11, 15, 17, 19, 21, 23, 25 };
  int corner_ids[8] = { 0, 2, 6, 8, 18, 20, 24, 26 };

  // setup - create data
  for (i = 0; i < numVert; i++)
    data[i] = 0.;
  data[center_ids[0]] = big_value;

  // setup - create test mesh & var
  fc_createDataset("temp.xxx", &dataset);
  rc = fc_createSimpleHexMesh(dataset, "simple hex mesh", 2, 2, 2, bottom, top,
			      &mesh);
  fail_unless(rc == FC_SUCCESS, "test aborted: failed to create simple mesh");
  fc_createVariable(mesh, var_name, &var);
  rc = fc_setVariableData(var, numVert, 1, FC_AT_VERTEX, FC_MT_SCALAR, 
                   FC_DT_DOUBLE, (void*)data);
  fail_unless(rc == FC_SUCCESS, "test aborted: failed to create test data");
  rc = fc_createMeshVertexBin(mesh, &bin);
  fail_unless(rc == FC_SUCCESS, "test aborted: failed to create vertex bin");

  // smoothing with really small radius should give back original data
  // case r < side_length
  rc = fc_geomSmoothVariable(1, &mesh, &bin, &var, side_length/10., 0,
                             &smooth_vars);
  fail_unless(rc == FC_SUCCESS, "failed to smooth");
  rc = fc_getVariableDataPtr(smooth_vars[0], (void**)&temp_data);
  fail_unless(rc == FC_SUCCESS, "failed to get smoothed var data");
  fail_unless(!memcmp(data, temp_data, sizeof(double)*numVert),
              "did not get the same data");
  fc_deleteVariable(smooth_vars[0]);
  free(smooth_vars);

  // smoothing with a really large radius should give uniform results
  // case r >> 2*side_length
  rc = fc_geomSmoothVariable(1, &mesh, &bin, &var, 100.*side_length, 0,
                             &smooth_vars);
  fail_unless(rc == FC_SUCCESS, "failed to smooth");
  rc = fc_getVariableDataPtr(smooth_vars[0], (void**)&temp_data);
  fail_unless(rc == FC_SUCCESS, "failed to get smoothed var data");
  for (i = 0; i < numVert; i++) 
    fail_unless(FC_DBL_EQUIV(temp_data[i], big_value/numVert), 
                "did not get expected data");
  fc_deleteVariable(smooth_vars[0]);
  free(smooth_vars);

  // do the other 3 cases
  for (i = 0; i < numCase; i++) {
    rc = fc_geomSmoothVariable(1, &mesh, &bin, &var, radius[i], 0,
                               &smooth_vars);
    fail_unless(rc == FC_SUCCESS, "failed to smooth");
    rc = fc_getVariableDataPtr(smooth_vars[0], (void**)&temp_data);
    fail_unless(rc == FC_SUCCESS, "failed to get smoothed var data");
    fail_unless(temp_data[center_ids[0]] == center_values[i],
                "mismatch of center values");
    for (j = 0; j < numFace; j++)
      fail_unless(temp_data[face_ids[j]] == face_values[i],
                  "mismatch of face values");
    for (j = 0; j < numEdge; j++)
      fail_unless(temp_data[edge_ids[j]] == edge_values[i],
                  "mismatch of edge values");
    for (j = 0; j < numCorner; j++)
      fail_unless(temp_data[corner_ids[j]] == corner_values[i],
                  "mismatch of corner values");
    fc_deleteVariable(smooth_vars[0]);
    free(smooth_vars);
  }

  fc_freeVertexBin(bin);
  fc_deleteDataset(dataset);
}
END_TEST

// test smoothing of element associated variable
// Given a 3 x 3 x 3 hex mesh (perfect cubes), set the variable on the
// center element to a large value, and 0 everywhere else. Then smoothing
// results should be easy to predict by choosing radius relative to
// the cube side lengths. 
START_TEST(element_smooth)
{
  FC_ReturnCode rc;
  int i, j;
  char var_name[20] = "arbitrary var";
  double side_length = 2.3;
  double big_value = 10.0;
  FC_Dataset dataset;
  FC_Mesh mesh;
  FC_Variable var, *smooth_vars;
  int numElem = 27;// numVert = 64, // 3 x 3 x 3 hex mesh
  FC_Coords bottom =  { 0., 0., 0. };
  FC_Coords top =     { 3.0*side_length, 3.0*side_length, 3.0*side_length };
  double data[27];
  double* temp_data;
  FC_VertexBin* bin;
  int numCase = 3;
  double radius[3] = { //side_length*0.8,  // case 0: r < s*sqrt(3)/2
                       //side_length*1.6,  // case 1: " < r < s*sqrt(11)/2
                       side_length*2.1,  // case 2: " < r < s*sqrt(19)/2
                       side_length*2.5,  // case 3: " < r < s*sqrt(3)*3/2
                       side_length*2.6 }; // case 4: r > s*sqrt(3)*3/2
  // each value is the big_value divided by the number of vertices that
  // fall within the radius.
  //                           case 0        case 1        case 2
  // note in last case, the distance to corners is large enough to
  // get 2nd order face neighbors.
  double center_values[3] = { big_value/7, big_value/19, big_value/27 };
  double face_values[3] = { big_value/6, big_value/14, big_value/19 };
  double edge_values[3] = { 0., big_value/10, big_value/14 };
  double corner_values[3] = { 0., 0., big_value/11 };
  //int numCenter = 1;
  int numFace = 6, numEdge = 12, numCorner = 8;
  int center_ids[1] = { 13 };
  int face_ids[6] = { 4, 10, 12, 14, 16, 22 };
  int edge_ids[12] = { 1, 3, 5, 7, 9, 11, 15, 17, 19, 21, 23, 25 };
  int corner_ids[8] = { 0, 2, 6, 8, 18, 20, 24, 26 };

  // setup - create data
  for (i = 0; i < numElem; i++)
    data[i] = 0.;
  data[center_ids[0]] = big_value;

  // setup - create test mesh & var
  fc_createDataset("temp.xxx", &dataset);
  rc = fc_createSimpleHexMesh(dataset, "simple hex mesh", 3, 3, 3, bottom, top,
			      &mesh);
  fail_unless(rc == FC_SUCCESS, "test aborted: failed to create simple mesh");
  fc_createVariable(mesh, var_name, &var);
  rc = fc_setVariableData(var, numElem, 1, FC_AT_ELEMENT, FC_MT_SCALAR, 
                   FC_DT_DOUBLE, (void*)data);
  fail_unless(rc == FC_SUCCESS, "test aborted: failed to create test data");
  rc = fc_createMeshVertexBin(mesh, &bin);
  fail_unless(rc == FC_SUCCESS, "test aborted: failed to create vertex bin");
  
  // smoothing with small radius should give back original data
  // case r < side_length*sqrt(11)/2
  rc = fc_geomSmoothVariable(1, &mesh, &bin, &var, 1.6*side_length, 0,
                             &smooth_vars);
  fail_unless(rc == FC_SUCCESS, "failed to smooth");
  rc = fc_getVariableDataPtr(smooth_vars[0], (void**)&temp_data);
  fail_unless(rc == FC_SUCCESS, "failed to get smoothed var data");
  fail_unless(!memcmp(data, temp_data, sizeof(double)*numElem),
              "did not get the same data");
  fc_deleteVariable(smooth_vars[0]);
  free(smooth_vars);
 
  // smoothing with a really large radius should give uniform results
  // case r >> 2*side_length
  rc = fc_geomSmoothVariable(1, &mesh, &bin, &var, 100.*side_length, 0,
                             &smooth_vars);
  fail_unless(rc == FC_SUCCESS, "failed to smooth");
  rc = fc_getVariableDataPtr(smooth_vars[0], (void**)&temp_data);
  fail_unless(rc == FC_SUCCESS, "failed to get smoothed var data");
  for (i = 0; i < numElem; i++) 
    fail_unless(FC_DBL_EQUIV(temp_data[i], big_value/numElem), 
                "did not get expected data");
  fc_deleteVariable(smooth_vars[0]);
  free(smooth_vars);

  // do the other 3 cases
  for (i = 0; i < numCase; i++) {
    rc = fc_geomSmoothVariable(1, &mesh, &bin, &var, radius[i], 0,
                               &smooth_vars);
    fail_unless(rc == FC_SUCCESS, "failed to smooth");
    rc = fc_getVariableDataPtr(smooth_vars[0], (void**)&temp_data);
    fail_unless(rc == FC_SUCCESS, "failed to get smoothed var data");
    fail_unless(temp_data[center_ids[0]] == center_values[i],
                "mismatch of center values");
    for (j = 0; j < numFace; j++)
      fail_unless(temp_data[face_ids[j]] == face_values[i],
                  "mismatch of face values");
    for (j = 0; j < numEdge; j++)
      fail_unless(temp_data[edge_ids[j]] == edge_values[i],
                  "mismatch of edge values");
    for (j = 0; j < numCorner; j++)
      fail_unless(temp_data[corner_ids[j]] == corner_values[i],
                  "mismatch of corner values");
    fc_deleteVariable(smooth_vars[0]);
    free(smooth_vars);
  }
  fc_freeVertexBin(bin);
  fc_deleteDataset(dataset);
}
END_TEST
  
// test smoothing of all types of variables
// (WSK 10/18/04 - this test could replace the previous two, but might
// as well keep them as they check in a different way)
// Given two 3 x 3 x 3 hex meshes (perfect cubes), set the variable on the
// first entity of each mesh to a large value, and 0 everywhere else.
START_TEST(entity_smooth)
{
  FC_ReturnCode rc;
  int i, j, k, m, n, p;
  double side_length = 2.3;
  double big_values[2] = { 7.0, 11.0 };
  char var_name[20] = "arbitrary var";
  //char smooth_name[30] = "arbitrary var_smoothed";
  FC_Dataset dataset;
  int numMesh = 2;
  FC_Mesh meshes[2], badMeshes[2], badMesh = { 999, 999 };
  FC_Variable vars[2], displs[2], badVars[2], badVar = { 999, 999 };
  FC_Variable *smooth_vars, *displ_smooth_vars, *ref_smooth_vars;
  int numElem = 27;// numVert = 64, // 3 x 3 x 3 hex mesh
  FC_Coords bottom1 =  { 0., 0., 0. };
  FC_Coords top1 =     { 3.0*side_length, 3.0*side_length, 3.0*side_length };
  FC_Coords bottom2, top2;
  int numEntity;
  double* coords[2];
  double datas[2][27*12], data_calc[27*12];
  double* temp_data, *temp_displ_data, *ref_data;
  FC_VertexBin* bins[2], *displBins[2], *badBins[2];
  int numRadius = 3;
  double radius[3] = { side_length*0.1,  // smaller than edge - no smoothing
                       side_length*2.1,  // smaller than mesh - some smoothing
                       side_length*100 }; // larger than mesh - complete smooth
  int numEntityType = 5;
  FC_AssociationType assocs[5] = { FC_AT_VERTEX, FC_AT_EDGE, FC_AT_FACE, 
                                   FC_AT_ELEMENT, FC_AT_WHOLE_MESH };
  int numDataType = 3;
  FC_DataType dataTypes[3] = { FC_DT_CHAR, FC_DT_INT, FC_DT_FLOAT };
  int numComp = 3;
  FC_Coords centers[2][5][27*12]; // index by mesh, entity type, entity ID
  int numID, *ids, id;

  // setup - create test meshes, 2nd one is shifted in x direction
  bottom2[0] = bottom1[0] + side_length*4.0; 
  top2[0] = top1[0] + side_length*4.0;
  for (i = 1; i < 3; i++) {
    bottom2[i] = bottom1[i];
    top2[i] = top1[i];
  }
  fc_createDataset("temp.xxx", &dataset);
  rc = fc_createSimpleHexMesh(dataset, "simple hex mesh1", 3, 3, 3, bottom1, 
			      top1, &meshes[0]);
  fail_unless(rc == FC_SUCCESS, "test aborted: failed to create simple mesh");
  rc = fc_createSimpleHexMesh(dataset, "simple hex mesh2", 3, 3, 3, bottom2, 
			      top2, &meshes[1]);
  fail_unless(rc == FC_SUCCESS, "test aborted: failed to create simple mesh");
  fc_getMeshCoordsPtr(meshes[0], &coords[0]);
  fc_getMeshCoordsPtr(meshes[1], &coords[1]);

  // create displacement variables
  rc = fc_getMeshCoordsAsVariable(meshes[0], &displs[0]);
  fail_unless(rc == FC_SUCCESS, "abort: failed to get coords as var");
  rc = fc_getMeshCoordsAsVariable(meshes[1], &displs[1]);
  fail_unless(rc == FC_SUCCESS, "abort: failed to get coords as var");

  // create bins
  for (i = 0; i < numMesh; i++) {
    rc = fc_createMeshVertexBin(meshes[i], &bins[i]);
    fail_unless(rc == FC_SUCCESS, "test aborted: failed to create vertex bin");
    rc = fc_createDisplacedMeshVertexBin(meshes[i], displs[i], &displBins[i]);
    fail_unless(rc == FC_SUCCESS, "test aborted: failed to create vertex bin");
  }

  // setup data (same for both meshes, all entity types)
  datas[0][0] = big_values[0];
  datas[1][0] = big_values[1];
  for (i = 1; i < 27*12; i++) {
    datas[0][i] = 0.;
    datas[1][i] = 0.;
  }

  // Calculate centers of all entities
  for (i = 0; i < numMesh; i++) {
    FC_ElementType elemType;
    int numVertPerElem, maxNumVertPerFace, *numVertPerFace;
    int *edgeConns, *faceConns, *elemConns;
    fc_getMeshElementType(meshes[i], &elemType);
    numVertPerElem = fc_getElementTypeNumVertex(elemType);
    fc_getMeshEdgeConnsPtr(meshes[i], &edgeConns);
    fc_getMeshFaceConnsPtr(meshes[i], &numVertPerFace, &maxNumVertPerFace,
			   &faceConns);
    fc_getMeshElementConnsPtr(meshes[i], &elemConns);

    for (j = 0; j < numEntityType; j++) {
      fc_getMeshNumEntity(meshes[i], assocs[j], &numEntity);
     
      for (k = 0; k < numEntity; k++) {
        // get ids
        if (assocs[j] == FC_AT_VERTEX) {
          numID = 1;
	  id = k;
          ids = &id;
        }
        else if (assocs[j] == FC_AT_EDGE) {
          numID = 2;
          ids = &edgeConns[k*2];
        }
        else if (assocs[j] == FC_AT_FACE) {
          numID = numVertPerFace[k];
          ids = &faceConns[k*maxNumVertPerFace];
        }
        else if (assocs[j] == FC_AT_ELEMENT) {
          numID = numVertPerElem;
	  ids = &elemConns[k*numVertPerElem];
        }
        else if (assocs[j] == FC_AT_WHOLE_MESH) {
          for (m = 0; m < 3; m++) {
            centers[0][j][k][m] = 1.5*side_length;
            centers[1][j][k][m] = 1.5*side_length;
          }
          centers[1][j][k][0] += 4.0;
          continue;
        } 

        // compute centers
        for (m = 0; m < 3; m++)
          centers[i][j][k][m] = 0;
        for (n = 0; n < numID; n++) 
          for (m = 0; m < 3; m++) 
            //printf("%d %d %d %d: %f\n", i, j, k, m, coords[i][3*ids[n]+m]);
            centers[i][j][k][m] += coords[i][3*ids[n] + m];
        for (m = 0; m < 3; m++)
          centers[i][j][k][m] = centers[i][j][k][m] / numID;
      }
    }
  }

  // test all of the radius' for all assocs on hex mesh
  for (i = 0; i < numEntityType; i++) {

    // create the vars
    for (j = 0; j < numMesh; j++) {
      fc_getMeshNumEntity(meshes[j], assocs[i], &numEntity);
      fc_createVariable(meshes[j], var_name, &vars[j]);
      rc = fc_setVariableData(vars[j], numEntity, 1, assocs[i], 
                          FC_MT_SCALAR, FC_DT_DOUBLE, (void*)datas[j]);
      fail_unless(rc == FC_SUCCESS, 
                  "test aborted: failed to create test var");
    }

    // test the smoothing
    
    for (j = 0; j < numRadius; j++) {
      //printf("i = %d, j = %d\n", i, j);

      // get result & compare for on all-mesh smoothing
      rc = fc_geomSmoothVariable(numMesh, meshes, bins, vars, radius[j], 0,
                                 &smooth_vars);
      fail_unless(rc == FC_SUCCESS, "failed to smooth");
      rc = fc_displacedGeomSmoothVariable(numMesh, meshes, displs,
					  displBins, vars, 2.*radius[j], 0,
					  &displ_smooth_vars);
      fail_unless(rc == FC_SUCCESS, "failed to displ smooth");

      // look at each new var
      for (k = 0; k < numMesh; k++) {
        rc = fc_getVariableDataPtr(smooth_vars[k], (void**)&temp_data);
        fail_unless(rc == FC_SUCCESS, "failed to get smoothed var data");
        rc = fc_getVariableDataPtr(displ_smooth_vars[k], (void**)&temp_displ_data);
        fail_unless(rc == FC_SUCCESS, "failed to get displ smoothed var data");
       
        //  calc answer 
        for (m = 0; m < numEntity; m++) {
          int isInList = 0;
          int numTotal = 0;
          data_calc[m] = 0.;
          for (n = 0; n < numMesh; n++) {
            fc_getEntitiesWithinSphere(meshes[n], bins[n], centers[k][i][m], 
                                       radius[j], assocs[i], &numID, &ids);
            numTotal += numID;
            for (p = 0; p < numID; p++) {
              data_calc[m] += datas[n][ids[p]];
              if (ids[p] == m)
                isInList = 1;
            }
            if (n == k && !isInList) {
              data_calc[m] += datas[n][m];
              numTotal++;
            }
            free(ids);
          }
          data_calc[m] = data_calc[m] / numTotal;
          //fprintf(stderr, "Expected %f but got %f\n", data_calc[m], temp_data[m]);          
        }
        
        // compare to calc answers
        fail_unless(!memcmp(data_calc, temp_data, numEntity*sizeof(double)),
                    "mismatch of data values");
        fail_unless(!memcmp(data_calc, temp_displ_data, numEntity*sizeof(double)),
                    "mismatch of displ_data values");

        // check boundary cases
        if (j == 0) { // no smoothing
          fail_unless(temp_data[0] == big_values[k], 
                      "mismatch of expected value");
          fail_unless(temp_displ_data[0] == big_values[k], 
                      "mismatch of displ expected value");
          for (m = 1; m < numEntity; m++) {
	    fail_unless(temp_data[m] == 0., "mismatch of expected values");
            fail_unless(temp_displ_data[m] == 0., "mismatch of displ expected values");
	  }
        }
        else if (j == 2) { // completely smooth
          double ave = (big_values[0] + big_values[1]) / (2*numEntity);
          for (m = 0; m < numEntity; m++) {
            fail_unless(temp_data[m] == ave, "mismatch of expected values");
            fail_unless(temp_displ_data[m] == ave, "mismatch of displ expected values");
	  }
        }
          
        // cleanup
        fc_deleteVariable(smooth_vars[k]);
      }
      
      // cleanup
      free(smooth_vars);
      free(displ_smooth_vars);

      // get result & compare for on-mesh only smoothing
      rc = fc_geomSmoothVariable(numMesh, meshes, bins, vars, radius[j], 1,
                                 &smooth_vars);
      fail_unless(rc == FC_SUCCESS, "failed to smooth");
      rc = fc_displacedGeomSmoothVariable(numMesh, meshes, displs,
					  displBins, vars, 2.*radius[j], 1,
					  &displ_smooth_vars);
      fail_unless(rc == FC_SUCCESS, "failed to displ smooth");

      // look at each new var
      for (k = 0; k < numMesh; k++) {
        rc = fc_getVariableDataPtr(smooth_vars[k], (void**)&temp_data);
        fail_unless(rc == FC_SUCCESS, "failed to get smoothed var data");
        rc = fc_getVariableDataPtr(displ_smooth_vars[k], (void**)&temp_displ_data);
        fail_unless(rc == FC_SUCCESS, "failed to get displ smoothed var data");
       
        //  calc answer 
        for (m = 0; m < numEntity; m++) {
          int isInList = 0;
          int numTotal = 0;
          data_calc[m] = 0.;
          for (n = k; n < k+1; n++) { // Just look on the current mesh
            fc_getEntitiesWithinSphere(meshes[n], bins[n], centers[k][i][m], 
                                       radius[j], assocs[i], &numID, &ids);
            numTotal += numID;
            for (p = 0; p < numID; p++) {
              data_calc[m] += datas[n][ids[p]];
              if (ids[p] == m)
                isInList = 1;
            }
            if (n == k && !isInList) {
              data_calc[m] += datas[n][m];
              numTotal++;
            }
            free(ids);
          }
          data_calc[m] = data_calc[m] / numTotal;
          //fprintf(stderr, "Expected %f but got %f\n", data_calc[m], temp_data[m]);          
        }
        
        // compare to calc answers
        fail_unless(!memcmp(data_calc, temp_data, numEntity*sizeof(double)),
                    "mismatch of data values");
        fail_unless(!memcmp(data_calc, temp_displ_data, numEntity*sizeof(double)),
                    "mismatch of displ_data values");

        // check boundary cases
        if (j == 0) { // no smoothing
          fail_unless(temp_data[0] == big_values[k], 
                      "mismatch of expected value");
          fail_unless(temp_displ_data[0] == big_values[k], 
                      "mismatch of displ expected value");
          for (m = 1; m < numEntity; m++) {
	    fail_unless(temp_data[m] == 0., "mismatch of expected values");
            fail_unless(temp_displ_data[m] == 0., "mismatch of displ expected values");
	  }
        }
        else if (j == 2) { // completely smooth
          double ave = big_values[k] / numEntity;
          for (m = 0; m < numEntity; m++) {
            fail_unless(temp_data[m] == ave, "mismatch of expected values");
            fail_unless(temp_displ_data[m] == ave, "mismatch of displ expected values");
	  }
        }
          
        // cleanup
        fc_deleteVariable(smooth_vars[k]);
      }
      
      // cleanup
      free(smooth_vars);
      free(displ_smooth_vars);

      // --- test bad args

      // numMesh
      rc = fc_geomSmoothVariable(0, meshes, bins, vars, radius[j], 0,
                                 &smooth_vars);
      fail_unless(rc != FC_SUCCESS, "should fail if numMesh < 1");
      fail_unless(smooth_vars == NULL, "fail should return null");
      rc = fc_displacedGeomSmoothVariable(0, meshes, displs, displBins,
					  vars, radius[j], 0, &smooth_vars);
      fail_unless(rc != FC_SUCCESS, "should fail if numMesh < 1");
      fail_unless(smooth_vars == NULL, "fail should return null");

      // bad meshes
      rc = fc_geomSmoothVariable(numMesh, NULL, bins, vars, radius[j], 0,
                                 &smooth_vars);
      fail_unless(rc != FC_SUCCESS, "should fail if no meshes");
      fail_unless(smooth_vars == NULL, "fail should return null");
      rc = fc_displacedGeomSmoothVariable(numMesh, NULL, displs, displBins,
					  vars, radius[j], 0, &smooth_vars);
      fail_unless(rc != FC_SUCCESS, "should fail if no meshes");
      fail_unless(smooth_vars == NULL, "fail should return null");
      badMeshes[0] = meshes[0];
      badMeshes[1] = badMesh;
      rc = fc_geomSmoothVariable(numMesh, badMeshes, bins, vars, radius[j], 0,
                                 &smooth_vars);
      fail_unless(rc != FC_SUCCESS, "should fail if a mesh is bad");
      fail_unless(smooth_vars == NULL, "fail should return null");
      rc = fc_displacedGeomSmoothVariable(numMesh, badMeshes, displs,
					  displBins, vars, radius[j], 0,
					  &smooth_vars);
      fail_unless(rc != FC_SUCCESS, "should fail if a mesh is bad");
      fail_unless(smooth_vars == NULL, "fail should return null");

      // bad displs (displaced version only)
      rc = fc_displacedGeomSmoothVariable(numMesh, meshes, NULL,
					  displBins, vars, radius[j], 0,
					  &smooth_vars);
      fail_unless(rc != FC_SUCCESS, "should fail if no displs");
      fail_unless(smooth_vars == NULL, "fail should return null");
      badVars[0] = displs[0];
      badVars[1] = badVar;
      rc = fc_displacedGeomSmoothVariable(numMesh, meshes, badVars,
					  displBins, vars, radius[j], 0,
					  &smooth_vars);
      fail_unless(rc != FC_SUCCESS, "should fail if a displ is bad");
      fail_unless(smooth_vars == NULL, "fail should return null");
      rc = fc_displacedGeomSmoothVariable(numMesh, meshes, vars,
					  displBins, vars, radius[j], 0,
					  &smooth_vars);
      fail_unless(rc != FC_SUCCESS, "should fail if a displs not valid displ");
      fail_unless(smooth_vars == NULL, "fail should return null");
  
      // bad bins 
      rc = fc_geomSmoothVariable(numMesh, meshes, NULL, vars, radius[j], 0,
                                 &smooth_vars);
      fail_unless(rc != FC_SUCCESS, "should fail if no bins");
      fail_unless(smooth_vars == NULL, "fail should return null");
      rc = fc_displacedGeomSmoothVariable(numMesh, meshes, displs, NULL, 
					  vars, radius[j], 0, &smooth_vars);
      fail_unless(rc != FC_SUCCESS, "should fail if no bins");
      fail_unless(smooth_vars == NULL, "fail should return null");
      badBins[0] = bins[1];
      badBins[1] = bins[0];
      rc = fc_geomSmoothVariable(numMesh, meshes, badBins, vars, radius[j], 0,
                                 &smooth_vars);
      fail_unless(rc != FC_SUCCESS, "should fail if bins for wrong meshes");
      fail_unless(smooth_vars == NULL, "fail should return null");
      badBins[0] = displBins[1];
      badBins[1] = displBins[0];
      rc = fc_displacedGeomSmoothVariable(numMesh, meshes, displs, badBins,
					  vars, radius[j], 0, &smooth_vars);
      fail_unless(rc != FC_SUCCESS, "should fail if bins for wrong meshes");
      fail_unless(smooth_vars == NULL, "fail should return null");

      // bad vars 
      rc = fc_geomSmoothVariable(numMesh, meshes, bins, NULL, radius[j], 0,
                                 &smooth_vars);
      fail_unless(rc != FC_SUCCESS, "should fail if no vars");
      fail_unless(smooth_vars == NULL, "fail should return null");
      rc = fc_displacedGeomSmoothVariable(numMesh, meshes, displs, displBins, NULL,
					  radius[j], 0, &smooth_vars);
      fail_unless(rc != FC_SUCCESS, "should fail if no vars");
      fail_unless(smooth_vars == NULL, "fail should return null");
      badVars[0] = vars[1];
      badVars[1] = vars[0];
      rc = fc_geomSmoothVariable(numMesh, meshes, bins, badVars, radius[j], 0,
                                 &smooth_vars);
      fail_unless(rc != FC_SUCCESS, "should fail if vars for wrong meshes");
      fail_unless(smooth_vars == NULL, "fail should return null");
      rc = fc_displacedGeomSmoothVariable(numMesh, meshes, displs, displBins,
					  badVars, radius[j], 0, &smooth_vars);
      fail_unless(rc != FC_SUCCESS, "should fail if vars for wrong meshes");
      fail_unless(smooth_vars == NULL, "fail should return null");

      // bad radius
      rc = fc_geomSmoothVariable(numMesh, meshes, bins, vars, 0, 0,
                                 &smooth_vars);
      fail_unless(rc != FC_SUCCESS, "should fail if radius is negative");
      fail_unless(smooth_vars == NULL, "fail should return null");
      rc = fc_displacedGeomSmoothVariable(numMesh, meshes, displs, displBins, 
					  vars, 0, 0, &smooth_vars);
      fail_unless(rc != FC_SUCCESS, "should fail if radius is negative");
      fail_unless(smooth_vars == NULL, "fail should return null");

      // bad smooth_vars
      rc = fc_geomSmoothVariable(numMesh, meshes, bins, vars, radius[j], 0,
                                 NULL);
      fail_unless(rc != FC_SUCCESS, "should fail if smooth_vars is NULL");
      rc = fc_displacedGeomSmoothVariable(numMesh, meshes, displs, displBins,
					  vars, radius[j], 0, NULL);
      fail_unless(rc != FC_SUCCESS, "should fail if smooth_vars is NULL");
    }
    
    // cleanup
    for (j = 0; j < numMesh; j++)
      fc_deleteVariable(vars[j]);
  }

  // --- more testing that didn't really fit into loops above

  // create a reference smoothed vars for further testing
  for (i = 0; i < numMesh; i++) {
    fc_createVariable(meshes[i], var_name, &vars[i]);
    rc = fc_setVariableData(vars[i], numElem, 1, FC_AT_ELEMENT,
                            FC_MT_SCALAR, FC_DT_DOUBLE, (void*)datas[i]);
    fail_unless(rc == FC_SUCCESS, 
                "test aborted: failed to create test var");
  }
  rc = fc_geomSmoothVariable(numMesh, meshes, bins, vars, radius[1], 0,
                             &ref_smooth_vars);
  fail_unless(rc == FC_SUCCESS, "failed to smooth");

  // test that other data types give same answer for elements, radius[1]  
  for (i = 0; i < numDataType; i++) {
    // create vars
    for (j = 0; j < numMesh; j++) {
      void* data_dt;
      if (dataTypes[i] == FC_DT_CHAR) {
        char *data_char = malloc(numElem*sizeof(char));
        data_char[0] = 'a';
        for (k = 1; k < numElem; k++)
          data_char[k] = 'z';
        data_dt = data_char;
      }
      else if (dataTypes[i] == FC_DT_INT) {
        int *data_int =  malloc(numElem*sizeof(int));
        data_int[0] = big_values[j];
        for (k = 1; k < numElem; k++)
          data_int[k] = 0;
        data_dt = data_int;
      }
      else if (dataTypes[i] == FC_DT_FLOAT) {
        float *data_float = malloc(numElem*sizeof(float));
        data_float[0] = big_values[j];
        for (k = 1; k < numElem; k++)
          data_float[k] = 0.;
        data_dt = data_float;
      }
      fc_createVariable(meshes[j], var_name, &vars[j]);
      rc = fc_setVariableDataPtr(vars[j], numElem, 1, FC_AT_ELEMENT, 
                                 FC_MT_SCALAR, dataTypes[i], data_dt);
      fail_unless(rc == FC_SUCCESS, 
                  "test aborted: failed to create test var");
    }

    // smooth 
    rc = fc_geomSmoothVariable(numMesh, meshes, bins, vars, radius[1], 0,
                               &smooth_vars);
    if (dataTypes[i] == FC_DT_CHAR)
      fail_unless(rc != FC_SUCCESS, "should fail to smooth char var");
    else {
      fail_unless(rc == FC_SUCCESS, "failed to smooth");
      for (j = 0; j < numMesh; j++) {
        fail_unless(!FC_HANDLE_EQUIV(smooth_vars[j], ref_smooth_vars[j]),
                    "handles should be different");
        rc = fc_getVariableDataPtr(smooth_vars[j], (void**)&temp_data);
        fail_unless(rc == FC_SUCCESS, "failed to get smoothed var data");
        rc = fc_getVariableDataPtr(ref_smooth_vars[j], (void**)&ref_data);
        fail_unless(rc == FC_SUCCESS, "failed to get ref smooth var data");
        fail_unless(!memcmp(temp_data, ref_data, numElem*sizeof(double)),
                    "mismatch of data values");
      }
    }
    free(smooth_vars);
  }

  // test that multicomponent gives same answer for elements, radius[1]  
  // create vars
  for (i = 0; i < numMesh; i++) {
    double* data_multi = malloc(numElem*numComp*sizeof(double));
    for (j = 0; j < numComp; j++)
      data_multi[j] = big_values[i];
    for (j = numComp; j < numElem*numComp; j++)
      data_multi[j] = 0.;
    fc_createVariable(meshes[i], var_name, &vars[i]);
    rc = fc_setVariableDataPtr(vars[i], numElem, numComp, FC_AT_ELEMENT, 
                               FC_MT_VECTOR, FC_DT_DOUBLE, data_multi);
    fail_unless(rc == FC_SUCCESS, 
                "test aborted: failed to create test var");
  }
  // smooth 
  rc = fc_geomSmoothVariable(numMesh, meshes, bins, vars, radius[1], 0,
                             &smooth_vars);
  fail_unless(rc == FC_SUCCESS, "failed to smooth");
  for (i = 0; i < numMesh; i++) {
    fail_unless(!FC_HANDLE_EQUIV(smooth_vars[i], ref_smooth_vars[i]),
                "handles should be different");
    rc = fc_getVariableDataPtr(smooth_vars[i], (void**)&temp_data);
    fail_unless(rc == FC_SUCCESS, "failed to get smoothed var data");
    rc = fc_getVariableDataPtr(ref_smooth_vars[i], (void**)&ref_data);
    fail_unless(rc == FC_SUCCESS, "failed to get ref smooth var data");
    for (j = 0; j < numElem; j++) 
      for (k = 0; k < numComp; k++) 
        fail_unless(temp_data[j*numComp+k] == ref_data[j],
                    "mismatch of data values");
  }
  free(smooth_vars);
  free(ref_smooth_vars);
  // create bins
  for (i = 0; i < numMesh; i++) {
    fc_freeVertexBin(bins[i]);
    fc_freeVertexBin(displBins[i]);
  }
 
  fc_deleteDataset(dataset);
}
END_TEST
  
// Since seq var smooth is just a wrapper, for geomSmooth (which was
// thoroughly tested above) we'll just test a particular case (and compare 
// against geomSmoothVariable answer)
START_TEST(seq_var_smooth)
{
  FC_ReturnCode rc;
  int i, j, k;
  double side_length = 2.3;
  double big_values[2] = { 7.0, 11.0 };
  FC_Dataset dataset;
  FC_Sequence sequence;
  int numMesh = 2, numStep = 5, temp_numStep;
  FC_Mesh meshes[2], badMeshes[2], badMesh = { 999, 999 };
  FC_Variable vars[2], *ref_smooth_vars;
  FC_Variable *seqVars[2], *displVars[2], *badSeqVars[2];
  FC_Variable **smooth_seq_vars, **displ_smooth_seq_vars;
  int numElem = 27, numVert = 64; // 3 x 3 x 3 hex mesh
  FC_Coords bottom1 =  { 0., 0., 0. };
  FC_Coords top1 =     { 3.0*side_length, 3.0*side_length, 3.0*side_length };
  FC_Coords bottom2, top2;
  double* coords[2];
  double datas[2][27*12], displ_buf[64*3];
  double* temp_data, *displ_data, *ref_data;
  FC_VertexBin* bins[2], **displBins[2], *badBins[2], **badDisplBins[2];
  double radius = side_length*2.1;
  double seq_coords[5] = { 0., 1., 2., 3., 4 };

  // setup - create test meshes & sequence
  bottom2[0] = bottom1[0] + side_length*4.0; 
  top2[0] = top1[0] + side_length*4.0;
  for (i = 1; i < 3; i++) {
    bottom2[i] = bottom1[i];
    top2[i] = top1[i];
  }
  fc_createDataset("temp.xxx", &dataset);
  fc_createSequence(dataset, "simple sequence", &sequence);
  rc = fc_setSequenceCoords(sequence, numStep, FC_DT_DOUBLE, seq_coords);
  fail_unless(rc == FC_SUCCESS, "test aborted: failed to create sequence");
  rc = fc_createSimpleHexMesh(dataset, "simple hex mesh1", 3, 3, 3, bottom1, 
			      top1, &meshes[0]);
  fail_unless(rc == FC_SUCCESS, "test aborted: failed to create simple mesh");
  rc = fc_createSimpleHexMesh(dataset, "simple hex mesh2", 3, 3, 3, bottom2, 
			      top2, &meshes[1]);
  fail_unless(rc == FC_SUCCESS, "test aborted: failed to create simple mesh");
  fc_getMeshCoordsPtr(meshes[0], &coords[0]);
  fc_getMeshCoordsPtr(meshes[1], &coords[1]);

  // setup data (same for both meshes, all entity types)
  datas[0][0] = big_values[0];
  datas[1][0] = big_values[1];
  for (i = 1; i < 27*12; i++) {
    datas[0][i] = 0.;
    datas[1][i] = 0.;
  }

  // create seq vars
  for (i = 0; i < numMesh; i++) {
    rc = fc_createSeqVariable(meshes[i], sequence, "displ",
			      &temp_numStep, &displVars[i]);
    fail_unless(rc == FC_SUCCESS, "failed to create seq displ");
    rc = fc_createSeqVariable(meshes[i], sequence, "abitrary var", 
                              &temp_numStep, &seqVars[i]);
    fail_unless(rc == FC_SUCCESS, "failed to create seq var");
    for (j = 0; j < numStep; j++) {
      for (k = 0; k < 3*numVert; k++)
	displ_buf[k] = i; // translate
      if (j > 0) { // make different data each step
        for (k = 0; k < numElem; k++)
          datas[i][k] += 1.;
      }
      rc = fc_setVariableData(displVars[i][j], numVert, 3, FC_AT_VERTEX,
			      FC_MT_VECTOR, FC_DT_DOUBLE, displ_buf);
      fail_unless(rc == FC_SUCCESS,
		  "abort: failed to create displ var");
      rc = fc_setVariableData(seqVars[i][j], numElem, 1, FC_AT_ELEMENT, 
                              FC_MT_SCALAR, FC_DT_DOUBLE, datas[i]);
      fail_unless(rc == FC_SUCCESS, 
                  "test aborted: failed to create test var");
    }
  }

  // create bins
  for (i = 0; i < numMesh; i++) {
    rc = fc_createMeshVertexBin(meshes[i], &bins[i]);
    fail_unless(rc == FC_SUCCESS, "test aborted: failed to create vertex bin");
    displBins[i] = malloc(numStep*sizeof(FC_VertexBin*));
    for (j = 0; j < numStep; j++) {
      rc = fc_createDisplacedMeshVertexBin(meshes[i], displVars[i][j],
					   &displBins[i][j]);
      fail_unless(rc == FC_SUCCESS, "abort: failed to create displ vertex bin");
    }
  }

  // Create reference smoothed vars to test against
  for (i = 0; i < numMesh; i++)
    vars[i] = seqVars[i][0];
  rc = fc_geomSmoothVariable(numMesh, meshes, bins, vars, radius, 0,
			     &ref_smooth_vars);
  fail_unless(rc == FC_SUCCESS, "failed to smooth");
 
  // do it: smooth 
  rc = fc_geomSmoothSeqVariable(numMesh, meshes, bins, numStep, seqVars, 
                                radius, 0, &smooth_seq_vars);
  fail_unless(rc == FC_SUCCESS, "failed to smooth");
  rc = fc_displacedGeomSmoothSeqVariable(numMesh, meshes, numStep, displVars,
					 displBins, numStep, seqVars, radius, 0,
					 &displ_smooth_seq_vars);
  fail_unless(rc == FC_SUCCESS, "failed to smooth on displaced mesh");

  // compare results
  for (i = 0; i < numMesh; i++) {
    for (j = 0; j < numStep; j++) {
      fail_unless(!FC_HANDLE_EQUIV(smooth_seq_vars[i][j], ref_smooth_vars[i]),
                  "handles should be different");
      rc = fc_getVariableDataPtr(smooth_seq_vars[i][j], (void**)&temp_data);
      fail_unless(rc == FC_SUCCESS, "failed to get smoothed var data");
      rc = fc_getVariableDataPtr(displ_smooth_seq_vars[i][j],
				 (void**)&displ_data);
      fail_unless(rc == FC_SUCCESS, "failed to get displaced smooth var data");
      rc = fc_getVariableDataPtr(ref_smooth_vars[i], (void**)&ref_data);
      fail_unless(rc == FC_SUCCESS, "failed to get ref smooth var data");
      for (k = 0; k < numElem; k++) {
        fail_unless(FC_DBL_EQUIV(temp_data[k], ref_data[k]+j), 
                    "mismatch of data values");
	fail_unless(FC_DBL_EQUIV(displ_data[k], ref_data[k]+j),
	    "mismatch of displ data values");
      }
    }
    fc_deleteSeqVariable(numStep, smooth_seq_vars[i]);
    free(smooth_seq_vars[i]);
    fc_deleteSeqVariable(numStep, displ_smooth_seq_vars[i]);
    free(displ_smooth_seq_vars[i]);
  }
  free(smooth_seq_vars);
  free(displ_smooth_seq_vars);

  // --- test bad args

  // bad numMesh
  rc = fc_geomSmoothSeqVariable(0, meshes, bins, numStep, seqVars, 
                                radius, 0, &smooth_seq_vars);
  fail_unless(rc != FC_SUCCESS, "should fail if numMesh < 1");
  fail_unless(smooth_seq_vars == NULL, "fail should return NULL");
  rc = fc_displacedGeomSmoothSeqVariable(0, meshes, numStep, displVars,
				displBins, numStep, seqVars, radius, 0,
					 &displ_smooth_seq_vars);
  fail_unless(rc != FC_SUCCESS, "should fail if numMesh < 1");
  fail_unless(smooth_seq_vars == NULL, "fail should return NULL");

  // bad meshes
  rc = fc_geomSmoothSeqVariable(numMesh, NULL, bins, numStep, seqVars, 
                                radius, 0, &smooth_seq_vars);
  fail_unless(rc != FC_SUCCESS, "should fail if no meshes");
  fail_unless(smooth_seq_vars == NULL, "fail should return nulls");
  rc = fc_displacedGeomSmoothSeqVariable(numMesh, NULL, numStep, displVars,
					 displBins, numStep, seqVars, 
                                radius, 0, &smooth_seq_vars);
  fail_unless(rc != FC_SUCCESS, "should fail if no meshes");
  fail_unless(smooth_seq_vars == NULL, "fail should return nulls");
  badMeshes[0] = meshes[0];
  badMeshes[1] = badMesh;
  rc = fc_geomSmoothSeqVariable(numMesh, badMeshes, bins, numStep, seqVars, 
                                radius, 0, &smooth_seq_vars);
  fail_unless(rc != FC_SUCCESS, "failed to smooth");
  fail_unless(smooth_seq_vars == NULL, "fail should return nulls");
  rc = fc_displacedGeomSmoothSeqVariable(numMesh, badMeshes, numStep,
				  displVars, displBins, numStep, seqVars, 
                                radius, 0, &smooth_seq_vars);
  fail_unless(rc != FC_SUCCESS, "failed to smooth");
  fail_unless(smooth_seq_vars == NULL, "fail should return nulls");

  // bad displ vars
  rc = fc_displacedGeomSmoothSeqVariable(numMesh, badMeshes, numStep+1,
				displVars, displBins, numStep, seqVars, 
                                radius, 0, &smooth_seq_vars);
  fail_unless(rc != FC_SUCCESS, "should fail if wrong numStep");
  fail_unless(smooth_seq_vars == NULL, "fail should return null");
  rc = fc_displacedGeomSmoothSeqVariable(numMesh, meshes, numStep, NULL,
					 displBins, numStep, seqVars, 
                                radius, 0, &smooth_seq_vars);
  fail_unless(rc != FC_SUCCESS, "should fail if no seq vars");
  fail_unless(smooth_seq_vars == NULL, "fail should return null");
  badSeqVars[0] = displVars[0];
  badSeqVars[1] = NULL;
  rc = fc_displacedGeomSmoothSeqVariable(numMesh, meshes, numStep, badSeqVars, 
				displBins, numStep, seqVars, 
                                radius, 0, &smooth_seq_vars);
  fail_unless(rc != FC_SUCCESS, "should fail if displ for wrong meshes");
  fail_unless(smooth_seq_vars == NULL, "fail should return null");
  rc = fc_displacedGeomSmoothSeqVariable(numMesh, meshes, numStep, seqVars, 
				displBins, numStep, seqVars, 
                                radius, 0, &smooth_seq_vars);
  fail_unless(rc != FC_SUCCESS, "should fail if not a proper displ var");
  fail_unless(smooth_seq_vars == NULL, "fail should return null");
   
  // bad bins 
  rc = fc_geomSmoothSeqVariable(numMesh, meshes, NULL, numStep, seqVars, 
                                radius, 0, &smooth_seq_vars);
  fail_unless(rc != FC_SUCCESS, "should fail if no bins");
  fail_unless(smooth_seq_vars == NULL, "fail should return null");
  rc = fc_displacedGeomSmoothSeqVariable(numMesh, meshes, numStep, displVars,
					 NULL, numStep, seqVars, 
                                radius, 0, &smooth_seq_vars);
  fail_unless(rc != FC_SUCCESS, "should fail if no bins");
  fail_unless(smooth_seq_vars == NULL, "fail should return null");
  badBins[0] = bins[1];
  badBins[1] = bins[0];
  badDisplBins[0] = displBins[0];
  badDisplBins[1] = NULL;
  rc = fc_geomSmoothSeqVariable(numMesh, meshes, badBins, numStep, seqVars, 
                                radius, 0, &smooth_seq_vars);
  fail_unless(rc != FC_SUCCESS, "should fail if bins for wrong meshes");
  fail_unless(smooth_seq_vars == NULL, "fail should return null");
  rc = fc_displacedGeomSmoothSeqVariable(numMesh, meshes, numStep, displVars,
					 badDisplBins, numStep, seqVars, 
                                radius, 0, &smooth_seq_vars);
  fail_unless(rc != FC_SUCCESS, "should fail if bad bins");
  fail_unless(smooth_seq_vars == NULL, "fail should return null");
  badDisplBins[0] = displBins[1];
  badDisplBins[1] = displBins[0];
  rc = fc_displacedGeomSmoothSeqVariable(numMesh, meshes, numStep, displVars,
					 badDisplBins, numStep, seqVars, 
                                radius, 0, &smooth_seq_vars);
  fail_unless(rc != FC_SUCCESS, "should fail if bad bins");
  fail_unless(smooth_seq_vars == NULL, "fail should return null");

  // bad vars
  rc = fc_geomSmoothSeqVariable(numMesh, meshes, bins, numStep+1, seqVars, 
                                radius, 0, &smooth_seq_vars);
  fail_unless(rc != FC_SUCCESS, "should fail if wrong numStep");
  fail_unless(smooth_seq_vars == NULL, "fail should return null");
  rc = fc_displacedGeomSmoothSeqVariable(numMesh, meshes, numStep, displVars,
					 displBins, numStep+1, seqVars, 
                                radius, 0, &smooth_seq_vars);
  fail_unless(rc != FC_SUCCESS, "should fail if wrong numStep");
  fail_unless(smooth_seq_vars == NULL, "fail should return null");
  rc = fc_geomSmoothSeqVariable(numMesh, meshes, bins, numStep, NULL, 
                                radius, 0, &smooth_seq_vars);
  fail_unless(rc != FC_SUCCESS, "should fail if no seq vars");
  fail_unless(smooth_seq_vars == NULL, "fail should return null");
  rc = fc_displacedGeomSmoothSeqVariable(numMesh, meshes, numStep, displVars,
					 displBins, numStep, NULL, 
                                radius, 0, &smooth_seq_vars);
  fail_unless(rc != FC_SUCCESS, "should fail if no seq vars");
  fail_unless(smooth_seq_vars == NULL, "fail should return null");
  badSeqVars[0] = seqVars[1];
  badSeqVars[1] = seqVars[0];
  rc = fc_geomSmoothSeqVariable(numMesh, meshes, bins, numStep, badSeqVars, 
                                radius, 0, &smooth_seq_vars);
  fail_unless(rc != FC_SUCCESS, "should fail if seq var for wrong meshes");
  fail_unless(smooth_seq_vars == NULL, "fail should return null");
  rc = fc_displacedGeomSmoothSeqVariable(numMesh, meshes, numStep, displVars,
					 displBins, numStep, badSeqVars, 
                                radius, 0, &smooth_seq_vars);
  fail_unless(rc != FC_SUCCESS, "should fail if seq var for wrong meshes");
  fail_unless(smooth_seq_vars == NULL, "fail should return null");

  // bad radius
  rc = fc_geomSmoothSeqVariable(numMesh, meshes, bins, numStep, seqVars, 
                                0, 0, &smooth_seq_vars);
  fail_unless(rc != FC_SUCCESS, "should fail if radius is zero");
  fail_unless(smooth_seq_vars == NULL, "fail should return null");
  rc = fc_displacedGeomSmoothSeqVariable(numMesh, meshes, numStep, displVars,
					 displBins, numStep, seqVars, 
                                0, 0, &smooth_seq_vars);
  fail_unless(rc != FC_SUCCESS, "should fail if radius is zero");
  fail_unless(smooth_seq_vars == NULL, "fail should return null");

  // bad smooth_seq-vars
  rc = fc_geomSmoothSeqVariable(numMesh, meshes, bins, numStep, seqVars, 
                                radius, 0, NULL);
  fail_unless(rc != FC_SUCCESS, "should fail if smooth_seq_vars is NULL");
  rc = fc_displacedGeomSmoothSeqVariable(numMesh, meshes, numStep, displVars,
					 displBins, numStep, seqVars, 
                                radius, 0, NULL);
  fail_unless(rc != FC_SUCCESS, "should fail if smooth_seq_vars is NULL");

  // cleanup
  for (i = 0; i < numMesh; i++) {
    fc_deleteSeqVariable(numStep, seqVars[i]);
    free(seqVars[i]);
    fc_deleteSeqVariable(numStep, displVars[i]);
    free(displVars[i]);
  }
  free(ref_smooth_vars);

  // free bins
  for (i = 0; i < numMesh; i++) {
    fc_freeVertexBin(bins[i]);
    for (j = 0; j < numStep; j++)
      fc_freeVertexBin(displBins[i][j]);
    free(displBins[i]);
  }

  rc = fc_deleteDataset(dataset);
  fail_unless(rc == FC_SUCCESS, "failed to delete dataset at end of test");
}
END_TEST
  
// Populate the Suite with the tests

Suite *geom_suite(void)
{
  Suite *suite = suite_create("Geometry");

  TCase *tc_general = tcase_create(" - General Geom ");
  TCase *tc_param = tcase_create(" - Parameterization ");
  TCase *tc_bound_box = tcase_create(" - Bounding Boxes ");
  TCase *tc_diameter = tcase_create(" - Set Diameter ");
  TCase *tc_centroid = tcase_create(" - Centroid ");
  TCase *tc_measures = tcase_create(" - Entity Measures ");
  TCase *tc_normals = tcase_create(" - Normals ");
  TCase *tc_deformation = tcase_create(" - Extract Deformation ");
  TCase *tc_proximity = tcase_create(" - Proximity ");
  TCase *tc_geom_neighs = tcase_create(" - Geometric Neighbors ");
  TCase *tc_smooth = tcase_create(" - Variable Smoothing ");

  suite_add_tcase(suite, tc_general);
  tcase_add_checked_fixture(tc_general, geom_setup, geom_teardown);
  tcase_add_test(tc_general, euclid_distances);
  tcase_add_test(tc_general, angle_between);
  tcase_add_test(tc_general, valid_displ);
  tcase_add_test(tc_general, coords_displ);

  suite_add_tcase(suite, tc_param);
  tcase_add_test(tc_param, quad_param);

  suite_add_tcase(suite, tc_bound_box);
  tcase_add_checked_fixture(tc_bound_box, geom_setup, geom_teardown);
  tcase_add_test(tc_bound_box, bb_ops);
  tcase_add_test(tc_bound_box, get_bb);
  tcase_add_test(tc_bound_box, bb_mesh);
  tcase_add_test(tc_bound_box, bb_projection);

  suite_add_tcase(suite, tc_diameter);
  tcase_add_checked_fixture(tc_diameter, geom_setup, geom_teardown);
  tcase_add_test(tc_diameter, get_diameter);

  suite_add_tcase(suite, tc_centroid);
  tcase_add_checked_fixture(tc_centroid, geom_setup, geom_teardown);
  tcase_add_test(tc_centroid, get_centroid);

  suite_add_tcase(suite, tc_measures);
  tcase_add_checked_fixture(tc_measures, geom_setup, geom_teardown);
  tcase_add_test(tc_measures, measure_cores);
  tcase_add_test(tc_measures, length_area_vol);
  tcase_add_test(tc_measures, subset_area);

  suite_add_tcase(suite, tc_normals);
  tcase_add_checked_fixture(tc_normals, geom_setup, geom_teardown);
  tcase_add_test(tc_normals, polygon_normal);
  tcase_add_test(tc_normals, surface_normals);

  suite_add_tcase(suite, tc_deformation);
  tcase_add_checked_fixture(tc_deformation, geom_setup, geom_teardown);
  tcase_add_test(tc_deformation, deform_vertex);
  tcase_add_test(tc_deformation, deform_centroid);

  suite_add_tcase(suite, tc_proximity);
  tcase_add_checked_fixture(tc_proximity, geom_setup, geom_teardown);
  tcase_add_test(tc_proximity, mesh_proximity);

  suite_add_tcase(suite, tc_geom_neighs);
  tcase_add_checked_fixture(tc_geom_neighs, geom_setup, geom_teardown);
  tcase_add_test(tc_geom_neighs, verts_in_bb);
  tcase_add_test(tc_geom_neighs, entities_in_bb);
  tcase_add_test(tc_geom_neighs, verts_in_sphere);
  tcase_add_test(tc_geom_neighs, elems_in_sphere);
  tcase_add_test(tc_geom_neighs, entities_in_sphere);

  suite_add_tcase(suite, tc_smooth);
  tcase_add_checked_fixture(tc_smooth, geom_setup, geom_teardown);
  tcase_add_test(tc_smooth, vertex_smooth);
  tcase_add_test(tc_smooth, element_smooth);
  tcase_add_test(tc_smooth, entity_smooth);
  tcase_add_test(tc_smooth, seq_var_smooth);


  return suite;
}
