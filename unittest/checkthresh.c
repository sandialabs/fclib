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
 * \file checkthresh.c
 * \brief Unit testing of \ref Threshold module
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/unittest/checkthresh.c,v $
 * $Revision: 1.16 $ 
 * $Date: 2006/10/19 03:14:53 $
 *
 * \description
 *
 *    Tests the threshold routines.
 *
 * \modifications
 *    - 11/4/04 WSK, Created.
 */

#include <stdlib.h>
#include <string.h>
#include <check.h>
#include "fc.h"
#include "fcP.h"
#include "checkall.h"

// **** test fixtures

static void thresh_setup(void) {
  FC_ReturnCode rc;
  if (isForking) {
    fc_setLibraryVerbosity(fc_messages);
    rc = fc_initLibrary();
    fail_unless(rc == FC_SUCCESS, "test aborted: failed to init library");
  }
}

static void thresh_teardown(void) {
  FC_ReturnCode rc;
  if (isForking) {
    rc = fc_finalLibrary();
    fail_unless(rc == FC_SUCCESS, "failed to final library at end of test");
  }
}

static int lt (double x, double y) { 
  return fc_lt(x, y, DBL_EPSILON, DBL_MIN); }
static int ltet(double x, double y) { 
  return fc_lteq(x, y, DBL_EPSILON, DBL_MIN); }
static int gt(double x, double y) { 
  return fc_gt(x, y, DBL_EPSILON, DBL_MIN); }
static int gtet(double x, double y) { 
  return fc_gteq(x, y, DBL_EPSILON, DBL_MIN); }
static int eq(double x, double y) { 
  return fc_eq(x, y, DBL_EPSILON, DBL_MIN); }
static int neq(double x, double y) { 
  return fc_neq(x, y, DBL_EPSILON, DBL_MIN); }
static int log_or(int a, int b) { return (a || b); };
static int log_and(int a, int b) { return (a && b); };
static int log_xor(int a, int b) { 
  if ((a && b) || (!a && !b)) 
    return 0;
  else 
    return 1;
}

START_TEST(do_threshold)
{
  FC_ReturnCode rc;
  int i, j, k;
  int temp_num;
  char name[20] = "banana", *temp_name;
  FC_Dataset dataset;
  FC_Mesh mesh, *returnMeshes;
  int numReturnMeshes;
  FC_Subset subset, subset2;
  double whole_data = 5.;
  int numVar = 5; // different assocs: vertex, edge, face, element, whole
  FC_Variable coordsVar, *temp_vars, variables[5], badVar = { 999, 999 };
  FC_Variable glbVar;
  FC_AssociationType assocs[5] = { FC_AT_VERTEX, FC_AT_EDGE, FC_AT_FACE,
				  FC_AT_ELEMENT, FC_AT_WHOLE_MESH };
  FC_AssociationType temp_assoc;
  double *data_ps[5];
  double means[5], sdevs[5];
  int numMember, *memberIDs, numEntity;
  int numOps = 6;
  char opstrings[6][4] =          {  "<", "<=", ">", ">=", "=", "!=" };
  int (*ops[6])(double, double) = {   lt, ltet,  gt, gtet,  eq,  neq };
  int numEdges[3];
  double lengths[3];
  double global_data = 9;

  // setup
  fc_loadDataset("../data/gen_hex.ex2", &dataset);
  rc = fc_createGlobalVariable(dataset, "global var", &glbVar);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create global var");
  rc = fc_setVariableData(glbVar, 1, 1, FC_AT_WHOLE_DATASET, FC_MT_SCALAR,
			  FC_DT_DOUBLE, &global_data);
  fail_unless(rc == FC_SUCCESS, "abort: failed to set global var data");
  rc = fc_getMeshByName(dataset, "hex mesh", &numReturnMeshes,&returnMeshes);
  fail_unless(rc == FC_SUCCESS, "aborted - cant get mesh by name");
  fail_unless(numReturnMeshes == 1, "failed to find unique mesh by name");
  mesh = returnMeshes[0];
  free(returnMeshes);
  fc_getMeshCoordsAsVariable(mesh, &coordsVar);
  fc_createComponentVariables(coordsVar, &temp_num, &temp_vars);
  variables[0] = temp_vars[0];
  free(temp_vars);
  fc_getEdgeLengths(mesh, &variables[1]);
  fc_getFaceAreas(mesh, &variables[2]);
  fc_getRegionVolumes(mesh, &variables[3]);
  rc = fc_createVariable(mesh, "whole var", &variables[4]);
  fail_unless(rc == FC_SUCCESS, "failed to set whole assoc var data");
  rc = fc_setVariableData(variables[4], 1, 1, FC_AT_WHOLE_MESH, FC_MT_SCALAR,
			  FC_DT_DOUBLE, (void*)&whole_data);
  fail_unless(rc == FC_SUCCESS, "failed to set whole assoc var data");
  for (i = 0; i < numVar; i++) {
    fc_getVariableMeanSdev(variables[i], &means[i], &sdevs[i]);
    fc_getVariableDataPtr(variables[i], (void**)&data_ps[i]);
  }
  
  // general tests
  for (i = 0; i < numVar; i++) {
    int saved_members[6];
    for (j = 0; j < numOps; j++) {
      rc = fc_createThresholdSubset(variables[i], opstrings[j], means[i],
				    name, &subset);
      fail_unless(rc == FC_SUCCESS, "failed to create threshold subset");
      fail_unless(fc_isSubsetValid(subset), "didn't create valid subset");
      fc_getSubsetName(subset, &temp_name);
      fail_unless(!strcmp(temp_name, name), "mismatch of name");
      free(temp_name);
      fc_getSubsetAssociationType(subset, &temp_assoc);
      fail_unless(temp_assoc == assocs[i], "mismatch of assoc");
      fc_getSubsetMembersAsArray(subset, &numMember, &memberIDs);
      saved_members[j] = numMember;
      for (k = 0; k < numMember; k++) 
	fail_unless(ops[j](data_ps[i][memberIDs[k]], means[i]), 
		    "members failed to meet criteria");
      free(memberIDs);
      fc_createSubsetComplement(subset, "evil subset", &subset2);
      fc_getSubsetMembersAsArray(subset2, &numMember, &memberIDs);
      for (k = 0; k < numMember; k++) 
	fail_unless(!ops[j](data_ps[i][memberIDs[k]], means[i]), 
		    "nonmembers should fail to meet criteria");
      free(memberIDs);
    }
    // Make sure numbers add up
    fc_getMeshNumEntity(mesh, assocs[i], &numEntity);
    fail_unless(saved_members[0] + saved_members[4] == saved_members[1],
		"LT + EQ should = LTEQ");
    fail_unless(saved_members[2] + saved_members[4] == saved_members[3],
		"GT + EQ should = GTEQ");
    fail_unless(saved_members[0] + saved_members[4] + saved_members[2] ==
		numEntity, "LT + GT + EQ should = Total");
    fail_unless(saved_members[0] + saved_members[2] == saved_members[5],
		"LT + GT should = NEQ");
  }
 
  // specific test - I know the edge lengths
  // Couldn't test 'eq' versions because of round off error
  lengths[0] = 8./12.;  // .66667 ...
  numEdges[0] = 6*11*12 + 12*6 + 12*11 + 12;
  lengths[1] = 5./6.;   // .83333 ...
  numEdges[1] = 6*11*12 + 6*11 + 6*12 + 6;
  lengths[2] = 10./11.; // .90909 ...
  numEdges[2] = 6*11*12 + 11*6 + 11*12 + 11;
  rc = fc_createThresholdSubset(variables[1], "<", 0.8, name, &subset);
  fail_unless(rc == FC_SUCCESS, "failed to create threshold subset");
  fc_getSubsetNumMember(subset, &numMember);
  fail_unless(numMember == numEdges[0], "unexpected number of members");
  fc_deleteSubset(subset);
  rc = fc_createThresholdSubset(variables[1], ">", .9, name, &subset);
  fail_unless(rc == FC_SUCCESS, "failed to create threshold subset");
  fc_getSubsetNumMember(subset, &numMember);
  fail_unless(numMember == numEdges[2], "unexpected number of members");
  fc_deleteSubset(subset);
 
  // --- test error conditions (bad args)
  rc = fc_createThresholdSubset(badVar, "<", 0.8, name, &subset);
  fail_unless(rc != FC_SUCCESS, "should fail if bad variable");
  fail_unless(FC_HANDLE_EQUIV(subset, FC_NULL_SUBSET), 
	      "failure should return null");
  rc = fc_createThresholdSubset(glbVar, "<", 0.8, name, &subset);
  fail_unless(rc != FC_SUCCESS, "should fail if global variable");
  fail_unless(FC_HANDLE_EQUIV(subset, FC_NULL_SUBSET), 
	      "failure should return null");
  rc = fc_createThresholdSubset(coordsVar, "<", 0.8, name, &subset);
  fail_unless(rc != FC_SUCCESS, "should fail if variable is not scalar");
  fail_unless(FC_HANDLE_EQUIV(subset, FC_NULL_SUBSET), 
	      "failure should return null");
  rc = fc_createThresholdSubset(variables[1], "badop", 0.8, name, &subset);
  fail_unless(rc != FC_SUCCESS, "shoudl fail if bad operation string");
  fail_unless(FC_HANDLE_EQUIV(subset, FC_NULL_SUBSET), 
	      "failure should return null");
  rc = fc_createThresholdSubset(variables[1], "<", 0.8, NULL, &subset);
  fail_unless(rc != FC_SUCCESS, "shoudl fail if NULL name");
  fail_unless(FC_HANDLE_EQUIV(subset, FC_NULL_SUBSET), 
	      "failure should return null");
  rc = fc_createThresholdSubset(variables[1], "<", 0.8, name, NULL);
  fail_unless(rc != FC_SUCCESS, "shoudl fail if NULL subset");

  rc = fc_deleteDataset(dataset);
  fail_unless(rc == FC_SUCCESS, "failed to delete dataset at end of tests");
}
END_TEST
  
// since threshold range just wraps up other functionality
// (fc_createThresholdSubset and fc_createSubsetIntersection)
// don't have to do exhaustive test
START_TEST(thresh_range)
{
  FC_ReturnCode rc;
  int i, j, k, m;
  char name[20] = "banana", *temp_name;
  FC_Dataset dataset;
  FC_Mesh mesh, *returnMeshes;
  int numReturnMeshes;
  FC_Subset subset, subset2;
  FC_Variable coordsVar, edge_lengths, badVar = { 999, 999 };
  FC_Variable glbVar;
  double *data_p, ref1, ref2;
  int numMember, *memberIDs;
  int numOps = 6, numLogOps = 3;
  char opstrings[6][4] =          {  "<", "<=", ">", ">=", "=", "!=" };
  int (*ops[6])(double, double) = {   lt, ltet,  gt, gtet,  eq,  neq };
  char logOpStrings[3][4] =    { "OR", "AND", "XOR" };
  int (*logOps[3])(int, int) = { log_or, log_and, log_xor };
  int numEdges[3], totalNumEdge;
  double lengths[3];
  double global_data = 9;

  // setup
  fc_loadDataset("../data/gen_hex.ex2", &dataset);
  rc = fc_createGlobalVariable(dataset, "global var", &glbVar);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create global var");
  rc = fc_setVariableData(glbVar, 1, 1, FC_AT_WHOLE_DATASET, FC_MT_SCALAR,
			  FC_DT_DOUBLE, &global_data);
  fail_unless(rc == FC_SUCCESS, "abort: failed to set global var data");
  rc = fc_getMeshByName(dataset, "hex mesh", &numReturnMeshes,&returnMeshes);
  fail_unless(rc == FC_SUCCESS, "aborted - cant get mesh by name");
  fail_unless(numReturnMeshes == 1, "failed to find unique mesh by name");
  mesh = returnMeshes[0];
  free(returnMeshes);
  fc_getMeshCoordsAsVariable(mesh, &coordsVar);
  fc_getEdgeLengths(mesh, &edge_lengths);
  fc_getVariableDataPtr(edge_lengths, (void**)&data_p);
  lengths[0] = 8./12.;  // .66667 ...
  numEdges[0] = 6*11*12 + 12*6 + 12*11 + 12;
  lengths[1] = 5./6.;   // .83333 ...
  numEdges[1] = 6*11*12 + 6*11 + 6*12 + 6;
  lengths[2] = 10./11.; // .90909 ...
  numEdges[2] = 6*11*12 + 11*6 + 11*12 + 11;
  ref1 = 0.8;
  ref2 = 0.9;
  fc_getMeshNumEdge(mesh, &totalNumEdge);
  fail_unless(numEdges[0] + numEdges[1] + numEdges[2] == totalNumEdge,
	      "abort: didn't calc numEdges correctly");

  // specific tests
  rc = fc_createThresholdRangeSubset(edge_lengths, ">", ref1, "&&",
				     "<", ref2, name, &subset);
  fail_unless(rc == FC_SUCCESS, "failed to create threshold subset");
  fc_getSubsetNumMember(subset, &numMember);
  fail_unless(numMember == numEdges[1], "unexpected number of members");
  fc_deleteSubset(subset);
  rc = fc_createThresholdRangeSubset(edge_lengths, "<", ref1, "&&",
				     ">", ref2, name, &subset);
  fail_unless(rc == FC_SUCCESS, "failed to create threshold subset");
  fc_getSubsetNumMember(subset, &numMember);
  fail_unless(numMember == 0, "unexpected number of members");
  fc_deleteSubset(subset);
  rc = fc_createThresholdRangeSubset(edge_lengths, ">", ref1, "||",
				     "<", ref2, name, &subset);
  fail_unless(rc == FC_SUCCESS, "failed to create threshold subset");
  fc_getSubsetNumMember(subset, &numMember);
  fail_unless(numMember == totalNumEdge, "unexpected number of members");
  fc_deleteSubset(subset);
  rc = fc_createThresholdRangeSubset(edge_lengths, "<", ref1, "||",
				     ">", ref2, name, &subset);
  fail_unless(rc == FC_SUCCESS, "failed to create threshold subset");
  fc_getSubsetNumMember(subset, &numMember);
  fail_unless(numMember == numEdges[0]+numEdges[2], 
	      "unexpected number of members");
  fc_deleteSubset(subset);
  rc = fc_createThresholdRangeSubset(edge_lengths, ">", ref1, "XOR",
				     "<", ref2, name, &subset);
  fail_unless(rc == FC_SUCCESS, "failed to create threshold subset");
  fc_getSubsetNumMember(subset, &numMember);
  fail_unless(numMember == numEdges[0]+numEdges[2], 
	      "unexpected number of members");
  fc_deleteSubset(subset);
  rc = fc_createThresholdRangeSubset(edge_lengths, "<", ref1, "XOR",
				     ">", ref2, name, &subset);
  fail_unless(rc == FC_SUCCESS, "failed to create threshold subset");
  fc_getSubsetNumMember(subset, &numMember);
  fail_unless(numMember == numEdges[0]+numEdges[2], 
	      "unexpected number of members");
  fc_deleteSubset(subset);

  // general tests
  for (i = 0; i < numLogOps; i++) {
    for (j = 0; j < numOps; j++) {
      for (k = 0; k < numOps; k++) {
	rc = fc_createThresholdRangeSubset(edge_lengths, opstrings[j], ref1,
				   logOpStrings[i], opstrings[k], ref2,
				   name, &subset);
	fail_unless(rc == FC_SUCCESS, "failed to create threshold subset");
	fail_unless(fc_isSubsetValid(subset), "didn't create valid subset");
	fc_getSubsetName(subset, &temp_name);
	fail_unless(!strcmp(temp_name, name), "mismatch of name");
	free(temp_name);
	fc_getSubsetMembersAsArray(subset, &numMember, &memberIDs);
	for (m = 0; m < numMember; m++) {
	  int test1, test2;
	  test1 = ops[j](data_p[memberIDs[m]], ref1);
	  test2 = ops[k](data_p[memberIDs[m]], ref2); 
	  fail_unless(logOps[i](test1, test2), 
		      "members failed to meet criteria");
	}
	free(memberIDs);
	fc_createSubsetComplement(subset, "evil subset", &subset2);
	fc_getSubsetMembersAsArray(subset2, &numMember, &memberIDs);
	for (m = 0; m < numMember; m++) {
	  int test1, test2;
	  test1 = ops[j](data_p[memberIDs[m]], ref1);
	  test2 = ops[k](data_p[memberIDs[m]], ref2); 
	  fail_unless(!logOps[i](test1, test2), 
		      "non-members failed to NOT meet criteria");
	}
	  free(memberIDs);
      }
    }
  }

  // --- test error conditions (bad args)
  rc = fc_createThresholdRangeSubset(badVar, ">", ref1, "&&",
				     "<", ref2, name, &subset);
  fail_unless(rc != FC_SUCCESS, "should fail if bad variable");
  fail_unless(FC_HANDLE_EQUIV(subset, FC_NULL_SUBSET), 
	      "failure should return null");
  rc = fc_createThresholdRangeSubset(glbVar, ">", ref1, "&&",
				     "<", ref2, name, &subset);
  fail_unless(rc != FC_SUCCESS, "should fail if global variable");
  fail_unless(FC_HANDLE_EQUIV(subset, FC_NULL_SUBSET), 
	      "failure should return null");
  rc = fc_createThresholdRangeSubset(coordsVar, ">", ref1, "&&",
				     "<", ref2, name, &subset); 
  fail_unless(rc != FC_SUCCESS, "should fail if variable is not scalar");
  fail_unless(FC_HANDLE_EQUIV(subset, FC_NULL_SUBSET), 
	      "failure should return null");
  rc = fc_createThresholdRangeSubset(coordsVar, "badOp", ref1, "&&",
				     "<", ref2, name, &subset); 
  fail_unless(rc != FC_SUCCESS, "should fail if bad operation string");
  fail_unless(FC_HANDLE_EQUIV(subset, FC_NULL_SUBSET), 
	      "failure should return null");
  rc = fc_createThresholdRangeSubset(coordsVar, ">", ref1, "&&",
				     "badOp", ref2, name, &subset); 
  fail_unless(rc != FC_SUCCESS, "should fail if bad operation string");
  fail_unless(FC_HANDLE_EQUIV(subset, FC_NULL_SUBSET), 
	      "failure should return null");
  rc = fc_createThresholdRangeSubset(coordsVar, ">", ref1, "badLogOp",
				     "<", ref2, name, &subset); 
  fail_unless(rc != FC_SUCCESS, "should fail if bad logical operation str"); 
  fail_unless(FC_HANDLE_EQUIV(subset, FC_NULL_SUBSET), 
	      "failure should return null");
  rc = fc_createThresholdRangeSubset(coordsVar, ">", ref1, "&&",
				     "<", ref2, NULL, &subset); 
  fail_unless(rc != FC_SUCCESS, "should fail if NULL name"); 
  fail_unless(FC_HANDLE_EQUIV(subset, FC_NULL_SUBSET), 
	      "failure should return null");
  rc = fc_createThresholdRangeSubset(coordsVar, ">", ref1, "&&",
				     "<", ref2, name, NULL); 
  fail_unless(rc != FC_SUCCESS, "should fail if NULL subset");

  rc = fc_deleteDataset(dataset);
  fail_unless(rc == FC_SUCCESS, "failed to delete dataset at end of tests");
}
END_TEST
  
// Populate the Suite with the tests

Suite *thresh_suite(void)
{
  Suite *suite = suite_create("Threshold");

  TCase *tc_threshold = tcase_create(" - Threshold ");

  suite_add_tcase(suite, tc_threshold);
  tcase_add_checked_fixture(tc_threshold, thresh_setup, thresh_teardown);
  tcase_add_test(tc_threshold, do_threshold);
  tcase_add_test(tc_threshold, thresh_range);

  return suite;
}
