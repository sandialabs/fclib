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
 * \file checktrack.c
 * \brief Unit testing of \ref FeatureTracking module.
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/unittest/checktrack.c,v $
 * $Revision: 1.13 $ 
 * $Date: 2006/09/22 04:40:45 $
 */
#include <string.h>
#include <stdlib.h>
#include <check.h>
#include "fc.h"
#include "fcP.h"
#include "checkall.h"

// **** global flags

// **** test fixtures

static void track_setup(void) {
  FC_ReturnCode rc;
  if (isForking) {
    fc_setLibraryVerbosity(fc_messages);
    rc = fc_initLibrary();
    fail_unless(rc == FC_SUCCESS, "test aborted: failed to init library");
  }
}

static void track_teardown(void) {
  FC_ReturnCode rc;
  if (isForking) {
    rc = fc_finalLibrary();
    fail_unless(rc == FC_SUCCESS, "failed to final library at end of test");
  }
}

START_TEST(default_scores)
{
  FC_ReturnCode rc;
  int i, j;
  FC_Dataset dataset;
  FC_Mesh mesh, *returnMeshes;
  int numReturnMeshes;
  int numSubset1 = 2, numSubset2 = 4;
  FC_Subset subsets1[2], subsets2[4];
  int numMemberSubset1[2] = { 3, 5 };
  int numMemberSubset2[4] = { 2, 3, 4, 5 };
  int membersSubset1[2][5] = { { 0, 1, 2 }, { 11, 12, 13, 14, 15 } };
  int membersSubset2[4][5] = { { 0, 1 }, { 3, 4, 5 }, { 2, 11, 12, 13 }, 
			       { 15, 16, 17, 18, 19 } };
  double good_scores[2][4] = { { 2, 0, 1, 0 }, { 0, 0, 3, 1 } };
  double **scores;

  // get a mesh
  rc = fc_loadDataset("../data/gen_hex.ex2", &dataset);
  fail_unless(rc == FC_SUCCESS, "abort: failed to load dataset");
  rc = fc_getMeshByName(dataset, "hex mesh", &numReturnMeshes,&returnMeshes);
  fail_unless(rc == FC_SUCCESS, "aborted - cant get mesh by name");
  fail_unless(numReturnMeshes == 1, "failed to find unique mesh by name");
  mesh = returnMeshes[0];
  free(returnMeshes);

  // create some subsets
  for (i = 0; i < numSubset1; i++) {
    rc = fc_createSubset(mesh, "temp subset", FC_AT_VERTEX, &subsets1[i]);
    fail_unless(rc == FC_SUCCESS, "abort: failed to create subset");
    for (j = 0; j < numMemberSubset1[i]; j++)
      fc_addMemberToSubset(subsets1[i], membersSubset1[i][j]);
  }
  for (i = 0; i < numSubset2; i++) {
    rc = fc_createSubset(mesh, "temp subset", FC_AT_VERTEX, &subsets2[i]);
    fail_unless(rc == FC_SUCCESS, "abort: failed to create subset");
    for (j = 0; j < numMemberSubset2[i]; j++)
      fc_addMemberToSubset(subsets2[i], membersSubset2[i][j]);
  }

  // test
  rc = fc_computeNumOverlapScore(numSubset1, subsets1, numSubset2, subsets2,
			       &scores);
  fail_unless(rc == FC_SUCCESS, "failed to compute scores");
  for (i = 0; i < numSubset1; i++) {
    for (j = 0; j < numSubset2; j++)
      fail_unless(scores[i][j] == good_scores[i][j], "mismatch of score");
    free(scores[i]);
  }
  free(scores);

  // test if numPrev or numNew are zero
  rc = fc_computeNumOverlapScore(numSubset1, subsets1, 0, NULL, &scores);
  fail_unless(rc == FC_SUCCESS, "failed to compute scores");
  fail_unless(scores == NULL, "if no scores, scores should be NULL");
  rc = fc_computeNumOverlapScore(0, NULL, numSubset2, subsets2, &scores);
  fail_unless(rc == FC_SUCCESS, "failed to compute scores");
  fail_unless(scores == NULL, "if no scores, scores should be NULL");

  // test bad args
  rc = fc_computeNumOverlapScore(-1, subsets1, numSubset2, subsets2, &scores);
  fail_unless(rc != FC_SUCCESS, "should fail if bad numPrev");
  fail_unless(scores == NULL, "fail should return NULLs");
  rc = fc_computeNumOverlapScore(numSubset1, NULL, numSubset2, subsets2, 
			       &scores);
  fail_unless(rc != FC_SUCCESS, "should fail if NULL prevSubsets");
  fail_unless(scores == NULL, "fail should return NULLs");
  rc = fc_computeNumOverlapScore(numSubset1, subsets1, -1, subsets2, &scores);
  fail_unless(rc != FC_SUCCESS, "should fail if bad numNew");
  fail_unless(scores == NULL, "fail should return NULLs");
  rc = fc_computeNumOverlapScore(numSubset1, subsets1, numSubset2, NULL,
			       &scores); 
  fail_unless(rc != FC_SUCCESS, "should fail if bad NULL newSubsets");
  fail_unless(scores == NULL, "fail should return NULLs");
  rc = fc_computeNumOverlapScore(numSubset1, subsets1, numSubset2, subsets2, 
			       NULL);
  fail_unless(rc != FC_SUCCESS, "should fail if NULL return arg");

  // all done
  rc = fc_deleteDataset(dataset);
  fail_unless(rc == FC_SUCCESS, "deleting dataset at end of test failed");
}
END_TEST

START_TEST(default_corresponds)
{
  FC_ReturnCode rc;
  int i, j;
  int numSubset1 = 2, numSubset2 = 4;
  // rint rounds *.5 up if * is odd and down if * is even
  double scores[2][4] = { { 1.1, 1.5, 3.4, -11.11 }, 
			  { 0.5, -0.5, -1.5, .99 } };
  double *scores_p[2] =  { scores[0], scores[1] } ;
  int good_corresponds[2][4] = { { 1, 2, 3, -11 }, 
				 { 0, 0, -2, 1 } };
  int **corresponds;

  // test
  rc = fc_computeSimpleCorrespond(numSubset1, numSubset2, scores_p,
					&corresponds);
  fail_unless(rc == FC_SUCCESS, "failed to compute corresponds");
  for (i = 0; i < numSubset1; i++) {
    for (j = 0; j < numSubset2; j++) 
      fail_unless(corresponds[i][j] == good_corresponds[i][j], 
		  "mismatch of score");
    free(corresponds[i]);
  }
  free(corresponds);
  
  // test if numPrev or numNew are zero
  rc = fc_computeSimpleCorrespond(numSubset1, 0, scores_p,
					&corresponds); 
  fail_unless(rc == FC_SUCCESS, "failed to compute correspondences");
  fail_unless(corresponds == NULL, 
	      "if no correspondences, corresponds should be NULL");
  rc = fc_computeSimpleCorrespond(0, numSubset2, scores_p,
					&corresponds); 
  fail_unless(rc == FC_SUCCESS, "failed to compute correspondences");
  fail_unless(corresponds == NULL, 
	      "if no correspondences, corresponds should be NULL");
  
  // test bad args
  rc = fc_computeSimpleCorrespond(-1, numSubset2, scores_p,
					&corresponds); 
  fail_unless(rc != FC_SUCCESS, "should fail if bad numPrev");
  fail_unless(corresponds == NULL, "fail should return NULLs");
  rc = fc_computeSimpleCorrespond(numSubset1, -1, scores_p, 
					&corresponds);
  fail_unless(rc != FC_SUCCESS, "should fail if bad numNew");
  fail_unless(corresponds == NULL, "fail should return NULLs");
  rc = fc_computeSimpleCorrespond(numSubset1, numSubset2, NULL,
					&corresponds); 
  fail_unless(rc != FC_SUCCESS, "should fail if NULL scores");
  fail_unless(corresponds == NULL, "fail should return NULLs");
  rc = fc_computeSimpleCorrespond(numSubset1, numSubset2, scores_p, 
					NULL);
  fail_unless(rc != FC_SUCCESS, "should fail if NULL return arg");
  scores_p[0] = NULL;
  rc = fc_computeSimpleCorrespond(numSubset1, numSubset2, scores_p,
					&corresponds); 
  fail_unless(rc != FC_SUCCESS, "should fail if bad scores");
  fail_unless(corresponds == NULL, "fail should return NULLs");
}
END_TEST

/// \FIX Add a 'split-merge' event to see if handled correctly
/// \FIX Change so that splitting & merging more than 2 things
START_TEST(default_addROIs)
{
  FC_ReturnCode rc;
  int i, j;
  FC_FeatureGroup *group;
  int numPrev;
  int temp_numStep, temp_numFeature, temp_numROI, *temp_stepIDs;
  FC_Subset *temp_rois;
  // setup these events: none, two, same two, merge, split, none, one 
  // Creating this structure 
  //
  //      stepID:  0    1    2    3    4    5    6
  //   -------------------------------------------
  //   featureID:       0 .. 0 ->   -> 3
  //                             \ /
  //                              2  
  //                             / \                 .
  //                    1 .. 1 ->   -> 4         5
  //
  int numStep = 7, numROI = 8, numFeature = 6;
  int numSegments[7] = { 0, 2, 2, 1, 2, 0, 1 };
  FC_Subset segments[7][2] = { {                    }, 
			       { { 0, 0 }, { 1, 1 } }, 
			       { { 2, 2 }, { 3, 3 } }, 
			       { { 4, 4 }           },
			       { { 5, 5 }, { 6, 6 } }, 
			       {                    }, 
			       { { 7, 7 }           }  };
  int corresponds[7][2][2] = { { },
			       { },
			       { { 1, 0 }, { 0, 1 } },
			       { { 1 }, { 1 } },
			       { { 1, 1 } },
			       { },
			       { } };
  int *corresponds_p[2], **corresponds_pp;
  int numROIPerFeature[6] = { 2, 2, 1, 1, 1, 1 };
  int stepIDsPerFeature[6][2] = { { 1, 2 }, { 1, 2 }, { 3 }, { 4 }, { 4 },
				  { 6 } };
  FC_Subset roisPerFeature[6][2] = { { segments[1][0], segments[2][0] },
				     { segments[1][1], segments[2][1] },
				     { segments[3][0] },
				     { segments[4][0] },
				     { segments[4][1] },
				     { segments[6][0] } };
  // Create a group
  rc = fc_createFeatureGroup(&group);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create feature group");

  // Add the rois
  for (i = 0; i < numStep; i++) {
    if (i == 0) 
      numPrev = 0;
    else 
      numPrev = numSegments[i-1];
    if (numPrev == 0 || numSegments[i] == 0) 
      corresponds_pp = NULL;
    else {
      corresponds_pp = corresponds_p;
      corresponds_p[0] = corresponds[i][0];
      corresponds_p[1] = corresponds[i][1];
    }
    rc = fc_addROIsUsingBasicEvents(i, numPrev, numSegments[i], segments[i],
				     corresponds_pp, group);
    fail_unless(rc == FC_SUCCESS, "failed to add rois");
  }

  // check the group
  fc_featureGroup_getNumStep(group, &temp_numStep);
  fail_unless(temp_numStep == numStep, "mismatch of numStep");
  fc_featureGroup_getNumFeature(group, &temp_numFeature);
  fail_unless(temp_numFeature == numFeature, "mismatch of numFeature");
  fc_featureGroup_getNumROI(group, &temp_numROI);
  fail_unless(temp_numROI == numROI, "mismatch of numROI");

  // check the features
  for (i = 0; i < numFeature; i++) {
    fc_getFeatureROIs(group, i, &temp_numROI, &temp_stepIDs, &temp_rois);
    fail_unless(temp_numROI == numROIPerFeature[i], "mismatch of numROI");
    for (j = 0; j < temp_numROI; j++) {
      fail_unless(temp_stepIDs[j] == stepIDsPerFeature[i][j],
		  "mismatch of step ID");
      fail_unless(FC_HANDLE_EQUIV(temp_rois[j], roisPerFeature[i][j]),
		  "mismatch of roi");
    }
    free(temp_stepIDs);
    free(temp_rois);
  }

  // test bad args

  // setup one more step
  numPrev = numSegments[numStep-1];
  corresponds_p[0] = corresponds[3][0]; // this is { { 1 } }
  // This would work
  //rc = fc_addROIsUsingBasicEvents(numStep, numPrev, numSegments[numStep-1],
  //			   segments[numStep-1], corresponds_p, group);
  //fail_unless(rc == FC_SUCCESS, "failed to add rois");

  rc = fc_addROIsUsingBasicEvents(-1, numPrev, numSegments[numStep-1],
				   segments[numStep-1], corresponds_p, group);
  fail_unless(rc != FC_SUCCESS, "should fail if bad numStep");
  rc = fc_addROIsUsingBasicEvents(numStep, -1, numSegments[numStep-1],
				   segments[numStep-1], corresponds_p, group);
  fail_unless(rc != FC_SUCCESS, "should fail if bad numPrev");
  rc = fc_addROIsUsingBasicEvents(numStep, numPrev+1, numSegments[numStep-1],
				   segments[numStep-1], corresponds_p, group);
  fail_unless(rc != FC_SUCCESS, "should fail if bad numPrev");
  rc = fc_addROIsUsingBasicEvents(numStep, numPrev, -1,
				   segments[numStep-1], corresponds_p, group);
  fail_unless(rc != FC_SUCCESS, "should fail if bad numNew");
  rc = fc_addROIsUsingBasicEvents(numStep, numPrev, numSegments[numStep-1],
				   NULL, corresponds_p, group);
  fail_unless(rc != FC_SUCCESS, "should fail if NULL segments & num > 0");
  rc = fc_addROIsUsingBasicEvents(numStep, numPrev, numSegments[numStep-1],
				   segments[numStep-1], NULL, group);
  fail_unless(rc != FC_SUCCESS, "should fail if NULL matrix & nums > 0");
   rc = fc_addROIsUsingBasicEvents(numStep, numPrev, numSegments[numStep-1],
				   segments[numStep-1], corresponds_p, NULL);
  fail_unless(rc != FC_SUCCESS, "should fail if bad group");

  fc_freeFeatureGroup(group);
}
END_TEST

START_TEST(track_tests);
{
  FC_ReturnCode rc;
  FC_Dataset dataset;
  FC_Mesh mesh, *returnMeshes;
  int numReturnMeshes;
  
  int i, j, k, i_case;
  int numCase = 4;
  FC_FeatureGroup *group;
  int temp_numStep, temp_numFeature, temp_numROI, *temp_stepIDs;
  FC_Subset *temp_rois;
  // setup these events: none, two, same two, merge, split, none, one 
  // Creating this structure 
  //
  //      stepID:  0    1    2    3    4    5    6
  //   -------------------------------------------
  //   featureID:       0 .. 0 ->   -> 3
  //                             \ /
  //                              2  
  //                             / \                 .
  //                    1 .. 1 ->   -> 4         5
  //
  int numStep = 7, numROI = 8, numFeature = 6;
  int numSegments[7] = { 0, 2, 2, 1, 2, 0, 1 };
  int numMember = 3; // there's going to be 3 members in each subset
  int memberIDs[7][2][3] = { { {}, {} }, 
			     { { 1, 2, 3 }, { 5, 6, 7 } }, 
			     { { 0, 1, 2 }, { 7, 8, 9 } }, 
			     { { 2, 5, 7 }, {} }, 
			     { { 0, 1, 2 }, { 7, 8, 9} }, 
			     { {}, {} }, 
			     { { 22, 23, 24 }, {} } };
  FC_Subset segments[7][2];
  FC_Subset *segments_p[7] = { segments[0], segments[1], segments[2], 
			       segments[3], segments[4], segments[5],
			       segments[6] };
  int numROIPerFeature[6] = { 2, 2, 1, 1, 1, 1 };
  int stepIDsPerFeature[6][2] = { { 1, 2 }, { 1, 2 }, { 3 }, { 4 }, { 4 },
  			  { 6 } };
  FC_Subset *roisPerFeature[6][2] = { { &segments[1][0], &segments[2][0] },
				      { &segments[1][1], &segments[2][1] },
				      { &segments[3][0] },
				      { &segments[4][0] },
				      { &segments[4][1] },
				      { &segments[6][0] } };

  // --- setup

  // get a mesh
  rc = fc_loadDataset("../data/gen_hex.ex2", &dataset);
  fail_unless(rc == FC_SUCCESS, "abort: failed to load dataset");
  rc = fc_getMeshByName(dataset, "hex mesh", &numReturnMeshes,&returnMeshes);
  fail_unless(rc == FC_SUCCESS, "aborted - cant get mesh by name");
  fail_unless(numReturnMeshes == 1, "failed to find unique mesh by name");
  mesh = returnMeshes[0];
  free(returnMeshes);

  // create some segments
  for (i = 0; i < numStep; i++) {
    for (j = 0; j < numSegments[i]; j++) {
      rc = fc_createSubset(mesh, "temp subset", FC_AT_VERTEX, &segments[i][j]);
      fail_unless(rc == FC_SUCCESS, "abort: failed to create subset");
      for (k = 0; k < numMember; k++)
	fc_addMemberToSubset(segments[i][j], memberIDs[i][j][k]);
    }
  }

  // --- Test all four ways of using track routines (should give same results)

  for (i_case = 0; i_case < numCase; i_case++) {

    // Create a group
    rc = fc_createFeatureGroup(&group);
    fail_unless(rc == FC_SUCCESS, "abort: failed to create feature group");
    
    // Do it - different for each case
    if (i_case == 0) {
      for (i = 0; i < numStep; i++) {
	rc = fc_trackStep_custom(i, numSegments[i], segments[i], group,
				 fc_computeNumOverlapScore,
				 fc_computeSimpleCorrespond, 
				 fc_addROIsUsingBasicEvents);
	fail_unless(rc == FC_SUCCESS, "failed to track step (custom)");
      }
    }
    else if (i_case == 1) {
      for (i = 0; i < numStep; i++) {
	rc = fc_trackStep(i, numSegments[i], segments[i], group);
 	fail_unless(rc == FC_SUCCESS, "failed to track step");
      }
    }
    else if (i_case == 2) {
      rc = fc_trackAllSteps_custom(numStep, numSegments, segments_p, group,
				 fc_computeNumOverlapScore,
				 fc_computeSimpleCorrespond, 
				 fc_addROIsUsingBasicEvents);
 	fail_unless(rc == FC_SUCCESS, "failed to track all steps (custom)");
    }
    else if (i_case == 3) {
      rc = fc_trackAllSteps(numStep, numSegments, segments_p, group);
 	fail_unless(rc == FC_SUCCESS, "failed to track all steps");
    }
    
    // check the group
    fc_featureGroup_getNumStep(group, &temp_numStep);
    fail_unless(temp_numStep == numStep, "mismatch of numStep");
    fc_featureGroup_getNumFeature(group, &temp_numFeature);
    fail_unless(temp_numFeature == numFeature, "mismatch of numFeature");
    fc_featureGroup_getNumROI(group, &temp_numROI);
    fail_unless(temp_numROI == numROI, "mismatch of numROI");
    
    // check the features
    for (i = 0; i < numFeature; i++) {
      fc_getFeatureROIs(group, i, &temp_numROI, &temp_stepIDs, &temp_rois);
      fail_unless(temp_numROI == numROIPerFeature[i], "mismatch of numROI");
      for (j = 0; j < temp_numROI; j++) {
	fail_unless(temp_stepIDs[j] == stepIDsPerFeature[i][j],
		    "mismatch of step ID");
	fail_unless(FC_HANDLE_EQUIV(temp_rois[j], *roisPerFeature[i][j]),
		  "mismatch of roi");
      }
      free(temp_stepIDs);
      free(temp_rois);
    }
    
    // cleanup
    fc_freeFeatureGroup(group);
  }

  // --- test error conditions

  // Create a group
  rc = fc_createFeatureGroup(&group);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create feature group");
  
  // fc_trackStep_custom() -- bad args
  rc = fc_trackStep_custom(-1, numSegments[0], segments[0], group,
			   fc_computeNumOverlapScore,
			   fc_computeSimpleCorrespond,
			   fc_addROIsUsingBasicEvents);
  fail_unless(rc != FC_SUCCESS, "should fail if bad step");
  rc = fc_trackStep_custom(1, numSegments[0], segments[0], group,
			   fc_computeNumOverlapScore,
			   fc_computeSimpleCorrespond,
			   fc_addROIsUsingBasicEvents);
  fail_unless(rc != FC_SUCCESS, "should fail if not next step ID");
  rc = fc_trackStep_custom(0, -1, segments[0], group,
			   fc_computeNumOverlapScore,
			   fc_computeSimpleCorrespond,
			   fc_addROIsUsingBasicEvents);
  fail_unless(rc != FC_SUCCESS, "should fail if bad numSegment");
  rc = fc_trackStep_custom(0, 1, NULL, group,
			   fc_computeNumOverlapScore,
			   fc_computeSimpleCorrespond,
			   fc_addROIsUsingBasicEvents);
  fail_unless(rc != FC_SUCCESS, 
	      "should fail if numSegment > 0 & null segments");
  rc = fc_trackStep_custom(0, numSegments[0], segments[0], NULL,
			   fc_computeNumOverlapScore,
			   fc_computeSimpleCorrespond,
			   fc_addROIsUsingBasicEvents);
  fail_unless(rc != FC_SUCCESS, "should fail if bad group");

  // fc_trackAllSteps_custom() -- bad args
  rc = fc_trackAllSteps_custom(-1, numSegments, segments_p, group,
			       fc_computeNumOverlapScore,
			       fc_computeSimpleCorrespond,
			       fc_addROIsUsingBasicEvents);
  fail_unless(rc != FC_SUCCESS, "should fail if bad step");
  rc = fc_trackAllSteps_custom(numStep, NULL, segments_p, group,
			       fc_computeNumOverlapScore,
			       fc_computeSimpleCorrespond,
			       fc_addROIsUsingBasicEvents);
  fail_unless(rc != FC_SUCCESS, "should fail if NULL numSegments");
  rc = fc_trackAllSteps_custom(numStep, numSegments, NULL, group,
			       fc_computeNumOverlapScore,
			       fc_computeSimpleCorrespond,
			       fc_addROIsUsingBasicEvents);
  fail_unless(rc != FC_SUCCESS, "should fail if NULL segments");
  rc = fc_trackAllSteps_custom(numStep, numSegments, segments_p, NULL,
			       fc_computeNumOverlapScore,
			       fc_computeSimpleCorrespond,
			       fc_addROIsUsingBasicEvents);
  fail_unless(rc != FC_SUCCESS, "should fail if NULL group");

  // fc_trackStep() -- bad args
  rc = fc_trackStep(-1, numSegments[0], segments[0], group);
  fail_unless(rc != FC_SUCCESS, "should fail if bad step");
  rc = fc_trackStep(1, numSegments[0], segments[0], group);
  fail_unless(rc != FC_SUCCESS, "should fail if not next step ID");
  rc = fc_trackStep(0, -1, segments[0], group);
  fail_unless(rc != FC_SUCCESS, "should fail if bad numSegment");
  rc = fc_trackStep(0, 1, NULL, group);
  fail_unless(rc != FC_SUCCESS, 
	      "should fail if numSegment > 0 & null segments");
  rc = fc_trackStep(0, numSegments[0], segments[0], NULL);
  fail_unless(rc != FC_SUCCESS, "should fail if bad group");

  // fc_trackAllSteps() -- bad args
  rc = fc_trackAllSteps(-1, numSegments, segments_p, group);
  fail_unless(rc != FC_SUCCESS, "should fail if bad step");
  rc = fc_trackAllSteps(numStep, NULL, segments_p, group);
  fail_unless(rc != FC_SUCCESS, "should fail if NULL numSegments");
  rc = fc_trackAllSteps(numStep, numSegments, NULL, group);
  fail_unless(rc != FC_SUCCESS, "should fail if NULL segments");
  rc = fc_trackAllSteps(numStep, numSegments, segments_p, NULL);
  fail_unless(rc != FC_SUCCESS, "should fail if NULL group");

  // cleanup
  fc_freeFeatureGroup(group);
  rc = fc_deleteDataset(dataset);
  fail_unless(rc == FC_SUCCESS, "deleting dataset at end of test failed");
}
END_TEST


START_TEST(self_corresponds)
{
  FC_ReturnCode rc;
  int i, j; 
  int temp_numGroup, *temp_numPerGroup, **temp_idsPerGroup;

  // Case 1: identity matrix - no grouping
  int num1 = 5;
  int A1[5][5] = { { 1, 0, 0, 0, 0 },
		   { 0, 1, 0, 0, 0 },
		   { 0, 0, 1, 0, 0 },
		   { 0, 0, 0, 1, 0 },
		   { 0, 0, 0, 0, 1 } };
  int numGroup1 = 5;
  int numPerGroup1[5] = { 1, 1, 1, 1, 1 };
  int idsPerGroup1[5][1] = { { 0 },
			     { 1 }, 
			     { 2 },
			     { 3 },
			     { 4 } };
  int** temp_A;
  // Case 2: all correspond - 1 big group
  int num2 = 5;
  int A2[5][5] = { { 1, 1, 1, 1, 1 },
		   { 0, 1, 1, 1, 1 },
		   { 0, 0, 1, 1, 1 },
		   { 0, 0, 0, 1, 1 },
		   { 0, 0, 0, 0, 1 } };
  int numGroup2 = 1;
  int numPerGroup2[1] = { 5 };
  int idsPerGroup2[1][5] = { { 0, 1, 2, 3, 4 } };

  // Case 3: some grouping at first & last
  int num3 = 5;
  int A3[5][5] = { { 1, 1, 0, 0, 0 },
		   { 0, 1, 0, 0, 0 },
		   { 0, 0, 1, 0, 0 },
		   { 0, 0, 0, 1, 1 },
		   { 0, 0, 0, 0, 1 } };
  int numGroup3 = 3;
  int numPerGroup3[3] = { 2, 1, 2 };
  int idsPerGroup3[3][2] = { { 0, 1 }, 
			     { 2 },
			     { 3, 4 } };

  // Case 4: some grouping in middle
  int num4 = 5;
  int A4[5][5] = { { 1, 0, 0, 0, 0 },
		   { 0, 1, 1, 0, 0 },
		   { 0, 0, 1, 1, 0 },
		   { 0, 0, 0, 1, 0 },
		   { 0, 0, 0, 0, 1 } };
  int numGroup4 = 3;
  int numPerGroup4[3] = { 1, 3, 1 };
  int idsPerGroup4[3][3] = { { 0 }, 
			     { 1, 2, 3 },
			     { 4 } };

  // Case 5: this case added because it was causing problems with early
  // version of the algorithm. A parent is not from diagonal OR first 
  // non-diagagonal entry of row (see 2nd row)
  int num5 = 5;
  int A5[5][5] = { { 1, 0, 0, 0, 1 },
		   { 0, 1, 0, 1, 1 },
		   { 0, 0, 1, 0, 0 },
		   { 0, 0, 0, 1, 0 },
		   { 0, 0, 0, 0, 1 } };
  int numGroup5 = 2;
  int numPerGroup5[2] = { 4, 1 };
  int idsPerGroup5[2][4] = { { 0, 1, 3, 4 }, { 2 } };

  // Case 6: this case added because it was causing problems with early
  // version of the algorithm. The algorithm was assigning obj 4 to
  // two groups, when it should have combined groups.
  int num6 = 5;
  int A6[5][5] = { { 1, 0, 1, 0, 0 },
		   { 0, 1, 0, 1, 0 },
		   { 0, 0, 1, 0, 1 },
		   { 0, 0, 0, 1, 1 },
		   { 0, 0, 0, 0, 1 } };
  int numGroup6 = 1;
  int numPerGroup6[1] = { 5 };
  int idsPerGroup6[1][5] = { { 0, 1, 2, 3, 4 } };

  // Case 7: make sure that single entry works
  int num7 = 1;
  int A7[1][1] = { { 1 } };
  int numGroup7 = 1;
  int numPerGroup7[1] = { 1 };
  int idsPerGroup7[1][1] = { { 0 } };

  // Note: stupid C! int[5][5] isn't the same as int**
  // Can't easily loop all of these tests, arggh!!!

  // Make temp_A with space for largest A
  temp_A = malloc(num1*sizeof(int**));
  for (i = 0; i < num1; i++) 
    temp_A[i] = malloc(num1*sizeof(int*));

  // Case 1 -
  for (i = 0; i < num1; i++)
    for (j = 0; j < num1; j++)
      temp_A[i][j] = A1[i][j];
  rc = fc_groupSelfCorrespond(num1, temp_A, &temp_numGroup,
			      &temp_numPerGroup, &temp_idsPerGroup);
  fail_unless(rc == FC_SUCCESS, "failed to group");
  fail_unless(temp_numGroup == numGroup1, "mismatch of numGroup1");
  for (i = 0; i < numGroup1; i++) {
    fail_unless(temp_numPerGroup[i] = numPerGroup1[i], 
		"mismatch of numPerGroup1");
    for (j = 0; j < numPerGroup1[i]; j++)
      fail_unless(temp_idsPerGroup[i][j] == idsPerGroup1[i][j],
		  "mismatch of idsPerGroup1");
  }
  free(temp_numPerGroup);
  for (i = 0; i < temp_numGroup; i++)
    free(temp_idsPerGroup[i]);
  free(temp_idsPerGroup);

  // Case 2 -
  for (i = 0; i < num2; i++)
    for (j = 0; j < num2; j++)
      temp_A[i][j] = A2[i][j];
  rc = fc_groupSelfCorrespond(num2, temp_A, &temp_numGroup,
			      &temp_numPerGroup, &temp_idsPerGroup);
  fail_unless(rc == FC_SUCCESS, "failed to group");
  fail_unless(temp_numGroup == numGroup2, "mismatch of numGroup2");
  for (i = 0; i < numGroup2; i++) {
    fail_unless(temp_numPerGroup[i] = numPerGroup2[i], 
		"mismatch of numPerGroup2");
    for (j = 0; j < numPerGroup2[i]; j++)
      fail_unless(temp_idsPerGroup[i][j] == idsPerGroup2[i][j],
		  "mismatch of idsPerGroup2");
  }
  free(temp_numPerGroup);
  for (i = 0; i < temp_numGroup; i++)
    free(temp_idsPerGroup[i]);
  free(temp_idsPerGroup);

  // Case 3 -
  for (i = 0; i < num3; i++)
    for (j = 0; j < num3; j++)
      temp_A[i][j] = A3[i][j];
  rc = fc_groupSelfCorrespond(num3, temp_A, &temp_numGroup,
			      &temp_numPerGroup, &temp_idsPerGroup);
  fail_unless(rc == FC_SUCCESS, "failed to group");
  fail_unless(temp_numGroup == numGroup3, "mismatch of numGroup3");
  for (i = 0; i < numGroup3; i++) {
    fail_unless(temp_numPerGroup[i] = numPerGroup3[i], 
		"mismatch of numPerGroup3");
    for (j = 0; j < numPerGroup3[i]; j++)
      fail_unless(temp_idsPerGroup[i][j] == idsPerGroup3[i][j],
		  "mismatch of idsPerGroup3");
  }
  free(temp_numPerGroup);
  for (i = 0; i < temp_numGroup; i++)
    free(temp_idsPerGroup[i]);
  free(temp_idsPerGroup);

  // Case 4 -
  for (i = 0; i < num4; i++)
    for (j = 0; j < num4; j++)
      temp_A[i][j] = A4[i][j];
  rc = fc_groupSelfCorrespond(num4, temp_A, &temp_numGroup,
			      &temp_numPerGroup, &temp_idsPerGroup);
  fail_unless(rc == FC_SUCCESS, "failed to group");
  fail_unless(temp_numGroup == numGroup4, "mismatch of numGroup4");
  for (i = 0; i < numGroup4; i++) {
    fail_unless(temp_numPerGroup[i] = numPerGroup4[i], 
		"mismatch of numPerGroup4");
    for (j = 0; j < numPerGroup4[i]; j++)
      fail_unless(temp_idsPerGroup[i][j] == idsPerGroup4[i][j],
		  "mismatch of idsPerGroup4");
  }
  free(temp_numPerGroup);
  for (i = 0; i < temp_numGroup; i++)
    free(temp_idsPerGroup[i]);
  free(temp_idsPerGroup);

  // Case 5 -
  for (i = 0; i < num5; i++)
    for (j = 0; j < num5; j++)
      temp_A[i][j] = A5[i][j];
  rc = fc_groupSelfCorrespond(num5, temp_A, &temp_numGroup,
			      &temp_numPerGroup, &temp_idsPerGroup);
  fail_unless(rc == FC_SUCCESS, "failed to group");
  fail_unless(temp_numGroup == numGroup5, "mismatch of numGroup5");
  for (i = 0; i < numGroup5; i++) {
    fail_unless(temp_numPerGroup[i] = numPerGroup5[i], 
		"mismatch of numPerGroup5");
    for (j = 0; j < numPerGroup5[i]; j++)
      fail_unless(temp_idsPerGroup[i][j] == idsPerGroup5[i][j],
		  "mismatch of idsPerGroup5");
  }
  free(temp_numPerGroup);
  for (i = 0; i < temp_numGroup; i++)
    free(temp_idsPerGroup[i]);
  free(temp_idsPerGroup);

  // Case 6 -
  for (i = 0; i < num6; i++)
    for (j = 0; j < num6; j++)
      temp_A[i][j] = A6[i][j];
  rc = fc_groupSelfCorrespond(num6, temp_A, &temp_numGroup,
			      &temp_numPerGroup, &temp_idsPerGroup);
  fail_unless(rc == FC_SUCCESS, "failed to group");
  fail_unless(temp_numGroup == numGroup6, "mismatch of numGroup6");
  for (i = 0; i < numGroup6; i++) {
    fail_unless(temp_numPerGroup[i] = numPerGroup6[i], 
		"mismatch of numPerGroup6");
    for (j = 0; j < numPerGroup6[i]; j++)
      fail_unless(temp_idsPerGroup[i][j] == idsPerGroup6[i][j],
		  "mismatch of idsPerGroup6");
  }
  free(temp_numPerGroup);
  for (i = 0; i < temp_numGroup; i++)
    free(temp_idsPerGroup[i]);
  free(temp_idsPerGroup);

  // Case 7 -
  for (i = 0; i < num7; i++)
    for (j = 0; j < num7; j++)
      temp_A[i][j] = A7[i][j];
  rc = fc_groupSelfCorrespond(num7, temp_A, &temp_numGroup,
			      &temp_numPerGroup, &temp_idsPerGroup);
  fail_unless(rc == FC_SUCCESS, "failed to group");
  fail_unless(temp_numGroup == numGroup7, "mismatch of numGroup7");
  for (i = 0; i < numGroup7; i++) {
    fail_unless(temp_numPerGroup[i] = numPerGroup7[i], 
		"mismatch of numPerGroup7");
    for (j = 0; j < numPerGroup7[i]; j++)
      fail_unless(temp_idsPerGroup[i][j] == idsPerGroup7[i][j],
		  "mismatch of idsPerGroup7");
  }
  free(temp_numPerGroup);
  for (i = 0; i < temp_numGroup; i++)
    free(temp_idsPerGroup[i]);
  free(temp_idsPerGroup);

  // cleanup
  for (i = 0; i < num1; i++)
    free(temp_A[i]);
  free(temp_A);
}
END_TEST

// Populate the Suite with the tests

Suite *track_suite(void)
{
  Suite *suite = suite_create("FeatureTracking");

  TCase *tc_default = tcase_create(" - Default Routines ");
  TCase *tc_track = tcase_create(" - Track Interface ");
  TCase *tc_related = tcase_create(" - Related Routines ");

  suite_add_tcase(suite, tc_default);
  tcase_add_checked_fixture(tc_default, track_setup, track_teardown);
  tcase_add_test(tc_default, default_scores);
  tcase_add_test(tc_default, default_corresponds);
  tcase_add_test(tc_default, default_addROIs);

  suite_add_tcase(suite, tc_track);
  tcase_add_checked_fixture(tc_track, track_setup, track_teardown);
  tcase_add_test(tc_track, track_tests);

  suite_add_tcase(suite, tc_related);
  tcase_add_checked_fixture(tc_related, track_setup, track_teardown);
  tcase_add_test(tc_related, self_corresponds);

  return suite;
}
