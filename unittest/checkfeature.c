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
 * \file checkfeature.c
 * \brief Unit testing of \ref Features module.
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/unittest/checkfeature.c,v $
 * $Revision: 1.30 $ 
 * $Date: 2006/08/30 19:20:05 $
 */
#include <string.h>
#include <stdlib.h>
#include <check.h>
#include "fc.h"
#include "fcP.h"
#include "checkall.h"

START_TEST(feature_basic)
{
  FC_ReturnCode rc;
  int i;
  int numROI = 3, numParent = 2, numChild = 4;
  int id = 77;
  int roiIDs[3] = { 12, 23, 41};
  int stepIDs[3] = { 2, 3, 4};
  int parentIDs[2] = { 74, 76 };
  int childIDs[4] = { 78, 79, 80, 81 };
  _FC_Feature* feature, *feature2;


  // --- test values of new feature

  rc = _fc_createFeature(id, &feature);
  fail_unless(rc == FC_SUCCESS, "create failed");
  fail_unless(feature->id == id, "mismatch of id");
  fail_unless(feature->numROI == 0 && feature->stepIDs == NULL && 
	      feature->roiIDs == NULL, "new feature should have no rois");
  fail_unless(feature->numParent == 0 && feature->parentIDs == NULL,
	      "new feature should have no parents");
  fail_unless(feature->numChild == 0 && feature->childIDs == NULL,
	      "new featuer should have no children");
  
  // --- test adding rois, parents & children

  // ROI IDs
  for (i = 0; i < numROI; i++) {
    rc = _fc_feature_addROI(feature, stepIDs[i], roiIDs[i]);
    fail_unless(rc == FC_SUCCESS, "addROI failed");
  }
  fail_unless(feature->numROI == numROI, "mismatch of numROI");
  for (i = 0; i < numROI; i++) {
    fail_unless(feature->stepIDs[i] == stepIDs[i], "mismatch in stepIDs");
    fail_unless(feature->roiIDs[i] == roiIDs[i], "mismatch in ROI_IDs");
  }

  // parent IDs
  for (i = 0; i < numParent; i++) {
    rc = _fc_feature_addParent(feature, parentIDs[i]);
    fail_unless(rc == FC_SUCCESS, "addParent failed");
  }
  fail_unless(feature->numParent == numParent, "mismatch of numParent");
  for (i = 0; i < numParent; i++)
    fail_unless(feature->parentIDs[i] == parentIDs[i], "mismatch in parentIDs");
 
  // child IDs
  for (i = 0; i < numChild; i++) {
    rc = _fc_feature_addChild(feature, childIDs[i]);
    fail_unless(rc == FC_SUCCESS, "addChild failed");
  }
  fail_unless(feature->numChild == numChild, "mismatch of numChild");
  for (i = 0; i < numChild; i++)
    fail_unless(feature->childIDs[i] == childIDs[i], "mismatch in childIDs");

  // --- test bad function calls

  // _fc_createFeature() -- bad args
  feature2 = (_FC_Feature*)1;
  rc = _fc_createFeature(-1, &feature2);
  fail_unless(rc != FC_SUCCESS, "should fail if id < 0");
  fail_unless(feature2 == NULL, "fail should return null");
  rc = _fc_createFeature(0, NULL);
  fail_unless(rc != FC_SUCCESS, "should fail if NULL feature pointer");

  // _fc_feature_addROI() -- bad args
  rc = _fc_feature_addROI(feature, -1, 0);
  fail_unless(rc != FC_SUCCESS, "should fail if id < 0");
  rc = _fc_feature_addROI(feature, 0, -1);
  fail_unless(rc != FC_SUCCESS, "should fail if id < 0");

  // _fc_feature_addParent() -- bad args
  rc = _fc_feature_addParent(feature, -1);
  fail_unless(rc != FC_SUCCESS, "should fail if id < 0");

  // _fc_feature_addChild() -- bad args
  rc = _fc_feature_addChild(feature, -1);
  fail_unless(rc != FC_SUCCESS, "should fail if id < 0");

  // Make sure feature hasn't changed
  fail_unless(feature->id == id, "mismatch of id");
  fail_unless(feature->numROI == numROI, "mismatch of numROI");
  for (i = 0; i < numROI; i++) {
    fail_unless(feature->stepIDs[i] == stepIDs[i], "mismatch in stepIDs");
    fail_unless(feature->roiIDs[i] == roiIDs[i], "mismatch in ROI_IDs");
  }
  fail_unless(feature->numParent == numParent, "mismatch of numParent");
  for (i = 0; i < numParent; i++)
    fail_unless(feature->parentIDs[i] == parentIDs[i], "mismatch in parentIDs");
  fail_unless(feature->numChild == numChild, "mismatch of numChild");
  for (i = 0; i < numChild; i++)
    fail_unless(feature->childIDs[i] == childIDs[i], "mismatch in childIDs");

  // --- all done - cleanup

  _fc_freeFeature(feature);
}
END_TEST

START_TEST(featureGroup_basic)
{
  FC_ReturnCode rc;
  int i, j;
  FC_FeatureGroup* group;
  _FC_Feature* feature1 = NULL, *feature2 = NULL;
  int stepIDs[3] = { 1, 4, 1 };
  int temp_numStep, temp_numROI, *temp_stepIDs;
  int temp_numFeature, *temp_featureIDs, temp_id;
  int temp_numParent, *temp_parentIDs, temp_numChild, *temp_childIDs;
  FC_Subset *temp_rois, temp_roi;
  _FC_Feature** temp_features;
  // Creating this structure (adding ROI's in order)
  //
  //      stepID:  0    1    2    3    4    5    6
  //   -------------------------------------------
  //   featureID:  0 .. 0 ->   -> 3
  //                        \ /
  //                         2  
  //                        / \                             .
  //                    1 ->   -> 4 .. 4         5
  //
  int numStep = 7, numROI = 8, numFeature = 6;
  FC_Subset rois[8] = { { 0, 0 }, { 1, 1 }, { 2, 2 }, { 3, 3 }, { 4, 4 },
			{ 5, 5 }, { 6, 6 }, { 7, 7 } };
  int numROIperStep[7] = { 1, 2, 1, 2, 1, 0, 1 };
  int featureIDsPerStep[7][2] = { { 0, -1 }, { 0,  1 }, { 2, -1 }, 
                           { 3,  4 }, { 4, -1 }, { -1, -1 }, { 5, -1 } };
  int numROIPerFeature[6] = { 2, 1, 1, 1, 2, 1 };
  int stepIDsPerFeature[6][3] = { { 0, 1 }, { 1 }, { 2 }, { 3 }, { 3, 4 }, 
				  { 6 } };
  FC_Subset roisPerFeature[6][3] = { { rois[0], rois[1] }, { rois[2] }, 
				     { rois[3] }, { rois[4] }, 
				     { rois[5], rois[6] }, { rois[7] } };
  int numParentPerFeature[6] = { 0, 0, 2, 1, 1, 0 };
  int parentIDsPerFeature[6][2] = { { }, { }, { 0, 1 }, { 2 }, { 2 }, { } };
  int numChildPerFeature[6] = { 1, 1, 2, 0, 0, 0 };
  int childIDsPerFeature[6][2] = { { 2 }, { 2 }, { 3, 4 }, { }, { }, { } };

  // --- queries on new (empty) feature group

  // create
  rc = fc_createFeatureGroup(&group);
  fail_unless(rc == FC_SUCCESS, "create failed");

  // Group query
  rc = fc_featureGroup_getNumStep(group, &temp_numStep);
  fail_unless(rc == FC_SUCCESS, "getNumStep failed on new group");
  fail_unless(temp_numStep == 0, "getNumStep on new group should return 0");
  rc = fc_featureGroup_getNumROI(group, &temp_numROI);
  fail_unless(rc == FC_SUCCESS, "getNumROI failed on new group");
  fail_unless(temp_numROI == 0,  "getNumROI on new group should return 0");
  rc = fc_featureGroup_getNumFeature(group, &temp_numFeature);
  fail_unless(rc == FC_SUCCESS, "getNumFeature on new group failed");
  fail_unless(temp_numFeature == 0, "should be 0 features in new group");

  // private group query
  rc = _fc_featureGroup_getFeatures(group, &temp_numFeature, &temp_features);
  fail_unless(rc == FC_SUCCESS, "getFeatures on new group failed");
  fail_unless(temp_numFeature == 0, 
	      "getFeatures on new group should return numFeature=0");
  fail_unless(temp_features == NULL, 
	      "getFeatures on new group should return NULL");

  // --- test private Helper routines

  // ROI -- use helper to add ROIs and check
  for (i = 0; i < 3; i++) {
    rc = _fc_featureGroup_addROI(group, i+1, stepIDs[i], rois[i], &temp_id);
    fail_unless(rc == FC_SUCCESS, "addROI failed");
  }
  fail_unless(group->numStep == 5, "unexpected numStep");
  fail_unless(group->numROIperStep[0] == 0, "unexpected numROIperStep[]");
  fail_unless(group->numROIperStep[1] == 2, "unexpected numROIperStep[]");
  fail_unless(group->numROIperStep[2] == 0, "unexpected numROIperStep[]");
  fail_unless(group->numROIperStep[3] == 0, "unexpected numROIperStep[]");
  fail_unless(group->numROIperStep[4] == 1, "unexpected numROIperStep[]");
  fail_unless(group->rois[0] == NULL, "unexpected rois[]");
  fail_unless(FC_HANDLE_EQUIV(group->rois[1][0], rois[0]), "unexpected rois[][]");
  fail_unless(FC_HANDLE_EQUIV(group->rois[1][1], rois[2]), "unexpected rois[][]");
  fail_unless(group->rois[2] == NULL, "unexpected rois[]");
  fail_unless(group->rois[3] == NULL, "unexpected rois[]");
  fail_unless(FC_HANDLE_EQUIV(group->rois[4][0], rois[1]), "unexpected rois[][]");
  fail_unless(group->roiFeatureIDs != NULL, "featureIDs should not be NULL");
  fail_unless(group->roiFeatureIDs[0] == NULL, "featureIDs[0] should be NUL");
  fail_unless(group->roiFeatureIDs[1][0] == 1, "featureIDs[1][0] should be 1");
  fail_unless(group->roiFeatureIDs[1][1] == 3, "featureIDs[1][1] should be 3");
  fail_unless(group->roiFeatureIDs[2] == NULL, "featureIDs[2] should be NUL");
  fail_unless(group->roiFeatureIDs[3] == NULL, "featureIDs[3] should be NUL");
  fail_unless(group->roiFeatureIDs[4][0] == 2, "featureIDs[4][0] should be 2");
	       
  // Feature -- use helper to add features and then check
  // (just passes around feature addresses so don't need to set them up)
  _fc_createFeature(2, &feature1);
  _fc_createFeature(12, &feature2);
  rc = _fc_featureGroup_addFeature(group, feature1);
  fail_unless(rc == FC_SUCCESS, "addFeature failed");
  rc = _fc_featureGroup_addFeature(group, feature2);
  fail_unless(rc == FC_SUCCESS, "addFeature failed");
  fail_unless(group->numFeature == 2, "unexpected numFeature");
  fail_unless(group->features[0] == feature1, "unexpected features[]");
  fail_unless(group->features[1] == feature2, "unexpected features[]");
  fail_unless(feature1->id == 0, "unexpected feature id");
  fail_unless(feature2->id == 1, "unexpected feature id");

  // close
  fc_freeFeatureGroup(group);

  //----------------------------------------------------------
  // test manipulations interface
  //----------------------------------------------------------
  // Creating this structure (adding ROI's in order)
  //
  //      stepID:  0    1    2    3    4    5    6
  //   -------------------------------------------
  //   featureID:  0 .. 0 ->   -> 3
  //                        \ /
  //                         2  
  //                        / \                             .
  //                    1 ->   -> 4 .. 4         5
  //
  rc = fc_createFeatureGroup(&group);
  fail_unless(rc == FC_SUCCESS, "failed to create feature group");
  rc = fc_featureGroup_newFeature(group, 0, rois[0], &temp_id);
  fail_unless(rc == FC_SUCCESS, "newFeature failed");
  fail_unless(temp_id == 0, "new feature should have id = 0");
  fail_unless(group->numStep == 1 && group->numROIperStep[0] == 1 && 
              group->roiFeatureIDs[0][0] == 0,
              "newFeature didn't do ROIs as expected");
  fail_unless(group->numFeature == 1 &&
              group->features[0]->id == 0 &&  
              group->features[0]->numROI == 1 &&
              group->features[0]->stepIDs[0] == 0 &&
              group->features[0]->roiIDs[0] == 0, 
              "newFeature didn't do features as expected");
  rc = fc_featureGroup_continueFeature(group, temp_id, 1, rois[1]);
  fail_unless(rc == FC_SUCCESS, "continueFeauture failed");
  fail_unless(group->numStep == 2 && 
              group->numROIperStep[0] == 1 && group->numROIperStep[1] == 1 && 
              group->roiFeatureIDs[0][0] == 0 &&
              group->roiFeatureIDs[1][0] == 0,
              "continueFeature didn't do ROIs as expected");
  fail_unless(group->numFeature == 1 && 
              group->features[0]->id == 0 && 
              group->features[0]->numROI == 2 &&
              group->features[0]->stepIDs[0] == 0 &&
              group->features[0]->roiIDs[0] == 0 && 
              group->features[0]->stepIDs[1] == 1 &&
              group->features[0]->roiIDs[0] == 0, 
              "continueFeature didn't do features as expected");
  // make a few more features to make interesting
  fc_featureGroup_newFeature(group, 1, rois[2], &temp_id); // id = 1
  fc_featureGroup_newFeature(group, 2, rois[3], &temp_id); // id = 2
  fc_featureGroup_newFeature(group, 3, rois[4], &temp_id); // id = 3
  fc_featureGroup_newFeature(group, 3, rois[5], &temp_id); // id = 4
  fc_featureGroup_continueFeature(group, temp_id, 4, rois[6]);
  fc_featureGroup_noFeaturesThisStep(group, 5);
  fail_unless(group->numStep == 6 && group->numROIperStep[5] == 0 &&
	      group->rois[5] == NULL && group->roiFeatureIDs[5] == 0,
	      "noFeaturesThisStep didn't do rois as expected");
  fc_featureGroup_newFeature(group, 6, rois[7], &temp_id); // id = 5 
  // 
  rc = fc_featureGroup_setParentChildPair(group, 0, 2);
  fail_unless(rc == FC_SUCCESS, "addChild failed");
  fc_featureGroup_setParentChildPair(group, 2, 3);
  fc_featureGroup_setParentChildPair(group, 1, 2);
  fc_featureGroup_setParentChildPair(group, 2, 4);
  fail_unless(group->features[0]->numChild == 1 &&
              group->features[0]->childIDs[0] == 2 &&
              group->features[0]->numParent == 0 &&
              group->features[1]->numChild == 1 &&
              group->features[1]->childIDs[0] == 2 &&
              group->features[1]->numParent == 0 &&
              group->features[2]->numChild == 2 &&
              group->features[2]->childIDs[0] == 3 &&
              group->features[2]->childIDs[1] == 4 &&
              group->features[2]->numParent == 2 &&
              group->features[2]->parentIDs[0] == 0 &&
              group->features[2]->parentIDs[1] == 1 &&
              group->features[3]->numChild == 0 &&
              group->features[3]->numParent == 1 &&
              group->features[3]->parentIDs[0] == 2,
              "addChilds didn't do as expected");

  //----------------------------------------------------------
  // test query interface
  //----------------------------------------------------------

  // Group query
  rc = fc_featureGroup_getNumStep(group, &temp_numStep);
  fail_unless(rc == FC_SUCCESS, "getNumStep failed");
  fail_unless(temp_numStep == numStep, "numStep mismatch");
  rc = fc_featureGroup_getNumFeature(group, &temp_numFeature);
  fail_unless(rc == FC_SUCCESS, "getNumFeature failed");
  fail_unless(temp_numFeature == numFeature, "should have numFeature features");
  rc = fc_featureGroup_getNumROI(group, &temp_numROI);
  fail_unless(rc == FC_SUCCESS, "getNumROI failed");
  fail_unless(temp_numROI == numROI, "numStep mismatch");

  // Query step of a group
  for (i = 0; i < numStep; i++) {
    rc = fc_featureGroup_getNumFeatureInStep(group, i, &temp_numFeature);
    fail_unless(rc == FC_SUCCESS, "getNumFeatureInStep failed");
    fail_unless(temp_numFeature == numROIperStep[i], 
		"mismatch in numFeaureInStep");
    rc = fc_featureGroup_getFeatureIDsInStep(group, i, &temp_numFeature, 
					     &temp_featureIDs);
    fail_unless(rc == FC_SUCCESS, "getFeatureIDsInStep failed");
    fail_unless(temp_numFeature == numROIperStep[i], 
		"getFeatureIDsInStep should return numROI[i]");
    if (temp_numFeature > 0) {
      for (j = 0; j < numROIperStep[i]; j++)
	fail_unless(temp_featureIDs[j] == featureIDsPerStep[i][j], 
		    "featureIDs mismatch");
      free(temp_featureIDs);
    }
    else 
      fail_unless(temp_featureIDs == NULL, "featureIDs mismatch");
    rc = fc_featureGroup_getROIsInStep(group, i, &temp_numROI, &temp_rois);
    fail_unless(rc == FC_SUCCESS, "getROIsInStep failed");
    fail_unless(temp_numROI == numROIperStep[i], "numStep mismatch");
    if (temp_numROI > 0) {
      for (j = 0; j < numROIperStep[i]; j++)
	fail_unless(FC_HANDLE_EQUIV(temp_rois[j], group->rois[i][j]), 
		    "rois mismatch");
      free(temp_rois);
    }
    else 
      fail_unless(temp_rois == NULL, "rois mismatch");
  }
  // Query features
  for (i = 0; i < numFeature; i++) {
    rc = fc_getFeatureNumROI(group, i, &temp_numROI);
    fail_unless(rc == FC_SUCCESS, "failed to get numROI in feature");
    fail_unless(temp_numROI == numROIPerFeature[i], "mismatch of numROI");
    rc = fc_getFeatureROIs(group, i, &temp_numROI, &temp_stepIDs, &temp_rois);
    fail_unless(rc == FC_SUCCESS, "failed to get rois in feature");
    fail_unless(temp_numROI == numROIPerFeature[i], "wrong numRoi in feature");
    for (j = 0; j < temp_numROI; j++) {
      fail_unless(temp_stepIDs[j] == stepIDsPerFeature[i][j],
		  "mismtach of step IDs");
      fail_unless(FC_HANDLE_EQUIV(temp_rois[j], roisPerFeature[i][j]),
		  "mismatch of roi");
      rc = fc_getFeatureROIAtStep(group, i, temp_stepIDs[j], &temp_roi);
      fail_unless(rc == FC_SUCCESS, "failed to get ROI at step");
      fail_unless(FC_HANDLE_EQUIV(temp_roi, temp_rois[j]),
		  "mismatch of roi");
    }
    free(temp_stepIDs);
    free(temp_rois);
    rc = fc_getFeatureNumParent(group, i, &temp_numParent);
    fail_unless(rc == FC_SUCCESS, "failed to get numParent in feature");
    fail_unless(temp_numParent == numParentPerFeature[i],
		"mismatch of numParent");
    rc = fc_getFeatureParentIDs(group, i, &temp_numParent, &temp_parentIDs);
    fail_unless(rc == FC_SUCCESS, "failed to get numParent in feature");
    fail_unless(temp_numParent == numParentPerFeature[i],
		"mismatch of numParent");
    if (temp_numParent > 0) {
      for (j = 0; j < temp_numParent; j++)
	fail_unless(temp_parentIDs[j] == parentIDsPerFeature[i][j],
		    "mismatch of parent IDs");
      free(temp_parentIDs);
    }
    else
      fail_unless(temp_parentIDs == NULL, "mismatch of parent IDs");
    rc = fc_getFeatureNumChild(group, i, &temp_numChild);
    fail_unless(rc == FC_SUCCESS, "failed to get numChild in feature");
    fail_unless(temp_numChild == numChildPerFeature[i],
		"mismatch of numChild");
    rc = fc_getFeatureChildIDs(group, i, &temp_numChild, &temp_childIDs);
    fail_unless(rc == FC_SUCCESS, "failed to get numChild in feature");
    fail_unless(temp_numChild == numChildPerFeature[i],
		"mismatch of numChild");
    if (temp_numChild > 0) {
      for (j = 0; j < temp_numChild; j++)
	fail_unless(temp_childIDs[j] == childIDsPerFeature[i][j],
		    "mismatch of child IDs");
      free(temp_childIDs);
    }
    else
      fail_unless(temp_childIDs == NULL, "mismatch of child IDs");
  }

  //----------------------------------------------------------
  // test private feature query interface
  //----------------------------------------------------------
  for (i = 0; i < numFeature; i++) {
    rc = _fc_featureGroup_getFeature(group, i, &feature1);
    fail_unless(group->features[i] == feature1, "feature mismatch");
  }
  rc = _fc_featureGroup_getFeatures(group, &temp_numFeature, &temp_features);
  fail_unless(rc == FC_SUCCESS, "getFeatures failed");
  fail_unless(temp_numFeature == numFeature, 
	      "getFeatures returned wrong numFeatuer");
  fail_unless(temp_features != NULL, "getFeatures failed to return features");
  for (i = 0; i < numFeature; i++) 
    fail_unless(temp_features[i] == group->features[i], "features mismatch");
  free(temp_features);

  //----------------------------------------------------------
  // test special cases
  //----------------------------------------------------------

  // fc_freeFeatureGroup() -- o.k. to free NULL
  fc_freeFeatureGroup(NULL); // shouldn't core dump

  // fc_getFeatureROIs() -- o.k. to pass NULL for stepIDs OR ROIs
  rc = fc_getFeatureROIs(group, 0, &temp_numROI, NULL, &temp_rois);
  fail_unless(rc == FC_SUCCESS, "should not fail if NULL stepIDs");
  fail_unless(temp_numROI == numROIPerFeature[0], "mismatch of numROI");
  for (i = 0; i < temp_numROI; i++)
    fail_unless(FC_HANDLE_EQUIV(temp_rois[i], roisPerFeature[0][i]),
		"mismatch of rois");
  free(temp_rois);
  rc = fc_getFeatureROIs(group, 0, &temp_numROI, &temp_stepIDs, NULL);
  fail_unless(rc == FC_SUCCESS, "should not fail if NULL ROIs");
  fail_unless(temp_numROI == numROIPerFeature[0], "mismatch of numROI");
  for (i = 0; i < temp_numROI; i++)
    fail_unless(temp_stepIDs[i] == stepIDsPerFeature[0][i],
      "mismatch of rois");
  free(temp_stepIDs);

  //----------------------------------------------------------
  // test error conditions
  //----------------------------------------------------------

  // fc_getFeatureROIAtStep() -- should fail if feature doesn't exist at 
  //that step
  rc = fc_getFeatureROIAtStep(group, 1, 0, &temp_roi);
  fail_unless(rc != FC_SUCCESS, "should fail if feature doensn't exist");
  fail_unless(FC_HANDLE_EQUIV(temp_roi, FC_NULL_SUBSET),
	      "fail should return NULL");

  // fc_createFeatureGroup -- bad args
  rc = fc_createFeatureGroup(NULL);
  fail_unless(rc != FC_SUCCESS, "should fail to create if NULL");

  // fc_featureGroup_newFeature() -- bad args
  rc = fc_featureGroup_newFeature(NULL, numStep, rois[0], &temp_id);
  fail_unless(rc != FC_SUCCESS, "should fail if bad group");
  fail_unless(temp_id == -1, "fail should return -1 id");
  rc = fc_featureGroup_newFeature(group, -1, rois[0], &temp_id);
  fail_unless(rc != FC_SUCCESS, "should fail if bad step id");
  fail_unless(temp_id == -1, "fail should return -1 id");
  rc = fc_featureGroup_newFeature(group, numStep, rois[0], NULL);
  fail_unless(rc != FC_SUCCESS, "should fail if NULL id to return");
  fail_unless(temp_id == -1, "fail should return -1 id");

  // fc_featureGroup_continueFeature() -- bad args
  rc = fc_featureGroup_continueFeature(NULL, numFeature-1, numStep, rois[0]);
  fail_unless(rc != FC_SUCCESS, "should fail if bad group");
  rc = fc_featureGroup_continueFeature(group, numFeature, numStep, rois[0]);
  fail_unless(rc != FC_SUCCESS, "should fail if bad feature ID");
  rc = fc_featureGroup_continueFeature(group, numFeature-1, -1, rois[0]);
  fail_unless(rc != FC_SUCCESS, "should fail if bad step ID");

  // fc_featureGroup_setParentChildPair() -- bad args
  rc = fc_featureGroup_setParentChildPair(NULL, 3, 4);
  fail_unless(rc != FC_SUCCESS, "should fail if bad group");
  rc = fc_featureGroup_setParentChildPair(group, numFeature, 4);
  fail_unless(rc != FC_SUCCESS, "should fail if bad parent ID");
  rc = fc_featureGroup_setParentChildPair(group, 3, numFeature);
  fail_unless(rc != FC_SUCCESS, "should fail if bad child ID");

  // fc_featureGroup_getNumStep() -- bad args
  rc = fc_featureGroup_getNumStep(NULL, &temp_numStep);
  fail_unless(rc != FC_SUCCESS, "should fail if bad group");
  fail_unless(temp_numStep == -1, "fail should return null");
  rc = fc_featureGroup_getNumStep(group, NULL);
  fail_unless(rc != FC_SUCCESS, "should fail if NULL return arg");

  // fc_featureGroup_getNumFeature() -- bad args
  rc = fc_featureGroup_getNumFeature(NULL, &temp_numFeature);
  fail_unless(rc != FC_SUCCESS, "should fail if bad group");
  fail_unless(temp_numFeature == -1, "fail should return null");
  rc = fc_featureGroup_getNumFeature(group, NULL);
  fail_unless(rc != FC_SUCCESS, "should fail if NULL return arg");

  // fc_featureGroup_getNumROI() -- bad args
  rc = fc_featureGroup_getNumROI(NULL, &temp_numROI);
  fail_unless(rc != FC_SUCCESS, "should fail if bad group");
  fail_unless(temp_numROI == -1, "fail should return null");
  rc = fc_featureGroup_getNumROI(group, NULL);
  fail_unless(rc != FC_SUCCESS, "should fail if NULL return arg");

  // fc_featureGroup_getNumFeatureInStep() -- bad args
  rc = fc_featureGroup_getNumFeatureInStep(NULL, numStep-1, &temp_numFeature); 
  fail_unless(rc != FC_SUCCESS, "should fail if bad group");
  fail_unless(temp_numFeature == -1, "fail should return null");
  rc = fc_featureGroup_getNumFeatureInStep(group, numStep, &temp_numFeature); 
  fail_unless(rc != FC_SUCCESS, "should fail if bad step ID");
  fail_unless(temp_numFeature == -1, "fail should return null");
  rc = fc_featureGroup_getNumFeatureInStep(group, numStep-1, NULL); 
  fail_unless(rc != FC_SUCCESS, "should fail if NULL return arg");

  // fc_featureGroup_getFeatureIDsInStep() -- bad args
  rc = fc_featureGroup_getFeatureIDsInStep(NULL, numStep-1, &temp_numFeature,
					   &temp_featureIDs); 
  fail_unless(rc != FC_SUCCESS, "should fail if bad group");
  fail_unless(temp_numFeature == -1 && temp_featureIDs == NULL, 
	      "fail should return null");
  rc = fc_featureGroup_getFeatureIDsInStep(group, numStep, &temp_numFeature,
					   &temp_featureIDs); 
  fail_unless(rc != FC_SUCCESS, "should fail if bad step ID");
  fail_unless(temp_numFeature == -1 && temp_featureIDs == NULL,
	      "fail should return null");
  rc = fc_featureGroup_getFeatureIDsInStep(group, numStep-1, NULL, 
					   &temp_featureIDs); 
  fail_unless(rc != FC_SUCCESS, "should fail if NULL return arg");
  fail_unless(temp_featureIDs == NULL, "fail should return null");
  rc = fc_featureGroup_getFeatureIDsInStep(group, numStep, &temp_numFeature,
					   NULL); 
  fail_unless(rc != FC_SUCCESS, "should fail if NULL return arg");
  fail_unless(temp_numFeature == -1, "fail should return null");

  // fc_featureGroup_getROIsInStep() -- bad args
  rc = fc_featureGroup_getROIsInStep(NULL, numStep-1, &temp_numROI,
					   &temp_rois); 
  fail_unless(rc != FC_SUCCESS, "should fail if bad group");
  fail_unless(temp_numROI == -1 && temp_rois == NULL, 
	      "fail should return null");
  rc = fc_featureGroup_getROIsInStep(group, numStep, &temp_numROI,
					   &temp_rois); 
  fail_unless(rc != FC_SUCCESS, "should fail if bad step ID");
  fail_unless(temp_numROI == -1 && temp_rois == NULL,
	      "fail should return null");
  rc = fc_featureGroup_getROIsInStep(group, numStep-1, NULL, 
					   &temp_rois); 
  fail_unless(rc != FC_SUCCESS, "should fail if NULL return arg");
  fail_unless(temp_rois == NULL, "fail should return null");
  rc = fc_featureGroup_getROIsInStep(group, numStep, &temp_numROI,
					   NULL); 
  fail_unless(rc != FC_SUCCESS, "should fail if NULL return arg");
  fail_unless(temp_numROI == -1, "fail should return null");

  // fc_getFeatureNumROI() -- bad args
  rc = fc_getFeatureNumROI(NULL, numFeature-1, &temp_numROI);
  fail_unless(rc != FC_SUCCESS, "should fail if bad group");
  fail_unless(temp_numROI == -1, "fail should return NULL");
  rc = fc_getFeatureNumROI(group, numFeature, &temp_numROI);
  fail_unless(rc != FC_SUCCESS, "should fail if bad feature ID");
  fail_unless(temp_numROI == -1, "fail should return NULL");
  rc = fc_getFeatureNumROI(group, numFeature-1, NULL);
  fail_unless(rc != FC_SUCCESS, "should fail if NULL return arg");

  // fc_getFeatureROIs() -- bad args
  rc = fc_getFeatureROIs(NULL, numFeature-1, &temp_numROI, &temp_stepIDs,
			 &temp_rois);
  fail_unless(rc != FC_SUCCESS, "should fail if bad group");
  fail_unless(temp_numROI == -1 && temp_stepIDs == NULL && temp_rois == NULL,
	      "fail should return NULL");
  rc = fc_getFeatureROIs(group, numFeature, &temp_numROI, &temp_stepIDs,
			 &temp_rois);
  fail_unless(rc != FC_SUCCESS, "should fail if bad feature ID");
  fail_unless(temp_numROI == -1 && temp_stepIDs == NULL && temp_rois == NULL,
	      "fail should return NULL");
  rc = fc_getFeatureROIs(group, numFeature-1, NULL, &temp_stepIDs,
			 &temp_rois);
  fail_unless(rc != FC_SUCCESS, "should fail if NULL return arg");
  fail_unless(temp_stepIDs == NULL && temp_rois == NULL, 
	      "fail should return NULL");
  rc = fc_getFeatureROIs(group, numFeature-1, &temp_numROI, NULL, NULL);
  fail_unless(rc != FC_SUCCESS, "should fail if both stepIDs & ROIs are NULL");
  fail_unless(temp_numROI == -1, "fail should return NULL");

  // fc_getFeatureROIAtStep() -- bad args
  rc = fc_getFeatureROIAtStep(NULL, 0, 0, &temp_roi);
  fail_unless(rc != FC_SUCCESS, "should fail if bad group");
  fail_unless(FC_HANDLE_EQUIV(temp_roi, FC_NULL_SUBSET),
	      "fail should return NULL");
  rc = fc_getFeatureROIAtStep(group, -1, 0, &temp_roi);
  fail_unless(rc != FC_SUCCESS, "should fail if bad feature ID");
  fail_unless(FC_HANDLE_EQUIV(temp_roi, FC_NULL_SUBSET),
	      "fail should return NULL");
  rc = fc_getFeatureROIAtStep(group, numFeature, 0, &temp_roi);
  fail_unless(rc != FC_SUCCESS, "should fail if out of range feature ID");
  fail_unless(FC_HANDLE_EQUIV(temp_roi, FC_NULL_SUBSET),
	      "fail should return NULL");
  rc = fc_getFeatureROIAtStep(group, 0, -1, &temp_roi);
  fail_unless(rc != FC_SUCCESS, "should fail if bad step ID");
  fail_unless(FC_HANDLE_EQUIV(temp_roi, FC_NULL_SUBSET),
	      "fail should return NULL");
  rc = fc_getFeatureROIAtStep(group, 0, numStep, &temp_roi);
  fail_unless(rc != FC_SUCCESS, "should fail if out of range step ID");
  fail_unless(FC_HANDLE_EQUIV(temp_roi, FC_NULL_SUBSET),
	      "fail should return NULL");
  rc = fc_getFeatureROIAtStep(group, 0, 0, NULL);
  fail_unless(rc != FC_SUCCESS, "should fail if bad group");

  // fc_getFeatureNumParent() -- bad args
  rc = fc_getFeatureNumParent(NULL, numFeature-1, &temp_numParent);
  fail_unless(rc != FC_SUCCESS, "should fail if bad group");
  fail_unless(temp_numParent == -1, "fail should return NULL");
  rc = fc_getFeatureNumParent(group, numFeature, &temp_numParent);
  fail_unless(rc != FC_SUCCESS, "should fail if bad feature ID");
  fail_unless(temp_numParent == -1, "fail should return NULL");
  rc = fc_getFeatureNumParent(group, numFeature-1, NULL);
  fail_unless(rc != FC_SUCCESS, "should fail if NULL return arg");
 
  // fc_getFeatureParentIDs() -- bad args
  rc = fc_getFeatureParentIDs(NULL, numFeature-1, &temp_numParent, 
			    &temp_parentIDs);
  fail_unless(rc != FC_SUCCESS, "should fail if bad group");
  fail_unless(temp_numParent == -1 && temp_parentIDs == NULL,
	      "fail should return NULL");
  rc = fc_getFeatureParentIDs(group, numFeature, &temp_numParent,
			      &temp_parentIDs);
  fail_unless(rc != FC_SUCCESS, "should fail if bad feature ID");
  fail_unless(temp_numParent == -1 && temp_parentIDs == NULL,
	      "fail should return NULL");
  rc = fc_getFeatureParentIDs(group, numFeature-1, NULL, &temp_parentIDs);
  fail_unless(rc != FC_SUCCESS, "should fail if NULL return arg");
  fail_unless(temp_parentIDs == NULL, "fail should return NULL");
  rc = fc_getFeatureParentIDs(group, numFeature-1, &temp_numParent, NULL);
  fail_unless(rc != FC_SUCCESS, "should fail if NULL return arg");
  fail_unless(temp_numParent == -1, "fail should return NULL");

  // fc_getFeatureNumChild() -- bad args
  rc = fc_getFeatureNumChild(NULL, numFeature-1, &temp_numChild);
  fail_unless(rc != FC_SUCCESS, "should fail if bad group");
  fail_unless(temp_numChild == -1, "fail should return NULL");
  rc = fc_getFeatureNumChild(group, numFeature, &temp_numChild);
  fail_unless(rc != FC_SUCCESS, "should fail if bad feature ID");
  fail_unless(temp_numChild == -1, "fail should return NULL");
  rc = fc_getFeatureNumChild(group, numFeature-1, NULL);
  fail_unless(rc != FC_SUCCESS, "should fail if NULL return arg");
 
  // fc_getFeatureChildIDs() -- bad args
  rc = fc_getFeatureChildIDs(NULL, numFeature-1, &temp_numChild, 
			    &temp_childIDs);
  fail_unless(rc != FC_SUCCESS, "should fail if bad group");
  fail_unless(temp_numChild == -1 && temp_childIDs == NULL,
	      "fail should return NULL");
  rc = fc_getFeatureChildIDs(group, numFeature, &temp_numChild,
			      &temp_childIDs);
  fail_unless(rc != FC_SUCCESS, "should fail if bad feature ID");
  fail_unless(temp_numChild == -1 && temp_childIDs == NULL,
	      "fail should return NULL");
  rc = fc_getFeatureChildIDs(group, numFeature-1, NULL, &temp_childIDs);
  fail_unless(rc != FC_SUCCESS, "should fail if NULL return arg");
  fail_unless(temp_childIDs == NULL, "fail should return NULL");
  rc = fc_getFeatureChildIDs(group, numFeature-1, &temp_numChild, NULL);
  fail_unless(rc != FC_SUCCESS, "should fail if NULL return arg");
  fail_unless(temp_numChild == -1, "fail should return NULL");

  //----------------------------------------------------------
  // test closing
  //----------------------------------------------------------

  fc_freeFeatureGroup(group);
}
END_TEST

// Test error conditions for IO - testing real output done in regression
// tests
START_TEST(featureGroup_io)
{
  FC_ReturnCode rc;
  FC_FeatureGroup* group;
  char* filename = "feature.out";
  
  rc = fc_createFeatureGroup(&group);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create group");
 
  // fc_print/writeFeatureGroup() -- bad args
  rc = fc_printFeatureGroup(NULL);
  fail_unless(rc != FC_SUCCESS, "should fail");
  rc = fc_writeFeatureGroup(NULL, filename);
  fail_unless(rc != FC_SUCCESS, "should fail");
  rc = fc_writeFeatureGroup(group, NULL);
  fail_unless(rc != FC_SUCCESS, "should fail");

  // fc_print/writeFeatureGraph() -- bad args
  rc = fc_printFeatureGraph(NULL);
  fail_unless(rc != FC_SUCCESS, "should fail");
  rc = fc_writeFeatureGraph(NULL, filename);
  fail_unless(rc != FC_SUCCESS, "should fail");
  rc = fc_writeFeatureGraph(group, NULL);
  fail_unless(rc != FC_SUCCESS, "should fail");

  // fc_print/writeFeatureGraph() -- shouldn't work on empty group
  rc = fc_printFeatureGraph(group);
  fail_unless(rc != FC_SUCCESS, "should fail");
  rc = fc_writeFeatureGraph(group, filename);
  fail_unless(rc != FC_SUCCESS, "should fail");

  fc_freeFeatureGroup(group);
}
END_TEST

// Populate the Suite with the tests

Suite *feature_suite(void)
{
  Suite *suite = suite_create("Feature");

  TCase *tc_feature = tcase_create(" - Feature Interface ");
  TCase *tc_featureGroup = tcase_create(" - Feature Group Interface ");

  suite_add_tcase(suite, tc_feature);
  tcase_add_test(tc_feature, feature_basic);
  
  suite_add_tcase(suite, tc_featureGroup);
  tcase_add_test(tc_featureGroup, featureGroup_basic);
  tcase_add_test(tc_featureGroup, featureGroup_io);

  return suite;
}
