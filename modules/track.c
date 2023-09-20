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
 * \file track.c 
 * \brief Implementation of \ref FeatureTracking module.
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/modules/track.c,v $
 * $Revision: 1.43 $ 
 * $Date: 2006/09/19 00:57:58 $
 */

// C library includes
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// fc library dependencies
#include "base.h"
#include "storage.h"
#include "library.h"
#include "subset.h"
#include "feature.h"

// this module
#include "track.h"

/**
 * \addtogroup FeatureTracking
 * \brief  Routines for assembling Regions Of Interest into Features.
 *
 * \description 
 *
 *   The feature tracking module provides routines for assembling Regions of
 *   interest from a series of different steps (usually time steps) into
 *   Features.
 *
 *   In each (time) step, tracking occurs in three steps.
 *   - 1. A score matrix is computed which indicates the amount of
 *        correspondence between each ROI in a new step with the ROI's
 *        in the previous step. 
 *   - 2. The score matrix is analyzed and a correspondence matrix is created
 *        which indicates, by a nonzero value, which of the new ROIs 
 *        correspond with the previous ROIs.
 *   - 3. The correspondence matrix and the ROIs are used to construct the 
 *        features.
 * 
 *   The routine \ref fc_trackAllSteps() takes in the ROIs over a number of
 *   steps and returns a Feature Group containing the features.
 *   The routine \ref fc_trackStep() allows the addition of
 *   a single step.
 *
 *   The default algorithm for computing the score score matrix
 *   (fc_computeNumOverlapScore) is based on the number of overlapping
 *   entities. The default algorithm for computing the correspondence matrix
 *   (fc_computeSimpleCorrespond) is just convert the score matrix to an
 *   int matrix.  The default algorithm for feature assignment
 *   (rc_addROIsUsingBasicEvents) is to start new features if there is a split
 *   or merge event (and of course when there is a new feature).
 *
 *   The routines \ref fc_trackAllSteps_custom() and 
 *   \ref fc_trackStep_custom()  allow the user to provide their own
 *   implementations of each of the three steps in the tracking process. A
 *   value of NULL indicates that the default implementation of that step be
 *   used.  The matrices are accessed by matrix[i][j] where i is the index into
 *   the previous step's ROI list and j is the index into the new step's ROI
 *   list.  Note that a modified Feature Group may be needed if you implement
 *   a different feature assignment step.
 *
 * \todo ? Should the score matrix and correspondence matrix really be separate
 *      routines?
 *
 * \modifications
 *   - FEB-18-2003  W Koegler  Created
 *   - APR-01-2003  W Koegler  Adapted to use feature groups
 *   - APR-14-2003  W Koegler  Integrated into FCLib
 *   - 07/01/04 WSK Changed to take FC_Subsets instead of Linked Lists
 */

/**
 *\ingroup FeatureTracking
 *\defgroup FeatureTrackingPrivateHelpers (Not available outside this module)
 */

void _fc_freeIntMatrix(int numPrev, int numNew, int** matrix);
void _fc_freeDoubleMatrix(int numPrev, int numNew, double** matrix);

/**
 * \ingroup  FeatureTracking
 * \brief      Track and create features for a new step.
 *
 * \description
 *  
 *    For a specific step, track and add the new ROIs to the feature group. 
 *
 * \modifications  
 *    - APR-14-2003  W Koegler  Documented
 */
FC_ReturnCode fc_trackStep(
  int stepIndex,           /**< input - the index of the step */
  int numNew,              /**< input - the number of new segments */
  FC_Subset* newSegments,  /**< input - the new segments to be added */
  FC_FeatureGroup* group   /**< input - feature group to store results */
)
{
	
  // log message
  fc_printfLogMessage("Tracking and creating features for a new step.");
	
  return fc_trackStep_custom(stepIndex, numNew, newSegments, group, NULL, NULL,
			     NULL);
}

/**
 * \ingroup  FeatureTracking
 * \brief      Track and create features.
 *
 * \description
 *
 *    Given the ROIs in each step over a number of steps, this function tracks
 *    the ROIs and creates features which are stored in the feature group.
 *
 * \modifications  
 *    - APR-14-2003  W Koegler  Documented
 */
FC_ReturnCode fc_trackAllSteps(
  int numStep,            /**< input - the number of steps */ 
  int* numSegments,       /**< input - array of number of segments per step */
  FC_Subset** segments,   /**< input - array of arrays of segments per step */
  FC_FeatureGroup* group  /**< input - feature group to store results */
) {

  // log message
  fc_printfLogMessage("Tracking and creating features."); 

  return fc_trackAllSteps_custom(numStep, numSegments, segments, group, 
				 NULL, NULL, NULL);
}

/**
 * \ingroup  FeatureTracking
 * \brief     Track and create features for a new step using custom tracking.
 *
 * \description
 * 
 *    This routine is a customizable version of fc_track_add_step().
 *    The use can specify a function for each of the three steps 
 *    of tracking 
 *      -# compute score matrix; 
 *      -# compute correspondence matrix; and 
 *      -# Assign ROIs to features. 
 *
 *    A NULL value indicates that the
 *    default implementation should be used for that step.
 *
 * \modifications  
 *    - APR-14-2003  W Koegler  Documented
 */
FC_ReturnCode fc_trackStep_custom(
  int stepIndex,           /**< input - the index of the step */
  int numNew,              /**< input - the number of new segments */
  FC_Subset* newSegments,  /**< input - the new segments to be added */
  FC_FeatureGroup* group,  /**< input - feature group to store results */
  int compute_scores(int, FC_Subset*, int, FC_Subset*, double***), 
        /**< input - user provided routine for computing the overlap table. 
	   If NULL is passed, will use the default routine */
  int compute_correspondence(int, int, double**, int***), 
        /**< input - user provided routine for computing the correspondence 
	   matrix. If NULL is passed, will use the default routine */
  int add_rois_to_features(int, int, int, FC_Subset*, int**, FC_FeatureGroup*)
        /**< input - user proved routine for assigning features based on 
	   correspondence matrix. If NULL is passed, will use the default 
	   routine */
)
{
  FC_ReturnCode rc;
  int numPrev;
  double **overlaps;
  int **matrix;
  FC_Subset *prevSegments;

  // check input
  if (stepIndex < 0 || numNew < 0 || (numNew > 0 && !newSegments) || !group) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // log message
  fc_printfLogMessage("Tracking and creating features for a new step "
		      "using custom tracking routines.");

  // assign any default routines
  if (compute_scores == NULL)
    compute_scores = &fc_computeNumOverlapScore;
  if (compute_correspondence == NULL)
    compute_correspondence = &fc_computeSimpleCorrespond;
  if (add_rois_to_features == NULL)
    add_rois_to_features = &fc_addROIsUsingBasicEvents;

  // Get the previous step's ROIs
  if (stepIndex == 0) {
    numPrev = 0;
    prevSegments = NULL;
  }
  else {
    rc = fc_featureGroup_getROIsInStep(group, stepIndex-1, &numPrev,
				       &prevSegments);
    if (rc != FC_SUCCESS) 
      return rc;
  }

  // first compute overlaps
  rc = compute_scores(numPrev, prevSegments, numNew, newSegments, &overlaps);
  free(prevSegments);
  if (rc != FC_SUCCESS) {
    _fc_freeDoubleMatrix(numPrev, numNew, overlaps);
    return rc;
  }

  // second, compute correspondence matrix
  rc = compute_correspondence(numPrev, numNew, overlaps, &matrix);
  _fc_freeDoubleMatrix(numPrev, numNew, overlaps);
  if (rc != FC_SUCCESS) {
    _fc_freeIntMatrix(numPrev, numNew, matrix);
    return rc;
  }

  // third & last, use correspondence matrix to decide on new features
  rc = add_rois_to_features(stepIndex, numPrev, numNew, newSegments, matrix,
			    group);
  _fc_freeIntMatrix(numPrev, numNew, matrix);

  return rc;
}

/**
 * \ingroup  FeatureTracking
 * \brief  Track and create features using custom tracking.
 *
 * \description
 *  
 *    This routine is a customizable version of fc_track().
 *    The use can specify a function for each of the three steps 
 *    of tracking 
 *         -# compute score matrix; 
 *         -# compute correspondence matrix; and 
 *         -# Assign ROIs to features. 
 *
 *    A NULL value indicates that the
 *    default implementation should be used for that step.
 *
 * \modifications  
 *    - APR-14-2003  W Koegler  Documented
 */
FC_ReturnCode fc_trackAllSteps_custom(
  int numStep,            /**< input - the number of steps */ 
  int* numSegments,       /**< input - array of number of segments per step */
  FC_Subset** segments,   /**< input - numStep arrays of segments */
  FC_FeatureGroup* group, /**< input - feature group to store results */
  int compute_scores(int, FC_Subset*, int, FC_Subset*, double***), 
        /**< input - user provided routine for computing the overlap table. 
	   If NULL is passed,  will use the default routine */
  int compute_correspondence(int, int, double**, int***), 
        /**< input - user provided routine for computing the correspondence 
	   matrix. If NULL is passed, will use the default routine */
  int add_rois_to_features(int, int, int, FC_Subset*, int**, FC_FeatureGroup*)
        /**< input - user proved routine for assigning features based on 
	   correspondence matrix. If NULL is passed, will use the default 
	   routine */
) {

  FC_ReturnCode rc;
  int i;

  // check input FIX add
  if (numStep < 1 || numSegments == NULL || segments == NULL || group == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_ERROR;
  }

  // log message
  fc_printfLogMessage("Tracking and creating features using custom tracking.");

  // add the steps
  for (i = 0; i < numStep; i++) {
    rc = fc_trackStep_custom(i, numSegments[i], segments[i], group, 
			     compute_scores, compute_correspondence, 
			     add_rois_to_features);
    if (rc != FC_SUCCESS)
      return FC_ERROR;
  }

  return FC_SUCCESS;
}

/**
 * \ingroup  FeatureTracking
 * \brief Compute Score Matrix based on overlapping entities.
 *
 * \description
 *
 *    Given a list of previous ROIs and a list of new ROIs, this routine
 *    computes the score matrix used to decide which ROIs belong in the same
 *    features. This default routine computes the scores as the number of mesh
 *    entities (e.g. vertices, elements) that the ROIs have in common. It
 *    requires that mesh entity IDs do not change between the steps (i.e.
 *    static topology).
 *
 *    The caller is responsible for freeing the dynamically allocated matrix.
 *
 *    Note that if numPrev or numNew are 0, no matrix will be generated.
 *
 *    \todo Should we enforce matching of entity types?
 *
 * \modifications  
 *    - APR-14-2003  W Koegler  Documented
 */
FC_ReturnCode fc_computeNumOverlapScore(
  int numPrev,           /**< input - number of segments in previous step*/
  FC_Subset* prevLists,  /**< input - previous ROIs */
  int numNew,            /**< input - number of segments in new step*/
  FC_Subset* newLists,   /**< input - new ROIs */
  double*** scores_p     /**< output - pointer to a numPrev x numNew matrix */
) {
  // assumes unchanging geometry & sorted arrays
  // returns number of entities in common
  FC_ReturnCode rc; 
  int i, j;
  double **scores;
  int* arrayA, numA, A;
  int* arrayB, numB, B;

  // default return
  if (scores_p)
    *scores_p = NULL;

  // check input
  if (numPrev < 0 || (numPrev > 0 && ! prevLists) || numNew < 0 ||
      (numNew > 0 && !newLists) || !scores_p)  {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // log message
  fc_printfLogMessage("Computing Score Matrix.");

  // early return if no overlaps to compute
  if (numPrev == 0 || numNew == 0) 
    return FC_SUCCESS;

  // create space for the overlap array
  scores = malloc(sizeof(double*)*numPrev);
  if (scores == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  for (i = 0; i < numPrev; i++) {
    scores[i] = malloc(sizeof(double)*numNew);
    if (scores[i] == NULL) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      return FC_MEMORY_ERROR;
    }
  }
  // for each combination of old and new, compute score
  for (i = 0; i < numPrev; i++) {
    for (j = 0; j < numNew; j++) {
      // walk through both lists at same time:
      // if id's match, increment overlap entry & advance both lists
      // else if old's id is less than new's, advance the old list
      // else advance the new list
      scores[i][j] = 0;
      rc = fc_getSubsetMembersAsArray(prevLists[i], &numA, &arrayA);
      if (rc != FC_SUCCESS)
        return rc;
      rc = fc_getSubsetMembersAsArray(newLists[j], &numB, &arrayB);
      if (rc != FC_SUCCESS)
        return rc;
      A = 0; B = 0;
      while (A < numA && B < numB) { // when hit end of either list
        if (arrayA[A] < arrayB[B])
          A++;
        else if (arrayA[A] > arrayB[B])
          B++; 
	else { // equal
          scores[i][j]++;
          A++;
          B++;
        }
     }
      free(arrayA);
      free(arrayB);
    }
  }

  *scores_p = scores;
  return FC_SUCCESS;
}

/**
 * \ingroup  FeatureTracking
 * \brief      Compute correspondence matrix
 *
 * \description 
 *
 *    Given the score matrix, create the correspondence matrix which has the
 *    information for assigning ROIs to Features. This default implementation
 *    just rounds the score matrix values to ints using the C rint() function
 *    to do the rounding. Rint usually behaves so that values of exactly *.5
 *    are rounded up when * is odd and rounded down when * is even. (Examples:
 *    0.5 => 0, 1.5 => 2, -0.5 => 0, -1.5 => -2.)
 *
 *    The caller is responsible for freeing the dynamically allocated matrix. 
 *
 *    Note that if numPrev or numNew are 0, no matrix will be generated.
 *
 * \modifications  
 *    - APR-14-2003  W Koegler  Documented
 */
FC_ReturnCode fc_computeSimpleCorrespond(
  int numPrev,     /**< input - number of segments in previous step*/
  int numNew,      /**< input - number of segments in new step*/
  double** scores, /**< input - the score matrix, numPrev x numNew */
  int*** matrix_p  /**< output - pointer to correspd. matrix, numPrev x numNew */
) {
  int i, j;

  // default return
  if (matrix_p)
    *matrix_p = NULL;

  // check input
  if (numPrev < 0 || numNew < 0 || (numPrev > 0 && numNew > 0 && !scores) ||
      !matrix_p)  {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
      return FC_INPUT_ERROR;
  }
  if (scores) {
    for (i = 0; i < numPrev; i++) {
      if (!scores[i]) {
	fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR)); 
	return FC_ERROR;
      }
    }
  }

  // log message
  fc_printfLogMessage("Computing correspondence matrix.");

  // early return if no overlaps to compute
  if (numPrev == 0 || numNew == 0) { 
    *matrix_p = NULL;
    return FC_SUCCESS;
  }

  // create space for the overlap array
  (*matrix_p) = malloc(sizeof(int*)*numPrev);
  if (*matrix_p == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  for (i = 0; i < numPrev; i++) {
    (*matrix_p)[i] = malloc(sizeof(int)*numNew);
    if ( (*matrix_p)[i] == NULL) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      return FC_MEMORY_ERROR;
    }
  }
  // for each entry, convert to an int
  for (i = 0; i < numPrev; i++) 
    for (j = 0; j < numNew; j++) 
      (*matrix_p)[i][j] = (int)rint(scores[i][j]);

  return FC_SUCCESS;
}

/**
 * \ingroup  FeatureTracking
 * \brief  Add ROIs to the features in a Feature Group
 *
 * \description
 *  
 *    Given the new ROIs and the correspondence matrix,
 *    create or extend the features in the feature Group.
 *    Use the devinition of basic events to decide on features (i.e.
 *    all events but continue events create new features).
 *
 * \modifications  
 *    - APR-14-2003  W Koegler  Documented
 */
FC_ReturnCode fc_addROIsUsingBasicEvents(
  int stepIndex,          /**< input - index of the new step */
  int numPrev,            /**< input - number of segments in prev step */
  int numNew,             /**< input - number of segments in new step */
  FC_Subset* newSegments, /**< input - an array of ROIs */
  int** overlaps,         /**< input - pointer to a numPrev X numNew */
  FC_FeatureGroup* group  /**< input - feature group to store results */
) {
  FC_ReturnCode rc;	 
  int i, j;
  int idx, id;
  int temp_numExist;     // store numExist before adding this step
  int* temp_existsIDs;   // store the existsIDs before adding this step
  int *numOldOverlaps, *numNewOverlaps;
  FC_Subset roi;

  // check input
  if (stepIndex < 0 || numPrev < 0 || numNew < 0 || 
      (numNew > 0 && !newSegments) || 
      (numPrev > 0 && numNew > 0 && !overlaps) || !group)  {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
      return FC_INPUT_ERROR;
  }
  if (overlaps) {
    for (i = 0; i < numPrev; i++) {
      if (!overlaps[i]) {
	fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR)); 
	return FC_ERROR;
      }
    }
  }

  // log message
  fc_printfLogMessage("Adding ROIs to the Feature Group for step %d.",
		      stepIndex);

  // if no new segments, no new features 
  if (numNew < 1) {
    rc = fc_featureGroup_noFeaturesThisStep(group, stepIndex);
    return rc;
  }

  // get previous step's existing features (should == numPrev)
  temp_numExist = 0;
  if (stepIndex > 0) {
    rc = fc_featureGroup_getFeatureIDsInStep(group, stepIndex-1, 
					     &temp_numExist, &temp_existsIDs);
    if (rc != FC_SUCCESS)
      return rc;
    if (temp_numExist != numPrev) {
      fc_printfErrorMessage("passed in numPrev does not agree with "
			    "feature group");
      free(temp_existsIDs);
      return FC_INPUT_ERROR;
    }
  }

  // if no existing features, all new segments are new features
  if (temp_numExist < 1) {
    for (i = 0; i < numNew; i++) {
      rc = fc_featureGroup_newFeature(group, stepIndex, newSegments[i], &id);
      if(rc != FC_SUCCESS)
        return rc;
    }
    return FC_SUCCESS;
  }

  // make numOverlap helper arrays
  numOldOverlaps = malloc(sizeof(int)*(temp_numExist));
  numNewOverlaps = malloc(sizeof(int)*numNew);
  if (numOldOverlaps == NULL || numNewOverlaps == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  for (i = 0; i < temp_numExist; i++)
    numOldOverlaps[i] = 0;
  for (i = 0; i < numNew; i++)
    numNewOverlaps[i] = 0;
  for (i = 0; i < temp_numExist; i++)
    for (j = 0; j < numNew; j++)
      if (overlaps[i][j] > 0) {
        numOldOverlaps[i]++;
        numNewOverlaps[j]++;
      }

  //temp_numExist = 0;
  // look for new features, merges, splits
  // test number of overlaps to see what to do
  for (i = 0; i < numNew; i++) {
    roi = newSegments[i];
    
    // if no overlaps with previous segments, it's a new feature
    if (numNewOverlaps[i] == 0) { 
      rc = fc_featureGroup_newFeature(group, stepIndex, roi, &id);
      if (rc != FC_SUCCESS)
        return rc;
    }
    
    // if there is one overlap it could be a split or just continue
    else if (numNewOverlaps[i] == 1) { 
      // find old segments index that causes this overlap
      for (j = 0; j < temp_numExist; j++)
        if (overlaps[j][i] > 0) { 
          idx = j;
          break;
        }
      if (numOldOverlaps[idx] == 1) {// this feature just continues
        rc = fc_featureGroup_continueFeature(group, temp_existsIDs[idx],
					     stepIndex, roi);
        if(rc != FC_SUCCESS)
          return rc;					 
      }
      else { // parent has split and this segment is a child
        rc = fc_featureGroup_newFeature(group, stepIndex, roi, &id);
        if(rc != FC_SUCCESS)
          return rc;
        rc = fc_featureGroup_setParentChildPair(group, temp_existsIDs[idx],
						id);
        if(rc != FC_SUCCESS)
          return rc;
      }
    }
    
    // if there are more than one overlaps, it's a merge product
    else { 
      // make a new feature and then add as child to previous features
      rc = fc_featureGroup_newFeature(group, stepIndex, roi, &id);
      if (rc != FC_SUCCESS)
        return rc;
      for (j = 0; j < temp_numExist; j++) 
        if (overlaps[j][i] > 0) {
          idx = (temp_existsIDs)[j];
          rc = fc_featureGroup_setParentChildPair(group, idx, id);
          if (rc != FC_SUCCESS)
            return rc;
        }
    }
  }

  // cleanup
  free(temp_existsIDs);
  free(numOldOverlaps);
  free(numNewOverlaps);

  return FC_SUCCESS;
}

/**
 * \ingroup  FeatureTracking
 * \brief  Group objects based on self-correspondance scores.
 *
 * \description 
 *
 *    Given a self-correspondance matrix (= for N objects, an N X N matrix
 *    where non-zero entry E(i,j) indicates that objects i and j 
 *    "correspond"), returns groups of entities that correspond.
 * 
 *    For example, if a matrix has only non-zero entries on the diagonal,
 *    N groups will be returned with 1 object per group. If all entries
 *    are non-zero, 1 group will be returned with all N objects.  
 *
 *    The caller is responsible for freeing the dynamically allocated
 *    returns: free(numPerGroup) and for (i=0; i<numGroup; i++) {
 *    free(idsPerGroup[i]); } free(idsPerGroup).
 *
 * \modifications  
 *    - 03/14/06 WSD Created
 */
FC_ReturnCode fc_groupSelfCorrespond(
  int num,           /**< input - number of objects */
  int** corresponds, /**< input - the corresponds, num X num, only
		      upper tri required (because it's symmetric) */
  int* numGroup,     /**< output - the number of groups */
  int** numPerGroup, /**< output - the number of objects per group,
			caller is responsible for freeing */
  int*** idsPerGroup /**< output - array of ids of objects per group,
			caller is responsible for freeing */
) {
  int i, j, k;
  int ret;
  int* parents, parentID, temp_num;
  FC_SortedIntArray* lists;

  // default return
  if (numGroup)
    *numGroup = -1;
  if (numPerGroup)
    *numPerGroup = NULL;
  if (idsPerGroup)
    *idsPerGroup = NULL;

  // check input
  if (num < 1 || !corresponds || !numGroup || !numPerGroup || !idsPerGroup)  {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
      return FC_INPUT_ERROR;
  }
  for (i = 0; i < num; i++) {
    if (!corresponds[i]) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR)); 
      return FC_ERROR;
    }
  }

  // log message
  fc_printfLogMessage("Grouping using self corresponds.");

  // initialize data structures
  parents = (int*)malloc(num*sizeof(int));
  lists = (FC_SortedIntArray*)malloc(num*sizeof(FC_SortedIntArray));
  if (!parents || !lists) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR)); 
    return FC_MEMORY_ERROR;
  }
  for (i = 0; i < num; i++) {
    parents[i] = -1; // means "none"
    fc_initSortedIntArray(&lists[i]);
  }

  // do it - keep track of ultimate "parent" for each obj
  //         accumulate groups in list on ultimate parents
  for (i = 0; i < num; i++) {
    // process diagonal entry (assumed to be nonzero)
    parentID = parents[i];
    if (parentID < 0) {
      parentID = i;
      ret = fc_addIntToSortedIntArray(&lists[parentID], i);
      if (ret < FC_SUCCESS)
	return ret;
    }
    // process the rest of entries on this row
    for (j = i+1; j < num; j++) {
      if (corresponds[i][j] > 0) { // only process non-zero entries
	if (parents[j] == -1) {
	  // add j to current parent's list
	  parents[j] = parentID;
	  ret = fc_addIntToSortedIntArray(&lists[parentID], j);
	  if (ret < FC_SUCCESS)
	    return ret;
	}
	else if (parents[j] > parentID) {	
	  // move j's parent's list into current parent's list
	  int tempParent = parents[j];
	  for (k = 0; k < lists[tempParent].numVal; k++)
	    parents[lists[tempParent].vals[k]] = parentID;
	  fc_addIntArrayToSortedIntArray(&lists[parentID], 
					 lists[tempParent].numVal,
					 lists[tempParent].vals, 1);
	  fc_freeSortedIntArray(&lists[tempParent]);
	}
	else if (parents[j] < parentID) {
	  // move current parent's list to j's parent's list
	  int tempParent = parentID;
	  parentID = parents[j];
	  for (k = 0; k < lists[tempParent].numVal; k++)
	    parents[lists[tempParent].vals[k]] = parentID;
	  fc_addIntArrayToSortedIntArray(&lists[tempParent],
					 lists[parentID].numVal,
					 lists[parentID].vals, 1);
	  fc_freeSortedIntArray(&lists[parentID]);
	}
	// if parents[j] == parentID, j already in list, do nothing
      }
    }
  }
  free(parents);

  // count groups
  temp_num = 0;
  for (i = 0; i < num; i++) {
    if (lists[i].numVal > 0)
      temp_num++;
  }

  // Convert group lists to arrays
  *numGroup = temp_num;
  *numPerGroup = (int*)malloc(temp_num*sizeof(int));
  *idsPerGroup = (int**)malloc(temp_num*sizeof(int*));
  if (!(*numPerGroup) || !(*idsPerGroup)) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR)); 
    return FC_MEMORY_ERROR;
  }
  temp_num = 0;
  for (i = 0; i < num; i++) {
    if (lists[i].numVal > 0) {
      (*numPerGroup)[temp_num] = lists[i].numVal;
      (*idsPerGroup)[temp_num] = lists[i].vals;
      temp_num++;
    }
  }
  // don't have to free individual lists 'cause they are passed off
  free(lists);

  return FC_SUCCESS;
}

/**
 * \ingroup  FeatureTrackingPrivateHelpers
 * \brief  Free int matrix
 *
 * \description
 *  
 *    This is a convenience utility to free matrices that are set up
 *    as 2D arrays.
 *
 * \modifications  
 *    - 11/15/04 WSK, changed name and no longer sets freed pointer
 *         to NULL.
 */
void _fc_freeIntMatrix(
  int numPrev,     /**< input - number of segments in previous step */
  int numNew,      /**< input - number of segments in new step*/
  int** matrix_p   /**< input - matrix to free */
) {
  int i;

  if (numPrev > 0 && numNew > 0 && matrix_p) {
    for (i = 0; i < numPrev; i++)
      if (matrix_p[i])
        free(matrix_p[i]);
    free(matrix_p);
  }

  return;
}

/**
 * \ingroup  FeatureTrackingPrivateHelpers
 * \brief  Free double matrix
 *
 * \description
 *  
 *    This is a convenience utility to free matrices that are set up
 *    as 2D arrays.
 *
 * \modifications  
 *    - 11/15/04 WSK, changed name and no longer sets freed pointer
 *         to NULL.
 */
void _fc_freeDoubleMatrix(
  int numPrev,        /**< input - number of segments in previous step*/
  int numNew,         /**< input - number of segments in new step*/
  double** matrix_p   /**< input - matrix to free */
) {
  int i;

  if (numPrev > 0 && numNew > 0 && matrix_p) {
    for (i = 0; i < numPrev; i++)
      if (matrix_p[i])
        free(matrix_p[i]);
    free(matrix_p);
  }

  return;
}

