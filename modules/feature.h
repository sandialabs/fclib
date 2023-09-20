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
 * \file feature.h 
 * \brief  Public declarations for \ref Features module.
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/modules/feature.h,v $
 * $Revision: 1.20 $ 
 * $Date: 2006/08/30 19:20:01 $
 */

#ifndef _FC_FEATURE_H_
#define _FC_FEATURE_H_

#ifdef __cplusplus
extern "C" {
#endif

/**
 * \ingroup  PrivateFeatures
 * \brief Feature object
 *
 * \description
 *
 *    This object is set up for use within and by a FC_FeatureGroup -- it's not
 *    really set up to be used on it's own.  Stores ids of ROI's (indices into
 *    a feature group's ROI list).  Stores ids of the child & parent IDs of the
 *    feature.
 *
 * \modifications  
 *   - FEB-18-2003  W Koegler  Created
 */
typedef struct {
  int id;         /**< this feature's id in the feature group */
  int numROI;     /**< number of ROIs in this feature */
  // the stepID & subset locate ROI in space-time
  int *stepIDs;   /**< the step (i.e. time) ids, numROI long */
  int *roiIDs;    /**< the ROI ids (for stepIDs[i]), numROI long */
  // the ids within the feature group
  int numParent;  /**< number of parent features */
  int *parentIDs; /**< parent ids */
  int numChild;   /**< number of child features */
  int *childIDs;  /**< children ids */
} _FC_Feature;

/**
 * \ingroup  Features
 * \brief Feature Group object
 *
 * \description
 *
 *    This object manages the ROIs and the features. An ROI is a space-time
 *    location, i.e., a step and a subset.  The feature group stores the subset
 *    handles of the ROIs (FC_Subset) but does not own them (they may
 *    disappear).  The group creates the features and is responsible for
 *    freeing them.
 *
 * \todo ? add reference to step values (a sequence perhaps)? 
 *
 * \modifications  
 *   - FEB-18-2003  W Koegler  Created
 */
typedef struct {
  // ROI information
  int numStep;         /**< Number of steps this group covers  */
  int* numROIperStep;  /**< Array of number of ROIs in each step */ 
  FC_Subset** rois;    /**< 2D array of ROIs (subsets), 
			  numStep x numROIperStep[numStep] */
  int** roiFeatureIDs; /**< 2D array - for each ROI, the id of feature 
			  containing  it */
  // Feature Information
  int numFeature;         /**< total number of features */ 
  _FC_Feature** features;  /**< 1D array of pointers to features */ 
} FC_FeatureGroup;

// Feature Group creation & destruction
FC_ReturnCode fc_createFeatureGroup(FC_FeatureGroup **group);
void fc_freeFeatureGroup(FC_FeatureGroup *group);

// Feature Group feature manipulation
FC_ReturnCode fc_featureGroup_newFeature(FC_FeatureGroup *group, int stepID, 
			       FC_Subset roi, int *featureID);
FC_ReturnCode fc_featureGroup_setParentChildPair(FC_FeatureGroup* group,
                               int parentID, int childID);
FC_ReturnCode fc_featureGroup_continueFeature(FC_FeatureGroup* group, 
                               int featureID, int stepID, FC_Subset roi);
FC_ReturnCode fc_featureGroup_noFeaturesThisStep(FC_FeatureGroup*group,
			       int stepID);
// Group Query
FC_ReturnCode fc_featureGroup_getNumStep(FC_FeatureGroup* group, int* numStep);
FC_ReturnCode fc_featureGroup_getNumFeature(FC_FeatureGroup *group, 
                               int *numFeature);
FC_ReturnCode fc_featureGroup_getNumROI(FC_FeatureGroup *group, int* numROI);
// Query Step of a group
FC_ReturnCode fc_featureGroup_getNumFeatureInStep(FC_FeatureGroup *group, 
                               int stepID, int *numFeature);
FC_ReturnCode fc_featureGroup_getFeatureIDsInStep(FC_FeatureGroup *group, 
                               int stepID, int *numFeature, int** featureIDs);
FC_ReturnCode fc_featureGroup_getROIsInStep(FC_FeatureGroup* group, 
                               int stepID, int* numROI, FC_Subset** rois);
// Feature Query
FC_ReturnCode fc_getFeatureNumROI(FC_FeatureGroup *group, int featureID, 
                               int* numROI);
FC_ReturnCode fc_getFeatureROIs(FC_FeatureGroup* group, int featureID, 
                               int* numROI, int** stepIDs, FC_Subset** rois);
FC_ReturnCode fc_getFeatureROIAtStep(FC_FeatureGroup* group, int featureID,
				     int stepID, FC_Subset* ROI);
FC_ReturnCode fc_getFeatureNumParent(FC_FeatureGroup *group, int featureID, 
                               int* numParent);
FC_ReturnCode fc_getFeatureParentIDs(FC_FeatureGroup* group, int featureID, 
			       int* numParent, int** parentIDs);
FC_ReturnCode fc_getFeatureNumChild(FC_FeatureGroup *group, int featureID, 
                               int* numChild);
FC_ReturnCode fc_getFeatureChildIDs(FC_FeatureGroup* group, int featureID, 
			       int* numChild, int** childIDs);

// Print
//-------------------------------
FC_ReturnCode fc_printFeatureGroup(FC_FeatureGroup *group);
FC_ReturnCode fc_writeFeatureGroup(FC_FeatureGroup *group, char *filename);
FC_ReturnCode fc_printFeatureGraph(FC_FeatureGroup *group);
FC_ReturnCode fc_writeFeatureGraph(FC_FeatureGroup *group, char *filename);


#ifdef __cplusplus
}
#endif

#endif  // _FC_FEATURE_H_
