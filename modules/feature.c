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
 * \file feature.c
 * \brief Implementation for \ref Features module.
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/modules/feature.c,v $
 * $Revision: 1.41 $ 
 * $Date: 2006/08/30 19:20:01 $
 */

// C library includes
#include <stdio.h>
#include <stdlib.h>
#include <string.h> // for memcpy

// fc library dependencies
#include "base.h"
#include "storage.h"
#include "library.h"
#include "subset.h"

// this module
#include "feature.h"
#include "featureP.h"

/**
 * \addtogroup  Features
 * \brief  Basic feature data structures and routines for manipulating them.
 *
 * \description 
 *
 *   <b> ** WARNING ** This interface will definitely change some day ** </b>
 *
 *   The feature module provides the basic feature data structures as
 *   well as routines for manipulating them. 
 *   - ROI - region of interest in a single time step (two pieces of data--the
 *       step ID and an FC_Subset--are used to represent this)
 *   - Features - an ordered collecting of ROIs from successive time steps that
 *     are considered to be the same structure, also knows its parents and 
 *     children features.
 *   - Feature Group - container for related features. It stores ROIs and
 *     creates and stores features.
 *
 *   Usage will usually be through the Feature Group interface so that it can
 *   properly manage ROIs and features. Generally, you add ROIs from successive
 *   time step to the Group and it will create the features.  Features that have
 *   an ROI for the latest time step are considered to be "open" (i.e. still
 *   exists), and conversely, features that do not have an ROI for the latest
 *   time step are considered to have terminated.
 *
 *  These routines purposefully do not check that the ROI's are valid
 *  FC_Subsets. This was done to make testing easier, but also makes
 *  the feature and feature group data structures more flexible.
 *
 * \todo: fundamental question, store ROIs as array of arrays indexed
 *      by step and then local ROI index, or as a single array with
 *      a helper array indexed by step & local ROI index? 
 *
 *      It is nice to have the 2 indices in the feature because
 *      can tell it's sequence extent directly from the stepIDs.
 *      But it might also be nice to have global ROI IDs.
 *      Could do both? can easily keep a total ROI count too
 *
 * \modifications
 *   - FEB-18-2003  W Koegler  Created
 *   - 07/02/04 WSK Removed FC_Roi type and replaced use with FC_Subset. 
 */

/**
 * \ingroup  Features
 * \defgroup PrivateFeatures (Private)
 */

/**
 * \ingroup  Features
 * \defgroup PrivateFeaturesHelpers Helpers (not available outside this module)
 */

/** \name Feature group creation and desctruction. */
//-------------------------------
//@{

/**
 * \ingroup  Features
 * \brief  Create a new Feature Group.
 *
 * \modifications  
 *    - MAR-17-2003  W Koegler  Created
 *    - 11/11/04 WSK Changed so that allocates room for itself.
 */
FC_ReturnCode fc_createFeatureGroup(
  FC_FeatureGroup **group   /**< output - the new feature group */
) {

  // default return
  if (group)
    *group = NULL;

  // check input
  if (!group) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Creating new feature group.");   
  
  // create
  *group = malloc(sizeof(FC_FeatureGroup));
  if (*group == NULL) {
   fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }

  // assign values
  (*group)->numStep = 0;
  (*group)->numROIperStep = NULL;
  (*group)->roiFeatureIDs = NULL;
  (*group)->rois = NULL;
  (*group)->numFeature = 0;
  (*group)->features = NULL;

  return FC_SUCCESS;
}

/**
 * \ingroup  Features
 * \brief  Free memory associated with a feature group
 *
 * \description  
 *
 *    This routine should only be called on a feature group created using
 *    fc_createFeatureGroup.  This frees the dynamically allocated memory
 *    inside the group as well as all features own by the group (but
 *    doesn't delete the ROI's).
 *
 * \modifications  
 *   - MAR-17-2003  W Koegler  Created
 *   - 11/11/04 WSK Changed so that it frees itself as well as it's
 *        dynamically allocated data.
 */
void fc_freeFeatureGroup(
  FC_FeatureGroup *group  /**< input - pointer to feature group */
) {
  int i;

  // early return
  if (group == NULL)
    return;

  // log message
  fc_printfLogMessage("Closing feature group."); 

  // free 2d array of roi handles
  if (group->numStep > 0 && group->numROIperStep && group->rois) {
    for (i = 0; i < group->numStep; i++) {
      if (group->rois[i]) 
        free(group->rois[i]);
    }
    free(group->rois);
  }

  // free 2d allocated roiFeatureIDs
  if (group->numStep > 0 && group->roiFeatureIDs) {
    for (i = 0; i < group->numStep; i++) 
      if (group->roiFeatureIDs[i])
        free(group->roiFeatureIDs[i]);
    free(group->roiFeatureIDs);
  }

  // free allocated features & their array
  if (group->numFeature > 0 && group->features) {
    for (i = 0; i < group->numFeature; i++) {
      if (group->features[i]) {
        _fc_freeFeature(group->features[i]);
      }
    }
    free(group->features);
  }

  // free other arrays
  if (group->numROIperStep)
    free(group->numROIperStep);

  // free itself
  free(group);
}

//@}

/** \name Feature group feature manipulation. */
//-------------------------------
//@{

/**
 * \ingroup  Features
 * \brief  Creates a new feature in the feature group.
 *
 * \description
 *
 *    Creates a new feature whose first ROI is the stepID and the given
 *    subset. If this new feature is a child of a previous feature,
 *    use fc_featureGroup_setParentChildPair() to set that relationship.
 *
 * \modifications  
 *   - MAR-17-2003  W Koegler  Created
 */
FC_ReturnCode fc_featureGroup_newFeature(
  FC_FeatureGroup *group,   /**< input - pointer to feature group */
  int timeID,               /**< input - time ID for the ROI */
  FC_Subset roi,            /**< input - the first ROI in feature */
  int *featureID            /**< output - the new feature's ID */
) {
  FC_ReturnCode rc;
  _FC_Feature* new_feature;
  int roiID;

  // default output
  if (featureID)
    *featureID = -1;

  // check input
  if (!group || timeID < 0 || /*!fc_isSubsetValid(roi) ||*/ !featureID) {
    if (featureID)
      *featureID = -1;
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
 
  // log message
  fc_printfLogMessage("Creating new feature in a feature group.");

  // get a new feature
  rc = _fc_createFeature(1, &new_feature);
  if (rc != FC_SUCCESS)
    return rc;
  
  // add feature to group
  rc = _fc_featureGroup_addFeature(group, new_feature); 
  if (rc != FC_SUCCESS)
    return rc;

  // get feature id
  *featureID = new_feature->id;

  // add ROI to the group & feature
  rc = _fc_featureGroup_addROI(group, *featureID, timeID, roi, &roiID);
  if (rc != FC_SUCCESS) 
    return rc;
  rc = _fc_feature_addROI(new_feature, timeID, roiID);
  if (rc != FC_SUCCESS) 
    return rc;

  return rc;
}

/**
 * \ingroup  Features
 * \brief  Create a Parent/Child Relationship in the feature Group.
 *
 * \description
 *  
 *    Creates a parent/child relationship in the feature Group.
 *    A pointer to the parent is added to the child feature; and
 *    a pointer to the child is added to the parent feature.
 *    Both features must already be in the feature group.
 *
 * \modifications  
 *   - MAR-25-2003  W Koegler  Created
 */
FC_ReturnCode fc_featureGroup_setParentChildPair(
  FC_FeatureGroup *group,   /**< input - pointer to feature group */
  int parentID,             /**< input - the feature's ID */
  int childID               /**< input - the id of the child to be added */
) {
  FC_ReturnCode rc;

  // check input - equal IDs not allowed (circular references!)
  if (!group || parentID == childID || childID < 0 || parentID < 0 ||
      parentID >= group->numFeature || childID >= group->numFeature) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
 
  // log message
  fc_printfLogMessage("Creating parent/child relationship in feature group."); 
 
  // add child
  rc = _fc_feature_addChild(group->features[parentID], childID);
  if (rc != FC_SUCCESS)
    return rc;
  rc = _fc_feature_addParent(group->features[childID], parentID);
 
  return rc;
}

/**
 * \ingroup  Features
 * \brief  Continues a feature
 *
 * \description
 * 
 *    Given a feature, this routine 'continues' the feature by adding the
 *    given ROI. 
 *
 * \modifications  
 *   - MAR-25-2003  W Koegler  Created
 */
FC_ReturnCode fc_featureGroup_continueFeature(
  FC_FeatureGroup *group,    /**< input - pointer to feature group */
  int featureID,             /**< input - ID of the feature to be continued */
  int timeID,                /**< input - time ID for the ROI */
  FC_Subset roi              /**< input - ROI to add to feature */
) {
  FC_ReturnCode rc;
  int roiID;

  // check input 
  if (!group || timeID < 0 || /*!fc_isSubsetValid(roi) || */ 
      featureID >= group->numFeature)  {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Continuing a feature by copying in a ROI.");   
  
  // add new ROI to group
  rc = _fc_featureGroup_addROI(group, featureID, timeID, roi, &roiID);
  if (rc != FC_SUCCESS)
    return rc;

  // Add ROI to feature
  rc = _fc_feature_addROI(group->features[featureID], timeID, roiID);

  return rc;
}

/**
 * \ingroup  Features
 * \brief  Declare that a step in the feature group has no features.
 *
 * \description
 *
 *    This routine is used to add a step to a feature group if there
 *    are no features in that step (steps are automatically added if
 *    there if you use fc_featureGroup_newFeature() or
 *    fc_featureGroup_continueFeature()). 
 *
 *    It is an error to call this routine on a step that already
 *    exists in the feature group.
 *
 * \modifications  
 *   - MAR-17-2003  W Koegler  Created
 */
FC_ReturnCode fc_featureGroup_noFeaturesThisStep(
  FC_FeatureGroup *group,   /**< input - pointer to feature group */
  int stepID                /**< input - step */
) {
  int i;
  int new_numStep, *temp_nums, **temp_featureIDs;
  FC_Subset** temp_rois;

  // check input
  if (!group || stepID < 0 || stepID < group->numStep) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
 
  // log message
  fc_printfLogMessage("Adding featureless step to a feature group.");

  // grow the arrays to accommodate this new step (fill intervening steps
  // with zeros).
  new_numStep = stepID + 1;
  temp_nums = (int*)realloc(group->numROIperStep, new_numStep*sizeof(int));
  temp_rois = (FC_Subset**)realloc(group->rois,
                                   new_numStep*sizeof(FC_Subset*));
  temp_featureIDs = (int**)realloc(group->roiFeatureIDs,
                                   new_numStep*sizeof(int*));
  if (!temp_nums || !temp_rois || !temp_featureIDs) {
    free(temp_nums);
    free(temp_rois);
    free(temp_featureIDs);
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  group->numROIperStep = temp_nums;
  group->rois = temp_rois;
  group->roiFeatureIDs = temp_featureIDs;
  for (i = group->numStep; i < new_numStep; i++) {
    temp_nums[i] = 0;
    temp_rois[i] = NULL;
    temp_featureIDs[i] = NULL;
  }
  group->numStep = new_numStep;

  return FC_SUCCESS;
}

//@}

/** \name Feature group group query. */
//-------------------------------
//@{

/**
 * \ingroup  Features
 * \brief  Get the number of steps in the Feature Group
 *
 * \modifications  
 *    - MAR-31-2003  W Koegler  Created
 */
FC_ReturnCode fc_featureGroup_getNumStep(
  FC_FeatureGroup *group,   /**< input - pointer to feature group */
  int* numStep_p            /**< output - the number of steps in this group */
) {
  // default output
  if (numStep_p)
    *numStep_p = -1;

  // check input
  if (!group || !numStep_p) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // log message
  fc_printfLogMessage("Getting number of steps in a feature group.");  
  // return number of features
  *numStep_p = group->numStep;

  return FC_SUCCESS;
}

/**
 * \ingroup  Features
 * \brief Get the number of Features in the Feature Group
 *
 * \modifications  
 *   - MAR-17-2003  W Koegler  Created
 */
FC_ReturnCode fc_featureGroup_getNumFeature(
  FC_FeatureGroup *group,   /**< input - pointer to feature group */
  int* numFeature_p         /**< output - the number of features in group */
) {
  // default output
  if (numFeature_p)
    *numFeature_p = -1;

  // check input
  if (!group || !numFeature_p) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // log message
  fc_printfLogMessage("Getting number of features in a feature group.");  

  // return number of features
  *numFeature_p = group->numFeature;

  return FC_SUCCESS;
}

//@}

/** \name Query step of a feature group. */
//-------------------------------
//@{

/**
 * \ingroup  Features
 * \brief  Get the number Features in a step of the Feature Group
 *
 * \description
 *  
 *    Returns the number of features that exist at a specified
 *    step in the feature group. Note: this is the same as the number
 *    of ROIs at the same step.
 *
 * \modifications  
 *   - APR-07-2003  W Koegler  Created
 */
FC_ReturnCode fc_featureGroup_getNumFeatureInStep(
  FC_FeatureGroup *group,   /**< input - pointer to feature group */
  int  stepID,              /**< input - the step id */
  int* numFeature_p         /**< output - the number of features in the step */
) {
  // default output
  if (numFeature_p)
    *numFeature_p = -1;

  // check input
  if (!group || stepID < 0 || stepID >= group->numStep || !numFeature_p) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Getting number of features in step %d of the feature "
                      " group.", stepID);    

  // return number of features
  *numFeature_p = group->numROIperStep[stepID];

  return FC_SUCCESS;
}

/**
 * \ingroup  Features
 * \brief  Get the IDs of the features in a step.
 *
 * \description
 *  
 *    Returns the number of features and an array of their IDs 
 *    which exist in the specified step of the feature group.
 *    This array needs to be freed by the calling routine.
 *
 * \modifications  
 *   - MAR-17-2003  W Koegler  Created
 */
FC_ReturnCode fc_featureGroup_getFeatureIDsInStep(
  FC_FeatureGroup *group,   /**< input - pointer to feature group */
  int  stepID,              /**< input - the step id */
  int* numFeature_p,        /**< output - the number of features in the step */
  int **featureIDs_p        /**< output - array of feature IDs */
) {
  // default output
  if (numFeature_p)
    *numFeature_p = -1;
  if (featureIDs_p)
    *featureIDs_p = NULL;
 
  // check input
  if (!group || stepID < 0 || stepID >= group->numStep || !numFeature_p ||
      !featureIDs_p) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // log message
  fc_printfLogMessage("Getting IDs of the existing feature.");  
 
  // get numStep and make array for ids
  *numFeature_p = group->numROIperStep[stepID];
  if (*numFeature_p <= 0) 
    *featureIDs_p = NULL;
  else {
    *featureIDs_p = malloc(sizeof(_FC_Feature*) * (*numFeature_p));
    if (*featureIDs_p == NULL ) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      return FC_MEMORY_ERROR;
    } 
    memcpy(*featureIDs_p, group->roiFeatureIDs[stepID], 
           sizeof(int)*(*numFeature_p));
  }

  return FC_SUCCESS;
}

/**
 * \ingroup  Features
 * \brief  Get the number of ROIs in the Feature Group
 *
 * \modifications  
 *    - 11/12/04 WSK Created.
 */
FC_ReturnCode fc_featureGroup_getNumROI(
  FC_FeatureGroup *group,  /**< input - pointer to feature group */
  int* numROI_p            /**< output - number of ROIs per step */
) {
  int i;

  // default output
  if (numROI_p)
    *numROI_p = -1;

  // check input
  if (!group || !numROI_p) { 
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // log message
  fc_printfLogMessage("Getting number of ROIs in a feature group.");  

  // sum up number of ROI's per step
  *numROI_p = 0;
  for (i = 0; i < group->numStep; i++)
    *numROI_p += group->numROIperStep[i];

  return FC_SUCCESS;
}

/**
 * \ingroup  Features
 * \brief  Get the ROIs in a step from the Feature Group
 *
 * \description
 *  
 *    Returns an array of subset handles that are the ROIs for the
 *    indicated step. Freeing the array is the responsibility of the caller.
 *
 * \modifications  
 *    - MAR-31-2003  W Koegler  Created
 *    - 11/11/04 WSK changed to just get a single steps ROIs
 */
FC_ReturnCode fc_featureGroup_getROIsInStep(
  FC_FeatureGroup *group,  /**< input - the feature group */
  int stepID,              /**< input - the step ID */
  int* numRoi_p,           /**< output - the number of ROIs */
  FC_Subset** rois_p       /**< output - arrays of ROIs */
) {
  int i;

  // default output
  if (numRoi_p)
    *numRoi_p = -1;
  if (rois_p)
    *rois_p = NULL;

  // check input
  if (!group || stepID < 0 || stepID >= group->numStep || !numRoi_p ||
      !rois_p) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // log message
  fc_printfLogMessage("Getting ROIs in step #%d a feature group.", stepID);

  // return rois
  *numRoi_p = group->numROIperStep[stepID];
  if (*numRoi_p > 0) {
    *rois_p = (FC_Subset*)malloc(sizeof(FC_Subset)*(*numRoi_p));
    if (*rois_p == NULL) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      return FC_MEMORY_ERROR;
    } 
    for (i = 0; i < group->numROIperStep[stepID]; i++)
      (*rois_p)[i] = group->rois[stepID][i];
  }

  return FC_SUCCESS;
}

//@}

/** \name Query features. */
//-------------------------------
//@{

/**
 * \ingroup  Features
 * \brief  Get the number of ROIs in a feature.
 *
 * \description
 *
 *    The number of ROI in a feature is equivalent to the number of steps
 *    over which the feature exists.
 *
 * \modifications  
 *    - 11/11/04 WSK, created.
 */
FC_ReturnCode fc_getFeatureNumROI(
  FC_FeatureGroup *group,  /**< input - the feature group */
  int featureID,           /**< input - the feature ID */
  int* numROI_p           /**< output - the number of steps */
) {
  // default output
  if (numROI_p)
    *numROI_p = -1;

  // check input
  if (!group || featureID < 0 || featureID >= group->numFeature || 
      !numROI_p) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // log message
  fc_printfLogMessage("Getting numStep in feature #%d a feature group.", 
                      featureID);

  // return rois
  *numROI_p = group->features[featureID]->numROI;

  return FC_SUCCESS;
}

/**
 * \ingroup  Features
 * \brief  Get the ROIs in a feature.
 *
 * \description
 *  
 *    Returns an array of stepIDs and an array of subset handles that are the 
 *    ROIs for the indicated feature. Freeing the arrays is the responsibility
 *    of the caller. Pass NULL for either array if you don't want it
 *    returned.
 *
 * \modifications  
 *    - 11/11/04 WSK, created.
 */
FC_ReturnCode fc_getFeatureROIs(
  FC_FeatureGroup *group,  /**< input - the feature group */
  int featureID,           /**< input - the feature ID */
  int* numRoi_p,           /**< output - the number of ROIs */
  int** stepIDs_p,         /**< output - (optional) array of step IDs */
  FC_Subset** rois_p       /**< output - (optional) arrays of ROIs */
) {
  int i;
  int stepID, roiID;

  // default output
  if (numRoi_p)
    *numRoi_p = -1;
  if (stepIDs_p)
    *stepIDs_p = NULL;
  if (rois_p)
    *rois_p = NULL;

  // check input
  if (!group || featureID < 0 || featureID >= group->numFeature || 
      !numRoi_p || (!stepIDs_p && !rois_p)) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // log message
  fc_printfLogMessage("Getting ROIs in feature #%d a feature group.", 
                      featureID);

  // return rois
  *numRoi_p = group->features[featureID]->numROI;
  if (*numRoi_p > 0) {
    if (stepIDs_p) {
      *stepIDs_p = (int*)malloc(sizeof(int)*(*numRoi_p));
      if (*stepIDs_p == NULL) {
        fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
        return FC_MEMORY_ERROR;
      }
      for (i = 0; i < *numRoi_p; i++)
        (*stepIDs_p)[i] = group->features[featureID]->stepIDs[i];
    }
    if (rois_p) {
      *rois_p = (FC_Subset*)malloc(sizeof(FC_Subset)*(*numRoi_p));
      if (*rois_p == NULL) {
        fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
        return FC_MEMORY_ERROR;
      } 
      for (i = 0; i < *numRoi_p; i++) {
        stepID = group->features[featureID]->stepIDs[i];
        roiID = group->features[featureID]->roiIDs[i];
        (*rois_p)[i] = group->rois[stepID][roiID];
      }
    }
  }

  return FC_SUCCESS;
}

/**
 * \ingroup  Features
 * \brief  Get the ROI for the specified step.
 *
 * \description
 *  
 *    Returns the ROI that exists in the feature for the specified step.
 *    If no ROI exists at the requested step, function returns an error.
 *
 * \modifications  
 *    - 11/122/04 WSK, created.
 */
FC_ReturnCode fc_getFeatureROIAtStep(
  FC_FeatureGroup *group,  /**< input - the feature group */
  int featureID,           /**< input - the feature ID */
  int stepID,              /**< input - the step ID */
  FC_Subset* roi_p         /**< output - ROI */
) {
  FC_ReturnCode rc;
  int i;
  int numROI, *stepIDs;
  FC_Subset* rois;

  // default output
  if (roi_p)
    *roi_p = FC_NULL_SUBSET;

  // check input
  if (!group || featureID < 0 || featureID >= group->numFeature ||
      stepID < 0 || stepID >= group->numStep || !roi_p) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Getting ROI at stepID %d in feature #%d", 
                      stepID, featureID);

  // return roi
  rc = fc_getFeatureROIs(group, featureID, &numROI, &stepIDs, &rois);
  if (rc != FC_SUCCESS)
    return rc;
  for (i = 0; i < numROI; i++) {
    if (stepIDs[i] == stepID) {
      *roi_p = rois[i];
      break;
    }
  }
  free(stepIDs);
  free(rois);

  if (i == numROI)
    return FC_ERROR;
  else
    return FC_SUCCESS;
}

/**
 * \ingroup  Features
 * \brief  Get the number of parents of a feature.
 *
 * \modifications  
 *    - 11/11/04 WSK, created.
 */
FC_ReturnCode fc_getFeatureNumParent(
  FC_FeatureGroup *group,  /**< input - the feature group */
  int featureID,           /**< input - the feature ID */
  int* numParent_p         /**< output - the number of parents */
) {
  // default output
  if (numParent_p)
    *numParent_p = -1;

  // check input
  if (!group || featureID < 0 || featureID >= group->numFeature || 
      !numParent_p) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // log message
  fc_printfLogMessage("Getting numParent in feature #%d a feature group.", 
                      featureID);

  // return rois
  *numParent_p = group->features[featureID]->numParent;

  return FC_SUCCESS;
}

/**
 * \ingroup  Features
 * \brief  Get the parents of a feature.
 *
 * \description  
 *
 *    Returns IDs (in the same feature group) of the parents of the feature.
 *    Freeing the arrays is the responsibility of the caller.
 *
 * \modifications  
 *    - 11/11/04 WSK, created.
 */
FC_ReturnCode fc_getFeatureParentIDs(
  FC_FeatureGroup *group,  /**< input - the feature group */
  int featureID,           /**< input - the feature ID */
  int* numParent_p,        /**< output - the number of parents */
  int** parentIDs_p        /**< output - the parent IDs */
) {
  int i;

  // default output
  if (numParent_p)
    *numParent_p = -1;
  if (parentIDs_p)
    *parentIDs_p = NULL;

  // check input
  if (!group || featureID < 0 || featureID >= group->numFeature || 
      !numParent_p || !parentIDs_p) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // log message
  fc_printfLogMessage("Getting parents of feature #%d a feature group.", 
                      featureID);

  // return rois
  *numParent_p = group->features[featureID]->numParent;
  if (*numParent_p > 0) {
    *parentIDs_p = (int*)malloc(sizeof(int)*(*numParent_p));
    if (*parentIDs_p == NULL) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      return FC_MEMORY_ERROR;
    }
    for (i = 0; i < *numParent_p; i++)
      (*parentIDs_p)[i] = group->features[featureID]->parentIDs[i];
  }

  return FC_SUCCESS;
}

/**
 * \ingroup  Features
 * \brief  Get the number of children of a feature.
 *
 * \modifications  
 *    - 11/11/04 WSK, created.
 */
FC_ReturnCode fc_getFeatureNumChild(
  FC_FeatureGroup *group,  /**< input - the feature group */
  int featureID,           /**< input - the feature ID */
  int* numChild_p          /**< output - the number of children */
) {
  // default output
  if (numChild_p)
    *numChild_p = -1;

  // check input
  if (!group || featureID < 0 || featureID >= group->numFeature || 
      !numChild_p) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // log message
  fc_printfLogMessage("Getting numChild in feature #%d a feature group.", 
                      featureID);

  // return rois
  *numChild_p = group->features[featureID]->numChild;

  return FC_SUCCESS;
}

/**
 * \ingroup  Features
 * \brief  Get the children of a feature.
 *
 * \description  
 *
 *    Returns IDs (in the same feature group) of the children of the feature.
 *    Freeing the arrays is the responsibility of the caller.
 *
 * \modifications  
 *    - 11/11/04 WSK, created.
 */
FC_ReturnCode fc_getFeatureChildIDs(
  FC_FeatureGroup *group,  /**< input - the feature group */
  int featureID,           /**< input - the feature ID */
  int* numChild_p,         /**< output - the number of children */
  int** childIDs_p         /**< output - the child IDs */
) {
  int i;

  // default output
  if (numChild_p)
    *numChild_p = -1;
  if (childIDs_p)
    *childIDs_p = NULL;

  // check input
  if (!group || featureID < 0 || featureID >= group->numFeature || 
      !numChild_p || !childIDs_p) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // log message
  fc_printfLogMessage("Getting children of feature #%d a feature group.", 
                      featureID);

  // return rois
  *numChild_p = group->features[featureID]->numChild;
  if (*numChild_p > 0) {
    *childIDs_p = (int*)malloc(sizeof(int)*(*numChild_p));
    if (*childIDs_p == NULL) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      return FC_MEMORY_ERROR;
    }
    for (i = 0; i < *numChild_p; i++)
      (*childIDs_p)[i] = group->features[featureID]->childIDs[i];
  }

  return FC_SUCCESS;
}

//@}

/** \name Printing. */
//-------------------------------
//@{

/**
 * \ingroup  Features
 * \brief  Print the contents of a feature group.
 *
 * \description
 *  
 *    Prints human readable information about the features
 *    in a featureGroup.
 *
 * \modifications  
 *    - APR-01-2003  W Koegler  Created
 */
FC_ReturnCode fc_printFeatureGroup(
  FC_FeatureGroup* group   /**< input - pointer to the feature group */
) {
  // check input
  if (!group) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // log message
  fc_printfLogMessage("Printing features.");  

  // do it
  return _fc_writeFeatureGroup(group, stdout);
}

/**
 * \ingroup  Features
 * \brief  Writes the contents of a feature group to file.
 *
 * \description
 *  
 *    Prints human readable information about the features
 *    in a featureGroup into a file.
 *
 * \todo Add a flag fore more info - like the association type
 *    and number of members in each subset.
 *
 * \modifications  
 *    - APR-01-2003  W Koegler  Created
 */
FC_ReturnCode fc_writeFeatureGroup(
  FC_FeatureGroup* group,   /**< input - pointer to the feature group */
  char *filename            /**< input - filename */
) {
  FC_ReturnCode rc;
  FILE* file;

  // check input
  if (!group || !filename) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }  

  // log message
  fc_printfLogMessage("Writing contents of a feature group to a file.");  

  // get the file pointer
  file = fopen(filename, "w");
  if (!file)
    return FC_FILE_IO_ERROR;

  // do it
  rc = _fc_writeFeatureGroup(group, file);

  // close file
  fclose(file);

  return rc;
}

/**
 * \ingroup  Features
 * \brief  Print a feature graph
 *
 * \description
 *  
 *    Prints code that the program 'dot' in ATT&T's graphviz
 *    package can convert to a pretty picture of the graph of the features.
 *
 * \modifications  
 *    - APR-09-2003  W Koegler  Created
 */
FC_ReturnCode fc_printFeatureGraph(
  FC_FeatureGroup* group   /**< input - pointer to the feature group */
) {
  // check input
  if (!group) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
 
  // log message
  fc_printfLogMessage("Printing a feature graph.");  

  // do it
  return _fc_writeFeatureGraph(group, stdout);
}

/**
 * \ingroup  Features
 * \brief  Print a feature graph
 *
 * \description
 *  
 *    Prints code that the program 'dot' in ATT&T's graphviz
 *    package can convert to a pretty picture of the graph of the features.
 *
 * \modifications  
 *    - APR-09-2003  W Koegler  Created
 */
FC_ReturnCode fc_writeFeatureGraph(
  FC_FeatureGroup* group,   /**< input - pointer to the feature group */
  char *filename            /**< input - filename */
) {
  FC_ReturnCode rc;
  FILE *file;
  int numFeature;
 
  // check input
  if (!group || !filename) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }  
  fc_featureGroup_getNumFeature(group, &numFeature);
  if (numFeature < 1) {
    fc_printfErrorMessage("Cannot write graph of empty feature group");
    return FC_INPUT_ERROR;
  }

  // log message
  fc_printfLogMessage("Writing a feature graph.");  

  // open file
  file = fopen(filename, "w");
  if (file == NULL)
    return FC_FILE_IO_ERROR;

  // do it
  rc = _fc_writeFeatureGraph(group, file);

  // close file
  fclose(file);

  return rc;
}

//@}

/**
 * \ingroup  PrivateFeatures
 * \brief Create new (empty) feature and assign ID.
 *
 * \modifications  
 *    - FEB-18-2003  W Koegler  Created
 *    - 11/11/04 WSK Changed so that allocates room for itself.
 *
 * \todo?: Change so have to specify an ROI?
 */
FC_ReturnCode _fc_createFeature(
  int id,               /**< input - the feature's ID */
  _FC_Feature **feature  /**< output - the new feature */
) {
  // default return
  if (feature)
    *feature = NULL;

  // check input
  if (id < 0 || !feature) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Creating new feature.");

  // create
  *feature = malloc(sizeof(_FC_Feature));
  if (*feature == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }

  // assign values
  (*feature)->id = id;
  (*feature)->numROI = 0;
  (*feature)->stepIDs = NULL;
  (*feature)->roiIDs = NULL;
  (*feature)->numParent = 0;
  (*feature)->parentIDs = NULL;
  (*feature)->numChild = 0;
  (*feature)->childIDs = NULL;

  return FC_SUCCESS;
}

/**
 * \ingroup  PrivateFeatures
 * \brief  Free memory associated with a feature.
 *
 * \description
 *
 *    This routine should only be called on features created using
 *    fc_createFeature.  This frees the dynamically allocated memory inside a
 *    feature, then the feature itself.
 *
 * \modifications  
 *   - FEB-18-2003  W Koegler  Created.
 *   - 11/11/04 WSK Changed so that it frees itself as well as it's
 *      dynamically allocated arrays.
 */
void _fc_freeFeature(
  _FC_Feature *feature  /**< input - pointer to the feature */
) {
  // early return
  if (feature == NULL)
    return;
  
  // log message
  fc_printfLogMessage("Freeing feature.");  

  // free allocated arrays
  if (feature->stepIDs)
    free(feature->stepIDs);
  if (feature->roiIDs)
    free(feature->roiIDs);
  if (feature->parentIDs) 
    free(feature->parentIDs);
  if (feature->childIDs)
    free(feature->childIDs);

  // free itself
  free(feature);
}

/**
 * \ingroup  PrivateFeatures
 * \brief  Add an ROI reference to a feature
 *
 * \description  
 *
 *    Adds an ROI to the end of a Feature's ROI list. The ROI is really in the
 *    feature group, here we store a reference to the ROI in the group. The
 *    reference is specified first by the step it is in, and then by its index
 *    in that step.
 *
 * \modifications  
 *   - FEB-18-2003  W Koegler  Created
 */
FC_ReturnCode _fc_feature_addROI(
  _FC_Feature *feature,   /**< input - pointer to the feature */
  int stepID,            /**< input - the index for the step the ROI is in */
  int ROI_ID             /**< input - the ROI's index (in the feature group)
                            within the specified step */
) {
  int i;
  int *temp_steps;
  int *temp_ROIs;

  // check input -- this is an internal routine so we don't check that
  // ROI actually exists on the group, or that ROI doesn't already exits
  // on this feature
  if (!feature || stepID < 0 || ROI_ID < 0) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // log message
  fc_printfLogMessage("Add ROI to a feature.");  

  // expand id arrays
  if (feature->numROI > 0) {
    temp_steps = malloc(sizeof(int)*(feature->numROI + 1));
    temp_ROIs = malloc(sizeof(int)*(feature->numROI + 1));
    if (temp_steps == NULL || temp_ROIs == NULL) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      return FC_MEMORY_ERROR;
    }
    for (i = 0; i < feature->numROI; i++) {
      temp_steps[i] = feature->stepIDs[i];
      temp_ROIs[i] = feature->roiIDs[i];
    }
    free(feature->stepIDs);
    free(feature->roiIDs);
    feature->stepIDs = temp_steps;
    feature->roiIDs = temp_ROIs;
  }
  else {
    feature->stepIDs = malloc(sizeof(int));
    feature->roiIDs = malloc(sizeof(int));
    if (feature->stepIDs == NULL || feature->roiIDs == NULL) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      return FC_MEMORY_ERROR;
    }
    feature->numROI = 0;
  }

  // new values
  feature->stepIDs[feature->numROI] = stepID;
  feature->roiIDs[feature->numROI] = ROI_ID;
  feature->numROI++;

  return FC_SUCCESS;
}

/**
 * \ingroup  PrivateFeatures
 * \brief  Add a parent reference to this feature
 *
 * \description
 *  
 *    Add the specified feature to this feature's parent list. The
 *    reference is the id of the parent in the feature group.
 *
 * \modifications  
 *   - FEB-18-2003  W Koegler  Created
 */
FC_ReturnCode _fc_feature_addParent(
  _FC_Feature *feature,   /**< input - pointer to the feature */
  int featureID          /**< input - the index of the parent feature */
) {
  int i;
  int *temp_ids;
  
  // check input & state (this only catches automatically allocated features)
  if (!feature || featureID < 0)  {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // log message
  fc_printfLogMessage("Adding parent to a feature."); 

  // expand id array
  temp_ids = malloc(sizeof(int)*(feature->numParent + 1));
  if (temp_ids == NULL ) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  } 
  for (i = 0; i < feature->numParent; i++) 
    temp_ids[i] = feature->parentIDs[i];
  free(feature->parentIDs);
  feature->parentIDs = temp_ids;
    
  // new values
  feature->parentIDs[feature->numParent] = featureID;
  feature->numParent++;
  
  return FC_SUCCESS;
}

/**
 * \ingroup  PrivateFeatures
 * \brief  Add a child reference to this feature
 *
 * \description
 *  
 *    Add the specified feature to this feature's child list.  The
 *    reference is the id of the parent in the feature group.
 *
 * \modifications  
 *   - FEB-18-2003  W Koegler  Created
 */
FC_ReturnCode _fc_feature_addChild(
  _FC_Feature *feature,   /**< input - pointer to the feature */
  int featureID          /**< input - the index of the child feature */
){
  int i;
  int *temp_ids;

  // check input
  if (!feature || featureID < 0) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // log message
  fc_printfLogMessage("Adding child to a feature."); 

  // expand id array
  temp_ids = malloc(sizeof(int)*(feature->numChild + 1));
  if (temp_ids == NULL ) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  } 
  for (i = 0; i < feature->numChild; i++) 
    temp_ids[i] = feature->childIDs[i];
  free(feature->childIDs);
  feature->childIDs = temp_ids;
    
  // new values
  feature->childIDs[feature->numChild] = featureID;
  feature->numChild++;

  return FC_SUCCESS;
}

/**
 * \ingroup  PrivateFeaturesHelpers
 * \brief  Helper routine to add a feature to the feature group
 *
 * \description
 *  
 *    Adds a feature to the feature Group. Helper routine to
 *    extend the internal structures of the group. It assigns
 *    the feature an ID and adds the feature to both the 
 *    features array. 
 *    The feature group now has ownership of the passed in 
 *    feature.
 *
 * \modifications  
 *    - MAR-20-2003  W Koegler  Created
 */
FC_ReturnCode _fc_featureGroup_addFeature(
  FC_FeatureGroup *group,        /**< input - pointer to feature group */
  _FC_Feature* feature            /**< input - pointer to a feature */
) {
  int numFeature;
  _FC_Feature** temp_features;

  // check input
  if (!group || !feature)  {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // log message
  fc_printfLogMessage("Adding feature to a feature group."); 
      
  // check feature group state // FIX?
  numFeature = group->numFeature;
  if (numFeature < 0)
    numFeature = 0;
  
  // NOTE: finish setting up feature before touching feature groups' info
  // (this also provides a check on the feature's fitness)

  // set feature ID
  feature->id = numFeature;

  // extend features array
  if (numFeature > 0) {
    temp_features = malloc(sizeof(_FC_Feature*) * (numFeature + 1));
    if (temp_features == NULL) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      return FC_MEMORY_ERROR;
    } 
    memcpy(temp_features, group->features, sizeof(_FC_Feature*)*numFeature);
    free(group->features);
    group->features = temp_features;
  }
  else 
    group->features = malloc(sizeof(_FC_Feature*));
    if (group->features == NULL ) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      return FC_MEMORY_ERROR;
    } 

  // add feature to features array
  group->features[numFeature] = feature;
  group->numFeature = numFeature + 1;

  return FC_SUCCESS;
}

/**
 * \ingroup  PrivateFeaturesHelpers
 * \brief      Helper routine to Add an ROI to the feature group
 *
 * \description
 *  
 *    Adds an ROI to the feature Group. Helper routine to
 *    extend the internal structures of the group. It uses 
 *    the ROI's stepID to place it in the group's structures 
 *    and sets its roiID. It also requires an argument specifying
 *    the id of the owning feature. 
 *    The feature group now has ownership of the passed in ROI.
 *
 * \modifications  
 *   - MAR-19-2003  W Koegler  Created
 */
FC_ReturnCode _fc_featureGroup_addROI(
  FC_FeatureGroup *group,   /**< input - pointer to feature group */
  int featureID,            /**< input - the ID of the feature that owns ROI */
  int stepID,               /**< input - the stepID of the ROI */
  FC_Subset roi,            /**< input - the ROI */
  int* roiID                /**< output - the ID of the roi within the step on 
                               the group */
) {
  int i;
  FC_Subset** temp_rois;
  FC_Subset* temp_step_rois;
  int numROI;
  int new_numStep;
  int* temp_numROIperStep;
  int** temp_roiFeatureIDs;
  int* temp_step_roiFeatureIDs;

  // default output
  if (roiID)
    *roiID = -1;

  // check input
  if (!group || featureID < 0 || stepID < 0 || /*!fc_isSubsetValid(roi) ||*/ 
      !roiID)  {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Adding ROI a feature group.");
  
  // NOTE: finish setting up ROI before touching feature groups' info
  // (this also provides a check on the ROI's fitness)

  // get next roi ID 
  if (stepID < group->numStep)
    *roiID = group->numROIperStep[stepID];
  else
    *roiID = 0; // first in the next step
  
  // if needed, extend numROI & outer level of roiFeatureIDs & rois (over steps)
  if (stepID >= group->numStep) {
    new_numStep = stepID + 1;
    // allocate arrays
    temp_numROIperStep = malloc(sizeof(int)*new_numStep);
    temp_roiFeatureIDs = malloc(sizeof(int*)*new_numStep);
    temp_rois = malloc(sizeof(FC_Subset*)*new_numStep);
    if (temp_numROIperStep == NULL || temp_rois == NULL || 
        temp_roiFeatureIDs == NULL ) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      return FC_MEMORY_ERROR;
    } 
    // copy existing values
    memcpy(temp_numROIperStep, group->numROIperStep, 
           sizeof(int)*group->numStep);
    memcpy(temp_roiFeatureIDs, group->roiFeatureIDs, 
           sizeof(int*)*group->numStep);
    memcpy(temp_rois, group->rois, sizeof(FC_Subset*)*group->numStep);
    // initialize new values
    for (i = group->numStep; i < new_numStep; i++) {
      temp_numROIperStep[i] = 0;
      temp_roiFeatureIDs[i] = NULL;
      temp_rois[i] = NULL;
    }
    // swap pointers
    free(group->numROIperStep);
    group->numROIperStep = temp_numROIperStep;
    free(group->roiFeatureIDs);
    group->roiFeatureIDs = temp_roiFeatureIDs;
    free(group->rois);
    group->rois = temp_rois;
    // set numStep
    group->numStep = new_numStep;
  }

  // extend the inner roiFeatureIDs and rois arrays for this specific step
  numROI = group->numROIperStep[stepID];
  if (numROI < 1) {
    group->roiFeatureIDs[stepID] = malloc(sizeof(int));
    group->rois[stepID] = malloc(sizeof(FC_Subset));
    if (group->roiFeatureIDs[stepID] == NULL || group->rois[stepID] == NULL) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      return FC_MEMORY_ERROR;
    } 
  }
  else {
    temp_step_roiFeatureIDs = malloc(sizeof(int)*(numROI+1));
    temp_step_rois = malloc(sizeof(FC_Subset)*(numROI+1));
    if (temp_step_roiFeatureIDs == NULL || temp_step_rois == NULL) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      return FC_MEMORY_ERROR;
    } 
    memcpy(temp_step_roiFeatureIDs, group->roiFeatureIDs[stepID], 
           sizeof(int)*numROI);
    for (i = 0; i < numROI; i++) 
      temp_step_rois[i] = group->rois[stepID][i];
    free(group->roiFeatureIDs[stepID]);
    group->roiFeatureIDs[stepID] = temp_step_roiFeatureIDs;
    free(group->rois[stepID]);
    group->rois[stepID] = temp_step_rois;
  }
  
  // add the new roi
  // FIX the featureID will be added when it is assigned
  group->roiFeatureIDs[stepID][numROI] = featureID;
  group->rois[stepID][numROI] = roi;
  group->numROIperStep[stepID]++;

  return FC_SUCCESS;
}

/**
 * \ingroup  PrivateFeatures
 * \brief  Get a feature from the feature group
 *
 * \description
 *  
 *    Given a feature's ID, get the feature from a feature group.
 *
 * \modifications  
 *    - MAR-17-2003  W Koegler  Created
 */
FC_ReturnCode _fc_featureGroup_getFeature(
  FC_FeatureGroup *group,   /**< input - pointer to feature group */
  int featureID,            /**< input - a feature's ID */
  _FC_Feature **feature      /**< output - pointer to the feature */
){
  // default output
  if (feature)
    *feature = NULL;

  // check input
  if (!group || !feature) {
     fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Getting feature from a feature group.");  
   
  
  // get return value
  if (featureID >= group->numFeature)  // FIX? should this error?
    *feature = NULL;
  else
    *feature = group->features[featureID];
 
  return FC_SUCCESS;
}

/**
 * \ingroup  PrivateFeatures
 * \brief Get the features in a feature group.
 *
 * \description
 *  
 *    Returns the number of features in the feature group and
 *    an array of pointers to the features in the feature group.
 *    This array needs to be freed by the calling routine.
 *    (But don't free the features, they belong to the group,
 *    free the group when completely done with the features.)
 *
 * \modifications  
 *    - MAR-17-2003  W Koegler  Created
 */
FC_ReturnCode _fc_featureGroup_getFeatures(
  FC_FeatureGroup *group,   /**< input - pointer to feature group */
  int* numFeature_p,        /**< output - number of features in the group */
  _FC_Feature ***features_p  /**< output - array of pointers to features */
) {
  int i;

  // default output
  if (numFeature_p)
    *numFeature_p = -1;
  if (features_p)
    *features_p = NULL;
 
  // check input
  if (!group || !numFeature_p || !features_p) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // log message
  fc_printfLogMessage("Getting features in a feature group.");  
 
  // make arrays
  *numFeature_p = group->numFeature;
  if (*numFeature_p > 0) {
    *features_p = malloc(sizeof(_FC_Feature*) * (*numFeature_p));
    if (*features_p == NULL ) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      return FC_MEMORY_ERROR;
    } 
    for (i = 0; i < *numFeature_p; i++)
      (*features_p)[i] = group->features[i];
  }
  else
    *features_p = NULL;

  return FC_SUCCESS;
}

/**
 * \ingroup  PrivateFeaturesHelpers
 * \brief  Print the contents of a feature group.
 *
 * \description
 *  
 *    Helper routine writes features to the specified FILE
 *    Prints human readable information about the features
 *    in a featureGroup.
 *
 * \modifications  
 *    - APR-01-2003  W Koegler  Created
 */
FC_ReturnCode _fc_writeFeatureGroup(
  FC_FeatureGroup* group,   /**< input - pointer to the feature group */
  FILE* file                /**< input - pointer to the file */
) {
  int i, j;
  int numParent, numChild, numROI, *parentIDs, *childIDs, *stepIDs;
  FC_Subset *ROIs;
  char *temp_name;

  if (!group || !file) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }  
  
  // log message
  fc_printfLogMessage("Printing contents of a feature group.");    

  // try the first one, if ok assume the rest will be ok
  if (fprintf(file, "Feature Group:\n") < 0)
    return FC_FILE_IO_ERROR;
  fprintf(file, "    numStep = %d, numFeature = %d\n", group->numStep,
          group->numFeature);
  fprintf(file, "\n");
  fflush(NULL);

  // for each feature
  for (i = 0; i < group->numFeature; i++) {

    // general feature info, numROI, parents & children
    fprintf(file, "Feature %d:\n", i);
    fflush(NULL);

    // parent ids
    fc_getFeatureParentIDs(group, i, &numParent, &parentIDs);
    fprintf(file, "    numParent = %d, IDs: ", numParent);
    if (numParent > 0) {
      fprintf(file, "%d", parentIDs[0]);
      for (j = 1; j < numParent; j++)
        fprintf(file, ", %d", parentIDs[j]);
    }
    fprintf(file, "\n");
    fflush(NULL);
    free(parentIDs);
      
    // child ids
    fc_getFeatureChildIDs(group, i, &numChild, &childIDs);
    fprintf(file, "    numChild  = %d, IDs: ", numChild);
    if (numChild > 0) {
      fprintf(file, "%d", childIDs[0]);
      for (j = 1; j < numChild; j++)
        fprintf(file, ", %d", childIDs[j]);
    }
    fprintf(file, "\n");
    fflush(NULL);
    free(childIDs);
      
    // rois
    fc_getFeatureROIs(group, i, &numROI, &stepIDs, &ROIs);
    fprintf(file, "    numROI = %d\n", numROI);
    fprintf(file, "    ROIs:\n");
    for (j = 0; j < group->features[i]->numROI; j++) {
      fc_getSubsetName(ROIs[j], &temp_name);
      fprintf(file, "        %d: stepID = %d, ROI = '%s'\n", j, stepIDs[j],
              temp_name);
      fflush(NULL);
      free(temp_name);
    }
    fprintf(file, "\n");
    fflush(NULL);
    free(stepIDs);
    free(ROIs);
  }

  return FC_SUCCESS;
}

/**
 * \ingroup  PrivateFeaturesHelpers
 * \brief  Write a feature graph to a file
 *
 * \description
 *  
 *    Prints code that the program 'dot' in ATT&T's graphviz
 *    package can convert to a pretty picture of the graph of the features.
 *
 * \modifications  
 *    - APR-01-2003  W Koegler  Created
 */
FC_ReturnCode _fc_writeFeatureGraph(
  FC_FeatureGroup* group,   /**< input - pointer to the feature group */
  FILE* file                /**< input - pointer to the file */
) {
  int i, j;
  int numFeature;
  char name[1024];
  char endName[1024];
  char childName[1024];
  char** stepNames;
  _FC_Feature* feature;

  // check input (group already checked?)
  if (!group || !file) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  fc_featureGroup_getNumFeature(group, &numFeature);
  if (numFeature < 1) {
    fc_printfErrorMessage("Cannot write graph of empty feature group");
    return FC_INPUT_ERROR;
  }

  // log message
  fc_printfLogMessage("Writing feature graph to a file.");  

  // NOTE: lots of things must be enclosed in quotes so there are
  //       lots of  \" 's, yuck!

  //---start graph
  // try the first one, if ok assume the rest will be ok
  fprintf(file, "digraph G {\n");
  if (fprintf(file, "orientation = portrait;\n") < 0)
    return FC_FILE_IO_ERROR;
  fprintf(file, "size = \"10,7.5\"; \n");
  fprintf(file, "ranksep = \"0.25\"; \n");
  fprintf(file, "node [height=0.25,width=0.25];\n");
  fflush(NULL);

  //---make step names
  stepNames = malloc(sizeof(char*)*group->numStep);
  if (stepNames == NULL ) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  } 
  for (i = 0; i < group->numStep; i++) {
    stepNames[i] = malloc(sizeof(char)*1024);
    if (stepNames[i] == NULL ) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      return FC_MEMORY_ERROR;
    } 
    sprintf(stepNames[i], "T%d", i);
  }

  //---write time line graph
  fprintf(file, "\n");
  fprintf(file, "{ node [shape=box]; \n");
  fprintf(file, "\"%s\"", stepNames[0]);
  fflush(NULL);
  for (i = 1; i < group->numStep; i++) {
    fprintf(file, " -> \"%s\"", stepNames[i]);
    fflush(NULL);
  }
  fprintf(file, " ; \n");
  fprintf(file, "} \n\n");
  fflush(NULL);

  //---for each feature
  for (i = 0; i < group->numFeature; i++) {
    // dereference feature once
    feature = group->features[i];

    // default names
    sprintf(name, "%d", i);
    sprintf(endName, "%s", name); // same as name
    fflush(NULL);

    //---if feature has no parent, distinguish it's node
    if (feature->numParent == 0)  
      fprintf(file, "\"%s\" [style=filled,color=grey];\n", name);
    
    //---set rank same as it's start step, will show up at same position 
    fprintf(file, "{ rank = same; \"%s\" ; \"%s\" ; }\n", 
           stepNames[feature->stepIDs[0]], name);
    
    //---if feature exists for more than 1 step, make an "end" node
    //   and connect to start node
    if (feature->numROI > 1) {
      sprintf(endName, "%s_end", name);
      // label="" makes it not get a label (default label is node name)
      fprintf(file, "\"%s\" [label=\"\", style=filled,color=black,"
             "shape=diamond,height=0.2,width=0.2] ;\n", endName);
      fprintf(file, "{ rank = same ; \"%s\"; \"%s\" ; }\n", 
             stepNames[feature->stepIDs[feature->numROI-1]], endName);
      fprintf(file, "\"%s\" -> \"%s\"\n", name, endName);
    }
    fflush(NULL);

    //---do child edges
    for (j = 0; j < feature->numChild; j++) {
      sprintf(childName, "%d", feature->childIDs[j]);
      fprintf(file, "\"%s\" -> \"%s\"\n", endName, childName);
      fflush(NULL);
    }
  }

  // end of digraph;
  fprintf(file, "}\n");
  fflush(NULL);

  // free resources
  for (i = 0; i < group->numStep; i++)
    free(stepNames[i]);
  free(stepNames);

  return FC_SUCCESS;
}
