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
 * \file featureP.h 
 * \brief  Private declarations for \ref Features module.
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/modules/featureP.h,v $
 * $Revision: 1.4 $ 
 * $Date: 2006/08/30 19:20:01 $
 */

#ifndef _FC_FEATURE_P_H_
#define _FC_FEATURE_P_H_

#ifdef __cplusplus
extern "C" {
#endif

// features

// create/free
FC_ReturnCode _fc_createFeature(int id, _FC_Feature **feature);
void _fc_freeFeature(_FC_Feature *feature);
// helpers
FC_ReturnCode _fc_feature_addROI(_FC_Feature *feature, int stepID, 
                               int ROI_ID); 
FC_ReturnCode _fc_feature_addParent(_FC_Feature *feature, int parentID);
FC_ReturnCode _fc_feature_addChild(_FC_Feature *feature, int childID);

// feature group

// helpers
FC_ReturnCode _fc_featureGroup_addFeature(FC_FeatureGroup *group, 
                               _FC_Feature* feature);
FC_ReturnCode _fc_featureGroup_addROI(FC_FeatureGroup *group, int featureID,
                               int stepID, FC_Subset roi, int* roiID);
// query
FC_ReturnCode _fc_featureGroup_getFeature(FC_FeatureGroup *group, 
			       int featureID, _FC_Feature **feature);
FC_ReturnCode _fc_featureGroup_getFeatures(FC_FeatureGroup *group, 
                               int* numFeature, _FC_Feature ***features);
// core of IO
FC_ReturnCode _fc_writeFeatureGroup(FC_FeatureGroup *group, FILE* file);
FC_ReturnCode _fc_writeFeatureGraph(FC_FeatureGroup *group, FILE* file);

#ifdef __cplusplus
}
#endif

#endif  // _FC_FEATURE_P_H_
