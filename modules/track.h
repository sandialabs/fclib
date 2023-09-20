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
 * \file track.h 
 * \brief Declaration of tracking routines
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/modules/track.h,v $
 * $Revision: 1.22 $ 
 * $Date: 2006/08/30 19:20:01 $
 */

#ifndef _FC_TRACK_H_
#define _FC_TRACK_H_

#ifdef __cplusplus
extern "C" {
#endif

// Tracking routines
FC_ReturnCode fc_trackStep(int stepIndex, int numNew, FC_Subset* newSegments,
                       FC_FeatureGroup* featureGroup);
FC_ReturnCode fc_trackAllSteps(int numStep, int* numSegments, 
                       FC_Subset** segments, FC_FeatureGroup* featureGroup);

// Customizable tracking routines
FC_ReturnCode fc_trackStep_custom(int stepIndex, int numNew, 
		       FC_Subset* newSegments, FC_FeatureGroup* featureGroup,
		       int compute_scores(int, FC_Subset*, int, FC_Subset*, 
					  double***),
                       int compute_correspondence(int, int, double**, int***),
		       int add_rois_to_features(int, int, int, FC_Subset*,
						int**, FC_FeatureGroup*));
FC_ReturnCode fc_trackAllSteps_custom(int numStep, int* numSegments, 
                       FC_Subset** segments, FC_FeatureGroup* featureGroup,
                       int compute_scores(int, FC_Subset*, int, FC_Subset*, 
					  double***),
                       int compute_correspondence(int, int, double**, int***),
		       int add_rois_to_features(int, int, int, FC_Subset*, 
                                                int**, FC_FeatureGroup*));

// Default routines for tracking
FC_ReturnCode fc_computeNumOverlapScore(int numPrev, FC_Subset* prevSegments, 
		       int numNew, FC_Subset* newSegments, double*** scores);
FC_ReturnCode fc_computeSimpleCorrespond(int numPrev, int numNew, 
                       double** scores, int*** matrix);
FC_ReturnCode fc_addROIsUsingBasicEvents(int stepIndex, int numPrev, 
		       int numNew, FC_Subset* newSegments, int** matrix, 
                       FC_FeatureGroup* featureGroup);

// Related routines

// Given a self-correspondance matrix for a group of objects, figure out
// how they should be grouped
FC_ReturnCode fc_groupSelfCorrespond(int num, int** corresponds, int* numGroup,
				     int** numPerGroup, int***idsPerGroup);


#ifdef __cplusplus
}
#endif

#endif // _FC_TRACK_H_
