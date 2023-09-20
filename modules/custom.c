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
 * \file custom.c
 * \brief Implementation of \ref Custom module.
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/modules/custom.c,v $
 * $Revision: 1.46 $ 
 * $Date: 2006/08/30 19:20:01 $
 */

// C library includes
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <stdio.h>

// fc library dependencies
#include "base.h"
#include "storage.h"
#include "library.h"
#include "dataset.h"
#include "mesh.h"
#include "subset.h"
#include "variable.h"
#include "fileio.h"
#include "geom.h"
#include "topo.h"
#include "feature.h"
#include "track.h"
#include "util.h"
#include "datasetP.h"
#include "meshP.h"

// this module
#include "custom.h"

/**
 * \addtogroup  Custom
 * \brief Sandbox for trying out new stuff.
 *
 * \description
 *
 *   This is a sandbox (i.e. playground) for
 *   custom characterizations that may someday be incorporated into
 *   the fc library.
 *
 * \modifications
 *   - 2003 W Koegler  Created.
 */

/** 
 * \ingroup Custom
 * \brief helper function, generates combinations 
 */
void _fc_generateCombination(int numElement, int ithCombination, int* comb) {
  int i, j;
  for (j = 0; j < numElement; j++)
    comb[j] = 0;
 
  j = 0;
  i = ithCombination;
  while (i > 0) {
    comb[j] = i % 2;
    i /= 2;
    j++;
  }

  return;
}

/**
 * \ingroup  Custom 
 * \brief Compute Score Matrix based on relative amount of overlap.
 *
 * \description
 *
 *    Given a list of previous ROIs and a list of new ROIs, this routine
 *    computes the score matrix used to decide which ROIs belong in the same
 *    features. This uses the number of overlapping entities, not actually
 *    entity size.
 *
 *    The caller is responsible for freeing the dynamically allocated matrix.
 *
 *    Note that if numPrev or numNew are 0, no matrix will be generated.
 *
 *    \todo Should we enforce matching of entity types?

 *   BestMatch Tracking algorithm
 *   Copyright (c) 1999-2000 Laboratory of Visiometrics and Modeling, 
 *   Rutgers University. All rights reserved.

 *
 * \modifications  
 *    - APR-14-2003  W Koegler  Documented
 */
FC_ReturnCode fc_computeBestMatchNumOverlapScore(
  int numPrev,           /**< input - number of segments in previous step*/
  FC_Subset* prevLists,  /**< input - previous ROIs */
  int numNew,            /**< input - number of segments in new step*/
  FC_Subset* newLists,   /**< input - new ROIs */
  double*** scores_p     /**< output - pointer to a numPrev x numNew matrix */
) {
  // assumes unchanging geometry & sorted arrays
  // returns number of entities in common
  FC_ReturnCode rc; 
  int i, j, k;
  double **scores;
  int* newNums;          // number of entities in each new ROI 
  int* prevNums;         // number of entities in each prev ROI
  double** numOverlaps;  

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
  fc_printfLogMessage("Computing Relative NumOverlap Score Matrix.");

  // early return if no overlaps to compute
  if (numPrev == 0 || numNew == 0) 
    return FC_SUCCESS;

  // Get the number of overlaps for all combinations of ROIs
  rc = fc_computeNumOverlapScore(numPrev, prevLists, numNew, newLists,
				 &numOverlaps);
  if (rc != FC_SUCCESS)
    return rc;

  // Get sizes of all ROIs
  newNums = malloc(numNew*sizeof(int));
  prevNums = malloc(numPrev*sizeof(int));
  if (newNums == NULL || prevNums == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  for (i = 0; i < numPrev; i++) {
    rc = fc_getSubsetNumMember(prevLists[i], &prevNums[i]);
    if (rc != FC_SUCCESS)
      return rc;
  }
  for (i = 0; i < numNew; i++) {
    rc = fc_getSubsetNumMember(newLists[i], &newNums[i]);
    if (rc != FC_SUCCESS)
      return rc;
  }

  // create space for the overlap array, initialize to zero
  scores = calloc(numPrev, sizeof(double*));
  if (scores == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  for (i = 0; i < numPrev; i++) {
    scores[i] = calloc(numNew, sizeof(double));
    if (scores[i] == NULL) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      return FC_MEMORY_ERROR;
    }
  }

  // compute forward scores for each previous ROI
  for (i = 0; i < numPrev; i++) {
    // find the new ROIs that this previous ROI overlaps with & keep ids
    int numOverlap = 0;
    int numCombo;
    int* overlapIDs = malloc(numNew*sizeof(int));
    for (j = 0; j < numNew; j++) {
      if (numOverlaps[i][j] > 0) {
	overlapIDs[numOverlap] = j;
	numOverlap++;
      }
    }

    // iterate through all possible combinations of overlaps
    numCombo = (int) pow(2.0, (double)numOverlap);
    for (k = 1; k < numCombo; k++) { // skip 1st comb which is all 0's
      double numIntersect = 0;
      double totalNum = 0;
      double meanNum;
      int cost;
      int* comboMatrix = malloc(numOverlap*sizeof(int));
      _fc_generateCombination(numOverlap, k, comboMatrix);
      for (j = 0; j < numOverlap; j++) {
	if (comboMatrix[j]) {
	  numIntersect += numOverlaps[i][overlapIDs[j]];
	  totalNum += newNums[overlapIDs[j]];
	}
      }
      meanNum = sqrt(prevNums[i]*totalNum);
      //?? cast to int -- just setting arbitrary significant digs
      cost = (int) (numIntersect/meanNum*1000.);
      // for each ROI in this combination, update scores
      for (j = 0; j < numOverlap; j++) {
	if (comboMatrix[j]) {
	  if (cost > scores[i][overlapIDs[j]])
	    scores[i][overlapIDs[j]] = cost;
	}
      }
      free(comboMatrix);
    }
    free(overlapIDs);
  }

  // compute backward scores for each new ROI
  for (j = 0; j < numNew; j++) {
    // find the prev ROIs that this new ROI overlaps with
    int numOverlap = 0;
    int numCombo;
    int* overlapIDs = malloc(numPrev*sizeof(int));
    for (i = 0; i < numPrev; i++) {
      if (numOverlaps[i][j] > 0) {
	overlapIDs[numOverlap] = i;
	numOverlap++;
      }
    }

    // iterate through all possible combinations of overlaps
    numCombo = (int) pow(2.0, (double)numOverlap);
    for (k = 1; k < numCombo; k++) { // skip 1st comb which is all 0's
      double numIntersect = 0;
      double totalNum = 0;
      double meanNum;
      int cost;
      int* comboMatrix = malloc(numOverlap*sizeof(int));
      _fc_generateCombination(numOverlap, k, comboMatrix);
      for (i = 0; i < numOverlap; i++) {
	if (comboMatrix[i]) {
	  numIntersect += numOverlaps[overlapIDs[i]][j];
	  totalNum += prevNums[overlapIDs[i]];
	}
      }
      meanNum = sqrt(newNums[j]*totalNum);
      //?? cast to int -- just setting arbitrary significant digs
      cost = (int) (numIntersect/meanNum*1000.);
      // for each ROI in this combination, update scores
      for (i = 0; i < numOverlap; i++) {
	if (comboMatrix[i]) {
	  if (cost > scores[overlapIDs[i]][j])
	    scores[overlapIDs[i]][j] = cost;
	}
      }
      free(comboMatrix);
    }
    free(overlapIDs);
  }

  // cleanup
  free(prevNums);
  free(newNums);

  // debug
  //printf("Score Matrix\n");
  //for (i = 0; i < numPrev; i++) {
  //for (j = 0; j < numNew; j++)
  //  printf("%f ", scores[i][j]);
  //printf("\n");
  //}

  *scores_p = scores;
  return FC_SUCCESS;
}

/**
 * \ingroup  Custom
 * \brief      Compute correspondence matrix based on best match.
 *
 * \description 
 *
 *    Given the score matrix, create the correspondence matrix which has the
 *    information for assigning ROIs to Features. Chooses best match based
 *    on highest scores where ties are permissible.
 *
 *    The caller is responsible for freeing the dynamically allocated matrix. 
 *
 *    Note that if numPrev or numNew are 0, no matrix will be generated.
 *
 * \modifications  
 *    - 03/04/05 WSK Created
 */
FC_ReturnCode fc_computeBestMatchCorrespond(
  int numPrev,     /**< input - number of segments in previous step*/
  int numNew,      /**< input - number of segments in new step*/
  double** scores, /**< input - the score matrix, numPrev x numNew */
  int*** matrix_p  /**< output - pointer to correspd. matrix, numPrev x numNew */
) {
  int i, j, k;

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
  fc_printfLogMessage("Computing Best Match correspondence matrix.");

  // early return if no overlaps to compute
  if (numPrev == 0 || numNew == 0) { 
    *matrix_p = NULL;
    return FC_SUCCESS;
  }

  // create space for the overlap array
  (*matrix_p) = calloc(numPrev, sizeof(int*));
  if (*matrix_p == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  for (i = 0; i < numPrev; i++) {
    (*matrix_p)[i] = calloc(numNew, sizeof(int));
    if ( (*matrix_p)[i] == NULL) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      return FC_MEMORY_ERROR;
    }
  }

  // for each entry, convert to an int
  for (i = 0; i < numPrev; i++) {
    for (j = 0; j < numNew; j++) {
      // if this is the max value in this row & column, keep
      int maxScore = scores[i][j];
      for (k = 0; k < numPrev; k++) {
	if (scores[k][j] > maxScore)
	  maxScore = scores[k][j];
      }
      for (k = 0; k < numNew; k++) {
	if (scores[i][k] > maxScore)
	  maxScore = scores[i][k];
      }
      if (FC_DBL_EQUIV(scores[i][j], maxScore)) {
	(*matrix_p)[i][j] = 1;
      }
    }
  }

  return FC_SUCCESS;
}
