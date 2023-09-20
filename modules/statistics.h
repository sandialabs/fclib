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
 * \file statistics.h
 * \brief Public declarations for \ref Statistics module.
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/modules/statistics.h,v $
 * $Revision: 1.39 $ 
 * $Date: 2006/09/19 00:57:58 $
 */

#ifndef _FC_STATISTICS_H_
#define _FC_STATISTICS_H_

#ifdef __cplusplus
extern "C" {
#endif


// For a variable
FC_ReturnCode fc_getVariableMinMax(FC_Variable var, double *min, 
                            int *min_index, double *max, int *max_index);
  
FC_ReturnCode fc_getVariableMeanSdev(FC_Variable var, double *mean, 
                            double *sdev);

FC_ReturnCode fc_getVariableSum(FC_Variable var, double* sum);

// For subset of a variable
FC_ReturnCode fc_getVariableSubsetMinMax(FC_Variable var, FC_Subset subset,
                            double *min, int *min_index, 
                            double *max, int *max_index);
  
FC_ReturnCode fc_getVariableSubsetMeanSdev(FC_Variable var, FC_Subset subset,
                            double *mean, double *sdev);

FC_ReturnCode fc_getVariableSubsetSum(FC_Variable var, FC_Subset subset,
                            double* sum);

// For a sequence variable
FC_ReturnCode fc_getSeqVariableMinMax(int numStep, FC_Variable *seqVar,
			    double *min, int *min_seq_id, int* min_entity_id,
			    double *max, int *max_seq_id, int* max_entity_id);
FC_ReturnCode fc_getSeqVariableMeanSdev(int numStep, FC_Variable *seqVar,
			    double *mean, double *sdev);

// For a sequence variable, per data point and component
FC_ReturnCode fc_getSeqVariableSeriesMinMax(int numStep, FC_Variable *seqvar,
                            FC_Variable *minvar, FC_Variable *maxvar,
                            FC_Variable *minindexvar, FC_Variable *maxindexvar);


FC_ReturnCode fc_getSeqVariableSeriesMeanSdev(int numStep, FC_Variable *seqvar, 
                            FC_Variable *meanvar, FC_Variable *sdevvar);

FC_ReturnCode fc_getSeqVariableSeriesSum(int numStep, FC_Variable *seqvar,
                            FC_Variable *sumvar);

// For a sequence
FC_ReturnCode fc_getSequenceMinMaxMono(FC_Sequence seq,
                            double *min, int *min_index,
                            double *max, int *max_index, int *mono);

FC_ReturnCode fc_getSequenceSpacingMinMaxMeanSdev(FC_Sequence seq,
                            double *min, int *min_index,
                            double *max, int *max_index,
                            double *mean, double *std);
                                                   
FC_ReturnCode fc_getClosestSequenceValue(FC_Sequence seq, char*  comparer,
                            double compval, double *seqval, int *step_index);


#ifdef __cplusplus
}
#endif

#endif // _FC_STATISTICS_H_
