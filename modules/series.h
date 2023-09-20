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
 * \file series.h
 * \brief Public declarations for \ref Series module.
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/modules/series.h,v $
 * $Revision: 1.47 $
 * $Date: 2006/09/19 00:57:58 $
 *
 */

#ifndef _FC_SERIES_H_
#define _FC_SERIES_H_

#ifdef __cplusplus
extern "C" {
#endif

//************************************************************
//these map a single sequence var into another single sequence var
//************************************************************
FC_ReturnCode fc_leadingWindowAverage(int numStep, FC_Variable *seqvar,
                                      int windowsize, char* varname,
                                      FC_Variable **newseqvar);

FC_ReturnCode fc_leadingWindowAverage_Time(int numStep, FC_Variable *seqvar,
                                      double time_val, char* varname,
                                      FC_Variable **newseqvar);

FC_ReturnCode fc_centeredWindowAverage(int numStep, FC_Variable *seqvar,
                                      int windowsize, char* varname,
                                      FC_Variable **newseqvar);

FC_ReturnCode fc_centeredWindowAverage_Time(int numStep, FC_Variable *seqvar,
                                      double time_val, char* varname,
                                      FC_Variable **newseqvar);

FC_ReturnCode fc_normalForm(int numStep, FC_Variable *seqvar,
			    char* varname,
			    FC_Variable **newseqvar);

FC_ReturnCode fc_linearInterpolation(int numStep, FC_Variable *seqvar,
				      FC_Sequence newseq,
				      char* varname,
				      FC_Variable **newseqvar);

FC_ReturnCode fc_firstDerivative_REA(int numStep, FC_Variable *seqvar,
				     double ueps,
				      char* varname,
				      FC_Variable **newseqvar);


//************************************************************
//these operate on a single sequence var and return a value or values
//based on that var but the return val is *not* a sequence var
//************************************************************
FC_ReturnCode fc_integral_TR(int numStep, FC_Variable *seqvar1,
					 char* varname, FC_Variable *newvar);


FC_ReturnCode fc_linearLeastSquares(int numStep, FC_Variable *seqvar,
				FC_Variable *newVar);


//************************************************************
//these map a pair of sequence vars into a single sequence var
//************************************************************



//************************************************************
//these operate on a pair of sequence vars and return a value or values
//based on those vars but the return val is *not* a sequence var
//************************************************************
FC_ReturnCode  fc_euclideanDistanceBetweenCurves(int numStep, 
					     FC_Variable *seqvar1,
					     FC_Variable *seqvar2,
					     char* varname,
					     FC_Variable *newvar);


FC_ReturnCode fc_dimensionlessAreaBetweenCurves(int numStep, 
					     FC_Variable *refseqvar,
					     FC_Variable *compseqvar,
					     char* varname,
					     FC_Variable *newvar);


//************************************************************
//to be put in
//1) convolution
//2) shift and scale in the (time)series variable



#ifdef __cplusplus
}
#endif

#endif // _FC_SERIES_H_
