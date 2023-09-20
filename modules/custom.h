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
 * \file custom.h
 * \brief Declarations for \ref Custom module.
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/modules/custom.h,v $
 * $Revision: 1.33 $ 
 * $Date: 2006/08/30 19:20:01 $
 */

#ifndef _FC_CUSTOM_H_
#define _FC_CUSTOM_H_

#ifdef __cplusplus
extern "C" {
#endif

// specialized algorithms for tracking. These need to be tested and
// perhaps they should be made the default algorithms.
FC_ReturnCode fc_computeBestMatchNumOverlapScore(int numPrev, 
		       FC_Subset* prevSegments, int numNew, 
		       FC_Subset* newSegments, double*** scores);
FC_ReturnCode fc_computeBestMatchCorrespond(int numPrev, int numNew, 
                       double** scores, int*** matrix);

void _fc_generateCombination(int numElement, int ithCombination, int* comb);

#ifdef __cplusplus
}
#endif

#endif  // end of _FC_CUSTOM_H_
