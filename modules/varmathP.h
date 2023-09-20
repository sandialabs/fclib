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
 * \file varmathP.h
 * \brief Private declarations for \ref VariableMath module.
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/modules/varmathP.h,v $
 * $Revision: 1.8 $
 * $Date: 2006/08/30 19:20:01 $
 */ 

#ifndef _FC_VARMATH_P_H_
#define _FC_VARMATH_P_H_

#ifdef __cplusplus
extern "C" {
#endif

// unary operations on arrays
FC_ReturnCode _fc_negateInts(int num, int* data, int* resultData);
FC_ReturnCode _fc_negateDoubles(int num, double* data, double* resultData);

// binary operations on arrays
FC_ReturnCode _fc_addInts(int numComp, int numDataPoint1, int* data1, 
		 int numDataPoint2, int* data2, int* resultData);
FC_ReturnCode _fc_subtractInts(int numComp, int numDataPoint1, int* data1, 
		 int numDataPoint2, int* data2, int* resultData);
FC_ReturnCode _fc_multiplyInts(int numComp, int numDataPoint1, int* data1, 
		 int numDataPoint2, int* data2, int* resultData);
FC_ReturnCode _fc_divideInts(int numComp, int numDataPoint1, int* data1, 
		 int numDataPoint2, int* data2, int* resultData);
FC_ReturnCode _fc_addDoubles(int numComp, int numDataPoint1, double* data1, 
		 int numDataPoint2, double* data2, double* resultData);
FC_ReturnCode _fc_subtractDoubles(int numComp, int numDataPoint1, 
                 double* data1, int numDataPoint2, double* data2, 
                 double* resultData);
FC_ReturnCode _fc_multiplyDoubles(int numComp, int numDataPoint1, 
                 double* data1, int numDataPoint2, double* data2, 
                 double* resultData);
FC_ReturnCode _fc_divideDoubles(int numComp, int numDataPoint1, double* data1, 
		 int numDataPoint2, double* data2, double* resultData);

// binary operator support
FC_ReturnCode _fc_binaryOperatorCore(int numComp, int numData1, void* data1,
		 FC_DataType datatype1, char* operation, int numData2,
		 void* data2, FC_DataType datatype2, FC_Variable origVar, 
                 char* new_var_name, FC_Variable * new_var);

// decompose vector support
void _fc_decomposeVector(int inputDim, double* inputVector, int refDim,
                 double* referenceVector, double* normalComponent, 
                 double* tangentComponent);

#ifdef __cplusplus
}
#endif

#endif  // _FC_VARMATH_P_H_
