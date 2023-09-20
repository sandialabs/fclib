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
 * \file varmath.h
 * \brief Public declarations for \ref VariableMath module.
 *
 * $Source: /home/Repositories/fcdmf/fclib/modules/varmath.h,v $
 * $Revision: 1.39 $
 * $Date: 2006/11/20 20:06:38 $
 */ 

#ifndef _FC_VARMATH_H_
#define _FC_VARMATH_H_

#ifdef __cplusplus
extern "C" {
#endif

//---- Built-in operations on vars
FC_ReturnCode fc_operatorVar(FC_Variable var1, char* operation,
                      char* new_var_name, FC_Variable* resultVar);
FC_ReturnCode fc_varOperatorVar(FC_Variable var1, char* operation, 
                      FC_Variable var2, char* new_var_name, 
                      FC_Variable* new_var);
FC_ReturnCode fc_varOperatorConst(FC_Variable var, char* operation,      
                      int constNumComp, void* constVector, 
                      FC_DataType constType, char* new_var_name, 
                      FC_Variable* new_var);
FC_ReturnCode fc_constOperatorVar(int constNumComp, void* constVector,
                      FC_DataType constType, char* operation, FC_Variable var, 
		      char* new_var_name, FC_Variable* new_var) ;

//---- Built-in operations on seq vars
FC_ReturnCode fc_operatorSeqVar(int numStep, FC_Variable *seqvar1,
		      char* operation, char* new_seq_var_name,
		      FC_Variable** new_seq_var);
FC_ReturnCode fc_seqVarOperatorSeqVar(int numstep, FC_Variable *seqvar1,
		      char* operation, FC_Variable *seqvar2, 
		      char* new_seq_var_name, FC_Variable** new_seq_var);
FC_ReturnCode fc_seqVarOperatorSeqConst(int numstep, FC_Variable *seqvar1,
                      char* operation, int constNumComp, void* constVectors,
                      FC_DataType constType, char* new_seq_var_name, 
                      FC_Variable** new_seq_var);
FC_ReturnCode fc_seqConstOperatorSeqVar(int numstep, int constNumComp,
                      void* constVectors, FC_DataType constType, 
                      char* operation, FC_Variable *seqvar2, 
                      char* new_seq_var_name, FC_Variable** new_seq_var);

//---- Built-in operations between vars and seq vars
FC_ReturnCode fc_seqVarOperatorVar(int numstep, FC_Variable *seqvar1,
		      char* operation, FC_Variable var2, 
		      char* new_seq_var_name, FC_Variable** new_seq_var);
FC_ReturnCode fc_varOperatorSeqVar(FC_Variable var1, char* operation, 
		      int numstep, FC_Variable *seqvar2, 
		      char* new_seq_var_name, FC_Variable** new_seq_var);
FC_ReturnCode fc_seqVarOperatorConst(int numstep, FC_Variable *seqvar1,
		      char* operation, int constNumComp, void* constVector, 
                      FC_DataType constType, 
		      char* new_seq_var_name, FC_Variable** new_seq_var);
FC_ReturnCode fc_constOperatorSeqVar(int constNumComp, void* constVector, 
                      FC_DataType constType, char* operation, 
  	              int numstep, FC_Variable *seqvar2, 
  	              char* new_seq_var_name, FC_Variable** new_seq_var);
   
//---- functions on vars
FC_ReturnCode fc_varUnaryFunction(FC_Variable var, 
                      double (*unaryFunction)(double),       
                      char* new_var_name, FC_Variable* new_var);
FC_ReturnCode fc_varBinaryFunctionVarVar(FC_Variable var1, FC_Variable var2,
                      double (*binaryFunction)(double, double),
                      char* new_var_name, FC_Variable* new_var);
FC_ReturnCode fc_varBinaryFunctionVarConst(FC_Variable var, 
		      int numConst, double* constValues,
                      double (*binaryFunction)(double, double),
                      char* new_var_name, FC_Variable* new_var);
FC_ReturnCode fc_varBinaryFunctionConstVar(int constLen, double* constValues,
                      FC_Variable var,
                      double (*binaryFunction)(double, double),
                      char* new_var_name, FC_Variable* new_var);

//---- functions on seq vars
FC_ReturnCode fc_seqVarUnaryFunction(int numStep, FC_Variable *seqvar, 
                      double (*unaryFunction)(double),       
                      char* new_var_name, FC_Variable** new_seq_var);
FC_ReturnCode fc_seqVarBinaryFunctionSeqVarSeqVar(int numStep,
		      FC_Variable *seqvar1, FC_Variable *seqvar2,
		      double (*binaryFunction)(double, double),
		      char* new_seq_var_name, FC_Variable** new_seq_var);

//---- functions on seq vars and constants
FC_ReturnCode fc_seqVarBinaryFunctionSeqVarConst(
                      int numStep, FC_Variable *seq_var,
		      int constNumComp, double* constValues, 
		      double (*binaryFunction)(double ddata1, double ddata2),
		      char* new_seq_var_name, FC_Variable** new_seq_var);
FC_ReturnCode fc_seqVarBinaryFunctionConstSeqVar(
		      int constNumComp, double* constValues,		      
		      int numStep, FC_Variable *seq_var,
		      double (*binaryFunction)(double ddata1, double ddata2),
		      char* new_seq_var_name, FC_Variable** new_seq_var);
FC_ReturnCode fc_seqVarBinaryFunctionSeqVarSeqConst(
		      int numStep, FC_Variable *seq_var,
		      int constNumComp, void* constVectors, FC_DataType  constType,
		      double (*binaryFunction)(double ddata1, double ddata2),
		      char* new_seq_var_name,  FC_Variable** new_seq_var);
FC_ReturnCode fc_seqVarBinaryFunctionSeqConstSeqVar(
                      int numStep,
		      int constNumComp,  void* constVectors, FC_DataType  constType,
		      FC_Variable *seq_var,		      
		      double (*binaryFunction)(double ddata1, double ddata2),
		      char* new_seq_var_name,
		      FC_Variable** new_seq_var);

//---- functions between vars and seq vars
FC_ReturnCode fc_seqVarBinaryFunctionSeqVarVar(int numStep,
		      FC_Variable *seqvar1, FC_Variable var2,
		      double (*binaryFunction)(double, double),
		      char* new_seq_var_name, FC_Variable** new_seq_var);
FC_ReturnCode fc_seqVarBinaryFunctionVarSeqVar(FC_Variable var1, int numStep,
		      FC_Variable *seqvar2,
		      double (*binaryFunction)(double, double),
		      char* new_seq_var_name, FC_Variable** new_seq_var);

//------ Other --- more component aware
FC_ReturnCode fc_createMagnitudeVariable(FC_Variable var, char* name, 
                      FC_Variable* magnitude_var);
// make a createDeterminant version for "tensor" variables?
FC_ReturnCode fc_createNormalTangentVariables(FC_Variable inputVector, 
                      int numDim, double* referenceVector, 
                      FC_Variable* normalComponent, 
                      FC_Variable* tangentComponent);
FC_ReturnCode fc_createNormalTangentVariables2(FC_Variable inputVector, 
                      FC_Variable referenceVector, 
                      FC_Variable* normalComponent, 
                      FC_Variable* tangentComponent);
FC_ReturnCode fc_createStressPressureVariables(FC_Variable* sigmaVars, 
		      FC_Variable* stressVar, FC_Variable* pressureVar);

 
#ifdef __cplusplus
}
#endif

#endif  // _FC_VARMATH_H_
