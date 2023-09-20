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
 * \file variable.h
 * \brief Public declarations for \ref Variable Module.
 *
 * $Source: /home/Repositories/fcdmf/fclib/modules/variable.h,v $
 * $Revision: 1.30 $ 
 * $Date: 2007/02/28 20:18:28 $
 */

#ifndef _FC_VARIABLE_H_
#define _FC_VARIABLE_H_

#ifdef __cplusplus
extern "C" {
#endif

// create new variable or sequence variable from scratch
//-------------------------------
FC_ReturnCode fc_createGlobalVariable(FC_Dataset dataset, char* varname, 
		       FC_Variable *variable);
FC_ReturnCode fc_createGlobalSeqVariable(FC_Dataset dataset, 
                       FC_Sequence sequence, char* seqvarname, 
		       int *numStep, FC_Variable **seqVar);
FC_ReturnCode fc_createVariable(FC_Mesh mesh, char *varname, 
                       FC_Variable *variable);
FC_ReturnCode fc_createSeqVariable(FC_Mesh mesh, FC_Sequence sequence, 
                       char* seqvarname, int* numStep, 
                       FC_Variable **seqVar);
FC_ReturnCode fc_setVariableData(FC_Variable variable, int numDataPoint, 
		       int numComponent, FC_AssociationType assoc, 
		       FC_MathType mathtype, FC_DataType datatype, void *data);
FC_ReturnCode fc_setVariableDataPtr(FC_Variable variable, int numDataPoint, 
		       int numComponent, FC_AssociationType assoc, 
		       FC_MathType mathtype, FC_DataType datatype, void *data_p);

// other ways to get new variables or sequence variables (see also varmath)
//-------------------------------
FC_ReturnCode fc_copyGlobalVariable(FC_Variable src_var, FC_Dataset dest_ds, 
		       char* newVarName, FC_Variable* new_var);
FC_ReturnCode fc_copyGlobalSeqVariable(int numStep, FC_Variable* src_vars, 
                       FC_Dataset dest_ds, FC_Sequence dest_seq, 
                       char* newSeqVarName, FC_Variable** new_vars);
FC_ReturnCode fc_copyVariable(FC_Variable src_var, FC_Mesh dest_mesh, 
		       char* newVarName, FC_Variable* new_var);
FC_ReturnCode fc_copyVariableFromRegionMesh(FC_Variable src_var, FC_Mesh dest_mesh, 
					    FC_Variable mapping, void* fillval,
					    char* newVarName, FC_Variable* new_var);
FC_ReturnCode fc_copySeqVariableStepFromRegionMesh( FC_Variable src_var,
						    FC_Mesh dest_mesh,
						    FC_Sequence dest_seq,
						    int targetStep,
						    FC_Variable mapping,
						    void* fillval,
						    char* newVarName,
						    FC_Variable** dest_seqvar); 
FC_ReturnCode fc_copyVariableToRegionMesh(FC_Variable src_var, FC_Mesh dest_mesh, 
					    FC_Variable mapping, void* fillval,
					    char* newVarName, FC_Variable* new_var);
FC_ReturnCode fc_copySeqVariable(int numStep, FC_Variable* src_vars, 
                       FC_Mesh dest_mesh, FC_Sequence dest_seq, 
                       char* newSeqVarName, FC_Variable** new_vars);
FC_ReturnCode fc_copyVariableWithNewAssociation(FC_Variable src_var,
		       FC_AssociationType newAssoc,
  		       char* newName, FC_Variable* new_var); 
FC_ReturnCode fc_copySeqVariableWithNewAssociation(int src_numStep,
		       FC_Variable* src_seqVar, 
                       FC_AssociationType newAssoc, char* newName,
  		       FC_Variable** new_seqVar);
FC_ReturnCode fc_convertVariablesToSeqVariable(int numStep, FC_Variable* vars, 
		       FC_Sequence sequence, char* new_name, 
		       FC_Variable** seqVar);
FC_ReturnCode fc_createComponentVariables(FC_Variable var, int* numComponent, 
		       FC_Variable** components);
FC_ReturnCode fc_createComponentSeqVariables(int numStep, FC_Variable* seqVar, 
		       int* numComponent, FC_Variable*** componentSeqVars);
FC_ReturnCode fc_mergeComponentVariables(int numComponent, 
                       FC_Variable* components, char* name, 
		       FC_MathType mathtype, FC_Variable* variable);
FC_ReturnCode fc_mergeComponentSeqVariables(int numComponent, int numStep, 
		       FC_Variable** componentSeqVars, char* name, 
                       FC_MathType mathtype, FC_Variable** seqVariable);

// get variables or sequence variables
//-------------------------------
FC_ReturnCode fc_getNumGlobalVariable(FC_Dataset ds, int *numGlbVar);
FC_ReturnCode fc_getGlobalVariables(FC_Dataset ds, int *numGlbVar, 
		       FC_Variable **glbVars);
FC_ReturnCode fc_getGlobalVariableByName(FC_Dataset ds, char *varName,
					 int *numGlbVar,
					 FC_Variable **glbVars);
FC_ReturnCode fc_getNumGlobalSeqVariable(FC_Dataset ds, int* numGlbSeqVar);
FC_ReturnCode fc_getGlobalSeqVariables(FC_Dataset ds, int *numGlbSeqVar, 
		       int** numStepPerSeq, FC_Variable ***glbSeqVars);
FC_ReturnCode fc_getGlobalSeqVariableByName(FC_Dataset ds, char *seqVarName, 
					    int * numGlbSeqVar,
					    int** numStepPerSeq,
					    FC_Variable ***glbSeqVar);
FC_ReturnCode fc_getNumVariable(FC_Mesh mesh, int *numVar);
FC_ReturnCode fc_getVariables(FC_Mesh mesh, int *numVar, 
		       FC_Variable **variables);
FC_ReturnCode fc_getVariableByName(FC_Mesh mesh, char *varName, 
				   int *numVar,
				   FC_Variable **variables);  
FC_ReturnCode fc_getVariableComponentsByName(FC_Mesh mesh, char *varName,
					     FC_MathType *type,
					     int *numVar,
					     FC_Variable **variables);
FC_ReturnCode fc_getOrGenerateUniqueVariableByName(FC_Mesh mesh, 
						   char *varName, 
						   FC_Variable *variable);  
FC_ReturnCode fc_getNumSeqVariable(FC_Mesh mesh, int* numSeqVar);
FC_ReturnCode fc_getSeqVariables(FC_Mesh mesh, int *numSeqVar, 
		       int** numStepPerSeq, FC_Variable ***seqVariables);
FC_ReturnCode fc_getSeqVariableByName(FC_Mesh mesh, char *seqVarName, 
				      int *numSeqVar,
				      int** numStepPerSeq,
				      FC_Variable ***seqVariables);
FC_ReturnCode fc_getOrGenerateUniqueSeqVariableByName(FC_Mesh mesh, 
						      char *seqVarName, 
						      int* numStep,
						      FC_Variable **seqVariable);
FC_ReturnCode fc_getSeqVariableComponentsByName(FC_Mesh mesh, char *varName,
						FC_MathType *type,
						int *numVar,
						int **numStepPerVar,
						FC_Variable ***seqVariables);


// change the name of a variable or sequence variable
//-------------------------------
FC_ReturnCode fc_changeVariableName(FC_Variable variable, char* newName);
FC_ReturnCode fc_changeSeqVariableName(int numStep, FC_Variable* variable, 
                       char* newName);

// release/delete variables or sequence variables
//-------------------------------
FC_ReturnCode fc_releaseVariable(FC_Variable variable);
FC_ReturnCode fc_releaseSeqVariable(int numStep, FC_Variable *seqVar);
FC_ReturnCode fc_deleteVariable(FC_Variable variable);
FC_ReturnCode fc_deleteSeqVariable(int numStep, FC_Variable *seqVar);

// get meta data
//-------------------------------
int fc_isVariableValid(FC_Variable variable);
int fc_isSeqVariableValid(int numStep, FC_Variable* seqVar);
int fc_isVariableGlobal(FC_Variable variable);
FC_ReturnCode fc_getVariableName(FC_Variable variable, char** varName);
FC_ReturnCode fc_getDatasetFromVariable(FC_Variable variable, FC_Dataset *ds);
FC_ReturnCode fc_getMeshFromVariable(FC_Variable variable, FC_Mesh *mesh);
FC_ReturnCode fc_getSequenceFromSeqVariable(int numStep, FC_Variable* seqVar,
		       FC_Sequence* sequence);
FC_ReturnCode fc_getVariableInfo(FC_Variable variable, int *numDataPoint, 
                       int* numComponent, FC_AssociationType *assoc, 
                       FC_MathType *mathtype, FC_DataType *datatype);
FC_ReturnCode fc_getVariableNumDataPoint(FC_Variable variable, int *count);
FC_ReturnCode fc_getVariableNumComponent(FC_Variable variable, 
                       int *numComponent);
FC_ReturnCode fc_getVariableAssociationType(FC_Variable variable, 
                       FC_AssociationType *assoc);  
FC_ReturnCode fc_getVariableMathType(FC_Variable variable, 
                       FC_MathType *mathtype);
FC_ReturnCode fc_getVariableDataType(FC_Variable variable, 
                       FC_DataType *datatype);

// get big data
//-------------------------------
FC_ReturnCode fc_getVariableDataPtr(FC_Variable variable, void **data_p);
FC_ReturnCode fc_getVariableDataAsDataType(FC_Variable variable, 
		       FC_DataType datatype, void **data_p); 
FC_ReturnCode fc_getSeqVariableDataPointSlice(int numStep, FC_Variable* seqVar,
		       int dataPointID, void** history);
FC_ReturnCode fc_getSeqVariableDataPointSliceAsSeqVariable(int numStep,
			 FC_Variable* seqVar, int dataPointID,
			 char* slicevarname, FC_Variable ** sliceVar);

// Print
FC_ReturnCode fc_printVariable(FC_Variable variable, char* label, 
                       int print_data);
FC_ReturnCode fc_printSeqVariable(int numStep, FC_Variable* seq_vars, 
		       char* label, int print_data);
FC_ReturnCode fc_printSeqVariableRange(int numStep, FC_Variable* seq_vars, 
		       int range_start, int range_stop,
		       char* label, int print_data);


#ifdef __cplusplus
}
#endif

#endif // _FC_VARIABLE_H_
