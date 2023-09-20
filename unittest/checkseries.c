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
 * \file checkseries.c
 * \brief Unit testing of \ref series module
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/unittest/checkseries.c,v $
 * $Revision: 1.63 $
 * $Date: 2006/10/19 03:14:53 $
 *
 * \description
 *
 *    Tests the series routines.
 *    Note: these tests use fc_createRegularSequence in the Sequence module
 *
 * \todo
 *    Decide if want to combine window and single series checks
 *    Decide if want to move comparison seq vectors into global.
 *    Waiting on this last to see if can reuse values.
 *
 * \modifications
 *    -1/12/05 ACG  Created.
 *    -1/13/05 ACG  changing test to made up dataset
 *    -1/31/05 ACG  moved helper function _fc_convertSeqVarType
 *                  into timeseries now named 
 *                  fc_createChangedDataTypeSeqVariable
 *    -2/10/05 ACG  and createChangedDataTypeSeqVariable is back. There is now
 *                  fc_getVariableDataAsDataType as a more elemental piece of
 *                  this
 *    -2/24/05 ACG  moved some setup things into the glob setup. 
 *                  this affects everythign except regtoseq since its up in
 *                  the air whether that will be kept or not. 
 *    -4/12/05 ACG  regtoseq is becoming a real check again since need to
 *                  add in ability to check two different sequence variables
 *                  by putting them on the same sequence. this will require
 *                  adding in a linear interpolation function. this test
 *                  will be updated once the linear interpolation is written.
 *    -4/13/05 ACG  adding some more glob sequences for different types
 *    -4/13/05 ACG  putting in linearInterpolation as a separate check
 *                  and leaving reg to seq alone. may rip out reg to seq....
 *    -4/21/05 ACG  moving createRegularSequence to sequence module. copying
 *                  test for that over to Sequence. will come back and
 *                  rip out the reg to seq test thats in here, once im
 *                  sure the move to Sequence module is ok
 *    -4/26/05 ACG  decided to leave in regtoseqcheck (though removing the
 *                  parts that checked creatign the regular sequence since
 *                  thats in the sequence module) even thoguh it just now
 *                  checks a private function since i use it to make vars
 *                  for testing alot.
 *    -6/1/05 ACG   name change from timeseries to series
 *
 */
#include <stdlib.h>
#include <check.h>
#include "fc.h"
#include "fcP.h"
#include "seriesP.h"
#include "checkall.h"

//globals
static FC_Dataset glob_dataset;
static FC_Mesh glob_mesh;
static int glob_meshDataPoint =3;
static int glob_meshDim = 1;
static double glob_meshcoords[3]={0.0,1.0,2.0}; 

static int glob_numStep =5;
static FC_DataType glob_seqDataType = FC_DT_FLOAT;
static FC_Sequence glob_fseq;
static float glob_ftimes[5] = {0.0,1.0,2.0,3.0,4.0};
static double glob_dtimes[5] = {0.0,1.0,2.0,3.0,4.0};
static int glob_itimes[5] = {0,1,2,3,4};
static char glob_ctimes[5] = {'a','b','c','d','e'}; 
static FC_Sequence glob_cseq;
static FC_Sequence glob_dseq;
static FC_Sequence glob_iseq;

static FC_Variable* glob_sv1;
static char* glob_sv1_name = "glob_sv1";
static FC_Variable* glob_svchar;
static char* glob_svchar_name = "glob_svchar";

static FC_DataType glob_dtype[3] = {FC_DT_FLOAT, FC_DT_INT, FC_DT_DOUBLE};

static FC_Variable bad_var = {999,999};
static FC_Mesh bad_mesh;
static FC_Sequence bad_seq, empty_seq;
static float glob_decrtimes[5] = {0.1,1.0,0.5,2.0,3.0};
static float glob_nonincrtimes[5] = {0.1,1.0,1.0,2.0,3.0};
static FC_Sequence glob_decrseq; //sequence itself is not monotonically increasing
//but it is a function
static FC_Sequence glob_nonincrseq; //seq vars on this arent functions

//known vals based on formula for glob_sv1
//includes test of soemthing with 0 sdev
static double glob_sv1_msdev[6][2]= {{0.2,0},{2.2,1.581139},
		     {2.2,1.581139},{4.2,3.162278},
		     {4.2,3.162278},{6.2,4.743416}};



//helper function to change type of a seq var for testing various data types
static FC_ReturnCode ag_createChangedSeqDataTypeSeqVariable(int numStep,
							 FC_Variable *seqvar,
							 FC_DataType
							 newdatatype,
							 char * newname,
							 FC_Variable
							 **newseqvar){
  //no checks since helper function
  //know the sequences i have all have 5 steps
  FC_ReturnCode rc;
  FC_Mesh mesh;
  FC_AssociationType assoc;
  FC_Sequence seq, *newseq;
  FC_MathType mathtype;
  FC_DataType seqdatatype,seqvardatatype;
  int numDataPoint, numComponent;
  void* data;
  int i;

  if (numStep != 5){
    fc_printfErrorMessage("prechosen seq all only have 5 steps");
    return FC_ERROR;
  }
    

  rc = fc_getVariableInfo(seqvar[0],&numDataPoint,
			  &numComponent,
			  &assoc,&mathtype,NULL);
  if (rc != FC_SUCCESS){
    return rc;
  }


  rc = fc_getMeshFromVariable(seqvar[0],&mesh);
  if (rc != FC_SUCCESS){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(rc));
    return rc;

  }

  rc = fc_getSequenceFromSeqVariable(numStep,seqvar,&seq);
  if (rc != FC_SUCCESS)
    return rc;

  rc = fc_getSequenceDataType(seq,&seqdatatype);
  if (rc != FC_SUCCESS)
    return rc;

  if (seqdatatype == newdatatype){
    newseq = &seq;
  }else {
    switch(newdatatype){
    case FC_DT_INT:
      newseq = &glob_iseq;
      break;
    case FC_DT_DOUBLE:
      newseq = &glob_dseq;
      break;
    case FC_DT_FLOAT:
      newseq = &glob_dseq;
      break;
    case FC_DT_CHAR:
      newseq = &glob_cseq;
      break;
    default:
      //bad
      newseq = NULL;
      break;
    }
  }
      
  rc = fc_createSeqVariable(mesh, *newseq, newname, &numStep,newseqvar);
  if (rc != FC_SUCCESS){
    fc_deleteSeqVariable(numStep, *newseqvar); 
    free(*newseqvar);
    *newseqvar=NULL;
    return rc;
  }

  rc = fc_getVariableDataType(seqvar[0],&seqvardatatype);
  if (rc != FC_SUCCESS)
    return rc;

  for (i = 0; i < numStep; i++){
    rc = fc_getVariableDataPtr(seqvar[i], &data);
    if (rc != FC_SUCCESS){
      fc_deleteSeqVariable(numStep,*newseqvar);
      free(*newseqvar);
      *newseqvar=NULL;
      return rc;
    }

    rc = fc_setVariableData((*newseqvar)[i], numDataPoint,  numComponent,
			    assoc, mathtype, seqvardatatype, data);
    if (rc != FC_SUCCESS){
      fc_deleteSeqVariable(numStep,*newseqvar);
      free(*newseqvar);
      *newseqvar=NULL;
      free(data);
      return rc;
    }
  }

  return FC_SUCCESS;
}


//helper function to change type of a seq var for testing various data types
static FC_ReturnCode ag_createChangedDataTypeSeqVariable(int numStep,
							 FC_Variable *seqvar,
							 FC_DataType
							 newdatatype,
							 char * newname,
							 FC_Variable
							 **newseqvar){
  //no checks since helper function
  FC_ReturnCode rc;
  int numDataPoint, numComponent;
  FC_AssociationType assoc;
  FC_MathType mathtype;
  int i;
  void* data;


  rc = fc_getVariableInfo(seqvar[0],&numDataPoint,
			  &numComponent,
			  &assoc,&mathtype,NULL);
  if (rc != FC_SUCCESS){
    return rc;
  }


  rc = _fc_makeNewSeqVar(numStep,seqvar,newname,newseqvar);
  if (rc != FC_SUCCESS){
      fc_deleteSeqVariable(numStep,*newseqvar);
      free(*newseqvar);
      *newseqvar=NULL;
    return rc;
  }


  for (i = 0; i < numStep; i++){
    rc = fc_getVariableDataAsDataType(seqvar[i],newdatatype, &data);
    if (rc != FC_SUCCESS){
      fc_deleteSeqVariable(numStep,*newseqvar);
      free(*newseqvar);
      *newseqvar=NULL;
      return rc;
    }

    rc = fc_setVariableDataPtr((*newseqvar)[i], numDataPoint,  numComponent,
			    assoc, mathtype, newdatatype, data);
    if (rc != FC_SUCCESS){
      fc_deleteSeqVariable(numStep,*newseqvar);
      free(*newseqvar);
      *newseqvar=NULL;
      free(data);
      return rc;
    }
  }

  return FC_SUCCESS;
}


static void series_setup(void) {
  FC_ReturnCode rc;
  int seqStep;
  int numComponent;
  float* componentData;
  char* charData;
  int i,j,k;

  if (isForking) {
    fc_setLibraryVerbosity(fc_messages);
    rc = fc_initLibrary();
    fail_unless(rc == FC_SUCCESS, "test aborted: failed to init library");
  }

  //make a data set to play in
  rc = fc_createDataset("glob dataset", &glob_dataset);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create dataset");

  //make a mesh
  rc = fc_createMesh(glob_dataset," glob mesh", &glob_mesh);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create mesh");
  rc = fc_setMeshCoords(glob_mesh, glob_meshDim, glob_meshDataPoint,
			glob_meshcoords);
  fail_unless(rc == FC_SUCCESS, "failed to set vertex coords");


  //create float sequence - main global sequence
  rc = fc_createSequence(glob_dataset,"glob seq",&glob_fseq);
  fail_unless(rc == FC_SUCCESS, "failed to create sequence");
  rc = fc_setSequenceCoords(glob_fseq,glob_numStep,glob_seqDataType,
			    glob_ftimes);
  fail_unless(rc == FC_SUCCESS, "failed to set sequence coords");

  //some seq vars based on this sequence
  numComponent = 2;
  rc = fc_createSeqVariable(glob_mesh,glob_fseq,glob_sv1_name,&seqStep,
		       &glob_sv1);
  fail_unless(rc == FC_SUCCESS, "failed to create sv1");
  rc = fc_createSeqVariable(glob_mesh,glob_fseq,glob_svchar_name,&seqStep,
		       &glob_svchar);
  fail_unless(rc == FC_SUCCESS, "failed to create svchar");

  componentData = (float*)malloc(glob_meshDataPoint*
				 numComponent*sizeof(float));
  charData = (char*)malloc(glob_meshDataPoint*numComponent*sizeof(char));
  for (i  = 0; i < glob_numStep; i++){ //like each time
    for (j = 0; j < glob_meshDataPoint; j++) { //like each node
      for (k = 0; k < numComponent; k++){
	componentData[numComponent*j+k] = glob_ftimes[i]  * 
	  ((float)j+(float)k)+0.2;
	charData[numComponent*j+k] = 'a';
      }
    }
    rc = fc_setVariableData(glob_sv1[i],glob_meshDataPoint,numComponent,
		       FC_AT_VERTEX,FC_MT_VECTOR,
		       FC_DT_FLOAT,componentData);
    fail_unless(rc == FC_SUCCESS, "failed to set data for glob_sv1");
    rc = fc_setVariableData(glob_svchar[i],glob_meshDataPoint,numComponent,
			    FC_AT_VERTEX,FC_MT_VECTOR,
			    FC_DT_CHAR,charData);
    fail_unless(rc == FC_SUCCESS, "failed to set data for glob_svchar");
  }
  free(componentData);
  free(charData);


  // additional sequences. will use these with varying
  // coordinates, so will make variables on these specific to the tests
  rc = fc_createSequence(glob_dataset,"glob_cseq",&glob_cseq);
  fail_unless(rc == FC_SUCCESS, "failed to create sequence");
  rc = fc_setSequenceCoords(glob_cseq,glob_numStep,FC_DT_CHAR,glob_ctimes);
  fail_unless(rc == FC_SUCCESS, "failed to set sequence coords");

  rc = fc_createSequence(glob_dataset,"glob_iseq",&glob_iseq);
  fail_unless(rc == FC_SUCCESS, "failed to create sequence");
  rc = fc_setSequenceCoords(glob_iseq,glob_numStep,FC_DT_INT,glob_itimes);
  fail_unless(rc == FC_SUCCESS, "failed to set sequence coords");

  rc = fc_createSequence(glob_dataset,"glob_dseq",&glob_dseq);
  fail_unless(rc == FC_SUCCESS, "failed to create sequence");
  rc = fc_setSequenceCoords(glob_dseq,glob_numStep,FC_DT_DOUBLE,glob_dtimes);
  fail_unless(rc == FC_SUCCESS, "failed to set sequence coords");

  rc = fc_createSequence(glob_dataset,"glob_decrseq",&glob_decrseq);
  fail_unless(rc == FC_SUCCESS, "failed to create sequence");
  rc = fc_setSequenceCoords(glob_decrseq,glob_numStep,FC_DT_FLOAT,
			    glob_decrtimes);
  fail_unless(rc == FC_SUCCESS, "failed to set sequence coords");

  rc = fc_createSequence(glob_dataset,"glob_nonincrseq",&glob_nonincrseq);
  fail_unless(rc == FC_SUCCESS, "failed to create sequence");
  rc = fc_setSequenceCoords(glob_nonincrseq,glob_numStep,FC_DT_FLOAT,
			    glob_nonincrtimes);
  fail_unless(rc == FC_SUCCESS, "failed to set sequence coords");


  //and some bad things
  rc = fc_createMesh(glob_dataset,"bad_mesh", &bad_mesh);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create bad mesh");
  rc = fc_setMeshCoords(bad_mesh, glob_meshDim, glob_meshDataPoint,
			glob_meshcoords);
  fail_unless(rc == FC_SUCCESS, "failed to set vertex coords");
  //bad seq is just bad cause its a different seq with the same values
  rc = fc_createSequence(glob_dataset,"bad_seq",&bad_seq);
  fail_unless(rc == FC_SUCCESS, "failed to create bad sequence");
  rc = fc_setSequenceCoords(bad_seq,glob_numStep,glob_seqDataType,
			    glob_ftimes);
  fail_unless(rc == FC_SUCCESS, "failed to set bad sequence coords");
  rc = fc_createSequence(glob_dataset,"empty_seq",&empty_seq);
  fail_unless(rc == FC_SUCCESS, "failed to create empty sequence");

}

static void series_teardown(void) {
  FC_ReturnCode rc;

  rc = fc_deleteSeqVariable(glob_numStep,glob_sv1);
  fail_unless(rc == FC_SUCCESS, "failed to delete glob_sv1 at end of tests");
  free(glob_sv1);
  rc = fc_deleteSeqVariable(glob_numStep,glob_svchar);
  fail_unless(rc == FC_SUCCESS, "failed to delete glob_svchar at end of tests");
  free(glob_svchar);
  rc = fc_deleteSequence(glob_fseq);
  fail_unless(rc == FC_SUCCESS, "failed to delete seq at end of tests");
  rc = fc_deleteSequence(glob_iseq);
  fail_unless(rc == FC_SUCCESS, "failed to delete seq at end of tests");
  rc = fc_deleteSequence(glob_cseq);
  fail_unless(rc == FC_SUCCESS, "failed to delete seq at end of tests");
  rc = fc_deleteSequence(glob_dseq);
  fail_unless(rc == FC_SUCCESS, "failed to delete seq at end of tests");
  rc = fc_deleteSequence(glob_decrseq);
  fail_unless(rc == FC_SUCCESS, "failed to delete seq at end of tests");
  rc = fc_deleteSequence(glob_nonincrseq);
  fail_unless(rc == FC_SUCCESS, "failed to delete seq at end of tests");
  rc = fc_deleteSequence(bad_seq);
  fail_unless(rc == FC_SUCCESS, "failed to delete badseq at end of tests");
  rc = fc_deleteSequence(empty_seq);
  fail_unless(rc == FC_SUCCESS, "failed to delete emptyseq at end of tests");
  rc = fc_deleteMesh(glob_mesh);
  fail_unless(rc == FC_SUCCESS, "failed to delete mesh at end of tests");
  rc = fc_deleteMesh(bad_mesh);
  fail_unless(rc == FC_SUCCESS, "failed to delete badmesh at end of tests");
  rc = fc_deleteDataset(glob_dataset);
  fail_unless(rc == FC_SUCCESS, "failed to delete dataset at end of tests");

  if (isForking) {
    rc = fc_finalLibrary();
    fail_unless(rc == FC_SUCCESS, "failed to final library at end of test");
  }
}

START_TEST(do_regtoseq)
{
  //no longer checks reg seq since thats part of sequence now.
  //all this checks is regvartoseqvar and thats a private function
  //ill keep the test in since i use it a lot to make tests

  //this does *not* use global data set and mesh
  FC_ReturnCode rc;
  FC_Dataset dataset;
  FC_Mesh mesh;

  //seq
  FC_Sequence seq, checkSeq, *returnSequences;
  FC_Variable *seqVar, **checkSeqVars;
  int newNumComponent, newNumDataPoint;
  int checkStep, *checkSteps, numCheckVars;
  int numReturnSequences;
  FC_AssociationType newAssoc;
  FC_MathType newMathtype;
  FC_DataType newDatatype;
  void* coordsp;

  //normal
  FC_Variable normal;
  FC_AssociationType assoc = FC_AT_VERTEX;
  FC_MathType mathtype = FC_MT_VECTOR;

  //mesh
  int numDataPoint =5;
  int meshDim = 1;
  double meshcoords[5]={0.1,1.1,2.1,3.1,4.1};
  int numComponent = 3;


  float* fcomponentData;
  int* icomponentData;
  double* dcomponentData;
  char* ccomponentData;

  int i,j,k;

  FC_DataType dtype[4]= {FC_DT_FLOAT, FC_DT_INT, FC_DT_DOUBLE, FC_DT_CHAR};
  //dont test unknown cause cant make var of that type

  //setup
  //make a data set to play in
  rc = fc_createDataset("temp dataset", &dataset);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create dataset");

  //make a mesh
  rc = fc_createMesh(dataset,"temp_mesh", &mesh);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create mesh");
  rc = fc_setMeshCoords(mesh, meshDim, numDataPoint, meshcoords);
  fail_unless(rc == FC_SUCCESS, "failed to set vertex coords");


  //make data
  fcomponentData = (float*)malloc(numDataPoint*numComponent*sizeof(float));
  icomponentData = (int*)malloc(numDataPoint*numComponent*sizeof(int));
  dcomponentData = (double*)malloc(numDataPoint*numComponent*sizeof(double));
  ccomponentData = (char*)malloc(numDataPoint*numComponent*sizeof(char));

  for (i = 0; i < numComponent*numDataPoint; i++){
    icomponentData[i] = i;
    fcomponentData[i] = (float)(icomponentData[i]);
    dcomponentData[i] = (double)(icomponentData[i]);
    ccomponentData[i] = 'a';
  }


  //make a sequence
  rc = fc_createRegularSequence(dataset,numDataPoint,1.2,1.3,
				"seq_name",&seq);
  fail_unless(rc == FC_SUCCESS, "couldnt create regular sequnce");

  //testing vals for it here
  rc = fc_getSequenceNumStep(seq,&checkStep);
  fail_unless(numDataPoint == checkStep, "created wrong number of steps");
  rc = fc_getSequenceByName(dataset,"seq_name",&numReturnSequences,&returnSequences);
  fail_unless(rc == FC_SUCCESS, "unable to get created seq by name");
  fail_unless(numReturnSequences == 1, "unable to find unique matching sequence");
  fail_unless(FC_HANDLE_EQUIV(returnSequences[0],seq),"mismatch of seq handles");
  free(returnSequences);

  rc = fc_getSequenceCoordsPtr (seq, &coordsp);
  fail_unless(rc == FC_SUCCESS, "unable to get created seq var coords");
  for (j = 0; j < numDataPoint; j++){
    fail_unless( FC_DBL_EQUIV(((double*)coordsp)[j],1.2+j*1.3),
		 "bad seq coords");
  }

  //now test reg var to seq var
  for (i = 0; i < 4; i++){
    //make normal var
    rc = fc_createVariable(mesh,"junk",&normal);
    fail_unless(rc == FC_SUCCESS, "failed to create test variable");
    switch(dtype[i]){
    case FC_DT_INT:
      rc = fc_setVariableData(normal,numDataPoint,numComponent,
			      assoc,mathtype,
			      dtype[i],icomponentData);
      fail_unless(rc == FC_SUCCESS, "failed to set test data");
      break;
    case FC_DT_DOUBLE:
      rc = fc_setVariableData(normal,numDataPoint,numComponent,
			      assoc,mathtype,
			      dtype[i],dcomponentData);
      fail_unless(rc == FC_SUCCESS, "failed to set test data");
      break;
    case FC_DT_FLOAT:
      rc = fc_setVariableData(normal,numDataPoint,numComponent,
			      assoc,mathtype,
			      dtype[i],fcomponentData);
      fail_unless(rc == FC_SUCCESS, "failed to set test data");
      break;
    case FC_DT_CHAR:
      rc = fc_setVariableData(normal,numDataPoint,numComponent,
			      assoc,mathtype,
			      dtype[i],ccomponentData);
      fail_unless(rc == FC_SUCCESS, "failed to set test data");
    default:
      //wont happen
      break;
    }


    //make a seq var from it
    rc = _fc_createSeqVariableFromVariable(numDataPoint,seq,normal,
					  "seq_var_name",
					  &seqVar);
    fail_unless(rc == FC_SUCCESS, 
		"failed to create seq var from regular var");
    

    //checks seq var
    rc = fc_getVariableInfo(seqVar[0],&newNumDataPoint,
			    &newNumComponent,
			    &newAssoc,&newMathtype,&newDatatype);
    fail_unless(rc == FC_SUCCESS, "unable to get info on new seq var");
    
    fail_unless(newNumDataPoint == 1,"mismatch of numDataPoint");
    fail_unless(newNumComponent == numComponent,"mismatch of numComponent");
    fail_unless(newAssoc == FC_AT_WHOLE_MESH,"mismatch of assoc");
    fail_unless(newMathtype == mathtype,"mismatch of mathtype");
    fail_unless(newDatatype == dtype[i],"mismatch of datatype");

    rc = fc_getSequenceFromSeqVariable(numDataPoint,seqVar,&checkSeq);
    fail_unless(rc == FC_SUCCESS, "unable to get seq of new seq var");
    fail_unless(FC_HANDLE_EQUIV(checkSeq,seq),"mismatch of seq handles");

    rc = fc_getSeqVariableByName(mesh,"seq_var_name",&numCheckVars,
				 &checkSteps, &checkSeqVars);
    fail_unless(rc == FC_SUCCESS, "unable to get created seq var by name");
    fail_unless(numCheckVars == 1, "wrong number of matching vars");
    fail_unless(FC_HANDLE_EQUIV((*checkSeqVars[0]),(*seqVar)),
		"mismatch of seq handles");
    free(checkSeqVars[0]);
    free(checkSeqVars);
    free(checkSteps);

    rc = fc_getSeqVariableByName(mesh,"seq_var_name",&numCheckVars,
				 &checkSteps, &checkSeqVars);
    fail_unless(rc == FC_SUCCESS, "unable to get created seq var by name");
    fail_unless(numCheckVars == 1, "wrong number of matching vars");
    fail_unless(FC_HANDLE_EQUIV((*checkSeqVars[0]),(*seqVar)),
		"mismatch of seq handles");
    free(checkSeqVars[0]);
    free(checkSeqVars);
    free(checkSteps);

    //check values for seq var
    for (j = 0; j < numDataPoint; j++){
      void* data;
      rc = fc_getVariableDataPtr(seqVar[j],&data);
      fail_unless(rc == FC_SUCCESS, "cant get seq var data");
      for (k = 0; k <numComponent; k++){
	switch(dtype[i]){
	case FC_DT_FLOAT:
	  fail_unless(FC_FLT_EQUIV(((float*)data)[k],
				   fcomponentData[j*numComponent+k]),
		      "mismatch of value");
	  break;
	case FC_DT_INT:
	  fail_unless(((int*)data)[k] == 
		      icomponentData[j*numComponent+k],
		      "mismatch of value");
	  break;
	case FC_DT_DOUBLE:
	  fail_unless(FC_DBL_EQUIV(((double*)data)[k],
				   dcomponentData[j*numComponent+k]),
		      "mismatch of value");
	  break;
	case FC_DT_CHAR:
	  fail_unless(((char*)data)[k] == 'a',"mismatch of value");
	  break;
	default:
	  //wont happen
	  break;
	}
      }
    }

   fc_deleteSeqVariable(numDataPoint,seqVar);
   free(seqVar);

   if (i != 3){
     //keep it for below 
     fc_deleteVariable(normal);
   }
  }
  
  free(fcomponentData);
  free(dcomponentData);
  free(icomponentData);
  free(ccomponentData);


  // --- test error conditions for seq var (bad args)
  rc = _fc_createSeqVariableFromVariable(0,seq,normal,
					      "seq_var_name",
					      &seqVar);
  fail_unless(rc != FC_SUCCESS, "should fail if bad input variable");
  fail_unless(&seqVar != NULL, "should return NULL seqvar when fails");

  rc = _fc_createSeqVariableFromVariable(numDataPoint,FC_NULL_SEQUENCE,normal,
					      "seq_var_name",
					      &seqVar);
  fail_unless(rc != FC_SUCCESS, "should fail if NULL input seq ");
  fail_unless(&seqVar != NULL, "should return NULL seqvar when fails");

  rc = _fc_createSeqVariableFromVariable(numDataPoint,seq,FC_NULL_VARIABLE,
					      "seq_var_name",
					      &seqVar);
  fail_unless(rc != FC_SUCCESS, "should fail if NULL input normal var ");
  fail_unless(&seqVar != NULL, "should return NULL seqvar when fails");

  rc = _fc_createSeqVariableFromVariable(numDataPoint,seq,bad_var,
					      "seq_var_name",
					      &seqVar);
  fail_unless(rc != FC_SUCCESS, "should fail if bad input normal var ");
  fail_unless(&seqVar != NULL, "should return NULL seqvar when fails");

  rc = _fc_createSeqVariableFromVariable(numDataPoint,seq,normal,
					NULL, &seqVar);
  fail_unless(rc != FC_SUCCESS, "should fail if NULL input seq var name ");
  fail_unless(&seqVar != NULL, "should return NULL seqvar when fails");

  rc = _fc_createSeqVariableFromVariable(numDataPoint,seq,normal,
					"seq_var_name", NULL);
  fail_unless(rc != FC_SUCCESS, "should fail if NULL input seq var");
  fail_unless(&seqVar != NULL, "should return NULL seqvar when fails");


  fc_deleteVariable(normal);
  fc_deleteSequence(seq);
  fc_deleteSeqVariable(numDataPoint,seqVar);
  free(seqVar);


  rc = fc_deleteMesh(mesh);
  fail_unless(rc == FC_SUCCESS, "failed to delete mesh at end of tests");
  rc = fc_deleteDataset(dataset);
  fail_unless(rc == FC_SUCCESS, "failed to delete dataset at end of tests");
}
END_TEST;



START_TEST(do_interpolate)
{
  //to messy to put in with other single things

  FC_ReturnCode rc;

  //seqvar
  int numComponent;
  int numDataPoint;
  int numStep;
  FC_AssociationType assoc;
  FC_MathType mathtype;
  FC_DataType datatype;

  FC_Sequence newSeq;
  int newNumStep;
  FC_DataType newSeqDatatype;

  //compare against
  FC_Variable *newSeqVar, *checkSeqVar, *convertSeqVar, **checkSeqVars;
  FC_AssociationType newAssoc;
  FC_MathType newMathtype;
  FC_DataType newDatatype;
  int newNumDataPoint, newNumComponent;
  int *checkSteps, numCheckVars;
  void *data; 
  void* newSeqCoords;
  
  int i,j,k,jj,kk;

  //use global dataset,mesh,seq, sv1
  //uses known vals from formula based on sv1 below

  //this first block iterates over new seq types and diff data types of the seq
  //var. it also tests some permutations of the ranges.
  //different seq types in the seq variable are handled separately below
  for (kk = 0; kk < 3; kk++){ //iterate over new seq types
    newSeqDatatype = glob_dtype[kk];
    switch(newSeqDatatype){
    case FC_DT_DOUBLE:
      //first lets get a new regular sequnce thats offset from glob_fseq
      //this seq ends with the same last point as the series we are
      // comparing it with and overlap with some of the points in them
      newNumStep = 7;
      rc = fc_createRegularSequence(glob_dataset,newNumStep,1.6,0.4,
				    "seq_name",&newSeq);
      fail_unless(rc == FC_SUCCESS, "couldnt create regular sequence");
      break;
    case FC_DT_INT:{
      //first lets get a new regular sequence thats offset from glob_fseq
      //this seq starts with the same first point as the series we are
      // comparing it with and overlaps with some of the points in them
      int localtimes[3] = {0,1,3};
      newNumStep = 3;
      rc = fc_createSequence(glob_dataset,"junkseq",&newSeq);
      fail_unless(rc == FC_SUCCESS, "failed to create sequence");
      rc = fc_setSequenceCoords(newSeq,newNumStep,FC_DT_INT,localtimes);
      fail_unless(rc == FC_SUCCESS, "failed to set sequence coords");
      break;
    }
    case FC_DT_FLOAT:{
      //first lets get a new regular sequence thats offset from glob_fseq
      //this seq starts with the same first point as the series we are
      // comparing it with and overlaps with some of the points in them
      float localtimes[5] = {0.0,1.2,1.8,3.0,3.2};
      newNumStep = 5;
      rc = fc_createSequence(glob_dataset,"junkseq",&newSeq);
      fail_unless(rc == FC_SUCCESS, "failed to create sequence");
      rc = fc_setSequenceCoords(newSeq,newNumStep,FC_DT_FLOAT,localtimes);
      fail_unless(rc == FC_SUCCESS, "failed to set sequence coords");
      break;
    }
    default:
      //wont happen
      break;
    }

    for (jj = 0; jj < 3; jj++){ //iterate over variable data types with
      //float seq
      numStep = glob_numStep;
      rc = ag_createChangedDataTypeSeqVariable(numStep,glob_sv1,
					       glob_dtype[jj],
					       "convert",&convertSeqVar);
      fail_unless(rc == FC_SUCCESS, "couldnt convert seqvar");
      
      rc = fc_getVariableInfo(convertSeqVar[0],&numDataPoint,&numComponent,
			      &assoc,&mathtype,&datatype);
      fail_unless(rc == FC_SUCCESS, "couldnt get seqvar info");

      /*
      printf("orig:\n");
      for (i = 0; i < numStep; i++){
	printf("step %10g: ",((float*)glob_ftimes)[i]);
	rc = fc_getVariableDataPtr(convertSeqVar[i],&data);
	fail_unless(rc == FC_SUCCESS, "can't get convertSeqVar data");
	for (j = 0; j < numDataPoint; j++){
	  for (k = 0; k < numComponent; k++){
	    switch(glob_dtype[jj]){
	    case FC_DT_FLOAT:
	      printf("%10g ",((float*)data)[numComponent*j+k]);
	      break;
	    case FC_DT_INT:
	      printf("%d ",((int*)data)[numComponent*j+k]);
	      break;
	    case FC_DT_DOUBLE:
	      printf("%10g ",((double*)data)[numComponent*j+k]);
	      break;
	    default:
	      break;
	    }
	  }
	}
	printf("\n");
      }
      */

      //linearInterpolation test
      rc = fc_linearInterpolation(numStep,convertSeqVar,newSeq,"junk",
				  &newSeqVar);
      fail_unless(rc == FC_SUCCESS, "failed to do linearInterpolation");
      fail_unless(fc_isSeqVariableValid(newNumStep,newSeqVar),
		  "didn't create valid newvar with correct number of steps");
      fc_getVariableInfo(newSeqVar[0],&newNumDataPoint,&newNumComponent,
			 &newAssoc,&newMathtype,&newDatatype);
      fail_unless(newNumDataPoint == numDataPoint,"mismatch of numDataPoint");
      fail_unless(newNumComponent == numComponent,"mismatch of numComponent");
      fail_unless(newAssoc == assoc,"mismatch of assoc");
      fail_unless(newMathtype == mathtype,"mismatch of mathtype");
      fail_unless(newDatatype == FC_DT_DOUBLE,"mismatch of datatype");

      rc = fc_getSeqVariableByName(glob_mesh,"junk",&numCheckVars,
				   &checkSteps,&checkSeqVars);
      fail_unless(rc == FC_SUCCESS, "failed to get interpolated var by name");
      fail_unless(numCheckVars == 1, "wrongnumber of matching vars");
      fail_unless(FC_HANDLE_EQUIV(newSeqVar[0],checkSeqVars[0][0]),
		  "mismatch of variable handles");
      free(checkSeqVars[0]);
      free(checkSeqVars);
      free(checkSteps);

      //validity test 
      //      printf("new: \n");

      rc = fc_getSequenceCoordsPtr(newSeq, &newSeqCoords);
      fail_unless(rc == FC_SUCCESS, "couldnt get seq coords");
      for (i = 0; i < newNumStep; i++){
	double val;
	switch (newSeqDatatype){
	case FC_DT_INT:
	  val =  (double)(((int*)newSeqCoords)[i]);
	  break;
	case FC_DT_FLOAT:
	  val =  (double)(((float*)newSeqCoords)[i]);
	  break;
	case FC_DT_DOUBLE:
	  val =  ((double*)newSeqCoords)[i];
	  break;
	default:
	  //wont happen
	  break;
	}

	//	printf("step %5g: ",val);
	rc = fc_getVariableDataPtr(newSeqVar[i],&data);
	fail_unless(rc == FC_SUCCESS, "can't get newSeqVar data");
	for (j = 0; j < newNumDataPoint; j++){
	  for (k = 0; k < newNumComponent; k++){
	    //known values
	    double temp;
	    if (glob_dtype[jj] == FC_DT_INT){
	      temp = val*(j+k);
	    }else {
	      temp = val*(j+k)+0.2;
	    }
	    //	    printf("(%5g, %5g) ",
	    //		   ((double*)data)[newNumComponent*j+k],temp);
	    fail_unless(FC_FLT_EQUIV(((double*)data)[newNumComponent*j+k],
				     temp),
			"bad match for interpolation value");
	  }
	}
	//	printf("\n");
      }

      //clean up for this test
      rc = fc_deleteSeqVariable(numStep,convertSeqVar);
      fail_unless(rc == FC_SUCCESS, "failed to delete seq var");
      free(convertSeqVar);
      rc = fc_deleteSeqVariable(newNumStep,newSeqVar);
      fail_unless(rc == FC_SUCCESS, "failed to delete seq var");
      free(newSeqVar);
    }
    rc = fc_deleteSequence(newSeq);
    fail_unless(rc == FC_SUCCESS, "failed to delete seq");
  }

  //  printf("\n");


  //now testing works ok for diff types of seq in the seq variables,
  //double and int.its easier to stick these separately rather than
  //put them in the loop above. float has already been checked above.
  //i dont think its necessary to check all possible permutations for this. 
  {
    //do a double and an int version in the sequence with the variable values
    //of glob_sv1, since we already know the results for that one
    //glob_dseq and glob_iseq have the same values.
    
    //comparison seq used for both 
    newNumStep = 7;
    rc = fc_createRegularSequence(glob_dataset,newNumStep,1.6,0.4,
				  "seq_name",&newSeq);
    fail_unless(rc == FC_SUCCESS, "couldnt create regular sequence");

    for (jj = 0; jj <2 ;jj++){// 0 = double, 1 = int
      numStep = glob_numStep;
      if (jj == 0){ //double seq
	rc = fc_createSeqVariable(glob_mesh,glob_dseq,"junk",&numStep,
				  &convertSeqVar);
	fail_unless(rc == FC_SUCCESS, "failed to create variable");
      }else { // int seq
	rc = fc_createSeqVariable(glob_mesh,glob_iseq,"junk",&numStep,
				  &convertSeqVar);
	fail_unless(rc == FC_SUCCESS, "failed to create variable");
      }
      
      numDataPoint = glob_meshDataPoint;
      for (i  = 0; i < numStep; i++){ //copy glob_sv1 data over
	void *data1; 
	rc = fc_getVariableDataPtr(glob_sv1[i],&data1); //glob_sv1 data is float
	fail_unless(rc == FC_SUCCESS, "failed to get variable data");
	rc = fc_setVariableData(convertSeqVar[i],numDataPoint,2,
				FC_AT_VERTEX,FC_MT_VECTOR,
				FC_DT_FLOAT,data1);
	fail_unless(rc == FC_SUCCESS, "failed to set variable data");
      }
      
      //linearInterpolation test
      rc = fc_linearInterpolation(numStep,convertSeqVar,newSeq,"junk_int",
				  &newSeqVar);
      fail_unless(rc == FC_SUCCESS, "failed to do linearInterpolation");
      fail_unless(fc_isSeqVariableValid(newNumStep,newSeqVar),
		  "didn't create valid newvar with correct number of steps");
      fc_getVariableInfo(newSeqVar[0],&newNumDataPoint,&newNumComponent,
			 &newAssoc,&newMathtype,&newDatatype);
      fail_unless(newNumDataPoint == numDataPoint,"mismatch of numDataPoint");
      fail_unless(newNumComponent == numComponent,"mismatch of numComponent");
      fail_unless(newAssoc == assoc,"mismatch of assoc");
      fail_unless(newMathtype == mathtype,"mismatch of mathtype");
      fail_unless(newDatatype == FC_DT_DOUBLE,"mismatch of datatype");
      
      rc = fc_getSeqVariableByName(glob_mesh,"junk_int",&numCheckVars,
				   &checkSteps,&checkSeqVars);
      fail_unless(rc == FC_SUCCESS, "failed to get interpolated var by name");
      fail_unless(numCheckVars == 1 , "wrong number of matching vars");
      fail_unless(FC_HANDLE_EQUIV(newSeqVar[0],checkSeqVars[0][0]),
				  "mismatch of variable handles");
      free(checkSeqVars[0]);
      free(checkSeqVars);
      free(checkSteps);
      
      //validity test 
      rc = fc_getSequenceCoordsPtr(newSeq, &newSeqCoords);
      fail_unless(rc == FC_SUCCESS, "couldnt get seq coords");
      
      for (i = 0; i < newNumStep; i++){
	double val = ((double*)newSeqCoords)[i];  //i know coords are double
	rc = fc_getVariableDataPtr(newSeqVar[i],&data);
	fail_unless(rc == FC_SUCCESS, "can't get newSeqVar data");
	for (j = 0; j < newNumDataPoint; j++){
	  for (k = 0; k < newNumComponent; k++){
	    //known values
	    double temp;
	    temp = val*(j+k)+0.2;
	    fail_unless(FC_FLT_EQUIV(((double*)data)[newNumComponent*j+k],
				     temp),"bad match for interpolation value");
	  }
	}
      }
      
      //clean up - keep comparison seq
      rc = fc_deleteSeqVariable(numStep,convertSeqVar);
      fail_unless(rc == FC_SUCCESS, "failed to delete seq var");
      free(convertSeqVar);
      rc = fc_deleteSeqVariable(newNumStep,newSeqVar);
      fail_unless(rc == FC_SUCCESS, "failed to delete seq var");
      free(newSeqVar);
    }

    //clean up
    rc = fc_deleteSequence(newSeq);
    fail_unless(rc == FC_SUCCESS, "failed to delete seq");
  }



  // --- test error conditions (bad args)
  rc = fc_linearInterpolation(999,glob_sv1,glob_fseq,"junk",&newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail if bad numstep");
  fail_unless(newSeqVar == NULL, "failure should return NULL seq var");
  rc = fc_linearInterpolation(glob_numStep,&bad_var,glob_fseq,
			      "junk",&newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail if bad input seq variable");
  fail_unless(newSeqVar == NULL, "failure should return NULL seq var");
  rc = fc_linearInterpolation(glob_numStep,NULL,glob_fseq,"junk",&newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail if null input seq variable");
  fail_unless(newSeqVar == NULL, "failure should return NULL seq var");
  rc = fc_linearInterpolation(numStep,glob_sv1,empty_seq,"junk",&newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail if empty input seq");
  fail_unless(newSeqVar == NULL, "failure should return NULL seq var");
  rc = fc_linearInterpolation(glob_numStep,glob_sv1,glob_fseq,
			      NULL,&newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail if NULL name");
  fail_unless(newSeqVar == NULL, "failure should return NULL seq var");
  rc = fc_linearInterpolation(glob_numStep,glob_sv1,glob_fseq,"junk",NULL);
  fail_unless(rc != FC_SUCCESS, "should fail if null return seq variable");

  //now check other failures 
  {
    int seqStep;
    FC_Variable *junk_sv;


    //need a good new seq again
    rc = fc_createRegularSequence(glob_dataset,7,1.6,0.4,
				  "seq_name",&newSeq);
    fail_unless(rc == FC_SUCCESS, "couldnt create regular sequence");

    //check failure for char variable data
    rc = fc_linearInterpolation(glob_numStep,glob_svchar,newSeq,"junk",
				&newSeqVar);
    fail_unless(rc != FC_SUCCESS,
		"should fail for seq var type == FC_DT_CHAR");
    fail_unless(newSeqVar == NULL, "failure should return NULL seq var");
    //no clean up for this test

    //check failure for char sequence in the variable
    //create seq var
    rc = fc_createSeqVariable(glob_mesh,glob_cseq,"junksv",&seqStep,
			      &junk_sv);
    fail_unless(rc == FC_SUCCESS, "failed to create sv");

    for (i  = 0; i < glob_numStep; i++){ //like each time
      rc = fc_setVariableData(junk_sv[i],1,1,FC_AT_WHOLE_MESH,FC_MT_SCALAR,
			      FC_DT_FLOAT,&(glob_ftimes[i]));
      fail_unless(rc == FC_SUCCESS, "failed to set data for junk_sv");
    }

    //check failure 
    rc = fc_linearInterpolation(glob_numStep,junk_sv,newSeq,
				"junk",&newSeqVar);
    fail_unless(rc != FC_SUCCESS,
		"should fail for seqdatatype == FC_DT_CHAR");
    fail_unless(newSeqVar == NULL, "failure should return NULL seq var");
    //clean up for this test
    rc = fc_deleteSeqVariable(glob_numStep,junk_sv);
    fail_unless(rc == FC_SUCCESS, "failed to delete junk_sv");
    free(junk_sv);


    //check failure for char sequence in the new sequence
    rc = fc_linearInterpolation(glob_numStep,glob_sv1,glob_cseq,"junk",
				&newSeqVar);
    fail_unless(rc != FC_SUCCESS,
		"should fail for seqdatatype == FC_DT_CHAR");
    fail_unless(newSeqVar == NULL, "failure should return NULL seq var");
    //no clean up for this test


    //check failure for decreasing sequence in variable
    //create seq var
    numStep = glob_numStep;
    rc = fc_createSeqVariable(glob_mesh,glob_decrseq,"junksv",&seqStep,
			      &junk_sv);
    fail_unless(rc == FC_SUCCESS, "failed to create sv");

    for (i  = 0; i < numStep; i++){ //like each time
      rc = fc_setVariableData(junk_sv[i],1,1,FC_AT_WHOLE_MESH,FC_MT_SCALAR,
			      FC_DT_FLOAT,&(glob_ftimes[i]));
      fail_unless(rc == FC_SUCCESS, "failed to set data for junk_sv");
    }

    //check failure 
    rc = fc_linearInterpolation(glob_numStep,junk_sv,newSeq,
				"junk",&newSeqVar);
    fail_unless(rc != FC_SUCCESS, "should fail for decreasing sequence");
    fail_unless(newSeqVar == NULL, "failure should return NULL seq var");
    //clean up for this test
    rc = fc_deleteSeqVariable(numStep,junk_sv);
    fail_unless(rc == FC_SUCCESS, "failed to delete junk_sv");
    free(junk_sv);

    //check failure for decreasing sequence in new seq
    rc = fc_linearInterpolation(glob_numStep,glob_sv1,glob_decrseq,
				"junk",&newSeqVar);
    fail_unless(rc != FC_SUCCESS, "should fail for decreasing sequence");
    fail_unless(newSeqVar == NULL, "failure should return NULL seq var");
    //no clean up for this test


    //check failure for unchanging sequence in the variable
    //create seq var
    numStep = glob_numStep;
    rc = fc_createSeqVariable(glob_mesh,glob_nonincrseq,"junksv",&seqStep,
			      &junk_sv);
    fail_unless(rc == FC_SUCCESS, "failed to create sv");

    for (i  = 0; i < numStep; i++){ //like each time
      rc = fc_setVariableData(junk_sv[i],1,1,FC_AT_WHOLE_MESH,FC_MT_SCALAR,
			      FC_DT_FLOAT,glob_ftimes);
      fail_unless(rc == FC_SUCCESS, "failed to set data for junk_sv");
    }

    //check failure 
    rc = fc_linearInterpolation(glob_numStep,junk_sv,newSeq,
				"junk",&newSeqVar);
    fail_unless(rc != FC_SUCCESS, "should fail for unchanging sequence");
    fail_unless(newSeqVar == NULL, "failure should return NULL seq var");
    //clean up for this test
    rc = fc_deleteSeqVariable(numStep,junk_sv);
    fail_unless(rc == FC_SUCCESS, "failed to delete junk_sv");
    free(junk_sv);

    //check failure for unchanging sequence in the new seq
    rc = fc_linearInterpolation(glob_numStep,glob_sv1,glob_nonincrseq,
				"junk",&newSeqVar);
    fail_unless(rc != FC_SUCCESS, "should fail for unchanging sequence");
    fail_unless(newSeqVar == NULL, "failure should return NULL seq var");
    //no clean up for this test

    //new seq bad boundaries
    {
      int localnumStep = 3;
      float timesa[3] = {-1.0,0.0,1.0};
      float timesb[3] = {1.0,2.0,10.0};
      FC_Sequence tempseq;

      rc = fc_createSequence(glob_dataset,"junk",&tempseq);
      fail_unless(rc == FC_SUCCESS, "failed to create sequence");
      rc = fc_setSequenceCoords(tempseq,localnumStep,FC_DT_FLOAT,
				timesa);
      fail_unless(rc == FC_SUCCESS, "failed to set sequence coords");

      //check failure for new seq begin too soon
      rc = fc_linearInterpolation(glob_numStep,glob_sv1,tempseq,
				"junk",&newSeqVar);
      fail_unless(rc != FC_SUCCESS, "should fail for new seq begin too soon");
      fail_unless(newSeqVar == NULL, "failure should return NULL seq var");
      //clean up for this test
      fc_deleteSequence(tempseq);

      rc = fc_createSequence(glob_dataset,"junk",&tempseq);
      fail_unless(rc == FC_SUCCESS, "failed to create sequence");
      rc = fc_setSequenceCoords(tempseq,localnumStep,FC_DT_FLOAT,
				timesb);
      fail_unless(rc == FC_SUCCESS, "failed to set sequence coords");

      //check failure for new seq end too late
      rc = fc_linearInterpolation(glob_numStep,glob_sv1,tempseq,
				"junk",&newSeqVar);
      fail_unless(rc != FC_SUCCESS, "should fail for new seq end too late");
      fail_unless(newSeqVar == NULL, "failure should return NULL seq var");
      //clean up for this test
      fc_deleteSequence(tempseq);
    }

    rc= fc_deleteSequence(newSeq);
  }

  //not checking for other mathtypes becuase all it does is copy that
  //field over when calling fc_setVariableData internally. i have
  //already checked that it does a copy correctly.

  //clean up
  newSeqVar = NULL;
  checkSeqVar = NULL;
  convertSeqVar = NULL;

}
END_TEST


START_TEST(do_window)
{
  //too messy to put non window things in here
  //this tests window ave
  FC_ReturnCode rc;

  //seqvar
  int numComponent;
  int numDataPoint;
  int numStep;
  FC_AssociationType assoc;
  FC_MathType mathtype;
  FC_DataType datatype;

  //compare against
  FC_Variable *newSeqVar, *checkSeqVar, *convertSeqVar, *origSeqVar, **checkSeqVars;
  FC_AssociationType newAssoc;
  FC_MathType newMathtype;
  FC_DataType newDatatype;
  FC_Sequence newSeq;
  int newNumDataPoint, newNumComponent;
  int checkStep, *checkSteps, numCheckVars;
  void *data, *data2;
  float *componentData;
  
  int window;
  double timewindow;
  int i,j,k,jj;

  //use global dataset,mesh,seq
  //want a new seq var that isnt linear to distinguish the centered and
  //leading cases
  numComponent = 2;
  rc = fc_createSeqVariable(glob_mesh,glob_fseq,"junk",&checkStep,
		       &origSeqVar);
  fail_unless(rc == FC_SUCCESS, "failed to create origSeqVar");
  componentData = (float*)malloc(glob_meshDataPoint*
				 numComponent*sizeof(float));
  for (i  = 0; i < glob_numStep; i++){ //like each time
    for (j = 0; j < glob_meshDataPoint; j++) { //like each node
      for (k = 0; k < numComponent; k++){
	componentData[numComponent*j+k] =
	  ((float)j+(float)k)*glob_ftimes[i]*glob_ftimes[i] +0.2;
      }
    }
    rc = fc_setVariableData(origSeqVar[i],glob_meshDataPoint,numComponent,
			    FC_AT_VERTEX,FC_MT_VECTOR,
			    FC_DT_FLOAT,componentData);
    fail_unless(rc == FC_SUCCESS, "failed to set data for origSeqVar");
  }
  free(componentData);



  for (jj = 0; jj < 3; jj++){
    FC_Variable *tempSeqVar;
    FC_Sequence lseq;

    numStep = glob_numStep;
    //change sequence type
    rc = ag_createChangedSeqDataTypeSeqVariable(numStep,origSeqVar,glob_dtype[jj],
					     "timeconvert",&tempSeqVar);
    fail_unless(rc == FC_SUCCESS, "couldnt convert seqvar");
    //change variable type
    rc = ag_createChangedDataTypeSeqVariable(numStep,tempSeqVar,glob_dtype[jj],
					     "convert",&convertSeqVar);
    fail_unless(rc == FC_SUCCESS, "couldnt convert seqvar");
    fc_deleteSeqVariable(numStep,tempSeqVar);
    free(tempSeqVar);
    rc = fc_getSequenceFromSeqVariable(numStep,convertSeqVar,&lseq);
    fail_unless(rc == FC_SUCCESS, "couldnt get sequence seqvar");

    rc = fc_getVariableInfo(convertSeqVar[0],&numDataPoint,&numComponent,
		       &assoc,&mathtype,&datatype);
    fail_unless(rc == FC_SUCCESS, "couldnt get seqvar info");

    //average test
    for (window =1; window < 4; window++){
      //only do window sequence for float, for others only window = 3
      if (datatype != FC_DT_FLOAT && window != 3){
	continue;
      }

      rc = fc_leadingWindowAverage(numStep,convertSeqVar,
				   window,"window",&newSeqVar);
      fail_unless(rc == FC_SUCCESS, "failed to do window average");
      fail_unless(fc_isSeqVariableValid(numStep,newSeqVar),
		  "didn't create valid seqvar");
      fc_getVariableInfo(newSeqVar[0],&newNumDataPoint,&newNumComponent,
			 &newAssoc,&newMathtype,&newDatatype);
      fail_unless(newNumDataPoint == numDataPoint,"mismatch of numDataPoint");
      fail_unless(newNumComponent == numComponent,"mismatch of numComponent");
      fail_unless(newAssoc == assoc,"mismatch of assoc");
      fail_unless(newMathtype == mathtype,"mismatch of mathtype");
      fail_unless(newDatatype == (datatype == FC_DT_FLOAT ?
				  FC_DT_FLOAT: FC_DT_DOUBLE),
		  "mismatch of datatype");
      fc_getSequenceFromSeqVariable(numStep,newSeqVar,&newSeq);
      fail_unless(FC_HANDLE_EQUIV(newSeq,lseq),
		  "mismatch of sequence handles");
      fc_getSeqVariableByName(glob_mesh,"window",&numCheckVars,&checkSteps,&checkSeqVars);
      fail_unless(numCheckVars == 1, "wrong number of matching vars");
      for (i = 0; i < numStep; i++){
	fail_unless(FC_HANDLE_EQUIV(newSeqVar[i],checkSeqVars[0][i]),
		    "mismatch of variable handles");
      }
      free(checkSeqVars[0]);
      free(checkSeqVars);
      free(checkSteps);

      //test correct vals and validity of results here
      for (i = 0; i < numStep; i++){
	rc = fc_getVariableDataPtr(convertSeqVar[i],&data);
	fail_unless(rc == FC_SUCCESS, "can't get convertSeqVar data");
	rc = fc_getVariableDataPtr(newSeqVar[i],&data2);
	fail_unless(rc == FC_SUCCESS, "can't get newSeqVar data");
	for (j = 0; j < numDataPoint; j++){
	  for (k = 0; k < numComponent; k++){
	    if (glob_dtype[jj] == FC_DT_FLOAT){
	      float currval = ((float*)data2)[numComponent*j+k];
	      float origval = ((float*)data)[numComponent*j+k];
	      //known vals from formula:
	      if (window == 1 || i == 0){ 
		fail_unless(FC_VALUE_EQUIV(currval,origval,
					   FLT_EPSILON,FLT_MIN),
			    "bad ave value for window");
	      } else if (i == 1 || window == 2 ) {
		fail_unless(FC_VALUE_EQUIV(currval,
					   origval+((j+k)*(float)(-2*glob_ftimes[i]+1))/2.0,
					   FLT_EPSILON,FLT_MIN),
			    "bad ave value for window");
	      } else{
		fail_unless(FC_VALUE_EQUIV(currval,
					   origval+((j+k)*(float)(-6*glob_ftimes[i]+5))/3.0,
					   FLT_EPSILON,FLT_MIN),
			    "bad ave value for window");
	      }
	    }else {
	      double currval = ((double*)data2)[numComponent*j+k];	      
	      double origval;
	      if (glob_dtype[jj] == FC_DT_INT){
		origval = (double)(((int*)data)[numComponent*j+k]);
	      }else {
		origval = ((double*)data)[numComponent*j+k];
	      }
 	      //known vals from formula. note window = 3 for these
	      //using FLT_VALS since cast from float
	      switch(i){
	      case 0:
		fail_unless(FC_VALUE_EQUIV(currval,origval,
					   FLT_EPSILON,FLT_MIN),
			    "bad ave value for window");
		break;
	      case 1:
		fail_unless(FC_VALUE_EQUIV(currval,
					   origval+(j+k)*(float)(-2*glob_ftimes[i]+1)/2.0,
					   FLT_EPSILON,FLT_MIN),
			    "bad ave value for window");
		break;
	      default:
		fail_unless(FC_VALUE_EQUIV(currval,
					   origval+(j+k)*(float)(-6*glob_ftimes[i]+5)/3.0,
					   FLT_EPSILON,FLT_MIN),
			    "bad ave value for window");
		break;
	      }
	    }
	  }
	}
	data = NULL;
	data2 = NULL;
      }
      fc_deleteSeqVariable(numStep,newSeqVar);
      free(newSeqVar);
    } // end window


    for (window = 0; window < 5; window++){
      timewindow = (double)(window)/2.0;
      //only do window sequence for float,
      //for others only time window = 1.5 (window  = 3)
      if (datatype != FC_DT_FLOAT && window != 3){
	continue;
      }

      rc = fc_leadingWindowAverage_Time(numStep,convertSeqVar,
				   timewindow,"window",&newSeqVar);
      fail_unless(rc == FC_SUCCESS, "failed to do window average");
      fail_unless(fc_isSeqVariableValid(numStep,newSeqVar),
		  "didn't create valid seqvar");
      fc_getVariableInfo(newSeqVar[0],&newNumDataPoint,&newNumComponent,
			 &newAssoc,&newMathtype,&newDatatype);
      fail_unless(newNumDataPoint == numDataPoint,"mismatch of numDataPoint");
      fail_unless(newNumComponent == numComponent,"mismatch of numComponent");
      fail_unless(newAssoc == assoc,"mismatch of assoc");
      fail_unless(newMathtype == mathtype,"mismatch of mathtype");
      fail_unless(newDatatype == (datatype == FC_DT_FLOAT ?
				  FC_DT_FLOAT: FC_DT_DOUBLE),
		  "mismatch of datatype");
      fc_getSequenceFromSeqVariable(numStep,newSeqVar,&newSeq);
      fail_unless(FC_HANDLE_EQUIV(newSeq,lseq),
		  "mismatch of sequence handles");
      fc_getSeqVariableByName(glob_mesh,"window",&numCheckVars,&checkSteps,&checkSeqVars);
      fail_unless(numCheckVars == 1, "wrong number of matching vars");
      for (i = 0; i < numStep; i++){
	fail_unless(FC_HANDLE_EQUIV(newSeqVar[i],checkSeqVars[0][i]),
		    "mismatch of variable handles");
      }
      free(checkSeqVars[0]);
      free(checkSeqVars);
      free(checkSteps);
      
      //test correct vals and validity of results here
      for (i = 0; i < numStep; i++){
	rc = fc_getVariableDataPtr(convertSeqVar[i],&data);
	fail_unless(rc == FC_SUCCESS, "can't get convertSeqVar data");
	rc = fc_getVariableDataPtr(newSeqVar[i],&data2);
	fail_unless(rc == FC_SUCCESS, "can't get newSeqVar data");
	for (j = 0; j < numDataPoint; j++){
	  for (k = 0; k < numComponent; k++){
	    switch(glob_dtype[jj]){
	    case FC_DT_FLOAT:
	      {
		float currval = ((float*)data2)[numComponent*j+k];
		float origval = ((float*)data)[numComponent*j+k];
		//known vals from formula:
		if (timewindow < 1 || i == 0){
		  fail_unless(FC_VALUE_EQUIV(currval,origval,
					     FLT_EPSILON,FLT_MIN),
			      "bad ave value for window");
		} else if (timewindow < 2 || i == 1){
		  fail_unless(FC_VALUE_EQUIV(currval,
					     origval+((j+k)*(float)(-2*glob_ftimes[i]+1))/2.0,
					     FLT_EPSILON,FLT_MIN),
			      "bad ave value for window");
		} else{ //window == 2 && i >= 2
		  fail_unless(FC_VALUE_EQUIV(currval,
					     origval+(j+k)*(float)(-6*glob_ftimes[i]+5)/3.0,
					     FLT_EPSILON,FLT_MIN),
			      "bad ave value for window");
		}
	      }
	      break;
	    default: 
	      {
		double currval = ((double*)data2)[numComponent*j+k];	      
		double origval;
		if (glob_dtype[jj] == FC_DT_INT){
		  origval = (double)(((int*)data)[numComponent*j+k]);
		}else {
		  origval = ((double*)data)[numComponent*j+k];
		}
		//known vals from formula: window = 3
		if (timewindow < 1 || i == 0){
		  fail_unless(FC_VALUE_EQUIV(currval,origval,
					     FLT_EPSILON,FLT_MIN),
			      "bad ave value for window");
		} else if (timewindow < 2 || i == 1){
		  fail_unless(FC_VALUE_EQUIV(currval,
					     origval+((j+k)*(float)(-2*glob_ftimes[i]+1))/2.0,
					     FLT_EPSILON,FLT_MIN),
			      "bad ave value for window");
		} else{ //window == 2 && i >= 2
		  fail_unless(FC_VALUE_EQUIV(currval,
					     origval+(j+k)*(float)(-6*glob_ftimes[i]+5)/3.0,
					     FLT_EPSILON,FLT_MIN),
			      "bad ave value for window");
		}
	      }
	    }
	  }
	}
	data = NULL;
	data2 = NULL;
      }
      fc_deleteSeqVariable(numStep,newSeqVar);
      free(newSeqVar);
    }

    for (window = 1; window < 7; window+=2){
      //only do window sequence for float, for others only window 3;
      if (datatype != FC_DT_FLOAT && window != 3){
	continue;
      }

      rc = fc_centeredWindowAverage(numStep,convertSeqVar,
				   window,"window",&newSeqVar);
      fail_unless(rc == FC_SUCCESS, "failed to do window average");
      fail_unless(fc_isSeqVariableValid(numStep,newSeqVar),
		  "didn't create valid seqvar");
      fc_getVariableInfo(newSeqVar[0],&newNumDataPoint,&newNumComponent,
			 &newAssoc,&newMathtype,&newDatatype);
      fail_unless(newNumDataPoint == numDataPoint,"mismatch of numDataPoint");
      fail_unless(newNumComponent == numComponent,"mismatch of numComponent");
      fail_unless(newAssoc == assoc,"mismatch of assoc");
      fail_unless(newMathtype == mathtype,"mismatch of mathtype");
      fail_unless(newDatatype == (datatype == FC_DT_FLOAT ?
				  FC_DT_FLOAT: FC_DT_DOUBLE),
		  "mismatch of datatype");
      fc_getSequenceFromSeqVariable(numStep,newSeqVar,&newSeq);
      fail_unless(FC_HANDLE_EQUIV(newSeq,lseq),
		  "mismatch of sequence handles");
      fc_getSeqVariableByName(glob_mesh,"window",&numCheckVars,&checkSteps,&checkSeqVars);
      fail_unless(numCheckVars == 1, "wrong number of matching vars");
      for (i = 0; i < numStep; i++){
	fail_unless(FC_HANDLE_EQUIV(newSeqVar[i],checkSeqVars[0][i]),
		    "mismatch of variable handles");
      }
      free(checkSeqVars[0]);
      free(checkSeqVars);
      free(checkSteps);

      //test correct vals and validity of results here
      for (i = 0; i < numStep; i++){
	rc = fc_getVariableDataPtr(convertSeqVar[i],&data);
	fail_unless(rc == FC_SUCCESS, "can't get convertSeqVar data");
	rc = fc_getVariableDataPtr(newSeqVar[i],&data2);
	fail_unless(rc == FC_SUCCESS, "can't get newSeqVar data");
	for (j = 0; j < numDataPoint; j++){
	  for (k = 0; k < numComponent; k++){
	    switch(glob_dtype[jj]){
	    case FC_DT_FLOAT:
	      {
		float currval = ((float*)data2)[numComponent*j+k];
		float origval = ((float*)data)[numComponent*j+k];
		//known vals from formula:
		if (window == 1 || i == 0 || i == (numStep -1)){
		  fail_unless(FC_VALUE_EQUIV(currval,origval,
					     FLT_EPSILON,FLT_MIN),
			      "bad ave value for window");
		} else if (window == 3  || i == 1 || i == (numStep -2)){
		  fail_unless(FC_VALUE_EQUIV(currval,
					     origval+2.0*(j+k)/3.0,
					     FLT_EPSILON,FLT_MIN),
			      "bad ave value for window");
		} else{ 
		  fail_unless(FC_VALUE_EQUIV(currval,origval+2.0*(j+k),
					     FLT_EPSILON,FLT_MIN),
			      "bad ave value for window");
		}
	      }
	      break;
	    default: 
	      {
		double currval = ((double*)data2)[numComponent*j+k];	      
		double origval;
		if (glob_dtype[jj] == FC_DT_INT){
		  origval = (double)(((int*)data)[numComponent*j+k]);
		}else {
		  origval = ((double*)data)[numComponent*j+k];
		}
		//known vals from formula. know window == 3
		if (i == 0 || i == (numStep-1)){
		  fail_unless(FC_VALUE_EQUIV(currval,origval,
					     FLT_EPSILON,FLT_MIN),
			      "bad ave value for window");
		} else{
		  fail_unless(FC_VALUE_EQUIV(currval,origval+2.0*(j+k)/3.0,
					     FLT_EPSILON,FLT_MIN),
			      "bad ave value for window");
		}
	      }
	    }
	  }
	}
	data = NULL;
	data2 = NULL;
      }
      fc_deleteSeqVariable(numStep,newSeqVar);
      free(newSeqVar);
    }

    for (window = 0; window < 6; window++){
      timewindow = (double)(window)/2.0;
      //only do window sequence for float, for others only time window = 1.5 (window  = 3)
      if (datatype != FC_DT_FLOAT && window != 3){
	continue;
      }

      rc = fc_centeredWindowAverage_Time(numStep,convertSeqVar,
				   timewindow,"window",&newSeqVar);
      fail_unless(rc == FC_SUCCESS, "failed to do window average");
      fail_unless(fc_isSeqVariableValid(numStep,newSeqVar),
		  "didn't create valid seqvar");
      fc_getVariableInfo(newSeqVar[0],&newNumDataPoint,&newNumComponent,
			 &newAssoc,&newMathtype,&newDatatype);
      fail_unless(newNumDataPoint == numDataPoint,"mismatch of numDataPoint");
      fail_unless(newNumComponent == numComponent,"mismatch of numComponent");
      fail_unless(newAssoc == assoc,"mismatch of assoc");
      fail_unless(newMathtype == mathtype,"mismatch of mathtype");
      fail_unless(newDatatype == (datatype == FC_DT_FLOAT ?
				  FC_DT_FLOAT: FC_DT_DOUBLE),
		  "mismatch of datatype");
      fc_getSequenceFromSeqVariable(numStep,newSeqVar,&newSeq);
      fail_unless(FC_HANDLE_EQUIV(newSeq,lseq),
		  "mismatch of sequence handles");
      fc_getSeqVariableByName(glob_mesh,"window",&numCheckVars,&checkSteps,&checkSeqVars);
      fail_unless(numCheckVars == 1, "wrong number of matching vars");
      for (i = 0; i < numStep; i++){
	fail_unless(FC_HANDLE_EQUIV(newSeqVar[i],checkSeqVars[0][i]),
		    "mismatch of variable handles");
      }
      free(checkSeqVars[0]);
      free(checkSeqVars);
      free(checkSteps);

      //test correct vals and validity of results here
      for (i = 0; i < numStep; i++){
	rc = fc_getVariableDataPtr(convertSeqVar[i],&data);
	fail_unless(rc == FC_SUCCESS, "can't get convertSeqVar data");
	rc = fc_getVariableDataPtr(newSeqVar[i],&data2);
	fail_unless(rc == FC_SUCCESS, "can't get newSeqVar data");
	for (j = 0; j < numDataPoint; j++){
	  for (k = 0; k < numComponent; k++){
	    switch(glob_dtype[jj]){
	    case FC_DT_FLOAT:
	      {
		float currval = ((float*)data2)[numComponent*j+k];
		float origval = ((float*)data)[numComponent*j+k];
		//known vals from formula:
		if (timewindow < 2 || i == 0 || i == (numStep-1)){
		  fail_unless(FC_VALUE_EQUIV(currval,origval,
					     FLT_EPSILON,FLT_MIN),
			      "bad ave value for window");
		} else if (timewindow < 4 || i == 1 || i ==(numStep -2)){
		  fail_unless(FC_VALUE_EQUIV(currval,
					     origval+2.0*(j+k)/3.0,
					     FLT_EPSILON,FLT_MIN),
			      "bad ave value for window");
		} else{ 
		  fail_unless(FC_VALUE_EQUIV(currval,
					     origval+2.0*(j+k),
					     FLT_EPSILON,FLT_MIN),
			      "bad ave value for window");
		}
	      }
	      break;
	    default: 
	      {
		double currval = ((double*)data2)[numComponent*j+k];	      
		double origval;
		if (glob_dtype[jj] == FC_DT_INT){
		  origval = (double)(((int*)data)[numComponent*j+k]);
		}else {
		  origval = ((double*)data)[numComponent*j+k];
		}
		//known vals from formula: window = 3
		if (timewindow < 2 || i == 0 || i == (numStep-1)){
		  fail_unless(FC_VALUE_EQUIV(currval,origval,
					     FLT_EPSILON,FLT_MIN),
			      "bad ave value for window");
		} else if (timewindow < 4 || i == 1 || i == (numStep-2)){
		  fail_unless(FC_VALUE_EQUIV(currval,
					     origval+2.0*(j+k)/3.0,
					     FLT_EPSILON,FLT_MIN),
			      "bad ave value for window");
		} else{ 
		  fail_unless(FC_VALUE_EQUIV(currval,
					     origval+2.0*(j+k),
					     FLT_EPSILON,FLT_MIN),
			      "bad ave value for window");
		}
	      }
	    }
	  }
	}
	data = NULL;
	data2 = NULL;
      }
      fc_deleteSeqVariable(numStep,newSeqVar);
      free(newSeqVar);
    }

    fc_deleteSeqVariable(numStep,convertSeqVar);
    free(convertSeqVar);
  }

  fc_deleteSeqVariable(numStep,origSeqVar);
  free(origSeqVar);



  // --- test error conditions (bad args)
  rc = fc_leadingWindowAverage(numStep,&bad_var,3,"window",&newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail if bad input seq variable");
  fail_unless(newSeqVar == NULL, "failure should return NULL seq var");
  rc = fc_leadingWindowAverage_Time(numStep,&bad_var,3,"window",&newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail if bad input seq variable");
  fail_unless(newSeqVar == NULL, "failure should return NULL seq var");
  rc = fc_centeredWindowAverage(numStep,&bad_var,3,"window",&newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail if bad input seq variable");
  fail_unless(newSeqVar == NULL, "failure should return NULL seq var");
  rc = fc_centeredWindowAverage_Time(numStep,&bad_var,3,"window",&newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail if bad input seq variable");
  fail_unless(newSeqVar == NULL, "failure should return NULL seq var");

  rc = fc_leadingWindowAverage(numStep,NULL,3,"window",&newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail if null input seq variable");
  fail_unless(newSeqVar == NULL, "failure should return NULL seq var");
  rc = fc_leadingWindowAverage_Time(numStep,NULL,3,"window",&newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail if null input seq variable");
  fail_unless(newSeqVar == NULL, "failure should return NULL seq var");
  rc = fc_centeredWindowAverage(numStep,NULL,3,"window",&newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail if null input seq variable");
  fail_unless(newSeqVar == NULL, "failure should return NULL seq var");
  rc = fc_centeredWindowAverage_Time(numStep,NULL,3,"window",&newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail if null input seq variable");
  fail_unless(newSeqVar == NULL, "failure should return NULL seq var");

  rc = fc_leadingWindowAverage(numStep,glob_sv1,0,"window",&newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail if window == 0 ");
  fail_unless(newSeqVar == NULL, "failure should return NULL seq var");
  rc = fc_leadingWindowAverage(numStep,glob_sv1,-1,"window",&newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail if window negative");
  fail_unless(newSeqVar == NULL, "failure should return NULL seq var");
  rc = fc_leadingWindowAverage_Time(numStep,glob_sv1,-1,"window",&newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail if window negative");
  fail_unless(newSeqVar == NULL, "failure should return NULL seq var");
  rc = fc_centeredWindowAverage(numStep,glob_sv1,-1,"window",&newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail if window negative");
  fail_unless(newSeqVar == NULL, "failure should return NULL seq var");
  rc = fc_centeredWindowAverage(numStep,glob_sv1,0, "window",&newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail if window = 0");
  fail_unless(newSeqVar == NULL, "failure should return NULL seq var");
  rc = fc_centeredWindowAverage(numStep,glob_sv1,4, "window",&newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail if window is even value ");
  fail_unless(newSeqVar == NULL, "failure should return NULL seq var");
  rc = fc_centeredWindowAverage_Time(numStep,glob_sv1,-1,"window",&newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail if window negative");
  fail_unless(newSeqVar == NULL, "failure should return NULL seq var");

  rc = fc_leadingWindowAverage(numStep,glob_sv1,3,NULL,&newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail if NULL name");
  fail_unless(newSeqVar == NULL, "failure should return NULL seq var");
  rc = fc_leadingWindowAverage_Time(numStep,glob_sv1,3,NULL,&newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail if NULL name");
  fail_unless(newSeqVar == NULL, "failure should return NULL seq var");
  rc = fc_centeredWindowAverage(numStep,glob_sv1,3,NULL,&newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail if NULL name");
  fail_unless(newSeqVar == NULL, "failure should return NULL seq var");
  rc = fc_centeredWindowAverage_Time(numStep,glob_sv1,3,NULL,&newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail if NULL name");
  fail_unless(newSeqVar == NULL, "failure should return NULL seq var");

  rc = fc_leadingWindowAverage(numStep,glob_sv1,3,"window",NULL);
  fail_unless(rc != FC_SUCCESS, "should fail if null return seq variable");
  rc = fc_leadingWindowAverage_Time(numStep,glob_sv1,3,"window",NULL);
  fail_unless(rc != FC_SUCCESS, "should fail if null return seq variable");
  rc = fc_centeredWindowAverage(numStep,glob_sv1,3,"window",NULL);
  fail_unless(rc != FC_SUCCESS, "should fail if null return seq variable");
  rc = fc_centeredWindowAverage_Time(numStep,glob_sv1,3,"window",NULL);
  fail_unless(rc != FC_SUCCESS, "should fail if null return seq variable");

  //check failure for invalid data types
  rc = fc_leadingWindowAverage(numStep,glob_svchar,3,"window",&newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail for datatype == FC_DT_CHAR");
  fail_unless(newSeqVar == NULL, "failure should return NULL seq var");
  rc = fc_leadingWindowAverage_Time(numStep,glob_svchar,3,"window",&newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail for datatype == FC_DT_CHAR");
  fail_unless(newSeqVar == NULL, "failure should return NULL seq var");
  rc = fc_centeredWindowAverage(numStep,glob_svchar,3,"window",&newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail for datatype == FC_DT_CHAR");
  fail_unless(newSeqVar == NULL, "failure should return NULL seq var");
  rc = fc_centeredWindowAverage_Time(numStep,glob_svchar,3,"window",&newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail for datatype == FC_DT_CHAR");
  fail_unless(newSeqVar == NULL, "failure should return NULL seq var");

  //seq types only matter for the time windows
  rc = ag_createChangedSeqDataTypeSeqVariable(numStep,glob_sv1,FC_DT_CHAR,
					      "convert",&convertSeqVar);
  fail_unless(rc == FC_SUCCESS, "couldnt convert seqvar");
  rc = fc_leadingWindowAverage_Time(numStep,convertSeqVar,3,"window",&newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail for seq datatype == FC_DT_CHAR");
  fail_unless(newSeqVar == NULL, "failure should return NULL seq var");
  rc = fc_centeredWindowAverage_Time(numStep,convertSeqVar,3,"window",&newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail for seq datatype == FC_DT_CHAR");
  fail_unless(newSeqVar == NULL, "failure should return NULL seq var");
  fc_deleteSeqVariable(numStep,convertSeqVar);
  free(convertSeqVar);

  //not checking for other mathtypes becuase all it does is copy that
  //field over when calling fc_setVariableData internally. i have
  //already checked that it does a copy correctly.

  //clean up
  newSeqVar = NULL;
  checkSeqVar = NULL;
  convertSeqVar = NULL;
}
END_TEST


START_TEST(do_comparison)
{
  //include in here tests of things on two seq var
  //right now this is euclidean dist
  //and dimensionlessAreaBetweenCurves
  FC_ReturnCode rc;

  //seq
  FC_Variable *seqVar1, *seqVar2, euc_var;
  FC_Variable *altseqVar1, *altseqVar2;
  FC_Variable *junk_sv;
  int numStep;
  int seqStep;

  //mesh
  int altnumdatapoint;

  //seqvar
  int numDataPoint;
  int numComponent, altnumcomponent;
  float *componentData1, *componentData2;
  float *comparisonData, *comparisonResult; //these are for calc the dimensionless area
  FC_AssociationType assoc = FC_AT_VERTEX;
  FC_AssociationType altassoc;
  FC_MathType mathtype, altmathtype;
  FC_DataType datatype, altdatatype;

  //bad
  FC_Variable *bad_seqVar;

  void *data;
  int i,j,k,ii;

  //setup - use global data set and mesh and seq
  //need own variables though

  //create component seq var
  mathtype = FC_MT_VECTOR;
  datatype = FC_DT_FLOAT;
  numComponent = 2;
  numDataPoint  = glob_meshDataPoint;
  numStep = glob_numStep;
  rc = fc_createSeqVariable(glob_mesh,glob_fseq,"componentseqvar1",&seqStep,
		       &seqVar1);
  fail_unless(rc == FC_SUCCESS, "failed to create 1st seqvar");
  rc = fc_createSeqVariable(glob_mesh,glob_fseq,"componentseqvar2",&seqStep,
		       &seqVar2);
  fail_unless(rc == FC_SUCCESS, "failed to create 2ndst seqvar");
  componentData1 = (float*)malloc(numDataPoint*numComponent*sizeof(float));
  componentData2 = (float*)malloc(numDataPoint*numComponent*sizeof(float));
  comparisonData = (float*)malloc(numStep*numDataPoint*numComponent*sizeof(float));
  comparisonResult = (float*)malloc(numStep*numDataPoint*numComponent*sizeof(float));
  for (i = 0; i < numDataPoint*numComponent; i++){
    comparisonResult[i] = 0.0;
  }
  for (i  = 0; i < numStep; i++){ //like each time
    for (j = 0; j < numDataPoint; j++) { //like each node
      for (k = 0; k < numComponent; k++){
	componentData1[numComponent*j+k] = (float)(j+k)+1.0;
	componentData2[numComponent*j+k] = (float)((i*i)*(j+k))+1.0;
	comparisonData[i*numDataPoint*numComponent+j*numComponent+k] = 
	  (componentData1[numComponent*j+k] - componentData2[numComponent*j+k])*
	  (componentData1[numComponent*j+k] - componentData2[numComponent*j+k]);
      }
    }

    //data for checks below
    rc = fc_setVariableData(seqVar1[i],numDataPoint,numComponent,
		       assoc,mathtype,
		       datatype,componentData1);
    fail_unless(rc == FC_SUCCESS, "failed to set data for 1st seqvar");
    rc = fc_setVariableData(seqVar2[i],numDataPoint,numComponent,
		       assoc,mathtype,
		       datatype,componentData2);
    fail_unless(rc == FC_SUCCESS, "failed to set data for 2nd seqvar");
  }

  //calcing for dimensionless area
  for (i = 0; i < numStep; i++){
    float val;
    for (j = 0; j < numDataPoint; j++) { 
      for (k = 0; k < numComponent; k++){
	 //trapezoidal rule for equally spaced abcissas.know h = 1	
	 val = comparisonData[i*numDataPoint*numComponent+j*numComponent+k];
	 comparisonResult[j*numComponent+k]+= (i == 0 || i == numStep-1) ?
	   0.5*val : val;
      }
    }
  }
  free(comparisonData);

  //will use these below, so dont free them yet
  //  free(componentData1);
  //  free(componentData2);

  //validity check
  for (ii = 0; ii < 3; ii++){
    rc = ag_createChangedDataTypeSeqVariable(numStep,seqVar1,glob_dtype[ii],
					     "convert1", &altseqVar1);
    rc = ag_createChangedDataTypeSeqVariable(numStep,seqVar2,glob_dtype[ii],
					     "convert2", &altseqVar2);

    //check euc_var
    rc = fc_euclideanDistanceBetweenCurves(numStep,altseqVar1,altseqVar2,
					"eucvar",&euc_var);
    fail_unless(rc == FC_SUCCESS, "failed to calc euclidean dist");

    fc_getVariableInfo(euc_var,&altnumdatapoint,&altnumcomponent,
		       &altassoc,&altmathtype,&altdatatype);
    //always returns doubles
    fail_unless(altdatatype == FC_DT_DOUBLE,"mismatch of datatype");

    rc = fc_getVariableDataPtr(euc_var,&data);
    fail_unless(rc == FC_SUCCESS, "failed to get euclidean data ptr");

    for (j = 0; j < altnumdatapoint; j++){
      for (k = 0; k < altnumcomponent; k++){
	double currval = ((double*)data)[altnumcomponent*j+k];
	//known value
	fail_unless(FC_VALUE_EQUIV(currval,(double)((j+k)*sqrt(299.0)),
				   DBL_EPSILON,DBL_MIN),
		    "bad value for euc dist");
      }
    }
    rc = fc_deleteVariable(euc_var);
    fail_unless(rc == FC_SUCCESS, "unable to delete euc_var");
    data = NULL;


    //check dist between same curves gives zero
    rc = fc_euclideanDistanceBetweenCurves(numStep,altseqVar1,altseqVar1,
					"eucvar",&euc_var);
    fail_unless(rc == FC_SUCCESS, "failed to calc euclidean dist");
    rc = fc_getVariableDataPtr(euc_var,&data);
    fail_unless(rc == FC_SUCCESS, "failed to get euclidean data ptr");
    for (j = 0; j < altnumdatapoint; j++){
      for (k = 0; k < altnumcomponent; k++){
	double currval = ((double*)data)[altnumcomponent*j+k];
	fail_unless(FC_DBL_EQUIV(currval,0.0),
		    "bad value for euc dist");
      }
    }
    rc = fc_deleteVariable(euc_var);
    fail_unless(rc == FC_SUCCESS, "unable to delete euc_var");
    data = NULL;


    //check dimensionlessArea
    rc = fc_dimensionlessAreaBetweenCurves(numStep,altseqVar1,altseqVar2,
					"eucvar",&euc_var);
    fail_unless(rc == FC_SUCCESS,
	"failed to calc dimensionless area between curves");

    fc_getVariableInfo(euc_var,&altnumdatapoint,&altnumcomponent,
		       &altassoc,&altmathtype,&altdatatype);
    //always returns doubles
    fail_unless(altdatatype == FC_DT_DOUBLE,"mismatch of datatype");

    rc = fc_getVariableDataPtr(euc_var,&data);
    fail_unless(rc == FC_SUCCESS, "failed to get dimensionless area data ptr");

    for (j = 0; j < altnumdatapoint; j++){
      for (k = 0; k < altnumcomponent; k++){
	double currval = ((double*)data)[altnumcomponent*j+k];
	//know max val is the last val and series size = 4
	double calcval = comparisonResult[altnumcomponent*j+k]/
	  (componentData1[altnumcomponent*j+k]*4.0);
	fail_unless(FC_DBL_EQUIV(currval,calcval),
		    "bad value for euc dist");
      }
    }

    rc = fc_deleteVariable(euc_var);
    fail_unless(rc == FC_SUCCESS, "unable to delete return var");

    data = NULL;

    //check dimensionless area between same curves gives zero
    rc = fc_dimensionlessAreaBetweenCurves(numStep,altseqVar1,altseqVar1,
					"eucvar",&euc_var);
    fail_unless(rc == FC_SUCCESS, "failed to calc euclidean dist");
    rc = fc_getVariableDataPtr(euc_var,&data);
    fail_unless(rc == FC_SUCCESS, "failed to get euclidean data ptr");
    for (j = 0; j < altnumdatapoint; j++){
      for (k = 0; k < altnumcomponent; k++){
	double currval = ((double*)data)[altnumcomponent*j+k];
	fail_unless(FC_DBL_EQUIV(currval,0.0),
		    "bad value for euc dist");
      }
    }
    rc = fc_deleteVariable(euc_var);
    fail_unless(rc == FC_SUCCESS, "unable to delete euc_var");
    data = NULL;

    //clean up for both tests
    fc_deleteSeqVariable(numStep,altseqVar1);
    free(altseqVar1);
    fc_deleteSeqVariable(numStep,altseqVar2);
    free(altseqVar2);
  }

  //clean up for dimensionless area calc
  free(comparisonResult);


  //additional special cases that work
  //see at end of test for special cases that always fail
  {
    FC_Sequence junk_seq;
    FC_Variable *junk_sv2;
    int numDataPointx = 1;
    int numComponentx = 6;
    float data1[6] = {1.0,2.0,3.0,4.0,5.0,6.0};
    float data2[6] = {2.0,3.0,4.0,5.0,6.0,7.0};
    float data3[6] = {-2.0,0.0,-4.0,-5.0,-6.0,-7.0};
    double *datax;
    int numStepx;

    //check euclidean distance works for single point sequence
    //but fails for dimensionless area.
    rc = fc_createRegularSequence(glob_dataset,1,1.2,1.3,
				  "junk_seq",&junk_seq);
    fail_unless(rc == FC_SUCCESS, "couldnt create regular sequnce");
    rc = fc_createSeqVariable(glob_mesh,junk_seq,"junk",&numStepx,
			      &junk_sv);
    fail_unless(rc == FC_SUCCESS, "failed to create variable");
    rc = fc_setVariableData(junk_sv[0],numDataPointx,numComponentx,
				FC_AT_WHOLE_MESH,FC_MT_VECTOR,
				FC_DT_FLOAT,data1);
    fail_unless(rc == FC_SUCCESS, "failed to set variable data");
    rc = fc_createSeqVariable(glob_mesh,junk_seq,"junk2",&numStepx,
			      &junk_sv2);
    fail_unless(rc ==  FC_SUCCESS,  "failed to create variable");
    rc = fc_setVariableData(junk_sv2[0],numDataPointx,numComponentx,
				FC_AT_WHOLE_MESH,FC_MT_VECTOR,
				FC_DT_FLOAT,data2);
    fail_unless(rc == FC_SUCCESS, "failed to set variable data");

    rc = fc_euclideanDistanceBetweenCurves(1,junk_sv,junk_sv2,
					"eucvar",&euc_var);
    fail_unless(rc == FC_SUCCESS, "should work for single data point sequence");
    rc = fc_getVariableDataPtr(euc_var,(void**) &datax);
    fail_unless(rc == FC_SUCCESS, "failed to get euclidean data ptr");
    for (i = 0; i < numDataPointx; i++){
      for (j = 0; j < numComponentx; j++){
	//known values - recal it returns a positive number
	fail_unless(FC_DBL_EQUIV(datax[i*numComponentx+j],1.0),
		    "bad val for euclidean distance");
      }
    }
    rc = fc_deleteVariable(euc_var);
    fail_unless(rc == FC_SUCCESS, "failed to delete eucvar");


    rc = fc_dimensionlessAreaBetweenCurves(1,junk_sv,junk_sv2,
					"eucvar",&euc_var);
    fail_unless(rc != FC_SUCCESS, "should fail for single data point sequence");
    fail_unless(FC_HANDLE_EQUIV(euc_var,FC_NULL_VARIABLE),
		" should return NULL variable when fail");

    //clean up for this test
    rc = fc_deleteSeqVariable(1,junk_sv);
    fail_unless(rc == FC_SUCCESS, "failed to delete junk_sv");
    free(junk_sv);
    rc = fc_deleteSeqVariable(1,junk_sv2);
    fail_unless(rc == FC_SUCCESS, "failed to delete junk_sv2");
    free(junk_sv2);
    rc = fc_deleteSequence(junk_seq);
    fail_unless(rc == FC_SUCCESS, "failed to delete junk_seq");



    //check failure for zero max ref in dimensionless area
    numStepx = glob_numStep;
    rc = fc_createSeqVariable(glob_mesh,glob_fseq,"junksv",&numStepx,
			      &junk_sv);
    fail_unless(rc == FC_SUCCESS, "failed to create sv");
    for (i  = 0; i < numStepx; i++){ 
      //know 3 and 2 are glob_sv1
      rc = fc_setVariableData(junk_sv[i],3,2,FC_AT_VERTEX,FC_MT_VECTOR,
			      FC_DT_FLOAT,data3);
      fail_unless(rc == FC_SUCCESS, "failed to set data for junk_sv");
    }
    //check failure 
    rc = fc_dimensionlessAreaBetweenCurves(numStepx,junk_sv,glob_sv1,
					       "eucvar",&euc_var);
    fail_unless(rc != FC_SUCCESS, "should fail for zero max in ref function");
    fail_unless(FC_HANDLE_EQUIV(euc_var,FC_NULL_VARIABLE),
		"should return NULL var when fail");
    //should work other way around
    rc = fc_dimensionlessAreaBetweenCurves(numStepx,glob_sv1,junk_sv,
					       "eucvar",&euc_var);
    fail_unless(rc == FC_SUCCESS,
		"should work for zero max only in comparison function");

    //clean up for this test
    rc = fc_deleteVariable(euc_var);
    fail_unless(rc == FC_SUCCESS, "failed to delete eucvar");
    rc = fc_deleteSeqVariable(numStepx,junk_sv);
    fail_unless(rc == FC_SUCCESS, "failed to delete junk_sv");
    free(junk_sv);
  }


  // --- test error conditions (bad args)
  // this is for both euc dist and dimensionless area
  rc = fc_euclideanDistanceBetweenCurves(numStep,&bad_var,seqVar2,
				      "eucvar",&euc_var);
  fail_unless(rc != FC_SUCCESS, "should fail if bad input seq variable 1");
  fail_unless(FC_HANDLE_EQUIV(euc_var,FC_NULL_VARIABLE),
	      " should return NULL variable when fail");
  rc = fc_dimensionlessAreaBetweenCurves(numStep,&bad_var,seqVar2,
				      "eucvar",&euc_var);
  fail_unless(rc != FC_SUCCESS, "should fail if bad input seq variable 1");
  fail_unless(FC_HANDLE_EQUIV(euc_var,FC_NULL_VARIABLE),
	      " should return NULL variable when fail");
  rc = fc_euclideanDistanceBetweenCurves(numStep,seqVar1, &bad_var,
				      "eucvar",&euc_var);
  fail_unless(rc != FC_SUCCESS, "should fail if bad input seq variable 2");
  fail_unless(FC_HANDLE_EQUIV(euc_var,FC_NULL_VARIABLE),
	      " should return NULL variable when fail");
  rc = fc_dimensionlessAreaBetweenCurves(numStep,seqVar1, &bad_var,
				      "eucvar",&euc_var);
  fail_unless(rc != FC_SUCCESS, "should fail if bad input seq variable 2");
  fail_unless(FC_HANDLE_EQUIV(euc_var,FC_NULL_VARIABLE),
	      " should return NULL variable when fail");
  rc = fc_euclideanDistanceBetweenCurves(numStep,NULL,seqVar2,
				      "eucvar",&euc_var);
  fail_unless(rc != FC_SUCCESS, "should fail if null input seq variable");
  fail_unless(FC_HANDLE_EQUIV(euc_var,FC_NULL_VARIABLE),
	      " should return NULL variable when fail");
  rc = fc_dimensionlessAreaBetweenCurves(numStep,NULL,seqVar2,
				      "eucvar",&euc_var);
  fail_unless(rc != FC_SUCCESS, "should fail if null input seq variable");
  fail_unless(FC_HANDLE_EQUIV(euc_var,FC_NULL_VARIABLE),
	      " should return NULL variable when fail");
  rc = fc_euclideanDistanceBetweenCurves(numStep,seqVar1,NULL,
				      "eucvar",&euc_var);
  fail_unless(rc != FC_SUCCESS, "should fail if null input seq variable");
  fail_unless(FC_HANDLE_EQUIV(euc_var,FC_NULL_VARIABLE),
	      " should return NULL variable when fail");
  rc = fc_dimensionlessAreaBetweenCurves(numStep,seqVar1,NULL,
				      "eucvar",&euc_var);
  fail_unless(rc != FC_SUCCESS, "should fail if null input seq variable");
  fail_unless(FC_HANDLE_EQUIV(euc_var,FC_NULL_VARIABLE),
	      " should return NULL variable when fail");
  rc = fc_euclideanDistanceBetweenCurves(numStep,seqVar1,seqVar2,NULL,
				      &euc_var);
  fail_unless(rc != FC_SUCCESS, "should fail if NULL name");
  fail_unless(FC_HANDLE_EQUIV(euc_var,FC_NULL_VARIABLE),
	      " should return NULL variable when fail");
  rc = fc_dimensionlessAreaBetweenCurves(numStep,seqVar1,seqVar2,NULL,
				      &euc_var);
  fail_unless(rc != FC_SUCCESS, "should fail if NULL name");
  fail_unless(FC_HANDLE_EQUIV(euc_var,FC_NULL_VARIABLE),
	      " should return NULL variable when fail");
  rc = fc_euclideanDistanceBetweenCurves(numStep,seqVar1,seqVar2,"euc_var",
				      NULL);
  fail_unless(rc != FC_SUCCESS, "should fail if null return seq variable");
  rc = fc_dimensionlessAreaBetweenCurves(numStep,seqVar1,seqVar2,"euc_var",
				      NULL);
  fail_unless(rc != FC_SUCCESS, "should fail if null return seq variable");


  //-- test bad cases - this is for both euc dist and dimensionless area

  //should fail if mismatching seqvars.
  //mismatching meshes
  rc = fc_createSeqVariable(bad_mesh,glob_fseq,"bad_seqVar",&seqStep,
		       &bad_seqVar);
  fail_unless(rc == FC_SUCCESS, "failed to create bad seqvar");
  for (i = 0; i < numStep; i++){
    rc = fc_setVariableData(bad_seqVar[i],numDataPoint,numComponent,
			  assoc,mathtype,
			  datatype,componentData1);
    fail_unless(rc == FC_SUCCESS, "failed to assign data to bad seqvar");
  }
  rc = fc_euclideanDistanceBetweenCurves(numStep,bad_seqVar,seqVar2,
				      "eucvar",&euc_var);
  fail_unless(rc != FC_SUCCESS, "should fail if mismatching meshes");
  fail_unless(FC_HANDLE_EQUIV(euc_var,FC_NULL_VARIABLE),
	      " should return NULL variable when fail");
  rc = fc_euclideanDistanceBetweenCurves(numStep,seqVar1,bad_seqVar,
				      "eucvar",&euc_var);
  fail_unless(rc != FC_SUCCESS, "should fail if mismatching meshes");
  fail_unless(FC_HANDLE_EQUIV(euc_var,FC_NULL_VARIABLE),
	      " should return NULL variable when fail");
  rc = fc_dimensionlessAreaBetweenCurves(numStep,bad_seqVar,seqVar2,
				      "eucvar",&euc_var);
  fail_unless(rc != FC_SUCCESS, "should fail if mismatching meshes");
  fail_unless(FC_HANDLE_EQUIV(euc_var,FC_NULL_VARIABLE),
	      " should return NULL variable when fail");
  rc = fc_dimensionlessAreaBetweenCurves(numStep,seqVar1,bad_seqVar,
				      "eucvar",&euc_var);
  fail_unless(rc != FC_SUCCESS, "should fail if mismatching meshes");
  fail_unless(FC_HANDLE_EQUIV(euc_var,FC_NULL_VARIABLE),
	      " should return NULL variable when fail");
  rc = fc_deleteSeqVariable(seqStep,bad_seqVar);
  fail_unless(rc == FC_SUCCESS, "unable to delete bad seqvar");
  free(bad_seqVar);

  //mismatching sequences
  rc = fc_createSeqVariable(glob_mesh,bad_seq,"bad_seqVar",&seqStep,
		       &bad_seqVar);
  fail_unless(rc == FC_SUCCESS, "failed to create bad seqvar");
  for (i = 0; i < glob_numStep; i++){
    rc = fc_setVariableData(bad_seqVar[i],numDataPoint,numComponent,
			    assoc,mathtype,
			    datatype,componentData2);
    fail_unless(rc == FC_SUCCESS, "failed to assign data to bad seqvar");
  }
  rc = fc_euclideanDistanceBetweenCurves(numStep,bad_seqVar,seqVar2,
				      "eucvar",&euc_var);
  fail_unless(rc != FC_SUCCESS, "should fail if mismatching sequences");
  fail_unless(FC_HANDLE_EQUIV(euc_var,FC_NULL_VARIABLE),
	      " should return NULL variable when fail");
  rc = fc_euclideanDistanceBetweenCurves(numStep,seqVar1,bad_seqVar,
				      "eucvar",&euc_var);
  fail_unless(rc != FC_SUCCESS, "should fail if mismatching sequences");
  fail_unless(FC_HANDLE_EQUIV(euc_var,FC_NULL_VARIABLE),
	      " should return NULL variable when fail");
  rc = fc_dimensionlessAreaBetweenCurves(numStep,bad_seqVar,seqVar2,
				      "eucvar",&euc_var);
  fail_unless(rc != FC_SUCCESS, "should fail if mismatching sequences");
  fail_unless(FC_HANDLE_EQUIV(euc_var,FC_NULL_VARIABLE),
	      " should return NULL variable when fail");
  rc = fc_dimensionlessAreaBetweenCurves(numStep,seqVar1,bad_seqVar,
				      "eucvar",&euc_var);
  fail_unless(rc != FC_SUCCESS, "should fail if mismatching sequences");
  fail_unless(FC_HANDLE_EQUIV(euc_var,FC_NULL_VARIABLE),
	      " should return NULL variable when fail");
  rc = fc_deleteSeqVariable(seqStep,bad_seqVar);
  fail_unless(rc == FC_SUCCESS, "unable to delete bad seqvar");
  free(bad_seqVar);


  //mismatching numDataPoint
  rc = fc_createSeqVariable(glob_mesh,glob_fseq,"bad_seqVar",&seqStep,
		       &bad_seqVar);
  fail_unless(rc == FC_SUCCESS, "failed to create to bad seqvar");
  // will check numDataPoint before the data anyway
  rc = fc_euclideanDistanceBetweenCurves(numStep,bad_seqVar,seqVar2,
				      "eucvar",&euc_var);
  fail_unless(rc != FC_SUCCESS, "should fail if mismatching numDataPoint");
  fail_unless(FC_HANDLE_EQUIV(euc_var,FC_NULL_VARIABLE),
	      " should return NULL variable when fail");
  rc = fc_euclideanDistanceBetweenCurves(numStep,seqVar1,bad_seqVar,
				      "eucvar",&euc_var);
  fail_unless(rc != FC_SUCCESS, "should fail if mismatching numDataPoint");
  fail_unless(FC_HANDLE_EQUIV(euc_var,FC_NULL_VARIABLE),
	      " should return NULL variable when fail");
  rc = fc_dimensionlessAreaBetweenCurves(numStep,bad_seqVar,seqVar2,
				      "eucvar",&euc_var);
  fail_unless(rc != FC_SUCCESS, "should fail if mismatching numDataPoint");
  fail_unless(FC_HANDLE_EQUIV(euc_var,FC_NULL_VARIABLE),
	      " should return NULL variable when fail");
  rc = fc_dimensionlessAreaBetweenCurves(numStep,seqVar1,bad_seqVar,
				      "eucvar",&euc_var);
  fail_unless(rc != FC_SUCCESS, "should fail if mismatching numDataPoint");
  fail_unless(FC_HANDLE_EQUIV(euc_var,FC_NULL_VARIABLE),
	      " should return NULL variable when fail");
  rc = fc_deleteSeqVariable(seqStep,bad_seqVar);
  fail_unless(rc == FC_SUCCESS, "unable to delete bad seqvar");
  free(bad_seqVar);


  //mismatching numComponent
  rc = fc_createSeqVariable(glob_mesh,glob_fseq,"bad_seqVar",&seqStep,
		       &bad_seqVar);
  fail_unless(rc == FC_SUCCESS, "failed to create to bad seqvar");
  for (i = 0; i < numStep; i++){
    rc = fc_setVariableData(bad_seqVar[i],numDataPoint,numComponent+1,
			    assoc,mathtype,
			    datatype,glob_ftimes);
    fail_unless(rc == FC_SUCCESS, "failed to assign data to bad seqvar");
  }
  rc = fc_euclideanDistanceBetweenCurves(numStep,bad_seqVar,seqVar1,
				      "eucvar",&euc_var);
  fail_unless(rc != FC_SUCCESS, "should fail if mismatching numComponent");
  fail_unless(FC_HANDLE_EQUIV(euc_var,FC_NULL_VARIABLE),
	      " should return NULL variable when fail");
  rc = fc_euclideanDistanceBetweenCurves(numStep,seqVar1,bad_seqVar,
				      "eucvar",&euc_var);
  fail_unless(rc != FC_SUCCESS, "should fail if mismatching numComponent");
  fail_unless(FC_HANDLE_EQUIV(euc_var,FC_NULL_VARIABLE),
	      " should return NULL variable when fail");
  rc = fc_dimensionlessAreaBetweenCurves(numStep,bad_seqVar,seqVar1,
				      "eucvar",&euc_var);
  fail_unless(rc != FC_SUCCESS, "should fail if mismatching numComponent");
  fail_unless(FC_HANDLE_EQUIV(euc_var,FC_NULL_VARIABLE),
	      " should return NULL variable when fail");
  rc = fc_dimensionlessAreaBetweenCurves(numStep,seqVar1,bad_seqVar,
				      "eucvar",&euc_var);
  fail_unless(rc != FC_SUCCESS, "should fail if mismatching numComponent");
  fail_unless(FC_HANDLE_EQUIV(euc_var,FC_NULL_VARIABLE),
	      " should return NULL variable when fail");
  rc = fc_deleteSeqVariable(seqStep,bad_seqVar);
  fail_unless(rc == FC_SUCCESS, "unable to delete bad seqvar");
  free(bad_seqVar);


  /*
  //mismatching assoc
     dont have a good way to check this because changing the
     assocation type also affects the numDataPoints, etc
  */


  //mismatching mathtype
  rc = fc_createSeqVariable(glob_mesh,glob_fseq,"bad_seqVar",&seqStep,
		       &bad_seqVar);
  fail_unless(rc == FC_SUCCESS, "failed to create to bad seqvar");
  //will check mathtype before the data anyway
  rc = fc_setVariableData(*bad_seqVar,numDataPoint,numComponent,
			  assoc,FC_MT_TENSOR,
		       datatype,componentData1);
  rc = fc_euclideanDistanceBetweenCurves(numStep,bad_seqVar,seqVar1,
				      "eucvar",&euc_var);
  fail_unless(rc != FC_SUCCESS, "should fail if mismatching mathtype");
  fail_unless(FC_HANDLE_EQUIV(euc_var,FC_NULL_VARIABLE),
	      " should return NULL variable when fail");
  rc = fc_euclideanDistanceBetweenCurves(numStep,seqVar1,bad_seqVar,
				      "eucvar",&euc_var);
  fail_unless(rc != FC_SUCCESS, "should fail if mismatching mathtype");
  fail_unless(FC_HANDLE_EQUIV(euc_var,FC_NULL_VARIABLE),
	      " should return NULL variable when fail");
  rc = fc_dimensionlessAreaBetweenCurves(numStep,bad_seqVar,seqVar1,
				      "eucvar",&euc_var);
  fail_unless(rc != FC_SUCCESS, "should fail if mismatching mathtype");
  fail_unless(FC_HANDLE_EQUIV(euc_var,FC_NULL_VARIABLE),
	      " should return NULL variable when fail");
  rc = fc_dimensionlessAreaBetweenCurves(numStep,seqVar1,bad_seqVar,
				      "eucvar",&euc_var);
  fail_unless(rc != FC_SUCCESS, "should fail if mismatching mathtype");
  fail_unless(FC_HANDLE_EQUIV(euc_var,FC_NULL_VARIABLE),
	      " should return NULL variable when fail");
  rc = fc_deleteSeqVariable(seqStep,bad_seqVar);
  fail_unless(rc == FC_SUCCESS, "unable to delete bad seqvar");
  free(bad_seqVar);


  //mismatching datatype 
  rc = fc_createSeqVariable(glob_mesh,glob_fseq,"bad_seqVar",&seqStep,
		       &bad_seqVar);
  fail_unless(rc == FC_SUCCESS, "failed to create to bad seqvar");
  //will check datatype before the data anyway
  rc = fc_setVariableData(*bad_seqVar,numDataPoint,numComponent,
			  assoc,mathtype,
			  FC_DT_INT,componentData1);
  rc = fc_euclideanDistanceBetweenCurves(numStep,bad_seqVar,seqVar1,
				      "eucvar",&euc_var);
  fail_unless(rc != FC_SUCCESS, "should fail if mismatching datatype");
  fail_unless(FC_HANDLE_EQUIV(euc_var,FC_NULL_VARIABLE),
	      " should return NULL variable when fail");
  rc = fc_euclideanDistanceBetweenCurves(numStep,seqVar1,bad_seqVar,
				      "eucvar",&euc_var);
  fail_unless(rc != FC_SUCCESS, "should fail if mismatching datatype");
  fail_unless(FC_HANDLE_EQUIV(euc_var,FC_NULL_VARIABLE),
	      " should return NULL variable when fail");
  rc = fc_dimensionlessAreaBetweenCurves(numStep,bad_seqVar,seqVar1,
				      "eucvar",&euc_var);
  fail_unless(rc != FC_SUCCESS, "should fail if mismatching datatype");
  fail_unless(FC_HANDLE_EQUIV(euc_var,FC_NULL_VARIABLE),
	      " should return NULL variable when fail");
  rc = fc_dimensionlessAreaBetweenCurves(numStep,seqVar1,bad_seqVar,
				      "eucvar",&euc_var);
  fail_unless(rc != FC_SUCCESS, "should fail if mismatching datatype");
  fail_unless(FC_HANDLE_EQUIV(euc_var,FC_NULL_VARIABLE),
	      " should return NULL variable when fail");
  rc = fc_deleteSeqVariable(seqStep,bad_seqVar);
  fail_unless(rc == FC_SUCCESS, "unable to delete bad seqvar");
  free(bad_seqVar);


  //check failure for invalid data types - char
  rc = fc_createSeqVariable(glob_mesh,glob_fseq,"bad_seqVar",&seqStep,
			    &bad_seqVar);
  fail_unless(rc == FC_SUCCESS, "unable to create bad seq var");
  rc = fc_setVariableData(*bad_seqVar,numDataPoint,numComponent,
		       assoc,mathtype, FC_DT_CHAR,componentData1);
  fail_unless(rc == FC_SUCCESS, "unable to set data for bad seq var");
  rc = fc_euclideanDistanceBetweenCurves(numStep,bad_seqVar,bad_seqVar,
				      "eucvar",&euc_var);
  fail_unless(rc != FC_SUCCESS, "should fail for datatype == FC_DT_CHAR");
  fail_unless(FC_HANDLE_EQUIV(euc_var,FC_NULL_VARIABLE),
	      " should return NULL variable when fail");
  rc = fc_dimensionlessAreaBetweenCurves(numStep,bad_seqVar,
					     bad_seqVar,
				      "eucvar",&euc_var);
  fail_unless(rc != FC_SUCCESS, "should fail for datatype == FC_DT_CHAR");
  fail_unless(FC_HANDLE_EQUIV(euc_var,FC_NULL_VARIABLE),
	      " should return NULL variable when fail");
  rc = fc_deleteSeqVariable(seqStep,bad_seqVar);
  fail_unless(rc == FC_SUCCESS, "unable to delete bad seqvar");
  free(bad_seqVar);


  //additional special cases that always fail....
  //additional failure checks for euclidean distance:
  //no failures - decided to accept :
  //     non-unique sequence (not a function)
  //     non-monotonically increasing sequence (is a function)


  //additional checks for dimensionlessArea:
  //check failure for decreasing sequence - nonmonotonically increasing
  //create seq var
  numStep = glob_numStep;
  rc = fc_createSeqVariable(glob_mesh,glob_decrseq,"junksv",&seqStep,
			    &junk_sv);
  fail_unless(rc == FC_SUCCESS, "failed to create sv");
  
  for (i  = 0; i < numStep; i++){ //like each time
    rc = fc_setVariableData(junk_sv[i],1,1,FC_AT_WHOLE_MESH,FC_MT_SCALAR,
			    FC_DT_FLOAT,&(glob_decrtimes[i]));
    fail_unless(rc == FC_SUCCESS, "failed to set data for junk_sv");
  }
  //check failure 
  rc = fc_dimensionlessAreaBetweenCurves(numStep,junk_sv,junk_sv,
				      "eucvar",&euc_var);
  fail_unless(rc != FC_SUCCESS, "should fail for decreasing sequence");
  fail_unless(FC_HANDLE_EQUIV(euc_var,FC_NULL_VARIABLE),
	      "should return NULL var when fail");
  //clean up for this test
  rc = fc_deleteSeqVariable(numStep,junk_sv);
  fail_unless(rc == FC_SUCCESS, "failed to delete junk_sv");
  free(junk_sv);


  //check failure for unchanging sequence  - not a fucntion
  //create seq var
  numStep = glob_numStep;
  rc = fc_createSeqVariable(glob_mesh,glob_nonincrseq,"junksv",&seqStep,
			    &junk_sv);
  fail_unless(rc == FC_SUCCESS, "failed to create sv");

  for (i  = 0; i < numStep; i++){ //like each time
    rc = fc_setVariableData(junk_sv[i],1,1,FC_AT_WHOLE_MESH,FC_MT_SCALAR,
			    FC_DT_FLOAT,&(glob_decrtimes[i]));
    fail_unless(rc == FC_SUCCESS, "failed to set data for junk_sv");
  }
  //check failure 
  rc = fc_dimensionlessAreaBetweenCurves(numStep,junk_sv,junk_sv,
				      "eucvar",&euc_var);
  fail_unless(rc != FC_SUCCESS, "should fail for unchanging sequence");
  fail_unless(FC_HANDLE_EQUIV(euc_var,FC_NULL_VARIABLE),
	      "should return NULL var when fail");
  //clean up for this test
  rc = fc_deleteSeqVariable(numStep,junk_sv);
  fail_unless(rc == FC_SUCCESS, "failed to delete junk_sv");
  free(junk_sv);


  //cleaning up
  free(componentData1);
  free(componentData2);

  rc = fc_deleteSeqVariable(glob_numStep,seqVar1);
  fail_unless(rc == FC_SUCCESS, "unable to delete seqvar1");
  free(seqVar1);
  rc = fc_deleteSeqVariable(glob_numStep,seqVar2);
  fail_unless(rc == FC_SUCCESS, "unable to delete seqvar2");
  free(seqVar2);

}
END_TEST


START_TEST(do_single)
{
  //too messy to put in with window
  //this tests normal form, integral
  FC_ReturnCode rc;

  //seqvar
  int numComponent;
  int numDataPoint;
  FC_AssociationType assoc;
  FC_MathType mathtype;
  FC_DataType datatype;
  int numStep;

  //compare against
  FC_Variable *newSeqVar, *checkSeqVar, *convertSeqVar, **checkSeqVars;
  FC_Variable newVar, *returnVars;
  int numReturnVars;
  FC_AssociationType newAssoc;
  FC_MathType newMathtype;
  FC_DataType newDatatype;
  FC_Sequence newSeq;
  int newNumDataPoint, newNumComponent;
  int *checkSteps, numCheckVars;
  void *data, *data2;

  
  int i,j,k,ii;

  //setup - use global dataset,mesh,seq, sv1
  for (ii= 0; ii<3; ii++){
    numStep = glob_numStep;
    rc = ag_createChangedDataTypeSeqVariable(numStep,glob_sv1,glob_dtype[ii],
					     "convert", &convertSeqVar);
    fail_unless(rc == FC_SUCCESS, "couldnt convert seqvar");

    rc = fc_getVariableInfo(convertSeqVar[0],&numDataPoint,&numComponent,
		       &assoc,&mathtype,&datatype);
    fail_unless(rc == FC_SUCCESS, "couldnt get seqvar info");


    //normal test
    rc = fc_normalForm(numStep,convertSeqVar,"junk",&newSeqVar);
    fail_unless(rc == FC_SUCCESS, "failed to do normal form");
    fail_unless(fc_isSeqVariableValid(numStep,newSeqVar),
		"didn't create valid seqvar");
    fc_getVariableInfo(newSeqVar[0],&newNumDataPoint,&newNumComponent,
		       &newAssoc,&newMathtype,&newDatatype);
    fail_unless(newNumDataPoint == numDataPoint,"mismatch of numDataPoint");
    fail_unless(newNumComponent == numComponent,"mismatch of numComponent");
    fail_unless(newAssoc == assoc,"mismatch of assoc");
    fail_unless(newMathtype == mathtype,"mismatch of mathtype");
    if (datatype == FC_DT_FLOAT){
      fail_unless(newDatatype == FC_DT_FLOAT,"mismatch of datatype");
    } else{
      fail_unless(newDatatype == FC_DT_DOUBLE,"mismatch of datatype");
    }

    fc_getSequenceFromSeqVariable(numStep,newSeqVar,&newSeq);
    fail_unless(FC_HANDLE_EQUIV(newSeq,glob_fseq),"mismatch of sequence handles");
    fc_getSeqVariableByName(glob_mesh,"junk",&numCheckVars,&checkSteps,&checkSeqVars);
    fail_unless(numCheckVars == 1, "wrong number of matching vars");
    for (i = 0; i < numStep; i++){
      fail_unless(FC_HANDLE_EQUIV(newSeqVar[i],checkSeqVars[0][i]),
		  "mismatch of variable handles");
    }
    free(checkSeqVars[0]);
    free(checkSeqVars);
    free(checkSteps);

    //test correct vals and validity of results here
    for (i = 0; i < numStep; i++){
      //      printf("\ntimestep %d\n",i);
      rc = fc_getVariableDataPtr(convertSeqVar[i],&data);
      fail_unless(rc == FC_SUCCESS, "can't get seqVar data");
      rc = fc_getVariableDataPtr(newSeqVar[i],&data2);
      fail_unless(rc == FC_SUCCESS, "can't get newSeqVar data");
      for (j = 0; j < numDataPoint; j++){
	//	printf("\tdatapoint %d ",j);
	for (k = 0; k < numComponent; k++){
	  switch(glob_dtype[ii]){
	  case FC_DT_FLOAT:
	    {
	      float currval = ((float*)data2)[numComponent*j+k];
	      float origval = ((float*)data)[numComponent*j+k];
	      double sdevx = glob_sv1_msdev[numComponent*j+k][1];
	      double meanx = glob_sv1_msdev[numComponent*j+k][0];
	      //known vals from formula:
	      if (FC_FLT_EQUIV(sdevx,0.0)){
		//		printf("(%f %f) ",origval,currval);
		fail_unless(FC_FLT_EQUIV(currval,origval),
			    "bad value for normal");
	      }else {
		float compval = (float)((origval-meanx)/sdevx);
		//		printf("(%10g %10g %10g %10g %10g) ",
		//              origval,currval,compval,
		//              FLT_EPSILON, FLT_MIN);
		//roundoff 1e-6 in comparison values
		fail_unless(FC_VALUE_EQUIV(fabs(currval-compval),
					   0.0,FLT_EPSILON,1e-6),
			    "bad value for normal");
	      }
	    }
	    break;
	  case FC_DT_INT:
	    {
	      double currval = ((double*)data2)[numComponent*j+k];
	      double origval = (double)(((int*)data)[numComponent*j+k]);
	      double sdevx = glob_sv1_msdev[numComponent*j+k][1];
	      //known val. subtract becuase of conversion to int
	      double meanx = glob_sv1_msdev[numComponent*j+k][0]-0.2;

	      //known vals from formula:
	      if (FC_DBL_EQUIV(sdevx,0.0)){
		//		printf("(%f %f) ",origval,currval);
		fail_unless(FC_DBL_EQUIV(currval,origval),
			    "bad value for normal");
	      }else {
		double compval = (double)((origval-meanx)/sdevx);
		//		printf("(%10g %10g %10g %10g %10g) ",
                //                     origval,currval,compval,
		//		       DBL_EPSILON, DBL_MIN);
		fail_unless(FC_VALUE_EQUIV(fabs(currval-compval),0.0,
					   DBL_EPSILON,1e-6),
			    "bad value for normal");
	      }
	    }
	    break;
	  default:
	    {
	      double currval = ((double*)data2)[numComponent*j+k];
	      double origval = ((double*)data)[numComponent*j+k];
	      double sdevx = glob_sv1_msdev[numComponent*j+k][1];
	      double meanx = glob_sv1_msdev[numComponent*j+k][0];
	    //known vals from formula:
	    if (FC_DBL_EQUIV(sdevx,0.0)){
	      //	      printf("(%f %f) ",origval,currval);
	      fail_unless(FC_DBL_EQUIV(currval,origval),
			  "bad value for normal");
	    }else {
	      double compval = (double)((origval-meanx)/sdevx);
	      //	      printf("(%10g %10g %10g %10g %10g) ",
              //                     origval,currval,compval,
	      //		     DBL_EPSILON, DBL_MIN);
	      fail_unless(FC_VALUE_EQUIV(fabs(currval-compval),0.0,
					 DBL_EPSILON,1e-6),
			  "bad value for normal");
	    }
	  }
	  break;
	  }
	}

      //	printf("\n");
      }
      data = NULL;
      data2 = NULL;
    }
    //clean up for this test
    rc = fc_deleteSeqVariable(numStep,newSeqVar);
    fail_unless(rc == FC_SUCCESS, "failed to delete newseqvar");
    free(newSeqVar);


    //integral test
    rc = fc_integral_TR(numStep,convertSeqVar,"junk",&newVar);
    fail_unless(rc == FC_SUCCESS, "failed to do integral");
    fail_unless(fc_isVariableValid(newVar),
		"didn't create valid newvar");
    fc_getVariableInfo(newVar,&newNumDataPoint,&newNumComponent,
		       &newAssoc,&newMathtype,&newDatatype);
    fail_unless(newNumDataPoint == numDataPoint,"mismatch of numDataPoint");
    fail_unless(newNumComponent == numComponent,"mismatch of numComponent");
    fail_unless(newAssoc == assoc,"mismatch of assoc");
    fail_unless(newMathtype == mathtype,"mismatch of mathtype");
    fail_unless(newDatatype == FC_DT_DOUBLE,"mismatch of datatype");

    rc = fc_getVariableByName(glob_mesh,"junk",&numReturnVars,&returnVars);
    fail_unless(rc == FC_SUCCESS, "failed to get integral val by name");
    fail_unless(numReturnVars == 1, "failed to get integral val by name");
    fail_unless(FC_HANDLE_EQUIV(newVar,returnVars[0]),
		"mismatch of variable handles");
    free(returnVars);


    //test correct vals and validity of results here
    rc = fc_getVariableDataPtr(newVar,&data);
    fail_unless(rc == FC_SUCCESS, "can't get newVar data");

    //known values from formula for glob_sv1
    for (j = 0; j < numDataPoint; j++){
      for (k = 0; k < numComponent; k++){
	double knownval = 8.0*(j+k)+0.8;
	double compval = ((double*)data)[numComponent*j+k];
	switch (glob_dtype[ii]){
	case FC_DT_INT:{
	  double kval = (double)((int)knownval);
	  //	  printf("compval = %10g knownval = %10g\n",compval,kval);
	  fail_unless(FC_VALUE_EQUIV(fabs(kval-compval),0.0,
				     DBL_EPSILON,1e-6),
		      "bad value for integral"); 
	  break;
	}
	case FC_DT_FLOAT:
	  //	  printf("compval = %10g knownval = %10g\n",knownval,compval);
	  fail_unless(FC_VALUE_EQUIV(fabs(knownval-compval),0.0,
				     FLT_EPSILON,1e-6),
		      "bad value for integral"); 
	  break;
	case FC_DT_DOUBLE:
	  //	  printf("compval = %10g knownval = %10g\n",knownval,compval);
	  fail_unless(FC_VALUE_EQUIV(fabs(knownval-compval),0.0,
				     DBL_EPSILON,1e-6),
		      "bad value for integral"); 
	  break;
	default:
	  fail_unless(0 == 1, "shouldnt reach this case: ERROR!");
	  break;
	}
      }
    }

    //i know i should test for if sequences are different type than just
    //float to see if it gets the val out right, but im not going to do
    //that right now. i do check for error conditions on the seq below
    // though
	  
    //clean up for this test
    rc = fc_deleteVariable(newVar);
    fail_unless(rc == FC_SUCCESS, "failed to delete newvar");


    //done testing for this iteration

    //clean up for this iteration
    fc_deleteSeqVariable(numStep,convertSeqVar);
    free(convertSeqVar);
  }
  checkSeqVar = NULL;
  convertSeqVar = NULL;

  // --- test error conditions (bad args)
  //test for normal form
  rc = fc_normalForm(numStep,&bad_var,"junk",&newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail if bad input seq variable");
  fail_unless(newSeqVar == NULL, "fail should return null seqvar");
  rc = fc_normalForm(numStep,NULL,"junk",&newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail if null input seq variable");
  fail_unless(newSeqVar == NULL, "fail should return null seqvar");
  rc = fc_normalForm(numStep,glob_sv1,NULL,&newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail if NULL name");
  fail_unless(newSeqVar == NULL, "fail should return null seqvar");
  rc = fc_normalForm(numStep,glob_sv1,"window",NULL);
  fail_unless(rc != FC_SUCCESS, "should fail if null return seq variable");

  //check failure for invalid data types
  rc = fc_normalForm(glob_numStep,glob_svchar,"window",&newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail for datatype == FC_DT_CHAR");
  fail_unless(newSeqVar == NULL, "fail should return null seqvar");

  //not checking for other mathtypes becuase all it does is copy that
  //field over when calling fc_setVariableData internally. i have
  //already checked that it does a copy correctly.

  //clean up for this test
  newSeqVar = NULL;


  //test for integral except for sequence part (bad args)
  rc = fc_integral_TR(numStep,&bad_var,"junk",&newVar);
  fail_unless(rc != FC_SUCCESS, "should fail if bad input seq variable");
  fail_unless(FC_HANDLE_EQUIV(newVar,FC_NULL_VARIABLE),
	      "should return NULL var when fail");
  rc = fc_integral_TR(numStep,NULL,"junk",&newVar);
  fail_unless(rc != FC_SUCCESS, "should fail if null input seq variable");
  fail_unless(FC_HANDLE_EQUIV(newVar,FC_NULL_VARIABLE),
	      "should return NULL var when fail");
  rc = fc_integral_TR(numStep,glob_sv1,NULL,&newVar);
  fail_unless(rc != FC_SUCCESS, "should fail if NULL name");
  fail_unless(FC_HANDLE_EQUIV(newVar,FC_NULL_VARIABLE),
	      "should return NULL var when fail");
  rc = fc_integral_TR(numStep,glob_sv1,"junk",NULL);
  fail_unless(rc != FC_SUCCESS, "should fail if null return variable");

  //check failure for invalid data types
  rc = fc_integral_TR(glob_numStep,glob_svchar,"junk",&newVar);
  fail_unless(rc != FC_SUCCESS, "should fail for datatype == FC_DT_CHAR");
  fail_unless(FC_HANDLE_EQUIV(newVar,FC_NULL_VARIABLE),
	      "should return NULL var when fail");

  //now check failures for bad sequences 
  {
    int seqStep;
    FC_Variable *junk_sv;
    FC_Sequence junk_seq;
    double data1[1] = {1.0};

    //check failure for char sequence
    //create seq var
    numStep = glob_numStep;
    rc = fc_createSeqVariable(glob_mesh,glob_cseq,"junksv",&seqStep,
			      &junk_sv);
    fail_unless(rc == FC_SUCCESS, "failed to create sv");

    for (i  = 0; i < numStep; i++){ //like each time
      rc = fc_setVariableData(junk_sv[i],1,1,FC_AT_WHOLE_MESH,FC_MT_SCALAR,
			      FC_DT_FLOAT,&(glob_ftimes[i]));
      fail_unless(rc == FC_SUCCESS, "failed to set data for junk_sv");
    }

    //check failure 
    rc = fc_integral_TR(numStep,junk_sv,"junk",&newVar);
    fail_unless(rc != FC_SUCCESS, "should fail for seqdatatype == FC_DT_CHAR");
    fail_unless(FC_HANDLE_EQUIV(newVar,FC_NULL_VARIABLE),
		"should return NULL var when fail");
    //clean up for this test
    rc = fc_deleteSeqVariable(numStep,junk_sv);
    fail_unless(rc == FC_SUCCESS, "failed to delete junk_sv");
    free(junk_sv);


    //check failure for decreasing sequence
    //create seq var
    numStep = glob_numStep;
    rc = fc_createSeqVariable(glob_mesh,glob_decrseq,"junksv",&seqStep,
			      &junk_sv);
    fail_unless(rc == FC_SUCCESS, "failed to create sv");

    for (i  = 0; i < numStep; i++){ //like each time
      rc = fc_setVariableData(junk_sv[i],1,1,FC_AT_WHOLE_MESH,FC_MT_SCALAR,
			      FC_DT_FLOAT,&(glob_ftimes[i]));
      fail_unless(rc == FC_SUCCESS, "failed to set data for junk_sv");
    }

    //check failure 
    rc = fc_integral_TR(numStep,junk_sv,"junk",&newVar);
    fail_unless(rc != FC_SUCCESS, "should fail for decreasing sequence");
    fail_unless(FC_HANDLE_EQUIV(newVar,FC_NULL_VARIABLE),
		"should return NULL var when fail");
    //clean up for this test
    rc = fc_deleteSeqVariable(numStep,junk_sv);
    fail_unless(rc == FC_SUCCESS, "failed to delete junk_sv");
    free(junk_sv);


    //check failure for unchanging sequence 
    //create seq var
    numStep = glob_numStep;
    rc = fc_createSeqVariable(glob_mesh,glob_nonincrseq,"junksv",&seqStep,
			      &junk_sv);
    fail_unless(rc == FC_SUCCESS, "failed to create sv");

    for (i  = 0; i < numStep; i++){ //like each time
      rc = fc_setVariableData(junk_sv[i],1,1,FC_AT_WHOLE_MESH,FC_MT_SCALAR,
			      FC_DT_FLOAT,glob_ftimes);
      fail_unless(rc == FC_SUCCESS, "failed to set data for junk_sv");
    }

    //check failure 
    rc = fc_integral_TR(numStep,junk_sv,"junk",&newVar);
    fail_unless(rc != FC_SUCCESS, "should fail for unchanging sequence");
    fail_unless(FC_HANDLE_EQUIV(newVar,FC_NULL_VARIABLE),
		"should return NULL var when fail");
    //clean up for this test
    rc = fc_deleteSeqVariable(numStep,junk_sv);
    fail_unless(rc == FC_SUCCESS, "failed to delete junk_sv");
    free(junk_sv);


    //check failure for single point sequence
    rc = fc_createRegularSequence(glob_dataset,1,1.2,1.3,
				  "junk_seq",&junk_seq);
    fail_unless(rc == FC_SUCCESS, "couldnt create regular sequnce");
    rc = fc_createSeqVariable(glob_mesh,junk_seq,"junk",&numStep,
			      &junk_sv);
    fail_unless(rc == FC_SUCCESS, "failed to create variable");
    rc = fc_setVariableData(junk_sv[0],1,1,
			    FC_AT_WHOLE_MESH,FC_MT_SCALAR,
			    FC_DT_DOUBLE,data1);
    fail_unless(rc == FC_SUCCESS, "failed to set variable data");
    rc = fc_integral_TR(1,junk_sv,"junk",&newVar);
    fail_unless(rc != FC_SUCCESS, "should fail for single data point sequence");
    fail_unless(FC_HANDLE_EQUIV(newVar,FC_NULL_VARIABLE),
		"should return NULL var when fail");
    //clean up for this test
    rc = fc_deleteSeqVariable(numStep,junk_sv);
    fail_unless(rc == FC_SUCCESS, "failed to delete junk_sv");
    free(junk_sv);
    rc = fc_deleteSequence(junk_seq);
    fail_unless(rc == FC_SUCCESS, "failed to delete junk_seq");
  }

  //not checking for other mathtypes becuase all it does is copy that
  //field over when calling fc_setVariableData internally. i have
  //already checked that it does a copy correctly.

  //clean up for this test - none

}
END_TEST


START_TEST(do_derivatives)
{

  // THIS IS NOT DONE YET
  //do numdp and numcomp and type checks

  //not putting in with single becuase there are restrictions
  //on number of pts int eh seq. also may put in other derivative
  //calculation and it would jst get too messy

  FC_ReturnCode rc;

  //seqvar
  int numComponent;
  int numDataPoint;
  FC_AssociationType assoc;
  FC_MathType mathtype;
  FC_DataType datatype;
  FC_Variable *seqVar;
  FC_Sequence seq;
  int numStep;
  int *itimes;
  float *ftimes, *fvals;
  double* dtimes, *dvals;

  //compare against
  FC_Variable *newSeqVar, **checkSeqVars;
  FC_AssociationType newAssoc;
  FC_MathType newMathtype;
  FC_DataType newDatatype;
  FC_Sequence newSeq;
  int newNumDataPoint, newNumComponent;
  int checkStep, *checkSteps, numCheckVars;
  //  void *data, *data2;
  void* data;

  int i,j; //,k,ii;
  
  //check legt vals
  numStep = 10;
  numDataPoint = 1;
  numComponent = 1;

  for (j = 0; j < 2; j++){
    //choosing by 0.1 since htese were the squirliest values for the reg check
    ftimes = (float*)malloc(numStep*sizeof(float));
    fvals = (float*)malloc(numDataPoint*numComponent*numStep*sizeof(float));
    for (i = 0; i < numStep; i++){
      ftimes[i] = 0.1*i;
      fvals[i] = ftimes[i]*ftimes[i];
      if (j ==1 )ftimes[i]*=-1; // j == 1 is decreasing
    }
    //create seq
    rc = fc_createSequence(glob_dataset,"origseq",&seq);
    fail_unless(rc == FC_SUCCESS, "failed to create sequence");
    rc = fc_setSequenceCoords(seq,numStep,FC_DT_FLOAT,
			      ftimes);
    fail_unless(rc == FC_SUCCESS, "failed to set sequence coords");
    
    //create seq var
    rc = fc_createSeqVariable(glob_mesh,seq,"origsv",&checkStep,
			      &seqVar);
    fail_unless(rc == FC_SUCCESS, "failed to create sv");
    
    assoc = FC_AT_WHOLE_MESH;
    datatype = FC_DT_FLOAT;
    mathtype = FC_MT_SCALAR;
    for (i  = 0; i < numStep; i++){ 
      rc = fc_setVariableData(seqVar[i],numDataPoint,numComponent,assoc,
			      mathtype,datatype,&fvals[i]);
      fail_unless(rc == FC_SUCCESS, "failed to set data for junk_sv");
    }
    
    rc = fc_firstDerivative_REA(numStep,seqVar,0,"junk",&newSeqVar);
    if (j == 0){
      fail_unless(rc == FC_SUCCESS, "failed to do first derivative REA");
      fail_unless(fc_isSeqVariableValid(numStep,newSeqVar),
		  "didn't create valid seqvar");
      fc_getVariableInfo(newSeqVar[0],&newNumDataPoint,&newNumComponent,
			 &newAssoc,&newMathtype,&newDatatype);
      fail_unless(newNumDataPoint == numDataPoint,"mismatch of numDataPoint");
      fail_unless(newNumComponent == numComponent,"mismatch of numComponent");
      fail_unless(newAssoc == assoc,"mismatch of assoc");
      fail_unless(newMathtype == mathtype,"mismatch of mathtype");
      fail_unless(newDatatype == FC_DT_FLOAT,"mismatch of datatype");
      fc_getSequenceFromSeqVariable(numStep,newSeqVar,&newSeq);
      fail_unless(FC_HANDLE_EQUIV(newSeq,seq),"mismatch of sequence handles");
      rc = fc_getSeqVariableByName(glob_mesh,"junk",&numCheckVars,&checkSteps,
				   &checkSeqVars);
      fail_unless(rc == FC_SUCCESS, "cannot get seq var by name");
      fail_unless(numCheckVars == 1, "wrong number of matching vars");
      fail_unless(checkSteps[0] == numStep, "mismatch of num steps");
      for (i = 0; i < numStep; i++){
	fail_unless(FC_HANDLE_EQUIV(newSeqVar[i],checkSeqVars[0][i]),
		    "mismatch of variable handles");
      }
      free(checkSeqVars[0]);
      free(checkSeqVars);
      free(checkSteps);
      
      //test correct vals and validity of results here
      //      printf("\n");
      for (i = 0; i < numStep; i++){
	rc = fc_getVariableDataPtr(newSeqVar[i],&data);
	fail_unless(rc == FC_SUCCESS, "can't get newVar data");
	//know only 1 dp 1 comp for now
	//	printf("%10.8g %10.8g\n",((float*)data)[0], 2.0*ftimes[i]);
	if (i < 2 || i > numStep-3){
	  fail_unless(FC_FLT_EQUIV(((float*)data)[0], 0.0),
		      "mismatch of return values");
	}else {
	  fail_unless(FC_VALUE_EQUIV(((float*)data)[0], 2.0*ftimes[i],10*FLT_EPSILON,FLT_MIN),
		      "mismatch of return values");
	}
      }
      rc = fc_deleteSeqVariable(numStep,newSeqVar);
      fail_unless(rc == FC_SUCCESS, "failed to delete newSeqVar");
    } else {
      fail_unless(rc != FC_SUCCESS, "should fail for decreasing seq");
    }
  
    //clean up for this test
    rc = fc_deleteSeqVariable(numStep,seqVar);
    fail_unless(rc == FC_SUCCESS, "failed to delete orig seq var");
    rc = fc_deleteSequence(seq);
    fail_unless(rc == FC_SUCCESS, "failed to delete seq");
    free(ftimes);
    free(fvals);
    free(seqVar);
    free(newSeqVar);
  }

  //***TO DO***  check diff num dp num comp

  //note: check for type double and int is below

  //---------special cases 
  //infinte is same as unchanging sequence (see below)
  //0 derivative
  //  printf("zero deriviative\n");
  numStep = 5;
  ftimes = (float*)malloc(numStep*sizeof(float));
  fvals = (float*)malloc(numStep*sizeof(float));
  for (i = 0; i < numStep; i++){
    ftimes[i] = 0.1*i;
    fvals[i] = 0.1;
  }

  //create seq
  rc = fc_createSequence(glob_dataset,"origseq",&seq);
  fail_unless(rc == FC_SUCCESS, "failed to create sequence");
  rc = fc_setSequenceCoords(seq,numStep,FC_DT_FLOAT,
			      ftimes);
  fail_unless(rc == FC_SUCCESS, "failed to set sequence coords");

  //create seq var
  rc = fc_createSeqVariable(glob_mesh,seq,"origsv",&checkStep,
			      &seqVar);
  fail_unless(rc == FC_SUCCESS, "failed to create sv");

  for (i  = 0; i < numStep; i++){ 
    rc = fc_setVariableData(seqVar[i],1,1,FC_AT_WHOLE_MESH,FC_MT_SCALAR,
			    FC_DT_FLOAT,&fvals[i]);
    fail_unless(rc == FC_SUCCESS, "failed to set data for seqvar");
  }

  rc = fc_firstDerivative_REA(numStep,seqVar,0,"junk",&newSeqVar);
  fail_unless(rc == FC_SUCCESS, "should work for zero derivative");
  for (i = 0; i < numStep; i++){
    rc = fc_getVariableDataPtr(newSeqVar[i],&data);
    fail_unless(rc == FC_SUCCESS, "can't get newVar data");
    //know only 1 dp 1 comp for now
    //    printf("%10.8g %10.8g\n",((float*)data)[0], 0.0); 
    fail_unless(FC_FLT_EQUIV(1.0+((float*)data)[0], 1.0),
		"mismatch of return values");
  }
  
  //clean up for this test
  rc = fc_deleteSeqVariable(numStep,seqVar);
  fail_unless(rc == FC_SUCCESS, "failed to delete orig seq var");
  rc = fc_deleteSeqVariable(numStep,newSeqVar);
  fail_unless(rc == FC_SUCCESS, "failed to delete newSeqVar");
  rc = fc_deleteSequence(seq);
  fail_unless(rc == FC_SUCCESS, "failed to delete seq");
  free(ftimes);
  free(fvals);
  free(seqVar);
  free(newSeqVar);


  // check int with datatype double vs float with datatype double
  //could check for more combos....
  for (j = 0; j < 2; j++){
    //    printf("checking mixed\n");
    numStep = 10;
    itimes = (int*)malloc(numStep*sizeof(int));
    ftimes = (float*)malloc(numStep*sizeof(float));
    dvals = (double*)malloc(numStep*sizeof(double));
    if (j == 0){
      for (i = 0; i < numStep; i++){
	itimes[i] = i;
	dvals[i] = (double)(itimes[i]*itimes[i]);
      }
    } else{
      for (i = 0; i < numStep; i++){
	ftimes[i] = i;
	dvals[i] = (double)(ftimes[i]*ftimes[i]);
      }
    }


    //create seq
    rc = fc_createSequence(glob_dataset,"origseq",&seq);
    fail_unless(rc == FC_SUCCESS, "failed to create sequence");
    if (j ==0){
      rc = fc_setSequenceCoords(seq,numStep,FC_DT_INT,
				itimes);
    } else {
      rc = fc_setSequenceCoords(seq,numStep,FC_DT_FLOAT,
				ftimes);
    }

    fail_unless(rc == FC_SUCCESS, "failed to set sequence coords");
    
    //create seq var
    rc = fc_createSeqVariable(glob_mesh,seq,"origsv",&checkStep,
			      &seqVar);
    fail_unless(rc == FC_SUCCESS, "failed to create sv");
    
    for (i  = 0; i < numStep; i++){ 
      rc = fc_setVariableData(seqVar[i],1,1,FC_AT_WHOLE_MESH,FC_MT_SCALAR,
			      FC_DT_DOUBLE,&dvals[i]);
      fail_unless(rc == FC_SUCCESS, "failed to set data for seqvar");
    }
    
    rc = fc_firstDerivative_REA(numStep,seqVar,0,"junk",&newSeqVar);
    //should work
    fail_unless(rc == FC_SUCCESS, "failed to do first derivative REA");
    fail_unless(fc_isSeqVariableValid(numStep,newSeqVar),
		"didn't create valid seqvar");
    fc_getVariableInfo(newSeqVar[0],&newNumDataPoint,&newNumComponent,
		       &newAssoc,&newMathtype,&newDatatype);
    fail_unless(newNumDataPoint == numDataPoint,"mismatch of numDataPoint");
    fail_unless(newNumComponent == numComponent,"mismatch of numComponent");
    fail_unless(newAssoc == assoc,"mismatch of assoc");
    fail_unless(newMathtype == mathtype,"mismatch of mathtype");
    if (j == 0){
      fail_unless(newDatatype == FC_DT_DOUBLE,"mismatch of datatype");
    } else {
      fail_unless(newDatatype == FC_DT_FLOAT,"mismatch of datatype");
    }
    fc_getSequenceFromSeqVariable(numStep,newSeqVar,&newSeq);
    fail_unless(FC_HANDLE_EQUIV(newSeq,seq),"mismatch of sequence handles");
    fc_getSeqVariableByName(glob_mesh,"junk",&numCheckVars,&checkSteps,&checkSeqVars);
    fail_unless(numCheckVars == 1, "wrong number of matching vars");
    fail_unless(checkSteps[0] == numStep, "mismatch of num steps");
    for (i = 0; i < numStep; i++){
      fail_unless(FC_HANDLE_EQUIV(newSeqVar[i],checkSeqVars[0][i]),
		  "mismatch of variable handles");
    }
    free(checkSeqVars[0]);
    free(checkSeqVars);
    free(checkSteps);

    //test correct vals and validity of results here
    //      printf("\n");
    for (i = 0; i < numStep; i++){
      rc = fc_getVariableDataPtr(newSeqVar[i],&data);
      fail_unless(rc == FC_SUCCESS, "can't get newVar data");
      //know only 1 dp 1 comp for now
      if (j == 0){
	//	printf("%10.8g %10d\n",((double*)data)[0], 2*itimes[i]);
	if (i < 2 || i > numStep-3){
	  fail_unless(FC_DBL_EQUIV(((double*)data)[0], 0.0),
		      "mismatch of return values");
	}else {
	  fail_unless(FC_DBL_EQUIV(((double*)data)[0],2*itimes[i]),
		      "mismatch of return values");
	}
      }else {
	//	printf("%10.8g %10.8g\n",((float*)data)[0], 2*ftimes[i]);
	if (i < 2 || i > numStep-3){
	  fail_unless(FC_DBL_EQUIV(((float*)data)[0], 0.0),
		      "mismatch of return values");
	}else {
	  fail_unless(FC_DBL_EQUIV(((float*)data)[0],2*ftimes[i]),
		      "mismatch of return values");
	}
      }
    }

    //clean up for this test
    rc = fc_deleteSeqVariable(numStep,newSeqVar);
    fail_unless(rc == FC_SUCCESS, "failed to delete new seq var");
    rc = fc_deleteSeqVariable(numStep,seqVar);
    fail_unless(rc == FC_SUCCESS, "failed to delete orig seq var");
    rc = fc_deleteSequence(seq);
    fail_unless(rc == FC_SUCCESS, "failed to delete seq");
    free(itimes);
    free(ftimes);
    free(dvals);
    free(seqVar);
    free(newSeqVar);
  }

  // --- test error conditions (bad args)
  //test for derivatives except for sequence part (bad args)
  numStep = 5;
  rc = fc_firstDerivative_REA(numStep,&bad_var,0,"junk",&newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail if bad input seq variable");
  fail_unless(newSeqVar == NULL,
	      "should return NULL var when fail");
  rc = fc_firstDerivative_REA(numStep,NULL,0,"junk",&newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail if null input seq variable");
  fail_unless(newSeqVar == NULL,
	      "should return NULL var when fail");
  rc = fc_firstDerivative_REA(numStep,glob_sv1,0,NULL,&newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail if NULL name");
  fail_unless(newSeqVar == NULL,
	      "should return NULL var when fail");
  rc = fc_firstDerivative_REA(numStep,glob_sv1,0,"junk",NULL);
  fail_unless(rc != FC_SUCCESS, "should fail if null return variable");

  //check failure for invalid data types
  rc = fc_firstDerivative_REA(glob_numStep,glob_svchar,0,"junk",&newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail for datatype == FC_DT_CHAR");
  fail_unless(newSeqVar == NULL,
	      "should return NULL var when fail");
  
  //------now check failures for bad sequences 
  {
    int seqStep;
    FC_Variable *junk_sv;

    //check failure for char sequence
    //create seq var
    numStep = glob_numStep;
    rc = fc_createSeqVariable(glob_mesh,glob_cseq,"junksv",&seqStep,
			      &junk_sv);
    fail_unless(rc == FC_SUCCESS, "failed to create sv");

    for (i  = 0; i < numStep; i++){ //like each time - data doesnt matter
      rc = fc_setVariableData(junk_sv[i],1,1,FC_AT_WHOLE_MESH,FC_MT_SCALAR,
			      FC_DT_FLOAT,&(glob_ftimes[i]));
      fail_unless(rc == FC_SUCCESS, "failed to set data for junk_sv");
    }

    //check failure 
    rc = fc_firstDerivative_REA(numStep,junk_sv,0,"junk",&newSeqVar);
    fail_unless(rc != FC_SUCCESS, "should fail for seqdatatype == FC_DT_CHAR");
    fail_unless(newSeqVar == NULL,
		"should return NULL var when fail");
    //clean up for this test
    rc = fc_deleteSeqVariable(numStep,junk_sv);
    fail_unless(rc == FC_SUCCESS, "failed to delete junk_sv");
    free(junk_sv);
  }


  // check 1) double = ok  and 2) failure for irregular sequence
  for (j = 0; j < 2; j++){
    //    printf ("\n IRREG SEQ\n");
    numStep = 10;
    dtimes = (double*)malloc(numStep*sizeof(double));
    dvals = (double*)malloc(numStep*sizeof(double));
    for (i = 0; i < numStep; i++){
      dtimes[i] = 0.1*i;
      dvals[i] = dtimes[i]*dtimes[i];
    }
    if (j == 1){ //irregular
      dtimes[3]+= FLT_EPSILON;
    }
    //    printf("times: ");
    //    for (i  = 0; i < numStep; i++){
    //      printf ("%10.8g ",dtimes[i]);
    //    }
    //    printf("\n");
    
    //create seq
    rc = fc_createSequence(glob_dataset,"origseq",&seq);
    fail_unless(rc == FC_SUCCESS, "failed to create sequence");
    rc = fc_setSequenceCoords(seq,numStep,FC_DT_DOUBLE,
			      dtimes);
    fail_unless(rc == FC_SUCCESS, "failed to set sequence coords");
    
    //create seq var
    rc = fc_createSeqVariable(glob_mesh,seq,"origsv",&checkStep,
			      &seqVar);
    fail_unless(rc == FC_SUCCESS, "failed to create sv");
    
    for (i  = 0; i < numStep; i++){ 
      rc = fc_setVariableData(seqVar[i],1,1,FC_AT_WHOLE_MESH,FC_MT_SCALAR,
			      FC_DT_DOUBLE,&dvals[i]);
      fail_unless(rc == FC_SUCCESS, "failed to set data for seqvar");
    }

    rc = fc_firstDerivative_REA(numStep,seqVar,0,"junk",&newSeqVar);
    if (j == 0){
      //should work
      fail_unless(rc == FC_SUCCESS, "should work for double series");
      fail_unless(rc == FC_SUCCESS, "failed to do first derivative REA");
      fail_unless(fc_isSeqVariableValid(numStep,newSeqVar),
		  "didn't create valid seqvar");
      fc_getVariableInfo(newSeqVar[0],&newNumDataPoint,&newNumComponent,
			 &newAssoc,&newMathtype,&newDatatype);
      fail_unless(newNumDataPoint == numDataPoint,"mismatch of numDataPoint");
      fail_unless(newNumComponent == numComponent,"mismatch of numComponent");
      fail_unless(newAssoc == assoc,"mismatch of assoc");
      fail_unless(newMathtype == mathtype,"mismatch of mathtype");
      fail_unless(newDatatype == FC_DT_DOUBLE,"mismatch of datatype");
      fc_getSequenceFromSeqVariable(numStep,newSeqVar,&newSeq);
      fail_unless(FC_HANDLE_EQUIV(newSeq,seq),"mismatch of sequence handles");
      fc_getSeqVariableByName(glob_mesh,"junk",&numCheckVars,&checkSteps,&checkSeqVars);
      fail_unless(numCheckVars == 1, "wrong number of matching vars");
      fail_unless(checkSteps[0] == numStep, "mismatch of num steps");
      for (i = 0; i < numStep; i++){
	fail_unless(FC_HANDLE_EQUIV(newSeqVar[i],checkSeqVars[0][i]),
		    "mismatch of variable handles");
      }
      free(checkSeqVars[0]);
      free(checkSeqVars);
      free(checkSteps);

      //test correct vals and validity of results here
      //      printf("\n");
      for (i = 0; i < numStep; i++){
	rc = fc_getVariableDataPtr(newSeqVar[i],&data);
	fail_unless(rc == FC_SUCCESS, "can't get newVar data");
	//know only 1 dp 1 comp for now
	//	printf("%10.8g %10.8g\n",((double*)data)[0], 2.0*dtimes[i]);
	if (i < 2 || i > numStep-3){
	  fail_unless(FC_DBL_EQUIV(((double*)data)[0], 0.0),
		      "mismatch of return values");
	}else {
	  fail_unless(FC_VALUE_EQUIV(((double*)data)[0],
				     2.0*dtimes[i],10*DBL_EPSILON,DBL_MIN),
		      "mismatch of return values");
	}
      }
      rc = fc_deleteSeqVariable(numStep,newSeqVar);
      fail_unless(rc == FC_SUCCESS, "failed to delete new seq var");
      free(newSeqVar);
    } else {
      fail_unless(rc != FC_SUCCESS, "should fail for irreg series");
      fail_unless(newSeqVar == NULL,
		  "should return NULL var when fail");
    }
  
    //clean up for this test
    rc = fc_deleteSeqVariable(numStep,seqVar);
    fail_unless(rc == FC_SUCCESS, "failed to delete orig seq var");
    rc = fc_deleteSequence(seq);
    fail_unless(rc == FC_SUCCESS, "failed to delete seq");
    free(dtimes);
    free(dvals);
    free(seqVar);
  }


  //check failure for unchanging sequence 
  numStep = 5;
  ftimes = (float*)malloc(numStep*sizeof(float));
  fvals = (float*)malloc(numStep*sizeof(float));
  for (i = 0; i < numStep; i++){
    ftimes[i] = 0.1;
    fvals[i] = 0.1*i;
  }

  //create seq
  rc = fc_createSequence(glob_dataset,"origseq",&seq);
  fail_unless(rc == FC_SUCCESS, "failed to create sequence");
  rc = fc_setSequenceCoords(seq,numStep,FC_DT_FLOAT,
			      ftimes);
  fail_unless(rc == FC_SUCCESS, "failed to set sequence coords");

  //create seq var
  rc = fc_createSeqVariable(glob_mesh,seq,"origsv",&checkStep,
			      &seqVar);
  fail_unless(rc == FC_SUCCESS, "failed to create sv");

  for (i  = 0; i < numStep; i++){ 
    rc = fc_setVariableData(seqVar[i],1,1,FC_AT_WHOLE_MESH,FC_MT_SCALAR,
			    FC_DT_FLOAT,&fvals[i]);
    fail_unless(rc == FC_SUCCESS, "failed to set data for seqvar");
  }

  rc = fc_firstDerivative_REA(numStep,seqVar,0,"junk",&newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail for unchaning seq");
  fail_unless(newSeqVar == NULL,
	      "should return NULL var when fail");
  
  //clean up for this test
  rc = fc_deleteSeqVariable(numStep,seqVar);
  fail_unless(rc == FC_SUCCESS, "failed to delete orig seq var");
  rc = fc_deleteSequence(seq);
  fail_unless(rc == FC_SUCCESS, "failed to delete seq");
  free(ftimes);
  free(fvals);
  free(seqVar);


  //check failure for seq had < 5 pts
  numStep = 4;
  ftimes = (float*)malloc(numStep*sizeof(float));
  fvals = (float*)malloc(numStep*sizeof(float));
  for (i = 0; i < numStep; i++){
    ftimes[i] = 0.1*i;
    fvals[i] = 0.1*i;
  }
  
  //create seq
  rc = fc_createSequence(glob_dataset,"origseq",&seq);
  fail_unless(rc == FC_SUCCESS, "failed to create sequence");
  rc = fc_setSequenceCoords(seq,numStep,FC_DT_FLOAT,
			      ftimes);
  fail_unless(rc == FC_SUCCESS, "failed to set sequence coords");

  //create seq var
  rc = fc_createSeqVariable(glob_mesh,seq,"origsv",&checkStep,
			      &seqVar);
  fail_unless(rc == FC_SUCCESS, "failed to create sv");

  for (i  = 0; i < numStep; i++){ 
    rc = fc_setVariableData(seqVar[i],1,1,FC_AT_WHOLE_MESH,FC_MT_SCALAR,
			    FC_DT_FLOAT,&fvals[i]);
    fail_unless(rc == FC_SUCCESS, "failed to set data for seqvar");
  }

  rc = fc_firstDerivative_REA(numStep,seqVar,0,"junk",&newSeqVar);
  fail_unless(rc != FC_SUCCESS, "should fail for seq too short");
  fail_unless(newSeqVar == NULL,
	      "should return NULL var when fail");
  
  //clean up for this test
  rc = fc_deleteSeqVariable(numStep,seqVar);
  fail_unless(rc == FC_SUCCESS, "failed to delete orig seq var");
  rc = fc_deleteSequence(seq);
  fail_unless(rc == FC_SUCCESS, "failed to delete seq");
  free(ftimes);
  free(fvals);
  free(seqVar);

  //not checking for other mathtypes becuase all it does is copy that
  //field over when calling fc_setVariableData internally. 
  //clean up for this test - none
}
END_TEST


START_TEST(do_squares)
{
  // checks linearleastsquares
  // too messy to put in with other things.
  // may split out more of the other stuff
  // as this goes on
  FC_ReturnCode rc;

  //seqvar
  FC_Variable *baseVar, *baseCharVar;
  int numComponent;
  int numDataPoint;
  FC_AssociationType assoc;
  FC_MathType mathtype;
  FC_DataType datatype;
  int numStep;
  float *componentData;

  //compare against
  FC_Variable *checkSeqVar, *convertSeqVar;
  FC_Variable newVar, *returnVars;
  int numReturnVars;
  FC_AssociationType newAssoc;
  FC_MathType newMathtype;
  FC_DataType newDatatype;
  FC_Sequence newSeq;
  int* newISeqCoords; 
  double* newDSeqCoords;
  int newNumDataPoint, newNumComponent;
  void *data;

  int i,j,ii;

  //make a single component seq var out of the 1st component of glob_sv1
  //with a little alteration....
  //so copied that code over
  rc = fc_createSeqVariable(glob_mesh,glob_fseq,glob_sv1_name,&numStep,
		       &baseVar);
  fail_unless(rc == FC_SUCCESS, "failed to create baseVar");
  rc = fc_createSeqVariable(glob_mesh,glob_fseq,glob_svchar_name,&numStep,
		       &baseCharVar);
  fail_unless(rc == FC_SUCCESS, "failed to create baseCharVar");

  numDataPoint = glob_meshDataPoint;

  componentData = (float*)malloc(numDataPoint*sizeof(float));
  for (i  = 0; i < glob_numStep; i++){ //like each time
    //first fill in same data as glob_sv1
    for (j = 0; j < glob_meshDataPoint; j++) { //like each node
      componentData[j] = glob_ftimes[i] * (float)j+0.2;
      //if its the last data point change the sign
      if (j == glob_meshDataPoint-1){
	componentData[j]*=-1;
      }
    }

    rc = fc_setVariableData(baseVar[i],numDataPoint,1,
		       FC_AT_VERTEX,FC_MT_SCALAR,
		       FC_DT_FLOAT,componentData);
    fail_unless(rc == FC_SUCCESS, "failed to set data for baseVar");
    rc = fc_setVariableData(baseCharVar[i],numDataPoint,1,
			    FC_AT_VERTEX,FC_MT_SCALAR,
			    FC_DT_CHAR,glob_ctimes);
    fail_unless(rc == FC_SUCCESS, "failed to set data for baseCharVar");
  }
  free(componentData);


  //now do same as in others, only from local vars rather than glob vars
  for (ii= 0; ii<3; ii++){
    numStep = glob_numStep;
    rc = ag_createChangedDataTypeSeqVariable(numStep,baseVar,glob_dtype[ii],
					     "convert", &convertSeqVar);
    fail_unless(rc == FC_SUCCESS, "couldnt convert seqvar");

    rc = fc_getVariableInfo(convertSeqVar[0],&numDataPoint,&numComponent,
		       &assoc,&mathtype,&datatype);
    fail_unless(rc == FC_SUCCESS, "couldnt get seqvar info");

    //linearLeastSquares test
    rc = fc_linearLeastSquares(numStep,convertSeqVar,&newVar);
    fail_unless(rc == FC_SUCCESS, "failed to do linearLeastSquares");
    fail_unless(fc_isVariableValid(newVar),
		"didn't create valid newvar");
    fc_getVariableInfo(newVar,&newNumDataPoint,&newNumComponent,
		       &newAssoc,&newMathtype,&newDatatype);
    fail_unless(newNumDataPoint == numDataPoint,"mismatch of numDataPoint");
    fail_unless(newNumComponent == 5,"mismatch of numComponent");
    fail_unless(newAssoc == assoc,"mismatch of assoc");
    fail_unless(newMathtype == FC_MT_VECTOR,"mismatch of mathtype");
    fail_unless(newDatatype == FC_DT_DOUBLE,"mismatch of datatype");

    rc = fc_getVariableByName(glob_mesh,"convert_LSQ",&numReturnVars,
			      &returnVars);
    fail_unless(rc == FC_SUCCESS, "failed to get integral val by name");
    fail_unless(numReturnVars == 1, "failed to get integral val by name");
    fail_unless(FC_HANDLE_EQUIV(newVar,returnVars[0]),
		"mismatch of variable handles");
    free(returnVars);

    //test correct vals and validity of results here
    rc = fc_getVariableDataPtr(newVar,&data);
    fail_unless(rc == FC_SUCCESS, "can't get newVar data");
    for (j = 0; j < numDataPoint; j++){
      double a,b,r2,sea,seb, intercept;
      a = ((double*)data)[5*j];
      b = ((double*)data)[5*j+1];
      r2 = ((double*)data)[5*j+2];
      sea = ((double*)data)[5*j+3];
      seb = ((double*)data)[5*j+4];

      //      printf("a %10g b %10g r2 %10g sea %10g seb %10g\n",
      //	     a,b,r2,sea,seb);
      //known values
      //dont know why should have to use FLT and not DBL here
      if (j !=2){
	fail_unless(FC_FLT_EQUIV(b,j), "bad slope");
	intercept =  (glob_dtype[ii] == FC_DT_INT) ? 0.0 : 0.2;
      }else {
	fail_unless(FC_FLT_EQUIV(b,-1*j), "bad slope");
	intercept =  (glob_dtype[ii] == FC_DT_INT) ? 0.0 : -0.2;
      }
      fail_unless(FC_VALUE_EQUIV(a,intercept,10*FLT_EPSILON,FLT_MIN),
		  "bad intercept");
      switch(j){
      case 0:
	//know ssyy == 0 , ssxy == 0 
	fail_unless(FC_FLT_EQUIV(r2,-1.0), "bad correlation coefficient");
	fail_unless(FC_FLT_EQUIV(sea,0.0), "bad std error a");
	fail_unless(FC_FLT_EQUIV(seb,0.0), "bad std error b");
	break;
      case 1:
	//know ssxx == ssxy == ssyy
	fail_unless(FC_FLT_EQUIV(r2,1.0), "bad correlation coefficient");
	fail_unless(FC_FLT_EQUIV(sea,0.0), "bad std error a");
	fail_unless(FC_FLT_EQUIV(seb,0.0), "bad std error b");
	break;
      case 2:
	//know ssyy = ssxy2/ssxx 
	fail_unless(FC_FLT_EQUIV(r2,1.0), "bad correlation coefficient");
	fail_unless(FC_FLT_EQUIV(sea,0.0), "bad std error a");
	fail_unless(FC_FLT_EQUIV(seb,0.0), "bad std error b");
	break;
      default:
	fail_unless(0 == 1,"Developer error - there is no other test....");
	break;
      }
    }

    //clean up for this test
    data = NULL;
    rc = fc_deleteVariable(newVar);
    fail_unless(rc == FC_SUCCESS, "failed to delete newvar");

    //done testing for this iteration

    //clean up for this iteration
    fc_deleteSeqVariable(numStep,convertSeqVar);
    free(convertSeqVar);
  }

  checkSeqVar = NULL;
  convertSeqVar = NULL;


  //now lets try something thats not a pure line
  fc_deleteSeqVariable(numStep,baseVar);
  free(baseVar);
  rc = fc_createSeqVariable(glob_mesh,glob_fseq,glob_sv1_name,&numStep,
		       &baseVar);
  fail_unless(rc == FC_SUCCESS, "failed to create baseVar");
  numDataPoint = 1;
  componentData = (float*)malloc(numDataPoint*sizeof(float));
  for (i  = 0; i < glob_numStep; i++){ //like each time
    componentData[0] = glob_ftimes[i];
    //if its the first or last time step, change it
    if (i == 0){
      componentData[0]--;
    }
    if (i == glob_numStep-1){
      componentData[0]++;
    }
    rc = fc_setVariableData(baseVar[i],numDataPoint,1,
		       FC_AT_WHOLE_MESH,FC_MT_SCALAR,
		       FC_DT_FLOAT,componentData);
    fail_unless(rc == FC_SUCCESS, "failed to set data for baseVar");
  }
  free(componentData);
  rc = fc_linearLeastSquares(numStep,baseVar,&newVar);
  fail_unless(rc == FC_SUCCESS, "failed to do linearLeastSquares");
  fail_unless(fc_isVariableValid(newVar),
	      "didn't create valid newvar");
  rc = fc_getVariableDataPtr(newVar,&data);
  fail_unless(rc == FC_SUCCESS, "can't get newVar data");
  {
    double a,b,r2,sea,seb, s;
    a = ((double*)data)[0];
    b = ((double*)data)[1];
    r2 = ((double*)data)[2];
    sea = ((double*)data)[3];
    seb = ((double*)data)[4];
    
    //    printf("a %10g b %10g r2 %10g sea %10g seb %10g\n",
    //	   a,b,r2,sea,seb);
    //known values ssxx = 10 ssxy = 14 ssyy = 20
    fail_unless(FC_FLT_EQUIV(b,1.4), "bad slope");
    fail_unless(FC_VALUE_EQUIV(a,-0.8,FLT_EPSILON,FLT_MIN),
		"bad intercept");
    fail_unless(FC_FLT_EQUIV(r2,0.98), "bad correlation coefficient");
    s = sqrt(0.133333333333);
    fail_unless(FC_FLT_EQUIV(sea,s*sqrt(0.6)), "bad std error a");
    fail_unless(FC_FLT_EQUIV(seb,s/sqrt(10.0)), "bad std error b");
  }
    
  //clean up for this test
  data = NULL;
  rc = fc_deleteVariable(newVar);
  fail_unless(rc == FC_SUCCESS, "failed to delete newvar");

  //combining checking with double and int sequence below with
  //other error checking

  // --- test error conditions (bad args)
  //test for sum of squares except for sequence part
  rc = fc_linearLeastSquares(numStep,glob_sv1,&newVar);
  fail_unless(rc != FC_SUCCESS, "should fail if multiple component data");
  fail_unless(FC_HANDLE_EQUIV(newVar,FC_NULL_VARIABLE),
	      "should return NULL var when fail");
  rc = fc_linearLeastSquares(numStep,&bad_var,&newVar);
  fail_unless(rc != FC_SUCCESS, "should fail if bad input seq variable");
  fail_unless(FC_HANDLE_EQUIV(newVar,FC_NULL_VARIABLE),
	      "should return NULL var when fail");
  rc = fc_linearLeastSquares(numStep,NULL,&newVar);
  fail_unless(rc != FC_SUCCESS, "should fail if null input seq variable");
  fail_unless(FC_HANDLE_EQUIV(newVar,FC_NULL_VARIABLE),
	      "should return NULL var when fail");
  rc = fc_linearLeastSquares(numStep,baseVar,NULL);
  fail_unless(rc != FC_SUCCESS, "should fail if null return variable");

  //check failure for invalid data types
  rc = fc_linearLeastSquares(glob_numStep,baseCharVar,&newVar);
  fail_unless(rc != FC_SUCCESS, "should fail for datatype == FC_DT_CHAR");
  fail_unless(FC_HANDLE_EQUIV(newVar,FC_NULL_VARIABLE),
	      "should return NULL var when fail");
  //clean up for this test 
  fc_deleteSeqVariable(numStep,baseVar);
  free(baseVar);
  fc_deleteSeqVariable(numStep,baseCharVar);
  free(baseCharVar);


  //check failure for char sequence
  numStep = glob_numStep;

  rc = fc_createSeqVariable(glob_mesh,glob_cseq,"junk",&numStep,
		       &baseVar);
  fail_unless(rc == FC_SUCCESS, "failed to create baseVar");
  numDataPoint = 1;
  for (i  = 0; i < numStep; i++){ //like each time
    rc = fc_setVariableData(baseVar[i],numDataPoint,1,
		       FC_AT_WHOLE_MESH,FC_MT_SCALAR,
		       FC_DT_FLOAT,&(glob_ftimes[i]));
    fail_unless(rc == FC_SUCCESS, "failed to set data for baseVar");
  }
  rc = fc_linearLeastSquares(numStep,baseVar,&newVar);
  fail_unless(rc != FC_SUCCESS, "should fail for char sequence");

  //clean up for this test
  data = NULL;
  rc = fc_deleteVariable(newVar);
  fail_unless(rc == FC_SUCCESS, "failed to delete newvar");
  rc = fc_deleteSeqVariable(numStep,baseVar);
  fail_unless(rc == FC_SUCCESS, "failed to delete basevar");
  free(baseVar);

  //check works ok for double sequence but return sea and seb = -1 
  //if num time steps < 3
  numStep = 2; //so cant use global double sequence
  newDSeqCoords = (double*)malloc(numStep*sizeof(double));
  for (j = 0; j < numStep; j++){
    newDSeqCoords[j] = (double)j;
  }
  rc = fc_createSequence(glob_dataset,"new seq",&newSeq); 
  fail_unless(rc == FC_SUCCESS, "failed to create sequence");
  rc = fc_setSequenceCoords(newSeq,numStep,FC_DT_DOUBLE,newDSeqCoords);
  fail_unless(rc == FC_SUCCESS, "failed to set sequence coords");
  free(newDSeqCoords);

  rc = fc_createSeqVariable(glob_mesh,newSeq,"junk",&numStep,
		       &baseVar);
  fail_unless(rc == FC_SUCCESS, "failed to create baseVar");
  numDataPoint = 1;
  for (i  = 0; i < numStep; i++){ //like each time
    rc = fc_setVariableData(baseVar[i],numDataPoint,1,
		       FC_AT_WHOLE_MESH,FC_MT_SCALAR,
		       FC_DT_FLOAT,&(glob_ftimes[i]));
    fail_unless(rc == FC_SUCCESS, "failed to set data for baseVar");
  }
  rc = fc_linearLeastSquares(numStep,baseVar,&newVar);
  fail_unless(rc == FC_SUCCESS, "failed to do linearLeastSquares");
  fail_unless(fc_isVariableValid(newVar),
	      "didn't create valid newvar");
  rc = fc_getVariableDataPtr(newVar,&data);
  fail_unless(rc == FC_SUCCESS, "can't get newVar data");
  {
    double a,b,r2,sea,seb;
    a = ((double*)data)[0];
    b = ((double*)data)[1];
    r2 = ((double*)data)[2];
    sea = ((double*)data)[3];
    seb = ((double*)data)[4];

    //    printf("a %10g b %10g r2 %10g sea %10g seb %10g\n",
    //	   a,b,r2,sea,seb);
    //known values 
    fail_unless(FC_FLT_EQUIV(b,1.0), "bad slope");
    fail_unless(FC_VALUE_EQUIV(a,0.0,FLT_EPSILON,FLT_MIN),
		"bad intercept");
    fail_unless(FC_FLT_EQUIV(r2,1.0), "bad correlation coefficient");
    fail_unless(FC_FLT_EQUIV(sea,-1.0), "bad std error a");
    fail_unless(FC_FLT_EQUIV(seb,-1.0), "bad std error b");
  }
  //clean up for this test
  data = NULL;
  rc = fc_deleteVariable(newVar);
  fail_unless(rc == FC_SUCCESS, "failed to delete newvar");
  rc = fc_deleteSeqVariable(numStep,baseVar);
  fail_unless(rc == FC_SUCCESS, "failed to delete basevar");
  free(baseVar);
  rc = fc_deleteSequence(newSeq);
  fail_unless(rc == FC_SUCCESS, "failed to delete newSeq");


  //check works ok for int sequnce but returns -1 for vertical line
  //cant use glob_iseq since need vertical line
  numStep = glob_numStep;
  rc = fc_createSequence(glob_dataset,"new seq",&newSeq); 
  fail_unless(rc == FC_SUCCESS, "failed to create sequence");
  newISeqCoords = (int*)malloc(numStep*sizeof(int));
  for (j = 0; j < numStep; j++){
    newISeqCoords[j] = 3.2;
  }
  rc = fc_setSequenceCoords(newSeq,numStep,FC_DT_INT,newISeqCoords);
  fail_unless(rc == FC_SUCCESS, "failed to set sequence coords");
  free(newISeqCoords);

  rc = fc_createSeqVariable(glob_mesh,newSeq,"junk",&numStep,
		       &baseVar);
  fail_unless(rc == FC_SUCCESS, "failed to create baseVar");
  numDataPoint = 1;
  for (i  = 0; i < numStep; i++){ //like each time
    rc = fc_setVariableData(baseVar[i],numDataPoint,1,
		       FC_AT_WHOLE_MESH,FC_MT_SCALAR,
		       FC_DT_FLOAT,&(glob_ftimes[i]));
    fail_unless(rc == FC_SUCCESS, "failed to set data for baseVar");
  }

  rc = fc_linearLeastSquares(numStep,baseVar,&newVar);
  fail_unless(rc == FC_SUCCESS, "failed to do linearLeastSquares");
  fail_unless(fc_isVariableValid(newVar),
	      "didn't create valid newvar");
  rc = fc_getVariableDataPtr(newVar,&data);
  fail_unless(rc == FC_SUCCESS, "can't get newVar data");
  {
    double a,b,r2,sea,seb;
    a = ((double*)data)[0];
    b = ((double*)data)[1];
    r2 = ((double*)data)[2];
    sea = ((double*)data)[3];
    seb = ((double*)data)[4];
    
    //    printf("a %10g b %10g r2 %10g sea %10g seb %10g\n",
    //	   a,b,r2,sea,seb);

    //known values ssxx = 0 
    fail_unless(FC_FLT_EQUIV(b,-1.0), "bad slope");
    fail_unless(FC_VALUE_EQUIV(a,-1.0,FLT_EPSILON,FLT_MIN),
		"bad intercept");
    fail_unless(FC_FLT_EQUIV(r2,-1.0), "bad correlation coefficient");
    fail_unless(FC_FLT_EQUIV(sea,-1.0), "bad std error a");
    fail_unless(FC_FLT_EQUIV(seb,-1.0), "bad std error b");
  }
    
  //clean up for this test
  data = NULL;
  rc = fc_deleteVariable(newVar);
  fail_unless(rc == FC_SUCCESS, "failed to delete newvar");
  rc = fc_deleteSeqVariable(numStep,baseVar);
  fail_unless(rc == FC_SUCCESS, "failed to delete basevar");
  free(baseVar);
  rc = fc_deleteSequence(newSeq);
  fail_unless(rc == FC_SUCCESS, "failed to delete newSeq");
}
END_TEST

  
// Populate the Suite with the tests

Suite *series_suite(void)
{
  Suite *suite = suite_create("Series");

  TCase *tc_regtoseq = tcase_create(" - RegToSeq ");
  TCase *tc_window = tcase_create(" - Window ");
  TCase *tc_interpolate = tcase_create(" - LinearInterpolate ");
  TCase *tc_derivatives = tcase_create(" - Derivatives ");
  TCase *tc_squares = tcase_create(" - LeastSquares ");
  TCase *tc_comparison = tcase_create(" - Comparison ");
  TCase *tc_single = tcase_create(" - SingleSeries ");

  //checking this even though it is not in the public interface
  suite_add_tcase(suite, tc_regtoseq);
  tcase_add_checked_fixture(tc_regtoseq, series_setup,
			    series_teardown);
  tcase_add_test(tc_regtoseq, do_regtoseq);

  suite_add_tcase(suite, tc_window);
  tcase_add_checked_fixture(tc_window, series_setup,
			    series_teardown);
  tcase_add_test(tc_window, do_window);

  suite_add_tcase(suite, tc_interpolate);
  tcase_add_checked_fixture(tc_interpolate, series_setup,
			    series_teardown);
  tcase_add_test(tc_interpolate, do_interpolate);

  suite_add_tcase(suite, tc_derivatives);
  tcase_add_checked_fixture(tc_derivatives, series_setup,
			    series_teardown);
  tcase_add_test(tc_derivatives, do_derivatives);

  suite_add_tcase(suite, tc_single);
  tcase_add_checked_fixture(tc_single, series_setup,
			    series_teardown);
  tcase_add_test(tc_single, do_single);

  suite_add_tcase(suite, tc_squares);
  tcase_add_checked_fixture(tc_squares, series_setup,
			    series_teardown);
  tcase_add_test(tc_squares, do_squares);

  //comparison must be after single because the comparison calculation
  //    relies on integral which is tested in single
  suite_add_tcase(suite, tc_comparison);
  tcase_add_checked_fixture(tc_comparison, series_setup,
			    series_teardown);
  tcase_add_test(tc_comparison, do_comparison);


  return suite;
}
