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
 * \file series.c
 * \brief Declarations for \ref Series Module
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/modules/series.c,v $
 * $Revision: 1.103 $
 * $Date: 2006/10/19 03:14:50 $
 *
 */

// C library includes
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include "math.h"

// fc library dependencies
#include "base.h"
#include "storage.h"
#include "library.h"
#include "sequence.h"
#include "mesh.h"
#include "variable.h"
#include "datasetP.h"
#include "variableP.h"
#include "util.h"
#include "varmath.h"
#include "statistics.h"

// this module
#include "series.h"
#include "seriesP.h"


/** 
 * \addtogroup Series
 * \brief Sequence-based analysis routines (e.g., time series) and
 *        helper functions.
 *
 * \description 
 *
 *    Sequence-based analysis routines (e.g., time series)
 *    and helper functions. These may
 *      - compare 2 sequence vars by generating a characteristic value
 *      representing the comparison (e.g., distance between them)
 *      - generate characteristic values of a single sequence var (e.g.,
 *        integral) 
 *      - map a sequence var into another sequence var  (e.g.,
 *        normal form)
 *
 * \modifications
 *      - 6/01/05 ACG this module was renamed from timeseries to series
 */



/** 
 * \ingroup Series
 * \defgroup PrivateSeries (Private)
 */

/** 
 * \ingroup PrivateSeries
 * \brief Makes a new seq var on the mesh that is on the same sequence
 *        as a given variable. this is a helper function.
 *        inputs are not checked
 *        
 *
 * \description
 *
 *        Makes a new seq var on the mesh that is on the same sequence
 *        as a given variable. this is a helper function.
 *        inputs are not checked.
 *
 *        Like fc_createSeqVariable, it actaully allocates the memory
 *        for the array here. Since this uses createseqvar, the user has to
 *        call a deleteSeqvar *and* a free on the seq var when
 *        deallocating it, unless delete gets changed.
 *        
 *        
 *
 * \modifications  
 *    - 1/26/05 ACG, Created
 */
FC_ReturnCode _fc_makeNewSeqVar(
  int numStep, /**< input - numSteps in orig seq */
  FC_Variable *seqvar, /**<input - orig seq var */
  char* varname, /**< output - name of new seq */
  FC_Variable **newseqvar /**< output - new seqvar */
  ){

  FC_ReturnCode rc;
  FC_Mesh mesh;
  FC_Sequence seq;
  int newSteps;

  //not checking input since helper function

  rc = fc_getMeshFromVariable(seqvar[0],&mesh);
  if (rc != FC_SUCCESS){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(rc));
    return rc;

  }

  rc = fc_getSequenceFromSeqVariable(numStep,seqvar,&seq);
  if (rc != FC_SUCCESS)
    return rc;

  //like the implementation of fc_createSeqVariable method,
  //actaully allocating the memory for the return variable here.
  //but i will free it in case of error, if i can
  rc = fc_createSeqVariable(mesh, seq, varname, &newSteps,newseqvar);
  if (rc != FC_SUCCESS){
    fc_deleteSeqVariable(newSteps, *newseqvar); 
    free(*newseqvar);
    *newseqvar=NULL;
    return rc;
  }

  return FC_SUCCESS;

}


/** 
 * \ingroup PrivateSeries
 * \brief Given a regular multi-datapoint variable, make a sequence
 *        variable from it.
 *
 * \description
 *
 *        Given a regular multi-datapoint variable, make a sequence
 *        variable from it. Each data point maps into a different
 *        step in the sequence (lenght of sequnce = num data point)
 *        in the order they are in the regular variable's data array.
 *        Values then are assigned the assocation FC_AT_WHOLE_MESH.
 *
 *        A sequence must be provided such that the length of the
 *        sequence must equal the number of datapoints.
 *        fc_createRegularSequence exists in sequence module for
 *        assistance with that.
 *
 *        Handles multi component data.
 *
 *        NOTE: becuase this uses createSeqVariable the user must call
 *        deleteseqvariable *and* free when deallocating it
 *        unless delete gets changed.
 *
 * \modifications  
 *    - 2/1/05 ACG, Moved here from checkstats.c and renamed
 *     from fc_seqVarFromNormalVar to fc_createSeqVariableFromRegularVariable
 *    - 2/2/05 ACG, changed to require user to pass in sequence, so
 *             that user can create several seq variables on the same 
 *             sequence by this method. seqvars have to be on the same
 *             sequence to use the seqvar comparison routines in this module.
 *             Removed word "regular" from title so as not to confuse
 *             with use of "regular" elsewhere in this module.
 *    -2/10/05 ACG moved to private, since it may or may not be
 *             something we want to keep since ordering isnt necessarily
 *             meaningful in the translations. nevertheless, 
 *             its worth hanging onto for the time being
 *    - 01/09/06 ACG changing/removing handling/cleanup in case of datatype error
 *             at WSK recommendation since it skews the coverage testing. Either
 *             there is already an explict type check or if those lines are reached
 *             you've got serious errors and the cleanup doesnt matter.
 */
FC_ReturnCode _fc_createSeqVariableFromVariable(
  int numStep, /**< input - numstep in the provided seq */
  FC_Sequence seq, /**< input  - handle to the provided sequence */
  FC_Variable normal,  /**< input - the regular variable */
  char* seq_var_name, /**< output - name of the new seq variable */
  FC_Variable** seqvar /**< output - handles to the new seq var */
){

  FC_ReturnCode rc;
  int i;
  FC_Mesh mesh;
  FC_MathType mathtype;
  FC_DataType datatype;
  FC_AssociationType assoc;
  int numComponent, numDataPoint, mystery;

  void *data, *dcounter;


  if (seqvar)
    *seqvar = NULL;

  if ((!fc_isVariableValid(normal)) || (!fc_isSequenceValid(seq))
      || seq_var_name == NULL){ 
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  };


  rc = fc_getVariableInfo (normal,&numDataPoint,&numComponent,
                          &assoc, &mathtype, &datatype);
  if (rc != FC_SUCCESS)
    return rc;

  if (numDataPoint != numStep){
    fc_printfErrorMessage("%s", "invalid sized seq for this regular variable");
    return FC_INPUT_ERROR;
  }

  fc_printfLogMessage("Creating seqvar from var");


  if (datatype != FC_DT_INT && datatype != FC_DT_FLOAT && 
      datatype != FC_DT_DOUBLE &&
      datatype != FC_DT_CHAR){
    fc_printfErrorMessage("Cant handle this datatype");
    return FC_INPUT_ERROR;
  }

  rc = fc_getMeshFromVariable(normal,&mesh);
  if (rc != FC_SUCCESS)
    return rc;

  rc = fc_createSeqVariable(mesh,seq,seq_var_name,
                            &mystery, seqvar);
  if (rc != FC_SUCCESS)
    return rc;

  rc = fc_getVariableDataPtr(normal,&data);
  if (rc != FC_SUCCESS){
    fc_deleteSeqVariable(numStep,*seqvar);
    free(*seqvar);
    *seqvar=NULL;
    return rc;
  }
  for (i = 0; i <numStep; i++){
    switch (datatype){
    case FC_DT_INT:
      dcounter = &(((int*)data)[i*numComponent]);
      break;
    case FC_DT_FLOAT:
      dcounter = &(((float*)data)[i*numComponent]);
      break;
    case FC_DT_DOUBLE:
      dcounter = &(((double*)data)[i*numComponent]);
      break;
    default: //its a char cuase i checked types above
      dcounter = &(((char*)data)[i*numComponent]);
      break;
    }


    rc = fc_setVariableData((*seqvar)[i],1,numComponent,
                            FC_AT_WHOLE_MESH,mathtype,
                            datatype,dcounter);
    if (rc != FC_SUCCESS){
      fc_deleteSeqVariable(numStep,*seqvar);
      free(*seqvar);
      *seqvar=NULL;
      return rc;
    }
  }

  data = NULL;

  return FC_SUCCESS;
}


//real analysis functions

/** \name Functions that map a single sequence var into another single sequence var */
//-------------------------------
//@{


/** 
 * \ingroup Series
 * \brief Return a handle to a seqvar that is the old seq var on a new
 *    sequence. Linear interpolation is done to get the values
 *    at the new sequence points.
 *
 * \description
 *
 *    Returns a handle to a seqvar that is the old seq var on a new
 *    sequence. Linear interpolation is done to get the values
 *    at the new sequence points. Each Data Point and Component
 *    are calculated separately. The return values of the variable
 *    will be Double.
 *
 *    The new and old sequences must be monotonically increasing.
 *
 *    This will return with an error if the sequences and the variable
 *    are not numerical. This will return with an error if the new
 *    sequence range is greater than that of the original sequence.
 *
 *    Like fc_createSeqVariable, it actaully allocates the memory
 *    for the array here - the user must call deleteSeqVariable
 *    *and* free unless delete gets changed.
 *
 * \todo see if this algorithm (or a trivially modified one) will work for
 *    monotonically decreasing sequences. 
 * \todo see about rewriting taking advantage of  fc_getSequenceClosestValue.
 *
 * \modifications  
 *    - 4/12/05 ACG, Created
 *    - 01/09/06 ACG changing/removing handling/cleanup in case of datatype 
 *        error at WSK recommendation since it skews the coverage testing. 
 *        Either there is already an explict type check or if those lines are 
 *        reached you've got serious errors and the cleanup doesnt matter.
 *    - 01/19/06 ACG clean up  
 */
FC_ReturnCode fc_linearInterpolation(
     int numStep,  /**< input - number of steps in the seq var  */
     FC_Variable *seqvar, /**< input - seq var averaging over  */
     FC_Sequence newseq, /**< input - sequence new variable will be over */
     char* varname, /**< input - name for output seq var  */
     FC_Variable **newseqvar /**< output - handles for output seq var  */
){

  FC_ReturnCode rc, rc1, rc2;
  _FC_VarSlot* varSlot;
  int numDataPoint;
  int numComponent;
  FC_Mesh mesh;
  FC_Sequence seq;
  FC_AssociationType assoc;
  FC_MathType mathtype;
  FC_DataType datatype, newdatatype;
  FC_DataType seqdatatype, newseqdatatype;
  int seqnumstep, newseqnumstep;
  void *seqcoords, *newseqcoords;
  double* writedata;
  double seqprev, seqnext, newseqprev, newseqnext;
  int iseq, inewseq;
  int numArrayVals;
  int j,k;


  if (newseqvar)
    *newseqvar = NULL;

  //will check below to see that its appropriate numerical type
  if ((!fc_isSeqVariableValid(numStep,seqvar)) || !newseqvar
      || !varname  || !fc_isSequenceValid(newseq)){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  };

  rc = fc_getVariableDataType(seqvar[0],&datatype);
  if (rc != FC_SUCCESS)
    return rc;

  if (!(datatype == FC_DT_FLOAT ||
        datatype == FC_DT_DOUBLE ||
        datatype == FC_DT_INT )
      ){
    return FC_INPUT_ERROR;
  } 

  rc = fc_getSequenceFromSeqVariable(numStep,seqvar,&seq);
  if (rc != FC_SUCCESS)
    return rc;

  rc1 = fc_getSequenceInfo(seq, &seqnumstep, &seqdatatype);
  rc2 = fc_getSequenceInfo(newseq, &newseqnumstep, &newseqdatatype);
  if (rc1 != FC_SUCCESS || rc2 != FC_SUCCESS) 
    return (rc1 != FC_SUCCESS ? rc1: rc2);     //if both are bad, will get first back. 

  if (seqnumstep < 1 || newseqnumstep < 1){
    return FC_INPUT_ERROR;
  }

  //check datatypes are right
  //already checked that the variable was numerical
  if (!((seqdatatype == FC_DT_FLOAT ||
         seqdatatype == FC_DT_DOUBLE ||
         seqdatatype == FC_DT_INT) &&
        (newseqdatatype == FC_DT_FLOAT ||
         newseqdatatype == FC_DT_DOUBLE ||
         newseqdatatype == FC_DT_INT ))){
    fc_printfErrorMessage("%s", "Invalid datatype");
    return FC_INPUT_ERROR;
  }
  newdatatype = FC_DT_DOUBLE; // return var type will be double

  //lets check monotonicty upfront
  rc1 = fc_getSequenceMinMaxMono(seq,&seqprev,NULL,&seqnext,NULL,&j);
  rc2 = fc_getSequenceMinMaxMono(newseq,&newseqprev,NULL,&newseqnext,NULL,&k);
  if (rc1 !=  FC_SUCCESS || rc2 != FC_SUCCESS)
    return (rc1 != FC_SUCCESS ? rc1: rc2);     //if both are bad, will get first back. 
    
  if (j != 1 && j != 2){
    fc_printfErrorMessage("%s", "Orig sequence not monotonically increasing");
    return FC_INPUT_ERROR;
  }

  if (k != 1 && k !=2){
    fc_printfErrorMessage("%s", "New sequence not monotonically increasing");
    return FC_INPUT_ERROR;
  }

  //now check BC - if the boundaries are off by DBL_epsilon its still ok
  if (newseqnext > seqnext && !FC_DBL_EQUIV(newseqnext,seqnext)){
    fc_printfErrorMessage("%s", "New sequence range end too high");
    return FC_INPUT_ERROR;
  }

  if (newseqprev < seqprev && !FC_DBL_EQUIV(newseqprev,seqprev)){
    fc_printfErrorMessage("%s", "New sequence range begin too low");
    return FC_INPUT_ERROR;
  }

  //finally ok! 
  varSlot = _fc_getVarSlot(seqvar[0]);
  fc_printfLogMessage("Linear interpolations for seq variable '%s'",
                      varSlot->header.name);

  rc = fc_getMeshFromVariable(seqvar[0],&mesh);
  if (rc != FC_SUCCESS){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(rc));
    return rc;
  }

  rc = fc_createSeqVariable(mesh, newseq, varname, &newseqnumstep, newseqvar);
  if (rc != FC_SUCCESS){
    fc_deleteSeqVariable(newseqnumstep,*newseqvar); 
    free(*newseqvar);
    *newseqvar=NULL;
    return rc;
  }

  rc = fc_getVariableInfo(seqvar[0], &numDataPoint, &numComponent, &assoc,
                           &mathtype, &datatype);
  if (rc != FC_SUCCESS){
    fc_deleteSeqVariable(newseqnumstep,*newseqvar); 
    free(*newseqvar);
    *newseqvar=NULL;
    return rc;
  }

  //make writeout arrays
  numArrayVals = numDataPoint*numComponent; 
  writedata = (double*)malloc(numArrayVals*sizeof(double));

  //get data arrays
  rc1 = fc_getSequenceCoordsPtr(seq, &seqcoords);
  rc2 = fc_getSequenceCoordsPtr(newseq, &newseqcoords);
  if (rc1 != FC_SUCCESS || rc2 != FC_SUCCESS){
    fc_deleteSeqVariable(newseqnumstep,*newseqvar); 
    free(*newseqvar);
    *newseqvar=NULL;
    return (rc1 != FC_SUCCESS ? rc1: rc2);     //if both are bad, will get first back. 
  }

  iseq = 0;
  seqnext= seqprev; // get it back to beginning
  for (inewseq = 0; inewseq < newseqnumstep; inewseq++) { //like each time
    //get current value
    newseqprev = newseqnext;
    switch(newseqdatatype){
    case FC_DT_INT:
      newseqnext = (double)((int*)newseqcoords)[inewseq];
      break;
    case FC_DT_FLOAT:
      newseqnext = (double)((float*)newseqcoords)[inewseq];
      break;
    default: // its a double cuase i checked this above
      newseqnext = ((double*)newseqcoords)[inewseq];
      break;
    }


    //find the range for this seq val in the old seq
    while( seqnext < newseqnext  && !FC_DBL_EQUIV(newseqnext,seqnext)){
      iseq++; 
      //know shouldnt overrun seqnumstep, but checking it anyway
      if (iseq > seqnumstep-1){
	fc_printfErrorMessage("Overran sequence bounds. Developer Error!");
	return FC_ERROR;
      }

      seqprev = seqnext;
      switch(seqdatatype){
      case FC_DT_INT:
        seqnext = (double)((int*)seqcoords)[iseq];
        break;
      case FC_DT_FLOAT:
        seqnext = (double)((float*)seqcoords)[iseq];
        break;
      default: //its a double cause i checked this above
        seqnext = ((double*)seqcoords)[iseq];
        break;
      }
    }

    //now we have the boundaries for the val we are looking for
    if (FC_DBL_EQUIV(newseqnext,seqnext)){ //then just copy over val
      double* datanext = NULL;
      rc = fc_getVariableDataAsDataType(seqvar[iseq], FC_DT_DOUBLE,
                                        (void**)(&datanext));
      if (rc != FC_SUCCESS){
        fc_deleteSeqVariable(newseqnumstep,*newseqvar);
        free(*newseqvar);
        *newseqvar=NULL;
        if (datanext) free(datanext);
        if (writedata) free(writedata);
        return rc;
      }

      //copy it over
      rc = fc_setVariableData((*newseqvar)[inewseq], numDataPoint,
                                numComponent, assoc, mathtype,
                                newdatatype,datanext);
      free(datanext);
      if (rc != FC_SUCCESS){
        fc_deleteSeqVariable(newseqnumstep,*newseqvar);
        free(*newseqvar);
        *newseqvar=NULL;
        if (writedata) free(writedata);
        return rc;
      }
    } else {       //interpolate
      double* dataprev = NULL;
      double* datanext = NULL;

      rc1 = fc_getVariableDataAsDataType(seqvar[iseq], FC_DT_DOUBLE,
					(void**)(&datanext));
      rc2 = fc_getVariableDataAsDataType(seqvar[iseq-1], FC_DT_DOUBLE,
					(void**)(&dataprev));
      if (rc1 != FC_SUCCESS || rc2 != FC_SUCCESS){
	fc_deleteSeqVariable(newseqnumstep,*newseqvar);
	free(*newseqvar);
	*newseqvar=NULL;
	if (datanext) free(datanext);
	if (dataprev) free(dataprev);
	if (writedata) free(writedata);
	return (rc1 != FC_SUCCESS ? rc1: rc2);     //if both are bad, will get first back. 
      }

      //filling in variable
      for (j = 0; j < numArrayVals; j++){
	double dt, dvar, diff;
	
	dt = seqnext-seqprev; 
	//dt shouldnt be zero since we checked it with the monotonicity
	if (FC_DBL_EQUIV(dt,0.0)){
	  printf("** Developer Alert! Should never reach this point!**\n");
	  fflush(NULL);
	  fc_printfErrorMessage("%s", "Zero time step ");
	  return FC_ERROR;
	}
	dvar = datanext[j]- dataprev[j];
	diff = newseqnext-seqprev;
	
	if (FC_DBL_EQUIV(diff,0.0) || FC_DBL_EQUIV(dvar,0.0)){
	  writedata[j] = dataprev[j];
	} else {
	  writedata[j] = dvar*(diff/dt)+dataprev[j];
	}
      }
      
      if (dataprev) free(dataprev);
      if (datanext) free(datanext);
        
      rc = fc_setVariableData((*newseqvar)[inewseq], numDataPoint,
			      numComponent, assoc, mathtype,
			      newdatatype,writedata);
      if (rc != FC_SUCCESS){
	fc_deleteSeqVariable(newseqnumstep,*newseqvar);
	free(*newseqvar);
	*newseqvar=NULL;
	if (writedata) free(writedata);
	return rc;
      }
    }
  }
                              
  if (writedata) free(writedata);

  return FC_SUCCESS;

}


/** 
 * \ingroup Series
 * \brief Return a handle to a seqvar that is the first derivative
 *         of an existing seqvar calulated using the Richardson
 *         Extrapolation Approximation 
 *
 * \description
 *
 *    Return a handle to a seqvar that is the first derivative
 *    of an existing seqvar calulated using the Richardson
 *    Extrapolation Approximation.The sequence must be regularly
 *    spaced.
 *
 *    f'(x) = (f(x-2h)-8f(x-h)+8f(x+h)-f(x+2h))/12h.
 *    where h is the regular spacing.
 *
 *    There must be at least 5 steps in the sequence. Return vals
 *    for the first 2 and last 2 pts in the series will be zero.
 *
 *    Does calc for each data point and each component separately.
 *    i.e. original seqvar has (numseriessteps*numdatapoint*numcomponent)
 *    dimensions; the result has numdatapoint*numcomponent) dimensions.
 *
 *    This will return with an error if the seqvar is not numerical.
 *    This will return with an error if the seq is not regularly 
 *    increasingly spaced. User can pass tolerance in. If user passes
 *    in zero it defaults to (type)_EPS
 *
 *    If either the sequence datatype or the value datatype is FC_DT_FLOAT,
 *    then the return value is datatype FC_DT_FLOAT; otherwise it is
 *    FC_DT_DOUBLE;
 *
 *    NOTE: This method is really not recommended for floats becuase of roundoff
 *    error.
 *
 *    NOTE: becuase this allocates memory for seq var, user must call
 *    deleteSeqVar *and* free, unless delete gets changed.
 *
 * \todo
 *    - THIS IS NOT CHEKCED FOR MULTIPLE DP, MULTIPLE COMPONENT
 *    - pass tolerance parameter IS NOT CHECKED
 *    - see if reg spacing criteria works ok for most real life cases.
 *    - can this be made more robust for floats (esp determining irregualr sequence ?)
 *    - is this the methodolgy we want ? (REA)
 *    - check out higher order derivative methods
 *
 * \modifications  
 *    - 07/05/05 ACG, Created
 *    - 07/27/05 ACG added user passes in tolerance
 *    - 01/09/06 ACG changing/removing handling/cleanup in case of datatype error
 *             at WSK recommendation since it skews the coverage testing. Either
 *             there is already an explict type check or if those lines are reached
 *             you've got serious errors and the cleanup doesnt matter.
 *    - 01/12/06 ACG cleaned up, esp wrt writeout arrays
 */

FC_ReturnCode fc_firstDerivative_REA(
     int numStep,  /**< input - number of steps in the seq var  */
     FC_Variable *seqvar, /**< input - seq var differnetiating */
     double ueps, /**< input - user pass in tolerance */
     char* varname, /**< input - name for averaged output var  */
     FC_Variable **newseqvar /**< output - handle for output var  */
     ){
  FC_ReturnCode rc;
  _FC_VarSlot* varSlot;
  FC_Sequence seq;
  int numDataPoint, numComponent, numArrayVals;
  FC_AssociationType assoc;
  FC_MathType mathtype;
  FC_DataType datatype, seqdatatype, newdatatype;
  double dh;
  float fh;
  void *calcdata, *writedata;
  double  ssmax, ssmin, ssmean, sssdev,tol;
  int i,j, count;

  if (newseqvar)
    *newseqvar = NULL;

  //will check below to see that its appropriate numerical type
  if ((!fc_isSeqVariableValid(numStep,seqvar)) ||  !varname || !newseqvar){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  };


  varSlot = _fc_getVarSlot(seqvar[0]);
  fc_printfLogMessage("Getting first derivative for seq variable '%s'",
                      varSlot->header.name);


  rc = fc_getVariableDataType(seqvar[0],&datatype);
  if (rc != FC_SUCCESS) 
    return rc;

  if (!(datatype == FC_DT_FLOAT ||
        datatype == FC_DT_DOUBLE ||
        datatype == FC_DT_INT )
      ){
    fc_printfErrorMessage("%s", "seqvar is not numerical");
    return FC_INPUT_ERROR;
  } 

  if (numStep < 5){
    fc_printfErrorMessage("%s", "seqvar is too short");
    return FC_INPUT_ERROR;
  }


  rc = fc_getSequenceFromSeqVariable(numStep,seqvar,&seq);
  if (rc != FC_SUCCESS) 
    return rc;

  rc = fc_getSequenceDataType(seq,&seqdatatype);
  if (rc != FC_SUCCESS) 
    return rc;

  if (seqdatatype != FC_DT_FLOAT && seqdatatype != FC_DT_DOUBLE &&
      seqdatatype != FC_DT_INT){
    fc_printfErrorMessage("seq data type non numerical");
    return FC_INPUT_ERROR;
  }

  newdatatype = (datatype == FC_DT_FLOAT || seqdatatype == FC_DT_FLOAT) ?
    FC_DT_FLOAT: FC_DT_DOUBLE;


  rc = fc_getSequenceSpacingMinMaxMeanSdev(seq,&ssmin,&i,
                                           &ssmax,&j,&ssmean,&sssdev);
  if (rc != FC_SUCCESS)
    return rc;

  //also need sequence to be monotonically increasing. rather than do mono check
  //since will use mean as spacing, just check that the mean is positive.
  // checks below will make sure that it is regular. these two
  //together will guarentee that the spacing is good enuf
  if (ssmean  < 0){
    fc_printfErrorMessage("%s", "sequence is not monotonically increasing");
    return FC_INPUT_ERROR;
  }
  
  //  printf("\nmin %10.8g max %10.8g  mean %10.8g sdev %10.8g max-min %10.8g max-mean %10.8g
  // min-mean %10.8g sc(max-min) %10.8g sc(max-mean) %10.8g sc(min-mean) %10.8g\n",
  //     ssmin,ssmax,ssmean,sssdev, ssmax-ssmin, ssmax-ssmean, ssmin-ssmean, (ssmax-ssmin)/ssmax,
  //     (ssmax-ssmean)/ssmax, (ssmin-ssmean)/ssmean);



  //check for regular spacing - may have to adjust these constraints later
  if (seqdatatype == FC_DT_FLOAT){
    tol = FLT_EPSILON;
    if (ueps > tol) tol = ueps;
    if (FC_VALUE_EQUIV(ssmin,ssmax,tol,FLT_MIN)  || FC_FLT_EQUIV(sssdev,0.0) ||
        (FC_VALUE_EQUIV(ssmin,ssmean,tol,FLT_MIN) && FC_VALUE_EQUIV(ssmax,ssmean,tol,FLT_MIN))){
      if (FC_FLT_EQUIV(ssmean,0.0)){
        fc_printfErrorMessage("%s", "cant calc with 0 mean spacing");
        return FC_INPUT_ERROR;
      }else {
        fh = ssmean;
      }
    } else{
      fc_printfErrorMessage("%s", "seqvar is not on regular sequence");
      return FC_INPUT_ERROR;
    }
  }else{ //its a double or a int, use dbl conditions in either case
    tol = DBL_EPSILON;
    if (ueps > tol) tol = ueps;
    if (FC_VALUE_EQUIV(ssmin,ssmax,tol,DBL_MIN) || FC_DBL_EQUIV(sssdev,0.0) ||
        (FC_VALUE_EQUIV(ssmin,ssmean,tol,DBL_MIN) && FC_VALUE_EQUIV(ssmax,ssmean,tol,DBL_MIN))){
      if (FC_DBL_EQUIV(ssmean,0.0)){
        fc_printfErrorMessage("%s", "cant calc with 0 mean spacing");
        return FC_INPUT_ERROR;
      }else {
        dh = ssmean;
      }
    } else{
      fc_printfErrorMessage("%s", "seqvar is not on regular sequence");
      return FC_INPUT_ERROR;
    }
  }


  rc = _fc_makeNewSeqVar(numStep, seqvar, varname,newseqvar); 
  if (rc != FC_SUCCESS)
    return rc;

  rc = fc_getVariableInfo(seqvar[0], &numDataPoint, &numComponent, &assoc,
                           &mathtype, &datatype);
  if (rc != FC_SUCCESS){
    fc_deleteSeqVariable(numStep,*newseqvar);
    free(*newseqvar);
    *newseqvar=NULL;
    return rc;
  }


  //first make calc arrays
  numArrayVals = numDataPoint*numComponent;
  if (datatype == FC_DT_FLOAT){
    calcdata = (float*)malloc(numArrayVals*sizeof(float));
  }else{
    calcdata = (double*)malloc(numArrayVals*sizeof(double));
  }
  if (newdatatype == FC_DT_FLOAT){
    writedata = (float*)malloc(numArrayVals*sizeof(float));
  }else{
    writedata = (double*)malloc(numArrayVals*sizeof(double));
  }
  if(calcdata == NULL || writedata == NULL){
    //dont clean up in case of memory error
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }

  if (datatype == FC_DT_FLOAT && newdatatype == FC_DT_FLOAT){
    for (i = 0; i < numArrayVals; i++){
      ((float*)calcdata)[i] = 0.0;
      ((float*)writedata)[i] = 0.0;
    }
  } else if (newdatatype == FC_DT_DOUBLE){ //they both have to be dbls then
    for (i = 0; i < numArrayVals; i++){
      ((double*)calcdata)[i] = 0.0;
      ((double*)writedata)[i] = 0.0;
    }
  } else{ //data must be a dbl and new must be a float
    for (i = 0; i < numArrayVals; i++){
      ((double*)calcdata)[i] = 0.0;
      ((float*)writedata)[i] = 0.0;
    }
  }
      

  //assign zeros to first two and last two steps
  j = 0;
  for (i = 0; i < 4; i++){
    if (i == 2){
      j = numStep - 2;
    }

    rc = fc_setVariableData((*newseqvar)[j], numDataPoint,
			    numComponent, assoc, mathtype,
			    newdatatype,writedata);
    j++;
    if (rc != FC_SUCCESS){
      fc_printfErrorMessage("%s", "cant set var data");
      fc_deleteSeqVariable(numStep,*newseqvar);
      free(*newseqvar);
      *newseqvar=NULL;
      free(writedata);
      free(calcdata);
      return rc;
    }
  }

  //now do all the others
  for (count = 2; count < numStep-2; count++){

    if (datatype == FC_DT_FLOAT && newdatatype == FC_DT_FLOAT){
      for (i = 0; i < numArrayVals; i++){
	((float*)calcdata)[i] = 0.0;
	((float*)writedata)[i] = 0.0;
      }
    } else if (newdatatype == FC_DT_DOUBLE){ //they both have to be dbls then
      for (i = 0; i < numArrayVals; i++){
	((double*)calcdata)[i] = 0.0;
	((double*)writedata)[i] = 0.0;
      }
    } else{ //data must be a dbl and new must be a float
      for (i = 0; i < numArrayVals; i++){
	((double*)calcdata)[i] = 0.0;
	((float*)writedata)[i] = 0.0;
      }
    }

    //now calc vals
    for (j = 0; j < 5; j++){ //for each of the 4 time vals that contrib to this quantity
      void *data;
      int coeff;

      switch(j){ //weighted add to the sum 
      case 0:
	coeff = 1;
	break;
      case 1:
	coeff = -8;
	break;
      case 3: 
	coeff = 8;
	break;
      case 4:
	coeff = -1;
	break;
      default:
	//continue. this actually happens for j = 2
	continue;
      }
      rc = fc_getVariableDataPtr(seqvar[count-2+j], &data); //get a time slice
      for (i = 0; i < numArrayVals; i++){ //for each datapoint
	switch(datatype){
	case FC_DT_FLOAT:
	  ((float*)calcdata)[i]+=coeff*((float*)data)[i];
	  break;
	case FC_DT_INT:
	  ((double*)calcdata)[i]+=coeff*((int*)data)[i];
	  break;
	case FC_DT_DOUBLE:
	  ((double*)calcdata)[i]+=coeff*((double*)data)[i];
	  break;
	default: //this wont happen. do nothing
	  break;
	}
      }
    } //now have finished the 4 vals that contrib to the sum


    //divide by 12h
    if (seqdatatype != FC_DT_FLOAT && datatype != FC_DT_FLOAT){
      for (i = 0; i < numArrayVals; i++){
        if (FC_DBL_EQUIV(((double*)calcdata)[i],0.0)){
	  ((double*)writedata)[i] = 0.0;
        } else {
	  ((double*)writedata)[i] = ((double*)calcdata)[i]/(12.0*dh);
        }
      }
    } else {
      //one or the other or both are a float 
      if (seqdatatype != FC_DT_FLOAT){
	fh = (float)dh;
      }

      if (datatype == FC_DT_FLOAT){
	for (i = 0; i < numArrayVals; i++){
	  float temp = ((float*)calcdata)[i];
	  if (FC_FLT_EQUIV(temp,0.0)){
	    ((float*)writedata)[i] = 0.0;
	  }else{
	    ((float*)writedata)[i]=temp/(12.0*fh);
	  }
	}
      } else{
	for (i = 0; i < numArrayVals; i++){
	  float temp = (float)(((double*)calcdata)[i]);
	  if (FC_FLT_EQUIV(temp,0.0)){
	    ((float*)writedata)[i] = 0.0;
	  }else{
	    ((float*)writedata)[i]=temp/(12.0*fh);
	  }
	}
      }
    }
    
    //assign to this time pt
    rc = fc_setVariableData((*newseqvar)[count], numDataPoint,
			    numComponent, assoc, mathtype,
                            newdatatype,writedata);
    if (rc != FC_SUCCESS){
      fc_deleteSeqVariable(numStep,*newseqvar);
      free(*newseqvar);
      *newseqvar=NULL;
      free(calcdata);
      free(writedata);
      return rc;
    }
  } //now have finished deriv for this time step

  free(calcdata);
  free(writedata);
  return FC_SUCCESS;

}


/** 
 * \ingroup Series
 * \brief Return a handle to a seqvar that is the window averaged
 *         version of an existing seqvar
 *
 * \description
 *
 *    Return a handle to a seqvar that is the window averaged
 *    version of an existing seqvar. Return value is an average
 *    of the current f(t) and the previous (window-1) values of
 *    f. Therefore choosing window = 1, returns the orig sequence.
 *    For early times if there are not enough values to satisfy the window, 
 *    it uses what it can. Time averaged version of this is
 *    \ref fc_leadingWindowAverage_Time.
 *
 *    Averaging is over each dimension of the variable separately
 *
 *    Floats and Doubles will return vars with same type; Ints will
 *    return Doubles. Someone can floor/ceil it later if they
 *    want integer values. This will return with an error if the
 *    seqvar is not numerical.
 *
 * \todo
 *
 *    - Weighted version ?
 *    - Think about making a wrapper and passing in ptrs to function
 *    for the particular calculation parts
 *
 * \modifications  
 *    - 12/21/04 ACG, Created
 *    - 01/09/06 ACG changing/removing handling/cleanup in case of datatype error
 *             at WSK recommendation since it skews the coverage testing. Either
 *             there is already an explict type check or if those lines are reached
 *    - 01/12/06 ACG cleaned up the way the writeout arrays were handled.
 */
FC_ReturnCode fc_leadingWindowAverage(
     int numStep,  /**< input - number of steps in the seq var  */
     FC_Variable *seqvar, /**< input - seq var averaging over  */
     int windowsize, /**< input - window size for average  */
     char* varname, /**< input - name for averaged output seq var  */
     FC_Variable **newseqvar /**< output - handles for output seq var  */
){


  FC_ReturnCode rc;
  _FC_VarSlot* varSlot;
  int numDataPoint;
  int numComponent;
  FC_AssociationType assoc;
  FC_MathType mathtype;
  FC_DataType datatype, newdatatype;
  double* curr;
  void* vwritedata;
  int i,j,count,numArrayVals;

  if (newseqvar)
    *newseqvar = NULL;

  //will check below to see that its appropriate numerical type
  if ((!fc_isSeqVariableValid(numStep,seqvar)) || (newseqvar == NULL)
      || varname == NULL || windowsize < 1){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  };

  rc = fc_getVariableDataType(seqvar[0],&datatype);
  if (rc != FC_SUCCESS)
    return rc;

  if (!(datatype == FC_DT_FLOAT ||
        datatype == FC_DT_DOUBLE ||
        datatype == FC_DT_INT )
      ){
    return FC_INPUT_ERROR;
  } 

  varSlot = _fc_getVarSlot(seqvar[0]);
  fc_printfLogMessage("Getting windowave for seq variable '%s'",
                      varSlot->header.name);

  //special case
  if (windowsize == 1){
    FC_Mesh mesh;
    FC_Sequence seq;

    rc = fc_getMeshFromVariable(seqvar[0],&mesh);
    if (rc != FC_SUCCESS){
      fc_printfErrorMessage("%s", fc_getReturnCodeText(rc));
      return rc;
    }

    rc = fc_getSequenceFromSeqVariable(numStep,seqvar,&seq);
    if (rc != FC_SUCCESS)
      return rc;

    rc = fc_copySeqVariable(numStep,seqvar,mesh,seq,varname,newseqvar); 
    if (rc != FC_SUCCESS)
      return rc;

    return FC_SUCCESS;
  }


  rc = _fc_makeNewSeqVar(numStep, seqvar, varname,newseqvar); 
  if (rc != FC_SUCCESS)
    return rc;

  rc = fc_getVariableInfo(seqvar[0], &numDataPoint, &numComponent, &assoc,
                           &mathtype, &datatype);
  if (rc != FC_SUCCESS){
    fc_deleteSeqVariable(numStep,*newseqvar);
    free(*newseqvar);
    *newseqvar=NULL;
    return rc;
  }


  //make writeout array and calc array
  numArrayVals = numDataPoint*numComponent;
  curr = (double*)malloc(numArrayVals*sizeof(double));
  if (datatype == FC_DT_FLOAT){
    newdatatype = FC_DT_FLOAT;
    vwritedata = (float*)malloc(numArrayVals*sizeof(float));    
  } else { //knows its an it or a double cuase that was checked above)
    newdatatype = FC_DT_DOUBLE;
    vwritedata = (double*)malloc(numArrayVals*sizeof(double));
  }
  if (curr == NULL || vwritedata == NULL){
    //no clean up in case of memory error
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }

  for (i = 0; i < numArrayVals; i++){
    curr[i] = 0.0;
  }

  count = 0;


  for (i = 0; i < numStep; i++) { //like each time
    int dropping = 0;
    void *data, *dropdata;

    if (count == windowsize){
      //if a value is ageing out, get it
      int idx = i - windowsize;
      dropping = 1;
      rc = fc_getVariableDataPtr(seqvar[idx], &dropdata);
      if (rc != FC_SUCCESS){
        data = NULL;
        dropdata = NULL;
        free(curr);
        fc_deleteSeqVariable(numStep,*newseqvar);
        free(*newseqvar);
        *newseqvar=NULL;
	free(vwritedata);
        return rc;
      }
      count--;
    }


    //get new val
    rc = fc_getVariableDataPtr(seqvar[i], &data);
    if (rc != FC_SUCCESS){
      data = NULL;
      dropdata = NULL;
      free(curr);
      fc_deleteSeqVariable(numStep,*newseqvar);
      free(*newseqvar);
      *newseqvar=NULL;
      free(vwritedata);
      return rc;
    }
    count++;

    for (j = 0; j < numArrayVals; j++){
      //if dropping subtract aged value off the curr total.
      //always add new value in.
      double temp;
      switch(datatype){
      case FC_DT_INT:
        temp = (double)((int*)data)[j];
        if (dropping)   temp -= (double)((int*)dropdata)[j];
        break;
      case FC_DT_FLOAT:
        temp = (double)((float*)data)[j];
        if (dropping) temp -= (double)((float*)dropdata)[j];
        break;
      default: //knows its a double cause that was checked above
        temp = ((double*)data)[j];
        if (dropping) temp -= ((double*)dropdata)[j];
        break;
      }

      curr[j] += temp;

      //divide at end so less roundoff; count not always fixed
      //(e.g., at early times)
      if (newdatatype == FC_DT_FLOAT){
        ((float*)vwritedata)[j] = (float)(curr[j]/(double)count);
      }else {
        ((double*)vwritedata)[j] = curr[j]/(double)count;
      }
    }

    rc = fc_setVariableData((*newseqvar)[i], numDataPoint,
                            numComponent, assoc, mathtype,
                            newdatatype,vwritedata);
    if (rc != FC_SUCCESS){
      //additional clean up
      fc_deleteSeqVariable(numStep,*newseqvar);

      free(*newseqvar);
      free(vwritedata);
      *newseqvar=NULL;
      free(curr);
      return rc;
    }
  }



  free(curr);
  free(vwritedata);

  return FC_SUCCESS;
}


/** 
 * \ingroup Series
 * \brief Return a handle to a seqvar that is the time window averaged
 *         version of an existing seqvar.
 *
 * \description
 *
 *    Return a handle to a seqvar that is the time window averaged
 *    version of an existing seqvar.
 *
 *    Return value f(t) = ave(x(t-window)....x(t)), i.e. it is
 *    an average of the current and past values for which the
 *    sequence time is <= (currenttime-timewindow).
 *    For early times, if there are not enough values to satisfy
 *    the window, it uses as many as it can. This is the
 *    window duration version of \ref fc_leadingWindowAverage. There
 *    is nothing that checks that the x variable is time of course,
 *    but it most likely is so I'm using that name. If the
 *    time window is zero, it returns a copy of the orig sequence.
 *
 *    This does not check that the sequence is regular, nor does it
 *    weight the vals in the average by their duration in the
 *    window etc. It is a blind average.
 *
 *    Averaging is over each dimension of the variable separately
 *
 *    Floats and Doubles will return vars with same type; Ints will
 *    return Doubles. Someone can floor/ceil it later if they
 *    want integer values. This will return with an error if the seqvar
 *    or seq is not numerical.
 *
 *
 * \todo
 *    - any interested in a weighted version ?
 *
 * \modifications  
 *    - 01/23/06 ACG, Created
 */
FC_ReturnCode fc_leadingWindowAverage_Time(
     int numStep,  /**< input - number of steps in the seq var  */
     FC_Variable *seqvar, /**< input - seq var averaging over  */
     double windowtime, /**< input - time window size for average  */
     char* varname, /**< input - name for averaged output seq var  */
     FC_Variable **newseqvar /**< output - handles for output seq var  */
){
  FC_ReturnCode rc;
  _FC_VarSlot* varSlot;
  int numDataPoint;
  int numComponent;
  FC_AssociationType assoc;
  FC_MathType mathtype;
  FC_DataType datatype, seqdatatype, newdatatype;
  FC_Sequence seq;

  void *seqcoords, *vwritedata;
  int i,j,rangeindex, numArrayVals;

  if (newseqvar)
    *newseqvar = NULL;

  //will check below to see that its appropriate numerical type
  if ((!fc_isSeqVariableValid(numStep,seqvar)) || (newseqvar == NULL)
      || varname == NULL || windowtime < 0){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  };

  rc = fc_getVariableDataType(seqvar[0],&datatype);
  if (rc != FC_SUCCESS)
    return rc;

  if (!(datatype == FC_DT_FLOAT ||
        datatype == FC_DT_DOUBLE ||
        datatype == FC_DT_INT )
      ){
    fc_printfErrorMessage("var data must be numerical");
    return FC_INPUT_ERROR;
  } 

  rc = fc_getSequenceFromSeqVariable(numStep,seqvar,&seq);
  if (rc != FC_SUCCESS)
    return rc;

  rc = fc_getSequenceDataType(seq,&seqdatatype);
  if (rc != FC_SUCCESS)
    return rc;

  if (!(seqdatatype == FC_DT_FLOAT ||
        seqdatatype == FC_DT_DOUBLE ||
        seqdatatype == FC_DT_INT )
      ){
    fc_printfErrorMessage("seq data must be numerical");
    return FC_INPUT_ERROR;
  } 

  varSlot = _fc_getVarSlot(seqvar[0]);
  fc_printfLogMessage("Getting windowave for seq variable '%s'",
                      varSlot->header.name);


  //special case
  if (FC_DBL_EQUIV(windowtime,0.0)){
    FC_Mesh mesh;
    rc = fc_getMeshFromVariable(seqvar[0],&mesh);
    if (rc != FC_SUCCESS){
      fc_printfErrorMessage("%s", fc_getReturnCodeText(rc));
      return rc;
    }

    rc = fc_copySeqVariable(numStep,seqvar,mesh,seq,varname,newseqvar); 
    if (rc != FC_SUCCESS)
      return rc;

    return FC_SUCCESS;
  }

  rc = fc_getSequenceCoordsPtr(seq,&seqcoords);
  if (rc != FC_SUCCESS)
    return rc;

  rc = _fc_makeNewSeqVar(numStep, seqvar, varname,newseqvar); 
  if (rc != FC_SUCCESS)
    return rc;

  rc = fc_getVariableInfo(seqvar[0], &numDataPoint, &numComponent, &assoc,
                           &mathtype, &datatype);
  if (rc != FC_SUCCESS){
    fc_deleteSeqVariable(numStep,*newseqvar);
    free(*newseqvar);
    *newseqvar=NULL;
    return rc;
  }

  //make writeout array
  numArrayVals = numDataPoint*numComponent;
  if (datatype == FC_DT_FLOAT){
    newdatatype = FC_DT_FLOAT;
    vwritedata = (float*)malloc(numArrayVals*sizeof(float));    
  } else { //knows its an int or a double cause that was checked above)
    newdatatype = FC_DT_DOUBLE;
    vwritedata = (double*)malloc(numArrayVals*sizeof(double));
  }
  if (vwritedata == NULL){
    //no clean up in case of memory error
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }

  rangeindex = 0;
  for (i = 0; i < numStep; i++) { //like each time
    int count = 0;
    double currtime, itertime;

    //get new time
    switch(seqdatatype){
    case FC_DT_INT:
      currtime = ((int*)(seqcoords))[i];
      itertime = ((int*)(seqcoords))[rangeindex];
      break;
    case FC_DT_FLOAT: 
      //we are going to use double for the time no matter
      //what cause its easier
      currtime = ((float*)seqcoords)[i];
      itertime = ((float*)seqcoords)[rangeindex];
      break;
    default: //know its a double
      currtime = ((double*)seqcoords)[i];
      itertime = ((double*)seqcoords)[rangeindex];
    }

    //move rangeindex if necessary
    while ((currtime - itertime) > windowtime){
      rangeindex++;
      switch(seqdatatype){
      case FC_DT_INT:
	itertime = ((int*)seqcoords)[rangeindex];
	break;
      case FC_DT_FLOAT:
	itertime = ((float*)seqcoords)[rangeindex];
	break;
      default: //know its a double
	itertime = ((double*)seqcoords)[rangeindex];
      }
    }

    for (j = 0; j < numArrayVals; j++){
      if (newdatatype == FC_DT_FLOAT){
	((float*)vwritedata)[j] = 0.0;
      }else {
	((double*)vwritedata)[j] = 0.0;
      }
    }

    for (j = rangeindex; j <= i ; j++){
      void *data;
      int k;

      rc = fc_getVariableDataPtr(seqvar[j], &data);
      if (rc != FC_SUCCESS){
	data = NULL;
	fc_deleteSeqVariable(numStep,*newseqvar);
	free(*newseqvar);
	*newseqvar=NULL;
	free(vwritedata);
	return rc;
      }

      for (k = 0; k < numArrayVals; k++){
	switch(datatype){
	case FC_DT_INT:
	  ((double*)vwritedata)[k] += ((int*)data)[k];
	  break;
	case FC_DT_FLOAT:
	  ((float*)vwritedata)[k] += ((float*)data)[k];
	  break;
	default: //knows its a double cause that was checked above
	  ((double*)vwritedata)[k] += ((double*)data)[k];
	  break;
	}
      }
      count++;
    } //end for - over all those in range

    for (j = 0 ; j < numArrayVals; j++){
      if (datatype == FC_DT_FLOAT){
        ((float*)vwritedata)[j] /= (float)count;
      }else {
        ((double*)vwritedata)[j] /= (double)count;
      }
    }

    rc = fc_setVariableData((*newseqvar)[i], numDataPoint,
                            numComponent, assoc, mathtype,
                            newdatatype,vwritedata);
    if (rc != FC_SUCCESS){
      //additional clean up
      fc_deleteSeqVariable(numStep,*newseqvar);

      free(*newseqvar);
      free(vwritedata);
      *newseqvar=NULL;
      return rc;
    }
  } // end for - timestep

  free(vwritedata);

  return FC_SUCCESS;
}


/** 
 * \ingroup Series
 * \brief Return a handle to a seqvar that is the window averaged
 *         version of an existing seqvar
 *
 * \description
 *
 *    Return a handle to a seqvar that is the window averaged
 *    version of an existing seqvar. Return value is an average
 *    of the the window number values of f centered about the
 *    current time. Therefore choosing window = 1 
 *    returns the orig sequence. For early/late times if there are not enough
 *    values to satisfy the window, it uses what it can and preserves
 *    the symettry of the window. Time averaged version of this is
 *    \ref fc_centeredWindowAverage_Time.
 *
 *    Averaging is over each dimension of the variable separately
 *
 *    Floats and Doubles will return vars with same type; Ints will
 *    return Doubles. Someone can floor/ceil it later if they
 *    want integer values. This will return with an error if the
 *    seqvar is not numerical. Windowsize must be odd.
 *
 *
 * \todo
 *
 *    - see notes in \ref fc_leadingWindowAverage
 *
 * \modifications  
 *    - 01/24/06 ACG, Created
 */

FC_ReturnCode fc_centeredWindowAverage(
     int numStep,  /**< input - number of steps in the seq var  */
     FC_Variable *seqvar, /**< input - seq var averaging over  */
     int windowsize, /**< input - window size for average  */
     char* varname, /**< input - name for averaged output seq var  */
     FC_Variable **newseqvar /**< output - handles for output seq var  */
){


  FC_ReturnCode rc;
  _FC_VarSlot* varSlot;
  int numDataPoint;
  int numComponent;
  FC_AssociationType assoc;
  FC_MathType mathtype;
  FC_DataType datatype, newdatatype;
  void* vwritedata;
  int i,j,numArrayVals;

  if (newseqvar)
    *newseqvar = NULL;

  //will check below to see that its appropriate numerical type
  if ((!fc_isSeqVariableValid(numStep,seqvar)) || (newseqvar == NULL)
      || varname == NULL || windowsize < 1 || (windowsize-1)%2 != 0){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  };

  rc = fc_getVariableDataType(seqvar[0],&datatype);
  if (rc != FC_SUCCESS)
    return rc;

  if (!(datatype == FC_DT_FLOAT ||
        datatype == FC_DT_DOUBLE ||
        datatype == FC_DT_INT )
      ){
    return FC_INPUT_ERROR;
  } 


  varSlot = _fc_getVarSlot(seqvar[0]);
  fc_printfLogMessage("Getting windowave for seq variable '%s'",
                      varSlot->header.name);

  //special case
  if (windowsize == 1){
    FC_Mesh mesh;
    FC_Sequence seq;

    rc = fc_getMeshFromVariable(seqvar[0],&mesh);
    if (rc != FC_SUCCESS){
      fc_printfErrorMessage("%s", fc_getReturnCodeText(rc));
      return rc;
    }

    rc = fc_getSequenceFromSeqVariable(numStep,seqvar,&seq);
    if (rc != FC_SUCCESS)
      return rc;

    rc = fc_copySeqVariable(numStep,seqvar,mesh,seq,varname,newseqvar); 
    if (rc != FC_SUCCESS)
      return rc;

    return FC_SUCCESS;
  }


  rc = _fc_makeNewSeqVar(numStep, seqvar, varname,newseqvar); 
  if (rc != FC_SUCCESS)
    return rc;

  rc = fc_getVariableInfo(seqvar[0], &numDataPoint, &numComponent, &assoc,
                           &mathtype, &datatype);
  if (rc != FC_SUCCESS){
    fc_deleteSeqVariable(numStep,*newseqvar);
    free(*newseqvar);
    *newseqvar=NULL;
    return rc;
  }

  //make writeout array
  numArrayVals = numDataPoint*numComponent;
  if (datatype == FC_DT_FLOAT){
    newdatatype = FC_DT_FLOAT;
    vwritedata = (float*)malloc(numArrayVals*sizeof(float));    
  } else { //knows its an int or a double cause that was checked above)
    newdatatype = FC_DT_DOUBLE;
    vwritedata = (double*)malloc(numArrayVals*sizeof(double));
  }
  if (vwritedata == NULL){
    //no clean up in case of memory error
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }


  //this is not the most efficient, but it is the easiest
  for (i = 0; i < numStep; i++) { //like each time
    int count = 0;
    int range = (windowsize-1)/2;
    if (i-range < 0) range = i;
    if (i+range > (numStep-1)) range = numStep-1-i;

    for (j = 0; j < numArrayVals; j++){
      if (newdatatype == FC_DT_FLOAT){
	((float*)vwritedata)[j] = 0.0;
      }else {
	((double*)vwritedata)[j] = 0.0;
      }
    }

    for (j = (i-range); j <= (i+range); j++){
      void *data;
      int k;

      rc = fc_getVariableDataPtr(seqvar[j], &data);
      if (rc != FC_SUCCESS){
	data = NULL;
	fc_deleteSeqVariable(numStep,*newseqvar);
	free(*newseqvar);
	*newseqvar=NULL;
	free(vwritedata);
	return rc;
      }

      for (k = 0; k < numArrayVals; k++){
	switch(datatype){
	case FC_DT_INT:
	  ((double*)vwritedata)[k] += ((int*)data)[k];
	  break;
	case FC_DT_FLOAT:
	  ((float*)vwritedata)[k] += ((float*)data)[k];
	  break;
	default: //knows its a double cause that was checked above
	  ((double*)vwritedata)[k] += ((double*)data)[k];
	  break;
	}
      }
      count++;
    } //end for - over all those in range

    for (j = 0 ; j < numArrayVals; j++){
      if (datatype == FC_DT_FLOAT){
        ((float*)vwritedata)[j] /= (float)count;
      }else {
        ((double*)vwritedata)[j] /= (double)count;
      }
    }

    rc = fc_setVariableData((*newseqvar)[i], numDataPoint,
                            numComponent, assoc, mathtype,
                            newdatatype,vwritedata);
    if (rc != FC_SUCCESS){
      //additional clean up
      fc_deleteSeqVariable(numStep,*newseqvar);

      free(*newseqvar);
      free(vwritedata);
      *newseqvar=NULL;
      return rc;
    }
  } // end for - timestep

  free(vwritedata);

  return FC_SUCCESS;
}


/** 
 * \ingroup Series
 * \brief Return a handle to a seqvar that is the time window averaged
 *         version of an existing seqvar.
 *
 * \description
 *
 *    Return a handle to a seqvar that is the time window averaged
 *    version of an existing seqvar.
 *
 *    Return value f(t) = ave(x(t-window/2)....x(t+window/2)), i.e. it is
 *    an average of the vals of f in the time window centered about
 *    the current time. For early/late times, if there are not enough
 *    values to satisfy the window, it uses as many as it can and
 *    preserves the symmetry of the window. This is the
 *    window duration version of \ref fc_centeredWindowAverage. There
 *    is nothing that checks that the x variable is time of course,
 *    but it most likely is so I'm using that name. If the
 *    time window is zero, it returns a copy of the orig sequence.
 *
 *    This does not check that the sequence is regular, nor does it
 *    weight the vals in the average by their duration in the
 *    window etc. It is a blind average.
 *
 *    Averaging is over each dimension of the variable separately
 *
 *    Floats and Doubles will return vars with same type; Ints will
 *    return Doubles. Someone can floor/ceil it later if they
 *    want integer values. This will return with an error if the seqvar
 *    or seq is not numerical.
 *
 *
 * \todo
 *    - any interested in a weighted version ?
 *
 * \modifications  
 *    - 01/24/06 ACG, Created
 */

FC_ReturnCode fc_centeredWindowAverage_Time(
     int numStep,  /**< input - number of steps in the seq var  */
     FC_Variable *seqvar, /**< input - seq var averaging over  */
     double windowtime, /**< input - time window size for average  */
     char* varname, /**< input - name for averaged output seq var  */
     FC_Variable **newseqvar /**< output - handles for output seq var  */
){
  FC_ReturnCode rc;
  _FC_VarSlot* varSlot;
  int numDataPoint;
  int numComponent;
  FC_Sequence seq;
  FC_AssociationType assoc;
  FC_MathType mathtype;
  FC_DataType datatype, seqdatatype,newdatatype;
  void *vwritedata, *seqcoords;
  int i,j,numArrayVals;

  if (newseqvar)
    *newseqvar = NULL;

  //will check below to see that its appropriate numerical type
  if ((!fc_isSeqVariableValid(numStep,seqvar)) || (newseqvar == NULL)
      || varname == NULL || windowtime < 0){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  };

  rc = fc_getVariableDataType(seqvar[0],&datatype);
  if (rc != FC_SUCCESS)
    return rc;

  if (!(datatype == FC_DT_FLOAT ||
        datatype == FC_DT_DOUBLE ||
        datatype == FC_DT_INT )
      ){
    fc_printfErrorMessage("var data must be numerical");
    return FC_INPUT_ERROR;
  } 

  rc = fc_getSequenceFromSeqVariable(numStep,seqvar,&seq);
  if (rc != FC_SUCCESS)
    return rc;

  rc = fc_getSequenceDataType(seq,&seqdatatype);
  if (rc != FC_SUCCESS)
    return rc;

  if (!(seqdatatype == FC_DT_FLOAT ||
        seqdatatype == FC_DT_DOUBLE ||
        seqdatatype == FC_DT_INT )
      ){
    fc_printfErrorMessage("seq data must be numerical");
    return FC_INPUT_ERROR;
  } 


  varSlot = _fc_getVarSlot(seqvar[0]);
  fc_printfLogMessage("Getting windowave for seq variable '%s'",
                      varSlot->header.name);


  //special case
  if (FC_DBL_EQUIV(windowtime,0.0)){
    FC_Mesh mesh;
    rc = fc_getMeshFromVariable(seqvar[0],&mesh);
    if (rc != FC_SUCCESS){
      fc_printfErrorMessage("%s", fc_getReturnCodeText(rc));
      return rc;
    }

    rc = fc_copySeqVariable(numStep,seqvar,mesh,seq,varname,newseqvar); 
    if (rc != FC_SUCCESS)
      return rc;

    return FC_SUCCESS;
  }

  rc = fc_getSequenceCoordsPtr(seq,&seqcoords);
  if (rc != FC_SUCCESS)
    return rc;

  rc = _fc_makeNewSeqVar(numStep, seqvar, varname,newseqvar); 
  if (rc != FC_SUCCESS)
    return rc;

  rc = fc_getVariableInfo(seqvar[0], &numDataPoint, &numComponent, &assoc,
                           &mathtype, &datatype);
  if (rc != FC_SUCCESS){
    fc_deleteSeqVariable(numStep,*newseqvar);
    free(*newseqvar);
    *newseqvar=NULL;
    return rc;
  }

  //make writeout array
  numArrayVals = numDataPoint*numComponent;
  if (datatype == FC_DT_FLOAT){
    newdatatype = FC_DT_FLOAT;
    vwritedata = (float*)malloc(numArrayVals*sizeof(float));    
  } else { //knows its an int or a double cause that was checked above)
    newdatatype = FC_DT_DOUBLE;
    vwritedata = (double*)malloc(numArrayVals*sizeof(double));
  }
  if (vwritedata == NULL){
    //no clean up in case of memory error
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }

  for (i = 0; i < numStep; i++) { //like each time
    double currtime, lowtime, hightime, prevtime;
    double timebounds,newtimebounds;
    int highindex, lowindex, previndex;
    int count;

    timebounds = windowtime/2.0;
    //get new time - currtime is alwasy a double cause its easier
    switch(seqdatatype){
    case FC_DT_INT:
      currtime = ((int*)(seqcoords))[i];
      break;
    case FC_DT_FLOAT:
      currtime = ((float*)(seqcoords))[i];
      break;
    default: //know its a double
      currtime = ((double*)(seqcoords))[i];
      break;
    }

    //find time range indicies
    lowtime = currtime;
    lowindex = i;
    prevtime = currtime;
    previndex = i;
    while((currtime-lowtime) <= timebounds){
      previndex = lowindex;
      prevtime = lowtime;

      if (lowindex == 0) break;
      lowindex--;
      switch(seqdatatype){
      case FC_DT_INT:
	lowtime = ((int*)(seqcoords))[lowindex];
	break;
      case FC_DT_FLOAT:
	lowtime = ((float*)(seqcoords))[lowindex];
	break;
      default: //know its a double
	lowtime = ((double*)(seqcoords))[lowindex];
	break;
      }
    }
    lowindex = previndex;
    lowtime = prevtime;

    hightime = currtime;
    highindex = i;
    prevtime = currtime;
    previndex = i;
    while((hightime-currtime) <= timebounds){
      previndex = highindex;
      prevtime = hightime;

      if (highindex == numStep -1) break;
      highindex++;
      switch(seqdatatype){
      case FC_DT_INT:
	hightime = ((int*)(seqcoords))[highindex];
	break;
      case FC_DT_FLOAT:
	hightime = ((float*)(seqcoords))[highindex];
	break;
      default: //know its a double
	hightime = ((double*)(seqcoords))[highindex];
	break;
      }
    }

    highindex = previndex;
    hightime = prevtime;

    //now we have valid high and low times and indicies, but
    //did we stop because we hit the bounds ?
    newtimebounds = timebounds;
    if (lowindex == 0){
      newtimebounds = currtime-lowtime;
    }
    if (highindex == (numStep-1)){
      newtimebounds = (newtimebounds < (hightime-currtime) ? 
		       newtimebounds : (hightime-currtime));
      
    }
    if (!FC_DBL_EQUIV(newtimebounds,timebounds) &&
	newtimebounds < timebounds){ //redo it
      //latter check just incase there is some weird round off thing
      timebounds = newtimebounds;

      //find time range indicies
      lowtime = currtime;
      lowindex = i;
      prevtime = currtime;
      previndex = i;
      while((currtime-lowtime) <= timebounds){
	previndex = lowindex;
	prevtime = lowtime;
	
	if (lowindex == 0) break;
	lowindex--;
	switch(seqdatatype){
	case FC_DT_INT:
	  lowtime = ((int*)(seqcoords))[lowindex];
	  break;
	case FC_DT_FLOAT:
	  lowtime = ((float*)(seqcoords))[lowindex];
	  break;
	default: //know its a double
	  lowtime = ((double*)(seqcoords))[lowindex];
	  break;
	}
      }
      lowindex = previndex;
      lowtime = prevtime;

      hightime = currtime;
      highindex = i;
      prevtime = currtime;
      previndex = i;
      while((hightime-currtime) <= timebounds){
	previndex = highindex;
	prevtime = hightime;
	
	if (highindex == numStep -1) break;
	highindex++;
	switch(seqdatatype){
	case FC_DT_INT:
	  hightime = ((int*)(seqcoords))[highindex];
	  break;
	case FC_DT_FLOAT:
	  hightime = ((float*)(seqcoords))[highindex];
	  break;
	default: //know its a double
	  hightime = ((double*)(seqcoords))[highindex];
	  break;
	}
      }

      highindex = previndex;
      hightime = prevtime;
    }

    //now weve got valid indicie bounds for this step
    for (j = 0; j < numArrayVals; j++){
      if (newdatatype == FC_DT_FLOAT){
	((float*)vwritedata)[j] = 0.0;
      }else {
	((double*)vwritedata)[j] = 0.0;
      }
    }

    count = 0;
    for (j = lowindex; j <= highindex ; j++){
      void *data;
      int k;

      rc = fc_getVariableDataPtr(seqvar[j], &data);
      if (rc != FC_SUCCESS){
	data = NULL;
	fc_deleteSeqVariable(numStep,*newseqvar);
	free(*newseqvar);
	*newseqvar=NULL;
	free(vwritedata);
	return rc;
      }

      for (k = 0; k < numArrayVals; k++){
	switch(datatype){
	case FC_DT_INT:
	  ((double*)vwritedata)[k] += ((int*)data)[k];
	  break;
	case FC_DT_FLOAT:
	  ((float*)vwritedata)[k] += ((float*)data)[k];
	  break;
	default: //knows its a double cause that was checked above
	  ((double*)vwritedata)[k] += ((double*)data)[k];
	  break;
	}
      }
      count++;
    } //end for - over all those in range

    for (j = 0 ; j < numArrayVals; j++){
      if (datatype == FC_DT_FLOAT){
        ((float*)vwritedata)[j] /= (float)count;
      }else {
        ((double*)vwritedata)[j] /= (double)count;
      }
    }

    rc = fc_setVariableData((*newseqvar)[i], numDataPoint,
                            numComponent, assoc, mathtype,
                            newdatatype,vwritedata);
    if (rc != FC_SUCCESS){
      //additional clean up
      fc_deleteSeqVariable(numStep,*newseqvar);

      free(*newseqvar);
      free(vwritedata);
      *newseqvar=NULL;
      return rc;
    }
  } // end for - timestep

  free(vwritedata);

  return FC_SUCCESS;
}


/** 
 * \ingroup Series
 * \brief Return a handle to a seqvar that is the normal form of
 *         an existing seqvar
 *
 * \description
 *
 *    Normal forms are done over each dimension and component
 *    separately. Normal form X(t)=  (x(t)-xmean(t))/sdev(x);
 *
 *    This will return with an error if the seqvar is not numerical.
 *    This will return the vals unchanged if the sdev for a sequence
 *    is zero.
 *
 *    Floats and Doubles will return vars with same type; Ints will
 *    return Doubles. Someone can floor/ceil it later if they
 *    want integer values.
 *
 *    NOTE: becuase this allocates memory for seq var, user must call
 *    deleteSeqVar *and* free, unless delete gets changed.
 *
 * \todo
 *    - ask new seqvar functions are written for varmath see if can
 *      replace any of the innards of this
 *    - Decide if want to change return if stdev is zero to error.
 *    Hesitating on that since want to get the legitimate ones back.
 *
 *
 *
 * \modifications  
 *    - 1/21/04 ACG, Created
 *    - 01/09/06 ACG changing/removing handling/cleanup in case of datatype error
 *             at WSK recommendation since it skews the coverage testing. Either
 *             there is already an explict type check or if those lines are reached
 *             you've got serious errors and the cleanup doesnt matter.
 *    - 01/12/06 ACG cleaned up the way the writeout arrays were handled.
 */
        
FC_ReturnCode fc_normalForm(
     int numStep,  /**< input - number of steps in the seq var  */
     FC_Variable *seqvar, /**< input - seq var changing to normal form  */
     char* varname, /**< input - name for averaged output seq var  */
     FC_Variable **newseqvar /**< output - handles for output seq var  */
){


  FC_ReturnCode rc;
  _FC_VarSlot* varSlot;
  int numDataPoint;
  int numComponent;
  FC_AssociationType assoc;
  FC_MathType mathtype;
  FC_DataType datatype, newdatatype;
  FC_Variable meanvar, sdevvar;
  void* vwritedata;
  int i,j,numArrayVals;
  double *meandata, *sdevdata;


  if (newseqvar)
    *newseqvar = NULL;

  //will check below to see that its appropriate numerical type
  if ((!fc_isSeqVariableValid(numStep,seqvar)) || (newseqvar == NULL)
      || varname == NULL){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  };

  rc = fc_getVariableDataType(seqvar[0],&datatype);
  if (rc != FC_SUCCESS)
    return rc;

  if (!(datatype == FC_DT_FLOAT ||
        datatype == FC_DT_DOUBLE ||
        datatype == FC_DT_INT )
      ){
    return FC_INPUT_ERROR;
  } 


  varSlot = _fc_getVarSlot(seqvar[0]);
  fc_printfLogMessage("Getting normal form seq variable '%s'",
                      varSlot->header.name);


  rc = fc_getSeqVariableSeriesMeanSdev(numStep,seqvar,&meanvar,&sdevvar);
  if (rc != FC_SUCCESS){
    fc_printfErrorMessage("%s","can't get mean and stdev for this variable");
    return FC_INPUT_ERROR;
  }

  //going to do the cast on these so dont have to allocate arrays for them
  rc = fc_getVariableDataPtr(meanvar,(void**)&meandata);
  if (rc != FC_SUCCESS){
    fc_deleteVariable(meanvar);
    fc_deleteVariable(sdevvar);
    return rc;
  }

  rc = fc_getVariableDataPtr(sdevvar,(void**)&sdevdata);
  if (rc != FC_SUCCESS){
    fc_deleteVariable(meanvar);
    fc_deleteVariable(sdevvar);
    return rc;
  }


  rc = _fc_makeNewSeqVar(numStep, seqvar, varname,newseqvar); 
  if (rc != FC_SUCCESS){
    fc_deleteVariable(meanvar);
    fc_deleteVariable(sdevvar);
    return rc;
  }

  rc = fc_getVariableInfo(seqvar[0], &numDataPoint, &numComponent, &assoc,
                           &mathtype, &datatype);
  if (rc != FC_SUCCESS){
    fc_deleteSeqVariable(numStep,*newseqvar);
    free(*newseqvar);
    *newseqvar=NULL;
    fc_deleteVariable(meanvar);
    fc_deleteVariable(sdevvar);
    return rc;
  }


  //make writeout array
  numArrayVals = numDataPoint*numComponent;
  if (datatype == FC_DT_FLOAT){
    newdatatype = FC_DT_FLOAT;
    vwritedata = (float*)malloc(numArrayVals*sizeof(float));    
  }else {     //its an int or a double cuase i checked it above
    newdatatype = FC_DT_DOUBLE;
    vwritedata = (double*)malloc(numArrayVals*sizeof(double));
  }
  if (vwritedata == NULL ){
    //no clean up in case of memory error
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }

  for (i = 0; i < numStep; i++) { //like each time
    double* data;
    rc = fc_getVariableDataAsDataType(seqvar[i], FC_DT_DOUBLE,
                                      (void**)(&data));
    if (rc != FC_SUCCESS){
      fc_deleteSeqVariable(numStep,*newseqvar);
      free(*newseqvar);
      *newseqvar=NULL;
      if (data) free(data);
      free(vwritedata);
      fc_deleteVariable(meanvar);
      fc_deleteVariable(sdevvar);
      return rc;
    }

    for (j = 0; j < numArrayVals; j++){
      if (newdatatype == FC_DT_FLOAT){
        ((float*)vwritedata)[j] = (FC_DBL_EQUIV(sdevdata[j],0.0) ?
                         (float)(data[j]):
                         (float)(data[j] - meandata[j])/sdevdata[j]);
      }else {
        ((double*)vwritedata)[j] = (FC_DBL_EQUIV(sdevdata[j],0.0) ?
                         data[j]:(data[j] - meandata[j])/sdevdata[j]);
      }
    }


    free(data);
    rc = fc_setVariableData((*newseqvar)[i], numDataPoint,
                            numComponent, assoc, mathtype,
                            newdatatype,vwritedata);


    if (rc != FC_SUCCESS){
      //addtional cleanup
      free(vwritedata);
      fc_deleteVariable(meanvar);
      fc_deleteVariable(sdevvar);
      fc_deleteSeqVariable(numStep,*newseqvar);
      free(*newseqvar);
      *newseqvar=NULL;
      return rc;
    }
  }

  free(vwritedata);
  fc_deleteVariable(meanvar);
  fc_deleteVariable(sdevvar);

  return FC_SUCCESS;
}


//@}


/** \name Functions that operate on a single sequence var and return a value or values based on that var but the return val is not a sequence var */
//-------------------------------
//@{


/** 
 * \ingroup Series
 * \brief Return a handle to a var that is the integral of an
 *        existing seqvar over its sequence.
 *
 * \description
 *
 *    Does calc for each data point and each component separately.
 *    i.e. original seqvar has (numseriessteps*numdatapoint*numcomponent)
 *    dimensions; the result has numdatapoint*numcomponent) dimensions.
 *    Uses the trapezoidal rule for the calculation. The sequence
 *    does not have to be equally spaced.
 *
 *    This will return with an error if the seqvar is not numerical.
 *    This will return with an error if the seq is not increasing.
 *    Must be at least 2 steps in the sequence.
 *
 *    Value is returned as datatype FC_DT_DOUBLE.
 *
 * \todo
 *  - should floats return floats ?
 *  - is it worth making a version that takes advantage of equally spaced
 *  sequences ?
 *
 * \modifications  
 *    - 2/22/05 ACG, Created
 *    - 3/14/05 ACG, tested
 *    - 7/06/05 ACG renamed for consistency with derviative
 *    - 01/09/06 ACG changing/removing handling/cleanup in case of datatype error
 *             at WSK recommendation since it skews the coverage testing. Either
 *             there is already an explict type check or if those lines are reached
 *             you've got serious erros and the cleanup doesnt matter.
 */
FC_ReturnCode fc_integral_TR(
     int numStep,  /**< input - number of steps in the seq var  */
     FC_Variable *seqvar, /**< input - seq var integrating */
     char* varname, /**< input - name for averaged output var  */
     FC_Variable *newvar /**< output - handle for output var  */
     ){
  FC_ReturnCode rc;
  _FC_VarSlot* varSlot;
  FC_Mesh mesh;
  FC_Sequence seq;
  int numDataPoint;
  int numComponent;
  FC_AssociationType assoc;
  FC_MathType mathtype;
  FC_DataType datatype,dt;
  int varsize;
  int i,j;
  void *data;
  double *currtot, *prevval, *nextval, deltat;
  void* seqcoords;

  if (newvar)
    *newvar = FC_NULL_VARIABLE;

  //check input
  if (newvar == NULL || varname == NULL ||
      (!fc_isSeqVariableValid(numStep,seqvar))){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }


  varSlot = _fc_getVarSlot(seqvar[0]);
  fc_printfLogMessage("Getting integral for seq variable '%s'",
                      varSlot->header.name);


  rc = fc_getVariableInfo(seqvar[0], &numDataPoint, &numComponent, &assoc,
                           &mathtype, &datatype);
  if (rc != FC_SUCCESS) return rc;

  if (!(datatype == FC_DT_FLOAT ||
        datatype == FC_DT_DOUBLE ||
        datatype == FC_DT_INT )
      ){
    return FC_INPUT_ERROR;
  } 


  if (numStep < 2){
    fc_printfErrorMessage("sequence doesnt have enough steps to integrate");
    return FC_INPUT_ERROR;
  }

  rc = fc_getSequenceFromSeqVariable(numStep,seqvar,&seq);
  if (rc != FC_SUCCESS)
    return rc;

  rc = fc_getSequenceCoordsPtr(seq,&seqcoords);
  if (rc != FC_SUCCESS)
    return rc;

  rc = fc_getSequenceDataType(seq,&dt);
  if (rc != FC_SUCCESS) 
    return rc;

  if (!(dt == FC_DT_FLOAT ||
        dt == FC_DT_DOUBLE ||
        dt == FC_DT_INT )
      ){
    return FC_INPUT_ERROR;
  } 


  // set up data arrays
  varsize = numDataPoint*numComponent;
  currtot = (double*)malloc(varsize*sizeof(double));
  nextval = (double*)malloc(varsize*sizeof(double));
  prevval = (double*)malloc(varsize*sizeof(double));
  if (currtot == NULL || nextval == NULL || prevval == NULL){
    //no clean up in case of memory error
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }

  //get first point 
  rc = fc_getVariableDataPtr(seqvar[0], &data);
  if (rc != FC_SUCCESS){
    data = NULL;
    free(currtot);
    free(nextval);
    free(prevval);
    return rc;
  }
  for (j = 0; j < varsize; j++){
    currtot[j] = 0;
    switch(datatype){
    case FC_DT_DOUBLE:
      prevval[j]= ((double*)data)[j];
      break;
    case FC_DT_INT:
      prevval[j]= (double)(((int*)data)[j]);
      break;
    default: //its a float cuase i checked it above
      prevval[j]= (double)(((float*)data)[j]);
      break;
    }
  }

  //do rest of them....
  for (i = 1; i < numStep; i++){
    switch (dt){
    case FC_DT_DOUBLE:
      deltat= ((double*)seqcoords)[i]-((double*)seqcoords)[i-1];
      break;
    case FC_DT_INT:
      deltat= (double)((int*)seqcoords)[i]-(double)((int*)seqcoords)[i-1];
      break;
    default: //its a float cuase i checked it above
      deltat= (double)((float*)seqcoords)[i]-
        (double)((float*)seqcoords)[i-1];
      break;
    }
    
    if (FC_DBL_EQUIV(deltat,0.0)){
      fc_printfErrorMessage("divide by zero");
      free(currtot);
      free(nextval);
      free(prevval);
      return FC_ERROR;
    }
    if (deltat < 0.0){
      fc_printfErrorMessage("sequence not increasing");
      free(currtot);
      free(nextval);
      free(prevval);
      return FC_ERROR; 
    }


    rc = fc_getVariableDataPtr(seqvar[i], &data);
    if (rc != FC_SUCCESS){
      data = NULL;
      free(currtot);
      free(nextval);
      free(prevval);
      return rc;
    }
    for (j = 0; j < varsize; j++){
      switch(datatype){
      case FC_DT_DOUBLE:
        nextval[j] = ((double*)data)[j];
        break;
      case FC_DT_INT:
        nextval[j] = (double)(((int*)data)[j]);
        break;
      default: //its a float cuase i checked this above
        nextval[j] = (double)(((float*)data)[j]);
        break;
      }
      currtot[j] += deltat*(nextval[j]+prevval[j])/2.0;
      prevval[j] = nextval[j];
    }
  }

  free(nextval);
  free(prevval);

  rc = fc_getMeshFromVariable(seqvar[0],&mesh);
  if (rc != FC_SUCCESS){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(rc));
    return rc;
  }

  rc = fc_createVariable(mesh,varname,newvar);
  if (rc != FC_SUCCESS){
    free(currtot);
    return rc;
  }

  rc = fc_setVariableDataPtr(*newvar,numDataPoint, numComponent,
                             assoc,mathtype,FC_DT_DOUBLE,currtot);
  if (rc != FC_SUCCESS){
    fc_deleteVariable(*newvar);
    return rc;
  }

  return FC_SUCCESS;
}


/**
 * \ingroup Series
 * \brief Given an N Data Point, Single Component sequence variable, 
 *  returns a handle to a var containing linear least squares parameters,
 *  a,b,r2,SE(a),SE(b), where y=a+bx. 
 *
 * \description
 *
 *    Given an N Data Point, Single Component sequence variable,
 *    return a handle to an N DataPoint, 5 component non-sequence
 *    variable containing least squares parameters
 *    a,b,r2,SE(a),SE(b), where y=a+bx. Therefore, the least-squares
 *    fit is doen for each data point's series independently
 *    The variable must be single component and it is understood
 *    that x=seq values and y = variable values. 
 *
 *    The output variable name is the name of the original variable,
 *    postpended with "_LSQ" with values in the order given above.
 *    It will be a Double.
 *
 *    This will return with an error if the seqvars are not numerical.
 *    if ssxx == 0, all values will be -1.
 *    if ssyy == 0, all values are calc, except r2 ==-1
 *    if numStep < 3, all values are calc, except SE(a) and SE(b) == -1
 *
 * \todo
 *    see about relaxing the single component constraint.
 *
 * \modifications  
 *    - 3/21/05 ACG, Created
 *    - 4/06/05 ACG, changed from calculating the sumofsquares to 
 *                   the actual linear least squares fit. will worry about
 *                   what to do for non-linear fits later
 *    - 01/09/06 ACG changing/removing handling/cleanup in case of datatype error
 *             at WSK recommendation since it skews the coverage testing. Either
 *             there is already an explict type check or if those lines are reached
 *             you've got serious erros and the cleanup doesnt matter.
 */
FC_ReturnCode fc_linearLeastSquares(
     int numStep,  /**< input - number of steps in the seq var  */
     FC_Variable *seqvar, /**< input - seq var */
     FC_Variable *newvar /**< output - handle for output var  */
     ){
  FC_ReturnCode rc;
  _FC_VarSlot* varSlot;
  FC_Mesh mesh;
  int numDataPoint;
  int numComponent;
  FC_AssociationType assoc;
  FC_MathType mathtype;
  FC_DataType datatype, dt;
  FC_Sequence seq;
  int sum_dim;
  int i,j;
  void* data;
  double *sumx, *sumy, *sumzz, *retdata;
  char *varname, *returnname;
  void* seqcoords;


  if (newvar)
    *newvar = FC_NULL_VARIABLE;

  //check input
  if (newvar == NULL ||
      (!fc_isSeqVariableValid(numStep,seqvar))){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }


  varSlot = _fc_getVarSlot(seqvar[0]);
  fc_printfLogMessage("Getting sumofSquares for seq variable '%s'",
                      varSlot->header.name);


  rc = fc_getVariableInfo(seqvar[0], &numDataPoint, &numComponent, &assoc,
                           &mathtype, &datatype);
  if (rc != FC_SUCCESS) return rc;

  if (!(datatype == FC_DT_FLOAT ||
        datatype == FC_DT_DOUBLE ||
        datatype == FC_DT_INT )
      ){
    return FC_INPUT_ERROR;
  } 

  if (numComponent != 1){
    return FC_INPUT_ERROR;
  }

  rc = fc_getSequenceFromSeqVariable(numStep,seqvar,&seq);
  if (rc != FC_SUCCESS)
    return rc;

  rc = fc_getSequenceCoordsPtr(seq,&seqcoords);
  if (rc != FC_SUCCESS)
    return rc;

  rc = fc_getSequenceDataType(seq,&dt);
  if (rc != FC_SUCCESS) 
    return rc;

  if (!(dt == FC_DT_FLOAT ||
        dt == FC_DT_DOUBLE ||
        dt == FC_DT_INT )
      ){
    return FC_INPUT_ERROR;
  } 

  sum_dim = numDataPoint*numComponent; //numcomponent==1
  sumx = (double*)malloc(sum_dim*sizeof(double));
  sumy = (double*)malloc(sum_dim*sizeof(double));
  sumzz = (double*)malloc(3*sum_dim*sizeof(double));
  retdata= (double*)malloc(5*sum_dim*sizeof(double));
  if (sumx == NULL || sumy == NULL || sumzz == NULL || retdata == NULL){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }

  for (i =0; i < sum_dim; i++){
    sumx[i] = 0;
    sumy[i] = 0;
    for (j = 0; j < 3; j++){
      sumzz[3*i+j] = 0.0;
    }
    for (j = 0; j < 5; j++){
      retdata[5*i+j] = 0.0;
    }
  }


  for (i = 0; i < numStep; i++){
    double currx, curry;
    switch(dt){
    case FC_DT_DOUBLE:
      currx = ((double*)seqcoords)[i];
      break;
    case FC_DT_INT:
      currx = (double)(((int*)seqcoords)[i]);
      break;
    default: // know its a float cause i checked it above
      currx = (double)(((float*)seqcoords)[i]);
      break;
    }

    rc = fc_getVariableDataPtr(seqvar[i], &data);
    if (rc != FC_SUCCESS){
      data = NULL;
      free(sumx);
      free(sumy);
      free(sumzz);
      free(retdata);
      return rc;
    }
   
    for (j = 0; j < sum_dim; j++){
      switch(datatype){
      case FC_DT_DOUBLE:
         curry = ((double*)data)[j];
        break;
      case FC_DT_INT:
         curry = (double)(((int*)data)[j]);
        break;
      default: // know its a float cause i checked it above
	curry = (double)(((float*)data)[j]);
      }
    
      sumx[j]+=currx;
      sumy[j]+=curry;
      sumzz[3*j+0]+=(currx*currx);
      sumzz[3*j+1]+=(currx*curry);
      sumzz[3*j+2]+=(curry*curry);
    }
  }

  //  printf("\n");

  //now calc least squares val
  for (j = 0; j < sum_dim; j++){
    double a,b,r2,sea,seb,s;
    int id, xx,xy,yy;

    //formulas from mathworld on least squares fitting
    sumx[j]/=(double)numStep;
    sumy[j]/=(double)numStep;
    sumzz[3*j+0]-=(sumx[j]*sumx[j]*(double)numStep);
    sumzz[3*j+1]-=(sumx[j]*sumy[j]*(double)numStep);
    sumzz[3*j+2]-=(sumy[j]*sumy[j]*(double)numStep);

    //    printf("sxx %10g sxy %10g syy %10g",sumzz[3*j],sumzz[3*j+1],sumzz[3*j+2]);

    xx = (FC_DBL_EQUIV(sumzz[3*j],0.0) ? 1:0);
    xy = (FC_DBL_EQUIV(sumzz[3*j+1],0.0) ? 1:0);
    yy = (FC_DBL_EQUIV(sumzz[3*j+2],0.0) ? 1:0);

    // if xx means xx == 0
    if (xx) id = 1;
    if (!xx && xy && yy) id = 2;
    if (!xx && xy && !yy) id = 3;
    if (!xx && !xy && yy) id = 4;
    if (!xx && !xy && !yy) id = 5;

    switch(id){
    case 1:
      //just returning junk
      a = -1;
      b = -1;
      r2 = -1;
      sea = -1; 
      seb = -1; 
      break;
    case 2: 
      a = sumy[j];
      b = 0;
      r2 = -1;
      sea = 0;
      seb = 0;
      break;
    case 3:  //sea and seb calc below
      a = sumy[j];
      b = 0;
      r2 = 0;
      if (numStep > 2){ //s cant be zero
        s = sqrt(sumzz[3*j+2]/(double)(numStep-2)); 
        //      printf(" s %10g",s); 
      }
      break;
    case 4: //sea and seb calc below
      b = sumzz[3*j+1]/sumzz[3*j+0]; 
      a = sumy[j]- b*sumx[j];
      r2 = -1;
      if (numStep >2){  //s cant be zero
        s = sqrt((sumzz[3*j+1]*sumzz[3*j+1]/sumzz[3*j])/(double)(numStep-2));
      }
      break;
    case 5: //sea and seb calc below
      b =  sumzz[3*j+1]/sumzz[3*j+0]; 
      a = sumy[j]- b*sumx[j];
      r2 = sumzz[3*j+1]*sumzz[3*j+1]/(sumzz[3*j+0]*sumzz[3*j+2]); 
      if (numStep > 2){
        //check if s == 0
        double p1, p2;
        p1 = sumzz[3*j+2]/(double)(numStep-2);
        p2 = (sumzz[3*j+1]/(double)(numStep-2))*(sumzz[3*j+1]/sumzz[3*j]);
        //      printf(" inside %10g eps %10g",p1-p2, DBL_EPSILON);
        if (FC_VALUE_EQUIV(p1,p2,50*DBL_EPSILON,DBL_MIN)){
          s = 0;
        }else {
          s = sqrt(p1-p2);
        }
      }
      break;
    default:
      printf("** Developer Alert! Should never reach this point!**\n");
      fflush(NULL);
      fc_printfErrorMessage("Invalid SS values");
      return FC_ERROR;
      break;
    }
      
    if (id ==3 || id == 4 || id == 5){
      if (numStep < 3){
        sea = -1;
        seb = -1;
      } else {
        //      printf(" s %10g",s); 
        if (FC_DBL_EQUIV(s,0.0)){ 
          sea = 0;
          seb = 0;
        }else {
          double acoeff, bcoeff;        
          acoeff = sqrt(1.0/(double)numStep + sumx[j]*sumx[j]/sumzz[3*j]);
          bcoeff = sqrt(sumzz[3*j]);
          //      printf(" coeffa %10g coeffb %10g",acoeff,bcoeff); 
          sea = s*acoeff;
          seb = s/bcoeff;
        }
      }
    }

    //    printf ("\n");
    retdata[5*j] = a;
    retdata[5*j+1] = b;
    retdata[5*j+2] = r2;
    retdata[5*j+3] = sea;
    retdata[5*j+4] = seb;
  }

  free(sumx);
  free(sumy);
  free(sumzz);

  rc = fc_getMeshFromVariable(seqvar[0],&mesh);
  if (rc != FC_SUCCESS){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(rc));
    free(retdata);
    return rc;
  }


  varname = varSlot->header.name;
  returnname = malloc((strlen(varname)+20)*sizeof(char));
  if (returnname == NULL) {
    //no clean up in case of memory error
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  sprintf(returnname, "%s_LSQ", varname);
  
  rc = fc_createVariable(mesh,returnname,newvar);
  if (rc != FC_SUCCESS){
    free(retdata);
    free(returnname);
    return rc;
  }

  free(returnname);

  rc = fc_setVariableDataPtr(*newvar,numDataPoint,5,
                             assoc,FC_MT_VECTOR,FC_DT_DOUBLE,retdata);
  if (rc != FC_SUCCESS){
    fc_deleteVariable(*newvar);
    return rc;
  }

  return FC_SUCCESS;
}

//@}

/** \name Functions that operate on a pair of sequence vars and return a value or values based on those vars but the return val is not a sequence var */
//-------------------------------
//@{

/** 
 * \ingroup Series
 * \brief Return a handle to a var that is the euclidean distance
 *        between two curves traced out by sequence variables where
 *        y = seq var data value, x = sequence value (e.g., time).
 *
 * \description
 *
 *    Return a handle to a var that is the euclidean distance
 *    between two curves traced out by sequence variables where
 *        y = seq var data value, x = sequence value (e.g., time).
 *    NOTE: this is not normalized in any way (ie, it is dependent upon the
 *    number of data points and the spacing)
 *    Does a distance calc for each data point and each component separately.
 *    i.e. original seqvar has (numseriessteps*numdatapoint*numcomponent)
 *    dimensions; the result has numdatapoint*numcomponent) dimensions.
 *
 *    Always returns datatype FC_DT_DOUBLE.
 *
 *    This will return with an error if the seqvars are not numerical.
 *    The seq vars must be of the same mathtype, datatype, and association.
 *    For now seq vars must be on same mesh and same sequence.
 *
 *    This does *not* check that the sequence is a function.
 *    This does *not* check that the sequence is monotonically increasing.
 * 
 *
 * \todo
 *
 *    - In future may want to see about issues if the time pts arent regular.
 *    - Also may want to relax constraint of same mesh and same sequence -
 *    there are issues here however (comparing things disjoint in time etc).
 *    - Had changed this to each component separately - see if there is a
 *    reason to want to combine them.
 *    - should it check if the sequence is a fucntion ?
 *    - should it check if the sequence is monotonically increasing ?
 *
 * \modifications  
 *   - 1/14/05 ACG, Created
 *   - 4/25/05 ACG changed to use new functions fc_calcSeqVarUnaryFunction and
 *                 fc_calcSeqVarBinaryFunction
 *   - 4/26/05 ACG because using fc_getSeqVariableTimeSum return type
 *                 must always be double.
 *   - 6/07/05 ACG replaced calls to fc_calcSeqVarUnaryFunction and
 *                 fc_calcSeqVarBinaryFunction with calls to
 *                 appropriate seqvar functions in varmath
 *   - 7/06/05 ACG name change
 */
FC_ReturnCode fc_euclideanDistanceBetweenCurves(
     int numStep,  /**< input - number of steps in the seq var  */
     FC_Variable *seqvar1, /**< input - 1st seq var  */
     FC_Variable *seqvar2, /**< input - 2nd seq var  */
     char* varname, /**< input - name for averaged output var  */
     FC_Variable *newvar /**< output - handle for output var  */
){
  FC_ReturnCode rc,rc1,rc2;
  _FC_VarSlot *varSlot1, *varSlot2;
  int numDataPoint1, numDataPoint2, numComponent1, numComponent2;
  FC_AssociationType assoc1, assoc2;
  FC_MathType mathtype1, mathtype2;
  FC_DataType datatype1, datatype2;
  FC_Sequence seq1, seq2;
  FC_Mesh mesh1, mesh2;

  FC_Variable *diffvar, *prodvar;
  FC_Variable sumvar;


  if (newvar)
    *newvar = FC_NULL_VARIABLE;

  if (newvar == NULL || varname == NULL ||
      !fc_isSeqVariableValid(numStep,seqvar1) ||
      !fc_isSeqVariableValid(numStep,seqvar2)){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  //check seqvars are compatible
  rc1 = fc_getVariableInfo(seqvar1[0], &numDataPoint1, &numComponent1, &assoc1,
                          &mathtype1, &datatype1);
  rc2 = fc_getVariableInfo(seqvar2[0], &numDataPoint2, &numComponent2, &assoc2,
                           &mathtype2, &datatype2);
  if (rc1 != FC_SUCCESS || rc2 != FC_SUCCESS)
    return (rc1 != FC_SUCCESS ? rc1: rc2);

  if ((numDataPoint1 != numDataPoint2) ||
      (numComponent1 != numComponent2) ||
      (assoc1 != assoc2) ||
      (mathtype1 != mathtype2) ||
      (datatype1 != datatype2)){
    fc_printfErrorMessage("Mismatch in seqvar info");
    return FC_INPUT_ERROR;
  }

  //check datatypes are right
  if (!(datatype1 == FC_DT_FLOAT ||
        datatype1 == FC_DT_DOUBLE ||
        datatype1 == FC_DT_INT )
      ){
    fc_printfErrorMessage("seqvars are not numerical");
    return FC_INPUT_ERROR;
  } 

  //more checking seq vars compatible
  rc1 = fc_getSequenceFromSeqVariable(numStep,seqvar1,&seq1);
  rc2 = fc_getSequenceFromSeqVariable(numStep,seqvar2,&seq2);  
  if (rc1 != FC_SUCCESS || rc2 != FC_SUCCESS)
    return (rc1 != FC_SUCCESS ? rc1: rc2);

  if (!(FC_HANDLE_EQUIV(seq1,seq2))){
    fc_printfErrorMessage("seq vars must be on same seq");
    return FC_INPUT_ERROR;
  }


  rc1 = fc_getMeshFromVariable(seqvar1[0],&mesh1);
  rc2 = fc_getMeshFromVariable(seqvar2[0],&mesh2);
  if (rc1 != FC_SUCCESS || rc2 != FC_SUCCESS)
    return (rc1 != FC_SUCCESS ? rc1: rc2);

  if (!(FC_HANDLE_EQUIV(mesh1,mesh2))){
    fc_printfErrorMessage("seq vars must be on same mesh");
    return FC_INPUT_ERROR;
  }

  varSlot1 = _fc_getVarSlot(seqvar1[0]);
  varSlot2 = _fc_getVarSlot(seqvar2[0]);
  fc_printfLogMessage("Getting euclidean distance for seq variables '%s' '%s'",
                      varSlot1->header.name,varSlot2->header.name);

  varSlot1 = NULL;
  varSlot2 = NULL;


  //first get the difference between the seq vars
  rc = fc_seqVarOperatorSeqVar(numStep,seqvar1,"-",seqvar2,
                                   "fc_junk_diff_var", &diffvar);
  if (rc != FC_SUCCESS){
    fc_printfErrorMessage("cannot get diff between two sequnces");
    return rc;
  }

  //now square it 
  rc = fc_seqVarOperatorSeqVar(numStep,diffvar,"*",diffvar,
                                   "fc_junk_diff_var", &prodvar);
  if (rc != FC_SUCCESS){
    fc_printfErrorMessage("cannot get diff between two sequnces");
    fc_deleteSeqVariable(numStep,diffvar);
    free(diffvar);
    return rc;
  }
  rc =  fc_deleteSeqVariable(numStep,diffvar);
  free(diffvar);
    
  //now sum each data point and each component over time
  // returns doubles
  rc = fc_getSeqVariableSeriesSum (numStep, prodvar, &sumvar);
  if (rc != FC_SUCCESS){
    fc_printfErrorMessage("cannot get series sum");
    fc_deleteSeqVariable(numStep,prodvar);
    free(prodvar);
    return rc;
  }
  rc = fc_deleteSeqVariable(numStep,prodvar);
  free(prodvar);

  //now need sqrt of the series sum
  rc = fc_varUnaryFunction(sumvar, sqrt, varname, newvar);
  if (rc != FC_SUCCESS){
    fc_printfErrorMessage("cannot get sqrt of var");
    fc_deleteVariable(sumvar);
    return rc;
  }

  return FC_SUCCESS;
}


/** 
 * \ingroup Series
 * \brief Return a handle to a var that is a dimensionless area
 *        between two curves traced out by sequence variables where
 *        y = seq var data value, x = sequence value (e.g., time).
 *
 * \description
 *
 *    Return a handle to a var that is a dimensionless area
 *    between two curves traced out by sequence variables where
 *    y = seq var data value, x = sequence value (e.g., time).
 *    This is defined:
 *    (1/(xrange*ref_fct_max))*
 *    (integral sqr(reference variable - comparison variable)

 *    This function is the formula used by Patrick Klein
 *    Note that that is the absolute max on the ref fct, not 
 *    max absolute value, so its subject to the placement
 *    of the functions on the y axis. Note that the variables 
 *    are not symmetric.
 *
 *    Does a calc for each data point and each component separately.
 *    i.e. original seqvar has (numseriessteps*numdatapoint*numcomponent)
 *    dimensions; the result has numdatapoint*numcomponent) dimensions.
 *
 *    This will return with an error if the seq is not monotonically
 *    increasing, if the seqvars are not numerical, or if the sequence
 *    has less than 2 points. The seq vars must be of the same mathtype,
 *    datatype, and association. For now seq vars must be on same mesh
 *    and same sequence. 
 *
 *    This will always return datatype FC_DT_DOUBLE
 *
 *
 * \todo
 *   - despite PK definition, should we go with max abs val on
 *     the ref fct, as opposed to just max val ? or some distance
 *     val that means things are not subject to where the functions
 *     are placed on the axis?
 * 
 * \modifications  
 *   - 4/26/05 ACG, Created
 *   - 6/07/05 ACG revised to use new seqvar varmath functions
 *   - 7/06/05 ACG name change
 *   - 1/20/06 ACG error checking revisions
 */
FC_ReturnCode fc_dimensionlessAreaBetweenCurves(
     int numStep,  /**< input - number of steps in the seq var  */
     FC_Variable *refseqvar, /**< input - reference seq var  */
     FC_Variable *compseqvar, /**< input - comparison seq var  */
     char* varname, /**< input - name for averaged output var  */
     FC_Variable *newvar /**< output - handle for output var  */
){
  FC_ReturnCode rc,rc1,rc2;
  _FC_VarSlot *varSlot1, *varSlot2;
  int numDataPoint1, numDataPoint2, numComponent1, numComponent2;
  FC_AssociationType assoc1, assoc2;
  FC_MathType mathtype1, mathtype2;
  FC_DataType datatype1, datatype2;
  FC_Sequence seq1, seq2;
  FC_Mesh mesh1, mesh2;

  FC_Mesh mesh;
  FC_Variable *diffvar, *prodvar, intvar;
  FC_Variable minvar,maxvar,minindexvar,maxindexvar;
  FC_DataType datatype;
  FC_AssociationType assoc;
  FC_MathType mathtype;
  int numDataPoint, numComponent;
  double *maxdata, *intdata, *writedata;
  double min,max;
  int minindex, maxindex,mono;
  double range;
  int i,j;

  if (newvar)
    *newvar = FC_NULL_VARIABLE;

  if (newvar == NULL || varname == NULL){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  if ((!fc_isSeqVariableValid(numStep,refseqvar)) || 
      (!fc_isSeqVariableValid(numStep,compseqvar))){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  //check seqvars are compatible
  rc1 = fc_getVariableInfo(refseqvar[0], &numDataPoint1, &numComponent1, &assoc1,
                          &mathtype1, &datatype1);
  rc2 = fc_getVariableInfo(compseqvar[0], &numDataPoint2, &numComponent2, &assoc2,
                           &mathtype2, &datatype2);
  if (rc1 != FC_SUCCESS || rc2 != FC_SUCCESS) 
    return (rc1 != FC_SUCCESS ? rc1: rc2);


  if ((numDataPoint1 != numDataPoint2) ||
      (numComponent1 != numComponent2) ||
      (assoc1 != assoc2) ||
      (mathtype1 != mathtype2) ||
      (datatype1 != datatype2)){
    fc_printfErrorMessage("seqvars not compatible");
    return FC_INPUT_ERROR;
  }

  //check datatypes are right
  if (!(datatype1 == FC_DT_FLOAT || datatype1 == FC_DT_DOUBLE ||
	datatype1 == FC_DT_INT )
      ){
    fc_printfErrorMessage("seqvar data type not numerical");
    return FC_INPUT_ERROR;
  } 

  //more checking seq vars compatible
  rc1 = fc_getSequenceFromSeqVariable(numStep,refseqvar,&seq1);
  rc2 = fc_getSequenceFromSeqVariable(numStep,compseqvar,&seq2);
  if (rc1 != FC_SUCCESS || rc2 != FC_SUCCESS) 
    return (rc1 != FC_SUCCESS ? rc1: rc2);


  if (!(FC_HANDLE_EQUIV(seq1,seq2))){
    fc_printfErrorMessage("seqvars must be on same sequence");
    return FC_INPUT_ERROR;
  }

  if (numStep < 2){
    fc_printfErrorMessage("seqvars must have at least 2 datapoints");
    return FC_INPUT_ERROR;
  }

  rc=fc_getSequenceMinMaxMono(seq1,&min,&minindex,&max,&maxindex,&mono);
  if (rc !=FC_SUCCESS){
    return rc;
  }

  if (mono != 1){
    fc_printfErrorMessage("sequence must be monotonically increasing");
    return FC_ERROR;
  }

  range = max-min;
  if (FC_DBL_EQUIV(range,0.0)){ //check range here, check max & max*range below
    fc_printfErrorMessage("cant get area between zero range sequences");
    return FC_INPUT_ERROR;
  }


  rc1 = fc_getMeshFromVariable(refseqvar[0],&mesh1);
  rc2 = fc_getMeshFromVariable(compseqvar[0],&mesh2);
  if (rc1 != FC_SUCCESS || rc2 != FC_SUCCESS) 
    return (rc1 != FC_SUCCESS ? rc1: rc2);

  if (!(FC_HANDLE_EQUIV(mesh1,mesh2))){
    return FC_INPUT_ERROR;
  }


  varSlot1 = _fc_getVarSlot(refseqvar[0]);
  varSlot2 = _fc_getVarSlot(compseqvar[0]);
  fc_printfLogMessage("Getting area between curves of seq variables '%s' '%s'",
                      varSlot1->header.name,varSlot2->header.name);

  varSlot1 = NULL;
  varSlot2 = NULL;

  rc = fc_getVariableInfo(refseqvar[0],&numDataPoint, &numComponent,
                          &assoc,&mathtype, &datatype);
  if (rc != FC_SUCCESS){
    fc_printfErrorMessage("cant get info on ref variable");
    return rc;
  }

  
  // now get max of ref fcnt and make sure it is non zero....
  rc =  fc_getSeqVariableSeriesMinMax(numStep, refseqvar,&minvar,
				      &maxvar, &minindexvar,&maxindexvar);
  if (rc != FC_SUCCESS){
    fc_printfErrorMessage("cant get min max of ref variable");
    return rc;
  }
  
  rc=fc_deleteVariable(minvar);
  rc=fc_deleteVariable(minindexvar);
  rc=fc_deleteVariable(maxindexvar);

  rc = fc_getVariableDataPtr(maxvar,(void**)&maxdata);
  if (rc != FC_SUCCESS){
    fc_printfErrorMessage("cant get max data");
    fc_deleteVariable(maxvar);
    return rc;
  }


  //its a shame to loop thru these twice but it will save
  //work below if there is a problem discoverable now.
  for (i = 0; i < numDataPoint; i++){
    for (j = 0; j < numComponent; j++){
      if (FC_DBL_EQUIV(maxdata[numComponent*i+j],0.0) || 
	  FC_DBL_EQUIV(maxdata[numComponent*i+j]*range,0.0)){
        fc_printfErrorMessage("cant normalize by zero max function or max*range");
        fc_deleteVariable(maxvar);
        return FC_ERROR;
      }
    }
  }

  //now get the difference between the seq vars
  rc = fc_seqVarOperatorSeqVar(numStep,refseqvar,"-",compseqvar,
                                   "fc_junk_diff_var", &diffvar);
  if (rc != FC_SUCCESS){
    fc_printfErrorMessage("cannot get diff between two sequnces");
    fc_deleteVariable(maxvar);
    return rc;
  }

  //now square it
  rc = fc_seqVarOperatorSeqVar(numStep,diffvar,"*",diffvar,
                                   "fc_junk_prod_var", &prodvar);
  if (rc != FC_SUCCESS){
    fc_printfErrorMessage("cannot get square");
    fc_deleteVariable(maxvar);
    fc_deleteSeqVariable(numStep,diffvar);
    free(diffvar);
    return rc;
  }
  rc =  fc_deleteSeqVariable(numStep,diffvar);
  free(diffvar);

  //now integrate it - returns doubles. no longer a seq var
  rc = fc_integral_TR(numStep,prodvar,"fc_junk_int_var",
                                &intvar);
  if (rc != FC_SUCCESS){
    fc_printfErrorMessage("cannot get integral");
    fc_deleteVariable(maxvar);
    fc_deleteSeqVariable(numStep,prodvar);
    free(prodvar);
    return rc;
  }
  rc =  fc_deleteSeqVariable(numStep,prodvar);
  free(prodvar);


  //now do calc
  //still have maxdataptr

  rc = fc_getVariableDataPtr(intvar,(void**)&intdata);
  if (rc != FC_SUCCESS){
    fc_printfErrorMessage("cant get integral data");
    fc_deleteVariable(maxvar);
    fc_deleteVariable(intvar);
    return rc;
  }

  writedata = (double*) malloc(numDataPoint*numComponent*sizeof(double));

  for (i = 0; i < numDataPoint; i++){
    for (j = 0; j < numComponent; j++){
      //have already checked for divide by zero
      //      printf("int %10g max %10g range %10g\n",intdata[numComponent*i+j],
      //	     maxdata[numComponent*i+j],range);
      writedata[numComponent*i+j] = intdata[numComponent*i+j]/
				      (maxdata[numComponent*i+j]*range);
    }
  }

  fc_deleteVariable(maxvar);
  fc_deleteVariable(intvar);

  //assign data to variable
  rc = fc_getMeshFromVariable(refseqvar[0],&mesh);
  if (rc != FC_SUCCESS){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(rc));
    free (writedata);
    return rc;
  }

  rc = fc_createVariable(mesh,varname,newvar);
  if (rc != FC_SUCCESS){
    free(writedata);
    return rc;
  }

  rc = fc_setVariableDataPtr(*newvar,numDataPoint, numComponent,
                             assoc,mathtype,FC_DT_DOUBLE,writedata);
  if (rc != FC_SUCCESS){
    if (writedata) free (writedata);
    fc_deleteVariable(*newvar);
    return rc;
  }

  return FC_SUCCESS;
}


//@}


