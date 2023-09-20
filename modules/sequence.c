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
 * \file sequence.c
 * \brief Implementation for \ref Sequence module.
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/modules/sequence.c,v $
 * $Revision: 1.54 $ 
 * $Date: 2006/09/22 21:38:57 $
 */

// C library includes
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// fc library dependencies
#include "base.h"
#include "storage.h"
#include "library.h"
#include "util.h"
#include "tableP.h"
#include "datasetP.h"
#include "meshP.h"
#include "variableP.h"
#include "fileioP.h"
#include "statistics.h"

// this module
#include "sequence.h"
#include "sequenceP.h"

/**
 * \addtogroup  Sequence
 * \brief  Operations on sequences.
 *
 * \description
 *
 *    (For an explanation of general data manipulations, see 
 *    \ref DataInterface.)
 */

/**
 * \ingroup Sequence
 * \defgroup PrivateSequence (Private)
 */

/** 
 * \ingroup PrivateSequence
 * \brief The number of allocated slots in the sequence table.
 */
static int seqTableSize = 0;
/** 
 * \ingroup PrivateSequence
 * \brief The sequence table.
 */
static _FC_SeqSlot **seqTable = NULL;
/**
 * \ingroup PrivateSequence
 * \brief The list of unused slots in the sequence table.
 */
static FC_SortedIntArray seqOpenSlots = { 0, 0, 0 };

/** \name Create new sequence from scratch. */
//-------------------------------
//@{

/**
 * \ingroup  Sequence
 * \brief  Create a new sequence.
 *
 * \description
 *
 *    Create a new, in-memory, sequence.
 *
 * \modifications  
 *    - Nancy Collins, Created.
 */
FC_ReturnCode fc_createSequence(
  FC_Dataset dataset,    /**< input - dataset to create sequence in */
  char *seqname,         /**< input - name of sequence */
  FC_Sequence *sequence  /**< output - new sequence handle */
) {
  _FC_DsSlot* dsSlot;
  _FC_SeqSlot* seqSlot;
  int *temp_slotIDs;
  
  // default return value
  if (sequence)
    *sequence = FC_NULL_SEQUENCE;
  
  // check input
  dsSlot = _fc_getDsSlot(dataset);
  if (dsSlot == NULL || !seqname || !sequence) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Creating new sequence '%s'", seqname);

  // get an open slot
  seqSlot = _fc_getNewSeqSlot();
  if (seqSlot == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  
  // table header information
  _FC_GET_HANDLE(*sequence, seqSlot);
  _fc_setSlotHeaderName(&(seqSlot->header), seqname);
  
  // set back and forth references
  seqSlot->ds = dataset;
  
  // update owning dataset with info about this sequence
  temp_slotIDs = realloc(dsSlot->seqIDs, sizeof(int) * (dsSlot->numSeq + 1));
  if (temp_slotIDs == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  temp_slotIDs[dsSlot->numSeq] = seqSlot->header.slotID;
  dsSlot->seqIDs = temp_slotIDs;
  dsSlot->numSeq++;
  
  return FC_SUCCESS;
}

/**
 * \ingroup  Sequence
 * \brief  Set the sequence coordinates by copy.
 *
 * \description
 *  
 *    Set the coordinates for the sequence. Must provide the number of steps,
 *    the data type, and a pointer to the 1-D array of coordinates.
 *    This routine makes a copy of the coords for the library,
 *    so the user still has responsibility of the original buffer.
 *    If you would like to give the coordinates array to the sequence instead
 *    of just copying them, use fc_setSequenceCoordsPtr().
 *
 *    You cannot set coords if they already exist.
 *
 * \modifications 
 *    - 2003-NOV-13  WSK  Created
 *    - 09/10/04 WSK, changed so cannot set coords more than once.
 */
FC_ReturnCode fc_setSequenceCoords(
  FC_Sequence sequence,  /**< Input - Sequence handle */
  int numStep,           /**< Input - Number of steps in the sequence */
  FC_DataType dataType,  /**< Input - The data type of the coordinate array
                                      (i.e. int, float, etc) */
  void *coords           /**< Input - Array of coordinate data */
) {
  int size;
  _FC_SeqSlot *seqSlot;
  
  // NOTE!: if you change stuff here, also change fc_setSetSequenceCoordsPtr()

  // check input
  seqSlot = _fc_getSeqSlot(sequence);
  if (seqSlot == NULL || numStep < 1 || !fc_isDataTypeValid(dataType) ||
      dataType == FC_DT_UNKNOWN || coords == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
 
  // a little more checking -- can't overwrite coords 
  if (seqSlot->coords != NULL || seqSlot->header.committed == 1) {
    fc_printfErrorMessage("Coords already exist on sequence '%s'",
                          seqSlot->header.name);
    return FC_ERROR;
  }

  // log message
  fc_printfLogMessage("Setting coords on sequence '%s'", seqSlot->header.name);

  // Meta data
  seqSlot->numStep = numStep;
  seqSlot->datatype = dataType;

  // Big data
  size = fc_sizeofDataType(dataType);
  size *= numStep;
  seqSlot->coords = malloc(size);
  if (seqSlot->coords == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      return FC_MEMORY_ERROR;
  }
  memcpy(seqSlot->coords, coords, size);

  return FC_SUCCESS;
}

/**
 * \ingroup  Sequence
 * \brief  Set the pointer to the mesh coordinates.
 *
 * \description
 *  
 *    Set the coordinates for the sequence. Must provide the number of steps,
 *    the data type, and a pointer to the 1-D array of coordinates. 
 *    This routine passes the coords array directly to the library,
 *    so the user no longer has responsibility of the original buffer.
 *    If you would like to keep the coords array, you should use 
 *    fc_setSequenceCoords() instead which makes a copy of the coords array.
 *
 *    You cannot set coords if they already exist.
 *
 * \modifications 
 *    - 10/01/04 WSK, created.
 */
FC_ReturnCode fc_setSequenceCoordsPtr(
  FC_Sequence sequence,  /**< Input - Sequence handle */
  int numStep,           /**< Input - Number of steps in the sequence */
  FC_DataType dataType,  /**< Input - The data type of the coordinate array
                                      (i.e. int, float, etc) */
  void *coords           /**< Input - Array of coordinate data */
) {
  _FC_SeqSlot *seqSlot;
  
  // NOTE!: if you change stuff here, also change fc_setSetSequenceCoords()

  // check input
  seqSlot = _fc_getSeqSlot(sequence);
  if (seqSlot == NULL || numStep < 1 || !fc_isDataTypeValid(dataType) ||
      dataType == FC_DT_UNKNOWN || coords == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
 
  // a little more checking -- can't overwrite coords 
  if (seqSlot->coords != NULL || seqSlot->header.committed == 1) {
    fc_printfErrorMessage("Coords already exist on sequence '%s'",
                          seqSlot->header.name);
    return FC_ERROR;
  }

  // log message
  fc_printfLogMessage("Setting coords on sequence '%s'", seqSlot->header.name);

  // Meta data
  seqSlot->numStep = numStep;
  seqSlot->datatype = dataType;

  // Big data
  seqSlot->coords = coords;

  return FC_SUCCESS;
}

//@}

/** \name Other ways to get new sequences. */
//-------------------------------
//@{

/**
 * \ingroup  Sequence
 * \brief  Copy a sequence
 *
 * \description
 *  
 *    Create a copy of a sequence, in the specified dataset, with the 
 *    given name.
 *
 * \todo
 *    Perhaps somehow flag when the destination dataset is the same
 *    as the source dataset. Also do this for other copies. Very useful
 *    for fc_copyVariable.
 *
 * \modifications 
 *   - 2003-NOV-11 WSK Created
 *   - 2008-MAR-26 CDU Updated to allow zero-length sequences
 */
FC_ReturnCode fc_copySequence(
  FC_Sequence src_seq, /**< Input - the source sequence to be copied */
  FC_Dataset dest_ds, /**< Input - the destination dataset to be copied into */
  char *newName,       /**< Input - the name of the new sequence, a NULL value
                            is a flag to use the name of the source sequence */
  FC_Sequence* new_seq /**< Output - handle to the new sequence */
) { 
  FC_ReturnCode rc;
  _FC_DsSlot* dest_dsSlot;
  _FC_SeqSlot* src_seqSlot;
  int numStep;
  FC_DataType dataType;
  void* coords;

  // default return
  if (new_seq)
    *new_seq = FC_NULL_SEQUENCE;

  // check input
  dest_dsSlot = _fc_getDsSlot(dest_ds);
  src_seqSlot = _fc_getSeqSlot(src_seq);
  if (dest_dsSlot == NULL || src_seqSlot == NULL || new_seq == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // test for NULL name
  if (newName == NULL)
    newName = src_seqSlot->header.name;

  // log message
  fc_printfLogMessage("Copying sequence '%s' to '%s'", 
                      src_seqSlot->header.name, newName);

  // get src data
  rc = fc_getSequenceInfo(src_seq, &numStep, &dataType);
  if (rc != FC_SUCCESS) {
    fc_printfErrorMessage("Failed to get info for sequence '%s'",
                          src_seqSlot->header.name);
    return rc;
  }
  rc = fc_getSequenceCoordsPtr(src_seq, &coords);
  if (rc != FC_SUCCESS) {
    fc_printfErrorMessage("Failed to get coords for sequence '%s'",
                          src_seqSlot->header.name);
    return rc;
  }

  // create new sequence 
  rc = fc_createSequence(dest_ds, newName, new_seq);
  if (rc != FC_SUCCESS) {
    fc_printfErrorMessage("Failed to create sequence '%s'", newName);
    return rc;
  }

  // Copy the coordinates
  if(numStep==0){
    //Some datasets have sequences that have 0 steps (eg, Cubit
    //will add a sequence named "time" that has 0 steps). Thus, we
    //allow users to copy the empty subset from one dataset to another.
    fc_printfWarningMessage("Copying an empty sequence named '%s'",
			    src_seqSlot->header.name );

  } else {
    rc = fc_setSequenceCoords(*new_seq, numStep, dataType, coords);
    if (rc != FC_SUCCESS) {
      fc_printfErrorMessage("Failed to set coords for sequence '%s'", newName);
      return rc;
    }
  }

  return rc;
}


/** 
 * \ingroup Sequence
 * \brief Given a seq, return a new sequence that is the old one shifted and scaled.
 *
 * \description
 *
 *    Given a seq, return a new sequence that is the old one shifted and
 *    scaled.  Possible usage: This function can be followed by
 *    fc_copySeqVariable in order to create a new seqVar that is an oldSeqVar
 *    normalized in time such that it starts at time zero.
 *    
 *    The shift is applied before the scale. 
 *    i.e., new_seq_val = (old_seq_val+shift)*(scale).
 *    If you don't want to shift, pass in 0 for the shift; if you
 *    dont want to scale pass in 1 for the scale.
 *
 *
 *    The datatype of the new sequence is FC_DT_DOUBLE regardless of
 *    the original sequence datatype.
 *
 *    This returns with an error if the orignal sequence is not numerical.
 *    Scaling by zero is not allowed.
 *
 * \todo change this so the user can specify the return type.
 *
 * \modifications
 *    6/14/05 ACG created. 
 *
 */
FC_ReturnCode fc_shiftAndScaleSequence(
  FC_Sequence src_seq, /**< Input - the source sequence to be shifted and scaled */
  FC_Dataset dest_ds, /**< Input - the destination dataset to be copied into */
  double shift, /**< input - val shifting seq by */ 
  double scale, /**< input - val scaling seq by */ 
  char *newName,       /**< Input - the name of the new sequence */
  FC_Sequence* new_seq /**< Output - handle to the new sequence */
) { 
  FC_ReturnCode rc;

  _FC_DsSlot* dest_dsSlot;
  _FC_SeqSlot* src_seqSlot;
  int numStep;
  FC_DataType seqdatatype;
  double* prev;
  void* coordsp;
  int i;

  // default return
  if (new_seq)
    *new_seq = FC_NULL_SEQUENCE;

  // check input
  dest_dsSlot = _fc_getDsSlot(dest_ds);
  src_seqSlot = _fc_getSeqSlot(src_seq);
  if (dest_dsSlot == NULL || src_seqSlot == NULL ||
      new_seq == NULL || newName == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // log message
  fc_printfLogMessage("shifting and scaling sequence '%s' to '%s'", 
                      src_seqSlot->header.name, newName);

  if (FC_DBL_EQUIV(scale,0.0)){
    fc_printfErrorMessage("scaling by zero is not allowed");
    return FC_INPUT_ERROR;
  }

  rc = fc_getSequenceInfo (src_seq, &numStep, &seqdatatype);
  if (rc != FC_SUCCESS){
    return rc;
  }

  if (!(seqdatatype == FC_DT_INT || seqdatatype == FC_DT_FLOAT ||
        seqdatatype == FC_DT_DOUBLE)){
    fc_printfErrorMessage("seq data type is not numerical");
    return FC_ERROR;
  }

  rc = fc_getSequenceCoordsPtr(src_seq,&coordsp);
  if (rc!= FC_SUCCESS){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_ERROR;
  }


  prev = (double*)malloc(numStep*sizeof(double));
  if (prev == NULL){
    fc_printfErrorMessage("%s","cant make data array");
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }

  for (i = 0; i < numStep; i++){
    switch(seqdatatype){
    case FC_DT_INT:
      prev[i] = (double)(((int*)coordsp)[i]);
      break;
    case FC_DT_FLOAT:
      prev[i] = (double)(((float*)coordsp)[i]);
      break;
    default: //know its a double
      prev[i] = ((double*)coordsp)[i];
      break;
    }

    //this isnt quite ideal - e.g., if magnitude of both shift and scale are small...
    if (FC_DBL_EQUIV(prev[i],-shift)){
      prev[i] = 0.0;
    }else{
      prev[i]+=shift;
      if (!FC_DBL_EQUIV(scale,1.0)){
        prev[i]*=scale;
      }
    }
  }

  // create new sequence 
  rc = fc_createSequence(dest_ds, newName, new_seq);
  if (rc != FC_SUCCESS) {
    fc_printfErrorMessage("Failed to create sequence '%s'", newName);
    free(prev);
    return rc;
  }

  // set new seq coords
  rc = fc_setSequenceCoordsPtr(*new_seq,numStep,FC_DT_DOUBLE,(void*)prev);
  if (rc!= FC_SUCCESS){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    free(prev);
    fc_deleteSequence(*new_seq);
    *new_seq = FC_NULL_SEQUENCE;
    return FC_ERROR;
  }
  
  return rc;
}


/** 
 * \ingroup Sequence
 * \brief Creates a sequence of specified length, starting at a given
 *        val with subsequent vals spaced by given stepsize
 *
 * \description
 *
 *        Creates a sequence of specified with starting at a given
 *        val with subsequent vals spaced by given stepsize. Stepsize
 *        cannot = 0 unless it requesting a single stepsize sequence.
 *        If want a single stepsize sequence the spacing requested
 *        doesnt matter.
 *
 *        Return seq data type is FC_DT_DOUBLE.
 *
 *        This routine is often used as a companion function to 
 *        functions in the timeseries module. This function can
 *        be used to create a common sequence for two sequence
 *        variables. The old variables can be put on the new sequence
 *        using fc_LinearInterpolation and then used in the comparison
 *        routines in the timeseries module.
 *
 *        This is also is a companion function to the private function
 *        _fc_createSeqVariableFromRegularVariable in the timeseries module
 *        used in some test modules.In this way 2 sequence variables can be 
 *        created from 2 different multi-datapoint variables and used in
 *        the comparison routines in the timeseries module, since they will
 *        be on the same sequence. 
 *
 * \todo should there be another version where you set the bounds rather than
 *       the stepsize to avoid roundoff issues ?
 *
 * \modifications  
 *    - 2/2/05 ACG, created from being split out from 
 *      fc_createSeqvariableFromVariable
 *    - 2/10/05 ACG changed so takes an offset and a step size.
 *    - 2/10/05 ACG moved to private becuase thats where the companion
 *              functions _fc_createSeqVariableFromRegularVariable is.
 *              This may be useful and perhaps put in Sequence module, on
 *              the other hand its sufficiently short that maybe it cannot
 *              justify its existence.
 *    - 4/12/05 ACG and now we need it for comparing variables on
 *              different sequences.
 *    - 4/21/05 ACG moved it here from its previous home in timeseries.
 *    - 4/27/05 ACG better handling of single data point sequences
 *
 */
FC_ReturnCode fc_createRegularSequence(
  FC_Dataset dset, /**< input - dataset for the sequence */
  int numStep, /**< input - numstep in created seq */
  double offset, /**< initial val in the sequence */ 
  double stepsize, /**< spacing of the sequence */ 
  char* seq_name, /**< input - name of the new sequence that will be created */
  FC_Sequence *seq /**< output - handle to the new sequence */
){

  FC_ReturnCode rc;
  _FC_DsSlot* dsSlot;
  double *times;
  int k;


  // default return
  if (seq)
    *seq = FC_NULL_SEQUENCE;

  dsSlot = _fc_getDsSlot(dset);
  if (dsSlot == NULL || !seq_name || !seq || numStep < 1){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  if (FC_DBL_EQUIV(stepsize,0.0) && numStep != 1){
    fc_printfErrorMessage("cant ask for more than 1 step but zero step size");
    return FC_INPUT_ERROR;
  }

  fc_printfLogMessage("Creating regular sequence");

  rc = fc_createSequence(dset,seq_name,seq);
  if (rc != FC_SUCCESS){
    fc_deleteSequence(*seq);
    *seq = FC_NULL_SEQUENCE;
    return rc;
  }


  times = (double*)malloc(numStep*sizeof(double));
  times[0] = offset;
  for (k = 1; k < numStep; k++){
    times[k] = times[k-1]+stepsize;
  }
      
  rc = fc_setSequenceCoords(*seq,numStep,FC_DT_DOUBLE,times);
  free(times);
  if (rc != FC_SUCCESS){
    fc_deleteSequence(*seq);
    *seq = FC_NULL_SEQUENCE;
    return rc;
  }

  return FC_SUCCESS;

}


/** 
 * \ingroup Sequence
 * \brief Given an array of sequences on the same dataset, creates a new 
 *        regularly spaced sequence whose range is the intersection
 *        of the sequences.
 *
 * \description
 *
 *        Given an array of sequences on the same dataset, creates a new 
 *        regularly spaced sequence whose range is the intersection
 *        of the sequences. Input sequences must be numerical, but dont
 *        have to be same datatype. Input sequences must be 
 *        monotonically increasing. New seq data type is FC_DT_DOUBLE.
 *        Seq is null if no valid range. If numStep == 0, then it
 *        returns a sequence with the greatest number of steps that
 *        an input sequence has on that range. 
 *
 *        This routine is often used as a companion function to 
 *        functions in the timeseries module. This function creates
 *        a common sequence for two or more sequence variables. The
 *        old variables can be put on the new sequence
 *        using fc_LinearInterpolation and then used in the comparison
 *        routines in the timeseries module.
 *
 * \todo
 *     see if can combine this with rc_getIntersectingRangeRegularSequence 
 *
 * \modifications  
 *    - 4/25/05 ACG created.
 *    - 4/27/05 ACG better handling of single data point sequences
 *    - 4/28/05 ACG explicitly set last point of seq to timemax to
 *              avoid roundoff issues as a result of createRegularSequence
 *              being based on stepsize rather than range. if we make a
 *              another version of createRegularSequence that is based on
 *              ranage, then use that version here.
 *    - 5/13/05 ACG works for non-monotonic sequences so that
 *              can use this to check for any overlap to begin with. of course
 *              if you are going to map a seq var onto a new seq that was
 *              calculated by this method, then the seqvar will have to
 *              be made into a monotonic seq to before the interpolation.
 *              This change utilizes the new stats function 
 *              fc_getSequenceMinMaxMonotonicity.
 *    - 5/16/05 ACG making this so has to be monotonically increasing.
 *              this makes things easier for doing the automatic
 *              numsteps option.
 *    - 5/17/05 ACG added default case of numStep == 0 returns max
 *              number of steps that the input sequences have on that
 *              range. in this way you dont lose any resolution that an
 *              input sequence has (well, you kind of do if it wasnt regularly
 *              spaced to begin with) and you wont gain any when you do
 *              a linear interpolation over the new sequence becuase its just
 *              more points over what was a linear segment in the original
 *              resolution anyway. 
 *    - 6/27/05 ACG updated to use new getClosestSequenceValue
 *
 */
FC_ReturnCode fc_createIntersectingRegularSequence(
  int numseqs,/**< input - num of seqs in orig sequence array */
  FC_Sequence *seqs,/**< input - orig sequence array */
  int numStep, /**< input - numstep in created seq */
  char* seqName, /**<input - name of created sequence */
  FC_Sequence *seq /**< output - handle to the new sequence */
  ){

  FC_ReturnCode rc;
  FC_Dataset dset, dsetcomp;
  _FC_DsSlot* dsSlot;
  double timemin, timemincomp, timemax, timemaxcomp, timespacing;
  void *data;
  int nstep;
  int mono,monocomp;
  FC_DataType datatype;
  int i;


  if (seq)
    *seq = FC_NULL_SEQUENCE;

  if (!seqName || !seq || numStep < 0 || numseqs < 1 || !seqs){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  fc_printfLogMessage("Creating Intersecting regular sequence");

  if (!fc_isSequenceValid(seqs[0])){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  rc = fc_getDatasetFromSequence (seqs[0], &dset);
  if (rc != FC_SUCCESS){
    fc_printfErrorMessage("Cannot get dataset from sequence");
    return FC_ERROR;
  }

  dsSlot = _fc_getDsSlot(dset);
  if (dsSlot == NULL){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  rc = fc_getSequenceInfo(seqs[0], &nstep,&datatype);
  if (rc != FC_SUCCESS){
    fc_printfErrorMessage("failed to get seq info");
    return FC_ERROR;
  }

  if (datatype != FC_DT_INT && datatype != FC_DT_FLOAT &&
      datatype != FC_DT_DOUBLE){
    fc_printfErrorMessage("Invalid sequence datatype");
    return FC_INPUT_ERROR;
  }

  for (i= 1; i < numseqs; i++){
    if (!fc_isSequenceValid(seqs[i])){
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
      return FC_INPUT_ERROR;
    }
    rc = fc_getDatasetFromSequence (seqs[i], &dsetcomp);
    if (rc != FC_SUCCESS){
      fc_printfErrorMessage("Cannot get dataset from sequence");
      return FC_ERROR;
    }
    
    if (!FC_HANDLE_EQUIV(dset,dsetcomp)){
      fc_printfErrorMessage("Input sequences not on same dataset");
      return FC_INPUT_ERROR;
    } 

    rc = fc_getSequenceInfo(seqs[i], &nstep,&datatype);
    if (rc != FC_SUCCESS){
      fc_printfErrorMessage("failed to get seq info");
      return FC_ERROR;
    }

    if (datatype != FC_DT_INT && datatype != FC_DT_FLOAT &&
        datatype != FC_DT_DOUBLE){
      fc_printfErrorMessage("Invalid sequence datatype");
      return FC_INPUT_ERROR;
    }
  }


  rc = fc_getSequenceMinMaxMono(seqs[0],&timemin,NULL,
                                &timemax,NULL,&mono);
  if (rc != FC_SUCCESS){
    fc_printfErrorMessage("failed to get minmax of sequence");
    return FC_ERROR;
  }

  if (mono != 1 && mono != 2){
    fc_printfErrorMessage("invalid monotonicity for this sequence");
    return FC_INPUT_ERROR;
  }


  for (i = 1; i < numseqs; i++){
    rc = fc_getSequenceMinMaxMono(seqs[i],&timemincomp,NULL,
                                  &timemaxcomp,NULL,&monocomp);
    if (rc != FC_SUCCESS){
      fc_printfErrorMessage("failed to get minmax of sequence");
      return FC_ERROR;
    }


    //leaving this construction in case switch to allowing decreasing later
    if (!(monocomp == 1 || monocomp ==2)){ //increasing or single point
      fc_printfErrorMessage("invalid monotonicity for this call");
      return FC_INPUT_ERROR;
    } else {
      if (mono == 2){
        //its gotta be single point sequence
        mono = monocomp;
      }
    }

    //want narrowest bounds
    if (timemincomp > timemin) timemin= timemincomp;
    if (timemaxcomp < timemax) timemax= timemaxcomp;

    if (timemin > timemax){
      //no interesecting bounds - bad input set, but dont want to error
      fc_printfErrorMessage("no intersecting range");
      return FC_ERROR;
    }
  }

  //now have valid range
  
  //if numStep == 0 need to calc numStep
  if (numStep == 0){
    if (FC_DBL_EQUIV(timemax,timemin)){
      numStep = 1; //just make single point sequence
    } else {
      double closemintime, closemaxtime;
      int closeminindex, closemaxindex;
      int npts;
      int holdnpts = -1;

      for (i = 0; i < numseqs; i++){
        rc = fc_getClosestSequenceValue(seqs[i],"<=",timemin,&closemintime,
                                        &closeminindex);
        if (rc != FC_SUCCESS){
          fc_printfErrorMessage("failed to get closest seq val");
          return FC_ERROR;
        }
        rc = fc_getClosestSequenceValue(seqs[i],">=",timemax,&closemaxtime,
                                        &closemaxindex);
        if (rc != FC_SUCCESS){
          fc_printfErrorMessage("failed to get closest seq val");
          return FC_ERROR;
        }
        
        //know its monotonically increasing
        npts = closemaxindex-closeminindex+1;
        if (holdnpts == -1 || npts > holdnpts){
          holdnpts = npts;
        }
      }
      numStep = holdnpts;
    }
  }


  //check for zero range sequences
  if (numStep == 1){
    if (FC_DBL_EQUIV(timemax,timemin)){
      //good - createRegularsequence should just return the single
      //data point
      timespacing = 0.0; //actually this value doesnt matter
    } else {
      fc_printfErrorMessage("cant have nonzero range with only one data point");
      return FC_ERROR;
    }
  } else {
    if (FC_DBL_EQUIV(timemin,timemax)){
      fc_printfErrorMessage("cant have zero range with more than one data point");
      return FC_ERROR;
    } else {
      //good
      timespacing= (timemax-timemin)/((double)numStep-1);
    }
  }
  //printf("new seq: tmin %10g tmax %10g spacing %10g\n",timemin, timemax,
  //timespacing);
  rc = fc_createRegularSequence(dset,numStep,timemin, timespacing,
                                  seqName, seq);
  if (rc != FC_SUCCESS){
    fc_printfErrorMessage("failed to create reg sequence");
    return FC_ERROR;
  }
  
  //just in case there are some roundoff issues, lets explicitly set that
  //last point to the timemax. (this wont be thread safe since playing 
  //with the value thats the return val, but i dont think thats a big issue
  //right now
  rc = fc_getSequenceCoordsPtr(*seq, &data);
  if (rc != FC_SUCCESS){
    fc_printfErrorMessage("failed to get seq data");
  }else {
    ((double*)data)[numStep-1] = timemax;
  }

  return FC_SUCCESS;
}

/** 
 * \ingroup Sequence
 * \brief Given an array of sequences on the same dataset, and new
 *        bounds, creates a new regularly spaced sequence whose range 
 *        is that of the given bounds and whose number of steps is
 *        the max number of steps that the input sequences have on that
 *        range.
 *
 * \description
 *
 *        Given an array of sequences on the same dataset, and new
 *        bounds, creates a new regularly spaced sequence whose range 
 *        is that of the given bounds and whose number of steps is
 *        the max number of steps that the input sequences have on that
 *        range. This is similar to fc_getIntersectingReqularSequence with
 *        input parameter numStep = 0, but where the range is specified, rather
 *        than the entire overlapping range.
 *
 *        Input sequences must be numerical, but dont
 *        have to be same datatype. Input sequences must be 
 *        monotonically increasing. New seq data type is FC_DT_DOUBLE.
 *        The bounds must be included within every input sequence.
 *
 *        This routine is often used as a companion function to 
 *        functions in the timeseries module. This function creates
 *        a common sequence for two or more sequence variables. The
 *        old variables can be put on the new sequence
 *        using fc_LinearInterpolation and then used in the comparison
 *        routines in the timeseries module.
 *
 * \todo
 *     see if can combine this with rc_getIntersectingRegularSequence 
 *
 * \modifications  
 *    - 6/27/05 ACG created.
 */
FC_ReturnCode fc_createIntersectingRangeRegularSequence(
  int numseqs,/**< input - num of seqs in orig sequence array */
  FC_Sequence *seqs,/**< input - orig sequence array */
  double seqmin,/**< input - min val for new sequence */
  double seqmax,/**< input - max val for new sequence */
  char* seqName, /**<input - name of created sequence */
  FC_Sequence *seq /**< output - handle to the new sequence */
  ){

  FC_ReturnCode rc;
  FC_Dataset dset, dsetcomp;
  _FC_DsSlot* dsSlot;
  double closemintime, closemaxtime, timespacing;
  int closeminindex, closemaxindex;
  void *data;
  int numStep, holdnpts;
  int mono;
  FC_DataType datatype;
  int i;

  if (seq)
    *seq = FC_NULL_SEQUENCE;

  if (!seqName || !seq || seqmin > seqmax || numseqs < 1 || !seqs){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }


  fc_printfLogMessage("Creating Intersecting Range regular sequence");

  //checks

  if (!fc_isSequenceValid(seqs[0])){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  };

  rc = fc_getSequenceInfo(seqs[0], &numStep,&datatype);
  if (rc != FC_SUCCESS){
    fc_printfErrorMessage("failed to get seq info");
    return FC_ERROR;
  }

  if (datatype != FC_DT_INT && datatype != FC_DT_FLOAT &&
      datatype != FC_DT_DOUBLE){
    fc_printfErrorMessage("Invalid sequence datatype");
    return FC_INPUT_ERROR;
  }

  rc = fc_getDatasetFromSequence (seqs[0], &dset);
  if (rc != FC_SUCCESS){
    fc_printfErrorMessage("Cannot get dataset from sequence");
    return FC_ERROR;
  }

  dsSlot = _fc_getDsSlot(dset);
  if (dsSlot == NULL){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }


  rc = fc_getSequenceMinMaxMono(seqs[0],NULL,NULL,
                                NULL,NULL,&mono);
  if (rc != FC_SUCCESS){
    fc_printfErrorMessage("failed to get minmax of sequence");
    return FC_ERROR;
  }

  if (mono != 1 && mono != 2){
    fc_printfErrorMessage("invalid monotonicity for this sequence");
    return FC_INPUT_ERROR;
  }

  for (i= 1; i < numseqs; i++){
    if (!fc_isSequenceValid(seqs[i])){
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
      return FC_INPUT_ERROR;
    }
    rc = fc_getDatasetFromSequence (seqs[i], &dsetcomp);
    if (rc != FC_SUCCESS){
      fc_printfErrorMessage("Cannot get dataset from sequence");
      return FC_ERROR;
    }
    
    if (!FC_HANDLE_EQUIV(dset,dsetcomp)){
      fc_printfErrorMessage("Input sequences not on same dataset");
      return FC_INPUT_ERROR;
    } 

    rc = fc_getSequenceInfo(seqs[i], &numStep,&datatype);
    if (rc != FC_SUCCESS){
      fc_printfErrorMessage("failed to get seq info");
      return FC_ERROR;
    }

    if (datatype != FC_DT_INT && datatype != FC_DT_FLOAT &&
        datatype != FC_DT_DOUBLE){
      fc_printfErrorMessage("Invalid sequence datatype");
      return FC_INPUT_ERROR;
    }

    rc = fc_getSequenceMinMaxMono(seqs[i],NULL,NULL,
                                  NULL,NULL,&mono);
    if (rc != FC_SUCCESS){
      fc_printfErrorMessage("failed to get minmax of sequence");
      return FC_ERROR;
    }

    if (mono != 1 && mono != 2){
      fc_printfErrorMessage("invalid monotonicity for this sequence");
      return FC_INPUT_ERROR;
    }
  }

  //finally ok 

  // now get npts and make sure bounds are valid for all exisitng seqs
  // know all are monotonically increasing
  holdnpts = -1;
  for (i= 0 ; i < numseqs; i++){
    int npts;
    rc = fc_getClosestSequenceValue(seqs[i],"<=",seqmin,&closemintime,
                                    &closeminindex);
    if (rc != FC_SUCCESS){
      fc_printfErrorMessage("failed to get closest seq val");
      return FC_ERROR;
    }

    if (closeminindex == -1){
      fc_printfErrorMessage("%s","seqmin is out of bounds for input sequence");
      return FC_ERROR;
    }

    rc = fc_getClosestSequenceValue(seqs[i],">=",seqmax,&closemaxtime,
                                    &closemaxindex);
    if (rc != FC_SUCCESS){
      fc_printfErrorMessage("failed to get closest seq val");
      return FC_ERROR;
    }

    if (closemaxindex == -1){
      fc_printfErrorMessage("%s","seqmax is out of bounds for input sequence");
      return FC_ERROR;
    }
        
    npts = closemaxindex-closeminindex+1;
    if (holdnpts == -1 || npts > holdnpts){
      holdnpts = npts;
    }
  }
  numStep = holdnpts;
  //now have valid range and numpts

  /*
  if (numStep == 1 && !FC_DBL_EQUIV(seqmax,seqmin)){
    //got a problem
    fc_printfErrorMessage("DEVELOEPR ERROR _ HOW CAN THIS BE ?");
    printf("DEVELOEPR ERROR _ HOW CAN THIS BE ?");
    return FC_ERROR;
  }
  */

  //check for zero range sequences
  if (FC_DBL_EQUIV(seqmax,seqmin)){
    //its realy a single pt sequence
    numStep = 1;
    timespacing = 0.0; //actually this value doesnt matter
  } else {
    timespacing= (seqmax-seqmin)/((double)numStep-1);
  }

  rc = fc_createRegularSequence(dset,numStep,seqmin, timespacing,
                                  seqName, seq);
  if (rc != FC_SUCCESS){
    fc_printfErrorMessage("failed to create reg sequence");
    return FC_ERROR;
  }
  
  //just in case there are some roundoff issues, lets explicitly set that
  //last point to the seqmax. (this wont be thread safe since playing 
  //with the value thats the return val, but i dont think thats a big issue
  //right now
  rc = fc_getSequenceCoordsPtr(*seq, &data);
  if (rc != FC_SUCCESS){
    fc_printfErrorMessage("failed to get seq data");
  }else {
    ((double*)data)[numStep-1] = seqmax;
  }

  return FC_SUCCESS;
}

//@}

/** \name Get existing sequences. */
//-------------------------------
//@{

/**
 * \ingroup  Sequence
 * \brief  Get the sequences in a dataset.
 *
 * \description
 *  
 *    This call returns an array of sequences in the dataset.
 *    The caller is responsible for freeing the array of sequence handles.
 *
 * \modifications
 *    - Nancy Collins, Created.
 *    - 11/21/03 RM, added default return numSeq = -1, and 
 *               seqNames = NULL if anything is wrong
 *    - 9/8/04 WSK, changed to return sequence handles instead of sequence
 *               names.
 */
FC_ReturnCode fc_getSequences(
 FC_Dataset dataset,   /**< input - dataset */
 int *numSeq,          /**< output - number of sequences returned */
 FC_Sequence **sequences  /**< output - array of sequences */
) {
  int i;
  _FC_DsSlot* dsSlot;
  _FC_SeqSlot* seqSlot;
  
  // default returns if something goes wrong
  if (numSeq != NULL) 
    *numSeq = -1;
  if (sequences)
    *sequences = NULL;  
  
  // check input
  dsSlot = _fc_getDsSlot(dataset);
  if (dsSlot == NULL || numSeq == NULL || sequences == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Getting sequences from dataset '%s'",
                      dsSlot->header.name);
  
  // get num sequences
  *numSeq = dsSlot->numSeq;
  
  // get sequences
  if (dsSlot->numSeq > 0) {
    *sequences = malloc(dsSlot->numSeq*sizeof(FC_Sequence));
    if (*sequences == NULL) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      return FC_MEMORY_ERROR;
    }
    for (i = 0; i < dsSlot->numSeq; i++) {
      seqSlot = _fc_getSeqSlotFromID(dsSlot->seqIDs[i]);
      _FC_GET_HANDLE((*sequences)[i], seqSlot);
    }
  }
  
  return FC_SUCCESS;
}

/**
 * \ingroup  Sequence
 * \brief  Get the number of sequences in a dataset.
 *
 * \modifications
 *    - 9/28/04 WSK, created
 */
FC_ReturnCode fc_getNumSequence(
 FC_Dataset dataset,   /**< input - dataset */
 int *numSeq           /**< output - number of sequences returned */
) {
  _FC_DsSlot* dsSlot;
  
  // default returns if something goes wrong
  if (numSeq != NULL) 
    *numSeq = -1;
  
  // check input
  dsSlot = _fc_getDsSlot(dataset);
  if (dsSlot == NULL || numSeq == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Getting number sequences from dataset '%s'",
                      dsSlot->header.name);
  
  // get num sequences
  *numSeq = dsSlot->numSeq;
  
  return FC_SUCCESS;
}

/**
 * \ingroup  Sequence
 * \brief  Get all sequences with a given name from a dataset.
 *
 * \description
 *  
 *    Returns all sequences on the dataset with the specified name. 
 *
 * \modifications  
 *   - Nancy Collins  Created.
 *   - 2003-MAY-30  W Koegler  Moved getting sequence meta data from
 *       _fc_setupDsSlot() to here.
 *   - 11/25/03 RM, added default return sequence = \ref FC_NULL_SEQUENCE 
 *   - 12/01/03 WSK Changed so it removes an uncommitted sequence
 *       from the owning dataset's seqID list.
 *   - 12/02/03 RM, added checking for NULL sequence
 *   - 12/12/03 RM, eliminated _fc_lookup_seqSlot() and _fc_namematch(),
 *       those input checking changed and added into the body of this function
 *   - 9/8/03 WSK, changed name and moved to dataset (to be next to
 *       fc_getSequences) 
 *   - 9/22/06 ACG returns array of matching sequences
 */
FC_ReturnCode fc_getSequenceByName(
 FC_Dataset dataset,    /**< input - handle for this dataset */
 char *seqName,         /**< input - name of sequence to get */
 int *numSequences,     /**< output - number of matching sequences */
 FC_Sequence **sequences  /**< output - handle for this sequence */ 
) {
  int i;
  _FC_DsSlot* dsSlot;
  _FC_SeqSlot *seqSlot;
  FC_Sequence sequence;

  FC_Sequence *lmatch = NULL;
  int numlmatch, maxnumlmatch;
  
  // default values returned if anything goes wrong
  if (numSequences)
    *numSequences = -1;
  if (sequences) 
    *sequences = NULL;
  
  // check input
  dsSlot = _fc_getDsSlot(dataset);
  if (dsSlot == NULL || seqName == NULL || !sequences || !numSequences ) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Getting sequence '%s'", seqName);

  numlmatch = 0;
  maxnumlmatch = 0;
  
  // look for existing slots with the matching name
  for (i = 0; i < dsSlot->numSeq; i++) {
    seqSlot = _fc_getSeqSlotFromID(dsSlot->seqIDs[i]);
    if (!strcmp(seqSlot->header.name, seqName)) {
      _FC_GET_HANDLE(sequence, seqSlot);

      if (lmatch == NULL){
	lmatch = (FC_Sequence*)malloc(sizeof(FC_Sequence));
	if (lmatch == NULL){
	  fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
	  return FC_MEMORY_ERROR;
	}
	numlmatch = 1;
	maxnumlmatch = 1;
	lmatch[0] = sequence;
      }else{
	if(numlmatch == maxnumlmatch){
	  FC_Sequence *temp;
	  temp = (FC_Sequence*)realloc(lmatch,(2*numlmatch)*sizeof(FC_Sequence));
	  if (temp == NULL){
	    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
	    return FC_MEMORY_ERROR;
	  }
	  maxnumlmatch*=2;
	  lmatch = temp;
	}
	lmatch[numlmatch] = sequence;
	numlmatch++;
      }
    }
  }

  *numSequences = numlmatch;
  *sequences = lmatch;

  return FC_SUCCESS;
}

//@}

/** \name Change the name of a sequence. */
//-------------------------------
//@{

/**
 * \ingroup Sequence
 * \brief  Change the name of a sequence.
 *
 * \modifications
 *   - 8/24/2005 WSD. Created.
 */
FC_ReturnCode fc_changeSequenceName(
  FC_Sequence sequence, /**< Input - the sequence. */
  char* newName       /**< Input - new name for the sequence. */
) {
  _FC_SeqSlot* seqSlot;

  // check input
  seqSlot = _fc_getSeqSlot(sequence);
  if (seqSlot == NULL || newName == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // log message
  fc_printfLogMessage("Changing name of sequence from '%s' to '%s'", 
                      seqSlot->header.name, newName);
  
  // Do it
  return _fc_setSlotHeaderName(&seqSlot->header, newName);
}

//@}

/** \name Release/delete sequences. */
//-------------------------------
//@{

/**
 * \ingroup  Sequence
 * \brief  Attempt to minimize sequence's memory usage.
 *
 * \description
 *  
 *    Call to try and release un-needed memory in a  sequence. 
 *    If the sequence has been saved to disk, the large data is 
 *    released. If the sequence has not been saved to disk, this 
 *    will do nothing.
 *
 * \modifications  
 *    - 04/19/04 WSK Created.
 */
FC_ReturnCode fc_releaseSequence(
 FC_Sequence sequence       /**< input - sequence handle */
) {
  _FC_SeqSlot* seqSlot;
  
  // special case: releasing a null handle is not an error
  if ( FC_HANDLE_EQUIV(sequence, FC_NULL_SEQUENCE)) {
    fc_printfWarningMessage("Releasing FC_NULL_SEQUENCE");
    return FC_SUCCESS;
  }

  // check input
  seqSlot = _fc_getSeqSlot(sequence);
  if (seqSlot == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Cleaning up sequence '%s'", seqSlot->header.name);

  // remove big data if this slot is committed
  if (seqSlot->header.committed == 1) {
    free(seqSlot->coords);
    seqSlot->coords = NULL;
  }

  return FC_SUCCESS;
}

/**
 * \ingroup  Sequence
 * \brief  Delete a sequence.
 *
 * \description
 * 
 *    Call when this sequence is no longer needed. Deleting a NULL
 *    sequence handle is not an error, but trying to delete a sequence if any
 *    meshes have a sequence variable on it is an error.  So, delete meshes
 *    before deleting sequences.
 *
 * \modifications  
 *    - Nancy Collins, Created.
 *    - 11/25/03 RM, added more input checking, if sequence is 
 *      FC_NULL_SEQUENCE) then do nothing and return 
 *      success (issue warning if fc_getLibraryVerbosity()). 
 *    - 12/01/03 WSK Changed so it removes an uncommitted sequence
 *      form the owning dataset's seqID list.
 *    - 9/9/04 WSK changed name to delete. Changed behavior so that
 *      delete always deletes no matter whether on disk or not.
 */
FC_ReturnCode fc_deleteSequence(
 FC_Sequence sequence     /**< input - sequence to be deleted */
) {
  int i, j;
  _FC_SeqSlot* seqSlot;
  _FC_DsSlot* dsSlot;
  
  // special case: deleting a null handle is not an error
  if ( FC_HANDLE_EQUIV(sequence, FC_NULL_SEQUENCE)) {
    fc_printfWarningMessage("Deleting FC_NULL_SEQUENCE");
    return FC_SUCCESS;
  }
  
  // check input
  seqSlot = _fc_getSeqSlot(sequence);
  if (seqSlot == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Deleting sequence '%s'", seqSlot->header.name);

  // setup
  dsSlot = _fc_getDsSlot(seqSlot->ds);

  // special case: don't delete if any mesh has a reference to
  for (i = 0; i < dsSlot->numMesh; i++) {
    _FC_MeshSlot* meshSlot = _fc_getMeshSlotFromID(dsSlot->meshIDs[i]);
    for (j = 0; j < meshSlot->numSeqVar; j++) {
      _FC_VarSlot* varSlot = _fc_getVarSlotFromID(meshSlot->seqVarIDs[j][0]);
      if (FC_HANDLE_EQUIV(sequence, varSlot->sequence)) {
        fc_printfErrorMessage("Can't delete sequence '%s' because mesh '%s' "
                              "has a reference to it", seqSlot->header.name,
                              meshSlot->header.name);
        return FC_ERROR;
      }
    }
  }
  
  // remove reference from parent dataset
  if (dsSlot->numSeq == 1) {
    free(dsSlot->seqIDs);
    dsSlot->seqIDs = NULL;
    dsSlot->numSeq = 0;
  }
  else {
    for (i = 0; i < dsSlot->numSeq; i++) {
      if (dsSlot->seqIDs[i] == seqSlot->header.slotID)
        break;
    }
    for (j = i+1; j < dsSlot->numSeq; j++) {
      dsSlot->seqIDs[j-1] = dsSlot->seqIDs[j];
    }
    dsSlot->numSeq--;
  }

  // clear & delete slot
  return _fc_deleteSeqSlot(sequence);
}

//@}

/** \name Get sequence metadata. */
//-------------------------------
//@{

/**
 * \ingroup  Sequence
 * \brief Check that the handle refers to a valid sequence.
 *
 * \description
 *
 *    Returns 1 (true) is the handle is valid, and 0 (false) if it is not.
 *
 * \modifications  
 *    - 05/25/04 WSK Created.
 */
int fc_isSequenceValid(
  FC_Sequence sequence  /**< input - sequence handle */
) {
  return _fc_isHandleValid(sequence.slotID, sequence.uID, seqTableSize,
                           (_FC_SlotHeader**)seqTable);
}

/**
 * \ingroup  Sequence
 * \brief  Return name of a sequence.
 *
 * \description
 *  
 *    Return the name associated with a sequence handle. Spaced
 *    is allocated for the string so it needs to be freed by the
 *    user when finished.
 *
 * \modifications  
 *    - Aug, 9, 2002.  W Koegler  Created
 *    - 11/25/03 RM, added checking for NULL dataset handle, 
 *         default return seqname = NULL if anything is wrong
 */
FC_ReturnCode fc_getSequenceName(
  FC_Sequence sequence,   /**< input - sequence handle */
  char **seqname          /**< output - sequence name, space allocated here
                             and must be freed by user when finished */
) {
  char *cp;
  _FC_SeqSlot* seqSlot;
  
  // default returns if something goes wrong
  if (seqname)
    *seqname = NULL;  
  
  // check input
  seqSlot = _fc_getSeqSlot(sequence);
  if (seqSlot == NULL || seqname == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Getting name for sequence '%s'", seqSlot->header.name);

  // copy name
  cp = seqSlot->header.name;
  *seqname = malloc(strlen(cp) + 1);
  if (*seqname == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  strcpy(*seqname, cp);
  
  return FC_SUCCESS;
}

/**
 * \ingroup  Sequence
 * \brief  Get the parent dataset of the sequence.
 *
 * \modifications  
 *    - Nancy Collins, Created.
 *    - 11/25/03 RM, added default return dataset = \ref FC_NULL_DATASET
 *    - 12/22/03 Rm, added checking for input/output dataset == NULL
 */
FC_ReturnCode fc_getDatasetFromSequence(
  FC_Sequence sequence,   /**< input - sequence handle */
  FC_Dataset *dataset     /**< output - dataset handle */
) {
  _FC_SeqSlot* seqSlot;
  
  // default values returned if anything goes wrong
  if (dataset != NULL) 
    *dataset = FC_NULL_DATASET;
  
  // check input
  seqSlot = _fc_getSeqSlot(sequence);
  if (seqSlot == NULL || dataset == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Getting dataset from sequence '%s'", 
                      seqSlot->header.name);
  
  // return requested value
  *dataset = seqSlot->ds;
  
  return FC_SUCCESS;
}

/**
 * \ingroup  Sequence
 * \brief  Get information about a sequence.
 *
 * \description
 *  
 *    Return the number of steps in a sequence, and the datatype of the
 *    sequences coordinates.
 *
 * \modifications 
 *    - Nancy Collins Created.
 *    - 2003-APR-18  W Koegler  Finishing sequence implementation
 *    - 11/25/03 RM, added default return numStep = -1, datatype = 
 *      \ref FC_DT_UNKNOWN
 */
FC_ReturnCode fc_getSequenceInfo(
  FC_Sequence sequence,   /**< input - sequence handle */
  int *numStep,           /**< output - number of steps in the sequence */
  FC_DataType *datatype   /**< output - data type of sequence coordinates */
) {
  _FC_SeqSlot* seqSlot;
  
  // default values returned if anything goes wrong
  if (numStep != NULL) 
    *numStep = -1 ;
  if (datatype != NULL)
    *datatype = FC_DT_UNKNOWN;
  
  // check input
  seqSlot = _fc_getSeqSlot(sequence);
  if (seqSlot == NULL || (!numStep && !datatype)) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Getting info for sequence '%s'", 
                      seqSlot->header.name);
  
  // fill in requested values
  if (numStep)
    *numStep = seqSlot->numStep;
  if (datatype)
    *datatype = seqSlot->datatype;
  
  return FC_SUCCESS;
}

/**
 * \ingroup  Sequence
 * \brief  Get number of steps in a sequence.
 *
 * \modifications 
 *    - 9/21/04 WSK, created.
 */
FC_ReturnCode fc_getSequenceNumStep(
  FC_Sequence sequence,   /**< input - sequence handle */
  int *numStep            /**< output - number of steps in the sequence */
) {
  _FC_SeqSlot* seqSlot;
  
  // default values returned if anything goes wrong
  if (numStep != NULL) 
    *numStep = -1 ;
  
  // check input
  seqSlot = _fc_getSeqSlot(sequence);
  if (seqSlot == NULL || !numStep) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Getting numStep for sequence '%s'", 
                      seqSlot->header.name);
  
  // fill in requested values
  if (numStep)
    *numStep = seqSlot->numStep;
  
  return FC_SUCCESS;
}

/**
 * \ingroup  Sequence
 * \brief  Get the data type of a sequence's coordinates.
 *
 * \modifications 
 *    - 9/21/04 WSK, created.
 */
FC_ReturnCode fc_getSequenceDataType(
  FC_Sequence sequence,   /**< input - sequence handle */
  FC_DataType *datatype   /**< output - data type of sequence coordinates */
) {
  _FC_SeqSlot* seqSlot;
  
  // default values returned if anything goes wrong
  if (datatype != NULL)
    *datatype = FC_DT_UNKNOWN;
  
  // check input
  seqSlot = _fc_getSeqSlot(sequence);
  if (seqSlot == NULL || !datatype) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Getting info for sequence '%s'", 
                      seqSlot->header.name);
  
  // fill in requested value
  if (datatype)
    *datatype = seqSlot->datatype;
  
  return FC_SUCCESS;
}

//@}

/** \name Get sequence big data. */
//-------------------------------
//@{

/**
 * \ingroup  Sequence
 * \brief  Get sequence coordinates
 *
 * \description
 *  
 *    Return a pointer to the sequence's coordinate array. This
 *    must be treated as read-only.
 *
 * \modifications 
 *   - Nancy Collins Created.
 *   - 2003-APR-18  W Koegler  Finishing sequence implementation
 *   - 11/25/03 RM, added checking for NULL sequence handle, 
 *         default return *coords = NULL, issue an error if 
 *         fc_getLibraryVerbosity()
 */
FC_ReturnCode fc_getSequenceCoordsPtr(
  FC_Sequence sequence,   /**< input - sequence handle */
  void **coords_p         /**< output - pointer to coord array */
  ) {
  FC_ReturnCode rc;
  _FC_SeqSlot* seqSlot;
  
  // default values returned if anything goes wrong
  if (coords_p)
    *coords_p = NULL;
  
  // check input
  seqSlot = _fc_getSeqSlot(sequence);
  if (seqSlot == NULL || coords_p == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Getting coords for sequence '%s'", 
                      seqSlot->header.name);
  
  // lazy fetching of sequence coords done here
  if (seqSlot->coords == NULL && seqSlot->numStep > 0) {
    rc = _fc_readSequenceCoords(sequence);
    if (rc != FC_SUCCESS) {
      fc_printfErrorMessage("Lazy loading of coords of sequence '%s' failed",
                            seqSlot->header.name);
      return rc;
    }
  }
  
  // set return values
  *coords_p = seqSlot->coords;
  
  return FC_SUCCESS;
}


/**
 * \ingroup  Sequence
 * \brief  Get a copy of the sequence's coords as the requested data type
 *
 * \description
 *  
 *    Returns the sequence's data as the type requested by the user.
 *    The user is responsible for freeing the array. The routine casts
 *    to the requested user type so there is a potential for lost information
 *    (e.g. if you cast a floating point number to int it is truncated, not
 *    rounded).
 *
 *    If you are worried about memory or casting, you should use
 *    fc_getVariableDataPtr() to get a pointer to the raw data array.
 *
 * \modifications
 *   - 06/09/2006 WSD Created.
 */
FC_ReturnCode fc_getSequenceCoordsAsDataType(
  FC_Sequence sequence,   /**< input - sequence handle */
  FC_DataType dataType,   /**< input - the data type for the returned coords */
  void **coords           /**< output - pointer to coord array */
  ) {
  FC_ReturnCode rc;
  int i;
  int numStep;
  FC_DataType dataType_orig;
  void* coords_orig;
  
  // default values returned if anything goes wrong
  if (coords)
    *coords = NULL;
  
  // check input
  if (!fc_isSequenceValid(sequence) || !fc_isDataTypeValid(dataType) ||
      dataType == FC_DT_UNKNOWN || coords == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Getting coords for sequence '%s' as type %s", 
		      _fc_getSeqSlot(sequence)->header.name,
		      fc_getDataTypeText(dataType));
  
  // get sequence info and coords
  rc = fc_getSequenceInfo(sequence, &numStep, &dataType_orig);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_getSequenceCoordsPtr(sequence, &coords_orig);
  if (rc != FC_SUCCESS)
    return rc;

  // allocate return array
  *coords = malloc(numStep*fc_sizeofDataType(dataType));
  if (!*coords) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  
  // special case - the same data type
  if (dataType == dataType_orig) {
    memcpy(*coords, coords_orig, numStep*fc_sizeofDataType(dataType));
    return FC_SUCCESS;
  }

  // Fun!
  if (dataType == FC_DT_CHAR) {
    if (dataType_orig == FC_DT_INT) 
      for (i = 0; i < numStep; i++)
        ((char*)(*coords))[i] = ((int*)coords_orig)[i];
    else if (dataType_orig == FC_DT_FLOAT) 
      for (i = 0; i < numStep; i++)
        ((char*)(*coords))[i] = ((float*)coords_orig)[i];
    else if (dataType_orig == FC_DT_DOUBLE) 
      for (i = 0; i < numStep; i++)
        ((char*)(*coords))[i] = ((double*)coords_orig)[i];
  }
  else if (dataType == FC_DT_INT) {
    if (dataType_orig == FC_DT_CHAR) 
      for (i = 0; i < numStep; i++)
        ((int*)(*coords))[i] = ((char*)coords_orig)[i];
    else if (dataType_orig == FC_DT_FLOAT) 
      for (i = 0; i < numStep; i++)
        ((int*)(*coords))[i] = ((float*)coords_orig)[i];
    else if (dataType_orig == FC_DT_DOUBLE) 
      for (i = 0; i < numStep; i++)
        ((int*)(*coords))[i] = ((double*)coords_orig)[i];
  }
  else if (dataType == FC_DT_FLOAT) {
    if (dataType_orig == FC_DT_CHAR) 
      for (i = 0; i < numStep; i++)
        ((float*)(*coords))[i] = ((char*)coords_orig)[i];
    else if (dataType_orig == FC_DT_INT) 
      for (i = 0; i < numStep; i++)
        ((float*)(*coords))[i] = ((int*)coords_orig)[i];
    else if (dataType_orig == FC_DT_DOUBLE) 
      for (i = 0; i < numStep; i++)
        ((float*)(*coords))[i] = ((double*)coords_orig)[i];
  }
  else if (dataType == FC_DT_DOUBLE) {
    if (dataType_orig == FC_DT_CHAR) 
      for (i = 0; i < numStep; i++)
        ((double*)(*coords))[i] = ((char*)coords_orig)[i];
    else if (dataType_orig == FC_DT_INT) 
      for (i = 0; i < numStep; i++)
        ((double*)(*coords))[i] = ((int*)coords_orig)[i];
    else if (dataType_orig == FC_DT_FLOAT) 
      for (i = 0; i < numStep; i++)
        ((double*)(*coords))[i] = ((float*)coords_orig)[i];
  }

  return FC_SUCCESS;
}

//@}

/** \name Print sequence. */
//-------------------------------
//@{

/**
 * \ingroup  Sequence
 * \brief Print sequence meta data to stdout, and optionally the coordinates.
 *
 * \modifications 
 *    - 11/20/03 WSK Created.
 *    - 10/31/07 ACG optional precision flag
 */
FC_ReturnCode fc_printSequence(
  FC_Sequence sequence, /**< input - sequence */
  char* label,          /**< input - label to prepend variable name with, 
                                     can be NULL */
  int print_coords      /**< input - 1 = print the coordinates. 
			   2 = print the in exponential notation with 5 digits of precision */
) {
  FC_ReturnCode rc;
  int i;
  _FC_SeqSlot* seqSlot;
  int numStep;
  FC_DataType dataType;
  void* dbuf;

  // special case: printing a null handle is not an error
  if (FC_HANDLE_EQUIV(sequence, FC_NULL_SEQUENCE)) {
    printf("Sequence: FC_NULL_SEQUENCE\n");
    fflush(NULL);
    return FC_SUCCESS;
  }

  // check input
  seqSlot = _fc_getSeqSlot(sequence);
  if (seqSlot == NULL) {
    fc_printfWarningMessage("Deleting FC_NULL_DATASET");
    return FC_INPUT_ERROR;
  }

  // log message
  fc_printfLogMessage("Printing sequence '%s'", seqSlot->header.name);

  // print label & var name
  if (label)
    printf("%s: ", label);
  else
    printf("Sequence: ");
  printf("'%s'\n", seqSlot->header.name);
  fflush(NULL);

  // print meta data
  rc = fc_getSequenceInfo(sequence, &numStep, &dataType);
  if (rc != FC_SUCCESS) {
    fc_printfErrorMessage("Failed to get info for sequence '%s'",
                          seqSlot->header.name);
    return rc;
  }
  printf("%snumStep = %d, dataType = %s\n", INDENT_STR, numStep,
         fc_getDataTypeText(dataType));
  fflush(NULL);

  // print coordinates
  if (print_coords) {
    rc = fc_getSequenceCoordsPtr(sequence, &dbuf);
    if (rc != FC_SUCCESS) {
      fc_printfErrorMessage("Failed to get coords for sequence '%s'",
                          seqSlot->header.name);
      return rc;
    }
    
    // print the coordinates
    printf("%sCoordinates:\n", INDENT_STR);
    fflush(NULL);
    switch (dataType) {
    case FC_DT_CHAR: 
      {
        char* cast_buf = (char *)dbuf;
        for (i = 0; i < numStep; i++)
          printf("%s%s%d:  %c\n", INDENT_STR, INDENT_STR, i, cast_buf[i]);
      }
      break;
    case FC_DT_INT:
      {
        int* cast_buf = (int *)dbuf;
        for (i = 0; i < numStep; i++) 
          printf("%s%s%d:  %d\n", INDENT_STR, INDENT_STR, i, cast_buf[i]);
      }
      break;
    case FC_DT_FLOAT:
      {
        float* cast_buf = (float *)dbuf;
        for (i = 0; i < numStep; i++) 
          printf("%s%s%d:  %.6g\n", INDENT_STR, INDENT_STR, i, cast_buf[i]);
      }
      break;
    case FC_DT_DOUBLE:
      {
        double* cast_buf = (double *)dbuf;
        for (i = 0; i < numStep; i++) 
	  switch(print_coords){
	  case 2:
	    printf("%s%s%d:  %5leg\n", INDENT_STR, INDENT_STR, i, cast_buf[i]);
	    break;
	  default:
	    printf("%s%s%d:  %.15g\n", INDENT_STR, INDENT_STR, i, cast_buf[i]);
	    break;
	  }
      }
      break;
    case FC_DT_UNKNOWN:
      printf("Cannot print data type '%s'\n", fc_getDataTypeText(dataType));
    }
    fflush(NULL);
  }

  return FC_SUCCESS;
}

//@}

/**
 * \ingroup  PrivateSequence
 * \brief  Initialize a sequence slot
 *
 * \modifications  
 *   APR-30-2003  W Koegler  Created.
 */
void _fc_initSeqSlot(
  _FC_SeqSlot* seqSlot    /**< input - the qtable slot to initialize */
) {
  if (seqSlot != NULL) {
    _fc_initSlotHeader(&seqSlot->header);
    seqSlot->ds = FC_NULL_DATASET;
    seqSlot->numStep = 0;
    seqSlot->datatype = FC_DT_UNKNOWN;
    seqSlot->coords = NULL;
  }   
}

/**
 * \ingroup  PrivateSequence
 * \brief  Release non-necessary resources in a sequence slot.
 *
 * \description
 *  
 *    Releases the coordinates.
 *
 * \modifications  
 *   2003-APR-17  W Koegler  Created
 */
void _fc_releaseSeqSlot(
  _FC_SeqSlot* seqSlotp    /**< input - qtable pointer */
) {

  // release big data
  if (seqSlotp->coords)
    free(seqSlotp->coords);
  seqSlotp->coords = NULL;
    
  return;
}

/**
 * \ingroup  PrivateSequence
 * \brief  Clear a sequence slot.
 *
 * \description
 *  
 *    Releases all dynamically allocated resources in a sequence,
 *    and reinitializes all members.
 *
 * \modifications  
 *   - 2003-APR-17  W Koegler  Created
 *   - 2003-APR-31  W Koegler  Made more comprehensive
 */
void _fc_clearSeqSlot(
 _FC_SeqSlot* seqSlotp       /**< input - qtable pointer */
) {
  if (seqSlotp != NULL) {
    // release big data
    _fc_releaseSeqSlot(seqSlotp);

    // free other dynamic arrays, if any

    // clear table header
    _fc_clearSlotHeader(&seqSlotp->header);

    // reinit
    _fc_initSeqSlot(seqSlotp);
  }
}

/**
 * \ingroup  PrivateSequence
 * \brief  Get the size of the sequence table.
 *
 * \modifications  
 *   - 05/21/04 WSK Created.
 */
int _fc_getSeqTableSize() {
  return seqTableSize;
}

/**
 * \ingroup  PrivateSequence
 * \brief   Get an empty sequence slot from the sequence table.
 *
 * \description
 *  
 *    Will return a new slot or NULL. A NULL means a memory error occurred.
 *
 * \modifications  
 *    - 2003-MAY-01  W Koegler  Created
 *    - 12/08/03 WSK Changed because changed XTables from arrays of slots
 *      to arrays of handles to slots. The _fc_new_XSlot() routines do
 *      not call a common helper routine anymore.
 *    - 5/11/04 WSK  Moved common code from _fc_getNewDsSlot(), 
 *      _fc_getNewSeqSlot(), etc to _fc_getNewSlot().
 */
_FC_SeqSlot* _fc_getNewSeqSlot() {
  _FC_SeqSlot* seqSlot;

  seqSlot = (_FC_SeqSlot*) _fc_getNewSlot(&seqTableSize, 
					  (_FC_SlotHeader***)&seqTable, 
					  &seqOpenSlots, sizeof(_FC_SeqSlot));
  _fc_initSeqSlot(seqSlot);

  return seqSlot;
}

/**
 * \ingroup  PrivateSequence
 * \brief Get a sequence slot from a sequence handle.
 *
 * \description
 *g 
 *    The slot id of the sequence handle is used to find the
 *    the slot in the sequence table. The identity is then
 *    verified by checking that the slot and the handle's
 *    uIDs are the same. On success, a pointer to the slot 
 *    is returned. Otherwise a NULL pointer is returned.
 *
 * \modifications 
 *   - 2003-OCT-19  WSK  Created to replace using _fc_checkhandle and GET_ID 
 *       for better type checking (using pointer's to slots is less bug
 *       prone than lot's of xxxTable[xxxID] references).
 *   - 05/25/04 WSK Changed to call new _fc_isHandleValid common to all tables.
 */
_FC_SeqSlot* _fc_getSeqSlot(
  FC_Sequence sequence    /**< input - sequence */
) {
  if (_fc_isHandleValid(sequence.slotID, sequence.uID, seqTableSize,
                        (_FC_SlotHeader**)seqTable))
    return seqTable[sequence.slotID];
  else
    return NULL;
}

/**
 * \ingroup  PrivateSequence
 *
 * \brief Get a sequence slot based on ID.
 *
 * \description
 * 
 *    The _FC_SeqSlot at the requested slotID is returned. This should
 *    only be used when you really know that it is the slot you want
 *    because there is no checking of the uID like with _fc_getSeqSlot.
 *    The only checking is that the slotID is exists in the table and
 *    that the uID is > 0 (i.e. it is not empty).
 *
 * \modifications 
 *    - 11/27/03 WSK Created so that external routines can index into
 *      tables without knowing about them.
 *    - 05/25/04 WSK Changed to call new _fc_getSlot() common to all tables.
 */
_FC_SeqSlot* _fc_getSeqSlotFromID(
  int seqID
) {
  return (_FC_SeqSlot*) _fc_getSlot(seqID, seqTableSize, 
                                    (_FC_SlotHeader**)seqTable);
}

/**
 * \ingroup  PrivateSequence
 * \brief Delete the sequence slot associated with the sequence handle.
 *
 * \description 
 *    
 *    This routine clears the dynamic data in the slot and then
 *    deletes the slot from the table.
 *
 * \modifications
 *   - 12/20/2005 WSD. Created.
 */
FC_ReturnCode _fc_deleteSeqSlot(
  FC_Sequence sequence    /**< input - sequence */
) {
  _FC_SeqSlot* seqSlot;
  
  seqSlot = _fc_getSeqSlot(sequence);
  if (seqSlot == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // do it
  _fc_clearSeqSlot(seqSlot);
  return _fc_deleteSlot(seqSlot->header.slotID, seqTableSize, 
			(_FC_SlotHeader**)seqTable, &seqOpenSlots);
}

/**
 * \ingroup  PrivateSequence
 * \brief  Print the contents of the Sequence Table (\ref seqTable) to stderr.
 *
 * \description
 *
 *    Prints the contents of the Sequence Table in a human readable form.
 *    It takes a string to use as a label for this print out of the table.
 *    Pass NULL if you don't want a label.
 *
 * \modifications 
 *    - Nancy Collins Created.
 *    - 2003-NOV-13  WSK  Fixed up.
 */
void _fc_printSeqTable(
  char *label  /**< Input - label for this version the table (without a
                           trailing \\n) pass NULL if no label is desired */
) {
  int i;
  _FC_SeqSlot* seqSlot;
  
  // print table heading
  fprintf(stderr, TABLE_DIVIDER_STRING);
  fprintf(stderr, "Sequence Table (seqTable):\n");
  if (label != NULL)
    fprintf(stderr, "%s\n", label);
  fprintf(stderr, TABLE_DIVIDER_STRING);
  fflush(NULL);

  // We're done if it's empty
  if (seqTableSize < 1) {
    fprintf(stderr, "(empty)\n");
    fprintf(stderr, "\n");
    fflush(NULL);
    return;
  }

  // print contents for each slot
  for (i = 0; i < seqTableSize; i++) {
    seqSlot = seqTable[i];
    fprintf(stderr, "%3d: %s\n", i, (seqSlot) ? "exists" : "NULL");
    if (seqSlot == NULL)
      continue;
    _fc_printSlotHeader(seqSlot->header);
    fprintf(stderr, "     FC_Dataset = { %d, %d }\n",
            seqSlot->ds.slotID, seqSlot->ds.uID);
    fprintf(stderr, "     numStep = %d, datatype = %s, coords = %s\n", 
            seqSlot->numStep, fc_getDataTypeText(seqSlot->datatype),
            (seqSlot->coords) ? "exists" : "NULL");
  }
  fprintf(stderr, "\n");
  fflush(NULL);

  return;
}

/**
 * \ingroup  PrivateSequence
 * \brief  Free all entries in the sequence table
 *
 * \description
 *
 *    This should only be called after making sure that
 *    all slots have been cleared.
 *
 * \modifications  
 *   - 05/21/04 WSK Created.
 */
void _fc_freeSeqTable() {
  int i;
  for (i = 0; i < seqTableSize; i++) {
    _fc_clearSeqSlot(seqTable[i]);
    free(seqTable[i]);
  }
  free(seqTable);
  fc_freeSortedIntArray(&seqOpenSlots);
  seqTableSize = 0;
  seqTable = NULL;
}

