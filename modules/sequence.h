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
 * \file sequence.h
 * \brief Public declarations for the \ref Sequence Module.
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/modules/sequence.h,v $
 * $Revision: 1.21 $ 
 * $Date: 2006/09/22 21:38:57 $
 */

#ifndef _FC_SEQUENCE_H_
#define _FC_SEQUENCE_H_

#ifdef __cplusplus
extern "C" {
#endif

// create new sequence from scratch
//-------------------------------
FC_ReturnCode fc_createSequence(FC_Dataset dataset, char *seqname, 
		          FC_Sequence *sequence);
FC_ReturnCode fc_setSequenceCoords(FC_Sequence sequence, int numStep,
                          FC_DataType datatype, void* coords);
FC_ReturnCode fc_setSequenceCoordsPtr(FC_Sequence sequence, int numStep,
                          FC_DataType datatype, void* coords_p);


// other ways to get new sequences
//-------------------------------
FC_ReturnCode fc_copySequence(FC_Sequence src_seq, FC_Dataset dest_ds, 
                          char* newSeqName, FC_Sequence* new_seq);
FC_ReturnCode fc_shiftAndScaleSequence(FC_Sequence src_seq, FC_Dataset
					 dataset, double shift, double scale,
					 char* newSeqName, 
					 FC_Sequence* newseq);
FC_ReturnCode fc_createRegularSequence(FC_Dataset dataset,
		       int numStep, double offset, double stepsize,
		       char* seq_name, FC_Sequence *seq);
FC_ReturnCode fc_createIntersectingRegularSequence(
			   int numseqs, FC_Sequence* seqs,
			   int numStep,
			   char* seqName,
			   FC_Sequence *sequence);
FC_ReturnCode fc_createIntersectingRangeRegularSequence(
			   int numseqs, FC_Sequence* seqs,
			   double seqmin, double seqmax,
			   char* seqName,
			   FC_Sequence *sequence);

// get sequences
//-------------------------------
FC_ReturnCode fc_getSequences(FC_Dataset dataset, int *numSeq, 
			      FC_Sequence **sequences);
FC_ReturnCode fc_getNumSequence(FC_Dataset dataset, int *numSeq);
FC_ReturnCode fc_getSequenceByName(FC_Dataset dataset, char *seqName, 
				   int *numSequences,
				   FC_Sequence **sequences);

// change the name of a sequence
//-------------------------------
FC_ReturnCode fc_changeSequenceName(FC_Sequence sequence, char* newName);

// release/delete sequence
//-------------------------------
FC_ReturnCode fc_releaseSequence(FC_Sequence sequence);
FC_ReturnCode fc_deleteSequence(FC_Sequence sequence);

// get Sequence metadata
//-------------------------------
int fc_isSequenceValid(FC_Sequence sequence);
FC_ReturnCode fc_getSequenceName(FC_Sequence sequence, char** seqName);
FC_ReturnCode fc_getDatasetFromSequence(FC_Sequence sequence, 
                          FC_Dataset *dataset);
FC_ReturnCode fc_getSequenceInfo(FC_Sequence sequence, int *numStep, 
                          FC_DataType *dataType);
FC_ReturnCode fc_getSequenceNumStep(FC_Sequence sequence, int* numStep);
FC_ReturnCode fc_getSequenceDataType(FC_Sequence sequence, 
                          FC_DataType *dataType);

// get big data
//-------------------------------
FC_ReturnCode fc_getSequenceCoordsPtr(FC_Sequence sequence, void **coords_p);
FC_ReturnCode fc_getSequenceCoordsAsDataType(FC_Sequence sequence, 
					     FC_DataType datatype, void **data_p);
// Print
//-------------------------------
FC_ReturnCode fc_printSequence(FC_Sequence sequence, char* label, 
                          int print_coords);


#ifdef __cplusplus
}
#endif

#endif // _FC_SEQUENCE_H_
