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
 * \file checksequence.c
 * \brief Unit tests for \ref Sequence module.
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/unittest/checksequence.c,v $
 * $Revision: 1.30 $ 
 * $Date: 2006/09/22 21:43:56 $
 *
 * \modifications
 *    9/10/04 WSK, split off of checklibrary
 */

#include <stdlib.h>
#include <string.h>
#include <check.h>
#include "fc.h"
#include "fcP.h"
#include "checkall.h"

// **** test fixtures

static void sequence_setup(void) {
  FC_ReturnCode rc;

  if (isForking) {
    fc_setLibraryVerbosity(fc_messages);
    //    fc_setLibraryVerbosity(FC_ERROR_MESSAGES);
    rc = fc_initLibrary();
    fail_unless(rc == FC_SUCCESS, "test aborted: failed to init library");
  }
}

static void sequence_teardown(void) {
  FC_ReturnCode rc;

  if (isForking) {
    rc = fc_finalLibrary();
    fail_unless(rc == FC_SUCCESS, "failed to final library at end of test");
  }
}

// Test that creating and deleting slots modifies table as expected
// Deleting should open up slots that get reused.
START_TEST(slot_new_delete)
{
  FC_ReturnCode rc;
  int i;
  int num = 7;
  // uID_incs = diff between current uID and startID
  int uID_start, uID_incs[7] = { 0, 7, 2, 8, 4, 9, 6 };
  FC_Dataset dataset;
  FC_Sequence sequences[7];
  
  // create parent datset
  rc = fc_createDataset("fake", &dataset);
  fail_unless(rc == FC_SUCCESS, "failed to create dataset");

  // Create num
  for (i = 0; i < num; i++) {
    rc = fc_createSequence(dataset, "fake", &sequences[i]);
    fail_unless(rc == FC_SUCCESS, "failed to create sequence");
  }
  uID_start = sequences[0].uID;

  // Delete "odd" slots
  for (i = 1; i < num; i+=2) {
    rc = fc_deleteSequence(sequences[i]);
    fail_unless(rc == FC_SUCCESS, "failed to delete sequence");
  }

  // Create some more (fills in the cracks)
  for (i = 1; i < num; i+=2) {
    rc = fc_createSequence(dataset, "fake", &sequences[i]);
    fail_unless(rc == FC_SUCCESS, "failed to create sequence");
  }

  // Check slotID and uID
  for (i = 0; i < num; i++) {
    fail_unless(sequences[i].slotID == i, "mismatch of slot id");
    fail_unless(sequences[i].uID == uID_start + uID_incs[i],
		"mismatch of uID");
  }

  // cleanup
  rc = fc_deleteDataset(dataset);
  fail_unless(rc == FC_SUCCESS, "failed to delete dataset at end of test");
}
END_TEST

// test creating and deleting sequences, also test querying dataset
// for sequences before, inbetween and after.
START_TEST(create_get_delete)
{
  FC_ReturnCode rc;
  int i;
  FC_Dataset dataset;
  int numSequence = 10, temp_numSequence;
  FC_Sequence sequences[10], *temp_sequences, temp_sequence;
  FC_Sequence badSequence = { 999, 999 };
  char names[10][100] = { "one", "two", "three", "four", "five", "six",
			  "seven", "eight", "nine", "ten" };
  char newName[100] = { "sparkly new" };

  // make a dataset to play in
  rc = fc_createDataset("temp dataset", &dataset);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create dataset");

  // create some sequences
  for (i = 0; i < numSequence; i++) {
    rc = fc_createSequence(dataset, names[i], &sequences[i]);
    fail_unless(rc == FC_SUCCESS, "failed to create sequence");
    fail_unless(!FC_HANDLE_EQUIV(sequences[i], FC_NULL_SEQUENCE),
		"created sequence should not = FC_NULL_SEQUENCE");
  }

  // get sequences (should be numSequence);
  rc = fc_getSequences(dataset, &temp_numSequence, &temp_sequences);
  fail_unless(rc == FC_SUCCESS, "failed to get sequences");
  fail_unless(temp_numSequence == numSequence, "mismatch of numSequence");
  for (i = 0; i < numSequence; i++) 
    fail_unless(FC_HANDLE_EQUIV(temp_sequences[i], sequences[i]), 
		"mismatch of sequence handles");
  free(temp_sequences);
  rc = fc_getNumSequence(dataset, &temp_numSequence);
  fail_unless(rc == FC_SUCCESS, "failed to get numSequence");
  fail_unless(temp_numSequence == numSequence, "mismatch of numSequence");
  for (i = 0; i < numSequence; i++) {
    rc = fc_getSequenceByName(dataset, names[i], &temp_numSequence,&temp_sequences);
    fail_unless(rc == FC_SUCCESS, "failed to get sequence by name");
    fail_unless(temp_numSequence == 1, "wrong number of matching sequences");
    fail_unless(FC_HANDLE_EQUIV(temp_sequences[0], sequences[i]),
		"mismatch of sequence handles");
    free(temp_sequences);
  }

  //temporary sequence for testing multiples
  rc = fc_createSequence(dataset, names[0], &temp_sequence);
  fail_unless(rc == FC_SUCCESS, "failed to create sequence");
  fail_unless(!FC_HANDLE_EQUIV(temp_sequence, FC_NULL_SEQUENCE),
	      "created sequence should not = FC_NULL_SEQUENCE");
  rc = fc_getSequenceByName(dataset, names[0], &temp_numSequence,&temp_sequences);
  fail_unless(rc == FC_SUCCESS, "failed to get sequence by name");
  fail_unless(temp_numSequence == 2, "wrong number of matching sequences");
  fail_unless(((FC_HANDLE_EQUIV(temp_sequences[0], sequences[0]) &&
	       FC_HANDLE_EQUIV(temp_sequences[1], temp_sequence)) ||
	      (FC_HANDLE_EQUIV(temp_sequences[1], sequences[0]) &&
	       FC_HANDLE_EQUIV(temp_sequences[0], temp_sequence))),
	      "mismatch of sequences handles");
  fc_deleteSequence(temp_sequence);
  free(temp_sequences);


  // change the name of the first sequence
  rc = fc_changeSequenceName(sequences[0], newName);
  fail_unless(rc == FC_SUCCESS, "failed to change sequence name");
  rc = fc_getSequenceByName(dataset, names[0], &temp_numSequence,&temp_sequences);
  fail_unless(rc == FC_SUCCESS, "failed to get sequence by name");
  fail_unless(temp_numSequence == 0, "should be no matching sequences for old name");
  fail_unless(temp_sequences == NULL, "should return null for no matching sequences");
  rc = fc_getSequenceByName(dataset, newName, &temp_numSequence,&temp_sequences);
  fail_unless(rc == FC_SUCCESS, "new name should work for find");
  fail_unless(temp_numSequence == 1, "wrong number of matching sequences");
  fail_unless(FC_HANDLE_EQUIV(temp_sequences[0], sequences[0]),
	      "mismatch of sequence handle");
  free(temp_sequences);

  // delete half of the sequences (alternate)
  for (i = 0; i < numSequence; i+=2) {
    rc = fc_deleteSequence(sequences[i]);
    fail_unless(rc == FC_SUCCESS, "failed to delete sequence");
  }

  // get sequences (should be numSequence/2);
  fc_getSequences(dataset, &temp_numSequence, &temp_sequences);
  fail_unless(rc == FC_SUCCESS, "failed to get sequences");
  fail_unless(temp_numSequence == numSequence/2, "mismatch of numSequence");
  for (i = 0; i < numSequence/2; i++) 
    fail_unless(FC_HANDLE_EQUIV(temp_sequences[i], sequences[i*2+1]), 
		"mismatch of sequence handles");
  free(temp_sequences);
  for (i = 0; i < numSequence/2; i++) {
    rc = fc_getSequenceByName(dataset, names[i*2+1], &temp_numSequence,&temp_sequences);
    fail_unless(rc == FC_SUCCESS, "failed to get sequence by name");
    fail_unless(temp_numSequence ==1, "wrong number of matching sequences");
    fail_unless(FC_HANDLE_EQUIV(temp_sequences[0], sequences[i*2+1]),
		"mismatch of sequence handles");
    free(temp_sequences);
  }

  // delete remaining sequences
  for (i = 1; i < numSequence; i+=2) { 
    rc = fc_deleteSequence(sequences[i]);
    fail_unless(rc == FC_SUCCESS, "failed to delete sequence");
  }
  
  // get sequences (should be none)
  rc = fc_getSequences(dataset, &temp_numSequence, &temp_sequences);
  fail_unless(rc == FC_SUCCESS, "failed to get sequences from empty library");
  fail_unless(temp_numSequence == 0 && temp_sequences == NULL,
	      "should return 0 if all sequences deleted");

  // create one more sequence for further testing
  rc = fc_createSequence(dataset, names[0], &sequences[0]);
  fail_unless(rc == FC_SUCCESS, "aborted: failed to create sequence sequence");

  // ---- test special cases

  // not an error to delete FC_NULL_SEQUENCE
  rc = fc_deleteSequence(FC_NULL_SEQUENCE);
  fail_unless(rc == FC_SUCCESS, "should not error to delete NULL sequence");

  // ---- test error conditions

  // bad args to fc_createSequence()
  temp_sequence = badSequence;
  rc = fc_createSequence(FC_NULL_DATASET, names[0], &temp_sequence);
  fail_unless(rc != FC_SUCCESS, 
	      "should have failed to create sequence with null database");
  fail_unless(FC_HANDLE_EQUIV(temp_sequence, FC_NULL_SEQUENCE),
	      "fail should return NULL sequence");
  temp_sequence = badSequence;
  rc = fc_createSequence(dataset, NULL, &temp_sequence);
  fail_unless(rc != FC_SUCCESS, 
	      "should have failed to create sequence with null name");
  fail_unless(FC_HANDLE_EQUIV(temp_sequence, FC_NULL_SEQUENCE),
	      "fail should return NULL sequence");
  rc = fc_createSequence(dataset, names[0], NULL);
  fail_unless(rc != FC_SUCCESS, 
	      "should have failed to create sequence with null handle");

  // bad args to fc_changeSequenceName()
  rc = fc_changeSequenceName(sequences[0], NULL);
  fail_unless(rc != FC_SUCCESS, "should fail if new name is NULL");
  rc = fc_getSequenceByName(dataset, names[0], &temp_numSequence,&temp_sequences);
  fail_unless(rc == FC_SUCCESS && temp_numSequence == 1 &&
	      FC_HANDLE_EQUIV(temp_sequences[0],sequences[0]),
	      "should not change name of sequence");
  free(temp_sequences);

  // bad args to fc_getSequences()
  temp_numSequence = 99;
  temp_sequences = (FC_Sequence*)1;
  rc = fc_getSequences(FC_NULL_DATASET, &temp_numSequence, &temp_sequences);
  fail_unless(rc != FC_SUCCESS, 
	      "should have failed to create sequence with null dataset");
  fail_unless(temp_numSequence == -1,
	      "fail should return -1 for numSequence"); 
  fail_unless(temp_sequences == NULL, "fail should return null");
  temp_sequences = (FC_Sequence*)1;
  rc = fc_getSequences(dataset, NULL, &temp_sequences);
  fail_unless(rc != FC_SUCCESS, 
	      "should have failed to get sequence with NULL numSequence");
  fail_unless(temp_sequences == NULL, "fail should return null");
  temp_numSequence = 999;
  rc = fc_getSequences(dataset, &temp_numSequence, NULL);
  fail_unless(rc != FC_SUCCESS, "should error to get just numSequence");
  fail_unless(temp_numSequence == -1, "fail should return null");

   
  // bad args to fc_getNumSequence()
  temp_numSequence = 99;
  rc = fc_getNumSequence(FC_NULL_DATASET, &temp_numSequence);
  fail_unless(rc != FC_SUCCESS, 
	      "should have failed to get numSequence with null dataset");
  fail_unless(temp_numSequence == -1,
	      "fail should return -1 for numSequence"); 
  rc = fc_getNumSequence(dataset, NULL);
  fail_unless(rc != FC_SUCCESS, 
	      "should have failed to get numSsequence with NULL numSequence");
  fail_unless(temp_sequences == NULL, "mismatch of sequences");
 
  // bad args to fc_getSequenceByName()
  rc = fc_getSequenceByName(FC_NULL_DATASET, names[0], &temp_numSequence, 
			    &temp_sequences);
  fail_unless(rc != FC_SUCCESS, 
	      "should have failed to get sequence with null dataset");
  fail_unless(temp_numSequence == -1, "fail should return -1");
  fail_unless(temp_sequences == NULL, "fail should return NULL");
  rc = fc_getSequenceByName(dataset, NULL, &temp_numSequence, 
			    &temp_sequences);
  fail_unless(rc != FC_SUCCESS, 
	      "should have failed to get sequence with null dataset");
  fail_unless(temp_numSequence == -1, "fail should return -1");
  fail_unless(temp_sequences == NULL, "fail should return NULL");
  rc = fc_getSequenceByName(dataset, names[0], NULL,
			    &temp_sequences);
  fail_unless(rc != FC_SUCCESS, 
	      "should have failed to get sequence with null arg");
  fail_unless(temp_sequences == NULL, "fail should return NULL");
  rc = fc_getSequenceByName(dataset, names[0], &temp_numSequence, NULL);
  fail_unless(rc != FC_SUCCESS, 
	      "should have failed to get sequence with null arg");
  fail_unless(temp_numSequence == -1, "fail should return -1");
  
  // bad args to fc_deleteSequence()
  rc = fc_deleteSequence(badSequence);
  fail_unless(rc != FC_SUCCESS, 
	      "should error to delete nonexistent sequence");

  // --- done

  // delete last sequence
  rc = fc_deleteSequence(sequences[0]);
  fail_unless(rc == FC_SUCCESS, "final delete sequence failed");

  // delete the dataset
  rc = fc_deleteDataset(dataset);
  fail_unless(rc == FC_SUCCESS, "failed to delete dataset at end of test");

}
END_TEST


// query meta data
START_TEST(metadata_query)
{
  FC_ReturnCode rc;
  char name[20] = "blue berry", *temp_name;
  FC_Dataset dataset, temp_dataset, badDataset = { 999, 999 };
  FC_Sequence sequence, badSequence = { 999, 999 };

  // make a dataset to play in
  rc = fc_createDataset("temp dataset", &dataset);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create dataset");

  // create a sequence to play with
  rc = fc_createSequence(dataset, name, &sequence);
  fail_unless(rc == FC_SUCCESS, "failed to created sequence");
 
  // test fc_isSequenceValid()
  fail_unless(fc_isSequenceValid(sequence), 
	      "failed to validate a valid sequence");

  // test fc_getSequenceName()
  temp_name = NULL;
  rc = fc_getSequenceName(sequence, &temp_name);
  fail_unless(rc == FC_SUCCESS, "failed to get sequence name");
  fail_unless(!strcmp(name, temp_name), "mismatch of name");
  // should have returned a copy, mess with temp_name & try again to make sure
  temp_name[0] = 'q';
  free(temp_name);
  rc = fc_getSequenceName(sequence, &temp_name);
  fail_unless(rc == FC_SUCCESS, "failed to get sequence name");
  fail_unless(!strcmp(name, temp_name), "mismatch of name");
  free(temp_name);

  // test fc_getDatasetFromSequence()
  temp_dataset = badDataset;
  rc = fc_getDatasetFromSequence(sequence, &temp_dataset);
  fail_unless(rc == FC_SUCCESS, "failed to get parent dataset");
  fail_unless(FC_HANDLE_EQUIV(temp_dataset, dataset),
	      "mismatch of parent dataset");

  // --- check with bad args

  // fc_isSequenceValid()
  fail_unless(!fc_isSequenceValid(badSequence), 
	      "badSequence should NOT be valid");
  fail_unless(!fc_isSequenceValid(FC_NULL_SEQUENCE), 
	      "FC_NULL_SEQUENCE should not be valid");

  // fc_getSequenceName()
  temp_name = (char*)1;
  rc = fc_getSequenceName(badSequence, &temp_name);
  fail_unless(rc != FC_SUCCESS, "badSequence should NOT return a name");
  fail_unless(temp_name == NULL, "fail should return NULL name");
  rc = fc_getSequenceName(sequence, NULL);
  fail_unless(rc != FC_SUCCESS, "should fail for NULL name");

  // fc_getDatasetFromSequence()
  temp_dataset = badDataset;
  rc = fc_getDatasetFromSequence(badSequence, &temp_dataset);
  fail_unless(rc != FC_SUCCESS, "badSequence should NOT return a dataset");
  fail_unless(FC_HANDLE_EQUIV(temp_dataset, FC_NULL_DATASET),
	      "failure should return NULL dataset");
  rc = fc_getDatasetFromSequence(sequence, NULL);
  fail_unless(rc != FC_SUCCESS, "should fail for NULL dataset");

  // --- all done

  // cleanup & one last test
  fc_deleteSequence(sequence);
  fail_unless(!fc_isSequenceValid(sequence), 
	      "handle should not be valid after a delete"); 

  // delete the dataset
  rc = fc_deleteDataset(dataset);
  fail_unless(rc == FC_SUCCESS, "failed to delete dataset at end of test");
}
END_TEST

// query the sequence coords
START_TEST(coords_query)
{
  FC_ReturnCode rc;
  int i;
  char name[20] = "blue berry", name2[20] = "banana nut";
  FC_Dataset dataset;
  FC_Sequence sequence, sequence2, empty_sequence, badSequence = { 999, 999 };
  _FC_SeqSlot *seqSlot;
  int numDataType = 4, numStep = 10, temp_numStep;
  FC_DataType dataTypes[4] = { FC_DT_CHAR, FC_DT_INT, 
			       FC_DT_FLOAT, FC_DT_DOUBLE };
  FC_DataType temp_dataType;
  char charCoords[10] = { 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j' };
  int intCoords[10] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
  float floatCoords[10];
  double doubleCoords[10];
  int sizes[4] = {sizeof(char), sizeof(int), sizeof(float), sizeof(double) };
  void *coords[4] = { charCoords, intCoords, floatCoords, doubleCoords };
  void *temp_coords, *coords_cpy;
  
  // make a dataset to play in
  rc = fc_createDataset("temp dataset", &dataset);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create dataset");

  // setup coords
  for (i = 0; i < numStep; i++) {
    floatCoords[i] = intCoords[i] + 0.1;
    doubleCoords[i] = intCoords[i] + 0.00000000001;
  }

  for (i = 0; i < numDataType; i++) {
    // create three sequences
    rc = fc_createSequence(dataset, name, &sequence);
    fail_unless(rc == FC_SUCCESS, "abort: failed to created sequence");
    rc = fc_createSequence(dataset, name2, &sequence2);
    fail_unless(rc == FC_SUCCESS, "abort: failed to created sequence");
    rc = fc_createSequence(dataset, "empty sequence", &empty_sequence);

    // check default values
    rc = fc_getSequenceInfo(sequence, &temp_numStep, &temp_dataType);
    fail_unless(rc == FC_SUCCESS, "failed to get sequence info");
    fail_unless(temp_numStep == 0 && temp_dataType == FC_DT_UNKNOWN,
		"empty sequence should have empty values");

    // sequence1: set/get coords by copy on and check
    temp_numStep = 999;
    temp_dataType = FC_DT_UNKNOWN;
    temp_coords = NULL;
    rc = fc_setSequenceCoords(sequence, numStep, dataTypes[i], coords[i]);
    fail_unless(rc == FC_SUCCESS, "failed to set coords");
    rc = fc_getSequenceInfo(sequence, &temp_numStep, &temp_dataType);
    fail_unless(rc == FC_SUCCESS, "failed to get info");
    fail_unless(temp_numStep == numStep, "mismatch of numStep");
    fail_unless(temp_dataType == dataTypes[i], "mismatch of datatype");
    rc = fc_getSequenceCoordsPtr(sequence, &temp_coords);
    fail_unless(rc == FC_SUCCESS, "failed to get coords");
    fail_unless(!memcmp(temp_coords, coords[i], numStep*sizes[i]),
		"mismatch of coords");

    // sequence2: set/get coords pointer and check
    coords_cpy = malloc(numStep*sizes[i]);
    memcpy(coords_cpy, coords[i], numStep*sizes[i]);
    temp_numStep = 999;
    temp_dataType = FC_DT_UNKNOWN;
    temp_coords = NULL;
    rc = fc_setSequenceCoordsPtr(sequence2, numStep, dataTypes[i],
				 coords_cpy); 
    fail_unless(rc == FC_SUCCESS, "failed to set coords");
    rc = fc_getSequenceInfo(sequence2, &temp_numStep, &temp_dataType);
    fail_unless(rc == FC_SUCCESS, "failed to get info");
    fail_unless(temp_numStep == numStep, "mismatch of numStep");
    fail_unless(temp_dataType == dataTypes[i], "mismatch of datatype");
    rc = fc_getSequenceCoordsPtr(sequence2, &temp_coords);
    fail_unless(rc == FC_SUCCESS, "failed to get coords");
    fail_unless(!memcmp(temp_coords, coords[i], numStep*sizes[i]),
		"mismatch of coords");

    // test set by copy vs. set the pointer
    seqSlot = _fc_getSeqSlot(sequence);
    fail_unless(seqSlot->coords != coords[i], "should have copied the coords");
    seqSlot = _fc_getSeqSlot(sequence2);
    fail_unless(seqSlot->coords == coords_cpy, "should have ptr to coords");

    // check other query functions that are simpler forms of fc_getSequenceInfo
    rc = fc_getSequenceNumStep(sequence, &temp_numStep);
    fail_unless(rc == FC_SUCCESS, "failed to get numStep");
    fail_unless(temp_numStep == numStep, "mismatch of numStep");
    rc = fc_getSequenceDataType(sequence, &temp_dataType);
    fail_unless(rc == FC_SUCCESS, "failed to get data type");
    fail_unless(temp_dataType == dataTypes[i], "mismatch of datatype");
 
    // ---- check special cases

    // fc_getSequenceInfo -- most arguments are optional
    rc = fc_getSequenceInfo(sequence, NULL, &temp_dataType);
    fail_unless(rc == FC_SUCCESS, "failed to get info with null numStep");
    fail_unless(temp_dataType == dataTypes[i], "mismatch of datatype");
    rc = fc_getSequenceInfo(sequence, &temp_numStep, NULL);
    fail_unless(rc == FC_SUCCESS, "failed to get info with null datatype");
    fail_unless(temp_numStep == numStep, "mismatch of numStep");

    // ---- check errors

    // fc_getSequenceInfo() -- bad args
    rc = fc_getSequenceInfo(badSequence, &temp_numStep, &temp_dataType);
    fail_unless(rc != FC_SUCCESS, "should fail with bad sequence");
    fail_unless(temp_numStep == -1 && temp_dataType == FC_DT_UNKNOWN,
		"fail should return NULLs");
    rc = fc_getSequenceInfo(sequence, NULL, NULL);
    fail_unless(rc != FC_SUCCESS, "can't have all args be NULL");

    // functions similar to fc_getSequenceInfo() -- bad args
    rc = fc_getSequenceNumStep(badSequence, &temp_numStep);
    fail_unless(rc != FC_SUCCESS, "should fail with bad sequence");
    fail_unless(temp_numStep == -1, "fail should return NULL");
    rc = fc_getSequenceNumStep(sequence, NULL);
    fail_unless(rc != FC_SUCCESS, "should fail with NULL numStep");
    rc = fc_getSequenceDataType(badSequence, &temp_dataType);
    fail_unless(rc != FC_SUCCESS, "should fail with bad sequence");
    fail_unless(temp_dataType == FC_DT_UNKNOWN, "fail should return NULL");
    rc = fc_getSequenceDataType(sequence, NULL);
    fail_unless(rc != FC_SUCCESS, "should fail NULL data type");

    // fc_getSequenceCoordsPtr() -- bad args
    rc = fc_getSequenceCoordsPtr(badSequence, &temp_coords);
    fail_unless(rc != FC_SUCCESS, "should fail with bad sequence");
    fail_unless(temp_coords == NULL, "fail should return nulls");
    rc = fc_getSequenceCoordsPtr(sequence, NULL);
    fail_unless(rc != FC_SUCCESS, "should fail if coords = NULL");
    
    // fc_setSequenceCoords() -- error to try to reset
    rc = fc_setSequenceCoords(sequence, numStep, dataTypes[i], coords[i]);
    fail_unless(rc != FC_SUCCESS, "shouldn't be able to reset coords");

    // fc_setSequenceCoords() -- bad args
    rc = fc_setSequenceCoords(badSequence, numStep, dataTypes[i], coords[i]);
    fail_unless(rc != FC_SUCCESS, "should fail with bad sequence");
    rc = fc_setSequenceCoords(empty_sequence, 0, dataTypes[i], coords[i]);
    fail_unless(rc != FC_SUCCESS, "should fail with numStep < 1");
    rc = fc_setSequenceCoords(empty_sequence, numStep, FC_DT_UNKNOWN, coords[i]);
    fail_unless(rc != FC_SUCCESS, "should fail with unknown datatype");
    rc = fc_setSequenceCoords(empty_sequence, numStep, -999, coords[i]);
    fail_unless(rc != FC_SUCCESS, "should fail with bad datatype");
    rc = fc_setSequenceCoords(empty_sequence, numStep, dataTypes[i], NULL);
    fail_unless(rc != FC_SUCCESS, "should fail with no coords");
    rc = fc_getSequenceInfo(empty_sequence, &temp_numStep, &temp_dataType);
    fail_unless(rc == FC_SUCCESS, "failed to get sequence info");
    fail_unless(temp_numStep == 0 && temp_dataType == FC_DT_UNKNOWN,
		"empty sequence should have empty values");

    // fc_setSequenceCoordsPtr() -- error to try to reset
    rc = fc_setSequenceCoordsPtr(sequence, numStep, dataTypes[i], coords[i]);
    fail_unless(rc != FC_SUCCESS, "shouldn't be able to reset coords");

    // fc_setSequenceCoordsPtr() -- bad args
    rc = fc_setSequenceCoordsPtr(badSequence, numStep, dataTypes[i], coords[i]);
    fail_unless(rc != FC_SUCCESS, "should fail with bad sequence");
    rc = fc_setSequenceCoordsPtr(empty_sequence, 0, dataTypes[i], coords[i]);
    fail_unless(rc != FC_SUCCESS, "should fail with numStep < 1");
    rc = fc_setSequenceCoordsPtr(empty_sequence, numStep, FC_DT_UNKNOWN, coords[i]);
    fail_unless(rc != FC_SUCCESS, "should fail with unknown datatype");
    rc = fc_setSequenceCoordsPtr(empty_sequence, numStep, -999, coords[i]);
    fail_unless(rc != FC_SUCCESS, "should fail with bad datatype");
    rc = fc_setSequenceCoordsPtr(empty_sequence, numStep, dataTypes[i], NULL);
    fail_unless(rc != FC_SUCCESS, "should fail with no coords");
    rc = fc_getSequenceInfo(empty_sequence, &temp_numStep, &temp_dataType);
    fail_unless(rc == FC_SUCCESS, "failed to get sequence info");
    fail_unless(temp_numStep == 0 && temp_dataType == FC_DT_UNKNOWN,
		"empty sequence should have empty values");

    // empty sequence should still be empty (all set's failed)
    rc = fc_getSequenceInfo(empty_sequence, &temp_numStep, &temp_dataType);
    fail_unless(rc == FC_SUCCESS, "failed to get sequence info");
    fail_unless(temp_numStep == 0 && temp_dataType == FC_DT_UNKNOWN,
		"empty sequence should have empty values");

    // --- all done
    // cleanup
    rc = fc_deleteSequence(sequence);
    fail_unless(rc == FC_SUCCESS, 
		"failed to delete sequence at end of tests");
    rc = fc_deleteSequence(sequence2);
    fail_unless(rc == FC_SUCCESS, 
		"failed to delete sequence at end of tests");
    rc = fc_deleteSequence(empty_sequence);
    fail_unless(rc == FC_SUCCESS, 
		"failed to delete sequence at end of tests");
  }

  // delete the dataset
  rc = fc_deleteDataset(dataset);
  fail_unless(rc == FC_SUCCESS, "failed to delete dataset at end of test");
}
END_TEST

// query the sequence coords
START_TEST(get_coords_as)
{
  FC_ReturnCode rc;
  int i, j;
  char names[4][20] = { "blue berry", "pineapple", "guava", "coconut" };
  FC_Dataset dataset;
  FC_Sequence sequences[4], badSequence = { 999, 999 };
  int numDataType = 4, numStep = 10;
  FC_DataType dataTypes[4] = { FC_DT_CHAR, FC_DT_INT, 
			       FC_DT_FLOAT, FC_DT_DOUBLE };
  char charCoords[10] = { 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j' };
  int intCoords[10] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
  float floatCoords[10];
  double doubleCoords[10];
  void *coords[4] = { charCoords, intCoords, floatCoords, doubleCoords };
  void *temp_coords;
  
  // make a dataset to play in
  rc = fc_createDataset("temp dataset", &dataset);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create dataset");

  // setup flt & double coords
  for (i = 0; i < numStep; i++) {
    floatCoords[i] = intCoords[i] + 0.1;
    doubleCoords[i] = intCoords[i] + 0.00000000001;
  }

  // create a sequence for each dataType
  for (i = 0; i < numDataType; i++) {
    rc = fc_createSequence(dataset, names[i], &sequences[i]);
    fail_unless(rc == FC_SUCCESS, "abort: failed to created sequence");
    rc = fc_setSequenceCoords(sequences[i], numStep, dataTypes[i], coords[i]);
    fail_unless(rc == FC_SUCCESS, "abort: failed to set coords");
  }

  // test getting char sequence as each type
  rc = fc_getSequenceCoordsAsDataType(sequences[0], FC_DT_CHAR, &temp_coords);
  fail_unless(rc == FC_SUCCESS, "failed to get coords as char");
  for (j = 0; j < numStep; j++) 
    fail_unless(((char*)temp_coords)[j] == charCoords[j],
		"mismatch of char to char");
  free(temp_coords);
  rc = fc_getSequenceCoordsAsDataType(sequences[0], FC_DT_INT, &temp_coords);
  fail_unless(rc == FC_SUCCESS, "failed to get coords as int");
  for (j = 0; j < numStep; j++) 
    fail_unless(((int*)temp_coords)[j] == (int)charCoords[j],
		"mismatch of char to int");
  free(temp_coords);
  rc = fc_getSequenceCoordsAsDataType(sequences[0], FC_DT_FLOAT, &temp_coords);
  fail_unless(rc == FC_SUCCESS, "failed to get coords as float");
  for (j = 0; j < numStep; j++) 
    fail_unless(((float*)temp_coords)[j] == (float)charCoords[j],
		"mismatch of char to float");
  free(temp_coords);
  rc = fc_getSequenceCoordsAsDataType(sequences[0], FC_DT_DOUBLE, 
				      &temp_coords);
  fail_unless(rc == FC_SUCCESS, "failed to get coords as double");
  for (j = 0; j < numStep; j++) 
    fail_unless(((double*)temp_coords)[j] == (double)charCoords[j],
		"mismatch of char to double");
  free(temp_coords);

  // test getting int sequence as each type
  rc = fc_getSequenceCoordsAsDataType(sequences[1], FC_DT_CHAR, &temp_coords);
  fail_unless(rc == FC_SUCCESS, "failed to get coords as char");
  for (j = 0; j < numStep; j++) 
    fail_unless(((char*)temp_coords)[j] == (char)intCoords[j],
		"mismatch of int to char");
  free(temp_coords);
  rc = fc_getSequenceCoordsAsDataType(sequences[1], FC_DT_INT, &temp_coords);
  fail_unless(rc == FC_SUCCESS, "failed to get coords as int");
  for (j = 0; j < numStep; j++) 
    fail_unless(((int*)temp_coords)[j] == intCoords[j],
		"mismatch of int to int");
  free(temp_coords);
  rc = fc_getSequenceCoordsAsDataType(sequences[1], FC_DT_FLOAT, &temp_coords);
  fail_unless(rc == FC_SUCCESS, "failed to get coords as float");
  for (j = 0; j < numStep; j++) 
    fail_unless(((float*)temp_coords)[j] == (float)intCoords[j],
		"mismatch of int to float");
  free(temp_coords);
  rc = fc_getSequenceCoordsAsDataType(sequences[1], FC_DT_DOUBLE, 
				      &temp_coords);
  fail_unless(rc == FC_SUCCESS, "failed to get coords as double");
  for (j = 0; j < numStep; j++) 
    fail_unless(((double*)temp_coords)[j] == (double)intCoords[j],
		"mismatch of int to double");
  free(temp_coords);

  // test getting float sequence as each type
  rc = fc_getSequenceCoordsAsDataType(sequences[2], FC_DT_CHAR, &temp_coords);
  fail_unless(rc == FC_SUCCESS, "failed to get coords as char");
  for (j = 0; j < numStep; j++) 
    fail_unless(((char*)temp_coords)[j] == (char)floatCoords[j],
		"mismatch of float to char");
  free(temp_coords);
  rc = fc_getSequenceCoordsAsDataType(sequences[2], FC_DT_INT, &temp_coords);
  fail_unless(rc == FC_SUCCESS, "failed to get coords as int");
  for (j = 0; j < numStep; j++) 
    fail_unless(((int*)temp_coords)[j] == (int)floatCoords[j],
		"mismatch of float to int");
  free(temp_coords);
  rc = fc_getSequenceCoordsAsDataType(sequences[2], FC_DT_FLOAT, &temp_coords);
  fail_unless(rc == FC_SUCCESS, "failed to get coords as float");
  for (j = 0; j < numStep; j++) 
    fail_unless(((float*)temp_coords)[j] == floatCoords[j],
		"mismatch of float to float");
  free(temp_coords);
  rc = fc_getSequenceCoordsAsDataType(sequences[2], FC_DT_DOUBLE, 
				      &temp_coords);
  fail_unless(rc == FC_SUCCESS, "failed to get coords as double");
  for (j = 0; j < numStep; j++) 
    fail_unless(((double*)temp_coords)[j] == (double)floatCoords[j],
		"mismatch of float to double");
  free(temp_coords);

  // test getting double sequence as each type
  rc = fc_getSequenceCoordsAsDataType(sequences[3], FC_DT_CHAR, &temp_coords);
  fail_unless(rc == FC_SUCCESS, "failed to get coords as char");
  for (j = 0; j < numStep; j++) 
    fail_unless(((char*)temp_coords)[j] == (char)doubleCoords[j],
		"mismatch of douible to char");
  free(temp_coords);
  rc = fc_getSequenceCoordsAsDataType(sequences[3], FC_DT_INT, &temp_coords);
  fail_unless(rc == FC_SUCCESS, "failed to get coords as int");
  for (j = 0; j < numStep; j++) 
    fail_unless(((int*)temp_coords)[j] == (int)doubleCoords[j],
		"mismatch of double to int");
  free(temp_coords);
  rc = fc_getSequenceCoordsAsDataType(sequences[3], FC_DT_FLOAT, &temp_coords);
  fail_unless(rc == FC_SUCCESS, "failed to get coords as float");
  for (j = 0; j < numStep; j++) 
    fail_unless(((float*)temp_coords)[j] == (float)doubleCoords[j],
		"mismatch of double to float");
  free(temp_coords);
  rc = fc_getSequenceCoordsAsDataType(sequences[3], FC_DT_DOUBLE, 
				      &temp_coords);
  fail_unless(rc == FC_SUCCESS, "failed to get coords as double");
  for (j = 0; j < numStep; j++) 
    fail_unless(((double*)temp_coords)[j] == doubleCoords[j],
		"mismatch of double to double");
  free(temp_coords);

  // ---- check errors
  
  // fc_getSequenceCoorsAsDataType() -- bad args
  rc = fc_getSequenceCoordsAsDataType(badSequence, FC_DT_CHAR, &temp_coords);
  fail_unless(rc != FC_SUCCESS, "should fail with bad sequence");
  fail_unless(temp_coords == NULL, "fail should return NULL");
  rc = fc_getSequenceCoordsAsDataType(sequences[0], FC_DT_UNKNOWN, &temp_coords);
  fail_unless(rc != FC_SUCCESS, "should fail with FC_DT_UNKNOWN");
  fail_unless(temp_coords == NULL, "fail should return NULL");
  rc = fc_getSequenceCoordsAsDataType(sequences[0], -99, &temp_coords);
  fail_unless(rc != FC_SUCCESS, "should fail with bad datatype");
  fail_unless(temp_coords == NULL, "fail should return NULL");

  // --- all done

  // delete the dataset
  rc = fc_deleteDataset(dataset);
  fail_unless(rc == FC_SUCCESS, "failed to delete dataset at end of test");
}
END_TEST

// copy and shiftandscale
START_TEST(copy_test)
{
  FC_ReturnCode rc;
  int i,j,k,l;
  char name[20] = "blue berry", copy_name[20] = "banana nut", *temp_name;
  char ss_name[20] = "ss_name";
  FC_Dataset dataset, dataset2;
  int numSequence;
  FC_Sequence sequence, copy_sequence, badSequence = { 999, 999 };
  FC_Sequence ss_sequence;
  int numDataType = 4, numStep = 10, temp_numStep;
  FC_DataType dataTypes[4] = { FC_DT_CHAR, FC_DT_INT, 
			       FC_DT_FLOAT, FC_DT_DOUBLE };
  FC_DataType temp_dataType;
  char charCoords[10] = { 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j' };
  int intCoords[10] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
  float floatCoords[10];
  double doubleCoords[10];
  int sizes[4] = {sizeof(char), sizeof(int), sizeof(float), sizeof(double) };
  void *coords[4] = { charCoords, intCoords, floatCoords, doubleCoords };
  void *temp_coords, *ss_coords;

  int numShiftVals = 3;
  int numScaleVals = 3;
  double shiftVals[3] = {0.0,1.0,-1.5};
  double scaleVals[3] = {1.0,0.4,-0.4}; //scale == 0 is tested in bad args
  double shift,scale;

  
  // setup coords
  for (i = 0; i < numStep; i++) {
    floatCoords[i] = intCoords[i] + 0.1;
    doubleCoords[i] = intCoords[i] + 0.00000000001;
  }

  // make a dataset to play in
  rc = fc_createDataset("temp dataset", &dataset);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create dataset");

  // make a second dataset
  rc = fc_createDataset("temp dataset2", &dataset2);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create dataset");

  for (i = 0; i < numDataType; i++) {
    // create a sequence & set coords
    rc = fc_createSequence(dataset, name, &sequence);
    fail_unless(rc == FC_SUCCESS, "abort: failed to created sequence");
    rc = fc_setSequenceCoords(sequence, numStep, dataTypes[i], coords[i]);
    fail_unless(rc == FC_SUCCESS, "failed to set coords");

    // fc_copySequence() to same dataset & check it
    rc = fc_copySequence(sequence, dataset, copy_name, &copy_sequence);
    fail_unless(rc == FC_SUCCESS, "failed to copy to same dataset");
    rc = fc_getSequenceName(copy_sequence, &temp_name);
    fail_unless(rc == FC_SUCCESS, "failed to get name of copy");
    fail_unless(!strcmp(temp_name, copy_name), "mismatch of name");
    free(temp_name);
    rc = fc_getSequenceInfo(copy_sequence, &temp_numStep, &temp_dataType);
    fail_unless(rc == FC_SUCCESS, "failed to get info");
    fail_unless(temp_numStep == numStep, "mismatch of numStep");
    fail_unless(temp_dataType == dataTypes[i], "mismatch of datatype");
    rc = fc_getSequenceCoordsPtr(sequence, &temp_coords);
    fail_unless(rc == FC_SUCCESS, "failed to get coords");
    fail_unless(!memcmp(temp_coords, coords[i], numStep*sizes[i]),
		"mismatch of coords");

    //fc_shiftAndScaleSequence to same dataset and check it
    for (j = 0; j < numShiftVals; j++){
      for (k = 0; k < numScaleVals; k++){
	shift = shiftVals[j];
	scale = scaleVals[k];

	rc = fc_shiftAndScaleSequence(sequence, dataset, shift,scale,ss_name, &ss_sequence);
	if (dataTypes[i] == FC_DT_CHAR || FC_DBL_EQUIV(scale,0.0)){
	  fail_unless(rc != FC_SUCCESS, "should fail to shift and scale char seq or"
		      " scale == 0.0");
	  fail_unless(FC_HANDLE_EQUIV(ss_sequence, FC_NULL_SEQUENCE),
		      "fail should return NULL");
	} else {
	  fail_unless(rc == FC_SUCCESS, "failed to shiftand scale to same dataset");
	  rc = fc_getSequenceName(ss_sequence, &temp_name);
	  fail_unless(rc == FC_SUCCESS, "failed to get name of ss_sequence");
	  fail_unless(!strcmp(temp_name, ss_name), "mismatch of name");
	  free(temp_name);
	  rc = fc_getSequenceInfo(ss_sequence, &temp_numStep, &temp_dataType);
	  fail_unless(rc == FC_SUCCESS, "failed to get info");
	  fail_unless(temp_numStep == numStep, "mismatch of numStep");
	  fail_unless(temp_dataType == FC_DT_DOUBLE, "mismatch of datatype");
	  rc = fc_getSequenceCoordsPtr(ss_sequence, &temp_coords);
	  fail_unless(rc == FC_SUCCESS, "failed to get coords");
	  rc = fc_getSequenceCoordsPtr(sequence, &ss_coords);
	  fail_unless(rc == FC_SUCCESS, "failed to get coords");
	  for (l = 0; l < numStep; l++){
	    double orig;
	    switch(dataTypes[i]){
	    case FC_DT_INT:
	      orig = (double)(((int*)ss_coords)[l]);
	      break;
	    case FC_DT_FLOAT:
	      orig = (double)(((float*)ss_coords)[l]);
	      break;
	    case FC_DT_DOUBLE:
	      orig = ((double*)ss_coords)[l];
	      break;
	    default:
	      fail_unless(0==1,"bad data type - shouldnt reach this point");
	      break;
	    }
	    fail_unless(FC_DBL_EQUIV(((double*)temp_coords)[l],(orig+shift)*scale),
			"mismatch of datavals");
	  }
	  fc_deleteSequence(ss_sequence);
	}
      }
    }


    // fc_copySequence() to same different dataset & check it
    rc = fc_copySequence(sequence, dataset2, copy_name, &copy_sequence);
    fail_unless(rc == FC_SUCCESS, "failed to copy to same dataset");
    rc = fc_getSequenceName(copy_sequence, &temp_name);
    fail_unless(rc == FC_SUCCESS, "failed to get name of copy");
    fail_unless(!strcmp(temp_name, copy_name), "mismatch of name");
    free(temp_name);
    rc = fc_getSequenceInfo(copy_sequence, &temp_numStep, &temp_dataType);
    fail_unless(rc == FC_SUCCESS, "failed to get info");
    fail_unless(temp_numStep == numStep, "mismatch of numStep");
    fail_unless(temp_dataType == dataTypes[i], "mismatch of datatype");
    rc = fc_getSequenceCoordsPtr(sequence, &temp_coords);
    fail_unless(rc == FC_SUCCESS, "failed to get coords");
    fail_unless(!memcmp(temp_coords, coords[i], numStep*sizes[i]),
		"mismatch of coords");

    //fc_shiftAndScale() to different dataset, but im not checking it.
    rc = fc_shiftAndScaleSequence(sequence, dataset2, 1.0,1.0, ss_name, &ss_sequence);
    if (dataTypes[i] == FC_DT_CHAR){
      fail_unless(rc != FC_SUCCESS, "should fail to shiftandscale char seq to new dataset");
      fail_unless(FC_HANDLE_EQUIV(ss_sequence, FC_NULL_SEQUENCE),
		  "fail should return NULL");
    }else {
      fail_unless(rc == FC_SUCCESS, "failed to shiftandscale to new dataset");
      rc = fc_getSequenceName(ss_sequence, &temp_name);
      fail_unless(rc == FC_SUCCESS, "failed to get name of ss sequence");
      fail_unless(!strcmp(temp_name, ss_name), "mismatch of name");
      free(temp_name);
      rc = fc_getSequenceInfo(ss_sequence, &temp_numStep, &temp_dataType);
      fail_unless(rc == FC_SUCCESS, "failed to get info");
      fail_unless(temp_numStep == numStep, "mismatch of numStep");
      fail_unless(temp_dataType == FC_DT_DOUBLE, "mismatch of datatype");
      fc_deleteSequence(ss_sequence);
    }
 
    // --- special cases

    // copy with NULL name will use original's name
    rc = fc_copySequence(sequence, dataset, NULL, &copy_sequence);
    fail_unless(rc == FC_SUCCESS, "should fail to copy bad sequence");
    rc = fc_getSequenceName(copy_sequence, &temp_name);
    fail_unless(rc == FC_SUCCESS, "failed to get name of copy");
    fail_unless(!strcmp(temp_name, name), "mismatch of name");
    free(temp_name);
    rc = fc_deleteSequence(copy_sequence);
    fail_unless(rc == FC_SUCCESS, "failed to delete copied sequence");
  
    // ---- check errors

    // fc_copySequence() -- bad args
    rc = fc_copySequence(badSequence, dataset, copy_name, &copy_sequence);
    fail_unless(rc != FC_SUCCESS, "should fail to copy bad sequence");
    fail_unless(FC_HANDLE_EQUIV(copy_sequence, FC_NULL_SEQUENCE),
		"fail should return NULL");
    rc = fc_copySequence(sequence, FC_NULL_DATASET, copy_name, &copy_sequence);
    fail_unless(rc != FC_SUCCESS, "should fail to copy to bad dataset");
    fail_unless(FC_HANDLE_EQUIV(copy_sequence, FC_NULL_SEQUENCE),
		"fail should return NULL");
    rc = fc_copySequence(sequence, dataset, copy_name, NULL);
    fail_unless(rc != FC_SUCCESS, "should fail to copy if NULL handle");

    //fc_shiftAndScaleSequence() -- bad args
    rc = fc_shiftAndScaleSequence(badSequence, dataset, 1.0,1.0, ss_name, &ss_sequence);
    fail_unless(rc != FC_SUCCESS, "should fail to copy bad sequence");
    fail_unless(FC_HANDLE_EQUIV(ss_sequence, FC_NULL_SEQUENCE),
		"fail should return NULL");
    rc = fc_shiftAndScaleSequence(sequence, FC_NULL_DATASET, 1.0,1.0,
				  ss_name, &ss_sequence);
    fail_unless(rc != FC_SUCCESS, "should fail to shiftand scale to bad dataset");
    fail_unless(FC_HANDLE_EQUIV(ss_sequence, FC_NULL_SEQUENCE),
		"fail should return NULL");
    rc = fc_shiftAndScaleSequence(sequence, dataset, 1.0,0.0,
				  ss_name, &ss_sequence);
    fail_unless(rc != FC_SUCCESS, "should fail to shiftandscale with zero scale");
    fail_unless(FC_HANDLE_EQUIV(ss_sequence, FC_NULL_SEQUENCE),
		"fail should return NULL");
    rc = fc_shiftAndScaleSequence(sequence, dataset, 1.0,1.0,
				  NULL, &ss_sequence);
    fail_unless(rc != FC_SUCCESS, "should fail to shiftandscale with NULL name");
    fail_unless(FC_HANDLE_EQUIV(ss_sequence, FC_NULL_SEQUENCE),
		"fail should return NULL");
    rc = fc_shiftAndScaleSequence(sequence, dataset, 1.0,1.0,
				  ss_name, NULL);
    fail_unless(rc != FC_SUCCESS, "should fail to shiftandscale with NULL handle");
    fail_unless(FC_HANDLE_EQUIV(ss_sequence, FC_NULL_SEQUENCE),
		"fail should return NULL");
  }

  // a little more checking
  rc = fc_getNumSequence(dataset2, &numSequence);
  fail_unless(rc == FC_SUCCESS, "shouldn't fail to get sequence");
  fail_unless(numSequence == numDataType, "mismatch of numSequence");
  rc = fc_deleteDataset(dataset2);
  fail_unless(rc == FC_SUCCESS, "failed to delete local dataset");
  rc = fc_getNumSequence(dataset, &numSequence);
  fail_unless(rc == FC_SUCCESS, "shouldn't fail to get sequence");
  fail_unless(numSequence == 2*numDataType, "mismatch of numSequence");  

  // delete the dataset
  rc = fc_deleteDataset(dataset);
  fail_unless(rc == FC_SUCCESS, "failed to delete dataset at end of test");
}
END_TEST

//test making a regular sequence and intersecting reg seq
START_TEST(regseq_test)
{
  //making this self contained since im not the author of the
  //other tests in this module.

  FC_ReturnCode rc;
  FC_Dataset dset,baddset;
  FC_Sequence seq, seqb, seqc, seqd, seqe, checkSeq;
  FC_Sequence empty_seq, badseq, charseq, intseq, floatseq, nonmonoseq;
  FC_Sequence *seqarray, *returnSequences;
  FC_DataType dt;
  int numReturnSequences;
  int checkStep;
  void* coordsp;
  char charcoords[3]={'a','b','c'};
  int intcoords[3]={1,2,3};
  float floatcoords[4]={1.1,2.2,2.3,3.3};
  float nonmonocoords[3]={1.1,3.3,2.2};
  int numStep =5;
  int lnumStep;
  int j;
  

  //setup
  //make a data set to play in
  rc = fc_createDataset("temp dataset", &dset);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create dataset");

  //first tests for regular sequence

  //make a sequence
  rc = fc_createRegularSequence(dset,numStep,1.2,1.3,
				"seq_name",&seq);
  fail_unless(rc == FC_SUCCESS, "couldnt create regular sequnce");

  //testing vals for it here
  rc = fc_getSequenceNumStep(seq,&checkStep);
  fail_unless(numStep == checkStep, "created wrong number of steps");
  rc = fc_getSequenceByName(dset,"seq_name",&numReturnSequences,&returnSequences);
  fail_unless(rc == FC_SUCCESS, "unable to get created seq by name");
  fail_unless(numReturnSequences == 1, "wrong number of matching sequences");
  fail_unless(FC_HANDLE_EQUIV(returnSequences[0],seq),"mismatch of seq handles");
  checkSeq = returnSequences[0];
  free(returnSequences);

  rc = fc_getSequenceDataType(checkSeq,&dt);
  fail_unless(rc == FC_SUCCESS, "unable to get created seq data type");
  fail_unless(dt = FC_DT_DOUBLE,"created seq is wrong datatype");
  rc = fc_getSequenceCoordsPtr (checkSeq, &coordsp);
  fail_unless(rc == FC_SUCCESS, "unable to get created seq var coords");
  for (j = 0; j < numStep; j++){
    fail_unless( FC_DBL_EQUIV(((double*)coordsp)[j],1.2+j*1.3),
		 "bad seq coords");
  }
  rc = fc_deleteSequence(seq);
  fail_unless(rc == FC_SUCCESS, "unable to delete sequence");

  //make a decreasing sequence
  rc = fc_createRegularSequence(dset,numStep,1.2,-1.3,
				"seq_name",&seq);
  fail_unless(rc == FC_SUCCESS, "couldnt create regular sequnce");

  //testing vals for it here
  rc = fc_getSequenceNumStep(seq,&checkStep);
  fail_unless(numStep == checkStep, "created wrong number of steps");
  rc = fc_getSequenceCoordsPtr (seq, &coordsp);
  fail_unless(rc == FC_SUCCESS, "unable to get created seq var coords");
  for (j = 0; j < numStep; j++){
    fail_unless( FC_DBL_EQUIV(((double*)coordsp)[j],1.2-j*1.3),
		 "bad seq coords");
  }
  rc = fc_deleteSequence(seq);
  fail_unless(rc == FC_SUCCESS, "unable to delete sequence");

  //check single point sequence  with non-zero step size (should be ok)
  rc = fc_createRegularSequence(dset,1,1.2,1.3,
				"seq_name",&seq);
  fail_unless(rc == FC_SUCCESS, "couldnt create regular single step sequence");
  rc = fc_getSequenceNumStep(seq,&checkStep);
  fail_unless(checkStep == 1, "created wrong number of steps");
  rc = fc_getSequenceCoordsPtr (seq, &coordsp);
  fail_unless(rc == FC_SUCCESS, "unable to get created seq var coords");
  for (j = 0; j < checkStep; j++){
    fail_unless( FC_DBL_EQUIV(((double*)coordsp)[j],1.2),
		 "bad seq coords");
  }
  rc = fc_deleteSequence(seq);
  fail_unless(rc == FC_SUCCESS, "unable to delete sequence");


  //check single point sequence with zero step size (should be ok)
  rc = fc_createRegularSequence(dset,1,1.2,0.0,
				"seq_name",&seq);
  fail_unless(rc == FC_SUCCESS, "couldnt create regular single step sequence");
  rc = fc_getSequenceNumStep(seq,&checkStep);
  fail_unless(checkStep == 1, "created wrong number of steps");
  rc = fc_getSequenceCoordsPtr (seq, &coordsp);
  fail_unless(rc == FC_SUCCESS, "unable to get created seq var coords");
  for (j = 0; j < checkStep; j++){
    fail_unless( FC_DBL_EQUIV(((double*)coordsp)[j],1.2),
		 "bad seq coords");
  }
  rc = fc_deleteSequence(seq);
  fail_unless(rc == FC_SUCCESS, "unable to delete sequence");

  //testing bad args
  rc = fc_createRegularSequence(FC_NULL_DATASET,numStep,0,1.3,
				"seq_name",&seq);
  fail_unless(rc != FC_SUCCESS, "should fail for NULL data set");
  fail_unless(&seq != NULL, "should return NULL seq when fails");

  rc = fc_createRegularSequence(dset,0,0,1.3,
				"seq_name",&seq);
  fail_unless(rc != FC_SUCCESS, "should fail for bad numdatapoint");
  fail_unless(&seq != NULL, "should return NULL seq when fails");

  rc = fc_createRegularSequence(dset,-1,0,1.3,
				"seq_name",&seq);
  fail_unless(rc != FC_SUCCESS, "should fail for bad numdatapoint");
  fail_unless(&seq != NULL, "should return NULL seq when fails");

  rc = fc_createRegularSequence(dset,numStep,0.0,0.0,
				"seq_name",&seq);
  fail_unless(rc != FC_SUCCESS, 
	      "should fail for stepsize = 0 when numstep > 1");
  fail_unless(&seq != NULL, "should return NULL seq when fails");

  rc = fc_createRegularSequence(dset,numStep,0,1.3,
				NULL,&seq);
  fail_unless(rc != FC_SUCCESS, "should fail for NULL seq name");
  fail_unless(&seq != NULL, "should return NULL seq when fails");

  rc = fc_createRegularSequence(dset,numStep,0,1.3,
				"seq_name",NULL);
  fail_unless(rc != FC_SUCCESS, "should fail for NULL outputvariable");
  fail_unless(&seq != NULL, "should return NULL seq when fails");


  //now tests for intersecting regular sequence 
  //make some seqs
  rc = fc_createRegularSequence(dset,numStep,1.2,1.3,
				"seq_name",&seq);
  fail_unless(rc == FC_SUCCESS, "couldnt create regular sequnce");
  rc = fc_createRegularSequence(dset,numStep,0.0,1.0,
				"seqb",&seqb);
  fail_unless(rc == FC_SUCCESS, "couldnt create regular sequnce");
  rc = fc_createRegularSequence(dset,10,3.5,0.3,
				"seqc",&seqc);
  fail_unless(rc == FC_SUCCESS, "couldnt create regular sequnce");
  //these two overlap seq at only one point
  rc = fc_createRegularSequence(dset,2,1.1,0.1,
				"seqd",&seqd);
  fail_unless(rc == FC_SUCCESS, "couldnt create regular sequnce");
  rc = fc_createRegularSequence(dset,1,1.5,0.0,
				"seqe",&seqe);
  fail_unless(rc == FC_SUCCESS, "couldnt create regular sequnce");

  rc = fc_createSequence(dset,"empty",&empty_seq);
  fail_unless(rc == FC_SUCCESS, "couldnt create empty sequnce");
  //different types
  rc= fc_createSequence(dset,"charseq",&charseq);
  fail_unless(rc == FC_SUCCESS, "couldnt create char sequence");
  rc = fc_setSequenceCoords(charseq,3,FC_DT_CHAR,(void*)charcoords);
  fail_unless(rc == FC_SUCCESS, "couldnt set char sequence coords");
  rc= fc_createSequence(dset,"intseq",&intseq);
  fail_unless(rc == FC_SUCCESS, "couldnt create int sequence");
  rc = fc_setSequenceCoords(intseq,3,FC_DT_INT,(void*)intcoords);
  fail_unless(rc == FC_SUCCESS, "couldnt set int sequence coords");
  rc= fc_createSequence(dset,"floatseq",&floatseq);
  fail_unless(rc == FC_SUCCESS, "couldnt create float sequence");
  rc = fc_setSequenceCoords(floatseq,4,FC_DT_FLOAT,(void*)floatcoords);
  fail_unless(rc == FC_SUCCESS, "couldnt set float sequence coords");
  rc= fc_createSequence(dset,"nonmonoseq",&nonmonoseq);
  fail_unless(rc == FC_SUCCESS, "couldnt create nonmono sequence");
  rc = fc_setSequenceCoords(nonmonoseq,3,FC_DT_FLOAT,(void*)nonmonocoords);
  fail_unless(rc == FC_SUCCESS, "couldnt set nonmono sequence coords");


  //validity tests
  //checking only one sequence - should be same range but with
  //possibly different num steps
  lnumStep = 4;
  seqarray = (FC_Sequence*)malloc(1*sizeof(FC_Sequence));
  seqarray[0] = seq;
  rc = fc_createIntersectingRegularSequence(1,seqarray, lnumStep,
					    "junk_intersection",
					    &checkSeq);
  fail_unless(rc == FC_SUCCESS,
	      "should be able to get single array intersecting sequence");
  rc = fc_getSequenceNumStep(checkSeq,&checkStep);
  fail_unless(checkStep == lnumStep, "created wrong number of steps");
  rc = fc_getSequenceByName(dset,"junk_intersection",&numReturnSequences,&returnSequences);
  fail_unless(rc == FC_SUCCESS, "unable to get created seq by name");
  fail_unless(numReturnSequences == 1, "wrong number of matching sequences");
  fail_unless(FC_HANDLE_EQUIV(checkSeq,returnSequences[0]),"mismatch of seq handles");  
  badseq = returnSequences[0];
  free(returnSequences);

  rc = fc_getSequenceDataType(checkSeq,&dt);
  fail_unless(rc == FC_SUCCESS, "unable to get created seq data type");
  fail_unless(dt = FC_DT_DOUBLE,"created seq is wrong datatype");
  rc = fc_getSequenceCoordsPtr (checkSeq, &coordsp);
  fail_unless(rc == FC_SUCCESS, "unable to get created seq var coords");
  for (j = 0; j < checkStep; j++){
    fail_unless( FC_DBL_EQUIV(((double*)coordsp)[j],
			      1.2+(1.3*(numStep-1)*j)/(double)(lnumStep-1)),
		 "bad seq coords");
  }
  //clean up for this test - partial
  rc = fc_deleteSequence(checkSeq);
  fail_unless(rc == FC_SUCCESS, "unable to delete sequence");


  //same but for automated number of pnts. this should return
  //the same seq since seq was created by a regular seq
  rc = fc_createIntersectingRegularSequence(1,seqarray, 0,
					    "junk_intersection",
					    &checkSeq);
  fail_unless(rc == FC_SUCCESS,
	      "should be able to get single array intersecting sequence");
  rc = fc_getSequenceNumStep(checkSeq,&checkStep);
  fail_unless(checkStep == numStep, "created wrong number of steps");

  rc = fc_getSequenceByName(dset,"junk_intersection",&numReturnSequences,&returnSequences);
  fail_unless(rc == FC_SUCCESS, "unable to get created seq by name");
  fail_unless(numReturnSequences == 1, "wrong number of matching sequences");
  fail_unless(FC_HANDLE_EQUIV(checkSeq,returnSequences[0]),"mismatch of seq handles");  
  badseq = returnSequences[0];
  free(returnSequences);

  rc = fc_getSequenceDataType(checkSeq,&dt);
  fail_unless(rc == FC_SUCCESS, "unable to get created seq data type");
  fail_unless(dt = FC_DT_DOUBLE,"created seq is wrong datatype");
  rc = fc_getSequenceCoordsPtr (checkSeq, &coordsp);
  fail_unless(rc == FC_SUCCESS, "unable to get created seq var coords");
  for (j = 0; j < checkStep; j++){
    fail_unless( FC_DBL_EQUIV(((double*)coordsp)[j],1.2+1.3*j),
		 "bad seq coords");
  }
  //clean up for this test - partial
  rc = fc_deleteSequence(checkSeq);
  fail_unless(rc == FC_SUCCESS, "unable to delete sequence");

  //and now try range with single array
  rc = fc_createIntersectingRangeRegularSequence(1,seqarray, 1.2,2.6,
					    "junk_intersection",
					    &checkSeq);
  fail_unless(rc == FC_SUCCESS,
	      "should be able to get single array intersecting range sequence");
  rc = fc_getSequenceNumStep(checkSeq,&checkStep);
  fail_unless(checkStep == 3, "created wrong number of steps");
  rc = fc_getSequenceByName(dset,"junk_intersection",&numReturnSequences,&returnSequences);
  fail_unless(rc == FC_SUCCESS, "unable to get created seq by name");
  fail_unless(numReturnSequences == 1, "wrong number of matching sequences");
  fail_unless(FC_HANDLE_EQUIV(checkSeq,returnSequences[0]),"mismatch of seq handles");  
  badseq = returnSequences[0];
  free(returnSequences);
  rc = fc_getSequenceDataType(checkSeq,&dt);
  fail_unless(rc == FC_SUCCESS, "unable to get created seq data type");
  fail_unless(dt = FC_DT_DOUBLE,"created seq is wrong datatype");
  rc = fc_getSequenceCoordsPtr (checkSeq, &coordsp);
  fail_unless(rc == FC_SUCCESS, "unable to get created seq var coords");
  fail_unless( FC_DBL_EQUIV(((double*)coordsp)[0],1.2), "bad seq coords");
  fail_unless( FC_DBL_EQUIV(((double*)coordsp)[1],1.9), "bad seq coords");
  fail_unless( FC_DBL_EQUIV(((double*)coordsp)[2],2.6), "bad seq coords");

  //clean up for this test - partial
  rc = fc_deleteSequence(checkSeq);
  fail_unless(rc == FC_SUCCESS, "unable to delete sequence");

  free(seqarray);


  //checking single overlap #1
  seqarray = (FC_Sequence*)malloc(2*sizeof(FC_Sequence));
  seqarray[0] = seq;
  seqarray[1] = seqd;
  rc = fc_createIntersectingRegularSequence(2,seqarray, 1,
					    "junk_intersection",
					    &checkSeq);
  fail_unless(rc == FC_SUCCESS,
	      "should be able to get single point intersecting sequence");
  rc = fc_getSequenceNumStep(checkSeq,&checkStep);
  fail_unless(checkStep == 1, "created wrong number of steps");
  rc = fc_getSequenceCoordsPtr (checkSeq, &coordsp);
  fail_unless(rc == FC_SUCCESS, "unable to get created seq var coords");
  for (j = 0; j < checkStep; j++){
    fail_unless( FC_DBL_EQUIV(((double*)coordsp)[j],1.2),
		 "bad seq coords");
  }
  //clean up for this test - partial
  rc = fc_deleteSequence(checkSeq);
  fail_unless(rc == FC_SUCCESS, "unable to delete sequence");

  //same but for automated number of pts
  rc = fc_createIntersectingRegularSequence(2,seqarray, 0,
					    "junk_intersection",
					    &checkSeq);
  fail_unless(rc == FC_SUCCESS,
	      "should be able to get single point intersecting sequence");
  rc = fc_getSequenceNumStep(checkSeq,&checkStep);
  fail_unless(checkStep == 1, "created wrong number of steps");
  rc = fc_getSequenceCoordsPtr (checkSeq, &coordsp);
  fail_unless(rc == FC_SUCCESS, "unable to get created seq var coords");
  for (j = 0; j < checkStep; j++){
    fail_unless( FC_DBL_EQUIV(((double*)coordsp)[j],1.2),
		 "bad seq coords");
  }
  //clean up for this test - partial
  rc = fc_deleteSequence(checkSeq);
  fail_unless(rc == FC_SUCCESS, "unable to delete sequence");


  //and now range test
  rc = fc_createIntersectingRangeRegularSequence(2,seqarray, 1.2,1.2,
					    "junk_intersection",
					    &checkSeq);
  fail_unless(rc == FC_SUCCESS,
	      "should be able to get single point intersecting sequence");
  rc = fc_getSequenceNumStep(checkSeq,&checkStep);
  fail_unless(checkStep == 1, "created wrong number of steps");
  rc = fc_getSequenceCoordsPtr (checkSeq, &coordsp);
  fail_unless(rc == FC_SUCCESS, "unable to get created seq var coords");
  for (j = 0; j < checkStep; j++){
    fail_unless( FC_DBL_EQUIV(((double*)coordsp)[j],1.2),
		 "bad seq coords");
  }
  //clean up for this test - partial
  rc = fc_deleteSequence(checkSeq);
  fail_unless(rc == FC_SUCCESS, "unable to delete sequence");
  //keep seqarray for next test
  
  //checking single overlap #2
  seqarray[1] = seqe;
  rc = fc_createIntersectingRegularSequence(2,seqarray, 1,
					    "junk_intersection",
					    &checkSeq);
  fail_unless(rc == FC_SUCCESS,
	      "should be able to get single point intersecting sequence");
  rc = fc_getSequenceNumStep(checkSeq,&checkStep);
  fail_unless(checkStep == 1, "created wrong number of steps");
  rc = fc_getSequenceCoordsPtr (checkSeq, &coordsp);
  fail_unless(rc == FC_SUCCESS, "unable to get created seq var coords");
  for (j = 0; j < checkStep; j++){
    fail_unless( FC_DBL_EQUIV(((double*)coordsp)[j],1.5),
		 "bad seq coords");
  }
  //clean up for this test - partial
  rc = fc_deleteSequence(checkSeq);
  fail_unless(rc == FC_SUCCESS, "unable to delete sequence");

  //same but for automated number of pts
  rc = fc_createIntersectingRegularSequence(2,seqarray, 0,
					    "junk_intersection",
					    &checkSeq);
  fail_unless(rc == FC_SUCCESS,
	      "should be able to get single point intersecting sequence");
  rc = fc_getSequenceNumStep(checkSeq,&checkStep);
  fail_unless(checkStep == 1, "created wrong number of steps");
  rc = fc_getSequenceCoordsPtr (checkSeq, &coordsp);
  fail_unless(rc == FC_SUCCESS, "unable to get created seq var coords");
  for (j = 0; j < checkStep; j++){
    fail_unless( FC_DBL_EQUIV(((double*)coordsp)[j],1.5),
		 "bad seq coords");
  }
  //clean up for this test - partial
  rc = fc_deleteSequence(checkSeq);
  fail_unless(rc == FC_SUCCESS, "unable to delete sequence");

  //cleanup for these two tests - total
  free (seqarray);


  //checking that asking for two points gives the bounds
  seqarray = (FC_Sequence*)malloc(3*sizeof(FC_Sequence));
  seqarray[0] = seq;
  seqarray[1] = seqb;
  seqarray[2] = seqc;

  rc = fc_createIntersectingRegularSequence(3,seqarray, 2,
					    "junk_intersection",
					    &checkSeq);
  fail_unless(rc == FC_SUCCESS,
	      "should be able to get 2 point intersecting sequence");
  rc = fc_getSequenceNumStep(checkSeq,&checkStep);
  fail_unless(checkStep == 2, "created wrong number of steps");
  rc = fc_getSequenceCoordsPtr (checkSeq, &coordsp);
  fail_unless(rc == FC_SUCCESS, "unable to get created seq var coords");
  //known values
  fail_unless( FC_DBL_EQUIV(((double*)coordsp)[0],3.5),
	       "bad seq coords");
  fail_unless( FC_DBL_EQUIV(((double*)coordsp)[1],4.0),
	       "bad seq coords");
  //clean up for this test - partial
  rc = fc_deleteSequence(checkSeq);
  fail_unless(rc == FC_SUCCESS, "unable to delete sequence");
  free (seqarray);


  //checking seq if different numerical types
  seqarray = (FC_Sequence*)malloc(3*sizeof(FC_Sequence));
  seqarray[0] = seq;
  seqarray[1] = intseq;
  seqarray[2] = floatseq;

  lnumStep = 4;
  rc = fc_createIntersectingRegularSequence(3,seqarray, lnumStep,
					    "junk_intersection",
					    &checkSeq);
  fail_unless(rc == FC_SUCCESS,
	      "should be able to get  intersecting sequence");
  rc = fc_getSequenceNumStep(checkSeq,&checkStep);
  fail_unless(checkStep == lnumStep, "created wrong number of steps");
  rc = fc_getSequenceCoordsPtr (checkSeq, &coordsp);
  fail_unless(rc == FC_SUCCESS, "unable to get created seq var coords");
  //known values
  for (j = 0; j < 4; j++){
    fail_unless( FC_DBL_EQUIV(((double*)coordsp)[j],
			      1.2+j*(3.0-1.2)/(double)(lnumStep -1)),
		 "bad seq coords");
  }
  rc = fc_deleteSequence(checkSeq);
  fail_unless(rc == FC_SUCCESS, "unable to delete sequence");

  //same but with automated number of pts - this also tests calc with
  //diff number of pts in the range of the seqs in the seqarray
  rc = fc_createIntersectingRegularSequence(3,seqarray, 0,
					    "junk_intersection",
					    &checkSeq);
  fail_unless(rc == FC_SUCCESS,
	      "should be able to get  intersecting sequence");
  rc = fc_getSequenceNumStep(checkSeq,&checkStep);
  //known value
  fail_unless(checkStep == 4, "created wrong number of steps");
  rc = fc_getSequenceCoordsPtr (checkSeq, &coordsp);
  fail_unless(rc == FC_SUCCESS, "unable to get created seq var coords");
  //known values
  for (j = 0; j < checkStep; j++){
    fail_unless( FC_DBL_EQUIV(((double*)coordsp)[j],
			      1.2+j*(3.0-1.2)/(double)(checkStep -1)),
		 "bad seq coords");
  }
  rc = fc_deleteSequence(checkSeq);
  fail_unless(rc == FC_SUCCESS, "unable to delete sequence");

  //range test
  rc = fc_createIntersectingRangeRegularSequence(3,seqarray, 1.2,1.3,
						 "junk_intersection",
						 &checkSeq);
  fail_unless(rc == FC_SUCCESS,
	      "should be able to get  intersecting sequence");
  rc = fc_getSequenceNumStep(checkSeq,&checkStep);
  //known value
  fail_unless(checkStep == 2, "created wrong number of steps");
  rc = fc_getSequenceCoordsPtr (checkSeq, &coordsp);
  fail_unless(rc == FC_SUCCESS, "unable to get created seq var coords");
  //known values
  fail_unless( FC_DBL_EQUIV(((double*)coordsp)[0],1.2),
			    "bad seq coords");
  fail_unless( FC_DBL_EQUIV(((double*)coordsp)[1],1.3),
			    "bad seq coords");
  rc = fc_deleteSequence(checkSeq);
  fail_unless(rc == FC_SUCCESS, "unable to delete sequence");

  //clean up
  free (seqarray);


  //test error conditions 
  //bad args
  rc = fc_createIntersectingRegularSequence(3,NULL, numStep,
					    "junk_intersection",
					    &checkSeq);
  fail_unless(rc != FC_SUCCESS, "should fail if null seq array");
  fail_unless(&checkSeq != NULL, "should return null seq if fail");

  rc = fc_createIntersectingRangeRegularSequence(3,NULL, 3.5,3.7,
					    "junk_intersection",
					    &checkSeq);
  fail_unless(rc != FC_SUCCESS, "should fail if null seq array");
  fail_unless(&checkSeq != NULL, "should return null seq if fail");

  seqarray = (FC_Sequence*)malloc(3*sizeof(FC_Sequence));
  seqarray[0] = seq;
  seqarray[1] = seqb;
  seqarray[2] = empty_seq;

  rc = fc_createIntersectingRegularSequence(3,seqarray, 3,
					    "junk_intersection",
					    &checkSeq);
  fail_unless(rc != FC_SUCCESS, "should fail if empty seq in the seqarray");
  fail_unless(&checkSeq != NULL, "should return null seq if fail");

  rc = fc_createIntersectingRangeRegularSequence(3,seqarray, 3.5,3.7,
					    "junk_intersection",
					    &checkSeq);
  fail_unless(rc != FC_SUCCESS, "should fail if empty seq in the seqarray");
  fail_unless(&checkSeq != NULL, "should return null seq if fail");


  seqarray[2] = FC_NULL_SEQUENCE;

  rc = fc_createIntersectingRegularSequence(3,seqarray, 3,
					    "junk_intersection",
					    &checkSeq);
  fail_unless(rc != FC_SUCCESS, "should fail if NULL seq in the seqarray");
  fail_unless(&checkSeq != NULL, "should return null seq if fail");

  rc = fc_createIntersectingRangeRegularSequence(3,seqarray, 3.5,3.7,
					    "junk_intersection",
					    &checkSeq);
  fail_unless(rc != FC_SUCCESS, "should fail if NULL seq in the seqarray");
  fail_unless(&checkSeq != NULL, "should return null seq if fail");

  seqarray[2] = seqc;

  rc = fc_createIntersectingRegularSequence(3,seqarray, -1,
					    "junk_intersection",
					    &checkSeq);
  fail_unless(rc != FC_SUCCESS, "should fail if numstep < 0");
  fail_unless(&checkSeq != NULL, "should return null seq if fail");


  rc = fc_createIntersectingRangeRegularSequence(3,seqarray, 3.7,3.5,
					    "junk_intersection",
					    &checkSeq);
  fail_unless(rc != FC_SUCCESS, "should fail if seqmax < seqmin");
  fail_unless(&checkSeq != NULL, "should return null seq if fail");


  rc = fc_createIntersectingRegularSequence(3,seqarray, 3,
					    NULL,
					    &checkSeq);
  fail_unless(rc != FC_SUCCESS, "should fail if NULL new seqvar name");
  fail_unless(&checkSeq != NULL, "should return null seq if fail");

  rc = fc_createIntersectingRangeRegularSequence(3,seqarray, 3.5,3.7,
					    NULL,
					    &checkSeq);
  fail_unless(rc != FC_SUCCESS, "should fail if NULL new seqvar name");
  fail_unless(&checkSeq != NULL, "should return null seq if fail");


  rc = fc_createIntersectingRegularSequence(3,seqarray, 3,
					    "junk_intersection",
					    NULL);
  fail_unless(rc != FC_SUCCESS, "should fail if NULL return new seqvar");
  fail_unless(&checkSeq != NULL, "should return null seq if fail");

  rc = fc_createIntersectingRangeRegularSequence(3,seqarray, 3.5,3.7,
					    "junk_intersection",
					    NULL);
  fail_unless(rc != FC_SUCCESS, "should fail if NULL return new seqvar");
  fail_unless(&checkSeq != NULL, "should return null seq if fail");

  //clean up for these tests
  free (seqarray);

  //bad input 
  //not on same datset
  rc = fc_createDataset("temp dataset", &baddset);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create dataset");

  rc = fc_createRegularSequence(baddset,numStep,1.2,1.3,
				"seq_name",&badseq);
  fail_unless(rc == FC_SUCCESS, "couldnt create regular sequnce");
  
  seqarray = (FC_Sequence*)malloc(2*sizeof(FC_Sequence));
  seqarray[0] = seq;
  seqarray[1] = badseq;
  rc = fc_createIntersectingRegularSequence(2,seqarray, numStep,
					    "junk_intersection",
					    &checkSeq);
  fail_unless(rc != FC_SUCCESS, "should fail if seqs on different datasets");
  fail_unless(&checkSeq != NULL, "should return null seq if fail");

  rc = fc_createIntersectingRangeRegularSequence(2,seqarray, 3.5,3.7,
					    "junk_intersection",
					    &checkSeq);
  fail_unless(rc != FC_SUCCESS, "should fail if seqs on different datasets");
  fail_unless(&checkSeq != NULL, "should return null seq if fail");

  //clean up for this test
  free (seqarray);
  rc = fc_deleteSequence(badseq);
  fail_unless(rc == FC_SUCCESS, "failed to delete bad seq");
  rc = fc_deleteDataset(baddset);
  fail_unless(rc == FC_SUCCESS, "failed to delete bad dataset");


  //overlapping sequences of more than 1 pt but only asking for 1 pt
  seqarray = (FC_Sequence*)malloc(2*sizeof(FC_Sequence));
  seqarray[0] = seq;
  seqarray[1] = seqb;
  rc = fc_createIntersectingRegularSequence(2,seqarray, 1,
					    "junk_intersection",
					    &checkSeq);
  fail_unless(rc != FC_SUCCESS,
      "should fail for non-zero range but only 1 datapoint");
  fail_unless(&checkSeq != NULL, "should return null seq if fail");

  //range bounds bad
  rc = fc_createIntersectingRangeRegularSequence(2,seqarray, -0.5,1.5,
					    "junk_intersection",
					    &checkSeq);
  fail_unless(rc != FC_SUCCESS,
	      "should fail for range bounds too low");
  fail_unless(&checkSeq != NULL, "should return null seq if fail");

  rc = fc_createIntersectingRangeRegularSequence(2,seqarray, 1.2,10.0,
					    "junk_intersection",
					    &checkSeq);
  fail_unless(rc != FC_SUCCESS,
	      "should fail for range bounds too high");
  fail_unless(&checkSeq != NULL, "should return null seq if fail");

  //clean up for this test
  free(seqarray);

  //overlapping sequences of only 1 pt but asking for more than 1 pt
  seqarray = (FC_Sequence*)malloc(2*sizeof(FC_Sequence));
  seqarray[0] = seq;
  seqarray[1] = seqe;
  rc = fc_createIntersectingRegularSequence(2,seqarray, 3,
					    "junk_intersection",
					    &checkSeq);
  fail_unless(rc != FC_SUCCESS,
      "should fail for 1pt overlap but asking for more than 1 datapoint");
  fail_unless(&checkSeq != NULL, "should return null seq if fail");

  //clean up for this test
  free(seqarray);


  //nonoverlapping sequences
  rc = fc_createRegularSequence(dset,numStep,20.0,1.3,
				"seq_name",&badseq);
  fail_unless(rc == FC_SUCCESS, "couldnt create regular sequnce");

  seqarray = (FC_Sequence*)malloc(2*sizeof(FC_Sequence));
  seqarray[0] = seq;
  seqarray[1] = badseq;
  rc = fc_createIntersectingRegularSequence(2,seqarray, numStep,
					    "junk_intersection",
					    &checkSeq);
  fail_unless(rc != FC_SUCCESS, "should fail if non-overlapping seqs");
  fail_unless(&checkSeq != NULL, "should return null seq if fail");

  //same but with automated numpts
  rc = fc_createIntersectingRegularSequence(2,seqarray, 0,
					    "junk_intersection",
					    &checkSeq);
  fail_unless(rc != FC_SUCCESS, "should fail if non-overlapping seqs");
  fail_unless(&checkSeq != NULL, "should return null seq if fail");

  //clean up for this test
  free (seqarray);
  rc = fc_deleteSequence(badseq);
  fail_unless(rc == FC_SUCCESS, "failed to delete bad seq");

  //char seq
  seqarray = (FC_Sequence*)malloc(2*sizeof(FC_Sequence));
  seqarray[0] = seq;
  seqarray[1] = charseq;
  rc = fc_createIntersectingRegularSequence(2,seqarray, numStep,
					    "junk_intersection",
					    &checkSeq);
  fail_unless(rc != FC_SUCCESS, "should fail for char seq seqs");
  fail_unless(&checkSeq != NULL, "should return null seq if fail");

  rc = fc_createIntersectingRangeRegularSequence(2,seqarray, 1.2,1.3,
					    "junk_intersection",
					    &checkSeq);
  fail_unless(rc != FC_SUCCESS, "should fail for char seq seqs");
  fail_unless(&checkSeq != NULL, "should return null seq if fail");
  //clean up for this test
  free(seqarray);


  //checking non-monotonic seq fails
  seqarray = (FC_Sequence*)malloc(3*sizeof(FC_Sequence));
  seqarray[0] = seq;
  seqarray[1] = intseq;
  seqarray[2] = nonmonoseq;

  lnumStep = 4;
  rc = fc_createIntersectingRegularSequence(3,seqarray, lnumStep,
					    "junk_intersection",
					    &checkSeq);
  fail_unless(rc != FC_SUCCESS,
	      "should fail for non monotonic sequence array");
  fail_unless(&checkSeq != NULL, "should return null seq if fail");


  rc = fc_createIntersectingRangeRegularSequence(3,seqarray,1.3,1.4, 
					    "junk_intersection",
					    &checkSeq);
  fail_unless(rc != FC_SUCCESS,
	      "should fail for non monotonic sequence array");
  fail_unless(&checkSeq != NULL, "should return null seq if fail");
  //clean up for this test
  free(seqarray);


  //checking decreasing seq fails - NOTE changing seq for this last test
  numStep = 4;
  lnumStep = 4;
  rc = fc_deleteSequence(seq);
  fail_unless(rc == FC_SUCCESS, "failed to delete seq");
  rc = fc_deleteSequence(seqb);
  fail_unless(rc == FC_SUCCESS, "failed to delete seq");
  rc = fc_createRegularSequence(dset,numStep,1.2,-1.3,
				"seq_name",&seq);
  fail_unless(rc == FC_SUCCESS, "couldnt create regular sequnce");
  rc = fc_createRegularSequence(dset,numStep,1.1,-1.3,
				"seq_name",&seqb);
  fail_unless(rc == FC_SUCCESS, "couldnt create regular sequnce");
  //checking non-monotonic seq fails
  seqarray = (FC_Sequence*)malloc(2*sizeof(FC_Sequence));
  seqarray[0] = seq;
  seqarray[1] = seqb;

  rc = fc_createIntersectingRegularSequence(2,seqarray, lnumStep,
					    "junk_intersection",
					    &checkSeq);
  fail_unless(rc != FC_SUCCESS,
	      "should fail for decreasing sequence array");
  fail_unless(&checkSeq != NULL, "should return null seq if fail");

  rc = fc_createIntersectingRangeRegularSequence(2,seqarray, -1.0,0.0,
					    "junk_intersection",
					    &checkSeq);
  fail_unless(rc != FC_SUCCESS,
	      "should fail for decreasing sequence array");
  fail_unless(&checkSeq != NULL, "should return null seq if fail");


  //clean up for this test
  rc = fc_deleteSequence(seq);
  fail_unless(rc == FC_SUCCESS, "failed to delete seq");
  rc = fc_deleteSequence(seqb);
  fail_unless(rc == FC_SUCCESS, "failed to delete seq");
  free(seqarray);


  //clean up
  rc = fc_deleteSequence(seqc);
  fail_unless(rc == FC_SUCCESS, "failed to delete seq");
  rc = fc_deleteSequence(seqd);
  fail_unless(rc == FC_SUCCESS, "failed to delete seq");
  rc = fc_deleteSequence(seqe);
  fail_unless(rc == FC_SUCCESS, "failed to delete seq");
  rc = fc_deleteSequence(charseq);
  fail_unless(rc == FC_SUCCESS, "failed to delete seq");
  rc = fc_deleteSequence(intseq);
  fail_unless(rc == FC_SUCCESS, "failed to delete seq");
  rc = fc_deleteSequence(floatseq);
  fail_unless(rc == FC_SUCCESS, "failed to delete seq");
  rc = fc_deleteSequence(nonmonoseq);
  fail_unless(rc == FC_SUCCESS, "failed to delete seq");
  rc = fc_deleteDataset(dset);
  fail_unless(rc == FC_SUCCESS, "failed to delete dataset at end of tests");
}
END_TEST

// test error input for fc_printSequence()
START_TEST(print)
{
  FC_ReturnCode rc;
  FC_Sequence badSequence = { 999, 999 };

  // fc_printSequence()
  rc = fc_printSequence(badSequence, "Bad Sequence!", 1);
  fail_unless(rc != FC_SUCCESS, "should fail if bad sequence");
}
END_TEST 

// release seqs
START_TEST(release_seq)
{
  FC_ReturnCode rc;
  FC_Dataset dataset;
  FC_Sequence seq, seq2, badSeq = { 999, 999 };
  FC_Sequence *returnSequences;
  int numReturnSequences;
  _FC_SeqSlot* seqSlot, *seqSlot2;
  void *coords;

  // setup
  fc_loadDataset("../data/gen_hex_seq.ex2", &dataset);

  rc = fc_getSequenceByName(dataset,"time",&numReturnSequences,&returnSequences);
  fail_unless(rc == FC_SUCCESS, "unable to get created seq by name");
  fail_unless(numReturnSequences == 1, "wrong number of matching sequences");
  seq = returnSequences[0];
  free(returnSequences);
  seqSlot = _fc_getSeqSlot(seq);
  fail_unless(seqSlot != NULL, "abort: failed to get seqslot");

  // test that lazily loaded coords is not there yet
  fail_unless(seqSlot->coords == NULL, "abort: should start with no coords");

  // make copy (uncommitted) sequence
  rc = fc_copySequence(seq, dataset, "new_seq", &seq2);
  fail_unless(rc == FC_SUCCESS, "failed to copy seq");
  seqSlot2 = _fc_getSeqSlot(seq2);
  fail_unless(seqSlot2 != NULL, "abort: failed to get seqslot2");

  // --- coords

  // load & unload sequence coords
  rc = fc_getSequenceCoordsPtr(seq, &coords); // force loading
  fail_unless(rc == FC_SUCCESS, "failed to get coords");
  fail_unless(seqSlot->coords != NULL, "abort: coords should have loaded");
  fail_unless(coords == seqSlot->coords, "mismatch of coords pointers");
  rc = fc_releaseSequence(seq);
  fail_unless(rc == FC_SUCCESS, "failed to release coords");
  fail_unless(seqSlot->coords == NULL, "coords should now be null");

  // uncommitted seq's shouldn't get deleted
  fail_unless(seqSlot2->coords != NULL, "should start with coords");
  fc_releaseSequence(seq); // make sure that copy doesn't just copy pointers
  rc = fc_releaseSequence(seq2);
  fail_unless(rc == FC_SUCCESS, "failed to release coords");
  fail_unless(seqSlot2->coords != NULL, "should not delete temp seq coords");

  // --- special cases
  
  // releasing a null handle does not cause an error
  rc = fc_releaseSequence(FC_NULL_SEQUENCE);
  fail_unless(rc == FC_SUCCESS, "NULL handle should not fail");
  
  // --- errors
  
  // fc_releaseSequence() --- bad args
  rc = fc_releaseSequence(badSeq);
  fail_unless(rc != FC_SUCCESS, "bad sequence should fail");
 
  // cleanup
  fc_deleteDataset(dataset);
}
END_TEST 


// *********************************************
// ***** Populate the Suite with the tests
// *********************************************

Suite *sequence_suite(void)
{
  Suite *suite = suite_create("Sequence");

  TCase *tc_private_seq = tcase_create(" - Private Sequence ");
  TCase *tc_sequence = tcase_create(" - Sequence Interface ");
  
  // private sequence
  suite_add_tcase(suite, tc_private_seq);
  tcase_add_checked_fixture(tc_private_seq, sequence_setup, sequence_teardown);
  tcase_add_test(tc_private_seq, slot_new_delete);

  // sequence interface 
  suite_add_tcase(suite, tc_sequence);
  tcase_add_checked_fixture(tc_sequence, sequence_setup, sequence_teardown);
  tcase_add_test(tc_sequence, create_get_delete);
  tcase_add_test(tc_sequence, metadata_query);
  tcase_add_test(tc_sequence, coords_query);
  tcase_add_test(tc_sequence, get_coords_as);
  tcase_add_test(tc_sequence, copy_test);
  tcase_add_test(tc_sequence, regseq_test);
  tcase_add_test(tc_sequence, release_seq);
  tcase_add_test(tc_sequence, print);

  return suite;
}
