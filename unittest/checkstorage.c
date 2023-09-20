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
 * \file checkstorage.c
 * \brief Unit testing of \ref SimpleDataObjects Module.
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/unittest/checkstorage.c,v $
 * $Revision: 1.12 $ 
 * $Date: 2006/08/30 19:20:05 $
 *
 * \modifications
 *    - 02/25/04 WSK, created.
 *    - 05/07/04 WSK moved checking of mesh structures to checkmesh.c,
 *      moved checking of linked lists and masks from checkutil to here.
 */

#include <stdlib.h>
#include <stdio.h>
#include <check.h>
#include <string.h>
#include "fc.h"
#include "storageP.h"
#include "checkall.h"

// sorted int arrays

START_TEST(sia_helpers) 
{
  FC_ReturnCode rc;
  FC_SortedIntArray sia;
  int i, flag, idx, midID = 3;
  int numVal = 7, sortedVals[7] = { -9, -2, 0, 3, 6, 11, 14 };
  int val = 99;
  int temp_length, *temp_array;

  // -- test _fc_expandIntArray()

  // length can't be negative
  temp_length = -1;
  temp_array = (int*)1;
  rc = _fc_expandIntArray(&temp_length, &temp_array);
  fail_unless(rc != FC_SUCCESS, "should fail if length < 0");
  fail_unless(temp_array == (int*)1, "should not change array");
  // if 0, get length of 1
  temp_length = 0;
  rc = _fc_expandIntArray(&temp_length, &temp_array);
  fail_unless(rc == FC_SUCCESS, "should not fial if length 0");
  fail_unless(temp_length == 1, "should be of size 1");
  fail_unless(temp_array != NULL && temp_array != (int*)1, 
	      "should have changed array");
  temp_array[0] = sortedVals[0];
  rc = _fc_expandIntArray(&temp_length, &temp_array);
  fail_unless(rc == FC_SUCCESS, "should not fail if length > 0");
  fail_unless(temp_length == 2, "should be of size 2");
  fail_unless(temp_array[0] == sortedVals[0], 
	      "first val should be the same");
  temp_array[1] = sortedVals[1];
  rc = _fc_expandIntArray(&temp_length, &temp_array);
  fail_unless(rc == FC_SUCCESS, "should not fail if length > 0");
  fail_unless(temp_length == 4, "should be of size 4");
  fail_unless(!memcmp(temp_array, sortedVals, 2*sizeof(int)), 
	      "vals should be the same");
  temp_array[2] = sortedVals[2];
  temp_array[3] = sortedVals[3];
  rc = _fc_expandIntArray(&temp_length, &temp_array);
  fail_unless(rc == FC_SUCCESS, "should not fail if length > 0");
  fail_unless(temp_length == 8, "should be of size 8");
  fail_unless(!memcmp(temp_array, sortedVals, 4*sizeof(int)), 
	      "vals should be the same");
  free(temp_array);

  // -- test _fc_lookupIntInSortedIntArray

  // test lookup = binary search algrithm (hack an sia)
  sia.numVal = numVal;
  sia.maxNumVal = numVal;
  sia.nextTicket = numVal;
  sia.vals = malloc(numVal*sizeof(int));
  sia.tickets = malloc(numVal*sizeof(int));
  fail_unless(sia.vals != NULL, "abort: memory allocation failure");
  fail_unless(sia.tickets != NULL, "abort: memory allocation failure");
  memcpy(sia.vals, sortedVals, numVal*sizeof(int));
  for(i=0; i<numVal; i++)
    sia.tickets[i]=i;
  for (i = 0; i < numVal; i++) {
    rc = _fc_lookupIntInSortedIntArray(&sia, sortedVals[i], &flag, &idx);
    fail_unless(rc == FC_SUCCESS, "failed to get index");
    fail_unless(flag == 1 && idx == i, "mismatch of found index");
  }
  for (i = 0; i < numVal; i++) {
    rc = _fc_lookupIntInSortedIntArray(&sia, sortedVals[i]-1, &flag, &idx);
    fail_unless(rc == FC_SUCCESS, "failed to get index");
    fail_unless(flag == 0 && idx == i, "mismatch of not found index");
  }
  for (i = 0; i < numVal; i++) {
    rc = _fc_lookupIntInSortedIntArray(&sia, sortedVals[i]+1, &flag, &idx);
    fail_unless(rc == FC_SUCCESS, "failed to get index");
    fail_unless(flag == 0 && idx == i+1, "mismatch of not found index");
  }

  // test bad args
  rc = _fc_lookupIntInSortedIntArray(NULL, 5, &flag, &idx);
  fail_unless(rc != FC_SUCCESS, "should fail for null input");
  fail_unless(flag == -1 && idx == -1, "fail should return nulls");
  rc = _fc_lookupIntInSortedIntArray(&sia, 5, NULL, &idx);
  fail_unless(rc != FC_SUCCESS, "should fail for null input");
  fail_unless(idx == -1, "fail should return nulls");
  rc = _fc_lookupIntInSortedIntArray(&sia, 5, &flag, NULL);
  fail_unless(rc != FC_SUCCESS, "should fail for null input");
  fail_unless(flag == -1, "fail should return nulls");

  // -- test _fc_addEntryToSortedIntArray & _fc_deleteEntryFromSortedIntArray

  // test add/delete from front
  rc = _fc_addEntryToSortedIntArray(&sia, 0, val);
  fail_unless(rc == FC_SUCCESS, "failed to add");
  fail_unless(sia.numVal == numVal+1, "mismatch of val");
  fail_unless(sia.vals[0] == val, "mismatch of val");
  fail_unless(!memcmp(sia.vals+1, sortedVals, numVal*sizeof(int)),
	      "mismatch of vals");
  rc = _fc_deleteEntryFromSortedIntArray(&sia, 0);
  fail_unless(rc == FC_SUCCESS, "failed to add");
  fail_unless(sia.numVal == numVal, "mismatch of val");
  fail_unless(!memcmp(sia.vals, sortedVals, numVal*sizeof(int)),
	      "mismatch of vals");

  // test add/delete from middle
  rc = _fc_addEntryToSortedIntArray(&sia, midID, val);
  fail_unless(rc == FC_SUCCESS, "failed to add");
  fail_unless(sia.numVal == numVal+1, "mismatch of val");
  fail_unless(sia.vals[midID] == val, "mismatch of val");
  fail_unless(!memcmp(sia.vals, sortedVals, midID*sizeof(int)),
	      "mismatch of vals");
  fail_unless(!memcmp(sia.vals+midID+1, sortedVals+midID, 
		      (numVal-midID)*sizeof(int)),
  	      "mismatch of vals");
  rc = _fc_deleteEntryFromSortedIntArray(&sia, midID);
  fail_unless(rc == FC_SUCCESS, "failed to add");
  fail_unless(sia.numVal == numVal, "mismatch of val");
  fail_unless(!memcmp(sia.vals, sortedVals, numVal*sizeof(int)),
	      "mismatch of vals");

  // test add/delete from back
  rc = _fc_addEntryToSortedIntArray(&sia, numVal, val);
  fail_unless(rc == FC_SUCCESS, "failed to add");
  fail_unless(sia.numVal == numVal+1, "mismatch of val");
  fail_unless(sia.vals[numVal] == val, "mismatch of val");
  fail_unless(!memcmp(sia.vals, sortedVals, numVal*sizeof(int)),
	      "mismatch of vals");
  rc = _fc_deleteEntryFromSortedIntArray(&sia, numVal);
  fail_unless(rc == FC_SUCCESS, "failed to add");
  fail_unless(sia.numVal == numVal, "mismatch of val");
  fail_unless(!memcmp(sia.vals, sortedVals, numVal*sizeof(int)),
	      "mismatch of vals");

  // test bad args
  rc = _fc_addEntryToSortedIntArray(NULL, 0, val);
  fail_unless(rc != FC_SUCCESS, "should fail w/ null args");
  rc = _fc_deleteEntryFromSortedIntArray(NULL, 0);
  fail_unless(rc != FC_SUCCESS, "should fail w/ null args");
  rc = _fc_addEntryToSortedIntArray(&sia, -1, val);
  fail_unless(rc != FC_SUCCESS, "should fail w/ negative index");
  rc = _fc_deleteEntryFromSortedIntArray(&sia, -1);
  fail_unless(rc != FC_SUCCESS, "should fail w/ negative index");
  // not "too big index" is not the same for the two routines
  rc = _fc_addEntryToSortedIntArray(&sia, numVal+1, val);
  fail_unless(rc != FC_SUCCESS, "should fail w/ index too big");
  rc = _fc_deleteEntryFromSortedIntArray(&sia, numVal);
  fail_unless(rc != FC_SUCCESS, "should fail w/ index too big");

  // cleanup
  free(sia.vals);
  free(sia.tickets);

  // test border cases: empty arrays w & w/o vals=NULL
  // & also that growing works
  sia.numVal = 0;
  sia.maxNumVal = 0;
  sia.nextTicket = 36;
  sia.vals = NULL;
  sia.tickets = NULL;
  rc = _fc_deleteEntryFromSortedIntArray(&sia, 0);
  fail_unless(rc != FC_SUCCESS, "should fail to delete from empty");
  rc = _fc_addEntryToSortedIntArray(&sia, 0, val);
  fail_unless(rc == FC_SUCCESS, "should add to empty");
  fail_unless(sia.numVal == 1, "mismatch of numVal");
  fail_unless(sia.vals[0] = val, "mismatch of vals");
  fail_unless(sia.tickets[0] == 36, "insert ticket didn't match");
  fail_unless(sia.nextTicket == 37, "current ticket didn't update");
  for (i = 0; i < 10; i++) { // add a bunch to test growing
    rc = _fc_addEntryToSortedIntArray(&sia, 0, val);
    fail_unless(rc == FC_SUCCESS, "should add to empty");
  }
  fail_unless(sia.numVal == 11, "mismatch of numVal");
  sia.numVal = 0;
  sia.vals[0] = val*2;
  rc = _fc_deleteEntryFromSortedIntArray(&sia, 0);
  fail_unless(rc != FC_SUCCESS, "should fail to delete from empty");
  rc = _fc_addEntryToSortedIntArray(&sia, 0, val);
  fail_unless(rc == FC_SUCCESS, "should add to empty");
  fail_unless(sia.numVal == 1, "mismatch of numVal");
  fail_unless(sia.vals[0] = val, "mismatch of vals");

  // cleanup again
  free(sia.vals);
  free(sia.tickets);
}
END_TEST

START_TEST(sia_valid) 
{
  FC_SortedIntArray sia;
  int flag;
  int numVal = 7;
  int sortedVals[7] = { -9, -2, 0, 3, 6, 11, 14 };
  int insertIDs[7]  = {  6,  5, 4, 3, 2,  1,  0 }; //Reverse order
  int next_ticket = 7;
  // hack the sia & test

  // valid cases
  // all zeros 
  sia.numVal = 0;
  sia.maxNumVal = 0;
  sia.nextTicket = 0;
  sia.vals = NULL;
  sia.tickets = NULL;
  flag = fc_isSortedIntArrayValid(&sia);
  fail_unless(flag == 1, "should return true");
  // no entries but has space
  sia.numVal = 0;
  sia.maxNumVal = numVal;
  sia.nextTicket = next_ticket;
  sia.vals = sortedVals;
  sia.tickets = insertIDs;
  flag = fc_isSortedIntArrayValid(&sia);
  fail_unless(flag == 1, "should return true");
  // has same number of entries and space
  sia.numVal = numVal;
  sia.maxNumVal = numVal;
  sia.nextTicket = next_ticket;
  sia.vals = sortedVals;
  sia.tickets = insertIDs;
  flag = fc_isSortedIntArrayValid(&sia);
  fail_unless(flag == 1, "should return true");
  // has fewer entires than space
  sia.numVal = 1;
  sia.maxNumVal = numVal;
  sia.nextTicket = next_ticket;
  sia.vals = sortedVals;
  sia.tickets = insertIDs;
  flag = fc_isSortedIntArrayValid(&sia);
  fail_unless(flag == 1, "should return true");
  
  // invalid cases
  // numVal is negative
  sia.numVal = -1;
  sia.maxNumVal = 0;
  sia.nextTicket = 0;
  sia.vals = NULL;
  sia.tickets = NULL;
  flag = fc_isSortedIntArrayValid(&sia);
  fail_unless(flag == 0, "should return false");
  // maxNumVal is negative
  sia.numVal = 0;
  sia.maxNumVal = -1;
  sia.nextTicket = 0;
  sia.vals = NULL;
  sia.tickets = NULL;
  flag = fc_isSortedIntArrayValid(&sia);
  fail_unless(flag == 0, "should return false");
  // space is zero, but there appears to be something there
  sia.numVal = 0;
  sia.maxNumVal = 0;
  sia.nextTicket = 0;
  sia.vals = sortedVals;
  sia.tickets = insertIDs;
  flag = fc_isSortedIntArrayValid(&sia);
  fail_unless(flag == 0, "should return false");
  // there should be something there, but there is not
  sia.numVal = numVal;
  sia.maxNumVal = numVal;
  sia.nextTicket = next_ticket;
  sia.vals = NULL;
  sia.tickets = NULL;
  flag = fc_isSortedIntArrayValid(&sia);
  fail_unless(flag == 0, "should return false");
  // there are more entries than space
  sia.numVal = numVal+1;
  sia.maxNumVal = numVal;
  sia.nextTicket = next_ticket;
  sia.vals = sortedVals;
  sia.tickets = insertIDs;
  flag = fc_isSortedIntArrayValid(&sia);
  fail_unless(flag == 0, "should return false");
}
END_TEST

START_TEST(sia_basic) 
{
  FC_ReturnCode rc;
  FC_SortedIntArray sia;
  int i, count;
  // have enough vals to cause multiple grows, but do not
  // have a multiple of 2 so can test squeeze
  int numVal = 7, unsortedVals[7] = { 3, 14, -9, 0, 6, 11, -2 };
  // note the vals are all at least 2 apart so I can easily, rigoursly
  // test the binary search algorithm by using vals[i]+/-1.
  int sortedVals[7] = { -9, -2, 0, 3, 6, 11, 14 };
  // have one "in front", "in middle", and "in back" of vals
  int numNoVal = 3, valsNotInArray[3] = { -99, 1, 99 };
  int temp_numVal, *temp_vals;

  // test fc_initSortedIntArray
  sia.numVal = sia.maxNumVal = -2;
  sia.vals = sortedVals;
  rc = fc_initSortedIntArray(&sia);
  fail_unless(rc == FC_SUCCESS, "bad return");
  fail_unless(sia.numVal == 0, "bad val from init");
  fail_unless(sia.maxNumVal == 0, "bad val from init");
  fail_unless(sia.vals == NULL, "bad val from init");

  // test query & delete of empty & NULL sia
  count = fc_isIntInSortedIntArray(&sia, 5);
  fail_unless(count == 0, "should not find vals in empty sia");
  rc = fc_getSortedIntArrayNumValue(&sia, &temp_numVal);
  fail_unless(rc == FC_SUCCESS, "should find numVal");
  fail_unless(temp_numVal == 0, "numVal should be zero for empty");
  rc = fc_getSortedIntArrayValues(&sia, &temp_numVal, &temp_vals);
  fail_unless(rc == FC_SUCCESS, "should find values");
  fail_unless(temp_numVal == 0, "numVal should be zero for empty");
  fail_unless(temp_vals == NULL, "vals should be zero for empty");
  count = fc_deleteIntFromSortedIntArray(&sia, 5);
  fail_unless(count == 0, "should not delete vals from empty sia");
 
  // test add single int
  for (i = 0; i < numVal; i++) {
    count = fc_addIntToSortedIntArray(&sia, unsortedVals[i]);
    fail_unless(count == 1, "create should create");
    fail_unless(sia.numVal == i+1, "sorted wrong");
  }
  fail_unless(sia.numVal == numVal, "bad add");
  fail_unless(!memcmp(sia.vals, sortedVals, numVal*sizeof(int)),
	      "sorted wrong");
 
  // test query
  for (i = 0; i < numVal; i++) {
    count = fc_isIntInSortedIntArray(&sia, unsortedVals[i]);
    fail_unless(count == 1, "should find vals");
  }
  for (i = 0; i < numNoVal; i++) {
    count = fc_isIntInSortedIntArray(&sia, valsNotInArray[i]);
    fail_unless(count == 0, "should not find not there vals");
  }
  rc = fc_getSortedIntArrayNumValue(&sia, &temp_numVal);
  fail_unless(rc == FC_SUCCESS, "should find numVal");
  fail_unless(temp_numVal == numVal, "numVal should agree");
  rc = fc_getSortedIntArrayValues(&sia, &temp_numVal, &temp_vals);
  fail_unless(rc == FC_SUCCESS, "should find values");
  fail_unless(temp_numVal == numVal, "numVal mismatch");
  fail_unless(!memcmp(temp_vals, sortedVals, numVal*sizeof(int)),
	      "values mismatch");
  free(temp_vals);
  fail_unless(sia.numVal == numVal, "queries altered contents");
  fail_unless(!memcmp(sia.vals, sortedVals, numVal*sizeof(int)),
	      "queries alterned contents");

  //Test getting the values back in inserted order
  rc = fc_getSortedIntArrayValuesInInsertOrder(&sia, &temp_numVal, &temp_vals);
  fail_unless(rc == FC_SUCCESS, "should find values");
  fail_unless(temp_numVal == numVal, "should have had same number of insert items");
  fail_unless(!memcmp(temp_vals, unsortedVals, numVal*sizeof(int)), "results differ when getting vals in insert order");
  free(temp_vals);



  // test try to create again (add same ones)
  for (i = 0; i < numVal; i++) {
    count = fc_addIntToSortedIntArray(&sia, unsortedVals[i]);
    fail_unless(count == 0, "create doesn't create if already there");
  }
  fail_unless(sia.numVal == numVal, "recreating went wrong");
  fail_unless(!memcmp(sia.vals, sortedVals, numVal*sizeof(int)),
	      "recreating went wrong");

  //Should still have same insert order
  rc = fc_getSortedIntArrayValuesInInsertOrder(&sia, &temp_numVal, &temp_vals);
  fail_unless(rc == FC_SUCCESS, "should find values");
  fail_unless(temp_numVal == numVal, "should have had same number of insert items");
  fail_unless(!memcmp(temp_vals, unsortedVals, numVal*sizeof(int)), "results differ from getting vals in insert order");
  free(temp_vals);


  // test delete single int
  for (i = 0; i < numNoVal; i++) {
    count = fc_deleteIntFromSortedIntArray(&sia, valsNotInArray[i]);
    fail_unless(count == 0, "should not delete not there vals");
  }
  fail_unless(sia.numVal == numVal, "recreating went wrong");
  //Should still have same insert order
  rc = fc_getSortedIntArrayValuesInInsertOrder(&sia, &temp_numVal, &temp_vals);
  fail_unless(rc == FC_SUCCESS, "should find values");
  fail_unless(temp_numVal == numVal, "should have had same number of insert items");
  fail_unless(!memcmp(temp_vals, unsortedVals, numVal*sizeof(int)), "results differ from getting vals in insert order");
  free(temp_vals);

  for (i = 0; i < numVal; i++) {

    count = fc_deleteIntFromSortedIntArray(&sia, unsortedVals[i]);
    fail_unless(count == 1, "should delete vals");
    fail_unless(sia.numVal == numVal-i-1, "didn't delete vals");

    //Should still have same insert order, just with the first n vals removed
    rc = fc_getSortedIntArrayValuesInInsertOrder(&sia, &temp_numVal, &temp_vals);
    fail_unless(rc == FC_SUCCESS, "should find values");
    fail_unless(temp_numVal == numVal-(i+1), "should have had same number of insert items");
    if(i!=numVal-1){
      fail_unless(!memcmp(temp_vals, &unsortedVals[i+1], (numVal-(i+1))*sizeof(int)), "results differ from getting vals in insert order");
    } else {
      fail_unless(temp_vals == NULL, "all deleted should have null unsorted list");
    }
    free(temp_vals);
  }
  fail_unless(sia.numVal == 0, "should be empty");
  fail_unless(sia.maxNumVal != 0, "but still have space");
  fail_unless(sia.vals != NULL, "but still have space");
  fail_unless(sia.tickets != NULL, "but still have space");

  // test query and delete of empty, but not NULL sia
  count = fc_deleteIntFromSortedIntArray(&sia, 5);
  fail_unless(count == 0, "should not delete vals from empty sia");
  count = fc_isIntInSortedIntArray(&sia, 5);
  fail_unless(count == 0, "should not find vals in empty sia");
  rc = fc_getSortedIntArrayNumValue(&sia, &temp_numVal);
  fail_unless(rc == FC_SUCCESS, "should find numVal");
  fail_unless(temp_numVal == 0, "empty sia should have no values");
  rc = fc_getSortedIntArrayValues(&sia, &temp_numVal, &temp_vals);
  fail_unless(rc == FC_SUCCESS, "should find values");
  fail_unless(temp_numVal == 0, "numVal mismatch");
  fail_unless(temp_vals == NULL, "values mismatch");
  
  // special cases to test memcopies
  i = fc_deleteIntFromSortedIntArray(&sia,4);
  fail_unless(i == 0, "bad return val from delete");
  i = fc_isIntInSortedIntArray(&sia,4);
  fail_unless(i == 0, "bad return val from lookup");
  i = fc_addIntToSortedIntArray(&sia,4);
  fail_unless(i == 1, "bad return val from lookup/create");
  fail_unless(sia.numVal == 1, "bad val from lookup/create");
  fail_unless(sia.vals[0] == 4, "bad val from lookup/create");
  i = fc_isIntInSortedIntArray(&sia,4);
  fail_unless(i == 1, "bad return val from lookup");
  fail_unless(sia.numVal == 1, "bad val from lookup");
  fail_unless(sia.vals[0] == 4, "bad val from lookup");
  i = fc_addIntToSortedIntArray(&sia,4);
  fail_unless(i == 0, "bad return val from lookup/create");
  fail_unless(sia.numVal == 1, "bad val from lookup/create");
  fail_unless(sia.vals[0] == 4, "bad val from lookup/create");

  i = fc_isIntInSortedIntArray(&sia,3);
  fail_unless(i == 0, "bad return val from lookup");
  fail_unless(sia.numVal == 1, "bad val from lookup");
  fail_unless(sia.vals[0] == 4, "bad val from lookup");
  i = fc_addIntToSortedIntArray(&sia,3);
  fail_unless(i == 1, "bad return val from lookup/create");
  fail_unless(sia.numVal == 2, "bad val from lookup");
  fail_unless(sia.vals[0] == 3, "bad val from lookup");
  fail_unless(sia.vals[1] == 4, "bad val from lookup");

  i = fc_addIntToSortedIntArray(&sia,5);
  fail_unless(i == 1, "bad return val from lookup/create");
  fail_unless(sia.numVal == 3, "bad val from lookup");
  fail_unless(sia.vals[0] == 3, "bad val from lookup");
  fail_unless(sia.vals[1] == 4, "bad val from lookup");
  fail_unless(sia.vals[2] == 5, "bad val from lookup");

  i = fc_deleteIntFromSortedIntArray(&sia,4);
  fail_unless(i == 1, "bad return val from delete");
  fail_unless(sia.numVal == 2, "bad val from delete");
  fail_unless(sia.vals[0] == 3, "bad val from delete");
  fail_unless(sia.vals[1] == 5, "bad val from delete");

  i = fc_addIntToSortedIntArray(&sia,4);
  i = fc_deleteIntFromSortedIntArray(&sia,5);
  fail_unless(sia.numVal == 2, "bad val from delete");
  fail_unless(sia.vals[0] == 3, "bad val from delete");
  fail_unless(sia.vals[1] == 4, "bad val from delete");
  i = fc_addIntToSortedIntArray(&sia,5);
  fail_unless(sia.numVal == 3, "bad val from lookup/create");
  fail_unless(sia.vals[0] == 3, "bad val from lookup/create");
  fail_unless(sia.vals[1] == 4, "bad val from lookup/create");
  fail_unless(sia.vals[2] == 5, "bad val from lookup/create");
  // delay cleanup until after bad args test
   
  // test bad args
  // fc_initSortedIntArray
  rc = fc_initSortedIntArray(NULL);
  fail_unless(rc != FC_SUCCESS, "should fail for null input");
  // fc_getSortedIntArrayNumValue
  rc = fc_getSortedIntArrayNumValue(NULL, &temp_numVal);
  fail_unless(rc != FC_SUCCESS, "should fail for null input");
  fail_unless(temp_numVal == -1, "fail should return nulls");
  rc = fc_getSortedIntArrayNumValue(&sia, NULL);
  fail_unless(rc != FC_SUCCESS, "should fail for null input");
  // fc_getSortedIntArrayValues
  rc = fc_getSortedIntArrayValues(NULL, &temp_numVal, &temp_vals);
  fail_unless(rc != FC_SUCCESS, "should fail for null input");
  fail_unless(temp_numVal == -1, "fail should return nulls");
  fail_unless(temp_vals == NULL, "fail should return nulls");
  rc = fc_getSortedIntArrayValues(&sia, NULL, &temp_vals);
  fail_unless(rc != FC_SUCCESS, "should fail for null input");
  fail_unless(temp_vals == NULL, "fail should return nulls");
  rc = fc_getSortedIntArrayValues(&sia, &temp_numVal, NULL);
  fail_unless(rc != FC_SUCCESS, "should fail for null input");
  fail_unless(temp_numVal == -1, "fail should return nulls");
  // fc_addIntToSortedIntArray
  count = fc_addIntToSortedIntArray(NULL, 5);
  fail_unless(count < FC_SUCCESS, "should fail for null input");
  // fc_isIntInSortedIntArray
  count = fc_isIntInSortedIntArray(NULL, 5);
  fail_unless(count < FC_SUCCESS, "should fail for null input");
  // fc_deleteIntFromSortedIntArray
  count = fc_deleteIntFromSortedIntArray(NULL, 5);
  fail_unless(count < FC_SUCCESS, "should fail for null input");
  // fc_printSortedIntArray
  rc = fc_printSortedIntArray(NULL);
  fail_unless(rc != FC_SUCCESS, "should fail for null input");
  // fc_getSortedIntArrayValuesInInsertOrder
  rc = fc_getSortedIntArrayValuesInInsertOrder(NULL, &temp_numVal, &temp_vals);
  fail_unless(rc != FC_SUCCESS, "should fail on null input");
  fail_unless(temp_numVal == -1, "fail should return nulls");
  fail_unless(temp_vals == NULL, "fail should return nulls");
  rc = fc_getSortedIntArrayValuesInInsertOrder(&sia, NULL, &temp_vals);
  fail_unless(rc != FC_SUCCESS, "should fail for null input");
  fail_unless(temp_vals == NULL, "fail should return nulls");
  rc = fc_getSortedIntArrayValuesInInsertOrder(&sia, &temp_numVal, NULL);
  fail_unless(rc != FC_SUCCESS, "should fail for null input");
  fail_unless(temp_numVal == -1, "fail should return nulls");
  // fc_freeSortedIntArray
  fc_freeSortedIntArray(NULL);


  // test freeing 
  fc_freeSortedIntArray(&sia);
  fail_unless(sia.numVal == 0, "bad val from free");
  fail_unless(sia.maxNumVal == 0, "bad val from free");
  fail_unless(sia.vals == NULL, "bad val from free");
}
END_TEST

START_TEST(sia_copy) 
{
  FC_ReturnCode rc;
  int i;
  FC_SortedIntArray sia_src, sia_dest;
  int numVal = 7, sortedVals[7] = { -9, -2, 0, 3, 6, 11, 14 };
  int temp_numVal, *temp_vals;

  // test copy of empty & NULL
  rc = fc_initSortedIntArray(&sia_src);
  fail_unless(rc == FC_SUCCESS, "abort: failed to init sia");
  rc = fc_copySortedIntArray(&sia_src, &sia_dest);
  fail_unless(rc == FC_SUCCESS, "failed to copy empty");
  fail_unless(sia_dest.numVal == 0 && 
	      sia_dest.maxNumVal == sia_dest.numVal,
	      "mismatch of numVal & maxNumVal");
  fail_unless(sia_dest.vals == NULL, "mismatch of vals");

  // test copy of non empty
  for (i = 0; i < numVal; i++)
    fc_addIntToSortedIntArray(&sia_src, sortedVals[i]);
  fail_unless(rc >= FC_SUCCESS, "abort: failed to add vals to sia");
  rc = fc_copySortedIntArray(&sia_src, &sia_dest);
  fail_unless(rc == FC_SUCCESS, "failed to copy empty");
  fail_unless(sia_dest.numVal == numVal && 
	      sia_dest.maxNumVal == sia_dest.numVal,
	      "mismatch of numVal & maxNumVal");
  fail_unless(!memcmp(sia_dest.vals, sortedVals, numVal*sizeof(int)), 
	      "mismatch of vals");
  //Make sure the insert order is still ok on both dest and src
  rc = fc_getSortedIntArrayValuesInInsertOrder(&sia_dest, &temp_numVal, &temp_vals);
  fail_unless(rc == FC_SUCCESS, "Failed to get in insert order");
  fail_unless(temp_numVal == numVal, "should have had same number of insert items");
  fail_unless(!memcmp(temp_vals, sortedVals, numVal*sizeof(int)), "results differ when getting vals in insert order");
  free(temp_vals);
  rc = fc_getSortedIntArrayValuesInInsertOrder(&sia_src, &temp_numVal, &temp_vals);
  fail_unless(rc == FC_SUCCESS, "Failed to get in insert order");
  fail_unless(temp_numVal == numVal, "should have had same number of insert items");
  fail_unless(!memcmp(temp_vals, sortedVals, numVal*sizeof(int)), "results differ when getting vals in insert order");
  free(temp_vals);
  fc_freeSortedIntArray(&sia_dest);

  // test copy of empty & not NULL
  for (i = 0; i < numVal; i++)
    fc_deleteIntFromSortedIntArray(&sia_src, sortedVals[i]);
  fail_unless(rc >= FC_SUCCESS, "abort: failed to add vals to sia");
  rc = fc_copySortedIntArray(&sia_src, &sia_dest);
  fail_unless(rc == FC_SUCCESS, "failed to copy empty");
  fail_unless(sia_dest.numVal == 0 && 
	      sia_dest.maxNumVal == sia_dest.numVal,
	      "mismatch of numVal & maxNumVal");
  fail_unless(sia_dest.vals == NULL, "mismatch of vals");

  // test bad args
  rc = fc_copySortedIntArray(NULL, &sia_dest);
  fail_unless(rc != FC_SUCCESS, "should fail for bad args");
  rc = fc_copySortedIntArray(&sia_src, NULL);
  fail_unless(rc != FC_SUCCESS, "should fail for bad args");

  // cleanup
  fc_freeSortedIntArray(&sia_src);
}
END_TEST

START_TEST(sia_squeeze) 
{
  FC_ReturnCode rc;
  int i;
  FC_SortedIntArray sia;
  // do not have a multiple of 2 so can test squeeze
  int numVal = 7, sortedVals[7] = { -9, -2, 0, 3, 6, 11, 14 };


  // test squeeze of empty & NULL
  rc = fc_initSortedIntArray(&sia);
  fail_unless(rc == FC_SUCCESS, "abort: failed to init sia");
  rc = fc_squeezeSortedIntArray(&sia);
  fail_unless(rc == FC_SUCCESS, "failed to squeeze");
  fail_unless(sia.numVal == 0, "mismatch of numVal");
  fail_unless(sia.maxNumVal == 0, "mismatch of maxNumVal");
  fail_unless(sia.vals == NULL, "mismatch of vals");
  fail_unless(sia.tickets == NULL, "mismatch of tickets");


  // test squeeze of non empty & extra room
  for (i = 0; i < numVal; i++) {
    rc = fc_addIntToSortedIntArray(&sia, sortedVals[i]);
    fail_unless(rc >= FC_SUCCESS, "abort: failed to add vals to sia");
  }
  fail_unless(sia.maxNumVal > sia.numVal, "abort: test not setup correctly");
  rc = fc_squeezeSortedIntArray(&sia);
  fail_unless(rc == FC_SUCCESS, "failed to squeeze");
  fail_unless(sia.numVal == numVal, "mismatch of numVal");
  fail_unless(sia.maxNumVal == numVal, "mismatch of maxNumVal");
  fail_unless(!memcmp(sia.vals, sortedVals, numVal*sizeof(int)),
	      "mismatch of vals");

  // test squeeze of empty & extra room
  fail_unless(sia.maxNumVal == sia.numVal && sia.numVal > 0, 
	      "abort: test not setup correctly");
  rc = fc_squeezeSortedIntArray(&sia);
  fail_unless(rc == FC_SUCCESS, "failed to squeeze");
  fail_unless(sia.numVal == numVal, "mismatch of numVal");
  fail_unless(sia.maxNumVal == numVal, "mismatch of maxNumVal");
  fail_unless(!memcmp(sia.vals, sortedVals, numVal*sizeof(int)),
	      "mismatch of vals");

  // test copy of empty & not NULL
  for (i = 0; i < numVal; i++)
    fc_deleteIntFromSortedIntArray(&sia, sortedVals[i]);
  fail_unless(sia.maxNumVal > 0 && sia.numVal == 0, 
	      "abort: test not setup correctly");
  rc = fc_squeezeSortedIntArray(&sia);
  fail_unless(rc == FC_SUCCESS, "failed to squeeze");
  fail_unless(sia.numVal == 0, "mismatch of numVal");
  fail_unless(sia.maxNumVal == 0, "mismatch of maxNumVal");
  fail_unless(sia.vals == NULL, "mismatch of vals");

  // test bad args
  rc = fc_squeezeSortedIntArray(NULL);
  fail_unless(rc != FC_SUCCESS, "should fail for null input");

  // cleanup not needed! (if squeeze is working)
}
END_TEST

START_TEST(sia_multiadd) 
{
  FC_SortedIntArray sia;
  int i, count;
  // have enough vals to cause multiple grows, but do not
  // have a multiple of 2 so can test squeeze
  int numVal = 7, unsortedVals[7] = { 3, 14, -9, 0, 6, 11, -2 };
  // this set has duplicates
  int numVal2 = 10, unsortedVals2[10] = { 3, 3, 3, 14, -9 , 0, 0, 6, 11, -2 };
  // note the vals are all at least 2 apart so I can easily, rigoursly
  // test the binary search algorithm by using vals[i]+/-1.
  int sortedVals[7] = { -9, -2, 0, 3, 6, 11, 14 };
  // have one "in front", "in middle", and "in back" of vals
  int numNoVal = 3, valsNotInArray[3] = { -99, 1, 99 };

  // test adding an array in multiple ways

  // test add unsorted
  // add to empty
  fc_initSortedIntArray(&sia);
  count = fc_addIntArrayToSortedIntArray(&sia, numVal, unsortedVals, 0);
  fail_unless(count == numVal, "failed to add unsorted array");
  fail_unless(sia.numVal == numVal, "mismatch of numVal");
  fail_unless(!memcmp(sia.vals, sortedVals, numVal*sizeof(int)),
	      "sorted wrong");
  // read the same ones
  count = fc_addIntArrayToSortedIntArray(&sia, numVal, unsortedVals, 0);
  fail_unless(count == 0, "should not add any already there");
  fail_unless(sia.numVal == numVal, "should not change numVal");
  fail_unless(!memcmp(sia.vals, sortedVals, numVal*sizeof(int)),
	      "should not change vals");
  // delete off first two & add back
  fc_deleteIntFromSortedIntArray(&sia, sia.vals[0]);
  fc_deleteIntFromSortedIntArray(&sia, sia.vals[0]);
  count = fc_addIntArrayToSortedIntArray(&sia, 2, sortedVals, 0);
  fail_unless(count == 2, "should add both");
  fail_unless(sia.numVal == numVal, "mismatch of numVal");
  fail_unless(!memcmp(sia.vals, sortedVals, numVal*sizeof(int)),
	      "sorted wrong");
  // delete off last two & add back
  fc_deleteIntFromSortedIntArray(&sia, sia.vals[sia.numVal-1]);
  fc_deleteIntFromSortedIntArray(&sia, sia.vals[sia.numVal-1]);
  count = fc_addIntArrayToSortedIntArray(&sia, 2, sortedVals+(numVal-2), 0);
  fail_unless(count == 2, "should add both");
  fail_unless(sia.numVal == numVal, "mismatch of numVal");
  fail_unless(!memcmp(sia.vals, sortedVals, numVal*sizeof(int)),
	      "sorted wrong");
  // delete off two in middle
  fc_deleteIntFromSortedIntArray(&sia, sia.vals[1]);
  fc_deleteIntFromSortedIntArray(&sia, sia.vals[1]);
  count = fc_addIntArrayToSortedIntArray(&sia, 2, sortedVals+1, 0);
  fail_unless(count == 2, "should add both");
  fail_unless(sia.numVal == numVal, "mismatch of numVal");
  fail_unless(!memcmp(sia.vals, sortedVals, numVal*sizeof(int)),
	      "sorted wrong");
  // delete off first, last, & middle & add all back
  fc_deleteIntFromSortedIntArray(&sia, sia.vals[0]);
  fc_deleteIntFromSortedIntArray(&sia, sia.vals[1]);
  fc_deleteIntFromSortedIntArray(&sia, sia.vals[sia.numVal-1]);
  count = fc_addIntArrayToSortedIntArray(&sia, numVal, unsortedVals, 0);
  fail_unless(count == 3, "should add 3");
  fail_unless(sia.numVal == numVal, "mismatch of numVal");
  fail_unless(!memcmp(sia.vals, sortedVals, numVal*sizeof(int)),
	      "sorted wrong");
  fc_freeSortedIntArray(&sia);

  // test add unsorted w/ duplicates
  // add to empty
  fc_initSortedIntArray(&sia);
  count = fc_addIntArrayToSortedIntArray(&sia, numVal2, unsortedVals2, 0);
  fail_unless(count == numVal, "failed to add unsorted array with duplicates");
  fail_unless(sia.numVal == numVal, "mismatch of numVal");
  for (i = 0; i < numVal; i++) 
    fail_unless(sia.vals[i] == sortedVals[i], "sorted wrong");
  // delete off first, last, & middle & add all back
  fc_deleteIntFromSortedIntArray(&sia, sia.vals[0]);
  fc_deleteIntFromSortedIntArray(&sia, sia.vals[1]);
  fc_deleteIntFromSortedIntArray(&sia, sia.vals[sia.numVal-1]);
  count = fc_addIntArrayToSortedIntArray(&sia, numVal2, unsortedVals2, 0);
  fail_unless(count == 3, "should add 3");
  fail_unless(sia.numVal == numVal, "mismatch of numVal");
  fail_unless(!memcmp(sia.vals, sortedVals, numVal*sizeof(int)),
	      "mismatch of vals");
  fc_freeSortedIntArray(&sia);

  // test add sorted
  // add to empty
  fc_initSortedIntArray(&sia);
  count = fc_addIntArrayToSortedIntArray(&sia, numVal, sortedVals, 1);
  fail_unless(count == numVal, "failed to add sorted array");
  fail_unless(sia.numVal == numVal, "mismatch of numVal");
  for (i = 0; i < numVal; i++) 
    fail_unless(sia.vals[i] == sortedVals[i], "sorted wrong");
  // delete off first, last, & middle & add all back
  fc_deleteIntFromSortedIntArray(&sia, sia.vals[0]);
  fc_deleteIntFromSortedIntArray(&sia, sia.vals[1]);
  fc_deleteIntFromSortedIntArray(&sia, sia.vals[sia.numVal-1]);
  count = fc_addIntArrayToSortedIntArray(&sia, numVal, sortedVals, 1);
  fail_unless(count == 3, "should add 3");
  fail_unless(sia.numVal == numVal, "mismatch of numVal");
  fail_unless(!memcmp(sia.vals, sortedVals, numVal*sizeof(int)),
	      "mismatch of vals");
  // delay cleanup until after bad args test

  // o.k. to add nothing
  count = fc_addIntArrayToSortedIntArray(&sia, 0, NULL, 0);
  fail_unless(count == 0, "should be o.k. to add null array");
  fail_unless(sia.numVal == numVal, "mismatch of numVal");
  fail_unless(!memcmp(sia.vals, sortedVals, numVal*sizeof(int)),
	      "mismatch of vals");

  // test bad args
  // fc_addIntArrayToSortedIntArray
  count = fc_addIntArrayToSortedIntArray(NULL, numNoVal, valsNotInArray, 0);
  fail_unless(count < FC_SUCCESS, "should fail for null input");
  count = fc_addIntArrayToSortedIntArray(&sia, -1, valsNotInArray, 0);
  fail_unless(count < FC_SUCCESS, "should fail for null input");
  count = fc_addIntArrayToSortedIntArray(&sia, numNoVal, NULL, 0);
  fail_unless(count < FC_SUCCESS, "should fail for null input");

  // cleanup
  fc_freeSortedIntArray(&sia);
}
END_TEST

START_TEST(sia_pop)
{
  FC_ReturnCode rc;
  int count;
  FC_SortedIntArray sia, emptySia;
  int numVal = 7, sortedVals[7] = { -9, -2, 0, 3, 6, 11, 14 };
  int temp_val;
  int temp_numVal, *temp_vals; 

  // setup
  rc = fc_initSortedIntArray(&sia);
  rc = fc_initSortedIntArray(&emptySia);
  fail_unless(rc == FC_SUCCESS, "abort: failed to init sia");
  count = fc_addIntArrayToSortedIntArray(&sia, numVal, sortedVals, 1);
  fail_unless(count == numVal, "abort: failed to add ints to sia");

  // query front & back of empty sia (it's allowable, but not recommended)
  rc = fc_getSortedIntArrayFront(&emptySia, &temp_val);
  fail_unless(rc == FC_SUCCESS, "failed to get front of empty");
  fail_unless(temp_val == -1, "mismatch of empty front");
  rc = fc_getSortedIntArrayBack(&emptySia, &temp_val);
  fail_unless(rc == FC_SUCCESS, "failed to get back of empty");
  fail_unless(temp_val == -1, "mismatch of empty back");

  // pop front & back of empty sia
  count = fc_popSortedIntArrayFront(&emptySia);
  fail_unless(count == 0, "pop of empty = 0");
  fail_unless(emptySia.numVal == 0, "should still be empty");
  count = fc_popSortedIntArrayBack(&emptySia);
  fail_unless(count == 0, "pop of empty = 0");
  fail_unless(emptySia.numVal == 0, "should still be empty");

  // query front & back of nonempty sia
  rc = fc_getSortedIntArrayFront(&sia, &temp_val);
  fail_unless(rc == FC_SUCCESS, "failed to get front");
  fail_unless(temp_val == sortedVals[0], "mismatch of front");
  rc = fc_getSortedIntArrayBack(&sia, &temp_val);
  fail_unless(rc == FC_SUCCESS, "failed to get front");
  fail_unless(temp_val == sortedVals[numVal-1], "mismatch of front");

  // pop front & back of nonempty sia
  count = fc_popSortedIntArrayFront(&sia);
  fail_unless(count == 1, "pop of non empty should be 1");
  fail_unless(sia.numVal = numVal-1, "pop should decrease numVal");
  count = fc_popSortedIntArrayBack(&sia);
  fail_unless(!memcmp(sia.vals, sortedVals+1, sia.numVal*sizeof(int)),
	      "mismatch of vals");
  fail_unless(count == 1, "pop of non empty should be 1");
  fail_unless(sia.numVal = numVal-2, "pop should decrease numVal");
  fail_unless(!memcmp(sia.vals, sortedVals+1, sia.numVal*sizeof(int)),
	      "mismatch of vals");
  
  //Check and make sure insert order is still correct
  rc = fc_getSortedIntArrayValuesInInsertOrder(&sia, &temp_numVal, &temp_vals);
  fail_unless(rc == FC_SUCCESS, "failed to get vals in insert order");
  fail_unless(temp_numVal == numVal-2, "Wrong number of vals");
  fail_unless(!memcmp(&temp_vals[0], &sortedVals[1], (numVal-2)*sizeof(int)), "unsorted vals don't line up");
  free(temp_vals);

  // cleanup
  fc_freeSortedIntArray(&sia);
}
END_TEST

START_TEST(sia_insert_order)
{
  FC_ReturnCode rc;
  int count;
  FC_SortedIntArray sia, emptySia;
  int numVal = 7, sortedVals[7] = { -9, -2, 0, 3, 6, 11, 14 };
  int temp_val;
  int temp_numVal, *temp_vals; 

  int i,j;

  typedef struct {
    int num;
    char op;
    int vals[10];
    int exp_num; //Expected number of sorted vals
    int exp_sort[20]; //Expected sorted vals
    int exp_insrt[20]; //Expected vals in insert order
  } commands_t;


  commands_t cmd[] = 
    { { 3, '+', {  0,2,4 },     3, {0,2,4},                 {0,2,4} },
      { 3, '+', {  1,3,5 },     6, {0,1,2,3,4,5},           {0,2,4,1,3,5} }, 
      { 3, '+', { 10,6,8 },     9, {0,1,2,3,4,5,6,8,10},    {0,2,4,1,3,5,10,6,8} },
      { 3, '+', {  9,6,7 },    11, {0,1,2,3,4,5,6,7,8,9,10},{0,2,4,1,3,5,10,6,8,9,7} }, //note: duplicate 6
      { 3, '-', {  2,0,7 },     8, {1,3,4,5,6,8,9,10},      {4,1,3,5,10,6,8,9,7} },
      { 5, '-', { 1,5,8,10,9 }, 3, {3,4,6},                 {4,3,6} },
      { 5, '+', { 1,8,8,7,1 },  6, {1,3,4,6,7,8},           {4,3,6,1,8,7} },
      { 0, 'x', {},             0, {}, {}}
    };

  rc = fc_initSortedIntArray(&sia);
  fail_unless(rc == FC_SUCCESS, "abort: failed to init sia");

  for(i=0; cmd[i].num; i++){
    switch(cmd[i].op){
    case '+':
      rc = fc_addIntArrayToSortedIntArray(&sia, cmd[i].num, cmd[i].vals,0);
      break;
    case '-':
      for(j=0; j<cmd[i].num; j++)
	rc = fc_deleteIntFromSortedIntArray(&sia, cmd[i].vals[j]);
      break; 
    default:
      fail_unless(0, "unexpected state in case statement");
    }

    //Check sorted version
    rc = fc_getSortedIntArrayValues(&sia, &temp_numVal, &temp_vals);
    fail_unless(rc == FC_SUCCESS, "failed to get values from sorted int array");
    fail_unless(temp_numVal == cmd[i].exp_num, "Got wrong number of values");
    fail_unless(!memcmp(temp_vals, cmd[i].exp_sort, cmd[i].exp_num*sizeof(int)), "Got wrong sorted values");
    free(temp_vals);

    //Check unsorted version 
    rc = fc_getSortedIntArrayValuesInInsertOrder(&sia, &temp_numVal, &temp_vals);
    fail_unless(rc == FC_SUCCESS, "failed to get values from sorted int array w/ insert order");
    fail_unless(temp_numVal == cmd[i].exp_num, "Got wrong number of values for insert order");
    fail_unless(!memcmp(temp_vals, cmd[i].exp_insrt, cmd[i].exp_num*sizeof(int)), "Got wrong values for insert order");
    free(temp_vals);
  }

  fc_freeSortedIntArray(&sia);
  

}
END_TEST



// sorted blob arrays

START_TEST(sba_helpers) 
{
  FC_ReturnCode rc;
  FC_SortedBlobArray sba;
  int i, flag, idx, midID = 3;
  int numVal = 7, sortedVals[7] = { -9, -2, 0, 3, 6, 11, 14 };
  int val = 99;
  int temp_length;
  void** temp_array;

  // -- test _fc_expandBlobArray() (blob = int*)

  // length can't be negative
  temp_length = -1;
  temp_array = (void**)1;
  rc = _fc_expandBlobArray(&temp_length, &temp_array);
  fail_unless(rc != FC_SUCCESS, "should fail if length < 0");
  fail_unless(temp_array == (void*)1, "should not change array");
  // if 0, get length of 1
  temp_length = 0;
  rc = _fc_expandBlobArray(&temp_length, &temp_array);
  fail_unless(rc == FC_SUCCESS, "should not fial if length 0");
  fail_unless(temp_length == 1, "should be of size 1");
  fail_unless(temp_array != NULL && temp_array != (void**)1,
	      "should have changed array");
  temp_array[0] = (void*)1;
  rc = _fc_expandBlobArray(&temp_length, &temp_array);
  fail_unless(rc == FC_SUCCESS, "should not fail if length > 0");
  fail_unless(temp_length == 2, "should be of size 2");
  fail_unless(temp_array[0] == (void*)1, 
	      "first val should be the same");
  temp_array[1] = (void*)2;
  rc = _fc_expandBlobArray(&temp_length, &temp_array);
  fail_unless(rc == FC_SUCCESS, "should not fail if length > 0");
  fail_unless(temp_length == 4, "should be of size 4");
  fail_unless(temp_array[0] == (void*)1, "1st val should be the same");
  fail_unless(temp_array[1] == (void*)2, "2nd val should be the same");
  temp_array[2] = (void*)3;
  temp_array[3] = (void*)4;
  rc = _fc_expandBlobArray(&temp_length, &temp_array);
  fail_unless(rc == FC_SUCCESS, "should not fail if length > 0");
  fail_unless(temp_length == 8, "should be of size 8");
  fail_unless(temp_array[0] == (void*)1, "1st val should be the same");
  fail_unless(temp_array[1] == (void*)2, "2nd val should be the same");
  fail_unless(temp_array[2] == (void*)3, "3rd val should be the same");
  fail_unless(temp_array[3] == (void*)4, "4th val should be the same");
  free(temp_array);

  // -- test _fc_lookupBlobInSortedBlobArray
 
  // test lookup = binary search algrithm (hack an sba)
  sba.numBlob = numVal;
  sba.maxNumBlob = numVal;
  sba.blobs = malloc(numVal*sizeof(void*));
  fail_unless(sba.blobs != NULL, "abort: memory allocation failure");
  for (i = 0; i < numVal; i++)
    sba.blobs[i] = &sortedVals[i];
  for (i = 0; i < numVal; i++) {
    rc = _fc_lookupBlobInSortedBlobArray(&sba, &sortedVals[i], fc_intCompare,
					 &flag, &idx);
    fail_unless(rc == FC_SUCCESS, "failed to get index");
    fail_unless(flag == 1 && idx == i, "mismatch of found index");
  }
  for (i = 0; i < numVal; i++) {
    val = sortedVals[i]-1;
    rc = _fc_lookupBlobInSortedBlobArray(&sba, &val, fc_intCompare,
					 &flag, &idx);
    fail_unless(rc == FC_SUCCESS, "failed to get index");
    fail_unless(flag == 0 && idx == i, "mismatch of not found index");
  }
  for (i = 0; i < numVal; i++) {
    val = sortedVals[i]+1;
    rc = _fc_lookupBlobInSortedBlobArray(&sba, &val, fc_intCompare,
					 &flag, &idx);
    fail_unless(rc == FC_SUCCESS, "failed to get index");
    fail_unless(flag == 0 && idx == i+1, "mismatch of not found index");
  }
  val = 99;

  // test bad args
  rc = _fc_lookupBlobInSortedBlobArray(NULL, &val, fc_intCompare, &flag, &idx);
  fail_unless(rc != FC_SUCCESS, "should fail for null input");
  fail_unless(flag == -1 && idx == -1, "fail should return nulls");
  rc = _fc_lookupBlobInSortedBlobArray(&sba, NULL, fc_intCompare, &flag, &idx);
  fail_unless(rc != FC_SUCCESS, "should fail for null input");
  fail_unless(flag == -1 && idx == -1, "fail should return nulls");
  rc = _fc_lookupBlobInSortedBlobArray(&sba, &val, NULL, &flag, &idx);
  fail_unless(rc != FC_SUCCESS, "should fail for null input");
  fail_unless(flag == -1 && idx == -1, "fail should return nulls");
  rc = _fc_lookupBlobInSortedBlobArray(&sba, &val, fc_intCompare, NULL, &idx);
  fail_unless(rc != FC_SUCCESS, "should fail for null input");
  fail_unless(idx == -1, "fail should return nulls");
  rc = _fc_lookupBlobInSortedBlobArray(&sba, &val, fc_intCompare, &flag, NULL);
  fail_unless(rc != FC_SUCCESS, "should fail for null input");
  fail_unless(flag == -1, "fail should return nulls");

  // -- test _fc_addEntryToSortedBlobArray & _fc_deleteEntryFromSortedBlobArray

  // test add/delete from front
  rc = _fc_addEntryToSortedBlobArray(&sba, 0, &val);
  fail_unless(rc == FC_SUCCESS, "failed to add");
  fail_unless(sba.numBlob == numVal+1, "mismatch of numBlob");
  fail_unless(sba.blobs[0] == &val, "mismatch of blob");
  for (i = 0; i < numVal; i++)
    fail_unless(sba.blobs[i+1] == &sortedVals[i], "mismathc of blobs");
  rc = _fc_deleteEntryFromSortedBlobArray(&sba, 0);
  fail_unless(rc == FC_SUCCESS, "failed to add");
  fail_unless(sba.numBlob == numVal, "mismatch of numBlob");
  for (i = 0; i < numVal; i++)
    fail_unless(sba.blobs[i] == &sortedVals[i], "mismathc of blobs");
  
  // test add/delete from middle
  rc = _fc_addEntryToSortedBlobArray(&sba, midID, &val);
  fail_unless(rc == FC_SUCCESS, "failed to add");
  fail_unless(sba.numBlob == numVal+1, "mismatch of blob");
  fail_unless(sba.blobs[midID] == &val, "mismatch of blob");
  for (i = 0; i < midID; i++)
    fail_unless(sba.blobs[i] == &sortedVals[i], "mismatch of blobs");
  for (i = midID; i < numVal; i++)
    fail_unless(sba.blobs[i+1] == &sortedVals[i], "mismath of blobs");
  rc = _fc_deleteEntryFromSortedBlobArray(&sba, midID);
  fail_unless(rc == FC_SUCCESS, "failed to add");
  fail_unless(sba.numBlob == numVal, "mismatch of blob");
  for (i = 0; i < numVal; i++)
    fail_unless(sba.blobs[i] == &sortedVals[i], "mismatch of blobs");

  // test add/delete from back
  rc = _fc_addEntryToSortedBlobArray(&sba, numVal, &val);
  fail_unless(rc == FC_SUCCESS, "failed to add");
  fail_unless(sba.numBlob == numVal+1, "mismatch of blob");
  fail_unless(sba.blobs[numVal] == &val, "mismatch of blob");
  for (i = 0; i < numVal; i++)
    fail_unless(sba.blobs[i] == &sortedVals[i], "mismatch of blobs");
  rc = _fc_deleteEntryFromSortedBlobArray(&sba, numVal);
  fail_unless(rc == FC_SUCCESS, "failed to add");
  fail_unless(sba.numBlob == numVal, "mismatch of blob");
  for (i = 0; i < numVal; i++)
    fail_unless(sba.blobs[i] == &sortedVals[i], "mismatch of blobs");

  // test bad args
  rc = _fc_addEntryToSortedBlobArray(NULL, 0, &val);
  fail_unless(rc != FC_SUCCESS, "should fail w/ null args");
  rc = _fc_addEntryToSortedBlobArray(&sba, 0, NULL);
  fail_unless(rc != FC_SUCCESS, "should fail w/ null args");
  rc = _fc_deleteEntryFromSortedBlobArray(NULL, 0);
  fail_unless(rc != FC_SUCCESS, "should fail w/ null args");
  rc = _fc_addEntryToSortedBlobArray(&sba, -1, &val);
  fail_unless(rc != FC_SUCCESS, "should fail w/ negative index");
  rc = _fc_deleteEntryFromSortedBlobArray(&sba, -1);
  fail_unless(rc != FC_SUCCESS, "should fail w/ negative index");
  // not "too big index" is not the same for the two routines
  rc = _fc_addEntryToSortedBlobArray(&sba, numVal+1, &val);
  fail_unless(rc != FC_SUCCESS, "should fail w/ index too big");
  rc = _fc_deleteEntryFromSortedBlobArray(&sba, numVal);
  fail_unless(rc != FC_SUCCESS, "should fail w/ index too big");

  // cleanup
  free(sba.blobs);

  // test border cases: empty arrays w & w/o blobs=NULL
  // & also that growing works
  sba.numBlob = 0;
  sba.maxNumBlob = 0;
  sba.blobs = NULL;
  rc = _fc_deleteEntryFromSortedBlobArray(&sba, 0);
  fail_unless(rc != FC_SUCCESS, "should fail to delete from empty");
  rc = _fc_addEntryToSortedBlobArray(&sba, 0, &val);
  fail_unless(rc == FC_SUCCESS, "should add to empty");
  fail_unless(sba.numBlob == 1, "mismatch of numBlob");
  fail_unless(sba.blobs[0] == &val, "mismatch of blobs");
  for (i = 0; i < 10; i++) { // add a bunch to test growing
    rc = _fc_addEntryToSortedBlobArray(&sba, 0, &val);
    fail_unless(rc == FC_SUCCESS, "should add to empty");
  }
  fail_unless(sba.numBlob == 11, "mismatch of numBlob");
  sba.numBlob = 0;
  sba.blobs[0] = &sortedVals[0];
  rc = _fc_deleteEntryFromSortedBlobArray(&sba, 0);
  fail_unless(rc != FC_SUCCESS, "should fail to delete from empty");
  rc = _fc_addEntryToSortedBlobArray(&sba, 0, &val);
  fail_unless(rc == FC_SUCCESS, "should add to empty");
  fail_unless(sba.numBlob == 1, "mismatch of numBlob");
  fail_unless(sba.blobs[0] == &val, "mismatch of blobs");

  // cleanup again
  free(sba.blobs);
}
END_TEST

START_TEST(sba_valid) 
{
  FC_SortedBlobArray sba;
  int flag, i;
  int numVal = 7, sortedVals[7] = { -9, -2, 0, 3, 6, 11, 14 };
  void* blobs[7];

  // hack the sba & test
  for (i = 0; i < numVal; i++)
    blobs[i] = &sortedVals[i];

  // valid cases
  // all zeros 
  sba.numBlob = 0;
  sba.maxNumBlob = 0;
  sba.blobs = NULL;
  flag = fc_isSortedBlobArrayValid(&sba);
  fail_unless(flag == 1, "should return true");
  // no entries but has space
  sba.numBlob = 0;
  sba.maxNumBlob = numVal;
  sba.blobs = blobs;
  flag = fc_isSortedBlobArrayValid(&sba);
  fail_unless(flag == 1, "should return true");
  // has same number of entries and space
  sba.numBlob = numVal;
  sba.maxNumBlob = numVal;
  sba.blobs = blobs;
  flag = fc_isSortedBlobArrayValid(&sba);
  fail_unless(flag == 1, "should return true");
  // has fewer entires than space
  sba.numBlob = 1;
  sba.maxNumBlob = numVal;
  sba.blobs = blobs;
  flag = fc_isSortedBlobArrayValid(&sba);
  
  // invalid cases
  // numBlob is negative
  sba.numBlob = -1;
  sba.maxNumBlob = 0;
  sba.blobs = NULL;
  flag = fc_isSortedBlobArrayValid(&sba);
  fail_unless(flag == 0, "should return false");
  // maxNumBlob is negative
  sba.numBlob = 0;
  sba.maxNumBlob = -1;
  sba.blobs = NULL;
  flag = fc_isSortedBlobArrayValid(&sba);
  fail_unless(flag == 0, "should return false");
  // space is zero, but there appears to be something there
  sba.numBlob = 0;
  sba.maxNumBlob = 0;
  sba.blobs = blobs;
  flag = fc_isSortedBlobArrayValid(&sba);
  fail_unless(flag == 0, "should return false");
  // there should be something there, but there is not
  sba.numBlob = numVal;
  sba.maxNumBlob = numVal;
  sba.blobs = NULL;
  flag = fc_isSortedBlobArrayValid(&sba);
  fail_unless(flag == 0, "should return false");
  // there are more entries than space
  sba.numBlob = numVal+1;
  sba.maxNumBlob = numVal;
  sba.blobs = blobs;
  flag = fc_isSortedBlobArrayValid(&sba);
  fail_unless(flag == 0, "should return false");
}
END_TEST

START_TEST(sba_basic) 
{
  FC_ReturnCode rc;
  FC_SortedBlobArray sba;
  int i, count;
  // have enough vals to cause multiple grows, but do not
  // have a multiple of 2 so can test squeeze
  int numVal = 7, unsortedVals[7] = { 3, 14, -9, 0, 6, 11, -2 };
  // The same unsorted vals, but they look like different blobs
  int unsortedVals2[7] = { 3, 14, -9, 0, 6, 11, -2 };
  // note the vals are all at least 2 apart so I can easily, rigoursly
  // test the binary search algorithm by using vals[i]+/-1.
  int sortedVals[7] = { -9, -2, 0, 3, 6, 11, 14 };
  // have one "in front", "in middle", and "in back" of vals
  int numNoVal = 3, valsNotInArray[3] = { -99, 1, 99 };
  int temp_numBlob;
  void* temp_blob, **temp_blobs;
 
  // test fc_initSortedBlobArray
  sba.numBlob = sba.maxNumBlob = -2;
  sba.blobs = (void**)1;
  rc = fc_initSortedBlobArray(&sba);
  fail_unless(rc == FC_SUCCESS, "bad return");
  fail_unless(sba.numBlob == 0, "bad val from init");
  fail_unless(sba.maxNumBlob == 0, "bad val from init");
  fail_unless(sba.blobs == NULL, "bad val from init");
 
  // test query & delete of empty & NULL sba
  count = fc_isBlobInSortedBlobArray(&sba, &sortedVals[0], fc_intCompare,
				     NULL);
  fail_unless(count == 0, "should not find vals in empty sba");
  rc = fc_getSortedBlobArrayNumBlob(&sba, &temp_numBlob);
  fail_unless(rc == FC_SUCCESS, "should find numVal");
  fail_unless(temp_numBlob == 0, "numVal should be zero for empty");
  rc = fc_getSortedBlobArrayBlobs(&sba, &temp_numBlob, &temp_blobs);
  fail_unless(rc == FC_SUCCESS, "should find values");
  fail_unless(temp_numBlob == 0, "numVal should be zero for empty");
  fail_unless(temp_blobs == NULL, "vals should be zero for empty");
  count = fc_deleteBlobFromSortedBlobArray(&sba, &sortedVals[0], fc_intCompare,
					   NULL);
  fail_unless(count == 0, "should not delete vals from empty sba");

  // test add single blob
  for (i = 0; i < numVal; i++) {
    count = fc_addBlobToSortedBlobArray(&sba, &unsortedVals[i], fc_intCompare);
    fail_unless(count == 1, "create should create");
    fail_unless(sba.numBlob == i+1, "sorted wrong");
  }
  fail_unless(sba.numBlob == numVal, "bad add");
  for (i = 0; i < numVal; i++)
    fail_unless(((int*)sba.blobs[i])[0] == sortedVals[i], "sorted wrong");
 
  // test query
  for (i = 0; i < numVal; i++) {
    count = fc_isBlobInSortedBlobArray(&sba, &unsortedVals[i], fc_intCompare,
				       &temp_blob);
    fail_unless(count == 1, "should find vals");
    fail_unless(temp_blob == &unsortedVals[i], "mismatch of found blob");
  }
  for (i = 0; i < numVal; i++) {
    count = fc_isBlobInSortedBlobArray(&sba, &unsortedVals2[i], fc_intCompare,
				       &temp_blob);
    fail_unless(count == 1, "should find vals");
    fail_unless(((int*)temp_blob)[0] == unsortedVals2[i], 
		"mismatch of found blob");
  }
  for (i = 0; i < numNoVal; i++) {
    count = fc_isBlobInSortedBlobArray(&sba, &valsNotInArray[i], fc_intCompare,
				       &temp_blob);
    fail_unless(count == 0, "should not find not there vals");
    fail_unless(temp_blob == NULL, "should be null if no find blob");
  }
  rc = fc_getSortedBlobArrayNumBlob(&sba, &temp_numBlob);
  fail_unless(rc == FC_SUCCESS, "should find numVal");
  fail_unless(temp_numBlob == numVal, "numVal should agree");
  rc = fc_getSortedBlobArrayBlobs(&sba, &temp_numBlob, &temp_blobs);
  fail_unless(rc == FC_SUCCESS, "should find values");
  fail_unless(temp_numBlob == numVal, "numVal mismatch");
  fail_unless(!memcmp(temp_blobs, sba.blobs, numVal*sizeof(int)),
        "values mismatch");
  free(temp_blobs);
  fail_unless(sba.numBlob == numVal, "queries altered contents");
  for (i = 0; i < numVal; i++)
    fail_unless(((int*)sba.blobs[i])[0] == sortedVals[i], "sorted wrong");

  // test try to add again (add same blob content, but not same blob)
  for (i = 0; i < numVal; i++) {
    count = fc_addBlobToSortedBlobArray(&sba, &unsortedVals2[i], fc_intCompare);
    fail_unless(count == 0, "create doesn't create if already there");
  }
  fail_unless(sba.numBlob == numVal, "recreating went wrong");
  for (i = 0; i < numVal; i++)
    fail_unless(((int*)sba.blobs[i])[0] == sortedVals[i],
		"recreating went wrong");
 
  // test delete single blob
  for (i = 0; i < numNoVal; i++) {
    count = fc_deleteBlobFromSortedBlobArray(&sba, &valsNotInArray[i],
					     fc_intCompare, &temp_blob);
    fail_unless(count == 0, "should not delete not there vals");
    fail_unless(temp_blob == 0, "if not there, should return null");
  }
  fail_unless(sba.numBlob == numVal, "deleting nothing went wrong");
  for (i = 0; i < numVal; i++) {
    count = fc_deleteBlobFromSortedBlobArray(&sba, &unsortedVals2[i],
					     fc_intCompare, &temp_blob);
    fail_unless(count == 1, "should delete vals");
    fail_unless(sba.numBlob == numVal-i-1, "didn't delete vals");
    fail_unless(((int*)temp_blob)[0] == unsortedVals2[i],
		"mismtach of deleted blob");
  }
  fail_unless(sba.numBlob == 0, "should be empty");
  fail_unless(sba.maxNumBlob != 0, "but still have space");
  fail_unless(sba.blobs != NULL, "but still have space");

  // test query and delete of empty, but not NULL sba
  count = fc_deleteBlobFromSortedBlobArray(&sba, &sortedVals[0], fc_intCompare,
					   NULL);
  fail_unless(count == 0, "should not delete vals from empty sba");
  count = fc_isBlobInSortedBlobArray(&sba, &sortedVals[0], fc_intCompare,
				     NULL);
  fail_unless(count == 0, "should not find vals in empty sba");
  rc = fc_getSortedBlobArrayNumBlob(&sba, &temp_numBlob);
  fail_unless(rc == FC_SUCCESS, "should find numVal");
  fail_unless(temp_numBlob == 0, "empty sba should have no values");
  rc = fc_getSortedBlobArrayBlobs(&sba, &temp_numBlob, &temp_blobs);
  fail_unless(rc == FC_SUCCESS, "should find values");
  fail_unless(temp_numBlob == 0, "numVal mismatch");
  fail_unless(temp_blobs == NULL, "values mismatch");
  
  // setup for bad args test
  for (i = 0; i < numVal; i++)
    fc_addBlobToSortedBlobArray(&sba, &unsortedVals[i], fc_intCompare);
  fail_unless(sba.numBlob == numVal, "abort: bad add");
  for (i = 0; i < numVal; i++)
    fail_unless(((int*)sba.blobs[i])[0] == sortedVals[i],
		"abort: bad add");
  
  // test bad args
  // fc_initSortedBlobArray
  rc = fc_initSortedBlobArray(NULL);
  fail_unless(rc != FC_SUCCESS, "should fail for null input");
  // fc_getSortedBlobArrayNumBlob
  rc = fc_getSortedBlobArrayNumBlob(NULL, &temp_numBlob);
  fail_unless(rc != FC_SUCCESS, "should fail for null input");
  fail_unless(temp_numBlob == -1, "fail should return nulls");
  rc = fc_getSortedBlobArrayNumBlob(&sba, NULL);
  fail_unless(rc != FC_SUCCESS, "should fail for null input");
  // fc_getSortedBlobArrayBlobs
  rc = fc_getSortedBlobArrayBlobs(NULL, &temp_numBlob, &temp_blobs);
  fail_unless(rc != FC_SUCCESS, "should fail for null input");
  fail_unless(temp_numBlob == -1, "fail should return nulls");
  fail_unless(temp_blobs == NULL, "fail should return nulls");
  rc = fc_getSortedBlobArrayBlobs(&sba, NULL, &temp_blobs);
  fail_unless(rc != FC_SUCCESS, "should fail for null input");
  fail_unless(temp_blobs == NULL, "fail should return nulls");
  rc = fc_getSortedBlobArrayBlobs(&sba, &temp_numBlob, NULL);
  fail_unless(rc != FC_SUCCESS, "should fail for null input");
  fail_unless(temp_numBlob == -1, "fail should return nulls");
  // fc_addIntToSortedBlobArray
  count = fc_addBlobToSortedBlobArray(NULL, &sortedVals[0], fc_intCompare);
  fail_unless(count < FC_SUCCESS, "should fail for null input");
  count = fc_addBlobToSortedBlobArray(&sba, NULL, fc_intCompare);
  fail_unless(count < FC_SUCCESS, "should fail for null input");
  count = fc_addBlobToSortedBlobArray(&sba, &sortedVals[0], NULL);
  fail_unless(count < FC_SUCCESS, "should fail for null input");
  // fc_isIntInSortedBlobArray
  count = fc_isBlobInSortedBlobArray(NULL, &sortedVals[0], fc_intCompare,
				    &temp_blob);
  fail_unless(count < FC_SUCCESS, "should fail for null input");
  fail_unless(temp_blob == NULL, "fail should return null");
  count = fc_isBlobInSortedBlobArray(&sba, NULL, fc_intCompare, &temp_blob);
  fail_unless(count < FC_SUCCESS, "should fail for null input");
  fail_unless(temp_blob == NULL, "fail should return null");
  count = fc_isBlobInSortedBlobArray(&sba, &sortedVals[0], NULL, &temp_blob);
  fail_unless(count < FC_SUCCESS, "should fail for null input");
  fail_unless(temp_blob == NULL, "fail should return null");
  // fc_deleteIntFromSortedBlobArray
  count = fc_deleteBlobFromSortedBlobArray(NULL, &sortedVals[0], fc_intCompare,
					  &temp_blob);
  fail_unless(count < FC_SUCCESS, "should fail for null input");
  fail_unless(temp_blob == NULL, "fail should return null");
  count = fc_deleteBlobFromSortedBlobArray(&sba, NULL, fc_intCompare,
					   &temp_blob); 
  fail_unless(count < FC_SUCCESS, "should fail for null input");
  fail_unless(temp_blob == NULL, "fail should return null");
  count = fc_deleteBlobFromSortedBlobArray(&sba, &sortedVals[0], NULL,
					   &temp_blob);
  fail_unless(count < FC_SUCCESS, "should fail for null input");
  fail_unless(temp_blob == NULL, "fail should return null");
  // fc_printSortedBlobArray
  //rc = fc_printSortedBlobArray(NULL);
  //fail_unless(rc != FC_SUCCESS, "should fail for null input");
  // fc_freeSortedBlobArray
  fc_freeSortedBlobArray(NULL);

  // test freeing 
  fc_freeSortedBlobArray(&sba);
  fail_unless(sba.numBlob == 0, "bad val from free");
  fail_unless(sba.maxNumBlob == 0, "bad val from free");
  fail_unless(sba.blobs == NULL, "bad val from free");
}
END_TEST

START_TEST(sba_squeeze) 
{
  FC_ReturnCode rc;
  int i;
  FC_SortedBlobArray sba;
  // do not have a multiple of 2 so can test squeeze
  int numVal = 7, sortedVals[7] = { -9, -2, 0, 3, 6, 11, 14 };

  // test squeeze of empty & NULL
  rc = fc_initSortedBlobArray(&sba);
  fail_unless(rc == FC_SUCCESS, "abort: failed to init sba");
  rc = fc_squeezeSortedBlobArray(&sba);
  fail_unless(rc == FC_SUCCESS, "failed to squeeze");
  fail_unless(sba.numBlob == 0, "mismatch of numBlob");
  fail_unless(sba.maxNumBlob == 0, "mismatch of maxNumBlob");
  fail_unless(sba.blobs == NULL, "mismatch of blobs");

  // test squeeze of non empty & extra room
  for (i = 0; i < numVal; i++) {
    rc = fc_addBlobToSortedBlobArray(&sba, &sortedVals[i], fc_intCompare);
    fail_unless(rc >= FC_SUCCESS, "abort: failed to add blobs to sba");
  }
  fail_unless(sba.maxNumBlob > sba.numBlob, "abort: test not setup correctly");
  rc = fc_squeezeSortedBlobArray(&sba);
  fail_unless(rc == FC_SUCCESS, "failed to squeeze");
  fail_unless(sba.numBlob == numVal, "mismatch of numBlob");
  fail_unless(sba.maxNumBlob == numVal, "mismatch of maxNumBlob");
  for (i = 0; i < numVal; i++)
    fail_unless(((int*)sba.blobs[i])[0] == sortedVals[i],
		"mismatch of blobs");

  // test squeeze of empty & extra room
  fail_unless(sba.maxNumBlob == sba.numBlob && sba.numBlob > 0, 
	      "abort: test not setup correctly");
  rc = fc_squeezeSortedBlobArray(&sba);
  fail_unless(rc == FC_SUCCESS, "failed to squeeze");
  fail_unless(sba.numBlob == numVal, "mismatch of numBlob");
  fail_unless(sba.maxNumBlob == numVal, "mismatch of maxNumBlob");
  for (i = 0; i < numVal; i++)
    fail_unless(((int*)sba.blobs[i])[0] == sortedVals[i],
		"mismatch of blobs");

  // test copy of empty & not NULL
  for (i = 0; i < numVal; i++)
    fc_deleteBlobFromSortedBlobArray(&sba, &sortedVals[i], fc_intCompare, NULL);
  fail_unless(sba.maxNumBlob > 0 && sba.numBlob == 0, 
	      "abort: test not setup correctly");
  rc = fc_squeezeSortedBlobArray(&sba);
  fail_unless(rc == FC_SUCCESS, "failed to squeeze");
  fail_unless(sba.numBlob == 0, "mismatch of numBlob");
  fail_unless(sba.maxNumBlob == 0, "mismatch of maxNumBlob");
  fail_unless(sba.blobs == NULL, "mismatch of blobs");

  // test bad args
  rc = fc_squeezeSortedBlobArray(NULL);
  fail_unless(rc != FC_SUCCESS, "should fail for null input");

  // cleanup not needed! (if squeeze is working)
}
END_TEST

// --------------------------------------------------------------
//              conversions 
// --------------------------------------------------------------

// test back and forth between sia and int array
START_TEST(sia_array)
{
  FC_ReturnCode rc;
  FC_SortedIntArray sia, temp_sia;
  int i;
  int numVal = 7;
  int unsortedVals[7] = { 3, 14, -9, 0, 6, 11, -2 };
  int sortedVals[7] = { -9, -2, 0, 3, 6, 11, 14 };
  int *array, *sorted_array; 
  int temp_numVal, *temp_array;

  // setup 
  fc_initSortedIntArray(&sia);
  rc = fc_addIntArrayToSortedIntArray(&sia, numVal, sortedVals, 1);
  fail_unless(rc >= FC_SUCCESS, "abort: failed to add ints to sia");
  array = malloc(numVal*sizeof(int));
  sorted_array = malloc(numVal*sizeof(int));
  fail_unless(array != NULL && sorted_array != NULL, 
	      "abort: failed to alloc space for arrays");
  for (i = 0; i < numVal; i++) {
    array[i] = unsortedVals[i];
    sorted_array[i] = sortedVals[i];
  }

  // convert sorted int array -> sia
  rc = fc_convertIntArrayToSortedIntArray(numVal, sorted_array, 1, &temp_sia);
  fail_unless(rc == FC_SUCCESS, "failed to convert int array to sia");
  fail_unless(temp_sia.numVal == numVal, "mismatch of numVal");
  fail_unless(!memcmp(temp_sia.vals, sortedVals, numVal*sizeof(int)),
	      "mismatch of vals");
  fc_freeSortedIntArray(&temp_sia);

  // convert int array -> sia
  rc = fc_convertIntArrayToSortedIntArray(numVal, array, 0, &temp_sia);
  fail_unless(rc == FC_SUCCESS, "failed to convert int array to sia");
  fail_unless(temp_sia.numVal == numVal, "mismatch of numVal");
  fail_unless(temp_sia.maxNumVal == numVal, "mismatch of maxNumVal");
  fail_unless(!memcmp(temp_sia.vals, sortedVals, numVal*sizeof(int)),
	      "mismatch of vals");
  fc_freeSortedIntArray(&temp_sia);

  // convert sia -> int array
  rc = fc_convertSortedIntArrayToIntArray(&sia, &temp_numVal, &temp_array);
  fail_unless(rc == FC_SUCCESS, "failed to convert sia to int aray");
  fail_unless(temp_numVal == numVal, "mismatch of numVal");
  fail_unless(!memcmp(temp_array, sortedVals, numVal*sizeof(int)),
	      "mismatch of vals");
  fail_unless(sia.numVal == 0, "numVal should go to zero");
  fail_unless(sia.maxNumVal == 0, "maxNumVal should go to zero");
  fail_unless(sia.vals == NULL, "vals should go to zero");
  free(temp_array);
}
END_TEST

// test back and forth between mask and array
START_TEST(array_mask)
{
  FC_ReturnCode rc;
  int i;
  int numMask = 10;
  int mask_good[10] =  { 0, 1, 0, 0, 1, 0, 0, 0, 0, 1 };
  int mask_empty[10] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
  int numArray_good = 3;
  int array_good[3] = { 1, 4, 9 };
  int bad_array1[3] = { 2, -1, 4 };
  int bad_array2[3] = { 2, 12, 4 };
  int numArray, *array, *mask;

  //--- test fc_createIntArrayFromMask

  // test bad input
  numArray = 12; array = (int*)1;
  rc = fc_createIntArrayFromMask(0, mask_good, &numArray, &array);
  fail_unless(rc != FC_SUCCESS, "should fail with bad input");
  fail_unless(numArray == -1, "failed should return null");
  fail_unless(array == 0, "failed should return null");
  numArray = 12; array = (int*)1;
  rc = fc_createIntArrayFromMask(numMask, NULL, &numArray, &array);
  fail_unless(rc != FC_SUCCESS, "should fail with bad input");
  fail_unless(numArray == -1, "failed should return null");
  fail_unless(array == 0, "failed should return null");
  numArray = 12; array = (int*)1;
  rc = fc_createIntArrayFromMask(numMask, mask_good, NULL, &array);
  fail_unless(rc != FC_SUCCESS, "should fail with bad input");
  fail_unless(array == 0, "failed should return null");
  numArray = 12; array = (int*)1;
  rc = fc_createIntArrayFromMask(numMask, mask_good, &numArray, NULL);
  fail_unless(rc != FC_SUCCESS, "should fail with bad input");
  fail_unless(numArray == -1, "failed should return null");

  // test good input - border case: empty mask
  array = (int*)1;
  rc = fc_createIntArrayFromMask(numMask, mask_empty, &numArray, &array);
  fail_unless(rc == FC_SUCCESS, "should not fail on empty mask");
  fail_unless(numArray == 0, "mismatch of numArray");
  fail_unless(array == NULL, "mismatch of array");
 
  // test good input
  rc = fc_createIntArrayFromMask(numMask, mask_good, &numArray, &array);
  fail_unless(rc == FC_SUCCESS, "should not fail");
  fail_unless(numArray == numArray_good, "mismatch of numArray");
  for (i = 0; i < numArray; i++)
    fail_unless(array[i] == array_good[i], "mismatch of ints");
  free(array);

  //--- test fc_createMaskFromIntArray

  // test bad input
  mask = (int*)1;
  rc = fc_createMaskFromIntArray(-1, array_good, numMask, &mask);
  fail_unless(rc != FC_SUCCESS, "should fail with bad input");
  fail_unless(mask == NULL, "failed should return null");
  mask = (int*)1;
  rc = fc_createMaskFromIntArray(numArray_good, NULL, numMask, &mask);
  fail_unless(rc != FC_SUCCESS, "should fail with bad input");
  fail_unless(mask == NULL, "failed should return null");
  mask = (int*)1;
  rc = fc_createMaskFromIntArray(numArray_good, array_good, 0, &mask);
  fail_unless(rc != FC_SUCCESS, "should fail with bad input");
  fail_unless(mask == NULL, "failed should return null");
  rc = fc_createMaskFromIntArray(numArray_good, array_good, numMask, NULL);
  fail_unless(rc != FC_SUCCESS, "should fail with bad input");

  // test bad input - int array members out of range
  rc = fc_createMaskFromIntArray(numArray_good, bad_array1, numMask, &mask);
  fail_unless(rc != FC_SUCCESS, "should fail with bad input");  
  fail_unless(mask == NULL, "failed should return null");
  rc = fc_createMaskFromIntArray(numArray_good, bad_array2, numMask, &mask);
  fail_unless(rc != FC_SUCCESS, "should fail with bad input");  
  fail_unless(mask == NULL, "failed should return null");

  // test good input: border case = numArray = 0
  rc = fc_createMaskFromIntArray(0, NULL, numMask, &mask);
  fail_unless(rc == FC_SUCCESS, "should not fail to create empty mask");
  for (i = 0; i < numMask; i++)
    fail_unless(mask[i] == 0, "mismatch of mask member");
  free(mask);

  // test good input: general case
  rc = fc_createMaskFromIntArray(numArray_good, array_good, numMask, &mask);
  fail_unless(rc == FC_SUCCESS, "should not fail to create empty mask");
  for (i = 0; i < numMask; i++)
    fail_unless(mask[i] == mask_good[i], "mismatch of mask member");
  free(mask);
}
END_TEST

// Populate the Suite with the tests

Suite *storage_suite(void)
{
  Suite *suite = suite_create("Storage");

  TCase *tc_sortedIntArrays = tcase_create(" - Sorted Int Arrays ");
  TCase *tc_sortedBlobArrays = tcase_create(" - Sorted Blob Arrays ");
  TCase *tc_conversions = tcase_create(" - Conversions ");
  
  suite_add_tcase(suite, tc_sortedIntArrays);
  tcase_add_test(tc_sortedIntArrays, sia_helpers);
  tcase_add_test(tc_sortedIntArrays, sia_valid);
  tcase_add_test(tc_sortedIntArrays, sia_basic);
  tcase_add_test(tc_sortedIntArrays, sia_copy);
  tcase_add_test(tc_sortedIntArrays, sia_squeeze);
  tcase_add_test(tc_sortedIntArrays, sia_multiadd);
  tcase_add_test(tc_sortedIntArrays, sia_pop);
  tcase_add_test(tc_sortedIntArrays, sia_insert_order);

  
  suite_add_tcase(suite, tc_sortedBlobArrays);
  tcase_add_test(tc_sortedBlobArrays, sba_helpers);
  tcase_add_test(tc_sortedBlobArrays, sba_valid);
  tcase_add_test(tc_sortedBlobArrays, sba_basic);
  tcase_add_test(tc_sortedBlobArrays, sba_squeeze);

  suite_add_tcase(suite, tc_conversions);
  tcase_add_test(tc_conversions, sia_array);
  tcase_add_test(tc_conversions, array_mask);
  
  return suite;
}
