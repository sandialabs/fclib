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
 * \file storage.c
 * \brief Implementation for \ref SimpleDataObjects Module.
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/modules/storage.c,v $
 * $Revision: 1.17 $ 
 * $Date: 2006/08/30 19:20:01 $
 */

// C library includes
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// fc library dependencies
#include "base.h"
#include "library.h"
#include "util.h"

// this module
#include "storage.h"
#include "storageP.h"

/**
 * \addtogroup  SimpleDataObjects
 * \brief Simple Data Objects
 *   
 * \description
 *
 *    The majority of these data objects were created to manipulate
 *    collections of ids.
 */

/**
 * \ingroup  SimpleDataObjects
 * \defgroup PrivateSimpleDataObjects (Private)
 */

/** 
 * \name Sorted Int Array Routines
 *
 * Sorted, exclusive, self-exanding int array.
 */
//-------------------------------
//@{

/**
 * \ingroup SimpleDataObjects
 * \brief Initialize a sorted int array.
 *
 * \description
 *
 *    The sia must already have been allocated
 *
 * \modifications
 *    - 06/06/06 ACG Created.
 */
FC_ReturnCode fc_initSortedIntArray(
  FC_SortedIntArray *sia /**< input/output - sortedIntArray */
)
{
  if (sia == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  fc_printfLogMessage("Initing sorted int array");

  sia->numVal = 0;
  sia->maxNumVal = 0;
  sia->vals = NULL;
  
  return FC_SUCCESS;
}

/**
 * \ingroup SimpleDataObjects
 * \brief Check to see if a sorted int array is valid.
 *
 * \description
 *
 *    Checks that the sorted int array has "sane" values. Returns
 *    1 if the sorted array looks valid and 0 if it does not.
 *
 * \modifications
 *    - 06/06/06 ACG Created.
 */
int fc_isSortedIntArrayValid(
  FC_SortedIntArray *sia /**< input - sortedIntArray */
)
{
  if (!sia || sia->numVal < 0 || sia->maxNumVal < 0 ||
      sia->numVal > sia->maxNumVal ||
      (sia->maxNumVal == 0 && sia->vals != NULL) ||
      (sia->maxNumVal > 0 && sia->vals == NULL)) {
    return 0;
  }
  else
    return 1;
}

/**
 * \ingroup SimpleDataObjects
 * \brief Copy (deep copy) a sorted int array.
 *
 * \description
 *
 *    The destination sia must already have been allocated but does not have to
 *    be initialized. This peforms a deep copy (i.e. internal data is copied).
 *
 * \modifications
 *    - 08/04/06 WSD Created.
 */
FC_ReturnCode fc_copySortedIntArray(
  FC_SortedIntArray *sia_src, /**< input - sortedIntArray to copy*/
  FC_SortedIntArray *sia_dest /**< output - new sortedIntArray */
)
{
  if (!fc_isSortedIntArrayValid(sia_src) || !sia_dest) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  fc_printfLogMessage("Copying sorted int array");

  sia_dest->numVal = sia_src->numVal;
  sia_dest->maxNumVal = sia_src->numVal;
  if (sia_src->numVal > 0) {
    sia_dest->vals = malloc(sia_src->numVal*sizeof(int));
    if (sia_dest->vals == NULL) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      return FC_MEMORY_ERROR;
    }
    memcpy(sia_dest->vals, sia_src->vals, sia_src->numVal*sizeof(int));
  }
  else {
    sia_dest->vals = NULL;
  }
  
  return FC_SUCCESS;
}

/**
 * \ingroup SimpleDataObjects
 * \brief Get number of values in sorted int array.
 *
 * \modifications
 *    - 08/03/06 WSD created.
 */
FC_ReturnCode fc_getSortedIntArrayNumValue(
  FC_SortedIntArray *sia,  /**< input - sortedIntArray */
  int* numValue            /**< output - number of values */
){
  // default returns
  if (numValue)
    *numValue = -1;

  // check arguments
  if (!fc_isSortedIntArrayValid(sia) || !numValue) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // do it
  *numValue = sia->numVal;

  return FC_SUCCESS;
}

/**
 * \ingroup SimpleDataObjects
 * \brief Get a copy of the values stored in the sorted int array.
 *
 * \description
 *
 *    Returns a copy of the values stores in the sorted int array.  The caller
 *    is responsible for freeing the return array.  If you know you never want
 *    to use the current values of the sorted int array again, you can instead
 *    convert by using \ref fc_convertSortedIntArrayToIntArray.
 *
 * \modifications
 *    - 08/04/06 WSD created.
 */
FC_ReturnCode fc_getSortedIntArrayValues(
  FC_SortedIntArray *sia,  /**< input - sortedIntArray */
  int* numValue,           /**< output - number of values */
  int** values             /**< output - values (needs to be freed) */
){
  // default returns
  if (numValue)
    *numValue = -1;
  if (values)
    *values = NULL;

  // check arguments
  if (!fc_isSortedIntArrayValid(sia) || !numValue || !values) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // easy case - no valus
  if (sia->numVal == 0) {
    *numValue = 0;
    *values = NULL;
    return FC_SUCCESS;
  }

  // do it
  *values = malloc(sia->numVal*sizeof(int));
  if (!(*values)) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  *numValue = sia->numVal;
  memcpy(*values, sia->vals, sia->numVal*sizeof(int));

  return FC_SUCCESS;
}

/**
 * \ingroup SimpleDataObjects
 * \brief Query for int in sorted int array.
 *
 * \description
 *
 *    Returns 1 if the int is in the array, 0 if it is not, and an
 *    FC_ReturnCode if an error occurs.
 *
 * \modifications
 *    - 06/06/06 ACG Created.
 *    - 07/13/06 WSD Separated binary search into it's own routine.
 *      Separated lookupInt into two seperate add and query functions.
 */
int fc_isIntInSortedIntArray(
  FC_SortedIntArray *sia, /**< input - sortedIntArray */
  int ix                  /**< input - int to look for */
){
  FC_ReturnCode rc;
  int foundInt, idx;

  // this call also checks args
  rc = _fc_lookupIntInSortedIntArray(sia, ix, &foundInt, &idx);
  if (rc < 0)
    return rc;

  return foundInt;
}

/**
 * \ingroup SimpleDataObjects
 * \brief Add an int to a sorted int array.
 *
 * \description
 *
 *    Add an int to a sorted int array. Returns 1 if the int is added.
 *    Returns 0 if the int was already there. Returns an FC_ReturnCode
 *    if there is an error.
 *
 * \modifications
 *    - 06/06/06 ACG Created.
 *    - 07/13/06 WSD Separated binary search into it's own routine.
 *      Separated lookupInt into two seperate add and query functions.
 */
int fc_addIntToSortedIntArray(
  FC_SortedIntArray *sia, /**< input/output - sortedIntArray */
  int ix                  /**< input - int to add */
){
  FC_ReturnCode rc;
  int foundInt, idx;

  // this call also checks args
  rc = _fc_lookupIntInSortedIntArray(sia, ix, &foundInt, &idx);
  if (rc < 0)
    return rc;

  // early return
  if (foundInt == 1)
    return 0;

  // do it (this routine allocates memory if needed)
  rc = _fc_addEntryToSortedIntArray(sia, idx, ix);
  if (rc != FC_SUCCESS)
    return rc;
  else
    return 1;
}

/**
 * \ingroup SimpleDataObjects
 * \brief Add an array of ints to a sorted int array.
 *
 * \description
 *
 *    Add an array of ints to a sorted int array. Returns the total number
 *    of ints added (if all ints already existed, would return zero).
 *    Returns an FC_ReturnCode if there is an error.
 *
 *    The extra flag, isSorted, is for performance. A value of 1
 *    means the incoming ints are already sorted (0 means they
 *    may not be).
 *
 *    The incoming int array can have repeated values.
 *
 * \modifications
 *    - 07/13/06 WSD Created.
 */
int fc_addIntArrayToSortedIntArray(
  FC_SortedIntArray *sia, /**< input/output - sortedIntArray */
  int num,                /**< input - number of ints to add */
  int* array,             /**< input - the ints to add */
  int isSorted            /**< input - flag, 1 = incoming ints are sorted */
){
  int i, j;
  FC_ReturnCode rc;
  int *sortedArray;
  int numAdded;
  int newNum, newMaxNum, *newVals;

  // check input
  if (!fc_isSortedIntArrayValid(sia) || num < 0 || (num > 0 && array == NULL)) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // early return - empty array
  if (num == 0)
    return 0;
  
  // sort, if necessary
  if (isSorted)
    sortedArray = array;
  else {
    sortedArray = (int*)malloc(num*sizeof(int));
    if (!sortedArray) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      return FC_MEMORY_ERROR;
    }
    memcpy(sortedArray, array, num*sizeof(int));
    rc = fc_sortIntArray(num, sortedArray);
    if (rc != FC_SUCCESS)
      return rc;
  }

  // make space for result (will replace current sia->vals)
  newMaxNum = sia->numVal + num;
  newVals = (int*)malloc(newMaxNum*sizeof(int));
  if (newVals == NULL){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }

  // make new vals
  numAdded = 0;
  newNum = 0;
  i = j = 0;
  while (i < sia->numVal && j < num) {
    if (sia->vals[i] < sortedArray[j]) {
      newVals[newNum] = sia->vals[i];
      newNum++;
      i++;
    }
    else if (sortedArray[j] < sia->vals[i]) {
      newVals[newNum] = sortedArray[j];
      newNum++;
      numAdded++;
      j++;
      while (j < num && sortedArray[j] == sortedArray[j-1]) // skip duplicates
	j++;
    }
    else {
      newVals[newNum] = sia->vals[i];
      newNum++;
      i++;
      j++;
      while (j < num && sortedArray[j] == sortedArray[j-1]) // skip duplicates
	j++;
    }
  }
  if (i < sia->numVal) {
    memcpy(newVals+newNum, sia->vals+i, (sia->numVal-i)*sizeof(int));
    newNum += sia->numVal-i;
  }
  else if (j < num) {
    while (j < num) {
      newVals[newNum] = sortedArray[j];
      newNum++;
      numAdded++;
      j++;
      while (j < num && sortedArray[j] == sortedArray[j-1]) // skip duplicates
	j++;
    }
  }

  // swap old vals for new
  if (numAdded == 0) {
    free(newVals);
  }
  else {
    free(sia->vals);
    sia->numVal = newNum;
    sia->maxNumVal = newMaxNum;
    sia->vals = newVals;
  }

  // cleanup after sort, if necessary
  if (!isSorted)
    free(sortedArray);

  return numAdded;
}

/**
 * \ingroup SimpleDataObjects
 * \brief Remove an int from a sorted int array.
 *
 * \description
 *
 *    Returns 1 if the int was found and deleted, returns 0 if the int was not
 *    in the array, and returns an error code if something went wrong. (Note:
 *    the total size of the array does not change, call \ref
 *    fc_squeezeSortedIntArray() to save space).
 *
 * \modifications
 *    - 06/06/06 ACG Created.
 *    - 07/13/06 WSD Separated binary search into it's own routine.
 */
int fc_deleteIntFromSortedIntArray(
  FC_SortedIntArray *sia,  /**< input/output - sortedIntArray */
  int ix                   /**< input - int to lookup */
){
  FC_ReturnCode rc;
  int foundInt, idx;

  // This call also checks sia
  rc = _fc_lookupIntInSortedIntArray(sia, ix, &foundInt, &idx);
  if (rc < 0)
    return rc;

  // early return
  if (foundInt == 0)
    return 0;

  // do it
  rc = _fc_deleteEntryFromSortedIntArray(sia, idx);
  if (rc != FC_SUCCESS)
    return rc;
  else
    return 1;
}

/**
 * \ingroup SimpleDataObjects
 * \brief Get first value in the sorted int array (WARNING: does not check for
 *        empty array)
 *
 * \description
 *
 *    Returns the first value in a non empty sorted int array. Do not use
 *    this funciton on empty sorted int arrays. Use \ref
 *    fc_getSortedIntArrayNumValue to determine if an array has values.
 *    To remove the front value, use \ref fc_popSortedIntArrayFront.
 *
 *    You can also manipulate the back of a sorted int array; see also 
 *    \ref fc_getSortedIntArrayBack and fc_popSortedIntArrayBack.
 *
 * \modifications
 *    - 08/03/06 WSD Created.
 */
int fc_getSortedIntArrayFront(
  FC_SortedIntArray *sia, /**< input - sortedIntArray */
  int* value              /**< input - first value in the array */
){
  // default return
  if (value)
    *value = -1;
  
  // check arguments
  if (!fc_isSortedIntArrayValid(sia)) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // do it
  if (sia->numVal > 0) {
    *value = sia->vals[0];
  }
  else {
    fc_printfWarningMessage("Getting the front of an empty sorted int array");
  }

  return FC_SUCCESS;
}

/**
 * \ingroup SimpleDataObjects
 * \brief Get last value in the sorted int array (WARNING: does not check for
 *        empty array)
 *
 * \description
 *
 *    Returns the last value in a non empty sorted int array. Do not use
 *    this funciton on empty sorted int arrays. Use \ref
 *    fc_getSortedIntArrayNumValue to determine if an array has values.
 *    To remove the back value, use \ref fc_popSortedIntArrayBack.
 *
 *    You can also manipulate the front of a sorted int array; see also 
 *    \ref fc_getSortedIntArrayFront and fc_popSortedIntArrayFront.
 *
 * \modifications
 *    - 08/03/06 WSD Created.
 */
int fc_getSortedIntArrayBack(
  FC_SortedIntArray *sia, /**< input - sortedIntArray */
  int* value              /**< input - last value in the array */
){
  // default return
  if (value)
    *value = -1;
  
  // check arguments
  if (!fc_isSortedIntArrayValid(sia)) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // do it
  if (sia->numVal > 0) {
    *value = sia->vals[sia->numVal-1];
  }
  else {
    fc_printfWarningMessage("Getting the front of an empty sorted int array");
  }

  return FC_SUCCESS;
}

/**
 * \ingroup SimpleDataObjects
 * \brief Delete the first value in a sorted int array.
 *
 * \description
 *
 *    Delete the first value in a sorted int array. The return value
 *    is an FC_ReturnCode if there is an error, 0 if the array is
 *    empty, or 1 if the value has been deleted.
 *
 *    If you want to know the value of the front, before popping you must
 *    call \ref fc_getSortedIntArrayFront.
 *
 *    You can also manipulate the back of a sorted int array; see also 
 *    \ref fc_getSortedIntArrayBack and fc_popSortedIntArrayBack.
 *
 * \modifications
 *    - 08/03/06 WSD Created.
 */
int fc_popSortedIntArrayFront(
  FC_SortedIntArray *sia /**< input/output - sortedIntArray */
){
  // check arguments
  if (!fc_isSortedIntArrayValid(sia)) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // do it (do minimum amt of work)
  if (sia->numVal > 0) {
    sia->numVal--;
    if (sia->numVal > 0) 
      memmove(sia->vals, sia->vals+1, sia->numVal*sizeof(int)); 
    return 1;
  }
  else {
    return 0;
  }
}

/**
 * \ingroup SimpleDataObjects
 * \brief Delete the last value in a sorted int array.
 *
 * \description
 *
 *    Delete the last value in a sorted int array. The return value
 *    is an FC_ReturnCode if there is an error, 0 if the array is
 *    empty, or 1 if the value has been deleted.
 *
 *    If you want to know the value of the back, before popping you must
 *    call \ref fc_getSortedIntArrayBack.
 *
 *    You can also manipulate the front of a sorted int array; see also 
 *    \ref fc_getSortedIntArrayFront and fc_popSortedIntArrayFront.
 *
 * \modifications
 *    - 08/03/06 WSD Created.
 */
int fc_popSortedIntArrayBack(
  FC_SortedIntArray *sia /**< input/output - sortedIntArray */
){
  // check arguments
  if (!fc_isSortedIntArrayValid(sia)) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // do it (do minimum amt of work)
  if (sia->numVal > 0) {
    sia->numVal--;
    return 1;
  }
  else {
    return 0;
  }
}

/**
 * \ingroup SimpleDataObjects
 * \brief Realloc the sorted int array so it has no extra space.
 *
 * \modifications
 *    - 06/06/06 ACG Created.
 */
FC_ReturnCode fc_squeezeSortedIntArray(
  FC_SortedIntArray *sia /**< input - sortedIntArray */
){
  int* temp;

  // check args
  if (!fc_isSortedIntArrayValid(sia)){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // early return
  if (sia->numVal == sia->maxNumVal) // already compressed
    return FC_SUCCESS;

  // special case - empty
  if (sia->numVal == 0) {
    free(sia->vals);
    sia->maxNumVal = 0;
    sia->vals = NULL;
    return FC_SUCCESS;
  }

  // do it
  temp = (int*)realloc(sia->vals,(sia->numVal)*sizeof(int));
  if (temp == NULL){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  sia->maxNumVal = sia->numVal;
  sia->vals = temp;

  return FC_SUCCESS;
}

/**
 * \ingroup SimpleDataObjects
 * \brief Free the sorted int array's values.
 *
 * \description
 *
 *    This routine frees the internals of a sorted int array, but not
 *    the sia itself. The user must dealloc the sia itself, unless it
 *    was automatically allocated.
 *
 * \modifications
 *    - 06/06/06 ACG Created.
 */
void fc_freeSortedIntArray(
  FC_SortedIntArray *sia /**< input/output - sortedIntArray */
){

  if (!fc_isSortedIntArrayValid(sia)){
    return; 
  }

  if (sia->vals){
    free(sia->vals);
    sia->numVal = 0;
    sia->maxNumVal = 0;
    sia->vals = NULL;
  }
}

/**
 * \ingroup SimpleDataObjects
 * \brief prints the sortedint array
 *
 * \description
 *
 *    Prints the contents of a sorted int array.
 *
 * \modifications
 *    - 06/06/06 ACG Created.
 *    - 06/07/06 ACG takes ptr to array as parameter
 */
FC_ReturnCode fc_printSortedIntArray(
  FC_SortedIntArray *sia /**< input - sortedIntArray */
){
  int i;

  // check args
  if (!fc_isSortedIntArrayValid(sia)){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  printf("Sorted Int Array:\n");
  printf("    numVal = %d\n", sia->numVal);
  printf("    maxNumVal = %d\n", sia->maxNumVal);
  printf("    vals =");
  if (sia->numVal == 0){
    printf(" (empty)\n");
  }
  else {
    for (i = 0; i < sia->numVal; i++){
      printf(" %d",sia->vals[i]);
    }
    printf("\n");
  }
  fflush(NULL);
  
  return FC_SUCCESS;
}

//@}

/** 
 * \name Sorted Blob Array Routines
 *
 * Sorted, exclusive, self-exanding array of pointers to data blobs. The array
 * is implemented as an array of pointers (void*) to the objects of interest
 * and do not manage the blobs themselves.  The user must also provide a
 * function that compares blobs and returns -1 if the first blob is "less than"
 * the 2nd, 0 if the blobs are equivalent, and 1 if the 2nd glob is greater
 * than first blob (the same as compare functions for qsort).
 */
//-------------------------------
//@{

/**
 * \ingroup SimpleDataObjects
 * \brief Initialize a sorted blob array.
 *
 * \description
 *
 *    The sba must already have been allocated
 *
 * \modifications
 *    - 08/07/06 WSK Created.
 */
FC_ReturnCode fc_initSortedBlobArray(
  FC_SortedBlobArray *sba /**< input/output - sortedBlobArray */
)
{
  if (sba == NULL){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  fc_printfLogMessage("Initing sorted blob array");

  sba->numBlob = 0;
  sba->maxNumBlob = 0;
  sba->blobs = NULL;
  
  return FC_SUCCESS;
}

/**
 * \ingroup SimpleDataObjects
 * \brief Check to see if a sorted blob array is valid.
 *
 * \description
 *
 *    Checks that the sorted blob array has "sane" values. Returns
 *    1 if the sorted array looks valid and 0 if it does not.
 *
 * \modifications
 *    - 06/06/06 ACG Created.
 */
int fc_isSortedBlobArrayValid(
  FC_SortedBlobArray *sba /**< input - sortedBlobArray */
)
{
  if (!sba || sba->numBlob < 0 || sba->maxNumBlob < 0 ||
      sba->numBlob > sba->maxNumBlob ||
      (sba->maxNumBlob == 0 && sba->blobs != NULL) ||
      (sba->maxNumBlob > 0 && sba->blobs == NULL)) {
    return 0;
  }
  else
    return 1;
}

/**
 * \ingroup SimpleDataObjects
 * \brief Get number of blobs in sorted blob array.
 *
 * \modifications
 *    - 08/08/06 WSD created.
 */
FC_ReturnCode fc_getSortedBlobArrayNumBlob(
  FC_SortedBlobArray *sba,  /**< input - sortedBlobArray */
  int* numBlob              /**< output - number of blobs */
){
  // default returns
  if (numBlob)
    *numBlob = -1;

  // check arguments
  if (!fc_isSortedBlobArrayValid(sba) || !numBlob) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // do it
  *numBlob = sba->numBlob;

  return FC_SUCCESS;
}

/**
 * \ingroup SimpleDataObjects
 * \brief Get a copy of the blob locations stored in the sorted blob array.
 *
 * \description
 *
 *    Returns a copy of the blob locations storee in the sorted blob array.
 *    The caller is responsible for freeing the returned array.
 *
 * \modifications
 *    - 08/08/06 WSD created.
 */
FC_ReturnCode fc_getSortedBlobArrayBlobs(
  FC_SortedBlobArray *sba,  /**< input - sortedBlobArray */
  int* numBlob,             /**< output - number of values */
  void*** blobs             /**< output - blob locations (needs to be freed) */
){
  // default returns
  if (numBlob)
    *numBlob = -1;
  if (blobs)
    *blobs = NULL;

  // check arguments
  if (!fc_isSortedBlobArrayValid(sba) || !numBlob || !blobs) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // easy case - no valus
  if (sba->numBlob == 0) {
    *numBlob = 0;
    *blobs = NULL;
    return FC_SUCCESS;
  }

  // do it
  *blobs = malloc(sba->numBlob*sizeof(void*));
  if (!(*blobs)) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  *numBlob = sba->numBlob;
  memcpy(*blobs, sba->blobs, sba->numBlob*sizeof(void*));

  return FC_SUCCESS;
}

/**
 * \ingroup SimpleDataObjects
 * \brief Query for blob in sorted blob array.
 *
 * \description
 *
 *    Returns 1 if the blob is in the array, 0 if it is not, and an
 *    FC_ReturnCode if an error occurs. Also returns pointer to
 *    the blob already in the array if requested (pass NULL if you
 *    don't want it).
 *
 * \modifications
 *    - 08/08/2006 WSD Created.
 */
int fc_isBlobInSortedBlobArray(
  FC_SortedBlobArray *sba, /**< input - sortedBlobArray */
  void* blob,              /**< input - blob to look for */
  int blobCmp(const void*, const void*), /**< input - the blob comparison 
					    function */
  void** foundBlob         /**< output - (optional) the found blob */
){
  FC_ReturnCode rc;
  int foundInt, idx;

  // default return
  if (foundBlob)
    *foundBlob = NULL;

  // this call also checks args
  rc = _fc_lookupBlobInSortedBlobArray(sba, blob, blobCmp, &foundInt, &idx);
  if (rc < 0)
    return rc;

  if (foundInt && foundBlob)
    *foundBlob = sba->blobs[idx];

  return foundInt;
}

/**
 * \ingroup SimpleDataObjects
 * \brief Add a blob to a sorted blob array.
 *
 * \description
 *
 *    Add a blob to a sorted blob array. Returns 1 if the blob is added.
 *    Returns 0 if the blob was already there. Returns an FC_ReturnCode
 *    if there is an error.
 *
 * \modifications
 *    - 08/08/2006 WSD Created.
 */
int fc_addBlobToSortedBlobArray(
  FC_SortedBlobArray *sba, /**< input/output - sortedBlobArray */
  void* blob,              /**< input - blob to add */
  int blobCmp(const void*, const void*) /**< input - the blob comparison 
					    function */
){
  FC_ReturnCode rc;
  int foundInt, idx;

  // this call also checks args
  rc = _fc_lookupBlobInSortedBlobArray(sba, blob, blobCmp, &foundInt, &idx);
  if (rc < 0)
    return rc;

  // early return
  if (foundInt == 1)
    return 0;

  // do it (this routine allocates memory if needed)
  rc = _fc_addEntryToSortedBlobArray(sba, idx, blob);
  if (rc != FC_SUCCESS)
    return rc;
  else
    return 1;
}

/**
 * \ingroup SimpleDataObjects
 * \brief Remove a blob from a sorted blob array.
 *
 * \description
 *
 *    Returns 1 if the blob was found and deleted, returns 0 if the blob was not
 *    in the array, and returns an error code if something went wrong. (Note:
 *    the total size of the array does not change, call \ref
 *    fc_squeezeSortedBlobArray() to save space). Also returns pointer to
 *    the blob deleted from the array if requested (pass NULL if you
 *    don't want it).
 *
 * \modifications
 *    - 08/08/2006 WSD Created.
 */
int fc_deleteBlobFromSortedBlobArray(
  FC_SortedBlobArray *sba, /**< input/output - sortedBlobArray */
  void* blob,              /**< input - blob to delete */
  int blobCmp(const void*, const void*), /**< input - the blob comparison 
					    function */
  void** deletedBlob       /**< output - (optional) the deleted blob */
){
  FC_ReturnCode rc;
  int foundInt, idx;

  // default return
  if (deletedBlob)
    *deletedBlob = NULL;

  // This call also checks sia
  rc = _fc_lookupBlobInSortedBlobArray(sba, blob, blobCmp, &foundInt, &idx);
  if (rc < 0)
    return rc;

  // early return
  if (foundInt == 0)
    return 0;

  // do it
  if (deletedBlob)
    *deletedBlob = sba->blobs[idx];
  rc = _fc_deleteEntryFromSortedBlobArray(sba, idx);
  if (rc != FC_SUCCESS)
    return rc;
  else
    return 1;
}

/**
 * \ingroup SimpleDataObjects
 * \brief Realloc the sorted blob array so it has no extra space.
 *
 * \modifications
 *    - 06/06/06 ACG Created.
 */
FC_ReturnCode fc_squeezeSortedBlobArray(
  FC_SortedBlobArray *sba /**< input - sortedBlobArray */
){
  void** temp;

  // check args
  if (!fc_isSortedBlobArrayValid(sba)){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // early return
  if (sba->numBlob == sba->maxNumBlob) // already compressed
    return FC_SUCCESS;

  // special case - empty
  if (sba->numBlob == 0) {
    free(sba->blobs);
    sba->maxNumBlob = 0;
    sba->blobs = NULL;
    return FC_SUCCESS;
  }

  // do it
  temp = (void**)realloc(sba->blobs, (sba->numBlob)*sizeof(void*));
  if (temp == NULL){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  sba->maxNumBlob = sba->numBlob;
  sba->blobs = temp;

  return FC_SUCCESS;
}

/**
 * \ingroup SimpleDataObjects
 * \brief Free the sorted blob array's values.
 *
 * \description
 *
 *    This routine frees the internals of a sorted blob array, but not
 *    the sba itself. The user must dealloc the sba itself, unless it
 *    was automatically allocated.
 *
 * \modifications
 *    - 06/06/06 ACG Created.
 */
void fc_freeSortedBlobArray(
  FC_SortedBlobArray *sba /**< input/output - sortedBlobArray */
){

  if (!fc_isSortedBlobArrayValid(sba)){
    return; 
  }

  if (sba->blobs){
    free(sba->blobs);
    sba->numBlob = 0;
    sba->maxNumBlob = 0;
    sba->blobs = NULL;
  }
}

//@}

/** 
 * \name Conversion routines
 *
 * Conversion implies that original data object will be consumed (no longer
 * available) in the process.
 */
//-------------------------------
//@{

/**
 * \ingroup SimpleDataObjects
 * \brief Convert int array (no duplicates) to sorted int array.
 *
 * \description
 *
 *    Convert the int array to a sorted int array. The input array is assumed
 *    to have no duplicates. If it does, use 
 *    \ref fc_addIntArrayToSortedIntArray. For performance, there is a flag 
 *    to specify if the array has already been sorted.
 *
 *    The sorted int array has to already have been allocated. The int array is
 *    "destroyed" in that the responsibility for freeing the array is
 *    transfered to the sorted int array. 
 *
 *    WARNING: Do not use this routine on an automatically allocated array
 *    (i.e. one declared as int array[8]). Instead use \ref
 *    fc_addIntArrayToSortedIntArray().
 *
 *    NOTE: It was decided not to handle duplicates because this is meant
 *    to be a quick conversion.
 *
 * \modifications
 *    - 08/04/06 WSD Created.
 */
int fc_convertIntArrayToSortedIntArray(
  int num,                /**< input - number of ints in array */
  int* array,             /**< input - the array of ints */
  int isSorted,           /**< input - sorted flag, 1 = sorted, 0 = not sorted */
  FC_SortedIntArray *sia  /**< input/output - sortedIntArray */
){
  FC_ReturnCode rc;

  // default return
  // ??? init the sia?
  
  // check arguments
  if (num < 0 || (num > 0 && array == NULL) || !sia) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // setup sia
  rc = fc_initSortedIntArray(sia);
  if (rc != FC_SUCCESS)
    return rc;

  // simple case - no entries
  if (num == 0)
    return FC_SUCCESS;

  // sort 
  if (!isSorted)
    fc_sortIntArray(num, array);

  // do it
  sia->numVal = num;
  sia->maxNumVal = num;
  sia->vals = array;

  return FC_SUCCESS;
}

/**
 * \ingroup SimpleDataObjects
 * \brief Convert sorted int array to array.
 *
 * \description
 *
 *    Convert the sorted int array to an int array. The sorted int array
 *    will be "destroyed" in that it's contents will be freed (and set back
 *    to zero). The sorted int array may still need to be freed it was
 *    malloc'd. 
 *
 * \modifications
 *    - 08/04/06 WSD Created.
 */
int fc_convertSortedIntArrayToIntArray(
  FC_SortedIntArray *sia, /**< input/output - sortedIntArray */
  int* num,               /**< input - number of ints in array */
  int** array             /**< input - the array of ints */
){
  FC_ReturnCode rc;

  // default return
  if (num)
    *num = -1;
  if (array)
    *array = NULL;
  
  // check arguments
  if (!sia || !num || !array) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // minimize space usage
  rc = fc_squeezeSortedIntArray(sia);
  if (rc != FC_SUCCESS)
    return rc;

  // do it
  *num = sia->numVal;
  *array = sia->vals;
  sia->numVal = 0;
  sia->maxNumVal = 0;
  sia->vals = NULL;

  return FC_SUCCESS;
}

/**
 * \ingroup  SimpleDataObjects
 * \brief Create an Int Array from a Mask
 *
 * \description
 *       
 *    Creates an Int Array from the given mask. The original mask
 *    is left untouched. The returned Int Array will be sorted.
 *
 * \modifications 
 *    - 07/01/04 WSK, Created
 */
FC_ReturnCode fc_createIntArrayFromMask(
  int numMask,    /**< input - number of entries in Mask */
  int* mask,       /**< input - mask array */ 
  int* numArray,  /**< output - number of entries in returned array */
  int **array     /**< output - returned int array, will be sorted */ 
) {
  int i, num;

  // default return
  if (numArray)
    *numArray = -1;
  if (array)
    *array = NULL;
  
  // input checking
  if (numMask < 1 || !mask || !numArray || !array) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }  

  // log message
  fc_printfLogMessage("Creating integer array from mask."); 

  // do it
  num = 0;
  for (i = 0; i < numMask; i++)
    if (mask[i] > 0)
      num++;
  *numArray = num;
  if (num > 0) {
    *array = malloc(num*sizeof(int));
    if (*array == NULL) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      return FC_MEMORY_ERROR;
    }
    num = 0;
    for (i = 0; i < numMask; i++) {
      if (mask[i] > 0) {
        (*array)[num] = i;
        num++;
      }
    }
  }

  return FC_SUCCESS;
}

/**
 * \ingroup  SimpleDataObjects
 * \brief Create a Mask from an Int Array
 *
 * \description
 *       
 *    Creates a Mask from the given Int Array; the length of the mask
 *    must also be provided. The original array is left untouched.
 *    numArray = 0 is a valid input that creates a mask with all
 *    entries being 0. If any entries of array are outside the range
 *    [0, numMask-1], the routine does not create the mask and
 *    returns an error.
 *
 * \modifications 
 *    - 07/01/04 WSK, Created
 */
FC_ReturnCode fc_createMaskFromIntArray(
  int numArray,  /**< input - number of entries in array */
  int *array,    /**< input - the array */ 
  int numMask,   /**< input - number of entries in Mask */
  int** mask     /**< output - the mask array */ 
) {
  int i;

  // default return
  if (mask)
    *mask = NULL;
  
  // input checking
  if (numArray < 0 || (numArray > 0 && !array) || numMask < 1 || !mask) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }  

  // log message
  fc_printfLogMessage("Creating mask from integer array."); 

  // do it
  *mask = calloc(numMask, sizeof(int));
  if (*mask == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  
  for (i = 0; i < numArray; i++) {
    if (array[i] < 0 || array[i] > numMask-1) {
      free(*mask);
      *mask = NULL;
      fc_printfErrorMessage("Input Error: array element has invalid value.");
      return FC_INPUT_ERROR;
    }
    (*mask)[array[i]] = 1;
  }

  return FC_SUCCESS;
}

//@}

/**
 * \ingroup PrivateSimpleDataObjects
 * \brief Expand an int array.
 *
 * \description
 *
 *    If length == 0, an array of length 1 is generated.
 *    If length > 0, then the array length is doubled. If length < 0,
 *    an error is generated.
 *
 * \modifications
 *    - 08/07/06 WSD Created.
 */
FC_ReturnCode _fc_expandIntArray(
  int* length,  /**< input/output - length of the array */
  int** array   /**< input/output - the array */
){
  int new_length;
  int* new_array;

  // test args (note, array doesn't have to be null if length = 0)
  if (!length || *length < 0 || !array) 
    return FC_INPUT_ERROR;

  // do it
  new_length = *length > 0 ? 2*(*length) : 1;
  if (*length == 0)  // malloc in case they didn't pass in NULL
    new_array = malloc(new_length*sizeof(int));
  else 
    new_array = realloc(*array, new_length*sizeof(int));
  if (!new_array) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
 
  // return
  *length = new_length;
  *array = new_array;
  return FC_SUCCESS;
}

/**
 * \ingroup PrivateSimpleDataObjects
 * \brief Lookup int in the sorted int array.
 *
 * \description
 *
 *    Look for an int in a sorted int array. A flag for whether the int
 *    is found is returned in foundInt (1 = found, 0 = not found). The index
 *    that the int is at (or would be at if it was in the array) is returned
 *    in idx.
 *
 *    This is implemented as binary search and should have average run
 *    time of O(log(sia->numVal)).
 *
 * \modifications
 *    - 06/06/06 ACG Created.
 *    - 07/13/06 WSD Separated binary search into it's own routine.
 */
FC_ReturnCode _fc_lookupIntInSortedIntArray(
  FC_SortedIntArray *sia,  /**< input - sortedIntArray */
  int i,                   /**< input - int to lookup */    
  int* foundInt,           /**< output - 1 if find int, or 0 */
  int* idx                 /**< output - the index if int is found, or
			      the index the int would have been at */
){
  int low, high, mid;

  // default returns
  if (foundInt)
    *foundInt = -1;
  if (idx)
    *idx = -1;

  // check arguments
  if (!fc_isSortedIntArrayValid(sia) || !foundInt || !idx) {
    return FC_INPUT_ERROR;
  }

  // early return
  if (sia->numVal < 1) {
    *foundInt = 0;
    *idx = 0;
    return FC_SUCCESS;
  }

  // binary search
  low = 0;
  high = sia->numVal-1;
  while (low <= high) {
    mid = (low+high)/2;
    if (i == sia->vals[mid]) {
      *foundInt = 1;
      *idx = mid;
      return FC_SUCCESS;
    } else if (i < sia->vals[mid]) { // val is before mid
      high = mid-1;
    } else {
      low = mid+1;
    }
  }

  // no match
  *foundInt = 0;
  *idx = (i < sia->vals[mid] ? mid : mid+1);
  return FC_SUCCESS;
}

/**
 * \ingroup PrivateSimpleDataObjects
 * \brief Add entry to a sorted int array at the specified location.
 *
 * \description
 *
 *    Use this after you know where to put the new entry (usually
 *    use \ref _fc_lookupIntInSortedIntArray).
 *
 *    NOTE! Make sure that use of this function does not violate the sorted
 *    and uniqueness of the sorted int array!
 *
 * \modifications
 *    - 08/04/06 WSD Created.
 */
FC_ReturnCode _fc_addEntryToSortedIntArray(
  FC_SortedIntArray *sia, /**< input/output - sortedIntArray */
  int idx,                /**< input - location to add the value */
  int value               /**< input - the value to add */
)
{
  // check args
  if (!fc_isSortedIntArrayValid(sia) || idx < 0 || idx > sia->numVal) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // make sure there is room
  if (sia->numVal >= sia->maxNumVal)
    _fc_expandIntArray(&sia->maxNumVal, &sia->vals);

  // do it
  if (idx < sia->numVal) {
    memmove((sia->vals)+idx+1,(sia->vals)+idx,
	    (sia->numVal-idx)*sizeof(int));
  }
  sia->vals[idx] = value;
  sia->numVal++;

  return FC_SUCCESS;
}

/**
 * \ingroup PrivateSimpleDataObjects
 * \brief Remove an entry at a specified location from a sorted int array.
 *
 * \description
 *
 *    Use this after you know which entry you want to delete (usually
 *    use \ref _fc_lookupIntInSortedIntArray).
 *
 *    This routine will fail if you try to delete an entry from an empty 
 *    array.
 *
 * \modifications
 *    - 08/04/06 WSD Created.
 */
FC_ReturnCode _fc_deleteEntryFromSortedIntArray(
  FC_SortedIntArray *sia, /**< input/output - sortedIntArray */
  int idx                 /**< input - location to remove value */
)
{
  // check args
  if (!fc_isSortedIntArrayValid(sia) || idx < 0 || idx > sia->numVal-1) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // do it
  if (idx+1 < sia->numVal) {
    memmove((sia->vals)+idx,(sia->vals)+idx+1,
	    (sia->numVal-(idx+1))*sizeof(int));
  }
  sia->numVal--;
  
  return FC_SUCCESS;
}

/**
 * \ingroup PrivateSimpleDataObjects
 * \brief Expand a blob array.
 *
 * \description
 *
 *    If length == 0, an array of length 1 is generated.
 *    If length > 0, then the array length is doubled. If length < 0,
 *    an error is generated.
 *
 * \modifications
 *    - 08/07/06 WSD Created.
 */
FC_ReturnCode _fc_expandBlobArray(
  int* length,   /**< input/output - length of the array */
  void*** array  /**< input/output - the array */
){
  int new_length;
  void** new_array;

  // test args (note, array doesn't have to be null if length = 0)
  if (!length || *length < 0 || !array) 
    return FC_INPUT_ERROR;

  // do it
  new_length = *length > 0 ? 2*(*length) : 1;
  if (*length == 0)  // malloc in case they didn't pass in NULL
    new_array = malloc(new_length*sizeof(void*));
  else 
    new_array = realloc(*array, new_length*sizeof(void*));
  if (!new_array) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
 
  // return
  *length = new_length;
  *array = new_array;
  return FC_SUCCESS;
}

/**
 * \ingroup PrivateSimpleDataObjects
 * \brief Lookup blob in the sorted blob array.
 *
 * \description
 *
 *    Look for an blob in a sorted blob array using the provided compare
 *    function. A flag for whether a matching blob is found is returned in
 *    foundBlob (1 = found, 0 = not found). The index that the int is at (or
 *    would be at if it was in the array) is returned in idx.
 *
 *    This is implemented as binary search and should have average run
 *    time of O(log(sia->numVal)).
 *
 * \modifications
 *    - 06/06/06 ACG Created.
 *    - 07/13/06 WSD Separated binary search into it's own routine.
 */
FC_ReturnCode _fc_lookupBlobInSortedBlobArray(
  FC_SortedBlobArray *sba,  /**< input - sortedIntArray */
  void* blob,               /**< input - blob to lookup */
  int blobCmp(const void* blob1, const void* blob2), /**< input - the
							  compare function */
  int* foundBlob,           /**< output - 1 if find blob, or 0 */
  int* idx                  /**< output - the index if blob is found, or
			       the index the blob would have been at */
){
  int low, high, mid, cmpVal;

  // default returns
  if (foundBlob)
    *foundBlob = -1;
  if (idx)
    *idx = -1;

  // check arguments
  //?FIX could blob be null?
  if (!fc_isSortedBlobArrayValid(sba) || !blob || !blobCmp || !foundBlob ||
      !idx) {
    return FC_INPUT_ERROR;
  }

  // early return
  if (sba->numBlob < 1) {
    *foundBlob = 0;
    *idx = 0;
    return FC_SUCCESS;
  }

  // binary search
  low = 0;
  high = sba->numBlob-1;
  while (low <= high) {
    mid = (low+high)/2;
    cmpVal = blobCmp((const void*)blob, (const void*)sba->blobs[mid]);
    if (cmpVal == 0) {
      *foundBlob = 1;
      *idx = mid;
      return FC_SUCCESS;
    } else if (cmpVal < 0) { // blob is before mid
      high = mid-1;
    } else {
      low = mid+1;
    }
  }

  // no match
  *foundBlob = 0;
  *idx = (cmpVal < 0 ? mid : mid+1);
  return FC_SUCCESS;
}

/**
 * \ingroup PrivateSimpleDataObjects
 * \brief Add entry to a sorted blob array at the specified location.
 *
 * \description
 *
 *    Use this after you know where to put the new entry (usually
 *    use \ref _fc_lookupBlobInSortedBlobArray).
 *
 *    NOTE! Make sure that use of this function does not violate the sorted
 *    and uniqueness of the sorted blob array!
 *
 * \modifications
 *    - 08/08/06 WSD Created.
 */
FC_ReturnCode _fc_addEntryToSortedBlobArray(
  FC_SortedBlobArray *sba, /**< input/output - sortedBlobArray */
  int idx,                 /**< input - location to add the blob */
  void* blob_p             /**< input - the pointer to the blob */
)
{
  // check args
  if (!fc_isSortedBlobArrayValid(sba) || idx < 0 || idx > sba->numBlob ||
      !blob_p) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // make sure there is room
  if (sba->numBlob >= sba->maxNumBlob)
    _fc_expandBlobArray(&sba->maxNumBlob, &sba->blobs);

  // do it
  if (idx < sba->numBlob) {
    memmove((sba->blobs)+idx+1,(sba->blobs)+idx,
	    (sba->numBlob-idx)*sizeof(void*));
  }
  sba->blobs[idx] = blob_p;
  sba->numBlob++;

  return FC_SUCCESS;
}

/**
 * \ingroup PrivateSimpleDataObjects
 * \brief Remove an entry at a specified location from a sorted blob array.
 *
 * \description
 *
 *    Use this after you know which entry you want to delete (usually
 *    use \ref _fc_lookupBlobInSortedBlobArray).
 *
 *    This routine will fail if you try to delete an entry from an empty 
 *    array.
 *
 * \modifications
 *    - 08/08/06 WSD Created.
 */
FC_ReturnCode _fc_deleteEntryFromSortedBlobArray(
  FC_SortedBlobArray *sba, /**< input/output - sortedBlobArray */
  int idx                  /**< input - location to remove value */
)
{
  // check args
  if (!fc_isSortedBlobArrayValid(sba) || idx < 0 || idx > sba->numBlob-1) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // do it
  if (idx+1 < sba->numBlob) {
    memmove((sba->blobs)+idx,(sba->blobs)+idx+1,
	    (sba->numBlob-(idx+1))*sizeof(void*));
  }
  sba->numBlob--;
  
  return FC_SUCCESS;
}

