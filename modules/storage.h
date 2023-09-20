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
 * \file storage.h
 * \brief Public declarations for \ref SimpleDataObjects Module.
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/modules/storage.h,v $
 * $Revision: 1.12 $ 
 * $Date: 2006/08/30 19:20:01 $
 */

#ifndef _FC_STORAGE_H_
#define _FC_STORAGE_H_

#ifdef __cplusplus
extern "C" {
#endif

// Simple data objects: int arrays, linked lists & masks

/** \name Typedefs */
//-------------------------------
//@{

/**
 * \ingroup   SimpleDataObjects
 * \brief Coordinates data structure
 *
 * \description
 * 
 *    Structure to hold coordinates. Maximum dimensionality is 3. 
 *    (This is essentially a double[3] declaration).
 */
typedef double FC_Coords[3];

/**
 * \ingroup   SimpleDataObjects
 * \brief Vector data structure
 *
 * \description
 * 
 *    Structure to hold a vector, assuming that it originates
 *    at the origin. Maximum dimensionality is 3. (This is
 *    essentially a double[3] declaration).
 *
 */
typedef double FC_Vector[3];

/**
 * \ingroup  SimpleDataObjects
 * \brief Sorted, exclusive, self-exanding int array
 *
 * \description
 *
 *   Use this object and it's functions to maintain a sorted int array.
 *   The array is self expanding and maintains a sorted, unique, set
 *   of integers. The array does not automatically shrink, but \ref 
 *   fc_squeezeSortedIntArray() can be used to minimize memory usage.
 *
 * \modifications 
 *    - 2006/06/06 ACG created
 *    - 2008/05/02 CDU Updated to support getting unsorted vals back
 */
typedef struct {
  int  numVal;     /**< number of values */
  int  maxNumVal;  /**< total allocated space (max possible number of values) */
  int* vals;      /**< int array holding the values */
  int  nextTicket; /**< the next unique ticket to assign to an inserted val */
  int* tickets;   /**< each val has an insert ticket so we can extract vals in order they were inserted */
} FC_SortedIntArray;

/**
 * \ingroup  SimpleDataObjects
 * \brief Sorted, exclusive, self-exanding array of data blobs
 *
 * \description
 *
 *   Use this object and it's functions to maintain a sorted array of data
 *   "blobs". It's implemented as an array of void pointers and the user must
 *   provide the size of the blobs and a comparison function. The array is self
 *   expanding and maintains a sorted, unique, set of the blobs integers. The
 *   array does not automatically shrink, but \ref fc_squeezeSortedBlobArray()
 *   can be used to minimize memory usage.
 *
 * \todo ?Store the compare function and have it provided during init? What
 *     Happens on free?
 *
 * \modifications 
 *    - 06/07/06 WSD created
 */
typedef struct {
  int numBlob;     /**< number of blobs */
  int maxNumBlob;  /**< total allocated space (max possible number of blobs) */
  void** blobs;    /**< array of blobs */
} FC_SortedBlobArray;

//@}

//sorted int array routines
FC_ReturnCode fc_initSortedIntArray(FC_SortedIntArray *sia);
int fc_isSortedIntArrayValid(FC_SortedIntArray* sia);
FC_ReturnCode fc_copySortedIntArray(FC_SortedIntArray *sia_src,
				   FC_SortedIntArray *sia_dest);
FC_ReturnCode fc_getSortedIntArrayNumValue(FC_SortedIntArray *sia, 
		        int* numVal);
FC_ReturnCode fc_getSortedIntArrayValues(FC_SortedIntArray *sia,
			int* numVal, int** values);
FC_ReturnCode fc_getSortedIntArrayValuesInInsertOrder(FC_SortedIntArray *sia,  
		        int* numValue, int** values);

int fc_isIntInSortedIntArray(FC_SortedIntArray *sia, int i);
int fc_addIntToSortedIntArray(FC_SortedIntArray *sia, int i);
int fc_addIntArrayToSortedIntArray(FC_SortedIntArray *sia, int num, 
				   int* array, int isSorted);
int fc_deleteIntFromSortedIntArray(FC_SortedIntArray *sia, int i);
FC_ReturnCode fc_getSortedIntArrayFront(FC_SortedIntArray *sia, int* value);
FC_ReturnCode fc_getSortedIntArrayBack(FC_SortedIntArray *sia, int* value); 
int fc_popSortedIntArrayFront(FC_SortedIntArray *sia);
int fc_popSortedIntArrayBack(FC_SortedIntArray *sia);
FC_ReturnCode fc_squeezeSortedIntArray(FC_SortedIntArray *sia);
void fc_freeSortedIntArray(FC_SortedIntArray *sia);
FC_ReturnCode fc_printSortedIntArray(FC_SortedIntArray *sia);

// sorted blob array routines
FC_ReturnCode fc_initSortedBlobArray(FC_SortedBlobArray *sba);
int fc_isSortedBlobArrayValid(FC_SortedBlobArray* sba);
FC_ReturnCode fc_getSortedBlobArrayNumBlob(FC_SortedBlobArray *sba, 
		        int* numBlob);
FC_ReturnCode fc_getSortedBlobArrayBlobs(FC_SortedBlobArray *sba,
			int* numBlob, void*** blobs);
int fc_isBlobInSortedBlobArray(FC_SortedBlobArray *sba, void* blob, 
		   int blobCmp(const void* blob1, const void* blob2),
		   void** found_blob);
int fc_addBlobToSortedBlobArray(FC_SortedBlobArray *sba, void* blob,
		   int blobCmp(const void* blob1, const void* blob2));
int fc_deleteBlobFromSortedBlobArray(FC_SortedBlobArray *sba, void* blob,
                   int blobCmp(const void* blob1, const void* blob2),
		   void** deleted_blob);
FC_ReturnCode fc_squeezeSortedBlobArray(FC_SortedBlobArray *sba);
void fc_freeSortedBlobArray(FC_SortedBlobArray *sba);

// routines to convert between simple data objects

// array/sia
FC_ReturnCode fc_convertIntArrayToSortedIntArray(int num, int* array, 
                        int isSorted, FC_SortedIntArray* sia);
FC_ReturnCode fc_convertSortedIntArrayToIntArray(FC_SortedIntArray* sia,
			int* num, int** array);

// array/mask
FC_ReturnCode fc_createIntArrayFromMask(int numMask, int* mask, int* numArray,
					int** array);
FC_ReturnCode fc_createMaskFromIntArray(int numArray, int* array, int numMask, 
					int** mask);


#ifdef __cplusplus
}
#endif

#endif // _FC_STORAGE_H_

