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
 * \file storageP.h
 * \brief Private declarations for \ref SimpleDataObjects Module.
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/modules/storageP.h,v $
 * $Revision: 1.6 $ 
 * $Date: 2006/09/19 00:57:58 $
 */

#ifndef _FC_STORAGE_P_H_
#define _FC_STORAGE_P_H_

#ifdef __cplusplus
extern "C" {
#endif

// sia helpers
FC_ReturnCode _fc_expandIntArray(int* length, int** array);
FC_ReturnCode _fc_lookupIntInSortedIntArray(FC_SortedIntArray *sia, int i,
					    int* foundInt, int* idx);
FC_ReturnCode _fc_addEntryToSortedIntArray(FC_SortedIntArray *sia, int idx,
					   int value);
FC_ReturnCode _fc_deleteEntryFromSortedIntArray(FC_SortedIntArray *sia, 
						int idx);

// sba helpers
FC_ReturnCode _fc_expandBlobArray(int* length, void*** array);
FC_ReturnCode _fc_lookupBlobInSortedBlobArray(FC_SortedBlobArray *sia, 
			 void* blob_p,
			 int blobCmp(const void* blob1_p, const void* blob2_p),
			 int* foundBlob, int* idx);
FC_ReturnCode _fc_addEntryToSortedBlobArray(FC_SortedBlobArray *sba, int idx,
					   void* blob_p);
FC_ReturnCode _fc_deleteEntryFromSortedBlobArray(FC_SortedBlobArray *sba, 
                         int idx);


#ifdef __cplusplus
}
#endif

#endif // _FC_STORAGE_P_H_

