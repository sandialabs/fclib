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
 * \file dataset.h
 * \brief Public declarations for the \ref Dataset Module.
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/modules/dataset.h,v $
 * $Revision: 1.15 $ 
 * $Date: 2006/09/22 06:39:53 $
 */

#ifndef _FC_DATASET_H_
#define _FC_DATASET_H_

#ifdef __cplusplus
extern "C" {
#endif

// create new dataset from scratch
//-------------------------------
FC_ReturnCode fc_createDataset(char *name, FC_Dataset *dataset);

// other ways to get new datasets
//-------------------------------
FC_ReturnCode fc_copyDataset(FC_Dataset src_ds, char* newDsName, 
			     FC_Dataset* new_ds);

// get datasets
//-------------------------------
FC_ReturnCode fc_getDatasets(int* numDataset, FC_Dataset** datasets);
FC_ReturnCode fc_getNumDataset(int* numDataset);
FC_ReturnCode fc_getDatasetByName(char* datasetName, int *numDatasets,
				  FC_Dataset** dataset);

// change the name of a dataset
//-------------------------------
FC_ReturnCode fc_changeDatasetName(FC_Dataset dataset, char* newName);

// release/delete dataset
//-------------------------------
FC_ReturnCode fc_releaseDataset(FC_Dataset dataset);
FC_ReturnCode fc_deleteDataset(FC_Dataset dataset);

// get Dataset metadata
//-------------------------------
int fc_isDatasetValid(FC_Dataset dataset);
FC_ReturnCode fc_getDatasetName(FC_Dataset dataset, char** dsName);

// Print
//-------------------------------
FC_ReturnCode fc_printDataset(FC_Dataset dataset, char* label);

#ifdef __cplusplus
}
#endif

#endif // _FC_DATASET_H_
