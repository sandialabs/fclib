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
 * \file fileio.h
 * \brief Public declarations for \ref FileIO module.
 *
 * \revisions
 *  - 10/15/2007 ACG removed saf
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/modules/fileio.h,v $
 * $Revision: 1.11 $ 
 * $Date: 2006/08/30 19:20:01 $
 */

#ifndef _FC_FILE_IO_H_
#define _FC_FILE_IO_H_

#ifdef __cplusplus
extern "C" {
#endif

/** \name File IO enumerated type. */
//-------------------------------
//@{

/**
 * \ingroup FileIO
 * \brief Types of file formats for dataset file IO
 *
 * \description
 *
 *    All possible file types are defined, but not all may be available
 *    in a specific installation.
 *
 *    NOTE: purposefully do not have this type stuff in base.h because
 *    then only this set of files (fileio) needs to be modified if
 *    different format types are available.
 *
 * \todo ? Add html pointers for each type?
 * \todo ? Is the NOTE above true? move this & helpers to base.h?
 *
 * \modifications
 *    - 7/13/05 WSD Created.
 */
typedef enum {
  FC_FT_NONE = 0,   /**< Type not specified */
  FC_FT_EXODUS,     /**< Exodus II */
  FC_FT_LSDYNA,     /**< LS-DYNA */
} FC_FileIOType;

//@}

// helpers
int fc_isFileIOTypeValid(FC_FileIOType fileType);
char* fc_getFileIOTypeText(FC_FileIOType fileType);
int fc_isFileIOTypeReadSupported(FC_FileIOType fileType);
int fc_isFileIOTypeWriteSupported(FC_FileIOType fileType);

// Load dataset
// FIX? someday add r/w option
FC_ReturnCode fc_loadDataset(char *filename, FC_Dataset *dataset);
FC_ReturnCode fc_loadDatasetWithFormat(char *filename, FC_FileIOType fileType,
			FC_Dataset *dataset);

// Write dataset
// FIX? add a clobber option?
FC_ReturnCode fc_writeDataset(FC_Dataset dataset, char* filename,
			FC_FileIOType fileType);
FC_ReturnCode fc_rewriteDataset(FC_Dataset dataset, char* filename,
			FC_FileIOType fileType);

// Other
FC_ReturnCode fc_getDatasetFileIOType(FC_Dataset dataset, 
                        FC_FileIOType *fileType);

#ifdef __cplusplus
}
#endif

#endif // _FC_FILE_IO_H_
