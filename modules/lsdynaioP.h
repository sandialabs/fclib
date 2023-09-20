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
 * \file lsdynaioP.h
 * \brief Private declarations for \ref LSDynaFileIO module.
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/modules/lsdynaioP.h,v $
 * $Revision: 1.17 $ 
 * $Date: 2006/09/15 22:09:36 $
 */

#ifndef _FC_LS_DYNA_IO_P_H_
#define _FC_LS_DYNA_IO_P_H_

#ifdef __cplusplus
extern "C" {
#endif

/**
 * \ingroup LSDynaFileIO
 * \brief Struct hold info about the file format for dataset.
 * 
 * \description
 *
 *    Need to keep word size and endian-ness.
 *
 *    It would be more aesthetically pleasing to put the vertIDs
 *    in each mesh, but if put here, it is easier to manage 
 *    cleaning up the dynamic memory. Ditto for conns_chunks
 *
 * \modifications
 *    - 7/17/2005 WSD Created.
 */
typedef struct {
  int wordlen;       /**< size of units to be read from database */
  int doEndianSwap;  /**< flag for whether to flip bits or not */
  int numTotalVert;   /**< number of global vertices */
  // Coords info is in here since it is really global and not per mesh
  int coordsFileNum;  /**< number postfix for file (don't use if 0). */
  long coords_offset; /**< offset to start of coords */
  long conns8_offset; /**< offset to start of conns for solids (3D) */
  long conns4_offset; /**< offset to start of conns for shells (2D) */
  long conns2_offset; /**< offset to start of conns for beams (1D) */
  int numMesh;        /**< the number of meshes (need to dealloc memory) */
  int** vertIDs;      /**< for each mesh, the global vertex IDs 
			 (index by local ID) */
  int* numElemChunks; /**< for each mesh, the number of chunks that the
			 conns/elem_data are stored in */
  int** startElemPerChunk; /** for each mesh, for each chunk, the offset to
			       the chunk (in # of elements--will usually have 
			       to multiply by # of whatevers for that element 
			       to step through conns/elem_data arrays) */
  int** numElemPerChunk; /** for each mesh, for each chunk, the number
			      of elements in that chunk */
} _FC_DsLSDynaIOInfo;

/**
 * \ingroup LSDynaFileIO
 * \brief Struct hold info about the file format for mesh.
 * 
 * \description
 *
 *    Need to keep offset to topology and global IDs.
 *
 *    NOTE: We also keep the mesh's part or material ID, but this is stored up
 *    a level on the mesh's fileInfo object so that it can be retained when
 *    writing to a different file.
 *
 * \modifications
 *    - 7/18/2005 WSD Created.
 */
typedef struct {
  int meshID;         /**< id to look up vertIDs */
} _FC_MeshLSDynaIOInfo;

/**
 * \ingroup LSDynaFileIO
 * \brief Struct hold info about the file format for variable.
 * 
 * \description
 *
 *    Need to keep filename, offset to data. Also info to decode "striped"
 *    data. (Striped = all data written for 1st element, then all data written
 *    for next element; non striped would have single data written for all
 *    elements, then next data written for all elements).
 *
 *    Set buffer size is o.k. because root names are limited to 75 char.
 *
 * \modifications
 *    - 7/18/2005 WSD Created.
 */
typedef struct {
  int dataFileNum;     /**< number postfix for file (don't use if 0). */
  long data_offset;    /**< offset to start of data */
  // for vertex variables
  int doDispl;         /**< flag for whether to calc displacement */
  // for element variables
  int isStriped;       /**< if non-zero, this is the stride size
			(i.e. the number of data values for each elem) */
  int stripeID;        /**< the ID of the value of interest (e.g. has
			to be between 0 and isStriped-1) */
} _FC_VarLSDynaIOInfo;

// Overall IO setup & take down
//FC_ReturnCode _fc_initLSDyna();
//FC_ReturnCode _fc_finalLSDyna();

// File reader
FC_ReturnCode _fc_loadLSDynaDataset(char *filename, FC_Dataset *dataset);

// Lazy readers for big data
//FC_ReturnCode _fc_readLSDynaSequenceCoords(FC_Sequence sequence);
FC_ReturnCode _fc_readLSDynaMeshCoords(FC_Mesh mesh);
FC_ReturnCode _fc_readLSDynaMeshElemConns(FC_Mesh mesh);
//FC_ReturnCode _fc_readLSDynaSubsetMembers(FC_Subset subset);
FC_ReturnCode _fc_readLSDynaVariableData(FC_Variable variable);

// File writer
// NO file write & never will be

// File close
FC_ReturnCode _fc_closeLSDynaDataset(FC_Dataset dataset);

// Misc
int _fc_hasSharedLSDynaVertices(FC_Dataset dataset);

#ifdef __cplusplus
}
#endif

#endif // _FC_LS_DYNA_IO_P_H_
