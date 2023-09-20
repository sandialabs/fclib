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
 * \file exodusioP.h
 * \brief Private declarations for \ref ExodusFileIO module.
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/modules/exodusioP.h,v $
 * $Revision: 1.11 $
 * $Date: 2006/11/28 19:09:14 $
 */

#ifndef _FC_EXODUS_IO_P_H_
#define _FC_EXODUS_IO_P_H_

#ifdef __cplusplus
extern "C" {
#endif

// Include any dependencies for exodus io
#include "netcdf.h"
#include "exodusII.h"

/**
 * \ingroup ExodusFileIO
 * \brief Struct hold info about the file format for dataset.
 *
 * \description
 *
 *    Need to keep exoid
 *
 *    It would be more aesthetically pleasing to put the vertIDs
 *    in each mesh, but if put here, it is easier to manage 
 *    cleaning up the dynamic memory.
 *
 * \modifications
 *    - 8/01/2005 WSD Created.
 */
typedef struct {
  int cpu_word_size;
  int io_word_size;
  int exoid;
  int numTotalVert; /**< number of global vertices */
  int numMesh;      /**< the number of meshes */
  int** vertIDs;    /**< for each mesh, the global vertex IDs (index by
		       local ID) */
} _FC_DsExodusIOInfo;

/**
 * \ingroup ExodusFileIO
 * \brief Struct hold info about the file format for mesh.
 *
 * \description
 *
 *    Holds indexing information
 *
 *    NOTE: We also keep the mesh's Exodus ID, but this is stored up a level on
 *    the mesh's fileInfo object so that it can be retained when writing to a
 *    different file.
 *
 * \modifications
 *    - 9/01/2005 WSD Created.
 */
typedef struct {
  //int blkid;   /**< exodus ID for the block */
  int meshID;  /**< id to look up vertIDs in dsIOInfo */
  int elemOffset; /**< Offset into global element numbering 
		        (elemID + offset = global) */
} _FC_MeshExodusIOInfo;

/**
 * \ingroup ExodusFileIO
 * \brief Struct hold info about the file format for subset.
 *
 * \description
 *
 *    Need to keep node block id & association.
 *    NOTE: Only side sets are read lazily -- it's easier (spacewise) to 
 *    just do the node sets during initial read.
 *
 * \modifications
 *    - 8/08/2005 WSD Created.
 */
typedef struct {
  int assoc_flag;  /**< 0 = node, 1 = side, 2 = elem */
  int nsssid;      /**< exodus ID for the node set/side/elem set */
} _FC_SubExodusIOInfo;

/**
 * \ingroup ExodusFileIO
 * \brief Struct hold info about the file format for variable
 *
 * \description
 *
 *    Need to keep type of var and it's index. global vars
 *    cannot be read lazily.
 *
 * \modifications
 *    - 9/02/2005 WSD Created.
 *    - 9/09/2007 ACG added attrIndex for basic var. I think
 *                this really could just use the varIndex since
 *                we know otherwise if its a seq or basic var
 *                (note that FC_AT_WHOLE_MESH is not an attr).
 *                Probably could even drop this if didn't
 *                do lazy loading
 */
typedef struct {
  char varType[2]; /**< "n" for nodal, or "e" for element */
  int varIndex;    /**< exodus index (starts from 1) */
  int attrIndex;    /**< exodus index (starts from 1) */
} _FC_VarExodusIOInfo;

// Helper - compare variables for sorting
int _fc_exo_seqVar_cmp(const void* a, const void* b);
int _fc_exo_seqVar_meshID_cmp(const void* a, const void* b);

// Helper - create names for var components
FC_ReturnCode _fc_createExoVarCompName(char* name, int numComp, 
                         FC_MathType mathtype, int compID, char** compName);

// Helper - drill down to some exoinfo quickly
FC_ReturnCode _fc_getMeshExodusGlobalNodalIDsPtr(FC_Mesh mesh, int** LUTable);
FC_ReturnCode _fc_getMeshExodusGlobalElementIDOffset(FC_Mesh mesh, int* offset);

// Overall IO setup & take down
FC_ReturnCode _fc_initExodus(void);
// don't need a final

// File reader
FC_ReturnCode _fc_loadExodusDataset(char *filename, FC_Dataset *dataset);

// Lazy readers for big data
FC_ReturnCode _fc_readExodusSequenceCoords(FC_Sequence sequence);
FC_ReturnCode _fc_readExodusMeshCoords(FC_Mesh mesh);
FC_ReturnCode _fc_readExodusMeshElemConns(FC_Mesh mesh);
FC_ReturnCode _fc_readExodusSubsetMembers(FC_Subset subset);
FC_ReturnCode _fc_readExodusVariableData(FC_Variable variable);
FC_ReturnCode _fc_readExodusOneStepVariableData(FC_Variable variable);
FC_ReturnCode _fc_readExodusAttributeData(FC_Variable variable);

// File writer
FC_ReturnCode _fc_writeExodusDataset(FC_Dataset dataset, char* filename);

// File close
FC_ReturnCode _fc_closeExodusDataset(FC_Dataset dataset);

#ifdef __cplusplus
}
#endif

#endif // _FC_EXODUS_IO_P_H_
