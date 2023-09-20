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
 * \file fileioP.h
 * \brief Private declarations for \ref FileIO module.
 *
 * \modifications
 *  - 10/15/07 ACG removed saf
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/modules/fileioP.h,v $
 * $Revision: 1.27 $ 
 * $Date: 2006/09/19 00:57:58 $
 */

#ifndef _FC_FILE_IO_P_H_
#define _FC_FILE_IO_P_H_

#ifdef __cplusplus
extern "C" {
#endif

#ifdef HAVE_LIBXML2
#include <libxml/tree.h>
#endif

// dataset file types
#ifdef HAVE_EXODUS
#include "exodusioP.h"
#undef TRUE
#endif
#include "lsdynaioP.h"


// Our "Ex" native format - used to boot strap example data
FC_ReturnCode _fc_readVertexFile(char* vertfilename, int* numVert, 
		      int* numDim, double** coords);
FC_ReturnCode _fc_readElementFile(char* elemfilename, int* numVert, 
		      int* numElem, FC_ElementType* elemType, char** name,
		      int** conns);
FC_ReturnCode _fc_readDataFile(char* datafilename, int* numPoint, 
		      int* numComponent, FC_AssociationType *assoc, 
		      FC_MathType* mathType, char** name, float** data);
FC_ReturnCode _fc_makeMeshFromFiles(FC_Dataset dataset, char* vertfile, 
                      char* elemfile, FC_Mesh* mesh);
FC_ReturnCode _fc_makeVariableFromFile(FC_Mesh mesh, char* datafile, 
		      FC_Variable* var);

#ifdef HAVE_LIBXML2
// aux file - under development
FC_ReturnCode _fc_writeAuxFileHeader(FILE* file);
FC_ReturnCode _fc_writeAuxFileFooter(FILE* file);
FC_ReturnCode _fc_writeAuxFileSubset(FILE* file, FC_Subset subset);
FC_ReturnCode _fc_writeAuxFileSubsetCore(FILE* file, char* subsetName, 
		      char* meshName, FC_AssociationType assoc,
		      int maxNumMember, int numMember,int* memberIDs);
FC_ReturnCode _fc_writeAuxFileTear(FILE* file, char* tearName, int numSubset,
		      char** subsetNames, double length);
FC_ReturnCode _fc_importAuxFileXMLDoc(char* filename, xmlDocPtr* doc); 
FC_ReturnCode _fc_getAuxFileSubsets(xmlDocPtr doc, int* numSubset, 
		      char*** subsetNames, char*** meshNames,
		      FC_AssociationType** assocs, int** maxNumMembers, 
                      int** numMembers, int*** memberIDs);
FC_ReturnCode _fc_getAuxFileTears(xmlDocPtr doc, int* numTear, 
                      char*** tearNames, int** numSubsetPerTear,
		      char**** subsetNames, double** lengths);

// PVLookmark BB file 
FC_ReturnCode _fc_writeBBFileHeader(FILE* bb_file);
FC_ReturnCode _fc_writeBBFileFooter(FILE* bb_file);
FC_ReturnCode _fc_writeBBFileBoundingBox(FILE* bb_file, char* name, int stepID,
		      char* comment_str, int numDim,
		      FC_Coords lowers, FC_Coords uppers);
FC_ReturnCode _fc_readBBFile(char* filename, int* numBB, char*** names,
		       int** stepIDs,  char***comments,
		       FC_Coords **lowers, FC_Coords **uppers);
#endif

// Misc helpers 
// Used to read .names file, common to both Exodus & LSDyna (for now)
/**
 * \ingroup PrivateFileIO
 * \brief Struct to hold int/string pairs.
 */
typedef struct {
  int id;
  char* string;
} _FC_IntStringPair;
FC_ReturnCode _fc_readNamesFile(char* datasetfilename, 
		      FC_SortedBlobArray* block_names, 
                      FC_SortedBlobArray* nodeset_names, 
		      FC_SortedBlobArray* sideset_names);
void _fc_freeNamesFileNames(FC_SortedBlobArray* block_names, 
		      FC_SortedBlobArray* nodeset_names, 
		      FC_SortedBlobArray* sideset_names);

/** 
 * \ingroup PrivateFileIO
 * \brief Struct to hold info about the file format for dataset.
 * 
 * \description
 *
 *    Holds information pertinent to file IO for a dataset.
 *
 *    NOTE: only the members of this strut should have dynamic data because
 *    it's easier to do cleanup here (via _fc_closeDataset). The infos for
 *    sequeces meshes, subsets and variables should NOT have dynamic data.
 *
 * \modifications
 *    - 7/14/2005 WSD created.
 */
typedef struct {
#ifdef HAVE_EXODUS
  _FC_DsExodusIOInfo exoInfo; /**< Holds Exodus specific file IO info */
#endif
  _FC_DsLSDynaIOInfo dynaInfo; /**< Holds LSDyna specific file IO info */
} _FC_DsFileIOInfo;

/** 
 * \ingroup PrivateFileIO
 * \brief Struct to hold info about the file format for sequence.
 * 
 * \description
 *
 *    Holds information pertinent to file IO for a sequence.
 *
 *    NOTE: purposefully do not have any dynamically allocated bits
 *    in here or we have to have a special free call (see _FC_DsFileIOInfo).
 *
 * \modifications
 *    - 7/18/2005 WSD created.
 */
typedef struct {
  // LSDyna - nothing
} _FC_SeqFileIOInfo;

/** 
 * \ingroup PrivateFileIO
 * \brief Struct to hold info about the file format for mesh.
 * 
 * \description
 *
 *    Holds information pertinent to file IO for a mesh.
 *
 *    NOTE: purposefully do not have any dynamically allocated bits
 *    in here or we have to have a special free call (see _FC_DsFileIOInfo).
 *
 * \modifications
 *    - 7/18/2005 WSD created.
 */
typedef struct {
  int id;  /**< Since some file formats use an abitrary number instead of
	      name to specify a mesh, this field is used carry that ID along 
	      in case of writing
	      to such formats (to keep ID consistent if rewriting to
	      such a format). Should be 1 or greater. */
#ifdef HAVE_EXODUS
  _FC_MeshExodusIOInfo exoInfo;  /**< Holds Exodus specific file IO info */
#endif
  _FC_MeshLSDynaIOInfo dynaInfo; /**< Holds LSDyna specific file IO info */
} _FC_MeshFileIOInfo;

/** 
 * \ingroup PrivateFileIO
 * \brief Struct to hold info about the file format for subset.
 * 
 * \description
 *
 *    Holds information pertinent to file IO for a subset.
 *
 *    NOTE: purposefully do not have any dynamically allocated bits
 *    in here or we have to have a special free call (see _FC_DsFileIOInfo).
 *
 * \modifications
 *    - 7/18/2005 WSD created.
 */
typedef struct {
#ifdef HAVE_EXODUS
  _FC_SubExodusIOInfo exoInfo;   /**< Holds Exodus specific file IO info */
#endif
  // LSDyna - does not have subsets
} _FC_SubFileIOInfo;

/** 
 * \ingroup PrivateFileIO
 * \brief Struct to hold info about the file format for variable.
 * 
 * \description
 *
 *    Holds information pertinent to file IO for a variable.
 *
 *    NOTE: purposefully do not have any dynamically allocated bits
 *    in here or we have to have a special free call (see _FC_DsFileIOInfo).
 *
 * \modifications
 *    - 7/18/2005 WSD created.
 */
typedef struct {
#ifdef HAVE_EXODUS
  _FC_VarExodusIOInfo exoInfo; /**< Holds Exodus specific file IO info */
#endif
  _FC_VarLSDynaIOInfo dynaInfo; /**< Holds LSDyna specific file IO info */
} _FC_VarFileIOInfo;

// Overall IO setup & take down
FC_ReturnCode _fc_initFileIO(void);
FC_ReturnCode _fc_finalFileIO(void);

// Lazy readers for big data
FC_ReturnCode _fc_readSequenceCoords(FC_Sequence sequence);
FC_ReturnCode _fc_readMeshCoords(FC_Mesh mesh);
FC_ReturnCode _fc_readMeshElemConns(FC_Mesh mesh);
FC_ReturnCode _fc_readSubsetMembers(FC_Subset subset);
FC_ReturnCode _fc_readVariableData(FC_Variable variable);

// helpers for closing dataset
FC_ReturnCode _fc_closeDataset(FC_Dataset dataset);

#ifdef __cplusplus
}
#endif

#endif // _FC_FILE_IO_P_H_
