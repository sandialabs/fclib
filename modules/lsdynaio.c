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
 * \file lsdynaio.c
 * \brief Implementation for \ref LSDynaFileIO module.
 *  
 * $Source: /usr/local/Repositories/fcdmf/fclib/modules/lsdynaio.c,v $
 * $Revision: 1.52 $ 
 * $Date: 2006/11/07 23:49:02 $
 */

// C library includes
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// fc library dependencies
#include "base.h"
#include "storage.h"
#include "library.h"
#include "dataset.h"
#include "sequence.h"
#include "mesh.h"
#include "variable.h"
#include "tableP.h"
#include "libraryP.h"
#include "datasetP.h"
#include "sequenceP.h"
#include "subsetP.h"
#include "meshP.h"
#include "variableP.h"
#include "util.h"

// this module
#include "lsdynaioP.h"

/**
 * \addtogroup LSDynaFileIO
 * \brief LS-Dyna file IO.
 *
 * \description
 * 
 *    FCLib does not provide a writer for LS-Dyna files. 
 *
 *    The file IO for d3plot binary files is based on the documentation
 *    in "LS-DYNA Database Binary Output Files" January 1994 Revised May 2002.
 *    HOWEVER, there were a lot of places where the documentation differed
 *    from actual LS-DYNA files. The ultimate resource should be the code
 *    in lsdynaio.c, as we are reasonable sure of its correctness. We
 *    have also tried to document these differences in comments within the
 *    code, and also below.
 *
 *    This reader is not completely general, and will not work on LS-DYNA d3plot
 *    files with SPH nodes, thick shell elements, fluid materials or CFD
 *    values.
 *
 *    Documentation differences:
 *    - In Control Data section and the State Data section, we have encountered
 *      datasets where the number of global variables (NGLBV) is 1, although 
 *      NGLBV is usually > 7. We have no idea which of the 7 normal global vars
 *      it is.
 *    - In Control Data section and the State Data section, NGLBV = <b>7</b> +
 *      6*(NUMMAT8 + ...) + NUMRW. (7 instead of 6).
 *    - In the Geometry Data section, the 2 node beam element connectivities 
 *      come before the 4 node shell element connectivities, not after.
 *    - In the User Material, Node and Element Identification Numbers section, 
 *      the total section length is increased by 6 +
 *      <b>3</b>*(NUMMAT8+NUMMAT4+NUMMAT2+NUMMATT). (The 3 was missing.)
 *    - In the User Material, Node and Element Identification Numbers section,
 *      some dataset seem to be padded and total section length is less than
 *      NARBS.
 *    - In the State Data section, the documentation implies that the
 *      temperature nodal variable comes first. It actually comes second. The
 *      order of nodal variables is: displacement, temperature, velocity, 
 *      acceleration.
 *
 *    Misc notes:
 *    - We never read the number of rigid wall global vars and instead
 *      calculate them using the formula for the total number of global vars.
 *    - The element data is written as all data for an element, then all
 *      data for the next element (= striped).
 *    - Note that while elements are in order based on element type, the
 *      parts/meshes are not in any particular order. \ref
 *      _fc_calcLSDynaElemChunks() creates "chunk" information that maps
 *      from partIDs back into the element data arrays. 
 *
 * \todo Combine stress tensor scalars into 1 variable.
 */

/**
 * \ingroup LSDynaFileIO
 * \brief Helper function to change Endianness.
 *
 * \description 
 *
 *    Will only handle 4 or 8 bytes.
 */
static void _fc_swapEndian(void* Addr, size_t size) {
  unsigned char* word = (unsigned char*)Addr;
  unsigned char temp[4];
  switch(size) {
  case 4:
    temp[0] = word[0];
    temp[1] = word[1];
    word[0] = word[3];
    word[1] = word[2];
    word[2] = temp[1];
    word[3] = temp[0];
    break;
  case 8:
    temp[0] = word[0];
    temp[1] = word[1];
    temp[2] = word[2];
    temp[3] = word[3];
    word[0] = word[7];
    word[1] = word[6];
    word[2] = word[5];
    word[3] = word[4];
    word[4] = temp[3];
    word[5] = temp[2];
    word[6] = temp[1];
    word[7] = temp[0];
  }
}

/**
 * \ingroup LSDynaFileIO
 * \brief Helper function to extract ints from void buffer.
 *
 * \description 
 *
 *    Will only handle 4 or 8 bytes (int and long long)
 */
static int _fc_getIntFromBuffer(void* ptr, size_t size, int idx) {
  long long *longptr;
  int *intptr;
  switch (size) {
  case 4: intptr =  (int*) ptr;        return intptr[idx];   
  case 8: longptr = (long long*) ptr;  return longptr[idx];
  }
  return -1;
}
    
/**
 * \ingroup LSDynaFileIO
 * \brief Helper function to extract floats from void buffer.
 *
 * \description 
 *
 *    Will only handle 4 or 8 bytes (float and double)
 */
static double _fc_getDoubleFromBuffer(void* ptr, size_t size, int idx) {
  float *floatptr;
  double *doubleptr;
  switch (size) {
  case 4: floatptr =  (float*) ptr;   return floatptr[idx];   
  case 8: doubleptr = (double*) ptr;  return doubleptr[idx];
  }
  return -1;
}

/**
 * \ingroup LSDynaFileIO
 * \brief Helper function to do endianSwap while reading from file.
 *
 * \description 
 *
 *    This function wraps fread and performs the requested read with the
 *    following additional functionalities:
 *      - checks to see if we've reach the end of the current file
 *        and opens the next file if so.
 *      - reorders bits if the reading platform has a different Endianness
 *        than the writing platform did.
 *
 *    The last argument is optional (can pass NULL) and only applies
 *    if an error condition is returned. If the flag is 1, then the routine
 *    could not find a next file (this is how we tell we are done in
 *    the main read routine).
 *   
 * \modifications
 *    - 8/5/2005 WSD Added automatically opening the next file.
 */
static FC_ReturnCode _fc_fread_endianSwap(
   void* ptr,    /**< output - pointer to already allocated buffer to read 
                    into */
   size_t size,  /**< input - the size of each object */
   size_t nmemb, /**< input - the number of objects */ 
   FILE** file,  /**< input/output - the current file to read from */
   char* basename, /**< input - the base file name */
   int* fileNum,   /**< input/output - pointer to the current file number */
   int doEndianSwap, /**< input - flag for whether to flip bits */
   int* noNextFile   /**< output - flag, 1 means error condition is 
                        because there is no next file */
) {
  size_t count, i;
  char filename[1024] = ""; // root names are limited to 75 char
  char test[size];

  // default return
  if (noNextFile)
    *noNextFile = 0;

  // read
  count = fread(ptr, size, nmemb, *file);

  // Take first entry & swap endian so can test for 'end of file' word
  if (count > 0) {
    memcpy(test, ptr, size);
    if (doEndianSwap)
      _fc_swapEndian(&test, size);
  }

  // If an error or see 'end of file' word, try opening another file
  if (count != nmemb || _fc_getDoubleFromBuffer(test, size, 0) == -999999) {
    fclose(*file);
    // check all the way up to 999, 'cause they are allowed to skip #'s
    while(*fileNum < 999) {
      (*fileNum)++;
      if (*fileNum < 100)
	sprintf(filename, "%s%02d", basename, *fileNum);
      else if (*fileNum < 1000)
	sprintf(filename, "%s%03d", basename, *fileNum);
      fc_printfLogMessage("Attempting to open next file '%s'", filename);
      *file = fopen(filename, "r");
      if (*file != NULL)
	break;
    }
    // If fail to find file, return special code.
    if (*file == NULL) { 
      if (noNextFile)
        *noNextFile = 1;
      // no error message on purpose since not always an error
      return FC_ERROR;
    }
    // Reread
    count = fread(ptr, size, nmemb, *file);
    if (count != nmemb) {
      fc_printfErrorMessage("Read failed. Expected count of %d, got %d", 
                            (int)nmemb, (int)count); 
      return FC_FILE_IO_ERROR;
    }
  }

  // swap bits
  if (doEndianSwap) {
    for (i = 0; i < nmemb; i++)
      _fc_swapEndian((char*)ptr + i*size, size);
  }

  return FC_SUCCESS;
}

/**
 * \ingroup LSDynaFileIO
 * \brief Helper to discover how elements are assigned to parts (chunked).
 *
 * \description
 *
 *    Connectivity information is stored per element, but the order of elements
 *    is kinda random. This routine uses the partID stored with each element to
 *    create "chunk" offsets and "chunk" lengths to make grabbing element data
 *    later for a specific part easier. (Chunk = contiguous group of elements
 *    from the same part; a part may have 1 or more chunks).
 *
 *    The results are ordered by partID, not by element type. So even though
 *    the first element is a solid, it might be in part 3; and a shell element
 *    might be in the first part.
 *   
 * \modifications
 *    - 9/2006 created WSD.
 */
static FC_ReturnCode _fc_calcLSDynaElemChunks(
  int numTotalMesh,
  int numTotalVert,
  FC_ElementType elemType,
  int numMesh,
  int numElem,
  FC_ElementType *elemTypes,
  int* numVerts,
  int* numElems,
  int** vertIDs,
  int* numElemChunks,
  int** startElemPerChunk,
  int** numElemPerChunk,
  int wordlen,
  FILE** file,
  char* filename,
  int* fileNum,
  int doEndianSwap
)
{
  // FIX bars don't follow numVertPerElem+1 convention
  // I think the below code is correct, but needs to be checked
  FC_ReturnCode rc;
  int i, j, k, m;
  int prevMeshID, count, meshCount, *vertLU;
  void* temp, *temp2;
  int numVertPerElem, numPlus;
  
  // setup
  numVertPerElem = fc_getElementTypeNumVertex(elemType);
  numPlus = numVertPerElem + 1;
  if (elemType == FC_ET_LINE)
    numPlus = 6; // atypical case, discard values at indices 2-4
 
  // Read the conns plus part IDs
  int* global_conns_plus = (int*)malloc(numPlus*numElem*sizeof(int));
  if (!global_conns_plus) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  if (wordlen == 4) {
    rc = _fc_fread_endianSwap(global_conns_plus, wordlen, 
			     numPlus*numElem, file,
			     filename, fileNum, doEndianSwap, NULL);
    if (rc != FC_SUCCESS) {
      fc_printfErrorMessage("Failed to read element conns");
      return rc;
    }
  }
  else {
    long long *temp_conns_plus = (long long*)malloc(numPlus*
				         numElem*sizeof(long long));
    rc = _fc_fread_endianSwap(temp_conns_plus, wordlen, 
			     numPlus*numElem, file,
			     filename, fileNum, doEndianSwap, NULL);
    if (rc != FC_SUCCESS) {
      fc_printfErrorMessage("Failed to read element conns");
      return rc;
    }
    for (i = 0; i < numPlus*numElem; i++)
      global_conns_plus[i] = temp_conns_plus[i];
    free(temp_conns_plus);
  }
  
  // Cannot assume that all elements for a part will be in 1 continguous block.
  // so, walk global_conns_plus once to determine blocks for each mesh 
  prevMeshID = -1;
  for (i = 0; i < numElem; i++) {
    // Note: MeshID = PartID - 1
    int meshID = global_conns_plus[i*numPlus + numPlus-1] - 1;
    if (meshID != prevMeshID) {
      elemTypes[meshID] = elemType;
      temp = realloc(numElemPerChunk[meshID], 
		     (1+numElemChunks[meshID])*sizeof(int));
      temp2 = realloc(startElemPerChunk[meshID], 
		     (1+numElemChunks[meshID])*sizeof(int));
      if (!temp || !temp2) {
	fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
	return FC_MEMORY_ERROR;
      }
      numElemPerChunk[meshID] = temp;
      startElemPerChunk[meshID] = temp2;
      numElemPerChunk[meshID][numElemChunks[meshID]] = 1;
      startElemPerChunk[meshID][numElemChunks[meshID]] = i;
      numElemChunks[meshID]++;
      prevMeshID = meshID;
    }
    else
      numElemPerChunk[meshID][numElemChunks[meshID]-1]++;
  }

  // this data structure is used anew in each loop below
  vertLU = (int*)malloc(numTotalVert*sizeof(int)); // global -> local
  if (!vertLU) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }

  // creating global -> local, and local to global vert LUTs
  // Also filling in numElems & numVerts (per part)
  meshCount = 0;
  for (i = 0; i < numTotalMesh; i++) { // loop over meshes
    // skip meshes that are not of the current element type
    if (elemTypes[i] != elemType)
      continue;
    meshCount++;
    // examine conns &  calc numElems[i], calc numVerts[i]
    // assign new ID numbers, collect globalIDs
    // initialize global -> local LUT
    for (j = 0; j < numTotalVert; j++) 
      vertLU[j] = -1;
    numElems[i] = 0;
    for (j = 0; j < numElemChunks[i]; j++) { // loop over chunks
      int current_offset = numPlus*startElemPerChunk[i][j];
      numElems[i] += numElemPerChunk[i][j];
      for (k = 0; k < numElemPerChunk[i][j]; k++) { // l. elems in chunks
	for (m = 0; m < numVertPerElem; m++) {
	  // global IDs start with 1, so have to subtract 1 off
	  vertLU[global_conns_plus[current_offset+m]-1] = 0;
	}
	current_offset += numPlus;
      }
    }
    numVerts[i] = 0;
    for (j = 0; j < numTotalVert; j++) {
      if (vertLU[j] > -1) {
        vertLU[j] = numVerts[i]; // new ID
        numVerts[i]++;
      }
    }
    vertIDs[i] = (int*)malloc(numVerts[i]*sizeof(int));
    if (!vertIDs[i]) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      return FC_MEMORY_ERROR;
    }
    count = 0;
    for (j = 0; j < numTotalVert; j++) {
      if (vertLU[j] > -1) {
        vertIDs[i][count] = j;
        count++;
      }
    } 
  }
  free(global_conns_plus);
  free(vertLU);

  // consistency check
  if (meshCount != numMesh) {
    fc_printfErrorMessage("didn't process the expected number of meshes");
    return FC_ERROR;
  }

  return FC_SUCCESS;
}

/**
 * \ingroup  LSDynaFileIO
 * \brief  Load all members of an LS-Dyna Dataset
 *
 * \description
 *
 *    Given a the name of the root LS-Dyna database file, this routine loads
 *    the dataset and makes all of it's members available for use.  All
 *    meta data data is loaded immediately and most big data is loaded 
 *    lazily (i.e. not fetched until needed).
 *
 * \todo To get node sets and side sets, will need to read input deck
 *    at same time ...
 *
 * \modifications  
 *   - 01/10/04 WSK, Created.
 */
FC_ReturnCode _fc_loadLSDynaDataset(
 char *filename,       /**< input - name of LS-Dyna file to open */
 FC_Dataset *dataset   /**< output - handle for this dataset */
) {
  FC_ReturnCode rc;
  int i, j;
  _FC_DsSlot *dsSlot;
  _FC_MeshSlot* meshSlot;
  _FC_VarSlot* varSlot;
  FC_Sequence sequence;
  FC_Mesh* meshes;
  FC_Variable **vars, *seqVar;
  char globalVarNames[5][20] = { "global_KE", "global_IE", "global_TE", 
                                  "global_velocity", "global_work" };
  int globalVarNumComps[5] = { 1, 1, 1, 3, 1 };
  char perMaterialVarNames[4][15] = { "mat_IE", "mat_KE", "mat_velocity",
                                      "mass" };
  int perMatVarComps[4] = { 1, 1, 3, 1 };                         
  char tempVarName[12] =  "temperature";
  char displVarName[13] = "displacement";
  char velVarName[9] =    "velocity";
  char accelVarName[13] = "acceleration";
  char elem8VarNames[7][15] = { "sigma-xx", "sigma-yy", "sigma-zz",
                               "sigma-xy", "sigma-yz", "sigma-zx", 
                               "plastic_strain" };
  // FIX! elem4 meshes have two variables named IE (prev one is FC_AT_WHOLE_MESH)
  char elem4VarNames[24][15] = { "mxx", "myy", "mxy", "qx", "qy", "nxx", "nyy",
				 "nxy", "thickness", "dv01", "dv02",
				 "in-eps-xx", "in-eps-yy", "in-eps-zz",
				 "in-eps-xy", "in-eps-yz", "in-eps-zx",
				 "in-eps-xx", "in-eps-yy", "in-eps-zz",
				 "in-eps-xy", "in-eps-yz", "in-eps-zx",
				 "IE" };
  char elem2VarNames[11][15] = { "axial_force", "s-shear", "t-shear",
				 "s-bending", "t-bending", "torsion",
				 "rs-shear_stress", "tr-shear_stress", 
				 "axial-stress", "plastic_strain",
				 "axial_strain" };
  char elemDeathVarName[11] = { "elem_death" };
  FILE* file;
  _FC_IntStringPair *blob;
  FC_SortedBlobArray block_names, nodeset_names, sideset_names;
  int wordlen, doEndianSwap; // file storage structure
  int doTempVar, doDisplVar, doVelVar, doAccelVar; // flags for vert vars
  int doElemDeathVar; // flags for elem vars
  int numArbWord; // space for arbitrary node & element ordering section
  int fileNum;  // number postfix for file
  int numDim, numTotalVert, numTotalElem8, numTotalElem2, numTotalElem4;
  int *numVerts, *numElems;
  int **vertIDs;
  int NEIPH; // # of additional values per solid element (elem8)
  int NEIPS; // # of additional values per integration point for shells (elem4)
  int MAXINT; // used to calc # of shell variables
  int ISTRN; // used to calc # of shell variables - i.e. are strains present
  int IOSHL1, IOSHL2, IOSHL3, IOSHL4; // used to calc # of shell variables
  int numMesh8, numMesh2, numMesh4, numTotalMesh;
  int numTotalGlbVar, numGlbVar, numRW, numPerMatGlbVar, numSeqVar;
  int numElemVar8, numElemVar2, numElemVar4;
  int numStep, maxStep;
  char tempName[1024];
  char title[8*10+1];
  unsigned char wordBuf[8*64];
  int intBuf;
  void *time_coords;
  FC_ElementType *elemTypes;
  int seqVarID;
  FC_DataType datatype;
  int *numElemChunks;
  int** startElemPerChunk;
  int** numElemPerChunk;

  // Note: sometimes part and material are used interchangeably here,
  // But partID is the order/indexing in the LS-DYNA file, indexing
  // from 1 (therefore meshID is partID-1), the material ID is
  // an arbitrary number. The materials section gives a mapping from
  // partID to matID.

  // Assumption: each part will be of only 1 element type, i.e.
  // there will never be a part with some elem8 and some elem4.
  // I don't know if this is a valid assumption.

  // default return value
  if (dataset != NULL)
    *dataset = FC_NULL_DATASET;
  
  // check input
  if ((filename == NULL) || (dataset == NULL)) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // make sure library is ready, return error if it's not
  if (!_fc_isLibraryInit()) {
    fc_printfErrorMessage("Cannot load a dataset until library is "
                          "initialized.");
    return FC_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Loading dataset '%s'", filename);

  //--- Set up

  // open first file
  fileNum = 0;
  file = fopen(filename, "r");
  if (file == NULL) {
    fc_printfErrorMessage("couldn't open file '%s'", filename);
    return FC_FILE_IO_ERROR;
  }
  
  // test byte order by checking that a specific value, ICODE, makes sense
  wordlen = 4;
  doEndianSwap = 0;
  fseek(file, 17*wordlen, SEEK_SET);
  rc = _fc_fread_endianSwap(wordBuf, wordlen, 1, &file, filename, &fileNum, 
                           doEndianSwap, NULL);
  intBuf = _fc_getIntFromBuffer(wordBuf, wordlen, 0);
  if (rc != FC_SUCCESS || !(intBuf == 2 || intBuf == 6)) {
    doEndianSwap = 1;
    fseek(file, 17*wordlen, SEEK_SET);
    rc = _fc_fread_endianSwap(&wordBuf, wordlen, 1, &file, filename, &fileNum, 
                             doEndianSwap, NULL);
    intBuf = _fc_getIntFromBuffer(wordBuf, wordlen, 0);
    if (rc != FC_SUCCESS || !(intBuf == 2 || intBuf == 6)) {
      wordlen = 8;
      doEndianSwap = 0;
      fseek(file, 17*wordlen, SEEK_SET);
      rc = _fc_fread_endianSwap(&wordBuf, wordlen, 1, &file, filename, &fileNum,
                               doEndianSwap, NULL);
      intBuf = _fc_getIntFromBuffer(wordBuf, wordlen, 0);
      if (rc != FC_SUCCESS || !(intBuf == 2 || intBuf == 6)) {
        doEndianSwap = 1;
        fseek(file, 17*wordlen, SEEK_SET);
        rc = _fc_fread_endianSwap(&wordBuf, wordlen, 1, &file, filename, 
                                 &fileNum, doEndianSwap, NULL);
        intBuf = _fc_getIntFromBuffer(wordBuf, wordlen, 0);
        if (rc != FC_SUCCESS || !(intBuf == 2 || intBuf == 6)) {
          fc_printfErrorMessage("Does not look like an LS-Dyna database");
          fclose(file);
          return FC_FILE_IO_ERROR;
        }
      }
    }
  }

  if (wordlen == 4)
    datatype = FC_DT_FLOAT;
  else
    datatype = FC_DT_DOUBLE;

  // NOTE: Uppercase variable names are LS-Dyna's names for the values
  // See the manual "LS-DYNA Database Binary Output Files"
  // Binary file is read in words, so we can seek right to point we want

  //--- Read the title

  // don't need to do endian stuff since chars 1 byte
  fseek(file, 0, SEEK_SET);
  fread(title, sizeof(char), 10*wordlen, file);
  title[10*wordlen] = '\0'; // make sure string is terminated
  // clip trailing spaces
  for (i = 10*wordlen-1; i >= 0; i--) {
    if (title[i] == ' ')
      title[i] = '\0';
    else
      break;
  }

  //--- Read entire control section (guaranteed to be in 1st file)

  fseek(file, 0, SEEK_SET);
  rc = _fc_fread_endianSwap(&wordBuf, wordlen, 64, &file, filename, &fileNum, 
                           doEndianSwap, NULL);
  if (rc != FC_SUCCESS) {
    fc_printfErrorMessage("Failed to read control section");
    fclose(file);
    return rc;
  }

  // NDIM (word 15) is combo flag for number of dimensions and MATTYPE
  // (NDIM == 4 means MATTYPE = 0)
  if (_fc_getIntFromBuffer(wordBuf, wordlen, 15) != 4) {
    fc_printfErrorMessage("Cannot read dataset unless it is 3D and has "
                          "unpacked ints (i.e. LS-Dyna parameter NDIM=4 "
                          "not %d)", _fc_getIntFromBuffer(wordBuf, wordlen, 15));
    fclose(file);
    return FC_ERROR;
  }
  numDim = 3;
  // NUMNP (word 16) - numTotalVert
  numTotalVert = _fc_getIntFromBuffer(wordBuf, wordlen, 16);
  // ICODE (word 17) - flag for type of generating code
  if (_fc_getIntFromBuffer(wordBuf, wordlen, 17) != 6) {
    fc_printfErrorMessage("This is not a NIKE3D, LS-DYNA3D or LS-NIKE3D "
                          "dataset (expected LS-Dyna parameter ICODE=6 not "
                          "%d)\n", _fc_getIntFromBuffer(wordBuf, wordlen, 17));
    fclose(file);
    return FC_ERROR;
  }
  // NGLBV (word 18) - number of global variables
  // (NOTE: numTotalGlbVar == 1 is an undocumented special case)
  numTotalGlbVar = _fc_getIntFromBuffer(wordBuf, wordlen, 18);
  if (numTotalGlbVar < 7 && numTotalGlbVar != 1) {
    fc_printfErrorMessage("Expected NGLBV (%d) to be either 1 or >= 7",
			  numTotalGlbVar); 
    return FC_ERROR;
  } 
  // IT (word 19) - flag for temperature variable
  doTempVar = _fc_getIntFromBuffer(wordBuf, wordlen, 19);
  // IU (word 20) - flag for current geometry
  doDisplVar = _fc_getIntFromBuffer(wordBuf, wordlen, 20);
  // IV (word 21) - flag for velocities
  doVelVar = _fc_getIntFromBuffer(wordBuf, wordlen, 21);
  // IA (word 22) - flag for accelerations
  doAccelVar = _fc_getIntFromBuffer(wordBuf, wordlen, 22);
  // NEL8 (word 23) - number of hexes, numTotalElem8
  numTotalElem8 = _fc_getIntFromBuffer(wordBuf, wordlen, 23);
  // NUMMAT8 (word 24) - number of materials used by elements = numMesh8
  // (a mesh per material)
  numMesh8 = _fc_getIntFromBuffer(wordBuf, wordlen, 24);
  // (word 25 & 26 are not used by dyna)
  // NV3D (word 27) - number of 3D element vars
  numElemVar8 = _fc_getIntFromBuffer(wordBuf, wordlen, 27);
  // NEL2 (word 28) - number of line elements
  numTotalElem2 = _fc_getIntFromBuffer(wordBuf, wordlen, 28);
  // NUMMAT2 (word 29) - number of materials used by 1D elements = numMesh2
  // (a mesh per material)
  numMesh2 = _fc_getIntFromBuffer(wordBuf, wordlen, 29);
  // NV1D (word 27) - number of 1D element vars
  numElemVar2 = _fc_getIntFromBuffer(wordBuf, wordlen, 30);
  // NEL4 (word 31) - number of quad elements
  numTotalElem4 = _fc_getIntFromBuffer(wordBuf, wordlen, 31);
  // NUMMAT4 (word 32) - number of materials used by 2D elements = numMesh4
  // (a mesh per material)
  numMesh4 = _fc_getIntFromBuffer(wordBuf, wordlen, 32);
  // NV2D (word 27) - number of 2D element vars
  numElemVar4 = _fc_getIntFromBuffer(wordBuf, wordlen, 33);
  // NEIPH (word 34) - number of additional vars for hex elements
  NEIPH = _fc_getIntFromBuffer(wordBuf, wordlen, 34);
  if (numElemVar8 > 0 && NEIPH != numElemVar8 - 7) {
    fc_printfErrorMessage("Unexpected NEIPH value (expected NEIPH "
                          "(=%d) + 7 to equal  NV3D (=%d)\n)", 
                          NEIPH, numElemVar8);
    return FC_ERROR;
  }
  // NEIPS (word 35) - number of additional values for integrations points
  NEIPS = _fc_getIntFromBuffer(wordBuf, wordlen, 35);
  // FIX? we could handle these ...
  if (_fc_getIntFromBuffer(wordBuf, wordlen, 35) != 0) {
    fc_printfErrorMessage("Expected 0 not %d 1D additional values for "
                          "integration points on shell elements (i.e. LS-Dyna"
                          " parameter NEIPS=0)\n", 
                          _fc_getIntFromBuffer(wordBuf, wordlen, 35));
    return FC_ERROR;
  }
  // MAXINT (word 36) - flag for MAXINT (number of integration points for shells)
  // and MDLOPT (element deletion table)
  MAXINT = _fc_getIntFromBuffer(wordBuf, wordlen, 36);
  if (MAXINT >= 0) {
    doElemDeathVar = 0;
  }
  else if (MAXINT < -10000) {
    doElemDeathVar = 1;
    MAXINT = abs(MAXINT) - 10000;
  }
  else {
    fc_printfErrorMessage("Expected MDLOPT to be 0 or 2, not 1 (can't handle "
			  "node deletion)");
    return FC_ERROR;
  }
  // NMSPH (word 37) - number of SPH nodes
  if (_fc_getIntFromBuffer(wordBuf, wordlen, 37) != 0) {
    fc_printfErrorMessage("Expected 0 not %d SPH nodes (i.e. LS-Dyna "
                          "parameter NMSPH=0)\n", 
                          _fc_getIntFromBuffer(wordBuf, wordlen, 37));
    return FC_ERROR;
  }
  // NGPSPH (word 38) - number of SPH materials
  if (_fc_getIntFromBuffer(wordBuf, wordlen, 38) != 0) {
    fc_printfErrorMessage("Expected 0 not %d SPH materials (i.e. LS-Dyna "
                          "parameter NGPSPH=0)\n",  
                          _fc_getIntFromBuffer(wordBuf, wordlen, 38));
    return FC_ERROR;
  }
  // NARBS (word 39) - additional storage space for arbitrary node and
  // element numbering
  numArbWord = _fc_getIntFromBuffer(wordBuf, wordlen, 39);
  // NELT (word 40) - number of 8 node thick shell elements
  if (_fc_getIntFromBuffer(wordBuf, wordlen, 40) != 0) {
    fc_printfErrorMessage("Expected 0 not %d thick shell elems (i.e. LS-Dyna "
                          "parameter NELT=0)\n",  
                          _fc_getIntFromBuffer(wordBuf, wordlen, 40));
    return FC_ERROR;
  }
  // (word 41 & 42 are about thick shell elements and are skipped) 
  // IOSHL(1) (word 43) - 6 stress components flag
  if (_fc_getIntFromBuffer(wordBuf, wordlen, 43) == 1000)
    IOSHL1 = 1;
  else
    IOSHL1 = 0;
  // IOSHL(2) (word 44) - 6 stress components flag
  if (_fc_getIntFromBuffer(wordBuf, wordlen, 44) == 1000)
    IOSHL2 = 1;
  else
    IOSHL2 = 0;
  // IOSHL(3) (word 45) - 6 stress components flag
  if (_fc_getIntFromBuffer(wordBuf, wordlen, 45) == 1000)
    IOSHL3 = 1;
  else
    IOSHL3 = 0;
  // IOSHL(4) (word 46) - 6 stress components flag
  if (_fc_getIntFromBuffer(wordBuf, wordlen, 46) == 1000)
    IOSHL4 = 1;
  else
    IOSHL4 = 0;
  // IALEMAT (word 47) - number of fluid material ids
  if (_fc_getIntFromBuffer(wordBuf, wordlen, 47) != 0) {
    fc_printfErrorMessage("Expected 0 not %d fluid materials (i.e. LS-Dyna "
                          "parameter IALEMAT=0)\n",  
                          _fc_getIntFromBuffer(wordBuf, wordlen, 47));
    return FC_ERROR;
  }
  // NCFDV1 (word 48) - bit flags for CFD nodal values
  if (_fc_getIntFromBuffer(wordBuf, wordlen, 48) != 0) {
    fc_printfErrorMessage("Expected 0 not %d bit flag for CFD values (i.e. "
                          "LS-Dyna parameter NCFDV1=0)\n",  
                          _fc_getIntFromBuffer(wordBuf, wordlen, 48));
    return FC_ERROR;
  }
  // NCFDV2 (word 49) - more bit flags for CFD nodal values
  if (_fc_getIntFromBuffer(wordBuf, wordlen,49) != 0) {
    fc_printfErrorMessage("Expected 0 not %d more bit flag for CFD values "
                          "(i.e. LS-Dyna parameter NCFDV2=0)\n",  
                          _fc_getIntFromBuffer(wordBuf, wordlen, 49));
    return FC_ERROR;
  }
  // (word 50-63 are not used by ls-dyna)

  // Checks & some more processing
  numTotalMesh = numMesh8 + numMesh2 + numMesh4;
  if (numTotalMesh < 1) {
    fc_printfErrorMessage("Expected at least one solid, shell or beam "
			  "element");
    return FC_ERROR;
  }
  if (numTotalElem4 > 0 && numElemVar4 -
      (MAXINT*(6*IOSHL1 + IOSHL2 + NEIPS) + 8*IOSHL3 + 4*IOSHL4) > 1) {
    ISTRN = 1;
  }
  else {
    ISTRN = 0;
  }

  //--- Create new dataset

  rc = fc_createDataset(title, dataset);
  if (rc != FC_SUCCESS) {
    fc_printfErrorMessage("Failed to create dataset");
    return rc;
  }
  // twiddle writable flag so we don't accidently clobber the original file
  dsSlot = _fc_getDsSlot(*dataset);
  dsSlot->header.committed = 1;
  dsSlot->fileType = FC_FT_LSDYNA;
  dsSlot->baseFile = (char*)malloc((strlen(filename)+1)*sizeof(char));
  if (!dsSlot->baseFile) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  strcpy(dsSlot->baseFile, filename);
  dsSlot->fileInfo.dynaInfo.wordlen = wordlen;  // save for later
  dsSlot->fileInfo.dynaInfo.doEndianSwap = doEndianSwap; // save for later
  dsSlot->fileInfo.dynaInfo.numTotalVert = numTotalVert; 
  dsSlot->writable = 0;

  //--- Read geometry section (assumed to all be in 1st file) & create meshes

  // Do coords

  // Save coords start location
  dsSlot->fileInfo.dynaInfo.coordsFileNum = fileNum;
  dsSlot->fileInfo.dynaInfo.coords_offset = ftell(file);
  fseek(file, wordlen*numDim*numTotalVert, SEEK_CUR);

  // Do conns
  // (NOTE: The LS-DYNA documentation says the 2D elements are before
  // the 1D elements, but in practice we have seen the 2D elements
  // after the 1D elements.)
 
  // Setup verts and elemChunks arrays, initialize elem chunks to zero
  elemTypes = (FC_ElementType*)calloc(numTotalMesh,sizeof(FC_ElementType));
  numVerts = (int*)malloc(numTotalMesh*sizeof(int)); // numVertPerMesh
  numElems = (int*)malloc(numTotalMesh*sizeof(int)); // numElemPerMesh
  vertIDs = (int**)malloc(numTotalMesh*sizeof(int*)); // local -> global
  numElemChunks = (int*)calloc(numTotalMesh, sizeof(int));
  startElemPerChunk= (int**)calloc(numTotalMesh, sizeof(int*));
  numElemPerChunk = (int**)calloc(numTotalMesh, sizeof(int*));
  if (!elemTypes || !numVerts || !numElems || !vertIDs || 
      !numElemChunks || !numElemPerChunk || !startElemPerChunk) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    free(elemTypes);
    free(numVerts);
    free(numElems);
    free(vertIDs);
    free(numElemChunks);
    free(startElemPerChunk);
    free(numElemPerChunk);
    return FC_MEMORY_ERROR;
  }

  // The same chunks data structure is used in all 3 calls - collects all info
  // Handle 3D elements (solids)
  if (numMesh8 > 0) {
    dsSlot->fileInfo.dynaInfo.conns8_offset = ftell(file);
    _fc_calcLSDynaElemChunks(numTotalMesh, numTotalVert, FC_ET_HEX, numMesh8,
                             numTotalElem8, elemTypes, numVerts, numElems,
                             vertIDs, numElemChunks, startElemPerChunk,
                             numElemPerChunk, wordlen, &file, filename,
                             &fileNum, doEndianSwap);
  }
  // Handle 1D elements (bars)
  if (numMesh2 > 0) {
    dsSlot->fileInfo.dynaInfo.conns2_offset = ftell(file);
    _fc_calcLSDynaElemChunks(numTotalMesh, numTotalVert, FC_ET_LINE, numMesh2,
                             numTotalElem2, elemTypes, numVerts, numElems,
                             vertIDs, numElemChunks, startElemPerChunk,
                             numElemPerChunk, wordlen, &file, filename,
                             &fileNum, doEndianSwap);
  }
  // Handle 2D elements (shells)
  if (numMesh4 > 0) {
    dsSlot->fileInfo.dynaInfo.conns4_offset = ftell(file);
    _fc_calcLSDynaElemChunks(numTotalMesh, numTotalVert, FC_ET_QUAD, numMesh4,
                             numTotalElem4, elemTypes, numVerts, numElems,
                             vertIDs, numElemChunks, startElemPerChunk,
                             numElemPerChunk, wordlen, &file, filename,
                             &fileNum, doEndianSwap);
  }
  
  // Make the meshes
  meshes = (FC_Mesh*)malloc(numTotalMesh*sizeof(FC_Mesh));
  if (!meshes) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  for (i = 0; i < numTotalMesh; i++) {
    //if (elemTypes[i] == 0) {
    //printf("Skipping something\n");
    //meshes[i] = FC_NULL_MESH;
    //continue;
    //}
    sprintf(tempName, "part_%d", i+1); 
    rc = fc_createMesh(*dataset, tempName, &meshes[i]);
    if (rc != FC_SUCCESS) {
      fc_printfErrorMessage("Failed to create mesh #%d '%s'", i, tempName);
      return rc;
    }
    meshSlot = _fc_getMeshSlot(meshes[i]);
    meshSlot->fileInfo.id = i+1;
    meshSlot->fileInfo.dynaInfo.meshID = i;
    meshSlot->dim = numDim;
    meshSlot->topodim = fc_getElementTypeTopoDim(elemTypes[i]);
    meshSlot->numVertex = numVerts[i];
    meshSlot->numElement = numElems[i];
    meshSlot->elemType = elemTypes[i];
    meshSlot->header.committed = 1;
  }

  // save necessary info to dataset
  dsSlot->fileInfo.dynaInfo.numMesh = numTotalMesh;
  dsSlot->fileInfo.dynaInfo.vertIDs = vertIDs;
  dsSlot->fileInfo.dynaInfo.numElemChunks = numElemChunks;
  dsSlot->fileInfo.dynaInfo.numElemPerChunk = numElemPerChunk;
  dsSlot->fileInfo.dynaInfo.startElemPerChunk = startElemPerChunk;
  // free other stuff except elemtypes
  free(numVerts);
  free(numElems);

  //--- Arbitrary IDs section (assumed to be in 1st file)

  //printf("***numArbWord = %d\n", numArbWord);
  if (numArbWord > 0) {
    int NSORT, doMatNums;
    int temp_count = 10 + numTotalVert + numTotalElem8 +
      numTotalElem2 + numTotalElem4; // to help figure out padding at end

    // Read first 10 constants
    rc = _fc_fread_endianSwap(&wordBuf, wordlen, 10, &file, filename, &fileNum,
			      doEndianSwap, NULL);
    if (rc != FC_SUCCESS) {
      fc_printfErrorMessage("Failed to read user numbers constants");
      return rc;
    }
    // Check NSORT - if negative, we do material numbers too
    NSORT = _fc_getIntFromBuffer(wordBuf, wordlen, 0);
    if (NSORT < 0)
      doMatNums = 1;
    else
      doMatNums = 0;
    // Read other constants for mats
    if (doMatNums) {
      rc = _fc_fread_endianSwap(&wordBuf, wordlen, 6, &file, filename,
				&fileNum, doEndianSwap, NULL);
      if (rc != FC_SUCCESS) {
	fc_printfErrorMessage("Failed to read user numbers constants");
	return rc;
      }
      // (NOTE: In the LS-DYNA documentation the 3 was missing)
      temp_count += 6 + 3*(numMesh8 + numMesh2 + numMesh4);
    }
    
    // user defined numbers for nodes
    for (i = 0; i < numTotalMesh; i++) {
      FC_Variable temp_var;
      rc = fc_createVariable(meshes[i], "vert_IDs", &temp_var);
      if (rc != FC_SUCCESS) {
	fc_printfErrorMessage("failed to create nodal ID var on mesh #%d",
			      i);
	return rc;
      }
      varSlot = _fc_getVarSlot(temp_var);
      meshSlot = _fc_getMeshSlot(meshes[i]);
      varSlot->numDataPoint = meshSlot->numVertex;
      varSlot->numComponent = 1;
      varSlot->assoc = FC_AT_VERTEX;
      varSlot->mathtype = FC_MT_SCALAR;
      varSlot->datatype = FC_DT_INT;
      varSlot->header.committed = 1;
      varSlot->fileInfo.dynaInfo.dataFileNum = fileNum;
      varSlot->fileInfo.dynaInfo.data_offset = ftell(file);
      varSlot->fileInfo.dynaInfo.doDispl = 0;
      varSlot->fileInfo.dynaInfo.isStriped = 0;
    }
    fseek(file, numTotalVert*wordlen, SEEK_CUR);
    // user defined numbers for solid elements 
    for (i = 0; i < numTotalMesh; i++) {
      FC_Variable temp_var;
      if (elemTypes[i] != FC_ET_HEX)
	continue;
      rc = fc_createVariable(meshes[i], "elem_IDs", &temp_var);
      if (rc != FC_SUCCESS) {
	fc_printfErrorMessage("failed to create element ID var on mesh #%d",
			      i);
	return rc;
      }
      varSlot = _fc_getVarSlot(temp_var);
      meshSlot = _fc_getMeshSlot(meshes[i]);
      varSlot->numDataPoint = meshSlot->numElement;
      varSlot->numComponent = 1;
      varSlot->assoc = FC_AT_ELEMENT;
      varSlot->mathtype = FC_MT_SCALAR;
      varSlot->datatype = FC_DT_INT;
      varSlot->header.committed = 1;
      varSlot->fileInfo.dynaInfo.dataFileNum = fileNum;
      varSlot->fileInfo.dynaInfo.data_offset = ftell(file);
      varSlot->fileInfo.dynaInfo.isStriped = 0;
    }
    fseek(file, numTotalElem8*wordlen, SEEK_CUR);
    // user defined numbers beam elements
    for (i = 0; i < numTotalMesh; i++) {
      FC_Variable temp_var;
      if (elemTypes[i] != FC_ET_LINE)
	continue;
      rc = fc_createVariable(meshes[i], "elem_IDs", &temp_var);
      if (rc != FC_SUCCESS) {
	fc_printfErrorMessage("failed to create element ID var on mesh #%d",
			      i);
	return rc;
      }
      varSlot = _fc_getVarSlot(temp_var);
      meshSlot = _fc_getMeshSlot(meshes[i]);
      varSlot->numDataPoint = meshSlot->numElement;
      varSlot->numComponent = 1;
      varSlot->assoc = FC_AT_ELEMENT;
      varSlot->mathtype = FC_MT_SCALAR;
      varSlot->datatype = FC_DT_INT;
      varSlot->header.committed = 1;
      varSlot->fileInfo.dynaInfo.dataFileNum = fileNum;
      varSlot->fileInfo.dynaInfo.data_offset = ftell(file);
      varSlot->fileInfo.dynaInfo.isStriped = 0;
    }
    fseek(file, numTotalElem2*wordlen, SEEK_CUR);
    // user defined numbers for shell elements
    for (i = 0; i < numTotalMesh; i++) {
      FC_Variable temp_var;
      if (elemTypes[i] != FC_ET_QUAD)
	continue;
      rc = fc_createVariable(meshes[i], "elem_IDs", &temp_var);
      if (rc != FC_SUCCESS) {
	fc_printfErrorMessage("failed to create element ID var on mesh #%d",
			      i);
	return rc;
      }
      varSlot = _fc_getVarSlot(temp_var);
      meshSlot = _fc_getMeshSlot(meshes[i]);
      varSlot->numDataPoint = meshSlot->numElement;
      varSlot->numComponent = 1;
      varSlot->assoc = FC_AT_ELEMENT;
      varSlot->mathtype = FC_MT_SCALAR;
      varSlot->datatype = FC_DT_INT;
      varSlot->header.committed = 1;
      varSlot->fileInfo.dynaInfo.dataFileNum = fileNum;
      varSlot->fileInfo.dynaInfo.data_offset = ftell(file);
      varSlot->fileInfo.dynaInfo.isStriped = 0;
    }
    fseek(file, numTotalElem4*wordlen, SEEK_CUR);

    // get mat IDs and use in mesh names
    if (doMatNums) {
      int numMat;
      void* tempWordBuf;
      
      // check NMMAT
      numMat =_fc_getIntFromBuffer(wordBuf, wordlen, 5);
      if (numMat != numTotalMesh) {
	fc_printfErrorMessage("Expected numMat (%d) to equal numMesh8 (%d) + "
			      "numMesh2 (%d) + numMesh4 (%d)", numMat, numMesh8, 
			      numMesh2, numMesh4);
	return FC_ERROR;
      }
      // read & use material IDs
      tempWordBuf = malloc(numMat*wordlen);
      if (tempWordBuf == NULL) {
	fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
	return FC_MEMORY_ERROR;
      }
      rc = _fc_fread_endianSwap(tempWordBuf, wordlen, numMat, &file, filename, 
				&fileNum, doEndianSwap, NULL);
      for (i = 0; i < numTotalMesh; i++) {
	if (FC_HANDLE_EQUIV(meshes[i], FC_NULL_MESH))
	  continue;
	meshSlot = _fc_getMeshSlot(meshes[i]);
	meshSlot->fileInfo.id = _fc_getIntFromBuffer(tempWordBuf, wordlen, i);
	sprintf(tempName, "Material_ID_%d", meshSlot->fileInfo.id);
	_fc_setSlotHeaderName(&meshSlot->header, tempName);
      }
      free(tempWordBuf);

      // skip some more stuff
      fseek(file, 2*numMat*wordlen, SEEK_CUR);
    }

    // NOTE: For some reason, some datasets seem to be padded ... skip padding
    if (temp_count != numArbWord)
      fseek(file, (numArbWord - temp_count)*wordlen, SEEK_CUR);
  }

  //--- Read the names file, if it exists, & rename meshes
  rc = _fc_readNamesFile(filename, &block_names, &nodeset_names, &sideset_names);
  if (rc != FC_SUCCESS)
    return rc;
  if (block_names.numBlob > 0) {
    // rename meshes
    for (i = 0; i < numTotalMesh; i++) {
      int key;
      if (FC_HANDLE_EQUIV(meshes[i], FC_NULL_MESH))
	continue;
      meshSlot = _fc_getMeshSlot(meshes[i]);
      if (!strncmp(meshSlot->header.name, "part", 4)) 
	key = atoi(&meshSlot->header.name[5]);
      else
	key = atoi(&meshSlot->header.name[12]);
      if (fc_isBlobInSortedBlobArray(&block_names, &key, fc_intCompare,
				     (void**)&blob)) {
	_fc_setSlotHeaderName(&meshSlot->header, blob->string);
      }
      // else - keep current name
    }
  }
  _fc_freeNamesFileNames(&block_names, &nodeset_names, &sideset_names);

  //--- Read state data section (might be in multiple files, but assume
  //    first state is in first file)

  // (Assume all of a single state will be in the same file).

  // NOTE: because the State Data will span multiple files, it's difficult
  // to figure out the total number of states a priori, so we'll just have
  // to grow arrays of things instead of a priori allocation.

  // NOTE: There is a special case not documented in LS-DYNA documentation,
  // where numTotalGlbVar = 1. In that case, there are no per material global 
  // vars and no per rigid wall global vars.
  // (p.s. some day might want to have a flag other than numTotalGlbVar)
  if (numTotalGlbVar == 1) {
    numGlbVar = numTotalGlbVar;
    numPerMatGlbVar = 0;
    numRW = 0;
  }
  else {
    // (NOTE: The LS-DYNA documentation has 6 global vars but there are 7)
    numGlbVar = 5; // would be 7, but we combine 3 velocity comps
    numPerMatGlbVar = 4; // would be 6, but we combine 3 vel comps
    numRW = numTotalGlbVar - 7 - 6*numTotalMesh;
  }

  // Total number of vars we need to turn into seqvars
  // Store vars as we make them, turn into seqvars at very end
  numSeqVar = numGlbVar + numRW 
              + numTotalMesh*(numPerMatGlbVar + doTempVar + doDisplVar +
                              doVelVar + doAccelVar + doElemDeathVar)
              + numMesh8*numElemVar8 
              + numMesh4*numElemVar4
              + numMesh2*numElemVar2;
  vars = (FC_Variable**)malloc(numSeqVar*sizeof(FC_Variable*));
  if (!vars) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  for (i = 0; i < numSeqVar; i++) {
    vars[i] = (FC_Variable*)malloc(sizeof(FC_Variable));
    if (vars[i] == NULL) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      return FC_MEMORY_ERROR;
    }
  }

  numStep = 0;
  maxStep = 1;
  time_coords = (void*)malloc(wordlen);
  if (!time_coords) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  while (1) {
    int noNextFile;

    // read time coord value,
    rc = _fc_fread_endianSwap(&wordBuf, wordlen, 1, &file, filename, &fileNum,
                               doEndianSwap,&noNextFile);
    if (rc != FC_SUCCESS) {
      if (noNextFile)
        break; // no more files, we're done! (the file was closed)
      else {
        fc_printfErrorMessage("Failed to read time coord");
        fclose(file);
        return FC_FILE_IO_ERROR; // error => exit function
      }
    }

    // expand arrays if needed
    if (numStep >= maxStep) {
      void* temp_p;
      FC_Variable* temp_vars;
      maxStep *= 2;
      temp_p = realloc(time_coords, maxStep*wordlen);
      if (temp_p == NULL) {
        fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
        return FC_MEMORY_ERROR;
      }
      time_coords = temp_p;
      for (i = 0; i < numSeqVar; i++) {
        temp_vars = (FC_Variable*)realloc(vars[i],maxStep*sizeof(FC_Variable));
        if (temp_vars == NULL) {
          fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
          return FC_MEMORY_ERROR;
        }
        vars[i] = temp_vars;
      }
    }

    // Assign time coord value
    if (wordlen == 4) 
      ((float*)time_coords)[numStep] = _fc_getDoubleFromBuffer(wordBuf, wordlen,0);
    else
      ((double*)time_coords)[numStep] = _fc_getDoubleFromBuffer(wordBuf,wordlen,0);

    // init seqVarID
    seqVarID = 0;

    // True Global variables
    for (i = 0 ; i < numGlbVar; i++) {
      rc = fc_createGlobalVariable(*dataset, globalVarNames[i],
                                   &vars[seqVarID][numStep]);
      if (rc != FC_SUCCESS) {
        fc_printfErrorMessage("failed to create global var '%s'", 
                              globalVarNames[i]);
          return rc;
      }
      varSlot = _fc_getVarSlot(vars[seqVarID][numStep]);
      varSlot->numDataPoint = 1;
      varSlot->numComponent = globalVarNumComps[i];
      varSlot->assoc = FC_AT_WHOLE_DATASET;
      if (globalVarNumComps[i] == 1)
        varSlot->mathtype = FC_MT_SCALAR;
      else
        varSlot->mathtype = FC_MT_VECTOR;
      varSlot->datatype = datatype;
      varSlot->header.committed = 1;
      varSlot->fileInfo.dynaInfo.dataFileNum = fileNum;
      varSlot->fileInfo.dynaInfo.data_offset = ftell(file);
      varSlot->fileInfo.dynaInfo.doDispl = 0;
      varSlot->fileInfo.dynaInfo.isStriped = 0;
      seqVarID++;
      fseek(file, globalVarNumComps[i]*wordlen, SEEK_CUR);
    }

    // Per material "global" variables - each mesh gets a set
    for (i = 0; i < numPerMatGlbVar; i++) {
      for (j = 0; j < numTotalMesh; j++) {
        if (!FC_HANDLE_EQUIV(meshes[j], FC_NULL_MESH)) {
          rc = fc_createVariable(meshes[j], perMaterialVarNames[i],
                                 &vars[seqVarID][numStep]);
          if (rc != FC_SUCCESS) {
            fc_printfErrorMessage("failed to create per material global var "
                                  "'%s' on mesh #%d", perMaterialVarNames[i],
                                  j);
            return rc;
          }
          varSlot = _fc_getVarSlot(vars[seqVarID][numStep]);
          varSlot->numDataPoint = 1;
          varSlot->numComponent = perMatVarComps[i];
          varSlot->assoc = FC_AT_WHOLE_MESH;
          if (perMatVarComps[i] == 1)
            varSlot->mathtype = FC_MT_SCALAR;
          else
            varSlot->mathtype = FC_MT_VECTOR;
          varSlot->datatype = datatype;
          varSlot->header.committed = 1;
          varSlot->fileInfo.dynaInfo.dataFileNum = fileNum;
          varSlot->fileInfo.dynaInfo.data_offset = ftell(file);
          varSlot->fileInfo.dynaInfo.doDispl = 0;
          varSlot->fileInfo.dynaInfo.isStriped = 0;
          seqVarID++;
        }
        fseek(file, perMatVarComps[i]*wordlen, SEEK_CUR);
      }
    }

    // Per rigid wall global variables (treat as global variable)
    for (i = 0 ; i < numRW; i++) {     
      // index from 1 so have same names as lsprepost has
      sprintf(tempName, "force_on_rw_%d", i+1); 
      //sprintf(tempName, "global_%d", i);
      
      rc = fc_createGlobalVariable(*dataset, tempName, &vars[seqVarID][numStep]);
      if (rc != FC_SUCCESS) {
        fc_printfErrorMessage("failed to create rigid wall var #%d", i);
        return rc;
      }
      varSlot = _fc_getVarSlot(vars[seqVarID][numStep]);
      varSlot->numDataPoint = 1;
      varSlot->numComponent = 1;
      varSlot->assoc = FC_AT_WHOLE_DATASET;
      varSlot->mathtype = FC_MT_SCALAR;
      varSlot->datatype = datatype;
      varSlot->header.committed = 1;
      varSlot->fileInfo.dynaInfo.dataFileNum = fileNum;
      varSlot->fileInfo.dynaInfo.data_offset = ftell(file);
      varSlot->fileInfo.dynaInfo.doDispl = 0;
      varSlot->fileInfo.dynaInfo.isStriped = 0;
      seqVarID++;
      fseek(file, wordlen, SEEK_CUR);
    }

    // Displacement, temperature, velocity & acceleration vectors per node
    // (NOTE: This ordering differs from that implied in the documentation
    // which implies that temperature is first.)
    for (i = 0; i < doDisplVar + doTempVar + doVelVar + doAccelVar; i++) {
      char *seqVarName;
      FC_MathType mathType;
      int temp_numDim;
      if (doDisplVar && i == 0) {
        seqVarName = displVarName;
	mathType = FC_MT_VECTOR;
	temp_numDim = numDim;
      }
      else if (doTempVar && i == doDisplVar) {
	seqVarName = tempVarName;
	mathType = FC_MT_SCALAR;
	temp_numDim = 1;
      }
      else if (doVelVar && i == doDisplVar + doTempVar) {
        seqVarName = velVarName;
	mathType = FC_MT_VECTOR;
	temp_numDim = numDim;
      }
      else if (doAccelVar && i == doDisplVar + doTempVar + doVelVar) {
        seqVarName = accelVarName;
	mathType = FC_MT_VECTOR;
	temp_numDim = numDim;
      }
      for (j = 0; j < numTotalMesh; j++) {
	if (FC_HANDLE_EQUIV(meshes[j], FC_NULL_MESH))
	  continue;
        rc = fc_createVariable(meshes[j], seqVarName,
                               &vars[seqVarID][numStep]);
        if (rc != FC_SUCCESS) {
          fc_printfErrorMessage("failed to vector var '%s' on mesh #%d", 
                                seqVarName, j);
          return rc;
        }
        varSlot = _fc_getVarSlot(vars[seqVarID][numStep]);
        meshSlot = _fc_getMeshSlot(meshes[j]);
        varSlot->numDataPoint = meshSlot->numVertex;
	varSlot->numComponent = temp_numDim;
        varSlot->assoc = FC_AT_VERTEX;
	varSlot->mathtype = mathType;
        varSlot->datatype = datatype;
        varSlot->header.committed = 1;
        varSlot->fileInfo.dynaInfo.dataFileNum = fileNum;
        varSlot->fileInfo.dynaInfo.data_offset = ftell(file);
        if (i == 0 && doDisplVar)
          varSlot->fileInfo.dynaInfo.doDispl = 1;
        else
          varSlot->fileInfo.dynaInfo.doDispl = 0;
        varSlot->fileInfo.dynaInfo.isStriped = 0;
        seqVarID++;
      }
      fseek(file, temp_numDim*numTotalVert*wordlen, SEEK_CUR);
    }

    // Element vars - These are written as all values for a data point,
    // then all values for the next data point, so we have to unstripe them
    // per element type

    // Element vars for elem8
    for (i = 0; i < numTotalMesh; i++) {
      if (elemTypes[i] != FC_ET_HEX)
	continue;
      for (j = 0; j < numElemVar8; j++) {
        char *seqVarName = tempName;
        if (j < 7)
          seqVarName = elem8VarNames[j];
        else // start with 1 to agree with LSPrepost naming
          sprintf(seqVarName, "elem_var_%d", j - 7 + 1);
        rc = fc_createVariable(meshes[i], seqVarName,
                               &vars[seqVarID][numStep]); 
        if (rc != FC_SUCCESS) {
          fc_printfErrorMessage("failed to create scalar var '%s' on "
                                "mesh #%d", seqVarName, i);
          return rc;
        }
        varSlot = _fc_getVarSlot(vars[seqVarID][numStep]);
        meshSlot = _fc_getMeshSlot(meshes[i]);
        varSlot->numDataPoint = meshSlot->numElement;
        varSlot->numComponent = 1;
        varSlot->assoc = FC_AT_ELEMENT;
        varSlot->mathtype = FC_MT_SCALAR;
        varSlot->datatype = datatype;
        varSlot->header.committed = 1;
        varSlot->fileInfo.dynaInfo.dataFileNum = fileNum;
        varSlot->fileInfo.dynaInfo.data_offset = ftell(file);
        varSlot->fileInfo.dynaInfo.isStriped = numElemVar8;
        varSlot->fileInfo.dynaInfo.stripeID = j;
        seqVarID++;
      }
    }
    fseek(file, numTotalElem8*numElemVar8*wordlen, SEEK_CUR);

    // Elem vars for elem2
    for (i = 0; i < numTotalMesh; i++) {
      if (elemTypes[i] != FC_ET_LINE)
	continue;
      for (j = 0; j < numElemVar2; j++) {
	char *seqVarName = tempName;
	if (j < 6)
	  seqVarName = elem2VarNames[j];
	else {
	  int IP = (j-6)/5 + 1; // integration point #
	  int IPVarID = (j-6) % 5;
	  sprintf(seqVarName, "%s-IP%d", elem2VarNames[6+IPVarID], IP);
	}
	rc = fc_createVariable(meshes[i], seqVarName,
                               &vars[seqVarID][numStep]); 
        if (rc != FC_SUCCESS) {
          fc_printfErrorMessage("failed to create scalar var '%s' on "
                                "mesh #%d", seqVarName, i);
          return rc;
        }
        varSlot = _fc_getVarSlot(vars[seqVarID][numStep]);
        meshSlot = _fc_getMeshSlot(meshes[i]);
        varSlot->numDataPoint = meshSlot->numElement;
        varSlot->numComponent = 1;
        varSlot->assoc = FC_AT_ELEMENT;
        varSlot->mathtype = FC_MT_SCALAR;
        varSlot->datatype = datatype;
        varSlot->header.committed = 1;
        varSlot->fileInfo.dynaInfo.dataFileNum = fileNum;
        varSlot->fileInfo.dynaInfo.data_offset = ftell(file);
        varSlot->fileInfo.dynaInfo.isStriped = numElemVar2;
        varSlot->fileInfo.dynaInfo.stripeID = j;
        seqVarID++;
      }
    }
    fseek(file, wordlen*numTotalElem2*numElemVar2, SEEK_CUR);

    // Elem vars for elem4
    for (i = 0; i < numTotalMesh; i++) {
      if (elemTypes[i] != FC_ET_QUAD)
	continue;
      for (j = 0; j < numElemVar4; j++) {
	char *seqVarName = tempName;
	if (j < 7)
	  sprintf(seqVarName, "mid-%s", elem8VarNames[j]);
	else if (j < 14)
	  sprintf(seqVarName, "in-%s", elem8VarNames[j-7]);
	else if (j < 21)
	  sprintf(seqVarName, "out-%s", elem8VarNames[j-14]);
	else if (j == numElemVar4-1)
	  seqVarName = elem4VarNames[23]; // IE is always last
	else
	  seqVarName = elem4VarNames[j-21];
	rc = fc_createVariable(meshes[i], seqVarName,
                               &vars[seqVarID][numStep]); 
        if (rc != FC_SUCCESS) {
          fc_printfErrorMessage("failed to create scalar var '%s' on "
                                "mesh #%d", seqVarName, i);
          return rc;
        }
        varSlot = _fc_getVarSlot(vars[seqVarID][numStep]);
        meshSlot = _fc_getMeshSlot(meshes[i]);
        varSlot->numDataPoint = meshSlot->numElement;
        varSlot->numComponent = 1;
        varSlot->assoc = FC_AT_ELEMENT;
        varSlot->mathtype = FC_MT_SCALAR;
        varSlot->datatype = datatype;
        varSlot->header.committed = 1;
        varSlot->fileInfo.dynaInfo.dataFileNum = fileNum;
        varSlot->fileInfo.dynaInfo.data_offset = ftell(file);
        varSlot->fileInfo.dynaInfo.isStriped = numElemVar4;
        varSlot->fileInfo.dynaInfo.stripeID = j;
        seqVarID++;
      }
    }
    fseek(file, wordlen*numTotalElem4*numElemVar4, SEEK_CUR);

    // element death var for all elems
    if (doElemDeathVar) {
      for (i = 0; i < numTotalMesh; i++) {
	if (FC_HANDLE_EQUIV(meshes[i], FC_NULL_MESH))
	  continue;
        rc = fc_createVariable(meshes[i], elemDeathVarName, 
                               &vars[seqVarID][numStep]);
        if (rc != FC_SUCCESS) {
          fc_printfErrorMessage("failed to create elem_death var on "
                                "mesh #%d", i);

          return rc;
        }
        varSlot = _fc_getVarSlot(vars[seqVarID][numStep]);
        meshSlot = _fc_getMeshSlot(meshes[i]);
        varSlot->numDataPoint = meshSlot->numElement;
        varSlot->numComponent = 1;
        varSlot->assoc = FC_AT_ELEMENT;
        varSlot->mathtype = FC_MT_SCALAR;
        varSlot->datatype = datatype;
        varSlot->header.committed = 1;
        varSlot->fileInfo.dynaInfo.dataFileNum = fileNum;
	if (elemTypes[i] == FC_ET_HEX)
	  varSlot->fileInfo.dynaInfo.data_offset = ftell(file);
	else if (elemTypes[i] == FC_ET_QUAD)
	  varSlot->fileInfo.dynaInfo.data_offset =
	                        ftell(file) + numTotalElem8*wordlen;
	else if (elemTypes[i] == FC_ET_LINE)
	  varSlot->fileInfo.dynaInfo.data_offset = 
	    ftell(file) + (numTotalElem8+numTotalElem4)*wordlen;
        varSlot->fileInfo.dynaInfo.isStriped = 0;
        seqVarID++;
      }
      fseek(file, (numTotalElem8+numTotalElem4+numTotalElem2)*wordlen, 
	    SEEK_CUR);
    }

    // increment step ID
    numStep++;

    // debug
    //printf("numSeqVar = %d, seqVarID = %d\n", numSeqVar, seqVarID);
  }

  // debug
  //for (i = 0; i < numStep; i++) 
  //printf("%d: %f\n", i, time_coords[i]);

  //--- Create sequence

  realloc(time_coords, numStep*wordlen);
  rc = fc_createSequence(*dataset, "time", &sequence);
  if (rc != FC_SUCCESS) {
    fc_printfErrorMessage("Failed to create sequence");
    return rc;
  }
  rc = fc_setSequenceCoordsPtr(sequence, numStep, datatype, time_coords); 
  if (rc != FC_SUCCESS)
    return rc;

  //--- Create sequence variables - convert the created vars to seq variables
  for (i = 0; i < numSeqVar; i++) {
    rc = fc_convertVariablesToSeqVariable(numStep, vars[i], sequence, 
                                          NULL, &seqVar);
    free(seqVar);
    if (rc != FC_SUCCESS) {
      fc_printfErrorMessage("Failed to convert variable #%d to sequence var",
                            i);
      return rc;
    }
  }

  //--- All done

  free(elemTypes);
  free(meshes);
  for (i = 0; i < numSeqVar; i++) 
    free(vars[i]);
  free(vars);

  return FC_SUCCESS;
}

/**
 * \ingroup  LSDynaFileIO
 * \brief  Read data from LS-DYNA file into a mesh
 *
 * \description
 *  
 *    Assuming that the given mesh had already been setup, then read the 
 *    mesh vertex coords.
 *
 * \modifications
 *   - 7/18/2005 WSD.
 */
FC_ReturnCode _fc_readLSDynaMeshCoords(
  FC_Mesh mesh     /**< input - mesh to be read into */
) {
  FC_ReturnCode rc;
  int i, j;
  _FC_MeshSlot* meshSlot;
  _FC_DsSlot* dsSlot;
  double* global_coords;
  int global_count;
  int local_count;
  FILE* file;
  int* vertIDs;
  char filename[1024]; // root names are limited to 75 char 
  int fileNum;

  // check input
  meshSlot = _fc_getMeshSlot(mesh);
  if (!meshSlot || !meshSlot->header.name) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Reading coords for  mesh '%s'", meshSlot->header.name);

  // setup
  dsSlot = _fc_getDsSlot(meshSlot->ds);
  vertIDs = dsSlot->fileInfo.dynaInfo.vertIDs[meshSlot->fileInfo.dynaInfo.meshID];
  global_count = meshSlot->dim * dsSlot->fileInfo.dynaInfo.numTotalVert;
  local_count = meshSlot->dim * meshSlot->numVertex;

  // open the file
  fileNum = dsSlot->fileInfo.dynaInfo.coordsFileNum;
  if (fileNum  == 0)
    strcpy(filename, dsSlot->baseFile);
  else if (fileNum < 100)
    sprintf(filename, "%s%02d", dsSlot->baseFile, fileNum);
  else
    sprintf(filename, "%s%d", dsSlot->baseFile, fileNum);
  file = fopen(filename, "r");
  if (file == NULL) {
    fc_printfErrorMessage("could't open file '%s'", filename);
    return FC_FILE_IO_ERROR;
  }

   // Set up
  // Read global coords 
  //FIX? I don't think there is any good way to avoid loading all coords
  // FIX fiddle with reading different sizes
  global_coords = (double*)malloc(global_count*sizeof(double));  
  if (!global_coords) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }  
  fseek(file, dsSlot->fileInfo.dynaInfo.coords_offset, SEEK_SET);
  if (dsSlot->fileInfo.dynaInfo.wordlen == 4) {
    float* temp_coords = (float*)malloc(global_count*sizeof(float));
    if (!temp_coords) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      return FC_MEMORY_ERROR;
    }
    rc = _fc_fread_endianSwap(temp_coords, dsSlot->fileInfo.dynaInfo.wordlen, 
                             global_count, &file, dsSlot->baseFile, &fileNum, 
                             dsSlot->fileInfo.dynaInfo.doEndianSwap, NULL);
    if (rc != FC_SUCCESS) {
      fc_printfErrorMessage("Failed to read mesh coords");
      return rc;
    }
    for (i = 0; i < global_count; i++)
      global_coords[i] = _fc_getDoubleFromBuffer(temp_coords, 
                                             dsSlot->fileInfo.dynaInfo.wordlen,
                                             i);
    free(temp_coords);
  }
  else {
    rc = _fc_fread_endianSwap(global_coords, dsSlot->fileInfo.dynaInfo.wordlen,
                             global_count, &file, dsSlot->baseFile, &fileNum, 
                             dsSlot->fileInfo.dynaInfo.doEndianSwap, NULL);
    if (rc != FC_SUCCESS) {
      fc_printfErrorMessage("Failed to read mesh coords");
      return rc;
    }
  }

  // close file
  fclose(file);

  // make local coords array
  meshSlot->coords = (double*)malloc(local_count*sizeof(double));
  if (!meshSlot->coords) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  // FIX! use a memcpy here, get rid of local_count
  for (i = 0; i < meshSlot->numVertex; i++)
    for (j = 0; j < meshSlot->dim; j++)
      meshSlot->coords[i*meshSlot->dim + j] = 
        global_coords[vertIDs[i]*meshSlot->dim + j]; 

  // cleanup
  free(global_coords);
  
  return FC_SUCCESS;
}

/**
 * \ingroup  LSDynaFileIO
 * \brief  Read data from LS-DYNA file into a mesh
 *
 * \description
 *  
 *    Assuming that the given mesh had already been setup, then read the 
 *    mesh element connectivities.
 *
 * \modifications
 *   - 7/18/2005 WSD.
 */
FC_ReturnCode _fc_readLSDynaMeshElemConns(
  FC_Mesh mesh     /**< input - mesh to be read into */
) {
  FC_ReturnCode rc;
  int i, j, k;
  _FC_MeshSlot* meshSlot;
  _FC_DsSlot* dsSlot;
  FILE* file;
  int* conns;
  void* temp_conns;
  int numVertPerElem, numPlus;
  int* vertLU, numVert;
  char filename[1024]; // root names are limited to 75 char 
  int fileNum;
  int meshID, numElem;
  long offset;

  // check input
  meshSlot = _fc_getMeshSlot(mesh);
  if (!meshSlot || !meshSlot->header.name) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Reading conns for  mesh '%s'", meshSlot->header.name);

  // setup
  meshID = meshSlot->fileInfo.dynaInfo.meshID;
  dsSlot = _fc_getDsSlot(meshSlot->ds);
  numVertPerElem = fc_getElementTypeNumVertex(meshSlot->elemType);
  numPlus = numVertPerElem + 1;
  if (meshSlot->elemType == FC_ET_LINE)
    numPlus = 6;
  if (meshSlot->elemType == FC_ET_HEX)
    offset = dsSlot->fileInfo.dynaInfo.conns8_offset;
  else if (meshSlot->elemType == FC_ET_QUAD)
    offset = dsSlot->fileInfo.dynaInfo.conns4_offset;
  else if (meshSlot->elemType == FC_ET_LINE)
    offset = dsSlot->fileInfo.dynaInfo.conns2_offset;
 
  // open the file
  fileNum = dsSlot->fileInfo.dynaInfo.coordsFileNum;
  if (fileNum  == 0)
    strcpy(filename, dsSlot->baseFile);
  else if (fileNum < 100)
    sprintf(filename, "%s%02d", dsSlot->baseFile, fileNum);
  else
    sprintf(filename, "%s%d", dsSlot->baseFile, fileNum);
  file = fopen(filename, "r");
  if (file == NULL) {
    fc_printfErrorMessage("could't open file '%s'", filename);
    return FC_FILE_IO_ERROR;
  }

  // setup - make room for the conns
  conns = malloc(meshSlot->numElement*numVertPerElem*sizeof(int));
  if (!conns) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      return FC_MEMORY_ERROR;
  }

  // Read conns in blocks (but will be in global verts + has partID)
  // copy into conns into new array w/o partID
  numElem = 0;
  for (i = 0; i < dsSlot->fileInfo.dynaInfo.numElemChunks[meshID]; i++) {
    int numElemPerChunk = dsSlot->fileInfo.dynaInfo.numElemPerChunk[meshID][i];
    temp_conns = malloc(numElemPerChunk*numPlus*
			dsSlot->fileInfo.dynaInfo.wordlen);
    if (!temp_conns) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      return FC_MEMORY_ERROR;
    }  
    fseek(file, offset + numPlus*dsSlot->fileInfo.dynaInfo.wordlen*
	  dsSlot->fileInfo.dynaInfo.startElemPerChunk[meshID][i], SEEK_SET);
    rc = _fc_fread_endianSwap(temp_conns, dsSlot->fileInfo.dynaInfo.wordlen, 
			     numElemPerChunk*numPlus, &file, 
			     dsSlot->baseFile, &fileNum, 
			     dsSlot->fileInfo.dynaInfo.doEndianSwap, NULL);
    if (rc != FC_SUCCESS) {
      fc_printfErrorMessage("Failed to read mesh conns");
      free(temp_conns);
      return rc;
    }

    if (dsSlot->fileInfo.dynaInfo.wordlen == 4) {
      int* cast_conns = (int*)temp_conns;
      for (j = 0; j < numElemPerChunk; j++) { 
	for (k = 0; k < numVertPerElem; k++) 
	  conns[numElem*numVertPerElem+k] = cast_conns[j*numPlus+k];
	numElem++;
      }
    }
    else {
      long long* cast_conns = (long long*)temp_conns;
      for (j = 0; j < numElemPerChunk; j++) { 
	for (k = 0; k < numVertPerElem; k++)
	  conns[numElem*numVertPerElem+k] = cast_conns[j*numPlus+k];
	numElem++;
      }
    }
    free(temp_conns);
  }
  if (numElem != meshSlot->numElement) {
    printf("****** oh no, mismatch of numElement ******!\n");
    fflush(NULL);
  }

  // Close file
  fclose(file);

  // Convert to local verts
  // Make vertLU
  vertLU = (int*)malloc(dsSlot->fileInfo.dynaInfo.numTotalVert*sizeof(int));
  if (!vertLU) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  for (i = 0; i < dsSlot->fileInfo.dynaInfo.numTotalVert; i++)
    vertLU[i] = -1;
  for (i = 0; i < meshSlot->numElement; i++) 
    for (j = 0; j < numVertPerElem; j++)
      vertLU[conns[i*numVertPerElem+j]-1] = 0;
  numVert = 0;
  for (i = 0; i <dsSlot->fileInfo.dynaInfo.numTotalVert; i++) {
    if (vertLU[i] > -1) {
      vertLU[i] = numVert;
      numVert++;
    }
  }
  if (numVert != meshSlot->numVertex) {
    printf("****** oh no, mismatch of numVertex ******!\n");
    fflush(NULL);
  }
  for (i = 0; i < numVertPerElem*meshSlot->numElement; i++)
    conns[i] = vertLU[conns[i]-1];
  free(vertLU);

  // done
  meshSlot->elemToVertConns = conns;
  return FC_SUCCESS;
}

/**
 * \ingroup  LSDynaFileIO
 * \brief  Read data from LS-DYNA file into a variable
 *
 * \description
 *  
 *    Assuming that the given mesh had already been setup, then read the 
 *    variable data.
 *
 * \modifications
 *   - 7/19/2005 WSD.
 */
FC_ReturnCode _fc_readLSDynaVariableData(
  FC_Variable variable   /**< input - variable to be read into */
) {
  FC_ReturnCode rc;
  int i, j;
  _FC_VarSlot* varSlot;
  _FC_DsSlot* dsSlot;
  _FC_MeshSlot* meshSlot;
  FILE* file;
  void* temp_data;
  size_t read_number;
  double* coords;
  char filename[1024]; // root names are limited to 75 char 
  int fileNum;
  int meshID, numElem;

  // check input
  varSlot = _fc_getVarSlot(variable);
  if (!varSlot || !varSlot->header.name) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Reading data for variable '%s'", varSlot->header.name);

  // setup
  dsSlot = _fc_getDsSlot(varSlot->ds);
  meshSlot = _fc_getMeshSlot(varSlot->mesh);
  if (meshSlot)
    meshID = meshSlot->fileInfo.dynaInfo.meshID;

  // open the file
  fileNum = varSlot->fileInfo.dynaInfo.dataFileNum;
  if (fileNum  == 0)
    strcpy(filename, dsSlot->baseFile);
  else if (fileNum < 100)
    sprintf(filename, "%s%02d", dsSlot->baseFile, fileNum);
  else
    sprintf(filename, "%s%d", dsSlot->baseFile, fileNum);
  file = fopen(filename, "r");
  if (file == NULL) {
    fc_printfErrorMessage("could't open file '%s'", filename);
    return FC_FILE_IO_ERROR;
  }
  
  // 4 different cases - vert, global, elem & striped elem

  // Vertex variables need to be unpacked - read all and then copy
  // just the ones for this mesh
  if (varSlot->assoc == FC_AT_VERTEX) {
    int* vertIDs = dsSlot->fileInfo.dynaInfo.vertIDs[meshID];
    read_number = dsSlot->fileInfo.dynaInfo.numTotalVert*varSlot->numComponent;
    temp_data = malloc(read_number*dsSlot->fileInfo.dynaInfo.wordlen);
    varSlot->data = malloc(varSlot->numDataPoint*varSlot->numComponent*
                           fc_sizeofDataType(varSlot->datatype));
    if (!temp_data || !varSlot->data) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      return FC_MEMORY_ERROR;
    }
    fseek(file, varSlot->fileInfo.dynaInfo.data_offset, SEEK_SET);
    rc = _fc_fread_endianSwap(temp_data, dsSlot->fileInfo.dynaInfo.wordlen, 
			     read_number, &file, dsSlot->baseFile, &fileNum, 
			     dsSlot->fileInfo.dynaInfo.doEndianSwap, NULL);
    if (rc != FC_SUCCESS) {
      fc_printfErrorMessage("Failed to read variable data");
      free(temp_data); 
      free(varSlot->data); varSlot->data = NULL;
      return rc;
    }
    if (varSlot->datatype == FC_DT_INT && 
	dsSlot->fileInfo.dynaInfo.wordlen == 4) {
      int* data = (int*)varSlot->data;
      int* tdata = (int*)temp_data;
      for (i = 0; i < varSlot->numDataPoint; i++) 
        memcpy(&data[i*varSlot->numComponent], 
               &tdata[vertIDs[i]*varSlot->numComponent],
               varSlot->numComponent*dsSlot->fileInfo.dynaInfo.wordlen);
    }
    else if (varSlot->datatype == FC_DT_INT && 
	     dsSlot->fileInfo.dynaInfo.wordlen == 8) {
      int* data = (int*)varSlot->data;
      long long* tdata = (long long*)temp_data;
      for (i = 0; i < varSlot->numDataPoint; i++) 
	for (j = 0; j < varSlot->numComponent; j++) 
	  data[i*varSlot->numComponent+j] = tdata[vertIDs[i]*varSlot->numComponent+j];
    
    }
    else if (varSlot->datatype == FC_DT_FLOAT) {
      float* data = (float*)varSlot->data;
      float* tdata = (float*)temp_data;
      for (i = 0; i < varSlot->numDataPoint; i++) 
        memcpy(&data[i*varSlot->numComponent], 
               &tdata[vertIDs[i]*varSlot->numComponent],
               varSlot->numComponent*dsSlot->fileInfo.dynaInfo.wordlen);
      if (varSlot->fileInfo.dynaInfo.doDispl) {
        rc = fc_getMeshCoordsPtr(varSlot->mesh, &coords);
        for (i = 0; i < varSlot->numDataPoint*varSlot->numComponent; i++)
          data[i] -= coords[i];
      }
    }
    else if (varSlot->datatype == FC_DT_DOUBLE) {
      double* data = (double*)varSlot->data;
      double* tdata = (double*)temp_data;
      for (i = 0; i < varSlot->numDataPoint; i++) 
        memcpy(&data[i*varSlot->numComponent], 
               &tdata[vertIDs[i]*varSlot->numComponent],
               varSlot->numComponent*dsSlot->fileInfo.dynaInfo.wordlen);
      if (varSlot->fileInfo.dynaInfo.doDispl) {
        rc = fc_getMeshCoordsPtr(varSlot->mesh, &coords);
        for (i = 0; i < varSlot->numDataPoint*varSlot->numComponent; i++)
          data[i] -= coords[i];
      }
    }
    free(temp_data);
  }
   
  // Global & whole vars can be read as is
  // assumes we will never see int global vars
  else if (varSlot->assoc == FC_AT_WHOLE_MESH || 
	   varSlot->assoc == FC_AT_WHOLE_DATASET) {
    read_number = varSlot->numDataPoint*varSlot->numComponent;
    temp_data = malloc(read_number*dsSlot->fileInfo.dynaInfo.wordlen);
    if (!temp_data) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      return FC_MEMORY_ERROR;
    }
    fseek(file, varSlot->fileInfo.dynaInfo.data_offset, SEEK_SET);
    rc = _fc_fread_endianSwap(temp_data, dsSlot->fileInfo.dynaInfo.wordlen, 
			     read_number, &file, dsSlot->baseFile, &fileNum, 
			     dsSlot->fileInfo.dynaInfo.doEndianSwap, NULL);
    if (rc != FC_SUCCESS) {
      fc_printfErrorMessage("Failed to read variable data");
      free(temp_data); 
      return rc;
    }
    varSlot->data = temp_data;
  }

  // unstriped element data can be read as is, but in chunks
  // assume to be single comopnent
  else if (varSlot->assoc == FC_AT_ELEMENT &&
	   varSlot->fileInfo.dynaInfo.isStriped == 0) {
    varSlot->data = malloc(varSlot->numDataPoint*
			   fc_sizeofDataType(varSlot->datatype));
    if (!varSlot->data) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      return FC_MEMORY_ERROR;
    }
    numElem = 0;
    for (i = 0; i < dsSlot->fileInfo.dynaInfo.numElemChunks[meshID]; i++) {
      int numElemPerChunk = dsSlot->fileInfo.dynaInfo.numElemPerChunk[meshID][i];
      void* current_data;
      if (varSlot->datatype == FC_DT_INT)
	current_data = (void*) &(((int*)varSlot->data)[numElem]);
      if (varSlot->datatype == FC_DT_FLOAT)
	current_data = (void*) &(((float*)varSlot->data)[numElem]);
      else
	current_data = (void*) &(((double*)varSlot->data)[numElem]);
      fseek(file, varSlot->fileInfo.dynaInfo.data_offset + 
	    dsSlot->fileInfo.dynaInfo.wordlen*
	    dsSlot->fileInfo.dynaInfo.startElemPerChunk[meshID][i], SEEK_SET);
      if (varSlot->datatype == FC_DT_INT && 
	  dsSlot->fileInfo.dynaInfo.wordlen == 8) {
	long long* tdata = (long long*)malloc(numElemPerChunk*
					      dsSlot->fileInfo.dynaInfo.wordlen);
	if (!tdata) {
	  fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
	  return FC_MEMORY_ERROR;
	}
	rc = _fc_fread_endianSwap(tdata, dsSlot->fileInfo.dynaInfo.wordlen,
				  numElemPerChunk, &file, dsSlot->baseFile, &fileNum,
				  dsSlot->fileInfo.dynaInfo.doEndianSwap, NULL);
	if (rc != FC_SUCCESS) {
	  fc_printfErrorMessage("Failed to read variable data");
	  free(varSlot->data);
	  varSlot->data = NULL;
	  return rc;
	}
	for (j = 0; j < numElemPerChunk; j++) {
	  ((int*)varSlot->data)[numElem+j] = tdata[j];
	}
	free(tdata);
      }
      else {
	rc = _fc_fread_endianSwap(current_data, dsSlot->fileInfo.dynaInfo.wordlen,
				  numElemPerChunk, &file, dsSlot->baseFile, &fileNum, 
				  dsSlot->fileInfo.dynaInfo.doEndianSwap, NULL);
	if (rc != FC_SUCCESS) {
	  fc_printfErrorMessage("Failed to read variable data");
	  free(varSlot->data);
	  varSlot->data = NULL;
	  return rc;
	}
      }
      numElem += numElemPerChunk;
    }
  }

  // striped element data must be unstriped, and unchunked
  // assumes single component
  // assumes we will never see striped int data
  else if (varSlot->assoc == FC_AT_ELEMENT && 
	   varSlot->fileInfo.dynaInfo.isStriped > 0) {
    varSlot->data = malloc(varSlot->numDataPoint*dsSlot->fileInfo.dynaInfo.wordlen);
    if (!varSlot->data) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      return FC_MEMORY_ERROR;
    }
    numElem = 0;
    for (i = 0; i < dsSlot->fileInfo.dynaInfo.numElemChunks[meshID]; i++) {
      int numElemPerChunk = dsSlot->fileInfo.dynaInfo.numElemPerChunk[meshID][i];
      read_number = numElemPerChunk*varSlot->fileInfo.dynaInfo.isStriped;
      temp_data = malloc(read_number*dsSlot->fileInfo.dynaInfo.wordlen);
      if (!temp_data) {
	fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
	return FC_MEMORY_ERROR;
      }
      fseek(file, varSlot->fileInfo.dynaInfo.data_offset + 
	    varSlot->fileInfo.dynaInfo.isStriped*dsSlot->fileInfo.dynaInfo.wordlen*
	    dsSlot->fileInfo.dynaInfo.startElemPerChunk[meshID][i], SEEK_SET);
      rc = _fc_fread_endianSwap(temp_data, dsSlot->fileInfo.dynaInfo.wordlen, 
			       read_number, &file, dsSlot->baseFile, &fileNum, 
			       dsSlot->fileInfo.dynaInfo.doEndianSwap, NULL);
      if (rc != FC_SUCCESS) {
	fc_printfErrorMessage("Failed to read variable data");
	free(temp_data); 
	free(varSlot->data); varSlot->data = NULL;
	return rc;
      }
      if (varSlot->datatype == FC_DT_FLOAT) {
	float* data = &(((float*)varSlot->data)[numElem]);
	float* tdata = (float*)temp_data;
	for (j = 0; j < numElemPerChunk; j++) 
	  data[j] = tdata[j*varSlot->fileInfo.dynaInfo.isStriped + 
			  varSlot->fileInfo.dynaInfo.stripeID];
      }
      else if (varSlot->datatype == FC_DT_DOUBLE) {
	double* data = &(((double*)varSlot->data)[numElem]);
	double* tdata = (double*)temp_data;
	for (j = 0; j < numElemPerChunk; j++) 
	  data[j] = tdata[j*varSlot->fileInfo.dynaInfo.isStriped + 
			  varSlot->fileInfo.dynaInfo.stripeID];
      }
      else {
	fc_printfErrorMessage("Unknown data type");
	return FC_ERROR;
      }
      numElem += numElemPerChunk;
      free(temp_data);
    }
  }

  // fifth case - what if not one of the other 4 cases
  else {
    fc_printfErrorMessage("**Algorithm error, should never reach this point");
    return FC_ERROR;
  }

  // close the file
  fclose(file);

  return FC_SUCCESS;
}

/**
 * \ingroup  LSDynaFileIO
 * \brief  Close a dataset.
 *
 * \description
 *  
 *    Do any necessary file related cleanup for the Dataset.
 *
 * \modifications  
 *   - 7/19/2005 WSD, created.
 */
FC_ReturnCode _fc_closeLSDynaDataset(
  FC_Dataset dataset /**< input - the dataset */
) {
  int i;
  _FC_DsSlot* dsSlot;

  // check input
  dsSlot = _fc_getDsSlot(dataset);
  if (dsSlot == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // free dynamically allocated memory in db
  for (i = 0; i < dsSlot->fileInfo.dynaInfo.numMesh; i++) {
    free(dsSlot->fileInfo.dynaInfo.vertIDs[i]);
    free(dsSlot->fileInfo.dynaInfo.numElemPerChunk[i]);
    free(dsSlot->fileInfo.dynaInfo.startElemPerChunk[i]);
  }
  free(dsSlot->fileInfo.dynaInfo.vertIDs);
  free(dsSlot->fileInfo.dynaInfo.numElemChunks);
  free(dsSlot->fileInfo.dynaInfo.numElemPerChunk);
  free(dsSlot->fileInfo.dynaInfo.startElemPerChunk);

  return FC_SUCCESS;
}

/**
 * \ingroup  LSDynaFileIO
 * \brief  Ask if any vertices are shared.
 *
 * \description 
 *
 *    Because LS-DYNA has a single global list of vertices and FCLib has local
 *    lists of vertices per mesh, it is possible for a single global vertex to
 *    show up in multiple meshes.
 *
 *    Returns true (>1), false (0), or -1 if the dataset is not an
 *    LS-DYNA Dataset.
 *
 * \modifications  
 *   - 2/3/2006 WSD, created.
 */
int _fc_hasSharedLSDynaVertices(
  FC_Dataset dataset /**< input - the dataset */
) {
  FC_ReturnCode rc;
  int i,j;
  _FC_DsSlot* dsSlot;
  int* counts, numShare;
  int numMesh, numVert;
  FC_Mesh* meshes;
  _FC_MeshSlot *meshSlot;

  // check input
  dsSlot = _fc_getDsSlot(dataset);
  if (dsSlot == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // setup
  rc = fc_getMeshes(dataset, &numMesh, &meshes);
  if (rc != FC_SUCCESS)
    return rc;

  // Early return - no meshes, no shared global verts
  if (numMesh < 1)
    return 0;

  // collect counts
  counts = (int*)calloc(dsSlot->fileInfo.dynaInfo.numTotalVert, sizeof(int));
  if (!counts) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  for (i = 0; i < numMesh; i++) {
    int meshID;
    meshSlot = _fc_getMeshSlot(meshes[i]);
    meshID = meshSlot->fileInfo.dynaInfo.meshID;
    if (meshID < 0 || meshID >= dsSlot->fileInfo.dynaInfo.numMesh)
      continue; // this isn't an lsdyna mesh
    fc_getMeshNumVertex(meshes[i], &numVert);
    if (rc != FC_SUCCESS)
      return rc;
    for (j = 0; j < numVert; j++) 
      counts[dsSlot->fileInfo.dynaInfo.vertIDs[meshID][j]]++;
  }
  free(meshes);

  // look at counts
  numShare = 0;
  for (i = 0; i < dsSlot->fileInfo.dynaInfo.numTotalVert; i++) {
    if (counts[i] > 1)
      numShare++;
    else if (counts[i] < 1)
      printf("***** ieeee ***\n");
  }
  return numShare;
}
