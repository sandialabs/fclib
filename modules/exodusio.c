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
 * \file exodusio.c
 * \brief Implementation for \ref ExodusFileIO module.
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/modules/exodusio.c,v $
 * $Revision: 1.53 $
 * $Date: 2007/08/29 07:03:58 $
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
#include "subset.h"
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
#include "exodusII_ext.h"
#include "exodusioP.h"

/**
 * \addtogroup ExodusFileIO
 * \brief Exodus II file IO.
 *
 * \description
 *   
 *   Be aware that if you open copy, and rewrite an exodus file
 *   using FCLib, THE BLOCK/NODESET/SIDESET IDS MAY CHANGE!
 *  
 *   Exodus restricts the lengths of names of objects to 32 chars.
 *
 *   This file format stores a single global vertex array, and all vertex
 *   associated variables (nodal variables) are stored on that single vertex
 *   array instead of per mesh as they are in FCLib. When the file is written,
 *   if a nodal variable does not exist on every mesh, the global array will be
 *   padded with zeros.
 *
 *   Only recently has Exodus II written arbitrary edge and face
 *   We still currently stick with the convention of writing
 *   "sidesets" which will be edge subsets for 2D topologies and face
 *   subsets for 3D topologies. 
 *
 *   Exodus II has global node, element, and sides sets. This reader
 *   breaks the subsets into per mesh subsets. They are not put back
 *   together again upon writing. The optional orientation array
 *   for sets not read in (and thus not written out).
 *
 *   Exodus II does not support multicomponent variables or attributes.
 *   When reading and writing, each multicomponent variable or attribute
 *   is broken into components and renamed varname_x,
 *   varname_y, varname_z for vectors or varname_c1, varname_c2 ... for 
 *   tensors.
 *
 *   This reader does not support time-changing mesh geometry.
 *
 *   During reading, Exodus info, records, number maps, order
 *   maps, and distribution factors are discarded. 
 *
 *   - Exodus variables are seq variables with FC_DT_DOUBLE (when writing, 
 *   other data types are converted to double).
 *   - Exodus properties are FC_AT_WHOLE_MESH and FC_DT_INT (when writing,
 *   other data types of FC_AT_WHOLE are discarded). NOTE: 'ID' is a name
 *   already used by Exodus.
 *   - Exodus attributes are non seq variables. We support these for 
 *   FC_AT_ELEMENT and FC_AT_VERTEX as element and nodal attributes with
 *   data type FC_DT_DOUBLE (when writing, all other data types on elements
 *   are converted to double, and all other assoc types
 *   are discarded).
 */

/**
 * \ingroup ExodusFileIO
 * \brief Compare algorithm for sorting seq vars
 *
 * \description
 *
 *    Compare two sequence variables (using only the first step).
 *
 *    Most important thing is that we sort first by association type.  Will
 *    also sort by name, numComp, and mathType, and dataType (numDataPoint can
 *    be ignored).
 *    
 *    Only used in this file -- helper for writing files. Sort by assoc first
 *    since different types treated much differently by exodus. 
 *
 * \todo
 *   - ACG is using this for non-seq vars as well. I may rename this.
 *     It is used in writing out non-seq-vars as attributes.
 *
 * \modifications
 *    - 9/14/2005 WSD Created.
 *    - 8/11/2006 WSD. Changed to use seqVars directly instead of via
 *        a linked list node object.
 */
int _fc_exo_seqVar_cmp(const void* a, const void* b) {
  const FC_Variable* seqVar_a = (const FC_Variable*)a;
  const FC_Variable* seqVar_b = (const FC_Variable*)b;
  _FC_VarSlot *varSlot_a = _fc_getVarSlot(seqVar_a[0]);
  _FC_VarSlot *varSlot_b = _fc_getVarSlot(seqVar_b[0]);
  
  // Compare by association type
  if (varSlot_a->assoc > varSlot_b->assoc)
    return 1;
  else if (varSlot_a->assoc < varSlot_b->assoc)
    return -1;
  else {

    // Compare by names
    // (Note: strcmp is not restricted to returning 1 and -1)
    if (strcmp(varSlot_a->header.name, varSlot_b->header.name) > 0)
      return 1;
    else if (strcmp(varSlot_a->header.name, varSlot_b->header.name) < 0)
      return -1;
    else {

      // Compare by numComp
      if (varSlot_a->numComponent > varSlot_b->numComponent)
        return 1;
      else if (varSlot_a->numComponent < varSlot_b->numComponent)
        return -1;
      else {

        // Compare by mathType
        if (varSlot_a->mathtype > varSlot_b->mathtype)
          return 1;
        else if (varSlot_a->mathtype < varSlot_b->mathtype)
          return -1;
        else {
          
          // Compare by dataType
          if (varSlot_a->datatype > varSlot_b->datatype)
            return 1;
          else if (varSlot_a->datatype < varSlot_b->datatype)
            return -1;
          else 
            return 0;
        }
      }
    }
  }
}

/**
 * \ingroup ExodusFileIO
 * \brief Compare algorithm for sorting seq vars including the meshID.
 *
 * \description
 *
 *    Compare two sequence variables (using only the first step).
 *
 *    Most important thing is that we sort first by association type and
 *    last by which mesh (i.e. exodus meshID). Will also sort by name,
 *    numComp, and mathType, and dataType (numDataPoint can be ignored).
 *    
 *    Only used in this file -- helper for writing files. Sort by assoc first
 *    since different types treated much differently by exodus.  Sort by meshID
 *    last because that is easiest way to see the groupings over all
 *    meshes. The order of other things that denote different variables doesn't
 *    matter.
 *
 * \modifications
 *    - 9/2/2005 WSD Created.
 *    - 9/14/2005 WSD. Core code moved to simpler compare. This is a more
 *        complex wrapper.
 *    - 8/11/2006 WSD. Changed to use seqVars directly instead of via
 *        a linked list node object.
 */
int _fc_exo_seqVar_meshID_cmp(const void* a, const void* b) {
  const FC_Variable* seqVar_a = (const FC_Variable*)a;
  const FC_Variable* seqVar_b = (const FC_Variable*)b;
  _FC_VarSlot *varSlot_a = _fc_getVarSlot(seqVar_a[0]);
  _FC_VarSlot *varSlot_b = _fc_getVarSlot(seqVar_b[0]);
  int cmp;
  
  // Compare everything else first
  cmp = _fc_exo_seqVar_cmp(a, b);
  if (cmp != 0)
    return cmp;
  else {
    // Compare by meshID
    _FC_MeshSlot *meshSlot_a = _fc_getMeshSlot(varSlot_a->mesh);
    _FC_MeshSlot *meshSlot_b = _fc_getMeshSlot(varSlot_b->mesh);
    if (meshSlot_a->fileInfo.exoInfo.meshID >
        meshSlot_b->fileInfo.exoInfo.meshID) 
      return 1;
    else if (meshSlot_a->fileInfo.exoInfo.meshID < 
             meshSlot_b->fileInfo.exoInfo.meshID)
      return -1;
    else
      return 0;
  }
}

/**
 * \ingroup ExodusFileIO
 * \brief Create the name for a component of a variable to use in Exodus file.
 *
 * \description
 *
 *   Given a name, the number of components, the math type and the component
 *   index, return a new string that is a name for that component.
 *
 *   - If it's a scalar: comp_name = name
 *   - If it's a vector: comp_name = name_x or name_y or name_z
 *   - If it's any other multi-component = name_c0 or name_c1 or ...
 *
 *   This functions returns malloc memory that must be freed by the caller.
 *
 * \todo Promote this to the fileio module so that all IO writers can
 *    use the same naming conventions.
 *
 * \modifications
 *   - 2006-11-27 WSK. Created because common code was cropping up
 *        in multiple places.
 *   - 10-11-1007 ACG restrict to MAX_STRLEN length
 *   
 */
FC_ReturnCode _fc_createExoVarCompName(
  char* name,      /**< input - the name of the variable */
  int numComp,     /**< input - the number of components of the variable */
  FC_MathType mathType, /**< input - the math type of the variable */
  int compID,      /**< input - the index of this particular component */
  char** compName  /**< output - the created component name (must be freed
                      by the caller). */
) 
{
  char *tempName;
  char vector_endings[3][2] = { "x", "y", "z" }; 

  // default return
  if (compName)
    *compName = NULL;

  // check input
  if (!name || numComp < 1 || compID < 0 || compID >= numComp || 
      !fc_isMathTypeValid(mathType) || mathType == FC_MT_UNKNOWN  || !compName) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // make space
  
  tempName = (char*)malloc(MAX_STR_LENGTH*sizeof(char));
  if (!tempName) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }

  // do it
  if (numComp == 1){
    strncpy(tempName, name, MAX_STR_LENGTH-1);
    tempName[MAX_STR_LENGTH-1] = '\0';
  }
  else {
    char clippedname[MAX_STR_LENGTH-4];
    snprintf(clippedname, MAX_STR_LENGTH-4,name);
    if (numComp < 4 && mathType == FC_MT_VECTOR){
      sprintf(tempName, "%s_%s", clippedname, vector_endings[compID]);
    }  else {
      sprintf(tempName, "%s_c%d", clippedname, compID);
    }
  }

  *compName = tempName;
  return FC_SUCCESS;
}

/**
 * \ingroup ExodusFileIO
 * \brief Return LU map for local to exodus global nodal IDs
 *
 * \description
 *  
 *    To use, index into the lookup table with the vertex ID. The returned
 *    value is the index you work use to access a global array of nodal
 *    values. Add 1 to get the actual Exodus ID.
 *
 *    Do NOT free the returned array.
 *
 *    Only works if the owning dataset is tied to an Exodus file.
 *
 * \modifications
 *    - 2006-01-23 WSD Created.
 */
FC_ReturnCode _fc_getMeshExodusGlobalNodalIDsPtr(
  FC_Mesh mesh,   /**< input - the mesh */
  int** LUTable   /**< output - ptr to local to global ID map,
		     Do NOT free! */
) {
  _FC_DsSlot* dsSlot;
  _FC_MeshSlot* meshSlot;
  
  // default return
  if (LUTable)
    *LUTable = NULL;

  // check input
  if (!fc_isMeshValid(mesh) || !LUTable) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  meshSlot = _fc_getMeshSlot(mesh);
  dsSlot = _fc_getDsSlot(meshSlot->ds);
  if (dsSlot->fileType != FC_FT_EXODUS) {
    fc_printfErrorMessage("Dataset is not of type %s",
			  fc_getFileIOTypeText(FC_FT_EXODUS));
    return FC_INPUT_ERROR;
  }

  // do it
  *LUTable = dsSlot->fileInfo.exoInfo.vertIDs[meshSlot->fileInfo.exoInfo.meshID];

  return FC_SUCCESS;
}

/**
 * \ingroup ExodusFileIO
 * \brief Return offset for local to exodus global element IDs
 *
 * \description
 *  
 *    To use, add the offset to the element ID. This value is the index you
 *    would use to access a global array of element values. Add 1 to get the
 *    actual Exodus ID.
 *
 *    Only works if the owning dataset is tied to an Exodus file.
 
 * \modifications
 *    - 2006-01-23 WSD Created.
 */
FC_ReturnCode _fc_getMeshExodusGlobalElementIDOffset(
  FC_Mesh mesh,   /**< input - the mesh */
  int* offset    /**< output - offset */
) {
  _FC_DsSlot* dsSlot;
  _FC_MeshSlot* meshSlot;
  
  // default return
  if (offset)
    *offset = -1;

  // check input
  if (!fc_isMeshValid(mesh) || !offset) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  meshSlot = _fc_getMeshSlot(mesh);
  dsSlot = _fc_getDsSlot(meshSlot->ds);
  if (dsSlot->fileType != FC_FT_EXODUS) {
    fc_printfErrorMessage("Dataset is not of type %s",
			  fc_getFileIOTypeText(FC_FT_EXODUS));
    return FC_INPUT_ERROR;
  }

  // do it
  *offset = meshSlot->fileInfo.exoInfo.elemOffset;

  return FC_SUCCESS;
}

/**
 * \ingroup ExodusFileIO
 * \brief Initialize for Exodus IO
 *
 * \description
 *
 *    Uses the setting of the library verbosity at the time of call.
 *
 * \modifications
 *    - 8/2/2005 WSD Created
 */
FC_ReturnCode _fc_initExodus(void) 
{
  // log message
  fc_printfLogMessage("Initializing Exodus II");

  // Set exodus messages
  switch (fc_getLibraryVerbosity()) {
  case FC_QUIET :             ex_opts(0); break;
  case FC_ERROR_MESSAGES :    // fall through
  case FC_WARNING_MESSAGES :  ex_opts(EX_VERBOSE); break;
  case FC_LOG_MESSAGES :      // fall through
  case FC_DEBUG_MESSAGES :    ex_opts(EX_VERBOSE | EX_DEBUG); break;
  }

  // can't check for errors
  return FC_SUCCESS;
}


/**
 * \ingroup ExodusFileIO
 * \brief Load an Exodus II dataset from a file.
 *
 * \description
 *
 *    Given a name of an Exodus II file, this routine loads the
 *    dataset and makes all of its members available for use.
 *    Members are loaded lazily: that is all of their metadata 
 *    is loaded immediately, but the big data arrays are not 
 *    loaded until needed (this happens automatically).
 *
 *  \todo
 *   - make elem sets lazy load
 *   - should we automatically take out the exodus ID named items (props)?
 *
 * \modifications
 *    - 8/1/05 WSD Created.
 *    - 9/4/07 ACG. Commencing reading elem and nodal attr
 *    - 10/3/07 ACG Nodal attr
 *    - 10/X/07 ACG Elem sets
 *    - 10/17/07 ACG removed names support
 */
FC_ReturnCode _fc_loadExodusDataset(
  char* filename,      /**< input - name of Exodus dataset to open. */
  FC_Dataset *dataset  /**< output - handle for this dataset. */
) {
  FC_ReturnCode rc;
  int i, j, k, m;
  _FC_DsSlot* dsSlot;
  _FC_SeqSlot* seqSlot;
  int error, exoid, cpu_word_size, io_word_size, *blk_ids;
  float version;
  FC_Sequence sequence;
  FC_Mesh *meshes;
  int numMesh, numTotalVert, numTotalElem, numDim, numNodeSet, numSideSet, numElemSet;
  int numStep, *numVerts, *numElems, numProp;
  int *numElemAttr;
  int numNodalAttr;
  int **hasVerts;
  char title[MAX_LINE_LENGTH+1], namebuf[MAX_STR_LENGTH];
  char *temp_name;
  float fdum;
  char* cdum;
  int numVarType = 3;
  char* exoVarTypes[3] = { "g", "n", "e" };
  ex_init_params exoParams;
  FC_AssociationType fcVarTypes[3] = { FC_AT_WHOLE_DATASET, 
                                       FC_AT_VERTEX, 
                                       FC_AT_ELEMENT };

  // default return value
  if (dataset != NULL)
    *dataset = FC_NULL_DATASET;

  // check input
  if (!filename || !dataset) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // FIX? move this from all io helpers into fileio?
  // make sure library is ready, return error if it's not
  if (!_fc_isLibraryInit()) {
    fc_printfErrorMessage("Cannot load a dataset until library is initialized.");
    return FC_ERROR;
  }

  // log message
  fc_printfLogMessage("Loading Exodus dataset '%s'", filename);

  // Open the file
  cpu_word_size = 8;
  io_word_size = 0;
  exoid = ex_open(filename, EX_READ, &cpu_word_size, &io_word_size, &version);
  if (exoid < 0) {
    fc_printfErrorMessage("Exodus failed to open database '%s'", filename);
    return FC_FILE_IO_ERROR;
  }


  // Get file parameters
  error = ex_get_init_ext(exoid, &exoParams);
  if (error != 0) {
    fc_printfErrorMessage("Exodus failed to read file parameters");
    return FC_FILE_IO_ERROR;
  }
  strncpy(title, exoParams.title, MAX_STR_LENGTH);
  numDim = exoParams.num_dim;
  numTotalVert = exoParams.num_nodes;
  numTotalElem = exoParams.num_elem;
  numMesh = exoParams.num_elem_blk;
  numNodeSet = exoParams.num_node_sets;
  numSideSet = exoParams.num_side_sets;
  numElemSet = exoParams.num_elem_sets;

  // create a dataset
  rc = fc_createDataset(title, dataset);
  if (rc != FC_SUCCESS) {
    fc_printfErrorMessage("Failed to create dataset");
    return rc;
  }
  // twiddle writable flag so we don't acciently clobber
  dsSlot = _fc_getDsSlot(*dataset);
  dsSlot->header.committed = 1;
  dsSlot->writable = 0;
  dsSlot->fileType = FC_FT_EXODUS;
  dsSlot->baseFile = (char*)malloc((strlen(filename)+1)*sizeof(char));
  if (!dsSlot->baseFile) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  strcpy(dsSlot->baseFile, filename);
  dsSlot->fileInfo.exoInfo.exoid = exoid;
  dsSlot->fileInfo.exoInfo.cpu_word_size = cpu_word_size;
  dsSlot->fileInfo.exoInfo.io_word_size = io_word_size;
  dsSlot->fileInfo.exoInfo.numTotalVert = numTotalVert;
  dsSlot->fileInfo.exoInfo.numMesh = numMesh;
  dsSlot->fileInfo.exoInfo.vertIDs = (int**)malloc(numMesh*sizeof(int*));
  if (!dsSlot->fileInfo.exoInfo.vertIDs) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }

  // Create meshes
  meshes = (FC_Mesh*)malloc(numMesh*sizeof(FC_Mesh));
  numVerts = (int*)malloc(numMesh*sizeof(int));
  numElems = (int*)malloc(numMesh*sizeof(int));
  numElemAttr = (int*)malloc(numMesh*sizeof(int));
  blk_ids = (int*)malloc(numMesh*sizeof(int));
  hasVerts = (int**)malloc(numMesh*sizeof(int**));
  if (!meshes || !numVerts || !numElems || !numElemAttr ||
      !blk_ids || !hasVerts) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    free(meshes);
    free(blk_ids);
   return FC_MEMORY_ERROR;
  }
  for (i = 0; i < numMesh; i++) {
    hasVerts[i] = (int*)malloc(numTotalVert*sizeof(int*));
    if (!hasVerts[i]) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      for (j = i - 1; j >= 0; j++)
        free(hasVerts[j]);
      free(meshes); free(blk_ids); free(hasVerts);
      free(numElems); free(numElemAttr);
      return FC_MEMORY_ERROR;
    }
  }
  error = ex_get_elem_blk_ids(exoid, blk_ids);
  if (error != 0) {
    fc_printfErrorMessage("Exodus failed to get elem blk ids");
    return FC_FILE_IO_ERROR;
  }
  numTotalElem = 0;  // keep a running count
  for (i = 0; i < numMesh; i++) {
    int *conns, *vertIDs;
    int numElem, numVert, numVertPerElem;
    char exoElemType[MAX_STR_LENGTH+1];
    FC_ElementType elemType;
    _FC_MeshSlot* meshSlot;

    // Read & process meta data
    error = ex_get_elem_block(exoid, blk_ids[i], exoElemType, &numElem,
                              &numVertPerElem, &numElemAttr[i]);
    if (error != 0) {
      fc_printfErrorMessage("Exodus failed to get elem block meta data");
      return FC_FILE_IO_ERROR;
    }
    if (numVertPerElem == 1)
      elemType = FC_ET_POINT;
    else if (numVertPerElem == 2)
      elemType = FC_ET_LINE;
    else if (numVertPerElem == 3 && !strncasecmp(exoElemType, "tri", 3))
      elemType = FC_ET_TRI;
    else if (numVertPerElem == 4 && !strncasecmp(exoElemType, "quad", 4))
      elemType = FC_ET_QUAD;
    else if (numVertPerElem == 4 && !strncasecmp(exoElemType, "tet", 3))
      elemType = FC_ET_TET;
    else if (numVertPerElem == 5)
      elemType = FC_ET_PYRAMID;
    else if (numVertPerElem == 6 && 
             (!strncasecmp(exoElemType, "prism", 5) || 
              !strncasecmp(exoElemType, "wedge", 5)))
      elemType = FC_ET_PRISM;
    else if (numVertPerElem == 8 && !strncasecmp(exoElemType, "hex", 3))
      elemType = FC_ET_HEX;
    else {
      fc_printfErrorMessage("Cannot handle exodus element type of '%s' "
                            "with %d verts per element", exoElemType,
                            numVertPerElem);
      return FC_ERROR;
    }

    // read & process conns
    conns = (int*)malloc(numElem*numVertPerElem*sizeof(int));
    if (!conns) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      return FC_MEMORY_ERROR;
    }      
    error = ex_get_elem_conn(exoid, blk_ids[i], conns);
    if (error != 0) {
      fc_printfErrorMessage("Exodus failed to get elem conns");
      return FC_FILE_IO_ERROR;
    }
    for (j = 0; j < numTotalVert; j++)
      hasVerts[i][j] = -1;
    for (j = 0; j < numElem*numVertPerElem; j++) 
      hasVerts[i][conns[j]-1] = 1; //exodus numbering starts at 1
    free(conns);
    numVert = 0;
    for (j = 0; j < numTotalVert; j++) {
      if (hasVerts[i][j] > -1) {
        hasVerts[i][j] = numVert; 
        numVert++;
      }
    }
    // hasVerts [mesh][exonumber-1] = either the vert id on this mesh (uniquely numbered
    // starting at zero) for all those participating verts in the mesh, -1 otherwise.
    vertIDs = (int*)malloc(numVert*sizeof(int));
    if (!vertIDs) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      return FC_MEMORY_ERROR;
    }
    for (j = 0; j < numTotalVert; j++)
      if (hasVerts[i][j] > -1) 
        vertIDs[hasVerts[i][j]] = j; 
    // vertIDs[local_per_mesh_number(beginning at zero)] = exonumber
    numVerts[i] = numVert;
    numElems[i] = numElem;

    // Create mesh & save info
    // Try to get exo name
    namebuf[0] = '\0';
    error = ex_get_name(exoid, EX_ELEM_BLOCK, blk_ids[i], namebuf);
    if (error != 0) {
      fc_printfErrorMessage("Exodus failed to get elem block name");
      return FC_FILE_IO_ERROR;
    }
    if (strcmp(namebuf,"")) {
      temp_name = namebuf;
    } else {
      temp_name = (char*)malloc(20*sizeof(char));
      sprintf(temp_name, "block_%d", blk_ids[i]);
    }
    rc = fc_createMesh(*dataset, temp_name, &meshes[i]);
    if (rc != FC_SUCCESS) {
      fc_printfErrorMessage("Failed to create mesh #%d '%s'", i, temp_name);
      return rc;
    }
    if (!strcmp(namebuf, "")) 
      free(temp_name);
    meshSlot = _fc_getMeshSlot(meshes[i]);
    meshSlot->fileInfo.id = blk_ids[i];
    meshSlot->fileInfo.exoInfo.meshID = i;
    //offset to go from fclib id to exodus id - num Elem up til this mesh
    meshSlot->fileInfo.exoInfo.elemOffset = numTotalElem;
    numTotalElem += numElem;
    dsSlot->fileInfo.exoInfo.vertIDs[i] = vertIDs;
    meshSlot->dim = numDim;
    meshSlot->numVertex = numVert;
    meshSlot->numElement = numElem;
    meshSlot->elemType = elemType;
    meshSlot->topodim = fc_getElementTypeTopoDim(elemType);
    meshSlot->header.committed = 1;
  }
 
  // Create node sets
  // NOTE: becuase this cannot be lazily loaded, the committed flag is
  // not set. This is so that this data will not be released in a release
  // call since it cannot be reloaded. This has the side effect though of 
  // allowing this subset to be altered from what is on the disk
  if (numNodeSet > 0) {
    int* nodeset_ids = (int*)malloc(numNodeSet*sizeof(int));
    if (!nodeset_ids) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      free(meshes);
      return FC_MEMORY_ERROR;
    }
    error = ex_get_node_set_ids(exoid, nodeset_ids);
    if (error != 0) {
      fc_printfErrorMessage("Exodus failed to get elem blk ids");
      return FC_FILE_IO_ERROR;
    }
    for (i = 0; i < numNodeSet; i++) {
      int numNode, numDistFactor, *nodes;
      
      // Read & process members
      error = ex_get_node_set_param(exoid, nodeset_ids[i], &numNode,
                                    &numDistFactor);
      if (error != 0) {
        fc_printfErrorMessage("Exodus failed to get node set param");
        free(meshes); free(nodeset_ids);
        return FC_FILE_IO_ERROR;
      }
      nodes = (int*)malloc(numNode*sizeof(int));
      if (!nodes) {
        fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
        free(meshes); free(nodeset_ids);
        return FC_MEMORY_ERROR;
      }
      error = ex_get_node_set(exoid, nodeset_ids[i], nodes);
      if (error != 0) {
        fc_printfErrorMessage("Exodus failed to get nodes in node set");
        free(meshes); free(nodeset_ids); free(nodes);
        return FC_FILE_IO_ERROR;
      }
      
      // Create subset & save info & data
      // (Each node set might span meshes and will be multiple subsets)
      // Have to save too much info or do too much work for lazy read,
      // so, do it all here
      namebuf[0] = '\0';
      error = ex_get_name(exoid, EX_NODE_SET, nodeset_ids[i], namebuf);
      if (error != 0) {
        fc_printfErrorMessage("Exodus failed to get node set name");
        return FC_FILE_IO_ERROR;
      }
      if (strcmp(namebuf, "")) {
        temp_name = namebuf;
      } else {
        temp_name = (char*)malloc(20*sizeof(char));
        sprintf(temp_name, "nodeset_%d", nodeset_ids[i]);
      }

      //split the subset into per-mesh subsets containing only the ids that
      //exist on that mesh
      for (j = 0; j < numMesh; j++) {
        FC_Subset subset = FC_NULL_SUBSET;
        for (k = 0; k < numNode; k++) {
          if (hasVerts[j][nodes[k]-1] > -1) {
            if (FC_HANDLE_EQUIV(subset, FC_NULL_SUBSET)) {
              rc = fc_createSubset(meshes[j],  temp_name, FC_AT_VERTEX,
                                   &subset);
              if (rc != FC_SUCCESS) {
                fc_printfErrorMessage("Failed to create nodeset #%d '%s' " 
                                      " on mesh #%d", i, temp_name, j);
                return rc;
              }
            }
            fc_addMemberToSubset(subset, hasVerts[j][nodes[k]-1]);
          }
        }
      }
      free(nodes);
      if (!strcmp(namebuf, ""))
        free(temp_name);
    } //loop over nodesets
    free(nodeset_ids);
  } //if numNodeSet


  // Create elem sets. not doing this as lazy load since
  // we have to process through all the elements anyway (see
  // how the count is retained in sideset). This means
  // we cant release this data -- see if that becomes
  // desireable
  // NOTE: becuase this cannot be lazily loaded, the committed flag is
  // not set. This is so that this data will not be released in a release
  // call since it cannot be reloaded. This has the side effect though of 
  // allowing this subset to be altered from what is on the disk
  if (numElemSet > 0) {
    int* elemset_ids = (int*)malloc(numElemSet*sizeof(int));
    if (!elemset_ids) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      free(meshes);
      return FC_MEMORY_ERROR;
    }
    error = ex_get_ids(exoid, EX_ELEM_SET, elemset_ids);
    if (error != 0) {
      fc_printfErrorMessage("Exodus failed to get elem set ids");
      return FC_FILE_IO_ERROR;
    }
    for (i = 0; i < numElemSet; i++) {
      int numElem, numDistFactor, *elems;
      
      // Read & process members
      error = ex_get_set_param(exoid, EX_ELEM_SET, elemset_ids[i], &numElem,
                                    &numDistFactor);
      if (error != 0) {
        fc_printfErrorMessage("Exodus failed to get elem set param");
        free(meshes); free(elemset_ids);
        return FC_FILE_IO_ERROR;
      }
      elems = (int*)malloc(numElem*sizeof(int));
      if (!elems) {
        fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
        free(meshes); free(elemset_ids);
        return FC_MEMORY_ERROR;
      }
      error = ex_get_set(exoid, EX_ELEM_SET, elemset_ids[i], elems, NULL);
      if (error != 0) {
        fc_printfErrorMessage("Exodus failed to get elems in elem set");
        free(meshes); free(elemset_ids); free(elems);
        return FC_FILE_IO_ERROR;
      }
      
      // Create subset & save info & data
      // (Each elem set might span meshes and will be multiple subsets)
      namebuf[0] = '\0';
      error = ex_get_name(exoid, EX_ELEM_SET, elemset_ids[i], namebuf);
      if (error != 0) {
	fc_printfErrorMessage("Exodus failed to get elem set name");
	return FC_FILE_IO_ERROR;
      }
      if (strcmp(namebuf, "")) {
        temp_name = namebuf;
      } else {
        temp_name = (char*)malloc(20*sizeof(char));
        sprintf(temp_name, "elemset_%d", elemset_ids[i]);
      }

      //split the subset into per-mesh subsets containing only the ids that
      //exist on that mesh
      for (j = 0; j < numMesh; j++) {
        FC_Subset subset = FC_NULL_SUBSET;
	_FC_MeshSlot* meshSlot = _fc_getMeshSlot(meshes[j]);
        for (k = 0; k < numElem; k++) {
	  if (elems[k] - 1 >=  meshSlot->fileInfo.exoInfo.elemOffset &&
	      elems[k] - 1 <  meshSlot->fileInfo.exoInfo.elemOffset + meshSlot->numElement) {
            if (FC_HANDLE_EQUIV(subset, FC_NULL_SUBSET)) {
              rc = fc_createSubset(meshes[j],  temp_name, FC_AT_ELEMENT,
                                   &subset);
              if (rc != FC_SUCCESS) {
                fc_printfErrorMessage("Failed to create elemset #%d '%s' " 
                                      " on mesh #%d", i, temp_name, j);
                return rc;
              }
            }
            fc_addMemberToSubset(subset,elems[k]-1-meshSlot->fileInfo.exoInfo.elemOffset);
	  } //elem on this mesh
	}
      }
      free(elems);
      if (!strcmp(namebuf, ""))
        free(temp_name);
    } //loop over elemsets
    free(elemset_ids);
  } //if numElemSet

  
  // Create side sets
  if (numSideSet > 0) {
    int* sideset_ids = (int*)malloc(numSideSet*sizeof(int));
    if (!sideset_ids) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      free(meshes);
      return FC_MEMORY_ERROR;
    }
    error = ex_get_side_set_ids(exoid, sideset_ids);
    if (error != 0) {
      fc_printfErrorMessage("Exodus failed to get elem blk ids");
      return FC_FILE_IO_ERROR;
    }
    for (i = 0; i < numSideSet; i++) {
      int numSide, numDistFactor, *elems, *sides;
      FC_Subset subsets[numMesh];
      int elemOffsets[numMesh];
      
      // Read & process members
      error = ex_get_side_set_param(exoid, sideset_ids[i], &numSide,
                                    &numDistFactor);
      if (error != 0) {
        fc_printfErrorMessage("Exodus failed to get sides in side set");
        free(meshes); free(sideset_ids);
        return FC_FILE_IO_ERROR;
      }
      elems = (int*)malloc(numSide*sizeof(int));
      sides = (int*)malloc(numSide*sizeof(int));
      if (!elems || !sides) {
        fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
        free(meshes); free(sideset_ids); free(elems); free(sides);
        return FC_MEMORY_ERROR;
      }
      error = ex_get_side_set(exoid, sideset_ids[i], elems, sides);
      if (error != 0) {
        fc_printfErrorMessage("Exodus failed to get sides in side set");
        free(meshes); free(sideset_ids); free(elems); free(sides);
        return FC_FILE_IO_ERROR;
      }
      free(sides); // don't need this right now
      
      // Create subset & save info
      // (Each side set might span meshes and will be multiple subsets)
      // Unlike node sets, will be lazily read.
      namebuf[0] = '\0';
      error = ex_get_name(exoid, EX_SIDE_SET, sideset_ids[i], namebuf);
      if (error != 0) {
        fc_printfErrorMessage("Exodus failed to get side set name");
        return FC_FILE_IO_ERROR;
      }
      if (strcmp(namebuf, "")) {
        temp_name = namebuf;
      } else {
        temp_name = (char*)malloc(20*sizeof(char));
        sprintf(temp_name, "sideset_%d", sideset_ids[i]);
      }
      for (j = 0; j < numMesh; j++) { 
        subsets[j] = FC_NULL_SUBSET;
        elemOffsets[j] = _fc_getMeshSlot(meshes[j])->fileInfo.exoInfo.elemOffset;
      }
      for (j = 0; j < numSide; j++) {
        _FC_SubSlot* subSlot;
        // Figure out which mesh this side is on  - compare exo id to offset which is the running count of the total up till that point
        for (k = numMesh-1; k >= 0; k--) 
          if (elems[j] - 1 >= elemOffsets[k])
            break;
        subSlot = _fc_getSubSlot(subsets[k]);
        // If first time we've seen a side for this mesh, make a subset
        if (!subSlot) {
          FC_AssociationType assoc;
          FC_ElementType elemType;
          int topoDim;

          rc = fc_getMeshElementType(meshes[k], &elemType);
          if (rc != FC_SUCCESS)
            return rc;
          topoDim = fc_getElementTypeTopoDim(elemType);
          if (topoDim == 3)
            assoc = FC_AT_FACE;
          else if (topoDim == 2)
            assoc = FC_AT_EDGE;
          else {
            fc_printfErrorMessage("Encountered unexpected element topology");
            return FC_ERROR;
          }
	  //create the subset but dont fill in the ids
          rc = fc_createSubset(meshes[k],  temp_name, assoc, &subsets[k]);
          if (rc != FC_SUCCESS) {
            fc_printfErrorMessage("Failed to create sideset #%d '%s' " 
                                  " on mesh #%d", i, temp_name, j);
            return rc;
          }
          subSlot = _fc_getSubSlot(subsets[k]);
          subSlot->fileInfo.exoInfo.assoc_flag = 1;
          subSlot->fileInfo.exoInfo.nsssid = sideset_ids[i];
	  //Because this data exists on disk and can be reloaded, set the committed flag
          subSlot->header.committed = 1;
        }
        // count the side
        subSlot->numMember++;
      }
      free(elems);
      if (!strcmp(namebuf, ""))
        free(temp_name);
    } // end of loop over side sets
    free(sideset_ids);
  }

  // Create vars
  // Get element block properties = FC_AT_WHOLE_MESH & FC_DT_INT vars
  error = ex_inquire(exoid, EX_INQ_EB_PROP, &numProp, &fdum, cdum);
  if (error != 0) {
    fc_printfErrorMessage("Exodus failed to get number of elem block props");
    return FC_FILE_IO_ERROR;
  }
  if (numProp > 0) {
    char** prop_names;
    prop_names = (char**)malloc(numProp*sizeof(char*));
    if (!prop_names) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      return FC_MEMORY_ERROR;
    }
    for (i = 0; i < numProp; i++) {
      prop_names[i] = (char*)malloc((MAX_STR_LENGTH)*sizeof(char));
      if (!prop_names[i]) {
        fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
        return FC_MEMORY_ERROR;
      }
    }
    error = ex_get_prop_names(exoid, EX_ELEM_BLOCK, prop_names);
    if (error != 0) {
      fc_printfErrorMessage("Exodus failed to get elem block prop names");
      return FC_FILE_IO_ERROR;
    }
    for (i = 0; i < numMesh; i++) {
      for (j = 0; j < numProp; j++) {
        int prop_value;
        FC_Variable var;
        error = ex_get_prop(exoid, EX_ELEM_BLOCK, blk_ids[i], prop_names[j],
                            &prop_value);
        if (error != 0) {
          fc_printfErrorMessage("Exodus failed to get prop '%s' for mesh #%d",
                                prop_names[j], i);
          return FC_FILE_IO_ERROR;
        }
        rc = fc_createVariable(meshes[i], prop_names[j], &var);
        if (rc != FC_SUCCESS)
          return rc;
	//Property data is set here and is not reloadable, so do not set committed flag
        rc = fc_setVariableData(var, 1, 1, FC_AT_WHOLE_MESH, FC_MT_SCALAR,
                                FC_DT_INT, (void*)&prop_value);
        if (rc != FC_SUCCESS)
          return rc;
      }
    }
    for (i = 0; i < numProp; i++)
      free(prop_names[i]);
    free(prop_names);
  }

  // Blk Attributes (currently only elem)
  for (i = 0; i < numMesh; i++){
    if (numElemAttr[i] > 0){
      char** attr_names = (char**)malloc(numElemAttr[i]*sizeof(char*));
      if (!attr_names){
        fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
        return FC_MEMORY_ERROR;
      }

      for (j = 0; j < numElemAttr[i]; j++){
	attr_names[j] = (char*)malloc((MAX_STR_LENGTH+1)*sizeof(char));
        if (!attr_names[j]){
	  fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
	  return FC_MEMORY_ERROR;
        }
      }

      error = ex_get_elem_attr_names(exoid, blk_ids[i], attr_names);
      if (error != 0){
	fc_printfErrorMessage("Exodus failed to get attribute names");
	return FC_FILE_IO_ERROR;
      }
      for (j = 0; j < numElemAttr[i]; j++){
      	//does not set var data
	FC_Variable var;
	_FC_VarSlot* varSlot;
	rc = fc_createVariable(meshes[i], attr_names[j], &var);
	if (rc != FC_SUCCESS)
	  return rc;
	varSlot = _fc_getVarSlot(var);
	strcpy(varSlot->fileInfo.exoInfo.varType, "e");
	//need attr index for reading data later (starts at 1)
	varSlot->fileInfo.exoInfo.attrIndex = j + 1;
	varSlot->stepID = -1;
        varSlot->numDataPoint = numElems[i];
	varSlot->numComponent = 1;
	varSlot->assoc = FC_AT_ELEMENT;
	varSlot->mathtype = FC_MT_SCALAR;
	varSlot->datatype = FC_DT_DOUBLE;
	//EB attributes exist on disk and are lazily loaded, so set committed flag
	varSlot->header.committed = 1;

	free(attr_names[j]);
      }
      free(attr_names);
    }
  }

  free(numElemAttr);
  free(blk_ids);


  //nodal attr (exist on all meshes)
  numNodalAttr = 0;
  error = ex_get_attr_param(exoid, EX_NODAL, 0, &numNodalAttr);
  if (error != 0){
    fc_printfErrorMessage("Exodus failed to get num nodal attributes");
    return FC_FILE_IO_ERROR;
  }
  
  if (numNodalAttr > 0){
    char** attr_names = (char**)malloc(numNodalAttr*sizeof(char*));
    if (!attr_names){
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      return FC_MEMORY_ERROR;
    }
    for (i = 0; i < numNodalAttr; i++){
      attr_names[i] = (char*)calloc((MAX_STR_LENGTH+1), sizeof(char));
      if (!attr_names[i]){
	fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
	return FC_MEMORY_ERROR;
      }
    }
    error = ex_get_attr_names (exoid, EX_NODAL, 0, attr_names);
    if (error != 0){
      fc_printfErrorMessage("Exodus failed to get nodal attribute names");
      return FC_FILE_IO_ERROR;
    }

    for (i = 0; i < numNodalAttr; i++){
      for (j = 0; j < numMesh; j++){
      	//does not set var data
	FC_Variable var;
	_FC_VarSlot* varSlot;
	rc = fc_createVariable(meshes[j], attr_names[i], &var);
	if (rc != FC_SUCCESS)
	  return rc;
	varSlot = _fc_getVarSlot(var);
	strcpy(varSlot->fileInfo.exoInfo.varType, "n");
	//need attr index for reading data later (starts at 1)
	varSlot->fileInfo.exoInfo.attrIndex = i + 1;
	varSlot->stepID = -1;
        varSlot->numDataPoint = numVerts[j];
	varSlot->numComponent = 1;
	varSlot->assoc = FC_AT_VERTEX;
	varSlot->mathtype = FC_MT_SCALAR;
	varSlot->datatype = FC_DT_DOUBLE;
	//Nodal attributes exist on disk and are lazily loaded so set committed flag
	varSlot->header.committed = 1;
      }
      free(attr_names[i]);
    }
    free(attr_names);
  }

  
  // Create the sequence
  error = ex_inquire(exoid, EX_INQ_TIME, &numStep, &fdum, cdum);
  if (error != 0) {
    fc_printfErrorMessage("Exodus failed to get numStep");
    return FC_FILE_IO_ERROR;
  }
  rc = fc_createSequence(*dataset, "time", &sequence);
  if (rc != FC_SUCCESS) {
    fc_printfErrorMessage("Failed to create sequence");
    return rc;
  }
  seqSlot = _fc_getSeqSlot(sequence);
  seqSlot->numStep = numStep;
  seqSlot->datatype = FC_DT_DOUBLE;
  //seq exists on disk and is lazily loaded so set committed flag
  seqSlot->header.committed = 1;

  // Create Seq Vars
  // FIX? automatically stitch up vectors? could do this after creating --
  // use fc_mergeComponents()
  for (i = 0; i < numVarType; i++) {
    char** var_names;
    int numVar;
    int* numDatas;
    int* elem_truth_tab = NULL;

    // get number of seq vars of this type
    error = ex_get_var_param(exoid, exoVarTypes[i], &numVar);
    if (error != 0) {
      fc_printfErrorMessage("Exodus failed to get num of vars");
      return FC_FILE_IO_ERROR;
    }
    if (numVar < 1)
      continue; // done if there are none

    // Setup
    if (i == 1)
      numDatas = numVerts;
    else if (i == 2)
      numDatas = numElems;
    var_names = (char**)malloc(numVar*sizeof(char*));
    for (j = 0; j < numVar; j++)
      var_names[j] = (char*)malloc((MAX_STR_LENGTH+1)*sizeof(char));
    error = ex_get_var_names(exoid, exoVarTypes[i], numVar, var_names);
    if (error != 0) {
      fc_printfErrorMessage("Exodus failed to get '%s' var names",
                            exoVarTypes[i]);
      return FC_FILE_IO_ERROR;
    }
    if (i == 2) { //no truth tab for nodal
      elem_truth_tab = (int*)malloc(numMesh*numVar*sizeof(int));
      if (!elem_truth_tab) {
        fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
        return FC_MEMORY_ERROR;
      }
      error = ex_get_elem_var_tab(exoid, numMesh, numVar, elem_truth_tab);
      if (error != 0) {
      fc_printfErrorMessage("Exodus failed to read element truth table");
      return FC_FILE_IO_ERROR;
      }
    }

    // Make seq vars 
    for (j = 0; j < numVar; j++) {
      double* data = NULL;

      // Globals
      // (go ahead and get global data, probably easier not to do lazy)
      if (i == 0) {                                     
        FC_Variable *seqVar;
        rc = fc_createGlobalSeqVariable(*dataset, sequence, var_names[j],
                                        &numStep, &seqVar);
        if (rc != FC_SUCCESS)
          return rc;
        data = (double*)malloc(numStep*sizeof(double));
        error = ex_get_glob_var_time(exoid, j+1, 1, -1, data);
        if (error != 0) {
          fc_printfErrorMessage("Exodus failed to get global var '%s'",
                                var_names[j]);
          return FC_FILE_IO_ERROR;
        }
        for (m = 0; m < numStep; m++) {
	  //not lazily loaded so do not set committed flag
          rc = fc_setVariableData(seqVar[m], 1, 1, FC_AT_WHOLE_DATASET, 
				  FC_MT_SCALAR, FC_DT_DOUBLE, (void**)&data[m]);
          if (rc != FC_SUCCESS)
            return rc;
        }
        free(seqVar);
      }
      // Non-globals
      else {
        for (k = 0; k < numMesh; k++) {
          FC_Variable *seqVar;
          
          // skip zero entries in truth table (don't really exist)
          if (elem_truth_tab && elem_truth_tab[k*numVar+j] == 0)
            continue;
          //does not set var data
          rc = fc_createSeqVariable(meshes[k], sequence, var_names[j], &numStep,
                                    &seqVar);
          if (rc != FC_SUCCESS)
            return rc;
          for (m = 0; m < numStep; m++) {
            _FC_VarSlot* varSlot = _fc_getVarSlot(seqVar[m]);
            strcpy(varSlot->fileInfo.exoInfo.varType, exoVarTypes[i]);
            varSlot->fileInfo.exoInfo.varIndex = j + 1;
            varSlot->numDataPoint = numDatas[k];
            varSlot->numComponent = 1;
            varSlot->assoc = fcVarTypes[i];
            varSlot->mathtype = FC_MT_SCALAR;
            varSlot->datatype = FC_DT_DOUBLE;
	    //non-globals are laily loaded so set committed flag
            varSlot->header.committed = 1;
          }
          free(seqVar);
        } // end of loop over meshes
      }
      free(var_names[j]);
      free(data);
    } // end of loop over vars for type i
    free(elem_truth_tab);
    free(var_names);
  } // end of creating vars
  
  // cleanup
  free(meshes);
  free(numVerts);
  free(numElems);
  for (i = 0; i < numMesh; i++)
    free(hasVerts[i]);
  free(hasVerts);
  
  return FC_SUCCESS;
}

/**
 * \ingroup ExodusFileIO
 * \brief Read data from Exodus file into a sequence
 *
 * \modifications
 *    - 8/8/2005 WSD Created.
 */
FC_ReturnCode _fc_readExodusSequenceCoords(
  FC_Sequence sequence     /**< input - sequence to be read into */
) {
  int error;
  _FC_SeqSlot* seqSlot;
  _FC_DsSlot* dsSlot;
  double* times;

  // check input
  seqSlot = _fc_getSeqSlot(sequence);
  if (!seqSlot || !seqSlot->header.name) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Reading coords for sequence '%s'", seqSlot->header.name);
  
  // setup
  dsSlot = _fc_getDsSlot(seqSlot->ds);
  
  // get the sequence coords
  times = malloc(seqSlot->numStep*sizeof(double));
  error = ex_get_all_times(dsSlot->fileInfo.exoInfo.exoid, times);
  if (error != 0) {
    fc_printfErrorMessage("Exodus failed to get time values");
    return FC_FILE_IO_ERROR;
  }
  
  // don't have to do anything!
  seqSlot->coords = times;
  return FC_SUCCESS;
}

/**
 * \ingroup ExodusFileIO
 * \brief Read data from Exodus file into a mesh
 *
 * \modifications
 *    - 8/8/2005 WSD Created.
 */
FC_ReturnCode _fc_readExodusMeshCoords(
  FC_Mesh mesh     /**< input - mesh to be read into */
) {
  int error;
  int i, j;
  _FC_MeshSlot* meshSlot;
  _FC_DsSlot* dsSlot;
  int* vertIDs; 
  double *global_coords[3];

  // check input
  meshSlot = _fc_getMeshSlot(mesh);
  if (!meshSlot || !meshSlot->header.name) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Reading coords for mesh '%s'", meshSlot->header.name);

  // setup
  dsSlot = _fc_getDsSlot(meshSlot->ds);
  vertIDs = dsSlot->fileInfo.exoInfo.vertIDs[meshSlot->fileInfo.exoInfo.meshID];
  
  // get global coords
  for (i = 0; i < meshSlot->dim; i++)
    global_coords[i] = (double*)malloc(dsSlot->fileInfo.exoInfo.numTotalVert*sizeof(double));
  error = ex_get_coord(dsSlot->fileInfo.exoInfo.exoid, global_coords[0],
                       global_coords[1], global_coords[2]);
  if (error != 0) {
    fc_printfErrorMessage("Exodus failed to get global coords");
    return FC_FILE_IO_ERROR;
  }

  // Make local coords array
  meshSlot->coords = (double*)malloc(meshSlot->numVertex*meshSlot->dim*sizeof(double));
  if (!meshSlot->coords) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  for (i = 0; i < meshSlot->numVertex; i++)
    for (j = 0; j < meshSlot->dim; j++)
      meshSlot->coords[i*meshSlot->dim + j] = global_coords[j][vertIDs[i]];
  // vertIDs[local_per_mesh_number(beginning at zero)] = exonumber

  // cleanup
  for (i = 0; i < meshSlot->dim; i++)
    free(global_coords[i]);

  return FC_SUCCESS;
}

/**
 * \ingroup ExodusFileIO
 * \brief Read data from Exodus file into a mesh
 *
 * \modifications
 *    - 8/8/2005 WSD Created.
 */
FC_ReturnCode _fc_readExodusMeshElemConns(
  FC_Mesh mesh     /**< input - mesh to be read into */
) {
  int error;
  int i;
  _FC_MeshSlot* meshSlot;
  _FC_DsSlot* dsSlot;
  int numVertPerElem, numVert, numTotalVert;
  int* hasVert, *conns;

  // check input
  meshSlot = _fc_getMeshSlot(mesh);
  if (!meshSlot || !meshSlot->header.name) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Reading coords for mesh '%s'", meshSlot->header.name);

  // setup
  dsSlot = _fc_getDsSlot(meshSlot->ds);
  numTotalVert = dsSlot->fileInfo.exoInfo.numTotalVert;
  numVertPerElem = fc_getElementTypeNumVertex(meshSlot->elemType);
  
  // read global conns (conns for just this blk, but has global ids)
  conns = (int*)malloc(meshSlot->numElement*numVertPerElem*sizeof(int));
  if (!conns) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }      
  error = ex_get_elem_conn(dsSlot->fileInfo.exoInfo.exoid, 
                           meshSlot->fileInfo.id, conns);
  if (error != 0) {
    fc_printfErrorMessage("Exodus failed to get elem conns");
    return FC_FILE_IO_ERROR;
  }
  
  // Make a vertLU table
  hasVert = (int*)malloc(numTotalVert*sizeof(int));
  if (!hasVert) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  for (i = 0; i < numTotalVert; i++)
    hasVert[i] = -1;
  for (i = 0; i < meshSlot->numElement*numVertPerElem; i++)
    hasVert[conns[i]-1] = 1;
  numVert = 0;
  for (i = 0; i < numTotalVert; i++) {
    if (hasVert[i] > 0) {
      hasVert[i] = numVert;
      numVert++;
    }
  }
  if (numVert != meshSlot->numVertex) {
    fc_printfErrorMessage("Programmer error!");
    return FC_ERROR;
  }

  // Convert global conns to local conns
  for (i = 0; i < meshSlot->numElement*numVertPerElem; i++) 
    conns[i] = hasVert[conns[i]-1];
  meshSlot->elemToVertConns = conns;

  // cleanup
  free(hasVert);

  return FC_SUCCESS;
}

/**
 * \ingroup ExodusFileIO
 * \brief Read data from Exodus file into a subset.
 *
 * \description
 *   
 *    Lazy read only for side sets. Node and elem sets are handled on initial load.
 *
 * \modifications
 *    - 8/30/2005 WSD Created.
 */
FC_ReturnCode _fc_readExodusSubsetMembers(
  FC_Subset subset     /**< input - subset to be read into */
) {
  FC_ReturnCode rc;
  int error, ret;
  int i;
  _FC_SubSlot* subSlot;
  _FC_MeshSlot* meshSlot;
  _FC_DsSlot* dsSlot;
  int elemOffset;
  int numSide, numDistFactor, topoDim;
  int *elemIDs, *sideIDs;
  int numEntityPerElem, *elemToEntityConns;
  FC_ElementType elemType;

  // check input
  subSlot = _fc_getSubSlot(subset);
  if (!subSlot || !subSlot->header.name) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  if (subSlot->fileInfo.exoInfo.assoc_flag != 1) {
    fc_printfErrorMessage("Programmer error? can only lazy read side sets ("
                          "edges or faces)");
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Reading members for subset '%s'", subSlot->header.name);

  // setup
  meshSlot = _fc_getMeshSlot(subSlot->mesh);
  dsSlot = _fc_getDsSlot(meshSlot->ds);
  elemOffset = meshSlot->fileInfo.exoInfo.elemOffset;
  rc = fc_getMeshElementType(subSlot->mesh, &elemType);
  if (rc != FC_SUCCESS)
    return rc;
  topoDim = fc_getElementTypeTopoDim(elemType);
  if (topoDim == 2) {
    numEntityPerElem = fc_getElementTypeNumEdge(elemType);
    rc = fc_getMeshElementToEdgeConnsPtr(subSlot->mesh, &elemToEntityConns);
    if (rc != FC_SUCCESS)
      return rc;
  }
  else if (topoDim == 3) {
    numEntityPerElem = fc_getElementTypeNumFace(elemType);
    rc = fc_getMeshElementToFaceConnsPtr(subSlot->mesh, &elemToEntityConns);
    if (rc != FC_SUCCESS)
      return rc;
  }
  else {
    fc_printfErrorMessage("Programmer error? unexpected topo dim");
    return FC_ERROR;
  }

  // Read members
  error = ex_get_side_set_param(dsSlot->fileInfo.exoInfo.exoid, 
                                subSlot->fileInfo.exoInfo.nsssid, &numSide,
                                &numDistFactor);
  if (error != 0) {
    fc_printfErrorMessage("Exodus failed to get sides in side set");
    return FC_FILE_IO_ERROR;
  }
  elemIDs = (int*)malloc(numSide*sizeof(int));
  sideIDs = (int*)malloc(numSide*sizeof(int));
  if (!elemIDs || !sideIDs) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    free(elemIDs); free(sideIDs);
    return FC_MEMORY_ERROR;
  }
  error = ex_get_side_set(dsSlot->fileInfo.exoInfo.exoid, 
                          subSlot->fileInfo.exoInfo.nsssid, elemIDs, sideIDs);
  if (error != 0) {
    fc_printfErrorMessage("Exodus failed to get sides in side set");
    free(elemIDs); free(sideIDs);
    return FC_FILE_IO_ERROR;
  }

  // Process just the sides for this mesh & collect appropriate entity IDs
  for (i = 0; i < numSide; i++) {
    if (elemIDs[i] - 1 >= elemOffset && 
        elemIDs[i] - 1 <  elemOffset + meshSlot->numElement) {
      int elemID = elemIDs[i] - 1 - elemOffset;
      ret = fc_addIntToSortedIntArray(&subSlot->sia, 
               elemToEntityConns[elemID*numEntityPerElem+sideIDs[i]-1]);
      if (ret < 0) {
        fc_printfErrorMessage("Failed to add ID to sia");
        return ret;
      }
    }
  }
  free(elemIDs);
  free(sideIDs);
  if (subSlot->sia.numVal != subSlot->numMember) {
    fc_printfErrorMessage("Error collecting mesh ids from global ids");
    return FC_ERROR;
  }

  return FC_SUCCESS;
}

/**
 * \ingroup ExodusFileIO
 * \brief Read data from Exodus file into a basic or one step
 *        of a seq variable
 *
 * \modifications
 *    - 9/2/2005 WSD Created.
 *    - 9/9/2007 ACG turned into wrapper to split into basic
 *               or one step of seq var
 */
FC_ReturnCode _fc_readExodusVariableData(
  FC_Variable variable     /**< input - variable to be read into */
) {
  _FC_VarSlot* varSlot;

  // check input
  varSlot = _fc_getVarSlot(variable);
  if (!varSlot || !varSlot->header.name) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  if (varSlot->stepID < 0) {
    return _fc_readExodusAttributeData(variable);
  } else {
    return _fc_readExodusOneStepVariableData(variable);
  }
}


/**
 * \ingroup ExodusFileIO
 * \brief Read data from Exodus file into one step of a seq variable
 *
 * \modifications
 *    - 9/2/2005 WSD Created.
 *    - 9/9/2007 ACG renamed to this from _fc_readExodusVariableData
 */
FC_ReturnCode _fc_readExodusOneStepVariableData(
  FC_Variable variable     /**< input - variable to be read into */
) {
  int error;
  int i;
  _FC_VarSlot* varSlot;
  _FC_MeshSlot* meshSlot;
  _FC_DsSlot* dsSlot;
  double* global_data, *data;
  int exoid, exoStepID, *vertIDs, numTotalVert;


  // check input
  varSlot = _fc_getVarSlot(variable);
  if (!varSlot || !varSlot->header.name) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  if (varSlot->numComponent > 1) {
    fc_printfErrorMessage("lazy read of exodus data can't handle multicomp");
    return FC_ERROR;
  }
  if (varSlot->stepID < 0) {
    fc_printfErrorMessage("Programmer error? shouldn't call this on non "
			  "sequence variables");
    return FC_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Reading data for variable '%s'",
                      varSlot->header.name);

  // setup
  dsSlot = _fc_getDsSlot(varSlot->ds);
  meshSlot = _fc_getMeshSlot(varSlot->mesh);
  exoid = dsSlot->fileInfo.exoInfo.exoid;
  exoStepID = varSlot->stepID + 1;
  numTotalVert = dsSlot->fileInfo.exoInfo.numTotalVert;
  vertIDs = dsSlot->fileInfo.exoInfo.vertIDs[meshSlot->fileInfo.exoInfo.meshID];

  // Read the step's data
  // For nodal variables, have to read data for all nodes 
  if (!strcmp("n", varSlot->fileInfo.exoInfo.varType)) {
    data = (double*)malloc(varSlot->numDataPoint*sizeof(double));
    global_data = (double*)malloc(numTotalVert*sizeof(double));
    if (!data || !global_data) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      return FC_MEMORY_ERROR;
    }
    error = ex_get_nodal_var(exoid, exoStepID, varSlot->fileInfo.exoInfo.varIndex,
                             numTotalVert, global_data);
    if (error != 0) {
      fc_printfErrorMessage("Exodus failed to get nodal var data");
      return FC_FILE_IO_ERROR;
    }
    for (i = 0; i < meshSlot->numVertex; i++)
      data[i] = global_data[vertIDs[i]];
    free(global_data);
  }
  // for Element variables, just get this block's data
  else if (!strcmp("e", varSlot->fileInfo.exoInfo.varType)) {
    data = (double*)malloc(varSlot->numDataPoint*sizeof(double));
    error = ex_get_elem_var(exoid, exoStepID, varSlot->fileInfo.exoInfo.varIndex,
                            meshSlot->fileInfo.id,
                            meshSlot->numElement, data);
    if (error != 0) {
      fc_printfErrorMessage("Exodus failed to get element var data");
      return FC_FILE_IO_ERROR;
    }
  }
  else {
    fc_printfErrorMessage("Unknown exoVarType encountered");
    return FC_ERROR;
  }

  // all done
  varSlot->data = data;

  return FC_SUCCESS;
}


/**
 * \ingroup ExodusFileIO
 * \brief Read attribute data from Exodus file into a basic variable
 *
 * \modifications
 *    - 09/04/2007 ACG Commencing. 
 *    - 10/03/2007 adding nodal attr
 */
FC_ReturnCode _fc_readExodusAttributeData(
  FC_Variable variable     /**< input - variable to be read into */
) {

  int error;
  _FC_VarSlot* varSlot;
  _FC_MeshSlot* meshSlot;
  _FC_DsSlot* dsSlot;
  double* data;
  int exoid;

  // check input
  varSlot = _fc_getVarSlot(variable);
  if (!varSlot || !varSlot->header.name) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR)); 
    return FC_INPUT_ERROR;
  }
  if (varSlot->numComponent > 1) {
    fc_printfErrorMessage("lazy read of exodus data can't handle multicomp");
    return FC_ERROR;
  }
  if (varSlot->stepID >= 0) {
    fc_printfErrorMessage("Programmer error? shouldn't call this on "
                          "sequence variables");
    return FC_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Reading data for variable '%s'",
                      varSlot->header.name);


  //only works on elem or nodal attr for now

  if ( strcmp("e", varSlot->fileInfo.exoInfo.varType) &&
       strcmp("n", varSlot->fileInfo.exoInfo.varType) ){
    return FC_SUCCESS;
  }

  // setup
  dsSlot = _fc_getDsSlot(varSlot->ds);
  meshSlot = _fc_getMeshSlot(varSlot->mesh);
  exoid = dsSlot->fileInfo.exoInfo.exoid;

  // for Element variables, just get this block's data
  if (!strcmp("e", varSlot->fileInfo.exoInfo.varType)) {
    data = (double*)malloc(varSlot->numDataPoint*sizeof(double));
    if (!data){
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      return FC_MEMORY_ERROR;
    }
    error = ex_get_one_elem_attr(exoid, 
				 meshSlot->fileInfo.id,
				 varSlot->fileInfo.exoInfo.attrIndex,
				 data);
    if (error != 0) {
      fc_printfErrorMessage("Exodus failed to get element var data");
      return FC_FILE_IO_ERROR;
    }

  } else if (!strcmp("n", varSlot->fileInfo.exoInfo.varType)) {
    // For nodal attr have to read data for all nodes 
    double* global_data;
    int numTotalVert;
    int *vertIDs;
    int i;
    
    numTotalVert = dsSlot->fileInfo.exoInfo.numTotalVert;
    vertIDs = dsSlot->fileInfo.exoInfo.vertIDs[meshSlot->fileInfo.exoInfo.meshID];

    data = (double*)malloc(varSlot->numDataPoint*sizeof(double));
    global_data = (double*)malloc(numTotalVert*sizeof(double));
    if (!data || !global_data) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      return FC_MEMORY_ERROR;
    }

    error = ex_get_one_attr(exoid, EX_NODAL, 0, varSlot->fileInfo.exoInfo.attrIndex,
                             global_data);
    if (error != 0) {
      fc_printfErrorMessage("Exodus failed to get nodal var data");
      return FC_FILE_IO_ERROR;
    }

    for (i = 0; i < meshSlot->numVertex; i++)
      data[i] = global_data[vertIDs[i]];
    free(global_data);

  } else {
    fc_printfErrorMessage("Unknown exoVarType encountered");
    return FC_ERROR;
  }

  // all done
  varSlot->data = data;

  return FC_SUCCESS;

}


/**
 * \ingroup ExodusFileIO
 * \brief Write a dataset
 *
 * \description
 *
 * \todo 
 *    - Committed flag for seq vars and now non-seqvars. see what the
 *      issue is there.
 *    - why dont we combine node sets by names? currently
 *      I will not combine elem sets.
 *
 *
 * \modifications
 *    - 08/01/2005 WSD Created.
 *    - 07/31/2007 ACG added elemblk attr writeouts for elem non seq vars
 *    - 08/10/2007 ACG added nodal attr writeouts for vertex non seq vars
 *    - 10/05/2007 ACG write out elem sets
 *    - 10/17/2007 ACG removed names file
 */
FC_ReturnCode _fc_writeExodusDataset(
  FC_Dataset dataset, /**< input - the dataset to be written */
  char* filename      /**< input - the file name (including path) to write
                         into. */
) {
  FC_ReturnCode rc;
  int i, j, k, m, n, p;
  int maxID, *blk_ids;
  _FC_DsSlot* dsSlot;
  FC_Mesh* meshes;
  FC_Sequence* sequences, sequence;
  FC_Subset* nodeSets;
  FC_Subset* sideSets;
  FC_Subset* elemSets;
  int numMesh, numTotalVert, numTotalElem, maxNumDim, numNodeSet,
    numSideSet, numElemSet;
  int numSequence, numStep;
  int *numVerts, *numElems, *vertOffsets, *elemOffsets, *numDims, *numElemAttrs;
  int error, exoid, cpu_word_size, io_word_size;
  char* temp_name;
  char* coord_names[3] = { "coord_x", "coord_y", "coord_z" };
  double* coords[3] = { NULL, NULL, NULL };
  FC_SortedBlobArray sba;
  FC_SortedBlobArray nattr_sba;
  int numNodalAttrs;
  ex_init_params exoParams;

  // check input
  dsSlot = _fc_getDsSlot(dataset);
  if (dsSlot == NULL || dsSlot->writable == 0 || 
      (dsSlot->numSeq == 0 && dsSlot->numMesh == 0) ||
      !filename) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // log message
  fc_printfLogMessage("Writing dataset '%s' to Exodus file '%s'", 
                      dsSlot->header.name, filename);

  // --- Create the database file  ---

  // FIX??  should this check for an existing file?
  cpu_word_size = 8;
  io_word_size = 8;
  exoid = ex_create(filename, EX_CLOBBER, &cpu_word_size, &io_word_size);
  if (exoid < 0) {
    fc_printfErrorMessage("Failed to create Exodus dataset '%s'",
                          dsSlot->header.name);
    return FC_FILE_IO_ERROR;
  }

  // save file info
  dsSlot->baseFile = (char*) malloc((strlen(filename)+1)*sizeof(char));
  strcpy(dsSlot->baseFile, filename);
  dsSlot->fileInfo.exoInfo.cpu_word_size = cpu_word_size;
  dsSlot->fileInfo.exoInfo.io_word_size = io_word_size;
  dsSlot->fileInfo.exoInfo.exoid = exoid;

  // Preprocessing before can initialize file parameters (need num of things)
  rc = fc_getMeshes(dataset, &numMesh, &meshes);
  if (rc != FC_SUCCESS)
    return rc;
  numTotalVert = 0;
  numTotalElem = 0;
  maxNumDim = 0;
  numNodeSet = 0;
  numSideSet = 0;
  numElemSet = 0;
  vertOffsets = malloc(numMesh*sizeof(int));
  elemOffsets = malloc(numMesh*sizeof(int));
  numVerts = malloc(numMesh*sizeof(int));
  numElems = malloc(numMesh*sizeof(int));
  numDims = malloc(numMesh*sizeof(int));
  numElemAttrs = calloc(numMesh, sizeof(int)); // init to zero
  nodeSets = NULL;
  sideSets = NULL;
  elemSets = NULL;
  if (!vertOffsets || !elemOffsets || !numVerts || !numElems || !numDims ||
      !numElemAttrs){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    free(vertOffsets); free(elemOffsets); free(numVerts); free(numElems);
    free(numDims); free(numElemAttrs);
    return FC_MEMORY_ERROR;
  }
  for (i = 0; i < numMesh; i++) {
    int numVert, numElem, numDim;
    int numSubset, numVar;
    FC_Subset* subsets;
    FC_Variable* vars;
    FC_ElementType elemType;

    // get mesh info
    rc = fc_getMeshInfo(meshes[i], NULL, &numDim, &numVert, &numElem,
                        &elemType);
    if (rc != FC_SUCCESS)
      return rc;

    // count numVert & calc offsets and maxNumDim
    numVerts[i] = numVert;
    numElems[i] = numElem;
    vertOffsets[i] = numTotalVert;
    elemOffsets[i] = numTotalElem;
    numDims[i] = numDim;
    numTotalVert += numVert;
    numTotalElem += numElem;
    if (numDim > maxNumDim)
      maxNumDim = numDim;

    // count (& save) subsets for nodesets and sidesets
    rc = fc_getSubsets(meshes[i], &numSubset, &subsets);
    if (rc != FC_SUCCESS)
      return rc;
    for (j = 0; j < numSubset; j++) {
      int numMember;
      FC_AssociationType assoc;
      FC_Subset *temp_subsets;
      int assoc_flag;  // -1 = neither, 0 = nodeset, 1 = sideset, 2 = elemset

      // process
      rc = fc_getSubsetInfo(subsets[j], &numMember, NULL, &assoc);
      if (rc != FC_SUCCESS)
        return rc;
      if (!fc_isAssociationTypeValid(assoc)) {
        fc_printfErrorMessage("subset did not have valid assoc");
        return FC_ERROR;
      }
      assoc_flag = -1;
      switch(assoc) {
      case FC_AT_VERTEX: assoc_flag = 0; break;
      case FC_AT_EDGE: // if elements are 2D, add to sidesets
        if (fc_getElementTypeTopoDim(elemType) == 2) assoc_flag = 1; break;
      case FC_AT_FACE: // if elements are 3D, add to sidesets
        if (fc_getElementTypeTopoDim(elemType) == 3) assoc_flag = 1; break;
      case FC_AT_ELEMENT: 
	assoc_flag = 2; break;
      case FC_AT_WHOLE_MESH: 
      case FC_AT_WHOLE_DATASET: case FC_AT_UNKNOWN:
        ; // -1
      }
      switch (assoc_flag){
      case 0:
	// add to nodesets
        temp_subsets = realloc(nodeSets, (numNodeSet+1)*sizeof(FC_Subset));
        if (!temp_subsets) {
          fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
          return FC_MEMORY_ERROR;
        }
        nodeSets = temp_subsets;
        nodeSets[numNodeSet] = subsets[j];
        numNodeSet++;
	break;
      case 1:
        temp_subsets = realloc(sideSets, (numSideSet+1)*sizeof(FC_Subset));
        if (!temp_subsets) {
          fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
          return FC_MEMORY_ERROR;
        }
        sideSets = temp_subsets;
        sideSets[numSideSet] = subsets[j];
        numSideSet++;
	break;
      case 2:
	//elem set
	temp_subsets = realloc(elemSets, (numElemSet+1)*sizeof(FC_Subset));
        if (!temp_subsets) {
          fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
          return FC_MEMORY_ERROR;
        }
	elemSets = temp_subsets;
	elemSets[numElemSet] = subsets[j];
	numElemSet++;
	break;
      default:
        rc = fc_getSubsetName(subsets[j], &temp_name);
        if (rc != FC_SUCCESS)
          return rc;
        fc_printfWarningMessage("Skipping subset '%s' because exodus cannot "
                                "handle subsets of type %s on a %dD mesh", 
                                temp_name, fc_getAssociationTypeText(assoc),
                                fc_getElementTypeTopoDim(elemType));
        free(temp_name);
      }
    } // end of loop over subsets
    free(subsets);

    // count vars that can be elem attributes (not nodal)
    rc = fc_getVariables(meshes[i], &numVar, &vars);
    if (rc != FC_SUCCESS)
      return rc;
    for (j = 0; j < numVar; j++) {
      int numComp;
      FC_AssociationType assoc;
      FC_DataType dataType;
      rc = fc_getVariableInfo(vars[j], NULL, &numComp, &assoc, NULL,
                              &dataType);
      if (rc != FC_SUCCESS)
        return rc;
      if (assoc == FC_AT_ELEMENT &&
	  (dataType == FC_DT_INT || dataType == FC_DT_FLOAT ||
	   dataType == FC_DT_DOUBLE || dataType == FC_DT_CHAR)){
        numElemAttrs[i]+=numComp;
      }
     } // end of loop over vars
    free(vars);
  } // end of loop over meshes
  if (maxNumDim < 1) {
        fc_printfErrorMessage("Exodus cannot write files with "
                              "dimensionality < 1");
    return FC_ERROR;
  }

  // Initialize file parameters
  rc = fc_getDatasetName(dataset, &temp_name);
  strncpy(exoParams.title, temp_name, MAX_STR_LENGTH);
  exoParams.num_dim = maxNumDim;
  exoParams.num_nodes = numTotalVert;
  exoParams.num_edge = 0;
  exoParams.num_edge_blk = 0;
  exoParams.num_face = 0;
  exoParams.num_face_blk = 0;
  exoParams.num_elem = numTotalElem;
  exoParams.num_elem_blk = numMesh;
  exoParams.num_node_sets = numNodeSet;
  exoParams.num_edge_sets = 0;
  exoParams.num_face_sets = 0;
  exoParams.num_side_sets = numSideSet;
  exoParams.num_elem_sets = numElemSet;
  //we drop all maps
  exoParams.num_node_maps = 0;
  exoParams.num_edge_maps = 0;
  exoParams.num_face_maps = 0;
  exoParams.num_elem_maps = 0;

  error = ex_put_init_ext(exoid, &exoParams);

  free(temp_name);
  if (error != 0) {
    fc_printfErrorMessage("Exodus failed to put file parameters");
    return FC_FILE_IO_ERROR;
  }

  // Save Vert LUs
  dsSlot->fileInfo.exoInfo.numTotalVert = numTotalVert;
  dsSlot->fileInfo.exoInfo.numMesh = numMesh;
  dsSlot->fileInfo.exoInfo.vertIDs = (int**)malloc(numMesh*sizeof(int*));
  if (!dsSlot->fileInfo.exoInfo.vertIDs) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;  
  }
  for (i = 0; i < numMesh; i++) {
    int* vertIDs = (int*)malloc(numVerts[i]*sizeof(int));
    if (!vertIDs) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      return FC_MEMORY_ERROR;  
    }
    for (j = 0; j < numVerts[i]; j++)
      vertIDs[j] = vertOffsets[i] + j;
    dsSlot->fileInfo.exoInfo.vertIDs[i] = vertIDs;
  }

  // --- Write nodal coordinates ---

  // Make space
  for (i = 0; i < maxNumDim; i++) {
    coords[i] = (double*)malloc(numTotalVert*sizeof(double));
    if (!coords[i]) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      return FC_MEMORY_ERROR;
    }
  }
  // create global, components of coords
  for (i = 0; i < numMesh; i++) {
    double* temp_coords;
    rc = fc_getMeshCoordsPtr(meshes[i], &temp_coords);
    if (rc != FC_SUCCESS)
      return rc;
    for (j = 0; j < numVerts[i]; j++) {
      for (k = 0; k < numDims[i]; k++)
        coords[k][vertOffsets[i] + j] = temp_coords[j*numDims[i]+k];
      for (k = numDims[i]; k < maxNumDim; k++)
        coords[k][vertOffsets[i] + j] = 0.;
    } 
  }
  // do it
  error = ex_put_coord(exoid, coords[0], coords[1], coords[2]);
  for (i = 0; i < maxNumDim; i++)
    free(coords[i]);
  if (error != 0) {
    fc_printfErrorMessage("Exodus failed to put coords");
    return FC_FILE_IO_ERROR;
  }
  error = ex_put_coord_names(exoid, coord_names);
  if (error != 0) {
    fc_printfErrorMessage("Exodus failed to put coord names");
    return FC_FILE_IO_ERROR;
  }

  // --- Write element blocks (i.e. the conns) ---

  // Create block ids if they do not already exits
  // To get unique ids, start w/ the max existing id + 1
  maxID = 0;
  for (i = 0; i < numMesh; i++) {
    _FC_MeshSlot* meshSlot = _fc_getMeshSlot(meshes[i]);
    if (meshSlot->fileInfo.id > maxID)
      maxID = meshSlot->fileInfo.id;
  }
  for (i = 0; i < numMesh; i++) {
    _FC_MeshSlot* meshSlot = _fc_getMeshSlot(meshes[i]);
    if (meshSlot->fileInfo.id == -1) {
      maxID++;
      meshSlot->fileInfo.id = maxID;
    }
  }
  blk_ids = (int*)malloc(numMesh*sizeof(numMesh));
  if (!blk_ids) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  for (i = 0; i < numMesh; i++) {
    _FC_MeshSlot* meshSlot = _fc_getMeshSlot(meshes[i]);
    blk_ids[i] = meshSlot->fileInfo.id;
  }


  for (i = 0; i < numMesh; i++) {
    _FC_MeshSlot* meshSlot = _fc_getMeshSlot(meshes[i]);
    char exoElemType[MAX_STR_LENGTH+1];
    int numVertPerElem, *temp_conns, *conns;
    meshSlot->fileInfo.exoInfo.meshID = i;
    meshSlot->fileInfo.exoInfo.elemOffset = elemOffsets[i];
    switch (meshSlot->elemType) {
    case FC_ET_POINT:
      if (maxNumDim == 1)
        sprintf(exoElemType, "POINT"); // FIX!?
      else if (maxNumDim == 2)
        sprintf(exoElemType, "CIRCLE");
      else if (maxNumDim == 3)
        sprintf(exoElemType, "SPHERE");
      break;
    case FC_ET_LINE:  sprintf(exoElemType, "BEAM");   break;
    case FC_ET_TRI:   sprintf(exoElemType, "TRI");    break;
    case FC_ET_QUAD:  sprintf(exoElemType, "QUAD");   break;
    case FC_ET_TET:   sprintf(exoElemType, "TETRA");  break;
    case FC_ET_PYRAMID: sprintf(exoElemType, "PYRAMID"); break;
    case FC_ET_PRISM: sprintf(exoElemType, "WEDGE");  break;
    case FC_ET_HEX:   sprintf(exoElemType, "HEX");    break;
    default:
      fc_printfErrorMessage("Elementype %s not supported",
                            fc_getElementTypeText(meshSlot->elemType));
      return FC_ERROR;
    }
    numVertPerElem = fc_getElementTypeNumVertex(meshSlot->elemType);
    error = ex_put_elem_block(exoid, meshSlot->fileInfo.id,
                              exoElemType, numElems[i],
                              numVertPerElem, numElemAttrs[i]);
    if (error != 0) {
      fc_printfErrorMessage("Exodus failed to put elem blk meta data");
      return FC_FILE_IO_ERROR;
    }
 
    rc = fc_getMeshElementConnsPtr(meshes[i], &temp_conns);
    if (rc != FC_SUCCESS)
      return rc;
    conns = malloc(numElems[i]*numVertPerElem*sizeof(int));
    if (!conns) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      return FC_MEMORY_ERROR;
    }
    for (j = 0; j < numElems[i]*numVertPerElem; j++) 
      conns[j] = temp_conns[j] + vertOffsets[i] + 1;
    error = ex_put_elem_conn(exoid, meshSlot->fileInfo.id, conns);
    free(conns);
    if (error != 0) {
      fc_printfErrorMessage("Exodus failed to put elem conns");
      return FC_FILE_IO_ERROR;
    }
    // restrict name to MAX_STR_LENGTH
    temp_name = malloc(MAX_STR_LENGTH*sizeof(char));
    strncpy(temp_name, meshSlot->header.name, MAX_STR_LENGTH);
    temp_name[MAX_STR_LENGTH-1] = '\0';
    if (strlen(meshSlot->header.name) >= MAX_STR_LENGTH)
      fc_printfWarningMessage("Mesh '%s's name is being clipped by Exodus "
                              "to '%s'", meshSlot->header.name, temp_name);
    error = ex_put_name(exoid, EX_ELEM_BLOCK, meshSlot->fileInfo.id,
                        temp_name);
    free(temp_name);
    if (error != 0) {
      fc_printfErrorMessage("Exodus failed to put elem block name");
      return FC_FILE_IO_ERROR;
    }
    meshSlot->header.committed = 1;
  }

  // --- Write side sets ---

  // FIX? change to use concat method which is faster

  for (i = 0; i < numSideSet; i++) {
    int meshID, sideset_id;
    int numMember, *members, *sides;
    int numEntityPerElem, *elemToEntityConns;
    int *numElemPerEntity, **elemParentsPerEntity;
    _FC_SubSlot* subSlot;
    _FC_MeshSlot* meshSlot;

    // setup
    subSlot = _fc_getSubSlot(sideSets[i]);
    meshSlot = _fc_getMeshSlot(subSlot->mesh);
    meshID = meshSlot->fileInfo.exoInfo.meshID;
    if (subSlot->assoc == FC_AT_EDGE) {
      numEntityPerElem = fc_getElementTypeNumEdge(meshSlot->elemType);
      rc = fc_getMeshElementToEdgeConnsPtr(subSlot->mesh, &elemToEntityConns);
      if (rc != FC_SUCCESS)
        return rc;
      rc = _fc_getMeshElementParentsOfEdgesPtr(subSlot->mesh, &numElemPerEntity,
                                            &elemParentsPerEntity);
      if (rc != FC_SUCCESS)
        return rc;
    }
    else if (subSlot->assoc == FC_AT_FACE) {
      numEntityPerElem = fc_getElementTypeNumFace(meshSlot->elemType);
      rc = fc_getMeshElementToFaceConnsPtr(subSlot->mesh, &elemToEntityConns);
      if (rc != FC_SUCCESS)
        return rc;
      rc = _fc_getMeshElementParentsOfFacesPtr(subSlot->mesh, &numElemPerEntity,
                                            &elemParentsPerEntity);
      if (rc != FC_SUCCESS)
        return rc;
    }
    else {
      fc_printfErrorMessage("Programmer error? should be edge or face subset");
      return rc;
    }
   
    // Get subset data - this is by copy
    rc = fc_getSubsetMembersAsArray(sideSets[i], &numMember, &members);
    if (rc != FC_SUCCESS)
      return rc;
    sideset_id = i+1;
    // Put side set params
    error = ex_put_side_set_param(exoid, sideset_id, numMember, 0);
    if (error != 0) {
      fc_printfErrorMessage("Exodus failed to put side set param");
      return FC_FILE_IO_ERROR;
    }
    // have to convert mesh side ids to global elem ids & local to elem side id
    // (members gets turned into element IDs)
    sides = (int*)malloc(numMember*sizeof(int));
    if (!sides) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      return FC_MEMORY_ERROR;
    }
    for (j = 0; j < numMember; j++) {
      int elemID, entityID;
      // get first element that has this entity, then figure out local ID
      elemID = elemParentsPerEntity[members[j]][0];
      for (k = 0; k < numEntityPerElem; k++) {
        if (elemToEntityConns[elemID*numEntityPerElem+k] == members[j]) {
          entityID = k;
          break;
        }
      }
      members[j] = elemID + meshSlot->fileInfo.exoInfo.elemOffset + 1;
      sides[j] = entityID + 1;
    }
    // put side set
    error = ex_put_side_set(exoid, sideset_id, members, sides);
    free(members);
    free(sides);
    if (error != 0) {
      fc_printfErrorMessage("Exodus failed to put side set");
      return FC_FILE_IO_ERROR;
    }
    // restrict name to MAX_STR_LENGTH
    temp_name = malloc(MAX_STR_LENGTH*sizeof(char));
    strncpy(temp_name, subSlot->header.name, MAX_STR_LENGTH);
    temp_name[MAX_STR_LENGTH-1] = '\0';
    if (strlen(subSlot->header.name) >= MAX_STR_LENGTH)
      fc_printfWarningMessage("Subset '%s's name is being clipped by Exodus"
                              "to '%s'", subSlot->header.name, temp_name);
    error = ex_put_name(exoid, EX_SIDE_SET, sideset_id, temp_name);
    free(temp_name);
    if (error != 0) {
      fc_printfErrorMessage("Exodus failed to put side set name");
      return FC_FILE_IO_ERROR;
    }
    subSlot->fileInfo.exoInfo.assoc_flag = 0;
    subSlot->fileInfo.exoInfo.nsssid = sideset_id;
    subSlot->header.committed = 1;
  }
  free(sideSets);

  // --- Write node sets ---

  // FIX? change to use concat method which is supposedly faster

  for (i = 0; i < numNodeSet; i++) {
    int meshID, nodeset_id;
    int numMember, *members;
    _FC_SubSlot* subSlot;
    _FC_MeshSlot* meshSlot;

    // setup
    subSlot = _fc_getSubSlot(nodeSets[i]);
    meshSlot = _fc_getMeshSlot(subSlot->mesh);
    meshID = meshSlot->fileInfo.exoInfo.meshID;
    // Get subset data - this is by copy
    rc = fc_getSubsetMembersAsArray(nodeSets[i], &numMember, &members);
    if (rc != FC_SUCCESS)
      return rc;
    nodeset_id = i+1;
    // Put node set params
    error = ex_put_node_set_param(exoid, nodeset_id, numMember, 0);
    if (error != 0) {
      fc_printfErrorMessage("Exodus failed to put node set param");
      return FC_FILE_IO_ERROR;
    }
    // have to convert mesh node ids to global node ids
    for (j = 0; j < numMember; j++)
      members[j] += vertOffsets[meshID] + 1;
    // put node set
    error = ex_put_node_set(exoid, nodeset_id, members);
    free(members);
    if (error != 0) {
      fc_printfErrorMessage("Exodus failed to put node set");
      return FC_FILE_IO_ERROR;
    }
    // restrict name to MAX_STR_LENGTH
    temp_name = malloc(MAX_STR_LENGTH*sizeof(char));
    strncpy(temp_name, subSlot->header.name, MAX_STR_LENGTH);
    temp_name[MAX_STR_LENGTH-1] = '\0';
    if (strlen(subSlot->header.name) >= MAX_STR_LENGTH)
      fc_printfWarningMessage("Subset '%s's name is being clipped by Exodus"
                              "to '%s'", subSlot->header.name, temp_name);
    error = ex_put_name(exoid, EX_NODE_SET, nodeset_id, temp_name);
    free(temp_name);
    if (error != 0) {
      fc_printfErrorMessage("Exodus failed to put node set name");
      return FC_FILE_IO_ERROR;
    }
    subSlot->fileInfo.exoInfo.assoc_flag = 0;
    subSlot->fileInfo.exoInfo.nsssid = nodeset_id;
    subSlot->header.committed = 1;
  }
  free(nodeSets);

  // write elem sets
  //  printf("num elem set + %d\n", numElemSet);
  for (i = 0; i < numElemSet; i++) {
    int meshID, elemset_id;
    int numMember, *members;
    _FC_SubSlot* subSlot;
    _FC_MeshSlot* meshSlot;

    // setup
    subSlot = _fc_getSubSlot(elemSets[i]);
    meshSlot = _fc_getMeshSlot(subSlot->mesh);
    meshID = meshSlot->fileInfo.exoInfo.meshID;
    // Get subset data - this is by copy
    rc = fc_getSubsetMembersAsArray(elemSets[i], &numMember, &members);
    if (rc != FC_SUCCESS)
      return rc;
    elemset_id = i+1;
    // Put elem set params
    error = ex_put_set_param(exoid, EX_ELEM_SET, elemset_id, numMember, 0);
    if (error != 0) {
      fc_printfErrorMessage("Exodus failed to put elem set param");
      return FC_FILE_IO_ERROR;
    }
      
    // have to convert mesh elem ids to global elem ids
    for (j = 0; j < numMember; j++){
      members[j] += meshSlot->fileInfo.exoInfo.elemOffset + 1;
    }
    // put elem set
    error = ex_put_set(exoid, EX_ELEM_SET, elemset_id, members, NULL);
    free(members);
    if (error != 0) {
      fc_printfErrorMessage("Exodus failed to put elem set");
      return FC_FILE_IO_ERROR;
    }
    // restrict name to MAX_STR_LENGTH
    temp_name = malloc(MAX_STR_LENGTH*sizeof(char));
    strncpy(temp_name, subSlot->header.name, MAX_STR_LENGTH);
    temp_name[MAX_STR_LENGTH-1] = '\0';
    if (strlen(subSlot->header.name) >= MAX_STR_LENGTH)
      fc_printfWarningMessage("Subset '%s's name is being clipped by Exodus"
                              "to '%s'", subSlot->header.name, temp_name);
    //    printf("putting elem set <%s>\n", temp_name);
    error = ex_put_name(exoid, EX_ELEM_SET, elemset_id, temp_name);
    free(temp_name);
    if (error != 0) {
      fc_printfErrorMessage("Exodus failed to put elem set name");
      return FC_FILE_IO_ERROR;
    }
    subSlot->fileInfo.exoInfo.assoc_flag = 2;
    subSlot->fileInfo.exoInfo.nsssid = elemset_id;
    subSlot->header.committed = 1;
  }
  free(elemSets);



  // --- Write vars

  // FC_AT_WHOLE_MESH, FC_DT_INT vars are properties
  //   (other FC_AT_WHOLE_MESH are discarded)
  // FC_AT_ELEMENT are attributes
  // FC_AT_VERTEX are attributes and will be handled later
  for (i = 0; i < numMesh; i++) {
    int numVar;
    FC_Variable *vars;
    char** elemAttrNames = NULL;
    int currElemAttr = 0;
    char* temp_name2;

    rc = fc_getVariables(meshes[i], &numVar, &vars);
    if (rc != FC_SUCCESS)
      return rc;
    if (numElemAttrs[i] > 0){
      elemAttrNames = (char**)malloc(numElemAttrs[i]*sizeof(char*));
      if (!elemAttrNames){
	return FC_MEMORY_ERROR;
      }
    }

    for (j = 0; j < numVar; j++) {
      FC_AssociationType assoc;
      FC_DataType dataType;
      FC_MathType mathType;
      int numComp, *prop_vals;
      rc = fc_getVariableInfo(vars[j], NULL, &numComp, &assoc, &mathType,
                              &dataType);
      if (rc != FC_SUCCESS)
        return rc;
      rc = fc_getVariableName(vars[j], &temp_name);
      if (rc != FC_SUCCESS)
	return rc;
      switch (assoc){
      case FC_AT_WHOLE_MESH:
	if (dataType != FC_DT_INT) {
          fc_printfWarningMessage("Variable '%s' on mesh %d discard because "
                                  "cannot write datatype of %s", temp_name,
                                  i, fc_getDataTypeText(dataType));
          free(temp_name);
	  break;
	}
        // 'ID' is reserved for Exodus
        if (!strcmp(temp_name, "ID")) {
          // Purposefully not adding warning message because we'd get one every
          // time we rewrote an Exodus dataset -WSD 11/24/2006
          free(temp_name);
	  break;
        }
        rc = fc_getVariableDataPtr(vars[j], (void**)&prop_vals);
        if (rc != FC_SUCCESS)
          return rc;
        for (k = 0; k < numComp; k++) {
          rc = _fc_createExoVarCompName(temp_name, numComp, mathType, k,
                                        &temp_name2);
          if (rc != FC_SUCCESS)
            return rc;
          error = ex_put_prop(exoid, EX_ELEM_BLOCK, blk_ids[i], temp_name2,
                              prop_vals[k]);
          if (error != 0) {
            fc_printfErrorMessage("Failed to put prop '%s' on mesh #%d",
                                  temp_name2, i);
            return FC_FILE_IO_ERROR;
          }
          free(temp_name2);
	}
        free(temp_name);
 	// end of if FC_AT_WHOLE_MESH
	break;
      case FC_AT_ELEMENT:
	{
	  void* elemAttrVals = NULL;
	  double* elemAttrData = NULL;

	  //note that numElemAttrs has taken into account dataType
	  //and numComponent above, since we needed the numattr
	  //to establish the elem block. cross ref with that if
	  //this code ever changes

	  if (dataType != FC_DT_DOUBLE  &&
	      dataType != FC_DT_FLOAT  &&
	      dataType != FC_DT_INT &&
	      dataType != FC_DT_CHAR){
	    fc_printfWarningMessage("Variable '%s' on mesh %d elems discard "
				    "because cannot write datatype of %s",
				    temp_name,
				    i, fc_getDataTypeText(dataType));
	    free(temp_name);
	    break;
	  }
	  // 'ID' reserved for Exodus 
	  if (!strcmp(temp_name, "ID")) {
	    //this actually shouldnt matter because ID will be a prop for this var
	    free(temp_name);
	    break;
	  }
	  elemAttrData = (double*)malloc(numElems[i]*sizeof(double));
	  if (!elemAttrData){
	    return FC_MEMORY_ERROR;
	  }
	  rc = fc_getVariableDataPtr(vars[j], &elemAttrVals);
	  if (rc != FC_SUCCESS){
	    free(temp_name);
	    return rc;
	  }

	  for (k = 0; k < numComp; k++) {
	    int kk;
	    for (kk = 0; kk < numElems[i]; kk++){
        	switch(dataType){
		case FC_DT_INT:     
		  elemAttrData[kk] = ((int*)elemAttrVals)[kk*numComp+k];
		  break;
		case FC_DT_FLOAT:   
		  elemAttrData[kk] = ((float*)elemAttrVals)[kk*numComp+k];
		  break;
		case FC_DT_DOUBLE:  
		  elemAttrData[kk] = ((double*)elemAttrVals)[kk*numComp+k];
		  break;
		case FC_DT_CHAR:  
		  elemAttrData[kk] = ((char*)elemAttrVals)[kk*numComp+k];
		  break;
		default:
		  //shouldnt happen
		  break;
		}
	    }
	    //attribute numbering starts at 1
	    error = ex_put_one_elem_attr(exoid, blk_ids[i],currElemAttr+1,
					 elemAttrData);
	    if (error != 0){
	      free(temp_name);
	      fc_printfErrorMessage("Failed to put one elem attr mesh #%d",
				    i);
	      return FC_FILE_IO_ERROR;
	    }
	    rc = _fc_createExoVarCompName(temp_name, numComp, mathType, k,
					  &temp_name2);
	    if (rc != FC_SUCCESS){
	      free(temp_name);
	      return rc;
	    }
	    elemAttrNames[currElemAttr] = temp_name2;
	    currElemAttr++;
	  } //numComp
	  free(temp_name);
	  free(elemAttrData);
	}
	break;
      case FC_AT_VERTEX:
	//these will be handled below
        free(temp_name);
	break;
      default:
        fc_printfWarningMessage("Variable '%s' on mesh %d discarded because "
                                "cannot write assoc type of %s", temp_name,
                                i, fc_getAssociationTypeText(assoc));
        free(temp_name);
	break;
      } // end of switch
    } // end of loop over vars
    free(vars);

    if (numElemAttrs[i] > 0){
      error = ex_put_elem_attr_names (exoid, blk_ids[i], elemAttrNames);
      for (j = 0; j < numElemAttrs[i]; j++){
	free(elemAttrNames[j]);
      }
      free(elemAttrNames);
      if (error != 0){
	fc_printfErrorMessage("Failed to put elem attr names on mesh #%d",
			      i);
	return FC_FILE_IO_ERROR;
      }
    }

  } //end of loop over meshes

  //we dont need the attributes anymore
  if (numElemAttrs) free(numElemAttrs);

  // --- Write sequence (there can only be 0 or 1) ---
  
  rc = fc_getSequences(dataset, &numSequence, &sequences);
  if (rc != FC_SUCCESS)
    return rc;
  if (numSequence > 0) {
    FC_DataType dataType;
    void* data;

    // Warn if too many sequences
    if (numSequence > 1) {
      rc = fc_getSequenceName(sequences[i], &temp_name);
      if (rc != FC_SUCCESS)
        return rc;
      fc_printfWarningMessage("Exodus can only have 1 sequence, using the "
                              "first one: '%s'", temp_name);
      free(temp_name);
      for (i = 1; i < numSequence; i++) {
        rc = fc_getSequenceName(sequences[i], &temp_name);
        if (rc != FC_SUCCESS)
          return rc;
        fc_printfWarningMessage("Sequence '%s' discarded (also all seq vars "
                                " on this seq discarded)", temp_name);
        free(temp_name);
      } 
    }
    
    // do it
    sequence = sequences[0]; // keep for later
    rc = fc_getSequenceInfo(sequence, &numStep, &dataType);
    if (rc != FC_SUCCESS)
      return rc;
    rc = fc_getSequenceCoordsPtr(sequence, &data);
    if (rc != FC_SUCCESS)
      return rc;
    for (i = 0; i < numStep; i++) {
      double time_coord;
      switch (dataType) {
      case FC_DT_CHAR:     time_coord = ((char*)data)[i];    break;
      case FC_DT_INT:      time_coord = ((int*)data)[i];     break;
      case FC_DT_FLOAT:    time_coord = ((float*)data)[i];   break;
      case FC_DT_DOUBLE:   time_coord = ((double*)data)[i];  break;
      case FC_DT_UNKNOWN:
        fc_printfErrorMessage("sequence has unknown datatype");
        return FC_ERROR;
      }
      error = ex_put_time(exoid, i+1, &time_coord);
      if (error != 0) {
        fc_printfErrorMessage("Failed to put time value %d (of %d)", i, numStep);
        return FC_FILE_IO_ERROR;
      }
    } // end of loop over steps
    // Sequence is not read lazily, so don't set committed flag
  }
  free(sequences);

  // --- Write global seq variables ---

  {
    int numGlbSeqVar, *numStepPerSeqVar;
    int numGlbVar;
    int* numComps = NULL; // a zero entry means was on another seq
    FC_Variable **glbSeqVars;
    
    // get the vars
    rc = fc_getGlobalSeqVariables(dataset, &numGlbSeqVar,
                                  &numStepPerSeqVar, &glbSeqVars);
    if (rc != FC_SUCCESS)
      return rc;
    
    // Get number of globar vars (split up multi-component vars,
    // ignore vars on other sequences)
    numGlbVar = 0;
    if (numGlbSeqVar > 0) {
      numComps = malloc(numGlbSeqVar*sizeof(int));
      if (!numComps) {
	fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
	return FC_MEMORY_ERROR;
      }
      for (i = 0; i < numGlbSeqVar; i++) {
	FC_Sequence temp_sequence;
	rc = fc_getSequenceFromSeqVariable(numStepPerSeqVar[i], glbSeqVars[i],
					   &temp_sequence);
	if (rc != FC_SUCCESS)
	  return rc;
	if (!FC_HANDLE_EQUIV(sequence, temp_sequence)) {
	  numComps[i] = 0;
	  fc_getVariableName(glbSeqVars[i][0], &temp_name);
	  fc_printfWarningMessage("Global var '%s is not on the first "
				  "sequence, skipping", temp_name);
	  free(temp_name);
	  continue;
	}
	rc = fc_getVariableNumComponent(glbSeqVars[i][0], &numComps[i]);
	if (rc != FC_SUCCESS)
	  return rc;
	numGlbVar += numComps[i];
      }
    }

    // Further processing only necessary if there are any
    if (numGlbVar > 0) {
      char** names;
     	
      // write number
      error = ex_put_var_param(exoid, "g", numGlbVar);
      if (error != 0) {
	fc_printfErrorMessage("Failed to put 'g' var param");
	return FC_FILE_IO_ERROR;
      }

      // compile names
      names = (char**)malloc(numGlbVar*sizeof(char*));
      if (!names) {
	fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
	return FC_MEMORY_ERROR;
      }
      numGlbVar = 0;
      for (i = 0; i < numGlbSeqVar; i++) {
	FC_MathType mathtype;
	if (numComps[i] == 0)
	  continue; // skip if on another sequence
	rc = fc_getVariableName(glbSeqVars[i][0], &temp_name);
	if (rc != FC_SUCCESS)
	  return rc;
	rc = fc_getVariableMathType(glbSeqVars[i][0], &mathtype);
	if (rc != FC_SUCCESS)
	  return rc;
        for (j = 0; j < numComps[i]; j++) {
          rc = _fc_createExoVarCompName(temp_name, numComps[i], mathtype, j,
                                        &names[numGlbVar]);
          if (rc != FC_SUCCESS)
            return rc;
          numGlbVar++;
        }
        free(temp_name);
      }
      
      // write number and names
      error = ex_put_var_names(exoid, "g", numGlbVar, names);
      if (error != 0) {
	fc_printfErrorMessage("Failed to put 'g' var names");
	return FC_FILE_IO_ERROR;
      }  
            
      // cleanup names
      for (i = 0; i < numGlbVar; i++)
	free(names[i]);
      free(names);
      
      // write data for each step
      FC_DataType dataType;
      void* data;
      double values[numGlbVar];
      for (j = 0; j < numStep; j++) {
	int current = 0;
	for (i = 0; i < numGlbSeqVar; i++) {
	  if (numComps[i] == 0) // skip if on another sequence
	    continue;
	  rc = fc_getVariableDataType(glbSeqVars[i][0], &dataType);
	  if (rc != FC_SUCCESS)
	    return rc;
	  rc = fc_getVariableDataPtr(glbSeqVars[i][j], &data);
	  if (rc != FC_SUCCESS)
	    return rc;
	  for (k = 0; k < numComps[i]; k++) {
	    if (dataType == FC_DT_CHAR)
	      values[current] = ((char*)data)[k];
	    else if (dataType == FC_DT_INT)
	      values[current] = ((int*)data)[k];
	    else if (dataType == FC_DT_FLOAT)
              values[current] = ((float*)data)[k];
	    else if (dataType == FC_DT_DOUBLE)
	      values[current] = ((double*)data)[k];
	    current++;
	  }
	}
	error = ex_put_glob_vars(exoid, j+1, numGlbVar, values);
	if (error != 0) {
	  fc_printfErrorMessage("Exodus failed to put global vars for "
				"step %d", j);
	  return FC_FILE_IO_ERROR;
	}
	// global vars not read lazily, so no committed flag
      }
    }

    // cleanup (outside because numGlbVar could be zero even if numGlbSeqVar not)
    free(numComps);
    for (i = 0; i < numGlbSeqVar; i++)
      free(glbSeqVars[i]);
    free(glbSeqVars);
    free(numStepPerSeqVar);
  }


  // Write nodal vars
  // count vars that can be nodal attributes 
  // FIXME: may want to check for ID, though that shouldnt matter
  // because it will be a prop 
  fc_initSortedBlobArray(&nattr_sba);
  for (i = 0; i < numMesh; i++) {
    int numVar;
    FC_Variable *vars;
    rc = fc_getVariables(meshes[i], &numVar, &vars);
    if (rc != FC_SUCCESS)
      return rc;
    for (j = 0; j < numVar; j++) {
      int numComp;
      FC_AssociationType assoc;
      FC_DataType dataType;
      rc = fc_getVariableInfo(vars[j], NULL, &numComp, &assoc, NULL,
                              &dataType);
      if (rc != FC_SUCCESS)
        return rc;
      if (assoc == FC_AT_VERTEX &&
	  (dataType == FC_DT_INT || dataType == FC_DT_FLOAT ||
	   dataType == FC_DT_DOUBLE || dataType == FC_DT_CHAR) ){	   
	FC_Variable *tempVar = (FC_Variable*) malloc(sizeof(FC_Variable));
	if (!tempVar){
	  fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
	  return FC_MEMORY_ERROR;
	}

	*tempVar = vars[j];
	rc = fc_addBlobToSortedBlobArray(&nattr_sba, tempVar,
					 _fc_exo_seqVar_meshID_cmp); 
	if (rc < FC_SUCCESS) {
	  fc_printfErrorMessage("failed to add nodal var to nattr_sba");
	  return FC_ERROR;
	}
	// Warn about "duplicates" which are dropped ...
	if (rc == 0) {
	  fc_printfWarningMessage("Encountered duplicate var, ignoring");
	  free(tempVar);
	}
      }
    } // end of loop over vars
    free(vars);
  } // end of loop over meshes

  numNodalAttrs = 0;
  if (nattr_sba.numBlob > 0) {
    FC_Variable ***sortedVars; //[varIndex][meshID] 
    FC_Variable *last_Var, *curr_Var; 
    int* numNodalComps;

    numNodalAttrs = 1;
    if (nattr_sba.numBlob > 1){
      last_Var = (FC_Variable*)nattr_sba.blobs[0];
      for ( i = 1; i < nattr_sba.numBlob; i++){
	//we already know the assoc types are all the same
	curr_Var = (FC_Variable*)nattr_sba.blobs[i];
	if (_fc_exo_seqVar_cmp(curr_Var, last_Var) != 0){
	  numNodalAttrs++; //note that this hasnt split things up into components
	}
	last_Var = curr_Var;
      }
    }

    // convert to array of var handles (w/ possible empty slots)
    // handing off responsibility for vars
    k = 0;
    curr_Var = (FC_Variable*)nattr_sba.blobs[k];
    sortedVars = (FC_Variable***)malloc(numNodalAttrs*sizeof(FC_Variable**));
    if (!sortedVars){
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      return FC_MEMORY_ERROR;
    }
    for (j = 0; j < numNodalAttrs; j++){
      sortedVars[j] = (FC_Variable**)calloc(numMesh,sizeof(FC_Variable*));
      if (!sortedVars[j]) {
	fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
	return FC_MEMORY_ERROR;
      }

      last_Var = curr_Var; // quaranteed to get 1st one
      while (k < nattr_sba.numBlob) {
	curr_Var = (FC_Variable*)nattr_sba.blobs[k];
	int meshID = _fc_getMeshSlot(_fc_getVarSlot(*curr_Var)->mesh)->fileInfo.exoInfo.meshID;
	if (_fc_exo_seqVar_cmp(curr_Var, last_Var) != 0)
	  break; // move onto next var
	sortedVars[j][meshID] = curr_Var;
	k++;
      }
    }

    // Write number and names of nodal variables
    {
      int current;
      int numVar = numNodalAttrs;
      int numExoVar = numVar; // currently includes multicomponents as solos                    
      char** names, **temp_names;

      // assemble names array                                                                
      // CHECK for multicomponent => have to write as separate components                    
      names = (char**)malloc(numVar*sizeof(char*)); // this might grow                       
      numNodalComps = (int*)malloc(numVar*sizeof(int)); // won't grow                          
      if (!names || !numNodalComps) {
        free(names); free(numNodalComps);
        fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
        return FC_MEMORY_ERROR;
      }
      current = 0;

      for (j = 0; j < numVar; j++) {
        int validID;
        int numComp;
        FC_MathType mathtype;

        // Find a valid seqVar (not all meshes may have a seqvar)                            
        validID = -1;
        for (k = 0; k < numMesh; k++) {
          if (sortedVars[j][k]) {
            validID = k;
            break;
          }
        }
        if (validID < 0) {
          fc_printfErrorMessage("Programmer error?");
          return FC_ERROR;
        }
	// Add the var's name (or names if multicomponent)                                   
        rc = fc_getVariableName(*(sortedVars[j][validID]),
                                &temp_name);
        if (rc != FC_SUCCESS)
          return rc;
        rc = fc_getVariableInfo(*(sortedVars[j][validID]), NULL,
                                &numComp, NULL, &mathtype, NULL);
        if (rc != FC_SUCCESS)
          return rc;
        numNodalComps[j] = numComp; // save for later usage                                   
        if (numComp > 1) {
          numExoVar += numComp - 1;  // -1 'cause 1 component already counted                
          temp_names = realloc(names, numExoVar*sizeof(char*));
          if (!temp_names) {
            fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
            return FC_MEMORY_ERROR;
          }
          names = temp_names;
        }
        for (k = 0; k < numComp; k++) {
          rc = _fc_createExoVarCompName(temp_name, numComp, mathtype, k,
                                        &names[current]);
          if (rc != FC_SUCCESS)
            return rc;
          current++;
        }
        free(temp_name);
      }

      error = ex_put_attr_param(exoid, EX_NODAL, 0, numExoVar);
      if (error != 0) {
	fc_printfErrorMessage("Failed to put nodal var param");
	return FC_FILE_IO_ERROR;
      }
      error = ex_put_attr_names(exoid, EX_NODAL, 0, names);
      if (error != 0) {
	fc_printfErrorMessage("Failed to put nodal var names");
	return FC_FILE_IO_ERROR;
      }

      // cleanup names array 
      for (j = 0; j < numExoVar; j++)
        free(names[j]);
      free(names);
    }

    // write nodal data
    {
      // write var for all nodes, per step                                                   
      double* values = malloc(numTotalVert*sizeof(double));
      if (!values) {
	fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
	return FC_MEMORY_ERROR;
      }
      void* data;
      FC_DataType dataType;
      int current = 0;
      for (k = 0; k < numNodalAttrs; k++){
	for (m = 0; m < numNodalComps[k]; m++) {
	  for (n = 0; n < numTotalVert; n++)
	    values[n] = 0; // initialize in case values don't exist                      
	  for (n = 0; n < numMesh; n++) {
	    // FIX add a warning for this - do it in names section above ??
	    if (!sortedVars[k][n])
	      continue; // skip if var doesn't exist == pad with zeros                   
	    fc_getVariableDataType(*sortedVars[k][n], &dataType);
	    fc_getVariableDataPtr(*sortedVars[k][n], &data);
	    for (p = 0; p < numVerts[n]; p++) {
	      if (dataType == FC_DT_INT)
		values[vertOffsets[n] + p] = ((int*)data)[p*numNodalComps[k]+m];
	      else if (dataType == FC_DT_FLOAT)
		values[vertOffsets[n] + p] = ((float*)data)[p*numNodalComps[k]+m];
	      else if (dataType == FC_DT_DOUBLE)
		values[vertOffsets[n] + p] = ((double*)data)[p*numNodalComps[k]+m];
	      else if (dataType == FC_DT_CHAR){
		values[vertOffsets[n] + p] = ((char*)data)[p*numNodalComps[k]+m];
	      }
	    }
	  }
	  error = ex_put_one_attr(exoid, EX_NODAL, 0, current+1, values);
	  if (error != 0) {
	    fc_printfErrorMessage("Exodus failed to put nodal attr");
	    return FC_FILE_IO_ERROR;
	  }
	  current++;
	}
      }
      free(values);
    }

    // Cleanup                                                                               
    free(numNodalComps);
   for (j = 0; j < numNodalAttrs; j++) {
      for (k = 0; k < numMesh; k++)
	free(sortedVars[j][k]);
      free(sortedVars[j]);
    }
    free(sortedVars);
  } //numBlob > 0

  // delete nattr_sba (internal data was handed off above)
  fc_freeSortedBlobArray(&nattr_sba);


  // --- Write seq variables ---
  // FIX? committed flag and exoinfo stuff not set so cannot release ...
  //      this is hard to fix due to mulitcomponent variables

  // Sort the seq vars by assoc type, then other stuff, & lastely by meshID.
  fc_initSortedBlobArray(&sba);
  for (i = 0; i < numMesh; i++) {
    int numSeqVar, *numStepPerSeq;
    FC_Variable** seqVars;
    rc = fc_getSeqVariables(meshes[i], &numSeqVar, &numStepPerSeq, &seqVars);
    if (rc != FC_SUCCESS)
      return rc;
    for (j = 0; j < numSeqVar; j++) {
      FC_Sequence temp_sequence;
      rc = fc_getSequenceFromSeqVariable(numStepPerSeq[j], seqVars[j],
                                         &temp_sequence);
      if (rc != FC_SUCCESS)
        return rc;
      if (!FC_HANDLE_EQUIV(sequence, temp_sequence)) {
        fc_printfWarningMessage("Encountered var not on the first sequence, "
                                "ignoring")
        continue;
      }
      rc = fc_addBlobToSortedBlobArray(&sba, seqVars[j],
				       _fc_exo_seqVar_meshID_cmp); 
      if (rc < FC_SUCCESS) {
        fc_printfErrorMessage("failed to add seqvar to sba");
        return FC_ERROR;
      }
      // Warn about "duplicates" which are dropped ...
      if (rc == 0) {
	free(seqVars[j]);
	fc_printfWarningMessage("Encountered duplicate seqvar, ignoring");
      }
    }
    free(numStepPerSeq);
    free(seqVars);
  }

  if (sba.numBlob > 0) {
    // FIX: Could make this code clearer - actually have each assoc type
    // as separate array?
    //int total possible is = 7;
    int numAssocType = 0;
    int numVarPerAssocType[7] = { 0, 0, 0, 0, 0, 0 }; // num uniquely name seqvars
                                                   // of that type
    int numExoVarPerAssocType[7];
    FC_AssociationType lastAssoc = -1, assocTypes[7];
    char exoVarTypes[7][2];
    FC_Variable*** sortedSeqVars[7]; // [assocType][varIndexPerType][meshID]
    int* numComps[7]; // [assocType][varIndexPerType]
    FC_Variable* last_seqVar = NULL, *curr_seqVar;

    // count number of vars of each assoc type
    for (i = 0; i < sba.numBlob; i++) {
      FC_AssociationType assoc;
      curr_seqVar = (FC_Variable*)sba.blobs[i];
      rc = fc_getVariableAssociationType(curr_seqVar[0], &assoc);
      if (rc != FC_SUCCESS)
        return rc;
      // if of same assoc, only increment if really a new var
      if (assoc == lastAssoc) {
        if (_fc_exo_seqVar_cmp(curr_seqVar, last_seqVar) != 0)
          numVarPerAssocType[numAssocType-1]++;
      }
      // else, change the current association type and count as first
      else {
        numAssocType++;
        numVarPerAssocType[numAssocType-1]++;
        assocTypes[numAssocType-1] = assoc;
        lastAssoc = assoc;
      }
      last_seqVar = curr_seqVar;
    }
    // convert to array of var handles (w/ possible empty slots)
    // handing off responsibility for seqVar arrays
    k = 0;
    curr_seqVar = (FC_Variable*)sba.blobs[k];
    for (i = 0; i < numAssocType; i++) {
      sortedSeqVars[i] = (FC_Variable***)malloc(numVarPerAssocType[i]*
                                                sizeof(FC_Variable**));
      if (!sortedSeqVars[i]) {
        fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
        return FC_MEMORY_ERROR;
      }
      for (j = 0; j < numVarPerAssocType[i]; j++) {
        sortedSeqVars[i][j] = (FC_Variable**)calloc(numMesh,
                                                    sizeof(FC_Variable*));
        if (!sortedSeqVars[i][j]) {
          fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
          return FC_MEMORY_ERROR;
        }
        last_seqVar = curr_seqVar; // quaranteed to get 1st one
	while (k < sba.numBlob) {
          curr_seqVar = (FC_Variable*)sba.blobs[k];
          int meshID = _fc_getMeshSlot(_fc_getVarSlot(curr_seqVar[0])->mesh)
                                                ->fileInfo.exoInfo.meshID;
          if (_fc_exo_seqVar_cmp(curr_seqVar, last_seqVar) != 0)
            break; // move onto next var
          sortedSeqVars[i][j][meshID] = curr_seqVar;
          k++;
        }
      }
    }
    // delete sba (internal data was handed off above)
    fc_freeSortedBlobArray(&sba);


    // Write number and names of variables
    for (i = 0; i < numAssocType; i++) {
      int current;
      int numVar = numVarPerAssocType[i]; // # of fclib variables
      int numExoVar = numVar; // includes multicomponents as solos
      char** names, **temp_names;

      // determine type
      switch (assocTypes[i]) {
      case FC_AT_UNKNOWN:  strcpy(exoVarTypes[i], "");   break;
      case FC_AT_VERTEX:   strcpy(exoVarTypes[i], "n");  break;
      case FC_AT_EDGE:     strcpy(exoVarTypes[i], "");   break;
      case FC_AT_FACE:     strcpy(exoVarTypes[i], "");   break;
      case FC_AT_ELEMENT:  strcpy(exoVarTypes[i], "e");  break;
      case FC_AT_WHOLE_MESH:    strcpy(exoVarTypes[i], "");  break;
      case FC_AT_WHOLE_DATASET:    strcpy(exoVarTypes[i], "");  break;
      default:
        fc_printfErrorMessage("invalid assocation type encountered");
        return FC_ERROR;
      }

      // assemble names array
      // CHECK for multicomponent => have to write as separate components
      names = (char**)malloc(numVar*sizeof(char*)); // this might grow
      numComps[i] = (int*)malloc(numVar*sizeof(int)); // won't grow
      if (!names || !numComps[i]) {
        free(names); free(numComps[i]);
        fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
        return FC_MEMORY_ERROR;
      }
      current = 0;
      for (j = 0; j < numVar; j++) {
        int validID;
        int numComp;
        FC_MathType mathtype;

        // Find a valid seqVar (not all meshes may have a seqvar)
        validID = -1;
        for (k = 0; k < numMesh; k++) {
          if (sortedSeqVars[i][j][k]) {
            validID = k;
            break;
          }
        }
        if (validID < 0) {
          fc_printfErrorMessage("Programmer error?");
          return FC_ERROR;
        }

        // Add the var's name (or names if multicomponent)
        rc = fc_getVariableName(sortedSeqVars[i][j][validID][0], 
                                &temp_name);
        if (rc != FC_SUCCESS)
          return rc;
        rc = fc_getVariableInfo(sortedSeqVars[i][j][validID][0], NULL, 
                                &numComp, NULL, &mathtype, NULL);
        if (rc != FC_SUCCESS)
          return rc;
        numComps[i][j] = numComp; // save for later useage
        if (numComp > 1) {
          numExoVar += numComp - 1;  // -1 'cause 1 component already counted
          temp_names = realloc(names, numExoVar*sizeof(char*));
          if (!temp_names) {
            fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
            return FC_MEMORY_ERROR;
          }
          names = temp_names;
        }
        for (k = 0; k < numComp; k++) {
          rc = _fc_createExoVarCompName(temp_name, numComp, mathtype, k,
                                        &names[current]);
          if (rc != FC_SUCCESS)
            return rc;
          current++;
        }
        free(temp_name);
      }
      numExoVarPerAssocType[i] = numExoVar; // save for later
 
      // Process names - write to exo or explain why can't
      if (!strcmp(exoVarTypes[i], "")) {
        fc_printfWarningMessage("Exodus can't handle variables with "
                                "association type '%s'", 
                                fc_getAssociationTypeText(assocTypes[i]));
        for (j = 0; j < numExoVar; j++)
          fc_printfWarningMessage("  Discarding %s variable '%s'",
                                  fc_getAssociationTypeText(assocTypes[i]),
                                  names[j]);
      }
      else {
        error = ex_put_var_param(exoid, exoVarTypes[i], numExoVar);
        if (error != 0) {
          fc_printfErrorMessage("Failed to put '%s' var param", exoVarTypes[i]);
          return FC_FILE_IO_ERROR;
        }
        error = ex_put_var_names(exoid, exoVarTypes[i], numExoVar, names);
        if (error != 0) {
          fc_printfErrorMessage("Failed to put '%s' var names", exoVarTypes[i]);
          return FC_FILE_IO_ERROR;
        }
      }

      // cleanup names array
      for (j = 0; j < numExoVar; j++)
        free(names[j]);
      free(names);

      // Write Element var truth table if elements
      if (assocTypes[i] == FC_AT_ELEMENT) {
        int* truth_table;
        truth_table = calloc(numMesh*numExoVar, sizeof(int));
        if (!truth_table) {
          fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
          return FC_MEMORY_ERROR;
        }
        current = 0;
        for (j = 0; j < numMesh; j++) {
          for (k = 0; k < numVar; k++) {
            for (m = 0; m < numComps[i][k]; m++) {
              if (sortedSeqVars[i][k][j])
                truth_table[current] = 1;
              current++;
            }
          }
        }
        error = ex_put_elem_var_tab(exoid, numMesh, numExoVar, truth_table);
        if (error != 0) {
          fc_printfErrorMessage("Failed to put element var truth table");
          return FC_FILE_IO_ERROR;
        }
        free(truth_table);
      }
    }

    // write data
    for (i = 0; i < numAssocType; i++) {
      int numVar = numVarPerAssocType[i];
        
      // skip types we aren't writing
      if (!strcmp(exoVarTypes[i], "")) 
        continue;
      
      // do global vars
      else if (!strcmp(exoVarTypes[i], "g")) {
        double* values;
        void* data;
        int validIDs[numVarPerAssocType[i]];
        FC_DataType dataType;
        
        // Find 1st valid seqVar for each global var
        // FIX? warn that we are only taking the first one? - will become
        // moot if implement proper global vars
        for (j = 0; j < numVar; j++) {
          validIDs[j] = -1;
          for (k = 0; k < numMesh; k++) {
            if (sortedSeqVars[i][j][k]) {
              validIDs[j] = k;
              break;
            }
          }      
          if (validIDs[j] < 0) {
            fc_printfErrorMessage("Programmer error?");
            return FC_ERROR;
          }
        }

        // make exo values array
        values = malloc(numExoVarPerAssocType[i]*sizeof(double));
        if (!values) {
          fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
          return FC_MEMORY_ERROR;
        }
        // For each step, write all global vars
        for (j = 0; j < numStep; j++) {
          int current = 0;
          for (k = 0; k < numVar; k++) {
            fc_getVariableDataType(sortedSeqVars[i][k][validIDs[k]][j], &dataType);
            fc_getVariableDataPtr(sortedSeqVars[i][k][validIDs[k]][j], &data);
            for (m = 0; m < numComps[i][k]; m++) {
              if (dataType == FC_DT_CHAR)
                values[current] = ((char*)data)[m];
              else if (dataType == FC_DT_INT)
                values[current] = ((int*)data)[m];
              else if (dataType == FC_DT_FLOAT)
                values[current] = ((float*)data)[m];
              else if (dataType == FC_DT_DOUBLE)
                values[current] = ((double*)data)[m];
              current++;
            }
          }
          error = ex_put_glob_vars(exoid, j+1, numExoVarPerAssocType[i], 
                                   values);
          if (error != 0) {
            fc_printfErrorMessage("Exodus failed to put global vars for "
                                  "step %d", j);
            return FC_FILE_IO_ERROR;
          }
          // global vars not read lazily, so no committed flag
        }
        free(values);
      }

      // do nodal vars
      // write var for all nodes, per step
      else if (!strcmp(exoVarTypes[i], "n")) {
        double* values = malloc(numTotalVert*sizeof(double));
        if (!values) {
          fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
          return FC_MEMORY_ERROR;
        }
        void* data;
        FC_DataType dataType;
        for (j = 0; j < numStep; j++) {
          int current = 0;
          for (k = 0; k < numVarPerAssocType[i]; k++) {
            for (m = 0; m < numComps[i][k]; m++) {
              for (n = 0; n < numTotalVert; n++) 
                values[n] = 0; // initialize in case values don't exist
              for (n = 0; n < numMesh; n++) {
                // FIX add a warning for this - do it in names section above
                if (!sortedSeqVars[i][k][n]) 
                  continue; // skip if var doesn't exist == pad with zeros
                fc_getVariableDataType(sortedSeqVars[i][k][n][j], &dataType);
                fc_getVariableDataPtr(sortedSeqVars[i][k][n][j], &data);
                for (p = 0; p < numVerts[n]; p++) {
                  if (dataType == FC_DT_CHAR)
                    values[vertOffsets[n] + p] = ((char*)data)[p*numComps[i][k]+m];
                  else if (dataType == FC_DT_INT)
                    values[vertOffsets[n] + p] = ((int*)data)[p*numComps[i][k]+m];
                  else if (dataType == FC_DT_FLOAT)
                    values[vertOffsets[n] + p] = ((float*)data)[p*numComps[i][k]+m];
                  else if (dataType == FC_DT_DOUBLE)
                    values[vertOffsets[n] + p] = ((double*)data)[p*numComps[i][k]+m];
                }
              }
              error = ex_put_nodal_var(exoid, j+1, current+1, numTotalVert, values);
              if (error != 0) {
                fc_printfErrorMessage("Exodus failed to put nodal vars for "
                                      "step %d", j);
                return FC_FILE_IO_ERROR;
              }
              current++;
            }
          }
        }
        free(values);
      }

      // do element vars
      // write var per step per block
      else if (!strcmp(exoVarTypes[i], "e")) {
        int numDP;
        FC_DataType dataType;
        double* values;
        void* data;
        for (j = 0; j < numStep; j++) {
          for (m = 0; m < numMesh; m++) {
            int current = 0;
            for (k = 0; k < numVarPerAssocType[i]; k++) {
              if (!sortedSeqVars[i][k][m])
                current += numComps[i][k];
              else {
                fc_getVariableNumDataPoint(sortedSeqVars[i][k][m][j], &numDP);
                fc_getVariableDataType(sortedSeqVars[i][k][m][j], &dataType);
                fc_getVariableDataPtr(sortedSeqVars[i][k][m][j], &data);
                values = malloc(numDP*sizeof(double));
                if (!values) {
                  fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
                  return FC_MEMORY_ERROR;
                }
                for (n = 0; n < numComps[i][k]; n++) {
                  for (p = 0; p < numDP; p++) {
                    values[p] = 12;
                    if (dataType == FC_DT_CHAR)
                      values[p] = ((char*)data)[p*numComps[i][k]+n];
                    else if (dataType == FC_DT_INT)
                      values[p] = ((int*)data)[p*numComps[i][k]+n];
                    else if (dataType == FC_DT_FLOAT)
                      values[p] = ((float*)data)[p*numComps[i][k]+n];
                    else if (dataType == FC_DT_DOUBLE)
                      values[p] = ((double*)data)[p*numComps[i][k]+n];
                  }
                  // time id, var id, block id
                  error = ex_put_elem_var(exoid, j+1, current+1, blk_ids[m], 
					  numDP, values);
                  if (error != 0) {
                    fc_printfErrorMessage("Exodus failed to put elem vars for "
                                          "step %d on block %d", j, blk_ids[m]);
                    return FC_FILE_IO_ERROR;
                  }
                  current++;
                }
                free(values);
              }
            }
          }
        }
      }
    } // end of write data
    
    // Cleanup 
    for (i = 0; i < numAssocType; i++) {
      free(numComps[i]);
      for (j = 0; j < numVarPerAssocType[i]; j++) {
        for (k = 0; k < numMesh; k++) 
          free(sortedSeqVars[i][j][k]);
        free(sortedSeqVars[i][j]);
      }
      free(sortedSeqVars[i]);
    }
    // Do NOT free sortedSeqVars - it's an automatic variable
    
  } // done with seq vars

  // cleanup
  free(blk_ids);
 
  // Flush the buffers
  error = ex_update(exoid);
  if (error != 0) {
    fc_printfErrorMessage("Exodus failed to update at end of writing");
    return FC_FILE_IO_ERROR;
  }

  // Move these up in the file? .... delete when no longer used
  free(numDims);
  free(meshes);
  free(vertOffsets);
  free(elemOffsets);
  free(numVerts);
  free(numElems);

  dsSlot->fileType = FC_FT_EXODUS;
  dsSlot->header.committed = 1;
  dsSlot->writable = 0;

  return FC_SUCCESS;
}

/**
 * \ingroup  ExodusFileIO
 * \brief  Close a dataset.
 *
 * \description
 *  
 *    Do any necessary file related cleanup for the Dataset.
 *
 * \modifications  
 *   - 8/1/2005 WSD.
 */
FC_ReturnCode _fc_closeExodusDataset(
  FC_Dataset dataset /**< input - the dataset */
) {
  int error;
  int i;
  _FC_DsSlot* dsSlot;

  // check input
  dsSlot = _fc_getDsSlot(dataset);
  if (dsSlot == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // free dynamically allocate memory in db
  for (i = 0; i < dsSlot->fileInfo.exoInfo.numMesh; i++)
    free(dsSlot->fileInfo.exoInfo.vertIDs[i]);
  free(dsSlot->fileInfo.exoInfo.vertIDs);

  // close any files associated with this dataset
  error = ex_close(dsSlot->fileInfo.exoInfo.exoid);
  if (error != 0) 
    fc_printfWarningMessage("Exodus database for dataset '%s' did not close "
                            "properly", dsSlot->header.name);

  return FC_SUCCESS;
}


