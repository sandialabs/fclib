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
 * \file fileio.c
 * \brief Implementation for \ref FileIO module.
 * 
 * $Source: /home/Repositories/fcdmf/fclib/modules/fileio.c,v $
 * $Revision: 1.40 $ 
 * $Date: 2006/11/07 23:49:02 $
 */

// C library includes
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h> // for stat() call to see if vald file

// aux file types
#ifdef HAVE_LIBXML2
//#include <libxml/tree.h> // already included from fileioP.h
#include <libxml/parser.h>
#endif

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

// this module (specific io types embedded in fileioP.h"
#include "fileio.h"
#include "fileioP.h"

/**
 * \addtogroup FileIO
 * \brief Generic file IO.
 *
 * \description
 *
 *    This file provides wrappers to the various file IO options to hide the
 *    specific implementations (the specific implementations should generally
 *    not be used). We used to support certian SAF conventions, but
 *    no longer. For details about a specific type of IO
 *    see \ref ExodusFileIO \ref LSDynaFileIO.
 *
 * \modifications
 *    - 07/13/2005 WSD created to separate file IO out more so that it is
 *      easier to add new types.
 *    - 10/15/2007 ACG removed saf
 */

/**
 * \ingroup FileIO
 * \defgroup PrivateFileIO (Private Generic IO)
 * \brief Generic file IO.
 */

/**
 * \ingroup FileIO
 * \defgroup ExodusFileIO (Private Exodus IO)
 */

/**
 * \ingroup FileIO
 * \defgroup LSDynaFileIO (Private LS-Dyna IO)
 */

/**
 * \ingroup  PrivateFileIO
 * \brief  Read FC vertex file.
 *
 * \description  
 *
 *    This routine is used to bootstrap example dataset files. The file
 *    contains the coordinates of a set of points. The coords are returned
 *    interleaved, e.g. Vert0Coord0 Vert0Coord1 Vert0Coord2 Vert1Coord0
 *    Vert1Coord1 Vert1Coord2 ... VertNCoord0 VertNCoord1 VertNCoord2.
 *
 *    Currently this assumes 3D coordinates.
 *
 * \todo Modify this and _fc_readDataFile() to take numDim/numComponent.
 *
 * \modifications  
 *   - 8/24/2005 WSD, created.
 */
FC_ReturnCode _fc_readVertexFile(
  char* vertfilename, /**< input - file with vertex coords. */
  int* numVert,       /**< output - number of vertices. */
  int* numDim,        /**< output - dimensionality of the coords. */
  double** coords     /**< output - the coordinates of the vertes. */
) {
  int i, count, temp_numVert, temp_numDim;
  FILE* vertfile;
  double* temp_coords;

  // default return
  if (numVert)
    *numVert = -1;
  if (numDim)
    *numDim = -1;
  if (coords)
    *coords = NULL;

  // check inputs
  if (!vertfilename || !numVert || !numDim || !coords)
    return FC_INPUT_ERROR;
  
  // Attempt to open the files
  vertfile = fopen(vertfilename, "r");
  if (!vertfile) {
    fc_printfErrorMessage("Could not open vert file '%s'", vertfilename);
    return FC_FILE_IO_ERROR;
  }

  // Read metada
  temp_numDim = 3; // ASSUMPTION!
  count = fscanf(vertfile, "%d", &temp_numVert);
  if (count != 1 || temp_numVert < 1) {
    fc_printfErrorMessage("Error reading number of verts in vert file '%s'",
			  vertfilename);
    fclose(vertfile);
    return FC_ERROR;
  }

  // Read coords
  temp_coords = malloc(temp_numVert*temp_numDim*sizeof(double));
  if (!temp_coords) {
    fclose(vertfile);
    return FC_MEMORY_ERROR;
  }
  for (i = 0; i < temp_numVert*temp_numDim; i++) {
    if (fscanf(vertfile, "%lf", &temp_coords[i]) != 1) {
      fc_printfErrorMessage("error reading coords for vert %d", i/temp_numDim);
      free(temp_coords);
      fclose(vertfile);
      return FC_FILE_IO_ERROR;
    }
  }

  // All dones
  fclose(vertfile);
  *numVert = temp_numVert;
  *numDim = temp_numDim;
  *coords = temp_coords;

  return FC_SUCCESS;
}

/**
 * \ingroup  PrivateFileIO
 * \brief  Read FC element file.
 *
 * \description
 *  
 *    This routine is used to bootstrap example dataset files. The file
 *    contains the element connetivities for a mesh. The conns are 
 *    returned interleaved, e.g. Elem0Vert0 Elem0Vert1 .. Elem0VertN
 *    Elem1Vert0 Elem1Vert1 .. Elem1VertN ... ElemMVert0 ElemMVert1 ..
 *    ElemMVertN.
 *
 * \modifications  
 *   - 8/24/2005 WSD, created.
 */
FC_ReturnCode _fc_readElementFile(
  char* elemfilename, /**< input - file with element conns */
  int* numVert,       /**< output - number of vertices */
  int* numElem,       /**< output - number of elements */
  FC_ElementType* elemType, /**< output - element type */
  char** name,  /**< output - name of mesh (user responsible for freeing */
  int** conns         /**< output - element conns */
) {
  int i, count;
  int numVertPerElem, temp_numVert, temp_numElem;
  FILE* elemfile;
  char temp_name[1024];
  FC_ElementType temp_elemType;
  int* temp_conns;

  // default return
  if (numVert)
    *numVert = -1;
  if (numElem)
    *numElem = -1;
  if (elemType)
    *elemType = FC_ET_UNKNOWN;
  if (name)
    *name = NULL;
  if (conns)
    *conns = NULL;

  // check inputs
  if (!elemfilename || !numVert || !numElem || !elemType || !name || !conns)
    return FC_INPUT_ERROR;
  
  // Attempt to open the file
  elemfile = fopen(elemfilename, "r");
  if (!elemfile) {
    fc_printfErrorMessage("Could not open elem file '%s'", elemfilename);
    return FC_FILE_IO_ERROR;
  }

  // Read metadata
  // FIX? don't have to put in mesh name?
  count = fscanf(elemfile, "%d%d", &temp_numVert, &temp_numElem);
  if (count != 2 || temp_numVert < 1 || temp_numElem < 1) {
    fc_printfErrorMessage("error reading number of verts and elems in "
			  "elem file '%s'", elemfilename);
    fclose(elemfile);
    return FC_ERROR;
  }
  fscanf(elemfile, "%s\n", temp_name);
  if (!strcmp(temp_name, "POINT"))
    temp_elemType = FC_ET_POINT;
  else if (!strcmp(temp_name, "LINE"))
    temp_elemType = FC_ET_LINE;
  else if (!strcmp(temp_name, "TRI"))
    temp_elemType = FC_ET_TRI;
  else if (!strcmp(temp_name, "QUAD"))
    temp_elemType = FC_ET_QUAD;
  else if (!strcmp(temp_name, "TET"))
    temp_elemType = FC_ET_TET;
  else if (!strcmp(temp_name, "PYRAMID"))
    temp_elemType = FC_ET_PYRAMID;
  else if (!strcmp(temp_name, "PRISM"))
    temp_elemType = FC_ET_PRISM;
  else if (!strcmp(temp_name, "HEX"))
    temp_elemType = FC_ET_HEX;
  else {
    fc_printfErrorMessage("unknown element type encountered '%s'",
			  temp_name);
    fclose(elemfile);
    return FC_ERROR;
  }
  numVertPerElem = fc_getElementTypeNumVertex(temp_elemType);
  fscanf(elemfile, "%[^\n]", temp_name); // now the mesh name

  // Read conns
  temp_conns = malloc(temp_numElem*numVertPerElem*sizeof(int));
  if (!temp_conns) {
    fclose(elemfile);
    return FC_MEMORY_ERROR;
  }
  for (i = 0; i < temp_numElem*numVertPerElem; i++) {
    if (fscanf(elemfile, "%d", &temp_conns[i]) != 1) {
      fc_printfErrorMessage("error reading conns for element %d",
			    i/numVertPerElem);
      free(temp_conns);
      fclose(elemfile);
      return FC_FILE_IO_ERROR;
    }
  }
  fclose(elemfile);

  // All done
  *name = malloc((strlen(temp_name)+1)*sizeof(char));
  if (!(*name)) {
    free(temp_conns);
    return FC_MEMORY_ERROR;
  }
  strcpy(*name, temp_name);
  *numVert = temp_numVert;
  *numElem = temp_numElem;
  *elemType = temp_elemType;
  *conns = temp_conns;
  
  return FC_SUCCESS;
}

/**
 * \ingroup  PrivateFileIO
 * \brief  Read FC data file.
 *
 * \description
 *  
 *    This routine is used to bootstrap example dataset files. The file
 *    contains data values for a mesh. The values are 
 *    returned interleaved, e.g. Point0Comp0 Point0Comp1 .. Point0CompN
 *    Point1Comp0 Point1Comp1 .. Point1CompN ... PointMComp0 PointMComp1 ..
 *    PointMCompN.
 *
 *    Only supports scalars, 3D vectors and 3 component symmetric tensors.
 *
 * \todo Make this read doubles? or make able to read multiple types
 *
 * \modifications  
 *   - 8/24/2005 WSD, created.
 *   - 5/11/2006 WSD modified to read vectors and 3 component symtensor.
 */
FC_ReturnCode _fc_readDataFile(
  char* datafilename, /**< input - file with data values */
  int* numPoint,      /**< output - number of data points */
  int* numComp,       /**< output - number of components */
  FC_AssociationType* assoc, /**< output - how data is associated with mesh */
  FC_MathType* mathType, /**< output - math type of data (scalar or vector) */
  char** name,  /**< output - name of variable (user responsible for freeing) */
  float** data        /**< output - data values */
) {
  int i, count;
  int temp_numPoint, temp_numComp;
  FILE* datafile;
  char temp_name[1024];
  float* temp_data;
  FC_AssociationType temp_assoc;
  FC_MathType temp_mathType;

  // default return
  if (numPoint)
    *numPoint = -1;
  if (numComp)
    *numComp = -1;
  if (assoc)
    *assoc = FC_AT_UNKNOWN;
  if (mathType)
    *mathType = FC_MT_UNKNOWN;
  if (name)
    *name = NULL;
  if (data)
    *data = NULL;

  // check inputs
  if (!datafilename || !numPoint || !numComp || !assoc || !mathType ||
      !name || !data)
    return FC_INPUT_ERROR;
  
  // Attempt to open the file
  datafile = fopen(datafilename, "r");
  if (!datafile) {
    fc_printfErrorMessage("Could not open data file '%s'", datafilename);
    return FC_FILE_IO_ERROR;
  }

  // Read metadata
  count = fscanf(datafile, "%d", &temp_numPoint);
  if (count != 1 || temp_numPoint < 1) {
    fc_printfErrorMessage("error reading number of points in data "
			  "file '%s'", datafilename);
    fclose(datafile);
    return FC_ERROR;
  }
  fscanf(datafile, "%s\n", temp_name);
  if (!strcmp(temp_name, "VERTEX"))
    temp_assoc = FC_AT_VERTEX;
  else if (!strcmp(temp_name, "ELEMENT"))
    temp_assoc = FC_AT_ELEMENT;
  else {
    fc_printfErrorMessage("unknown assoc type encountered '%s'",
			  temp_name);
    fclose(datafile);
    return FC_ERROR;
  }
  fscanf(datafile, "%s\n", temp_name);
  if (!strcmp(temp_name, "SCALAR")) {
    temp_numComp = 1; 
    temp_mathType = FC_MT_SCALAR;
  }
  else if (!strcmp(temp_name, "VECTOR")) {
    temp_numComp = 3;  // ASSUMPTION!
    temp_mathType = FC_MT_VECTOR;
  }
  else if (!strcmp(temp_name, "SYMTENSOR")) {
    temp_numComp = 3;  // ASSUMPTION!
    temp_mathType = FC_MT_SYMTENSOR;
  }
  else {
    fc_printfErrorMessage("unknown math type encountered '%s'",
			  temp_name);
    fclose(datafile);
    return FC_ERROR;
  }
  fscanf(datafile, "%[^\n]", temp_name); // now the mesh name

  // Read conns
  temp_data = malloc(temp_numPoint*temp_numComp*sizeof(float));
  if (!temp_data) {
    fclose(datafile);
    return FC_MEMORY_ERROR;
  }
  for (i = 0; i < temp_numPoint*temp_numComp; i++) {
    if (fscanf(datafile, "%f", &temp_data[i]) != 1) {
      fc_printfErrorMessage("error reading data for point %d",
			    i/temp_numComp);
      free(temp_data);
      fclose(datafile);
      return FC_FILE_IO_ERROR;
    }
  }
  fclose(datafile);

  // All done
  *name = malloc((strlen(temp_name)+1)*sizeof(char));
  if (!(*name)) {
    free(temp_data);
    return FC_MEMORY_ERROR;
  }
  *numComp = temp_numComp;
  *numPoint = temp_numPoint;
  *assoc = temp_assoc;
  *mathType = temp_mathType;
  strcpy(*name, temp_name);
  *data = temp_data;
  
  return FC_SUCCESS;
}

/**
 * \ingroup  PrivateFileIO
 * \brief  Read in vert and elem files and make a mesh.
 *
 * \description
 *  
 *    This routine is used to bootstrap example dataset files.
 *
 * \modifications  
 *   - 8/24/2005 WSD, created.
 */
FC_ReturnCode _fc_makeMeshFromFiles(
  FC_Dataset dataset, /**< input - dataset to add mesh too. */
  char* vertfilename, /**< input - file with vertex coords */
  char* elemfilename, /**< input - file with element conns */
  FC_Mesh* mesh       /**< output - the created mesh. */
) {
  FC_ReturnCode rc;
  int numVert, numElem, temp_numVert, numDim;
  char* name;
  FC_ElementType elemType;
  double* coords;
  int* conns;

  // default return
  if (mesh)
    *mesh = FC_NULL_MESH;

  // check inputs
  if (!fc_isDatasetValid(dataset) || !vertfilename || !elemfilename || !mesh)
    return FC_INPUT_ERROR;
  
  // Get info we need
  rc = _fc_readVertexFile(vertfilename, &numVert, &numDim, &coords);
  if (rc != FC_SUCCESS)
    return rc;
  rc = _fc_readElementFile(elemfilename, &temp_numVert, &numElem, &elemType,
			   &name, &conns);
  if (rc != FC_SUCCESS) {
    free(coords);
    return rc;
  }

  // a little more checking
  if (temp_numVert != numVert) {
    fc_printfErrorMessage("The number of verts in vertfile '%s' and "
			  "elemfile '%s' do not agree (%d != %d)",
			  vertfilename, elemfilename, numVert, temp_numVert);
    free(coords);
    free(conns);
    free(name);
    return FC_ERROR;
  }

  // Make the mesh
  rc = fc_createMesh(dataset, name, mesh);
  free(name);
  if (rc != FC_SUCCESS) {
    free(coords);
    free(conns);
    return rc;
  }
  rc = fc_setMeshCoordsPtr(*mesh, numDim, numVert, coords);
  if (rc != FC_SUCCESS) {
    free(coords);
    free(conns);
    return rc;
  }
  rc = fc_setMeshElementConnsPtr(*mesh, elemType, numElem, conns);
  if (rc != FC_SUCCESS) {
    free(conns);
    return rc;
  }

  return FC_SUCCESS;
}

/**
 * \ingroup  PrivateFileIO
 * \brief  Read in data file and add to a mesh.
 *
 * \description
 *  
 *    This routine is used to bootstrap example dataset files.
 *
 * \modifications  
 *   - 8/24/2005 WSD, created.
 */
FC_ReturnCode _fc_makeVariableFromFile(
  FC_Mesh mesh,          /**< input - mesh to add data too. */
  char* datafilename,    /**< input - file with data values */
  FC_Variable* variable  /**< output - the created variable. */
) {
  FC_ReturnCode rc;
  int numPoint, numComp, numEntity;
  char* name;
  FC_AssociationType assoc;
  FC_MathType mathType;
  float* data;

  // default return
  if (variable)
    *variable = FC_NULL_VARIABLE;

  // check inputs
  if (!fc_isMeshValid(mesh) || !datafilename || !variable)
    return FC_INPUT_ERROR;
  
  // Get info we need
  rc = _fc_readDataFile(datafilename, &numPoint, &numComp, &assoc, &mathType, 
			&name, &data);
  if (rc != FC_SUCCESS)
    return rc;

  // A little more checking
  rc = fc_getMeshNumEntity(mesh, assoc, &numEntity);
  if (rc != FC_SUCCESS)
    return rc;
  if (numEntity != numPoint) {
    fc_printfErrorMessage("The number of mesh entities in datafile '%s' "
			  "does not agree with the mesh (%d != %d)",
			  datafilename, numPoint, numEntity);
    free(data);
    free(name);
    return FC_ERROR;
  }

  // Make the variable
  rc = fc_createVariable(mesh, name, variable);
  free(name);
  if (rc != FC_SUCCESS) {
    free(data);
    return rc;
  }
  rc = fc_setVariableDataPtr(*variable, numPoint, numComp, assoc, mathType, 
			     FC_DT_FLOAT, data);
  if (rc != FC_SUCCESS) {
    free(data);
    return rc;
  }

  return FC_SUCCESS;
}

#ifdef HAVE_LIBXML2

/**
 * \ingroup PrivateFileIO
 * \brief Write header of an already opened xml file
 *
 * \todo ? reference dataset? date?
 *
 * \modifications
 *    - 06/05/2006 WSD Created.
 */
FC_ReturnCode _fc_writeAuxFileHeader(
  FILE* file /**< input - writable file handle */
) {
  fprintf(file, "<?xml version=\"1.0\"?>\n");
  fprintf(file, "<FC_AuxFile>\n");
  return FC_SUCCESS;
}

/**
 * \ingroup PrivateFileIO
 * \brief Write footer into xml file
 *
 * \modifications
 *    - 06/05/2006 WSD Created.
 */
FC_ReturnCode _fc_writeAuxFileFooter(
  FILE* file /**< input - writable file handle */
) {
  fprintf(file, "</FC_AuxFile>\n");
  return FC_SUCCESS;
}

/**
 * \ingroup PrivateFileIO
 * \brief Write a subset element to a subset file
 *
 * \modifications
 *    - 06/05/2006 WSD Created.
 */
FC_ReturnCode _fc_writeAuxFileSubset(
  FILE* file,       /**< input - writable file handle */
  FC_Subset subset  /**< input - the subset */
) {
  FC_ReturnCode rc;
  FC_Mesh mesh;
  char* meshName, *subsetName;
  FC_AssociationType assoc;
  int numMember, maxNumMember, *memberIDs;
  
  // check input
  if (!fc_isSubsetValid(subset))
    return FC_ERROR;

  // get info
  rc = fc_getMeshFromSubset(subset, &mesh);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_getMeshName(mesh, &meshName);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_getSubsetName(subset, &subsetName);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_getSubsetInfo(subset, &numMember, &maxNumMember, &assoc);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_getSubsetMembersAsArray(subset, &numMember, &memberIDs);
  if (rc != FC_SUCCESS)
    return rc;

  // do it
  rc = _fc_writeAuxFileSubsetCore(file, subsetName, meshName, assoc, 
				  maxNumMember, numMember, memberIDs);
  free(meshName);
  free(subsetName);
  free(memberIDs);
  return rc;
}

/**
 * \ingroup PrivateFileIO
 * \brief Write a subset to a subset file
 *
 * \modifications
 *    - 06/05/2006 WSD Created.
 */
FC_ReturnCode _fc_writeAuxFileSubsetCore(
  FILE* file,                 /**< input - writable file handle */
  char* subsetName,           /**< input - name of subset */
  char* meshName,             /**< input - name of parent mesh */
  FC_AssociationType assoc,   /**< input - association type of subset */
  int maxNumMember,           /**< input - max possible number of members */
  int numMember,              /**< input - number of members */
  int* memberIDs              /**< input - array of member IDs */
) {
  int i;
  char* typeName = "";

  // check input - minimal amount because we don't expect external use)
  if (!file || !subsetName || !meshName)
    return FC_INPUT_ERROR;

  // setup
  switch(assoc) {
  case FC_AT_UNKNOWN: typeName = "unknown"; break;
  case FC_AT_VERTEX:  typeName = "vertex";  break;
  case FC_AT_EDGE:    typeName = "edge";    break;
  case FC_AT_FACE:    typeName = "face";    break;
  case FC_AT_ELEMENT: typeName = "element"; break;
  case FC_AT_WHOLE_MESH:   typeName = "whole";   break;
  case FC_AT_WHOLE_DATASET: 
    return FC_ERROR;
  }

  // do it - write an XML element
  fprintf(file, "<FC_Subset\n");
  fprintf(file, "  Name=\"%s\"\n", subsetName);
  fprintf(file, "  ParentMeshName=\"%s\"\n", meshName);
  fprintf(file, "  Type=\"%s\"\n", typeName);
  fprintf(file, "  MaxNumMember=\"%d\"\n", maxNumMember);
  fprintf(file, "  NumMember=\"%d\"\n", numMember);
  fprintf(file, "  MemberIDs=\"%d", memberIDs[0]);
  for (i = 1; i < numMember; i++)
    fprintf(file, " %d", memberIDs[i]);
  fprintf(file, "\"\n");
  fprintf(file, "/>\n");
  
  return FC_SUCCESS;
}

/**
 * \ingroup PrivateFileIO
 * \brief Write a tear to a subset file
 *
 * \description
 *
 *    The tear is required to have at least one subset. All subset names
 *    have to be non-NULL.
 *
 * \todo How many digits of a double should we print out?
 *
 * \modifications
 *    - 06/05/2006 WSD Created.
 */
FC_ReturnCode _fc_writeAuxFileTear(
  FILE* file,           /**< input - writable file handle */
  char* tearName,       /**< input - name of tear */
  int numSubset,        /**< input - number of subsets that make up tear */
  char** subsetNames,   /**< input - names of the subsets */
  double length         /**< input - length of tear */
) {
  int i;

  // check input - minimal amount because we don't expect external use)
  if (!file || !tearName || numSubset < 1 || !subsetNames)
    return FC_INPUT_ERROR;
  for (i = 0; i < numSubset; i++)
    if (!subsetNames[i])
      return FC_INPUT_ERROR;

  // do it - write an XML element
  fprintf(file, "<FC_Tear\n");
  fprintf(file, "  Name=\"%s\"\n", tearName);
  fprintf(file, "  Length=\"%.11g\"\n", length);
  fprintf(file, "  NumSubset=\"%d\">\n", numSubset);
  for (i = 0; i < numSubset; i++)
    fprintf(file, "  <SubsetName>%s</SubsetName>\n", subsetNames[i]);
  fprintf(file, "</FC_Tear>\n");
  
  return FC_SUCCESS;
}

/**
 * \ingroup PrivateFileIO
 * \brief Import an aux file to create an xmlDoc.
 *
 * \description
 *
 *    This reads the xml file and creats a libxml2 xmlDoc object.
 *    The user is responsible for freeing the xmlDoc using xmlFreeDoc().
 *
 * \todo Can we control error reporting from libxml2?
 *
 * \modifications
 *    - 06/27/2006 WSD Created.
 */
FC_ReturnCode _fc_importAuxFileXMLDoc(
  char* filename,      /**< input - file name */
  xmlDocPtr* xmldoc    /**< output - pointer to the xml Doc object */
) {
  xmlDocPtr doc;
  xmlNode* root = NULL;
   
  // default return
  if (xmldoc)
    *xmldoc = NULL;

  // test input
  if (!filename || !xmldoc)
    return FC_INPUT_ERROR;

  // parse the file
  doc = xmlParseFile(filename);
  if (doc == NULL) 
    return FC_FILE_IO_ERROR;

  // must be an aux doc
  root = xmlDocGetRootElement(doc);
  if (!root || !root->name || xmlStrcmp(root->name,(const xmlChar*)"FC_AuxFile")) {
    xmlFreeDoc(doc);
    fc_printfErrorMessage("XML document is not of type FC_AuxFile");
    return FC_ERROR;
  }

  // return
  *xmldoc = doc;
  return FC_SUCCESS;
}

/**
 * \ingroup PrivateFileIO
 * \brief Get subsets from an aux file.
 *
 * \description
 *
 *    This creates a bunch of data which the user has the responsibility
 *    of freeing. The xml file has to have already been parsed into
 *    an xmlDocPtr.
 *
 *    NOTE: we're assuming that xmlChar* and char* can be cast back and
 *    forth w/o worry. That is, libxml2 uses UTF-8 and we cast to
 *    char* hoping that our words are the proper size.
 *
 * \modifications
 *    - 06/27/2006 WSD Created.
 */
FC_ReturnCode _fc_getAuxFileSubsets(
  xmlDocPtr doc,    /**< input - xml Doc ptr */
  int* numSubset,      /**< output - the number of subset */
  char*** subsetNames, /**< output - the names of each subset */
  char*** meshNames,   /**< output - the names of the parent meshes */
  FC_AssociationType** assocs, /**< output - the association types */
  int** maxNumMembers, /**< output - the max possible number of members */
  int** numMembers,    /**< output - the number of members */
  int*** memberIDs     /**< output - the member IDs */
) {
  int i;
  xmlNode* root = NULL, *cur_node;
  char* prop, *token;
  int num = 0, len = 0;
  char** sNames = NULL, **mNames = NULL;
  FC_AssociationType* asss = NULL;
  int* maxs = NULL, *nums = NULL, **ids = NULL;
  void* temp1, *temp2, *temp3, *temp4, *temp5, *temp6;

  // test input
  if (!doc || !numSubset || !subsetNames || !meshNames || !assocs ||
      !maxNumMembers || !numMembers || !memberIDs )
    return FC_INPUT_ERROR;

  // must be an aux doc
  root = xmlDocGetRootElement(doc);
  if (!root || !root->name || xmlStrcmp(root->name,(const xmlChar*)"FC_AuxFile")) {
    fc_printfErrorMessage("XML document is not of type FC_AuxFile");
    return FC_ERROR;
  }

  // walk the tree looking for subsets
  for (cur_node = root->children; cur_node != NULL; cur_node = cur_node->next)
  {
    if (cur_node->type == XML_ELEMENT_NODE && 
	!xmlStrcmp(cur_node->name, (const xmlChar*)"FC_Subset")) {
      if (num == len) {
	if (len == 0)
	  len = 2;
	else
	  len *= 2;
	temp1 = realloc(sNames, len*sizeof(char*));
	temp2 = realloc(mNames, len*sizeof(char*));
	temp3 = realloc(asss, len*sizeof(FC_AssociationType));
	temp4 = realloc(maxs, len*sizeof(int));
	temp5 = realloc(nums, len*sizeof(int));
	temp6 = realloc(ids, len*sizeof(int*));
	if (!temp1 || !temp2 || !temp3 || !temp4 || !temp5 || !temp6)
	  return FC_MEMORY_ERROR;
	sNames = temp1;
	mNames = temp2;
	asss = temp3;
	maxs = temp4;
	nums = temp5;
	ids = temp6;
      }
      sNames[num] = (char*) xmlGetProp(cur_node, (const xmlChar*)"Name");
      mNames[num] = (char*) xmlGetProp(cur_node, 
				       (const xmlChar*)"ParentMeshName");
      prop = (char*) xmlGetProp(cur_node, (const xmlChar*)"Type");
      if (!strcmp(prop, "vertex"))
	asss[num] = FC_AT_VERTEX;
      else if (!strcmp(prop, "edge"))
	asss[num] = FC_AT_EDGE;
      else if (!strcmp(prop, "face"))
	asss[num] = FC_AT_FACE;
      else if (!strcmp(prop, "element"))
	asss[num] = FC_AT_ELEMENT;
      else if (!strcmp(prop, "whole"))
	asss[num] = FC_AT_WHOLE_MESH;
      else if (!strcmp(prop,"unknown"))
	asss[num] = FC_AT_UNKNOWN;
      else {
	fc_printfWarningMessage("nonvalid assocation type '%s' encountered",
				prop);
	asss[num] = FC_AT_UNKNOWN;
      }
      free(prop);
      prop = (char*) xmlGetProp(cur_node, (const xmlChar*)"MaxNumMember");
      maxs[num] = atoi(prop);
      free(prop);
      prop = (char*) xmlGetProp(cur_node, (const xmlChar*)"NumMember");
      nums[num] = atoi(prop);
      free(prop);
      ids[num] = calloc(nums[num], sizeof(int));
      if (!ids[num])
	return FC_MEMORY_ERROR;
      prop = (char*) xmlGetProp(cur_node, (const xmlChar*)"MemberIDs");
      token = strtok(prop," ,\t\n");
      i = 0;
      while (token != NULL && i < nums[num]) {
	ids[num][i] = atoi(token);
	token = strtok(NULL, " ,\t\n");
	i++;
      }
      free(prop);
      if (i > nums[num])
	fc_printfWarningMessage("MemberIDs attribute had more members than specified");
      num++;
    }
  }

  // return values
  *numSubset = num;
  *subsetNames = sNames;
  *meshNames = mNames;
  *assocs = asss;
  *maxNumMembers = maxs;
  *numMembers = nums;
  *memberIDs = ids;

  return FC_SUCCESS;
}
/**
 * \ingroup PrivateFileIO
 * \brief Get tears from an aux file.
 *
 * \description
 *
 *    This creates a bunch of data which the user has the responsibility
 *    of freeing. The xml file has to have already been parsed into
 *    an xmlDocPtr.
 *
 *    NOTE: we're assuming that xmlChar* and char* can be cast back and
 *    forth w/o worry. That is, libxml2 uses UTF-8 and we cast to
 *    char* hoping that our words are the proper size.
 *
 * \modifications
 *    - 07/03/2006 WSD Created.
 */
FC_ReturnCode _fc_getAuxFileTears(
  xmlDocPtr doc,        /**< input - xml Doc ptr */
  int* numTear,            /**< output - the number of tears */
  char*** tearNames,       /**< output - the names of each tear */
  int** numSubsetPerTear,  /**< output - the number of subsets in each tear */
  char**** subsetNames,    /**< output - the names of the subsets */
  double** lengths         /**< output - the lengths of the tears */ 
) {
  xmlNode* root = NULL, *cur_node, *cur_child_node;
  char* prop;
  int num = 0, len = 0, count;
  char** tNames = NULL, ***sNames = NULL;
  double* lens = NULL;
  int *numSubsets = NULL;
  void* temp1, *temp2, *temp3, *temp4;

  // test input
  if (!doc || !numTear || !tearNames || ! numSubsetPerTear || !subsetNames)
    return FC_INPUT_ERROR;

  // must be an aux doc
  root = xmlDocGetRootElement(doc);
  if (!root || !root->name || xmlStrcmp(root->name,(const xmlChar*)"FC_AuxFile")) {
    fc_printfErrorMessage("XML document is not of type FC_AuxFile");
    return FC_ERROR;
  }
  
  // walk the tree looking for subsets
  for (cur_node = root->children; cur_node != NULL; 
       cur_node = cur_node->next) {
    if (cur_node->type == XML_ELEMENT_NODE && 
	!xmlStrcmp(cur_node->name, (const xmlChar*)"FC_Tear")) {
      if (num == len) {
	if (len == 0)
	  len = 2;
	else
	  len *= 2;
	temp1 = realloc(tNames, len*sizeof(char*));
	temp2 = realloc(sNames, len*sizeof(char**));
	temp3 = realloc(numSubsets, len*sizeof(int));
	temp4 = realloc(lens, len*sizeof(double));
	if (!temp1 || !temp2 || !temp3 || !temp4)
	  return FC_MEMORY_ERROR;
	tNames = temp1;
	sNames = temp2;
	numSubsets = temp3;
	lens = temp4;
      }
      tNames[num] = (char*) xmlGetProp(cur_node, (const xmlChar*)"Name");
      prop = (char*) xmlGetProp(cur_node, (const xmlChar*)"NumSubset");
      numSubsets[num] = atoi(prop);
      free(prop);
      sNames[num] = (char**) calloc(numSubsets[num], sizeof(char*));
      if (!sNames[num])
	return FC_MEMORY_ERROR;
      prop = (char*) xmlGetProp(cur_node, (const xmlChar*)"Length");
      lens[num] = atof(prop);
      free(prop);

      count = 0;
      for (cur_child_node = cur_node->children; 
	   cur_child_node != NULL; cur_child_node = cur_child_node->next) {
	if (cur_child_node->type == XML_ELEMENT_NODE &&
	    !xmlStrcmp(cur_child_node->name, (const xmlChar*)"SubsetName")) {
	  if (count == numSubsets[num]) {
	    fc_printfWarningMessage("Tear %s has more subset subset names than"  
				    " expected (%d != %d)", tNames[num],
				    count+1, numSubsets[num]);
	  }
	  else {
	    prop = (char*) xmlNodeGetContent(cur_child_node);
	    sNames[num][count] = prop;
	  }
	  count++;
	}
      } // loop over child nodes

      num++;
    }
  }

  // return values
  *numTear = num;
  *tearNames = tNames;
  *subsetNames = sNames;
  *numSubsetPerTear = numSubsets;
  *lengths = lens;

  return FC_SUCCESS;
}

/**
 * \ingroup PrivateFileIO
 * \brief Write header of an already opened bb file
 *
 * \modifications
 *    - 06/29/2006 WSD Created.
 */
FC_ReturnCode _fc_writeBBFileHeader(
  FILE* bb_file /**< input - writable file handle */
) {
  //fprintf(bb_file, "<?xml version=\"1.0\"?>\n");
  fprintf(bb_file, "<BoundingBoxFile>\n");
  return FC_SUCCESS;
}

/**
 * \ingroup PrivateFileIO
 * \brief Write footer into xml file
 *
 * \modifications
 *    - 06/29/2006 WSD Created.
 */
FC_ReturnCode _fc_writeBBFileFooter(
				    FILE* bb_file /**< */
) {
  fprintf(bb_file, "</BoundingBoxFile>\n");
  return FC_SUCCESS;
}

/**
 * \ingroup PrivateFileIO
 * \brief Write a bounding box element
 *
 * \description
 *   
 *    If numDim < 3, the extra dimensions are padded with zeros.
 *
 * \modifications
 *    - 06/29/2006 WSD Created.
 */
FC_ReturnCode _fc_writeBBFileBoundingBox(
  FILE* bb_file, /**< input - writable file handle */
  char* name,    /**< input - the name of the bounding box */
  int stepID,    /**< input - the step ID */
  char* comment_str, /**< input - comments */
  int numDim,    /**< input - number of dimensions */
  FC_Coords lowers, /**< input - the lower corner coords */
  FC_Coords uppers  /**< input - the upper corner coords */
) {
  int i;
  char axisName[3] = { 'X', 'Y', 'Z' };

  fprintf(bb_file, "  <BoundingBox\n");
  fprintf(bb_file, "    Name=\"%s\" Timestep=\"%d\"\n", name, stepID);
  fprintf(bb_file, "    Comments=\"%s\"\n", comment_str);
  for (i = 0; i < numDim; i++)
    fprintf(bb_file, "    Min%c=\"%g\" Max%c=\"%g\"\n", axisName[i], lowers[i],
	    axisName[i], uppers[i]);
  for (i = numDim; i < 3; i++)
    fprintf(bb_file, "    Min%c=\"0\" Max%c=\"0\"\n", axisName[i], axisName[i]);
  fprintf(bb_file, "  />\n");
  return FC_SUCCESS;
}

/**
 * \ingroup PrivateFileIO
 * \brief Read an xml bounding box file
 *
 * \description
 *
 *    This creates a bunch of data which the user has the responsibility
 *    of freeing.
 *
 *    NOTE: we're assuming that xmlChar* and char* can be cast back and
 *    forth w/o worry. That is, libxml2 uses UTF-8 and we cast to
 *    char* hoping that our words are the proper size.
 *
 * \modifications
 *    - 06/29/2006 WSD Created.
 */
FC_ReturnCode _fc_readBBFile(
  char* filename,       /**< input - file name */
  int* numBB,           /**< output - the number of bounding boxes */
  char*** bbNames,      /**< output - the names of the bbs */
  int** bbStepIDs,      /**< output - step IDs */
  char*** bbComments,   /**< output - comment strings */
  FC_Coords** bbLowers, /**< output - lower corners of the bbs */
  FC_Coords** bbUppers  /**< output - upper corners of the bb */
) {
  int i;
  xmlDocPtr doc;
  xmlNode* root = NULL, *cur_node;
  char* prop;
  int num = 0, len = 0;
  char** names = NULL, **comments = NULL;
  int* stepIDs = NULL;
  FC_Coords* lowers = NULL, *uppers = NULL;
  void* temp1, *temp2, *temp3, *temp4, *temp5;
  char* MinAttNames[3] = { "MinX", "MinY", "MinZ" };
  char* MaxAttNames[3] = { "MaxX", "MaxY", "MaxZ" };

  // test input
  if (!filename || !numBB || !bbNames || !bbComments || !bbStepIDs || 
      !bbLowers || !bbUppers)
    return FC_INPUT_ERROR;

  // parse the doc
  doc = xmlParseFile(filename);
  if (doc == NULL) 
    return FC_FILE_IO_ERROR;

  // must be an bounding box file
  root = xmlDocGetRootElement(doc);
  if (!root || !root->name || 
      xmlStrcmp(root->name,(const xmlChar*)"BoundingBoxFile")) {
    xmlFreeDoc(doc);
    fc_printfErrorMessage("XML document is not of type BoundingBoxFile");
    return FC_ERROR;
  }

  // walk the tree looking for bounding boxes
  for (cur_node = root->children; cur_node != NULL; cur_node = cur_node->next)
  {
    if (cur_node->type == XML_ELEMENT_NODE && 
	!xmlStrcmp(cur_node->name, (const xmlChar*)"BoundingBox")) {
      if (num == len) {
	if (len == 0)
	  len = 2;
	else
	  len *= 2;
	temp1 = realloc(names, len*sizeof(char*));
	temp2 = realloc(comments, len*sizeof(char*));
	temp3 = realloc(stepIDs, len*sizeof(int));
	temp4 = realloc(lowers, len*sizeof(FC_Coords));
	temp5 = realloc(uppers, len*sizeof(FC_Coords));
	if (!temp1 || !temp2 || !temp3 || !temp4 || !temp5)
	  return FC_MEMORY_ERROR;
	names = temp1;
	comments = temp2;
	stepIDs = temp3;
	lowers = temp4;
	uppers = temp5;
      }
      names[num] = (char*) xmlGetProp(cur_node, (const xmlChar*)"Name");
      comments[num] = (char*) xmlGetProp(cur_node, 
					 (const xmlChar*)"Comments");
      prop = (char*) xmlGetProp(cur_node, (const xmlChar*)"Timestep");
      stepIDs[num] = atoi(prop);
      free(prop);
      for (i = 0; i < 3; i++) {
	prop = (char*) xmlGetProp(cur_node, (const xmlChar*)MinAttNames[i]);
	lowers[num][i] = atof(prop);
	free(prop);
	prop = (char*) xmlGetProp(cur_node, (const xmlChar*)MaxAttNames[i]);
	uppers[num][i] = atof(prop);
	free(prop);
      }
      num++;
    }
  }

  // done with doc
  xmlFreeDoc(doc);

  // return values
  *numBB = num;
  *bbNames = names;
  *bbComments = comments;
  *bbStepIDs = stepIDs;
  *bbLowers = lowers;
  *bbUppers = uppers;

  return FC_SUCCESS;
}

#endif // HAVE_LIBXML2

/**
 * \ingroup PrivateFileIO
 * \brief Read a .names file
 *
 * \description
 *
 *    Originally for Exodus files, but being coopted for LSDyna files.  File
 *    format is 'block' or 'nodeset' or 'sideset' on a single line. All lines
 *    between these descriptors are of the listed type.
 *
 *    The returned sorted blob arrays hold pointers to \ref _FC_IntStringPair
 *    objects. The string on each object will need to be deleted and each
 *    object will need to be deleted. A helper function 
 *    \ref _fc_freeNamesFileNames() can be used.
 *
 *    It is not an error to not find a names file and the sorted blob
 *    arrays will come back empty.
 *
 * \modifications
 *    - 6/19/2006 WSD Created to replaced duplicated code in exodusio and
 *      lsdynaio modules.
 *    - 8/10/2006 WSD Changed to use sorted blob arrays instead of linked lists.
 */
FC_ReturnCode _fc_readNamesFile(
  char* datasetfilename,     /**< input - path to and name of dataset file */
  FC_SortedBlobArray* block_names,    /**< output- the block names,
				    caller is responsible for freeing internal data */
  FC_SortedBlobArray* nodeset_names,  /**< output - the nodeset names,
				    caller is responsible for freeing internal data */
  FC_SortedBlobArray* sideset_names   /**< output - the sideset names,
				    caller is responsible for freeing internal data */
) 
{
  int rc;
  char* dirname, *basefilename, *namesfilename;
  FILE* namesfile;
  size_t len;
  FC_SortedBlobArray *sba_p;
  _FC_IntStringPair *blob;
  char strbuf[1024];  // FIX!!!! overrun possible!

  // default return - leave blob arrays alone

  // checkinputs
  if (!datasetfilename || !block_names || !nodeset_names || !sideset_names ){
    return FC_INPUT_ERROR;
  }

  // init lists (so will return empties if no file found)
  rc = fc_initSortedBlobArray(block_names);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_initSortedBlobArray(nodeset_names);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_initSortedBlobArray(sideset_names);
  if (rc != FC_SUCCESS)
    return rc;

  // attempt to open the names file
  dirname = fc_getDirname(datasetfilename);
  basefilename = fc_getBasenameWOExtension(datasetfilename, 1);
  len = strlen(dirname) + 1 + strlen(basefilename) + 7;
  namesfilename = (char*)malloc(len*sizeof(char));
  if (!namesfilename) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  if (dirname != NULL)
    sprintf(namesfilename, "%s/%s.names", dirname, basefilename);
  else
    sprintf(namesfilename, "%s.names", basefilename);
  free(dirname);
  free(basefilename);
  namesfile = fopen(namesfilename, "r");
  free(namesfilename);

  // early return - no names file
  if (!namesfile)
    return FC_SUCCESS;

  // scan file
  while(fscanf(namesfile, "%s", strbuf) > 0) {
    if (!strcmp(strbuf, "blocks")) 
      sba_p = block_names; // switch to block list
    else if(!strcmp(strbuf, "sidesets"))
      sba_p = sideset_names; // switch to sidesets list
    else if(!strcmp(strbuf, "nodesets"))
      sba_p = nodeset_names; // switch to nodests list
    else { // add to current list
      blob = (_FC_IntStringPair*)malloc(sizeof(_FC_IntStringPair));
      if (blob == NULL) {
	fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
	return FC_MEMORY_ERROR;
      }
      blob->id = atoi(strbuf);
      fscanf(namesfile, "%[ \t]", strbuf); // discard spaces
      if (fscanf(namesfile, "%[^\n]", strbuf) < 1) { // read rest of line
	fc_printfErrorMessage("Error reading names from '%s'",
			      namesfilename);
	return FC_FILE_IO_ERROR;
      }
      blob->string = (char*)malloc((strlen(strbuf)+1)*sizeof(char));
      if (blob->string == NULL) {
	fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
	return FC_MEMORY_ERROR;
      }
      strcpy(blob->string, strbuf);
      fc_addBlobToSortedBlobArray(sba_p, blob, fc_intCompare);
    }
  }

  // all done
  fclose(namesfile);
  return FC_SUCCESS;
}

/**
 * \ingroup PrivateFileIO
 * \brief Delete the contents returned from _fc_readNamesFile.
 *
 * \description
 *
 *    Deletes the blob data in the sorted blob arrays return from
 *    \ref _fc_readNamesFile. It is assumed that the blobs are
 *    malloc'd \ref _FC_IntStringPair objects with malloc'd strings.
 *
 * \modifications
 *    - 8/10/2006 WSD Created.
 */
void _fc_freeNamesFileNames(
  FC_SortedBlobArray* block_names,    /**< input/output- the block names */
  FC_SortedBlobArray* nodeset_names,  /**< input/output - the nodeset names */
  FC_SortedBlobArray* sideset_names   /**< input/output - the sideset names */
) 
{
  int i, j;
  FC_SortedBlobArray *sba_p;
  _FC_IntStringPair *blob;

  for (i = 0; i < 3; i++) {
    switch(i) {
    case 0: sba_p = block_names;   break;
    case 1: sba_p = nodeset_names; break;
    case 2: sba_p = sideset_names; break;
    }
    if (fc_isSortedBlobArrayValid(sba_p)) {
      for (j = 0; j < sba_p->numBlob; j++) {
	blob = sba_p->blobs[j];
	free(blob->string);
	free(blob);
      }
      fc_freeSortedBlobArray(sba_p);
    }
  }
}

/** \name Queries of File IO enumerated type. */
//-------------------------------
//@{

/**
 * \ingroup   FileIO
 * \brief  Verify that the given value is valid for this enum
 *
 * \modifications   
 *    - 7/14/2005 WSD, Created
 */
int fc_isFileIOTypeValid(
  FC_FileIOType fileType   
) {
  switch(fileType) {
  case FC_FT_NONE:     // fall through
  case FC_FT_EXODUS:   // fall through
  case FC_FT_LSDYNA:   // fall through
  // no default case so will get warning if new values added to enum
    return 1;
  }
  return 0;
}

/**
 * \ingroup   FileIO
 * \brief  Return name of the FC_FileIOType's value.
 *
 * \modifications   
 *    - 7/14/2005 WSD, Created
 */
char *fc_getFileIOTypeText(
  FC_FileIOType fileType   
) {
  switch(fileType) {
  case FC_FT_NONE:    return "FC_FT_NONE";
  case FC_FT_EXODUS:  return "FC_FT_EXODUS";
  case FC_FT_LSDYNA:  return "FC_FT_LSDYNA";
  // no default case so will get warning if new values added to enum
  }
  return "invalid value for enum FC_FileIOType";
}

/**
 * \ingroup   FileIO
 * \brief  Verify that a file IO type can be read.
 *
 * \description
 *
 *    The library can be build with or without support for various file types
 *    and some types are supported only for reading.  Use this routine to
 *    determine if the requested file type is supported for reading.
 *    
 *    Note that a file IO type of FC_FT_NONE is considered supported
 *    since it requires no support.
 *
 * \todo Think more about what FC_FT_NONE should return.
 *
 * \modifications   
 *    - 2/1/2006 WSD, Created
 */
int fc_isFileIOTypeReadSupported(
  FC_FileIOType fileType   
) {
  switch(fileType) {
  case FC_FT_NONE:     
    return 1;  // requires no support, so supported by default
  case FC_FT_EXODUS:  
#ifdef HAVE_EXODUS
    return 1;
#else
    return 0;
#endif
  case FC_FT_LSDYNA:
    return 1;
  // no default case so will get warning if new values added to enum
  }
  return 0;
}

/**
 * \ingroup   FileIO
 * \brief  Verify that a file IO type can be written.
 *
 * \description
 *
 *    The library can be build with or without support for various file types
 *    and some types are supported only for reading.  Use this routine to
 *    determine if the requested file type is supported for writing.
 *    
 *    Note that a file IO type of FC_FT_NONE is considered supported
 *    since it requires no support.
 *
 * \todo Think more about what FC_FT_NONE should return.
 *
 * \modifications   
 *    - 2/1/2006 WSD, Created
 */
int fc_isFileIOTypeWriteSupported(
  FC_FileIOType fileType   
) {
  switch(fileType) {
  case FC_FT_NONE:     
    return 0;  // can't do, so not supported
  case FC_FT_EXODUS:  
#ifdef HAVE_EXODUS
    return 1;
#else
    return 0;
#endif
  case FC_FT_LSDYNA:
    return 0;
  // no default case so will get warning if new values added to enum
  }
  return 0;
}

//@}

/** \name File IO operations on datasets. */
//-------------------------------
//@{

/** 
 * \ingroup FileIO
 * \brief Get type of file format.
 *
 * \modifications
 *    - 7/14/2005 WSD Created.
 */
FC_ReturnCode fc_getDatasetFileIOType(
  FC_Dataset dataset,      /**< input - dataset handle */
  FC_FileIOType* fileType  /**< output - file format type */
) {
  _FC_DsSlot* dsSlot;

  // default return value
  if (fileType)
    *fileType = FC_FT_NONE;

  // check input
  dsSlot = _fc_getDsSlot(dataset);
  if (dsSlot == NULL || fileType == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Getting file type for dataset '%s'", 
		      dsSlot->header.name);
  
  // do it
  *fileType = dsSlot->fileType;
  return FC_SUCCESS;
}

/**
 * \ingroup  FileIO
 * \brief  Load a dataset and it's members from a file.
 *
 * \description 
 *
 *    Given a the name of a database file, this routine loads the dataset and
 *    makes all of it's members available for use.
 *
 *    This routine automatically picks file format based on the name of the
 *    file. If this is unsatisfactory, use fc_loadDatasetWithFormat.  See
 *    FC_FileIOType for supported formats.
 *
 *    For some formats the members are loaded lazily: that is all of
 *    their metadata data is loaded immediately, but the big data arrays are
 *    not loaded until needed (this happens automatically). (this was
 *    true with saf and is only partially true with exodus)
 *
 * \modifications
 *   - 7/13/2004 WSD. Old method contents moved to fc_loadDatasetWithFormat.
 *     This is a convenience function w/o the fileType arg. 
 */
FC_ReturnCode fc_loadDataset(
 char *filename,         /**< input - name of dataset to open */
 FC_Dataset *dataset     /**< output - handle for this dataset */
) {
  return fc_loadDatasetWithFormat(filename, FC_FT_NONE, dataset);
}

/**
 * \ingroup  FileIO
 * \brief  Load a dataset and it's members from a file.
 *
 * \description 
 *
 *    Given a the name of a database file, this routine loads the dataset and
 *    makes all of it's members available for use.
 *
 *    If a specific filetype is passed in, the load will fail if the file is
 *    not of that type. If FC_FT_NONE is passed in, this routine automatically
 *    determines the dataset type.  See FC_FileIOType for supported formats.
 *
 *    For some formats, the members are loaded lazily: that is all of
 *    their metadata data is loaded immediately, but the big data arrays are
 *    not loaded until needed (this happens automatically). (This was
 *    true with saf, and is ony partially true with exodus)
 *
 * \modifications
 *   - Nancy Collins  Created.
 *   - Aug 5, 2002  W Koegler  Now calls _fc_setupDsSlot() so that 
 *       fc_openDataset() and  use_dataset() can use the same code to setup 
 *       the ftable slot. (RM comment: use_dataset() does not exist anymore.)
 *   - 11/20/03 RM changed behavior: when library is not initialized - 
 *       return an error.
 *   - 11/25/03 RM, always return default NULL dataset handle if anything
 *       goes wrong.
 *   - 09/07/04 WSK, changed name from fc_openDataset to fc_loadDataset.
 *   - 01/13/05 WSK, this code now wraps calls to database type specific
 *       codes. Most of the old contents have been moved to
 *       _fc_loadSAFDataset(). 
 *   - 7/13/02005 WSD Moved to fileio.c.
 *   - 7/13/2004 WSD Added flag for file format type. Name changed.
 *   - 4/17/2007 CDU Added check to see if file is a regular file
 */
FC_ReturnCode fc_loadDatasetWithFormat(
 char *filename,         /**< input - name of dataset to open */
 FC_FileIOType fileType, /**< input - File Type (pass FC_FT_NONE if unknonw */
 FC_Dataset *dataset     /**< output - handle for this dataset */
) {
  struct stat f_stat;
  int res;
  FC_ReturnCode rc;
  
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
    fc_printfErrorMessage("Cannot load a dataset until library is initialized.");
    return FC_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Loading dataset '%s'", filename);

  // Make sure this is a readable file in case lower io libs 
  // don't check (ie, exodus)
  res=stat(filename, &f_stat);
  if(res || !S_ISREG(f_stat.st_mode)){
    fc_printfErrorMessage("File %s is not a regular file.", filename);
    return FC_ERROR;
  }

  // Pick a type 
  // FIX make this smarter?
  if (fileType == FC_FT_NONE) {
    char *basename = fc_getBasename(filename);
    char *extension = fc_getExtension(filename, 1);
    if (!strcasecmp(basename, "d3plot")) // test this first
      fileType = FC_FT_LSDYNA;
    else if (extension) {
      if (!strcasecmp(extension, "e") || !strcasecmp(extension, "e2") ||
	  !strcasecmp(extension, "exo"))
	fileType = FC_FT_EXODUS;
    }
    free(basename);
    free(extension);
  }

  switch (fileType) {
  case FC_FT_EXODUS:
#ifdef HAVE_EXODUS
    return _fc_loadExodusDataset(filename, dataset);
#else
    fc_printfErrorMessage("Exodus module not available");
    return FC_FILE_IO_ERROR;
#endif
  case FC_FT_LSDYNA:
    return _fc_loadLSDynaDataset(filename, dataset);
  case FC_FT_NONE: // fall through
  default:
    // Try every type until we get one that works
#ifdef HAVE_EXODUS
    rc = _fc_loadExodusDataset(filename, dataset);
    if (rc == FC_SUCCESS)
      return rc;
#endif
    rc = _fc_loadLSDynaDataset(filename, dataset);
    if (rc == FC_SUCCESS)
      return rc;
    // Couldn't find one 
    fc_printfErrorMessage("Could not open database '%s'", 
			  filename);
    return FC_FILE_IO_ERROR;
  }
}

/**
 * \ingroup  FileIO
 * \brief Write a dataset
 *
 * \description  
 *
 *    Currently, only a created dataset can be written; the function cannot
 *    update changes to a dataset read from disk.  To save changes to
 *    a dataset that has been read to disk, or already written, use
 *    fc_rewriteDataset(). 
 *
 *    Empty subsets do not get written.
 *
 * \modifications  
 *   - 2003-NOV-11 WSK Created.
 *   - 7/13/2005 WSD Moved to fileio.c.
 *   - 7/13/2005 WSD Made into a wrapper (old guts are now _fc_writeSAFDataset)
 *   - 10/15/2007 ACG exodus is now default file type as saf is removed
 */
FC_ReturnCode fc_writeDataset(
  FC_Dataset dataset,      /**< the dataset to write */
  char* filename,          /**< the filename to write to */
  FC_FileIOType fileType   /**< the file type to write */
) { 
  _FC_DsSlot* dsSlot;

  // check input
  dsSlot = _fc_getDsSlot(dataset);
  if (dsSlot == NULL || dsSlot->writable == 0 || 
      (dsSlot->numSeq == 0 && dsSlot->numMesh == 0) ||
      !fc_isFileIOTypeValid(fileType)) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // log message
  fc_printfLogMessage("Writing dataset '%s'", dsSlot->header.name);

  // if no type, choose one
  if (fileType == FC_FT_NONE)
    fileType = FC_FT_EXODUS;

  // Call format specific writers
  switch (fileType) {
  case FC_FT_EXODUS:
#ifdef HAVE_EXODUS
    return _fc_writeExodusDataset(dataset, filename);
#else
    fc_printfErrorMessage("Exodus module not available");
    return FC_FILE_IO_ERROR;
#endif
  case FC_FT_LSDYNA: 
    fc_printfErrorMessage("This installation does not support writing files "
			  "in LS-Dyna format\n");
    return FC_ERROR;
  case FC_FT_NONE: // fall through
    fc_printfErrorMessage("Should never reach this point! should have "
			  "chosen a default format");
    return FC_ERROR;
  }

  // error if reach this point
  fc_printfErrorMessage("Programmer error!: Should never reach this point\n");
  return FC_ERROR;
}

/**
 * \ingroup  FileIO
 * \brief Rewrite a dataset.
 *
 * \description  
 *
 *    This will force the writing of a dataset that is already associated with
 *    a file on disk, e.g. to rewrite an Exodus dataset or write an Exodus
 *    dataset into a differnt type of file. Calling on a dataset that has
 *    never been associated with a file is o.k.
 *
 *    Empty subsets do not get written.
 *
 * \todo ?Passing NULL for the name will make it use the original
 *    dataset's filename and type -- essentially rewriting it?
 * \todo ?Replace fc_writeDataset with this version?
 * \todo Write to filename.tmp then move filename.tmp to filename --
 *    then if rewriting a file, old version less likely to get
 *    clobbered in case of an error.
 *
 * \modifications  
 *   - 3/25/2006 WSD Created.
 */
FC_ReturnCode fc_rewriteDataset(
  FC_Dataset dataset,      /**< the dataset to write */
  char* filename,          /**< the filename to write to */
  FC_FileIOType fileType   /**< the file type to write */
) { 
  FC_ReturnCode rc;
  _FC_DsSlot* dsSlot;

  // check input
  dsSlot = _fc_getDsSlot(dataset);
  if (dsSlot == NULL || (dsSlot->numSeq == 0 && dsSlot->numMesh == 0) ||
      !fc_isFileIOTypeValid(fileType)) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // if not writable, make writable
  if (dsSlot->writable == 0) {
    int i, j, k;
    int numSeq, numGlbVar, numGlbSeqVar, numMesh;
    int numSub, numVar, numSeqVar, numMember, *numSteps;
    FC_Sequence *sequences;
    _FC_SeqSlot* seqSlot;
    FC_Mesh* meshes;
    _FC_MeshSlot* meshSlot;
    FC_Subset* subsets;
    _FC_SubSlot* subSlot;
    FC_Variable *glbVars, **glbSeqVars, *variables, **seqVars;
    _FC_VarSlot* varSlot;
    double* coords;
    int *conns, *ids;
    void* data;

    // Force loading of all data & set to uncommited
    fc_getGlobalVariables(dataset, &numGlbVar, &glbVars);
    for (i = 0; i < numGlbVar; i++) {
      fc_getVariableDataPtr(glbVars[i], &data);
      varSlot = _fc_getVarSlot(glbVars[i]);
      varSlot->header.committed = 0;
    }
    free(glbVars);
    fc_getGlobalSeqVariables(dataset, &numGlbSeqVar, &numSteps, &glbSeqVars);
    for (i = 0; i < numGlbSeqVar; i++) {
      for (j = 0; j < numSteps[i]; j++) {
        fc_getVariableDataPtr(glbSeqVars[i][j], &data);
        varSlot = _fc_getVarSlot(glbSeqVars[i][j]);
        varSlot->header.committed = 0;
      }
      free(glbSeqVars[i]);
    }
    free(numSteps);
    free(glbSeqVars);
    fc_getSequences(dataset, &numSeq, &sequences);
    for (i = 0; i < numSeq; i++) {
      fc_getSequenceCoordsPtr(sequences[i], &data);
      seqSlot = _fc_getSeqSlot(sequences[i]);
      seqSlot->header.committed = 0;
    }
    fc_getMeshes(dataset, &numMesh, &meshes);
    for (i = 0; i < numMesh; i++) {
      fc_getMeshCoordsPtr(meshes[i], &coords);
      fc_getMeshElementConnsPtr(meshes[i], &conns);
      fc_getSubsets(meshes[i], &numSub, &subsets);
      for (j = 0; j < numSub; j++) {
	fc_getSubsetMembersAsArray(subsets[j], &numMember, &ids);
	free(ids);
	subSlot = _fc_getSubSlot(subsets[j]);
	subSlot->header.committed = 0;
      }
      free(subsets);
      fc_getVariables(meshes[i], &numVar, &variables);
      for (j = 0; j < numVar; j++) {
	fc_getVariableDataPtr(variables[j], &data);
	varSlot = _fc_getVarSlot(variables[j]);
	varSlot->header.committed = 0;
      }
      free(variables);
      fc_getSeqVariables(meshes[i], &numSeqVar, &numSteps, &seqVars);
      for (j = 0; j < numSeqVar; j++) {
	for (k = 0; k < numSteps[j]; k++) { 
	  fc_getVariableDataPtr(seqVars[j][k], &data);
	  varSlot = _fc_getVarSlot(seqVars[j][k]);
	  varSlot->header.committed = 0;
	}
	free(seqVars[j]);
      }
      free(numSteps);
      free(seqVars);
      meshSlot = _fc_getMeshSlot(meshes[i]);
      meshSlot->header.committed = 0;
    }
    free(sequences);
    free(meshes);
    dsSlot->header.committed = 0;

    // reset to look like a created dataset
    rc = _fc_closeDataset(dataset);
    if (rc != FC_SUCCESS)
      return rc;
    dsSlot->writable = 1;
  }
 
  return fc_writeDataset(dataset, filename, fileType);
}

//@}

/**
 * \ingroup PrivateFileIO
 * \brief Initial for File IO
 *
 * \description
 *  
 *    This sets up stuff for FILE IO
 *
 * \todo ? could have flags so that specific type of file IO is not initialized
 *   until use--but there could possibly be stuff that you'd want to do up 
 *   front ....
 *
 * \modifications  
 *   - 7/19/2005 WSD, created to better abstract reading routines.
 */
FC_ReturnCode _fc_initFileIO(void) {
  FC_ReturnCode rc, rc_keep = FC_SUCCESS;

#ifdef HAVE_LIBXML2
  LIBXML_TEST_VERSION
#endif

#ifdef HAVE_EXODUS
  rc = _fc_initExodus();
  if (rc != FC_SUCCESS) {
    fc_printfErrorMessage("Unable to init Exodus File IO");
    rc_keep = rc;
  }
#endif
  // nothing for LSDyna
  return rc_keep;
}

/**
 * \ingroup PrivateFileIO
 * \brief Finalize for File IO
 *
 * \description
 *  
 *    This finalizes stuff for FILE IO
 *
 * \todo? could have flags so that specific type of file IO is only finalize
 *   if it was inited ... see Fix on _fc_initFileIO().
 *
 * \modifications  
 *   - 7/19/2005 WSD, created to better abstract reading routines.
 */
FC_ReturnCode _fc_finalFileIO(void) {
  FC_ReturnCode rc, rc_keep = FC_SUCCESS;
  rc = FC_SUCCESS; // so the compiler won't complaid

#ifdef HAVE_LIBXML2
  xmlCleanupParser();
#endif

  // no Exodus
  // no LSDyna

  return rc_keep;
}

/**
 * \ingroup  PrivateFileIO
 * \brief  Read coords from file into a sequence.
 *
 * \description
 *  
 *    Reads step coordinates from file into a sequence. This is a wrapper
 *    that chooses the proper implentation based on the file type.
 *
 * \modifications  
 *   - 7/13/2005 WSD, created to better abstract reading routines.
 */
FC_ReturnCode _fc_readSequenceCoords(
  FC_Sequence sequence  /**< input - sequence to read into */
) {
  _FC_SeqSlot* seqSlot = _fc_getSeqSlot(sequence);
  _FC_DsSlot* dsSlot;

  // check input
  if (!seqSlot || !seqSlot->header.name) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // log message
  fc_printfLogMessage("Reading coords for sequence '%s'", 
		      seqSlot->header.name);

  // Call format specific routines
  dsSlot = _fc_getDsSlot(seqSlot->ds);
  switch(dsSlot->fileType) {
  case FC_FT_EXODUS:
#ifdef HAVE_EXODUS
    return _fc_readExodusSequenceCoords(sequence);
#else
    fc_printfErrorMessage("Exodus module not available");
    return FC_ERROR;
#endif
  case FC_FT_LSDYNA:
    return FC_SUCCESS; // do nothing, no lazy capability
  case FC_FT_NONE:
    fc_printfErrorMessage("Called without a file type");
    return FC_ERROR;
  }

  // shouldn't reach this point
  fc_printfErrorMessage("Unrecognized file type");
  return FC_ERROR;
}

/**
 * \ingroup  PrivateFileIO
 * \brief  Read coords from file into a mesh
 *
 * \description
 *  
 *    Reads vertex coordinates from file into a mesh. This is a wrapper
 *    that chooses the proper implentation based on the file type.
 *
 * \modifications  
 *   - 7/13/2005 WSD, created to better abstract reading routines.
 */
FC_ReturnCode _fc_readMeshCoords(
  FC_Mesh mesh  /**< input - mesh to read into */
) {
  _FC_MeshSlot* meshSlot;
  _FC_DsSlot* dsSlot;

  // check input
  meshSlot = _fc_getMeshSlot(mesh);
  if (!meshSlot || !meshSlot->header.name) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // log message
  fc_printfLogMessage("Reading coords for mesh '%s'", 
		      meshSlot->header.name);

  // Call format specific routines
  dsSlot = _fc_getDsSlot(meshSlot->ds);
  switch (dsSlot->fileType) {
  case FC_FT_EXODUS:
#ifdef HAVE_EXODUS
    return _fc_readExodusMeshCoords(mesh);
#else
    fc_printfErrorMessage("Exodus module not available");
    return FC_ERROR;
#endif
  case FC_FT_LSDYNA: 
    return _fc_readLSDynaMeshCoords(mesh);
  case FC_FT_NONE:
    fc_printfErrorMessage("Called without a file type");
    return FC_ERROR;
  }

  // shouldn't reach this point
  fc_printfErrorMessage("Unrecognized file type");
  return FC_ERROR;
}

/**
 * \ingroup  PrivateFileIO
 * \brief  Read element conns from file into a mesh
 *
 * \description
 *  
 *    Reads element connectivities from file into a mesh. This is a wrapper
 *    that chooses the proper implentation based on the file type.
 *
 * \modifications  
 *   - 7/13/2005 WSD, created to better abstract reading routines.
 */
FC_ReturnCode _fc_readMeshElemConns(
  FC_Mesh mesh  /**< input - mesh to read into */
) {
  _FC_MeshSlot* meshSlot;
  _FC_DsSlot* dsSlot;

  // check input
  meshSlot = _fc_getMeshSlot(mesh);
  if (!meshSlot || !meshSlot->header.name) {
    fc_printfErrorMessage("INPUT ERROR: invalid mesh");
    return FC_INPUT_ERROR;
  }

  // log message
  fc_printfLogMessage("Reading element conns for mesh '%s'", 
		      meshSlot->header.name);

  // Call format specific routines
  dsSlot = _fc_getDsSlot(meshSlot->ds);
  switch (dsSlot->fileType) {
  case FC_FT_EXODUS:
#ifdef HAVE_EXODUS
   return _fc_readExodusMeshElemConns(mesh);
#else
    fc_printfErrorMessage("Exodus module not available");
    return FC_ERROR;
#endif
  case FC_FT_LSDYNA: 
    return _fc_readLSDynaMeshElemConns(mesh);
  case FC_FT_NONE:
    fc_printfErrorMessage("Called without a file type");
    return FC_ERROR;
  }

  // shouldn't reach this point
  fc_printfErrorMessage("Unrecognized file type");
  return FC_ERROR;
}

/**
 * \ingroup  PrivateFileIO
 * \brief  Read subset members from file into a subset
 *
 * \description
 *  
 *    Reads subset members from file into a subset. This is a wrapper that
 *     chooses the proper implentation based on the file type.
 *
 * \modifications  
 *   - 7/13/2005 WSD, created to better abstract reading routines.
 */
FC_ReturnCode _fc_readSubsetMembers(
  FC_Subset subset  /**< input - subset to read into */
) {
  _FC_SubSlot* subSlot = _fc_getSubSlot(subset);
  _FC_DsSlot* dsSlot;

  // check input
  if (!subSlot || !subSlot->header.name) {
    fc_printfErrorMessage("INPUT ERROR: invalid subSlot");
    return FC_INPUT_ERROR;
  }

  // log message
  fc_printfLogMessage("Reading members for subset '%s'", 
		      subSlot->header.name);

  // Call format specific routine
  dsSlot = _fc_getDsSlot(_fc_getMeshSlot(subSlot->mesh)->ds);
  switch (dsSlot->fileType) {
  case FC_FT_EXODUS:
#ifdef HAVE_EXODUS
    return _fc_readExodusSubsetMembers(subset);
#else
    fc_printfErrorMessage("Exodus module not available");
    return FC_ERROR;
#endif
  case FC_FT_LSDYNA:
    return FC_SUCCESS; // no subsets in lsdyna
  case FC_FT_NONE:
    fc_printfErrorMessage("Called without a file type");
    return FC_ERROR;
  }

  // shouldn't reach this point
  fc_printfErrorMessage("Unrecognized file type");
  return FC_ERROR;
}

/**
 * \ingroup  PrivateFileIO
 * \brief  Read data from file into a variable
 *
 * \description
 *  
 *    Reads data from file into a variable. This is a wrapper that chooses the
 *    proper implentation based on the file type.
 *
 * \modifications  
 *   - 7/13/2005 WSD, created to better abstract reading routines.
 */
FC_ReturnCode _fc_readVariableData(
  FC_Variable variable  /**< input - variable to read into */
) {
  _FC_VarSlot* varSlot = _fc_getVarSlot(variable);
  _FC_DsSlot* dsSlot;
  
  // check input
  if (!varSlot || !varSlot->header.name) {
    fc_printfErrorMessage("INPUT ERROR: invalid variable");
    return FC_INPUT_ERROR;
  }

  // log message
  fc_printfLogMessage("Reading data for variable '%s'", 
		      varSlot->header.name);

  // Call format specific routines
  dsSlot = _fc_getDsSlot(varSlot->ds);
  switch (dsSlot->fileType) {
  case FC_FT_EXODUS:
#ifdef HAVE_EXODUS
    return _fc_readExodusVariableData(variable);
#else
    fc_printfErrorMessage("Exodus module not available");
    return FC_ERROR;
#endif
  case FC_FT_LSDYNA: 
    return _fc_readLSDynaVariableData(variable);
  case FC_FT_NONE:
    fc_printfErrorMessage("Called without a file type");
    return FC_ERROR;
  }

  // shouldn't reach this point
  fc_printfErrorMessage("Unrecognized file type");
  return FC_ERROR;
}

/**
 * \ingroup  PrivateFileIO
 * \brief  Close a dataset.
 *
 * \description 
 * 
 *    Close the dataset file (if needed) and do any necessary file related
 *    cleanup for the Dataset. This is a wrapper that calls format specific
 *    implementations.
 *
 * \modifications  
 *   - 7/14/2005 WSD, created to better abstract file IO routines.
 */
FC_ReturnCode _fc_closeDataset(
  FC_Dataset dataset /**< input - the dataset */
) {
  FC_ReturnCode rc;
  _FC_DsSlot* dsSlot;

  // check input
  dsSlot = _fc_getDsSlot(dataset);
  if (dsSlot == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // Call format specific routines
  dsSlot = _fc_getDsSlot(dataset);
  switch (dsSlot->fileType) {
  case FC_FT_EXODUS:
#ifdef HAVE_EXODUS
    rc =_fc_closeExodusDataset(dataset); break;
#else
    fc_printfErrorMessage("Exodus module not available");
    return FC_ERROR;
#endif
  case FC_FT_LSDYNA: 
    rc = _fc_closeLSDynaDataset(dataset); break;
  case FC_FT_NONE:
    fc_printfErrorMessage("Called without a file type");
    return FC_ERROR;
  }
  if (rc != FC_SUCCESS)
    return rc;

  // Generic stuff
  dsSlot->fileType = FC_FT_NONE;
  free(dsSlot->baseFile);
  dsSlot->baseFile = NULL;

  return FC_SUCCESS;
}

