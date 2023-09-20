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
 * \file checkfileio.c
 * \brief Unit tests for \ref FileIO Module.
 * 
 * $Source: /usr/local/Repositories/fcdmf/fclib/unittest/checkfileio.c,v $
 * $Revision: 1.45 $
 * $Date: 2006/11/28 19:09:15 $
 *
 * \modifications
 *    - 08/01/2005 WSD. Created.
 *    - 10/15/2007 ACG removed saf
 */

#include <stdlib.h>
#include <string.h>
#include <check.h>
#include "fc.h"
#include "fcP.h"
#include "checkall.h"

#define BAD_FILE_NAME "sdfssf"
#define WRONG_FILE_NAME "checkfileio.c"
#define EXODUS_FILE_NAME "temp.ex2"

// *** global vars 

static int numMeshType = 8;
static char* vertFiles[8] = { "../data/gen_point_line_tri_quad.vert",
			      "../data/gen_point_line_tri_quad.vert",
			      "../data/gen_point_line_tri_quad.vert",
			      "../data/gen_point_line_tri_quad.vert",
			      "../data/gen_tet_prism_hex.vert",
			      "../data/gen_pyramid.vert",
			      "../data/gen_tet_prism_hex.vert",
			      "../data/gen_tet_prism_hex.vert" };
static char* elemFiles[8] = { "../data/gen_point.elem",
			      "../data/gen_line.elem",
			      "../data/gen_tri.elem",
			      "../data/gen_quad.elem",
			      "../data/gen_tet.elem",
			      "../data/gen_pyramid.elem",
			      "../data/gen_prism.elem",
			      "../data/gen_hex.elem" };
static int topoDims[8] = { 0, 1, 2, 2, 3, 3, 3, 3 };

static int numAssocType = 5;
static FC_AssociationType assocTypes[5] = { FC_AT_VERTEX, 
					    FC_AT_EDGE, 
					    FC_AT_FACE,
					    FC_AT_ELEMENT, 
					    FC_AT_WHOLE_MESH };

// ***** Fixtures for the rest of the file io tests

static void fileio_setup(void) 
{
  FC_ReturnCode rc;
  
  if (isForking) {
    fc_setLibraryVerbosity(fc_messages);
    rc = fc_initLibrary();
    fail_unless(rc == FC_SUCCESS, "test aborted: failed to init library");
  }
}

static void fileio_teardown(void) 
{
  FC_ReturnCode rc;
  
  if (isForking) {
    rc = fc_finalLibrary();
    fail_unless(rc == FC_SUCCESS, "failed to final library at end of test");
  }
}

// **** test bootstrap - that is, read generated data files

START_TEST(vert_files)
{
  FC_ReturnCode rc;
  int i;
  int temp_numVert, temp_numDim;
  double* temp_coords;
  int numFile = 1;
  char* filenames[1] = { "../data/gen_point_line_tri_quad.vert" };
  int numVerts[1] = { 156 };
  int numDims[1] = { 3 };
  double startCoords[1] = { 0 };
  double endCoords[1] = { 2.3925586960645e-05 };

  // Should fail if file doesn't exist
  rc = _fc_readVertexFile("nonsense", &temp_numVert, &temp_numDim, 
			  &temp_coords);
  fail_unless(rc != FC_SUCCESS, "should fail to open file");
  fail_unless(temp_numVert == -1 && temp_numDim == -1 && temp_coords == NULL, 
	      "fail should return null values");

  // Should fail if bad args
  rc = _fc_readVertexFile(NULL, &temp_numVert, &temp_numDim, &temp_coords);
  fail_unless(rc != FC_SUCCESS, "should fail if file is NULL");
  fail_unless(temp_numVert == -1 && temp_numDim == -1 && temp_coords == NULL, 
	      "fail should return null values");
  rc = _fc_readVertexFile(filenames[0], NULL, &temp_numDim, &temp_coords);
  fail_unless(rc != FC_SUCCESS, "should fail if numvert is NULL");
  fail_unless(temp_numDim == -1 && temp_coords == NULL, 
	      "fail should return null values");
  rc = _fc_readVertexFile(filenames[0], &temp_numVert, NULL, &temp_coords);
  fail_unless(rc != FC_SUCCESS, "should fail if numDim is NULL");
  fail_unless(temp_numVert == -1 && temp_coords == NULL, 
	      "fail should return null values");
  rc = _fc_readVertexFile(filenames[0], &temp_numVert, &temp_numDim, NULL);
  fail_unless(rc != FC_SUCCESS, "should fail if coords is NULL");
  fail_unless(temp_numVert == -1 && temp_numDim == -1, 
	      "fail should return null values");
  
  // Test all types of vert files
  for (i = 0; i < numFile; i++) {
    rc = _fc_readVertexFile(filenames[i], &temp_numVert, &temp_numDim, 
			    &temp_coords);
    fail_unless(rc == FC_SUCCESS, "should not fail to read file");
    fail_unless(temp_numVert == numVerts[i], "mismatch of numVert");
    fail_unless(temp_numDim == numDims[i], "mismatch of numDim");
    fail_unless(temp_coords[0] == startCoords[i], "mismatch of 1st coord");
    fail_unless(temp_coords[temp_numDim*temp_numVert-1] == endCoords[i],
		"mismatch of final coord");
    free(temp_coords);
  }
}
END_TEST

START_TEST(elem_files)
{
  FC_ReturnCode rc;
  int i;
  int temp_numVert, temp_numElem;
  FC_ElementType temp_elemType;
  char* temp_name;
  int* temp_conns;
  int numFile = 8;
  char* filenames[8] = { "../data/gen_point.elem", "../data/gen_line.elem",
			 "../data/gen_tri.elem", "../data/gen_quad.elem",
			 "../data/gen_tet.elem", "../data/gen_pyramid.elem",
			 "../data/gen_prism.elem", "../data/gen_hex.elem" };
  int numVerts[8] = { 156, 156, 156, 156, 1092, 1884, 1092, 1092 };
  int numElems[8] = { 156, 396, 264, 132, 3960, 4752, 1584, 792 };
  FC_ElementType elemTypes[8] = { FC_ET_POINT, FC_ET_LINE, FC_ET_TRI,
				  FC_ET_QUAD, FC_ET_TET, FC_ET_PYRAMID,
				  FC_ET_PRISM, FC_ET_HEX };
  char* names[8] = { "point mesh", "line mesh", "tri mesh", "quad mesh",
		     "tet mesh", "pyramid mesh", "prism mesh", "hex mesh" };
  int startConns[8] = { 0, 0, 0, 0, 0, 0, 0, 0 };
  int endConns[8] = { 155, 142, 154, 154, 1090, 1883, 1090, 1090 };
  int numVperE;

  // Should fail if file doesn't exist
  rc = _fc_readElementFile("nonsense", &temp_numVert, &temp_numElem, 
			   &temp_elemType, &temp_name, &temp_conns);
  fail_unless(rc != FC_SUCCESS, "should fail to open file");
  fail_unless(temp_numVert == -1 && temp_numElem == -1 && 
	      temp_elemType == FC_ET_UNKNOWN && temp_name == NULL &&
	      temp_conns == NULL, "fail should return null values");

  // Should fail if bad args
  rc = _fc_readElementFile(NULL, &temp_numVert, &temp_numElem,
			   &temp_elemType, &temp_name, &temp_conns);
  fail_unless(rc != FC_SUCCESS, "should fail if file is NULL");
  fail_unless(temp_numVert == -1 && temp_numElem == -1 && 
	      temp_elemType == FC_ET_UNKNOWN && temp_name == NULL &&
	      temp_conns == NULL, "fail should return null values");
  rc = _fc_readElementFile(filenames[0], NULL, &temp_numElem,
			   &temp_elemType, &temp_name, &temp_conns);
  fail_unless(rc != FC_SUCCESS, "should fail if numVert is NULL");
  fail_unless(temp_numElem == -1 && 
	      temp_elemType == FC_ET_UNKNOWN && temp_name == NULL &&
	      temp_conns == NULL, "fail should return null values");
  rc = _fc_readElementFile(filenames[0], &temp_numVert, NULL,
			   &temp_elemType, &temp_name, &temp_conns);
  fail_unless(rc != FC_SUCCESS, "should fail if numElem is NULL");
  fail_unless(temp_numVert == -1 &&  
	      temp_elemType == FC_ET_UNKNOWN && temp_name == NULL &&
	      temp_conns == NULL, "fail should return null values");
  rc = _fc_readElementFile(filenames[0], &temp_numVert, &temp_numElem,
			   NULL, &temp_name, &temp_conns);
  fail_unless(rc != FC_SUCCESS, "should fail if elemType is NULL");
  fail_unless(temp_numVert == -1 && temp_numElem == -1 && 
	      temp_name == NULL &&
	      temp_conns == NULL, "fail should return null values");
  rc = _fc_readElementFile(filenames[0], &temp_numVert, &temp_numElem,
			   &temp_elemType, NULL, &temp_conns);
  fail_unless(rc != FC_SUCCESS, "should fail if name is NULL");
  fail_unless(temp_numVert == -1 && temp_numElem == -1 && 
	      temp_elemType == FC_ET_UNKNOWN &&
	      temp_conns == NULL, "fail should return null values");
  rc = _fc_readElementFile(filenames[0], &temp_numVert, &temp_numElem,
			   &temp_elemType, &temp_name, NULL);
  fail_unless(rc != FC_SUCCESS, "should fail if conns is NULL");
  fail_unless(temp_numVert == -1 && temp_numElem == -1 && 
	      temp_elemType == FC_ET_UNKNOWN && temp_name == NULL,
	      "fail should return null values");

  // Test all types of elems files
  for (i = 0; i < numFile; i++) {
    rc = _fc_readElementFile(filenames[i], &temp_numVert, &temp_numElem,
			     &temp_elemType, &temp_name, &temp_conns);
    fail_unless(rc == FC_SUCCESS, "should not fail to read file");
    fail_unless(temp_numVert == numVerts[i], "mismatch of numVert");
    fail_unless(temp_numElem == numElems[i], "mismatch of numElem");
    fail_unless(temp_elemType == elemTypes[i], "mismatch of elemType");
    fail_unless(!strcmp(temp_name, names[i]), "mismatch of name");
    fail_unless(temp_conns[0] == startConns[i], "mismatch of 1st conn");
    numVperE = fc_getElementTypeNumVertex(temp_elemType);
    fail_unless(temp_conns[temp_numElem*numVperE-1] == endConns[i],
		"mismatch of final conn");
    free(temp_name);
    free(temp_conns);
  }
}
END_TEST

START_TEST(data_files)
{
  FC_ReturnCode rc;
  int i, j;
  int temp_numPoint, temp_numComp;
  FC_AssociationType temp_assoc;
  FC_MathType temp_mathType;
  char* temp_name;
  float* temp_values;
  int numFile = 5;
  char* filenames[5] = { "../data/gen_point_line_tri_quad.data", 
			 "../data/gen_tri_perElem.data",
			 "../data/gen_small_tri_displ.data",
			 "../data/gen_small_tri_vstress.data",
			 "../data/gen_small_tri_estress.data" };
  int numPoints[5] = { 156, 264, 200, 200, 383 };
  int numComps[5] = { 1, 1, 3, 3, 3 };
  FC_AssociationType assocs[5] = { FC_AT_VERTEX, FC_AT_ELEMENT,
				   FC_AT_VERTEX, FC_AT_VERTEX,
				   FC_AT_ELEMENT };
  FC_MathType mathTypes[5] = { FC_MT_SCALAR, FC_MT_SCALAR,
			       FC_MT_VECTOR, FC_MT_SYMTENSOR,
			       FC_MT_SYMTENSOR };
  char* names[5] = { "temperature", "temperature", "displacement",
		     "stress", "stress" };
  float startValues[5][3] = { { 60.033909 }, 
			      { 60.067757 },
			      { 0.121143, 0.51015, 1.08402 },
			      { 0.121143, 0.51015, 1.08402 },
			      { 0.101315, 0.775898, 1.039923 } };
  float endValues[5][3] = { { 60.033909 }, 
			    { 60.114097 },
			    { 0.749443, 0.270375, 1.0831 }, 
			    { 0.749443, 0.270375, 1.0831 }, 
			    { 0.671951, 0.313660, 1.048657 } };

  // Should fail if file doesn't exist
  rc = _fc_readDataFile("nonsense", &temp_numPoint, &temp_numComp, &temp_assoc,
			&temp_mathType, &temp_name, &temp_values);
  fail_unless(rc != FC_SUCCESS, "should fail to open file");
  fail_unless(temp_numPoint == -1 && temp_numComp == -1 &&
	      temp_assoc == FC_AT_UNKNOWN && temp_mathType == FC_ET_UNKNOWN &&
	      temp_name == NULL && temp_values == NULL, 
	      "fail should return null values");

  // Should fail if bad args
  rc = _fc_readDataFile(NULL, &temp_numPoint, &temp_numComp, &temp_assoc, 
			&temp_mathType, &temp_name, &temp_values);
  fail_unless(rc != FC_SUCCESS, "should fail if filename is NULL");
  fail_unless(temp_numPoint == -1 && temp_numComp == -1 &&
	      temp_assoc == FC_AT_UNKNOWN && temp_mathType == FC_ET_UNKNOWN &&
	      temp_name == NULL && temp_values == NULL, 
	      "fail should return null values");
  rc = _fc_readDataFile(filenames[0], NULL, &temp_numComp, &temp_assoc, 
			&temp_mathType, &temp_name, &temp_values);
  fail_unless(rc != FC_SUCCESS, "should fail if numPoint is NULL");
  fail_unless(temp_numComp == -1 && temp_assoc == FC_AT_UNKNOWN && 
	      temp_mathType == FC_ET_UNKNOWN && temp_name == NULL && 
	      temp_values == NULL, "fail should return null values");
  rc = _fc_readDataFile(filenames[0], &temp_numPoint, NULL, &temp_assoc, 
			&temp_mathType, &temp_name, &temp_values);
  fail_unless(rc != FC_SUCCESS, "should fail if numComp is NULL");
  fail_unless(temp_numPoint == -1 && temp_assoc == FC_AT_UNKNOWN && 
	      temp_mathType == FC_ET_UNKNOWN && temp_name == NULL && 
	      temp_values == NULL, "fail should return null values");
  rc = _fc_readDataFile(filenames[0], &temp_numPoint, &temp_numComp, NULL, 
			&temp_mathType, &temp_name, &temp_values);
  fail_unless(rc != FC_SUCCESS, "should fail if assoc is NULL");
  fail_unless(temp_numPoint == -1 && temp_numComp == -1 &&
	      temp_mathType == FC_ET_UNKNOWN && temp_name == NULL && 
	      temp_values == NULL, "fail should return null values");
  rc = _fc_readDataFile(filenames[0], &temp_numPoint, &temp_numComp, 
			&temp_assoc, NULL, &temp_name, &temp_values);
  fail_unless(rc != FC_SUCCESS, "should fail if mathType is NULL");
  fail_unless(temp_numPoint == -1 && temp_numComp == -1 &&
	      temp_assoc == FC_AT_UNKNOWN && temp_name == NULL && 
	      temp_values == NULL, "fail should return null values");
  rc = _fc_readDataFile(filenames[0], &temp_numPoint, &temp_numComp, 
			&temp_assoc, &temp_mathType, NULL, &temp_values);
  fail_unless(rc != FC_SUCCESS, "should fail if name is NULL");
  fail_unless(temp_numPoint == -1 && temp_numComp == -1 &&
	      temp_assoc == FC_AT_UNKNOWN && temp_mathType == FC_ET_UNKNOWN &&
	      temp_values == NULL, "fail should return null values");
  rc = _fc_readDataFile(filenames[0], &temp_numPoint, &temp_numComp, 
			&temp_assoc, &temp_mathType, &temp_name, NULL);
  fail_unless(rc != FC_SUCCESS, "should fail if values is NULL");
  fail_unless(temp_numPoint == -1 && temp_numComp == -1 &&
	      temp_assoc == FC_AT_UNKNOWN && temp_mathType == FC_ET_UNKNOWN &&
	      temp_name == NULL, "fail should return null values");

  // Test all types of elems files
  for (i = 0; i < numFile; i++) {
    rc = _fc_readDataFile(filenames[i], &temp_numPoint, &temp_numComp, 
			  &temp_assoc, &temp_mathType, &temp_name, 
			  &temp_values);
    fail_unless(rc == FC_SUCCESS, "should not fail to read file");
    fail_unless(temp_numPoint == numPoints[i], "mismatch of numPoint");
    fail_unless(temp_numComp == numComps[i], "mismatch of numComp");
    fail_unless(temp_assoc == assocs[i], "mismatch of assoc");
    fail_unless(temp_mathType == mathTypes[i], "mismatch of mathType");
    fail_unless(!strcmp(temp_name, names[i]), "mismatch of name");
    for (j = 0; j < numComps[i]; j++) {
      fail_unless(temp_values[j] == startValues[i][j], 
		  "mismatch of 1st value");
      fail_unless(temp_values[(temp_numPoint-1)*numComps[i]+j] == endValues[i][j],
		  "mismatch of final value");
    }
    free(temp_name);
    free(temp_values);
  }
}
END_TEST

// _fc_makeMeshFromFile is a wrapper - going to assume that testing 1
// case is good enough.
START_TEST(make_mesh)
{
  FC_ReturnCode rc;
  char* vertfile = "../data/gen_point_line_tri_quad.vert";
  char* elemfile = "../data/gen_tri.elem";
  char* mismatched_elemfile = "../data/gen_hex.elem";
  FC_Dataset dataset, bad_dataset = { 100, 100 };
  FC_Mesh mesh;
  // f prefix is for "from the files"
  char* name,* fname;
  int numVert, fnumVert, numElem, fnumElem, numDim, fnumDim, numVperE;
  FC_ElementType elemType, felemType;
  double* coords, *fcoords;
  int* conns, *fconns;

  // Create a dataset to add to
  rc = fc_createDataset("happy", &dataset);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create dataset");

  // Read files & save data for later comparison
  rc = _fc_readVertexFile(vertfile, &fnumVert, &fnumDim, &fcoords);
  fail_unless(rc == FC_SUCCESS, "abort: failed to read vert file");
  rc = _fc_readElementFile(elemfile, &fnumVert, &fnumElem, &felemType,
			   &fname, &fconns);
  fail_unless(rc == FC_SUCCESS, "abort: failed to read elem file"); 

  // Test bad inputs - should fail
  rc = _fc_makeMeshFromFiles(bad_dataset, vertfile, elemfile, &mesh);
  fail_unless(rc != FC_SUCCESS, "should fail if bad dataset");
  fail_unless(FC_HANDLE_EQUIV(mesh, FC_NULL_MESH), "fail should return null");
  rc = _fc_makeMeshFromFiles(dataset, "nonsense", elemfile, &mesh);
  fail_unless(rc != FC_SUCCESS, "should fail if bad vertfile");
  fail_unless(FC_HANDLE_EQUIV(mesh, FC_NULL_MESH), "fail should return null");
  rc = _fc_makeMeshFromFiles(dataset, vertfile, "nonsense", &mesh);
  fail_unless(rc != FC_SUCCESS, "should fail if bad elem file");
  fail_unless(FC_HANDLE_EQUIV(mesh, FC_NULL_MESH), "fail should return null");
  rc = _fc_makeMeshFromFiles(dataset, vertfile, elemfile, NULL);
  fail_unless(rc != FC_SUCCESS, "should fail if no mesh to return");
  
  // Test incompatible input - vertfile & elemfile need to agree
  rc = _fc_makeMeshFromFiles(dataset, vertfile, mismatched_elemfile, &mesh);
  fail_unless(rc != FC_SUCCESS, "should fail if mismatch files");
  fail_unless(FC_HANDLE_EQUIV(mesh, FC_NULL_MESH), "fail should return null");

  // Test good input
  rc = _fc_makeMeshFromFiles(dataset, vertfile, elemfile, &mesh);
  fail_unless(rc == FC_SUCCESS, "should succeed");
  rc = fc_getMeshName(mesh, &name);
  fail_unless(rc == FC_SUCCESS, "failed to get mesh name");
  rc = fc_getMeshInfo(mesh, NULL, &numDim, &numVert, &numElem, &elemType);
  fail_unless(rc == FC_SUCCESS, "failed to get mesh info");
  fc_getMeshCoordsPtr(mesh, &coords);
  fc_getMeshElementConnsPtr(mesh, &conns);
  fail_unless(numVert == fnumVert, "mismatch of numVert");
  fail_unless(numDim == fnumDim, "mismatch of numDim");
  fail_unless(!memcmp(coords, fcoords, numVert*numDim*sizeof(double)),
	      "mismatch of coords");
  fail_unless(numElem == fnumElem, "mismatch of numElem");
  fail_unless(elemType == felemType, "mismatch of elemType");
  numVperE = fc_getElementTypeNumVertex(elemType);
  fail_unless(!strcmp(name, fname), "msimatch of name");
  fail_unless(!memcmp(conns, fconns, numElem*numVperE*sizeof(int)),
	      "mismatch of conns");
    
  // All done
  free(name);
  free(fname);
  free(fcoords);
  free(fconns);
  rc = fc_deleteDataset(dataset);
  fail_unless(rc == FC_SUCCESS, "failed to delete dataset at end of test");
}
END_TEST

// _fc_makeVarFromFile is a wrapper - going to assume that testing 1
// case is good enough.
START_TEST(make_var)
{
  FC_ReturnCode rc;
  char* vertfile = "../data/gen_point_line_tri_quad.vert";
  char* elemfile = "../data/gen_tri.elem";
  char* datafile = "../data/gen_point_line_tri_quad.data";
  char* mismatched_datafile = "../data/gen_tet_prism_hex.data";
  FC_Dataset dataset;
  FC_Mesh mesh, bad_mesh = { 100, 100 };
  FC_Variable variable;
  // f prefix is for "from the files"
  char* name, *fname;
  int numPoint, fnumPoint, numComp, fnumComp;
  FC_AssociationType assoc, fassoc;
  FC_MathType mathType, fmathType;
  FC_DataType dataType, fdataType = FC_DT_FLOAT;
  float* data, *fdata;

  // Create a dataset & mesh to add to
  rc = fc_createDataset("happy", &dataset);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create dataset");
  rc = _fc_makeMeshFromFiles(dataset, vertfile, elemfile, &mesh);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create mesh");

  // Read file & save data for later comparison
  rc = _fc_readDataFile(datafile, &fnumPoint, &fnumComp, &fassoc, &fmathType,
			&fname, &fdata);
  fail_unless(rc == FC_SUCCESS, "abort: failed to read data file");

  // Test bad inputs - should fail
  rc = _fc_makeVariableFromFile(bad_mesh, datafile, &variable);
  fail_unless(rc != FC_SUCCESS, "should fail if bad mesh");
  fail_unless(FC_HANDLE_EQUIV(variable, FC_NULL_VARIABLE), 
	      "fail should return null");
  rc = _fc_makeVariableFromFile(mesh, "nonsense", &variable);
  fail_unless(rc != FC_SUCCESS, "should fail if bad data file");
  fail_unless(FC_HANDLE_EQUIV(variable, FC_NULL_VARIABLE), 
	      "fail should return null");
  rc = _fc_makeVariableFromFile(mesh, datafile, NULL);
  fail_unless(rc != FC_SUCCESS, "should fail if bad mesh");

  // Test incompatible input - data file needs to agree with vert or elem file
  rc = _fc_makeVariableFromFile(mesh, mismatched_datafile, &variable);
  fail_unless(rc != FC_SUCCESS, "should fail if mismatched data file");
  fail_unless(FC_HANDLE_EQUIV(variable, FC_NULL_VARIABLE), 
	      "fail should return null");

  // Test good input
  rc = _fc_makeVariableFromFile(mesh, datafile, &variable);
  fail_unless(rc == FC_SUCCESS, "should suceed");
  rc = fc_getVariableName(variable, &name);
  fail_unless(rc == FC_SUCCESS, "failed to get variable name");
  rc = fc_getVariableInfo(variable, &numPoint, &numComp, &assoc, &mathType,
			  &dataType);
  fail_unless(rc == FC_SUCCESS, "failed to get variable info");
  fc_getVariableDataPtr(variable, (void**)&data);
  fail_unless(numPoint == fnumPoint, "mismatch of numPoint");
  fail_unless(numComp == fnumComp, "mismatch of numComp");
  fail_unless(assoc == fassoc, "mismatch of assoc");
  fail_unless(mathType == fmathType, "mismatch of mathType");
  fail_unless(dataType == fdataType, "mismatch of dataType");
  fail_unless(!strcmp(name, fname), "mismatch of name");
  fail_unless(!memcmp(data, fdata, numPoint*numComp*sizeof(float)), 
	      "mismatch of data");
  
  // All done
  free(name);
  free(fname);
  free(fdata);
  rc = fc_deleteDataset(dataset);
  fail_unless(rc == FC_SUCCESS, "failed to delete dataset at end of test");
 }
END_TEST

// **** helpers 	 
	  	 
START_TEST(namesfile) 	 
 { 	 
   FC_ReturnCode rc; 	 
   int i, j, ifile; 	 
   FILE *file; 	 
   int numFile = 3; 	 
   char* gooddsfilenames[3] = { "gen_namestest1.notreal", 	 
                                "gen_namestest2.notreal", 	 
                                "./gen_namestest1.notreal" }; 	 
   char* baddsfilename = "icantpossiblyexistcani.notreal"; 	 
   char* namesfilenames[2] = { "gen_namestest1.names", "gen_namestest2.names" }; 	 
   int numType = 3; 	 
   char* typeNames[3] = { "blocks", "nodesets", "sidesets" }; 	 
   int numPerType[3] = { 4, 3, 5 }; 	 
   char* names[3][5] = { { "hearts", "diamonds", "clubs", "spades" }, 	 
                         { "red with space", "yellow   ", "blue" }, 	 
                         { "one", "two", "three", "four", "five" } }; 	 
   int IDs[3][5] = { { 3, 11, 9, 7 }, 	 
                     { 100, 101, 102 }, 	 
                     { 5, 4, 3, 2, 1 } }; 	 
   // if sorted ... 	 
   int maps[3][5] = { { 0, 3, 2, 1 }, 	 
                      { 0, 1, 2 }, 	 
                      { 4, 3, 2, 1, 0 } }; 	 
   FC_SortedBlobArray temp_names[3]; 	 
   _FC_IntStringPair *blob; 	 
  	 
   // generate test data 	 
   file = fopen(namesfilenames[0], "w"); 	 
   fail_unless(file != NULL, "abort: failed to open file for write"); 	 
   for (i = 0; i < numType; i++) { 	 
     fprintf(file, "%s\n", typeNames[i]); 	 
     for (j = 0; j < numPerType[i]; j++) 	 
       fprintf(file, "%-6d %s\n", IDs[i][j], names[i][j]); 	 
   } 	 
   fclose(file); 	 
   // write in reverse order - will be interprested the same! 	 
   file = fopen(namesfilenames[1], "w"); 	 
   fail_unless(file != NULL, "abort: failed to open file for write"); 	 
   for (i = numType-1; i >= 0; i--) { 	 
     fprintf(file, "%s\n", typeNames[i]); 	 
     for (j = numPerType[i]-1; j >= 0; j--) 	 
       fprintf(file, "%-6d %s\n", IDs[i][j], names[i][j]); 	 
   } 	 
   fclose(file); 	 
  	 
   // it's o.k. to call on nonexistent file 	 
   rc = _fc_readNamesFile(baddsfilename, &temp_names[0], &temp_names[1], 	 
                          &temp_names[2]); 	 
   fail_unless(rc == FC_SUCCESS, "failed to read names file"); 	 
   for (i = 0; i < numType; i++) 	 
     fail_unless(temp_names[i].numBlob == 0, "if no file, should return empties"); 	 
  	 
   // test files 	 
   for (ifile = 0; ifile < numFile; ifile++) { 	 
     rc = _fc_readNamesFile(gooddsfilenames[ifile], &temp_names[0], 	 
                            &temp_names[1], &temp_names[2]); 	 
     fail_unless(rc == FC_SUCCESS, "failed to read names file"); 	 
     for (i = 0; i < numType; i++) { 	 
       fail_unless(temp_names[i].numBlob == numPerType[i], 	 
                   "mismatch of num of that type"); 	 
       for (j = 0; j < numPerType[i]; j++) { 	 
         blob = temp_names[i].blobs[j]; 	 
         fail_unless(blob->id == IDs[i][maps[i][j]], 	 
                     "mismatch of block ID"); 	 
         fail_unless(!strcmp(blob->string, names[i][maps[i][j]]), 	 
                     "mismatch of name"); 	 
       } 	 
     } 	 
  	 
     // test free 	 
     _fc_freeNamesFileNames(&temp_names[0], &temp_names[1], &temp_names[2]); 	 
   } 	 
 } 	 
END_TEST
// **** test init/final for io

START_TEST(fileiotype)
{
  int i;
  int numInput = 5;
  int inputs[5] = { FC_FT_NONE, FC_FT_EXODUS, FC_FT_LSDYNA, 
		    FC_FT_NONE - 1, FC_FT_LSDYNA + 1 };
  int validOuts[5] = { 1, 1, 1, 0, 0 };
  char* textOuts[5] = { "FC_FT_NONE", "FC_FT_EXODUS", "FC_FT_LSDYNA",
			"invalid value for enum FC_FileIOType",
			"invalid value for enum FC_FileIOType" };
  int isReadable[5] = { 1, 1, 1, 0, 0 };
  int isWritable[5] = { 0, 1, 0, 0, 0 };

#ifndef HAVE_EXODUS
  isReadable[1] = 0;
  isWritable[1] = 0;
#endif
  
  for (i = 0; i < numInput; i++) {
    fail_unless(fc_isFileIOTypeValid(inputs[i]) == validOuts[i],
		"isValid mismatch");
    fail_unless(!strcmp(fc_getFileIOTypeText(inputs[i]), textOuts[i]),
		"getText mismatch");
    fail_unless(fc_isFileIOTypeReadSupported(inputs[i]) == isReadable[i],
		"isFileIOTypeReadSupported mismatch");
    fail_unless(fc_isFileIOTypeWriteSupported(inputs[i]) == isWritable[i],
		"isFileIOTypeWriteSupported mismatch");
  }
}
END_TEST

START_TEST(init_final)
{
  if (isForking) {
    FC_ReturnCode rc;
    
    // Init
    fc_setLibraryVerbosity(fc_messages);
    rc = fc_initLibrary();
    fail_unless(rc == FC_SUCCESS, "failed to init library");
    
    // FIX? - check io initialization done properly
    
    // Final
    rc = fc_finalLibrary();
    fail_unless(rc == FC_SUCCESS, "failed to final library");
  }
  else {
    printf("**CK_NOFORK: Can't test init/final for IO\n");
  }
}
END_TEST

// NOTE: We have to depend on a specific file IO for this test, and
// Exodus is the least worst choice.
// Should really have after testing Exodus IO, but conceptually makes
// sense to have in this test case.
#ifdef HAVE_EXODUS
START_TEST(rewrite_test)
{
  FC_ReturnCode rc;
  int i, j, k, m;
  FC_Dataset origDataset, newDataset;
  FC_Sequence sequence, *temp_sequences;
  FC_Mesh *temp_meshes;
  FC_Subset *temp_subsets;
  FC_Variable *temp_glbVars, **temp_glbSeqVars, *temp_vars, **temp_seqVars;
  char origfilename[20] = "temp.ex2", newfilename[20] = "temp2.ex2";
  char rmcommand[20] = "rm -f temp2.ex2"; // make sure this is same as newfilename
  char* temp_name;
  int numSeq = 1, numGlbVar = 2, numGlbSeqVar = 2, numMesh = 2;
  int numSubset = 2, numVar = 2, numSeqVar = 2;
  int numStep = 3;
  // 1 hex, 8 hexess
  int numDim = 3, dims[2][3] = { { 1, 1, 1 }, { 2, 2, 2 } };
  int numVertPerElem = 8, numVerts[2] = { 8, 27 }, numElems[2] = { 1, 8 };
  FC_Coords lowers = { 0, 0, 0 }, uppers = { 1, 1, 1 };
  int temp_numSeq, temp_numGlbVar, temp_numGlbSeqVar, temp_numMesh;
  int temp_numSubset, temp_numVar, temp_numSeqVar;
  int temp_numStep, *temp_numStepPerSeqVar, temp_numDataPoint, temp_numComp;
  int temp_topodim, temp_numDim, temp_numVert, temp_numElem;
  int temp_numMember, temp_maxNumMember, *temp_memberIDs;
  FC_DataType dataType = FC_DT_DOUBLE, temp_dataType;
  FC_MathType mathType = FC_MT_SCALAR, temp_mathType;
  FC_AssociationType assocs[2] = { FC_AT_VERTEX, FC_AT_ELEMENT };
  FC_AssociationType temp_assoc;
  FC_ElementType elemType = FC_ET_HEX, temp_elemType;
  int numMembers[2] = { 2, 3 }, memberIDs[2][3] = { { 0, 2 }, { 0, 2, 4 } };
  char ds_name[20] = "blue berry";
  char seq_name[20] = "time";
  char glbVar_names[2][20] = { "glbVar1", "glbVar2" };
  char glbSeqVar_names[2][20] = { "glbSeqVar1", "glbSeqVar2" };
  char mesh_names[2][20] = { "mesh1", "mesh2" };
  char subset_names[2][20] = { "subset1", "subset2" };
  char var_names[2][20] = { "var1", "var2" };
  char seqVar_names[2][20] = { "seqVar1", "seqVar2" };
  double data[100]; // 100 is arbitrary, should be way big enough
  int conns[2][8*8], *temp_conns;
  double coords[2][3*27], *temp_coords;
  void* temp_data;

  // call setup fixture manually (didn't want for all tests in this tcase)
  fileio_setup();

  // setup
  for (i = 0; i < 100; i ++) 
    data[i] = i/10.;

  // Make simple dataset with 2 of everything on it (except 1 sequence)
  // and only vertex not elem subset for mesh 0
  rc = fc_createDataset(ds_name, &origDataset);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create dataset");
  // sequence
  rc = fc_createSequence(origDataset, seq_name, &sequence);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create sequence");
  rc = fc_setSequenceCoords(sequence, numStep, dataType, data);
  fail_unless(rc == FC_SUCCESS, "abort: failed to set sequence coords");
  // global vars
  for (i = 0; i < numGlbVar; i++) {
    FC_Variable glbVar;
    rc = fc_createGlobalVariable(origDataset, glbVar_names[i], &glbVar);
    fail_unless(rc == FC_SUCCESS, "abort: failed to create global variable");
    rc = fc_setVariableData(glbVar, 1, 1, FC_AT_WHOLE_DATASET,
			    mathType, dataType, data);
    fail_unless(rc == FC_SUCCESS, "abort: failed to set global var data");
  }
  // global seq vars
  for (i = 0; i < numGlbSeqVar; i++) {
    FC_Variable *glbSeqVar;
    rc = fc_createGlobalSeqVariable(origDataset, sequence, glbSeqVar_names[i], 
                                    &temp_numStep, &glbSeqVar);
    fail_unless(rc == FC_SUCCESS, "abort: failed to create global seqVar");
    for (j = 0; j < temp_numStep; j++) {
      rc = fc_setVariableData(glbSeqVar[j], 1, 1, FC_AT_WHOLE_DATASET,
			      mathType, dataType, data);
      fail_unless(rc == FC_SUCCESS, "abort: failed to set global seqVar data");
    }
    free(glbSeqVar);
  }
  // meshes
  for (i = 0; i < numMesh; i++) {
    FC_Mesh mesh;
    rc = fc_createSimpleHexMesh(origDataset, mesh_names[i], dims[i][0], 
                                dims[i][1], dims[i][2], lowers, uppers, &mesh);
    fail_unless(rc == FC_SUCCESS, "abort: failed to create simpl hex mesh");
    rc = fc_getMeshElementConnsPtr(mesh, &temp_conns);
    fail_unless(rc == FC_SUCCESS, "abort: failed to get mesh conns");
    memcpy(conns[i], temp_conns, numElems[i]*numVertPerElem*sizeof(int));
    rc = fc_getMeshCoordsPtr(mesh, &temp_coords);
    fail_unless(rc == FC_SUCCESS, "abort: failed to get mesh coords");
    memcpy(coords[i], temp_coords, numVerts[i]*numDim*sizeof(double));
    // subsets
    for (j = 0; j < numSubset; j++) {
      FC_Subset subset;
      rc = fc_createSubset(mesh, subset_names[j], assocs[j], &subset);
      fail_unless(rc == FC_SUCCESS, "abort: failed to create subset");
      if (assocs[j] == FC_AT_ELEMENT && i == 0){//just one membe
	rc = fc_addMemberToSubset(subset, 0);
      } else {
	for (k = 0; k < numMembers[j]; k++) {
	  rc = fc_addMemberToSubset(subset, memberIDs[j][k]);
	  fail_unless(rc == FC_SUCCESS, "abort: failed to add members");
	}
	//otherwise dont add the subset
      }
    }
    // vars
    for (j = 0; j < numVar; j++) {
      int numDataPoint;
      FC_Variable var;
      rc = fc_createVariable(mesh, var_names[j], &var);
      fail_unless(rc == FC_SUCCESS, "abort: failed to create variable");
      //its either a nodal or an elem var
      if (assocs[j] == FC_AT_VERTEX)
        numDataPoint = numVerts[i];
      else
        numDataPoint = numElems[i];
      rc = fc_setVariableData(var, numDataPoint, 1, assocs[j],
                              mathType, dataType, data);
      fail_unless(rc == FC_SUCCESS, "abort: failed to set var data");
    }
    // seqVars
    for (j = 0; j < numSeqVar; j++) {
      int numDataPoint;
      FC_Variable *seqVar;
      rc = fc_createSeqVariable(mesh, sequence, seqVar_names[j], 
                                &temp_numStep, &seqVar);
      fail_unless(rc == FC_SUCCESS, "abort: failed to create seqVariable");
      if (assocs[j] == FC_AT_VERTEX)
        numDataPoint = numVerts[i];
      else
        numDataPoint = numElems[i];
      for (k = 0; k < temp_numStep; k++) {
        rc = fc_setVariableData(seqVar[k], numDataPoint, 1, assocs[j], 
                                mathType, dataType, data);
        fail_unless(rc == FC_SUCCESS, "abort: failed to set seqVar data");
      }
      free(seqVar);
    }
  }

  
  for (m = 0; m < 2; m++) { // two different cases

    //** Case 1 - rewrite in core dataset
    if (m == 0) {
      rc = fc_writeDataset(origDataset, origfilename, FC_FT_EXODUS);
      fail_unless(rc == FC_SUCCESS, "abort: failed to write dataset");
    }
    //** Case 2 - rewrite loaded dataset
    else {
      rc = fc_loadDataset(origfilename, &origDataset);
      fail_unless(rc == FC_SUCCESS, "abort: failed to load dataset");
    }
      
    // rewrite - (first make sure that a file doesn't already exist)
    system(rmcommand);
    rc = fc_loadDataset(newfilename, &newDataset);
    fail_unless(rc != FC_SUCCESS, "should fail, we just deleted it");
    rc = fc_rewriteDataset(origDataset, newfilename, FC_FT_EXODUS);
    fail_unless(rc == FC_SUCCESS, "failed to rewrite dataset");
    // delete origDataset
    fc_deleteDataset(origDataset);

    // reload
    rc = fc_loadDataset(newfilename, &newDataset);
    fail_unless(rc == FC_SUCCESS, "failed to load rewritten dataset");
    
    // test contents
    // test name
    rc = fc_getDatasetName(newDataset, &temp_name);
    fail_unless(rc == FC_SUCCESS, "failed to get name of copy");
    fail_unless(!strcmp(temp_name, ds_name), "mismatch of name");
    free(temp_name);
    // test sequences
    rc = fc_getSequences(newDataset, &temp_numSeq, &temp_sequences);
    fail_unless(rc == FC_SUCCESS, "should not fail to get sequences");
    fail_unless(temp_numSeq == numSeq, "mismatch of numSequence");
    for (i = 0; i < numSeq; i++) {
      rc = fc_getSequenceName(temp_sequences[i], &temp_name);
      fail_unless(rc == FC_SUCCESS, "failed to get sequence name");
      fail_unless(!strcmp(temp_name, seq_name), "seq name mismatch");
      free(temp_name);
      rc = fc_getSequenceInfo(temp_sequences[i], &temp_numStep, &temp_dataType);
      fail_unless(rc == FC_SUCCESS, "failed to get sequence info");
      fail_unless(temp_numStep == numStep, "mismatch of numStep");
      fail_unless(temp_dataType == dataType, "mismatch of dataType");
      rc = fc_getSequenceCoordsPtr(temp_sequences[i], &temp_data);
      fail_unless(rc == FC_SUCCESS, "failed to get sequence coords");
      fail_unless(!memcmp(temp_data, data, numStep*sizeof(double)),
		  "mismatch of data");
    }
    free(temp_sequences);
    // test global vars - FIX some day - aren't written to Exodus
    rc = fc_getGlobalVariables(newDataset, &temp_numGlbVar, &temp_glbVars);
    fail_unless(rc == FC_SUCCESS, "should not fail to get global vars");
    fail_unless(temp_numGlbVar == 0, "mismatch of numGlbVar");
    // test global seq vars
    rc = fc_getGlobalSeqVariables(newDataset, &temp_numGlbSeqVar,
				  &temp_numStepPerSeqVar, &temp_glbSeqVars);
    fail_unless(rc == FC_SUCCESS, "should not fail to get glbSeqVars");
    fail_unless(temp_numGlbSeqVar == numGlbSeqVar, "mismatch of numGlbSeqVar");
    for (i = 0; i < numGlbSeqVar; i++) {
      fail_unless(temp_numStepPerSeqVar[i] == numStep, "mismatch of numStep");
      rc = fc_getVariableName(temp_glbSeqVars[i][0], &temp_name);
      fail_unless(rc == FC_SUCCESS, "failed to get glbSeqVar name");
      fail_unless(!strcmp(temp_name, glbSeqVar_names[i]), 
		  "glbSeqVar name mismatch");
      free(temp_name);
      rc = fc_getVariableInfo(temp_glbSeqVars[i][0], &temp_numDataPoint,
			      &temp_numComp, &temp_assoc, &temp_mathType, 
			      &temp_dataType);
      fail_unless(rc == FC_SUCCESS, "failed to get glbSeqVar info");
      fail_unless(temp_numDataPoint == 1, "mismatch of numDataPoints");
      fail_unless(temp_numComp == 1, "mismatch of numComp");
      fail_unless(temp_assoc == FC_AT_WHOLE_DATASET, "mismatch of assoc");
      fail_unless(temp_mathType == mathType, "mismatch of mathType");
      fail_unless(temp_dataType == dataType, "mismatch of dataType");
      for (j = 0; j < temp_numStepPerSeqVar[i]; j++) {
	rc = fc_getVariableDataPtr(temp_glbSeqVars[i][j], &temp_data);
	fail_unless(rc == FC_SUCCESS, "failed to get glbSeqVar data");
	fail_unless(!memcmp(temp_data, data, sizeof(double)),
		    "mismatch of glbSeqVar data");
      }
      free(temp_glbSeqVars[i]);
    }
    free(temp_numStepPerSeqVar);
    free(temp_glbSeqVars);
    // test meshes
    rc = fc_getMeshes(newDataset, &temp_numMesh, &temp_meshes);
    fail_unless(rc == FC_SUCCESS, "should not fail to get meshes");
    fail_unless(temp_numMesh == numMesh, "mismatch of numMesh");
    for (i = 0; i < numMesh; i++) {
      rc = fc_getMeshName(temp_meshes[i], &temp_name);
      fail_unless(rc == FC_SUCCESS, "failed to get mesh name");
      fail_unless(!strcmp(temp_name, mesh_names[i]), "mesh name mismatch");
      free(temp_name);
      rc = fc_getMeshInfo(temp_meshes[i], &temp_topodim, &temp_numDim,
			  &temp_numVert, &temp_numElem, &temp_elemType);
      fail_unless(rc == FC_SUCCESS, "failed to get mesh info");
      fail_unless(temp_topodim == numDim, "mismatch of topodim");
      fail_unless(temp_numDim == numDim, "mismatch of numdim");
      fail_unless(temp_numVert == numVerts[i], "mismatch of numVert");
      fail_unless(temp_numElem == numElems[i], "mismatch of numElem");
      fail_unless(temp_elemType == elemType, "mismatch of elemType");
      rc = fc_getMeshCoordsPtr(temp_meshes[i], &temp_coords);
      fail_unless(rc == FC_SUCCESS, "failed to get mesh coords");
      fail_unless(!memcmp(temp_coords, coords[i], numDim*numVerts[i]*sizeof(double)),
		  "mismatch of coords");
      rc = fc_getMeshElementConnsPtr(temp_meshes[i], &temp_conns);
      fail_unless(rc == FC_SUCCESS, "failed to get mesh coords");
      fail_unless(!memcmp(temp_conns, conns[i], 
			  numVertPerElem*numElems[i]*sizeof(int)),
		  "mismatch of conns");

      // test subsets
      rc = fc_getSubsets(temp_meshes[i], &temp_numSubset, &temp_subsets);
      fail_unless(rc == FC_SUCCESS, "should not fail to get subsets");
      fail_unless(temp_numSubset == numSubset, "mismatch of numSubset");
      //one is nodal the other is elem
      {
	int foundNodal = 0;
	int foundElem = 0;
	int compid = 0;

	for (j = 0; j < numSubset; j++) {
	  rc = fc_getSubsetName(temp_subsets[j], &temp_name);
	  fail_unless(rc == FC_SUCCESS, "failed to get subset name");
	  rc = fc_getSubsetInfo(temp_subsets[j], &temp_numMember,
				&temp_maxNumMember, &temp_assoc);
	  fail_unless(rc == FC_SUCCESS, "failed to get subset info");
	  switch (temp_assoc){
	  case FC_AT_VERTEX:
	    foundNodal = 1;
	    fail_unless(temp_maxNumMember == numVerts[i], 
			"mismatch of maxNumMember");
	    compid = (assocs[j] == FC_AT_VERTEX ? j: !j);
	    break;
	  case FC_AT_ELEMENT:
	    foundElem = 1;
	    fail_unless(temp_maxNumMember == numElems[i], 
			"mismatch of maxNumMember");
	    compid = (assocs[j] == FC_AT_ELEMENT ? j: !j);
	    break;
	  default:
	    fail_unless( 0 == 1, "wrong return assoc for subset");
	    break;
	  }
	  fail_unless(!strcmp(temp_name, subset_names[compid]), "subset name mismatch");
	  free(temp_name);
	  if (i == 0 && temp_assoc == FC_AT_ELEMENT){
	    fail_unless( temp_numMember == 1, "wrong num members for subset");
	  } else {
	    fail_unless(temp_numMember == numMembers[compid], "mismatch of numMember");
	    rc = fc_getSubsetMembersAsArray(temp_subsets[j], &temp_numMember, 
					    &temp_memberIDs);
	    fail_unless(rc == FC_SUCCESS, "failed to get member ids");
	    
	    fail_unless(!memcmp(temp_memberIDs, memberIDs[compid], numMembers[compid]*sizeof(int)),
			"mismatch of memberIDs");
	    free(temp_memberIDs);
	  }
	}
	free(temp_subsets);
	if (i == 0 ){
	  fail_unless(foundNodal, "missing return assoc for subsets");
	} else {
	  fail_unless(foundNodal && foundElem, "missing return assoc for subsets");
	}
      }


      // test vars -  - will see the "ID" var that exodus creates,
      // any elem attr and node attr as vars so add one of each.
      rc = fc_getVariables(temp_meshes[i], &temp_numVar, &temp_vars);
      fail_unless(rc == FC_SUCCESS, "should not fail to get vars");
      fail_unless(temp_numVar == 3, "mismatch of numvar");
      //one must be ID, others are elem or nodal var
      {
	int foundID = 0;
	int foundElem = 0;
	int foundNodal = 0;
	int elemVarIndex = (assocs[0] == FC_AT_ELEMENT ? 0:1);
	int nodalVarIndex = !elemVarIndex;
	for (j = 0; j < temp_numVar; j++) {
	  rc = fc_getVariableName(temp_vars[j], &temp_name);
	  if (!strcmp(temp_name, "ID")){
	    free(temp_name);
	    rc = fc_getVariableInfo(temp_vars[j], &temp_numDataPoint,
				    &temp_numComp, &temp_assoc, &temp_mathType, 
				    &temp_dataType);
	    fail_unless(rc == FC_SUCCESS, "failed to get var info");
	    fail_unless(temp_numDataPoint == 1, "mismatch of numDataPoint");
	    fail_unless(temp_numComp == 1, "mismatch of numComp");
	    fail_unless(temp_assoc == FC_AT_WHOLE_MESH, "mismatch of assoc");
	    fail_unless(temp_mathType == FC_MT_SCALAR, "mismatch of mathtype");
	    fail_unless(temp_dataType == FC_DT_INT, "mismatch of datatype");
	    rc = fc_getVariableDataPtr(temp_vars[j], &temp_data);
	    fail_unless(rc == FC_SUCCESS, "failed to get var data");
	    fail_unless(((int*)temp_data)[0] == i+1, "mismatch of whole var data");
	    foundID = 1;
	  } else {
	    rc = fc_getVariableInfo(temp_vars[j], &temp_numDataPoint,
				    &temp_numComp, &temp_assoc, &temp_mathType, 
				    &temp_dataType);
	    fail_unless(rc == FC_SUCCESS, "failed to get var info");
	    fail_unless(temp_numComp == 1, "mismatch of numComp");
	    fail_unless(temp_mathType == FC_MT_SCALAR, "mismatch of mathtype");
	    fail_unless(temp_dataType == FC_DT_DOUBLE, "mismatch of datatype");
	    rc = fc_getVariableDataPtr(temp_vars[j], &temp_data);
	    fail_unless(rc == FC_SUCCESS, "failed to get var data");
	    for (k = 0; k < temp_numDataPoint; k++){
	      fail_unless(((double*)temp_data)[i] == data[i], "mismatch of elem var data");
	    }
	    if (!strcmp(temp_name, var_names[elemVarIndex])){
	      fail_unless(temp_numDataPoint == numElems[i], "mismatch of numDataPoint");
	      fail_unless(temp_assoc == FC_AT_ELEMENT, "mismatch of assoc");
	      foundElem = 1;
	    } else if (!strcmp(temp_name, var_names[nodalVarIndex])){
	      fail_unless(temp_numDataPoint == numVerts[i], "mismatch of numDataPoint");
	      fail_unless(temp_assoc == FC_AT_VERTEX, "mismatch of assoc");
	      foundNodal = 1;
	    } else {
	      fail_unless( 1 , "invalid var");
	    }
	    free(temp_name);
	  }
	}
	fail_unless(foundID && foundElem, "invalid vars");
	fail_unless(foundNodal, "invalid vars");
      }
      free(temp_vars);

      // test seq vars - are only on last mesh & last sequence
      rc = fc_getSeqVariables(temp_meshes[i], &temp_numSeqVar,
			      &temp_numStepPerSeqVar, &temp_seqVars);
      fail_unless(rc == FC_SUCCESS, "should not fail to get seqVars");
      fail_unless(temp_numSeqVar == numSeqVar, "mismatch of numSeqVar");
      for (j = 0; j < numSeqVar; j++) {
	int numDataPoint;
	fail_unless(temp_numStepPerSeqVar[j] == numStep, "mismatch of numStep");
	rc = fc_getVariableName(temp_seqVars[j][0], &temp_name);
	fail_unless(rc == FC_SUCCESS, "failed to get seqVar name");
	fail_unless(!strcmp(temp_name, seqVar_names[j]), "seqVar name mismatch");
	free(temp_name);
	rc = fc_getVariableInfo(temp_seqVars[j][0], &temp_numDataPoint,
				&temp_numComp, &temp_assoc, &temp_mathType, 
				&temp_dataType);
	fail_unless(rc == FC_SUCCESS, "failed to get seqVar info");
	if (assocs[j] == FC_AT_VERTEX)
	  numDataPoint = numVerts[i];
	else
	  numDataPoint = numElems[i];
	fail_unless(temp_numDataPoint == numDataPoint, 
		    "mismatch of numDataPoint");
	fail_unless(temp_numComp == 1, "mismatch of numComp");
	fail_unless(temp_assoc == assocs[j], "mismatch of assoc");
	fail_unless(temp_mathType == mathType, "mismatch of mathType");
	fail_unless(temp_dataType == dataType, "mismatch of dataType");
	for (k = 0; k < temp_numStepPerSeqVar[j]; k++) {
	  rc = fc_getVariableDataPtr(temp_seqVars[j][k], &temp_data);
	  fail_unless(rc == FC_SUCCESS, "failed to get seqVar data");
	  fail_unless(!memcmp(temp_data, data, numDataPoint*sizeof(double)),
		      "mismatch of seqVar data");
	}
	free(temp_seqVars[j]);
      }
      free(temp_numStepPerSeqVar);
      free(temp_seqVars);
    }
    free(temp_meshes);
  
    // all done 
    rc = fc_deleteDataset(newDataset);
    fail_unless(rc == FC_SUCCESS, "failed to delete rewritten dataset");
  } // end of loop over cases

  // call teardown fixture manually (didn't want tof all tests in this tcase)
  fileio_teardown();
}
END_TEST
#endif // HAVE_EXODUS

// ***** Exodus

#ifdef HAVE_EXODUS

// Test passing bad args to private functions (because these usually
// called via a wrapper so no other way to test)
START_TEST(exo_private_bad_args)
{
  FC_ReturnCode rc;
  FC_Dataset dataset, badDataset = { 100, 100 };
  FC_Sequence badSequence = { 100, 100 };
  FC_Mesh badMesh = { 100, 100 };
  FC_Subset badSubset = { 100, 100 };
  FC_Variable badVariable = { 100, 100 };


  // _fc_loadExodusDataset()
  rc = _fc_loadExodusDataset(NULL, &dataset);
  fail_unless(rc != FC_SUCCESS, "null should fail");
  fail_unless(FC_HANDLE_EQUIV(dataset, FC_NULL_DATASET),
	      "fail should return null");
  rc = _fc_loadExodusDataset(EXODUS_FILE_NAME, NULL);
  fail_unless(rc != FC_SUCCESS, "null should fail");

  // _fc_readExodusXBigData()
  rc = _fc_readExodusSequenceCoords(badSequence);
  fail_unless(rc != FC_SUCCESS, "bad seq should fail");
  rc = _fc_readExodusMeshCoords(badMesh);
  fail_unless(rc != FC_SUCCESS, "bad mesh should fail");
  rc = _fc_readExodusMeshElemConns(badMesh);
  fail_unless(rc != FC_SUCCESS, "bad mesh should fail");
  rc = _fc_readExodusSubsetMembers(badSubset);
  fail_unless(rc != FC_SUCCESS, "bad subset should fail");
  rc = _fc_readExodusVariableData(badVariable);
  fail_unless(rc != FC_SUCCESS, "bad var should fail");
  rc = _fc_readExodusAttributeData(badVariable);
  fail_unless(rc != FC_SUCCESS, "bad var should fail");
  rc = _fc_readExodusOneStepVariableData(badVariable);
  fail_unless(rc != FC_SUCCESS, "bad var should fail");

  // _fc_writeExodusDataset()
  // create a temp dataset for testing
  rc = fc_createDataset("temp", &dataset);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create dataset");
  rc = _fc_writeExodusDataset(badDataset, EXODUS_FILE_NAME);
  fail_unless(rc != FC_SUCCESS, "bad dataset should fail");
  rc = _fc_writeExodusDataset(dataset, NULL);
  fail_unless(rc != FC_SUCCESS, "null name should fail");
  fc_deleteDataset(dataset);

  // _fc_closeExodusDataset()
  rc = _fc_closeExodusDataset(badDataset);
  fail_unless(rc != FC_SUCCESS, "bad datset should fail");
}
END_TEST

START_TEST(exo_helpers)
{
  FC_ReturnCode rc;
  int i, j;
  int numCase = 4;
  int numComps[4] = { 1, 3, 2, 4 };
  FC_MathType mathTypes[4] = { FC_MT_SCALAR, FC_MT_VECTOR, FC_MT_SYMTENSOR,
                               FC_MT_TENSOR };
  char name[10] = "fruity", *compName, tempName[1024];
  char vector_endings[3][2] = { "x", "y", "z" };
  char longname[MAX_STR_LENGTH+5];

  // test _fc_createExoVarCompName()
  for (i = 0; i < numCase; i++) {
    for (j = 0; j < numComps[i]; j++) {
      if (mathTypes[i] == FC_MT_SCALAR)
        strcpy(tempName, name);
      else if (mathTypes[i] == FC_MT_VECTOR)
        sprintf(tempName, "%s_%s", name, vector_endings[j]);
      else
        sprintf(tempName, "%s_c%d", name, j);
      rc = _fc_createExoVarCompName(name, numComps[i], mathTypes[i], j,
                                    &compName);
      fail_unless(rc == FC_SUCCESS, "failed to create comp name");
      fail_unless(!strcmp(compName, tempName), "mismatch of comp name");
      free(compName);
    }
  }

  compName = NULL;
  for (i = 0; i < MAX_STR_LENGTH+5; i++){
    longname[i] = 'a';
  }
  longname[MAX_STR_LENGTH+4] = '\0';
  rc = _fc_createExoVarCompName(longname, numComps[0], mathTypes[0],
				0, &compName);
  fail_unless(rc == FC_SUCCESS, "should work for long name");
  fail_unless(strlen(compName) <= MAX_STR_LENGTH, "should truncate length");
  free(compName);

  // test bad args to _fc_createExoVarCompName()
  rc = _fc_createExoVarCompName(NULL, 1, FC_MT_SCALAR, 0, &compName);
  fail_unless(rc != FC_SUCCESS, "should fail if null input"); 
  fail_unless(compName == NULL, "fail should return NULL");
  rc = _fc_createExoVarCompName(name, 0, FC_MT_SCALAR, 0, &compName);
  fail_unless(rc != FC_SUCCESS, "should fail if numComp < 1"); 
  fail_unless(compName == NULL, "fail should return NULL");
  rc = _fc_createExoVarCompName(name, 1, FC_MT_UNKNOWN, 0, &compName);
  fail_unless(rc != FC_SUCCESS, "should fail if unknown mathType"); 
  fail_unless(compName == NULL, "fail should return NULL");
  rc = _fc_createExoVarCompName(name, 1, -99, 0, &compName);
  fail_unless(rc != FC_SUCCESS, "should fail not possible mathType"); 
  fail_unless(compName == NULL, "fail should return NULL");
  rc = _fc_createExoVarCompName(name, 1, FC_MT_SCALAR, -1, &compName);
  fail_unless(rc != FC_SUCCESS, "should fail if comp id is out of range"); 
  fail_unless(compName == NULL, "fail should return NULL");
  rc = _fc_createExoVarCompName(name, 1, FC_MT_SCALAR, 1, &compName);
  fail_unless(rc != FC_SUCCESS, "should fail if comp id is out of range"); 
  fail_unless(compName == NULL, "fail should return NULL");
  rc = _fc_createExoVarCompName(name, 1, FC_MT_SCALAR, 1, NULL);
  fail_unless(rc != FC_SUCCESS, "should fail if Null output"); 
}
END_TEST

START_TEST(exo_dataset)
{
  FC_ReturnCode rc;
  FC_Dataset dataset = { 999, 999 };
  char* ds_name = "dataset";

  // Should fail to load non existent file or file of wrong type
  rc = fc_loadDatasetWithFormat(BAD_FILE_NAME, FC_FT_EXODUS, &dataset);
  fail_unless(rc != FC_SUCCESS, "should fail to load nonexistant file");
  fail_unless(FC_HANDLE_EQUIV(dataset, FC_NULL_DATASET),
	      "fail should return null handle");
  rc = fc_loadDatasetWithFormat(WRONG_FILE_NAME, FC_FT_EXODUS, &dataset);
  fail_unless(rc != FC_SUCCESS, "should fail to load wrong type of file");
  fail_unless(FC_HANDLE_EQUIV(dataset, FC_NULL_DATASET),
	      "fail should return null handle");
   
  // Should fail if write called with bad dataset
  rc = fc_writeDataset(dataset, EXODUS_FILE_NAME, FC_FT_EXODUS);
  fail_unless(rc != FC_SUCCESS, "should fail to write bad dataset");

  // create a dataset
  rc = fc_createDataset(ds_name, &dataset);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create dataset");
  
  // should fail to write because it is empty
  rc = fc_writeDataset(dataset, EXODUS_FILE_NAME, FC_FT_EXODUS);
  fail_unless(rc != FC_SUCCESS, "should fail to write empty dataset");

  // all done
  rc = fc_deleteDataset(dataset);
  fail_unless(rc == FC_SUCCESS, "failed to delete dataset at end");
}
END_TEST

START_TEST(exo_mesh)
{
  FC_ReturnCode rc;
  int i, j, k;
  FC_Dataset dataset;
  int temp_numMesh;
  FC_Mesh mesh, *meshes;
  char* ds_name = "dataset", *mesh_name;
  int numVert, numElem, numDim, topoDim, numVperE;
  int* conns, *temp_conns;
  double* coords, *temp_coords;
  char* temp_name;
  int temp_topodim, temp_numDim, temp_numVert, temp_numElem;
  FC_ElementType elemType, temp_elemType;

  // Test writing & reading a single mesh data w/variety of dims
  for (i = 0; i < numMeshType; i++) {
    for (j = 1; j <= 3; j++) { // j => numDim

      // skip if dim too small for this type
      if (topoDims[i] > j)
	continue;

      // Setup      
      rc = _fc_readVertexFile(vertFiles[i], &numVert, &numDim, &coords);
      fail_unless(rc == FC_SUCCESS, "abort: failed to read vert file");
      rc = _fc_readElementFile(elemFiles[i], &temp_numVert, &numElem, &elemType,
			       &mesh_name, &conns);
      fail_unless(rc == FC_SUCCESS, "abort: failed to read elem file");
      fail_unless(temp_numVert == numVert, "abort: numVerts don't agree");
      topoDim = fc_getElementTypeTopoDim(elemType);
      numVperE = fc_getElementTypeNumVertex(elemType);
      if (j != numDim) {
	temp_coords = malloc(numVert*j*sizeof(double));
	for (k = 0; k < numVert; k++)
	  memcpy(&temp_coords[k*j], &coords[k*numDim], j*sizeof(double));
	free(coords);
	coords = temp_coords;
      }
      numDim = j;

      // Create the dataset & mesh
      rc = fc_createDataset(ds_name, &dataset);
      fail_unless(rc == FC_SUCCESS, "abort: failed to create dataset");
      rc = fc_createMesh(dataset, mesh_name, &mesh);
      fail_unless(rc == FC_SUCCESS, "abort: failed to create mesh");
      rc = fc_setMeshCoords(mesh, numDim, numVert, coords);
      fail_unless(rc == FC_SUCCESS, "abort: failed to set coords");
      rc = fc_setMeshElementConns(mesh, elemType, numElem, conns);
      fail_unless(rc == FC_SUCCESS, "abort: failed to set conns");

      // test writing
      rc = fc_writeDataset(dataset, EXODUS_FILE_NAME, FC_FT_EXODUS);
      fail_unless(rc == FC_SUCCESS, "should be able to write dataset");
      
      // test close & reopen
      rc = fc_deleteDataset(dataset);
      fail_unless(rc == FC_SUCCESS, "should close dataset w/o error");
      rc = fc_loadDatasetWithFormat(EXODUS_FILE_NAME, FC_FT_EXODUS, &dataset);
      fail_unless(rc == FC_SUCCESS, "should load dataset");
      
      // check dataset
      rc = fc_getDatasetName(dataset, &temp_name);
      fail_unless(rc == FC_SUCCESS, "should get dataset name");
      fail_unless(!strcmp(temp_name, ds_name), "ds name mismatch");
      free(temp_name);
      rc = fc_getMeshes(dataset, &temp_numMesh, &meshes);
      fail_unless(rc == FC_SUCCESS, "failed to get meshes");
      fail_unless(temp_numMesh == 1, "should be 1 mesh");
      
      // check mesh
      rc = fc_getMeshName(meshes[0], &temp_name);
      fail_unless(rc == FC_SUCCESS, "should get mesh name");
      fail_unless(!strcmp(temp_name, mesh_name), "mesh name mismatch");
      free(temp_name);
      rc = fc_getMeshInfo(meshes[0], &temp_topodim, &temp_numDim, &temp_numVert,
			  &temp_numElem, &temp_elemType);
      fail_unless(rc == FC_SUCCESS, "failed to get mesh info");
      fail_unless(temp_topodim == topoDim, "mismatch of topodim");
      fail_unless(temp_numDim == numDim, "mismatch of dim");
      fail_unless(temp_numVert == numVert, "mismatch of numVert");
      fail_unless(temp_numElem == numElem, "mismatch of numElem");
      fail_unless(temp_elemType == elemType, "mismatch of elemType");
      rc = fc_getMeshCoordsPtr(meshes[0], &temp_coords);
      fail_unless(rc == FC_SUCCESS, "failed to get mesh coords");
      fail_unless(!memcmp(temp_coords, coords, numVert*numDim*sizeof(double)),
		  "mismatch of coords");
      rc = fc_getMeshElementConnsPtr(meshes[0], &temp_conns);
      fail_unless(rc == FC_SUCCESS, "failed to get mesh conns");
      fail_unless(!memcmp(temp_conns, conns, numElem*numVperE*sizeof(int)),
		  "mismatch of conns");

      // all done
      rc = fc_deleteDataset(dataset);
      fail_unless(rc == FC_SUCCESS, "failed to close after testing");
      free(mesh_name);
      free(meshes);
      free(coords);
      free(conns);
    }
  }
}
END_TEST

START_TEST(exo_multi_mesh)
{
  FC_ReturnCode rc;
  int i, j, k;
  FC_Dataset dataset;
  int temp_numMesh;
  FC_Mesh mesh, *meshes;
  char* ds_name = "dataset", *mesh_names[8];
  int numVerts[8], numElems[8], numDim, numVperEs[8];
  int* conns[8], *temp_conns;
  double* coords[8], *temp_coords;
  char* temp_name;
  int temp_topodim, temp_numDim, temp_numVert, temp_numElem;
  FC_ElementType elemTypes[8], temp_elemType;
  
  // Two cases: 1st has all types of meshes, all in 3D.
  //            2nd has mix of 1D, 2D & 3D meshes (will all end up 3d)

  // *** Case 1: all types of meshes, all in 3D

  // Create dataset
  rc = fc_createDataset(ds_name, &dataset);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create dataset");

  // Add the meshes
  for (i = 0; i < numMeshType; i++) {
    
    // Setup
    rc = _fc_readVertexFile(vertFiles[i], &numVerts[i], &numDim, &coords[i]);
    fail_unless(rc == FC_SUCCESS, "abort: failed to read vert file");
    fail_unless(numDim == 3, "abort: expect all to have numDim = 3");
    rc = _fc_readElementFile(elemFiles[i], &temp_numVert, &numElems[i], 
			     &elemTypes[i], &mesh_names[i], &conns[i]);
    fail_unless(rc == FC_SUCCESS, "abort: failed to read elem file");
    fail_unless(temp_numVert == numVerts[i], "abort: numVerts don't agree");
    numVperEs[i] = fc_getElementTypeNumVertex(elemTypes[i]);

    // Create the mesh
    rc = fc_createMesh(dataset, mesh_names[i], &mesh);
    fail_unless(rc == FC_SUCCESS, "abort: failed to create mesh");
    rc = fc_setMeshCoords(mesh, numDim, numVerts[i], coords[i]);
    fail_unless(rc == FC_SUCCESS, "abort: failed to set coords");
    rc = fc_setMeshElementConns(mesh, elemTypes[i], numElems[i], conns[i]);
    fail_unless(rc == FC_SUCCESS, "abort: failed to set conns");
  }

  // test writing
  rc = fc_writeDataset(dataset, EXODUS_FILE_NAME, FC_FT_EXODUS);
  fail_unless(rc == FC_SUCCESS, "should be able to write dataset");
    
  // test close & reopen
  rc = fc_deleteDataset(dataset);
  fail_unless(rc == FC_SUCCESS, "should close dataset w/o error");
  rc = fc_loadDatasetWithFormat(EXODUS_FILE_NAME, FC_FT_EXODUS, &dataset);
  fail_unless(rc == FC_SUCCESS, "should load dataset");
  
  // check dataset
  rc = fc_getMeshes(dataset, &temp_numMesh, &meshes);
  fail_unless(rc == FC_SUCCESS, "failed to get meshes");
  fail_unless(temp_numMesh == numMeshType, "mismatch of numMesh");
  
  // check the meshes
  for (i = 0; i < numMeshType; i++) {
    rc = fc_getMeshName(meshes[i], &temp_name);
    fail_unless(rc == FC_SUCCESS, "should get mesh name");
    fail_unless(!strcmp(temp_name, mesh_names[i]), "mesh name mismatch");
    free(temp_name);
    rc = fc_getMeshInfo(meshes[i], &temp_topodim, &temp_numDim, &temp_numVert,
			&temp_numElem, &temp_elemType);
    fail_unless(rc == FC_SUCCESS, "failed to get mesh info");
    fail_unless(temp_topodim == topoDims[i], "mismatch of topodim");
    fail_unless(temp_numDim == numDim, "mismatch of dim");
    fail_unless(temp_numVert == numVerts[i], "mismatch of numVert");
    fail_unless(temp_numElem == numElems[i], "mismatch of numElem");
    fail_unless(temp_elemType == elemTypes[i], "mismatch of elemType");
    rc = fc_getMeshCoordsPtr(meshes[i], &temp_coords);
    fail_unless(rc == FC_SUCCESS, "failed to get mesh coords");
    fail_unless(!memcmp(temp_coords, coords[i], numVerts[i]*numDim*sizeof(double)),
		"mismatch of coords");
    rc = fc_getMeshElementConnsPtr(meshes[i], &temp_conns);
    fail_unless(rc == FC_SUCCESS, "failed to get mesh conns");
    fail_unless(!memcmp(temp_conns, conns[i], numElems[i]*numVperEs[i]*sizeof(int)),
		"mismatch of conns");
  }


  // done with dataset, but keep the arrays of info
  rc = fc_deleteDataset(dataset);
  fail_unless(rc == FC_SUCCESS, "failed to close after testing");
  free(meshes);

  // *** Case 2: mix of 1D, 2D & 3D meshes (will all end up 3d)

  // Create dataset
  rc = fc_createDataset(ds_name, &dataset);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create dataset");

  // Add the meshes
  for (i = 0; i < numMeshType; i++) {
    
    // Setup - truncate point mesh to numDim = 1, and line, tri & quad
    // meshes to numDim = 2
    if (i == 0) 
      temp_numDim = 1;
    else if (i > 0 && i < 4)
      temp_numDim = 2;
    else
      temp_numDim = 3;
    if (temp_numDim < numDim) { // numDim = 3!
      temp_coords = malloc(numVerts[i]*temp_numDim*sizeof(double));
      for (j = 0; j < numVerts[i]; j++)
	memcpy(&temp_coords[j*temp_numDim], &coords[i][j*numDim],
	       temp_numDim*sizeof(double));
      free(coords[i]);
      coords[i] = temp_coords;
    }

    // Create the mesh
    rc = fc_createMesh(dataset, mesh_names[i], &mesh);
    fail_unless(rc == FC_SUCCESS, "abort: failed to create mesh");
    rc = fc_setMeshCoords(mesh, temp_numDim, numVerts[i], coords[i]);
    fail_unless(rc == FC_SUCCESS, "abort: failed to set coords");
    rc = fc_setMeshElementConns(mesh, elemTypes[i], numElems[i], conns[i]);
    fail_unless(rc == FC_SUCCESS, "abort: failed to set conns");
  }

  // test writing
  rc = fc_writeDataset(dataset, EXODUS_FILE_NAME, FC_FT_EXODUS);
  fail_unless(rc == FC_SUCCESS, "should be able to write dataset");
    
  // test close & reopen
  rc = fc_deleteDataset(dataset);
  fail_unless(rc == FC_SUCCESS, "should close dataset w/o error");
  rc = fc_loadDatasetWithFormat(EXODUS_FILE_NAME, FC_FT_EXODUS, &dataset);
  fail_unless(rc == FC_SUCCESS, "should load dataset");
  
  // check dataset
  rc = fc_getMeshes(dataset, &temp_numMesh, &meshes);
  fail_unless(rc == FC_SUCCESS, "failed to get meshes");
  fail_unless(temp_numMesh == numMeshType, "mismatch of numMesh");
  
  // check the meshes
  for (i = 0; i < numMeshType; i++) {
    rc = fc_getMeshName(meshes[i], &temp_name);
    fail_unless(rc == FC_SUCCESS, "should get mesh name");
    fail_unless(!strcmp(temp_name, mesh_names[i]), "mesh name mismatch");
    free(temp_name);
    rc = fc_getMeshInfo(meshes[i], &temp_topodim, &temp_numDim, &temp_numVert,
			&temp_numElem, &temp_elemType);
    fail_unless(rc == FC_SUCCESS, "failed to get mesh info");
    fail_unless(temp_topodim == topoDims[i], "mismatch of topodim");
    fail_unless(temp_numDim == 3, "mismatch of dim");
    fail_unless(temp_numVert == numVerts[i], "mismatch of numVert");
    fail_unless(temp_numElem == numElems[i], "mismatch of numElem");
    fail_unless(temp_elemType == elemTypes[i], "mismatch of elemType");
    rc = fc_getMeshCoordsPtr(meshes[i], &temp_coords);
    fail_unless(rc == FC_SUCCESS, "failed to get mesh coords");
    if (i == 0) {
      for (j = 0; j < numVerts[i]; j++) {
	fail_unless(temp_coords[j*numDim] == coords[i][j], 
		    "mismatch of coords");
	for (k = 1; k < numDim; k++)
	  fail_unless(temp_coords[j*numDim+k] == 0., "mismatch of coords");
      }
    }
    else if (i > 0 && i < 4) {
      for (j = 0; j < numVerts[i]; j++) {
	fail_unless(!memcmp(&temp_coords[j*numDim], &coords[i][j*2], 
			    2*sizeof(double)), "mismatch of coords");
	fail_unless(temp_coords[j*numDim+2] == 0., "mismatch of coords");
      }
    }
    else
      fail_unless(!memcmp(temp_coords, coords[i], numVerts[i]*numDim*sizeof(double)),
    	  "mismatch of coords");
    rc = fc_getMeshElementConnsPtr(meshes[i], &temp_conns);
    fail_unless(rc == FC_SUCCESS, "failed to get mesh conns");
    fail_unless(!memcmp(temp_conns, conns[i], numElems[i]*numVperEs[i]*sizeof(int)),
		"mismatch of conns");
  }

  // all done
  rc = fc_deleteDataset(dataset);
  fail_unless(rc == FC_SUCCESS, "failed to close after testing");
  free(meshes);
  for (i = 0; i < numMeshType; i++) {
    free(mesh_names[i]);
    free(coords[i]);
    free(conns[i]);
  }
}
END_TEST

START_TEST(exo_node_set)
{
  FC_ReturnCode rc;
  int i, j, k, m, mm;
  int temp_numMesh, temp_numSubset, temp_numMember, *temp_members;
  FC_Dataset dataset;
  FC_Mesh mesh, *temp_meshes;
  FC_Subset subset, *temp_subsets;
  FC_AssociationType temp_assoc;
  char* ds_name = "dataset";
  char* sub_names[3] = { "First 5 verts", "random verts", "empty" };
  char* temp_name;
  int numMembers[3] = { 5, 7, 0 };
  int members[3][7] = { { 0, 1, 2, 3, 4 }, { 2, 7, 11, 13, 21, 23, 25 }, {}  };

  // Test writing single mesh or two meshes with 1 or two node sets
  for (i = 1; i <= 2; i++) { // number of meshes
    for (j = 1; j <= 3; j++) { // number of nodeset per mesh
      
      // Create dataset, meshes & subset
      rc = fc_createDataset(ds_name, &dataset);
      fail_unless(rc == FC_SUCCESS, "abort: failed to create dataset");
      for (k = 0; k < i; k++) { 
	rc = _fc_makeMeshFromFiles(dataset, vertFiles[2+k], elemFiles[2+k], 
				   &mesh);
	fail_unless(rc == FC_SUCCESS, "abort: failed to create mesh");
	for (m = 0; m < j; m++) {
	  rc = fc_createSubset(mesh, sub_names[m], FC_AT_VERTEX, &subset);
	  fail_unless(rc == FC_SUCCESS, "abort: failed to create subset");
	  if (numMembers[m] > 0){
	    rc = fc_addArrayMembersToSubset(subset, numMembers[m], members[m]);
	    fail_unless(rc == FC_SUCCESS, 
			"abort: failed to add members to subset");
	  }
	}
      }

      // test writing
      rc = fc_writeDataset(dataset, EXODUS_FILE_NAME, FC_FT_EXODUS);
      fail_unless(rc == FC_SUCCESS, "should be able to write dataset");

      // test close & reopen
      rc = fc_deleteDataset(dataset);
      fail_unless(rc == FC_SUCCESS, "should close dataset w/o error");
      rc = fc_loadDatasetWithFormat(EXODUS_FILE_NAME, FC_FT_EXODUS, &dataset);
      fail_unless(rc == FC_SUCCESS, "should load dataset");

      // check dataset, mesh & subsets
      rc = fc_getMeshes(dataset, &temp_numMesh, &temp_meshes);
      fail_unless(temp_numMesh == i, "mismatch of numMesh");
      for (k = 0; k < i; k++) {
	int found[3] = {0,0,0};
	rc = fc_getSubsets(temp_meshes[k], &temp_numSubset, &temp_subsets);
	if (j == 3){
	  //empty subsets put on all meshes and not merged
	  fail_unless(temp_numSubset = (j-1)+temp_numMesh, "mismatch of num subset");
	} else {
	  fail_unless(temp_numSubset = j, "mismatch of num subset");
	}
	// check subsets
	for (m = 0; m < temp_numSubset; m++) {
	  subset = temp_subsets[m];
	  rc = fc_getSubsetName(subset, &temp_name);
	  fail_unless(rc == FC_SUCCESS, "should get subset name");
	  for (mm = 0; mm < 3; mm++){
	    if (!strcmp(temp_name, sub_names[mm])){
	      found[mm]++;
	      rc = fc_getSubsetInfo(subset, &temp_numMember, NULL, &temp_assoc);
	      fail_unless(rc == FC_SUCCESS, "failed to get subset info");
	      fail_unless(temp_numMember == numMembers[mm], 
			  "mismatch of numMember");
	      fail_unless(temp_assoc == FC_AT_VERTEX, "mismatch of assoc");
	      if (temp_numMember > 0){
		rc = fc_getSubsetMembersAsArray(subset, &temp_numMember,
						&temp_members);
		fail_unless(rc == FC_SUCCESS, "failed to get members");
		fail_unless(!memcmp(temp_members, members[mm],
				    numMembers[mm]*sizeof(int)), "mismatch of members");
		free(temp_members);
	      }
	    }
	  }
	  free(temp_name);
	}
	for (m = 0; m < j &&  j!= 3; m++){
	  fail_unless(found[m] == 1, "mismatch of meshes");
	}
	for (m = j; m < 2; m++)
	  fail_unless(found[m] == 0, "mismatch of meshes");
	if (j == 3)
	  fail_unless(found[2] == temp_numMesh, "mismatch of meshes");
	else 
	  fail_unless(found[2] == 0, "mismatch of meshes");
	free(temp_subsets);
      }
      free(temp_meshes);

      // all done
      rc = fc_deleteDataset(dataset);
      fail_unless(rc == FC_SUCCESS, "failed to close after testing");
    }
  }
}
END_TEST


START_TEST(exo_seq_node_set)
{
  FC_ReturnCode rc;
  int i, j, k, m, ii, jj;
  int temp_numMesh, temp_numSubset, temp_numMember, *temp_members;
  FC_Dataset dataset;
  FC_Mesh mesh, *temp_meshes;
  FC_Subset subset, *temp_subsets, *seqSubset;
  FC_AssociationType temp_assoc;
  FC_Sequence sequence;
  int numStep = 5;
  char* ds_name = "dataset";
  char* sub_names[2] = { "First 5 verts", "random verts" };
  char* seq_name = "sequence";
  char* temp_name;
  int numMembers[2] = { 5, 7 };
  int members[2][7] = { { 0, 1, 2, 3, 4 }, { 2, 7, 11, 13, 21, 23, 25 } };
  int data[5];
  int mask[2][6];

  for (i = 0; i < 5; i++) 
    data[i] = i/10.;
    

  // Test writing single mesh or two meshes with 1 or two node sets
  for (i = 1; i <= 2; i++) { // number of meshes
    for (j = 1; j <= 2; j++) { // number of nodeset per mesh
      
      // Create dataset, meshes & subset
      rc = fc_createDataset(ds_name, &dataset);
      fail_unless(rc == FC_SUCCESS, "abort: failed to create dataset");
      // sequence
      rc = fc_createSequence(dataset, seq_name, &sequence);
      fail_unless(rc == FC_SUCCESS, "abort: failed to create sequence");
      rc = fc_setSequenceCoords(sequence, numStep, FC_DT_FLOAT, data);
      fail_unless(rc == FC_SUCCESS, "abort: failed to set sequence coords");

      for (k = 0; k < i; k++) { 
	rc = _fc_makeMeshFromFiles(dataset, vertFiles[2+k], elemFiles[2+k], 
				   &mesh);
	fail_unless(rc == FC_SUCCESS, "abort: failed to create mesh");

	temp_subsets = (FC_Subset*)malloc(numStep*sizeof(FC_Subset));
	for (m = 0; m < j; m++) {
	  rc = fc_createSubset(mesh, sub_names[m], FC_AT_VERTEX, &subset);
	  fail_unless(rc == FC_SUCCESS, "abort: failed to create subset");
	  rc = fc_addArrayMembersToSubset(subset, numMembers[m], members[m]);
	  fail_unless(rc == FC_SUCCESS, 
		      "abort: failed to add members to subset");
	  for (ii = 0; ii < numStep; ii++){
	    rc = fc_copySubset(subset, mesh, "junk", &(temp_subsets[ii]));
	    fail_unless(rc == FC_SUCCESS, 
			"abort: failed to copy subset");
	  }
	  rc = fc_convertSubsetsToSeqSubset(numStep, temp_subsets, sequence, sub_names[m], 
					    &seqSubset);
	  fail_unless(rc == FC_SUCCESS, 
		      "abort: failed to create seq subset");
	  free(seqSubset);
	}
	free(temp_subsets);
      }

      // test writing
      rc = fc_writeDataset(dataset, EXODUS_FILE_NAME, FC_FT_EXODUS);
      fail_unless(rc == FC_SUCCESS, "should be able to write dataset");

      // test close & reopen
      rc = fc_deleteDataset(dataset);
      fail_unless(rc == FC_SUCCESS, "should close dataset w/o error");
      rc = fc_loadDatasetWithFormat(EXODUS_FILE_NAME, FC_FT_EXODUS, &dataset);
      fail_unless(rc == FC_SUCCESS, "should load dataset");

      // check dataset, mesh & subsets
      rc = fc_getMeshes(dataset, &temp_numMesh, &temp_meshes);
      fail_unless(temp_numMesh == i, "mismatch of numMesh");
      for (k = 0; k < i; k++) {
	for (m = 0; m < 2; m++){
	  for (ii = 0; ii < 6; ii++){
	    mask[m][ii] = 0;
	  }
	}

	rc = fc_getSubsets(temp_meshes[k], &temp_numSubset, &temp_subsets);
	fail_unless(temp_numSubset == j*numStep+j, "mismatch of numSubset");
	// check subsets
	for (m = 0; m < temp_numSubset; m++) {
	  int found = 0;
	  subset = temp_subsets[m];
	  rc = fc_getSubsetName(subset, &temp_name);
	  fail_unless(rc == FC_SUCCESS, "should get subset name");
	  //we should have the one orig with no numbering and then numStep with numbering per j
	  for (ii = 0; ii < j; ii++){
	    if (!strncmp(temp_name, sub_names[ii], strlen(sub_names[ii]))){
	      for (jj = 0; jj < numStep+1; jj++){
		char* comp_name = (char*)malloc((strlen(sub_names[ii])+40)*sizeof(char));
		if (jj < numStep){
		  strcpy(comp_name,sub_names[ii]);
		} else {
		  snprintf(comp_name, strlen(sub_names[ii])+40,"%s_%d", sub_names[ii],jj);
		}
		if (!strcmp(temp_name,comp_name)){
		  fail_unless(mask[ii][jj] == 0, "already found subset");
		  mask[ii][jj] = 1;
		  rc = fc_getSubsetInfo(subset, &temp_numMember, NULL, &temp_assoc);
		  fail_unless(rc == FC_SUCCESS, "failed to get subset info");
		  fail_unless(temp_numMember == numMembers[ii], 
			      "mismatch of numMember");
		  fail_unless(temp_assoc == FC_AT_VERTEX, "mismatch of assoc");
		  rc = fc_getSubsetMembersAsArray(subset, &temp_numMember,
						  &temp_members);
		  fail_unless(rc == FC_SUCCESS, "failed to get members");
		  fail_unless(!memcmp(temp_members, members[ii],
				      numMembers[ii]*sizeof(int)), "mismatch of members"); 
		  free(temp_members);
		  free(comp_name);
		  found = 1;
		  break;
		}
		free(comp_name);
	      } //jj
	    } //if
	    if (found) break;
	  } //ii
	  free(temp_name);
	}
	free(temp_subsets);
      }
      free(temp_meshes);
    }
  }


  // all done
  rc = fc_deleteDataset(dataset);
  fail_unless(rc == FC_SUCCESS, "failed to close after testing");
}
END_TEST


//checks both side and elem sets
START_TEST(exo_elem_side_set)
{
  FC_ReturnCode rc;
  int i, j, k;
  int temp_numMesh, temp_numSubset, temp_numMember, *temp_members;
  FC_Dataset dataset;
  FC_Mesh mesh, *temp_meshes;
  FC_Subset subset, *temp_subsets;
  char* ds_name = "dataset";
  FC_AssociationType assoc, temp_assoc;
  char* sub_names[4] = { "First 5 'sides'", "random 'sides'" , "emptyside", "emptyelem"};
  char* temp_name;
  int numMembers[4] = { 5, 7, 0 };
  //note these numbers also work for elem sets as well
  int members[4][7] = { { 0, 1, 2, 3, 4 }, { 2, 7, 11, 13, 21, 23, 25 }, {} };
  int worked[8][5] = { { 1, -1, -1, 1, 0 },   // point
		       { 1, -1, -1, 1, 0 },   // line
		       { 1,  1, -1, 1, 0 },   // tri
		       { 1,  1, -1, 1, 0 },   // quad
		       { 1,  0,  1, 1, 0 },   // tet
		       { 1,  0,  1, 1, 0 },   // pyramid
		       { 1,  0,  1, 1, 0 },   // prism
		       { 1,  0,  1, 1, 0 } }; // hex
  
  // ** Case 1 -- Single mesh single side set - test every kind sideset/mesh combo

  fflush(NULL);
  for (i = 0; i < numMeshType; i++) {
    for (j = 1; j < numAssocType; j++) { // skip vertex, done in prv test

      // skip inappropriate combinatations 
      if ( ((i == 0 || i == 1) && (j == 1 || j == 2)) || 
         ((i == 2 || i == 3) && j == 2) )
	continue;
					 
      // Create dataset, mesh & subset
      rc = fc_createDataset(ds_name, &dataset);
      fail_unless(rc == FC_SUCCESS, "abort: failed to create dataset");
      rc = _fc_makeMeshFromFiles(dataset, vertFiles[i], elemFiles[i], 
				 &mesh);
      fail_unless(rc == FC_SUCCESS, "abort: failed to create mesh");
      //note that this also works for elem set as well, even though still
      //using the side names and numbers for the data
      rc = fc_createSubset(mesh, sub_names[0], assocTypes[j], &subset);
      fail_unless(rc == FC_SUCCESS, "abort: failed to create subset");
      if (j < 4)
	rc = fc_addArrayMembersToSubset(subset, numMembers[0], members[0]);
      else
	rc = fc_addMemberToSubset(subset, 0);
      fail_unless(rc == FC_SUCCESS, 
		  "abort: failed to add members to subset");
      
      // test writing
      rc = fc_writeDataset(dataset, EXODUS_FILE_NAME, FC_FT_EXODUS);
      fail_unless(rc == FC_SUCCESS, "should be able to write dataset");
      
      // test close & reopen
      rc = fc_deleteDataset(dataset);
      fail_unless(rc == FC_SUCCESS, "should close dataset w/o error");
      
      rc = fc_loadDatasetWithFormat(EXODUS_FILE_NAME, FC_FT_EXODUS, &dataset);
      fail_unless(rc == FC_SUCCESS, "should load dataset");

      // check subsets
      rc = fc_getMeshes(dataset, &temp_numMesh, &temp_meshes);
      fail_unless(temp_numMesh == 1, "mismatch of numMesh");
      rc = fc_getSubsets(temp_meshes[0], &temp_numSubset, &temp_subsets);
      fail_unless(temp_numSubset == worked[i][j], "mismatch of numSubset");
      if (temp_subsets) {
	rc = fc_getSubsetName(temp_subsets[0], &temp_name);
	fail_unless(rc == FC_SUCCESS, "should get subset name");
	fail_unless(!strcmp(temp_name, sub_names[0]), "subset name mismatch");
	free(temp_name);
	rc = fc_getSubsetInfo(temp_subsets[0], &temp_numMember, NULL, &temp_assoc);
	fail_unless(rc == FC_SUCCESS, "failed to get subset info");
	fail_unless(temp_numMember == numMembers[0], 
		    "mismatch of numMember");
	fail_unless(temp_assoc == assocTypes[j], "mismatch of assoc");
	rc = fc_getSubsetMembersAsArray(temp_subsets[0], &temp_numMember,
					&temp_members);
	fail_unless(rc == FC_SUCCESS, "failed to get members");
	fail_unless(!memcmp(temp_members, members[0],
			    numMembers[0]*sizeof(int)), "mismatch of members");
	free(temp_members);
      }

      // all done
      free(temp_meshes);
      free(temp_subsets);
      rc = fc_deleteDataset(dataset);
      fail_unless(rc == FC_SUCCESS, "failed to close after testing");
      
    }
  }

  fflush(NULL);

  // ** Case 2 -- Multiple meshes & multiple side sets (edges on 2D meshes) and empty side and elem set

  // Create dataset, meshes & subset
  rc = fc_createDataset(ds_name, &dataset);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create dataset");
  assoc = FC_AT_EDGE;
  for (i = 2; i < 4; i++) { // tri & quad
    rc = _fc_makeMeshFromFiles(dataset, vertFiles[i], elemFiles[i],
			       &mesh);
    fail_unless(rc == FC_SUCCESS, "abort: failed to create mesh");
    for (j = 0; j < 4; j++) {
      rc = fc_createSubset(mesh, sub_names[j], assoc, &subset);
      fail_unless(rc == FC_SUCCESS, "abort: failed to create subset");
      if (numMembers[j] > 0){
	rc = fc_addArrayMembersToSubset(subset, numMembers[j], members[j]);
	fail_unless(rc == FC_SUCCESS, 
		    "abort: failed to add members to subset");
      }
    }
  }
 
  // test writing
  rc = fc_writeDataset(dataset, EXODUS_FILE_NAME, FC_FT_EXODUS);
  fail_unless(rc == FC_SUCCESS, "should be able to write dataset");

  // test close & reopen
  rc = fc_deleteDataset(dataset);
  fail_unless(rc == FC_SUCCESS, "should close dataset w/o error");
  rc = fc_loadDatasetWithFormat(EXODUS_FILE_NAME, FC_FT_EXODUS, &dataset);
  fail_unless(rc == FC_SUCCESS, "should load dataset");

  // check subsets
  rc = fc_getMeshes(dataset, &temp_numMesh, &temp_meshes);
  fail_unless(temp_numMesh == 2, "mismatch of numMesh");
  for (i = 0; i < 2; i++) {
    int found[4] = {0,0,0,0};
    rc = fc_getSubsets(temp_meshes[i], &temp_numSubset, &temp_subsets);
    //empty side and elem sets will be put on all meshes
    fail_unless(temp_numSubset == 6, "mismatch of numSubset");
    for (j = 0; j < 6; j++) {
      rc = fc_getSubsetName(temp_subsets[j], &temp_name);
      fail_unless(rc == FC_SUCCESS, "should get subset name");
      //get nummesh*nummesh versions of empty subset since they was
      //initially put on all meshes and then will be copied onto all meshes
      //as part of the load
      for (k = 0; k < 4; k++){
	if (!strcmp(temp_name, sub_names[k])){
	  rc = fc_getSubsetInfo(temp_subsets[j], &temp_numMember, NULL, &temp_assoc);
	  fail_unless(rc == FC_SUCCESS, "failed to get subset info");
	  fail_unless(temp_numMember == numMembers[k], 
		      "mismatch of numMember");
	  fail_unless(temp_assoc == assoc, "mismatch of assoc");
	  if (temp_numMember > 0){
	    rc = fc_getSubsetMembersAsArray(temp_subsets[j], &temp_numMember,
					    &temp_members);
	    fail_unless(rc == FC_SUCCESS, "failed to get members");
	    fail_unless(!memcmp(temp_members, members[k],
			    numMembers[k]*sizeof(int)), "mismatch of members");
	    free(temp_members);
	  }
	  found[k]++;
	}
      }
      free(temp_name);
    }
    fail_unless(found[0] == 1, "wrong num copies of subset");
    fail_unless(found[1] == 1, "wrong num copies of subset");
    fail_unless(found[2] == 2, "wrong num copies of empty subset");
    fail_unless(found[3] == 2, "wrong num copies of empty subset");
    free(temp_subsets);
  }
  free(temp_meshes);

  // all done
  rc = fc_deleteDataset(dataset);
  fail_unless(rc == FC_SUCCESS, "failed to close after testing");
}
END_TEST


//checks both side and elem sets
START_TEST(exo_seq_elem_side_set)
{
  FC_ReturnCode rc;
  int i, j, m, ii, jj, mm;
  int temp_numMesh, temp_numSubset, temp_numMember, *temp_members;
  FC_Dataset dataset;
  FC_Mesh mesh, *temp_meshes;
  FC_Subset subset, *temp_subsets, *seqSubset;
  char* ds_name = "dataset";
  FC_AssociationType  temp_assoc, assoc;
  char* sub_names[2] = { "First 5 'sides'", "random 'sides'" };
  char* temp_name;
  int numMembers[2] = { 5, 7 };
  FC_Sequence sequence;
  int numStep = 5;
  char* seq_name = "sequence";
  int data[5];
  int mask[2][6];

  //note these numbers also work for elem sets as well
  int members[2][7] = { { 0, 1, 2, 3, 4 }, { 2, 7, 11, 13, 21, 23, 25 } };
  int worked[8][5] = { { 1, -1, -1, 1, 0 },   // point
		       { 1, -1, -1, 1, 0 },   // line
		       { 1,  1, -1, 1, 0 },   // tri
		       { 1,  1, -1, 1, 0 },   // quad
		       { 1,  0,  1, 1, 0 },   // tet
		       { 1,  0,  1, 1, 0 },   // pyramid
		       { 1,  0,  1, 1, 0 },   // prism
		       { 1,  0,  1, 1, 0 } }; // hex
  
  // ** Case 1 -- Single mesh single side set - test every kind sideset/mesh combo

  fflush(NULL);

  for (i = 0; i < 5; i++) 
    data[i] = i/10.;

  for (i = 0; i < numMeshType; i++) {
    for (j = 1; j < numAssocType; j++) { // skip vertex, done in prv test

      // skip inappropriate combinatations 
      if ( ((i == 0 || i == 1) && (j == 1 || j == 2)) || 
         ((i == 2 || i == 3) && j == 2) )
	continue;

      // Create dataset, mesh & subset
      rc = fc_createDataset(ds_name, &dataset);
      fail_unless(rc == FC_SUCCESS, "abort: failed to create dataset");
      // sequence
      rc = fc_createSequence(dataset, seq_name, &sequence);
      fail_unless(rc == FC_SUCCESS, "abort: failed to create sequence");
      rc = fc_setSequenceCoords(sequence, numStep, FC_DT_FLOAT, data);
      fail_unless(rc == FC_SUCCESS, "abort: failed to set sequence coords");

      rc = _fc_makeMeshFromFiles(dataset, vertFiles[i], elemFiles[i], 
				 &mesh);
      fail_unless(rc == FC_SUCCESS, "abort: failed to create mesh");
      //note that this also works for elem set as well, even though still
      //using the side names and numbers for the data
      rc = fc_createSubset(mesh, sub_names[0], assocTypes[j], &subset);
      fail_unless(rc == FC_SUCCESS, "abort: failed to create subset");
      if (j < 4)
	rc = fc_addArrayMembersToSubset(subset, numMembers[0], members[0]);
      else
	rc = fc_addMemberToSubset(subset, 0);
      fail_unless(rc == FC_SUCCESS, 
		  "abort: failed to add members to subset");
      temp_subsets = (FC_Subset*)malloc(numStep*sizeof(FC_Subset));
      for (ii = 0; ii < numStep; ii++){
	rc = fc_copySubset(subset, mesh, "junk", &(temp_subsets[ii]));
	fail_unless(rc == FC_SUCCESS, 
		    "abort: failed to copy subset");
      }
      rc = fc_convertSubsetsToSeqSubset(numStep, temp_subsets, sequence, sub_names[0], 
					&seqSubset);
      fail_unless(rc == FC_SUCCESS, 
		  "abort: failed to create seq subset");
      free(seqSubset);
      free(temp_subsets);

      // test writing
      rc = fc_writeDataset(dataset, EXODUS_FILE_NAME, FC_FT_EXODUS);
      fail_unless(rc == FC_SUCCESS, "should be able to write dataset");
      
      // test close & reopen
      rc = fc_deleteDataset(dataset);
      fail_unless(rc == FC_SUCCESS, "should close dataset w/o error");
      
      rc = fc_loadDatasetWithFormat(EXODUS_FILE_NAME, FC_FT_EXODUS, &dataset);
      fail_unless(rc == FC_SUCCESS, "should load dataset");

      // check subsets
      rc = fc_getMeshes(dataset, &temp_numMesh, &temp_meshes);
      fail_unless(temp_numMesh == 1, "mismatch of numMesh");
      rc = fc_getSubsets(temp_meshes[0], &temp_numSubset, &temp_subsets);
      fail_unless(temp_numSubset == worked[i][j]*numStep+worked[i][j], "mismatch of numSubset");
      if (temp_subsets) {
	// check subsets
	for (mm = 0; mm < 2; mm++){
	  for (m = 0; m < numStep+1; m++){
	    mask[mm][m] = 0;
	  }
	}
	for (m = 0; m < temp_numSubset; m++) {
	  subset = temp_subsets[m];
	  rc = fc_getSubsetName(subset, &temp_name);
	  fail_unless(rc == FC_SUCCESS, "should get subset name");
	  //we should have the one orig with no numbering and then numStep with numbering
	  if (!strncmp(temp_name, sub_names[0], strlen(sub_names[0]))){
	    for (jj = 0; jj < numStep+1; jj++){
	      char* comp_name = (char*)malloc((strlen(sub_names[0])+40)*sizeof(char));
	      if (jj < numStep){
		strcpy(comp_name,sub_names[0]);
	      } else {
		snprintf(comp_name, strlen(sub_names[0])+40,"%s_%d", sub_names[0],jj);
	      }
	      if (!strcmp(temp_name,comp_name)){
		fail_unless(mask[0][jj] == 0, "already found subset");
		mask[0][jj] = 1;
		rc = fc_getSubsetInfo(subset, &temp_numMember, NULL, &temp_assoc);
		fail_unless(rc == FC_SUCCESS, "failed to get subset info");
		fail_unless(temp_assoc == assocTypes[j], "mismatch of assoc");
		fail_unless(temp_numMember == numMembers[0], 
			      "mismatch of numMember");
		rc = fc_getSubsetMembersAsArray(subset, &temp_numMember,
						&temp_members);
		fail_unless(rc == FC_SUCCESS, "failed to get members");
		fail_unless(!memcmp(temp_members, members[0],
				    numMembers[0]*sizeof(int)), "mismatch of members"); 
		free(temp_members);
		free(comp_name); 
		break;
	      }
	      free(comp_name); 
	    } //jj
	  } //if
	  free(temp_name);
	}
	free(temp_subsets);
      }

      // all done
      free(temp_meshes);
      rc = fc_deleteDataset(dataset);
      fail_unless(rc == FC_SUCCESS, "failed to close after testing");
      
    }
  }

  fflush(NULL);

  // ** Case 2 -- Multiple meshes & multiple side sets (edges on 2D meshes)

  // Create dataset, meshes & subset
  rc = fc_createDataset(ds_name, &dataset);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create dataset");
  rc = fc_createSequence(dataset, seq_name, &sequence);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create sequence");
  rc = fc_setSequenceCoords(sequence, numStep, FC_DT_FLOAT, data);
  fail_unless(rc == FC_SUCCESS, "abort: failed to set sequence coords");

  assoc = FC_AT_EDGE;
  for (i = 2; i < 4; i++) { // tri & quad
    rc = _fc_makeMeshFromFiles(dataset, vertFiles[i], elemFiles[i],
			       &mesh);
    fail_unless(rc == FC_SUCCESS, "abort: failed to create mesh");
    for (j = 0; j < 2; j++) {
      rc = fc_createSubset(mesh, sub_names[j], assoc, &subset);
      fail_unless(rc == FC_SUCCESS, "abort: failed to create subset");
      rc = fc_addArrayMembersToSubset(subset, numMembers[j], members[j]);
      fail_unless(rc == FC_SUCCESS, 
		  "abort: failed to add members to subset");
      temp_subsets = (FC_Subset*)malloc(numStep*sizeof(FC_Subset));
      for (ii = 0; ii < numStep; ii++){
	rc = fc_copySubset(subset, mesh, "junk", &(temp_subsets[ii]));
	fail_unless(rc == FC_SUCCESS, 
		    "abort: failed to copy subset");
      }
      rc = fc_convertSubsetsToSeqSubset(numStep, temp_subsets, sequence, sub_names[0], 
					&seqSubset);
      fail_unless(rc == FC_SUCCESS, 
		  "abort: failed to create seq subset");
      free(seqSubset);
      free(temp_subsets);
    }
  }
 
  // test writing
  rc = fc_writeDataset(dataset, EXODUS_FILE_NAME, FC_FT_EXODUS);
  fail_unless(rc == FC_SUCCESS, "should be able to write dataset");

  // test close & reopen
  rc = fc_deleteDataset(dataset);
  fail_unless(rc == FC_SUCCESS, "should close dataset w/o error");
  rc = fc_loadDatasetWithFormat(EXODUS_FILE_NAME, FC_FT_EXODUS, &dataset);
  fail_unless(rc == FC_SUCCESS, "should load dataset");

  // check subsets
  rc = fc_getMeshes(dataset, &temp_numMesh, &temp_meshes);
  fail_unless(temp_numMesh == 2, "mismatch of numMesh");
  for (i = 0; i < 2; i++) {
    rc = fc_getSubsets(temp_meshes[i], &temp_numSubset, &temp_subsets);
    fail_unless(temp_numSubset == 2*numStep+2, "mismatch of numSubset");
    if (temp_subsets) {
      // check subsets
      for (m = 0; m < 2; m++){
	for (mm = 0; mm < numStep+1; mm++){
	  mask[m][mm] = 0;
	}
      }
      for (m = 0; m < temp_numSubset; m++) {
	int found = 0;
	subset = temp_subsets[m];
	rc = fc_getSubsetName(subset, &temp_name);
	fail_unless(rc == FC_SUCCESS, "should get subset name");
	//we should have the one orig with no numbering and then numStep with numbering
	//for each orig subset
	for (mm = 0; mm < 2; mm++){
	  if (!strncmp(temp_name, sub_names[mm], strlen(sub_names[mm]))){
	    for (jj = 0; jj < numStep+1; jj++){
	      char* comp_name = (char*)malloc((strlen(sub_names[mm])+40)*sizeof(char));
	      if (jj < numStep){
		strcpy(comp_name,sub_names[mm]);
	      } else {
		snprintf(comp_name, strlen(sub_names[mm])+40,"%s_%d", sub_names[mm],jj);
	      }
	      if (!strcmp(temp_name,comp_name)){
		fail_unless(mask[mm][jj] == 0, "already found subset");
		mask[mm][jj] = 1;
		rc = fc_getSubsetInfo(subset, &temp_numMember, NULL, &temp_assoc);
		fail_unless(rc == FC_SUCCESS, "failed to get subset info");
		fail_unless(temp_assoc == assoc, "mismatch of assoc");
		fail_unless(temp_numMember == numMembers[mm], 
			    "mismatch of numMember");
		rc = fc_getSubsetMembersAsArray(subset, &temp_numMember,
						&temp_members);
		fail_unless(rc == FC_SUCCESS, "failed to get members");
		fail_unless(!memcmp(temp_members, members[mm],
				    numMembers[mm]*sizeof(int)), "mismatch of members"); 
		free(temp_members);
		free(comp_name); 
		found = 1;
		break;
	      }
	      free(comp_name); 
	    } //jj
	  } //if
	  if (found) break;
	} //mm
	free(temp_name);
      } //m
      free(temp_subsets);
    }
  } // i 

  free(temp_meshes);

  // all done
  rc = fc_deleteDataset(dataset);
  fail_unless(rc == FC_SUCCESS, "failed to close after testing");
}
END_TEST


// FC can't write global subsets, but may encounter them: So write an exo file
// (via exodus API) with 1 global nodeset & 1 global elem set & 1 global
// sideset & test read.  Also test that vertices and elems shared between meshes 
// are handles properly
// Also test that names w/o names file are as expected.
START_TEST(exo_global_subsets)
{
  FC_ReturnCode rc;
  int error, i, j, k;
  int exoid, cpu_word_size = 8, io_word_size = 8;
  int numDim = 3, numVert = 12, numElem = 2, numMesh = 2;
  int numNodeSet = 1, numSideSet = 1, numElemSet = 1, numMembers[3] = { 5, 2, 2 };
  int meshIDs[2] = { 5, 7 }, subsetIDs[3] = { 3, 11, 13 }; //where are these numbers from?
  int nodes[6] = { 1, 5, 6, 8, 12 }; // 1st & last (to catch off by 1 bug)
  int elems[2] = { 1, 2 }, sides[2] = { 1, 1 };
  int numElemPerBlock[2] = { 1, 1 }, numVertPerElem = 8;
  double coords[3][12];
  char *coord_names[3] = { "xcoord", "ycoord", "zcoord" };
  int conns[2][8] = { { 1, 2, 4, 3, 5, 6, 8, 7 }, 
		      { 5, 6, 8, 7, 9, 10, 12, 11 } };
  char str_buf[1024];
  FC_Dataset dataset;
  FC_Mesh *temp_meshes;
  FC_Subset *temp_subsets;
  int numSubset = 3;
  int temp_numMesh, temp_topodim, temp_numDim, temp_numVert, temp_numElem;
  int temp_numSubset, *temp_conns, temp_numMember, *temp_members;
  char* temp_name;
  FC_ElementType temp_elemType;
  FC_AssociationType temp_assoc;
  double* temp_coords;
  ex_init_params exoParams;

  // Setup - create coords
  for (k = 0; k < 3; k++) {
    for (j = 0; j < 2; j++) {
      for (i = 0; i < 2; i++) {
	int idx = k*2*2 + j*2 + i;
	coords[0][idx] = i;
	coords[1][idx] = j;
	coords[2][idx] = k;
      }
    }
  }

  // Create exodus file with 2 joined hexes, each it's own mesh
  // and 1 side set spanning each, 1 node set spanning each, 1 elem
  // set spanning each
  exoid = ex_create(EXODUS_FILE_NAME, EX_CLOBBER, &cpu_word_size,
		    &io_word_size); 
  fail_unless(exoid > 0, "abort: failed to created exodus file");
  strncpy(exoParams.title, "simple", MAX_STR_LENGTH);
  exoParams.num_dim = numDim;
  exoParams.num_nodes = numVert;
  exoParams.num_edge = 0;
  exoParams.num_edge_blk = 0;
  exoParams.num_face = 0;
  exoParams.num_face_blk = 0;
  exoParams.num_elem = numElem;
  exoParams.num_elem_blk = numMesh;
  exoParams.num_node_sets = numNodeSet;
  exoParams.num_edge_sets = 0;
  exoParams.num_face_sets = 0;
  exoParams.num_side_sets = numSideSet;
  exoParams.num_elem_sets = numElemSet;
  exoParams.num_node_maps = 0;
  exoParams.num_edge_maps = 0;
  exoParams.num_face_maps = 0;
  exoParams.num_elem_maps = 0;
  error = ex_put_init_ext(exoid, &exoParams);
  fail_unless(error == 0, "abort: failed to init file");
  error = ex_put_coord(exoid, coords[0], coords[1], coords[2]);
  fail_unless(error == 0, "abort: failed to put coords");
  error = ex_put_coord_names(exoid, coord_names);
  fail_unless(error == 0, "abort: failed to put coord names");
  for (i = 0; i < numMesh; i++) {
    error = ex_put_elem_block(exoid, meshIDs[i], "HEX", numElemPerBlock[i],
			      numVertPerElem, 0);
    fail_unless(error == 0, "abort: Failed to put element block");
  }
  for (i = 0; i < numMesh; i++) {
    error = ex_put_elem_conn(exoid, meshIDs[i], conns[i]);
    fail_unless(error == 0, "abort: Failed to put element block conns");
  }
  error = ex_put_node_set_param(exoid, subsetIDs[0], numMembers[0], 0);
  fail_unless(error == 0, "abort: Failed to put node set param");
  error = ex_put_node_set(exoid, subsetIDs[0], nodes);
  fail_unless(error == 0, "abort: Failed to put node set nodes");
  error = ex_put_side_set_param(exoid, subsetIDs[1], numMembers[1], 0);
  fail_unless(error == 0, "abort: Failed to put side set param");
  error = ex_put_side_set(exoid, subsetIDs[1], elems, sides);
  fail_unless(error == 0, "abort: Failed to put side set sides");
  error = ex_put_set_param(exoid, EX_ELEM_SET, subsetIDs[2], numMembers[2], 0);
  fail_unless(error == 0, "abort: Failed to put elem set param");
  error = ex_put_set(exoid, EX_ELEM_SET, subsetIDs[2], elems, NULL);
  fail_unless(error == 0, "abort: Failed to put elem set param");
  error = ex_close(exoid);
  fail_unless(error == 0, "abort: Failed to close exodus file");
  
  // Read the dataset
  rc = fc_loadDatasetWithFormat(EXODUS_FILE_NAME, FC_FT_EXODUS, &dataset);
  fail_unless(rc == FC_SUCCESS, "should load dataset");

  // check the mesh & subsets
  rc = fc_getMeshes(dataset, &temp_numMesh, &temp_meshes);
  fail_unless(rc == FC_SUCCESS, "failed to get meshes");
  fail_unless(temp_numMesh == numMesh, "mismatch of numMesh");
  for (i = 0; i < numMesh; i++) {
    int vertOffset = 0;
    if (i == 1)
      vertOffset = 4;

    // check mesh
    rc = fc_getMeshName(temp_meshes[i], &temp_name);
    fail_unless(rc == FC_SUCCESS, "should get mesh name");
    sprintf(str_buf, "block_%d", meshIDs[i]);
    fail_unless(!strcmp(temp_name, str_buf), "mesh name mismatch");
    free(temp_name);
    rc = fc_getMeshInfo(temp_meshes[i], &temp_topodim, &temp_numDim, 
			&temp_numVert, &temp_numElem, &temp_elemType);
    fail_unless(rc == FC_SUCCESS, "failed to get mesh info");
    fail_unless(temp_topodim == numDim, "mismatch of topodim");
    fail_unless(temp_numDim == numDim, "mismatch of dim");
    fail_unless(temp_numVert == numVertPerElem, "mismatch of numVert");
    fail_unless(temp_numElem == numElemPerBlock[i], "mismatch of numElem");
    fail_unless(temp_elemType == FC_ET_HEX, "mismatch of elemType");
    rc = fc_getMeshCoordsPtr(temp_meshes[i], &temp_coords);
    fail_unless(rc == FC_SUCCESS, "failed to get mesh coords");
    for (j = 0; j < temp_numVert; j++)
      for (k = 0; k < numDim; k++)
	fail_unless(temp_coords[j*numDim+k] == coords[k][j+vertOffset], 
		    "mismatch of coords");
    rc = fc_getMeshElementConnsPtr(temp_meshes[i], &temp_conns);
    fail_unless(rc == FC_SUCCESS, "failed to get mesh conns");
    for (j = 0; j < temp_numVert; j++)
      fail_unless(temp_conns[j] == conns[0][j] - 1, "mismatch of conns");

    // check subsets
    rc = fc_getSubsets(temp_meshes[i], &temp_numSubset, &temp_subsets);
    fail_unless(temp_numSubset == numSubset, "mismatch of numSubset");
    for (j = 0; j < numSubset; j++) {
      rc = fc_getSubsetName(temp_subsets[j], &temp_name);
      fail_unless(rc == FC_SUCCESS, "should get subset name");
      rc = fc_getSubsetInfo(temp_subsets[j], &temp_numMember, NULL,
			    &temp_assoc);
      fail_unless(rc == FC_SUCCESS, "failed to get subset info");
      rc = fc_getSubsetMembersAsArray(temp_subsets[j], &temp_numMember,
				      &temp_members);
      fail_unless(rc == FC_SUCCESS, "failed to get members");
      switch (temp_assoc){
      case FC_AT_VERTEX:
	sprintf(str_buf, "nodeset_%d", subsetIDs[0]);
	fail_unless(temp_numMember == numMembers[0]-1, "mismatch of numMember");
	for (k = 0; k < temp_numMember; k++) {
	  fail_unless(temp_members[k] == nodes[k+i]-vertOffset-1, 
		      "mismatch of members");
	}
	break;
      case FC_AT_FACE:
	sprintf(str_buf, "sideset_%d", subsetIDs[1]);
	fail_unless(temp_numMember == numMembers[1]-1, "mismatch of numMember");
	for (k = 0; k < temp_numMember; k++) {
	  fail_unless(temp_members[k] == sides[i]-1, "mismatch of members");
	}
	break;
      case FC_AT_ELEMENT:
	sprintf(str_buf, "elemset_%d", subsetIDs[2]);
	fail_unless(temp_numMember == 1, "mismatch of numMember");
	//elems begin numbering again
	fail_unless(temp_members[0] == 0, "mismatch of members"); 
	break;
      default:
	fail_unless(0 == 1, "wrong subset assoc");
      }
      fail_unless(!strcmp(temp_name, str_buf), "subset name mismatch");
      free(temp_name);
      free(temp_members);
    }
    free(temp_subsets);
  }
  free(temp_meshes);

  // all done
  rc = fc_deleteDataset(dataset);
  fail_unless(rc == FC_SUCCESS, "failed to close after testing"); 
}
END_TEST 

START_TEST(exo_whole_var)
{
  FC_ReturnCode rc;
  int i, j, k, m;
  FC_Dataset dataset;
  int numMesh = 3, numPerMesh = 3;
  FC_Mesh mesh, *temp_meshes;
  FC_Variable var, *temp_vars;
  char* ds_name = "dataset";
  int numDataType = 4;
  FC_DataType dataTypes[4] = { FC_DT_CHAR, FC_DT_INT, 
                               FC_DT_FLOAT, FC_DT_DOUBLE };
  char var_names[4][1024], temp_var_name[1024];
  char char_data[3*3+3] = { 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j',
                            'k', 'l' };
  int int_data[3*3+3];
  float flt_data[3*3+3];
  double dbl_data[3*3+3];
  void* datas[4] = { char_data, int_data, flt_data, dbl_data };
  void *data_p;
  // Not all data types work
  int numDataTypeWork = 1;
  // LU index by order seen after read -> order written 
  int dataTypeWorkLU[1] = { 1 }; // not counting the "ID" var
  int temp_numMesh, temp_numVar;
  int temp_numDataPoint, temp_numComp;
  char* temp_name;
  FC_AssociationType temp_assoc;
  FC_MathType temp_mathtype;
  FC_DataType temp_datatype;
  void* temp_data;
  FC_MathType mathTypes[3] = { FC_MT_SCALAR, FC_MT_VECTOR, FC_MT_TENSOR };
  int numComps[3] = { 1, 3, 2 }, totalNumComp = 1+3+2;

  // Exodus always creates a FC_AT_WHOLE_MESH & FC_DT_INT var "ID"

  // setup (don't start w/ zero, too boring for 1 comp whole vars)
  for (i = 0; i < 3*3+3; i++) {
    int_data[i] = (i+1);
    flt_data[i] = (i+1)/10.;
    dbl_data[i] = (i+1)/100.;
  }
 
  // ** Case 1: one of each data type on a single mesh, scalars

  // Create dataset, sequence, hex mesh & seqVars
  rc = fc_createDataset(ds_name, &dataset);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create dataset");
  rc = _fc_makeMeshFromFiles(dataset, vertFiles[7], elemFiles[7], 
			     &mesh);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create mesh");
  for (i = 0; i < numDataType; i++) {
    sprintf(var_names[i], "variable as %s",
	    fc_getDataTypeText(dataTypes[i]));
    rc = fc_createVariable(mesh, var_names[i], &var);
    fail_unless(rc == FC_SUCCESS, "abort: failed to create variable");
    rc = fc_setVariableData(var, 1, 1, FC_AT_WHOLE_MESH,
                            FC_MT_SCALAR, dataTypes[i], datas[i]);
    fail_unless(rc == FC_SUCCESS, "abort: failed to set data on variable");
    //fc_printVariable(var, "hiya", 1);
  }

  // test writing
  rc = fc_writeDataset(dataset, EXODUS_FILE_NAME, FC_FT_EXODUS);
  fail_unless(rc == FC_SUCCESS, "should be able to write dataset");
  
  // test close & reopen
  rc = fc_deleteDataset(dataset);
  fail_unless(rc == FC_SUCCESS, "should close dataset w/o error");    
  rc = fc_loadDatasetWithFormat(EXODUS_FILE_NAME, FC_FT_EXODUS, &dataset);
  fail_unless(rc == FC_SUCCESS, "should load dataset");

  // check variables
  rc = fc_getMeshes(dataset, &temp_numMesh, &temp_meshes);
  fail_unless(rc == FC_SUCCESS, "should find meshes");
  fail_unless(temp_numMesh == 1, "mismatch of numMesh");
  mesh = temp_meshes[0];
  free(temp_meshes);
  rc = fc_getVariables(mesh, &temp_numVar, &temp_vars);
  fail_unless(rc == FC_SUCCESS, "should find vars");
  fail_unless(temp_numVar == numDataTypeWork + 1, "mismatch of numVar");
  // check the "ID" var made by Exodus
  rc = fc_getVariableName(temp_vars[0], &temp_name);
  fail_unless(!strcmp(temp_name, "ID"), "mismatch of name");
  free(temp_name);
  rc = fc_getVariableInfo(temp_vars[0], &temp_numDataPoint,
                          &temp_numComp, &temp_assoc, &temp_mathtype,
                          &temp_datatype);
  fail_unless(rc == FC_SUCCESS, "failed to get variable info");
  fail_unless(temp_numDataPoint == 1, "mismatch of numDataPoint");
  fail_unless(temp_numComp == 1, "mismatch of numComp");
  fail_unless(temp_assoc == FC_AT_WHOLE_MESH, "mismatch of assoc");
  fail_unless(temp_mathtype == FC_MT_SCALAR, "mismatch of mathtype");
  fail_unless(temp_datatype == FC_DT_INT, "mismatch of datatype");
  rc = fc_getVariableDataPtr(temp_vars[0], &temp_data);
  fail_unless(rc == FC_SUCCESS, "failed to get variable  data");
  fail_unless(((int*)temp_data)[0] == 1, "mismatch of ID data");
  for (i = 0; i < numDataTypeWork; i++) {
    FC_Variable temp_var = temp_vars[i+1];
    // check info
    rc = fc_getVariableName(temp_var, &temp_name);
    fail_unless(!strcmp(temp_name, var_names[dataTypeWorkLU[i]]), 
                "mismatch of name");
    free(temp_name);
    rc = fc_getVariableInfo(temp_var, &temp_numDataPoint, &temp_numComp,
                            &temp_assoc, &temp_mathtype, &temp_datatype);
    fail_unless(rc == FC_SUCCESS, "failed to get variable info");
    fail_unless(temp_numDataPoint == 1,
		"mismatch of numDataPoint");
    fail_unless(temp_numComp == 1, "mismatch of numComp");
    fail_unless(temp_assoc == FC_AT_WHOLE_MESH, "mismatch of assoc");
    fail_unless(temp_mathtype == FC_MT_SCALAR, "mismatch of mathtype");
    fail_unless(temp_datatype == dataTypes[dataTypeWorkLU[i]],
                "mismatch of datatype");
    // check data
    rc = fc_getVariableDataPtr(temp_var, &temp_data);
    fail_unless(rc == FC_SUCCESS, "failed to get variable  data");
    fail_unless(!memcmp(temp_data, datas[dataTypeWorkLU[i]],
                        fc_sizeofDataType(dataTypeWorkLU[i])),
                "mismatch of data");
  }
  free(temp_vars);

  // all done
  rc = fc_deleteDataset(dataset);
  fail_unless(rc == FC_SUCCESS, "failed to close after testing"); 
 
  // ** Case 2: multiple of each data type on multiple meshes, scalars

  // Create dataset, sequence, hex meshes & seqVars
  rc = fc_createDataset(ds_name, &dataset);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create dataset");
  for (m = 0; m < numMesh; m++) {
    rc = _fc_makeMeshFromFiles(dataset, vertFiles[7], elemFiles[7], 
			       &mesh);
    fail_unless(rc == FC_SUCCESS, "abort: failed to create mesh");
    for (i = 0; i < numDataType; i++) {
      for (j = 0; j < numPerMesh; j++) {
        if (i == 0)
          data_p = &char_data[m*3+j];
        else if (i == 1)
          data_p = &int_data[m*3+j];
        else if (i == 2)
          data_p = &flt_data[m*3+j];
        else
          data_p = &dbl_data[m*3+j];
	sprintf(temp_var_name, "%s %d", var_names[i], j);
	rc = fc_createVariable(mesh, temp_var_name, &var);
	fail_unless(rc == FC_SUCCESS, "abort: failed to create variable");
        rc = fc_setVariableData(var, 1, 1, FC_AT_WHOLE_MESH,
                                FC_MT_SCALAR, dataTypes[i], data_p);
        fail_unless(rc == FC_SUCCESS, "abort: failed to set data on variable");
        //fc_printVariable(var, "hiya", 1);
      }
    }
  }

  // test writing
  rc = fc_writeDataset(dataset, EXODUS_FILE_NAME, FC_FT_EXODUS);
  fail_unless(rc == FC_SUCCESS, "should be able to write dataset");

  // test close & reopen
  rc = fc_deleteDataset(dataset);
  fail_unless(rc == FC_SUCCESS, "should close dataset w/o error");    
  rc = fc_loadDatasetWithFormat(EXODUS_FILE_NAME, FC_FT_EXODUS, &dataset);
  fail_unless(rc == FC_SUCCESS, "should load dataset");

  // check variables
  rc = fc_getMeshes(dataset, &temp_numMesh, &temp_meshes);
  fail_unless(rc == FC_SUCCESS, "should find meshes");
  fail_unless(temp_numMesh == numMesh, "mismatch of numMesh");
  for (m = 0; m < numMesh; m++) {
    rc = fc_getVariables(temp_meshes[m], &temp_numVar, &temp_vars);
    fail_unless(rc == FC_SUCCESS, "should find seq vars");
    fail_unless(temp_numVar == numDataTypeWork*numPerMesh+1, 
		"mismatch of numVar");
    // check the "ID" var made by Exodus
    rc = fc_getVariableName(temp_vars[0], &temp_name);
    fail_unless(!strcmp(temp_name, "ID"), "mismatch of name");
    free(temp_name);
    rc = fc_getVariableInfo(temp_vars[0], &temp_numDataPoint,
                            &temp_numComp, &temp_assoc, &temp_mathtype,
                            &temp_datatype);
    fail_unless(rc == FC_SUCCESS, "failed to get variable info");
    fail_unless(temp_numDataPoint == 1, "mismatch of numDataPoint");
    fail_unless(temp_numComp == 1, "mismatch of numComp");
    fail_unless(temp_assoc == FC_AT_WHOLE_MESH, "mismatch of assoc");
    fail_unless(temp_mathtype == FC_MT_SCALAR, "mismatch of mathtype");
    fail_unless(temp_datatype == FC_DT_INT, "mismatch of datatype");
    rc = fc_getVariableDataPtr(temp_vars[0], &temp_data);
    fail_unless(rc == FC_SUCCESS, "failed to get variable  data");
    fail_unless(((int*)temp_data)[0] == m+1, "mismatch of ID data");
    for (i = 0; i < numDataTypeWork; i++) {
      for (j = 0; j < numPerMesh; j++) {
        FC_Variable temp_var = temp_vars[i*numPerMesh+j+1];
        if (i == 0)
          data_p = &int_data[m*3+j];
	// check info
	rc = fc_getVariableName(temp_var, &temp_name);
	sprintf(temp_var_name, "%s %d", var_names[dataTypeWorkLU[i]], j);
	fail_unless(!strcmp(temp_name, temp_var_name), "mismatch of name");
	free(temp_name);
	rc = fc_getVariableInfo(temp_var, &temp_numDataPoint,
				&temp_numComp, &temp_assoc, &temp_mathtype,
				&temp_datatype);
	fail_unless(rc == FC_SUCCESS, "failed to get variable info");
	fail_unless(temp_numDataPoint == 1,
		    "mismatch of numDataPoint");
	fail_unless(temp_numComp == 1, "mismatch of numComp");
	fail_unless(temp_assoc == FC_AT_WHOLE_MESH, "mismatch of assoc");
	fail_unless(temp_mathtype == FC_MT_SCALAR, "mismatch of mathtype");
	fail_unless(temp_datatype == dataTypes[dataTypeWorkLU[i]], 
                    "mismatch of datatype");
	// check data
        rc = fc_getVariableDataPtr(temp_var, &temp_data);
        fail_unless(rc == FC_SUCCESS, "failed to get variable step data");
        fail_unless(!memcmp(temp_data, data_p,
                            fc_sizeofDataType(dataTypeWorkLU[i])),
                    "mismatch of data");
      }
    }
    free(temp_vars);
  }
  free(temp_meshes);

  // all done
  rc = fc_deleteDataset(dataset);
  fail_unless(rc == FC_SUCCESS, "failed to close after testing");
 
  // ** Case 3: multiple of each data type on multiple meshes, scalars
  //       BUT don't have complete coverage (leave one off of middle mesh)

  // Create dataset, sequence, hex meshes & seqVars
  rc = fc_createDataset(ds_name, &dataset);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create dataset");
  for (m = 0; m < numMesh; m++) {
    rc = _fc_makeMeshFromFiles(dataset, vertFiles[7], elemFiles[7], 
			       &mesh);
    fail_unless(rc == FC_SUCCESS, "abort: failed to create mesh");
    for (i = 0; i < numDataType; i++) {
      for (j = 0; j < numPerMesh; j++) {
	if (m == 1 && j == 1)
	  continue; // drop one of each kind of variable on middle mesh
        if (i == 0)
          data_p = &char_data[m*3+j];
        else if (i == 1)
          data_p = &int_data[m*3+j];
        else if (i == 2)
          data_p = &flt_data[m*3+j];
        else
          data_p = &dbl_data[m*3+j];
	sprintf(temp_var_name, "%s %d", var_names[i], j);
	rc = fc_createVariable(mesh, temp_var_name, &var);
	fail_unless(rc == FC_SUCCESS, "abort: failed to create variable");
        rc = fc_setVariableData(var, 1, 1, FC_AT_WHOLE_MESH,
                                FC_MT_SCALAR, dataTypes[i], data_p);
        fail_unless(rc == FC_SUCCESS, "abort: failed to set data on variable");
        //fc_printVariable(var, "hiya", 1);
      }
    }
  }
 
  // test writing
  rc = fc_writeDataset(dataset, EXODUS_FILE_NAME, FC_FT_EXODUS);
  fail_unless(rc == FC_SUCCESS, "should be able to write dataset");

  // test close & reopen
  rc = fc_deleteDataset(dataset);
  fail_unless(rc == FC_SUCCESS, "should close dataset w/o error");    
  rc = fc_loadDatasetWithFormat(EXODUS_FILE_NAME, FC_FT_EXODUS, &dataset);
  fail_unless(rc == FC_SUCCESS, "should load dataset");

  // check variables
  rc = fc_getMeshes(dataset, &temp_numMesh, &temp_meshes);
  fail_unless(rc == FC_SUCCESS, "should find meshes");
  fail_unless(temp_numMesh == numMesh, "mismatch of numMesh");
  for (m = 0; m < numMesh; m++) {
    rc = fc_getVariables(temp_meshes[m], &temp_numVar, &temp_vars);
    fail_unless(rc == FC_SUCCESS, "should find seq vars");
    fail_unless(temp_numVar == numDataTypeWork*numPerMesh+1, 
		"mismatch of numVar");
    // check the "ID" var made by Exodus
    rc = fc_getVariableName(temp_vars[0], &temp_name);
    fail_unless(!strcmp(temp_name, "ID"), "mismatch of name");
    free(temp_name);
    rc = fc_getVariableInfo(temp_vars[0], &temp_numDataPoint,
                            &temp_numComp, &temp_assoc, &temp_mathtype,
                            &temp_datatype);
    fail_unless(rc == FC_SUCCESS, "failed to get variable info");
    fail_unless(temp_numDataPoint == 1, "mismatch of numDataPoint");
    fail_unless(temp_numComp == 1, "mismatch of numComp");
    fail_unless(temp_assoc == FC_AT_WHOLE_MESH, "mismatch of assoc");
    fail_unless(temp_mathtype == FC_MT_SCALAR, "mismatch of mathtype");
    fail_unless(temp_datatype == FC_DT_INT, "mismatch of datatype");
    rc = fc_getVariableDataPtr(temp_vars[0], &temp_data);
    fail_unless(rc == FC_SUCCESS, "failed to get variable  data");
    fail_unless(((int*)temp_data)[0] == m+1, "mismatch of ID data");
    for (i = 0; i < numDataTypeWork; i++) {
      for (j = 0; j < numPerMesh; j++) {
        FC_Variable temp_var = temp_vars[i*numPerMesh+j+1];
        if (i == 0)
          data_p = &int_data[m*3+j];
	// check info
	rc = fc_getVariableName(temp_var, &temp_name);
	sprintf(temp_var_name, "%s %d", var_names[dataTypeWorkLU[i]], j);
	fail_unless(!strcmp(temp_name, temp_var_name), "mismatch of name");
	free(temp_name);
	rc = fc_getVariableInfo(temp_var, &temp_numDataPoint,
				&temp_numComp, &temp_assoc, &temp_mathtype,
				&temp_datatype);
	fail_unless(rc == FC_SUCCESS, "failed to get variable info");
	fail_unless(temp_numDataPoint == 1,
		    "mismatch of numDataPoint");
	fail_unless(temp_numComp == 1, "mismatch of numComp");
	fail_unless(temp_assoc == FC_AT_WHOLE_MESH, "mismatch of assoc");
	fail_unless(temp_mathtype == FC_MT_SCALAR, "mismatch of mathtype");
	fail_unless(temp_datatype == dataTypes[dataTypeWorkLU[i]], 
                    "mismatch of datatype");
	// check data
        rc = fc_getVariableDataPtr(temp_var, &temp_data);
        fail_unless(rc == FC_SUCCESS, "failed to get variable step data");
        if (m == 1 && j == 1) // missing ones are set to zero
          fail_unless(((int*)temp_data)[0] == 0, "mismatch of data");
        else 
          fail_unless(!memcmp(temp_data, data_p,
                              fc_sizeofDataType(dataTypeWorkLU[i])),
                      "mismatch of data");
      }
    }
    free(temp_vars);
  }
  free(temp_meshes);

  // all done
  rc = fc_deleteDataset(dataset);
  fail_unless(rc == FC_SUCCESS, "failed to close after testing");

  // ** Case 4: multiple of each data type on multiple meshes, scalars
  //       BUT each has different mathtype

  // Create dataset, sequence, hex meshes & seqVars
  rc = fc_createDataset(ds_name, &dataset);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create dataset");
  for (m = 0; m < numMesh; m++) {
    rc = _fc_makeMeshFromFiles(dataset, vertFiles[7], elemFiles[7], 
			       &mesh);
    fail_unless(rc == FC_SUCCESS, "abort: failed to create mesh");
    for (i = 0; i < numDataType; i++) {
      for (j = 0; j < numPerMesh; j++) {
        if (i == 0)
          data_p = &char_data[m*3+j];
        else if (i == 1)
          data_p = &int_data[m*3+j];
        else if (i == 2)
          data_p = &flt_data[m*3+j];
        else
          data_p = &dbl_data[m*3+j];
	sprintf(temp_var_name, "%s %d", var_names[i], j);
	rc = fc_createVariable(mesh, temp_var_name, &var);
	fail_unless(rc == FC_SUCCESS, "abort: failed to create variable");
        rc = fc_setVariableData(var, 1, numComps[j], FC_AT_WHOLE_MESH,
                                mathTypes[j], dataTypes[i], data_p);
        fail_unless(rc == FC_SUCCESS, "abort: failed to set data on variable");
        //fc_printVariable(var, "hiya", 1);
      }
    }
  }
 
  // test writing
  rc = fc_writeDataset(dataset, EXODUS_FILE_NAME, FC_FT_EXODUS);
  fail_unless(rc == FC_SUCCESS, "should be able to write dataset");

  // test close & reopen
  rc = fc_deleteDataset(dataset);
  fail_unless(rc == FC_SUCCESS, "should close dataset w/o error");    
  rc = fc_loadDatasetWithFormat(EXODUS_FILE_NAME, FC_FT_EXODUS, &dataset);
  fail_unless(rc == FC_SUCCESS, "should load dataset");

  // check variables
  rc = fc_getMeshes(dataset, &temp_numMesh, &temp_meshes);
  fail_unless(rc == FC_SUCCESS, "should find meshes");
  fail_unless(temp_numMesh == numMesh, "mismatch of numMesh");
  for (m = 0; m < numMesh; m++) {
    int idx;
    char vector_endings[3][2] = { "x", "y", "z" };
    rc = fc_getVariables(temp_meshes[m], &temp_numVar, &temp_vars);
    fail_unless(rc == FC_SUCCESS, "should find seq vars");
    fail_unless(temp_numVar == numDataTypeWork*totalNumComp+1, 
		"mismatch of numVar");
    // check the "ID" var made by Exodus
    rc = fc_getVariableName(temp_vars[0], &temp_name);
    fail_unless(!strcmp(temp_name, "ID"), "mismatch of name");
    free(temp_name);
    rc = fc_getVariableInfo(temp_vars[0], &temp_numDataPoint,
                            &temp_numComp, &temp_assoc, &temp_mathtype,
                            &temp_datatype);
    fail_unless(rc == FC_SUCCESS, "failed to get variable info");
    fail_unless(temp_numDataPoint == 1, "mismatch of numDataPoint");
    fail_unless(temp_numComp == 1, "mismatch of numComp");
    fail_unless(temp_assoc == FC_AT_WHOLE_MESH, "mismatch of assoc");
    fail_unless(temp_mathtype == FC_MT_SCALAR, "mismatch of mathtype");
    fail_unless(temp_datatype == FC_DT_INT, "mismatch of datatype");
    rc = fc_getVariableDataPtr(temp_vars[0], &temp_data);
    fail_unless(rc == FC_SUCCESS, "failed to get variable  data");
    fail_unless(((int*)temp_data)[0] == m+1, "mismatch of ID data");
    idx = 1;
    for (i = 0; i < numDataTypeWork; i++) {
      for (j = 0; j < numPerMesh; j++) {
        for (k = 0; k < numComps[j]; k++) {
          FC_Variable temp_var = temp_vars[idx];
          if (i == 0)
            data_p = &int_data[m*3+j];
          // check info
          rc = fc_getVariableName(temp_var, &temp_name);
          if (mathTypes[j] == FC_MT_SCALAR)
            sprintf(temp_var_name, "%s %d", var_names[dataTypeWorkLU[i]], j);
          else if (mathTypes[j] == FC_MT_VECTOR)
            sprintf(temp_var_name, "%s %d_%s", var_names[dataTypeWorkLU[i]], j,
                    vector_endings[k]);
          else
            sprintf(temp_var_name, "%s %d_c%d", var_names[dataTypeWorkLU[i]],
                    j, k);
          fail_unless(!strcmp(temp_name, temp_var_name), "mismatch of name");
          free(temp_name);
          rc = fc_getVariableInfo(temp_var, &temp_numDataPoint,
                                  &temp_numComp, &temp_assoc, &temp_mathtype,
                                  &temp_datatype);
          fail_unless(rc == FC_SUCCESS, "failed to get variable info");
          fail_unless(temp_numDataPoint == 1,
                      "mismatch of numDataPoint");
          fail_unless(temp_numComp == 1, "mismatch of numComp");
          fail_unless(temp_assoc == FC_AT_WHOLE_MESH, "mismatch of assoc");
          fail_unless(temp_mathtype == FC_MT_SCALAR, "mismatch of mathtype");
          fail_unless(temp_datatype == dataTypes[dataTypeWorkLU[i]], 
                      "mismatch of datatype");
          // check data
          rc = fc_getVariableDataPtr(temp_var, &temp_data);
          fail_unless(rc == FC_SUCCESS, "failed to get variable step data");
          fail_unless(!memcmp(temp_data, &((int*)data_p)[k],
                              fc_sizeofDataType(dataTypeWorkLU[i])),
                      "mismatch of data");
          idx++;
        }
      }
    }
    free(temp_vars);
  }
  free(temp_meshes);

  // all done
  rc = fc_deleteDataset(dataset);
  fail_unless(rc == FC_SUCCESS, "failed to close after testing");
}
END_TEST

START_TEST(exo_seq)
{
  FC_ReturnCode rc;
  int i;
  int temp_numSeq, temp_numStep;
  FC_Dataset dataset;
  FC_Mesh mesh;
  FC_Sequence sequence, *temp_sequences, sequence2;
  FC_DataType temp_dataType;
  char* ds_name = "dataset";
  char* seq_name = "sequence";
  char* temp_name;
  int numStep = 6;
  char char_times[6];
  int int_times[6];
  float flt_times[6];
  double dbl_times[6];
  double dbl_times2[6];
  void* times[4], *temp_times;
  int numDataType = 4;
  FC_DataType dataTypes[4] = { FC_DT_CHAR, FC_DT_INT, FC_DT_FLOAT,
			       FC_DT_DOUBLE };
  
  // setup
  for (i = 0; i < numStep; i++) {
    char_times[i] = (char)i;
    int_times[i] = i;
    flt_times[i] = i;
    dbl_times[i] = i;
    dbl_times2[i] = i/2.;
  }
  times[0] = char_times;
  times[1] = int_times;
  times[2] = flt_times;
  times[3] = dbl_times;
  
  // ** Case 1 - a single sequence of each type (with a mesh, but no vars)

  for (i = 0; i < numDataType; i++) {
    
    // Create dataset, mesh & sequence
    rc = fc_createDataset(ds_name, &dataset);
    fail_unless(rc == FC_SUCCESS, "abort: failed to create dataset");
    rc = _fc_makeMeshFromFiles(dataset, vertFiles[3], elemFiles[3], 
			       &mesh);
    fail_unless(rc == FC_SUCCESS, "abort: failed to create mesh");
    rc = fc_createSequence(dataset, seq_name, &sequence);
    fail_unless(rc == FC_SUCCESS, "abort: failed to create seqeunce");
    rc = fc_setSequenceCoords(sequence, numStep, dataTypes[i], times[i]);
    fail_unless(rc == FC_SUCCESS, "abort: failed to set seq coords");
    
    // test writing
    rc = fc_writeDataset(dataset, EXODUS_FILE_NAME, FC_FT_EXODUS);
    fail_unless(rc == FC_SUCCESS, "should be able to write dataset");
    
    // test close & reopen
    rc = fc_deleteDataset(dataset);
    fail_unless(rc == FC_SUCCESS, "should close dataset w/o error");    
    rc = fc_loadDatasetWithFormat(EXODUS_FILE_NAME, FC_FT_EXODUS, &dataset);
    fail_unless(rc == FC_SUCCESS, "should load dataset");
    
    // check sequence
    rc = fc_getSequences(dataset, &temp_numSeq, &temp_sequences);
    fail_unless(rc == FC_SUCCESS, "failed to get sequences");
    fail_unless(temp_numSeq == 1, "mismatch of numSequence");
    rc = fc_getSequenceName(temp_sequences[0], &temp_name);
    fail_unless(rc == FC_SUCCESS, "should get sequence name");
    fail_unless(!strcmp(temp_name, "time"), "mismatch of name");
    free(temp_name);
    rc = fc_getSequenceInfo(temp_sequences[0], &temp_numStep, &temp_dataType);
    fail_unless(rc == FC_SUCCESS, "failed to get info");
    fail_unless(temp_numStep == numStep, "mismatch of numStep");
    fail_unless(temp_dataType == FC_DT_DOUBLE, "mismatch of datatype");
    rc = fc_getSequenceCoordsPtr(temp_sequences[0], &temp_times);
    fail_unless(rc == FC_SUCCESS, "failed to get sequence coords");
    fail_unless(!memcmp(temp_times, dbl_times, numStep*sizeof(double)),
		"mismatch of coords");
    free(temp_sequences);

    // all done
    rc = fc_deleteDataset(dataset);
    fail_unless(rc == FC_SUCCESS, "failed to close after testing"); 
  }

  // ** Case 2 - two sequence (with a mesh but no vars)
  //    (2nd sequence will get clobbered)

  // Create dataset, mesh & sequences
  rc = fc_createDataset(ds_name, &dataset);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create dataset");
  rc = _fc_makeMeshFromFiles(dataset, vertFiles[3], elemFiles[3], 
			     &mesh);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create mesh");
  rc = fc_createSequence(dataset, "another sequence", &sequence2);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create seqeunce");
  rc = fc_setSequenceCoords(sequence2, numStep, FC_DT_DOUBLE, dbl_times2);
  fail_unless(rc == FC_SUCCESS, "abort: failed to set seq coords");
  rc = fc_createSequence(dataset, seq_name, &sequence);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create seqeunce");
  rc = fc_setSequenceCoords(sequence, numStep, FC_DT_DOUBLE, dbl_times);
  fail_unless(rc == FC_SUCCESS, "abort: failed to set seq coords");
  
  // test writing
  rc = fc_writeDataset(dataset, EXODUS_FILE_NAME, FC_FT_EXODUS);
  fail_unless(rc == FC_SUCCESS, "should be able to write dataset");
  
  // test close & reopen
  rc = fc_deleteDataset(dataset);
  fail_unless(rc == FC_SUCCESS, "should close dataset w/o error");    
  rc = fc_loadDatasetWithFormat(EXODUS_FILE_NAME, FC_FT_EXODUS, &dataset);
  fail_unless(rc == FC_SUCCESS, "should load dataset");
  
  // check sequence - only 1!
  rc = fc_getSequences(dataset, &temp_numSeq, &temp_sequences);
  fail_unless(rc == FC_SUCCESS, "failed to get sequences");
  fail_unless(temp_numSeq == 1, "mismatch of numSequence");
  rc = fc_getSequenceName(temp_sequences[0], &temp_name);
  fail_unless(rc == FC_SUCCESS, "should get sequence name");
  fail_unless(!strcmp(temp_name, "time"), "mismatch of name");
  free(temp_name);
  rc = fc_getSequenceInfo(temp_sequences[0], &temp_numStep, &temp_dataType);
  fail_unless(rc == FC_SUCCESS, "failed to get info");
  fail_unless(temp_numStep == numStep, "mismatch of numStep");
  fail_unless(temp_dataType == FC_DT_DOUBLE, "mismatch of datatype");
  rc = fc_getSequenceCoordsPtr(temp_sequences[0], &temp_times);
  fail_unless(rc == FC_SUCCESS, "failed to get sequence coords");
  fail_unless(!memcmp(temp_times, dbl_times2, numStep*sizeof(double)),
	      "mismatch of coords");
  free(temp_sequences);
  
  // all done
  rc = fc_deleteDataset(dataset);
  fail_unless(rc == FC_SUCCESS, "failed to close after testing"); 
}
END_TEST

// greater than the first one, same as the 2nd one, less than the 3rd
START_TEST(exo_var_cmp)
{
  FC_ReturnCode rc;
  int i, j, k, m, n, p, i2, j2, k2, m2, n2, p2;
  int num = 3;
  FC_Dataset dataset;
  FC_Mesh meshes[3];
  FC_Variable variables[3][3][3][3][3][3]; 
  _FC_MeshSlot* meshSlot;
  _FC_VarSlot* varSlot;
  char* ds_name = "dataset";
  char* var_names[3] = { "apple", "truth", "zebra" } ;
  int numComps[3] = { 1, 2, 3 };
  FC_AssociationType assocs[3] = { FC_AT_VERTEX, FC_AT_EDGE, FC_AT_ELEMENT };
  FC_MathType mathtypes[3] = { FC_MT_SCALAR, FC_MT_VECTOR, FC_MT_TENSOR };
  FC_DataType datatypes[3] = { FC_DT_INT, FC_DT_FLOAT, FC_DT_DOUBLE };
  FC_AssociationType temp_assoc;
  FC_MathType temp_mathtype;
  FC_DataType temp_datatype;
  
  // setup - make sure assocs, mathtypes & datatypes are sorted
  for (i = 2; i > 0; i--) {
    if (assocs[i] < assocs[i-1]) {
      temp_assoc = assocs[i-1];
      assocs[i-1] = assocs[i];
      assocs[i] = temp_assoc;
    }
  }
  for (i = 2; i > 0; i--) {
    if (mathtypes[i] < mathtypes[i-1]) {
      temp_mathtype = mathtypes[i-1];
      mathtypes[i-1] = mathtypes[i];
      mathtypes[i] = temp_mathtype;
    }
  }
  for (i = 2; i > 0; i--) {
    if (datatypes[i] < datatypes[i-1]) {
      temp_datatype = datatypes[i-1];
      datatypes[i-1] = datatypes[i];
      datatypes[i] = temp_datatype;
    }
  }
  fail_unless(assocs[2] > assocs[1] && assocs[1] > assocs[0], 
	      "abort: assocs not sorted");
  fail_unless(mathtypes[2] > mathtypes[1] && mathtypes[1] > mathtypes[0], 
	      "abort: mathtypes not sorted");
  fail_unless(datatypes[2] > datatypes[1] && datatypes[1] > datatypes[0], 
	      "abort: datatypes not sorted");

  // setup - Make some fake vars
  rc = fc_createDataset(ds_name, &dataset);
  fail_unless(rc == FC_SUCCESS, "failed to create dataset");
  for (p = 0; p < num; p++) { // meshID
    rc = fc_createMesh(dataset, "mesh", &meshes[p]);
    fail_unless(rc == FC_SUCCESS, "failed to create mesh");
    meshSlot = _fc_getMeshSlot(meshes[p]);
    meshSlot->fileInfo.exoInfo.meshID = p;
    for (i = 0; i < num; i++) { // assoc type
      for (j = 0; j < num; j++) { // name
	for (k = 0; k < num; k++) { // numComp
	  for (m = 0; m < num; m++) { // mathtype
	    for (n = 0; n < num; n++) { // datatype
	      rc = fc_createVariable(meshes[p], var_names[j],
				     &variables[i][j][k][m][n][p]); 
	      fail_unless(rc == FC_SUCCESS, "failed to create variable");
	      varSlot = _fc_getVarSlot(variables[i][j][k][m][n][p]);
	      varSlot->assoc = assocs[i];
	      varSlot->numComponent = numComps[k];
	      varSlot->mathtype = mathtypes[m];
	      varSlot->datatype = datatypes[n];
	    }
	  }
	}
      }
    }
  }

  // test all possible combinations
  for (i = 0; i < num; i++) {
    for (i2 = 0; i2 < num; i2++) {
      for (j = 0; j < num; j++) {
	for (j2 = 0; j2 < num; j2++) {
	  for (k = 0; k < num; k++) {
	    for (k2 = 0; k2 < num; k2++) {
	      for (m = 0; m < num; m++) {
		for (m2 = 0; m2 < num; m2++) {
		  for (n = 0; n < num; n++) {
		    for (n2 = 0; n2 < num; n2++) {
		      for (p = 0; p < num; p++) {
			for (p2 = 0; p2 < num; p2++) {
			  //printf("i = %d, j = %d, k = %d, m = %d, n = %d, p = %d\n",
			  //i, j, k, m, n, p);
			  FC_Variable *a, *b;
			  int cmp1, cmp2;
			  a = &variables[i][j][k][m][n][p];
			  b = &variables[i2][j2][k2][m2][n2][p2];
			  cmp1 = _fc_exo_seqVar_cmp(a, b);
			  cmp2 = _fc_exo_seqVar_meshID_cmp(a, b);
			  if (i < i2) 
			    fail_unless(cmp1 == -1 && cmp2 == -1, 
					"problem with assoc cmp");
			  else if (i > i2)
			    fail_unless(cmp1 == 1 && cmp2 == 1, 
					"problem with assoc cmp");
			  else {
			    if (j < j2)
			      fail_unless(cmp1 == -1 && cmp2 == -1, 
					  "problem with name cmp");
			    else if (j > j2)
			      fail_unless(cmp1 == 1 && cmp2 == 1, 
					  "problem with name cmp");
			    else {
			      if (k < k2)
				fail_unless(cmp1 == -1 && cmp2 == -1, 
					    "problem with numComp cmp");
			      else if (k > k2)
				fail_unless(cmp1 == 1 && cmp2 == 1, 
					    "problem with numComp cmp");
			      else {
				if (m < m2)
				  fail_unless(cmp1 == -1 && cmp2 == -1, 
					      "problem with mathtype cmp");
				else if (m > m2)
				  fail_unless(cmp1 == 1 && cmp2 == 1, 
					      "problem with mathtype cmp");
				else {
				  if (n < n2)
				    fail_unless(cmp1 == -1 && cmp2 == -1, 
						"problem with datatype cmp");
				  else if (n > n2)
				    fail_unless(cmp1 == 1 && cmp2 == 1, 
						"problem with datatype cmp");
				  else {
				    fail_unless(cmp1 == 0, 
						"equal except for meshID problem");
				    if (p < p2)
				      fail_unless(cmp2 == -1, "problem with meshID");
				    else if (p > p2)
				      fail_unless(cmp2 == 1, "problem with meshID");
				    else
				      fail_unless(cmp2 == 0, "equal should return 0");
				  }
				}
			      }
			    }
			  }
			}
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
  
  // all done
  rc = fc_deleteDataset(dataset);
  fail_unless(rc == FC_SUCCESS, "failed to close after testing"); 
}
END_TEST

START_TEST(exo_global_seq_var)
{
  FC_ReturnCode rc;
  int i, j, k;
  FC_Dataset dataset;
  FC_Mesh mesh;
  FC_Sequence sequence;
  int numDataType = 4;
  FC_DataType dataTypes[4] = { FC_DT_CHAR, FC_DT_INT,
                               FC_DT_FLOAT, FC_DT_DOUBLE };
  int numMathType = 3;
  FC_MathType mathTypes[3] = { FC_MT_SCALAR, FC_MT_VECTOR,
                               FC_MT_SYMTENSOR };
  int numComps[3] = { 1, 3, 6 }, maxNumComp = 6, numTotalComp = 1+3+6;
  FC_Variable *glbSeqVar, **temp_glbSeqVars;
  char* ds_name = "dataset";
  char* seq_name = "sequence";
  char var_names[numDataType][1024], temp_var_name[1024];
  int numStep = 5;
  double times[5];
  char char_data[6*5];
  int int_data[6*5];
  float flt_data[6*5];
  double dbl_data[6*5];
  void *data_p;
  int temp_numGlbSeqVar, *temp_numSteps;
  int temp_numDataPoint, temp_numComp;
  char* temp_name;
  FC_AssociationType temp_assoc;
  FC_MathType temp_mathtype;
  FC_DataType temp_datatype;
  void* temp_data;
  int idx;
  char vector_endings[3][2] = { "x", "y", "z" };

  // setup
  for (i = 0; i < numStep; i++)
    times[i] = i/10.;
  // get unique data starting at different point inside of data
  // for example, timestep 0: &data[0], timestep 1: &data[1] ...
  for (i = 0; i < maxNumComp*numStep; i++) {
    char_data[i] = 'a' + i;
    int_data[i] = 1 + i;
    flt_data[i] = 1 + i/10.;
    dbl_data[i] = 1 + i/100.;
  }

  // ** Case 1: One of each data type, scalars

  // Create dataset, mesh & sequence
  rc = fc_createDataset(ds_name, &dataset);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create dataset");
  rc = _fc_makeMeshFromFiles(dataset, vertFiles[3], elemFiles[3], 
                             &mesh);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create mesh");
  rc = fc_createSequence(dataset, seq_name, &sequence);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create sequence");
  rc = fc_setSequenceCoords(sequence, numStep, FC_DT_DOUBLE, times);
  fail_unless(rc == FC_SUCCESS, "abort: failed to set sequence coords");
  // create glbvars
  for (i = 0; i < numDataType; i++) {
    sprintf(var_names[i], "glbVar %s", fc_getDataTypeText(dataTypes[i]));
    rc = fc_createGlobalSeqVariable(dataset, sequence, var_names[i], &numStep, 
                                    &glbSeqVar);
    fail_unless(rc == FC_SUCCESS, "abort: failed to create glb seq var");
    for (j = 0; j < numStep; j++) {
      if (i == 0)
        data_p = &char_data[j];
      else if (i == 1)
        data_p = &int_data[j];
      else if (i == 2)
        data_p = &flt_data[j];
      else 
        data_p = &dbl_data[j];
      rc = fc_setVariableData(glbSeqVar[j], 1, 1, FC_AT_WHOLE_DATASET,
                              FC_MT_SCALAR, dataTypes[i], data_p);
      fail_unless(rc == FC_SUCCESS, "abort: failed to set data on variable");
    }
    //fc_printSeqVariable(numStep, glbSeqVar, "", 1);
    free(glbSeqVar);
  }

  // test writing
  rc = fc_writeDataset(dataset, EXODUS_FILE_NAME, FC_FT_EXODUS);
  fail_unless(rc == FC_SUCCESS, "should be able to write dataset");

  // test close & reopen
  rc = fc_deleteDataset(dataset);
  fail_unless(rc == FC_SUCCESS, "should close dataset w/o error");    
  rc = fc_loadDatasetWithFormat(EXODUS_FILE_NAME, FC_FT_EXODUS, &dataset);
  fail_unless(rc == FC_SUCCESS, "should load dataset");
 
  // check variables
  rc = fc_getGlobalSeqVariables(dataset, &temp_numGlbSeqVar, &temp_numSteps,
                                &temp_glbSeqVars);
  fail_unless(rc == FC_SUCCESS, "should find glb seq vars");
  fail_unless(temp_numGlbSeqVar == numDataType, "mismatch of numGlbSeqVar");
  for (i = 0; i < numDataType; i++) {
    // check info
    fail_unless(temp_numSteps[i] == numStep, "mismatch of numStep");
    rc = fc_getVariableName(temp_glbSeqVars[i][0], &temp_name);
    fail_unless(!strcmp(temp_name, var_names[i]), "mismatch of name");
    free(temp_name);
    rc = fc_getVariableInfo(temp_glbSeqVars[i][0], &temp_numDataPoint,
                            &temp_numComp, &temp_assoc, &temp_mathtype,
                            &temp_datatype);
    fail_unless(rc == FC_SUCCESS, "failed to get variable info");
    fail_unless(temp_numDataPoint == 1, "mismatch of numDataPoint");
    fail_unless(temp_numComp == 1, "mismatch of numComp");
    fail_unless(temp_assoc == FC_AT_WHOLE_DATASET, "mismatch of assoc");
    fail_unless(temp_mathtype == FC_MT_SCALAR, "mismatch of mathtype");
    fail_unless(temp_datatype == FC_DT_DOUBLE, "mismatch of datatype");
    // check data
    for (j = 0; j < numStep; j++) {
      rc = fc_getVariableDataPtr(temp_glbSeqVars[i][j], &temp_data);
      fail_unless(rc == FC_SUCCESS, "failed to get variable step data");
      if (i == 0) 
        fail_unless(((double*)temp_data)[0] == char_data[j],
                      "mismatch of char data");
      else if (i == 1) 
        fail_unless(((double*)temp_data)[0] == int_data[j],
                    "mismatch of int data");
      else if (i == 2)
        fail_unless(((double*)temp_data)[0] == flt_data[j],
			  "mismatch of float data");
      else
        fail_unless(((double*)temp_data)[0] == dbl_data[j],
                    "mismatch of double data");
    }
    free(temp_glbSeqVars[i]);
  }
  free(temp_numSteps);
  free(temp_glbSeqVars);
 
  // all done
  rc = fc_deleteDataset(dataset);
  fail_unless(rc == FC_SUCCESS, "failed to close after testing"); 

  // ** Case 2: different math types, all doubles

  // Create dataset, mesh & sequence
  rc = fc_createDataset(ds_name, &dataset);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create dataset");
  rc = _fc_makeMeshFromFiles(dataset, vertFiles[3], elemFiles[3], 
                             &mesh);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create mesh");
  rc = fc_createSequence(dataset, seq_name, &sequence);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create sequence");
  rc = fc_setSequenceCoords(sequence, numStep, FC_DT_DOUBLE, times);
  fail_unless(rc == FC_SUCCESS, "abort: failed to set sequence coords");
  // create glbvars
  for (i = 0; i < numMathType; i++) {
    sprintf(var_names[i], "glbVar %s", fc_getMathTypeText(mathTypes[i]));
    rc = fc_createGlobalSeqVariable(dataset, sequence, var_names[i], &numStep, 
                                    &glbSeqVar);
    fail_unless(rc == FC_SUCCESS, "abort: failed to create glb seq var");
    for (j = 0; j < numStep; j++) {
      data_p = &dbl_data[j*numComps[i]];
      rc = fc_setVariableData(glbSeqVar[j], 1, numComps[i], FC_AT_WHOLE_DATASET,
                              mathTypes[i], FC_DT_DOUBLE, data_p);
      fail_unless(rc == FC_SUCCESS, "abort: failed to set data on variable");
    }
    //fc_printSeqVariable(numStep, glbSeqVar, "", 1);
    free(glbSeqVar);
  }
 
  // test writing
  rc = fc_writeDataset(dataset, EXODUS_FILE_NAME, FC_FT_EXODUS);
  fail_unless(rc == FC_SUCCESS, "should be able to write dataset");

  // test close & reopen
  rc = fc_deleteDataset(dataset);
  fail_unless(rc == FC_SUCCESS, "should close dataset w/o error");    
  rc = fc_loadDatasetWithFormat(EXODUS_FILE_NAME, FC_FT_EXODUS, &dataset);
  fail_unless(rc == FC_SUCCESS, "should load dataset");
  
  // check variables
  rc = fc_getGlobalSeqVariables(dataset, &temp_numGlbSeqVar, &temp_numSteps,
                                &temp_glbSeqVars);
  fail_unless(rc == FC_SUCCESS, "should find glb seq vars");
  fail_unless(temp_numGlbSeqVar == numTotalComp, "mismatch of numGlbSeqVar");
  idx = 0;
  for (i = 0; i < numMathType; i++) {
    for (j = 0; j < numComps[i]; j++) {
      // check info
      fail_unless(temp_numSteps[idx] == numStep, "mismatch of numStep");
      rc = fc_getVariableName(temp_glbSeqVars[idx][0], &temp_name);
      if (mathTypes[i] == FC_MT_SCALAR)
        sprintf(temp_var_name, "%s", var_names[i]);
      else if (mathTypes[i] == FC_MT_VECTOR) 
        sprintf(temp_var_name, "%s_%s", var_names[i], vector_endings[j]);
      else 
        sprintf(temp_var_name, "%s_c%d", var_names[i], j);
      fail_unless(!strcmp(temp_name, temp_var_name), "mismatch of name");
      free(temp_name);
      rc = fc_getVariableInfo(temp_glbSeqVars[idx][0], &temp_numDataPoint,
                              &temp_numComp, &temp_assoc, &temp_mathtype,
                              &temp_datatype);
      fail_unless(rc == FC_SUCCESS, "failed to get variable info");
      fail_unless(temp_numDataPoint == 1, "mismatch of numDataPoint");
      fail_unless(temp_numComp == 1, "mismatch of numComp");
      fail_unless(temp_assoc == FC_AT_WHOLE_DATASET, "mismatch of assoc");
      fail_unless(temp_mathtype == FC_MT_SCALAR, "mismatch of mathtype");
      fail_unless(temp_datatype == FC_DT_DOUBLE, "mismatch of datatype");
      // check data
      for (k = 0; k < numStep; k++) {
        rc = fc_getVariableDataPtr(temp_glbSeqVars[idx][k], &temp_data);
        fail_unless(rc == FC_SUCCESS, "failed to get variable step data");
        fail_unless(((double*)temp_data)[0] == dbl_data[j+numComps[i]*k],
                    "mismatch of double data");
      }
      free(temp_glbSeqVars[idx]);
      idx++;
    }
  }
  free(temp_numSteps);
  free(temp_glbSeqVars);

  // all done
  rc = fc_deleteDataset(dataset);
  fail_unless(rc == FC_SUCCESS, "failed to close after testing");
}
END_TEST

// FIX? might be able to speed up a little by moving last 3 cases to a
// separate test and don't bother to made edge & face variables.
START_TEST(exo_seq_var)
{
  FC_ReturnCode rc;
  int i, j, k, m, n, p;
  FC_Dataset dataset;
  FC_Sequence sequence;
  int numMesh = 3, numPerMesh = 3;
  FC_Mesh mesh, *temp_meshes;
  FC_Variable *seqVar, **temp_seqVars;
  char* ds_name = "dataset";
  char* seq_name = "sequence";
  char var_names[numAssocType][1024], temp_var_name[1024];
  int numStep = 6;
  double times[6];
  int numEntities[numAssocType], maxNumEntity = 2945; // for a hex
  double data[4*2945];
  char char_data[4*2945];
  int int_data[4*2945];
  float float_data[4*2945];
  void *data_p;
  // Not all assoc types work
  int numAssocTypeWork = 2;
  // LU index by order seen after read to order written 
  int assocTypeWorkLU[2] = { 0, 3 };
  int temp_numMesh, temp_numSeqVar, *temp_numSteps;
  int temp_numDataPoint, temp_numComp;
  char* temp_name;
  FC_AssociationType temp_assoc;
  FC_MathType temp_mathtype;
  FC_DataType temp_datatype;
  void* temp_data;
  FC_MathType mathTypes[3] = { FC_MT_SCALAR, FC_MT_VECTOR, FC_MT_TENSOR };
  int numComps[3] = { 1, 3, 2 }, totalNumComp = 1+3+2;
  FC_DataType dataTypes[3] = { FC_DT_CHAR, FC_DT_INT, FC_DT_FLOAT };

  // setup
  for (i = 0; i < numStep; i++)
    times[i] = i/10.;
  // get unique data starting at different point inside of data
  // for example, timestep 0: &data[0], timestep 1: &data[1] ...
  for (i = 0; i < 4*maxNumEntity; i++) {
    data[i] = 1 + i/100.;
    char_data[i] = (char)data[i];
    int_data[i] = (int)data[i];
    float_data[i] = (float)data[i];
  }

  // ** Case 1: one of each association type on a single mesh, scalars

  // Create dataset, sequence, hex mesh & seqVars
  rc = fc_createDataset(ds_name, &dataset);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create dataset");
  rc = fc_createSequence(dataset, seq_name, &sequence);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create seqeunce");
  rc = fc_setSequenceCoords(sequence, numStep, FC_DT_DOUBLE, times);
  fail_unless(rc == FC_SUCCESS, "abort: failed to set seq coords");
  rc = _fc_makeMeshFromFiles(dataset, vertFiles[7], elemFiles[7], 
			     &mesh);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create mesh");
  // save numEnities for repeated use
  for (i = 0; i < numAssocType; i++) {
    fc_getMeshNumEntity(mesh, assocTypes[i], &numEntities[i]);
    if (numEntities[i] > maxNumEntity)
      maxNumEntity = numEntities[i];
  }
  for (i = 0; i < numAssocType; i++) {
    sprintf(var_names[i], "variable on %s",
	    fc_getAssociationTypeText(assocTypes[i]));
    rc = fc_createSeqVariable(mesh, sequence, var_names[i], &numStep, &seqVar);
    fail_unless(rc == FC_SUCCESS, "abort: failed to create variable");
    for (j = 0; j < numStep; j++) {
      rc = fc_setVariableData(seqVar[j], numEntities[i], 1, assocTypes[i],
			      FC_MT_SCALAR, FC_DT_DOUBLE, 
			      &data[i*numStep+j]);
      fail_unless(rc == FC_SUCCESS, "abort: failed to set data on variable");
    }
    free(seqVar);
  }

  // test writing
  rc = fc_writeDataset(dataset, EXODUS_FILE_NAME, FC_FT_EXODUS);
  fail_unless(rc == FC_SUCCESS, "should be able to write dataset");
  
  // test close & reopen
  rc = fc_deleteDataset(dataset);
  fail_unless(rc == FC_SUCCESS, "should close dataset w/o error");    
  rc = fc_loadDatasetWithFormat(EXODUS_FILE_NAME, FC_FT_EXODUS, &dataset);
  fail_unless(rc == FC_SUCCESS, "should load dataset");

  // check variables
  rc = fc_getMeshes(dataset, &temp_numMesh, &temp_meshes);
  fail_unless(rc == FC_SUCCESS, "should find meshes");
  fail_unless(temp_numMesh == 1, "mismatch of numMesh");
  mesh = temp_meshes[0];
  free(temp_meshes);
  rc = fc_getSeqVariables(mesh, &temp_numSeqVar, &temp_numSteps,
			  &temp_seqVars);
  fail_unless(rc == FC_SUCCESS, "should find seq vars");
  fail_unless(temp_numSeqVar == numAssocTypeWork, "mismatch of numSeqVar");
  for (i = 0; i < temp_numSeqVar; i++) {
    // check info
    fail_unless(temp_numSteps[i] == numStep, "mismatch of numStep");
    rc = fc_getVariableName(temp_seqVars[i][0], &temp_name);
    fail_unless(!strcmp(temp_name, var_names[assocTypeWorkLU[i]]), 
		"mismatch of name");
    free(temp_name);
    rc = fc_getVariableInfo(temp_seqVars[i][0], &temp_numDataPoint,
			    &temp_numComp, &temp_assoc, &temp_mathtype,
			    &temp_datatype);
    fail_unless(rc == FC_SUCCESS, "failed to get variable info");
    fail_unless(temp_numDataPoint == numEntities[assocTypeWorkLU[i]],
		"mismatch of numDataPoint");
    fail_unless(temp_numComp == 1, "mismatch of numComp");
    fail_unless(temp_assoc == assocTypes[assocTypeWorkLU[i]],
		"mismatch of assoc");
    fail_unless(temp_mathtype == FC_MT_SCALAR, "mismatch of mathtype");
    fail_unless(temp_datatype == FC_DT_DOUBLE, "mismatch of datatype");
    // check data
    for (j = 0; j < numStep; j++) {
      rc = fc_getVariableDataPtr(temp_seqVars[i][j], &temp_data);
      fail_unless(rc == FC_SUCCESS, "failed to get variable step data");
      fail_unless(!memcmp(temp_data, &data[assocTypeWorkLU[i]*numStep+j],
			  temp_numDataPoint*sizeof(double)),
		  "mismatch of data");
    }
    free(temp_seqVars[i]);
  }
  free(temp_numSteps);
  free(temp_seqVars);

  // all done
  rc = fc_deleteDataset(dataset);
  fail_unless(rc == FC_SUCCESS, "failed to close after testing"); 

  // ** Case 2: multiple of each association type on multiple meshes, scalars

  // Create dataset, sequence, hex meshes & seqVars
  rc = fc_createDataset(ds_name, &dataset);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create dataset");
  rc = fc_createSequence(dataset, seq_name, &sequence);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create seqeunce");
  rc = fc_setSequenceCoords(sequence, numStep, FC_DT_DOUBLE, times);
  fail_unless(rc == FC_SUCCESS, "abort: failed to set seq coords");
  for (m = 0; m < numMesh; m++) {
    rc = _fc_makeMeshFromFiles(dataset, vertFiles[7], elemFiles[7], 
			       &mesh);
    fail_unless(rc == FC_SUCCESS, "abort: failed to create mesh");
    for (i = 0; i < numAssocType; i++) {
      for (j = 0; j < numPerMesh; j++) {
	sprintf(temp_var_name, "%s %d", var_names[i], j);
	rc = fc_createSeqVariable(mesh, sequence, temp_var_name, &numStep, 
				  &seqVar);
	fail_unless(rc == FC_SUCCESS, "abort: failed to create variable");
	for (k = 0; k < numStep; k++) {
	  int uid = ((m*numAssocType + i)*numPerMesh + j)*numStep + k;
	  rc = fc_setVariableData(seqVar[k], numEntities[i], 1, assocTypes[i],
				  FC_MT_SCALAR, FC_DT_DOUBLE, &data[uid]);
	  fail_unless(rc == FC_SUCCESS, "abort: failed to set data on variable");
	}
	free(seqVar);
      }
    }
  }
  
  // test writing
  rc = fc_writeDataset(dataset, EXODUS_FILE_NAME, FC_FT_EXODUS);
  fail_unless(rc == FC_SUCCESS, "should be able to write dataset");

  // test close & reopen
  rc = fc_deleteDataset(dataset);
  fail_unless(rc == FC_SUCCESS, "should close dataset w/o error");    
  rc = fc_loadDatasetWithFormat(EXODUS_FILE_NAME, FC_FT_EXODUS, &dataset);
  fail_unless(rc == FC_SUCCESS, "should load dataset");

  // check variables
  rc = fc_getMeshes(dataset, &temp_numMesh, &temp_meshes);
  fail_unless(rc == FC_SUCCESS, "should find meshes");
  fail_unless(temp_numMesh == numMesh, "mismatch of numMesh");
  for (m = 0; m < numMesh; m++) {
    rc = fc_getSeqVariables(temp_meshes[m], &temp_numSeqVar, &temp_numSteps,
			    &temp_seqVars);
    fail_unless(rc == FC_SUCCESS, "should find seq vars");
    fail_unless(temp_numSeqVar == numAssocTypeWork*numPerMesh, 
		"mismatch of numSeqVar");
    for (i = 0; i < numAssocTypeWork; i++) {
      for (j = 0; j < numPerMesh; j++) {
	int idx = i*numPerMesh+j;
	// check info
	fail_unless(temp_numSteps[idx] == numStep, "mismatch of numStep");
	rc = fc_getVariableName(temp_seqVars[idx][0], &temp_name);
	sprintf(temp_var_name, "%s %d", var_names[assocTypeWorkLU[i]], j);
	fail_unless(!strcmp(temp_name, temp_var_name), "mismatch of name");
	free(temp_name);
	rc = fc_getVariableInfo(temp_seqVars[idx][0], &temp_numDataPoint,
				&temp_numComp, &temp_assoc, &temp_mathtype,
				&temp_datatype);
	fail_unless(rc == FC_SUCCESS, "failed to get variable info");
	fail_unless(temp_numDataPoint == numEntities[assocTypeWorkLU[i]],
		    "mismatch of numDataPoint");
	fail_unless(temp_numComp == 1, "mismatch of numComp");
	fail_unless(temp_assoc == assocTypes[assocTypeWorkLU[i]],
		    "mismatch of assoc");
	fail_unless(temp_mathtype == FC_MT_SCALAR, "mismatch of mathtype");
	fail_unless(temp_datatype == FC_DT_DOUBLE, "mismatch of datatype");
	// check data
	for (k = 0; k < numStep; k++) {
	  int uid = ((m*numAssocType + assocTypeWorkLU[i])*numPerMesh + 
                     j)*numStep + k;
	  rc = fc_getVariableDataPtr(temp_seqVars[idx][k], &temp_data);
	  fail_unless(rc == FC_SUCCESS, "failed to get variable step data");
	  fail_unless(!memcmp(temp_data, &data[uid],
			      temp_numDataPoint*sizeof(double)),
		      "mismatch of data");
	}
	free(temp_seqVars[idx]);
      }
    }
    free(temp_numSteps);
    free(temp_seqVars);
  }
  free(temp_meshes);

  // all done
  rc = fc_deleteDataset(dataset);
  fail_unless(rc == FC_SUCCESS, "failed to close after testing"); 

  // ** Case 3: multiple of each association type on multiple meshes, scalars
  //       BUT don't have complete coverage (leave one off of middle mesh)

  // Create dataset, sequence, hex meshes & seqVars
  rc = fc_createDataset(ds_name, &dataset);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create dataset");
  rc = fc_createSequence(dataset, seq_name, &sequence);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create seqeunce");
  rc = fc_setSequenceCoords(sequence, numStep, FC_DT_DOUBLE, times);
  fail_unless(rc == FC_SUCCESS, "abort: failed to set seq coords");
  for (m = 0; m < numMesh; m++) {
    rc = _fc_makeMeshFromFiles(dataset, vertFiles[7], elemFiles[7], 
			       &mesh);
    fail_unless(rc == FC_SUCCESS, "abort: failed to create mesh");
    for (i = 0; i < numAssocType; i++) {
      for (j = 0; j < numPerMesh; j++) {
	if (m == 1 && j == 1)
	  continue; // drop one of each kind of variable on middle mesh
	sprintf(temp_var_name, "%s %d", var_names[i], j);
	rc = fc_createSeqVariable(mesh, sequence, temp_var_name, &numStep, 
				  &seqVar);
	fail_unless(rc == FC_SUCCESS, "abort: failed to create variable");
	for (k = 0; k < numStep; k++) {
	  int uid = ((m*numAssocType + i)*numPerMesh + j)*numStep + k;
	  rc = fc_setVariableData(seqVar[k], numEntities[i], 1, assocTypes[i],
				  FC_MT_SCALAR, FC_DT_DOUBLE, &data[uid]);
	  fail_unless(rc == FC_SUCCESS, "abort: failed to set data on variable");
	}
 	free(seqVar);
      }
    }
  }
  
  // test writing
  rc = fc_writeDataset(dataset, EXODUS_FILE_NAME, FC_FT_EXODUS);
  fail_unless(rc == FC_SUCCESS, "should be able to write dataset");

  // test close & reopen
  rc = fc_deleteDataset(dataset);
  fail_unless(rc == FC_SUCCESS, "should close dataset w/o error");    
  rc = fc_loadDatasetWithFormat(EXODUS_FILE_NAME, FC_FT_EXODUS, &dataset);
  fail_unless(rc == FC_SUCCESS, "should load dataset");

  // check variables
  rc = fc_getMeshes(dataset, &temp_numMesh, &temp_meshes);
  fail_unless(rc == FC_SUCCESS, "should find meshes");
  fail_unless(temp_numMesh == numMesh, "mismatch of numMesh");
  for (m = 0; m < numMesh; m++) {
    rc = fc_getSeqVariables(temp_meshes[m], &temp_numSeqVar, &temp_numSteps,
			    &temp_seqVars);
    fail_unless(rc == FC_SUCCESS, "should find seq vars");
    if (m == 1)
      fail_unless(temp_numSeqVar == numAssocTypeWork*numPerMesh-1,
		  "mismatch of numSeqVar");
    else 
      fail_unless(temp_numSeqVar == numAssocTypeWork*numPerMesh, 
		  "mismatch of numSeqVar");
    for (i = 0; i < numAssocTypeWork; i++) {
      for (j = 0; j < numPerMesh; j++) {
	int idx = i*numPerMesh+j;
	if (m == 1 && assocTypes[assocTypeWorkLU[i]] == FC_AT_ELEMENT 
	    && j == 1)
	  j++; // skips the case that doesn't exist
	// check info
	fail_unless(temp_numSteps[idx] == numStep, "mismatch of numStep");
	rc = fc_getVariableName(temp_seqVars[idx][0], &temp_name);
	sprintf(temp_var_name, "%s %d", var_names[assocTypeWorkLU[i]], j);
	fail_unless(!strcmp(temp_name, temp_var_name), "mismatch of name");
	free(temp_name);
	rc = fc_getVariableInfo(temp_seqVars[idx][0], &temp_numDataPoint,
				&temp_numComp, &temp_assoc, &temp_mathtype,
				&temp_datatype);
	fail_unless(rc == FC_SUCCESS, "failed to get variable info");
	fail_unless(temp_numDataPoint == numEntities[assocTypeWorkLU[i]],
		    "mismatch of numDataPoint");
	fail_unless(temp_numComp == 1, "mismatch of numComp");
	fail_unless(temp_assoc == assocTypes[assocTypeWorkLU[i]],
		    "mismatch of assoc");
	fail_unless(temp_mathtype == FC_MT_SCALAR, "mismatch of mathtype");
	fail_unless(temp_datatype == FC_DT_DOUBLE, "mismatch of datatype");
	// check data
	for (k = 0; k < numStep; k++) {
	  int uid = ((m*numAssocType + assocTypeWorkLU[i])*numPerMesh + 
                     j)*numStep + k;
	  rc = fc_getVariableDataPtr(temp_seqVars[idx][k], &temp_data);
	  fail_unless(rc == FC_SUCCESS, "failed to get variable step data");
	  if (m == 1 && i == 0 && j == 1) { // the middle vertex variable
	    for (n = 0; n < temp_numDataPoint; n++) {
	      fail_unless(((double*)temp_data)[n] == 0., 
			  "mismatch of zero data");
	    }
	  }
	  else {
	    fail_unless(!memcmp(temp_data, &data[uid],
	       		temp_numDataPoint*sizeof(double)),
			"mismatch of data");
	  }
	}
	free(temp_seqVars[idx]);
      }
    }
    free(temp_numSteps);
    free(temp_seqVars);
  }
  free(temp_meshes);

  // all done
  rc = fc_deleteDataset(dataset);
  fail_unless(rc == FC_SUCCESS, "failed to close after testing"); 

  // ** Case 4: multiple of each association type on multiple meshes--
  //       BUT each has different mathtype

  // Create dataset, sequence, hex meshes & seqVars
  rc = fc_createDataset(ds_name, &dataset);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create dataset");
  rc = fc_createSequence(dataset, seq_name, &sequence);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create seqeunce");
  rc = fc_setSequenceCoords(sequence, numStep, FC_DT_DOUBLE, times);
  fail_unless(rc == FC_SUCCESS, "abort: failed to set seq coords");
  for (m = 0; m < numMesh; m++) {
    rc = _fc_makeMeshFromFiles(dataset, vertFiles[7], elemFiles[7], 
			       &mesh);
    fail_unless(rc == FC_SUCCESS, "abort: failed to create mesh");
    for (i = 0; i < numAssocType; i++) {
      for (j = 0; j < numPerMesh; j++) {
	sprintf(temp_var_name, "%s %d", var_names[i], j);
	rc = fc_createSeqVariable(mesh, sequence, temp_var_name, &numStep, 
				  &seqVar);
	fail_unless(rc == FC_SUCCESS, "abort: failed to create variable");
	for (k = 0; k < numStep; k++) {
	  int uid = ((m*numAssocType + i)*numPerMesh + j)*numStep + k;
	  rc = fc_setVariableData(seqVar[k], numEntities[i], numComps[j], 
				  assocTypes[i], mathTypes[j], FC_DT_DOUBLE, 
				  &data[uid]);
	  fail_unless(rc == FC_SUCCESS, "abort: failed to set data on variable");
	}
	free(seqVar);
      }
    }
  }

  // test writing
  rc = fc_writeDataset(dataset, EXODUS_FILE_NAME, FC_FT_EXODUS);
  fail_unless(rc == FC_SUCCESS, "should be able to write dataset");

  // test close & reopen
  rc = fc_deleteDataset(dataset);
  fail_unless(rc == FC_SUCCESS, "should close dataset w/o error");    
  rc = fc_loadDatasetWithFormat(EXODUS_FILE_NAME, FC_FT_EXODUS, &dataset);
  fail_unless(rc == FC_SUCCESS, "should load dataset");

  //FIXME: do a variable test where clipping occurs
  // check variables
  rc = fc_getMeshes(dataset, &temp_numMesh, &temp_meshes);
  fail_unless(rc == FC_SUCCESS, "should find meshes");
  fail_unless(temp_numMesh == numMesh, "mismatch of numMesh");
  for (m = 0; m < numMesh; m++) {
    int idx;
    char vector_endings[3][2] = { "x", "y", "z" };
    rc = fc_getSeqVariables(temp_meshes[m], &temp_numSeqVar, &temp_numSteps,
			    &temp_seqVars);
    fail_unless(rc == FC_SUCCESS, "should find seq vars");
    fail_unless(temp_numSeqVar == numAssocTypeWork*totalNumComp, 
		"mismatch of numSeqVar");
    idx = 0;
    for (i = 0; i < numAssocTypeWork; i++) {
      for (j = 0; j < numPerMesh; j++) {
	for (n = 0; n < numComps[j]; n++) {
	  // check info
	  fail_unless(temp_numSteps[idx] == numStep, "mismatch of numStep");
	  rc = fc_getVariableName(temp_seqVars[idx][0], &temp_name);
	  if (j == 0)
	    sprintf(temp_var_name, "%s %d", var_names[assocTypeWorkLU[i]], j);
	  else if (j == 1)
	    sprintf(temp_var_name, "%s %d_%s", var_names[assocTypeWorkLU[i]],
		    j, vector_endings[n]);
	  else
	    sprintf(temp_var_name, "%s %d_c%d", var_names[assocTypeWorkLU[i]],
		    j, n);
	  fail_unless(!strcmp(temp_name, temp_var_name), "mismatch of name");
	  free(temp_name);
	  rc = fc_getVariableInfo(temp_seqVars[idx][0], &temp_numDataPoint,
				&temp_numComp, &temp_assoc, &temp_mathtype,
				  &temp_datatype);
	  fail_unless(rc == FC_SUCCESS, "failed to get variable info");
	  fail_unless(temp_numDataPoint == numEntities[assocTypeWorkLU[i]],
		      "mismatch of numDataPoint");
	  fail_unless(temp_numComp == 1, "mismatch of numComp");
	  fail_unless(temp_assoc == assocTypes[assocTypeWorkLU[i]],
		      "mismatch of assoc");
	  fail_unless(temp_mathtype == FC_MT_SCALAR, "mismatch of mathtype");
	  fail_unless(temp_datatype == FC_DT_DOUBLE, "mismatch of datatype");
	  // check data
	  for (k = 0; k < numStep; k++) {
	    int uid = ((m*numAssocType + assocTypeWorkLU[i])*numPerMesh + 
                       j)*numStep + k;
	    rc = fc_getVariableDataPtr(temp_seqVars[idx][k], &temp_data);
	    fail_unless(rc == FC_SUCCESS, "failed to get variable step data");
	    if (j == 0) // the scalar ones
	      fail_unless(!memcmp(temp_data, &data[uid],
				  temp_numDataPoint*sizeof(double)),
			  "mismatch of scalar data");
	    else {
	      for (p = 0; p < temp_numDataPoint; p++)
		fail_unless(((double*)temp_data)[p] == 
			    data[uid + p*numComps[j] + n],
			    "mismatch of non scalar data");
	    }
	  }
	  free(temp_seqVars[idx]);
	  idx++;
	}
      }
    }
    free(temp_numSteps);
    free(temp_seqVars);
  }
  free(temp_meshes);

  // all done
  rc = fc_deleteDataset(dataset);
  fail_unless(rc == FC_SUCCESS, "failed to close after testing"); 

  // ** Case 5: multiple of each association type on multiple meshes, scalars
  //      BUT each has different data type

  // Create dataset, sequence, hex meshes & seqVars
  rc = fc_createDataset(ds_name, &dataset);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create dataset");
  rc = fc_createSequence(dataset, seq_name, &sequence);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create seqeunce");
  rc = fc_setSequenceCoords(sequence, numStep, FC_DT_DOUBLE, times);
  fail_unless(rc == FC_SUCCESS, "abort: failed to set seq coords");
  for (m = 0; m < numMesh; m++) {
    rc = _fc_makeMeshFromFiles(dataset, vertFiles[7], elemFiles[7], 
			       &mesh);
    fail_unless(rc == FC_SUCCESS, "abort: failed to create mesh");
    for (i = 0; i < numAssocType; i++) {
      for (j = 0; j < numPerMesh; j++) {
	sprintf(temp_var_name, "%s %d", var_names[i], j);
	rc = fc_createSeqVariable(mesh, sequence, temp_var_name, &numStep, 
				  &seqVar);
	fail_unless(rc == FC_SUCCESS, "abort: failed to create variable");
	for (k = 0; k < numStep; k++) {
	  int uid = ((m*numAssocType + i)*numPerMesh + j)*numStep + k;
	  if (j == 0)
	    data_p = &char_data[uid];
	  else if (j == 1)
	    data_p = &int_data[uid];
	  else if (j == 2)
	    data_p = &float_data[uid];
	  rc = fc_setVariableData(seqVar[k], numEntities[i], 1, assocTypes[i],
				  FC_MT_SCALAR, dataTypes[j], data_p);
	  fail_unless(rc == FC_SUCCESS, "abort: failed to set data on variable");
	}
	free(seqVar);
      }
    }
  }

  // test writing
  rc = fc_writeDataset(dataset, EXODUS_FILE_NAME, FC_FT_EXODUS);
  fail_unless(rc == FC_SUCCESS, "should be able to write dataset");

  // test close & reopen
  rc = fc_deleteDataset(dataset);
  fail_unless(rc == FC_SUCCESS, "should close dataset w/o error");    
  rc = fc_loadDatasetWithFormat(EXODUS_FILE_NAME, FC_FT_EXODUS, &dataset);
  fail_unless(rc == FC_SUCCESS, "should load dataset");

  // check variables
  rc = fc_getMeshes(dataset, &temp_numMesh, &temp_meshes);
  fail_unless(rc == FC_SUCCESS, "should find meshes");
  fail_unless(temp_numMesh == numMesh, "mismatch of numMesh");
  for (m = 0; m < numMesh; m++) {
    rc = fc_getSeqVariables(temp_meshes[m], &temp_numSeqVar, &temp_numSteps,
			    &temp_seqVars);
    fail_unless(rc == FC_SUCCESS, "should find seq vars");
    fail_unless(temp_numSeqVar == numAssocTypeWork*numPerMesh, 
		"mismatch of numSeqVar");
    for (i = 0; i < numAssocTypeWork; i++) {
      for (j = 0; j < numPerMesh; j++) {
	int idx = i*numPerMesh+j;
	// check info
	fail_unless(temp_numSteps[idx] == numStep, "mismatch of numStep");
	rc = fc_getVariableName(temp_seqVars[idx][0], &temp_name);
	sprintf(temp_var_name, "%s %d", var_names[assocTypeWorkLU[i]], j);
	fail_unless(!strcmp(temp_name, temp_var_name), "mismatch of name");
	free(temp_name);
	rc = fc_getVariableInfo(temp_seqVars[idx][0], &temp_numDataPoint,
				&temp_numComp, &temp_assoc, &temp_mathtype,
				&temp_datatype);
	fail_unless(rc == FC_SUCCESS, "failed to get variable info");
	fail_unless(temp_numDataPoint == numEntities[assocTypeWorkLU[i]],
		    "mismatch of numDataPoint");
	fail_unless(temp_numComp == 1, "mismatch of numComp");
	fail_unless(temp_assoc == assocTypes[assocTypeWorkLU[i]],
		    "mismatch of assoc");
	fail_unless(temp_mathtype == FC_MT_SCALAR, "mismatch of mathtype");
	fail_unless(temp_datatype == FC_DT_DOUBLE, "mismatch of datatype");
	// check data
	for (k = 0; k < numStep; k++) {
	  int uid = ((m*numAssocType + assocTypeWorkLU[i])*numPerMesh + 
                     j)*numStep + k;
	  rc = fc_getVariableDataPtr(temp_seqVars[idx][k], &temp_data);
	  fail_unless(rc == FC_SUCCESS, "failed to get variable step data");
	  for (n = 0; n < temp_numDataPoint; n++) {
	    if (j == 0) 
	      fail_unless(((double*)temp_data)[n] == char_data[uid+n],
			  "mismatch of char data");
	    else if (j == 1) 
	      fail_unless(((double*)temp_data)[n] == int_data[uid+n],
			  "mismatch of int data");
	    else if (j == 2)
	      fail_unless(((double*)temp_data)[n] == float_data[uid+n],
			  "mismatch of float data");
	  }
	}
	free(temp_seqVars[idx]);
      }
    }
    free(temp_numSteps);
    free(temp_seqVars);
  }
  free(temp_meshes);

  // all done
  rc = fc_deleteDataset(dataset);
  fail_unless(rc == FC_SUCCESS, "failed to close after testing"); 
}
END_TEST


// This is a basically a rewrite of exo_seq_var, without
// order reliance, but for var, however
// this is w/o FC_AT_WHOLE_MESH assoc since that is checked
// elsewhere. this checks attributes, the other checks properties.
START_TEST(exo_var)
{
  FC_ReturnCode rc;
  int i, j, k, m, n, p;
  FC_Dataset dataset;
  int numMesh = 3, numPerMesh = 3;
  FC_Mesh mesh, *temp_meshes;
  FC_Variable *temp_Vars;
  char* ds_name = "dataset";
  char var_names[numAssocType][1024], temp_var_name[1024];
  int numEntities[numAssocType], maxNumEntity = 2945; // for a hex
  double data[4*2945];
  char char_data[4*2945];
  int int_data[4*2945];
  float float_data[4*2945];
  void *data_p;
  // Not all assoc types work
  int numAssocTypeWork = 2;
  int assocTypeWorkLU[2] = { 0, 3 }; 
  int foundAssocTypes[numMesh][numAssocType][numPerMesh];
  int temp_numMesh, temp_numVar;
  int temp_numDataPoint, temp_numComp;
  char* temp_name;
  FC_AssociationType temp_assoc;
  FC_MathType temp_mathtype;
  FC_DataType temp_datatype;
  void* temp_data;
  FC_MathType mathTypes[3] = { FC_MT_SCALAR, FC_MT_VECTOR, FC_MT_TENSOR };
  int numComps[3] = { 1, 3, 2 }, totalNumComp = 1+3+2;
  FC_DataType dataTypes[3] = { FC_DT_CHAR, FC_DT_INT, FC_DT_FLOAT };

  // get unique data starting at different point inside of data
  // for example, timestep 0: &data[0], timestep 1: &data[1] ...
  for (i = 0; i < 4*maxNumEntity; i++) {
    data[i] = 1 + i/100.;
    char_data[i] = (char)data[i];
    int_data[i] = (int)data[i];
    float_data[i] = (float)data[i];
  }

  // ** Case 1: one of each association type on a single mesh, scalars

  // Create dataset, hex mesh & Vars
  rc = fc_createDataset(ds_name, &dataset);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create dataset");
  rc = _fc_makeMeshFromFiles(dataset, vertFiles[7], elemFiles[7], 
			     &mesh);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create mesh");
  // save numEnities for repeated use
  for (i = 0; i < numAssocType; i++) {
    fc_getMeshNumEntity(mesh, assocTypes[i], &numEntities[i]);
    if (numEntities[i] > maxNumEntity)
      maxNumEntity = numEntities[i];
  }
  for (i = 0; i < numAssocType; i++) {
    FC_Variable var;
    foundAssocTypes[0][i][0] = 0;

    if (assocTypes[i] == FC_AT_WHOLE_MESH){
      continue; // this is checked elsewhere
    }
    sprintf(var_names[i], "variable on %s",
	    fc_getAssociationTypeText(assocTypes[i]));
    rc = fc_createVariable(mesh, var_names[i], &var);
    fail_unless(rc == FC_SUCCESS, "abort: failed to create variable");
    //NOTE that FC_AT_WHOLE_MESH is only supported in exodus as integer
    rc = fc_setVariableData(var, numEntities[i], 1, assocTypes[i],
			    FC_MT_SCALAR, FC_DT_DOUBLE, 
			    &data[i]);
    fail_unless(rc == FC_SUCCESS, "abort: failed to set data on variable");
  }

  // test writing
  rc = fc_writeDataset(dataset, EXODUS_FILE_NAME, FC_FT_EXODUS);
  fail_unless(rc == FC_SUCCESS, "should be able to write dataset");
  
  // test close & reopen
  rc = fc_deleteDataset(dataset);
  fail_unless(rc == FC_SUCCESS, "should close dataset w/o error");    
  rc = fc_loadDatasetWithFormat(EXODUS_FILE_NAME, FC_FT_EXODUS, &dataset);
  fail_unless(rc == FC_SUCCESS, "should load dataset");

  // check variables
  rc = fc_getMeshes(dataset, &temp_numMesh, &temp_meshes);
  fail_unless(rc == FC_SUCCESS, "should find meshes");
  fail_unless(temp_numMesh == 1, "mismatch of numMesh");
  mesh = temp_meshes[0];
  free(temp_meshes);
  //delete the ID variable that exodus creates
  rc = fc_getVariableByName(mesh, "ID", &temp_numVar, &temp_Vars);
  fail_unless(rc == FC_SUCCESS, "cant get var by name");
  fail_unless(temp_numVar == 1, "should be 1 ID var");
  rc = fc_deleteVariable(temp_Vars[0]);
  fail_unless(rc == FC_SUCCESS, "cant delete ID var");
  free(temp_Vars);
  rc = fc_getVariables(mesh, &temp_numVar, &temp_Vars);
  fail_unless(rc == FC_SUCCESS, "should find vars");
  fail_unless(temp_numVar == numAssocTypeWork, "mismatch of num vars");
  for (i = 0; i < temp_numVar; i++) {
    // check info
    rc = fc_getVariableName(temp_Vars[i], &temp_name);
    rc = fc_getVariableInfo(temp_Vars[i], &temp_numDataPoint,
			    &temp_numComp, &temp_assoc, &temp_mathtype,
			    &temp_datatype);
    fail_unless(rc == FC_SUCCESS, "failed to get variable info");
    for ( j = 0; j < numAssocTypeWork; j++){
      //only one of each assoc
      if (temp_assoc == assocTypes[assocTypeWorkLU[j]]){
	foundAssocTypes[0][assocTypeWorkLU[j]][0]++;
	fail_unless(!strcmp(temp_name, var_names[assocTypeWorkLU[j]]), 
		    "mismatch of name");
	fail_unless(temp_numDataPoint == numEntities[assocTypeWorkLU[j]],
		    "mismatch of numDataPoint");
	fail_unless(temp_numComp == 1, "mismatch of numComp");
	fail_unless(temp_mathtype == FC_MT_SCALAR, "mismatch of mathtype");
	fail_unless(temp_datatype == FC_DT_DOUBLE, "mismatch of datatype");
	// check data
	rc = fc_getVariableDataPtr(temp_Vars[i], &temp_data);
	fail_unless(rc == FC_SUCCESS, "failed to get variable step data");
	fail_unless(!memcmp(temp_data, &data[assocTypeWorkLU[j]],
			    temp_numDataPoint*sizeof(double)),
		    "mismatch of data");
	break;
      }
    }
    free(temp_name);
  }
  free(temp_Vars);

  for (i = 0; i < numAssocTypeWork; i++){
    fail_unless(foundAssocTypes[0][assocTypeWorkLU[i]][0] == 1,
		"did not find var of correct assoc type");
  }


  // all done
  rc = fc_deleteDataset(dataset);
  fail_unless(rc == FC_SUCCESS, "failed to close after testing"); 


  // ** Case 2: multiple of each association type on multiple meshes, scalars

  // Create dataset, hex meshes & Vars
  rc = fc_createDataset(ds_name, &dataset);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create dataset");
  for (m = 0; m < numMesh; m++) {
    rc = _fc_makeMeshFromFiles(dataset, vertFiles[7], elemFiles[7], 
			       &mesh);
    fail_unless(rc == FC_SUCCESS, "abort: failed to create mesh");
    for (i = 0; i < numAssocType; i++) {
      if (assocTypes[i] == FC_AT_WHOLE_MESH){
	continue; // this is checked elsewhere
      }

      for (j = 0; j < numPerMesh; j++) {
	FC_Variable var;
	foundAssocTypes[m][i][j] = 0;

	sprintf(temp_var_name, "%s %d", var_names[i], j);
	rc = fc_createVariable(mesh, temp_var_name, &var);
	fail_unless(rc == FC_SUCCESS, "abort: failed to create variable");
	int uid = (m*numAssocType + i)*numPerMesh+j;
	rc = fc_setVariableData(var, numEntities[i], 1, assocTypes[i],
				FC_MT_SCALAR, FC_DT_DOUBLE, &data[uid]);
	fail_unless(rc == FC_SUCCESS, "abort: failed to set data on variable");
      }
    }
  }
  
  // test writing
  rc = fc_writeDataset(dataset, EXODUS_FILE_NAME, FC_FT_EXODUS);
  fail_unless(rc == FC_SUCCESS, "should be able to write dataset");

  // test close & reopen
  rc = fc_deleteDataset(dataset);
  fail_unless(rc == FC_SUCCESS, "should close dataset w/o error");    
  rc = fc_loadDatasetWithFormat(EXODUS_FILE_NAME, FC_FT_EXODUS, &dataset);
  fail_unless(rc == FC_SUCCESS, "should load dataset");

  // check variables
  rc = fc_getMeshes(dataset, &temp_numMesh, &temp_meshes);
  fail_unless(rc == FC_SUCCESS, "should find meshes");
  fail_unless(temp_numMesh == numMesh, "mismatch of numMesh");
  for (m = 0; m < numMesh; m++) {
    //delete the ID variable that exodus creates
    rc = fc_getVariableByName(temp_meshes[m], "ID", &temp_numVar, &temp_Vars);
    fail_unless(rc == FC_SUCCESS, "cant get var by name");
    fail_unless(temp_numVar == 1, "should be 1 ID var");
    rc = fc_deleteVariable(temp_Vars[0]);
    fail_unless(rc == FC_SUCCESS, "cant delete ID var");
    free(temp_Vars);
    rc = fc_getVariables(temp_meshes[m], &temp_numVar, &temp_Vars);
    fail_unless(rc == FC_SUCCESS, "should find vars");
    fail_unless(temp_numVar == numAssocTypeWork*numPerMesh, 
		"mismatch of numVar");
    for (i = 0; i < temp_numVar; i++) {
      // check info
      rc = fc_getVariableName(temp_Vars[i], &temp_name);
      rc = fc_getVariableInfo(temp_Vars[i], &temp_numDataPoint,
			      &temp_numComp, &temp_assoc, &temp_mathtype,
			      &temp_datatype);
      fail_unless(rc == FC_SUCCESS, "failed to get variable info");
      for ( j = 0; j < numAssocTypeWork; j++){
	if (temp_assoc == assocTypes[assocTypeWorkLU[j]]){
	  int uid;
	  // there are numPerMesh and they end with different digits
	  //this doesnt check the ending digit....
	  fail_unless(!strncmp(temp_name, var_names[assocTypeWorkLU[j]],
			       strlen(var_names[assocTypeWorkLU[j]])),
		      "mismatch of name");
	  uid = atoi(&(temp_name[strlen(temp_name)-1])); //does not detect errors
	  fail_unless((uid < numPerMesh && uid >= 0), "bad number");
	  foundAssocTypes[m][assocTypeWorkLU[j]][uid]++;
	  fail_unless(temp_numDataPoint == numEntities[assocTypeWorkLU[j]],
		      "mismatch of numDataPoint");
	  fail_unless(temp_numComp == 1, "mismatch of numComp");
	  fail_unless(temp_mathtype == FC_MT_SCALAR, "mismatch of mathtype");
	  fail_unless(temp_datatype == FC_DT_DOUBLE, "mismatch of datatype");
	  // check data
	  rc = fc_getVariableDataPtr(temp_Vars[i], &temp_data);
	  fail_unless(rc == FC_SUCCESS, "failed to get variable data");
	  uid += (m*numAssocType + assocTypeWorkLU[j])*numPerMesh;
	  fail_unless(!memcmp(temp_data, &data[uid],
			      temp_numDataPoint*sizeof(double)),
		      "mismatch of data");
	  break;
	}
      }
      free(temp_name);
    }
    free(temp_Vars);

    for (i = 0; i < numAssocType; i++){
      int shouldwork = 0;
      if (assocTypes[i] == FC_AT_WHOLE_MESH){
	continue;
      }
      for (k = 0; k < numAssocTypeWork; k++){
	if (assocTypeWorkLU[k] == i){
	  shouldwork = 1;
	  break;
	}
      }
      for (j = 0; j < numPerMesh; j++){
	fail_unless(foundAssocTypes[m][i][j] == shouldwork,
		    "prob with found/unfound var");
      }
    }
  }


  free(temp_meshes);

  // all done
  rc = fc_deleteDataset(dataset);
  fail_unless(rc == FC_SUCCESS, "failed to close after testing"); 

  // ** Case 3: multiple of each association type on multiple meshes, scalars
  //       BUT don't have complete coverage (leave one off of middle mesh)

  // Create dataset, hex meshes & Vars
  rc = fc_createDataset(ds_name, &dataset);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create dataset");
  for (m = 0; m < numMesh; m++) {
    rc = _fc_makeMeshFromFiles(dataset, vertFiles[7], elemFiles[7], 
			       &mesh);
    fail_unless(rc == FC_SUCCESS, "abort: failed to create mesh");
    for (i = 0; i < numAssocType; i++) {
      if (assocTypes[i] == FC_AT_WHOLE_MESH){
	continue; // this is checked elsewhere
      }

      for (j = 0; j < numPerMesh; j++) {
	FC_Variable var;
	foundAssocTypes[m][i][j] = 0;

	if (m != 1){
	  sprintf(temp_var_name, "%s %d", var_names[i], j);
	  rc = fc_createVariable(mesh, temp_var_name, &var);
	  fail_unless(rc == FC_SUCCESS, "abort: failed to create variable");
	  int uid = (m*numAssocType + i)*numPerMesh+j;
	  rc = fc_setVariableData(var, numEntities[i], 1, assocTypes[i],
				  FC_MT_SCALAR, FC_DT_DOUBLE, &data[uid]);
	  fail_unless(rc == FC_SUCCESS, "abort: failed to set data on variable");
	}
      }
    }
  }
  
  // test writing
  rc = fc_writeDataset(dataset, EXODUS_FILE_NAME, FC_FT_EXODUS);
  fail_unless(rc == FC_SUCCESS, "should be able to write dataset");

  // test close & reopen
  rc = fc_deleteDataset(dataset);
  fail_unless(rc == FC_SUCCESS, "should close dataset w/o error");    
  rc = fc_loadDatasetWithFormat(EXODUS_FILE_NAME, FC_FT_EXODUS, &dataset);
  fail_unless(rc == FC_SUCCESS, "should load dataset");

  // check variables
  rc = fc_getMeshes(dataset, &temp_numMesh, &temp_meshes);
  fail_unless(rc == FC_SUCCESS, "should find meshes");
  fail_unless(temp_numMesh == numMesh, "mismatch of numMesh");
  for (m = 0; m < numMesh; m++) {
    //delete the ID variable that exodus creates
    rc = fc_getVariableByName(temp_meshes[m], "ID", &temp_numVar, &temp_Vars);
    fail_unless(rc == FC_SUCCESS, "cant get var by name");
    fail_unless(temp_numVar == 1, "should be 1 ID var");
    rc = fc_deleteVariable(temp_Vars[0]);
    fail_unless(rc == FC_SUCCESS, "cant delete ID var");
    free(temp_Vars);
    rc = fc_getVariables(temp_meshes[m], &temp_numVar, &temp_Vars);
    fail_unless(rc == FC_SUCCESS, "should find vars");
    if (m == 1){
      int numVertexVars = numPerMesh;       //get vertex ones by default
      fail_unless(temp_numVar == numVertexVars, "mismatch of numVar");
    } else {
      fail_unless(temp_numVar == numAssocTypeWork*numPerMesh, 
		  "mismatch of numVar");
    }
    for (i = 0; i < temp_numVar; i++) {
      // check info
      rc = fc_getVariableName(temp_Vars[i], &temp_name);
      rc = fc_getVariableInfo(temp_Vars[i], &temp_numDataPoint,
			      &temp_numComp, &temp_assoc, &temp_mathtype,
			      &temp_datatype);
      fail_unless(rc == FC_SUCCESS, "failed to get variable info");
      for ( j = 0; j < numAssocTypeWork; j++){
	if (temp_assoc == assocTypes[assocTypeWorkLU[j]]){
	  int uid;
	  // there are numPerMesh and they end with different digits
	  //this doesnt check the ending digit....
	  fail_unless(!strncmp(temp_name, var_names[assocTypeWorkLU[j]],
			       strlen(var_names[assocTypeWorkLU[j]])),
		      "mismatch of name");
	  uid = atoi(&(temp_name[strlen(temp_name)-1])); //does not detect errors
	  fail_unless((uid < numPerMesh && uid >= 0), "bad number");
	  foundAssocTypes[m][assocTypeWorkLU[j]][uid]++;
	  fail_unless(temp_numDataPoint == numEntities[assocTypeWorkLU[j]],
		      "mismatch of numDataPoint");
	  fail_unless(temp_numComp == 1, "mismatch of numComp");
	  fail_unless(temp_mathtype == FC_MT_SCALAR, "mismatch of mathtype");
	  fail_unless(temp_datatype == FC_DT_DOUBLE, "mismatch of datatype");
	  // check data
	  rc = fc_getVariableDataPtr(temp_Vars[i], &temp_data);
	  fail_unless(rc == FC_SUCCESS, "failed to get variable data");
	  if (m!= 1){
	    uid += (m*numAssocType + assocTypeWorkLU[j])*numPerMesh;
	    fail_unless(!memcmp(temp_data, &data[uid],
				temp_numDataPoint*sizeof(double)),
			"mismatch of data");
	  } else {
	    //data should be all zeros by default
	    for (k = 0; k < temp_numDataPoint; k++){
	      double temp = ((double*)temp_data)[k];
	      fail_unless(FC_DBL_EQUIV(temp,0), "mismatch of data");
	    }
	  }
	  break;
	}
      }
      free(temp_name);
    }
    free(temp_Vars);

    for (i = 0; i < numAssocType; i++){
      int shouldwork = 0;
      if (assocTypes[i] == FC_AT_WHOLE_MESH){
	continue;
      }
      if (m != 1){
	for (k = 0; k < numAssocTypeWork; k++){
	  if (assocTypeWorkLU[k] == i){
	    shouldwork = 1;
	    break;
	  }
	}
      } else {
	if (assocTypes[i] == FC_AT_VERTEX){
	  shouldwork = 1;
	} 
      }

      for (j = 0; j < numPerMesh; j++){
	fail_unless(foundAssocTypes[m][i][j] == shouldwork,
		    "prob with found/unfound var");
      }
    }
  }

  free(temp_meshes);

  // all done
  rc = fc_deleteDataset(dataset);
  fail_unless(rc == FC_SUCCESS, "failed to close after testing"); 

  // ** Case 4: multiple of each association type on multiple meshes--
  //       BUT each has different mathtype

  // Create dataset, hex meshes & Vars
  rc = fc_createDataset(ds_name, &dataset);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create dataset");
  for (m = 0; m < numMesh; m++) {
    rc = _fc_makeMeshFromFiles(dataset, vertFiles[7], elemFiles[7], 
			       &mesh);
    fail_unless(rc == FC_SUCCESS, "abort: failed to create mesh");
    for (i = 0; i < numAssocType; i++) {
      if (assocTypes[i] == FC_AT_WHOLE_MESH){
	continue; // this is checked elsewhere
      }
      
      for (j = 0; j < numPerMesh; j++) {
	FC_Variable var;

	sprintf(temp_var_name, "%s %d", var_names[i], j);
	rc = fc_createVariable(mesh, temp_var_name, &var);
	fail_unless(rc == FC_SUCCESS, "abort: failed to create variable");
	int uid = (m*numAssocType + i)*3+j;	
	rc = fc_setVariableData(var, numEntities[i], numComps[j], 
				assocTypes[i], mathTypes[j], FC_DT_DOUBLE, 
				&data[uid]);
	fail_unless(rc == FC_SUCCESS, "abort: failed to set data on variable");
      }
    }
  }
  
  // test writing
  rc = fc_writeDataset(dataset, EXODUS_FILE_NAME, FC_FT_EXODUS);
  fail_unless(rc == FC_SUCCESS, "should be able to write dataset");

  // test close & reopen
  rc = fc_deleteDataset(dataset);
  fail_unless(rc == FC_SUCCESS, "should close dataset w/o error");    
  rc = fc_loadDatasetWithFormat(EXODUS_FILE_NAME, FC_FT_EXODUS, &dataset);
  fail_unless(rc == FC_SUCCESS, "should load dataset");

  // check variables
  rc = fc_getMeshes(dataset, &temp_numMesh, &temp_meshes);
  fail_unless(rc == FC_SUCCESS, "should find meshes");
  fail_unless(temp_numMesh == numMesh, "mismatch of numMesh");
  for (m = 0; m < numMesh; m++) {
    char vector_endings[3][2] = { "x", "y", "z" };

    //delete the ID variable that exodus creates
    rc = fc_getVariableByName(temp_meshes[m], "ID", &temp_numVar, &temp_Vars);
    fail_unless(rc == FC_SUCCESS, "cant get var by name");
    fail_unless(temp_numVar == 1, "should be 1 ID var");
    rc = fc_deleteVariable(temp_Vars[0]);
    fail_unless(rc == FC_SUCCESS, "cant delete ID var");
    free(temp_Vars);

    rc = fc_getVariables(temp_meshes[m], &temp_numVar, &temp_Vars);
    fail_unless(rc == FC_SUCCESS, "should find vars");
    fail_unless(temp_numVar == numAssocTypeWork*totalNumComp, 
		"mismatch of numVar");

    //changing the way we are iterating through them
    //because it is becoming too unweidly. we know
    //which ones we should have, now ask for them by name
    for (i = 0; i < numAssocTypeWork; i++){
      for (j = 0; j < 3; j++){
	for (k = 0; k < numComps[j]; k++){
	  FC_Variable* lVars;
	  if (j == 0){
	    sprintf(temp_var_name, "%s %d", var_names[assocTypeWorkLU[i]], j);
	  } else{
	    if (j == 1){
	      sprintf(temp_var_name, "%s %d_%s", var_names[assocTypeWorkLU[i]],
		      j, vector_endings[k]);
	    } else {
	      sprintf(temp_var_name, "%s %d_c%d", var_names[assocTypeWorkLU[i]],
		      j, k);
	    }
	  }
	  rc = fc_getVariableByName(temp_meshes[m], temp_var_name,
				    &temp_numVar, &lVars);
	  fail_unless(rc == FC_SUCCESS, "cant get var by name");
	  fail_unless(temp_numVar == 1, "wrong num return vars");
	  rc = fc_getVariableInfo(lVars[0], &temp_numDataPoint,
				  &temp_numComp, &temp_assoc, &temp_mathtype,
				  &temp_datatype);
	  fail_unless(rc == FC_SUCCESS, "failed to get variable info");
	  //everyone ends up this way
	  fail_unless(temp_numComp == 1, "mismatch of numComp");
	  fail_unless(temp_mathtype == FC_MT_SCALAR, "mismatch of mathtype");
	  fail_unless(temp_datatype == FC_DT_DOUBLE, "mismatch of datatype");
	  fail_unless(temp_numDataPoint == numEntities[assocTypeWorkLU[i]],
		      "mismatch of numDataPoint");
	  // check data
	  rc = fc_getVariableDataPtr(lVars[0], &temp_data);
	  fail_unless(rc == FC_SUCCESS, "failed to get variable data");

	  int uid = (m*numAssocType + assocTypeWorkLU[i])*3+j;	
	  if (j == 0) // the scalar ones
	    fail_unless(!memcmp(temp_data, &data[uid],
				temp_numDataPoint*sizeof(double)),
			"mismatch of scalar data");
	  else {
	    for (p = 0; p < temp_numDataPoint; p++)
	      fail_unless(((double*)temp_data)[p] == 
			  data[uid + p*numComps[j] + k],
			  "mismatch of non scalar data");
	  }
	  if (lVars) free(lVars);
	}
      }
    }
    free(temp_Vars);
  }


  free(temp_meshes);

  // all done
  rc = fc_deleteDataset(dataset);
  fail_unless(rc == FC_SUCCESS, "failed to close after testing"); 


  // ** Case 5: multiple of each association type on multiple meshes, scalars
  //      BUT each has different data type

  // Create dataset, hex meshes & Vars
  rc = fc_createDataset(ds_name, &dataset);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create dataset");
  for (m = 0; m < numMesh; m++) {
    rc = _fc_makeMeshFromFiles(dataset, vertFiles[7], elemFiles[7], 
			       &mesh);
    fail_unless(rc == FC_SUCCESS, "abort: failed to create mesh");
    for (i = 0; i < numAssocType; i++) {
      if (assocTypes[i] == FC_AT_WHOLE_MESH){
	continue; // this is checked elsewhere
      }
      
      for (j = 0; j < numPerMesh; j++) {
	FC_Variable var;

	foundAssocTypes[m][i][j] = 0;

	sprintf(temp_var_name, "%s %d", var_names[i], j);
	rc = fc_createVariable(mesh, temp_var_name, &var);
	fail_unless(rc == FC_SUCCESS, "abort: failed to create variable");
	int uid = (m*numAssocType + i)*3+j;	
	if (j == 0)
	  data_p = &char_data[uid];
	else if (j == 1)
	  data_p = &int_data[uid];
	else if (j == 2)
	  data_p = &float_data[uid];
	rc = fc_setVariableData(var, numEntities[i], 1, assocTypes[i],
				FC_MT_SCALAR, dataTypes[j], data_p);
	fail_unless(rc == FC_SUCCESS, "abort: failed to set data on variable");
      }
    }
  }
  
  // test writing
  rc = fc_writeDataset(dataset, EXODUS_FILE_NAME, FC_FT_EXODUS);
  fail_unless(rc == FC_SUCCESS, "should be able to write dataset");

  // test close & reopen
  rc = fc_deleteDataset(dataset);
  fail_unless(rc == FC_SUCCESS, "should close dataset w/o error");    
  rc = fc_loadDatasetWithFormat(EXODUS_FILE_NAME, FC_FT_EXODUS, &dataset);
  fail_unless(rc == FC_SUCCESS, "should load dataset");

  // check variables
  rc = fc_getMeshes(dataset, &temp_numMesh, &temp_meshes);
  fail_unless(rc == FC_SUCCESS, "should find meshes");
  fail_unless(temp_numMesh == numMesh, "mismatch of numMesh");
  for (m = 0; m < numMesh; m++) {
    //delete the ID variable that exodus creates
    rc = fc_getVariableByName(temp_meshes[m], "ID", &temp_numVar, &temp_Vars);
    fail_unless(rc == FC_SUCCESS, "cant get var by name");
    fail_unless(temp_numVar == 1, "should be 1 ID var");
    rc = fc_deleteVariable(temp_Vars[0]);
    fail_unless(rc == FC_SUCCESS, "cant delete ID var");
    free(temp_Vars);
    rc = fc_getVariables(temp_meshes[m], &temp_numVar, &temp_Vars);
    fail_unless(rc == FC_SUCCESS, "should find vars");
    fail_unless(temp_numVar == numAssocTypeWork*numPerMesh, 
		"mismatch of numVar");
    for (i = 0; i < temp_numVar; i++) {
      // check info
      rc = fc_getVariableName(temp_Vars[i], &temp_name);
      rc = fc_getVariableInfo(temp_Vars[i], &temp_numDataPoint,
			      &temp_numComp, &temp_assoc, &temp_mathtype,
			      &temp_datatype);
      fail_unless(rc == FC_SUCCESS, "failed to get variable info");
      for ( j = 0; j < numAssocTypeWork; j++){
	if (temp_assoc == assocTypes[assocTypeWorkLU[j]]){
	  int uid, nuid;
	  // there are numPerMesh and they end with different digits
	  //this doesnt check the ending digit....
	  fail_unless(!strncmp(temp_name, var_names[assocTypeWorkLU[j]],
			       strlen(var_names[assocTypeWorkLU[j]])),
		      "mismatch of name");
	  uid = atoi(&(temp_name[strlen(temp_name)-1])); //does not detect errors
	  fail_unless((uid < numPerMesh && uid >= 0), "bad number");
	  foundAssocTypes[m][assocTypeWorkLU[j]][uid]++;
	  fail_unless(temp_numDataPoint == numEntities[assocTypeWorkLU[j]],
		      "mismatch of numDataPoint");
	  fail_unless(temp_numComp == 1, "mismatch of numComp");
	  fail_unless(temp_mathtype == FC_MT_SCALAR, "mismatch of mathtype");
	  fail_unless(temp_datatype == FC_DT_DOUBLE, "mismatch of datatype");
	  // check data
	  rc = fc_getVariableDataPtr(temp_Vars[i], &temp_data);
	  fail_unless(rc == FC_SUCCESS, "failed to get variable data");
	  nuid = (m*numAssocType + assocTypeWorkLU[j])*3+uid;
	  for (n = 0; n < temp_numDataPoint; n++) {
	    if (uid == 0) 
	      fail_unless(((double*)temp_data)[n] == char_data[nuid+n],
			  "mismatch of char data");
	    else if (uid == 1) 
	      fail_unless(((double*)temp_data)[n] == int_data[nuid+n],
			  "mismatch of int data");
	    else if (uid == 2)
	      fail_unless(((double*)temp_data)[n] == float_data[nuid+n],
			  "mismatch of float data");
	  }
	}
      }
      free(temp_name);
    }
    free(temp_Vars);

    for (i = 0; i < numAssocType; i++){
      int shouldwork = 0;
      if (assocTypes[i] == FC_AT_WHOLE_MESH){
	continue;
      }
      for (k = 0; k < numAssocTypeWork; k++){
	if (assocTypeWorkLU[k] == i){
	  shouldwork = 1;
	  break;
	}
      }
      for (j = 0; j < numPerMesh; j++){
	fail_unless(foundAssocTypes[m][i][j] == shouldwork,
		    "prob with found/unfound var");
      }
    }
  }


  free(temp_meshes);

  // all done
  rc = fc_deleteDataset(dataset);
  fail_unless(rc == FC_SUCCESS, "failed to close after testing"); 

}
END_TEST


// Make sure functions return the proper info
START_TEST(exo_global_ids)
{
  FC_ReturnCode rc;
  int i, j;
  FC_Dataset dataset;
  FC_Mesh meshes[8];
  char* ds_name = "dataset", *mesh_name;
  int numVerts[8], numElems[8], numDim, temp_numVert;
  int* conns;
  double* coords;
  FC_ElementType elemType;
  _FC_DsSlot* dsSlot;
  int globalNodeID, globalElemID;

  // Just need to make a temporary dataset & write it to make it Exodus-y
  // Needs to have multiple meshes

  // Create dataset
  rc = fc_createDataset(ds_name, &dataset);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create dataset");
  dsSlot = _fc_getDsSlot(dataset);
  fail_unless(dsSlot != NULL, "abort: failed to get dsSlot");

  // Add the meshes
  for (i = 0; i < numMeshType; i++) {
    
    // Setup
    rc = _fc_readVertexFile(vertFiles[i], &numVerts[i], &numDim, &coords);
    fail_unless(rc == FC_SUCCESS, "abort: failed to read vert file");
    fail_unless(numDim == 3, "abort: expect all to have numDim = 3");
    rc = _fc_readElementFile(elemFiles[i], &temp_numVert, &numElems[i], 
			     &elemType, &mesh_name, &conns);
    fail_unless(rc == FC_SUCCESS, "abort: failed to read elem file");
    fail_unless(temp_numVert == numVerts[i], "abort: numVerts don't agree");

    // Create the mesh
    rc = fc_createMesh(dataset, mesh_name, &meshes[i]);
    free(mesh_name);
    fail_unless(rc == FC_SUCCESS, "abort: failed to create mesh");
    rc = fc_setMeshCoordsPtr(meshes[i], numDim, numVerts[i], coords);
    fail_unless(rc == FC_SUCCESS, "abort: failed to set coords");
    rc = fc_setMeshElementConnsPtr(meshes[i], elemType, numElems[i], conns);
    fail_unless(rc == FC_SUCCESS, "abort: failed to set conns");
  }

  // write it
  rc = fc_writeDataset(dataset, EXODUS_FILE_NAME, FC_FT_EXODUS);
  fail_unless(rc == FC_SUCCESS, "abort: faile to write dataset");
    
  // test close & reopen
  //rc = fc_deleteDataset(dataset);
  //fail_unless(rc == FC_SUCCESS, "should close dataset w/o error");
  //rc = fc_loadDatasetWithFormat(EXODUS_FILE_NAME, FC_FT_EXODUS, &dataset);
  //fail_unless(rc == FC_SUCCESS, "should load dataset");

  // Testing 
  globalNodeID = 0;
  globalElemID = 0;
  for (i = 0; i < numMeshType; i++) {
    // setup
    int meshExoID, *LUT, offset;
    _FC_MeshSlot *meshSlot = _fc_getMeshSlot(meshes[i]);
    fail_unless(meshSlot != NULL, "abort: failed to get meshSlot");
    meshExoID = meshSlot->fileInfo.exoInfo.meshID;
    // This doesn't have to be true, but we assume it is in the testing below
    // With a newly created dataset, global IDs will increase monotonically
    fail_unless(meshExoID == i, "abort: exoMeshID is not as expected");
    
    // node ids
    rc = _fc_getMeshExodusGlobalNodalIDsPtr(meshes[i], &LUT);
    fail_unless(rc == FC_SUCCESS, "failed to get mesh exo global node ids");
    fail_unless(LUT == dsSlot->fileInfo.exoInfo.vertIDs[meshExoID],
		"mismatch of LUT");
    for (j = 0; j < numVerts[i]; j++) {
      fail_unless(LUT[j] == globalNodeID, "mismatch of global node ID");
      globalNodeID++;
    }

    // elem ids
    rc = _fc_getMeshExodusGlobalElementIDOffset(meshes[i], &offset);
    fail_unless(rc == FC_SUCCESS, "failed to get mesh exo global elem offset");
    fail_unless(offset == meshSlot->fileInfo.exoInfo.elemOffset,
		"mismatch of global mesh offset");
    fail_unless(offset == globalElemID, "mismatch of global elem ID");
    globalElemID += numElems[i];
  }

  // all done
  rc = fc_deleteDataset(dataset);
  fail_unless(rc == FC_SUCCESS, "failed to close after testing");
}
END_TEST

#endif // HAVE_EXODUS

// Test a dataset with temperature
// We are confident of the given values.
// Test metadata on everything, test big data only on selection:
// 1st of that thing, middle thing of that thing, last of that thing..
// some of the middle things nest: e.g. only test big var data on
// meshes that we testing big data on.
// NOTE: double numbers were generated with %.17g so that there were
// extra digits past the significant digits so that FC_DBL_EQUIV worked.
START_TEST(lsdyna_small)
{
  int i, j, k, m, n;
  int midID = 2;  
  int numSeq = 1, numGlbVar = 0, numGlbSeqVar = 5, numMesh = 1;
  int numSubset = 0, numVar = 2, numSeqVar = 15;
  int numWholeVar = 4, numVertVar = 3; //numElemVar = 8
  int temp_num, *temp_numStepPerSeqVar;
  FC_ReturnCode rc;
  FC_Dataset dataset;
  FC_Sequence* sequences;
  FC_Mesh* meshes;
  FC_Subset* subsets;
  FC_Variable *glbVars, **glbSeqVars, *vars, **seqVars;
  char* dataset_file_name = "../data/StressTests/Compression/d3plot";
  char dataset_name[1024] = "MAXSTRESSTEST-2";
  char* temp_name;
  char sequence_name[1024] = "time";
  int numStep = 12;
  FC_DataType dataType = FC_DT_FLOAT;
  float seqCoords[3] = { 0, 0.00019869821, 0.0010003619 };
  FC_DataType temp_dataType;
  float* temp_seqCoords;
  char meshNames[16][1024] = { "part_1" };
  int topodim = 3, numDim = 3, numVertPerElem;
  int numVertPerMesh[1] = { 216 };
  int numElemPerMesh[1] = { 125 };
  FC_ElementType elemType = FC_ET_HEX;
  // index by vert case, then dim (verts 0, 2, & 215)
  double coords[3][3] = { { 1, 3, -5 },
                          { -0.90268743, 3, -2.8146491 },
                          { 5, -5, 3 } };
  //index by elem case, then numVertPerElem (elems 0, 2, 124)
  int conns[3][8] = { { 0, 1, 2, 3, 4, 5, 6, 7 },
                      { 12, 13, 14, 15, 16, 17, 18, 19 },
                      { 9, 11, 10, 8, 28, 29, 45, 43 } };
  int temp_topodim, temp_numDim, temp_numVert, temp_numElem;
  FC_ElementType temp_elemType;
  int *temp_conns;
  double* temp_coords;
  char varNames[2][1024] = { "vert_IDs", "elem_IDs" };
  char glbSeqVarNames[5][1024] = { "global_KE",
                                   "global_IE",
                                   "global_TE",
                                   "global_velocity",
                                   "global_work" };
  int numCompPerGlbSeqVar[5] = { 1, 1, 1, 3, 1 };
  // index by numvar, then time case, then dim
  float global_data[5][3][3] = {
    { { 0 }, { 0.036911439 }, { 0.3661226 } },
    { { 9.9999997e-21 }, { 61.572666 }, { 1546.3472 } },
    { { 9.9999997e-21 }, { 61.609577 }, { 1546.7133 } },
    { { 0, 0, 0 }, 
      { 0.00045430628, -0.011901795, 0.00025582098 },
      { 0.00070099818, -0.79212922, 0.00049107615 } },
    { { 0 }, { 61.572666 }, { 1546.3472 } } };
  char seqVarNames[15][1024] = { "mat_IE",
				 "mat_KE",
				 "mat_velocity",
				 "mass",
				 "displacement",
				 "velocity",
				 "acceleration",
				 "sigma-xx",
				 "sigma-yy",
				 "sigma-zz",
				 "sigma-xy",
				 "sigma-yz",
				 "sigma-zx",
				 "plastic_strain",
				 "elem_death" };
  int numCompPerSeqVar[15] = { 1, 1, 3, 1, 3, 3, 3, 1, 1, 1,
                               1, 1, 1, 1, 1 }; 
  // index by numWholeVar, then by time case, dim
  float whole_data[4][3][3] = {
    { { 0 }, { 0.039879926 }, { 0.22684088 } },
    { { 0 }, { 0.00044252575 }, { 0.00075191067 } },
    { { 0, 0, 0.73499942 },
      { -0.0023735808, 0.00024508088, 0.73499942 },
      { -0.59749424, 0.0005285865, 0.73499942 } },
    { { 0 }, { 0.3513217 }, { 9.1221275 } } };
// index by mesh case, then numVertVar, then time case
// then entity case, and finally dim
float vert_data[3][3][3][3] = {
  //seqVar4
  { { { 0, 0, 0 },
      { 0, 0, 0 },
      { 0, 0, 0 } },
    { { 2.1696091e-05, -0.00049233437, -0.00011014938 },
      { -2.0265579e-05, -0.00048708916, -6.5326691e-05 },
      { 0, 0, 0 } },
    { { 0.00011134148, -0.0024626255, -0.00055932999 },
      { -0.00010222197, -0.002468586, -0.00033426285 },
      { 0, 0, 0 } } },
  //seqVar5
  { { { 0, 0, 0 },
      { 0, 0, 0 },
      { 0, 0, 0 } },
    { { 0.065397196, -0.17186408, -0.12185585 },
      { -0.061299477, 0.12308237, -0.030312598 },
      { 0, 0, 0 } },
    { { -0.01205952, -1.3162524, -0.29368019 },
      { -0.026158731, -1.0940681, -0.051936805 },
      { 0, 0, 0 } } },
  //seqVar6
  { { { 0, 0, 0 },
      { 0, 0, 0 },
      { 0, 0, 0 } },
    { { 3592.1484, -5.999723, -6183.8442 },
      { -1785.0477, -29903.186, -3625.8408 },
      { 0, 0, 0 } },
    { { -966.27118, -64484.707, -9764.917 },
      { -3286.9871, -73619.188, -483.51315 },
      { 0, 0, 0 } } } };
float elem_data[8][3][3] = {
  //seqVar7
  { { 0, 0, 0 },
    { 74.809372, 35.035095, -54.616531 },
    { 269.1846, 75.245438, -454.29208 } },
  //seqVar8
  { { 0, 0, 0 },
    { -1931.4797, -1932.2974, -2342.1147 },
    { -9670.7754, -9674.5312, -11896.861 } },
  //seqVar9
  { { 0, 0, 0 },
    { 35.434586, 74.95739, -53.808445 },
    { 76.483086, 270.08322, -455.07498 } },
  //seqVar10
  { { 0, 0, 0 },
    { -0.51165181, -19.970253, 0.63499516 },
    { -2.476459, -95.063782, 0.59348369 } },
  //seqVar11
  { { 0, 0, 0 },
    { 19.833843, 0.18871289, 0.29217809 },
    { 92.615768, 1.948981, -0.94667685 } },
  //seqVar12
  { { 0, 0, 0 },
    { 0.53380167, -0.74469697, -0.20847093 },
    { 1.6090398, -2.63253, -0.90011197 } },
  //seqVar13
  { { 0, 0, 0 },
    { 2.0240002e-05, 2.0251884e-05, 2.7233664e-05 },
    { 0.00010363303, 0.00010367855, 0.0001423199 } },
  //seqVar14
  { { 1, 1, 1 },
    { 1, 1, 1 },
    { 1, 1, 1 } } };
  int temp_numComp;
  FC_AssociationType temp_assoc;
  FC_MathType temp_mathType;
  float* temp_data;
  int* temp_int_data;

  // *** Test load dataset

  rc = fc_loadDataset(dataset_file_name, &dataset);
  fail_unless(rc == FC_SUCCESS, "failed to load dataset");

  // *** Test dataset

  rc = fc_getDatasetName(dataset, &temp_name);
  fail_unless(rc == FC_SUCCESS, "failed to get dataset name");
  fail_unless(!strcmp(temp_name, dataset_name), 
	      "mismatch of dataset name");
  free(temp_name);

  // *** Test sequences

  rc = fc_getSequences(dataset, &temp_num, &sequences);
  fail_unless(rc == FC_SUCCESS, "failed to get sequences");
  fail_unless(temp_num == numSeq, "mismatch of numSequence");
  for (i = 0; i < numSeq; i++) {

    // test sequence meta data
    rc = fc_getSequenceName(sequences[i], &temp_name);
    fail_unless(rc == FC_SUCCESS, "failed to get sequence name");
    fail_unless(!strcmp(temp_name, sequence_name), 
		"mismatch of sequence name");
    free(temp_name);
    rc = fc_getSequenceInfo(sequences[i], &temp_num, &temp_dataType);
    fail_unless(rc == FC_SUCCESS, "failed to get sequence info");
    fail_unless(temp_num == numStep, "mismatch of numStep");
    fail_unless(temp_dataType == dataType, 
		"mismatch of sequence dataType");

   // test sequence big data
    rc = fc_getSequenceCoordsPtr(sequences[i], (void**)&temp_seqCoords);
    fail_unless(rc == FC_SUCCESS, "failed to get sequence coords");
    for (j = 0; j < 3; j++) {
      int timeID = 0;
      if (j == 1)
	timeID = midID;
      else if (j == 2)
	timeID = numStep-1;
      fail_unless(FC_FLT_EQUIV(seqCoords[j], temp_seqCoords[timeID]),
		  "mismatch of seq coords");
    }
  }
  free(sequences);

  // *** Test global vars

  rc = fc_getGlobalVariables(dataset, &temp_num, &glbVars);
  fail_unless(rc == FC_SUCCESS, "failed to get global vars");
  fail_unless(temp_num == numGlbVar, "mismatch of numGlbVar");

  // *** Test global seq vars
  
  rc = fc_getGlobalSeqVariables(dataset, &temp_num, 
                                &temp_numStepPerSeqVar, &glbSeqVars);
  fail_unless(rc == FC_SUCCESS, "failed to get global seq vars");
  fail_unless(temp_num == numGlbSeqVar, "mismatch of numGlbSeqVar");
  for (i = 0; i < numGlbSeqVar; i++) {

    // test glb seq variable meta data
    rc = fc_getVariableName(glbSeqVars[i][0], &temp_name);
    fail_unless(rc == FC_SUCCESS, "failed to get glb seq variable name");
    fail_unless(!strcmp(temp_name, glbSeqVarNames[i]), 
                "mismatch of glb seq variable name");
    free(temp_name);
    fail_unless(temp_numStepPerSeqVar[i] == numStep, 
                "mismatch of numStepPerSeqVar");
    rc = fc_getVariableInfo(glbSeqVars[i][0], &temp_num, &temp_numComp,
                            &temp_assoc, &temp_mathType, &temp_dataType);
    fail_unless(rc == FC_SUCCESS, "failed to get glb seq variable info");
    fail_unless(temp_numComp == numCompPerGlbSeqVar[i], "mismatch of numComp");
    fail_unless(temp_num == 1, "mismatch of numDataPoint");
    fail_unless(temp_assoc == FC_AT_WHOLE_DATASET, "mismatch of assoc");
    if (temp_numComp == 3)
      fail_unless(temp_mathType == FC_MT_VECTOR, 
                  "mismatch of glb seq var mathType");
    else
      fail_unless(temp_mathType == FC_MT_SCALAR, 
                  "mismatch of glb seq var mathType");
    fail_unless(temp_dataType == dataType, "mismatch of glb seq var dataType");
    
    // test glb seq variable big data
    for (j = 0; j < 3; j++) {
      int timeID = 0;
      if (j == 1)
        timeID = midID;
      else if (j == 2)
        timeID = numStep-1;
      rc = fc_getVariableDataPtr(glbSeqVars[i][timeID], (void**)&temp_data);
      fail_unless(rc == FC_SUCCESS, "failed to get variable data");
      for (k = 0; k < temp_numComp; k++)
        fail_unless(FC_FLT_EQUIV(temp_data[k], global_data[i][j][k]),
                    "mismatch of global data");
    }
    
    // done with this glb seq variable
    fc_deleteSeqVariable(temp_numStepPerSeqVar[i], glbSeqVars[i]);
    free(glbSeqVars[i]);
  }
  free(temp_numStepPerSeqVar);
  free(glbSeqVars);

  // *** Test meshes

  rc = fc_getMeshes(dataset, &temp_num, &meshes);
  fail_unless(rc == FC_SUCCESS, "failed to get meshes");
  fail_unless(temp_num == numMesh, "mismatch of numMesh");
  fail_unless(numMesh == 1, "abort: following tests expect 1 mesh");
  for (i = 0; i < numMesh; i++) {

    // test mesh meta data
    rc = fc_getMeshName(meshes[i], &temp_name);
    fail_unless(rc == FC_SUCCESS, "failed to get mesh name");
    fail_unless(!strcmp(temp_name, meshNames[i]), "mismatch of mesh name");
    free(temp_name);
    
    rc = fc_getMeshInfo(meshes[i], &temp_topodim, &temp_numDim,
			&temp_numVert, &temp_numElem, &temp_elemType);
    fail_unless(rc == FC_SUCCESS, "failed to get mesh info");
    fail_unless(temp_topodim == topodim, "mismatch of topodim");
    fail_unless(temp_numDim == numDim, "mismatch of numDim");
    fail_unless(temp_numVert == numVertPerMesh[i], "mismatch of numVert");
    fail_unless(temp_numElem == numElemPerMesh[i], "mismatch of numElem");
    fail_unless(temp_elemType == elemType, "mismatch of elemType");
 
    // test mesh big data
    rc = fc_getMeshCoordsPtr(meshes[i], &temp_coords);
    fail_unless(rc == FC_SUCCESS, "failed to get mesh coords");
    rc = fc_getMeshElementConnsPtr(meshes[i], &temp_conns);
    fail_unless(rc == FC_SUCCESS, "failed to get mesh elem conns");
    numVertPerElem = fc_getElementTypeNumVertex(temp_elemType);
    for (j = 0; j < 3; j++) { // do 1st entity, middle entity & last entity
      int vertID = 0;
      int elemID = 0;
      if (j == 1) {
        vertID = midID;
        elemID = midID;
      }
      else if (j == 2) {
        vertID = numVertPerMesh[i]-1;
        elemID = numElemPerMesh[i]-1;
      }
      for (k = 0; k < numDim; k++)
	fail_unless(FC_FLT_EQUIV(temp_coords[vertID*numDim+k], 
				 coords[j][k]),
		    "mismatch of coords");
      fail_unless(!memcmp(&temp_conns[elemID*numVertPerElem],
                          conns[j], numVertPerElem*sizeof(int)),
		"mismatch of conns");
    }
   
    // *** Test subsets

    rc = fc_getSubsets(meshes[i], &temp_num, &subsets);
    fail_unless(rc == FC_SUCCESS, "failed to get subsets");
    fail_unless(temp_num == numSubset, "mismatch of numSubset");

    // *** Test variables

    rc = fc_getVariables(meshes[i], &temp_num, &vars);
    fail_unless(rc == FC_SUCCESS, "failed to get variables");
    fail_unless(temp_num == numVar, "mismatch of numVar");
    for (j = 0; j < numVar; j++) {
      
      // test var meta data
      rc = fc_getVariableName(vars[j], &temp_name);
      fail_unless(rc == FC_SUCCESS, "failed to get variable name");
      fail_unless(!strcmp(temp_name, varNames[j]),
		  "mismatch of variable name");
      free(temp_name);
      rc = fc_getVariableInfo(vars[j], &temp_num, &temp_numComp,
			      &temp_assoc, &temp_mathType, &temp_dataType);
      fail_unless(rc == FC_SUCCESS, "failed to get seq variable info");
      fail_unless(temp_numComp == 1, "mismatch of numComp");
      if (j == 0) {
	fail_unless(temp_num == numVertPerMesh[i], "mismatch of numDataPoint");
	fail_unless(temp_assoc == FC_AT_VERTEX, "mismatch of assoc");
      }
      else {
	fail_unless(temp_num == numElemPerMesh[i], "mismatch of numDataPoint");
	fail_unless(temp_assoc == FC_AT_ELEMENT, "mismatch of assoc");
      }
      fail_unless(temp_mathType == FC_MT_SCALAR, "mismatch of var mathType");
      fail_unless(temp_dataType == FC_DT_INT, "mismatch of var dataType");

      // test variable big data
      rc = fc_getVariableDataPtr(vars[j], (void**)&temp_int_data);
      fail_unless(rc == FC_SUCCESS, "failed to get variable data");
      for (k = 0; k < temp_num; k++) 
	fail_unless(temp_int_data[k] == k + 1, "data mismatch");
    }
    free(vars);

    // *** Test seq variables

    rc = fc_getSeqVariables(meshes[i], &temp_num, &temp_numStepPerSeqVar, 
			    &seqVars);
    fail_unless(rc == FC_SUCCESS, "failed to get seq variables");
    fail_unless(temp_num == numSeqVar, "mismatch of numSeqVar");
    for (j = 0; j < numSeqVar; j++) {
       
      // test seq variable meta data
      rc = fc_getVariableName(seqVars[j][0], &temp_name);
      fail_unless(rc == FC_SUCCESS, "failed to get seq variable name");
      fail_unless(!strcmp(temp_name, seqVarNames[j]), 
		  "mismatch of seq variable name");
      free(temp_name);
      fail_unless(temp_numStepPerSeqVar[j] == numStep, 
		  "mismatch of numStepPerSeqVar");
      rc = fc_getVariableInfo(seqVars[j][0], &temp_num, &temp_numComp,
			      &temp_assoc, &temp_mathType, &temp_dataType);
      fail_unless(rc == FC_SUCCESS, "failed to get seq variable info");;
      fail_unless(temp_numComp == numCompPerSeqVar[j], "mismatch of numComp");
      if (j < numWholeVar) {
	fail_unless(temp_num == 1, "mismatch of numDataPoint");
	fail_unless(temp_assoc == FC_AT_WHOLE_MESH, "mismatch of assoc");
      }
      else if (j < numWholeVar + numVertVar) {
	fail_unless(temp_num == numVertPerMesh[i], "mismatch of numDataPoint");
	fail_unless(temp_assoc == FC_AT_VERTEX, "mismatch of assoc");
      }
      else {
	fail_unless(temp_num == numElemPerMesh[i], "mismatch of numDataPoint");
	fail_unless(temp_assoc == FC_AT_ELEMENT, "mismatch of assoc");
      }
      if (temp_numComp == 3)
	fail_unless(temp_mathType == FC_MT_VECTOR, 
		    "mismatch of seq var mathType");
      else
	fail_unless(temp_mathType == FC_MT_SCALAR, 
		    "mismatch of seq var mathType");
      fail_unless(temp_dataType == dataType, "mismatch of seq var dataType");
      
      // test seq variable big data
      for (k = 0; k < 3; k++) {
	int varLU;
	int timeID = 0;
	if (k == 1)
	  timeID = midID;
	else if (k == 2)
	  timeID = numStep-1;
	rc = fc_getVariableDataPtr(seqVars[j][timeID], (void**)&temp_data);
	fail_unless(rc == FC_SUCCESS, "failed to get variable data");
	//printf("%d %d %d: temp = %g, known = %g\n", i, j, k,
	//     temp_data[0], whole_data[meshLU][j][k]);
	if (j < numWholeVar) {
          for (n = 0; n < temp_numComp; n++)
            fail_unless(FC_FLT_EQUIV(temp_data[n], whole_data[j][k][n]),
                        "mismatch of whole data");
	}
	else {
	  for (m = 0; m < 3; m++) {
	    int entityID = 0;
	    if (m == 1)
	      entityID = midID;
	    else if (m == 2)
	      entityID = temp_num-1;
	    if (j < numWholeVar + numVertVar) {
	      varLU = j - numWholeVar;
              for (n = 0; n < temp_numComp; n++)
                fail_unless(FC_FLT_EQUIV(temp_data[entityID*temp_numComp+n],
                                         vert_data[varLU][k][m][n]),
                            "mismatch of vert data");
	    }
            else { 
              varLU = j - numWholeVar - numVertVar;
              fail_unless(FC_FLT_EQUIV(temp_data[entityID],
                                       elem_data[varLU][k][m]),
                          "Mismatch elem data");
            }
          }
        }
      }
      
      // done with this seq variable
      fc_deleteSeqVariable(temp_numStepPerSeqVar[j], seqVars[j]);
      free(seqVars[j]);
    }
    free(temp_numStepPerSeqVar);
    free(seqVars);

    // done with this mesh
    fc_deleteMesh(meshes[i]);
  }
  free(meshes);

  // All done
  rc = fc_deleteDataset(dataset);
  fail_unless(rc == FC_SUCCESS, "failed to delete dataset at end of test");
}
END_TEST

// **** xml subset file IO

#ifdef HAVE_LIBXML2

START_TEST(aux_readwrite) {
  FC_ReturnCode rc;
  int i, j;
  xmlDocPtr doc;
  char* filename = "auxfiletest.fcx";
  FILE* file;
  int numType = 6;
  FC_AssociationType assocs[6] = { FC_AT_UNKNOWN, FC_AT_VERTEX,
				   FC_AT_EDGE, FC_AT_FACE,
				   FC_AT_ELEMENT, FC_AT_WHOLE_MESH };
  //int numMesh = 2;
  char *subsetNames[6] = { "apple", "banana nut", "c*w", 
			   "pig", "mug", "donut" }; 
  char *meshNames[2] = { "R.V.", "a cute tiny car" };
  int maxNumMembers[2] = { 25, 1001 };
  int numMembers[6] = { 2, 4, 1, 3, 0, 100 };
  int memberIDs[6][100] = { { 0, 1 }, { 3, 5, 7, 9 }, { 22 },
			    { 0, 1, 2 }, { }, { } };
  int temp_numSubset;
  char** temp_subsetNames, **temp_meshNames;
  FC_AssociationType* temp_assocs;
  int* temp_maxNumMembers, *temp_numMembers, **temp_memberIDs;
  int numTear = 3;
  char* tearNames[3] = { "Teensy BabyTear", "MommyTear", "DaddyTear" };
  int numSubsetPerTear[3] = { 1, 2, 3 };
  double lengths[3] = { .1, 33.3, 4.5 };
  int tearSubsetMap[3][3] = { { 2 }, { 1, 5 }, { 0, 3, 4 } };
  char** temp_tearNames, ***temp_subsetNamesPerTear;
  int temp_numTear, *temp_numSubsetPerTear;
  double* temp_lengths;

  // setup
  i = 5;
  for (j = 0; j < numMembers[i]; j++)
    memberIDs[i][j] = j;

  // test write
  file = fopen(filename, "w");
  fail_unless(file != NULL, "abort: failed to open file for writing");
  rc = _fc_writeAuxFileHeader(file);
  fail_unless(rc == FC_SUCCESS, "failed to write subset file header");
  for (i = 0; i < numType; i++) {
    if (i < numType/2)
      j = 0;
    else 
      j = 1;
    rc = _fc_writeAuxFileSubsetCore(file, subsetNames[i], meshNames[j],
				    assocs[i], maxNumMembers[j], 
				    numMembers[i], memberIDs[i]);
    fail_unless(rc == FC_SUCCESS, "failed to write subset");
  }
  for (i = 0; i < numTear; i++) {
    char* names[3];
    for (j = 0; j < numSubsetPerTear[i]; j++)
      names[j] = subsetNames[tearSubsetMap[i][j]];
    rc = _fc_writeAuxFileTear(file, tearNames[i], numSubsetPerTear[i],
			      names, lengths[i]);
    fail_unless(rc == FC_SUCCESS, "failed to write tear");
  }
  rc = _fc_writeAuxFileFooter(file);
  fail_unless(rc == FC_SUCCESS, "failed to write subset file footer");
  fclose(file);

  // test read
  rc = _fc_importAuxFileXMLDoc(filename, &doc);
  fail_unless(rc == FC_SUCCESS, "failed to import aux file");
  rc = _fc_getAuxFileSubsets(doc, &temp_numSubset, &temp_subsetNames,
		       &temp_meshNames, &temp_assocs, &temp_maxNumMembers,
		       &temp_numMembers, &temp_memberIDs); 
  fail_unless(rc == FC_SUCCESS, "failed to get subsets");
  fail_unless(temp_numSubset == numType, "mismatch of numSubset");
  for (i = 0; i < numType; i++) {
    if (i < numType/2)
      j = 0;
    else
      j = 1;
    fail_unless(!strcmp(temp_subsetNames[i], subsetNames[i]),
		"mismatch of name");
    fail_unless(!strcmp(temp_meshNames[i], meshNames[j]),
		"mismatch of mesh name");
    fail_unless(temp_assocs[i] == assocs[i], "mismatch of assoc");
    fail_unless(temp_maxNumMembers[i] == maxNumMembers[j],
		"mismatch of maxNumMember");
    fail_unless(temp_numMembers[i] == numMembers[i], "mismatch of numMember"); 
    fail_unless(!memcmp(temp_memberIDs[i], memberIDs[i],
			numMembers[i]*sizeof(int)), "mismatch of member IDs");
    free(temp_subsetNames[i]);
    free(temp_meshNames[i]);
    free(temp_memberIDs[i]);
  }
  free(temp_subsetNames);
  free(temp_meshNames);
  free(temp_assocs);
  free(temp_maxNumMembers);
  free(temp_numMembers);
  free(temp_memberIDs);
  rc = _fc_getAuxFileTears(doc, &temp_numTear, &temp_tearNames,
			   &temp_numSubsetPerTear, &temp_subsetNamesPerTear,
			   &temp_lengths);
  fail_unless(rc == FC_SUCCESS, "failed to get tears");
  fail_unless(temp_numTear == numTear, "mismatch of numTear");
  for (i = 0; i < numTear; i++) {
    fail_unless(!strcmp(temp_tearNames[i], tearNames[i]), "mismatch of name");
    fail_unless(temp_numSubsetPerTear[i] == numSubsetPerTear[i],
		"mismatch of numSubsetPerTear");
    for (j = 0; j < numSubsetPerTear[i]; j++) { 
      fail_unless(!strcmp(temp_subsetNamesPerTear[i][j],
			  subsetNames[tearSubsetMap[i][j]]), 
			  "mismatch of subsetNamePerTear");
      free(temp_subsetNamesPerTear[i][j]);
    }
    fail_unless(temp_lengths[i] == lengths[i], "mismatch of length");
    free(temp_tearNames[i]);
    free(temp_subsetNamesPerTear[i]);
  }
  free(temp_tearNames);
  free(temp_numSubsetPerTear);
  free(temp_subsetNamesPerTear);
  free(temp_lengths);
  xmlFreeDoc(doc);
}
END_TEST

START_TEST(bb_readwrite) {
  FC_ReturnCode rc;
  int i, j;
  char* filename = "fake.bb";
  FILE* file;
  int numBB = 4;
  char *names[4] = { "apple", "banana", "c*w_pig", "mug donut" };
  char *comments[4] = { "roses are red, violets are blue", "sugar",
			"is sweet", "and so is Splenda" };
  int stepIDs[4] = { -1, 1, 3, 2 };
  int numDims[4] = { 3, 3, 2, 1 };
  FC_Coords lowers[4] = { { 0, 0, 0 }, { .1, .2, 250 }, 
			  { -2, 1.6667, 250 }, { -2, -29, 29 } };
  FC_Coords uppers[4] = { { 1, 1, 1 }, { .2, .4, 1000 }, 
			  { -1, 2, 1000 }, { 0, 29, 29 } };
  int temp_numBB;
  char** temp_names, **temp_comments;
  int* temp_stepIDs;
  FC_Coords* temp_lowers, *temp_uppers;

  // test write
  file = fopen(filename, "w");
  fail_unless(file != NULL, "abort: failed to open file for writing");
  rc = _fc_writeBBFileHeader(file);
  fail_unless(rc == FC_SUCCESS, "failed to write bb file header");
  for (i = 0; i < numBB; i++) {
    rc = _fc_writeBBFileBoundingBox(file, names[i], stepIDs[i], comments[i],
				    numDims[i], lowers[i], uppers[i]);
    fail_unless(rc == FC_SUCCESS, "failed to write bb");
  }
  rc = _fc_writeBBFileFooter(file);
  fail_unless(rc == FC_SUCCESS, "failed to write bb file footer");
  fclose(file);

  // test read
  rc = _fc_readBBFile(filename, &temp_numBB, &temp_names, &temp_stepIDs,
		       &temp_comments, &temp_lowers, &temp_uppers); 
  fail_unless(rc == FC_SUCCESS, "failed to parse file");
  fail_unless(temp_numBB == numBB, "mismatch of numBB");
  for (i = 0; i < numBB; i++) {
    fail_unless(!strcmp(temp_names[i], names[i]),
		"mismatch of name");
    fail_unless(!strcmp(temp_comments[i], comments[i]),
		"mismatch of comments");
    fail_unless(temp_stepIDs[i] == stepIDs[i], "mismatch of stepID");
    for (j = 0; j < numDims[i]; j++) {
      fail_unless(temp_lowers[i][j] == lowers[i][j], "mismatch of lowers");
      fail_unless(temp_uppers[i][j] == uppers[i][j], "mismatch of uppers");
    }
    for (j = numDims[i]; j < 3; j++) {
      fail_unless(temp_lowers[i][j] == 0, "mismatch of lowers");
      fail_unless(temp_uppers[i][j] == 0, "mismatch of uppers");
    }
    free(temp_names[i]);
    free(temp_comments[i]);
  }
  free(temp_names);
  free(temp_comments);
  free(temp_stepIDs);
  free(temp_lowers);
  free(temp_uppers);
}
END_TEST

#endif // HAVE_LIBXML2

// *********************************************
// ***** Populate the Suite with the tests
// *********************************************

Suite *fileio_suite(void)
{
  Suite *suite = suite_create("FileIO");

  TCase *tc_bootstrap = tcase_create(" - File IO Bootstrap ");
  TCase *tc_bootstrap2 = tcase_create(" - File IO Bootstrap2 ");
  TCase *tc_helpers = tcase_create(" - File IO Helpers ");
  TCase *tc_overall = tcase_create(" - File IO Overall ");
#ifdef HAVE_EXODUS
  TCase *tc_exodus = tcase_create(" - Exodus IO ");
#endif
  TCase *tc_lsdyna = tcase_create(" - LSDyna IO ");
#ifdef HAVE_LIBXML2
  TCase *tc_xmlsubset = tcase_create(" - XML Subset IO ");
  TCase *tc_xmlbb = tcase_create(" - BoundingBox IO ");
#endif

  // File IO Bootstrap - read generated data files
  suite_add_tcase(suite, tc_bootstrap);
  tcase_add_test(tc_bootstrap, vert_files);
  tcase_add_test(tc_bootstrap, elem_files);
  tcase_add_test(tc_bootstrap, data_files);


  // File IO Bootstrap 2 - conveniences functions 
  suite_add_tcase(suite, tc_bootstrap2);
  tcase_add_checked_fixture(tc_bootstrap2, fileio_setup, fileio_teardown);
  tcase_add_test(tc_bootstrap2, make_mesh);
  tcase_add_test(tc_bootstrap2, make_var);

  // misc helpers
  suite_add_tcase(suite, tc_helpers);
  tcase_add_test(tc_helpers, namesfile);


  // all types of file io - (no fixtures)
  suite_add_tcase(suite, tc_overall);
  tcase_add_test(tc_overall, fileiotype);
  tcase_add_test(tc_overall, init_final);

#ifdef HAVE_EXODUS
  tcase_add_test(tc_overall, rewrite_test);
#endif

#ifdef HAVE_EXODUS
  // exodus file io
  suite_add_tcase(suite, tc_exodus);
  tcase_add_checked_fixture(tc_exodus, fileio_setup, fileio_teardown);
  tcase_add_test(tc_exodus, exo_private_bad_args);
  tcase_add_test(tc_exodus, exo_helpers);
  tcase_add_test(tc_exodus, exo_dataset);
  tcase_add_test(tc_exodus, exo_mesh);
  tcase_add_test(tc_exodus, exo_multi_mesh);
  tcase_add_test(tc_exodus, exo_node_set);
  tcase_add_test(tc_exodus, exo_seq_node_set);
  tcase_add_test(tc_exodus, exo_elem_side_set);
  tcase_add_test(tc_exodus, exo_seq_elem_side_set);
  tcase_add_test(tc_exodus, exo_global_subsets);
  tcase_add_test(tc_exodus, exo_whole_var);
  tcase_add_test(tc_exodus, exo_seq);
  tcase_add_test(tc_exodus, exo_var_cmp);
  tcase_add_test(tc_exodus, exo_global_seq_var);
  tcase_add_test(tc_exodus, exo_seq_var);
  tcase_add_test(tc_exodus, exo_var);
  tcase_add_test(tc_exodus, exo_global_ids);
#endif

  // lsdyna file io 
  suite_add_tcase(suite, tc_lsdyna);
  tcase_add_checked_fixture(tc_lsdyna, fileio_setup, fileio_teardown);
  tcase_add_test(tc_lsdyna, lsdyna_small);

#ifdef HAVE_LIBXML2  
  // xml subset files
  suite_add_tcase(suite, tc_xmlsubset);
  tcase_add_test(tc_xmlsubset, aux_readwrite);

  // xml bb files
  suite_add_tcase(suite, tc_xmlbb);
  tcase_add_test(tc_xmlbb, bb_readwrite); 
#endif

  return suite;
}
