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
 * \file checkbase.c
 * \brief Unit testing of \ref DataTypes, \ref DataInterface, and
 *     \ref Handles Modules.
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/unittest/checkbase.c,v $
 * $Revision: 1.30 $ 
 * $Date: 2006/10/19 03:14:52 $
 *
 * \modifications
 *    - 02/25/04 WSK, created.
 *    - 05/07/04 WSK moved checking of mesh structures to checkmesh.c,
 *      moved checking of linked list and masks from checkutil to here.
 */

#include <stdlib.h>
#include <stdio.h>
#include <check.h>
#include <string.h>
#include "base.h"
#include "checkall.h"

// itests for all enumerated types

START_TEST(verbosityLevel) 
{
  int i;
  int numInput = 7;
  int inputs[7] = { FC_QUIET, FC_ERROR_MESSAGES, FC_WARNING_MESSAGES,
		    FC_LOG_MESSAGES, FC_DEBUG_MESSAGES, FC_QUIET - 1,
		    FC_DEBUG_MESSAGES + 1 };
  int validOuts[7] = { 1, 1, 1, 1, 1, 0, 0 };
  char* textOuts[7] = { "FC_QUIET", "FC_ERROR_MESSAGES", 
			"FC_WARNING_MESSAGES", "FC_LOG_MESSAGES",
			"FC_DEBUG_MESSAGES", 0, 0 };

  for (i = 0; i < numInput; i++) {
    fail_unless(fc_isVerbosityLevelValid(inputs[i]) == validOuts[i], 
		"isValid mismatch");
    if (textOuts[i])
      fail_unless(!strcmp(fc_getVerbosityLevelText(inputs[i]), textOuts[i]),
		  "getText mismatch");
    else
      fail_unless(fc_getVerbosityLevelText(inputs[i]) == NULL,
		  "getText mismatch");
  } 
}
END_TEST

START_TEST(returnCode) 
{
  int i;
  int numInput = 7;
  int inputs[7] = { FC_SUCCESS, FC_ERROR, FC_MEMORY_ERROR, FC_INPUT_ERROR, 
		    FC_FILE_IO_ERROR, FC_SUCCESS + 1, FC_FILE_IO_ERROR - 1 };
  int validOuts[7] = { 1, 1, 1, 1, 1, 0, 0 };
  char* textOuts[7] = { "FC_SUCCESS", "FC_ERROR", "FC_MEMORY_ERROR",
			"FC_INPUT_ERROR", "FC_FILE_IO_ERROR", 0, 0 };

  for (i = 0; i < numInput; i++) {
    fail_unless(fc_isReturnCodeValid(inputs[i]) == validOuts[i], 
		"isValid mismatch");
    if (textOuts[i])
      fail_unless(!strcmp(fc_getReturnCodeText(inputs[i]), textOuts[i]),
		  "getText mismatch");
    else
      fail_unless(fc_getReturnCodeText(inputs[i]) == NULL,
		  "getText mismatch");
  } 
}
END_TEST

START_TEST(elementType) 
{
  int i, j;
  int numInput = 12;
  int inputs[12] = { FC_ET_UNKNOWN, FC_ET_POINT, FC_ET_LINE, FC_ET_TRI, 
		     FC_ET_QUAD, FC_ET_TET, FC_ET_PYRAMID, FC_ET_PRISM, 
		     FC_ET_HEX, FC_ET_MIXED, FC_ET_UNKNOWN - 1,
		     FC_ET_MIXED + 1 };
  int validOuts[12] = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0 };
  char* textOuts[12] = { "FC_ET_UNKNOWN", "FC_ET_POINT", "FC_ET_LINE", "FC_ET_TRI", 
			 "FC_ET_QUAD", "FC_ET_TET", "FC_ET_PYRAMID", "FC_ET_PRISM", 
			 "FC_ET_HEX", "FC_ET_MIXED", 0, 0 };
  int topoDims[12] = { -1, 0, 1, 2, 2, 3, 3, 3, 3, -1, -1, -1 };
  int numVerts[12] = { -1, 1, 2, 3, 4, 4, 5, 6, 8, -1, -1, -1 };
  int numEdges[12] = { -1, 0, 1, 3, 4, 6, 8, 9, 12, -1, -1, -1 };
  int numFaces[12] = { -1, 0, 0, 1, 1, 4, 5, 5, 6, -1, -1, -1 };
  FC_ElementType faceTypes[12] = { FC_ET_UNKNOWN, FC_ET_UNKNOWN, FC_ET_UNKNOWN,
				   FC_ET_TRI, FC_ET_QUAD, FC_ET_TRI,
				   FC_ET_MIXED, FC_ET_MIXED, FC_ET_QUAD,
				   FC_ET_UNKNOWN, FC_ET_UNKNOWN };
  int numInput2 = 9;
  int inputs2[9] = { FC_AT_UNKNOWN, FC_AT_VERTEX, FC_AT_EDGE, FC_AT_FACE,
		     FC_AT_ELEMENT, FC_AT_WHOLE_MESH, FC_AT_WHOLE_DATASET,
		     FC_AT_UNKNOWN - 1, FC_AT_WHOLE_DATASET + 1 };
  // zero entries will get set for each i
  int numOuts2[9] = { -1, 0, 0, 0, 1, -1, -1, -1, -1 }; 
  
  for (i = 0; i < numInput; i++) {
    fail_unless(fc_isElementTypeValid(inputs[i]) == validOuts[i], 
		"isValid mismatch");
    if (textOuts[i]) 
      fail_unless(!strcmp(fc_getElementTypeText(inputs[i]), textOuts[i]),
		  "getText mismatch");
    else
      fail_unless(fc_getElementTypeText(inputs[i]) == NULL,
		  "getText mismatch");
    fail_unless(fc_getElementTypeTopoDim(inputs[i]) == topoDims[i],
		"getTopoDim mismatch");
    fail_unless(fc_getElementTypeNumVertex(inputs[i]) == numVerts[i],
		"getNumVertex mismatch");
    fail_unless(fc_getElementTypeNumEdge(inputs[i]) == numEdges[i],
		"getNumEdge mismatch");
    fail_unless(fc_getElementTypeNumFace(inputs[i]) == numFaces[i],
		"getNumFace mismatch");
    numOuts2[1] = numVerts[i];
    numOuts2[2] = numEdges[i];
    numOuts2[3] = numFaces[i];
    for (j = 0; j < numInput2; j++) {
      fail_unless(fc_getElementTypeNumEntity(inputs[i], inputs2[j]) ==
		numOuts2[j], "getNumEntity mismatch");
    }
    fail_unless(fc_getElementTypeFaceType(inputs[i]) == faceTypes[i],
		"getFaceType mismatch");
  } 
}
END_TEST

START_TEST(dataType) 
{
  int i;
  int numInput = 7;
  int inputs[7] = { FC_DT_UNKNOWN, FC_DT_CHAR, FC_DT_INT, FC_DT_FLOAT,
		    FC_DT_DOUBLE, FC_DT_UNKNOWN - 1, FC_DT_DOUBLE + 1 };
  int validOuts[7] = { 1, 1, 1, 1, 1, 0, 0 };
  char* textOuts[7] = { "FC_DT_UNKNOWN", "FC_DT_CHAR", "FC_DT_INT",
			"FC_DT_FLOAT", "FC_DT_DOUBLE", 0, 0 };
  int sizes[7] = { -1, sizeof(char), sizeof(int), sizeof(float),
		   sizeof(double), -1, -1 };
  
  for (i = 0; i < numInput; i++) {
    fail_unless(fc_isDataTypeValid(inputs[i]) == validOuts[i], 
		"isValid mismatch");
    if (textOuts[i])
      fail_unless(!strcmp(fc_getDataTypeText(inputs[i]), textOuts[i]),
		  "getText mismatch");
    else
      fail_unless(fc_getDataTypeText(inputs[i]) == NULL,
		  "getText mismatch");
    fail_unless(fc_sizeofDataType(inputs[i]) == sizes[i],
		"sizeOf mismatch");
  } 
}
END_TEST

START_TEST(mathType)
{
  int i;
  int numInput = 7;
  int inputs[7] = { FC_MT_UNKNOWN, FC_MT_SCALAR, FC_MT_VECTOR, FC_MT_SYMTENSOR,
		    FC_MT_TENSOR, FC_MT_UNKNOWN - 1, FC_MT_TENSOR + 1 };
  int validOuts[7] = { 1, 1, 1, 1, 1, 0, 0 };
  char* textOuts[7] = { "FC_MT_UNKNOWN", "FC_MT_SCALAR", "FC_MT_VECTOR",
			"FC_MT_SYMTENSOR", "FC_MT_TENSOR", 0, 0 };
  
  for (i = 0; i < numInput; i++) {
    fail_unless(fc_isMathTypeValid(inputs[i]) == validOuts[i], 
		"isValid mismatch");
    if (textOuts[i])
      fail_unless(!strcmp(fc_getMathTypeText(inputs[i]), textOuts[i]),
		  "getText mismatch");
    else
      fail_unless(fc_getMathTypeText(inputs[i]) == NULL,
		  "getText mismatch");
  } 
}
END_TEST

START_TEST(associationType) 
{
  int i;
  int numInput = 9;
  int inputs[9] = { FC_AT_UNKNOWN, FC_AT_VERTEX, FC_AT_EDGE, FC_AT_FACE,
		    FC_AT_ELEMENT, FC_AT_WHOLE_MESH, FC_AT_WHOLE_DATASET,
		    FC_AT_UNKNOWN - 1, FC_AT_WHOLE_DATASET + 1 };
  int validOuts[9] = { 1, 1, 1, 1, 1, 1, 1, 0, 0 };
  char* textOuts[9] = { "FC_AT_UNKNOWN", "FC_AT_VERTEX", "FC_AT_EDGE",
			"FC_AT_FACE", "FC_AT_ELEMENT", "FC_AT_WHOLE_MESH",
			"FC_AT_WHOLE_DATASET", 0, 0 };
  
  for (i = 0; i < numInput; i++) {
    fail_unless(fc_isAssociationTypeValid(inputs[i]) == validOuts[i], 
		"isValid mismatch");
    if (textOuts[i]) 
      fail_unless(!strcmp(fc_getAssociationTypeText(inputs[i]), textOuts[i]),
		  "getText mismatch");
    else
      fail_unless(fc_getAssociationTypeText(inputs[i]) == NULL,
		  "getText mismatch");
  } 
}
END_TEST

// Populate the Suite with the tests

Suite *base_suite(void)
{
  Suite *suite = suite_create("Base");

  TCase *tc_enums = tcase_create(" - Enumerated types ");
  
  suite_add_tcase(suite, tc_enums);
  tcase_add_test(tc_enums, verbosityLevel);
  tcase_add_test(tc_enums, returnCode);
  tcase_add_test(tc_enums, elementType);
  tcase_add_test(tc_enums, dataType);
  tcase_add_test(tc_enums, mathType);
  tcase_add_test(tc_enums, associationType);
  
  return suite;
}
