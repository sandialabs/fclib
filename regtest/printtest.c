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
 * \file printtest.c
 * \brief Test various fc_printX() routines.
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/regtest/printtest.c,v $
 * $Revision: 1.25 $ 
 * $Date: 2006/10/19 03:14:51 $
 *
 * \modifications
 *    - 12/8/04 W Koegler. Created
 */

#include <string.h> // for memcpy
#include <fc.h>
#include <fcP.h>

int main(void)
{
  FC_ReturnCode rc;
  int i, j;
  
  // 
  // initialize the fcdmf library
  //
  rc = fc_setLibraryVerbosity(FC_WARNING_MESSAGES);
  if (rc != FC_SUCCESS)
    fc_exitIfError(rc);
  rc = fc_initLibrary();
  if (rc != FC_SUCCESS)
    fc_exitIfError(rc);
  
  //
  // Print NULL handles (printing of real objects is handled
  // by running fcdump).
  //
  { 
    printf("Printing NULL handles:\n");
    fflush(NULL);
    fc_printDataset(FC_NULL_DATASET, "This is a NULL handle");
    fc_printSequence(FC_NULL_SEQUENCE, "This is a NULL handle", 1);
    fc_printMesh(FC_NULL_MESH, "This is a NULL handle", 1, 1, 1, 1);
    fc_printSubset(FC_NULL_SUBSET, "This is a NULL handle", 1);
    fc_printVariable(FC_NULL_VARIABLE, "This is a NULL handle", 1);
    printf("\n");
    fflush(NULL);
  }

  //
  // Print the internal data tables
  //
  { 
    FC_Dataset dataset;

    printf("Printing the internal tables (clean slate):\n");
    fflush(NULL);
    _fc_printDsTable("Should be empty");
    _fc_printSeqTable("Should be empty");
    _fc_printMeshTable("Should be empty");
    _fc_printSubTable("Should be empty");
    _fc_printVarTable("Should be empty");
    printf("\n");
    fflush(NULL);

    printf("Printing the internal tables (after loading dataset):\n");
    fflush(NULL);
    fc_loadDataset("../data/gen_multivar_seq.ex2", &dataset);
    _fc_printDsTable(NULL);
    _fc_printSeqTable(NULL);
    _fc_printMeshTable(NULL);
    _fc_printSubTable(NULL);
    _fc_printVarTable(NULL);
    printf("\n");
    fflush(NULL);

    fc_deleteDataset(dataset);

    printf("Printing the internal tables (after delete dataset):\n");
    fflush(NULL);
    _fc_printDsTable(NULL);
    _fc_printSeqTable(NULL);
    _fc_printMeshTable(NULL);
    _fc_printSubTable(NULL);
    _fc_printVarTable(NULL);
    printf("\n");
    fflush(NULL);
  }

  //
  // Additional subset printing routines
  //
  {
    FC_Dataset dataset;
    FC_Mesh mesh,*meshes;
    int numMeshes;
    FC_Subset subsets[10];
    FC_Variable variables[5][4];
    int numAssoc = 5, numDataType = 4, numEntity;
    FC_AssociationType assocs[5] = { FC_AT_VERTEX, FC_AT_EDGE, FC_AT_FACE,
				     FC_AT_ELEMENT, FC_AT_WHOLE_MESH };
    FC_DataType dataTypes[4] = { FC_DT_CHAR, FC_DT_INT, FC_DT_FLOAT,
				 FC_DT_DOUBLE };
    char* names[10] = { "empty verts", "empty edges", "empty faces", 
			"empty elements", "empty whole", "first 5 odd verts",
			"first 5 odd edges", "first 5 odd faces", 
			"first 5 odd elements", "the whole thing" };
    int numMember = 5, memberIDs[5] = { 1, 3, 5, 7, 9 };
    char charData[10] = { 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j' };
    int intData[10];
    float floatData[10];
    double doubleData[10];
    void* data, *goodDatas[4] = { charData, intData, floatData, doubleData };
    char temp_name[1024];

    // setup
    fc_loadDataset("../data/gen_prism.ex2", &dataset);
    fc_getMeshByName(dataset, "prism mesh", &numMeshes,&meshes);
    if (numMeshes != 1){
      exit(FC_INPUT_ERROR);
    }
    mesh= meshes[0];
    free(meshes);
    for (i = 0; i < 2*numAssoc; i++) 
      fc_createSubset(mesh, names[i], assocs[i%5], &subsets[i]);
    for (i = 0; i < numAssoc-1; i++) 
      fc_addArrayMembersToSubset(subsets[numAssoc+i], numMember, memberIDs);
    fc_addMemberToSubset(subsets[2*numAssoc-1], 0);
    for (i = 0; i < 2*numMember; i++) {
      intData[i] = i;
      floatData[i] = i + i/10.;
      doubleData[i] = i + i/1000.;
    }
    for (i = 0; i < numAssoc; i++) {
      fc_getMeshNumEntity(mesh, assocs[i], &numEntity);
      for (j = 0; j < numDataType; j++) {
	sprintf(temp_name, "%s_%s", fc_getAssociationTypeText(assocs[i]),
		fc_getDataTypeText(dataTypes[j]));
	fc_createVariable(mesh, temp_name, &variables[i][j]);
	data = calloc(numEntity, fc_sizeofDataType(dataTypes[j]));
	if (assocs[i] == FC_AT_WHOLE_MESH)
	  memcpy(data, goodDatas[j], fc_sizeofDataType(dataTypes[j]));
	else
	  memcpy(data, goodDatas[j],
		 2*numMember*fc_sizeofDataType(dataTypes[j]));
	fc_setVariableDataPtr(variables[i][j], numEntity, 1, assocs[i],
			      FC_MT_SCALAR, dataTypes[j], data);
      }
    }

    // do subset of mesh
    printf("Printing empty subsets of a mesh:\n");
    fflush(NULL);
    for (i = 0; i < numAssoc; i++)
      fc_printSubsetOfMesh(subsets[i], mesh);
    printf("\n");
    printf("Printing subsets of a mesh:\n");
    fflush(NULL);
    for (i = 0; i < numAssoc; i++)
      fc_printSubsetOfMesh(subsets[numAssoc + i], mesh);
    printf("\n");
    fflush(NULL);

    // do subset of variables
    printf("Printing empty subsets of variable:\n");
    fflush(NULL);
    for (i = 0; i < numAssoc; i++)
      fc_printSubsetOfVariable(subsets[i], variables[i][0]);
    printf("\n");
    printf("Printing subsets of all data types of variables:\n");
    fflush(NULL);
    for (i = 0; i < numAssoc; i++) 
      for (j = 0; j < numDataType; j++)
	fc_printSubsetOfVariable(subsets[numAssoc + i], variables[i][j]);
    printf("\n");
    fflush(NULL);

    // cleanup 
    fc_deleteDataset(dataset);
  }

  //
  // Print sorted int arrays
  //
  {
    FC_SortedIntArray sia;
    int numVal = 7, vals[7] = { -1, 0, 1, 4, 7, 9, 101 };

    // print empty - maxVal = 0
    printf("Printing a newly inited sorted int array:\n");
    fc_initSortedIntArray(&sia);
    fc_printSortedIntArray(&sia);
    printf("\n");
    fflush(NULL);

    // print full
    printf("Printing a populated sorted int array:\n");
    for (i = 0; i < numVal; i++)
      fc_addIntToSortedIntArray(&sia, vals[i]);
    fc_printSortedIntArray(&sia);
    printf("\n");
    fflush(NULL);

    // print empty - maxVal > 0
    printf("Printing an empty, but previously populated, sorted int array:\n");
    for (i = 0; i < numVal; i++)
      fc_deleteIntFromSortedIntArray(&sia, vals[i]);
    fc_printSortedIntArray(&sia);
    printf("\n");
    fflush(NULL);
  }

  //
  // Print & write features
  //
  {
    FC_Dataset dataset;
    FC_Mesh mesh, *meshes;
    int numMeshes, numReturnSeqVars, *numStepPerReturnSeqVar;
    FC_Variable **returnSeqVars;
    FC_Variable* seqVars = NULL;
    int numStep, numSegment;
    FC_Subset *segments = NULL, subset;
    char database_name[1024] = "../data/gen_gaussians.ex2";
    char mesh_name[1024] = "grid"; 
    char var_name[1024] = "temperature";
    char step_name[1024];
    double threshold = 30.;
    FC_FeatureGroup *group;
    char *group_file = "group.out";
    char *graph_file = "graph.out";
    
    fc_createFeatureGroup(&group);

    // test print/write on empty group
    printf("Printing & writing empty feature group:\n");
    fflush(NULL);
    fc_printFeatureGroup(group);
    fc_writeFeatureGroup(group, group_file);
   
    // Fill up the group - regions greater than threshold
    fc_loadDataset(database_name, &dataset);
    fc_getMeshByName(dataset, mesh_name, &numMeshes,&meshes);
    if (numMeshes != 1){
      exit(FC_INPUT_ERROR);
    }
    mesh = meshes[0];
    free(meshes);
    fc_getSeqVariableByName(mesh, var_name, &numReturnSeqVars,
			    &numStepPerReturnSeqVar,
			    &returnSeqVars);
    if (numReturnSeqVars != 1){
      exit(FC_INPUT_ERROR);
    }
    numStep = numStepPerReturnSeqVar[0];
    seqVars = returnSeqVars[0];
    free(returnSeqVars);
    free(numStepPerReturnSeqVar);

    for (i = 0; i < numStep; i++) {
      sprintf(step_name, "Step%d", i);
      fc_createThresholdSubset(seqVars[i], ">", threshold, step_name, 
			       &subset);
      fc_segment(subset, 0, &numSegment, &segments);
      fc_trackStep(i, numSegment, segments, group);
      free(segments);
    }
    free(seqVars);

    // test print/write on full group
    printf("Printing & writing feature group:\n");
    fflush(NULL);
    fc_printFeatureGroup(group);
    fc_writeFeatureGroup(group, group_file);

    // test print/write of feature graph (won't work on empty)
    printf("Printing & writing feature graph:\n");
    fflush(NULL);
    fc_printFeatureGraph(group);
    fc_writeFeatureGraph(group, graph_file);
    
    // cleanup
    fc_freeFeatureGroup(group);
    fc_deleteDataset(dataset);
  }

 
  // 
  // finalize the fcdmf library
  //
  fc_finalLibrary();
  
  exit(FC_SUCCESS);
}
