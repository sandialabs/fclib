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
 * \file ex2_bigDataAccess.c
 * \brief Example of accessing big data.
 *
 * \description
 *    
 *    This is an annotated example that walks through the use of basic aspects
 *    of the Feature Characterization Library (FCLib).  It covers big data
 *    access (mesh coordinates, variable data values).
 *
 *    Please see ex1_intro.c first for an example of basic data access
 *    and using built-in FCLib characterizations. See also
 *    ex3_createCharactrztn.c for details bout creating your own
 *    characterization. 
 *
 *    NOTE: If you are looking at this file through the Doxygen interface,
 *    all objects and functions are linked to their documentation--so 
 *    follow the links for more details!
 *
 *    NOTE: Don't forget to try running the program while looking at this file.
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/examples/ex2_bigDataAccess.c,v $
 * $Revision: 1.6 $ 
 * $Date: 2006/09/26 21:31:32 $
 *
 * \modifications
 *    - 12/06/04 WSK, created.
 *    - 10/12/05 WSD, moved making a characterization to it's own file
 *      because gcc 4.0 didn't like function declarations within functions.
 */

#include <string.h>
#include <fc.h>

int main(void)
{
  FC_ReturnCode rc;
  int i, j;
  int numStep, topodim, dim, numVertex, numElement, numVertexPerElement;
  int numMember, maxNumMember, numDataPoint, numComponent, *memberIDs;
  int numMeshes;
  int numReturnSeqVars, *numStepPerReturnSeqVar, numReturnSubsets;
  int *conns;
  double *coords;
  void* seqCoords, *data;
  FC_Dataset dataset;
  FC_Sequence sequence;
  FC_Mesh mesh, *meshes;
  FC_Subset subset, *returnSubsets;
  FC_Variable *seqVariable, **returnSeqVars;
  FC_DataType dataType;
  FC_ElementType elemType;
  FC_AssociationType assocType;
  FC_MathType mathType;
  char dataset_name[1024];
  char mesh_name[1024];
  char seqVariable_name[1024];
  char sequence_name[1024];
  char subset_name[1024];
  char* temp_name;

  printf("Running ex2_bigDataAccess ...\n");
  printf("\n");

  // --- Init library and get Data Object handles

  // The code to get the data object handles is taken directly from
  // ex1_intro.c which we assume the reader has read.
  // Get the hex mesh, the temperature sequence variable, and the sequence
  rc = fc_setLibraryVerbosity(FC_WARNING_MESSAGES);
  if (rc != FC_SUCCESS)
    exit(rc);
  rc = fc_initLibrary();
  if (rc != FC_SUCCESS)
    exit(rc);
  strcpy(dataset_name, "../data/gen_multivar_seq.ex2");
  rc = fc_loadDataset(dataset_name, &dataset);
  if (rc != FC_SUCCESS)
    exit(rc);
  strcpy(mesh_name, "hex mesh");
  rc = fc_getMeshByName(dataset, mesh_name, &numMeshes, &meshes);
  if (rc != FC_SUCCESS)
    exit(rc);
  if (numMeshes != 1){
    rc = FC_INPUT_ERROR;
    fc_exitIfErrorPrintf(rc, "Failed to find mesh '%s'",mesh_name);
  }
  mesh = meshes[0];
  free(meshes);
  strcpy(subset_name, "every 5th element");
  rc = fc_getSubsetByName(mesh, subset_name, &numReturnSubsets, &returnSubsets);
  if (rc != FC_SUCCESS)
    exit(rc);
  if (numReturnSubsets != 1)
    exit(FC_INPUT_ERROR);
  subset = returnSubsets[0];
  free(returnSubsets);
  strcpy(seqVariable_name, "temperature per vertex");
  rc = fc_getSeqVariableByName(mesh, seqVariable_name, &numReturnSeqVars,
			       &numStepPerReturnSeqVar,&returnSeqVars);
  if (rc != FC_SUCCESS || numReturnSeqVars != 1)
    exit(rc);
  numStep = numStepPerReturnSeqVar[0];
  seqVariable = returnSeqVars[0];
  free(numStepPerReturnSeqVar);
  free(returnSeqVars);
  rc = fc_getSequenceFromSeqVariable(numStep, seqVariable, &sequence);
  if (rc != FC_SUCCESS)
    exit(rc);
  fc_getSequenceName(sequence, &temp_name);
  strcpy(sequence_name, temp_name);
  free(temp_name);

  // --- NOTE: differences between meta data and big data

  // For the major data objects, the data is separated into two categories:
  // meta data and big data. In general the big data tends to be big arrays
  // and the meta data provides the info for decoding the big data.
  // Meta data for data objects can be gotten using the fc_getXInfo() routines.
  // Meta data and most data returned from routines is copied from the library
  // and the caller is responsible for freeing any arrays. Big data is different
  // in the the library passes a pointer to its data and the caller is NOT
  // supposed to free it. Any routines that return a pointer to an array that
  // is NOT to be freed has 'Ptr' appended to it (e.g. fc_getMeshCoordsPtr()).
  // These distinctions will pointed out again in the example below.

  // Aside: Lazy loading of big data. If you load load from a file, all meta
  // data is loaded into memory, but the big data is not loaded until it is
  // actually needed. 

  // --- Accessing Sequence big data : sequence coords
  
  // A sequence is a 1D series of values that are the coordinates of the
  // sequence "space" (usually time). These can be stored as char, int,
  // float or double. The following code gets the meta data and uses
  // it to decode the big data. Note that seqCoords is declared as (void*)
  // and the 'Ptr' in the name of the big data call: fc_getSequenceCoordsPtr().
  printf("Sequence '%s':\n", sequence_name);
  fc_getSequenceInfo(sequence, &numStep, &dataType);
  printf("    Meta data: numStep = %d, dataType = %s\n", numStep,
         fc_getDataTypeText(dataType));
  printf("    Big data:\n");
  // Use the dataType to cast coords array to appropriate type
  fc_getSequenceCoordsPtr(sequence, &seqCoords);
  printf("        Coords:\n");
  if (dataType == FC_DT_CHAR) {
    for (i = 0; i < numStep; i++)
      printf("            %c\n", ((char*)seqCoords)[i]);
  }
  else if (dataType == FC_DT_INT) {
    for (i = 0; i < numStep; i++)
      printf("            %d\n", ((int*)seqCoords)[i]);
  }
  else if (dataType == FC_DT_FLOAT) {
    for (i = 0; i < numStep; i++)
      printf("            %f\n", ((float*)seqCoords)[i]);
  }
  else if (dataType == FC_DT_DOUBLE) {
    for (i = 0; i < numStep; i++)
      printf("            %f\n", ((double*)seqCoords)[i]);
  }
  else
    exit(-1);
  // Don't free seqCoords!! But do clear pointer if we're done with it.
  seqCoords = NULL;
  printf("\n");

  // --- Accessing Mesh big data : mesh coords & connectivity

  // The big data for a mesh are the spatial coordinate of the vertices
  // and the connectivity array which describes how the vertices make up
  // the elements.
  printf("Mesh '%s':\n", mesh_name);
  fc_getMeshInfo(mesh, &topodim, &dim, &numVertex, &numElement, &elemType);
  printf("    Meta data: topodim = %d, dim = %d\n", topodim, dim);
  printf("               numVertex = %d, numElement = %d, elemType = %s\n",
         numVertex, numElement, fc_getElementTypeText(elemType));
  printf("    Big data:\n");
  // The vertex coordinates are always doubles and are ordered as the coords
  // for the first vertex, then the coords for the next vertex, etc.
  // The number values that make up each set of coordinates is 'dim', the
  // dimensionality of the space the dataset occupies.
  fc_getMeshCoordsPtr(mesh, &coords);
  printf("        Vertex Coords (only for first 10 vertices):\n");
  for (i = 0; i < 10; i++) { // would normally loop over numVertex
    printf("            %d: ", i);
    for (j = 0; j < dim; j++)
      printf("%f ", coords[i*dim + j]); 
    printf("\n");
  }
  printf("               ...\n");
  // The connectivities are always ints (they are the IDs of the vertices that
  // make up the elements) and ordered as the vertices that make up the first
  // element, then the vertices that make up the second element, etc.
  // To decode, we also need to know the number of vertices in a element
  // which can be queried from the element type. 
  fc_getMeshElementConnsPtr(mesh, &conns);
  numVertexPerElement = fc_getElementTypeNumVertex(elemType);
  printf("        Elements Conns (only for first 10 elements):\n");
  for (i = 0; i < 10; i++) { // would normally loop over numElement
    printf("            %d: ", i);
    for (j = 0; j < numVertexPerElement; j++)
      printf("%d ", conns[i*numVertexPerElement + j]); 
    printf("\n");
  }
  printf("               ...\n");
  // Don't free coords or conns! But do clear pointer if done with them.
  coords = NULL;
  conns = NULL;
  printf("\n");

  // --- Accessing Variable big data : data values

  // A variable has a number of data points corresponding to the number of
  // entities it is associated with (e.g. a variable associated with the
  // vertices of a mesh must have numDataPoint = numVertex). The number of
  // components is how many values per data point make up the variable.
  // For example, a 3D vector variable has 3 components. The mathtype 
  // gives the user a little more information about how the components are
  // ordered and what they mean. For example, a 3D vector and a 2D
  // symtensor would both have numComponent = 3 (the symmetric tensor only
  // needs to store the upper part of the matrix) so the math type is
  // needed to properly use the data.
  //
  // The data is ordered as all data values for the first data point,
  // then all data values for the next data point, etc.
  printf("Sequence Variable '%s', Step #0:\n", seqVariable_name);
  fc_getVariableInfo(seqVariable[0], &numDataPoint, &numComponent, &assocType,
                     &mathType, &dataType);
  printf("    Meta data: numDataPoint = %d, numComponent = %d, assoc = %s\n",
         numDataPoint, numComponent, fc_getAssociationTypeText(assocType));
  printf("                mathType = %s, dataType = %s\n", 
         fc_getMathTypeText(mathType), fc_getDataTypeText(dataType));
  printf("    Big data:\n");
  // Use the dataType to cast data array to appropriate type
  // (Note, it'd be more efficient to test for data type outside of the loop,
  // but it makes the sample more cluttered).
  fc_getVariableDataPtr(seqVariable[0], &data);
  printf("        Data (only for first 10 data points):\n");
  for (i = 0; i < 10; i++) { // would normally loop over numDataPoint
    printf("            %d: ", i);
    for (j = 0; j < numComponent; j++) {
      if (dataType == FC_DT_CHAR)
        printf("%c ", ((char*)data)[i*numComponent + j]);
      else if (dataType == FC_DT_INT)
        printf("%d ", ((int*)data)[i*numComponent + j]);
      else if (dataType == FC_DT_FLOAT)
        printf("%f ", ((float*)data)[i*numComponent + j]);
      else if (dataType == FC_DT_DOUBLE)
        printf("%f ", ((double*)data)[i*numComponent + j]);
      else
        exit(-1);
    } 
    printf("\n");
  }
  printf("               ...\n");
  // Don't free data!! But do clear pointer if we're done with it.
  data = NULL;
  printf("\n");

  // --- Accessing Subset "big data" : member IDs

  // Accessing the "big data" for a subset is a little different than
  // for the other data types because subsets are a lot less static
  // than the other data types (meshes can't change the number of
  // vertices they have, but you can change how many vertices in a subset).
  // The "big data" for subsets is the IDs of the members of the
  // subset, but it is not treated as big data in that it is always
  // returned as a copy instead of a pointer. That is, you get a snapshot
  // of the state of a subset (if the subset changes, your snapshot is
  // out of date) and you are responsible for freeing the arrays.
  printf("Subset '%s':\n", subset_name);
  fc_getSubsetInfo(subset, &numMember, &maxNumMember, &assocType);
  printf("    Meta data: numMember = %d, maxNumMember = %d, assoc = %s\n",
         numMember, maxNumMember, fc_getAssociationTypeText(assocType)); 
  printf("    \"Big data\":\n");
  // Note: no 'Ptr' in this function name, so remember to free!
  fc_getSubsetMembersAsArray(subset, &numMember, &memberIDs);
  printf("        Member IDs (only the first 10):\n");
  for (i = 0; i < 10; i++)  // would normall loop over numMember
    printf("            %d: %d\n", i, memberIDs[i]);
  printf("               ...\n");
  // Special case: for subsets we do free the "big data"
  free(memberIDs);
  printf("\n");
  
  // Depending on what you want to do with the IDs after you get them, you
  // may also retrieve them as a linked list (helpful if growing a list) or
  // a mask (quicker for indexing).

  // --- All done

  // Memory cleanup - Arrays of small data and data objects not managed
  // by the library.
  free(seqVariable);

  // When you are all done it is a good idea to finalize the library.
  // This will automatically delete all datasets (and deleting a
  // dataset automatically deletes all of its sequences & meshes & etc.)
  fc_finalLibrary();

  // Bye!
  printf("... finished running ex2_bigDataAccess\n");
  exit(0);
}
