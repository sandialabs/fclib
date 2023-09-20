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
 * \file dataset_generator.c
 * \brief Generates datasets (for tests).
 * 
 * $Source: /usr/local/Repositories/fcdmf/fclib/data/dataset_generator.c,v $
 * $Revision: 1.9 $
 * $Date: 2006/11/07 23:49:01 $
 *
 * \description
 * 
 *   This is a utility to generate a variety of files for testing purposes.
 *   Each call to this program creates a single dataset file.
 *  
 *   It requires the name of the output file and two files per mesh
 *
 *     - 'vertex file' - lists the number of vertices & their 3D coordinates 
 *     - 'element file' - lists the number of vertices & the number of 
 *         elements, the type of element, and an ordered list of vertex indices
 *         (connectivities).
 *
 *   Optionally, multiple variables can be provided for each mesh via a file
 *   for each
 *  
 *     - 'variable file' - lists number of data points, whether they are 
 *         associated with vertices or elements, the name of the variable, 
 *         and then the variable values.
 *
 *   There are switches to create abitrary data or to copy variables
 *   to different types so that datasets are created with a variety of
 *   data for testing.
 *
 *   Usage - 
 *
 *   this option must come first...
 * 
 *     -m 2                 : number of meshes
 *
 *   required options...
 *   (# - indicates which mess associated with, numbering from 1)
 *   - -o output.output.ex2     : exodus output filename
 *   - -v # input.vert              : file containing vertex coordinates
 *   - -e # input.elem              : file containing element connectivities
 * 
 *   optional...
 *   - -d # input.data   : file containing data -- Can have multiple of these
 *   - -Md               : make char, int, float & double versions of variable,
 *                           requires use of -d flag
 *   - -Mi               : make vertex & element indices into int variables
 *   - -Ma               : make scalar, vector & symtensor versions of
 *                           variable, requires use of -d flag
 *   - -Mt               : make a time varying veresion of the variable 
 *                           (seqVar), requires use of -d flag
 *   - -Mq               : make sequences with different data types - char,
 *                           int, float & double, requires use of -d flag,
 *                           implies -Mt
 *   - -Ms               : make default scalar data (temperature per vertex)
 *   - -Mv               : make default vector data (velocity per vertex)
 *   - -Mst              : make default symmetric tensor data (stress per element)
 *   - -Msub             : make default subsets
 *
 *   Example -
 *
 *   dataset_generator -m 2 -o gook.ex2 -v 1 gook.vert -e 1 gook.elem \
 *   -v 2 goop.vert -e 2 goop.elem -Mi
 *   
 *
 *   Input file formats - 
 *
 *   Vertex File (*.vert) - 
 *   - int                    // the number of vertices
 *   - double double double   // the coordinates (x,y,z) of 1st vertex
 *   - double double double   // the coordinates (x,y,z) of 2nd vertex
 *   ...
 *   ...
 *   ...
 *   etc
 *
 *   Element File (*.elem) -
 *   - int int     //the number of vertices, the number of elements
 *   - string      // type: POINT, LINE, TRI, QUAD, TET, PYRAMID, PRISM, or HEX
 *   - string      // the mesh name
 *   - int int int (int ... ) // the indices of the verts making up 1st element
 *   - int int int (int ... ) // the indices of the verts making up 2nd element
 *   ...
 *   ...
 *   ...
 *   etc
 *   (Note: the index order is important, indices start from 0)
 *
 *   Data File (*.data) -
 *   - int         // the # of data points
 *   - string      // element association: choices are VERTEX or ELEMENT
 *   - string      // The type of data: choices are SCALAR, VECTOR or SYMTENSOR
 *   - string      // the variable name
 *   - double (double ...) // data values associated with 1st vert (or element)
 *   - double (double ...) // data values associated with 2nd vert (or element)
 *   ...
 *   ...
 *   ...
 *   etc
 *     
 *   simple_example.vert -
 *   - 4
 *   - 1.1 1.1 1.1
 *   - 1.1 1.2 1.2
 *   - 1.2 1.1 1.3
 *   - 1.2 1.2 1.4
 *   
 *   simple_example.elem -
 *   - 4 2
 *   - TRI
 *   - triangle mesh
 *   - 0 3 1
 *   - 2 3 0
 *   
 *   simple_example.data - 
 *   - 4 
 *   - VERTEX
 *   - SCALAR
 *   - temperature
 *   - 12.35
 *   - 13.6
 *   - 14.62
 *   - 25.12
 *   
 *   gives roughly -
 *   - (1.1,1.2) T=13.6
 *   - 1------3 (1.2,1.2) T=25.12
 *   - |...../|
 *   - |..../.|
 *   - |.../..|
 *   - |../...|
 *   - |./....|
 *   - |/.....|
 *   - 0------2 (1.2,1.1) T=14.62
 *   - (1.1,1.1)
 *   - T=12.35
 *
 * \modifications
 * - nsc 20sept01
 * - wsk 28mar02  made more general
 * - wsk 07may02  added multi mesh capability
 * - wsk 02Aug02  extended to elements of type POINT
 * - wsk 08Aug02  added mesh's name to the elem file spec
 * - wsk 03Dec02  Changed the basic structure of writing variables from a 
 *         complicated while loop to some straightforward tests. Also, many 
 *         things were componentized into subroutines.
 * - wsk 05Dec02  Adding new flags to create test saf files for comprehensive 
 *         testing of different types of variables
 * - wsk 16Dec02  Adding new flag to create test saf file containing time seris
 * - 05/27/04  WSK, added flag to create some subsets
 * - 08/22/04 WSD, removed SAF depedence--made it use fclib data structures.
 * - 02/06/06 WSD, the dataset names are now the file name stripped
 *       of 'gen_' and '.saf'.
 * - 01/10/06 ACG added flag for option for exodus data files rather than saf. 
 *            This will eventually be dropped when we drop saf.
 * - 10/04/07 ACG now defaults to exodus
 * - 10/15/2007 ACG removed saf
 */

#include <string.h>
#include <fc.h>
#include <fcP.h>

// Assumptions (FIX? extend file formats?)
// coords are doubles and are 3D (i.e x, y & z)
// data type for data files is assumed to be floats
// only one type of element in a mesh

#define MAX_STR_LEN 1024

int main(int argc, char **argv)
{
  FC_ReturnCode rc;
  int i, j;
  int m;                     // loop counter for meshes
  int v;                     // loop counter for data files
  int numStep = 11;          // Number of time steps to be generated
  int numMesh;               // number of meshes
  int numSeq = 1;            // number of sequences, if they exist
  FC_Dataset dataset;        // The dataset we are creating
  FC_Sequence *sequences;    // the sequences
  FC_Mesh mesh;              // a mesh
  char* datasetFileName = NULL;  // Name of the output dataset file
  char* datasetName;
  char **vertFiles;          // Names of files with vertex coords, numMesh
  char **elemFiles;          // Names files with element connectivities, numMesh
  int *numDataFiles;         // Number of data files per mesh
  char ***dataFiles;         // Pointer to array of data files
  int generate_scalar = 0;   // flag, 1 = generate default scalar variable
  int generate_vector = 0;   // flag, 1 = generate default vector variable
  int generate_tensor = 0;   // flag, 1 = generate default symtensor variable
  int generate_indices = 0;  // flag, 1 = generate vertex & element indices
  int generate_datatypes = 0;// flag, 1 = gen. all datatypes for each data file
  int generate_algtypes = 0; // flag, 1 = gen. all algtypes for each data file
  int generate_seqvar = 0;   // flag, 1 = gen. seq var from the data
  int generate_seqdatatypes = 0; // flag, 1 = gen. sequences with diff datatypes
  int generate_subsets = 0;  // flag, 1 = gen. subsets for each mesh

  //-------------------------------------------------------------------
  //  Process command line args
  //-------------------------------------------------------------------

  // If no arguments, print usage and exit
  if (argc == 1) {
  usage:
    printf("%s options...\n", argv[0]);
    printf("this option must come first:\n");
    printf("  -m 2               : number of meshes\n");
    printf("required options...\n");
    printf("  (# indicates which mesh associated with, numbered from 1)\n");
    printf("  -o output.exo     : exo output filename\n");
    printf("  -v # input.vert   : file containing vertex coordinates\n");
    printf("  -e # input.elem   : file containing element's vertex ids\n");
    printf("optional...\n");
    printf("  (-M options do not generate physically meaningful data, they\n");
    printf("   also operate globally - e.g. for all meshes or all variables)\n");
    printf("  -d # input.data   : file containing variable data (mutiple\n");
    printf("                        allowed each with their own -d flag)\n");
    printf("  -Mi               : make vertex & element indices into int variables (not possible for exodus)\n");
    printf("  -Md               : make char, int, float & double versions of\n");
    printf("                        variable, requires use of -d flag (char and some int not possible for exodus)\n");
    printf("  -Ma               : make scalar, vector & symtensor versions of\n");
    printf("                        variable, requires use of -d flag\n");
    printf("  -Mt               : make a time varying version of the data (seqVar)\n");
    printf("                        requires use of -d flag\n");
    printf("  -Mq               : make sequences with time as char, int, flaot\n");
    printf("                        and double, requires use of -d flag\n");
    printf("  -Ms               : make default scalar variable (temperature per vert)\n");
    printf("  -Mv               : make default vector variable (velocity per vert)\n");
    printf("  -Mst              : make default symmetric tensor variable (stress\n");
    printf("                        per element)\n");
    printf("  -Msub             : make default subsets\n");
    exit(0);
  }

  // Handle number of meshes
  numMesh = -1;
  if (!strcmp(argv[1], "-m")) 
    numMesh = atoi(argv[2]);
  if (numMesh < 1) goto usage;

  // setup filenames -- make space & initialize
  elemFiles = (char**)malloc(numMesh*sizeof(char*));
  vertFiles = (char**)malloc(numMesh*sizeof(char*));
  numDataFiles = (int*)malloc(numMesh*sizeof(int));
  dataFiles = (char***)malloc(numMesh*sizeof(char**));
  for (i = 0; i < numMesh; i++) {
    numDataFiles[i] = 0;
    dataFiles[i] = NULL;   // these are created from arguments
  }

  // loop over remaining arguments
  for (i = 3; i < argc; i++) {
    
    if ((!strcmp(argv[i], "help")) ||
        (!strcmp(argv[i], "-h")) ||
        (!strcmp(argv[i], "-help")) ||
        (!strcmp(argv[i], "-")) ||
        (!strcmp(argv[i], "--"))) 
      goto usage;
    
    else if (!strcmp(argv[i], "-o")) {
      if (i+1 >= argc || argv[i+1][0] == '-') goto usage;
      datasetFileName = argv[i+1];
      i++;
    }
  
    else if (!strcmp(argv[i], "-v")) {
      if (i+2 >= argc || argv[i+1][0] == '-') goto usage;
      m = atoi(argv[i+1]) - 1;
      vertFiles[m] = argv[i+2];
      i += 2;
    }
    
    else if (!strcmp(argv[i], "-e")) {
      if (i+2 >= argc || argv[i+1][0] == '-') goto usage;
      m = atoi(argv[i+1]) - 1;
      elemFiles[m] = argv[i+2];
      i += 2;
    }

    else if (!strcmp(argv[i], "-d")) {
      if (i+2 >= argc || argv[i+1][0] == '-') goto usage;
      m = atoi(argv[i+1]) - 1;
      dataFiles[m] = (char**)realloc(dataFiles[m], sizeof(char*) * 
				     (numDataFiles[m]+1));
      if (!dataFiles[m])
	fc_exitIfError(FC_MEMORY_ERROR);
      dataFiles[m][numDataFiles[m]] = argv[i+2];
      numDataFiles[m]++;
      i += 2;
    }
    
    else if (!strcmp(argv[i], "-Md")) 
      generate_datatypes = 1;
    
    else if (!strcmp(argv[i], "-Mi")) 
      generate_indices = 1;
    
    else if (!strcmp(argv[i], "-Ma")) 
      generate_algtypes = 1;
    
    else if (!strcmp(argv[i], "-Ms")) 
      generate_scalar = 1;
    
    else if (!strcmp(argv[i], "-Mt")) 
      generate_seqvar = 1;
    
    else if (!strcmp(argv[i], "-Mq")) {
      generate_seqvar = 1;
      generate_seqdatatypes = 1;
    }
    
    else if (!strcmp(argv[i], "-Mv")) 
      generate_vector = 1;
    
    else if (!strcmp(argv[i], "-Mst")) 
      generate_tensor = 1;
    
    else if (!strcmp(argv[i], "-Msub")) 
      generate_subsets = 1;
    
    else 
      fprintf(stderr, "%s: ignored argument `%s'\n", argv[0], argv[i]);
  }
  
  // check that required options were processed
  if (datasetFileName == NULL || !strcmp(datasetFileName, "")) goto usage;
  for (i = 0; i < numMesh; i++) {
    // a coordinate and a topo file are required for each mesh
    if (!strcmp(vertFiles[i],"") || !strcmp(elemFiles[i], "") ) goto usage;
    // cannot generate variables of all datatypes if no data is given
    if (generate_datatypes && numDataFiles[i] < 1) goto usage;
    // cannot generate variables of all algtypes if no data is given
    if (generate_algtypes && numDataFiles[i] < 1) goto usage;
    // cannot generate times if no data is given
    if (generate_seqvar && numDataFiles[i] < 1) goto usage;
  }

  // Create the dataset name, clip 'gen_' and '.ex2'
  datasetName = (char*)malloc((strlen(datasetFileName)+1)*sizeof(char));
  strcpy(datasetName, &datasetFileName[4]);
  datasetName[strlen(datasetName)-4] = '\0';

  //-------------------------------------------------------------------
  //  Initialize FCLib
  //-------------------------------------------------------------------
  // FIX? verbosity
  rc = fc_initLibrary();
  fc_exitIfErrorPrintf(rc, "failed to init fclib");

  //ACG
  rc = fc_setLibraryVerbosity(FC_WARNING_MESSAGES);
  

    
  //-------------------------------------------------------------------
  //  Create a new dataset
  //-------------------------------------------------------------------
  rc = fc_createDataset(datasetName, &dataset);
  fc_exitIfErrorPrintf(rc, "failed to create dataset '%s'", datasetFileName);

  //-------------------------------------------------------------------
  //  Create sequences
  //-------------------------------------------------------------------
  if (generate_seqvar) {
    float* float_times;
    void* times[4];
    char* seqNames[4];
    FC_DataType dataTypes[4];
    char* default_name = "time";
    char* type_names[4] = { "time (char)",
			    "time (int)",
			    "time (float)",
			    "time (double)" };
    
    // Make room for sequences handles (will need these later)
    if (generate_seqdatatypes)
      numSeq = 4;
    else
      numSeq = 1;
    sequences = (FC_Sequence*)malloc(numSeq*sizeof(FC_Sequence));

    // Make the float coords - will always use
    float_times = malloc(sizeof(float)*numStep);
    if (!float_times || !sequences)
      fc_exitIfError(FC_MEMORY_ERROR);
    for (i = 0; i < numStep; i++) 
      float_times[i] = i/10.; 

    // Set pointers (will get changed if generate_seqdatatypes)
    times[0] = float_times;
    seqNames[0] = default_name;
    dataTypes[0] = FC_DT_FLOAT;
    
    // If doing all datatypes, make multiple coords & set pointer
    if (generate_seqdatatypes) {
      char* char_times;
      int* int_times;
      double* double_times;
      FC_DataType types[4] = { FC_DT_CHAR, FC_DT_INT, FC_DT_FLOAT,
			       FC_DT_DOUBLE };

      // Make coords other than float
      char_times = malloc(numStep*sizeof(char));
      int_times = malloc(numStep*sizeof(int));
      double_times = malloc(numStep*sizeof(double));
      if (!char_times || !int_times || !double_times)
	fc_exitIfError(FC_MEMORY_ERROR);
      for (i = 0; i < numStep; i++) {
	char_times[i] = (char)(97+i);
	int_times[i] = i;
	double_times[i] = i + i/1000.;
      }
      times[0] = char_times;
      times[1] = int_times;
      times[2] = float_times;
      times[3] = double_times;
      for (i = 0; i < numSeq; i++) {
	seqNames[i] = type_names[i];
	dataTypes[i] = types[i];
      }
    }
       
    // Make the sequences
    for (i = 0; i < numSeq; i++) {
      rc = fc_createSequence(dataset, seqNames[i], &sequences[i]);
      fc_exitIfErrorPrintf(rc, "failed to create sequence");
      rc = fc_setSequenceCoordsPtr(sequences[i], numStep, dataTypes[i], 
				   times[i]);
      fc_exitIfErrorPrintf(rc, "failed to set sequence coords");
    }
  } // done with generate_seqvar
  

  //-------------------------------------------------------------------
  //  Process each mesh
  //-------------------------------------------------------------------
  for (m = 0; m < numMesh; m++) {
    
    //-------------------------------------------------------------------
    //  Create the mesh
    //-------------------------------------------------------------------
    
    rc = _fc_makeMeshFromFiles(dataset, vertFiles[m], elemFiles[m], &mesh);
    fc_exitIfErrorPrintf(rc, "failed to make mesh #%d\n", m);
    
    // Generate subsets if requested
    if (generate_subsets) {
      FC_Subset subset;
      int numElem;
      
      // the first 20 vertices
      rc = fc_createSubset(mesh, "first 20 vertices", FC_AT_VERTEX, &subset);
      fc_exitIfErrorPrintf(rc, "failed to create vertex subset");
      for (i = 0; i < 20; i++) {
	rc = fc_addMemberToSubset(subset, i);
	fc_exitIfErrorPrintf(rc, "failed to add member to vertex subset");
      }
      
      // every fifth element
      rc = fc_createSubset(mesh, "every 5th element", FC_AT_ELEMENT, &subset);
      fc_exitIfErrorPrintf(rc, "failed to create element subset");
      rc = fc_getMeshNumElement(mesh, &numElem);
      for (i = 0; i < numElem; i+=5) {
	rc = fc_addMemberToSubset(subset, i);
	fc_exitIfErrorPrintf(rc, "failed to add member to element subset");
      }
    } // end of generate_subsets
    
    
    //-------------------------------------------------------------------
    //  Create variable for each data file
    //-------------------------------------------------------------------    
    for (v = 0; v < numDataFiles[m]; v++) {
      FC_Variable variable;
      
      rc = _fc_makeVariableFromFile(mesh, dataFiles[m][v], &variable);
      fc_exitIfErrorPrintf(rc, "failed to create variable #%d on mesh #%d", 
			   v, m);    
      
      //-------------------------------------------------------------------
      // if the generate flags are used
      //-------------------------------------------------------------------
      
      if (generate_datatypes || generate_algtypes || generate_seqvar) {
	char baseName[1024], newName[1024];
	char* name;
	int numPoint, numComp;
	FC_AssociationType assoc;
	FC_MathType mathType;
	FC_DataType dataType;
	float* data;
	
	// Query the variable
	rc = fc_getVariableInfo(variable, &numPoint, &numComp, &assoc,
				&mathType, &dataType);
	fc_exitIfErrorPrintf(rc, "failed to get variable info");
	rc = fc_getVariableDataPtr(variable, (void**)&data);
	fc_exitIfErrorPrintf(rc, "failed to get variable data");

	// Rename the variable to give it a unique name
	fc_getVariableName(variable, &name);
	if (assoc == FC_AT_VERTEX)
	  sprintf(baseName, "%s per vertex", name);
	else
	  sprintf(baseName, "%s per element", name);
	free(name);
	if (generate_datatypes || generate_algtypes) {
	  sprintf(newName, "%s (original float scalar)", baseName);
	  rc = fc_changeVariableName(variable, newName);
	}
	else
	  rc = fc_changeVariableName(variable, baseName);
	fc_exitIfErrorPrintf(rc, "failed to change variable name");	

	// If generating all datatypes
	if (generate_datatypes) {
	  // Float already written - write char, int & double
	  FC_Variable temp_variable;
	  char* char_data;
	  int* int_data;
	  double* double_data;
	  void* datas[3];
	  FC_DataType dataTypes[3] = { FC_DT_CHAR, FC_DT_INT, FC_DT_DOUBLE };
	  char* name_endings[3] = { "(converted to chars)",
				    "(converted to ints)",
				    "(converted to doubles)" };
	  int temp_int;
	  
	  // Make datas
	  char_data = malloc(numPoint*numComp*sizeof(char));
	  int_data = malloc(numPoint*numComp*sizeof(int));
	  double_data = malloc(numPoint*numComp*sizeof(double));
	  if (!char_data || !int_data || !double_data) 
	    fc_exitIfError(FC_MEMORY_ERROR);
	  datas[0] = char_data;
	  datas[1] = int_data;
	  datas[2] = double_data;
	  // Convert (restric to basic character set - numerals)
	  for (i = 0; i < numPoint*numComp; i++) {
	    temp_int = (int)data[i];
	    if (temp_int < 33 || temp_int > 126  || // not a basic character 
		(temp_int >= 48 && temp_int <= 57) )  // also exclude numerals
	      char_data[i] = '&';
	    else
	      char_data[i] = data[i];
	    int_data[i] = temp_int;
	    double_data[i] = data[i];
	  }
	  
	  // Create the variables
	  for (i = 0; i < 3; i++) {
	    sprintf(newName, "%s %s", baseName, name_endings[i]);
	    rc = fc_createVariable(mesh, newName, &temp_variable);
	    fc_exitIfErrorPrintf(rc, "failed to created diff var datatype");
	    rc = fc_setVariableDataPtr(temp_variable, numPoint, numComp, assoc,
				       mathType, dataTypes[i], datas[i]);
	    fc_exitIfErrorPrintf(rc, "failed to set data for diff var datatype");
          }
	}
	
	// If generating all algtypes
	if (generate_algtypes) {
	  // Scalar already done - write vector and symtensor
	  FC_Variable temp_variable;
	  float* vector_data;
	  float* tensor_data;
	  int numComps[2] = { 3, 3 };
	  void* datas[2];
	  FC_MathType mathTypes[2] = { FC_MT_VECTOR, FC_MT_SYMTENSOR };
	  char* name_endings[2] = { "(expanded to a vector)",
				    "(expanded to a symtensor)" };
	  
	  // Make datas
	  vector_data = malloc(numPoint*numComps[0]*sizeof(float));
	  tensor_data = malloc(numPoint*numComps[1]*sizeof(float));
	  if (!vector_data || !tensor_data) 
	    fc_exitIfError(FC_MEMORY_ERROR);
	  datas[0] = vector_data;
	  datas[1] = tensor_data;
	  // Create extra components
	  for (i = 0; i < numPoint; i++) {
	    for (j = 0; j < numComps[0]; j++)
	      vector_data[i*numComps[0]+j] = data[i] + (double)j;
	    for (j = 0; j < numComps[1]; j++)
	      tensor_data[i*numComps[1]+j] = data[i]*(j+1.);
	  }
	  
	  // Create the variables
	  for (i = 0; i < 2; i++) {
	    sprintf(newName, "%s %s", baseName, name_endings[i]);
	    rc = fc_createVariable(mesh, newName, &temp_variable);
	    fc_exitIfErrorPrintf(rc, "failed to created diff var algtype");
	    rc = fc_setVariableDataPtr(temp_variable, numPoint, numComps[i], 
				       assoc, mathTypes[i], dataType, datas[i]);
	    fc_exitIfErrorPrintf(rc, "failed to set data for diff var algtype");
          }
	}
	
	// If generating any seq vars
	if (generate_seqvar) {
	  FC_Variable vars[numStep];
	  FC_Variable* seqVar;
	  
	  // Make the first (or only) seq var
	  vars[0] = variable;
	  for (i = 1; i < numStep; i++) {
	    float* temp_data = malloc(numPoint*numComp*sizeof(float));
	    if (!temp_data)
	      fc_exitIfError(FC_MEMORY_ERROR);
	    for (j = 0; j < numPoint*numComp; j++)
	      temp_data[j] = (data[j])*(1. + i/10.); // * 1.t
	    rc = fc_createVariable(mesh, baseName, &vars[i]);
	    fc_exitIfErrorPrintf(rc, "failed to create step of seq var");
	    rc = fc_setVariableDataPtr(vars[i], numPoint, numComp, assoc,
				       mathType, dataType, temp_data);
	    fc_exitIfErrorPrintf(rc, "failed to set data on step of seq var");
	  }
	  rc = fc_convertVariablesToSeqVariable(numStep, vars, sequences[0],
						NULL, &seqVar);
	  fc_exitIfErrorPrintf(rc, "failed to convert steps to seq var");
	  
	  // if generating sequences of different types - make copy onto
	  // other sequences
	  if (generate_seqdatatypes) {
	    char* name_endings[4] = { "(for char seq)" ,
				      "(for int seq)",
				      "(for float seq)", 
				      "(for double seq)" };

	    // rename 1st seq var
	    sprintf(newName, "%s %s", baseName, name_endings[0]);
	    rc = fc_changeSeqVariableName(numStep, seqVar, newName);
	    fc_exitIfErrorPrintf(rc, "failed to rename seqVar");

	    // Exact same data on other sequences
	    for (j = 1; j < numSeq; j++) {
	      FC_Variable* temp_seqVar;
	      sprintf(newName, "%s %s", baseName, name_endings[j]);
	      rc = fc_copySeqVariable(numStep, seqVar, mesh, sequences[j],
				      newName, &temp_seqVar);
	      fc_exitIfErrorPrintf(rc, "failed to copy seqVar");
	      free(temp_seqVar);      
	    }
	  } // end of generate_seqdatatypes
	} // end of generate_seqvar

      } // end of genaration flags built on var

    } // end of loop over data files

    //-------------------------------------------------------------------
    //  If requested, generate scalar per vertex variable
    //-------------------------------------------------------------------
    if (generate_scalar) {
      FC_Variable variable;
      int numVert;
      float* data; 
      
      // Get number of verts
      fc_getMeshNumVertex(mesh, &numVert);
      
      // random values between 62 & 82
      data = malloc(numVert*sizeof(float));
      for (i = 0; i < numVert; i++)
	data[i] = 72.0 + (10.0 * (((float)rand() / RAND_MAX) - 0.5));
      
      // Create variable
      rc = fc_createVariable(mesh, "temperature", &variable);
      fc_exitIfErrorPrintf(rc, "failed to create variable");
      rc = fc_setVariableDataPtr(variable, numVert, 1, FC_AT_VERTEX,
				 FC_MT_SCALAR, FC_DT_FLOAT, data);
      fc_exitIfErrorPrintf(rc, "failed to set data for generated scalar");
    }
    
    //-------------------------------------------------------------------
    //  If requested, generate vector data per vertex variable
    //-------------------------------------------------------------------
    if (generate_vector) {
      FC_Variable variable;
      int numVert, numComp = 3;
      float* data; 
      
      // Get number of verts
      fc_getMeshNumVertex(mesh, &numVert);
      
      // create vectors
      data = malloc(numVert*numComp*sizeof(float));
      for (i = 0; i < numVert; i++) {
	data[3*i]   = i*0.1; // x component
	data[3*i+1] = i*0.2; // y component
	data[3*i+2] = i*0.3; // z component
      }
      
      // Create variable
      rc = fc_createVariable(mesh, "velocity", &variable);
      fc_exitIfErrorPrintf(rc, "failed to create variable");
      rc = fc_setVariableDataPtr(variable, numVert, numComp, FC_AT_VERTEX,
				 FC_MT_VECTOR, FC_DT_FLOAT, data);
      fc_exitIfErrorPrintf(rc, "failed to set data for generated vector");
    }
    
    //-------------------------------------------------------------------
    // If requested, tensor per element variable
    //-------------------------------------------------------------------
    if (generate_tensor) {
      FC_Variable variable;
      int numElem, numComp = 3;
      float* data;
      FC_ElementType elemType;
      
      // Get number of elems
      fc_getMeshNumElement(mesh, &numElem);
      fc_getMeshElementType(mesh, &elemType);
      
      data = malloc(numElem*numComp*sizeof(float));
      for (i = 0; i < numElem; i++) {
	if (elemType == FC_ET_TRI) { // holdover from borrowing code from SAF
	  // stress for tris
	  data[3*i]   = i*0.1; // xx component
	  data[3*i+1] = i*0.2; // xy component
	  data[3*i+2] = i*0.3; // yy component
	}
	else {
	  // stress for everybody else
	  data[3*i]   = i+0.0; // xx component
	  data[3*i+1] = i+0.1; // xy component
	  data[3*i+2] = i+0.2; // yy component
	}
      }
      
      // Create variable
      rc = fc_createVariable(mesh, "stress tensor", &variable);
      fc_exitIfErrorPrintf(rc, "failed to create variable");
      rc = fc_setVariableDataPtr(variable, numElem, numComp, FC_AT_ELEMENT,
				 FC_MT_SYMTENSOR, FC_DT_FLOAT, data);
      fc_exitIfErrorPrintf(rc, "failed to set data for generated tensor");
    }
      
    //-------------------------------------------------------------------
    // If requested, generate vertex id & element id variables
    //-------------------------------------------------------------------
    if (generate_indices) {
      FC_Variable variable;
      int numVert, numElem;
      int* ids;
      
      // Get number of verts
      fc_getMeshNumVertex(mesh, &numVert);
      fc_getMeshNumElement(mesh, &numElem);
      
      // fill up ids
      ids = malloc((numVert > numElem ? numVert : numElem)*sizeof(int));
      for (i = 0; i < (numVert > numElem ? numVert : numElem); i++)
	ids[i] = i;
      
      // Create variables
      rc = fc_createVariable(mesh, "vertex ids", &variable);
      fc_exitIfErrorPrintf(rc, "failed to create variable");
      rc = fc_setVariableData(variable, numVert, 1, FC_AT_VERTEX,
			      FC_MT_SCALAR, FC_DT_INT, ids);
      fc_exitIfErrorPrintf(rc, "failed to set data for vertex ids");
      rc = fc_createVariable(mesh, "element ids", &variable);
      fc_exitIfErrorPrintf(rc, "failed to create variable");
      rc = fc_setVariableData(variable, numElem, 1, FC_AT_ELEMENT,
			      FC_MT_SCALAR, FC_DT_INT, ids);
      fc_exitIfErrorPrintf(rc, "failed to set data element ids");
      free(ids);
    }
    
  } // end of loop over meshes
    
  //-------------------------------------------------------------------
  // Write dataset
  //-------------------------------------------------------------------

  rc = fc_writeDataset(dataset, datasetFileName, FC_FT_EXODUS);

  fc_exitIfErrorPrintf(rc, "failed to write dataset");
  
  //-------------------------------------------------------------------
  // Finalize access to the library
  //-------------------------------------------------------------------
  rc = fc_finalLibrary();
  fc_exitIfErrorPrintf(rc, "failed to final library");

  // clean up

  exit(0);
}
