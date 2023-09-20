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
 * \file varSmooth.c
 * \brief Create a new dataset in which the desired variable has
 *        been smoothed.
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/tools/varSmooth.c,v $
 * $Revision: 1.24 $ 
 * $Date: 2006/11/17 22:51:25 $
 *
 * \description
 *    Usage: varSmooth [options] dataset "variable name" radius new_file
 *
 *    Create a new dataset in which the desired variable has been
 *    smoothed using the given kernel radius. The new dataset will
 *    be of the same file type as the original dataset.
 *
 *    Example: varSmooth data.ex2 "temperature" 1.2 data_smooth.ex2
 *
 *    Note that smoothing a non-time-varying variable but applying a
 *    time-varying displacment will create a time-varying smooted variable.
 *
 * \todo  Add switch so user can specify output file type. 
 *
 * \modifications
 *   - 10/11/04 WSK, created. 
 *   - 05/08/06 WSD, added displaced option (ended up greatly reworking
 *       internals).
 */

#include <string.h>
#include "fc.h"
 
int main(int argc, char** argv) {
  FC_ReturnCode rc;
  int i, j;
  int numMesh, numBasicVar, numSeqVar, numStep = 0, temp_numStep;
  FC_Dataset ds_orig;
  FC_Mesh* meshes;
  FC_Variable* basic_vars = NULL, *displs = NULL;  
  FC_Variable** seq_vars = NULL, **seqDispls = NULL;   
  FC_Sequence sequence;
  char* file_name = NULL, *new_file_name;
  char* var_name;
  char* displ_name = NULL;
  FC_VerbosityLevel verbose_level = FC_QUIET; 
  double radius;
  int doLocalSmooth = 0;
  FC_Variable* smooth_vars = NULL, **smooth_seq_vars = NULL;
  int numSpecialSeqVar = 0; // smoothing basic var w/ seq displ ...
  FC_Variable **special_seq_vars = NULL;
  int doReplace = 0;
  int doDispl = 0; // value of 1 = basic var, 2 = seq var
  FC_FileIOType fileType;
  int *basic_meshIDs, *seq_meshIDs, *meshUsedFlags;
  
  // handle arguments
  if (argc < 5) {
  usage:
    printf("usage: %s [options] dataset \"var name\" radius new_file_name\n", 
           argv[0]);
    printf("options: \n");
    printf("   -d var_name     : calculate on the displaced dataset using this "
           "displacement\n");
    printf("                     variable\n");
    printf("   --replace       : replace the original variable with the smoothed "
           "variable\n");
    printf("   --doLocalSmooth : if true, only points on local mesh are used for smoothing\n"); 
    printf("                     if false, points on all meshes are used for smoothing \n");
    printf("   -h              : print this help message\n");
    printf("   -v              : verbose: print warning and error messages\n");
    printf("   -V              : very verbose: print log and error messages\n");
    printf("\n");
    printf("Create a new dataset which will be a copy of the original dataset\n");
    printf("plus an additional variable which is the geometrically smoothed\n");
    printf("version of the specified variable. The smoothing will be done using\n"); 
    printf("the given kernel radius and the new variable's name will be the \n");
    printf("original variable's name postfixed with '_smoothed'. (If the\n");
    printf("'--replace' flag is used, the smoothed variable will replace the\n"); 
    printf("original variable and keep it's name.)\n");
    fflush(NULL);
    exit(-1);
  }
  for (i = 1; i < argc; i++) {
    if (!strcmp(argv[i], "-d")) {
      doDispl = 1;
      i++;
      displ_name = argv[i];
    }
    else if (!strcmp(argv[i], "--replace") || !strcmp(argv[i], "-replace")) {
      doReplace = 1;
    }
    else if (!strcmp(argv[i], "--doLocalSmooth") || !strcmp(argv[i], "-doLocalSmooth")) {
      doLocalSmooth = 1;
    }
    else if (!strcmp(argv[i], "-v")) {
      verbose_level = FC_WARNING_MESSAGES;
    }
    else if (!strcmp(argv[i], "-V")) {
      verbose_level = FC_LOG_MESSAGES;
    }
    else if (!strncmp(argv[i], "-h", 2) || !strncmp(argv[i], "-", 1))
      goto usage;
    else {
      if (i+3 >= argc)
        goto usage;
      file_name = argv[i];
      var_name = argv[i+1];
      radius = atof(argv[i+2]);
      new_file_name = argv[i+3];
      i+=3;
    }
  }
  if (!file_name || !var_name || radius <= 0. || !new_file_name)
    goto usage;
  
  // init library and load dataset & meshes
  rc = fc_setLibraryVerbosity(verbose_level);
  fc_exitIfError(rc);
  rc = fc_initLibrary();
  fc_exitIfError(rc);
  rc = fc_loadDataset(file_name, &ds_orig);
  fc_exitIfError(rc);
  rc = fc_getMeshes(ds_orig, &numMesh, &meshes);
  fc_exitIfError(rc);
   
  // loop over meshes, looking for var and or seq var to be smoothed
  basic_meshIDs = (int*)malloc(numMesh*sizeof(int));
  seq_meshIDs = (int*)malloc(numMesh*sizeof(int));
  meshUsedFlags = (int*)calloc(numMesh, sizeof(int)); // init to zeros
  basic_vars = (FC_Variable*)malloc(numMesh*sizeof(FC_Variable));
  seq_vars = (FC_Variable**)malloc(numMesh*sizeof(FC_Variable*));
  if (!basic_meshIDs || !seq_meshIDs || !meshUsedFlags || !basic_vars || !seq_vars)
    fc_exitIfError(FC_MEMORY_ERROR);
  numBasicVar = 0;
  numSeqVar = 0;
  for (i = 0; i < numMesh; i++) {
    // try to access as basic var & seq var
    rc = fc_getOrGenerateUniqueVariableByName(meshes[i], var_name,
					      &basic_vars[numBasicVar]);
    fc_exitIfError(rc);
    if (FC_HANDLE_EQUIV(basic_vars[numBasicVar], FC_NULL_VARIABLE)){
      fc_printfWarningMessage("Problems finding unique variable '%s'",
			      var_name);
    } else {
      basic_meshIDs[numBasicVar] = i;
      meshUsedFlags[i] = 1;
      numBasicVar++;
    }

    rc = fc_getOrGenerateUniqueSeqVariableByName(meshes[i], var_name,
						 &temp_numStep,
						 &seq_vars[numSeqVar]);
    fc_exitIfError(rc);
    if (!seq_vars[numSeqVar]){
      fc_printfWarningMessage("Problems finding unique variable '%s'",
			      var_name);
    } else {
      if (numSeqVar == 0) {
	numStep = temp_numStep;
	rc = fc_getSequenceFromSeqVariable(numStep, seq_vars[numSeqVar],
					   &sequence);
	fc_exitIfErrorPrintf(rc, "Could not get squence from seq variable");
      }
      else if (temp_numStep != numStep)
	fc_exitIfErrorPrintf(FC_ERROR, "Expected the sequence variable to have "
			     "the same number of steps on all meshes");
      seq_meshIDs[numSeqVar] = i;
      meshUsedFlags[i] = 1;
      numSeqVar++;
    }
  }
  if (numBasicVar < 1 && numSeqVar < 1)
    fc_exitIfErrorPrintf(FC_ERROR, "Failed to find variable '%s' in dataset "
                         "'%s'", var_name, file_name)

  // if displacing, make sure we have displacement variables
  if (doDispl) {
    int hasBasicDispl, hasSeqDispl;
    displs = (FC_Variable*)malloc(numMesh*sizeof(FC_Variable));
    seqDispls = (FC_Variable**)malloc(numMesh*sizeof(FC_Variable*));
    if (!displs || ! seqDispls) 
      fc_exitIfError(FC_MEMORY_ERROR);
    // not catching errors on purpose
    for (i = 0; i < numMesh; i++) {
      rc = fc_getOrGenerateUniqueVariableByName(meshes[i], displ_name,
						&displs[i]);
      fc_exitIfError(rc);

      rc = fc_getOrGenerateUniqueSeqVariableByName(meshes[i], displ_name,
						   &temp_numStep,
						   &seqDispls[i]);
      fc_exitIfError(rc);


      if (seqDispls[i]) {
	if (numStep == 0) {
	  numStep = temp_numStep;
	  rc = fc_getSequenceFromSeqVariable(numStep, seqDispls[i],
					     &sequence);
	  fc_exitIfErrorPrintf(rc, "Could not get sequence from seq variable");
	}
	else if (temp_numStep != numStep)
	  fc_exitIfErrorPrintf(FC_ERROR, "Expected all seq variables to "
			       "have same number of steps");
      }
    }
    // check that we have displacement vars
    hasBasicDispl = 0;
    for (i = 0; i < numMesh; i++) {
      if (meshUsedFlags[i]) {
	if (fc_isValidDisplacementVariable(meshes[i], displs[i]))
	  hasBasicDispl = 1;
	else {
	  hasBasicDispl = 0;
	  break;
	}
      }
    }
    hasSeqDispl = 0;
    for (i = 0; i < numMesh; i++) {
      if (meshUsedFlags[i]) {
	if (seqDispls[i] &&
	    fc_isValidDisplacementVariable(meshes[i], seqDispls[i][0])){
	  hasSeqDispl = 1;
	} else {
	  hasSeqDispl = 0;
	  break;
	}
      }
    }
    if (!hasBasicDispl && !hasSeqDispl) 
      fc_exitIfErrorPrintf(FC_ERROR, "Could not find valid displacment "
			   "variable '%s'", displ_name);
    if (hasSeqDispl)
      doDispl = 2;
    if (hasBasicDispl && hasSeqDispl)
      fc_printfWarningMessage("Displacement variable '%s' exists as both a "
			      "basic var and a seq var (will ignore basic)",
			      displ_name);
  }

  // do the smoothing - 6 total cases: 3 displ cases (none, basic, seq) X
  // 2 var cases (basic, seq). The first 4 are handle in the first if.
  if (doDispl < 2) {
    FC_VertexBin *vertexBins[numMesh];
    FC_Mesh temp_meshes[numMesh];
    FC_VertexBin *temp_bins[numMesh];
    FC_VertexBin **temp_seqBins[numMesh];
    FC_Variable temp_displs[numMesh];
    FC_Variable *temp_seqDispls[numMesh];

    // Make the bins
    for (i = 0; i < numMesh; i++) {
      vertexBins[i] = NULL;
      if (meshUsedFlags[i]) {
	if (doDispl)
	  rc = fc_createDisplacedMeshVertexBin(meshes[i], displs[i],
					      &vertexBins[i]);
	else
	  rc = fc_createMeshVertexBin(meshes[i], &vertexBins[i]);
	fc_exitIfErrorPrintf(rc, "Failed to create mesh vertex bin");
      }
    }
    
    // first do basic vars
    if (numBasicVar > 0) {
      // create smoothed variable
      for (i = 0; i < numBasicVar; i++) {
	temp_meshes[i] = meshes[basic_meshIDs[i]];
	temp_bins[i] = vertexBins[basic_meshIDs[i]];
	if (doDispl)
	  temp_displs[i] = displs[basic_meshIDs[i]];
      }
      if (doDispl)
	rc = fc_displacedGeomSmoothVariable(numBasicVar, temp_meshes,
					    temp_displs, temp_bins, basic_vars,
					    radius, doLocalSmooth, &smooth_vars);
      else
	rc = fc_geomSmoothVariable(numBasicVar, temp_meshes, temp_bins, 
				 basic_vars, radius, doLocalSmooth, &smooth_vars);
      fc_exitIfErrorPrintf(rc, "Failed to smooth variable");
    }
    
    // then do the seq vars by step
    if (numSeqVar > 0) {
      // created smoothed sequence variable
      for (i = 0; i < numSeqVar; i++) {
	temp_meshes[i] = meshes[seq_meshIDs[i]];
	temp_bins[i] = vertexBins[seq_meshIDs[i]];
	temp_seqBins[i] = malloc(numStep*sizeof(FC_VertexBin*));
	for (j = 0; j < numStep; j++)
	  temp_seqBins[i][j] = vertexBins[seq_meshIDs[i]];
	if (doDispl) {
	  temp_displs[i] = displs[seq_meshIDs[i]];
	  temp_seqDispls[i] = seqDispls[seq_meshIDs[i]];
	}
      }
      if (doDispl)
	rc = fc_displacedGeomSmoothSeqVariable(numSeqVar, temp_meshes, numStep,
					       temp_seqDispls, temp_seqBins,
					       numStep, seq_vars, radius,
					       doLocalSmooth, &smooth_seq_vars);
      else
	rc = fc_geomSmoothSeqVariable(numSeqVar, temp_meshes,
				      temp_bins, numStep, seq_vars,
				      radius, doLocalSmooth, &smooth_seq_vars);
      fc_exitIfErrorPrintf(rc, "Failed to smooth sequence variable");
      for (i = 0; i < numSeqVar; i++)
	free(temp_seqBins[i]);
    }
    
    // cleanup
    for (i = 0; i < numMesh; i++) {
      if (vertexBins[i])
	fc_freeVertexBin(vertexBins[i]);
    }
  }
  else if (doDispl == 2) { // displacements are seq variable
    FC_VertexBin *vertexBins[numMesh];
    FC_Mesh temp_meshes[numMesh];
    FC_VertexBin *temp_bins[numMesh];
    FC_Variable temp_seq_vars[numMesh];
    FC_Variable temp_displs[numMesh];
    FC_Variable *temp_smooth_vars[numStep];
    FC_Variable *temp_smooth_seq_vars[numStep];
    for (j = 0; j < numStep; j++) {
      temp_smooth_vars[j] = NULL;
      temp_smooth_seq_vars[j] = NULL;
    }

    for (j = 0; j < numStep; j++) {

      // Make the bins
      for (i = 0; i < numMesh; i++) {
	vertexBins[i] = NULL;
	if (meshUsedFlags[i]) {
	  rc = fc_createDisplacedMeshVertexBin(meshes[i], seqDispls[i][j],
					       &vertexBins[i]);
	  fc_exitIfErrorPrintf(rc, "Failed to create mesh vertex bin");
	}
      }
    
      // first do basic vars
      if (numBasicVar > 0) {
	// create smoothed variable
	for (i = 0; i < numBasicVar; i++) {
	  temp_meshes[i] = meshes[basic_meshIDs[i]];
	  temp_bins[i] = vertexBins[basic_meshIDs[i]];
	  temp_displs[i] = seqDispls[basic_meshIDs[i]][j];
	}
	rc = fc_displacedGeomSmoothVariable(numBasicVar, temp_meshes,
					    temp_displs, temp_bins, basic_vars,
					    radius, doLocalSmooth, &temp_smooth_vars[j]);
	fc_exitIfErrorPrintf(rc, "Failed to smooth variable");
      }
      
      // then do the seq vars by step
      if (numSeqVar > 0) {
	// created smoothed sequence variable
	for (i = 0; i < numSeqVar; i++) {
	  temp_meshes[i] = meshes[seq_meshIDs[i]];
	  temp_bins[i] = vertexBins[seq_meshIDs[i]];
	  temp_displs[i] = seqDispls[seq_meshIDs[i]][j];
	  temp_seq_vars[i] = seq_vars[i][j];
	}
	rc = fc_displacedGeomSmoothVariable(numSeqVar, temp_meshes, 
					    temp_displs, temp_bins,
					    temp_seq_vars, radius, doLocalSmooth,
					    &temp_smooth_seq_vars[j]);
	fc_exitIfErrorPrintf(rc, "Failed to smooth sequence variable");
      }
      
      // cleanup
      for (i = 0; i < numMesh; i++) {
	if (vertexBins[i])
	  fc_freeVertexBin(vertexBins[i]);
      }
    }

    // stich into seq var
    if (numBasicVar > 0) {
      FC_Variable temp_seqVar[numStep];
      numSpecialSeqVar = numBasicVar;
      special_seq_vars = malloc(numBasicVar*sizeof(FC_Variable*));
      if (!special_seq_vars)
	fc_exitIfError(FC_MEMORY_ERROR);
      for (i = 0; i < numBasicVar; i++) {
	for (j = 0; j < numStep; j++)
	  temp_seqVar[j] = temp_smooth_vars[j][i];
	rc = fc_convertVariablesToSeqVariable(numStep, temp_seqVar, sequence,
					      NULL, &special_seq_vars[i]);
	fc_exitIfErrorPrintf(rc, "Failed to convert series of variables to a "
			     "seq variable");
      }
    }
    if (numSeqVar > 0) {
      FC_Variable temp_seqVar[numStep];
      smooth_seq_vars = malloc(numSeqVar*sizeof(FC_Variable*));
      for (i = 0; i < numSeqVar; i++) {
	for (j = 0; j < numStep; j++)
	  temp_seqVar[j] = temp_smooth_seq_vars[j][i];
	rc = fc_convertVariablesToSeqVariable(numStep, temp_seqVar, sequence,
					      NULL, &smooth_seq_vars[i]);
	fc_exitIfErrorPrintf(rc, "Failed to convert series of variables to a "
			     "seq variable");
      }
    }

    // cleanup
    for (j = 0; j < numStep; j++) {
      free(temp_smooth_vars[j]);
      free(temp_smooth_seq_vars[j]);
    }
  }

  // Handle replace request
  if (doReplace) {
    // change name of new vars
    if (smooth_vars) { // may not exist
      for (i = 0; i < numBasicVar; i++) {
	rc = fc_changeVariableName(smooth_vars[i], var_name);
	fc_exitIfErrorPrintf(rc, "Failed to change variable name");
      }
    }
    for (i = 0; i < numSpecialSeqVar; i++) {
      rc = fc_changeSeqVariableName(numStep, special_seq_vars[i],
				    var_name);
      fc_exitIfErrorPrintf(rc, "Failed to change seq variable name");
    }
    for (i = 0; i < numSeqVar; i++) {
      rc = fc_changeSeqVariableName(numStep, smooth_seq_vars[i],
				    var_name);
      fc_exitIfErrorPrintf(rc, "Failed to change seq variable name");
    }
    // delete the old vars
    for (i = 0; i < numBasicVar; i++) 
      fc_deleteVariable(basic_vars[i]);
    for (i = 0; i < numSeqVar; i++)
      fc_deleteSeqVariable(numStep, seq_vars[i]);
  }
      
  // write the new dataset
  rc = fc_getDatasetFileIOType(ds_orig, &fileType);
  fc_exitIfErrorPrintf(rc, "Failed to get file IO type");
  rc = fc_rewriteDataset(ds_orig, new_file_name, fileType);
  fc_exitIfErrorPrintf(rc, "Failed to write dataset");
  
  // cleanup
  free(meshes);
  free(basic_meshIDs);
  free(seq_meshIDs);
  free(meshUsedFlags);
  free(basic_vars);
  free(smooth_vars);
  for (i = 0; i < numSpecialSeqVar; i++)
    free(special_seq_vars[i]);
  free(special_seq_vars);
  for (i = 0; i < numSeqVar; i++) {
    free(seq_vars[i]);
    free(smooth_seq_vars[i]);
  }
  free(seq_vars);
  free(smooth_seq_vars);
  if (doDispl) {
    free(displs);
    for (i = 0; i < numMesh; i++) {
      if (seqDispls[i]) free(seqDispls[i]);
    }
    free(seqDispls);
  }

  // done
  fc_finalLibrary();
    
  exit(0);
}
