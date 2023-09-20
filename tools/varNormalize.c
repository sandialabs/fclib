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
 * \file varNormalize.c
 * \brief Create a new dataset in which the desired variable has
 *        been normalized.
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/tools/varNormalize.c,v $
 * $Revision: 1.5 $ 
 * $Date: 2007/03/14 01:47:22 $
 *
 * \description
 *
 *    Usage: varNormalize [options] dataset "variable name" new_file
 *
 *    Create a new dataset in which the desired variable has been
 *    normalized to be the ratio of the original values to a maximum value.
 *
 *    If you specify "von mises stress" or "pressure" as the var and
 *    it does not exist in the file, it will be calculated.
 *
 *    The maximum value(s) can be provided or "discovered"--the discovered
 *    values are calculated from variable itself. The default execution is
 *    that a single value is applied to all meshes and all steps. If a value
 *    is not provided on the command line, the maximum value of the variable
 *    over all meshes and time steps is used.
 *
 *    Command line arguments can switch the behavior to applying per mesh
 *    values, per time step values, or both. For example, if per mesh behavior
 *    is requested, the maximum value over all timesteps is found for each mesh
 *    and applied to only that mesh.
 *
 *    To apply a provided value to all meshes and steps, you can provide
 *    the value on the command line (-max #). If you want to provide per
 *    mesh values, you give the name of a file containing the mesh names and
 *    their values. The format of the file is pair of lines: the first line
 *    has the name of the mesh exactly (NO trailing spaces), the second line
 *    is the value you would like applied to that mesh.
 *
 *    Currently, this tool does not allow providing per time max values.
 *
 *    The new dataset will be of the same file type as the original dataset.
 *
 * \todo 
 *    - Be able to specify per step "provided" max values (both per mesh 
 *      and overall). For the per mesh, could extend the current file
 *      format to have multiple vals on the val line.
 *    - The max stress tool calculates the tension and compression
 *      separately and reports min and max for those separately.
 *      Will we want that here?
 *    - I believe this does not check for the normalized var already
 *      exisitng by name - fix this.
 *    - What happens for divide by zero? check this
 *    - The von mises stress is fairly well tested via reliance on values
 *      obtained by maxStress, but that is not true of pressure. Decide how
 *      necessary it is to rigorously test this givent hat the var is built the
 *      same in either case and that the vals come from the unittested function
 *      fc_createStressAndPressureVariables.
 *    - for sparse keep track of remaining vars per sequence and delete sequences
 *      with no seqvars
 *
 * \modifications
 *   - 11/17/06 WSD, created
 *   - 03/12/07 ACG option to calculate von mises stress and pressure vars. Much
 *                  of this was lifted from the maxstress tool. If you make a
 *                  change in either file, please cross reference.
 *   - 10/11/07 ACG changed append of "normalized" to "nmz" because of
 *                  var name length limitation in exodus
 */

#include <string.h>
#include "fc.h"

static FC_ReturnCode createStressOrPressureVar(int, FC_Mesh*, char**, char*);

// the endings of the sigma matrix in the order we plan to process them
static char sigmaPostfix[6][3] = { "xx", "yy", "zz", "xy", "yz", "zx" };
static char* von_mises_stress_name = "von mises stress";
static char* pressure_name = "pressure";


FC_ReturnCode createStressOrPressureVar(int numMesh, FC_Mesh *meshes, 
				       char **mesh_names, char *var_name){
  //lifting this from maxstress
  int numStep;
  FC_Variable **sigmas[6]; //index as sigmas[compID][meshID][stepID]
  FC_Sequence sequence;
  FC_ReturnCode rc;
  int i,j,k;

  for (i = 0; i < 6; i++) {
    char charBuf[1028];
    sprintf(charBuf, "sigma-%s", sigmaPostfix[i]);
    sigmas[i] = (FC_Variable**)malloc(numMesh*sizeof(FC_Variable*));
    if (!sigmas[i])
      fc_exitIfError(FC_MEMORY_ERROR);
    for (j = 0; j < numMesh; j++) {
      rc = fc_getOrGenerateUniqueSeqVariableByName(meshes[j], charBuf,
						   &numStep,
						   &sigmas[i][j]);
      fc_exitIfErrorPrintf(rc, "Could not find variable '%s' on mesh '%s'",
			   charBuf, mesh_names[j]);
      if (sigmas[i][j] == 0){
	fc_exitIfErrorPrintf(FC_INPUT_ERROR,
			     "Problems finding (unique) variable '%s' on mesh '%s'",
			     charBuf, mesh_names[j]);
      }
    }
  }

  rc = fc_getSequenceFromSeqVariable(numStep, sigmas[0][0], &sequence);
  fc_exitIfErrorPrintf(rc, "Failed to get sequence from first seq variable");


  for (i = 0; i < numMesh; i++){
    FC_Variable stress[numStep], pressure[numStep];
    FC_Variable *newseq_var;

    for (j = 0; j < numStep; j++) {
      FC_Variable sigs[6];
      for (k = 0; k < 6; k++)
	sigs[k] = sigmas[k][i][j];
      
      // get stress & pressure
      rc = fc_createStressPressureVariables(sigs, &stress[j], &pressure[j]);
      fc_exitIfErrorPrintf(rc, "failed to compute stress & pressure");
    }

    //make one into a seq var and kill the other one
    if (!strcmp(var_name,von_mises_stress_name)){
      rc = fc_convertVariablesToSeqVariable (numStep, stress, sequence,
					     von_mises_stress_name,
					     &newseq_var);
      for (k = 0; k < numStep; k++){
	fc_deleteVariable(pressure[k]);
      }
    } else{
      rc = fc_convertVariablesToSeqVariable (numStep, pressure, sequence,
					    pressure_name, &newseq_var);
      for (k = 0; k < numStep; k++){
	fc_deleteVariable(stress[k]);
      }
    }

    //we dont need to keep our array of handles (the var itself stays)
    free(newseq_var);
  }

  for (i = 0; i < 6; i++){
    for (j = 0; j < numMesh; j++){
      free(sigmas[i][j]);
    }
    free(sigmas[i]);
  }


  return FC_SUCCESS;
}


 
int main(int argc, char** argv) {
  FC_ReturnCode rc;
  int i, j;
  int numMesh, numBasicVar, numSeqVar, numStep = 0;
  FC_Dataset ds_orig;
  FC_Mesh* meshes;
  char** meshNames;
  FC_Variable* basic_vars = NULL, *new_basic_vars = NULL;
  FC_Variable** seq_vars = NULL, **new_seq_vars = NULL;
  FC_Sequence sequence;
  char* file_name = NULL, *new_file_name;
  char* var_name, *new_var_name;
  FC_VerbosityLevel verbose_level = FC_QUIET; 
  FC_FileIOType fileType;
  int doPerMesh = 0; // 0 = same value for all meshes
  int doPerStep = 0; // 0 = same value for all steps
  int *basic_meshIDs, *seq_meshIDs;
  int maxWasProvided = 0; // 1 = a max value was provided
  double providedMaxVal;
  int maxFileWasProvided = 0; // 1 = a file of max values was provided
  char* maxFileName = NULL;
  int* haveMaxValPerMesh;
  double* providedMaxValsPerMesh;
  int doSparse = 0;

  // handle arguments
  if (argc < 4) {
  usage:
    printf("usage: %s [options] dataset \"var name\" new_file_name\n", 
           argv[0]);
    printf("options: \n");
    printf("  --permesh  : Use a different max value for each mesh\n"); 
    printf("  --perstep  : Use a different max value for each step\n");
    printf("  --max #    : Provide a single max value to use for all meshes and steps.\n"); 
    printf("  -f file    : A file to read multiple provided max values from.\n");
    printf("  --sparse   : Create sparse output--contains only the new variable\n");
    printf("  -h         : Print this help message\n");
    printf("  -v         : Verbose: print warning and error messages\n");
    printf("  -V         : Very verbose: print log and error messages\n");
    printf("\n");
    printf("Create a new dataset which will be a copy of the original dataset\n");
    printf("plus an additional variable which is a normalized version of the\n");
    printf("specified variable. The user can provide the max value(s) to normalize\n");
    printf("with using the --max or -f arguments. Otherwise, a max value is calucated\n");
    printf("from the variable field itself. The -f flag requires the --permesh flag.\n");
    printf("Both the --permesh and --perstep flags can be used at the same time.\n");
    printf("Currently, the user cannot provide data for --perstep execution.\n");
    printf("\n");
    printf("The format for the file of max values is to have the name of a mesh on\n");
    printf("a single line (no extra spaces), and then the max value for that mesh\n");
    printf("on the next line. Then the next mesh name on the next line, and its\n");
    printf("value on the line after that.\n");
    printf("\n");   
    printf("If the variable name is either von mises stress or pressure, it will\n");
    printf("be calculated if it does not exist.\n");
    fflush(NULL);
    exit(-1);
  }
  for (i = 1; i < argc; i++) {
    if (!strcmp(argv[i], "--permesh")) {
      doPerMesh = 1;
    }
    else if (!strcmp(argv[i], "--perstep")) {
      doPerStep = 1;
    }
    else if (!strcmp(argv[i], "--max")) {
      if (i+1 >= argc)
        goto usage;
      maxWasProvided = 1;
      providedMaxVal = atof(argv[i+1]);
      i++;
    }
    else if (!strcmp(argv[i], "-f")) {
      if (i+1 >= argc)
        goto usage;
      maxFileWasProvided = 1;
      maxFileName = argv[i+1];
      i++;
    }
    else if (!strcmp(argv[i], "--sparse")) {
      doSparse = 1;
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
      if (i+2 >= argc)
        goto usage;
      file_name = argv[i];
      var_name = argv[i+1];
      new_file_name = argv[i+2];
      i+=2;
    }
  }
  if (!file_name || !var_name || !new_file_name || 
      (maxFileWasProvided && !maxFileName))
    goto usage;
  if (maxWasProvided && maxFileWasProvided) 
  fc_exitIfErrorPrintf(FC_ERROR, "You cannot provide a single max value and "
                       "a file of per mesh max values at the same time");

  // init library and load dataset & meshes
  rc = fc_setLibraryVerbosity(verbose_level);
  fc_exitIfError(rc);
  rc = fc_initLibrary();
  fc_exitIfError(rc);
  rc = fc_loadDataset(file_name, &ds_orig);
  fc_exitIfError(rc);

  // Read the meshes & store their names
  rc = fc_getMeshes(ds_orig, &numMesh, &meshes);
  fc_exitIfError(rc);
  meshNames = (char**)malloc(numMesh*sizeof(char**));
  for (i = 0; i < numMesh; i++) {
    rc = fc_getMeshName(meshes[i], &meshNames[i]);
    fc_exitIfErrorPrintf(rc, "Failed to get mesh name");
  }
 
  // read the values from provided file & associate with a mesh
  // Requires uniquely named meshes
  if (maxFileWasProvided) {
    char temp_name[1028], temp_valstr[1028];
    FILE* file;
    providedMaxValsPerMesh = (double*)malloc(numMesh*sizeof(double));
    haveMaxValPerMesh = (int*)calloc(numMesh, sizeof(int));
    if (!providedMaxValsPerMesh)
      fc_exitIfError(FC_MEMORY_ERROR);
    file = fopen(maxFileName, "r");
    if (!file)
      fc_exitIfErrorPrintf(FC_FILE_IO_ERROR, "Could not open file '%s'",
                           maxFileName);
    while(fscanf(file, "%[^\n]\n", temp_name) > 0) {
      if (fscanf(file, "%[^\n]\n", temp_valstr) < 1) {
        fc_exitIfErrorPrintf(FC_FILE_IO_ERROR, "Could not find the value "
                             "for mesh '%s'", temp_name);
      }
      for (i = 0; i < numMesh; i++) {
        if (!strcmp(temp_name, meshNames[i])) {
          haveMaxValPerMesh[i] = 1;
          providedMaxValsPerMesh[i] = atof(temp_valstr);
          break;
        }
      }
    }
    fclose(file);
    // require that we have something for every mesh
    for (i = 0; i < numMesh; i++) {
      if (haveMaxValPerMesh[i] != 1) {
        fc_exitIfErrorPrintf(FC_ERROR, "In file '%s', a max val was not "
                             "provided for mesh '%s'", maxFileName, 
                             meshNames[i]);
      }
    }
  }


  if ((!strcmp(var_name, von_mises_stress_name)) || 
      (!strcmp(var_name,pressure_name))){
    rc = createStressOrPressureVar(numMesh, meshes, meshNames,var_name);
    fc_exitIfErrorPrintf(rc, "failed to create von mises stress or pressure vars");
  }

 
  // --- loop over meshes, looking for var and or seq var to be normalized
  basic_vars = (FC_Variable*)malloc(numMesh*sizeof(FC_Variable));
  seq_vars = (FC_Variable**)malloc(numMesh*sizeof(FC_Variable*));
  basic_meshIDs = (int*)malloc(numMesh*sizeof(int));
  seq_meshIDs = (int*)malloc(numMesh*sizeof(int));
  if ( !basic_vars || !seq_vars || !basic_meshIDs || !seq_meshIDs)
    fc_exitIfError(FC_MEMORY_ERROR);
  numBasicVar = 0;
  numSeqVar = 0;
  for (i = 0; i < numMesh; i++) {
    int tempNumStep;
    // try to access as basic var
    rc = fc_getOrGenerateUniqueVariableByName(meshes[i], var_name,
					      &basic_vars[numBasicVar]);
    fc_exitIfError(rc);
    if (!FC_HANDLE_EQUIV(basic_vars[numBasicVar], FC_NULL_VARIABLE)){
      basic_meshIDs[numBasicVar] = i;
      numBasicVar++;
    }
    
    // try to access as seq vars

    rc = fc_getOrGenerateUniqueSeqVariableByName(meshes[i], var_name,
						 &tempNumStep,
						 &seq_vars[numSeqVar]);
    fc_exitIfError(rc);
    if (seq_vars[numSeqVar] != NULL){
      if (numSeqVar == 0) {
	//this is the first one
        numStep = tempNumStep;
	rc = fc_getSequenceFromSeqVariable(numStep, seq_vars[numSeqVar],
					   &sequence);
	fc_exitIfErrorPrintf(rc, "Could not get squence from seq variable");
      } else if (tempNumStep != numStep) {
        fc_exitIfErrorPrintf(FC_ERROR, "Expected the sequence variable to "
                             "have the same number of steps on all meshes");
      }
      seq_meshIDs[numSeqVar] = i;
      numSeqVar++;
    }
  }
  if (numBasicVar < 1 && numSeqVar < 1) {
    fc_exitIfErrorPrintf(FC_ERROR, "Failed to find variable '%s' in dataset "
                         "'%s'", var_name, file_name);
  }

  // --- Make new var name - append "_nmz"
  new_var_name = (char*)malloc((strlen(var_name)+12)*sizeof(char));
  if (!new_var_name)
    fc_exitIfError(FC_MEMORY_ERROR);
  sprintf(new_var_name, "%s_nmz", var_name);
  
  // --- Do basic var
  if (numBasicVar > 0) {
    double maxVal;
    double* maxVals;
    new_basic_vars = (FC_Variable*)malloc(numBasicVar*sizeof(FC_Variable));
    if (maxWasProvided) {
      maxVal = providedMaxVal;
    }
    else if (maxFileWasProvided) {
      maxVals = (double*)malloc(numBasicVar*sizeof(double));
      if (!maxVals)
        fc_exitIfError(FC_MEMORY_ERROR);
      for (i = 0; i < numBasicVar; i++) {
        maxVals[i] = providedMaxValsPerMesh[basic_meshIDs[i]];
      }
    }
    else {
      maxVals = (double*)malloc(numBasicVar*sizeof(double));
      if (!maxVals)
        fc_exitIfError(FC_MEMORY_ERROR);
      for (i = 0; i < numBasicVar; i++) {
        rc = fc_getVariableMinMax(basic_vars[i], NULL, NULL, &maxVals[i], NULL);
        fc_exitIfErrorPrintf(rc, "Failed to get variable min max");
      }
      if (!doPerMesh) {
        for (i = 0; i < numBasicVar; i++) {
          if (i == 0)
            maxVal = maxVals[i];
          else if (maxVals[i] > maxVal)
            maxVal = maxVals[i];
        }
      } 
    }
    for (i = 0; i < numBasicVar; i++) {
      double temp_max;
      if (doPerMesh)
        temp_max = maxVals[i];
      else
        temp_max = maxVal;
      printf("Max value for basic var on mesh '%s' = %g\n",
             meshNames[basic_meshIDs[i]], temp_max);
      fflush(NULL);
      rc = fc_varOperatorConst(basic_vars[i], "/", 1, &temp_max, FC_DT_DOUBLE,
                               new_var_name, &new_basic_vars[i]);
      fc_exitIfErrorPrintf(rc, "Failed to divide var by const value");
    }
    if (!maxWasProvided)
      free(maxVals);
  }

  // --- Do seq var
  if (numSeqVar > 0) {
    new_seq_vars = (FC_Variable**)malloc(numSeqVar*sizeof(FC_Variable*));
    if (doPerStep) {
      double* seqMaxVal = (double*)malloc(numStep*sizeof(double));
      double** seqMaxVals = (double**)malloc(numSeqVar*sizeof(double*));
      if (!seqMaxVals || !seqMaxVals)
        fc_exitIfError(FC_MEMORY_ERROR);
      for (i = 0; i < numSeqVar; i++) {
        seqMaxVals[i] = (double*)malloc(numStep*sizeof(double));
        if (!seqMaxVals[i])
          fc_exitIfError(FC_MEMORY_ERROR);
      }
      for (i = 0; i < numSeqVar; i++) {
        for (j = 0; j < numStep; j++) {
          rc = fc_getVariableMinMax(seq_vars[i][j], NULL, NULL,
                                    &seqMaxVals[i][j], NULL);
          fc_exitIfErrorPrintf(rc, "Failed to get seq var step's min max");
        }
      }
      if (!doPerMesh) {
        for (j = 0; j < numStep; j++) {
          for (i = 0; i < numSeqVar; i++) {
            if (i == 0) 
              seqMaxVal[j] = seqMaxVals[i][j];
            else if (seqMaxVals[i][j] > seqMaxVal[j])
              seqMaxVal[j] = seqMaxVals[i][j];
          }
        }
      }
      for (i = 0; i < numSeqVar; i++) {
        double* temp_seqMax;
        if (doPerMesh)
          temp_seqMax = seqMaxVals[i];
        else
          temp_seqMax = seqMaxVal;
        printf("Max values for seq var on mesh '%s' = \n",
               meshNames[seq_meshIDs[i]]);
        for (j = 0; j < numStep; j++)
          printf("  %d: %g\n", j, temp_seqMax[j]);
        rc = fc_seqVarOperatorSeqConst(numStep, seq_vars[i], "/", 
                                       1, temp_seqMax, FC_DT_DOUBLE, 
                                       new_var_name, &new_seq_vars[i]);
        fc_exitIfErrorPrintf(rc, "Failed to divide seq var by const values");
      }
      for (i = 0; i < numSeqVar; i++)
        free(seqMaxVals[i]);
      free(seqMaxVals);
      free(seqMaxVal);
    }
    else {
      double maxVal;
      double* maxVals;
      if (maxWasProvided) {
        maxVal = providedMaxVal;
      }
      else if (maxFileWasProvided) {
        maxVals = (double*)malloc(numSeqVar*sizeof(double));
        if (!maxVals)
          fc_exitIfError(FC_MEMORY_ERROR);
        for (i = 0; i < numSeqVar; i++) {
          maxVals[i] = providedMaxValsPerMesh[seq_meshIDs[i]];
        }
      }
      else {
        maxVals = (double*)malloc(numSeqVar*sizeof(double));
        if (!maxVals)
          fc_exitIfError(FC_MEMORY_ERROR);
        for (i = 0; i < numSeqVar; i++) {
          rc = fc_getSeqVariableMinMax(numStep, seq_vars[i], NULL, NULL, NULL,
                                       &maxVals[i], NULL, NULL);
          fc_exitIfErrorPrintf(rc, "Failed to get seq variable min max");
        }
        if (!doPerMesh) {
          for (i = 0; i < numSeqVar; i++) {
            if (i == 0)
              maxVal = maxVals[i];
            else if (maxVals[i] > maxVal)
              maxVal = maxVals[i];
          }
        }
      }
      for (i = 0; i < numSeqVar; i++) {
        double temp_max;
        if (doPerMesh)
          temp_max = maxVals[i];
        else
          temp_max = maxVal;
        printf("Max value for seq var on mesh '%s' = %g\n",
               meshNames[seq_meshIDs[i]], temp_max);
        rc = fc_seqVarOperatorConst(numStep, seq_vars[i], "/", 
                                    1, &temp_max, FC_DT_DOUBLE, 
                                    new_var_name, &new_seq_vars[i]);
        fc_exitIfErrorPrintf(rc, "Failed to divide seq var by const value");
      }
      if (!maxWasProvided)
        free(maxVals);
    }
  }

  // --- If requested, make output sparse by deleting all but the
  // new variable. note that in exodus we always keep the sequence
  // even if there are no more seq vars on that sequence
  if (doSparse) {
    int temp_num, *temp_numSteps;
    FC_Subset *temp_subsets;
    FC_Variable *temp_vars, **temp_seqVars; 
    // because vars are in order by mesh, we don't have to search
    int i_var = 0; // current basic var
    int i_seqVar = 0; // current seq var

    fc_getGlobalVariables(ds_orig, &temp_num, &temp_vars);
    for (i = 0; i < temp_num; i++)
      fc_deleteVariable(temp_vars[i]);
    free(temp_vars);
    fc_getGlobalSeqVariables(ds_orig, &temp_num, &temp_numSteps,
                             &temp_seqVars);
    for (i = 0; i < temp_num; i++) {
      fc_deleteSeqVariable(temp_numSteps[i], temp_seqVars[i]);
      free(temp_seqVars[i]);
    }
    free(temp_numSteps);
    free(temp_seqVars);
    for (i = 0; i < numMesh; i++) {
      fc_getSubsets(meshes[i], &temp_num, &temp_subsets);
      for (j = 0; j < temp_num; j++) 
        fc_deleteSubset(temp_subsets[j]);
      free(temp_subsets);
      fc_getVariables(meshes[i], &temp_num, &temp_vars);
      for (j = 0; j < temp_num; j++) {
        if (numBasicVar > 0 && 
            FC_HANDLE_EQUIV(temp_vars[j], new_basic_vars[i_var])) {
          continue;
        }
        fc_deleteVariable(temp_vars[j]);
      }
      free(temp_vars);
      fc_getSeqVariables(meshes[i], &temp_num, &temp_numSteps, &temp_seqVars);
      for (j = 0; j < temp_num; j++) {
        if (numSeqVar > 0 &&
            FC_HANDLE_EQUIV(temp_seqVars[j][0], new_seq_vars[i_seqVar][0])) {
          ; // nothing
        }
        else {
          fc_deleteSeqVariable(temp_numSteps[j], temp_seqVars[j]);
        }
        free(temp_seqVars[j]);
      }
      free(temp_numSteps);
      free(temp_seqVars);
      if (numBasicVar > 0 && basic_meshIDs[i_var] == i)
        i_var++;
      if (numSeqVar > 0 && seq_meshIDs[i_seqVar] == i)
        i_seqVar++;
    }
  }

  // --- Write the new dataset
  rc = fc_getDatasetFileIOType(ds_orig, &fileType);
  fc_exitIfErrorPrintf(rc, "Failed to get file IO type");
  if (fileType == FC_FT_LSDYNA)
    fileType = FC_FT_EXODUS;
  rc = fc_rewriteDataset(ds_orig, new_file_name, fileType);
  fc_exitIfErrorPrintf(rc, "Failed to write dataset");

  // --- done

  // cleanup
  for (i = 0; i < numMesh; i++)
    free(meshNames[i]);
  free(meshNames);
  free(meshes);
  free(basic_vars);
  free(new_basic_vars);
  for (i = 0; i < numSeqVar; i++) {
    free(seq_vars[i]);
    free(new_seq_vars[i]);
  }
  free(seq_vars);
  free(new_seq_vars);

  // final library
  fc_finalLibrary();
 
  exit(0);
}
