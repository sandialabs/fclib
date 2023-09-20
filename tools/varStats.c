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
 * \file varStats.c
 * \brief Report variable statistics on a dataset
 *
 * $Source: /home/Repositories/fcdmf/fclib/tools/varStats.c,v $
 * $Revision: 1.24 $ 
 * $Date: 2006/11/07 23:49:03 $
 *
 * \description
 *    Usage: varStats [options] dataset [var_name]
 *
 *    Reports variable statistics. If a variable name is given, it reports
 *    the stats for only that variable, otherwise it reports stats for
 *    all variables.
 *
 * \modifications
 *   - 04/20/04 WSK, created.
 */

#include <string.h>
#include "fc.h"

#define INDENT_STR "    "

static void get_and_print_stats(FC_Variable var, char* message, 
                                char* indent_str) {
  FC_ReturnCode rc;
  double min, max, ave, stdev, sum;
  int minID, maxID;
  
  printf("%s%s\n", indent_str, message);
  fflush(NULL);

  rc = fc_getVariableMinMax(var, &min, &minID, &max, &maxID);
  if (rc == FC_SUCCESS) {
    printf("%s%smin = %g, at ID = %d\n", indent_str, INDENT_STR, min, minID);
    printf("%s%smax = %g, at ID = %d\n", indent_str, INDENT_STR, max, maxID);
  }
  else
    printf("%s%sFailed to get min/max statistics\n", indent_str, INDENT_STR);
  fflush(NULL);

  rc = fc_getVariableMeanSdev(var, &ave, &stdev);
  if (rc == FC_SUCCESS) 
    printf("%s%save = %g +/- %g\n", indent_str, INDENT_STR, ave, stdev);
  else
    printf("%s%sFailed to get mean/stdev statistics\n", indent_str, INDENT_STR);
  fflush(NULL);

  rc = fc_getVariableSum(var, &sum);
  if (rc == FC_SUCCESS) 
    printf("%s%ssum = %g\n", indent_str, INDENT_STR, sum);
  else
    printf("%s%sFailed to get sum statistics\n", indent_str, INDENT_STR);
  fflush(NULL);

  return;  
}

int main(int argc, char** argv) {
  FC_ReturnCode rc;
  int i, j, k, m;
  int numMesh, numVar, numSeqVar, *numStepPerVar, numComp;
  FC_Dataset ds;
  FC_Mesh *meshes;
  FC_Variable *vars, **seqVars;
  char* file_name = NULL;
  char* mesh_name;
  char* var_name = NULL;
  FC_VerbosityLevel verbose_level = FC_QUIET;
  int doOnlyOne = 0; // flag for whether to do one (1) or all (0);
  char message[1024]; // FIX? fixed buffer size
  FC_MathType mathtype;

  // handle arguments
  if (argc < 2 || argc > 5) {
  usage:
    printf("usage: %s [options] dataset [var_name]\n", 
           argv[0]);
    printf("options: \n");
    printf("   -h         : print this help message\n");
    printf("   -v         : verbose: print warning and error messages\n");
    printf("   -V         : very verbose: prints log and error messages\n");
    printf("\n");
    printf("Prints statistics of the variables in the dataset. If a\n");
    printf("variable name is provided, prints only a summary for that "
           "variable.\n");
    fflush(NULL);
    exit(-1);
  }
  for (i = 1; i < argc; i++) {
    if (!strcmp(argv[i], "-v")) {
      verbose_level = FC_WARNING_MESSAGES;
    }
    else if (!strcmp(argv[i], "-V")) {
      verbose_level = FC_LOG_MESSAGES;
    }
    else if (!strncmp(argv[i], "-h", 2) || !strncmp(argv[i], "-", 1))
      goto usage;
    else {
      file_name = argv[i];
      if (i+1 < argc) {
        i++;
        doOnlyOne = 1;
        var_name = argv[i];
      }
    }
  }
  if (!file_name)
    goto usage;
    
  
  // init library and load dataset 
  rc = fc_setLibraryVerbosity(verbose_level);
  fc_exitIfError(rc);
  rc = fc_initLibrary();
  fc_exitIfError(rc);
  rc = fc_loadDataset(file_name, &ds);
  fc_exitIfErrorPrintf(rc, "Failed to load dataset '%s'", file_name);

  // Get number of meshes and mesh names
  rc = fc_getMeshes(ds, &numMesh, &meshes);
  fc_exitIfError(rc);

  // header
  printf("Variable statistics summary for dataset '%s'\n", file_name);
  fflush(NULL);
  if (numMesh < 0) {
    printf("No meshes were found\n");
    fflush(NULL);
    exit(0);
  }
  else {
    printf("%d meshes\n", numMesh);
    fflush(NULL);
  }

  // loop over meshes
  for (i = 0; i < numMesh; i++) {

    fc_getMeshName(meshes[i], &mesh_name);

    //------- do basic vars -----------------

    // get var names
    if (doOnlyOne) {
      numVar = 1;
      vars = malloc(sizeof(FC_Variable));
      rc = fc_getOrGenerateUniqueVariableByName(meshes[i],
						var_name,
						&vars[0]);
      fc_exitIfError(rc);
      if (FC_HANDLE_EQUIV(vars[0], FC_NULL_VARIABLE)){
	fc_printfWarningMessage("Problems finding unique variable '%s'",
				var_name);
	numVar = 0;
      }
    }
    else {
      rc = fc_getVariables(meshes[i], &numVar, &vars);
      fc_exitIfError(rc);
    }

    // do var statistics
    for (j = 0; j < numVar; j++) {
      // print var info
      printf("\n");
      fflush(NULL);
      sprintf(message, "Var %d (on Mesh %d: '%s')", j, i, mesh_name);
      fc_printVariable(vars[j], message, 0);
      
      // Handle based on mathtype
      rc = fc_getVariableInfo(vars[j], NULL, NULL, NULL, &mathtype, NULL);
      fc_exitIfError(rc);

      // print stats
      if (mathtype == FC_MT_UNKNOWN) {
        printf("%sCan't do statistics for variables with mathtype %s\n",
               INDENT_STR, fc_getMathTypeText(mathtype));
        fflush(NULL);
      }
      else if (mathtype == FC_MT_SCALAR)
        get_and_print_stats(vars[j], "Data Statistics:", INDENT_STR);
      else {
        FC_Variable magnitude, *components;
        printf("%sData Statistics:\n", INDENT_STR);
        fflush(NULL);
        if (mathtype == FC_MT_VECTOR) {
          rc = fc_createMagnitudeVariable(vars[j], "temp", &magnitude);
          fc_exitIfError(rc);
          get_and_print_stats(magnitude, "Magnitude of Vector:",
                              INDENT_STR INDENT_STR);
          fc_deleteVariable(magnitude);
        }
        rc = fc_createComponentVariables(vars[j], &numComp, &components);
        fc_exitIfError(rc);
        for (k = 0; k < numComp; k++) { 
          sprintf(message, "Component %d:", k);
          get_and_print_stats(components[k], message, INDENT_STR INDENT_STR);
          fc_deleteVariable(components[k]);
        }
        free(components);
      }

      // cleanup
      fc_deleteVariable(vars[j]);
    }
    
    // cleanup 
    free(vars);
    
    //------- do sequence vars -----------------

    // get var names
    if (doOnlyOne) {
      numSeqVar = 1;
      numStepPerVar = malloc(sizeof(int));
      seqVars = malloc(sizeof(FC_Variable*));

      rc = fc_getOrGenerateUniqueSeqVariableByName(meshes[i],
						   var_name,
						   &numStepPerVar[0],
						   &seqVars[0]);
      fc_exitIfError(rc);
      if (!seqVars[0]){
	fc_printfWarningMessage("Problems finding unique variable '%s'",
				var_name);
        numSeqVar = 0;
	numStepPerVar = 0;
      }
    } else {
      rc = fc_getSeqVariables(meshes[i], &numSeqVar,&numStepPerVar,
                              &seqVars);
      fc_exitIfError(rc);
    }

    // do var statistics
    for (j = 0; j < numSeqVar; j++) {
      // FIX? get sequence coords to print with stat data?
      
      // print var info
      printf("\n");
      fflush(NULL);
      sprintf(message, "SeqVar %d (on Mesh %d: '%s')", j, i, mesh_name);
      fc_printVariable(seqVars[j][0], message, 0);
      
      // Handle based on mathtype
      rc = fc_getVariableInfo(seqVars[j][0], NULL, NULL, NULL, &mathtype, 
                              NULL);
      fc_exitIfError(rc);
      if (mathtype == FC_MT_UNKNOWN) {
        printf("%sCan't do statistics for variables with mathtype %s\n",
               INDENT_STR, fc_getMathTypeText(mathtype));
        fflush(NULL);
      }
      
      // print statistics for each step
      for (k = 0; k < numStepPerVar[j]; k++) {
        if (mathtype == FC_MT_SCALAR) {
          sprintf(message, "Data Statistics for Step %d:", k);
          get_and_print_stats(seqVars[j][k], message, INDENT_STR);
        }
        else {
          FC_Variable magnitude, *components;
          printf("%sData Statistics for Step %d\n", INDENT_STR, k);
          fflush(NULL);
          if (mathtype == FC_MT_VECTOR) {
            rc = fc_createMagnitudeVariable(seqVars[j][k], "temp", &magnitude);
            fc_exitIfError(rc);
            sprintf(message, "Magnitude of Vector:");
            get_and_print_stats(magnitude, message, INDENT_STR INDENT_STR);
            fc_deleteVariable(magnitude);
          }
          rc = fc_createComponentVariables(seqVars[j][k], &numComp, &components);
          fc_exitIfError(rc);
          for (m = 0; m < numComp; m++) {
            sprintf(message, "Component %d:", m);
            get_and_print_stats(components[m], message, INDENT_STR INDENT_STR);
            fc_deleteVariable(components[m]);
          }
          free(components);
        }
      }
      
      // cleanup
      fc_deleteSeqVariable(numStepPerVar[j], seqVars[j]);
      free(seqVars[j]);
    }
    
    // cleanup
    free(numStepPerVar);
    free(seqVars);
    
    // extra reporting
    if (doOnlyOne && numVar == 0 && numSeqVar == 0) {
      printf("Variable '%s' not found on Mesh '%s'\n", var_name, 
             mesh_name);
      fflush(NULL);
    }
    
    // delete mesh
    fc_deleteMesh(meshes[i]);
    free(mesh_name);
  }
  free(meshes);
    
  // shut down
  fc_deleteDataset(ds);
  fc_finalLibrary();

  exit(0);
}
