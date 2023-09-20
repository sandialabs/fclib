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
 * \file windowVars.c
 * \brief window averages the seq vars on a dataset
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/tools/windowVars.c,v $
 * $Revision: 1.12 $
 * $Date: 2006/09/28 06:34:39 $
 *
 * \description
 *    Usage: windowVars [options] dataset [var_name]
 *
 *    Reports windowaverages. If a variable name is given, it reports
 *    the averages for only that variable, otherwise it reports averages for
 *    all variables.
 *
 *    Code modeled on varStats.
 *
 *
 * \modifications
 *   - 12/21/04 ACG, created.
 *   -  1/14/05 ACG, cleaned out testing stuff since now in checktimeseries.c
 */

#include <string.h>
#include "fc.h"

#define INDENT_STR "    "

int main(int argc, char** argv) {
  FC_ReturnCode rc;
  int i, j,ii,jj,k;
  int numMesh, numSeqVar, *numStepPerVar;
  FC_Dataset ds;
  FC_Mesh *meshes;
  FC_Variable **seqVars, *newseqvar;
  void* data;
  void* data2;
  char* file_name = NULL;
  char* mesh_name;
  char* var_name = NULL;
  FC_VerbosityLevel verbose_level = FC_QUIET;
  int doOnlyOne = 0; // flag for whether to do one (1) or all (0);
  char message[1024]; // FIX? fixed buffer size
  int window;

  // handle arguments
  if (argc < 3 || argc > 6) {
  usage:
    printf("usage: %s [options] dataset window_size [var_name] n", 
           argv[0]);
    printf("options: \n");
    printf("   -h         : print this help message\n");
    printf("   -v         : verbose: print warning and error messages\n");
    printf("   -V         : very verbose: prints log and error messages\n");
    printf("\n");
    printf("Prints window averages of the seq variables in the dataset. If a\n");
    printf("variable name is provided, prints only averages for that "
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
      if (i+1 < argc){
        file_name = argv[i];
        window = atoi(argv[i+1]);
        i+=2;
      }
      if (i < argc) {
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

    if (doOnlyOne){ // get var names
      numSeqVar = 1;
      numStepPerVar = malloc(sizeof(int));
      seqVars = malloc(sizeof(FC_Variable*));

      rc = fc_getOrGenerateUniqueSeqVariableByName(meshes[i],
						   var_name,
						   &numStepPerVar[0],
						   &seqVars[0]);
      fc_exitIfErrorPrintf(rc, "Problems finding unique variable '%s'",
			   var_name);
      if (!seqVars[0]){
	fc_printfWarningMessage("Problems finding unique variable '%s'",
				var_name);
	numSeqVar = 0;
        numStepPerVar = 0;
        printf("No seq variable %s on mesh %s\n",var_name,mesh_name);
        fflush(NULL);
      }
    }else {
      rc = fc_getSeqVariables(meshes[i], &numSeqVar,
                              &numStepPerVar,&seqVars);
      fc_exitIfError(rc);
      if (numSeqVar ==0){
        printf("No seq variables on mesh %s\n",mesh_name);
        fflush(NULL);
      }
    }

    for (jj = 0; jj < numSeqVar; jj++) {
      // print var info
      int numDataPoint, numComponent;
      FC_AssociationType assoc, newassoc;
      FC_MathType mathtype, newmathtype;
      FC_DataType datatype, newdatatype;
      int newnumDataPoint, newnumComponent;

      printf("\n");
      fflush(NULL);
      sprintf(message, "SeqVar %d (on Mesh %d: '%s')", jj, i, mesh_name);
      fc_printVariable(seqVars[jj][0], message, 0);

      rc = fc_getVariableInfo(seqVars[jj][0], &numDataPoint, &numComponent,
                              &assoc, &mathtype, &datatype);

      rc = fc_leadingWindowAverage(numStepPerVar[jj],seqVars[jj],
                                   window,"window",&newseqvar);

      if (rc != FC_SUCCESS){
        fc_printfErrorMessage("failed to compute leading window average");
        exit(-1);
      }

      rc = fc_getVariableInfo(newseqvar[0], &newnumDataPoint, &newnumComponent,
                              &newassoc, &newmathtype, &newdatatype);

      for (ii  = 0; ii < numStepPerVar[jj]; ii++){ //like each time
        fc_getVariableDataPtr(newseqvar[ii],&data);
        fc_getVariableDataPtr(seqVars[jj][ii],&data2);

        printf("step %d\n", ii);
        fflush(NULL);
        for (j = 0; j < newnumDataPoint; j++) { //like each node
          printf("%d: ", j);
          for (k = 0; k < newnumComponent; k++) { //like dimensionality
            if (datatype == FC_DT_INT){
              printf("(%d ", (((int*)data2)[j*newnumComponent + k]));
            } else {
              if (datatype == FC_DT_FLOAT){
                printf("(%f ",(((float*)data2)[j*newnumComponent + k]));
              }else {
                printf("(%f ",(((double*)data2)[j*newnumComponent + k]));
              }
            }
            if (newdatatype == FC_DT_FLOAT){
              printf("%f) ", (((float*)data)[j*newnumComponent + k]));
            }else{
              printf("%f) ", (((double*)data)[j*newnumComponent + k]));
            }
          }
          printf("\n");
          fflush(NULL);
        }
        data = NULL;
        data2 = NULL;
      }
      fc_deleteSeqVariable(numStepPerVar[jj],newseqvar);
      free(newseqvar);
    }

    // cleanup
    free(numStepPerVar);
    free(seqVars);
    
    // delete mesh
    fc_deleteMesh(meshes[i]);
    free(mesh_name);
  }
    
  // shut down
  fc_deleteDataset(ds);
  fc_finalLibrary();
  
  exit(0);
}
