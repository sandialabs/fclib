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
 * \file loadHistory.c
 * \brief Calculate the load history on a subset of a mesh.
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/tools/loadHistory.c,v $
 * $Revision: 1.10 $ 
 * $Date: 2006/09/28 06:34:39 $
 *
 * \description
 *    Usage: loadHistory [options] filename meshname subsetname varname
 *
 *    Computes a curve where each point is the sum of the variable over
 *    the subset at that time. If the variable has multiple components,
 *    this is done for each component. Also, if the variable is a vector,
 *    this is done for the magnitude.
 *
 * \todo What if variable is tensor? Currently program stops with error.
 *   Could still do the components...
 *
 * \modifications
 *   - 01/19/05 WSK, created.
 */

#include <string.h>
#include <time.h>
#include "fc.h"
 
int main(int argc, char** argv) {
  FC_ReturnCode rc;
  time_t start_time = time(NULL); // the time analysis starts
  int i, j;
  int numStep, numComp;
  FC_Dataset dataset;
  FC_Sequence sequence;
  FC_Mesh mesh, *returnMeshes;
  FC_Subset subset, *returnSubsets;
  FC_Variable *seqVar, **compSeqVars;
  int numReturnMeshes, numReturnSubsets;
  char* file_name = NULL;
  char* mesh_name = NULL;
  char* subset_name = NULL;
  char* var_name = NULL;
  char* output_root_name = "loadHistory"; // this is default, can change
  int str_len;
  char* datfile_name;
  char* gnufile_name;
  char* plotfile_name;
  char** comp_names;
  FC_VerbosityLevel verbose_level = FC_QUIET;
  FC_AssociationType subset_assoc, var_assoc;
  double **sums;
  FILE *file;
  void* times;
  FC_DataType datatype;

  // handle arguments
  if (argc < 5) {
  usage:
    printf("\n");
    printf("usage: %s [options] filename meshname subsetname varname\n", 
           argv[0]);
    printf("options: \n");
    printf("   -o         : root name for output files (default is "
           "loadHistory)\n");
    printf("   -h         : print this help message\n");
    printf("   -v         : print warning and error messages\n");   
    printf("   -v         : print log and error messages\n");   
    printf("\n");
    printf("Creates a file, '*.dat', that contains the sum of the variable "
           "over the subset\n" );
    printf("(decomposed into components if the variable is a vector) for "
           "each time\n");
    printf("step. It also creates a file, '*.gnu', which produces a plot of "
           "the results. \n");
    printf("\n");
    fflush(NULL);
    exit(-1);
  }
  for (i = 1; i < argc; i++) {
    if (!strcmp(argv[i], "-o")) {
      output_root_name = argv[i+1];
      i++;
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
      if (i + 3 >= argc)
        goto usage;
      file_name = argv[i];
      mesh_name = argv[i+1];
      subset_name = argv[i+2];
      var_name = argv[i+3];
      i+=4;
    }
  }
  if (!file_name)
    goto usage;

  // init library and load dataset
  rc = fc_setLibraryVerbosity(verbose_level);
  fc_exitIfError(rc);
  rc = fc_initLibrary();
  fc_exitIfError(rc);
  rc = fc_loadDataset(file_name, &dataset);
  fc_exitIfErrorPrintf(rc, "Failed to load dataset '%s'", file_name);

  // get mesh, subset and sequence variable
  rc = fc_getMeshByName(dataset, mesh_name, &numReturnMeshes,
			&returnMeshes);
  fc_exitIfErrorPrintf(rc, "Failed to get mesh '%s' by name", mesh_name);
  if (numReturnMeshes != 1){
    fc_exitIfErrorPrintf(FC_INPUT_ERROR, 
			 "Failed to find (unique) mesh '%s' - found %d matches",
			 mesh_name,numReturnMeshes);
  } 
  mesh = returnMeshes[0];
  free(returnMeshes);
  rc = fc_getSubsetByName(mesh, subset_name, &numReturnSubsets, &returnSubsets);
  fc_exitIfErrorPrintf(rc, "Failed to get subset '%s' by name", subset_name);
  if (numReturnSubsets != 1){
    fc_exitIfErrorPrintf(FC_INPUT_ERROR,
			 "Failed to find (unique) subset '%s' - found %d matches",
			 subset_name, numReturnSubsets);
  }
  subset = returnSubsets[0];
  free(returnSubsets);
  rc = fc_getOrGenerateUniqueSeqVariableByName(mesh, var_name,
					       &numStep, &seqVar);
  fc_exitIfErrorPrintf(rc, "Failed to find or create sequence variable '%s'\n",
		       var_name);
  if (!seqVar){
    fc_exitIfErrorPrintf(FC_INPUT_ERROR,
			 "Failed to find or create sequence variable '%s'\n",
			 var_name);
  }

  rc = fc_getSequenceFromSeqVariable(numStep, seqVar, &sequence);
  fc_exitIfErrorPrintf(rc, "Failed to get sequence from seqVariable");
  rc = fc_getSequenceDataType(sequence, &datatype);
  fc_exitIfErrorPrintf(rc, "Failed to get data from sequence");
  rc = fc_getSequenceCoordsPtr(sequence, &times);
  fc_exitIfErrorPrintf(rc, "Failed to get data from sequence");

  // Make sure that subset and seqVar are on same association
  rc = fc_getSubsetAssociationType(subset, &subset_assoc);
  rc = rc + fc_getVariableAssociationType(seqVar[0], &var_assoc);
  if (subset_assoc != var_assoc || rc != FC_SUCCESS) 
    fc_exitIfErrorPrintf(FC_ERROR, "Subset and seqVar must have the same "
                         "assoc (subset = '%s', var = '%s')",
                         fc_getAssociationTypeText(subset_assoc),
                         fc_getAssociationTypeText(var_assoc));

  // If a vector seqVar, decompose seqVar into components & get magnitude
  rc = fc_getVariableNumComponent(seqVar[0], &numComp);
  fc_exitIfErrorPrintf(rc, "Could not get number of components of variable");
  if (numComp == 1) {
    compSeqVars = (FC_Variable**)malloc(sizeof(FC_Variable*));
    compSeqVars[0] = (FC_Variable*)malloc(numStep*sizeof(FC_Variable));
    for (i = 0; i < numStep; i++)
      compSeqVars[0][i] = seqVar[i];
  }
  else {                        
    void* temp;
    rc = fc_createComponentSeqVariables(numStep, seqVar, &numComp,
                                        &compSeqVars);
    fc_exitIfErrorPrintf(rc, "Could not split var into components");
    temp = realloc(compSeqVars, (numComp+1)*sizeof(FC_Variable*));
    if (!temp) 
      fc_exitIfError(FC_MEMORY_ERROR);
    compSeqVars = temp;
    compSeqVars[numComp] = (FC_Variable*)malloc(numStep*sizeof(FC_Variable));
    if (!compSeqVars[numComp]) 
      fc_exitIfError(FC_MEMORY_ERROR);
    for (i = 0; i < numStep; i++) {
      char* tempName = (char*)malloc((strlen(var_name)+20)*sizeof(char));
      sprintf(tempName, "%s_magnitude", var_name);
      rc = fc_createMagnitudeVariable(seqVar[i], tempName, 
                                      &compSeqVars[numComp][i]);
      free(tempName);
      fc_exitIfErrorPrintf(rc, "could not create magnitude variable");
    }
    numComp++;
  }

  // Collect Sum over subset for all components and steps
  sums = (double**)malloc(numComp*sizeof(double*));
  for (i = 0; i < numComp; i++)
    sums[i] = (double*)malloc(numStep*sizeof(double));
  for (i = 0; i < numComp; i++) {
    for (j = 0; j < numStep; j++) {
      rc = fc_getVariableSubsetSum(compSeqVars[i][j], subset, &sums[i][j]);
      fc_exitIfErrorPrintf(rc, "Get sum failed for step %d of component %d",
                           j, i);
    }
  }

  // Get names of things
  // get component names
  comp_names = (char**)malloc(numComp*sizeof(char*));
  for (i = 0; i < numComp; i++) {
    rc = fc_getVariableName(compSeqVars[i][0], &comp_names[i]);
    fc_exitIfErrorPrintf(rc, "failed to get name of component %d", i);
  }

  // setup output names
  str_len = strlen(output_root_name)+1;
  datfile_name = (char*)malloc((str_len + 4)*sizeof(char));
  sprintf(datfile_name, "%s.dat", output_root_name);
  gnufile_name = (char*)malloc((str_len + 4)*sizeof(char));
  sprintf(gnufile_name, "%s.gnu", output_root_name);
  plotfile_name = (char*)malloc((str_len + 4)*sizeof(char));
  sprintf(plotfile_name, "%s.eps", output_root_name);

  // create dat file
  file = fopen(datfile_name, "w");
  fprintf(file, "# Data file generated by FCLib's %s\n", argv[0]);
  fprintf(file, "# %s", ctime(&start_time));
  fprintf(file, "#\n");
  fprintf(file, "# Database file: %s\n", file_name);
  fprintf(file, "# Mesh = '%s'\n", mesh_name);
  fprintf(file, "# Subset = '%s'\n", subset_name);
  fprintf(file, "# SeqVariable = '%s'\n", var_name);
  fprintf(file, "#\n");
  fprintf(file, "# numComp = %d, numStep = %d\n", numComp, numStep);
  fprintf(file, "# time");
  fflush(NULL);
  for (i = 0; i < numComp; i++) 
    fprintf(file, " %s", comp_names[i]);
  fprintf(file, "\n");
  for (i = 0; i < numStep; i++) {
    switch (datatype) {
    case FC_DT_INT:     fprintf(file, "%11d", ((int*)times)[i]);       break;
    case FC_DT_FLOAT:   fprintf(file, "%11.6g", ((float*)times)[i]);   break;
    case FC_DT_DOUBLE:  fprintf(file, "%11.6g", ((double*)times)[i]);  break;
    default:
      printf("Warning: could not decode time coords, writing time IDs instead\n");
      fprintf(file, "%11d", i);
    }
    for (j = 0; j < numComp; j++) 
      fprintf(file, " %11.6g", sums[j][i]);
    fprintf(file, "\n");
  }
  fflush(NULL);
  fclose(file);

  // create gnuplot file
  file = fopen(gnufile_name, "w");
  if (file == NULL)
    fc_exitIfErrorPrintf(FC_FILE_IO_ERROR, "Couldn't open file");
  fprintf(file, "# Gnuplot script generated by FCLib's %s\n", argv[0]);
  fprintf(file, "# %s", ctime(&start_time));
  fprintf(file, "#\n");
  fprintf(file, "# References data file '%s'\n", datfile_name);
  fprintf(file, "# To run: 'gnuplot %s' \n", gnufile_name); 
  fprintf(file, "#\n");

  fprintf(file, "# set up\n");
  fprintf(file, "set terminal postscript eps color\n");
  fprintf(file, "set output '%s'\n", plotfile_name);
  fprintf(file, "set multiplot\n");
  fprintf(file, "set size 1, %g\n", 1./numComp); // all plots same size
  fprintf(file, "set xlabel 'Time'\n");
  fprintf(file, "\n");
  fflush(NULL);

  //fprintf(file, "set title \"Title\"\n");
  for (i = 0; i < numComp; i++) {
    fprintf(file, "set origin 0, %g\n", 1.- (i+1.)/numComp);
    fprintf(file, "set ylabel '%s'\n", comp_names[i]);
    fprintf(file, "plot '%s' using 1:%d notitle with linespoints\n", 
            datfile_name, 2+i);
  }
  fflush(NULL);
  fclose(file);

  // all done
  free(seqVar);
  for (i = 0; i < numComp; i++) 
    free(compSeqVars[i]);
  free(compSeqVars);
  for (i = 0; i < numComp; i++)
    free(sums[i]);
  free(sums);
  free(datfile_name);
  free(gnufile_name);
  free(plotfile_name);
  fc_finalLibrary();

  exit(0);
}
