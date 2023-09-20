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
 * \file fc2ensight.c
 * \brief Conversion program
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/tools/fc2ensight.c,v $
 * $Revision: 1.20 $ 
 * $Date: 2006/10/19 03:14:52 $
 *
 * \description 
 *
 *   Opens a dataset & then writes ensight gold files with same root file 
 *   name but different extensions in directory named rootname_ens in
 *   the current directory. (E.g. starting with can.saf you get
 *   can_ens/can.case, can_ens/can.geo, etc.)
 *
 *   Currently necessary because ensight cannot read the type of saf
 *   files produced by fclib. Also, ensight automatically clips dead
 *   elements when reading LsDyna files so this is a helpful workaround
 *   if you don't want that to happen.
 *
 *   Subsets are promoted to parts. Global variabls are not written
 *   unless explicitly asked for with the -g option.
 *
 * \modifications
 *    - 4/13/05  WSK  Created to replace the saf only version.
 *
 * \todo? prefix subset names so it's obvious they are a subset?
 * \todo haven't test cases with subsets or variables with association whole
 * \todo a bunch more stuff indicated in comments throughout code.
 * \todo - do we need this anymore? 10/15/07 ACG. I am not
 *    doing anything to this file becuase I dont knwo what the issues are
 */

#include <string.h>
#include <fc.h>

static char* getEnsightElementType(FC_ElementType elemType) {
  switch(elemType) {
  case FC_ET_POINT:    return "point";
  case FC_ET_LINE:     return "bar2";
  case FC_ET_TRI:      return "tria3";
  case FC_ET_QUAD:     return "quad4";
  case FC_ET_TET:      return "tetra4";
  case FC_ET_PYRAMID:  return "pyramid5";
  case FC_ET_PRISM:    return "penta6";
  case FC_ET_HEX:      return "hexa8";
  case FC_ET_MIXED:    return "mixed";
  case FC_ET_UNKNOWN:  return "unknown";
  }
  return "unknown";
}

// all str lengths are w/o trailing '\0'
int main(int argc, char **argv)
{
  FC_VerbosityLevel verbose_level = FC_QUIET; 
  int i, j, k, m, n, p;
  int flag;
  FC_ReturnCode rc;
  char* orig_file_name;
  char* temp_string, *temp_string2;
  char* path_base_file_name; // path + base file name
  char* base_file_name;      // points into path_base_file_name
  int path_base_length;      // length of path_base_file_name
  int base_length;           // length of base_file_name
  int path_length;           // used to cut off path -> index to start of name
  char* case_file_name; // name for the case file
  char* geom_file_name; // name for the case file
  char* var_file_name;   // name for the case file
  char* seqVar_file_name; // name for the case file
  FILE* case_file;
  FILE* geom_file;
  FILE* var_file;
  FILE* seqVar_file;
  int doGlobals = 0;

  FC_Dataset dataset;
  int numSequence;
  FC_Sequence* sequences;
  int numMesh;
  FC_Mesh* meshes;
  int* numSubsetPerMesh;
  FC_Subset** subsetsPerMesh;  
  int numVar, numSeqVar, *numStepPerSeqVar, numStep;
  FC_Variable var, *variables, *seqVar, **seqVariables;

  FC_ElementType *meshElemTypes;    // numMesh long
  FC_AssociationType **subsetAssocsPerMesh; // numMesh long, numSubsetPerMesh
  int** numSubsetVertsPerMesh, ***subsetVertIDsPerMesh; // numMesh, numSubsetPerMesh
  int** numMemberPerMesh, ***memberIDsPerMesh; // numMesh x numSubsetPerMesh
  int numTotalVar, numTotalSeqVar;
  char **varNames, **seqVarNames; // numTotalVar, numTotalSeqVar long
  int *varNumComps, *seqVarNumComps; // numTotalVar, numTotalSeqVar long
  FC_AssociationType *varAssocs, *seqVarAssocs; // "
  FC_MathType *varMathTypes, *seqVarMathTypes;  // "
  FC_DataType *varDataTypes, *seqVarDataTypes;  // "
  int *seqVarNumSteps; // numTotalSeqVar long
  int *seqVarSeqIDs; // numTotalSeqVar long

  int partID;

  //------------------------------------------------------------------
  //  process command line args
  //------------------------------------------------------------------  
  // FIX? to specify output base_file_name?
  if (argc < 2) {
  usage:
    printf("usage:  %s [options] dataset_file\n", argv[0]);
    printf("options: \n");
    printf("    -g    : output the global values (default is not to)\n");
    printf("    -h    : print this help message\n");
    printf("    -v    : verbose: print warning and error messages\n");
    printf("    -V    : very verbose: print log and error message\n");
    printf("  Creates ensight files in the directory the command is\n");
    printf("executed from. The created files will have the same base\n");
    printf("name as the dataset file.\n"); 
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
    else if (!strcmp(argv[i], "-g")) {
      doGlobals = 1;
    }
    else {
      if (i != argc-1)
        goto usage;
      orig_file_name = argv[i];
    }
  }
  
  //------------------------------------------------------------------
  //  determine base file name
  //------------------------------------------------------------------
  // make a directory and put all files in it
  temp_string = (char*)malloc((strlen(orig_file_name) + 1)*sizeof(char));       
  strcpy(temp_string, orig_file_name);
  // take off outer extension
  for (i = strlen(temp_string); i >= 0; i--) {
    if (temp_string[i] == '/') { // there was no extension
      break;
    }
    if (temp_string[i] == '.') {
      temp_string[i] = '\0';
      break;
    }
  } // if no early breaks, then no extension  
  // take off leading directories (2nd pointer points into first string)
  temp_string2 = temp_string;
  for (i = strlen(temp_string); i >= 0; i--) {
    if (temp_string[i] == '/') {
      temp_string2 = &temp_string[i+1];
      break;
    }
  }
  // make base_name w & w/o path
  base_length = strlen(temp_string2);
  path_length = base_length + 5;
  path_base_length = path_length + base_length;
  path_base_file_name = (char*)malloc((path_base_length+1)*sizeof(char));
  sprintf(path_base_file_name, "%s_ens/%s", temp_string2, temp_string2);
  base_file_name = &path_base_file_name[path_length];
  free(temp_string);
  // debug
  //printf("path_base_file_name = '%s'\n", path_base_file_name);
  //printf("base_file_name = '%s'\n", base_file_name);
  // make directory
  temp_string = (char*)malloc((base_length + 20)*sizeof(char));
  sprintf(temp_string, "mkdir %s_ens", base_file_name);
  flag = system(temp_string);
  if (flag != 0) 
    fc_exitIfErrorPrintf(FC_ERROR, "Cannot create directory '%s_ens'", 
                         base_file_name);
  free(temp_string);

  //---------------------------------------------------------------------
  //  Initialize FCLib & open dataset
  //---------------------------------------------------------------------
  rc = fc_setLibraryVerbosity(verbose_level);
  fc_exitIfErrorPrintf(rc, "Failed to set library verbosity level");
  rc = fc_initLibrary();
  fc_exitIfErrorPrintf(rc, "Failed to initialize library");
  rc = fc_loadDataset(orig_file_name, &dataset);
  fc_exitIfErrorPrintf(rc, "Failed to load dataset '%s'", orig_file_name);

  //---------------------------------------------------------------------
  //  Get handles to sequences, meshes & subsets
  //---------------------------------------------------------------------
  rc = fc_getSequences(dataset, &numSequence, &sequences);
  fc_exitIfErrorPrintf(rc, "Failed to get sequences");
  rc = fc_getMeshes(dataset, &numMesh, &meshes);
  fc_exitIfErrorPrintf(rc, "Failed to get meshes");
  numSubsetPerMesh = (int*)malloc(numMesh*sizeof(int));
  subsetsPerMesh = (FC_Subset**)malloc(numMesh*sizeof(FC_Subset*));
  if (!numSubsetPerMesh || !subsetsPerMesh)
    fc_exitIfError(FC_MEMORY_ERROR);
  for (i = 0; i < numMesh; i++) {
    rc = fc_getSubsets(meshes[i], &numSubsetPerMesh[i], &subsetsPerMesh[i]);
    fc_exitIfErrorPrintf(rc, "Failed to get subsets");
  }
  // FIX remove any empty subsets from list
 
  //---------------------------------------------------------------------
  //  Make master list of all variable and sequence variable names 
  //  (Also collect info about vars)
  //---------------------------------------------------------------------
  // FIX! this assumes that all variables of the same name will be
  //        of the same type (i.e. nodal or elemental).
  numTotalVar = numTotalSeqVar = 0;
  varNames = seqVarNames = NULL;
  varNumComps = seqVarNumComps = NULL;
  varAssocs = seqVarAssocs = NULL;
  varMathTypes = seqVarMathTypes = NULL;
  varDataTypes = seqVarDataTypes = NULL;
  seqVarNumSteps = NULL;
  seqVarSeqIDs = NULL;
  for (i = 0; i < numMesh; i++) {
    rc = fc_getVariables(meshes[i], &numVar, &variables);
    fc_exitIfErrorPrintf(rc, "Failed to get variables");
    for (j = 0; j < numVar; j++) {
      rc = fc_getVariableName(variables[j], &temp_string);
      fc_exitIfErrorPrintf(rc, "Failed to get variable name");
      for (k = 0; k < numTotalVar; k++) {
        if (!strcmp(temp_string, varNames[k])) {
          free(temp_string);
          break; // found a match
        }
      }
      if (k == numTotalVar) { // didn't find match
        char** temp_char = realloc(varNames, (numTotalVar+1)*sizeof(char*));
        int *temp_comp = realloc(varNumComps, (numTotalVar+1)*sizeof(int));
        FC_AssociationType *temp_assoc = realloc(varAssocs, (numTotalVar+1)*
                                                 sizeof(FC_AssociationType));
        FC_MathType* temp_math = realloc(varMathTypes, (numTotalVar+1)*
                                         sizeof(FC_MathType));
        FC_DataType* temp_data = realloc(varDataTypes, (numTotalVar+1)*
                                         sizeof(FC_DataType));
        if (!temp_char || !temp_comp || !temp_assoc || !temp_math || !temp_data)
          fc_exitIfError(FC_MEMORY_ERROR);
        varNames = temp_char;
        varNumComps = temp_comp;
        varAssocs = temp_assoc;
        varMathTypes = temp_math;
        varDataTypes =temp_data;
        varNames[k] = temp_string;
        rc = fc_getVariableInfo(variables[j], NULL, &varNumComps[k],
                                &varAssocs[k], &varMathTypes[k],
                                &varDataTypes[k]);
        numTotalVar++;
      }
    }
    free(variables);
  }
  for (i = 0; i < numMesh; i++) {
    rc = fc_getSeqVariables(meshes[i], &numSeqVar, &numStepPerSeqVar, 
                            &seqVariables);
    fc_exitIfErrorPrintf(rc, "Failed to get seq variables");
    for (j = 0; j < numSeqVar; j++) {
      rc = fc_getVariableName(seqVariables[j][0], &temp_string);
      fc_exitIfErrorPrintf(rc, "Failed to get seq variable name");
      for (k = 0; k < numTotalSeqVar; k++) {
        if (!strcmp(temp_string, seqVarNames[k])) {
          free(temp_string);
          break; // found a match
        }
      }
      if (k == numTotalSeqVar) { // didn't find match, so add
        FC_Sequence sequence;
        int temp_ID = -1;
        char** temp_char = realloc(seqVarNames, (numTotalSeqVar+1)*sizeof(char*));
        int *temp_comp = realloc(seqVarNumComps, (numTotalSeqVar+1)*sizeof(int));
        FC_AssociationType *temp_assoc = realloc(seqVarAssocs, (numTotalSeqVar+1)*
                                                 sizeof(FC_AssociationType));
        FC_MathType* temp_math = realloc(seqVarMathTypes, (numTotalSeqVar+1)*
                                         sizeof(FC_MathType));
        FC_DataType* temp_data = realloc(seqVarDataTypes, (numTotalSeqVar+1)*
                                         sizeof(FC_DataType));
        int* temp_numStep = realloc(seqVarNumSteps, (numTotalSeqVar+1)*sizeof(int));
        int* temp_seqID = realloc(seqVarSeqIDs, (numTotalSeqVar+1)*sizeof(int));
        if (!temp_char || !temp_comp || !temp_assoc || !temp_math || 
            !temp_data || !temp_numStep || !temp_seqID)
          fc_exitIfError(FC_MEMORY_ERROR);
        seqVarNames = temp_char;
        seqVarNumComps = temp_comp;
        seqVarAssocs = temp_assoc;
        seqVarMathTypes = temp_math;
        seqVarDataTypes = temp_data;
        seqVarNames[k] = temp_string;
        seqVarNumSteps = temp_numStep;
        seqVarSeqIDs = temp_seqID;
        rc = fc_getVariableInfo(seqVariables[j][0], NULL, &seqVarNumComps[k],
                                &seqVarAssocs[k], &seqVarMathTypes[k],
                                &seqVarDataTypes[k]);
        seqVarNumSteps[k] = numStepPerSeqVar[j];
        rc = fc_getSequenceFromSeqVariable(numStepPerSeqVar[j],
                                           seqVariables[j], &sequence);
        for (m = 0; m < numSequence; m++) {
          if (FC_HANDLE_EQUIV(sequence, sequences[m])) {
            temp_ID = m;
            break;
          }
        }
        if (temp_ID < 0) 
          fc_exitIfErrorPrintf(FC_ERROR, "failed to find matching sequence");
        seqVarSeqIDs[k] = temp_ID;
        numTotalSeqVar++;
      }
    }
    free(numStepPerSeqVar);
    for (j = 0; j < numSeqVar; j++)
      free(seqVariables[j]);
    free(seqVariables);
  }
  // debug
  //printf("vars\n");
  //for (i = 0; i < numTotalVar; i++)
  //printf("%d: %s\n", i, varNames[i]);
  //printf("seq vars\n");
  //for (i = 0; i < numTotalSeqVar; i++)
  //printf("%d: %s\n", i, seqVarNames[i]);

  //------------------------------------------------------------------
  //  create case file - keep this file open throughout entire process
  //------------------------------------------------------------------
  // create file names
  temp_string = malloc((path_base_length + 6)*sizeof(char));
  case_file_name = malloc((base_length + 6)*sizeof(char));
  if (!temp_string || !case_file_name)
    fc_exitIfError(FC_MEMORY_ERROR);
  sprintf(temp_string, "%s.case", path_base_file_name);
  sprintf(case_file_name, "%s.case", base_file_name);

  // open
  case_file = fopen(temp_string, "w");
  if (!case_file)
    fc_exitIfErrorPrintf(FC_FILE_IO_ERROR, "Failed to open file '%s'", 
                         temp_string);
  free(temp_string);

  // header
  fprintf(case_file, "#\n");
  fprintf(case_file, "# Generated by fc2ensight\n");
  fprintf(case_file, "# Case File: %s\n", case_file_name);
  fprintf(case_file, "\nFORMAT\n\n");
  fprintf(case_file, "type:    ensight gold\n");
  free(case_file_name);
  fflush(NULL);

  fprintf(case_file, "\nGEOMETRY\n\n");
  fflush(NULL);

  //------------------------------------------------------------------
  //  create geom file
  //------------------------------------------------------------------
  // create file names
  temp_string = malloc((path_base_length + 5)*sizeof(char));
  geom_file_name = malloc((base_length + 5)*sizeof(char));\
  if (!temp_string || ! geom_file_name)
    fc_exitIfError(FC_MEMORY_ERROR);
  sprintf(temp_string, "%s.geo", path_base_file_name);
  sprintf(geom_file_name, "%s.geo", base_file_name);

  // --- Aside: Write geom file name to case file
  fprintf(case_file, "model:   %s\n", geom_file_name);
  fflush(NULL);
    
  // open geom file
  geom_file = fopen(temp_string, "w");
  if (!geom_file)
    fc_exitIfErrorPrintf(FC_FILE_IO_ERROR, "Failed to open file '%s'",
                         temp_string);
  free(temp_string);
  
  // geom file header
  fprintf(geom_file, "%s\n", base_file_name);
  fprintf(geom_file, "\n");
  fprintf(geom_file, "node id given\n");
  fprintf(geom_file, "element id given\n");
  fflush(NULL);

  // Write part info (put some in arrays to keep for handling vars later)
  meshElemTypes = malloc(numMesh*sizeof(FC_ElementType));
  subsetAssocsPerMesh = malloc(numMesh*sizeof(FC_AssociationType*));
  numSubsetVertsPerMesh = malloc(numMesh*sizeof(int*));
  numMemberPerMesh = malloc(numMesh*sizeof(int*));
  subsetVertIDsPerMesh = malloc(numMesh*sizeof(int**));
  memberIDsPerMesh = malloc(numMesh*sizeof(int**));
  for (i = 0; i < numMesh; i++) {
    subsetAssocsPerMesh[i] = malloc(numSubsetPerMesh[i]*sizeof(FC_AssociationType));
    numSubsetVertsPerMesh[i] = malloc(numSubsetPerMesh[i]*sizeof(int));
    numMemberPerMesh[i] = malloc(numSubsetPerMesh[i]*sizeof(int));
    subsetVertIDsPerMesh[i] = malloc(numSubsetPerMesh[i]*sizeof(int*));
    memberIDsPerMesh[i] = malloc(numSubsetPerMesh[i]*sizeof(int*));
  }
  partID = 0;
  for (i = 0; i < numMesh; i++) {
    int numDim, numVert, numElem, numVertPerElem;
    int* conns;
    double* coords;
    
    // mesh query
    rc = fc_getMeshName(meshes[i], &temp_string);
    fc_exitIfErrorPrintf(rc, "Failed to get mesh name");
    rc = fc_getMeshInfo(meshes[i], NULL, &numDim, &numVert, &numElem, 
                        &meshElemTypes[i]);
    fc_exitIfErrorPrintf(rc, "Failed to get mesh info");
    numVertPerElem = fc_getElementTypeNumVertex(meshElemTypes[i]);
    rc = fc_getMeshCoordsPtr(meshes[i], &coords);
    fc_exitIfErrorPrintf(rc, "Failed to get mesh coords");
    rc = fc_getMeshElementConnsPtr(meshes[i], &conns);
    fc_exitIfErrorPrintf(rc, "Failed to get mesh elem conns");
    
    // special cases
    // FIX! Could handle this, just have to do a lot of work
    if (meshElemTypes[i] == FC_ET_MIXED)
      fc_exitIfErrorPrintf(rc, "Cannot yet handle mixed element types");
    
    // write metadata
    fprintf(geom_file, "part\n");
    fprintf(geom_file, "%10d\n", partID+1);  // index from 1, not 0
    fprintf(geom_file, "%s\n", temp_string);
    fflush(NULL);
    free(temp_string);
    
    // write vertex ids & coords
    fprintf(geom_file, "coordinates\n");
    fprintf(geom_file, "%10d\n", numVert);
    fflush(NULL);
    for (j = 0; j < numVert; j++) 
      fprintf(geom_file, "%10d\n", j);
    // print coords as all x's, then all y's, then all z's 
    for (j = 0; j < numDim; j++) 
      for (k = 0; k < numVert; k++)
        fprintf(geom_file, "%12.5e\n", coords[k*numDim+j]);
    for (j = numDim; j < 3; j++)
      for (k = 0; k < numVert; k++)
        fprintf(geom_file, "%12.5e\n", 0.);
    fflush(NULL);

    // write element ids & conns
    fprintf(geom_file, "%s\n", getEnsightElementType(meshElemTypes[i]));
    fprintf(geom_file, "%10d\n", numElem);
    fflush(NULL);
    for (j = 0; j < numElem; j++) 
      fprintf(geom_file, "%10d\n", j);
    // have to add one because ensight indices start at 1 not 0
    for (j = 0; j < numElem; j++) {
      for (k = 0; k < numVertPerElem; k++)
        fprintf(geom_file, "%10d", conns[j*numVertPerElem+k] + 1);
      fprintf(geom_file, "\n");
    }
    fflush(NULL);

    // finished a part
    partID++;

    for (j = 0; j < numSubsetPerMesh[i]; j++) {
      int numMember, *memberIDs;
      FC_AssociationType assoc;
      int numSubsetVert, *subsetVertIDs, numVertPerSubsetElem;
      FC_ElementType subsetElemType;
      int *vertLU;

      // subset query
      rc = fc_getSubsetName(subsetsPerMesh[i][j], &temp_string);
      fc_exitIfErrorPrintf(rc, "Failed to get subset name");
      rc = fc_getSubsetAssociationType(subsetsPerMesh[i][j], &assoc);
      fc_exitIfErrorPrintf(rc, "Failed to get subset assoc");
      
      rc = fc_getSubsetMembersAsArray(subsetsPerMesh[i][j], &numMember,
                                      &memberIDs);
      fc_exitIfErrorPrintf(rc, "Failed to get subset members");
      rc = fc_getMeshEntityElementType(meshes[i], assoc, memberIDs[0],
                                       &subsetElemType);
      numVertPerSubsetElem = fc_getElementTypeNumVertex(subsetElemType);
      fc_exitIfErrorPrintf(rc, "Failed to get entity element type");
      rc = fc_changeMeshEntityType(meshes[i], assoc, numMember, memberIDs,
                                   FC_AT_VERTEX, 1, &numSubsetVert,
                                   &subsetVertIDs);
      fc_exitIfErrorPrintf(rc, "Failed to get subset vertices");
      vertLU = (int*)malloc(numVert*sizeof(int));
      if (!vertLU)
        fc_exitIfError(FC_MEMORY_ERROR);
      for (k = 0; k < numSubsetVert; k++)
        vertLU[subsetVertIDs[k]] = k;

      // special cases
      // FIX! Could handle this, just have to do a lot of work
      if (subsetElemType == FC_ET_MIXED)
        fc_exitIfErrorPrintf(rc, "Cannot yet handle mixed element types");
      if (assoc == FC_AT_WHOLE_MESH) {
        int num, *IDs;
        rc = fc_changeMeshEntityType(meshes[i], assoc, numMember, memberIDs,
                                     FC_AT_ELEMENT, 1, &num, &IDs);
        fc_exitIfErrorPrintf(rc, "Failed to get subset whole");
        free(memberIDs);
        numMember = num;
        memberIDs = IDs;
        assoc = FC_AT_ELEMENT;
      }

      // save for later
      subsetAssocsPerMesh[i][j] = assoc;
      numSubsetVertsPerMesh[i][j] = numSubsetVert;
      subsetVertIDsPerMesh[i][j] = subsetVertIDs;
      numMemberPerMesh[i][j] = numMember;
      memberIDsPerMesh[i][j] = memberIDs;
      
      // write metadata
      fprintf(geom_file, "part\n");
      fprintf(geom_file, "%10d\n", partID+1); // index from 1 not 0
      fprintf(geom_file, "%s\n", temp_string);
      fflush(NULL);
      free(temp_string);

      // write vertex ids & coords
      fprintf(geom_file, "coordinates\n");
      fprintf(geom_file, "%10d\n", numSubsetVert);
      for (k = 0; k < numSubsetVert; k++)
        fprintf(geom_file, "%10d\n", subsetVertIDs[k]); 
      // print coords as all x's, then all y's, then all z's
      for (k = 0; k < numDim; k++)
        for (n = 0; n < numSubsetVert; n++)
          fprintf(geom_file, "%12.5e\n", coords[subsetVertIDs[n]*numDim + k]);
      for (k = numDim; k < 3; k++)
        for (n = 0; n < numSubsetVert; n++)
          fprintf(geom_file, "%12.5e\n", 0.);
      fflush(NULL);

      // write element ids & conns
      fprintf(geom_file, "%s\n", getEnsightElementType(subsetElemType));
      fprintf(geom_file, "%10d\n", numMember);
      fflush(NULL);
      for (k = 0; k < numMember; k++)
        fprintf(geom_file, "%10d\n", memberIDs[k]);
      fflush(NULL);
      // have to add one because ensight indices start at 1 not 0
      if (assoc == FC_AT_VERTEX) {
        for (k = 0; k < numMember; k++)
          fprintf(geom_file, "%10d\n", k+1);
      } 
      else {
        int* tempConns;
        if (assoc == FC_AT_EDGE) {
          rc = fc_getMeshEdgeConnsPtr(meshes[i], &tempConns);
          fc_exitIfErrorPrintf(rc, "Failed to get edge conns");
        }
        else if (assoc == FC_AT_FACE) {
          rc = fc_getMeshFaceConnsPtr(meshes[i], NULL, &numVertPerSubsetElem,
                                      &tempConns);
          fc_exitIfErrorPrintf(rc, "Failed to get face conns");
        }
        else if (assoc == FC_AT_ELEMENT) {
          tempConns = conns;
        }
        for (k = 0; k < numMember; k++) {
          for (n = 0; n < numVertPerSubsetElem; n++)
            fprintf(geom_file, "%10d", 
                    vertLU[tempConns[memberIDs[k]*numVertPerSubsetElem+n]] + 1);
          fprintf(geom_file, "\n");
        }
      }
      fflush(NULL);
      
      // finished a part
      free(vertLU);
      partID++;
    }
    fc_releaseMesh(meshes[i]);
  }
  fclose(geom_file);

  // --- Aside: write to case file, beginning of variable section
  fprintf(case_file, "\nVARIABLE\n\n");
  fflush(NULL);

  //---------------------------------------------------------------------
  //  Write variable files
  //---------------------------------------------------------------------
  for (i = 0; i < numTotalVar; i++) {
    // make name w/o spaces
    temp_string = malloc((strlen(varNames[i])+1)*sizeof(char));
    strcpy(temp_string, varNames[i]);
    for (j = 0; j < (int)strlen(temp_string); j++)
      if (temp_string[j] == ' ')
        temp_string[j] = '_';

    // FIX? for now, punt if variables are on edges or face -- but we
    //        could put those on edge or face subsets if we really want to
    if (varAssocs[i] == FC_AT_EDGE || varAssocs[i] == FC_AT_FACE) {
      printf("Skipping variable '%s' because ensight can't handle assoc of "
             " type %s\n", varNames[i], 
             fc_getAssociationTypeText(varAssocs[i]));
      fflush(NULL);
    }
    else if (varDataTypes[i] != FC_DT_INT && varDataTypes[i] != FC_DT_FLOAT &&
        varDataTypes[i] != FC_DT_DOUBLE) {
      printf("Skipping variable '%s' because ensight can't handle data of "
             " type %s\n", varNames[i], fc_getDataTypeText(varDataTypes[i]));
      fflush(NULL);
    }
    else if (varMathTypes[i] != FC_MT_SCALAR && varMathTypes[i] != FC_MT_VECTOR) {
      printf("Skipping variable '%s' because this function can't yet handle "
             "data of math type %s\n", varNames[i], 
             fc_getMathTypeText(varMathTypes[i]));
      fflush(NULL);
    }
    // constants
    else if (varAssocs[i] == FC_AT_WHOLE_MESH && varNumComps[i] == 1) {
      if (doGlobals) {
	for (j = 0; j < numMesh; j++) {
	  FC_Variable *returnVariables;
	  int numReturnVariables;
	  double* data;

	  rc = fc_getMeshName(meshes[j], &temp_string2);
	  fc_exitIfErrorPrintf(rc, "Failed to get mesh name");
	  for (k = 0; k < (int)strlen(temp_string2); k++)
	    if (temp_string2[k] == ' ')
	      temp_string2[k] = '_';
	  rc = fc_getVariableByName(meshes[j], varNames[i],
				    &numReturnVariables,
				    &returnVariables);
	  fc_exitIfErrorPrintf(rc, "Failed to get variable by name");
	  if (numReturnVariables != 1){
	    fc_exitIfErrorPrintf(FC_INPUT_ERROR,"Failed to get unique variable '%s' - found %d matches",
				 varNames[i],numReturnVariables);
	  }
	  var = returnVariables[0];
	  free(returnVariables);

	  rc = fc_getVariableDataAsDataType(var, FC_DT_DOUBLE, (void**)&data);
	  fc_exitIfErrorPrintf(rc, "Failed to get variable data");

	  // --- Aside: write to case file
	  fprintf(case_file, "constant per case:   %s_on_%s ",
		  temp_string, temp_string2);
	  fprintf(case_file, "%g\n", data[0]);
	  fflush(NULL);
	  
	  free(temp_string2);
	  free(data);
	  fc_releaseVariable(var);
	}
      }
    } 
    // stuff that will require a file
    else if (varAssocs[i] == FC_AT_VERTEX || varAssocs[i] == FC_AT_ELEMENT) {
      int trunc_numComp = varNumComps[i];
      // check vector variables
      if (varMathTypes[i] == FC_MT_VECTOR && varNumComps[i] < 3) {
        printf("Warning, vector variable '%s' has less than three components, "
               "the remaining components will be padded with zeros\n",
               varNames[i]);
        fflush(NULL);
      }
      else if (varMathTypes[i] == FC_MT_VECTOR && varNumComps[i] > 3) {
        printf("Warning, vector variable '%s' has %d components, but will be "
               "truncated to just three components\n", varNames[i], 
               varNumComps[i]);
        fflush(NULL);
        trunc_numComp = 3;
      }
      
      // names
      temp_string2 = malloc((path_base_length + strlen(varNames[i]) + 2)*
                            sizeof(char));
      var_file_name = malloc((base_length + strlen(varNames[i]) + 2)*
                             sizeof(char));
      sprintf(temp_string2, "%s.%s", path_base_file_name, temp_string);
      sprintf(var_file_name, "%s.%s", base_file_name, temp_string);

      // --- Aside: Write var file info to case file
      if (varMathTypes[i] == FC_MT_SCALAR)
        fprintf(case_file, "scalar ");
      else if (varMathTypes[i] == FC_MT_VECTOR)
        fprintf(case_file, "vector ");
      if (varAssocs[i] == FC_AT_VERTEX)
        fprintf(case_file, "per node:    ");
      else if (varAssocs[i] == FC_AT_ELEMENT)
        fprintf(case_file, "per element:    ");
      fprintf(case_file, "%s  %s\n", temp_string, var_file_name);
      fflush(NULL);
      free(var_file_name);
 
      // open var file
      var_file = fopen(temp_string2, "w");
      if (!var_file)
        fc_exitIfErrorPrintf(FC_FILE_IO_ERROR, "Failed to open file '%s'",
                             temp_string);
      free(temp_string2);
      
      // var file header
      fprintf(var_file, "%s\n", temp_string);
      fflush(NULL);

      // Write part info
      partID = 0;
      for (j = 0; j < numMesh; j++) {
	FC_Variable *returnVariables;
	int numReturnVariables;
        int numDataPoint;
        double* data;
        fc_getVariableByName(meshes[j], varNames[i], 
				  &numReturnVariables,
				  &returnVariables);
        // skip this mesh & its subsets if this variable doesn't exist
	if (numReturnVariables != 1){
	  if (numReturnVariables > 1){
	    fc_printfWarningMessage("Problems finding (unique) variable '%s' - found %d matches",
				    varNames[i],numReturnVariables);
	    if (returnVariables) free(returnVariables);
	  }
          partID += 1 + numSubsetPerMesh[j];
          continue;
	}
	var = returnVariables[0];
	free(returnVariables);
        rc = fc_getVariableNumDataPoint(var, &numDataPoint);
        fc_exitIfErrorPrintf(rc, "Failed to get numDataPoint");
        rc = fc_getVariableDataAsDataType(var, FC_DT_DOUBLE, (void**)&data);
        fc_exitIfErrorPrintf(rc, "Failed to get variable data");
        // mesh part 
        fprintf(var_file, "part\n");
        fprintf(var_file, "%10d\n", partID + 1);  // index from 1
        fflush(NULL);
        if (varAssocs[i] == FC_AT_VERTEX)
          fprintf(var_file, "coordinates\n");
        else if (varAssocs[i] == FC_AT_ELEMENT)
          fprintf(var_file, "%s\n", getEnsightElementType(meshElemTypes[j]));
        for (k = 0; k < trunc_numComp; k++) 
          for (n = 0; n < numDataPoint; n++)
            fprintf(var_file, "%12.5e\n", data[n*varNumComps[i]+k]);
        if (varMathTypes[i] == FC_MT_VECTOR && varNumComps[i] < 3) {
          for (k = trunc_numComp; k < 3; k++)
            for (n = 0; n < numDataPoint; n++)
              fprintf(var_file, "%12.5e\n", 0.);
        }
        fflush(NULL);
        partID++;
        // subset parts
        for (k = 0; k < numSubsetPerMesh[j]; k++) {
          if (varAssocs[i] == FC_AT_VERTEX) {
            fprintf(var_file, "part\n");
            fprintf(var_file, "%10d\n", partID + 1);  // index from 1
            fprintf(var_file, "coordinates\n");
            for (m = 0; m < trunc_numComp; m++) 
              for (n = 0; n < numSubsetVertsPerMesh[j][k]; n++)
                fprintf(var_file, "%12.5e\n", 
                        data[subsetVertIDsPerMesh[j][k][n]*varNumComps[i]+m]);
            if (varMathTypes[i] == FC_MT_VECTOR && varNumComps[i] < 3) {
              for (m = trunc_numComp; m < 3; m++)
                for (n = 0; n < numSubsetVertsPerMesh[j][k]; n++)
                  fprintf(var_file, "%12.5e\n", 0.);
            }
          }
          else if (varAssocs[i] == FC_AT_ELEMENT &&
                   subsetAssocsPerMesh[j][k] == FC_AT_ELEMENT) {
            fprintf(var_file, "part\n");
            fprintf(var_file, "%10d\n", partID + 1);  // index from 1
            fprintf(var_file, "%s\n", getEnsightElementType(meshElemTypes[j]));
            for (m = 0; m < trunc_numComp; m++)
              for (n = 0; n < numMemberPerMesh[j][k]; n++)
                fprintf(var_file, "%12.5e\n",
                        data[memberIDsPerMesh[j][k][n]*varNumComps[i]+m]);
            if (varMathTypes[i] == FC_MT_VECTOR && varNumComps[i] < 3) {
              for (m = trunc_numComp; m < 3; m++)
                for (n = 0; n < numMemberPerMesh[j][k]; n++)
                  fprintf(var_file, "%12.5e\n", 0.);
            }
          }
          fflush(NULL);     
          partID++;
        } // end of loop over subsets
        free(data);
        fc_releaseVariable(var);
      } // end of loop over meshes
      fclose(var_file);
    } // end of stuff that requires a file
    else {
      printf("Skipping variable '%s' because ensight can't handle it:\n", 
             varNames[i]);
      printf(" numComp = %d, assoc = %s, mathType = %s, dataType = %s\n",
             varNumComps[i], fc_getAssociationTypeText(varAssocs[i]),
             fc_getMathTypeText(varMathTypes[i]),
             fc_getDataTypeText(varDataTypes[i]));
      fflush(NULL);
    }
    free(temp_string);
  } // loop over variables

  //---------------------------------------------------------------------
  //  Write seq variable files
  //---------------------------------------------------------------------
  for (i = 0; i < numTotalSeqVar; i++) {
    // make name w/o spaces
    temp_string = malloc((strlen(seqVarNames[i])+1)*sizeof(char));
    strcpy(temp_string, seqVarNames[i]);
    for (j = 0; j < (int)strlen(temp_string); j++)
      if (temp_string[j] == ' ')
        temp_string[j] = '_';

    // FIX? for now, punt if variables are on edges or face -- but we
    //        could put those on edge or face subsets if we really want to
    if (seqVarAssocs[i] == FC_AT_EDGE || seqVarAssocs[i] == FC_AT_FACE) {
      printf("Skipping seqVariable '%s' because ensight can't handle assoc of "
             " type %s\n", seqVarNames[i], 
             fc_getAssociationTypeText(seqVarAssocs[i]));
      fflush(NULL);
    }
    else if (seqVarDataTypes[i] != FC_DT_INT && seqVarDataTypes[i] != FC_DT_FLOAT &&
        seqVarDataTypes[i] != FC_DT_DOUBLE) {
      printf("Skipping seqVariable '%s' because ensight can't handle data of "
             " type %s\n", seqVarNames[i], fc_getDataTypeText(seqVarDataTypes[i]));
      fflush(NULL);
    }
    else if (seqVarMathTypes[i] != FC_MT_SCALAR && seqVarMathTypes[i] != FC_MT_VECTOR) {
      printf("Skipping seqVariable '%s' because this function can't yet handle "
             "data of math type %s\n", seqVarNames[i], 
             fc_getMathTypeText(seqVarMathTypes[i]));
      fflush(NULL);
    }
    // constants
    else if (seqVarAssocs[i] == FC_AT_WHOLE_MESH && seqVarNumComps[i] == 1) {
      if (doGlobals) {
	for (j = 0; j < numMesh; j++) {
	  FC_Variable **returnSeqVars;
	  int numReturnVars,*returnSteps;
	  double* data;
	  rc = fc_getMeshName(meshes[j], &temp_string2);
	  fc_exitIfErrorPrintf(rc, "Failed to get mesh name");
	  for (k = 0; k < (int)strlen(temp_string2); k++)
	    if (temp_string2[k] == ' ')
	      temp_string2[k] = '_';
	  rc = fc_getSeqVariableByName(meshes[j], seqVarNames[i],
				       &numReturnVars,
				       &returnSteps,&returnSeqVars);
	  fc_exitIfErrorPrintf(rc, "Failed to get seqVariable by name");
	  if (numReturnVars != 1){
	    fc_exitIfErrorPrintf(FC_INPUT_ERROR, 
				 "Failed to get (unique) seqVariable '%s' by name - found %d matches",
				 seqVarNames[i],numReturnVars);
	  }
	  numStep = returnSteps[0];
	  seqVar = returnSeqVars[0];
	  free(returnSteps);
	  free(returnSeqVars);
	  
	  // --- Aside: write to case file
	  fprintf(case_file, "constant per case:   %d  %s_on_%s ", 
		  seqVarSeqIDs[i]+1, temp_string, temp_string2);
	  fflush(NULL);
	  for (k = 0; k < numStep; k++) {
	    // put in some line breaks, ensight chokes if too long a line
	    if (k != 0 && k%5 == 0)
	      fprintf(case_file, "\n");
	    rc = fc_getVariableDataAsDataType(seqVar[k], FC_DT_DOUBLE, 
					      (void**)&data);
	    fc_exitIfErrorPrintf(rc, "Failed to get seqVariable data");
	    fprintf(case_file, " %g", data[0]);
	    fc_releaseVariable(seqVar[k]);
	    free(data);
	  }
	  fprintf(case_file, "\n");
	  fflush(NULL);
	  
	  // cleanup
	  free(temp_string2);
	  free(seqVar);
	}
      }
    } 
    // stuff that will require a file
    else if (seqVarAssocs[i] == FC_AT_VERTEX || 
             seqVarAssocs[i] == FC_AT_ELEMENT) {
      int temp_int;
      int trunc_numComp = seqVarNumComps[i];
      int number_len; // length of the step number field in name
      char* number_p;  // char pointer to start of number field in name
      char controlStr[50];
      // check vector seqVariables
      if (seqVarMathTypes[i] == FC_MT_VECTOR && seqVarNumComps[i] < 3) {
        printf("Warning, vector seqVariable '%s' has less than three components, "
               "the remaining components will be padded with zeros\n",
               seqVarNames[i]);
        fflush(NULL);
      }
      else if (seqVarMathTypes[i] == FC_MT_VECTOR && seqVarNumComps[i] > 3) {
        printf("Warning, vector seqVariable '%s' has %d components, but will be "
               "truncated to just three components\n", seqVarNames[i], 
               seqVarNumComps[i]);
        fflush(NULL);
        trunc_numComp = 3;
      }
      
      // names
      number_len = 1;
      temp_int = seqVarNumSteps[i];
      while(temp_int /= 10)
        number_len++;
      temp_string2 = malloc((path_base_length + strlen(seqVarNames[i]) + 
                             number_len + 3)*sizeof(char));
      seqVar_file_name = malloc((base_length + strlen(seqVarNames[i]) + 
                                 number_len + 3)*sizeof(char));
      sprintf(temp_string2, "%s.%s_", path_base_file_name, temp_string);
      sprintf(seqVar_file_name, "%s.%s_", base_file_name, temp_string);
      number_p = &seqVar_file_name[base_length + strlen(seqVarNames[i]) + 2];
      for (j = 0; j < number_len; j++)
        number_p[j] = '*';
      number_p[number_len] = '\0';
      number_p = &temp_string2[path_base_length + strlen(seqVarNames[i]) + 2];
      sprintf(controlStr, "%%0%dd", number_len);
      //printf("controlStr = %s\n", controlStr);

      // --- Aside: Write seqVar file info to case file
      if (seqVarMathTypes[i] == FC_MT_SCALAR)
        fprintf(case_file, "scalar ");
      else if (seqVarMathTypes[i] == FC_MT_VECTOR)
        fprintf(case_file, "vector ");
      if (seqVarAssocs[i] == FC_AT_VERTEX)
        fprintf(case_file, "per node:    ");
      else if (seqVarAssocs[i] == FC_AT_ELEMENT)
        fprintf(case_file, "per element:    ");
      fprintf(case_file, "%d  %s  %s\n", seqVarSeqIDs[i]+1, temp_string, 
              seqVar_file_name);
      fflush(NULL);
      free(seqVar_file_name);

      // loop over steps
      for (p = 0; p < seqVarNumSteps[i]; p++) {
        
        // open seqVar file
        sprintf(number_p, controlStr, p+1);
        seqVar_file = fopen(temp_string2, "w");
        if (!seqVar_file)
          fc_exitIfErrorPrintf(FC_FILE_IO_ERROR, "Failed to open file '%s'",
                               temp_string2);
      
        // seqVar file header
        fprintf(seqVar_file, "%s\n", temp_string);
        fflush(NULL);

        // Write part info
        partID = 0;
        for (j = 0; j < numMesh; j++) {
	  FC_Variable **returnVars;
	  int numReturnVars, *returnVarsSteps;
          int numDataPoint;
          double* data;
          fc_getSeqVariableByName(meshes[j], seqVarNames[i],
				       &numReturnVars,
				       &returnVarsSteps, &returnVars);
          // skip this mesh & its subsets if this seqVariable doesn't exist
	  if (numReturnVars != 1){
	    if (numReturnVars > 1){
	      fc_printfWarningMessage("Problems finding (unique) variable '%s' - found %d matches",
				      seqVarNames[i],numReturnVars);
	      for (k = 0; k < numReturnVars; k++){
		free(returnVars[k]);
	      }
	      free(returnVars);
	      free(returnVarsSteps);
	    }
            partID += 1 + numSubsetPerMesh[j];
            continue;
	  }
	  seqVar = returnVars[0];
	  numStep = returnVarsSteps[0];
	  free(returnVars);
	  free(returnVarsSteps);

          rc = fc_getVariableNumDataPoint(seqVar[p], &numDataPoint);
          fc_exitIfErrorPrintf(rc, "Failed to get numDataPoint");
          rc = fc_getVariableDataAsDataType(seqVar[p], FC_DT_DOUBLE, 
                                            (void**)&data);
          fc_exitIfErrorPrintf(rc, "Failed to get seqVariable data");
          // mesh part 
          fprintf(seqVar_file, "part\n");
          fprintf(seqVar_file, "%10d\n", partID + 1);  // index from 1
          if (seqVarAssocs[i] == FC_AT_VERTEX)
            fprintf(seqVar_file, "coordinates\n");
          else if (seqVarAssocs[i] == FC_AT_ELEMENT)
            fprintf(seqVar_file, "%s\n", getEnsightElementType(meshElemTypes[j]));
          for (k = 0; k < trunc_numComp; k++) 
            for (n = 0; n < numDataPoint; n++)
              fprintf(seqVar_file, "%12.5e\n", data[n*seqVarNumComps[i]+k]);
          if (seqVarMathTypes[i] == FC_MT_VECTOR && seqVarNumComps[i] < 3) {
            for (k = trunc_numComp; k < 3; k++)
              for (n = 0; n < numDataPoint; n++)
                fprintf(seqVar_file, "%12.5e\n", 0.);
          }
          fflush(NULL);
          partID++;
          // subset parts
          for (k = 0; k < numSubsetPerMesh[j]; k++) {
            if (seqVarAssocs[i] == FC_AT_VERTEX) {
              fprintf(seqVar_file, "part\n");
              fprintf(seqVar_file, "%10d\n", partID + 1);  // index from 1
              fprintf(seqVar_file, "coordinates\n");
              for (m = 0; m < trunc_numComp; m++) 
                for (n = 0; n < numSubsetVertsPerMesh[j][k]; n++)
                  fprintf(seqVar_file, "%12.5e\n", 
                          data[subsetVertIDsPerMesh[j][k][n]*seqVarNumComps[i]+m]);
              if (seqVarMathTypes[i] == FC_MT_VECTOR && seqVarNumComps[i] < 3) {
                for (m = trunc_numComp; m < 3; m++)
                  for (n = 0; n < numSubsetVertsPerMesh[j][k]; n++)
                    fprintf(seqVar_file, "%12.5e\n", 0.);
              }
            }
            else if (seqVarAssocs[i] == FC_AT_ELEMENT &&
                     subsetAssocsPerMesh[j][k] == FC_AT_ELEMENT) {
              fprintf(seqVar_file, "part\n");
              fprintf(seqVar_file, "%10d\n", partID + 1);  // index from 1
              fprintf(seqVar_file, "%s\n", getEnsightElementType(meshElemTypes[j]));
              for (m = 0; m < trunc_numComp; m++)
                for (n = 0; n < numMemberPerMesh[j][k]; n++)
                  fprintf(seqVar_file, "%12.5e\n",
                          data[memberIDsPerMesh[j][k][n]*seqVarNumComps[i]+m]);
              if (seqVarMathTypes[i] == FC_MT_VECTOR && seqVarNumComps[i] < 3) {
                for (m = trunc_numComp; m < 3; m++)
                  for (n = 0; n < numMemberPerMesh[j][k]; n++)
                    fprintf(seqVar_file, "%12.5e\n", 0.);
              }
            }       
            fflush(NULL);
            partID++;
          } // end of loop over subsets
          free(data);
          fc_releaseVariable(seqVar[p]);
          free(seqVar);
        } // end of loop over meshes
        fclose(seqVar_file);
      } // end of loop over steps
      free(temp_string2);
    } // end of stuff that requires a file
    else {
      printf("Skipping seqVariable '%s' because ensight can't handle it:\n", 
             seqVarNames[i]);
      printf(" numComp = %d, assoc = %s, mathType = %s, dataType = %s\n",
             seqVarNumComps[i], fc_getAssociationTypeText(seqVarAssocs[i]),
             fc_getMathTypeText(seqVarMathTypes[i]),
             fc_getDataTypeText(seqVarDataTypes[i]));
      fflush(NULL);
    }
    free(temp_string);
  } // loop over seqVariables

  //---------------------------------------------------------------------
  //  Write sequences to case file
  //---------------------------------------------------------------------
  if (numSequence > 0) {

    fprintf(case_file, "\nTIME\n\n");
    fflush(NULL);

    for (i = 0; i < numSequence; i++) {
      FC_DataType dataType;
      void* coords;

      rc = fc_getSequenceInfo(sequences[i], &numStep, &dataType);
      fc_exitIfErrorPrintf(rc, "Failed to get sequence numStep");
      if (dataType != FC_DT_INT && dataType != FC_DT_FLOAT &&
          dataType != FC_DT_DOUBLE)
        fc_exitIfErrorPrintf(rc, "Cannot handle sequende with data type "
                             "%s\n", fc_getDataTypeText(dataType));
      rc = fc_getSequenceCoordsPtr(sequences[i], &coords);
      fc_exitIfErrorPrintf(rc, "Failed to get sequence coords");
      
      fprintf(case_file, "time set:          %d\n", i + 1);
      fprintf(case_file, "number of steps:       %d\n", numStep);
      fprintf(case_file, "filename start number: 1\n");
      fprintf(case_file, "filename increment:    1\n");
      fprintf(case_file, "time values:\n");
      fflush(NULL);
      for (j = 0; j < numStep; j++) {
        switch(dataType) {
        case FC_DT_INT: 
          fprintf(case_file, "     %d\n", ((int*)coords)[j]);     break;
        case FC_DT_FLOAT:
          fprintf(case_file, "     %g\n", ((float*)coords)[j]);   break;
        case FC_DT_DOUBLE:
          fprintf(case_file, "     %g\n", ((double*)coords)[j]);  break;
        default:
          printf("Can't reach this!!!!\n");
        }
      }
      fflush(NULL);
    }
  }

  //---------------------------------------------------------------------
  //  Finish the case file
  //---------------------------------------------------------------------
  fclose(case_file);

  //---------------------------------------------------------------------
  //  Finalize FCLib
  //---------------------------------------------------------------------
  rc = fc_finalLibrary();
  fc_exitIfErrorPrintf(rc, "Failed to finalize library");

  //------------------------------------------------------------------
  //   cleanup & exit
  //------------------------------------------------------------------
  free(path_base_file_name);
  free(sequences);
  free(meshes);
  free(numSubsetPerMesh);
  for (i = 0; i < numMesh; i++)
    free(subsetsPerMesh[i]);
  free(subsetsPerMesh);
  exit(0);
}
