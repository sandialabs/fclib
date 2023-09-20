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
 * \file volumeSegments.c
 * \brief calculate volumes of segments of dataset based on thresholding
 *        a variable.
 *
 * $Source: /home/Repositories/fcdmf/fclib/tools/volumeSegments.c,v $
 * $Revision: 1.15 $
 * $Date: 2006/11/07 23:49:03 $
 *
 * \description
 *    Usage: volumeSegments dataset "variable name" "op" threshold_value
 *    doStrict 
 *
 *    "op: is the threshold operator such as ">", "<=",
 *    "==", etc. Quotes are required around the op string, but can be
 *    dropped around the variable name if there are no spaces or 
 *    punctuation in the name. "doStrict" is used in the call to
 *    fc_copySubsetWithNewAssocation and more info on it can be found
 *    there.
 *
 *    Example: volumeSegments data.ex2 "cool factor" ">" 0.5 1
 *
 *    Prints to stdout the volumes of segments in the dataset
 *    for which the given variable is above the given threshold value.
 *    The format of the output is one line per segment (per
 *    step, if its a sequence variable) with form
 *
 *    Volume of segmentname = vol
 *
 *    where, segmentname is the name of the segment
 *    and vol is the volume of the segment, not of the bounding box
 *    of the segment. Steps are designated by step number, if appropriate.
 *    Total volume covered and percent volume covered are also given.
 *    Outputs error message if mesh dim < 3D.
 *   
 *
 * \modifications
 *   - 10/04/04 AG based on threshBoundBox & varStats
 *   - 11/22/04 AG revised for updated fclib
 *   - 11/23/04 AG renamed to volumeSegments and created new volumeSubsets.
 *   - 11/30/04 AG replaced with new getSubsetVolume and added doStrict flag
 *   - 06/08/2006 WSD add option to do with displacements
 */

#include <string.h>
#include "fc.h"

static int isThresholdable(FC_Variable);
static int calcSegmentsVolumes(FC_Subset,int, FC_Variable* point_to_displ_var);

int main(int argc, char** argv) {
  FC_ReturnCode rc;
  int i;
  int numMesh;
  int doStrict;
  FC_Dataset ds;
  FC_Mesh* meshes;
  double threshold;
  char* file_name = NULL;
  char* mesh_name;
  char* var_name;
  char* displ_name = NULL;
  char* op;
  FC_VerbosityLevel verbose_level = FC_QUIET;
  int doDispl = 0;


  // handle arguments
  if (argc < 5) {
  usage:
    printf("usage: %s [options] dataset \"var name\" \"op\" thresh_value doStrict\n", argv[0]);
    printf("options: \n");
    printf("   -d var_name  : calc the displaced volumes using this displacement\n");
    printf("                  variable\n");
    printf("   -h           : print this help message\n");
    printf("   -v           : verbose: print warning and error messages\n");
    printf("   -V           : very verbose: print log and error messages\n");
    printf("\n");
    printf("Prints to stdout the volumes of segments in the dataset\n");
    printf("for which the given variable satisfies the threshold operation.\n");
    printf("\"op\" is an operator such as \">\", \"<=\", \"==\", etc.\n");
    printf("Quotes are required around the op string, but can be dropped\n");
    printf("around the variable name if it contains no spaces or\n");
    printf("punctuation. \"doStrict\" is used in the internal call to\n");
    printf("fc_copySubsetWithNewAssociation and more info can be found there.\n");
    printf("\n");
    printf("Example: %s data.ex2 \"cool factor\" \">\" 0.5 1\n", argv[0]);
    printf("\n");
    printf("The format of the output is one line per segment\n");
    printf("(per step if it is a sequence variable) with form:\n");
    printf("\tVolume of segment segmentname = vol\n");
    printf("where segmentname is the name of the segment,\n");
    printf("and vol is the volume of the segment, not of the bounding\n");
    printf("box of the segment. Steps are indicated by step index, if\n");
    printf("appropriate. Total volume covered and percent volume are also\n");
    printf("given.\n");
    printf("\n");
    printf("Outputs error message if dimension < 3D.\n");
    fflush(NULL);
    exit(-1);
  }

  
  for (i = 1; i < argc; i++) {
    if (!strcmp(argv[i], "-d")) {
      doDispl = 1;
      i++;
      displ_name = argv[i];
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
      if (i+4 >= argc)
        goto usage;
      file_name = argv[i];
      var_name = argv[i+1];
      op = argv[i+2];
      threshold = atof(argv[i+3]);
      doStrict = atoi(argv[i+4]);
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
  rc = fc_loadDataset(file_name, &ds);
  fc_exitIfErrorPrintf(rc,"Can't load dataset\n");

  // get number of meshes and setup var arrays
  rc = fc_getMeshes(ds, &numMesh, &meshes);


  printf("Subset volume summary for dataset '%s'\n",file_name);
  fflush(NULL);
  
  if (numMesh < 1){
    printf("No meshes were found\n");
    fflush(NULL);
    exit(0);
  } else {
    printf("%d meshes:\n", numMesh);
    fflush(NULL);
    
    // loop over meshes
    for (i = 0; i < numMesh; i++) {
      int meshTopoDim = 0, numStep;
      FC_Variable var_handle;
      FC_Variable displ, *seqDispl = NULL;

      fc_getMeshName(meshes[i], &mesh_name);
      //see if want to do areas at some point
      rc = fc_getMeshTopodim(meshes[i], &meshTopoDim);
      if (meshTopoDim < 3){
        printf("Can't calculate volume on mesh '%s' dimensionality too small = %d\n",
               mesh_name, meshTopoDim);
        fflush(NULL);
        free(mesh_name);
        continue;
      }

      if (doDispl) {
	rc = fc_getOrGenerateUniqueVariableByName(meshes[i], displ_name,
						  &displ);
	fc_exitIfError(rc);	
	if (!FC_HANDLE_EQUIV(displ, FC_NULL_VARIABLE) && 
	    !fc_isValidDisplacementVariable(meshes[i], displ)) {
	  fc_exitIfErrorPrintf(FC_ERROR, "'%s' is not a valid displ variable",
			       displ_name);
	}

	rc = fc_getOrGenerateUniqueSeqVariableByName(meshes[i], displ_name,
						     &numStep, &seqDispl);
	fc_exitIfError(rc);
	if (seqDispl && !fc_isValidDisplacementVariable(meshes[i],seqDispl[0])) {
	  fc_exitIfErrorPrintf(FC_ERROR, "'%s' is not a valid displ variable",
			       displ_name);
	}

	if (FC_HANDLE_EQUIV(displ, FC_NULL_VARIABLE) && seqDispl == NULL) {
	  fc_exitIfErrorPrintf(FC_ERROR, "Could not find displ variable '%s'", 
			       displ_name);
	}
      }
      
      printf("Results for mesh '%s':\n",mesh_name);
      fflush(NULL);
      //check if basic var
      rc = fc_getOrGenerateUniqueVariableByName(meshes[i], var_name,
						&var_handle);
      fc_exitIfError(rc);
      if (FC_HANDLE_EQUIV(var_handle, FC_NULL_VARIABLE)){
	  fc_printfWarningMessage("Problems finding unique variable '%s'",
				  var_name);
      }
      
      if (fc_isVariableValid(var_handle)){
        if (isThresholdable(var_handle)){
          //      printf("Subset on mesh %s:\n",mesh_name);
          char* sub_name = "temp_subset";
          FC_Subset subs_handle;
          
          rc = fc_createThresholdSubset(var_handle, op, threshold, sub_name,
                                        &subs_handle); //gives me all the ids
          fc_exitIfError(rc);
          
          //now get contiguous segments
	  if (doDispl) {
	    if (seqDispl) {
	      int j;
	      for (j = 0; j < numStep; j++) {
		printf("Step %d:\n", j);
		rc = calcSegmentsVolumes(subs_handle,doStrict,&seqDispl[j]);
		fc_exitIfError(rc);
	      }
	    }
	    else {
	      rc = calcSegmentsVolumes(subs_handle,doStrict,&displ);
	      fc_exitIfError(rc);
	    }
	  }
	  else {
	    rc = calcSegmentsVolumes(subs_handle,doStrict,NULL);
	    fc_exitIfError(rc);
	  }
        }
      } else {
        //check if seq var
        FC_Variable* seq_var_handles;
        int temp_numStep;
	int j;

	rc = fc_getOrGenerateUniqueSeqVariableByName(meshes[i], var_name,
						     &temp_numStep,
						     &seq_var_handles);
        fc_exitIfError(rc);
	if (!seq_var_handles){
	  fc_printfErrorMessage("Problems finding unique variable '%s'",
				var_name);
	}

        if (fc_isSeqVariableValid(temp_numStep,seq_var_handles)){
          if(isThresholdable(seq_var_handles[0])){
            //do step by step
            for (j=0; j < temp_numStep; j++){
              char* sub_name = "temp_subset";
              FC_Subset subs_handle;
              
              printf("Step %d:\n",j);
              fflush(NULL);
              rc = fc_createThresholdSubset(seq_var_handles[j], op, threshold,
                                            sub_name, &subs_handle); //gives me all the ids
              fc_exitIfError(rc);
              
              //now get contiguous segments
	      if (doDispl) {
		if (temp_numStep != numStep)
		  fc_exitIfErrorPrintf(FC_ERROR, 
				       "the var and displ var have different numStep");
		if (seqDispl) 
		  rc = calcSegmentsVolumes(subs_handle, doStrict, &seqDispl[j]);
		else 
		  rc = calcSegmentsVolumes(subs_handle, doStrict, &displ);
	      }
	      else 
		rc = calcSegmentsVolumes(subs_handle, doStrict, NULL);
	      fc_exitIfError(rc);
            }
          }
        } else {
          printf("Cannot get variable '%s'\n", var_name);
          fflush(NULL);
          continue;
	}
        free(seq_var_handles);
      }

      fc_deleteMesh(meshes[i]);
      if (seqDispl) free(seqDispl);
      free(mesh_name);
    }
  }
  free(meshes);
  fc_finalLibrary();
  exit(0);
}

int calcSegmentsVolumes(FC_Subset subs_handle, int doStrict, FC_Variable* displ_p){

  int i;
  FC_Subset* segments;
  FC_Mesh mesh; 
  int numSegments;
  double coveredvol = 0.0;
  double meshvol = -1.0;

  //get segments
  FC_ReturnCode rc = fc_segment(subs_handle, 0, &numSegments, &segments);
  fc_exitIfError(rc);

  for (i = 0; i < numSegments; i++){
    double vol = 0.0;
    char* seg_name;
    rc = fc_getSubsetName(segments[i],&seg_name);
    fc_exitIfError(rc);
    if (displ_p) 
      rc = fc_getDisplacedSubsetVolume(segments[i], doStrict, *displ_p,  &vol);
    else
      rc = fc_getSubsetVolume(segments[i],doStrict,&vol);
    fc_exitIfError(rc);
    printf("\tVolume of segment '%s' = %f\n",seg_name,vol);
    fflush(NULL);
    coveredvol+=vol;
    free(seg_name);
    fc_deleteSubset(segments[i]);
  }
  free(segments);

  rc = fc_getMeshFromSubset(subs_handle,&mesh);
  fc_exitIfError(rc);
  if (displ_p)
    rc = fc_getDisplacedMeshVolume(mesh, *displ_p, &meshvol);
  else
    rc = fc_getMeshVolume(mesh,&meshvol);
  fc_exitIfError(rc);
  printf("\tPercent coverage = %f\n",(coveredvol/meshvol));
  fflush(NULL);

  return FC_SUCCESS;
}

int isThresholdable(FC_Variable var_handle){
  int numComp;
  FC_ReturnCode rc = fc_getVariableNumComponent(var_handle,&numComp);
  fc_exitIfError(rc);
  if (numComp != 1){
    char* var_name;
    rc = fc_getVariableName(var_handle,&var_name);
    fc_exitIfError(rc);
    fc_exitIfErrorPrintf(FC_ERROR, "Cannot threshold variable '%s' because it has more than one component\n", var_name);
  }
  return 1;
}

