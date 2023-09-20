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
 * \file subsetter.c
 * \brief given a var (elem or nodal) name and some condition, 
 *        writes out a subsets (or promoted subsets) satisfying that
 *        condition. This is intended to be used with seq var, where the
 *        subsets are promoted to meshes and the reassmebler is used to
 *        reconstitute the data onto the orig big mesh
 *
 * \description
 *    subsetter dataset "ouputtype" "variable name" "op" threshold_value
 *    doStrict 
 *
 *    outputtype is by default subsets with option of -M to promote the
 *    segmented threshold subset to new meshes. In either case, the original
 *    mesh will be written out as well. Note that in the subset to mesh
 *    promotion the numbering will change, but there will be a elem and vert
 *    numbering map variable associated with the new mesh.
 *
 *    The same val of doStrict is used
 *    for subset to mesh promotion in the case of vertex vars, but we
 *    may not want this.
 *
 *    "op: is the threshold operator such as ">", "<=",
 *    "==", etc. Quotes are required around the op string, but can be
 *    dropped around the variable name if there are no spaces or 
 *    punctuation in the name. "doStrict" is used in the call to
 *    fc_copySubsetWithNewAssocation and more info on it can be found
 *    there.
 *
 *    Multicomponent vars will be discovered and merged and their
 *    magnitude is used for the thresh comparison.
 *
 *    New subsets or meshes are named <meshname>_T_s# in the case
 *    of variables and <meshname>_T#_s# where # is the timestep
 *    number in the case of sequence variables and the s# is the segment
 *    number
 *
 *    var can be a normal var, but in the case of a seq var, this
 *    can be paired with region reassembler for reconsitituting
 *
 *    output files are:
 *    if dont promote: subsetter_subsets.ex2 - has mesh and subsets on it (no vars)
 *    if promote: 1) subsetter_geomonly.ex2 - for the main mesh with geom only
 *                2) subsetter_subsets.ex2  - for the promoted meshes (with all seq vars
 *                                            at that step and with elem/vertex mapping
 *                                            vars for reconsitituting)
 *
 *
 *
 * \todo
 *   - check doStrict
 *
 * \modifications
 *   - 12/02/07 ACG 
 *   - 03/28/08 ACG - picking up again now that nodal subset writeout is fixed. Note - not using
 *                    seq subset since segmenting the subset - this way I can set the names. 
 *   - 04/12/08 ACG - rearranging a bit. 
 *   - 06/XX/08 ACG - various fixes becuase of filesize issues
 *   - 07/13/08 ACG - put sequences on geom only for reconstituting later
 *    -08/10/08 ACG - cleaning up for pairing with reassembler
 */

#include <string.h>
#include "exodusII.h"
#include "fc.h"

static FC_Dataset ds;
static FC_Dataset new_ds;
static FC_Dataset geomonly_ds;

static char* op;
static double threshold;
static int doStrict;
static int doPromote = 0; //promote subsets to meshes
static int threshType = -1;
static FC_VerbosityLevel verbose_level = FC_QUIET;
static char* file_name = NULL;
static char* var_name = NULL;
static char* ffilebase = NULL;

static int doFeature = 0;
static int* numROIPerStep;
static FC_Subset** ROIsPerStep;
static FC_FeatureGroup *featureGroup;

static void getArgs(int, char**);
static int createThresholdableVar(FC_Variable, FC_Variable*);
static int processOneVar(FC_Mesh mesh, FC_Variable var, int step, int currnSubsets);

enum {ABSOLUTE=0, SCALEMIN, SCALEMAX};

int main(int argc, char** argv) {
  FC_ReturnCode rc;
  int i;
  int numMesh;
  FC_Mesh* meshes;

  getArgs(argc, argv);

  // init library and load dataset
  rc = fc_setLibraryVerbosity(verbose_level);
  fc_exitIfError(rc);
  rc = fc_initLibrary();
  fc_exitIfError(rc);
  rc = fc_loadDataset(file_name, &ds);
  fc_exitIfErrorPrintf(rc,"Can't load dataset\n");

  //set up new dataset -- is this really the only way I can
  //get a dataset I can drop things from and write?
  rc = fc_createDataset("copied_ds", &new_ds); //will have to handle case where there are no subsets created
  fc_exitIfErrorPrintf(rc, "Failed to create new dataset");
  rc = fc_createDataset("geomonly_ds", &geomonly_ds); //this will have only the orig geom. have to
  //write to a separte ds so that it doesnt gain vertex vars in the eventual exodus writeout
  fc_exitIfErrorPrintf(rc, "Failed to create new dataset");

  // get number of meshes and setup var arrays
  rc = fc_getMeshes(ds, &numMesh, &meshes);

  printf("Subset summary for dataset '%s'\n",file_name);
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
      FC_Variable var_handle;
      char* mesh_name;
      fc_getMeshName(meshes[i], &mesh_name);
      printf("Results for mesh '%s':\n",mesh_name);
      fflush(NULL);

      //check if basic var
      rc = fc_getOrGenerateUniqueVariableByName(meshes[i], var_name,
						&var_handle);
      fc_exitIfError(rc);
      if (FC_HANDLE_EQUIV(var_handle, FC_NULL_VARIABLE)){
	  fc_printfWarningMessage("Problems finding unique basic variable '%s'",
				  var_name);
      }
      
      if (fc_isVariableValid(var_handle)){
	//make sure its the right type
	FC_AssociationType assoc;

	rc = fc_getVariableAssociationType(var_handle, &assoc);
	fc_exitIfError(rc);
	
	if (assoc == FC_AT_ELEMENT || assoc == FC_AT_VERTEX ){
	  processOneVar( meshes[i], var_handle, -1, 0);
	}
      } else {
        //check if seq var
        FC_Variable* seq_var_handles;
        int temp_numStep;
	int j;
	int nSubsets = 0;
	rc = fc_getOrGenerateUniqueSeqVariableByName(meshes[i], var_name,
						     &temp_numStep,
						     &seq_var_handles);
        fc_exitIfError(rc);
	if (!seq_var_handles){
	  fc_printfErrorMessage("Problems finding unique seq variable '%s'",
				var_name);
	}

        if (fc_isSeqVariableValid(temp_numStep,seq_var_handles)){
	  //make sure its the right type
	  FC_AssociationType assoc;
	  
	  rc = fc_getVariableAssociationType(seq_var_handles[0], &assoc);
	  fc_exitIfError(rc);

	  if (assoc == FC_AT_ELEMENT || assoc == FC_AT_VERTEX ){ 
	    if (doFeature){
	      numROIPerStep = (int*)malloc(temp_numStep*sizeof(int));
	      ROIsPerStep = (FC_Subset**)malloc(temp_numStep*sizeof(FC_Subset*));
	    }

	    //do step by step
	    for (j=0; j < temp_numStep; j++){
	      int crSubsets;
	      crSubsets =  processOneVar(meshes[i], seq_var_handles[j], 
					 j, nSubsets);
	      nSubsets+=crSubsets;
	    } //steps
	    if (doFeature){
	      char featurefilename[100];
	      int sone = strlen(ffilebase);
	      if (sone > 20) sone = 20;
	      int stwo = strlen(mesh_name);
	      if (stwo > 20) stwo = 20;
	      int sthree = strlen(var_name);
	      if (sthree > 20) sthree = 20;
	      char formatstring[20];
	      snprintf(formatstring,20,"%s%d%s%d%s%d%s","%",sone,"s_%",stwo,"s_%",sthree,"s");
	      snprintf(featurefilename, 100, formatstring,ffilebase,mesh_name,var_name);
	      featurefilename[100-1] = '\0';

	      rc = fc_createFeatureGroup(&featureGroup);
	      fc_exitIfErrorPrintf(rc, "Can't create feature group");
	      rc = fc_trackAllSteps(temp_numStep, numROIPerStep,
				    ROIsPerStep, featureGroup);
	      fc_exitIfErrorPrintf(rc, "Can't track steps");
	      //now can clean up
	      free(numROIPerStep);
	      for (j=0; j < temp_numStep; j++){
		free(ROIsPerStep[j]);
	      }
	      free(ROIsPerStep);
	      fc_printFeatureGroup(featureGroup);
	      fc_writeFeatureGraph(featureGroup,featurefilename);
	      fc_freeFeatureGroup(featureGroup);
	    }
	  } //assoc
        } else {
          printf("Cannot get variable '%s'\n", var_name);
          fflush(NULL);
          continue;
	}
        free(seq_var_handles);
      }
      free(mesh_name);
    } //nummesh
  }
  free(meshes);

  {
    //if there are anymeshes on the geom only ds, then copy over the
    //sequences as well for reconstituing later
    
    FC_Sequence *seqs,newseq;
    int numseq;

    rc = fc_getSequences(ds, &numseq, &seqs);
    if (rc != FC_SUCCESS){
      printf("Cannot get sequences from orig ds\n");
    } else {
      for (i = 0; i < numseq; i++){
	rc = fc_copySequence(seqs[i],geomonly_ds,NULL,&newseq);
	if (rc != FC_SUCCESS){
	  printf("Cannot get copy sequence from orig ds\n");
	}
      }
      free(seqs);
    }
  }

  /*
  printf("-----------------------------\n");
  fc_printDataset(ds,"test - orig ds");
  rc = fc_getMeshes(ds, &numMesh, &meshes);
  for (i = 0; i < numMesh; i++){
    fc_printMesh(meshes[i],"",0,0,0,0);
  }
  if (meshes) free(meshes);
  fc_printDataset(new_ds,"test - new ds");
  rc = fc_getMeshes(new_ds, &numMesh, &meshes);
  for (i = 0; i < numMesh; i++){
    fc_printMesh(meshes[i],"",0,0,0,0);
  }
  if (meshes) free(meshes);
  fc_printDataset(geomonly_ds,"test - geom only ds");
  rc = fc_getMeshes(geomonly_ds, &numMesh, &meshes);
  for (i = 0; i < numMesh; i++){
    fc_printMesh(meshes[i],"",0,0,0,0);
  }
  if (meshes) free(meshes);
  printf("-----------------------------\n");
  */

  if (doPromote){
    rc = fc_writeDataset(geomonly_ds, "subsetter_geomonly.ex2", FC_FT_EXODUS);
  }
  rc = fc_writeDataset(new_ds, "subsetter_subsets.ex2", FC_FT_EXODUS);

  fc_exitIfError(rc);
  rc = fc_deleteDataset(ds);
  fc_exitIfError(rc);

  fc_finalLibrary();
  exit(0);
}



int createThresholdableVar(FC_Variable var_handle, FC_Variable* thresh_var_handle){
  int numComp;
  FC_ReturnCode rc = fc_getVariableNumComponent(var_handle,&numComp);
  fc_exitIfError(rc);
  if (numComp != 1){
    // use the magnitude instead
    rc = fc_createMagnitudeVariable(var_handle, "mag", thresh_var_handle);
    fc_exitIfErrorPrintf(rc, "Cannot make magnitude var\n");
  } else {
    *thresh_var_handle = var_handle;
  }
  return 1;
}


static int processOneVar(FC_Mesh mesh, 
			 FC_Variable var_handle, 
			 int step, int currnSubsets){
  //if its not a seqvar step = -1
  //currnSubsets is number of segments currently written out for this mesh (if this
  //the first time one will be written out, we also write out the mesh to the new
  //dataset
  FC_ReturnCode rc;
  FC_Variable thresh_var;
  char* mesh_name;

  if (createThresholdableVar(var_handle,&thresh_var)){
    int numMember;
    FC_Subset subs_handle;
    char sub_name[MAX_STR_LENGTH];
    double currthreshold;

    fc_getMeshName(mesh, &mesh_name);
    if (step == -1 ){
      snprintf(sub_name, MAX_STR_LENGTH, "%s_T",mesh_name);
    } else {
      snprintf(sub_name, MAX_STR_LENGTH, "%s_T%d",mesh_name,step);
    }
    sub_name[MAX_STR_LENGTH-1] = '\0';

    switch (threshType){
    case ABSOLUTE:
      currthreshold = threshold;
      break;
    default:
      {
	double minv, maxv;
	rc = fc_getVariableMinMax(var_handle, &minv, NULL, &maxv, NULL);
	fc_exitIfErrorPrintf(rc, "Cannot get min/max\n");
	currthreshold  = (threshType ==SCALEMIN)? minv: maxv;
	currthreshold *= threshold; //threshold is a scale val
	//	printf("max = %f min = %f scale = %f abs = %f\n", 
	//	       maxv, minv, threshold, currthreshold);
	break;
      }
    }

    rc = fc_createThresholdSubset(thresh_var, op, currthreshold, sub_name,
				  &subs_handle); //gives me all the ids
    fc_exitIfError(rc);
    
    rc = fc_getSubsetNumMember(subs_handle, &numMember);
    fc_exitIfError(rc);
    if (numMember < 1){
      if (step != -1 && doFeature){
	//for feature tracking - can only do if seqvar
	numROIPerStep[step] = 0;
	ROIsPerStep[step] = NULL;
      }
      rc = fc_deleteSubset(subs_handle);
      fc_exitIfError(rc);
      free(mesh_name);
      return 0;
    } else {
      //segment
      int numSegments;
      FC_Subset *segments;
      int num_dim;
      int i,j;

      rc = fc_segment(subs_handle, 0, &numSegments, &segments); //names work out ok
      fc_exitIfError(rc);
      //no longer need the thresh subset
      rc = fc_deleteSubset(subs_handle);
      fc_exitIfError(rc);

      if (step != -1){
	printf("Step %d:\n",step);
      }
      //write it out to the new dataset
      if (currnSubsets == 0){
	//will want this geom only later
	FC_Mesh copiedMesh;
	if (doPromote){ 
	  //copy it to a diff file than the promotions
	  //write out the whole mesh geom 
	  rc = fc_copyMesh(mesh, geomonly_ds, mesh_name, 0,0,0,0, &copiedMesh);
	  fc_exitIfErrorPrintf(rc, "Can't copy over original mesh geom");  
	} else {
	  //copy it to the same file as the non-promoted subsets since the subsets will be on this mesh
	  rc = fc_copyMesh(mesh, new_ds, mesh_name, 0,0,0,0, &copiedMesh);
	  fc_exitIfErrorPrintf(rc, "Can't copy over original mesh geom");  
	}
      }

      rc = fc_getMeshDim(mesh, &num_dim);
      fc_exitIfErrorPrintf(rc,"Can't get mesh dim");
      for (i = 0; i < numSegments; i++){
	char *segname;

	rc = fc_getSubsetName(segments[i], &segname);
	fc_exitIfErrorPrintf(rc, "Can't get subset name");
	fc_printSubset(segments[i],"thresh subset",0);

	if (doPromote){
	  FC_Mesh promotedMesh;
	  FC_Variable vertmap, elemmap;
	  //copy over the vars and also the seq vars, but for this timestep only
	  //first get the vars
	  rc = fc_createSubsetMesh(segments[i], new_ds, doStrict, 1,0,0,0,
				   segname, &promotedMesh, &vertmap, &elemmap);
	  fc_exitIfErrorPrintf(rc, "cant promote subset");		
	  //depending on doStrict, this might actually be a null mesh
	  //note that we will still have a segment, since that is a vertex subset
	  if (FC_HANDLE_EQUIV(promotedMesh, FC_NULL_MESH)){
	    fc_deleteMesh(promotedMesh);
	  }
	  //now get the seq vars for this timestep only
	  if (step > 0){
	    int numSeqVar;
	    int *numStepPerSeq;
	    FC_Variable** seqVars;
	    FC_Variable junkVar;
	    int inneri;
	    rc = fc_getSeqVariables(mesh, &numSeqVar, &numStepPerSeq, &seqVars);
	    fc_exitIfErrorPrintf(rc, "cant get seqvars on orig mesh");		
	    for (inneri = 0; inneri < numSeqVar; inneri++){
	      FC_AssociationType assoc;
	      rc = fc_getVariableAssociationType(seqVars[inneri][0], &assoc);
	      fc_exitIfErrorPrintf(rc, "cant get var assoc type");		
              switch(assoc){
	      case FC_AT_VERTEX:
		rc = fc_copyVariableToRegionMesh(seqVars[inneri][step], promotedMesh,
						 vertmap, NULL, NULL, &junkVar) ;
		break;
	      case FC_AT_ELEMENT:
		rc = fc_copyVariableToRegionMesh(seqVars[inneri][step], promotedMesh,
						 elemmap, NULL, NULL, &junkVar) ;
		break;
	      default:
		fc_printfErrorMessage("Can't copy var to region var for assoc type");
		break;
	      }
	    }
	    for (j = 0; j < numSeqVar; j++)
	      free(seqVars[j]);
	    if (seqVars) free(seqVars);
	    if (numStepPerSeq) free (numStepPerSeq);
          }
	  fc_printMesh(promotedMesh,"promoted mesh",0,0,0,0);
	} else {
	  int numCopiedMeshes;
	  FC_Mesh* copiedMeshes;
	  FC_Subset copiedSubset;

	  rc = fc_getMeshByName(new_ds, mesh_name, &numCopiedMeshes,
				&copiedMeshes);
	  fc_exitIfError(rc);
	  if (numCopiedMeshes != 1){
	    fc_exitIfErrorPrintf(FC_ERROR,
				 "Cant find unique new mesh on which to put this subset");
	  }
	  rc = fc_copySubset(segments[i], copiedMeshes[0], segname, &copiedSubset);
	  fc_exitIfErrorPrintf(rc, "cant copy subset");
	  free(copiedMeshes);
	}
	free(segname);
      } //for

      //feature tracking
      if (step == -1 || !doFeature) {
	free(segments);
      } else {
	//for feature tracking
	numROIPerStep[step] = numSegments;
	ROIsPerStep[step] = segments;
      }
      free(mesh_name);
      return numSegments;
    }

  }
  return 0;
}


void getArgs(int argc, char** argv){
  // handle arguments
  int i;

  if (argc < 5) {
  usage:
    printf("usage: %s [options] dataset \"var name\" \"op\" thresh_value/perc doStrict\n", argv[0]);
    printf("options: \n");
    printf("   -h           : print this help message\n");
    printf("   -v           : verbose: print warning and error messages\n");
    printf("   -V           : very verbose: print log and error messages\n");
    printf("   -M           : promote thresholded subsets to meshes (for reassembler later)\n");
    printf("   -f ffile     : turn on feature tracking. will write per-mesh output files with ffilename base\n");
    printf("   -a/smax/smin (required) : is the thresh value absolute (a) or a scale (times) of the max or min per timestep\n");
    printf("\n");
    printf("Prints to stdout the subsets in the dataset for which\n");
    printf("the given nodal or elem variable satisfies the threshold operation.\n");
    printf("\"op\" is an operator such as \">\", \"<=\", \"==\", etc.\n");
    printf("Quotes are required around the op string, but can be dropped\n");
    printf("around the variable name if it contains no spaces or\n");
    printf("punctuation. \"doStrict\" is used in the internal call to\n");
    printf("fc_copySubsetWithNewAssociation and more info can be found there.\n");
    printf("\n");
    printf("Example: %s -M -f \"./junk/sub\" -a ../data/gen_multivar_seq.ex2 \"temperature per vertex\" \">\" 300 1\n", argv[0]);
    printf("\n");
    printf("\n");
    fflush(NULL);
    exit(-1);
  }

  //defaults
  doStrict = -1;

  for (i = 1; i < argc; i++) {
    if (!strcmp(argv[i], "-v")) {
      verbose_level = FC_WARNING_MESSAGES;
    }
    else if (!strcmp(argv[i], "-V")) {
      verbose_level = FC_LOG_MESSAGES;
    }
    else if (!strcmp(argv[i], "-M")){
      doPromote = 1;
    }
    else if (!strcmp(argv[i], "-f")){
      doFeature = 1;
      if (i+1 > argc)
	goto usage;
      ffilebase = argv[++i];
    }
    else if (!strcmp(argv[i], "-a") || !(strcmp(argv[i], "-smin")) || !(strcmp(argv[i], "-smax"))){
      //putting these together for clarity
      if (threshType != -1){
	goto usage;
      }
      if (!strcmp(argv[i], "-a")){
	threshType = ABSOLUTE;
      } else if (!strcmp(argv[i], "-smin")){
	threshType = SCALEMIN;
      } else {
	threshType = SCALEMAX;
      }
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
  if (!file_name || !var_name || (doStrict == -1) || threshType == -1)
    goto usage;
}
