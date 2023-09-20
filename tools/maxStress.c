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
 * \file maxStress.c
 * \brief Report regions of maximum stress, either in compression or tension.
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/tools/maxStress.c,v $
 * $Revision: 1.7 $ 
 * $Date: 2007/04/04 01:07:14 $
 *
 * \description
 *    Usage: maxStress [options] dataset [mesh_name]
 *
 *  Report the minimum and maximum Von Mises stress in regions in compression
 *  and in regions of tension. Also those max and min normalized to the
 *  max stress over time on a per mesh basis. Regions in compression are regions
 *  where the pressure is less than 0 and regions in tension are regions where the
 *  pressure is greater than 0. The Von Mises stress and the pressure are computed
 *  from the symmetric stress tensor.
 *
 *  WARNING: Right now assumes that the stress tensor (symmetric) will
 *  be read in as separate one-component variables named: sigma-xx,
 *  sigma-yy, sigma-zz, sigma-xy, sigma-yz, sigma-zx.
 *
 * \modifications
 *   - 08/08/2006 WSD Created.
 *   - 03/13/2007 ACG not a modification, but the innards wrt creation of the stress and pressure
 *                variables have been duplicated into varNormalize. If you make a change in either
 *                file please cross reference.
 *   - 04/02/2007 ACG normalization added
 */

#include <string.h>
#include <math.h>
#include "fc.h"

// the endings of the sigma matrix in the order we plan to process them
static char sigmaPostfix[6][3] = { "xx", "yy", "zz", "xy", "yz", "zx" };

int main(int argc, char** argv) {
  FC_ReturnCode rc;
  int i, j, k;
  int numMesh = 0, numStep;
  FC_Dataset ds;
  FC_Sequence sequence;
  FC_Mesh *meshes;
  char* file_name = NULL;
  char** mesh_names = NULL;
  FC_VerbosityLevel verbose_level = FC_QUIET;
  int firstt = 1, firstc = 1;
  double max_cmax, max_tmax, min_cmin, min_tmin;
  int cminMeshID, tminMeshID, cmaxMeshID, tmaxMeshID;
  int cminStepID, tminStepID, cmaxStepID, tmaxStepID;
  FC_Variable **sigmas[6]; // index as sigmas[compID][meshID][stepID]

  // handle arguments
  if (argc < 2) {
  usage:
    printf("usage: %s [options] dataset [mesh_name]\n", argv[0]);
    printf("options: \n");
    printf("   -h            : print this help message\n");
    printf("   -v            : print warning and error messages\n");   
    printf("   -V            : print log and error messages\n");   
    printf("\n");
    printf("UPDATE ME\n");
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
      // optional mesh names
      for (j = i+1; j < argc; j++) {
	void* tmp;
	if (!strncmp(argv[j], "-", 1))
	  break;
	tmp = (char**)realloc(mesh_names, (numMesh+1)*sizeof(char*));
	if (!tmp)
	  fc_exitIfError(FC_MEMORY_ERROR);
	mesh_names = tmp;
	mesh_names[numMesh] = (char*)malloc((strlen(argv[j])+1)*sizeof(char));
	strcpy(mesh_names[numMesh], argv[j]);
	numMesh++;
	i++;
      }
    }
  }
  if (!file_name)
    goto usage;

  // init library
  rc = fc_setLibraryVerbosity(verbose_level);
  fc_exitIfError(rc);
  rc = fc_initLibrary();
  fc_exitIfError(rc);

  // load dataset
  rc = fc_loadDataset(file_name, &ds);
  fc_exitIfErrorPrintf(rc, "Failed to load dataset '%s'", file_name);

  // Get meshes
  if (numMesh > 0) {
    meshes = (FC_Mesh*)malloc(numMesh*sizeof(FC_Mesh));
    if (!meshes) 
      fc_exitIfError(FC_MEMORY_ERROR);
    for (i = 0; i < numMesh; i++) {
      FC_Mesh *returnMeshes;
      int numReturnMeshes;

      rc = fc_getMeshByName(ds, mesh_names[i], &numReturnMeshes,
		       &returnMeshes);
      fc_exitIfErrorPrintf(rc, "Failed to get mesh by name");
      if (numReturnMeshes != 1){
	fc_exitIfErrorPrintf(FC_INPUT_ERROR,
			     "Problems finding (unique) mesh '%s' - found %d matches",
			     mesh_names[i],numReturnMeshes);
      }
      meshes[i] = returnMeshes[0];
      free(returnMeshes);
    }
  }
  else {
    rc = fc_getMeshes(ds, &numMesh, &meshes);
    fc_exitIfError(rc);
    mesh_names = (char**)malloc(numMesh*sizeof(char*));
    if (!mesh_names) {
      fc_exitIfError(FC_MEMORY_ERROR);
    }
    for (i = 0; i < numMesh; i++) {
      rc = fc_getMeshName(meshes[i], &mesh_names[i]);
      fc_exitIfError(rc);
    }
  }

  // get the variables
  // FIX? check that numStep is the same?
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
      if (!sigmas[i][j]){
	fc_exitIfErrorPrintf(FC_INPUT_ERROR,
			     "Problems finding (unique) variable '%s' on mesh '%s'",
			     charBuf, mesh_names[j]);
      }
    }
  }
  rc = fc_getSequenceFromSeqVariable(numStep, sigmas[0][0], &sequence);
  fc_exitIfErrorPrintf(rc, "Failed to get sequence from first seq variable");

  // header
  printf("Maximum stresses for dataset '%s'\n", file_name);
  fflush(NULL);
  if (numMesh < 0) {
    printf("No meshes were found\n");
    fflush(NULL);
    exit(0);
  }
  else 
    printf("%d meshes\n", numMesh);
  fflush(NULL);

  // loop over meshes
  for (i = 0; i < numMesh; i++) {
    FC_Variable stress[numStep], pressure[numStep];
    // FC_Variable *temp_seqVar; // for debuggin
    FC_Subset tense, compress;
    int tmin_id, tmax_id;
    double tmin, tmax;
    int cmin_id, cmax_id;
    double cmin, cmax;
    int numTense, numCompress;
    double series_smax;

    // print mesh name
    printf("\n");
    printf("Mesh %d: '%s'\n", i, mesh_names[i]);
    fflush(NULL);


    //make complete vars first inorder to get min max for normalization
    for (j = 0; j < numStep; j++) {
      FC_Variable sigs[6];
      double currmax;

      for (k = 0; k < 6; k++)
	sigs[k] = sigmas[k][i][j];
      
      // get stress & pressure
      rc = fc_createStressPressureVariables(sigs, &stress[j], &pressure[j]);
      fc_exitIfErrorPrintf(rc, "failed to compute stress & pressure");
      rc = fc_getVariableMinMax (stress[j], NULL, NULL, &currmax, NULL);
      if (j == 0){
	series_smax = currmax;
      } else {
	if (currmax >  series_smax){
	  series_smax = currmax;
	}
      }
    }

 
    for (j = 0; j < numStep; j++) {
      printf("  Step %d:\n", j);

      // get compressions & tension subsets
      rc = fc_createThresholdSubset(pressure[j], "<", 0, "compression elems", 
				    &compress);
      fc_getSubsetNumMember(compress, &numCompress);
      rc = fc_createThresholdSubset(pressure[j], ">", 0, "tension elems", 
				    &tense);
      fc_getSubsetNumMember(tense, &numTense);
      
      // get stats
      if (numTense == 0) {
	printf("    No tension elements\n");
      }
      else {
	fc_getVariableSubsetMinMax(stress[j], tense, &tmin, &tmin_id, 
				   &tmax, &tmax_id);
	printf("    Min_tension normalized_to_mesh_max_stress:     %9.6g %9.6g at ID %d\n", tmin, tmin/series_smax, tmin_id);
	printf("    Max_tension normalized_to_mesh_max_stress:     %9.6g %9.6g at ID %d\n", tmax, tmax/series_smax, tmax_id);

	// save extremes for summary
	if (firstt) {
	  min_tmin = tmin;
	  tminMeshID = i;
	  tminStepID = j;
	  max_tmax = tmax;
	  tmaxMeshID = i;
	  tmaxStepID = j;
	  firstt = 0;
	}
	else {
	  if (tmin < min_tmin) {
	    min_tmin = tmin;
	    tminMeshID = i;
	    tminStepID = j;
	  }
	  if (tmax > max_tmax) {
	    max_tmax = tmax;
	    tmaxMeshID = i;
	    tmaxStepID = j;
	  }
	}
      } 
      
      if (numCompress == 0) {
	printf("    No compression elements\n");
      } 
      else {
	fc_getVariableSubsetMinMax(stress[j], compress, &cmin, &cmin_id, 
				   &cmax, &cmax_id);
	printf("    Min_compression normalized_to_mesh_max_stress: %9.6g %9.6g at ID %d\n", cmin, cmin/series_smax, cmin_id);
	printf("    Max_compression normalized_to_mesh_max_stress: %9.6g %9.6g at ID %d\n", cmax, cmax/series_smax, cmax_id);
	// save extremes for summary
	if (firstc) {
	  min_cmin = cmin;
	  cminMeshID = i;
	  cminStepID = j;
	  max_cmax = cmax;
	  cmaxMeshID = i;
	  cmaxStepID = j;
	  firstc = 0;
	}
	else {
	  if (cmin < min_cmin) {
	    min_cmin = cmin;
	    cminMeshID = i;
	    cminStepID = j;
	  }
	  if (cmax > max_cmax) {
	    max_cmax = cmax;
	    cmaxMeshID = i;
	    cmaxStepID = j;
	  }
	}
      }
    }

    // for debugging
    /*
    // convert variables to seq variables so we can see in output
    rc = fc_convertVariablesToSeqVariable(numStep, stress, sequence, NULL,
					  &temp_seqVar);
    fc_exitIfErrorPrintf(rc, "Failed to conver variables to seqVar");
    free(temp_seqVar);
    rc = fc_convertVariablesToSeqVariable(numStep, pressure, sequence, NULL,
					  &temp_seqVar);
    fc_exitIfErrorPrintf(rc, "Failed to conver variables to seqVar");
    free(temp_seqVar);
    */

    // cleanup
    fc_deleteMesh(meshes[i]);
  }

  // Summary
  if (numMesh > 1) {

    printf("\n");
    printf("Summary:\n");

    if (firstt) {
      printf("  No tensioned elements found\n");
    }
    else {
      printf("  Min tension:     %g at step %d on mesh %d:'%s'\n", 
	     min_tmin, tminStepID, tminMeshID, mesh_names[tminMeshID]);
      printf("  Max tension:     %g at step %d on mesh %d:'%s'\n", 
	     max_tmax, tmaxStepID, tmaxMeshID, mesh_names[tmaxMeshID]);
    }
    if (firstc) {
      printf("  No compressed elements found\n");
    }
    else {
      printf("  Min compression: %g at step %d on mesh %d:'%s'\n", 
	     min_cmin, cminStepID, cminMeshID, mesh_names[cminMeshID]);
      printf("  Max compression: %g at step %d on mesh %d:'%s'\n", 
	     max_cmax, cmaxStepID, cmaxMeshID, mesh_names[cmaxMeshID]);
    }
  }

  // for debuggin
  // write new dataset with new data
  //fc_rewriteDataset(ds, "maxStress.ex2", FC_FT_EXODUS);

  // shut down
  for (i = 0; i < numMesh; i++)
    free(mesh_names[i]);
  free(mesh_names);
  for (i = 0; i < 6; i++) {
    for (j = 0; j < numMesh; j++)
      free(sigmas[i][j]);
    free(sigmas[i]);
  }
  fc_deleteDataset(ds);
  fc_finalLibrary();

  exit(0);
}
