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
 * \file extents.c
 * \brief Report the extents of a dataset and its meshes.
 *
 * $Source: /home/Repositories/fcdmf/fclib/tools/extents.c,v $
 * $Revision: 1.10 $ 
 * $Date: 2006/11/17 22:51:25 $
 *
 * \description
 *    Usage: extents [options] dataset [mesh_name]
 *
 *    Reports the bounding box of each mesh in a dataset and the bounding box
 *    of the entired dataset.
 *
 * \todo ?Add a flag to get overall average of edges?
 *
 * \modifications
 *   - 05/24/2006 WSD Created.
 */

#include <string.h>
#include "fc.h"
 
// print extents in format: "  %s = [ %g, %g, %g ] - [ %g, %g, %g ]\n"
static void print_extent(char* label, int numDim, FC_Coords lowpoint, 
			 FC_Coords highpoint) {
  int i;
  printf("  %s = [ %g", label, lowpoint[0]);
  for (i = 1; i < numDim; i++)
    printf(", %g", lowpoint[i]);
  printf(" ] - [ %g", highpoint[0]);
  for (i = 1; i < numDim; i++)
    printf(", %g", highpoint[i]);
  printf(" ]\n");
}

// combine bounding boxes
// accumulate results in the combined bb
static FC_ReturnCode combine_bb(int numDim, FC_Coords low, FC_Coords high, 
				FC_Coords *lowc, FC_Coords *highc) {
  int i;
  FC_Coords temp_low, temp_high;

  // assume low & high are real, but lowc & highc might not be
  // if lowc,highc is not real copy low/high
  for (i = 0; i < 3; i++) {
    if ((*highc)[i] < (*lowc)[i]) {
      (*lowc)[i] = low[i];
      (*highc)[i] = high[i];
    }
  }

  // do it
  for (i = 0; i < numDim; i++) {
    temp_low[i] = (*lowc)[i];
    temp_high[i] = (*highc)[i];
  }
  return fc_combineBoundingBoxes(numDim, low, high, temp_low, temp_high, 
				 lowc, highc);
}

int main(int argc, char** argv) {
  FC_ReturnCode rc;
  int i, j;
  int numMesh, numDim, ds_numDim = 0, *numSteps, maxNumStep;
  FC_Dataset ds;
  FC_Mesh *meshes;
  FC_Variable *displs, **seqDispls;   
  char* file_name = NULL;
  char* mesh_name = NULL;
  char* displ_name = NULL;
  FC_VerbosityLevel verbose_level = FC_QUIET;
  int doOnlyOne = 0; // flag for whether to do one (1) or all (0);
  int doDispl = 0;
  int foundDispl = 0, foundSeqDispl = 0;
  int* foundDispls, *foundSeqDispls, numNoDispl;
  FC_Coords ds_low, ds_high, displ_ds_low, displ_ds_high;
  FC_Coords *seqDispl_ds_lows, *seqDispl_ds_highs;

  // handle arguments
  if (argc < 2) {
  usage:
    printf("usage: %s [options] dataset [mesh_name]\n", argv[0]);
    printf("options: \n");
    printf("   -d var_name   : calc the displaced values using this displacement\n");
    printf("                   variable\n");
    printf("   -h            : print this help message\n");
    printf("   -v            : print warning and error messages\n");   
    printf("   -V            : print log and error messages\n");   
    printf("\n");
    printf("Prints summary of edge lengths on meshes in the dataset.\n");
    printf("If a mesh name is provided, prints only a summary for that "
           "mesh.\n");
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
      file_name = argv[i];
      if (i+1 < argc) {
        i++;
        doOnlyOne = 1;
        mesh_name = argv[i];
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
  if (doOnlyOne) {
    FC_Mesh *returnMeshes;
    int numReturnMeshes;
    numMesh = 1;
    meshes = malloc(sizeof(FC_Mesh));
    rc = fc_getMeshByName(ds, mesh_name, &numReturnMeshes,&returnMeshes);
    fc_exitIfErrorPrintf(rc, "Can't get mesh by name");
    if (numReturnMeshes != 1){
      fc_exitIfErrorPrintf(FC_INPUT_ERROR, "Problems finding (unique) mesh '%s' - found %d matches",
             mesh_name, numReturnMeshes);
    }
    meshes[0] = returnMeshes[0];
    free(returnMeshes);
  }
  else {
    rc = fc_getMeshes(ds, &numMesh, &meshes);
    fc_exitIfError(rc);
  }

  // header
  printf("Extents summary for dataset '%s'\n", file_name);
  fflush(NULL);
  if (numMesh < 0) {
    printf("No meshes were found\n");
    fflush(NULL);
    exit(0);
  }
  else 
    printf("%d meshes\n", numMesh);
  fflush(NULL);

  // if displacing, make sure we have displacement variables
  if (doDispl) {
    foundDispls = (int*)malloc(numMesh*sizeof(int));
    foundSeqDispls = (int*)malloc(numMesh*sizeof(int));
    displs = (FC_Variable*)malloc(numMesh*sizeof(FC_Variable));
    numSteps = (int*)malloc(numMesh*sizeof(int));
    seqDispls = (FC_Variable**)malloc(numMesh*sizeof(FC_Variable*));
    if (!foundDispls || !foundSeqDispls || !displs || !seqDispls ||
	!numSteps)
      fc_exitIfError(FC_MEMORY_ERROR);
		    
    // try to get handles to vars
    numNoDispl = 0;
    for (i = 0; i < numMesh; i++) {
      foundDispls[i] = 0;
      foundSeqDispls[i] = 0;
      rc = fc_getOrGenerateUniqueVariableByName(meshes[i],
						displ_name,
						&displs[i]);
      fc_exitIfError(rc);

      if (fc_isValidDisplacementVariable(meshes[i], displs[i])) {
	foundDispl = 1;
	foundDispls[i] = 1;
      }

      rc = fc_getOrGenerateUniqueSeqVariableByName(meshes[i], displ_name,
						   &numSteps[i],
						   &seqDispls[i]);
      fc_exitIfError(rc);
      
      if (seqDispls[i] && 
	  fc_isValidDisplacementVariable(meshes[i], seqDispls[i][0])){
	foundSeqDispl = 1;
	foundSeqDispls[i] = 1;
      }
      if (foundDispls[i] + foundSeqDispls[i] == 0) 
	numNoDispl++;
    }

 // check
    if (numNoDispl == numMesh) {
      fc_exitIfErrorPrintf(FC_ERROR, "Displacement variable '%s' was not found",
			   displ_name);
    }
    else if (numNoDispl > 0) {
      int first = 1;
      printf("WARNING: the displacement variable does not exist on meshes ");
      first = 0;
      for (i = 0; i < numMesh; i++) {
	if (foundDispls[i] + foundSeqDispls[i] == 0) {
	  if (first)
	    printf("%d", i);
	  else
	    printf(", %d", i);
	}
      }
      printf("\n");
    }
  }

  // initialize overall dataset bb's
  if (doDispl) {
    maxNumStep = 0;
    for (i = 0; i < numMesh; i++) {
      if (numSteps[i] > maxNumStep)
	maxNumStep = numSteps[i];
    }
    if (maxNumStep > 0) {
      seqDispl_ds_lows = (FC_Coords*)malloc(maxNumStep*sizeof(FC_Coords));
      seqDispl_ds_highs = (FC_Coords*)malloc(maxNumStep*sizeof(FC_Coords));	
      if (!seqDispl_ds_lows || !seqDispl_ds_highs)
	fc_exitIfError(FC_MEMORY_ERROR);
      for (i = 0; i < maxNumStep; i++) {
	for (j = 0; j < 3; j++) {
	  seqDispl_ds_lows[i][j] = -1;
	  seqDispl_ds_highs[i][j] = -2;
	}
      }
    }
  }
  for (i = 0; i < 3; i++) {
    ds_low[i] = -1;
    ds_high[i] = -2;
    displ_ds_low[i] = -1;
    displ_ds_high[i] = -2;
  }

  // loop over meshes
  for (i = 0; i < numMesh; i++) {
    FC_Coords lowpoint, highpoint;

    // print mesh name
    rc = fc_getMeshName(meshes[i], &mesh_name);
    fc_exitIfError(rc);
    printf("\n");
    printf("Mesh %d: '%s'\n", i, mesh_name);
    free(mesh_name);
    fflush(NULL);
 
    // Get bounding box of undisplaced mesh
    rc = fc_getMeshBoundingBox(meshes[i], &numDim, &lowpoint, &highpoint);
    fc_exitIfErrorPrintf(rc, "Failed to get mesh bounding box");
    print_extent("Extent", numDim, lowpoint, highpoint);
    if (numDim > ds_numDim)
      ds_numDim = numDim;
    rc = combine_bb(numDim, lowpoint, highpoint, &ds_low, &ds_high);
    fc_exitIfError(rc);

    if (doDispl) { 
      // if didn't find displs
      if (foundDispls[i] + foundSeqDispls[i] == 0) {
	printf("    (Displacement var not found, cannot calc extents)\n");
	fflush(NULL);
	break;
      }

      // do basic displ var
      if (foundDispls[i]) {
	rc = fc_getDisplacedMeshBoundingBox(meshes[i], displs[i], &numDim,
					    &lowpoint, &highpoint);
	fc_exitIfErrorPrintf(rc, "Failed to get displaced mesh bounding box");
	print_extent("Displ extent", numDim, lowpoint, highpoint);
	rc = combine_bb(numDim, lowpoint, highpoint, &displ_ds_low, &displ_ds_high);
	fc_exitIfError(rc);
      }

      // do a seq displ variable	
      if (foundSeqDispls[i]) {
	FC_Coords combinedlow, combinedhigh;
	for (j = 0; j < 3; j++) {
	  combinedlow[j] = -1;
	  combinedhigh[j] = -2;
	}
	for (j = 0; j < numSteps[i]; j++) {
	  char strbuf[1028];
	  rc = fc_getDisplacedMeshBoundingBox(meshes[i], seqDispls[i][j], 
					      &numDim, &lowpoint, &highpoint);
	  fc_exitIfErrorPrintf(rc, "Failed to get displaced mesh bounding box");
	  sprintf(strbuf, "Displ extent at Step %d", j);
	  print_extent(strbuf, numDim, lowpoint, highpoint);
	  rc = combine_bb(numDim, lowpoint, highpoint, &combinedlow, 
			  &combinedhigh);
	  fc_exitIfError(rc);
	  rc = combine_bb(numDim, lowpoint, highpoint, &seqDispl_ds_lows[j],
			  &seqDispl_ds_highs[j]);
	  fc_exitIfError(rc);
	}
	print_extent("Displ extent for all Steps", numDim, combinedlow,
		     combinedhigh);
      }
    }

    fc_deleteMesh(meshes[i]);
  }

  // Summary
  if (numMesh > 1) {
    printf("\n");
    printf("Entire Dataset:\n");
    print_extent("Extent", ds_numDim, ds_low, ds_high);
    if (doDispl) {
      if (foundDispl) {
	print_extent("Displ extent", ds_numDim, displ_ds_low, displ_ds_high);
      }
      if (foundSeqDispl) {
	FC_Coords combinedlow, combinedhigh;
	for (i = 0; i < 3; i++) {
	  combinedlow[i] = -1;
	  combinedhigh[i] = -2;
	}
	for (i = 0; i < maxNumStep; i++) {
	  char strbuf[1028];
	  sprintf(strbuf, "Displ extent at Step %d", i);
	  print_extent(strbuf, ds_numDim, seqDispl_ds_lows[i],
		       seqDispl_ds_highs[i]);
	  rc = combine_bb(ds_numDim, seqDispl_ds_lows[i], seqDispl_ds_highs[i],
			  &combinedlow, &combinedhigh);
	  fc_exitIfError(rc);
	}
	print_extent("Displ extent for all Steps", ds_numDim, combinedlow,
		     combinedhigh);
      }
    }
  }
 

  // cleanup
  if (doDispl) {
    free(foundDispls);
    free(foundSeqDispls);
    if (maxNumStep > 0) {
      free(seqDispl_ds_lows);
      free(seqDispl_ds_highs);
    }
    for (i = 0; i < numMesh; i++){
      if (seqDispls[i]) free(seqDispls[i]);
    }
    free(numSteps);
    free(displs);
    free(seqDispls);
  }

  // shut down
  fc_deleteDataset(ds);
  fc_finalLibrary();

  exit(0);
}
