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
 * \file edgeStats.c
 * \brief Report edge length statistics on a dataset
 *
 * $Source: /home/Repositories/fcdmf/fclib/tools/edgeStats.c,v $
 * $Revision: 1.22 $ 
 * $Date: 2006/11/07 23:49:03 $
 *
 * \description
 *    Usage: edgeStats [options] dataset [mesh_name]
 *
 *    Reports edge length statistics. If a mesh is given, it reports
 *    the stats for only that mesh, otherwise it reports stats for
 *    all meshes.
 *
 * \todo ?Add a flag to get overall average of edges?
 *
 * \modifications
 *   - 04/20/04 WSK, created.
 */

#include <string.h>
#include "fc.h"
 
//  pretty output depends on var already being checked to be a valid
// displacement variable
static void get_and_print_lengths(FC_Mesh mesh, FC_Variable var,
				  char* message) {
  FC_ReturnCode rc;
  FC_Variable lengths;
  double min, max, ave, stdev;

  printf("    %s\n", message);
  fflush(NULL);

  // get lengths  
  if (FC_HANDLE_EQUIV(var, FC_NULL_VARIABLE))
    rc = fc_getEdgeLengths(mesh, &lengths);
  else
    rc = fc_getDisplacedEdgeLengths(mesh, var, &lengths);
  if (rc != FC_SUCCESS) {
    printf("        Failed to get edge lengths\n");
    fflush(NULL);
    return;
  }

  // do stats
  rc = fc_getVariableMinMax(lengths, &min, NULL, &max, NULL);
  if (rc == FC_SUCCESS) {
    printf("        min = %g\n", min);
    printf("        max = %g\n", max);
  }
  else 
    printf("        Failed to get min/max lengths\n");
  fflush(NULL);
  rc = fc_getVariableMeanSdev(lengths, &ave, &stdev);
  if (rc == FC_SUCCESS) 
    printf("        ave = %g +/- %g\n", ave, stdev);
  else
    printf("        Failed to get mean/stdev lengths\n");
  fflush(NULL);

  // cleanup
  fc_deleteVariable(lengths);

  return;
}

int main(int argc, char** argv) {
  FC_ReturnCode rc;
  int i, j;
  int numMesh, numElem, numEdge, *numSteps;
  FC_Dataset ds;
  FC_Mesh *meshes;
  FC_Variable *displs, **seqDispls;   
  char* file_name = NULL;
  char* mesh_name = NULL;
  char* displ_name = NULL;
  char message_buf[1028];
  FC_VerbosityLevel verbose_level = FC_QUIET;
  int doOnlyOne = 0; // flag for whether to do one (1) or all (0);
  FC_ElementType elemType;
  int doDispl = 0;
  int* foundDispls, *foundSeqDispls, numNoDispl;

  // handle arguments
  if (argc < 2) {
  usage:
    printf("usage: %s [options] dataset [mesh_name]\n", argv[0]);
    printf("options: \n");
    printf("   -d var_name   : calc the displaced values using this displacement\n");
    printf("                   variable\n");
    printf("   -h            : print this help message\n");
    printf("   -v            : print warning and error messages\n");   
    printf("   -v            : print log and error messages\n");   
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
    int numReturnMeshes;
    FC_Mesh *returnMeshes;

    numMesh = 1;
    meshes = malloc(sizeof(FC_Mesh));
    rc = fc_getMeshByName(ds, mesh_name, &numReturnMeshes,&returnMeshes);
    fc_exitIfErrorPrintf(rc,"Failed to get mesh by name");
    if (numReturnMeshes != 1){
      fc_exitIfErrorPrintf(FC_INPUT_ERROR, "Failed to find (unique) mesh '%s' - found %d matches",
             mesh_name,numReturnMeshes);
    }
    meshes[0] = returnMeshes[0];
    free(returnMeshes);
  }
  else {
    rc = fc_getMeshes(ds, &numMesh, &meshes);
    fc_exitIfError(rc);
  }

  // header
  printf("Edge length summary for dataset '%s'\n", file_name);
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
	! numSteps)
      fc_exitIfError(FC_MEMORY_ERROR);
		    
    // try to get handles to vars
    numNoDispl = 0;
    for (i = 0; i < numMesh; i++) {
      foundDispls[i] = 0;
      foundSeqDispls[i] = 0;

      rc = fc_getOrGenerateUniqueVariableByName(meshes[i], displ_name,
						&displs[i]);
      fc_exitIfError(rc);						
      if (fc_isValidDisplacementVariable(meshes[i], displs[i])) 
	foundDispls[i] = 1;

      rc = fc_getOrGenerateUniqueSeqVariableByName(meshes[i], displ_name,
						   &numSteps[i], 
						   &seqDispls[i]);
      fc_exitIfError(rc);						
      if (seqDispls[i] && 
	  fc_isValidDisplacementVariable(meshes[i], seqDispls[i][0]))
	foundSeqDispls[i] = 1;

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

  // loop over meshes
  for (i = 0; i < numMesh; i++) {

    // describe mesh
    rc = fc_getMeshName(meshes[i], &mesh_name);
    rc = fc_getMeshInfo(meshes[i], NULL, NULL, NULL, &numElem, &elemType);
    rc = fc_getMeshNumEdge(meshes[i], &numEdge);
    fc_exitIfError(rc);

    // print mesh info
    printf("\n");
    printf("Mesh %d: '%s'\n", i, mesh_name);
    free(mesh_name);
    printf("    numElement = %d,  elemType = %s\n", numElem,
           fc_getElementTypeText(elemType));
    printf("    numEdge = %d\n", numEdge);
    fflush(NULL);
 
    // analyze edges, if any
    if (numEdge > 0) {
      if (!doDispl)
	get_and_print_lengths(meshes[i], FC_NULL_VARIABLE, "Edge lengths:");
      else {
	// do basic var
	if (foundDispls[i])
	  get_and_print_lengths(meshes[i], displs[i], "Displaced edge lengths:");

	// do a seq variable	
	if (foundSeqDispls[i]) {
	  for (j = 0; j < numSteps[i]; j++) {
	    sprintf(message_buf, "Displaced edge lengths for step %d", j);
	    get_and_print_lengths(meshes[i], seqDispls[i][j], message_buf);
	  }
	}
	
	// didn't find one
	if (foundDispls[i] + foundSeqDispls[i] == 0) {
	  printf("    (Displacement var not found, cannot calc edge lengths)\n");
	  fflush(NULL);
	}
      }
    }

    fc_deleteMesh(meshes[i]);
  }

  // cleanup
  if (doDispl) {
    free(foundDispls);
    free(foundSeqDispls);
    for (i = 0; i < numMesh; i++){
      free(seqDispls[i]); //dont delete the var
    }
    free(seqDispls);
    free(numSteps);
    free(displs);
  }

  // shut down
  fc_deleteDataset(ds);
  fc_finalLibrary();

  exit(0);
}
