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
 * \file geomExtract.c
 * \brief Extract the requested geometry to a new dataset.
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/tools/geomExtract.c,v $
 * $Revision: 1.5 $ 
 * $Date: 2006/09/28 06:34:39 $
 *
 * \description
 *    Usage: geomExtract [options] dataset [meshes]
 *
 *    Extract the geometry of a single step of a dataset as a
 *    stateless Exodus file.
 *
 * \todo Add handling of dead elements. Clip them out?
 *
 * \modifications
 *   - 09/18/2006 WSD Created.
 */

#include <string.h>
#include "fc.h"

int main(int argc, char** argv) {
  FC_ReturnCode rc;
  int i, j;
  char* dataset_file_name = NULL;
  char* output_file_name = NULL;
  char* displ_var_name = NULL;
  char** meshNames = NULL;
  FC_VerbosityLevel verbose_level = FC_QUIET;
  int numMesh = 0, numStep, temp_numStep;
  FC_Dataset dataset, new_ds;
  FC_Mesh* meshes;
  int doDispl = 0;
  FC_Variable **displs;  // displs[meshID][stepID]
  int stepID = -1;  // the id of the timestep of interest
  char* temp_name;

  // --- handle arguments

  if (argc < 2 ) {
  usage:
    printf("usage: %s [options] dataset [mesh_names]\n", 
           argv[0]);
    printf("options: \n");
    printf("   -h               : print this help message\n");
    printf("   -v               : verbose: print warning and error messages\n");
    printf("   -V               : very verbose: prints log and error messages\n");
    printf("   -o dataset_name  : name of output dataset (default is extracted.g)\n");
    printf("   -d var_name      : the displacement variable\n");
    printf("   -t time_step     : The index of the timestep at which to find tears\n");
    printf("                      (indexing from 0). The default is to use the last one.\n");
    printf("\n");
    printf("Extracts all meshes (or the specific meshes provided as arguments)\n");
    printf("into a stateless Exodus file. If a displacement variable and a time step\n");
    printf("are provided, this will extract the displaced geometry at that step.\n");
    printf("\n");
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
    else if (!strncmp(argv[i], "-h", 2))
      goto usage;
    else if (!strcmp(argv[i], "-o")) {
      i++;
      output_file_name = argv[i];
    }
    else if (!strcmp(argv[i], "-d")) {
      doDispl = 1;
      i++;
      displ_var_name = argv[i];
    }
    else if (!strcmp(argv[i], "-t")) {
      i++;
      stepID = atoi(argv[i]);
    }
    else {
      dataset_file_name = argv[i];
      // optional mesh names
      for (j =i+1; j < argc; j++) {
	void* tmp;
	if (argv[j][0] == '-') 
	  break;
	tmp = (char**)realloc(meshNames, (numMesh+1)*sizeof(char*));
	if (!tmp)
	  fc_exitIfError(FC_MEMORY_ERROR);
	meshNames = tmp;
	meshNames[numMesh] = (char*)malloc((strlen(argv[j])+1)*sizeof(char));
	strcpy(meshNames[numMesh], argv[j]);
	numMesh++;
	i++;
      }
    }
  }

  // testing args
  if (!dataset_file_name)
    goto usage;
  if (output_file_name == NULL)
    output_file_name = "extracted.g";

  // --- setup

  // init library and load dataset 
  rc = fc_setLibraryVerbosity(verbose_level);
  fc_exitIfError(rc);
  rc = fc_initLibrary();
  fc_exitIfError(rc);
  rc = fc_loadDataset(dataset_file_name, &dataset);
  fc_exitIfErrorPrintf(rc, "Failed to load dataset '%s'", dataset_file_name);
  

  // Make sure meshes exist (or get all meshes) {
  if (numMesh > 0) {
    meshes = (FC_Mesh*)malloc(numMesh*sizeof(FC_Mesh));
    if (!meshes)
      fc_exitIfError(FC_MEMORY_ERROR);
    for (i = 0; i < numMesh; i++) {
      FC_Mesh *returnMeshes;
      int numReturnMeshes;
      rc = fc_getMeshByName(dataset, meshNames[i], &numReturnMeshes,
			    &returnMeshes);
      fc_exitIfErrorPrintf(rc, "failed to get mesh by name");
      if (numReturnMeshes != 1){
	fc_exitIfErrorPrintf(FC_INPUT_ERROR,
			     "Failed to find (unique) mesh '%s' - found %d matches",
			     meshNames[i],numReturnMeshes);
      }
      meshes[i] = returnMeshes[0];
      free(returnMeshes);
    }
  }
  else {
    rc = fc_getMeshes(dataset, &numMesh, &meshes);
    fc_exitIfError(rc);
    if (numMesh < 1)
      fc_exitIfErrorPrintf(FC_ERROR, "No meshes in dataset");
    meshNames = (char**)malloc(numMesh*sizeof(char*));
    if (!meshNames)
      fc_exitIfError(FC_MEMORY_ERROR);
    for (i = 0; i < numMesh; i++) {
      rc = fc_getMeshName(meshes[i], &meshNames[i]);
      fc_exitIfErrorPrintf(rc, "Failed to get mesh name");
    }
  }

  // Make sure required variable fields exist on the meshes
  if (doDispl) {
    displs = (FC_Variable**)malloc(numMesh*sizeof(FC_Variable*));
    if (!displs)
      fc_exitIfError(FC_MEMORY_ERROR);
    for (i = 0; i < numMesh; i++) {
      rc = fc_getOrGenerateUniqueSeqVariableByName(meshes[i], displ_var_name,
						   &temp_numStep,
						   &displs[i]);
      fc_exitIfErrorPrintf(rc, "Failed to find displacement variable '%s'",
			   displ_var_name);
      if (!displs[i]){
	fc_exitIfErrorPrintf(FC_INPUT_ERROR,
			     "Failed to find (unique) displacement variable '%s'", displ_var_name);
      }


      if (i == 0)
	numStep = temp_numStep;
      else {
	if (temp_numStep != numStep) {
	  fc_exitIfErrorPrintf(FC_ERROR, "the displacement vars do not have "
			       "the same number of steps on each mesh");
	}
      }
    }
    // setup the stepID
    if (stepID == -1) {
      stepID = numStep-1;
    }
    else {
      if (stepID < 0 || stepID > numStep-1)
	fc_exitIfErrorPrintf(FC_ERROR, "stepID (%d) is out of range (0-%d)",
			     stepID, numStep-1);
    }
  }

  // --- Set up new dataset
  
  rc = fc_getDatasetName(dataset, &temp_name);
  fc_exitIfErrorPrintf(rc, "Failed to get dataset name");
  rc = fc_createDataset(temp_name, &new_ds);
  fc_exitIfErrorPrintf(rc, "Failed to create new dataset");
  free(temp_name);

  // --- Copy mesh geometries
  
  for (i = 0; i < numMesh; i++) {
    FC_Mesh new_mesh;
    FC_ElementType elemType;
    int topodim, numDim, numVert, numElem, *conns;
    double *coords;

    // Create new mesh
    rc = fc_getMeshName(meshes[i], &temp_name);
    fc_exitIfErrorPrintf(rc, "failed to get mesh name");
    rc = fc_createMesh(new_ds, temp_name, &new_mesh);
    fc_exitIfErrorPrintf(rc, "failed to create new mesh");
    free(temp_name);

    // Get mesh info & data (make a copy of coords & hand off)
    rc = fc_getMeshInfo(meshes[i], &topodim, &numDim, &numVert, &numElem,
			&elemType);
    fc_exitIfErrorPrintf(rc, "failed to get mesh info");
    if (doDispl) {
      rc = fc_getDisplacedMeshCoords(meshes[i], displs[i][stepID], &coords);
      fc_exitIfErrorPrintf(rc, "failed to get displaced mesh coords");
    }
    else {
      double* temp_coords;
      rc = fc_getMeshCoordsPtr(meshes[i], &temp_coords);
      fc_exitIfErrorPrintf(rc, "failed to get mesh coords");
      coords = malloc(numDim*numVert*sizeof(double));
      memcpy(coords, temp_coords, numDim*numVert*sizeof(double));
    }
    rc = fc_getMeshElementConnsPtr(meshes[i], &conns);
    fc_exitIfErrorPrintf(rc, "failed to get mesh conns");

    // Add mesh coords
    rc = fc_setMeshCoordsPtr(new_mesh, numDim, numVert, coords);
    fc_exitIfErrorPrintf(rc, "failed to set mesh coords");

    // Add Mesh conns
    rc = fc_setMeshElementConns(new_mesh, elemType, numElem, conns);
    fc_exitIfErrorPrintf(rc, "failed to set mesh conns");
  }

  // --- Write new dataset

  rc = fc_writeDataset(new_ds, output_file_name, FC_FT_EXODUS);
  fc_exitIfErrorPrintf(rc, "failed to write new dataset %s", 
		       output_file_name);

  // --- cleanup

  if (doDispl) {
    for (i = 0; i < numMesh; i++) 
      free(displs[i]);
    free(displs);
  }
  for (i = 0; i < numMesh; i++) 
    free(meshNames[i]);
  free(meshNames);
  free(meshes);
  fc_deleteDataset(dataset);
  fc_finalLibrary();

  exit(0);

}
