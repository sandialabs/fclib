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
 * \file gaps_generator.c
 * \brief A simple program to create a dataset suitable for gap testing.
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/data/gaps_generator.c,v $
 * $Revision: 1.9 $
 * $Date: 2006/08/31 23:57:23 $
 *
 * \description
 *
 *   The meshes will be just a copy of the hex mesh, but moved around in
 *   space. The primary mesh (mesh0) will be stationary. Additional meshes 
 *   will start out touching the primary mesh and then be moved/deformed
 *   over time to cause gaps.
 *     - mesh1: no alignment of verts and edges, translate as (1, 1, 1)
 *     - mesh2: alignment of verts with edges, but not verts with verts,
 *              translate as (1, 1, 0)
 *     - mesh3: alignment of verts, it shears in the x direction proportional
 *              to y coordinate (i.e. top will pull away, but bottom stay
 *              connected)
 *     - mesh4: overlaps more than one side of the mesh, translation as
 *              (0, -1, 0)
 *
 * \modifications:
 *    - 06/12/2006 WSD Created. Used to be internal to gaplines for testing.
 *    - 06/21/2006 WSD Added mesh 3 because was seeing errors w/ a new
 *         algorithm being tested to replace current one in gaplines.
 *    - 08/25/2006 WSD Adding shear case to test normal/tangent resolving.
 *    - 08/25/2006 Changed to build a simple hex mesh instead of loading one.
 *    - 08/25/2006 Added another mesh to test case of more 1 set of sides.
 *    - 08/31/2006 WSD Added death variable.
 */
#include <string.h>
#include "fc.h"
#include "fcP.h"

int main(void) {
  FC_ReturnCode rc;
  int i, j, k;
  FC_Dataset dataset;
  FC_Sequence sequence;
  int numElemPerDim[3] = { 11, 12, 6 };
  FC_Coords lowpoint  = { 0., 0., 0. }, highpoint = { 10., 8., 5. };
  double elemWidths[3];
  int numMesh = 4;
  FC_Mesh meshes[5];
  FC_Variable *displs[5];
  FC_Subset temp_subset;
  FC_Mesh temp_mesh;
  char* mesh_names[5] = { "mesh0", "mesh1", "mesh2", "mesh3", "mesh4" };
  char* displ_name = "displacement";
  double *coords[5], *temp_coords;
  int topodim, numDim, numVertex, numElement, numStep, temp_numStep;
  int* conns, *temp_conns, numVperE;
  FC_Coords move; 
  FC_ElementType elemType;
  double* time_coords;
  double* data;
  int numDead = 18;
  int deadIDs[18] = {   0,   1,   2,    11,  12,  13,  132, 133, 134,
		      143, 144, 145,   264, 265, 266,  275, 276, 277 };

  // --- Setup

  // Init library
  rc = fc_setLibraryVerbosity(FC_WARNING_MESSAGES);
  fc_exitIfError(rc);
  rc = fc_initLibrary();
  fc_exitIfError(rc);

  // create dataset
  fc_createDataset("gapping dataset", &dataset);

  // create a sequence
  numStep = 11;
  time_coords = malloc(numStep*sizeof(double));
  for (i = 0; i < numStep; i++)
    time_coords[i] = i + .1*i;
  fc_createSequence(dataset, "time", &sequence);
  fc_setSequenceCoordsPtr(sequence, numStep, FC_DT_DOUBLE, 
			  (void*)time_coords);

  // make element widths
  for (i = 0; i < 3; i++)
    elemWidths[i] = (highpoint[i] - lowpoint[i])/numElemPerDim[i];
  
  // --- Make meshes 0, 1, 2 & 3 which all have the same mesh conn

  // Mesh0
  // create simple hex mesh
  fc_createSimpleHexMesh(dataset, mesh_names[0], numElemPerDim[0],
			 numElemPerDim[1], numElemPerDim[2], lowpoint,
			 highpoint, &meshes[0]);
  fc_getMeshInfo(meshes[0], &topodim, &numDim, &numVertex, &numElement, 
		 &elemType);
  fc_getMeshCoordsPtr(meshes[0], &coords[0]);
  fc_getMeshElementConnsPtr(meshes[0], &conns);

  // Mesh1 - verts fall inside face 
  // make a second mesh which is on top of the first and moved about halfway
  // over in x and y direction
  move[0] = highpoint[0]/2. + 0.123456789; // trying to prevent elements
  move[1] = highpoint[1]/2. + 0.123456789; // from lining up too nicely
  move[2] = highpoint[2];
  coords[1] = (double*)malloc(numVertex*numDim*sizeof(double));
  for (i = 0; i < numVertex; i++) 
    for (j = 0; j < numDim; j++)
      coords[1][i*numDim + j] = coords[0][i*numDim + j] + move[j];
  fc_createMesh(dataset, mesh_names[1], &meshes[1]);
  fc_setMeshCoordsPtr(meshes[1], numDim, numVertex, coords[1]); // don't copy!
  fc_setMeshElementConns(meshes[1], elemType, numElement, conns); // copy!

  // Mesh2 - verts fall on face edges
  // make a third mesh which is to the back and moved about halfway over
  move[0] = highpoint[0]/2. + 0.123456789; 
  move[1] = highpoint[1];
  move[2] = 0;
  coords[2] = (double*)malloc(numVertex*numDim*sizeof(double));
  for (i = 0; i < numVertex; i++)
    for (j = 0; j < numDim; j++)
      coords[2][i*numDim + j] = coords[0][i*numDim + j] + move[j];
  fc_createMesh(dataset, mesh_names[2], &meshes[2]);
  fc_setMeshCoordsPtr(meshes[2], numDim, numVertex, coords[2]); // don't copy!
  fc_setMeshElementConns(meshes[2], elemType, numElement, conns); // copy!

  // Mesh3 - verts align
  // make a fourth mesh which is directly to the right 
  move[0] = highpoint[0];
  move[1] = 0;
  move[2] = 0;
  coords[3] = (double*)malloc(numVertex*numDim*sizeof(double));
  for (i = 0; i < numVertex; i++) 
    for (j = 0; j < numDim; j++)
      coords[3][i*numDim + j] = coords[0][i*numDim+j] + move[j];
  fc_createMesh(dataset, mesh_names[3], &meshes[3]);
  fc_setMeshCoordsPtr(meshes[3], numDim, numVertex, coords[3]); // don't copy!
  fc_setMeshElementConns(meshes[3], elemType, numElement, conns); // copy!

  // create displacements -
  // mesh0 - 0's
  // mesh1 - translate in 3D (vector = 1,1,1)
  // mesh2 - translate in 2D (vector = 1,1,0)
  // mesh3 - shear in x direction, proportional to y coord
  for (i = 0; i < numMesh; i++) 
    fc_createSeqVariable(meshes[i], sequence, displ_name, &temp_numStep,
			 &displs[i]);
  for (i = 0; i < numStep; i++) {
    // mesh 0
    data = calloc(numVertex*numDim, sizeof(double));
    fc_setVariableData(displs[0][i], numVertex, numDim, FC_AT_VERTEX,
		       FC_MT_VECTOR, FC_DT_DOUBLE, (void*)data);
    // mesh 3
    for (j = 0; j < numDim*numVertex; j += numDim)
      data[j] = i*coords[3][j+1]/10.;
    fc_setVariableData(displs[3][i], numVertex, numDim, FC_AT_VERTEX,
		       FC_MT_VECTOR, FC_DT_DOUBLE, (void*)data);
    // mesh 2
    for (j = 0; j < numDim*numVertex; j += numDim)
      data[j] = i/10.;
    for (j = 1; j < numDim*numVertex; j += numDim)
      data[j] = i/10.;
    fc_setVariableData(displs[2][i], numVertex, numDim, FC_AT_VERTEX,
		       FC_MT_VECTOR, FC_DT_DOUBLE, (void*)data);
    // mesh 1
    for (j = 2; j < numDim*numVertex; j += numDim)
      data[j] = i/10.;
    fc_setVariableData(displs[1][i], numVertex, numDim, FC_AT_VERTEX,
		       FC_MT_VECTOR, FC_DT_DOUBLE, (void*)data);
    free(data);
  }

  // --- Make mesh 4 - has different conn from other meshes

  // Mesh4 - 2 overlapping sides, will translate in -y direction
  // Copy the front lower corner set of elements into a new mesh & translate
  fc_createSubset(meshes[0], "temp", FC_AT_ELEMENT, &temp_subset);
  for (i = 0; i < numElemPerDim[0]; i++) {
    fc_addMemberToSubset(temp_subset, i);
    fc_addMemberToSubset(temp_subset, numElemPerDim[0] + i);
    fc_addMemberToSubset(temp_subset, numElemPerDim[0]*2 + i);
    fc_addMemberToSubset(temp_subset, numElemPerDim[0]*numElemPerDim[1] + i);
    fc_addMemberToSubset(temp_subset, numElemPerDim[0]*numElemPerDim[1]*2 +
			 i);
    }
  fc_createSubsetMesh(temp_subset, dataset, 1, "temp", &temp_mesh);
  fc_deleteSubset(temp_subset);
  fc_getMeshInfo(temp_mesh, &topodim, &numDim, &numVertex, &numElement,
		 &elemType);
  fc_getMeshCoordsPtr(temp_mesh, &temp_coords);
  fc_getMeshElementConnsPtr(temp_mesh, &temp_conns);
  numVperE = fc_getElementTypeNumVertex(elemType);
  conns = (int*)malloc(numElement*numVperE*sizeof(int));
  coords[4] = (double*)malloc(numVertex*numDim*sizeof(double));
  memcpy(conns, temp_conns, numElement*numVperE*sizeof(int));
  memcpy(coords[4], temp_coords, numVertex*numDim*sizeof(double));
  for (i = 0; i < numVertex; i++) {
    for (j = 1; j < numDim; j++)
      coords[4][i*numDim + j] -= elemWidths[j];
  }
  // munge coords so we get face of different orientation:
  // switch vertices 66 & 67
  // this tests the ave face calc in gaplines
  for (i = 0; i < numElement*numVperE; i++) {
    if (conns[i] == 66)
      conns[i] = 67;
    else if (conns[i] == 67)
      conns[i] = 66;
  }
  for (i = 0; i < numDim; i++) {
    double temp_double = coords[4][66*numDim+i];
    coords[4][66*numDim+i] = coords[4][67*numDim+i];
    coords[4][67*numDim+i] = temp_double;
  }
  fc_createMesh(dataset, mesh_names[4], &meshes[4]);
  fc_setMeshCoordsPtr(meshes[4], numDim, numVertex, coords[4]); // don't copy!
  fc_setMeshElementConns(meshes[4], elemType, numElement, conns); // copy!
  fc_deleteMesh(temp_mesh);

  // create displacements
  fc_createSeqVariable(meshes[4], sequence, displ_name, &temp_numStep, 
		       &displs[4]);
  for (i = 0; i < numStep; i++) {
    data = calloc(numVertex*numDim, sizeof(double));
    for (j = 1; j < numDim*numVertex; j += numDim)
      data[j] = -1*i/10.;
    fc_setVariableDataPtr(displs[4][i], numVertex, numDim, FC_AT_VERTEX,
			  FC_MT_VECTOR, FC_DT_DOUBLE, (void*)data);
  }

  // --- Create element death variable for each mesh
  
  for (i = 0; i < numMesh; i++) {
    int* death_data;
    FC_Variable *temp_seqVar;
    fc_getMeshNumElement(meshes[i], &numElement);
    fc_createSeqVariable(meshes[i], sequence, "isDead", &temp_numStep,
			 &temp_seqVar);
    for (j = 0; j < numStep; j++) {
      int stop_index = j;
      if (j == numStep - 1)
	stop_index = numDead; // all die on last step
      death_data = calloc(numElement, sizeof(int));
      for (k = 0; k < stop_index; k++) {
	if (deadIDs[k] < numElement)
	  death_data[deadIDs[k]] = 1;
      }
      fc_setVariableDataPtr(temp_seqVar[j], numElement, 1, FC_AT_ELEMENT,
			    FC_MT_SCALAR, FC_DT_INT, (void*)death_data);
    }
    free(temp_seqVar);
  }

  // --- Write the dataset

  // Write the dataset to a new file
  fc_writeDataset(dataset, "gen_gaps.ex2", FC_FT_EXODUS);

  // --- All done

  // Final library
  fc_finalLibrary();

  exit(0);  
}
