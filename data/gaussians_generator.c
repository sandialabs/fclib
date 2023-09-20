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
 * \file gaussians_generator.c
 * \brief A simple program to create dataset with interesting features.
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/data/gaussians_generator.c,v $
 * $Revision: 1.12 $
 * $Date: 2006/09/19 00:57:58 $
 *
 * \description
 *
 *   Originally written to create datasets for IDA 2005 paper to as a simple
 *   example of features. The features are hot spots modeled by a
 *   four gaussians, one in each quadrant. Each gaussian is scaled over
 *   time by a linear factor.
 *   
 *   The temperature field is vertex centered data on 100x100 grid of 
 *   vertices connected into quads. Another variable, displacement,
 *   maps the temperature to a height above the original grid and
 *   is a vector sequence variable. 
 *
 *   Also added some arbitrary subsets to make this dataset is more useful for
 *   testing. And added an element variable "damage" that is the temperature
 *   data interpolated onto the elements.
 *
 *   \modifications
 *      - 02/2005 WSK Created.
 *      - 04/11/05 WSK Added to repository.
 */

#include <fc.h>

// This version goes to 1 as x -> mu and zero as x get very far from mu
static double normal_dist_1D(double x, double mu, double sigma) {
  return exp(-1./2.*pow((x-mu)/sigma,2));
}

static double normal_dist_2D(double* X, double* mu, double sigma) {
  return normal_dist_1D(X[0], mu[0], sigma) * 
         normal_dist_1D(X[1], mu[1], sigma);
}

int main(void) {
  // generic
  FC_ReturnCode rc;
  int i, j, k, idx;

  // Names and handles
  char dataset_file_name_exo[] = "gen_gaussians.ex2";
  char dataset_name[] = "gaussians";
  char sequence_name[] = "time";
  char mesh_name[] = "grid";
  char var_name[] = "temperature";
  char var2_name[] = "damage";
  char vert_subset_name[] = "bottom edge";
  char elem_subset_name[] = "SW quadrant";
  FC_Dataset dataset;
  FC_Sequence sequence;
  FC_Mesh mesh;
  FC_Variable *seqVar, *seqVar2, *displs;
  FC_Subset elemSubset, vertSubset;

  // Time information
  int numStep = 11;
  double time_origin = 0;
  double step_size = 1;
  double *time_coords;
   
  // The underlying grid is assumed to be 2D, embedded in 3D space
  int topoDim = 2;
  int numDim = 3;
  int numVertPerDim[2] = { 101, 101 }, numElem, numVert;
  FC_Coords origin = { 0, 0, 0 };
  FC_Coords extent = { 100, 100, 0 };
  double* coords;
  int* conns, *conn_p;

  // The data
  double *data, *data2, *displ_data;

  // Gaussians
  int numPt = 4;     // number of gaussians
  FC_Coords centers[4] =  // where
    { { 25, 25, 0 }, 
      { 25, 75, 0 },
      { 75, 25, 0 },
      { 75, 75, 0 } };
  float sigmas[4] =      // widths
    { 15, 15, 15, 15 };
  float initial_scales[4] = // initial scale
    { 50, 80, 50, 10 };
  float rates[4] =       // rate that scale changes at
    { 4, -15, 8, 6 };
 
  // Setup the dataset
  rc = fc_initLibrary();
  fc_exitIfErrorPrintf(rc, "failed to initialize fc library");
  rc = fc_setLibraryVerbosity(FC_WARNING_MESSAGES);
  rc = fc_createDataset(dataset_name, &dataset);
  fc_exitIfErrorPrintf(rc, "Failed to create dataset");

  // Create the sequence
  time_coords = (double*)malloc(numStep*sizeof(double));
  for (i = 0; i < numStep; i++)
    time_coords[i] = time_origin + i*step_size;
  rc = fc_createSequence(dataset, sequence_name, &sequence);
  fc_exitIfErrorPrintf(rc, "failed to create sequence");
  rc = fc_setSequenceCoordsPtr(sequence, numStep, FC_DT_DOUBLE, time_coords);
  fc_exitIfErrorPrintf(rc, "failed to set time values");

  // Create underlying mesh
  // Fastest varying index is in X direction, then Y
  numVert = numVertPerDim[0] * numVertPerDim[1];
  coords = (double*)malloc(numVert*numDim*sizeof(double));
  // make coordinates
  for (i = 0; i < numVertPerDim[0]; i++) {
    for (j = 0; j < numVertPerDim[1]; j++) { 
      idx = j*numVertPerDim[0] + i;
      coords[idx*numDim] = (extent[0] - origin[0])*i/(numVertPerDim[0]-1);
      coords[idx*numDim + 1] = (extent[1] - origin[1])*j/(numVertPerDim[1]-1);
      coords[idx*numDim + 2] = 0.;
    }
  }
  // make conns
  numElem = (numVertPerDim[0] - 1) * (numVertPerDim[1] - 1);
  conns = (int*)malloc(numElem*4*sizeof(int));
  conn_p = conns;
  for (j = 0; j < numVertPerDim[1]-1; j++) {
    for (i = 0; i < numVertPerDim[0]-1; i++) {
      conn_p[0] = j*numVertPerDim[0] + i;
      conn_p[1] = j*numVertPerDim[0] + i+1;
      conn_p[2] = (j+1)*numVertPerDim[0] + i+1;
      conn_p[3] = (j+1)*numVertPerDim[0] + i;
      conn_p += 4;
    }
  }
  rc = fc_createMesh(dataset, mesh_name, &mesh);
  fc_exitIfErrorPrintf(rc, "Failed to create mesh");
  rc = fc_setMeshCoordsPtr(mesh, numDim, numVert, coords);
  fc_exitIfErrorPrintf(rc, "Failed to set mesh coords");
  rc = fc_setMeshElementConnsPtr(mesh, FC_ET_QUAD, numElem, conns);
  fc_exitIfErrorPrintf(rc, "Failed to set mesh conns");

  // Create the sequence variables - temperature & damage
  rc = fc_createSeqVariable(mesh, sequence, var_name, &numStep,
			    &seqVar);
  fc_exitIfErrorPrintf(rc, "Failed to create sequence variable");
  rc = fc_createSeqVariable(mesh, sequence, var2_name, &numStep,
			    &seqVar2);
  for (i = 0; i < numStep; i++) {
    // create vert data for this step
    data = malloc(numVert*sizeof(double));
    data2 = malloc(numElem*sizeof(double));
    for (j = 0; j < numVert; j++) {
      double *X = &coords[j*numDim];
      data[j] = 0;
      for (k = 0; k < numPt; k++) {
	data[j] += (initial_scales[k] + i*rates[k]) *
	  normal_dist_2D(X, centers[k], sigmas[k]);
      }
    }
    // create elem data for this step -- average vert values
    for (j = 0; j < numVertPerDim[0] - 1; j++) {
      for (k = 0; k < numVertPerDim[1] - 1; k++) {
	data2[j*(numVertPerDim[0]-1)+k] = data[j*numVertPerDim[0]+k];
	data2[j*(numVertPerDim[0]-1)+k] += data[(j+1)*numVertPerDim[0]+k];
	data2[j*(numVertPerDim[0]-1)+k] += data[j*numVertPerDim[0]+k+1];
	data2[j*(numVertPerDim[0]-1)+k] += data[(j+1)*numVertPerDim[0]+k+1];
	data2[j*(numVertPerDim[0]-1)+k] /= 4.;
      }
    }
    rc = fc_setVariableDataPtr(seqVar[i], numVert, 1, FC_AT_VERTEX,
			       FC_MT_SCALAR, FC_DT_DOUBLE, data);
    fc_exitIfErrorPrintf(rc, "failed to set data for step %d of vert "
			 "variable", i);
    rc = fc_setVariableDataPtr(seqVar2[i], numElem, 1, FC_AT_ELEMENT,
			       FC_MT_SCALAR, FC_DT_DOUBLE, data2);
    fc_exitIfErrorPrintf(rc, "failed to set data for step %d of elem "
			 "variable", i);
  }
  
  // If 2D, create a displacement var from the seq var for pretty pictures
  if (topoDim == 2) {
    rc = fc_createSeqVariable(mesh, sequence, "displacement", &numStep,
			      &displs);
    fc_exitIfErrorPrintf(rc, "failed to create displacement variable");
    for (i = 0; i < numStep; i++) {
      rc = fc_getVariableDataPtr(seqVar[i], (void**)&data);
      displ_data = malloc(numVert*3*sizeof(double));
      for (j = 0; j < numVert; j++) {
	displ_data[j*3] = 0;
	displ_data[j*3 + 1] = 0;
	displ_data[j*3 + 2] = data[j]/2.;
      }
      rc = fc_setVariableDataPtr(displs[i], numVert, 3, FC_AT_VERTEX,
				 FC_MT_VECTOR, FC_DT_DOUBLE, displ_data);
    }
  }

  // Add some subsets to make more interesting
  rc = fc_createSubset(mesh, vert_subset_name, FC_AT_VERTEX, &vertSubset);
  fc_exitIfErrorPrintf(rc, "failed to create vertex subset");
  for (i = 0; i < numVertPerDim[0]; i++) {
    rc = fc_addMemberToSubset(vertSubset, i);
    fc_exitIfErrorPrintf(rc, "failed to add member to vertex subset");
  }
  rc = fc_createSubset(mesh, elem_subset_name, FC_AT_ELEMENT, &elemSubset);
  fc_exitIfErrorPrintf(rc, "failed to create element subset");
  for (i = 0; i < (numVertPerDim[0]-1)/2; i++) {
    for (j = 0; j < (numVertPerDim[0]-1)/2; j++) {
      rc = fc_addMemberToSubset(elemSubset, i*(numVertPerDim[0]-1)+j);
      fc_exitIfErrorPrintf(rc, "failed to add member to element subset");
    }
  }

  // Some reporting
  //printf("'%s':'%s':'%s'\n", dataset_file_name, mesh_name, var_name);
  //printf("%d X %d quad elements\n", numVertPerDim[0]-1, numVertPerDim[1]-1);
  //printf("origin = [ %g, %g ]\n", origin[0], origin[1]);
  //printf("extent = [ %g, %g ]\n", extent[0], extent[1]);
  //printf("numPoint = %d\n", numPt);
  //printf("         Mu-x     Mu-y    sigma  scale_0 scale_rate\n");
  //for (i = 0; i < numPt; i++) {
  //printf("%3d: %8g %8g %8g %8g %8g\n", i, centers[i][0], centers[i][1], 
  //   sigmas[i], initial_scales[i], rates[i]);
  //}

  // Write the dataset 
  rc = fc_writeDataset(dataset, dataset_file_name_exo, FC_FT_EXODUS);
  fc_exitIfErrorPrintf(rc, "failed to write dataset");
  
  rc = fc_finalLibrary();
  fc_exitIfErrorPrintf(rc, "Failed to finalize fc library");

  exit(0);
}
