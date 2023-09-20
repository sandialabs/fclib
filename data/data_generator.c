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
 * \file data_generator.c 
 * \brief  Generates a selected type of data file.
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/data/data_generator.c,v $
 * $Revision: 1.19 $
 * $Date: 2006/08/30 19:20:00 $
 * 
 *  \description
 *
 *    This program was created to generate test data.  See the documentation
 *    for saf_generator.c for the format of the data file produced and the
 *    formats of the input files.
 *
 *    This program reads in a vertex file (and optionally an element file) and
 *    generates a requested type of data file. The choices are 1) scalar
 *    temperature (per vertex or per element), 2) displacement vector (3D, per
 *    vertex), or 3) symmetric tensor stress (3 components, per element).
 *    
 *    The scalar temperature field has the form of the sum of two gaussians. If
 *    only the vertex file is provided, the values will be per vertex. If the
 *    element file is also provided the values will be for each elemement.
 *    The user can also specify whether the gaussians are 2D (using only
 *    the x and y coords) or 3D.
 *
 *    The displacement vector are the values of the vertex coordinates. This
 *    means that if the displacement vector is applied to the coordinates, the
 *    mesh will appear to have doubled in size.
 *
 *    The stress tensor components will also be the values of the vertex
 *    coordinates or element centroids, and probably have no physical meaning.
 *
 *    For the temperature calculations, the program will have problems if all
 *    the vertices lie in a single plane or line or point. It automatically
 *    exits if fewer than four vertices are specified.  Details: the
 *    coordinates are normalized to a domain ranging form 0 to 1 in each
 *    axis. Then the gaussians are calculated from with predetermined
 *    parameters that looked "good" to me. That is there are two easily
 *    identifiable humps. If the domain is not square, the humps will be
 *    distorted.
 */
  
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <fc.h>
#include <fcP.h>

// Assumptions
// coordinates will be in 3D

// Global Variables
// these variables identify the gaussian information
double mu1[3] = { 0.3, 0.7, 0.3 };
double mu2[3] = { 0.8, 0.2, 0.8 };
double sigma1 = 0.2;
double sigma2 = 0.1;
double base = 60;
double scale1 = 7.746;  // scale per dimension
double scale2 = 5.4772;
//double scale1 = 60;  // scale for 2d
//double scale2 = 30;

// functions
double normal_dist(double x, double mu, double sigma);
double calc_data(int numDim, double* X);

int main(int argc, char **argv)
{
  FC_ReturnCode rc;
  int i, j, k;
  int numVert, numElem, numDim, numDataPoint, numComp, temp_num;
  int gauss_numDim;
  FC_ElementType elemType;
  int doElement = 0;   // 1 means data is per element, 0 data is per vertex
  int doTemp = 0;     
  int doDispl = 0;
  int doStress = 0;
  double* coords, *centroids, *data;
  int* conns;
  char* data_file_name = NULL;
  char* vert_file_name = NULL;
  char* elem_file_name = NULL;
  char* data_name, *assoc_name, *mathtype_name;
  char* temp_name;   // Temp string
  FILE* data_file;

  //------------------------------------------------------------------
  //  argument handling
  //------------------------------------------------------------------
  // If no arguments, print usage and exit
  if (argc < 4) {
  usage:
    printf("usage: %s <temp2D|temp3D|displ|stress> -v mesh.vert\n", argv[0]);
    printf("\n");
    printf("Where the first is argument is the type of data to generate (temperature \n");
    printf("based on 2D gaussians, temperature based on 3D gaussians, displacement vector,\n");
    printf("or stress symmetric tensor), and mesh.vert is the file containing the vertex\n");
    printf("coordinates. The displacement vector and the 3 component stress tensor\n");
    printf("are just set equal to the coordinates of the vertices, or the centroids of the\n");
    printf("elements (displacement vector will not be created per element).\n");
    printf("\n");
    printf("options:\n");
    printf("   -o mesh.data       : output filename (default is stdout)\n");
    printf("   -e mesh.elem       : force the data to be per elem (default is per vertex\n");
    printf("                        and provide the connectivity file.\n");
    exit(0);
  }
  for (i = 1; i < argc; i++) {

    if ((!strcmp(argv[i], "help")) ||
        (!strcmp(argv[i], "-h")) ||
        (!strcmp(argv[i], "-help")) ||
        (!strcmp(argv[i], "-")) ||
        (!strcmp(argv[i], "--"))) 
      goto usage;  

    else if (!strcmp(argv[i], "-o")) {
      if (i+1 >= argc || argv[i+1][0] == '-') goto usage;
      data_file_name = argv[i+1];
      i++;
    }

    else if (!strcmp(argv[i], "-v")) {
      if (i+1 >= argc || argv[i+1][0] == '-') goto usage;
      vert_file_name = argv[i+1];
      i++;
    }

    else if (!strcmp(argv[i], "-e")) {
      if (i+1 >= argc || argv[i+1][0] == '-') goto usage;
      elem_file_name = argv[i+1];
      doElement = 1;
      i++;
    }

    else {
      if (!strcmp(argv[i], "temp2D")) {
	gauss_numDim = 2;
	doTemp = 1;
      }
      else if (!strcmp(argv[i], "temp3D")) {
	gauss_numDim = 3;
	doTemp = 1;
      }
      else if (!strcmp(argv[i], "displ"))
	doDispl = 1;
      else if (!strcmp(argv[i], "stress"))
	doStress = 1;
      else
	goto usage;
    }
  }

  // check that required options were processed
  if ( !data_file_name || !vert_file_name || (doElement && !elem_file_name) ||
       (doTemp + doDispl + doStress != 1) || (doDispl && doElement) )
    goto usage;

  //------------------------------------------------------------------
  //  read vertex file
  //------------------------------------------------------------------
  rc = _fc_readVertexFile(vert_file_name, &numVert, &numDim, &coords);
  fc_exitIfErrorPrintf(rc, "failed to read vertex file"); 
  if (numVert < 4) 
    fc_exitIfErrorPrintf(FC_ERROR, "Vert file must have at least 4 verts\n");

  //------------------------------------------------------------------
  //  optionally read element file & create centroids
  //------------------------------------------------------------------
  if (doElement) {
    int numVperE;

    // read
    rc = _fc_readElementFile(elem_file_name, &temp_num, &numElem, &elemType, 
			     &temp_name, &conns);
    fc_exitIfErrorPrintf(rc, "failed to read element file"); 
    if (temp_num != numVert) 
      fc_exitIfErrorPrintf(FC_ERROR, 
			   "numVert in vertex and elem files do not agree");
    free(temp_name);

    // create centroids
    numVperE = fc_getElementTypeNumVertex(elemType);
    centroids = (double*)calloc(numElem*numDim, sizeof(double));
    for(i = 0; i < numElem; i++) {
      double point[3] = { 0., 0., 0 };
      for (j = 0; j < numVperE; j++) {
	int vertId = conns[i*numVperE + j];
	for (k = 0; k < numDim; k++)
	  point[k] += coords[vertId*numDim+k];
      }
      for (j = 0; j < numDim; j++)
	centroids[i*numDim+j] = point[j]/numVperE;
    }

    // cleanup
    free(conns);
  }

  //------------------------------------------------------------------
  //  set metadata & calculate the requested type of data ...
  //  (free big data as you go)
  //------------------------------------------------------------------
  if (doElement) {
    numDataPoint = numElem;
    assoc_name = "ELEMENT";  }
  else {
    numDataPoint = numVert;
    assoc_name = "VERTEX";
  }

  //------------------------------------------------------------------
  //  displacement vector data
  //------------------------------------------------------------------
  if (doDispl) {
    // metadata
    mathtype_name = "VECTOR";
    data_name = "displacement";
    numComp = numDim;
    // data
    data = coords;
  }
  
  //------------------------------------------------------------------
  //  stress symtensor data
  //------------------------------------------------------------------
  else if (doStress) {
    // metadata
    mathtype_name = "SYMTENSOR";
    data_name = "stress";
    numComp = 3;
    // data
    if (doElement) {
      data = centroids;
      free(coords);
    }
    else 
      data = coords;
  }  
  
  //------------------------------------------------------------------
  //  gaussian scalar data
  //------------------------------------------------------------------
  else{
    double lowerbounds[3], upperbounds[3], range[3];
    double* points;
    
    // metadata
    mathtype_name = "SCALAR";
    data_name = "temperature";
    numComp = 1;
    
    // get bounding box info of the verts
    for (j = 0; j < numDim; j++) {
      lowerbounds[j] = coords[j];
      upperbounds[j] = coords[j];
    }
    for (i = 1; i < numVert; i++) {
      for (j = 0; j < numDim; j++) {
	if (coords[i*numDim+j] < lowerbounds[j])
	  lowerbounds[j] = coords[i*numDim+j];
	else if (coords[i*numDim+j] > upperbounds[j])
	  upperbounds[j] = coords[i*numDim+j];
      }
    }
    for (j = 0; j < numDim; j++) 
      range[j] = upperbounds[j] - lowerbounds[j];

    // intermediate data - normalize points (numDim) to bounding box
    if (doElement) {
      points = centroids;
      free(coords);
    }
    else
      points = coords;
    for (i = 0; i < numDataPoint; i++) 
      for (j = 0; j < numDim; j++)
	points[i*numDim+j] = (points[i*numDim+j] - lowerbounds[j])/range[j];
    
    // calc data (numComp)
    data = (double*)malloc(numDataPoint*sizeof(double));
    if (!data)
      fc_exitIfError(FC_MEMORY_ERROR);
    for (i = 0; i < numDataPoint; i++) 
      data[i] = calc_data(gauss_numDim, &points[i*numDim]);
    free(points);
  }

  //------------------------------------------------------------------
  //  Write the data file
  //------------------------------------------------------------------
  data_file = fopen(data_file_name, "w");
  if (data_file == NULL) {
    printf("Can't open file %s\n", data_file_name);
    return -1;
  }
  fprintf(data_file, "%d\n", numDataPoint);
  fprintf(data_file, "%s\n", assoc_name);
  fprintf(data_file, "%s\n", mathtype_name);
  fprintf(data_file, "%s\n", data_name);
  for (i = 0; i < numDataPoint; i++) {
    fprintf(data_file, "%f", data[numComp*i]);
    for (j = 1; j < numComp; j++)
      fprintf(data_file, " %f", data[numComp*i+j]);
    fprintf(data_file, "\n");
  }  
  fclose(data_file);

  //------------------------------------------------------------------
  //  cleanup & exit
  //------------------------------------------------------------------
  free(data);
  exit(0);
}

// calculating sum of two normal distributions PLUS baseline
// Gaussian parameters are declared globally
double calc_data(int numDim, double* X) {
  int i;

  double answer1 = 1;
  double answer2 = 1;
  for (i = 0; i < numDim; i++) {
    answer1 *= scale1*normal_dist(X[i], mu1[i], sigma1);
    answer2 *= scale2*normal_dist(X[i], mu2[i], sigma2);
  }
 
  return base + answer1 + answer2;
}

double normal_dist(double x, double mu, double sigma) {
  return 1/sqrt(2.*M_PI*sigma) * exp(-1./2.*pow((x-mu)/sigma,2));
}
