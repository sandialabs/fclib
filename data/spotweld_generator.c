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
 * \file spotweld_generator.c
 * \brief  A simple program to create datasets for verification testing
 *    of analyzeSpotWelds.
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/data/spotweld_generator.c,v $
 * $Revision: 1.21 $
 * $Date: 2006/08/30 19:20:00 $
 *
 * \description
 *
 *   - Case 0: no translation
 *   - Case 1: pure z translation (normal)
 *   - Case 2: pure x translation (tangential)
 *   - Case 3: arbitrary translation I
 *   - Case 4: arbitrary translation II
 *
 * A spotweld is defined as a connection between a vertex on one surface
 * and a face on another surface (at a specific attachement point).
 *
 * Two sets of displaced coords are created. One with rotation and one without
 * For each case, a 'bar' is created connecting the attachement point on the 
 * quad and the vertex to make viewing easier.
 
 * \modifications
 *  - 2003-OCT-16 WSK Created.
 *  - 06/22/2004 WSK Changed to write nodeset and sideset as subsets
 *    instead of sets.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fc.h>

#define MAX_STR_LEN 1024
#define NUM_CASE 5
#define NUM_MESH 3
#define NUM_SUBSET 2
#define NUM_VAR 2
#define NUM_STEP 16

// Helper routine to calculate arbirary rotation of coordinates
static void rotate(double* x_p, double* y_p, double* z_p) {
  double x = *x_p, y = *y_p, z = *z_p;
  double new_x, new_y, new_z;

  // FIX calculate these from degrees using PI = d*pi/180
  // arbitrary rotations about each axis:
  double THETA_X = 0.366519143; // 21 degrees
  double THETA_Y = 0.296705973; // 17 degrees
  double THETA_Z = 0.191986218; // 11 degrees
  double sinx = sin(THETA_X), siny = sin(THETA_Y), sinz = sin(THETA_Z);
  double cosx = cos(THETA_X), cosy = cos(THETA_Y), cosz = cos(THETA_Z);

  new_x = x*cosy*cosz - y*cosy*sinz + z*siny;
  new_y = x*(sinx*siny*cosz + cosx*sinz) + y*(-sinx*siny*sinz + cosx*cosz) -
    z*sinx*cosy;
  new_z = x*(-cosx*siny*cosz + sinx*sinz) + y*(cosx*siny*sinz + sinx*cosz) +
    z*cosx*cosy;

  *x_p = new_x;
  *y_p = new_y;
  *z_p = new_z;
}

int main(void)
{
  int i, j, k, m, case_i;  // iterators
  FC_ReturnCode rc;
  int numStep = NUM_STEP, temp_numStep;  // number of time steps
  int numCase = NUM_CASE;  // number of spot welds
  int numMesh = NUM_MESH;    // number of set types, order is: node, face, bar
  int numDispl = NUM_VAR; // number of types of displacment vars
  char dbFileName_exo[MAX_STR_LEN] = "gen_spotweldmodel.ex2";
  char dbName[MAX_STR_LEN] = "spotweldmodel";
  char setNames[NUM_MESH][MAX_STR_LEN] = { "Nodes", "Surfaces", "Connectors" };
  char nodeNames[NUM_CASE][MAX_STR_LEN] = { "nodelist_0",
                                            "nodelist_1", 
                                            "nodelist_2",
                                            "nodelist_3",
                                            "nodelist_4" };
  char faceNames[NUM_CASE][MAX_STR_LEN] = { "surface_hex8_quadface4_0",
                                            "surface_hex8_quadface4_1",
                                            "surface_hex8_quadface4_2",
                                            "surface_hex8_quadface4_3",
                                            "surface_hex8_quadface4_4" };
  // FIX analyze spot welds expects the names to be postpended with '_vec'
  //char displNames[NUM_VAR][MAX_STR_LEN] = { "displ", "displ_rot" };
  char displNames[NUM_VAR][MAX_STR_LEN] = { "displ", // just trans
                                              "displ_rot" }; // rotated
  FC_ElementType celltype[NUM_MESH] = { FC_ET_POINT,
				       FC_ET_QUAD,
				       FC_ET_LINE };
  int nodeConns[5]; // nodes set connectivity
  int faceConns[4*5]; // faces set connectivity
  int barConns[2*5]; // bar connectivity
  double refNodeCoords[3] =     { 4.4, 4.5, 0. }; // node coordinates
  double refFaceCoords[4 * 3]  = { 2.8, 1.,  0.,
                                   9.3, 2.5, 0.,
                                   9.,  8.,  0.,
                                   0.4, 3.5, 0. }; // face coords in quadrant
  double refFaceCoords0[4 * 3] = {  1.,  1., 0.,
                                   -1.,  1., 0.,
                                   -1., -1., 0.,
                                    1., -1., 0. }; // face coords for origin
  double nodeCoords[NUM_CASE*3]; // for each case, node coords
  double faceCoords[NUM_CASE*4*3]; // for each case, face's vertex coords
  double barCoords[NUM_CASE*2*3]; // for each case, bar coords
  double newNodeCoords[3];  // buffer for node displaced & rotated coordinates
  double nodeDispls[NUM_CASE*3]; // buffer for node displacements 
  double newFaceCoords[4*3]; // buffer for face displaced & rotated coords
  double faceDispls[NUM_CASE*4*3]; // buffer for face displacements
  double newBarCoords[2*3]; // buffer for bar displaced & rotated coords
  double barDispls[NUM_CASE*2*3];  // buff for bar displacements
  float times[NUM_STEP];
  FC_Dataset dataset;
  FC_Mesh meshes[NUM_MESH]; // set order is: node, face, bar
  FC_Subset temp_subset;
  FC_Variable *seqVars[NUM_MESH][NUM_VAR]; 
  FC_Sequence sequence;

  // temp test rotation
  //double x = 1, y = 1, z = 1;
  //for (i = 0; i < 5; i++) {
  //rotate(&x, &y, &z);
  //printf ("%g %g %g\n", x, y, z);
  //}

  //---------------------------------------------------------------------
  //  Initialize library
  //---------------------------------------------------------------------
  rc = fc_initLibrary();
  fc_exitIfErrorPrintf(rc, "problem initalizing library");
  fc_setLibraryVerbosity(FC_WARNING_MESSAGES);
  
  //-------------------------------------------------------------------
  //  Create new ex2 dataset
  //-------------------------------------------------------------------
  rc = fc_createDataset(dbName, &dataset);
  fc_exitIfErrorPrintf(rc, "Error: could not create dataset '%s'",
		       dbFileName_exo);

  //-----------------------------------------------------------------  
  // Set up initial conns for all meshes
  //----------------------------------------------------------------- 
  for (i = 0; i < 1*numCase; i++)
    nodeConns[i] = i;
  for (i = 0; i < 4*numCase; i++)
    faceConns[i] = i;
  for (i = 0; i < 2*numCase; i++)
    barConns[i] = i;

  //-----------------------------------------------------------------  
  // Set up initial coords for all meshes
  //----------------------------------------------------------------- 
  // Case 0 - copy reference coords 0
  for (i = 0; i < 3; i++) 
    nodeCoords[i] = 0; // the origin
  memcpy(&faceCoords[0], refFaceCoords0, sizeof(double)*3*4);
  
  // Put a case in each quadrant of the XY plane
  // Copy reference coords for remaining cases
  for (i = 1; i < numCase; i++) {
    memcpy(&nodeCoords[i*3], refNodeCoords, sizeof(double)*3);
    memcpy(&faceCoords[i*3*4], refFaceCoords, sizeof(double)*3*4);
  }
  // Case 1 - leave in 1st quadrant.
  // Case 2 - make x coords negative
  nodeCoords[2*3] *= -1.;
  for (i = 0; i < 4; i++) 
    faceCoords[2*4*3 + i*3] *= -1.;
  // Case 3 - make y coords negative
  nodeCoords[3*3 + 1] *= -1.;
  for (i = 0; i < 4; i++)
    faceCoords[3*4*3 + i*3 + 1] *= -1.;
  // Case 4 - make x & y coords negative
  nodeCoords[4*3] *= -1.;
  nodeCoords[4*3 + 1] *= -1.;
  for (i = 0; i < 4; i++) {
    faceCoords[4*4*3 + i*3] *= -1.;
    faceCoords[4*4*3 + i*3 + 1] *= -1.;
  }

  // initially, ends of bars are the same as the node coords
  for (case_i = 0; case_i < numCase; case_i++) {
    for (i = 0; i < 2; i++)
      for (j = 0; j < 3; j++)
        barCoords[case_i*2*3 + i*3 + j] = nodeCoords[case_i*3+j];
  }

  //-----------------------------------------------------------------  
  // Create Meshes for the nodes, faces and bars
  //----------------------------------------------------------------- 
  {
    int numComp = 3;
    int numElem = numCase;
    int numVertPerElem[NUM_MESH] = { 1, 4, 2 };
    int* conns[NUM_MESH] = { nodeConns, faceConns, barConns };
    double* coords[NUM_MESH] = { nodeCoords, faceCoords, barCoords};
    
    // create meshes
    for (j = 0; j < numMesh; j++) {
      rc = fc_createMesh(dataset, setNames[j], &meshes[j]);
      fc_exitIfErrorPrintf(rc, "Failed to create mesh");
      rc = fc_setMeshCoords(meshes[j], numComp, numVertPerElem[j]*numElem, 
			    coords[j]);
      fc_exitIfErrorPrintf(rc, "Failed to set mesh coords");
      rc = fc_setMeshElementConns(meshes[j], celltype[j], numElem, conns[j]);
      fc_exitIfErrorPrintf(rc, "Failed to set mesh conns");
    }
  }
 
  //-----------------------------------------------------------------  
  // Create subsets for the nodes and faces for each case
  //-----------------------------------------------------------------
  for (case_i = 0; case_i < numCase; case_i++) {
    rc = fc_createSubset(meshes[0], nodeNames[case_i], celltype[0],
			 &temp_subset);
    fc_exitIfErrorPrintf(rc, "failed to create node subset");
    rc = fc_addMemberToSubset(temp_subset, case_i);
    fc_exitIfErrorPrintf(rc, "failed to add member to node subset");
    rc = fc_createSubset(meshes[1], faceNames[case_i], celltype[1],
			 &temp_subset);
    fc_exitIfErrorPrintf(rc, "failed to create face subset");
    rc = fc_addMemberToSubset(temp_subset, case_i);
    fc_exitIfErrorPrintf(rc, "failed to add member to face subset");
  }
 
  //-----------------------------------------------------------------  
  // Create the sequence & seq variables
  //----------------------------------------------------------------- 
  for (i = 0; i < numStep; i++)
    times[i] = i;
  rc = fc_createSequence(dataset, "time", &sequence);
  fc_exitIfErrorPrintf(rc, "failed to create sequence");
  rc = fc_setSequenceCoords(sequence, numStep, FC_DT_FLOAT, times);
  fc_exitIfErrorPrintf(rc, "failed to set time values");
  for (i = 0; i < numMesh; i++) {
    for (j = 0; j < numDispl; j++) {
      rc = fc_createSeqVariable(meshes[i], sequence, displNames[j], &temp_numStep,
				&seqVars[i][j]); 
      fc_exitIfErrorPrintf(rc, "failed to creat seq variable");
    }
  }

  //-----------------------------------------------------------------  
  // Calculate displacements & fill in seq variables
  // (note that 2nd point of bars has same coordinates as the node)
  //----------------------------------------------------------------- 
  for (i = 0; i < numStep; i++) {
    for (j = 0; j < numDispl; j++) { // do both displ & displ_rot
      
      for (case_i = 0; case_i < numCase; case_i++) {
	
	// start with initial coords
	// note, only need to do 1st of bar coords, 2nd will be same as node
	memcpy(newNodeCoords, &nodeCoords[case_i*3], sizeof(double)*3);
	memcpy(newBarCoords, &nodeCoords[case_i*3], sizeof(double)*3);
	memcpy(newFaceCoords, &faceCoords[case_i*3*4], sizeof(double)*3*4);
	
	// Node translations (don't translate quads or bars)
	if (case_i == 0) {
	  // do nothing
	}
	if (case_i == 1) {
	  // just change z 
	  newNodeCoords[2] += i;
	}
	else if (case_i == 2) {
	  // just change x 
	  newNodeCoords[0] += i;
	}
	else if (case_i == 3) {
	  // change x, y & z the same amount
	  for (k = 0; k < 3; k++)
	    newNodeCoords[k] += i;
	}
	else if (case_i == 4) {
	  // arbitrary data based on i
	  newNodeCoords[0] += 0.2*i;
	  newNodeCoords[1] += 0.5*i;
	  newNodeCoords[2] += 0.7*i;
	}
	
        if (case_i == 0) {
          for (k = 0; k < 3; k++) {
            nodeDispls[k] = 0.;
            for (m = 0; m < 2; m++)
              barDispls[m*3+k] = 0.;
            for (m = 0; m < 4; m++)
              faceDispls[m*3+k] = 0.;
          }
        }
        else { // all cases but 0 are translated
          // 2nd time through, rotate i times
          if (j == 1) {
            for (k = 0; k < i; k++) {
              rotate(&newNodeCoords[0], &newNodeCoords[1], &newNodeCoords[2]);
              rotate(&newBarCoords[0], &newBarCoords[1], &newBarCoords[2]);
              for (m = 0; m < 4; m++)
                rotate(&newFaceCoords[m*3], &newFaceCoords[m*3+1],
                       &newFaceCoords[m*3+2]);
            }
          }
          
          // calculate the displacements
          for (k = 0; k < 3; k++) {
            nodeDispls[case_i*3+k] = newNodeCoords[k] - nodeCoords[case_i*3+k];
            barDispls[case_i*2*3+k] = newBarCoords[k] - nodeCoords[case_i*3+k];
            barDispls[case_i*2*3+3+k] = nodeDispls[k];
            for (m = 0; m < 4; m++)
              faceDispls[case_i*4*3 + m*3 + k] = newFaceCoords[m*3+k] - 
                faceCoords[case_i*4*3 + m*3 + k];
          }
        }
      }

      // create displacement vars on each set
      fc_setVariableData(seqVars[0][j][i], 5*1, 3, FC_AT_VERTEX,
			 FC_MT_VECTOR, FC_DT_DOUBLE, nodeDispls);
      fc_setVariableData(seqVars[1][j][i], 5*4, 3, FC_AT_VERTEX,
			 FC_MT_VECTOR, FC_DT_DOUBLE, faceDispls);
      fc_setVariableData(seqVars[2][j][i], 5*2, 3, FC_AT_VERTEX,
			 FC_MT_VECTOR, FC_DT_DOUBLE, barDispls);
    }
  }

  //-------------------------------------------------------------------
  // Write Dataset & Finalize access to the library
  //-------------------------------------------------------------------
  rc = fc_writeDataset(dataset, dbFileName_exo, FC_FT_EXODUS);
  fc_exitIfErrorPrintf(rc, "failed to write dataset");
  
  rc = fc_finalLibrary();
  fc_exitIfErrorPrintf(rc, "Failed to finalize fc library");

  //-----------------------------------------------------------------  
  // cleanup & exit
  //-----------------------------------------------------------------
      // FIX? need we clean up info arrays?
  exit(0);
}
