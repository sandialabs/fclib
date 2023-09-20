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
 * \file elemdeath_ex1.c
 * \brief Example of using some of the elemdeath module and related calls
 *
 * \description
 *
 *   Example of using some of the elemdeath module and related calls.
 *   Loops through 4 different cases of dead elements in a 3x3 hex mesh
 *   and shows how to find various topological information about the
 *   dead element regions wrt the mesh sides.
 *
 *   The code prints out the side info of the original mesh, and info
 *   about 4 different dead region cases.
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/examples/elemdeath_ex1.c,v $
 * $Revision: 1.27 $
 * $Date: 2006/08/30 19:20:00 $
 *
 *  \todo
 *   - may replace this with the sidesDeath code, or that code applied
 *     to these meshes
 *   - rewrite to take advantage of new calc of sides in shape, once
 *     i decide what will happen to thin shape and if there will be a new shape struct
 *
 *  \modifications
 *   - 12/05/2005 ACG created from surfacedemo.c in my sandbox
 *   - 01/04/2006 ACG added new methods from shape. Got rid of eval method to
 *     clean this up.
 *   - 06/04/2006 ACG updated to use new getShapeSides interface
 *   - 06/16/2006 ACG cleaned up to use shape and remove thin shape
 */
#include "fc.h"

/**
 * \brief - makes dead elements in the mesh
 *
 * \description
 *  makes dead elements in the mesh
 */
static void makeDeadVar(
 FC_Mesh mesh, /**< input - mesh */
 int iteration, /**< input - which iteration this is */
 int numelempertype, /**< input - num element in this mesh */
 FC_Variable* deadvar /**< output - variable where val of 1 = dead elem */
){
  FC_ReturnCode rc;
  int numDataPoint,j;
  int *data;
 
  rc = fc_createVariable(mesh,"deadvar",deadvar);
  fc_exitIfErrorPrintf(rc, "failed to create deadvar");
    
  //now set dead data
  numDataPoint = numelempertype;
  data = (int*)malloc(numDataPoint*sizeof(int));
  switch(iteration){
  case 0: //empty set
    printf("CASE 0: Empty set. \n");
    for (j = 0; j < numDataPoint; j++)
      data[j] = 0;
    break;
  case 1: //complete set
    printf("CASE 1: Complete set. \n");
    for (j = 0; j < numDataPoint; j++)
      data[j] = 1;
    break;
  case 2:
    printf("CASE 2: Bottom along the corner, top across the middle. \n");
    for (j = 0; j < numDataPoint; j++)
      data[j] = 0;
    //known test case
    data[0] = 1;
    data[1] = 1;
    data[2] = 1;
    data[7] = 1;
    data[16] = 1;
    data[25] = 1;
    break;
  case 3:
    printf("CASE 3: Through the middle\n");
    for (j = 0; j < numDataPoint; j++)
      data[j] = 0;
    //known test case
    data[4] = 1;
    data[13] = 1;
    data[22] = 1;
    break;
  default: //shouldnt occur
    break;
  }
  fflush(NULL);
  
  //set var data
  rc = fc_setVariableData(*deadvar,numDataPoint,1,FC_AT_ELEMENT,
			  FC_MT_SCALAR, FC_DT_INT,data);
  fc_exitIfErrorPrintf(rc,"failed to set var data");
  free(data);
};


/**
 * \brief main
 *
 * \description
 *  Example of using some of the elemdeath module and related calls.
 *  Loops through 4 different cases of dead elements in a 3x3 hex mesh
 *  and shows how to find various topological information about the
 *  dead element regions wrt the mesh sides.
 *
 */
int main(void){
  FC_ReturnCode rc;
  int i,j;
  FC_Dataset dataset;
  FC_Mesh mesh;

  //coord info for the mesh
  FC_Coords lowers = { 0., 0., 0. };
  FC_Coords uppers = { 1., 1., 0.5 };
  int numElemPerType = 27;

  FC_VerbosityLevel verbose_level = FC_QUIET;

  int numMeshShapes;
  FC_Shape* meshShapes;
  FC_Subset meshskin;

   // init library and load dataset
  rc = fc_setLibraryVerbosity(verbose_level);
  fc_exitIfError(rc);
  rc = fc_initLibrary();
  fc_exitIfError(rc);

  // setup - create test mesh
  rc = fc_createDataset("dset", &dataset);
  fc_exitIfErrorPrintf(rc,"abort: failed to create dataset");

  rc = fc_createSimpleHexMesh(dataset,"hex_mesh",3,3,3,lowers,uppers,&mesh);
  fc_exitIfErrorPrintf(rc, "abort: failed to creat mesh");

  rc = fc_getMeshShapes(mesh,30,0,&numMeshShapes,
			    &meshShapes);
  fc_exitIfErrorPrintf(rc, "failed to getMeshShapes");

  for (i = 0; i < numMeshShapes; i++){
    fc_printShape(&meshShapes[i],"",0);
  }
  //we know there is only 1 mesh shape
  if (numMeshShapes != 1){
    fc_printfErrorMessage("There should only be 1 mesh shape");
    exit(-1);
  }

  //now get the skin (could also get this by merging the meshSideFaces)
  //becuase we will be looking at decayed mesh skin segments
  rc= fc_getMeshSkin(mesh,&meshskin);
  if (rc != FC_SUCCESS){
    fc_printfErrorMessage("cant get mesh skin");
    exit(1);
  }

  printf("-------------------------------------------------------------\n");
  
  //now try 4 different dead cases:
  for (i = 0 ; i < 4; i++){
    int numdeadregions;
    FC_Variable deadvar;
    FC_Subset deadsubset,*deadregions; 

    //set the dead elements
    makeDeadVar(mesh,i,numElemPerType,&deadvar);
    
    //now get dead subset
    rc = fc_createThresholdSubset(deadvar,"=",1,"temp_dead", &deadsubset);
    fc_exitIfErrorPrintf(rc, "failed to create dead subset");
    fc_deleteVariable(deadvar);
    
    //segment it by connectivity in case its not contiguous
    //to get the dead regions
    rc = fc_segment(deadsubset,1,&numdeadregions,&deadregions);
    fc_exitIfErrorPrintf(rc, "failed to segment dead subset");
    fc_deleteSubset(deadsubset);
    
    if (numdeadregions  == 1){
      printf("There is 1 dead region:\n");
    }else {
      printf("There are %d dead regions:\n",numdeadregions);
    }
    fflush(NULL);

    //for each dead region, see what surfaces of each meshshape
    //it intersects with
    for (j = 0; j < numdeadregions; j++){
      FC_Subset *decayedMeshSkin_seg; // segments the intersection of
                                      // the dead region with the mesh skin
      int numdecayedseg; //number of the above
      FC_Subset deadskin;
      int** decayflags;
      int l,ll,lll,llll;

      //note - may change later what im doing here, in face maybe
      //shapeDeath may become this example
      //segment the intersection of the dead region with the mesh skin
      rc = fc_getDecayedSkin(deadregions[j],&meshskin,&deadskin);
      fc_exitIfErrorPrintf(rc, "failed to get decayed skin");
      if (FC_HANDLE_EQUIV(deadskin,FC_NULL_SUBSET)){
	numdecayedseg = 0;
      }else{
	rc = fc_segment(deadskin,1,&numdecayedseg,&decayedMeshSkin_seg);
	fc_exitIfErrorPrintf(rc, "failed to get decayed mesh skin segments");
	fc_deleteSubset(deadskin);
      }

      decayflags = (int**)malloc(numdecayedseg*sizeof(int*));
      if (decayflags == NULL){
	fc_printfErrorMessage("MEMORY ERROR\n");
	exit(-1);
      }
				 
      //for each dead segment get which mesh part sides it intersects
      printf("DeadRegion %d intersects the Mesh surface in ",j);
      printf("%d decayedMeshSkin segments\n", numdecayedseg);
      for (l = 0; l < numdecayedseg; l++){
	int count;
	//	fc_printSubset(decayedMeshSkin_seg[l],"",0);
	rc = fc_getShapeSidesDecayType(decayedMeshSkin_seg[l],
				       &meshShapes[0],
			      &decayflags[l]);
	fc_exitIfErrorPrintf(rc, "failed to get sides decay");
	
	count = 0;
	printf("Segment %d sides decay: ", l);
	for(ll = 0; ll < meshShapes[0].numSides;ll++){
	  if (decayflags[l][ll] == 1 ){
	    printf("%d ",ll);
	    count++;
	  }
	}
	if (!count) printf("none");
	printf("\n");

	count = 0;
	printf("Segment %d sides erosion: ", l);
	for(ll = 0; ll < meshShapes[0].numSides;ll++){
	  if (decayflags[l][ll] == 2){
	    printf("%d ",ll);
	    count++;
	  }
	}
	if (!count) printf("none");
	printf("\n");
      }

      for (l = 0; l < numdecayedseg; l++){
	fc_deleteSubset(decayedMeshSkin_seg[l]);
      }
      free(decayedMeshSkin_seg);

      //check the sides decay flag against the side adjacencies
      for (l = 0; l < numdecayedseg; l++){
	for (ll = l+1; ll < numdecayedseg ; ll++){
	  for(lll = 0; lll < meshShapes[0].numSides; lll++){
	    for(llll = 0; llll < meshShapes[0].numSides; llll++){
	      if (lll != llll){
		// printf("comparing flags[%d][%d] flags[%d][%d]\n",
		//      l,lll,ll,llll);
		if (decayflags[l][lll] == 1  && decayflags[ll][llll] == 1 &&
		    meshShapes[l].adjmatrix[lll][llll] == 0){
		  printf("Deadregion %d intersects the mesh in non-adjacent sides: %d %d\n",j,lll,llll);
		}
	      }
	    }
	  }
	}
      }
      
      //clean up
      for(l = 0; l < numdecayedseg; l++){
	free(decayflags[l]);
      }
      free(decayflags);

      fc_deleteSubset(deadregions[j]);
    } // j - dead regions
    printf("\n");
    free(deadregions);
  }//k - cases

  //clean up
  for(i = 0; i < numMeshShapes; i++){
    fc_freeShape(&meshShapes[i]);
  }
  free(meshShapes);

  fc_deleteMesh(mesh);
  fc_deleteSubset(meshskin);
  
  return 1;
}

