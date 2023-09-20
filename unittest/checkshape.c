/* 
   Copyright (2000) Sandia Corporation. Under the terms of Contract
  DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
  certain rights in this software.
 
  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are met:
 
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the
      distribution.

    * Neither the name of Sandia nor the names of any contributors may
      be used to endorse or promote products derived from this software
      without specific prior written permission.
  
    * Modified source versions must be plainly marked as such, and must
      not be misrepresented as being the original software.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
  ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHORS OR CONTRIBUTORS BE LIABLE FOR 
  ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT 
  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
  OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH 
  DAMAGE.
*/

/**
 * \file checkshape.c
 * \brief Unit testing of \ref Shape module
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/unittest/checkshape.c,v $
 * $Revision: 1.82 $
 * $Date: 2006/09/29 22:03:42 $
 *
 * \description
 *
 *    Tests the shape routines.
 *
 * \modifications
 *    - 12/19/05 ACG created. moved some tests here from topo
 *    - 02/08/06 ACG reorganized test for matrix and simple becuase
 *      of logical flaw in old adjacency matrix calculation
 *    - 06/04/06 ACG commenting everything out while the redo of shape
 *       sides is going on
 *    - 06/23/06 ACG as part of redo of shpae, things are rather complicated.
 *      I am splitting up some checks becuase of this.
 *      \ref _fc_getSubsetSides is now a private function and that
 *      info will be obtained thru the wrapper functions 
 *      \ref fc_getSubsetShape and \ref fc_getMeshShape.
 *      The bad tests of these are in the the ShapeSides-bad test below.
 *      The good tests of these are in the the ShapeMatrix test below.
 *      Good tests of subsetSides are not done explicitly, because
 *      that function is good if the Shape functions are good.
 */

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <check.h>
#include "fc.h"
#include "fcP.h"
#include "shapeP.h"
#include "checkall.h"

// **** global variables


// **** test fixtures

static void shape_setup(void) {
  FC_ReturnCode rc;
  if (isForking) {
    fc_setLibraryVerbosity(fc_messages);
    rc = fc_initLibrary();
    fail_unless(rc == FC_SUCCESS, "test aborted: failed to init library");
  }
}

static void shape_teardown(void) {
  FC_ReturnCode rc;
  if (isForking) {
    rc = fc_finalLibrary();
    fail_unless(rc == FC_SUCCESS, "failed to final library at end of test");
  }
}


//this tests shape helper functions
START_TEST(helpers)
{
  FC_ReturnCode rc;
  FC_Dataset dataset;
  FC_Mesh mesh;

  FC_Coords lowers = { 0., 0., 0. };
  FC_Coords uppers = { 1., 0.75, 0.5 };

  FC_Shape *shapes,copyshape,copyshape2;
  FC_Subset *holdfaces, *holdelems;
  int numSides;
  int numShapes;

  int i,j;

  // setup - create test dataset
  rc = fc_createDataset("temp.xxx", &dataset);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create dataset");

  
  rc = fc_createSimpleHexMesh(dataset, "junk", 2, 2, 3, lowers, uppers,
			      &mesh);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create simple hex mesh");

  rc = fc_getMeshShapes(mesh,80,0,&numShapes,&shapes);
  fail_unless(rc == FC_SUCCESS, "abort: failed to get MeshSides");
  fail_unless(numShapes == 1, "abort: bad num of shapes");
  fail_unless(shapes[0].numSides == 6, "abort: bad num of MeshSides");

  //---copy shape ----
  
  //make a copy of the shape
  rc = fc_copyShape(&shapes[0],&copyshape);
  fail_unless(rc == FC_SUCCESS, "cant copy shape");
  fail_unless(copyshape.numSides == shapes[0].numSides, 
	      "wrong num of sides");
  //minimal check for subsets since it calls copy
  for (i = 0; i <copyshape.numSides; i++){
    int num1, num2;
    rc = fc_getSubsetNumMember(copyshape.faces[i],&num1);
    fail_unless(rc == FC_SUCCESS, "cant get subset nummembers");
    rc = fc_getSubsetNumMember(shapes[0].faces[i],&num2);
    fail_unless(rc == FC_SUCCESS, "cant get subset nummembers");
    fail_unless(num1 == num2, "bad shape copy");
    rc = fc_getSubsetNumMember(copyshape.elems[i],&num1);
    fail_unless(rc == FC_SUCCESS, "cant get subset nummembers");
    rc = fc_getSubsetNumMember(shapes[0].elems[i],&num2);
    fail_unless(rc == FC_SUCCESS, "cant get subset nummembers");
    fail_unless(num1 == num2, "bad shape copy");
    for (j = 0; j < copyshape.numSides; j++){
      fail_unless(copyshape.adjmatrix[i][j] ==
		  shapes[0].adjmatrix[i][j], 
		  "bad adj matrix copy");
    }
  }

  // --- init shape ---
  // have to be careful how to test this since i lose the handles to things
  rc = fc_initShape(&copyshape2);
  fail_unless(rc == FC_SUCCESS, "can't init shape");
  fail_unless(copyshape2.numSides == 0, "bad val for init shape");
  fail_unless(copyshape2.faces == NULL, "bad val for init shape");
  fail_unless(copyshape2.elems == NULL, "bad val for init shape");
  fail_unless(copyshape2.adjmatrix == NULL, "bad val for init shape");

  //give it some values
  copyshape2.numSides = copyshape.numSides;
  copyshape2.faces = copyshape.faces;
  copyshape2.elems = copyshape.elems;
  copyshape2.adjmatrix = copyshape.adjmatrix;
  rc = fc_initShape(&copyshape2);
  fail_unless(rc == FC_SUCCESS, "can't init shape");
  fail_unless(copyshape2.numSides == 0, "bad val for init shape");
  fail_unless(copyshape2.faces == NULL, "bad val for init shape");
  fail_unless(copyshape2.elems == NULL, "bad val for init shape");
  fail_unless(copyshape2.adjmatrix == NULL, "bad val for init shape");

  //make sure those vals werent dealloc
  fail_unless(copyshape.numSides == shapes[0].numSides,
	      "init other shape shouldnt affect the orig shape - dealloc problem");
  for (i = 0; i <copyshape.numSides; i++){
    int num1, num2;
    rc = fc_getSubsetNumMember(copyshape.faces[i],&num1);
    fail_unless(rc == FC_SUCCESS, "cant get subset nummembers - dealloc problem");
    rc = fc_getSubsetNumMember(shapes[0].faces[i],&num2);
    fail_unless(rc == FC_SUCCESS, "cant get subset nummembers");
    fail_unless(num1 == num2, "dealloc problem");
    rc = fc_getSubsetNumMember(copyshape.elems[i],&num1);
    fail_unless(rc == FC_SUCCESS, "cant get subset nummembers - dealloc problem ");
    rc = fc_getSubsetNumMember(shapes[0].elems[i],&num2);
    fail_unless(rc == FC_SUCCESS, "cant get subset nummembers");
    fail_unless(num1 == num2, "bad shape copy - dealloc problem");
    for (j = 0; j < copyshape.numSides; j++){
      fail_unless(copyshape.adjmatrix[i][j] ==
		  shapes[0].adjmatrix[i][j], 
		  "bad adj matrix - dealloc problem");
    }
  }


  //test clear while also doing the clear
  //if is deep copy, clearing the copy wont affect the orig
  numSides = copyshape.numSides;
  holdfaces = (FC_Subset*) malloc(numSides*sizeof(FC_Subset));
  holdelems = (FC_Subset*) malloc(numSides*sizeof(FC_Subset));
  fail_unless(holdfaces && holdelems, "memory error");
  for (j = 0; j < numSides; j++){
    holdfaces[j] = copyshape.faces[j];
    holdelems[j] = copyshape.elems[j];
  }
  rc = fc_clearShape(&copyshape);
  fail_unless(rc == FC_SUCCESS, "cant clear the shape");
  fail_unless(copyshape.numSides == 0, "bad return after clear shape");
  fail_unless(copyshape.adjmatrix == NULL, "bad return after clear shape");
  fail_unless(copyshape.faces == NULL, "bad return after clear shape");
  fail_unless(copyshape.elems == NULL, "bad return after clear shape");
  fail_unless(shapes[0].numSides != 0, "bad return after freeing copy - problem with deep copy ");
  fail_unless(shapes[0].adjmatrix != NULL, "bad return after free shape - problem with deep copy");
  fail_unless(shapes[0].faces != NULL, "bad return after free shape - problem with deep copy");
  fail_unless(shapes[0].elems != NULL, "bad return after free shape - problem with deep copy");
  for (j = 0; j < numSides; j++){
    fail_unless(fc_isSubsetValid(holdfaces[j]),
		"subset should still be valid");
    rc = fc_deleteSubset(holdfaces[j]);
    fail_unless(rc == FC_SUCCESS, "should be able to delete this subset");
    fail_unless(fc_isSubsetValid(shapes[0].faces[j]),
		"subset should still be valid - problem with deep copy");
    fail_unless(fc_isSubsetValid(holdelems[j]),
		"subset should still be valid");
    rc = fc_deleteSubset(holdelems[j]);
    fail_unless(rc == FC_SUCCESS, "should be able to delete this subset");
    fail_unless(fc_isSubsetValid(shapes[0].elems[j]),
		"subset should still be valid - problem with deep copy");
  }
  free(holdfaces);
  free(holdelems);
  fail_unless(rc == FC_SUCCESS, "cant free the copy");

  //bad args for copy 
  copyshape.numSides = 3;
  rc = fc_copyShape(NULL,&copyshape);
  fail_unless(rc != FC_SUCCESS, "should fail for NULL inputshape ");
  fail_unless(copyshape.numSides == 0, "shape should return as init for fail");
  fail_unless(copyshape.faces == NULL, "shape should return as init for fail");
  fail_unless(copyshape.elems == NULL, "shape should return as init for fail"); 
  fail_unless(copyshape.adjmatrix == NULL, "shape should return as init for fail"); 

  rc = fc_copyShape(&shapes[0],NULL);
  fail_unless(rc != FC_SUCCESS, "should fail for NULL return arg ");
  
  //bad args for print
  rc = fc_printShape(NULL," ",1);
  fail_unless(rc != FC_SUCCESS, "should fail for NULL inputshape ");
  
  //bad args for free
  rc = fc_freeShape(NULL);
  fail_unless(rc != FC_SUCCESS, "should fail for NULL inputshape ");
  
  //bad args for clear
  rc = fc_clearShape(NULL);
  fail_unless(rc != FC_SUCCESS, "should fail for NULL inputshape ");
  
  //bad args for init
  rc = fc_initShape(NULL);
  fail_unless(rc != FC_SUCCESS, "should fail for NULL inputshape ");

  
  //test free while also doing the free
  numSides = shapes[0].numSides;
  holdfaces = (FC_Subset*) malloc(numSides*sizeof(FC_Subset));
  holdelems = (FC_Subset*) malloc(numSides*sizeof(FC_Subset));
  fail_unless(holdfaces && holdelems, "memory error");
  for (j = 0; j < numSides; j++){
    holdfaces[j] = shapes[0].faces[j];
    holdelems[j] = shapes[0].elems[j];
  }

  rc = fc_freeShape(&shapes[0]);
  fail_unless(rc == FC_SUCCESS, "cant free the shape");
  fail_unless(shapes[0].numSides == 0, "bad return after free shape");
  fail_unless(shapes[0].adjmatrix == NULL, "bad return after free shape");
  fail_unless(shapes[0].faces == NULL, "bad return after free shape");
  fail_unless(shapes[0].elems == NULL, "bad return after free shape");
  for (j = 0; j < numSides; j++){
    fail_unless(!fc_isSubsetValid(holdfaces[j]),"subset should not be valid");
    fail_unless(!fc_isSubsetValid(holdelems[j]),"subset should not be valid");
  }
  
  free(holdfaces);
  free(holdelems);

  free(shapes);
}
END_TEST



// Simple assumes that fucntionality of sides, shapes and matrix work.
// This is testing of simple functions (used to be testing of simple shapes
// right now this tests 1) areas 2) create descending area order 
// 3) reordershape 4) largeandnonadjsides order 5) largeandopposingsides order
START_TEST(orders)
{
  FC_ReturnCode rc;
  FC_Dataset dataset;
  FC_Mesh mesh;

  FC_Coords lowers = { 0., 0., 0. };
  FC_Coords uppers = { 1., 0.75, 0.5 };

  FC_Shape *shapes;
  FC_Shape copyshape;
  FC_Subset holdsubset;
  int numShapes;
  FC_Subset badSubset={999,999}, emptyf, test;
  double checkareas[6];
  double *areas;
  int *neworder;

  int i,j,ii;

  // setup - create test dataset
  rc = fc_createDataset("temp.xxx", &dataset);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create dataset");

  //do it twice, once with 2 sides equal -- see change in uppers at bottom of loop
  for (ii = 0; ii < 2; ii++){
    // Make the hex mesh;
    rc = fc_createSimpleHexMesh(dataset, "junk", 2, 2, 3, lowers, uppers,
				&mesh);
    fail_unless(rc == FC_SUCCESS, "abort: failed to create simple hex mesh");
    
    rc = fc_getMeshShapes(mesh,80,0,&numShapes,&shapes);
    fail_unless(rc == FC_SUCCESS, "abort: failed to get MeshSides");
    fail_unless(numShapes == 1, "abort: bad num of shapes");
    fail_unless(shapes[0].numSides == 6, "abort: bad num of MeshSides");
    
    //----- areas ---------

    //get results. have to do this at runtime since the numbering
    //is arbitary
    for (i = 0 ; i < shapes[0].numSides; i++){
      FC_Coords cent1, cent2;
      int dim;
      checkareas[i] = -1;
      rc = fc_getSubsetCentroid(shapes[0].faces[i], &dim, &cent1);
      fail_unless(rc == FC_SUCCESS, "cant get centroid");
      for (j = 0; j < shapes[0].numSides; j++){
	if (i ==j) continue;
	rc = fc_getSubsetCentroid(shapes[0].faces[j], &dim, &cent2);
	fail_unless(rc == FC_SUCCESS, "cant get centroid");
	if (FC_DBL_EQUIV(cent1[0],cent2[0]) &&
	    FC_DBL_EQUIV(cent1[1],cent2[1]) &&
	    ! FC_DBL_EQUIV(cent1[2],cent2[2])){
	  checkareas[i] = uppers[0]*uppers[1];
	  break;
	}
	if (FC_DBL_EQUIV(cent1[0],cent2[0]) && 
	    FC_DBL_EQUIV(cent1[2],cent2[2]) &&
	    ! FC_DBL_EQUIV(cent1[1],cent2[1])){
	  checkareas[i] = uppers[0]*uppers[2];
	  break;
	}
	if (FC_DBL_EQUIV(cent1[1],cent2[1]) &&
	    FC_DBL_EQUIV(cent1[2],cent2[2]) &&
	    ! FC_DBL_EQUIV(cent1[0],cent2[0])){
	  checkareas[i] = uppers[1]*uppers[2];
	  break;
	}
      }
    }

    rc = fc_getShapeSidesAreas(&shapes[0], &areas);
    fail_unless(rc == FC_SUCCESS, "abort: failed to get areas");

    //check results
    for (i = 0; i <shapes[0].numSides; i++){
      fail_unless(FC_DBL_EQUIV(areas[i],checkareas[i]), "wrong areas");
    }

    //----- descending area order ---------

    rc = fc_createDescendingAreaOrder(&shapes[0],&neworder);
    fail_unless(rc == FC_SUCCESS, 
		"abort: failed to get descending area order");
    for (i = 0; i < shapes[0].numSides-1; i++){
      fail_unless(areas[neworder[i]]  >= areas[neworder[i+1]],
		  "bad descending order");
    }
    free(areas);

    
    //----- reorder ------------

    rc = fc_copyShape(&shapes[0],&copyshape);
    fail_unless(rc == FC_SUCCESS, "cant copy shape");

    rc = fc_reorderShape(neworder,&shapes[0]);
    fail_unless(rc == FC_SUCCESS, "cant reorder shape");

    for (i = 0; i <copyshape.numSides; i++){
      int *orig;
      int *reord;
      int numorig, numreord;
      rc = fc_getSubsetMembersAsArray(copyshape.faces[neworder[i]],
				      &numorig,&orig);
      fail_unless(rc == FC_SUCCESS, "cant get subset members");
      rc = fc_getSubsetMembersAsArray(shapes[0].faces[i],
				      &numreord,&reord);
      fail_unless(rc == FC_SUCCESS, "cant get subset members");
      fail_unless(numorig == numreord, "bad subset reord");
      for (j = 0; j < numorig; j++){
	fail_unless(orig[j] == reord[j], "bad subset reord");
      }
      free(orig);
      free(reord);
      rc = fc_getSubsetMembersAsArray(copyshape.elems[neworder[i]],
				      &numorig,&orig);
      fail_unless(rc == FC_SUCCESS, "cant get subset members");
      rc = fc_getSubsetMembersAsArray(shapes[0].elems[i],
				      &numreord,&reord);
      fail_unless(rc == FC_SUCCESS, "cant get subset members");
      fail_unless(numorig == numreord, "bad subset reord");
      for (j = 0; j < numorig; j++){
	fail_unless(orig[j] == reord[j], "bad subset reord");
      }
      free(orig);
      free(reord);
      for (j = 0; j < copyshape.numSides; j++){
	fail_unless(copyshape.adjmatrix[neworder[i]][neworder[j]] ==
		    shapes[0].adjmatrix[i][j], "bad adj matrix reorder");
      }
    }
		    
    fc_freeShape(&copyshape);
    free(neworder);

    //----- thin shape order---------
    //for largeandnonadjside order:
    //ii = 0 - legit case
    //       -- failure cases:
    //       - then internally changed so two largest sides are adj      
    //       - then internally changed so there is only one largest side,
    //         and the next two nonadj sides are of equal size
    //       - no non adj sides 
    //       - only 1 side 
    //ii = 1 - the 4 largest sides are of equal area, should fail
    //
    //for largeandopposingside order:
    //ii = 0 - legitcase
    //       - then change so the opposing side is adj, should then find
    //         next pair (this checks both adj and opposing angle)
    //       - no nonadj sides
    //       - no opposing sides
    //
    if (0){
      //get the areas of the current incarnation of the shape
      rc = fc_getShapeSidesAreas(&shapes[0], &areas);
      fail_unless(rc == FC_SUCCESS, "abort: failed to get areas");

      fc_printShape(&shapes[0],"",0);
      printf("\n before call order\n");
      for (i = 0; i <shapes[0].numSides; i++){
	printf("side %d area %10g\n",i,areas[i]);
      }
      free(areas);
    }

    if (ii){
      int *order;

      rc = fc_getShapeSidesAreas(&shapes[0], &areas);
      fail_unless(rc == FC_SUCCESS, "abort: failed to get areas");

      rc = fc_createDescendingAreaOrder(&shapes[0], &order);
      fail_unless(rc == FC_SUCCESS, "abort: failed to get order");

      rc = fc_createLargeAndNonAdjacentSidesOrder(&shapes[0],&neworder);
      fail_unless(rc == FC_SUCCESS,
		  "should work for legit case that doesnt satisfy the large non-adjside criteria");
      fail_unless(neworder == NULL, "should return NULL order");

      //will put something in here like this if put in uniqueness...
      /*
      rc = fc_createLargeAndOpposingSidesOrder(&shapes[0],10,&neworder);
      fail_unless(rc == FC_SUCCESS, "cant get opposing side");
      for (i = 0; i < shapes[0].numSides; i++){
	printf("%d ",neworder[i]);
      }
      fail_unless(((neworder[0] == order[0] && neworder[1] == order[1]) ||
		   (neworder[0] == order[1] && neworder[1] == order[0])),
		  "wrong answer for large and opposing sides");
      free(neworder);
      */

      free(areas);
      free(order);
    }else{
      int *order;
      int a, b, c, d;

      rc = fc_getShapeSidesAreas(&shapes[0], &areas);
      fail_unless(rc == FC_SUCCESS, "abort: failed to get areas");

      rc = fc_createDescendingAreaOrder(&shapes[0], &order);
      fail_unless(rc == FC_SUCCESS, "abort: failed to get order");

      rc = fc_createLargeAndNonAdjacentSidesOrder(&shapes[0],&neworder);
      fail_unless(rc == FC_SUCCESS, "should be thin shape");
      fail_unless(neworder != NULL, "should return an order");
      
      //known case
      fail_unless(FC_DBL_EQUIV(areas[neworder[0]],areas[neworder[1]]),
		  "for this case the two sides of interest should be the same area");
      fail_unless(shapes[0].adjmatrix[neworder[0]][neworder[1]] == 0,
		  "two sides should be non-adjacent");
      fail_unless(areas[neworder[0]] > areas[neworder[2]],
		  "0th side is not the largest side");
      //check that all sides are represented;
      for (i = 0; i < shapes[0].numSides; i++){  
	int found = 0;
	for (j = 0; j < shapes[0].numSides; j++){
	  if (neworder[j] == i){
	    found = 1;
	    break;
	  }
	}
	fail_unless(found,"missing index in neworder");
      }
      free(neworder);


      rc = fc_createLargeAndOpposingSidesOrder(&shapes[0],10,&neworder);
      fail_unless(rc == FC_SUCCESS, "cant get opposing side");
      fail_unless(((neworder[0] == order[0] && neworder[1] == order[1]) ||
		   (neworder[0] == order[1] && neworder[1] == order[0])),
		  "wrong answer for large and opposing sides");
      free(neworder);


      // now mutate the shape so that the two largest sides
      // are adj, and make sure there is another non adj side
      a = order[0];
      b = order[1];
      c = order[2];
      shapes[0].adjmatrix[a][b] = 1;
      shapes[0].adjmatrix[b][a] = 1;
      shapes[0].adjmatrix[a][c] = 0;
      shapes[0].adjmatrix[c][a] = 0;

      rc = fc_createLargeAndNonAdjacentSidesOrder(&shapes[0],&neworder);
      fail_unless(rc == FC_SUCCESS,
		  "should work for legit case that doesnt satisfy the large non-adjside criteria");
      fail_unless(neworder == NULL, "should return NULL order");


      rc = fc_createLargeAndOpposingSidesOrder(&shapes[0],10,&neworder);
      fail_unless(rc == FC_SUCCESS, "cant get opposing side");
      fail_unless(((neworder[0] == order[2] && neworder[1] == order[3]) ||
		   (neworder[0] == order[3] && neworder[1] == order[2])),
		  "wrong answer for large and opposing sides");
      free(neworder);

      //and fix the shape back
      shapes[0].adjmatrix[a][b] = 0;
      shapes[0].adjmatrix[b][a] = 0;
      shapes[0].adjmatrix[a][c] = 1;
      shapes[0].adjmatrix[c][a] = 1;

      // now change so no non adj sides
      for (i = 0; i < shapes[0].numSides; i++){
	for (j = 0; j < shapes[0].numSides; j++){
	  shapes[0].adjmatrix[i][j] = 1;
	}
      }

      rc = fc_createLargeAndOpposingSidesOrder(&shapes[0],10,&neworder);
      fail_unless(rc == FC_SUCCESS, "cant get opposing side");
      fail_unless(neworder == NULL," should return NULL for no largeandopposing sides");

      //and fix the shape back
      for (i = 0; i < shapes[0].numSides; i++){
	for (j = 0; j < shapes[0].numSides; j++){
	  shapes[0].adjmatrix[i][j] = 1;
	}
	shapes[0].adjmatrix[i][i] = 0;
      }
      shapes[0].adjmatrix[order[0]][order[1]] = 0;
      shapes[0].adjmatrix[order[1]][order[0]] = 0;
      shapes[0].adjmatrix[order[2]][order[3]] = 0;
      shapes[0].adjmatrix[order[3]][order[2]] = 0;
      shapes[0].adjmatrix[order[4]][order[5]] = 0;
      shapes[0].adjmatrix[order[5]][order[4]] = 0;


      // now change so no opposing sides
      // use just the first 3 sides, if any of them are the
      // same size then just make them adj
      shapes[0].numSides = 3;
      for (i = 0; i < 3; i++){
	for (j = 0; j < 3; j++){
	  if (FC_FLT_EQUIV(areas[i],areas[j])){
	    shapes[0].adjmatrix[i][j] = 1;
	  }
	}
      }

      rc = fc_createLargeAndOpposingSidesOrder(&shapes[0],10,&neworder);
      fail_unless(rc == FC_SUCCESS, "cant get opposing side");
      fail_unless(neworder == NULL," should return NULL for no largeandopposing sides");

      //and fix the shape back
      shapes[0].numSides = 6;
      for (i = 0; i < shapes[0].numSides; i++){
	for (j = 0; j < shapes[0].numSides; j++){
	  shapes[0].adjmatrix[i][j] = 1;
	}
	shapes[0].adjmatrix[i][i] = 0;
      }
      shapes[0].adjmatrix[order[0]][order[1]] = 0;
      shapes[0].adjmatrix[order[1]][order[0]] = 0;
      shapes[0].adjmatrix[order[2]][order[3]] = 0;
      shapes[0].adjmatrix[order[3]][order[2]] = 0;
      shapes[0].adjmatrix[order[4]][order[5]] = 0;
      shapes[0].adjmatrix[order[5]][order[4]] = 0;


      //mutations with subset
      {
	FC_Subset tempSubset;

	//now mutate so that there is only one largest side,
	//and the next two nonadj sides are of equal size
	a = order[0];
	b = order[1];
	c = order[2];
	d = order[3];

	tempSubset = shapes[0].faces[b];
	shapes[0].faces[b] = shapes[0].faces[order[shapes[0].numSides-1]];
	shapes[0].adjmatrix[a][c] = 0;
	shapes[0].adjmatrix[c][a] = 0;
	shapes[0].adjmatrix[a][d] = 0;
	shapes[0].adjmatrix[d][a] = 0;

	rc = fc_createLargeAndNonAdjacentSidesOrder(&shapes[0],&neworder);
	fail_unless(rc == FC_SUCCESS,
		    "should work for legit case that doesnt satisfy the large non-adjside criteria");
	fail_unless(neworder == NULL, "should return NULL order");

	//dont fix the shape back, but make the
	//everything adj to the first side
	for (i = 1; i < shapes[0].numSides; i++){
	  shapes[0].adjmatrix[order[0]][order[i]] = 1;
	  shapes[0].adjmatrix[order[i]][order[0]] = 1;
	}


	rc = fc_createLargeAndNonAdjacentSidesOrder(&shapes[0],&neworder);
	fail_unless(rc == FC_SUCCESS,
		    "should work for legit case that doesnt satisfy the large non-adjside criteria");
	fail_unless(neworder == NULL, "should return NULL order");


	//now fix the shape back
	for (i = 0; i < shapes[0].numSides; i++){
	  for (j = 0; j < shapes[0].numSides; j++){
	    shapes[0].adjmatrix[i][j] = 1;
	  }
	}
	for (i = 0; i < shapes[0].numSides; i++){
	  shapes[0].adjmatrix[i][i] = 0;
	}
	shapes[0].adjmatrix[order[0]][order[1]] = 0;
	shapes[0].adjmatrix[order[1]][order[0]] = 0;
	shapes[0].adjmatrix[order[2]][order[3]] = 0;
	shapes[0].adjmatrix[order[3]][order[2]] = 0;
	shapes[0].adjmatrix[order[4]][order[5]] = 0;
	shapes[0].adjmatrix[order[5]][order[4]] = 0;

	shapes[0].faces[b] = tempSubset;


	//only 1 side
	a = shapes[0].numSides;
	shapes[0].numSides = 1;
	  
	rc = fc_createLargeAndNonAdjacentSidesOrder(&shapes[0],&neworder);
	fail_unless(rc == FC_SUCCESS,
		    "should work for legit case that doesnt satisfy the large non-adjside criteria");
	fail_unless(neworder == NULL, "should return NULL order");

	//fix shape back
	shapes[0].numSides = a;
      } //end mutating with subset

      free(areas);
      free(order);
    }

    //clean up for next iteration
    if (ii ==0){
      for (j = 0; j < numShapes; j++){
	fc_freeShape(&shapes[j]);
      }
      free(shapes);
      fc_deleteMesh(mesh);

      //set two sides equal. also makes sure that 1st one isnt always the largest one
      uppers[0] = uppers[2];
    }
  }
  //recall havent cleaned up meshes or shapes on last mesh yet...

  //empty side test
  rc = fc_createSubset(mesh,"emptyfaces", FC_AT_FACE, &emptyf);
  fail_unless(rc == FC_SUCCESS, "cant createsubset");

  test = shapes[0].faces[1];
  shapes[0].faces[1] = emptyf;

  rc = fc_getShapeSidesAreas(&shapes[0],&areas);
  fail_unless(rc == FC_SUCCESS, "should work for empty side");
  fail_unless(FC_DBL_EQUIV(areas[1],0.0), "wrong return val for empty side area");
  free(areas);
  fc_deleteSubset(emptyf);

  //bad args for side areas
  shapes[0].faces[1] = FC_NULL_SUBSET;
  rc = fc_getShapeSidesAreas(&shapes[0],&areas);
  fail_unless(rc != FC_SUCCESS, "should fail for NULL side");
  fail_unless(areas == NULL, "should return NULL when fail");
  shapes[0].faces[1] = test; //fix the shape

  rc = fc_getShapeSidesAreas(NULL,&areas);
  fail_unless(rc != FC_SUCCESS, "should fail for NULL shape");
  fail_unless(areas == NULL, "should return NULL when fail");

  rc = fc_getShapeSidesAreas(&shapes[0],NULL);
  fail_unless(rc != FC_SUCCESS, "should fail for NULL return arg");


  //bad args for reorder
  //lets get that copy back for reorder
  rc = fc_copyShape(&shapes[0],&copyshape);
  fail_unless(rc == FC_SUCCESS, "should be able to copy");
  rc = fc_reorderShape(NULL,&copyshape);
  fail_unless(rc != FC_SUCCESS, "should fail for NULL input array ");
  //if reorder fails, shape should be unchanged.
  //copying over check from above to test this
  for (i = 0; i <shapes[0].numSides; i++){
    int *orig;
    int *reord;
    int numorig, numreord;
    rc = fc_getSubsetMembersAsArray(copyshape.faces[i],&numorig,&orig);
    fail_unless(rc == FC_SUCCESS, "cant get subset members");
    rc = fc_getSubsetMembersAsArray(shapes[0].faces[i],
				    &numreord,&reord);
    fail_unless(rc == FC_SUCCESS, "cant get subset members");
    fail_unless(numorig == numreord, "bad subset reord");
    for (j = 0; j < numorig; j++){
      fail_unless(orig[j] == reord[j], "bad subset reord");
    }
    free(orig);
    free(reord);
    rc = fc_getSubsetMembersAsArray(copyshape.elems[i],&numorig,&orig);
    fail_unless(rc == FC_SUCCESS, "cant get subset members");
    rc = fc_getSubsetMembersAsArray(shapes[0].elems[i],
				    &numreord,&reord);
    fail_unless(rc == FC_SUCCESS, "cant get subset members");
    fail_unless(numorig == numreord, "bad subset reord");
    for (j = 0; j < numorig; j++){
      fail_unless(orig[j] == reord[j], "bad subset reord");
    }
    fail_unless(orig[0] == reord[0], "bad subset reord");
    free(orig);
    free(reord);
    for (j = 0; j < copyshape.numSides; j++){
      fail_unless(shapes[0].adjmatrix[i][j] ==
		  copyshape.adjmatrix[i][j], "bad adj matrix reorder");
    }
  }

  //now get an order array and make it bad
  rc = fc_createDescendingAreaOrder(&shapes[0],&neworder);
  fail_unless(rc == FC_SUCCESS, "abort: failed to get descending area order");
  neworder[0] = -1;
  rc = fc_reorderShape(neworder,&copyshape);
  fail_unless(rc != FC_SUCCESS, "should fail for badorder array ");
  neworder[0] = neworder[1];
  rc = fc_reorderShape(neworder,&copyshape);
  fail_unless(rc != FC_SUCCESS, "should fail for badorder array ");
  free(neworder);
  //and one more check of the copyshape
  for (i = 0; i <shapes[0].numSides; i++){
    int *orig;
    int *reord;
    int numorig, numreord;
    rc = fc_getSubsetMembersAsArray(copyshape.faces[i],&numorig,&orig);
    fail_unless(rc == FC_SUCCESS, "cant get subset members");
    rc = fc_getSubsetMembersAsArray(shapes[0].faces[i],
				    &numreord,&reord);
    fail_unless(rc == FC_SUCCESS, "cant get subset members");
    fail_unless(numorig == numreord, "bad subset reord");
    for (j = 0; j < numorig; j++){
      fail_unless(orig[j] == reord[j], "bad subset reord");
    }
    free(orig);
    free(reord);
    rc = fc_getSubsetMembersAsArray(copyshape.elems[i],&numorig,&orig);
    fail_unless(rc == FC_SUCCESS, "cant get subset members");
    rc = fc_getSubsetMembersAsArray(shapes[0].elems[i],
				    &numreord,&reord);
    fail_unless(rc == FC_SUCCESS, "cant get subset members");
    fail_unless(numorig == numreord, "bad subset reord");
    for (j = 0; j < numorig; j++){
      fail_unless(orig[j] == reord[j], "bad subset reord");
    }
    free(orig);
    free(reord);
    for (j = 0; j < copyshape.numSides; j++){
      fail_unless(shapes[0].adjmatrix[i][j] ==
		  copyshape.adjmatrix[i][j], "bad adj matrix reorder");
    }
  }
  fc_freeShape(&copyshape);


  rc = fc_createDescendingAreaOrder(&shapes[0],&neworder);
  fail_unless(rc == FC_SUCCESS, "abort: failed to get descending area order");
  rc = fc_reorderShape(neworder,NULL);
  fail_unless(rc != FC_SUCCESS, "should fail for NULL shape ");
  free(neworder);

  //--- bad args for descending --- //
  rc = fc_createDescendingAreaOrder(NULL,&neworder);
  fail_unless(rc != FC_SUCCESS, "should fail NULL shape");
  fail_unless(neworder == NULL, "should return NULL when fail");
  rc = fc_createDescendingAreaOrder(&shapes[0],NULL);
  fail_unless(rc != FC_SUCCESS, "should fail NULL returnarg");

  //mess with the shape - making bad sides
  holdsubset = shapes[0].faces[0];
  shapes[0].faces[0] = FC_NULL_SUBSET;
  rc = fc_createDescendingAreaOrder(&shapes[0],&neworder);
  fail_unless(rc != FC_SUCCESS, "should fail for NULL subset");
  shapes[0].faces[0] = badSubset;
  rc = fc_createDescendingAreaOrder(&shapes[0],&neworder);
  fail_unless(rc != FC_SUCCESS, "should fail for bad subset");
  shapes[0].faces[0] = holdsubset;

  //--- bad args for large nonadj side ---//
  rc = fc_createLargeAndNonAdjacentSidesOrder(NULL,&neworder);
  fail_unless(rc != FC_SUCCESS, "should fail NULL shape");
  fail_unless(neworder == NULL, "should return NULL when fail");
  rc = fc_createLargeAndNonAdjacentSidesOrder(&shapes[0],NULL);
  fail_unless(rc != FC_SUCCESS, "should fail NULL return arg");

  //mess with the shape - no non adj sides
  rc = fc_copyShape(&shapes[0],&copyshape);
  fail_unless(rc == FC_SUCCESS, "should be able to copy");
  for (i = 0; i < copyshape.numSides; i++){
    for (j = i+1; j < copyshape.numSides; j++){
      copyshape.adjmatrix[i][j] =1;
      copyshape.adjmatrix[j][i] =1;
    }
  }
  fc_freeShape(&copyshape);

  //mess with the shape - making bad sides
  holdsubset = shapes[0].faces[0];
  shapes[0].faces[0] = FC_NULL_SUBSET;
  rc = fc_createLargeAndNonAdjacentSidesOrder(&shapes[0],&neworder);
  fail_unless(rc != FC_SUCCESS, "should fail for NULL subset");
  shapes[0].faces[0] = badSubset;
  rc = fc_createLargeAndNonAdjacentSidesOrder(&shapes[0],&neworder);
  fail_unless(rc != FC_SUCCESS, "should fail for bad subset");
  shapes[0].faces[0] = holdsubset;

  // --- bad args for large opposing side ---//
  rc = fc_createLargeAndOpposingSidesOrder(NULL,90,&neworder);
  fail_unless(rc != FC_SUCCESS, "should fail NULL shape");
  fail_unless(neworder == NULL, "should return NULL when fail");
  rc = fc_createLargeAndOpposingSidesOrder(&shapes[0],-1,&neworder);
  fail_unless(rc != FC_SUCCESS, "should fail for bad angle");
  fail_unless(neworder == NULL, "should return NULL when fail");
  rc = fc_createLargeAndOpposingSidesOrder(&shapes[0],181,&neworder);
  fail_unless(rc != FC_SUCCESS, "should fail for bad angle");
  fail_unless(neworder == NULL, "should return NULL when fail");
  rc = fc_createLargeAndOpposingSidesOrder(&shapes[0],90,NULL);
  fail_unless(rc != FC_SUCCESS, "should fail NULL return arg");

  //mess with the shape - making bad sides
  holdsubset = shapes[0].faces[0];
  shapes[0].faces[0] = FC_NULL_SUBSET;
  rc = fc_createLargeAndOpposingSidesOrder(&shapes[0],90,&neworder);
  fail_unless(rc != FC_SUCCESS, "should fail for NULL subset");
  shapes[0].faces[0] = badSubset;
  rc = fc_createLargeAndOpposingSidesOrder(&shapes[0],90,&neworder);
  fail_unless(rc != FC_SUCCESS, "should fail for bad subset");
  shapes[0].faces[0] = holdsubset;

  //clean up the rest
  for (j = 0; j < numShapes; j++){
    fc_freeShape(&shapes[j]);
  }
  free(shapes);

  fc_deleteMesh(mesh);
  fc_deleteDataset(dataset);
}
END_TEST



// Tests 1)createSideAdjacencyMatrix 2)calcDistanceMatrix 3) segmentadjmatrix
// 4) createSubsetShapes 5) createMeshShapes 
// in all cases i check that the num of faces and elem is correct per side
// and that the adj matrix has the right disribution. in all but one
// case i check that that the faces are indeed children of the elems per side;
// in the remaining case i will check that sides do line up properly with the 
// adj matrix -- that last is not done yet (check if ive done this....)
// also do i want checking of the private functions since they are key
START_TEST(shapematrix)
{
  FC_ReturnCode rc;
  FC_Dataset dataset;
  FC_Mesh mesh;

  FC_Subset subset, testSubset,emptySubset;
  FC_Subset *facesx;
  FC_Subset badSubset={999,999};
  
  //generating meshes and subsets to look at
  FC_Coords lowers = { 0., 0., 0. };
  FC_Coords uppers = { 1., 1., 1. };

  //generating meshes and subsets to look at
  int elemIDS[2];
  int numElem;

  //comparing the shape adj matrix to the calc
  //adj matrix, also producing sides for the calc adj matrix
  //also for checking shape
  FC_Shape *shapes;
  int numShapes;
  int numids, *ids;
  _FC_INTPAIR *checkfaces, *checkelems;
  int numcheckfaces, numcheckelems;

  //adj matrix and dist matrix
  int **adjmatrix_0, **distmatrix;

  // segmentadjmatrix
  int numSeg;
  int *sidesPerSeg;
  int **segIDS;

  //other
  int i,j,k;


  // setup - create test dataset
  rc = fc_createDataset("temp.xxx", &dataset);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create dataset");

  rc = fc_createSimpleHexMesh(dataset, "junk", 2, 2, 2, lowers, uppers,
				&mesh);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create simple hex mesh");


  //special case - empty subset
  rc = fc_createSubset(mesh,"temp",FC_AT_FACE,&emptySubset);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create subset");
  rc = fc_createSideAdjacencyMatrix(1,&emptySubset,0,&adjmatrix_0);
  fail_unless(rc != FC_SUCCESS, "should fail for empty subset");
  fail_unless(adjmatrix_0 == NULL, "bad return val for empty subset");
  fc_deleteSubset(emptySubset);

  //special case - full mesh. 
  rc = fc_getMeshShapes(mesh,85,1,&numShapes,&shapes);
  fail_unless(rc == FC_SUCCESS, "abort: failed to get mesh Sides");
  fail_unless(numShapes == 1, "wrong number of mesh shapes");
  fail_unless(shapes[0].numSides == 6, "wrong number of mesh sides");
  rc = fc_getMeshSkin(mesh,&testSubset);
  fail_unless(rc == FC_SUCCESS, "cant get mesh skin");
  rc = fc_getSubsetMembersAsArray(testSubset,&numids,&ids);
  fail_unless(rc == FC_SUCCESS, "cant get mesh skin as array");
  numcheckfaces = numids;
  checkfaces = (_FC_INTPAIR*)malloc(numcheckfaces*sizeof(_FC_INTPAIR));
  numcheckelems = 8;
  checkelems = (_FC_INTPAIR*)malloc(numcheckelems*sizeof(_FC_INTPAIR));
  fail_unless(checkfaces != NULL && checkelems != NULL, "memory error");
  for (i = 0; i < numcheckfaces; i++){
    checkfaces[i].ID1 = ids[i];
    checkfaces[i].ID2 = 0;
  }
  free(ids);
  for (i = 0; i < numcheckelems; i++){ //since its a full mesh numbers from 0 to numelem
    checkelems[i].ID1 = i;
    checkelems[i].ID2 = 0;
  }
  for (i = 0; i < shapes[0].numSides; i++){
    int numfaces, *faces;
    int numelems, *elems;
    rc = fc_getSubsetMembersAsArray(shapes[0].faces[i],&numfaces,&faces);
    fail_unless(rc == FC_SUCCESS, "cant get face ids");
    fail_unless(numfaces == 4, "wrong num faces members");
    for (j = 0; j < 4; j++){
      for (k = 0; k < numcheckfaces; k++){
	if (faces[j] == checkfaces[k].ID1){
	  checkfaces[k].ID2++;
	  break;
	}
      }
    }
    rc = fc_getSubsetMembersAsArray(shapes[0].elems[i],&numelems,&elems);
    fail_unless(rc == FC_SUCCESS, "cant get elem ids");
    fail_unless(numelems == 4, "wrong num elem members");
    for (j = 0; j < numelems; j++){
      for (k = 0; k < numcheckelems; k++){
	if (elems[j] == checkelems[k].ID1){
	  checkelems[k].ID2++;
	  break;
	}
      }
    }

    //make sure that every face is a child of an element for this side.
    {
      FC_SortedIntArray facelist;
      int numidsx, *idsx;
      int jj, kk;

      rc = fc_initSortedIntArray(&facelist);
      fail_unless(rc == FC_SUCCESS, "cant init sorted int array");

      for (jj = 0; jj < numelems; jj++){
	rc = fc_getMeshEntityChildren(mesh,FC_AT_ELEMENT,
				      elems[jj],
				      FC_AT_FACE, &numidsx,&idsx);
	fail_unless(rc == FC_SUCCESS, "cant get elem children");
	for (kk = 0; kk < numidsx; kk++){
	  fc_addIntToSortedIntArray(&facelist,idsx[kk]);
	  fail_unless(rc == FC_SUCCESS, "cant add face to SIA");
	}
	free(idsx);
      }
      for (jj = 0; jj < numfaces; jj++){ 
	int ret = fc_isIntInSortedIntArray(&facelist,faces[jj]);
	fail_unless(rc == FC_SUCCESS, "cant lookup int in SIA");
	fail_unless(ret, "face isnt on an elem of this side");
      }
	
      fc_freeSortedIntArray(&facelist);
    }

    free(faces);
    free(elems);
  }

  //each skin face should be in the return set once
  for (i = 0; i < numcheckfaces; i++){
    fail_unless(checkfaces[i].ID2 == 1, "wrong return face");
  }
  //each skin elem should be in the return set 3 times for this case
  for (i = 0; i < numcheckelems; i++){
    fail_unless(checkelems[i].ID2 == 3, "wrong return elem");
  }
  fc_deleteSubset(testSubset);
  free(checkfaces);
  free(checkelems);
  
  rc = fc_createSideAdjacencyMatrix(shapes[0].numSides,shapes[0].faces,
				    0,&adjmatrix_0);
  fail_unless(rc == FC_SUCCESS, "abort: failed to get adj matrix");
  
  rc = fc_calcDistanceMatrix(shapes[0].numSides,adjmatrix_0,&distmatrix);
  fail_unless(rc == FC_SUCCESS, "abort: failed to get distmatrix");


  //checks both adj matricies against each other as well
  for (i = 0; i < shapes[0].numSides; i++){
    int count = 0;
    for (j = 0; j < shapes[0].numSides; j++){
      fail_unless(adjmatrix_0[i][j] == adjmatrix_0[j][i],
		  "matrix should be symmetric");
      fail_unless(distmatrix[i][j] == distmatrix[j][i],
		  "matrix should be symmetric");
      fail_unless(adjmatrix_0[i][j] == shapes[0].adjmatrix[i][j],
		  "matricies should be the same");
      if (i == j){
	fail_unless(adjmatrix_0[i][j] == 0, "matrix should be 0 on diagonal");
	fail_unless(distmatrix[i][j] == 0, "matrix should be 0 on diagonal");
      }else{
	fail_unless(distmatrix[i][j] == (adjmatrix_0[i][j]? 1: 2),
		    "bad val for dist matrix");
      }
      if (adjmatrix_0[i][j] == 1) count++;
    }
    fail_unless(count == 4, "wrong number of adj sides");
  }

  rc = fc_segmentAdjacencyMatrix(shapes[0].numSides,shapes[0].adjmatrix,
				 &numSeg,&sidesPerSeg,&segIDS);
  fail_unless(rc == FC_SUCCESS, "cant segment adj matrix");
  fail_unless(numSeg ==1, "wrong num of seg");
  fail_unless(sidesPerSeg[0] == shapes[0].numSides, "wrong num of seg");
  //all sides should be in this seg
  for (i = 0; i < shapes[0].numSides; i++){
    int found = 0;
    for (j = 0; j < shapes[0].numSides; j++){
      if (segIDS[0][i] == j){
	found = 1;
	break;
      }
    }
    fail_unless(found, "wrong side in segment id list");
  }

  for (i = 0; i < shapes[0].numSides; i++){
    free(adjmatrix_0[i]);
    free(distmatrix[i]);
  }
  free(adjmatrix_0);
  free(distmatrix);

  for (i = 0; i < numSeg; i++){
    free(segIDS[i]);
  }
  free(sidesPerSeg);
  free(segIDS);

  fc_freeShape(&shapes[0]);
  free(shapes);

  //special case - 1 block
  numElem = 1;
  elemIDS[0] = 0;

  rc = fc_createSubset(mesh,"temp",FC_AT_ELEMENT,&subset);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create subset");

  rc = fc_addArrayMembersToSubset(subset,numElem,elemIDS);
  fail_unless(rc == FC_SUCCESS, "abort: failed to add members to subset");

  rc = fc_getSubsetShapes(subset,85,1,&numShapes,&shapes);
  fail_unless(rc == FC_SUCCESS, "abort: failed to get  subset Sides");
  fail_unless(numShapes == 1, "wrong number of mesh shapes");
  fail_unless(shapes[0].numSides == 6, "wrong number of mesh sides");
  rc = fc_getSubsetSkin(subset,&testSubset);
  fail_unless(rc == FC_SUCCESS, "cant get mesh skin");
  rc = fc_getSubsetMembersAsArray(testSubset,&numids,&ids);
  fail_unless(rc == FC_SUCCESS, "cant get mesh skin as array");
  numcheckfaces = numids;
  checkfaces = (_FC_INTPAIR*)malloc(numcheckfaces*sizeof(_FC_INTPAIR));
  fail_unless(checkfaces != NULL, "memory error");
  for (i = 0; i < numcheckfaces; i++){
    checkfaces[i].ID1 = ids[i];
    checkfaces[i].ID2 = 0;
  }
  free(ids);
  for (i = 0; i < shapes[0].numSides; i++){
    int numfaces, *faces;
    int numelems, *elems;
    rc = fc_getSubsetMembersAsArray(shapes[0].faces[i],&numfaces,&faces);
    fail_unless(rc == FC_SUCCESS, "cant get face ids");
    fail_unless(numfaces == 1, "wrong num faces members");
    for (k = 0; k < numcheckfaces; k++){
      if (faces[0] == checkfaces[k].ID1){
	  checkfaces[k].ID2++;
	  break;
      }
    }
    rc = fc_getSubsetMembersAsArray(shapes[0].elems[i],&numelems,&elems);
    fail_unless(rc == FC_SUCCESS, "cant get elem ids");
    fail_unless(numelems == 1, "wrong num elem members");
    fail_unless(elems[0] == elemIDS[0], "wrong elem id");

    //make sure that every face is a child of an element for this side.
    {
      int numidsx, *idsx;
      int jj, kk;

      rc = fc_getMeshEntityChildren(mesh,FC_AT_ELEMENT,elems[0],
				    FC_AT_FACE, &numidsx,&idsx);
      fail_unless(rc == FC_SUCCESS, "cant get elem children");
      for (jj = 0; jj < numfaces; jj++){
	int found = 0;
	for (kk = 0; kk < numidsx; kk++){
	  if (faces[jj] == idsx[kk]){
	    found = 1;
	    break;
	  }
	}
	fail_unless(found, "face isnt on an elem of this side");
      }
      free(idsx);
    }

    free(faces);
    free(elems);
  }  
  //each skin face should be in the return set once
  for (i = 0; i < numcheckfaces; i++){
    fail_unless(checkfaces[i].ID2 == 1, "wrong return face");
  }

  fc_deleteSubset(testSubset);
  free(checkfaces);

  rc = fc_createSideAdjacencyMatrix(shapes[0].numSides,shapes[0].faces,
				    1,&adjmatrix_0);
  fail_unless(rc == FC_SUCCESS, "abort: failed to get adj matrix");


  for (i = 0; i < shapes[0].numSides; i++){
    int count = 0;
    for (j = 0; j < shapes[0].numSides; j++){
      fail_unless(adjmatrix_0[i][j] == adjmatrix_0[j][i],
		  "matrix should be symmetric");
      fail_unless(shapes[0].adjmatrix[i][j] == adjmatrix_0[i][j],
		  "matricies should be the same");
      if (i == j){
	fail_unless(adjmatrix_0[i][j] == 0, "matrix should be 0 on diagonal");
      }
      if (adjmatrix_0[i][j] == 1) count++;
    }
    fail_unless(count == 4, "wrong number of adj sides");
  }


  rc = fc_segmentAdjacencyMatrix(shapes[0].numSides,shapes[0].adjmatrix,
				 &numSeg,&sidesPerSeg,&segIDS);
  fail_unless(rc == FC_SUCCESS, "cant segment adj matrix");
  fail_unless(numSeg ==1 , "wrong num of seg");
  fail_unless(sidesPerSeg[0] == shapes[0].numSides , "wrong num of seg");
  //all sides should be in this seg
  for (i = 0; i < shapes[0].numSides; i++){
    int found = 0;
    for (j = 0; j < shapes[0].numSides; j++){
      if (segIDS[0][i] == j){
	found = 1;
	break;
      }
    }
    fail_unless(found, "wrong side in segment id list");
  }


  for (i = 0; i < numSeg; i++){
    free(segIDS[i]);
  }
  free(sidesPerSeg);
  free(segIDS);

  for (i = 0; i < shapes[0].numSides; i++){
    free(adjmatrix_0[i]);
  }
  free(adjmatrix_0);
  fc_freeShape(&shapes[0]);
  free(shapes);

  fc_deleteSubset(subset);

  //special case - adj 2 blocks
  numElem = 2;
  elemIDS[0] = 0;
  elemIDS[1] = 2;

  rc = fc_createSubset(mesh,"temp",FC_AT_ELEMENT,&subset);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create subset");

  rc = fc_addArrayMembersToSubset(subset,numElem,elemIDS);
  fail_unless(rc == FC_SUCCESS, "abort: failed to add members to subset");

  rc = fc_getSubsetShapes(subset,85,1,&numShapes,&shapes);
  fail_unless(rc == FC_SUCCESS, "abort: failed to get  subset Sides");
  fail_unless(numShapes == 1, "wrong number of shapes");
  fail_unless(shapes[0].numSides == 6, "wrong number of mesh sides");
  rc = fc_getSubsetSkin(subset,&testSubset);
  fail_unless(rc == FC_SUCCESS, "cant get subset skin");
  rc = fc_getSubsetMembersAsArray(testSubset,&numids,&ids);
  fail_unless(rc == FC_SUCCESS, "cant get mesh skin as array");
  numcheckfaces = numids;
  checkfaces = (_FC_INTPAIR*)malloc(numcheckfaces*sizeof(_FC_INTPAIR));
  numcheckelems = numElem;
  checkelems = (_FC_INTPAIR*)malloc(numcheckelems*sizeof(_FC_INTPAIR));
  fail_unless(checkfaces != NULL && checkelems != NULL, "memory error");
  for (i = 0; i < numcheckfaces; i++){
    checkfaces[i].ID1 = ids[i];
    checkfaces[i].ID2 = 0;
  }
  free(ids);
  for (i = 0; i < numcheckelems; i++){
    checkelems[i].ID1 = elemIDS[i];
    checkelems[i].ID2 = 0;
  }
  for (i = 0; i < shapes[0].numSides; i++){
    int numfaces, *faces;
    int numelems, *elems;
    rc = fc_getSubsetMembersAsArray(shapes[0].faces[i],&numfaces,&faces);
    fail_unless(rc == FC_SUCCESS, "cant get face ids");
    rc = fc_getSubsetMembersAsArray(shapes[0].elems[i],&numelems,&elems);
    fail_unless(rc == FC_SUCCESS, "cant get elem ids");
    switch (numelems){
    case 2:
      fail_unless(elems[0] == elemIDS[0], "wrong elem id"); //its an ordered array
      fail_unless(elems[1] == elemIDS[1], "wrong elem id");
      fail_unless(numfaces == 2, "wrong num of faces");
      break;
    case 1:
      fail_unless(elems[0] == elemIDS[0] || elems[0] == elemIDS[1], "wrong elem id"); 
      fail_unless(numfaces == 1, "wrong num of faces");
      break;
    default:
      fail_unless(0, "wrong num of elems");
    }

    for (j = 0; j < numfaces; j++){
      for (k = 0; k < numcheckfaces; k++){
	if (faces[j] == checkfaces[k].ID1){
	  checkfaces[k].ID2++;
	  break;
	}
      }
    }

    for (j = 0; j < numelems; j++){
      for (k = 0; k < numcheckelems; k++){
	if (elems[j] == checkelems[k].ID1){
	  checkelems[k].ID2++;
	  break;
	}
      }
    }

    //make sure that every face is a child of an element for this side.
    {
      FC_SortedIntArray facelist;
      int numidsx, *idsx;
      int jj, kk;

      rc = fc_initSortedIntArray(&facelist);
      fail_unless(rc == FC_SUCCESS, "cant init sorted int array");

      for (jj = 0; jj < numelems; jj++){
	rc = fc_getMeshEntityChildren(mesh,FC_AT_ELEMENT,
				      elems[jj],
				      FC_AT_FACE, &numidsx,&idsx);
	fail_unless(rc == FC_SUCCESS, "cant get elem children");
	for (kk = 0; kk < numidsx; kk++){
	  fc_addIntToSortedIntArray(&facelist,idsx[kk]);
	  fail_unless(rc == FC_SUCCESS, "cant add face to SIA");
	}
	free(idsx);
      }

      for (jj = 0; jj < numfaces; jj++){
	int ret = fc_isIntInSortedIntArray(&facelist,faces[jj]);
	fail_unless(rc == FC_SUCCESS, "cant lookup int in SIA");
	fail_unless(ret, "face isnt on an elem of this side");
      }
	
      fc_freeSortedIntArray(&facelist);
    }

    free(faces);
    free(elems);
  }  
  //each skin face should be in the return set once
  for (i = 0; i < numcheckfaces; i++){
    fail_unless(checkfaces[i].ID2 == 1, "wrong return face");
  }
  //each skin elem should be in the return set 5 times
  for (i = 0; i < numcheckelems; i++){
    fail_unless(checkelems[i].ID2 == 5, "wrong return elem");
  }
  fc_deleteSubset(testSubset);
  free(checkfaces);
  free(checkelems);

  rc = fc_createSideAdjacencyMatrix(shapes[0].numSides,shapes[0].faces,
				    1,&adjmatrix_0);
  fail_unless(rc == FC_SUCCESS, "abort: failed to get adj matrix");


  for (i = 0; i < shapes[0].numSides; i++){
    int count = 0;
    for (j = 0; j < shapes[0].numSides; j++){
      fail_unless(adjmatrix_0[i][j] == adjmatrix_0[j][i],
		  "matrix should be symmetric");
      fail_unless(shapes[0].adjmatrix[i][j] == adjmatrix_0[i][j],
		  "matricies should be the same");
      if (i == j){
	fail_unless(adjmatrix_0[i][j] == 0, "matrix should be 0 on diagonal");
      }
      if (adjmatrix_0[i][j] == 1) count++;
    }
    fail_unless(count == 4, "wrong number of adj sides");
  }


  rc = fc_segmentAdjacencyMatrix(shapes[0].numSides,shapes[0].adjmatrix,
				 &numSeg,&sidesPerSeg,&segIDS);
  fail_unless(rc == FC_SUCCESS, "cant segment adj matrix");
  fail_unless(numSeg ==1 , "wrong num of seg");
  fail_unless(sidesPerSeg[0] == shapes[0].numSides , "wrong num of seg");
  //all sides should be in this seg
  for (i = 0; i < shapes[0].numSides; i++){
    int found = 0;
    for (j = 0; j < shapes[0].numSides; j++){
      if (segIDS[0][i] == j){
	found = 1;
	break;
      }
    }
    fail_unless(found, "wrong side in segment id list");
  }


  for (i = 0; i < numSeg; i++){
    free(segIDS[i]);
  }
  free(sidesPerSeg);
  free(segIDS);

  for (i = 0; i < shapes[0].numSides; i++){
    free(adjmatrix_0[i]);
  }
  free(adjmatrix_0);
  fc_freeShape(&shapes[0]);
  free(shapes);

  fc_deleteSubset(subset);

  //this is the big test for adj - takes into consideration diffences
  // in the segdim flag

  //this configuration also checks single element blocks
  //known vals
  elemIDS[0] = 0;
  elemIDS[1] = 3;

  rc = fc_createSubset(mesh,"temp",FC_AT_ELEMENT,&subset);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create subset");

  rc = fc_addArrayMembersToSubset(subset,numElem,elemIDS);
  fail_unless(rc == FC_SUCCESS, "abort: failed to add members to subset");

  for (k = 0; k < 2; k++){
    //k = segdim - determines whether those sides in the
    //same plane are considered as the same side or not. in
    //either case they are considered the same shape.
    int num4 = 0;
    int num6 = 0;
    int num8 = 0;

    int *neighbors;

    rc = fc_getSubsetShapes(subset,85,k,&numShapes,&shapes);
    fail_unless(rc == FC_SUCCESS, "abort: failed to get subset shapes");
    fail_unless(numShapes == 1, "bad num of shapes");
    //more checks for this below...

    rc = fc_createSideAdjacencyMatrix(shapes[0].numSides,shapes[0].faces,
				      k,&adjmatrix_0);
    fail_unless(rc == FC_SUCCESS, "abort: failed to get adj matrix");

    neighbors = (int*)malloc(shapes[0].numSides*sizeof(int));
    fail_unless(neighbors != NULL, "memory error");
    for (i = 0; i < shapes[0].numSides; i++){
      neighbors[i] = 0;
    }
    for (i = 0; i < shapes[0].numSides; i++){
      for (j = 0; j < shapes[0].numSides; j++){
	fail_unless(adjmatrix_0[i][j] == adjmatrix_0[j][i],
		    "adj matrix should be symm"); 
	fail_unless(adjmatrix_0[i][j] == shapes[0].adjmatrix[i][j],
		    "adj matrix should be the same"); 
	if (i == j){
	  fail_unless(adjmatrix_0[i][j] == 0,
		      "adj matrix should be 0 on diagonal"); 
	}
	if (shapes[0].adjmatrix[i][j] == 1){
	  neighbors[i]++;
	}
      }
    }

    if (k == 0){
      //for this case check that adjacency matrix
      //sides should line up as they should.
      //(checks more generally the number of
      //      adj but not their identities below
      // this may
      //not be entirely rigourous, but i dont want
      //to rely on knowing the ids of the faces
      for (i = 0; i < shapes[0].numSides; i++){
	switch (neighbors[i]){
	case 8:
	  for (j = 0; j < shapes[0].numSides; j++){
	    int numelemsi,*elemsi;
	    rc= fc_getSubsetMembersAsArray(shapes[0].elems[i],
					   &numelemsi,&elemsi);
	    fail_unless(rc == FC_SUCCESS, "cant get subset members");
	    fail_unless(numelemsi == 2, "mismatch of neighbors and num elements");
	    free(elemsi);
	    
	    if (i == j || neighbors[j] == 8){
	      fail_unless(shapes[0].adjmatrix[i][j] == 0, 
			  "wrong side adjacency");
	    }else{
	      fail_unless(shapes[0].adjmatrix[i][j] == 1, 
			  "wrong side adjacency");
	    }
	  }
	  break;
	case 6:
	  for (j = 0; j < shapes[0].numSides; j++){
	    //symmetry checks the sides with 8
	    if (neighbors[j] == 6 && i != j){
	      fail_unless(shapes[0].adjmatrix[i][j] == 1, 
			  "wrong side adjacency");
	    }else{
	      if (neighbors[j] == 4){
		int numelemsi,*elemsi;
		int numelemsj,*elemsj;
		rc= fc_getSubsetMembersAsArray(shapes[0].elems[i],
					       &numelemsi,&elemsi);
		fail_unless(rc == FC_SUCCESS, "cant get subset members");
		fail_unless(numelemsi == 1, "wrong num elems for this side");
		rc= fc_getSubsetMembersAsArray(shapes[0].elems[j],
					       &numelemsj,&elemsj);
		fail_unless(rc == FC_SUCCESS, "cant get subset members");
		//its only adj to one side with 4, and each
		//side with 4 is adj to only one side with 6.
		//this isnt entirely rigourous, because it
		//doesnt uniquely id the side with 4.
		if (elemsi[0] != elemsj[0]){
		  //not adj to any 4 neighbor ones on the other element
		  fail_unless(shapes[0].adjmatrix[i][j] == 0, 
			      "wrong side adjacency");
		}
		//the remainins side must be a side with 4 becuase there
		//are no other sides
		free(elemsi);
		free(elemsj);
	      }
	    }
	    break;
	  case 4:
	    for (j = 0; j < shapes[0].numSides; j++){
	      //symmetry checks the sides with 6 and 8
	      if (neighbors[j] == 4 && i !=j ){
		int numelemsi,*elemsi;
		int numelemsj,*elemsj;
		rc= fc_getSubsetMembersAsArray(shapes[0].elems[i],
					       &numelemsi,&elemsi);
		fail_unless(rc == FC_SUCCESS, "cant get subset members");
		fail_unless(numelemsi == 1, "wrong num elems for this side");
		rc= fc_getSubsetMembersAsArray(shapes[0].elems[j],
					       &numelemsj,&elemsj);
		fail_unless(rc == FC_SUCCESS, "cant get subset members");
		if (elemsi[0] == elemsj[0]){
		  fail_unless(shapes[0].adjmatrix[i][j] == 1, 
			      "wrong side adjacency");
		}else{
		  fail_unless(shapes[0].adjmatrix[i][j] == 0, 
			      "wrong side adjacency");
		}
		free(elemsi);
		free(elemsj);
	      }
	    }
	  }
	  break;
	default:
	  fail_unless(0,"wrong num neighbors");
	}
      }
    }

    for (i = 0; i < shapes[0].numSides; i++){
      switch(neighbors[i]){
      case 4:
	num4++;
	break;
      case 6:
	num6++;
	break;
      case 8:
	num8++;
	break;
      }
    }

    switch(k){
    case 0:
      fail_unless(shapes[0].numSides == 10, "bad num of subset Sides");
      fail_unless(num4 == 4, " bad adj matrix");
      fail_unless(num6 == 4, " bad adj matrix");
      fail_unless(num8 == 2, " bad adj matrix");

      for (i = 0; i < shapes[0].numSides; i++){
	int numelems,*elems;
	int numfaces, *faces;

	rc= fc_getSubsetMembersAsArray(shapes[0].elems[i],&numelems,&elems);
	fail_unless(rc == FC_SUCCESS, "cant get subset members");
	rc= fc_getSubsetMembersAsArray(shapes[0].faces[i],&numfaces,&faces);
	fail_unless(rc == FC_SUCCESS, "cant get subset members");
	switch (numelems){
	case 2:
	  fail_unless(elems[0] == elemIDS[0], "wrong elem id"); //its an ordered array
	  fail_unless(elems[1] == elemIDS[1], "wrong elem id");
	  fail_unless(numfaces == 2, "wrong num of faces");
	  break;
	case 1:
	  fail_unless(elems[0] == elemIDS[0] || elems[0] == elemIDS[1],
		      "wrong elem id"); 
	  fail_unless(numfaces == 1, "wrong num of faces");
	  break;
	default:
	  fail_unless(0, "wrong num elems");
	}
	free(elems);
	free(faces);
      } 
      break;
    case 1:
      fail_unless(shapes[0].numSides == 12, "bad num of subset Sides");
      fail_unless(num4 == 8, " bad adj matrix");
      fail_unless(num6 == 4, " bad adj matrix");
      for (i = 0; i < shapes[0].numSides; i++){
	int numelems,*elems;
	int numfaces, *faces;

	rc= fc_getSubsetMembersAsArray(shapes[0].elems[i],&numelems,&elems);
	fail_unless(rc == FC_SUCCESS, "cant get subset members");
	rc= fc_getSubsetMembersAsArray(shapes[0].faces[i],&numfaces,&faces);
	fail_unless(rc == FC_SUCCESS, "cant get subset members");
	fail_unless(numelems == 1 , "wrong num of subset members");
	fail_unless(elems[0] == elemIDS[0] || elems[0] == elemIDS[1],
		      "wrong elem id"); 
	fail_unless(numfaces == 1, "wrong num of faces");
	free(elems);
	free(faces);
      } 
      break;
    }

    free(neighbors);

    rc = fc_segmentAdjacencyMatrix(shapes[0].numSides,shapes[0].adjmatrix,
				   &numSeg,&sidesPerSeg,&segIDS);
    fail_unless(rc == FC_SUCCESS, "cant segment adj matrix");
    fail_unless(numSeg ==1 ,"wrong num of seg");
    fail_unless(sidesPerSeg[0] == shapes[0].numSides , "wrong num of seg");
    //all sides should be in this seg
    for (i = 0; i < shapes[0].numSides; i++){
      int found = 0;
      for (j = 0; j < shapes[0].numSides; j++){
	if (segIDS[0][i] == j){
	  found = 1;
	  break;
	}
      }
      fail_unless(found, "wrong side in segment id list");
    }

    for (i = 0; i < numSeg; i++){
      free(segIDS[i]);
    }
    free(sidesPerSeg);
    free(segIDS);

    for (i = 0; i < shapes[0].numSides; i++){
      free(adjmatrix_0[i]);
    }
    free(adjmatrix_0);

    fc_freeShape(&shapes[0]);
    free(shapes);
  }
  fc_deleteSubset(subset);


  // this is for infinite dist things and segment, primarily
  numElem = 2;
  elemIDS[0] = 0;
  elemIDS[1] = 7;

  rc = fc_createSubset(mesh,"temp",FC_AT_ELEMENT,&subset);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create subset");

  rc = fc_addArrayMembersToSubset(subset,numElem,elemIDS);
  fail_unless(rc == FC_SUCCESS, "abort: failed to add members to subset");

  rc = fc_getSubsetShapes(subset,85,1,&numShapes,&shapes);
  fail_unless(rc == FC_SUCCESS, "abort: failed to get subset shapes");
  fail_unless(numShapes == 2, "bad num of shapes");
  for (j = 0; j < numShapes; j++){
    for (i = 0; i < shapes[j].numSides; i++){
      int numelems,*elems;
      int numfaces, *faces;
      int numCheckSubsets;
      char* temp_name;
      FC_Subset *checkSubsets;

      temp_name = malloc(50*sizeof(char));
      fail_unless(temp_name != NULL, "memory error");

      snprintf(temp_name, 50,
	       "temp_Shape%d_FaceSurface%d",j,i);
      rc = fc_getSubsetByName(mesh,temp_name, &numCheckSubsets,&checkSubsets);
      fail_unless(rc == FC_SUCCESS, "cant get subset by name");
      fail_unless(numCheckSubsets == 1, "wrong number of matching subsets");
      fail_unless(FC_HANDLE_EQUIV(checkSubsets[0],shapes[j].faces[i]),
		  "subset name mismatch");
      free(checkSubsets);
      snprintf(temp_name, 50,
	       "temp_Shape%d_ElemSurface%d",j,i);
      rc = fc_getSubsetByName(mesh,temp_name,&numCheckSubsets,&checkSubsets);
      fail_unless(rc == FC_SUCCESS, "cant get subset by name");
      fail_unless(numCheckSubsets == 1, "wrong number of matching subsets");
      fail_unless(FC_HANDLE_EQUIV(checkSubsets[0],shapes[j].elems[i]),
		  "subset name mismatch");
      free(checkSubsets);
      free(temp_name);
      rc= fc_getSubsetMembersAsArray(shapes[j].elems[i],&numelems,&elems);
      fail_unless(rc == FC_SUCCESS, "cant get subset members");
      rc= fc_getSubsetMembersAsArray(shapes[j].faces[i],&numfaces,&faces);
      fail_unless(rc == FC_SUCCESS, "cant get subset members");
      fail_unless(numelems == 1 , "wrong num of subset members");
      fail_unless(elems[0] == elemIDS[0] || elems[0] == elemIDS[1],
		  "wrong elem id"); 
      if (!j){   // using k to make sure the two shapes have different elem
	k = elems[0];
      }else{
	fail_unless(elems[0] != k, "repeat of elem in other shape");
      }
      fail_unless(numfaces == 1, "wrong num of faces");

      //make sure that every face is a child of an element for this side.
      {
	int numidsx, *idsx;
	int jj, kk;
	
	rc = fc_getMeshEntityChildren(mesh,FC_AT_ELEMENT,elems[0],
				      FC_AT_FACE, &numidsx,&idsx);
	fail_unless(rc == FC_SUCCESS, "cant get elem children");
	for (jj = 0; jj < numfaces; jj++){
	  int found = 0;
	  for (kk = 0; kk < numidsx; kk++){
	    if (faces[jj] == idsx[kk]){
	      found = 1;
	      break;
	    }
	  }
	  fail_unless(found, "face isnt on an elem of this side");
	}
	free(idsx);
      }

      free(elems);
      free(faces);
    } 
  }

  facesx = (FC_Subset*)malloc(12*sizeof(FC_Subset));
  fail_unless(facesx != NULL, "memory error");
  for (i = 0; i < numShapes; i++){
    for (j = 0; j < 6; j++){ //numbering will only work
                            //if they have same num of sides
      rc = fc_getSubsetNumMember(shapes[i].faces[j],&k);
      fail_unless(rc == FC_SUCCESS, "abort: cant get skin num members");
      fail_unless(k == 1, "wrong num members");
      facesx[i*6+j] = shapes[i].faces[j];
    }
  } //this doesnt rigourously check the numbers of the faces, but
  //if they are wrong then the upcoming stuff will be wrong. 

  //lets make 1 big adj matrix for this item
  rc = fc_createSideAdjacencyMatrix(12,facesx,1,&adjmatrix_0);
  fail_unless(rc == FC_SUCCESS, "abort: failed to get adj matrix");

  //everybody has only 4 adjacencies
  for (i = 0; i < 12; i++){
    int count = 0;
    fail_unless(adjmatrix_0[i][i] == 0, "diag elem should be zero");
    for (j = 0; j < 12; j++){
      fail_unless(adjmatrix_0[i][j] == adjmatrix_0[j][i],
		  "adj matrix should be symm"); 
      fail_unless((adjmatrix_0[i][j] == 0 || adjmatrix_0[i][j] == 1),
		  "wrong val for adj matrix"); 
      if (adjmatrix_0[i][j] == 1) count++;
    }
    fail_unless(count == 4, "wrong num of neighbors");
  }
  //that isnt a rigourous test, but if that messed up,
  //then it will mess up the the segmentadjmatrix and the distance matrix

  rc = fc_calcDistanceMatrix(12,adjmatrix_0,&distmatrix);
  fail_unless(rc == FC_SUCCESS, "abort: failed to get distmatrix");
  for (i = 0; i < 12; i++){
    int num1 = 0;
    int num2 = 0;
    int numinf = 0;
    fail_unless(distmatrix[i][i] == 0, "bad val for dist matrix");
    for (j = 0; j < 12; j++){
      switch(distmatrix[i][j]){
      case 0:
	break;
      case 1:
	num1++;
	break;
      case 2:
	num2++;
	break;
      case -1:
	numinf++;
	break;
      default:
	fail_unless(0,"bad val for distance matrix");
      }
    }
    fail_unless(num1 == 4,"bad val for distance matrix");
    fail_unless(num2 == 1,"bad val for distance matrix");
    fail_unless(numinf == 6,"bad val for distance matrix");
  }

  rc = fc_segmentAdjacencyMatrix(12,adjmatrix_0,
				 &numSeg,&sidesPerSeg,&segIDS);
  fail_unless(rc == FC_SUCCESS, "cant segment adj matrix");
  fail_unless(numSeg == 2 ,"wrong num of seg");
  for (i = 0; i < numSeg; i++){
    fail_unless(sidesPerSeg[i] == 6 , "wrong num of seg");
  }
  for (i = 0; i < 12;i++){ 
    int found = 0;
    for (j = 0; j < 2; j++){
      for(k = 0; k < 6; k++){
	if (i == segIDS[j][k]){
	  found =1;
	  break;
	}
      }
    }
    fail_unless(found, "missing id in segIDS");
  }
  //check the id transfer is generally ok
  for (i = 0; i < 2; i++){ 
    for (j = 0; j < 6; j++){
      for (k = 0; k < 6; k++){
	fail_unless(shapes[i].adjmatrix[j][k] ==
		    adjmatrix_0[segIDS[i][j]][segIDS[i][k]],
		    "mismatch in adj matrix after segmenting");
      }
    }
  }
  for (i = 0; i < numSeg; i++){
    free(segIDS[i]);
  }
  free(sidesPerSeg);
  free(segIDS);

  for (j = 0; j < numShapes; j++){
    fc_freeShape(&shapes[j]);
  }
  free(shapes);

  for (i = 0; i < 12; i++){
    free(adjmatrix_0[i]);
    free(distmatrix[i]);
  }
  free(adjmatrix_0);
  free(distmatrix);
  free(facesx);

  fc_deleteSubset(subset);      

  //---------------bad args for matrix things --------------------------/ 
  //---------------bad shape tests are in their own method--------------/ 
  
  //get a shape back - consecutive blocks
  numElem = 2;
  elemIDS[0] = 0;
  elemIDS[1] = 2;

  rc = fc_createSubset(mesh,"temp",FC_AT_ELEMENT,&subset);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create subset");

  rc = fc_addArrayMembersToSubset(subset,numElem,elemIDS);
  fail_unless(rc == FC_SUCCESS, "abort: failed to add members to subset");

  rc = fc_getSubsetShapes(subset,85,1,&numShapes,&shapes);
  fail_unless(rc == FC_SUCCESS, "abort: failed to get  subset Sides");
  fail_unless(numShapes == 1, "wrong number of shapes");
  fail_unless(shapes[0].numSides == 6, "wrong number of mesh sides");
  
  //bad args for side adj matrix
  rc = fc_createSideAdjacencyMatrix(shapes[0].numSides,shapes[0].faces,-1,&adjmatrix_0);
  fail_unless(rc != FC_SUCCESS, "should fail for bad shared dim");
  fail_unless(adjmatrix_0 == NULL, "should return NULL when fail");
  rc = fc_createSideAdjacencyMatrix(0,shapes[0].faces,0,&adjmatrix_0);
  fail_unless(rc != FC_SUCCESS, "should fail for no subsets");
  fail_unless(adjmatrix_0 == NULL, "should return NULL when fail");
  rc = fc_createSideAdjacencyMatrix(shapes[0].numSides,shapes[0].faces,2,&adjmatrix_0);
  fail_unless(rc != FC_SUCCESS, "should fail for bad shared dim");
  fail_unless(adjmatrix_0 == NULL, "should return NULL when fail");
  rc = fc_createSideAdjacencyMatrix(shapes[0].numSides,shapes[0].elems,0,&adjmatrix_0);
  fail_unless(rc != FC_SUCCESS, "should fail for wrong assoc");
  fail_unless(adjmatrix_0 == NULL, "should return NULL when fail");
  rc = fc_createSideAdjacencyMatrix(shapes[0].numSides,NULL,0,&adjmatrix_0);
  fail_unless(rc != FC_SUCCESS, "should fail for NULL subset");
  fail_unless(adjmatrix_0 == NULL, "should return NULL when fail");
  rc = fc_createSideAdjacencyMatrix(shapes[0].numSides,&badSubset,0,&adjmatrix_0);
  fail_unless(rc != FC_SUCCESS, "should fail for bad subset");
  fail_unless(adjmatrix_0 == NULL, "should return NULL when fail");
  rc = fc_createSideAdjacencyMatrix(shapes[0].numSides,shapes[0].faces,0,NULL);
  fail_unless(rc != FC_SUCCESS, "should fail for NULL return arg");

  //bad args for dist 
  rc = fc_calcDistanceMatrix(0,adjmatrix_0,&distmatrix);
  fail_unless(rc != FC_SUCCESS, "should fail for 0 subsets");
  fail_unless(distmatrix == NULL, "should return NULL when fail");
  rc = fc_calcDistanceMatrix(2,NULL,&distmatrix);
  fail_unless(rc != FC_SUCCESS, "should fail for NULL adj matrix");
  fail_unless(distmatrix == NULL, "should return NULL when fail");
  rc = fc_calcDistanceMatrix(2,adjmatrix_0,NULL);
  fail_unless(rc != FC_SUCCESS, "should fail for NULL return arg");

  //bad args for segment
  rc = fc_segmentAdjacencyMatrix(0,shapes[0].adjmatrix,
				 &numSeg,&sidesPerSeg,&segIDS);
  fail_unless(rc != FC_SUCCESS, "should fail for 0 subsets");
  fail_unless(numSeg == -1, "should return -1 when fail");
  fail_unless(sidesPerSeg == NULL, "should return NULL when fail");
  fail_unless(segIDS == NULL, "should return NULL when fail");

  rc = fc_segmentAdjacencyMatrix(shapes[0].numSides,NULL,
				 &numSeg,&sidesPerSeg,&segIDS);
  fail_unless(rc != FC_SUCCESS, "should fail for NULL adjmatrix");
  fail_unless(numSeg == -1, "should return -1 when fail");
  fail_unless(sidesPerSeg == NULL, "should return NULL when fail");
  fail_unless(segIDS == NULL, "should return NULL when fail");

  rc = fc_segmentAdjacencyMatrix(shapes[0].numSides,shapes[0].adjmatrix,
				 NULL,&sidesPerSeg,&segIDS);
  fail_unless(rc != FC_SUCCESS, "should fail for NULL arg");
  fail_unless(sidesPerSeg == NULL, "should return NULL when fail");
  fail_unless(segIDS == NULL, "should return NULL when fail");

  rc = fc_segmentAdjacencyMatrix(shapes[0].numSides,shapes[0].adjmatrix,
				 &numSeg,NULL,&segIDS);
  fail_unless(rc != FC_SUCCESS, "should fail for NULL arg");
  fail_unless(numSeg == -1, "should return -1 when fail");
  fail_unless(segIDS == NULL, "should return NULL when fail");

  rc = fc_segmentAdjacencyMatrix(shapes[0].numSides,shapes[0].adjmatrix,
				 &numSeg,&sidesPerSeg,NULL);
  fail_unless(rc != FC_SUCCESS, "should fail for NULL arg");
  fail_unless(numSeg == -1, "should return -1 when fail");
  fail_unless(sidesPerSeg == NULL, "should return NULL when fail");

  //now send in a bad adj matrix  - note only
  //segmetn checks for bad adj matrix, the calc in dsitance matrix
  //is supposed to be general and handle non symm and weighted
  //matricies, thoguh i dont guarentee that

  //obviously not all possible bad cases, but this is representative
  shapes[0].adjmatrix[0][0] = 1;
  rc = fc_segmentAdjacencyMatrix(shapes[0].numSides,shapes[0].adjmatrix,
				 &numSeg,&sidesPerSeg,NULL);
  fail_unless(rc != FC_SUCCESS, "should fail for bad adjmatrix");
  fail_unless(numSeg == -1, "should return -1 when fail");
  fail_unless(sidesPerSeg == NULL, "should return NULL when fail");
  fail_unless(segIDS == NULL, "should return NULL when fail");

  shapes[0].adjmatrix[0][0] = 0;
  shapes[0].adjmatrix[0][1] = (!shapes[0].adjmatrix[0][1]);
  rc = fc_segmentAdjacencyMatrix(shapes[0].numSides,shapes[0].adjmatrix,
				 &numSeg,&sidesPerSeg,NULL);
  fail_unless(rc != FC_SUCCESS, "should fail for bad adjmatrix");
  fail_unless(numSeg == -1, "should return -1 when fail");
  fail_unless(sidesPerSeg == NULL, "should return NULL when fail");
  fail_unless(segIDS == NULL, "should return NULL when fail");
  shapes[0].adjmatrix[0][1] = (!shapes[0].adjmatrix[0][1]);

  fc_freeShape(&shapes[0]);
  free(shapes);
  fc_deleteSubset(subset);


  fc_deleteDataset(dataset);
}
END_TEST


//this assumes that shape methods and reorder work
//(must be called after matrixshape, badargs, and simple).
// tests createScrewOrder 
START_TEST(screw){
  FC_ReturnCode rc;
  FC_Dataset dataset;
  FC_Mesh mesh, invalidMesh;
  FC_Shape *rawshapes;
  int numRawshapes;
  int *screworder, *areaorder;

  FC_Coords lowers = { 0., 0., 0. };
  FC_Coords uppers = { 1., 1., 1. };

  int i;

  // create test screw mesh
  rc = fc_createDataset("screw dataset", &dataset);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create dataset");
  rc = _fc_makeMeshFromFiles(dataset, "../data/screw.vert",
			     "../data/screw.elem", &mesh);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create mesh");

  //check inner methods explicitly
  //meshshapes is already tested
  rc = fc_getMeshShapes(mesh,80,1,&numRawshapes,&rawshapes);
  fail_unless(rc == FC_SUCCESS, "abort: failed to get shapes"); 
  fail_unless(numRawshapes == 1, "wrong number of shapes");
  fail_unless(rawshapes[0].numSides == 5, "wrong number of sides");

  rc = fc_createScrewShapeOrder(&rawshapes[0],&screworder);
  fail_unless(rc == FC_SUCCESS, "abort: failed to get screworder");
  rc = fc_createDescendingAreaOrder(&rawshapes[0],&areaorder);
  fail_unless(rc == FC_SUCCESS, "abort: failed to get area order");
  fail_unless(areaorder[screworder[0]] < areaorder[screworder[4]],
	      "head area should be smaller than base area");

  //we know the adj matrix is right, compare the order to it,
  //dont have to examine the whole adj matrix to test
  fail_unless(rawshapes[0].adjmatrix[screworder[0]][screworder[1]],
	      "head connectivity is wrong");
  fail_unless(rawshapes[0].adjmatrix[screworder[4]][screworder[3]],
	      "base connectivity is wrong");
  fail_unless(rawshapes[0].adjmatrix[screworder[2]][screworder[3]] &&
	      rawshapes[0].adjmatrix[screworder[2]][screworder[1]],
	      "mid plane connectivity is wrong");

  free(screworder);
  free(areaorder);
  for (i = 0; i < numRawshapes; i++){
    fc_freeShape(&rawshapes[i]);
  }
  free(rawshapes);

  //non screw mesh 
  rc = fc_createSimpleHexMesh(dataset, "junk", 1, 1, 1, lowers, uppers,
			      &invalidMesh);
  rc = fc_getMeshShapes(invalidMesh,80,1,&numRawshapes,&rawshapes);
  fail_unless(rc == FC_SUCCESS, "abort: failed to get shapes"); 
  fail_unless(numRawshapes == 1, "wrong number of shapes");
  fail_unless(rawshapes[0].numSides == 6, "wrong number of sides");
  rc = fc_createScrewShapeOrder(&rawshapes[0],&screworder);
  fail_unless(rc == FC_SUCCESS, "abort: should work for nonscrewmesh");
  fail_unless(screworder == NULL, "abort: return null order for nonscrewmesh");

  for (i = 0; i < numRawshapes; i++){
    fc_freeShape(&rawshapes[i]);
  }
  free(rawshapes);
  
  //bad args
  //now get the shape back
  rc = fc_getMeshShapes(mesh,80,1,&numRawshapes,&rawshapes);
  fail_unless(rc == FC_SUCCESS, "abort: failed to get shapes"); 
  fail_unless(numRawshapes == 1, "wrong number of shapes");
  fail_unless(rawshapes[0].numSides == 5, "wrong number of sides");

  //bad args for screw order
  rc = fc_createScrewShapeOrder(NULL,&screworder);
  fail_unless(rc != FC_SUCCESS, "should fail NULL shape");
  fail_unless(screworder == NULL, "should return NULL when fail");
  rc = fc_createScrewShapeOrder(&rawshapes[0],NULL);
  fail_unless(rc != FC_SUCCESS, "should fail NULL returnarg");

  //now change numside
  rawshapes[0].numSides = 3;
  rc = fc_createScrewShapeOrder(&rawshapes[0],&screworder);
  fail_unless(rc == FC_SUCCESS, 
	      "should work for nonscrew becuase of numsides");
  fail_unless(screworder == NULL, "should return NULL when fail");
  //fix it back
  rawshapes[0].numSides = 5;

  //mess with the adj matrix
  //im goint to use the areaorder matrix so i 
  //hang onto it even though its a screw order
  rc = fc_createScrewShapeOrder(&rawshapes[0],&screworder);
  fail_unless(rc == FC_SUCCESS, "cant get screworder");
  rawshapes[0].adjmatrix[screworder[0]][screworder[1]] = 0;
  free(screworder);
  rc = fc_createScrewShapeOrder(&rawshapes[0],&screworder);
  fail_unless(rc != FC_SUCCESS, "should fail for bad adjmatrix");
  fail_unless(screworder == NULL, "should return NULL when fail");

  for (i = 0; i < numRawshapes; i++){
    fc_freeShape(&rawshapes[i]);
  }
  free(rawshapes);

  fc_deleteDataset(dataset);

}
END_TEST


//ends assumes shape, screworder and matrix work (for the innards)
START_TEST(ends){
  FC_ReturnCode rc;
  FC_Dataset dataset;
  FC_Mesh mesh, hexMesh;
  FC_Shape *shapes, newshape;
  int numShapes;
  int numEnds;
  int *endArray;
  int *screworder;

  FC_Coords lowers = { 0., 0., 0. };
  FC_Coords uppers = { 1., 1., 1. };

  // create test screw mesh
  rc = fc_createDataset("screw dataset", &dataset);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create dataset");
  rc = _fc_makeMeshFromFiles(dataset, "../data/screw.vert",
			     "../data/screw.elem", &mesh);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create mesh");

  rc = fc_getMeshShapes(mesh,80,1,&numShapes,&shapes);
  fail_unless(rc == FC_SUCCESS, "abort: failed to get shapes"); 
  fail_unless(numShapes == 1, "wrong number of shapes");
  fail_unless(shapes[0].numSides == 5, "wrong number of sides");

  rc = fc_createScrewShapeOrder(&shapes[0],&screworder);
  fail_unless(rc == FC_SUCCESS, "abort: failed to get screworder");

  
  rc = fc_getShapeEnds(&shapes[0],&numEnds,&endArray);
  fail_unless(rc == FC_SUCCESS, "cant get screw ends");
  fail_unless(numEnds == 2, "wrong number of ends for screw");
  fail_unless(endArray[0] == screworder[4], "smallest screw end should be screw bottom"); 
  fail_unless(endArray[1] == screworder[0], "largest screw end should be screw top"); 
  free(endArray);
  free(screworder);


  //----------------special cases-------------------------------------//
  //single side
  newshape.numSides = 1;
  newshape.faces = &(shapes[0].faces[0]);
  newshape.elems = &(shapes[0].elems[0]);
  newshape.adjmatrix = (int**)malloc(sizeof(int*));
  fail_unless(newshape.adjmatrix != NULL, "memory error");
  newshape.adjmatrix[0] = (int*)malloc(sizeof(int));
  fail_unless(newshape.adjmatrix[0] != NULL, "memory error");
  newshape.adjmatrix[0][0] = 0;
  rc = fc_getShapeEnds(&newshape,&numEnds,&endArray);
  fail_unless(rc == FC_SUCCESS, "should work for 1 side ends");
  fail_unless(numEnds == 0, "1 side object should have 0 ends");
  fail_unless(endArray == NULL, "1 side should result in NULL end array");

  fc_freeShape(&shapes[0]);
  free(shapes);
  //newshape we have already freed the subsets
  free(newshape.adjmatrix[0]);
  free(newshape.adjmatrix);


  //no ends
  rc = fc_createSimpleHexMesh(dataset, "junk", 1, 1, 1, lowers, uppers,
			      &hexMesh);
  fail_unless(rc == FC_SUCCESS, "cant make hex mesh");
  rc = fc_getMeshShapes(hexMesh,80,0,&numShapes,&shapes);
  fail_unless(rc == FC_SUCCESS, "abort: failed to get shapes"); 
  fail_unless(numShapes == 1, "wrong number of shapes");
  fail_unless(shapes[0].numSides == 6, "wrong number of sides");

  rc = fc_getShapeEnds(&shapes[0],&numEnds,&endArray);
  fail_unless(rc == FC_SUCCESS, "should work for 0 side ends");
  fail_unless(numEnds == 0, "1 side object should have 0 ends");
  fail_unless(endArray == NULL, "1 side should result in NULL end array");


  //not checking results for shared_dim since thats just a passed parameter
  //to the adj matrix and its tested there

  //keep shape for bad args

  //bad args
  rc = fc_getShapeEnds(NULL,&numEnds,&endArray);
  fail_unless(rc != FC_SUCCESS, "should fail for NULL shape");
  fail_unless(numEnds == -1, "should return -1 when fail");
  fail_unless(endArray == NULL, "should return NULL when fail");

  rc = fc_getShapeEnds(&shapes[0],NULL,&endArray);
  fail_unless(rc != FC_SUCCESS, "should fail for NULL return arg");
  fail_unless(endArray == NULL, "should return NULL when fail");

  rc = fc_getShapeEnds(&shapes[0],&numEnds,NULL);
  fail_unless(rc != FC_SUCCESS, "should fail for NULL return arg");
  fail_unless(numEnds == -1, "should return -1 when fail");

  //make the shape bad
  shapes[0].numSides = 0;
  rc = fc_getShapeEnds(&shapes[0],&numEnds,&endArray);
  fail_unless(rc != FC_SUCCESS, "should fail for NULL no sides ");
  fail_unless(numEnds == -1, "should return -1 when fail");
  fail_unless(endArray == NULL, "should return NULL when fail");
  //fix it
  shapes[0].numSides = 6;

  //bad adj matrix
  shapes[0].adjmatrix[0][1] = 1;
  shapes[0].adjmatrix[1][0] = 0;
  rc = fc_getShapeEnds(&shapes[0],&numEnds,&endArray);
  fail_unless(rc != FC_SUCCESS, "should fail for bad adj matrix ");
  fail_unless(numEnds == -1, "should return -1 when fail");
  fail_unless(endArray == NULL, "should return NULL when fail");

  fc_freeShape(&shapes[0]);
  free(shapes);

  fc_deleteDataset(dataset);

}
END_TEST


// assumes that shape building works althouught this is kind of circular
// since the shape uses the normal. however i think it will
// all work out if i either test getOutwardFaceNormal explcitly
// or compare to known answers
START_TEST(normaltest){
  FC_ReturnCode rc;
  FC_Dataset dataset;
  FC_Mesh hexMesh;
  //the nums used in the hex mesh are used for comparison testing below
  FC_Coords lowers = { 0., 0., 0. };
  FC_Coords uppers = { 2., 1., 1. };

  FC_Shape *shapes;
  int numShapes;
  FC_Coords normal;
  FC_Vector std;
  int *numVertPerFace, maxNumVertPerFace;
  int *faceConns;
  double *mesh_coords;
  int i, j, k;

  rc = fc_createDataset("temp.xxx", &dataset);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create dataset");

  //the nums used in the hex mesh are used for comparison testing below
  rc = fc_createSimpleHexMesh(dataset, "junk", 2, 1, 1,
			      lowers, uppers, &hexMesh);
  fail_unless(rc == FC_SUCCESS, "cant make hex mesh");

  rc = fc_getMeshShapes(hexMesh,80,1,&numShapes,&shapes);
  fail_unless(rc == FC_SUCCESS, "abort: failed to get shapes"); 
  fail_unless(numShapes == 1, "wrong number of shapes");
  fail_unless(shapes[0].numSides == 6, "wrong number of sides");

  rc = fc_getMeshFaceConnsPtr(hexMesh, &numVertPerFace, &maxNumVertPerFace,
			      &faceConns);
  rc = fc_getMeshCoordsPtr(hexMesh,&mesh_coords);

  // --- normal shape --- //
  for (i = 0; i < shapes[0].numSides; i++){
    FC_Coords coords[4];
    FC_Coords v1, v2;
    int compval[3];
    int numface, *facearr;
    int multiplier = 1;

    rc = fc_getPlanarSideNormal(&shapes[0],i, &normal);
    fail_unless(rc == FC_SUCCESS, "cant get normal");

    rc = fc_getSubsetMembersAsArray(shapes[0].faces[i],&numface, &facearr);
    fail_unless(rc == FC_SUCCESS, "cant get elems as array");
    fail_unless(numface > 0, "no faces for this side");
    fail_unless(numVertPerFace[facearr[0]] == 4, "wrong num vert for this example");

    for (j = 0; j < numVertPerFace[facearr[0]]; j++){
      for (k = 0; k < 3; k++){
	coords[j][k] = mesh_coords[faceConns[facearr[0]*maxNumVertPerFace+j]*3+k];
      }
      if ((!(int)coords[j][0] && !(int)coords[j][1] && !(int)coords[j][2])){
	multiplier = -1; //known direction of normal
      }
    }

    for (j = 0; j < 3; j++){
      v1[j] = coords[2][j] - coords[1][j];
      v2[j] = coords[0][j] - coords[1][j];
    }
    
    compval[0] = multiplier*abs((int)(v1[1]*v2[2]-v1[2]*v2[1]));
    compval[1] = multiplier*abs((int)(v1[2]*v2[0]-v1[0]*v2[2]));
    compval[2] = multiplier*abs((int)(v1[0]*v2[1]-v1[1]*v2[0]));


    fail_unless((((int)normal[0] == compval[0]) &&
		 ((int)normal[1] == compval[1]) &&
		 ((int)normal[2] == compval[2])),
		"bad normal");
    
    free(facearr);
    rc = fc_getSideNormalMeanSdev(&shapes[0],i,
				  &normal,&std);
    fail_unless(rc == FC_SUCCESS, "cant get normal mean sdev");

    fail_unless((((int)normal[0] == compval[0]) &&
		 ((int)normal[1] == compval[1]) &&
		 ((int)normal[2] == compval[2])),
		"bad normal");
  }

  //keep the shape....

  //---- combo together two sides into one --- //
  {
    //have to manually manipulate the sides becuase want all this
    //checked before reshape. for ease, lets find 2 adj
    //sides each of which has 2 faces.
    FC_Shape copyShape;
    int numVals, *vals;
    FC_Coords normal1, normal2;
    int origside;
    int compside;

    rc = fc_copyShape(&shapes[0],&copyShape);
    fail_unless(rc == FC_SUCCESS, "cant copy shape");
    origside = -1;
    for (i = 0; i < shapes[0].numSides; i++){
      fc_getSubsetNumMember(copyShape.faces[i],&k);
      if (k == 2){
	origside = i;
	break;
      }
    }
    fail_unless(origside != -1, "can't find an orig side");

    compside = -1;
    for (i = 0; i < shapes[0].numSides; i++){
      if (copyShape.adjmatrix[0][i] == 1){
	fc_getSubsetNumMember(copyShape.faces[i],&k);
	if (k == 2){
	  compside = i;
	  break;
	}
      }
    }
    fail_unless(compside != -1, "can't find an adj side to the orig side");

    rc = fc_getSubsetMembersAsArray(copyShape.faces[compside],&numVals,&vals);
    fail_unless(rc == FC_SUCCESS, "cant get subset members as array");
    rc = fc_addArrayMembersToSubset(copyShape.faces[origside],numVals,vals);
    fail_unless(rc == FC_SUCCESS, "cant add array members to subset");
    free(vals);
    rc = fc_getSubsetMembersAsArray(copyShape.elems[compside],&numVals,&vals);
    fail_unless(rc == FC_SUCCESS, "cant get subset members as array");
    rc = fc_addArrayMembersToSubset(copyShape.elems[origside],numVals,vals);
    fail_unless(rc == FC_SUCCESS, "cant add array members to subset");
    free(vals);

    //we know at this point that sidenormalmeansdev works for the 
    //planar sides individually
    rc = fc_getSideNormalMeanSdev(&shapes[0],origside,&normal1,&std);
    fail_unless(rc == FC_SUCCESS, "cant get normal mean sdev");
    rc = fc_getSideNormalMeanSdev(&shapes[0],compside,&normal2,&std);
    fail_unless(rc == FC_SUCCESS, "cant get normal mean sdev");

    rc = fc_getSideNormalMeanSdev(&copyShape,origside,&normal,&std);
    fail_unless(rc == FC_SUCCESS, "cant get normal mean sdev");
    for (i = 0; i < 3; i++){
      if (fabs(normal1[i]) > 0){
	fail_unless(FC_FLT_EQUIV(normal[i], 0.5*normal1[i]), "wrong val for normal");
	fail_unless(FC_FLT_EQUIV(std[i],sqrt(1.0/3.0)), "wrong val for std");
      }else if (fabs(normal2[i]) > 0){
	fail_unless(FC_FLT_EQUIV(normal[i],0.5*normal2[i]), "wrong val for normal");
	fail_unless(FC_FLT_EQUIV(std[i],sqrt(1.0/3.0)), "wrong val for std");
      }else{
	fail_unless(FC_FLT_EQUIV(normal[i],0), "wrong val for normal");
	fail_unless(FC_FLT_EQUIV(std[i],0), "wrong val for std");
      }
    }
    fc_freeShape(&copyShape);
  }

  fc_freeShape(&shapes[0]);
  free(shapes);

  //------now try with one giant side--------/
  rc = fc_getMeshShapes(hexMesh,100,1,&numShapes,&shapes);
  fail_unless(rc == FC_SUCCESS, "abort: failed to get shapes"); 
  fail_unless(numShapes == 1, "wrong number of shapes");
  fail_unless(shapes[0].numSides == 1, "wrong number of sides");
  rc = fc_getSideNormalMeanSdev(&shapes[0],0,
				&normal,&std);
  fail_unless(rc == FC_SUCCESS, "cant get normal mean sdev");
  fail_unless(FC_FLT_EQUIV(normal[0],0.0) &&
	      FC_FLT_EQUIV(normal[1],0.0) &&
	      FC_FLT_EQUIV(normal[2],0.0),
	      "bad normal");
  //known vals from it being rectangular and 2 divisions in x, 1 in y and z
  fail_unless(FC_FLT_EQUIV(std[0],sqrt(2.0/9.0)) &&
	      FC_FLT_EQUIV(std[1],sqrt(4.0/9.0)) &&
	      FC_FLT_EQUIV(std[2],sqrt(4.0/9.0)),
	      "bad std");


  fc_freeShape(&shapes[0]);
  free(shapes);

  //now get our normal shape back
  rc = fc_getMeshShapes(hexMesh,80,1,&numShapes,&shapes);
  fail_unless(rc == FC_SUCCESS, "abort: failed to get shapes"); 
  fail_unless(numShapes == 1, "wrong number of shapes");
  fail_unless(shapes[0].numSides == 6, "wrong number of sides");

  //bad args
  rc = fc_getPlanarSideNormal(NULL, 1, &normal);
  fail_unless(rc != FC_SUCCESS, "should fail for NULL shape");
  for (i = 0; i < 3; i++){
    fail_unless(normal[i] == 0, "should return 0 when fail");
  }

  rc = fc_getSideNormalMeanSdev(NULL, 1, &normal, &std);
  fail_unless(rc != FC_SUCCESS, "should fail for NULL shape");
  for (i = 0; i < 3; i++){
    fail_unless(normal[i] == 0, "should return 0 when fail");
    fail_unless(std[i] == 0, "should return 0 when fail");
  }

  rc = fc_getPlanarSideNormal(&shapes[0], -1, &normal);
  fail_unless(rc != FC_SUCCESS, "should fail for bad side id");
  for (i = 0; i < 3; i++){
    fail_unless(normal[i] == 0, "should return 0 when fail");
  }

  rc = fc_getSideNormalMeanSdev(&shapes[0], -1, &normal, &std);
  fail_unless(rc != FC_SUCCESS, "should fail for bad side id");
  for (i = 0; i < 3; i++){
    fail_unless(normal[i] == 0, "should return 0 when fail");
    fail_unless(std[i] == 0, "should return 0 when fail");
  }

  rc = fc_getPlanarSideNormal(&shapes[0], 7, &normal);
  fail_unless(rc != FC_SUCCESS, "should fail for bad side id");
  for (i = 0; i < 3; i++){
    fail_unless(normal[i] == 0, "should return 0 when fail");
  }

  rc = fc_getSideNormalMeanSdev(&shapes[0], 7, &normal, &std);
  fail_unless(rc != FC_SUCCESS, "should fail for bad side id");
  for (i = 0; i < 3; i++){
    fail_unless(normal[i] == 0, "should return 0 when fail");
    fail_unless(std[i] == 0, "should return 0 when fail");
  }

  rc = fc_getPlanarSideNormal(&shapes[0], 1, NULL);
  fail_unless(rc != FC_SUCCESS, "should fail for NULL return arg");
  for (i = 0; i < 3; i++){
    fail_unless(normal[i] == 0, "should return 0 when fail");
  }

  rc = fc_getSideNormalMeanSdev(&shapes[0], 1, NULL, &std);
  fail_unless(rc != FC_SUCCESS, "should fail for NULL return arg");
  for (i = 0; i < 3; i++){
    fail_unless(std[i] == 0, "should return 0 when fail");
  }

  rc = fc_getSideNormalMeanSdev(&shapes[0], 1, &normal, NULL);
  fail_unless(rc != FC_SUCCESS, "should fail for NULL return arg");
  for (i = 0; i < 3; i++){
    fail_unless(normal[i] == 0, "should return 0 when fail");
  }

  for (i = 0; i < numShapes; i++){
    fc_freeShape(&shapes[i]);
  }
  free(shapes);

  fc_deleteDataset(dataset);
}
END_TEST

// This testes the bad cases of shape:
// 1) subsetShape 2) meshShape 
// Things are sufficiently complex that i prefer to 
// separate out the bad cases from the good ones.
// Also the good cases are fairly well linked with the 
// matrix cals, so i will put the good tests there
//
// 08/19/08 ACG subsetSides has been moved into
// getSubsetShapes, so it is no longer tested independently.
START_TEST(shapesides_bad){
  FC_ReturnCode rc;
  FC_Dataset dataset;
  FC_Subset subset;
  FC_Mesh mesh;

  //generating meshes and subsets to look at
  FC_Coords lowers = { 0., 0., 0. };
  FC_Coords uppers = { 1., 1., 1. };

  //generating meshes and subsets to look at
  int elemIDS[2];
  int numElem;

  FC_Subset emptySubset;
  FC_Mesh badMesh = {999,999};
  FC_Subset badSubset= {999,999};
  FC_Subset wrongassoc;


  //for shape
  FC_Shape *shapes;
  int numShapes;

  // setup - create test dataset
  rc = fc_createDataset("temp.xxx", &dataset);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create dataset");

  rc = fc_createSimpleHexMesh(dataset, "junk", 2, 2, 2, lowers, uppers,
				&mesh);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create simple hex mesh");
  rc = fc_createSubset(mesh,"temp",FC_AT_FACE,&emptySubset);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create subset");


  numElem = 2;
  elemIDS[0] = 0;
  elemIDS[1] = 2;


  rc = fc_createSubset(mesh,"temp",FC_AT_ELEMENT,&subset);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create subset");

  rc = fc_addArrayMembersToSubset(subset,numElem,elemIDS);
  fail_unless(rc == FC_SUCCESS, "abort: failed to add members to subset");

  rc = fc_getSubsetShapes(FC_NULL_SUBSET,30,1,&numShapes,&shapes);
  fail_unless(rc != FC_SUCCESS, "Should fail for NULL subset");
  fail_unless(numShapes == -1, "fail should return -1");
  fail_unless(shapes == NULL, "fail should return NULL");

  rc = fc_getSubsetShapes(emptySubset,30,1,&numShapes,&shapes);
  fail_unless(rc != FC_SUCCESS, "Should fail for empty subset");
  fail_unless(numShapes == -1, "fail should return -1");
  fail_unless(shapes == NULL, "fail should return NULL");

  rc = fc_getSubsetShapes(badSubset,30,1,&numShapes,&shapes);
  fail_unless(rc != FC_SUCCESS, "Should fail for bad subset");
  fail_unless(numShapes == -1, "fail should return -1");
  fail_unless(shapes == NULL, "fail should return NULL");

  rc = fc_getMeshShapes(badMesh,30,1,&numShapes,&shapes);
  fail_unless(rc != FC_SUCCESS, "Should fail for bad mesh");
  fail_unless(numShapes == -1, "fail should return -1");
  fail_unless(shapes == NULL, "fail should return NULL");

  rc = fc_getSubsetShapes(subset,-1,1,&numShapes,&shapes);
  fail_unless(rc != FC_SUCCESS, "Should fail for bad angle");
  fail_unless(numShapes == -1, "fail should return -1");
  fail_unless(shapes == NULL, "fail should return NULL");

  rc = fc_getMeshShapes(mesh,-1,1,&numShapes,&shapes);
  fail_unless(rc != FC_SUCCESS, "Should fail for bad angle");
  fail_unless(numShapes == -1, "fail should return -1");
  fail_unless(shapes == NULL, "fail should return NULL");

  rc = fc_getSubsetShapes(subset,181,1,&numShapes,&shapes);
  fail_unless(rc != FC_SUCCESS, "Should fail for bad angle");
  fail_unless(numShapes == -1, "fail should return -1");
  fail_unless(shapes == NULL, "fail should return NULL");

  rc = fc_getMeshShapes(mesh,181,1,&numShapes,&shapes);
  fail_unless(rc != FC_SUCCESS, "Should fail for bad angle");
  fail_unless(numShapes == -1, "fail should return -1");
  fail_unless(shapes == NULL, "fail should return NULL");

  rc = fc_getSubsetShapes(subset,30,-1,&numShapes,&shapes);
  fail_unless(rc != FC_SUCCESS, "Should fail for bad neighbordim");
  fail_unless(numShapes == -1, "fail should return -1");
  fail_unless(shapes == NULL, "fail should return NULL");

  rc = fc_getMeshShapes(mesh,30,-1,&numShapes,&shapes);
  fail_unless(rc != FC_SUCCESS, "Should fail for bad neighbordim");
  fail_unless(numShapes == -1, "fail should return -1");
  fail_unless(shapes == NULL, "fail should return NULL");

  rc = fc_getSubsetShapes(subset,30,2,&numShapes,&shapes);
  fail_unless(rc != FC_SUCCESS, "Should fail for bad neighbordim");
  fail_unless(numShapes == -1, "fail should return -1");
  fail_unless(shapes == NULL, "fail should return NULL");

  rc = fc_getMeshShapes(mesh,30,2,&numShapes,&shapes);
  fail_unless(rc != FC_SUCCESS, "Should fail for bad neighbordim");
  fail_unless(numShapes == -1, "fail should return -1");
  fail_unless(shapes == NULL, "fail should return NULL");

  rc = fc_getSubsetShapes(subset,30,1,NULL,&shapes);
  fail_unless(rc != FC_SUCCESS, "Should fail for NULL arg");
  fail_unless(shapes == NULL, "fail should return NULL");
  
  rc = fc_getSubsetShapes(subset,30,1,&numShapes,NULL);
  fail_unless(rc != FC_SUCCESS, "Should fail for NULL arg");
  fail_unless(numShapes == -1, "fail should return -1");

  //and lets try a subset with the wrong assoc
  rc = fc_createSubset(mesh,"junk2",FC_AT_FACE,&wrongassoc);
  fail_unless(rc == FC_SUCCESS, "cant make new assoc to test");
  rc = fc_addMemberToSubset(wrongassoc,1);
  fail_unless(rc == FC_SUCCESS, "cant add member to subset");

  rc = fc_getSubsetShapes(wrongassoc,30,1,&numShapes,&shapes);
  fail_unless(rc != FC_SUCCESS, "Should fail for NULL arg");
  fail_unless(shapes == NULL, "fail should return NULL");
  fc_deleteSubset(wrongassoc);

  //test meshshapes on wrong meshes below

  fc_deleteSubset(subset);
  fc_deleteDataset(dataset);


  { //bad mesh types - reusing mesh set up from something WSD had in one
    //of her tests, doesnt matter what the vals are, i just want some
    //meshes...note: the hex mesh should work, but im not testing
    //against it here, but its too hard to rip it out of the conns list.

    int i;
    FC_Mesh hexmesh;

    int numDim = 3;
    int numElemType = 4;
    FC_ElementType elemTypes[4] = { FC_ET_POINT,  FC_ET_LINE,
				    FC_ET_QUAD,   FC_ET_HEX };
    
    int numVertPerType[4] = { 4, 4, 16, 64}, numElemPerType[4] = { 4, 3, 9, 27 };
    int conns[4][64*27] = { { 0,       1,       2,     3},
			    { 0, 1,    1, 2,    2, 3},
			    { 0, 1, 5, 4,    1, 2, 6, 5,    2, 3, 7, 6,
			      4, 5, 9, 8,    5, 6, 10, 9,   6, 7, 11, 10,
			      8, 9, 13, 12,  9, 10, 14,13,  10, 11,15, 14
			    },
			    { } };
    double* coords;
    //lowers and uppers havent changed
    //    FC_Coords lowers = { 0., 0., 0. };
    //    FC_Coords uppers = { 1., 1., 1. };


    // setup - create test dataset
    rc = fc_createDataset("temp.xxx", &dataset);
    fail_unless(rc == FC_SUCCESS, "abort: failed to create dataset");
    
    // Make the hex mesh, then reuse coords for other meshes  
    rc = fc_createSimpleHexMesh(dataset, "hex", 3, 3, 3, lowers, uppers,
			      &hexmesh);
    fail_unless(rc == FC_SUCCESS, "abort: failed to create simple hex mesh");
    fc_getMeshCoordsPtr(hexmesh, &coords);
    for (i = 0; i < numElemType-1; i++) {
      FC_Mesh meshx;
      rc = fc_createMesh(dataset, "testmesh", &meshx);
      fail_unless(rc == FC_SUCCESS, "abort: failed to created mesh");
      rc = fc_setMeshCoords(meshx, numDim, numVertPerType[i], coords);
      fail_unless(rc == FC_SUCCESS, "failed to set vertex coords");
      rc = fc_setMeshElementConns(meshx, elemTypes[i], numElemPerType[i], 
				  conns[i]);
      fail_unless(rc == FC_SUCCESS, "failed to set element conns");
      
      rc = fc_createSubset(meshx, "junk", FC_AT_ELEMENT, &subset);
      fail_unless(rc == FC_SUCCESS, "failed to create subset");
      rc = fc_addMemberToSubset(subset,1);
      fail_unless(rc == FC_SUCCESS, "cant add member to subset");

      rc = fc_getSubsetShapes(subset,30,1,&numShapes,&shapes);
      fail_unless(rc != FC_SUCCESS, "Should fail for bad mesh type");
      fail_unless(numShapes == -1, "fail should return -1");
      fail_unless(shapes == NULL, "fail should return NULL");

      rc = fc_getMeshShapes(meshx,30,1,&numShapes,&shapes);
      fail_unless(rc != FC_SUCCESS, "Should fail for bad mesh type");
      fail_unless(numShapes == -1, "fail should return -1");
      fail_unless(shapes == NULL, "fail should return NULL");

      //dont have to free the shapes becuase tey never get created

      fc_deleteSubset(subset);
      fc_deleteMesh(meshx);
    }
    fc_deleteMesh(hexmesh);
    fc_deleteDataset(dataset);
  }

}
END_TEST


START_TEST(reduce){
  FC_ReturnCode rc;
  FC_Dataset dataset;
  FC_Mesh mesh;

  FC_Coords lowers = { 0., 0., 0. };
  FC_Coords uppers = { 1., 0.75, 0.5 };

  FC_Subset temp_subset;
  FC_Shape *shapes,temp_shape;
  int numShapes, tempnum, *temparr;
  int numcombomem, numcheckmem, comboside;
  int numReducedSides;
  char *name1, *name2;
  int i,j;

  // setup - create test dataset
  rc = fc_createDataset("temp.xxx", &dataset);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create dataset");
  
  rc = fc_createSimpleHexMesh(dataset, "junk", 2, 2, 3, lowers, uppers,
			      &mesh);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create simple hex mesh");
  
  rc = fc_getMeshShapes(mesh,80,0,&numShapes,&shapes);
  fail_unless(rc == FC_SUCCESS, "abort: failed to get MeshSides");
  fail_unless(numShapes == 1, "abort: bad num of shapes");
  fail_unless(shapes[0].numSides == 6, "abort: bad num of MeshSides");


  // ----- reduceShapeNumSides ---- //

  //we cant try all possible cases, but we will try several
  //the tests are set up specifically for the cases given
  //and if parameters change (numReducedSides) or the
  //shape orders, it will break the tests

  //1 side merge 
  //  fc_printShape(&shapes[0],"before",1);
  numReducedSides = shapes[0].numSides-1;
  rc = fc_reduceShapeNumSides(&shapes[0],numReducedSides,&temp_shape);
  fail_unless(rc == FC_SUCCESS, "cant reduce numsides");
  fail_unless(temp_shape.numSides == numReducedSides, "wrong num sides");
  //  fc_printShape(&temp_shape,"after",1);

  //determine the side that becomes the combo side
  //this will change if the algorithm changes
  comboside = -1;
  for (i = 0; i < shapes[0].numSides-2; i++){
    if (shapes[0].adjmatrix[i][shapes[0].numSides-1] == 1){
      comboside = i;
      break;
    }
  }
  fail_unless(comboside != -1, 
	      "abort: can't determine combo side to check");
  //  printf("combo side = %d\n",comboside);
  //faces will have no duplicates
  rc = fc_copySubset(shapes[0].faces[comboside],mesh,
		     "temp",&temp_subset);
  fail_unless(rc == FC_SUCCESS, "cant copy subset");
  rc = fc_getSubsetMembersAsArray(shapes[0].faces[shapes[0].numSides-1],
				  &tempnum, &temparr);
  fail_unless(rc == FC_SUCCESS, "cant get subset members");
  rc = fc_addArrayMembersToSubset(temp_subset,tempnum,temparr);
  free(temparr);
  fail_unless(rc == FC_SUCCESS, "cant add array members to subset");
  rc = fc_getSubsetNumMember(temp_shape.faces[comboside],
				  &numcheckmem);
  fail_unless(rc == FC_SUCCESS, "cant get subset members");
  rc = fc_getSubsetNumMember(temp_subset, &numcombomem);
  fail_unless(rc == FC_SUCCESS, "cant get subset members");
  fail_unless(numcheckmem == numcombomem, "wrong number of combo mem");
  //since already checked equal members, just have to check superset in either order
  rc = fc_isSubsetSuperset(temp_subset,
			   temp_shape.faces[comboside],&i);
  fail_unless(rc == FC_SUCCESS, "cant determine superset");
  fail_unless(i == 1, "wrong members in combo side");
  fail_unless(rc == FC_SUCCESS, "cant get subset members");
  fc_deleteSubset(temp_subset);

  //elements will have duplicates
  rc = fc_copySubset(shapes[0].elems[comboside],mesh,
		     "temp",&temp_subset);
  fail_unless(rc == FC_SUCCESS, "cant copy subset");
  rc = fc_getSubsetMembersAsArray(shapes[0].elems[shapes[0].numSides-1],
				  &tempnum, &temparr);
  fail_unless(rc == FC_SUCCESS, "cant get subset members");
  rc = fc_addArrayMembersToSubset(temp_subset,tempnum,temparr);
  free(temparr);
  fail_unless(rc == FC_SUCCESS, "cant add array members to subset");
  rc = fc_getSubsetNumMember(temp_shape.elems[comboside],
				  &numcheckmem);
  fail_unless(rc == FC_SUCCESS, "cant get subset members");
  rc = fc_getSubsetNumMember(temp_subset, &numcombomem);
  fail_unless(rc == FC_SUCCESS, "cant get subset members");
  fail_unless(numcheckmem == numcombomem, "wrong number of combo mem");
  //since already checked equal members, just have to check superset in either order
  rc = fc_isSubsetSuperset(temp_subset,
			   temp_shape.elems[comboside],&i);
  fail_unless(rc == FC_SUCCESS, "cant determine superset");
  fail_unless(i == 1, "wrong members in combo side");
  fail_unless(rc == FC_SUCCESS, "cant get subset members");
  fc_deleteSubset(temp_subset);

  //make sure all the other sides havent changed and that all names are correct
  for (i = 0; i < numReducedSides; i++){
    if (i != comboside){
      int check1, check2;
      rc = fc_isSubsetSuperset(shapes[0].faces[i],temp_shape.faces[i],&check1);
      fail_unless(rc == FC_SUCCESS, "cant determine superset");
      rc = fc_isSubsetSuperset(temp_shape.faces[i],shapes[0].faces[i],&check2);
      fail_unless(rc == FC_SUCCESS, "cant determine superset");
      fail_unless(check1 && check2, "non combo sides shouldnt have changed");
      rc = fc_isSubsetSuperset(shapes[0].elems[i],temp_shape.elems[i],&check1);
      fail_unless(rc == FC_SUCCESS, "cant determine superset");
      rc = fc_isSubsetSuperset(temp_shape.elems[i],shapes[0].elems[i],&check2);
      fail_unless(rc == FC_SUCCESS, "cant determine superset");
      fail_unless(check1 && check2, "non combo sides shouldnt have changed");
    }
    
    rc = fc_getSubsetName(shapes[0].faces[0],&name1);
    fail_unless(rc == FC_SUCCESS, "cant get subset name");
    rc = fc_getSubsetName(temp_shape.faces[0],&name2);
    fail_unless(rc == FC_SUCCESS, "cant get subset name");
    fail_unless(!strcmp(name1,name2), "incorrect name for subset");
    free(name1);
    free(name2);
    rc = fc_getSubsetName(shapes[0].elems[0],&name1);
    fail_unless(rc == FC_SUCCESS, "cant get subset name");
    rc = fc_getSubsetName(temp_shape.elems[0],&name2);
    fail_unless(rc == FC_SUCCESS, "cant get subset name");
    fail_unless(!strcmp(name1,name2), "incorrect name for subset");
    free(name1);
    free(name2);
  }

  for (i = 0; i < temp_shape.numSides; i++){
    for (j = 0; j < temp_shape.numSides; j++){
      if (i == j){
	fail_unless(temp_shape.adjmatrix[i][j] == 0,
		    "wrong adj matrix");
      }else if (shapes[0].adjmatrix[i][j] == 1){
	fail_unless(temp_shape.adjmatrix[i][j] == 1,
		    "wrong adj matrix");
      }else{
	//they may not have been adj before they could be
	//now if the combo made them adj
	if (i == comboside){
	  fail_unless(temp_shape.adjmatrix[i][j] == shapes[0].adjmatrix[shapes[0].numSides-1][j],
		      "wrong adj matrix");
	}else if (j == comboside){
	  fail_unless(temp_shape.adjmatrix[i][j] == shapes[0].adjmatrix[i][shapes[0].numSides-1],
		      "wrong adj matrix");
	}else{
	  fail_unless(temp_shape.adjmatrix[i][j] == shapes[0].adjmatrix[i][j],
		      "wrong adj matrix");
	}
      }
    }
  }
  
  rc = fc_freeShape(&temp_shape);
  fail_unless(rc == FC_SUCCESS, "problem freeing shape");


  //two side merge where it requires two passes
  //first lets mutate the shape - make the last side only
  //adjacent to the side prev to it
  //  fc_printShape(&shapes[0],"before 1",1);
  for (i = 0; i < shapes[0].numSides; i++){
    shapes[0].adjmatrix[i][shapes[0].numSides-1] = 0;
    shapes[0].adjmatrix[shapes[0].numSides-1][i] = 0;
  }
  shapes[0].adjmatrix[shapes[0].numSides-2][shapes[0].numSides-1] = 1;
  shapes[0].adjmatrix[shapes[0].numSides-1][shapes[0].numSides-2] = 1;

  //  fc_printShape(&shapes[0],"before 2",1);
  numReducedSides = shapes[0].numSides-2;
  rc = fc_reduceShapeNumSides(&shapes[0],numReducedSides,&temp_shape);
  fail_unless(rc == FC_SUCCESS, "cant reduce numsides");
  fail_unless(temp_shape.numSides == numReducedSides, "wrong num sides");
  //  fc_printShape(&temp_shape,"after",1);

  //determine the side that becomes the combo side
  //this will change if the algorithm changes
  comboside = -1;
  for (i = 0; i < shapes[0].numSides-2; i++){
    if (shapes[0].adjmatrix[i][shapes[0].numSides-2] == 1){
      comboside = i;
      break;
    }
  }
  fail_unless(comboside != -1, 
	      "abort: can't determine combo side to check");
  //  printf("combo side = %d\n",comboside);
  //faces will have no duplicates
  rc = fc_copySubset(shapes[0].faces[comboside],mesh,
		     "temp",&temp_subset);
  fail_unless(rc == FC_SUCCESS, "cant copy subset");
  rc = fc_getSubsetMembersAsArray(shapes[0].faces[shapes[0].numSides-1],
				  &tempnum, &temparr);
  fail_unless(rc == FC_SUCCESS, "cant get subset members");
  rc = fc_addArrayMembersToSubset(temp_subset,tempnum,temparr);
  free(temparr);
  rc = fc_getSubsetMembersAsArray(shapes[0].faces[shapes[0].numSides-2],
				  &tempnum, &temparr);
  fail_unless(rc == FC_SUCCESS, "cant get subset members");
  rc = fc_addArrayMembersToSubset(temp_subset,tempnum,temparr);
  free(temparr);
  fail_unless(rc == FC_SUCCESS, "cant add array members to subset");
  rc = fc_getSubsetNumMember(temp_shape.faces[comboside],
				  &numcheckmem);
  fail_unless(rc == FC_SUCCESS, "cant get subset members");
  rc = fc_getSubsetNumMember(temp_subset, &numcombomem);
  fail_unless(rc == FC_SUCCESS, "cant get subset members");
  fail_unless(numcheckmem == numcombomem, "wrong number of combo mem");
  //since already checked equal members, just have to check superset in either order
  rc = fc_isSubsetSuperset(temp_subset,
			   temp_shape.faces[comboside],&i);
  fail_unless(rc == FC_SUCCESS, "cant determine superset");
  fail_unless(i == 1, "wrong members in combo side");
  fail_unless(rc == FC_SUCCESS, "cant get subset members");
  fc_deleteSubset(temp_subset);

  //elements will have duplicates
  rc = fc_copySubset(shapes[0].elems[comboside],mesh,
		     "temp",&temp_subset);
  fail_unless(rc == FC_SUCCESS, "cant copy subset");
  rc = fc_getSubsetMembersAsArray(shapes[0].elems[shapes[0].numSides-1],
				  &tempnum, &temparr);
  fail_unless(rc == FC_SUCCESS, "cant get subset members");
  rc = fc_addArrayMembersToSubset(temp_subset,tempnum,temparr);
  free(temparr);
  rc = fc_getSubsetMembersAsArray(shapes[0].elems[shapes[0].numSides-2],
				  &tempnum, &temparr);
  fail_unless(rc == FC_SUCCESS, "cant get subset members");
  rc = fc_addArrayMembersToSubset(temp_subset,tempnum,temparr);
  free(temparr);
  fail_unless(rc == FC_SUCCESS, "cant add array members to subset");
  rc = fc_getSubsetNumMember(temp_shape.elems[comboside],
				  &numcheckmem);
  fail_unless(rc == FC_SUCCESS, "cant get subset members");
  rc = fc_getSubsetNumMember(temp_subset, &numcombomem);
  fail_unless(rc == FC_SUCCESS, "cant get subset members");
  fail_unless(numcheckmem == numcombomem, "wrong number of combo mem");
  //since already checked equal members, just have to check superset in either order
  rc = fc_isSubsetSuperset(temp_subset,
			   temp_shape.elems[comboside],&i);
  fail_unless(rc == FC_SUCCESS, "cant determine superset");
  fail_unless(i == 1, "wrong members in combo side");
  fail_unless(rc == FC_SUCCESS, "cant get subset members");
  fc_deleteSubset(temp_subset);

  //make sure all the other sides havent changed and that all names are correct
  for (i = 0; i < numReducedSides; i++){
    if (i != comboside){
      int check1, check2;
      rc = fc_isSubsetSuperset(shapes[0].faces[i],temp_shape.faces[i],&check1);
      fail_unless(rc == FC_SUCCESS, "cant determine superset");
      rc = fc_isSubsetSuperset(temp_shape.faces[i],shapes[0].faces[i],&check2);
      fail_unless(rc == FC_SUCCESS, "cant determine superset");
      fail_unless(check1 && check2, "non combo sides shouldnt have changed");
      rc = fc_isSubsetSuperset(shapes[0].elems[i],temp_shape.elems[i],&check1);
      fail_unless(rc == FC_SUCCESS, "cant determine superset");
      rc = fc_isSubsetSuperset(temp_shape.elems[i],shapes[0].elems[i],&check2);
      fail_unless(rc == FC_SUCCESS, "cant determine superset");
      fail_unless(check1 && check2, "non combo sides shouldnt have changed");
    }
    rc = fc_getSubsetName(shapes[0].faces[0],&name1);
    fail_unless(rc == FC_SUCCESS, "cant get subset name");
    rc = fc_getSubsetName(temp_shape.faces[0],&name2);
    fail_unless(rc == FC_SUCCESS, "cant get subset name");
    fail_unless(!strcmp(name1,name2), "incorrect name for subset");
    free(name1);
    free(name2);
    rc = fc_getSubsetName(shapes[0].elems[0],&name1);
    fail_unless(rc == FC_SUCCESS, "cant get subset name");
    rc = fc_getSubsetName(temp_shape.elems[0],&name2);
    fail_unless(rc == FC_SUCCESS, "cant get subset name");
    fail_unless(!strcmp(name1,name2), "incorrect name for subset");
    free(name1);
    free(name2);
  }

  for (i = 0; i < temp_shape.numSides; i++){
    for (j = 0; j < temp_shape.numSides; j++){
      if (i == j){
	fail_unless(temp_shape.adjmatrix[i][j] == 0,
		    "wrong adj matrix");
      }else if (shapes[0].adjmatrix[i][j] == 1){
	fail_unless(temp_shape.adjmatrix[i][j] == 1,
		    "wrong adj matrix");
      }else{
	//they may not have been adj before they could be
	//now if the combo made them adj
	if (i == comboside){
	  fail_unless(temp_shape.adjmatrix[i][j] == 
		      (shapes[0].adjmatrix[shapes[0].numSides-1][j] ||
		       shapes[0].adjmatrix[shapes[0].numSides-2][j]),
		      "wrong adj matrix");
	}else if (j == comboside){
	  fail_unless(temp_shape.adjmatrix[i][j] == 
		      (shapes[0].adjmatrix[i][shapes[0].numSides-1] ||
		       shapes[0].adjmatrix[i][shapes[0].numSides-2]),
		      "wrong adj matrix");
	}else{
	  fail_unless(temp_shape.adjmatrix[i][j] == shapes[0].adjmatrix[i][j],
		      "wrong adj matrix");
	}
      }
    }
  }
  
  rc = fc_freeShape(&temp_shape);
  fail_unless(rc == FC_SUCCESS, "problem freeing shape");

  //bad case where there is no adj side
  //mutate the shape once more
  for (i = 0; i < shapes[0].numSides; i++){
    shapes[0].adjmatrix[i][shapes[0].numSides-1] = 0;
    shapes[0].adjmatrix[shapes[0].numSides-1][i] = 0;
  }

  numReducedSides = shapes[0].numSides-1;
  rc = fc_reduceShapeNumSides(&shapes[0],numReducedSides,&temp_shape);
  fail_unless(rc != FC_SUCCESS, "should fail for non-adj side in the shape");

  //get a pure shape[0] for now
  for (i = 0; i < numShapes; i++){
    fc_freeShape(&shapes[i]);
  }
  free(shapes);
  rc = fc_getMeshShapes(mesh,80,0,&numShapes,&shapes);
  fail_unless(rc == FC_SUCCESS, "abort: failed to get MeshSides");
  fail_unless(numShapes == 1, "abort: bad num of shapes");
  fail_unless(shapes[0].numSides == 6, "abort: bad num of MeshSides");

  //bad args for reduce numsides
  temp_shape.numSides =3;
  rc = fc_reduceShapeNumSides(NULL,2,&temp_shape);
  fail_unless(rc != FC_SUCCESS, "should fail for NULL shape");
  fail_unless(temp_shape.numSides == 0, "shape should return as init for fail");
  fail_unless(temp_shape.faces == NULL, "shape should return as init for fail");
  fail_unless(temp_shape.elems == NULL, "shape should return as init for fail"); 
  fail_unless(temp_shape.adjmatrix == NULL, "shape should return as init for fail"); 

  temp_shape.numSides =3;
  rc = fc_reduceShapeNumSides(&shapes[0],-1,&temp_shape);
  fail_unless(rc != FC_SUCCESS, "should fail for bad numSides");
  fail_unless(temp_shape.numSides == 0, "shape should return as init for fail");
  fail_unless(temp_shape.faces == NULL, "shape should return as init for fail");
  fail_unless(temp_shape.elems == NULL, "shape should return as init for fail"); 
  fail_unless(temp_shape.adjmatrix == NULL, "shape should return as init for fail"); 

  temp_shape.numSides =3;
  rc = fc_reduceShapeNumSides(&shapes[0],0,&temp_shape);
  fail_unless(rc != FC_SUCCESS, "should fail for bad numSides");
  fail_unless(temp_shape.numSides == 0, "shape should return as init for fail");
  fail_unless(temp_shape.faces == NULL, "shape should return as init for fail");
  fail_unless(temp_shape.elems == NULL, "shape should return as init for fail"); 
  fail_unless(temp_shape.adjmatrix == NULL, "shape should return as init for fail"); 

  temp_shape.numSides =3;
  rc = fc_reduceShapeNumSides(&shapes[0],shapes[0].numSides,&temp_shape);
  fail_unless(rc != FC_SUCCESS, "should fail for bad numSides");
  fail_unless(temp_shape.numSides == 0, "shape should return as init for fail");
  fail_unless(temp_shape.faces == NULL, "shape should return as init for fail");
  fail_unless(temp_shape.elems == NULL, "shape should return as init for fail"); 
  fail_unless(temp_shape.adjmatrix == NULL, "shape should return as init for fail"); 

  rc = fc_reduceShapeNumSides(&shapes[0],2,NULL);
  fail_unless(rc != FC_SUCCESS, "should fail for NULL arg");

  for (i = 0; i < numShapes; i++){
    fc_freeShape(&shapes[i]);
  }
  free(shapes);
}
END_TEST


START_TEST(reshape){
  FC_ReturnCode rc;
  FC_Dataset dataset;
  FC_Mesh mesh;
  FC_Shape *shapes, temp_shape, *testshapes;
  int numShapes, numTestShapes, *origorder, *sidecheck;
  //known vals are given for checking
  double angles[2] = {80,70}; //these have to be between 10 &90 for testing
  int i,j,k,ii;

  // create test screw mesh
  rc = fc_createDataset("screw dataset", &dataset);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create dataset");
  rc = _fc_makeMeshFromFiles(dataset, "../data/screw.vert",
			     "../data/screw.elem", &mesh);
  fail_unless(rc == FC_SUCCESS, "abort: failed to create mesh");


  // ---  screw shape --- //
  rc = fc_getMeshShapes(mesh,80,1,&numShapes,&shapes);
  fail_unless(rc == FC_SUCCESS, "abort: failed to get shapes"); 
  fail_unless(numShapes == 1, "wrong number of shapes");
  fail_unless(shapes[0].numSides == 5, "wrong number of sides");

  rc = fc_createScrewShapeOrder(&shapes[0],&origorder);
  fail_unless(rc == FC_SUCCESS, "failed to get screw order");

  for (ii = 0; ii < 2; ii++){
    int *temporder;

    rc = fc_reshapeByAngle(&shapes[0],angles[ii],1,&temp_shape);
    fail_unless(rc == FC_SUCCESS, "failed to reshape by angle");
    fail_unless(temp_shape.numSides == shapes[0].numSides, 
		"reshape should return the same number of sides");    

    //check sides here but they might not be in the same order
    rc = fc_createScrewShapeOrder(&temp_shape,&temporder);
    fail_unless(rc == FC_SUCCESS, "failed to get screw order");
    
    for (i = 0; i < shapes[0].numSides; i++){
      for (j = 0; j < 2; j++){
	FC_Subset origsubset, tempsubset;
	int numorigmember, numtempmember;
	
	if (j == 0){
	  origsubset = shapes[0].faces[origorder[j]];
	  tempsubset = temp_shape.faces[temporder[j]];
	}else{
	  origsubset = shapes[0].elems[origorder[j]];
	  tempsubset = temp_shape.elems[temporder[j]];
	}
	rc = fc_getSubsetNumMember(origsubset, &numorigmember);
	fail_unless(rc == FC_SUCCESS, "cant get subset members");
	rc = fc_getSubsetNumMember(tempsubset, &numtempmember);
	fail_unless(rc == FC_SUCCESS, "cant get subset members");
	fail_unless(numorigmember == numtempmember,
		    "wrong number of members");
	//since already checked equal members, just have to check
	//superset in either order
	rc = fc_isSubsetSuperset(tempsubset,origsubset,&k);
	fail_unless(rc == FC_SUCCESS, "cant determine superset");
	fail_unless(k == 1, "wrong members");
      }
	
      for (j = 0; j < shapes[0].numSides; j++){
	fail_unless(shapes[0].adjmatrix[origorder[i]][origorder[j]] ==
		    temp_shape.adjmatrix[temporder[i]][temporder[j]],
		    "wrong adj matrix");
      }
    }

    free(temporder);
    fc_freeShape(&temp_shape);
  }
  free(origorder);

  // ---  becomes 1 side --- //
  rc = fc_getMeshShapes(mesh,90,1,&numTestShapes,&testshapes);
  fail_unless(rc == FC_SUCCESS, "abort: failed to get shapes"); 
  fail_unless(numTestShapes == 1, "wrong number of shapes");
  fail_unless(testshapes[0].numSides == 1, "wrong number of sides");  

  rc = fc_reshapeByAngle(&shapes[0],90,1,&temp_shape);
  fail_unless(rc == FC_SUCCESS, "failed to reshape by angle");
  fail_unless(temp_shape.numSides == 1,
		"reshape should return the same number of sides");    

  rc = fc_isSubsetSuperset(temp_shape.faces[0],testshapes[0].faces[0],&k);
  fail_unless(rc == FC_SUCCESS, "cant determine superset");
  fail_unless(k == 1, "wrong members");

  rc = fc_isSubsetSuperset(testshapes[0].faces[0],temp_shape.faces[0],&k);
  fail_unless(rc == FC_SUCCESS, "cant determine superset");
  fail_unless(k == 1, "wrong members");

  fc_freeShape(&testshapes[0]);
  free(testshapes);
  fc_freeShape(&temp_shape);

  // --- more than 5 sides -- //
  rc = fc_getMeshShapes(mesh,9,1,&numTestShapes,&testshapes);
  fail_unless(rc == FC_SUCCESS, "abort: failed to get shapes"); 
  fail_unless(numTestShapes == 1, "wrong number of shapes");
  fail_unless(testshapes[0].numSides == 42, "wrong number of sides");  

  rc = fc_reshapeByAngle(&shapes[0],9,1,&temp_shape);
  fail_unless(rc == FC_SUCCESS, "failed to reshape by angle");
  fail_unless(temp_shape.numSides == testshapes[0].numSides,
		"reshape should return the same number of sides");    

  sidecheck = (int*)malloc(testshapes[0].numSides*sizeof(int));
  for (i = 0; i < testshapes[0].numSides; i++){
    sidecheck[i] = 0;
  }

  rc = fc_isSubsetSuperset(testshapes[0].faces[39],
			   temp_shape.faces[1], &k);  
  fail_unless(rc == FC_SUCCESS, "cant determine superset");
  for (i = 0; i < testshapes[0].numSides; i++){
    for (j = 0; j < testshapes[0].numSides; j++){
      k = 0;
      rc = fc_isSubsetSuperset(testshapes[0].faces[i],
			       temp_shape.faces[j], &k);
      fail_unless(rc == FC_SUCCESS, "cant determine superset");
      if (k){
	rc = fc_isSubsetSuperset(temp_shape.faces[j], 
				 testshapes[0].faces[i],&k);
	fail_unless(rc == FC_SUCCESS, "cant determine superset");
	fail_unless(k, "subsets with the same face dont match");
	rc = fc_isSubsetSuperset(temp_shape.elems[j], 
				 testshapes[0].elems[i],&k);
	fail_unless(rc == FC_SUCCESS, "cant determine superset");
	fail_unless(k, "elem subsets with the same face dont match");
	rc = fc_isSubsetSuperset(testshapes[0].elems[i],
				 temp_shape.elems[j], &k);
	fail_unless(rc == FC_SUCCESS, "cant determine superset");
	fail_unless(k, "elem subsets with the same face dont match");
	sidecheck[j] = 1;
	break;
      }
    }
  }
  for (i = 0; i < testshapes[0].numSides; i++){
    fail_unless(sidecheck[i] == 1, "missing side in the check");
  }

  fc_freeShape(&temp_shape);
  fc_freeShape(&testshapes[0]);
  free(testshapes);
  free(sidecheck);

  //bad args for reshape by angle
  temp_shape.numSides =3;
  rc = fc_reshapeByAngle(NULL,30,2,&temp_shape);
  fail_unless(rc != FC_SUCCESS, "should fail for NULL shape");
  fail_unless(temp_shape.numSides == 0, "shape should return as init for fail");
  fail_unless(temp_shape.faces == NULL, "shape should return as init for fail");
  fail_unless(temp_shape.elems == NULL, "shape should return as init for fail"); 
  fail_unless(temp_shape.adjmatrix == NULL, "shape should return as init for fail"); 

  temp_shape.numSides =3;
  rc = fc_reshapeByAngle(&shapes[0],-1,2,&temp_shape);
  fail_unless(rc != FC_SUCCESS, "should fail for bad angle");
  fail_unless(temp_shape.numSides == 0, "shape should return as init for fail");
  fail_unless(temp_shape.faces == NULL, "shape should return as init for fail");
  fail_unless(temp_shape.elems == NULL, "shape should return as init for fail"); 
  fail_unless(temp_shape.adjmatrix == NULL, "shape should return as init for fail"); 

  temp_shape.numSides =3;
  rc = fc_reshapeByAngle(&shapes[0],181,2,&temp_shape);
  fail_unless(rc != FC_SUCCESS, "should fail for bad angle");
  fail_unless(temp_shape.numSides == 0, "shape should return as init for fail");
  fail_unless(temp_shape.faces == NULL, "shape should return as init for fail");
  fail_unless(temp_shape.elems == NULL, "shape should return as init for fail"); 
  fail_unless(temp_shape.adjmatrix == NULL, "shape should return as init for fail"); 

  temp_shape.numSides =3;
  rc = fc_reshapeByAngle(&shapes[0],30,-1,&temp_shape);
  fail_unless(rc != FC_SUCCESS, "should fail for bad segdim");
  fail_unless(temp_shape.numSides == 0, "shape should return as init for fail");
  fail_unless(temp_shape.faces == NULL, "shape should return as init for fail");
  fail_unless(temp_shape.elems == NULL, "shape should return as init for fail"); 
  fail_unless(temp_shape.adjmatrix == NULL, "shape should return as init for fail"); 

  temp_shape.numSides =3;
  rc = fc_reshapeByAngle(&shapes[0],30,3,&temp_shape);
  fail_unless(rc != FC_SUCCESS, "should fail for bad segdim");
  fail_unless(temp_shape.numSides == 0, "shape should return as init for fail");
  fail_unless(temp_shape.faces == NULL, "shape should return as init for fail");
  fail_unless(temp_shape.elems == NULL, "shape should return as init for fail"); 
  fail_unless(temp_shape.adjmatrix == NULL, "shape should return as init for fail"); 

  rc = fc_reshapeByAngle(&shapes[0],30,1,NULL);
  fail_unless(rc != FC_SUCCESS, "should fail for NULL arg");

  for (i = 0; i < numShapes; i++){
    fc_freeShape(&shapes[i]);
  }
  free(shapes);
}
END_TEST



// Populate the Suite with the tests

Suite *shape_suite(void)
{
  //see notes above re assumptions of what works, which
  //means order is important here.

  Suite *suite = suite_create("Shape");
  TCase *tc_shapesides_bad = tcase_create(" - ShapeSides Bad ");  
  TCase *tc_shapematrix = tcase_create(" - ShapeMatrix ");
  TCase *tc_basic = tcase_create(" - Basic ");
  TCase *tc_normal = tcase_create(" - SideNormal ");
  TCase *tc_reshape_all = tcase_create(" - Reshape_all ");


  suite_add_tcase(suite, tc_shapesides_bad);
  tcase_add_checked_fixture(tc_shapesides_bad, shape_setup, shape_teardown);
  tcase_add_test(tc_shapesides_bad, shapesides_bad);

  suite_add_tcase(suite, tc_shapematrix);
  tcase_add_checked_fixture(tc_shapematrix, shape_setup, shape_teardown);
  tcase_add_test(tc_shapematrix, shapematrix);


  //normal has to be tested before order (meanstdev is used in opposing sides)
  suite_add_tcase(suite, tc_normal);
  tcase_add_checked_fixture(tc_normal, shape_setup, shape_teardown);
  tcase_add_test(tc_normal, normaltest);

  suite_add_tcase(suite, tc_basic);
  tcase_add_checked_fixture(tc_basic, shape_setup, shape_teardown);
  tcase_add_test(tc_basic, helpers);
  tcase_add_test(tc_basic, orders);
  tcase_add_test(tc_basic, screw);
  tcase_add_test(tc_basic, ends);

  //reshape has to be tested after order, uses order
  //for checking results
  suite_add_tcase(suite, tc_reshape_all);
  tcase_add_checked_fixture(tc_reshape_all, shape_setup, shape_teardown);
  tcase_add_test(tc_reshape_all, reduce);
  tcase_add_test(tc_reshape_all, reshape);

  return suite;
}
