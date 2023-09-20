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
 * \file geom.c
 * \brief Implementation of \ref GeometricRelations module.
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/modules/geom.c,v $
 * $Revision: 1.193 $ 
 * $Date: 2007/07/17 06:08:55 $
 */

// C library includes
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

// fc library dependencies
#include "base.h"
#include "storage.h"
#include "library.h"
#include "dataset.h"
#include "mesh.h"
#include "variable.h"
#include "subset.h"
#include "topo.h"
#include "util.h"
#include "meshP.h"
#include "subsetP.h"
#include "topoP.h"

// this module
#include "geom.h"
#include "geomP.h"

/**
 * \addtogroup GeometricRelations Geometric Relations
 * 
 * \brief  Relationships and values which depend on the coordinates 
 *         in space. 
 *
 * \description
 *
 *    The geometric functions compute relationships and values which depend
 *    on the coordinates of mesh vertices.  For example, computing
 *    the distances between vertices, the element and vertex normals of 
 *    surfaces, etc.
 */

/**
 *\ingroup GeometricRelations
 *\defgroup PrivateGeometricRelations (Private)
 */

/** \name General geometry */
//-------------------------------
//@{

/**
 * \ingroup  GeometricRelations
 * \brief  Compute the distance between two vertices
 *
 * \description 
 *
 *    Computes the euclidean distance between the two sets of vertex 
 *    coordinates
 *
 *    It's cheaper to compute the squared euclidean distance than 
 *    the euclidean distance. If you want to compare distances, it's
 *    faster to compute and compare the squared euclidean distance
 *    (see \ref fc_calcSquaredEuclideanDistance()).
 *
 * \modifications 
 *   - 07-AUG-2002 Natalie Rooney Created.
 *   - 28-AUG-2002 Natalie Rooney Changing inputs to be Vectors instead 
 *   of float arrays.
 *   - 2003-AUG-05  W Koegler  Moved from displacement.c to here,
 *   made internal calculations use doubles, minor name change.
 *   - 6/14/2006 WSD Guts moved to fc_calcSquaredEuclideanDistance().
 */
FC_ReturnCode fc_calcEuclideanDistance( 
  FC_Coords v1Coords,   /**< input - coordinates of one vertex */
  FC_Coords v2Coords,   /**< input - coordinates of the other vertex */
  int dim,              /**< input - number of dimensions in the coords */
  double* distance      /**< output - distance between the two vertices */
) {
  FC_ReturnCode rc;
  double dist2;

  // default output
  if (distance)
    *distance = -1;

  // check input
  if (distance == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // log message
  fc_printfLogMessage("Calculating euclidean distance.");
  
  // do it
  rc = fc_calcSquaredEuclideanDistance(v1Coords, v2Coords, dim, &dist2);
  if (rc != FC_SUCCESS)
    return rc;
  *distance = sqrt(dist2);

  return FC_SUCCESS;
}

/**
 * \ingroup  GeometricRelations
 * \brief  Compute the squared distance between two vertices
 *
 * \description 
 *
 *    Computes the squared euclidean distance between the two sets of vertex
 *    coordinates.
 *
 *    It's cheaper to compute the squared euclidean distance than 
 *    the euclidean distance. If you want to compare distances, it's
 *    faster to compute and compare the squared euclidean distance.
 *
 * \modifications 
 *   - 6/14/2006 WSD Created. Contains old guts of fc_calcEuclideanDistance.
 */
FC_ReturnCode fc_calcSquaredEuclideanDistance( 
  FC_Coords v1Coords,   /**< input - coordinates of one vertex */
  FC_Coords v2Coords,   /**< input - coordinates of the other vertex */
  int dim,              /**< input - number of dimensions in the coords */
  double* distance2     /**< output - distance between the two vertices */
) {
  int i;
  double sumOfSquares, difference;

  // default output
  if (distance2)
    *distance2 = -1;
  
  // check input
  if (v1Coords == NULL || v2Coords == NULL || dim < 1 || dim > 3 ||
      distance2 == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // log message
  fc_printfLogMessage("Calculating squared euclidean distance.");

  sumOfSquares = 0.;
  for (i = 0; i < dim; i++ ) {
    difference = v1Coords[i] - v2Coords[i];
    sumOfSquares += difference * difference;
  }

  *distance2 = sumOfSquares;

  return FC_SUCCESS;
}


/** 
 * \ingroup  GeometricRelations
 * \brief  Determine angle between two vectors. 
 *
 * \description
 *  
 *    Determine angle between two vectors. Because of the return val
 *    of acos, vals for the angle are limited to 0->180 degrees.
 *
 * \modifications 
 *    - 01/16/06 ACG revised from old surface_grow when switched from
 *       old recursive algorithm to new WSK method for getSurfaceSides.
 *    - 01/16/06 ACG moved here from shape.
 *    - 06/01/06 ACG changing comparison from DBL_EQUIV to VALUE_EQUIV
 *                 with a greater epsilon becuase some cases that
 *                 should have been rounded off were not.
 */
FC_ReturnCode fc_calcAngleBetweenVectors(   
  FC_Vector vector1,        /**< input - 1st normal */
  FC_Vector vector2,        /**< input - 2nd normal */
  double *diff_angle          /**< output - angle for normal comparison */ 
) {

  double tempA, tempB, magA, magB, dotprod, diffrad;
  int j;


  //default return values
  if (diff_angle)
    *diff_angle = -1;

  if(vector1 == NULL || vector2 == NULL || !diff_angle){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  fc_printfLogMessage("Getting angle between two vectors");

  tempA = 0;
  tempB = 0;
  dotprod = 0;
  for (j = 0; j < 3; j++){
    tempA+=vector1[j]*vector1[j];
    tempB+=vector2[j]*vector2[j];
    dotprod+=vector1[j]*vector2[j];
  }
  magA = sqrt(tempA);
  magB = sqrt(tempB);

  if (FC_DBL_EQUIV(magA,0.0) || FC_DBL_EQUIV(magB,0.0) ||
      FC_DBL_EQUIV((magA*magB),0.0)){
    fc_printfErrorMessage("Zero magitude input vector");
    return FC_INPUT_ERROR;
  }

  dotprod/=(magA*magB);

  if (FC_VALUE_EQUIV(dotprod,-1,(1000*DBL_EPSILON),DBL_MIN)){ //round
    dotprod = -1;
  }else if (FC_VALUE_EQUIV(dotprod,1,(1000*DBL_EPSILON),DBL_MIN)){ //round
    dotprod = 1;
  }

  if (dotprod < -1 || dotprod > 1){
    fc_printfErrorMessage("Invalid val for dot prod arg to acos: %20g\n",dotprod);
    fc_printfErrorMessage("magA = %10g magB = %10g vector1 = (%10g,%10g,%10g) vector2 = (%10g,%10g,%10g) dotprod = %10g\n",magA, magB,vector1[0],vector1[1],vector1[2],vector2[0],vector2[1],vector2[2],dotprod);
    return FC_ERROR;
  }

  diffrad = acos(dotprod); //will go between 0 and pi rad
  *diff_angle = diffrad*180./FC_PI;

  return FC_SUCCESS; 
}


//@}

/** \name General geometry helpers */
//-------------------------------
//@{


//@}

/** \name Displaced coordinates */
//-------------------------------
//@{

/**
 * \ingroup GeometricRelations
 * \brief Test to see if variable can be used as displacment variable.
 *
 * \description
 *
 *   Test a variable to see if it can be used as a displacement variable on the
 *   given mesh. The variable does NOT have to be on the given mesh, but has to
 *   match enough to be used. The variable must have the same number of as the
 *   mesh has vertices and the number of components must equal the
 *   dimensionality of the mesh.
 *
 *   Returns true (1) if the variable can be used as displacements on the
 *   mesh, false (0) if the variable is not appropriate to use on the
 *   mesh (including the case if either the mesh or variable are not valid).
 *
 * \todo Should we be more strick. Should the var have to be on verts? 
 *       Should the var have the correct math type?
 *  
 * \modifications
 *    - 02/14/2006 WSD, Created.
 */
int fc_isValidDisplacementVariable(
  FC_Mesh mesh,         /**< input - the mesh that the displacement would be  
			   applied to */
  FC_Variable variable  /**< input - the variable to check */
) {
  int numDim, numVertex, numDataPoint, numComp;

  // test input
  if (!fc_isMeshValid(mesh) || !fc_isVariableValid(variable)) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return 0;
  }

  // collect info
  fc_getMeshInfo(mesh, NULL, &numDim, &numVertex, NULL, NULL);
  fc_getVariableInfo(variable, &numDataPoint, &numComp, NULL, NULL, NULL);

  // testing
  if (numDim != numComp || numVertex != numDataPoint)
    return 0;

  // passed tests
  return 1;
}

/**
 * \ingroup  GeometricRelations
 * \brief Create displaced vertex coordinates
 *
 * \description
 *  
 *   Create vertex coordinates which are the sum of the vertex coordinates
 *   and a displacement variable. Unlike the call 'fc_getMeshCoords', the 
 *   caller is responsible for freeing the created coordinates.
 *
 *   The variable to use as displacements must be associated with vertices
 *   and have the same dimensionality as the vertices.
 *
 * \modifications 
 *    - 2003-SEP-04  W Koegler  Moved code from getDisplacedSurfaceNormals here
 *    - 01/29/04 WSK, Changed to return double array instead of array of
 *      _FC_Vertex. Now can be a public function and more generally useful.
 */
FC_ReturnCode fc_getDisplacedMeshCoords(
  FC_Mesh mesh,               /**< input - mesh */
  FC_Variable displ,          /**< input - variable to displace coordinates by */
  double **displaced_coords   /**< output - array of coords in XYZXYZ... order */
) {
  FC_ReturnCode rc;
  int i;
  FC_DataType datatype;
  int numData, dataDim;
  double* coords;
  double* temp_displ_coords;
  void* data;
  
  // default return
  if (displaced_coords)
    *displaced_coords = NULL;

  // check input
  if (!fc_isMeshValid(mesh) || 
      !fc_isValidDisplacementVariable(mesh, displ) ||
      displaced_coords == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // log message
  fc_printfLogMessage("Computing displaced vertex coordinates.");

  // get info & big data
  rc = fc_getVariableInfo(displ, &numData, &dataDim, NULL, NULL, &datatype);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_getMeshCoordsPtr(mesh, &coords);
  if (rc != FC_SUCCESS)
    return rc;
  fc_getVariableDataPtr(displ, &data);
  if (rc != FC_SUCCESS)
    return rc;

  // make the displacement coords
  temp_displ_coords = (double*)malloc(sizeof(double)*numData*dataDim);
  if ( temp_displ_coords == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  } 
  switch (datatype) {
  case FC_DT_INT: {
    int *cast_data = (int*)data;
    for (i = 0; i < numData*dataDim; i++)
      temp_displ_coords[i] = coords[i] + cast_data[i];
  } break;
  case FC_DT_FLOAT: {
    float *cast_data = (float*)data;
    for (i = 0; i < numData*dataDim; i++)
      temp_displ_coords[i] = coords[i] + cast_data[i];
  } break;
  case FC_DT_DOUBLE: {
    double *cast_data = (double*)data;
    for (i = 0; i < numData*dataDim; i++)
      temp_displ_coords[i] = coords[i] + cast_data[i];
  } break;
  default:
    free(temp_displ_coords);
    return FC_ERROR;
  }  
  *displaced_coords = temp_displ_coords;

  return FC_SUCCESS;
}

//@}

/** \name Bounding boxes */
//-------------------------------
//@{

/**
 * \ingroup  GeometricRelations
 * \brief Determine if a bounding box is valid.
 *
 * \description
 *
 *    For a bounding box to be valid, the lowpoint cannot be "higher" in
 *    coordinates space than the highpoint. This is determined by a
 *    component-wise comparison of coordinates. A zero width bounding
 *    box (lowpoint = highpoint) is considered valid.
 *
 *    Returns 1 (true) for valid, 0 (false) for not valid.
 *
 * \modifications   
 *   - 10/28/04 WSK, created.
 */
int fc_isBoundingBoxValid(
  int numDim,          /**< input - the number of dimensions of bbox */    
  FC_Coords lowpoint,  /**< input - lower left corner of bbox */
  FC_Coords highpoint  /**< input - upper right corner of bbox */
) {
  int i;

  if (numDim < 1 || numDim > 3 || !lowpoint || !highpoint)
    return 0; // not valid
  for (i = 0; i < numDim; i++) {
    if (fc_gtd(lowpoint[i], highpoint[i]))
      return 0; // not valid
  }
 
  // passed all tests
  return 1; // valid
}

/**
 * \ingroup  GeometricRelations
 * \brief Determine if two bounding boxes overlap in space
 *
 * \description
 *
 *    Returns a flag which will be 0 if the bounding boxes do not overlap
 *    and greater than 1 if they do. The flag also indicates the type of
 *    overlap where 1 = bb1 is inside of bb2, 2 = bb2 is inside of bb1,
 *    3 = bb1 is the same as bb2, and 4 = partial overlap.
 *
 *    Note that a bb can still be inside of another if they share a boundary.
 *    If the bounding boxes do intersect, the bounding box of the overlap can
 *    be returned.
 *
 * \modifications   
 *   - 10/27/04 WSK, created.
 */
FC_ReturnCode fc_getBoundingBoxesOverlap(
  int numDim,            /**< input - the number of dimensions of bbox */    
  FC_Coords lowpoint1,   /**< input - lower left corner of bbox 1 */
  FC_Coords highpoint1,  /**< input - upper right corner of bbox 1 */
  FC_Coords lowpoint2,   /**< input - lower left corner of bbox 2 */
  FC_Coords highpoint2,  /**< input - upper right corner of bbox 2 */
  int* overlap_flag,     /**< output - 0 = no overlap, 1 = they
                            do overlap, -1 = error */
  FC_Coords *overlap_lowpoint,    /**< output - (optional) lower left 
                                       corner of overlap bbox */
  FC_Coords *overlap_highpoint    /**< output - (optional) upper right
                                       corner of overlap bbox */
) {
  int i;
  int oneIsTwo, oneInsideTwo, twoInsideOne;
 
  // default returns
  if (overlap_flag)
    *overlap_flag = -1;
  if (overlap_lowpoint) 
    for (i = 0; i < 3; i++)
      (*overlap_lowpoint)[i] = -1.;
  if (overlap_highpoint) 
    for (i = 0; i < 3; i++)
      (*overlap_highpoint)[i] = -2.; // makes bb invalid

  // check input
  if (numDim < 1 || numDim > 3 || !lowpoint1 || !highpoint1 || !lowpoint2 || 
      !highpoint2 || !overlap_flag ||
      !fc_isBoundingBoxValid(numDim, lowpoint1, highpoint1) ||
      !fc_isBoundingBoxValid(numDim, lowpoint2, highpoint2) ) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Getting bounding box overlap");   

  // look for no overlap - if either high1 < low2 or low1 > high2 
  *overlap_flag = 1;
  for (i = 0; i < numDim; i++) {
    if (fc_ltd(highpoint1[i], lowpoint2[i]) ||
	fc_gtd(lowpoint1[i], highpoint2[i])) {
      *overlap_flag = 0; // has to overlap in all dimensions
      break;
    }
  }

  // if they don't overlap, we're done
  if (*overlap_flag == 0)
    return FC_SUCCESS;


  // figure out what kind of overlap
  oneIsTwo = 0;
  oneInsideTwo = 0;
  twoInsideOne = 0;
  for (i = 0; i < numDim; i++) {
    int lowDblEquiv = FC_DBL_EQUIV(lowpoint1[i], lowpoint2[i]);
    int highDblEquiv = FC_DBL_EQUIV(highpoint1[i], highpoint2[i]); 
    if (lowDblEquiv && highDblEquiv)
      oneIsTwo++;
    else if (( lowpoint1[i]  >= lowpoint2[i]  || lowDblEquiv ) && 
	     ( highpoint1[i] <= highpoint2[i] || highDblEquiv ))
      oneInsideTwo++;
    else if (( lowpoint2[i]  >= lowpoint1[i]  || lowDblEquiv ) &&
             ( highpoint2[i] <= highpoint1[i] || highDblEquiv ))
      twoInsideOne++;
  }
  
  // set overlap flag & set overlap bb if requested
  if (oneIsTwo == numDim || oneInsideTwo == numDim) {
    if (oneIsTwo == numDim)
      *overlap_flag = 3;
    else
      *overlap_flag = 1;
    if (overlap_lowpoint) 
      for (i = 0; i < numDim; i++) 
        (*overlap_lowpoint)[i] = lowpoint1[i];
    if (overlap_highpoint)
      for (i = 0; i < numDim; i++)
        (*overlap_highpoint)[i] = highpoint1[i];
  }
  else if (twoInsideOne == numDim) {
    *overlap_flag = 2;
    if (overlap_lowpoint) 
      for (i = 0; i < numDim; i++) 
        (*overlap_lowpoint)[i] = lowpoint2[i];
    if (overlap_highpoint)
      for (i = 0; i < numDim; i++)
        (*overlap_highpoint)[i] = highpoint2[i];
  }
  else {
    *overlap_flag = 4;
    if (overlap_lowpoint) {
      for (i = 0; i < numDim; i++) {
        if (lowpoint1[i] > lowpoint2[i])
          (*overlap_lowpoint)[i] = lowpoint1[i];
        else
          (*overlap_lowpoint)[i] = lowpoint2[i];
      }
    }
    if (overlap_highpoint) {
      for (i = 0; i < numDim; i++) {
        if (highpoint1[i] < highpoint2[i])
          (*overlap_highpoint)[i] = highpoint1[i];
        else
          (*overlap_highpoint)[i] = highpoint2[i];
      }
    }
  }

  // if numDim < 3, fill in rest of FC_Coords with zeros
  if (overlap_lowpoint) 
    for (i = numDim; i < 3; i++)
      (*overlap_lowpoint)[i] = 0.;
  if (overlap_highpoint)
    for (i = numDim; i < 3; i++)
      (*overlap_highpoint)[i] = 0.;
    
  return FC_SUCCESS;
}

/**
 * \ingroup  GeometricRelations
 * \brief Determine the bounding box of two bounding boxes.
 *
 * \description
 *
 *    Returns the bounding box that contains the given bounding boxes.
 *
 * \modifications   
 *   - 10/27/04 WSK, created.
 */
FC_ReturnCode fc_combineBoundingBoxes(
  int numDim,            /**< input - the number of dimensions of bbox */    
  FC_Coords lowpoint1,   /**< input - lower left corner of bbox 1 */
  FC_Coords highpoint1,  /**< input - upper right corner of bbox 1 */
  FC_Coords lowpoint2,   /**< input - lower left corner of bbox 2 */
  FC_Coords highpoint2,  /**< input - upper right corner of bbox 2 */
  FC_Coords *combine_lowpoint,    /**< output - (optional) lower left 
                                       corner of combined bbox */
  FC_Coords *combine_highpoint    /**< output - (optional) upper right
                                       corner of combined bbox */
) {
  int i;
 
  // default returns
  if (combine_lowpoint) 
    for (i = 0; i < 3; i++)
      (*combine_lowpoint)[i] = -1.;
  if (combine_highpoint) 
    for (i = 0; i < 3; i++)
      (*combine_highpoint)[i] = -2.; // makes bb invalid

  // check input
  if (!fc_isBoundingBoxValid(numDim, lowpoint1, highpoint1) ||
      !fc_isBoundingBoxValid(numDim, lowpoint2, highpoint2) ||
      !combine_lowpoint || !combine_highpoint) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Getting bounding box of bounding boxes");   

  // return the bb if asked for
  for (i = 0; i < numDim; i++) {
    if (lowpoint1[i] < lowpoint2[i])
      (*combine_lowpoint)[i] = lowpoint1[i];
    else
      (*combine_lowpoint)[i] = lowpoint2[i];
    if (highpoint1[i] > highpoint2[i])
      (*combine_highpoint)[i] = highpoint1[i];
    else
      (*combine_highpoint)[i] = highpoint2[i];
  }
  for (i = numDim; i < 3; i++) {
    (*combine_lowpoint)[i] = 0.;
    (*combine_highpoint)[i] = 0.;
  }
  
  return FC_SUCCESS;
}

/**
 * \ingroup  GeometricRelations
 * \brief Compute mesh bounding box
 *
 * \description
 *
 *    The smallest bounding box (sides parallel to the coordinate
 *    axes) that fits around the entire mesh. It is described by two points:
 *    lowpoint is the smallest coordinates of the box, and highpoint is
 *    the largest coordinates of the box.
 *
 *    Mesh entities may lie (some will by definition) on the boundary of 
 *    the bounding box.
 *
 * \modifications   
 *   - Nancy Collins  Created.
 *   - July 3, 2002  W Koegler  Adapted to modified _FC_Vertex\
 *   - July 9, 2002  W Koegler  Adapted to modified struct Vector
 *   - 2003-AUG-05  W Koegler  changed to use doubles internally
 *   - 04/19/04 WSK, Added return argument for number of dimensions
 */
FC_ReturnCode fc_getMeshBoundingBox(
  FC_Mesh mesh,            /**< input - mesh handle */
  int *numDim,             /**< output - the number of dimensions of bbox */
  FC_Coords *lowpoint,     /**< output - lower left corner of bbox */
  FC_Coords *highpoint     /**< output - upper right corner of bbox */
) {
  FC_ReturnCode rc;
  int i, j;
  int numVertex;
  double *coords_p, temp;
  // must initalize all entries of coords to 0 in case numDim < 3
  FC_Coords low = { 0., 0., 0. }, high = { 0., 0., 0. };

  // set default output
  if (numDim)
    *numDim = -1;
  if (lowpoint) 
    for (i = 0; i < 3; i++)
      (*lowpoint)[i] = -1.;
  if (highpoint) 
    for (i = 0; i < 3; i++)
      (*highpoint)[i] = -2.; // makes bb invalid

  // check input
  if (!fc_isMeshValid(mesh) || !numDim || !lowpoint || !highpoint) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // log message
  fc_printfLogMessage("Computing mesh bounding box.");  

  // setup
  rc = fc_getMeshInfo(mesh, NULL, numDim, &numVertex, NULL, NULL);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_getMeshCoordsPtr(mesh, &coords_p);
  if (rc != FC_SUCCESS)
    return rc;
  
  // initialize  - start with 1st point
  for (i = 0; i < *numDim; i++) 
    low[i] = high[i] = coords_p[i];
  
  // find min & max in each dim
  for (i = 1; i < numVertex; i++) {
    for (j = 0; j < *numDim; j++) {
      temp = coords_p[i*(*numDim)+j];
      if (temp < low[j])
        low[j] = temp;
      if (temp > high[j])
        high[j] = temp;
    }
  }

  for (i = 0; i < *numDim; i++) {
    (*lowpoint)[i] = low[i];
    (*highpoint)[i] = high[i];
  }

  return FC_SUCCESS;
}

/**
 * \ingroup  GeometricRelations
 * \brief Compute displaced mesh bounding box
 *
 * \description
 *
 *    See fc_getMeshBoundingBox() for a description of bounding boxes.
 *
 *    Gets the bounding box of the mesh after displacing the coordinates
 *    with the given variable.
 *
 * \modifications   
 *   - 10/27/04 WSK, Created.
 */
FC_ReturnCode fc_getDisplacedMeshBoundingBox(
  FC_Mesh mesh,            /**< input - mesh handle */
  FC_Variable coord_displ, /**< input - coord displacements */
  int *numDim,             /**< output - the number of dimensions of bbox */
  FC_Coords *lowpoint,     /**< output - lower left corner of bbox */
  FC_Coords *highpoint     /**< output - upper right corner of bbox */
) {
  FC_ReturnCode rc;
  int i, j;
  int numVertex;
  double *coords, temp;
  // must initalize all entries of coords to 0 in case numDim < 3
  FC_Coords low = { 0., 0., 0. }, high = { 0., 0., 0. };

  // set default output
  if (numDim)
    *numDim = -1;
  if (lowpoint) 
    for (i = 0; i < 3; i++)
      (*lowpoint)[i] = -1.;
  if (highpoint) 
    for (i = 0; i < 3; i++)
      (*highpoint)[i] = -2.; // makes bb invalid

  // check input
  if (!fc_isMeshValid(mesh) || 
      !fc_isValidDisplacementVariable(mesh, coord_displ) || 
      !numDim || !lowpoint || !highpoint) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // log message
  fc_printfLogMessage("Computing displaced mesh bounding box.");  

  // setup
  rc = fc_getMeshInfo(mesh, NULL, numDim, &numVertex, NULL, NULL);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_getDisplacedMeshCoords(mesh, coord_displ, &coords);
  if (rc != FC_SUCCESS)
    return rc;
  
  // initialize  - start with 1st point
  for (i = 0; i < *numDim; i++) 
    low[i] = high[i] = coords[i];
  
  // find min & max in each dim
  for (i = 1; i < numVertex; i++) {
    for (j = 0; j < *numDim; j++) {
      temp = coords[i*(*numDim)+j];
      if (temp < low[j])
        low[j] = temp;
      if (temp > high[j])
        high[j] = temp;
    }
  }
  free(coords);

  // set return values
  for (i = 0; i < *numDim; i++) {
    (*lowpoint)[i] = low[i];
    (*highpoint)[i] = high[i];
  }

  return FC_SUCCESS;
}

/**
 * \ingroup  GeometricRelations
 * \brief Compute bounding box of subset.
 *
 * \description
 *
 *    See fc_getMeshBoundingBox() for a description of bounding boxes.
 *
 *    Returns with an error if nothing is marked.
 *
 * \modifications   
 *   - 02/05/04 WSK, created.
 *   - 04/19/04 WSK, Added return argument for number of dimensions
 *   - 06/15/04 RM, marked list changed to subset
 *                  added return FC_ERROR in switch statement if wrong assoc 
 */
FC_ReturnCode fc_getSubsetBoundingBox(
 FC_Subset subset,        /**< input - subset handle   */  
 int *numDim,             /**< output - the number of dimensions of bbox */    
 FC_Coords *lowpoint,     /**< output - lower left corner of bbox */
 FC_Coords *highpoint     /**< output - upper right corner of bbox */
)
{
  FC_ReturnCode rc;
  int i, j;
  int numVertex, *vertexIDs, numMember, *memberIDs;
  double *coords_p, temp;
  // must initialize all entries of coords to 0 in case numDim < 3
  FC_Coords low = { 0., 0., 0. }, high = { 0., 0., 0. };
  FC_Mesh mesh;
  FC_AssociationType assoc;

  // default output
  if (numDim)
    *numDim = -1;
  if (lowpoint) 
    for (i = 0; i < 3; i++)
      (*lowpoint)[i] = -1.;
  if (highpoint) 
    for (i = 0; i < 3; i++)
      (*highpoint)[i] = -2.; // makes bb invalid

  // check input 
  if (!fc_isSubsetValid(subset) || !numDim || !lowpoint || !highpoint) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }  
  
  // log message
  fc_printfLogMessage("Computing subset bounding box");    
  
  // setup I
  rc = fc_getMeshFromSubset(subset, &mesh);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_getSubsetInfo(subset, &numMember, NULL, &assoc);
  if (rc != FC_SUCCESS)
    return rc;

  // Error case: if subset has zero members - return error
  if (numMember == 0 ) {
    fc_printfErrorMessage("Subset has 0 members");
    return FC_ERROR;
  }
  
  // special case -- whole association
  if (assoc == FC_AT_WHOLE_MESH) 
    return fc_getMeshBoundingBox(mesh, numDim, lowpoint, highpoint);

  // get vertices of the members
  rc = fc_getSubsetMembersAsArray(subset, &numMember, &memberIDs);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_changeMeshEntityType(mesh, assoc, numMember, memberIDs, 
                               FC_AT_VERTEX, 1, &numVertex, &vertexIDs);
  free(memberIDs);

  // setup II
  rc = fc_getMeshDim(mesh, numDim);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_getMeshCoordsPtr(mesh, &coords_p);
  if (rc != FC_SUCCESS)
    return rc;

  // initialize - start with 1st point
  for (i = 0; i < *numDim; i++)
    low[i] = high[i] = coords_p[vertexIDs[0]*(*numDim)+i];

  // find min & max in each dim
  for (i = 1; i < numVertex; i++) {
    for (j = 0; j < *numDim; j++) {
      temp = coords_p[vertexIDs[i]*(*numDim)+j];
      if (temp < low[j])
        low[j] = temp;
      if (temp > high[j])
        high[j] = temp;
    }
  }
  free(vertexIDs);

  for (i = 0; i < *numDim; i++) {
    (*lowpoint)[i] = low[i];
    (*highpoint)[i] = high[i];
  }

  return FC_SUCCESS;
}

/**
 * \ingroup  GeometricRelations
 * \brief Compute bounding box of subset of a displaced mesh.
 *
 * \description
 *
 *    See fc_getMeshBoundingBox() for a description of bounding boxes.
 *
 *    Gets the bounding box of the subset after displacing the coordinates of
 *    the parent mesh with the given variable.
 *
 *    Returns with an error if nothing is marked.
 *
 * \modifications   
 *   - 10/27/04 WSK, Created.
 */
FC_ReturnCode fc_getDisplacedSubsetBoundingBox(
  FC_Subset subset,        /**< input - subset handle   */  
  FC_Variable coord_displ, /**< input - coord displacements */
  int *numDim,             /**< output - the number of dimensions of bbox */    
  FC_Coords *lowpoint,     /**< output - lower left corner of bbox */
  FC_Coords *highpoint     /**< output - upper right corner of bbox */
)
{
  FC_ReturnCode rc;
  FC_Mesh mesh;
  int i, j;
  int numVertex, *vertexIDs, numMember, *memberIDs;
  double *coords, temp;
  // must initalize all entries of coords to 0 in case numDim < 3
  FC_Coords low = { 0., 0., 0. }, high = { 0., 0., 0. };
  FC_AssociationType assoc;

  // default output
  if (numDim)
    *numDim = -1;
  if (lowpoint) 
    for (i = 0; i < 3; i++)
      (*lowpoint)[i] = -1.;
  if (highpoint) 
    for (i = 0; i < 3; i++)
      (*highpoint)[i] = -2.; // makes bb invalid

  // check input 
  if (!fc_isSubsetValid(subset) || !numDim || !lowpoint || !highpoint) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }  
  rc = fc_getMeshFromSubset(subset, &mesh);
  if (rc != FC_SUCCESS)
    return rc;
  if (!fc_isValidDisplacementVariable(mesh, coord_displ)) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }  

  // log message
  fc_printfLogMessage("Computing displaced subset bounding box");    
  
  // setup I
  rc = fc_getSubsetInfo(subset, &numMember, NULL, &assoc);
  if (rc != FC_SUCCESS)
    return rc;

  // Error case: if subset has zero members - return error
  if (numMember == 0 ) {
    fc_printfErrorMessage("Subset has 0 members");
    return FC_ERROR;
  }
  
  // --- special case -- whole association
  if (assoc == FC_AT_WHOLE_MESH) 
    return fc_getDisplacedMeshBoundingBox(mesh, coord_displ, numDim, lowpoint, 
                                          highpoint);

  // get vertices of the members
  rc = fc_getSubsetMembersAsArray(subset, &numMember, &memberIDs);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_changeMeshEntityType(mesh, assoc, numMember, memberIDs, 
                               FC_AT_VERTEX, 1, &numVertex, &vertexIDs);
  free(memberIDs);

  // setup II
  rc = fc_getMeshDim(mesh, numDim);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_getDisplacedMeshCoords(mesh, coord_displ, &coords);
  if (rc != FC_SUCCESS)
    return rc;

  // initialize - start with 1st point
  for (i = 0; i < *numDim; i++)
    low[i] = high[i] = coords[vertexIDs[0]*(*numDim)+i];

  // find min & max in each dim
  for (i = 1; i < numVertex; i++) {
    for (j = 0; j < *numDim; j++) {
      temp = coords[vertexIDs[i]*(*numDim)+j];
      if (temp < low[j])
        low[j] = temp;
      if (temp > high[j])
        high[j] = temp;
    }
  }
  free(vertexIDs);
  free(coords);

  // set return values
  for (i = 0; i < *numDim; i++) {
    (*lowpoint)[i] = low[i];
    (*highpoint)[i] = high[i];
  }

  return FC_SUCCESS;
}

/**
 * \ingroup  GeometricRelations
 * \brief Create a mesh that represents a bounding box.
 *
 * \description
 *
 *    Returns a mesh, made of line elements, that represents
 *    a bounding box. This is useful for visualizing bounding boxes.
 *
 * \modifications   
 *   - 02/19/04 WSK, created.
 */
FC_ReturnCode fc_createBoundingBoxMesh(
  FC_Dataset dest_ds,     /* input - the dataset to add mesh to */
  char* bb_name,          /* input - the name of the new mesh */
  int dim,                /* input - the dimensionality of the bounding box */
  FC_Coords lowerCoords,  /* input - the lower coords of the bounding box */
  FC_Coords upperCoords,  /* input - the upper coords of the bounding box */
  FC_Mesh* bb_mesh        /* output - the new mesh */
) {
  FC_ReturnCode rc;
  int i;
  int numVertex[3] = { 2, 4, 8 }; 
  int numElement[3] = { 1, 4, 12 };
  int conns[3][12*2] = { { 0, 1 },
			 { 0, 1,
			   1, 2,
			   2, 3,
			   3, 0 },
			 { 0, 1,
			   1, 2,
			   2, 3,
			   3, 0,
			   4, 5,
			   5, 6,
			   6, 7,
			   7, 4,
			   0, 4,
			   1, 5,
			   2, 6,
			   3, 7 } };
  double coords[8*3];

  // default return
  if (bb_mesh)
    *bb_mesh = FC_NULL_MESH;

  // check input
  if (!fc_isDatasetValid(dest_ds) || !bb_name || dim < 1 || dim > 3 ||
      !lowerCoords || !upperCoords || !bb_mesh) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // log message
  fc_printfLogMessage("Creating a mesh representing bounding box.");   

  rc = fc_createMesh(dest_ds, bb_name, bb_mesh);
  if (rc != FC_SUCCESS)
    return rc;

  if (dim == 1) {
    coords[0] = lowerCoords[0];
    coords[1] = upperCoords[0];
  }
  else if (dim == 2) {
    // x coords
    coords[0] =   lowerCoords[0];
    coords[1*2] = upperCoords[0];
    coords[2*2] = upperCoords[0];
    coords[3*2] = lowerCoords[0];
    // y coords
    coords[1] =     lowerCoords[1];
    coords[1*2+1] = lowerCoords[1];
    coords[2*2+1] = upperCoords[1];
    coords[3*2+1] = upperCoords[1];
  }
  else if (dim == 3) {
    // x coords
    coords[0] =   lowerCoords[0];
    coords[1*3] = upperCoords[0];
    coords[2*3] = upperCoords[0];
    coords[3*3] = lowerCoords[0];
    coords[4*3] = lowerCoords[0];
    coords[5*3] = upperCoords[0];
    coords[6*3] = upperCoords[0];
    coords[7*3] = lowerCoords[0];
    // y coords
    coords[1] =     lowerCoords[1];
    coords[1*3+1] = lowerCoords[1];
    coords[2*3+1] = upperCoords[1];
    coords[3*3+1] = upperCoords[1];
    coords[4*3+1] = lowerCoords[1];
    coords[5*3+1] = lowerCoords[1];
    coords[6*3+1] = upperCoords[1];
    coords[7*3+1] = upperCoords[1];
    // z coords
    for (i = 0; i < 4; i++)
      coords[i*3+2] = lowerCoords[2];
    for (i = 4; i < 8; i++)
      coords[i*3+2] = upperCoords[2];
  }

  rc = fc_setMeshCoords(*bb_mesh, dim, numVertex[dim-1], coords);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_setMeshElementConns(*bb_mesh, FC_ET_LINE, numElement[dim-1], 
			      conns[dim-1]);
  if (rc != FC_SUCCESS)
    return rc;
  
  return FC_SUCCESS;
}

//@}

/** \name Projected bounding boxes */
//-------------------------------
//@{

/**
 * \ingroup GeometricRelations
 * \brief Returns indication if the bounding boxes of 2 elements, when projected
 *  onto the same planar surface (given a projection direction) haven
 *  any intersection.
 *
 * \description
 * 
 *  Returns indication if the bounding boxes of 2 elements, when projected
 *  onto the same planar surface (given a projection direction) have
 *  any intersection.
 *
 *  The intent of this method is to find overlap of projections of elements,
 *  but the mathematics of that can be very complex for general cases. Instead
 *  we are going with the more trivial bounding box approach and hoping
 *  that gives us enough information.
 *
 *  Returns 1 if there is any overlap and 0 if there is none. Single point
 *  or single edge overlap is considered not to be overlap, since that
 *  is only boundary overlap, except for cases where an element projects into
 *  a line or point. Returns -1 in case of error.
 *
 *  Note: I have the proj normal as FC_Coords rather than a vector. I have
 *  debated chanign this to a FC_Vector instead.
 *
 * \todo
 *   - see if the BB version gives the resolution necessary
 *   - see if the exceptions for the cases of points and lines
 *     are really the exceptions we want to make.
 *   - check results for 2 dimensions -- SHOULD THIS EVEN WORK FOR 2D elems ?
 *   - can reduce the allocation of intermediate results if turns
 *     out we wont actaully need to keep them (some intermediates
 *     im keeping while i decide on the algorithm)
 *
 *
 * \modifications 
 *   - 10/21/05 ACG created (but not put into production)
 *   - 10/24/05 ACG changed from face occlusion to projected element
 *     bounding box overlap to simplify algorithmic necessity of handling
 *     possibly concave simple polygon intersection.
 *   - 11/21/05 ACG putting into production now that lower level
 *      type changes are completed
 *   - 03/15/06 ACG changed projnormal type from double* to FC_Coords
 *   - 03/28/06 ACG fixed bug in proj matrix - i had it antisymmetirc
 *                  rather than symmetric
 */
FC_ReturnCode fc_projectedElementBoundingBoxIntersection(
  FC_Mesh mesh,          /**< input - mesh */
  int elemID_1,       /**< input -  elementID 1 */
  int elemID_2,      /**< input -  elementID 2 */  
  FC_Coords proj_normal,  /**< input - the proj normal */
  int* intersectval   /**< output - result 1 or 0 */
){
  FC_ReturnCode rc;
  int *elem_id_array;
  FC_ElementType elemtype;
  FC_Coords *refcoords, *projrefcoords, *compcoords, *projcompcoords; 
  double *mesh_coords;
  double minmax[2][6];
  double projmatrix[3][3];
  double tempval, val;

  //  FC_Vector temp_normal;
  int topodim,numDim;
  int line[2][3];
  int numelem, numvertex,e_vertex;
  int i,j,k;

  if (intersectval)
    *intersectval = -1;

  if (!fc_isMeshValid(mesh) || !intersectval){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  

  // log message
  fc_printfLogMessage("calc elem bounding box intersection.");   

  rc = fc_getMeshInfo (mesh,&topodim,&numDim,&numvertex,&numelem,&elemtype);
  if (rc != FC_SUCCESS)
    return rc;
 //Note: verify that ID is same as array index
  if (elemID_1 < 0 || elemID_1 > numelem ||
      elemID_2 < 0 || elemID_2 > numelem){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  tempval = 0.;
  for (i = 0; i < 3; i++){
    tempval+=proj_normal[i]*proj_normal[i];
  }

  if (tempval < 0){
    fc_printfErrorMessage("bad val for normal");
    return FC_INPUT_ERROR;
  }
  val = sqrt(fabs(tempval));

  //make sure normal has mag 1
  if (!(FC_FLT_EQUIV(val,1.0))){
    fc_printfErrorMessage("normal does not have magnitude 1");
    return FC_INPUT_ERROR;
  }

  //get numverticies for the elem
  e_vertex = fc_getElementTypeNumVertex(elemtype);
  if (e_vertex < 1){
    fc_printfErrorMessage("bad number of verticies for these elements");
    return FC_INPUT_ERROR;
  }

  rc = fc_getMeshElementConnsPtr(mesh,&elem_id_array);
  if (rc != FC_SUCCESS)
    return FC_ERROR;

  rc = fc_getMeshCoordsPtr(mesh,&mesh_coords);
  if (rc != FC_SUCCESS)
    return FC_ERROR;


  refcoords = (FC_Coords*)malloc(e_vertex*sizeof(FC_Coords));
  compcoords = (FC_Coords*)malloc(e_vertex*sizeof(FC_Coords));
  projrefcoords = (FC_Coords*)malloc(e_vertex*sizeof(FC_Coords));
  projcompcoords = (FC_Coords*)malloc(e_vertex*sizeof(FC_Coords));

  if (!refcoords || !compcoords || !projrefcoords || !projcompcoords){
    if (refcoords) free(refcoords);
    if (projrefcoords) free(projrefcoords);
    if (compcoords) free(compcoords);
    if (compcoords) free(projcompcoords);
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;			  
  }

  //  printf("\n");
  for (j = 0; j < e_vertex; j++){
    for (k = 0; k < numDim; k++)
      //       refcoords[j][k] = mesh_coords[conns[elemID_1*e_vertex+j]*numDim+k];
      refcoords[j][k] = mesh_coords[elem_id_array[elemID_1*e_vertex +j]*numDim+k];
    if (numDim < 3) refcoords[j][2] = 0;
  }

  for (j = 0; j < e_vertex; j++){
    for (k = 0; k < numDim; k++)
      compcoords[j][k] = mesh_coords[elem_id_array[elemID_2*e_vertex+j]*numDim+k];
    if (numDim < 3) compcoords[j][2] = 0;
  }


  // project each of the points onto the plane defined by the proj normal
  // its is ax+by+cz = 0 for the plane that intersects the origin where
  // a,b,and c are the components of the proj normal
  projmatrix[0][0] = proj_normal[1]*proj_normal[1]+
    proj_normal[2]*proj_normal[2];
  projmatrix[1][1] = proj_normal[0]*proj_normal[0]+
    proj_normal[2]*proj_normal[2];
  projmatrix[2][2] = proj_normal[0]*proj_normal[0]+
    proj_normal[1]*proj_normal[1];
  projmatrix[0][1] = -proj_normal[0]*proj_normal[1];
  projmatrix[0][2] = -proj_normal[0]*proj_normal[2];
  projmatrix[1][0] = projmatrix[0][1];
  projmatrix[2][0] = projmatrix[0][2];
  projmatrix[1][2] = -proj_normal[1]*proj_normal[2];
  projmatrix[2][1] = projmatrix[1][2];


  //proj each point of the elems onto the projection plane. dont
  //actaully need to keep these, but holding on to them
  //for now in case i change the algorithm
  //and get the bounding boxes for the resulting elements
  for (i = 0; i < 2; i++){
    for (j = 0; j <3; j++){
      minmax[i][2*j] = DBL_MAX;
      minmax[i][2*j+1] = -1*DBL_MAX;
    }
  }


  for (i = 0; i < e_vertex; i++){ //for each vertex
    for (j = 0; j < 3; j++){ //for each component
      projrefcoords[i][j] = 
	projmatrix[j][0] * refcoords[i][0]+
	projmatrix[j][1] * refcoords[i][1]+
	projmatrix[j][2] * refcoords[i][2];
      projcompcoords[i][j] = 
	projmatrix[j][0] * compcoords[i][0]+
	projmatrix[j][1] * compcoords[i][1]+
	projmatrix[j][2] * compcoords[i][2];
      minmax[0][2*j] = (projrefcoords[i][j] < minmax[0][2*j] ?
			projrefcoords[i][j] : minmax[0][2*j]);
      minmax[0][2*j+1] = (projrefcoords[i][j] > minmax[0][2*j+1] ?
			projrefcoords[i][j] : minmax[0][2*j+1]);
      minmax[1][2*j] = (projcompcoords[i][j] < minmax[1][2*j] ?
			projcompcoords[i][j] : minmax[1][2*j]);
      minmax[1][2*j+1] = (projcompcoords[i][j] > minmax[1][2*j+1] ?
			projcompcoords[i][j] : minmax[1][2*j+1]);
    }
  }


  //	can we capitalize on having projected them into the same plane ?
  //    (these are really flat elements now)
  //    rot the plane or something and then not have to consider "z" ???
  //    is it ok to check only one way ??

  //first check if it reduces to a line or point case
  for (i =0; i < 2; i++){
    for (j = 0; j < 3; j++){
      line[i][j] = (FC_DBL_EQUIV(minmax[i][2*j],minmax[i][2*j+1])? 1: 0);
    }
  }
  //if its a line in a direction, then we will consider it
  //as overlapping, if the edge aligns, otherwise need more
  //than just edge overlap
  k = 1;
  for (i = 0; i < 3; i++){
    if (line[0][i] || line [1][i]){ //allow edge overlap
      if ((minmax[0][2*i] > minmax[1][2*i+1]) ||
	  (minmax[0][2*i+1] < minmax[1][2*i])){
	k = 0;
	break;
      }
    }else {
      if (!(minmax[0][2*i] < minmax[1][2*i+1] &&
	  minmax[0][2*i+1] > minmax[1][2*i])){
	k = 0;
	break;
      }
    }
  }

  *intersectval = k;
  free(refcoords);
  free(projrefcoords);
  free(compcoords);
  free(projcompcoords);
  return FC_SUCCESS;
}

/**
 * \ingroup GeometricRelations
 * \brief Given a subset of elements ids and a reference subset of element ids,
 *  projects those elements onto the same plane and calculates which elements of
 *  the reference subset have bounding boxes in the projected plane that are
 *  intersected by the bounding boxes of the comparison elements in the 
 *  projected plane.
 *
 * \description
 * 
 *  Given a subset of elements ids and a reference subset of element ids,
 *  projects those elements onto the same plane and calculates which elements of
 *  the reference subset have bounding boxes in the projected plane that are
 *  intersected by the bounding boxes of the comparison elements in the 
 *  projected plane. 
 *
 *  This is intended for getting an estimate of amount of the tie surface
 *  that has dead elements above it, in some direction. Neither subset
 *  can be NULL. The reference subset cannot be empty. The comp
 *  subset can be empty.
 *
 *  Name of output subset is that of the input ref subset, postpended with
 *  "_Intersect"
 *
 *  Uses fc_projectedElementBoundingBoxIntersection.
 *
 * \modifications 
 *   - 10/25/05 ACG created (but not put into production)
 *   - 11/21/05 ACG putting into production now that lower level type changes are completed
 *   - 03/15/06 ACG changed projnormal type from double* to FC_Coords
 *   - 08/07/06 WSD changed orig linked lists to arrays, and the out linked
 *          list to a sorted int array.
 *   - 10/10/06 ACG changed so once an item satisfies overlap it is removed 
 *          from comparison
 *   - 10/16/06 removed mesh from arg list to make more consistent with other
 *          subset calls
 */
FC_ReturnCode fc_projectedSubsetElementBoundingBoxIntersection(
  FC_Subset compelemIDs, /**< input - subset of comp elem ids */
  FC_Subset refelemIDs, /**< input - subset of reference elem ids */
  FC_Coords proj_normal,  /**< input - proj normal */
  FC_Subset* intersectingrefIDs /**< output - subset of ref elem ids
				    satisfying the calculation */
){
  FC_ReturnCode rc;
  FC_Mesh mesh1, mesh2;
  int *orig1, *orig2;
  FC_SortedIntArray ref2, out2;
  FC_AssociationType assoc;
  int count1, count2;
  int i, j;

  if (intersectingrefIDs)
    *intersectingrefIDs = FC_NULL_SUBSET;

  if (!fc_isSubsetValid(compelemIDs) ||!fc_isSubsetValid(refelemIDs) ||
      !intersectingrefIDs){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  for (i = 0; i < 2; i++){
    FC_Subset temp;
    temp = i? compelemIDs: refelemIDs;
    rc = fc_getSubsetAssociationType(temp, &assoc);
    if (rc != FC_SUCCESS)
      return rc;
    
    if (assoc != FC_AT_ELEMENT){
      fc_printfErrorMessage("Wrong association type");
      return FC_INPUT_ERROR;
    }
  }
  //normal is checked in inner method since its just passed through

  rc = fc_getMeshFromSubset(compelemIDs,&mesh1);
  if (rc != FC_SUCCESS){
    fc_printfErrorMessage("Can't get mesh from subset");
    return rc;
  }

  rc = fc_getMeshFromSubset(refelemIDs,&mesh2);
  if (rc != FC_SUCCESS){
    fc_printfErrorMessage("Can't get mesh from subset");
    return rc;
  }

  if (!FC_HANDLE_EQUIV(mesh1,mesh2)){
    fc_printfErrorMessage("Subsets must be on same mesh");
    return FC_INPUT_ERROR;
  }

  rc = fc_getSubsetMembersAsArray(compelemIDs,&count1,&orig1);
  if (rc != FC_SUCCESS){
    fc_printfErrorMessage("Can't get subset memebers");
    return rc;
  }

  rc = fc_getSubsetMembersAsArray(refelemIDs,&count2,&orig2);
  if (rc != FC_SUCCESS){
    fc_printfErrorMessage("Can't get subset memebers");
    return rc;
  }

  if (count2 == 0){
    fc_printfErrorMessage("refsubset cannot be empty");
    free(orig1);
    free(orig2);
    return FC_INPUT_ERROR;
  }

  rc = fc_convertIntArrayToSortedIntArray(count2,orig2,1,&ref2);
  if (rc != FC_SUCCESS){
    fc_printfErrorMessage("Can't convert ref array");
    free(orig1);
    free(orig2);
  }

  //compare each item in the orig1 list to all (remaining) items in the ref2 list
  fc_initSortedIntArray(&out2);
  for (i = 0; i < count1 && ref2.numVal; i++) {
    int elem1 = orig1[i];
    FC_SortedIntArray rem; 
    fc_initSortedIntArray(&rem);
    for (j = 0; j < ref2.numVal; j++) {
      int elem2 = ref2.vals[j];
      int intersectval;
      //      printf("\t(%d) and elem %d\n",count2++,elem2);
      rc = fc_projectedElementBoundingBoxIntersection(mesh1,elem1,elem2,
						      proj_normal,&intersectval);
      if (rc != FC_SUCCESS){
	fc_printfErrorMessage("Can't determine intersection");
	free(orig1);
	fc_freeSortedIntArray(&out2);
	fc_freeSortedIntArray(&ref2);
	fc_freeSortedIntArray(&rem);
	return rc;
      }
      if (intersectval == 1){
	int ret = fc_addIntToSortedIntArray(&out2, elem2);
	if (ret < 0 ){
	  fc_printfErrorMessage("Can't add id to sorted int array");
	  free(orig1);
	  fc_freeSortedIntArray(&out2);
	  fc_freeSortedIntArray(&ref2);
	  fc_freeSortedIntArray(&rem);
	  return ret; //its an error code
	}
	ret = fc_addIntToSortedIntArray(&rem, elem2);
	if (ret < 0 ){
	  fc_printfErrorMessage("Can't add id to sorted int array");
	  free(orig1);
	  fc_freeSortedIntArray(&out2);
	  fc_freeSortedIntArray(&ref2);
	  fc_freeSortedIntArray(&rem);
	  return ret; //its an error code
	}
      }
    } //j
    for (j = 0; j < rem.numVal; j++){
      int ret = fc_deleteIntFromSortedIntArray(&ref2,rem.vals[j]);
      if (ret != 1){
	fc_printfErrorMessage("Can't delete id from sorted int array");
	free(orig1);
	fc_freeSortedIntArray(&out2);
	fc_freeSortedIntArray(&ref2);
	fc_freeSortedIntArray(&rem);
	return FC_ERROR;
      }
    }
    fc_freeSortedIntArray(&rem);
  } //i
  free(orig1);
  fc_freeSortedIntArray(&ref2);
  //rem is freed

  {// allowed to send back empty subsets
    _FC_SubSlot *subSlot;
    char *subsetName, *tempName;
    subSlot = _fc_getSubSlot(refelemIDs);
    subsetName = subSlot->header.name;
    tempName = malloc((strlen(subsetName)+20)*sizeof(char));
    if (tempName == NULL){
       fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
       fc_freeSortedIntArray(&out2);
       return FC_MEMORY_ERROR;
    }
    sprintf(tempName,"%s_Intersect",subsetName);
    rc = fc_createSubset(mesh1,tempName,FC_AT_ELEMENT,intersectingrefIDs);
    free(tempName);
    if (rc != FC_SUCCESS){
      fc_printfErrorMessage("cannot create new element subset");
      fc_freeSortedIntArray(&out2);
      return rc;
    }
    rc = fc_addArrayMembersToSubset(*intersectingrefIDs,out2.numVal,
				       out2.vals);
    fc_freeSortedIntArray(&out2);
    if (rc != FC_SUCCESS){
      fc_printfErrorMessage("cannot add items to subset");
      return rc;
    }
  }

  return FC_SUCCESS;
}

//@}

/** \name Diameter (max distance between any points in a point set) */
//-------------------------------
//@{

/**
 * \ingroup  GeometricRelations
 * \brief Get the diameter of a mesh.
 * 
 * The diameter of a set of points is the maximum distance between any two
 * points in the set. (It is also the diameter of the mimimum area circle
 * circumscribing the set of points). This routine returns the diameter
 * of the set of vertices that make up the mesh.
 *
 * If requested, this routine will also return the IDs of the vertices
 * that determined the diameter. There may be multiple point pairs
 * that have the same distance, so the number of pairs is returned
 * and an array which lists the first pair, then the next, etc.
 *
 * \modifications
 *    - 6/27/05 WSD Created.
 */
FC_ReturnCode fc_getMeshDiameter(
  FC_Mesh mesh,  /**< input - A mesh */
  double* diameter,  /**< ouput - diameter of the vertex set */
  int* numPair,      /**< output - number of pairs of vertices that are
		                   diameter distance from each other */
  int** pairIDs      /**< output - IDs of the vertex pairs */
) {
  FC_ReturnCode rc;
  FC_Subset subset;

  // log message
  fc_printfLogMessage("Computing mesh diameter.");

  // set default returns & test output args
  rc = _fc_getDiameterOutputArgsTest(diameter, numPair, pairIDs);
  if (rc != FC_SUCCESS)
    return rc;

  // test input
  if (!fc_isMeshValid(mesh)) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // Create subset that is whole
  rc = fc_createSubset(mesh, "temp", FC_AT_WHOLE_MESH, &subset);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_addMemberToSubset(subset, 0);
  if (rc != FC_SUCCESS)
    return rc;

  // Send to core getSubsetDiameter code
  rc = _fc_getSubsetDiameterCore(subset, NULL, diameter, numPair, pairIDs);

  // Cleanup & return
  fc_deleteSubset(subset);
  return rc;
}

/**
 * \ingroup  GeometricRelations
 * \brief Get the diameter of a displaced mesh.
 * 
 * This routine returns the diameter of the set of displaced vertices that make
 * up the mesh. For more details see fc_getMeshDiameter().
 *
 * \modifications
 *    - 6/27/05 WSD Created.
 */
FC_ReturnCode fc_getDisplacedMeshDiameter(
  FC_Mesh mesh,       /**< input - A mesh */
  FC_Variable coord_displ,  /**< input - the coordinate displacements */ 
  double* diameter,   /**< ouput - diameter of the vertex set */
  int* numPair,       /**< output - number of pairs of vertices that are
		                    diameter distance from each other */
  int** pairIDs       /**< output - IDs of the vertex pairs */
) {
  FC_ReturnCode rc;
  FC_Subset subset;

  // log message
  fc_printfLogMessage("Computing displaced mesh diameter.");

  // set default returns & test output args
  rc = _fc_getDiameterOutputArgsTest(diameter, numPair, pairIDs);
  if (rc != FC_SUCCESS)
    return rc;

  // test input
  if (!fc_isMeshValid(mesh) || 
      !fc_isValidDisplacementVariable(mesh, coord_displ)) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // Create subset that is whole
  rc = fc_createSubset(mesh, "temp", FC_AT_WHOLE_MESH, &subset);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_addMemberToSubset(subset, 0);
  if (rc != FC_SUCCESS)
    return rc;

  // Send to core getSubsetDiameter code
  rc = _fc_getSubsetDiameterCore(subset, &coord_displ, diameter, numPair, pairIDs);

  // Cleanup & return
  fc_deleteSubset(subset);
  return rc;
}

/**
 * \ingroup  GeometricRelations
 * \brief Get the diameter of a subset.
 * 
 * This routine returns the diameter of the set of vertices that make up the
 * members of the subset. For more details see fc_getMeshDiameter().
 *
 * The function will return with an error if the subset is empty.
 *
 * \modifications
 *    - 6/24/05 WSD Created.
 */
FC_ReturnCode fc_getSubsetDiameter(
  FC_Subset subset,  /**< input - A subset */
  double* diameter,  /**< ouput - diameter of the vertex set */
  int* numPair,      /**< output - number of pairs of vertices that are
		                   diameter distance from each other */
  int** pairIDs      /**< output - IDs of the vertex pairs */
) {
  FC_ReturnCode rc;
 
  // log message
  fc_printfLogMessage("Computing subset diameter.");

  // set default returns & test output args
  rc = _fc_getDiameterOutputArgsTest(diameter, numPair, pairIDs);
  if (rc != FC_SUCCESS)
    return rc;

  // test input
  if (!fc_isSubsetValid(subset)) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // do it
  return _fc_getSubsetDiameterCore(subset, NULL, diameter, numPair, pairIDs);
}

/**
 * \ingroup  GeometricRelations
 * \brief Get the diameter of a displaced subset.
 * 
 * This routine returns the diameter of the set of displaced vertices that make
 * up the members of the subset.  For more details see fc_getMeshDiameter().
 *
 * \modifications
 *    - 6/24/05 WSD Created.
 */
FC_ReturnCode fc_getDisplacedSubsetDiameter(
  FC_Subset subset,  /**< input - A subset */
  FC_Variable coord_displ,  /**< input - the owning mesh's coordinate 
			                 displacements */ 
  double* diameter,  /**< ouput - diameter of the vertex set */
  int* numPair,      /**< output - number of pairs of vertices that are
		                   diameter distance from each other */
  int** pairIDs      /**< output - IDs of the vertex pairs */
) {
  FC_ReturnCode rc;
  FC_Mesh mesh;
 
  // log message
  fc_printfLogMessage("Computing displaced subset diameter.");

  // set default returns & test output args
  rc = _fc_getDiameterOutputArgsTest(diameter, numPair, pairIDs);
  if (rc != FC_SUCCESS)
    return rc;

  // test input
  fc_getMeshFromSubset(subset, &mesh);
  if (!fc_isSubsetValid(subset) || 
      !fc_isValidDisplacementVariable(mesh, coord_displ)) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // do it
  return _fc_getSubsetDiameterCore(subset, &coord_displ, diameter, numPair, pairIDs);
}

/**
 * \ingroup  GeometricRelations
 * \brief Get the diameter of a group of meshes.
 * 
 * The diameter of a set of points is the maximum distance between any two
 * points in the set. (It is also the diameter of the mimimum area circle
 * circumscribing the set of points). This routine returns the diameter
 * of the set of vertices that make up a group of meshes.
 *
 * If requested, this routine will also return the IDs of the vertices
 * that determined the diameter. There may be multiple point pairs
 * that have the same distance, so the number of pairs is returned
 * and two arrays which identify the pairs. The pairIDs list the vertex
 * IDs (local to a mesh) for the first pair, then the next pair.
 * The pairMeshes is the originating mesh for each vertex in pairIDs.
 *
 * \modifications
 *    - 3/9/06 WSD Created.
 */
FC_ReturnCode fc_getMeshesDiameter(
  int numMesh, /**< input - number of meshes */
  FC_Mesh* meshes,  /**< input - Array of meshes */
  double* diameter,  /**< ouput - diameter of the vertex set */
  int* numPair,      /**< output - number of pairs of vertices that are
		                   diameter distance from each other */
  int** pairIDs,     /**< output - IDs of the vertex pairs */
  FC_Mesh** pairMeshes  /**< output - the meshes of the pairs */
) {
  FC_ReturnCode rc;
  int i;
  FC_Subset* subsets, *pairSubs;

  // log message
  fc_printfLogMessage("Computing meshes diameter.");

  // set default returns & test output args
  rc = _fc_getDiameterOutputArgsTestMeshes(diameter, numPair, pairIDs,
					   pairMeshes);
  if (rc != FC_SUCCESS)
    return rc;

  // test input
  if (numMesh < 1 || !meshes) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  for (i = 0; i < numMesh; i++) {
    if (!fc_isMeshValid(meshes[i])) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
      return FC_INPUT_ERROR;
    }
  }

  // Create subsets that are whole
  subsets = (FC_Subset*)malloc(numMesh*sizeof(FC_Subset));
  if (!subsets) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  for (i = 0; i < numMesh; i++) {
    rc = fc_createSubset(meshes[i], "temp", FC_AT_WHOLE_MESH, &subsets[i]);
    if (rc != FC_SUCCESS)
      return rc;
    rc = fc_addMemberToSubset(subsets[i], 0);
    if (rc != FC_SUCCESS)
      return rc;
  }

  // Send to core getSubsetsDiameter code
  rc = _fc_getSubsetsDiameterCore(numMesh, subsets, NULL, diameter, numPair, 
				  pairIDs, &pairSubs);

  // convert pairSubs to pairMeshes
  if (rc == FC_SUCCESS && numPair && *numPair > 0) {
    FC_Mesh *subMeshes;
    subMeshes = (FC_Mesh*)malloc(2*(*numPair)*sizeof(FC_Mesh));
    if (!meshes) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      return FC_MEMORY_ERROR;
    }
    for (i = 0; i < 2*(*numPair); i++) {
      rc = fc_getMeshFromSubset(pairSubs[i], &subMeshes[i]);
      if (rc != FC_SUCCESS)
      return rc;
    }
    *pairMeshes = subMeshes;
    free(pairSubs);
  }

  // Cleanup & return
  for (i = 0; i < numMesh; i++)
    fc_deleteSubset(subsets[i]);
  free(subsets);

  return rc;
}

/**
 * \ingroup  GeometricRelations
 * \brief Get the diameter of a group of displaced meshes.
 * 
 * This routine returns the diameter of the set of displaced vertices that make
 * up a group of meshes. For more details see fc_getMeshesDiameter().
 *
 * \modifications
 *    - 3/9/06 WSD Created.
 */
FC_ReturnCode fc_getDisplacedMeshesDiameter(
  int numMesh, /**< input - number of meshes */
  FC_Mesh* meshes,      /**< input - Array of meshes*/
  FC_Variable* displs,  /**< input - the coordinate displacements */ 
  double* diameter,   /**< ouput - diameter of the vertex set */
  int* numPair,       /**< output - number of pairs of vertices that are
		                    diameter distance from each other */
  int** pairIDs,      /**< output - IDs of the vertex pairs */
  FC_Mesh** pairMeshes  /**< output - the meshes of the pairs */
) {
  FC_ReturnCode rc;
  int i;
  FC_Subset* subsets, *pairSubs;

  // log message
  fc_printfLogMessage("Computing displaced meshes diameter.");

  // set default returns & test output args
  rc = _fc_getDiameterOutputArgsTestMeshes(diameter, numPair, pairIDs,
					   pairMeshes);
  if (rc != FC_SUCCESS)
    return rc;

  // test input
  if (numMesh < 1 || !meshes || !displs) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  for (i = 0; i < numMesh; i++) {
    if (!fc_isMeshValid(meshes[i]) || 
	!fc_isValidDisplacementVariable(meshes[i], displs[i])) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
      return FC_INPUT_ERROR;
    }
  }

  // Create subsets that are whole
  subsets = (FC_Subset*)malloc(numMesh*sizeof(FC_Subset));
  if (!subsets) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  for (i = 0; i < numMesh; i++) {
    rc = fc_createSubset(meshes[i], "temp", FC_AT_WHOLE_MESH, &subsets[i]);
    if (rc != FC_SUCCESS)
      return rc;
    rc = fc_addMemberToSubset(subsets[i], 0);
    if (rc != FC_SUCCESS)
      return rc;
  }

  // Send to core getSubsetDiameter code
  rc = _fc_getSubsetsDiameterCore(numMesh, subsets, displs, diameter, 
				  numPair, pairIDs, &pairSubs);

  // convert pairSubs to pairMeshes
  if (rc == FC_SUCCESS && numPair && *numPair > 0) {
    FC_Mesh* subMeshes;
    subMeshes = (FC_Mesh*)malloc(2*(*numPair)*sizeof(FC_Mesh));
    if (!subMeshes) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      return FC_MEMORY_ERROR;
    }
    for (i = 0; i < 2*(*numPair); i++) {
      rc = fc_getMeshFromSubset(pairSubs[i], &subMeshes[i]);
      if (rc != FC_SUCCESS)
      return rc;
    }
    *pairMeshes = subMeshes;
    free(pairSubs);
  }

  // Cleanup & return
  for (i = 0; i < numMesh; i++)
    fc_deleteSubset(subsets[i]);
  free(subsets);

  return rc;
}

/**
 * \ingroup  GeometricRelations
 * \brief Get the diameter of a group of subsets.
 * 
 * This routine returns the diameter of the set of vertices that make up the
 * members of a group of subsets. For more details see fc_getMeshesDiameter().
 *
 * \modifications
 *    - 3/9/06 WSD Created.
 */
FC_ReturnCode fc_getSubsetsDiameter(
  int numSubset,     /**< input - number of subsets */
  FC_Subset* subsets,  /**< input - Array of subsets */
  double* diameter,  /**< ouput - diameter of the vertex set */
  int* numPair,      /**< output - number of pairs of vertices that are
		                   diameter distance from each other */
  int** pairIDs,     /**< output - IDs of the vertex pairs */
  FC_Subset** pairSubsets /**< output - the subsets of the pairs */
) {
  FC_ReturnCode rc;
  int i;

  // log message
  fc_printfLogMessage("Computing subsets diameter.");

  // set default returns & test output args
  rc = _fc_getDiameterOutputArgsTestSubsets(diameter, numPair, pairIDs,
				     pairSubsets);
  if (rc != FC_SUCCESS)
    return rc;

  // test input
  if (numSubset < 1 || !subsets) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  for (i = 0; i < numSubset; i++) {
    if (!fc_isSubsetValid(subsets[i])) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
      return FC_INPUT_ERROR;
    }
  }

  // do it
  return _fc_getSubsetsDiameterCore(numSubset, subsets, NULL, diameter,
				    numPair, pairIDs, pairSubsets);
}

/**
 * \ingroup  GeometricRelations
 * \brief Get the diameter of a group of displaced subsets.
 * 
 * This routine returns the diameter of the set of displaced vertices that make
 * up the members of a gorup of subsets.  For more details see
 * fc_getMeshesDiameter().
 *
 * \modifications
 *    - 3/9/06 WSD Created.
 */
FC_ReturnCode fc_getDisplacedSubsetsDiameter(
  int numSubset,     /**< input - number of subsets */
  FC_Subset* subsets,  /**< input - Array of subset */
  FC_Variable* displs,  /**< input - the owning mesh's coordinate 
			                 displacements */ 
  double* diameter,  /**< ouput - diameter of the vertex set */
  int* numPair,      /**< output - number of pairs of vertices that are
		                   diameter distance from each other */
  int** pairIDs,     /**< output - IDs of the vertex pairs */
  FC_Subset** pairSubsets /**< output - the subsets of the pairs */
) {
  FC_ReturnCode rc;
  int i;
  FC_Mesh mesh;
 
  // log message
  fc_printfLogMessage("Computing displaced subsets diameter.");

  // set default returns & test output args
  rc = _fc_getDiameterOutputArgsTestSubsets(diameter, numPair, pairIDs,
					    pairSubsets);
  if (rc != FC_SUCCESS)
    return rc;

  // test input
  if (numSubset < 1 || !subsets || !displs) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  for (i = 0; i < numSubset; i++) {
    fc_getMeshFromSubset(subsets[i], &mesh);
    if (!fc_isSubsetValid(subsets[i]) || 
	!fc_isValidDisplacementVariable(mesh, displs[i])) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
      return FC_INPUT_ERROR;
    }
  }

  // do it
  return _fc_getSubsetsDiameterCore(numSubset, subsets, displs, diameter, 
				    numPair, pairIDs, pairSubsets);
}

//@}

/** \name Diameter helpers */
//-------------------------------
//@{

/** 
 * \ingroup  PrivateGeometricRelations
 * \brief Helper for 'getDiameter' routines.
 *
 * Sets the default return values of output args and tests that output
 * arg values are sane.
 *
 * \modifications
 *    - 6/27/05 WSD Created.
 */
FC_ReturnCode _fc_getDiameterOutputArgsTest(
  double* diameter,  /**< ouput - diameter of the vertex set */
  int* numPair,      /**< output - number of pairs of vertices that are
		                   diameter distance from each other */
  int** pairIDs      /**< output - IDs of the vertex pairs */
) {
  // default return
  if (diameter)
    *diameter = -1;
  if (numPair)
    *numPair = -1;
  if (pairIDs)
    *pairIDs = NULL;
  
  // test output args
  if (!diameter) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  if ((numPair && !pairIDs) || (!numPair && pairIDs)) {
    fc_printfErrorMessage("Did not provide both numPair & pairIDs for "
			  "returned vertex IDs");
    return FC_INPUT_ERROR;
  }
  
  return FC_SUCCESS;
}

/** 
 * \ingroup  PrivateGeometricRelations
 * \brief Helper for 'getDiameter' routines.
 *
 * Sets the default return values of output args and tests that output
 * arg values are sane.
 *
 * \modifications
 *    - 3/13/06 WSD Created.
 */
FC_ReturnCode _fc_getDiameterOutputArgsTestMeshes(
  double* diameter,  /**< ouput - diameter of the vertex set */
  int* numPair,      /**< output - number of pairs of vertices that are
		                   diameter distance from each other */
  int** pairIDs,     /**< output - IDs of the vertex pairs */
  FC_Mesh** pairMeshes /**< output - owning meshes for the verts */
) {
  // default return
  if (diameter)
    *diameter = -1;
  if (numPair)
    *numPair = -1;
  if (pairIDs)
    *pairIDs = NULL;
  if (pairMeshes)
    *pairMeshes = NULL;
  
  // test output args
  if (!diameter) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  if ( (numPair && (!pairIDs || !pairMeshes)) || 
       (!numPair && (pairIDs || pairMeshes)) ) {
    fc_printfErrorMessage("Did not provide both numPair & pairIDs for "
			  "returned vertex IDs");
    return FC_INPUT_ERROR;
  }
  
  return FC_SUCCESS;
}

/** 
 * \ingroup  PrivateGeometricRelations
 * \brief Helper for 'getDiameter' routines.
 *
 * Sets the default return values of output args and tests that output
 * arg values are sane.
 *
 * \modifications
 *    - 3/13/06 WSD Created.
 */
FC_ReturnCode _fc_getDiameterOutputArgsTestSubsets(
  double* diameter,  /**< ouput - diameter of the vertex set */
  int* numPair,      /**< output - number of pairs of vertices that are
		                   diameter distance from each other */
  int** pairIDs,     /**< output - IDs of the vertex pairs */
  FC_Subset** pairSubsets /**< output - owning subsets for the verts */
) {
  // default return
  if (diameter)
    *diameter = -1;
  if (numPair)
    *numPair = -1;
  if (pairIDs)
    *pairIDs = NULL;
  if (pairSubsets)
    *pairSubsets = NULL;
  
  // test output args
  if (!diameter) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  if ( (numPair && (!pairIDs || !pairSubsets)) || 
       (!numPair && (pairIDs || pairSubsets)) ) {
    fc_printfErrorMessage("Did not provide both numPair & pairIDs for "
			  "returned vertex IDs");
    return FC_INPUT_ERROR;
  }
  
  return FC_SUCCESS;
}

/** 
 * \ingroup  PrivateGeometricRelations
 * \brief Core code of 'getDiameter' routines.
 *
 * Implements brute force comparison of all vertices.
 *
 * \todo This routine will be very slow for point sets. The best optimization
 * is to find the convex hull of the set and then brute force search only the
 * points of the hull polygon. However, the 3D convex hull algorithm
 * is not easy to implement so I punted -WSD.
 *
 * \todo This routine should probably be replaced by
 *    _fc_getSubsetsDiameterCore() so that the core code isn't duplicated.
 *
 * \modifications
 *    - 6/27/05 WSD Created.
 */
FC_ReturnCode _fc_getSubsetDiameterCore(
  FC_Subset subset,    /**< input - a subset */
  FC_Variable* displ,  /**< input - pointer to displ (NULL = no displ) */
  double* diameter,    /**< ouput - diameter of the vertex set */
  int* numPair,        /**< output - number of pairs of vertices that are
		                     diameter distance from each other */
  int** pairIDs        /**< output - IDs of the vertex pairs */
) {
  FC_ReturnCode rc;
  int i, j;
  int numDim, numVertex, *vertexIDs, numMember, *memberIDs;
  double *coords_p, temp;
  FC_Mesh mesh;
  FC_AssociationType assoc;
  int doPairs = 0, num = 0, maxNum = 0, *ids;
  double diameter2 = -1; // squared diameter;
  
  // setup I
  if (numPair && pairIDs)
    doPairs = 1;
  rc = fc_getMeshFromSubset(subset, &mesh);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_getSubsetInfo(subset, &numMember, NULL, &assoc);
  if (rc != FC_SUCCESS)
    return rc;
  
  // get vertices of the members
  rc = fc_getSubsetMembersAsArray(subset, &numMember, &memberIDs);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_changeMeshEntityType(mesh, assoc, numMember, memberIDs, 
                               FC_AT_VERTEX, 1, &numVertex, &vertexIDs);
  free(memberIDs);

  // Special case: if 0 vertex - return error!
  if (numVertex == 0) {
    fc_printfErrorMessage("Subset has 0 members, diameter is undefined.");
    return FC_INPUT_ERROR;
  }
  // Special case: if 1 vertex - return 0
  if (numVertex == 1) {
    free(vertexIDs);
    *diameter = 0.;
    if (doPairs)
      *numPair = 0;
    return FC_SUCCESS;
  }

  // setup II
  rc = fc_getMeshDim(mesh, &numDim);
  if (rc != FC_SUCCESS)
    return rc;
  if (displ) {
    rc = fc_getDisplacedMeshCoords(mesh, *displ, &coords_p);
    if (rc != FC_SUCCESS)
      return rc;
  } 
  else {
    rc = fc_getMeshCoordsPtr(mesh, &coords_p);
    if (rc != FC_SUCCESS)
      return rc;
  }

  // Start off with space for at least 1 pair
  if (doPairs) {
    ids = (int*)malloc(2*sizeof(int));
    num = 0;
    maxNum = 1;
  }
  // Do it
  for (i = 0; i < numVertex; i++) {
    for (j = i+1; j < numVertex; j++) {
      fc_calcSquaredEuclideanDistance(&coords_p[vertexIDs[i]*numDim], 
				      &coords_p[vertexIDs[j]*numDim],
				      numDim, &temp);
      if (temp > diameter2) {
	diameter2 = temp;
	if (doPairs) {
	  ids[0] = vertexIDs[i];
	  ids[1] = vertexIDs[j];
	  num = 1;
	}
      }
      else if (doPairs && FC_DBL_EQUIV(temp, diameter2)) {
	if (num+1 > maxNum) {
	  ids = (int*)realloc(ids, 4*maxNum*sizeof(int));
	  if (!ids)
	    return FC_MEMORY_ERROR;
	  maxNum *= 2;
	}
	ids[num*2] = vertexIDs[i];
	ids[num*2+1] = vertexIDs[j];
	num++;
      }
    }
  }
  if (displ)
    free(coords_p);
  free(vertexIDs);

  // Return?
  *diameter = sqrt(diameter2);
  if (doPairs) {
    *numPair = num;
    *pairIDs = (int*)realloc(ids, 2*num*sizeof(int)); // squeeze
  }

  return FC_SUCCESS;
}

/** 
 * \ingroup  PrivateGeometricRelations
 * \brief Core code of 'getDiameter' routines.
 *
 * Implements brute force comparison of all vertices.
 *
 * \todo This routine will be very slow for point sets. The best optimization
 * is to find the convex hull of the set and then brute force search only the
 * points of the hull polygon. However, the 3D convex hull algorithm
 * is not easy to implement so I punted -WSD.
 *
 * \modifications
 *    - 6/27/05 WSD Created.
 *    - 3/9/06 WSD Changed to do multiple subsets.
 */
FC_ReturnCode _fc_getSubsetsDiameterCore(
  int numSubset,       /**< input - number of subsets */
  FC_Subset* subsets,  /**< input - a subset */
  FC_Variable* displs, /**< input - array of displs (NULL = no displ) */
  double* diameter,    /**< ouput - diameter of the vertex set */
  int* numPair,        /**< output - number of pairs of vertices that are
		                     diameter distance from each other */
  int** pairIDs,       /**< output - IDs of the vertex pairs */
  FC_Subset** pairSubsets /**< output - Subsets of the vertex pairs */
) {
  FC_ReturnCode rc;
  int i, j, i_s, j_s;
  int numDim, numVerts[numSubset], *vertIDs[numSubset], numTotalVert;
  int numMember, *memberIDs;
  double *coords[numSubset], temp;
  FC_Mesh meshes[numSubset];
  FC_Subset* subs = NULL, *temp_subs;
  FC_AssociationType assoc;
  int doPairs = 0, num = 0, maxNum = 1, *ids = NULL, *temp_ids;
  double diameter2 = -1; // squared diameter
  
  // assume args are o.k.
  // FIX: assume numDim is the same for all subsets,
  // pass as an argument!?

  // setup I
  if (numPair && pairSubsets && pairIDs)
    doPairs = 1;

  // Get vertIDs per subset
  numTotalVert = 0;
  for (i = 0; i < numSubset; i++) {
    rc = fc_getMeshFromSubset(subsets[i], &meshes[i]);
    if (rc != FC_SUCCESS)
      return rc;
    rc = fc_getMeshDim(meshes[i], &numDim); // FIX
    if (rc != FC_SUCCESS)
      return rc; 
    rc = fc_getSubsetInfo(subsets[i], &numMember, NULL, &assoc);
    if (rc != FC_SUCCESS)
      return rc;
    // get vertices of the members
    rc = fc_getSubsetMembersAsArray(subsets[i], &numMember, &memberIDs);
    if (rc != FC_SUCCESS)
      return rc;
    if (assoc == FC_AT_VERTEX) {
      numVerts[i] = numMember;
      vertIDs[i] = memberIDs;
    }
    else {
      rc = fc_changeMeshEntityType(meshes[i], assoc, numMember, memberIDs, 
				   FC_AT_VERTEX, 1, &numVerts[i], &vertIDs[i]);
      free(memberIDs);
      if (rc != FC_SUCCESS)
	return rc;
    }
    numTotalVert += numVerts[i];
  }

  // Special case: if 0 vertex - return error!
  if (numTotalVert == 0) {
    fc_printfErrorMessage("Subsets have 0 members, diameter is undefined.");
    return FC_INPUT_ERROR;
  }
  // Special case: if 1 vertex - return 0
  if (numTotalVert == 1) {
    for (i = 0; i < numSubset; i++)
      free(vertIDs[i]);
    *diameter = 0.;
    if (doPairs)
      *numPair = 0;
    return FC_SUCCESS;
  }

  // Get the coords
  for (i = 0; i < numSubset; i++) {
    if (displs) 
      rc = fc_getDisplacedMeshCoords(meshes[i], displs[i], &coords[i]);
    else 
      rc = fc_getMeshCoordsPtr(meshes[i], &coords[i]);
    if (rc != FC_SUCCESS)
      return rc;
  }

  // Start off with space for at least 1 pair
  if (doPairs) {
    ids = (int*)malloc(2*sizeof(int));
    subs = (FC_Subset*)malloc(2*sizeof(FC_Subset));
    num = 0;
    maxNum = 1;
  }
  // Note: diameter2 starts at -1
  // Do within each subset
  for (i_s = 0; i_s < numSubset; i_s++) {
    for (i = 0; i < numVerts[i_s]; i++) {
      for (j = i+1; j < numVerts[i_s]; j++) {
	fc_calcSquaredEuclideanDistance(&coords[i_s][vertIDs[i_s][i]*numDim], 
					&coords[i_s][vertIDs[i_s][j]*numDim],
					numDim, &temp);
	if (temp > diameter2) {
	  diameter2 = temp;
	  if (doPairs) {
	    ids[0] = vertIDs[i_s][i];
	    ids[1] = vertIDs[i_s][j];
	    subs[0] = subsets[i_s];
	    subs[1] = subsets[i_s];
	    num = 1;
	  }
	}
	else if (doPairs && FC_DBL_EQUIV(temp, diameter2)) {
	  if (num+1 > maxNum) {
	    temp_ids = (int*)realloc(ids, 4*maxNum*sizeof(int));
	    temp_subs = (FC_Subset*)realloc(subs, 4*maxNum*sizeof(FC_Subset));
	    if (!temp_ids || !temp_subs)
	      return FC_MEMORY_ERROR;
	    ids = temp_ids;
	    subs = temp_subs;
	    maxNum *= 2;
	  }
	  ids[num*2] = vertIDs[i_s][i];
	  ids[num*2+1] = vertIDs[i_s][j];
	  subs[num*2] = subsets[i_s];
	  subs[num*2+1] = subsets[i_s];
	  num++;
	}
      }
    }
  }
  // check each vertex in one subset against all verts in other subsets
  for (i_s = 0; i_s < numSubset; i_s++) {
    for (i = 0; i < numVerts[i_s]; i++) {
      for (j_s = i_s+1; j_s < numSubset; j_s++) {
	for (j = 0; j < numVerts[j_s]; j++) {
	  fc_calcSquaredEuclideanDistance(&coords[i_s][vertIDs[i_s][i]*numDim], 
					  &coords[j_s][vertIDs[j_s][j]*numDim],
					  numDim, &temp);
	  if (temp > diameter2) {
	    diameter2 = temp;
	    if (doPairs) {
	      ids[0] = vertIDs[i_s][i];
	      ids[1] = vertIDs[j_s][j];
	      subs[0] = subsets[i_s];
	      subs[1] = subsets[j_s];
	      num = 1;
	    }
	  }
	  else if (doPairs && FC_DBL_EQUIV(temp, diameter2)) {
	    if (num+1 > maxNum) {
	      temp_ids = (int*)realloc(ids, 4*maxNum*sizeof(int));
	      temp_subs = (FC_Subset*)realloc(subs, 4*maxNum*sizeof(FC_Subset));
	      if (!temp_ids || !temp_subs)
		return FC_MEMORY_ERROR;
	      maxNum*= 2;
	      ids = temp_ids;
	      subs = temp_subs;
	    }
	    ids[num*2] = vertIDs[i_s][i];
	    ids[num*2+1] = vertIDs[j_s][j];
	    subs[num*2] = subsets[i_s];
	    subs[num*2+1] = subsets[j_s];
	    num++;
	  }
	}
      }
    }
  }
  if (displs) 
    for (i = 0; i < numSubset; i++)
      free(coords[i]);
  for (i = 0; i < numSubset; i++)
    free(vertIDs[i]);
  
  // Return?
  *diameter = sqrt(diameter2);
  if (doPairs) {
    *numPair = num;
    *pairSubsets = (FC_Subset*)realloc(subs, 2*num*sizeof(FC_Subset));
    *pairIDs = (int*)realloc(ids, 2*num*sizeof(int)); // squeeze
  }

  return FC_SUCCESS;
}


//@}

/** \name Centroid */
//-------------------------------
//@{

/**
 * \ingroup  GeometricRelations
 * \brief Compute centroid of mesh.
 *
 * \description
 *  
 *    The centroid is the geometric center of the mesh and is calculated by
 *    averaging the coordinates in each dimension. Also called center
 *    of mass or center of gravity.
 *
 * \modifications   
 *   - Wendy Koegler  Created.
 *   - July 3, 2002 W. Koegler removed dependence on var
 *   - July 3, 2002  W Koegler  Adapted to modified _FC_Vertex
 *   - July 9, 2002  W Koegler  Adapted to modified struct Vector
 *   - 2003-AUG-05  W Koegler  Changed to use doubles internally
 *   - 10/29/04 WSK, added 'mesh' to name & added dim to returned arg.
 */
FC_ReturnCode fc_getMeshCentroid(
  FC_Mesh mesh,        /**< input - mesh id */
  int *dim,            /**< output - dimensionality of coords */
  FC_Coords *center    /**< output - coordinates of centroid */
) {
  FC_ReturnCode rc;
  int i, j;
  int numDim, numVertex;
  double *coords_p;
  double coords_sum[3];
  
  // default
  if (dim)
    *dim = -1;
  if (center) {
    for (i = 0; i < 3; i++)
      (*center)[i] = -1.;
  }

  // check input
  if (!fc_isMeshValid(mesh) || !dim || !center) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // log message
  fc_printfLogMessage("Computing mesh centroid.");  
  
  // get vertices 
  rc = fc_getMeshInfo(mesh, NULL, &numDim, &numVertex, NULL, NULL);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_getMeshCoordsPtr(mesh, &coords_p);
  if (rc != FC_SUCCESS) 
    return rc;
  
  // calculate coordinate sums
  for (i = 0; i < 3; i++)
    coords_sum[i] = 0.;
  
  for (i = 0; i < numVertex; i++)
    for (j = 0; j < numDim; j++)
      coords_sum[j] += coords_p[i*numDim + j];
  
  // calculate centroids
  *dim = numDim;
  for (i = 0; i < numDim; i++)
    (*center)[i] = coords_sum[i]/numVertex;
  for (i = numDim; i < 3; i++)
    (*center)[i] = 0.;
  
  return FC_SUCCESS;
}

/**
 * \ingroup  GeometricRelations
 * \brief Compute centroid of displaced mesh.
 *
 * \description
 *  
 *    The centroid is the geometric center of the mesh and is calculated by
 *    averaging the coordinates in each dimension.
 *
 * \modifications   
 *   - 10/29/04 WSK, created.
 */
FC_ReturnCode fc_getDisplacedMeshCentroid(
  FC_Mesh mesh,             /**< input - mesh */
  FC_Variable coord_displ,  /**< input - coord displacements */
  int *dim,                 /**< output - dimensionality of coords */
  FC_Coords *center         /**< output - coordinates of centroid */
) {
  FC_ReturnCode rc;
  int i, j;
  int numVertex, numDim;
  double *coords;
  double coords_sum[3];
  
  // default output
  if (dim)
    *dim = -1;
  if (center) {
    for (i = 0; i < 3; i++)
      (*center)[i] = -1.;
  }

  // check input
  if (!fc_isMeshValid(mesh) || 
      !fc_isValidDisplacementVariable(mesh, coord_displ) || 
      !dim || !center) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // log message
  fc_printfLogMessage("Computing displaced mesh centroid.");  
  
  // get vertices 
  rc = fc_getMeshInfo(mesh, NULL, &numDim, &numVertex, NULL, NULL);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_getDisplacedMeshCoords(mesh, coord_displ, &coords);
  if (rc != FC_SUCCESS) 
    return rc;
  
  // initialize sums
  for (i = 0; i < 3; i++)
    coords_sum[i] = 0.;
  
  // collect sums
  for (i = 0; i < numVertex; i++)
    for (j = 0; j < numDim; j++)
      coords_sum[j] += coords[i*numDim + j];
  free(coords);
  
  // calculate centroids
  *dim = numDim;
  for (i = 0; i < numDim; i++)
    (*center)[i] = coords_sum[i]/numVertex;
  for (i = numDim; i < 3; i++)
    (*center)[i] = 0.;

  return FC_SUCCESS;
}

/**
 * \ingroup  GeometricRelations
 * \brief Compute variable weighted centroid of a mesh.
 *
 * \description  
 *
 *    Computes the centroid of a mesh weighted by the variable values
 *    (equivalent to center of mass where the variable is the masses). The
 *    weights variable must be a scalar variable on the vertices (and should have
 *    the same number of vertices as the mesh if it is from a different mesh).
 *
 * \modifications   
 *    - Wendy Koegler  Created.
 *    - July 3, 2002  W Koegler  Created, work in progress
 *    - July 3, 2002  W Koegler  Adapted to modified _FC_Vertex
 *    - July 9, 2002  W Koegler  Adapted to modified struct Vector
 *    - 10/29/04 WSK added mesh to input args and dim to output args.
 *        Fixed up.
 */
FC_ReturnCode fc_getVariableWeightedMeshCentroid(
  FC_Mesh mesh,            /**< input - mesh */
  FC_Variable weightsVar,  /**< input - weights */
  int* dim,                /**< output - dimensionality of centroid */
  FC_Coords *center        /**< output - coordinates of centroid */
) {
  FC_ReturnCode rc;
  int i, j;
  int numDim, numVertex, numDataPoint, data_dim;
  FC_AssociationType assoc; // data association
  FC_DataType datatype; // data type
  double coords_sum[3];
  double* coords_p, weight_sum, weight;
  void* data_p;

  // default output
  if (dim)
    *dim = -1;
  if (center) {
    for (i = 0; i < 3; i++)
      (*center)[i] = -1.;
  }

  // check input
  if (!fc_isMeshValid(mesh) || !fc_isVariableValid(weightsVar) || !dim || 
      !center) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // get mesh & variable information
  rc = fc_getMeshInfo(mesh, NULL, &numDim, &numVertex, NULL, NULL);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_getVariableInfo(weightsVar, &numDataPoint, &data_dim, &assoc,  
			  NULL, &datatype);
  if (rc != FC_SUCCESS)
    return rc;
  
   // more checking
  if (assoc != FC_AT_VERTEX || numDataPoint != numVertex || data_dim != 1 ||
      !fc_isDataTypeValid(datatype) || datatype == FC_DT_UNKNOWN || 
      datatype == FC_DT_CHAR) {
    fc_printfErrorMessage("Inappropriate variable for weights");
    return FC_ERROR;
  }

  // log message
  fc_printfLogMessage("Computing variable weighted centroid of mesh.");  

  // get vertices and data
  rc = fc_getMeshCoordsPtr(mesh, &coords_p);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_getVariableDataPtr(weightsVar, &data_p);
  if (rc != FC_SUCCESS)
    return rc;


  // init coord sums
  weight_sum = 0.;
  for (i = 0; i < 3; i++)
    coords_sum[i] = 0;

  // collect sums
  for (i = 0; i < numVertex; i++) {
    switch(datatype) {
    case FC_DT_INT :    weight = ((int*)data_p)[i];     break;
    case FC_DT_FLOAT:   weight = ((float*)data_p)[i];   break;
    case FC_DT_DOUBLE:  weight = ((double*)data_p)[i];  break;
    default:
      ; // nothing
    }
    weight_sum += weight;
    for (j = 0; j < numDim; j++) 
      coords_sum[j] += coords_p[i*numDim + j] * weight;
  }

  // calculate weight centroid
  *dim = numDim;
  for (i = 0; i < numDim; i++)
    (*center)[i] = coords_sum[i]/weight_sum;
  for (i = numDim; i < 3; i++)
    (*center)[i] = 0.;

  return FC_SUCCESS;
}

/**
 * \ingroup  GeometricRelations
 * \brief Compute variable weighted centroid of a displaced mesh.
 *
 * \description  
 *
 *    Computes the centroid of a mesh weighted by the variable values
 *    (equivalent to center of mass where the variable is the masses). The
 *    weights variable must be a scalar variable on the vertices (and should have
 *    the same number of vertices as the mesh if it is from a different mesh).
 *
 * \modifications   
 *    - 10/29/04 WSK created.
 */
FC_ReturnCode fc_getVariableWeightedDisplacedMeshCentroid(
  FC_Mesh mesh,            /**< input - mesh */
  FC_Variable coord_displ,  /**< input - coord displacements */
  FC_Variable weightsVar,  /**< input - weights */
  int* dim,                /**< output - dimensionality of centroid */
  FC_Coords *center        /**< output - coordinates of centroid */
) {
  FC_ReturnCode rc;
  int i, j;
  int numDim, numVertex, numDataPoint, data_dim;
  FC_AssociationType assoc; // data association
  FC_DataType datatype; // data type
  double coords_sum[3];
  double* coords, weight_sum, weight;
  void* data_p;

  // default output
  if (dim)
    *dim = -1;
  if (center) {
    for (i = 0; i < 3; i++)
      (*center)[i] = -1.;
  }

  // check input
  if (!fc_isMeshValid(mesh) || 
      !fc_isValidDisplacementVariable(mesh, coord_displ) ||
      !fc_isVariableValid(weightsVar) || !dim || !center) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // get mesh & weights variable information
  rc = fc_getMeshInfo(mesh, NULL, &numDim, &numVertex, NULL, NULL);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_getVariableInfo(weightsVar, &numDataPoint, &data_dim, &assoc,  
			  NULL, &datatype);
  if (rc != FC_SUCCESS)
    return rc;
  
   // more checking
  if (assoc != FC_AT_VERTEX || numDataPoint != numVertex || data_dim != 1 ||
      !fc_isDataTypeValid(datatype) || datatype == FC_DT_UNKNOWN || 
      datatype == FC_DT_CHAR) {
    fc_printfErrorMessage("Inappropriate variable for weights");
    return FC_ERROR;
  }

  // log message
  fc_printfLogMessage("Computing variable weighted centroid of mesh.");  

  // get vertices and data
  rc = fc_getDisplacedMeshCoords(mesh, coord_displ, &coords);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_getVariableDataPtr(weightsVar, &data_p);
  if (rc != FC_SUCCESS)
    return rc;

  // init coord sums
  weight_sum = 0.;
  for (i = 0; i < 3; i++)
    coords_sum[i] = 0;

  // collect sums
  for (i = 0; i < numVertex; i++) {
    switch(datatype) {
    case FC_DT_INT :    weight = ((int*)data_p)[i];     break;
    case FC_DT_FLOAT:   weight = ((float*)data_p)[i];   break;
    case FC_DT_DOUBLE:  weight = ((double*)data_p)[i];  break;
    default:
      ; // nothing
    }
    weight_sum += weight;
    for (j = 0; j < numDim; j++) {
      coords_sum[j] += coords[i*numDim + j] * weight;
    }
  }
  free(coords);

  // calculate weight centroid
  *dim = numDim;
  for (i = 0; i < numDim; i++)
    (*center)[i] = coords_sum[i]/weight_sum;
  for (i = numDim; i < 3; i++)
    (*center)[i] = 0.;

  return FC_SUCCESS;
}

/**
 * \ingroup  GeometricRelations
 * \brief Compute centroid of a subset.
 *
 * \description
 *  
 *    The centroid is the geometric center of the subset and is calculated by
 *    averaging the coordinates in each dimension. Also called center
 *    of mass or center of gravity.
 *
 * \modifications   
 *    - 10/29/04 WSK Created.
 */
FC_ReturnCode fc_getSubsetCentroid(
  FC_Subset subset,    /**< input - subset */
  int *dim,            /**< output - dimensionality of coords */
  FC_Coords *center    /**< output - coordinates of centroid */
) {
  FC_ReturnCode rc;
  int i, j;
  FC_Mesh mesh;
  FC_AssociationType assoc;
  int numDim, numVertex, *vertexIDs, numMember, *memberIDs;
  double *coords_p;
  double coords_sum[3];
  
  // default
  if (dim)
    *dim = -1;
  if (center) {
    for (i = 0; i < 3; i++)
      (*center)[i] = -1.;
  }

  // check input
  if (!fc_isSubsetValid(subset) || !dim || !center) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // log message
  fc_printfLogMessage("Computing subset centroid.");  
  
  // setup I
  rc = fc_getMeshFromSubset(subset, &mesh);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_getSubsetInfo(subset, &numMember, NULL, &assoc);
  if (rc != FC_SUCCESS)
    return rc;

  // Error case: if subset has zero members - return error
  if (numMember == 0 ) {
    fc_printfErrorMessage("Subset has 0 members");
    return FC_ERROR;
  }
  
  // special case -- whole association
  if (assoc == FC_AT_WHOLE_MESH) 
    return fc_getMeshCentroid(mesh, dim, center);
    
  // setup II
  rc = fc_getMeshDim(mesh, &numDim);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_getMeshCoordsPtr(mesh, &coords_p);
  if (rc != FC_SUCCESS)
    return rc;
  
  // get vertices of the members
  rc = fc_getSubsetMembersAsArray(subset, &numMember, &memberIDs);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_changeMeshEntityType(mesh, assoc, numMember, memberIDs, 
                               FC_AT_VERTEX, 1, &numVertex, &vertexIDs);
  free(memberIDs);

  // calculate coordinate sums
  for (i = 0; i < 3; i++)
    coords_sum[i] = 0.;
  
  for (i = 0; i < numVertex; i++)
    for (j = 0; j < numDim; j++)
      coords_sum[j] += coords_p[vertexIDs[i]*numDim + j];
  free(vertexIDs);
  
  // calculate centroids
  *dim = numDim;
  for (i = 0; i < numDim; i++)
    (*center)[i] = coords_sum[i]/numVertex;
  for (i = numDim; i < 3; i++)
    (*center)[i] = 0.;
  
  return FC_SUCCESS;
}

/**
 * \ingroup  GeometricRelations
 * \brief Compute centroid of a displaced subset.
 *
 * \description
 *  
 *    The centroid is the geometric center of the subset and is calculated by
 *    averaging the coordinates in each dimension. Also called center
 *    of mass or center of gravity.
 *
 * \modifications   
 *    - 10/29/04 WSK Created.
 */
FC_ReturnCode fc_getDisplacedSubsetCentroid(
  FC_Subset subset,    /**< input - subset */
  FC_Variable coord_displ,  /**< input - coord displacements */
  int *dim,            /**< output - dimensionality of coords */
  FC_Coords *center    /**< output - coordinates of centroid */
) {
  FC_ReturnCode rc;
  int i, j;
  FC_Mesh mesh;
  FC_AssociationType assoc;
  int numDim, numVertex, *vertexIDs, numMember, *memberIDs;
  double *coords;
  double coords_sum[3];
  
  // default
  if (dim)
    *dim = -1;
  if (center) {
    for (i = 0; i < 3; i++)
      (*center)[i] = -1.;
  }

  // check input
  if (!fc_isSubsetValid(subset) || !dim || !center) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  rc = fc_getMeshFromSubset(subset, &mesh);
  if (rc != FC_SUCCESS)
    return rc;
  if (!fc_isValidDisplacementVariable(mesh, coord_displ)) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // log message
  fc_printfLogMessage("Computing displaced subset centroid.");  
  
  // setup I
  rc = fc_getSubsetInfo(subset, &numMember, NULL, &assoc);
  if (rc != FC_SUCCESS)
    return rc;

  // Error case: if subset has zero members - return error
  if (numMember == 0 ) {
    fc_printfErrorMessage("Subset has 0 members");
    return FC_ERROR;
  }
  
  // special case -- whole association
  if (assoc == FC_AT_WHOLE_MESH) 
    return fc_getDisplacedMeshCentroid(mesh, coord_displ, dim, center);
    
  // setup II
  rc = fc_getMeshDim(mesh, &numDim);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_getDisplacedMeshCoords(mesh, coord_displ, &coords);
  if (rc != FC_SUCCESS)
    return rc;
  
  // get vertices of the members
  rc = fc_getSubsetMembersAsArray(subset, &numMember, &memberIDs);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_changeMeshEntityType(mesh, assoc, numMember, memberIDs, 
                               FC_AT_VERTEX, 1, &numVertex, &vertexIDs);
  free(memberIDs);

  // calculate coordinate sums
  for (i = 0; i < 3; i++)
    coords_sum[i] = 0.;
  
  for (i = 0; i < numVertex; i++)
    for (j = 0; j < numDim; j++)
      coords_sum[j] += coords[vertexIDs[i]*numDim + j];
  free(vertexIDs);
  free(coords);
  
  // calculate centroids
  *dim = numDim;
  for (i = 0; i < numDim; i++)
    (*center)[i] = coords_sum[i]/numVertex;
  for (i = numDim; i < 3; i++)
    (*center)[i] = 0.;
  
  return FC_SUCCESS;
}

/**
 * \ingroup  GeometricRelations
 * \brief Compute centroid of a subset.
 *
 * \description
 *  
 *    The centroid is the geometric center of the subset and is calculated by
 *    averaging the coordinates in each dimension. Also called center
 *    of mass or center of gravity.
 *
 * \modifications   
 *    - 10/29/04 WSK Created.
 */
FC_ReturnCode fc_getVariableWeightedSubsetCentroid(
  FC_Subset subset,    /**< input - subset */
  FC_Variable weightsVar,  /**< input - weights */
  int *dim,            /**< output - dimensionality of coords */
  FC_Coords *center    /**< output - coordinates of centroid */
) {
  FC_ReturnCode rc;
  int i, j;
  FC_Mesh mesh;
  FC_AssociationType var_assoc, subset_assoc;
  FC_DataType datatype; // data type
  int numDim, numVertexInMesh, numDataPoint, data_dim;
  int numVertex, *vertexIDs, numMember, *memberIDs;
  double coords_sum[3];
  double *coords_p, weight_sum, weight;
  void* data_p;

  // default
  if (dim)
    *dim = -1;
  if (center) {
    for (i = 0; i < 3; i++)
      (*center)[i] = -1.;
  }

  // check input
  if (!fc_isSubsetValid(subset) || !fc_isVariableValid(weightsVar) || !dim 
      || !center) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // get mesh & variable information
  rc = fc_getMeshFromSubset(subset, &mesh);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_getMeshInfo(mesh, NULL, &numDim, &numVertexInMesh, NULL, NULL);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_getVariableInfo(weightsVar, &numDataPoint, &data_dim, &var_assoc,  
			  NULL, &datatype);
  if (rc != FC_SUCCESS)
    return rc;
  
  // more checking
  if (var_assoc != FC_AT_VERTEX || numDataPoint != numVertexInMesh || 
      data_dim != 1 || !fc_isDataTypeValid(datatype) || 
      datatype == FC_DT_UNKNOWN || datatype == FC_DT_CHAR) {
    fc_printfErrorMessage("Inappropriate variable for weights");
    return FC_ERROR;
  }

  // log message
  fc_printfLogMessage("Computing variable weighted subset centroid.");  
  
  // setup I
  rc = fc_getSubsetInfo(subset, &numMember, NULL, &subset_assoc);
  if (rc != FC_SUCCESS)
    return rc;

  // Error case: if subset has zero members - return error
  if (numMember == 0 ) {
    fc_printfErrorMessage("Subset has 0 members");
    return FC_ERROR;
  }
  
  // special case -- whole association
  if (subset_assoc == FC_AT_WHOLE_MESH) 
    return fc_getVariableWeightedMeshCentroid(mesh, weightsVar, dim, center);
    
  // setup II
  rc = fc_getMeshCoordsPtr(mesh, &coords_p);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_getVariableDataPtr(weightsVar, &data_p);
  if (rc != FC_SUCCESS)
    return rc;
  
  // get vertices of the members
  rc = fc_getSubsetMembersAsArray(subset, &numMember, &memberIDs);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_changeMeshEntityType(mesh, subset_assoc, numMember, memberIDs, 
                               FC_AT_VERTEX, 1, &numVertex, &vertexIDs);
  free(memberIDs);

  // init coordinate sums
  weight_sum = 0.;
  for (i = 0; i < 3; i++)
    coords_sum[i] = 0.;
  
  // collect sums
  for (i = 0; i < numVertex; i++) {
    switch(datatype) {
    case FC_DT_INT :    weight = ((int*)data_p)[vertexIDs[i]];     break;
    case FC_DT_FLOAT:   weight = ((float*)data_p)[vertexIDs[i]];   break;
    case FC_DT_DOUBLE:  weight = ((double*)data_p)[vertexIDs[i]];  break;
    default:
      ; // nothing
    }
    weight_sum += weight;
    for (j = 0; j < numDim; j++)
      coords_sum[j] += coords_p[vertexIDs[i]*numDim + j] * weight;
  }
  free(vertexIDs);
  
  // calculate centroids
  *dim = numDim;
  for (i = 0; i < numDim; i++)
    (*center)[i] = coords_sum[i]/weight_sum;
  for (i = numDim; i < 3; i++)
    (*center)[i] = 0.;
  
  return FC_SUCCESS;
}

/**
 * \ingroup  GeometricRelations
 * \brief Compute centroid of a displaced subset.
 *
 * \description
 *  
 *    The centroid is the geometric center of the subset and is calculated by
 *    averaging the coordinates in each dimension. Also called center
 *    of mass or center of gravity.
 *
 * \modifications   
 *    - 10/29/04 WSK Created.
 */
FC_ReturnCode fc_getVariableWeightedDisplacedSubsetCentroid(
  FC_Subset subset,    /**< input - subset */
  FC_Variable coord_displ, /**< input - coord displacements */
  FC_Variable weightsVar,  /**< input - weights */
  int *dim,            /**< output - dimensionality of coords */
  FC_Coords *center    /**< output - coordinates of centroid */
) {
  FC_ReturnCode rc;
  int i, j;
  FC_Mesh mesh;
  FC_AssociationType var_assoc, subset_assoc;
  FC_DataType datatype; // data type
  int numDim, numVertexInMesh, numDataPoint, data_dim;
  int numVertex, *vertexIDs, numMember, *memberIDs;
  double coords_sum[3];
  double *coords, weight_sum, weight;
  void* data_p;

  // default
  if (dim)
    *dim = -1;
  if (center) {
    for (i = 0; i < 3; i++)
      (*center)[i] = -1.;
  }

  // check input
  if (!fc_isSubsetValid(subset) || !fc_isVariableValid(weightsVar) ||
      !dim || !center) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  rc = fc_getMeshFromSubset(subset, &mesh);
  if (rc != FC_SUCCESS)
    return rc;
  if (!fc_isValidDisplacementVariable(mesh, coord_displ)) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // get mesh & weights variable information
  rc = fc_getMeshInfo(mesh, NULL, &numDim, &numVertexInMesh, NULL, NULL);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_getVariableInfo(weightsVar, &numDataPoint, &data_dim, &var_assoc,  
			  NULL, &datatype);
  if (rc != FC_SUCCESS)
    return rc;
  
  // more checking
  if (var_assoc != FC_AT_VERTEX || numDataPoint != numVertexInMesh || 
      data_dim != 1 || !fc_isDataTypeValid(datatype) || 
      datatype == FC_DT_UNKNOWN || datatype == FC_DT_CHAR) {
    fc_printfErrorMessage("Inappropriate variable for weights");
    return FC_ERROR;
  }

  // log message
  fc_printfLogMessage("Computing variable weighted subset centroid.");  
  
  // setup I
  rc = fc_getSubsetInfo(subset, &numMember, NULL, &subset_assoc);
  if (rc != FC_SUCCESS)
    return rc;

  // Error case: if subset has zero members - return error
  if (numMember == 0 ) {
    fc_printfErrorMessage("Subset has 0 members");
    return FC_ERROR;
  }
  
  // special case -- whole association
  if (subset_assoc == FC_AT_WHOLE_MESH) 
    return fc_getVariableWeightedDisplacedMeshCentroid(mesh, coord_displ,
						   weightsVar, dim, center);
    
  // setup II
  rc = fc_getDisplacedMeshCoords(mesh, coord_displ, &coords);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_getVariableDataPtr(weightsVar, &data_p);
  if (rc != FC_SUCCESS)
    return rc;
  
  // get vertices of the members
  rc = fc_getSubsetMembersAsArray(subset, &numMember, &memberIDs);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_changeMeshEntityType(mesh, subset_assoc, numMember, memberIDs, 
                               FC_AT_VERTEX, 1, &numVertex, &vertexIDs);
  free(memberIDs);

  // init coordinate sums
  weight_sum = 0.;
  for (i = 0; i < 3; i++)
    coords_sum[i] = 0.;
  
  // collect sums
  for (i = 0; i < numVertex; i++) {
    switch(datatype) {
    case FC_DT_INT :    weight = ((int*)data_p)[vertexIDs[i]];     break;
    case FC_DT_FLOAT:   weight = ((float*)data_p)[vertexIDs[i]];   break;
    case FC_DT_DOUBLE:  weight = ((double*)data_p)[vertexIDs[i]];  break;
    default:
      ; // nothing
    }
    weight_sum += weight;
    for (j = 0; j < numDim; j++)
      coords_sum[j] += coords[vertexIDs[i]*numDim + j] * weight;
  }
  free(vertexIDs);
  free(coords);
  
  // calculate centroids
  *dim = numDim;
  for (i = 0; i < numDim; i++)
    (*center)[i] = coords_sum[i]/weight_sum;
  for (i = numDim; i < 3; i++)
    (*center)[i] = 0.;
  
  return FC_SUCCESS;
}

//@}

/** \name Mesh entities' measurements */
//-------------------------------
//@{


/**
 * \ingroup  GeometricRelations
 * \brief Return a variable with the edge lengths of the mesh
 *
 * \description
 *
 *    If the mesh is made of 1D elements (i.e. lines)only lines, the length
 *    variable will be associated with the elements. Otherwise it will be 
 *    associated with the edges of the mesh.
 *
 *    This routine will return with an error if the elements of the mesh are 
 *    not at least 2D (e.g. it is an error to ask for lengths of a point mesh)
 *
 * \modifications  
 *    - 02/03/04 WSK, Created to be more generally useful than
 *       fc_createEdgeMesh
 *    - 08/23/04 WSK, Fixed bug: was failing for line mesh.
 */
FC_ReturnCode fc_getEdgeLengths(
  FC_Mesh mesh,              /**< input - mesh */
  FC_Variable *edge_lengths  /**< output - edge lengths variable handle */
) {
  FC_ReturnCode rc;
  int i;
  int dim, numEdge;
  double *lengths;
  FC_AssociationType assoc;
  FC_ElementType elemType;
  double* coords;
  int* edgeConns;
  
  // default return
  if (edge_lengths) 
    *edge_lengths = FC_NULL_VARIABLE;

  // check input
  if (!fc_isMeshValid(mesh) || edge_lengths == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
    
  // log message
  fc_printfLogMessage("Getting the edge lengths of the mesh.");  

  // setup & more checking
  rc = fc_getMeshInfo(mesh, NULL, &dim, NULL, NULL, &elemType);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_getMeshCoordsPtr(mesh, &coords);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_getMeshNumEdge(mesh, &numEdge);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_getMeshEdgeConnsPtr(mesh, &edgeConns);
  if (rc != FC_SUCCESS)
    return rc;
  if (numEdge < 1 && elemType != FC_ET_LINE) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // compute edge lengths
  lengths = (double *)malloc(sizeof(double) * numEdge);
  if (lengths == NULL ) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  for (i = 0; i < numEdge; i++) 
    fc_calcEuclideanDistance(&coords[edgeConns[i*2]*dim],
                             &coords[edgeConns[i*2 + 1]*dim], 
                             dim, &lengths[i]);

  // create new variable
  if (elemType == FC_ET_LINE)
    assoc = FC_AT_ELEMENT;
  else
    assoc = FC_AT_EDGE;
  rc = fc_createVariable(mesh, "edge_lengths", edge_lengths);
  if (rc != FC_SUCCESS) {
    free(lengths);
    return rc;
  }
  rc = fc_setVariableDataPtr(*edge_lengths, numEdge, 1, assoc, FC_MT_SCALAR, 
                             FC_DT_DOUBLE, lengths);
  // don't free lengths!
  
  return rc;
}

/**
 * \ingroup  GeometricRelations
 * \brief Return a variable with the edge lengths of the displaced mesh
 *
 * \description
 *
 *    If the mesh is made of 1D elements (i.e. lines)only lines, the length
 *    variable will be associated with the elements. Otherwise it will be 
 *    associated with the edges of the mesh.
 *
 *    This routine will return with an error if the elements of the mesh are 
 *    not at least 2D (e.g. it is an error to ask for lengths of a point mesh)
 *
 * \modifications  
 *    - 11/01/04 WSK created.
 */
FC_ReturnCode fc_getDisplacedEdgeLengths(
  FC_Mesh mesh,              /**< input - mesh */
  FC_Variable coord_displ,   /**< input - coordinated displacements */
  FC_Variable *edge_lengths  /**< output - edge lengths variable handle */
) {
  FC_ReturnCode rc;
  int i, j, k;
  int dim, numEdge;
  double *lengths;
  FC_AssociationType assoc;
  FC_ElementType elemType;
  int *edgeConns;
  double* coords;
  FC_Coords points[2] = { { 0., 0., 0.}, { 0., 0., 0.} };
  
  // default return
  if (edge_lengths) 
    *edge_lengths = FC_NULL_VARIABLE;

  // check input
  if (!fc_isMeshValid(mesh) || 
      !fc_isValidDisplacementVariable(mesh, coord_displ) ||
      edge_lengths == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
    
  // log message
  fc_printfLogMessage("Getting the edge lengths of the displaced mesh.");  

  // setup & more checking
  rc = fc_getMeshInfo(mesh, NULL, &dim, NULL, NULL, &elemType);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_getMeshNumEdge(mesh, &numEdge);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_getMeshEdgeConnsPtr(mesh, &edgeConns);
  if (rc != FC_SUCCESS)
    return rc;
  if (numEdge < 1 && elemType != FC_ET_LINE) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  rc = fc_getDisplacedMeshCoords(mesh, coord_displ, &coords);
  if (rc != FC_SUCCESS) 
    return rc;
  
  // compute edge lengths
  lengths = (double *)malloc(sizeof(double) * numEdge);
  if (lengths == NULL ) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  for (i = 0; i < numEdge; i++) {
    for (j = 0; j < 2; j++) 
      for (k = 0; k < dim; k++)
	points[j][k] = coords[dim*edgeConns[i*2+j] + k];
    fc_calcEuclideanDistance(points[0], points[1], dim, &lengths[i]);
  }
  free(coords);

  // create new variable
  if (elemType == FC_ET_LINE)
    assoc = FC_AT_ELEMENT;
  else
    assoc = FC_AT_EDGE;
  rc = fc_createVariable(mesh, "displaced_edge_lengths", edge_lengths);
  if (rc != FC_SUCCESS) {
    free(lengths);
    return rc;
  }
  rc = fc_setVariableDataPtr(*edge_lengths, numEdge, 1, assoc, FC_MT_SCALAR, 
                             FC_DT_DOUBLE, lengths);
  // don't free lengths!
  
  return rc;
}

/**
 * \ingroup  GeometricRelations
 * \brief Return a variable with the face areas of the mesh
 *
 * \description
 *
 *    If the mesh is made of 2D elements (i.e. tris and quads), the
 *    area variable will be associated with the elements. Otherwise it
 *    will be associated with the faces of the mesh.
 *
 *    This routine will return with an error if the elements of the mesh are 
 *    not at least 2D (e.g. it is an error to ask for areas of an edge mesh). 
 *
 *    This routine assumes that faces are planar.
 *
 * \modifications  
 *    - 08/23/04 WSK, Created
 */
FC_ReturnCode fc_getFaceAreas(
  FC_Mesh mesh,            /**< input - mesh */
  FC_Variable *face_areas  /**< output - face areas variable handle */
) {
  FC_ReturnCode rc;
  int i, j, k;
  int dim, numFace;
  double *areas;
  FC_AssociationType assoc;
  FC_ElementType elemType, *faceTypes;
  FC_Coords coords[4];    // buffer for coords
  double* mesh_coords;    // coords for entire mesh
  int maxNumVertPerFace, *numVertPerFaces, *faceConns;
  
  // default return
  if (face_areas) 
    *face_areas = FC_NULL_VARIABLE;

  // check input
  if (!fc_isMeshValid(mesh) || face_areas == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
    
  // log message
  fc_printfLogMessage("Getting the face areas of the mesh.");  

  // setup & more checking
  rc = fc_getMeshInfo(mesh, NULL, &dim, NULL, NULL, &elemType);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_getMeshCoordsPtr(mesh, &mesh_coords);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_getMeshNumFace(mesh, &numFace);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_getMeshFaceConnsPtr(mesh, &numVertPerFaces, &maxNumVertPerFace,
			      &faceConns);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_getMeshFaceTypesPtr(mesh, &faceTypes);
  if (numFace == 0 && !(elemType == FC_ET_TRI || elemType == FC_ET_QUAD)) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // compute face areas
  areas = (double *)malloc(sizeof(double) * numFace);
  if (areas == NULL ) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  } 
  for (i = 0; i < numFace; i++) {
    for (j = 0; j < numVertPerFaces[i]; j++)
      for (k = 0; k < dim; k++) 
        coords[j][k] = mesh_coords[faceConns[i*maxNumVertPerFace+j]*dim+k];
    switch (faceTypes[i]) {
    case FC_ET_TRI:   _fc_calcTriArea(dim, coords, &areas[i]);   break;
    case FC_ET_QUAD:  _fc_calcQuadArea(dim, coords, &areas[i]);  break;
    default:
      fc_printfErrorMessage("Unexpected face type %s encountered",
                            fc_getElementTypeText(faceTypes[i]));
      return FC_ERROR;
    }
  }

  // create new variable
  if (elemType == FC_ET_TRI || elemType == FC_ET_QUAD)
    assoc = FC_AT_ELEMENT;
  else
    assoc = FC_AT_FACE;
  rc = fc_createVariable(mesh, "face_areas", face_areas);
  if (rc != FC_SUCCESS) {
    free(areas);
    return rc;
  }
  rc = fc_setVariableDataPtr(*face_areas, numFace, 1, assoc, FC_MT_SCALAR, 
                          FC_DT_DOUBLE, areas);
  // don't free areas! 
  
  return rc;
}

/**
 * \ingroup  GeometricRelations
 * \brief Return a variable with the face areas of the displaced mesh
 *
 * \description
 *
 *    If the mesh is made of 2D elements (i.e. tris and quads), the
 *    area variable will be associated with the elements. Otherwise it
 *    will be associated with the faces of the mesh.
 *
 *    This routine will return with an error if the elements of the mesh are 
 *    not at least 2D (e.g. it is an error to ask for areas of an edge mesh). 
 *
 *    This routine assumes that faces are planar.
 *
 * \modifications  
 *    - 11/1/04 WSK, Created
 */
FC_ReturnCode fc_getDisplacedFaceAreas(
  FC_Mesh mesh,             /**< input - mesh */
  FC_Variable coord_displ,  /**< input - coordinated displacements */
  FC_Variable *face_areas   /**< output - face areas variable handle */
) {
  FC_ReturnCode rc;
  int i, j, k;
  int dim, numFace;
  double *areas;
  FC_AssociationType assoc;
  FC_ElementType elemType, *faceTypes;
  double* coords;
  int maxNumVertPerFace, *numVertPerFaces, *faceConns;
  // buffer for coords of face
  FC_Coords points[4] = { { 0., 0., 0.}, { 0., 0., 0.},
			  { 0., 0., 0.}, { 0., 0., 0.} };  
  
  // default return
  if (face_areas) 
    *face_areas = FC_NULL_VARIABLE;

  // check input
  if (!fc_isMeshValid(mesh) || 
      !fc_isValidDisplacementVariable(mesh, coord_displ) ||
      face_areas == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
    
  // log message
  fc_printfLogMessage("Getting the face areas of the displaced mesh.");  

  // setup & more checking
  rc = fc_getMeshInfo(mesh, NULL, &dim, NULL, NULL, &elemType);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_getMeshNumFace(mesh, &numFace);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_getMeshFaceConnsPtr(mesh, &numVertPerFaces, &maxNumVertPerFace,
			      &faceConns);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_getMeshFaceTypesPtr(mesh, &faceTypes);
  if (rc != FC_SUCCESS)
    return rc;
  if (numFace == 0 && !(elemType == FC_ET_TRI || elemType == FC_ET_QUAD)) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  rc = fc_getDisplacedMeshCoords(mesh, coord_displ, &coords);
  if (rc != FC_SUCCESS)
    return rc;
  
  // compute face areas
  areas = (double *)malloc(sizeof(double) * numFace);
  if (areas == NULL ) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  } 
  for (i = 0; i < numFace; i++) {
    for (j = 0; j < numVertPerFaces[i]; j++)
      for (k = 0; k < dim; k++) 
        points[j][k] = coords[dim*faceConns[i*maxNumVertPerFace+j] + k];
    switch (faceTypes[i]) {
    case FC_ET_TRI:   _fc_calcTriArea(dim, points, &areas[i]);   break;
    case FC_ET_QUAD:  _fc_calcQuadArea(dim, points, &areas[i]);  break;
    default:
      fc_printfErrorMessage("Unexpected face type %s encountered",
                            fc_getElementTypeText(faceTypes[i]));
      return FC_ERROR;
    }
  }
  free(coords);

  // create new variable
  if (elemType == FC_ET_TRI || elemType == FC_ET_QUAD)
    assoc = FC_AT_ELEMENT;
  else
    assoc = FC_AT_FACE;
  rc = fc_createVariable(mesh, "displaced_face_areas", face_areas);
  if (rc != FC_SUCCESS) {
    free(areas);
    return rc;
  }
  rc = fc_setVariableDataPtr(*face_areas, numFace, 1, assoc, FC_MT_SCALAR, 
			     FC_DT_DOUBLE, areas);
  // don't free areas! 
  
  return rc;
}

/**
 * \ingroup  GeometricRelations
 * \brief Return a variable with the region volumes of the mesh.
 *
 * \description
 *
 *    This will return with an error if the elements of the mesh are not 3D 
 *    (e.g. it is an error to ask for volumes of a quad mesh).
 *
 *    This routine assumes that the faces of the regions are planar.
 *
 * \modifications  
 *    - 08/23/04 WSK, Created
 */
FC_ReturnCode fc_getRegionVolumes(
  FC_Mesh mesh,            /**< input - mesh */
  FC_Variable *region_volumes  /**< output - region volumes variable handle */
) {
  FC_ReturnCode rc;
  int i, j, k;
  int dim, numElement, numVertPerElem;
  double *vols;
  FC_ElementType elemType;
  FC_Coords points[8];    // buffer for coords
  double* coords;
  int* conns;

  // default return
  if (region_volumes) 
    *region_volumes = FC_NULL_VARIABLE;
  
  // check input
  if (!fc_isMeshValid(mesh) || region_volumes == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Getting the region volumes of the mesh.");  
  
  // setup & more checking
  rc = fc_getMeshInfo(mesh, NULL, &dim, NULL, &numElement, &elemType);
  if (rc != FC_SUCCESS)
    return rc;
  if (elemType != FC_ET_TET && elemType != FC_ET_PYRAMID && 
      elemType != FC_ET_PRISM && elemType != FC_ET_HEX) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  rc = fc_getMeshCoordsPtr(mesh, &coords);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_getMeshElementConnsPtr(mesh, &conns);
  if (rc != FC_SUCCESS)
    return rc;
  numVertPerElem = fc_getElementTypeNumVertex(elemType);
  
  // compute element volumes
  vols = (double *)malloc(sizeof(double) * numElement);
  if (vols == NULL ) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  } 
  for (i = 0; i < numElement; i++) {
    for (j = 0; j < numVertPerElem; j++)
      for (k = 0; k < dim; k++) 
        points[j][k] = coords[conns[i*numVertPerElem+j]*dim+k];
    switch (elemType) {
    case FC_ET_TET:      _fc_calcTetVolume(points, &vols[i]);      break;
    case FC_ET_PYRAMID:  _fc_calcPyramidVolume(points, &vols[i]);  break;
    case FC_ET_PRISM:    _fc_calcPrismVolume(points, &vols[i]);    break;
    case FC_ET_HEX:      _fc_calcHexVolume(points, &vols[i]);      break;
    default:
      fc_printfErrorMessage("Unexpected element type %s encountered",
                            fc_getElementTypeText(elemType));
      return FC_ERROR;
    }
  }

  // create new variable
  rc = fc_createVariable(mesh, "region_volumes", region_volumes);
  if (rc != FC_SUCCESS) {
    free(vols);
    return rc;
  }
  rc = fc_setVariableDataPtr(*region_volumes, numElement, 1, FC_AT_ELEMENT,
                             FC_MT_SCALAR, FC_DT_DOUBLE, vols);
  // don't free vols!
  
  return rc;
}

/**
 * \ingroup  GeometricRelations
 * \brief Return a variable with the region volumes of the displaced mesh.
 *
 * \description
 *
 *    This will return with an error if the elements of the mesh are not 3D 
 *    (e.g. it is an error to ask for volumes of a quad mesh).
 *
 *    This routine assumes that the faces of the regions are planar.
 *
 * \modifications  
 *    - 11/01/04 WSK, Created
 */
FC_ReturnCode fc_getDisplacedRegionVolumes(
  FC_Mesh mesh,             /**< input - mesh */
  FC_Variable coord_displ,  /**< input - coordinated displacements */
  FC_Variable *region_volumes  /**< output - region volumes variable handle */
) {
  FC_ReturnCode rc;
  int i, j, k;
  int dim, numElement, numVertPerElem;
  double *vols;
  FC_ElementType elemType;
  double* coords;
  int* conns;
  // buffer for coords of the region
  FC_Coords points[8] = { { 0., 0., 0.}, { 0., 0., 0.},
			  { 0., 0., 0.}, { 0., 0., 0.},
			  { 0., 0., 0.}, { 0., 0., 0.},
			  { 0., 0., 0.}, { 0., 0., 0.} };   

  // default return
  if (region_volumes) 
    *region_volumes = FC_NULL_VARIABLE;
  
  // check input
  if (!fc_isMeshValid(mesh) || 
      !fc_isValidDisplacementVariable(mesh, coord_displ) ||
      region_volumes == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Getting the region volumes of the displaced mesh.");  
  
  // setup & more checking
  rc = fc_getMeshInfo(mesh, NULL, &dim, NULL, &numElement, &elemType);
  if (rc != FC_SUCCESS)
    return rc;
  if (elemType != FC_ET_TET && elemType != FC_ET_PYRAMID && 
      elemType != FC_ET_PRISM && elemType != FC_ET_HEX) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  rc = fc_getDisplacedMeshCoords(mesh, coord_displ, &coords);
  if (rc != FC_SUCCESS)
    return rc;

  rc = fc_getMeshElementConnsPtr(mesh, &conns);
  if (rc != FC_SUCCESS)
    return rc;
  numVertPerElem = fc_getElementTypeNumVertex(elemType);
  
  // compute element volumes
  vols = (double *)malloc(sizeof(double) * numElement);
  if (vols == NULL ) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  } 
  for (i = 0; i < numElement; i++) {
    for (j = 0; j < numVertPerElem; j++)
      for (k = 0; k < dim; k++) 
        points[j][k] = coords[dim*conns[i*numVertPerElem+j] + k];
    switch (elemType) {
    case FC_ET_TET:      _fc_calcTetVolume(points, &vols[i]);      break;
    case FC_ET_PYRAMID:  _fc_calcPyramidVolume(points, &vols[i]);  break;
    case FC_ET_PRISM:    _fc_calcPrismVolume(points, &vols[i]);    break;
    case FC_ET_HEX:      _fc_calcHexVolume(points, &vols[i]);      break;
    default:
      fc_printfErrorMessage("Unexpected element type %s encountered",
                            fc_getElementTypeText(elemType));
      return FC_ERROR;
    }
  }
  free(coords);

  // create new variable
  rc = fc_createVariable(mesh, "region_volumes", region_volumes);
  if (rc != FC_SUCCESS) {
    free(vols);
    return rc;
  }
  rc = fc_setVariableDataPtr(*region_volumes, numElement, 1, FC_AT_ELEMENT,
                             FC_MT_SCALAR, FC_DT_DOUBLE, vols);
  // don't free vols!
  
  return rc;
}


/**
 * \ingroup  GeometricRelations
 * \brief Return a variable with the volume of the mesh.
 *
 * \description
 *
 *    This will return with an error and volume = -1 if the elements
 *    of the mesh are not 3D (e.g. it is an error to ask for volumes
 *    of a quad mesh).
 *
 *    This routine uses fc_getRegionVolumes and therefore has the same
 *    assumptions fc_getRegionVolumes does.
 *
 * \modifications  
 *    - 11/23/04 ACG, Created
 */
FC_ReturnCode fc_getMeshVolume(
  FC_Mesh mesh,            /**< input - mesh */
  double *volume  /**< output - volume of mesh */
) {

  FC_ReturnCode rc;
  FC_Variable elemvols;
  double vol;

  // log message
  fc_printfLogMessage("Calculating volume of mesh");


  // default return
  if (volume)
    *volume = -1;

  // test input
  if (!volume) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  //this will test if mesh is valid
  rc = fc_getRegionVolumes(mesh,&elemvols);
  if (rc != FC_SUCCESS){
    fc_deleteVariable(elemvols);
    return rc;
  }
  rc = fc_getVariableSum(elemvols,&vol);
  if (rc != FC_SUCCESS){
    fc_deleteVariable(elemvols);
    return rc;
  }

  fc_deleteVariable(elemvols);
  *volume = vol;
  return FC_SUCCESS;
}


/**
 * \ingroup  GeometricRelations
 * \brief Compute displaced mesh volume
 *
 * \description
 *
 *    Gets the volume of the mesh after displacing the coordinates
 *    with the given variable.
 *
 *    This will return with an error and vol = -1 if the elements of
 *    the mesh are not 3D (e.g. it is an error to ask for volumes of
 *    a quad mesh).
 *
 *    This routine uses fc_getDisplacedRegionVolumes and therefore has the same
 *    assumptions fc_getDisplacedRegionVolumes does.
 *
 *
 * \modifications   
 *   - 12/03/04 ACG, Created.
 */
FC_ReturnCode fc_getDisplacedMeshVolume(
  FC_Mesh mesh,            /**< input - mesh handle */
  FC_Variable coord_displ, /**< input - coord displacements */
  double *volume  /**< output - volume of mesh */
) {

  FC_ReturnCode rc;
  FC_Variable elemvols;
  double vol;

  // default return
  if (volume)
    *volume = -1;

  // test input
  if (!fc_isMeshValid(mesh) || 
      !fc_isValidDisplacementVariable(mesh, coord_displ) || !volume) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // log message
  fc_printfLogMessage("Calculating volume of displaced mesh");

  //this will test if mesh and the coord_displ are valid
  rc = fc_getDisplacedRegionVolumes(mesh,coord_displ,&elemvols);
  if (rc != FC_SUCCESS){
    fc_deleteVariable(elemvols);
    return rc;
  }
  rc = fc_getVariableSum(elemvols,&vol);
  if (rc != FC_SUCCESS){
    fc_deleteVariable(elemvols);
    return rc;
  }

  fc_deleteVariable(elemvols);
  *volume = vol;
  return FC_SUCCESS;
}



/**
 * \ingroup  GeometricRelations
 * \brief Returns a variable with the volume of the subset
 *
 * \description
 *
 *    This will return with an error if the elements of the mesh are not 3D 
 *    (e.g. it is an error to ask for volumes of a subset on a quad mesh).
 *
 *    This routine uses fc_getRegionVolumes and therefore has the same
 *    assumptions fc_getRegionVolumes does.
 *
 *    This routine uses fc_copySubsetWithNewAssociation to make
 *    a temporary intermediate subset with association FC_AT_ELEMENT.
 *    the "doStrict" input variable is for this call. See documentation
 *    there for further information.
 * 
 *
 * \modifications  
 *    - 11/30/04 ACG, Created
 */
FC_ReturnCode fc_getSubsetVolume(
  FC_Subset subs, /**< input - subset */
  int doStrict, /**< input - flag for dropping incomplete entities */
  double *volume  /**< output - volume of subset */
) {
  int k;
  FC_Mesh mesh;
  FC_ReturnCode rc;
  FC_Variable elemvols;
  FC_Subset  newSubset;
  void* actual_vols;
  int numDataPoint;
  char newSubsetName[20] = "junkSubset";
  int subsize;
  int* memarray;
  double vol = 0.0;

  // default return
  if (volume)
    *volume = -1;

  // test input
  if (!fc_isSubsetValid(subs) || !volume ||
      (doStrict != 0 && doStrict != 1)) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // log message
  fc_printfLogMessage("Computing subset volume.");

  //need parent 
  rc = fc_getMeshFromSubset(subs,&mesh);
  if (rc != FC_SUCCESS){
    //    fc_printfErrorMessage("cant get mesh from subset");
    return rc;
  }

  // get volumes of regions of mesh
  // know these have assoc fc_AT_ELEMENT
  rc = fc_getRegionVolumes(mesh,&elemvols);
  if (rc != FC_SUCCESS){
    //    fc_printfErrorMessage("cant get region volumes");
    fc_deleteVariable(elemvols);
    return rc;
  }

  rc = fc_getVariableDataPtr(elemvols,&actual_vols);
  if (rc != FC_SUCCESS){
    //    fc_printfErrorMessage("cant get actual vols");
    fc_deleteVariable(elemvols);
    return rc;
  }

  //just want to know how many there are for safety check below
  rc = fc_getVariableNumDataPoint(elemvols, &numDataPoint);
  if (rc != FC_SUCCESS){
    fc_deleteVariable(elemvols);
    return rc;
  }

  rc = fc_copySubsetWithNewAssociation
    (subs,newSubsetName,FC_AT_ELEMENT,doStrict,&newSubset);
  if (rc!= FC_SUCCESS){
    fc_deleteVariable(elemvols);
    return rc;
  }

  //next steps to match up subset elements with mesh region vols
  rc = fc_getSubsetMembersAsArray(newSubset,&subsize,&memarray);
  if (rc != FC_SUCCESS){
    //    fc_printfErrorMessage("cant get subset members as array");
    fc_deleteVariable(elemvols);
    rc = fc_deleteSubset(newSubset);
    free(memarray);
    return rc;
  }

  vol = 0.0;
  for (k = 0; k < subsize; k++)
    vol += ((double *)(actual_vols))[memarray[k]];

  fc_deleteVariable(elemvols);
  rc = fc_deleteSubset(newSubset);
  free(memarray);
  *volume = vol;
  return FC_SUCCESS;
}

/**
 * \ingroup  GeometricRelations
 * \brief Returns a variable with the volume of the subset on the
 *        displaced mesh
 *
 * \description
 *
 *    This will return with an error if the elements of the mesh are not 3D 
 *    (e.g. it is an error to ask for volumes of a subset on a quad mesh).
 *
 *    This routine uses fc_getDisplacedRegionVolumes and therefore has the same
 *    assumptions fc_getDisplacedRegionVolumes does.
 *
 *    This routine uses fc_copySubsetWithNewAssociation to make
 *    a temporary intermediate subset with association FC_AT_ELEMENT.
 *    the "doStrict" input variable is for this call. See documentation
 *    there for further information.
 * 
 *
 * \modifications  
 *    - 12/06/04 ACG, Created
 */
FC_ReturnCode fc_getDisplacedSubsetVolume(
  FC_Subset subs, /**< input - subset */
  int doStrict, /**< input - flag for dropping incomplete entities */
  FC_Variable coord_displ, /**< input - coord displacements */
  double *volume  /**< output - volume of subset */
) {
  int k;
  FC_Mesh mesh;
  FC_ReturnCode rc;
  FC_Variable elemvols;
  FC_Subset  newSubset;
  void* actual_vols;
  int numDataPoint;
  char newSubsetName[20] = "junkSubset";
  int subsize;
  int* memarray;
  double vol = 0.0;

  // default return
  if (volume)
    *volume = -1;

  // test input
  if (!fc_isSubsetValid(subs) || (doStrict != 0 && doStrict != 1) ||
      !volume) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  rc = fc_getMeshFromSubset(subs, &mesh);
  if (rc != FC_SUCCESS)
    return rc;
  if (!fc_isValidDisplacementVariable(mesh, coord_displ)) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // log message
  fc_printfLogMessage("Computing displaced subset volume.");

  // get volumes of displaced regions of mesh
  rc = fc_getDisplacedRegionVolumes(mesh, coord_displ, &elemvols);
  if (rc != FC_SUCCESS) {
    fc_deleteVariable(elemvols);
    return rc;
  }

  rc = fc_getVariableDataPtr(elemvols,&actual_vols);
  if (rc != FC_SUCCESS) {
    //    fc_printfErrorMessage("cant get actual vols");
    fc_deleteVariable(elemvols);
    return rc;
  }

  // just want to know how many there are for safety check below
  rc = fc_getVariableNumDataPoint(elemvols, &numDataPoint);
  if (rc != FC_SUCCESS) {
    fc_deleteVariable(elemvols);
    return rc;
  }

  rc = fc_copySubsetWithNewAssociation(subs, newSubsetName, FC_AT_ELEMENT, 
				       doStrict, &newSubset);
  if (rc!= FC_SUCCESS) {
    fc_deleteVariable(elemvols);
    return rc;
  }

  // next steps to match up subset elements with displaced mesh region vols
  rc = fc_getSubsetMembersAsArray(newSubset, &subsize, &memarray);
  if (rc != FC_SUCCESS){
    //    fc_printfErrorMessage("cant get subset members as array");
    fc_deleteVariable(elemvols);
    rc = fc_deleteSubset(newSubset);
    free(memarray);
    return rc;
  }

  vol = 0.0;
  for (k = 0; k < subsize; k++)
    vol += ((double *)(actual_vols))[memarray[k]];

  fc_deleteVariable(elemvols);
  rc = fc_deleteSubset(newSubset);
  free(memarray);
  *volume = vol;
  return FC_SUCCESS;
}


/**
 * \ingroup  GeometricRelations
 * \brief Returns a variable with the area of the subset
 *
 * \description
 *
 *    Returns the area of subset. The subset must have assoc
 *    FC_AT_FACE. Fails for NULL Subset, Returns 0 for empty subset.
 *
 * \todo
 *    - make displaced version
 *
 * \modifications  
 *    - 07/18/06 ACG, Created from old fc_getSidesAreas
 */
FC_ReturnCode fc_getSubsetArea(
  FC_Subset subset, /**< input - subset */
  double *area  /**< output - area of subset */
) {
  FC_Mesh mesh;
  FC_ReturnCode rc;
  FC_Variable faceareas;
  double *areaarray;
  double retarea = 0.0;
  FC_AssociationType assoc;
  int numMem, *ids;
  int i;

  //default return
  if (area)
    *area = -1;

  if (!fc_isSubsetValid(subset) || !area ){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }


  rc = fc_getSubsetAssociationType(subset, &assoc);
  if (rc != FC_SUCCESS)
    return rc;
  if ( assoc != FC_AT_FACE){
    fc_printfErrorMessage("Wrong association type");
    return FC_INPUT_ERROR;
  }

  // log message
  fc_printfLogMessage("Computing Subset Area");


  rc = fc_getSubsetMembersAsArray(subset, &numMem, &ids);
  if (rc!= FC_SUCCESS){
    return rc;
  }
  if (numMem == 0){
    *area = 0.0;
    return FC_SUCCESS;
  }

  //get areas for all faces of the mesh
  rc = fc_getMeshFromSubset(subset,&mesh);
  if (rc!= FC_SUCCESS){
    free(ids);
    return rc;
  }

  rc = fc_getFaceAreas(mesh,&faceareas);
  if (rc!= FC_SUCCESS){
    fc_printfErrorMessage("cant get face areas");
    free(ids);
    return rc;
  }

  rc = fc_getVariableDataPtr(faceareas,(void **) &areaarray);
  if (rc!= FC_SUCCESS){
    free(ids);
    return rc;
  }


  for (i = 0; i < numMem; i++){
    retarea += areaarray[ids[i]];
  }
  free(ids);
  fc_deleteVariable(faceareas);
  *area = retarea;

  return FC_SUCCESS;
}




//@}

/** \name Normals */
//-------------------------------
//@{

/**
 * \ingroup  GeometricRelations
 * \brief Compute surface normals on all elements in a mesh
 *
 * \description
 *  
 *    Compute surface normals on all elements in a mesh. Currently
 *    only correct for meshes with planar elements (like triangles
 *    and planar quads).
 *
 * \modifications  
 *    - Steph Boyles  Created.
 *    - 2002-JUL-09  W Koegler  Adapted to modified struct Vector
 *    - 2003-SEP-04  W Koegler  Created: _fc_SurfaceNormal was split between
 *    this function and _fc_calcSurfaceNormal for better modularity.
 *
 * \todo needs a flag if the normals are already computed
 * \todo expand to 3D elements -> could return surface normals of faces?
 *      ... but what would their orientation be?
 */
FC_ReturnCode fc_getSurfaceNormals(
  FC_Mesh mesh,          /**< input - mesh to compute surface normals of */
  FC_Variable *new_var   /**< output - new variable w/ surface normals as data */
) {
  FC_ReturnCode rc;
  int i, j, k;
  double *normals;
  FC_Vector temp_normal;  // normal vector to an element
  FC_Coords temp_vertices[4]; // array of vertex coords to pass around
  int numElem, dim, numVertPerElem;
  FC_ElementType elemType;
  double* coords;
  int* conns;
  
  // default return
  if (new_var)
    *new_var = FC_NULL_VARIABLE;
  
  // check input
  if (!fc_isMeshValid(mesh) || !new_var) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // setup and more checking
  rc = fc_getMeshInfo(mesh, NULL, &dim, NULL, &numElem, &elemType);
  if (rc != FC_SUCCESS)
    return rc;
  if (elemType != FC_ET_TRI && elemType != FC_ET_QUAD) {
    fc_printfErrorMessage("%s, can't handle element type %s", 
        fc_getReturnCodeText(FC_INPUT_ERROR),fc_getElementTypeText(elemType));
    return FC_INPUT_ERROR;
  }
  rc = fc_getMeshCoordsPtr(mesh, &coords);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_getMeshElementConnsPtr(mesh, &conns);
  if (rc != FC_SUCCESS)
    return rc;
  numVertPerElem = fc_getElementTypeNumVertex(elemType);
  
  // log message
  fc_printfLogMessage("Computing surface normals.");
  
  // make space for the data array 
  normals = (double *)malloc(sizeof(double) * dim * numElem);
  if (normals == NULL ) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  } 
  
  // compute Surface normals, save in buffer
  for (i = 0; i < numElem; i++) {

    // make list of pointers to this elements vertices
    for (j = 0; j < numVertPerElem; j++) {
      for (k = 0; k < dim; k++)
        temp_vertices[j][k] = coords[conns[i*numVertPerElem+j]*dim+k];
      for (k = dim; k < 3; k++)
        temp_vertices[j][k] = 0;
    }

    // calculate normal
    _fc_calcSurfaceNormal(numVertPerElem, temp_vertices, &temp_normal);
    
    // save results in element info and an array to make into a var
    for (j = 0; j < dim; j++) 
      normals[i*dim+j] = temp_normal[j];
  }
  
  // add new variable to current mesh
  fc_createVariable(mesh, "surfaceNormals", new_var);
  fc_setVariableDataPtr(*new_var, numElem, dim, FC_AT_ELEMENT, FC_MT_VECTOR, 
                        FC_DT_DOUBLE, normals);
  // don't free normals
  
  return FC_SUCCESS;
}

/**
 * \ingroup  GeometricRelations
 * \brief Compute surface normals on all elements in a mesh
 *
 * \description
 *
 *    Compute surface normals on all elements in a mesh after the mesh has been
 *    deformed by adding a displacement vector to the coordinates. See
 *    fc_getSurfaceNormals.
 *
 * \modifications 
 *   - 2003-JUL-01  W Koegler  Created from fc_getSurfaceNormals 
 *   - 2003-SEP-04  W Koegler  Moved calculation of displaced coords to its
 *   own function, fc_getDisplacedVertices. 
 */
FC_ReturnCode fc_getDisplacedSurfaceNormals(
  FC_Mesh mesh,          /**< input - mesh to compute surface normals of */
  FC_Variable displ,     /**< input - coord displacment variable */
  FC_Variable *new_var   /**< output - new variable w/ surface normals as data */
) 
{
  FC_ReturnCode rc;
  int i, j, k;
  double *normals;
  FC_Vector temp_normal;         // normal vector an the element
  int numElem, dim, numVertPerElem;
  FC_ElementType elemType;
  double *displ_coords; // mesh's displaced coords
  int* conns;
  FC_Coords temp_vertices[4];

  // default return
  if (new_var)
    *new_var = FC_NULL_VARIABLE;
  
  // check input
  if (!fc_isMeshValid(mesh) || 
      !fc_isValidDisplacementVariable(mesh, displ) || 
      !new_var) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // log message
  fc_printfLogMessage("Getting displaced surface normals.");  
  
  // get this meshes element info
  rc = fc_getMeshInfo(mesh, NULL, &dim, NULL, &numElem, &elemType);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_getMeshElementConnsPtr(mesh, &conns);
  if (rc != FC_SUCCESS)
    return rc;
  numVertPerElem = fc_getElementTypeNumVertex(elemType);
 
  // get the displaced vertex list
  fc_getDisplacedMeshCoords(mesh, displ, &displ_coords);
  
  // make space for the data array 
  normals = (double *)malloc(sizeof(double) * dim * numElem);
  if (normals == NULL ) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  } 

  // compute Surface normals, save in buffer
  for (i = 0; i < numElem; i++) {
    // make list of pointers to this elements vertices
    for (j = 0; j < numVertPerElem; j++) {
      for (k = 0; k < dim; k++)
        temp_vertices[j][k] = displ_coords[conns[i*numVertPerElem+j]*dim + k];
      for (k = dim; k < 3; k++)
        temp_vertices[j][k] = 0.;
    }

    // calculate normal  
    _fc_calcSurfaceNormal(numVertPerElem, temp_vertices, &temp_normal);
    
    // save in array to make var
    for (j = 0; j < dim; j++) 
      normals[i*dim+j] = temp_normal[j];
  }
  free(displ_coords);
  
  // add new variable to current mesh
  fc_createVariable(mesh, "displacedSurfaceNormals", new_var);
  fc_setVariableDataPtr(*new_var, numElem, dim, FC_AT_ELEMENT, FC_MT_VECTOR, 
                        FC_DT_DOUBLE, normals);
  // don't free normals
 
  return FC_SUCCESS;
}

//@}

/** \name Deformation */
//-------------------------------
//@{

/** 
 * \ingroup GeometricRelations
 * \brief Calculation the deformation of a mesh.
 *
 * \description
 *
 *   The displacement variable that is used to define changing mesh geometry 
 *    is a combination of mesh translation, rotation, and deformation.  This 
 *   routine attempts to extract the deformation component.  It assumes that 
 *   deformation is zero in the first step.
 *
 *   This routines assumes that there is no rotation, and then uses the
 *   supplied reference vertex to calculate the translation component. The 
 *   deformation component is then the displacement minus the translation.
 *   This routine will not perform well if the reference vertex is in a region
 *   that deforms.
 *
 * \todo Want to add two more routines: 1) takes 3 points & so might handles
 *   rotation, 2) use optimization to find transformation matrix that
 *   minimizes deformation - most robust and automatic.
 *
 * \modifications
 *   - 04/19/05 WSK Created.
 */
FC_ReturnCode fc_calcMeshDeformationUsingRefVertex(
  FC_Mesh mesh,        /**< input - mesh handle */
  int vertexID,        /**< input - vertex to use as translation reference */
  int numStep,  /**< input - number of steps in displacement variable */
  FC_Variable* displ,  /**< input - the displacement variable */
  FC_Variable** deform /**< output - the deformation variable */
) {
  FC_ReturnCode rc;
  int i, j, k;
  int numDim, numVertex, temp_numStep;
  FC_Sequence sequence;
  double* displ_coords, *current, orig[3];
  double* displ_data, *deform_data;
  
  // default output
  if (deform)
    *deform = NULL;

  // check input
  if (!fc_isMeshValid(mesh) || vertexID < 0 || 
      !fc_isSeqVariableValid(numStep, displ) ||
      !fc_isValidDisplacementVariable(mesh, displ[0])) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // setup
  rc = fc_getVariableInfo(displ[0], &numVertex, &numDim, NULL, NULL, NULL);
  if (rc != FC_SUCCESS)
    return rc;
  
  // more checking
  if (vertexID >= numVertex) {
    fc_printfErrorMessage("The vertex ID is out of range");
    return FC_INPUT_ERROR;
  }
  
  // create the deformation variable
  rc = fc_getSequenceFromSeqVariable(numStep, displ, &sequence);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_createSeqVariable(mesh, sequence, "deformation", &temp_numStep,
			    deform);
  if (rc != FC_SUCCESS) {
    fc_printfErrorMessage("Could not create sequence variable");
    return rc;
  }

  // Set the data
  for (i = 0; i < numStep; i++) {
    // get the displacement of the vertex
    // (might be slightly more efficient to calc this on our own)
    rc = fc_getDisplacedMeshCoords(mesh, displ[i], &displ_coords);
    if (rc != FC_SUCCESS)
      return rc;
    current = &displ_coords[vertexID*numDim];
    if (i == 0) 
      for (j = 0; j < numDim; j++)
	orig[j] = current[j];
    
    // create and set deform data
    rc = fc_getVariableDataAsDataType(displ[i], FC_DT_DOUBLE,
				      (void**)&displ_data);
    if (rc != FC_SUCCESS) {
      fc_printfErrorMessage("Could not get displ data as double");
      return rc;
    }
    deform_data = (double*)malloc(numDim*numVertex*sizeof(double));
    if (!deform_data) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      return FC_MEMORY_ERROR;
    }
    // IMPORTANT! To prevent bad roundoff errors, check that values are 
    // different before subtracting. If they are not different, answer is 0!
    // (Using FC_VALUE_EQUIV instead of FC_DBL_EQUIV because of roundoff in 
    // centroid calculation, have to be more flexible)
    for (j = 0; j < numVertex; j++) {
      for (k = 0; k < numDim; k++) {
	if (FC_VALUE_EQUIV(displ_data[j*numDim+k] + orig[k], current[k],
			   100*DBL_EPSILON, DBL_MIN))
	  deform_data[j*numDim+k] = 0.;
	else
	  deform_data[j*numDim+k] = displ_data[j*numDim+k] + orig[k] - current[k];
      }
    }
    free(displ_coords);
    rc = fc_setVariableDataPtr((*deform)[i], numVertex, numDim, FC_AT_VERTEX, 
			       FC_MT_VECTOR, FC_DT_DOUBLE, deform_data);
    free(displ_data);
  }
  
  return FC_SUCCESS;
}

/** 
 * \ingroup GeometricRelations
 * \brief Calculation the deformation of a mesh.
 *
 * \description
 *
 *   The displacement variable that is used to define changing mesh geometry 
 *   is a combination of mesh translation, rotation, and deformation. This 
 *   routine attempts to extract the deformation component.  It assumes that
 *   deformation is zero in the first step.
 *
 *   This routines assumes that there is no rotation, and then uses the
 *   centroid of the mesh to calculate the translation component. The 
 *   deformation component is then the displacement minus the translation.
 *   This routine will not perform well if deformation noticeably affects
 *   the centroid of the mesh.
 *
 *  \modifications
 *    - 04/19/05 WSK Created.
 */
FC_ReturnCode fc_calcMeshDeformationUsingCentroid(
  FC_Mesh mesh,        /**< input - mesh handle */
  int numStep,  /**< input - number of steps in displacement variable */
  FC_Variable* displ,  /**< input - the displacement variable */
  FC_Variable** deform /**< output - the deformation variable */
) {
  FC_ReturnCode rc;
  int i, j, k;
  int numDim, numVertex, temp_numStep;
  FC_Sequence sequence;
  double* displ_data, *deform_data;
  FC_Coords orig, current;
  
  // default output
  if (deform)
    *deform = NULL;

  // check input
  if (!fc_isMeshValid(mesh) || !fc_isSeqVariableValid(numStep, displ) ||
      !fc_isValidDisplacementVariable(mesh, displ[0])) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // setup
  rc = fc_getVariableInfo(displ[0], &numVertex, &numDim, NULL, NULL, NULL);
  if (rc != FC_SUCCESS)
    return rc;
  
  // create the deformation variable
  rc = fc_getSequenceFromSeqVariable(numStep, displ, &sequence);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_createSeqVariable(mesh, sequence, "deformation", &temp_numStep,
			    deform);
  if (rc != FC_SUCCESS) {
    fc_printfErrorMessage("Could not create sequence variable");
    return rc;
  }

  // Set the data
  for (i = 0; i < numStep; i++) {
    // get the displacement of the centroid
    rc = fc_getDisplacedMeshCentroid(mesh, displ[i], &numDim, &current);
    if (rc != FC_SUCCESS) {
      fc_printfErrorMessage("Could not get centroid");
      return rc;
    }
    if (i == 0) {
      for (j = 0; j < numDim; j++) 
	orig[j] = current[j];
    }

    // create and set deform data
    rc = fc_getVariableDataAsDataType(displ[i], FC_DT_DOUBLE,
				      (void**)&displ_data);
    if (rc != FC_SUCCESS) {
      fc_printfErrorMessage("Could not get displ data as double");
      return rc;
    }
    deform_data = (double*)malloc(numDim*numVertex*sizeof(double));
    if (!deform_data) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      return FC_MEMORY_ERROR;
    }
    // IMPORTANT! To prevent bad roundoff errors, check that values are 
    // different before subtracting. If they are not different, answer is 0!
    // (Using FC_VALUE_EQUIV instead of FC_DBL_EQUIV because of roundoff in 
    // centroid calculation, have to be more flexible)
    for (j = 0; j < numVertex; j++) {
      for (k = 0; k < numDim; k++) {
	if (FC_VALUE_EQUIV(displ_data[j*numDim+k] + orig[k], current[k],
			   100*DBL_EPSILON, DBL_MIN))
	  deform_data[j*numDim+k] = 0.;
	else
	  deform_data[j*numDim+k] = displ_data[j*numDim+k] + orig[k] - current[k];
      }
    }
    rc = fc_setVariableDataPtr((*deform)[i], numVertex, numDim, FC_AT_VERTEX, 
			       FC_MT_VECTOR, FC_DT_DOUBLE, deform_data);
    free(displ_data);
  }

  return FC_SUCCESS;
}

//@}

/** 
 *\name Proximity
 *
 * Use these routines to test if meshes and subsets are close to each other.
 * The actual test is to find pairs of vertices (one from each mesh/subset)
 * that are less than or equal to a specified distance cutoff.
 */
//-------------------------------
//@{

/** 
 * \ingroup GeometricRelations
 * \brief Return the proximal vertices between two meshes.
 *
 * \description
 *   
 *    This returns all vertex pairs that are less than or equal
 *    to min_dist apart. Vertices may show up in more than one pair.
 *    (In the extrememe case if you put in min_dist = infinity, you
 *    would end up with each vertex being connected to all vertices
 *    of the other mesh.)
 *
 *    The routine can be forced to return "early", after the first
 *    match is found, if you use NULL for the last two arguments (the
 *    vertex arrays). In this case, the returned numPair is
 *    really a flag for whether a pair was found or not and will
 *    never be greater than 1 even if there are more than one
 *    possible pairs.
 *
 *    Note, passing in min_dist = 0 will find overlapping vertices.
 *
 *    Note, meshes may be closer together than any of the vertex
 *    pairs (imagine 1 vertex of one mesh very close to middle of
 *    a face of another mesh. the vert to face distance is less
 *    than distance to the face's verts).
 *
 *    Will return if error if mesh1 == mesh2.
 *
 * \modifications
 *    - 02/13/2006 WSD Created.
 *    - 03/14/2006 WSD guts common to all proximity routines to a common
 *         core routine.
 *    - 03/16/2006 WSD changed behavior so that if NULLs are passed
 *         in for the vertIDs, will return after 1st match.
 */
FC_ReturnCode fc_getMeshesProximity(
  FC_Mesh mesh1, /**< input - the first mesh */
  FC_Mesh mesh2, /**< input - the second mesh */
  double min_dist, /**< input - the distance cutoff */
  int* numPair_p, /**< output - the number of vertex pairs found */
  int** mesh1VertIDs_p, /**< output - ids of the mesh1 verts (order matters),
			 user is responsible for freeing */
  int** mesh2VertIDs_p /**< output - ids of the mesh2 verts (order matters),
			 user is responsible for freeing */
) {
  FC_ReturnCode rc;
  FC_Subset subset1, subset2;

  // argument testing & default returns
  rc = _fc_getSubsetsProximityArgsTest(&mesh1, NULL, NULL, &mesh2,
				       NULL, NULL, min_dist, numPair_p,
				       mesh1VertIDs_p, mesh2VertIDs_p);
  if (rc != FC_SUCCESS)
    return rc;

  // log message
  fc_printfLogMessage("Computing meshes proximity");

  // convert to subsets (WHOLE_MESH)
  rc = fc_createSubset(mesh1, "temp", FC_AT_WHOLE_MESH, &subset1);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_addMemberToSubset(subset1, 0);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_createSubset(mesh2, "temp", FC_AT_WHOLE_MESH, &subset2);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_addMemberToSubset(subset2, 0);
  if (rc != FC_SUCCESS)
    return rc;

  // do it
  rc = _fc_getSubsetsProximityCore(subset1, NULL, subset2, NULL, min_dist, 
				   numPair_p, mesh1VertIDs_p, mesh2VertIDs_p);

  // cleanup
  fc_deleteSubset(subset1);
  fc_deleteSubset(subset2);

  return rc;
}

/** 
 * \ingroup GeometricRelations
 * \brief Return the proximal vertices between two displaced meshes.
 *
 * \description
 *   
 *    This returns all vertex pairs that are less than or equal
 *    to min_dist apart. Vertices may show up in more than one pair.
 *    (In the extrememe case if you put in min_dist = infinity, you
 *    would end up with each vertex being connected to all vertices
 *    of the other mesh.)
 *
 *    The routine can be forced to return "early", after the first
 *    match is found, if you use NULL for the last two arguments (the
 *    vertex arrays). In this case, the returned numPair is
 *    really a flag for whether a pair was found or not and will
 *    never be greater than 1 even if there are more than one
 *    possible pairs.
 *
 *    Note, passing in min_dist = 0 will find overlapping vertices.
 *
 *    Note, meshes may be closer together than any of the vertex
 *    pairs (imagine 1 vertex of one mesh very close to middle of
 *    a face of another mesh. the vert to face distance is less
 *    than distance to the face's verts).
 *
 *    Will return if error if mesh1 == mesh2 AND displ1 == displ2.
 *    (mesh1 can be the same as mesh2 if displ1 != displ2)
 *
 * \modifications
 *    - 02/14/2006 WSD Created.
 *    - 03/14/2006 WSD guts common to all proximity routines to a common
 *         core routine.
 *    - 03/16/2006 WSD changed behavior so that if NULLs are passed
 *         in for the vertIDs, will return after 1st match.
 */
FC_ReturnCode fc_getDisplacedMeshesProximity(
  FC_Mesh mesh1, /**< input - the first mesh */
  FC_Variable displ1, /**< input - the displacements on the first mesh */
  FC_Mesh mesh2, /**< input - the second mesh */
  FC_Variable displ2, /**< input - the displacements on the second mesh */
  double min_dist, /**< input - the distance cutoff */
  int* numPair_p, /**< output - the number of vertex pairs found */
  int** mesh1VertIDs_p, /**< output - ids of the mesh1 verts (order matters),
			 user is responsible for freeing */
  int** mesh2VertIDs_p /**< output - ids of the mesh2 verts (order matters),
			 user is responsible for freeing */
) {
  FC_ReturnCode rc;
  FC_Subset subset1, subset2;

  // argument testing & default returns
  rc = _fc_getSubsetsProximityArgsTest(&mesh1, NULL, &displ1, &mesh2,
				       NULL, &displ2, min_dist, numPair_p,
				       mesh1VertIDs_p, mesh2VertIDs_p);
  if (rc != FC_SUCCESS)
    return rc;

  // log message
  fc_printfLogMessage("Computing displaced meshes proximity");

  // convert to subsets that are whole
  rc = fc_createSubset(mesh1, "temp", FC_AT_WHOLE_MESH, &subset1);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_addMemberToSubset(subset1, 0);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_createSubset(mesh2, "temp", FC_AT_WHOLE_MESH, &subset2);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_addMemberToSubset(subset2, 0);
  if (rc != FC_SUCCESS)
    return rc;

  // do it
  rc = _fc_getSubsetsProximityCore(subset1, &displ1, subset2, &displ2,
				   min_dist, numPair_p, mesh1VertIDs_p,
				   mesh2VertIDs_p);

  // cleanup
  fc_deleteSubset(subset1);
  fc_deleteSubset(subset2);

  return rc;
}

/** 
 * \ingroup GeometricRelations
 * \brief Return the proximal vertices between two subsets.
 *
 * \description
 *   
 *    This returns all vertex pairs that are less than or equal
 *    to min_dist apart. Vertices may show up in more than one pair.
 *    (In the extrememe case if you put in min_dist = infinity, you
 *    would end up with each vertex being connected to all vertices
 *    of the other subset.)
 *
 *    The routine can be forced to return "early", after the first
 *    match is found, if you use NULL for the last two arguments (the
 *    vertex arrays). In this case, the returned numPair is
 *    really a flag for whether a pair was found or not and will
 *    never be greater than 1 even if there are more than one
 *    possible pairs.
 *
 *    Note, passing in min_dist = 0 will find overlapping vertices,
 *    or if subsets are on the same mesh, will find common vertices.
 *
 *    Note, subsets may be closer together than any of the vertex
 *    pairs (imagine 1 vertex of one mesh very close to middle of
 *    a face of another mesh. the vert to face distance is less
 *    than distance to the face's verts).
 *
 *    Will return if error if subset1 == subset2.
 *
 * \modifications
 *    - 02/14/2006 WSD Created.
 *    - 03/14/2006 WSD guts common to all proximity routines to a common
 *         core routine.
 *    - 03/16/2006 WSD changed behavior so that if NULLs are passed
 *         in for the vertIDs, will return after 1st match.
 */
FC_ReturnCode fc_getSubsetsProximity(
  FC_Subset subset1, /**< input - the first subset */
  FC_Subset subset2, /**< input - the second subset */
  double min_dist, /**< input - the distance cutoff */
  int* numPair_p, /**< output - the number of vertex pairs found */
  int** mesh1VertIDs_p, /**< output - ids of the mesh1 verts (order matters),
			 user is responsible for freeing */
  int** mesh2VertIDs_p /**< output - ids of the mesh2 verts (order matters),
			 user is responsible for freeing */
) {
  FC_ReturnCode rc;

  // argument testing & default returns
  rc = _fc_getSubsetsProximityArgsTest(NULL, &subset1, NULL, NULL,
				       &subset2, NULL, min_dist, numPair_p,
				       mesh1VertIDs_p, mesh2VertIDs_p);
  if (rc != FC_SUCCESS)
    return rc;
 
  // log message
  fc_printfLogMessage("Computing subsets proximity");

  // do it
  return _fc_getSubsetsProximityCore(subset1, NULL, subset2, NULL, min_dist,
				     numPair_p, mesh1VertIDs_p,
				     mesh2VertIDs_p);
}

/** 
 * \ingroup GeometricRelations
 * \brief Return the proximal vertices between two displaced subsets.
 *
 * \description
 *   
 *    This returns all vertex pairs that are less than or equal
 *    to min_dist apart. Vertices may show up in more than one pair.
 *    (In the extrememe case if you put in min_dist = infinity, you
 *    would end up with each vertex being connected to all vertices
 *    of the other subset.)
 *
 *    The routine can be forced to return "early", after the first
 *    match is found, if you use NULL for the last two arguments (the
 *    vertex arrays). In this case, the returned numPair is
 *    really a flag for whether a pair was found or not and will
 *    never be greater than 1 even if there are more than one
 *    possible pairs.
 *
 *    Note, passing in min_dist = 0 will find overlapping vertices,
 *    or if subsets are on the same mesh, will find common vertices.
 *
 *    Note, subsets may be closer together than any of the vertex
 *    pairs (imagine 1 vertex of one mesh very close to middle of
 *    a face of another mesh. the vert to face distance is less
 *    than distance to the face's verts).
 *
 *    Will return if error if subset1 == subset2 AND displ1 == displ2.
 *    (subset1 can be the same as subset2 if displ1 != displ2)
 *
 * \modifications
 *    - 02/14/2006 WSD Created.
 *    - 03/14/2006 WSD guts common to all proximity routines to a common
 *         core routine.
 *    - 03/16/2006 WSD changed behavior so that if NULLs are passed
 *         in for the vertIDs, will return after 1st match.
 */
FC_ReturnCode fc_getDisplacedSubsetsProximity(
  FC_Subset subset1, /**< input - the first subset */
  FC_Variable displ1, /**< input - the displacements on the mesh of the
		       first subset */
  FC_Subset subset2, /**< input - the second subset */
  FC_Variable displ2, /**< input - the displacements on the mesh of the
		       second subset */
  double min_dist, /**< input - the distance cutoff */
  int* numPair_p, /**< output - the number of vertex pairs found */
  int** mesh1VertIDs_p, /**< output - ids of the mesh1 verts (order matters),
			 user is responsible for freeing */
  int** mesh2VertIDs_p /**< output - ids of the mesh2 verts (order matters),
			 user is responsible for freeing */
) {
  FC_ReturnCode rc;

  // argument testing & default returns
  rc = _fc_getSubsetsProximityArgsTest(NULL, &subset1, &displ1, NULL,
				       &subset2, &displ2, min_dist, numPair_p,
				       mesh1VertIDs_p, mesh2VertIDs_p);
  if (rc != FC_SUCCESS)
    return rc;

  // log message
  fc_printfLogMessage("Computing displaced subsets proximity");

  return _fc_getSubsetsProximityCore(subset1, &displ1, subset2, &displ2,
				     min_dist, numPair_p,
				     mesh1VertIDs_p, mesh2VertIDs_p);
}

//@}

/** \name Proximity helpers */
//-------------------------------
//@{

/** 
 * \ingroup PrivateGeometricRelations
 * \brief Common argument testing for all proximity routines.
 *
 * \description
 *
 * \modifications
 *    - 03/14/2006 WSD Created.
 *    - 03/16/2006 WSD changed behavior so that if NULLs are passed
 *         in for the vertIDs, will return after 1st match.
 */
FC_ReturnCode _fc_getSubsetsProximityArgsTest(
  FC_Mesh* mesh1,      /**< input - the first mesh */
  FC_Subset* subset1,  /**< input - the first subset */
  FC_Variable* displ1, /**< input - the displs on the first mesh/subset */ 
  FC_Mesh* mesh2,      /**< input - the first mesh */
  FC_Subset* subset2,  /**< input - the second subset */
  FC_Variable* displ2, /**< input - the disps on the 2nd mesh/subset */
  double min_dist,     /**< input - the distance cutoff */
  int* numPair,        /**< output - the number of vertex pairs found */
  int** mesh1VertIDs,  /**< output - ids of the mesh1 verts (order matters),
			  user is responsible for freeing */
  int** mesh2VertIDs  /**< output - ids of the mesh2 verts (order matters),
			 user is responsible for freeing */
) {
  int numDim1, numDim2;
  FC_Mesh temp_mesh1, temp_mesh2;

  // default returns
  if (numPair)
    *numPair = -1;
  if (mesh1VertIDs)
    *mesh1VertIDs = NULL;
  if (mesh2VertIDs)
    *mesh2VertIDs = NULL;

  // test input
  // It is o.k. to pass in mesh1VertIDs = NULL && mesh2VertIDs = NULL
  if (!numPair || (!mesh1VertIDs && mesh2VertIDs) ||
      (mesh1VertIDs && !mesh2VertIDs)) { 
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  if ( (mesh1 && !fc_isMeshValid(*mesh1)) ||
       (mesh2 && !fc_isMeshValid(*mesh2)) ||
       (subset1 && !fc_isSubsetValid(*subset1)) ||
       (subset2 && !fc_isSubsetValid(*subset2)) ||
       (displ1 && !fc_isVariableValid(*displ1)) ||
       (displ2 && !fc_isVariableValid(*displ2)) ||
       fc_lt(min_dist, 0., DBL_EPSILON, DBL_MIN) ) { 
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  // If have 1, have to have the other, and can't be equal
  // (NOTE displ1 == displ2 implies mesh1 = mesh2, so not testing that
  if ( (mesh1 && !mesh2) || (!mesh1 && mesh2) ||
       (subset1 && !subset2) || (!subset1 && subset2) ||
       (displ1 && !displ2) || (!displ1 && displ2) || 
       (mesh1 && !displ1 && FC_HANDLE_EQUIV(*mesh1, *mesh2)) ||
       (mesh1 && displ1 && FC_HANDLE_EQUIV(*displ1, *displ2)) || 
       (subset1 && !displ1 && FC_HANDLE_EQUIV(*subset1, *subset2)) || 
       (subset1 && displ1 && (FC_HANDLE_EQUIV(*subset1, *subset2) &&
			      FC_HANDLE_EQUIV(*displ1, *displ2))) ) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  // some setup for further testing
  if (mesh1) 
    temp_mesh1 = *mesh1;
  else
    fc_getMeshFromSubset(*subset1, &temp_mesh1);
  if (mesh2)
    temp_mesh2 = *mesh2;
  else
    fc_getMeshFromSubset(*subset2, &temp_mesh2);
  // more tests
  if ( (displ1 && !fc_isValidDisplacementVariable(temp_mesh1, *displ1)) ||
       (displ2 && !fc_isValidDisplacementVariable(temp_mesh2, *displ2))) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  fc_getMeshDim(temp_mesh1, &numDim1);
  fc_getMeshDim(temp_mesh2, &numDim2);
  if (numDim1 != numDim2) {
    fc_printfErrorMessage("Meshes must have same numDim (%d != %d)",
			  numDim1, numDim2);
    return FC_INPUT_ERROR;
  }

  // made it past tests
  return FC_SUCCESS;
}

/** 
 * \ingroup PrivateGeometricRelations
 * \brief Return the proximal vertices between two displaced subsets.
 *
 * \description
 *   
 *    This returns all vertex pairs that are less than or equal
 *    to min_dist apart. Vertices may show up in more than one pair.
 *    (In the extrememe case if you put in min_dist = infinity, you
 *    would end up with each vertex being connected to all vertices
 *    of the other subset.)
 *
 *    The routine can be forced to return "early", after the first
 *    match is found, if you use NULL for the last two arguments (the
 *    vertex arrays). In this case, the returned numPair is
 *    really a flag for whether a pair was found or not and will
 *    never be greater than 1 even if there are more than one
 *    possible pairs.
 *
 *    Note, passing in min_dist = 0 will find overlapping vertices,
 *    or if subsets are on the same mesh, will find common vertices.
 *
 *    Note, subsets may be closer together than any of the vertex
 *    pairs (imagine 1 vertex of one mesh very close to middle of
 *    a face of another mesh. the vert to face distance is less
 *    than distance to the face's verts).
 *
 *    Will return if error if subset1 == subset2 AND displ1 == displ2.
 *    (subset1 can be the same as subset2 if displ1 != displ2)
 *
 * \todo (performance) Restrict the brute force search to vertices 
 *     in the bounding box overlap
 *
 * \todo (performance) make the pairs grow by greater than 1 each time
 *
 * \modifications
 *    - 02/14/2006 WSD Created.
 *    - 03/16/2006 WSD changed behavior so that if NULLs are passed
 *         in for the vertIDs, will return after 1st match.
 */
FC_ReturnCode _fc_getSubsetsProximityCore(
  FC_Subset subset1, /**< input - the first subset */
  FC_Variable* displ1, /**< input - the displacements on the mesh of the
		       first subset */
  FC_Subset subset2, /**< input - the second subset */
  FC_Variable* displ2, /**< input - the displacements on the mesh of the
		       second subset */
  double min_dist, /**< input - the distance cutoff */
  int* numPair_p, /**< output - the number of vertex pairs found */
  int** mesh1VertIDs_p, /**< output - ids of the mesh1 verts (order matters),
			 user is responsible for freeing */
  int** mesh2VertIDs_p /**< output - ids of the mesh2 verts (order matters),
			 user is responsible for freeing */
) {
  FC_ReturnCode rc;
  int i, j;
  FC_Mesh mesh1, mesh2;
  int numDim, temp_numDim, overlap_flag;
  int numMember1, numMember2, *memberIDs1, *memberIDs2;
  int numVert1, numVert2, *vertIDs1, *vertIDs2;
  FC_AssociationType assoc1, assoc2;
  FC_Coords lowers1, uppers1, lowers2, uppers2, lowers3, uppers3;
  double *coords1, *coords2;
  double min_dist2, distance2; // squared versions of min_dist & distance
  int numPair = 0, *mesh1VertIDs = NULL, *mesh2VertIDs = NULL;
  int *temp1, *temp2;
  int doDispl = 0, doPairs = 0;

  // assume args all tested & defaults set

  // setup
  if (displ1 && displ2)
    doDispl = 1;
  if (mesh1VertIDs_p && mesh2VertIDs_p)
    doPairs = 1;
  fc_getMeshFromSubset(subset1, &mesh1);
  fc_getMeshFromSubset(subset2, &mesh2);

  // plus one more test:
  fc_getMeshDim(mesh1, &numDim);
  fc_getMeshDim(mesh2, &temp_numDim);
  if (numDim != temp_numDim) {
    fc_printfErrorMessage("Meshes must have same numDim (%d != %d)",
			  numDim, temp_numDim);
    return FC_INPUT_ERROR;
  }

  // Compute bounding boxes +/- min_dist/2
  if (doDispl) {
    rc = fc_getDisplacedSubsetBoundingBox(subset1, *displ1, &numDim, &lowers1, 
					  &uppers1);
    if (rc != FC_SUCCESS)
      return rc;
    rc = fc_getDisplacedSubsetBoundingBox(subset2, *displ2, &numDim, &lowers2, 
					  &uppers2);
    if (rc != FC_SUCCESS)
      return rc;
  }
  else {
    rc = fc_getSubsetBoundingBox(subset1, &numDim, &lowers1, &uppers1);
    if (rc != FC_SUCCESS)
      return rc;
    rc = fc_getSubsetBoundingBox(subset2, &numDim, &lowers2, &uppers2);
    if (rc != FC_SUCCESS)
      return rc;
  }
  for (i = 0; i < numDim; i++) {
    lowers1[i] -= min_dist/2.;
    lowers2[i] -= min_dist/2.;
    uppers1[i] += min_dist/2.;
    uppers2[i] += min_dist/2.;
  }
  rc = fc_getBoundingBoxesOverlap(numDim, lowers1, uppers1, lowers2, uppers2, 
				  &overlap_flag, &lowers3, &uppers3);
  if (rc != FC_SUCCESS)
    return rc;

  // early return if no overlap
  if (overlap_flag == 0) {
    *numPair_p = 0;
    return FC_SUCCESS;
  }

  // convert subset members to verts (if needed)
  fc_getSubsetAssociationType(subset1, &assoc1);
  fc_getSubsetAssociationType(subset2, &assoc2);
  fc_getSubsetMembersAsArray(subset1, &numMember1, &memberIDs1);
  fc_getSubsetMembersAsArray(subset2, &numMember2, &memberIDs2);
  if (assoc1 == FC_AT_VERTEX) {
    numVert1 = numMember1;
    vertIDs1 = memberIDs1;
  }
  else {
    fc_changeMeshEntityType(mesh1, assoc1, numMember1, memberIDs1,
			    FC_AT_VERTEX, 0, &numVert1, &vertIDs1);
    free(memberIDs1);
  }
  if (assoc2 == FC_AT_VERTEX) {
    numVert2 = numMember2;
    vertIDs2 = memberIDs2;
  }
  else {
    fc_changeMeshEntityType(mesh2, assoc2, numMember2, memberIDs2,
			    FC_AT_VERTEX, 0, &numVert2, &vertIDs2);
    free(memberIDs2);
  }

  // get coords
  if (doDispl) {
    rc = fc_getDisplacedMeshCoords(mesh1, *displ1, &coords1);
    if (rc != FC_SUCCESS)
      return rc;
    rc = fc_getDisplacedMeshCoords(mesh2, *displ2, &coords2);
    if (rc != FC_SUCCESS)
      return rc;
  }
  else {
    rc = fc_getMeshCoordsPtr(mesh1, &coords1);
    if (rc != FC_SUCCESS)
      return rc;
    rc = fc_getMeshCoordsPtr(mesh2, &coords2);
    if (rc != FC_SUCCESS)
      return rc;
  }

  // FIX make this better 
  // 1) have smarter realloc of pairs - grow by *2
  // 2) could just test vertices in the overlap
  // 3) set up two sets of loops - switch from the first
  //    after the first pair is found so don't have the
  //    doPairs test in every loop
  
  // If find a pair, add to pairs
  // Note: early return if not doing pairs
  min_dist2 = min_dist*min_dist;
  for (i = 0; i < numVert1; i++) {
    int vertID1 = vertIDs1[i];
    for (j = 0; j < numVert2; j++) {
      int vertID2 = vertIDs2[j];
      fc_calcSquaredEuclideanDistance(&coords1[vertID1*numDim], 
				      &coords2[vertID2*numDim], numDim,
				      &distance2);
      if (fc_lteq(distance2, min_dist2, DBL_EPSILON, DBL_MIN)) {
	if (!doPairs) {
	  *numPair_p = 1;
	  free(vertIDs1);
	  free(vertIDs2);
	  if (doDispl) {
	    free(coords1);
	    free(coords2);
	  }
	  return FC_SUCCESS;
	}
	temp1 = realloc(mesh1VertIDs, (numPair+1)*sizeof(int));
	temp2 = realloc(mesh2VertIDs, (numPair+1)*sizeof(int));
	if (!temp1 || !temp2) {
	  free(mesh1VertIDs); free(mesh2VertIDs);
	  fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
	  return FC_MEMORY_ERROR;
	}
	mesh1VertIDs = temp1;
	mesh2VertIDs = temp2;
	mesh1VertIDs[numPair] = vertID1;
	mesh2VertIDs[numPair] = vertID2;
	numPair++;
      }
    }
  }
  free(vertIDs1);
  free(vertIDs2);
  if (doDispl) {
    free(coords1);
    free(coords2);
  }

  // set return values
  *numPair_p = numPair;
  if (doPairs) { // in case numPair is zero
    *mesh1VertIDs_p = mesh1VertIDs;
    *mesh2VertIDs_p = mesh2VertIDs;
  }

  return FC_SUCCESS;
}

//@}

/** 
 *\name Range queries and variable smoothing 
 *
 * There are no subset versions of these routines since the searching
 * data structure is per mesh. To find the entities of a subset within a range:
 * do the range search on the mesh, then use fc_createSubsetIntersection()
 * to "AND" that with the original subset. The result are the members of
 * the subset that fall within the range.
 */
//-------------------------------
//@{

/**
 * \ingroup  GeometricRelations
 * \brief Create spatial bins for vertices for geometric queries
 *
 * \description
 * 
 *     Create an FC_VertexBin to use for geometric queries.  The user is
 *     responsible for freeing the created FC_VertexBin by calling
 *     fc_freeVertexBin().
 *
 *     This is a wrapper around _fc_createVertexBinCore() which does all the
 *     work.
 *   
 * \todo Check that this works with 1D & 2D coords
 *
 * \modifications 
 *    - 03/08/04 WSK, created.
 *    - 02/20/06 WSD, moved core code into _fc_createVertexBinCore so could
 *         easily extend to dislaced and subset cases.
 */
FC_ReturnCode fc_createMeshVertexBin(
  FC_Mesh mesh,           /**< input - mesh handle */
  FC_VertexBin** bin_p    /**< output - bin */
) 
{
  FC_ReturnCode rc;
  int numDim, numVertex;
  double* coords;
  FC_Coords lowers, uppers;

  // default output
  if (bin_p)
    *bin_p = NULL;

  // Check input
  if (!fc_isMeshValid(mesh) || bin_p == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Creating spatial bins for vertices.");  
  
  // Get info from mesh
  rc = fc_getMeshInfo(mesh, NULL, &numDim, &numVertex, NULL, NULL);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_getMeshCoordsPtr(mesh, &coords);
  if (rc != FC_SUCCESS)
    return rc;

  // do it
  rc = fc_getMeshBoundingBox(mesh, &numDim, &lowers, &uppers);
  if (rc != FC_SUCCESS) 
    return rc;
  rc = _fc_createVertexBinCore(numVertex, numDim, coords, lowers, uppers,
			       bin_p);
  if (rc != FC_SUCCESS)
    return rc;
  (*bin_p)->mesh = mesh;
  (*bin_p)->displ = FC_NULL_VARIABLE;
  (*bin_p)->numDim = numDim;

  return FC_SUCCESS;
}

/**
 * \ingroup  GeometricRelations
 * \brief Create spatial bins for vertices for geometric queries
 *
 * \description
 * 
 *     Create an FC_VertexBin to use for geometric queries.  This version
 *     creates bins of the displaced vertices of the mesh. The user is
 *     responsible for freeing the created FC_VertexBin by calling
 *     fc_freeVertexBin().
 *
 *     This is a wrapper around _fc_createVertexBinCore() which does all the
 *     work.
 *   
 * \todo Check that this works with 1D & 2D coords
 *
 * \modifications 
 *    - 03/08/04 WSK, created.
 *    - 02/20/06 WSD, moved core code into _fc_createVertexBinCore so could
 *         easily extend to dislaced and subset cases.
 */
FC_ReturnCode fc_createDisplacedMeshVertexBin(
  FC_Mesh mesh,          /**< input - mesh handle */
  FC_Variable displ,     /**< input - displacement variable */
  FC_VertexBin** bin_p   /**< output - bin */
) 
{
  FC_ReturnCode rc;
  int numDim, numVertex;
  double* coords;
  FC_Coords lowers, uppers;

  // default output
  if (bin_p)
    *bin_p = NULL;

  // Check input
  if (!fc_isMeshValid(mesh) || 
      !fc_isValidDisplacementVariable(mesh, displ) || !bin_p) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Creating spatial bins for displaced vertices.");  
  
  // Get info from mesh
  rc = fc_getMeshInfo(mesh, NULL, &numDim, &numVertex, NULL, NULL);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_getDisplacedMeshCoords(mesh, displ, &coords);
  if (rc != FC_SUCCESS)
    return rc;

  // do it
  rc = fc_getDisplacedMeshBoundingBox(mesh, displ, &numDim, &lowers, &uppers);
  if (rc != FC_SUCCESS) 
    return rc;
  rc = _fc_createVertexBinCore(numVertex, numDim, coords, lowers, uppers,
			       bin_p);
  free(coords);
  if (rc != FC_SUCCESS)
    return rc;
  (*bin_p)->mesh = mesh;
  (*bin_p)->displ = displ;
  (*bin_p)->numDim = numDim;

  return FC_SUCCESS;
}

/**
 * \ingroup  GeometricRelations
 * \brief Clear memory in an FC_VertexBin
 *
 * \description
 * 
 *     Releases allocated memory in an FC_VertexBin.
 *
 * \modifications 
 *    - 03/08/04 WSK, created.
 */
void fc_freeVertexBin(
  FC_VertexBin* bin    /**< input/output - bin */
) 
{
  int i, j, k;

  // Check input
  if (bin == NULL)
    return;

  // log message
  fc_printfLogMessage("Freeing vertex bins.");       

  // Clear dynamically allocated storage
  for (i = 0; i < bin->numBins[0]; i++) {
    for (j = 0; j < bin->numBins[1]; j++) {
      for (k = 0; k < bin->numBins[2]; k++) 
        free(bin->vertexIDs[i][j][k]);
      free(bin->numVertex[i][j]);
      free(bin->vertexIDs[i][j]);
    }
    free(bin->numVertex[i]);
    free(bin->vertexIDs[i]);
  }
  free(bin->numVertex);
  free(bin->vertexIDs);
  bin->numVertex = NULL;
  bin->vertexIDs = NULL;

  free(bin);
}

/**
 * \ingroup  GeometricRelations
 * \brief Find all entities within a bounding box.
 *
 * \description 
 *
 *     Takes two points describing a bounding box and returns IDs of entities
 *     in a mesh that lie completely within the bounding box (and on it's
 *     boundary), regardless of the topology the mesh. To use this routine, you
 *     must first call fc_createMeshVertexBin or
 *     fc_createDisplacedMeshVertexBin to make the vertex bins used for
 *     searching
 *
 *     The returned ID array is created by this routine, but it is the
 *     responsibility of the caller to free it.
 *
 * \todo add argument so user can request entities compeletely within
 *     range, or entities that are partially in range.
 * \todo Change routines to return a subset instead 
 *
 * \modifications 
 *    - 10/27/04, WSK, Created.
 */
FC_ReturnCode fc_getEntitiesWithinBoundingBox(
  FC_Mesh mesh,         /**< input - mesh handle */
  FC_VertexBin* bin,    /**< input - vertex bin for the mesh */
  FC_Coords lowpoint,   /**< input - lower left corner of bounding box */
  FC_Coords highpoint,  /**< input - upper left corner of bounding box */
  FC_AssociationType assoc, /**< input - type of entity */
  int* numEntity,       /**< output - number of entities within sphere */
  int **entityIDs       /**< output - array of entity IDs */
)
{
  FC_ReturnCode rc;
  int numVertex, *vertexIDs;            // number & IDs in sphere

  // default return
  if (numEntity)
    *numEntity = -1;
  if (entityIDs)
    *entityIDs = NULL;
  
  // check input
  if (!fc_isMeshValid(mesh) || !bin || !FC_HANDLE_EQUIV(mesh, bin->mesh) ||
      !fc_isBoundingBoxValid(bin->numDim, lowpoint, highpoint) ||
      !fc_isAssociationTypeValid(assoc) || assoc == FC_AT_UNKNOWN || 
      assoc == FC_AT_WHOLE_DATASET || !numEntity || !entityIDs) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Finding entities within a bounding box.");    

  // get the vertices within radius
  rc = _fc_getVertexIDsWithinBoundingBoxUsingVertexBin(mesh, bin, lowpoint,
                                   highpoint, &numVertex, &vertexIDs);
  if (rc != FC_SUCCESS)
    return rc;

  // if we asked for vertices, we are done!
  if (assoc == FC_AT_VERTEX) {
    *numEntity = numVertex;
    *entityIDs = vertexIDs;
    return FC_SUCCESS;
  }

  // convert to requested entity type
  rc = fc_changeMeshEntityType(mesh, FC_AT_VERTEX, numVertex, vertexIDs,
			       assoc, 1, numEntity, entityIDs);
  if (rc != FC_SUCCESS)
    return rc;

  // cleanup
  free(vertexIDs);

  return FC_SUCCESS;
}

/**
 * \ingroup  GeometricRelations
 * \brief Find all entities within a sphere.
 *
 * \description 
 *
 *     Takes a point and a radius, and returns IDs of entities in a mesh that
 *     lie completely within the radius of the point (i.e. the elements whose
 *     vertices all fall within the radius of the point), regardless of the
 *     topology the mesh. To use this routine, you must first call
 *     fc_createMeshVertexBin or fc_createDisplacedMeshVertexBin to make the
 *     vertex bins used for searching.
 *
 *     The returned ID array is created by this routine, but it is the
 *     responsibility of the caller to free it.
 *
 * \todo add argument so user can request entities compeletely within
 *     range, or entities that are partially in range.
 * \todo Changed routines to return a subset instead 
 *
 * \modifications 
 *    - 10/14/04 WSK, Created to replace
 *          fc_getVertexIDsWithinSphereUsingBertexBin and
 *          fc_getElementIDsWithinSphereUsingVertexBin with a more general
 *          purpose routine.
 */
FC_ReturnCode fc_getEntitiesWithinSphere(
  FC_Mesh mesh,       /**< input - mesh handle */
  FC_VertexBin* bin,  /**< input - vertex bin for the mesh */
  FC_Coords center,   /**< input - center of the sphere */
  double radius,      /**< input - radius of the sphere */
  FC_AssociationType assoc, /**< input - type of entity */
  int* numEntity,     /**< output - number of entities within sphere */
  int **entityIDs    /**< output - array of entity IDs */
)
{
  FC_ReturnCode rc;
  int numVertex, *vertexIDs;            // number & IDs in sphere

  // default return
  if (numEntity)
    *numEntity = -1;
  if (entityIDs)
    *entityIDs = NULL;
  
  // check input
  if (!fc_isMeshValid(mesh) || !bin || !FC_HANDLE_EQUIV(mesh, bin->mesh) ||
      !center || radius < 0 || !fc_isAssociationTypeValid(assoc) || 
      assoc == FC_AT_UNKNOWN || assoc == FC_AT_WHOLE_DATASET || 
      !numEntity || !entityIDs) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Finding entities within a sphere.");    

  // get the vertices within radius
  rc = _fc_getVertexIDsWithinSphereUsingVertexBin(mesh, bin, center,
                                   radius, &numVertex, &vertexIDs);
  if (rc != FC_SUCCESS)
    return rc;

  // if we asked for vertices, we are done!
  if (assoc == FC_AT_VERTEX) {
    *numEntity = numVertex;
    *entityIDs = vertexIDs;
    return FC_SUCCESS;
  }

  // convert to requested entity type
  rc = fc_changeMeshEntityType(mesh, FC_AT_VERTEX, numVertex, vertexIDs,
			       assoc, 1, numEntity, entityIDs);
  if (rc != FC_SUCCESS)
    return rc;

  // cleanup
  free(vertexIDs);

  return FC_SUCCESS;
}

/**
 * \ingroup  GeometricRelations
 * \brief Smooth a variable field
 *
 * \description
 * 
 *    This routine creates a "smoothed" version of the original variable where
 *    each data point in the smoothed variable is the average of all data
 *    points within the given radius of the original data point. It uses the
 *    vertex bins created using fc_createMeshVertexBin(). Non-scalar variables
 *    are averaged component-wise. (Note that only entities that completely
 *    fall within the radius are used for the smoothing.) The smoothed
 *    variable will always have data type of double.
 *
 *    Because the vertex bins are per mesh, if you want to smooth a variable
 *    over all meshes, you must provide handles for the variable on all
 *    of the meshes and the appropriate vertex bins).
 *
 * \modifications 
 *    - 03/09/04 WSK, Created.
 *    - 04/26/04 WSK, added handling for smoothing element associated vars
 */
FC_ReturnCode fc_geomSmoothVariable(
  int numMesh,         /**< input - the number of meshes the variable is on */
  FC_Mesh* meshes,     /**< input - the meshes */ 
  FC_VertexBin** bins, /**< input - for each mesh, its vertex bin */
  FC_Variable* vars,   /**< input - for each mesh, the variable to smooth */
  double radius,       /**< input - the kernel radius */
  int doLocalSmooth,   /**< input - if true, only points on local mesh are used for smoothing 
                                    if false, points on all meshes are used for smoothing */
  FC_Variable** smooth_vars /**< output - the smoothed vars */
) {
  return _fc_geomSmoothCore(numMesh, meshes, 0, NULL, bins, vars,
			    radius, doLocalSmooth, smooth_vars);
}

/**
 * \ingroup  GeometricRelations
 * \brief Smooth a variable field on a displaced mesh
 *
 * \description 
 *
 *    This routine creates a "smoothed" version of the original variable where
 *    each data point in the smoothed variable is the average of all data
 *    points within the given radius of the original data point. It uses the
 *    vertex bins created using fc_createDisplacedMeshVertexBin(). Non-scalar variables
 *    are averaged component-wise. (Note that only entities that completely
 *    fall within the radius are used for the smoothing.) The smoothed
 *    variable will always have data type of double.
 *
 *    Because the vertex bins are per mesh, if you want to smooth a variable
 *    over all meshes, you must provide handles for the variable on all
 *    of the meshes and the appropriate vertex bins).
 *
 * \modifications 
 *    - 03/09/04 WSK, Created.
 *    - 04/26/04 WSK, added handling for smoothing element associated vars
 */
FC_ReturnCode fc_displacedGeomSmoothVariable(
  int numMesh,         /**< input - the number of meshes the variable is on */
  FC_Mesh* meshes,     /**< input - the meshes */ 
  FC_Variable* displs, /**< input - for each mesh, the displacement variable */
  FC_VertexBin** displBins, /**< input - for each mesh, its displaced vertex 
			       bin */
  FC_Variable* vars,   /**< input - for each mesh, the variable to smooth */
  double radius,       /**< input - the kernel radius */
  int doLocalSmooth,       /**< input - flag for whether smoothing is done over one mesh or all */
  FC_Variable** smooth_vars /**< output - the smoothed vars */
) {
  return _fc_geomSmoothCore(numMesh, meshes, 1, displs, displBins, vars,
			    radius, doLocalSmooth, smooth_vars);
}

/**
 * \ingroup  GeometricRelations
 * \brief Smooth a sequence variable field
 *
 * \description 
 *    
 *    Wrapper that calls fc_geomSmoothVariables() under the hood.
 *
 * \modifications 
 *    - 10/11/04 WSK Created.
 */
FC_ReturnCode fc_geomSmoothSeqVariable(
  int numMesh,         /**< input - the number of meshes to provide values */
  FC_Mesh* meshes,     /**< input - the meshes */ 
  FC_VertexBin** bins, /**< input - for each mesh, its vertex bin */
  int numStep,         /**< input - the number of steps for each seq variable */
  FC_Variable** seqVars, /**< input - for each mesh, the seq vars to smooth */
  double radius,       /**< input - the kernel radius */
  int doLocalSmooth,       /**< input - flag for whether smoothing is done over one mesh or all */
  FC_Variable*** smooth_seq_vars /**< output: the smoothed seq vars */
) {
  FC_ReturnCode rc;
  int i, j;
  FC_Mesh temp_mesh;
  FC_Sequence sequence;
  FC_Variable *temp_smooth;
  FC_Variable **smooth_vars, *temp_seqVar;
  
  // default output
  if (smooth_seq_vars)
    *smooth_seq_vars = NULL;

  if (numMesh < 1 || !meshes || !bins || numStep < 0 || !seqVars ||
      radius <= 0. || !smooth_seq_vars) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  for (i = 0; i < numMesh; i++) {
    if (!fc_isMeshValid(meshes[i]) || !bins[i] ||
        !fc_isSeqVariableValid(numStep, seqVars[i]) ||
	fc_isVariableGlobal(seqVars[i][0])) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
      return FC_INPUT_ERROR;
    }
    fc_getMeshFromVariable(seqVars[i][0], &temp_mesh);
    if (!FC_HANDLE_EQUIV(bins[i]->mesh, meshes[i]) || 
        !FC_HANDLE_EQUIV(temp_mesh, meshes[i]) ) {
      fc_printfErrorMessage("bins or vars do not match meshes");
      return FC_INPUT_ERROR;
    }
  }

  // make room for smoothed vars (this will become *smooth_seq_vars)
  smooth_vars = (FC_Variable**)malloc(numMesh*sizeof(FC_Variable*));
  if (!smooth_vars) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  for (i = 0; i < numMesh; i++) {
    smooth_vars[i] = (FC_Variable*)malloc(numStep*sizeof(FC_Variable));
    if (!smooth_vars[i]) {
      for (j = 0; j < i-1; j++)
        free(smooth_vars[j]);
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      return FC_MEMORY_ERROR;
    }
  }

  // smooth each time step (have to rearrange handles going in & out)
  for (i = 0; i < numStep; i++) {
    FC_Variable temp_vars[numMesh];
    for (j = 0; j < numMesh; j++)
      temp_vars[j] = seqVars[j][i];
    rc = fc_geomSmoothVariable(numMesh, meshes, bins, temp_vars,
                               radius, doLocalSmooth, &temp_smooth);
    if (rc != FC_SUCCESS)
      return rc;
    for (j = 0; j < numMesh; j++)
      smooth_vars[j][i] = temp_smooth[j];
    free(temp_smooth);
  }
  
  // convert the smoothed vars to seq vars
  rc = fc_getSequenceFromSeqVariable(numStep, seqVars[0], &sequence);
  if (rc != FC_SUCCESS)
    return rc;
  for (i = 0; i < numMesh; i++) {
    rc = fc_convertVariablesToSeqVariable(numStep, smooth_vars[i],
                                          sequence, NULL, &temp_seqVar);
    if (rc != FC_SUCCESS)
      return rc;
    free(temp_seqVar);
  }

  // our array of basic handles is the same as the array of seqVar handles
  *smooth_seq_vars = smooth_vars;

  return FC_SUCCESS;
}

/**
 * \ingroup  GeometricRelations
 * \brief Smooth each step of a sequence variable on a displaced mesh.
 *
 * \description 
 *    
 *    Wrapper that calls fc_displacedGeomSmoothVariables() under the hood.
 *    NOTE that unlike fc_geomSmoothSeqVariable, one has to input an
 *    array of displaced vertex bins for each step of the seqVariable.
 *
 * \modifications 
 *    - 2/28/06 WSD Created.
 */
FC_ReturnCode fc_displacedGeomSmoothSeqVariable(
  int numMesh,         /**< input - the number of meshes to provide values */
  FC_Mesh* meshes,     /**< input - the meshes */ 
  int numStep1,        /**< input - the number of steps for each displ variable */
  FC_Variable** seqDispls, /**< input - for each mesh, the displ seq var */
  FC_VertexBin*** displBins, /**< input - for each mesh, its displaced vertex 
				bins */
  int numStep2,        /**< input - the number of steps for each seq variable */
  FC_Variable** seqVars, /**< input - for each mesh, the seq vars to smooth */
  double radius,       /**< input - the kernel radius */
  int doLocalSmooth,       /**< input - flag for whether smoothing is done over one mesh or all */
  FC_Variable*** smooth_seq_vars /**< output: the smoothed seq vars */
) {
  FC_ReturnCode rc;
  int i, j;
  FC_Mesh temp_mesh;
  FC_Sequence sequence;
  FC_Variable *temp_smooth;
  FC_Variable **smooth_vars, *temp_seqVar;
  
  // default output
  if (smooth_seq_vars)
    *smooth_seq_vars = NULL;

  if (numMesh < 1 || !meshes || numStep1 < 0 || !seqDispls || !displBins || 
      numStep2 < 0 || numStep2 != numStep1 || !seqVars || radius <= 0. || 
      !smooth_seq_vars) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  for (i = 0; i < numMesh; i++) {
    if (!fc_isMeshValid(meshes[i]) || 
	!fc_isSeqVariableValid(numStep1, seqDispls[i]) || !displBins[i] ||
        !fc_isSeqVariableValid(numStep2, seqVars[i]) ||
	fc_isVariableGlobal(seqVars[i][0])) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
      return FC_INPUT_ERROR;
    }
    for (j = 0; j < numStep1; j++) {
      if (!displBins[i][j] || !FC_HANDLE_EQUIV(displBins[i][j]->mesh, meshes[i]) || 
	  !fc_isValidDisplacementVariable(meshes[i], seqDispls[i][j])) {
	fc_printfErrorMessage("bad bin or not a valid displ var");
	return FC_INPUT_ERROR;
      }
    }
    fc_getMeshFromVariable(seqVars[i][0], &temp_mesh);
    if (!FC_HANDLE_EQUIV(temp_mesh, meshes[i])) {
      fc_printfErrorMessage("vars do not match meshes");
      return FC_INPUT_ERROR;
    }
  }

  // make room for smoothed vars (this will become *smooth_seq_vars)
  smooth_vars = (FC_Variable**)malloc(numMesh*sizeof(FC_Variable*));
  if (!smooth_vars) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  for (i = 0; i < numMesh; i++) {
    smooth_vars[i] = (FC_Variable*)malloc(numStep1*sizeof(FC_Variable));
    if (!smooth_vars[i]) {
      for (j = 0; j < i-1; j++)
        free(smooth_vars[j]);
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      return FC_MEMORY_ERROR;
    }
  }

  // smooth each time step (have to rearrange handles going in & out)
  for (i = 0; i < numStep1; i++) {
    FC_Variable temp_displs[numMesh];
    FC_Variable temp_vars[numMesh];
    FC_VertexBin* temp_bins[numMesh];
    for (j = 0; j < numMesh; j++) {
      temp_displs[j] = seqDispls[j][i];
      temp_vars[j] = seqVars[j][i];
      temp_bins[j] = displBins[j][i];
    }
    rc = fc_displacedGeomSmoothVariable(numMesh, meshes, temp_displs, temp_bins,
					temp_vars, radius, doLocalSmooth, &temp_smooth);
    if (rc != FC_SUCCESS)
      return rc;
    for (j = 0; j < numMesh; j++)
      smooth_vars[j][i] = temp_smooth[j];
    free(temp_smooth);
  }
  
  // convert the smoothed vars to seq vars
  rc = fc_getSequenceFromSeqVariable(numStep1, seqVars[0], &sequence);
  if (rc != FC_SUCCESS)
    return rc;
  for (i = 0; i < numMesh; i++) {
    rc = fc_convertVariablesToSeqVariable(numStep1, smooth_vars[i],
                                          sequence, NULL, &temp_seqVar);
    if (rc != FC_SUCCESS)
      return rc;
    free(temp_seqVar);
  }

  // our array of basic handles is the same as the array of seqVar handles
  *smooth_seq_vars = smooth_vars;

  return FC_SUCCESS;
}

//@}

/** \name Range query helpers */
//-------------------------------
//@{

/**
 * \ingroup  PrivateGeometricRelations
 * \brief Create spatial bins for vertices for geometric queries
 *
 * \description
 * 
 *     Create an FC_VertexBin to use for geometric queries.
 *     The user is responsible for freeing the created FC_VertexBin
 *     by calling fc_freeVertexBin().
 *
 *     The algorithm attempts to create bins of roughly equal width
 *     in each dimension, and that hold on average about 1 vertex each.
 *     Given the Widths of the bounding box of the mesh (Wx,Wy,Wz)
 *     and the number of vertices (numVertex) and our requirement
 *     of one vertex per bin we can write a conservation of vertices:
 *
 *     Wx * Wy * Wz ~ (Wbin)^3 * numVertex
 *
 *     where Wbin is the bin width for each dimension. However, the bin
 *     widths are restricted so that the number of bins in the
 *     bounding box of the mesh is an integer. Given that 
 * 
 *     Wbin-x * nx = Wx
 *
 *     We first choose nx (number of bins in x direction) as
 *
 *     nx = (truncate to int) Wx * cubic_root(numVertex / (Wx*Wy*Wz))
 *     
 *     (except if nx < 1, then nx = 1) and then calculate the actual bin 
 *     width
 *
 *     Wbin-x = Wx / nx
 *
 *     And repeat for z and y.
 *
 * \todo - Check that this works with 2D coords
 *
 * \modifications 
 *    - 03/08/04 WSK, created.
 */
FC_ReturnCode _fc_createVertexBinCore(
  int numVertex, /**< input - number of verts to go into the vertex bin */
  //  int* vertIDs,     /**< input - dimensionsality of * coords */
  int numDim,           /**< input - dimensionsality of * coords */
  double* coords,       /**< input - coords of vertices */
  FC_Coords lowers,     /**< input - the bounding box of the point set */
  FC_Coords uppers,     /**< input - the bounding box of the point set */
  FC_VertexBin** bin_p  /**< output - bin */
) 
{
  //FC_ReturnCode rc;
  int i, j, I, J, K;
  double volume, density;   // number of points per unit volume, but cube root
  int icoords[3] = { 0, 0, 0 };          // i j k coords of a bin
  void* temp;
  FC_VertexBin* bin;

  // default output
  if (bin_p)
    *bin_p = NULL;

  // Check input
  if (bin_p == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // log message
  fc_printfLogMessage("Creating spatial bins for vertices.");

  // Setup basics
  bin = calloc(1, sizeof(FC_VertexBin));
  if (bin == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  for (i = 0; i < numDim; i++) {
    bin->lowerCoords[i] = lowers[i];
    bin->upperCoords[i] = uppers[i];
    bin->widths[i] = bin->upperCoords[i] - bin->lowerCoords[i];
  }

  // determine number of bins (at least 1 bin in each dimension)
  volume = 1;
  for (i = 0; i < numDim; i++)
    volume *= bin->widths[i];
  density = pow(numVertex/volume, 1./3.);
  for (i = 0; i < numDim; i++) {
    bin->numBins[i] = (int)(bin->widths[i]*density); // truncate
    if (bin->numBins[i] < 1) {
      bin->numBins[i] = 1;
      bin->binWidths[i] = bin->widths[i];
    }
    else
      bin->binWidths[i] = bin->widths[i]/bin->numBins[i];
  }

  // make space for bins
  bin->numVertex = malloc(bin->numBins[0]*sizeof(int**));
  bin->vertexIDs = malloc(bin->numBins[0]*sizeof(int***));
  if ( bin->numVertex == NULL ||  bin->vertexIDs == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  } 
  for (i = 0; i < bin->numBins[0]; i++) {
    bin->numVertex[i] = malloc(bin->numBins[1]*sizeof(int*));
    bin->vertexIDs[i] = malloc(bin->numBins[1]*sizeof(int**));
    if ( bin->numVertex[i] == NULL ||  bin->vertexIDs[i] == NULL) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      return FC_MEMORY_ERROR;
    } 
    for (j = 0; j < bin->numBins[1]; j++) {
      bin->numVertex[i][j] = calloc(bin->numBins[2], sizeof(int));
      bin->vertexIDs[i][j] = calloc(bin->numBins[2], sizeof(int*));
      if ( bin->numVertex[i][j] == NULL ||  bin->vertexIDs[i][j] == NULL) {
        fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
        return FC_MEMORY_ERROR;
      } 
    }
  }

  // Put vertices in the bins
  // Usually cells on boundaries of bins land in next bin up
  // (except for highest bin)
  for (i = 0; i < numVertex; i++) {
    for (j = 0; j < numDim; j++) { 
      icoords[j] = (int) ((coords[i*numDim+j] - bin->lowerCoords[j]) / 
        bin->binWidths[j]);
      // check boundaries
      if (icoords[j] < 0)
        icoords[j] = 0;
      else if (icoords[j] == bin->numBins[j])
        icoords[j] = bin->numBins[j] - 1;
    }
    for (j = numDim; j < 3; j++)
      icoords[j] = 0;
    I = icoords[0];
    J = icoords[1];
    K = icoords[2];
    temp = realloc(bin->vertexIDs[I][J][K], 
                   (1 + bin->numVertex[I][J][K])*sizeof(int));
    if (temp == NULL) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      return FC_MEMORY_ERROR;
    }
    bin->vertexIDs[I][J][K] = temp;
    bin->vertexIDs[I][J][K][bin->numVertex[I][J][K]] = i;
    bin->numVertex[I][J][K]++;
  }
  
  *bin_p = bin;
  return FC_SUCCESS;
}

/**
 * \ingroup  PrivateGeometricRelations
 * \brief Find all vertices within a bounding box
 *
 * \description
 * 
 *     Takes two points describing a bounding box and returns the IDs of all 
 *     vertices which are within the bounding box (including those on its
 *     boundary), regardless of the connectivity information in the mesh. 
 *     To use this routine, you must first call fc_createMeshVertexBin to make the 
 *     vertex bins used for searching.
 *
 *     The returned ID array is created by this routine, but it is the
 *     responsibility of the caller to free it.
 *
 * \modifications 
 *    - 10/27/04 WSK, Created.
 *    - 2/21/04 WSD can now handle displaced mesh
 */
FC_ReturnCode _fc_getVertexIDsWithinBoundingBoxUsingVertexBin(
  FC_Mesh mesh,         /**< input - mesh handle */
  FC_VertexBin* bin,    /**< input - vertex bin for mesh */
  FC_Coords lowpoint,   /**< input - lower corner of bb */
  FC_Coords highpoint,  /**< input - upper corner of bb */
  int* numVertex_p,     /**< output - number of vertices within radius */
  int **vertexIDs_p     /**< output - array of vertex IDs */
) 
{
  FC_ReturnCode rc;
  int i, j, k, m, n;
  int numVert, maxNumVert, *vertIDs, *temp;
  int lowerICoords[3] = { 0, 0, 0 };
  int upperICoords[3] = { 0, 0, 0 };
  double lower, upper;
  int intersects = 1; // true;
  double* coords;
  int do_displ;

  // default return
  if (numVertex_p)
    *numVertex_p = -1;
  if (vertexIDs_p)
    *vertexIDs_p = NULL;
  
  // check input
  if (!fc_isMeshValid(mesh) || !bin || !FC_HANDLE_EQUIV(mesh, bin->mesh) || 
      !fc_isBoundingBoxValid(bin->numDim, lowpoint, highpoint) || 
      !numVertex_p || !vertexIDs_p) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // log message
  fc_printfLogMessage("Finding vertices within given bounding box.");  
  
  // determine i,j,k bounding box for bins = intersection with bounding box
  for (i = 0; i < bin->numDim; i++) {
    lower = lowpoint[i];
    upper = highpoint[i];
    // special case - no intersection
    if (lower > bin->upperCoords[i] || upper < bin->lowerCoords[i]) {
      intersects = 0; // false
      break;
    }
    else {
      lowerICoords[i] = (int) ((lower - bin->lowerCoords[i])/bin->binWidths[i]);
      upperICoords[i] = (int) ((upper - bin->lowerCoords[i])/bin->binWidths[i]);
      if (lowerICoords[i] < 0)
        lowerICoords[i] = 0;
      else if (lowerICoords[i] >= bin->numBins[i])
        lowerICoords[i] = bin->numBins[i] - 1;
      if (upperICoords[i] < 0)
        upperICoords[i] = 0;
      else if (upperICoords[i] >= bin->numBins[i])
        upperICoords[i] = bin->numBins[i] - 1;
    }
  }

  // early return
  if (intersects == 0) {
    *numVertex_p = 0;
    return FC_SUCCESS; 
  }

  // Setup
  if (fc_isVariableValid(bin->displ)) {
    do_displ = 1;
    rc = fc_getDisplacedMeshCoords(mesh, bin->displ, &coords);
  }
  else {
    do_displ = 0;
    rc = fc_getMeshCoordsPtr(mesh, &coords);
  }
  if (rc != FC_SUCCESS)
    return rc;
  numVert = 0;
  maxNumVert = 0;
  vertIDs = NULL;

  // just check vertices in bins in i,j,k bounding box
  for (i = lowerICoords[0]; i <= upperICoords[0]; i++) {
    for (j = lowerICoords[1]; j <= upperICoords[1]; j++) {
      for (k = lowerICoords[2]; k <= upperICoords[2]; k++) {
        for (m = 0; m < bin->numVertex[i][j][k]; m++) {
	  int vertID = bin->vertexIDs[i][j][k][m];
	  intersects = 1;
          for (n = 0; n < bin->numDim; n++) {
	    if (coords[vertID*bin->numDim+n] < lowpoint[n] ||
		coords[vertID*bin->numDim+n] > highpoint[n]) {
	      intersects = 0;
	      break;
	    }
	  }
          if (intersects) {
	    if (numVert == maxNumVert) {
	      if (maxNumVert == 0)
		maxNumVert = 2;
	      else
		maxNumVert *= 2;
	      temp = realloc(vertIDs, maxNumVert*sizeof(int));
	      if (!temp) {
		free(vertIDs);
		fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
		return FC_MEMORY_ERROR;
	      }
	      vertIDs = temp;
	    }
            vertIDs[numVert] = bin->vertexIDs[i][j][k][m];
            numVert++;
          }
        }
      }
    }
  }
  if (do_displ)
    free(coords);

  // organize results
  if (numVert > 0) 
    fc_sortIntArray(numVert, vertIDs);
  *numVertex_p = numVert;
  *vertexIDs_p = vertIDs;

  return FC_SUCCESS;
}

/**
 * \ingroup  PrivateGeometricRelations
 * \brief Find all vertices within a given radius of a point.
 *
 * \description
 * 
 *     Takes a point and a radius, and returns the IDs of all vertices
 *     which are within the described sphere, regardless of the element
 *     (connectivity) information in the mesh. To use this routine,
 *     you must first call fc_createMeshVertexBin to make the vertex
 *     bins used for searching.
 *
 *     The returned ID array is created by this routine, but it is the
 *     responsibility of the caller to free it.
 *
 *     This routine will produce the same results as
 *     fc_markGeomNeighbors. It should be used when you plan to
 *     perform multiple queries.
 *
 * \modifications 
 *    - 03/09/04 WSK, Created.
 *    - 06/30/04 WSK Changed interface to return int array. 
 */
FC_ReturnCode _fc_getVertexIDsWithinSphereUsingVertexBin(
  FC_Mesh mesh,       /**< input - mesh handle */
  FC_VertexBin* bin,  /**< input - vertex bins */
  FC_Coords center,   /**< input - center of sphere */
  double radius,      /**< input - radius of sphere */
  int* numVertex_p,   /**< output - number of vertices within radius */
  int **vertexIDs_p   /**< output - array of vertex IDs */
) 
{
  FC_ReturnCode rc;
  int i, j, k, m, n;
  int numVert, maxNumVert, *vertIDs, *temp;
  double r2, t2, d;
  int lowerICoords[3] = { 0, 0, 0 };
  int upperICoords[3] = { 0, 0, 0 };
  double lower, upper;
  int intersects = 1; // true;
  double* coords;
  int do_displ;

  // default return
  if (numVertex_p)
    *numVertex_p = -1;
  if (vertexIDs_p)
    *vertexIDs_p = NULL;
  
  // check input
  if (!fc_isMeshValid(mesh) || !bin || !FC_HANDLE_EQUIV(mesh, bin->mesh) || 
      !center || radius < 0. || !numVertex_p || !vertexIDs_p) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // log message
  fc_printfLogMessage("Finding vertices within a sphere");  
  
  // algorithm - Most efficient algorithm would only need to search bins
  // that intersect boundary of the sphere - do not take  vertices from bins
  // completely outside of sphere and take all vertices from bins completely
  // inside of sphere. But determining intersection of sphere with the bins
  // is difficult. So instead, here we search all bins that intersect or
  // are within the bounding box of the sphere. 

  // determine i,j,k bounding box for bins = intersection with sphere's bounding box
  for (i = 0; i < bin->numDim; i++) {
    lower = center[i] - radius;
    upper = center[i] + radius;
    // special case - no intersection
    if (lower > bin->upperCoords[i] || upper < bin->lowerCoords[i]) {
      intersects = 0; // false
      break;
    }
    else {
      lowerICoords[i] = (int) ((lower - bin->lowerCoords[i])/bin->binWidths[i]);
      upperICoords[i] = (int) ((upper - bin->lowerCoords[i])/bin->binWidths[i]);
      if (lowerICoords[i] < 0)
        lowerICoords[i] = 0;
      else if (lowerICoords[i] >= bin->numBins[i])
        lowerICoords[i] = bin->numBins[i] - 1;
      if (upperICoords[i] < 0)
        upperICoords[i] = 0;
      else if (upperICoords[i] >= bin->numBins[i])
        upperICoords[i] = bin->numBins[i] - 1;
    }
  }

  // early return
  if (intersects == 0) {
    *numVertex_p = 0;
    return FC_SUCCESS; 
  }

  // Setup
  if (fc_isVariableValid(bin->displ)) {
    do_displ = 1;
    rc = fc_getDisplacedMeshCoords(mesh, bin->displ, &coords);
  }
  else {
    do_displ = 0;
    rc = fc_getMeshCoordsPtr(mesh, &coords);
  }
  if (rc != FC_SUCCESS)
    return rc;
  numVert = 0;
  maxNumVert = 0;
  vertIDs = NULL;

  // just check vertices in bins in i,j,k bounding box
  r2 = radius * radius; 
  for (i = lowerICoords[0]; i <= upperICoords[0]; i++) {
    for (j = lowerICoords[1]; j <= upperICoords[1]; j++) {
      for (k = lowerICoords[2]; k <= upperICoords[2]; k++) {
        for (m = 0; m < bin->numVertex[i][j][k]; m++) {
          t2 = 0;
          for (n = 0; n < bin->numDim; n++) {
            d = center[n] - coords[bin->vertexIDs[i][j][k][m]*bin->numDim+n];
            t2 += d * d;
          }
          if (t2 <= r2) {
	    if (numVert == maxNumVert) {
	      if (maxNumVert == 0)
		maxNumVert = 2;
	      else
		maxNumVert *= 2;
	      temp = realloc(vertIDs, maxNumVert*sizeof(int));
	      if (!temp) {
		free(vertIDs);
		fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
		return FC_MEMORY_ERROR;
	      }
	      vertIDs = temp;
	    }
            vertIDs[numVert] = bin->vertexIDs[i][j][k][m];
            numVert++;
          }
        }
      }
    }
  }
  if (do_displ)
    free(coords);

  // organize results
  if (numVert > 0)
    fc_sortIntArray(numVert, vertIDs);
  *numVertex_p = numVert;
  *vertexIDs_p = vertIDs;

  return FC_SUCCESS;
}

/**
 * \ingroup  PrivateGeometricRelations
 * \brief Smooth a variable field
 *
 * \description
 * 
 *    This routine creates a "smoothed" version of the original variable where
 *    each data point in the smoothed variable is the average of all data
 *    points within the given radius of the original data point. It uses the
 *    vertex bins created using fc_createMeshVertexBin(). Non-scalar variables
 *    are averaged component-wise. (Note that only entities that completely
 *    fall within the radius are used for the smoothing.) The smoothed
 *    variable will always have data type of double.
 *
 *    Because the vertex bins are per mesh, if you want to smooth a variable
 *    over all meshes, you must provide handles for the variable on all
 *    of the meshes and the appropriate vertex bins).
 *
 * \modifications 
 *    - 03/09/04 WSK, Created.
 *    - 04/26/04 WSK, added handling for smoothing element associated vars
 */
FC_ReturnCode _fc_geomSmoothCore(
  int numMesh,         /**< input - the number of meshes the variable is on */
  FC_Mesh* meshes,     /**< input - the meshes */ 
  int doDispl,         /**< input - flag for whether to do displacment */
  FC_Variable* displs, /**< input - for each mesh, the displacement variable */
  FC_VertexBin** bins, /**< input - for each mesh, its vertex bin */
  FC_Variable* vars,   /**< input - for each mesh, the variable to smooth */
  double radius,       /**< input - the kernel radius */
  int doLocalSmooth,       /**< input - flag for whether smoothing is done over one mesh or all */
  FC_Variable** smooth_vars /**< output - the smoothed vars */
) {
  FC_ReturnCode rc;
  int i, j, k, m, n;
  char* var_name, *temp_name;
  int dim, numComponent, temp_dim, temp_numComponent;
  int numDataPoints[numMesh];
  FC_Mesh temp_mesh;
  FC_AssociationType assoc, temp_assoc;
  FC_MathType mathtype, temp_mathtype;
  FC_DataType datatype, temp_datatype;
  int numID, *ids;
  double* coords[numMesh];
  int maxNumVertPerFace[numMesh], *numVertPerFaces[numMesh];
  int* edgeConns[numMesh], *faceConns[numMesh], *elemConns[numMesh];
  int num;
  double *sums;
  void *old_datas[numMesh];
  double *new_datas[numMesh];

  // default output
  if (smooth_vars)
    *smooth_vars = NULL;

  // check input
  if (numMesh < 1 || !meshes || (doDispl && !displs) || !bins || !vars || 
      radius <= 0. ||  !smooth_vars) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  for (i = 0; i < numMesh; i++) {
    if (doDispl && !fc_isValidDisplacementVariable(meshes[i], displs[i])) {
      fc_printfErrorMessage("displ[%d] is not a valid displ var", i);
      return FC_INPUT_ERROR;
    }
    fc_getMeshFromVariable(vars[i], &temp_mesh);
    if (!fc_isMeshValid(meshes[i]) || !bins[i] || 
	!fc_isVariableValid(vars[i]) || 
	!FC_HANDLE_EQUIV(temp_mesh, bins[i]->mesh) || 
        !FC_HANDLE_EQUIV(temp_mesh, meshes[i]) ) {
      fc_printfErrorMessage("bins or vars do not match meshes");
      return FC_INPUT_ERROR;
    }
  }
  
  // log message
  if (doDispl) {
    fc_printfLogMessage("Smoothing a variable on displaced meshes");
  }
  else {
    fc_printfLogMessage("Smoothing a variable.");
  }
  
  // get info & more checking
  rc = fc_getMeshDim(meshes[0], &dim);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_getVariableInfo(vars[0], &numDataPoints[0], &numComponent,
                          &assoc, &mathtype, &datatype);
  if (rc != FC_SUCCESS)
    return rc;
  if (datatype == FC_DT_CHAR) {
    fc_printfErrorMessage("Can't handle FC_DT_CHAR data");
    return FC_ERROR;
  }
  // All variables must have same dim, numComponent, assoc, mathtype, datatype
  for (i = 1; i < numMesh; i++) {
    fc_getMeshDim(meshes[i], &temp_dim);
    fc_getVariableInfo(vars[i], &numDataPoints[i], &temp_numComponent,
                       &temp_assoc, &temp_mathtype, &temp_datatype);
    if (temp_dim != dim || temp_numComponent != numComponent ||
        temp_assoc != assoc || temp_mathtype != mathtype ||
        temp_datatype != datatype) {
      fc_printfErrorMessage("Variables do not match");
      return FC_INPUT_ERROR;
    }
  }
  
  // For each mesh, get coords & variable data, make room for new data
  for (i = 0; i < numMesh; i++) {
    if (doDispl) 
      rc = fc_getDisplacedMeshCoords(meshes[i], displs[i], &coords[i]);
    else
      rc = fc_getMeshCoordsPtr(meshes[i], &coords[i]);
    if (rc != FC_SUCCESS)
      return rc;
    if (assoc == FC_AT_EDGE)
      fc_getMeshEdgeConnsPtr(meshes[i], &edgeConns[i]);
    else if (assoc == FC_AT_FACE)
      fc_getMeshFaceConnsPtr(meshes[i], &numVertPerFaces[i],
			     &maxNumVertPerFace[i], &faceConns[i]);
    else if (assoc == FC_AT_ELEMENT)
      fc_getMeshElementConnsPtr(meshes[i], &elemConns[i]);
    fc_getVariableDataPtr(vars[i], &old_datas[i]);
    new_datas[i] = (double*)malloc(sizeof(double)*
                                   numDataPoints[i]*numComponent);
    if (new_datas[i] == NULL) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      return FC_MEMORY_ERROR;
    }
  }

  // for each data point 
  // 1) find center
  // 2) find all other points on all meshes that fall within radius 
  //   ( 2A) make sure have current point)
  // 3) sum all contributions & divide by number of data points
  sums = (double*)malloc(sizeof(double)*numComponent);
  if (sums == NULL)  {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  for (i = 0; i < numMesh; i++) {
    for (j = 0; j < numDataPoints[i]; j++) {
      FC_Coords center;

      // Find center
      if (assoc == FC_AT_VERTEX) {
        for (k = 0; k < dim; k++)
          center[k] = coords[i][j*dim + k];
      }
      else if (assoc == FC_AT_WHOLE_MESH) {
        int numVert;
        fc_getMeshNumVertex(meshes[i], &numVert);
        for (k = 0; k < dim; k++)
          center[k] = 0;
        for (k = 0; k < numVert; k++)
          for (m = 0; m < dim; m++)
            center[m] += coords[i][k*dim + m];
        for (k = 0; k < dim; k++) 
          center[k] /= numVert;
      }
      else {
        // Centroid of this edge, face or element
        int numVert, *vertIDs;
        if (assoc == FC_AT_EDGE) {
          numVert = 2;
          vertIDs = &edgeConns[i][j*2];
        }
        else if (assoc == FC_AT_FACE) {
          numVert = numVertPerFaces[i][j];
          vertIDs = &faceConns[i][j*maxNumVertPerFace[i]];
        }
        else if (assoc == FC_AT_ELEMENT) {
	  FC_ElementType elemType;
	  fc_getMeshElementType(meshes[i], &elemType);
          numVert = fc_getElementTypeNumVertex(elemType);
          vertIDs = &elemConns[i][j*numVert];
        }
        for (k = 0; k < dim; k++)
          center[k] = 0;
        for (k = 0; k < numVert; k++)
          for (m = 0; m < dim; m++)
            center[m] += coords[i][vertIDs[k]*dim + m];
        for (k = 0; k < dim; k++) 
          center[k] /= numVert;
      }
 
      // find geometric neighbors on each mesh & add their values
      num = 0;
      for (k = 0; k < numComponent; k++)
        sums[k] = 0;
      for (k = (doLocalSmooth?i:0); k < (doLocalSmooth?i+1:numMesh); k++) {
        rc = fc_getEntitiesWithinSphere(meshes[k], bins[k], center, radius, 
                                        assoc, &numID, &ids);
        if (rc != FC_SUCCESS)
          return rc;

        // make sure we include the current entity
        if (k == i) {
          int found = 0;
          for (m = 0; m < numID; m++) {
            if (ids[m] == j) {
              found = 1;
              break;
            }
          }
          if (!found) {
            int* temp = realloc(ids, (numID+1)*sizeof(int));
            if (!temp) {
              fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
              return FC_MEMORY_ERROR;
            }
            ids = temp;
            ids[numID] = j;
            numID++;
          }
        }

        // add to our cumulative sums
        for (m = 0; m < numID; m++) {
          for (n = 0; n < numComponent; n++) {
            switch(datatype) {
            case FC_DT_INT:
              sums[n] += ((int*)(old_datas[k]))[ids[m]*numComponent+n];
              break;
            case FC_DT_FLOAT:
              sums[n] += ((float*)(old_datas[k]))[ids[m]*numComponent+n];
              break;
            case FC_DT_DOUBLE:
              sums[n] += ((double*)(old_datas[k]))[ids[m]*numComponent+n];
              break;
            default:
              return FC_ERROR;
            }
          }
          num++;
        }
        free(ids);
      }
      for (k = 0; k < numComponent; k++)
        new_datas[i][j*numComponent+k] = sums[k]/num;
    }
  }
  free(sums);
  if (doDispl) {
    for (i = 0; i < numMesh; i++)
      free(coords[i]);
  }

  // create new vars
  *smooth_vars = (FC_Variable*)malloc(numMesh*sizeof(FC_Variable));
  if (*smooth_vars == NULL) {
    fc_printfErrorMessage("memory allocation error");
    return FC_MEMORY_ERROR;
  }
  for (i = 0; i < numMesh; i++) {
    rc = fc_getVariableName(vars[i], &var_name);
    temp_name = (char*)malloc((strlen(var_name) + 10)*sizeof(char));
    sprintf(temp_name, "%s_smoothed", var_name);
    free(var_name);
    rc = fc_createVariable(meshes[i], temp_name, &(*smooth_vars)[i]);
    free(temp_name);
    if (rc != FC_SUCCESS)
      return rc;
    rc = fc_setVariableDataPtr((*smooth_vars)[i], numDataPoints[i], 
               numComponent, assoc, mathtype, FC_DT_DOUBLE, new_datas[i]);
    if (rc != FC_SUCCESS)
      return rc;
    // don't free new_datas[i]
  }
  
  return FC_SUCCESS;
}

//@}

/** \name Normals helpers */
//-------------------------------
//@{

/**
 * \ingroup  PrivateGeometricRelations
 * \brief  Compute unit normal of surface
 *
 * \description
 *
 *    This function implements Newell's method of calculating the normal of
 *    a planar polygon. It has the advantage of providing a decent estimate
 *    of a normal if the polygon is not planar.
 *
 *    Order of the points is important. Looking down on the plane,
 *    a counterclockwise traversal of the three points will yield
 *    a normal pointing up towards the viewer.
 *
 *    This routine currently assumes 3D coordinates. If you are less
 *    than 3D, pad with zeros. This can't do more than 3D.
 *
 *    Reference for Newell's method: Hill, Francis S, "Computer Graphics", 
 *    Prentice Hall, 1990, pg 258.
 *
 * \modifications  
 *    - 2003-SEP-04  W Koegler  Created: _fc_SurfaceNormal was split between
 *      this function and fc_getSurfaceNormals for better modularity.
 *    - 11/3/04 WSK Renamed & rewritten to use Newell's method (so will
 *      find estimate of surface normal of non-planar quads).
 */
void _fc_calcSurfaceNormal(
  int numPoint,       /**< input - number of points (at least 3) */
  FC_Coords *coords,  /**< input - array of point coordinates (at least 3) */
  FC_Vector *normal   /**< output - surface normal */
) {
  int i;
  double *point1, *point2; // temp pointers to points
  FC_Coords temp_norm;
  double magnitude;
 
  // log message
  fc_printfLogMessage("Computing a normal vector.");   
 
  // init 
  temp_norm[0] = 0.;
  temp_norm[1] = 0.;
  temp_norm[2] = 0.;

  // Newell method
  for (i = 0; i < numPoint; i++) {
    point1 = coords[i];
    if (i == numPoint - 1)
      point2 = coords[0];
    else
      point2 = coords[i+1];
    temp_norm[0] += (point1[1] - point2[1])*(point1[2] + point2[2]);  
    temp_norm[1] += (point1[2] - point2[2])*(point1[0] + point2[0]);  
    temp_norm[2] += (point1[0] - point2[0])*(point1[1] + point2[1]);  
  }

  // Convert to unit vector
  magnitude = 0;
  for (i = 0; i < 3; i++)
    magnitude += temp_norm[i]*temp_norm[i];
  magnitude = sqrt(magnitude);
  for (i = 0; i < 3; i++)
    (*normal)[i] = temp_norm[i]/magnitude;
}

//@}

/** \name Math helpers */
//-------------------------------
//@{

/**
 * \ingroup  PrivateGeometricRelations
 * \brief Calculate determinant of a 2x2 matrix
 *
 * \description
 *
 *    Input is by rows (top row, then bottom row).
 *
 * \modifications  
 *    - 08/23/04 WSK created.
 */
double _fc_calcDeterminant2x2(double a1, double a2,
                              double b1, double b2) {
  return a1*b2 - a2*b1;
}

/**
 * \ingroup  PrivateGeometricRelations
 * \brief Calculate determinant of a 3x3 matrix
 *
 * \description
 *
 *    Input is by rows (top row, middle row, etc).
 *
 * \modifications  
 *    - 08/23/04 WSK created.
 */
double _fc_calcDeterminant3x3(double a1, double a2, double a3,
                              double b1, double b2, double b3,
                              double c1, double c2, double c3) {
  return (  a1 * _fc_calcDeterminant2x2(b2, b3, c2, c3)
          - b1 * _fc_calcDeterminant2x2(a2, a3, c2, c3)
          + c1 * _fc_calcDeterminant2x2(a2, a3, b2, b3) );
}

//@}

/** \name Parameterization helpers */
//-------------------------------
//@{

/**
 * \ingroup  PrivateGeometricRelations
 * \brief Compute parameterization of a point on a quad
 *
 * \description
 *  
 *   Given a point on a quad, computes the parameterization location
 *   (in u & v) of that point. Combined with fc_calcQuadLocation, it
 *   can be used to find the location of a point on a quad that
 *   has been moved or distorted. (If the point is not on the surface
 *   of the quad, the closest point on the surface of the quad is
 *   parameterized.)
 *
 *   I.e. given an initial quad and a point in xyz space, the points'
 *   relative position in the quad can determined as the parameters u and
 *   v. If the quad's vertices' coordinates change. The new xzy coordinates
 *   of the point can be determined from u and v and the new quad's
 *   coordinates. See fc_calcQuadLocation.
 *
 *   Although this routine returns doubles, they are only gauranted
 *   to have as many significant digits as a float.
 *
 * \modifications   
 *    - 2003-SEP-04  Wendy Koegler  Created.
 *    - 2006-09-07 WSD Clarified exit conditions.
 */
void _fc_calcQuadParams(
  FC_Coords *vertices,     /**< input - the coordinates of the vertices of the quad */
  FC_Coords point,         /**< input - 3d coordinates of a point */
  int* numParam,           /**< output - the number of params found (0 or 2) */
  double** params          /**< output - the parameters found */
) {
  int i, j;
  int iter, MAX_ITER = 10;
  int converged;
  double diff[2], current[2], current_m[2];
  double fcol[2]; // solution vector
  double rcol[2], scol[2];  // columns of the Jacobian matrix
  double determinant;  // matrix determinant
  double weights[4];
  double derivs[8];
  double proj_point[3];
  FC_Vector normal;
  int indices[2];
  int maxID; // the index of the max component (which we will drop)
  double maxComponent;
  
  // default values
  *numParam = 0;
  *params = NULL;

  // get normal
  _fc_calcSurfaceNormal(4, vertices, &normal);

  // project point to plane of quad
  {
    double temp_point[3];
    double temp;

    // vector from point to "origin" of plane
    for (i = 0; i < 3; i++)
      temp_point[i] = point[i] - vertices[0][i];
    // compute dot product
    temp = 0;
    for (i = 0; i < 3; i++)
      temp += normal[i]*temp_point[i];
    // ??
    for (i = 0; i < 3; i++)
      proj_point[i] = point[i] - temp*normal[i];
  }

  // choose indices
  maxID = 2;  
  maxComponent = fabs(normal[2]);
  for (i = 0; i < 2; i++) {
    if (fabs(normal[i]) > maxComponent) {
      maxID = i;
      maxComponent = fabs(normal[i]);
    }
  }
  j = 0;
  for (i = 0; i < 3; i++) {
    if (i != maxID) {
      indices[j] = i;
      j++;
    }
  }

  // Use Newton's method to solve for parametric coordinates

  // initial guess - starting at (0, 0) cause that can be a problem point
  current[0] = current[1] = 0;

  converged = 0;
  for (iter = 0; iter < MAX_ITER; iter++) {
    // intermediate calcs
    current_m[0] = 1. - current[0];
    current_m[1] = 1. - current[1];
    weights[0] = current_m[0] * current_m[1];
    weights[1] =   current[0] * current_m[1];
    weights[2] =   current[0] *   current[1];
    weights[3] = current_m[0] *   current[1];
    derivs[0] = -1 * current_m[1];
    derivs[1] =      current_m[1];
    derivs[2] =        current[1];
    derivs[3] = -1 *   current[1];
    derivs[4] = -1 * current_m[0];
    derivs[5] = -1 *   current[0];
    derivs[6] =        current[0];
    derivs[7] =      current_m[0];

    // calculate newton functions
    for (i = 0; i < 2; i++) {
      fcol[i] = -1 * proj_point[indices[i]];
      rcol[i] = scol[i] = 0.0;
    }

    for (i = 0; i < 4; i++) {
      for (j = 0; j < 2; j++) { // overdetermined system, only use 2 of 3 dims
        fcol[j] += vertices[i][indices[j]] * weights[i];
        rcol[j] += vertices[i][indices[j]] * derivs[i];
        scol[j] += vertices[i][indices[j]] * derivs[i+4];
      }
    }

    determinant = rcol[0]*scol[1] - scol[0]*rcol[1];
    // FIX: check, a zero determinant is bad
    diff[0] = (fcol[0]*scol[1] - scol[0]*fcol[1])/determinant;
    diff[1] = (rcol[0]*fcol[1] - fcol[0]*rcol[1])/determinant;

    // Test for exit, and set up next guess
    if (FC_DBL_EQUIV(current[0] - diff[0], current[0]) &&
	FC_DBL_EQUIV(current[1] - diff[1], current[1])) {
      converged = 1;
      break;
    }
    
    // guess for next round
    current[0] -= diff[0];
    current[1] -= diff[1];
  } // end of iteration

  // debug
  //printf("***iter = %d, converged = %d\n", iter, converged);
  //printf("***diff        = %.20g %.20g\n", diff[0], diff[1]);
  //printf("***params      = %.20g %.20g\n", current[0], current[1]);
  //printf("***params-diff = %.20g %.20g\n", current[0]-diff[0], current[1]-diff[1]);

  // HACK - The exit condition for the loop is very stringent
  // here we have a more relaxed condition
  if (fabs(diff[0]) < FLT_EPSILON && fabs(diff[1]) < FLT_EPSILON) {
    converged = 1;
  }

  if (converged) {
    *numParam = 2;
    *params = malloc(sizeof(double)*2);
    if ( *params == NULL) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      exit(FC_MEMORY_ERROR);
    } 
    (*params)[0] = current[0];
    (*params)[1] = current[1];
  }
}

/**
 * \ingroup  PrivateGeometricRelations
 * \brief Compute parameterization of a point on a quad
 *
 * \description
 *  
 *   Given a point on a quad, computes the parameterization location
 *   (in u & v) of that point. Combined with fc_calcQuadLocation, it
 *   can be used to find the location of a point on a quad that
 *   has been moved or distorted.
 *
 *   I.e. given an initial quad and a point in xyz space, the points'
 *   relative position in the quad can determined as the parameters u and
 *   v. If the quad's vertices' coordinates change. The new xzy coordinates
 *   of the point can be determined from u and v and the new quad's
 *   coordinates. See fc_calcQuadLocation.
 *
 *   Assumes that the quad is planar. Code is based on implementation in
 *   VTK (visualization toolkit). Somebody should probably check the
 *   theoretics ....
 *
 * \modifications   
 *    - 2003-SEP-04  Wendy Koegler  Created.
 */
void _fc_calcQuadLocation(
  FC_Coords *vertices,     /**< input -  the coordinates of the vertices of the quad */
  int numParam,            /**< input - the number of params  */
  double* params,          /**< input - the parameters of the point */
  FC_Coords* point         /**< output - 3d coordinates of the point */
) {
  int i, j;
  double params_m[2];
  double weights[4];
  
  // FIX? default value?

  if (numParam != 2)
    return;

  // intermediate calcs
  params_m[0] = 1. - params[0];
  params_m[1] = 1. - params[1];
  weights[0] = params_m[0] * params_m[1];
  weights[1] =   params[0] * params_m[1];
  weights[2] =   params[0] *   params[1];
  weights[3] = params_m[0] *   params[1];

  for (i = 0; i < 3; i++)
    (*point)[i] = 0;

  for (i = 0; i < 4; i++) 
    for (j = 0; j < 3; j++) 
      (*point)[j] += vertices[i][j] * weights[i];
}

//@}

/** \name Mesh entities' measurement helpers */
//-------------------------------
//@{

/**
 * \ingroup  PrivateGeometricRelations
 * \brief Compute surface area of a triangle
 *
 * \description
 *    
 *    Provide array of three vertex coords and get the area back.
 *
 * \modifications
 *    - 8/23/04 WSK Created
 *    - 6/30/06 WSD Saw roundoff error w/ gcc 4 with Kahan's algorithm,
 *        switched back to Heron's which works fine.
 */
void _fc_calcTriArea(
  int numDim,           /**< input - the dimensionality of coordinates */
  FC_Coords *vertices,  /**< input - the coordinates of the tri */
  double* area          /**< output - surface area */
) {
  double a2, b2, c2;
  //double a, b, c, da, db, dc;

  // get lengths of sides of triangle
  fc_calcSquaredEuclideanDistance(vertices[0], vertices[1], numDim, &a2);
  fc_calcSquaredEuclideanDistance(vertices[1], vertices[2], numDim, &b2);
  fc_calcSquaredEuclideanDistance(vertices[2], vertices[0], numDim, &c2);

  // A version of Heron's formula
  *area = 0.25 * sqrt(fabs(4.0*a2*c2 - (a2-b2+c2)*(a2-b2+c2)));

  // Not used: found lacking on 6/30/06 by WSD. Switch back to old formula
  // This formula is from pdf found on the web:
  // http://www.cs.berkeley.edu/~wkahan/VtetLang.pdf
  // "What has the Volume of a Tetrahedron to do with Computer
  // Programming Languages?" 
  // by Prof. W. Kahan, Math Dept & EE & CS Dept at Berkeley
  //
  // It is supposed to be more accurate in terms of roundoff error
  // than the typical Heron's formula and other methods. 
  /*
  if (b > c)
    da = (b - a) + c;
  else 
    da = (c - a) + b;
  if (a > c)
    db = (a - b) + c;
  else
    db = (c - b) + a;
  if (a > b)
    dc = (a - c) + b;
  else
    dc = (b - c) + a;
  *area = 0.25 * sqrt((a+b+c)*da*db*dc);
  */

  return;
}

/**
 * \ingroup  PrivateGeometricRelations
 * \brief Compute surface area of a quadrilateral
 *
 * \description
 *    
 *    Provide array of vertex coords and get the area back. Note that
 *    if the coords are 2D, the last entry of FC_Coords should be 0.
 *
 *    This routine works for planar convex and planar nonconvex quads.
 *    (Will return results for non-planar quads, but behavior may not be what
 *    you want).
 *
 *    Basic procedure: divide the quad into two triangles and sum the
 *    areas of those triangles. There are two possible ways to divide a quad
 *    into triangles. This choice is important for non-convex quads and it
 *    should be chosen so that the shared edge falls within the quad. 
 *    However, figuring out which way to divide quads is difficult. Instead,
 *    the procedure is to calculate the quad area for both divisions and
 *    then to choose the smaller area (the wrong division for a non-convex
 *    quad would double count area outside of the quad).
 *
 * \modifications   
 *    - 8/23/04 WSK Created
 */
void _fc_calcQuadArea(
  int numDim,           /**< input - the dimensionality of coordinates */
  FC_Coords *vertices,  /**< input - the coordinates of the tri */
  double* area          /**< output - surface area */
) {
  int i, j, k;
  FC_Coords coords[4][3];
  int lookup[4][3] = { { 0, 1, 2 }, { 0, 2, 3 }, { 0, 1, 3 }, { 1, 2, 3 } };
  double area_tri1, area_tri2, area1, area2;

  // setup triangles
  for (i = 0; i < 4; i++) {
    for (j = 0; j < 3; j++) {
      for (k = 0; k < numDim; k++) {
        coords[i][j][k] = vertices[lookup[i][j]][k];
      }
    }
  }

  _fc_calcTriArea(numDim, coords[0], &area_tri1);
  _fc_calcTriArea(numDim, coords[1], &area_tri2);
  area1 = area_tri1 + area_tri2;
  _fc_calcTriArea(numDim, coords[2], &area_tri1);
  _fc_calcTriArea(numDim, coords[3], &area_tri2);
  area2 = area_tri1 + area_tri2;

  // return smaller area
  if (area1 < area2)
    *area = area1;
  else
    *area = area2;
  
  return;
}

/**
 * \ingroup  PrivateGeometricRelations
 * \brief Compute volume of a tetrahedron
 *
 * \description
 *    
 *    Provide array of four vertex coords and get the volume back.
 *
 *    Implementation taken from vtk: vtkTetra::ComputeVolume.
 *
 * \modifications   
 *    - 8/23/04 WSK Created
 */
void _fc_calcTetVolume(
  FC_Coords *vertices,  /**< input - the coordinates of the tri */
  double* volume        /**< output - surface area */
) {
  *volume = _fc_calcDeterminant3x3(vertices[1][0] - vertices[0][0],
                                   vertices[2][0] - vertices[0][0],
                                   vertices[3][0] - vertices[0][0],
                                   vertices[1][1] - vertices[0][1],
                                   vertices[2][1] - vertices[0][1],
                                   vertices[3][1] - vertices[0][1],
                                   vertices[1][2] - vertices[0][2],
                                   vertices[2][2] - vertices[0][2],
                                   vertices[3][2] - vertices[0][2]) / 6.0;
  return;
}

/**
 * \ingroup  PrivateGeometricRelations
 * \brief Compute volume of a pyramid
 *
 * \description
 *    
 *    Provide array of five vertex coords and get the volume back
 *    (the tip of the pyramid is the last last vertex).
 *
 *    Implementation. Divides the pyramid into 2 tets and sum their volumes.
 *    Only accurate if square face is planar.
 *
 * \modifications   
 *    - 8/24/04 WSK Created
 */
void _fc_calcPyramidVolume(
  FC_Coords *vertices,  /**< input - the coordinates of the tri */
  double* volume        /**< output - region volume */
) {
  // Note, if the square face is not concave, this works no matter which
  // way you split the tet's because if you split "wrong" you get a positve
  // and negative tet volume that add up to the right answer. I think.
  int i, j, k;
  int numDim = 3;
  FC_Coords coords[4];
  int lookup[2][4] = { { 0, 1, 3, 4 }, { 1, 2, 3, 4 } };
  double temp_vol;

  *volume = 0.0;
  for (i = 0; i < 2; i++) {
    // get coords for tet i
    for (j = 0; j < 4; j++)
      for (k = 0; k < numDim; k++) 
        coords[j][k] = vertices[lookup[i][j]][k];
    // add volume for tet i
    _fc_calcTetVolume(coords, &temp_vol);
    *volume += temp_vol;
  }

  return;
}

/**
 * \ingroup  PrivateGeometricRelations
 * \brief Compute volume of a prism
 *
 * \description
 *    
 *    Provide array of six vertex coords and get the volume back.
 *
 *    Implementation. Divides the prism into 3 tets and sums their volumes.
 *    Probably only accurate if square faces are convex and planar.
 *
 * \modifications   
 *    - 8/24/04 WSK Created
 */
void _fc_calcPrismVolume(
  FC_Coords *vertices,  /**< input - the coordinates of the tri */
  double* volume        /**< output - region volume */
) {
  int i, j, k;
  int numDim = 3;
  FC_Coords coords[4];
  int lookup[3][4] = { { 0, 1, 2, 3 }, { 1, 4, 2, 3 }, { 2, 4, 5, 3 } };
  double temp_vol;

  *volume = 0.0;
  for (i = 0; i < 3; i++) {
    // get coords for tet i
    for (j = 0; j < 4; j++)
      for (k = 0; k < numDim; k++) 
        coords[j][k] = vertices[lookup[i][j]][k];
    // add volume for tet i
    _fc_calcTetVolume(coords, &temp_vol);
    *volume += temp_vol;
  }

  return;
}

/**
 * \ingroup  PrivateGeometricRelations
 * \brief Compute volume of a hex
 *
 * \description
 *    
 *    Provide array of vertex coords and get the volume back.
 *
 *    Probably only works if the hex is convex.
 *
 *    Implementation. Divides the hex into 5 tets and sum their volumes.
 *    Probably only accurate if square faces are convex and planar.
 *
 * \modifications   
 *    - 8/23/04 WSK Created
 */
void _fc_calcHexVolume(
  FC_Coords *vertices,  /**< input - the coordinates of the tri */
  double* volume        /**< output - region volume */
) {
  int i, j, k;
  int numDim = 3;
  FC_Coords coords[4];
  int lookup[5][4] = { { 0, 1, 2, 5 }, { 0, 2, 3, 7 }, { 0, 5, 7, 4 },
                       { 2, 7, 5, 6 }, { 0, 5, 2, 7 } };
  double temp_vol;

  *volume = 0.0;
  for (i = 0; i < 5; i++) {
    // get coords for tet i
    for (j = 0; j < 4; j++)
      for (k = 0; k < numDim; k++) 
        coords[j][k] = vertices[lookup[i][j]][k];
    // add volume for tet i
    _fc_calcTetVolume(coords, &temp_vol);
    *volume += temp_vol;
  }

  return;
}

//@}


/**
 * \ingroup GeometricRelations
 * \brief Deteremine where a ray intersects with a triangle
 *
 * \description
 * 
 *   Given a ray (ie, point and direction) and a triangle, determine if and where
 *   the ray intersects the triangle. This implementation is taken from the
 *   well-known algorithm from Moller and Trumbore's 1997 "Fast, Minimum Storage
 *   Ray-Triangle Intersection" paper. 
 *
 *    returns 1 if an intersection occurred, or 0 if they do not intersect
 *
 * \modifications
 *    - 03/20/2008 CDU Initial version
 *
 */
int fc_getIntersectionBetweenRayAndTriangle(
     double ray_origin[3],    /**< input - coordinate for where the ray originates */
     double ray_direction[3], /**< input - normalized direction of the ray */
     double tri_vert0[3],     /**< input - Triangle vertex coordinate 1 */ 
     double tri_vert1[3],     /**< input - Triangle vertex coordinate 2 */ 
     double tri_vert2[3],     /**< input - Triangle vertex coordinate 3 */ 
     double *t,               /**< output - Distance from origin to the point of intersection (can be negative if opposite direction of ray) */
     double *u,               /**< output - Intersection coordinate within triangle */
     double *v){              /**< output - Intersecrion coordinate within triangle */
  
#define RTI_DEBUG (0)
#define RTI_DBG(a,b) { if(RTI_DEBUG) printf(" (%.3lf) %s\n", b, (a)?"hit":"miss"); }

  double edge1[3], edge2[3], tvec[3], pvec[3], qvec[3];
  double det,inv_det;
  double tt,uu,vv;
  const double RTI_EPSILON =  0.000001;

  #define CROSS(dest,v1,v2) \
    dest[0]=v1[1]*v2[2]-v1[2]*v2[1];	   \
    dest[1]=v1[2]*v2[0]-v1[0]*v2[2];	   \
    dest[2]=v1[0]*v2[1]-v1[1]*v2[0];
  #define DOT(v1,v2) (v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2])
  #define SUB(dest,v1,v2)	       \
    dest[0]=v1[0]-v2[0];	       \
    dest[1]=v1[1]-v2[1];	       \
    dest[2]=v1[2]-v2[2]; 

#if RTI_DEBUG
  { //Debug code for comparing intersections
    int i;
    printf("Origin:");
    for(i=0; i<3;i++) printf(" %.3lf", ray_origin[i]);
    printf("  direction:");
    for(i=0; i<3;i++) printf(" %.3lf", ray_direction[i]);
    printf("  p0:");
    for(i=0; i<3;i++) printf(" %.3lf", tri_vert0[i]);
    printf("  p1:");
    for(i=0; i<3;i++) printf(" %.3lf", tri_vert1[i]);
    printf("  p2:");
    for(i=0; i<3;i++) printf(" %.3lf", tri_vert2[i]);
  }
#endif

  //find vectors for two edges sharing vert0 
  SUB(edge1, tri_vert1, tri_vert0);
  SUB(edge2, tri_vert2, tri_vert0);

  //begin calculating determinant - also used to calculate U parameter
  CROSS(pvec, ray_direction, edge2);

  //if determinant is near zero, ray lies in plane of triangle 
  det = DOT(edge1, pvec);

  //the non-culling branch
  if (det > -RTI_EPSILON && det < RTI_EPSILON){
    RTI_DBG(0,0.0);
    return 0;
  }
  inv_det = 1.0 / det;

  //calculate distance from vert0 to ray origin */
  SUB(tvec, ray_origin, tri_vert0);

  //calculate U parameter and test bounds
  uu = DOT(tvec, pvec) * inv_det;
  if (uu < 0.0 || uu > 1.0){
    RTI_DBG(0,0.0);
    return 0;
  }

  //prepare to test V parameter
  CROSS(qvec, tvec, edge1);

  //calculate V parameter and test bounds
  vv = DOT(ray_direction, qvec) * inv_det;
  if (vv < 0.0 || uu + vv > 1.0){
    RTI_DBG(0,0.0);
    return 0;
  }

  //calculate t, ray intersects triangle
  tt = DOT(edge2, qvec) * inv_det;

  if(t) *t=tt;
  if(u) *u=uu;
  if(v) *v=vv;

  RTI_DBG(1,tt);
  return 1;

  #undef RTI_DEBUG
  #undef RTI_DBG
  #undef CROSS
  #undef DOT
  #undef SUB
}




/**
 * \ingroup PrivateGeometricRelations
 * \brief Internal function for determining whether a point is located within an element
 *
 * \description
 * 
 *   This internal function determines whether a point resides inside of an element
 *   or not. The function requires a number of inputs, including information about
 *   the mesh's faces and vertex coordinates (which can be obtained with
 *   fc_getMeshCoordsPtr() and fc_getMeshFaceConnsPtr()).
 *
 *   This function obtains all of the element's faces and then breaks each face into
 *   a collection of triangles. A ray-triangle intersection computation is performed
 *   with each triangle. If the point is inside the element, there will only be
 *   one intersection with a triangle. This function tries to detect when an
 *   intersection point happens to fall on the edge between two triangles (by verifying
 *   the distance is the same). In the case where the input point is lies in one
 *   of the element's faces, it is counted as being contained within the element.
 *
 *   note: it is possible for multiple elements to contain an element if 
 *         (a) the elements overlap, or (b) the
 *         point lies on a vertex or face the two elements have in common.
 *
 * \modifications
 *    - 03/20/2008 CDU Initial version
 *
 */
FC_ReturnCode  _fc_doesElementContainPoint(
   FC_Mesh mesh,    /**< input - The mesh this element belong to */
   int element_id,  /**< input - The id of the mesh element being examined */
   double point[3], /**< input - The xyz coordinate being evaluated */
   int dim,         /**< input - The dimension of the coordinates (must be 3) */
   int max_num_vert_per_face, /**< input - maximum number of vertices a face can have */
   int *num_vert_per_face,    /**< input - number of vertices each face has */
   int *face_conns,           /**< input - face connectivity list */
   double *vert_coords,       /**< input - coordinates for each vertex */
   int *hit_type              /**< output - Result (0=does not contain, 1=containts, 2=contains, but point is on a face) */
){

  int i,j;
  int num_faces;
  int *face_ids;


  int found_positive_hit, found_negative_hit;
  double t_val_positive, t_val_negative;
  int face_offset;
  double t;
  int hit;

  //In theory, we should do the test multiple times with randomly
  //generated directions. This helps weed out corner cases where
  //an edge/vertex happens to be in just the wrong spot. However, if we
  //compare the distance of all he intersections, and consider t=0.0 to
  //be a special case, I believe we can live with a single direction.
  double ray_direction[3] = {1.0, 0.0, 0.0};


  FC_ReturnCode rc;

  //Check inputs. Note: Mesh is checked in following step
  if(  (element_id<0) || (!point) || (!ray_direction) || 
       (max_num_vert_per_face<=0) || (!num_vert_per_face) ||
       (!vert_coords) || (!hit_type) ) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  if(dim!=3){
    fc_printfErrorMessage("Containment functions only work with 3-Dimensional data");
    return FC_INPUT_ERROR;
  }

  *hit_type = 0; //clear the output

 
  //Need to get all the faces for this element
  rc = fc_getMeshEntityChildren(mesh, 
				FC_AT_ELEMENT, element_id, //Parent: element 
				FC_AT_FACE,                //Children: faces
				&num_faces, &face_ids);
  if(rc!=FC_SUCCESS) return rc;

  found_positive_hit = found_negative_hit = 0;

  //Scan through all faces of element
  for(i=0; (i<num_faces); i++){

    if(num_vert_per_face[ face_ids[i] ] < 3){
      fc_printfErrorMessage("Element face must have at least 3 vertices for hit test");
      free(face_ids);
      return FC_ERROR;
    }

    face_offset = max_num_vert_per_face * face_ids[i]; //Shorthand

    //Break each face into a set of 1 or more triangles
    for(j=0; (j<num_vert_per_face[ face_ids[i] ] - 2); j++){
      hit =  fc_getIntersectionBetweenRayAndTriangle(
		     point, //Point is ray's origin
		     ray_direction,
                     &vert_coords[ dim * face_conns[ face_offset + 0   ] ], //triangle point 0
		     &vert_coords[ dim * face_conns[ face_offset + j+1 ] ], //triangle point 1
		     &vert_coords[ dim * face_conns[ face_offset + j+2 ] ], //triangle point 2
		     &t,NULL,NULL); //Only care about t, the distance to triangle

      if(hit){

	if((fc_eqd(t,0.0)) || (fc_eqd(t,-0.0)) ){
	  *hit_type = 2; //point is on a face
	  free(face_ids);
	  return FC_SUCCESS;

	} else if (fc_gtd(t,0.0)){

	  if(found_positive_hit){

	    if(t!=t_val_positive){
	      //Hit more than one place. Bail out.
	      free(face_ids);
	      return FC_SUCCESS;
	    }
	    //Otherwise, this must have been a common edge/vertex, keep going
	  } else {
	    //First unique hit, record distance so we can count intersections
	    t_val_positive = t;
	    found_positive_hit = 1;
	  }
	} else { //negative hit
	  
	  if(found_negative_hit){
	    if(t!=t_val_negative){
	      //Hit more than one place. Bail out.
	      free(face_ids);
	      return FC_SUCCESS;
	    }
	    //Otherwise, this must have been a common edge/vertex, keep going
	  } else {
	    //First unique hit, record distance so we can count intersecrions
	    t_val_negative = t;
	    found_negative_hit = 1;
	  }
	}
      }
    }
  }
  *hit_type = found_positive_hit && found_negative_hit;
  free(face_ids);
  return FC_SUCCESS;

}

/**
 * \ingroup GeometricRelations
 * \brief Determines whether a point is located within an element or not
 *
 * \description
 * 
 *   This function determines whether a point resides inside, outside, or exactly
 *   on the boundary of an element. Internally this function breaks an element's
 *   faces into a set of triangles and does ray-triangle intersections to
 *   determine whether the element contains the point.
 *
 *   If multiple elements are being analyzed for a single point, it is 
 *   more efficient to utilize the fc_getElementsThatContainPoint() or
 *   fc_getElementsFromSubsetThatContainPoint() functions.
 *   
 *
 * \modifications
 *    - 03/20/2008 CDU Initial version
 *
 */
FC_ReturnCode fc_doesElementContainPoint(
  FC_Mesh m,        /**< input - The mesh that the element belongs to */
  int element_id,   /**< input - The id of the element in this mesh that is to be examined */ 
  double *point,    /**< input - The xyz coordinate being evaluated */
  int *contain_type /**< output - Result (0=does not contain, 1=contains, 2=contains, but point is on surface) */
){

  FC_ReturnCode rc;
  FC_ElementType elem_type;
  int     dim, num_element, vert_per_elem;
  int    *num_vert_per_face; //ro
  int     max_num_vert_per_face;
  int    *face_conns;  //ro
  double *vert_coords; //ro

  if((!point) || (!contain_type)){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  rc = fc_getMeshInfo(m, NULL, &dim, NULL, &num_element, &elem_type); 
  if(rc!=FC_SUCCESS) return rc;

  if((dim!=3) || (elem_type<FC_ET_TET)){
    fc_printfErrorMessage("Containment functions only work with 3D coordinates and elements>FC_ET_TRIANGLE");
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  if((element_id<0)||(element_id>=num_element)){
    fc_printfErrorMessage("%s: Element id %d is out of bounds for this mesh, which has %d elements", 
			  fc_getReturnCodeText(FC_INPUT_ERROR), element_id, num_element);
    return FC_INPUT_ERROR;
  }

  vert_per_elem = fc_getElementTypeNumVertex(elem_type);

  rc = fc_getMeshCoordsPtr(m, &vert_coords);
  if(rc!=FC_SUCCESS) return rc;

  rc = fc_getMeshFaceConnsPtr(m, &num_vert_per_face, &max_num_vert_per_face, &face_conns);
  if(rc!=FC_SUCCESS) return rc;


  //Test an element
  rc = _fc_doesElementContainPoint(m, element_id, point, dim, 
		      max_num_vert_per_face, num_vert_per_face, 
		      face_conns, vert_coords, contain_type);

  return rc;
 

}



/**
 * \ingroup GeometricRelations
 * \brief Determines which elements (individually) enclose a point
 *
 * \description
 * 
 *   This function determines if a point resides inside or exactly
 *   on the boundary of each of the input elements. 
 *
 *   note: Geometry must be 3-dimensional
 *   note: Elements must be tets or higher
 *
 *   If only a single element is being examined, it may be easier to use
 *   the fc_doesElementContainPoint() function.
 *
 *   If a subset of the elements of interest is available, it may be more 
 *   convenient to use fc_getElementsFromSubsetThatContainPoint().
 *
 * \modifications
 *    - 03/20/2008 CDU Initial version
 *
 */
FC_ReturnCode fc_getElementsThatContainPoint(
  FC_Mesh m,          /**< input - The mesh that the element belongs to */
  int *element_mask,  /**< input - A mask for the mesh that specifiy which elements to consider (NULL=all) */
  double point[3],    /**< input - The xyz coordinate being evaluated */
  int *num_elem_ids,  /**< output - The number of elements that were found to contain the point */
  int **elem_ids      /**< output - The IDs of the elements that were found to contain the point */
){
 
  int i;
  FC_ReturnCode rc;
  FC_ElementType elem_type;
  int     dim, num_element, vert_per_elem;
  int    *num_vert_per_face; //ro
  int     max_num_vert_per_face;
  int    *face_conns;  //ro
  double *vert_coords; //ro

  int hit_type;

  int  res_num_elem_ids;
  int *res_elem_ids;

  if(!point || !num_elem_ids){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  rc = fc_getMeshInfo(m, NULL, &dim, NULL, &num_element, &elem_type); 
  if(rc!=FC_SUCCESS) return rc;

  if((dim!=3) || (elem_type<FC_ET_TET)){
    fc_printfErrorMessage("Containment functions only work with 3D coordinates and elements>FC_ET_TRIANGLE");
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  vert_per_elem = fc_getElementTypeNumVertex(elem_type);

  rc = fc_getMeshCoordsPtr(m, &vert_coords);
  if(rc!=FC_SUCCESS) return rc;

  rc = fc_getMeshFaceConnsPtr(m, &num_vert_per_face, &max_num_vert_per_face, &face_conns);
  if(rc!=FC_SUCCESS) return rc;


  //Allocate space for largest possible result
  res_num_elem_ids=0;
  res_elem_ids=(int *)malloc(num_element*sizeof(int));


  //Scan all elements and see which ones we fit in
  for(i=0; i<num_element; i++){

    if(element_mask && (!element_mask[i])) continue; //Bail out if 

    //Test an element
    rc = _fc_doesElementContainPoint(m, i, point, dim, 
			max_num_vert_per_face, num_vert_per_face, 
			face_conns, vert_coords, &hit_type);
		   
    if(rc!=FC_SUCCESS){      
      free(res_elem_ids);
      return rc;
    }

    //Add this element to the list if it contains point
    if(hit_type>0){
      res_elem_ids[res_num_elem_ids++] = i;
    }
  }

  //Pass the results back to the user
  if(res_num_elem_ids){
    if(elem_ids) //Only pass back if user asks
      *elem_ids = realloc(res_elem_ids, res_num_elem_ids*sizeof(int));
  } else {
    //No results, free by hand
    if(elem_ids)
      *elem_ids = NULL;
    free(res_elem_ids); 
  }
  *num_elem_ids = res_num_elem_ids;

  return FC_SUCCESS;
}

/**
 * \ingroup GeometricRelations
 * \brief Determines which elements (individually) enclose a point
 *
 * \description
 * 
 *   This function determines if a point resides inside or exactly
 *   ont the boundary of each input element listed in a subset.
 *
 *   note: Geometry must be 3-dimensional
 *   note: Elements must be tets or higher
 *
 *   If only a single element is being examined, it may be easier to use
 *   the fc_doesElementContainPoint() function.
 *
 *
 * \modifications
 *    - 03/20/2008 CDU Initial version
 *
 */
FC_ReturnCode fc_getElementsThatContainPointFromSubset(
  FC_Subset subset,  /**< input - The subset of elements that need to be inspected */
  double point[3],   /**< input - The xyz coordinate being evaluated */
  int *num_elem_ids, /**< output - The number of elements that were found to contain the point */ 
  int **elem_ids     /**< output - The IDs of the elements that were found to contain the point */
){

  FC_ReturnCode rc;
  FC_Mesh m;
  FC_AssociationType assoc;
  int  num_subset_members;
  int  mask_length;
  int *mask;
  int  tmp_num_elem_ids;
  int *tmp_elem_ids;

  if((!point)||(!num_elem_ids)){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  rc = fc_getSubsetInfo(subset, &num_subset_members, NULL, &assoc);
  if(rc!=FC_SUCCESS) return rc;


  //Set default return values
  *num_elem_ids = 0;
  if(elem_ids) *elem_ids=NULL;

  if(assoc!=FC_AT_ELEMENT){
    fc_printfErrorMessage("Association type must be FC_AT_ELEMENT");
    return FC_INPUT_ERROR;
  }

  if(!num_subset_members){
    //No items in subset, result will be empty
    return FC_SUCCESS;
  }    

  rc = fc_getMeshFromSubset(subset, &m);
  if(rc!=FC_SUCCESS) return rc;

  rc = fc_getSubsetMembersAsMask(subset, &mask_length, &mask);
  if(rc!=FC_SUCCESS) return rc;
  

  rc = fc_getElementsThatContainPoint(m, mask, point, &tmp_num_elem_ids, &tmp_elem_ids);
  free(mask);
  if(rc!=FC_SUCCESS) return rc;
  
  //Pass back to user
  *num_elem_ids = tmp_num_elem_ids;
  if(elem_ids) *elem_ids = tmp_elem_ids;

  return FC_SUCCESS;
}
