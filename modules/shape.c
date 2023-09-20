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
 * \file shape.c
 * \brief Implementation for \ref Shape modle.
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/modules/shape.c,v $
 * $Revision: 1.140 $
 * $Date: 2006/10/20 19:36:48 $
 *
 *
 * \modifications 
 *    - 12/19/05 ACG created. Moving side stuff here from topo. 
 *    - 02/08/06 ACG changes becuase of error in adjMatrix calc.
 *                   see notes of this date in individual methods.
 */

// C library includes
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

// fc library dependencies
#include "base.h"
#include "storage.h"
#include "library.h"
#include "mesh.h"
#include "meshP.h"
#include "subset.h"
#include "subsetP.h"
#include "variable.h"
#include "util.h"
#include "geom.h"
#include "geomP.h"
#include "topo.h"

// this module
#include "shape.h"
#include "shapeP.h"

/**
 * \addtogroup Shape
 * \brief Information about Shapes.
 * 
 * \description
 * 
 *   Information about Shapes (shapes of meshes etc). These may have
 *   a topological flavor to them, but since they arent heirarchical
 *   (e.g., children-parent), they are being located here rather
 *   than \ref TopologyRelations. This is meant to work in conjunction with 
 *   \ref ElemDeath e.g., given some information about the shape
 *   of a mesh is there some information about a dead element region
 *   that would be of interest, such as the region cutting through the shape ? 
 */


/**
 * \ingroup Shape
 * \defgroup PrivateShape (Private)
 */


/** 
 * \ingroup  PrivateShape
 * \brief inits the sfi array
 *
 * \description
 * 
 *  inits the sfi array
 *
 * \modifications
 *   - 06/07/06 ACG created
 */
FC_ReturnCode _fc_initSFI_Array(
  _FC_SFI_ARRAY *sfi /**< input - sorted face info array */
  ){
  sfi->numVals = 0;
  sfi->maxVals = 0;
  sfi->vals = NULL;
  
  return FC_SUCCESS;
}

/** 
 * \ingroup  PrivateShape
 * \brief add an sfi to the sfi array
 *
 * \description
 * 
 *  adds an sfi to the sfi array.
 *  returns 1 if it added it,
 *  0 if there is one already there.
 *
 * \modifications
 *   - 06/07/06 ACG created
 */
int _fc_addToSFI_Array(
  _FC_SFI_ARRAY *sfi,
  _FC_FACEINFO *fi
){
  int low, high, mid;
  int indexx;
  _FC_FACEINFO **temp;


  if (sfi->vals == NULL){ //empty
    sfi->vals = (_FC_FACEINFO**)malloc(sizeof(_FC_FACEINFO*));
    if (sfi->vals == NULL){
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      return FC_MEMORY_ERROR;
    }
    sfi->numVals = 1;
    sfi->maxVals = 1;
    sfi->vals[0] = fi; //gets the ptr, now responsible for the memory
    return 1;
  }
  
  //binary search
  low = 0;
  high = sfi->numVals-1;
  while(low <= high){
    mid = (low+high)/2;
    //      printf("low %d high %d mid %d\n",low,high,mid);
    if ((sfi->vals[mid])->faceID == fi->faceID){
      return 0;
    } else if ((sfi->vals[mid])->faceID < fi->faceID){
      low = mid+1;
    } else{
      high = mid -1;
    }
  }

  //either goes below or above mid. realloc if necessary
  indexx = ((sfi->vals[mid])->faceID > fi->faceID ? mid:mid+1);
  if (sfi->numVals >= sfi->maxVals){ //grow
    temp = (_FC_FACEINFO**)realloc(sfi->vals,(2*sfi->maxVals)*sizeof(_FC_FACEINFO*));
    if (temp == NULL){
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      return FC_MEMORY_ERROR;
    }
    sfi->maxVals*=2;
    sfi->vals = temp;
  }
  //memcpy didnt work right when i moved more than 1 element up
  memmove((sfi->vals)+indexx+1,(sfi->vals)+indexx,
	  (sfi->numVals-indexx)*sizeof(_FC_FACEINFO*));
  sfi->vals[indexx] = fi; //gets the ptr, now responsible for the memory
  sfi->numVals++;
  return 1;
  
}

/** 
 * \ingroup  PrivateShape
 * \brief looks for a fi index in sfi array
 *
 * \description
 * 
 *  looks for a fi index in sfi array. returns a pointer to the item
 *  if its there, NULL if not
 *
 * \modifications
 *   - 06/07/06 ACG created
 */
FC_ReturnCode _fc_lookupInSFI_Array(
  _FC_SFI_ARRAY *sfi,
  int i,
  _FC_FACEINFO** fi
){

  int low, high, mid;

  *fi = NULL;

  if (sfi->vals == NULL){ //empty
    return FC_SUCCESS;
  }

  //binary search
  
  low = 0;
  high = sfi->numVals-1;
  while(low <= high){
    mid = (low+high)/2;
    //      printf("low %d high %d mid %d\n",low,high,mid);
    if ((sfi->vals[mid])->faceID == i){
      *fi = sfi->vals[mid];
      return FC_SUCCESS;
    } else if ((sfi->vals[mid])->faceID < i){
      low = mid+1;
    } else{
      high = mid -1;
    }
  }

  return FC_SUCCESS;
}

/** 
 * \ingroup  PrivateShape
 * \brief removes an item from the sfi array
 *
 * \description
 * 
 *  removes an item from the sfi array. returns a pointer to it.
 *  the user is responsible for freeing it when he is done with it
 *
 * \modifications
 *   - 06/07/06 ACG created
 */
FC_ReturnCode _fc_removeFromSFI_Array(
  _FC_SFI_ARRAY *sfi,
  int i,
   _FC_FACEINFO** fi
){
  int low, high, mid;

  *fi = NULL;

  if (sfi->vals == NULL){ //empty
    return FC_SUCCESS;
  }

  //binary search
  
  low = 0;
  high = sfi->numVals-1;
  while(low <= high){
    mid = (low+high)/2;
    //      printf("low %d high %d mid %d\n",low,high,mid);
    if ((sfi->vals[mid])->faceID == i){
      *fi = sfi->vals[mid];
      //shuffle everything down
      memmove((sfi->vals)+mid,(sfi->vals)+mid+1,
	      (sfi->numVals-mid)*sizeof(_FC_FACEINFO*));
      sfi->numVals--;
      return FC_SUCCESS;
    } else if ((sfi->vals[mid])->faceID < i){
      low = mid+1;
    } else{
      high = mid -1;
    }
  }

  return FC_SUCCESS;
}

/** 
 * \ingroup  PrivateShape
 * \brief reallocs the sfi array
 *
 * \description
 * 
 *  reallocs the sfi array
 *
 * \modifications
 *   - 06/07/06 ACG created
 */
int _fc_reallocSFI_Array(
 _FC_SFI_ARRAY *sfi
 ){
  _FC_FACEINFO **temp = (_FC_FACEINFO**)realloc(sfi->vals,
						 (sfi->numVals)*sizeof(_FC_FACEINFO*));
  if (temp == NULL){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  sfi->maxVals = sfi->numVals;
  sfi->vals = temp;

  return 1;
}

/** 
 * \ingroup  PrivateShape
 * \brief prints the sfi array
 *
 * \description
 * 
 *  prints the sfi array
 *
 * \modifications
 *   - 06/07/06 ACG created
 */
void _fc_printSFI_Array(
  _FC_SFI_ARRAY *sfi
){
  int i;
  printf("SFI_Array has %d vals\n",sfi->numVals);
  for (i = 0; i < sfi->numVals; i++){
    printf("faceid %d elemid %d normal %10g %10g %10g\n",
	   (sfi->vals[i])->faceID,
	   (sfi->vals[i])->elemID,
	   (sfi->vals[i])->normal[0],
	   (sfi->vals[i])->normal[1],
	   (sfi->vals[i])->normal[2]);
  }
}


/** 
 * \ingroup  PrivateShape
 * \brief frees the sfi array
 *
 * \description
 * 
 *  frees the sfi array
 *
 * \modifications
 *   - 06/07/06 ACG created
 */
void _fc_freeSFI_Array(
  _FC_SFI_ARRAY *sfi
 ){
  int i;
  for(i = 0; i < sfi->numVals; i++){
    //free each face info
    _FC_FACEINFO *fi = sfi->vals[i];
    free(fi);
  }
  //free the array if it isnt already
  if (sfi->vals) free(sfi->vals);
  sfi->vals = NULL;
  sfi->numVals = 0;
  sfi->maxVals = 0;
}


/** 
 * \ingroup  PrivateShape
 * \brief return sorted array of ptrs to faceinfos for the subset skin
 *
 * \description 
 *
 *    This is a variation of fc_getSubsetSkin that also returns
 *    the element and the normal giving rise to those faces. Its
 *    easier to get these all at the same time.
 *
 *    There seems to be an issue in skin on how to handle an empty subset.
 *    Here I will return empty subset as an error.
 *
 *    Note: Since this is based on skin, in theory it handles
 *    not just element-faces but other topodim as well. In those
 *    cases I don't think there is a guarentee that the face-equivalent
 *    item is unique in the set.
 *
 *    This is a helper for the revised _fc_getSubsetShapes. It is
 *    a combination of the old 
 *    _fc_makeSubsetSkinFaceListWithElementData and the
 *    the old _fc_makeFaceListWithNormalsData. This was 
 *    then renamed from getDetailedSkinData
 *
 * \todo
 *   - replaced the linked list with a numFace allocated array to
 *     simulate a hash table to speed things up. will want a real hash table
 *
 * \modifications 
 *   - 06/02/06 ACG created from the elem data and the normal data
 *                  to increase efficiency of sides.
 *   - 06/05/06 ACG trying replacing the linked list with a numFace allocated
 *                  array to simulate a hash table to speed things up.
 *   - 06/07/06 ACG replacing the skinlist with a sorted array of
 *                  _FC_FACEINFO since it will need to search 
 *                  it later
 *   - 08/19/06 ACG renamed from getDetailedSkinData.
 */
FC_ReturnCode _fc_buildFaceInfoListFromSubset(
   FC_Subset subset, /**< input - subset */
   _FC_SFI_ARRAY *skinlist /**< output - sorted array of faceinfos */
  ){

  FC_ReturnCode rc;
  int i, j;
  FC_AssociationType assoc_orig, assoc_skin;
  int numMember, *memberIDs, maxNumSkin, *counts, numChild, *childIDs;
  int topodim_orig;
  FC_Mesh mesh;
  FC_ElementType elemType;

  int *temp_hasharray; //this wil be replaced by a real hashtable later


  //---default return values
  //  if (skinlist)
  //    *skinlist = NULL;

  // check input
  if (!fc_isSubsetValid(subset) ||
      !skinlist) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // log message
  fc_printfLogMessage("Creating faceinfo list");

  // Get info from subset and a little more checking
  rc = fc_getSubsetNumMember(subset, &numMember);
  if (rc != FC_SUCCESS)
    return rc;

  if (numMember < 1) {
    fc_printfErrorMessage("Can't create skin of empty subset");
    return FC_INPUT_ERROR;
  }

  rc = fc_getSubsetInfo(subset, NULL, NULL, &assoc_orig); 
  if (rc != FC_SUCCESS) 
    return rc;


  if (!(assoc_orig == FC_AT_ELEMENT || assoc_orig == FC_AT_WHOLE_MESH)){
    fc_printfErrorMessage("Can't handle subset with association type %s",
                          fc_getAssociationTypeText(assoc_orig));
    return FC_INPUT_ERROR;
  }

  rc = fc_getMeshFromSubset(subset, &mesh);
  if (rc != FC_SUCCESS) 
    return rc;

  // determine topodim of the members
  rc = fc_getMeshEntityElementType(mesh, assoc_orig, 0, &elemType);
  if (rc != FC_SUCCESS){
    free(memberIDs);
    return rc;
  }
  topodim_orig = fc_getElementTypeTopoDim(elemType);

  if (topodim_orig != 3){
    fc_printfErrorMessage("Invalid topodim %d (must be 3)\n",topodim_orig);
    return FC_INPUT_ERROR;
  }
  assoc_skin = FC_AT_FACE;

  // Get members -- If subset is whole, get all elements
  if (assoc_orig == FC_AT_WHOLE_MESH) {
    assoc_orig = FC_AT_ELEMENT;
    rc = fc_getMeshNumElement(mesh, &numMember);
    if (rc != FC_SUCCESS) 
      return rc;

    memberIDs = malloc(numMember*sizeof(int));
    if (!memberIDs) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      return FC_MEMORY_ERROR;
    }
    for (i = 0; i < numMember; i++)
      memberIDs[i] = i;
  }
  else {
    rc = fc_getSubsetMembersAsArray(subset, &numMember, &memberIDs);
    if (rc != FC_SUCCESS)
      return rc;
  }

  rc = fc_getMeshNumEntity(mesh, assoc_skin, &maxNumSkin);
  if (rc != FC_SUCCESS){
    free(memberIDs);
    return rc;
  }
  counts = calloc(maxNumSkin, sizeof(int));
  if (!counts) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }

  //adding all the faces for each element with the element id as the data.
  //if a face arises again from a different element, then remove it
  temp_hasharray = (int*)malloc(maxNumSkin*sizeof(int));
  if (temp_hasharray == NULL){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  for (i = 0; i < maxNumSkin; i++){
    temp_hasharray[i] = -1;
  }
  
  //  printf("determining the skin faces and elems\n");
  for (i = 0; i < numMember; i++) {
    fc_getMeshEntityChildren(mesh, assoc_orig, memberIDs[i], assoc_skin,
			     &numChild, &childIDs);
    for (j = 0; j < numChild; j++){
      switch(++counts[childIDs[j]]){
      case 1:{ //add it
	temp_hasharray[childIDs[j]] = memberIDs[i];
	break;
      }
      case 2:{ //remove it
	temp_hasharray[childIDs[j]] = -1;
	break;
      }
      default:  //its gone already
	//      printf("skipping face %d element %d\n",childIDs[j],memberIDs[i]);
	break;
      }
    }
    free(childIDs);
  }
  free(memberIDs);
  free(counts);
  //  printf("done determining the skin faces and elems\n");

  //  printf("buiding the SFI array and determining the normals\n");
  //since its already sorted there is a faster way to add things, but we will
  //use the interface for now

  rc = _fc_initSFI_Array(skinlist);
  for (i = maxNumSkin-1; i >= 0; i--){
    if (temp_hasharray[i] != -1){
      int test = 0;
      _FC_FACEINFO *face_data = (_FC_FACEINFO*)malloc(sizeof(_FC_FACEINFO));
      face_data->faceID = i;
      face_data->elemID = temp_hasharray[i];
      //      printf("adding face %d element %d\n",face_data->faceID,
      //	     face_data->elemID);
      rc = _fc_getOutwardFaceNormal(mesh, i, temp_hasharray[i],&(face_data->normal));
      if (rc != FC_SUCCESS){
	return FC_ERROR;
      }

      test = _fc_addToSFI_Array(skinlist,face_data);
      if (!test){
	free(temp_hasharray);
	if (face_data) free(face_data);
	//skinlist will be freed outside
	return rc;
      }
    }
  }
  //  printf("done building the skinlist and determining the normals\n");
  free(temp_hasharray);

  return FC_SUCCESS;
}


/** 
 * \ingroup  PrivateShape
 * \brief given a current shape, return sorted array of ptrs to faceinfos 
 *        for the shape sides
 *
 * \description 
 *
 *    Similiar to \ref _fc_buildFaceInfoListFromShape, but in that
 *    case you have the subset of the entire shape elements,
 *    including its innards. This method is instead called
 *    when you are going to reshape, and therefore have only
 *    the information about the sides.
 *
 * \todo
 *   - will this work ?
 *   - untested
 *
 * \modifications 
 *   - 08/19/06 ACG created
 */
FC_ReturnCode _fc_buildFaceInfoListFromShape(
  FC_Shape* shape, /**< input -shape */
   _FC_SFI_ARRAY *skinlist /**< output - sorted array of faceinfos */
  ){

  FC_ReturnCode rc;
  FC_Mesh mesh;
  int i,k;
  
  if (!shape || shape->numSides < 1){
    return FC_INPUT_ERROR;
  }


  rc = fc_getMeshFromSubset(shape->faces[0],&mesh);
  if (rc != FC_SUCCESS){
    fc_printfErrorMessage("Can't get mesh from subset");
    return rc;
  }

  rc = _fc_initSFI_Array(skinlist);
  if (rc != FC_SUCCESS){
    fc_printfErrorMessage("Can't init SFI array");
    return rc;
  }

  //each face can only be a child of one element making up the sides
  for (i = 0; i < shape->numSides; i++){
    int *farr, *earr, numf, nume;
    FC_SortedIntArray facearr;
    FC_SortedIntArray elemarr;

    rc = fc_getSubsetMembersAsArray(shape->faces[i],&numf,&farr);
    if (rc != FC_SUCCESS){
      fc_printfErrorMessage("Can't get face members");
      _fc_freeSFI_Array(skinlist);
      return rc;
    }

    rc = fc_getSubsetMembersAsArray(shape->elems[i],&nume,&earr);
    if (rc != FC_SUCCESS){
      fc_printfErrorMessage("Can't get elem members");
      free(farr);
      _fc_freeSFI_Array(skinlist);
      return rc;
    }

    rc = fc_convertIntArrayToSortedIntArray(numf,farr,1,&facearr);
    if (rc != FC_SUCCESS){
      fc_printfErrorMessage("Can't convert face array");
      free(farr);
      free(earr);
      _fc_freeSFI_Array(skinlist);
      return rc;
    }

    rc = fc_convertIntArrayToSortedIntArray(nume,earr,1,&elemarr);
    if (rc != FC_SUCCESS){
      fc_printfErrorMessage("Can't convert elem array");
      fc_freeSortedIntArray(&facearr);
      free(earr);
      _fc_freeSFI_Array(skinlist);
      return rc;
    }

    while(elemarr.numVal){
      int curre;
      int numChild, *childIDs;
      rc = fc_getSortedIntArrayBack(&elemarr,&curre);
      if (rc != FC_SUCCESS){
	fc_printfErrorMessage("Can't get last elem");
	fc_freeSortedIntArray(&facearr);
	fc_freeSortedIntArray(&elemarr);
	_fc_freeSFI_Array(skinlist);
	return rc;
      }

      rc = fc_getMeshEntityChildren(mesh,FC_AT_ELEMENT,curre,
				    FC_AT_FACE,&numChild,&childIDs);
      if (rc != FC_SUCCESS){
	fc_printfErrorMessage("Can't get elem children");
	fc_freeSortedIntArray(&facearr);
	fc_freeSortedIntArray(&elemarr);
	_fc_freeSFI_Array(skinlist);
	return rc;
      }

      for (k = 0; k < numChild;k++){
	int ret;
	if (fc_deleteIntFromSortedIntArray(&facearr,childIDs[k])){
	  _FC_FACEINFO *face_data = (_FC_FACEINFO*)malloc(sizeof(_FC_FACEINFO));
	  face_data->faceID = childIDs[k];
	  face_data->elemID = curre;
	  rc = _fc_getOutwardFaceNormal(mesh, face_data->faceID,
					face_data->elemID,
					&(face_data->normal));
	  if (rc != FC_SUCCESS){
	    fc_printfErrorMessage("Can't get normal for this face");
	    fc_freeSortedIntArray(&facearr);
	    fc_freeSortedIntArray(&elemarr);
	    _fc_freeSFI_Array(skinlist);
	    return FC_ERROR;
	  }

	  ret = _fc_addToSFI_Array(skinlist,face_data);
	  if (!ret){
	    fc_printfErrorMessage("Can't add to SFI array");
	    fc_freeSortedIntArray(&facearr);
	    fc_freeSortedIntArray(&elemarr);
	    _fc_freeSFI_Array(skinlist);
	    if (face_data) free(face_data);
	    return FC_ERROR;
	  }
	}
      }//loop over children of this element
      free(childIDs);
      fc_popSortedIntArrayBack(&elemarr);
    } //loop over elements of this side

    //now all the faces and elems should have been moved to the skinlist
    if (facearr.numVal && elemarr.numVal){
      fc_printfErrorMessage("there are still remaining faces or elems for this side");
      fc_freeSortedIntArray(&facearr);
      fc_freeSortedIntArray(&elemarr);
      _fc_freeSFI_Array(skinlist);
      return FC_ERROR;
    }

    fc_freeSortedIntArray(&facearr);
    fc_freeSortedIntArray(&elemarr);

  }//loop over sides

  return FC_SUCCESS;
}

/** 
 * \ingroup  PrivateShape
 * \brief Returns outward normal vector for a given face on a given element 
 *
 * \description 
 *
 *  Returns outward normal vector for a given face on a given element 
 *  Uses assumption of fc_calcSurfaceNormals that this is only
 *  correct for planar elements (like triangles
 *  and planar quads), but checking FaceType rather than ElemType.
 *
 * \todo
 *   - faster search
 *   - check mem allocation, esp in inner _fc_calcSurfaceNormal
 *   - does this need to be explicitly unittested? it is tested inso
 *     far as getSubsetSides and getPlanarSideNormal have been tested
 *
 * \modifications 
 *   - 10/07/05 ACG created
 *   - 12/19/2005 ACG moved from \ref TopologyRelations to \ref Shape
 */
FC_ReturnCode _fc_getOutwardFaceNormal(
  FC_Mesh mesh,        /**< input - mesh */
  int faceID,          /**< input - faceID */                                  
  int elemID,          /**< input - elemID */                                  
  FC_Vector* normal    /**< output - normal vector */
) {
  FC_ReturnCode rc;
  FC_ElementType elemType, *faceTypes;
  FC_Coords coords[4]; //restricted to tri and quad
  FC_Vector temp_normal;
  int numelem, numface, orient, numDim, numFacePerElem;
  int i,j,k;
  double* mesh_coords;
  int maxNumVertPerFace, *numVertPerFace, *faceConns, *elemToFaceConns;
  int *elemFaceOrients;

  //default
  if (normal){
    for (i = 0; i < 3; i++){
      (*normal)[i] = 0;
    }
  }

  fc_printfLogMessage("calc normal for face %d on elem %d\n", faceID, elemID);

  if (!fc_isMeshValid(mesh) || !normal){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  rc = fc_getMeshInfo(mesh, NULL, &numDim, NULL, &numelem, &elemType);
  if (rc != FC_SUCCESS) 
    return rc;

  rc = fc_getMeshNumFace(mesh,&numface);
  if (rc != FC_SUCCESS) 
    return rc;

  numFacePerElem = fc_getElementTypeNumFace(elemType);

  if (elemID < 0 || elemID > numelem || faceID < 0 || faceID > numface){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  rc = fc_getMeshFaceConnsPtr(mesh, &numVertPerFace, &maxNumVertPerFace, 
                              &faceConns);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_getMeshFaceTypesPtr(mesh, &faceTypes);
  if (rc != FC_SUCCESS) 
    return rc;

  if (faceTypes[faceID] != FC_ET_TRI && faceTypes[faceID] != FC_ET_QUAD){
    fc_printfErrorMessage("%s, can't handle face type %s",
                          fc_getReturnCodeText(FC_INPUT_ERROR),
                          fc_getElementTypeText(faceTypes[faceID]));
   return FC_INPUT_ERROR;
  }

  rc = fc_getMeshElementToFaceConnsPtr(mesh, &elemToFaceConns);
  if (rc != FC_SUCCESS) 
    return rc;
  rc = fc_getMeshElementFaceOrientsPtr(mesh, &elemFaceOrients);
  if (rc != FC_SUCCESS) 
    return rc;
  rc = fc_getMeshCoordsPtr(mesh, &mesh_coords);
  if (rc != FC_SUCCESS) 
    return rc;

  orient = 0;
  for (i = 0; i < numFacePerElem; i++){
    if (elemToFaceConns[elemID*numFacePerElem+i] == faceID){
      orient = elemFaceOrients[elemID*numFacePerElem+i];
      break;
    }
  }
  if (orient == 0){
    fc_printfErrorMessage("face %d not in face list for element %d",
                          faceID,elemID);
    return FC_INPUT_ERROR;
  }
  
  for (j = 0; j < numVertPerFace[faceID]; j++) {
    for (k = 0; k < numDim; k++)
     coords[j][k] = mesh_coords[faceConns[faceID*maxNumVertPerFace+j]*numDim+k];
    for (k = numDim; k < 3; k++)
      coords[j][k] = 0;
  }

  _fc_calcSurfaceNormal(numVertPerFace[faceID], coords, &temp_normal);
  //  printf("surface normal (%g %g %g) ",
  //     temp_normal[0],temp_normal[1], temp_normal[2]);

  //  normal = (FC_Vector*)malloc(sizeof(FC_Vector));
  for (i = 0; i < 3; i++)
    (*normal)[i] = orient*temp_normal[i];

  //  printf(" after orient (%g %g %g)\n",(*normal)[0],(*normal)[1],
  //     (*normal)[2]);
  
  return FC_SUCCESS;
}


/** 
 * \ingroup  PrivateShape
 * \brief makes the shape parts from the faceinfolist
 *
 * \description 
 *    makes the shape parts from the faceinfolist. This was the major
 *    guts of the old _fc_getSubsetSides (now defunct as a result),
 *    but I'm pulling it out in order to be able to do a reangle of
 *    an existing shape.
 *
 * \todo
 *   - will need to standardize on devouring the faceinfo list, 
 *     in case of input error etc
 *   - just writing this now, it isnt tested
 *
 * \modifications 
 *   - 10/04/05 ACG created as surfacesegment. trying to keep this
 *                  as close as possible
 *                  in algorithm to fc_segment
 *   - 11/30/05 ACG changed name from surfaceSegment to getSubsetSides to
 *                   reflect new terminology.
 *   - 01/13/06 ACG rewrote using new WSK algorithm in \ref fc_segment
 *   - 06/06/06 ACG replacing the currlist (built from the working list) with
 *                  two growing int arrays which are then handed to
 *                  the subsets. this drops out one of the (sorting) 
 *                  linked lists.
 *   - 06/07/06 ACG replacing the working list (built from the face list) with
 *                  a growing unsorted array of _FC_FACEINFO ptrs.
 *                  this drops out another sorted linked list.
 *   - 06/07/06 ACG replacing the faceinfo list (built from buildFaceInfoList)
 *                  with a growing sorted array of _FC_FACEINFO ptrs. this
 *                  drops out another sorted linked list. We wont be
 *                  adding anything to it, but the search will be faster
 *   - 06/16/06 ACG warning about creases
 *   - 06/28/06 ACG moved to private
 *   - 08/16/06 ACG created from splitting out of _fc_getSubsetSides.
 *                  renamed to createShapePartsFromFaceInfoList 
 */
FC_ReturnCode _fc_createShapePartsFromFaceInfoList(
  FC_Mesh mesh, /**< input -- the mesh */				   
  _FC_SFI_ARRAY *in_faceinfolist, /**< input - the starting array
				     with all the facedata, if this
				     method is successful,  will
				     end empty */
  char **basename,             /**< input - base name for the sides */
  double comp_angle,           /**< input - if two neighboring faces
                                differ by more than this angle then
                                they are considered to be on 
                                different surfaces.*/
  int shared_dim,              /**< input - the minimum dimensionality
				  of shared part of neighbors of faces */
  int *numSides,               /**< output - number of new subsets */
  FC_Subset** newFaceSubsets,  /**< output - array of segmented
				  regions (by face)*/
  FC_Subset** newElemSubsets,  /**< output - array of segmented
				  regions (by element)*/
  int ***adjmatrix              /**< output - adjmatrix */
  ){

  FC_ReturnCode rc;

  char* temp_name;

  _FC_INTPAIR *edgePairs = NULL;
  int numEdgePairs = 0;
  int maxEdgePairs = 0;
  int *sideArray = NULL; //temporary till I get a hash table - edge face ID to side table
  int numMeshFace = 0;

  int numSegments; //these are temporary till we set up the new Shape Object
  FC_Subset* faceSegments = NULL;
  FC_Subset* elemSegments = NULL;

  int **intmatrix = NULL;
  int i,j;
  
  //---default return values
  if (numSides)
    *numSides = -1;
  if (newFaceSubsets)
    *newFaceSubsets = NULL;
  if (newElemSubsets)
    *newElemSubsets = NULL;
  if (adjmatrix)
    *adjmatrix = NULL;

  
  // check input
  if (!in_faceinfolist || !numSides ||
      shared_dim < 0 || shared_dim > 1  || !adjmatrix ||
      !newFaceSubsets || ! newElemSubsets || comp_angle < 0 ||
      comp_angle > 180) { 
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // log message
  fc_printfLogMessage("getting shapeinfo from faceinfolist");

  if (in_faceinfolist->numVals == 0){ //return if empty subset
    *numSides = 0;
    return FC_SUCCESS;
  }

  // setup space for names of segments
  temp_name = malloc((strlen(*basename)+30)*sizeof(char));
  if (temp_name == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }

  //set up space for edgepairs
  edgePairs = (_FC_INTPAIR*)malloc(sizeof(_FC_INTPAIR));
  if (edgePairs == NULL){
    free(temp_name);
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }

  numEdgePairs = 0;
  maxEdgePairs = 1;
  //set up space for sideArray, to be replaced by a hashtable
  rc = fc_getMeshNumFace(mesh,&numMeshFace);
  if (rc != FC_SUCCESS){
    free(temp_name);
    return rc;
  }

  sideArray = (int*)malloc(numMeshFace*sizeof(int));
  if (sideArray == NULL){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  for (i=0 ; i < numMeshFace; i++){
    sideArray[i] = -1;
  }

  //  printf("segmenting faceinfo list\n");
  //---segment the faceinfo list
  numSegments = 0;
  while(in_faceinfolist->numVals){ // each pass will create a new side
    //first set up data structures....

    //info list which will contain members of the current side
    //being built. this list needs to be added to, popped off
    //either end and order doesnt matter (neither here or later).
    //contains info on the face, assoc elem, and normal.
    //therefore choosing faceinfo array. allocating space for
    //1 since we know there will be at least one in this face
    //(that is - the one we pop off the facelist
    //directly below)
    _FC_FACEINFO** working_list  = (_FC_FACEINFO**)malloc(sizeof(_FC_FACEINFO*));
    int numworkinglist = 0;
    int maxworkinglist = 1;

    //id arrays for the current side. chose (unsorted) int arrays 
    //becuase its the most efficient way to later add the IDs to the
    ///2 subsets(faces and elems) - add the whole array at once
    //(subset will do a sorting and uniqueness of the array).
    //they dont need to be sorted or searched. start by
    //allocating space for 1, since we know there wil be at least
    //1 side in this face (that is - the one we pop off the facelist
    //directly below
    int *currside_faces = (int*)malloc(sizeof(int)); 
    int *currside_elems = (int*)malloc(sizeof(int));
							
    int numcurrsideIDs = 0;
    int maxcurrsideIDs = 1;

    if (!working_list || ! currside_faces || ! currside_elems){
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      return FC_MEMORY_ERROR;
    }

    //     printf("working on side %d\n",numSegments);

    //initialize - pop off the last one and addit to the working list
    _FC_FACEINFO *temp_node = in_faceinfolist->vals[in_faceinfolist->numVals-1];
    in_faceinfolist->numVals--;
    working_list[0]= temp_node; //we know we have space for it
    numworkinglist = 1;

    // "grow" connected region from seed by checking neighbors
    while(numworkinglist){
      int currID;
      FC_Vector currNormal; 

      // remove last entity from working list and put it on
      // the 3 lists with side info
      _FC_FACEINFO *currnode = working_list[numworkinglist-1]; 
      numworkinglist--;

      sideArray[currnode->faceID] = numSegments; //fill in its side in the array (hashtable)
      if (numcurrsideIDs == maxcurrsideIDs){//grow the lists
	int* temp = (int*)realloc(currside_faces,(2*maxcurrsideIDs)*sizeof(int));
	if (temp == NULL){
	  fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
	  return FC_MEMORY_ERROR;
	}
	currside_faces = temp;
	temp = (int*)realloc(currside_elems,(2*maxcurrsideIDs)*sizeof(int));
	if (temp == NULL){
	  fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
	  return FC_MEMORY_ERROR;
	}
	currside_elems = temp;
	maxcurrsideIDs*=2;
      }
      currside_faces[numcurrsideIDs] = currnode->faceID;
      currside_elems[numcurrsideIDs] = currnode->elemID;
      numcurrsideIDs++;
      
      //we only need these if there is a neighbor, but it
      //allows us to free the currnode
      currID = currnode->faceID;
      for (i = 0; i < 3; i++){
	currNormal[i] = currnode->normal[i];
      }
      free(currnode); //free the actual info, not the ptr.

      if (in_faceinfolist->numVals){ //we have neighbors to check (otherwise
	// just loop till empty out working list
	int *neighborIDs;
	int numNeighbor;

	//	printf("checking neighbors of currid = %d\n",currID);
	// Move neighbors of that entity that satisfy the angle criteria
	// that are in the orig list to the working list
	rc = fc_getMeshEntityNeighbors(mesh, currID, FC_AT_FACE, shared_dim,
				       &numNeighbor, &neighborIDs);
	if (rc != FC_SUCCESS)
	  //clean up
	  return rc;
	
	for (i = 0; i < numNeighbor; i++) {
	  double angle;
	  _FC_FACEINFO *node;


	  rc = _fc_lookupInSFI_Array(in_faceinfolist,neighborIDs[i],&node);
	  if (rc != FC_SUCCESS){
	    fc_printfErrorMessage("cant check for key in list");
	    //clean up
	    return rc;
	  }
	  
	  if (node == NULL){ //this face is not part of the skin
	    //	    printf("face %d in not ont he skin\n",neighborIDs[i]);
	    continue;
	  }

	  //does this neighbor satisfy the angle criteria ?
	  rc = fc_calcAngleBetweenVectors(currNormal,
					  node->normal,
					  &angle);
	  if (rc != FC_SUCCESS){
	    fc_printfErrorMessage("cant get angle between vectors");
	    //clean up
	    return rc;
	  }

	  if (angle > comp_angle){ 
	    //this face is not part of this side, but it is an edgepair
	    void *temp;
	    //	    printf("face %d in not on the side, but is an edge\n",node->faceID);
	    if (numEdgePairs + 1 > maxEdgePairs) {
	      maxEdgePairs *= 2;
	      temp = realloc(edgePairs, maxEdgePairs*sizeof(_FC_INTPAIR));
	      if (temp == NULL) {
		fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
		return FC_MEMORY_ERROR;  
	      }
	      edgePairs = temp;
	    }
	    
	    edgePairs[numEdgePairs].ID1 = currID;
	    edgePairs[numEdgePairs].ID2 = node->faceID;
	    numEdgePairs++;
	    continue;
	  }
	  
	  //its part of the skin and its on this side
	  //	  printf("face %d is on the side\n",node->faceID);
	  rc = _fc_removeFromSFI_Array(in_faceinfolist,neighborIDs[i],&node);
	  if (rc != FC_SUCCESS){
	    fc_printfErrorMessage("cant remove node");
	    //clean up
	    return rc;
	  }

	  //now add that data to the workinglist
	  if (numworkinglist == maxworkinglist){ //grow the array
	    _FC_FACEINFO** temp = (_FC_FACEINFO**)realloc(working_list,
							  (2*maxworkinglist)*sizeof(_FC_FACEINFO*));
	    if (temp == NULL){
	      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
	      return FC_MEMORY_ERROR;
	    }
	    working_list = temp;
	    maxworkinglist*=2;
	  }
	  working_list[numworkinglist] = node; //gets the ptr
	  numworkinglist++;
	  
	} //end for
	free(neighborIDs);
      } //end if facelist != NULL
    } // end while working list
    free(working_list); //this is the array itself, there are no members now

    //now weve got faces and elem arrays for this side 
    { // create subsets from arrays of IDs
      void *temp;

      //      printf("making subsets from arrays\n");
      temp = realloc(faceSegments, (numSegments+1)*sizeof(FC_Subset));
      if (temp == NULL) {
	fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
	return FC_MEMORY_ERROR;
      }
      faceSegments = temp;
      sprintf(temp_name, "%s_FaceSurface%d", *basename, numSegments);
      rc = fc_createSubset(mesh, temp_name, FC_AT_FACE, &faceSegments[numSegments]);
      if (rc != FC_SUCCESS) {
	fc_printfErrorMessage("cannot create newFaceSubset[%d]", numSegments);
	//clean up
	return FC_ERROR; 
      }

      rc = fc_addArrayMembersToSubset(faceSegments[numSegments], numcurrsideIDs,
				      currside_faces);
      if (rc != FC_SUCCESS) {
	fc_printfErrorMessage("cannot add ids to facesegment");
	//clean up
	return FC_ERROR; 
      }
      
      temp = realloc(elemSegments, (numSegments+1)*sizeof(FC_Subset));
      if (temp == NULL) {
	fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
	return FC_MEMORY_ERROR;
      }
      elemSegments = temp;
      sprintf(temp_name, "%s_ElemSurface%d", *basename, numSegments);
      rc = fc_createSubset(mesh, temp_name, FC_AT_ELEMENT, &elemSegments[numSegments]);
      if (rc != FC_SUCCESS) {
	fc_printfErrorMessage("cannot create newElemSubset[%d]", numSegments);
	//clean up
	return FC_ERROR; 
      }

      rc = fc_addArrayMembersToSubset(elemSegments[numSegments], numcurrsideIDs,
				      currside_elems);
      if (rc != FC_SUCCESS) {
	fc_printfErrorMessage("cannot add ids to facesegment");
	//clean up
	return FC_ERROR; 
      }

      //clean up currside info
      free(currside_faces);
      free(currside_elems);

      numSegments++;
    } //created subsets
  } // while orig_list

  free(in_faceinfolist->vals); //free the array of ptrs, already freed what they are pting to


  //adjmatrix calculation
  intmatrix = (int**)malloc(numSegments*sizeof(int*));
  for (i = 0; i < numSegments; i++){
    intmatrix[i] = (int*)calloc(numSegments,sizeof(int));
  }
  //  printf("internal calc of adjmatrix\n");

  for (i = 0; i < numEdgePairs; i++){
    //      printf("ID %d side %d ID %d side %d\n",edgePairs[i].ID1,sideArray[edgePairs[i].ID1],
    //	    edgePairs[i].ID2,sideArray[edgePairs[i].ID2]);
    if (sideArray[edgePairs[i].ID1] == sideArray[edgePairs[i].ID2]){
      fc_printfWarningMessage("Faces %d %d are edges but both on same side %d - Might be a crease\n",
			      edgePairs[i].ID1,edgePairs[i].ID2,sideArray[edgePairs[i].ID1]);
    }else{
      intmatrix[sideArray[edgePairs[i].ID1]][sideArray[edgePairs[i].ID2]] = 1;
      intmatrix[sideArray[edgePairs[i].ID2]][sideArray[edgePairs[i].ID1]] = 1;
    }
  }

  if (0){
    //print adj matrix
    printf("before leaving shape parts\n");
    for (i = 0; i < numSegments; i++){	
      for (j = 0; j < numSegments; j++){	
	printf("%d ", intmatrix[i][j]);
      }
      printf("\n");
    }
  }

  // cleanup (orig_list is empty)
  free(temp_name);
  free(sideArray);
  free(edgePairs);

  //return values 
  *numSides = numSegments;
  *newFaceSubsets = faceSegments;
  *newElemSubsets = elemSegments;
  *adjmatrix = intmatrix;

  return FC_SUCCESS;
}



/** \name Creating General Shapes from Subsets and Meshes */
//----------------------------------------------------------
//@{

/**
 * \ingroup Shape
 * \brief given a subset, returns an array of FC_Shapes
 *
 * \description 
 *
 *  Given a subset returns an array of FC_Shapes.
 *  Vals for comp_angle are limited to 0->180 degrees because
 *  of the return val of acos. The choice is tricky because 
 *  for curved sides, setting this too low will result in
 *  each face along the curve will be considerd a new side, but
 *  setting it too high means that some things that should be
 *  distinct sides will be combined into one.
 *
 *  The face and elem subsets will be named as the name of the 
 *  original subset, postpended with  "Shape#" followed by
 *  "_FaceSurface#" or "_ElemSurface#" where # will be the index into the
 *  array of subsets for that shape.
 *
 *  In a given shape the indicies of the face and elem arrays correspond,
 *  that is side 0 consists of the faces in the subset face[0] and the
 *  elements in the subset elem[0]. There is no correspondence between
 *  the order of the faces and the order of the elems, within a subset,
 *  that is the first face in side 0 does not necessarily come from
 *  the first elem in side 0. In fact, the subsets may have different
 *  lengths e.g, if one element corresponds to more than one face
 *  on a side (which can happen depending on the angle chosen).
 *
 *  Note that pieces that are touching will come out as
 *  1 distinct shape with a lot of edges as opposed to
 *  two distinct shapes. Segment will do the same thing.
 *
 *  Fails if NULL subset. Returns 0 shapes if empty subset.
 *
 *  Notes:
 *
 *    - this does not use displacements.
 *    - vals for comp_angle are limited to 0->180 degrees because
 *      of the return val of acos. Similary, if the mesh angle
 *      is not in that range, it will return with an error.
 *    - becuase of the calculation of the normals, this method is limited
 *      to elements of topodim 3 who have faces of type \ref FC_ET_TRI or
 *      \ref FC_ET_QUAD
 *    - internally, this uses edge info between sides to build the adj
 *      matrix. In the case of a crease in the shape, two faces may
 *      appear to be on different sides at some point, but later
 *      discovered to be on the same side via neighbors that are
 *      without a crease between them. I believe this all works out
 *      fine, but I dont guarentee it. In the building of the
 *      adjacenecy matrix, the user will be warned about pairs
 *      of faces like this.
 *
 * \todo
 *  - want a hashtable to do the adjmatrix calc (see innards)
 *  - note, this does not clean up memory well in case of weird error
 *     (that shouldnt occur) deep within the list building
 *
 * \modifications
 * - 06/04/06 ACG created
 * - 08/19/06 ACG _fc_getSubsetSides (created 10/04/05 ACG) merged
 *                in now that that facelist processing has been
 *                separated out. (see \ref _fc_buildFaceInfoListFromSubset)
 * - 09/18/06 ACG better subset names
 */
FC_ReturnCode fc_getSubsetShapes(
  FC_Subset in_subset,     /**< input - subset */
  double comp_angle,       /**< input - if two neighboring faces
                               differ by more than this angle then
                               they are considered to be on 
                               different surfaces.*/
  int shared_dim,          /**< input - the minimum dimensionality of shared 
                                part of neighbors of faces */
  int *numShapes,          /**< output - number of shapes */
  FC_Shape **shapes        /**< output shapes array */
){

  FC_ReturnCode rc;
  FC_Mesh mesh;
  FC_AssociationType assoc;

  char *subsetName;
  int numMem, junk;

  _FC_SFI_ARRAY faceinfo_list;

  int numSides;              
  FC_Subset* newFaceSubsets; 
  FC_Subset* newElemSubsets; 
  int **adjmatrix;

  int numParts;
  int *sidesPerPart;
  int **partSideIDS;
  FC_Shape *lshapes;
  
  int i,j,k;

  //defaults
  if (numShapes)
    *numShapes = -1;
  if (shapes)
    *shapes = NULL;

  if (!numShapes || !shapes  || !fc_isSubsetValid(in_subset) ||
  shared_dim < 0 || shared_dim > 1 || comp_angle < 0 ||
  comp_angle > 180){
    fc_printfErrorMessage("%s",fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }


  fc_printfLogMessage("get subset shapes");
  
  // setup and more checking
  rc = fc_getSubsetInfo (in_subset, &numMem, &junk,&assoc);
  if (rc != FC_SUCCESS)
    return rc;

  if (assoc != FC_AT_ELEMENT && assoc != FC_AT_WHOLE_MESH){
    fc_printfErrorMessage("subset must have assoc FC_AT_ELEMENT or FC_AT_WHOLE_MESH");
    return FC_INPUT_ERROR;
  }

  if (numMem == 0){ //return if empty subset
    *numShapes = 0;
    return FC_SUCCESS;
  }


  rc = fc_getMeshFromSubset(in_subset, &mesh);
  if (rc != FC_SUCCESS)
    return rc;

  //let's get all the info we need about the skin
  rc = _fc_buildFaceInfoListFromSubset(in_subset,&faceinfo_list);
  if (rc != FC_SUCCESS){
    fc_printfErrorMessage("cant build face info list\n");
    return rc;
  }

  rc = fc_getSubsetName(in_subset,&subsetName);
  if (rc != FC_SUCCESS){
    _fc_freeSFI_Array(&faceinfo_list);
    free(subsetName);
    return rc;
  }

  rc = _fc_createShapePartsFromFaceInfoList(mesh, &faceinfo_list,&subsetName,
                                            comp_angle,shared_dim,&numSides,
                                            &newFaceSubsets,&newElemSubsets,
					    &adjmatrix);
  if (rc != FC_SUCCESS){
    free(subsetName);
    _fc_freeSFI_Array(&faceinfo_list); //see note about
    //deciding how to handle the freeing of this
  }


  rc = fc_segmentAdjacencyMatrix(numSides,adjmatrix,
				 &numParts, &sidesPerPart, &partSideIDS);
  if (rc != FC_SUCCESS){
    fc_printfErrorMessage("Can't segment adj matrix sides");
    free(subsetName);
    return rc;
  }

  lshapes = (FC_Shape*)malloc(numParts*sizeof(FC_Shape));
  if (lshapes == NULL){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }

  for(i = 0; i < numParts; i++){
    lshapes[i].numSides = sidesPerPart[i];
    lshapes[i].faces = (FC_Subset*)malloc(sidesPerPart[i]*sizeof(FC_Subset));
    lshapes[i].elems = (FC_Subset*)malloc(sidesPerPart[i]*sizeof(FC_Subset));
    lshapes[i].adjmatrix = (int**)malloc(sidesPerPart[i]*sizeof(int*));
    for (j = 0; j < sidesPerPart[i]; j++){
      char* temp_name;
      temp_name = malloc((strlen(subsetName)+40)*sizeof(char));
      if (temp_name == NULL) {
	fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
	return FC_MEMORY_ERROR;
      }
      snprintf(temp_name, 
	       strlen(subsetName)+40,
	       "%s_Shape%d_FaceSurface%d",
	       subsetName, i,j);
      lshapes[i].adjmatrix[j] = (int*)malloc(sidesPerPart[i]*sizeof(int));
      lshapes[i].faces[j] = newFaceSubsets[partSideIDS[i][j]];
      fc_changeSubsetName(lshapes[i].faces[j],temp_name);
      snprintf(temp_name, 
	       strlen(subsetName)+40,
	       "%s_Shape%d_ElemSurface%d", 
	       subsetName,i,j);
      lshapes[i].elems[j] = newElemSubsets[partSideIDS[i][j]];
      fc_changeSubsetName(lshapes[i].elems[j],temp_name);
      free(temp_name);
      for (k = 0; k < sidesPerPart[i]; k++){
	lshapes[i].adjmatrix[j][k] = adjmatrix[partSideIDS[i][j]][partSideIDS[i][k]];
      }
    }
  }


  free(subsetName);
  for (i = 0; i < numSides; i++){
    free(adjmatrix[i]);
  }
  free(adjmatrix);
  free(newFaceSubsets); //just the array
  free(newElemSubsets); //just the array

  for (i = 0 ; i <numParts; i++){
    free(partSideIDS[i]);
  }
  free(partSideIDS);
  free(sidesPerPart);

  *numShapes = numParts;
  *shapes = lshapes;
  
  return FC_SUCCESS;
}


/**
 * \ingroup Shape
 * \brief given a mesh, returns an array of FC_Shapes
 *
 * \description 
 *
 *  given a mesh, returns an array of \ref FC_Shape. This is the
 *  mesh version of \ref fc_getSubsetShapes. See notes there.
 *
 * \modifications
 * - 6/04/06 ACG created
 * - 7/16/06 ACG removed "skin" in name
 *
 */
FC_ReturnCode fc_getMeshShapes(
  FC_Mesh mesh,         /**< input - mesh */
  double comp_angle,           /**< input - if two neighboring faces
                                differ by more than this angle then
                                they are considered to be on 
                                different surfaces.*/
  int shared_dim,              /**< input - the minimum dimensionality of shared 
                                part of neighbors of faces */
  int *numShapes,               /**< output - number of shapes */
  FC_Shape **shapes        /**< output shapes array */
){
  
  FC_ReturnCode rc;
  FC_Subset temp_subset;
  char *meshname;

  //---default return values
  if (numShapes)
    *numShapes = -1;
  if (shapes)
    *shapes = NULL;

  //check input - some checked in sides
  if(!numShapes || !shapes){
    fc_printfErrorMessage("%s",fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // log message
  fc_printfLogMessage("getting mesh shapes");

  //does all error checking except is mesh valid inside _fc_getSubsetShapes
  if (!fc_isMeshValid(mesh)){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // create a temp subset which is the whole mesh & call fc_getSubsetShapes
  rc = fc_getMeshName(mesh,&meshname);
  if (rc != FC_SUCCESS) 
    return rc;

  rc = fc_createSubset(mesh, meshname, FC_AT_WHOLE_MESH, &temp_subset);
  if (rc != FC_SUCCESS) {
     fc_deleteSubset(temp_subset);
     free(meshname);
     return rc;
  }
  free(meshname);

  rc = fc_addMemberToSubset(temp_subset, 0);
  if (rc != FC_SUCCESS) {
     fc_deleteSubset(temp_subset);
     return rc;
  }

   // do it
  rc = fc_getSubsetShapes(temp_subset, comp_angle,shared_dim,numShapes,
			      shapes);
  fc_deleteSubset(temp_subset);
  return rc;
}


/**
 * \ingroup Shape
 * \brief given a shape and a new angle, creates a new shape
 *        where the sides are distinguished byt he new angle
 *
 * \description 
 *
 *        given a shape and a new angle, creates a new shape
 *        where the sides are distinguished by the new angle.
 *        The shared dim is used to determine neighbors for comparison
 *        when you build the sides, but since you are limited
 *        only to those sides that are in the orginal shape, you should
 *        keep this to the orignal number used. (I have no way of checking that)
 *
 *        This function is intended to be used 1) when you get too
 *        few sides or 2) if you have too many sides, as an alternate
 *        crtieria for merging sides than \ref fc_reduceShapeNumSides uses.
 *
 *        The face and elem subsets will be named "Reshape"
 *        postpended with  "Shape#" followed by
 *        "_FaceSurface#" or "_ElemSurface#" where # will be the index into the
 *        array of subsets for that shape.
 *
 * \todo
 *  - check if segdim is a problem or not
 *  - better new shape name
 *       
 * \modifications     
 *  - 08/19/06 ACG creating
 */
FC_ReturnCode fc_reshapeByAngle(
 FC_Shape *inshape, /**< input - orig shape */
 double newangle, /**< input - new angle used for the reshaping */
 int seg_dim, /**< input - seg dim, should be that used for building the orig shape */
 FC_Shape *newshape /**< output - the new shape */
){

  FC_ReturnCode rc;
  FC_Mesh mesh;

  _FC_SFI_ARRAY faceinfo_list;
  int numSides;
  FC_Subset* newFaceSubsets;
  FC_Subset* newElemSubsets;
  int **adjmatrix;
  char* newname = "Reshape";

  // default return - freeing shape is not an option
  if (newshape){
    fc_initShape(newshape);
  }


  if (!inshape || !newshape || newangle < 0 ||
      newangle > 180 || seg_dim < 0 || seg_dim > 2 || 
      !(inshape->numSides)){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  //log message
  fc_printfLogMessage("reshape by angle");

  //assume that since we have a shape already, the topology etc is correct.
  rc = fc_getMeshFromSubset(inshape->faces[0],&mesh);
  if (rc != FC_SUCCESS){
    fc_printfErrorMessage("cant get mesh for shape");
    return rc;
  }

  rc = _fc_buildFaceInfoListFromShape(inshape,&faceinfo_list);
  if (rc != FC_SUCCESS){
    fc_printfErrorMessage("can't build face info list\n");
    return rc;
  }

  rc = _fc_createShapePartsFromFaceInfoList(mesh,&faceinfo_list,
					    &newname,newangle,
					    seg_dim,
					    &numSides,
					    &newFaceSubsets,
					    &newElemSubsets,
					    &adjmatrix);
  if (rc != FC_SUCCESS){
    fc_printfErrorMessage("can't create shape parts from face info list\n");
    _fc_freeSFI_Array(&faceinfo_list);
    return rc;
  }

  //unlike subsetshapes we want to return this as one shape
  newshape->numSides = numSides;
  newshape->faces = newFaceSubsets;
  newshape->elems = newElemSubsets;
  newshape->adjmatrix = adjmatrix; 

  //dont free anything since theyve been assigned to the shape
  return FC_SUCCESS;
}



/**
 * \ingroup Shape
 * \brief given a shape and a desired number of sides, fewer
 *        than the current number of sides, creates a new shape
 *        with the desired number of sides
 *
 * \description 
 *
 *        Given a shape and a desired number of sides, fewer
 *        than the current number of sides, creates a new shape
 *        with the desired number of sides. This function is 
 *        intended to be used to merge small, perhaps curved
 *        faces into a major side. For instance, if on a
 *        rectangular prism with curved edges you set the
 *        angle to determine the division between edges as
 *        large, then it may resolve into one continuous side;
 *        on the other hand, if you set it to small, then all 
 *        faces along the curve may resolve into distinct sides.
 *        If you instead know that you want to end with 6 sides,
 *        then call shape with a small angle, resulting in too many
 *        sides, and then call this function to reduce the shape 
 *        to a 6 sided shape, retaining the 6 largest sides.
 *
 *        The innards iterate assigning each reduced out
 *        side to an arbitrary side of distance 1 as given 
 *        by the adj matrix. The reducedout faces therefore arent
 *        elegantly distributed into the new sides, and you
 *        may end up with odd adjacency matricies as a result.
 *
 *        The kept sides retain their same names and ordering.
 *
 *        See \ref fc_reshapeByAngle
 *
 * \todo
 *  - avoid the shape copy at the end
 *  - intent is to later have another function that class this one that
 *    does the shape call and somehow decides what the ideal
 *    number and IDs of sides is. this could be done for basic shapes
 *    if there is a great discrepancy in the areas of the resulting
 *    sides - e.g., 6 large sides and a zillion small ones
 *       
 * \modifications     
 * - 08/09/06 ACG creating
 * - 09/27/06 ACG changing algorithm
 */
FC_ReturnCode fc_reduceShapeNumSides(
  FC_Shape *inshape, /**< input - shape to be reduced */
  int reducedNumSides, /**< input - number of sides to be in reduced shape */
  FC_Shape *reducedshape /**< output - the reduced shape */
  ){

  FC_ReturnCode rc, rc1, rc2;
  FC_Mesh mesh;
  int **distmatrix;

  //dont need to be sorted, just taking advantage of the selfexpansion
  FC_SortedIntArray *combosides, *newlyassignedsides;
  FC_SortedIntArray remainingsides;
  FC_Shape newshape;
  int i,j,k,l;

  //defaults
  if (reducedshape){
    fc_initShape(reducedshape);
  }

  if (!inshape || reducedNumSides < 1 || 
      reducedNumSides >= inshape->numSides ||
      !reducedshape){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  //log message
  fc_printfLogMessage("reduceShapeNumSides");

  rc = fc_getMeshFromSubset(inshape->faces[0],&mesh);
  if (rc != FC_SUCCESS){
    fc_printfErrorMessage("Cant get mesh from subset");
    return rc;
  }

  //first copy the sides that we are keeping over
  combosides = (FC_SortedIntArray*)malloc(reducedNumSides*
					  sizeof(FC_SortedIntArray));
  newlyassignedsides = (FC_SortedIntArray*)malloc(reducedNumSides*
					  sizeof(FC_SortedIntArray));
  if (!combosides || ! newlyassignedsides){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }

  for (i = 0; i < reducedNumSides; i++){
    rc1 = fc_initSortedIntArray(&combosides[i]);
    rc2 = fc_initSortedIntArray(&newlyassignedsides[i]);
    if (rc1 !=FC_SUCCESS || rc2 != FC_SUCCESS){
      fc_printfErrorMessage("Can't init SIA");
      for (j = 0; j < i; j++){
	fc_freeSortedIntArray(&combosides[j]);
	fc_freeSortedIntArray(&newlyassignedsides[j]);
      }
      free(combosides);
      free(newlyassignedsides);
      return (rc1 != FC_SUCCESS ? rc1 : rc2);
    }
    rc = fc_addIntToSortedIntArray(&combosides[i],i);
    if (rc < 0){
      fc_printfErrorMessage("Can't add int to SIA");
      for (j = 0; j < i; j++){
	fc_freeSortedIntArray(&combosides[j]);
	fc_freeSortedIntArray(&newlyassignedsides[j]);
      }
      free(combosides);
      free(newlyassignedsides);
      return rc;
    }
  }

  //get ids for left over faces
  rc = fc_initSortedIntArray(&remainingsides);
  if (rc !=FC_SUCCESS){
    fc_printfErrorMessage("Cant init sia");
    for (j = 0; j < reducedNumSides; j++){
      fc_freeSortedIntArray(&combosides[j]);   
      fc_freeSortedIntArray(&newlyassignedsides[j]);
    }
    free(combosides);
    free(newlyassignedsides);
    return rc;
  }

  for (i = reducedNumSides; i < inshape->numSides; i++){
    rc = fc_addIntToSortedIntArray(&remainingsides,i);
    if (rc < 0){
      fc_printfErrorMessage("cant add int to SIA");
      for (j = 0; j < reducedNumSides; j++){
	fc_freeSortedIntArray(&combosides[j]);
	fc_freeSortedIntArray(&newlyassignedsides[j]);
      }
      free(combosides);
      free(newlyassignedsides);
      fc_freeSortedIntArray(&remainingsides);
      return rc;
    }
  }

  rc = fc_calcDistanceMatrix(inshape->numSides, inshape->adjmatrix,
			     &distmatrix);
  if (rc != FC_SUCCESS){
    fc_printfErrorMessage("Can't get dist matrix");
    for (j = 0; j < reducedNumSides; j++){
      fc_freeSortedIntArray(&combosides[j]);
      fc_freeSortedIntArray(&newlyassignedsides[j]);
    }
    free(combosides);
    free(newlyassignedsides);
    fc_freeSortedIntArray(&remainingsides);
    return rc;
  }

  while(remainingsides.numVal){
    int assigned = 0;
    
    //find all sides that are distance 1 from one ofthe established sides
    //only update the established sides after each pass
    for (i = 0; i < remainingsides.numVal; i++){
      int foundside = 0;
      int comboid;
      int currside = remainingsides.vals[i];
      for (j = 0; j < reducedNumSides; j++){
	for (k = 0; k < combosides[j].numVal; k++){
	  int compside = combosides[j].vals[k];
	  comboid = j;
	  if (distmatrix[currside][compside] ==1){
	    foundside = 1;
	    break;
	  }
	}
	if (foundside)
	  break;
      }
      if (foundside){
	rc1 = fc_addIntToSortedIntArray(&newlyassignedsides[j],currside);
	rc2 = fc_deleteIntFromSortedIntArray(&remainingsides,currside);
	if (rc1 < 0 || rc2 < 0){
	  fc_printfErrorMessage("cant add/delete int to/from SIA");
	  for (j = 0; j < reducedNumSides; j++){
	    fc_freeSortedIntArray(&combosides[j]);
	    fc_freeSortedIntArray(&newlyassignedsides[j]);
	  }
	  free(combosides);
	  free(newlyassignedsides);
	  fc_freeSortedIntArray(&remainingsides);
	  for (j = 0; j < inshape->numSides; j++){
	    free(distmatrix[j]);
	  }
	  free(distmatrix);
	  return (rc1 < 0 ? rc1 : rc2);
	}
	assigned++;
      }
    }
    
    if (!assigned){
      fc_printfErrorMessage("None of the remaining sides are distance 1 from any of the current sides");
      for (j = 0; j < reducedNumSides; j++){
	fc_freeSortedIntArray(&combosides[j]);
	fc_freeSortedIntArray(&newlyassignedsides[j]);
      }
      free(combosides);
      free(newlyassignedsides);
      fc_freeSortedIntArray(&remainingsides);
      for (j = 0; j < inshape->numSides; j++){
	free(distmatrix[j]);
      }
      free(distmatrix);
      return FC_INPUT_ERROR;
    }

    //copy the newly assigned sides into the established sides
    for (i = 0; i < reducedNumSides; i++){
      for (j = 0; j < newlyassignedsides[i].numVal; j++){
	rc = fc_addIntToSortedIntArray(&combosides[i],newlyassignedsides[i].vals[j]);
	if (rc < 0){
	  fc_printfErrorMessage("cant add int to SIA");
	  for (j = 0; j < reducedNumSides; j++){
	    fc_freeSortedIntArray(&combosides[j]);
	    fc_freeSortedIntArray(&newlyassignedsides[j]);
	  }
	  free(combosides);
	  free(newlyassignedsides);
	  fc_freeSortedIntArray(&remainingsides);
	  for (j = 0; j < inshape->numSides; j++){
	    free(distmatrix[j]);
	  }
	  free(distmatrix);
	  return rc;
	}
      }
      fc_freeSortedIntArray(&newlyassignedsides[i]);
    }
  } //while

  //clean up some things
  for (j = 0; j < inshape->numSides; j++){
    free(distmatrix[j]);
  }
  free(distmatrix);
  free(newlyassignedsides); //individual ones already freed
  fc_freeSortedIntArray(&remainingsides);

  //now build the new shape
  newshape.numSides = reducedNumSides;
  newshape.faces = (FC_Subset*)malloc(reducedNumSides*sizeof(FC_Subset));
  newshape.elems = (FC_Subset*)malloc(reducedNumSides*sizeof(FC_Subset));
  newshape.adjmatrix = (int**)malloc(reducedNumSides*sizeof(int*));
  if (!newshape.faces || !newshape.elems || !newshape.adjmatrix){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }

  for (i = 0; i < reducedNumSides; i++){
    //allocate them all now so can use free shape as a set later
    newshape.adjmatrix[i] = (int*)malloc(reducedNumSides*sizeof(int));
    if (!newshape.adjmatrix[i]){
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      return FC_MEMORY_ERROR;
    }
    for (j = 0; j < reducedNumSides; j++){
      newshape.adjmatrix[i][j] = 0;
    }
  }

  for (i = 0; i < reducedNumSides; i++){
    int nume, numf;
    int *arrf, *arre;
    char *temp;

    //assign the first item in each combo side to the new side
    rc1 = fc_copySubset(inshape->faces[combosides[i].vals[0]],mesh,
			"temp",&(newshape.faces[i]));
    rc2 = fc_copySubset(inshape->elems[combosides[i].vals[0]],mesh,
			"temp",&(newshape.elems[i]));
    if (rc1 != FC_SUCCESS || rc2 != FC_SUCCESS){
      fc_printfErrorMessage("Cant copy subsets");
      for (j = 0; j <= i; j++){
	fc_freeSortedIntArray(&combosides[j]);
      }
      free(combosides);
      fc_freeShape(&newshape);
      return (rc1 != FC_SUCCESS ? rc1:rc2);
    }
    
    //rename them to what they were in that original position
    //if cant rename, just warn
    rc = fc_getSubsetName(inshape->faces[i],&temp);
    if (rc != FC_SUCCESS || temp == NULL){
      fc_printfWarningMessage("fc_reduceShapeNumSides can't get name for face side %d", i);
    }else{
      rc = fc_changeSubsetName(newshape.faces[i],temp);
      if (rc!= FC_SUCCESS){
	fc_printfWarningMessage("fc_reduceShapeNumSides can't rename face side %d", i);
      }
    }
    if (temp) free(temp);
	  
    rc = fc_getSubsetName(inshape->elems[i],&temp);
    if (rc!= FC_SUCCESS || temp == NULL){
      fc_printfWarningMessage("fc_reduceShapeNumSides can't get name for elem side %d", i);
    }else{
      rc = fc_changeSubsetName(newshape.elems[i],temp);
      if (rc!= FC_SUCCESS){
	fc_printfWarningMessage("fc_reduceShapeNumSides can't rename elem side %d", i);
      }
    }
    if (temp) free(temp);

    //now copy all other members over
    for (j = 1; j < combosides[i].numVal; j++){
      int currside = combosides[i].vals[j];
      rc1 = fc_getSubsetMembersAsArray(inshape->faces[currside],
				       &numf,&arrf);
      rc2 = fc_getSubsetMembersAsArray(inshape->elems[currside],
				       &nume,&arre);
      if (rc1 != FC_SUCCESS || rc2 != FC_SUCCESS){
	fc_printfErrorMessage("Can't get subset members");
	for (j = 0; j < reducedNumSides; j++){
	  fc_freeSortedIntArray(&combosides[j]);
	}
	free(combosides);
	fc_freeShape(&newshape);
	if (arrf) free (arrf);
	if (arre) free (arre);
	return (rc1 != FC_SUCCESS ? rc1 : rc2);
      }

      //make sure wont fail if add duplicates
      rc1 = fc_addArrayMembersToSubset(newshape.faces[i],
				       numf,arrf);
      rc1 = fc_addArrayMembersToSubset(newshape.elems[i],
				       nume,arre);
      free(arrf);
      free(arre);
      if (rc1 != FC_SUCCESS || rc2 != FC_SUCCESS){
	fc_printfErrorMessage("Can't add subset members");
	for (j = 0; j < reducedNumSides; j++){
	  fc_freeSortedIntArray(&combosides[j]);
	}
	free(combosides);
	fc_freeShape(&newshape);
	return (rc1 != FC_SUCCESS ? rc1 : rc2);
      }
    } //combosides.numVal
  } // reducedNumSides


  //now all the sides have been assigned and the faces and the elems put in
  //the adj matrix is all zeros at this point
  for (i = 0; i < reducedNumSides; i++){
    for (j = i+1; j < reducedNumSides; j++){
      if (inshape->adjmatrix[i][j]){
	//if they were adj before, they still are
	newshape.adjmatrix[i][j] = 1;
	newshape.adjmatrix[j][i] = 1;
      }else{
	//they may now be adj becuase of the reducing
	for (k = 0; k < combosides[i].numVal; k++){
	  for (l = 0; l < combosides[j].numVal; l++){
	    int sideA = combosides[i].vals[k];
	    int sideB = combosides[j].vals[l];
	    if (inshape->adjmatrix[sideA][sideB]){
	      newshape.adjmatrix[i][j] = 1;
	      newshape.adjmatrix[j][i] = 1;
	      break;
	    }
	  }
	  if (newshape.adjmatrix[i][j]){
	    //we found a match inside
	    break;
	  }
	}
      }
    }
  }

  //cleanup
  for (i = 0; i < reducedNumSides; i++){
    fc_freeSortedIntArray(&combosides[i]);
  }
  free(combosides);

  //handle this shape better
  rc = fc_copyShape(&newshape,reducedshape);
  fc_freeShape(&newshape);

  return FC_SUCCESS;
}


//@}


/** \name Shape Helpers */
//-------------------------------
//@{

/**
 * \ingroup Shape
 * \brief Copy a shape.
 *
 * \description
 *
 *    This is a deep copy; it allocates memory and the created shape will need
 *    to be freed by the caller.
 *
 * \modifications
 *  - 06/18/06 ACG created
 */
FC_ReturnCode fc_copyShape(
  FC_Shape *inshape, /**< input - ptr to shape */
  FC_Shape *outshape /**< output - ptr to the shape */
){

  FC_ReturnCode rc;
  FC_Shape shape;
  int i,j;

  //defaults
  if (outshape){
    fc_initShape(outshape);
  }


  if (!inshape || !outshape ){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // log message
  fc_printfLogMessage("copying shape");

  shape.numSides = inshape->numSides;
  shape.faces = (FC_Subset*)malloc(shape.numSides*sizeof(FC_Subset));
  shape.elems = (FC_Subset*)malloc(shape.numSides*sizeof(FC_Subset));
  shape.adjmatrix = (int**)malloc(shape.numSides*sizeof(int*));
  if (!shape.faces || ! shape.elems || ! shape.adjmatrix){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  for (i = 0; i < shape.numSides; i++){
    FC_Mesh mesh; //these really should all be onthe same mesh
    char* subname;

    rc = fc_getMeshFromSubset(inshape->faces[i],&mesh);
    if (rc != FC_SUCCESS){
      fc_printfErrorMessage("cant get mesh from subset");
      fc_freeShape(&shape);
      return rc;
    }
    rc = fc_getSubsetName(inshape->faces[i],&subname);
    if (rc != FC_SUCCESS){
      fc_printfErrorMessage("cant get name of subset");
      fc_freeShape(&shape);
      return rc;
    }
    //is it ok to use the same name ?
    rc =  fc_copySubset(inshape->faces[i],mesh,subname,&(shape.faces[i]));
    if (rc != FC_SUCCESS){
      fc_printfErrorMessage("cant copy subset");
      fc_freeShape(&shape);
      return rc;
    }
    free(subname);

    rc = fc_getMeshFromSubset(inshape->elems[i],&mesh);
    if (rc != FC_SUCCESS){
      fc_printfErrorMessage("cant get mesh from subset");
      fc_freeShape(&shape);
      return rc;
    }
    rc = fc_getSubsetName(inshape->elems[i],&subname);
    if (rc != FC_SUCCESS){
      fc_printfErrorMessage("cant get name of subset");
      fc_freeShape(&shape);
      return rc;
    }
    //is it ok to use the same name ?
    rc =  fc_copySubset(inshape->elems[i],mesh,subname,&(shape.elems[i]));
    if (rc != FC_SUCCESS){
      fc_printfErrorMessage("cant copy subset");
      fc_freeShape(&shape);
      return rc;
    }
    
    free(subname);

    shape.adjmatrix[i] = (int*)malloc(shape.numSides*sizeof(int));
    if (!shape.adjmatrix[i]){
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      return FC_MEMORY_ERROR;
    }
  }
    
  for (i = 0; i < shape.numSides; i++){
    for (j = 0; j < shape.numSides; j++){
      shape.adjmatrix[i][j] = inshape->adjmatrix[i][j];
    }
  }

  outshape->numSides = shape.numSides;
  outshape->faces = shape.faces;
  outshape->elems = shape.elems;
  outshape->adjmatrix = shape.adjmatrix;

  return FC_SUCCESS;
}


/**
 * \ingroup Shape
 * \brief given a new ordering of the sides, reorders a shape 
 *
 * \description
 *
 *  Given a new ordering of the sides, reorders a shape.
 *  Example: given the orig sides (set0,set1,set2,set3,set4) with
 *  new order = (2,4,1,0,3) then the sides will be
 *  reordered to (set2,set4,set1,set0,set3). That is
 *  what was in place 2 has been moved to place 0, as
 *  opposed to moving what was in place 0 to place 2.
 *
 *  This keeps the original names of the subsets to allow
 *  referring to subsets by name to still be consistent. If you
 *  want to reorder the names, you must do so manually.
 *
 * \modifications
 *  - 06/18/06 ACG created
 *  - 06/19/06 ACG changed reorder so the indicies
 *    are now what the new order is going to be, rather
 *    than where the old order is moving to, because
 *    its more intuitive that way.
 */
FC_ReturnCode fc_reorderShape(
  int *neworder, /**< input - array of new order of indicies */
  FC_Shape *shape /**< output - ptr to the shape */
){

  FC_Subset *temp, *temp2; 
  int numSides;
  int** adjmatrix;
  int i,j;

  //no defaults

  if (!shape || !neworder){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // log message
  fc_printfLogMessage("reordering shape");

  numSides = shape->numSides;

  //make sure neworder array makes sense
  int *test = calloc(numSides, sizeof(int));
  if (test == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }

  for (i = 0; i < numSides; i++){
    if (neworder[i] < 0 || neworder[i] >= numSides){
      free(test);
      return FC_INPUT_ERROR;
    }
    test[neworder[i]] = 1;
  }

  for (i = 0; i < numSides; i++){
    if (!test[i]){
      free(test);
      return FC_INPUT_ERROR;
    }
  }
  
  free(test);

  //doesnt swap ptrs to subsets, since the subsets themselves are just handles.
  temp = (FC_Subset*)malloc(numSides*sizeof(FC_Subset));
  temp2 = (FC_Subset*)malloc(numSides*sizeof(FC_Subset));
  adjmatrix = (int**)malloc(numSides*sizeof(int*));
  if (temp == NULL || temp2 == NULL || adjmatrix == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  for (i = 0; i < numSides; i++){
    temp[i] = shape->faces[neworder[i]];
    temp2[i] = shape->elems[neworder[i]];
    adjmatrix[i] = (int*)malloc(numSides*sizeof(int));
    if (adjmatrix[i] == NULL){
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      return FC_MEMORY_ERROR;
    }
  }

  free(shape->faces);
  free(shape->elems);
  shape->faces = temp;
  shape->elems = temp2;

  for (i =0; i < numSides; i++){
    for (j =0; j < numSides; j++){
      adjmatrix[i][j] = shape->adjmatrix[neworder[i]][neworder[j]];
    }
  }
  for (i =0; i < numSides; i++){
    free(shape->adjmatrix[i]);
  }
  free(shape->adjmatrix);

  shape->adjmatrix = adjmatrix;
  
  return FC_SUCCESS;
}


/**
 * \ingroup Shape
 * \brief frees a Shape
 *
 * \description 
 *
 *    Deallocate contents of an \ref FC_Shape. 
 *
 * \modifications
 * - 06/04/06 ACG created
 * - 07/13/06 ACG renamed from freeShape
 * - 07/16/06 removed "skin" from name
 */
FC_ReturnCode fc_freeShape(
  FC_Shape *fcss /**< input - the shape */
){
  int i;
  FC_ReturnCode rc;

  if (!fcss){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  for (i = 0; i < fcss->numSides; i++){
    rc = fc_deleteSubset(fcss->faces[i]);
    if (rc != FC_SUCCESS){
      fc_printfErrorMessage("Can't free face %d",i);
      return rc;
    }
    rc = fc_deleteSubset(fcss->elems[i]);
    if (rc != FC_SUCCESS){
      fc_printfErrorMessage("Can't free elems %d",i);
      return rc;
    }
    free(fcss->adjmatrix[i]);
  }
  free(fcss->faces);
  free(fcss->elems);
  free(fcss->adjmatrix);
  fcss->faces= NULL;
  fcss->elems= NULL;
  fcss->adjmatrix= NULL;
  fcss->numSides = 0;

  return FC_SUCCESS;
}


/**
 * \ingroup Shape
 * \brief clears the allocated shape memory, but leaves the
 *        subsets alive
 *
 * \description 
 *
 *    Clears the allocated shape memory, but leaves
 *    the subsets alive. The user will need handles to
 *    them to eventually delete them.
 *
 * \modifications
 * - 08/22/06 ACG created
 */
FC_ReturnCode fc_clearShape(
  FC_Shape *fcss /**< input - the shape */
){
  int i;

  if (!fcss){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  for (i = 0; i < fcss->numSides; i++){
    free(fcss->adjmatrix[i]);
  }
  free(fcss->faces);
  free(fcss->elems);
  free(fcss->adjmatrix);
  fcss->faces= NULL;
  fcss->elems= NULL;
  fcss->adjmatrix= NULL;
  fcss->numSides = 0;

  return FC_SUCCESS;
}


/**
 * \ingroup Shape
 * \brief inits a Shape
 *
 * \description 
 *
 *  Blanks out a shape's fields w/o doing deallocation.
 *  This is used in cases where \ref fc_freeShape is
 *  not feasible (e.g., can't just free a newly allocated
 *  shape becuase the numSides may be anything).
 *
 *  In most cases the user wont need to init a shape.
 *  The shape creation methods will handle the newly
 *  instantiated shape pointer properly and after that
 *  the user will be working with the created shape.
 *
 * \modifications
 * - 09/28/06 ACG created
 */
FC_ReturnCode fc_initShape(
  FC_Shape *fcss /**< input - the shape */
){
  if (!fcss){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  fcss->faces= NULL;
  fcss->elems= NULL;
  fcss->adjmatrix= NULL;
  fcss->numSides = 0;

  return FC_SUCCESS;
}


/**
 * \ingroup Shape
 *
 * \brief prints a Shape
 *
 * \description
 *
 *  prints an \ref FC_Shape
 *
 * \todo 
 *  - test
 *
 * \modifications
 * - 06/05/06 ACG created
 * - 07/13/06 ACG dropped 'FC' in name
 * - 07/16/06 ACG dropped "skin" in name
 * - 07/18/06 ACG shape are is ptr
 */
FC_ReturnCode fc_printShape(
  FC_Shape *fcss, /**< input - the shape */
  char* label,
  int verb
){
  FC_ReturnCode rc;
  int i,j;

  if (!fcss){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  if (label)
    printf(label);

  printf("Shape has %d sides\n",fcss->numSides);
  if (verb){
    printf("Faces:\n");
    for (i = 0; i < fcss->numSides; i++){
      rc = fc_printSubset(fcss->faces[i],"",verb);
    }
    printf("Elems:\n");
    for (i = 0; i < fcss->numSides; i++){
      rc = fc_printSubset(fcss->elems[i],"",verb);
    }
  }

  printf("AdjMatrix:\n");
  for (i = 0; i < fcss->numSides; i++){
    for (j = 0; j < fcss->numSides; j++){
      printf("%d ",fcss->adjmatrix[i][j]);
    }
    printf("\n");
  }
  printf("\n");

  return FC_SUCCESS;
}
//@}


/** \name Side Orders */
//----------------------------------------------------------
//@{
/**
 * \ingroup Shape
 * \brief given a shape, returns the reordering array for a screw
 *
 * \description 
 *
 *    Given a shape, returns the reordering array for a screw.
 *    It returns the sides in the following order:
 *    head plane, head cylinder, middle plane, base cylinder, base plane
 *    (I've made up these names)
 *
 *    If the shape does not satisfy the screw criteria, the method
 *    returns FC_SUCCESS but a NULL order array. These cases are
 *    shapes without 5 sides or those whose adjacency 
 *    matricies don't fit the screw adjmatrix pattern. Some non-screw
 *    things in theory could slip through the testing cracks. It is up
 *    to the user to make sure he uses this in the right way.
 *
 *
 * \modifications
 *  - 06/18/06 ACG created
 *  - 07/19/06 ACG changed so if the shape is not a screw it returns 
 *                 FC_SUCCESS but an empty order array. Only returns
 *                 FC_ERROR if there is an actual error.
 *  - 08/06/06 ACG more of the above - for if not 5 sided figure.
 */
FC_ReturnCode fc_createScrewShapeOrder(
 FC_Shape *shape, /**< input - shape */
 int **neworder       /**< output - order array */
 ){

  //head plane has 1 adj side
  //base plane has 1 adj side
  //head cylinder is adj to the head plane and the mid plane
  //base cylinder is adj to the base plane and the mid plane
  //mid plane is adj to the base cylinder and the head cylinder
  //head plane has greater area than base plane

  FC_ReturnCode rc;

  int endcount = 0;
  int othercount = 0;
  int endoptions[5], otheroptions[5];
  double areas[2];
  int *order;
  int nSides;
  int i,j;



  //defaults
  if (neworder)
    *neworder = NULL;

  if (!shape || !neworder){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // log message
  fc_printfLogMessage("Getting screwshape order");

  nSides = shape->numSides;
  if (nSides != 5){
    //success but null return
    return FC_SUCCESS;
  }

  order = (int*)malloc(5*sizeof(int));
  if (!order){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }

  for (i = 0; i < nSides; i++){
    int count = 0;

    //can initialize these here because endcount and othercount
    //cant exceed i
    order[i] = -1;
    endoptions[i] = -1;
    otheroptions[i] = -1;

    for (j = 0; j < nSides; j++){ 
      //      printf(" %d",shape->adjmatrix[i][j]);
      //check matrix for symmetry
      if (shape->adjmatrix[i][j] != shape->adjmatrix[j][i]){
	fc_printfErrorMessage("adj matrix not symmetric");
	free(order);
	return FC_INPUT_ERROR; 
      }
      if (shape->adjmatrix[i][j]) count++; //count the adjacencies
    }
    //    printf("\n");
    if (count == 1){
      //its the base plane or top plane
      endoptions[endcount++] = i;
      //      printf("end - %d\n",i);
    }else if (count == 2){
      //its the top cylinder, base cylinder, or middle plane
      otheroptions[othercount++] = i;
      //      printf("other - %d\n",i);
    }else {
      //this is an invalid shape
      //      fc_printfErrorMessage("Invalid screw shape");
      free(order);
      return FC_SUCCESS; //not a screw, null return
    }
  }
      
  if (endcount != 2 || othercount != 3){
    //    fc_printfErrorMessage("Invalid screw shape");
    free(order);
    return FC_SUCCESS; // not a screw
  }

  //determine ends
  for(i = 0; i < endcount; i++){ //endcount is 2
    rc = fc_getSubsetArea(shape->faces[endoptions[i]],&(areas[i]));
    if (rc !=FC_SUCCESS){
      fc_printfErrorMessage("cant get areas");
      free(order);
      return rc;
    }
  }

  i = (areas[0] > areas[1] ? 0:1);
  order[0] = endoptions[i];
  order[4] = endoptions[!i];


  //determine cylinders and midplane
  for (i = 0; i < othercount; i++){
    if (shape->adjmatrix[otheroptions[i]][order[0]] == 1){
      order[1] = otheroptions[i];
    }else if (shape->adjmatrix[otheroptions[i]][order[4]] == 1){
      order[3] = otheroptions[i];
    }else {
      order[2] = otheroptions[i];
    }
  }

  //check midplane adjacencies
  if ((shape->adjmatrix[order[2]][order[1]] != 1) ||
      shape->adjmatrix[order[2]][order[3]] != 1){
    //    fc_printfErrorMessage("invalid screw shape");
    free(order);
    return FC_SUCCESS; //not a screw, null return
  }

  *neworder = order;

  return FC_SUCCESS;
}




/** 
 * \ingroup  Shape
 * \brief given a shape returns an array of the side ids
 *        in descending face area order
 *
 * \description
 * 
 *        Given a shape returns an array of the side ids
 *        in descending face area order. This can then be fed into
 *        \ref fc_reorderShape
 *
 * \modifications 
 *   - 06/18/06 ACG, created
 *   - 06/19/06 ACG updated to reflect change in the
 *     ordering of the array expected by reorderShape, unittested
 *   - 09/27/06 ACG changed from ascending to descending 
 */
FC_ReturnCode fc_createDescendingAreaOrder(
  FC_Shape *shape, /**< input - shape */
  int **neworder       /**< output - ids of sides */
  ) {

  FC_ReturnCode rc;
  double *areas;
  int numSides;
  int *order;
  int i,j;

  //defaults
  if (neworder)
    *neworder = NULL;

  if ( !shape ||  !neworder){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // log message
  fc_printfLogMessage("Getting descending area order");

  numSides = shape->numSides;

  rc= fc_getShapeSidesAreas(shape,&areas);
  if (rc != FC_SUCCESS){
    fc_printfErrorMessage("cant get sides areas");
    return rc;
  }
    
  order = (int*)malloc(numSides*sizeof(int));
  if (order == NULL){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  
  for (i =0; i < numSides; i++){
    order[i] = i;
  }
  

  //order sides by area
  //sort sides areas. yes i know, its bubblesort.
  {
    double tempa;
    int tempi;
    int test; 
    
    for(i = numSides - 1; i > 0; i--){
      test=0;
      for(j = 0; j < i; j++){
	if (areas[j] > areas[j+1]){
	  tempa = areas[j]; 
	  areas[j] = areas[j+1];
	  areas[j+1] = tempa;
	  tempi = order[j]; 
	  order[j] = order[j+1];
	  order[j+1] = tempi;
	  test=1;
	}
      }
      if(test==0) break; //exit if sorted
    }
  }
  free(areas);

  //reverse the order
  //flip order around for convenience
  for (i = 0; i < numSides/2; i++){
    int a = order[i];
    order[i] = order[numSides-1-i];
    order[numSides-1-i] = a;
  }

  *neworder = order;

  return FC_SUCCESS;
}


/** 
 * \ingroup  Shape
 * \brief As an approximation to obtaining the major
 *        sides of a thin shape, given a shape returns
 *        an array of the side ids in the ordered desired
 *        for shape with a large side and a non-adjacent side,
 *        subject to certain constraints.
 *
 * \description
 *        As an approximation to obtaining the major
 *        sides of a thin shape, given a shape returns
 *        an array of the side ids in the ordered desired
 *        for shape with a large side and a non-adjacent side,
 *        subject to certain constraints.
 * 
 *        A thin shape is a shape that in some 
 *        dimension is narrow and in a roughly perpendicular
 *        direction has a pair of large opposing sides,
 *        which are the major sides that the user is 
 *        interested in. The thin shape order would then have
 *        those major sides in the 0th and 1st
 *        position, with all other sides in random order.
 *        This can then be fed into \ref fc_reorderShape.
 *
 *        This method uses simplifying assumtions to get the
 *        major sides, which may not be correct in all cases.
 *        Specifically, this method attempts to find the largest
 *        side and its largest non-adjacent side and call
 *        them the major sides subject to the following
 *        constraints (in order):
 *        - If there are three or more sides that are the largest,
 *          it is not a thin shape (e.g., rectangular prism)
 *        - If the two largest sides are of equal area, returns them
 *          in random order if they are nonadj. If they are adjacent,
 *          then it is not a thin shape.
 *        - If the largest non-adj side to the unique largest side
 *          is equal in size to another non-adj side, then it is not a thin shape.
 *        - If there is no side that is non-adjacent to the largest
 *          side then it is not a thin shape.
 *
 *        Note that tests for equals are FC_FLT_EQUIV in order
 *        to keep the interface small. If you want to change
 *        those conditions, you have to write your own method.
 * 
 *        Note that there can be undesirable results for shapes
 *        that have no sharp corners and therefore what may be
 *        expected to be the thin side actually wraps around the shape
 *        as one large continuous side and is therefore the largest side.
 *        This case is attempted to be handled by
 *        \ref fc_createLargeAndOpposingSidesOrder.
 *
 *        Returns with an error if there is something
 *        fundamentally wrong. Returns FC_SUCCESS but empty
 *        order array if the shape does not satisfy the
 *        thin shape criteria.
 *
 * \todo
 *  - allows user to send in criteria for equals comparison?
 *
 *
 * \modifications 
 *   - 07/18/06 ACG created
 *   - 07/19/06 ACG changed so returns FC_SUCCESS and null order if it
 *                  is not an appropriate shape. Only error if there
 *                  is actually an error. 
 *   - 08/08/06 ACG renamed from createThinShapeOrder.
 *   - 09/19/06 ACG changed equal comparison from dbl to float
 */
FC_ReturnCode fc_createLargeAndNonAdjacentSidesOrder(
  FC_Shape *shape, /**< input - shape */
  int **neworder       /**< output - ids of sides */
  ) {

  FC_ReturnCode rc;
  double *areas;
  int numSides;
  int *order;
  int sidemax,nextside,testside;
  int i,count;

  //defaults
  if (neworder)
    *neworder = NULL;

  if (!shape || !neworder){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // log message
  fc_printfLogMessage("Getting large nonadjside area order");

  numSides = shape->numSides;

  //must be at least 3 sides for this to do anything
  if (numSides < 3){
    fc_printfErrorMessage("shape must have at least 3 sides to have a nonadj side");
    //    return FC_INPUT_ERROR;
    return FC_SUCCESS; //null return
  }

  rc = fc_getShapeSidesAreas(shape,&areas);
  if (rc != FC_SUCCESS){
    fc_printfErrorMessage("cant get sides areas");
    return rc;
  }

  rc = fc_createDescendingAreaOrder(shape,&order);
  if (rc != FC_SUCCESS){
    fc_printfErrorMessage("cant get descending order");
    free(areas);
    return rc;
  }

  sidemax = -1;
  nextside = -1;
  testside = -1;

  if (FC_FLT_EQUIV(areas[order[0]], areas[order[1]])){ //two largest sides are of equal area
    if (FC_FLT_EQUIV(areas[order[1]], areas[order[2]])){
      fc_printfErrorMessage("3 or more largest sides - cant get thin shape order");
      free(areas);
      free(order);
      return FC_SUCCESS; // not a thin shape, null return
    }else{
      //only two largest sides but if the two largest sides are adj - bad
      if (shape->adjmatrix[order[0]][order[1]] == 1){
	fc_printfErrorMessage("two largest sides adj - cant get thin shape order");
	free(areas);
	free(order);
	return FC_SUCCESS; // not a thin shape, null return
      }else{
	//return them in random order
	sidemax = order[0];
	nextside = order[1];
      }
    }
  }else{
    //two largest sides are not of equal area

    sidemax = order[0];
    //find the largest non-adjacent side. if that isnt unique - bad

    nextside = -1;
    testside = -1;
    for (i = 1; i < shape->numSides; i++){
      if (shape->adjmatrix[sidemax][order[i]] == 0){
	if (nextside == -1){
	  nextside = order[i];
	}else{
	  testside = order[i];
	  break;
	}
      }
    }
  }

  //sidemax cannot be -1
  if (nextside == -1){
    fc_printfErrorMessage("no non-adj side - cant get thin shape order");
    free(areas);
    free(order);
    return FC_SUCCESS; // not a thin shape, null return
  }
  

  if (testside != -1){
    if (FC_FLT_EQUIV(areas[nextside],areas[testside])){
      fc_printfErrorMessage("two largest non-adj sides are of equal area - cant get thin shape order");
      free(areas);
      free(order);
      return FC_SUCCESS; // not a thin shape, null return
    }
  }

  //now we have our sides fill in all other ones randomly
  free(order);
  free(areas);

  order = (int*)malloc((shape->numSides)*sizeof(int));
  order[0] = sidemax;
  order[1] = nextside;
  count = 2;
  for (i = 0; i < shape->numSides; i++){
    if (i != sidemax && i!=nextside){
      order[count]= i;
      count++;
    }
  }
  
  *neworder = order;
  
  return FC_SUCCESS;
}



/** 
 * \ingroup  Shape
 * \brief As an approximation to obtaining the major
 *        sides of a thin shape, given a shape returns
 *        an array of the side ids in the ordered desired
 *        for a shape with a large and opposing side.
 *
 * \description
 * 
 *        As an approximation to obtaining the major
 *        sides of a thin shape, given a shape returns
 *        an array of the side ids in the ordered desired
 *        for a shape with a large and opposing side.
 *        A thin shape is a shape that in some
 *        dimension is narrow and in a roughly perpendicular
 *        direction has a pair of large opposing sides,
 *        which are the major
 *        sides that the user is interested in. The thin shape
 *        order would then have those major sides in the 0th and 1st
 *        position, with all other sides in random order.
 *        This can then be fed into \ref fc_reorderShape.
 *
 *        This method uses simplifying assumtions to get the
 *        major sides, which may not be correct in all cases.
 *        This method will iterate through the sides, starting
 *        with the largest one and attempt to find the
 *        the largest non-adjacent side for which the angle
 *        between the mean of the normals of each face for each side
 *        is in the range (180- the input angle) to 180.
 *        Currently, there is no restriction on the size of the
 *        sides relative to one another.
 *
 *        Note that there can be undesirable results for
 *        some shapes - see \ref fc_getSideNormalMeanSdev
 *        for comments onthe limitations of the mean of
 *        the normals and angle calculation. See
 *        also \ref fc_createLargeAndNonAdjacentSidesOrder for
 *        another approximation to the ideal thin shape order.
 *
 *        Note that there are no chekcs for unique sides like
 *        there are in the largeandnonadj side method. (that
 *        is, when there are multiple side options with the
 *        same area)
 *
 *        Returns with an error if there is something
 *        fundamentally wrong. Returns FC_SUCCESS but empty
 *        order array if the shape does not satisfy the
 *        large opposing sides criteria.
 *
 * \todo
 *   - there are no checks on uniqueness in here like
 *     there are in non-adj side, prob want to put those in.
 *   - make sure both sides are large? add an arg that also
 *     limits the ration of the sides relative to each other ?
 *
 * \modifications 
 *   - 08/04/06 ACG created
 */
FC_ReturnCode fc_createLargeAndOpposingSidesOrder(
  FC_Shape *shape, /**< input - shape */
  double maxangle, /**< input - opposing angle range*/
  int **neworder       /**< output - ids of sides */
  ) {

  FC_ReturnCode rc;
  FC_Vector* mean;
  int numSides;
  int *order;
  int currside, compside;
  double angle;
  int found = 0;
  int i,j;

  //defaults
  if (neworder)
    *neworder = NULL;

  if (!shape ||  !neworder || maxangle < 0 || maxangle > 180){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // log message
  fc_printfLogMessage("Getting large opposing sides order");

  numSides = shape->numSides;

  //must be at least 3 sides for this to do anything
  if (numSides < 3){
    fc_printfErrorMessage("shape must have at least 3 sides to have opposing sides/thin shape");
    //    return FC_INPUT_ERROR;
    return FC_SUCCESS;//null return
  }

  rc = fc_createDescendingAreaOrder(shape,&order);
  if (rc != FC_SUCCESS){
    fc_printfErrorMessage("cant get ascending area order");
    return rc;
  }

  mean = (FC_Vector*)malloc((shape->numSides)*sizeof(FC_Vector));
  if (!mean){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }

  for (i = 0; i < shape->numSides; i++){
    FC_Vector junk;
    rc = fc_getSideNormalMeanSdev(shape,i,&mean[i],&junk);
    if (rc != FC_SUCCESS){
      free(mean);
      return rc;
    }
  }


  for (i = 0; i < numSides; i++){
    currside = order[i];
    for (j = i+1; j < numSides; j++){
      compside = order[j];
      if (!shape->adjmatrix[currside][compside]){
	rc = fc_calcAngleBetweenVectors(mean[currside],mean[compside],
					&angle);
	if (rc != FC_SUCCESS){
	  fc_printfErrorMessage("Can't get angle between vectors");
	  free(mean);
	  free(order);
	  return rc;
	}
	if (angle >= (180-maxangle)){
	  found = 1;
	  break;
	}
      }
    }
    if (found){
      break;
    }
  }

  free(mean);

  if (found){
    int count = 2;
    order[0] = currside;
    order[1] = compside;
    //now we have our sides fill in all other ones randomly
    for (i = 0; i < shape->numSides; i++){
      if (i != currside && i != compside){
	order[count] = i;
	count++;
      }
    }
   
    *neworder = order;
  }else{
    //not a large opposing sides shape, null return
    free(order);
  }

  return FC_SUCCESS;
}




//@}


/** \name Adjacencies */
//-------------------------------
//@{


/**
 * \ingroup Shape
 * \brief Given an array of subsets specifying the faces of each
 *      side of a shape, returns a matrix specifying the side adjacencies
 *
 * \description
 *
 *      Given an array n of subsets, specifying the side faces 
 *      returns an nxn matrix of 1s and 0s depending on whether or
 *      not the sides are adjacenct. Have opted to make the diagonal elements 0.
 *      Note that the shapecreating methods automatically create adjacency
 *      matricies in a somewhat different way.
 *
 *      If any of the subsets are empty or NULL, this method returns with an error.
 *
 * \modifications
 *    - 02/07/06 ACG created.to handle the screw case.
 *    - 02/08/06 ACG unittest added.
 *    - 02/08/06 ACG found another flaw in the logic of old adjacencyMatrix test.
 *      Therefore removed old version entirely. 
 */
FC_ReturnCode fc_createSideAdjacencyMatrix(
  int numSubsets,       /**< input - number of input subsets */ 
  FC_Subset *subsets,   /**< input - array of face subsets */
  int shared_dim,       /**< input - shared dim for determining neighbors */
  int ***intmatrix /**< output -  1/0 matrix */
 ){
  FC_ReturnCode rc;
  FC_Mesh mesh;
  FC_AssociationType assoc;
  int i,j,k;
  
  //default return
  if (intmatrix)
    *intmatrix = NULL;

  if (!numSubsets || !subsets || !intmatrix || !(shared_dim == 1 || shared_dim == 0)){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  //check that all the subsets are valid
  for (i = 0; i < numSubsets; i++){
    if (!fc_isSubsetValid(subsets[i])){
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
      return FC_INPUT_ERROR;
    }
    rc = fc_getSubsetAssociationType(subsets[i], &assoc);
    if (rc != FC_SUCCESS) 
      return rc;

    if (assoc != FC_AT_FACE){
      fc_printfErrorMessage("Wrong assoc type");
      return FC_INPUT_ERROR;
    }

    rc= fc_getSubsetNumMember(subsets[i], &k);
    if (rc != FC_SUCCESS) 
      return rc;

    if (k < 1){
      fc_printfErrorMessage("sides must have at least one member");
      return FC_INPUT_ERROR;
    }
    
  }

  // log message
  fc_printfLogMessage("Computing Side Adjacency Matrix");

  rc = fc_getMeshFromSubset(subsets[0],&mesh);
  if (rc !=FC_SUCCESS){
    fc_printfErrorMessage("Can't get mesh from subset");
    return rc;
  }


  // create space for the adj array, initialize to zero
  *intmatrix = calloc(numSubsets, sizeof(int*));
  if (*intmatrix == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  for (i = 0; i < numSubsets; i++) {
    (*intmatrix)[i] = calloc(numSubsets, sizeof(int));
    if ((*intmatrix)[i] == NULL) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      return FC_MEMORY_ERROR;
    }
  }

  for (i = 0; i < numSubsets; i++){
    int numNbr, *nbrIDs;
    FC_Subset neighborfaces;

    //    printf("determining adj for subset %d\n",i);

    (*intmatrix)[i][i] = 0;


    //get the face neighbors for this subset
    //    printf("getting subset neighbors\n");
    rc = fc_getSubsetNeighbors(subsets[i], 1, shared_dim, &assoc, &numNbr, &nbrIDs);
    if (rc != FC_SUCCESS || assoc != FC_AT_FACE){
      fc_printfErrorMessage("Failed to get subset neighbors");
      for (k = 0; k < numSubsets; k++){
	free((*intmatrix)[k]);
      }
      free((*intmatrix));
      (*intmatrix) = NULL;
      return rc;
    }

    rc = fc_createSubset(mesh,"temp",FC_AT_FACE,&neighborfaces);
    if (rc !=FC_SUCCESS){
      fc_printfErrorMessage("cant create subset");
      for (k = 0; k < numSubsets; k++){
	free((*intmatrix)[k]);
      }
      free((*intmatrix));
      (*intmatrix) = NULL;
      free(nbrIDs);
      return rc;
    }
    
    //    printf("Adding array members to subset\n");
    rc = fc_addArrayMembersToSubset(neighborfaces,numNbr,nbrIDs);
    if (rc !=FC_SUCCESS){
      fc_printfErrorMessage("cant create subset");
      for (k = 0; k < numSubsets; k++){
	free((*intmatrix)[k]);
      }
      free((*intmatrix));
      (*intmatrix) = NULL;
      fc_deleteSubset(neighborfaces);
      free(nbrIDs);
      return rc;
    }
    free(nbrIDs);
    //    printf("done adding array members to subset\n");

    for (j = i+1; j < numSubsets; j++){ //symmetric result
      int retval;
      //      printf("checking intersection\n");
      rc = fc_doSubsetsIntersect(neighborfaces,subsets[j], &retval);
      if (rc != FC_SUCCESS){
        fc_printfErrorMessage("Failed to get subset intersection");
	for (k = 0; k < numSubsets; k++){
	  free((*intmatrix)[k]);
	}
	free((*intmatrix));
	(*intmatrix) = NULL;
	fc_deleteSubset(neighborfaces);
	return rc;
      }
      (*intmatrix)[i][j] = retval;
      (*intmatrix)[j][i] = retval;
    }
    //    printf("done checking intersection\n");

    fc_deleteSubset(neighborfaces);
  }

  return FC_SUCCESS;
}

/**
 * \ingroup Shape
 * \brief Given an adjacency matrix, returns segments of it based on
 *        connectivities 
 *
 *
 * \description
 *  
 *        Given an adjacency matrix, returns segments of it based on
 *        connectivities 
 *
 *        Note that pieces that are touching will come out as
 *        1 distinct shape with a lot of edges as opposed to
 *        two distinct shapes. Segment will do the same thing.
 *
 *
 * \modifications
 *    - 06/01/06 ACG created.
 *    - 06/23/06 ACG unittested, added defaults and input checking
 **/
FC_ReturnCode fc_segmentAdjacencyMatrix(
 int matrixdim, /**< input - dim of the adj matrix */
 int **adjmatrix, /**< input -  adj matrix */
 int *numsegments, /**< output - number of distinct shapes */
 int **sidespersegment, /**< output - number of sides for each distinct shape */
 int ***sideids /**<output - array of sideids(orig matrix numbers) for each distinct shape */
){
  int numSides;

  int **partsideids; /* array of arrays, 1 for each part */
  int *sidesperpart;
  int numparts;
  int maxparts;

  int *search;
  int i,j,k;

  //defaults
  if (numsegments)
    *numsegments = -1;
  if (sidespersegment)
    *sidespersegment = NULL;
  if (sideids)
    *sideids = NULL;


  if (matrixdim < 1 || !adjmatrix || !numsegments || !sidespersegment ||
      !sideids){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // log message
  fc_printfLogMessage("Segmenting Adj Matrix");

  //test that the matrix is symm and only has 0 & 1
  for (i = 0; i < matrixdim; i++){
    if (adjmatrix[i][i] != 0){
      fc_printfErrorMessage("adj matrix diag elem must be zero");
      return FC_ERROR;
    }
    for (j = 0; j < matrixdim; j++){
      if (adjmatrix[i][j] != adjmatrix[j][i] || 
	  (adjmatrix[i][j] != 1 && adjmatrix[i][j] != 0)){
	fc_printfErrorMessage("adj matrix must be symm and elems must be 0 or 1");
	return FC_ERROR;
      }
    }
  }

  numSides = matrixdim;

  partsideids = (int**)malloc(sizeof(int*));
  sidesperpart = (int*)malloc(sizeof(int));
  search = (int*)malloc(numSides*sizeof(int));  
  if (partsideids == NULL || sidesperpart == NULL || search == NULL){
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      return FC_MEMORY_ERROR;
  }
  for(i = 0; i < numSides; i++){
    search[i] = 1;
  }
  numparts = 0;
  maxparts = 1;

  for (i = 0; i < numSides; i++){
    if (search[i]){//start a new part
      int* currpartids;
      int currmax, currnum;

      currpartids = (int*)malloc(sizeof(int));
      if(currpartids == NULL){ 
	fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
	return FC_MEMORY_ERROR;
      }
      currpartids[0] = i;
      currnum = 1;
      currmax = 1;
      search[i] = 0;


      for (j = 0; j < currnum; j++){
	//find all the adj of this member
	for (k = 0; k < numSides; k++){
	  if (search[k] && adjmatrix[currpartids[j]][k] == 1){
	    //add to the currpart ids
	    if (currnum == currmax){
	      int *temp;
	      temp = (int*)realloc(currpartids,(2*currmax)*sizeof(int));
	      if (temp == NULL){
		fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
		return FC_MEMORY_ERROR;
	      }
	      currmax*=2;
	      currpartids = temp;
	    }
	    currpartids[currnum] = k; //add this side to the current part
	    currnum++;
	    search[k] = 0; //take it out of the options
	  }
	}
      }

      //now hand it off
      if (numparts == maxparts){
	int **temp;
	int *temp2;
	temp = (int**)realloc(partsideids,(2*maxparts)*sizeof(int*));
	if (temp == NULL){
	  fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
	  return FC_MEMORY_ERROR;
	}
	temp2 = (int*)realloc(sidesperpart,(2*maxparts)*sizeof(int));
	if (temp == NULL){
	  fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
	  return FC_MEMORY_ERROR;
	}
	maxparts*=2;
	partsideids = temp;
	sidesperpart = temp2;
      }
      partsideids[numparts] = currpartids;
      sidesperpart[numparts] = currnum;
      numparts++;
    }
  }

  *numsegments = numparts;
  *sidespersegment =sidesperpart;
  *sideids = partsideids;

  free(search);

  return FC_SUCCESS;
}


/**
 * \ingroup Shape
 * \brief  Given an nxn adjacency matrix returns the (shortestpath) 
 *   distance matrix
 *
 * \description
 *  
 *  Given an nxn adjacency matrix returns the (shortest path)
 *  distance matrix. This is not limited to cases where the adj matrix represents
 *  the adj of the sides, but it is intended to be used in that way.
 *  Does not take advantage of the symmetry of the adj matrix,
 *  so this could be used for non symm matricies. Uses Floyd's algorithm
 *
 *  Note: -1 is used for infinite distance entries. 
 *
 *  Note: in theory, this should work for weighted entries; I haven't
 *  tried this however, since this method isnt in here for general
 *  matrix calculations, but as a precursor for determining
 *  side pairs of interest.
 *
 * \modifications
 *    - 12/21/05 ACG created.
 *    - 01/30/05 ACG changed return val so infinite distance entries
 *      return value = -1.
 */
FC_ReturnCode fc_calcDistanceMatrix(
  int nsides, /**< input - nsides (nxn matrix) */
  int **adjmatrix, /**< input - adj matrix */
  int ***distmatrix /**< output - distance matrix */
){
  int** imatrix_a; //intermediatematrix
  int** imatrix_b; //intermediatematrix
  int ***curr, ***next, ***temp;
  int i,j,k;


  //default return
  if (distmatrix)
    *distmatrix = NULL;

  if (!nsides || !adjmatrix || !distmatrix){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // log message
  fc_printfLogMessage("Computing Side Distance Matrix");

  // create space for the adj array, initialize to zero
  imatrix_a = calloc(nsides, sizeof(int*));
  imatrix_b = calloc(nsides, sizeof(int*));
  if (imatrix_a == NULL || imatrix_b == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  for (i = 0; i < nsides; i++) {
    imatrix_a[i] = calloc(nsides, sizeof(int));
    imatrix_b[i] = calloc(nsides, sizeof(int));
    if (imatrix_a[i] == NULL || imatrix_b[i] == NULL) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      return FC_MEMORY_ERROR;
    }
  }

  //intial vals
  for (i = 0; i < nsides; i++){
    imatrix_a[i][i] = 0;
    for (j = 0; j < nsides; j++){ 
      if (i == j){
	imatrix_a[i][i] = 0;
      } else {
	if (adjmatrix[i][j] == 1){
	  imatrix_a[i][j] = 1;
	} else {
	  imatrix_a[i][j] = -1; 
	  //	  imatrix_a[i][j] = INT_MAX; 
	}
      }
    }
  }

  next = &imatrix_a; 
  curr = &imatrix_b; 
  for (k = 0; k < nsides; k++){
    temp = next;
    next = curr;
    curr = temp;
    for (i = 0; i < nsides; i++){
      for (j = 0; j < nsides; j++){
	int a = (*curr)[i][j];
	int b;
	//	if ((*curr)[i][k] == INT_MAX || (*curr)[k][j] == INT_MAX){
	//	  b = INT_MAX;
	if ((*curr)[i][k] == -1 || (*curr)[k][j] == -1){
	  b = -1;
	} else {
	  b = (*curr)[i][k]+(*curr)[k][j];
	}

	if (a == -1){
	  (*next)[i][j] = b;
	}else if (b == -1){
	  (*next)[i][j] = a;
	}else {
	  (*next)[i][j] = a < b ? a : b;
	}
	//	  (*next)[i][j] = a < b ? a : b;
      }
    }
  }

  for (i = 0; i < nsides; i++){
    free((*curr)[i]);
  }
  free(*curr);

  *distmatrix = *next;
  return FC_SUCCESS;
}


//@}




/** \name Misc */
//-------------------------------
//@{


/** 
 * \ingroup  Shape
 * \brief given a shape returns the shape's ends.
 *
 * \description
 * 
 *    Given a shape an array containing the ids of the shape's ends. 
 *    The ends are those sides who are connected
 *    to only 1 other side in the shape (Determined by consideration of the
 *    adjacency matrix). The ids are those of the indicies of the sides
 *    in the order they are specified in the face array. Order of the
 *    indicies in the return array are in order from smallest 
 *    (area) to largest.
 *
 *    Note: nothing checks that the shape actually makes any sense.
 *     
 *    Returns -1 for numEnds and NULL array in case of error.
 *    If numSides == 0, it is considered an error.
 *    If numSides == 1 it is not an error and returns zero ends. 
 *
 * \todo
 *   - decide if this should be expanded to ask for any adj matrix value
 *
 * \modifications 
 *   - 02/07/06 ACG, created
 *   - 06/18/06 ACG, created - revision of old getSidesEnds
 *   - 06/28/06 ACG  now the official shape ends
 */
FC_ReturnCode fc_getShapeEnds(
  FC_Shape *shape, /**< input - shape */
  int *numEnds,        /**< output - num of ends */
  int **endArray       /**< output - ids of ends */
) {

  int* eholding;
  int nends;
  int numSides;
  int i,j;

  //defaults
  if (numEnds)
    *numEnds = -1;
  if (endArray)
    *endArray = NULL;

  if (!shape || !numEnds || !endArray){ 
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // log message
  fc_printfLogMessage("Getting shape ends B");

  numSides = shape->numSides;
  if (!numSides){
    fc_printfErrorMessage("shape has no sides");
    return FC_ERROR;
  }

  eholding = calloc(numSides,sizeof(int));
  if (!eholding) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    //not cleaning up matrix in case of memory error
    return FC_MEMORY_ERROR;
  }

  nends = 0;
  for (i = 0; i < numSides; i++){
    int conn = 0;

    if (shape->adjmatrix[i][i] != 0){
      fc_printfErrorMessage("non zero diagonal adj matrix");
      free(eholding);
      return FC_ERROR;
    }
    for (j = 0; j < numSides; j++){
      if (shape->adjmatrix[i][j] != shape->adjmatrix[j][i]){
	fc_printfErrorMessage("non symmetric adj matrix");
	free(eholding);
	return FC_ERROR;
      }
      if (shape->adjmatrix[i][j] !=0 ){
	conn+=1;
      }
    }
    if (conn == 1){
      eholding[nends++]=i;
    }
  }

  switch (nends){
  case 0:
    free(eholding);
    eholding = NULL;
    break;
  case 1:
    //do nothing
    break;
  default: 
    {
      FC_ReturnCode rc;
      double *areas;

      areas = malloc(nends*sizeof(FC_Subset));
      if (!areas){
	fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
	//not cleaning up in case of memory error
	return FC_MEMORY_ERROR;
      }
    
      for (i = 0; i < nends; i++){
	rc= fc_getSubsetArea(shape->faces[eholding[i]],&areas[i]);
	if (rc != FC_SUCCESS){
	  fc_printfErrorMessage("cant get subset areas");
	  free(eholding);
	  free(areas);
	  return rc;
	}
      }
    
    
      //order sides by area
      //sort sides areas. yes i know, its bubblesort.
      {
	double tempa;
	int tempi;
	int test; 
	
	for(i = nends - 1; i > 0; i--){
	  test=0;
	  for(j = 0; j < i; j++){
	    if (areas[j] > areas[j+1]){
	      tempa = areas[j]; 
	      areas[j] = areas[j+1];
	      areas[j+1] = tempa;
	      tempi = eholding[j]; 
	      eholding[j] = eholding[j+1];
	      eholding[j+1] = tempi;
	      test=1;
	    }
	  }
	  if(test==0) break; //exit if sorted
	}
      }
      free(areas);
      break;
    }
  }

  *numEnds = nends;
  //note: eholding was allocated with numSides members, not nends.
  //this is ok becuase free will free the whole array
  *endArray = eholding;

  return FC_SUCCESS;
}


/**
 * \ingroup Shape
 * \brief  Given a Shape, returns the (non-displaced) area 
 *  for each side
 *
 * \description
 *  
 *  Given a Shape, returns the (non-displaced) area 
 *  for each side. Checks to see if all the sides are
 *  on the same mesh.
 *
 *  If there is an error in calculating this for any side, then
 *  this returns NULL for the entire array.
 *
 * \todo
 *   - may want to do a displaced version
 *
 * \modifications
 *   - 01/03/06 ACG
 *   - 07/17/06 ACG born out of old getSidesAreas
 */
FC_ReturnCode fc_getShapeSidesAreas(
  FC_Shape *shape, /**< input - ptr to shape */
  double **area /**< output - area array */
){
  FC_Mesh meshlast, meshnext;
  FC_ReturnCode rc;
  FC_Variable faceareas;
  double *areaarray;
  FC_AssociationType assoc;
  int i,j;

  //default return
  if (area)
    *area = NULL;

  if (!shape || !area ){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  //check that all the subsets are valid
  for (i = 0; i < shape->numSides; i++){
    if (!fc_isSubsetValid(shape->faces[i])){
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
      return FC_INPUT_ERROR;
    }
    rc = fc_getSubsetAssociationType(shape->faces[i], &assoc);
    if (rc != FC_SUCCESS)
      return rc;
    if ( assoc != FC_AT_FACE){
      fc_printfErrorMessage("Wrong association type");
      return FC_INPUT_ERROR;
    }
  }

  // log message
  fc_printfLogMessage("Computing ShapeSide Areas");

  //get areas for all faces of the mesh
  rc = fc_getMeshFromSubset(shape->faces[0],&meshlast);
  if (rc!= FC_SUCCESS){
    return rc;
  }

  rc = fc_getFaceAreas(meshlast,&faceareas);
  if (rc!= FC_SUCCESS){
    fc_printfErrorMessage("cant get face areas");
    return rc;
  }

  rc = fc_getVariableDataPtr(faceareas,(void **) &areaarray);
  if (rc!= FC_SUCCESS){
    return rc;
  }

  *area = calloc(shape->numSides, sizeof(double));
  if (*area == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }

  for (i = 0; i < shape->numSides; i++){
    int numMem, *ids;

    rc = fc_getMeshFromSubset(shape->faces[i],&meshnext);
    if (rc!= FC_SUCCESS){
      free(*area);
      *area = NULL;
      return rc;
    }

    if (!FC_HANDLE_EQUIV(meshnext, meshlast)){
      fc_printfErrorMessage("all subsets must be on the same mesh");
      free(*area);
      *area = NULL;
      return FC_INPUT_ERROR;
    }

    rc = fc_getSubsetMembersAsArray(shape->faces[i], &numMem, &ids);
    if (rc!= FC_SUCCESS){
      free(*area);
      *area = NULL;
      return rc;
    }

    for (j = 0; j < numMem; j++){
      (*area)[i] += areaarray[ids[j]];
    }
    free(ids);

    meshlast = meshnext;
  }


  fc_deleteVariable(faceareas);

  return FC_SUCCESS;
}


/**
 * \ingroup Shape
 * \brief Given a shape and the ID of a planar side, returns
 *        the outward normal to that face. 
 *
 * \description
 *  
 *        Given a shape and the ID of a planar side, returns
 *        the outward normal to that face. Does not actually test to
 *        see if the face is planar. 
 *
 * \todo
 *    - decide if need more unittest
 *
 * \modifications
 *    - 03/16/06 ACG created.
 *    - 08/08/06 ACG changed interface so pass sides in terms
 *                   of shape rahter than faces and elems.
 */
FC_ReturnCode fc_getPlanarSideNormal(
  FC_Shape *shape, /**< input - the shape */				     
  int sideid, /**< input - the shape side id */
  FC_Coords *normal /**< output - normal */
  ){
  
  FC_ReturnCode rc1, rc2, rc;
  FC_Mesh mesh1,mesh2;
  FC_AssociationType assoc1, assoc2;
  FC_Vector retvec;
  int faceID,elemID;
  int numItems;
  int *arr;
  int i,j;

  //default
  if (normal){
    for (i = 0; i < 3; i++){
      (*normal)[i] = 0;
    }
  }
    
  if (!shape || (sideid < 0) || (sideid > (shape->numSides)-1) ||
      !fc_isSubsetValid(shape->faces[sideid]) || 
      !fc_isSubsetValid(shape->elems[sideid]) || !normal ){
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
      return FC_INPUT_ERROR;
  }

  rc1 = fc_getSubsetAssociationType(shape->faces[sideid], &assoc1);
  rc2 = fc_getSubsetAssociationType(shape->elems[sideid], &assoc2);
  if (rc1 != FC_SUCCESS || rc2 != FC_SUCCESS) 
    return(rc1 != FC_SUCCESS ? rc1 : rc2);

  if (assoc1 != FC_AT_FACE || assoc2 != FC_AT_ELEMENT){
    fc_printfErrorMessage("Wrong assoc type");
    return FC_INPUT_ERROR;
  }

  rc1 = fc_getMeshFromSubset(shape->faces[sideid],&mesh1);
  rc2 = fc_getMeshFromSubset(shape->elems[sideid],&mesh2);
  if (rc1 != FC_SUCCESS || rc2 != FC_SUCCESS) 
    return(rc1 != FC_SUCCESS ? rc1 : rc2);

  if (!FC_HANDLE_EQUIV(mesh1,mesh2)){
    fc_printfErrorMessage("faces and elems must be on same mesh");
    return FC_INPUT_ERROR;
  }

  
  //since its planar just pick the first face
  rc = fc_getSubsetMembersAsArray(shape->faces[sideid],&numItems,&arr);
  if (rc != FC_SUCCESS){
    fc_printfErrorMessage("Can't get faces");
    return rc;
  }
  if (numItems < 1){
    fc_printfErrorMessage("must have at least 1 face");
    return FC_INPUT_ERROR;
  }
  faceID = arr[0];
  free(arr);

  //find the elem that gave rise to it - can be only one.
  rc = fc_getSubsetMembersAsArray(shape->elems[sideid],&numItems,&arr);
  if (rc != FC_SUCCESS){
    fc_printfErrorMessage("Can't get elems");
    return rc;
  }

  elemID = -1;
  for (i = 0; i < numItems; i++){
    int numChild, *childIDs;
    rc = fc_getMeshEntityChildren(mesh1, FC_AT_ELEMENT, arr[i], FC_AT_FACE,
                             &numChild, &childIDs);
    if (rc != FC_SUCCESS){
      fc_printfErrorMessage("Cant get elem->face data");
      free(arr);
      return rc;
    }
    for (j = 0; j < numChild; j++){
      if (childIDs[j] == faceID){
	elemID = arr[i];
	break;
      }
    }
    free(childIDs);

    if (elemID != -1){
      break;
    }
  }
  free(arr);
   
  if (elemID == -1){
    fc_printfErrorMessage("Cant find elem that gave rise to this face");
    return FC_INPUT_ERROR;
  }

  rc =  _fc_getOutwardFaceNormal(mesh1,faceID,elemID, &retvec);
  if (rc != FC_SUCCESS){
    fc_printfErrorMessage("Cant get outward face normal");
    return FC_INPUT_ERROR;
  }

  for (i = 0; i < 3; i++){
    (*normal)[i] = retvec[i];
  }

  return FC_SUCCESS;
}

/**
 * \ingroup Shape
 * \brief Given a shape and a sideid, returns
 *        the outward normal and std in each direction to that face. 
 *
 * \description
 *  
 *        Given a shape and a sideid, returns
 *        the mean of the components of the outward normals
 *        and std in each direction to that face. Note that the
 *        mean is not a normal and also that for the same shape,
 *        this can retrun differnet values depending onthe sizes
 *        and shapes of the individual faces of the overall side.
 *        The hope is that in most cases, if one is using this to compare
 *        the results of different sides inthe same mesh, that the
 *        elements are sufficiently regular that a meaningful result
 *        is obtained.
 *
 * \modifications
 *    - 08/04/06 ACG created.
 *    - 08/08/06 ACG changed interface so takes a shape rather
 *                   than the faces and elems of the side
 */
FC_ReturnCode fc_getSideNormalMeanSdev(
  FC_Shape *shape, /**< input - shape */
  int sideid, /**< input - side id */
  FC_Vector *mean, /**< output - mean */
  FC_Vector *std /**< output - std in each dir separately */
  ){
  
  FC_ReturnCode rc1, rc2, rc;
  FC_Mesh mesh1, mesh2;
  FC_AssociationType assoc1, assoc2;
  FC_Vector sum, sumsq;
  FC_Vector compmean, sdev, temp;
  int faceID,elemID;
  int numFaces,numElems;
  int *arrf, *arre;
  FC_SortedIntArray sia;
  int i,j,k;

  //default
  if (mean){
    for (i = 0; i < 3; i++){
      (*mean)[i] = 0;
    }
  }

  if (std){
    for (i = 0; i < 3; i++){
      (*std)[i] = 0;
    }
  }

    
  if (!shape || (sideid < 0) || (sideid > ((shape->numSides)-1)) ||
      !fc_isSubsetValid(shape->faces[sideid]) || 
      !fc_isSubsetValid(shape->elems[sideid]) || !mean  || !std){
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
      return FC_INPUT_ERROR;
  }

  rc1 = fc_getSubsetAssociationType(shape->faces[sideid], &assoc1);
  rc2 = fc_getSubsetAssociationType(shape->elems[sideid], &assoc2);
  if (rc1 != FC_SUCCESS || rc2 != FC_SUCCESS) 
    return(rc1 != FC_SUCCESS ? rc1 : rc2);

  if (assoc1 != FC_AT_FACE || assoc2 != FC_AT_ELEMENT){
    fc_printfErrorMessage("Wrong assoc type");
    return FC_INPUT_ERROR;
  }

  rc1 = fc_getMeshFromSubset(shape->faces[sideid],&mesh1);
  rc2 = fc_getMeshFromSubset(shape->elems[sideid],&mesh2);
  if (rc1 != FC_SUCCESS || rc2 != FC_SUCCESS) 
    return(rc1 != FC_SUCCESS ? rc1 : rc2);

  if (!FC_HANDLE_EQUIV(mesh1,mesh2)){
    fc_printfErrorMessage("faces and elems must be on same mesh");
    return FC_INPUT_ERROR;
  }

  
  rc = fc_getSubsetMembersAsArray(shape->faces[sideid],&numFaces,&arrf);
  if (rc != FC_SUCCESS){
    fc_printfErrorMessage("Can't get faces");
    return rc;
  }
  if (numFaces < 1){
    fc_printfErrorMessage("must have at least 1 face");
    free(arrf);
    return FC_INPUT_ERROR;
  }
  rc = fc_getSubsetMembersAsArray(shape->elems[sideid],&numElems,&arre);
  if (rc != FC_SUCCESS){
    fc_printfErrorMessage("Can't get elems");
    free(arrf);
    free(arre);
    return rc;
  }
  
  rc = fc_initSortedIntArray(&sia);
  if (rc != FC_SUCCESS){
    fc_printfErrorMessage("Can't init sorted int array");
    free(arrf);
    free(arre);
    return rc;
  }

  //when we get the array from the subset, we know it is sorted
  i = fc_addIntArrayToSortedIntArray(&sia,numFaces,arrf,1);
  free(arrf);

  if (!i){
    fc_printfErrorMessage("Can't add array to sia int array");
    free(arre);
    fc_freeSortedIntArray(&sia);
    return FC_ERROR;
  }
  

  for (i = 0; i < 3; i++){
    sum[i] = 0;
    sumsq[i] = 0;
  }

  //for each elem, get its kids and match them up with the faces
  //each face can come from only 1 elem.
  //figure its faster to get the elems kids once and check thru a sorted
  //face list than to try to go thru the face list and iterate thur all
  //the elems
  //then get the mean of the normals
  for (i = 0; i < numElems; i++){
    int numChild, *childIDs;
    elemID = arre[i];
    rc = fc_getMeshEntityChildren(mesh1, FC_AT_ELEMENT, arre[i], FC_AT_FACE,
				  &numChild, &childIDs);
    if (rc != FC_SUCCESS){
      fc_printfErrorMessage("Cant get elem->face data");
      free(arre);
      fc_freeSortedIntArray(&sia);
      return rc;
    }

    for (j = 0; j < numChild; j++){
      FC_Vector retvec;
      //is this face in the list ?
      faceID = childIDs[j];
      if (fc_deleteIntFromSortedIntArray(&sia,faceID)){
	rc =  _fc_getOutwardFaceNormal(mesh1,faceID,elemID, &retvec);
	if (rc != FC_SUCCESS){
	  fc_printfErrorMessage("Cant get outward face normal");
	  free(arre);
	  fc_freeSortedIntArray(&sia);
	  return rc;
	}

	for (k = 0; k < 3; k++){
	  sum[k]+=retvec[k];
	  sumsq[k]+=(retvec[k]*retvec[k]);
	}
      }
    }

    free(childIDs);
  }

  free(arre);
  //make sure all the faces have been removed
  if (sia.numVal){
    fc_printfErrorMessage("Mismatch of faces and elems");
    fc_freeSortedIntArray(&sia);    
    return FC_INPUT_ERROR;
  }
  fc_freeSortedIntArray(&sia);    

  for (i = 0; i <3; i++){
    compmean[i] = sum[i]/(double)numFaces;

    if (numFaces > 1){
      temp[i] = (sumsq[i] - (sum[i] * sum[i] /(double) numFaces)) /
	(double)(numFaces - 1);
      sdev[i] = sqrt(fabs(temp[i]));
      // might be negative if roundoff error
    }else{
      sdev[i] = 0.0;
    }
  }

  for (i = 0; i < 3; i++){
    (*mean)[i] = compmean[i];
    (*std)[i] = sdev[i];
  }

  return FC_SUCCESS;
}

//@}



