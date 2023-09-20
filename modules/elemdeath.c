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
 * \file elemdeath.c
 * \brief Implementation of \ref ElemDeath module.
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/modules/elemdeath.c,v $
 * $Revision: 1.64 $
 * $Date: 2006/10/19 03:14:50 $
 */

// C library includes
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <stdio.h>

// fc library dependencies
#include "base.h"
#include "storage.h"
#include "library.h"
#include "variable.h"
#include "subset.h"
#include "geom.h"
#include "topo.h"
#include "threshold.h"
#include "shape.h"
#include "sequence.h"
#include "feature.h"
#include "track.h"
#include "util.h"

// this module
#include "elemdeathP.h"
#include "elemdeath.h"

/**
 * \addtogroup ElemDeath
 * \brief Dead element analysis routines and helper functions
 *
 * \description
 *
 *   In many simulations, changing mesh topology is approximated by 
 *   allowing elements to "die". The mesh topology stays the same, but any
 *   elements that are labeld "dead" no longer participate in the simulation.
 *   Dead elements can be used to model rips, tears and other changes in the
 *   mesh.
 *
 *   The input of most of the routines in this module is a subset representing
 *   a dead element region, and is assumed to have the association type of
 *   FC_AT_ELEMENT. It is also assumed that the coordinates of vertices within
 *   an dead element region cannot be trusted (the vertices on the boundary of
 *   the dead element region may be o.k. if they are still on live elements).
 *   
 *   A dead element region does not have to be a single topological segment,
 *   but most of the results are more easily interpreted if this is true.
 */

/**
 * \ingroup ElemDeath
 * \defgroup PrivateElemDeath (Private)
 */


/** \name Decay of Subsets */
//----------------------------------------------------------
//@{


/**
 * \ingroup ElemDeath
 * \brief Get exposed skin.
 *
 * \description 
 *
 *   Given a subset of elements, this routine returns a subset of the entities
 *   that would become "exposed" (i.e. become part of the mesh skin) if the
 *   given elements were removed from the parent mesh. The exposed skin is the
 *   boundary between the dead region and the rest of the mesh. It is useful as
 *   a subset of the dead region because all members are still "real" unlike
 *   entities that are totally within a dead region which may no longer have
 *   physical meaning.
 *
 *   Since calculating the skin of the parent mesh can be very computationally
 *   expension, there is an optional argument to pass in an already computed
 *   mesh skin.
 *
 *   If there are no exposed entities, this routine returns successfully but
 *   with FC_NULL_SUBSET (instead of a subset with no members).
 *
 *   This function will fail if there are no members in the input subset.
 *
 * \modifications
 *   - 6/27/2005 WSD. Created.
 */
FC_ReturnCode fc_getExposedSkin(
  FC_Subset subset,          /**< input - A subset of "dead" elements */
  FC_Subset* opt_meshSkin,   /**< input - (optional) pointer to mesh skin 
				(or NULL) */
  FC_Subset* exposed         /**< output - The exposed skin */
) {
  FC_ReturnCode rc;
  FC_AssociationType in_assoc;
  int numMember, maxNumMember;
  FC_Mesh mesh, temp_mesh;
  FC_Subset mesh_skin, subset_skin, union_subset;
  
  // default returns
  if (exposed)
    *exposed = FC_NULL_SUBSET;

  // Test input
  if (!fc_isSubsetValid(subset) || !exposed ||
      (opt_meshSkin != NULL && !fc_isSubsetValid(*opt_meshSkin))) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  rc = fc_getSubsetInfo(subset, &numMember, &maxNumMember, &in_assoc);
  if (rc != FC_SUCCESS)
    return rc;
  if (in_assoc != FC_AT_ELEMENT && in_assoc != FC_AT_WHOLE_MESH) {
    fc_printfErrorMessage("Subset must have element or whole association");
    return FC_INPUT_ERROR;
  }

  // A little setup and more testing
  rc = fc_getMeshFromSubset(subset, &mesh);
  if (rc != FC_SUCCESS)
    return rc;
  if (opt_meshSkin) {
    rc = fc_getMeshFromSubset(*opt_meshSkin, &temp_mesh);
    if (rc != FC_SUCCESS)
      return rc;
    if (!FC_HANDLE_EQUIV(mesh, temp_mesh)) {
      fc_printfErrorMessage("Subset and mesh sking must be on the same mesh");
      return FC_INPUT_ERROR;
    }
  }

  // Error case - No dead elements // FIX? is this really an error?
  if (numMember == 0) {
    fc_printfErrorMessage("Cannot handle empty subset");
    return FC_INPUT_ERROR;
  }
    
  // Special case - all elements are dead
  if (in_assoc == FC_AT_WHOLE_MESH ||  numMember == maxNumMember) {
    // return NULL handle
    return FC_SUCCESS;
  }

  // Setup - get the skin's of the mesh and the subset
  if (opt_meshSkin)
    mesh_skin = *opt_meshSkin;
  else {
    rc = fc_getMeshSkin(mesh, &mesh_skin);
    if (rc != FC_SUCCESS)
      return rc;
  }
  rc = fc_getSubsetSkin(subset, &subset_skin);
  if (rc != FC_SUCCESS)
    return rc;

  // Exposed enties are subset_skin - any skin shared with orig mesh
  // Union of mesh & subset skins will give us overlap
  rc = fc_createSubsetIntersection(mesh_skin, "AND", subset_skin, "temp",
				   &union_subset);
  if (!opt_meshSkin)
    fc_deleteSubset(mesh_skin);
  if (rc != FC_SUCCESS) {
    fc_deleteSubset(subset_skin);
    return rc;
  }
  // Xor of overlap and subset skin should give skin NOT overlapping
  rc = fc_createSubsetIntersection(subset_skin, "XOR", union_subset,
				   "temp", exposed);
  fc_deleteSubset(subset_skin);
  fc_deleteSubset(union_subset);

  // done
  return rc;
}

/**
 * \ingroup ElemDeath
 * \brief Get decayed skin
 *
 * \description 
 *
 *   Given a dead element region and the skin of a mesh or subset, 
 *   returns a subset which is the intersection of the dead
 *   element region and the skin, referred to as the "decayedSkin".
 *   The return subset is named "DecayedSkin". 
 *
 *   The dead region must have association FC_AT_ELEMENT, the
 *   skin must have association FC_AT_FACE, and the return subset
 *   has the assocation FC_AT_FACE. Both input subsets must be on the
 *   same mesh. If you want the elements of the decayed face, call
 *   \ref fc_createSubsetIntersection directly on the dead region and
 *   subset giving rise to the input skin.
 *
 *   This method does not check that the initial subsets actually
 *   describe a dead region and a skin, so this can in fact be used to
 *   calculate the intersection of any element subset and some other skin.
 *   It is being placed in elemdeath however, because 
 *   this function was developed for a specific dead element need and
 *   and using the deadelement terminology in this function name
 *   is beneficial.
 *
 *   Notes:
 *    - if the deadregion is NULL or the input skin is empty,
 *      returns with error (this is to be consistent with
 *      \ref fc_getExposedSkin)
 *    - if the input skin is NULL, it creates the mesh skin
 *    - if the dead region is empty, returns empty subset.
 *    - if there are no intersections, returns a subset with no members.
 *
 *
 * \modifications
 *   - 11/30/2005 ACG Created via moving here from surface demo code.
 *   - 12/05/2005 ACG Added the unittest (see todo above).
 *   - 03/01/2006 ACG added explicit type checking for hte input subsets
 *   - 03/02/2006 ACG changing name for DecayedMeshSkin to DecayedSkin since
 *                    the skin doesnt have to be a meshskin
 *   - 07/16/2006 ACG This was fc_getDecayedSkinSegments, but removing 
 *                    the segmenting
 */
FC_ReturnCode fc_getDecayedSkin(
  FC_Subset deadregion, /**< input - dead region elements*/
  FC_Subset* opt_meshSkin, /**< input - optional ptr to mesh skin or NULL */
  FC_Subset* decayedSkin /**< output - skin items ids comprising the dead skin (not segmented) */
  ){
  FC_ReturnCode rc;
  FC_Mesh mesh1, mesh2;
  FC_Subset deadregionskin,meshskin,createmeshskin,deadfaces;
  FC_AssociationType assocdead, assocskin;
  int n_dead, n_skin;
  int junk;


  // default returns
  if (decayedSkin)
    *decayedSkin = FC_NULL_SUBSET;

  // Test input. 
  if (!fc_isSubsetValid(deadregion) || !decayedSkin ||
      (opt_meshSkin != NULL && !fc_isSubsetValid(*opt_meshSkin))){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  rc = fc_getSubsetInfo(deadregion,&n_dead,&junk,&assocdead);
  if (rc != FC_SUCCESS){
    fc_printfErrorMessage("cant get info of dead region");    
    return rc;
  }

  if (assocdead != FC_AT_ELEMENT){
    fc_printfErrorMessage("deadregion must have element assoc");
    return FC_INPUT_ERROR;
  }
  rc = fc_getMeshFromSubset(deadregion,&mesh1);
  if (rc != FC_SUCCESS)
    return rc;
  if (opt_meshSkin) {
    rc = fc_getMeshFromSubset(*opt_meshSkin, &mesh2);
    if (rc != FC_SUCCESS)
      return rc;
    if (!FC_HANDLE_EQUIV(mesh1, mesh2)) {
      fc_printfErrorMessage("deadregion and skin must be on the same mesh");
      return FC_INPUT_ERROR;
    }
  }

  createmeshskin = FC_NULL_SUBSET;
  if (opt_meshSkin){
    meshskin = *opt_meshSkin;
  } else {
    rc = fc_getMeshSkin(mesh1, &createmeshskin);
    if (rc != FC_SUCCESS)
      return rc;
    meshskin = createmeshskin;
  }
  rc = fc_getSubsetInfo(meshskin,&n_skin,&junk,&assocskin);
  if (rc != FC_SUCCESS){
    fc_printfErrorMessage("cant get info of skin");    
    fc_deleteSubset(createmeshskin);
    return rc;
  }
  if (assocskin != FC_AT_FACE){
    fc_printfErrorMessage("skin must have face assoc");
    fc_deleteSubset(createmeshskin);
    return FC_INPUT_ERROR;
  }
  //have to check if the skin is empty before checking if the dead set is
  if (n_skin == 0){
    fc_printfErrorMessage("cant work on an empty skin");
    fc_deleteSubset(createmeshskin);
    return FC_INPUT_ERROR;
  }
  if (n_dead == 0){
    return FC_SUCCESS;
  }

  fc_printfLogMessage("get decayed skin ");

  //get the skin of the dead region
  rc = fc_getSubsetSkin(deadregion,&deadregionskin);
  if (rc != FC_SUCCESS){
    fc_deleteSubset(createmeshskin);
    fc_printfErrorMessage("cant get dead region skin");
    return rc;
  } 
  
  //intersect to get the dead faces
  rc = fc_createSubsetIntersection(deadregionskin,"AND",meshskin,
				   "DeadSkin",&deadfaces);
  fc_deleteSubset(deadregionskin);
  fc_deleteSubset(createmeshskin);
  if (rc != FC_SUCCESS){
    fc_printfErrorMessage("cant get dead mesh faces");
    return rc;
  } 

  *decayedSkin = deadfaces;

  return FC_SUCCESS;
};


/**
 * \ingroup ElemDeath
 * \brief Get identifier telling if deadSubset intersects, or entirely 
 *        erodes a subset
 *
 * \description
 * 
 *   Given a dead subset, returns an identifier
 *   telling if the deadSubset intersects, or entirely erodes the 
 *   comparisonSubset.
 *
 *   Identifier values:
 *     -  0 if the deadSubset doesnt intersect the compSubset at all.
 *     -  1 if the deadSubset intersects the compSubset.
 *     -  2 if the deadSubset erodes the whole compSubset.
 *     -  -1 in case of error
 *
 *   The deadSubset and the compSubset both have to have the same assocation.
 *   It must either be FC_AT_FACE or FC_AT_ELEMENT. They must also be on the
 *   same mesh. This function is currently a wrapper around
 *   \ref fc_doSubsetsIntersect and \ref fc_isSubsetSuperset and is subject to the
 *   constraints of those methods.
 *
 *   Fails if the compsubset is empty or null of if the deadsubset is null.
 *   The dead set can be empty and returns 0. (The inner calls would support
 *   permuations of empty and null, but the logic of the call does not).
 *
 *   In general, it is expect that this will be used to compare a dead
 *   region and the side of a mesh, but there is nothing that checks that
 *   so this can in fact be used more generally, and also misused at your
 *   peril. It is being placed in elemdeath however, because this function 
 *   was developed for a specific dead element need and using the deadelement
 *   terminology in this function name is beneficial. If you want to do this
 *   comparison for all the sides of a shape, use 
 *   \ref fc_getShapeSidesDecayType. 
 *
 * \todo
 *   - make more efficient by not having to call both the intersection and the
 *     superset code.
 *
 * \modifications
 *   - 02/14/2006 ACG fc_getSidesDecay created, which replaced
 *                    fc_getDecayedMeshSkinSideIntersections.
 *   - 07/17/2006 ACG new version born of fc_getSidesDecay.
 */
FC_ReturnCode fc_getSubsetDecayType(
   FC_Subset deadSubset, /**< input - deadSubset */
   FC_Subset compSubset, /**< input - comparison subset */
   int* sidedecayflag /**<output - decay flag  */
   ){

  FC_ReturnCode rc;
  FC_AssociationType assocdead, assoccomp;
  int ret = -1;
  int i;

  //default returns
  if (sidedecayflag)
    *sidedecayflag = -1;

  if (!fc_isSubsetValid(deadSubset) || !fc_isSubsetValid(compSubset) ||
      !sidedecayflag){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  if (FC_HANDLE_EQUIV(compSubset,FC_NULL_SUBSET)){
    fc_printfErrorMessage("comp subset cannot be NULL");
    return FC_INPUT_ERROR;
  } 

  rc = fc_getSubsetNumMember(compSubset,&i);
  if (rc != FC_SUCCESS){
    fc_printfErrorMessage("cant get num members of comp subset");
    return rc;
  } 
  if (!i){
    fc_printfErrorMessage("comp subset cannot be empty");
    return FC_INPUT_ERROR;
  }

  //check assoc, mesh is checked in internal calls. 
  rc = fc_getSubsetAssociationType(deadSubset,&assocdead);
  if (assocdead != FC_AT_ELEMENT && assocdead!= FC_AT_FACE){
    fc_printfErrorMessage("deadSubset must have face or element assoc");
    return FC_INPUT_ERROR;
  }
  rc = fc_getSubsetAssociationType(compSubset,&assoccomp);
  if (assoccomp != assocdead){
    fc_printfErrorMessage("compSubset must have smae assoc as deadSubset");
    return FC_INPUT_ERROR;
  }

  rc = fc_doSubsetsIntersect(deadSubset,compSubset,&ret);
  if (rc != FC_SUCCESS){
    return rc;
  }
  if (ret == 1){
    int ind2;
    
    rc = fc_isSubsetSuperset(deadSubset,compSubset,&ind2);
    if (rc != FC_SUCCESS){
      return rc;
    }
    ret = (ind2 == 1 ? 2:1);
  }

  *sidedecayflag  = ret;
  return FC_SUCCESS;
};


//@}


/** \name Decay of Shapes */
//-------------------------------
//@{

/**
 * \ingroup ElemDeath
 * \brief Get the decayed Skin
 *
 * \description 
 *
 *   Analogous method to \ref fc_getDecayedSkin but where
 *   the meshskin is that defined by a shape. If you want
 *   the decay information on a side-by-side basis, use
 *   \ref fc_getDecayedShapeSides.
 *   
 * \modifications
 *   - 07/16/2006 ACG created
 *                    
 */
FC_ReturnCode fc_getDecayedShapeSkin(
  FC_Subset deadregion, /**< input - dead region elements*/
  FC_Shape  *shape, /**< input - ptr to shape */
  FC_Subset* decayedSkin /**< output - decayedSkin */
  ){
  FC_ReturnCode rc;
  FC_Subset *sidesdecay;
  int i,j;

  // default returns
  if (decayedSkin)
    *decayedSkin = FC_NULL_SUBSET;

  // Test input. 
  if (!fc_isSubsetValid(deadregion) || !shape || !shape->numSides ||
      !decayedSkin){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  rc = fc_getDecayedShapeSides(deadregion,shape,&sidesdecay);
  if (rc != FC_SUCCESS){
    fc_printfErrorMessage("cant get decayed shape skin");
    return rc;
  }


  for (i = 1; i < shape->numSides; i++){
    int num, *arr;
    if (!FC_HANDLE_EQUIV(sidesdecay[i],FC_NULL_SUBSET)){
      rc = fc_getSubsetMembersAsArray(sidesdecay[i],&num,&arr);
      if (rc != FC_SUCCESS){
	fc_printfErrorMessage("cant get subset members as array");
	fc_deleteSubset(sidesdecay[0]);
	for (j = i; j < shape->numSides; j++){
	  fc_deleteSubset(sidesdecay[j]);
	}
	free(sidesdecay);
	return rc;
      }
      rc = fc_addArrayMembersToSubset(sidesdecay[0],num,arr);
      free(arr);
      if (rc != FC_SUCCESS){
	fc_printfErrorMessage("cant add aray members to subset");
	fc_deleteSubset(sidesdecay[0]);
	for (j = i; j < shape->numSides; j++){
	  fc_deleteSubset(sidesdecay[j]);
	}
	free(sidesdecay);
	return rc;
      }
    }
  }

  *decayedSkin = sidesdecay[0];
  free(sidesdecay);
  return FC_SUCCESS;
}


/**
 * \ingroup ElemDeath
 * \brief Get the decayed item of each of a shapes sides
 *
 * \description 
 *
 *   Given a dead element region and a shape,
 *   returns an array of subsets, which is the intersection 
 *   of the dead and the faces of each side. This is basically
 *   \ref fc_getDecayedSkin called for each side in turn, so
 *   see the text there. 
 *
 *   If you want to get thsi information for a single side,
 *   call \ref fc_getDecayedSkin and send it the single
 *   side face subset.
 *
 *   If the inner calls fail for any of the sides, the whole thing fails.
 *   
 *
 * \modifications
 *   - 03/02/2006 ACG fc_getDecayedSidesSegments created
 *   - 07/16/2006 ACG This is born from fc_getDecayedSkinSegments,
 *                    but without the segmenting.
 */
FC_ReturnCode fc_getDecayedShapeSides(
  FC_Subset deadregion, /**< input - dead region elements*/
  FC_Shape  *shape, /**< input - shape */
  FC_Subset** decayedSides /**< output - decayedSkin subsets,
			     one for each face*/
  ){
  FC_ReturnCode rc;
  FC_Subset *sidesdecay;
  int i,j;

  // default returns
  if (decayedSides!= NULL)
    *decayedSides = NULL;

  // Test input. 
  if (!fc_isSubsetValid(deadregion) || !shape || !shape->numSides ||
      !decayedSides){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  sidesdecay = (FC_Subset*)malloc((shape->numSides) * sizeof(FC_Subset));
  if (!sidesdecay){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }

  for (i = 0; i < shape->numSides; i++){
    rc = fc_getDecayedSkin(deadregion,&(shape->faces[i]),&sidesdecay[i]);
    if (rc != FC_SUCCESS){
      fc_printfErrorMessage("cant get decayedSkin of face %d",i);
      for (j = 0; j < i; j++){
	fc_deleteSubset(sidesdecay[j]);
      }
      free(sidesdecay);
      return rc;
    }
  }

  *decayedSides = sidesdecay;
  return FC_SUCCESS;
}

/**
 * \ingroup ElemDeath
 * \brief Get identifier telling if deadSubset intersects or entirely 
 *        erodes the sides of a shape.
 *
 * \description
 * 
 *   Given a dead subset, returns an array of identifiers,
 *   one for each side of a shape, telling if the deadSubset
 *   intersects or entirely erodes the sides of a shape.
 *
 *   Identifier values:
 *     -  0 if the deadSubset doesnt intersect the sides at all.
 *     -  1 if the deadSubset intersects the side
 *     -  2 if the deadSubset erodes the whole side.
 *     -  -1 in case of error
 *
 *   The dead subset can have assoc FC_AT_FACE or FC_AT_ELEM.
 *   The dead subset and shape have to be on the same mesh.
 *   This function is essentially a wrapper around 
 *   \ref fc_getSubsetDecayType. If you want this functionally
 *   for only 1 side of the Shape, then call \ref fc_getSubsetDecayType.
 *
 *   Right now, if any side fails, the whole thing fails, since
 *   you wouldnt have an empty or NULL side in a shape and therefore
 *   this is probably indicative of a more serious error in your shape.
 *
 *   While this is expected to be used for a dead region and a shape,
 *   there is nothing that checks that it is a dead region and
 *   so it can be used more generally. It is being placed in
 *   elemdeath however, because this function 
 *   was developed for a specific dead element need and using 
 *   the deadelement
 *   terminology in this function name is beneficial. 
 *
 * \todo
 *   - do we want this to succeed if any one side is bad?
 *
 * \modifications
 *   - 02/14/2006 ACG fc_getSidesDecay created, which replaced
 *                    fc_getDecayedMeshSkinSideIntersections.
 *   - 07/17/2006 ACG new version born of fc_getSidesDecay.
 */
FC_ReturnCode fc_getShapeSidesDecayType(
   FC_Subset deadSubset, /**< input - deadSubset */
   FC_Shape *shape, /**< input - ptr to shape */
   int** sidedecayflag /**<output - decay flag array */
   ){

  FC_ReturnCode rc;
  FC_AssociationType assocdead;
  int *flags;
  int usefaces = 0;
  int i;

  //default returns
  if (sidedecayflag)
    *sidedecayflag = NULL;

  //other things chekced inside
  if (!fc_isSubsetValid(deadSubset) || !shape ||
      !sidedecayflag){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  flags = (int*)malloc(shape->numSides*sizeof(int));
  if (!flags){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;  
  }

  //check assoc, mesh is checked in internal calls. 
  rc = fc_getSubsetAssociationType(deadSubset,&assocdead);
  if (assocdead != FC_AT_ELEMENT && assocdead!= FC_AT_FACE){
    fc_printfErrorMessage("deadSubset must have face or element assoc");
    return FC_INPUT_ERROR;
  }
  if (assocdead == FC_AT_FACE){
    usefaces = 1;
  } //otherwise use elems

  for (i = 0; i < shape->numSides; i++){
    FC_Subset tempSubset;
    if (usefaces){
      tempSubset = shape->faces[i];
    }else{
      tempSubset = shape->elems[i];
    }
    rc = fc_getSubsetDecayType(deadSubset,tempSubset,&flags[i]);
    if (rc != FC_SUCCESS){
      fc_printfErrorMessage("Can't get decay type for side %d",i);
      free(flags);
      return rc;
    }
  }

  *sidedecayflag = flags;

  return FC_SUCCESS;
};

    

//@}


/** \name Dead region segmenting mesh */
//-------------------------------
//@{


/**
 * \ingroup ElemDeath
 * \brief returns the segments (and their number)
 *  that result from segmenting the mesh by the input subset. 
 *
 * \description
 *
 *  Returns the segments (and their number) that result from segmenting
 *  the mesh by the input subset. If you dont want the segment array,
 *  then newSubsets and varname are optional parameters. If you do
 *  want the segments, they will be named as the varname
 *  parameter pospended with "_Seg#" where # will be the index
 *  into the array of returned subsets.
 *
 *  Returns -1 in case of error. If you pass in an empty
 *  subset -- this is not an error, unlike fc_exposedSkin.
 *
 *  One can pass in subsets that have nothing
 *  to do with dead cells as well, of course, but we anticipate
 *  that this will be most often used with a dead region. If
 *  the number of returned segments > the original number of
 *  segments then the subset has
 *  broken the mesh. If the number of returned segments = 0 then
 *  the subset has consumed the mesh, and therefore may also
 *  be considered to have broken the mesh.
 *
 *  This is the mesh equivalent of \ref fc_subsetSegmentsSubset
 *
 * \todo
 *   - test this for things other than elem
 *
 * \modifications
 *   - 8/2/2005 ACG. Created.
 *   - 8/5/05 ACG. reincarnation of fc_subsetBreaksMesh
 *   - 8/11/05 ACG made empty subset an ok case because the user
 *     prob wants to compare the number of segments to the orignal
 *     number to see if there are more. an empty subset is a valid
 *     case, esp if one is doing this over time.
 *   - 9/22/05 ACG added shared segmentation dim
 */
FC_ReturnCode fc_subsetSegmentsMesh(
   FC_Subset subset, /**< input - subset testing on */
   int shared_segdim, /**< input - shared dim for segmentation */
   char* varname, /**< input - prefix for returned segments (optional) */
   int *numSubset, /**< output - number of resulting segments */
   FC_Subset **newSubsets /**< output - array of resulting segments 
			     (optional) */
){
  FC_ReturnCode rc;
  FC_Subset tempSubset,*tempsubarray;
  //  FC_AssociationType in_assoc;
  //  int numMember, maxNumMember;
  int i;
  char* tempname;

  // default returns
  if (numSubset)
    *numSubset = -1;
  if (newSubsets)
    *newSubsets = NULL;


  // Test input. (shared_segdim checked in fc_segment)
  if (!fc_isSubsetValid(subset) || !numSubset || 
      (!newSubsets && varname) || (newSubsets && !varname))
       //couldnt get xor (^) to work
    {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // no - dont limit this to elems
  //NOTE: next two test sections copied from getExposedSkin above 
  //  rc = fc_getSubsetInfo(subset, &numMember, &maxNumMember, &in_assoc);
  //  if (rc != FC_SUCCESS)
  //    return rc;
  //
  //  if (in_assoc != FC_AT_ELEMENT && in_assoc != FC_AT_WHOLE_MESH) {
  //    fc_printfErrorMessage("Subset must have element or whole association");
  //    return FC_INPUT_ERROR;
  //  }

  fc_printfLogMessage("subsetSegmentsMesh");

  
  tempname = varname ? varname : "temp_subset";

  //this works for empty subset
  rc = fc_createSubsetComplement(subset,tempname, &tempSubset);
  if (rc != FC_SUCCESS){
    return rc;
  }

  if (newSubsets){ //even though this is more code, onyl 1 if check 
    //to determine the arg and do the delete, vs doing a temp
    //pts variable like i did for the name above and then
    //having to check again for the delete
    rc = fc_segment(tempSubset,shared_segdim, numSubset,newSubsets);
    if (rc != FC_SUCCESS){
      fc_deleteSubset(tempSubset);
      return rc;
    }
  } else {
    rc = fc_segment(tempSubset,shared_segdim, numSubset,&tempsubarray);
    if (rc != FC_SUCCESS){
      fc_deleteSubset(tempSubset);
      return rc;
    }
    for (i = 0; i < *numSubset; i++){
      fc_deleteSubset(tempsubarray[i]);
    }
    free(tempsubarray);
  }

  //clean up
  fc_deleteSubset(tempSubset);

  return FC_SUCCESS;
}


/**
 * \ingroup ElemDeath
 * \brief returns the segments (and their number)
 *  that result from segmenting a subset by another subset. 
 *
 * \description
 *
 *  Returns the segments (and their number) that result from segmenting
 *  a subset by another subset. If you dont want the segment array,
 *  then newSubsets and varname are optional parameters. If you do
 *  want the segments, they will be named as the varname
 *  parameter postpended with "_Seg#" where # will be the index
 *  into the array of returned subsets. Note: this is not
 *  the additional number of segments that result, but the total - that
 *  is,
 *
 *  Returns -1 in case of error. If the outer subset is empty or
 *  MULL this returns with an error. If the inner subset is empty
 *  or NULL this is not an error, unlike fc_exposedSkin. Both
 *  subsets must have the same association and be on the same mesh.
 *  See \ref fc_segment for restrictions on the segmentation dimension.
 *
 *  One can pass in subsets that have nothing
 *  to do with dead cells as well, of course, but we anticipate
 *  that this will be most often used with a dead region. If
 *  the number of returned segments > the original number of
 *  segments then the inner subset has broken the outer.
 *  If the number of returned segments = 0 then
 *  the inner subset has consumed the outer, and therefore may also
 *  be considered to have broken the outer.
 *
 *  This is the subset equivalent of \ref fc_subsetSegmentsMesh
 *
 * \todo
 *   - test this for things other than elem
 *
 * \modifications
 *   - 2/23/06 ACG created
 */
FC_ReturnCode fc_subsetSegmentsSubset(
   FC_Subset subset_inner, /**< input - subset that cuases the segmentation */
   FC_Subset subset_outer, /**< input - subset being segmented*/
   int shared_segdim, /**< input - shared dim for segmentation */
   char* varname, /**< input - prefix for returned segments (optional) */
   int *numSubset, /**< output - number of resulting segments */
   FC_Subset **newSubsets /**< output - array of resulting segments 
			     (optional) */
){
  FC_ReturnCode rc1, rc2;
  FC_Subset tempSubset,*tempsubarray;
  FC_Mesh mesh_inner, mesh_outer;
  FC_AssociationType assoc_inner, assoc_outer;
  int numMember_inner, numMember_outer;
  int maxNumMember_inner, maxNumMember_outer;
  int i;
  char* tempname;

  // default returns
  if (numSubset)
    *numSubset = -1;
  if (newSubsets)
    *newSubsets = NULL;


  // Test input. (shared_segdim checked in fc_segment)
  if (!fc_isSubsetValid(subset_inner) ||!fc_isSubsetValid(subset_outer) ||
      !numSubset || (!newSubsets && varname) ||(newSubsets && !varname)){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  rc1 = fc_getMeshFromSubset(subset_inner,&mesh_inner);
  rc2 = fc_getMeshFromSubset(subset_outer,&mesh_outer);
  if (rc1 != FC_SUCCESS)
    return rc1;
  else if (rc2 != FC_SUCCESS)
    return rc2;
  if (!FC_HANDLE_EQUIV(mesh_outer, mesh_inner)){
    fc_printfErrorMessage("Subsets must be on same mesh");
    return FC_INPUT_ERROR;
  }

  rc1 = fc_getSubsetInfo(subset_inner, &numMember_inner, &maxNumMember_inner, &assoc_inner);
  rc2 = fc_getSubsetInfo(subset_outer, &numMember_outer, &maxNumMember_outer, &assoc_outer);
  if (rc1 != FC_SUCCESS)
    return rc1;
  else if (rc2 != FC_SUCCESS)
    return rc2;

  if (FC_HANDLE_EQUIV(subset_outer, FC_NULL_SUBSET) || numMember_outer == 0){
    fc_printfErrorMessage("Outer subset cannot be empty or NULL ");
    return FC_INPUT_ERROR;
  }

  //  if (assoc_inner != FC_AT_ELEMENT || assoc_outer != FC_AT_ELEMENT) {
  //    fc_printfErrorMessage("Subsets must have element assoc ");
  //    return FC_INPUT_ERROR;
  //  }

  if (assoc_inner !=  assoc_outer) {
    fc_printfErrorMessage("Subsets must have same assoc ");
    return FC_INPUT_ERROR;
  }

  tempname = varname ? varname : "temp_subset";


  rc1 = fc_copySubset(subset_outer,mesh_outer,tempname, &tempSubset);
  if (rc1 != FC_SUCCESS)
    return rc1;
  if (!FC_HANDLE_EQUIV(subset_inner,FC_NULL_SUBSET) && !numMember_inner == 0){
    int* array;

    rc1 = fc_getSubsetMembersAsArray(subset_inner, &numMember_inner, &array);
    if (rc1 != FC_SUCCESS)
      return rc1;
    rc1 = fc_deleteArrayMembersFromSubset(tempSubset, numMember_inner, array);
    if (rc1 != FC_SUCCESS)
      return rc1;
    free(array);
  }

  if (newSubsets){ //even though this is more code, onyl 1 if check 
    //to determine the arg and do the delete, vs doing a temp
    //pts variable like i did for the name above and then
    //having to check again for the delete
    rc1 = fc_segment(tempSubset,shared_segdim, numSubset,newSubsets);
    if (rc1 != FC_SUCCESS){
      fc_deleteSubset(tempSubset);
      return rc1;
    }
  } else {
    rc1 = fc_segment(tempSubset,shared_segdim, numSubset,&tempsubarray);
    if (rc1 != FC_SUCCESS){
      fc_deleteSubset(tempSubset);
      return rc1;
    }
    for (i = 0; i < *numSubset; i++){
      fc_deleteSubset(tempsubarray[i]);
    }
    free(tempsubarray);
  }

  //clean up
  fc_deleteSubset(tempSubset);

  return FC_SUCCESS;
}

/**
 * \ingroup ElemDeath
 * \brief given an input subset and parameters for
 *  specifying neghbors, returns IDs of neighbors
 *  for which adding any single one individually to the subset will result in
 *  a greater segmentation of the mesh than results from
 *  segmenting based on the orginal subset
 *
 * \description
 *
 *  given an input subset and parameters for
 *  specifying neghbors, returns IDs of neighbors
 *  for which adding any single one individually to the subset will result in
 *  a greater segmentation of the mesh than results from
 *  segmenting based on the orginal subset.
 *
 *  Returns -1 in case of error. Passing in an empty subset returns
 *  -1.
 *
 *  Notes: 
 *   - No sense passing in anything more than level 1 since
 *     we are considering no element combos at this point.
 *   - This uses recursive innards that test on big sets of elem
 *     first and only narrow it down to smaller set if the big
 *     sets segment first. i.e., It is assumed that if a whole
 *     group doesnt result in any different segmenting, that the
 *     items of the group taken individually wont rest in different
 *     segmenting. There is the possibility that there could be
 *     offsetting effects that occur within the bigger set that
 *     will cuase this to be a faulty assumption, particularly
 *     when the dead region is near the boundary of the
 *     mesh. However, if this
 *     is used primarily for cases where we are looking for the
 *     mesh to be just about to break, then we expect that the
 *     num of elements that will satisfy the condition is small
 *     and so this should be ok. 
 *     
 *
 * \todo
 *
 *  - prob need a version that doesnt make the assumption that
 *    the recusion implies - esp for dead regions near
 *    boundary of mesh.
 *  - should this return IDS if it segments into fewer (ie
 *    subsumes the mesh, or causes to segments to combine ?). Note
 *    may cause problems with considering things as groups as effects
 *    might counter each other. 
 *  - may want a version that does combos of elem. the recursive
 *    innards, would not support that.
 *  - should we return num of resulting segments ?
 *  - should we make a time version or no ?
 *  - note: will need to upate if change return vals of recurive innard
 *  - do we want to explicity test the shareddim earlier in the call, and
 *    the unittest ratehr than them just being caught in the innards ?
 *
 * \modifications
 *   - 8/24/2005 ACG. Created.
 *   - 8/31/2005 ACG recusive version
 *   - 9/20/2005 ACG changed so uses new recursive innards taht retrun segmenting nodes
 *     instead of non-segmenting nodes so dont have to invert the list at the end.
 *   - 9/20/2005 ACG name change to reflect that will only be dealing with greater segmentation
 *   - 9/22/2005 ACG added segmentation dim
 *   - 10/9/05  ACG limited it to neighbor depth level 1
 */
FC_ReturnCode fc_subsetPlusNeighborGreaterSegmentationMesh(
  FC_Subset subset, /**< input - subset testing on */
  int shared_neighbordim, /**< input - minimum dimensionality of shared part of neighbors */
  int shared_segdim, /**< input - shared dim for segmentation */
  int *numNbr, /**< output - number of neigbors that will individually cause greater segmentation */
  int **nbrIDs /**< output - IDS of neighbors that will individually cause greater segmentation */
  ){
  FC_ReturnCode rc;
  FC_AssociationType in_assoc;
  int numMember, maxNumMember;
  int numSubset, natNbr, *natNbrIDs, nmaxnsegids, *maxnsegids;
  int i;
  int level = 1;

  // default returns
  if (numNbr)
    *numNbr = -1;
  if (nbrIDs)
    *nbrIDs = NULL;


  // Test input (segment will check shared_segdim)
  if (!fc_isSubsetValid(subset) ||  level < 1 || shared_neighbordim < 0 || 
      shared_neighbordim > 2 || !numNbr || !nbrIDs) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  rc = fc_getSubsetInfo(subset, &numMember, &maxNumMember, &in_assoc);
  if (rc != FC_SUCCESS)
    return rc;

  if (in_assoc != FC_AT_ELEMENT && in_assoc != FC_AT_WHOLE_MESH) {
    fc_printfErrorMessage("Subset must have element or whole association");
    return FC_INPUT_ERROR;
  }


  if (numMember == 0) {
    fc_printfErrorMessage("Cannot handle empty subset");
    return FC_INPUT_ERROR;
  }


  fc_printfLogMessage("subsetPlusNeighborsSegmentsMesh");

  //first determine the orignal segmentation
  rc = fc_subsetSegmentsMesh(subset,shared_segdim,NULL,&numSubset,NULL);
  if (rc != FC_SUCCESS){
    fc_printfErrorMessage("cant get original segmentation");
    return rc;
  }

  //now get the neighbors 
  rc = fc_getSubsetNeighbors(subset,level,shared_neighbordim,
			      &in_assoc,
			      &natNbr,&natNbrIDs);
  if (rc != FC_SUCCESS){
    fc_printfErrorMessage("cant get neighbors");
    return rc;
  }
  
  if (natNbr == 0){ //no neighbors - treat this as valid
    *numNbr = 0;
    return FC_SUCCESS;
  }

  rc = _fc_recursive_IndividuallySegmentingNodes(subset,numSubset,shared_segdim,
				    natNbr,natNbrIDs,
				    &nmaxnsegids, &maxnsegids);
  if (rc != FC_SUCCESS){
    fc_printfErrorMessage("cant get recursive segmenting nodes");
    if (natNbrIDs) free(natNbrIDs);
    return rc;
  }

  //clean up neighbor arrays
  if (natNbrIDs) free(natNbrIDs);

  //return vals from recursive call are all the vals we want to return
  //dont have to check them anymore

  //now set real return values
  *numNbr = nmaxnsegids;
  if (nmaxnsegids > 0){
    *nbrIDs = (int*)malloc(nmaxnsegids*sizeof(int));
    if (*nbrIDs == NULL){ //no lcean up in case of memory error
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      return FC_MEMORY_ERROR;
    }

    for (i = 0; i < nmaxnsegids;i++){
      (*nbrIDs)[i] = maxnsegids[i];
    }
  }
  if (maxnsegids) free(maxnsegids);

  return FC_SUCCESS;
}


//@}


/** \name Characterizations of Dead Element Regions */
//-------------------------------
//@{


/**
 * \ingroup ElemDeath
 * \brief Get tear length.
 *
 * \description 
 *
 *   This function is just for convenience - it calls fc_getExposedSkin()
 *   and then fc_getSubsetDiameter() (or fc_getDisplacedSubsetDiameter() if
 *   the displacement variable is provided).
 *
 * \todo Unit test this.
 *
 * \modifications
 *   - 7/6/2005 WSD. Created.
 */
FC_ReturnCode fc_calcTearLength(
  FC_Subset subset, /**< input - A subset of "dead" elements */
  FC_Subset* opt_meshSkin,   /**< input - (optional) pointer to mesh skin 
			      (or NULL) */
  FC_Variable* displ_coords, /**< input - (optional) pointer to the owning
				mesh's coordinate displacements (or NULL) */
  double* length    /**< ouput - diameter of the vertex set */
) {
  FC_ReturnCode rc;
  FC_Subset exposed;
  rc = fc_getExposedSkin(subset, opt_meshSkin, &exposed);
  if (rc != FC_SUCCESS) {
    *length = -1;
    return rc;
  }
  if (displ_coords)
    rc = fc_getDisplacedSubsetDiameter(exposed, *displ_coords, length, NULL,
				       NULL);
  else
    rc = fc_getSubsetDiameter(exposed, length, NULL, NULL);
  fc_deleteSubset(exposed);
  return rc;
}


//@}


/**
 * \ingroup PrivateElemDeath
 * \brief helper function that recursively determines 
 *  node ids that individually will result in greater segmentation than 
 *  some comparison value
 *  
 *
 * \description
 * 
 *  helper function that recursively determines 
 *  node ids that will individually result in greater segmentation than 
 *  some comparison value
 * 
 * Note:
 *   - This uses recursion to test on big sets of elem
 *     first and only narrow it down to smaller set if the big
 *     sets segment first. i.e., It is assumed that if a whole
 *     group doesnt result in a greater segmenting, that the
 *     items of the group taken individually wont result in a greater
 *     segmenting. There is the possibility that there could be
 *     offsetting effects that occur within the bigger set that
 *     will cuase this to be a faulty assumption, particularly
 *     when the dead region is close to the boundary of the
 *     mesh. However, if this
 *     is used primarily for cases where we are looking for the
 *     mesh to be just about to break, then we expect that the
 *     num of elements that will satisfy the condition is small
 *     and so this should be ok. 
 *
 * \todo
 * - will we want less segmentation too.
 *
 * \modifications
 *  - 9/20/05 ACG created in response to thinking that its better to
 *    return the segmenting rather than the nonsegmenting nodeset so
 *    dont have to do inversion of the set later. replaces 
 *    _fc_maximalNonSegmentingNodeset
 *  - 9/22/05 ACG added shared segdim
 *  - 10/09/05 ACG explicitly check that 1 element set doesnt result in 
 *    less segmentation (rather than just equal). will continue to
 *    loop for larger sets resulting in less segmentation (for instance
 *    if whole set of neighbors subsumes the live regions, still
 *    want to delve down into individual bits)
 */
FC_ReturnCode _fc_recursive_IndividuallySegmentingNodes(
  FC_Subset subset,
  int compare,
  int shared_segdim,
  int numin,
  int* inset,
  int *numout,
  int** outset
){
  FC_ReturnCode rc;
  FC_Mesh mesh;
  FC_Subset testSubset;
  int numSegment;
  int i;

  //defaults
  *numout = 0;
  if (outset)
    *outset = NULL;

  // Test input (shared_segdim checked in fc_segment)
  if (!fc_isSubsetValid(subset) || !numout || !outset ||
      numin < 0 || (!inset  && numin != 0)){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }


  fc_printfLogMessage("recursive segmenting nodes");
  
  //  printf("\tChecking set (%d) ",numin);
  //  for (i = 0; i < numin; i++){
  //    printf("%d ", inset[i]);
  //  }
  //  printf("\n");
  
  if (numin == 0){
    return FC_SUCCESS;
  }

  //now lets get a copy of the subset to play with
  rc = fc_getMeshFromSubset(subset,&mesh);
  if (rc != FC_SUCCESS){
    fc_printfErrorMessage("cant get mesh");
    return rc;
  }

  rc = fc_copySubset(subset,mesh,"test_subset",&testSubset);
  if (rc != FC_SUCCESS){
    fc_printfErrorMessage("cant copy subset");
    return rc;
  }
  
  //now see if adding if the whole set does anything
  rc = fc_addArrayMembersToSubset(testSubset,numin,inset);
  if (rc != FC_SUCCESS){
    fc_printfErrorMessage("cant add ids to subset");
    return rc;
  }


  rc =  fc_subsetSegmentsMesh(testSubset,shared_segdim,NULL,
			      &numSegment,NULL);
  if (rc != FC_SUCCESS){
    fc_printfErrorMessage("failed to segment");
    return rc;
  }

  rc = fc_deleteSubset(testSubset);

  //compare is passed in so we dont have to recalc it each time, but it better
  //be correct
  if (numSegment == compare){ 
    //    printf("elements as a group do not divide any more (or less) (%d). maximal recursion - returning 0 elem\n",numSegment);
    //then adding the vals doesnt cause any greater segmentation
    *numout = 0;
  } else {
    if (numin == 1){
      if (numSegment <= compare){
	//     printf("single element %d doesnt divide more. maximal recursion - returning 1 elem\n",inset[0]);
	*numout = 0;
      }else {
	//     printf("single element %d divides. maximal recursion - returning 1 elem\n",inset[0]);
	*numout = 1;
	*outset = (int*)malloc(sizeof(int));
	if (*outset == NULL){
	  fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
	  return FC_MEMORY_ERROR;
	}
	(*outset)[0] = inset[0];
      }
    }else{
      int idiv = numin/2;
      int *seta, *setb;
      int *reta, *retb, nreta, nretb;

      //      printf("this does divide - recursing\n");

      seta = malloc(idiv*sizeof(int));
      setb = malloc((numin-idiv)*sizeof(int));
      if ((seta == NULL && idiv > 0) || (setb == NULL && numin-idiv > 0)){
	fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
	return FC_MEMORY_ERROR;
      }
      for (i = 0; i < idiv; i++){
	seta[i] = inset[i];
      }
      for (i = 0; i < numin-idiv; i++){
	setb[i] = inset[i+idiv];
      }

      //      printf("\tcalling set A\n");
      rc = _fc_recursive_IndividuallySegmentingNodes(subset,compare,shared_segdim,
				       idiv,seta,&nreta,&reta);
      if (rc != FC_SUCCESS){
	fc_printfErrorMessage("failed to get recursive segmenting nodes");
	free(seta);
	free(setb);
	return rc;
      }

      //      printf("\tcalling set B\n");
      rc = _fc_recursive_IndividuallySegmentingNodes(subset,compare,shared_segdim,
				       numin-idiv,setb,&nretb,&retb);
      if (rc != FC_SUCCESS){
	fc_printfErrorMessage("failed to get recursive segmenting nodes");
	free(seta);
	free(setb);
	return rc;
      }

      free(seta);
      free(setb);
      
      *numout = nreta + nretb;
      if (*numout > 0){
	*outset = (int*)malloc(*numout*sizeof(int));
	if (*outset == NULL){
	  //no cleanup in case of memory error
	  fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
	  return FC_MEMORY_ERROR;
	}
	for (i = 0; i < nreta; i++){
	  (*outset)[i] = reta[i];
	}
	for (i = 0; i < nretb; i++){
	  (*outset)[i+nreta] = retb[i];
	}
	//      if (*numout != 0) printf("\trecur: %d in %d out\n",numin,*numout);
	if (reta) free(reta); //they will be NULL otherwise
	if (retb) free(retb);
      }
    }
  }

  return FC_SUCCESS;
}



