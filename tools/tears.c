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
 * \file tears.c
 * \brief Report tear characterizations on given dataset.
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/tools/tears.c,v $
 * $Revision: 1.31 $
 * $Date: 2006/11/10 19:38:18 $
 *
 * \description
 *    Usage: tears [options] dataset [var_name]
 *
 *    A tear is a region of dead elements. This routine reports for the last
 *    time step:
 *
 *      - The number of cells in a tear.
 *
 *      - The volume of the tear on the original mesh (dead elements may not
 *        behave proproperly so dead element volume cannot be calculated).
 *
 *      - A characteristic tear length which is the largest distance between
 *        any two vertices that define the surface of the dead region. If a
 *        displacement variable is provided, it will report the tear length of 
 *        the displaced geometry in addition to the tear length in the time 
 *        zero geometry.
 *
 *      - (Optional) Damage-weighted volumes. The user has to supply two
 *        "death" criteria, the first one for dead elements and the second
 *        for damaged elements. The tears are the unions of these two criteria.
 *        If the two sets overlap, the overlap elements are considered to
 *        be in the dead set.
 *
 *
 *   SHAPE INFORMATION
 *
 *     <DL>
 *     <DT> Shape Information: </DT>      
 *          <DD>
 *          <DL>
 *          <DT>For each mesh with a tear, shape information is printed.
 *          Example:</DT>
 *           <DD>
 *
 *           Shape (0:0) (5 sides):                         \n
 *           Thin shape approx 1: (0 4)                     \n
 *           Thin shape approx 2: (0 4)                     \n
 *           Screw shape: NA                                \n
 *           Ends: NA                                       \n
 *           Sides in area order (descending) ( 0 4 1 3 2)  \n
 *                                                          \n
 *           where Shape is labeled by (meshid:shapeid)     \n   
 *          </DD>
 *          </DL>
 *          </DD>
 *
 *    <DT> Shape-tear information: </DT>
 *         <DD>
 *         <DL> 
 *         <DT> For each tear, shape information is printed.
 *         Example:</DT>
 *            <DD>
 *
 *            Subtear 14 Intersections with Shape (0:0) (TUNNEL,NONADJSIDES,MAJOR): ( 17) ( 5)
 *            </DD>
 *         <DT> Summary info on a per shape basis is also provided, with the
 *           types of tears and (major/overall) ratio of each.
 *           Example:</DT>
 *           <DD>
 *
 *           Tears 18: BREAKS (0/0) TUNNELS (10/0) PITS (8/13)
 *           </DD>
 *         </DL>
 *         </DD>
 *     </DL>
 *   
 *
 *
 *    SHAPE DETAILS
 *
 *    A shape is a contiguous, distinct starting live region in a mesh.
 *    For instance, a mesh made of N separate screws will have N
 *    shapes. A shape is distinguished by its sides and the adjacency
 *    information of those sides. Only skin information is kept for
 *    each shape (no info on innards).
 *
 *    The number of sides is very dependent on the angle chosen (optional
 *    input argument) to distinguish which faces get merged into which side.
 *    If it is set too large, shapes with gradually curved
 *    edges (like a block wiht curved edges) may come out with too few sides.
 *    On the other hand, if the angle is set too low, a curved side in a
 *    shape (like a screw) may come out a large number of very small sides.
 *    There is currently no writeout to see what constitutes a side,
 *    though some meaningful information is given as described below.
 *
 *    Some approximations are made which attempt to distinguish sides
 *    of interest:
 *      - thin shape: tries to find sides of interest in a shape that
 *         is thin in some dimension
 *             <ul>
 *             <li>approx 1 returns the largest side and its largest
 *               non-adj side. this is fairly good for finding
 *               the two largest surfaces of something like a can
 *               with a curved inner and outer shell
 *             <li> approx 2 returns the largest side that has an
 *                  opposing side. An opposing side is defined as
 *                  one that is non-adj to the original side with
 *                  an average normal that points in a roughly opposing
 *                  direction  -- within OPPOSINGSIDEANGLE degrees of
 *                  180 (currently a constant 10 degrees) -- to that of
 *                  the original side.
 *             </ul>
 *      - screw: traditional 5 sided screw, listed in order from
 *               the top of the head down to the small base
 *      - ends: any side only adjacent to one other side
 *      - descending area order: sides listed in decreasing area
 *    NA is used to mark any case where the approx is not valid.
 *
 *    Information is given on each tear as to the number of intersections it
 *    has with a shape and which sides are intersected. It is further catagoried
 *    by:
 *     types: (these catagories are shape independent)
 *     - BREAK - breaks the shape into more than one piece
 *     - TUNNEL - intersects the shape (skin) in more than 1 place 
 *     - PIT - intersects in a single place
 *     subtypes: (these catagories are shape dependent)
 *     - SINGLESIDE - intersects only a single side (but may be multiple times)
 *     - NONADJSIDE - intersects at least two non adj sides
 *     - ADJSIDES - intersects only adj sides (but may be in multiple places)
 *     class: (these catagories are shape and thinshape assumption dependent)
 *     - MAJOR - intersects at least one major side fulfilling either thin shape assumption
 *     - MINOR - intersects no majors sides
 *
 *     Note that:
 *     - The shape dependent catagories can change depending
 *       on the shape angle - that is, the number and adj of sides may change
 *       if the angle change, so a SINGLESIDE PIT under one angle may become
 *       a NON-ADJSIDE PIT.
 *     - A BREAK segments the shape into more than one piece, so a dead chunk 
 *       taken out of a shape will be classified as a PIT, not a BREAK. This
 *       means:
 *       <ul>
 *       <li> A ring shape with a chunk taken clear out of it will be classified 
 *            by the algorithm as a PIT, even though a human might consider
 *            that as a BREAK.
 *       <li> Similarly a chunk taken out of a mesh may effectively act like
 *            a TUNNEL to a human if that dead region of the mesh bumps up
 *            against another mesh. That is, there now is a hole in the overall
 *            object. Again these are only classified by the algorithm as PITs.
 *     - I only consider cases where the tear intersects the shape - 
 *       totally internal pits or breaks will not be discovered (although
 *       I do use the innards for the break calculation, that calc is never done
 *       unless the tear intersects the shape (skin).
 *
 *
 * \todo 
 *   - ?Do feature tracking ?
 *   - shape stuff not unittested
 *   - shape stuff not valground
 *   - See ACG notes in the code about some thoughts on if should move some things around based on
 *             storage implcations. Think about if this should also affect the FCLib interfaces
 *             (how unweildly is it to do some of this stuff?)   
 *
 * \modifications
 *   - 07/06/2005 WSD, created.
 *   - 05/25/2006 Changed to support multiple simulataneous definitions of
 *       dead elements.
 *   - 06/07/2006 WSD, add capability to specify specific meshes on command
 *     line. 
 *   - 06/08/2006 WSD added damage weighted volume calcs.
 *   - 07/19/2006 WSD added flag to do a specific time step.
 *   - 08/26-27/2006 ACG added optional shape information
 *   - 09/08/2006 ACG fixed break logic - actually need the innards to
 *                determine a break. 
 */

// misc notes - going to combine tears  based on undisplaced geometry

#include <string.h>
#include <stdio.h>
#include "fc.h"
#include "fcP.h" // temporary until writeBB stuff gets made public

enum{ //defualt paramters for shape, now current params till i
  //make them aruments
  SHAPE_ANGLE = 40,
  OPPOSINGSIDES_ANGLE = 10,
  SHARED_DIM = 0,   
};

enum{ //these have to do with the intersection with the skin and
  //are not dependent on the definition of the shape sides
  BREAK,TUNNEL,PIT,
};

enum{ // these are shape dependent
  SINGLESIDE, NONADJSIDES, ADJSIDES,
};

typedef struct{
  char* name;
  char* op;
  double val;
} deathvarinfo;

typedef struct{
  int meshID;
  int shapeID;
  int type; //tunnel, broken, pit
  int subtype; //single, nonadj, adj
  int major; // if it involves major thinshape sides
  int numIntersections;
  FC_SortedIntArray* intersectedSides;
} ShapeIntersection;

typedef struct {
  int numDim;
  int numCell;
  double volume;              // region area/volume (can't do displaced)
  double diameter;            // region diameter (can't do displaced)
  double exp_diameter;        // exposed diameter
  double displ_exp_diameter;  // displaced exposed diameter
  FC_Coords lowers;           // region bounding box
  FC_Coords uppers;
  FC_Coords exp_lowers;       // exposed bounding box
  FC_Coords exp_uppers;
  FC_Coords displ_exp_lowers; // displaced exposed bounding box
  FC_Coords displ_exp_uppers;
  // damage stuff - this is the minimum set needed to calculate everything
  // (calculation is done during print_tear_characterization())
  int dead_numCell;
  int damaged_numCell;
  double dead_vol;
  double dead_logvol;
  double damaged_vol;
  double damaged_vol_x_damage;
  double damaged_vol_x_logdamage;
  double damaged_log_vol_x_damage;
} Tear_Characterization;

typedef struct {
  int tearID;
  FC_Subset region;  // dead element region - coords only valid of not displ'd
  FC_Subset exposed; // exposed surface = part of dead region shared w/ mesh
  int meshID;  // into meshes & meshnames array
  int stepID;  // this is silly?
  Tear_Characterization data;
  ShapeIntersection* shape_intersection; 
} Tear;

typedef struct {
  int tearID;
  int numTear;
  Tear** tears;
  Tear_Characterization data;
} SuperTear;

typedef struct{
  int meshID; //wish i didnt have repeat these ids, but its easiest
  int shapeID;
  int numSides; 
  int *nonadjorder; //has 2 members
  int *opposingorder; //has 2 members
  int *screworder;  //has numSides members
  int *areaorder; //has numSides members
  int numEnds;
  int *ends; //has numEnds members
}shape_order;

typedef struct{
  //now that we need the innards for breaking, should i delay the building
  //of the shape till later? right now keeping both the innards and all the
  //sides info until after do all the shapeintersection calcs on that mesh
  int meshID; //will move ids into shape later
  int shapeID;
  FC_Shape *shape;
  FC_Subset shapeinnards; //need this only for breaking
  shape_order *order;
}shape_w_order;

static int parseDeathVarSet(int argc, char** argv, int *idx, int *numvars,
			    deathvarinfo** deathvars) {
  long lcheck;
  double dcheck;
  int nvars;

  char* end_ptr;
  int i,j;

  if (*numvars !=0){
    fc_printfErrorMessage("Invalid deathvar syntax: already have death vars");
    return FC_ERROR;
  }
   

  i = *idx;
  if (argc < *idx+1){
    fc_printfErrorMessage("Invalid deathvar syntax: cant get num vars");
    return FC_ERROR; //cant get the num of vars
  }

  lcheck =  strtol(argv[i], &end_ptr,10);
  if (*end_ptr == '\0'){ //its a number
    nvars = (int)lcheck;
  }else { //not a number
    return FC_ERROR;
  }

  if (nvars == 0){
    fc_printfErrorMessage("Invalid deathvar syntax: no death vars");
    return FC_ERROR;
  }

  if (argc < i+3*nvars+1){
    fc_printfErrorMessage("Invalid deathvar syntax: not enough args for num vars");
    return FC_ERROR;
  }

  *deathvars = (deathvarinfo*)malloc(nvars*sizeof(deathvarinfo));
  if (!deathvars){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }

  i++;
  for(j = 0; j < nvars; j++){
    (*deathvars)[j].name = argv[i++];
    (*deathvars)[j].op = argv[i++];
    dcheck =  strtod(argv[i++], &end_ptr);
    if (*end_ptr == '\0'){ //its a number
      (*deathvars)[j].val = dcheck;
    } else{
      fc_printfErrorMessage("Invalid deathvar syntax: non-numerical value for cutoff");
      return FC_ERROR;
    }
  }

  *numvars = nvars;		     
  *idx = i-1;

  return FC_SUCCESS;
}

// Sort by diameter - in descending order
// FIX? sort by volume? by displaced exposed diameter?
static int cmpSuperTearsByLength(const void* n, const void* m) {
  double a = ((const SuperTear*)n)->data.diameter;
  double b = ((const SuperTear*)m)->data.diameter;
  if (a > b)
    return -1;
  else if (a < b)
    return 1;
  else
    return 0;
}

// Expand the tears array - make sure there is enough
// Note, if the Tear struct get's any more members, may want to change from 
// Tear* to Tear** so the array takes less space when emptyish.
static FC_ReturnCode grow_tears(int numNewTear, int numCurrentTear, 
				int* maxNumTear, Tear** tears) {
  Tear* temp_tears;
  int temp_maxNumTear;

  // internal, so no input testing

  // no new tears, do nothing
  if (numNewTear < 1)
    return FC_SUCCESS;

  // if we are going, to run out of room, expand the array
  if (numNewTear + numCurrentTear > *maxNumTear) {
    temp_maxNumTear = 2*(numNewTear + numCurrentTear);
    temp_tears = realloc(*tears, temp_maxNumTear*sizeof(Tear));
    if (!temp_tears) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      return FC_MEMORY_ERROR;
    }
    *maxNumTear = temp_maxNumTear;
    *tears = temp_tears;
  }
  // else, do nothing

  // done
  return FC_SUCCESS;
}



static void printShapeIntersection(FILE* output_file, int tearID,
				   ShapeIntersection *si){
  int i,j;

  if (!si) return;
  
  fprintf(output_file,"  Subtear %4d Intersections with Shape (%d:%d) ",tearID,
	  si->meshID,si->shapeID);

  switch(si->type){
  case (BREAK):
    fprintf(output_file,"(BREAK,");
    break;
  case (TUNNEL):
    fprintf(output_file,"(TUNNEL,");
    break;
  case (PIT):
    fprintf(output_file,"(PIT,");
    break;
  default:
    //shouldnt happen
    break;
  }
  switch(si->subtype){
  case (SINGLESIDE):
    fprintf(output_file,"SINGLESIDE,");
    break;
  case (NONADJSIDES):
    fprintf(output_file,"NONADJSIDES,");
    break;
  case (ADJSIDES):
    fprintf(output_file,"ADJSIDES,");
    break;
  }
  if (si->major){
    fprintf(output_file,"MAJOR)");
  }else{
    fprintf(output_file,"MINOR)");
  }

  fprintf(output_file,": ");

  for (i = 0; i < si->numIntersections; i++){
    fprintf(output_file,"(");
    for (j = 0; j < si->intersectedSides[i].numVal; j++){
      fprintf(output_file," %d",si->intersectedSides[i].vals[j]);
    }
    fprintf(output_file,") ");
  }
  fprintf(output_file, "\n");
}

// print tear data
// use "area" or "volume" where appropriate
// treat the "total" object differently
static void print_tear_characterization(FILE* output_file,
					Tear_Characterization *data,
					int isTotal, int topoDim,
					int doDispl, int doDamageWeight,
					double dead_damage) {
  char* volTypeStr[3] = { "area+volume", "area", "volume" };
  int strID = topoDim-1;
  if (topoDim < 0)
    strID = 0;

  fprintf(output_file, "  numCell = %d\n", data->numCell);
  if (topoDim > 1) 
    fprintf(output_file, "  region %s = %g\n", volTypeStr[strID],
	    data->volume);
  fprintf(output_file, "  region diameter  = %g\n", data->diameter);
  fprintf(output_file, "  exposed diameter = %g\n",
	  data->exp_diameter);
  if (doDispl) {
    fprintf(output_file, "  displ exposed diameter = %g\n", 
	    data->displ_exp_diameter);
  }
  // don't print bb is this is a total, they have no meaning
  if (!isTotal) {
    fprintf(output_file, "  region bb  = [ %g, %g, %g ] - [ %g, %g, %g ]\n",
	    data->lowers[0], data->lowers[1],
	    data->lowers[2], data->uppers[0],
	    data->uppers[1], data->uppers[2]);
    fprintf(output_file, "  exposed bb = [ %g, %g, %g ] - [ %g, %g, %g ]\n",
	    data->exp_lowers[0], data->exp_lowers[1],
	    data->exp_lowers[2], data->exp_uppers[0],
	    data->exp_uppers[1], data->exp_uppers[2]);
    if (doDispl) {
      fprintf(output_file, "  displ exposed bb = [ %g, %g, %g ] - [ %g, %g, %g ]\n",
	      data->displ_exp_lowers[0], 
	      data->displ_exp_lowers[1],
	      data->displ_exp_lowers[2], 
	      data->displ_exp_uppers[0],
	      data->displ_exp_uppers[1], 
	      data->displ_exp_uppers[2]);
    }
  }
  if (doDamageWeight) {
    double dead_vol_x_damage = data->dead_vol*dead_damage;
    double dead_vol_x_logdamage = data->dead_vol*log10(dead_damage);
    double dead_log_vol_x_damage = data->dead_logvol + 
      data->dead_numCell*log10(dead_damage);
    fprintf(output_file, "  --- Damage weighted %ss ---\n", volTypeStr[strID]);
    fprintf(output_file, "  dead numCell = %d\n", data->dead_numCell);
    fprintf(output_file, "  dead %s = %g\n", volTypeStr[strID], 
	    data->dead_vol);
    fprintf(output_file, "  dead log(%s) = %g\n", volTypeStr[strID], 
	    data->dead_logvol);
    fprintf(output_file, "  dead %s x damage = %g\n",volTypeStr[strID],
	    dead_vol_x_damage);
    fprintf(output_file, "  dead %s x log(damage) = %g\n",volTypeStr[strID],
	    dead_vol_x_logdamage);
    fprintf(output_file, "  dead log(%s x damage) = %g\n",volTypeStr[strID],
	    dead_log_vol_x_damage);
    fprintf(output_file, "  damaged numCell = %d\n", data->damaged_numCell);
    fprintf(output_file, "  damaged %s = %g\n", volTypeStr[strID],
	    data->damaged_vol);
    fprintf(output_file, "  damaged %s x damage = %g\n", volTypeStr[strID],
	    data->damaged_vol_x_damage);
    fprintf(output_file, "  damaged %s x log(damage) = %g\n",
	    volTypeStr[strID], data->damaged_vol_x_logdamage);
    fprintf(output_file, "  damaged log(%s x damage) = %g\n",
	    volTypeStr[strID], data->damaged_log_vol_x_damage);
    fprintf(output_file, "  total %s x damage = %g\n", volTypeStr[strID],
	    dead_vol_x_damage + data->damaged_vol_x_damage);
    fprintf(output_file, "  total %s x log(damage) = %g\n", volTypeStr[strID], 
	    dead_vol_x_logdamage + data->damaged_vol_x_logdamage);
    fprintf(output_file, "  total log(%s x damage) = %g\n", volTypeStr[strID],
	    dead_log_vol_x_damage + data->damaged_log_vol_x_damage);
  }
}


static void printShapeOrder(FILE* output_file,shape_order *order){
  if (!order) return;

  fprintf(output_file,"Shape (meshid:shapeid) (%d:%d) (numSides %d):\n",
	  order->meshID,order->shapeID,order->numSides); 
  fprintf(output_file,"  Thin shape approx 1: ");
  if (order->nonadjorder){
    fprintf(output_file,"(%d %d)\n",
           order->nonadjorder[0],order->nonadjorder[1]);
  }else{
    fprintf(output_file,"NA\n");
  }

  fprintf(output_file,"  Thin shape approx 2: ");
  if (order->opposingorder){
    fprintf(output_file,"(%d %d)\n",
           order->opposingorder[0],order->opposingorder[1]);
  }else{
    fprintf(output_file,"NA\n");
  }

  fprintf(output_file,"  Screw shape: ");
  if (order->screworder){
    int i;
   fprintf(output_file,"(");
    for (i = 0 ; i < order->numSides; i++){
      fprintf(output_file," %d",order->screworder[i]);
    }
    fprintf(output_file,")\n");
  }else{
    fprintf(output_file,"NA\n");
  }

  fprintf(output_file,"  Ends: ");
  if (order->ends){
    int i;
    fprintf(output_file,"(");
    for (i = 0 ; i < order->numEnds; i++){
      fprintf(output_file," %d",order->ends[i]);
    }
    fprintf(output_file,")\n");
  }else{
    fprintf(output_file,"NA\n");
  }

  if (order->areaorder){
    int i;
    fprintf(output_file,"  Sides in area order (descending) (");
    for (i = 0 ; i < order->numSides; i++){
      fprintf(output_file," %d",order->areaorder[i]);
    }
    fprintf(output_file,")\n");
  }
}


static void freeShapeOrder(shape_order *order){
  if (!order) return;

  order->meshID = -1;
  order->shapeID = -1;
  if (order->nonadjorder)  free(order->nonadjorder);
  if (order->screworder) free(order->screworder);
  if (order->opposingorder) free(order->opposingorder);
  if (order->ends) free(order->ends);
  if (order->areaorder) free(order->areaorder);
}

static void freeShapeWOrder_ShapeOnly(shape_w_order *wshape){
  if (!wshape) return;

  wshape->meshID = -1;
  wshape->shapeID = -1;
  fc_deleteSubset(wshape->shapeinnards);
  fc_freeShape(wshape->shape);
}


static void freeShapeWOrder(shape_w_order *wshape){
  freeShapeOrder(wshape->order);
  freeShapeWOrder_ShapeOnly(wshape);
}


static void createShapeWOrder(int mid, int sid,
			      FC_Shape* shape, 
			      FC_Subset innards,
			      shape_w_order* wshape){
  //  FC_ReturnCode rc;

  double angle = OPPOSINGSIDES_ANGLE;
  shape_order *order = (shape_order*)malloc(sizeof(shape_order));
  if (!order){
    fc_exitIfError(FC_MEMORY_ERROR);
  }
  
  wshape->meshID = mid;
  wshape->shapeID = sid;
  wshape->shape = shape;
  wshape->order = order;
  wshape->shapeinnards = innards;
  wshape->order->meshID = wshape->meshID;
  wshape->order->shapeID = wshape->shapeID;
  wshape->order->numSides = wshape->shape->numSides;

  //we dont care if these orders fail
  fc_createLargeAndNonAdjacentSidesOrder(shape,&(wshape->order->nonadjorder));
  fc_createLargeAndOpposingSidesOrder(shape,angle,
				      &(wshape->order->opposingorder));
  fc_createScrewShapeOrder(shape,&(wshape->order->screworder));
  fc_getShapeEnds(shape,&(wshape->order->numEnds),
		  &(wshape->order->ends));
  fc_createDescendingAreaOrder(shape,&(wshape->order->areaorder));
}

static FC_ReturnCode freeShapeIntersection(ShapeIntersection *si){
  int i;

  if (!si) return FC_SUCCESS;

  for (i = 0; i < si->numIntersections; i++){
    fc_freeSortedIntArray(&(si->intersectedSides[i]));
  }

  si->numIntersections = 0;
  free(si->intersectedSides);

  return FC_SUCCESS;
}

static FC_ReturnCode calcShapeIntersection(int numShapes,
					   shape_w_order *wshapes,
					   int shareddim,
					   Tear *tear){
  FC_ReturnCode rc;
  ShapeIntersection *si = NULL;
  int i,j,m,n;


  //each segment can intersect only 1 shape
  //if its entirely interior to a shape, then
  //it will intersect 0 shapes

  for (i = 0; i < numShapes; i++){
    FC_Subset decayedSkin, *livesegments;
    int numlive = 0;
    FC_Subset *decayedSidesSegments;
    int numDecayedSidesSegments;
    FC_Shape *currshape = wshapes[i].shape;
    shape_order *order = wshapes[i].order;
    FC_Subset currinnards = wshapes[i].shapeinnards;

    int numMem;
    int broken = 0;

    rc = fc_getDecayedShapeSkin(tear->region,currshape,
				&decayedSkin);
    if (rc != FC_SUCCESS){
      fc_printfErrorMessage("Cant get decayedSkin. skipping");
      continue;
    }

    if (FC_HANDLE_EQUIV(decayedSkin,FC_NULL_SUBSET)){
      //i hate having to check this
      printf("\tNULL SUBSET continuing\n");
      continue;
    }

    rc = fc_getSubsetNumMember(decayedSkin,&numMem);
    fc_exitIfErrorPrintf(rc, "can't get subset members");
    if (!numMem){
      //nothing in this shape
      fc_deleteSubset(decayedSkin);
      continue;
    }


    //segment should work if its an empty set
    //segment the dead skin
    //not using the shared dim becuase that was just for the shape
    rc = fc_segment(decayedSkin,0, &numDecayedSidesSegments,
		    &decayedSidesSegments);
    fc_deleteSubset(decayedSkin);
    if (rc!= FC_SUCCESS){
      fc_printfErrorMessage("Can't get Decayed sides segments. Skipping");
      continue;
    }

    if (!decayedSidesSegments){
      //tear doesnt intersect this shape
      continue;
    }

    //have to try to segment the whole innards to see if there is a break
    rc = fc_subsetSegmentsSubset(tear->region,currinnards,shareddim,
				 "live",&numlive,&livesegments);
    if (rc != FC_SUCCESS){
      fc_printfErrorMessage("Can't get segmentation on live shape");
      //just keep going
    }else{
      if (numlive > 1){
	//	printf("found a break - into %d segments of sizes: ",numlive);
	for (j = 0; j < numlive; j++){
	  //	  int num = -1;
	  //	  rc = fc_getSubsetNumMember(livesegments[j],&num);
	  //	  printf("%d ",num);
	  fc_deleteSubset(livesegments[j]);
	}
	//	printf("\n");
	free(livesegments);
	broken = 1;
      }
    }


    //tear does intersect this shape, set up for the intersection
    si = (ShapeIntersection*)malloc(sizeof(ShapeIntersection));
    if (!si){
      fc_exitIfError(FC_MEMORY_ERROR);
    }

    if (broken){
      si->type = BREAK;
    }else{
      if (numDecayedSidesSegments > 1){
	si->type = TUNNEL;
      }else {
	si->type = PIT;
      }
    }

    si->shapeID = wshapes[i].shapeID;
    si->meshID = wshapes[i].meshID;
    si->numIntersections = numDecayedSidesSegments;
    si->intersectedSides = (FC_SortedIntArray*)malloc((si->numIntersections)*
						      sizeof(FC_SortedIntArray));
    if (!(si->intersectedSides)){
      fc_exitIfError(FC_MEMORY_ERROR);
    }

    for (m = 0; m < si->numIntersections; m++){
      //each intersection of this tear with the skin
      int sidedecay;

      fc_initSortedIntArray(&(si->intersectedSides[m]));
      for (n = 0; n < currshape->numSides; n++){
	rc = fc_doSubsetsIntersect(decayedSidesSegments[m],
				   currshape->faces[n],
				   &sidedecay);
	fc_exitIfErrorPrintf(rc,"can't get subset intersection");
	
	if (sidedecay){
	  fc_addIntToSortedIntArray(&(si->intersectedSides[m]),n);
	}
      }
      
      //clean up
      fc_deleteSubset(decayedSidesSegments[m]);
    } //numintersections

    //now use intersections to determine major, subtype
    si->major = 0;
    for (m = 0; m < si->numIntersections; m++){
      if (order->nonadjorder){
	if (fc_isIntInSortedIntArray(&(si->intersectedSides[m]),
				     order->nonadjorder[0]) ||
	    fc_isIntInSortedIntArray(&(si->intersectedSides[m]),
				     order->nonadjorder[1])){
	  si->major = 1;
	  break;
	}
      }
      if (order->opposingorder){
	if (fc_isIntInSortedIntArray(&(si->intersectedSides[m]),
				     order->opposingorder[0]) ||
	    fc_isIntInSortedIntArray(&(si->intersectedSides[m]),
				     order->opposingorder[1])){
	  si->major = 1;
	  break;
	}
      }
    }

    //subtypes
    si->subtype = SINGLESIDE; //this will be updated below
    {
      FC_SortedIntArray uniqueintersections;
      fc_initSortedIntArray(&uniqueintersections);
      for (m = 0; m < si->numIntersections; m++){
	for (n = 0; n < si->intersectedSides[m].numVal; n++){
	  fc_addIntToSortedIntArray(&uniqueintersections,
				    si->intersectedSides[m].vals[n]);
	}
      }
      
      switch (uniqueintersections.numVal){
      case 0:
      case -1:
	fc_exitIfErrorPrintf(FC_ERROR,"Developer error: must be at least one intersection");
	break;
      case 1:
	si->subtype = SINGLESIDE;
	break;
      default:
	si->subtype = ADJSIDES;
	for (m = 0; m < uniqueintersections.numVal; m++){
	  for (n = m+1; n < uniqueintersections.numVal; n++){
	    if (currshape->adjmatrix[uniqueintersections.vals[m]][uniqueintersections.vals[n]]
		== 0){
	      si->subtype = NONADJSIDES;
	      break;
	    }
	  }
	  if (si->subtype == NONADJSIDES){
	    break;
	  }
	}
      }

      fc_freeSortedIntArray(&uniqueintersections);
    }


    //clean up
    free(decayedSidesSegments);

    //if weve gotten this far, we have an intersection
    break;

  } //numshapes

  tear->shape_intersection = si;
  return FC_SUCCESS;
}


// Expand the order array - make sure there is enough
static FC_ReturnCode grow_orders(int numNewShapes, int numCurrentShapes, 
				int* maxNumShapes, shape_order*** order) {
  shape_order** temp_order; //array of ptrs
  int temp_maxNumShapes;

  // internal, so no input testing

  // no new shapes do nothing
  if (numNewShapes < 1)
    return FC_SUCCESS;

  // if we are going, to run out of room, expand the array
  if (numNewShapes + numCurrentShapes > *maxNumShapes) {
    temp_maxNumShapes = 2*(numNewShapes + numCurrentShapes);
    temp_order = (shape_order**)realloc(*order,
					temp_maxNumShapes*
					sizeof(shape_order*));
    if (!temp_order) {
      fc_exitIfError(FC_MEMORY_ERROR);
    }
    *maxNumShapes = temp_maxNumShapes;
    *order = temp_order;
  }
  // else, do nothing

  // done
  return FC_SUCCESS;
}

int main(int argc, char** argv) {
  FC_ReturnCode rc;
  int i, j, k;
  char* dataset_file_name = NULL;
  char* output_file_name = NULL;
  char* displ_var_name = NULL;
  char** meshNames = NULL;
  FC_VerbosityLevel verbose_level = FC_QUIET;
  int numMesh = 0, numStep, temp_numStep;
  FC_Dataset dataset;
  FC_Mesh* meshes;
  FILE* output_file, *bb_file, *aux_file;
  int numDeathVar = 0;
  deathvarinfo* deathVarInfos;
  int numTear = 0, maxNumTear = 0;
  Tear *tears = NULL;
  int doDispl = 0, doCombine = 0;
  // This assumes elem_death and the displacements are seq vars 
  FC_Variable ***deathVars;  // deathVars[meshID][deathVarID][stepID]
  FC_Variable **displs;  // displs[meshID][stepID]
  int numSuperTear = 0;
  SuperTear* superTears = NULL;
  SuperTear superDuperTear; // grand total
  double min_dist;
  char name_buf[1028];
  int* topoDims;     // topoDim per mesh
  int globalTopoDim; // globalTopoDim = -1 means diff dims on diff meshes
  int doDamageWeight = 0;
  double dead_damage;
  int deathID = 0;  // assuming death criterion will be first
  int damageID = 1; // and damage criterion is second
  int stepID = -1;  // the id of the timestep of interest

  int useShape = 0; //default to not do shape calc
  int setangle = 0; //check for if the user set an angle 
  double shapeangle = SHAPE_ANGLE;
  int shareddim = SHARED_DIM;
  shape_order **shapeorder = NULL; //array of ptrs to already allocated shape_orders
  int numShapes = 0;
  int maxNumShapes = 0;

  // --- handle arguments

  if (argc < 2 ) {
  usage:
    printf("usage: %s [options] dataset [mesh_names]\n", 
           argv[0]);
    printf("options: \n");
    printf("   -h              : print this help message\n");
    printf("   -v              : verbose: print warning and error messages\n");
    printf("   -V              : very verbose: prints log and error messages\n");
    printf("   -o              : name of output file (default is stdout)\n");
    printf("                     value (default is \"<=\") \n");
    printf("   -d var_name     : calc the displaced values using this displacement\n");
    printf("                     variable\n");
    printf("   -m min_dist     : combine tears from different parts if they are less\n");
    printf("                     than (or equal to) a minimum distance apart\n");
    printf("   -n # <triplets> : Specify the rules for finding dead elements.\n");
    printf("                     The first parameter is the number of triplets. Each\n"); 
    printf("                     triplet is the name of variable to test, the comparison\n");
    printf("                     operator, and the value. The default is \n");
    printf("                     '-n 1 elem_death \"<=\" 0', meaning all elements\n");
    printf("                     where the variable 'elem_death' is less than or\n");
    printf("                     equal to zero are treated as being dead.\n");
    printf("   -t time_step    : The index of the timestep at which to find tears\n");
    printf("                     (indexing from 0). The default is to use the last one.\n");
    printf("   -w dead_damage  : Calculate damage weighted volumes. You must specify\n");
    printf("                     exactly 2 death criteria using the -n flag. The first\n");
    printf("                     must be for the dead elements and the second for the\n");
    printf("                     damaged elements. The supplied dead_damage value is\n");
    printf("                     used as the damage value for the dead elements.\n");
    printf("   -ss             : Include shape calculation, defaults to no shape calc \n");
    printf("                     shape (defaults to 40, must be between 0 and 180).\n");
    printf("   -sa angle       : Angle to be used for determining sides of a \n");
    printf("                     shape (defaults to 40, must be between 0 and 180).\n");
    printf("\n");
    printf("examples:\n");
    printf("   tears -d \"displ var\" data.ex2 - use default values but also calculate\n");
    printf("                       the displaced version.\n");
    printf("   tears -n 2 elem_death \"<=\" 0 damage \">\" .1 data.ex2 - tears are\n");
    printf("                       elements where elem_death is greater < 0 or damage\n");
    printf("                       is greater than .1\n");
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
    else if (!strcmp(argv[i], "-m")) {
      doCombine = 1;
      i++;
      min_dist = atof(argv[i]);
    }
    else if (!strcmp(argv[i], "-n")) {
      i++;
      rc = parseDeathVarSet(argc, argv, &i, &numDeathVar, &deathVarInfos);
      if (rc != FC_SUCCESS)
	goto usage;
    }
    else if (!strcmp(argv[i], "-t")) {
      i++;
      stepID = atoi(argv[i]);
    }
    else if (!strcmp(argv[i], "-w")) {
      doDamageWeight = 1;
      i++;
      dead_damage = atof(argv[i]);
    } else if (!strcmp(argv[i], "-ss")) {
      useShape = 1;
    } else if (!strcmp(argv[i], "-sa")) {
      char* end_ptr;
      i++;
      if (argc < i+1){
	goto usage;
      }
      shapeangle = strtod(argv[i], &end_ptr);
      if (*end_ptr == '\0'){ //its a number
        if (shapeangle < 0 || shapeangle > 180){
          fc_printfErrorMessage("Invalid angle: must be between 0-180");
          goto usage;
        }
      }else{
        fc_printfErrorMessage("Invalid angle: must be a number");
        goto usage;
      }
      setangle = 1;
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
  if (doCombine && min_dist < 0) 
    fc_exitIfErrorPrintf(FC_ERROR, "min_dist must be 0 or greater");
  if (doDamageWeight && numDeathVar != 2)
    fc_exitIfErrorPrintf(FC_ERROR, "you must specify 2 death criteria to "
			 "do damage weighting");
  if(setangle && !useShape)
    fc_exitIfErrorPrintf(FC_ERROR, "you must set use shape if you want to "
			 "specify an angle");


  // default death criteria
  if (numDeathVar < 1) {
    numDeathVar = 1;
    deathVarInfos = (deathvarinfo*)malloc(sizeof(deathvarinfo));
    if (deathVarInfos == NULL)
      fc_exitIfError(FC_MEMORY_ERROR);
    deathVarInfos[0].name = "elem_death";
    deathVarInfos[0].op = "<=";
    deathVarInfos[0].val = 0;
  }

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

      rc = fc_getMeshByName(dataset, meshNames[i], 
			    &numReturnMeshes,&returnMeshes);
      fc_exitIfErrorPrintf(rc, "Failed to get mesh by name");
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

  // Collect & test topo dims
  // 1) topodims need to be 2 or 3
  // 2) If combining, all topodims must be the same
  topoDims = (int*)malloc(numMesh*sizeof(int));
  if (!topoDims)
    fc_exitIfError(FC_MEMORY_ERROR);
  for (i = 0; i < numMesh; i++) {
    rc = fc_getMeshTopodim(meshes[i], &topoDims[i]);
    fc_exitIfErrorPrintf(rc, "Failed to get topo dims");
    if (topoDims[i] != 2 && topoDims[i] != 3)
      fc_exitIfErrorPrintf(FC_ERROR, "Meshes must be of topo dim 2 or 3");
    if (i == 0)
      globalTopoDim = topoDims[i];
    else if (topoDims[i] != topoDims[0]) {
      globalTopoDim = -1;
      if (doCombine)
	fc_exitIfErrorPrintf(FC_ERROR, "If combining tears across meshes, all"
			     "meshes must have the same topoDim");
    }
  }

  // Make sure required variable fields exist on the meshes
  // and that death vars are on elements
  deathVars = (FC_Variable***)malloc(numMesh*sizeof(FC_Variable**));
  if (!deathVars)
    fc_exitIfError(FC_MEMORY_ERROR);
  for (i = 0; i < numMesh; i++) {
    deathVars[i] = (FC_Variable**)malloc(numDeathVar*sizeof(FC_Variable*));
    if (!deathVars[i])
      fc_exitIfError(FC_MEMORY_ERROR);
    for (j = 0; j < numDeathVar; j++) {
      FC_AssociationType assoc;
      rc = fc_getOrGenerateUniqueSeqVariableByName(meshes[i], deathVarInfos[j].name,
						   &numStep, &deathVars[i][j]);
      fc_exitIfErrorPrintf(rc, "Failed to find element death variable '%s'",
			   deathVarInfos[j].name);
      if (!deathVars[i][j]){
	fc_exitIfErrorPrintf(FC_INPUT_ERROR,
			     "Failed to find (unique) element death variable '%s'",
			     deathVarInfos[j].name);
      }

      fc_getVariableAssociationType(deathVars[i][j][0], &assoc);
      if (assoc != FC_AT_ELEMENT) {
	fc_exitIfErrorPrintf(FC_ERROR, "Death var '%s' should be of assoc "
			     "%s not %s", deathVarInfos[j].name,
			     fc_getAssociationTypeText(FC_AT_ELEMENT),
			     fc_getAssociationTypeText(assoc));
      }
    }
  }
  if (doDispl) {
    displs = (FC_Variable**)malloc(numMesh*sizeof(FC_Variable*));
    if (!displs)
      fc_exitIfError(FC_MEMORY_ERROR);
    for (i = 0; i < numMesh; i++) {
      rc = fc_getOrGenerateUniqueSeqVariableByName(meshes[i], displ_var_name,
						   &temp_numStep, &displs[i]);

      fc_exitIfErrorPrintf(rc, "Failed to find element displacement variable '%s'",
			   displ_var_name);			   
      if (!displs[i]){
	fc_exitIfErrorPrintf(FC_INPUT_ERROR, 
			     "Failed to find (unique) displacement variable '%s'",
			     displ_var_name);
      }

      //this line was outside the loop but i think it was mean to be in
      if (temp_numStep != numStep)
	fc_exitIfErrorPrintf(FC_ERROR, "Element death var and displacment var"
			     "  must have the same number of steps");
    }
  }

  // set the stepID
  if (stepID == -1) {
    stepID = numStep-1;
  }
  else {
    if (stepID < 0 || stepID > numStep-1)
      fc_exitIfErrorPrintf(FC_ERROR, "stepID (%d) is out of range (0-%d)",
			   stepID, numStep-1);
  }

  // get the output file
  if (output_file_name) {
    output_file = fopen(output_file_name, "w");
    if (!output_file) 
      fc_exitIfErrorPrintf(FC_FILE_IO_ERROR, "Failed to open file '%s' to "
                           "write output to", output_file_name);
  }
  else
    output_file = stdout;

  // header
  fprintf(output_file, "Tear characterizations for dataset '%s'\n", 
          dataset_file_name);
  fprintf(output_file, "Tears criteria: '%s' %s %g", deathVarInfos[0].name,
	  deathVarInfos[0].op, deathVarInfos[0].val);
  for (i = 1; i < numDeathVar; i++) 
    fprintf(output_file, "  &&  '%s' %s %g", deathVarInfos[i].name,
	    deathVarInfos[i].op, deathVarInfos[i].val);
  fprintf(output_file, "\n");
  fprintf(output_file, "Time step index: %d\n", stepID);
  if (doDamageWeight) {
    fprintf(output_file, "Damage weighted volumes of dead elements calculated "
	    "using damage = %g\n", dead_damage);
  }
  fprintf(output_file, "%d mesh(es)\n", numMesh);
  fflush(NULL);
  
  // --- find (& characterize) dead element regions

  // for each mesh, look for tears in the specified time step
  for (i = 0; i < numMesh; i++) {
    int numSegment;
    FC_Subset threshPerDeathVar[numDeathVar], threshAll;
    FC_Subset *segsPerDeathVar[numDeathVar], *segments;

    shape_w_order *localwshapes = NULL;
    int numlocalshapes = 0;

    // Find dead regions in last step
    rc = fc_createSubset(meshes[i], "threshold", FC_AT_ELEMENT, &threshAll);
    fc_exitIfErrorPrintf(rc,"failed to create empty threshold subset");
    for (j = 0; j < numDeathVar; j++) {
      FC_Subset tempSubset;
      rc = fc_createThresholdSubset(deathVars[i][j][stepID],
				    deathVarInfos[j].op, deathVarInfos[j].val,
				    "temp", &threshPerDeathVar[j]);
      fc_exitIfErrorPrintf(rc, "Failed to theshold using op '%s' and value '%g'",
			   deathVarInfos[j].op, deathVarInfos[j].val);
      rc = fc_createSubsetIntersection(threshPerDeathVar[j], "OR", threshAll, 
				       "threshold", &tempSubset);
      fc_exitIfErrorPrintf(rc, "Failed to combine subsets");

      fc_deleteSubset(threshAll);
      threshAll = tempSubset;
    }
    rc = fc_segment(threshAll, 0, &numSegment, &segments);
    fc_deleteSubset(threshAll);
    fc_exitIfErrorPrintf(rc, "Failed to segment dead regions");

    // if doing damage weighting, split the segments by death var
    // Note: if elem passes both death criteria, dead claims elem since we 
    // are assuming damage variable has no meaning on dead elements
    if (doDamageWeight && numSegment > 0) {
      segsPerDeathVar[0] = (FC_Subset*)malloc(numSegment*sizeof(FC_Subset));
      segsPerDeathVar[1] = (FC_Subset*)malloc(numSegment*sizeof(FC_Subset));
      if (!segsPerDeathVar[0] || !segsPerDeathVar[1])
	fc_exitIfError(FC_MEMORY_ERROR);
      for (k = 0; k < numSegment; k++) {
	FC_Subset tempSubset, overlap;
	rc = fc_createSubsetIntersection(threshPerDeathVar[deathID], "AND", 
					 segments[k], "deadSegment", 
					 &segsPerDeathVar[deathID][k]);
	fc_exitIfErrorPrintf(rc, "Failed to AND subsets");
	rc = fc_createSubsetIntersection(threshPerDeathVar[damageID], "AND",
			       segments[k], "damageSegment", &tempSubset);
	fc_exitIfErrorPrintf(rc, "Failed to AND subsets");
	rc = fc_createSubsetIntersection(tempSubset, "AND",
					 segsPerDeathVar[deathID][k],
					 "temp", &overlap);
	fc_exitIfErrorPrintf(rc, "Failed to AND subsets");
	rc = fc_createSubsetIntersection(tempSubset, "XOR", overlap,
					 "damageSegment-deadSegment",
					 &segsPerDeathVar[damageID][k]);
	fc_exitIfErrorPrintf(rc, "Failed to AND subsets");
      }
    }
    for (j = 0; j < numDeathVar; j++) {
      fc_deleteSubset(threshPerDeathVar[j]);
    }

    // Midway report 
    fprintf(output_file, "Mesh %d: '%s' has %d dead element region(s)", 
	    i, meshNames[i], numSegment);
    fflush(output_file);

    //do shapes here and write out as part of the midway report
    if (useShape && numSegment > 0){
      //in order to check breaking, we need the innards of the shape and
      //therefore have to segment first
      int numinnards = 0;
      int numElem = 0;
      FC_Subset* innards = NULL;
      FC_Subset whole;

      rc = fc_getMeshNumElement(meshes[i],&numElem);
      fc_exitIfErrorPrintf(rc, "Failed get num elem");

      rc = fc_createSubset (meshes[i],"whole", FC_AT_ELEMENT,&whole);
      fc_exitIfErrorPrintf(rc, "Failed to create new subset whole");

      for (j = 0; j < numElem; j++){
	rc = fc_addMemberToSubset (whole,j);
	fc_exitIfErrorPrintf(rc, "Failed to add member to subset whole");
      }

      rc = fc_segment(whole,shareddim,&numinnards,&innards);
      fc_deleteSubset(whole);
      if (rc!= FC_SUCCESS){
	fc_printfErrorMessage("Can't get shapes for mesh %d. Skipping shape calculation.",i);
	numlocalshapes = 0;
      }else{
	//now make each innard into a shape
	numlocalshapes = numinnards;

	localwshapes = (shape_w_order*)malloc(numlocalshapes*sizeof(shape_w_order));
	if (!localwshapes){
	  fc_exitIfError(FC_MEMORY_ERROR);
	}

	for (j = 0; j < numinnards; j++){
	  //now make each innard into a separate shape w order
	  FC_Shape *singleshape;
	  int numsingleshape;

	  rc = fc_getSubsetShapes(innards[j],shapeangle,shareddim,&numsingleshape,&singleshape);
	  if (rc != FC_SUCCESS){
	    fc_printfErrorMessage("Can't get shapes for mesh %d. Skipping shape calculation.",i);
	    for (k = 0; k < numinnards; k++){
	      fc_deleteSubset(innards[k]);
	    }
	    numinnards = 0;
	    //	    free(innards); will free outside

	    //free any new shapes we already have for this mesh
	    for (k = 0; k < j; k++){
	      freeShapeWOrder(&localwshapes[k]);
	    }
	    free(localwshapes);
	    numlocalshapes = 0;
	  }else if (numsingleshape != 1){
	    fc_printfErrorMessage("Wrong shape division for mesh %d. Skipping shape calculation.",i);
	    for (k = 0; k < numinnards; k++){
	      fc_deleteSubset(innards[k]);
	    }
	    numinnards = 0;
	    //	    free(innards); will free outside
	    for (k = 0; k < j; k++){
	      freeShapeWOrder(&localwshapes[k]);
	    }
	    free(localwshapes);
	    numlocalshapes = 0;
	  }else{
	    createShapeWOrder(i,j,&(singleshape[0]),innards[j],&(localwshapes[j]));
	    //	    free(singleshape);
	    //never delete the shapes, later wshapes will free the memory
	  }
	} //each innard (shape)
	if (innards) free(innards);
      }
    }
    printf("\n");


    // Do as much characterization as possible of tears & save info
    if (numSegment > 0) {
      FC_Subset meshSkin;
      FC_Variable area_or_vol, logVol;
      FC_Variable volXdamage, logdamage, volXlogdamage, logVolXdamage;

      // expand the tears array
      rc = grow_tears(numSegment, numTear, &maxNumTear, &tears);
      fc_exitIfErrorPrintf(rc, "Failed to grow tears array");

      // get mesh skin
      rc = fc_getMeshSkin(meshes[i], &meshSkin);
      fc_exitIfErrorPrintf(rc, "Failed to get mesh skin");

      // get face areas or region volumes
      if (topoDims[i] == 2) {
	rc = fc_getFaceAreas(meshes[i], &area_or_vol);
	fc_exitIfErrorPrintf(rc, "Failed to get face areas");
      }
      else if (topoDims[i] == 3) {
	rc = fc_getRegionVolumes(meshes[i], &area_or_vol);
	fc_exitIfErrorPrintf(rc, "Failed to get region volumes");
      }

      // calc variables for damage weighting
      if (doDamageWeight) {
	FC_Variable damage = deathVars[i][damageID][stepID];
	rc = fc_varUnaryFunction(area_or_vol, log10, "logVol", &logVol);
	fc_exitIfErrorPrintf(rc, "Failed to tak log of var");
	rc = fc_varUnaryFunction(damage, log10, "logdamage", &logdamage);
	fc_exitIfErrorPrintf(rc, "Failed to take log of var");
	rc = fc_varOperatorVar(area_or_vol, "*", damage,
			       "volXdamage", &volXdamage);
	fc_exitIfErrorPrintf(rc, "Failed to multiply vars");
	rc = fc_varOperatorVar(area_or_vol, "*", logdamage, 
			       "volXlogdamage", &volXlogdamage);
	fc_exitIfErrorPrintf(rc, "Failed to multipl vars");
	rc = fc_varUnaryFunction(volXdamage, log10, "logVolXdamage",
				 &logVolXdamage);
	fc_exitIfErrorPrintf(rc, "Failed to take log of var");
	fc_deleteVariable(logdamage);
      }

      // Process each tear
      // Loop over segments - each is a new tear
      for (j = 0; j < numSegment; j++) {
	Tear *tear = &tears[numTear];
	Tear_Characterization *tearData = &tears[numTear].data;

	int numExposedMember;

	// basic info
	tear->tearID = numTear;
	tear->region = segments[j];
	tear->meshID = i;
	tear->stepID = stepID;
	
	// get exposed surface - the "live" part of tear
	rc = fc_getExposedSkin(segments[j], &meshSkin, &tear->exposed);
	fc_exitIfErrorPrintf(rc, "Failed to get exposed skin");
	
	fc_getSubsetNumMember(tear->exposed, &numExposedMember);
	if (numExposedMember < 1) {
	  fc_printfWarningMessage("There are no members in an exposed skin"
				  "subset; some exposed measurements are NA");
	}

	// rename the subsets for friendly output
	sprintf(name_buf, "region%d", numTear);
	fc_changeSubsetName(tear->region, name_buf);
	sprintf(name_buf, "exposed%d", numTear);
	fc_changeSubsetName(tear->exposed, name_buf);

        // Number of cells in region
        rc = fc_getSubsetNumMember(segments[j], &tearData->numCell);
        fc_exitIfErrorPrintf(rc, "Failed to get number of members in region");
        
        // area/volume of region - nondisplaced only
	rc = fc_getVariableSubsetSum(area_or_vol, segments[j], 
				     &tearData->volume);
	fc_exitIfErrorPrintf(rc, "Failed to get area of region");

	// bounding boxes
	// bounding box of region - nondisplaced only
	rc = fc_getSubsetBoundingBox(segments[j], &tearData->numDim, 
				     &tearData->lowers, &tearData->uppers);
	fc_exitIfErrorPrintf(rc, "Failed to get region bounding box");
	// bounding boxes of exposed 
	if (numExposedMember > 0) {
	  rc = fc_getSubsetBoundingBox(tear->exposed, &tearData->numDim, 
				       &tearData->exp_lowers, 
				       &tearData->exp_uppers);
	  fc_exitIfErrorPrintf(rc, "Failed to get exposed bounding box");
	  if (doDispl) {
	    rc = fc_getDisplacedSubsetBoundingBox(tear->exposed,
			   displs[i][stepID],  &tearData->numDim, 
			   &tearData->displ_exp_lowers, 
                           &tearData->displ_exp_uppers);
	    fc_exitIfErrorPrintf(rc, "Failed to get displaced exposed bb");
	  }
	}
	else {
	  for (k = 0; k < 3; k++) {
	    tearData->exp_lowers[k] = -2;
	    tearData->exp_uppers[k] = -1;
	    tearData->displ_exp_lowers[k] = -2;
	    tearData->displ_exp_uppers[k] = -1;
	  }
	}
        
	// diameters
        // diameter of region - nondisplaced only
        rc = fc_getSubsetDiameter(segments[j], &tearData->diameter, NULL, 
				  NULL);
        fc_exitIfErrorPrintf(rc, "Failed to get region diameter");
	// diameters of exposed
	if (numExposedMember > 0) {
	  rc = fc_getSubsetDiameter(tear->exposed, &tearData->exp_diameter,
				    NULL, NULL);
	  fc_exitIfErrorPrintf(rc, "Failed to get exposed diameter");
	  if (doDispl) {
	    rc = fc_getDisplacedSubsetDiameter(tear->exposed, 
					       displs[i][stepID], 
					       &tearData->displ_exp_diameter,
					       NULL, NULL);
	    fc_exitIfErrorPrintf(rc, "Failed to get displaced exposed diameter");
	  }
	}
	else {
	  tearData->exp_diameter = -1;
	  tearData->displ_exp_diameter = -1;
	}
	
	// damage weighted volumes
	if (doDamageWeight) {
	  FC_Subset deadSubset = segsPerDeathVar[deathID][j];
	  FC_Subset damageSubset = segsPerDeathVar[damageID][j];
	  fc_getSubsetNumMember(deadSubset, &tearData->dead_numCell);
	  fc_getSubsetNumMember(damageSubset, &tearData->damaged_numCell);
	  if (tearData->dead_numCell < 1) {
	    tearData->dead_vol = 0;
	    tearData->dead_logvol = 0;
	  }
	  else {
	    rc = fc_getVariableSubsetSum(area_or_vol, deadSubset,
					 &tearData->dead_vol); 
	    fc_exitIfErrorPrintf(rc, "Failed to get volume of dead part of tear");
	    rc = fc_getVariableSubsetSum(logVol, deadSubset,
					 &tearData->dead_logvol);
	    fc_exitIfErrorPrintf(rc, "Failed to get sum of log dead vol");
	  }
	  if (tearData->damaged_numCell < 1) {
	    tearData->damaged_vol = 0;
	    tearData->damaged_vol_x_damage = 0;
	    tearData->damaged_vol_x_logdamage = 0;
	    tearData->damaged_log_vol_x_damage = 0;
	  }
	  else {
	    rc = fc_getVariableSubsetSum(area_or_vol, damageSubset,
					 &tearData->damaged_vol);
	    fc_exitIfErrorPrintf(rc, "Failed to get volume of damage");
	    rc = fc_getVariableSubsetSum(volXdamage, damageSubset,
					 &tearData->damaged_vol_x_damage); 
	    fc_exitIfErrorPrintf(rc, "Failed to get sum of vol x damage");
	    rc = fc_getVariableSubsetSum(volXlogdamage, damageSubset,
					 &tearData->damaged_vol_x_logdamage);
	    fc_exitIfErrorPrintf(rc, "Failed to get sum of vol x logdamage");
	    rc = fc_getVariableSubsetSum(logVolXdamage,damageSubset,
					 &tearData->damaged_log_vol_x_damage);
	    fc_exitIfErrorPrintf(rc, "Failed to get sum of log (vol x damage)");
	  }
	}

	//will return null shape_intersection to the tear if no shapes
	rc = calcShapeIntersection(numlocalshapes,localwshapes,shareddim,tear);

	// done with this tear
	numTear++;
      }

      //add only the orders to the order array,
      //get rid of the shapes. 
      rc = grow_orders(numlocalshapes, numShapes,
		       &maxNumShapes, &shapeorder);

      for (j = 0; j < numlocalshapes; j++){
	shapeorder[numShapes+j] = localwshapes[j].order;
	freeShapeWOrder_ShapeOnly(&(localwshapes[j])); //frees the localshape
	localwshapes[j].order = NULL;
      }
      numShapes+=numlocalshapes;
      if (localwshapes)  free(localwshapes); //free the array

      // cleanup
      fc_deleteSubset(meshSkin);
      fc_deleteVariable(area_or_vol);
      free(segments);
      if (doDamageWeight) {
	free(segsPerDeathVar[deathID]);
	free(segsPerDeathVar[damageID]);
	fc_deleteVariable(logVol);
	fc_deleteVariable(volXdamage);
	fc_deleteVariable(volXlogdamage);
	fc_deleteVariable(logVolXdamage);
      }
    } //if numSegment > 0

    // get rid of extra data structures on meshes, maybe delete!
    if (numSegment == 0)
      fc_deleteMesh(meshes[i]);
    else
      fc_releaseMesh(meshes[i]);
  }
  for (i = 0; i < numMesh; i++) {
    for (j = 0; j < numDeathVar; j++)
      free(deathVars[i][j]);
    free(deathVars[i]);
  }
  free(deathVars);
  free(meshes);

  // Report how many tears found
  // Early exit: if no tears, we're done
  if (numTear == 0) {
    fprintf(output_file, "No tears found\n");
    if (output_file_name)
      fclose(output_file);
    fc_finalLibrary();
    exit(0);
  }
  else {
    fprintf(output_file, "Found %d dead element regions\n", numTear);
    fflush(NULL);
  }

  // --- Combine tears (if requested)

  // NOTE: whether combining or not, will convert from Tear to SuperTear type
  if (!doCombine) {
    // report that we are combining
    fprintf(output_file, "Combining of dead elem regions not requested\n");
    fflush(NULL);

    numSuperTear = numTear;
    superTears = (SuperTear*)malloc(numTear*sizeof(SuperTear));
    for (i = 0; i < numTear; i++) {
      superTears[i].tearID = tears[i].tearID;
      superTears[i].numTear = 1;
      superTears[i].tears = (Tear**)malloc(sizeof(Tear*));
      superTears[i].tears[0] = &tears[i];
    }
  }
  else {
    // NOTE: this is done with time zero geometry
    int** proxMatrix;
    int* numPerGroup, **idsPerGroup, numGroup, count;

    // report that we are combining
    fprintf(output_file, "Combining dead elem regions that have verts "
	    "within '%g' of each other ...\n", min_dist);
    fflush(NULL);
	    
    proxMatrix = (int**)malloc(numTear*sizeof(int*));
    for (i = 0; i < numTear; i++) 
      proxMatrix[i] = (int*)malloc(numTear*sizeof(int));
    // get the numPairs (proxities) - store in matrix 
    // FIX? is next thing what we really want to do?
    // Don't test if on same mesh, that was already done with segment!
    for (i = 0; i < numTear; i++) {
      for (j = i+1; j < numTear; j++) {
	if (tears[i].meshID == tears[j].meshID)
	  proxMatrix[i][j] = 0;
	else {
	  rc = fc_getSubsetsProximity(tears[i].region, tears[j].region, 
				      min_dist, &proxMatrix[i][j], NULL, NULL);
	  fc_exitIfErrorPrintf(rc, "failed to get subsets' proximity");
	}
      }
    }

    // debug - look at matrix
    /* 
    printf("ProxMatrix\n");
    printf("     ");
    for (i = 0; i < numTear; i++)
      printf("%3d ", i);
    printf("\n");
    printf("-----");
    for (i = 0; i < numTear; i++)
      printf("----");
    printf("\n");
    for (i = 0; i < numTear; i++) {
      printf("%3d: ", i);
      for (j = 0; j < numTear; j++) {
	if (j <= i)
	  printf("%3s ", "-");
	else
	  printf("%3d ", proxMatrix[i][j]);
      }
      printf("\n");
    }
    */

    // determine groups
    rc = fc_groupSelfCorrespond(numTear, proxMatrix, &numGroup,
				&numPerGroup, &idsPerGroup);
    fc_exitIfErrorPrintf(rc, "failed to group by proximity");
    for (i = 0; i < numTear; i++)
      free(proxMatrix[i]);
    free(proxMatrix);

    // make the super tears 
    superTears = malloc(numGroup*sizeof(SuperTear));
    numSuperTear = 0;
    count = 0; // for a quick check
    for (i = 0; i < numGroup; i++) {
      superTears[numSuperTear].tearID = numSuperTear;
      superTears[numSuperTear].numTear = numPerGroup[i];
      count += numPerGroup[i];
      superTears[numSuperTear].tears = malloc(numPerGroup[i]*sizeof(Tear*));
      for (j = 0; j < numPerGroup[i]; j++)
	superTears[numSuperTear].tears[j] = &tears[idsPerGroup[i][j]];
      numSuperTear++;
      free(idsPerGroup[i]);
    }
    free(numPerGroup);
    free(idsPerGroup);
    if (count != numTear) 
      fprintf(output_file, 
	      "WARNING! Algorithm error. Please report to FCLib developers.\n"
	      "WARNING! Total number of tears in groups (%d) does not match\n"
	      "WARNING! original number of tears(%d)\n", count, numTear);
  }

  // report number of super tear tears
  fprintf(output_file, "Found %d tears\n", numSuperTear);
  fflush(NULL);

  // --- Create data for the super tears

  // create data for the super tears
  for (i = 0; i < numSuperTear; i++) {
    int numDim;
    FC_Coords temp_lowers, temp_uppers;
    Tear_Characterization *superTearData;

    // copy 1st tear's data -- if numTear = 1, we are done
    superTears[i].data = superTears[i].tears[0]->data;
    if (superTears[i].numTear == 1)
      continue;

    // shortcuts for following work
    numDim = superTears[i].data.numDim;
    superTearData = &superTears[i].data;

    // for some quantities, can just add other tears' data
    for (j = 1; j < superTears[i].numTear; j++) {
      Tear_Characterization *tearData = &superTears[i].tears[j]->data;
      superTearData->numCell += tearData->numCell;
      superTearData->volume += tearData->volume;
      for (k = 0; k < numDim; k++) {
	temp_lowers[k] = superTearData->lowers[k];
	temp_uppers[k] = superTearData->uppers[k];
      }
      fc_combineBoundingBoxes(numDim, temp_lowers, temp_uppers,
			      tearData->lowers, tearData->uppers,
			      &superTearData->lowers,
			      &superTearData->uppers);
      for (k = 0; k < numDim; k++) {
	temp_lowers[k] = superTearData->exp_lowers[k];
	temp_uppers[k] = superTearData->exp_uppers[k];
      }
      fc_combineBoundingBoxes(numDim, temp_lowers, temp_uppers,
			      tearData->exp_lowers, tearData->exp_uppers,
			      &superTearData->exp_lowers,
			      &superTearData->exp_uppers);
      if (doDispl) {
	for (k = 0; k < numDim; k++) {
	  temp_lowers[k] = superTearData->displ_exp_lowers[k];
	  temp_uppers[k] = superTearData->displ_exp_uppers[k];
	}
	fc_combineBoundingBoxes(numDim, temp_lowers, temp_uppers,
				tearData->displ_exp_lowers,
				tearData->displ_exp_uppers,
				&superTearData->displ_exp_lowers,
				&superTearData->displ_exp_uppers);
      }
      if (doDamageWeight) {
	superTearData->dead_numCell += tearData->dead_numCell;
	superTearData->damaged_numCell += tearData->damaged_numCell;
	superTearData->dead_vol += tearData->dead_vol;
	superTearData->dead_logvol += tearData->dead_logvol;
	superTearData->damaged_vol += tearData->damaged_vol;
	superTearData->damaged_vol_x_damage += tearData->damaged_vol_x_damage;
	superTearData->damaged_vol_x_logdamage += tearData->damaged_vol_x_logdamage;
	superTearData->damaged_log_vol_x_damage += tearData->damaged_log_vol_x_damage;
      }
    }
    // for some quantities, have to recalculate using all tears
    if (superTears[i].numTear > 1) {
      int exposedsExist = 1;
      FC_Subset tear_regions[superTears[i].numTear];
      FC_Subset tear_exposeds[superTears[i].numTear];
      FC_Variable tear_displs[superTears[i].numTear];
      for (j = 0; j < superTears[i].numTear; j++) {
	if (FC_HANDLE_EQUIV(superTears[i].tears[j]->exposed, FC_NULL_SUBSET)) {
	  exposedsExist = 0;
	  break;
	}
      }
      for (j = 0; j < superTears[i].numTear; j++) {
	Tear* tear = superTears[i].tears[j];
	tear_regions[j] = tear->region;
	tear_exposeds[j] = tear->exposed;
	if (doDispl)
	  tear_displs[j] = displs[tear->meshID][tear->stepID];
      }
      rc = fc_getSubsetsDiameter(superTears[i].numTear, tear_regions,
				 &superTearData->diameter, NULL, NULL, NULL);
      if (exposedsExist) {
	rc = fc_getSubsetsDiameter(superTears[i].numTear, tear_exposeds,
				   &superTearData->exp_diameter, NULL,
				   NULL, NULL);
	if (doDispl) {
	  rc = fc_getDisplacedSubsetsDiameter(superTears[i].numTear,
					      tear_exposeds, tear_displs,
					  &superTearData->displ_exp_diameter,
					      NULL, NULL, NULL);
	  fc_exitIfErrorPrintf(rc, "Failed to get displaced exposed diameter");
	}
      }
      else {
	for (k = 0; k < 3; k++) {
	  superTearData->exp_lowers[k] = -2;
	  superTearData->exp_uppers[k] = -1;
	  superTearData->displ_exp_lowers[k] = -2;
	  superTearData->displ_exp_uppers[k] = -1;
	}
      }
    }
  }

  // midway cleanup - don't need displs any more
  if (doDispl) {
    for (i = 0; i < numMesh; i++)
      free(displs[i]);
    free(displs);
  }

  // --- sort

  // sort the super tears
  fprintf(output_file, "Sorting tears by region diameter (largest first) ...\n");
  fflush(NULL);
  qsort((void*)superTears, (size_t)numSuperTear, sizeof(SuperTear), 
	cmpSuperTearsByLength);
  //for (i = 0; i < numSuperTear; i++) 
  //printf("%d: superTear %d, diameter = %f, bbchord = %f\n", i,
  //   superTears[i].tearID, superTears[i].data.diameter,
  //   sqrt(pow(superTears[i].data.uppers[0]-superTears[i].data.lowers[0],2)+
  //	pow(superTears[i].data.uppers[1]-superTears[i].data.lowers[1],2)+
  //	pow(superTears[i].data.uppers[2]-superTears[i].data.lowers[2],2)));

  // --- calc totals of all super tears - the super duper tear

  // only some values are valid to super sum
  // and some value have slightly different meaning as total
  superDuperTear.tearID = -1;
  superDuperTear.numTear = superTears[0].numTear;
  superDuperTear.tears = NULL;
  superDuperTear.data = superTears[0].data;
  for (i = 1; i < numSuperTear; i++) {
    Tear_Characterization *superTearData = &superTears[i].data;
    superDuperTear.numTear += superTears[i].numTear;
    superDuperTear.data.numCell += superTearData->numCell;
    superDuperTear.data.volume += superTearData->volume;
    superDuperTear.data.diameter += superTearData->diameter;
    superDuperTear.data.exp_diameter += superTearData->exp_diameter;
    superDuperTear.data.displ_exp_diameter += superTearData->displ_exp_diameter;
    // bounding boxes don't mean anything
    if (doDamageWeight) {
      superDuperTear.data.dead_numCell += superTearData->dead_numCell;
      superDuperTear.data.damaged_numCell += superTearData->damaged_numCell;
      superDuperTear.data.dead_vol += superTearData->dead_vol;
      superDuperTear.data.dead_logvol += superTearData->dead_logvol;
      superDuperTear.data.damaged_vol += superTearData->damaged_vol;
      superDuperTear.data.damaged_vol_x_damage += superTearData->damaged_vol_x_damage;
      superDuperTear.data.damaged_vol_x_logdamage += superTearData->damaged_vol_x_logdamage;
      superDuperTear.data.damaged_log_vol_x_damage += superTearData->damaged_log_vol_x_damage;
    }
  }
  
  // --- report tears

  // print results per super tear
  for (i = 0; i < numSuperTear; i++) {
    SuperTear *superTear = &superTears[i];
    fprintf(output_file, "Tear %d:\n", i);
    fprintf(output_file, "  numDeadElementRegions = %d\n", superTear->numTear);
    // debug
    //fprintf(output_file, "  tearIDs = %d", superTear->tears[0]->tearID);
    //for (j = 1; j < superTear->numTear; j++)
    //fprintf(output_file, ", %d", superTear->tears[j]->tearID);
    //fprintf(output_file, "\n");
    fprintf(output_file, "  meshIDs = %d", superTear->tears[0]->meshID);
    for (j = 1; j < superTear->numTear; j++)
      fprintf(output_file, ", %d", superTear->tears[j]->meshID);
    fprintf(output_file, "\n");
    fprintf(output_file, "  meshNames = '%s'", 
	    meshNames[superTear->tears[0]->meshID]);
    for (j = 1; j < superTear->numTear; j++)
      fprintf(output_file, ", '%s'", meshNames[superTear->tears[j]->meshID]);
    fprintf(output_file, "\n");
    for (j = 0; j < superTear->numTear; j++){
      printShapeIntersection(output_file,superTear->tears[j]->tearID,
			     superTear->tears[j]->shape_intersection);
    }
    print_tear_characterization(output_file, &superTear->data, 0,
				topoDims[superTear->tears[0]->meshID],
				doDispl, doDamageWeight, dead_damage);
    fflush(NULL);
  }

  // --- report totals

  // print totals
  fprintf(output_file, "Totals:\n");
  fprintf(output_file, "  numDeadElementRegions = %d\n",
	  superDuperTear.numTear);
  print_tear_characterization(output_file, &superDuperTear.data, 1,
			      globalTopoDim, doDispl, doDamageWeight, 
			      dead_damage);
  fflush(NULL);


  // -- shape centric writeout
  if (useShape){
    int tot = 0;
    int totpits = 0;
    int totbreaks = 0;
    int tottunnels = 0;

    fprintf(output_file,"\nShape intersections:\n");
    for (i = 0; i < numShapes; i++){
      int pits = 0;
      int tunnels = 0;
      int breaks = 0;
      int pitsm = 0;
      int tunnelsm = 0;
      int breaksm = 0;
      int total = 0;
      printShapeOrder(output_file,shapeorder[i]);
      for (j = 0; j < numSuperTear; j++){
	SuperTear *superTear = &superTears[j];
	for (k = 0; k < superTear->numTear; k++){
	  ShapeIntersection *si = superTear->tears[k]->shape_intersection;
	  if (si && (shapeorder[i]->meshID == si ->meshID) &&
	      (shapeorder[i]->shapeID == si->shapeID)){
	    fprintf(output_file,"  SuperTear %4d ", superTear->tearID);
	    printShapeIntersection(output_file,superTear->tears[k]->tearID, si);
	    total++;
	    switch(si->type){
	    case (BREAK):
	      breaks++;
	      if (si->major) breaksm++;
	      break;
	    case (PIT):
	      pits++;
	      if (si->major) pitsm++;
	      break;
	    case (TUNNEL):
	      tunnels++;
	      if (si->major) tunnelsm++;
	      break;
	    }
	  }
	}
      } //num super tear
      fprintf(output_file,"  Shape Tears (%d): BREAKS (%d/%d) TUNNELS (%d/%d) PITS (%d/%d)\n",
	      total,breaksm,breaks,tunnelsm,tunnels,pitsm,pits);
      tot+=total;
      totpits+=pits;
      totbreaks+=breaks;
      tottunnels+=tunnels;
    } //num shapes
    fprintf(output_file,"Total Shape Tears %d: BREAKS (%d) TUNNELS (%d) PITS (%d)\n",
	      tot,totbreaks,tottunnels,totpits);

  }
      

  // --- print aux file
  sprintf(name_buf, "tears.fcx");
  //printf("Creating aux file '%s':\n", name_buf);
  aux_file = fopen(name_buf, "w");
  if (!aux_file)
    fc_exitIfErrorPrintf(FC_FILE_IO_ERROR, "Failed to open %s for writing",
			 name_buf);
  _fc_writeAuxFileHeader(aux_file);
  for (i = 0; i < numSuperTear; i++) {
    sprintf(name_buf, "Tear%d", i);
    char** temp_subsetNames;
    temp_subsetNames = (char**)malloc(superTears[i].numTear*sizeof(char*));
    for (j = 0; j < superTears[i].numTear; j++) 
      fc_getSubsetName(superTears[i].tears[j]->region, &temp_subsetNames[j]);
    _fc_writeAuxFileTear(aux_file, name_buf, superTears[i].numTear,
			 temp_subsetNames, superTears[i].data.diameter);
    for (j = 0; j < superTears[i].numTear; j++)
      free(temp_subsetNames[j]);
    free(temp_subsetNames);
  }
  for (i = 0; i < numTear; i++) {
    _fc_writeAuxFileSubset(aux_file, tears[i].region);
  }
  _fc_writeAuxFileFooter(aux_file);

  // --- print bounding boxes

  // DEBUG
  /*
  sprintf(name_buf, "tear-regions-origs.bb");
  printf("Creating bounding box file '%s':\n", name_buf);
  bb_file = fopen(name_buf, "w");
  if (!bb_file)
    fc_exitIfErrorPrintf(FC_FILE_IO_ERROR, "Failed to open %s for writing",
                         name_buf);
  for (i = 0; i < numTear; i++) {
    print_bb_line(bb_file, superTears[i].tears[0]->stepID,
		  tears[i].data.numDim,
		  tears[i].data.lowers, tears[i].data.uppers);
  }
  fclose(bb_file);
  */

  sprintf(name_buf, "tear-regions.bb");
  printf("Creating bounding box file '%s':\n", name_buf);
  bb_file = fopen(name_buf, "w");
  if (!bb_file)
    fc_exitIfErrorPrintf(FC_FILE_IO_ERROR, "Failed to open %s for writing", 
			 name_buf);
  _fc_writeBBFileHeader(bb_file);
  for (i = 0; i < numSuperTear; i++) {
    char comment_buf[1028];
    sprintf(comment_buf, "Diameter=%g", superTears[i].data.diameter);
    sprintf(name_buf, "Tear%d", i);
    _fc_writeBBFileBoundingBox(bb_file, name_buf, superTears[i].tears[0]->stepID,
		  comment_buf, superTears[i].data.numDim,
		  superTears[i].data.lowers, superTears[i].data.uppers);
  }
  _fc_writeBBFileFooter(bb_file);
  fclose(bb_file);

  sprintf(name_buf, "tear-exposeds.bb");
  printf("Creating bounding box file '%s':\n", name_buf);
  bb_file = fopen(name_buf, "w");
  if (!bb_file)
    fc_exitIfErrorPrintf(FC_FILE_IO_ERROR, "Failed to open %sfor writing",
			 name_buf);
  _fc_writeBBFileHeader(bb_file);
  for (i = 0; i < numSuperTear; i++) {
    char comment_buf[1028];
    sprintf(name_buf, "Tear%d", i);
    sprintf(comment_buf, "Diameter=%g", superTears[i].data.diameter);
    _fc_writeBBFileBoundingBox(bb_file, name_buf, superTears[i].tears[0]->stepID,
		  comment_buf, superTears[i].data.numDim,
		  superTears[i].data.exp_lowers, 
		  superTears[i].data.exp_uppers);
  }
  _fc_writeBBFileFooter(bb_file);
  fclose(bb_file);

  if (doDispl) {
    sprintf(name_buf, "displ-tear-exposeds.bb");
    printf("Creating displaced bounding box file '%s':\n", name_buf);
    bb_file = fopen(name_buf, "w");
    if (!bb_file)
      fc_exitIfErrorPrintf(FC_FILE_IO_ERROR, "Failed to open %s for writing",
			   name_buf);
    _fc_writeBBFileHeader(bb_file);
    for (i = 0; i < numSuperTear; i++) {
      char comment_buf[1028];
      sprintf(name_buf, "Tear%d", i);
      sprintf(comment_buf, "Diameter=%g DisplDiameter=%g", 
	      superTears[i].data.diameter,
	      superTears[i].data.displ_exp_diameter);
      _fc_writeBBFileBoundingBox(bb_file, name_buf, superTears[i].tears[0]->stepID,
		       comment_buf, superTears[i].data.numDim,
		       superTears[i].data.displ_exp_lowers,
		       superTears[i].data.displ_exp_uppers);
    }
    _fc_writeBBFileFooter(bb_file);
    fclose(bb_file);
  }

  // --- cleanup

  //clean up shape orders
  for (j = 0; j < numShapes; j++){
    freeShapeOrder(shapeorder[j]);
  }
  free(shapeorder); //free the array

  // shut down output file
  if (output_file_name)
    fclose(output_file);

  // cleanup memory
  for (i = 0; i < numMesh; i++)
    free(meshNames[i]);
  free(meshNames);
  for (i = 0; i < numSuperTear; i++){
    for (j = 0; j < superTears[i].numTear; j++){
      freeShapeIntersection(superTears[i].tears[j]->shape_intersection);
      superTears[i].tears[j]->shape_intersection = NULL;
    }
    free(superTears[i].tears);
  }
  free(superTears);
  free(tears);
  fc_deleteDataset(dataset);
  fc_finalLibrary();

  exit(0);

}

