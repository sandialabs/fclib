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
 * \file screwBreaks.c
 * \brief screw breaking analysis
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/tools/screwBreaks.c,v $
 * $Revision: 1.26 $
 * $Date: 2006/10/16 20:04:39 $
 *
 * \description
 *
 *    DEAD VARS:
 *    User specifies elements to be considered dead
 *    by a number of variables and threshold values
 *    for each specified on the command line. The dead
 *    set is thus the union of all the specified sets.
 *
 *    BREAKAGE CONDITIONS:
 *    Tool reports when the screws first break or when
 *    any screw side entirely erodes due to the dead set. 
 *    The eroded case means the screw has separated
 *    from whatever that side was attached to (this
 *    is typically a tie surface, which is therefore
 *    a special case breakage condition, where the
 *    screw has not split into pieces.). Also gives
 *    total number of broken screws.
 *
 *    BREAKAGE RATIO:
 *    Tool also reports the breakage ratio. This is a rough
 *    estimate of the accumulated cross sectional area of
 *    the dead elements in the column above the screw base.
 *    The screw base is the smallest side which is
 *    connected to only one other side - in a normal shaped
 *    screw this is the bottom circular section. In a 
 *    regular cylinder this will be either one of the
 *    two circular sections. This calculation considers
 *    the non-displaced areas and locations of the elements.
 *
 *    Explicitly, the breakage ratio is determined by
 *    1) projecting the dead elements in the screw and
 *    the base surface elements both onto a plane with
 *    normal = normal of the screw base. 2) determining which
 *    elements in the screw base are overlapped by the projected
 *    elements - this is done by comparing the bounding
 *    boxes of the projected dead elements and the screw
 *    base elements and is therefore an overestimation of
 *    overlap 3) summing the area of the overlapped screw
 *    base faces 4) calculating the ratio relative to
 *    the complete screw base area. Be careful how
 *    you use this measure - it is an overestimation and
 *    does not take into consideration how far the
 *    dead elements are relative to one another. 
 *
 *    RESULT OUTPUT:
 *    Results are presented in the following fashion:
 *    Every output line lists the Mesh, ScrewID and the timestep #,
 *    to make it easy to grep through an output file. Screw breakage
 *    ratios are reported for all time steps for each screw, up
 *    until the screw is broken (or is eroded).
 *
 *    SCREWS:
 *    The screwmeshes must consist entirely of distinct screws,
 *    but there can be multiple screws in each screwmesh.
 *    Screws are defined by segmenting the screwmesh -- there
 *    is nothing to check that there is anything screw-like
 *    about the screws, other than there must be at least one
 *    distinct "end" which is used in the calculation of the
 *    breakage ratio.
 *
 *    Screws are identified by an arbitary identifier
 *    which will be the same from run to run, as long
 *    as the mesh is identical - renumberings of the
 *    mesh verticies, etc, even if it is topologically
 *    the same may result in a different screw identifier.
 *
 *    Screwmeshes should be input by the same name they have
 *    if you run the fcdump tool on the dataset.(This
 *    name may be different than that given it by different
 *    tools - i.e., may be a different name than appears
 *    in ensight).
 *
 *    OTHER:
 *    Seg_dim is the segmenation_dimension for determining
 *    breakage, defaults to 0. The segmentation dimension 
 *    is used to specify the conditions for which live regions
 *    in each screw are considered separate. If 2 regions are
 *    to be considered separate if no element in one region
 *    shares a face with any element in another region (but they
 *    may share edges and verticies), then seg_dim =2. If the
 *    elements in the 2 regions cannot share faces or edges, but
 *    may share verticies, choose seg_dim = 1. If the elements
 *    in the two regions cannot share even verticies, 
 *    then seg_dim = 0.
 *
 * \todo
 *    - Note: if we decide to drop the breakage ratio, then this
 *      would be generally useful, since the constraint on the
 *      shapes having an identifiable base is only necessary for 
 *      calculating the BR. We could potentially drop that
 *      since we now have the damage variable as something we
 *      can loop through to get closeness to breaking.
 *    - also print out screw verticies or bounding box, so can
 *      ID screws from run to run.
 *
 * \modifications
 *   - 03/09/2006 ACG created. Modification of screw.c.
 *   - 03/14/2006 ACG changed flag for seg_dim, removed
 *                    writeout of nonbroken screws.   
 *   - 03/14/2006 ACG added breakage ratio
 *   - 03/20/2006 ACG commented out local normal function
 *                    since now have fc_getPlanarSideNormal
 *   - 03/22/2006 ACG user can specify more than 1 variable
 *                    and cutoff for the dead set.
 *   - 03/22/2006 ACG removed commented out screwbasenormal function
 *   - 03/24/2006 ACG switching to more general shape screws - 
 *                    using getShapeEnds rather than getScrewSides
 *   - 03/24/2006 ACG struct for screw and deathvars for clarity
 *   - 03/24/2006 ACG put proj function into its own subroutine.
 *                    now projects to
 *                    smallest end, no matter what the shape.
 *   - 03/26/2006 ACG moved some stuff into subroutines for clarity. Trying
 *                    to do that with screwBuild but ive got a memory leak...
 *   - 03/26/2006 ACG added loop to check all sides erosion, untested. Moved
 *                    broken, eroded into their own subroutines for symmetry
 *   - 03/27/2006 ACG added printscrew - note IDS start at zero per mesh and
 *                    therefore match ensight format id, but not the lsdyna ids
 *                    (ie, cant match up element ids looking at d3plot file)
 *   - 03/27/2006 ACG WSK fixed memory leak
 *   - 03/28/2006 ACG verified projection results, which resulted in bug fix
 *                    on geom
 *   - 03/28/2006 ACG moved to tools from sandbox (was screw_tool.c)
 *   - 03/29/2006 ACG documentation changes
 *   - 03/29/2006 ACG removed option to specify screwmesh by number becuase
 *                    that might get confused with the material ID. now must
 *                    specify mesh bysame name it has when you do fcdump.
 *   - 05/09/2006 ACG Now outputing info on each screw for all timesteps
 *                    for each screw, up until the screw is broken, even if
 *                    there are no dead elements in the screw. (WSK wanted this)
 *   - 05/09/2006 ACG handles multiple screw meshes (specified on teh command line)
 *   - 05/22/2006 WSD In no mesh names are supplied on the command line, this
 *                    will attemp to find all meshes with 'screw' or 'scrw'
 *                    in the name (case insensitive).
 *   - 05/22/2006 WSD Fiddled with output. Now there is always exactly 1 line
 *                    of output per screw per timestep. At the first break
 *                    it does not report BR any more. And after the first break
 *                    it know prints a line "*** still broken ***".
 *   - 06/02/2006 ACG updated for new fc_getSubsetSies interface 
 *   - 06/12/2006 WSD Prints values of sequence. Prints screw bounding boxes.
 *   - 06/28/2006 ACG uses new shape interface
*/

#include <string.h>
#include "fc.h"
#include <stdio.h>


typedef struct {
  //keeping everything for future use, although now we will just want the "base".
  FC_Subset whole;
  FC_Shape shape;
  int baseID;
  double baseArea;
  FC_Coords baseNormal;
  int broken;
  int eroded;
} screwinfo;

typedef struct{
  char* name;
  char* op;
  double val;
  FC_Variable *elemdeath;
} deathvarinfo;


static int parser(int argc, char**argv, 
		  char** file_name, int* numscrewmeshes, char ***screwmesh_names,
		  FC_VerbosityLevel *verbose_level, 
		  int* numdeathvars, deathvarinfo** deathvars,
		  int *seg_dim);
static int hasScrewInStr(char* name);
static int parseDeathVarSet(int argc, char** argv, int *idx, int *numvars,
			    deathvarinfo** deathvars);
static int printDeathVarSet(int numvars, deathvarinfo *deathvars);
static int segmentMesh(FC_Mesh mesh, int shared_dim,
		       int *numsegs, FC_Subset **segs);
static int screwInit(screwinfo *screw);
static int printScrew(screwinfo screw, char* screwmeshname, int screwID, int verb);
static int buildScrews(FC_Mesh mesh, int *numScrews, screwinfo** screws);
static int getDeathThresholdSubset(FC_Mesh mesh,int numdeathvars,
				   deathvarinfo* deathvars,
				   int timestep,
				   FC_Subset* threshold);
static double projection(FC_Subset deadscrew, screwinfo screw);
static int broken(FC_Subset deadscrew, int seg_dim, screwinfo* screw);
static int eroded(FC_Subset deadscrew, screwinfo* screw);


static int parser(int argc, char** argv, char** file_name,
		  int* numscrewmeshes, char ***screwmesh_names,
		  FC_VerbosityLevel* verbose_level,
		  int *numdeathvars, deathvarinfo** deathvars,
		  int *seg_dim){
  FC_ReturnCode rc;
  char** names;
  int i, j, nmesh;

  // handle arguments
  if (argc < 2){
  usage:
    printf("Usage: %s  [options] file [screwmesh_name(s)]\n", 
	   argv[0]);
    printf("options: \n");
    printf("   -h            : print this help message\n");
    printf("   -v            : verbose: print warning and error messages\n");
    printf("   -V            : very verbose: prints log and error messages\n");
    printf("   -n            : number of vars that will specify element death\n");
    printf("                   followed by var_name, op, value triads. Default is\n");
    printf("                   Only 1 triad of \"elem_death\" <= 0\n");
    printf("   -b seg_dim    : segmentation_dimension for determining breakage\n");
    printf("                   (default is 0)\n");
    printf("\n");
    printf("Examples:\n");
    printf("   ./screwBreaks -b 1 -n 2 elem_death \"<=\" 0 elem_var_8 \">\" 0.1 ~wkoegle/Work/fcdmf/Data/screw-preload/d3plot Material_ID_1 \n");
    printf("   ./screwBreaks ~wkoegle/Work/fcdmf/Data/hatch_run1/d3plot\n");
    printf("   ./screwBreaks ~wkoegle/Work/fcdmf/Data/run30-EMMI-screws/d3plot Material_ID_25 Material_ID_49 > junk; grep \"Screw 0\" junk\n");
    printf("   ./screwBreaks -n 1 damage \">=\" 0.5 ../data/gen_screws.ex2 \n");
    printf("\n");
    printf("Output: Prints Breakage Ratio and notification of when screws first break.\n");
    printf("\n");
    printf("Note: Screws can be normal screw shape or cylindrical etc. Only constraint is\n");
    printf("that they have an identifiable base for the Breakage Ratio calculation.\n");
    exit(-1);
  }


  nmesh = 0;
  names = NULL;
  for (i = 1; i < argc; i++) {
    if (!strcmp(argv[i], "-v")) {
      *verbose_level = FC_WARNING_MESSAGES;
    }
    else if (!strcmp(argv[i], "-V")) {
      *verbose_level = FC_LOG_MESSAGES;
    }
    else if (!strncmp(argv[i], "-h", 2))
      goto usage;
    else if (!strcmp(argv[i], "-n")) {      
      i++;
      rc = parseDeathVarSet(argc,argv,&i, numdeathvars,deathvars);
      if (rc!= FC_SUCCESS){
	goto usage;
      }
    }
    else if (!strcmp(argv[i], "-b")) {
      i++;
      if (i < argc){
	*seg_dim = atoi(argv[i]);
      }else{
	goto usage;
      }
    }else{ 
      *file_name = argv[i];
      // optional mesh names
      for (j = i+1; j < argc; j++) {
	void* tmp;
	if (!strncmp(argv[j],"-",1))
	  break;
	tmp = (char**)realloc(names, (nmesh+1)*sizeof(char*));
	if (!tmp)
	  fc_exitIfError(FC_MEMORY_ERROR);
	names = tmp;
	names[nmesh] = (char*)malloc((strlen(argv[j])+1)*sizeof(char));
	strcpy(names[nmesh], argv[j]);
	nmesh++;
	i++;
      }
    }
  }
  *screwmesh_names = names;
  *numscrewmeshes = nmesh;

  // test args
  if (!(*file_name))
    goto usage;
  if (*seg_dim < 0 || *seg_dim > 2){
    fc_printfErrorMessage("Invalid seg dim %d\n",*seg_dim);
    exit(-1);
  }

  if (*numdeathvars == 0){ //default vals
    *numdeathvars = 1;
    (*deathvars) = (deathvarinfo*)malloc(sizeof(deathvarinfo));
    if (deathvars == NULL){
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
      exit(-1);
    }
    (*deathvars)[0].name = "elem_death";
    (*deathvars)[0].op = "<=";
    (*deathvars[0]).val = 0;
  }

  if (0){
    printDeathVarSet(*numdeathvars,*deathvars);
  }

  return FC_SUCCESS;
}

// if any case of 'screw' or 'scrw' in name, return 1 (true) 
static int hasScrewInStr(char* name) {
  int i;
  int len = strlen(name);

  // scan string
  for (i = 0; i < len-3; i++) {
    if (name[i] == 's' || name[i] == 'S') {
      if (!strncasecmp(&name[i], "screw", 5) ||
	  !strncasecmp(&name[i], "scrw", 4)) {
	return 1;
      }
    }
  }

  // didnt' find one
  return 0;
}

static int parseDeathVarSet(int argc, char** argv, int *idx, int *numvars, 
			    deathvarinfo **deathvars){
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


static int printDeathVarSet(int numvars, deathvarinfo *deathvars){
  int i; 

  printf("Deathvars (%d):\n",numvars);
  for (i = 0; i < numvars; i++){
    printf("%s %s %12.4f\n",deathvars[i].name,deathvars[i].op,
	   deathvars[i].val);
  }
  return FC_SUCCESS;
}

static int segmentMesh(FC_Mesh mesh, int shared_dim,
		       int *numsegs, FC_Subset **segs){
  FC_ReturnCode rc;
  FC_Subset whole;
  int numMem;
  int k;

  rc = fc_getMeshNumElement(mesh, &numMem);
  fc_exitIfErrorPrintf(rc, "cant get mesh numelem");

  rc = fc_createSubset(mesh, "whole", FC_AT_ELEMENT, &whole);
  fc_exitIfErrorPrintf(rc, "cant create whole mesh");

  for (k = 0; k < numMem; k++){
    rc = fc_addMemberToSubset(whole,k);
    fc_exitIfErrorPrintf(rc, "cant add member to subset");
  }

  rc = fc_segment(whole, shared_dim,numsegs,segs);
  fc_exitIfErrorPrintf(rc, "cant segment mesh");

  fc_deleteSubset(whole);

  return FC_SUCCESS;
}


static int screwInit(screwinfo *screw){
  (*screw).whole = FC_NULL_SUBSET;
  //shape gets allocation defaults
  (*screw).baseID = -1;
  (*screw).baseArea = 0;
  (*screw).baseNormal[0] = 0;
  (*screw).baseNormal[1] = 0;
  (*screw).baseNormal[2] = 0;
  (*screw).broken = 0;
  (*screw).eroded = 0;

  return FC_SUCCESS;
}


static int printScrew(screwinfo screw, char* screwmeshname, int screwID, int verb){
  int i;

  printf("Mesh %20s Screw %d: numSides = %d, baseID = %d base area = %6.2f proj norm = (%3.2f, %3.2f, %3.2f)\n",
	 screwmeshname, screwID,screw.shape.numSides,screw.baseID,screw.baseArea,
	 screw.baseNormal[0],screw.baseNormal[1],screw.baseNormal[2]);
  if (verb){
    for (i = 0; i  < screw.shape.numSides; i++){
      if (i == screw.baseID){
	printf("Side %d:\n",i);
	fc_printSubset(screw.shape.faces[i],"faces",1);
	fc_printSubset(screw.shape.elems[i],"elems",1);
      }
    }
  }
  return FC_SUCCESS;
}


static int buildScrews(FC_Mesh mesh, int *nScrews, screwinfo** scr){
  FC_ReturnCode rc;
  FC_Subset *screwsegs;
  int numScrews;
  screwinfo *screws;

  int screwsideangle = 80;
  int screwshareddim = 0;
  int i;

  
  //segment mesh to get the screws
  //cant use the screw function becuase 1) we want the whole segments for
  //the projection 2) this does the plug shaped screws as well
  rc = segmentMesh(mesh,0, &numScrews, &screwsegs);
  fc_exitIfErrorPrintf(rc, "Can't segment mesh to get screws");
  //  printf("num screws = %d\n",numScrews);

  screws = malloc(numScrews*sizeof(screwinfo)); 
  if (screws == NULL){
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    exit(-1);
  }

  for (i = 0; i < numScrews; i++){
    int numEnds, *endArray;
    FC_Shape *shapes;
    int numShapes;

    screwInit(&(screws[i]));

    screws[i].whole = screwsegs[i];
    //this will get them in some random order
    rc = fc_getSubsetShapes(screws[i].whole,screwsideangle,screwshareddim,
				&numShapes,&shapes);
    if (rc != FC_SUCCESS || numShapes <= 1){
      fc_exitIfErrorPrintf(rc, "Can't get screw shape - aborting");
    }
    
    screws[i].shape = shapes[0];
    free(shapes); //but not the shape

    //now get the smallest end for the projection
    rc = fc_getShapeEnds(&screws[i].shape, &numEnds,&endArray);
    if (rc != FC_SUCCESS){
      fc_exitIfErrorPrintf(rc, "Can't get screwsubset ends - aborting");
    }

    // early exit
    // Question: we are rejecting the whole mesh, but we could just dump 
    // non-screw segemnts? YES we could, but not today AG 6/28/06
    if (numEnds <= 1) {
      // FIX: probably have within screw data to cleanup
      // need a screwFinal
      fc_freeShape(&screws[i].shape);
      free(screws);
      free(screwsegs);
      return 1;
    } 

    screws[i].baseID = endArray[0];
    free(endArray);

    rc = fc_getSubsetArea(screws[i].shape.faces[screws[i].baseID],
			  &screws[i].baseArea);
    fc_exitIfErrorPrintf(rc,"Can't get screwbase area");

    rc = fc_getPlanarSideNormal(&screws[i].shape,screws[i].baseID,
				&(screws[i].baseNormal));
    fc_exitIfErrorPrintf(rc,"Can't get screwbase normal");


  }
  
  *nScrews = numScrews;
  *scr = screws;

  free(screwsegs);

  return FC_SUCCESS;
}


static int getDeathThresholdSubset(FC_Mesh mesh, int numdeathvars,
				   deathvarinfo* deathvars, 
				   int timestep,
				   FC_Subset* threshold){
  FC_ReturnCode rc;
  int k;

  rc = fc_createSubset(mesh,"threshold",FC_AT_ELEMENT,threshold);
  fc_exitIfErrorPrintf(rc,"failed to create empty threshold subset");
  //get dead regions for this timestep
  for (k = 0; k < numdeathvars; k++){
    FC_Subset deathset;
    int *arr, nummem;

    rc = fc_createThresholdSubset(deathvars[k].elemdeath[timestep], 
				  deathvars[k].op,
				  deathvars[k].val,
				  "temp", &deathset);
    fc_exitIfErrorPrintf(rc, "Failed to theshold using op '%s' and value '%g'",
			 deathvars[k].op, deathvars[k].val);
    
    rc=fc_getSubsetMembersAsArray(deathset,&nummem,&arr);
    fc_exitIfErrorPrintf(rc, "cant get dead subset members");
    rc = fc_addArrayMembersToSubset(*threshold,nummem,arr);
    fc_exitIfErrorPrintf(rc, "cant add dead subset members to threshold subset");
    
    free(arr);
    fc_deleteSubset(deathset);
  }

  return FC_SUCCESS;
}


static double projection(FC_Subset deadscrew, screwinfo screw){
  FC_ReturnCode rc;
  FC_Subset maskedScrewBaseIDS;
  FC_Subset decayedBase;
  double ret = 0.0;
  int num1;

  //  fc_printSubset(deadscrew,"deadset",0);
  //  fc_printSubset(screw.elems[screw.baseID],"screwbase",1);
  rc = fc_getSubsetNumMember(deadscrew,&num1);
  fc_exitIfErrorPrintf(rc, "cant get subset num members");
  if (num1){
    rc = fc_projectedSubsetElementBoundingBoxIntersection(deadscrew,
							  screw.shape.elems[screw.baseID],
							  screw.baseNormal,
							  &maskedScrewBaseIDS);
    fc_exitIfErrorPrintf(rc, "Failed to get projected BBox Intersection");
    
    //Now we want the areas of the maskedIDS relative to the area of the
    //full base.
    rc = fc_getDecayedSkin(maskedScrewBaseIDS,&screw.shape.faces[screw.baseID],
			   &decayedBase); //go from elem->faces
    fc_exitIfErrorPrintf(rc, "Failed to get decayed skin");
    if (!FC_HANDLE_EQUIV(decayedBase,FC_NULL_SUBSET)){
      rc = fc_getSubsetArea(decayedBase,&ret);
      fc_exitIfErrorPrintf(rc, "Failed to get areas");
      fc_deleteSubset(decayedBase);
    }
    fc_deleteSubset(maskedScrewBaseIDS);
  }

  return ret;
}


static int broken(FC_Subset deadscrew, int seg_dim, screwinfo *screw){
  FC_ReturnCode rc;
  FC_Subset *livesegments;
  int numliveseg;
  int k;

  //segment live regions
  rc = fc_subsetSegmentsSubset(deadscrew, screw->whole,
			       seg_dim,"live",&numliveseg, &livesegments);
  fc_exitIfErrorPrintf(rc, "Failed in subsetSegmentsSubset");
	
  if (numliveseg > 1){
    screw->broken = 1;
  }
	
  for (k = 0; k < numliveseg; k++){
    fc_deleteSubset(livesegments[k]);
  }
  free(livesegments);

  return screw->broken;
}


static int eroded(FC_Subset deadscrew, screwinfo *screw){
  FC_ReturnCode rc;
  int *flag;
  int k;

  rc = fc_getShapeSidesDecayType(deadscrew,&(screw->shape),&flag);
  fc_exitIfErrorPrintf(rc, "Can't get decayed shape ");
  for (k = 0; k < screw->shape.numSides; k++){
    if (flag[k] ==2){
      screw->eroded = 1;
      break;
    }
  }
  free(flag);
  return screw->eroded;
}



/*******************************************************************
 *
 *   MAIN
 *
 ******************************************************************/
int main(int argc, char** argv) {
  FC_ReturnCode rc;
  FC_Dataset dataset;
  FC_Sequence *sequences, sequence;
  FC_Mesh* meshes, *screwmeshes;

  screwinfo *screws;
  int numMesh, numSeq, numStep;
  int numScrews;

  int seg_dim = 0;
  char* file_name = NULL;
  char** screwmesh_names;
  int numscrewmeshes = 0;
  int numdeathvars = 0;
  deathvarinfo *deathvars;
  int tscrews = 0;
  int tscrewsbroken = 0;

  FC_VerbosityLevel verbose_level = FC_QUIET;
  int i,j,imesh;
  
  rc = parser(argc,argv, &file_name, &numscrewmeshes,&screwmesh_names,
	      &verbose_level, &numdeathvars, &deathvars, &seg_dim);
  fc_exitIfErrorPrintf(rc,"can't parse args");

  // init library and load dataset 
  rc = fc_setLibraryVerbosity(verbose_level);
  fc_exitIfError(rc);
  rc = fc_initLibrary();
  fc_exitIfError(rc);
  rc = fc_loadDataset(file_name, &dataset);
  if (rc!=FC_SUCCESS){
    printf("Failed to load dataset %s\n",file_name);
    fc_exitIfErrorPrintf(rc, "Failed to load dataset '%s'", file_name);
  }
  
  // Get number of meshes and mesh names
  rc = fc_getMeshes(dataset, &numMesh, &meshes);
  fc_exitIfError(rc);

  // Get the screw meshes
  screwmeshes = (FC_Mesh*)malloc(numMesh*sizeof(FC_Mesh)); // numMesh is max
  if (!screwmeshes)
    fc_exitIfError(FC_MEMORY_ERROR);
  if (numscrewmeshes > 0) {
    for (imesh = 0; imesh < numscrewmeshes; imesh++) {
      FC_Mesh *returnMeshes;
      int numReturnMeshes;
      rc = fc_getMeshByName(dataset, screwmesh_names[imesh], 
			    &numReturnMeshes,&returnMeshes);
      fc_exitIfErrorPrintf(rc, "Failed to get mesh by name");
      if (numReturnMeshes != 1){
	fc_exitIfErrorPrintf(rc, "Problems finding (unique) mesh '%s' - found %d matches",
			     screwmesh_names[imesh],numReturnMeshes);
      }
      screwmeshes[imesh] = returnMeshes[0];
      free(returnMeshes);
    }
  }
  else {
    char* temp_name;
    screwmesh_names = (char**)malloc(numMesh*sizeof(char*));
    if (!screwmesh_names)
      fc_exitIfError(FC_MEMORY_ERROR);
    for (imesh = 0; imesh < numMesh; imesh++) {
      rc = fc_getMeshName(meshes[imesh], &temp_name);
      if (hasScrewInStr(temp_name)) {
	screwmeshes[numscrewmeshes] = meshes[imesh];
	screwmesh_names[numscrewmeshes] = temp_name;
	numscrewmeshes++;
      }
      else {
	free(temp_name);
      }
    }
    if (numscrewmeshes < 1) 
      fc_exitIfErrorPrintf(FC_ERROR, 
	 "No meshes with 'screw' or 'scrw' in the name were found");
  }
  if (numscrewmeshes < numMesh) {
    realloc(screwmeshes, numscrewmeshes*sizeof(FC_Mesh)); // shrink to actual
  }

  // get the sequence
  rc = fc_getSequences(dataset, &numSeq, &sequences);
  fc_exitIfErrorPrintf(rc, "failed to get sequences");
  if (numSeq < 1) {
    fc_exitIfErrorPrintf(FC_ERROR, "expecting at least one sequence");
  }
  else if (numSeq > 1) {
    fc_printfWarningMessage("Expected only 1 sequence not %d, will use "
			    "first one", numSeq);
  }
  sequence = sequences[0];
  free(sequences);
  
  // header
  printf("Screw characterizations for dataset '%s'\n", file_name);
  if (numMesh < 0) {
    printf("No meshes found\n");
    fc_printfErrorMessage("No meshes were found\n");
    exit(-1);
  }

  printf("Sequence values:\n");
  {
    double* coords;
    rc = fc_getSequenceNumStep(sequence, &numStep);
    fc_exitIfErrorPrintf(rc, "failed to get sequence numStep");
    rc = fc_getSequenceCoordsAsDataType(sequence, FC_DT_DOUBLE, 
					(void**)&coords);
    fc_exitIfErrorPrintf(rc, "failed to get coords as double");
    for (i = 0; i < numStep; i++)
      printf("  Step %d: %g\n", i, coords[i]);
    fflush(NULL);
    free(coords);
  }

  printf("%d mesh(es):\n", numscrewmeshes);
  fflush(NULL);
  
  for (imesh = 0; imesh < numscrewmeshes; imesh++){
    char* screwmesh_name = screwmesh_names[imesh];
    int numBrokePerStep[numStep];
    int* brokenScrewIDs[numStep];

    // initialize summary info
    for (i = 0; i < numStep; i++) {
      numBrokePerStep[i] = 0;
      brokenScrewIDs[i] = NULL;
    }

    for (i = 0; i < numdeathvars; i++){
      rc = fc_getOrGenerateUniqueSeqVariableByName(screwmeshes[imesh], deathvars[i].name,
						   &numStep, &deathvars[i].elemdeath);
      fc_exitIfErrorPrintf(rc, "Failed to find element death variable '%s'",
			   deathvars[i].name);
      if (!deathvars[i].elemdeath){
	fc_exitIfErrorPrintf(FC_INPUT_ERROR,
			     "Failed to find element death variable '%s'",
			     deathvars[i].name);
      }
    }
    
    rc = buildScrews(screwmeshes[imesh],&numScrews,&screws);
    if (rc == 1) {
      printf("WARNING: Mesh %s does not look like screws and will be ignored\n",
	     screwmesh_names[imesh]);
      continue;
    }
    fc_exitIfErrorPrintf(rc,"Cant build the screws");
    

    printf("Mesh %d: '%s' has %d screws\n", imesh, screwmesh_name, numScrews);
    printf("Screw bounding boxes at Step 0:\n");
    for (i = 0; i < numScrews; i++) {
      int numDim;
      FC_Coords lowers, uppers;
      rc = fc_getSubsetBoundingBox(screws[i].whole, &numDim, &lowers, &uppers);
      fc_exitIfErrorPrintf(rc, "Failed to get screw bounding box");
      printf("  Screw %d: [ %g", i, lowers[0]);
      for (j = 1; j < numDim; j++)
	printf(", %g", lowers[j]);
      printf(" ] - [ %g", uppers[0]);
      for (j = 1; j < numDim; j++)
	printf(", %g", uppers[j]);
      printf(" ]\n");
    }
    fflush(NULL);
    
    //------------------dead analysis----------------------------------//
    for (j = 0; j < numStep; j++){
      FC_Subset threshold;
      int num;
      
      // break out if no unbroken or eroded screws left to check
      num = 0;
      for (i = 0; i < numScrews; i++){
	if (!(screws[i].broken || screws[i].eroded)){
	  num = 1;
	  break;
	}
      }
      if (num == 0){
	break;
      }

      rc = getDeathThresholdSubset(screwmeshes[imesh],numdeathvars,
				   deathvars,j,&threshold);
      fc_exitIfErrorPrintf(rc,"cant get death subset");
      //    fc_printSubset(threshold,"communal",0);
    
      for (i = 0; i < numScrews; i++){
	FC_Subset deadscrew;
	double projArea;
	int broke = 0, erode = 0;
	
	if (!screws[i].broken && !screws[i].eroded){ 
	
	  //get the dead part of the screw at this time step
	  rc = fc_createSubsetIntersection(threshold,"AND",screws[i].whole,"dead",
					   &deadscrew);
	  fc_exitIfErrorPrintf(rc, "Cant get dead screw bits");
	  
	  rc = fc_getSubsetNumMember(deadscrew,&num);
	  fc_exitIfErrorPrintf(rc, "Cant get num dead elements");
	  
	  
	  //-------------------broken-------------------------------------//
	  broke = broken(deadscrew, seg_dim, &screws[i]);
	  
	  //-------------------eroded-------------------------------------//
	  if (!screws[i].broken) {
	    erode = eroded(deadscrew,&screws[i]);
	  }

	  //-------------------projection-------------------------------------//
	  if (!screws[i].broken && !screws[i].eroded) {
	    projArea = projection(deadscrew,screws[i]);
	  }

	  fc_deleteSubset(deadscrew);
	}

	// print results
	printf("Mesh %-20s Screw %-3d Step %-3d ", screwmesh_name, i, j);
	if (broke) {
	  printf("First broken");
	}
	else if (erode) {
	  printf("First broken, side eroded");
	}
	else if (screws[i].broken || screws[i].eroded) {
	  printf("*** still broken ***");
	} 
	else {
	  printf("BR =%6.2f (%7.2f/%7.2f)",(projArea/screws[i].baseArea),projArea,
		 screws[i].baseArea);
	}
	printf("\n");
	
	// collect summary info
	if (broke || erode) {
	  void* temp;
	  temp = realloc(brokenScrewIDs[j], (numBrokePerStep[j]+1)*sizeof(int));
	  if (!temp) 
	    fc_exitIfError(FC_MEMORY_ERROR);
	  brokenScrewIDs[j] = temp;
	  brokenScrewIDs[j][numBrokePerStep[j]] = i;
	  numBrokePerStep[j]++;
	}

	//diagnostic
	if (0){
	  if (screws[i].broken || screws[i].eroded){
	    printScrew(screws[i],screwmesh_name,i,1);
	  }
	}

      } //for each screw
      fc_deleteSubset(threshold);
    } //for each timestep

    //final writeout    
    i = 0;
    for (j = 0; j < numScrews; j++){
      if (screws[j].broken || screws[j].eroded){
	i++;
      }
    }
    printf("Mesh %-20s Broken/Total screws: %d/%d\n",screwmesh_name, i,numScrews);
    tscrews+=numScrews;
    tscrewsbroken +=i;

    // print summary info
    for (i = 0; i < numStep; i++) {
      if (numBrokePerStep[i] > 0) {
	printf("  %d screw(s) broke at Step %d: %d", numBrokePerStep[i],
	       i, brokenScrewIDs[i][0]);
	for (j = 1; j < numBrokePerStep[i]; j++)
	  printf(", %d", brokenScrewIDs[i][j]);
	printf("\n");
	free(brokenScrewIDs[i]);
      }
    }

    //cleanup after this mesh
    for (i = 0; i < numdeathvars; i++){ //release the deathvar on this mesh
      fc_deleteSeqVariable(numStep,deathvars[i].elemdeath);
      free(deathvars[i].elemdeath);
    }
    
    for (i = 0; i < numScrews; i++){
      fc_freeShape(&screws[i].shape);
      fc_deleteSubset(screws[i].whole);
    }
    free(screws);
  }
  printf("All %d Mesh(es) %10s Broken/Total screws: %d/%d\n",numscrewmeshes," ",
	 tscrewsbroken,tscrews);

  //clean up after all meshes
  free(meshes);
  for (i = 0; i < numscrewmeshes; i++)
    free(screwmesh_names[i]);
  free(screwmesh_names);
  // deathvars[i].elemdeath has already been freed
  free(deathvars);      //dont need to free name or op since assigned from argv 

  // shut down
  fc_deleteDataset(dataset);
  fc_finalLibrary();

  exit(0);

}

