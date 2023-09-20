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
 * \file gaplines.c
 * \brief Look at gaps by creatings lines connecting the surfaces of 
 *    two meshes.
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/tools/gaplines.c,v $
 * $Revision: 1.20 $ 
 * $Date: 2006/09/28 06:34:39 $
 *
 * \description
 *
 *   First have to figure out where gaps might happen -- i.e. find surfaces of
 *   different measures that initially are abutted. Then create lines
 *   (initially of 0 length) between the surfaces. If the lines lengthen when
 *   the meshes are deformed, then you have a gap.
 *
 *   (In future, thinking of creating elements of zero volume so that we can
 *   better characterize number of volumes of gaps.)
 *
 *   Input
 *     - Two meshes & their displacement seq vars
 *     - Min distance that points can be from each other and still be
 *       considered to be touching.
 *
 *   Assumptions: 
 *     - Meshes do not interpenetrate (more than some very small amount). 
 *     - Faces are planar
 *     - Can determine adjacent surfaces from initial geometry information.
 *
 *   Algorithm:
 *     - Get outside surfaces of each mesh - i.e. skin
 *     - For each mesh, find all vertices that are a) close to a skin face on
 *       the other mesh (close means than if you project the vertex to the 
 *       face, it lies within or on boundary of face AND the distance along 
 *       that projection is less than the min distance) and b) are in a skin 
 *       face themselves that have the same orientation as the face on the 
 *       other mesh. Keep track of the vertex, the faces on both meshes, and 
 *       the parameterized location of the projected vertex within the face on
 *       the other mesh.
 *     - Create a new mesh which is lines from each vertex to the paramterized
 *       point in the corresponding face on the other mesh. (Line lengths 
 *       should all be less than min distance chosen in 1.)
 *     - New! using Shapes, create subsets of the line mesh for each side-side
 *       grouping of gap lines.  
 *     - Create displacement variable for the line mesh.
 *     - Create variable which is the average normal of the two faces the
 *       line connects
 *
 *   Output:
 *     - Report #, and min, max, average and standard dev of line lengths per 
 *       step for each set of lines and the overall set of lines. Also report 
 *       normal and tangent components.
 *     - Write a gaplines dataset
 *
 * /modifications:
 *    - 01/05 Created, WSK.
 *    - 02/07/06 Had to update to changes in library.
 *    - 06/12/2006 WSD Making into an official tool.
 *    - 08/25/2006 WSD Adding normal & tangent break down of gapline line
 *      lengths (with respect to face normals). Also added break down by sides
 *      using FC_Shape.
 *    - 08/31/2006 WSD Added ability to drop gap lines connected to dead
 *      elements. 
 */
#include <string.h>
#include "fc.h"
#include "fcP.h"

// hold info about death vars
typedef struct{
  char* name;
  char* op;
  double val;
} deathvarinfo;

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

// get max edge length
static FC_ReturnCode getMeshMaxEdgeLength(FC_Mesh mesh, double* max_len) {
  FC_ReturnCode rc;
  FC_Variable lengths;
  
  rc = fc_getEdgeLengths(mesh, &lengths);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_getVariableMinMax(lengths, NULL, NULL, max_len, NULL);
  return rc;
}

// gets the first the orient of the first elem for that face
// assuming that skins will have unique parent element
// other faces are place holders
static FC_ReturnCode calcFaceOrients(FC_Mesh mesh, int numFace, 
				     int* faceIDs, int** faceOrients) {
  FC_ReturnCode rc;
  int i, j;
  int numTotalFace;
  int numFacePerElem, *elemToFaceConns;
  int* elemFaceOrients;
  int* numElemPerFace, **elemParentsPerFace;
  FC_ElementType elemtype;

  // gather info
  rc = fc_getMeshNumFace(mesh, &numTotalFace);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_getMeshElementType(mesh, &elemtype);
  if (rc != FC_SUCCESS)
    return rc;
  numFacePerElem = fc_getElementTypeNumFace(elemtype);
  rc = fc_getMeshElementToFaceConnsPtr(mesh, &elemToFaceConns);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_getMeshElementFaceOrientsPtr(mesh, &elemFaceOrients);
  if (rc != FC_SUCCESS)
    return rc;
  rc = _fc_getMeshElementParentsOfFacesPtr(mesh, &numElemPerFace,
					   &elemParentsPerFace);
  if (rc != FC_SUCCESS)
    return rc;

  // do it
  *faceOrients = (int*)malloc(numTotalFace*sizeof(int));
  if (!(*faceOrients))
    return FC_MEMORY_ERROR;
  for (i = 0; i < numFace; i++) {
    int elemID = elemParentsPerFace[faceIDs[i]][0];

    for (j = 0; j < numFacePerElem; j++) {
      if (elemToFaceConns[elemID*numFacePerElem+j] == faceIDs[i]) {
	(*faceOrients)[faceIDs[i]] = elemFaceOrients[elemID*numFacePerElem+j];
	break;
      }
    }
  }

  return FC_SUCCESS;
}


// Given a mesh of lines and a displacement variable, calculate the lines
// as vectors by subtracting one displaced point from the other
static FC_ReturnCode getLinesAsVectors(FC_Mesh mesh, FC_Variable displVar,
				       FC_Variable *lineVectors) {
  FC_ReturnCode rc;
  int i, j;
  int numElem, numDim, numVperE = 2; // assuming lines
  double* displCoords;
  int* conns;
  double* vectors;

  rc = fc_getMeshNumElement(mesh, &numElem);
  if (rc != FC_SUCCESS)
    return rc;  
  rc = fc_getMeshDim(mesh, &numDim);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_getDisplacedMeshCoords(mesh, displVar, &displCoords);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_getMeshElementConnsPtr(mesh, &conns);
  if (rc != FC_SUCCESS)
    return rc;

  vectors = malloc(numElem*numDim*sizeof(double));
  if (!vectors) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  for (i = 0; i < numElem; i++) {
    int vert1 = conns[i*numVperE];
    int vert2 = conns[i*numVperE+1];
    for (j = 0; j < numDim; j++) {
      vectors[i*numDim+j] = displCoords[vert1*numDim+j] -
	                    displCoords[vert2*numDim+j]; 
    }
  }
  free(displCoords);

  rc = fc_createVariable(mesh, "line vectors", lineVectors);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_setVariableDataPtr(*lineVectors, numElem, numDim, FC_AT_ELEMENT,
			     FC_MT_VECTOR, FC_DT_DOUBLE, (void*)vectors);
  return rc;
}

// find the shape and side that a face is in
static FC_ReturnCode findShapeSide(int faceID, int numShape, FC_Shape* shapes,
				   int* shapeID, int* sideID) {
  int i, j;
  for (i = 0; i < numShape; i++) {
    for (j = 0; j < shapes[i].numSides; j++) {
      if (fc_isMemberInSubset(shapes[i].faces[j], faceID)) {
	*shapeID = i;
	*sideID = j;
	return FC_SUCCESS;
      }
    }
  }

  *shapeID = -1;
  *sideID = -1;
  return FC_ERROR;
}

// In one go, remove vals from sia1 and sia2 that are the same,
// and stick them in sia3
static void getRemoveSIAIntersection(FC_SortedIntArray *sia1,
				     FC_SortedIntArray *sia2, 
				     FC_SortedIntArray *sia3) {
  int i, j;
  i = 0;
  j = 0;
  while (i < sia1->numVal && j < sia2->numVal) {
    if (sia1->vals[i] < sia2->vals[j]) {
      i++;
    }
    else if (sia1->vals[i] > sia2->vals[j]) {
      j++;
    }
    else {
      fc_addIntToSortedIntArray(sia3, sia1->vals[i]);
      _fc_deleteEntryFromSortedIntArray(sia1, i);
      _fc_deleteEntryFromSortedIntArray(sia2, j);
      // don't increment i or j
    }
  }
}

// determine if line connecting elements is dead or not
// only one of the elements has to be dead
static int isItDead(int numDeathVar, deathvarinfo* deathvarinfos, int id1,
		    double** data1, int id2, double** data2) {
  int i;
  int isItDead1 = 0, isItDead2 = 0;  // start off assuming dead

  for (i = 0; i < numDeathVar; i++) {
    int (*compare)(double, double);
    char* compare_str = deathvarinfos[i].op;

    // determine comparison
    // break out at the first dead classification
    if (!strcmp(compare_str, "=") || !strcmp(compare_str, "=="))
      compare = fc_eqd;
    else if (!strcmp(compare_str, "!=") || !strcmp(compare_str, "=!"))
      compare = fc_neqd;
    else if (!strcmp(compare_str, "<"))
      compare = fc_ltd;
    else if (!strcmp(compare_str, "<=") || !strcmp(compare_str, "=<"))
    compare = fc_lteqd;
    else if (!strcmp(compare_str, ">"))
      compare = fc_gtd;
    else if (!strcmp(compare_str, ">=") || !strcmp(compare_str, "=>"))
      compare = fc_gteqd;
    else {
      fc_printfErrorMessage("Unknown comparison string: '%s'", compare_str);
      return FC_INPUT_ERROR;
    }
    
    if (compare(data1[i][id1], deathvarinfos[i].val)) {
      isItDead1 = 1;
      break;
    }
    if (compare(data2[i][id2], deathvarinfos[i].val)) {
      isItDead2 = 1;
      break;
    }
  }

  if (isItDead1 || isItDead2)
    return 1;
  else
    return 0;
}


// print the stats of a specific seqVar over time
// pass in NULL for subset, or pointer to subset
static void printStats(char* title, int numStep, FC_Variable* lengthSeqVar,
		       FC_Subset* inputSubset, FC_Subset* alivePerTime,
		       double* seq_coords) {
  int i;
  FC_Subset temp_subset, subsets[numStep];

  // If not input subset, make a subset that is the entire mesh
  if (inputSubset)
    temp_subset = *inputSubset;
  else {
    int numMember;
    FC_Mesh mesh;
    fc_getMeshFromVariable(lengthSeqVar[0], &mesh);
    fc_getMeshNumElement(mesh, &numMember);
    fc_createSubset(mesh, "temp", FC_AT_ELEMENT, &temp_subset);
    for (i = 0; i < numMember; i++)
      fc_addMemberToSubset(temp_subset, i);
  }

  // If death is possible, find the ones that are still alive
  if (alivePerTime) {
    for (i = 0; i < numStep; i++)
      fc_createSubsetIntersection(temp_subset, "&&", alivePerTime[i],
				  "temp", &subsets[i]);
  }
  else {
    for (i = 0; i < numStep; i++)
      subsets[i] = temp_subset;
  }

  // print the stats
  printf("----------------------------------------------------------------\n"); 
  printf("    Step      |      |    %s\n", title);
  printf("--------------|------|------------------------------------------\n");
  printf("ID     Value  |  num |    min        max        mean      stdev\n");
  printf("----------------------------------------------------------------\n");
  for (i = 0; i < numStep; i++) {
    int num;
    double min, max, mean, stdev;
    fc_getSubsetNumMember(subsets[i], &num);
    fc_getVariableSubsetMinMax(lengthSeqVar[i], subsets[i], &min, NULL, 
			       &max, NULL);
    fc_getVariableSubsetMeanSdev(lengthSeqVar[i], subsets[i], &mean, &stdev);
    printf("%-2d %10f %6d %10f %10f %10f %10f\n", i, seq_coords[i], num,
	   min, max, mean, stdev);
    fflush(NULL);
  }

  if (!inputSubset)
    fc_deleteSubset(temp_subset);
  if (alivePerTime) {
    for (i = 0; i < numStep; i++)
      fc_deleteSubset(subsets[i]);
  }
}

int main(int argc, char** argv) {
  FC_ReturnCode rc;
  FC_VerbosityLevel verbose_level = FC_QUIET;
  int i, j, k, m, n;
  int numDeathVar = 0;
  deathvarinfo* deathVarInfos;
  FC_Dataset dataset;
  FC_Sequence sequence, sequences[2];
  FC_Mesh meshes[2], gap_mesh;
  FC_Variable *displs[2], **deathVars[2];
  FC_Variable *gap_displ, *aveFaceNorms, *isDeadVar = NULL;
  FC_Variable *lengths, *norms, *tangs;
  FC_Subset face_skins[2], vert_skins[2];
  FC_Subset *alivePerTime;
  char* file_name = NULL;
  char* mesh_names[2] = { NULL, NULL };
  char* displ_name = NULL;
  char charBuf[1024];
  int numDims[2], numDim, numStep;
  double max_lens[2]; // maximum edge length
  double max_dists[2]; // maximum possible distance between point in quad 
                       // and the vert (in the quad) furthest from point
                       // see also documentation at the calc below
  double *coords[2];
  int *numFacePerVert[2], **faceParentsPerVert[2];
  int *numElemPerFace[2], **elemParentsPerFace[2];
  int *numVertexPerFace[2], *faceToVertConns[2], stride[2];
  int numVertIDs[2], *vertIDs[2], numFaceIDs[2], *faceIDs[2];
  int *faceOrients[2];
  int numKeeps[2], *keepVertIDs[2], *keepThisFaceIDs[2], *keepOtherFaceIDs[2];
  double **keepParams[2];
  FC_VertexBin *bins[2];
  double min_dist = 0.001; // critical distance to identify "joined" meshes
  double* seq_coords;
  double angle = 80; // faces with angles greater than this are on diff faces
  int numShapes[2];
  FC_Shape *shapes[2];
  int numSide;
  FC_Subset *sides;

  // handle arguments
  if (argc < 6) {
  usage:
    printf("usage: %s [options] dataset mesh1name mesh2name displvarname min_dist\n",
	   argv[0]);
    printf("   -n # <triplets> : Specify the rules for finding dead elements.\n");
    printf("                     The first parameter is the number of triplets. Each\n"); 
    printf("                     triplet is the name of variable to test, the comparison\n");
    printf("                     operator, and the value. The default is \n");
    printf("                     '-n 1 elem_death \"<=\" 0', meaning all elements\n");
    printf("                     where the variable 'elem_death' is less than or\n");
    printf("                     equal to zero are treated as being dead.\n");
    printf("   -h              : print this help message\n");
    printf("   -v              : verbose: print calculation status and print warning\n");
    printf("                     and error messages\n");
    printf("   -V              : very verbose: print log and error messages\n");
    printf("\n");
    printf("Report gaps that form between the specified meshes assuming that surfaces\n");
    printf("of interest are those that are less than min_dist apart. This program will\n");
    printf("also create a dataset file 'gaplines.ex2' containing the gap lines.\n"); 
    printf("If dead elements are specified, gaplines attached to dead elements will \n");
    printf("be dropped for stat calcs.\n");
    printf("\n");
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
    else if (!strcmp(argv[i], "-n")) {
      i++;
      rc = parseDeathVarSet(argc, argv, &i, &numDeathVar, &deathVarInfos);
      if (rc != FC_SUCCESS)
	goto usage;
    }
    else {
      if (i+4 >= argc)
        goto usage;
      file_name = argv[i];
      mesh_names[0] = argv[i+1];
      mesh_names[1] = argv[i+2];
      displ_name = argv[i+3];
      min_dist = atof(argv[i+4]);
      i+=4;
    }
  }
  
  if (!file_name || !mesh_names[0] || !mesh_names[1] || !displ_name || 
      min_dist < 0.)
     goto usage;

  // Init library & load dataset
  rc = fc_setLibraryVerbosity(verbose_level);
  fc_exitIfError(rc);
  rc = fc_initLibrary();
  fc_exitIfError(rc);
  rc = fc_loadDataset(file_name, &dataset);
  fc_exitIfErrorPrintf(rc, "Failed to load dataset '%s'", file_name);

  // Get sequence, meshes, displacement seqVars & death seqVars
  for (i = 0; i < 2; i++) {
    FC_AssociationType temp_assoc;
    FC_Mesh *returnMeshes;
    int numReturnMeshes;

    rc = fc_getMeshByName(dataset, mesh_names[i], &numReturnMeshes,
			  &returnMeshes);
    fc_exitIfErrorPrintf(rc,"Failed to get mesh by name");
    if (numReturnMeshes != 1){
      fc_exitIfErrorPrintf(FC_INPUT_ERROR,
			   "Problems finding (unique) mesh '%s' - found %d matches",
			   mesh_names[i], numReturnMeshes);
    }
    meshes[i] = returnMeshes[0];
    free(returnMeshes);

    rc = fc_getMeshDim(meshes[i], &numDims[i]);
    fc_exitIfErrorPrintf(rc, "Failed to get numDim of mesh '%s'",
			 mesh_names[i]);
    rc = fc_getOrGenerateUniqueSeqVariableByName(meshes[i], displ_name,
						 &numStep, &displs[i]);
    fc_exitIfErrorPrintf(rc, "Failed to find or generate seqVar '%s' in mesh '%s'", 
			 displ_name, mesh_names[i]);
    
    if (!fc_isValidDisplacementVariable(meshes[i], displs[i][0]))
      fc_exitIfErrorPrintf(rc, "'%s' is not a valid displacement variable",
			   displ_name);
    rc = fc_getSequenceFromSeqVariable(numStep, displs[i], &sequences[i]);
    fc_exitIfErrorPrintf(rc, "failed to get sequence");
    if (numDeathVar > 0) {
      deathVars[i] = (FC_Variable**)malloc(numDeathVar*sizeof(FC_Variable*));
      if (!deathVars[i])
	fc_exitIfError(FC_MEMORY_ERROR);
      for (j = 0; j < numDeathVar; j++) {
	rc = fc_getOrGenerateUniqueSeqVariableByName(meshes[i],
						     deathVarInfos[j].name,
						     &numStep,
						     &deathVars[i][j]);
	fc_exitIfErrorPrintf(rc, "Failed to find or generate seqVar '%s' in mesh '%s'", 
			     deathVarInfos[j].name, mesh_names[i]);

	fc_getVariableAssociationType(deathVars[i][j][0], &temp_assoc);
	if (temp_assoc != FC_AT_ELEMENT) {
	  fc_exitIfErrorPrintf(FC_ERROR, "Death var '%s' should be of assoc "
			       "%s not %s'", deathVarInfos[j].name,
			       fc_getAssociationTypeText(FC_AT_ELEMENT),
			       fc_getAssociationTypeText(temp_assoc));
	}
      }
    }
  }
  if (numDims[0] != numDims[1]) 
    fc_exitIfErrorPrintf(FC_ERROR, "expect meshes to have same numDim "
			 "(%d!=%d)", numDims[0], numDims[1]);
  numDim = numDims[0];
  if (!FC_HANDLE_EQUIV(sequences[0], sequences[1])) 
    fc_exitIfErrorPrintf(FC_ERROR, 
			 "expected displacements to be on the same sequence");
  sequence = sequences[0];

  // Close everything else to 1) to make room, and 2) when we write the 
  // dataset we only get gapline related stuff
  {
    int temp_numSequence, temp_numMesh, temp_numSubset;
    int temp_numVar, temp_numSeqVar, *temp_numSteps;
    FC_Sequence *temp_sequences;
    FC_Mesh *temp_meshes;
    FC_Subset *temp_subsets;
    FC_Variable *temp_vars, **temp_seqVars;
    rc = fc_getSequences(dataset, &temp_numSequence, &temp_sequences);
    fc_exitIfError(rc);
    for (i = 0; i < temp_numSequence; i++) {
      if (!FC_HANDLE_EQUIV(temp_sequences[i], sequence))
	fc_deleteSequence(temp_sequences[i]);
    }
    free(temp_sequences);
    rc = fc_getMeshes(dataset, &temp_numMesh, &temp_meshes);
    fc_exitIfError(rc);
    for (i = 0; i < temp_numMesh; i++) {
      if (!FC_HANDLE_EQUIV(temp_meshes[i], meshes[0]) &&
	  !FC_HANDLE_EQUIV(temp_meshes[i], meshes[1]))
	fc_deleteMesh(temp_meshes[i]);
    }
    free(temp_meshes);
    for (i = 0; i < 2; i++) {
      fc_getSubsets(meshes[i], &temp_numSubset, &temp_subsets);
      for (j = 0; j < temp_numSubset; j++)
	fc_deleteSubset(temp_subsets[j]);
      free(temp_subsets);
      fc_getVariables(meshes[i], &temp_numVar, &temp_vars);
      for (j = 0; j < temp_numVar; j++)
	fc_deleteVariable(temp_vars[j]);
      free(temp_vars);
      fc_getSeqVariables(meshes[i], &temp_numSeqVar, &temp_numSteps, 
			 &temp_seqVars);
      for (j = 0; j < temp_numSeqVar; j++) {
	int hasMatch = 0;
	if (FC_HANDLE_EQUIV(temp_seqVars[j][0], displs[i][0]))
	  hasMatch = 1;
	else {
	  for (k = 0; k < numDeathVar; k++) {
	    if (FC_HANDLE_EQUIV(temp_seqVars[j][0], deathVars[i][k][0])) {
	      hasMatch = 1;
	      break;
	    }
	  }
	}
	if (!hasMatch) 
	  fc_deleteSeqVariable(temp_numSteps[j], temp_seqVars[j]);
	free(temp_seqVars[j]);
      }
      free(temp_numSteps);
      free(temp_seqVars);
    }
  }

  // ---- done with setup

  // Header
  printf("Generating gaplines for:\n");
  printf("Dataset: '%s'\n", file_name);
  printf("Meshes: '%s' and '%s'\n", mesh_names[0], mesh_names[1]);
  printf("Displ: '%s'\n", displ_name);
  printf("Min Dist: %g\n", min_dist);
  if (numDeathVar > 0) {
    printf("Death Criteria: '%s' %s %g", deathVarInfos[0].name,
	   deathVarInfos[0].op, deathVarInfos[0].val);
    for (i = 1; i < numDeathVar; i++) 
      printf(" or '%s' %s %g", deathVarInfos[i].name,
	     deathVarInfos[i].op, deathVarInfos[i].val);
    printf("\n");
  }
  fflush(NULL);

  // skin the meshes (get face skin and vertex skin)
  // also get shapes
  if (verbose_level > FC_QUIET) {
    printf("Skinning meshes ...\n");
    fflush(NULL);
  }
  for (i = 0; i < 2; i++) {
    rc = fc_getMeshCoordsPtr(meshes[i], &coords[i]);
    fc_exitIfError(rc);
    rc = getMeshMaxEdgeLength(meshes[i], &max_lens[i]);
    fc_exitIfErrorPrintf(rc, "failed to get mesh max edge length");
    // Imagine a 4 sided quad with edge lengths of max_lens[i]
    // squash it - the two vertices farthest apart are 2*max_lens[i] apart
    // Imagine a point min_dist above one of those two vertices
    // max_dists[i] is distance from that point to other vert.
    max_dists[i] = sqrt(4*max_lens[i]*max_lens[i] + min_dist*min_dist);
    rc = fc_getMeshSkin(meshes[i], &face_skins[i]);
    fc_exitIfError(rc);
    rc = fc_getSubsetMembersAsArray(face_skins[i], &numFaceIDs[i], &faceIDs[i]);
    fc_exitIfError(rc);
    sprintf(charBuf, "vert skin %d", i);
    rc = fc_copySubsetWithNewAssociation(face_skins[i], charBuf,
					 FC_AT_VERTEX, 1, &vert_skins[i]);
    fc_exitIfError(rc);
    rc = fc_getSubsetMembersAsArray(vert_skins[i], &numVertIDs[i], &vertIDs[i]);
    fc_exitIfError(rc);
    rc = fc_getMeshFaceConnsPtr(meshes[i], &numVertexPerFace[i], &stride[i],
				&faceToVertConns[i]);
    fc_exitIfError(rc);
    rc = _fc_getMeshFaceParentsOfVerticesPtr(meshes[i], &numFacePerVert[i],
					     &faceParentsPerVert[i]);
    fc_exitIfError(rc);
    rc = fc_getMeshShapes(meshes[i], angle, 0, &numShapes[i], &shapes[i]);
    fc_exitIfError(rc);
    calcFaceOrients(meshes[i], numFaceIDs[i], faceIDs[i], &faceOrients[i]);
    // FIX could move this to later ... don't need until constructing gap mesh
    rc = _fc_getMeshElementParentsOfFacesPtr(meshes[i], &numElemPerFace[i],
					     &elemParentsPerFace[i]);
    fc_exitIfError(rc);
    rc = fc_createMeshVertexBin(meshes[i], &bins[i]);
    fc_exitIfError(rc);
  }

  // Create lists of vert to face pairs
  if (verbose_level > FC_QUIET) {
    printf("Finding vert to face pairs ...\n");
    fflush(NULL);
  }
  for (i = 0; i < 2; i++) {
    int otherMesh = (i == 0 ? 1 : 0);
    numKeeps[i] = 0;
    keepVertIDs[i] = malloc(numVertIDs[i]*sizeof(int));
    keepThisFaceIDs[i] = malloc(numVertIDs[i]*sizeof(int));
    keepOtherFaceIDs[i] = malloc(numVertIDs[i]*sizeof(int));
    keepParams[i] = malloc(numVertIDs[i]*sizeof(double*));
    if (!keepVertIDs[i] || !keepThisFaceIDs[i] || !keepOtherFaceIDs[i] ||
	!keepParams[i]) {
      printf("Memory problem\n");
      exit(-1);
    }
    for (j = 0; j < numVertIDs[i]; j++) {
      // look for a containing face on other mesh
      // (assumption - can find such a face via nearest vertex
      //   this is not true, but true enough; cases this doesn't hold
      //   for are just discarded)
      int thisVertID = vertIDs[i][j];
      int numOtherFace, *otherFaceIDs;
      rc = fc_getEntitiesWithinSphere(meshes[otherMesh], bins[otherMesh],
				      &coords[i][thisVertID*numDim], 
				      max_dists[otherMesh], FC_AT_FACE, 
				      &numOtherFace, &otherFaceIDs);
      fc_exitIfErrorPrintf(rc, "failed to get entities within a sphere");
      for (k = 0; k < numOtherFace; k++) {
	int numParam;
	double *params;
	double distance2, min_dist2 = min_dist*min_dist; // squared distances
	//int otherFaceID = faceParentsPerVert[otherMesh][otherVertID][k];
	int otherFaceID = otherFaceIDs[k];
	int thisFaceID;
	FC_Vector otherNormal;
	FC_Coords otherFaceCoords[4], temp_newPoint;
	int foundMatch = 0;
	// skip if not in the skin
	if (!fc_isMemberInSubset(face_skins[otherMesh], otherFaceID))
	  continue;
	if (numVertexPerFace[otherMesh][otherFaceID] != 4)
	  printf("*** expected quads!\n");
	// setup
	for (m = 0; m < numVertexPerFace[otherMesh][otherFaceID]; m++) {
	  int totherVertID = faceToVertConns[otherMesh][otherFaceID*4+m];
	  memcpy(otherFaceCoords[m], &coords[otherMesh][numDim*totherVertID],
		 numDim*sizeof(double));
	}
	_fc_calcSurfaceNormal(numVertexPerFace[otherMesh][otherFaceID],
			      otherFaceCoords, &otherNormal);
	// check that normal is oposite some face in the first mesh
	for (m = 0; m < numFacePerVert[i][thisVertID]; m++) {
	  thisFaceID = faceParentsPerVert[i][thisVertID][m];
	  FC_Coords thisFaceCoords[4];
	  FC_Vector thisNormal;
	  double anglediff;
	  if (!fc_isMemberInSubset(face_skins[i], thisFaceID)) 
	    continue;
	  for (n = 0; n < numVertexPerFace[i][thisFaceID]; n++)
	    memcpy(thisFaceCoords[n],
		   &coords[i][numDim*faceToVertConns[i][thisFaceID*4+n]],
		   numDim*sizeof(double));
	  _fc_calcSurfaceNormal(numVertexPerFace[i][thisFaceID],
				thisFaceCoords, &thisNormal);
	  fc_calcAngleBetweenVectors(thisNormal, otherNormal, &anglediff);
	  //printf("angle = %g,faceorients = %d %d\n", anglediff,
	  // faceOrients[i][thisFaceID], faceOrients[otherMesh][otherFaceID]);
	  if ((faceOrients[i][thisFaceID] ==
	       faceOrients[otherMesh][otherFaceID] && anglediff > 175) ||
	      (faceOrients[i][thisFaceID] == 
	       -1*faceOrients[otherMesh][otherFaceID] && anglediff < 5)) {
	    foundMatch = 1;
	    break;
	  }
	}
	if (!foundMatch)
	  continue;
	// check to see if it's inside face - within tolerance of DBL_EPSILON
	_fc_calcQuadParams(otherFaceCoords, &coords[i][thisVertID*numDim], 
			   &numParam, &params);
	if (!params || 
	    params[0] < -1*FLT_EPSILON || params[0] > 1 + FLT_EPSILON || 
	    params[1] < -1*FLT_EPSILON || params[1] > 1 + FLT_EPSILON ) {
	  free(params);
	  continue; // not inside the quad
	}
	_fc_calcQuadLocation(otherFaceCoords, numParam, params, 
			     &temp_newPoint);
	// check distance
	fc_calcSquaredEuclideanDistance(&coords[i][thisVertID*numDim], 
					temp_newPoint, numDim, &distance2);
	// Below is friendlier version of "if (distance2 > min_dist2) {"
	if (distance2 - min_dist2 > DBL_EPSILON) {
	  free(params);
	  continue;
	}
	// passed the test, add to list and stop looking
	//printf("keepFace = %d\n", otherFaceID);
	keepVertIDs[i][numKeeps[i]] = thisVertID;
	keepThisFaceIDs[i][numKeeps[i]] = thisFaceID;
	keepOtherFaceIDs[i][numKeeps[i]] = otherFaceID;
	keepParams[i][numKeeps[i]] = params;
	numKeeps[i]++;
	break;
      }
      free(otherFaceIDs);
    }
    realloc(keepVertIDs[i], numKeeps[i]*sizeof(int));
    realloc(keepThisFaceIDs[i], numKeeps[i]*sizeof(int));
    realloc(keepOtherFaceIDs[i], numKeeps[i]*sizeof(int));
    realloc(keepParams[i], numKeeps[i]*sizeof(double*));
  }
  for (i = 0; i < 2; i++) {
    fc_deleteSubset(face_skins[i]);
    fc_deleteSubset(vert_skins[i]);
  }

  // Early exit - We are done if no pairs were found
  if (numKeeps[0] + numKeeps[1] < 1) {
    printf("\n");
    printf("Number of gap lines found = 0\n");
    exit(0);
  }
    
  //debug
  //printf("numKeeps[0] = %d, numKeeps[1] = %d\n", numKeeps[0], numKeeps[1]);
  /*
  { // for debugging - make subsets of thisFaces & otherFaces
    FC_Subset thisSubsets[2];
    FC_Subset otherSubsets[2];
    for (i = 0; i < 2; i++) {
      char buf[1048];
      sprintf(buf, "%s-thisFaces", mesh_names[i]);
      fc_createSubset(meshes[i], buf, FC_AT_FACE, &thisSubsets[i]);
      sprintf(buf, "%s-otherFaces", mesh_names[i]);
      fc_createSubset(meshes[i], buf, FC_AT_FACE, &otherSubsets[i]);
    }
    for (i = 0; i < 2; i++) {
      int otherMesh =(i == 0 ? 1 : 0);
      for (j = 0; j < numKeeps[i]; j++) {
	fc_addMemberToSubset(thisSubsets[i], keepThisFaceIDs[i][j]);
	fc_addMemberToSubset(otherSubsets[otherMesh], keepOtherFaceIDs[i][j]);
      }
    }
  }
  */

  // Create the line mesh & its displacements & isDead var
  if (verbose_level > FC_QUIET) {
    printf("Creating gaplines mesh ...\n");
    fflush(NULL);
  }  
  {
    int vindex, eindex;
    int numElement = numKeeps[0] + numKeeps[1];
    int numVertex = numElement*2;
    int temp_numStep;
    double* new_coords = malloc(numDim*numVertex*sizeof(double));
    int* new_conns = malloc(2*numElement*sizeof(int));
    double *displ_datas[2], **death_datas[2];
    double* new_displ_data, *new_ave_face_data;
    int *new_isDead_data;
    // create coords & conns
    vindex = 0;
    for (i = 0; i < 2; i++) {
      int otherMesh = (i == 0 ? 1 : 0);
      for (j = 0; j < numKeeps[i]; j++) {
	int otherFaceID = keepOtherFaceIDs[i][j];
	FC_Coords otherFaceCoords[4], newPoint;
	// First the vertex
	memcpy(new_coords+vindex*numDim, coords[i]+keepVertIDs[i][j]*numDim, 
	       numDim*sizeof(double));
	new_conns[vindex] = vindex;
	vindex++;
	// Then the parameterized point on the face
	for (k = 0; k < 4; k++) {
	  int totherVertID = faceToVertConns[otherMesh][otherFaceID*4+k];
	  memcpy(otherFaceCoords[k], &coords[otherMesh][numDim*totherVertID],
		 numDim*sizeof(double));
	}
	_fc_calcQuadLocation(otherFaceCoords, 2, keepParams[i][j], 
			     &newPoint);
	memcpy(new_coords+vindex*numDim, newPoint, numDim*sizeof(double));
	new_conns[vindex] = vindex;
	vindex++;
      }
    }
    // create mesh
    sprintf(charBuf, "gap lines");
    rc = fc_createMesh(dataset, charBuf, &gap_mesh);
    fc_exitIfError(rc);
    rc = fc_setMeshCoordsPtr(gap_mesh, numDim, numVertex, new_coords);
    fc_exitIfError(rc);
    rc = fc_setMeshElementConnsPtr(gap_mesh, FC_ET_LINE, numElement,
				   new_conns);
    fc_exitIfError(rc);
    // create displacement var and face normals var
    rc = fc_createSeqVariable(gap_mesh, sequence, displ_name, &temp_numStep,
			      &gap_displ);
    fc_exitIfError(rc);
    rc = fc_createSeqVariable(gap_mesh, sequence, "ave face normals", 
			 &temp_numStep, &aveFaceNorms);
    fc_exitIfError(rc);
    if (numDeathVar > 0) {
      rc = fc_createSeqVariable(gap_mesh, sequence, "isDead",
			   &temp_numStep, &isDeadVar);
      fc_exitIfError(rc);
    }
    for (i = 0; i < numStep; i++) {
      new_displ_data = calloc(numDim*numVertex, sizeof(double));
      new_ave_face_data = malloc(numDim*numElement*sizeof(double));
      if (!new_displ_data || !new_ave_face_data)
	fc_exitIfError(FC_MEMORY_ERROR);
      for (j = 0; j < 2; j++) {
	rc = fc_getVariableDataAsDataType(displs[j][i], FC_DT_DOUBLE,
					  (void**)&displ_datas[j]);
	fc_exitIfError(rc);
      }
      if (numDeathVar > 0) {
	new_isDead_data = malloc(numElement*sizeof(int));
	if (!new_isDead_data)
	  fc_exitIfError(FC_MEMORY_ERROR);
	for (j = 0; j < 2; j++) {
	  death_datas[j] = (double**)malloc(numDeathVar*sizeof(double*));
	  if (!death_datas[j])
	    fc_exitIfError(FC_MEMORY_ERROR);
	  for (k = 0; k < numDeathVar; k++) {
	    fc_getVariableDataAsDataType(deathVars[j][k][i], FC_DT_DOUBLE,
					 (void**)&death_datas[j][k]);
	    fc_exitIfError(rc);
	  }
	}
      }
      vindex = 0;
      eindex = 0;
      for (j = 0; j < 2; j++) {
	int otherMesh = (j == 0 ? 1 : 0);
	for (k = 0; k < numKeeps[j]; k++) {
	  int thisFaceID = keepThisFaceIDs[j][k];
	  int otherFaceID = keepOtherFaceIDs[j][k];
	  // assume skin face (only 1 parent element)
	  int thisElemID = elemParentsPerFace[j][thisFaceID][0];
	  int otherElemID = elemParentsPerFace[otherMesh][otherFaceID][0];
	  FC_Coords thisFaceCoords[4], otherFaceCoords[4], newPoint;
	  FC_Vector thisNormal, otherNormal;
	  // First the vertex
	  memcpy(&new_displ_data[vindex*numDim],
		 &displ_datas[j][keepVertIDs[j][k]*numDim], 
		 numDim*sizeof(double));
	  vindex++;
	  // Then the parameterized point on the face
	  for (m = 0; m < 4; m++) {
	    int otherVertID = faceToVertConns[otherMesh][otherFaceID*4+m];
	    for (n = 0; n < numDim; n++)
	      otherFaceCoords[m][n] = coords[otherMesh][numDim*otherVertID+n]
		+ displ_datas[otherMesh][numDim*otherVertID+n];
	  }
	  _fc_calcQuadLocation(otherFaceCoords, 2, keepParams[j][k],
			       &newPoint);
	  for (m = 0; m < numDim; m++)
	    new_displ_data[vindex*numDim+m] = newPoint[m] - 
	                                      new_coords[vindex*numDim+m];
	  vindex++;
	  // Also save the face normals
	  for (m = 0; m < 4; m++) {
	    int thisVertID = faceToVertConns[j][thisFaceID*4+m];
	    for (n = 0; n < numDim; n++)
	      thisFaceCoords[m][n] = coords[j][numDim*thisVertID+n]
		+ displ_datas[j][numDim*thisVertID+n];
	  }
	  _fc_calcSurfaceNormal(4, thisFaceCoords, &thisNormal);
	  _fc_calcSurfaceNormal(4, otherFaceCoords, &otherNormal);
	  if (faceOrients[j][thisFaceID] == faceOrients[otherMesh][otherFaceID]) {
	    for (m = 0; m < numDim; m++)
	      new_ave_face_data[eindex*numDim+m] = (thisNormal[m] - otherNormal[m])/2.;
	  }
	  else {
	    for (m = 0; m < numDim; m++)
	      new_ave_face_data[eindex*numDim+m] = (thisNormal[m] + otherNormal[m])/2.;
	  }
	  // also calc & save isDead
	  if (numDeathVar > 0) {
	    new_isDead_data[eindex] = isItDead(numDeathVar, deathVarInfos,
					       thisElemID, death_datas[j],
					       otherElemID, death_datas[otherMesh]);
	    //int isDead = 1;
	    //for (m = 0; m < numDeathVar; m++) {
	    //}

	    // FIX! do a real calc!
	    //if (death_datas[j][0][thisElemID] == 1 || 
	    //death_datas[otherMesh][0][otherElemID] == 1) 
	    // new_isDead_data[eindex] = 1;
	    //else
	    //new_isDead_data[eindex] = 0;
	  }
	  eindex++;
	}
      }
      rc = fc_setVariableDataPtr(gap_displ[i], numVertex, numDim, FC_AT_VERTEX,
				 FC_MT_VECTOR, FC_DT_DOUBLE, 
				 (void*)new_displ_data);
      fc_exitIfError(rc);
      rc = fc_setVariableDataPtr(aveFaceNorms[i], numElement, numDim, 
				 FC_AT_ELEMENT, FC_MT_VECTOR, FC_DT_DOUBLE, 
			    (void*)new_ave_face_data);
      fc_exitIfError(rc);
      free(displ_datas[0]);
      free(displ_datas[1]);
      if (numDeathVar > 0) {
	rc = fc_setVariableDataPtr(isDeadVar[i], numElement, 1, FC_AT_ELEMENT,
			      FC_MT_SCALAR, FC_DT_INT, (void*)new_isDead_data);
	fc_exitIfError(rc);
	for (j = 0; j < 2; j++) {
	  for (k = 0; k < numDeathVar; k++)
	    free(death_datas[j][k]);
	  free(death_datas[j]);
	}
      }
    }
  }
  free(displs[0]);
  free(displs[1]);
  // FIX free death seq vars arrays (but not the vars)

  // Create more variables
  {
    FC_Variable templengths[numStep], tempnorms[numStep], temptangs[numStep];
    for (i = 0; i < numStep; i++) {
      FC_Variable tempvector, tempnorm, temptang;
      rc = fc_getDisplacedEdgeLengths(gap_mesh, gap_displ[i], &templengths[i]);
      fc_exitIfError(rc);
      rc = getLinesAsVectors(gap_mesh, gap_displ[i], &tempvector);
      fc_exitIfError(rc);
      rc = fc_createNormalTangentVariables2(tempvector, aveFaceNorms[i],
					    &tempnorm, &temptang);
      fc_exitIfError(rc);
      rc = fc_varUnaryFunction(tempnorm, fabs, "normal", &tempnorms[i]);
      fc_exitIfError(rc);
      rc = fc_varUnaryFunction(temptang, fabs, "tanget", &temptangs[i]);
      fc_exitIfError(rc);
      fc_deleteVariable(tempvector);
      fc_deleteVariable(tempnorm);
      fc_deleteVariable(temptang);
    }
    rc = fc_convertVariablesToSeqVariable(numStep, templengths, sequence, 
				     "length", &lengths);
    fc_exitIfError(rc);
    rc = fc_convertVariablesToSeqVariable(numStep, tempnorms, sequence, 
				     "length_norm", &norms);
    fc_exitIfError(rc);
    rc = fc_convertVariablesToSeqVariable(numStep, temptangs, sequence, 
				     "length_tang", &tangs);
    fc_exitIfError(rc);
  }

  /*
  { // for debugging - make gaplines from this and from other
    int eindex = 0;
    for (i = 0; i < 2; i++) {
      FC_Mesh temp_mesh;
      FC_Subset temp_subset;
      if (i == 0) 
	sprintf(charBuf, "thisGaplines");
      else
	sprintf(charBuf, "otherGaplines");
      fc_createSubset(gap_mesh, charBuf, FC_AT_ELEMENT, &temp_subset);
      for (j = 0; j < numKeeps[i]; j++) {
	fc_addMemberToSubset(temp_subset, eindex);
	eindex++;
      }
      // HACK for exodus not writing element subsets
      fc_createSubsetMesh(temp_subset, dataset, 1, charBuf, &temp_mesh);
    }
  }
  */

  // Create per side gaplines subset
  {
    int eindex;
    FC_SortedIntArray **sias[2];
    // create initialize lists
    for (i = 0; i < 2; i++) { 
      sias[i] = malloc(numShapes[i]*sizeof(FC_SortedIntArray*));
      if (sias[i] == NULL)
	fc_exitIfError(FC_MEMORY_ERROR);
      for (j = 0; j < numShapes[i]; j++) {
	sias[i][j] = calloc(shapes[i][j].numSides, sizeof(FC_SortedIntArray));
	if (sias[i][j] == NULL)
	  fc_exitIfError(FC_MEMORY_ERROR);
      }
    }
    // add gapline ids to appropriate lists
    eindex = 0;
    for (i = 0; i < 2; i++) {
      int otherMesh = (i == 0 ? 1 : 0);
      for (j = 0; j < numKeeps[i]; j++) {
	int shapeID, sideID;
	findShapeSide(keepThisFaceIDs[i][j], numShapes[i], shapes[i], 
		      &shapeID, &sideID);
	fc_addIntToSortedIntArray(&sias[i][shapeID][sideID], eindex);
	findShapeSide(keepOtherFaceIDs[i][j], numShapes[otherMesh],
	      shapes[otherMesh], &shapeID, &sideID);
	fc_addIntToSortedIntArray(&sias[otherMesh][shapeID][sideID], eindex);
	eindex++;
      }
    }
    // make subset of non empty lists
    numSide = 0;
    sides = NULL;
    // loop over first mesh's shapes
    for (i = 0; i < numShapes[0]; i++) {
      for (j = 0; j < shapes[0][i].numSides; j++) {

	// skip if no entries
	if (sias[0][i][j].numVal < 1)
	  continue;

	// loop over 2nd mesh's shapes
	for (k = 0; k < numShapes[1]; k++) {
	  for (m = 0; m < shapes[1][k].numSides; m++) {
	    FC_SortedIntArray temp_sia = { 0, 0, 0 };

	    // skip if no entries
	    if (sias[1][k][m].numVal < 1)
	      continue;

	    // look for intersection
	    getRemoveSIAIntersection(&sias[0][i][j], &sias[1][k][m], 
				     &temp_sia);
	    if (temp_sia.numVal > 0) {
	      FC_Subset temp_subset;
	      char *subset_name = malloc((strlen(mesh_names[0]) +
					  strlen(mesh_names[1]) + 
					  100)*sizeof(char));
	      sprintf(subset_name, "%s_shape%d_side%d-%s_shape%d_side%d",
		      mesh_names[0], i, j, mesh_names[1], k, m);
	      rc = fc_createSubset(gap_mesh, subset_name, FC_AT_ELEMENT,
			      &temp_subset);
	      fc_exitIfError(rc);
	      free(subset_name);
	      fc_addArrayMembersToSubset(temp_subset, temp_sia.numVal, 
					 temp_sia.vals);
	      fc_freeSortedIntArray(&temp_sia);
	      sides = realloc(sides, (numSide+1)*sizeof(FC_Subset));
	      if (!sides) 
		fc_exitIfError(FC_MEMORY_ERROR);
	      sides[numSide] = temp_subset;
	      numSide++;
	    }
	  }
	} // done' with 2nd mesh's shapes
      }
    } // done with 1 mesh's shapes
    // delete lists
    for (i = 0; i < 2; i++) { 
      for (j = 0; j < numShapes[i]; j++) {
	for (k = 0; k < shapes[i][j].numSides; k++) 
	  fc_freeSortedIntArray(&sias[i][j][k]);
	free(sias[i][j]);
      }
      free(sias[i]);
    }
  } 

  // make alive subsets for entire gap mesh
  alivePerTime = NULL;
  if (numDeathVar > 0) {
    //FC_Mesh temp_mesh;
    alivePerTime = (FC_Subset*)malloc(numStep*sizeof(FC_Subset));
    for (i = 0; i < numStep; i++) {
      rc = fc_createThresholdSubset(isDeadVar[i], "==", 0, "temp",
				    &alivePerTime[i]);
      fc_exitIfErrorPrintf(rc, "failed to create alive subset");
      // FIX delete these subsets
      // debug hack
      //sprintf(charBuf, "alive at step %d", i);
      //rc = fc_createSubsetMesh(alivePerTime[i], dataset, 1, charBuf, 
      //		       &temp_mesh);
      //fc_exitIfErrorPrintf(rc, "could not create mesh versin of alive subset");
    }
  }

  // Report min/max/ave/stdev lengths of lines for each step
  rc = fc_getSequenceCoordsAsDataType(sequence, FC_DT_DOUBLE,
				      (void**)&seq_coords); 
  fc_exitIfErrorPrintf(rc, "Failed to get sequence coords");

  printf("\n");
  printf("Number of gap lines found = %d\n", numKeeps[0] + numKeeps[1]);
  printf("Number of sets of sides involved = %d\n", numSide);
  for (i = 0; i < numSide; i++) {
    int numMember;
    char *temp_name;
    fc_getSubsetName(sides[i], &temp_name);
    fc_getSubsetNumMember(sides[i], &numMember);
    printf("\n");
    printf("Stats for set %d ('%s'):\n", i, temp_name);
    free(temp_name);
    printf("numGapline = %d\n", numMember);
    printStats("Gap Length", numStep, lengths, &sides[i], alivePerTime,
	       seq_coords);
    printStats("Normal Component of Gap Length", numStep, norms,
	       &sides[i], alivePerTime, seq_coords);
    printStats("Tangent Component of Gap Length", numStep, tangs,
	       &sides[i], alivePerTime, seq_coords);
  }
  if (numSide > 1) {
    printf("\n");
    printf("Overall stats:\n");
    printf("numGapline = %d\n", numKeeps[0] + numKeeps[1]);
    fflush(NULL);
    printStats("Gap Length", numStep, lengths, NULL, alivePerTime, seq_coords);
    printStats("Normal Component of Gap Length", numStep, norms, NULL,
	       alivePerTime, seq_coords);
    printStats("Tangent Component of Gap Length", numStep, tangs, NULL,
	       alivePerTime, seq_coords);
  }
  fflush(NULL);

  // delete some more stuff so that it doesn't show up in output
  for (i = 0; i < 2; i++) {
    for (j = 0; j < numShapes[i]; j++)
      fc_freeShape(&shapes[i][j]);
    free(shapes[i]);
  }
  fc_deleteSeqVariable(numStep, aveFaceNorms);

  // FIX! HACK! Since exodus can't write element subsets, create mesh 
  // copies of sides
  for (i = 0; i < numSide; i++) {
    FC_Mesh temp_mesh;
    char *temp_name;
    fc_getSubsetName(sides[i], &temp_name);
    rc = fc_createSubsetMesh(sides[i], dataset, 1, temp_name, &temp_mesh);
    free(temp_name);
    fc_exitIfErrorPrintf(rc, 
                       "Could not create mesh versions of the side subsets");
  }

  // Write the gaplines dataset
  if (verbose_level > FC_QUIET) {
    printf("Writing gaplines dataset ...\n");
    fflush(NULL);
  }
  {
    char *temp_name;
    fc_getDatasetName(dataset, &temp_name);
    sprintf(charBuf, "gap lines for '");
    strncat(charBuf, temp_name, 1000);
    free(temp_name);
    strncat(charBuf, "'", 1);
    fc_changeDatasetName(dataset, charBuf);
    fc_rewriteDataset(dataset, "gaplines.ex2", FC_FT_EXODUS);
  }

  // Final library
  fc_finalLibrary();

  exit(0);  
}
