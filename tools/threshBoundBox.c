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
 * \file threshBoundBox.c
 * \brief Generate bounding boxes in a dataset based on thresholding
 *        a variable.
 *
 * $Source: /home/Repositories/fcdmf/fclib/tools/threshBoundBox.c,v $
 * $Revision: 1.39 $ 
 * $Date: 2006/11/07 23:49:03 $
 *
 * \description
 *    Usage: threshBoundBox dataset "variable name" "op" threshold_value
 *    where "op: is the threshold operator such as ">", "<=", "==", etc. 
 *    Quotes are required around the op string, but can be dropped around
 *    the variable name if there are no spaces or punctuation in the name.
 *
 *    Example: threshBoundBox data.ex2 "cool kid" > 0.5
 *
 *    Prints to stdout the bounding boxes of regions in the dataset
 *    for which the given variable is above the given threshold value.
 *    The quotes around the variable name can be dropped if it
 *    contains no spaces or punctuation.
 *    printf("The format of the output is one line per bounding box with form
 *
 *    t xmin ymin zmin xmax ymax zmax\n");
 *
 *    where t is the time step index, (xmin ymin zmin) is the corner
 *    of the bounding box with the smallest coordinate values, and
 *    xmax ymax zmax) is the corner with the largest coordinates 
 *    of the bounding box. If there are no time steps, t = -1.
 *
 *    Note: Widths of the bounding box can be zero. For example, if the
 *    bounding box contains a single point.
 *
 * \modifications
 *   - 02/05/04 WSK, created. 
 */

#include <string.h>
#include "fc.h"

// the core of the code
// if displs == NULL don't do displacements
int find_and_print_bbs(int stepID, int num, FC_Variable* displs, 
		       FC_Variable* vars, char* condition,
                       double reference, int doCollapse);
 
int main(int argc, char** argv) {
  FC_ReturnCode rc;
  int i, j;
  int numMesh, numBasicVar, numSeqVar, numStep, temp_numStep;
  FC_Dataset ds;
  FC_Mesh* meshes, *seq_meshes;
  FC_Variable* vars;      // array of vars (numMesh long)
  FC_Variable* displs;    // array of displacement vars (numMesh long)
  FC_Variable** seq_vars; // array of seq vars (numMesh long)
  FC_Variable** seq_displs; // array of displacement seq vars (numMesh long)
  double threshold;
  char* file_name = NULL;
  char* displ_var_name = NULL;
  char* var_name;
  char* op;
  FC_VerbosityLevel verbose_level = FC_QUIET; 
  int doSmoothing = 0;   // flag for doing smoothing
  int doCollapse = 0;    // flag for collapsing bounding boxes
  double radius;
  int doLocalSmooth = 0;
  FC_Variable* smooth_vars; // array of vars (numMesh long)
  FC_AssociationType assoc;
  int numComp;
  int doDispl = 0;
  
  
  // handle arguments
  if (argc < 5) {
  usage:
    printf("usage: %s [options] dataset \"var name\" \"op\" thresh_value \n", 
           argv[0]);
    printf("options: \n");
    printf("   -c              : combine any overlapping bounding boxes into a\n");
    printf("                     single bounding box.\n");
    printf("   -d var_name     : calc the displaced values using this displacement\n");
    printf("                     variable\n");
    printf("   -r radius       : smooth the variable by averaging each data\n");
    printf("                     point with all points within a radius r of it.\n");
    printf("   --doLocalSmooth : if true, only points on local mesh are used for smoothing\n"); 
    printf("                     if false, points on all meshes are used for smoothing \n");
    printf("   -h              : print this help message\n");
    printf("   -v              : verbose: print warning and error messages\n");
    printf("   -V              : very verbose: print log and error messages\n");
    printf("\n");
    printf("Prints to stdout the bounding boxes of regions in the dataset\n");
    printf("for which the given variable satisfies the threshold operation.\n");
    printf("\"op\" is an operator such as \">\", \"<=\", \"==\", etc.\n");
    printf("Quotes are required around the op string, but can be dropped\n");
    printf("dropped around the variable name if it contains no spaces or\n");
    printf("punctuation. \n");
    printf("\n");
    printf("Example: %s data.ex2 \"cool factor\" > 0.5\n", argv[0]);
    printf("\n");
    printf("The format of the output is one line per bounding box with form\n");
    printf("\n");
    printf("  t xmin ymin zmin xmax ymax zmax\n");
    printf("\n");
    printf("where t is the time step index, (xmin ymin zmin) is the corner\n");
    printf("of the bounding box with the smallest coordinate values, and\n");
    printf("(xmax ymax zmax) is the corner with the largest coordinate \n");
    printf("values. If there are no time steps, t = -1.\n");
    fflush(NULL);
    exit(-1);
  }
  for (i = 1; i < argc; i++) {
    if (!strcmp(argv[i], "-c")) {
      doCollapse = 1;
    }
    else if (!strcmp(argv[i], "-d")) {
      doDispl = 1;
      i++;
      displ_var_name = argv[i];
    }
    else if (!strcmp(argv[i], "-r")) {
      if (i+1 >= argc)
        goto usage;
      doSmoothing = 1;
      radius = atof(argv[i+1]);
      i++;
    }
    else if (!strcmp(argv[i], "--doLocalSmooth") || !strcmp(argv[i], "-doLocalSmooth")) {
      doLocalSmooth = 1;
    }
    else if (!strcmp(argv[i], "-v")) {
      verbose_level = FC_WARNING_MESSAGES;
    }
    else if (!strcmp(argv[i], "-V")) {
      verbose_level = FC_LOG_MESSAGES;
    }
    else if (!strncmp(argv[i], "-h", 2) || !strncmp(argv[i], "-", 1))
      goto usage;
    else {
      if (i+3 >= argc)
        goto usage;
      file_name = argv[i];
      var_name = argv[i+1];
      op = argv[i+2];
      threshold = atof(argv[i+3]);
      i+=3;
    }
  }
  if (!file_name)
    goto usage;
  
  // init library and load dataset
  rc = fc_setLibraryVerbosity(verbose_level);
  fc_exitIfError(rc);
  rc = fc_initLibrary();
  fc_exitIfError(rc);
  rc = fc_loadDataset(file_name, &ds);
  fc_exitIfError(rc);

  // get number of meshes and setup var arrays
  rc = fc_getMeshes(ds, &numMesh, &meshes);
  fc_exitIfError(rc);
  seq_meshes = (FC_Mesh*)malloc(numMesh*sizeof(FC_Mesh));
  vars = (FC_Variable*)malloc(numMesh*sizeof(FC_Variable));
  seq_vars = (FC_Variable**)malloc(numMesh*sizeof(FC_Variable*));

  // loop over meshes
  for (i = 0; i < numMesh; i++) {
    // keep another copy for seq vars
    seq_meshes[i] = meshes[i];

    // try to access as basic var & seq var,
    rc = fc_getOrGenerateUniqueVariableByName(meshes[i], var_name, &vars[i]);
    fc_exitIfError(rc);

    rc = fc_getOrGenerateUniqueSeqVariableByName(meshes[i], var_name, &numStep,
						 &seq_vars[i]);
    fc_exitIfError(rc);
  }

  // Remove mesh/var and mesh/seq_var entries that don't exist
  numBasicVar = 0;
  numSeqVar = 0;
  for (i = 0; i < numMesh; i++) {
    if (FC_HANDLE_EQUIV(vars[i], FC_NULL_VARIABLE)) {
      for (j = i + 1; j < numMesh; j++) {
        meshes[j-1] = meshes[j];
        vars[j-1] = vars[j]; 
      }
    }
    else 
      numBasicVar++;
    if (seq_vars[i] == NULL) {
      for (j = i + 1; j < numMesh; j++) {
        seq_meshes[j-1] = seq_meshes[j];
        seq_vars[j-1] = seq_vars[j];
      }
    }
    else 
      numSeqVar++;
  }
  if (numBasicVar == 0) {
    free(meshes);
    free(vars);
  }
  if (numSeqVar == 0) {
    free(seq_meshes);
    free(seq_vars);
  }
  
  // did we find anything?
  if (numBasicVar < 1 && numSeqVar < 1) 
    fc_exitIfErrorPrintf(FC_ERROR, "Failed to find variable '%s' in dataset "
                         "'%s'", var_name, file_name);

  // make sure variable makes sense
  if (numBasicVar > 0) {
    for (i = 0; i < numBasicVar; i++) {
      fc_getVariableInfo(vars[i], NULL, &numComp, &assoc, NULL, NULL);
      if (numComp != 1)
        fc_exitIfErrorPrintf(FC_ERROR, "Cannot threshold variable '%s' because"
                             " it has more than one component\n", var_name);
    }
  }
  // make sure seq variable makes sense
  if (numSeqVar > 0) {
    FC_Sequence sequence, temp_sequence;
    fc_getSequenceFromSeqVariable(numStep, seq_vars[0], &sequence);
    fc_exitIfError(rc);
    for (i = 0; i < numSeqVar; i++) {
      fc_getSequenceFromSeqVariable(numStep, seq_vars[i], &temp_sequence);
      if (!FC_HANDLE_EQUIV(temp_sequence, sequence)) 
        fc_exitIfErrorPrintf(FC_ERROR, "Abort because variable '%s' was found "
                             "on more than one sequence", var_name);
      fc_getVariableInfo(seq_vars[i][0], NULL, &numComp, &assoc, NULL, NULL);
      if (numComp != 1)
        fc_exitIfErrorPrintf(FC_ERROR, "Cannot threshold variable '%s' because"
                             " it has more than one component\n", var_name);
    }
  }

  // if displacement vars requested, make sure we have them
  displs = NULL;
  if (doDispl) {
    if (numBasicVar > 0) {
      displs = (FC_Variable*)malloc(numBasicVar*sizeof(FC_Variable));
      if (!displs)
	fc_exitIfError(FC_MEMORY_ERROR);
      for (i = 0; i < numBasicVar; i++) {
	rc = fc_getOrGenerateUniqueVariableByName(meshes[i], displ_var_name,
						  &displs[i]);
	fc_exitIfError(rc);
	if (FC_HANDLE_EQUIV(displs[i], FC_NULL_VARIABLE)){
	  fc_exitIfErrorPrintf(FC_INPUT_ERROR,
			       "Cannot find unique displacement variable '%s'",
			       displ_var_name);
	}
      }
    }
    if (numSeqVar > 0) {
      seq_displs = (FC_Variable**)malloc(numSeqVar*sizeof(FC_Variable*));
      if (!seq_displs)
	fc_exitIfError(FC_MEMORY_ERROR);
      for (i = 0; i < numSeqVar; i++) {

	rc = fc_getOrGenerateUniqueSeqVariableByName(seq_meshes[i], displ_var_name, 
						     &temp_numStep, &seq_displs[i]);
	if (!seq_displs[i]){
	  fc_exitIfErrorPrintf(FC_INPUT_ERROR, 
			       "Cannot find (unique) displacement variable '%s'",
			       displ_var_name);
	}
	
	if (temp_numStep != numStep) 
	  fc_exitIfErrorPrintf(FC_ERROR, "displacement variable and variable "
			       "to be smoothed must have same numStep");
      }
    }
  }

  // do the basic vars first
  if (numBasicVar > 0) {
    FC_VertexBin* bins[numBasicVar];

    // setup vars to use
    if (!doSmoothing) {
      smooth_vars = vars;
    }
    else { // (doSmoothing) 
      // create bins
      for (i = 0; i < numBasicVar; i++) {
	if (doDispl)
	  rc = fc_createDisplacedMeshVertexBin(meshes[i], displs[i], &bins[i]);
	else
	  rc = fc_createMeshVertexBin(meshes[i], &bins[i]);
        fc_exitIfError(rc);
      }
      // do smooth
      if (doDispl)
	rc = fc_displacedGeomSmoothVariable(numBasicVar, meshes, displs, bins,
					    vars, radius, doLocalSmooth, &smooth_vars);
      else
	rc = fc_geomSmoothVariable(numBasicVar, meshes, bins, vars, radius, 
				   doLocalSmooth, &smooth_vars);
      fc_exitIfError(rc);
      // smoothing cleanup
      for (i = 0; i < numBasicVar; i++)
        fc_freeVertexBin(bins[i]);
    }
    
    // do it
    rc = find_and_print_bbs(-1, numBasicVar, displs, smooth_vars, op, threshold,
                            doCollapse); 
    fc_exitIfError(rc);
    
    // cleanup
    if (doSmoothing) {
      for (i = 0; i < numBasicVar; i++)
	fc_deleteVariable(smooth_vars[i]);
      free(smooth_vars);
    }
    if (doDispl) {
      for (i = 0; i < numBasicVar; i++) 
	fc_deleteVariable(displs[i]);
      free(displs);
    }
    for (i = 0; i < numBasicVar; i++)
      fc_deleteVariable(vars[i]);
    free(vars);
    free(meshes); // don't close meshes, still in use via seq_meshes
  }

  // then do the seq vars by step
  if (numSeqVar > 0) {
    FC_VertexBin* bins[numSeqVar];
    for (i = 0; i < numStep; i++) {
      FC_Variable temp_vars[numSeqVar];
      FC_Variable* temp_displs, temp_displs_buf[numSeqVar];

      // get the vars for just this step
      for (j = 0; j < numSeqVar; j++) 
        temp_vars[j] = seq_vars[j][i];
      temp_displs = NULL;
      if (doDispl) {
	for (j = 0; j < numSeqVar; j++) 
	  temp_displs_buf[j] = seq_displs[j][i];
	temp_displs = temp_displs_buf;
      }

      if (!doSmoothing) {
	smooth_vars = temp_vars;
      }
      else { // (doSmoothing)
	// create bins
	if (doDispl) {
	  for (j = 0; j < numSeqVar; j++) {
	    rc = fc_createDisplacedMeshVertexBin(seq_meshes[j], temp_displs[j], 
						 &bins[j]);
	    fc_exitIfErrorPrintf(rc, "failed to create displ mesh vertex bin");
	  }
	}
	else if (i == 0) { // only have to do first time
	  for (j = 0; j < numSeqVar; j++) {
	    rc = fc_createMeshVertexBin(seq_meshes[j], &bins[j]);
	    fc_exitIfErrorPrintf(rc, "failed to create mesh vertex bin");
	  }
	}
	// do smooth
	if (doDispl) 
	  rc = fc_displacedGeomSmoothVariable(numSeqVar, seq_meshes,
					      temp_displs, bins, temp_vars,
					      radius, doLocalSmooth, &smooth_vars);
	else
	  rc = fc_geomSmoothVariable(numSeqVar, seq_meshes, bins, temp_vars, 
				     radius, doLocalSmooth, &smooth_vars); 
	fc_exitIfError(rc);
	// smoothing cleanup (for nonDispl case, only on last step)
	if (doDispl || i == numStep -1) { 
	  for (j = 0; j < numSeqVar; j++)
	    fc_freeVertexBin(bins[j]);
	}
      }
      
      // do it
      rc = find_and_print_bbs(i, numSeqVar, temp_displs, smooth_vars, op, 
			      threshold, doCollapse);
      fc_exitIfError(rc);

      // cleanup
      if (doSmoothing) {
	for (j = 0; j < numSeqVar; j++) 
	  fc_deleteVariable(smooth_vars[j]);
	free(smooth_vars);
      }
    }
    if (doDispl) {
      for (i = 0; i < numSeqVar; i++) {
	fc_deleteSeqVariable(numStep, seq_displs[i]);
	free(seq_displs[i]);
      }
      free(seq_displs);
    }
    for (i = 0; i < numSeqVar; i++) {
      fc_deleteSeqVariable(numStep, seq_vars[i]);
      free(seq_vars[i]);
    }
    free(seq_vars);
    free(seq_meshes);
  }
  
  // all done
  fc_finalLibrary();
  
  exit(0);
}

// Find the thresholded regions, combine bounding boxes if requested,
// & print output.
int find_and_print_bbs(int stepID, int num, FC_Variable* displs, 
		       FC_Variable* vars, char* condition, 
                       double threshold, int doCollapse) {
  FC_ReturnCode rc; 
  int i, j, k, m;
  int numDim, numBox, numSegment;
  FC_Coords *lowpoints, *highpoints;
  FC_Coords temp_lowpoint, temp_highpoint;
  FC_Coords combined_lowpoint, combined_highpoint;
  FC_Subset* segments;        // segments
  int used_temp, overlap_flag;
  FC_Subset subset;
  char *subsetName = "temp_subset";
  int doDispl = 0;

  // set doDispl flag
  if (displs)
    doDispl = 1;

  // get the bounding boxes of the thresholded regions
  numBox = 0;
  lowpoints = highpoints = NULL;
  for (i = 0; i < num; i++) { 
    // get segments
    rc = fc_createThresholdSubset(vars[i], condition, threshold, subsetName, 
                                  &subset);                             
    if (rc != FC_SUCCESS)
      return rc;
    rc = fc_segment(subset, 0, &numSegment, &segments);
    if (rc != FC_SUCCESS)
      return rc;
    fc_deleteSubset(subset);

    // If none, move on
    if (numSegment < 1)
      continue;

    // Extend bb array
    lowpoints = realloc(lowpoints, (numBox + numSegment)*sizeof(FC_Coords));
    highpoints = realloc(highpoints, (numBox + numSegment)*sizeof(FC_Coords));
    if (!lowpoints || !highpoints) 
      fc_exitIfError(FC_MEMORY_ERROR);

    for (j = 0; j < numSegment; j++) {
      // get bb
      if (doDispl) 
	rc = fc_getDisplacedSubsetBoundingBox(segments[j], displs[i],
					      &numDim, &temp_lowpoint,
					      &temp_highpoint);
      else
	rc = fc_getSubsetBoundingBox(segments[j], &numDim, &temp_lowpoint, 
				     &temp_highpoint);
      if (rc != FC_SUCCESS)
        return rc;
      fc_deleteSubset(segments[j]);

      // potentially combine with a previous bb (from a prev vars)
      used_temp = 0;
      if (doCollapse) {
        for (k = 0; k < numBox; k++) {
          fc_getBoundingBoxesOverlap(numDim, lowpoints[k], highpoints[k],
                                     temp_lowpoint, temp_highpoint, 
                                     &overlap_flag, NULL, NULL);
          if (overlap_flag > 0) {
            fc_combineBoundingBoxes(numDim, lowpoints[k], highpoints[k],
                                    temp_lowpoint, temp_highpoint,
                                    &combined_lowpoint, &combined_highpoint);
            for (m = 0; m < 3; m++) {
              lowpoints[k][m] = combined_lowpoint[m];
              highpoints[k][m] = combined_highpoint[m];
            }
            used_temp = 1;
            break;
          }
        }
      }

      // if we didn't combine, add to bb list
      if (!used_temp) {
        for (k = 0; k < 3; k++) {
          lowpoints[numBox][k] = temp_lowpoint[k];
          highpoints[numBox][k] = temp_highpoint[k];
        }
        numBox++;
      }
    }
    free(segments);
  }
  
  // print the bounding boxes - have to wait in case we combine
  for (i = 0; i < numBox; i++) {
    printf("%d %f %f %f %f %f %f\n", stepID, 
           lowpoints[i][0], lowpoints[i][1], lowpoints[i][2], 
           highpoints[i][0], highpoints[i][1], highpoints[i][2]);
    fflush(NULL);
  }

  // cleanup
  free(lowpoints);
  free(highpoints);
  
  return FC_SUCCESS;
}
