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
 * \file ex3_createCharactrztn.c
 * \brief Example of creating a new characterization. 
 *
 * \description
 *    
 *    This annotated example presents a simple example of creating a new
 *    characterization.
 *
 *    Please see ex1_intro.c and ex2_bigDataAccess.c first for examples of 
 *    data access and using built-in FCLib characterizations.
 *
 *    NOTE: If you are looking at this file through the Doxygen interface,
 *    all objects and functions are linked to their documentation--so 
 *    follow the links for more details!
 *
 *    NOTE: Don't forget to try running the program while looking at this file.
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/examples/ex3_createCharactrztn.c,v $
 * $Revision: 1.7 $ 
 * $Date: 2006/09/26 21:31:32 $
 *
 * \modifications
 *    - 10/12/05 WSD, created. Split out of example 2 because gcc 4.0
 *      didn't like function declarations within functions.
 */

#include <string.h>
#include <fc.h>

// --- Write a characterization

// Now that we know how to access the data object's big data (see
// ex2_bigDataAcess.c, let's write our own characterization. 

// Let's say that we want to find the average temperature in each
// Feature for each time step, but the FCLib function,
// fc_getVariableSubsetMeanSdev(), does not quite give us what we want.
// Let's say we want to weight the average so that values on the boundaries
// of the feature do not count as much. One way we could do this is to
// weigh values on a entity by how many of the entity's neighbors are
// also in the feature. (I admit that this is kind of a contrived example.)

// Here we create the new characterization and later we will apply it.
static FC_ReturnCode myCharacterizationFunkyAve(FC_Variable var, 
			  FC_Subset subset, double* funky_average) 
{
  // In C, have to declare all of your variables up front
  FC_ReturnCode rc;                // holds fclib error flags
  int i, j;                        // loop variables
  FC_Mesh mesh;                    // will hold the subset's parent mesh
  int numDataPoint, numComponent;  // variable info
  FC_DataType dataType;            // more var info
  FC_AssociationType assoc_var, assoc_subset; // more var & subset info
  int maxNumMember, numMember;     // subset info 
  void *datap;                     // var data
  int *membersMask;                // array maxNumMember long where 
  // nonzero entries indicate member
  int numNeighbor, *neighborIDs;   // neighbors
  double value, sum;               // helpers to do calc
  int totalNum, tempNum;           // "
  
  // Set a default return value in case anything goes wrong
  if (funky_average)
    *funky_average = -1;
  
  // Check that the input arguments seem sane. fc_isXValid() routines
  // exist for all data types and return true if the handle refers
      // to something that really exists in the library.
  if (!fc_isVariableValid(var) || !fc_isSubsetValid(subset) ||
      !funky_average) {
    return FC_INPUT_ERROR;
  }
  
  // Get info from var & subset
  rc = fc_getVariableInfo(var, &numDataPoint, &numComponent, &assoc_var,
			  NULL,  &dataType);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_getSubsetInfo(subset, &numMember, &maxNumMember, &assoc_subset);
  if (rc != FC_SUCCESS)
    return rc;
  
  // Do a little more checking - we refuse the right to serve variables
  // and subsets that do not meet our expectations!
  if (numComponent != 1 || numDataPoint < 1)  {
    fc_printfErrorMessage("Variable can only have 1 component not %d", 
			  numComponent);
    return FC_INPUT_ERROR;
  }
  if (!fc_isDataTypeValid(dataType) || dataType == FC_DT_UNKNOWN) {
    fc_printfErrorMessage("Can not operate on data type %s", 
			  fc_getDataTypeText(dataType));
    return FC_INPUT_ERROR;
  }
  if (assoc_var != assoc_subset || numDataPoint != maxNumMember) {
    fc_printfErrorMessage("Mismatch of entities");
    return FC_INPUT_ERROR;
  }
  if (numMember < 1) {
    fc_printfErrorMessage("Subset cannot be empty");
    return FC_INPUT_ERROR;
  }
  
  // Get variable big data & the subset members (we get the subset member
  // as a mask this time because we will want to quickly lookup whether
  // an entity is in the subset).
  rc = fc_getVariableDataPtr(var, &datap);
  if (rc != FC_SUCCESS)
    return rc;
  rc = fc_getSubsetMembersAsMask(subset, &maxNumMember, &membersMask);
  if (rc != FC_SUCCESS)
    return rc;
  // We are also going to need the parent mesh.
  rc = fc_getMeshFromSubset(subset, &mesh);
  if (rc != FC_SUCCESS)
    return rc;
  
  // Basic algorithm for the calculation:
  // 1) step through each member of subset
  //    a) figure out how many neighbors of the member are also in the 
  //       subset
  //    b) decode the void* data array to get the value on the member
  //    c) increment the weighted sum and the total number of values count
  // 2) Calculate average.
  
  // Do the calculation. 
  
  // Initialize our number of values and the sum of the values to zero.
  totalNum = 0;
  sum = 0.;
  
  // The members are only the entities who's IDs in the mask are equal 
  // to 1.
  for (i = 0; i < maxNumMember; i++) {
    if (membersMask[i] > 0) {
      
      // This routine comes from the topo module and it returns the IDs of
      // the neighbors of an entity. We then count how many of these
      // neighbors are in the subset. Don't forget to free the neighborIDs!
      tempNum = 0;
      rc = fc_getMeshEntityNeighbors(mesh, i, assoc_subset, 0, 
				     &numNeighbor, &neighborIDs);
      if (rc != FC_SUCCESS)
	return rc;
      for (j = 0; j < numNeighbor; j++) {
	if (membersMask[neighborIDs[j]] > 0)
	  tempNum++;
      }
      free(neighborIDs);
      
      // Get the value of this member
      switch(dataType) {
      case FC_DT_INT:     value = ((int*)datap)[i];     break;
      case FC_DT_FLOAT:   value = ((float*)datap)[i];   break;
      case FC_DT_DOUBLE:  value = ((double*)datap)[i];  break;
      default: ; // should never reach here 
      }
      
      // Increment the weighted sum and the total number of values
      sum += tempNum * value;
      totalNum += tempNum;
    }
  }
  
  // cleanup (free allocated memory)
  free(membersMask);
  
  // Calculate the average and stuff into return variable.
  *funky_average = sum / totalNum;
  
  return FC_SUCCESS;
}; // End of myCharacterizationFunkyAve()


// --- Use the characterization

// We will apply our characterization to the same host spots we found in
// ex1_intro.c. First, we'll find the hot spots on the hex mesh
// and track them. Then we will apply myCharacterizationFunkyAve() to
// each ROI of each feature.

int main(void)
{
  FC_ReturnCode rc;
  int i, j;
  int numStep, numFeature, numMeshes, numROI, *stepIDs;
  int numReturnSeqVars, *numStepPerReturnSeqVar, numReturnSubsets;
  double average, funky_average, **funky_averages;
  void* seqCoords;
  FC_Dataset dataset;
  FC_Sequence sequence;
  FC_Mesh mesh, *meshes;
  FC_Subset subset, *returnSubsets,thresholdSubset, *ROIs;
  FC_Variable *seqVariable, **returnSeqVars;
  FC_DataType dataType;
  FC_FeatureGroup *featureGroup;
  char char_buf[1024];
  char dataset_name[1024];
  char mesh_name[1024];
  char seqVariable_name[1024];
  char sequence_name[1024];
  char subset_name[1024];
  char* temp_name;

  printf("Running ex3_createCharactrztn ...\n");
  printf("\n");

  // !! Read the documentation for myCharacterizationFunkyAve() before
  // !! continuing !!

  // --- Init library and get Data Object handles
  
  // The code to get the data object handles is taken directly from
  // ex1_intro.c which we assume the reader has read.
  // Get the hex mesh, the temperature sequence variable, and the sequence
  rc = fc_setLibraryVerbosity(FC_WARNING_MESSAGES);
  if (rc != FC_SUCCESS)
    exit(rc);
  rc = fc_initLibrary();
  if (rc != FC_SUCCESS)
    exit(rc);
  strcpy(dataset_name, "../data/gen_multivar_seq.ex2");
  rc = fc_loadDataset(dataset_name, &dataset);
  if (rc != FC_SUCCESS)
    exit(rc);
  strcpy(mesh_name, "hex mesh");
  rc = fc_getMeshByName(dataset, mesh_name, &numMeshes,&meshes);
  if (rc != FC_SUCCESS)
    exit(rc);
  if (numMeshes != 1){
    rc = FC_INPUT_ERROR;
    fc_exitIfErrorPrintf(rc, "Failed to find mesh '%s'",mesh_name);
  }
  mesh = meshes[0];
  free(meshes);
  strcpy(subset_name, "every 5th element");
  rc = fc_getSubsetByName(mesh, subset_name, &numReturnSubsets,
			  &returnSubsets);
  if (rc != FC_SUCCESS)
    exit(rc);
  if (numReturnSubsets != 1)
    exit(FC_INPUT_ERROR);
  subset = returnSubsets[0];
  free(returnSubsets);
  strcpy(seqVariable_name, "temperature per vertex");
  rc = fc_getSeqVariableByName(mesh, seqVariable_name, &numReturnSeqVars,
			       &numStepPerReturnSeqVar,&returnSeqVars);
  if (rc != FC_SUCCESS || numReturnSeqVars != 1)
    exit(rc);
  numStep = numStepPerReturnSeqVar[0];
  seqVariable = returnSeqVars[0];
  free(numStepPerReturnSeqVar);
  free(returnSeqVars);
  rc = fc_getSequenceFromSeqVariable(numStep, seqVariable, &sequence);
  if (rc != FC_SUCCESS)
    exit(rc);
  fc_getSequenceName(sequence, &temp_name);
  strcpy(sequence_name, temp_name);
  free(temp_name);
  
  // Find the hot spots on the hex mesh and track them. (This is a little 
  // different from what we did in ex1_intro.c in that we'll track at
  // each step instead of all at once.)
  printf("Finding & Tracking hot spots:\n");
  rc = fc_createFeatureGroup(&featureGroup);
  if (rc != FC_SUCCESS)
    exit(rc);
  for (i = 0; i < numStep; i++) {
    sprintf(char_buf, "Step%d_ThresholdSubset", i);
    rc = fc_createThresholdSubset(seqVariable[i], ">", 100., 
				  char_buf, &thresholdSubset);
    if (rc != FC_SUCCESS)
      exit(rc);
    rc = fc_segment(thresholdSubset, 0, &numROI, &ROIs);
    if (rc != FC_SUCCESS)
      exit(rc);
    rc = fc_trackStep(i, numROI, ROIs, featureGroup);
    if (rc != FC_SUCCESS)
      exit(rc);
    fc_deleteSubset(thresholdSubset);
    free(ROIs);
  }
  fc_printFeatureGroup(featureGroup);

  // Now let's compare the normal average characterization with our
  // new funky average. Note that the funky average is always a little
  // greater than the normal average. This is what would be expected
  // as the hottest part of the feature (for this particular dataset)
  // lies in the center of the feature and we are weighing center
  // values more strongly in the funky average.
  printf("Comparing FCLib's average with the new funky average:\n");
  rc = fc_featureGroup_getNumFeature(featureGroup, &numFeature);
  if (rc != FC_SUCCESS)
    exit(rc);
  for (i = 0; i < numFeature; i++) {
    rc = fc_getFeatureROIs(featureGroup, i, &numROI, &stepIDs, &ROIs);
    if (rc != FC_SUCCESS)
      exit(rc);
    printf("Feature #%d:\n", i);
    printf("    numROI = %d\n", numROI);
    for (j = 0; j < numROI; j++) {
      rc = fc_getVariableSubsetMeanSdev(seqVariable[stepIDs[j]], ROIs[j], 
					&average, NULL);
      if (rc != FC_SUCCESS)
	exit(rc);
      rc = myCharacterizationFunkyAve(seqVariable[stepIDs[j]], ROIs[j], 
				      &funky_average);
      if (rc != FC_SUCCESS)
	exit(rc);
      printf("    %d (Step #%d): FCLib's average  = %6g\n", j, stepIDs[j],
	     average);
      printf("                 My funky average = %6g\n", funky_average); 
    }
    free(stepIDs);
    free(ROIs);
  }
  printf("\n");
  
  // Let's save the funky averages to print out differently later.
  // We'll make room to store an average for each feature for each
  // time step, but since most features do not exist over all
  // time steps, we are going to have some holes. Because we happen
  // to know that our averages will all be positive, if we 
  // initialize the storage space with -1's we'll be able to
  // see the holes more easily later.
  funky_averages = malloc(numFeature*sizeof(double*));
  for (i = 0; i < numFeature; i++) {
    funky_averages[i] = malloc(numStep*sizeof(double));
    for (j = 0; j < numStep; j++) 
      funky_averages[i][j] = -1;
  }
  for (i = 0; i < numFeature; i++) {
    fc_getFeatureROIs(featureGroup, i, &numROI, &stepIDs, &ROIs);
    for (j = 0; j < numROI; j++) {
      myCharacterizationFunkyAve(seqVariable[stepIDs[j]], ROIs[j], 
				 &funky_averages[i][stepIDs[j]]);
    }
    free(stepIDs);
    free(ROIs);
  }
  
  // As our final example in this annotated example, we'll print out the 
  // funky averages of the features in a plotable format.
  fc_getSequenceDataType(sequence, &dataType);
  fc_getSequenceCoordsPtr(sequence, &seqCoords);
  printf("Funky Average for each feature: (-1 = no value for this step)\n");
  printf("----");
  for (i = 0; i < numFeature; i++) 
    printf("---------------");
  printf("\n");
  printf("%4s", "Time");
  for (i = 0; i < numFeature; i++) {
    sprintf(char_buf, "Feature%d", i);
    printf("%15s", char_buf);
  }
  printf("\n");
  printf("----");
  for (i = 0; i < numFeature; i++) 
    printf("---------------");
  printf("\n");
  for (i = 0; i < numStep; i++) {
    if (dataType == FC_DT_CHAR)
      printf("%4c", ((char*)seqCoords)[i]);
    else if (dataType == FC_DT_INT)
      printf("%4d", ((int*)seqCoords)[i]);
    else if (dataType == FC_DT_FLOAT)
      printf("%4.2f", ((float*)seqCoords)[i]);
    else if (dataType == FC_DT_DOUBLE)
      printf("%4.2f", ((double*)seqCoords)[i]);
    for (j = 0; j < numFeature; j++) {
      printf("%15.6g", funky_averages[j][i]);
    }
    printf("\n");
  }
  printf("\n");
  
  // Cleanup the dynamically allocated memory.
  for (i = 0; i < numFeature; i++)
    free(funky_averages[i]);
  free(funky_averages);

  // --- All done

  // Memory cleanup - Arrays of small data and data objects not managed
  // by the library.
  free(seqVariable);
  //writeout feature graph
  fc_writeFeatureGraph(featureGroup,"ex3_FeatureGraph");
  fc_freeFeatureGroup(featureGroup);

  // When you are all done it is a good idea to finalize the library.
  // This will automatically delete all datasets (and deleting a
  // dataset automatically deletes all of its sequences & meshes & etc.)
  fc_finalLibrary();

  // Bye!
  printf("... finished running ex3_createCharactrztn\n");
  exit(0);
}
