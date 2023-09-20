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
 * \file ex1_intro.c
 * \brief Example of basic use of the library.
 *
 * \description
 *    
 *    This is an annotated example that walks through the use of basic aspects
 *    of the Feature Characterization Library (FCLib).  It covers basic data
 *    access (getting handles and querying meta data) and a simple example of
 *    applying feature tracking and FCLib characterizations.
 *
 *    See ex2_bigDataAccess.c for details of accessing the big data arrays and
 *    ex3_createCharactrztn.c for details about creating your own 
 *    characterization.
 *
 *    NOTE: If you are looking at this file through the Doxygen interface,
 *    all objects and functions are linked to their documentation--so 
 *    follow the links for more details!
 *
 *    NOTE: Don't forget to try running the program while looking at this file.
 *
 * $Source: /home/Repositories/fcdmf/fclib/examples/ex1_intro.c,v $
 * $Revision: 1.5 $ 
 * $Date: 2006/09/26 03:57:16 $
 *
 * \modifications
 *    - 11/18/04 WSK, created.
 */

#include <string.h>
#include <fc.h>

int main(void)
{
  FC_ReturnCode rc;
  int i, j;
  int numSequence, numMesh, numSubset, numVariable, numSeqVariable;
  int numReturnSeqVariables, *numStepPerReturnSeqVar;
  int numCompareMesh;
  int numStep, topodim, dim, numVertex, numElement, *numStepPerSeqVar;
  int numMember, maxNumMember, numDataPoint, numComponent;
  int numSegment, *numROIPerStep, numFeature, numROI, *stepIDs;
  int numParent, *parentIDs, numChild, *childIDs;
  double threshold, maxTemperature;
  FC_Coords centroid;
  FC_Dataset dataset;
  FC_Sequence *sequences, sequence;
  FC_Mesh *meshes, mesh, *comparemeshes;
  FC_Subset *subsets, thresholdSubset, *segments, **ROIsPerStep, *ROIs;
  FC_Variable *variables, **seqVariables, *seqVariable, **returnSeqVariables;
  FC_DataType dataType;
  FC_ElementType elemType;
  FC_AssociationType assocType;
  FC_MathType mathType;
  FC_FeatureGroup *featureGroup;
  char char_buf[1024];
  char dataset_name[1024];
  char mesh_name[1024];
  char seqVariable_name[1024];
  char *temp_name;

  printf("Running ex1_intro ...\n");
  printf("\n");

  // --- Setting up the library

  // The very first step is to initialize the library. 
  // The only thing you can do before that is to set the library verbosity
  // which instructs all library functions on how much to print to stdout.
  // Here we've chosen to have the library report error and warning messages.
  rc = fc_setLibraryVerbosity(FC_WARNING_MESSAGES);
  if (rc != FC_SUCCESS)
    exit(rc);
  rc = fc_initLibrary();
  if (rc != FC_SUCCESS)
    exit(rc);

  // Most function calls in FCLib return an FC_ReturnCode object (which is an 
  // enumerated type). You should always check the return code to make sure
  // that the function completed as expected and to handle errors
  // appropriately. An error may indicated something horrible has gone wrong
  // (memory allocation error) or it may indicate that you passed incompatible
  // arguments to a function.  

  // This is an example of what you can expect to see if an error message is
  // generated: Calling fc_initLibrary() more than once is an error, so we
  // should see an error message ...
  printf("About to call fc_initLibrary() a second time ...\n");
  printf("(Message format = 'MessageType:function[file:line#]:message')\n");
  rc = fc_initLibrary(); // this generates the message
  // ... and it should also return an error FC_ErrorCode.
  printf("... the type of error you get is: %s\n", 
         fc_getReturnCodeText(rc));
  printf("\n");

  // FC_ReturnCode and many other small data objects (e.g. FC_ElementType,
  // FC_DataType, etc) in FCLib are enumerated types. All enumerated types have
  // a function of the form fc_get<Type>TypeText() that returns a string
  // that prints the name of the enumerated value.

  // --- A note about "returned" values.

  // Technically in C, functions return a single value, which for most cases
  // in the library will be an FC_ReturnCode value. This means that most
  // of the "returned" values from FCLib functions are not returned at
  // all, but are values assigned to arguments passed into the functions.
  // For the remainder of the documentation, the used of "returned" will
  // usually mean a value returned via the arguments of the function.
  // In the Doxygen documentation, these arguments are labeled as "output"
  // arguments.

  // --- Playing with Datasets, Meshes & Variables (& Sequences & Subsets)

  // The major data objects of the library are Datasets, Sequences, Meshes,
  // Subsets and Variables. These objects are managed by the library and are
  // manipulated by the user via opaque handles (FC_Dataset, FC_Mesh, etc).
  // The data objects are organized hierarchically with Datasets containing
  // Sequences and Meshes, and Meshes containing Subsets and Variables.
  // In general, data objects are added to and deleted from the library by 
  // create and delete calls. A delete call will remove the object, and all
  // of the objects it contains, from the library.
 
  // The Dataset object is a container for the other objects and is usually
  // associated with a database file. Our next step is to load a dataset by
  // providing the file name. Here is an example of a value "returned" by
  // the functions--we pass in a pointer to where we want the dataset handle
  // to be stored.
  strcpy(dataset_name, "../data/gen_multivar_seq.ex2");
  rc = fc_loadDataset(dataset_name, &dataset);
  if (rc != FC_SUCCESS)
    exit(rc);

  // Loading a dataset creates all of the data objects that it contains for
  // you, but you will still need to get their handles. There are a variety
  // of ways to do this. Here we ask the dataset for all of its sequences
  // and meshes (next we will query their metadata).
  rc = fc_getSequences(dataset, &numSequence, &sequences);
  if (rc != FC_SUCCESS)
    exit(rc);
  rc = fc_getMeshes(dataset, &numMesh, &meshes);
  if (rc != FC_SUCCESS)
    exit(rc);
  printf("Dataset '%s' has %d sequences and %d meshes\n", dataset_name, 
         numSequence, numMesh);
  printf("\n");

  // A sequence is a container for coordinates orthogonal to the space
  // coordinates of the dataset. The most common type of sequence is
  // a time series. Here we query each sequence for name & metadata.
  printf("Sequences:\n");
  for (i = 0; i < numSequence; i++) {
    fc_getSequenceName(sequences[i], &temp_name);
    fc_getSequenceInfo(sequences[i], &numStep, &dataType);
    printf("Sequence #%d: name = '%s'\n", i, temp_name);
    printf("              numStep = %d, dataType = %s\n", numStep,
           fc_getDataTypeText(dataType));
    // make sure you free the sequence name when you are done with it
    free(temp_name);
  }
  printf("\n");

  // Again we see that the values that are "returned" by library calls are 
  // actually function arguments that pass a pointer to the thing to be
  // returned; the second and third arguments of fc_getSequenceInfo() are
  // output arguments. If one of the output arguments is an array, as is 
  // temp_name above, then the user is usually responsible for freeing the
  // array which was  dynamically allocated. This is documented for each 
  // function.

  // As with sequences, each mesh can be queried for name & meta data
  printf("Meshes:\n");
  for (i = 0; i < numMesh; i++) {
    fc_getMeshName(meshes[i], &temp_name);
    fc_getMeshInfo(meshes[i], &topodim, &dim, &numVertex, &numElement, 
                   &elemType);
    printf("Mesh #%d: name = '%s'\n", i, temp_name);
    printf("          topodim = %d, dim = %d\n", topodim, dim);
    printf("          numVertex = %d, numElement = %d, elemType = %s\n", 
           numVertex, numElement, fc_getElementTypeText(elemType)); 
    // make sure you free the mesh name when you are done with it
    free(temp_name);
  }
  printf("\n");

  // Alternately, if you known the name of the mesh of interest, you can
  // ask for that directly from the dataset. This will return all meshes
  // of the same name. In this case we check to make sure there is only
  // one mesh of the asked for name.
  strcpy(mesh_name, "hex mesh");
  rc = fc_getMeshByName(dataset, mesh_name, &numCompareMesh,&comparemeshes);
  if (rc != FC_SUCCESS)
    exit(rc);
  if (numCompareMesh != 1){
    printf("Did not find a unique mesh '%s'\n",mesh_name);
    exit(FC_ERROR);
  }
  mesh = comparemeshes[0];
  free(comparemeshes);

  // This is how you could figure out which of the meshes gotten from the
  // dataset using fc_getMeshes() is the same as the one you just got
  // by asking for it by name. Because the data objects are only handles
  // (that the library uses to look up the data), directly comparing the
  // handles will only tell you if you are looking at the same handle--
  // you need to use FC_HANDLE_EQUIV() to ask the library if the
  // two handles are referring to the same object. 
   for (i = 0; i < numMesh; i++) {
    if (FC_HANDLE_EQUIV(mesh, meshes[i]))
      printf("Mesh #%d is the same mesh as '%s'\n", i, mesh_name);
    else
      printf("Mesh #%d is NOT the same mesh as '%s'\n", i, mesh_name);
  }
  printf("\n");

  // Datasets contain sequences and meshes. Meshes in turn contain subsets
  // and variables.
  rc = fc_getSubsets(mesh, &numSubset, &subsets);
  if (rc != FC_SUCCESS)
    exit(rc);
  rc = fc_getVariables(mesh, &numVariable, &variables);
  if (rc != FC_SUCCESS)
    exit(rc);
  rc = fc_getSeqVariables(mesh, &numSeqVariable, &numStepPerSeqVar,
			  &seqVariables);
  if (rc != FC_SUCCESS)
    exit(rc);

  // A subset is a collection of mesh entities (mesh entities are vertices,
  // edges, faces, or elements of the mesh). FC_AssociationType is an
  // enumerated type that specifies the type of mesh entity. (maxNumMember
  // is equal to the number of that type of entity in the parent mesh).
  printf("Subsets in mesh '%s':\n", mesh_name);
  for (i = 0; i < numSubset; i++) {
    fc_getSubsetName(subsets[i], &temp_name);
    fc_getSubsetInfo(subsets[i], &numMember, &maxNumMember, &assocType);
    printf("Subset #%d: name = '%s'\n", i, temp_name);
    printf("            numMember = %d, maxNumMember = %d\n", numMember, 
           maxNumMember);
    printf("            assocType = %s\n", 
           fc_getAssociationTypeText(assocType)); 
    // make sure you free the subset name when you are done with it
    free(temp_name);
  }
  printf("\n");

  // Variables are the field data that is stored on a mesh (e.g. a temperature
  // field). In FCLib, variables come in two flavors based on whether they
  // vary over a sequence or not. In all cases the basic unit of operation
  // is the variable handle, FC_Variable. A variable that does not vary
  // over a sequence has a single variable handle, while a variable that
  // does vary over a sequence is an array of variable handles. This second
  // type of variables is distinguished from the first by calling it a
  // sequence variable. 
  printf("Variables in mesh '%s':\n", mesh_name);
  for (i = 0; i < numVariable; i++) {
    fc_getVariableName(variables[i], &temp_name);
    fc_getVariableInfo(variables[i], &numDataPoint, &numComponent, &assocType,
                       &mathType, &dataType);
    printf("Variable #%d: name = '%s'\n", i, temp_name);
    printf("              numDataPoint = %d, numComponent = %d\n", 
           numDataPoint, numComponent);
    printf("              assocType = %s\n", 
           fc_getAssociationTypeText(assocType)); 
    // make sure you free the variable name when you are done with it
    free(temp_name);
  }
  printf("\n");

  // A sequence variable is the pair of an int giving the length of the 
  // handle array and the handle array itself. Functions specifically
  // for sequence variables need both of these. Many other functions
  // you will have to apply to sequence variables manually by looping
  // over the steps and applying them--as we will see later.
  // Here we only bother to get info on the first step since the metadata
  // should be the same for all steps of the sequence variable.
  printf("Sequence Variables in mesh '%s':\n", mesh_name);
  for (i = 0; i < numSeqVariable; i++) {
    fc_getVariableName(seqVariables[i][0], &temp_name);
    fc_getVariableInfo(seqVariables[i][0], &numDataPoint, &numComponent, 
                       &assocType, &mathType, &dataType);
    printf("SeqVariable #%d: name = '%s'\n", i, temp_name);
    printf("                 numStep = %d\n", numStepPerSeqVar[i]);
    printf("                 numDataPoint = %d, numComponent = %d\n", 
           numDataPoint, numComponent);
    printf("                 assocType = %s\n", 
           fc_getAssociationTypeText(assocType)); 
    // make sure you free the seqVariable name when you are done with it
    free(temp_name);
  }
  printf("\n");

  // Again, if we know the name of the sequence variable we are interested
  // in (or any of the other major data objects) we can ask for it by name.
  strcpy(seqVariable_name, "temperature per vertex");
  rc = fc_getSeqVariableByName(mesh, seqVariable_name, &numReturnSeqVariables,
			       &numStepPerReturnSeqVar, &returnSeqVariables);
  if (rc != FC_SUCCESS)
    exit(rc);
  if (numReturnSeqVariables != 1){
    fc_exitIfErrorPrintf(FC_INPUT_ERROR, "Can't find (unique) variable '%s'",
			  seqVariable_name);
  }
  seqVariable = returnSeqVariables[0];
  numStep = numStepPerReturnSeqVar[0];
  free(returnSeqVariables);
  free(numStepPerReturnSeqVar);
  
  // Usually there is only one sequence in a dataset, but there can be
  // more and a mesh may have sequence variables that are on different
  // sequences. You can get the handle of the sequence that a sequence
  // variable is on by using fc_getSequenceFromSeqVariable() as we do
  // here for the first sequence variable.
  rc = fc_getSequenceFromSeqVariable(numStep, seqVariable, &sequence);
  if (rc != FC_SUCCESS)
    exit(rc);

  // Before moving on, let's cleanup up some allocated memory that we
  // no longer want. We are only deleting the arrays of handles--the
  // subset and variables still exist and we can get them again.
  free(sequences);
  free(meshes);
  free(subsets);
  free(variables);
  free(numStepPerSeqVar);
  for (i = 0; i < numSeqVariable; i++)
    free(seqVariables[i]);
  free(seqVariables);

  // --- Performing a characterization

  // Let's do a characterization by examining hot spots in the temperature
  // sequence variable on the hex mesh. First we'll find the hot spots in each
  // time step. Then we'll track them so we know their behavior over time.
  // Finally we'll characterize them by doing some basic statistics
  // and also finding their centroids.

  // Let's define a hot spot as a region with a temperature greater than 100.
  // To find these regions in the first step we first threshold the first step
  // of our sequence variable to create a subset that has ALL of the
  // mesh entities that have values greater than 100. Then we separate
  // these entities into separate connected components. The second argument
  // of fc_segment is a flag to tell the segmentation algorithm how to
  // judge connectedness.
  printf("Finding hot spots:\n");
  threshold = 100.;
  rc = fc_createThresholdSubset(seqVariable[0], ">", 100., 
                                "Step0_ThresholdSubset", &thresholdSubset);
  if (rc != FC_SUCCESS)
    exit(rc);
  rc = fc_segment(thresholdSubset, 0, &numSegment, &segments);
  if (rc != FC_SUCCESS)
    exit(rc);
  printf("Found %d hot spots in step#0 of sequence variable '%s'\n",
         numSegment, seqVariable_name);
  
  // We created an intermediate subset that we won't need again so let's
  // delete it.  Note that further calls with this handle will report
  // errors.
  rc = fc_deleteSubset(thresholdSubset);
  if (rc != FC_SUCCESS)
    exit(rc);

  // The segment (a single connected component) will from now on be
  // called an ROI (region of interest). Technically, an ROI is the subset 
  // handle plus the stepID. Let's store the current ROI and find others in the
  // rest of the seqVariable's steps.
  numROIPerStep = (int*)malloc(numStep*sizeof(int));
  ROIsPerStep = (FC_Subset**)malloc(numStep*sizeof(FC_Subset*));
  numROIPerStep[0] = numSegment;
  ROIsPerStep[0] = segments;
  for (i = 1; i < numStep; i++) {
    sprintf(char_buf, "Step%d_ThresholdSubset", i);
    rc = fc_createThresholdSubset(seqVariable[i], ">", 100., 
                                  char_buf, &thresholdSubset);
    if (rc != FC_SUCCESS)
      exit(rc);
    rc = fc_segment(thresholdSubset, 0, &numSegment, &segments);
    if (rc != FC_SUCCESS)
      exit(rc);
    printf("Found %d hot spots in step#%d of sequence variable '%s'\n",
           numSegment, i, seqVariable_name);
    rc = fc_deleteSubset(thresholdSubset);
    if (rc != FC_SUCCESS)
      exit(rc);
    numROIPerStep[i] = numSegment;
    ROIsPerStep[i] = segments;
  }
  printf("\n");
  
  // Although we've found all of the ROIs, there is still some
  // important information missing: how do we know which ROIs from
  // different time steps are really the same hot spot? ROI#0 from step#1
  // is not necessary the same hot spot as ROI#0 from step#2. Also, at step#5
  // we suddenly only have 1 hot spot when there used to be two. Did one hot
  // spot disappear? Did the two hot spots merge into one? Or is it a 
  // completely new hot spot? To figure this out we introduce feature tracking.
  // Features are structures that exist over time and are therefore essentially
  // a series of ROIs. In FCLib, features are stored in an object called
  // a feature group. fc_trackAllSteps() is used to create the feature
  // group from the ROIs. 
  rc = fc_createFeatureGroup(&featureGroup);
  if (rc != FC_SUCCESS)
    exit(rc);
  rc = fc_trackAllSteps(numStep, numROIPerStep, ROIsPerStep, featureGroup);
  if (rc != FC_SUCCESS)
    exit(rc);

  // Intemediate cleanup: All of the ROIs are now organized within the feature
  // group so we don't need our old organization of the data any more.
  // Again, we are only cleaning up handles, the ROIs still exists.
  free(numROIPerStep);
  for (i = 0; i < numStep; i++)
    free(ROIsPerStep[i]);
  free(ROIsPerStep);

  // Now you can ask the group about features
  rc = fc_featureGroup_getNumFeature(featureGroup, &numFeature);
  if (rc != FC_SUCCESS)
    exit(rc);
  printf("Found %d features:\n", numFeature);

  // Now report the characterizations for each feature - a feature is
  // accessed by the feature Group and the feature ID.
  for (i = 0; i < numFeature; i++) {
    rc = fc_getFeatureROIs(featureGroup, i, &numROI, &stepIDs, &ROIs);
    if (rc != FC_SUCCESS)
      exit(rc);
    rc = fc_getFeatureParentIDs(featureGroup, i, &numParent, &parentIDs);
    if (rc != FC_SUCCESS)
      exit(rc);
    rc = fc_getFeatureChildIDs(featureGroup, i, &numChild, &childIDs);
    if (rc != FC_SUCCESS)
      exit(rc);
    printf("Feature #%d:\n", i);
    printf("    numParent = %d, IDs = ", numParent);
    if (numParent > 0) {
      printf("%d", parentIDs[0]);
      for (j = 1; j < numParent; j++)
        printf(", %d", parentIDs[1]);
    }
    printf("\n");
    printf("    numChild  = %d, IDs = ", numChild);
    if (numChild > 0) {
      printf("%d", childIDs[0]);
      for (j = 1; j < numChild; j++)
        printf(", %d", childIDs[1]);
    }
    printf("\n");
    printf("    numStep = %d\n", numROI);

    // Do the characterizations, these are per step (ROI) of the feature
    //   1) get ROI maximum temperature
    //   2) get ROI centroid
    for (j = 0; j < numROI; j++) {
      rc = fc_getVariableSubsetMinMax(seqVariable[stepIDs[j]], ROIs[j], NULL,
                                      NULL, &maxTemperature, NULL);
      if (rc != FC_SUCCESS)
        exit(rc);
      rc = fc_getSubsetCentroid(ROIs[j], &dim, &centroid);
      if (rc != FC_SUCCESS)
        exit(rc);
      printf("    Step#%d: maxTemp = %f\n", stepIDs[j], maxTemperature);
      printf("             centroid = (%f, %f, %f)\n", centroid[0], 
             centroid[1], centroid[2]);
    }

    // cleanup
    free(stepIDs);
    free(ROIs);
    free(parentIDs);
    free(childIDs);
  }
  printf("\n");

  // cleanup
  fc_freeFeatureGroup(featureGroup);
  free(seqVariable);

  // --- All done

  // When you are all done it is a good idea to finalize the library.
  // This will automatically delete all datasets (and deleting a
  // dataset automatically deletes all of its sequences & meshes & etc.)
  fc_finalLibrary();

  // Bye!
  printf("... finished running ex1_intro\n");
  exit(0);
}
