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
 * \file  analyzeSpotWelds.c
 * \brief Analysis of the spot weld behavior in the results of a Presto 
 *        Simulation
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/tools/analyzeSpotWelds.c,v $
 * $Revision: 1.53 $ 
 * $Date: 2006/09/28 06:34:39 $
 *
 * \description 
 *
 *   This program generates a file with information about how many spot welds
 *   failed and when. This program also generates a parameter, Pfail, for each
 *   spot weld. Spotwelds fail when Pfail reaches one. Usage and validation
 *   documentation is available from FCLib developers.
 *
 *
 * \todo ?: more flexible coding of using sequences - what if more than 1?
 *       what if none?
 *  
 * \modifications
 *   - APR-08-2003  W Koegler  Created
 *   - 2003-JUN-25  W Koegler  Use new Sierra input file parsing code
 *       to remove hard coding dependencies
 *   - 2003-JUN-30  W Koegler  Have input filenames an arguments instead
 *       of hard coding.  Removed hard coding of sequence name.
 *   - 2003-JUL-10  W Koegler  Added Pfail calculations when taking into
 *       account coordinate displacements
 *   - 2003-JUL-10  W Koegler  Moved "broken" classifying into FCLib's custom.c
 *   - 2003-AUG-18  W Koegler  Changed code so that groups spot welds by
 *       'type' which means they share the same normal and tangent functions.
 *   - 2003-SEP-03  W Koegler  Added ability to use spot weld scale factors
 *       if they are present (just in case). Spot welds that use the same
 *       normal and tangential functions, but with different scale factors
 *       are of different types.
 *   - 2003-SEP-08  W Koegler  Did a lot of cleaning up and added analysis
 *       based on distance between the node and side sets of a spot weld.
 *       The force based analysis is really an approximation of this.
 *   - 2003-NOV-05  WSK  
 *       - Bug fix: was only printing 3 significant digits to the .dat file,
 *         now prints 6.
 *       - Added code to general a failure parameter for a group (average)
 *       - Fine tuning of .gnu file for pretty plots
 *   - 01/28/04 WSK, Changed from float to double
 *   - 06/22/04 WSK, nodesets/sidesets are now subsets instead of meshes.
 *   - 12/13/04 WSK, removed old force version since displacement version
 *       is technically more correct (force version was a hack).
 */

#include <stdarg.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <fc.h>

/** Converts a void* array to a double* array, given the length and the 
    FC_DataType. */
double* get_doubles(int num, FC_DataType dataType, void* data);

/** Holds information about a group a spot welds. */
typedef struct {
  char name[1024];
  // type identity depends on following 
  int normalId;
  int tangentId;
  char* normal_name;
  char* tangent_name;
  int decay_cycles;
  double normal_scale;
  double tangent_scale;
  double exponent;      /**< failure exponent */ 
 // other things to know about this type (gleaned from member functions)
  double ucritnorm;     // maximum displacement along normal 
  double ucrittang;     // maximum displacement along tangent 
  // the spotwelds
  int numSpotWeld;
  FC_SierraSpotWeldInfo** spotwelds;
  char** spotweld_names;
  char** nodeset_names;
  char** sideset_names;
  FC_Subset* nodesets;
  FC_Subset* sidesets;
  int numStep;
  FC_Variable** node_displs; // nodeset coord displacements per time step
  FC_Variable** side_displs; // sideset coord displacements per time step
  // analysis results from analysis
  double** fails;  /**< failure parameters per screw per step */
  double* group_fails; // normalized failure parameter over group, per step
  int* fail_stepIds; // Per spot weld, the step it failed at (-1 = didn't fail)
  int numBrokenTotal; // total number of broken spot welds
  int* numBrokens;   // number of broken spot welds per step
  double* fail_times; // per spot weld, estimate failure time if 
                      // fail_stepId > -1, else the max fail value seen
} SpotWeldType;

// generate names for spotWeldTypes, more document. with function definition
void nameSpotWeldTypes(int numType, SpotWeldType* types);

int main(int argc, char **argv)
{
  int i, j, k;  
  FC_ReturnCode rc;     // int to catch return codes
  time_t start_time = time(NULL); // the time analysis starts
  char* inputfile_name = NULL; // the simulation's input file
  char* database_name = NULL; // the simulation's results file
  char* output_root_name = NULL; // root name for output files
  char* outputfile_name = NULL;  // name for data file (.dat)
  char* gnuplotfile_name = NULL; // name of gnuplot script file (.gnu)
  char* scriptfile_name = NULL;  // name of script file that will launch
                                 // gnuplot and image magic (.script)
  char** plotfile_names = NULL;  // names for the plot files (.eps)
  char** jpgfile_names = NULL;   // names for jpg plot files (.jpg)
  char* temp_string;
  FC_VerbosityLevel beVerbose = FC_QUIET;   // verbosity flag
  FILE* outfile = NULL;
  FC_SierraInfo *sierra_info;
  char* displ_name;
  int numSpotWeld;
  int numType;
  SpotWeldType* types;
  FC_Dataset dataset;
  FC_Sequence *sequences, sequence;
  FC_Mesh* meshes;
  int numStep, numMesh;
  double* times;
  int numBrokenTotal; // value of analysis that was done
  int *numBrokens;  // pointer to the analysis that was done
  int count; 
  int time_index;
  int time_col; // time column (indexing from 1)
  int col;  // current column (indexing from 1)

  //*****************************************************************
  //*** End of variable declarations
  //*****************************************************************

  //*****************************************************************
  //*** Handle arguments
  //*****************************************************************

  if (argc < 7) { 
    usage:
      printf("\n");
      printf("summary:\n");
      printf("   analyzeSpotWelds performs an an analysis of Presto spotwelds (grouped\n"); 
      printf("   by spot weld type) and produces a data output file <output_root>.dat,\n");
      printf("   a gnuplot file <output.gnu>, and script file <output.script> that\n"); 
      printf("   can be used to create a series of eps (if gnuplot is available)\n");
      printf("   and jpg (if ImageMagick if available) plots <output>.*.[eps|jpg].\n");
      printf("\n");
      printf("usage:\n");
      printf("   analyzeSpotWelds -i <sierra_input_file> -s <database> -o <output_root>\n");
      printf("\noptions:\n");
      printf("   -v                         verbose output\n");
      printf("   -V                         very verbose output\n");
      printf("\n");
      printf("\n");
      printf("example:\n");
      printf("   analyzeSpotWelds -i problem.i -s problem.ex2 -o spotwelds\n");
      printf("\n");
      fflush(NULL);
      exit(-1);
  }
  for (i = 1; i < argc; i++) {
    if (!strcmp(argv[i], "-i") && i+1 < argc) {
      i++;
      inputfile_name = argv[i];
    }
    else if(!strcmp(argv[i], "-o") && i+1 < argc) {
      i++;
      output_root_name = argv[i];
    }
    else if(!strcmp(argv[i], "-s") && i+1 < argc) {
      i++;
      database_name = argv[i];
    }
    else if(!strcmp(argv[i], "-v")) {
      beVerbose = FC_WARNING_MESSAGES;
    }
    else if(!strcmp(argv[i], "-V")) {
      beVerbose = FC_LOG_MESSAGES;
    }
    else
      goto usage;
  }
  // check that we got everything we need
  if ( (inputfile_name == NULL) || (database_name == NULL) || 
       (output_root_name == NULL) )
    goto usage;

  //*** Setup

  if (beVerbose) {
    printf("Parsing input file...\n");
    fflush(NULL);
  }

  // parse the input file
  rc = fc_parseSierraInputFile(inputfile_name, &sierra_info);
  fc_exitIfErrorPrintf(rc, "Could not open input file '%s'", inputfile_name);

  // number of spot welds from sierra info
  numSpotWeld = sierra_info->numSpotWeld;
  if (numSpotWeld < 1)
    fc_exitIfErrorPrintf(FC_ERROR, "No spot welds were found");

  // displacement variable's names from sierra info
  displ_name = malloc(sizeof(char)*(strlen(sierra_info->displ_alias)+7));
  if (displ_name == NULL)
    fc_exitIfErrorPrintf(FC_MEMORY_ERROR, "Memory allocation failure");
  sprintf(displ_name, "%s", sierra_info->displ_alias);

  // distribute spotwelds into types
  numType = 0;
  types = NULL;
  for (i = 0; i < numSpotWeld; i++) {
    void* temp;
    FC_SierraSpotWeldInfo* sw = sierra_info->spotwelds[i];

    // find a matching type
    for (j = 0; j < numType; j++) {
      if (types[j].normalId == sw->normalId &&
          types[j].tangentId == sw->tangentId &&
          types[j].decay_cycles == sw->decay_cycles &&
          FC_FLT_EQUIV(types[j].normal_scale, sw->normal_scale) &&
          FC_FLT_EQUIV(types[j].tangent_scale, sw->tangent_scale) &&
          FC_FLT_EQUIV(types[j].exponent, sw->exponent))
        break;
    }
    // or create a new type
    if (j >= numType) {
      temp = realloc(types, sizeof(SpotWeldType)*(numType+1));
      if (temp == NULL)
        fc_exitIfErrorPrintf(FC_MEMORY_ERROR, "Could not allocate space for types");
      types = temp;
      types[j].normalId = sw->normalId;
      types[j].tangentId = sw->tangentId;
      types[j].normal_name = sierra_info->functions[sw->normalId]->name;
      types[j].tangent_name = sierra_info->functions[sw->tangentId]->name;
      types[j].decay_cycles = sw->decay_cycles;
      types[j].normal_scale = sw->normal_scale;
      types[j].tangent_scale = sw->tangent_scale;
      types[j].exponent = sw->exponent;
      types[j].ucritnorm = sierra_info->functions[sw->normalId]->abscissa_max;
      types[j].ucrittang = sierra_info->functions[sw->tangentId]->abscissa_max;
      types[j].numSpotWeld = 0;
      types[j].spotwelds = NULL;
      numType++;
    }

    // add spotweld
    temp = realloc(types[j].spotwelds, 
                   sizeof(FC_SierraSpotWeldInfo*)*(types[j].numSpotWeld+1));
    if (temp == NULL)
      fc_exitIfErrorPrintf(FC_MEMORY_ERROR, "Memory allocation failure");
    types[j].spotwelds = temp;
    types[j].spotwelds[types[j].numSpotWeld] = sw;
    types[j].numSpotWeld++;
  }

  // generate names for the spot weld types
  nameSpotWeldTypes(numType, types);

  // for each type, setup spot weld names, node and side set names
  for (i = 0; i < numType; i++) {
    SpotWeldType *type = &types[i];

    type->spotweld_names = (char**)malloc(sizeof(char*)*numSpotWeld);
    type->nodeset_names = (char**)malloc(sizeof(char*)*numSpotWeld);
    type->sideset_names = (char**)malloc(sizeof(char*)*numSpotWeld);
    if (type->spotweld_names == NULL || type->nodeset_names == NULL || 
        type->sideset_names == NULL)
      fc_exitIfErrorPrintf(FC_MEMORY_ERROR, "Memory allocation failure");

    for (j = 0; j < type->numSpotWeld; j++) {
      FC_SierraSpotWeldInfo* sw = type->spotwelds[j];
      type->spotweld_names[j] = sw->name;
      type->nodeset_names[j] = sw->nodeset_name;
      type->sideset_names[j] = sw->surface_name;
    }
  }

  //*****************************************************************
  //*** read dataset
  //*****************************************************************

  if (beVerbose) { 
    printf("Reading dataset...\n");
    fflush(NULL);
  }

  // init the library
  rc = fc_setLibraryVerbosity(beVerbose);
  fc_exitIfError(rc);  
  rc = fc_initLibrary();
  fc_exitIfErrorPrintf(rc, "Failed to initialize library");

  // load the dataset
  rc = fc_loadDataset(database_name, &dataset);
  fc_exitIfErrorPrintf(rc, "Could not load database '%s'", database_name);

  // get the sequence 
  {
    int numSeq;
    rc = fc_getSequences(dataset, &numSeq, &sequences);
    if (numSeq < 1) {
      fc_exitIfErrorPrintf(FC_ERROR, "Expected time series data");
    }
    else if (numSeq > 1) {
      fprintf(stderr, "Warning: Expected 1 time series, not %d: using "
              "first one\n", numSeq);
      fflush(NULL);
    }
    sequence = sequences[0];
    free(sequences);
  }

  // get sequence info & data
  {
    FC_DataType seq_dataType;
    void* seq_coords = NULL;

    rc = fc_getSequenceInfo(sequence, &numStep, &seq_dataType);
    fc_exitIfErrorPrintf(rc, "Failed to get sequence info");
    rc = fc_getSequenceCoordsPtr(sequence, &seq_coords); 
    fc_exitIfErrorPrintf(rc, "failed to get sequence coords");
    times = get_doubles(numStep, seq_dataType, seq_coords);
    if (times == NULL)
      fc_exitIfErrorPrintf(FC_MEMORY_ERROR, "Memory allocation failure");
  }

  // get the meshes
  {
    rc = fc_getMeshes(dataset, &numMesh, &meshes);
    if (rc != FC_SUCCESS)
      fc_exitIfErrorPrintf(rc, "Failed to get meshes");
    if (numMesh < 1)
      fc_exitIfErrorPrintf(FC_ERROR, "Could not find any meshes in dataset");
  }

  // for each spotweld type, get the node & side sets
  for (i = 0; i < numType; i++) {
    SpotWeldType *type = &types[i];

    type->nodesets = malloc(sizeof(FC_Mesh)*type->numSpotWeld);
    type->sidesets = malloc(sizeof(FC_Mesh)*type->numSpotWeld);
    if (type->nodesets == NULL || type->sidesets == NULL)
      fc_exitIfErrorPrintf(FC_MEMORY_ERROR, "Memory allocation failure");

    // FIX? can this be prettied up?
    // try getting subsets on each mesh, error if not exactly 1
    for (j = 0; j < type->numSpotWeld; j++) {
      int numNodeSet = 0, numSideSet = 0;
      int num_tempNodesets, num_tempSidesets;
      FC_Subset *temp_nodesets, *temp_sidesets;
      for (k = 0; k < numMesh; k++) {
        fc_getSubsetByName(meshes[k], type->nodeset_names[j], &num_tempNodesets,
			   &temp_nodesets); 
        fc_getSubsetByName(meshes[k], type->sideset_names[j], &num_tempSidesets,
			   &temp_sidesets);
        if (num_tempNodesets == 1){
          type->nodesets[j] = temp_nodesets[0];
          numNodeSet++;
        }
        if (num_tempSidesets == 1){
          type->sidesets[j] = temp_sidesets[0];
          numSideSet++;
        }
	if (temp_nodesets) free(temp_nodesets);
	if (temp_sidesets) free(temp_sidesets);
      }
      if (numNodeSet != 1) 
        fc_exitIfErrorPrintf(FC_ERROR, "problems finding nodeset %s",
                             type->nodeset_names[j]);
      if (numSideSet != 1) 
        fc_exitIfErrorPrintf(FC_ERROR, "problems finding sideset %s",
                             type->sideset_names[j]);
    }
  }
  free(meshes);
  
  // get variables: 
  for (i = 0; i < numType; i++) {
    SpotWeldType *type = &types[i];
    type->numStep = numStep;

    type->side_displs = malloc(sizeof(FC_Variable*) * type->numSpotWeld);
    if (type->side_displs == NULL)
      fc_exitIfErrorPrintf(FC_MEMORY_ERROR, "Memory allocation failure");
    type->node_displs = malloc(sizeof(FC_Variable*) * type->numSpotWeld);
    if (type->node_displs == NULL)
      fc_exitIfErrorPrintf(FC_MEMORY_ERROR, "Memory allocation failure");

    for (j = 0; j < type->numSpotWeld; j++) {
      int temp_numStep;
      FC_Mesh nodeset_mesh, sideset_mesh;
      
      // get parent mesh
      rc = fc_getMeshFromSubset(type->nodesets[j], &nodeset_mesh);
      fc_exitIfErrorPrintf(rc, "Failed to get nodeset's parent mesh");
      rc = fc_getMeshFromSubset(type->sidesets[j], &sideset_mesh);
      fc_exitIfErrorPrintf(rc, "Failed to get sideset's parent mesh");
      
      // get displs on the node sets
      rc = fc_getOrGenerateUniqueSeqVariableByName(nodeset_mesh,
						    displ_name,
						    &temp_numStep,
						    &(type->node_displs[j]));
      fc_exitIfErrorPrintf(rc, "Failed to find get displ variable on nodeset by name");
      if (!(type->node_displs[j])){
	fc_exitIfErrorPrintf(FC_INPUT_ERROR, "Failed to find get displ variable on nodeset by name");
      }
      
      if (temp_numStep != numStep)
        fc_exitIfErrorPrintf(FC_ERROR, "Mismatch between number of steps for"
                             " the sequence (%d) and the variable (%d)\n", 
                             numStep, temp_numStep);

      // get displs on the side sets
      rc = fc_getOrGenerateUniqueSeqVariableByName(sideset_mesh,
						   displ_name,
						   &temp_numStep,
						   &(type->side_displs[j]));
      fc_exitIfErrorPrintf(rc, "Failed to get displ variable on sideset by name");
      if (!(type->side_displs[j])){
	fc_exitIfErrorPrintf(FC_INPUT_ERROR, "Failed to find get displ variable on sideset by name");
      }
      if (temp_numStep != numStep)
        fc_exitIfErrorPrintf(FC_ERROR, "Mismatch between number of steps for "
                             "the sequence (%d) and the variable (%d)\n", 
                             numStep, temp_numStep);
    }
  }
  free(displ_name);

  //*****************************************************************
  //*** Do the analysis - calculated fails and fail_stepIds
  //*****************************************************************

  if (beVerbose) {
    printf("Performing spot weld analysis...\n");
    fflush(NULL);
  }

  // remember values for all spot welds
  numBrokenTotal = 0;
  numBrokens = calloc(numStep, sizeof(int));
  if (numBrokens == NULL)
    fc_exitIfErrorPrintf(FC_MEMORY_ERROR, "Memory allocation failure");
  
  // Do the analysis - get parameters and spot weld failures
  for (i = 0; i < numType; i++) {
    SpotWeldType* type = &types[i];
    
    // make room
    type->numBrokenTotal = 0;
    type->fails = malloc(sizeof(double*)*type->numSpotWeld);
    type->group_fails = calloc(numStep, sizeof(double)); // init to 0
    type->fail_stepIds = malloc(sizeof(int)*type->numSpotWeld);
    type->numBrokens = calloc(numStep, sizeof(int));
    if (type->fails == NULL || type->group_fails == NULL || 
        type->fail_stepIds == NULL || type->numBrokens == NULL)
      fc_exitIfErrorPrintf(FC_MEMORY_ERROR, "Memory allocation failure");
    type->fail_times = malloc(sizeof(double)*type->numSpotWeld);
    if (type->fail_times == NULL)
      fc_exitIfErrorPrintf(FC_MEMORY_ERROR, "Memory allocation failure");

    for (j = 0; j < type->numSpotWeld; j++) {
      int id;
      rc = fc_evalPrestoSpotWeld(type->spotwelds[j], type->nodesets[j],
                              type->sidesets[j], numStep, type->node_displs[j],
                              type->side_displs[j], &type->fails[j], 
                              &type->fail_stepIds[j]);
      // FIX? do we want this to shut things down?
      fc_exitIfErrorPrintf(rc, "Presto spot weld evaluation failed");
     
      // estimate failure time using linear interpolation
      id = type->fail_stepIds[j];
      if (id == 0)
        type->fail_times[j] = 0.;
      else if (id > 0) {
        double fails_diff = type->fails[j][id] - type->fails[j][id-1];
        double time_diff = times[id] - times[id-1];
        type->fail_times[j] = times[id-1] + 
          (1. - type->fails[j][id-1])*time_diff/fails_diff;
      }
      else { // (id < 0) 
        type->fail_times[j] = type->fails[j][0];
        for (k = 1; k < numStep; k++)
          if (type->fails[j][k] > type->fail_times[j])
            type->fail_times[j] = type->fails[j][k];
      }

      // accumulate broken info for spot weld type
      if (type->fail_stepIds[j] > -1) {
        type->numBrokenTotal++;
        for (k = type->fail_stepIds[j]; k < numStep; k++)
          type->numBrokens[k]++;
      }

      // accumulate sum of fails for group_fail
      for (k = 0; k < numStep; k++) {
        if (type->fails[j][k] < 1)
          type->group_fails[k] += type->fails[j][k];
        else
          type->group_fails[k] += 1.0;
      }
    } // end of loop over spotwelds

    // normalize group_fail
    for (j = 0; j < numStep; j++)
      type->group_fails[j] /= type->numSpotWeld;

    // accumulate broken info for all spot welds
    numBrokenTotal += type->numBrokenTotal;
    for (j = 0; j < numStep; j++)
      numBrokens[j] += type->numBrokens[j];
  }
   
  //*****************************************************************
  //*** Write analysis results
  //*****************************************************************

  if (beVerbose) {
    printf("Writing analysis results file...\n");
    fflush(NULL);
  }

  // open output file
  outputfile_name = malloc(sizeof(char)*(strlen(output_root_name)+5));
  sprintf(outputfile_name, "%s.dat", output_root_name);
  outfile = fopen(outputfile_name, "w");
  if (outfile == NULL)
    fc_exitIfErrorPrintf(FC_FILE_IO_ERROR, "Could not open file '%s' to write "
                   "results", outputfile_name);

  // Write header 
  fprintf(outfile, "# Data file generated by FCLib's analyzeSpotWelds\n");
  fprintf(outfile, "# %s", ctime(&start_time));
  fprintf(outfile, "#\n");
  fprintf(outfile, "# Input file: %s\n", inputfile_name);
  fprintf(outfile, "# Results file: %s\n", database_name);
  fprintf(outfile, "#\n");
  fflush(NULL);

  // Summary information
  fprintf(outfile, "# Summary:\n");
  fprintf(outfile, "# %d of %d total spot welds failed\n", numBrokenTotal, 
          numSpotWeld);
  count = 0;
  for (i = 0; i < numType; i++)
    if (types[i].numBrokenTotal > 0)
      count++;
  fprintf(outfile, "# %d of %d types of spot weld had failures\n", count, numType);
  fprintf(outfile, "# %d time steps (t = %g - %g)\n", numStep,
          times[0], times[numStep-1]);
  fprintf(outfile, "#\n");
  fflush(NULL);

  // NumBroken reported for each spot weld type
  fprintf(outfile, "# Breakdown by spot weld type:\n");
  fflush(NULL);
  for (i = 0; i < numType; i++) {
    SpotWeldType *type = &types[i];
    time_index = 0;
    while (time_index < numStep && type->numBrokens[time_index] <= 0)
      time_index ++;
    fprintf(outfile, "# %d:  %d of %d '%s' spot welds failed", i, 
            type->numBrokenTotal, type->numSpotWeld, type->name);
    if (type->numBrokenTotal > 0)
      fprintf(outfile, " beginning at time index %d (t = %g)\n", time_index, 
              times[time_index]);
    else
      fprintf(outfile, "\n");
    fflush(NULL);
  }
  fprintf(outfile, "#\n");
  fflush(NULL);

  // Write break times 
  fprintf(outfile, "# The estimated time of spot weld failure (NA = didn't fail)\n");
  fflush(NULL);
  for (i = 0; i < numType; i++) {
    SpotWeldType* type = &types[i];
    fprintf(outfile, "# type %d: '%s'\n", i, type->name);
    for (j = 0; j < type->numSpotWeld; j++) {
      fprintf(outfile, "#   %-13s ", type->spotweld_names[j]);
      if (type->fail_stepIds[j] < 0) 
        fprintf(outfile, "NA\n");
      else 
        fprintf(outfile, "%g\n", type->fail_times[j]);
    }
   fflush(NULL);
  }
  fprintf(outfile, "#\n");
  fflush(NULL);

  // Write plotable info
  fprintf(outfile, "# Data:\n");
  fprintf(outfile, "#   (Listed as time index, time value, total cumulative broken spot welds,\n");
  fprintf(outfile, "#   and then for each spot weld type: cumulative broken, group failure\n");
  fprintf(outfile, "#   parameter, and the failure parameters for each of the spot welds in\n");
  fprintf(outfile, "#   this type\n");
  fprintf(outfile, "#\n");
  // header line
  fprintf(outfile, "%-6s %11s %9s", "#Index", "Time", "totBroken");
  for (i = 0; i < numType; i++) {
    fprintf(outfile, " %9s %11s", "numBroken", "group_fail");
    for (j = 0; j < types[i].numSpotWeld; j++) 
      fprintf(outfile, " %11s", types[i].spotweld_names[j]);
  }
  fprintf(outfile, "\n");
  fflush(NULL);
  // data lines
  for (i = 0; i < numStep; i++) {
    fprintf(outfile, "%-6d %11.6g %9d", i, times[i], numBrokens[i]);
    for (k = 0; k < numType; k++) {
      fprintf(outfile, " %9d %11.6g", types[k].numBrokens[i], types[k].group_fails[i]);
      for (j = 0; j < types[k].numSpotWeld; j++)
        fprintf(outfile, " %11.6g", types[k].fails[j][i]);
    }
    fprintf(outfile, "\n");
    fflush(NULL);
  }

  // close file
  fclose(outfile);
  
  //*****************************************************************
  //*** Write gnuplot script file
  //*****************************************************************

  if (beVerbose) {
    printf("Writing gnuplot script file...\n");
    fflush(NULL);
  }

  // create gnuplot name and plot names 
  {
    int root_size = strlen(output_root_name) + 1;
    gnuplotfile_name = malloc(sizeof(char)*(root_size + 4));
    if (gnuplotfile_name == NULL)
      fc_exitIfErrorPrintf(FC_MEMORY_ERROR, "Could not allocate gnuplotfile_name");
    sprintf(gnuplotfile_name, "%s.gnu", output_root_name);
    
    plotfile_names = malloc(sizeof(char*)*(numType+1));
    jpgfile_names = malloc(sizeof(char*)*(numType+1));
    if (plotfile_names == NULL || jpgfile_names == NULL)
      fc_exitIfErrorPrintf(FC_MEMORY_ERROR, "Could not allocate plotfile_names");
    plotfile_names[0] = malloc(sizeof(char)*(root_size + 8));
    jpgfile_names[0] = malloc(sizeof(char)*(root_size + 8));
    if (plotfile_names[0] == NULL || jpgfile_names[0] == NULL)
      fc_exitIfErrorPrintf(FC_MEMORY_ERROR, "could not allocate plot or jpg file_names");
    sprintf(plotfile_names[0], "%s.all.eps", output_root_name);
    sprintf(jpgfile_names[0], "%s.all.jpg", output_root_name);
    root_size += (int)log10(numType) + 1;  // maximum chars need for numType
    for (i = 0; i < numType; i++) {
      plotfile_names[i+1] = malloc(sizeof(char)*(root_size + 5));
      jpgfile_names[i+1] = malloc(sizeof(char)*(root_size + 5));
      if (plotfile_names[i+1] == NULL || jpgfile_names[i+1] == NULL)
        fc_exitIfErrorPrintf(FC_MEMORY_ERROR, "Could not allocate plot or jpg file_names[%d]", i+1);
      sprintf(plotfile_names[i+1], "%s.%d.eps", output_root_name, i);
      sprintf(jpgfile_names[i+1], "%s.%d.jpg", output_root_name, i);
    }
  }

  // open output file
  outfile = fopen(gnuplotfile_name, "w");
  if (outfile == NULL)
    fc_exitIfErrorPrintf(FC_FILE_IO_ERROR, "Could not open file '%s' to write"
                         " results", gnuplotfile_name);
    
  // Write header 
  fprintf(outfile, "# Gnuplot script generated by FCLib's analyzeSpotWelds\n");
  fprintf(outfile, "# %s", ctime(&start_time));
  fprintf(outfile, "#\n");
  fprintf(outfile, "# References data file '%s'\n", outputfile_name);
  fprintf(outfile, "# To run: 'gnuplot %s' will generate '%s.<id>.eps'\n", 
          gnuplotfile_name, output_root_name);
  fprintf(outfile, "#   where <id> is 'all' or the index of the spot weld type\n");
  fprintf(outfile, "# To create jpgs: 'convert -density 144 %s.<id>.eps "
          "%s.<id>.jpg'\n", output_root_name, output_root_name);
  fprintf(outfile, "\n");
  fflush(NULL);

  // General setup
  fprintf(outfile, "# set up for all plots\n");
  fprintf(outfile, "set terminal postscript eps color\n");
  fprintf(outfile, "set xlabel 'Time'\n");
  fprintf(outfile, "set y2label 'Number of Broken Spot Welds'\n");
  fprintf(outfile, "set ytics nomirror\n");
  fprintf(outfile, "set y2tics\n");
  fprintf(outfile, "# uncomment below to move legend to left side of plots\n");
  fprintf(outfile, "#set key left Left reverse\n");
  fprintf(outfile, "\n");
  fflush(NULL);

  // As overall summary, plot group failure parameters and number of broken
  time_col = 2;
  col = time_col + 1;
  fprintf(outfile, "# plot all data\n");
  fprintf(outfile, "set output '%s'\n", plotfile_names[0]);
  fprintf(outfile, "set title \"Presto Spot Weld Analysis for All "
          "Spot Weld Groups\\n"
          "Total Number of Broken Spot Welds = %d (of %d)\"\n", 
          numBrokenTotal, numSpotWeld);
  fprintf(outfile, "set ylabel 'Normalized Group Failure Parameter (no units)'\n");
  if (numBrokenTotal == 0) 
    fprintf(outfile, "set y2range [0:1]\n");
  else
    fprintf(outfile, "set y2range [0:%d]\n", numBrokenTotal);
  fprintf(outfile, "plot [] [0:1] '%s' using %d:%d axes x1y2 title "
          "'NumBroken (on secondary y axis)' with linespoints", 
          outputfile_name, time_col, col);
  col++;
  for (i = 0; i < numType; i++) {
    SpotWeldType *type = &types[i];

    col++;
    fprintf(outfile, ", \\\n");
    fprintf(outfile, "'%s' using %d:%d title '%s' with linespoints", 
            outputfile_name, time_col, col, type->name);
    col++;
    col += types[i].numSpotWeld;
  }
  fprintf(outfile, "\n");
  fflush(NULL);

  // plot for each type of spot weld, failure parameters and broken
  col = time_col + 2; 
  for (i = 0; i < numType; i++) {
    SpotWeldType *type = &types[i];

    fprintf(outfile, "\n");
    fprintf(outfile, "# plot data for spot weld type index #%d (of %d)\n",
            i, numType);
    fprintf(outfile, "set output '%s'\n", plotfile_names[i+1]);
    fprintf(outfile, "set title \"Presto Spot Weld Analysis for Group '%s'\\n"
            "Number of Broken Spot Welds = %d (of %d)\"\n", type->name, 
            type->numBrokenTotal, type->numSpotWeld);
    fprintf(outfile, "set ylabel 'Failure Parameter (no units)'\n");
    if (type->numBrokenTotal == 0) 
      fprintf(outfile, "set y2range [0:1]\n");
    else
      fprintf(outfile, "set y2range [0:%d]\n", type->numBrokenTotal);
    fprintf(outfile, "plot [] [0:1] '%s' using %d:%d axes x1y2 title "
            "'NumBroken (on secondary y axis)' with linespoints", 
            outputfile_name, time_col, col);
    col += 2;
    for (j = 0; j < types[i].numSpotWeld; j++) {
      fprintf(outfile, ", \\\n");
      fprintf(outfile, "'%s' using %d:%d title '%s' with linespoints", 
              outputfile_name, time_col, col, type->spotweld_names[j]);
      col++;
    }
    fprintf(outfile, "\n");
    fflush(NULL);
  }

  // close file
  fclose(outfile);

  //*****************************************************************
  //*** Write script file
  //*****************************************************************

  if (beVerbose) {
    printf("Writing image creation script file ...\n");
    fflush(NULL);
  }

  // open script file
  scriptfile_name = malloc(sizeof(char)*(strlen(output_root_name)+8));
  sprintf(scriptfile_name, "%s.script", output_root_name);
  outfile = fopen(scriptfile_name, "w");
  if (outfile == NULL)
    fc_exitIfErrorPrintf(FC_FILE_IO_ERROR, "Could not open file '%s' to write "
                   "script", scriptfile_name);
 
  // Write header 
  fprintf(outfile, "# Script to make plots, generated by FCLib's analyzeSpotWelds\n");
  fprintf(outfile, "# %s", ctime(&start_time));
  fprintf(outfile, "#\n");
  fprintf(outfile, "# Depends on files '%s' and '%s'\n", outputfile_name, 
          gnuplotfile_name);
  fprintf(outfile, "# Uses programs 'gnuplot' and 'convert' (an ImageMagick program)\n");
  fprintf(outfile, "#\n");
  fflush(NULL);

  // write code to make graphs using gnuplot
  temp_string = malloc(sizeof(char)*(strlen(gnuplotfile_name)+12));
  fprintf(outfile, "gnuplot %s\n", gnuplotfile_name);
  fflush(NULL);
  free(temp_string);
  
  // write code to convert eps to jpg
  for (i = 0; i < numType + 1; i++) {
    temp_string = malloc(sizeof(char)*(2*strlen(plotfile_names[i]) + 40));
    fprintf(outfile, "convert -density 144 %s %s\n", plotfile_names[i],
            jpgfile_names[i]);
    fflush(NULL);
    free(temp_string);
  }

  // close file
  fclose(outfile);

  // make it executable
  temp_string = malloc(sizeof(char)*(strlen(scriptfile_name) + 20));
  sprintf(temp_string, "chmod +x %s\n", scriptfile_name);
  system(temp_string);
  free(temp_string);


  //*****************************************************************
  //*** Cleanup
  //*****************************************************************
  
  free(numBrokens);
  free(times);
  for (i = 0; i < numType; i++) {
    free(types[i].spotwelds);
    free(types[i].spotweld_names);
    free(types[i].nodeset_names);
    free(types[i].sideset_names);
    free(types[i].nodesets);
    free(types[i].sidesets);
    for (j = 0; j < types[i].numSpotWeld; j++)
      free(types[i].node_displs[j]);
    free(types[i].node_displs);
    for (j = 0; j < types[i].numSpotWeld; j++)
      free(types[i].side_displs[j]);
    free(types[i].side_displs);
    for (j = 0; j < types[i].numSpotWeld; j++)
      free(types[i].fails[j]);
    free(types[i].fails);
    free(types[i].group_fails);
    free(types[i].fail_stepIds);
    free(types[i].numBrokens);
    free(types[i].fail_times);
    //for (j = 0; j < types[i].numSpotWeld; j++) {
    //free(types[i].spotweld_names[j]);
    //free(types[i].nodeset_names[j]);
    //free(types[i].sideset_names[j]);
    //}
   }
  free(types);
  for (i = 0; i < numType + 1; i++) {
    free(plotfile_names[i]);
    free(jpgfile_names[i]);
  }
  free(plotfile_names);
  free(jpgfile_names);
  free(outputfile_name);
  free(gnuplotfile_name);
  free(scriptfile_name);

  // final the library
  fc_freeSierraInfo(sierra_info);
  // FIX? move this sooner?
  fc_deleteDataset(dataset);
  fc_finalLibrary();

  if (beVerbose) {
    printf("Done.\n");
    fflush(NULL);
  }

  exit(0); 
}

// convert void* array to doubles
double* get_doubles(int num, FC_DataType dataType, void* data) {
  int i;
  double* array;

  // make room
  array = malloc(sizeof(double)*num);
  
  switch(dataType) {
  case FC_DT_DOUBLE: {
    memcpy(array, data, num*sizeof(double));
  } break;
  case FC_DT_INT: {
    int* cast_data = (int*)data;
    for (i = 0; i < num; i++)
      array[i] = cast_data[i];
  } break;
  case FC_DT_FLOAT: {
    float* cast_data = (float*)data;
    for (i = 0; i < num; i++)
      array[i] = cast_data[i];
  } break;
  default:
    fprintf(stderr, "Warning:In helper function 'get_doubles', cannot convert "
            "data type '%s' to double\n", fc_getDataTypeText(dataType));
    fflush(NULL);
    free(array);
    array = NULL;
  }

  return array;
}

/** This function determines the shortest distinguishable name for each
 *  spot weld type. The name is a string listing the properties 
 *  of the spot weld types that differ from the other spotwelds in the
 *  array.
 *
 *  Procedure:
 *  -# Determine a base name for each type from the function names.
 *     Strips 'normal' and 'tangential' from the names of the functions to 
 *     see if the root names are the same. If they are, then the root name 
 *     is used. Other wise 'normal_name/tangent_name' is used.
 *  -# Check to see if names all differ. If they don't, switch
 *     all names back to the 'normal_name/tangent_name' form of the names. 
 *  -# Check to see what other parameters differ, and append 
 *     their values to all type's names.
 *
 *  \todo Might want the program to be even smarter and only append
 *        the different properties when they are necessary.
 *        For example, right now, if any two types have the same
 *        name and they have different scale factors, all types
 *        will get scale factors appended. It would be more concise
 *        to have the scale factors only added to the names of
 *        the two types that were equal.
 */
void nameSpotWeldTypes(int numType, SpotWeldType* types) {
  int i, j;
  int nameSame = 0;
  int decayDiff = 0;
  int normScaleDiff = 0;
  int tangScaleDiff = 0;
  int expDiff = 0;
  SpotWeldType *type0 = &types[0];
  
  // make base name
  for (i = 0; i < numType; i++) {
    SpotWeldType *type = &types[i];
    char temp_normal[1024] = "normal";
    char temp_tangent[1024] = "tangential";

    // First, check to see if they are equal
    if (!strcasecmp(type->normal_name, type->tangent_name)) {
      strcpy(temp_normal, type->normal_name);
      strcpy(temp_tangent, type->tangent_name);
    }
    else {
      int start;          // index to where strings first differ
      int norm_len = strlen(type->normal_name);
      int tang_len = strlen(type->tangent_name);
      int min_len;
      int found_norm;  // holds number of characters of norm string to skip
      int found_tang;  // holds number of characters of tangent str to skip

      if (norm_len < tang_len) 
        min_len = norm_len;
      else 
        min_len = tang_len;

      // Walk from start, until they differ
      for (start = 0; start < min_len; start++) {
        if (strncasecmp(&(type->normal_name[start]), 
                        &(type->tangent_name[start]), 1))
          break;
      }

      // check for variations of normal and tangential
      if (!strncasecmp(&(type->normal_name[start]), "normal", 6))
        found_norm = 6;
      else if (!strncasecmp(&(type->normal_name[start]), "norm", 4))
        found_norm = 4;
      if (!strncasecmp(&(type->tangent_name[start]), "tangential", 10))
        found_tang = 10;
      else if (!strncasecmp(&(type->tangent_name[start]), "tangent", 7))
        found_tang = 7;
      else if (!strncasecmp(&(type->tangent_name[start]), "tang", 4))
        found_tang = 4;

      if (found_norm > 0 && found_tang > 0) {
        // check for delimiters before and after (only need to take one)
        if (start > 0 && 
            (type->normal_name[start-1] == '_' || 
             type->normal_name[start-1] == '-')   ) {
          start--;
          found_norm++;
          found_tang++;
        }
        else {
          if (start + found_norm < norm_len  &&
              (type->normal_name[start + found_norm] == '_' ||
               type->normal_name[start + found_norm] == '-')   ) 
          found_norm++;
          if (start + found_tang < tang_len &&
              (type->tangent_name[start + found_tang] == '_' ||
               type->tangent_name[start + found_tang] == '-')  )
            found_tang++;
        }
        
        // copy all but in between bit
        if (start > 0) {
          strncpy(temp_normal, type->normal_name, start);
          strncpy(temp_tangent, type->tangent_name, start);
        }
        temp_normal[start] = '\0';
        temp_tangent[start] = '\0';
        if (start + found_norm < norm_len)
          strcat(temp_normal, &(type->normal_name[start + found_norm]));
        if (start + found_tang < tang_len)
          strcat(temp_tangent, &(type->tangent_name[start + found_tang]));

        //printf("norm_len= %d\n", norm_len);
        //printf("tang_len= %d\n", tang_len);
        //printf("start = %d\n", start);
        //printf("found_norm = %d\n", found_norm);
        //printf("found_tang = %d\n", found_tang);
        //printf("normal_name = %s\n", type->normal_name);
        //printf("temp_norm = %s\n", temp_normal);
        //printf("tangent_name = %s\n", type->tangent_name);
        //printf("temp_tang = %s\n", temp_tangent);
      }
    }

    // write base
    if (!strcmp(temp_normal, temp_tangent))
      sprintf(type->name, "%s", temp_normal);
    else
      sprintf(type->name, "%s/%s", type->normal_name, type->tangent_name);
  }
    
  // Check to see if any type names are the same
  for (i = 0; i < numType; i++) {
    for (j = i+1; j < numType; j++) {
      if (!strcmp(types[i].name, types[j].name)) {
        nameSame = 1;
        break;
      }
    }
    if (nameSame)
      break;
  }

  // if any are the same, switch back to long names
  if (nameSame) {
    for (i = 0; i < numType; i++) {
      SpotWeldType *type = &types[i];
      sprintf(type->name, "%s/%s", type->normal_name, type->tangent_name);
    }
    return;
  }
  
  // figure out what differs
  for (i = 1; i < numType; i++) {
    SpotWeldType *type = &types[i];
    if (decayDiff == 0 && type->decay_cycles != type0->decay_cycles)
      decayDiff = 1;
    if (normScaleDiff == 0 && type->normal_scale != type0->normal_scale)
      normScaleDiff = 1;
    if (tangScaleDiff == 0 && type->tangent_scale != type0->tangent_scale)
      tangScaleDiff = 1;
    if (expDiff == 0 && type->exponent != type0->exponent)
      expDiff = 1;
  }

  // append names of things that differ
  if (decayDiff || normScaleDiff || tangScaleDiff || expDiff) {
    for (i = 0; i < numType; i++) {
      SpotWeldType *type = &types[i];
      if (normScaleDiff) 
        sprintf(&(type->name[strlen(type->name)]), "/normal_scale_%g", 
                type->normal_scale);
      if (tangScaleDiff) 
        sprintf(&(type->name[strlen(type->name)]), "/tangent_scale_%g", 
                type->tangent_scale);
      if (decayDiff) 
        sprintf(&(type->name[strlen(type->name)]), "/decay_cycles_%d", 
                type->decay_cycles);
      if (expDiff) 
        sprintf(&(type->name[strlen(type->name)]), "/exponent_%g", 
                type->exponent);
    }
  }
}

