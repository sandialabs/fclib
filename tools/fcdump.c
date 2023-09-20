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
 * \file fcdump.c
 * \brief Read a dataset and spit out it's information to standard output
 *
 * $Source: /home/Repositories/fcdmf/fclib/tools/fcdump.c,v $
 * $Revision: 1.31 $ 
 * $Date: 2007/04/17 20:57:26 $
 *
 * \modifications
 *   - nsc 10dec01
 *   - updated for new interface, 18apr02
 *   - wsk Jul 2, 02  Moved coordinates & connectivities out of variable loop
 *                Also some general fixing up
 *   - 2003-NOV-14  WSK New & improved, made look more like scovdump.
 *                  Flag lets you see or hide big data. Better support 
 *                  of sequence variables.
 *   - 05/25/04 WSK At some point in the past, add -V flag for two
 *       levels of verbosity. Now removing most error messages local to
 *       this file because they will show up with -v.
 *   - 05/27/04 WSK Added subsets.
 *   - 01/13/05 WSK added ability to read ls-dyna files. Right now the
 *     model is to iterate over readers until one works (or all don't work).
 *   - 04/17/07 CDU added filters
 */

#include <string.h>
#include <fc.h>




//Structure for letting the user specify what they're looking for
typedef struct {
  char *mesh_name;  //"" means any mesh
  char *var_name;   //"" means any variable
  int   step_start; //-1 means all steps
  int   step_stop;  //-1 means all steps
} filter_t;

typedef enum itemType_e {
  NONE=0x00, ALL=0xFF,
  DATASET=0x01, SEQUENCE=0x02, GLOBAL_VAR=0x04, GLOBAL_SEQVAR=0x08,
  MESH=0x10,      SUBSET=0x20,        VAR=0x40,        SEQVAR=0x80,
} item_t;

typedef struct {
  char       *name;
  item_t      id;
} item_obj_t;

item_obj_t search_items[] = {
  { "ds",      DATASET },        { "dataset",                DATASET },
  { "seq",     SEQUENCE},        { "sequence",               SEQUENCE},        
  { "gvar",    GLOBAL_VAR },     { "globalvariable",         GLOBAL_VAR },
  { "gseqvar", GLOBAL_SEQVAR },  { "globalsequencevariable", GLOBAL_SEQVAR },
  { "m",       MESH },           { "mesh",                   MESH},
  { "sub",     SUBSET },         { "subset",                 SUBSET },
  { "var",     VAR},             { "variable",               VAR},
  { "seqvar",  SEQVAR },         { "sequencevariable",       SEQVAR },
  { NULL,      NONE } }; //Null terminates list

typedef struct {
  char *name;
  int   spot;
} index_names_t;

index_names_t index_names[] = {
  { "first",     -1},
  { "start",     -1},
  { "begin",     -1},
  { "last",      -2},
  { "stop",      -2},
  { "end",       -2},
  { NULL,        -1} }; //Null terminates list



static int isVariableInFilterList(char *meshName, char *varName, 
				  filter_t *f, int f_num,
				  int *start, int *stop){
  int i;
  if(f_num==0) { //no filters=all selected
    if(start) *start = -1;
    if(stop)  *stop  = -2;
    return 1; 
  }
  for(i=0;i<f_num;i++){
    if((!strcmp(f[i].mesh_name, meshName) && !strlen(f[i].var_name))          ||   //all of mesh
       (!strlen(f[i].mesh_name)           && !strcmp(f[i].var_name, varName)) ||   //all of var 
       (!strcmp(f[i].mesh_name, meshName) && !strcmp(f[i].var_name, varName))    ){//both match
      if(start) *start = f[i].step_start;
      if(stop)  *stop  = f[i].step_stop;
      return 1; 
    }
  }
  return 0;
}

static int getIndexFromString(char *s){
  index_names_t *idx;

  if(!s) return -2; //Null string, return end

  for(idx=index_names; idx->name; idx++)
    if(!strcasecmp(idx->name, s))
      return idx->spot;

  //Not found, must be a number
  return atoi(s);  
}

static void dumpHelp(void){

  printf("\n");
  printf("fcdump is a tool in the Feature Characterization Library (FCLib) that prints out\n");
  printf("information about a dataset. Users can specify filters on the command line to\n");
  printf("locate specific items of interest in the dataset.\n");
  printf("\n");        
  printf("Usage:\n");
  printf("  fcdump [-options]  <-s item1> <-s item2>  <-f filter1> <-f filter2>  dataset\n");
  printf("\nOptions:\n");
  printf("   -a           : print all data values\n");
  printf("   -A           : print all data values, using exponent notation\n");
  printf("   -h, -help    : print this help message\n");
  printf("   -v,          : verbose: print warning and error messages\n");
  printf("   -V,          : very verbose: print log and error messages\n");
  printf("\nCustomizations:\n");
  printf("   <-s item>    : show only specific items, where item is either\n");
  printf("       dataset                         - dataset info\n");
  printf("       seq                             - sequence info\n");
  printf("       gvar                            - global variable info\n");
  printf("       gseqvar                         - global sequence variable info\n");
  printf("       mesh                            - mesh info\n");
  printf("       subset                          - mesh subset info\n");
  printf("       var                             - mesh variable info\n");
  printf("       seqvar                          - mesh sequence variable info\n");
  printf("\n");
  printf("   <-f filter>  : only show filtered values, where a filter is\n");
  printf("       variable_name                   - a variable in all meshes\n");
  printf("       variable_name[s]                - step s of a variable\n");
  printf("       variable_name[s1:s2]            - steps s1-s2 of a variable\n");
  printf("       /mesh_name                      - all info about a mesh\n");
  printf("       /mesh_name/variable_name        - a specific variable in a mesh\n");
  printf("       /mesh_name/variable_name[s]     - step s of a variable in a mesh\n");
  printf("       /mesh_name/variable_name[s1:s2] - steps s1-s2 of a variable in a mesh\n");
  printf("\n");
  printf("      Note: sequence variable ranges can be singular, forwards, or backwards:\n");
  printf("        [3]                            - step  3\n");
  printf("        [3:5]                          - steps 3,4,5\n");
  printf("        [5:3]                          - steps 5,4,3\n");
  printf("        [last]                         - last sequence\n");
  printf("        [last:first]                   - full sequence in reverse order\n"); 
  printf("\n");
  fflush(NULL);
  exit(0);
}

int main(int argc, char **argv)
{
  FC_ReturnCode rc;
  int i, j, k;
  int printBigData = 0;
  FC_VerbosityLevel verbosity = FC_QUIET;
  int needed;

  char* filename = NULL;
  int numSeq, numGlbVar, numGlbSeqVar, numMesh;
  int numSubset, numVar, numSeqVar, *numStepPerVar;
  char temp_string[1024];
  FC_Dataset dataset;
  FC_Sequence* sequences = NULL, temp_seq;
  FC_Mesh *meshes;
  FC_Subset *subsets;
  FC_Variable *glbVars, **glbSeqVars, *variables, **seqVars;
  int seqID;
  char *varName;
  unsigned int specific_items=0;
  filter_t *filters=NULL;
  int       filters_num=0;
  char *meshName;
  

  //---------------------------------------------------------------------
  //  process command line args
  //---------------------------------------------------------------------
  if (argc < 2) {
    dumpHelp();
  }
  
  for (i = 1; i < argc; i++) {
    if ((!strcmp(argv[i], "help")) ||
        (!strcmp(argv[i], "-help")) ||
        (!strcmp(argv[i], "-h")) ||
        (!strcmp(argv[i], "-")) ||
        (!strcmp(argv[i], "--"))) 
      dumpHelp();
    else if (!strcmp(argv[i], "-a")) {
      printBigData = 1;
    }
    else if (!strcmp(argv[i], "-A")) {
      printBigData = 2;
    }
    else if (!strcmp(argv[i], "-f")) {
      char *t1,*t2, *t3, *t4, *t5=NULL;
      int  step_start=-1, step_stop = -2;
      filter_t *ff;
      i++;
      if(i==argc) dumpHelp(); 

      t1 = t2 = t3 = argv[i];
      filters = (filter_t *) realloc(filters, (filters_num+1)*sizeof(filter_t));
      ff = &filters[filters_num];

      //--Parse variable name to see if there is a range--
      // Supports: [10]   10  
      //           [3:5]  3,4,5
      //           [:5]   0,1,2,3,4,5
      //           [5:]   5,6,7... to end
      //Scan for [
      while((*t3!='\0')){
	if(*t3=='['){
	  *t3='\0';
	  t4 = t3+1;
	  //Scan for : or ]
	  while(*t4!='\0'){
	    if(*t4==']') 
	      *t4='\0';
	    else {
	      if(*t4==':'){
		*t4='\0';
		t5=t4+1;
	      }
	      t4++;
	    }
	  }

	  //Assign start and stop values
	  step_start = getIndexFromString(t3+1);      
	  step_stop  = (t5==NULL)  ? step_start :   //No second value
                              	     (*t5=='\0') ? -2 :      //Used a [:],  no value given
                                                   getIndexFromString(t5); //Used a [:x]  given a value
	 
	} else
	  t3++;
      }

      //Check for /mesh and /mesh/variable
      if(t1[0]=='/'){
	//Given at least a mesh
	t2 = ++t1;
	while((*t2!='\0') &&(*t2!='/'))
	  t2++;
	if(*t2=='/') { //is a /mesh/variable name
	  *t2='\0';
	  t2++;
	}
	ff->mesh_name = strdup(t1);
	ff->var_name  = strdup(t2);    
      } else { //Just a variable name
	ff->mesh_name = "";
	ff->var_name  = strdup(t1);	
      }
      //Everyone gets a start/stop
      ff->step_start = step_start;
      ff->step_stop  = step_stop;

      filters_num++;

    } else if (!strcmp(argv[i], "-s")) {
      i++;
      if(i==argc) dumpHelp(); 
      j=0;
      while((search_items[j].name) && 
	    (strcasecmp(search_items[j].name, argv[i])))
	j++;
      if(!search_items[j].name) dumpHelp();
      specific_items |= search_items[j].id;

    } else if (!strcmp(argv[i], "-v")) {
      verbosity = FC_WARNING_MESSAGES;
    }
    else if (!strcmp(argv[i], "-V")) {
      verbosity = FC_LOG_MESSAGES;
    }
    else 
      filename = argv[i];
  }
  if (!filename)
    dumpHelp();


  //Default is to display everything
  if(!specific_items) specific_items=ALL;


  //---------------------------------------------------------------------
  //  Initialize the library
  //---------------------------------------------------------------------
  rc = fc_setLibraryVerbosity(verbosity);
  fc_exitIfError(rc);
  rc = fc_initLibrary();
  fc_exitIfError(rc);
  
  //---------------------------------------------------------------------
  //  Load the file
  //---------------------------------------------------------------------
  rc = fc_loadDataset(filename, &dataset);
  fc_exitIfErrorPrintf(rc, "Failed to load dataset '%s'", filename);
  if(specific_items & DATASET){
    fc_printDataset(dataset, NULL);
  }
  //---------------------------------------------------------------------
  //  Query sequences
  //---------------------------------------------------------------------
  rc = fc_getSequences(dataset, &numSeq, &sequences);
  fc_exitIfError(rc);

  if(specific_items & SEQUENCE){
  
    // loop over sequences
    for (i = 0; i < numSeq; i++) {

      //See if we're listing this sequence
      rc = fc_getSequenceName(sequences[i], &varName);
      fc_exitIfError(rc);
      if(!isVariableInFilterList("", varName, filters, filters_num, NULL, NULL)){
	free(varName);
	continue;
      }
      free(varName);
      
      // print the sequence
      sprintf(temp_string, "\nSequence %d", i);
      rc = fc_printSequence(sequences[i], temp_string, printBigData);
      fc_exitIfError(rc);

    }
  }

  //---------------------------------------------------------------------
  //  Query Global variables
  //---------------------------------------------------------------------
  if(specific_items & GLOBAL_VAR){

    rc = fc_getGlobalVariables(dataset, &numGlbVar, &glbVars);
    fc_exitIfError(rc);
  
    // loop over global variables
    for (i = 0; i < numGlbVar; i++) {

      //See if this global variable is being listed
      rc = fc_getVariableName(glbVars[i], &varName);
      fc_exitIfError(rc);
      if(!isVariableInFilterList("", varName, filters, filters_num, NULL, NULL)){
	fc_deleteVariable(glbVars[i]);
	free(varName);
	continue;
      }
      free(varName);

      // print the global variable
      sprintf(temp_string, "\nGlobal Variable %d", i);
      rc = fc_printVariable(glbVars[i], temp_string, printBigData);
      fc_exitIfError(rc);
      
      // delete variable
      fc_deleteVariable(glbVars[i]);
    } // end of loop over global variables
    free(glbVars);
  }

  //---------------------------------------------------------------------
  //  Query Global Sequence Variables
  //---------------------------------------------------------------------
  // Note: GSV's have a sequence associated with them. Since we're
  //       deleting variables as we go, we have to access the GSV and
  //       delete them here, even if we don't use them. Otherwise,
  //       when we delete the sequence at the end and do a 
  //       deleteDataset(), we'll vars without sequences to be deleted.
  rc = fc_getGlobalSeqVariables(dataset, &numGlbSeqVar,
				&numStepPerVar, &glbSeqVars);
  fc_exitIfError(rc);

     
  // loop over global seq vars
  for (i = 0; i < numGlbSeqVar; i++) {
    int start,stop;

    if(specific_items & GLOBAL_SEQVAR){
      //See if this global seq var is in the filter list
      rc = fc_getVariableName(glbSeqVars[i][0], &varName);
      fc_exitIfError(rc);

      if(isVariableInFilterList("", varName, filters, filters_num, &start, &stop)){
	// print the seq variable
	rc = fc_getSequenceFromSeqVariable(numStepPerVar[i], glbSeqVars[i],
					   &temp_seq);
	fc_exitIfError(rc);
	seqID = -1;
	for (k = 0; k < numSeq; k++) {
	  if (FC_HANDLE_EQUIV(sequences[k], temp_seq))
	    seqID = k;
	}
	sprintf(temp_string, "\nGlobal SeqVar %d (Seq %d)", i, seqID);
	rc = fc_printSeqVariableRange(numStepPerVar[i], glbSeqVars[i], 
				      start, stop,
				      temp_string, printBigData);
	fc_exitIfError(rc);
      }
      free(varName);
    }

    // delete global seq vars
    fc_deleteSeqVariable(numStepPerVar[i], glbSeqVars[i]);
    free(glbSeqVars[i]);
  }
  free(numStepPerVar);
  free(glbSeqVars);    

    


  //---------------------------------------------------------------------
  //  Query meshes
  //---------------------------------------------------------------------
  rc = fc_getMeshes(dataset, &numMesh, &meshes);
  fc_exitIfError(rc);



  
  // loop over meshes
  for (i = 0; i < numMesh; i++) {

    //Get info up front so we can make filter decision before printing
    rc = fc_getMeshName(meshes[i], &meshName);
    fc_exitIfError(rc);


    //See if we're looking for any vars in this mesh
    needed=0;

    //1: Check Subsets
    if((!needed) && (specific_items & SUBSET)){
      rc = fc_getSubsets( meshes[i], &numSubset, &subsets);
      fc_exitIfError(rc);
      for(j=0; j<numSubset; j++){
	if(!needed){
	  rc = fc_getSubsetName(subsets[j], &varName);
	  fc_exitIfError(rc);
	  needed = isVariableInFilterList(meshName, varName, filters, filters_num, NULL, NULL);
	  free(varName);
	}
      }
      free(subsets);
    }

    //2: Check Variables
    if((!needed) && (specific_items & VAR)){
      rc = fc_getVariables(   meshes[i], &numVar, &variables);
      fc_exitIfError(rc);
      for(j=0; j<numVar; j++){
	if(!needed){
	  rc = fc_getVariableName(variables[j], &varName);
	  fc_exitIfError(rc);
	  needed = isVariableInFilterList(meshName, varName, filters, filters_num, NULL, NULL);
	  free(varName);
	}
      }
      free(variables);
    }

    //3: Check SeqVariables
    if((!needed) && (specific_items & SEQVAR)) {
      rc = fc_getSeqVariables(meshes[i], &numSeqVar, &numStepPerVar, &seqVars);
      fc_exitIfError(rc);

      for(j=0; j<numSeqVar; j++){
	if(!needed){
	  rc = fc_getVariableName(seqVars[j][0], &varName);
	  fc_exitIfError(rc);
	  needed = isVariableInFilterList(meshName, varName, filters, filters_num, NULL, NULL);
	  free(varName);
	}
	free(seqVars[j]);
      }
      free(numStepPerVar);
      free(seqVars);
    }

    //Skip this mesh if it doesn't have something the user wants
    if(!needed){
      free(meshName);
      fc_deleteMesh(meshes[i]);
      continue; //move on to next mesh
    }

     
    // print the mesh
    if(specific_items & MESH){
      sprintf(temp_string, "\nMesh %d", i);
      rc = fc_printMesh(meshes[i], temp_string, printBigData, printBigData, 
			printBigData, printBigData);
      fc_exitIfError(rc);
    }
    
    //---------------------------------------------------------------------
    //  Query subsets
    //---------------------------------------------------------------------
    if(specific_items & SUBSET){

      rc = fc_getSubsets(meshes[i], &numSubset, &subsets);
      fc_exitIfError(rc);
    
      // loop over subsets
      for (j = 0; j < numSubset; j++) {
	
	//See if this subset is in the filter list
	rc = fc_getSubsetName(subsets[j], &varName);
	fc_exitIfError(rc);
	if(!isVariableInFilterList(meshName, varName, filters, filters_num, NULL, NULL)){
	  fc_deleteSubset(subsets[j]);
	  free(varName);
	  continue;
	}
	free(varName);

	// print the subset
	sprintf(temp_string, "\nSubset %d (Mesh %d)", j, i);
	rc = fc_printSubset(subsets[j], temp_string, printBigData);
	fc_exitIfError(rc);
	
	// delete subset
	fc_deleteSubset(subsets[j]);
      } // end of loop over subsets
      free(subsets);
    }

    //---------------------------------------------------------------------
    //  Query variables
    //---------------------------------------------------------------------
    if(specific_items & VAR){
    
      rc = fc_getVariables(meshes[i], &numVar, &variables);
      fc_exitIfError(rc);

      // loop over variables
      for (j = 0; j < numVar; j++) {
	rc = fc_getVariableName(variables[j], &varName);
	fc_exitIfError(rc);
	if(!isVariableInFilterList(meshName, varName, filters, filters_num, NULL, NULL)){
	  fc_deleteVariable(variables[j]);
	  free(varName);
	  continue;
	}
	free(varName);

	// print the variable
	sprintf(temp_string, "\nVariable %d (Mesh %d)", j, i);
	rc = fc_printVariable(variables[j], temp_string, printBigData);
	fc_exitIfError(rc);
	
	// delete variable
	fc_deleteVariable(variables[j]);
      } // end of loop over variables
      free(variables);
    }


    //---------------------------------------------------------------------
    //  Query Sequence Variables
    //---------------------------------------------------------------------
    if(specific_items & SEQVAR){

      rc = fc_getSeqVariables(meshes[i], &numSeqVar, &numStepPerVar, &seqVars);
      fc_exitIfError(rc);
 
      // loop over seq vars
      for (j = 0; j < numSeqVar; j++) {
	int start, stop;

	//See if we need to do this one
	rc = fc_getVariableName(seqVars[j][0], &varName);
	fc_exitIfError(rc);
	if(!isVariableInFilterList(meshName, varName, filters, filters_num, &start,&stop)){
	  fc_deleteSeqVariable(numStepPerVar[j], seqVars[j]);
	  free(seqVars[j]);	
	  free(varName); 
	  continue;
	}     
	free(varName);
	

	// print the seq variable
	rc = fc_getSequenceFromSeqVariable(numStepPerVar[j], seqVars[j],
					   &temp_seq);
	fc_exitIfError(rc);
	seqID = -1;
	for (k = 0; k < numSeq; k++) {
	  if (FC_HANDLE_EQUIV(sequences[k], temp_seq))
	    seqID = k;
	}
	sprintf(temp_string, "\nSeqVar %d (Mesh %d, Seq %d)", j, i, seqID);
	rc = fc_printSeqVariableRange(numStepPerVar[j], seqVars[j],
				      start, stop,
				      temp_string, printBigData);
	fc_exitIfError(rc);
      
	// delete seq vars
	fc_deleteSeqVariable(numStepPerVar[j], seqVars[j]);
	free(seqVars[j]);

      }
      free(numStepPerVar);
      free(seqVars);
    }

    // delete mesh
    fc_deleteMesh(meshes[i]);
    free(meshName);                  
  }

  free(meshes);

  printf("\n");
  fflush(NULL);

  // delete the sequences
  for (i = 0; i < numSeq; i++)
    fc_deleteSequence(sequences[i]);
  free(sequences);
  
  fc_deleteDataset(dataset);
  fc_finalLibrary();  

  exit(FC_SUCCESS); 
}
