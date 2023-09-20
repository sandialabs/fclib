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
 * \file regionreassembler.c
 * \brief if have run subsetter with promote, works on out exodus
 *        file to copy var back to the main mesh using the map
 *
 *  reads in geom file and subsequent file where we have promoted
 *  meshes with well known naming convention (has timestep and segment id)
 *  and mappings for the orig mesh. reconsititutes all the vars from these
 *  promoted meshed onto the main mesh and fills in all the missing data points
 *  with value 0;
 *
 *  note: we cannot distinguish what were vars vs seqvars
 *  on the original mesh, so everything is treated as a step in a seqvar
 *
 *  writes the reassembled output to the named output file or to
 *  reassmebled.ex2
 *
 * \modifications
 *   - 04/26/08 ACG 
 *   - 06/20/08 ACG revised becuase of changes to subsetter
 *   - 07/13/08 ACG revised becuase of changes to subsetter
 *   - 08/11/08 ACG hardening
 */

#include <string.h>
#include "exodusII.h"
#include "fc.h"

static FC_Dataset ds;
static FC_Dataset rds;

static FC_VerbosityLevel verbose_level = FC_QUIET;
static char* geom_file_name = NULL;
static char* reduced_file_name = NULL;
static char* outfile = NULL;
// this is the fill val. can change this to be asked for and/or per variable 
static int* fillval = 0; 

static void getArgs(int, char**);
static int getMatchingMesh(char* origname,FC_Mesh **mainmeshes, int numMainMesh, 
			   int *imainmesh, int *timestep);

int main(int argc, char** argv) {
  FC_ReturnCode rc;
  int i,j;
  int numMesh, numMainMesh;
  FC_Mesh *meshes, *mainmeshes;
  FC_Sequence *seqs;
  int numseq, numstep;

  getArgs(argc, argv);

  // init library and load dataset
  rc = fc_setLibraryVerbosity(verbose_level);
  fc_exitIfError(rc);
  rc = fc_initLibrary();
  fc_exitIfError(rc);
  rc = fc_loadDataset(geom_file_name, &ds);
  fc_exitIfErrorPrintf(rc,"Can't load dataset\n");

  rc = fc_loadDataset(reduced_file_name, &rds);
  fc_exitIfErrorPrintf(rc,"Can't load dataset\n");

  // get meshes
  rc = fc_getMeshes(ds, &numMainMesh, &mainmeshes);
  if (numMainMesh < 1){
    printf("No mainmeshes were found\n");
    fflush(NULL);
    exit(0);
  } else {
    printf("%d main meshes\n", numMainMesh);
    fflush(NULL);
  }

  rc = fc_getMeshes(rds, &numMesh, &meshes);
  if (numMesh < 1){
    printf("No meshes were found\n");
    fflush(NULL);
    exit(0);
  } else {
    printf("%d reduced meshes\n", numMesh);
    fflush(NULL);
  }

  //there must be a unique seq on the main mesh dataset for reconstituting
  rc = fc_getSequences(ds,&numseq,&seqs);
  fc_exitIfErrorPrintf(rc, "Can't get sequence upon which to reconstitute");
  if (numseq != 1){
    fc_exitIfErrorPrintf(FC_ERROR, "Can't get unique sequence upon which to reconstitute");
  }
  rc = fc_getSequenceNumStep(seqs[0],&numstep);
  fc_exitIfErrorPrintf(rc, "Can't get sequence numstep");


  // loop over meshes
  for (i = 0; i < numMesh; i++) {
    FC_Variable *reducedvars;
    int numreducedvars;
    int imainmesh;
    int timestep;
    char* origname;
    char* mainname;

    //    printf("mesh %d\n",i);
    
    rc = fc_getMeshName(meshes[i], &origname);
    if (rc != FC_SUCCESS){
      fc_printfErrorMessage("cant get mesh by name");
      continue;
    }

    getMatchingMesh(origname,&mainmeshes,numMainMesh, 
		    &imainmesh, &timestep);
    if (imainmesh == -1 || timestep == -1){
      printf("mesh <%s> does not match a main mesh\n",origname);
      free(origname);
      continue;
    } 

    rc = fc_getMeshName(mainmeshes[imainmesh], &mainname);
    if (rc != FC_SUCCESS){
      fc_printfErrorMessage("cant get mesh by name");
      continue;
    }
    
    printf("mesh <%s> matches mainmesh <%s> at timestep <%d>\n", origname, mainname, timestep);
    free(mainname);
    free(origname);

    //if we dont have a seq var with this name make one.
    //then put the var data in the existing var at the correct step.
    
    //get only basic vars
    rc = fc_getVariables(meshes[i],&numreducedvars,&reducedvars);
    if (rc != FC_SUCCESS){
      fc_printfErrorMessage("cant get vars from mesh");
      continue;
    }

    for (j = 0; j < numreducedvars; j++){
      char* currvarname;
      FC_Variable *dest_seqvar;
      FC_AssociationType assoc;
      FC_Variable *mapvars;
      int nummapvar;
      
      rc = fc_getVariableName(reducedvars[j],&currvarname);
      if (rc != FC_SUCCESS){
	fc_printfErrorMessage("cant get var name");
	continue;
      }
      //      printf("Considering var <%s>\n",currvarname);
      
      //skip for exodus if var is ID, or if it is a mapping var)
      if (strcmp(currvarname, "ID") == 0 ||
	  strcmp(currvarname,"FC_ELEMENTMAP") == 0 ||
	  strcmp(currvarname,"FC_VERTEXMAP") == 0){
	fc_printfWarningMessage("not copying over var '%s'",currvarname);
	free(currvarname);
	continue;
      }

      rc = fc_getVariableAssociationType(reducedvars[j], &assoc);
      if (rc != FC_SUCCESS){
	fc_exitIfErrorPrintf(rc,"can't get var assoc type");
      }
      switch(assoc){
      case FC_AT_ELEMENT:
	rc = fc_getVariableByName(meshes[i],"FC_ELEMENTMAP",&nummapvar,&mapvars);
	fc_exitIfErrorPrintf(rc,"Can't get elem map var");
	break;
      case FC_AT_VERTEX:
	rc = fc_getVariableByName(meshes[i],"FC_VERTEXMAP",&nummapvar,&mapvars);
	fc_exitIfErrorPrintf(rc,"Can't get vertex map var");
	break;
      default:
	fc_printfWarningMessage("Wont copy basic vars of type %d", assoc);
	free(currvarname);
	continue;
      }
      if (nummapvar != 1){
	fc_printfErrorMessage("Can't get map var");
	exit(0);
      }

      printf("Reassembling var <%s> at timestep %d\n",currvarname,timestep);
      rc = fc_copySeqVariableStepFromRegionMesh(reducedvars[j],
						mainmeshes[imainmesh],
						seqs[0],
						timestep,
						mapvars[0],
						fillval,
						currvarname,
						&dest_seqvar);
      free(dest_seqvar);
      free(currvarname);
      free(mapvars);
    } //reduced vars
    free(reducedvars);
  } //nummesh
  
  free(meshes);
  free(mainmeshes);
  free(seqs);


  if (outfile != NULL){
    rc = fc_rewriteDataset(ds, outfile, FC_FT_EXODUS);
  } else {
    //    printf("trying to write to a different file\n");
    rc = fc_rewriteDataset(ds, "reassembled.ex2", FC_FT_EXODUS);
  }
  fc_exitIfError(rc);

  rc = fc_deleteDataset(ds);
  rc = fc_deleteDataset(rds);
  fc_exitIfError(rc);

  fc_finalLibrary();
  exit(0);
}


int getMatchingMesh(char* origname,FC_Mesh **mainmeshes, int numMainMesh, 
		    int *imainmesh, int *timestep){

  FC_ReturnCode rc;
  int i;

  //now looking for names like mainmeshname_T#_Seg#
  //where T is the timestep and Seg# is segment at that timestep.    
  //if it has the right naming convention then this is a mesh we want
  
  //defaults
  *imainmesh = -1;
  *timestep = -1;

  for (i = 0; i < numMainMesh; i++){
    char* mainname;
    rc = fc_getMeshName((*mainmeshes)[i],&mainname);
    if (rc != FC_SUCCESS){
      fc_printfErrorMessage("cant get mesh by name");
      continue;
    }

    //    printf("comparing <%s>(%d) with <%s>(%d)\n",origname,strlen(origname),
    //	   mainname, strlen(mainname));
    //origname must be at least mainname_T#_Seg#
    if (strncmp(origname,mainname,strlen(mainname)) == 0 &&
	strlen(origname) >= strlen(mainname)+7 &&
	origname[strlen(mainname)] == '_' &&
	origname[strlen(mainname)+1] == 'T'){

      char holdnum[10]; //FIXME: this cant be fixed
      int lindex = strlen(mainname)+2;
      int currdigit = 0;

      //      printf("matched base now looking for number\n");
      while(strlen(origname) >= (uint)(lindex+5)){
	if (origname[lindex] != '_'){
	  //FIXME -- is there a regex thing I can use?
	  //	  printf("eval <%c>\n",origname[lindex]);
	  holdnum[currdigit++] = origname[lindex++];
	  if (currdigit > 8){ //leave room for the \0
	    fc_exitIfErrorPrintf(FC_ERROR,"Programmer error - cannot handle timestep size");
	  }
	} else {
	  //FIXME - should check rest of name
	  //now we have a match
	  char* endptr;
	  long lcheck;
	  holdnum[currdigit] = '\0';
	  
	  lcheck = strtol(holdnum, &endptr, 10);
	  if (*endptr == '\0'){ // its a number
	    *timestep = (int)lcheck;
	    //	    printf("timestep = <%d>\n",*timestep);
	    *imainmesh = i;
	    free(mainname);
	    return 1;
	  } //else its not a number
	  break;
	}
      }
    } //checked against this meshname
    
    free(mainname);
    
  } //loop over main meshes

  return 0;
}


void getArgs(int argc, char** argv){
  // handle arguments
  int i;

  if (argc < 2) {
  usage:
    printf("usage: %s [options] geomdataset reducedmeshesdataset \n", argv[0]);
    printf("options: \n");
    printf("   -h           : print this help message\n");
    printf("   -v           : verbose: print warning and error messages\n");
    printf("   -V           : very verbose: print log and error messages\n");
    printf("   -o outfile   : outfile (exodus) \n");
    printf("\n");
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
    else if (!strcmp(argv[i], "-o")){
      if (i+1 > argc)
	goto usage;
      outfile = argv[++i];
    }
    else if (!strncmp(argv[i], "-h", 2) || !strncmp(argv[i], "-", 1))
      goto usage;
    else {
      if (i+2 > argc)
        goto usage;
      geom_file_name = argv[i];
      reduced_file_name = argv[i+1];
      i+=2;
    }
  }
  if (!geom_file_name || !reduced_file_name )
    goto usage;
}
