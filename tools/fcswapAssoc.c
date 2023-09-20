/**
 * \file fcswapAssoc.c
 * \brief Given a dataset, convert element data to vertex data, or vice versa. 
 * 
 * $Source: /usr/local/Repositories/fcdmf/fclib/tools/fcswapAssoc.c,v $
 * $Revision: 1.1 $
 * $Date: 2007/03/07 21:31:39 $
 *
 * \description
 * Usage: fcswapAssociation [options] input_file
 *   options:
 *    -o output_file_name
 *    -d [e2v|v2e]             -- direction: element-to-vertex or vertex-to-element
 *    -m mesh_name             -- specific mesh to work on
 *    -n variable_name         -- specific variable to work on
 *    -r new_label_extension   -- string appended to new variables
 *    -R                       -- remove original variable from output
 *    -v                       -- print warning messages
 *    -V                       -- print log messages
 *    -h                       -- print this help message
 *
 *        Given a dataset, create a new dataset whose variables have
 *        all been converted from element to vertex associations or
 *        vice versa. Through command line switches, the user can
 *        specify a specific variable/mesh to work on, an extension to
 *        append on to the new variable's name, and whether the original
 *        variable should be deleted. New variables are always type
 *        FC_DT_DOUBLE.
 *
 *        Note: Issues with the ExodusII file format limit just how
 *        selective the user can be with the -m -n options. If a variable
 *        is converted in one mesh, it is converted in all meshes.
 *
 *
 * \modifications
 *   - 03/06/2007 CDU created
 */
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <fc.h>


//=============================================================================
static int checkVarAndCreateName(FC_Variable v_in, 
				 FC_AssociationType search_assocType,
				 char  *target_var_name,
				 int    replace_var,
				 char  *mod_extension,
				 char **ret_new_name){
  
  char *temp_name, *new_name;
  FC_AssociationType v_atype;
  FC_DataType        v_dtype;
  FC_ReturnCode rc;

  *ret_new_name = NULL;

  //Get the variable's name to see if it matches
  rc = fc_getVariableName( v_in, &temp_name);
  fc_exitIfError(rc);
  rc = fc_getVariableInfo( v_in, NULL,NULL, &v_atype, NULL, &v_dtype);
  fc_exitIfError(rc);

  if( (v_atype != search_assocType)           ||     //Only elements
      !((v_dtype == FC_DT_DOUBLE) ||              //Only doubles, floats, ints
	(v_dtype == FC_DT_FLOAT)  ||
	(v_dtype == FC_DT_INT)      )      ||
      ((strlen(target_var_name)!=0) &&       //user specified a name and
       (strcmp(target_var_name,temp_name)))){ // name different than this
    *ret_new_name = NULL;
    free(temp_name);
    return 0;
  }

  //If replacing, then use name we were given
  if(replace_var){
    *ret_new_name = temp_name;
    return 1;
  }

  //Let someone else rename if we weren't given a modifier
  if(!strlen(mod_extension)){
    *ret_new_name = NULL;
    free(temp_name);
    return 1;
  }
  //Create a name for the new guy
  new_name = (char *)malloc(strlen(temp_name)+strlen(mod_extension)+1);
  if(!new_name)
    fc_exitIfErrorPrintf(FC_MEMORY_ERROR,"Memory Error");
  sprintf(new_name, "%s%s",temp_name, mod_extension);
  //fc_printfLogMessage("New name will be '%s'\n",new_name);
  //printf("New name will be '%s'\n",new_name);
  free(temp_name);
  *ret_new_name = new_name;
  return 1;
}

//=============================================================================
static void dumpHelp(void){
  printf("Usage: fcswapAssociation [options] input_file\n");
  printf("  options:\n");
  printf("     -o output_file_name\n");
  printf("     -d [e2v|v2e]             -- direction: element-to-vertex or vertex-to-element\n");
  printf("     -m mesh_name             -- specific mesh to work on\n");
  printf("     -n variable_name         -- specific variable to work on\n");
  printf("     -r new_label_extension   -- string appended to new variables\n");
  printf("     -R                       -- remove original variable from output\n");
  printf("     -v                       -- print warning messages\n");
  printf("     -V                       -- print log messages\n");
  printf("     -h                       -- print this help message\n");

  exit(0);
}


//=============================================================================
static void checkArgs(int argc,            char **argv,
		      FC_VerbosityLevel   *verbose_level,
		      FC_FileIOType       *new_fileType,
		      char                **file_in_name,
		      char                **file_out_name,
		      char                **mesh_name,
		      char                **var_name,
		      FC_AssociationType   *search_assocType,
		      FC_AssociationType   *target_assocType,
		      int                  *replace_var,
		      char                **mod_name){


  FC_FileIOType     ft[] = { FC_FT_EXODUS, FC_FT_EXODUS };
  char         *ft_ext[] = { ".exo",       ".ex2" };
  int           ft_num   = 2;
  char *infile_base;
  char *infile_ext;
  char *outfile_ext;
  int   type_set=0;
  int i,ich;

  if(argc<2) dumpHelp();


  //Parse input parameters
  while( (ich = getopt(argc, argv, "hvVm:n:o:t:d:r:R")) != EOF){
    switch(ich){

    case 'r': *mod_name      = optarg;              break;
    case 'R': *replace_var   = 1;                   break;
    case 'm': *mesh_name     = optarg;              break;
    case 'n': *var_name      = optarg;              break;
    case 'v': *verbose_level = FC_WARNING_MESSAGES; break;
    case 'V': *verbose_level = FC_LOG_MESSAGES;     break;

    case 'o': //Output file name
      *file_out_name = optarg;
      break;

    case 'd': //Direction
      if(!strcasecmp(optarg,"e2v"))      *search_assocType = FC_AT_ELEMENT;
      else if(!strcasecmp(optarg,"v2e")) *search_assocType = FC_AT_VERTEX;
      else dumpHelp();
      break;
	
    case 'h': 
    default: 
      dumpHelp(); break;
    }
  }
  if(optind < argc){
    *file_in_name = argv[optind];
  }

  //We at least need an input file name
  if(!strlen(*file_in_name)) dumpHelp();


  

  //Check the output file name
  if(strlen(*file_out_name)){
    //User gave us something. Make sure it's right type
    outfile_ext = fc_getExtension(*file_out_name,1);
    i=0;
    while((i<ft_num) && (strcasecmp(outfile_ext, &(ft_ext[i][1]))))
      i++;
    if( (i==ft_num) ||                                   //Unknown extension? 
        ((type_set) && (*new_fileType != ft[i])) ){ //Different extension than filename
      //The file name doesn't reflect the data. We need to append with
      //the proper extension
      i=0;
      while((i<ft_num) && (ft[i]!= *new_fileType))
	i++;
     
      *file_out_name = (char *)realloc(*file_out_name, 
				       strlen(*file_out_name)+strlen(ft_ext[i])+1);
      strcat(*file_out_name, ft_ext[i]);

    } else if(!type_set){
      //The user is setting the output type via file name
      *new_fileType = ft[i];
    }
    free(outfile_ext);

  } else {
    //User didn't give us a name. Create a simple one in local directory
    infile_base = fc_getBasenameWOExtension(*file_in_name, 1);

    //Figure out extension.. is there a call for this?
    i=0;
    while((i<ft_num) && (ft[i]!=*new_fileType))
      i++;
    //Note: Error check is handled when in arg case above

    infile_ext  = ft_ext[i];
    *file_out_name = (char *)malloc(strlen(infile_base) +
				    strlen("_mod")      +
				    strlen(infile_ext)  + 1);
    strcpy(*file_out_name, infile_base);
    strcat(*file_out_name, "_mod");
    strcat(*file_out_name, infile_ext);

    free(infile_base);
  }

  //Come up with the direction
  if      (*search_assocType == FC_AT_ELEMENT) *target_assocType = FC_AT_VERTEX;
  else if (*search_assocType == FC_AT_VERTEX)  *target_assocType = FC_AT_ELEMENT;
  else {
    printf("Unknown search association type?\n");
    dumpHelp();
  }
     
}


//=============================================================================
int main(int argc, char **argv){

  char *file_in_name = "";///datasets/can_crush.ex2";
  char *file_out_name = "";
  char *search_mesh_name     = "";
  char *search_variable_name = "";
  char *mod_extension = "";     //What we append to vars we find
  int   replace_var = 0;        //Replace the variable instead of duplicating it
  FC_ReturnCode rc;
  FC_VerbosityLevel verbose_level = FC_QUIET;
  FC_FileIOType new_fileType = FC_FT_EXODUS;
  FC_Dataset dataset;
  FC_Mesh *meshes;
  FC_Variable *vars, **seqVars;
  FC_AssociationType search_assocType = FC_AT_ELEMENT;
  FC_AssociationType target_assocType = FC_AT_VERTEX;
  int numMesh,numSeqVar, *numStepPerSeqVar;
  int numVar;
  char *temp_name, *new_name;
  int i,j;


  //Parse the args to figure out what we're supposed to do
  checkArgs(argc, argv, &verbose_level, &new_fileType, 
	    &file_in_name, &file_out_name,
	    &search_mesh_name, 
	    &search_variable_name,
	    &search_assocType,
	    &target_assocType,
	    &replace_var,
	    &mod_extension);


  //Setup the library and load in the dataset
  rc = fc_setLibraryVerbosity(verbose_level);        fc_exitIfError(rc);
  rc = fc_initLibrary();                             fc_exitIfError(rc);
  rc = fc_loadDataset(file_in_name, &dataset);       fc_exitIfError(rc);


  //Get the meshes that match a name, or all meshes
  if(!strcmp(search_mesh_name, "")){
    rc = fc_getMeshes(dataset, &numMesh, &meshes);
  } else {
    rc = fc_getMeshByName(dataset, search_mesh_name, &numMesh, &meshes); 
  }
  fc_exitIfError(rc);

  //Go through all the meshes, looking for var's with FC_DT_ELEMENT
  for(i=0; i< numMesh; i++){

    fc_getMeshName(meshes[i], &temp_name);
    //printf("Working on Mesh '%s'\n", temp_name);
    free(temp_name);

    //-- Check All variables --------------------------------------------------
    rc = fc_getVariables(meshes[i], &numVar, &vars); 
    fc_exitIfError(rc);
    
    for(j=0;j<numVar;j++){
      if(checkVarAndCreateName(vars[j], 
			       search_assocType, search_variable_name, 
			       replace_var, mod_extension, 
			       &new_name)){
	FC_Variable new_var;
	char *old_name;

	rc = fc_getVariableName(vars[j], &old_name);
	fc_exitIfError(rc);
	  
	//printf("  converting var %s to %s\n",old_name, (new_name)?new_name : "(library's choice)");

	rc = fc_copyVariableWithNewAssociation( vars[j],  
						target_assocType, new_name, 
						&new_var);
	fc_exitIfError(rc);

	if(replace_var) 
	  fc_deleteVariable(vars[j]);

	//Cleanup
	free(old_name);
	free(new_name);
      }
    }
    free(vars);


    //-- Check All SeqVariables -----------------------------------------------
    rc = fc_getSeqVariables(meshes[i], 
			    &numSeqVar, &numStepPerSeqVar,
			    &seqVars);
    fc_exitIfError(rc);

    for(j=0; j<numSeqVar; j++){

      if(checkVarAndCreateName(seqVars[j][0],
			       search_assocType, search_variable_name, 
			       replace_var, mod_extension, 
			       &new_name)){
	FC_Variable* new_seqVar;
	char *old_name;

	rc = fc_getVariableName(seqVars[j][0], &old_name);
	fc_exitIfError(rc)
	
	  //printf("  converting seqVar %s to %s\n",old_name, (new_name)?new_name : "(library's choice)");


	rc = fc_copySeqVariableWithNewAssociation(numStepPerSeqVar[j], seqVars[j],
						  target_assocType,
						  new_name, //or null
						  &new_seqVar);

	fc_exitIfError(rc);

	if(replace_var) 
	  fc_deleteSeqVariable(numStepPerSeqVar[j], seqVars[j]);

	//Cleanup
	free(old_name);
	free(new_name);
	free(new_seqVar);
      } 
      free(seqVars[j]);
    }
    free(numStepPerSeqVar);
    free(seqVars);

  }
  free(meshes);
 
  //Write out the results
  rc = fc_rewriteDataset(dataset, file_out_name, new_fileType);
 

  //Clean out everything
  fc_finalLibrary();
  exit(0);

}

