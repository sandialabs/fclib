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
 * \file sierra.c
 * \brief Implementation of the \ref Sierra module.
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/modules/sierra.c,v $
 * $Revision: 1.63 $ 
 * $Date: 2006/09/19 00:57:58 $
 */

// C library includes
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

// fc library dependencies
#include "base.h"
#include "storage.h"
#include "library.h"
#include "mesh.h"
#include "subset.h"
#include "variable.h"
#include "meshP.h"
#include "util.h"
#include "geom.h"
#include "geomP.h"
#include "varmathP.h"

// this module
#include "sierra.h"
#include "sierraP.h"

/**
 * \addtogroup Sierra
 * \brief Data structures and routines for dealing with Sierra generated data
 *
 * \description
 *
 *   Data structures and routines for interacting with Sierra generated 
 *   data. (currently focused on Presto spot welds.)
 *
 * \modifications
 *    - 2003-SEP-23  W Koegler  Created.
 */

/**
 *\ingroup Sierra
 *\defgroup PrivateSierra (Private)
 */

/** 
 * \ingroup Sierra
 * \defgroup PrivateSierraHelpers Helpers (not available outside this module)
 *
 * \description
 *
 *   These routines help parse the Sierra input file.
 *
 *   The Sierra input file is generally organized as nested blocks
 *   bounded by "BEGIN" and "END" strings. The helper parse routines
 *   parse specific blocks (e.g. _fc_parseSierra for Sierra blocks,
 *   _fc_parseSpotWeld for spot weld blocks, etc.). 
 *   Their job is to 
 *   1) keep track of 'begin' and 'end's so they know when to exit,
 *   2) copy information specific to the block to the FC_SierraInfo struct,
 *   3) call appropriate parser functions for encountered blocks.
 *
 *   If creating a new helper, it is recommended to copy an existing one
 *   to get all functions it is supposed to do.
 *
 *   If one helper needs fixing up or changing for general capabilities,
 *   you should fix up ALL helpers.
 *
 *   Helper functions expect to start with the file handle at a point where 
 *   it is pointing at the block's user specified name.
 *
 *   Helper functions consume one more 'end' than 'begin' (they're caller 
 *   saw their 'begin', but they see their own 'end'), so the caller of a
 *   helper function needs to decrement it's 'end' count when the
 *   helper returns.
 *
 * \todo The following Fixes apply to all helper parse routines, and perhaps
 *   also fc_parseSierraInputFile.
 *   - If a read in string doesn't match any known thing, discard the rest 
 *     of the line? Could speed things up but probably not necessary.
 *   - Routines can't deal with continuation lines (end with "/$" or "/#")
 *   - Add a generic read keyword value pair to centralize the code for 
 *     dealing with the variations 'keyword = value', 'keyword= value', 
 *     'keyword =value', and 'keyword==value' (note the spaces around '='). 
 *      BUT there is at least one case when value is more than one word!!! 
 *     (function = piecewise linear)
 *   - We don't have any early exits, (i.e if an error) & general error 
 *     checking
 */

/**
 * \ingroup PrivateSierra
 * \brief Global (to sierra.c only) char buffer.
 */
static char string[1024]; // FIX? more or limit scanf?
/**
 * \ingroup PrivateSierra
 * \brief Global (to sierra.c only) char buffer.
 */
static char temp_string[1024];

// These helper parse routines are not available outside of this module
// (not declared in any .h files).
// They should be handed off so next string to read is name of the block

// routines that will call other parse routines
FC_ReturnCode _fc_parseSierra(FILE* inputfile, FC_SierraInfo* sierra);
FC_ReturnCode _fc_parsePrestoProcedure(FILE* inputfile, FC_SierraInfo* sierra);
FC_ReturnCode _fc_parsePrestoRegion(FILE* inputfile, FC_SierraInfo* sierra);

// "leaf" routines (i.e. won't call any more parse routines);
FC_ReturnCode _fc_parseDefinitionForFunction(FILE* inputfile, FC_SierraInfo* sierra);
FC_ReturnCode _fc_parseResultsOutput(FILE* inputfile, FC_SierraInfo* sierra);
FC_ReturnCode _fc_parseSpotWeld(FILE* inputfile, FC_SierraInfo* sierra);

// general use routines
// FIX? reading '=' pairs is the same everywhere and error prone to repeat 
//FC_ReturnCode _fcreadKeywordValue(FILE* inputfile, char* keyword, char* current_string, char** value);
// FIX? replace scanf with scanf_skip_comments ?
// FIX? check scanf's that they really scanned something?

/**
 * \ingroup  PrivateSierra
 * \brief Initialize an FC_SierraInfo struct.
 *
 * \description
 *  
 *   Initialize an already allocated FC_SierraInfo struct with default
 *   values. 
 *
 * \modifications  
 *   - 2003-JUN-24  W Koegler  Created.
 *   - 12/23/04 WSK made private now that fc_parseSierraInputFile creates
 *      creates the object.
 */
FC_ReturnCode _fc_initSierraInfo(
  FC_SierraInfo* sierra /**< input - sierra info object */
) 
{
  // check input
  if (!sierra) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }

  // initialize
  strcpy(sierra->name, "");
  strcpy(sierra->displ_alias, "displacement");
  strcpy(sierra->f_ext_alias, "force_external");
  sierra->numFunction = 0;
  sierra->functions = NULL;
  sierra->numSpotWeld = 0;
  sierra->spotwelds = NULL;

  return FC_SUCCESS;
}

/**
 * \ingroup  Sierra
 * \brief Free memory associated with a FC_SierraInfo struct.
 *
 * \description
 *  
 *    This frees allocated members within the struct as well as
 *    the struct itself.
 *
 * \modifications  
 *   12/23/04 WSK created. Replaces fc_clearSierraInfo().
 */
void fc_freeSierraInfo(
  FC_SierraInfo* sierra /**< input - sierra info object */
) 
{
  int i;

  // check input
  if (!sierra)
    return;

  // log message
  fc_printfLogMessage("Freeing sierra info.");

  // release dynamically allocated members (and their members ...)
  if (sierra->numFunction > 0) {
    for (i = 0; i < sierra->numFunction; i++) {
      free(sierra->functions[i]->values);
      free(sierra->functions[i]);
    }
    free(sierra->functions);
  }
  if (sierra->numSpotWeld > 0) {
    for (i = 0; i < sierra->numSpotWeld; i++) 
      free(sierra->spotwelds[i]);
    free(sierra->spotwelds);
  }

  // free this
  free(sierra);
}

/**
 * \ingroup  Sierra
 * \brief Print an FC_SierraInfo struct.
 *
 * \description
 *  
 *   Print the contents of an FC_SierraInfo struct.
 *
 * \modifications  
 *   2003-JUN-24  W Koegler  Created.
 */
FC_ReturnCode fc_printSierraInfo(
  FC_SierraInfo* sierra /**< input - sierra info object */
) {
  int i, j;

  // check input
  if (!sierra) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
 
  printf("Sierra Problem '%s'\n", sierra->name);
  fflush(NULL);

  // variables
  printf("  Alias for 'displacement' = '%s'\n", sierra->displ_alias);
  printf("  Alias for 'force_external' = '%s'\n", sierra->f_ext_alias);
  fflush(NULL);

  // functions
  printf("  Found %d functions:\n", sierra->numFunction);
  for (i = 0; i < sierra->numFunction; i++) {
    FC_SierraFuncInfo *function = sierra->functions[i];
    printf("    %d: '%s'\n", i, function->name);
    printf("      type = '%s'\n", function->type);
    printf("      abscissa label = '%s'\n", function->abscissa_label);
    printf("      ordinate label = '%s'\n", function->ordinate_label);
    printf("      number of values = %d\n", function->numValue);
    printf("      values = %f", function->values[0]);
    if (function->numValue > 1)
      printf("  %f\n",function->values[1]);
    else
      printf("\n");
    if (function->numValue > 2) { 
      for (j = 1; j < function->numValue/2; j++)
        printf("               %f  %f\n", function->values[j*2], 
               function->values[j*2+1]);
    }
    printf("      min abscissa = %f\n", function->abscissa_min);
    printf("      max abscissa = %f\n", function->abscissa_max);
    if (function->numValue > 1) {
      printf("      min ordinate = %f\n", function->ordinate_min);
      printf("      max ordinate = %f\n", function->ordinate_max);
    }
    else {
      printf("      min ordinate = (not applicable)\n");
      printf("      max ordinate = (not applicable)\n");
    }
  }
  fflush(NULL);

  // spot welds
  printf("  Found %d spot welds:\n", sierra->numSpotWeld);
  for (i = 0; i < sierra->numSpotWeld; i++) {
    FC_SierraSpotWeldInfo *spotweld = sierra->spotwelds[i];
    printf("    %d: '%s'\n", i, spotweld->name);
    printf("      node set = '%s' (name = '%s')\n", spotweld->nodeset,
           spotweld->nodeset_name);
    printf("      surface  = '%s'  (name = '%s')\n", spotweld->surface,
           spotweld->surface_name);
    printf("      normal displ. function (ID = %d)     = '%s'\n", 
           spotweld->normalId, spotweld->normal);
    printf("      normal displ. scale factor = %f\n", spotweld->normal_scale);
    printf("      tangential displ. function (ID = %d) = '%s'\n", 
           spotweld->tangentId, spotweld->tangent);
    printf("      tangent displ. scale factor = %f\n", spotweld->tangent_scale);
    printf("      failure envelope exponent = %f\n", spotweld->exponent);
    printf("      failure decay cycles = %d\n", spotweld->decay_cycles);
  }
  fflush(NULL);

  return FC_SUCCESS;
}

/**
 * \ingroup  Sierra
 * \brief Parse a Sierra Input File and put info into FC_SierraInfo struct.
 *
 * \description
 *  
 *   Parses a Sierra input file and returns a new FC_SierraInfo object
 *   containibg the sierra info. It's memory can be released using
 *   fc_freeSierraInfo().
 *
 *   Currently the parser only returns the information needed to analyze
 *   Presto spot welds. However, the parsing routines were written to
 *   (hopefully) be easily expanded when additional info is required.
 *   See the private documentation for _fc_parseSierra for more
 *   details about the parsing process.
 *
 * \modifications  
 *    - 2003-JUN-24  W Koegler  Created.
 *    - 12/23/04 WSK Changed so it creates the FC_SierraInfo object instead
 *      of taking it as an object.
 */
FC_ReturnCode fc_parseSierraInputFile(
  char* filename,         /**< input - the sierra input file */
  FC_SierraInfo** sierra  /**< output - sierra info */ 
) 
{
  int i, j;
  FC_ReturnCode rc, rc_keep;
  int level = 0;
  FILE* inputfile;

  // default return
  if (sierra)
    *sierra = NULL;

  // check input
  if (!filename || !sierra) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(rc));
    return FC_INPUT_ERROR;
  }

  // open the input file
  inputfile = fopen(filename, "r");
  if (!inputfile) {
        fc_printfErrorMessage("sierra input parser could not open file '%s'\n", 
                 filename);
    return FC_FILE_IO_ERROR;
  }

  // setup sierra
  *sierra = malloc(sizeof(FC_SierraInfo));
  rc = _fc_initSierraInfo(*sierra);
  if (rc != FC_SUCCESS)
    return rc;

  // Loop over lines in the file
  while (fscanf(inputfile, "%s", string) > 0) { 
    
    // skip comment line 
    if (string[0] == '#' || string[0] == '$') 
      fscanf(inputfile, "%[^\n]%1[\n]", string, temp_string);
    
    // end block 
    else if (!strcasecmp(string, "end")) {
      level--;
      fscanf(inputfile, "%[^\n]%1[\n]", string, temp_string);
    }
    
    // begin block 0
    else if (!strcasecmp(string, "begin")) {
      level++;

      // determine type of begin
      fscanf(inputfile, "%s", string);
      
      // sierra block
      if (!strcasecmp(string, "sierra")) {
        rc = _fc_parseSierra(inputfile, *sierra);
        if (rc != FC_SUCCESS)
          return rc;
        level--;
      }
    }
  }
  
  // close the input file
  fclose(inputfile);

  if (level != 0) 
    fc_printfWarningMessage("parser: begins & ends do not match.");

  // Make extra info for spot welds
  rc_keep = FC_SUCCESS; // keep track of any errors
  for (i = 0; i < (*sierra)->numSpotWeld; i++) {
    FC_SierraSpotWeldInfo* spotweld = (*sierra)->spotwelds[i];

    // make connections to the functions
    spotweld->normalId = -1;
    spotweld->tangentId = -1;
    spotweld->normal_p = NULL;
    spotweld->tangent_p = NULL;
    for (j = 0; j < (*sierra)->numFunction; j++) {
      FC_SierraFuncInfo* function = (*sierra)->functions[j];
      if (!strcasecmp(spotweld->normal, function->name)) {
        spotweld->normalId = j;
        spotweld->normal_p = function;
      }
      if (!strcasecmp(spotweld->tangent, function->name)) {
        spotweld->tangentId = j;
        spotweld->tangent_p = function;
      }
    }
    if (spotweld->normalId == -1 || spotweld->tangentId == -1) {      
      rc_keep = FC_ERROR; // don't return, try to set up other spotwelds
      fc_printfWarningMessage("function for spotweld '%s' not found\n", 
              spotweld->name);
    }

    // make names used in database (e.g. 'nodelist_12' -> 'Node set ID 12'
    // and 'surface_3002' -> 'Side set ID 3003')
    sprintf(spotweld->nodeset_name, "nodelist_%s", &(spotweld->nodeset[9]));
    sprintf(spotweld->surface_name, "surface_hex8_quadface4_%s", &(spotweld->surface[8]));
  }
  
  return rc_keep;
}

/**
 * \ingroup  PrivateSierraHelpers
 * \brief Parse a Sierra block.
 *
 * \description
 *
 *   This routines parses a Sierra block. It calls helper parsers
 *   for function definition blocks and the presto procedure block.
 *
 *   What follows is applicable to all helper parse routines.
 * 
 *   The Sierra input file is generally organized as nested blocks
 *   bounded by "BEGIN" and "END" strings. The helper parse routines
 *   parse specific blocks (e.g. _fc_parseSierra for Sierra blocks,
 *   _fc_parseSpotWeld for spot weld blocks, etc.). 
 *   Their job is to 
 *   1) keep track of 'begin' and 'end's so they know when to exit,
 *   2) copy information specific to the block to the FC_SierraInfo struct,
 *   3) call appropriate parser functions for encountered blocks.
 *
 *   If creating a new helper, it is recommended to copy an existing one
 *   to get all functions it is supposed to do.
 *
 *   If one helper needs fixing up or changing for general capabilities,
 *   you should fix up ALL helpers.
 *
 *   Helper functions expect to start with the file handle at a point where 
 *   it is pointing at the block's user specified name.
 *
 *   Helper functions consume one more 'end' than 'begin' (they're caller 
 *   saw their 'begin', but they see their own 'end'), so the caller of a
 *   helper function needs to decrement it's 'end' count when the
 *   helper returns.
 *
 * \modifications  
 *   2003-JUN-24  W Koegler  Created.
 */
FC_ReturnCode _fc_parseSierra(FILE* inputfile, FC_SierraInfo* sierra) {
  FC_ReturnCode rc; // return code
  int level = 1;
  int parent_level = level - 1;

  // get block name
  fscanf(inputfile, "%s", sierra->name);

  // loop over lines until matching end found
  while (level > parent_level  &&  fscanf(inputfile, "%s", string) > 0) {
    
    // skip comment line 
    if (string [0] == '#' || string[0] == '$')
      fscanf(inputfile, "%[^\n]%1[\n]", string, temp_string);
    
    // end block
    else if (!strcasecmp(string, "end")) {
      level--;
      fscanf(inputfile, "%[^\n]%1[\n]", string, temp_string);
    }

    // begin block 
    else if (!strcasecmp(string, "begin")) {
      level++;

      // get type of begin
      fscanf(inputfile, "%s", string);
      
      // definition block 1
      if (!strcasecmp(string, "definition")) {
        fscanf(inputfile, "%s", string);
        if (!strcasecmp(string, "for")) {
          fscanf(inputfile, "%s", string);
          if (!strcasecmp(string, "function")) {
            rc = _fc_parseDefinitionForFunction(inputfile, sierra);
            if (rc != FC_SUCCESS)
              return rc;
            level--;
          }
        }
      } // definition block 1
            
      // procedure block
      else if (!strcasecmp(string, "presto")) {
        fscanf(inputfile, "%s", string);
        if (!strcasecmp(string, "procedure")) {
          rc = _fc_parsePrestoProcedure(inputfile, sierra);
          if (rc != FC_SUCCESS)
            return rc;
          level--;
        }
      } // procedure block
    } // begin block
  } // while loop

  return FC_SUCCESS;
}

/**
 * \ingroup  PrivateSierraHelpers
 * \brief Parse a function definition block
 *
 * \description
 *
 *   This routines parses a function definition block. It creates
 *   a FC_SierraFuncInfo struct and adds it to the FC_SierraInfo struct's
 *   array of functions. It reads in the sierra info and also calculates the
 *   min and max values of the ordinate and abscissa. It does not call any
 *   other helper parse routines.
 *
 *   See _fc_parseSierra for generic information about all helper parse 
 *   routines.
 *
 * \modifications  
 *   2003-JUN-24  W Koegler  Created.
 *
 *   2003-AUG-20  W Koegler  Fixed bug with comments within values block
 */
FC_ReturnCode _fc_parseDefinitionForFunction(FILE* inputfile, FC_SierraInfo* sierra) {
  int level = 1;
  int parent_level = level - 1;
  int numFunction;
  FC_SierraFuncInfo* function;
  void* temp;
  
  // make room for new function
  numFunction = sierra->numFunction;
  temp = realloc(sierra->functions, sizeof(FC_SierraFuncInfo*)*(numFunction + 1));
  if (temp == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  sierra->functions = temp;
  sierra->functions[numFunction] = malloc(sizeof(FC_SierraFuncInfo));
  if (sierra->functions[numFunction] == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  sierra->numFunction++;

  // initialize function
  function = sierra->functions[sierra->numFunction-1];
  strcpy(function->abscissa_label, "");
  strcpy(function->ordinate_label, "");
  function->numValue = 0;
  function->values = NULL;

  // get block name
  fscanf(inputfile, "%s", function->name);
  
  // loop over lines until matching end found
  while (level > parent_level  && fscanf(inputfile, "%s", string) > 0) {

    // skip comment line
    if (string[0] == '#' || string[0] == '$')
      fscanf(inputfile, "%[^\n]%1[\n]", string, temp_string);

    // end
    else if (!strcasecmp(string, "end")) {
      level--;
      fscanf(inputfile, "%[^\n]%1[\n]", string, temp_string);
    }

    // type
    else if (!strncasecmp(string, "type", 4)) {
      int len = strlen(string);
      if (len == 4) {
        fscanf(inputfile, "%s", string);
        if (string[0] == '=' && strlen(string) > 1) // probably '=type'
          strcpy(function->type, &string[1]);
        else
          fscanf(inputfile, "%s", function->type);
      }
      else if (len == 5) // probably was 'type='
        fscanf(inputfile, "%s", function->type);
      else if (len > 5) // probably was 'type=type'
        strcpy(function->type, &string[5]);
      if (!strcasecmp(function->type, "piecewise")) {
        fscanf(inputfile, "%s", string);
        sprintf(&function->type[strlen(function->type)], " %s", string);
      }
    } // end type

    // abscissa label
    else if (!strncasecmp(string, "abscissa", 8)) {
      int len = strlen(string);
      if (len == 8) {
        fscanf(inputfile, "%s", string);
        if (string[0] == '=' && strlen(string) > 1) // probably '=name'
          strcpy(function->abscissa_label, &string[1]);
        else
          fscanf(inputfile, "%s", function->abscissa_label);
      }
      else if (len == 9) // probably was 'abscissa='
        fscanf(inputfile, "%s", function->abscissa_label);
      else if (len > 9) // probably was 'abscissa=name'
        strcpy(function->abscissa_label, &string[9]);
    }  // end abscissa label

    // ordinate label
    else if (!strncasecmp(string, "ordinate", 8)) {
      int len = strlen(string);
      if (len == 8) {
        fscanf(inputfile, "%s", string);
        if (string[0] == '=' && strlen(string) > 1) // probably '=name'
          strcpy(function->ordinate_label, &string[1]);
        else
          fscanf(inputfile, "%s", function->ordinate_label);
      }
      else if (len == 9) // probably was 'ordinate='
        fscanf(inputfile, "%s", function->ordinate_label);
      else if (len > 9) // probably was 'ordinate=name'
        strcpy(function->ordinate_label, &string[9]);
    }  // end ordinate label

    // begin block
    else if (!strcasecmp(string, "begin")) {
      level++;

      // get type of begin
      fscanf(inputfile, "%s", string);
    
      // values block
      if (!strcasecmp(string, "values")) {
        while (fscanf(inputfile, "%s", string) > 0) {

          // skip comment line
          if (string[0] == '#' || string[0] == '$')
            fscanf(inputfile, "%[^\n]%1[\n]", string, temp_string);

          // end;
          else if (!strcasecmp(string, "end")) {
            level--;
            break;
          }

          // read the values and save the max
          else {
            double flt = atof(string);
            temp = realloc(function->values, 
                                 sizeof(double)*(function->numValue + 1));
            if (temp == NULL) {
              fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
              return FC_MEMORY_ERROR;
            }
            function->values = temp;
            function->values[function->numValue] = flt;
            function->numValue++;

            if (function->numValue%2 == 1) { // odd = abscissa (x)
              if (function->numValue == 1) {
                function->abscissa_min = flt;
                function->abscissa_max = flt;
              }
              else {
                if (flt < function->abscissa_min)
                  function->abscissa_min = flt;
                else if (flt > function->abscissa_max)
                  function->abscissa_max = flt;
              }
            }
            else {// even's
              if (function->numValue == 2) {
                function->ordinate_min = flt;
                function->ordinate_max = flt;
              }
              else {
                if (flt < function->ordinate_min)
                  function->ordinate_min = flt;
                else if (flt > function->ordinate_max)
                  function->ordinate_max = flt;
              } 
            }
          }
        }
      } // values block
    } // begin block
  } // while loop

  return FC_SUCCESS;
}

/**
 * \ingroup  PrivateSierraHelpers
 * \brief Parse a Presto procedure block
 *
 * \description
 *
 *   This routines parses a Presto procedure block. It calls a helper
 *   parser for the Presto region block.
 *
 *   See _fc_parseSierra for generic information about all helper parse 
 *   routines.
 *
 * \modifications  
 *   2003-JUN-24  W Koegler  Created.
 */
FC_ReturnCode _fc_parsePrestoProcedure(FILE* inputfile, FC_SierraInfo* sierra) {
  FC_ReturnCode rc; // return code
  int level = 1;
  int parent_level = level - 1;

  // get block name
  fscanf(inputfile, "%s", string); // discard

  // loop over lines until matching end found
  while (level > parent_level  &&  fscanf(inputfile, "%s", string) > 0) {
    
    // skip comment line 
    if (string [0] == '#' || string[0] == '$')
      fscanf(inputfile, "%[^\n]%1[\n]", string, temp_string);
    
    // end block 1
    else if (!strcasecmp(string, "end")) {
      level--;
      fscanf(inputfile, "%[^\n]%1[\n]", string, temp_string);
    }

    // begin block 
    else if (!strcasecmp(string, "begin")) {
      level++;

      // get type of begin
      fscanf(inputfile, "%s", string);
      
      // region block
      if (!strcasecmp(string, "presto")) {
        fscanf(inputfile, "%s", string);
        if (!strcasecmp(string, "region")) {
          rc = _fc_parsePrestoRegion(inputfile, sierra);
          if (rc != FC_SUCCESS)
            return rc;
          level--;
        }
      } // region block
    } // begin block
  } // while loop

  return FC_SUCCESS;
}

/**
 * \ingroup  PrivateSierraHelpers
 * \brief Parse a Presto region block
 *
 * \description
 *
 *   This routines parses a Presto region block. It calls helper
 *   parsers for the results block and spot weld blocks.
 *
 *   See _fc_parseSierra for generic information about all helper parse 
 *   routines.
 *
 * \modifications  
 *   2003-JUN-24  W Koegler  Created.
 */
FC_ReturnCode _fc_parsePrestoRegion(FILE* inputfile, FC_SierraInfo* sierra) {
  FC_ReturnCode rc; // return code
  int level = 1;
  int parent_level = 0;

  // get block name
  fscanf(inputfile, "%s", string); // discard

  // loop over lines until matching end found
  while (level > parent_level  &&  fscanf(inputfile, "%s", string) > 0) {
    
    // skip comment line 
    if (string [0] == '#' || string[0] == '$')
      fscanf(inputfile, "%[^\n]%1[\n]", string, temp_string);
    
    // end block 1
    else if (!strcasecmp(string, "end")) {
      level--;
      fscanf(inputfile, "%[^\n]%1[\n]", string, temp_string);
    }

    // begin block 
    else if (!strcasecmp(string, "begin")) {
      level++;

      // get type of begin
      fscanf(inputfile, "%s", string);
      
     // results block
      if (!strcasecmp(string, "results")) {
        fscanf(inputfile, "%s", string);
        if (!strcasecmp(string, "output")) {
          rc = _fc_parseResultsOutput(inputfile, sierra);
          if (rc != FC_SUCCESS)
            return rc;
          level--;
        }
      } // results block

      // spot weld block
      else if (!strcasecmp(string, "spot")) {
        fscanf(inputfile, "%s", string);
        if (!strcasecmp(string, "weld")) {
          rc = _fc_parseSpotWeld(inputfile, sierra);
          if (rc != FC_SUCCESS)
            return rc;
          level--;
        }
      } // spot weld block
    } // begin block
  } // while loop

  return FC_SUCCESS;
}

/**
 * \ingroup  PrivateSierraHelpers
 * \brief Parse a results output block
 *
 * \description
 *
 *   This routines parses a results output block. It saves the
 *   aliases for the variables "displacements" and "force_external". 
 *   It does not call any other helper parse routines.
 *
 *   See _fc_parseSierra for generic information about all helper parse 
 *   routines.
 *
 * \modifications  
 *   2003-JUN-24  W Koegler  Created.
 */
FC_ReturnCode _fc_parseResultsOutput(FILE* inputfile, FC_SierraInfo* sierra) {
  int level = 1;
  int parent_level = level - 1;

  // get block name
  fscanf(inputfile, "%s", string);
  
  // loop over lines until matching end found
  while (level > parent_level  && fscanf(inputfile, "%s", string) > 0) {

    startOfWhile: 

    // skip comment line
    if (string[0] == '#' || string[0] == '$')
      fscanf(inputfile, "%[^\n]%1[\n]", string, temp_string);

    // end
    else if (!strcasecmp(string, "end")) {
      level--;
      fscanf(inputfile, "%[^\n]%1[\n]", string, temp_string);
    }

    // begin block
    else if (!strcasecmp(string, "begin")) {
      level++;
    }

    // FIX instead of grabbing specific vars, collect all (& node/elem etc)
    // grab the alias for displacement & force_external
    else if (!strcasecmp(string, "nodal") || !strcasecmp(string, "element") 
             || !strcasecmp(string, "global")) {
      fscanf(inputfile, "%s", string);
      if (!strncasecmp(string, "variables", 9)) {
        // get the variable name
        int len = strlen(string);
        if (len == 9) {
          fscanf(inputfile, "%s", string);
          if (string[0] == '='  &&  strlen(string) > 1) // '=variable_name'
            strcpy(string, &string[1]);
          else
            fscanf(inputfile, "%s", string);
        }
        else if (len == 10) // we probably scanned in 'variables='
          fscanf(inputfile, "%s", string);
        else if (len > 10)  // we probably scanned in 'variables=variable_name'
          strcpy(string, &string[10]);

        // if the next string is not 'as' we've read in the next line!
        fscanf(inputfile, "%s", temp_string);
        if (strcasecmp(temp_string, "as") && level > parent_level) {
          strcpy(string, temp_string);
          goto startOfWhile;  // start another loop without a new scanf
        }

        // get the alias name
        if (!strcasecmp(string, "force_external"))
          fscanf(inputfile, "%s", sierra->f_ext_alias);
        else if (!strcasecmp(string, "displacement"))
          fscanf(inputfile, "%s", sierra->displ_alias);
      }
    } // force external
  } // while loop

  return FC_SUCCESS;
}

/**
 * \ingroup  PrivateSierraHelpers
 * \brief Parse a spot weld block
 *
 * \description
 *
 *   This routines parses a spot weld block. It creates a FC_SierraSpotWeldInfo
 *   struct and adds it to the FC_SierraInfo struct's array of spot welds.
 *   It reads int the sierra info (the sierra parser will create the
 *   function indices after all spot welds and functions have been parsed).
 *   It does not call any other helper parse routines.
 *
 *   See _fc_parseSierra for generic information about all helper parse 
 *   routines.
 *
 * \modifications  
 *   2003-JUN-24  W Koegler  Created.
 */
FC_ReturnCode _fc_parseSpotWeld(FILE* inputfile, FC_SierraInfo* sierra) {
  int level = 1;
  int parent_level = level - 1;
  int numSpotWeld;
  FC_SierraSpotWeldInfo* spotweld;
  void *temp;

  // make room for new spotweld
  numSpotWeld = sierra->numSpotWeld;
  temp = realloc(sierra->spotwelds, 
                 sizeof(FC_SierraSpotWeldInfo*)*(numSpotWeld + 1));
  if (temp == NULL)
    return FC_MEMORY_ERROR;
  sierra->spotwelds = temp;
  sierra->spotwelds[numSpotWeld] = malloc(sizeof(FC_SierraSpotWeldInfo));
  if (sierra->spotwelds[numSpotWeld] == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  sierra->numSpotWeld++;

  // initialize new spotweld
  spotweld = sierra->spotwelds[sierra->numSpotWeld-1];
  spotweld->normal_scale = 1.0F;
  spotweld->tangent_scale = 1.0F;
  spotweld->nodeset_name[0] = '\0';
  spotweld->surface_name[0] = '\0';

  // get block name
  fscanf(inputfile, "%s", spotweld->name);
  
  // loop over lines until matching end found
  while (level > parent_level  && fscanf(inputfile, "%s", string) > 0) {

    // skip comment line
    if (string[0] == '#' || string[0] == '$')
      fscanf(inputfile, "%[^\n]%1[\n]", string, temp_string);

    // end
    else if (!strcasecmp(string, "end")) {
      level--;
      fscanf(inputfile, "%[^\n]%1[\n]", string, temp_string);
    }

    // begin block
    else if (!strcasecmp(string, "begin")) {
      level++;
    }

    // node set name
    else if (!strcasecmp(string, "node")) {
      fscanf(inputfile, "%s", string);
      if (!strncasecmp(string, "set", 3)) {
        int len = strlen(string);
        if (len == 3) {
          fscanf(inputfile, "%s", string);
          if (string[0] == '='  &&  strlen(string) > 1) // probably '=name'
            strcpy(spotweld->nodeset, &string[1]);
          else
            fscanf(inputfile, "%s", spotweld->nodeset);
        }
        else if (len == 4) // probably was 'set='
          fscanf(inputfile, "%s", spotweld->nodeset);
        else if (len > 4) // probably was 'set=name' 
          strcpy(spotweld->nodeset, &string[4]);
      }
    }

    // surface name
    else if (!strncasecmp(string, "surface", 7)) {
      int len = strlen(string);
      if (len == 7) {
        fscanf(inputfile, "%s", string);
        if (string[0] == '=' && strlen(string) > 1) // probably '=name'
          strcpy(spotweld->surface, &string[1]);
        else
          fscanf(inputfile, "%s", spotweld->surface);
      }
      else if (len == 8) // probably was 'surface='
        fscanf(inputfile, "%s", spotweld->surface);
      else if (len > 8) // probably was 'surface=name'
        strcpy(spotweld->surface, &string[8]);
    }

    // normal displacement function name or scale factor
    else if (!strcasecmp(string, "normal")) {
      fscanf(inputfile, "%s", string);
      if (!strcasecmp(string, "displacement")) {
        fscanf(inputfile, "%s", string);

        // function name
        if (!strncasecmp(string, "function", 8)) {
          int len = strlen(string);
          if (len == 8) {
            fscanf(inputfile, "%s", string);
            if (string[0] == '=' && strlen(string) > 1) // probably '=name'
              strcpy(spotweld->normal, &string[1]);
            else
              fscanf(inputfile, "%s", spotweld->normal);
          }
          else if (len == 9) // probably was 'function='
            fscanf(inputfile, "%s", spotweld->normal);
          else if (len > 9) // probably was 'function=name'
            strcpy(spotweld->normal, &string[9]);
        } // end function name

        // scale factor
        else if (!strcasecmp(string, "scale")) {
          fscanf(inputfile, "%s", string);
          if (!strncasecmp(string, "factor", 6)) {
            int len = strlen(string);
            if (len == 6) {
              fscanf(inputfile, "%s", string);
              if (string[0] == '=' && strlen(string) > 1) // probably '=num'
                spotweld->normal_scale = atof(&string[1]);
              else
                fscanf(inputfile, "%lf", &spotweld->normal_scale);
            }
            else if (len == 7) // probably was 'scale='
              fscanf(inputfile, "%lf", &spotweld->normal_scale);
            else if (len > 7) // probably was 'scale=num'
              spotweld->normal_scale = atof(&string[7]);
          } 
        } // end scale factor
      }
    } // normal displacement function name or scale factor

    // tangential displacement function name or scale factor
    else if (!strcasecmp(string, "tangential")) {
      fscanf(inputfile, "%s", string);
      if (!strcasecmp(string, "displacement")) {
        fscanf(inputfile, "%s", string);

        // function name
        if (!strncasecmp(string, "function", 8)) {
          int len = strlen(string);
          if (len == 8) {
            fscanf(inputfile, "%s", string);
            if (string[0] == '=' && strlen(string) > 1) // probably '=name'
              strcpy(spotweld->tangent, &string[1]);
            else
              fscanf(inputfile, "%s", spotweld->tangent);
          }
          else if (len == 9) // probably was 'function='
            fscanf(inputfile, "%s", spotweld->tangent);
          else if (len > 9) // probably was 'function=name'
            strcpy(spotweld->tangent, &string[9]);
        } // end function name

        // scale factor
        else if (!strcasecmp(string, "scale")) {
          fscanf(inputfile, "%s", string);
          if (!strncasecmp(string, "factor", 6)) {
            int len = strlen(string);
            if (len == 6) {
              fscanf(inputfile, "%s", string);
              if (string[0] == '=' && strlen(string) > 1) // probably '=num'
                spotweld->tangent_scale = atof(&string[1]);
              else
                fscanf(inputfile, "%lf", &spotweld->tangent_scale);
            }
            else if (len == 7) // probably was 'scale='
              fscanf(inputfile, "%lf", &spotweld->tangent_scale);
            else if (len > 7) // probably was 'scale=num'
              spotweld->tangent_scale = atof(&string[7]);
          } 
        } // end scale factor
      }
    } // tangential displacement function name or scale factor

    // failure envelope exponent & decay cycles
    else if (!strcasecmp(string, "failure")) {
      fscanf(inputfile, "%s", string);

      // envelope exponent
      if (!strcasecmp(string, "envelope")) {
        fscanf(inputfile, "%s", string);
        if (!strncasecmp(string, "exponent", 8)) {
          int len = strlen(string);
          if (len == 8) {
            fscanf(inputfile, "%s", string);
            if (string[0] == '=' && strlen(string) > 1) // probably '=number'
              spotweld->exponent = atof(&string[1]);
            else
              fscanf(inputfile, "%lf", &spotweld->exponent);
          }
          else if (len == 9) // probably was 'function='
            fscanf(inputfile, "%lf", &spotweld->exponent);
          else if (len > 9) // probably was 'function=name'
            spotweld->exponent = atof(&string[9]);
        }
      } // end envelope exponent

      // decay cycles
      if (!strcasecmp(string, "decay")) {
        fscanf(inputfile, "%s", string);
        if (!strncasecmp(string, "cycles", 6)) {
          int len = strlen(string);
          if (len == 6) {
            fscanf(inputfile, "%s", string);
            if (string[0] == '=' && strlen(string) > 1) // probably '=number'
              spotweld->decay_cycles = atoi(&string[1]);
            else
              fscanf(inputfile, "%d", &spotweld->decay_cycles);
          }
          else if (len == 7) // probably was 'cycles='
            fscanf(inputfile, "%d", &spotweld->decay_cycles);
          else if (len > 7) // probably was 'cycles=num'
            spotweld->decay_cycles = atoi(&string[7]);
        }
      } // end envelope exponent
      
    } // failure envelope exponent

  } // while loop

  return FC_SUCCESS;
}

/**
 * \ingroup  Sierra
 * \brief Compute the failure parameter for a spot weld.
 *
 * \description
 *
 *    The failure parameter is calculated based on the distance between the
 *    node of the node set and the point of attachment on the sideset of the
 *    spotweld. The attachment remains intact as long as ufail remains less
 *    than zero
 *
 *    ufail = (unorm/umaxnorm)^p + (utang/umaxtang)^p 
 *
 *    where unorm and utang are the normal and tangential components of the
 *    distance between the node and side sets, and umaxnorm and umaxtang are
 *    the maximum displacements in the functions defining the spotweld. The
 *    caller is responsible for freeing the created params array.
 *
 *    The index of the time step at which the spot weld is considered to have
 *    failed is also returned. (A value of -1 indicates that the spot weld did
 *    not fail.) This is when ufail >= 1.
 *
 * \modifications 
 *    - 2003-SEP-08  W Koegler  Changed so that it also does the broken
 *      characterization.
 *    - 06/22/2004 WSK Changed to take subsets instead of meshes for
 *      the node and side sets.
 */
FC_ReturnCode fc_evalPrestoSpotWeld(
  FC_SierraSpotWeldInfo* spotweld, /**< input - Sierra spot weld info */
  FC_Subset node_set,  /**< input - node set specifying the spot weld */
  FC_Subset side_set,  /**< input - side set specifying the surface */
  int numStep,             /**< input - the number of time steps */
  FC_Variable* node_displs, /**< input - array of coord displacements on the
                                     node set */ 
  FC_Variable* side_displs, /**< input - array of coord displacements on the 
                                           side set */
  double** ufails_p, /**< output - created array of failure parameters per step */
  int *broken_id     /**< output - the step at which the spot weld fails 
                                (-1 = no failure */
) {
  FC_ReturnCode rc;
  int i, j, k;
  int numMember, nodesetID, sidesetID, *temp_int_array, dim, temp_num;
  int side_numVert, *side_vertexIDs;
  FC_AssociationType assoc;
  FC_ElementType elemtype;
  FC_Mesh node_mesh, side_mesh, temp_mesh;
  double umaxnorm, umaxtang, exponent;
  double unorm, utang; 
  double* ufails;
  int numParam;
  double* params;
  
  // default return values
  if (ufails_p)
    *ufails_p = NULL;
  if (broken_id)
    *broken_id = -1;
  
  // check input
  if (!spotweld || numStep < 1 || !node_displs || !side_displs || 
      !ufails_p || !broken_id) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // get parent meshes for subsets and make sure have same dim
  fc_getMeshFromSubset(node_set, &node_mesh);
  fc_getMeshFromSubset(side_set, &side_mesh);
  fc_getMeshDim(node_mesh, &dim);
  fc_getMeshDim(side_mesh, &temp_num);
  if (dim != temp_num) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }
  
  // make sure node displs are on the same mesh as node set
  for (i = 0; i < numStep; i++) {
    fc_getMeshFromVariable(node_displs[i], &temp_mesh);
    if (!FC_HANDLE_EQUIV(temp_mesh, node_mesh)) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
      return FC_INPUT_ERROR;
    }
  }
  
  // make sure side displaces are on the same mesh as side set
  for (i = 0; i < numStep; i++) {
    fc_getMeshFromVariable(side_displs[i], &temp_mesh);
    if (!FC_HANDLE_EQUIV(temp_mesh, side_mesh)) {
      fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
      return FC_INPUT_ERROR;
    }
  }
  
  // Get member ID of node set and make sure is o.k. node set
  rc = fc_getSubsetMembersAsArray(node_set, &numMember, &temp_int_array);
  if (rc != FC_SUCCESS)
    return rc;
  if (numMember != 1)
    return FC_ERROR;
  nodesetID = temp_int_array[0];
  free(temp_int_array);
  rc = fc_getSubsetInfo(node_set, NULL, NULL, &assoc);
  if (rc != FC_SUCCESS)
    return rc;
  if (assoc != FC_AT_VERTEX)
    return FC_ERROR;
  
  // Get member ID and vertexIDs of side set
  rc = fc_getSubsetMembersAsArray(side_set, &numMember, &temp_int_array);
  if (rc != FC_SUCCESS)
    return rc;
  if (numMember != 1)
    return FC_ERROR;
  sidesetID = temp_int_array[0];
  free(temp_int_array);
  rc = fc_getSubsetInfo(side_set, NULL, NULL, &assoc);
  if (rc != FC_SUCCESS)
    return rc;
  if (assoc == FC_AT_ELEMENT) {
    int* conns;
    rc = fc_getMeshElementConnsPtr(side_mesh, &conns);
    if (rc != FC_SUCCESS)
      return rc;
    fc_getMeshElementType(side_mesh, &elemtype);
    side_numVert = fc_getElementTypeNumVertex(elemtype);
    side_vertexIDs = &conns[sidesetID*side_numVert]; // don't free
  }
  else if (assoc == FC_AT_FACE) {
    FC_ElementType* side_faceTypes;
    int maxNumVertPerFace, *numVertPerFace, *faceConns;
    rc = fc_getMeshFaceTypesPtr(side_mesh, &side_faceTypes);
    if (rc != FC_SUCCESS)
      return rc;
    rc = fc_getMeshFaceConnsPtr(side_mesh, &numVertPerFace, &maxNumVertPerFace,
                                &faceConns);
    if (rc != FC_SUCCESS)
      return rc;
    elemtype = side_faceTypes[sidesetID];
    side_numVert = numVertPerFace[sidesetID];
    side_vertexIDs = &faceConns[sidesetID*maxNumVertPerFace]; // don't free
  }
  else {
    fc_printfErrorMessage("Sorry, fc_evalPrestoSpotWeld does not "
                          "handle sidesets with association %s yet\n",
                          fc_getAssociationTypeText(assoc));
    return FC_ERROR;
  }
  // FIX take this out when this routine can do other things
  if (elemtype != FC_ET_QUAD) {
    fc_printfErrorMessage("Sorry, fc_evalPrestoSpotWeld only does "
                          "quad sidesets at the moment\n");
    return FC_ERROR;
  }
  
  // get critical displs from spot weld info
  umaxnorm = spotweld->normal_p->abscissa_max;
  umaxtang = spotweld->tangent_p->abscissa_max;
  exponent = spotweld->exponent;
  
  // determine the params to locate attachment site in side set
  // FIX: make this work for a tri too!
  {
    FC_Coords side_vertCoords[4]; // array of side vertex coords, max 4 for quads
    FC_Coords attach_orig; // coordinates of attachment on side set at t0
    double* node_mesh_coords;
    double* side_mesh_coords;
   
    rc = fc_getMeshCoordsPtr(node_mesh, &node_mesh_coords);
    if (rc != FC_SUCCESS)
      return rc;
    rc = fc_getMeshCoordsPtr(side_mesh, &side_mesh_coords);
    if (rc != FC_SUCCESS)
      return rc;
    
    for (i = 0; i < side_numVert; i++) {
      for (j = 0; j < dim; j++)
        side_vertCoords[i][j] = side_mesh_coords[side_vertexIDs[i]*dim+j];
      for (j = dim; j < 3; j++)
        side_vertCoords[i][j] = 0.;
    }
    
    // FC_Coords always have 3 dimensions, pad with zeros
    for (i = 0; i < dim; i++) 
      attach_orig[i] = node_mesh_coords[nodesetID*dim+i];
    for (i = dim; i < 3; i++)
      attach_orig[i] = 0.;
    
    _fc_calcQuadParams(side_vertCoords, attach_orig, &numParam, &params);
    if (numParam != 2 )
      return FC_ERROR;
    // FIX? remove when we think things are working, or add reasonable err?
    if (params[0] < 0. || params[0] > 1. || params[1] < 0. || params[1] > 1.) {
      free(params);
      return FC_ERROR;
    }
  }
  
  // create room for return values
  ufails = malloc(sizeof(double)*numStep);
  if (ufails == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    free(params);
    return FC_MEMORY_ERROR;
  }
  for (i = 0; i < numStep; i++)
    ufails[i] = -1.;
  
  // For each step
  for (i = 0; i < numStep; i++) {
    int vertexID;
    double* node_displ_coords;  // displaced vertex coords for node mesh
    double* side_displ_coords;  // displaced vertex coords for side mesh
    FC_Coords side_vertCoords[4]; // array of coords for side vertices (max 4)
    FC_Coords node_point, side_point;
    double dist[3], norm[3];
    
    // get node coordinates
    rc = fc_getDisplacedMeshCoords(node_mesh, node_displs[i], 
                                   &node_displ_coords);
    if (rc != FC_SUCCESS) {
      free(params);
      free(ufails);
      return rc;
    }
    for (j = 0; j < 3; j++)
      node_point[j] = node_displ_coords[dim*nodesetID+j];
    free(node_displ_coords);
    
    // get attachment on side set coordinates & normal of sideset
    rc = fc_getDisplacedMeshCoords(side_mesh, side_displs[i], 
                                   &side_displ_coords);
    if (rc != FC_SUCCESS) {
      free(params);
      free(ufails);
      return rc;
    }
    for (j = 0; j < side_numVert; j++) {
      vertexID = side_vertexIDs[j];
      for (k = 0; k < dim; k++) 
        side_vertCoords[j][k] = side_displ_coords[vertexID*dim + k];
      for (k = dim; k < 3; k++)
        side_vertCoords[j][k] = 0.;
    }
    _fc_calcQuadLocation(side_vertCoords, numParam, params, &side_point);
    _fc_calcSurfaceNormal(side_numVert, side_vertCoords, &norm);
    free(side_displ_coords);
    
    // create a vector from side point to node point
    for (j = 0; j < 3; j++)
      dist[j] = node_point[j] - side_point[j];
    
    // decompose distance vector along sideset's normal
    _fc_decomposeVector(3, dist, dim, norm, &unorm, &utang);
    
    // calculate ufail
    // fabs needed - unorm and untang might be negative and exponent odd
    ufails[i] = pow(fabs(unorm/umaxnorm), exponent) + 
      pow(fabs(utang/umaxtang), exponent);
  }
  free(params);
  
  // assign return values
  *ufails_p = ufails;
  
  // decide if and when the spot weld breaks 
  // (this is when ufail first goes above 1)
  for (i = 0; i < numStep; i++) {
    if (ufails[i] > 1.0 || FC_VALUE_EQUIV(ufails[i], 1.0, 10*DBL_EPSILON,
					  DBL_MIN)) {
      *broken_id = i; 
      break;
    }
  }
  
  return FC_SUCCESS;
}



