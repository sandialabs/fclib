/**
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
 * \file exodump.c
 * \brief Exodus file explorer
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/tools/exodump.c,v $
 * $Revision: 1.29 $
 * $Date: 2007/08/02 06:58:44 $
 *
 * \description
 *    Program to write metadats from an Exodus II file out to stdout. This is
 *    to help developers to see underlying Exodus data model. It is not a
 *    comprehesive writeout.
 *
 *    Usage: exodump [-options] filename
 *
 *    In the writeout numbers that are exodus ids are indicated with (ID = numberhere),
 *    as opposed to non-exodus numbers that are written out as a convenience  -- e.g., 
 *    numbering lists of records, and lists of connectivities for an entity where the id of the 
 *    entity is not returned as part of the call (e.g., block conns call return an array of
 *    all conns/enttiy in sequential order but without explict return for the entity id
 *    (just explict conn id)) -- are indicated with #numberhere. In the latter case, these mostly
 *    start with 0, as they are iterations through a loop, however in the case of node,face,edge,elems
 *    in a list, they start with 1 since the exodus numbering for those will start with 1.
 *    
 *     
 *
 * \todo
 *   - documentation says there is no edge order map.
 *     is there a face order map?
 *   - could rewrite this to use the new interface entirely
 *   - is there anyway to know if there is a set extra list or not
 *   - waiting for elem set prop to be fixed
 *
 * \modifications
 *    - 08/02/2005 WSD Created.
 *    - 04/16/2007 ACG adding nodeset and sideset variable and 
 *                      truth table writeout
 *    - 04/18/2007 ACG commencing adding optional 4.46 items
 *    - 07/17-18/2007 ACG elem order map out and back in. error behavior 
 *                      may have changed at some point along the line
 *                      but it is consistent with the exodus test reader
 *    - 07/17/2007 ACG adding nodal attributes (but no set ones). 
 *    - 10/04/2007 ACG adding element subsets
 *    - 10/04/2007 ACG exodus 4.58 API calls are no longer conditional
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "netcdf.h"
#include "exodusII.h"
#include "exodusII_ext.h"


enum{NODESET=0, ELEMSET, FACESET, EDGESET, SIDESET};
enum{ELEMBLK=0, FACEBLK, EDGEBLK};
enum{NODE, ELEM, FACE, EDGE}; //generic types (for maps)
#define NUMSETTYPES 5
#define NUMBLKTYPES 3 
#define NUMTYPES 4

int main(int argc, char** argv)
{
  char* filename;
  int exoid, error;
  int num_dim, num_time_step;
  int num_node, num_elem, num_face, num_edge;
  int num_blk[NUMBLKTYPES];
  int num_set[NUMSETTYPES];
  int num_blk_prop[NUMBLKTYPES];
  int num_set_prop[NUMSETTYPES];
  int num_glob_var, num_nodal_var;
  int num_blk_var[NUMBLKTYPES];
  int num_set_var[NUMSETTYPES];
  int num_map[4]; 
  int i, j, k, ii, kk;
  int CPU_word_size, IO_word_size;
  int* ids;
  float version;
  char title[MAX_LINE_LENGTH+1], *coord_names[3], **prop_names;
  char namebuf[MAX_STR_LENGTH+1];
  char *qa_records[100][4], **infos, **var_names;
  int num_qa_rec, num_info;

  float fdum;
  char* cdum;
  int* map;
  int* truth_tab;
  int getconns = 0;

  ex_init_params ex_params;

  (void) kk; //for compilation of pre-446

  // handle arguments
  if (argc < 2 || argc > 3 || (argc == 3 && (strcmp(argv[1],"-a")))) {
    printf("Usage: %s [-options] <exodus file>\n", argv[0]);
    printf("   -a         : print all block connectivities, EB attributes,\n");
    printf("                set entity lists, orientation arrays\n");
    fflush(NULL);
    exit(-1);
  }
  getconns = 0;
  if (argc == 2){ 
    filename = argv[1];
  } else {
    filename = argv[2];
    getconns = 1;    
  }

  //init
  num_dim = 0;
  num_node = 0;
  num_elem = 0;
  num_face = 0;
  num_edge = 0;
  for (i = 0; i < NUMTYPES; i++){
    num_map[i] = 0;
  }
  for (i = 0; i < NUMBLKTYPES; i++){
    num_blk[i] = 0;
    num_blk_prop[i] = 0;
    num_blk_var[i] = 0;
  }
  for (i = 0; i < NUMSETTYPES; i++){
    num_set[i] = 0;
    num_set_prop[i] = 0;
    num_set_var[i] = 0;
  }
  num_nodal_var = 0;
  num_glob_var = 0;

  // ** init exodus

  ex_opts(EX_VERBOSE);

  // ** open Exodus II file

  // open file
  //  CPU_word_size = 0;
  CPU_word_size = 8;
  IO_word_size = 0;
  exoid = ex_open(filename, EX_READ, &CPU_word_size, &IO_word_size, &version);
  if (exoid < 0) {
    fprintf(stderr, "Error: failed to open file '%s'\n", filename);
    fflush(NULL);
    exit(-1);
  }
  //print results
  printf("Database: '%s'\n", filename);
  printf("  CPU_word_size = %d, IO_word_size = %d, version = %f\n",
         CPU_word_size, IO_word_size, version);
  fflush(NULL);

  // ** Read database parameters

  // Get info  
  error = ex_get_init_ext(exoid,&ex_params);
  if (error != 0) {
    fprintf(stderr, "Error: failed to read database parameters\n");
    fflush(NULL);
    exit(-1);
  }

  strncpy(title,ex_params.title,MAX_LINE_LENGTH+1);
  num_dim = ex_params.num_dim;
  num_node = ex_params.num_nodes;
  num_elem = ex_params.num_elem;
  num_face = ex_params.num_face;
  num_edge = ex_params.num_edge;
  num_map[NODE] = ex_params.num_node_maps;
  num_map[ELEM] = ex_params.num_elem_maps;
  num_map[FACE] = ex_params.num_face_maps;
  num_map[EDGE] = ex_params.num_edge_maps;
  num_blk[ELEMBLK] = ex_params.num_elem_blk;
  num_blk[FACEBLK] = ex_params.num_face_blk;
  num_blk[EDGEBLK] = ex_params.num_edge_blk;
  num_set[NODESET] = ex_params.num_node_sets;
  num_set[ELEMSET] = ex_params.num_elem_sets;
  num_set[FACESET] = ex_params.num_face_sets;
  num_set[EDGESET] = ex_params.num_edge_sets;
  num_set[SIDESET] = ex_params.num_side_sets;

  error = ex_inquire(exoid, EX_INQ_QA, &num_qa_rec, &fdum, cdum);
  if (error != 0) {
    fprintf(stderr, "Error: failed to read # of QA records\n");
    fflush(NULL);
    exit(-1);
  }
  error = ex_inquire(exoid, EX_INQ_INFO, &num_info, &fdum, cdum);
  if (error != 0) {
    fprintf(stderr, "Error: failed to read # of info records\n");
    fflush(NULL);
    exit(-1);
  }

  for (i = 0; i < NUMBLKTYPES; i++){
    int propinq = 0;
    char* paraminq;

    if (i != ELEMBLK) continue;

    switch (i){
    case ELEMBLK:
      propinq = EX_INQ_EB_PROP;
      paraminq = "e";
      break;
    case FACEBLK:
      propinq = EX_INQ_FACE_PROP;
      paraminq = "f";
      break;
    case EDGEBLK:
      propinq = EX_INQ_EDGE_PROP;
      paraminq = "l";
      break;
    default:
      fprintf(stderr, "Error: bad block type\n");
      fflush(NULL);
      exit(-1);
    }

    error = ex_inquire(exoid, propinq, &num_blk_prop[i], &fdum, cdum);
    if (error != 0) {
      fprintf(stderr, "Error: failed to read # of block properties\n");
      fflush(NULL);
      exit(-1);
    }
    error = ex_get_var_param(exoid, paraminq, &num_blk_var[i]);
    if (error != 0) {
      fprintf(stderr, "Error: failed to get # of blk vars\n");
      fflush(NULL);
      exit(-1);
    }
  }


  for (i = 0; i < NUMSETTYPES; i++){
    int propinq = 0;
    char *paraminq;
    switch (i){
    case NODESET:
      propinq = EX_INQ_NS_PROP;
      paraminq = "m";
      break;
    case ELEMSET:
      propinq = EX_INQ_ELS_PROP;
      paraminq = "t";
      break;
    case FACESET:
      propinq = EX_INQ_FS_PROP;
      paraminq = "a";
      break;
    case EDGESET:
      propinq = EX_INQ_ES_PROP;
      paraminq = "d";
      break;
    case SIDESET:
      propinq = EX_INQ_SS_PROP;
      paraminq = "s";
      break;
    default:
      fprintf(stderr, "Error: bad set type\n");
      fflush(NULL);
      exit(-1);
    }

    error = ex_inquire(exoid, propinq, &num_set_prop[i], &fdum, cdum);
    if (error != 0) {
      fprintf(stderr, "Error: failed to read # of set properties\n");
      fflush(NULL);
      exit(-1);
    }
    error = ex_get_var_param(exoid, paraminq, &num_set_var[i]);
    if (error != 0) {
      fprintf(stderr, "Error: failed to get # of blk vars\n");
      fflush(NULL);
      exit(-1);
    }
  }
      

  error = ex_get_var_param(exoid, "g", &num_glob_var);
  if (error != 0) {
    fprintf(stderr, "Error: failed to get # of global vars\n");
    fflush(NULL);
    exit(-1);
  }
  error = ex_get_var_param(exoid, "n", &num_nodal_var);
  if (error != 0) {
    fprintf(stderr, "Error: failed to get # of nodal vars\n");
    fflush(NULL);
    exit(-1);
  }
  error = ex_inquire(exoid, EX_INQ_TIME, &num_time_step, &fdum, cdum);
  if (error != 0) {
    fprintf(stderr, "Error: failed to get # of time steps\n");
    fflush(NULL);
    exit(-1);
  }

  // Print it
  printf("  Title = '%s'\n", title);
  printf("  num_dim = %d\n", num_dim);
  printf("  num_node = %d, num_elem = %d, num_face = %d, num_edge = %d\n",
         num_node, num_elem, num_face, num_edge);
  printf("  num_node_map = %d, num_elem_map = %d, num_face_map = %d, num_edge_map = %d\n",
         num_map[NODE], num_map[ELEM], num_map[FACE], num_map[EDGE]);
  printf("  num_elem_blk = %d, num_face_blk = %d, num_edge_blk = %d\n", 
	 num_blk[ELEMBLK], num_blk[FACEBLK], num_blk[EDGEBLK]);
  printf("  num_elem_blk_prop = %d, num_face_blk_prop = %d, num_edge_blk_prop = %d\n", 
	 num_blk_prop[ELEMBLK], num_blk_prop[FACEBLK], num_blk_prop[EDGEBLK]);
  printf("  num_node_set = %d, num_elem_set = %d, num_face_set = %d, num_edge_set = %d, num_side_set = %d\n", 
	 num_set[NODESET], num_set[ELEMSET], num_set[FACESET], num_set[EDGESET], num_set[SIDESET]);
  printf("  num_node_set_prop = %d, num_elem_set_prop = %d, num_face_set_prop = %d, num_edge_set_prop = %d, num_side_set_prop = %d\n", 
	 num_set_prop[NODESET], num_set_prop[ELEMSET], num_set_prop[FACESET], num_set_prop[EDGESET], num_set_prop[SIDESET]);
  printf("  num_glob_var = %d, num_nodal_var = %d\n",
         num_glob_var, num_nodal_var)
;  printf("  num_elem_blk_var = %d, num_face_blk_var = %d, num_edge_blk_var = %d\n", 
	 num_blk_var[ELEMBLK], num_blk_var[FACEBLK], num_blk_var[EDGEBLK]);
  printf("  num_node_set_var = %d, num_elem_set_var = %d, num_face_set_var = %d, num_edge_set_var = %d, num_side_set_var = %d\n", 
	 num_set_var[NODESET], num_set_var[ELEMSET], num_set_var[FACESET], num_set_var[EDGESET], num_set_var[SIDESET]);
  printf("  num_time_step = %d\n", num_time_step);
  printf("  num_qa_rec = %d, num_info = %d\n", num_qa_rec, num_info);

  printf("\n");
  fflush(NULL);

  // ** Look at QA records

  if (num_qa_rec > 0) {
    for (i = 0; i < num_qa_rec; i++) 
      for (j = 0; j < 4; j++)
        qa_records[i][j] = (char*)malloc((MAX_LINE_LENGTH+1)*sizeof(char));
    error = ex_get_qa(exoid, qa_records);
    if (error != 0) {
      fprintf(stderr, "Error: failed to read qa records\n");
      fflush(NULL);
      exit(-1);
    }
    for (i = 0; i < num_qa_rec; i++) {
      printf("QA Record #%d:\n", i);
      for (j = 0; j < 4; j++)
        printf("  %s\n", qa_records[i][j]);
      printf("\n");
      fflush(NULL);
    }
    for (i = 0; i < num_qa_rec; i++) 
      for (j = 0; j < 4; j++)
        free(qa_records[i][j]);
  }

  // ** Look at info records

  if (num_info > 0) {
    infos = (char**)malloc(num_info*sizeof(char*));
    for (i = 0; i < num_info; i++)
      infos[i] = (char*)malloc((MAX_LINE_LENGTH+1)*sizeof(char));
    error = ex_get_info(exoid, infos);
    if (error != 0) {
      fprintf(stderr, "Error: failed to read info records\n");
    fflush(NULL);
    exit(-1);
    }
    for (i = 0; i < num_info; i++) {
      printf("Info Record #%d:\n", i);
      printf(" %s\n", infos[i]);
      printf("\n");
      fflush(NULL);
    }
    for (i = 0; i < num_info; i++)
      free(infos[i]);
    free(infos);
  }

  // ** Look at the coords

  for (i = 0; i < num_dim; i++)
    coord_names[i] = (char*)malloc((MAX_STR_LENGTH+1)*sizeof(char));
  error = ex_get_coord_names(exoid, coord_names);
  if (error != 0) {
    fprintf(stderr, "Error: failed to read coord names\n");
    fflush(NULL);
    exit(-1);
  }
  printf("Coordinate names:\n");
  for (i = 0; i < num_dim; i++)
    printf("  %s\n", coord_names[i]);
  printf("\n");
  fflush(NULL);
  for (i = 0; i < num_dim; i++)
    free(coord_names[i]);

  //node, elem, face, and edge maps
  for (i = 0; i < NUMTYPES; i++){
    int obj_type;
    char* str;
    switch (i){
    case NODE:
      obj_type = EX_NODE_MAP;
      str = "Node Number";
      break;
    case ELEM:
      obj_type = EX_ELEM_MAP;
      str = "Elem Number";
      break;
    case FACE:
      obj_type = EX_FACE_MAP;
      str = "Face Number";
      break;
    case EDGE:
      obj_type = EX_EDGE_MAP;
      str = "Edge Number";
      break;
    default:
      fprintf(stderr, "Error: bad type for map\n");
      fflush(NULL);
      exit(-1);
    }

    if (num_map[i] != 0){
      ids = (int*)malloc(num_map[i]*sizeof(int));
      error = ex_get_ids(exoid, obj_type, ids);
      if (error != 0) {
	fprintf(stderr, "Error: failed to read map ids\n");
	fflush(NULL);
	exit(-1);
      }
      
      for (j = 0; j < num_map[i]; j++){
	error = ex_get_num_map(exoid, obj_type, ids[j], map);
	free(map);
	if (error < 0) {
	  fprintf(stderr, "Error: failed to read %s map %d \n",str,j);
	  fflush(NULL);
	  exit(-1);
	}
	else if (error > 0)
	  printf("%s map (ID = %d) does not exist\n",str,ids[j]);
	else 
	  printf("%s map (ID = %d) exists\n",str,ids[j]);
	printf("\n");
	fflush(NULL);
      }
    }

    //only look at the order map if its elements - can have this w/o a number map
    if (i == ELEM){
      map = (int*)malloc(num_elem*sizeof(int));
      error = ex_get_map(exoid, map);
      //      printf ("\nafter ex_get_map, error = %3d\n", error);
      if (error < 0) {
	fprintf(stderr, "Error: failed to read elem order map\n");
	fflush(NULL);
	exit(-1);
      } else if (error > 0){
	printf("Warning reading elem order map\n");
      } else {
	printf("Elem order map exists\n");
	//FCLib does nothing with the order map
	//    for (i=0; i< num_elem; i++){
	//      printf ("elem_map(%d) = %d \n", i, map[i]);
	//    }
	printf("\n");
      }
      fflush(NULL);
      free(map);
    }
  } // numtypes


  for (i = 0; i < NUMBLKTYPES; i++){
    int obj_type;
    char* str;
    switch (i){
    case ELEMBLK:
      obj_type = EX_ELEM_BLOCK;
      str = "ElemBlk";
      break;
    case FACEBLK:
      obj_type = EX_FACE_BLOCK;
      str = "FaceBlk";
      break;
    case EDGEBLK:
      obj_type = EX_EDGE_BLOCK;
      str = "EdgeBlk";
      break;
    default:
      fprintf(stderr, "Error: bad block type\n");
      fflush(NULL);
      exit(-1);
    }

    if (num_blk[i] == 0){
      continue;
    }

    ids = (int*)malloc(num_blk[i]*sizeof(int));
    error = ex_get_ids(exoid, obj_type, ids);
    if (error != 0) {
      fprintf(stderr, "Error: failed to read block ids\n");
      fflush(NULL);
      exit(-1);
    }
    for (j = 0; j < num_blk[i]; j++) {
      int num_entries_in_blk, num_nodes_per_entry, num_attr_per_entry;
      int num_edges_per_entry, num_faces_per_entry;
      char blk_type_name[MAX_STR_LENGTH+1];

      error = ex_get_block(exoid, obj_type, ids[j], blk_type_name, &num_entries_in_blk,
			   &num_nodes_per_entry, &num_edges_per_entry,
			   &num_faces_per_entry, &num_attr_per_entry);
      if (error != 0) {
	fprintf(stderr, "Error: failed to read block\n");
	fflush(NULL);
	exit(-1);
      }
      namebuf[0] = '\0';

      error = ex_get_name(exoid, obj_type, ids[j], namebuf);
      if (error != 0) {
	fprintf(stderr, "Error: failed to get block name\n");
	fflush(NULL);
	exit(-1);
      }

      printf("%s (ID = %d): '%s'\n", str, ids[j], namebuf);
      switch(i){
      case ELEMBLK:
	printf("  num_elem_in_blk = %d, num_nodes_per_elem = %d, num_faces_per_elem = % d, num_edges_per_elem = % d, elem_type = %s\n",
	       num_entries_in_blk, num_nodes_per_entry, num_faces_per_entry, num_edges_per_entry, blk_type_name);
	break;
      case FACEBLK:
	printf("  num_faces_in_blk = %d, num_nodes_per_face = %d, face_type = %s\n",
	       num_entries_in_blk, num_nodes_per_entry, blk_type_name);
	break;
      case EDGEBLK:
	printf("  num_edges_in_blk = %d, num_nodes_per_edge = %d, edge_type = %s\n",
	       num_entries_in_blk, num_nodes_per_entry, blk_type_name);
	break;
      default:
	fprintf(stderr, "Error: bad block type\n");
	fflush(NULL);
	exit(-1);
      }

      printf("  Attributes: %d:\n", num_attr_per_entry);
      fflush(NULL);

      if (num_attr_per_entry > 0) {
	char** attr_names = (char**)malloc(num_attr_per_entry*sizeof(char*));
	double* attr_vals = (double*)malloc(num_entries_in_blk*sizeof(double));
	if (!attr_names || !attr_vals){
	  exit(-1);
	}
	for (k = 0; k < num_attr_per_entry; k++){
	  attr_names[k] = (char*)malloc((MAX_STR_LENGTH+1)*sizeof(char*));
	  if (!attr_names[k]){
	    exit(-1);
	  }
	}
	error = ex_get_attr_names(exoid, obj_type, ids[j], attr_names);
	if (error != 0) {
	  fprintf(stderr, "Error: failed to get blk attr names\n");
	  fflush(NULL);
	  exit(-1);
	}
	for (k = 0; k < num_attr_per_entry; k++) {
	  printf("   #%d: '%s'\n", k, attr_names[k]);
	  if (getconns){
	    error = ex_get_one_attr(exoid, obj_type, ids[j], k+1, attr_vals);
	    if (error != 0) {
	      fprintf(stderr, "Error: failed to get blk attr\n");
	      fflush(NULL);
	      exit(-1);
	    }
	    for (kk = 0; kk < num_entries_in_blk; kk++){
	      printf("\t %d: %g\n", kk, attr_vals[kk]);
	      fflush(NULL);
	    }
	  } //getconns
	} //numattr
	for (k = 0 ; k < num_attr_per_entry; k++){
	  free(attr_names[k]);
	}
	free(attr_names);
	free(attr_vals);
      }
      printf("\n");
      fflush(NULL);

      // ** block conns (internal ids, so the full sequential numbering is represented)
      if (getconns){      
	int *node_connect;
	int *face_connect;
        int *edge_connect;
	node_connect = (int*) malloc(num_entries_in_blk*num_nodes_per_entry*sizeof(int));
	face_connect = NULL;
	edge_connect = NULL;
	if (i == ELEMBLK){
	  face_connect = (int*) malloc(num_entries_in_blk*num_faces_per_entry*sizeof(int));
	  edge_connect = (int*) malloc(num_entries_in_blk*num_edges_per_entry*sizeof(int));	
	}
	error = ex_get_conn(exoid, obj_type, ids[j], node_connect, edge_connect, face_connect);
	if (error != 0) {
	  fprintf(stderr, "Error: failed to get block conns\n");
	  fflush(NULL);
	  exit(-1);
	}
	printf(" node conns:\n");
	for (k = 0; k < num_entries_in_blk; k++){
	  printf("  Entry #%d (", k+1); //exodus numbering starts with 1 
	  for (kk = 0; kk < num_nodes_per_entry; kk++){
	    printf("%d ", node_connect[k*num_nodes_per_entry+kk]);
	  }
	  printf(")\n");
	}
	fflush(NULL);
	if ( i == ELEMBLK ){
	  if (num_faces_per_entry > 0){
	    printf(" face conns:\n");
	    for (k = 0; k < num_entries_in_blk; k++){ //elems
	      printf("  Entry #%d (", k+1); //exodus numbering starts with 1
	      for (kk = 0; kk < num_faces_per_entry; kk++){
		printf("%d ", face_connect[k*num_faces_per_entry+kk]);
	      }
	      printf(")\n");
	    }
	  } else {
	    printf(" no face conns\n");
	  } 
	  fflush(NULL);
	  if (num_edges_per_entry > 0){
	    printf(" edge conns:\n");
	    for (k = 0; k < num_entries_in_blk; k++){ //elems
	      printf("  Entry #%d (", k+1); //exodus numbering starts with 1
	      for (kk = 0; kk < num_edges_per_entry; kk++){
		printf("%d ", edge_connect[k*num_edges_per_entry+kk]);
	      }
	      printf(")\n");
	    }
	  } else {
	    printf(" no edge conns\n");
	  } 
	  fflush(NULL);
	}
	if (edge_connect) free(edge_connect);
	if (face_connect) free(face_connect);
	if (node_connect) free(node_connect);
      }
      printf("\n");
    }

    // ** block properties

    if (num_blk_prop[i] > 0) {
      prop_names = (char**)malloc(num_blk_prop[i]*sizeof(char*));
      for (j = 0; j < num_blk_prop[i]; j++)
	prop_names[j] = (char*)malloc((MAX_STR_LENGTH+1)*sizeof(char));
      error = ex_get_prop_names(exoid, obj_type, prop_names);
      if (error != 0) {
	fprintf(stderr, "Error: failed to get blk prop names\n");
	fflush(NULL);
	exit(-1);
      }
      for (j = 0; j < num_blk_prop[i]; j++) {
	printf("%s Property  #%d: '%s'\n", str, j, prop_names[j]);
	for (k = 0; k < num_blk[i]; k++) {
	  int prop_value;
	  error = ex_get_prop(exoid, obj_type, ids[k], prop_names[j],
			      &prop_value);
	  if (error != 0) {
	    fprintf(stderr, "Error: failed to get blk prop\n");
	    exit(-1);
	  }
	  printf("  %s (ID = %d): %d\n", str, ids[k], prop_value);
	}
	printf("\n");
	fflush(NULL);
      }

      for (j = 0; j < num_blk_prop[i]; j++)
	free(prop_names[j]);
      free(prop_names);
    }

    // done with block ids
    free(ids);
  }


  for (i = 0; i < NUMSETTYPES; i++){
    int obj_type;
    char* str;
    switch (i){
    case NODESET:
      str = "NodeSet";
      obj_type = EX_NODE_SET;
      break;
    case ELEMSET:
      str = "ElemSet";
      obj_type = EX_ELEM_SET;
      break;
    case FACESET:
      str = "FaceSet";
      obj_type = EX_FACE_SET;
      break;
    case EDGESET:
      str = "EdgeSet";
      obj_type = EX_EDGE_SET;
      break;
    case SIDESET:
      str = "SideSet";
      obj_type = EX_SIDE_SET;
      break;
    default:
      fprintf(stderr, "Error: bad set type\n");
      fflush(NULL);
      exit(-1);
    }

    if (num_set[i] == 0){
      continue;
    }

    ids = (int*)malloc(num_set[i]*sizeof(int));
    error = ex_get_ids(exoid, obj_type, ids);
    if (error != 0) {
      fprintf(stderr, "Error: failed to read set ids\n");
      fflush(NULL);
      exit(-1);
    }
    for (j = 0; j < num_set[i]; j++) {
      int num_entries_in_set, num_df_in_set;
      int *set_entity_list, *set_extra_list;
      error = ex_get_set_param(exoid, obj_type, ids[j], &num_entries_in_set,
                                    &num_df_in_set);
      if (error != 0) {
        fprintf(stderr, "Error: failed to get set params\n");
        exit(-1);
      }
      namebuf[0] = '\0';
      error = ex_get_name(exoid, obj_type, ids[j], namebuf);
      if (error != 0) {
        fprintf(stderr, "Error: failed to get set name\n");
        fflush(NULL);
        exit(-1);
      }
      printf("%s (ID = %d): '%s'\n", str, ids[j], namebuf);
      printf("  num_entries_in_set = %d, num_df_in_set = %d\n",
             num_entries_in_set, num_df_in_set);
      printf("\n");
      fflush(NULL);

      if (getconns){
	set_entity_list = (int*)malloc(num_entries_in_set*sizeof(int));
	if (i != NODESET && i!= ELEMSET){
	  //FIXME - is there anyway to know if there is an extra list or not?
	  set_extra_list = (int*)malloc(num_entries_in_set*sizeof(int));
	} else {
	  set_extra_list = NULL;
	}
	error = ex_get_set(exoid, obj_type, ids[j], set_entity_list, set_extra_list);
	if (error < 0 ) { 
	  fprintf(stderr, "Error: failed to get set entities\n");
	  fflush(NULL);
	  exit(-1);
	} else if (error > 0){
	  fprintf(stderr, "Warning: failed to get set entities\n");
	  fflush(NULL);
	} else {
	  switch (i){
	  case NODESET:
	    printf("  node_id:\n");
	    for (k = 0; k < num_entries_in_set; k++){
	      printf("    %d\n", set_entity_list[k]);
	    }
	    break;
	  case ELEMSET:
	    printf("  elem_id:\n");
	    for (k = 0; k < num_entries_in_set; k++){
	      printf("    %d\n", set_entity_list[k]);
	    }
	    break;
	  case SIDESET:
	    printf("  elem_id     side_index:\n");
	    for (k = 0; k < num_entries_in_set; k++){
	      printf("    %d            %d\n", set_entity_list[k],
		     set_extra_list[k]);
	    }
	    break;
	  default:
	    printf("  entity_id   orientation:\n");
	    for (k = 0; k < num_entries_in_set; k++){
	      printf("    %d            %d\n", set_entity_list[k],
		     set_extra_list[k]);
	    }
	  }
	}
	if (set_entity_list) free(set_entity_list);
	if (set_extra_list) free(set_extra_list);
      }
      printf("\n");
      fflush(NULL);
    }

    
    // ** set properties

    if (num_set_prop[i] > 0){
      //FIXME: commented out for elem_set until that gets fixed.
      if (i == ELEMSET) {
      	printf("Not printing out elemset prop at this time\n");
      	break;
      }

      prop_names = (char**)malloc(num_set_prop[i]*sizeof(char*));
      for (j = 0; j < num_set_prop[i]; j++)
        prop_names[j] = (char*)malloc((MAX_STR_LENGTH+1)*sizeof(char));
      error = ex_get_prop_names(exoid, obj_type, prop_names);
      if (error != 0) {
        fprintf(stderr, "Error: failed to get set prop names %d\n",i);
        fflush(NULL);
	exit(-1);
      } 
      for (j = 0; j < num_set_prop[i]; j++) {
	printf("%s Property #%d: '%s'\n", str, j, prop_names[j]);
	for (k = 0; k < num_set[i]; k++) {
	  int prop_value;
	  //this will return with EX_WARN if the set is empty.
	  //that is defined in exodusII_int.h
	  error = ex_get_prop(exoid, obj_type, ids[k], prop_names[j],
			      &prop_value);
	  if (error < 0) {
	    fprintf(stderr, "Error: failed to get set props\n");
	    fflush(NULL);
	    exit(-1);
	  } else if (error > 0 ){
	    fprintf(stderr, "Warning: failed to get set props");
	  } else {
	    printf("  %s (ID = %d): %d\n", str, ids[k], prop_value);
	  }
	}
	printf("\n");
	fflush(NULL);
      }
      for (j = 0; j < num_set_prop[i]; j++)
	free(prop_names[j]);
      free(prop_names);
    }

    // done with set ids
    free(ids);
  }


  //get nodal attr 
  {
    int num_attrs = 0;
    error = ex_get_attr_param(exoid, EX_NODAL, 0, &num_attrs);
    printf ("NODAL Attributes: %d\n", num_attrs);
    if (num_attrs > 0) {
      char** attr_names = (char**)malloc(num_attrs*sizeof(char*));
      double* attr_vals = (double*)malloc(num_node*sizeof(double));
      for (j=0; j< num_attrs; j++) {
	attr_names[j] = (char *)calloc ((MAX_STR_LENGTH+1), sizeof(char));
      }
      error = ex_get_attr_names (exoid, EX_NODAL, 0, attr_names);
      if (error != 0) {
	fprintf(stderr, "Error: failed to get nodal attr names\n");
	fflush(NULL);
	exit(-1);
      }
      for (j = 0; j < num_attrs; j++) {
	printf("   #%d: '%s'\n", j, attr_names[j]);
	if (getconns){
	  error = ex_get_one_attr(exoid, EX_NODAL, 0, j+1, attr_vals);
	  if (error != 0) {
	    fprintf(stderr, "Error: failed to get blk attr\n");
	    fflush(NULL);
	    exit(-1);
	  }
	  for (k = 0; k < num_node; k++){
	    printf("\t %d: %g\n", k, attr_vals[k]);
	    fflush(NULL);
	  }
	} //getconns
      } //numattr
      for (k = 0 ; k < num_attrs; k++){
	free(attr_names[k]);
      }
      free(attr_names);
      free(attr_vals);
    }
    printf("\n");
    fflush(NULL);
  }

  // ** global variable names
  
  if (num_glob_var > 0) {
    var_names = (char**)malloc(num_glob_var*sizeof(char*));
    for (i = 0; i < num_glob_var; i++)
      var_names[i] = (char*)malloc((MAX_STR_LENGTH+1)*sizeof(char));
    error = ex_get_var_names(exoid, "g", num_glob_var, var_names);
    if (error != 0) {
      fprintf(stderr, "Error: Failed to get global var names\n");
      fflush(NULL);
      exit(-1);
    }
    printf("Global Variables:\n");
    for (i = 0; i < num_glob_var; i++)
      printf(" #%d: %s\n", i, var_names[i]);
    for (i = 0; i < num_glob_var; i++)
      free(var_names[i]);
    free(var_names);
    printf("\n");
    fflush(NULL);
  }

  // ** nodal variable names
  
  if (num_nodal_var > 0) {
    var_names = (char**)malloc(num_nodal_var*sizeof(char*));
    for (i = 0; i < num_nodal_var; i++)
      var_names[i] = (char*)malloc((MAX_STR_LENGTH+1)*sizeof(char));
    error = ex_get_var_names(exoid, "n", num_nodal_var, var_names);
    if (error != 0) {
      fprintf(stderr, "Error: Failed to get nodal var names\n");
      fflush(NULL);
      exit(-1);
    }
    printf("Nodal Variables:\n");
    for (i = 0; i < num_nodal_var; i++)
      printf(" #%d: %s\n", i, var_names[i]);
    for (i = 0; i < num_nodal_var; i++)
      free(var_names[i]);
    free(var_names);
    printf("\n");
    fflush(NULL);
  }

  // ** blk vars
  for (ii = 0; ii < NUMBLKTYPES; ii++){
    char* varinq;
    char* str;
    switch (ii){
    case ELEMBLK:
      varinq = "e";
      str = "ElemBlk";
      break;
    case FACEBLK:
      varinq = "f";
      str = "FaceBlk";
      break;
    case EDGEBLK:
      varinq = "l";
      str = "EdgeBlk";
      break;
    default:
      fprintf(stderr, "Error: bad block type\n");
      fflush(NULL);
      exit(-1);
    }


    // blk variable names
    if (num_blk_var[ii] > 0) {
      var_names = (char**)malloc(num_blk_var[ii]*sizeof(char*));
      for (i = 0; i < num_blk_var[ii]; i++)
	var_names[i] = (char*)malloc((MAX_STR_LENGTH+1)*sizeof(char));
      error = ex_get_var_names(exoid, varinq, num_blk_var[ii], var_names);
      if (error != 0) {
	fprintf(stderr, "Error: Failed to get  blk var names\n");
	fflush(NULL);
	exit(-1);
      }
      printf("%s Variables:\n",str);
      for (i = 0; i < num_blk_var[ii]; i++)
	printf(" #%d: %s\n", i, var_names[i]);
      for (i = 0; i < num_blk_var[ii]; i++)
	free(var_names[i]);
      free(var_names);
      printf("\n");
      fflush(NULL);
    }

    // ** blk var truth table
    if (num_blk_var[ii] > 0) {
      truth_tab = (int*)malloc(num_blk[ii]*num_blk_var[ii]*sizeof(int));
      error = ex_get_var_tab(exoid, varinq, num_blk[ii], num_blk_var[ii], truth_tab);
      if (error != 0) {
	fprintf(stderr, "Error: Failed to get blk truth table\n");
	fflush(NULL);
	exit(-1);
      }

      // note cols can get off here if num blocks is very large
      printf("%s Variable Truth Table:\n",str);
      printf("  %14s"," ");
      for (i = 0; i < num_blk[ii]; i++)
	printf("%8s#%d ", str, i+1); //exodus numbering starts with 1
      printf("\n");
      for (i = 0; i < num_blk_var[ii]; i++) {
	printf("  %8sVar#%d ", str, i);
	for (j = 0; j < num_blk[ii]; j++)
	  printf(" %4s%2d%3s ", " ", truth_tab[j*num_blk_var[ii]+i], " ");
	printf("\n");
      }
      fflush(NULL);
      free(truth_tab);
    }
    printf("\n");
  }


  // **set vars
  for (ii = 0; ii < NUMSETTYPES; ii++){
    char* varinq;
    char* str;
    switch (ii){
    case NODESET:
      str = "NodeSet";
      varinq = "m";
      break;
    case ELEMSET:
      str = "ElemSet";
      varinq = "t";
      break;
    case FACESET:
      str = "FaceSet";
      varinq = "a";
      break;
    case EDGESET:
      str = "EdgeSet";
      varinq = "d";
      break;
    case SIDESET:
      str = "SideSet";
      varinq = "s";
      break;
    default:
      fprintf(stderr, "Error: bad set type\n");
      fflush(NULL);
      exit(-1);
    }
  
    //vars
    if (num_set_var[ii] > 0) {
      printf("\n");
      var_names = (char**)malloc(num_set_var[ii]*sizeof(char*));
      for (i = 0; i < num_set_var[ii]; i++)
	var_names[i] = (char*)malloc((MAX_STR_LENGTH+1)*sizeof(char));
      error = ex_get_var_names(exoid, varinq, num_set_var[ii], var_names);
      if (error != 0) {
	fprintf(stderr, "Error: Failed to get set var names\n");
	fflush(NULL);
	exit(-1);
      }
      printf("%s Variables:\n",str);
      for (i = 0; i < num_set_var[ii]; i++)
	printf("  %d: %s\n", i, var_names[i]);
      for (i = 0; i < num_set_var[ii]; i++)
	free(var_names[i]);
      free(var_names);
      printf("\n");
      fflush(NULL);
    }


    // truth table
    if (num_set_var[ii] > 0) {
      truth_tab = (int*)malloc(num_set[ii]*num_set_var[ii]*sizeof(int));
      error = ex_get_var_tab(exoid, varinq, num_set[ii], num_set_var[ii], truth_tab);
      if (error != 0) {
	fprintf(stderr, "Error: Failed to get set truth table\n");
	fflush(NULL);
	exit(-1);
      }
      printf("%s Variable Truth Table:\n",str);
      printf("  %11s"," ");
      for (i = 0; i < num_set[ii]; i++)
	printf(" %6s#%d ", str, i);
      printf("\n");
      for (i = 0; i < num_set_var[ii]; i++) {
	printf("  %6s Var # %d ", str, i);
	for (j = 0; j < num_set[ii]; j++)
	  printf(" %3s%d%3s ", " ", truth_tab[j*num_set_var[ii]+i], " ");
	printf("\n");
      }
      fflush(NULL);
      free(truth_tab);
    }
  }

  exit(0);
}
  


 


