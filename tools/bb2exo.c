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
 * \file bb2exo.c
 * \brief Convert bounding box file to exodus dataset
 *
 * $Source: /home/Repositories/fcdmf/fclib/tools/bb2exo.c,v $
 * $Revision: 1.6 $ 
 * $Date: 2006/11/07 23:49:03 $
 *
 * \description
 *
 *   Reads a bounding box file (see \ref threshBoundBox.c) and creates
 *   a new Exodus dataset with a mesh for each bounding box.
 *
 *   The naming convention of the bb meshes is: name_step#, where
 *   the name and step# are the same as the attributes on the bb element.
 *   (The step# is appended to help pick bb since the exo file doesn't
 *   have time.)
 *
 * \todo Should probably make a version 'bb2fc' where can choose final
 *    file type with '-t' like w/ fcconvert.
 *
 * \modifications
 *    - 03/01/2006 WSD, created.
 *    - 06/29/2006 WSD, converted to new bb file format.
 */

#include <string.h>
#include <stdio.h>
#include "fc.h"
#include "fcP.h" // temporary until readBB stuff gets made public

int main(int argc, char** argv) {
  FC_ReturnCode rc;
  int i;
  char* bbfile_name = NULL;
  char* base_name;
  char* exofile_name;
  char mesh_name_buf[1024];
  FC_VerbosityLevel verbose_level = FC_QUIET;
  FC_Dataset dataset;
  FC_Mesh mesh;
  int numBB, *stepIDs;
  FC_Coords *lows, *highs;
  char** bbNames, **comments;

  // handle arguments
  if (argc < 2) {
  usage:
    printf("Will create Exodus dataset (with name file.ex2)\n");
    printf("usage: bb2exo [options] file.bb\n");
    printf("options: \n");
    printf("   -h            : print this help message\n");
    printf("   -v            : verbose: print warning and error messages\n");
    printf("   -V            : very verbose: prints log and error messages\n");
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
    else {
      if (i + 1 > argc)
	goto usage;
      bbfile_name = argv[i];
    }
  }
  if (!bbfile_name)
    goto usage;

  // init library 
  rc = fc_setLibraryVerbosity(verbose_level);
  fc_exitIfError(rc);
  rc = fc_initLibrary();
  fc_exitIfError(rc);

  // read the bb file
  rc = _fc_readBBFile(bbfile_name, &numBB, &bbNames, &stepIDs, &comments,
		      &lows, &highs);
  fc_exitIfErrorPrintf(rc, "Failed to read bb file '%s'", bbfile_name);
  if (numBB < 1) {
    fc_exitIfErrorPrintf(FC_ERROR, "No bounding boxes found in bb file '%s'", 
			 bbfile_name);
  }
    
  // create name of exodus file
  base_name = fc_getBasenameWOExtension(bbfile_name, 1);
  exofile_name = (char*)malloc((strlen(base_name)+5)*sizeof(char));
  if (!base_name || !exofile_name) 
    fc_exitIfError(FC_MEMORY_ERROR);
  sprintf(exofile_name, "%s.ex2", base_name);
  free(base_name);

  // create a dataset to put bb meshes into
  rc = fc_createDataset("Bounding Boxes", &dataset);
  fc_exitIfErrorPrintf(rc, "Failed to create dataset");

  // FIX? short by stepID?

  // For each line in bb, create a bbmesh
  for (i = 0; i < numBB; i++) {
    // FIX unhardcode size of name
    if (stepIDs[i] < 0)
      sprintf(mesh_name_buf, "%s_nostep", bbNames[i]);
    else
      sprintf(mesh_name_buf, "%s_step%d", bbNames[i], stepIDs[i]);
    rc = fc_createBoundingBoxMesh(dataset, mesh_name_buf, 3, lows[i], highs[i],
				  &mesh);
    fc_exitIfErrorPrintf(rc, "Failed to create bounding box mesh with "
			 "([%g,%g,%g],[%g,%g,%g])", lows[i][0], lows[i][1],
			 lows[i][2], highs[i][0], highs[i][1], highs[i][2]);
    free(bbNames[i]);
    free(comments[i]);
  }
  free(bbNames);
  free(comments);
  free(stepIDs);
  free(lows);
  free(highs);

  // Write the dataset
  rc = fc_writeDataset(dataset, exofile_name, FC_FT_EXODUS);
  fc_exitIfErrorPrintf(rc, "Failed to write exodus dataset");

  // all done
  fc_finalLibrary();

  exit(0);
}
  
