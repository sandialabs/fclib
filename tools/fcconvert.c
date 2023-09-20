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
 * \file fcconvert.c
 * \brief Convert a dataset to a different file type.
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/tools/fcconvert.c,v $
 * $Revision: 1.7 $ 
 * $Date: 2006/10/13 18:13:20 $
 *
 * \description
 *    Usage is: fcconvert -f new_type orig_file
 *
 * \todo The SAF writer does not fail gracefully (it should ignore
 *    things it doesn't understand
 *
 * \modifications
 *    01/30/06 Created Wendy Doyle
 *    03/25/06 Changed to use new fc_writeAsDataset instead of doing a copy.
 *    10/13/06 Let it "convert" to the same type. This converts to fclib's 
 *      conventions for that type. (E.g. for Exodus, it splits the global 
 *      verts into the meshes, perhaps with duplication).
 *    10/15/07 ACG removed saf
 */

#include <string.h>
#include <fc.h>

int main(int argc, char** argv) {
  int i;
  FC_ReturnCode rc;
  int makeName = 1;
  char* orig_dsName = NULL, *new_dsName = NULL;
  char* extension;
  char extensions[1][20] = { "ex2" };
  FC_Dataset orig_ds;
  FC_FileIOType new_fileType = FC_FT_NONE;
  FC_VerbosityLevel verbose_level = FC_QUIET; 
  
  // handle arguments
  if (argc < 4) {
  usage:
    printf("Usage: fcconvert [options] -t new_type orig_file\n");
    printf("\n");
    printf("options: \n");
    printf("   -o new_name  : name for the converted dataset\n");
    printf("   -h           : print this help message\n");
    printf("   -v           : verbose: print warning and error messages\n");
    printf("   -V           : very verbose: prints log and error messages\n");
    printf("\n");
    printf("Converts the dataset to the requested file type. Possible types are\n");
    printf("EXO (ExodusII). Unless a new name is provided, the new dataset's\n");
    printf("name will have the same root but a different extension. For safety,\n");
    printf("fcconvert tries to prevent you from writing to the same file by comparing\n");
    printf("path strings. This can be circumvented by using roundabout paths and\n");
    printf("symbolic links, so care should still be taken to not overwrite the original\n");
    printf("dataset.\n");
    exit(0);
  }
  for (i = 1; i < argc; i++) {
    if (!strcmp(argv[i], "-o")) {
      i++;
      if (i >= argc)
	goto usage;
      new_dsName = argv[i];
      makeName = 0;
    }
    else if (!strcmp(argv[i], "-t")) {
      i++;
      if (i >= argc)
	goto usage;
      if (!strcmp(argv[i], "EXO")) {
	new_fileType = FC_FT_EXODUS;
	extension = extensions[0];
      }
      else {
	printf("Error: unknown file type");
	goto usage;
      }
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
      if (i+1 > argc)
	goto usage;
      orig_dsName = argv[i];
    }
  }
  if (!orig_dsName || new_fileType == FC_FT_NONE)
    goto usage;

  // Make sure writing file type is supported
  if (!fc_isFileIOTypeWriteSupported(new_fileType)) 
    fc_exitIfErrorPrintf(FC_ERROR, "Request file IO type is not supported");

  // Make new name, if necessary
  if (makeName) {
    int len;
    char* root = fc_getBasenameWOExtension(orig_dsName, 1);
    len = strlen(root) + strlen(extension) + 2;
    new_dsName = malloc(len*sizeof(char));
    if (!new_dsName)
      fc_exitIfError(FC_MEMORY_ERROR);
    sprintf(new_dsName, "%s.%s", root, extension);
  }
  // Don't allow write to same file
  // FIX!? Test isn't completely strict since you could have different
  // paths to the same file and they would look like different names
  if (!strcmp(orig_dsName, new_dsName)) 
    fc_exitIfErrorPrintf(FC_ERROR, "Cannot write to the same file. Use -o "
                         "to specify the new dataset filename.");

  // init library & load original dataset
  rc = fc_setLibraryVerbosity(verbose_level);
  fc_exitIfError(rc);  
  rc = fc_initLibrary();
  fc_exitIfError(rc);
  rc = fc_loadDataset(orig_dsName, &orig_ds);
  fc_exitIfErrorPrintf(rc, "Loading dataset '%s'", orig_dsName);

  // write new dataset
  rc = fc_rewriteDataset(orig_ds, new_dsName, new_fileType);
  fc_exitIfError(rc);

  // cleanup
  if (makeName)
    free(new_dsName);
  rc = fc_deleteDataset(orig_ds);
  fc_exitIfError(rc);
  rc = fc_finalLibrary();
  fc_exitIfError(rc);

  exit(0);
}
