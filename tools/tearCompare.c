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
 * \file tearCompare.c
 * \brief Report how tears from two different datasets compare.
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/tools/tearCompare.c,v $
 * $Revision: 1.7 $ 
 * $Date: 2006/08/30 19:20:04 $
 *
 * \description
 *
 *    Usage: tearCompare tears1.aux tears2.aux
 *
 *    The underlying models of the two datasets must be the same (i.e. have
 *    the same mesh topology).
 *
 *    The report assumes that you are interested in how the 2nd set of tears
 *    matches up to the 1st. If that is not quite what you want, switch of the
 *    order of the arguments.
 *
 * \todo Produce a graph of tear relationships.
 *
 * \modifications
 *   - 07/03/2006 WSD, created.
 */

#include <string.h>
#include <stdio.h>
#include "fc.h"
#include "fcP.h" // temporary until writeBB stuff gets made public

// assumes ids are sorted
// FIX sort ids before sending in
// calculates number of intersecting, not just if intersent
static int doSubsetsIntersect(char* meshname1, int numMember1, int* memberIDs1,
			      char* meshname2, int numMember2, int* memberIDs2)
{
  int i, j, count = 0;;
  // mesh names have to match
  if (strcmp(meshname1, meshname2))
    return 0;
  // find the first matching ID
  i = 0;
  j = 0;
  while (i < numMember1 && j < numMember2) {
    if (memberIDs1[i] == memberIDs2[j]) {
      count++;
      i++;
      j++;
    }
    else if (memberIDs1[i] < memberIDs2[j])
      i++;
    else
      j++;
  }

  return count;
}

int main(int argc, char** argv) {
  FC_ReturnCode rc;
  FC_VerbosityLevel verbose_level = FC_QUIET;
  int i, j, k, m, n;
  char* aux_filenames[2] = { NULL, NULL };
  int doSpecificTear = 0, tearID;
  char* tear_name = NULL;
  xmlDocPtr xmlDocs[2];
  int numSubset[2], *maxNumMembers[2], *numMembers[2], **memberIDs[2];
  char **subsetNames[2], **meshNames[2];
  FC_AssociationType *assocs[2];
  int numTear[2], *numSubsetPerTear[2];
  char **tearNames[2], ***subsetNamesPerTear[2];
  double *tearLengths[2];
  int **tearCorrespond, **tearToSubsetMap[2];
  int count;

  // --- handle arguments

  if (argc < 3 ) {
  usage:
    printf("usage: %s [options] auxFile1 auxFile2\n", 
           argv[0]);
    printf("options: \n");
    printf("   -t \"tear name\"  : report results only for a specific tear from the 1st set\n");
    printf("                     of tears\n");
    printf("   -h              : print this help message\n");
    printf("   -v              : verbose: print warning and error messages\n");
    printf("   -V              : very verbose: prints log and error messages\n");
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
    else if (!strcmp(argv[i], "-t")) {
      doSpecificTear = 1;
      i++;
      tear_name = argv[i];
    } 
    else {
      if (i+1 > argc)
	goto usage;
      aux_filenames[0] = argv[i];
      aux_filenames[1] = argv[i+1];
      i++;
    }
  }

  // testing args
  if (!aux_filenames[0] || !aux_filenames[1])
    goto usage;

  // --- setup

  // set library verbosity (don't need to init library)
  rc = fc_setLibraryVerbosity(verbose_level);
  fc_exitIfError(rc);

  // --- try to import the two files
  for (i = 0; i < 2; i++) {
    rc = _fc_importAuxFileXMLDoc(aux_filenames[i], &xmlDocs[i]);
    fc_exitIfErrorPrintf(rc, "failed to import '%s'", aux_filenames[i]);
  }

  // --- get the subsets & tears & specific tear, if requested
  for (i = 0; i < 2; i++) {
    rc = _fc_getAuxFileSubsets(xmlDocs[i], &numSubset[i], &subsetNames[i],
			       &meshNames[i], &assocs[i], &maxNumMembers[i],
			       &numMembers[i], &memberIDs[i]);
    fc_exitIfError(rc);
    rc = _fc_getAuxFileTears(xmlDocs[i], &numTear[i], &tearNames[i],
			     &numSubsetPerTear[i], &subsetNamesPerTear[i],
			     &tearLengths[i]);
    fc_exitIfError(rc);
    if (i == 0 && doSpecificTear) {
      tearID = -1;
      for (j = 0; j < numTear[i]; j++) {
	if (!strcmp(tearNames[i][j], tear_name)) {
	  tearID = j;
	  break;
	}
      }
      if (tearID == -1) 
	fc_exitIfErrorPrintf(FC_ERROR, "Could not find tear '%s' in file '%s'",
			     tear_name, aux_filenames[i]);
    }
  }

  // --- done with xmlDocs, free them
  for (i = 0; i < 2; i++)
    xmlFreeDoc(xmlDocs[i]);

  // Match up subsets to tears
  // FIX? make more efficient?
  for (i = 0; i < 2; i++) {
    tearToSubsetMap[i] = (int**)malloc(numTear[i]*sizeof(int*));
    for (j = 0; j < numTear[i]; j++) {
      tearToSubsetMap[i][j] = (int*)malloc(numSubsetPerTear[i][j]*sizeof(int));
      for (k = 0; k < numSubsetPerTear[i][j]; k++) {
	for (m = 0; m < numSubset[i]; m++) {
	  if (!strcmp(subsetNames[i][m], subsetNamesPerTear[i][j][k])) {
	    tearToSubsetMap[i][j][k] = m;
	    break;
	  }
	}
      }
    }
  }

  // Match subsets between tears
  // FIX? make more efficient?
  // FIX? assumes sorted int array
  tearCorrespond = (int**)malloc(numTear[0]*sizeof(int*));
  if (!tearCorrespond)
    fc_exitIfError(FC_MEMORY_ERROR);
  for (i = 0; i < numTear[0]; i++) {
    tearCorrespond[i] = (int*)calloc(numTear[1], sizeof(int));
    if (!tearCorrespond[i])
      fc_exitIfError(FC_MEMORY_ERROR);
  }
  for (i = 0; i < numTear[0]; i++) {
    for (j = 0; j < numTear[1]; j++) {
      for (m = 0; m < numSubsetPerTear[0][i]; m++) {
	for (n = 0; n < numSubsetPerTear[1][j]; n++) {
	  int subsetID0 = tearToSubsetMap[0][i][m];
	  int subsetID1 = tearToSubsetMap[1][j][n];
	  tearCorrespond[i][j] = doSubsetsIntersect(meshNames[0][subsetID0],
			  numMembers[0][subsetID0], memberIDs[0][subsetID0],
			  meshNames[1][subsetID1], numMembers[1][subsetID1],
			  memberIDs[1][subsetID1]);
	  if (tearCorrespond[i][j] > 0)
	    break;
	}
	if (tearCorrespond[i][j] > 0)
	  break;
      }
    }
  }

  // --- report
  
  // header - the first dataset is the "known" tears that 2nd dataset tears
  // are being matched to
  if (doSpecificTear) 
    printf("Comparing specific tear\n");
  else
    printf("Comparing tears\n");
  printf("A = 1st set of tears: '%s'\n", aux_filenames[0]);
  printf("B = 2nd set of tears: '%s'\n", aux_filenames[1]);
  printf("A numTear = %d\n", numTear[0]);
  printf("B numTear = %d\n", numTear[1]);

  // do "known" tears from 1st set
  for (i = 0; i < numTear[0]; i++) {
    double totalTearLength = 0;

    // skip output if doing specific tear
    if (doSpecificTear && i != tearID)
      continue;

    count = 0;
    for (j = 0; j < numTear[1]; j++) {
      if (tearCorrespond[i][j] > 0) {
	count++;
	totalTearLength += tearLengths[1][j];
      }
    }
    printf("'%s' from A has %d corresponding tear(s) in B",
	   tearNames[0][i], count);
    if (count == 0)
      printf("\n");
    else {
      printf(":\n");
      for (j = 0; j < numTear[1]; j++) {
	if (tearCorrespond[i][j] > 0) 
	  printf("  '%s (numOverlap = %d, %.1f%%, %.1f%%)'\n", 
		 tearNames[1][j], tearCorrespond[i][j], 
		 ((double)tearCorrespond[i][j])/numMembers[0][i]*100,
		 ((double)tearCorrespond[i][j])/numMembers[1][j]*100);
      }
    }
    if (doSpecificTear)
      printf("Total length of tears on B = %g\n", totalTearLength);
  }
  count = 0; // count number that do correspond
  for (j = 0; j < numTear[1]; j++) {
    for (i = 0; i < numTear[0]; i++) {
      if (tearCorrespond[i][j] > 0) {
	count++;
	break;
      }
    }
  }
  // do "leftover" tears from second set
  if (!doSpecificTear) {
    count = numTear[1] - count; // calc number that do NOT correspond
    printf("%d tear(s) from B do not correspond to tears in A", 
	   count);
    if (count == 0)
    printf("\n");
    else {
      printf(":\n");
      for (j = 0; j < numTear[1]; j++) {
	count = 0;
	for (i = 0; i < numTear[0]; i++)
	  if (tearCorrespond[i][j] > 0) 
	    count++;
	if (count == 0) {
	  printf("  '%s'\n", tearNames[1][j]);
	}
      }
    }
  }

  // --- all done

  // cleanup
  for (i = 0; i < numTear[0]; i++)
    free(tearCorrespond[i]);
  free(tearCorrespond);
  for (i = 0; i < 2; i++) {
    for (j = 0; j < numTear[i]; j++)
      free(tearToSubsetMap[i][j]);
    free(tearToSubsetMap[i]);
  }
  for (i = 0; i < 2; i++) {
    for (j = 0; j < numSubset[i]; j++) {
      free(subsetNames[i][j]);
      free(meshNames[i][j]);
      free(memberIDs[i][j]);
    }
    free(subsetNames[i]);
    free(meshNames[i]);
    free(maxNumMembers[i]);
    free(numMembers[i]);
    free(memberIDs[i]);
    free(assocs[i]);
    for (j = 0; j < numTear[i]; j++) {
      free(tearNames[i][j]);
      for (k = 0; k < numSubsetPerTear[i][j]; k++)
	free(subsetNamesPerTear[i][j][k]);
      free(subsetNamesPerTear[i][j]);
    }
    free(tearNames[i]);
    free(subsetNamesPerTear[i]);
    free(numSubsetPerTear[i]);
    free(tearLengths[i]);
  }

  fc_finalLibrary();

  exit(0);
}
