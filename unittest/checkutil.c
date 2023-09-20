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
 * \file checkutil.c
 * \brief Unit tests for the \ref Utilities module.
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/unittest/checkutil.c,v $
 * $Revision: 1.37 $ 
 * $Date: 2006/08/30 19:20:05 $
 *
 * \modifications
 *   - 05/07/04  WSK, moved linked list and mask array tests to checkbase.c 
 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <float.h>
#include <math.h>

#include <check.h>

#include "base.h"
#include "util.h"
#include "checkall.h"

// **** utility tests

// test FC_VALUE_EQUIV, FC_FLT_EQUIV, FC_DBL_EQUIV
START_TEST(value_equiv)
{
  FC_ReturnCode rc;
  float zeroF = 0.F, oneF = 1.F;
  float aboutZeroF = FLT_MIN/10.F, almostZeroF = FLT_MIN;
  float aboutOneF = 1.F + FLT_EPSILON/2.F, almostOneF = 1.F + 2.F*FLT_EPSILON;
  double zeroD = 0., oneD = 1.;
  double aboutZeroD = DBL_MIN/10., almostZeroD = DBL_MIN;
  double aboutOneD = 1. + DBL_EPSILON/2., almostOneD = 1. + 2.*DBL_EPSILON;

  // *** general value equiv tests

  // positive tests close to zero
  fail_unless(FC_VALUE_EQUIV(zeroF, zeroF, FLT_EPSILON, FLT_MIN), 
	      "failed to find 0. = 0.");
  fail_unless(FC_VALUE_EQUIV(zeroF, aboutZeroF, FLT_EPSILON, FLT_MIN),
	      "failed to find zeroF = aboutZeroF");
  fail_unless(FC_VALUE_EQUIV(aboutZeroF, zeroF, FLT_EPSILON, FLT_MIN),
	      "failed to find aboutZeroF = zeroF ");
  // negative tests close to zero
  fail_unless(!FC_VALUE_EQUIV(almostZeroF, zeroF, FLT_EPSILON, FLT_MIN),
	      "failed to find almostZeroF different from zeroF");
  fail_unless(!FC_VALUE_EQUIV(zeroF, almostZeroF, FLT_EPSILON, FLT_MIN),
	      "failed to find zeroF different from almostZeroF");
  // positive tests away from zero
  fail_unless(FC_VALUE_EQUIV(oneF, oneF, FLT_EPSILON, FLT_MIN),
              "failed to find 1. = 1. ");
  fail_unless(FC_VALUE_EQUIV(aboutOneF, oneF, FLT_EPSILON, FLT_MIN),
              "failed to find 1. + FLT_EPSILON/2 = 1.");
  fail_unless(FC_VALUE_EQUIV(oneF, aboutOneF, FLT_EPSILON, FLT_MIN),
              "failed to find 1. = 1. + FLT_EPSILON/2.");
  // negative tests away from zero
  fail_unless(!FC_VALUE_EQUIV(almostOneF, oneF, FLT_EPSILON, FLT_MIN),
              "failed to find 1. + FLT_EPSILON*2 different from 1.");
  fail_unless(!FC_VALUE_EQUIV(oneF, almostOneF, FLT_EPSILON, FLT_MIN),
              "failed to find 1. different from 1. + FLT_EPSILON*2.");

  // what if we do ints!
  fail_unless(FC_VALUE_EQUIV(1, 1, 1, 0.1), "using ints shouldn't hurt us");
  // strange things
  rc = FC_VALUE_EQUIV('a', 'b', 3, 4);
  fail_unless(rc != 0 || rc != 1, "character input shouldn't hurt us");

  // *** general float equiv tests

  // positive tests close to zero
  fail_unless(FC_FLT_EQUIV(zeroF, zeroF), "failed to find 0. = 0.");
  fail_unless(FC_FLT_EQUIV(zeroF, aboutZeroF),
              "failed to find zeroF = aboutZeroF");
  fail_unless(FC_FLT_EQUIV(aboutZeroF, zeroF),
              "failed to find aboutZeroF = zeroF");
  // negative tests close to zero
  fail_unless(!FC_FLT_EQUIV(almostZeroF, zeroF),
              "failed to find almostZeroF different from zeroF");
  fail_unless(!FC_FLT_EQUIV(zeroF, almostZeroF),
              "failed to find zeroF different from almostZeroF");
  // positive tests away from zero
  fail_unless(FC_FLT_EQUIV(oneF, oneF),
              "failed to find oneF = oneF");
  fail_unless(FC_FLT_EQUIV(aboutOneF, oneF),
              "failed to find aboutOneF = oneF");
  fail_unless(FC_FLT_EQUIV(oneF, aboutOneF),
              "failed to find oneF = aboutOneF");
  // negative tests away from zero
  fail_unless(!FC_FLT_EQUIV(almostOneF, oneF),
              "failed to find almostOneF different from oneF");
  fail_unless(!FC_FLT_EQUIV(oneF, almostOneF),
              "failed to find oneF different from almostOneF");
  
  // *** general double equiv tests

  // positive tests close to zero
  fail_unless(FC_DBL_EQUIV(zeroD, zeroD), "failed to find 0. = 0.");
  fail_unless(FC_DBL_EQUIV(zeroD, aboutZeroD),
              "failed to find zeroD = aboutZeroD");
  fail_unless(FC_DBL_EQUIV(aboutZeroD, zeroD),
              "failed to find aboutZeroD = zeroD");
  // negative tests close to zero
  fail_unless(!FC_DBL_EQUIV(almostZeroD, zeroD),
              "failed to find almostZeroD different from zeroD");
  fail_unless(!FC_DBL_EQUIV(zeroD, almostZeroD),
              "failed to find zeroD different from almostZeroD");
  // positive tests away from zero
  fail_unless(FC_DBL_EQUIV(oneD, oneD),
              "failed to find oneD = oneD");
  fail_unless(FC_DBL_EQUIV(aboutOneD, oneD),
              "failed to find aboutOneD = oneD");
  fail_unless(FC_DBL_EQUIV(oneD, aboutOneD),
              "failed to find oneD = aboutOneD");
  // negative tests away from zero
  fail_unless(!FC_DBL_EQUIV(almostOneD, oneD),
              "failed to find almostOneD different from oneD");
  fail_unless(!FC_DBL_EQUIV(oneD, almostOneD),
              "failed to find oneD different from almostOneD");

  // **** tests added because of bug for cases with values far apart but both
  //      less than MIN
  fail_unless(FC_VALUE_EQUIV(1.1, 1.2, 1.e-6, 2.), 
	      "failed to find values both less than MIN the same");
  fail_unless(FC_FLT_EQUIV(FLT_MIN/10., FLT_MIN/20.), 
	      "failed to find close to zero the same as even closer to zero");
  fail_unless(FC_DBL_EQUIV(DBL_MIN/10., DBL_MIN/20.), 
	      "failed to find close to zero the same as even closer to zero");

  // **** test for divide by zero possible bug introduced after fixing
  //      previous bug
  fail_unless(!FC_VALUE_EQUIV(0, 1, FLT_EPSILON, FLT_MIN), 
	      "0 should not equal 1");
  fail_unless(!FC_VALUE_EQUIV(1, 0, FLT_EPSILON, FLT_MIN), 
	      "1 should not equal 0");
  fail_unless(!FC_FLT_EQUIV(0, 1), "0 should not equal 1");
  fail_unless(!FC_FLT_EQUIV(1, 0), "1 should not equal 0");
  fail_unless(!FC_DBL_EQUIV(0, 1), "0 should not equal 1");
  fail_unless(!FC_DBL_EQUIV(1, 0), "1 should not equal 0");

  // **** test for another divide by zero bug possible if non zero entry was
  //      negative
  fail_unless(!FC_VALUE_EQUIV(0, -1, FLT_EPSILON, FLT_MIN), 
	      "0 should not equal 1");
  fail_unless(!FC_VALUE_EQUIV(-1, 0, FLT_EPSILON, FLT_MIN), 
	      "1 should not equal 0");
  fail_unless(!FC_FLT_EQUIV(0, -1), "0 should not equal 1");
  fail_unless(!FC_FLT_EQUIV(-1, 0), "1 should not equal 0");
  fail_unless(!FC_DBL_EQUIV(0, -1), "0 should not equal 1");
  fail_unless(!FC_DBL_EQUIV(-1, 0), "1 should not equal 0");
}
END_TEST

// test fc_eq, fc_neq, fc_lt, fc_lteq, fc_gt, fc_gteq
// just testing float values for above should be adequate
// test fc_eqf, fc_neqf, fc_ltf, fc_lteqf, fc_gtf, fc_gteqf
// test fc_eqd, fc_neqd, fc_ltd, fc_lteqd, fc_gt, fc_gteqd
// some D tests shouldn't pass w/ F test values!
START_TEST(comparisons)
{
  int i, j;
  float fvals1[3], fvals2[3];
  double dvals1[3], dvals2[3];
  float aboutZeroF = FLT_MIN/10.F, almostZeroF = FLT_MIN;
  float aboutOneF = 1.F + FLT_EPSILON/2.F, almostOneF = 1.F + 2.F*FLT_EPSILON;
  double aboutZeroD = DBL_MIN/10., almostZeroD = DBL_MIN;
  double aboutOneD = 1. + DBL_EPSILON/2., almostOneD = 1. + 2.*DBL_EPSILON;

  // *** compare equal values near zero
  fvals1[0] = fvals2[0] = 0;
  fvals1[1] = fvals2[1] = aboutZeroF;
  dvals1[0] = dvals2[0] = 0;
  dvals1[1] = dvals2[1] = aboutZeroD;
  for (i = 0; i < 2; i++) {
    for (j = 0; j < 2; j++) {
      // eq
      fail_unless(fc_eq(fvals1[i], fvals2[j], FLT_EPSILON, FLT_MIN), 
		  "should be true");
      fail_unless(fc_eqf(fvals1[i], fvals2[j]), "should be true");
      fail_unless(fc_eqd(dvals1[i], dvals2[j]), "should be true");
      // neq
      fail_unless(!fc_neq(fvals1[i], fvals2[j], FLT_EPSILON, FLT_MIN), 
		  "should be false");
      fail_unless(!fc_neqf(fvals1[i], fvals2[j]), "should be false");
      fail_unless(!fc_neqd(dvals1[i], dvals2[j]), "should be false");
      // lt
      fail_unless(!fc_lt(fvals1[i], fvals2[j], FLT_EPSILON, FLT_MIN), 
		  "should be false");
      fail_unless(!fc_ltf(fvals1[i], fvals2[j]), "should be false");
      fail_unless(!fc_ltd(dvals1[i], dvals2[j]), "should be false");
      // lteq
      fail_unless(fc_lteq(fvals1[i], fvals2[j], FLT_EPSILON, FLT_MIN), 
		  "should be true");
      fail_unless(fc_lteqf(fvals1[i], fvals2[j]), "should be true");
      fail_unless(fc_lteqd(dvals1[i], dvals2[j]), "should be true");
      // gt
      fail_unless(!fc_gt(fvals1[i], fvals2[j], FLT_EPSILON, FLT_MIN), 
		  "should be false");
      fail_unless(!fc_gtf(fvals1[i], fvals2[j]), "should be false");
      fail_unless(!fc_gtd(dvals1[i], dvals2[j]), "should be false");
      // gteq
      fail_unless(fc_gteq(fvals1[i], fvals2[j], FLT_EPSILON, FLT_MIN), 
		  "should be true");
      fail_unless(fc_gteqf(fvals1[i], fvals2[j]), "should be true");
      fail_unless(fc_gteqd(dvals1[i], dvals2[j]), "should be true");
    }
  }

  // *** compare equal values away from zero
  fvals1[0] = fvals2[0] = 1;
  fvals1[1] = fvals2[1] = aboutOneF;
  dvals1[0] = dvals2[0] = 1;
  dvals1[1] = dvals2[1] = aboutOneD;
  for (i = 0; i < 2; i++) {
    for (j = 0; j < 2; j++) {
      // eq
      fail_unless(fc_eq(fvals1[i], fvals2[j], FLT_EPSILON, FLT_MIN), 
		  "should be true");
      fail_unless(fc_eqf(fvals1[i], fvals2[j]), "should be true");
      fail_unless(fc_eqd(dvals1[i], dvals2[j]), "should be true");
      // neq
      fail_unless(!fc_neq(fvals1[i], fvals2[j], FLT_EPSILON, FLT_MIN), 
		  "should be false");
      fail_unless(!fc_neqf(fvals1[i], fvals2[j]), "should be false");
      fail_unless(!fc_neqd(dvals1[i], dvals2[j]), "should be false");
      // lt
      fail_unless(!fc_lt(fvals1[i], fvals2[j], FLT_EPSILON, FLT_MIN), 
		  "should be false");
      fail_unless(!fc_ltf(fvals1[i], fvals2[j]), "should be false");
      fail_unless(!fc_ltd(dvals1[i], dvals2[j]), "should be false");
      // lteq
      fail_unless(fc_lteq(fvals1[i], fvals2[j], FLT_EPSILON, FLT_MIN), 
		  "should be true");
      fail_unless(fc_lteqf(fvals1[i], fvals2[j]), "should be true");
      fail_unless(fc_lteqd(dvals1[i], dvals2[j]), "should be true");
      // gt
      fail_unless(!fc_gt(fvals1[i], fvals2[j], FLT_EPSILON, FLT_MIN), 
		  "should be false");
      fail_unless(!fc_gtf(fvals1[i], fvals2[j]), "should be false");
      fail_unless(!fc_gtd(dvals1[i], dvals2[j]), "should be false");
      // gteq
      fail_unless(fc_gteq(fvals1[i], fvals2[j], FLT_EPSILON, FLT_MIN), 
		  "should be true");
      fail_unless(fc_gteqf(fvals1[i], fvals2[j]), "should be true");
      fail_unless(fc_gteqd(dvals1[i], dvals2[j]), "should be true");
    }
  }

  // *** compare smaller value to bigger value
  fvals1[0] = 0;
  fvals2[0] = 1;
  fvals1[1] = 0;
  fvals2[1] = almostZeroF;
  fvals1[2] = 1;
  fvals2[2] = almostOneF;
  dvals1[0] = 0;
  dvals2[0] = 1;
  dvals1[1] = 0;
  dvals2[1] = almostZeroD;
  dvals1[2] = 1;
  dvals2[2] = almostOneD;
  for (i = 0; i < 3; i++) {
    // eq
    fail_unless(!fc_eq(fvals1[i], fvals2[i], FLT_EPSILON, FLT_MIN), 
		"should be false");
    fail_unless(!fc_eqf(fvals1[i], fvals2[i]), "should be false");
    fail_unless(!fc_eqd(dvals1[i], dvals2[i]), "should be false");
    // neq
    fail_unless(fc_neq(fvals1[i], fvals2[i], FLT_EPSILON, FLT_MIN), 
		"should be true");
    fail_unless(fc_neqf(fvals1[i], fvals2[i]), "should be true");
    fail_unless(fc_neqd(dvals1[i], dvals2[i]), "should be true");
    // lt
    fail_unless(fc_lt(fvals1[i], fvals2[i], FLT_EPSILON, FLT_MIN), 
		"should be true");
    fail_unless(fc_ltf(fvals1[i], fvals2[i]), "should be true");
    fail_unless(fc_ltd(dvals1[i], dvals2[i]), "should be true");
    // lteq
    fail_unless(fc_lteq(fvals1[i], fvals2[i], FLT_EPSILON, FLT_MIN), 
		"should be true");
    fail_unless(fc_lteqf(fvals1[i], fvals2[i]), "should be true");
    fail_unless(fc_lteqd(dvals1[i], dvals2[i]), "should be true");
    // gt
    fail_unless(!fc_gt(fvals1[i], fvals2[i], FLT_EPSILON, FLT_MIN), 
		"should be false");
    fail_unless(!fc_gtf(fvals1[i], fvals2[i]), "should be false");
    fail_unless(!fc_gtf(dvals1[i], dvals2[i]), "should be false");
    // gteq
    fail_unless(!fc_gteq(fvals1[i], fvals2[i], FLT_EPSILON, FLT_MIN), 
		"should be false");
    fail_unless(!fc_gteqf(fvals1[i], fvals2[i]), "should be false");
    fail_unless(!fc_gteqd(dvals1[i], dvals2[i]), "should be false");
  }

  // *** compare bigger value to smaller value
  fvals1[0] = 1;
  fvals2[0] = 0;
  fvals1[1] = almostZeroF;
  fvals2[1] = 0;
  fvals1[2] = almostOneF;
  fvals2[2] = 1;
  dvals1[0] = 1;
  dvals2[0] = 0;
  dvals1[1] = almostZeroD;
  dvals2[1] = 0;
  dvals1[2] = almostOneD;
  dvals2[2] = 1;
  for (i = 0; i < 3; i++) {
    // eq
    fail_unless(!fc_eq(fvals1[i], fvals2[i], FLT_EPSILON, FLT_MIN), 
		"should be false");
    fail_unless(!fc_eqf(fvals1[i], fvals2[i]), "should be false");
    fail_unless(!fc_eqd(dvals1[i], dvals2[i]), "should be false");
    // neq
    fail_unless(fc_neq(fvals1[i], fvals2[i], FLT_EPSILON, FLT_MIN), 
		"should be true");
    fail_unless(fc_neqf(fvals1[i], fvals2[i]), "should be true");
    fail_unless(fc_neqd(dvals1[i], dvals2[i]), "should be true");
    // lt
    fail_unless(!fc_lt(fvals1[i], fvals2[i], FLT_EPSILON, FLT_MIN), 
		"should be false");
    fail_unless(!fc_ltf(fvals1[i], fvals2[i]), "should be false");
    fail_unless(!fc_ltd(dvals1[i], dvals2[i]), "should be false");
    // lteq
    fail_unless(!fc_lteq(fvals1[i], fvals2[i], FLT_EPSILON, FLT_MIN), 
		"should be false");
    fail_unless(!fc_lteqf(fvals1[i], fvals2[i]), "should be false");
    fail_unless(!fc_lteqd(dvals1[i], dvals2[i]), "should be false");
    // gt
    fail_unless(fc_gt(fvals1[i], fvals2[i], FLT_EPSILON, FLT_MIN), 
		"should be true");
    fail_unless(fc_gtf(fvals1[i], fvals2[i]), "should be true");
    fail_unless(fc_gtd(dvals1[i], dvals2[i]), "should be true");
    // gteq
    fail_unless(fc_gteq(fvals1[i], fvals2[i], FLT_EPSILON, FLT_MIN), 
		"should be true");
    fail_unless(fc_gteqf(fvals1[i], fvals2[i]), "should be true");
    fail_unless(fc_gteqd(dvals1[i], dvals2[i]), "should be true");
  }
}
END_TEST

// int array

START_TEST(sort_int_array) 
{
  FC_ReturnCode rc;
  int i;
  int num = 6;
  int unsorted[6] = { 25, 0, -1, -21, 0, 2 };
  int sorted[6] = { -21, -1, 0, 0, 2, 25 };

  // test bad inputs
  rc = fc_sortIntArray(0, unsorted);
  fail_unless(rc != FC_SUCCESS, "should fail to sort 0 entries");
  rc = fc_sortIntArray(num, NULL);
  fail_unless(rc != FC_SUCCESS, "should fail to sort NULL array");
  
  // test
  fc_sortIntArray(num, unsorted);
  for (i = 0; i < num; i++)
    fail_unless(unsorted[i] == sorted[i], "mismatch of ints");
}
END_TEST

START_TEST(grow_int_array)
{
  FC_ReturnCode rc;
  int i;
  int num = 0, numValue = 6;
  int *array, values[6] = {1, 2, 3, 245, -1, 7};

  // add to an empty, but nonnull array
  rc = fc_addValueToIntArray(&array, &num, values[0]);
  fail_unless(rc == FC_SUCCESS, "should work");
  fail_unless(num == 1, "num should be 1 after adding first entry");
  fail_unless(array[0] == values[0], 
	      "1st array value should equal added value");

  // add to an empty, null array
  free(array);
  array = NULL;
  num = 0;
  rc =fc_addValueToIntArray(&array, &num, values[0]);
  fail_unless(rc == FC_SUCCESS, "should work");
  fail_unless(num == 1, 
	      "num should be 1 after adding first entry to null array");
  fail_unless(array[0] == values[0], "1st array value wrong in null array");

  // add remaining
  for (i = 1; i < numValue; i++) {
    rc = fc_addValueToIntArray(&array, &num, values[i]);
    fail_unless(rc == FC_SUCCESS, "should work");
  }
  fail_unless(num == numValue, "array doesn not have proper # of entries");
  fail_unless(!memcmp(array, values, sizeof(int)*numValue), 
              "arrays are not equal");

  // test bad args
  rc = fc_addValueToIntArray(NULL, &num, values[i]);
  fail_unless(rc != FC_SUCCESS, "should work");
  fail_unless(num == numValue, "array doesn not have proper # of entries");
  fail_unless(!memcmp(array, values, sizeof(int)*numValue), 
              "arrays are not equal");
  rc = fc_addValueToIntArray(&array, NULL, values[i]);
  fail_unless(rc != FC_SUCCESS, "should work");
  fail_unless(num == numValue, "array doesn not have proper # of entries");
  fail_unless(!memcmp(array, values, sizeof(int)*numValue), 
              "arrays are not equal");

  free(array);
  
}
END_TEST

//******* path extraction routines

// Assemble the different cases automatically
#define NUMDIR 3
#define NUMBASE 2
#define NUMEXT 3
#define NUMSLASH 2
static char* g_dirs[NUMDIR] = { "/usr/local", "local", NULL };
static char* g_bases[NUMBASE] = { "bin", ".bin" };
static char* g_exts[NUMEXT] = { "", "ext", "mor" };
static char* g_slashes[NUMSLASH] = { NULL, "/" };

START_TEST(path_extracts)
{
  int i, j, k, m, n;
  char path[1024], path_keep[1024];
  char buffer[1024];
  char answer[1024];
  char* temp_str;

  for (i = 0; i < NUMDIR; i++) { // dirnames
    for (j = 0; j < NUMBASE; j++) { // basenames
      for (k = 0; k < NUMEXT; k++) { // w or w/o extension
	for (m = 0; m < NUMSLASH; m++) { // w or w/o trailing slash
	  // assemble path
	  path[0] = '\0';
	  if (g_dirs[i]) {
	    sprintf(buffer, "%s/", g_dirs[i]);
	    strcat(path, buffer);
	  }
	  strcat(path, g_bases[j]);
	  if (k == 1) {
	    sprintf(buffer, ".%s", g_exts[1]);
	    strcat(path, buffer);
	  }
	  else if (k == 2) {
	    sprintf(buffer, ".%s.%s", g_exts[1], g_exts[2]);
	    strcat(path, buffer);
	  }
	  if (g_slashes[m])
	    strcat(path, g_slashes[m]);
	  strcpy(path_keep, path);
	  //printf("\npath = %s\n", path);

	  // dirname
	  temp_str = fc_getDirname(path);
	  if (g_dirs[i])
	    fail_unless(!strcmp(temp_str, g_dirs[i]), "mismatch of dirname");
	  else
	    fail_unless(!strcmp(temp_str, "."), "mismatch of dirname");
	  free(temp_str);
	  fail_unless(!strcmp(path, path_keep), "path was altered!");

	  // basename
	  temp_str = fc_getBasename(path);
	  if (k == 0) 
	    strcpy(answer, g_bases[j]);
	  else if (k == 1)
	    sprintf(answer, "%s.%s", g_bases[j], g_exts[1]);
	  else if (k == 2) 
	    sprintf(answer, "%s.%s.%s", g_bases[j], g_exts[1], g_exts[2]);
	  fail_unless(!strcmp(temp_str, answer), "mismatch of basename");
	  free(temp_str);
	  fail_unless(!strcmp(path, path_keep), "path was altered!");

	  // basename w/o extension
	  for (n = 0; n < 2; n++) {
	    temp_str = fc_getBasenameWOExtension(path, n);
	    if (k < 2 || n == 0)
	      fail_unless(!strcmp(temp_str, g_bases[j]), 
			  "mismatch of basename w/o ext");
	    else {
	      sprintf(answer, "%s.%s", g_bases[j], g_exts[1]);
	      fail_unless(!strcmp(temp_str, answer),
			  "mismatch of basename w/o ext");
	    }
	    free(temp_str);
	    fail_unless(!strcmp(path, path_keep), "path was altered!");
	  }

	  // extension
	  for (n = 0; n < 2; n++) {
	    temp_str = fc_getExtension(path, n);
	    if (k <= 1)
	      fail_unless(!strcmp(temp_str, g_exts[k]), "mismatch of ext");
	    else if (n == 0) {
	      sprintf(answer, "%s.%s", g_exts[1], g_exts[2]);
	      fail_unless(!strcmp(temp_str, answer), "mismatchof extension");
	    }
	    else
	      fail_unless(!strcmp(temp_str, g_exts[2]), "mismatch of ext");
	    free(temp_str);
	    fail_unless(!strcmp(path, path_keep), "path was altered!");
	  }

	  // check what happens if NULLs passed in
	  temp_str = fc_getDirname(NULL);
	  fail_unless(temp_str == NULL, "mismtach of NULL");
	  temp_str = fc_getBasename(NULL);
	  fail_unless(temp_str == NULL, "mismtach of NULL");
	  temp_str = fc_getBasenameWOExtension(NULL, 0);
	  fail_unless(temp_str == NULL, "mismtach of NULL");
	  temp_str = fc_getExtension(NULL, 0);
	  fail_unless(temp_str == NULL, "mismtach of NULL");
	}
      }
    }
  }
}
END_TEST

START_TEST(test_pi)
{
  fail_unless(FC_DBL_EQUIV(acos(0)*2, FC_PI),
	      "value of FC_PI does not agreed w/ acos cals");
}
END_TEST

START_TEST(replace_chars)
{
  char in_string[100] = "bill bull byll";
  char out_string[100] = "bill_bull_byll";
  char* temp_string;

  // --- test good input

  // case 1: change spaces to underscores
  temp_string = fc_replaceChars(in_string, ' ', '_');
  fail_unless(temp_string != NULL, "failed to replace spaces");
  fail_unless(!strcmp(temp_string, out_string), "mismatch");
  free(temp_string);
  // case 2: char doesn't exist in string - no change
  temp_string = fc_replaceChars(in_string, '&', '_');
  fail_unless(temp_string != NULL, "failed to replace spaces");
  fail_unless(!strcmp(temp_string, in_string), "shouldn't change");
  free(temp_string);
 
  // --- test bad input
  
  temp_string = fc_replaceChars(NULL, ' ', '_');
  fail_unless(temp_string == NULL, "null input should return null output");
}
END_TEST

// ***** Populate the Suite with the tests

Suite *util_suite(void)
{
  Suite *suite = suite_create("Utilities");

  TCase *tc_intArray = tcase_create(" - Int Arrays ");
  TCase *tc_floats = tcase_create(" - Floating Point ");
  TCase *tc_path = tcase_create(" - Path Processing");
  TCase *tc_misc = tcase_create(" - Misc Utilities");

  suite_add_tcase(suite, tc_floats);
  tcase_add_test(tc_floats, value_equiv);
  tcase_add_test(tc_floats, comparisons);

  suite_add_tcase(suite, tc_intArray);
  tcase_add_test(tc_intArray, grow_int_array);
  tcase_add_test(tc_intArray, sort_int_array);

  suite_add_tcase(suite, tc_path);
  tcase_add_test(tc_path, path_extracts);

  suite_add_tcase(suite, tc_misc);
  tcase_add_test(tc_misc, test_pi);
  tcase_add_test(tc_misc, replace_chars);

  return suite;
}
