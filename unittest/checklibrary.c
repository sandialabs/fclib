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
 * \file checklibrary.c
 * \brief Unit tests for \ref Library module.
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/unittest/checklibrary.c,v $
 * $Revision: 1.7 $ 
 * $Date: 2006/08/30 19:20:05 $
 *
 * \modifications
 *    3/30/04 WSK, added tests for library verbosity
 */

#include <stdlib.h>
#include <check.h>
#include "fc.h"
#include "fcP.h"
#include "checkall.h"

// *********************************************
// ***** General library interface tests
// *********************************************

// Normal init & then final
START_TEST(init_final)
{
  FC_ReturnCode rc;

  if (isForking) {
    fc_setLibraryVerbosity(fc_messages);

    // o.k. to final uninited library
    rc = fc_finalLibrary();
    fail_unless(rc == FC_SUCCESS, "should be o.k. to final uninited library");
    
    // normal init
    rc = fc_initLibrary();
    fail_unless(rc == FC_SUCCESS, "FCLIB did not init normally");
    
    // bad, can't reinit init library
    rc = fc_initLibrary();
    fail_unless(rc != FC_SUCCESS, "second init should have failed");

    // normal final
    rc = fc_finalLibrary();
    fail_unless(rc == FC_SUCCESS, "FCLIB did not final normally");

    // o.k. to final a finaled library
    rc = fc_finalLibrary();
    fail_unless(rc == FC_SUCCESS, "should be o.k. to final a finaled library");
    
    // bad, can't reinit a finaled library
    rc = fc_initLibrary();
    fail_unless(rc != FC_SUCCESS, "can't reinit a finaled library");
  }

  else 
    printf("**CK_NOFORK: Can't test out of order fc_initLibrary()/fc_finalLibrary() calls\n");
}
END_TEST

START_TEST(get_set_verbosity)
{
  // Note, not setting to LOG_MESSAGES because it clutters test output
  // the routines are simple enough that can adequately test with other values
  FC_ReturnCode rc;
  int verbose, default_verbose = FC_QUIET;
  char message[1024];

  if (isForking) {
    // tests before initing library-- should be o.k.
    rc = fc_setLibraryVerbosity(FC_QUIET);
    fail_unless(rc == FC_SUCCESS, "failed to set verbosity on uninitend libary");
    verbose = fc_getLibraryVerbosity();
    fail_unless(verbose == default_verbose, 
		"failed to get correct verbosity from uninitend libary");
    
    // init library  
    rc = fc_initLibrary();
    fail_unless(rc == FC_SUCCESS, "test aborted: fclib failed to init");
  }
  else
    printf("**CK_NOFORK: Can't test setting verbosity before library inited\n");

  // test with init library
  verbose = fc_getLibraryVerbosity();
  sprintf(message, "default library verbosity should be %s", 
	  fc_getVerbosityLevelText(default_verbose));
  fail_unless(verbose == default_verbose, message);
  rc = fc_setLibraryVerbosity(-99);
  fail_unless(rc != FC_SUCCESS, "should fail to set nonsense verbosity");
  verbose = fc_getLibraryVerbosity();
  fail_unless(verbose == FC_QUIET, "library verbosity should still be default");
  rc = fc_setLibraryVerbosity(FC_ERROR_MESSAGES);
  fail_unless(rc == FC_SUCCESS, "should not fail to set verbosity");
  verbose = fc_getLibraryVerbosity();
  fail_unless(verbose == FC_ERROR_MESSAGES, 
	      "verbosity should have changed to set value");

  if (isForking) {
    // final libray
    rc = fc_finalLibrary();
    fail_unless(rc == FC_SUCCESS,"test aborted: fclib failed to final");
    
    // tests after finaling library--should be o.k.
    rc = fc_setLibraryVerbosity(FC_QUIET);
    fail_unless(rc == FC_SUCCESS, "failed to set verbosity on finaled libary");
    verbose = fc_getLibraryVerbosity();
    fail_unless(verbose == FC_QUIET, "verbosity should have changed to set value");
  }
  else {
    fc_setLibraryVerbosity(FC_QUIET);
    printf("**CK_NOFORK: Can't test setting verbosity after library finaled\n");
  }
}
END_TEST
  
// *********************************************
// ***** Populate the Suite with the tests
// *********************************************

Suite *library_suite(void)
{
  Suite *suite = suite_create("Library");

  TCase *tc_library = tcase_create(" - Library Interface ");

  // general library init/final tests
  suite_add_tcase(suite, tc_library);
  tcase_add_test(tc_library, init_final);
  tcase_add_test(tc_library, get_set_verbosity);

  return suite;
}
