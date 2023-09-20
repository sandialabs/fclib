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
 * \file checkcustom.c
 * \brief Unit testing of \ref Custom module
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/unittest/checkcustom.c,v $
 * $Revision: 1.14 $ 
 * $Date: 2006/08/30 19:20:05 $
 *
 */

#include <stdlib.h>
#include <check.h>
#include <math.h>
#include "fc.h"
#include "fcP.h"
#include "checkall.h"

// **** test fixtures

static void custom_setup(void) {
  FC_ReturnCode rc;
  if (isForking) {
    fc_setLibraryVerbosity(fc_messages);
    rc = fc_initLibrary();
    fail_unless(rc == FC_SUCCESS, "test aborted: failed to init library");
  }
}

static void custom_teardown(void) {
  FC_ReturnCode rc;
  if (isForking) {
    rc = fc_finalLibrary();
    fail_unless(rc == FC_SUCCESS, "failed to final library at end of test");
  }
}

START_TEST(simple_test)
{

  fail_unless(1, "should always pass this test");

}
END_TEST

// Populate the Suite with the tests

Suite *custom_suite(void)
{
  Suite *suite = suite_create("Custom");

  // Example doesn't really have any tests, it's just an example
  // to remind how to set them up
  TCase *tc_example = tcase_create(" - Example ");
  TCase *tc_realtests = tcase_create(" - Real Tests ");

  suite_add_tcase(suite, tc_example);
  tcase_add_checked_fixture(tc_example, custom_setup, custom_teardown);
  tcase_add_test(tc_example, simple_test);

  suite_add_tcase(suite, tc_realtests);
  tcase_add_checked_fixture(tc_realtests, custom_setup, custom_teardown);
  //tcase_add_test(tc_tealtests, test);

  return suite;
}
