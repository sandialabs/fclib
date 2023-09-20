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
 * \file checktable.c
 * \brief Unit tests for \ref PrivateTable module.
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/unittest/checktable.c,v $
 * $Revision: 1.59 $ 
 * $Date: 2006/08/30 19:20:05 $
 *
 * \modifications
 *    3/30/04 WSK, added tests for library verbosity
 */

#include <stdlib.h>
#include <string.h>
#include <check.h>
#include "fc.h"
#include "fcP.h"
#include "checkall.h"

// *********************************************
// ***** General library interface tests
// *********************************************

// test slot header interface
START_TEST(slot_header)
{
  FC_ReturnCode rc;
  _FC_SlotHeader header;
  int slotID = 25, uID = 95;
  char name[10] = "hello";

  // fill header with stupid values
  header.slotID = slotID;
  header.uID = uID;
  header.name = (char*)1;
  header.committed = -99;

  // test init
  _fc_initSlotHeader(&header);
  fail_unless(header.slotID == slotID && header.uID == uID,
	      "init should NOT change slotID or uID");
  fail_unless(header.name == NULL && header.committed ==  0, 
	      "init should change flags to 0");

  // fill flags with real values  
  header.committed = 1;

  // test adding name -- bad arguments
  rc = _fc_setSlotHeaderName(NULL, name);
  fail_unless(rc != FC_SUCCESS, "should fail with NULL header");
  rc = _fc_setSlotHeaderName(&header, NULL);
  fail_unless(rc != FC_SUCCESS, "should fail with NULL name");  

  // test adding name -- good arguments
  rc = _fc_setSlotHeaderName(&header, name);
  fail_unless(rc == FC_SUCCESS, "failed to set name");
  fail_unless(!strcmp(header.name, name), "name mismatch");
  fail_unless(header.slotID == slotID && header.uID == uID,
	      "setName should NOT change slotID or uID");
  fail_unless(header.committed ==  1, 
	      "setName should NOTE change flags");

  // test clear
  _fc_clearSlotHeader(&header);
  fail_unless(header.slotID == slotID, "clear should NOT change slotID");
  fail_unless(header.uID == -1, "clear should change uID to -1");
  fail_unless(header.name == NULL && header.committed ==  0, 
	      "clear should change flags to 0");
}
END_TEST
  
// test stable interface
// Since these are very private routines, not testing bad arguments
START_TEST(table)
{
  FC_ReturnCode rc;
  int i, j;
  _FC_SlotHeader **table1 = NULL, **table2 = NULL, *slot;
  int size1 = 0, size2 = 0, temp_num;
  FC_SortedIntArray openSlots1 = { 0, 0, 0 };
  FC_SortedIntArray openSlots2 = { 0, 0, 0 };
  int uID = 0;
  int num = 20;   // make sure loop captures a couple of doublings
  // generic handle & object for testing (like FC_Mesh & _FC_MeshSlot)
  struct { int slotID; int uID; } handle;
  struct { _FC_SlotHeader header; } object; 

  //Make a bunch of of new slots in two different tables
  for (i = 0; i < num; i++) {
    for (j = 0; j < 2; j++) {
      slot = NULL;
      if (j == 0) 
	slot = _fc_getNewSlot(&size1, &table1, &openSlots1,
			      sizeof(_FC_SlotHeader));
      else
	slot = _fc_getNewSlot(&size2, &table2, &openSlots2, 
			      sizeof(_FC_SlotHeader));
      uID++;
      fail_unless(slot != NULL, "didn't get a new slot");
      fail_unless(slot->slotID == i, "should be first slot");
      fail_unless(slot->uID == uID, "should have next uID");
      fail_unless(slot->name == NULL && slot->committed ==  0, 
		  "newSlot should have 0s");
    }
  }

  // Test openSlots part of data struture (only need to check 1 since the same)
  for (i = 0; i < openSlots1.numVal; i++) 
    fail_unless(openSlots1.vals[i] == i+num, "mismatch of open slots");

  // Test _fc_getSlot
  for (i = 0; i < num; i++) {
    slot = _fc_getSlot(i, size1, table1);
    fail_unless(slot == table1[i], "mismatch in table1");
    slot = _fc_getSlot(i, size2, table2);
    fail_unless(slot == table2[i], "mismatch in table2");
  }

  // Test deleting/getting new slots in middle of table
  rc = _fc_deleteSlot(num/2, size1, table1, &openSlots1);
  fail_unless(rc == FC_SUCCESS, "failed to delete slot");
  fail_unless(table1[num/2] == NULL, "deleted slot should be null");
  slot = _fc_getNewSlot(&size1, &table1, &openSlots1, sizeof(_FC_SlotHeader));
  uID++;
  fail_unless(slot == table1[num/2], "should have gotten slot in middle");
  fail_unless(slot->slotID == num/2, "should get slotID in middle");
  fail_unless(slot->uID == uID, "should have next uID");
  fail_unless(slot->name == NULL && slot->committed ==  0, 
	      "newSlot should have 0s");
  temp_num = 0;
  for (i = 0; i < size1; i++) 
    if (table1[i] != NULL && table1[i]->uID > 0)
      temp_num++;
  fail_unless(temp_num == num, "total number should be the same");

  // Test _fc_isHandleValid
  for (i = 0; i < num; i++) {
    int isValid;
    isValid = _fc_isHandleValid(i, table1[i]->uID, size1, table1);
    fail_unless(isValid, "should be valid");
    isValid = _fc_isHandleValid(i, uID+1, size1, table1);
    fail_unless(!isValid, "bad uID should not be valid");
  }

  // Test _FC_GET_HANDLE
  for (i = 0; i < num; i++) {
    object.header = *table1[i];
    _FC_GET_HANDLE(handle, &object);
    fail_unless(handle.slotID == table1[i]->slotID, "mismatch of slotID");
    fail_unless(handle.uID == table1[i]->uID, "mismatch of uID");
  }

  // cleanup
  for (i = 0; i < size1; i++)
    free(table1[i]);
  free(table1);
  fc_freeSortedIntArray(&openSlots1);
  for (i = 0; i < size2; i++)
    free(table2[i]);
  free(table2);
  fc_freeSortedIntArray(&openSlots2);
}
END_TEST

// *********************************************
// ***** Populate the Suite with the tests
// *********************************************

Suite *table_suite(void)
{
  Suite *suite = suite_create("Table");

  TCase *tc_table = tcase_create(" - Generic Table Interface ");

  // general table interface
  suite_add_tcase(suite, tc_table);
  tcase_add_test(tc_table, slot_header);
  tcase_add_test(tc_table, table);

  return suite;
}
