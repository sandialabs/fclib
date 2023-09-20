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
 * \file checkall.h
 * \brief Header for all checks.
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/unittest/checkall.h,v $
 * $Revision: 1.12 $ 
 * $Date: 2006/08/30 19:20:05 $
 *
 * \modifications
 *    - 12/20/04 WSK, created.
 */

#ifndef _CHECK_ALL_H_
#define _CHECK_ALL_H_

// global flags

int isForking;  // whether the tests are forked or not
FC_VerbosityLevel fc_messages; // library verbosity

// Declation of suites (the default case e.g. check)
Suite *base_suite(void);
Suite *util_suite(void);
Suite *storage_suite(void);
Suite *table_suite(void);
Suite *library_suite(void);
Suite *dataset_suite(void);
Suite *sequence_suite(void);
Suite *mesh_suite(void);
Suite *subset_suite(void);
Suite *variable_suite(void);
Suite *fileio_suite(void);
Suite *topo_suite(void);
Suite *geom_suite(void);
Suite *stats_suite(void);
Suite *thresh_suite(void);
Suite *series_suite(void);
Suite *elemdeath_suite(void);
Suite *shape_suite(void);
Suite *varmath_suite(void);
Suite *feature_suite(void);
Suite *track_suite(void);
Suite *custom_suite(void);

// Delcaration of suites for if using external data (e.g. check2)
Suite *fileio2_suite(void);

#endif // _CHECK_ALL_H_
