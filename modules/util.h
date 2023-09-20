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
 * \file util.h
 * \brief Public declarations for \ref Utilities module.
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/modules/util.h,v $
 * $Revision: 1.29 $ 
 * $Date: 2006/08/30 19:20:01 $
 */

#ifndef _FC_UTIL_H_
#define _FC_UTIL_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "math.h" // for FC_VALUE_EQUIV macro
#include "float.h" // for FC_FLOAT_EQUIV and FC_DOUBLE_EQUIV

// For the Equiv tests below, the a==b might seem to be unnecessary
// but it catchs the a=b=0 case so we never have divide by zero.

/** \name Macros for floating point equality comparison. */
//-------------------------------
//@{

/**
 * \ingroup  PublicFloatingPoint
 * \brief Test equality of floating point numbers
 *
 * \description 
 *
 *   This macro returns true if both values are equal (that is the
 *   scaled difference between them is less than epsilon) or if they
 *   are both essentially zero (absolute value less than min). This 
 *   macro should only be used with floats and doubles (and long doubles).
 *
 *   Good guess values for EPS and MIN come from float.h 
 *   See discussion in 'Subtraction and Comparing Floats'
 *   and 'Precision of Mixed Type Operations'
 *
 * \modifications
 *   2003-JUL-23  W Koegler  Created.
 *
 */
#define FC_VALUE_EQUIV(a,   /* floating point number */ \
                       b,   /* floating point number */ \
                       EPS, /* smallest difference between numbers, scaled */ \
                       MIN  /* smallest number distinct from zero */ ) \
( ( (a) == (b) ) ||                                                    \
  ( (fabs(a) < (MIN)) && (fabs(b) < (MIN)) ) ||                        \
  ( fabs( ((a)-(b)) / (fabs(a) > fabs(b) ? (a) : (b)) ) < (EPS) )      \
)

/**
 * \ingroup  PublicFloatingPoint
 * \brief Test equality of two values with float precision
 *
 * \description 
 *
 *   This is a version of FC_VALUE_EQUIV specific for comparing values
 *   which have the precision of a float.
 *   This macro returns true if both values are equal (that is the
 *   scaled difference between them is less than FLT_EPSILON) or if they
 *   are both essentially zero (absolute value less than FLT_MIN). 
 *
 *   See discussion in 'Subtraction and Comparing Floats'
 *   and 'Precision of Mixed Type Operations'
 *
 * \modifications
 *   2003-JUL-23  W Koegler  Created.
 *
 */
#define FC_FLT_EQUIV(a,   /* floating point number */ \
                     b    /* floating point number */ ) \
( ( (a) == (b) ) ||                                             \
  ( (fabs(a) < FLT_MIN) && (fabs(b) < FLT_MIN) ) ||             \
  ( fabs( ((a)-(b)) / (fabs(a) > fabs(b) ? (a) : (b)) ) < FLT_EPSILON ) \
)

/**
 * \ingroup  PublicFloatingPoint
 * \brief Test equality of two values with double precision
 *
 * \description
 * 
 *   This is a version of FC_VALUE_EQUIV specific for comparing values
 *   which have the precision of a double.
 *   This macro returns true if both values are equal (that is the
 *   scaled difference between them is less than DBL_EPSILON) or if they
 *   are both essentially zero (absolute value less than DBL _MIN). 
 *
 *   See discussion in 'Subtraction and Comparing Floats'
 *   and 'Precision of Mixed Type Operations'
 *
 * \modifications
 *   2003-JUL-23  W Koegler  Created.
 */
#define FC_DBL_EQUIV(a,   /* floating point number */ \
                     b    /* floating point number */ ) \
( ( (a) == (b) ) ||                                             \
  ( (fabs(a) < DBL_MIN) && (fabs(b) < DBL_MIN) ) ||             \
  (fabs( ((a)-(b)) / (fabs(a) > fabs(b) ? (a) : (b)) ) < DBL_EPSILON )  \
)

//@}

// floating point comparisons
int fc_eq(double x, double y, double eps,double min);
int fc_neq(double x, double y, double eps,double min);
int fc_lt(double x, double y, double eps,double min);
int fc_lteq(double x, double y, double eps,double min);
int fc_gt(double x, double y, double eps,double min);
int fc_gteq(double x, double y, double eps,double min);
int fc_eqf(double x, double y);
int fc_neqf(double x, double y);
int fc_ltf(double x, double y);
int fc_lteqf(double x, double y);
int fc_gtf(double x, double y);
int fc_gteqf(double x, double y);
int fc_eqd(double x, double y);
int fc_neqd(double x, double y);
int fc_ltd(double x, double y);
int fc_lteqd(double x, double y);
int fc_gtd(double x, double y);
int fc_gteqd(double x, double y);

// path processing
char* fc_getDirname(char* path);
char* fc_getBasename(char* path);
char* fc_getBasenameWOExtension(char* path, int ext_flag);
char* fc_getExtension(char* path, int ext_flag);

/** 
 * \ingroup MiscUtilities
 * \brief Macro that gives value of PI.
*/
#define FC_PI 3.14159265358979323846

// misc
int fc_intCompare(const void* n, const void* m);
FC_ReturnCode fc_sortIntArray(int num, int* array);
FC_ReturnCode fc_addValueToIntArray(int **array, int* num, int new_value);
char* fc_replaceChars(char* in_string, char toBeReplaced, char replaceWith);

#ifdef __cplusplus
}
#endif

#endif // _FC_UTIL_H_
