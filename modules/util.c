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
 * \file util.c
 * \brief Implementation for \ref Utilities module.
 * 
 * $Source: /usr/local/Repositories/fcdmf/fclib/modules/util.c,v $
 * $Revision: 1.56 $ 
 * $Date: 2006/08/30 19:20:01 $
 */

// C library includes
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h>

// fc library dependencies
#include "base.h"
#include "library.h"

// this module
#include "util.h"

/**
 * \addtogroup Utilities
 * \brief Support routines.
 *
 * \description
 *
 *    This module containes various utilities for general computing purposes.
 *    So far these include a few functions to process path strings to 
 *    get base file names, extension, etc. (see below) and also macros
 *    and routines for scientifically dealing with floating point numbers
 *    (see \ref PublicFloatingPoint).
 */

/**
 * \ingroup Utilities
 * \defgroup PublicFloatingPoint Floating Point Comparisons
 * \brief Some useful constants and floats comparisons routines.
 *
 * \description 
 *
 *   It is important to be aware that errors can creep
 *   into calculations using non-integers (floating point numbers) because
 *   these numbers cannot necessarily be stored exactly in the computer,
 *   and because of roundoff errors and ignorance of the number
 *   of significant digits. Generally computations are done with more
 *   precision than is needed to minimize roundoff errors. However, errors
 *   from some numeric operations can still cause havoc. One particular
 *   operation is subtraction (see Subtraction and Comparing Floats in
 *   this chapter).
 *
 *   Some useful constants are defined in the C library file, float.h:
 *   
 *   *_EPSILON (where * = FLT, DBL, or LDBL) 
 *   is defined as the difference between 1.0 and the smallest number
 *   greater than 1.0 than can be represented as that type. This gives
 *   a measure of the granularity of the each type. FLT_EPSILON is
 *   on the order of 10^-7 and DBL_EPSILON is on the order of 10^-16.
 *
 *   *_DIG (where * = FLT, DBL, or LDBL) is the maximum possible number of 
 *   (decimal) significant digits that can be stored in that type.
 *
 *   *_MAX (where * = FLT, DBL, or LDBL) is the maximum number representable
 *   by that type (order of 10^38 for FLT and 10^308 for DBL). If you
 *   try to increase this number, you get Inf (infinity).
 *
 *   *_MIN (where * = FLT, DBL, or LDBL) is the minimum normalized number
 *   for that type (order of 10^-38 for FLT and 10^-308 for DBL). 
 *   This is not the minimum number representable (that number is
 *   *_MIN times *_EPSILON and is on the order of 10^-45 for FLT and
 *   10^-324 for DBL). Instead, this is the minimum number that can still
 *   have the full number of significant digits for that type.
 *
 *   (Note that because ints can be stored exactly, there are no corresponding
 *   constants in float.h.)
 *
 * \modifications
 *   2002-JUN-24  W Koegler  Touched up.
 *
 * -----------------------------------------------------------------------------
 *
 * Subtraction and Comparing Floats 
 *
 *
 * \modifications
 *   2002-JUN-24  W Koegler  Touched up.
 *
 * \description 
 *
 *   Subtraction is an arithmetic operation that
 *   can lead to the loss in the number of significant digits. For example, 
 *   2.00 - 1.99 (both have 3 significant digits) is equal to 0.01 (one
 *   significant digit), not 0.010. These means that the results of a
 *   subtraction operation can be highly sensitive to roundoff errors.
 *   For instance, if 2.00 is stored as 1.995 in the computer, and 
 *   1.99 is stored as 1.994, then the difference is 0.001.  Or if 2.00
 *   is stored as 2.004 and 1.99 is stored as 1.985, then the difference
 *   is 0.019. The answer with roundoff error can range from twice as
 *   large as the real answer to 10 times too small. 
 *
 *   This problem (loss of significant digits)
 *   becomes more pronounced as the two numbers
 *   being operated on approach equality. Therefore it is often
 *   advisable to test such an operation before deciding what to do
 *   with the results. For the operation 
 *
 *   a - b = c ,
 *
 *   if
 *
 *   |(c/a)| < epsilon 
 *
 *   (where |x| is the absolute value of x), and epsilon is some very small 
 *   number,
 *   then a and b are approaching equality and we may want to do
 *   something special with the results. Often this means assuming
 *   a = b, and c = 0.
 *
 *   One practical application is testing whether a and b are equal
 *
 *   if |(a-b)/a| < epsilon, then a = b.
 *
 *   CAVEAT-Unfortunately, if either a or b are zero, we have problems.
 *   If a = 0, we have a divide by zero and our test will always
 *   fail because x/0 = Inf > epsilon. If b = 0, then a/a = 1 and again
 *   our test will always fail. Sometimes we will be comparing
 *   a really really small number to 0, and expect them to be
 *   equal, but the above test will always fail!!!
 *
 *   If one value is zero, if they are equal the other one will
 *   be close to zero. Close to zero would be defined as
 *   |b| < almost_zero. Choice of almost_zero will be problem dependent,
 *   but a good guess would be FLT_MIN and DBL_MIN from float.h.
 *
 *   So now our procedure is to first test if a or b is equal to zero.
 *   If that is true, our test becomes
 *
 *   if |a-b| < almost_zero, then a = b = 0
 *
 *   and if neither a or b are equal to zero, then we can use our original test
 *
 *   if |(a-b)/a| < epsilon, then a = b.
 *
 *   Another practical application is to set the result of subtracting
 *   two similar numbers to zero. This will help prevent propagation of
 *   roundoff error in following calculations.
 *
 *   If |c/a| < epsilon, then c = 0.
 *
 *   Choosing epsilon will depend on the particular circumstances. It
 *   will depend on how strict you want to be, how many significant
 *   digits you should expect the values to have, how much roundoff
 *   error may have crept in, etc.
 *
 *   Some helpful constants can be found in the C library file, float.h
 *   (see the introduction to this chapter). Below are are some 
 *   specific examples.
 *
 *   A constant times *_EPSILON may be appropriate. Another possibility
 *   is to use 
 *
 *   epsilon = 10^(-1 * FLT_DIG)   (or DBL_DIG).
 *
 *   Typically *_EPSILON (where * = FLT, DBL or LDBL) will be too stringent 
 *   as round off errors can accumulate in previous calculations. 
 *   For example, a set of ci's, ci = ai + bi, i = 0 .. 10, where calculated
 *   with ai and bi provided as randomly generated floats. 
 *   The largest value of |(c - a - b)/c|
 *   seen was 4.7e-8. However, more operations are done such as
 *   |((c^2 - a^2 - 2*A*B - b^2)/c)/c|, the largest value seen was 1.3e-7.
 * 
 * --------------------------------------------------------------------------
 *
 *  Precision of Mixed Type Operations 
 *
 *
 * \modifications
 *   2002-JUN-24  W Koegler  Created. 
 *
 * \description
 *
 *   Here a '?' is symbol for any operation (+, -, *, etc.)
 * 
 *   INT ? FLT -> FLT
 *
 *   INT ? DBL -> DBL
 *
 *   FLT ? DBL -> FLT
 *   
 *   However, precision is not explicitly carried with floating point numbers.
 *   So when that number is passed from the current context, it may suddenly
 *   acquire more significance. (E.g. is the product of a float and a
 *   double is stored in a double, that double really has the precision of
 *   a float, but who can tell? Code should check the types and store
 *   the result in a float.)   
 *
 *   A full description of error propagation analysis is beyond the
 *   scope of this documentation.
 *
 *   \todo Add more documentation? 
 *
 * --------------------------------------------------------------------------
 *
 * I/O of Floating Point Numbers 
 *
 *
 * \modifications
 *   2002-JUN-24  W Koegler  Created.
 *
 * \description 
 *   
 *   If you want to keep the precision of a double, you should write
 *   all of it's digits.
 *
 *   When writing human readable I/O where you care about the number
 *   of significant digits, use %g and specify the number of
 *   significant digits as '.x' where x is the number of significant
 *   digits. Be aware that %g will drop trailing zeros. %g also
 *   automatically switches between decimal or scientific notation
 *   depending on a variety of factors. If you care about significant
 *   digits, never use %f, but %e might be o.k.
 *
 *   Example: "%11.6g" is used in the library to print floats.
 *   The '.6' part says, 'use 6 significant digits' and the '11'
 *   says 'try to fit it in 11 characters' (if it doesn't fit it just
 *   takes more room) (use '-11' to get left justified). Lets say we have
 *
 *   float a = 1.1, b = 222.222, and c = 0.000987643;
 *
 *   then 'printf("%11.6g %11.6g %11.6g)' will print (spacing may not 
 *   be correct)
 *
 *
 *   12345678911 12345678911 12345678911 (count spaces 1 to 11)
 *
 *           1.1     222.222 9.87643E-04
 *  
 *   and 'printf("%11.6f %11.6f %11.6f)' will print
 *   
 *      1.10000x  222.222xxx    0.000987         
 *
 *   where the 'x' are digits outside of the significant digit range
 *   and are probably garbage. Note that "%11.15g" is used in the
 *   library to print doubles. This will try to fit 15 significant
 *   digits into 11 characters, which is fine if there are lots
 *   of trailing zeros. If there are not trailing zeros, it will
 *   take at least 15 characters of space as needed.
 *
 *   For portable regression testing, you should not write out more than
 *   the number of significant digits because there is likely to
 *   be differences in the non significant digits. However, differences
 *   in the significant digits may also show up if there is more round
 *   off error on another machine.
 */

/**
 * \ingroup Utilities
 * \defgroup FilePathUtilities File Path Processing
 * \brief Various routines for decomposing file paths.
 */

/**
 * \ingroup Utilities
 * \defgroup MiscUtilities Misc Utilities
 */

/** \name General purpose floating point comparisons. */
//-------------------------------
//@{

/**
 * \ingroup PublicFloatingPoint Floating Point Comparisons
 * \brief Test if two doubles are equal.
 *
 * \description
 *
 *    This performs a comparison that takes into account floating point
 *    error. See FC_VALUE_EQUIV and \ref PublicFloatingPoint for an explanation
 *    of what eps and min mean.
 *
 * \modifications
 *    - 01/04/05 WSK moved from threshold here.
 */ 
int fc_eq(
  double x,    /**< first value */
  double y,    /**< second value */
  double eps,  /**< smallest distinquishable difference between values */
  double min   /**< smallest value distinquishable from zero */
) {
  return FC_VALUE_EQUIV(x, y, eps, min);
}

/**
 * \ingroup PublicFloatingPoint Floating Point Comparisons
 * \brief Test if two doubles are not equal.
 *
 * \description
 *
 *    This performs a comparison that takes into account floating point
 *    error. See FC_VALUE_EQUIV and \ref PublicFloatingPoint for an explanation
 *    of what eps and min mean.
 *
 * \modifications
 *    - 01/04/05 WSK moved from threshold here.
 */ 
int fc_neq(
  double x,    /**< first value */
  double y,    /**< second value */
  double eps,  /**< smallest distinquishable difference between values */
  double min   /**< smallest value distinquishable from zero */
) {
  return !FC_VALUE_EQUIV(x, y, eps, min);
}

/**
 * \ingroup PublicFloatingPoint Floating Point Comparisons
 * \brief Test if the first double is less than the second double.
 *
 * \description
 *
 *    This performs a comparison that takes into account floating point
 *    error. See FC_VALUE_EQUIV and \ref PublicFloatingPoint for an explanation
 *    of what eps and min mean.
 *
 * \modifications
 *    - 01/04/05 WSK moved from threshold here.
 */ 
int fc_lt(
  double x,    /**< first value */
  double y,    /**< second value */
  double eps,  /**< smallest distinquishable difference between values */
  double min   /**< smallest value distinquishable from zero */
) {
  if (x < y && !FC_VALUE_EQUIV(x, y, eps, min))
    return 1;
  else
    return 0;
}

/**
 * \ingroup PublicFloatingPoint Floating Point Comparisons
 * \brief Test if the first double is less than or equal to the second double.
 *
 * \description
 *
 *    This performs a comparison that takes into account floating point
 *    error. See FC_VALUE_EQUIV and \ref PublicFloatingPoint for an explanation
 *    of what eps and min mean.
 *
 * \modifications
 *    - 01/04/05 WSK moved from threshold here.
 */ 
int fc_lteq(
  double x,    /**< first value */
  double y,    /**< second value */
  double eps,  /**< smallest distinquishable difference between values */
  double min   /**< smallest value distinquishable from zero */
) {
  if (x < y || FC_VALUE_EQUIV(x, y, eps, min))
    return 1;
  else
    return 0;
}

/**
 * \ingroup PublicFloatingPoint Floating Point Comparisons
 * \brief Test if the first double is greater than the second double.
 *
 * \description
 *
 *    This performs a comparison that takes into account floating point
 *    error. See FC_VALUE_EQUIV and \ref PublicFloatingPoint for an explanation
 *    of what eps and min mean.
 *
 * \modifications
 *    - 01/04/05 WSK moved from threshold here.
 */ 
int fc_gt(
  double x,    /**< first value */
  double y,    /**< second value */
  double eps,  /**< smallest distinquishable difference between values */
  double min   /**< smallest value distinquishable from zero */
) {
  if (x > y && !FC_VALUE_EQUIV(x, y, eps, min))
    return 1;
  else
    return 0;
}

/**
 * \ingroup PublicFloatingPoint Floating Point Comparisons
 * \brief Test if the first double is greater than or equal to the second value.
 *
 * \description
 *
 *    This performs a comparison that takes into account floating point
 *    error. See FC_VALUE_EQUIV and \ref PublicFloatingPoint for an explanation
 *    of what eps and min mean.
 *
 * \modifications
 *    - 01/04/05 WSK moved from threshold here.
 */ 
int fc_gteq(
  double x,    /**< first value */
  double y,    /**< second value */
  double eps,  /**< smallest distinquishable difference between values */
  double min   /**< smallest value distinquishable from zero */
) {
  if (x > y || FC_VALUE_EQUIV(x, y, eps, min))
    return 1;
  else
    return 0;
}

//@}

/** \name Float versions of floating point comparisons. */
//-------------------------------
//@{

/**
 * \ingroup PublicFloatingPoint Floating Point Comparisons
 * \brief Test if two doubles are equal, assuming float precision.
 *
 * \description
 *
 *    This performs a comparison that takes into account floating point
 *    error. See \ref PublicFloatingPoint for more information.
 *
 * \modifications
 *    - 06/21/06 WSD Created.
 */ 
int fc_eqf(
  double x,    /**< first value */
  double y     /**< second value */
) {
  return FC_VALUE_EQUIV(x, y, FLT_EPSILON, FLT_MIN);
}

/**
 * \ingroup PublicFloatingPoint Floating Point Comparisons
 * \brief Test if two doubles are not equal, assuming float precision.
 *
 * \description
 *
 *    This performs a comparison that takes into account floating point
 *    error. See \ref PublicFloatingPoint for more information.
 *
 * \modifications
 *    - 06/21/06 WSD Created.
 */ 
int fc_neqf(
  double x,    /**< first value */
  double y    /**< second value */
) {
  return !FC_VALUE_EQUIV(x, y, FLT_EPSILON, FLT_MIN);
}

/**
 * \ingroup PublicFloatingPoint Floating Point Comparisons
 * \brief Test if the first double is less than the second double,
 *    assuming float precision.
 *
 * \description
 *
 *    This performs a comparison that takes into account floating point
 *    error. See \ref PublicFloatingPoint for more information.
 *
 * \modifications
 *    - 06/21/06 WSD Created.
 */ 
int fc_ltf(
  double x,    /**< first value */
  double y     /**< second value */
) {
  if (x < y && !FC_VALUE_EQUIV(x, y, FLT_EPSILON, FLT_MIN))
    return 1;
  else
    return 0;
}

/**
 * \ingroup PublicFloatingPoint Floating Point Comparisons
 * \brief Test if the first double is less than or equal to the second double,
 *    assuming float precision.
 *
 * \description
 *
 *    This performs a comparison that takes into account floating point
 *    error. See \ref PublicFloatingPoint for more information.
 *
 * \modifications
 *    - 06/21/06 WSD Created.
 */ 
int fc_lteqf(
  double x,    /**< first value */
  double y     /**< second value */
) {
  if (x < y || FC_VALUE_EQUIV(x, y, FLT_EPSILON, FLT_MIN))
    return 1;
  else
    return 0;
}

/**
 * \ingroup PublicFloatingPoint Floating Point Comparisons
 * \brief Test if the first double is greater than the second double,
 *    assuming float precision.
 *
 * \description
 *
 *    This performs a comparison that takes into account floating point
 *    error. See \ref PublicFloatingPoint for more information.
 *
 * \modifications
 *    - 06/21/06 WSD Created.
 */ 
int fc_gtf(
  double x,    /**< first value */
  double y     /**< second value */
) {
  if (x > y && !FC_VALUE_EQUIV(x, y, FLT_EPSILON, FLT_MIN))
    return 1;
  else
    return 0;
}

/**
 * \ingroup PublicFloatingPoint Floating Point Comparisons
 * \brief Test if the first double is greater than or equal to the second
 *     value, assmuming float precision. 
 *
 * \description
 *
 *    This performs a comparison that takes into account floating point
 *    error. See \ref PublicFloatingPoint for more information.
 *
 * \modifications
 *    - 06/21/06 WSD Created.
 */ 
int fc_gteqf(
  double x,    /**< first value */
  double y     /**< second value */
) {
  if (x > y || FC_VALUE_EQUIV(x, y, FLT_EPSILON, FLT_MIN))
    return 1;
  else
    return 0;
}

//@}

/** \name Float versions of floating point comparisons. */
//-------------------------------
//@{

/**
 * \ingroup PublicFloatingPoint Floating Point Comparisons
 * \brief Test if two doubles are equal, assuming double precision.
 *
 * \description
 *
 *    This performs a comparison that takes into account floating point
 *    error. See \ref PublicFloatingPoint for more information.
 *
 * \modifications
 *    - 06/21/06 WSD Created.
 */ 
int fc_eqd(
  double x,    /**< first value */
  double y     /**< second value */
) {
  return FC_VALUE_EQUIV(x, y, DBL_EPSILON, DBL_MIN);
}

/**
 * \ingroup PublicFloatingPoint Floating Point Comparisons
 * \brief Test if two doubles are not equal, assuming double precision.
 *
 * \description
 *
 *    This performs a comparison that takes into account floating point
 *    error. See \ref PublicFloatingPoint for more information.
 *
 * \modifications
 *    - 06/21/06 WSD Created.
 */ 
int fc_neqd(
  double x,    /**< first value */
  double y    /**< second value */
) {
  return !FC_VALUE_EQUIV(x, y, DBL_EPSILON, DBL_MIN);
}

/**
 * \ingroup PublicFloatingPoint Floating Point Comparisons
 * \brief Test if the first double is less than the second double,
 *    assuming double precision.
 *
 * \description
 *
 *    This performs a comparison that takes into account floating point
 *    error. See \ref PublicFloatingPoint for more information.
 *
 * \modifications
 *    - 06/21/06 WSD Created.
 */ 
int fc_ltd(
  double x,    /**< first value */
  double y     /**< second value */
) {
  if (x < y && !FC_VALUE_EQUIV(x, y, DBL_EPSILON, DBL_MIN))
    return 1;
  else
    return 0;
}

/**
 * \ingroup PublicFloatingPoint Floating Point Comparisons
 * \brief Test if the first double is less than or equal to the second double,
 *    assuming double precision.
 *
 * \description
 *
 *    This performs a comparison that takes into account floating point
 *    error. See \ref PublicFloatingPoint for more information.
 *
 * \modifications
 *    - 06/21/06 WSD Created.
 */ 
int fc_lteqd(
  double x,    /**< first value */
  double y     /**< second value */
) {
  if (x < y || FC_VALUE_EQUIV(x, y, DBL_EPSILON, DBL_MIN))
    return 1;
  else
    return 0;
}

/**
 * \ingroup PublicFloatingPoint Floating Point Comparisons
 * \brief Test if the first double is greater than the second double,
 *    assuming double precision.
 *
 * \description
 *
 *    This performs a comparison that takes into account floating point
 *    error. See \ref PublicFloatingPoint for more information.
 *
 * \modifications
 *    - 06/21/06 WSD Created.
 */ 
int fc_gtd(
  double x,    /**< first value */
  double y     /**< second value */
) {
  if (x > y && !FC_VALUE_EQUIV(x, y, DBL_EPSILON, DBL_MIN))
    return 1;
  else
    return 0;
}

/**
 * \ingroup PublicFloatingPoint Floating Point Comparisons
 * \brief Test if the first double is greater than or equal to the second
 *     value, assmuming float precision. 
 *
 * \description
 *
 *    This performs a comparison that takes into account floating point
 *    error. See \ref PublicFloatingPoint for more information.
 *
 * \modifications
 *    - 06/21/06 WSD Created.
 */ 
int fc_gteqd(
  double x,    /**< first value */
  double y     /**< second value */
) {
  if (x > y || FC_VALUE_EQUIV(x, y, DBL_EPSILON, DBL_MIN))
    return 1;
  else
    return 0;
}

//@}

/**
 * \ingroup FilePathUtilities
 * \brief Get the dirname in a path string.
 *
 * \description
 *
 *    Returns the part of the path before the last slash (ignoring trailing
 *    slashes). For example "/usr/local/bin" -> "/usr/local" or
 *    "../../cherry.txt" -> "../.." or "grumpy.xslt" -> ".".
 *
 *    User is responsible for freeing the returned string.
 *  
 *    If no path is provided, this routine returns NULL.
 *
 * \modifications
 *    - 7/14/2005 WSD Created.
 */ 
char* fc_getDirname(char* path) {
  char* temp_path;
  char* temp_str;
  char* dir;

  // early return
  if (path == NULL)
    return NULL;
  
  // Use C library function to get dir
  temp_path = (char*)malloc((strlen(path)+1)*sizeof(char));
  strcpy(temp_path, path);
  dir = dirname(temp_path);

  // copy, cleanup & return
  temp_str = (char*)malloc((strlen(dir)+1)*sizeof(char));
  strcpy(temp_str, dir);
  free(temp_path);
  return temp_str;
}

/**
 * \ingroup FilePathUtilities
 * \brief Get the basename in a path string.
 *
 * \description
 *
 *    Returns the part of the path after the last slash (ignoring trailing
 *    slashes). For example "/usr/local/bin" -> "bin" or
 *    "../../cherry.txt" -> "cherry.txt" or "grumpy.xslt" -> "grumpy.xslt".
 *
 *    User is responsible for freeing the returned string.
 *
 *    If no path is provided, this routine returns NULL.
 *
 * \modifications
 *    - 7/14/2005 WSD Created.
 */ 
char* fc_getBasename(char* path) {
  char* temp_path;
  char* temp_str;
  char* base;

  // early return
  if (path == NULL)
    return NULL;
  
  // Use C library function to get base 
  temp_path = (char*)malloc((strlen(path)+1)*sizeof(char));
  strcpy(temp_path, path);
  base = basename(temp_path);

  // copy, cleanup & return
  temp_str = (char*)malloc((strlen(base)+1)*sizeof(char));
  strcpy(temp_str, base);
  free(temp_path);
  return temp_str;
}

/**
 * \ingroup FilePathUtilities
 * \brief Get the root of a basename in a path string.
 *
 * \description
 *
 *    Returns the part of the path after the last slash (ignoring trailing
 *    slashes) minus any extension. For example "/usr/local/bin" -> "bin" or
 *    "../../cherry.txt" -> "cherry" or "grumpy.xslt" -> "grumpy".
 *
 *    The user can specify whether the extension starts with the
 *    first "." (0) or the last "." (1). For example, "bundle.tar.gz"
 *    would be "bundle" or "bundle.tar". A leading "." is not treated
 *    as an extension ".". For example. ".cshrc" -> ".cshrc";
 *
 *    User is responsible for freeing the returned string.
 *
 *    If no path is provided, this routine returns NULL.
 *
 * \modifications
 *    - 7/14/2005 WSD Created.
 */ 
char* fc_getBasenameWOExtension(
  char* path,   /**< the path */
  int ext_flag  /**< 0 = extension starts at 1st ".", 1 = extentsion starts at
		  last "." */
) {
  size_t i;
  char* temp_str;

  // early return
  if (path == NULL)
    return NULL;
  
  // get base
  temp_str = fc_getBasename(path);

  // clip off extentions
  if (ext_flag == 0) {
    for (i = 1; i < strlen(temp_str); i++) {
      if (temp_str[i] == '.') {
	temp_str[i] = '\0';
	break;
      }
    }
  }
  else {
    for (i = strlen(temp_str); i > 0; i--) {
      if (temp_str[i] == '.') {
	temp_str[i] = '\0';
	break;
      }
    }
  }

  return temp_str;
}

/**
 * \ingroup FilePathUtilities
 * \brief Get the extension of path
 *
 * \description
 *
 *    Returns the extension of a file if it exists. For example
 *    "/usr/local/bin" -> NULL or "../../cherry.txt" -> "txt" or
 *    "grumpy.xslt" -> "xslt".
 *
 *    The user can specify whether the extension starts with the
 *    first "." (0) or the last "." (1). For example, "bundle.tar.gz"
 *    would be "tar.gz" or ".gz". A leading "." is not treated
 *    as an extension ".". For example. ".cshrc" -> "";
 *
 *    NOTE: If there is no extension, a string is still returned, it
 *    just happens to be empty (and you still have to free it).
 *
 *    User is responsible for freeing the returned string.
 *
 * \modifications
 *    - 7/14/2005 WSD Created.
 */ 
char* fc_getExtension(
  char* path,   /**< input - the path */
  int ext_flag  /**< input - 0 = extension starts at 1st ".", 1 = extentsion 
		   starts at last "." */
) {
  size_t i;
  char* temp_str;
  char* ext;
  int start = -1;

  // early return
  if (path == NULL)
    return NULL;

  // get base
  temp_str = fc_getBasename(path);
  
  // find start of extension (".")
  if (ext_flag == 0) {
    for (i = 1; i < strlen(temp_str); i++) {
      if (temp_str[i] == '.') {
	start = i;
	break;
      }
    }
  }
  else {
    for (i = strlen(temp_str); i > 0; i--) {
      if (temp_str[i] == '.') {
	start = i;
	break;
      }
    }
  }

  // make the extension
  if (start < 0) {
    ext = (char*)malloc(sizeof(char));
    ext[0] = '\0';
  }
  else {
    ext = (char*)malloc((strlen(temp_str)-start)*sizeof(char));
    strcpy(ext, &temp_str[start+1]);
  }
  free(temp_str);

  return ext;
}

/**
 * \ingroup MiscUtilities
 * \brief Compare two ints.
 *
 * \description
 *
 *    Compares two ints and returns 1 if the first is greater,
 *    0 if they are equal, and -1 if the 2nd is greater.
 *    Return values are what routines like qsort expect.
 *
 * \modifications
 *    - 10/10/2005 WSD. Promoted (used to be inside fc_sortIntArray()).
 */
int fc_intCompare(const void* n, const void* m) 
{
  int a = ((const int*)n)[0];
  int b = ((const int*)m)[0];
  if (a > b)
    return 1;
  else if (a < b)
    return -1;
  else
    return 0;
}

/**
 * \ingroup MiscUtilities
 * \brief Sort an int array
 *
 * \description
 *
 *    Just a wrapper around qsort and fc_intCompare.
 *
 * \modifications
 *    - 07/01/04 WSK Created.
 */
FC_ReturnCode fc_sortIntArray(
  int num, /**< input - number of entries in array */
  int* array /**< input - the array */
) 
{
  if (num < 1 || array == NULL) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_INPUT_ERROR));
    return FC_INPUT_ERROR;
  }  

  qsort((void*)array, (size_t)num, sizeof(int), fc_intCompare);
  
  return FC_SUCCESS;
}

/**
 * \ingroup MiscUtilities
 * \brief  Add a new value to an array of ints.
 *
 * \description
 *  
 *    Given an array, the number of elements in the array,
 *    and a new value, append the value to the array.
 *    The array is assumed to be dynamically generated and the
 *    user is responsibly for freeing it when it is no longer
 *    in use. The new value is appended to the end and the count is
 *    increased by 1.
 *    A return value not equal to 0 indicates a failed malloc
 *    and the array and the arguments should remain unchanged.
 *
 * \modifications  
 *    - FEB-13-2002  W Koegler  Created.
 */
FC_ReturnCode fc_addValueToIntArray(int **array_p, int* num_p, int new_value) {
  int *temp_array;

  // check input
  if (!array_p || !num_p || *num_p < 0)
    return FC_ERROR;

  // copy values & replace array
  temp_array = malloc(sizeof(int)*(*num_p + 1));
  if (!temp_array) {
    fc_printfErrorMessage("%s", fc_getReturnCodeText(FC_MEMORY_ERROR));
    return FC_MEMORY_ERROR;
  }
  if (*num_p > 0) {
    memcpy(temp_array, *array_p, sizeof(int)*(*num_p));
    free(*array_p);
  }
  *array_p = temp_array;

  // add new value
  (*array_p)[*num_p] = new_value;
  (*num_p)++;

  return FC_SUCCESS;
}

/**
 * \ingroup MiscUtilities
 * \brief Replace occurances of a char in a string.
 *
 * \description
 *    
 *    Return a new string (to be freed by the caller) where occurences
 *    of a specified char are replaced with another specified char.
 *
 *    For example, fc_replaceChar("hi there!", ' ', '_') returns "hi_there!");
 */
char* fc_replaceChars(
  char* in_string,     /**< input - input string */
  char toBeReplaced,   /**< input - char to be replaced */
  char replaceWith     /**< input - char to be inserted */
) {
  size_t i, len;
  char* temp_string;

  // check input
  if (!in_string) 
    return NULL;

  // do it
  len = strlen(in_string);
  temp_string = (char*)malloc((len+1)*sizeof(char));
  if (!temp_string)
    return NULL;
  strcpy(temp_string, in_string);
  for (i = 0; i < len; i++)
    if (temp_string[i] == toBeReplaced)
      temp_string[i] = replaceWith;

  // done
  return temp_string;
}

