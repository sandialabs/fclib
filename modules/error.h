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
 * \file error.h
 * \brief Declarations for \ref ErrorHandling Module
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/modules/error.h,v $
 * $Revision: 1.25 $ 
 * $Date: 2006/08/30 19:20:01 $
 */

#ifndef _FC_ERROR_P_H_
#define _FC_ERROR_P_H_

#ifdef __cplusplus
extern "C" {
#endif

// NOTE! If error.c ever created, move the group comment below to error.c
/** 
 * \addtogroup ErrorHandling
 * \brief Data structures and routines to ease error handling.
 */

/** 
 * \ingroup ErrorHandling
 * \brief If the return code is an error, print diagnostic and exit.
 *
 * \description
 *
 *    This is a helper macro to ease writing command line programs; it
 *    groups the test for an error return code with the handling: exit with
 *    the same value as the return code. 
 *
 *    Example usage: fc_exitIfError(rc) = 
 *
 *    "Exiting with error FC_MEMORY_ERROR (-2) at [program.c:276]"
 *  
 *    Sometimes this macro has compiling problems if it's used in an if-else
 *    block. To fix, surround the call with { }.
 *
 * \modifications
 *    - 05/25/04 WSK Created
 *    - 11/28/2005 WSD. Added printing of type of error (used to be silent).
 */
#define fc_exitIfError(rc) \
  { \
    if ( (rc) != FC_SUCCESS ) { \
      fprintf(stderr, "Exiting with error %s (%d) at [%s:%d]\n", \
              fc_getReturnCodeText(rc), (rc), __FILE__, __LINE__);  \
      fflush(NULL); \
      exit(rc); \
    } \
  }

/*
 * \ingroup ErrorHandling
 * \brief If the return code is an error, print message and exit.
 *
 * \description
 * 
 *    This is a helper function to ease writing command line programs; it
 *    groups the test for an error return code with the handling: print message
 *    and exit.  The message can consist of a format string and arguments, just
 *    like printf (but you don't need the trailing '\\n').  The message will be
 *    'Error: Your message'.  It will exit with the same value as the return
 *    code.
 *
 *    Example usage: fc_exitIfError(rc, "Mesh #%d incomplete", i) = 
 *
 *    "Exiting with error FC_ERROR (-1) at [program.c:75]: Mesh #2 incomplete"
 *
 *    Sometimes this macro has compiling problems if it's used in an if-else
 *    block. To fix, surround the call with { }.
 *
 * \modifications
 *    - 2003-NOV-11  WSK  Created
 *    - 05/25/04 WSK changed name and changed to include the error condition 
 *        checking. Also deleted this function's unnecessary counterpart 'warn'
 *    - 11/28/2005 WSD. Changed to a macro and added printing of type of error.
 */
#define fc_exitIfErrorPrintf(rc, message, ...) \
  { \
    if ( (rc) != FC_SUCCESS ) { \
      fprintf(stderr, "Exiting with error %s (%d) at [%s:%d]: " message "\n", \
              fc_getReturnCodeText(rc), (rc), __FILE__, __LINE__, ##__VA_ARGS__); \
      fflush(NULL); \
      exit(rc); \
    } \
  }

#ifdef __cplusplus
}
#endif

#endif // _FC_ERROR_P_H_


