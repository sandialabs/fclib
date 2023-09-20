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
 * \file library.h
 * \brief Public declarations for the \ref Library Module.
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/modules/library.h,v $
 * $Revision: 1.13 $ 
 * $Date: 2006/08/30 19:20:01 $
 */

#ifndef _FC_LIBRARY_H_
#define _FC_LIBRARY_H_

#ifdef __cplusplus
extern "C" {
#endif

// Library verbosity
//-------------------------------
FC_ReturnCode fc_setLibraryVerbosity(FC_VerbosityLevel verbosity);
FC_VerbosityLevel fc_getLibraryVerbosity(void);

// Library init/final
//-------------------------------
FC_ReturnCode fc_initLibrary(void);
FC_ReturnCode fc_finalLibrary(void);

/** 
 * \name Message functions
 *
 * Use these to generate various levels of messages. The messages
 * that are actually sent depend on the verbosity level of the library.
 */
//-------------------------------
//@{
 
/** 
 * \ingroup Library
 * \brief Print an error message to stderr.
 *
 * \description
 *
 *    If the library verbosity indicates that error messages should
 *    be printed, this Macro prints the message to stderr.
 *    The message can consist of a format string and arguments,
 *    just like printf (but you don't need the trailing '\\n').
 *    The output will be 
 *
 *    'Error:function_name[filename:line#]: message'.
 *
 * \note
 *    This is a 'variadic macro' (i.e. it can have variable number of
 *    of arguments. I'm not sure how well supported it is
 *    over different compilers. __VA_ARGS__ is a predefined macro
 *    that maps the '...' in the argument list. The '##'
 *    is supposed to delete the preceding ',' if __VA_ARGS__ is empty.
 * 
 * \modifications 
 *    - 05/19/04  WSK, Created
 */
#define fc_printfErrorMessage(message, ...)  \
  { \
    if (fc_getLibraryVerbosity() >= FC_ERROR_MESSAGES)  { \
      fprintf(stderr, "FC Error:%s[%s:%d]: " message "\n", \
              __func__, __FILE__, __LINE__, ##__VA_ARGS__); \
      fflush(NULL); \
    } \
  }

/** 
 * \ingroup Library
 * \brief Print a warning message to stdout.
 *
 * \description
 *
 *    If the library verbosity indicates that warning messages should
 *    be printed, this Macro prints the message to stderr.
 *    The message can consist of a format string and arguments,
 *    just like printf (but you don't need the trailing '\\n').
 *    The output will be 
 *
 *    'FC Warning:function_name[filename:line#]: message'.
 *
 * \note
 *    This is a 'variadic macro' (i.e. it can have variable number of
 *    of arguments. I'm not sure how well supported it is
 *    over different compilers. __VA_ARGS__ is a predefined macro
 *    that maps the '...' in the argument list. The '##'
 *    is supposed to delete the preceding ',' if __VA_ARGS__ is empty.
 * 
 * \modifications 
 *    - 05/19/04  WSK, Created
 */
#define fc_printfWarningMessage(message, ...)  \
  { \
    if (fc_getLibraryVerbosity() >= FC_WARNING_MESSAGES) { \
      fprintf(stderr, "FC Warning:%s[%s:%d]: " message "\n", \
              __func__, __FILE__, __LINE__, ##__VA_ARGS__); \
      fflush(NULL); \
    } \
  }

/** 
 * \ingroup Library
 * \brief Print a log message to stdout.
 *
 * \description
 *
 *    If the library verbosity indicates that log messages should
 *    be printed, this Macro prints the message to stdout.
 *    The message can consist of a format string and arguments,
 *    just like printf (but you don't need the trailing '\\n').
 *    The output will be 
 *
 *    'FC Log:function_name[filename:line#]: message'.
 *
 * \note
 *    This is a 'variadic macro' (i.e. it can have variable number of
 *    of arguments. I'm not sure how well supported it is
 *    over different compilers. __VA_ARGS__ is a predefined macro
 *    that maps the '...' in the argument list. The '##'
 *    is supposed to delete the preceding ',' if __VA_ARGS__ is empty.
 * 
 * \modifications 
 *    - 05/19/04  WSK, Created
 */
#define fc_printfLogMessage(message, ...)  \
  { \
    if (fc_getLibraryVerbosity() >= FC_LOG_MESSAGES) { \
      fprintf(stdout, "FC Log:%s[%s:%d]: " message "\n", \
              __func__, __FILE__, __LINE__, ##__VA_ARGS__); \
      fflush(NULL); \
    } \
  }

/** 
 * \ingroup Library
 * \brief Print a debug message to stdout.
 *
 * \description
 *
 *    If the library verbosity indicates that debug messages should
 *    be printed, this Macro prints the message to stdout.
 *    The message can consist of a format string and arguments,
 *    just like printf (but you don't need the trailing '\\n').
 *    The output will be 
 *
 *    'FC Debug:function_name[filename:line#]: message'.
 *
 * \note
 *    This is a 'variadic macro' (i.e. it can have variable number of
 *    of arguments. I'm not sure how well supported it is
 *    over different compilers. __VA_ARGS__ is a predefined macro
 *    that maps the '...' in the argument list. The '##'
 *    is supposed to delete the preceding ',' if __VA_ARGS__ is empty.
 * 
 * \modifications 
 *    - 05/19/04  WSK, Created
 */
#define fc_printfDebugMessage(message, ...)  \
  { \
    if (fc_getLibraryVerbosity() >= FC_DEBUG_MESSAGES) { \
      fprintf(stdout, "FC Debug:%s[%s:%d]: " message "\n", \
              __func__, __FILE__, __LINE__, ##__VA_ARGS__); \
      fflush(NULL); \
    } \
  }

  // thinking about this
  /*
#define fc_returnErrorMessage(rc) \
  { \
    fc_printfErrorMessage("%s", fc_getReturnCodeText(rc)); \
    fflush(NULL); \
    return rc; \
 }
  */

//@}

#ifdef __cplusplus
}
#endif

#endif // _FC_LIBRARY_H_
