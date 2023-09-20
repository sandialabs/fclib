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
 * \file errortesthelper.c
 * \brief Test exit functions in \ref ErrorHandling Module
 *
 * \description
 *    
 *    Different cases try different exit/error possibilities.
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/regtest/errortesthelper.c,v $
 * $Revision: 1.4 $ 
 * $Date: 2006/09/08 23:39:29 $
 *
 * \modifications
 *    - 12/8/04 W Koegler. Created
 */

#include <string.h>
#include <fc.h>

int main(int argc, char** argv)
{
  FC_ReturnCode rc;

  // 
  // initialize the fcdmf library
  //
  rc = fc_setLibraryVerbosity(FC_WARNING_MESSAGES);
  if (rc != FC_SUCCESS)
    fc_exitIfError(rc);
  rc = fc_initLibrary();
  if (rc != FC_SUCCESS)
    fc_exitIfError(rc);

  //
  // Do the case indicated by args
  //
  if (argc != 2) {
    usage: 
      printf("ERROR: Expected a single arg with an acceptable case number\n");
      exit(-1);
  }
  if (!strcmp(argv[1], "case1")) {
    fc_exitIfError(FC_MEMORY_ERROR);
  }
  else if (!strcmp(argv[1], "case2")) {
    fc_exitIfErrorPrintf(FC_FILE_IO_ERROR, "Counting %d, %d & %d", 12, 44, 16);
  }
  else
    goto usage;

  // 
  // finalize the fcdmf library
  //
  fc_finalLibrary();
  
  exit(FC_SUCCESS);
}
