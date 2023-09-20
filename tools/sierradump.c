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
 * \file sierradump.c
 * \brief Read a sierra input file and report it's content to standard output
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/tools/sierradump.c,v $
 * $Revision: 1.17 $ 
 * $Date: 2006/08/30 19:20:04 $
 *
 * \modifications
 *   - 2003-JUN-24  W Koegler Created to test new Sierra parser 
 */

// C library includes
#include <stdio.h>
#include <stdlib.h>

// fclib includes
#include "base.h"
#include "sierra.h"

int main(int argc, char** argv) 
{
  FC_ReturnCode rc; // return code
  char* filename;
  FC_SierraInfo* sierra;

  // Process arguments
  if (argc != 2) {
    printf("usage: sierradump <sierra_input_file>\n");
    fflush(NULL);
    exit(-1);
  }
  filename = argv[1];

  // do it
  rc = fc_parseSierraInputFile(filename, &sierra);
  if (rc != FC_SUCCESS) {
    printf("fc_parseSierraInputFile returned error: '%s'\n", 
           fc_getReturnCodeText(rc));
    fflush(NULL);
    exit(rc);
  }
  rc = fc_printSierraInfo(sierra);
  if (rc != FC_SUCCESS) {
    printf("fc_printSierraInfo returned error: '%s'\n", 
           fc_getReturnCodeText(rc));
    fflush(NULL);
    exit(rc);
  }

  // cleanup
  fc_freeSierraInfo(sierra);

  exit(0);
}




