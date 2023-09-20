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
 * \file library.c
 * \brief Implementation for \ref Library module.
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/modules/library.c,v $
 * $Revision: 1.19 $ 
 * $Date: 2006/08/30 19:20:01 $
 */

// C library includes
#include <stdio.h>
#include <stdlib.h>

// fc library dependencies
#include "base.h"
#include "storage.h"
#include "dataset.h"
#include "datasetP.h"
#include "sequenceP.h"
#include "meshP.h"
#include "subsetP.h"
#include "variableP.h"

// this module
#include "library.h"
#include "libraryP.h"

/**
 * \addtogroup  Library
 * \brief  Operations on the library.
 *
 * \description
 *
 *    (For an explanation of general data manipulations, see 
 *    \ref DataInterface.)
 */

/**
 * \ingroup Library
 * \defgroup PrivateLibrary (Private)
 */ 

/** 
 * \ingroup PrivateLibrary
 * \brief Global variable for initialization state of the library.
 *
 * \description
 *
 *   1 = initialized, 0 = not initialized.
 */
static int library_init = 0;

/**
 * \ingroup PrivateLibrary
 * \brief  Global variable for Verbosity level of the library.    
 */
static FC_VerbosityLevel library_verbosity = FC_QUIET;

/** \name Set and query the library verbosity level. */
//-------------------------------
//@{

/** 
 * \ingroup Library
 * \brief Set the library verbosity level.
 *
 * \description
 * 
 *    See \ref FC_VerbosityLevel for explanation of possible settings.
 *
 * \modifications
 *    - 3/12/04 RM created.
 */
FC_ReturnCode fc_setLibraryVerbosity( 
 FC_VerbosityLevel verbosity   /**< input - desire verbosity */
)
{
  int flag;

  // test args
  flag = fc_isVerbosityLevelValid(verbosity);
  if (flag < 1) {
    fc_printfWarningMessage("Invalid verbosity level (%d)", verbosity);
    return FC_INPUT_ERROR;
  }

  // do it
  fc_printfLogMessage("Setting library verbosity to %s",
		      fc_getVerbosityLevelText(verbosity)); 
  library_verbosity = verbosity;

  return FC_SUCCESS;
}

/** 
 * \ingroup Library
 * \brief Get the current verbosity level of the library.
 *
 * \description
 * 
 *    See \ref FC_VerbosityLevel for explanation of possible settings.
 *
 * \modifications
 *   -  3/12/04 RM created.
 */
FC_VerbosityLevel fc_getLibraryVerbosity(void) 
{
  return  library_verbosity;
}

//@}
 
/** \name Initialize/finalize the library. */
//-------------------------------
//@{

/**
 * \ingroup  Library
 * \brief  Setup Feature Characterization Library for use.
 *
 * \description  
 *
 *    Executes library initialization.  Must be called before using any
 *    functions except for fc_getLibraryVerbosity and fc_setLibraryVerbosity
 *    (\ref fc_setLibraryVerbosity should be called before \ref
 *    fc_initLibrary() if you want any to insure that anything setup by the
 *    library (e.g.  the SAF i/o libraries) are initialized with your requested
 *    verbosity level).
 *
 *    Calling fc_initLibrary() more than once is an error.
 *
 * \modifications  
 *    - Nancy Collins, Created.
 *    - 06/03/04 WSK, modified to insure that library cannot be reinited.
 */
FC_ReturnCode fc_initLibrary(void) 
{
  FC_ReturnCode rc;
  static int been_here = 0;
  
  // reiniting library is an error
  if (been_here) {
    if (library_init) {
      fc_printfErrorMessage("Library already initialized");
    }
    else 
      fc_printfErrorMessage("Library cannot be reinitialized once finalized");
    return FC_ERROR;
  }

  // init the library
  fc_printfLogMessage("Initializing FC Library");

  // init IO
  rc = _fc_initFileIO();
  if (rc != FC_SUCCESS) {
    fc_printfErrorMessage("Failed to initialize File IO");
    return rc;
  }
  
  library_init = 1;
  been_here = 1;

  return FC_SUCCESS;
}

/**
 * \ingroup  Library
 * \brief   Cleanup when done using Feature Characterization Library.
 *
 * \description
 *  
 *    Executes library finalization.  Should be called after
 *    last access to library. Will delete all datasets.
 *
 *    It is not an error to final an uninitialized library, but you
 *    can not reinitialize a finalized library.
 *
 * \modifications  
 *    - Nancy Collins, Created.
 *    - 12/08/03 WSK Made sure it deletes databases and frees
 *            all of the tables.
 */
FC_ReturnCode fc_finalLibrary( 
 void
) {
  int i;
  FC_ReturnCode rc, rc_keep;
  _FC_DsSlot* dsSlot;
  FC_Dataset dataset;
  char* temp_name;

  // finalizing uninitialized library does not generate an error
  if (library_init == 0) {
    fc_printfWarningMessage("Library already finalized");
    return FC_SUCCESS;
  }
  
  // log message
  fc_printfLogMessage("Finalizing FC Library");

  // delete all datasets
  rc_keep = FC_SUCCESS;
  for (i = 0; i < _fc_getDsTableSize(); i++) {
    dsSlot = _fc_getDsSlotFromID(i);
    if (dsSlot) {
      _FC_GET_HANDLE(dataset, dsSlot);
      fc_getDatasetName(dataset, &temp_name);
      rc = fc_deleteDataset(dataset);
      if (rc != FC_SUCCESS) {
        rc_keep = rc;
        fc_printfWarningMessage("Failed to delete dataset '%s'", temp_name);
      }
      free(temp_name);
    }
  }

  // free tables
  _fc_freeDsTable();
  _fc_freeSeqTable();
  _fc_freeMeshTable();
  _fc_freeSubTable();
  _fc_freeVarTable();

  // finalize file IO
  rc = _fc_finalFileIO();
  if (rc != FC_SUCCESS) 
    fc_printfWarningMessage("Failed to finalize FileIO");
  
  // library is considered finalized even if error generated
  library_init = 0;

  // if any errors, make sure you return one of them
  if (rc == FC_SUCCESS)
    return rc_keep;
  else
    return rc;
}

//@}

/**
 * \ingroup  PrivateLibrary
 * \brief  Check initialization state of library.
 *
 * \description
 *  
 *    Returns true (1) is library is available for use.
 *
 * \modifications  
 *   - 05/20/04  W Koegler  Created.
 */
int _fc_isLibraryInit(void) 
{
  return library_init;
}
