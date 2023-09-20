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
 * \file fcP.h
 * \brief Master header file for private members and functions of FCLib.
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/modules/fcP.h,v $
 * $Revision: 1.51 $ 
 * $Date: 2006/08/30 19:20:01 $
 *
 * \description
 *
 *   Include this file to get access to all private FCLib functions (must
 *   include fc.h first).
 * 
 * \modifications 
 *   - ?? Created.
 *   - make July 10, 2002  W Koegler  Updated so it works.
 */

#ifndef _FC_P_H_
#define _FC_P_H_

// Data structures (?order important)
#include "storageP.h"
#include "tableP.h"
#include "libraryP.h"
#include "datasetP.h"
#include "sequenceP.h"
#include "meshP.h"
#include "subsetP.h"
#include "variableP.h"

// File IO
#include "fileioP.h"

// Characterizations & building blocks
#include "geomP.h"
#include "sierraP.h"
#include "topoP.h"
#include "varmathP.h"
#include "featureP.h"

#endif // _FC_P_H_
