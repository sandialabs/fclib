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
 * \file fc.h
 * \brief Master header file for FCLib.
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/modules/fc.h,v $
 * $Revision: 1.59 $ 
 * $Date: 2006/08/30 19:20:01 $
 *
 * \description
 *
 *   Include this file to get access to all public FCLib functions.
 *
 * \modifications 
 *   - ?? Created 
 *   - July 10, 2002  W Koegler  Prettied up.
 *   - APR-14-2003  W Koegler  Added 'feature' and 'track'
 */

#ifndef _FC_H_
#define _FC_H_

// C library requirements
#include <stdlib.h>
#include <stdio.h>

// Data structures (order important)
#include "base.h"
#include "storage.h"
#include "library.h"
#include "dataset.h"
#include "sequence.h"
#include "mesh.h"
#include "subset.h"
#include "variable.h"
#include "feature.h"
#include "shape.h"


// File IO
#include "fileio.h"

// Characterizations & building blocks
#include "error.h"
#include "geom.h"
#include "sierra.h"
#include "statistics.h"
#include "threshold.h"
#include "elemdeath.h"
#include "series.h"
#include "topo.h"
#include "util.h"
#include "track.h"
#include "varmath.h"
#include "custom.h"

#endif // _FC_H_
