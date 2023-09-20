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
 * \file seriesP.h
 * \brief Private declarations for \ref Series module.
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/modules/seriesP.h,v $
 * $Revision: 1.5 $
 * $Date: 2006/08/30 19:20:01 $
 */

#ifndef _FC_SERIES_P_H_
#define _FC_SERIES_P_H_

#ifdef __cplusplus
extern "C" {
#endif

//??


//2/10/04 ACG - currently these next two are just used in checkStats.c, but
//        I anticipate these beign used in checkTimeseries.c as well. 
//        the appropriateness of this location should be revisted
//        as we see how the code evolves.. 
//4/13/05 ACG moved createRegularSequence back into timeseries since
//        will use it inconjunction with linearInterpolation. it may
//        move to sequence later. this also means note above now
//        refers just to the createSeqVariableFromVariable only now.
//6/1/05 ACG renamed Timeseries to Series

//this is used in series, checkseries, checkstats, and checkvarmath
FC_ReturnCode _fc_createSeqVariableFromVariable(
     int numStep, FC_Sequence seq, FC_Variable regular,
     char* seq_var_name, FC_Variable** seqvar);


//this is used both in series and checkseries
FC_ReturnCode _fc_makeNewSeqVar(int numStep, FC_Variable *seqvar,
				     char* varname, FC_Variable **newseqvar);




#ifdef __cplusplus
}
#endif

#endif // _FC_SERIES_H_
