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
 * \file sierra.h
 * \brief Public declarations for the \ref Sierra module.
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/modules/sierra.h,v $
 * $Revision: 1.23 $ 
 * $Date: 2006/08/30 19:20:01 $
 */

#ifndef _FC_SIERRA_H_
#define _FC_SIERRA_H_

#ifdef __cplusplus
extern "C" {
#endif

/**
 * \ingroup Sierra
 * \brief Container for Sierra function information
 * 
 * \description
 *
 *   This struct holds information specific to Sierra functions.
 *   The abscissa (x) and ordinate (y) labels are optional and default
 *   values are "" (empty string). The function values are stored
 *   in a single array as 'xyxyxyxy...'. Minimum and maximum values
 *   of the abscissa and ordinate are also provided. Note that
 *   a function can have a single value, in which case the values
 *   of the min and max ordinate are undefined.
 *
 * \modifications
 *   - 2003-JUN-24  W Koegler  Created.
 *   - 2003-SEP-03  W Koegler  Expanded to support all function variables.
 */
typedef struct {
  // From the Sierra input deck
  char name[100];   /**< name of the function */
  char type[100];   /**< type is either 'constant' or 'piecewise linear' */
  char abscissa_label[100]; /**< optional label for abscissa (default = "") */
  char ordinate_label[100]; /**< optional label for abscissa (default = "") */
  int numValue;  /**< function has either 1 or an even number of values */
  double* values; /**< values are stored in array ordered as 'xyxyxy...' */
  // Calculated values
  double abscissa_min;  /**< minimum value of the abscissa (x min) */
  double ordinate_min;  /**< minimum value of the ordinate (y min) */
  double abscissa_max;  /**< maximum value of the abscissa (x max) */
  double ordinate_max;  /**< maximum value of the ordinate (y max) */
} FC_SierraFuncInfo;

/**
 * \ingroup  Sierra
 * \brief Container for Sierra spot weld information
 * 
 * \description
 *
 *   This struct holds information specific to Sierra spot welds.
 *   Indices into the owning FC_SierraInfo struct's function list
 *   are also provided for the normal and tangential functions.
 *   (E.g. given FC_SierraInfo sierra, the normal function of the
 *   ith spotweld can be found as 
 *   sierra->functions[sierra->spotwelds[i]->normalId].)
 *
 * \modifications
 *   - 2003-JUN-24  W Koegler  Created.
 *   - 2003-SEP-02  W Koegler  Added all possible spot weld variables.
 *   - 2007-OCT-15 ACG removed refs to saf.
 */
typedef struct {
  // From the Sierra input deck
  char name[100];    /**< name of the spot weld */
  char nodeset[100]; /**< name of node set (welded to the surface) */
  char surface[100]; /**< name of surface (welded to the node set) */
  char normal[100];  /**< name of normal displacement function */
  char tangent[100]; /**< name of tangential displacement function */
  double normal_scale; /**< scale factor for normal displacement function */ 
  double tangent_scale; /**< scale factor for tangential displacement function */ 
  double exponent;    /**< value of the failure envelope exponent */
  int decay_cycles;  /**< number of failure decay cycles */
  // Computed after parsing to make access easier
  int normalId;  /**< index of the normal function in SierraInfo function list */
  FC_SierraFuncInfo* normal_p; /**< pointer to normal function */
  int tangentId; /**< index of tangential function in SierraInfo function list */
  FC_SierraFuncInfo* tangent_p; /**< pointer to normal function */
  char nodeset_name[100]; /**< nodeset name */ 
  char surface_name[100]; /**< surface name */ 
} FC_SierraSpotWeldInfo;

/**
 * \ingroup  Sierra
 * \brief Container for Sierra information
 * 
 * \description
 *
 *   This struct holds information found in a Sierra input file.
 *   Sierra is a framework that provides management of a variety of
 *   physics simulation codes and their data. The input file contains
 *   instructions to perform a simulation.
 *
 *   Currently, the SierraInfo struct holds only the information needed 
 *   to analyze Presto spot welds. This includes the name of the
 *   the Sierra problem, alias for the result variables "displacement"
 *   and "force_external", the number of functions and pointers to
 *   FC_SierraFunctionInfo structs, and the number of spots welds
 *   and pointers to FC_SpotWeldInfo structs.
 *
 *   Basic handling routines (init, clear, print) are provided only
 *   for the SierraInfo struct, and not its subelements as they
 *   (the subelements) are not expected to be handled independently
 *   of the SierraInfo struct. 
 *
 *   Knowledge of this struct and it's member structs is necessary to
 *   access the contained data.
 *
 * \modifications
 *   - 2003-JUN-24  W Koegler Created.
 *
 * \todo
 *   Instead of a name for each results variable, keep a list
 *   of variable/alias pairs. Also, keep track of type of variable
 *   i.e. node, element, or global.
 */
typedef struct {
  char name[100];  /**< name of the sierra problem */
  char displ_alias[100];  /**< alias of the displacement variable */
  char f_ext_alias[100];  /**< alias of the force_external variable */
  int numFunction;        /**< number of functions */
  FC_SierraFuncInfo** functions; /**< array of pointers to functions */
  int numSpotWeld;        /**< number of spot welds */
  FC_SierraSpotWeldInfo** spotwelds; /**< array of pointers to spot welds */
} FC_SierraInfo;

// Basic SierraInfo handling routines
// This not a general parse, only enough info for the following characterizations
FC_ReturnCode fc_parseSierraInputFile(char* filename, FC_SierraInfo** sierra);
void fc_freeSierraInfo(FC_SierraInfo* sierra);
FC_ReturnCode fc_printSierraInfo(FC_SierraInfo* sierra);

// Presto characterizations

// based on the spot weld model of presto
FC_ReturnCode fc_evalPrestoSpotWeld(FC_SierraSpotWeldInfo *spotweld,
                    FC_Subset node_set, FC_Subset side_set,
                    int numStep, FC_Variable* node_displs, 
                    FC_Variable* side_displs,
                    double** params, int* break_id);

#ifdef __cplusplus
}
#endif

#endif  // _FC_SIERRA_H_
