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
 * \file base.c
 * \brief Implementation for \ref DataTypes, \ref DataInterface, and
 *     \ref Handles Modules.
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/modules/base.c,v $
 * $Revision: 1.74 $ 
 * $Date: 2006/10/19 03:14:50 $
 */

// C library includes
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// this module
#include "base.h"
#include "library.h"

/**
 * \addtogroup  DataTypes
 * \brief Data types of the library.
 *   
 * \description
 *
 *    The possible values of many types of data are represented with
 *    enumerations. This module contains the enum declarations and some basic
 *    functions for manipulating them.
 */

/**
 * \addtogroup  DataInterface
 * \brief  Primary interface between computational routines and 
 *         actual data. 
 *
 * \description
 *  
 *    The routines in this section are the primary interface between the
 *    computational routines in the Feature Characterization library and
 *    the actual data.  
 *
 *    There are five major data entities: datasets, sequences, meshes, subsets,
 *    and variables which follow this basic hierarchy (the library owns
 *    everything):
 *
 *     - Library
 *     - Dataset
 *     - Sequence
 *     - Mesh
 *     - Subset
 *     - Variable & Sequence Variable 
 *
 *    The data interface routines allow the calling code to get, query, copy,
 *    and delete each of these types of objects. In addition, you can
 *    load and write datasets. The objects are
 *    manipulated via opaque handles (see \ref Handles). These handles can be 
 *    assigned directly using '=', but cannot be checked for equality using 
 *    '==', so you must use \ref FC_HANDLE_EQUIV instead.
 *
 *    As a quick review, the library manages multiple datasets and each
 *    dataset can have multiple meshes. A sequence is the coordinates of 
 *    a parameter space orthogonal to the space of the mesh. The most common 
 *    type of sequence would be the time values of a time series. 
 *    A subset is a set of mesh subentities (e.g. vertices or elements).
 *    A variable is a function, such as temperature, over a mesh. It is
 *    associated with some subentity of a mesh such as vertices or
 *    elements. A sequence variable is actually an array of variable handles
 *    associated with a sequence, with one variable per step of the
 *    sequence.
 *
 *    It is very important to note that the library makes a distinction
 *    between meta data and "big data", and that access to these is treated
 *    very differently.
 *    Big data are the really large arrays of data that we want to 
 *    avoid duplicating or moving around. Currently, the coordinate arrays for
 *    the meshes and the sequences, and the data from the variables, are
 *    considered big data and everything else is meta data.
 *
 *    When users ask for metadata they get copies of the data and can
 *    manipulate it however they want, and they are responsible for freeing any
 *    arrays. On the other hand, when users ask for big data, they get a
 *    pointer to the original data (the names of such routines will end with
 *    'Ptr'). Users should treat this pointer as read only.  Users should never
 *    free big data directly, but should used the release and delete routines
 *    to try free up big data or to delete the
 *    entity that the big data is in.
 */

/**
 * \addtogroup  Handles
 * \brief Data objects handled through opaque handles.
 *   
 * \description
 *
 *    Handles are used to look up the appropriate data structure in
 *    the library's internals. Handles can be assigned 
 *    to (copies contents) with '=', but cannot be tested with '==',
 *    use \ref FC_HANDLE_EQUIV().
 *
 *    All of the handles have the same form
 *    (a struct containing the slot id and the sequence id), the library
 *    declares them as different types to clarify the user API and to add 
 *    more compile time type checking. But because they have the same
 *    form, general purpose routines and macros (like checking that a
 *    handle is valid) can be written. 
 *
 *    See \ref DataInterface for an explanation of using the data objects.
 */

/** 
 * \name Check that an enum has a valid value
 * 
 * These routines return true (1) if the int passed is a valid value for the
 * enumerated type, and false (0) if it is not. Enum values of the type 
 * FC_X_UNKNOWN are considered valid because they are possible values of the
 * enum. In many cases you may want to check for a "good" value and should
 * check "fc_isXValid(value) && value != FC_X_UNKNOWN".
 */
//-------------------------------
//@{

/**
 * \ingroup   DataTypes
 * \brief  Verify that the given value is valid for this enum
 *
 * \modifications   
 *    - 8/18/04 WSK, Created
 */
int fc_isVerbosityLevelValid(
  FC_VerbosityLevel verbosity   
) 
{
  switch(verbosity) {
  case FC_QUIET:             // fall through
  case FC_ERROR_MESSAGES:    // fall through
  case FC_WARNING_MESSAGES:  // fall through
  case FC_LOG_MESSAGES:      // fall through
  case FC_DEBUG_MESSAGES:    // fall through
  // no default case so will get warning if new values added to enum
    return 1;
  }
  return 0;
}

/**
 * \ingroup   DataTypes
 * \brief  Verify that the given value is valid for this enum
 *
 * \modifications   
 *    - 8/18/04 WSK, Created
 */
int fc_isReturnCodeValid(
  FC_ReturnCode rc      /**< input - enum to return printable string for */
) 
{
  switch(rc) {
  case FC_SUCCESS:        // fall through
  case FC_ERROR:          // fall through
  case FC_MEMORY_ERROR:   // fall through
  case FC_INPUT_ERROR:    // fall through
  case FC_FILE_IO_ERROR:  // fall through
  // no default case so will get warning if new values added to enum
    return 1;
  }
  return 0;
}

/**
 * \ingroup   DataTypes
 * \brief  Verify that the given value is valid for this enum
 *
 * \modifications   
 *    - 8/18/04 WSK, Created
 */
int fc_isElementTypeValid(
  FC_ElementType elemType   
) 
{
  switch(elemType) {
  case FC_ET_UNKNOWN:  // fall through
  case FC_ET_POINT:    // fall through
  case FC_ET_LINE:     // fall through
  case FC_ET_TRI:      // fall through
  case FC_ET_QUAD:     // fall through
  case FC_ET_TET:      // fall through
  case FC_ET_PYRAMID:  // fall through
  case FC_ET_PRISM:    // fall through
  case FC_ET_HEX:      // fall through
  case FC_ET_MIXED:    // fall through
  // no default case so will get warning if new values added to enum
    return 1;
  }
  return 0;
}

/**
 * \ingroup   DataTypes
 * \brief  Verify that the given value is valid for this enum
 *
 * \modifications   
 *    - 8/18/04 WSK, Created
 */
int fc_isDataTypeValid(
  FC_DataType dataType
) 
{
  switch (dataType) {
  case FC_DT_UNKNOWN:  // fall through
  case FC_DT_CHAR:     // fall through
  case FC_DT_INT:      // fall through
  case FC_DT_FLOAT:    // fall through
  case FC_DT_DOUBLE:   // fall through
  // no default case so will get warning if new values added to enum
    return 1;
  }
  return 0;
}

/**
 * \ingroup   DataTypes
 * \brief  Verify that the given value is valid for this enum
 *
 * \modifications   
 *    - 8/18/04 WSK, Created
 */
int fc_isMathTypeValid(
 FC_MathType mathtype
) 
{
  switch (mathtype) {
  case FC_MT_UNKNOWN:    // fall through
  case FC_MT_SCALAR:     // fall through
  case FC_MT_VECTOR:     // fall through
  case FC_MT_SYMTENSOR:  // fall through
  case FC_MT_TENSOR:     // fall through
  // no default case so will get warning if new values added to enum
    return 1;
  }
  return 0;
}

/**
 * \ingroup   DataTypes
 * \brief  Verify that the given value is valid for this enum
 *
 * \modifications   
 *    - 8/18/04 WSK, Created
 */
int fc_isAssociationTypeValid(
  FC_AssociationType assoc
) 
{
  switch (assoc) {
  case FC_AT_UNKNOWN:  // fall through
  case FC_AT_VERTEX:   // fall through
  case FC_AT_EDGE:     // fall through
  case FC_AT_FACE:     // fall through
  case FC_AT_ELEMENT:  // fall through
  case FC_AT_WHOLE_MESH:    // fall through
  case FC_AT_WHOLE_DATASET: // fall through
  // no default case so will get warning if new values added to enum
    return 1;
  }
  return 0;
}

//@}

/** 
 * \name Get the name of an enum value as text
 *
 * These routines return a text string holding the name of the enum value.
 * If the enum value is invalid, they return NULL. The returned string 
 * should NOT be freed.
 */
//-------------------------------
//@{

/**
 * \ingroup   DataTypes
 * \brief  Return the name of the \ref FC_VerbosityLevel's value.
 *
 * \description
 *
 *    Returns a string representing the name of the enum value or
 *    NULL is the enum value is invalid. The string should NOT be freed.
 * 
 * \modifications   
 *    - 3/30/04 WSK, Created
 */
char *fc_getVerbosityLevelText(
  FC_VerbosityLevel verbosity   
) 
{
  switch(verbosity) {
  case FC_QUIET:             return "FC_QUIET";
  case FC_ERROR_MESSAGES:    return "FC_ERROR_MESSAGES";
  case FC_WARNING_MESSAGES:  return "FC_WARNING_MESSAGES";
  case FC_LOG_MESSAGES:      return "FC_LOG_MESSAGES";
  case FC_DEBUG_MESSAGES:    return "FC_DEBUG_MESSAGES";
  // no default case so will get warning if new values added to enum
  }
  return NULL;
}

/**
 * \ingroup  DataTypes
 * \brief  Return the name of the \ref FC_ReturnCode's value.
 *
 * \description
 *  
 *    Returns a string representing the name of the enum value or
 *    NULL is the enum value is invalid. The string should NOT be freed.
 * 
 * \modifications 
 *    - 2003-SEP-22  W Koegler  Created
 */
char *fc_getReturnCodeText(
  FC_ReturnCode rc      /**< input - enum to return printable string for */
) 
{
    switch(rc) {
      case FC_SUCCESS:          return "FC_SUCCESS";
      case FC_ERROR:            return "FC_ERROR";
      case FC_MEMORY_ERROR:     return "FC_MEMORY_ERROR";
      case FC_INPUT_ERROR:      return "FC_INPUT_ERROR";
      case FC_FILE_IO_ERROR:    return "FC_FILE_IO_ERROR";
      // no default case so will get warning if new values added to enum
    }

    return NULL;
}

/**
 * \ingroup  DataTypes
 * \brief  Return the name of the \ref FC_ElementType's value.
 *
 * \description
 *
 *    Returns a string representing the name of the enum value or
 *    NULL is the enum value is invalid. The string should NOT be freed.
 *
 * \modifications   
 *    - Nancy Collins  Created.
 *    - 2003-NOV-17  WSK  Changed to return the enum values' name instead
 *      of more human readable word.
 */
char *fc_getElementTypeText(
  FC_ElementType elemType   
) 
{
  switch(elemType) {
  case FC_ET_UNKNOWN:    return "FC_ET_UNKNOWN";
  case FC_ET_POINT:      return "FC_ET_POINT";
  case FC_ET_LINE:       return "FC_ET_LINE";
  case FC_ET_TRI:        return "FC_ET_TRI";
  case FC_ET_QUAD:       return "FC_ET_QUAD";
  case FC_ET_TET:        return "FC_ET_TET";
  case FC_ET_PYRAMID:    return "FC_ET_PYRAMID";
  case FC_ET_PRISM:      return "FC_ET_PRISM";
  case FC_ET_HEX:        return "FC_ET_HEX";
  case FC_ET_MIXED:      return "FC_ET_MIXED";
  // no default case so will get warning if new values added to enum
  }
  return NULL;
}

/**
 * \ingroup  DataTypes
 * \brief  Return the name of the \ref FC_DataType's value.
 *
 * \description
 *
 *    Returns a string representing the name of the enum value or
 *    NULL is the enum value is invalid. The string should NOT be freed.
 * 
 * \modifications 
 *    - Nancy Collins Created.
 *    - 2003-NOV-17  WSK  Changed to return the enum values' name instead
 *      of more human readable word.
 */
char *fc_getDataTypeText(
 FC_DataType dataType
) 
{
  switch (dataType) {
  case FC_DT_UNKNOWN:       return "FC_DT_UNKNOWN";
  case FC_DT_CHAR:          return "FC_DT_CHAR";
  case FC_DT_INT:           return "FC_DT_INT";
  case FC_DT_FLOAT:         return "FC_DT_FLOAT";
  case FC_DT_DOUBLE:        return "FC_DT_DOUBLE";
  // no default case so will get warning if new values added to enum
  }
  return NULL;
}

/**
 * \ingroup  DataTypes
 * \brief  Return the name of the FC_MathType's value.
 *
 * \description
 *
 *    Returns a string representing the name of the enum value or
 *    NULL is the enum value is invalid. The string should NOT be freed.
 * 
 * \modifications 
 *    - Nancy Collins Created.
 *    - 2003-NOV-17  WSK  Changed to return the enum values' name instead
 *      of more human readable word.
 */
char *fc_getMathTypeText(
 FC_MathType mathtype
) 
{
  switch (mathtype) {
  case FC_MT_UNKNOWN:       return "FC_MT_UNKNOWN";
  case FC_MT_SCALAR:        return "FC_MT_SCALAR";
  case FC_MT_VECTOR:        return "FC_MT_VECTOR";
  case FC_MT_SYMTENSOR:     return "FC_MT_SYMTENSOR";
  case FC_MT_TENSOR:        return "FC_MT_TENSOR";
  // no default case so will get warning if new values added to enum
  }
  return NULL;
}

/**
 * \ingroup  DataTypes
 * \brief  Return the name of the FC_AssociationType's value.
 * 
 * \description
 *
 *    Returns a string representing the name of the enum value or
 *    NULL is the enum value is invalid. The string should NOT be freed.
 * 
 * \modifications 
 *    -Nancy Collins Created.
 *    - 2003-NOV-17  WSK  Changed to return the enum values' name instead
 *      of more human readable word.
 */
char *fc_getAssociationTypeText(
  FC_AssociationType assoc
) 
{
  switch (assoc) {
  case FC_AT_UNKNOWN:        return "FC_AT_UNKNOWN";
  case FC_AT_VERTEX:         return "FC_AT_VERTEX";
  case FC_AT_EDGE:           return "FC_AT_EDGE";
  case FC_AT_FACE:           return "FC_AT_FACE";
  case FC_AT_ELEMENT:        return "FC_AT_ELEMENT";
  case FC_AT_WHOLE_MESH:     return "FC_AT_WHOLE_MESH"; 
  case FC_AT_WHOLE_DATASET:  return "FC_AT_WHOLE_DATASET"; 
  // no default case so will get warning if new values added to enum
  }
  return NULL;
}

//@}

/** \name Get information about an FC_ElementType */
//-------------------------------
//@{

/**
 * \ingroup  DataTypes
 * \brief  Returns number of topological dimensions for a type of element.
 *
 * \description
 *
 *    Returns the topological dimensionality of the type of element.
 *    If the type is unknown (FC_ET_UNKNOWN), -1 is returned.
 *
 * \modifications  
 *    - Nancy Collins Created.
 *    - 2003-NOV-17 WSK Changed to return topodim instead of
 *      passing back out in argument list.
 */
int fc_getElementTypeTopoDim(
  FC_ElementType elemType
) 
{ 
  switch (elemType) {
  case FC_ET_MIXED:    // fall through to FC_ET_UNKNOWN
  case FC_ET_UNKNOWN:  return -1;
  case FC_ET_POINT:    return 0;
  case FC_ET_LINE:     return 1;
  case FC_ET_TRI:      // fall through to FC_ET_QUAD
  case FC_ET_QUAD:     return 2;
  case FC_ET_TET:      // fall through to FC_ET_HEX
  case FC_ET_PYRAMID:  // fall through to FC_ET_HEX
  case FC_ET_PRISM:    // fall through to FC_ET_HEX
  case FC_ET_HEX:      return 3;
  // no default case so will get warning if new values added to enum
  }
  fc_printfWarningMessage(
     "Warning: invalid value for enum FC_ElementType"); 
  return -1;
}

/**
 * \ingroup  DataTypes
 * \brief  Returns number of vertices in a type of element.
 *
 * \description
 *
 *    Returns the number of vertices of the type of element.
 *    If the type is unknown (FC_ET_UNKNOWN), -1 is returned.
 *
 * \modifications  
 *    - Nancy Collins  Created.
 *    - 2003-NOV-17 WSK Changed to return numVertex instead of
 *      passing back out in argument list.
 */
int fc_getElementTypeNumVertex(
  FC_ElementType elemType
) 
{ 
  switch (elemType) {
  case FC_ET_MIXED:   // fall through to FC_ET_UNKNOWN
  case FC_ET_UNKNOWN:  return -1;
  case FC_ET_POINT:    return 1;   
  case FC_ET_LINE:     return 2;
  case FC_ET_TRI:      return 3; 
  case FC_ET_QUAD:     // fall through to FC_ET_TET  
  case FC_ET_TET:      return 4;
  case FC_ET_PYRAMID:  return 5;
  case FC_ET_PRISM:    return 6;
  case FC_ET_HEX:      return 8;
  // no default case so will get warning if new values added to enum
  }
  fc_printfWarningMessage(
     "Warning: invalid value for enum FC_ElementType"); 
  return -1;
}

/**
 * \ingroup  DataTypes
 * \brief  Returns number of edges in a type of element.
 *
 * \description
 *
 *    Returns the number of edges in a type of element.
 *    If the type is unknown (FC_ET_UNKNOWN), -1 is returned.
 *
 * \modifications  
 *    - Nancy Collins  Created.
 *    - 2003-NOV-17 WSK Changed to return numVertex instead of
 *      passing back out in argument list.
 */
int fc_getElementTypeNumEdge(
  FC_ElementType elemType
)
{ 
  switch (elemType) {
  case FC_ET_MIXED:     // fall through to FC_ET_UNKNOWN
  case FC_ET_UNKNOWN:   return -1;
  case FC_ET_POINT:     return  0;
  case FC_ET_LINE:      return  1;  
  case FC_ET_TRI:       return  3; 
  case FC_ET_QUAD:      return  4;   
  case FC_ET_TET:       return  6;
  case FC_ET_PYRAMID:   return  8; 
  case FC_ET_PRISM:     return  9;  
  case FC_ET_HEX:       return 12; 
  // no default case so will get warning if new values added to enum
  }
  fc_printfWarningMessage(
     "Warning: invalid value for enum FC_ElementType"); 
  return -1;
}

/**
 * \ingroup  DataTypes
 * \brief  Returns number of faces in a type of element.
 * 
 * \description
 *
 *    Returns the number of faces for a type of element.
 *    If the type is unknown (FC_ET_UNKNOWN), -1 is returned.
 *
 * \modifications  
 *    - Nancy Collins Created.
 *    - 2003-NOV-17 WSK Changed to return numFace instead of
 *      passing back out in argument list.
 */
int fc_getElementTypeNumFace(
  FC_ElementType elemType
) 
{ 
  switch (elemType) {
  case FC_ET_MIXED:    // fall through to FC_ET_UNKNOWN
  case FC_ET_UNKNOWN:  return -1;
  case FC_ET_POINT:    // fall through to FC_ET_LINE
  case FC_ET_LINE:     return 0;  
  case FC_ET_TRI:      // fall through to FC_ET_QUAD
  case FC_ET_QUAD:     return 1;
  case FC_ET_TET:      return 4;
  case FC_ET_PYRAMID:  // fall through to FC_ET_PRISM
  case FC_ET_PRISM:    return 5;
  case FC_ET_HEX:      return 6;
  // no default case so will get warning if new values added to enum
  }
  fc_printfWarningMessage(
     "Warning: invalid value for enum FC_ElementType"); 
  return -1;
}

/**
 * \ingroup  DataTypes
 * \brief  Returns number of the entity type in a type of element.
 * 
 * \description
 *
 *    Returns the number of the entity type for a type of element.  If the type
 *    is unknown (FC_ET_UNKNOWN), -1 is returned.  If the entity type is
 *    inappropriate (e.g. FC_AT_WHOLE_MESH on any type of element, or e.g. 
 *    FC_AT_FACE on for FC_ET_LINE) it returns -1. Note, passing
 *    FC_AT_ELEMENT will return 1.
 *
 *    For entity specific behavior see fc_getElementTypeNumVertex(),
 *    fc_getElementTypeNumEdge(), and fc_getElementTypeNumFace().
 *
 * \modifications  
 *    - 10/24/2005 WSD Created.
 */
int fc_getElementTypeNumEntity(
  FC_ElementType elemType,  /**< input - Type of element */
  FC_AssociationType assoc  /**< Type of entity */
) 
{ 
  switch (assoc) {
  case FC_AT_VERTEX:         return fc_getElementTypeNumVertex(elemType);
  case FC_AT_EDGE:           return fc_getElementTypeNumEdge(elemType);
  case FC_AT_FACE:           return fc_getElementTypeNumFace(elemType);
  case FC_AT_ELEMENT:        return 1;
  case FC_AT_UNKNOWN:        // fall through
  case FC_AT_WHOLE_MESH:     return -1;
  case FC_AT_WHOLE_DATASET:  return -1;
  // no default case so will get warning if new values added to enum
  }
  fc_printfWarningMessage(
     "Warning: invalid value for enum FC_AssociationType"); 
  return -1;
}

/**
 * \ingroup  DataTypes
 * \brief  Returns the element type of the faces of the element.
 * 
 * \description
 *
 *    For FC_ET_PYRAMID and FC_ET_PRISM the face type is a 
 *    mix of FC_ET_TRI and FC_ET_QUAD, so the type
 *    FC_ET_MIXED is return.
 *
 *    For FC_ET_POINT and FC_ET_LINE which don't have faces,
 *    the type returned is FC_ET_UNKNOWN.
 *
 * \modifications  
 *    - 02/05/04 WSK, Created.
 */
FC_ElementType fc_getElementTypeFaceType(
  FC_ElementType elemType
) 
{ 
  switch (elemType) {
  case FC_ET_MIXED:    // fall through to FC_ET_UNKNOWN
  case FC_ET_UNKNOWN:  return FC_ET_UNKNOWN;
  case FC_ET_POINT:    // fall through to FC_ET_LINE
  case FC_ET_LINE:     return FC_ET_UNKNOWN;  
  case FC_ET_TRI:      return FC_ET_TRI;
  case FC_ET_QUAD:     return FC_ET_QUAD;
  case FC_ET_TET:      return FC_ET_TRI;
  case FC_ET_PYRAMID:  // fall through to FC_ET_PRISM
  case FC_ET_PRISM:    return FC_ET_MIXED;
  case FC_ET_HEX:      return FC_ET_QUAD;
  // no default case so will get warning if new values added to enum
  }
  fc_printfWarningMessage(
     "Warning: invalid value for enum FC_ElementType"); 
  return FC_ET_UNKNOWN;
}

//@}

/** \name Get information about an FC_DataType */
//-------------------------------
//@{

/**
 * \ingroup  DataTypes
 * \brief  Returns the number of bytes in an FC_DataType, like sizeof()
 *
 * \description
 *
 *    Returns the size (number of bytes) of a data type.
 *    If the type is unknown (FC_DT_UNKNOWN), -1 is returned.
 *
 * \modifications  
 *    - Nancy Collins  Created.
 *    - 2003-NOV-17 WSK Changed to return the number of bytes instead
 *      of passing back out in argument list. Also asks system for size
 *      instead of hard coding.
 */
int fc_sizeofDataType(
  FC_DataType dataType
) 
{
  switch (dataType) {
  case FC_DT_UNKNOWN: return -1;
  case FC_DT_CHAR:    return sizeof(char);
  case FC_DT_INT:     return sizeof(int);
  case FC_DT_FLOAT:   return sizeof(float);
  case FC_DT_DOUBLE:  return sizeof(double);
  // no default case so will get warning if new values added to enum
  }
  fc_printfWarningMessage(
     "Warning: invalid value for enum FC_DataType"); 
  return -1;
}

//@}

/** 
 * \ingroup  Handles
 * \brief Null (default) dataset handle
 */
FC_Dataset FC_NULL_DATASET = { -1, -1 };
/** 
 * \ingroup  Handles
 * \brief Null (default) sequence handle
 */
FC_Sequence FC_NULL_SEQUENCE = { -1, -1 };
/** 
 * \ingroup  Handles
 * \brief Null (default) mesh handle
 */
FC_Mesh FC_NULL_MESH = { -1, -1 };
/** 
 * \ingroup  Handles
 * \brief Null (default) variable handle
 */
FC_Variable FC_NULL_VARIABLE = { -1, -1 };
/** 
 * \ingroup  Handles
 * \brief Null (default) subset handle
 */
FC_Subset FC_NULL_SUBSET = { -1, -1 };
