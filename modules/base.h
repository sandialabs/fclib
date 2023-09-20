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
 * \file base.h
 * \brief Public declarations for \ref DataTypes, \ref DataInterface, and
 *     \ref Handles Modules.
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/modules/base.h,v $
 * $Revision: 1.49 $ 
 * $Date: 2006/10/19 03:14:50 $
 */

#ifndef _FC_BASE_H_
#define _FC_BASE_H_

#ifdef __cplusplus
extern "C" {
#endif

/** \name  DataType Enums */
//-------------------------------
//@{

/**
 * \ingroup DataTypes
 * \brief  Possible verbosity levels
 *
 * \description
 *
 *    Generally, error messages are for conditions that halt execution 
 *    (e.g. "Error: cannot operate on NULL input"); warning messages 
 *    are for unexpected conditions that can still be handled (e.g. "Warning:
 *    expected 1 input but found 3, using first one"); log messages are for
 *    normal events (e.g. "Log: Variable 'temp' created"); and debug messages 
 *    report state, internal variables and other details for developers 
 *    (e.g. "Debug: 12 entries in DsTable").
 *
 *    In FCLib, a verbosity level implies all levels below it. For instance,
 *    asking for warning message will also report error messages. The message 
 *    hierarchy, from less to more verbose, is:
 * 
 *    FC_QUIET < FC_ERROR_MESSAGES < FC_WARNING_MESSAGES < FC_LOG_MESSAGES < 
 *    FC_DEBUG_MESSAGES
 */
typedef enum 
{
  FC_QUIET = 0,        /**< no messages */
  FC_ERROR_MESSAGES,   /**< error messages   */
  FC_WARNING_MESSAGES, /**< warning and error messages */
  FC_LOG_MESSAGES,     /**< log, warning and error messages  */ 
  FC_DEBUG_MESSAGES    /**< debug, log, warning and error messages */
} FC_VerbosityLevel;	

/**
 * \ingroup  DataTypes
 * \brief Possible function return values
 *
 * \description
 * 
 *   Enumerated constants are provided for return values. Generally,
 *   a function which returns '0' or 'FC_SUCCESS' is considered to
 *   have succeeded, and return values less than 1 indicates that
 *   an error occurred. FC_ERROR (-1) is a generic error category.
 *
 * \modifications 
 *    - 2003-SEP-22  W Koegler  Created
 */
typedef enum 
{
  FC_SUCCESS = 0,
  FC_ERROR = -1,
  FC_MEMORY_ERROR = -2,
  FC_INPUT_ERROR = -3,
  FC_FILE_IO_ERROR = -4
} FC_ReturnCode;

/**
 * \ingroup   DataTypes
 * \brief Possible element types
 *
 * \description
 */
typedef enum 
{
  FC_ET_UNKNOWN = 0,
  FC_ET_POINT, 
  FC_ET_LINE, 
  FC_ET_TRI, 
  FC_ET_QUAD, 
  FC_ET_TET, 
  FC_ET_PYRAMID,
  FC_ET_PRISM, 
  FC_ET_HEX,
  // arbitrary polygon?
  FC_ET_MIXED /**< This is needed to describe faces on pyramids & prisms */
} FC_ElementType;

/**
 * \ingroup   DataTypes
 * \brief Possible data types.
 */
typedef enum 
{
  FC_DT_UNKNOWN = 0,
  FC_DT_CHAR, 
  // FC_DT_UCHAR,   
  // FC_DT_SHORT, 
  // FC_DT_USHORT,
  FC_DT_INT, 
  // FC_DT_UINT,     
  // FC_DT_LONG, 
  // FC_DT_ULONG,
  // FC_DT_LONGLONG,    /* 64 bits */
  // FC_DT_ULONGLONG,   /* 64 bits */
  FC_DT_FLOAT, 
  FC_DT_DOUBLE 
  // FC_DT_LONGDOUBLE,  /* 128 bits */
} FC_DataType;

/**
 * \ingroup   DataTypes
 * \brief Possible ways of organizing data values.
 *
 * \description
 * 
 *    Specifies the organization and order of the components of a
 *    variable. For instance, FC_MT_SCALAR denotes a scalar
 *    field which will have a single value per data point.
 *    FC_MT_VECTOR denotes a vector field with n (numComponent)
 *    ordered components such as X, Y and Z for a vector
 *    with 3 components.
 */
typedef enum 
{
  FC_MT_UNKNOWN = 0,
  FC_MT_SCALAR, 
  FC_MT_VECTOR, 
  FC_MT_SYMTENSOR,  /**< not fully supported by most library routines */
  FC_MT_TENSOR      /**< not fully supported by most library routines */
} FC_MathType;

/**
 * \ingroup   DataTypes
 * \brief Possible association of the data to the mesh
 *
 * \description
 * 
 *     Specifies which type of sub-entity of a mesh something is
 *     associated with. E.g a variable with a vertex association
 *     will have data for each vertex.
 *
 * \todo ?Surely we can come up with a better name? entity type?
 */
typedef enum 
{
  FC_AT_UNKNOWN = 0,
  FC_AT_VERTEX, 
  FC_AT_EDGE,      /**< Cannot exist in 0D meshes, and is the same
		    as elements in 1D meshes */
  FC_AT_FACE,      /**< Cannot exist in 0D or 1D meshes, and is the
		    same as elements in 2D meshes */
  FC_AT_ELEMENT,   /**< Will be 0D, 1D, 2D or 3D depending on element type */
  FC_AT_WHOLE_MESH,    /**< Refers to the entire mesh */
  FC_AT_WHOLE_DATASET  /**< Refers to the entire dataset */
} FC_AssociationType;

//@}

// check that enum has a valid value
int fc_isVerbosityLevelValid(FC_VerbosityLevel verbosity);
int fc_isReturnCodeValid(FC_ReturnCode rc);
int fc_isElementTypeValid(FC_ElementType elemtype);
int fc_isDataTypeValid(FC_DataType datatype);
int fc_isMathTypeValid(FC_MathType mathtype);
int fc_isAssociationTypeValid(FC_AssociationType assoc);

// printable text for enums
char *fc_getVerbosityLevelText(FC_VerbosityLevel verbosity);
char *fc_getReturnCodeText(FC_ReturnCode rc);
char *fc_getElementTypeText(FC_ElementType elemtype);
char *fc_getDataTypeText(FC_DataType datatype);
char *fc_getMathTypeText(FC_MathType mathtype);
char *fc_getAssociationTypeText(FC_AssociationType assoc);

// information about FC_ElementType
int fc_getElementTypeTopoDim(FC_ElementType elemtype);
int fc_getElementTypeNumVertex(FC_ElementType elemtype);
int fc_getElementTypeNumEdge(FC_ElementType elemtype);
int fc_getElementTypeNumFace(FC_ElementType elemtype);
int fc_getElementTypeNumEntity(FC_ElementType elemtype, 
		 FC_AssociationType assoc);
FC_ElementType fc_getElementTypeFaceType(FC_ElementType elemtype);

// information about FC_DataType
int fc_sizeofDataType(FC_DataType datatype);

/**
 * \ingroup  Handles
 * \brief Dataset Handle.
 *
 * \description
 * 
 *    Handle used by library routines to refer to a specific dataset.
 *    See \ref Dataset for functions.
 */
typedef struct {
    int slotID;    /**< Index into the dataset table. */
    int uID;       /**< Unique identifier. */
} FC_Dataset;

/**
 * \ingroup  Handles
 * \brief Sequence Handle
 *
 * \description
 * 
 *    Handle used by library routines to refer to a specific
 *    sequence (e.g. time series). See \ref Sequence for functions. 
 */
typedef struct {
    int slotID;    /**< Index into the sequence table. */
    int uID;       /**< Unique identifier. */
} FC_Sequence;

/**
 * \ingroup  Handles
 * \brief Mesh Handle
 *
 * \description
 * 
 *    Handle used by library routines to refer to a specific mesh.
 *    See \ref Mesh for functions.
 */
typedef struct {
    int slotID;    /**< Index into the mesh table. */
    int uID;       /**< Unique identifier. */
} FC_Mesh;

/**
 * \ingroup  Handles
 * \brief Subset Handle
 *
 * \description
 * 
 *    Handle used by library routines to refer to a specific subset.
 *    See \ref Subset for functions.
 *
 * \modifications 
 *    - 05/11/04 WSK Created.
 */
typedef struct {
  int slotID;    /**< Index into the subset table. */
  int uID;       /**< Unique identifier. */
} FC_Subset;

/**
 * \ingroup  Handles
 * \brief Variable Handle
 *
 * \description
 * 
 *    Handle used by library routines to refer to a specific
 *    variable (or a specific step of a sequence variable, e.g.
 *    pressure at time ti). See \ref Variable for functions.
 */
typedef struct {
    int slotID;    /**< Index into the variable table. */
    int uID;       /**< Unique identifier. */
} FC_Variable;

/** 
 * \ingroup  Handles
 * \brief  Check two handles for equality
 *
 * \description
 *
 *    Evaluates to true (1) if handles are equal, and false(0) 
 *    if they are not. Works for all handle types but note
 *    that null handles of all types will be equal.
 *
 * \todo ?NULL handles of different types could have different negative
 *       negative sequence numbers to prevent null handles of different
 *       types from being equal.
 */
#define FC_HANDLE_EQUIV(handle1, handle2) \
          (((handle1).slotID == (handle2).slotID) && \
          ((handle1).uID == (handle2).uID))

// "Empty" handles

extern FC_Dataset FC_NULL_DATASET;
extern FC_Sequence FC_NULL_SEQUENCE;
extern FC_Mesh FC_NULL_MESH;
extern FC_Variable FC_NULL_VARIABLE;
extern FC_Subset FC_NULL_SUBSET;

#ifdef __cplusplus
}
#endif

#endif // _FC_BASE_H_

