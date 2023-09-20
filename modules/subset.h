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
 * \file subset.h
 * \brief Public declarations for \ref Subset module.
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/modules/subset.h,v $
 * $Revision: 1.59 $ 
 * $Date: 2006/09/26 21:31:32 $
 */

#ifndef _FC_SUBSET_H_
#define _FC_SUBSET_H_

#ifdef __cplusplus
extern "C" {
#endif

// create new subset from scratch
//-------------------------------
FC_ReturnCode fc_createSubset(FC_Mesh mesh, char *subsetName, 
                     FC_AssociationType subsetAssoc, FC_Subset *subset);

// other ways to get new subsets
//-------------------------------
FC_ReturnCode fc_copySubset(FC_Subset subset, FC_Mesh destinationMesh, 
                     char* newSubsetName, FC_Subset *newSubset);
FC_ReturnCode fc_copySubsetWithNewAssociation(FC_Subset subset, 
		      char* newSubsetName, FC_AssociationType new_assoc, 
                      int doStrict, FC_Subset* new_subset);
FC_ReturnCode fc_createSubsetComplement(FC_Subset subset, char* newSubsetName, 
                     FC_Subset *newSubset);
FC_ReturnCode fc_createSubsetIntersection(FC_Subset subset1, char* operation, 
                     FC_Subset subset2, char* newSubsetName, 
		     FC_Subset *newSubset);

// get subsets
//-------------------------------
FC_ReturnCode fc_getSubsets(FC_Mesh mesh, int *numSubset, FC_Subset **subsets);
FC_ReturnCode fc_getNumSubset(FC_Mesh, int* numSubset);
FC_ReturnCode fc_getSubsetByName(FC_Mesh mesh, char* subsetName, int *numSubset,
                     FC_Subset** subset);

// Subset comparison
//-------------------------------
FC_ReturnCode fc_doSubsetsIntersect(FC_Subset subset1, FC_Subset subset2, 
		     int* indicator);
FC_ReturnCode fc_isSubsetSuperset(FC_Subset subset_super, FC_Subset subset_sub,
		     int* indicator);

// change the name of a subset
//-------------------------------
FC_ReturnCode fc_changeSubsetName(FC_Subset subset, char* newName);

// add members to a subset  
//---------------------------------
FC_ReturnCode fc_addMemberToSubset(FC_Subset subset, int member);
FC_ReturnCode fc_addArrayMembersToSubset(FC_Subset subset, int numMember, 
					 int *membersArray);
FC_ReturnCode fc_addMaskMembersToSubset(FC_Subset subset, int maskLength, 
					int *mask);
// delete members from a subset
//---------------------------------
FC_ReturnCode fc_deleteMemberFromSubset(FC_Subset subset, int memberID); 
FC_ReturnCode fc_deleteArrayMembersFromSubset(FC_Subset subset, int numElement,
                                    int *membersArray);
FC_ReturnCode fc_deleteMaskMembersFromSubset(FC_Subset subset, 
					     int maskLength,  int *mask);

// release/delete subset
//-------------------------------
FC_ReturnCode fc_releaseSubset(FC_Subset subset);
FC_ReturnCode fc_deleteSubset(FC_Subset subset);

// Get Subset metadata
//-------------------------
int fc_isSubsetValid(FC_Subset subset);
FC_ReturnCode fc_getSubsetName(FC_Subset subset, char** subsetName);
FC_ReturnCode fc_getMeshFromSubset(FC_Subset subset, FC_Mesh *mesh);
FC_ReturnCode fc_getSubsetInfo(FC_Subset subset, int *numMember, 
			       int *maxNumMember, FC_AssociationType *assoc);
FC_ReturnCode fc_getSubsetNumMember(FC_Subset subset, int *numMember);
FC_ReturnCode fc_getSubsetMaxNumMember(FC_Subset subset, int *maxNumMember);
FC_ReturnCode fc_getSubsetAssociationType(FC_Subset subset,
			       FC_AssociationType *assoc);
//  get/query subset members
//---------------------------------
int fc_isMemberInSubset(FC_Subset subset, int memberID);
FC_ReturnCode fc_getSubsetMembersAsArray (FC_Subset subset, int *numElement, 
					  int **membersArray);
FC_ReturnCode fc_getSubsetMembersAsMask(FC_Subset subset, int *maskLength, 
					int **mask);

// Print
//---------------------------------
FC_ReturnCode fc_printSubset(FC_Subset subset, char* label, int print_members);
FC_ReturnCode fc_printSubsetOfMesh(FC_Subset, FC_Mesh mesh);
FC_ReturnCode fc_printSubsetOfVariable(FC_Subset, FC_Variable variable);

#ifdef __cplusplus
}
#endif

#endif // _FC_SUBSET_H_

