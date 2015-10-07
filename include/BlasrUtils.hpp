// Copyright (c) 2014-2015, Pacific Biosciences of California, Inc.
//
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted (subject to the limitations in the
// disclaimer below) provided that the following conditions are met:
//
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//
//  * Redistributions in binary form must reproduce the above
//    copyright notice, this list of conditions and the following
//    disclaimer in the documentation and/or other materials provided
//    with the distribution.
//
//  * Neither the name of Pacific Biosciences nor the names of its
//    contributors may be used to endorse or promote products derived
//    from this software without specific prior written permission.
//
// NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE
// GRANTED BY THIS LICENSE. THIS SOFTWARE IS PROVIDED BY PACIFIC
// BIOSCIENCES AND ITS CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
// WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
// OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
// USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
// OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
// SUCH DAMAGE.

// Author: Mark Chaisson


#ifndef _BLASR_INC_UTILS_HPP_
#define _BLASR_INC_UTILS_HPP_

#include "BlasrHeaders.h"

//----------------------MODIFY ALIGNMENTS--------------------------//
//FIXME: refactor class SequenceIndexDatabase
void AssignRefContigLocation(T_AlignmentCandidate &alignment, 
                             SequenceIndexDatabase<FASTQSequence> &seqdb,
                             DNASequence &genome);

//FIXME: refactor class SequenceIndexDatabase
void AssignRefContigLocations(vector<T_AlignmentCandidate*> &alignmentPtrs,
                              SequenceIndexDatabase<FASTQSequence> &seqdb,
                              DNASequence &genome);

template<typename T_RefSequence>
//FIXME: refactor class SequenceIndexDatabase
void AssignGenericRefContigName(vector<T_AlignmentCandidate*> &alignmentPtrs,
                                T_RefSequence &genome);

//FIXME: move to class ReadAlignments
void StoreRankingStats(vector<T_AlignmentCandidate*> &alignments,
                       VarianceAccumulator<float> &accumPValue, 
                       VarianceAccumulator<float> &accumWeight);

//FIXME: mapQV should be assigned when alignments are created.
void AssignMapQV(vector<T_AlignmentCandidate*> &alignmentPtrs);

//FIXME: move to class ReadAlignments
void ScaleMapQVByClusterSize(T_AlignmentCandidate &alignment, 
                             MappingParameters &params);

void StoreMapQVs(SMRTSequence &read,
                 vector<T_AlignmentCandidate*> &alignmentPtrs, 
                 MappingParameters &params);


//--------------------SEARCH & CHECK ALIGNMENTS-------------------//
//FIXME: move to class ReadAlignments
template<typename T_Sequence>
bool CheckForSufficientMatch(T_Sequence &read, 
                             vector<T_AlignmentCandidate*> &alignmentPtrs,
                             MappingParameters &params);

//FIXME: move to class ReadAlignments
int FindMaxLengthAlignment(vector<T_AlignmentCandidate*> alignmentPtrs,
                           int &maxLengthIndex);

//FIXME: move to class T_AlignmentCandidate
void SumMismatches(SMRTSequence &read,
                   T_AlignmentCandidate &alignment,
                   int mismatchScore,
                   int fullIntvStart, int fullIntvEnd,
                   int &sum);

//FIXME: move to class T_AlignmentCandidate
/// \returns whether two alignments overlap by more than minPcercentOverlap%
bool AlignmentsOverlap(T_AlignmentCandidate &alnA,
                       T_AlignmentCandidate &alnB,
                       float minPercentOverlap);

/// \Partition overlapping alignments.
void PartitionOverlappingAlignments(vector<T_AlignmentCandidate*> &alignmentPtrs,
                                    vector<set<int> > &partitions,
                                    float minOverlap);


//--------------------FILTER ALIGNMENTS---------------------------//
//FIXME: move to class T_AlignmentCandidate and ReadAlignments
int RemoveLowQualitySDPAlignments(int readLength, 
                                  vector<T_AlignmentCandidate*> &alignmentPtrs, 
                                  MappingParameters &params);

//FIXME: move to class ReadAlignments
template<typename T_Sequence>
int RemoveLowQualityAlignments(T_Sequence &read,
                               vector<T_AlignmentCandidate*> &alignmentPtrs,
                               MappingParameters &params);

//FIXME: move to class ReadAlignments
int RemoveOverlappingAlignments(vector<T_AlignmentCandidate*> &alignmentPtrs, 
                                MappingParameters &params);

// FIXME: move to class ReadAlignments
// Delete all alignments from index startIndex in vector, inclusive.
void DeleteAlignments(vector<T_AlignmentCandidate*> &alignmentPtrs,
                      int startIndex=0);

//--------------------REFINE ALIGNMENTS---------------------------//
template<typename T_RefSequence, typename T_Sequence>
void RefineAlignment(vector<T_Sequence*> &bothQueryStrands,
                     T_RefSequence &genome,
                     T_AlignmentCandidate  &alignmentCandidate,
                     MappingParameters &params,
                     MappingBuffers &mappingBuffers);


template<typename T_RefSequence, typename T_Sequence>
void RefineAlignments(vector<T_Sequence*> &bothQueryStrands,
                      T_RefSequence &genome,
                      vector<T_AlignmentCandidate*> &alignmentPtrs, 
                      MappingParameters &params, 
                      MappingBuffers &mappingBuffers);


//--------------------PRINT ALIGNMENTS---------------------------//
vector<T_AlignmentCandidate*>
SelectAlignmentsToPrint(vector<T_AlignmentCandidate*> alignmentPtrs,
                        MappingParameters & params,
                        const int & associatedRandInt);

//
// The full read is not the subread, and does not have masked off characters.
//
void PrintAlignment(T_AlignmentCandidate &alignment, 
                    SMRTSequence &fullRead, 
                    MappingParameters &params, 
                    AlignmentContext &alignmentContext, 
                    ostream &outFile, 
#ifdef USE_PBBAM
                    PacBio::BAM::BamWriter * bamWriterPtr
#endif
                    );

// Print all alignments in vector<T_AlignmentCandidate*> alignmentPtrs
void PrintAlignments(vector<T_AlignmentCandidate*> alignmentPtrs,
                     SMRTSequence &read,
                     MappingParameters &params, ostream &outFile, 
                     AlignmentContext alignmentContext,
#ifdef USE_PBBAM
                     PacBio::BAM::BamWriter * bamWriterPtr,
#endif
                     MappingSemaphores & semaphores);

void PrintAlignmentPtrs(vector <T_AlignmentCandidate*> & alignmentPtrs,
                        ostream & out = cout);

// Print all alignments for subreads in allReadAlignments.
// Input: 
//   allReadAlignments - contains a set of subreads, each of which
//                       is associated with a group of alignments.
//   alignmentContext  - an alignment context of each subread used
//                       for printing in SAM format.
//   params            - mapping parameters.
// Output:
//   outFilePtr        - where to print alignments for subreads.
//   unalignedFilePtr  - where to print sequences for unaligned subreads.
void PrintAllReadAlignments(ReadAlignments & allReadAlignments,
                            AlignmentContext & alignmentContext,
                            ostream & outFilePtr,
                            ostream & unalignedFilePtr,
                            MappingParameters & params,
#ifdef USE_PBBAM
                            vector<SMRTSequence> & subreads, 
                            PacBio::BAM::BamWriter * bamWriterPtr,
#endif
                            MappingSemaphores & semaphores);

#include "BlasrUtilsImpl.hpp"
#endif
