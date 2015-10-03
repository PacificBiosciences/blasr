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
#ifndef __BLASR_ALIGN_HPP_
#define __BLASR_ALIGN_HPP_

#include "BlasrHeaders.h"
#include "BlasrMiscs.hpp"

//------------------MAP READS---------------------------------//
template<typename T_Sequence, typename T_RefSequence, typename T_SuffixArray, typename T_TupleCountTable>
void MapRead(T_Sequence &read, T_Sequence &readRC, T_RefSequence &genome, 
             T_SuffixArray &sarray, 
             BWT &bwt,
             SeqBoundaryFtr<FASTQSequence> &seqBoundary, 
             T_TupleCountTable &ct,
             SequenceIndexDatabase<FASTQSequence> &seqdb,
             MappingParameters &params,
             MappingMetrics    &metrics,
             vector<T_AlignmentCandidate*> &alignmentPtrs, 
             MappingBuffers &mappingBuffers,
             MappingIPC *mapData,
             MappingSemaphores & semaphores);
/*
void MapReads(MappingData<T_SuffixArray, T_GenomeSequence, T_Tuple> *mapData);
*/

//------------------MAKE ALIGNMENTS---------------------------//
template<typename T_TargetSequence, typename T_QuerySequence, typename TDBSequence>
void AlignIntervals(T_TargetSequence &genome, T_QuerySequence &read, T_QuerySequence &rcRead,
                    WeightedIntervalSet &weightedIntervals,
                    int mutationCostMatrix[][5], 
                    int ins, int del, int sdpTupleSize,
                    int useSeqDB, SequenceIndexDatabase<TDBSequence> &seqDB,
                    vector<T_AlignmentCandidate*> &alignments,
                    MappingParameters &params,
                    MappingBuffers &mappingBuffers,
                    int procId=0);
  
template<typename T_RefSequence, typename T_Sequence>
void PairwiseLocalAlign(T_Sequence &qSeq, T_RefSequence &tSeq, 
                        int k, 
                        MappingParameters &params, T_AlignmentCandidate &alignment, 
                        MappingBuffers &mappingBuffers,
                        AlignmentType alignType=Global);

// Extend target aligned sequence of the input alignement to both ends
// by flankSize bases. Update alignment->tAlignedSeqPos, 
// alignment->tAlignedSeqLength and alignment->tAlignedSeq.
void FlankTAlignedSeq(T_AlignmentCandidate * alignment,
                      SequenceIndexDatabase<FASTQSequence> &seqdb,
                      DNASequence & genome,
                      int flankSize);

// Align a subread of a SMRT sequence to target sequence of an alignment.
// Input:
//   subread         - a subread of a SMRT sequence.
//   unrolledRead    - the full SMRT sequence.
//   alignment       - an alignment.
//   passDirection   - whether or not the subread has the 
//                     same direction as query of the alignment.
//                     0 = true, 1 = false. 
//   subreadInterval - [start, end) interval of the subread in the 
//                     SMRT read.
//   subreadIndex    - index of the subread in allReadAlignments.
//   params          - mapping paramters.
// Output:
//   allReadAlignments - where the sequence and alignments of the 
//                       subread are saved.
//   threadOut         - an out stream for debugging the current thread. 
void AlignSubreadToAlignmentTarget(ReadAlignments & allReadAlignments,
        SMRTSequence & subread, SMRTSequence & unrolledRead,
        T_AlignmentCandidate * alignment,
        int passDirection, ReadInterval & subreadInterval,
        int subreadIndex, 
        MappingParameters & params, 
        MappingBuffers & mappingBuffers,
        ostream & threadOut);

#include "BlasrAlignImpl.hpp"
#endif
