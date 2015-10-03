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
#ifndef _BLASR_MISCS_HPP_
#define _BLASR_MISCS_HPP_

#include "BlasrHeaders.h"

//-------------------------Fetch Reads----------------------------//
template<typename T_Sequence>
bool GetNextReadThroughSemaphore(ReaderAgglomerate &reader, 
                                 MappingParameters &params,
                                 T_Sequence &read, 
                                 AlignmentContext &context, 
                                 int & associatedRandInt,
                                 MappingSemaphores & semaphores);

//---------------------MAKE & CHECK READS-------------------------//
//FIXME: move to SMRTSequence
bool ReadHasMeaningfulQualityValues(FASTQSequence &sequence);

//FIXME: Move to SMRTSequence
// Given a SMRT sequence and a subread interval, make the subread.
// Input: 
//   smrtRead         - a SMRT sequence 
//   subreadInterval  - a subread interval 
//   params           - mapping parameters
// Output: 
//   subreadSequence - the constructed subread
void MakeSubreadOfInterval(SMRTSequence & subreadSequence, 
                           SMRTSequence & smrtRead, 
                           ReadInterval & subreadInterval, 
                           MappingParameters & params);

//FIXME: Move to SMRTSequence
// Given a SMRT sequence and one of its subreads, make the 
// reverse complement of the subread in the coordinate of the
// reverse complement sequence of the SMRT sequence.
// Input: 
//   smrtRead          - a SMRT read
//   subreadSequence   - a subread of smrtRead     
// Output:
//   subreadSequenceRC - the reverse complement of the subread
//                       in the coordinate of the reverse 
//                       complement of the SMRT read.
void MakeSubreadRC(SMRTSequence & subreadSequenceRC,
                   SMRTSequence & subreadSequence,
                   SMRTSequence & smrtRead);


//-------------------------MISC-----------------------------------//
int CountZero(unsigned char *ptr, int length);

//FIXME: delete unused function
bool FirstContainsSecond(DNALength aStart, DNALength aEnd, DNALength bStart, DNALength bEnd);

//FIXME: delete unused function
// Assume the number of mismatches in a row follow a geometric distribution.
void GeometricDistributionSummaryStats(float pSuccess,
                                       float &mean, float &variance);

//FIXME: delete unused function
int ComputeExpectedWaitingBases(float mean, float variance, float certainty);

//FIXME: delete unused function
float ComputePMatch(float accuracy, int anchorLength);

#include "BlasrMiscsImpl.hpp"

#endif
