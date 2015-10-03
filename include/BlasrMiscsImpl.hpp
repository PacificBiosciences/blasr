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

#ifndef _BLASR_MISCS_IMPL_HPP_
#define _BLASR_MISCS_IMPL_HPP_

template<typename T_Sequence>
bool GetNextReadThroughSemaphore(ReaderAgglomerate &reader, 
                                 MappingParameters &params,
                                 T_Sequence &read, 
                                 string & readGroupId, 
                                 int & associatedRandInt,
                                 MappingSemaphores & semaphores)
{
    // Wait on a semaphore
    if (params.nProc > 1) {
#ifdef __APPLE__
        sem_wait(semaphores.reader);
#else
        sem_wait(&semaphores.reader);
#endif
    }

    bool returnValue = true;
    //
    // CCS Reads are read differently from other reads.  Do static casting here
    // of this.
    //
    if (reader.GetNext(read, associatedRandInt) == 0) {
        returnValue = false;
    }

    //
    // Set the read group id before releasing the semaphore, since other
    // threads may change the reader object to a new read group before
    // sending this alignment out to printing. 
    readGroupId = reader.readGroupId;

    if (params.nProc > 1) {
#ifdef __APPLE__
        sem_post(semaphores.reader);
#else
        sem_post(&semaphores.reader);
#endif
    }
    return returnValue;
}



bool ReadHasMeaningfulQualityValues(FASTQSequence &sequence) 
{
    if (sequence.qual.Empty() == true) {
        return 0;
    }
    else {
        int numZero=0, numNonZero=0;
        if (sequence.qual.data == NULL) {
            return false;
        }
        numZero = CountZero(sequence.qual.data, sequence.length);
        numNonZero = sequence.length - numZero;
        int subNumZero = 0, subNonZero = 0;

        if (sequence.substitutionQV.data == NULL) {
            return false;
        }
        subNumZero = CountZero(sequence.substitutionQV.data, sequence.length);
        subNonZero = sequence.length - subNumZero;

        if (numZero < 0.5*numNonZero and subNumZero < 0.5 * subNonZero) {
            return true;
        }
        else {
            return false;
        }
    }
}

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
                           MappingParameters & params)
{
    int start = subreadInterval.start;
    int end   = subreadInterval.end;

    assert(smrtRead.length >= subreadSequence.length);
    smrtRead.MakeSubreadAsMasked(subreadSequence, start, end); 

    if (!params.preserveReadTitle) {
        smrtRead.SetSubreadTitle(subreadSequence, 
                subreadSequence.SubreadStart(),
                subreadSequence.SubreadEnd());
    }
    else {
        subreadSequence.CopyTitle(smrtRead.title);
    }
    subreadSequence.zmwData = smrtRead.zmwData;
}

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
                   SMRTSequence & smrtRead) 
{
    assert(smrtRead.length >= subreadSequence.length);
    // Reverse complement sequence of the subread.
    subreadSequence.MakeRC(subreadSequenceRC);
    // Update start and end positions of subreadSequenceRC in the 
    // coordinate of reverse compelement sequence of the SMRT read.
    subreadSequenceRC.SubreadStart(smrtRead.length - subreadSequence.SubreadEnd());
    subreadSequenceRC.SubreadEnd  (smrtRead.length - subreadSequence.SubreadStart());
    subreadSequenceRC.zmwData     = smrtRead.zmwData;
}

int CountZero(unsigned char *ptr, int length) 
{
    int i;
    int nZero = 0;
    for (i = 0; i < length; i++) {
        if (ptr[i] == 0) { ++nZero; }
    }
    return nZero;
}

bool FirstContainsSecond(DNALength aStart, DNALength aEnd, DNALength bStart, DNALength bEnd) 
{
    return ((bStart > aStart and bEnd <= aEnd) or
            (bStart >= aStart and bEnd < aEnd));
}

float ComputePMatch(float accuracy, int anchorLength) 
{
    assert(anchorLength >= 0);
    if (anchorLength == 0) { 
        return 0;
    }
    else {
        return pow(accuracy,anchorLength);
    }
}

//
// Assume the number of mismatches in a row follow a geometric distribution.
//
void GeometricDistributionSummaryStats(float pSuccess,
                                       float &mean, float &variance) 
{
    mean = 1/pSuccess;
    variance = (1-pSuccess)/ (pow(pSuccess,2));
}



int ComputeExpectedWaitingBases(float mean, float variance, float certainty) 
{
    float nStdDev;
    assert(FindQNorm(certainty, nStdDev) != false);
    return mean + sqrt(variance) * nStdDev;
}

#endif
