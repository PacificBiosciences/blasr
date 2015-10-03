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

#ifndef __BLASR_READ_ALIGNMENTS__
#define __BLASR_READ_ALIGNMENTS__

#include <string>
#include <iostream>
#include <vector>
#include "SMRTSequence.hpp"
#include "datastructures/alignment/AlignmentCandidate.hpp"

using namespace std;

class ReadAlignments {
public:
/*
   This class stores the alignments from a read.  A read may be
   aligned in several different modes:
   1. Fullread    - Treat the read as a unit from start to end
   2. Subread     - Align each subread independently
   3. CCSDeNovo   - Only align the CCS sequence from a read
   4. CCSAllPass  - Align the de novo ccs sequences and then the
   subreads to where the denovo ccs aligned.
   5. CCSFullPass - Same as allpass, except using only complete
   subreads.
   6. ZmwSubreads - Align subreads of each zmw to where the longest 
   subread of the zmw aligned to.

   The alignments are a raggad array of n sequences; n is 1 for cases 
   1 and 3, the number of subreads for cases 2 and 4, and the number
   of full length passes for case 5.

   A ReadAligments class must only have alignments for a single type
   of read in it.

*/
    vector<vector<T_AlignmentCandidate*> > subreadAlignments;
    vector<SMRTSequence> subreads;
    AlignMode alignMode;
    SMRTSequence read;

    inline int GetNAlignedSeq();

    inline bool AllSubreadsHaveAlignments();

    inline void Clear();

    inline void Resize(int nSeq);

    inline void CheckSeqIndex(int seqIndex);
        
    inline void SetSequence(int seqIndex, SMRTSequence &seq);

    inline void AddAlignmentForSeq(int seqIndex, T_AlignmentCandidate *alignmentPtr);

    inline void AddAlignmentsForSeq(int seqIndex, vector<T_AlignmentCandidate*> &seqAlignmentPtrs);

    // Copy all T_AlignmentCandidate objects (to which subreadAlignment[seqIndex]
    // is pointing) to newly created objects, and then return pointers to the new
    // objects.
    inline vector<T_AlignmentCandidate*> CopySubreadAlignments(int seqIndex);

    inline void Print(ostream &out=cout);

    inline ~ReadAlignments();
};


inline int ReadAlignments::GetNAlignedSeq() {
    return subreadAlignments.size();
}

inline bool ReadAlignments::AllSubreadsHaveAlignments() {
    int i, nAlignedSeq;
    nAlignedSeq = subreadAlignments.size();
    for (i = 0; i < nAlignedSeq; i++) {
        if (subreadAlignments[i].size() == 0) {
            return false;
        }
    }
    return true;
}

inline void ReadAlignments::Clear() {
    int i;
    int nAlignedSeq;
    for (i = 0, nAlignedSeq = subreadAlignments.size(); i < nAlignedSeq; i++) {
        int nAlignments;
        int a;
        for (a = 0, nAlignments = subreadAlignments[i].size(); a < nAlignments; a++) {
            delete subreadAlignments[i][a];
        }
        subreadAlignments[i].clear();
    }

    for (i = 0, nAlignedSeq = subreads.size(); i< nAlignedSeq; i++) {
        subreads[i].Free();
    }
    subreadAlignments.clear();
    read.Free();
}

inline void ReadAlignments::Resize(int nSeq) {
    subreadAlignments.resize(nSeq);
    subreads.resize(nSeq);
}

inline void ReadAlignments::CheckSeqIndex(int seqIndex) {
    if ( seqIndex < 0 or seqIndex >= int(subreads.size()) ) {
        cout << "ERROR, adding a sequence to an unallocated position." 
            << endl;
        assert(0);
    }
}

inline void ReadAlignments::SetSequence(int seqIndex, SMRTSequence &seq) {
    CheckSeqIndex(seqIndex);
    subreads[seqIndex] = seq;
}

inline void ReadAlignments::AddAlignmentForSeq(int seqIndex, T_AlignmentCandidate *alignmentPtr) {
    CheckSeqIndex(seqIndex);
    subreadAlignments[seqIndex].push_back(alignmentPtr);
}

inline void ReadAlignments::AddAlignmentsForSeq(int seqIndex, vector<T_AlignmentCandidate*> &seqAlignmentPtrs) {
    CheckSeqIndex(seqIndex);
    subreadAlignments[seqIndex].insert(subreadAlignments[seqIndex].end(), seqAlignmentPtrs.begin(), seqAlignmentPtrs.end());
}

inline vector<T_AlignmentCandidate*> ReadAlignments::CopySubreadAlignments(int seqIndex) {
    vector<T_AlignmentCandidate*> ret;
    for (int i=0; i<int(subreadAlignments[seqIndex].size()); i++) {
        T_AlignmentCandidate * q = new T_AlignmentCandidate();
        *q = *(subreadAlignments[seqIndex][i]);
        ret.push_back(q);
    }
    return ret;
}

inline void ReadAlignments::Print(ostream &out=cout) { 
    out << "A ReadAlignments object with " 
        << subreadAlignments.size()
        << " groups of subread alignments." << endl;
    for (int i = 0; i < int(subreadAlignments.size()); i++) {
        out << "  subreadAlignment group [" << i << "/" 
            << subreadAlignments.size() << "] has "
            << subreadAlignments[i].size() << " alignments." << endl;
        for(int j = 0; j < int(subreadAlignments[i].size()); j++) {
            out << "    [" << i << "][" << j << "/" 
                << subreadAlignments[i].size() << "]" << endl;
            subreadAlignments[i][j]->Print(out);
        } 
    }
    out << "  read: ";
    read.Print(out);
    out << endl << endl;
}

inline ReadAlignments::~ReadAlignments() {
    read.Free();
}

#endif
