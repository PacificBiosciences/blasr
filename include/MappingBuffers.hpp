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
#ifndef __BLASR_MAPPING_BUFFERS__
#define __BLASR_MAPPING_BUFFERS__

#include <vector>
#include "tuples/DNATuple.hpp"
#include "tuples/TupleList.hpp"
#include "algorithms/alignment/sdp/SDPFragment.hpp"
#include "algorithms/anchoring/BasicEndpoint.hpp"
#include "datastructures/anchoring/ClusterList.hpp"
#include "datastructures/anchoring/MatchPos.hpp"

using namespace std;


//
// Define a list of buffers that are meant to grow to high-water
// marks, and not shrink down past that.   The memory is reused rather
// than having multiple calls to new.
//
class MappingBuffers {
public:
    vector<int> hpInsScoreMat, insScoreMat;
    vector<int> kbandScoreMat;
    vector<Arrow> hpInsPathMat, insPathMat;
    vector<Arrow> kbandPathMat;
    vector<int>   scoreMat;
    vector<Arrow> pathMat;
    vector<int>  affineScoreMat;
    vector<Arrow> affinePathMat;
    vector<ChainedMatchPos> matchPosList;
    vector<ChainedMatchPos> rcMatchPosList;
    vector<BasicEndpoint<ChainedMatchPos> > globalChainEndpointBuffer;
    vector<Fragment> sdpFragmentSet, sdpPrefixFragmentSet, sdpSuffixFragmentSet;
    TupleList<PositionDNATuple> sdpCachedTargetTupleList;
    TupleList<PositionDNATuple> sdpCachedTargetPrefixTupleList;
    TupleList<PositionDNATuple> sdpCachedTargetSuffixTupleList;
    std::vector<int> sdpCachedMaxFragmentChain;
    vector<double> probMat;
    vector<double> optPathProbMat;
    vector<float>  lnSubPValueMat;
    vector<float>  lnInsPValueMat;
    vector<float>  lnDelPValueMat;
    vector<float>  lnMatchPValueMat;
    vector<int>    clusterNumBases;
    ClusterList    clusterList;
    ClusterList    revStrandClusterList;

    void Reset(void);
};


inline void MappingBuffers::Reset(void) {
    vector<int>().swap(hpInsScoreMat);
    vector<int>().swap(insScoreMat);
    vector<int>().swap(kbandScoreMat);
    vector<Arrow>().swap(hpInsPathMat);
    vector<Arrow>().swap(insPathMat);
    vector<Arrow>().swap(kbandPathMat);
    vector<int>().swap(scoreMat);
    vector<Arrow>().swap(pathMat);
    vector<ChainedMatchPos>().swap(matchPosList);
    vector<ChainedMatchPos>().swap(rcMatchPosList);
    vector<BasicEndpoint<ChainedMatchPos> >().swap(globalChainEndpointBuffer);
    vector<Fragment>().swap(sdpFragmentSet);
    vector<Fragment>().swap(sdpPrefixFragmentSet);
    vector<Fragment>().swap(sdpSuffixFragmentSet);
    sdpCachedTargetTupleList.Reset();
    sdpCachedTargetPrefixTupleList.Reset();
    sdpCachedTargetSuffixTupleList.Reset();
    vector<int>().swap(sdpCachedMaxFragmentChain);
    vector<double>().swap(probMat);
    vector<double>().swap(optPathProbMat);
    vector<float>().swap(lnSubPValueMat);
    vector<float>().swap(lnInsPValueMat);
    vector<float>().swap(lnDelPValueMat);
    vector<float>().swap(lnMatchPValueMat);
    vector<int>().swap(clusterNumBases);
}

#endif
