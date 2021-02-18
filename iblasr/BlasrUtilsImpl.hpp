// Author: Mark Chaisson
#pragma once

#include "BlasrAlign.hpp"

//----------------------MODIFY ALIGNMENTS--------------------------//
void AssignRefContigLocation(T_AlignmentCandidate &alignment,
                             SequenceIndexDatabase<FASTQSequence> &seqdb, DNASequence &genome)
{
    //
    // If the sequence database is used, the start position of
    // the alignment is relative to the start of the chromosome,
    // not the entire index.  Subtract off the start position of
    // the chromosome to get the true position.
    //
    DNALength forwardTPos;
    int seqDBIndex;
    if (alignment.tStrand == 0) {
        forwardTPos = alignment.tAlignedSeqPos;
        seqDBIndex = seqdb.SearchForIndex(forwardTPos);
        alignment.tAlignedSeqPos -= seqdb.seqStartPos[seqDBIndex];
    } else {
        //
        // Flip coordinates into forward strand in order to find the boundaries
        // of the contig, then reverse them in order to find offset.
        //

        // Find the reverse complement coordinate of the index of the last aligned base.
        assert(alignment.tAlignedSeqLength > 0);
        forwardTPos =
            genome.MakeRCCoordinate(alignment.tAlignedSeqPos + alignment.tAlignedSeqLength - 1);
        seqDBIndex = seqdb.SearchForIndex(forwardTPos);

        //
        // Find the reverse comlement coordinate of the last base of this
        // sequence.  This would normally be the start of the next contig
        // -1 to get the length, but since an 'N' is added between every
        // pair of sequences, this is -2.
        //
        DNALength reverseTOffset;
        reverseTOffset = genome.MakeRCCoordinate(seqdb.seqStartPos[seqDBIndex + 1] - 2);
        alignment.tAlignedSeqPos -= reverseTOffset;
    }
}

void AssignRefContigLocations(std::vector<T_AlignmentCandidate *> &alignmentPtrs,
                              SequenceIndexDatabase<FASTQSequence> &seqdb, DNASequence &genome)
{

    UInt i;
    for (i = 0; i < alignmentPtrs.size(); i++) {
        T_AlignmentCandidate *aref = alignmentPtrs[i];
        AssignRefContigLocation(*aref, seqdb, genome);
    }
}

template <typename T_RefSequence>
void AssignGenericRefContigName(std::vector<T_AlignmentCandidate *> &alignmentPtrs,
                                T_RefSequence &genome)
{
    UInt i;
    for (i = 0; i < alignmentPtrs.size(); i++) {
        T_AlignmentCandidate *aref = alignmentPtrs[i];
        aref->tName = genome.title;
    }
}

void StoreRankingStats(std::vector<T_AlignmentCandidate *> &alignments,
                       VarianceAccumulator<float> &accumPValue,
                       VarianceAccumulator<float> &accumWeight)
{
    int i;
    for (i = 0; i < int(alignments.size()); i++) {
        alignments[i]->pvalVariance = accumPValue.GetVariance();
        alignments[i]->pvalNStdDev = accumPValue.GetNStdDev(alignments[i]->clusterScore);
        alignments[i]->weightVariance = accumWeight.GetVariance();
        alignments[i]->weightNStdDev = accumWeight.GetNStdDev(alignments[i]->clusterWeight);
    }
}

void AssignMapQV(std::vector<T_AlignmentCandidate *> &alignmentPtrs)
{
    int i;
    int mapQV = 1;
    if (alignmentPtrs.size() > 1 and alignmentPtrs[0]->score == alignmentPtrs[1]->score) {
        // the top two alignments have the same score, don't consider them as mapped.
        mapQV = 0;
    }

    for (i = 0; i < int(alignmentPtrs.size()); i++) {
        alignmentPtrs[i]->mapQV = mapQV;
    }
}

void ScaleMapQVByClusterSize(T_AlignmentCandidate &alignment, MappingParameters &params)
{
    if (alignment.numSignificantClusters > int(params.nCandidates)) {
        alignment.mapQV = Phred((1 - InversePhred(alignment.mapQV)) *
                                ((float)params.nCandidates / alignment.numSignificantClusters));
    } else if (alignment.numSignificantClusters == 0) {
        alignment.mapQV = 0;
    }
}

void StoreMapQVs(SMRTSequence &read, std::vector<T_AlignmentCandidate *> &alignmentPtrs,
                 MappingParameters &params)
{
    //
    // Only weight alignments for mapqv against eachother if they are overlapping.
    //
    int a;
    std::vector<std::set<int> >
        partitions;  // Each set contains alignments that overlap on the read.
    DistanceMatrixScoreFunction<DNASequence, FASTQSequence> distScoreFn;
    distScoreFn.del = params.deletion;
    distScoreFn.ins = params.insertion;
    // bug 24363, set affineOpen and affineExtend for distScoreFn
    distScoreFn.affineOpen = params.affineOpen;
    distScoreFn.affineExtend = params.affineExtend;
    distScoreFn.InitializeScoreMatrix(SMRTLogProbMatrix);
    IDSScoreFunction<DNASequence, FASTQSequence> idsScoreFn;
    idsScoreFn.ins = params.insertion;
    idsScoreFn.del = params.deletion;
    idsScoreFn.affineExtend = params.affineExtend;
    idsScoreFn.affineOpen = params.affineOpen;
    idsScoreFn.substitutionPrior = params.substitutionPrior;
    idsScoreFn.globalDeletionPrior = params.globalDeletionPrior;

    //
    // Rescore the alignment so that it uses probabilities.
    //
    for (a = 0; a < int(alignmentPtrs.size()); a++) {
        if (params.ignoreQualities == false) {
            // bug 24363, pass -affineAlign to compute correct alignment score.
            alignmentPtrs[a]->probScore =
                -ComputeAlignmentScore(*alignmentPtrs[a], alignmentPtrs[a]->qAlignedSeq,
                                       alignmentPtrs[a]->tAlignedSeq, idsScoreFn,
                                       params.affineAlign) /
                10.0;
        } else {
            alignmentPtrs[a]->probScore =
                -ComputeAlignmentScore(*alignmentPtrs[a], alignmentPtrs[a]->qAlignedSeq,
                                       alignmentPtrs[a]->tAlignedSeq, distScoreFn,
                                       params.affineAlign) /
                10.0;
        }
    }
    PartitionOverlappingAlignments(alignmentPtrs, partitions,
                                   params.minFractionToBeConsideredOverlapping);

    int p;
    std::set<int>::iterator partIt, partEnd;

    //
    // For each partition, store where on the read it begins, and where
    // it ends.
    //
    std::vector<int> partitionBeginPos, partitionEndPos;
    partitionBeginPos.resize(partitions.size());
    partitionEndPos.resize(partitions.size());
    std::fill(partitionBeginPos.begin(), partitionBeginPos.end(), -1);
    std::fill(partitionEndPos.begin(), partitionEndPos.end(), -1);
    std::vector<char> assigned;
    assigned.resize(alignmentPtrs.size());
    std::fill(assigned.begin(), assigned.end(), false);

    for (p = 0; p < int(partitions.size()); p++) {
        partEnd = partitions[p].end();
        int alnStart, alnEnd;

        if (partitions[p].size() > 0) {
            partIt = partitions[p].begin();
            alignmentPtrs[*partIt]->GetQInterval(alnStart, alnEnd);
            partitionBeginPos[p] = alnStart;
            partitionEndPos[p] = alnEnd;
            ++partIt;
            partEnd = partitions[p].end();
            for (; partIt != partEnd; ++partIt) {
                //  Comment out because all reads are now in the forward strand.
                //  alignmentPtrs[*partIt]->GetQInterval(alnStart, alnEnd, convertToForwardStrand);
                alignmentPtrs[*partIt]->GetQInterval(alnStart, alnEnd);
                if (alnEnd - alnStart > partitionEndPos[p] - partitionBeginPos[p]) {
                    partitionBeginPos[p] = alnStart;
                    partitionEndPos[p] = alnEnd;
                }
            }
        }
    }

    //
    // For each partition, determine the widest parts of the read that
    // are aligned in the partition.  All alignments will be extended to
    // the end of the widest parts of the partition.
    //
    const static bool convertToForwardStrand = true;

    UInt i;

    //
    // For now, just use the alignment score as the probability score.
    // Although it is possible to use the full forward probability, for
    // the most part it is pretty much the same as the Vitterbi
    // probability, but it takes a lot longer to compute.
    //

    //
    // Now estimate what the alignment scores would be if they were
    // extended past the ends of their current alignment.
    //

    for (p = 0; p < int(partitions.size()); p++) {
        partEnd = partitions[p].end();
        int alnStart, alnEnd;
        for (partIt = partitions[p].begin(); partitions[p].size() > 0 and partIt != partEnd;
             ++partIt) {
            int mismatchSum = 0;
            alignmentPtrs[*partIt]->GetQInterval(alnStart, alnEnd, convertToForwardStrand);
            if (alnStart - partitionBeginPos[p] > MAPQV_END_ALIGN_WIGGLE or
                partitionEndPos[p] - alnEnd > MAPQV_END_ALIGN_WIGGLE) {
                // bug 24363, use updated SumMismatches to compute mismatch score when
                // no QV is available.
                SumMismatches(read, *alignmentPtrs[*partIt], 15, partitionBeginPos[p],
                              partitionEndPos[p], params, mismatchSum);
            }
            //
            // Random sequence can be aligned with about 50% similarity due
            // to optimization, so weight the qv sum
            //
            alignmentPtrs[*partIt]->probScore += -(mismatchSum)*0.5;
        }
    }

    //
    // Determine mapqv by summing qvscores in partitions

    float mapQVDenominator = 0;
    for (p = 0; p < int(partitions.size()); p++) {
        std::set<int>::iterator nextIt;
        if (partitions[p].size() == 0) {
            continue;
        }
        int index = *partitions[p].begin();

        mapQVDenominator = alignmentPtrs[index]->probScore;

        if (partitions[p].size() > 1) {
            partIt = partitions[p].begin();
            partEnd = partitions[p].end();
            ++partIt;

            for (; partIt != partEnd; ++partIt) {
                index = *partIt;
                mapQVDenominator = LogSumOfTwo(mapQVDenominator, alignmentPtrs[index]->probScore);
            }
        }

        for (partIt = partitions[p].begin(); partIt != partitions[p].end(); ++partIt) {
            //
            // If only one alignment is found, assume maximum mapqv.
            //
            assigned[*partIt] = true;
            if (partitions[p].size() == 1) {
                alignmentPtrs[*partIt]->mapQV = MAX_PHRED_SCORE;
            }

            //
            // Look for overflow.
            //
            else if (alignmentPtrs[*partIt]->probScore - mapQVDenominator < -20) {
                alignmentPtrs[*partIt]->mapQV = 0;
            } else {
                double log10 = log(10);
                double sub = alignmentPtrs[*partIt]->probScore - mapQVDenominator;
                double expo = exp(log10 * sub);
                double diff = 1.0 - expo;
                int phredValue;

                if (expo == 0) {
                    phredValue = 0;
                } else if (diff == 0) {
                    phredValue = MAX_PHRED_SCORE;
                } else {
                    phredValue = Phred(diff);
                }
                if (phredValue > MAX_PHRED_SCORE) {
                    phredValue = MAX_PHRED_SCORE;
                }

                alignmentPtrs[*partIt]->mapQV = phredValue;
                assigned[*partIt] = true;
            }

            if (params.scaleMapQVByNumSignificantClusters) {
                ScaleMapQVByClusterSize(*alignmentPtrs[*partIt], params);
            }
        }
    }

    for (i = 0; i < assigned.size(); i++) {
        assert(assigned[i]);
    }
}

//--------------------SEARCH & CHECK ALIGNMENTS-------------------//
template <typename T_Sequence>
bool CheckForSufficientMatch(T_Sequence &read, std::vector<T_AlignmentCandidate *> &alignmentPtrs,
                             MappingParameters &params)
{
    (void)(read);
    if (alignmentPtrs.size() > 0 and alignmentPtrs[0]->score < params.maxScore) {
        return true;
    } else {
        return false;
    }
}

int FindMaxLengthAlignment(std::vector<T_AlignmentCandidate *> alignmentPtrs, int &maxLengthIndex)
{
    int i;
    int maxLength = 0;
    maxLengthIndex = -1;

    for (i = 0; i < int(alignmentPtrs.size()); i++) {
        int qStart, qEnd;
        alignmentPtrs[i]->GetQInterval(qStart, qEnd);
        if (qEnd - qStart > maxLength) {
            maxLengthIndex = i;
            maxLength = qEnd - qStart;
        }
    }
    return (maxLength != -1);
}

void SumMismatches(SMRTSequence &read, T_AlignmentCandidate &alignment, int mismatchScore,
                   int fullIntvStart, int fullIntvEnd, MappingParameters &params, int &sum)
{
    int alnStart, alnEnd;
    alignment.GetQIntervalOnForwardStrand(alnStart, alnEnd);
    int p;
    sum = 0;
    if (not params.ignoreQualities and read.substitutionQV.Empty() == false) {
        for (p = fullIntvStart; p < alnStart; p++) {
            sum += read.substitutionQV[p];
        }
        for (p = alnEnd; p < fullIntvEnd; p++) {
            sum += read.substitutionQV[p];
        }
    } else {
        // bug 24363, compute mismatch score when QV is not available.
        sum += mismatchScore * ((alnStart - fullIntvStart) + (fullIntvEnd - alnEnd));
    }
}

bool AlignmentsOverlap(T_AlignmentCandidate &alnA, T_AlignmentCandidate &alnB,
                       float minPercentOverlap)
{
    int alnAStart, alnAEnd, alnBStart, alnBEnd;
    bool useForwardStrand = true;
    alnA.GetQInterval(alnAStart, alnAEnd, useForwardStrand);
    alnB.GetQInterval(alnBStart, alnBEnd, useForwardStrand);
    // Look if one alignment encompasses the other
    int ovp = 0;
    if (alnAStart <= alnBStart and alnAEnd >= alnBEnd) {
        return true;
    } else if (alnBStart <= alnAStart and alnBEnd >= alnAEnd) {
        return true;
        //ovp = alnAEnd - alnAStart;
    } else {
        //
        // Look to see if the alignments overlap
        //

        if (alnAEnd >= alnBStart and alnAEnd <= alnBEnd) {
            ovp = alnAEnd - alnBStart;
        } else if (alnAStart >= alnBStart and alnAStart <= alnBEnd) {
            ovp = alnBEnd - alnAStart;
        }
    }

    // float ovpPercent = (2.0*ovp) / ((alnAEnd - alnAStart) + (alnBEnd - alnBStart));
    float ovpPercent = 0;
    if (alnAEnd - alnAStart > 0 and alnBEnd - alnBStart > 0) {
        // overlap percentage: maximum overlap percent in A and B.
        ovpPercent = std::max(float(ovp) / float(alnAEnd - alnAStart),
                              float(ovp) / float(alnBEnd - alnBStart));
    }

    // returns true when an overlap is found.
    return (ovpPercent > minPercentOverlap);
}

void PartitionOverlappingAlignments(std::vector<T_AlignmentCandidate *> &alignmentPtrs,
                                    std::vector<std::set<int> > &partitions, float minOverlap)
{
    if (alignmentPtrs.size() == 0) {
        partitions.clear();
        return;
    }

    std::set<int>::iterator setIt, setEnd;
    int i, p;
    bool overlapFound = false;
    for (i = 0; i < int(alignmentPtrs.size()); i++) {
        overlapFound = false;
        for (p = 0; p < int(partitions.size()) and overlapFound == false; p++) {
            setEnd = partitions[p].end();
            for (setIt = partitions[p].begin();
                 setIt != partitions[p].end() and overlapFound == false; ++setIt) {
                if (AlignmentsOverlap(*alignmentPtrs[i], *alignmentPtrs[*setIt], minOverlap) or
                    ((alignmentPtrs[i]->QAlignStart() <= alignmentPtrs[*setIt]->QAlignStart()) and
                     (alignmentPtrs[i]->QAlignEnd() > alignmentPtrs[*setIt]->QAlignEnd()))) {
                    partitions[p].insert(i);
                    overlapFound = true;
                }
            }
        }
        //
        // If this alignment does not overlap any other, create a
        // partition with it as the first element.
        //
        if (overlapFound == false) {
            partitions.push_back(std::set<int>());
            partitions[partitions.size() - 1].insert(i);
        }
    }
}

//--------------------FILTER ALIGNMENTS---------------------------//
int RemoveLowQualitySDPAlignments(int readLength,
                                  std::vector<T_AlignmentCandidate *> &alignmentPtrs,
                                  MappingParameters &params)
{
    // Just a hack.  For now, assume there is at least 1 match per 50 bases.
    int totalBasesMatched = 0;
    int a;
    for (a = 0; a < int(alignmentPtrs.size()); a++) {
        int b;
        for (b = 0; b < int(alignmentPtrs[a]->blocks.size()); b++) {
            totalBasesMatched += alignmentPtrs[a]->blocks[b].length;
        }
        int expectedMatches = params.sdpTupleSize / 50.0 * readLength;
        if (totalBasesMatched < expectedMatches) {
            delete alignmentPtrs[a];
            alignmentPtrs[a] = NULL;
        }
    }
    int packedAlignmentIndex = 0;
    for (a = 0; a < int(alignmentPtrs.size()); a++) {
        if (alignmentPtrs[a] != NULL) {
            alignmentPtrs[packedAlignmentIndex] = alignmentPtrs[a];
            packedAlignmentIndex++;
        }
    }
    alignmentPtrs.resize(packedAlignmentIndex);
    return packedAlignmentIndex;
}

template <typename T_Sequence>
int RemoveLowQualityAlignments(T_Sequence &read, std::vector<T_AlignmentCandidate *> &alignmentPtrs,
                               MappingParameters &params)
{
    PB_UNUSED(read);
    if (params.verbosity > 0) {
        std::cout << "checking at least " << alignmentPtrs.size()
                  << " alignments to see if they are accurate." << std::endl;
    }
    for (size_t i = 0; i < MIN(static_cast<size_t>(params.nCandidates), alignmentPtrs.size());
         i++) {
        if (params.verbosity > 0) {
            std::cout << "Quality check  " << i << " " << alignmentPtrs[i]->score << std::endl;
        }
        if (alignmentPtrs[i]->blocks.size() == 0 or alignmentPtrs[i]->score > params.maxScore) {
            //
            // Since the alignments are sorted according to alignment
            // score, once one of the alignments is too low of a score,
            // all remaining alignments are also too low, and should be
            // removed as well.  Do that all at once.
            //
            if (alignmentPtrs[i]->blocks.size() == 0 and params.verbosity > 0) {
                std::cout << "Removing empty alignment " << alignmentPtrs[i]->qName << std::endl;
            }
            if (params.verbosity > 0) {
                std::cout << alignmentPtrs[i]->qName << " alignment " << i
                          << " is too low of a score." << alignmentPtrs[i]->score << std::endl;
            }
            for (size_t deletedIndex = i; deletedIndex < alignmentPtrs.size(); deletedIndex++) {
                delete alignmentPtrs[deletedIndex];
                alignmentPtrs[deletedIndex] = NULL;
            }
            alignmentPtrs.erase(i + alignmentPtrs.begin(), alignmentPtrs.end());
            break;
        } else {
            if (params.verbosity > 0) {
                std::cout << "Keeping alignment " << i << " " << alignmentPtrs[i]->qPos << " "
                          << alignmentPtrs[i]->qLength << " " << alignmentPtrs[i]->tName << " "
                          << alignmentPtrs[i]->tPos << " " << alignmentPtrs[i]->tLength
                          << " from score: " << alignmentPtrs[i]->score << std::endl;
            }
        }
    }
    return alignmentPtrs.size();
}

//FIXME: move to class ReadAlignments
int RemoveOverlappingAlignments(std::vector<T_AlignmentCandidate *> &alignmentPtrs,
                                MappingParameters &params)
{
    std::vector<unsigned char> alignmentIsContained;
    alignmentIsContained.resize(alignmentPtrs.size());
    std::fill(alignmentIsContained.begin(), alignmentIsContained.end(), false);

    int j;
    int numContained = 0;
    int curNotContained = 0;

    if (alignmentPtrs.size() > 0) {
        UInt i;
        for (i = 0; i < alignmentPtrs.size() - 1; i++) {
            T_AlignmentCandidate *aref = alignmentPtrs[i];
            if (aref->pctSimilarity < params.minPctSimilarity) {
                continue;
            }
            for (j = i + 1; j < int(alignmentPtrs.size()); j++) {
                //
                // Make sure this alignment isn't already removed.
                //
                if (alignmentIsContained[j]) {
                    continue;
                }

                //
                // Only check for containment if the two sequences are from the same contig.
                //
                if (alignmentPtrs[i]->tIndex != alignmentPtrs[j]->tIndex) {
                    continue;
                }

                //
                // Check for an alignment that is fully overlapping another
                // alignment.
                if (aref->GenomicTBegin() <= alignmentPtrs[j]->GenomicTBegin() and
                    aref->GenomicTEnd() >= alignmentPtrs[j]->GenomicTEnd() and
                    alignmentPtrs[i]->tIndex == alignmentPtrs[j]->tIndex) {
                    //
                    // Alignment i is contained in j is only true if it has a worse score.
                    //
                    if (aref->score <= alignmentPtrs[j]->score) {
                        alignmentIsContained[j] = true;
                    }
                    if (params.verbosity >= 2) {
                        std::cout << "alignment " << i << " is contained in " << j << std::endl;
                        std::cout << aref->tAlignedSeqPos << " " << alignmentPtrs[j]->tAlignedSeqPos
                                  << " " << aref->tAlignedSeqPos + aref->tAlignedSeqLength << " "
                                  << alignmentPtrs[j]->tAlignedSeqPos +
                                         alignmentPtrs[j]->tAlignedSeqLength
                                  << std::endl;
                    }
                } else if (alignmentPtrs[j]->GenomicTBegin() <= aref->GenomicTBegin() and
                           alignmentPtrs[j]->GenomicTEnd() >= aref->GenomicTEnd() and
                           alignmentPtrs[i]->tIndex == alignmentPtrs[j]->tIndex) {
                    if (params.verbosity >= 2) {
                        std::cout << "ALIGNMENT " << j << " is contained in " << i << std::endl;
                        std::cout << alignmentPtrs[j]->tAlignedSeqPos << " " << aref->tAlignedSeqPos
                                  << " "
                                  << alignmentPtrs[j]->tAlignedSeqPos +
                                         alignmentPtrs[j]->tAlignedSeqLength
                                  << " " << aref->tAlignedSeqPos + aref->tAlignedSeqLength
                                  << std::endl;
                    }
                    if (alignmentPtrs[j]->score <= aref->score) {
                        alignmentIsContained[i] = true;
                    }
                }
            }
        }
        for (i = 0; i < alignmentPtrs.size(); i++) {
            T_AlignmentCandidate *aref = alignmentPtrs[i];
            if (alignmentIsContained[i]) {
                delete alignmentPtrs[i];
                alignmentPtrs[i] = NULL;
                numContained++;
            } else {
                alignmentPtrs[curNotContained] = aref;
                ++curNotContained;
            }
        }
        alignmentPtrs.resize(alignmentPtrs.size() - numContained);
    }
    return alignmentPtrs.size();
}

// Delete all alignments from index startIndex in vector, inclusive.
void DeleteAlignments(std::vector<T_AlignmentCandidate *> &alignmentPtrs, int startIndex)
{
    int i;
    for (i = startIndex; i < int(alignmentPtrs.size()); i++) {
        delete alignmentPtrs[i];
    }
    alignmentPtrs.resize(0);
}

//--------------------REFINE ALIGNMENTS---------------------------//
template <typename T_RefSequence, typename T_Sequence>
void RefineAlignment(std::vector<T_Sequence *> &bothQueryStrands, T_RefSequence &genome,
                     T_AlignmentCandidate &alignmentCandidate, MappingParameters &params,
                     MappingBuffers &mappingBuffers)
{
    (void)(genome);
    FASTQSequence qSeq;
    DNASequence tSeq;
    DistanceMatrixScoreFunction<DNASequence, FASTQSequence> distScoreFn(
        SMRTDistanceMatrix, params.deletion, params.insertion);

    DistanceMatrixScoreFunction<DNASequence, FASTQSequence> distScoreFn2(
        SMRTDistanceMatrix, params.indel, params.indel);

    QualityValueScoreFunction<DNASequence, FASTQSequence> scoreFn;
    IDSScoreFunction<DNASequence, FASTQSequence> idsScoreFn;
    idsScoreFn.InitializeScoreMatrix(SMRTDistanceMatrix);
    scoreFn.del = params.indel;
    scoreFn.ins = params.indel;
    idsScoreFn.ins = params.insertion;
    idsScoreFn.del = params.deletion;
    idsScoreFn.affineExtend = params.affineExtend;
    idsScoreFn.affineOpen = params.affineOpen;
    idsScoreFn.substitutionPrior = params.substitutionPrior;
    idsScoreFn.globalDeletionPrior = params.globalDeletionPrior;

    if (params.doGlobalAlignment) {
        // global and placeGapConsistently can not set at the same time
        assert(not params.placeGapConsistently);
        SMRTSequence subread;
        subread.ReferenceSubstring(*bothQueryStrands[0], bothQueryStrands[0]->SubreadStart(),
                                   (bothQueryStrands[0]->SubreadLength()));

        int drift = ComputeDrift(alignmentCandidate);
        T_AlignmentCandidate refinedAlignment;

        KBandAlign(subread, alignmentCandidate.tAlignedSeq, SMRTDistanceMatrix, params.insertion,
                   params.deletion, drift, mappingBuffers.scoreMat, mappingBuffers.pathMat,
                   refinedAlignment, idsScoreFn, Global);
        refinedAlignment.RemoveEndGaps();
        ComputeAlignmentStats(refinedAlignment, subread.seq, alignmentCandidate.tAlignedSeq.seq,
                              distScoreFn2);
        //idsScoreFn);

        alignmentCandidate.blocks = refinedAlignment.blocks;
        alignmentCandidate.gaps = refinedAlignment.gaps;
        alignmentCandidate.tPos = refinedAlignment.tPos;
        alignmentCandidate.qPos = refinedAlignment.qPos + bothQueryStrands[0]->SubreadStart();
        alignmentCandidate.score = refinedAlignment.score;
        subread.Free();
    } else if (params.useGuidedAlign) {
        T_AlignmentCandidate refinedAlignment;
        int lastBlock = alignmentCandidate.blocks.size() - 1;

        if (alignmentCandidate.blocks.size() > 0) {

            /*
             * Refine the alignment without expanding past the current
             * boundaries of the sequences that are already aligned.
             */

            //
            // NOTE** this only makes sense when
            // alignmentCandidate.blocks[0].tPos == 0. Otherwise the length
            // of the sequence is not correct.
            //
            tSeq.Copy(
                alignmentCandidate.tAlignedSeq, alignmentCandidate.tPos,
                (alignmentCandidate.blocks[lastBlock].tPos +
                 alignmentCandidate.blocks[lastBlock].length - alignmentCandidate.blocks[0].tPos));

            //      qSeq.ReferenceSubstring(alignmentCandidate.qAlignedSeq,
            qSeq.ReferenceSubstring(*bothQueryStrands[static_cast<size_t>(
                                        (alignmentCandidate.qStrand == Forward) ? 0 : 1)],
                                    alignmentCandidate.qAlignedSeqPos + alignmentCandidate.qPos,
                                    (alignmentCandidate.blocks[lastBlock].qPos +
                                     alignmentCandidate.blocks[lastBlock].length));

            assert(not(params.affineAlign and params.placeGapConsistently));

            if (!params.ignoreQualities &&
                ReadHasMeaningfulQualityValues(alignmentCandidate.qAlignedSeq)) {
                if (params.affineAlign) {
                    AffineGuidedAlign(qSeq, tSeq, alignmentCandidate, idsScoreFn, params.bandSize,
                                      mappingBuffers, refinedAlignment, Global, false);
                } else {
                    GuidedAlign(qSeq, tSeq, alignmentCandidate, idsScoreFn,
                                params.guidedAlignBandSize, mappingBuffers, refinedAlignment,
                                Global, false);
                }
            } else {
                if (params.affineAlign) {
                    AffineGuidedAlign(qSeq, tSeq, alignmentCandidate, distScoreFn, params.bandSize,
                                      mappingBuffers, refinedAlignment, Global, false);
                } else {
                    GuidedAlign(qSeq, tSeq, alignmentCandidate, distScoreFn,
                                params.guidedAlignBandSize, mappingBuffers, refinedAlignment,
                                Global, false);
                }
            }
            ComputeAlignmentStats(refinedAlignment, qSeq.seq, tSeq.seq, distScoreFn2,
                                  params.affineAlign);
            //
            // Copy the refine alignment, which may be a subsequence of the
            // alignmentCandidate into the alignment candidate.
            //

            // First copy the alignment block and gap (the description of
            // the base by base alignment).

            alignmentCandidate.blocks.clear();
            alignmentCandidate.blocks = refinedAlignment.blocks;

            alignmentCandidate.CopyStats(refinedAlignment);

            alignmentCandidate.gaps = refinedAlignment.gaps;
            alignmentCandidate.score = refinedAlignment.score;
            alignmentCandidate.nCells = refinedAlignment.nCells;

            // Next copy the information that describes what interval was
            // aligned.  Since the reference sequences of the alignment
            // candidate have been modified, they are reassigned.
            alignmentCandidate.tAlignedSeq.Free();
            alignmentCandidate.tAlignedSeq.TakeOwnership(tSeq);
            alignmentCandidate.ReassignQSequence(qSeq);
            alignmentCandidate.tAlignedSeqPos += alignmentCandidate.tPos;
            alignmentCandidate.qAlignedSeqPos += alignmentCandidate.qPos;

            //
            // tPos and qPos are the positions within the interval where the
            // alignment begins. The refined alignment has adifferent tPos
            // and qPos from the alignment candidate.
            alignmentCandidate.tPos = refinedAlignment.tPos;
            alignmentCandidate.qPos = refinedAlignment.qPos;

            // The lengths of the newly aligned sequences may differ, update those.
            alignmentCandidate.tAlignedSeqLength = tSeq.length;
            alignmentCandidate.qAlignedSeqLength = qSeq.length;
        }
    } else {
        assert(not params.placeGapConsistently);

        //
        // This assumes an SDP alignment has been performed to create 'alignmentCandidate'.

        //
        // Recompute the alignment using a banded smith waterman to
        // get rid of any spurious effects of usign the seeded gaps.
        //

        //
        // The k-banded alignment is over a subsequence of the first
        // (sparse dynamic programming, SDP) alignment.  The SDP
        // alignment is over a large window that may contain the
        // candidate sequence.  The k-band alignment is over a tighter
        // region.

        int drift = ComputeDrift(alignmentCandidate);

        //
        // Rescore the alignment with a banded alignment that has a
        // better model of sequencing error.
        //

        if (alignmentCandidate.blocks.size() == 0) {
            alignmentCandidate.score = 0;
            return;
        }
        int lastBlock = alignmentCandidate.blocks.size() - 1;

        //
        // Assign the sequences that are going to be realigned using
        // banded alignment.  The SDP alignment does not give that great
        // of a score, but it does do a good job at finding a backbone
        // alignment that closely defines the sequence that is aligned.
        // Reassign the subsequences for alignment with a tight bound
        // around the beginning and ending of each sequence, so that
        // global banded alignment may be performed.
        //

        //
        // This section needs to be cleaned up substantially.  Right now it
        // copies a substring from the ref to a temp, then from the temp
        // back to the ref.  It may be possible to just keep one pointer per
        // read to the memory that was allocated, then allow the seq
        // parameter to float around.  The reason for all the copying is
        // that in case there is a compressed version of the genome the
        // seqences must be transformed before alignment.
        //

        if (alignmentCandidate.qIsSubstring) {
            qSeq.ReferenceSubstring(*bothQueryStrands[0],  // the original sequence
                                    alignmentCandidate.qPos + alignmentCandidate.qAlignedSeqPos,
                                    alignmentCandidate.blocks[lastBlock].qPos +
                                        alignmentCandidate.blocks[lastBlock].length);
        } else {
            qSeq.ReferenceSubstring(
                alignmentCandidate.qAlignedSeq,  // the subsequence that the alignment points to
                alignmentCandidate.qPos + alignmentCandidate.qAlignedSeqPos,
                alignmentCandidate.blocks[lastBlock].qPos +
                    alignmentCandidate.blocks[lastBlock].length -
                    alignmentCandidate.blocks[0].qPos);
        }

        tSeq.Copy(alignmentCandidate.tAlignedSeq,  // the subsequence the alignment points to
                  alignmentCandidate.tPos,         // ofset into the subsequence
                  alignmentCandidate.blocks[lastBlock].tPos +
                      alignmentCandidate.blocks[lastBlock].length -
                      alignmentCandidate.blocks[0].tPos);

        T_AlignmentCandidate refinedAlignment;

        //
        // When the parameter bandSize is 0, set the alignment band size
        // to the drift off the diagonal, plus a little more for wiggle
        // room.  When the parameteris nonzero, use that as a fixed band.
        //
        int k;
        if (params.bandSize == 0) {
            k = std::abs(drift) * 1.5;
        } else {
            k = params.bandSize;
        }
        if (params.verbosity > 0) {
            std::cout << "drift: " << drift << " qlen: " << alignmentCandidate.qAlignedSeq.length
                      << " tlen: " << alignmentCandidate.tAlignedSeq.length << " k: " << k
                      << std::endl;
            std::cout << "aligning in " << k << " * " << alignmentCandidate.tAlignedSeq.length
                      << " " << k * alignmentCandidate.tAlignedSeq.length << std::endl;
        }
        if (k < 10) {
            k = 10;
        }

        alignmentCandidate.tAlignedSeqPos += alignmentCandidate.tPos;

        VectorIndex lastSDPBlock = alignmentCandidate.blocks.size() - 1;

        if (alignmentCandidate.blocks.size() > 0) {
            alignmentCandidate.tAlignedSeqLength = (alignmentCandidate.blocks[lastSDPBlock].tPos +
                                                    alignmentCandidate.blocks[lastSDPBlock].length -
                                                    alignmentCandidate.blocks[0].tPos);
        } else {
            alignmentCandidate.tAlignedSeqLength = 0;
        }

        alignmentCandidate.tPos = 0;
        alignmentCandidate.qAlignedSeqPos += alignmentCandidate.qPos;

        if (alignmentCandidate.blocks.size() > 0) {
            alignmentCandidate.qAlignedSeqLength = (alignmentCandidate.blocks[lastSDPBlock].qPos +
                                                    alignmentCandidate.blocks[lastSDPBlock].length -
                                                    alignmentCandidate.blocks[0].qPos);
        } else {
            alignmentCandidate.qAlignedSeqLength = 0;
        }
        alignmentCandidate.qPos = 0;

        alignmentCandidate.blocks.clear();
        alignmentCandidate.tAlignedSeq.Free();
        alignmentCandidate.tAlignedSeq.TakeOwnership(tSeq);
        alignmentCandidate.ReassignQSequence(qSeq);

        if (params.verbosity >= 2) {
            std::cout << "refining target: " << std::endl;
            alignmentCandidate.tAlignedSeq.PrintSeq(std::cout);
            std::cout << "refining query: " << std::endl;
            static_cast<DNASequence *>(&alignmentCandidate.qAlignedSeq)->PrintSeq(std::cout);
            std::cout << std::endl;
        }
        PairwiseLocalAlign(qSeq, tSeq, k, params, alignmentCandidate, mappingBuffers, Fit);
    }
}

template <typename T_RefSequence, typename T_Sequence>
void RefineAlignments(std::vector<T_Sequence *> &bothQueryStrands, T_RefSequence &genome,
                      std::vector<T_AlignmentCandidate *> &alignmentPtrs, MappingParameters &params,
                      MappingBuffers &mappingBuffers)
{

    UInt i;
    for (i = 0; i < alignmentPtrs.size(); i++) {
        RefineAlignment(bothQueryStrands, genome, *alignmentPtrs[i], params, mappingBuffers);
    }
    //
    // It's possible the alignment references change their order after running
    // the local alignments.  This is made into a parameter rather than resorting
    // every time so that the performance gain by resorting may be measured.
    //
    if (params.sortRefinedAlignments) {
        std::sort(alignmentPtrs.begin(), alignmentPtrs.end(), SortAlignmentPointersByScore());
    }
}

std::vector<T_AlignmentCandidate *> SelectAlignmentsToPrint(
    std::vector<T_AlignmentCandidate *> alignmentPtrs, MappingParameters &params,
    const int &associatedRandInt)
{
    if (params.placeRandomly) {
        assert(params.hitPolicy.IsRandombest());
    }

    if (alignmentPtrs.size() == 0) {
        return std::vector<T_AlignmentCandidate *>({});
    }

    std::sort(alignmentPtrs.begin(), alignmentPtrs.end(), SortAlignmentPointersByScore());

    // Apply filter criteria and hit policy.
    // Shallow copy AlignmentCandidate pointers.
    std::vector<T_AlignmentCandidate *> filtered;
    for (auto ptr : alignmentPtrs) {
        if (params.filterCriteria.Satisfy(ptr)) {
            filtered.push_back(ptr);
            if (int(filtered.size()) == params.nBest) break;
        }
    }

    return params.hitPolicy.Apply(filtered, false, associatedRandInt);
}

// The full read is not the subread, and does not have masked off characters.
void PrintAlignment(T_AlignmentCandidate &alignment, SMRTSequence &fullRead,
                    MappingParameters &params, AlignmentContext &alignmentContext,
                    std::ostream &outFile
#ifdef USE_PBBAM
                    ,
                    SMRTSequence &subread, PacBio::BAM::IRecordWriter *bamWriterPtr
#endif
                    )
{
    try {
        // Before printing alignments, make sure query is always Forward.
        // If query was Reverse, reverse complement both query and target
        // and recomputes all coordinates.
        alignment.MakeQueryForward();

        if (params.printFormat == StickPrint) {
            PrintAlignmentStats(alignment, outFile);
            StickPrintAlignment(alignment, (DNASequence &)alignment.qAlignedSeq,
                                (DNASequence &)alignment.tAlignedSeq, outFile,
                                alignment.qAlignedSeqPos, alignment.tAlignedSeqPos);
        } else if (params.printFormat == SAM) {
            SAMOutput::PrintAlignment(alignment, fullRead, outFile, alignmentContext,
                                      params.samQVList, params.clipping, params.cigarUseSeqMatch,
                                      params.allowAdjacentIndels);
        } else if (params.printFormat == BAM) {
#ifdef USE_PBBAM
            BAMOutput::PrintAlignment(alignment, fullRead, subread, *bamWriterPtr, alignmentContext,
                                      params.samQVList, params.clipping, params.cigarUseSeqMatch,
                                      params.allowAdjacentIndels);
#else
            REQUIRE_PBBAM_ERROR();
#endif
        } else if (params.printFormat == CompareXML) {
            XMLOutput::Print(alignment, (DNASequence &)alignment.qAlignedSeq,
                             (DNASequence &)alignment.tAlignedSeq, outFile,
                             alignment.qAlignedSeqPos, alignment.tAlignedSeqPos);
        } else if (params.printFormat == Vulgar) {
            PrintAlignmentStats(alignment, outFile);
            VulgarOutput::Print(alignment, outFile);
        } else if (params.printFormat == CompareSequencesParsable) {
            CompareSequencesOutput::Print(alignment, alignment.qAlignedSeq, alignment.tAlignedSeq,
                                          outFile);
        } else if (params.printFormat == Interval) {
            if (alignment.blocks.size() > 0) {
                IntervalOutput::Print(alignment, outFile);
            }
        } else if (params.printFormat == SummaryPrint) {
            if (alignment.blocks.size() > 0) {
                SummaryOutput::Print(alignment, outFile);
            }
        }
    } catch (const std::ostream::failure &f) {
        std::cout << "ERROR writing to output file. The output drive may be full, or you  "
                  << std::endl;
        std::cout << "may not have proper write permissions." << std::endl;
        std::exit(EXIT_FAILURE);
    }
}

// Print all alignments in std::vector<T_AlignmentCandidate*> alignmentPtrs
void PrintAlignments(std::vector<T_AlignmentCandidate *> alignmentPtrs, SMRTSequence &read,
                     MappingParameters &params, std::ostream &outFile,
                     AlignmentContext alignmentContext,
#ifdef USE_PBBAM
                     SMRTSequence &subread, PacBio::BAM::IRecordWriter *bamWriterPtr,
#endif
                     MappingSemaphores &semaphores)
{
    if (params.nProc > 1) {
#ifdef __APPLE__
        sem_wait(semaphores.writer);
#else
        sem_wait(&semaphores.writer);
#endif
    }
    for (int i = 0; i < int(alignmentPtrs.size()); i++) {
        T_AlignmentCandidate *aref = alignmentPtrs[i];

        if (aref->blocks.size() == 0) {

            //
            // If the SDP alignment finds nothing, there will be no
            // blocks.  This may happen if the sdp block size is larger
            // than the anchor size found with the suffix array.  When no
            // blocks are found there is no alignment, so zero-out the
            // score and continue.
            //
            aref->score = 0;
            if (params.verbosity > 0) {
                std::cout << "Zero blocks found for " << aref->qName << " " << aref->qAlignedSeqPos
                          << " " << aref->tAlignedSeqPos << std::endl;
            }
            continue;
        }

        //
        // Configure some of the alignment context before printing.
        //
        if (i > 0 and params.placeRandomly == false) {
            alignmentContext.isPrimary = false;
        } else {
            alignmentContext.isPrimary = true;
        }

        if (params.printSAM or params.printBAM) {
            DistanceMatrixScoreFunction<DNASequence, FASTASequence> editdistScoreFn(
                EditDistanceMatrix, 1, 1);
            T_AlignmentCandidate &alignment = *alignmentPtrs[i];
            alignmentContext.editDist = ComputeAlignmentScore(
                alignment, alignment.qAlignedSeq, alignment.tAlignedSeq, editdistScoreFn);
        }

        PrintAlignment(*alignmentPtrs[i], read, params, alignmentContext, outFile
#ifdef USE_PBBAM
                       ,
                       subread, bamWriterPtr
#endif
                       );
    }

    if (params.nProc > 1) {
#ifdef __APPLE__
        sem_post(semaphores.writer);
#else
        sem_post(&semaphores.writer);
#endif
    }
}

void PrintAlignmentPtrs(std::vector<T_AlignmentCandidate *> &alignmentPtrs, std::ostream &out)
{
    for (int alignmentIndex = 0; alignmentIndex < int(alignmentPtrs.size()); alignmentIndex++) {
        out << "[" << alignmentIndex << "/" << alignmentPtrs.size() << "]" << std::endl;
        T_AlignmentCandidate *alignment = alignmentPtrs[alignmentIndex];
        alignment->Print(out);
    }
    out << std::endl;
}

void PrintUnaligned(const SMRTSequence &unalignedRead, std::ostream &unalignedFilePtr,
                    const bool noPrintUnalignedSeqs)
{
    if (noPrintUnalignedSeqs) {
        std::string s = unalignedRead.GetTitle();
        SMRTTitle st(s);
        if (st.isSMRTTitle)
            unalignedFilePtr << st.ToString() << std::endl;
        else
            //size_t pos = s.rfind("/");
            //if (pos != string::npos)
            //    unalignedFilePtr << s.substr(0, pos) << std::endl;
            //else
            unalignedFilePtr << s << std::endl;
    } else
        unalignedRead.PrintSeq(unalignedFilePtr);
}

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
void PrintAllReadAlignments(ReadAlignments &allReadAlignments, AlignmentContext &alignmentContext,
                            std::ostream &outFilePtr, std::ostream &unalignedFilePtr,
                            MappingParameters &params, std::vector<SMRTSequence> &subreads,
#ifdef USE_PBBAM
                            PacBio::BAM::IRecordWriter *bamWriterPtr,
#endif
                            MappingSemaphores &semaphores)
{
    int subreadIndex;
    int nAlignedSubreads = allReadAlignments.GetNAlignedSeq();

    //
    // Initialize the alignemnt context with information applicable to SAM output.
    //
    alignmentContext.alignMode = allReadAlignments.alignMode;
    for (subreadIndex = 0; subreadIndex < nAlignedSubreads; subreadIndex++) {
        if (allReadAlignments.subreadAlignments[subreadIndex].size() > 0) {
            alignmentContext.numProperlyAlignedSubreads++;
        }
    }

    if (alignmentContext.numProperlyAlignedSubreads ==
        int(allReadAlignments.subreadAlignments.size())) {
        alignmentContext.allSubreadsProperlyAligned = true;
    }
    alignmentContext.nSubreads = nAlignedSubreads;

    for (subreadIndex = 0; subreadIndex < nAlignedSubreads; subreadIndex++) {
        alignmentContext.subreadIndex = subreadIndex;
        if (subreadIndex < nAlignedSubreads - 1 and
            allReadAlignments.subreadAlignments[subreadIndex + 1].size() > 0) {
            alignmentContext.nextSubreadPos =
                allReadAlignments.subreadAlignments[subreadIndex + 1][0]->QAlignStart();
            alignmentContext.nextSubreadDir =
                allReadAlignments.subreadAlignments[subreadIndex + 1][0]->qStrand;
            alignmentContext.rNext =
                allReadAlignments.subreadAlignments[subreadIndex + 1][0]->tName;
            alignmentContext.hasNextSubreadPos = true;
        } else {
            alignmentContext.nextSubreadPos = 0;
            alignmentContext.nextSubreadDir = 0;
            alignmentContext.rNext = "";
            alignmentContext.hasNextSubreadPos = false;
        }
        SMRTSequence *sourceSubread = &(allReadAlignments.subreads[subreadIndex]);
        if (subreads.size() == allReadAlignments.subreads.size()) {
            sourceSubread = &subreads[subreadIndex];
        }
        if (allReadAlignments.subreadAlignments[subreadIndex].size() > 0) {
            PrintAlignments(allReadAlignments.subreadAlignments[subreadIndex],
                            allReadAlignments.subreads[subreadIndex],
                            // for these alignments
                            params, outFilePtr,  //*mapData->outFilePtr,
                            alignmentContext,
#ifdef USE_PBBAM
                            *sourceSubread, bamWriterPtr,
#endif
                            semaphores);
        } else {
            //
            // Print the unaligned sequences.
            //
            if (params.printUnaligned == true) {
                if (params.nProc == 1) {
                    PrintUnaligned(*sourceSubread, unalignedFilePtr, params.noPrintUnalignedSeqs);
                } else {
#ifdef __APPLE__
                    sem_wait(semaphores.unaligned);
#else
                    sem_wait(&semaphores.unaligned);
#endif
                    PrintUnaligned(*sourceSubread,  //subreads[subreadIndex],
                                   unalignedFilePtr, params.noPrintUnalignedSeqs);
#ifdef __APPLE__
                    sem_post(semaphores.unaligned);
#else
                    sem_post(&semaphores.unaligned);
#endif
                }  // End of nproc > 1.
            }      // End of printing  unaligned sequences.
        }          // End of finding no alignments for the subread with subreadIndex.
    }              // End of printing and processing alignmentContext for each subread.
}
