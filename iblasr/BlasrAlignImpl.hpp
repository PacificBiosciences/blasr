// Author: Mark Chaisson
#pragma once

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
             MappingSemaphores & semaphores)
{

    bool matchFound;
    WeightedIntervalSet topIntervals(params.nCandidates);
    int numKeysMatched=0, rcNumKeysMatched=0;
    int expand = params.minExpand;
    metrics.clocks.total.Tick();
    int nTotalCells = 0;
    int forwardNumBasesMatched = 0, reverseNumBasesMatched = 0;
    do {
        matchFound = false;
        mappingBuffers.matchPosList.clear();
        mappingBuffers.rcMatchPosList.clear();
        alignmentPtrs.clear();
        topIntervals.clear();
        params.anchorParameters.expand = expand;

        metrics.clocks.mapToGenome.Tick();

        if (params.useSuffixArray) {
            params.anchorParameters.lcpBoundsOutPtr = mapData->lcpBoundsOutPtr;
            numKeysMatched   =
                MapReadToGenome(genome, sarray, read,   params.lookupTableLength, mappingBuffers.matchPosList,
                        params.anchorParameters);

            //
            // Only print values for the read in forward direction (and only
            // the first read).
            //
            mapData->lcpBoundsOutPtr = NULL;
            if (!params.forwardOnly) {
                rcNumKeysMatched =
                    MapReadToGenome(genome, sarray, readRC, params.lookupTableLength, mappingBuffers.rcMatchPosList,
                            params.anchorParameters);
            }
        }
        else if (params.useBwt){
            numKeysMatched   = MapReadToGenome(bwt, read, read.SubreadStart(), read.SubreadEnd(),
                    mappingBuffers.matchPosList, params.anchorParameters, forwardNumBasesMatched);
            if (!params.forwardOnly) {
                rcNumKeysMatched = MapReadToGenome(bwt, readRC, readRC.SubreadStart(), readRC.SubreadEnd(),
                        mappingBuffers.rcMatchPosList, params.anchorParameters, reverseNumBasesMatched);
            }
        }

        //
        // Look to see if only the anchors are printed.
        if (params.anchorFileName != "") {
            int i;
            if (params.nProc > 1) {
#ifdef __APPLE__
                sem_wait(semaphores.writer);
#else
                sem_wait(&semaphores.writer);
#endif
            }
            *mapData->anchorFilePtr << read.title << endl;
            for (i = 0; i < mappingBuffers.matchPosList.size(); i++) {
                *mapData->anchorFilePtr << mappingBuffers.matchPosList[i] << endl;
            }
            *mapData->anchorFilePtr << readRC.title << " (RC) " << endl;
            for (i = 0; i < mappingBuffers.rcMatchPosList.size(); i++) {
                *mapData->anchorFilePtr << mappingBuffers.rcMatchPosList[i] << endl;
            }

            if (params.nProc > 1) {
#ifdef __APPLE__
                sem_post(semaphores.writer);
#else
                sem_post(&semaphores.writer);
#endif
            }
        }

        metrics.totalAnchors += mappingBuffers.matchPosList.size() + mappingBuffers.rcMatchPosList.size();
        metrics.clocks.mapToGenome.Tock();

        metrics.clocks.sortMatchPosList.Tick();
        SortMatchPosList(mappingBuffers.matchPosList);
        SortMatchPosList(mappingBuffers.rcMatchPosList);
        metrics.clocks.sortMatchPosList.Tock();

        PValueWeightor lisPValue(read, genome, ct.tm, &ct);
        MultiplicityPValueWeightor lisPValueByWeight(genome);

        LISSumOfLogPWeightor<T_GenomeSequence,vector<ChainedMatchPos> > lisPValueByLogSum(genome);

        LISSizeWeightor<vector<ChainedMatchPos> > lisWeightFn;

        IntervalSearchParameters intervalSearchParameters;
        intervalSearchParameters.globalChainType = params.globalChainType;
        intervalSearchParameters.advanceHalf = params.advanceHalf;
        intervalSearchParameters.warp        = params.warp;
        intervalSearchParameters.fastMaxInterval = params.fastMaxInterval;
        intervalSearchParameters.aggressiveIntervalCut = params.aggressiveIntervalCut;
        intervalSearchParameters.verbosity = params.verbosity;

        //
        // If specified, only align a band from the anchors.
        //
        DNALength squareRefLength = read.length * 1.25 + params.limsAlign;
        if (params.limsAlign != 0) {
            int fi;
            for (fi = 0; fi < mappingBuffers.matchPosList.size(); fi++) {
                if (mappingBuffers.matchPosList[fi].t >= squareRefLength) { break; }
            }
            if (fi < mappingBuffers.matchPosList.size()) {
                mappingBuffers.matchPosList.resize(fi);
            }
        }

        metrics.clocks.findMaxIncreasingInterval.Tick();

        //
        // For now say that something that has a 50% chance of happening
        // by chance is too high of a p value. This is probably many times
        // the size.
        //
        intervalSearchParameters.maxPValue = log(0.5);
        intervalSearchParameters.aboveCategoryPValue = -300;
        VarianceAccumulator<float> accumPValue;
        VarianceAccumulator<float> accumWeight;
        VarianceAccumulator<float> accumNBases;

        mappingBuffers.clusterList.Clear();
        mappingBuffers.revStrandClusterList.Clear();

        //
        // Remove anchors that are fully encompassed by longer ones.  This
        // speeds up limstemplate a lot.
        //

        RemoveOverlappingAnchors(mappingBuffers.matchPosList);
        RemoveOverlappingAnchors(mappingBuffers.rcMatchPosList);

        if (params.pValueType == 0) {
            int original = mappingBuffers.matchPosList.size();

            int numMerged = 0;
            if (params.printDotPlots) {
                ofstream dotPlotOut;
                string dotPlotName = string(read.title) + ".anchors";
                CrucialOpen(dotPlotName, dotPlotOut, std::ios::out);
                int mp;
                for (mp = 0; mp < mappingBuffers.matchPosList.size(); mp++ ){
                    dotPlotOut << mappingBuffers.matchPosList[mp].q << " " << mappingBuffers.matchPosList[mp].t << " " << mappingBuffers.matchPosList[mp].l << " " << endl;
                }
                dotPlotOut.close();
            }
            /*
               This is an optimization that is being tested out that places a grid over the
               area where there are anchors, and then finds an increasing maximally weighted
               path through the grid.  The weight of a cell in the grid is the sum of the
               number of anchors in it.  All other anchors are to be removed.  This will likely
               only work for LIMSTemplate sequences, or other sequences with little structural
               variation.
               FindBand(mappingBuffers.matchPosList,
               refCopy, read, 100);
               */
            FindMaxIncreasingInterval(Forward,
                    mappingBuffers.matchPosList,
                    // allow for indels to stretch out the mapping of the read.
                    (DNALength) ((read.SubreadLength()) * (1 + params.indelRate)), params.nCandidates,
                    seqBoundary,
                    lisPValue,//lisPValue2,
                    lisWeightFn,
                    topIntervals, genome, read, intervalSearchParameters,
                    &mappingBuffers.globalChainEndpointBuffer,
                    mappingBuffers.clusterList,
                    accumPValue, accumWeight, accumNBases, read.title);
            // Uncomment when the version of the weight functor needs the sequence.

            mappingBuffers.clusterList.ResetCoordinates();

            FindMaxIncreasingInterval(Reverse, mappingBuffers.rcMatchPosList,
                    (DNALength) ((read.SubreadLength()) * (1 + params.indelRate)), params.nCandidates,
                    seqBoundary,
                    lisPValue,//lisPValue2
                    lisWeightFn,
                    topIntervals, genome, readRC, intervalSearchParameters,
                    &mappingBuffers.globalChainEndpointBuffer,
                    mappingBuffers.revStrandClusterList,
                    accumPValue, accumWeight, accumNBases, read.title);
        }
        else if (params.pValueType == 1) {
            FindMaxIncreasingInterval(Forward,
                    mappingBuffers.matchPosList,
                    // allow for indels to stretch out the mapping of the read.
                    (DNALength) ((read.SubreadLength()) * (1 + params.indelRate)), params.nCandidates,
                    seqBoundary,
                    lisPValueByWeight, // different from pvaltype == 2 and 0
                    lisWeightFn,
                    topIntervals, genome, read, intervalSearchParameters,
                    &mappingBuffers.globalChainEndpointBuffer,
                    mappingBuffers.clusterList,
                    accumPValue, accumWeight, accumNBases,
                    read.title);


            mappingBuffers.clusterList.ResetCoordinates();
            FindMaxIncreasingInterval(Reverse, mappingBuffers.rcMatchPosList,
                    (DNALength) ((read.SubreadLength()) * (1 + params.indelRate)), params.nCandidates,
                    seqBoundary,
                    lisPValueByWeight, // different from pvaltype == 2 and 0
                    lisWeightFn,
                    topIntervals, genome, readRC, intervalSearchParameters,
                    &mappingBuffers.globalChainEndpointBuffer,
                    mappingBuffers.revStrandClusterList,
                    accumPValue, accumWeight, accumNBases,
                    read.title);
        }
        else if (params.pValueType == 2) {
            FindMaxIncreasingInterval(Forward,
                    mappingBuffers.matchPosList,
                    // allow for indels to stretch out the mapping of the read.
                    (DNALength) ((read.SubreadLength()) * (1 + params.indelRate)), params.nCandidates,
                    seqBoundary,
                    lisPValueByLogSum, // different from pvaltype == 1 and 0
                    lisWeightFn,
                    topIntervals, genome, read, intervalSearchParameters,
                    &mappingBuffers.globalChainEndpointBuffer,
                    mappingBuffers.clusterList,
                    accumPValue, accumWeight, accumNBases,
                    read.title);

            mappingBuffers.clusterList.ResetCoordinates();
            FindMaxIncreasingInterval(Reverse, mappingBuffers.rcMatchPosList,
                    (DNALength) ((read.SubreadLength()) * (1 + params.indelRate)), params.nCandidates,
                    seqBoundary,
                    lisPValueByLogSum, // different from pvaltype == 1 and 0
                    lisWeightFn,
                    topIntervals, genome, readRC, intervalSearchParameters,
                    &mappingBuffers.globalChainEndpointBuffer,
                    mappingBuffers.revStrandClusterList,
                    accumPValue, accumWeight, accumNBases,
                    read.title);
        }

        mappingBuffers.clusterList.numBases.insert(mappingBuffers.clusterList.numBases.end(),
                mappingBuffers.revStrandClusterList.numBases.begin(),
                mappingBuffers.revStrandClusterList.numBases.end());

        mappingBuffers.clusterList.numAnchors.insert(mappingBuffers.clusterList.numAnchors.end(),
                mappingBuffers.revStrandClusterList.numAnchors.begin(),
                mappingBuffers.revStrandClusterList.numAnchors.end());

        metrics.clocks.findMaxIncreasingInterval.Tock();

        //
        // Print verbose output.
        //
        WeightedIntervalSet::iterator topIntIt, topIntEnd;
        topIntEnd = topIntervals.end();
        if (params.verbosity > 0) {
            int topintind = 0;
            cout << " intv: index start end qstart qend seq_boundary_start seq_boundary_end pvalue " << endl;
            for (topIntIt = topIntervals.begin();topIntIt != topIntEnd ; ++topIntIt) {
                cout << " intv: " << topintind << " " << (*topIntIt).start << " "
                    << (*topIntIt).end << " "
                    << (*topIntIt).qStart << " " << (*topIntIt).qEnd << " "
                    << seqBoundary((*topIntIt).start) << " " << seqBoundary((*topIntIt).end) << " "
                    << (*topIntIt).pValue << endl;
                if (params.verbosity > 2) {
                    int m;
                    for (m = 0; m < (*topIntIt).matches.size(); m++) {
                        cout << " (" << (*topIntIt).matches[m].q << ", " << (*topIntIt).matches[m].t << ", " << (*topIntIt).matches[m].l << ") ";
                    }
                    cout << endl;
                }
                ++topintind;
            }
        }

        //
        // Allocate candidate alignments on the stack.  Each interval is aligned.
        //
        alignmentPtrs.resize(topIntervals.size());
        UInt i;
        for (i = 0; i < alignmentPtrs.size(); i++ ) {
            alignmentPtrs[i] = new T_AlignmentCandidate;
        }
        metrics.clocks.alignIntervals.Tick();
        AlignIntervals( genome, read, readRC,
                topIntervals,
                SMRTDistanceMatrix,
                params.indel, params.indel,
                params.sdpTupleSize,
                params.useSeqDB, seqdb,
                alignmentPtrs,
                params,
                mappingBuffers,
                params.startRead );

        /*    cout << read.title << endl;
              for (i = 0; i < alignmentPtrs.size(); i++) {
              cout << alignmentPtrs[i]->clusterScore << " " << alignmentPtrs[i]->score << endl;
              }
              */
        StoreRankingStats(alignmentPtrs, accumPValue, accumWeight);

        std::sort(alignmentPtrs.begin(), alignmentPtrs.end(), SortAlignmentPointersByScore());
        metrics.clocks.alignIntervals.Tock();

        //
        // Evalutate the matches that are found for 'good enough'.
        //

        matchFound = CheckForSufficientMatch(read, alignmentPtrs, params);

        //
        // When no proper alignments are found, the loop will resume.
        // Delete all alignments because they are bad.
        //
        if (expand < params.maxExpand and matchFound == false) {
            DeleteAlignments(alignmentPtrs, 0);
        }

        //
        // Record some metrics that show how long this took to run per base.
        //

        if (alignmentPtrs.size() > 0) {
            metrics.RecordNumAlignedBases(read.length);
            metrics.RecordNumCells(alignmentPtrs[0]->nCells);
        }

        if (matchFound == true) {
            metrics.totalAnchorsForMappedReads += mappingBuffers.matchPosList.size() + mappingBuffers.rcMatchPosList.size();
        }
        ++expand;
    } while ( expand <= params.maxExpand and matchFound == false);
    metrics.clocks.total.Tock();
    UInt i;
    int totalCells = 0;
    for (i = 0; i< alignmentPtrs.size(); i++) {
        totalCells += alignmentPtrs[i]->nCells;
    }
    metrics.clocks.AddCells(totalCells);
    int totalBases = 0;
    for (i = 0; i < alignmentPtrs.size(); i++) {
        totalBases += alignmentPtrs[i]->qLength;
    }
    metrics.clocks.AddBases(totalBases);
    //
    //  Some of the alignments are to spurious regions. Delete the
    //  references that have too small of a score.
    //

    int effectiveReadLength = 0;
    for (i = 0; i< read.length; i++) {
        if (read.seq[i] != 'N') effectiveReadLength++;
    }
    if (params.sdpFilterType == 0) {
        RemoveLowQualityAlignments(read, alignmentPtrs, params);
    }
    else if (params.sdpFilterType == 1) {
        RemoveLowQualitySDPAlignments(effectiveReadLength, alignmentPtrs, params);
    }

    //
    // Now remove overlapping alignments.
    //

    vector<T_Sequence*> bothQueryStrands;
    bothQueryStrands.resize(2);
    bothQueryStrands[Forward] = &read;
    bothQueryStrands[Reverse] = &readRC;


    //
    // Possibly use banded dynamic programming to refine the columns
    // of an alignment and the alignment score.
    //
    if (params.refineAlignments) {
        RefineAlignments(bothQueryStrands, genome, alignmentPtrs, params, mappingBuffers);
        RemoveLowQualityAlignments(read,alignmentPtrs,params);
        RemoveOverlappingAlignments(alignmentPtrs, params);
    }

    //
    // Look to see if the number of anchors found for this read match
    // what is expected given the expected distribution of number of
    // anchors.
    //

    if (alignmentPtrs.size() > 0) {
        int clusterIndex;
        //
        // Compute some stats on the read.  For now this is fixed but will
        // be updated on the fly soon.
        //
        float meanAnchorBasesPerRead, sdAnchorBasesPerRead;
        float meanAnchorsPerRead, sdAnchorsPerRead;

        int lookupValue;
        //
        // If a very short anchor size was used, or very long min match
        // size there may be no precomputed distributions for it.
        // Handle this by bounding the min match by the smallest and
        // largest values for which there are precomputed statistics.

        int boundedMinWordMatchLength = min(max(params.minMatchLength, anchorMinKValues[0]), anchorMinKValues[1]);

        //
        // Do a similar bounding for match length and accuracy.
        //
        int boundedMatchLength  = min(max((int) alignmentPtrs[0]->qAlignedSeq.length, anchorReadLengths[0]), anchorReadLengths[1]);
        int boundedPctSimilarity = min(max((int)alignmentPtrs[0]->pctSimilarity, anchorReadAccuracies[0]), anchorReadAccuracies[1]);

        lookupValue = LookupAnchorDistribution(boundedMatchLength, boundedMinWordMatchLength, boundedPctSimilarity,
                meanAnchorsPerRead, sdAnchorsPerRead, meanAnchorBasesPerRead, sdAnchorBasesPerRead);

        float minExpAnchors = meanAnchorsPerRead -  sdAnchorsPerRead;
        //
        // The number of standard deviations is just trial and error.
        float minExpAnchorBases = meanAnchorBasesPerRead -  2 * sdAnchorBasesPerRead;
        if (lookupValue < 0 or minExpAnchorBases < 0) {
            minExpAnchorBases = 0;
        }
        int numSignificantClusters = 0;
        int totalSignificantClusterSize = 0;
        int maxClusterSize = 0;
        int maxClusterIndex = 0;
        int numAlnAnchorBases, numAlnAnchors, scaledMaxClusterSize;
        alignmentPtrs[0]->ComputeNumAnchors(boundedMinWordMatchLength, numAlnAnchors, numAlnAnchorBases);
        int totalAnchorBases = 0;
        if (numAlnAnchorBases > meanAnchorBasesPerRead + sdAnchorBasesPerRead) {
            numSignificantClusters = 1;
        }
        else {
            if (alignmentPtrs[0]->score < params.maxScore) {
                for (clusterIndex = 0; clusterIndex < mappingBuffers.clusterList.numBases.size(); clusterIndex++) {
                    if (mappingBuffers.clusterList.numBases[clusterIndex] > maxClusterSize) {
                        maxClusterSize = mappingBuffers.clusterList.numBases[clusterIndex];
                        maxClusterIndex = clusterIndex;
                    }
                }
                int scaledExpectedClusterSize = maxClusterSize / ((float)numAlnAnchorBases) * minExpAnchorBases;
                for (clusterIndex = 0; clusterIndex < mappingBuffers.clusterList.numBases.size(); clusterIndex++) {
                    bool isSignificant = false;
                    if (mappingBuffers.clusterList.numBases[clusterIndex] >= scaledExpectedClusterSize) {
                        //          cout << mappingBuffers.clusterList.numBases[clusterIndex] << " " << scaledExpectedClusterSize << " " << meanAnchorBasesPerRead << " " << sdAnchorBasesPerRead << endl;
                        ++numSignificantClusters;
                        totalSignificantClusterSize += meanAnchorBasesPerRead;
                        isSignificant = true;
                    }
                    //
                    // The following output block is useful in debugging mapqv
                    // calculation.   It should be uncommented and examined when
                    // mapqvs do not look correct.
                    //
                    totalAnchorBases +=  mappingBuffers.clusterList.numBases[clusterIndex];
                }
            }

            if (lookupValue == 0) {
                int scaledMaxClusterSize;
                alignmentPtrs[0]->ComputeNumAnchors(params.minMatchLength, numAlnAnchors, numAlnAnchorBases);
                scaledMaxClusterSize = (  ((float)numAlnAnchorBases )/ meanAnchorBasesPerRead) * maxClusterSize;
            }
        }

        for (i = 0; i < alignmentPtrs.size(); i++) {
            alignmentPtrs[i]->numSignificantClusters = numSignificantClusters;
        }
        if (mapData->clusterFilePtr != NULL and topIntervals.size() > 0 and alignmentPtrs.size() > 0) {
            WeightedIntervalSet::iterator intvIt = topIntervals.begin();
            if (params.nProc > 1) {
#ifdef __APPLE__
                sem_wait(semaphores.hitCluster);
#else
                sem_wait(&semaphores.hitCluster);
#endif
            }

            *mapData->clusterFilePtr << (*intvIt).size << " " << (*intvIt).pValue << " " << (*intvIt).nAnchors << " "
                << read.length << " " << alignmentPtrs[0]->score << " " << alignmentPtrs[0]->pctSimilarity << " "
                << " " << minExpAnchors << " " << alignmentPtrs[0]->qAlignedSeq.length << endl;

            if (params.nProc > 1) {
#ifdef __APPLE__
                sem_post(semaphores.hitCluster);
#else
                sem_post(&semaphores.hitCluster);
#endif
            }
        }

    }

    //
    // Assign the query name and strand for each alignment.
    //

    for (i = 0; i < alignmentPtrs.size(); i++) {
        T_AlignmentCandidate *aref = alignmentPtrs[i];
        if (aref->tStrand == 0) {
            aref->qName = read.GetName();
        }
        else {
            aref->qName = readRC.GetName();
        }
    }

    AssignRefContigLocations(alignmentPtrs, seqdb, genome);
}

template<typename T_Sequence>
void MapRead(T_Sequence &read, T_Sequence &readRC,
             vector<T_AlignmentCandidate*> &alignmentPtrs,
             MappingBuffers &mappingBuffers,
             MappingIPC *mapData,
             MappingSemaphores & semaphores)
{
    DNASuffixArray sarray;
    TupleCountTable<T_GenomeSequence, DNATuple> ct;
    SequenceIndexDatabase<FASTQSequence> seqdb;
    T_GenomeSequence    genome;
    BWT *bwtPtr = mapData->bwtPtr;
    mapData->ShallowCopySuffixArray(sarray);
    mapData->ShallowCopyReferenceSequence(genome);
    mapData->ShallowCopySequenceIndexDatabase(seqdb);
    mapData->ShallowCopyTupleCountTable(ct);
    SeqBoundaryFtr<FASTQSequence> seqBoundary(&seqdb);

    return
        MapRead(read, readRC,
                genome,           // possibly multi fasta file read into one sequence
                sarray, *bwtPtr,  // The suffix array, and the bwt-fm index structures
                seqBoundary,      // Boundaries of contigs in the
                // genome, alignments do not span
                // the ends of boundaries.
                ct,               // Count table to use word frequencies in the genome to weight matches.
                seqdb,            // Information about the names of
                // chromosomes in the genome, and
                // where their sequences are in the genome.
                mapData->params,// A huge list of parameters for
                // mapping, only compile/command
                // line values set.
                mapData->metrics, // Keep track of time/ hit counts,
                // etc.. Not fully developed, but
                // should be.
                alignmentPtrs,    // Where the results are stored.
                mappingBuffers,   // A class of buffers for structurs
                // like dyanmic programming
                // matrices, match lists, etc., that are not
                // reallocated between calls to
                // MapRead.  They are cleared though.
                mapData,          // Some values that are shared
                // across threads.
                semaphores);
}

template<typename T_TargetSequence, typename T_QuerySequence, typename TDBSequence>
void AlignIntervals(T_TargetSequence &genome, T_QuerySequence &read, T_QuerySequence &rcRead,
        WeightedIntervalSet &weightedIntervals,
        int mutationCostMatrix[][5],
        int ins, int del, int sdpTupleSize,
        int useSeqDB, SequenceIndexDatabase<TDBSequence> &seqDB,
        vector<T_AlignmentCandidate*> &alignments,
        MappingParameters &params,
        MappingBuffers &mappingBuffers,
        int procId) {

    vector<T_QuerySequence*> forrev;
    forrev.resize(2);
    forrev[Forward] = &read;
    forrev[Reverse] = &rcRead;

    //
    // Use an edit distance scoring function instead of IDS.  Although
    // the IDS should be more accurate, it is more slow, and it is more
    // important at this stage to have faster alignments than accurate,
    // since all alignments are rerun using GuidedAlignment later on.
    //
    DistanceMatrixScoreFunction<DNASequence, FASTQSequence> distScoreFn(SMRTDistanceMatrix, params.insertion, params.deletion);
    DistanceMatrixScoreFunction<DNASequence, FASTQSequence> distScoreFn2(SMRTDistanceMatrix, ins, ins);

    //
    // Assume there is at least one interval.
    //
    if (weightedIntervals.size() == 0)
        return;

    WeightedIntervalSet::iterator intvIt = weightedIntervals.begin();
    int alignmentIndex = 0;


    do {

        T_AlignmentCandidate *alignment = alignments[alignmentIndex];
        alignment->clusterWeight= (*intvIt).size; // totalAnchorSize == size
        alignment->clusterScore = (*intvIt).pValue;

        //
        // Advance references.  Intervals are stored in reverse order, so
        // go backwards in the list, and alignments are in forward order.
        // That should probably be changed.
        //
        ++alignmentIndex;

        //
        // Try aligning the read to the genome.
        //
        DNALength matchIntervalStart, matchIntervalEnd;
        matchIntervalStart = (*intvIt).start;
        matchIntervalEnd   = (*intvIt).end;

        bool readOverlapsContigStart    = false;
        bool readOverlapsContigEnd      = false;
        int  startOverlappedContigIndex = 0;
        int  endOverlappedContigIndex   = 0;
        if (params.verbosity > 0) {
            cout << "aligning interval : " << read.length << " " << (*intvIt).start << " "
                << (*intvIt).end  << " " << (*intvIt).qStart << " " << (*intvIt).qEnd
                << " " << matchIntervalStart << " to " << matchIntervalEnd << " "
                << params.approximateMaxInsertionRate << " "  << endl;
        }
        assert(matchIntervalEnd >= matchIntervalStart);

        //
        // If using a sequence database, check to make sure that the
        // boundaries of the sequence windows do not overlap with
        // the boundaries of the reads.  If the beginning is before
        // the boundary, move the beginning up to the start of the read.
        // If the end is past the end boundary of the read, similarly move
        // the window boundary to the end of the read boundary.

        DNALength tAlignedContigStart = 0;
        int seqDBIndex = 0;


        //
        // Stretch the alignment interval so that it is close to where
        // the read actually starts.
        //
        DNALength subreadStart = read.SubreadStart();
        DNALength subreadEnd   = read.SubreadEnd();
        if ((*intvIt).GetStrandIndex() == Reverse) {
            subreadEnd   = read.MakeRCCoordinate(read.SubreadStart()) + 1;
            subreadStart = read.MakeRCCoordinate(read.SubreadEnd()-1);
        }

        DNALength lengthBeforeFirstMatch = ((*intvIt).qStart - subreadStart) * params.approximateMaxInsertionRate ;
        DNALength lengthAfterLastMatch   = (subreadEnd - (*intvIt).qEnd) * params.approximateMaxInsertionRate;
        if (matchIntervalStart < lengthBeforeFirstMatch  or params.doGlobalAlignment) {
            matchIntervalStart = 0;
        }
        else {
            matchIntervalStart -= lengthBeforeFirstMatch;
        }

        if (genome.length < matchIntervalEnd + lengthAfterLastMatch or params.doGlobalAlignment) {
            matchIntervalEnd = genome.length;
        }
        else {
            matchIntervalEnd += lengthAfterLastMatch;
        }

        DNALength intervalContigStartPos, intervalContigEndPos;
        if (useSeqDB) {
            //
            // The sequence db index is the one where the actual match is
            // contained. The matchIntervalStart might be before the sequence
            // index boundary due to the extrapolation of alignment start by
            // insertion rate.  If this is the case, bump up the
            // matchIntervalStart to be at the beginning of the boundary.
            // Modify bounds similarly for the matchIntervalEnd and the end
            // of a boundary.
            //
            seqDBIndex = seqDB.SearchForIndex((*intvIt).start);
            intervalContigStartPos = seqDB.seqStartPos[seqDBIndex];
            if (intervalContigStartPos > matchIntervalStart) {
                matchIntervalStart = intervalContigStartPos;
            }
            intervalContigEndPos = seqDB.seqStartPos[seqDBIndex+1] - 1;
            if (intervalContigEndPos < matchIntervalEnd) {
                matchIntervalEnd = intervalContigEndPos;
            }
            alignment->tName    = seqDB.GetSpaceDelimitedName(seqDBIndex);
            alignment->tLength  = intervalContigEndPos - intervalContigStartPos;
            //
            // When there are multiple sequences in the database, store the
            // index of this sequence.  This lets one compare the contigs
            // that reads are mapped to, for instance.
            //
            alignment->tIndex   = seqDBIndex;
        }
        else {
            alignment->tLength     = genome.length;
            alignment->tName       = genome.GetName();
            intervalContigStartPos = 0;
            intervalContigEndPos   = genome.length;
            //
            // When there are multiple sequences in the database, store the
            // index of this sequence.  This lets one compare the contigs
            // that reads are mapped to, for instance.
            //
        }
        alignment->qName = read.title;
        //
        // Look to see if a read overhangs the beginning of a contig.
        //
        if (params.verbosity > 2) {
            cout << "Check for prefix/suffix overlap on interval: " << (*intvIt).qStart << " ?> " << (*intvIt).start - intervalContigStartPos <<endl;
        }
        if ( (*intvIt).qStart > (*intvIt).start - intervalContigStartPos) {
            readOverlapsContigStart = true;
            startOverlappedContigIndex = seqDBIndex;
        }

        //
        // Look to see if the read overhangs the end of a contig.
        //
        if (params.verbosity > 2) {
            cout << "Check for suffix/prefix overlap on interval, read overhang: " << read.length - (*intvIt).qEnd << " ?> " << matchIntervalEnd - (*intvIt).end  <<endl;
        }
        if (read.length - (*intvIt).qEnd > matchIntervalEnd - (*intvIt).end) {
            if (params.verbosity > 2) {
                cout << "read overlaps genome end." << endl;
            }
            readOverlapsContigEnd = true;
            endOverlappedContigIndex = seqDBIndex;
        }
        int alignScore;
        alignScore = 0;

        alignment->tAlignedSeqPos     = matchIntervalStart;
        alignment->tAlignedSeqLength  = matchIntervalEnd - matchIntervalStart;
        if ((*intvIt).GetStrandIndex() == Forward) {
            alignment->tAlignedSeq.Copy(genome, alignment->tAlignedSeqPos, alignment->tAlignedSeqLength);
            alignment->tStrand = Forward;
        }
        else {
            DNALength rcAlignedSeqPos = genome.MakeRCCoordinate(alignment->tAlignedSeqPos + alignment->tAlignedSeqLength - 1);
            genome.CopyAsRC(alignment->tAlignedSeq, rcAlignedSeqPos, alignment->tAlignedSeqLength);
            // Map forward coordinates into reverse complement.

            intervalContigStartPos    = genome.MakeRCCoordinate(intervalContigStartPos) + 1;
            intervalContigEndPos      = genome.MakeRCCoordinate(intervalContigEndPos - 1);
            swap(intervalContigStartPos, intervalContigEndPos);
            alignment->tAlignedSeqPos = rcAlignedSeqPos;
            alignment->tStrand        = Reverse;
        }

        // Configure the part of the query that is aligned.  The entire
        // query should always be aligned.
        alignment->qAlignedSeqPos    = 0;
        alignment->qAlignedSeq.ReferenceSubstring(read);
        alignment->qAlignedSeqLength = alignment->qAlignedSeq.length;
        alignment->qLength           = read.length;
        alignment->qStrand           = 0;

        if (params.verbosity > 1) {
            cout << "aligning read " << endl;
            static_cast<DNASequence*>(&(alignment->qAlignedSeq))->PrintSeq(cout);
            cout << endl << "aligning reference" << endl;
            static_cast<DNASequence*>(&(alignment->tAlignedSeq))->PrintSeq(cout);
            cout << endl;
        }

        //
        // The type of alignment that is performed depends on the mode
        // blasr is running in.  If it is running in normal mode, local
        // aligment is performed and guided by SDP alignment.  When
        // running in overlap mode, the alignments are forced to the ends
        // of reads.
        //

        int intervalSize = 0;
        int m;
        //
        // Check to see if the matches to the genome are sufficiently
        // dense to allow them to be used instead of having to redo
        // sdp alignment.
        //

        // First count how much of the read matches the genome exactly.
        for (m = 0; m < intvIt->matches.size(); m++) { intervalSize += intvIt->matches[m].l;}

        int subreadLength = forrev[(*intvIt).GetStrandIndex()]->SubreadEnd() - forrev[(*intvIt).GetStrandIndex()]->SubreadStart();
        if ((1.0*intervalSize) / subreadLength < params.sdpBypassThreshold and !params.emulateNucmer) {
            //
            // Not enough of the read maps to the genome, need to use
            // sdp alignment to define the regions of the read that map.
            //
            if (params.refineBetweenAnchorsOnly) {

                //
                // Run SDP alignment only between the genomic anchors,
                // including the genomic anchors as part of the alignment.
                //
                int m;

                vector<ChainedMatchPos> *matches;
                vector<ChainedMatchPos> rcMatches;
                Alignment anchorsOnly;
                DNASequence tAlignedSeq;
                FASTQSequence qAlignedSeq;
                //
                // The strand bookkeeping is a bit confusing, so hopefully
                // this will set things straight.
                //
                // If the alignment is forward strand, the coordinates of the
                // blocks are relative to the forward read, starting at 0, not
                // the subread start.
                // If the alignment is reverse strand, the coordinates of the
                // blocks are relative to the reverse strand, starting at the
                // position of the subread on the reverse strand.
                //
                // The coordinates of the blocks in the genome are always
                // relative to the forward strand on the genome, starting at
                // 0.
                //

                //
                // The first step to refining between anchors only is to make
                // the anchors relative to the tAlignedSeq.

                matches = (vector<ChainedMatchPos>*) &(*intvIt).matches;
                tAlignedSeq = alignment->tAlignedSeq;
                qAlignedSeq = alignment->qAlignedSeq;

                if (alignment->tStrand == 0) {
                    for (m = 0; m < matches->size(); m++) {
                        (*matches)[m].t -= alignment->tAlignedSeqPos;
                        (*matches)[m].q -= alignment->qAlignedSeqPos;
                    }
                }
                else  {
                    //
                    // Flip the entire alignment if it is on the reverse strand.
                    DNALength rcAlignedSeqPos = genome.MakeRCCoordinate(alignment->tAlignedSeqPos + alignment->tAlignedSeqLength - 1);
                    for (m = 0; m < matches->size(); m++) {
                        (*matches)[m].t -= rcAlignedSeqPos;
                        (*matches)[m].q -= alignment->qAlignedSeqPos;
                    }

                    alignment->tAlignedSeq.CopyAsRC(tAlignedSeq);
                    rcMatches.resize((*intvIt).matches.size());
                    //
                    // Make the reverse complement of the match list.
                    //

                    // 1. Reverse complement the coordinates.
                    for (m = 0; m < (*intvIt).matches.size(); m++) {
                        int revCompIndex = rcMatches.size() - m - 1;
                        rcMatches[revCompIndex].q = read.MakeRCCoordinate((*intvIt).matches[m].q + (*intvIt).matches[m].l - 1);
                        rcMatches[revCompIndex].t = tAlignedSeq.MakeRCCoordinate((*intvIt).matches[m].t + (*intvIt).matches[m].l - 1);
                        rcMatches[revCompIndex].l = (*intvIt).matches[m].l;
                    }
                    matches = &rcMatches;
                }

                /*
                   Uncomment to get a dot plot
                   ofstream matchFile;
                   matchFile.open("matches.txt");
                   matchFile << "q t l " << endl;
                   for (m = 0; matches->size() > 0 and m < matches->size() - 1; m++) {
                   matchFile << (*matches)[m].q << " " << (*matches)[m].t << " " << (*matches)[m].l << endl;
                   }
                   */
                DNASequence tSubSeq;
                FASTQSequence qSubSeq;
                for (m = 0; matches->size() > 0 and m < matches->size() - 1; m++) {
                    Block block;
                    block.qPos = (*matches)[m].q;
                    block.tPos = (*matches)[m].t;
                    block.length = (*matches)[m].l;

                    //
                    // Find the lengths of the gaps between anchors.
                    //
                    int tGap, qGap;
                    tGap = (*matches)[m+1].t - ((*matches)[m].t + (*matches)[m].l);
                    qGap = (*matches)[m+1].q - ((*matches)[m].q + (*matches)[m].l);
                    float gapRatio = (1.0*tGap)/qGap;

                    if (tGap > 0 and qGap > 0) {
                        DNALength tPos, qPos;
                        tPos = block.tPos + block.length;
                        qPos = block.qPos + block.length;
                        tSubSeq.ReferenceSubstring(tAlignedSeq, tPos, tGap);
                        qSubSeq.ReferenceSubstring(alignment->qAlignedSeq, qPos, qGap);
                        Alignment alignmentInGap;
                        int alignScore;

                        /*
                           The following code is experimental code for trying to do
                           something like affine gap alignment in long gaps.  It
                           would eventually be used in cDNA alignment to align
                           between exons, but for now is being tested here by using
                           it to align when there is a big gap between anchors.
                           */
                        if (params.separateGaps == true and
                                qSubSeq.length > 0 and tSubSeq.length > 0 and
                                ( (1.0*qSubSeq.length)/tSubSeq.length  < 0.25 )) {
                            alignScore = OneGapAlign(qSubSeq, tSubSeq, distScoreFn, mappingBuffers, alignmentInGap);
                        }
                        else {
                            /*
                               This is the 'normal/default' way to align between
                               gaps.  It is more well tested than OneGapAlign.
                               */
                            alignScore = SDPAlign(qSubSeq, tSubSeq, distScoreFn, params.sdpTupleSize,
                                    params.sdpIns, params.sdpDel, params.indelRate*2,
                                    alignmentInGap, mappingBuffers, Global,
                                    params.detailedSDPAlignment,
                                    params.extendFrontAlignment,
                                    params.recurseOver,
                                    params.fastSDP);
                        }

                        //
                        // Now, splice the fragment alignment into the current
                        // alignment.
                        //
                        if (alignmentInGap.blocks.size() > 0) {
                            int b;
                            //
                            // Configure this block to be relative to the beginning
                            // of the aligned substring.
                            //
                            for (b = 0; b < alignmentInGap.size(); b++) {
                                alignmentInGap.blocks[b].tPos += tPos + alignmentInGap.tPos;
                                alignmentInGap.blocks[b].qPos += qPos + alignmentInGap.qPos;
                                assert(alignmentInGap.blocks[b].tPos < alignment->tAlignedSeq.length);
                                assert(alignmentInGap.blocks[b].qPos < alignment->qAlignedSeq.length);
                            }
                        }
                        // Add the original block
                        alignment->blocks.push_back(block);
                        anchorsOnly.blocks.push_back(block);
                        // Add the blocks for the refined alignment
                        alignment->blocks.insert(alignment->blocks.end(),
                                alignmentInGap.blocks.begin(),
                                alignmentInGap.blocks.end());
                    }
                }

                // Add the last block
                m = (*matches).size() - 1;
                Block block;
                block.qPos = (*matches)[m].q;
                block.tPos = (*matches)[m].t;

                assert(block.tPos <= alignment->tAlignedSeq.length);
                assert(block.qPos <= alignment->qAlignedSeq.length);

                block.length = (*matches)[m].l;
                alignment->blocks.push_back(block);
                anchorsOnly.blocks.push_back(block);

                //
                // By convention, blocks start at 0, and the
                // alignment->tPos,qPos give the start of the alignment.
                // Modify the block positions so that they are offset by 0.
                alignment->tPos = alignment->blocks[0].tPos;
                alignment->qPos = alignment->blocks[0].qPos;
                int b;
                int blocksSize = alignment->blocks.size();
                for (b = 0; b < blocksSize ; b++) {
                    assert(alignment->tPos <= alignment->blocks[b].tPos);
                    assert(alignment->qPos <= alignment->blocks[b].qPos);
                    alignment->blocks[b].tPos -= alignment->tPos;
                    alignment->blocks[b].qPos -= alignment->qPos;
                }
                for (b = 0; b < anchorsOnly.blocks.size(); b++) {
                    anchorsOnly.blocks[b].tPos -= alignment->tPos;
                    anchorsOnly.blocks[b].qPos -= alignment->qPos;
                }
                anchorsOnly.tPos = alignment->tPos;
                anchorsOnly.qPos = alignment->qPos;
                ComputeAlignmentStats(*alignment, alignment->qAlignedSeq.seq, alignment->tAlignedSeq.seq,
                        distScoreFn);

                tAlignedSeq.Free();
                qAlignedSeq.Free();
                tSubSeq.Free();
                qSubSeq.Free();
            }
            else {
                alignScore = SDPAlign(alignment->qAlignedSeq, alignment->tAlignedSeq, distScoreFn,
                        sdpTupleSize, params.sdpIns, params.sdpDel, params.indelRate*3,
                        *alignment, mappingBuffers,
                        Local,
                        params.detailedSDPAlignment,
                        params.extendFrontAlignment,
                        params.recurseOver,
                        params.fastSDP);
                ComputeAlignmentStats(*alignment, alignment->qAlignedSeq.seq, alignment->tAlignedSeq.seq,
                        distScoreFn);
            }
        }
        else {
            //
            // The anchors used to anchor the sequence are sufficient to extend the alignment.
            //
            int m;
            for (m = 0; m < (*intvIt).matches.size(); m++ ){
                Block block;
                block.qPos = (*intvIt).matches[m].q - alignment->qAlignedSeqPos;
                block.tPos = (*intvIt).matches[m].t - alignment->tAlignedSeqPos;
                block.length = (*intvIt).matches[m].l;
                alignment->blocks.push_back(block);
            }
        }

        //
        //  The anchors/sdp alignments may leave portions of the read
        //  unaligned at the beginning and end.  If the parameters
        //  specify extending alignments, try and align extra bases at
        //  the beginning and end of alignments.
        if (params.extendAlignments) {

            //
            // Modify the alignment so that the start and end of the
            // alignment strings are at the alignment boundaries.
            //
            // Since the query sequence is pointing at a subsequence of the
            // read (and is always in the forward direction), just reference
            // a new portion of the read.
            alignment->qAlignedSeqPos = alignment->qAlignedSeqPos + alignment->qPos;
            alignment->qAlignedSeqLength = alignment->QEnd();
            alignment->qAlignedSeq.ReferenceSubstring(read, alignment->qAlignedSeqPos, alignment->qAlignedSeqLength );
            alignment->qPos = 0;

            //
            // Since the target sequence may be on the forward or reverse
            // strand, a copy of the subsequence is made, and the original
            // sequence free'd.
            //
            DNASequence tSubseq;
            alignment->tAlignedSeqPos = alignment->tAlignedSeqPos + alignment->tPos;
            alignment->tAlignedSeqLength = alignment->TEnd();
            tSubseq.Copy(alignment->tAlignedSeq, alignment->tPos, alignment->tAlignedSeqLength);
            alignment->tPos = 0;

            alignment->tAlignedSeq.Free();
            alignment->tAlignedSeq.TakeOwnership(tSubseq);

            DNALength maximumExtendLength = 500;

            if (alignment->blocks.size() > 0 ) {
                int lastAlignedBlock = alignment->blocks.size() - 1;
                DNALength lastAlignedQPos  = alignment->blocks[lastAlignedBlock].QEnd() + alignment->qPos + alignment->qAlignedSeqPos;
                DNALength lastAlignedTPos  = alignment->blocks[lastAlignedBlock].TEnd() + alignment->tPos + alignment->tAlignedSeqPos;
                T_AlignmentCandidate extendedAlignmentForward, extendedAlignmentReverse;
                int forwardScore, reverseScore;

                SMRTSequence  readSuffix;
                DNALength     readSuffixLength;
                DNASequence   genomeSuffix;
                DNALength     genomeSuffixLength;

                SMRTSequence   readPrefix;
                DNALength     readPrefixLength;
                DNASequence   genomePrefix;
                DNALength     genomePrefixLength;

                //
                // Align the entire end of the read if it is short enough.
                //
                readSuffixLength = min(read.length - lastAlignedQPos, maximumExtendLength);
                if (readSuffixLength > 0) {
                    readSuffix.ReferenceSubstring(read, lastAlignedQPos, readSuffixLength);
                }
                else {
                    readSuffix.length = 0;
                }

                //
                // Align The entire end of the genome up to the maximum extend length;
                //
                genomeSuffixLength = min(intervalContigEndPos - lastAlignedTPos, maximumExtendLength);
                if (genomeSuffixLength > 0) {
                    if (alignment->tStrand == Forward) {
                        genomeSuffix.Copy(genome, lastAlignedTPos, genomeSuffixLength);
                    }
                    else {
                        static_cast<DNASequence*>(&genome)->CopyAsRC(genomeSuffix, lastAlignedTPos, genomeSuffixLength);
                    }
                }
                else {
                    genomeSuffix.length = 0;
                }
                forwardScore = 0;
                if (readSuffix.length > 0 and genomeSuffix.length > 0) {
                    forwardScore = ExtendAlignmentForward(readSuffix, 0,
                            genomeSuffix, 0,
                            params.extendBandSize,
                            // Reuse buffers to speed up alignment
                            mappingBuffers.scoreMat,
                            mappingBuffers.pathMat,
                            // Do the alignment in the forward direction.
                            extendedAlignmentForward,
                            distScoreFn,
                            1, // don't bother attempting
                            // to extend the alignment
                            // if one of the sequences
                            // is less than 1 base long
                            params.maxExtendDropoff);
                }

                if ( forwardScore < 0 ) {
                    //
                    // The extended alignment considers the whole genome, but
                    // should be modified to be starting at the end of where
                    // the original alignment left off.
                    //
                    if (params.verbosity > 0) {
                        cout << "forward extended an alignment of score " << alignment->score << " with score " << forwardScore << " by " << extendedAlignmentForward.blocks.size() << " blocks and length " << extendedAlignmentForward.blocks[extendedAlignmentForward.blocks.size()-1].qPos << endl;
                    }
                    extendedAlignmentForward.tAlignedSeqPos = lastAlignedTPos;

                    extendedAlignmentForward.qAlignedSeqPos = lastAlignedQPos;

                    genomeSuffix.length = extendedAlignmentForward.tPos + extendedAlignmentForward.TEnd();
                    alignment->tAlignedSeq.Append(genomeSuffix);
                    alignment->qAlignedSeq.length += extendedAlignmentForward.qPos + extendedAlignmentForward.QEnd();
                    assert(alignment->qAlignedSeq.length <= read.length);
                    alignment->AppendAlignment(extendedAlignmentForward);
                }

                DNALength firstAlignedQPos = alignment->qPos + alignment->qAlignedSeqPos;
                DNALength firstAlignedTPos = alignment->tPos + alignment->tAlignedSeqPos;

                readPrefixLength = min(firstAlignedQPos, maximumExtendLength);
                if (readPrefixLength > 0) {
                    readPrefix.ReferenceSubstring(read, firstAlignedQPos-readPrefixLength, readPrefixLength);
                }
                else {
                    readPrefix.length = 0;
                }

                genomePrefixLength = min(firstAlignedTPos - intervalContigStartPos, maximumExtendLength);
                if (genomePrefixLength > 0) {
                    if (alignment->tStrand == 0) {
                        genomePrefix.Copy(genome, firstAlignedTPos - genomePrefixLength, genomePrefixLength);
                    }
                    else {
                        static_cast<DNASequence*>(&genome)->MakeRC(genomePrefix, firstAlignedTPos - genomePrefixLength, genomePrefixLength);
                    }
                }
                reverseScore = 0;
                if (readPrefix.length > 0 and genomePrefix.length > 0) {
                    reverseScore = ExtendAlignmentReverse(readPrefix, readPrefix.length-1,
                            genomePrefix, genomePrefixLength - 1,
                            params.extendBandSize, //k
                            mappingBuffers.scoreMat,
                            mappingBuffers.pathMat,
                            extendedAlignmentReverse,
                            distScoreFn,
                            1, // don't bother attempting
                            // to extend the alignment
                            // if one of the sequences
                            // is less than 1 base long
                            params.maxExtendDropoff);
                }

                if (reverseScore < 0 ) {
                    //
                    // Make alignment->tPos relative to the beginning of the
                    // extended alignment so that when it is appended, the
                    // coordinates match correctly.
                    if (params.verbosity > 0) {
                        cout << "reverse extended an alignment of score " << alignment->score << " with score " << reverseScore << " by " << extendedAlignmentReverse.blocks.size() << " blocks and length " << extendedAlignmentReverse.blocks[extendedAlignmentReverse.blocks.size()-1].qPos << endl;
                    }
                    extendedAlignmentReverse.tAlignedSeqPos = firstAlignedTPos - genomePrefixLength;
                    extendedAlignmentReverse.qAlignedSeqPos = firstAlignedQPos - readPrefixLength;
                    extendedAlignmentReverse.AppendAlignment(*alignment);

                    genomePrefix.Append(alignment->tAlignedSeq, genomePrefix.length - alignment->tPos);
                    alignment->tAlignedSeq.Free();
                    alignment->tAlignedSeq.TakeOwnership(genomePrefix);

                    alignment->blocks = extendedAlignmentReverse.blocks;

                    alignment->tAlignedSeqPos = extendedAlignmentReverse.tAlignedSeqPos;
                    alignment->tPos = extendedAlignmentReverse.tPos;


                    alignment->qAlignedSeqPos     = extendedAlignmentReverse.qAlignedSeqPos;
                    alignment->qAlignedSeq.length = readPrefix.length + alignment->qAlignedSeq.length;
                    alignment->qPos               = extendedAlignmentReverse.qPos;
                    alignment->qAlignedSeq.seq    = readPrefix.seq;
                    //
                    // Make sure the two ways of accounting for aligned sequence
                    // length are in sync.  This needs to go.
                    //
                    if (alignment->blocks.size() > 0) {
                        int lastBlock = alignment->blocks.size() - 1;
                        alignment->qAlignedSeqLength = alignment->qAlignedSeq.length;
                        alignment->tAlignedSeqLength = alignment->tAlignedSeq.length;
                    }
                    else {
                        alignment->qAlignedSeqLength = alignment->qAlignedSeq.length = 0;
                        alignment->tAlignedSeqLength = alignment->tAlignedSeq.length = 0;
                    }
                } // end of if (reverseScore < 0 )
                readSuffix.Free();
                readPrefix.Free();
                genomePrefix.Free();
                genomeSuffix.Free();
            }
            tSubseq.Free();
        }

        if (params.verbosity > 0) {
            cout << "interval align score: " << alignScore << endl;
            StickPrintAlignment(*alignment,
                    (DNASequence&) alignment->qAlignedSeq,
                    (DNASequence&) alignment->tAlignedSeq,
                    cout,
                    0, alignment->tAlignedSeqPos);

        }
        ComputeAlignmentStats(*alignment,
                alignment->qAlignedSeq.seq,
                alignment->tAlignedSeq.seq,
                distScoreFn2);
        //SMRTDistanceMatrix, ins, del );


        intvIt++;
    } while (intvIt != weightedIntervals.end());
}


template<typename T_RefSequence, typename T_Sequence>
void PairwiseLocalAlign(T_Sequence &qSeq, T_RefSequence &tSeq,
        int k,
        MappingParameters &params, T_AlignmentCandidate &alignment,
        MappingBuffers &mappingBuffers,
        AlignmentType alignType) {
    //
    // Perform a pairwise alignment between qSeq and tSeq, but choose
    // the pairwise alignment method based on the parameters.  The
    // options for pairwise alignment are:
    //  - Affine KBanded alignment: usually used for sequences with no
    //                              quality information.
    //  - KBanded alignment: For sequences with quality information.
    //                       Gaps are scored with quality values.
    //
    QualityValueScoreFunction<DNASequence, FASTQSequence> scoreFn;
    scoreFn.del = params.indel;
    scoreFn.ins = params.indel;

    DistanceMatrixScoreFunction<DNASequence, FASTASequence> distScoreFn2(
            SMRTDistanceMatrix, params.indel, params.indel);

    IDSScoreFunction<DNASequence, FASTQSequence> idsScoreFn;
    idsScoreFn.ins = params.insertion;
    idsScoreFn.del = params.deletion;
    idsScoreFn.substitutionPrior = params.substitutionPrior;
    idsScoreFn.globalDeletionPrior = params.globalDeletionPrior;
    idsScoreFn.InitializeScoreMatrix(SMRTDistanceMatrix);

    int kbandScore;
    int qvAwareScore;
    if (params.ignoreQualities || qSeq.qual.Empty() || !ReadHasMeaningfulQualityValues(qSeq) ) {

        kbandScore = AffineKBandAlign(qSeq, tSeq, SMRTDistanceMatrix,
                params.indel+2, params.indel - 3, // homopolymer insertion open and extend
                params.indel+2, params.indel - 1, // any insertion open and extend
                params.indel, // deletion
                k*1.2,
                mappingBuffers.scoreMat, mappingBuffers.pathMat,
                mappingBuffers.hpInsScoreMat, mappingBuffers.hpInsPathMat,
                mappingBuffers.insScoreMat, mappingBuffers.insPathMat,
                alignment, Global);

        alignment.score = kbandScore;
        if (params.verbosity >= 2) {
            cout << "align score: " << kbandScore << endl;
        }
    }
    else {


        if (qSeq.insertionQV.Empty() == false) {
            qvAwareScore = KBandAlign(qSeq, tSeq, SMRTDistanceMatrix,
                    params.indel+2, // ins
                    params.indel+2, // del
                    k,
                    mappingBuffers.scoreMat, mappingBuffers.pathMat,
                    alignment, idsScoreFn, alignType);
            if (params.verbosity >= 2) {
                cout << "ids score fn score: " << qvAwareScore << endl;
            }
        }
        else {
            qvAwareScore = KBandAlign(qSeq, tSeq, SMRTDistanceMatrix,
                    params.indel+2, // ins
                    params.indel+2, // del
                    k,
                    mappingBuffers.scoreMat, mappingBuffers.pathMat,
                    alignment, scoreFn, alignType);
            if (params.verbosity >= 2) {
                cout << "qv score fn score: " << qvAwareScore << endl;
            }
        }
        alignment.sumQVScore = qvAwareScore;
        alignment.score = qvAwareScore;
        alignment.probScore = 0;
    }
    // Compute stats and assign a default alignment score using an edit distance.
    ComputeAlignmentStats(alignment, qSeq.seq, tSeq.seq, distScoreFn2);

    if (params.scoreType == 1) {
        alignment.score = alignment.sumQVScore;
    }
}

// Extend target aligned sequence of the input alignement to both ends
// by flankSize bases. Update alignment->tAlignedSeqPos,
// alignment->tAlignedSeqLength and alignment->tAlignedSeq.
void FlankTAlignedSeq(T_AlignmentCandidate * alignment,
                      SequenceIndexDatabase<FASTQSequence> &seqdb,
                      DNASequence & genome,
                      int flankSize) {
    assert(alignment != NULL and alignment->tIsSubstring);

    UInt forwardTPos, newTAlignedSeqPos, newTAlignedSeqLen;
    // New aligned start position relative to this chromosome, with
    // the same direction as alignment->tStrand.
    newTAlignedSeqPos = UInt((alignment->tAlignedSeqPos > UInt(flankSize))?
            (alignment->tAlignedSeqPos - flankSize): 0);
    newTAlignedSeqLen = min(alignment->tAlignedSeqPos + alignment->tAlignedSeqLength +
            flankSize, alignment->tLength) - newTAlignedSeqPos;

    if (alignment->tStrand ==0) {
        forwardTPos = newTAlignedSeqPos;
    } else {
        forwardTPos = alignment->tLength - newTAlignedSeqPos - 1;
    }

    // Find where this chromosome is in the genome.
    int seqIndex = seqdb.GetIndexOfSeqName(alignment->tName);
    assert(seqIndex != -1);
    UInt newGenomePos = seqdb.ChromosomePositionToGenome(seqIndex, forwardTPos);

    if (alignment->tIsSubstring == false) {
        alignment->tAlignedSeq.Free();
    }
    alignment->tAlignedSeqPos = newTAlignedSeqPos;
    alignment->tAlignedSeqLength = newTAlignedSeqLen;
    if (alignment->tStrand == 0) {
        alignment->tAlignedSeq.ReferenceSubstring(genome, newGenomePos, newTAlignedSeqLen);
    } else {
        // Copy and then reverse complement.
        genome.MakeRC(alignment->tAlignedSeq,
                newGenomePos + 1 - alignment->tAlignedSeqLength,
                alignment->tAlignedSeqLength);
        alignment->tIsSubstring = false;
    }
}

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
        ostream & threadOut) {
    assert(passDirection == 0 or passDirection == 1);
    //
    // Determine where in the genome the subread has mapped.
    //
    DNASequence alignedForwardRefSequence, alignedReverseRefSequence;

    if (alignment->tStrand == 0) {
        // This needs to be changed -- copy copies RHS into LHS,
        // CopyAsRC copies LHS into RHS
        alignedForwardRefSequence.Copy(alignment->tAlignedSeq);
        alignment->tAlignedSeq.CopyAsRC(alignedReverseRefSequence);
    }
    else {
        alignment->tAlignedSeq.CopyAsRC(alignedForwardRefSequence);
        alignedReverseRefSequence.Copy(alignment->tAlignedSeq);
    }

    IDSScoreFunction<DNASequence, FASTQSequence> idsScoreFn;
    idsScoreFn.ins  = params.insertion;
    idsScoreFn.del  = params.deletion;
    idsScoreFn.InitializeScoreMatrix(SMRTDistanceMatrix);
    idsScoreFn.globalDeletionPrior = params.globalDeletionPrior;
    idsScoreFn.substitutionPrior   = params.substitutionPrior;

    DistanceMatrixScoreFunction<DNASequence, FASTQSequence> distScoreFn2(
            SMRTDistanceMatrix, params.indel, params.indel);
    //
    // Determine the strand to align the subread to.
    //
    T_AlignmentCandidate exploded;
    bool sameAlignmentPassDirection = (alignment->tStrand == passDirection);
    bool computeProbIsFalse = false;
    DNASequence & alignedRefSequence = (sameAlignmentPassDirection?
            alignedForwardRefSequence:alignedReverseRefSequence);
    //
    // In the original code, parameters: bandSize=10, alignType=Global,
    // sdpTupleSize=4 (instead of 12, Local and 6) were used when
    // alignment & pass have different directions.
    //
    int explodedScore = GuidedAlign(subread, alignedRefSequence,
            idsScoreFn, 12, params.sdpIns, params.sdpDel,
            params.indelRate, mappingBuffers, exploded,
            Local, computeProbIsFalse, 6);

    if (params.verbosity >= 3) {
        threadOut << "zmw " << unrolledRead.zmwData.holeNumber
            << ", subreadIndex " << subreadIndex
            << ", passDirection " << passDirection
            << ", subreadInterval [" << subreadInterval.start
            << ", " << subreadInterval.end << ")" << endl
            << "StickPrintAlignment subread-reference alignment which has"
            << " the " << (sameAlignmentPassDirection?"same":"different")
            << " direction as the ccs-reference (or the "
            << "longestSubread-reference) alignment. " << endl
            << "subread: " << endl;
        static_cast<DNASequence*>(&subread)->PrintSeq(threadOut);
        threadOut << endl;
        threadOut << "alignedRefSeq: " << endl;
        static_cast<DNASequence*>(&alignedRefSequence)->PrintSeq(threadOut);
        StickPrintAlignment(exploded, (DNASequence&) subread,
                (DNASequence&) alignedRefSequence,
                threadOut, exploded.qAlignedSeqPos,
                exploded.tAlignedSeqPos);
    }

    if (exploded.blocks.size() > 0) {
        DistanceMatrixScoreFunction<DNASequence, FASTQSequence> distScoreFn(
                SMRTDistanceMatrix, params.indel, params.indel);
        ComputeAlignmentStats(exploded, subread.seq,
                alignedRefSequence.seq,
                distScoreFn2);
        //SMRTDistanceMatrix, params.indel, params.indel);
        if (exploded.score <= params.maxScore) {
            //
            // The coordinates of the alignment should be
            // relative to the reference sequence (the specified chromosome,
            // not the whole genome).
            //
            exploded.qStrand = 0;
            exploded.tStrand = sameAlignmentPassDirection?0:1;
            exploded.qLength = unrolledRead.length;
            exploded.tLength = alignment->tLength;
            exploded.tAlignedSeq.Copy(alignedRefSequence);
            exploded.tAlignedSeqPos = (passDirection == 0)?
                (alignment->tAlignedSeqPos):
                (exploded.tLength - alignment->tAlignedSeqPos
                 - alignment->tAlignedSeqLength);
            exploded.tAlignedSeqLength = alignment->tAlignedSeqLength;

            exploded.qAlignedSeq.ReferenceSubstring(subread);
            exploded.qAlignedSeqPos = subreadInterval.start;
            exploded.qAlignedSeqLength = subreadInterval.end - subreadInterval.start;
            exploded.mapQV = alignment->mapQV;
            exploded.tName = alignment->tName;
            exploded.tIndex = alignment->tIndex;

            stringstream namestrm;
            namestrm << "/" << subreadInterval.start
                << "_" << subreadInterval.end;
            exploded.qName = string(unrolledRead.title) + namestrm.str();

            //
            // Don't call AssignRefContigLocation as the coordinates
            // of the alignment is already relative to the chromosome coordiantes.
            //
            // Save this alignment for printing later.
            //
            T_AlignmentCandidate *alignmentPtr = new T_AlignmentCandidate;
            *alignmentPtr = exploded;
            allReadAlignments.AddAlignmentForSeq(subreadIndex, alignmentPtr);
        } // End of exploded score <= maxScore.
        if (params.verbosity >= 3) {
            threadOut << "exploded score: " << exploded.score << endl
                << "exploded alignment: "<< endl;
            exploded.Print(threadOut);
            threadOut << endl;
        }
    } // End of exploded.blocks.size() > 0.
}
