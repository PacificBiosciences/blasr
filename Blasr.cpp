// Author: Mark Chaisson

#include "iblasr/BlasrAlign.hpp"
#include "iblasr/BlasrMiscs.hpp"
#include "iblasr/BlasrUtils.hpp"
#include "iblasr/RegisterBlasrOptions.h"

#if CMAKE_BUILD
#include "BlasrVersion.h"
#endif

//#define USE_GOOGLE_PROFILER
#ifdef USE_GOOGLE_PROFILER
#include "gperftools/profiler.h"
#endif

// Declare global structures that are shared between threads.
MappingSemaphores semaphores;
std::ostream *outFilePtr = NULL;
#ifdef USE_PBBAM
PacBio::BAM::IRecordWriter *bamWriterPtr = NULL;  // use IRecordWriter for both SAM ands BAM
#endif

HDFRegionTableReader *regionTableReader = NULL;
ReaderAgglomerate *reader = NULL;

// Add comment to version history for each version change !
//
// Version history
//
// 5.0 - a new major version number
// 5.1 - transiotion to POSIX notation - double sashes before multi-character flags
// 5.2 - --sam no longer supported
// 5.3 - --sam supported via pbbam/IRecordWriter
//
const std::string GetMajorVersion() { return PacBio::BlasrMajorVersion(); }

// version format is 3 numbers separated by dots, with an optional
// SHA1 separated by a dash: Version.Subversion.Patchversion[-SHA1]
const std::string GetVersion(void) { return PacBio::BlasrVersion() + PacBio::BlasrGitSha1(); }

/// Checks whether a smrtRead meets the following criteria
/// (1) is within the search holeNumber range specified by params.holeNumberRanges.
/// (2) its length greater than params.maxReadlength
/// (3) its read score (rq) is greater than params.minRawSubreadScore
/// (4) its qual is greater than params.minAvgQual.
/// Change stop to false if
/// HoleNumber of the smrtRead is greater than the search holeNumber range.
bool IsGoodRead(const SMRTSequence &smrtRead, MappingParameters &params, bool &stop)
{
    if (params.holeNumberRangesStr.size() > 0 and
        not params.holeNumberRanges.contains(smrtRead.HoleNumber())) {
        // Stop processing once the specified zmw hole number is reached.
        // Eventually this will change to just seek to hole number, and
        // just align one read anyway.
        if (smrtRead.HoleNumber() > params.holeNumberRanges.max()) {
            stop = true;
            return false;
        }
        return false;
    }
    //
    // Discard reads that are too small, or not labeled as having any
    // useable/good sequence.
    //
    if (smrtRead.highQualityRegionScore < params.minRawSubreadScore or
        (params.maxReadLength != 0 and smrtRead.length > UInt(params.maxReadLength)) or
        (int(smrtRead.length) < params.minReadLength)) {
        return false;
    }

    if (smrtRead.qual.Empty() != false and smrtRead.GetAverageQuality() < params.minAvgQual) {
        return false;
    }
    return true;
}

// Make primary intervals (which are intervals of subreads to align
// in the first round) from none BAM file using region table.
void MakePrimaryIntervals(RegionTable *regionTablePtr, SMRTSequence &smrtRead,
                          std::vector<ReadInterval> &subreadIntervals,
                          std::vector<int> &subreadDirections, int &bestSubreadIndex,
                          MappingParameters &params)
{
    std::vector<ReadInterval> adapterIntervals;
    //
    // Determine endpoints of this subread in the main read.
    //
    if (params.useRegionTable == false) {
        //
        // When there is no region table, the subread is the entire
        // read.
        //
        ReadInterval wholeRead(0, smrtRead.length);
        // The set of subread intervals is just the entire read.
        subreadIntervals.push_back(wholeRead);
    } else {
        //
        // Grab the subread & adapter intervals from the entire region table to
        // iterate over.
        //
        assert(regionTablePtr->HasHoleNumber(smrtRead.HoleNumber()));
        subreadIntervals = (*regionTablePtr)[smrtRead.HoleNumber()].SubreadIntervals(
            smrtRead.length, params.byAdapter);
        adapterIntervals = (*regionTablePtr)[smrtRead.HoleNumber()].AdapterIntervals();
    }

    // The assumption is that neighboring subreads must have the opposite
    // directions. So create directions for subread intervals with
    // interleaved 0s and 1s.
    CreateDirections(subreadDirections, subreadIntervals.size());

    //
    // Trim the boundaries of subread intervals so that only high quality
    // regions are included in the intervals, not N's. Remove intervals
    // and their corresponding dirctions, if they are shorter than the
    // user specified minimum read length or do not intersect with hq
    // region at all. Finally, return index of the (left-most) longest
    // subread in the updated vector.
    //
    int longestSubreadIndex =
        GetHighQualitySubreadsIntervals(subreadIntervals,   // a vector of subread intervals.
                                        subreadDirections,  // a vector of subread directions.
                                        smrtRead.lowQualityPrefix,  // hq region start pos.
                                        smrtRead.length - smrtRead.lowQualitySuffix,  // hq end pos.
                                        params.minSubreadLength);  // minimum read length.

    bestSubreadIndex = longestSubreadIndex;
    if (params.concordantTemplate == "longestsubread") {
        // Use the (left-most) longest full-pass subread as
        // template for concordant mapping
        int longestFullSubreadIndex =
            GetLongestFullSubreadIndex(subreadIntervals, adapterIntervals);
        if (longestFullSubreadIndex >= 0) {
            bestSubreadIndex = longestFullSubreadIndex;
        }
    } else if (params.concordantTemplate == "typicalsubread") {
        // Use the 'typical' full-pass subread as template for
        // concordant mapping.
        int typicalFullSubreadIndex =
            GetTypicalFullSubreadIndex(subreadIntervals, adapterIntervals);
        if (typicalFullSubreadIndex >= 0) {
            bestSubreadIndex = typicalFullSubreadIndex;
        }
    } else if (params.concordantTemplate == "mediansubread") {
        // Use the 'median-length' full-pass subread as template for
        // concordant mapping.
        int medianFullSubreadIndex =
            GetMedianLengthFullSubreadIndex(subreadIntervals, adapterIntervals);
        if (medianFullSubreadIndex >= 0) {
            bestSubreadIndex = medianFullSubreadIndex;
        }
    } else {
        assert(false);
    }
}

// Make primary intervals (which are intervals of subreads to align
// in the first round) for BAM file, -concordant,
void MakePrimaryIntervals(std::vector<SMRTSequence> &subreads,
                          std::vector<ReadInterval> &subreadIntervals,
                          std::vector<int> &subreadDirections, int &bestSubreadIndex)
{
    MakeSubreadIntervals(subreads, subreadIntervals);
    CreateDirections(subreadDirections, subreadIntervals.size());
    bestSubreadIndex = GetIndexOfConcordantTemplate(subreadIntervals);
}

/// Scan the next read from input.  This may either be a CCS read, unrolled (Polymerase) read,
/// or regular read (though this may be aligned in whole, or by
/// subread).
/// \params[in] reader: FASTA/FASTQ/BAX.H5/CCS.H5/BAM file reader
/// \params[in] regionTablePtr: RGN.H5 region table pointer.
/// \params[in] params: mapping parameters.
/// \params[out] smrtRead: to save smrt sequence.
/// \params[out] ccsRead: to save ccs sequence.
/// \params[out] readIsCCS: read is CCSSequence.
/// \params[out] readGroupId: associated read group id
/// \params[out] associatedRandInt: random int associated with this zmw,
///              required to for generating deterministic random
///              alignments regardless of nproc.
/// \params[out] stop: whether or not stop mapping remaining reads.
/// \returns whether or not to skip mapping reads of this zmw.
bool FetchReads(ReaderAgglomerate *reader, RegionTable *regionTablePtr, SMRTSequence &smrtRead,
                CCSSequence &ccsRead, std::vector<SMRTSequence> &subreads,
                MappingParameters &params, bool &readIsCCS, std::string &readGroupId,
                int &associatedRandInt, bool &stop)
{
    if ((reader->GetFileType() != FileType::PBBAM and
         reader->GetFileType() != FileType::PBDATASET) or
        not params.concordant) {
        if (reader->GetFileType() == FileType::HDFCCS ||
            reader->GetFileType() == FileType::HDFCCSONLY) {
            if (GetNextReadThroughSemaphore(*reader, params, ccsRead, readGroupId,
                                            associatedRandInt, semaphores) == false) {
                stop = true;
                return false;
            } else {
                readIsCCS = true;
                smrtRead.Copy(ccsRead);
                ccsRead.SetQVScale(params.qvScaleType);
                smrtRead.SetQVScale(params.qvScaleType);
            }
            assert(ccsRead.zmwData.holeNumber == smrtRead.zmwData.holeNumber and
                   ccsRead.zmwData.holeNumber == ccsRead.unrolledRead.zmwData.holeNumber);
        } else {
            if (GetNextReadThroughSemaphore(*reader, params, smrtRead, readGroupId,
                                            associatedRandInt, semaphores) == false) {
                stop = true;
                return false;
            } else {
                smrtRead.SetQVScale(params.qvScaleType);
            }
        }

        //
        // Only normal (non-CCS) reads should be masked.  Since CCS reads store the raw read, that is masked.
        //
        bool readHasGoodRegion = true;
        if (params.useRegionTable and params.useHQRegionTable) {
            if (readIsCCS) {
                readHasGoodRegion =
                    MaskRead(ccsRead.unrolledRead, ccsRead.unrolledRead.zmwData, *regionTablePtr);
            } else {
                readHasGoodRegion = MaskRead(smrtRead, smrtRead.zmwData, *regionTablePtr);
            }
            //
            // Store the high quality start and end of this read for masking purposes when printing.
            //
            int hqStart, hqEnd;
            int score;
            LookupHQRegion(smrtRead.zmwData.holeNumber, *regionTablePtr, hqStart, hqEnd, score);
            smrtRead.lowQualityPrefix = hqStart;
            smrtRead.lowQualitySuffix = smrtRead.length - hqEnd;
            smrtRead.highQualityRegionScore = score;
        } else {
            smrtRead.lowQualityPrefix = 0;
            smrtRead.lowQualitySuffix = 0;
        }

        if (not IsGoodRead(smrtRead, params, stop) or stop) return false;

        return readHasGoodRegion;
    } else {
        subreads.clear();
        std::vector<SMRTSequence> reads;
        if (GetNextReadThroughSemaphore(*reader, params, reads, readGroupId, associatedRandInt,
                                        semaphores) == false) {
            stop = true;
            return false;
        }

        for (const SMRTSequence &smrtRead : reads) {
            if (IsGoodRead(smrtRead, params, stop)) {
                subreads.push_back(smrtRead);
            }
        }
        if (subreads.size() != 0) {
            smrtRead.MadeFromSubreadsAsPolymerase(subreads);
            return true;
        } else {
            return false;
        }
    }
}

void MapReadsNonCCS(MappingData<T_SuffixArray, T_GenomeSequence, T_Tuple> *mapData,
                    MappingBuffers &mappingBuffers, SMRTSequence &smrtRead,
                    SMRTSequence &smrtReadRC, std::vector<SMRTSequence> &subreads,
                    MappingParameters &params, const int &associatedRandInt,
                    ReadAlignments &allReadAlignments, std::ofstream &threadOut)
{
    DNASuffixArray sarray;
    TupleCountTable<T_GenomeSequence, DNATuple> ct;
    SequenceIndexDatabase<FASTQSequence> seqdb;
    T_GenomeSequence genome;
    BWT *bwtPtr;

    mapData->ShallowCopySuffixArray(sarray);
    mapData->ShallowCopyReferenceSequence(genome);
    mapData->ShallowCopySequenceIndexDatabase(seqdb);
    mapData->ShallowCopyTupleCountTable(ct);

    bwtPtr = mapData->bwtPtr;
    SeqBoundaryFtr<FASTQSequence> seqBoundary(&seqdb);

    std::vector<ReadInterval> subreadIntervals;
    std::vector<int> subreadDirections;
    int bestSubreadIndex;

    if ((mapData->reader->GetFileType() != FileType::PBBAM and
         mapData->reader->GetFileType() != FileType::PBDATASET) or
        not params.concordant) {
        MakePrimaryIntervals(mapData->regionTablePtr, smrtRead, subreadIntervals, subreadDirections,
                             bestSubreadIndex, params);
    } else {
        MakePrimaryIntervals(subreads, subreadIntervals, subreadDirections, bestSubreadIndex);
    }

    // Flop all directions if direction of the longest subread is 1.
    if (bestSubreadIndex >= 0 and bestSubreadIndex < int(subreadDirections.size()) and
        subreadDirections[bestSubreadIndex] == 1) {
        UpdateDirections(subreadDirections, true);
    }

    int startIndex = 0;
    int endIndex = subreadIntervals.size();

    if (params.concordant) {
        // Only the longest subread will be aligned in the first round.
        // VR , change the comment
        startIndex = std::max(startIndex, bestSubreadIndex);
        endIndex = std::min(endIndex, bestSubreadIndex + 1);

        if (params.verbosity >= 1) {
            std::cout << "Concordant template subread index: " << bestSubreadIndex << ", "
                      << smrtRead.HoleNumber() << "/" << subreadIntervals[bestSubreadIndex]
                      << std::endl;
        }
    }

    //
    // Make room for alignments.
    //
    allReadAlignments.Resize(subreadIntervals.size());
    allReadAlignments.alignMode = Subread;

    for (int intvIndex = startIndex; intvIndex < endIndex; intvIndex++) {
        SMRTSequence subreadSequence, subreadSequenceRC;
        MakeSubreadOfInterval(subreadSequence, smrtRead, subreadIntervals[intvIndex], params);
        MakeSubreadRC(subreadSequenceRC, subreadSequence, smrtRead);

        //
        // Store the sequence that is being mapped in case no hits are
        // found, and missing sequences are printed.
        //
        allReadAlignments.SetSequence(intvIndex, subreadSequence);

        std::vector<T_AlignmentCandidate *> alignmentPtrs;
        mapData->metrics.numReads++;

        assert(subreadSequence.zmwData.holeNumber == smrtRead.zmwData.holeNumber);

        //
        // Try default and fast parameters to map the read.
        //
        MapRead(subreadSequence, subreadSequenceRC,
                genome,           // possibly multi fasta file read into one sequence
                sarray, *bwtPtr,  // The suffix array, and the bwt-fm index structures
                seqBoundary,      // Boundaries of contigs in the
                // genome, alignments do not span
                // the ends of boundaries.
                ct,     // Count table to use word frequencies in the genome to weight matches.
                seqdb,  // Information about the names of
                // chromosomes in the genome, and
                // where their sequences are in the genome.
                params,  // A huge list of parameters for
                // mapping, only compile/command
                // line values set.
                mapData->metrics,  // Keep track of time/ hit counts,
                // etc.. Not fully developed, but
                // should be.
                alignmentPtrs,   // Where the results are stored.
                mappingBuffers,  // A class of buffers for structurs
                // like dyanmic programming
                // matrices, match lists, etc., that are not
                // reallocated between calls to
                // MapRead.  They are cleared though.
                mapData,  // Some values that are shared
                // across threads.
                semaphores);

        //
        // No alignments were found, sometimes parameters are
        // specified to try really hard again to find an alignment.
        // This sets some parameters that use a more sensitive search
        // at the cost of time.
        //

        if ((alignmentPtrs.size() == 0 or alignmentPtrs[0]->pctSimilarity < 80) and
            params.doSensitiveSearch) {
            MappingParameters sensitiveParams = params;
            sensitiveParams.SetForSensitivity();
            MapRead(subreadSequence, subreadSequenceRC, genome, sarray, *bwtPtr, seqBoundary, ct,
                    seqdb, sensitiveParams, mapData->metrics, alignmentPtrs, mappingBuffers,
                    mapData, semaphores);
        }

        //
        // Store the mapping quality values.
        //
        if (alignmentPtrs.size() > 0 and alignmentPtrs[0]->score < params.maxScore and
            params.storeMapQV) {
            StoreMapQVs(subreadSequence, alignmentPtrs, params);
        }

        //
        // Select alignments for this subread.
        //
        std::vector<T_AlignmentCandidate *> selectedAlignmentPtrs =
            SelectAlignmentsToPrint(alignmentPtrs, params, associatedRandInt);
        allReadAlignments.AddAlignmentsForSeq(intvIndex, selectedAlignmentPtrs);

        //
        // Move reference from subreadSequence, which will be freed at
        // the end of this loop to the smrtRead, which exists for the
        // duration of aligning all subread of the smrtRead.
        //
        for (size_t a = 0; a < alignmentPtrs.size(); a++) {
            if (alignmentPtrs[a]->qStrand == 0) {
                alignmentPtrs[a]->qAlignedSeq.ReferenceSubstring(
                    smrtRead, alignmentPtrs[a]->qAlignedSeq.seq - subreadSequence.seq,
                    alignmentPtrs[a]->qAlignedSeqLength);
            } else {
                alignmentPtrs[a]->qAlignedSeq.ReferenceSubstring(
                    smrtReadRC, alignmentPtrs[a]->qAlignedSeq.seq - subreadSequenceRC.seq,
                    alignmentPtrs[a]->qAlignedSeqLength);
            }
        }
        // Fix for memory leakage bug due to undeleted Alignment Candidate objectts which wasn't selected
        // for printing
        // delete all AC which are in complement of SelectedAlignmemntPtrs vector
        // namely (SelectedAlignmentPtrs/alignmentPtrs)
        for (size_t ii = 0; ii < alignmentPtrs.size(); ii++) {
            int found = 0;
            for (size_t jj = 0; jj < selectedAlignmentPtrs.size(); jj++) {
                if (alignmentPtrs[ii] == selectedAlignmentPtrs[jj]) {
                    found = 1;
                    break;
                }
            }
            if (found == 0) delete alignmentPtrs[ii];
        }
        subreadSequence.Free();
        subreadSequenceRC.Free();
    }  // End of looping over subread intervals within [startIndex, endIndex).

    if (params.verbosity >= 3) allReadAlignments.Print(threadOut);

    // If not concordant , all done

    if (params.concordant) {
        allReadAlignments.read = smrtRead;
        allReadAlignments.alignMode = ZmwSubreads;

        if (startIndex >= 0 && startIndex < int(allReadAlignments.subreadAlignments.size())) {
            std::vector<T_AlignmentCandidate *> selectedAlignmentPtrs =
                allReadAlignments.CopySubreadAlignments(startIndex);

            for (int alignmentIndex = 0; alignmentIndex < int(selectedAlignmentPtrs.size());
                 alignmentIndex++) {
                FlankTAlignedSeq(selectedAlignmentPtrs[alignmentIndex], seqdb, genome,
                                 params.flankSize);
            }

            for (int intvIndex = 0; intvIndex < int(subreadIntervals.size()); intvIndex++) {
                if (intvIndex == startIndex) continue;
                int passDirection = subreadDirections[intvIndex];
                int passStartBase = subreadIntervals[intvIndex].start;
                int passNumBases = subreadIntervals[intvIndex].end - passStartBase;
                if (passNumBases <= params.minReadLength) {
                    continue;
                }

                mapData->metrics.numReads++;
                SMRTSequence subread;
                subread.ReferenceSubstring(smrtRead, passStartBase, passNumBases);
                subread.CopyTitle(smrtRead.title);
                // The unrolled alignment should be relative to the entire read.
                if (params.clipping == SAMOutput::subread) {
                    SMRTSequence maskedSubread;
                    MakeSubreadOfInterval(maskedSubread, smrtRead, subreadIntervals[intvIndex],
                                          params);
                    allReadAlignments.SetSequence(intvIndex, maskedSubread);
                    maskedSubread.Free();
                } else {
                    allReadAlignments.SetSequence(intvIndex, smrtRead);
                }

                for (size_t alnIndex = 0; alnIndex < selectedAlignmentPtrs.size(); alnIndex++) {
                    T_AlignmentCandidate *alignment = selectedAlignmentPtrs[alnIndex];
                    if (alignment->score > params.maxScore) break;
                    AlignSubreadToAlignmentTarget(allReadAlignments, subread, smrtRead, alignment,
                                                  passDirection, subreadIntervals[intvIndex],
                                                  intvIndex, params, mappingBuffers, threadOut);
                    if (params.concordantAlignBothDirections) {
                        AlignSubreadToAlignmentTarget(allReadAlignments, subread, smrtRead,
                                                      alignment, ((passDirection == 0) ? 1 : 0),
                                                      subreadIntervals[intvIndex], intvIndex,
                                                      params, mappingBuffers, threadOut);
                    }
                }  // End of aligning this subread to each selected alignment.
                subread.Free();
            }  // End of aligning each subread to where the template subread aligned to.
            for (size_t alignmentIndex = 0; alignmentIndex < selectedAlignmentPtrs.size();
                 alignmentIndex++) {
                if (selectedAlignmentPtrs[alignmentIndex])
                    delete selectedAlignmentPtrs[alignmentIndex];
            }
        }  // End of if startIndex >= 0 and < subreadAlignments.size()
    }      // End of if params.concordant
}

//
// invoked for mapping entire ZMW as a single entity
// either for CCS reads : all subreads of a ZMW collapsed/merged into a single read
// or Polymerase reads  : all subreads of a ZMW stitched into a single read
//
void MapReadsCCS(MappingData<T_SuffixArray, T_GenomeSequence, T_Tuple> *mapData,
                 MappingBuffers &mappingBuffers, SMRTSequence &smrtRead, SMRTSequence &smrtReadRC,
                 CCSSequence &ccsRead, const bool readIsCCS, MappingParameters &params,
                 const int &associatedRandInt, ReadAlignments &allReadAlignments,
                 std::ofstream &threadOut)
{
    DNASuffixArray sarray;
    TupleCountTable<T_GenomeSequence, DNATuple> ct;
    SequenceIndexDatabase<FASTQSequence> seqdb;
    T_GenomeSequence genome;
    BWT *bwtPtr;

    mapData->ShallowCopySuffixArray(sarray);
    mapData->ShallowCopyReferenceSequence(genome);
    mapData->ShallowCopySequenceIndexDatabase(seqdb);
    mapData->ShallowCopyTupleCountTable(ct);

    bwtPtr = mapData->bwtPtr;
    SeqBoundaryFtr<FASTQSequence> seqBoundary(&seqdb);

    //
    // The read must be mapped as a whole, even if it contains subreads.
    //
    std::vector<T_AlignmentCandidate *> alignmentPtrs;
    mapData->metrics.numReads++;
    smrtRead.SubreadStart(0).SubreadEnd(smrtRead.length);
    smrtReadRC.SubreadStart(0).SubreadEnd(smrtRead.length);

    MapRead(smrtRead, smrtReadRC, genome, sarray, *bwtPtr, seqBoundary, ct, seqdb, params,
            mapData->metrics, alignmentPtrs, mappingBuffers, mapData, semaphores);

    //
    // Store the mapping quality values.
    //
    if (alignmentPtrs.size() > 0 and alignmentPtrs[0]->score < params.maxScore and
        params.storeMapQV) {
        StoreMapQVs(smrtRead, alignmentPtrs, params);
    }

    //
    // Select de novo ccs-reference alignments for subreads to align to.
    //
    std::vector<T_AlignmentCandidate *> selectedAlignmentPtrs =
        SelectAlignmentsToPrint(alignmentPtrs, params, associatedRandInt);

    //
    // Just one sequence is aligned.  There is one primary hit, and
    // all other are secondary.
    //

    //
    // Here unrolled reads are aligned
    //
    if (readIsCCS == false or params.useCcsOnly) {
        // if -noSplitSubreads or -useccsdenovo.
        //
        // Record some information for proper SAM Annotation.
        //
        allReadAlignments.Resize(1);
        allReadAlignments.AddAlignmentsForSeq(0, selectedAlignmentPtrs);
        if (params.useCcsOnly) {
            allReadAlignments.alignMode = CCSDeNovo;
        } else {
            allReadAlignments.alignMode = Fullread;
        }
        allReadAlignments.SetSequence(0, smrtRead);
    }
    //
    // Here CCS reads are aligned
    //
    else if (readIsCCS) {  // if -useccsall or -useccs
        // Flank alignment candidates to both ends.
        for (size_t alignmentIndex = 0; alignmentIndex < selectedAlignmentPtrs.size();
             alignmentIndex++) {
            FlankTAlignedSeq(selectedAlignmentPtrs[alignmentIndex], seqdb, genome,
                             params.flankSize);
        }

        //
        // Align the ccs subread to where the denovo sequence mapped (explode).
        //
        CCSIterator ccsIterator;
        FragmentCCSIterator fragmentCCSIterator;
        CCSIterator *subreadIterator;

        //
        // Choose a different iterator over subreads depending on the
        // alignment mode.  When the mode is allpass, include the
        // framgents that are not necessarily full pass.
        //
        if (params.useAllSubreadsInCcs) {
            //
            // Use all subreads even if they are not full pass
            fragmentCCSIterator.Initialize(&ccsRead, mapData->regionTablePtr);
            subreadIterator = &fragmentCCSIterator;
            allReadAlignments.alignMode = CCSAllPass;
        } else {
            // Use only full pass reads.
            ccsIterator.Initialize(&ccsRead);
            subreadIterator = &ccsIterator;
            allReadAlignments.alignMode = CCSFullPass;
        }

        allReadAlignments.Resize(subreadIterator->GetNumPasses());

        int passDirection, passStartBase, passNumBases;
        SMRTSequence subread;

        //
        // The read was previously set to the smrtRead, which was the
        // de novo ccs sequence.  Since the alignments of exploded
        // reads are reported, the unrolled read should be used as the
        // reference when printing.
        //
        allReadAlignments.read = ccsRead.unrolledRead;
        subreadIterator->Reset();
        int subreadIndex;

        //
        // Realign all subreads to selected reference locations.
        //
        for (subreadIndex = 0; subreadIndex < subreadIterator->GetNumPasses(); subreadIndex++) {
            int retval = subreadIterator->GetNext(passDirection, passStartBase, passNumBases);
            assert(retval == 1);
            if (passNumBases <= params.minReadLength) {
                continue;
            }

            ReadInterval subreadInterval(passStartBase, passStartBase + passNumBases);

            subread.ReferenceSubstring(ccsRead.unrolledRead, passStartBase, passNumBases - 1);
            subread.CopyTitle(ccsRead.title);
            // The unrolled alignment should be relative to the entire read.
            allReadAlignments.SetSequence(subreadIndex, ccsRead.unrolledRead);

            //
            // Align this subread to all the positions that the de novo
            // sequence has aligned to.
            //
            for (size_t alignmentIndex = 0; alignmentIndex < selectedAlignmentPtrs.size();
                 alignmentIndex++) {
                T_AlignmentCandidate *alignment = selectedAlignmentPtrs[alignmentIndex];
                if (alignment->score > params.maxScore) break;
                AlignSubreadToAlignmentTarget(allReadAlignments, subread, ccsRead.unrolledRead,
                                              alignment, passDirection, subreadInterval,
                                              subreadIndex, params, mappingBuffers, threadOut);
            }  // End of aligning this subread to where the de novo ccs has aligned to.
            subread.Free();
        }  // End of alignining all subreads to where the de novo ccs has aligned to.
    }      // End of if readIsCCS and !params.useCcsOnly

    // Fix for memory leakage due to undeleted Alignment Candidate objectts not selected
    // for printing
    // delete all AC which are in complement of SelectedAlignmemntPtrs vector
    // namely (SelectedAlignmentPtrs/alignmentPtrs)
    for (size_t ii = 0; ii < alignmentPtrs.size(); ii++) {
        int found = 0;
        for (size_t jj = 0; jj < selectedAlignmentPtrs.size(); jj++) {
            if (alignmentPtrs[ii] == selectedAlignmentPtrs[jj]) {
                found = 1;
                break;
            }
        }
        if (found == 0) delete alignmentPtrs[ii];
    }
}

void MapReads(MappingData<T_SuffixArray, T_GenomeSequence, T_Tuple> *mapData)
{
    //
    // Step 1, initialize local pointers to map data
    // for programming shorthand.
    //
    MappingParameters params = mapData->params;

    DNASuffixArray sarray;
    TupleCountTable<T_GenomeSequence, DNATuple> ct;
    SequenceIndexDatabase<FASTQSequence> seqdb;
    T_GenomeSequence genome;

    mapData->ShallowCopySuffixArray(sarray);
    mapData->ShallowCopyReferenceSequence(genome);
    mapData->ShallowCopySequenceIndexDatabase(seqdb);
    mapData->ShallowCopyTupleCountTable(ct);

    SeqBoundaryFtr<FASTQSequence> seqBoundary(&seqdb);

    int numAligned = 0;

    SMRTSequence smrtRead, smrtReadRC;
    SMRTSequence unrolledReadRC;
    CCSSequence ccsRead;

    // Print verbose logging to pid.threadid.log for each thread.
    std::ofstream threadOut;
    if (params.verbosity >= 3) {
        std::stringstream ss;
        ss << getpid() << "." << pthread_self();
        std::string threadLogFileName = ss.str() + ".log";
        threadOut.open(threadLogFileName.c_str(), std::ios::out | std::ios::app);
    }

    //
    // Reuse the following buffers during alignment.  Since these keep
    // storage contiguous, hopefully this will decrease memory
    // fragmentation.
    //
    MappingBuffers mappingBuffers;
    while (true) {
        // Fetch reads from a zmw
        bool readIsCCS = false;
        AlignmentContext alignmentContext;
        // Associate each sequence to read in with a determined random int.
        int associatedRandInt = 0;
        bool stop = false;
        std::vector<SMRTSequence> subreads;
        bool readsOK =
            FetchReads(mapData->reader, mapData->regionTablePtr, smrtRead, ccsRead, subreads,
                       params, readIsCCS, alignmentContext.readGroupId, associatedRandInt, stop);
        if (stop) break;
        if (not readsOK) continue;

        if (params.verbosity > 1) {
            std::cout << "aligning read: " << std::endl;
            smrtRead.PrintSeq(std::cout);
        }

        smrtRead.MakeRC(smrtReadRC);

        // important
        // 1. CCS and unrolled mode are mutually exclusive
        // 2. Reverse Complement Read is generated fort CCS only
        //
        if (readIsCCS) {
            ccsRead.unrolledRead.MakeRC(unrolledReadRC);
        }

        //
        // When aligning subreads separately, iterate over each subread, and
        // print the alignments for these.
        //
        ReadAlignments allReadAlignments;
        allReadAlignments.read = smrtRead;

        // currently 3 ways of mapping
        // regular, CCS , and Polymerase (unrolled)
        //
        // for regular subreads MapReadsNonCCS
        // for mapping ZMW as a whole (CCS or Polymerase) MapReadsCCS
        // For the future , change the name of functions  to be more desriptive
        // noSplitSubreads is in essense unrolled - Polymerase read mode
        //
        if (readIsCCS == false and params.mapSubreadsSeparately) {
            // (not readIsCCS and not -noSplitSubreads)
            MapReadsNonCCS(mapData, mappingBuffers, smrtRead, smrtReadRC, subreads, params,
                           associatedRandInt, allReadAlignments, threadOut);
        }       // End of if (readIsCCS == false and params.mapSubreadsSeparately).
        else {  // if (readIsCCS or (not readIsCCS and -noSplitSubreads) )
            MapReadsCCS(mapData, mappingBuffers, smrtRead, smrtReadRC, ccsRead, readIsCCS, params,
                        associatedRandInt, allReadAlignments, threadOut);
        }  // End of if not (readIsCCS == false and params.mapSubreadsSeparately)

        PrintAllReadAlignments(allReadAlignments, alignmentContext, *mapData->outFilePtr,
                               *mapData->unalignedFilePtr, params, subreads,
#ifdef USE_PBBAM
                               bamWriterPtr,
#endif
                               semaphores);

        allReadAlignments.Clear();
        smrtReadRC.Free();
        smrtRead.Free();

        if (readIsCCS) {
            ccsRead.Free();
            unrolledReadRC.Free();
        }
        numAligned++;
        if (numAligned % 100 == 0) {
            mappingBuffers.Reset();
        }
    }  // End of while (true).
    smrtRead.Free();
    smrtReadRC.Free();
    unrolledReadRC.Free();
    ccsRead.Free();

    if (params.nProc > 1) {
#ifdef __APPLE__
        sem_wait(semaphores.reader);
        sem_post(semaphores.reader);
#else
        sem_wait(&semaphores.reader);
        sem_post(&semaphores.reader);
#endif
    }
}

int main(int argc, char *argv[])
{
    //
    // Configure parameters for refining alignments.
    //
    MappingParameters params;

    CommandLineParser clp;
    clp.SetHelp(BlasrHelp(params));
    clp.SetConciseHelp(BlasrConciseHelp());
    clp.SetProgramSummary(BlasrSummaryHelp());
    clp.SetProgramName("blasr");
    clp.SetVersion(GetVersion());

    // Register Blasr options.
    RegisterBlasrOptions(clp, params);

    // Parse command line args.
    clp.ParseCommandLine(argc, argv, params.readsFileNames);

    std::string commandLine;
    clp.CommandLineToString(argc, argv, commandLine);

    if (params.printVerboseHelp) {
        std::cout << BlasrHelp(params) << std::endl;
        std::exit(EXIT_SUCCESS);  // Not a failure.
    }
    if (argc < 3) {
        std::cout << BlasrConciseHelp();
        std::exit(EXIT_FAILURE);  // A failure.
    }

    int a, b;
    for (a = 0; a < 5; a++) {
        for (b = 0; b < 5; b++) {
            if (a != b) {
                SMRTDistanceMatrix[a][b] += params.mismatch;
            } else {
                SMRTDistanceMatrix[a][b] += params.match;
            }
        }
    }

    if (params.scoreMatrixString != "") {
        if (StringToScoreMatrix(params.scoreMatrixString, SMRTDistanceMatrix) == false) {
            std::cout << "ERROR. The string " << std::endl
                      << params.scoreMatrixString << std::endl
                      << "is not a valid format.  It should be a quoted, space separated string of "
                      << std::endl
                      << "integer values.  The matrix: " << std::endl
                      << "    A  C  G  T  N" << std::endl
                      << " A  1  2  3  4  5" << std::endl
                      << " C  6  7  8  9 10" << std::endl
                      << " G 11 12 13 14 15" << std::endl
                      << " T 16 17 18 19 20" << std::endl
                      << " N 21 22 23 24 25" << std::endl
                      << " should be specified as \"1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 "
                         "19 20 21 22 23 24 25\""
                      << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }

    std::cerr << "[INFO] " << GetTimestamp() << " [blasr] started." << std::endl;
    params.MakeSane();

    //
    // The random number generator is used for subsampling for debugging
    // and testing consensus and selecting hits when hit policy is random
    // or randombest.
    //
    if (params.useRandomSeed == true) {
        InitializeRandomGenerator(params.randomSeed);
    } else {
        InitializeRandomGeneratorWithTime();
    }

    //
    // Various aspects of timing are stored here.  However this isn't
    // quite finished.
    //
    MappingMetrics metrics;

    std::ofstream fullMetricsFile;
    if (params.fullMetricsFileName != "") {
        CrucialOpen(params.fullMetricsFileName, fullMetricsFile, std::ios::out);
        metrics.SetStoreList();
    }

    //
    // If reading a separate region table, there is a 1-1 correspondence
    // between region table and bas file.
    //
    if (params.readSeparateRegionTable) {
        if (FileOfFileNames::IsFOFN(params.regionTableFileName)) {
            FileOfFileNames::FOFNToList(params.regionTableFileName, params.regionTableFileNames);
        } else {
            params.regionTableFileNames.push_back(params.regionTableFileName);
        }
    }

    if (params.regionTableFileNames.size() != 0 and
        params.regionTableFileNames.size() != params.queryFileNames.size()) {
        std::cout << "Error, there are not the same number of region table files as input files."
                  << std::endl;
        std::exit(EXIT_FAILURE);
    }

    // If reading a separate ccs fofn, there is a 1-1 corresponence
    // between ccs fofn and base file.
    if (params.readSeparateCcsFofn) {
        if (FileOfFileNames::IsFOFN(params.ccsFofnFileName)) {
            FileOfFileNames::FOFNToList(params.ccsFofnFileName, params.ccsFofnFileNames);
        } else {
            params.ccsFofnFileNames.push_back(params.ccsFofnFileName);
        }
    }
    if (params.ccsFofnFileNames.size() != 0 and
        params.ccsFofnFileNames.size() != params.queryFileNames.size()) {
        std::cout << "Error, there are not the same number of ccs files as input files."
                  << std::endl;
        std::exit(EXIT_FAILURE);
    }

    SequenceIndexDatabase<FASTASequence> seqdb;
    SeqBoundaryFtr<FASTASequence> seqBoundary(&seqdb);

    //
    // Initialize the sequence index database if it used. If it is not
    // specified, it is initialized by default when reading a multiFASTA
    // file.
    //
    if (params.useSeqDB) {
        std::ifstream seqdbin;
        CrucialOpen(params.seqDBName, seqdbin);
        seqdb.ReadDatabase(seqdbin);
    }

    //
    // Make sure the reads file exists and can be opened before
    // trying to read any of the larger data structures.
    //

    FASTASequence fastaGenome;
    T_Sequence genome;
    FASTAReader genomeReader;

    //
    // The genome is in normal FASTA, or condensed (lossy homopolymer->unipolymer)
    // format.  Both may be read in using a FASTA reader.
    //
    if (!genomeReader.Init(params.genomeFileName)) {
        std::cout << "Could not open genome file " << params.genomeFileName << std::endl;
        std::exit(EXIT_FAILURE);
    }

    if (params.printSAM or params.printBAM) {
        genomeReader.computeMD5 = true;
    }
    //
    // If no sequence title database is supplied, initialize one when
    // reading in the reference, and consider a seqdb to be present.
    //
    if (!params.useSeqDB) {
        genomeReader.ReadAllSequencesIntoOne(fastaGenome, &seqdb);
        params.useSeqDB = true;
    } else {
        genomeReader.ReadAllSequencesIntoOne(fastaGenome);
    }
    genomeReader.Close();
    //
    // The genome may have extra spaces in the fasta name. Get rid of those.
    //
    for (int t = 0; t < fastaGenome.titleLength; t++) {
        if (fastaGenome.title[t] == ' ') {
            fastaGenome.titleLength = t;
            fastaGenome.title[t] = '\0';
            break;
        }
    }

    genome.seq = fastaGenome.seq;
    genome.length = fastaGenome.length;
    genome.title = fastaGenome.title;
    genome.deleteOnExit = false;
    genome.titleLength = fastaGenome.titleLength;
    genome.ToUpper();

    DNASuffixArray sarray;
    TupleCountTable<T_GenomeSequence, DNATuple> ct;

    std::ofstream outFile;
    outFile.exceptions(std::ostream::failbit);
    std::ofstream unalignedOutFile;
    BWT bwt;

    if (params.useBwt) {
        if (bwt.Read(params.bwtFileName) == 0) {
            std::cout << "ERROR! Could not read the BWT file. " << params.bwtFileName << std::endl;
            std::exit(EXIT_FAILURE);
        }
    } else {
        if (!params.useSuffixArray) {
            //
            // There was no explicit specification of a suffix
            // array on the command line, so build it on the fly here.
            //
            genome.ToThreeBit();
            std::vector<int> alphabet;
            sarray.InitThreeBitDNAAlphabet(alphabet);
            sarray.LarssonBuildSuffixArray(genome.seq, genome.length, alphabet);
            if (params.minMatchLength > 0) {
                if (params.anchorParameters.useLookupTable == true) {
                    if (params.lookupTableLength > params.minMatchLength) {
                        params.lookupTableLength = params.minMatchLength;
                    }
                    sarray.BuildLookupTable(genome.seq, genome.length, params.lookupTableLength);
                }
            }
            genome.ConvertThreeBitToAscii();
            params.useSuffixArray = 1;
        } else if (params.useSuffixArray) {
            if (sarray.Read(params.suffixArrayFileName)) {
                if (params.minMatchLength != 0) {
                    params.listTupleSize = std::min(8, params.minMatchLength);
                } else {
                    params.listTupleSize = sarray.lookupPrefixLength;
                }
                if (params.minMatchLength < int(sarray.lookupPrefixLength)) {
                    std::cerr << "WARNING. The value of -minMatch " << params.minMatchLength
                              << " is less than the smallest searched length of "
                              << sarray.lookupPrefixLength << ".  Setting -minMatch to "
                              << sarray.lookupPrefixLength << "." << std::endl;
                    params.minMatchLength = sarray.lookupPrefixLength;
                }
            } else {
                std::cout << "ERROR. " << params.suffixArrayFileName
                          << " is not a valid suffix array. " << std::endl
                          << " Make sure it is generated with the latest version of sawriter."
                          << std::endl;
                std::exit(EXIT_FAILURE);
            }
        }
    }

    if (params.minMatchLength < int(sarray.lookupPrefixLength)) {
        std::cerr << "WARNING. The value of -minMatch " << params.minMatchLength
                  << " is less than the smallest searched length of " << sarray.lookupPrefixLength
                  << ".  Setting -minMatch to " << sarray.lookupPrefixLength << "." << std::endl;
        params.minMatchLength = sarray.lookupPrefixLength;
    }

    //
    // It is required to have a tuple count table
    // for estimating the background frequencies
    // for word matching.
    // If one is specified on the command line, simply read
    // it in.  If not, this is operating under the mode
    // that everything is computed from scratch.
    //
    TupleMetrics saLookupTupleMetrics;
    if (params.useCountTable) {
        std::ifstream ctIn;
        CrucialOpen(params.countTableName, ctIn, std::ios::in | std::ios::binary);
        ct.Read(ctIn);
        saLookupTupleMetrics = ct.tm;

    } else {
        saLookupTupleMetrics.Initialize(params.lookupTableLength);
        ct.InitCountTable(saLookupTupleMetrics);
        ct.AddSequenceTupleCountsLR(genome);
    }

    TitleTable titleTable;
    if (params.useTitleTable) {
        std::ofstream titleTableOut;
        CrucialOpen(params.titleTableName, titleTableOut);
        //
        // When using a sequence index database, the title table is simply copied
        // from the sequencedb.
        //
        if (params.useSeqDB) {
            titleTable.Copy(seqdb.names, seqdb.nSeqPos - 1);
            titleTable.ResetTableToIntegers(seqdb.names, seqdb.nameLengths, seqdb.nSeqPos - 1);
        } else {
            //
            // No seqdb, so there is just one sequence. Still the user specified a title
            // table, so just the first sequence in the fasta file should be used.
            //
            titleTable.Copy(&fastaGenome.title, 1);
            titleTable.ResetTableToIntegers(&genome.title, &genome.titleLength, 1);
            fastaGenome.titleLength = strlen(genome.title);
        }
        titleTable.Write(titleTableOut);
    } else {
        if (params.useSeqDB) {
            //
            // When using a sequence index database, but not the titleTable,
            // it is necessary to truncate the titles at the first space to
            // be compatible with the way other alignment programs interpret
            // fasta titles.  When printing the title table, there is all
            // sorts of extra storage space, so the full line is stored.
            //
            seqdb.SequenceTitleLinesToNames();
        }
    }

    std::ostream *outFilePtr = &std::cout;
    std::ofstream outFileStrm;
    std::ofstream unalignedFile;
    std::ostream *unalignedFilePtr = NULL;
    std::ofstream metricsOut, lcpBoundsOut;
    std::ofstream anchorFileStrm;
    std::ofstream clusterOut, *clusterOutPtr;

    if (params.anchorFileName != "") {
        CrucialOpen(params.anchorFileName, anchorFileStrm, std::ios::out);
    }

    if (params.clusterFileName != "") {
        CrucialOpen(params.clusterFileName, clusterOut, std::ios::out);
        clusterOutPtr = &clusterOut;
        clusterOut << "total_size p_value n_anchors read_length align_score read_accuracy "
                      "anchor_probability min_exp_anchors seq_length"
                   << std::endl;
    } else {
        clusterOutPtr = NULL;
    }

    if (params.outFileName != "") {
        if (not params.printBAM) {
            CrucialOpen(params.outFileName, outFileStrm, std::ios::out);
            outFilePtr = &outFileStrm;
        }  // otherwise, use bamWriter and initialize it later
    }

    if (params.printHeader) {
        switch (params.printFormat) {
            case (SummaryPrint):
                SummaryOutput::PrintHeader(*outFilePtr);
                break;
            case (Interval):
                IntervalOutput::PrintHeader(*outFilePtr);
                break;
            case (CompareSequencesParsable):
                CompareSequencesOutput::PrintHeader(*outFilePtr);
                break;
        }
    }

    if (params.printUnaligned == true) {
        CrucialOpen(params.unalignedFileName, unalignedFile, std::ios::out);
        unalignedFilePtr = &unalignedFile;
    }

    if (params.metricsFileName != "") {
        CrucialOpen(params.metricsFileName, metricsOut);
    }

    if (params.lcpBoundsFileName != "") {
        CrucialOpen(params.lcpBoundsFileName, lcpBoundsOut);
        //    lcpBoundsOut << "pos depth width lnwidth" << std::endl;
    }

    //
    // Configure the mapping database.
    //

    MappingData<T_SuffixArray, T_GenomeSequence, T_Tuple> *mapdb =
        new MappingData<T_SuffixArray, T_GenomeSequence, T_Tuple>[ params.nProc ];

    int procIndex;
    pthread_attr_t *threadAttr = new pthread_attr_t[params.nProc];
    //  MappingSemaphores semaphores;
    //
    // When there are multiple processes running along, sometimes there
    // are semaphores to worry about.
    //

    if (params.nProc > 1) {
        semaphores.InitializeAll();
    }
    for (procIndex = 0; procIndex < params.nProc; procIndex++) {
        pthread_attr_init(&threadAttr[procIndex]);
    }

    //
    // Start the mapping jobs.
    //
    if (params.subsample < 1) {
        InitializeRandomGeneratorWithTime();
        reader = new ReaderAgglomerate(params.subsample);
    } else {
        reader = new ReaderAgglomerate(params.startRead, params.stride);
    }
    //  In case the input is fasta, make all bases in upper case.
    reader->SetToUpper();

    regionTableReader = new HDFRegionTableReader;
    RegionTable regionTable;
    //
    // Store lists of how long it took to map each read.
    //
    metrics.clocks.SetStoreList(true);
    if (params.useCcs) {
        reader->UseCCS();
    }

    std::string commandLineString;  // Restore command.
    clp.CommandLineToString(argc, argv, commandLineString);

    if (params.printSAM or params.printBAM) {
        std::string so = "UNKNOWN";          // sorting order;
        std::string version = GetVersion();  //blasr version;
        SAMHeaderPrinter shp(so, seqdb, params.queryFileNames, params.queryReadType,
                             params.samQVList, "BLASR", version, commandLineString);
        std::string headerString = shp.ToString();  // SAM/BAM header
        if (params.printSAM) {
            // this is not going to be executed since sam is printed via bam
            *outFilePtr << headerString;
        } else if (params.printBAM) {
// here both bam and sam are handled
#ifdef USE_PBBAM
            PacBio::BAM::BamHeader header = PacBio::BAM::BamHeader(headerString);
            // Create bam header
            // Both file name and SAMHeader are required in order to create a BamWriter.
            // sam_via_bam changes
            if (params.sam_via_bam) {
                bamWriterPtr = new PacBio::BAM::SamWriter(params.outFileName, header);
            } else {
                bamWriterPtr = new PacBio::BAM::BamWriter(params.outFileName, header);
            }
#else
            REQUIRE_PBBAM_ERROR();
#endif
        }
    }

    for (size_t readsFileIndex = 0; readsFileIndex < params.queryFileNames.size();
         readsFileIndex++) {
        params.readsFileIndex = readsFileIndex;
        //
        // Configure the reader to use the correct read and region
        // file names.
        //
        reader->SetReadFileName(params.queryFileNames[params.readsFileIndex]);

        // if PBBAM , need to construct scrap file name and check if exist
        //
        // Initialize using already set file names.
        //

        // unrolled Need to pass unrolled option
        // unrolled If not PBDATASET also need to construct scrap file name and
        // test if it exists in the same directory, if not exit with error message
        //
        int initReturnValue;

        if (((reader->GetFileType() == FileType::PBDATASET) ||
             (reader->GetFileType() == FileType::PBBAM)) and
            not params.mapSubreadsSeparately) {

            if (reader->GetFileType() == FileType::PBBAM) {
                reader->SetScrapsFileName(params.scrapsFileNames[params.readsFileIndex]);
            }
            initReturnValue = reader->Initialize(true, params.polymeraseMode);
        } else {
            initReturnValue = reader->Initialize();
        }
        if (initReturnValue <= 0) {
            std::cerr << "WARNING! Could not open file "
                      << params.queryFileNames[params.readsFileIndex] << std::endl;
            continue;
        }

        // Check whether use ccs only.
        if (reader->GetFileType() == FileType::HDFCCSONLY) {
            params.useAllSubreadsInCcs = false;
            params.useCcs = params.useCcsOnly = true;
        }

        std::string changeListIdString;
        reader->hdfBasReader.GetChangeListID(changeListIdString);
        ChangeListID changeListId(changeListIdString);
        params.qvScaleType = DetermineQVScaleFromChangeListID(changeListId);
        if (reader->FileHasZMWInformation() and params.useRegionTable) {
            if (params.readSeparateRegionTable) {
                if (regionTableReader->Initialize(
                        params.regionTableFileNames[params.readsFileIndex]) == 0) {
                    std::cout << "ERROR! Could not read the region table "
                              << params.regionTableFileNames[params.readsFileIndex] << std::endl;
                    std::exit(EXIT_FAILURE);
                }
                params.useRegionTable = true;
            } else {
                if (reader->HasRegionTable()) {
                    if (regionTableReader->Initialize(
                            params.queryFileNames[params.readsFileIndex]) == 0) {
                        std::cout << "ERROR! Could not read the region table "
                                  << params.queryFileNames[params.readsFileIndex] << std::endl;
                        std::exit(EXIT_FAILURE);
                    }
                    params.useRegionTable = true;
                } else {
                    params.useRegionTable = false;
                }
            }
        } else {
            params.useRegionTable = false;
        }

        //
        //  Check to see if there is a region table. If there is a separate
        //  region table, use that (over the region table in the bas
        // file).  If there is a region table in the bas file, use that,
        // without having to specify a region table on the command line.
        //
        if (params.useRegionTable) {
            regionTable.Reset();
            regionTableReader->ReadTable(regionTable);
            regionTableReader->Close();
        }

        //
        // Check to see if there is a separate ccs fofn. If there is a separate
        // ccs fofn, use that over the one in the bas file.
        //
        //if (params.readSeparateCcsFofn and params.useCcs) {
        //  if (reader->SetCCS(params.ccsFofnFileNames[params.readsFileIndex]) == 0) {
        //    std::cout << "ERROR! Could not read the ccs file "
        //         << params.ccsFofnFileNames[params.readsFileIndex] << std::endl;
        //    std::exit(EXIT_FAILURE);
        //  }
        // }

        if (reader->GetFileType() != FileType::HDFCCS and
            reader->GetFileType() != FileType::HDFBase and
            reader->GetFileType() != FileType::HDFPulse and
            reader->GetFileType() != FileType::PBBAM and
            reader->GetFileType() != FileType::PBDATASET and params.concordant) {
            std::cerr << "WARNING! Option concordant is only enabled when "
                      << "input reads are in PacBio bax/pls.h5, bam or "
                      << "dataset xml format." << std::endl;
            params.concordant = false;
        }

#ifdef USE_GOOGLE_PROFILER
        char *profileFileName = getenv("CPUPROFILE");
        if (profileFileName != NULL) {
            ProfilerStart(profileFileName);
        } else {
            ProfilerStart("google_profile.txt");
        }
#endif

        assert(initReturnValue > 0);
        if (params.nProc == 1) {
            mapdb[0].Initialize(&sarray, &genome, &seqdb, &ct, params, reader, &regionTable,
                                outFilePtr, unalignedFilePtr, &anchorFileStrm, clusterOutPtr);
            mapdb[0].bwtPtr = &bwt;
            if (params.fullMetricsFileName != "") {
                mapdb[0].metrics.SetStoreList(true);
            }
            if (params.lcpBoundsFileName != "") {
                mapdb[0].lcpBoundsOutPtr = &lcpBoundsOut;
            } else {
                mapdb[0].lcpBoundsOutPtr = NULL;
            }

            MapReads(&mapdb[0]);
            metrics.Collect(mapdb[0].metrics);
        } else {
            pthread_t *threads = new pthread_t[params.nProc];
            for (procIndex = 0; procIndex < params.nProc; procIndex++) {
                //
                // Initialize thread-specific parameters.
                //

                mapdb[procIndex].Initialize(&sarray, &genome, &seqdb, &ct, params, reader,
                                            &regionTable, outFilePtr, unalignedFilePtr,
                                            &anchorFileStrm, clusterOutPtr);
                mapdb[procIndex].bwtPtr = &bwt;
                if (params.fullMetricsFileName != "") {
                    mapdb[procIndex].metrics.SetStoreList(true);
                }
                if (params.lcpBoundsFileName != "") {
                    mapdb[procIndex].lcpBoundsOutPtr = &lcpBoundsOut;
                } else {
                    mapdb[procIndex].lcpBoundsOutPtr = NULL;
                }

                if (params.outputByThread) {
                    std::ofstream *outPtr = new std::ofstream;
                    mapdb[procIndex].outFilePtr = outPtr;
                    std::stringstream outNameStream;
                    outNameStream << params.outFileName << "." << procIndex;
                    mapdb[procIndex].params.outFileName = outNameStream.str();
                    CrucialOpen(mapdb[procIndex].params.outFileName, *outPtr, std::ios::out);
                }
                pthread_create(&threads[procIndex], &threadAttr[procIndex],
                               (void *(*)(void *))MapReads, &mapdb[procIndex]);
            }
            for (procIndex = 0; procIndex < params.nProc; procIndex++) {
                pthread_join(threads[procIndex], NULL);
            }
            for (procIndex = 0; procIndex < params.nProc; procIndex++) {
                metrics.Collect(mapdb[procIndex].metrics);
                if (params.outputByThread) {
                    delete mapdb[procIndex].outFilePtr;
                }
            }
            if (threads) {
                delete[] threads;
                threads = NULL;
            }
        }
        reader->Close();
    }

    if (!reader) {
        delete reader;
        reader = NULL;
    }

    fastaGenome.Free();
#ifdef USE_GOOGLE_PROFILER
    ProfilerStop();
#endif

    if (mapdb != NULL) {
        delete[] mapdb;
    }
    if (threadAttr != NULL) {
        delete[] threadAttr;
    }
    seqdb.FreeDatabase();
    if (regionTableReader) {
        delete regionTableReader;
    }
    if (params.metricsFileName != "") {
        metrics.PrintSummary(metricsOut);
    }
    if (params.fullMetricsFileName != "") {
        metrics.PrintFullList(fullMetricsFile);
    }
    if (params.outFileName != "") {
        if (params.printBAM) {
#ifdef USE_PBBAM
            assert(bamWriterPtr);
            try {
                if (!params.sam_via_bam) {  // no need to flush for SAM , but need to understand why
                    bamWriterPtr->TryFlush();
                }
                delete bamWriterPtr;
                bamWriterPtr = NULL;
            } catch (const std::exception &e) {
                std::cout << "Error, could not flush bam records to bam file." << std::endl;
                std::exit(EXIT_FAILURE);
            }
#else
            REQUIRE_PBBAM_ERROR();
#endif
        } else {
            outFileStrm.close();
        }
    }
    std::cerr << "[INFO] " << GetTimestamp() << " [blasr] ended." << std::endl;
    return 0;
}
