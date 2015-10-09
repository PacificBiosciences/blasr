/*
 * =====================================================================================
 *
 *       Filename:  SAMFilter.cpp
 *
 *    Description:  Filter SAM Hits according to 
 *                  filteration criteria
 *                     minPctSimilarity, minAccuracy, 
 *                     minLength, holeNumbers
 *                  and multiple-hit policy
 *                     random    : a random hit
 *                     all       : all hits
 *                     allbest   : all hits with the best score
 *                     randombest: a random hit selected from all the hits 
 *                                 that have the best score
 *
 *        Version:  1.0
 *        Created:  03/19/2013 01:19:43 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Yuan Li (yli), yli@pacificbiosciences.com
 *        Company:  Pacific Biosciences
 *
 * =====================================================================================
 */
#include <climits>

#include <iostream>

#include "FASTASequence.hpp"
#include "FASTAReader.hpp"
#include "CommandLineParser.hpp"
#include "ChangeListID.hpp"
#include "utils/TimeUtils.hpp"
#include "utils/RangeUtils.hpp"
#include "utils/SMRTReadUtils.hpp"
#include "algorithms/alignment/DistanceMatrixScoreFunction.hpp"
#include "algorithms/alignment/AlignmentUtils.hpp"
#include "algorithms/alignment/StringToScoreMatrix.hpp"
#include "sam/SAMReader.hpp"
#include "format/SAMPrinter.hpp"
#include "datastructures/alignment/AlignmentCandidate.hpp"
#include "datastructures/alignment/FilterCriteria.hpp"
#include "metagenome/TitleTable.hpp"
#include "datastructures/alignment/SAMToAlignmentCandidateAdapter.hpp"
#include "GFFFile.hpp"
#include "defs.h"
#include "RegisterFilterOptions.h"

//#define USE_GOOGLE_PROFILER
#ifdef USE_GOOGLE_PROFILER
#include "gperftools/profiler.h"
#endif

char VERSION[] = "v0.1.0";
char PERFORCE_VERSION_STRING[] = "$Change: 134995 $";
// By default negative score is better.
ScoreSign scoreSign = ScoreSign::NEGATIVE;

// Compare SAMAlignment objects by qName, score and 
// target positions.
bool byQNameScoreTStart(const SAMAlignment & a, 
                        const SAMAlignment & b) {
    if (a.qName == b.qName) {
        if (a.score == b.score) 
            return a.pos < b.pos;
        return Score(a.score, scoreSign).WorseThan(Score(b.score, scoreSign));
    }
    return (a.qName < b.qName);
}

// Compare SAMAlignment objects by rName and qName 
bool byRNameQName(const SAMAlignment & a, 
                  const SAMAlignment & b) {
    if (a.rName == b.rName) {
        return a.qName < b.qName;
    }
    return (a.rName < b.rName);
}

// Get the next group of SAM alignments that have the same qName from
// allSAMAlignments[groupBegin ... groupEnd)
// Note that allSAMAlignments is already sorted by qName, score and tPos.
void GetNextSAMAlignmentGroup(vector<SAMAlignment> & allSAMAlignments, 
                              unsigned int groupBegin, 
                              unsigned int & groupEnd) {
    assert(groupBegin < allSAMAlignments.size());
    groupEnd = groupBegin + 1;
    string queryName = allSAMAlignments[groupBegin].qName;
    while(groupEnd < allSAMAlignments.size()) {
        if (allSAMAlignments[groupEnd].qName == queryName) 
            groupEnd ++;
        else break;
    }
}

// Get the best SAM alignments whose alignment score are the best. 
// Assume that alignments in allSAMAlignments[groupBegin, groupEnd)
// all have the same queryName and are sorted by score and tPos 
// asscendingly: worst, ...., best
void GetBestSAMAlignmentsInGroup(vector<SAMAlignment> & allSAMAlignments, 
                                 const unsigned int & groupBegin, 
                                 const unsigned int & groupEnd, 
                                 unsigned int & bestBegin,
                                 unsigned int & bestEnd) {
    assert(groupEnd  <= allSAMAlignments.size() and
           groupBegin < groupEnd);

    bestEnd = groupEnd;
    bestBegin = groupEnd - 1;
    int groupBestScore = allSAMAlignments[bestBegin].score;
    string queryName = allSAMAlignments[bestBegin].qName;
    while (bestBegin >= groupBegin and bestBegin < groupEnd) {
        assert(allSAMAlignments[bestBegin].qName == queryName);
        if (allSAMAlignments[bestBegin].score == groupBestScore)
            bestBegin -= 1;
        else break;
    }
    bestBegin += 1;
}


// Apply hit policy to a group of SAM alignments and return indices
// of the selected alignments.
vector<unsigned int> ApplyHitPolicy(HitPolicy & hitPolicy, 
                                    vector<SAMAlignment> & allSAMAlignments, 
                                    const unsigned int & groupBegin, 
                                    const unsigned int & groupEnd) {
    vector<unsigned int> hitIndices;
    if (hitPolicy.IsAll()) {
        for(unsigned int i = groupBegin; i < groupEnd; i++){
            hitIndices.push_back(i);
        }
    } else if (hitPolicy.IsRandom()) {
        hitIndices.push_back(rand()%(groupEnd - groupBegin) + groupBegin);
    } else {
        unsigned int bestBegin, bestEnd;
        GetBestSAMAlignmentsInGroup(allSAMAlignments, groupBegin, groupEnd,
                                    bestBegin, bestEnd);
        if (hitPolicy.IsAllbest()) {
            for(unsigned int i = bestBegin; i < bestEnd; i++){
                hitIndices.push_back(i);
            }
        } else if (hitPolicy.IsRandombest()) {
            hitIndices.push_back(rand()%(bestEnd-bestBegin) + bestBegin);
        } else if (hitPolicy.IsLeftmost()) {
            hitIndices.push_back(bestBegin);
        } else {
            assert(false);
        }
    }
    return hitIndices;
}

// Convert references[...].title in reference.fasta to their corresponding
// indices in the title table.
void ConvertTitlesToTitleTableIndices(vector<FASTASequence> & references,
        string & titleTableName) {
    TitleTable tt;
    tt.Read(titleTableName);
    for(int i = 0; i < references.size(); i++) {
        string title = references[i].GetTitle();
        int idx = -1;
        if (tt.Lookup(title, idx)) {
            stringstream ss;
            ss << idx; 
            references[i].CopyTitle(ss.str());
        } else {
            cout << "ERROR, reference " << title << " does not exist "
                 << " in the title table " << titleTableName << ". The "
                 << "reference fasta and the title table do not match."
                 << endl;
            exit(1);
        }
    }
    tt.Free();
}

// Return true if the alignment can only map to an adapter specified 
// in the adapter GFF file. 
// A sample record in adapter GFF file:
// ref000001   .   adapter 10955   10999   0.00    +   .   xxxx
// ref000001   .   adapter 32886   32930   0.00    +   .   xxxx 
// Note that the first field (e.g., 'ref000001') is id of sequence  
// in a reference repository, not sequence name, so we need to 
// reconstruct the mapping between sequence id and sequence name.
bool CheckAdapterOnly(GFFFile & adapterGffFile, //Adapter gff file
    AlignmentCandidate<> & alignment, // An alignment
    map<string, int> & refNameToIndex) { 
    // Map target sequence name to its index in reference repository.
    if (refNameToIndex.find(alignment.tName) == refNameToIndex.end()) {
        // This should not happen ...
        cout << "ERROR, could not find alignment target name "
             << alignment.tName << " in the reference file." << endl;
        exit(1);
    }
    int refNameIndex = refNameToIndex[alignment.tName];
    char buf [16];
    sprintf(buf, "ref%06d", refNameIndex + 1);
    // Reconstruct ref id in the format "ref00000?".
    string refNameId(buf);
    int FUZZY_OVERLAP = 20;
    for(int eindex = 0; eindex < adapterGffFile.entries.size();
            eindex++) { 
        GFFEntry & entry = adapterGffFile.entries[eindex];
        // Convert each GFF record from 1-based inclusive to 
        // 0-based exclusive.
        if (entry.type == "adapter" and 
            (entry.name == alignment.tName or
             entry.name == refNameId)) {
            UInt estart = entry.start - 1;
            UInt eend = entry.end;
            if (entry.strand == '-') {
                UInt tmp = estart;
                estart = alignment.tLength - 1 - eend;
                eend = alignment.tLength - 1 - tmp;
            }
            if (not (eend < alignment.GenomicTBegin() or
                 estart > alignment.GenomicTEnd())) {
                int lengthUnion = max(eend, alignment.GenomicTEnd()) -
                                  min(estart, alignment.GenomicTBegin());
                if (lengthUnion < eend - estart + FUZZY_OVERLAP) {
                    return true;
                }
            }
        }
    }
    return false;
}

int main(int argc, char* argv[]) {
#ifdef USE_GOOGLE_PROFILER
    char *profileFileName = getenv("CPUPROFILE");
    if (profileFileName != NULL) {
      ProfilerStart(profileFileName);
    }
    else {
      ProfilerStart("google_profile.txt");
    }
#endif

    // Register inputs and outputs.
    string samFileName, refFileName, outFileName;

    CommandLineParser clp;
    clp.RegisterStringOption("file.sam", &samFileName,
                             "Input SAM file.");
    clp.RegisterStringOption("reference.fasta", &refFileName,
                             "Reference used to generate reads.");
    clp.RegisterStringOption("out.sam", &outFileName,
                             "Output SAM file.");
    clp.RegisterPreviousFlagsAsHidden();

    // Register filter criteria options.
    int minAlnLength = 50;
    float minPctSimilarity = 70, minPctAccuracy = 70;
    string hitPolicyStr = "randombest";
    bool useScoreCutoff = false;
    int  scoreCutoff = INF_INT;
    int  scoreSignInt = -1;
    RegisterFilterOptions(clp, minAlnLength, minPctSimilarity, 
                          minPctAccuracy, hitPolicyStr, useScoreCutoff,
                          scoreSignInt, scoreCutoff);

    int seed = 1; 
    clp.RegisterIntOption("seed", &seed,
            "(1)  Seed for random number generator.\n"
            "If seed is 0, then use current time as seed.",
            CommandLineParser::Integer);

    string holeNumberStr;
    Ranges holeNumberRanges;
    clp.RegisterStringOption("holeNumbers", &holeNumberStr,
            "A string of comma-delimited hole number ranges to output hits, "
            "such as '1,2,10-12'. "
            "This requires hit titles to be in SMRT read title format.");

    bool parseSmrtTitle = false;
    clp.RegisterFlagOption("smrtTitle", &parseSmrtTitle,
            "Use this option when filtering alignments generated by "
            "programs other than blasr, e.g. bwa-sw or gmap. "
            "  Parse read coordinates from the SMRT read title. " 
            "The title is in the format /name/hole/coordinates, where"
            " coordinates are in the format \\d+_\\d+, and represent "
            "the interval of the read that was aligned.");
    /* This experimental option can be useful for metagenomics, in which case
     * there are hundreds of sequences in the target, of which many titles are
     * long and may contain white spaces (e.g., ' ', '\t'). 
     * In order to save disc space and avoid the (possibly) none unique mapping
     * between full and short reference names, one may call blasr with 
     * -titleTable option to represent all target sequences in the output
     * by their indices in the title table.*/

    string titleTableName = "";
    clp.RegisterStringOption("titleTable", &titleTableName,
            "Use this experimental option when filtering alignments generated by "
            "blasr with -titleTable titleTableName, in which case "
            "reference titles in SAM are represented by their "
            "indices (e.g., 0, 1, 2, ...) in the title table.");

    string adapterGffFileName = "";
    clp.RegisterStringOption("filterAdapterOnly", &adapterGffFileName,
            "Use this option to remove reads which can only map to adapters " 
            "specified in the GFF file.");

    bool verbose = false;
    clp.RegisterFlagOption("v", &verbose, "Be verbose.");

    clp.SetExamples(
            "Because SAM has optional tags that have different meanings"
            " in different programs, careful usage is required in order "
            "to have proper output.  The \"xs\" tag in bwa-sw is used to "
            "show the suboptimal score, but in PacBio SAM (blasr) it is "
            "defined as the start in the query sequence of the alignment.\n"
            "When \"-smrtTitle\" is specified, the xs tag is ignored, but "
            "when it is not specified, the coordinates given by the xs and "
            "xe tags are used to define the interval of a read that is "
            "aligned.  The CIGAR string is relative to this interval.");

    clp.ParseCommandLine(argc, argv);

    // Set random number seed. 
    if (seed == 0) {
        srand(time(NULL));
    } else {
        srand(seed);
    }
    
    scoreSign = (scoreSignInt == -1)?ScoreSign::NEGATIVE:ScoreSign::POSITIVE;
    Score s(static_cast<float>(scoreCutoff), scoreSign);
    FilterCriteria filterCriteria(minAlnLength, minPctSimilarity, 
                                  minPctAccuracy, true, s);
    filterCriteria.Verbose(verbose);
    HitPolicy hitPolicy(hitPolicyStr, scoreSign);
                                  
    string errMsg;
    if (not filterCriteria.MakeSane(errMsg)) {
        cout << errMsg << endl;
        exit(1);
    }

    // Parse hole number ranges. 
    if (holeNumberStr.size() != 0) {
        if (not holeNumberRanges.setRanges(holeNumberStr)) {
            cout << "Could not parse hole number ranges: "
                 << holeNumberStr << "." << endl;
            exit(1);
        } 
    }

    // Open output file.
    ostream * outFilePtr = &cout;
	ofstream outFileStrm;
	if (outFileName != "") {
		CrucialOpen(outFileName, outFileStrm, std::ios::out);
		outFilePtr = &outFileStrm;
	}
    
    GFFFile adapterGffFile;
    if (adapterGffFileName != "")
        adapterGffFile.ReadAll(adapterGffFileName);
    
    SAMReader<SAMFullReferenceSequence, SAMReadGroup, SAMAlignment> samReader;
    FASTAReader fastaReader;

    //
    // Initialize samReader and fastaReader.
    //
    samReader.Initialize(samFileName);
    fastaReader.Initialize(refFileName);

    //
    // Configure the file log.
    //
    string command;
    CommandLineParser::CommandLineToString(argc, argv, command);
    string log = "Filter sam hits.";
    string program = "samFilter";
    string versionString = VERSION;
    AppendPerforceChangelist(PERFORCE_VERSION_STRING, versionString);

    //
    // Read necessary input.
    //
    vector<FASTASequence> references;
    fastaReader.ReadAllSequences(references);

    // If the SAM file is generated by blasr with -titleTable,
    // then references in the SAM are represented by 
    // their corresponding indices in the title table.
    // In that case, we need to convert reference titles in fasta file
    // to their corresponding indices in the title table, such that
    // references in both SAM and fasta files are represented
    // by title table indices and therefore can match.
    if (titleTableName != "") {
        ConvertTitlesToTitleTableIndices(references, titleTableName);
    }
 
    AlignmentSet<SAMFullReferenceSequence, SAMReadGroup, SAMAlignment> alignmentSet;
    vector<string> allHeaders = samReader.ReadHeader(alignmentSet); 

    // Process SAM Header.
    string commandLineString;
    clp.CommandLineToString(argc, argv, commandLineString);
    allHeaders.push_back("@PG\tID:SAMFILTER\tVN:" + versionString + \
                         "\tCL:" + program + " " + commandLineString);
    for (int i = 0; i < allHeaders.size(); i++) {
        outFileStrm << allHeaders[i] << endl;
    }

    //
    // The order of references in vector<FASTASequence> references and
    // AlignmentSet<, , >alignmentSet.references can be different.
    // Rearrange alignmentSet.references such that they are ordered in
    // exactly the same way as vector<FASTASequence> references.
    //
    alignmentSet.RearrangeReferences(references);

    // Map reference name obtained from SAM file to indices
    map<string, int> refNameToIndex;
    for (int i = 0; i < references.size(); i++) {
        string refName = alignmentSet.references[i].GetSequenceName();
        refNameToIndex[refName] = i;
    }

    //
    // Store the alignments.
    //
    SAMAlignment samAlignment;
    int alignIndex = 0; 

    //
    // For 150K, each chip produces about 300M sequences 
    // (not including quality values and etc.).
    // Let's assume that the sam file and reference data can 
    // fit in the memory. 
    // Need to scale for larger sequal data in the future.
    //
    vector<SAMAlignment> allSAMAlignments;
    while (samReader.GetNextAlignment(samAlignment)) {
        if (samAlignment.rName == "*") {
            continue;
        }

        if (parseSmrtTitle and holeNumberStr.size() != 0) {
            string movieName;
            int thisHoleNumber;
            if (not ParsePBIReadName(samAlignment.qName, 
                                     movieName, 
                                     thisHoleNumber)) {
                cout << "ERROR, could not parse SMRT title: "
                     << samAlignment.qName << "." << endl;
                exit(1);
            }
            if (not holeNumberRanges.contains(UInt(thisHoleNumber))) {
                if (verbose) 
                    cout << thisHoleNumber << " is not in range." << endl; 
                continue;
            }
        }

        if (samAlignment.cigar.find('P') != string::npos) {
            cout << "WARNING. Could not process SAM record with 'P' in "
                 << "its cigar string." << endl;
            continue;
        }

        vector<AlignmentCandidate<> > convertedAlignments;
        SAMAlignmentsToCandidates(samAlignment, 
                references, refNameToIndex,
                convertedAlignments, parseSmrtTitle, false);
        
        if (convertedAlignments.size() > 1) {
            cout << "WARNING. Ignore multiple segments." << endl;
            continue;
        }

        for (int i = 0; i < 1; i++) {
            AlignmentCandidate<> & alignment = convertedAlignments[i];

            //score func does not matter
            DistanceMatrixScoreFunction<DNASequence, DNASequence> distFunc; 
            ComputeAlignmentStats(alignment, alignment.qAlignedSeq.seq, 
                                  alignment.tAlignedSeq.seq, distFunc);
                                  
            // Check whether this alignment can only map to adapters in 
            // the adapter GFF file.
            if (adapterGffFileName != "" and 
                CheckAdapterOnly(adapterGffFile, alignment, refNameToIndex)) {
                if (verbose)
                    cout << alignment.qName << " filter adapter only."
                         << endl;
                continue;
            }

            // Assign score to samAlignment.
            samAlignment.score = samAlignment.as;

            if (not filterCriteria.Satisfy(static_cast<AlignmentCandidate<> *>(&alignment))) {
                continue;
            }
            allSAMAlignments.push_back( samAlignment ); 

            alignment.FreeSubsequences();
        }
        ++alignIndex;
    }

    // Sort all SAM alignments by qName, score and target position.
    sort(allSAMAlignments.begin(), allSAMAlignments.end(), 
         byQNameScoreTStart);

    unsigned int groupBegin = 0;
    unsigned int groupEnd = -1;
    vector<SAMAlignment> filteredSAMAlignments;
    while(groupBegin < allSAMAlignments.size()) {
        // Get the next group of SAM alignments which have the same qName
        // from allSAMAlignments[groupBegin ... groupEnd)
        GetNextSAMAlignmentGroup(allSAMAlignments, groupBegin, groupEnd);
        vector<unsigned int> hitIndices = ApplyHitPolicy(
                hitPolicy, allSAMAlignments, groupBegin, groupEnd);
        for(unsigned int i = 0; i < hitIndices.size(); i++) {
            filteredSAMAlignments.push_back(allSAMAlignments[hitIndices[i]]);
        }
        groupBegin = groupEnd;
    }

    // Sort all SAM alignments by reference name and query name
    sort(filteredSAMAlignments.begin(), filteredSAMAlignments.end(), 
         byRNameQName);

    for(unsigned int i = 0; i < filteredSAMAlignments.size(); i++) {
        filteredSAMAlignments[i].PrintSAMAlignment(outFileStrm);
    }

	if (outFileName != "") {
		outFileStrm.close();
	}
#ifdef USE_GOOGLE_PROFILER
  ProfilerStop();
#endif
    return 0;
}
