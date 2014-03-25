#include <string>
#include <iostream>
#include <vector>
#include <set>
#include <sstream>
#include <pthread.h>
#include <stdlib.h>
#include <time.h>
#include <signal.h>
#include <execinfo.h>

#include "MappingIPC.h"
#include "MappingSemaphores.h"

#include "CCSSequence.hpp"
#include "SMRTSequence.hpp"
#include "FASTASequence.hpp"
#include "FASTAReader.hpp"
#include "SeqUtils.hpp"
#include "defs.h"
#include "utils.h"

#include "tuples/DNATuple.hpp"
#include "tuples/HashedTupleList.h" // partial port
#include "algorithms/compare/CompareStrings.hpp"
#include "algorithms/alignment/AffineKBandAlign.hpp"
#include "algorithms/alignment/GuidedAlign.hpp"
#include "algorithms/alignment/AffineGuidedAlign.hpp"
#include "algorithms/alignment/FullQVAlign.hpp"
#include "algorithms/alignment/ExtendAlign.hpp"
#include "algorithms/alignment/OneGapAlignment.hpp"
#include "algorithms/alignment/AlignmentUtils.hpp"
#include "algorithms/alignment/QualityValueScoreFunction.h" // not ported
#include "algorithms/alignment/IDSScoreFunction.hpp"
#include "algorithms/alignment/DistanceMatrixScoreFunction.hpp"
#include "algorithms/alignment/StringToScoreMatrix.h" // not ported
#include "algorithms/alignment/AlignmentFormats.h" // not ported
#include "algorithms/anchoring/LISPValue.h" // not ported
#include "algorithms/anchoring/LISPValueWeightor.h" // not ported
#include "algorithms/anchoring/LISSizeWeightor.h" // not ported
#include "algorithms/anchoring/LISQValueWeightor.h" // not ported
#include "algorithms/anchoring/FindMaxInterval.h" // not ported
#include "algorithms/anchoring/MapBySuffixArray.h" // not ported
#include "datastructures/anchoring/ClusterList.h" // not ported
#include "algorithms/anchoring/ClusterProbability.h" // not ported
#include "algorithms/anchoring/BWTSearch.h" // not ported
#include "datastructures/metagenome/SequenceIndexDatabase.h" // not ported
#include "datastructures/metagenome/TitleTable.h" // not ported
#include "suffixarray/SharedSuffixArray.hpp"
#include "suffixarray/SuffixArrayTypes.hpp"
#include "tuples/TupleCountTable.hpp" 
#include "datastructures/anchoring/WeightedInterval.h" // not ported
#include "anchoring/AnchorParameters.hpp"
#include "datastructures/alignment/AlignmentCandidate.hpp"
#include "datastructures/alignment/AlignmentBlock.h" // partial port
#include "datastructures/alignment/AlignmentContext.hpp"
#include "datastructures/mapping/MappingMetrics.h" // not ported
#include "datastructures/reads/ReadInterval.h" // not ported
#include "utils/FileOfFileNames.h" // not ported
#include "utils/RegionUtils.h" // not ported
#include "utils/TimeUtils.h" // not ported
#include "qvs/QualityTransform.hpp"
#include "files/ReaderAgglomerate.hpp"
#include "files/CCSIterator.h" // not ported
#include "files/FragmentCCSIterator.h" // not ported
#include "HDFRegionTableReader.hpp"
#include "datastructures/bwt/BWT.h" // not ported
#include "datastructures/sequence/PackedDNASequence.h" // not ported
#include "CommandLineParser.h" // not ported
#include "qvs/QualityValue.hpp"
#include "statistics/VarianceAccumulator.h" // not ported
#include "statistics/pdfs.h" // not ported
#include "statistics/cdfs.h" // not ported
#include "statistics/statutils.h" // not ported
#include "statistics/LookupAnchorDistribution.h" // not ported

#define MAX_PHRED_SCORE 254
#define MAPQV_END_ALIGN_WIGGLE 5

// #define USE_GOOGLE_PROFILER

#ifdef USE_GOOGLE_PROFILER
#include "gperftools/profiler.h"
#endif

using namespace std;

MappingSemaphores semaphores;

/*
 * Declare global structures that are shared between threads.
 */

ostream *outFilePtr = NULL;
//ofstream clOut, nsOut, mcOut;

HDFRegionTableReader *regionTableReader = NULL;

typedef SMRTSequence T_Sequence;
typedef FASTASequence T_GenomeSequence;
typedef DNASuffixArray T_SuffixArray;
typedef DNATuple T_Tuple;

typedef LISPValueWeightor<T_GenomeSequence, DNATuple, vector<ChainedMatchPos> >  PValueWeightor;
typedef LISSMatchFrequencyPValueWeightor<T_GenomeSequence, DNATuple, vector<ChainedMatchPos> >  MultiplicityPValueWeightor;

ReaderAgglomerate *reader = NULL;

typedef MappingData<T_SuffixArray, T_GenomeSequence, T_Tuple> MappingIPC;

class ClusterInformation {
public:
  int maxClusterSize;
  float meanAnchorBasesPerRead;
  float sdAnchorBasesPerRead;
  int score;
  float pctSimilarity;
  int readLength;
  float nStdDev ;
  int numSignificant;
  int numForward, numReverse;
};


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
  int GetNAlignedSeq() {
    return subreadAlignments.size();
  }

  bool AllSubreadsHaveAlignments() {
    int i, nAlignedSeq;
    nAlignedSeq = subreadAlignments.size();
    for (i = 0; i < nAlignedSeq; i++) {
      if (subreadAlignments[i].size() == 0) {
        return false;
      }
    }
    return true;
  }

  void Clear() {
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
      subreads[i].FreeIfControlled();
      if (subreads[i].title != NULL) {
        delete[] subreads[i].title;
        subreads[i].title = NULL;
      }
    }
    subreadAlignments.clear();
    read.Free();
  }

  void Resize(int nSeq) {
    subreadAlignments.resize(nSeq);
    subreads.resize(nSeq);
  }

  void CheckSeqIndex(int seqIndex) {
    if ( seqIndex < 0 or seqIndex >= subreads.size() ) {
      cout << "ERROR, adding a sequence to an unallocated position." 
           << endl;
      assert(0);
    }
  }

  void SetSequence(int seqIndex, SMRTSequence &seq) {
    CheckSeqIndex(seqIndex);
    subreads[seqIndex] = seq;
  }

  void AddAlignmentForSeq(int seqIndex, T_AlignmentCandidate *alignmentPtr) {
    CheckSeqIndex(seqIndex);
    subreadAlignments[seqIndex].push_back(alignmentPtr);
  }

  void AddAlignmentsForSeq(int seqIndex, vector<T_AlignmentCandidate*> &seqAlignmentPtrs) {
    CheckSeqIndex(seqIndex);
    subreadAlignments[seqIndex].insert(subreadAlignments[seqIndex].end(), seqAlignmentPtrs.begin(), seqAlignmentPtrs.end());
  }

  // Copy all T_AlignmentCandidate objects (to which subreadAlignment[seqIndex]
  // is pointing) to newly created objects, and then return pointers to the new
  // objects.
  vector<T_AlignmentCandidate*> CopySubreadAlignments(int seqIndex) {
    vector<T_AlignmentCandidate*> ret;
    for (int i=0; i<subreadAlignments[seqIndex].size(); i++) {
      T_AlignmentCandidate * q = new T_AlignmentCandidate();
      *q = *(subreadAlignments[seqIndex][i]);
      ret.push_back(q);
    }
    return ret;
  }

  void Print(ostream &out=cout) { 
    out << "A ReadAlignments object with " 
        << subreadAlignments.size()
        << " groups of subread alignments." << endl;
        for (int i = 0; i < subreadAlignments.size(); i++) {
            out << "  subreadAlignment group [" << i << "/" 
                << subreadAlignments.size() << "] has "
                << subreadAlignments[i].size() << " alignments." << endl;
            for(int j = 0; j < subreadAlignments[i].size(); j++) {
                out << "    [" << i << "][" << j << "/" 
                    << subreadAlignments[i].size() << "]" << endl;
                subreadAlignments[i][j]->Print(out);
            } 
        }
        /* subreads may have been freed or not initialized when 
             * being printed. 
        for (int i = 0; i < subreads.size(); i++) {
            out << "  subread [" << i << "/" << subreads.size()
                << "]: ";
            subreads[i].Print(out);
        } */
        out << "  read: ";
        read.Print(out);
        out << endl << endl;
    }

  ~ReadAlignments() {
    read.FreeIfControlled();
  }
};


class SortAlignmentPointersByScore {
public:
  int operator()(T_AlignmentCandidate *lhs, T_AlignmentCandidate* rhs) {
    if (lhs->score == rhs->score) {
      return lhs->tPos + lhs->tAlignedSeqPos < rhs->tPos + rhs->tAlignedSeqPos;
    }
    else {
      return lhs->score < rhs->score;
    }
  }
};

class SortAlignmentPointersByMapQV {
public:
  int operator()(T_AlignmentCandidate *lhs, T_AlignmentCandidate* rhs) {
    if (lhs->mapQV == rhs->mapQV) {
      if (lhs->score == rhs->score) {
        return lhs->tPos + lhs->tAlignedSeqPos < rhs->tPos + rhs->tAlignedSeqPos;
      }
      else {
        return lhs->score < rhs->score;
      }
    }
    else {
      return lhs->mapQV > rhs->mapQV;
    }
  }
};

string GetMajorVersion() {
  return "1.3.1";
}

void GetVersion(string &version) {
  string perforceVersionString("$Change$");
  version = GetMajorVersion();
  if (perforceVersionString.size() > 12) {
    version.insert(version.size(), ".");
    version.insert(version.size(), perforceVersionString, 9, perforceVersionString.size() - 11);
  }
}


void MakeSAMHDString(string &hdString) {
  stringstream out;
  out << "@HD\t" << "VN:" << GetMajorVersion();
  hdString = out.str();
}

void MakeSAMPGString(string &commandLineString, string &pgString) {
  stringstream out;
  string version;
  GetVersion(version);
  out << "@PG\t" << "ID:BLASR\tVN:"<< version<<"\tCL:"<<commandLineString;
  pgString = out.str();
}


void ParseChipIdFromMovieName(string &movieName, string &chipId) {
  vector<string> movieNameFields;
  Tokenize(movieName, "_", movieNameFields);
  if (movieNameFields.size() == 1) {
    chipId = movieNameFields[0];
  }
  else if (movieNameFields.size() > 4) {
    chipId = movieNameFields[3];
  }
  else {
    chipId = "NO_CHIP_ID";
  }
}

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

  void Reset() {
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
};

void SetHelp(string &str) {
  stringstream helpStream;
  helpStream << "   Options for blasr " << endl
             << "   Basic usage: 'blasr reads.{fasta,bax.h5} genome.fasta [-options] " << endl
             << " option\tDescription (default_value)." << endl << endl
             << " Input Files." << endl
             << "   reads.fasta is a multi-fasta file of reads.  While any fasta file is valid input, " <<endl
             << "               it is preferable to use plx.h5 or bax.h5 files because they contain" << endl
             << "               more rich quality value information." << endl << endl
             << "   reads.bax.h5|reads.plx.h5 Is the native output format in Hierarchical Data Format of " <<endl
             << "               SMRT reads. This is the preferred input to blasr because rich quality" << endl
             << "               value (insertion,deletion, and substitution quality values) information is " << endl
             << "               maintained.  The extra quality information improves variant detection and mapping"<<endl
             << "               speed." << endl <<endl
             << "   -sa suffixArrayFile"<< endl
             << "               Use the suffix array 'sa' for detecting matches" << endl 
             << "               between the reads and the reference.  The suffix" << endl 
             << "               array has been prepared by the sawriter program." << endl << endl 
             << "   -ctab tab "<<endl
             << "               A table of tuple counts used to estimate match significance.  This is " << endl
             << "               by the program 'printTupleCountTable'.  While it is quick to generate on " << endl
             << "               the fly, if there are many invocations of blasr, it is useful to"<<endl
             << "               precompute the ctab." <<endl << endl
             << "   -regionTable table" << endl
             << "               Read in a read-region table in HDF format for masking portions of reads." << endl
             << "               This may be a single table if there is just one input file, " << endl
             << "               or a fofn.  When a region table is specified, any region table inside " << endl
             << "               the reads.plx.h5 or reads.bax.h5 files are ignored."<< endl
//             << "   -ccsFofn ccsFofn" << endl
//             << "               Read in a ciruclar consensus sequence (ccs) file in HDF format. " << endl
//             << "               This may be a single ccs.h5 file if there is just one input file, " << endl
//             << "               or a fofn. When a ccs file (or fofn) is specified, any ccs data inside the " << endl
//             << "               read.bax.h5 file (or files) are ignored." << endl
             << endl 
             << " Options for modifying reads. There is ancilliary information about substrings of reads " << endl
             << "               that is stored in a 'region table' for each read file.  Because " << endl
             << "               HDF is used, the region table may be part of the .bax.h5 or .plx.h5 file," << endl
             << "               or a separate file.  A contiguously read substring from the template " << endl
             << "               is a subread, and any read may contain multiple subreads. The boundaries " << endl
             << "               of the subreads may be inferred from the region table either directly or " <<endl
             << "               by definition of adapter boundaries.  Typically region tables also" << endl
             << "               contain information for the location of the high and low quality regions of"<<endl
             << "               reads.  Reads produced by spurrious reads from empty ZMWs have a high"<<endl
             << "               quality start coordinate equal to high quality end, making no usable read." <<endl
             << "   -useccs   " << endl
             << "               Align the circular consensus sequence (ccs), then report alignments" << endl
             << "               of the ccs subreads to the window that the ccs was mapped to.  Only " << endl
             << "               alignments of the subreads are reported." << endl  
             << "   -useccsall"<<endl 
             << "               Similar to -useccs, except all subreads are aligned, rather than just" << endl
             << "               the subreads used to call the ccs.  This will include reads that only"<<endl
             << "               cover part of the template." << endl
             << "   -useccsdenovo" << endl
             << "               Align the circular consensus, and report only the alignment of the ccs"<<endl
             << "               sequence." << endl
             << "   -noSplitSubreads (false)" <<endl
             << "               Do not split subreads at adapters.  This is typically only " << endl
             << "               useful when the genome in an unrolled version of a known template, and " << endl
             << "               contains template-adapter-reverse_template sequence." << endl
             << "   -ignoreRegions(false)" << endl
             << "               Ignore any information in the region table." << endl
             << "   -ignoreHQRegions (false)Ignore any hq regions in the region table." << endl
             << endl
             << " Alignment Output." << endl
             << "   -bestn n (10)" <<endl
             << "               Report the top 'n' alignments." << endl
             << "   -sam        Write output in SAM format." << endl
             << "   -clipping [none|hard|soft] (none)" << endl
             << "               Use no/hard/soft clipping for SAM output."<< endl
             << "   -out out (terminal)  " << endl
             << "               Write output to 'out'." << endl
             << "   -unaligned file" << endl
             << "               Output reads that are not aligned to 'file'" << endl
             << "   -m t           " << endl
             << "               If not printing SAM, modify the output of the alignment." << endl
             << "                t=" << StickPrint <<   " Print blast like output with |'s connecting matched nucleotides." << endl 
             << "                  " << SummaryPrint << " Print only a summary: score and pos." << endl 
             << "                  " << CompareXML <<   " Print in Compare.xml format." << endl 
             << "                  " << Vulgar <<       " Print in vulgar format (deprecated)." << endl
             << "                  " << Interval <<     " Print a longer tabular version of the alignment." << endl 
             << "                  " << CompareSequencesParsable  << " Print in a machine-parsable format that is read by compareSequences.py." << endl
             << "   -noSortRefinedAlignments (false) " << endl
             << "               Once candidate alignments are generated and scored via sparse dynamic "<< endl
             << "               programming, they are rescored using local alignment that accounts " << endl
             << "               for different error profiles." <<endl
             << "               Resorting based on the local alignment may change the order the hits are returned." << endl
             << "   -allowAdjacentIndels " << endl
             << "               When specified, adjacent insertion or deletions are allowed. Otherwise, adjacent " << endl
             << "               insertion and deletions are merged into one operation.  Using quality values " << endl
             << "               to guide pairwise alignments may dictate that the higher probability alignment "<<endl
             << "               contains adjacent insertions or deletions.  Current tools such as GATK do not permit" << endl
             << "               this and so they are not reported by default." << endl
             << "   -header" <<endl
             << "               Print a header as the first line of the output file describing the contents of each column."<<endl
             << "   -titleTable tab (NULL) " << endl
             << "               Construct a table of reference sequence titles.  The reference sequences are " << endl
             << "               enumerated by row, 0,1,...  The reference index is printed in alignment results" << endl
             << "               rather than the full reference name.  This makes output concise, particularly when" <<endl
             << "               very verbose titles exist in reference names."<<endl
             << "   -minPctIdentity p (0)"<< endl
             << "               Only report alignments if they are greater than p percent identity."<<endl
             << "   -unaligned file" << endl
             << "               Output reads that are not aligned to 'file'." << endl
             << "   -holeNumbers LIST " << endl
             << "               When specified, only align reads whose ZMW hole numbers are in LIST." << endl
             << "               LIST is a comma-delimited string of ranges, such as '1,2,3,10-13'." << endl
             << "               This option only works when reads are in base or pulse h5 format." << endl
             << "   -placeRepeatsRandomly (false)" << endl
             << "               When there are multiple positions to map a read with equal alignment scores, place the" << endl
             << "               read randomly at one of them." <<endl
             << endl 
             << " Options for anchoring alignment regions. This will have the greatest effect on speed and sensitivity." << endl
             << "   -minMatch m (12) " << endl
             << "               Minimum seed length.  Higher minMatch will speed up alignment, " << endl
             << "               but decrease sensitivity." << endl
//             << "   -maxExpand M (1)" << endl
//             << "               Perform no more than M iterations of searches through the suffix " << endl
//             << "               array for matches. At each iteration, all matches of length LCPi-M" << endl
//             << "               are found, where LCPi is the length of the longest common prefix " << endl
//             << "               between the string at i and anywhere in the genome."<<endl
//             << "               The number of matches grows as M increases, and can become very large with M > 3." << endl
             << "   -maxMatch l (inf)" << endl
             << "               Stop mapping a read to the genome when the lcp length reaches l.  " << endl
             << "               This is useful when the query is part of the reference, for example when " <<endl
             << "               constructing pairwise alignments for de novo assembly."<<endl
             << "   -maxLCPLength l (inf)" << endl
             << "               The same as -maxMatch." << endl
             << "   -maxAnchorsPerPosition m (inf) " << endl
             << "               Do not add anchors from a position if it matches to more than 'm' locations in the target." << endl
//             << "   -advanceHalf (false) " << endl
//             << "               A trick for speeding up alignments at the cost of sensitivity.  If " << endl
//             << "               a cluster of anchors of size n, (a1,...,an) is found, normally anchors " << endl
//             << "               (a2,...an) of size n-1 is also clustered to make sure a1 did not decrease the " << endl
//             << "               cluster score.  When advanceHalf is specified, clustering begins at a_(n/2)."<<endl<< endl
             << "   -advanceExactMatches E (0)" << endl
             << "               Another trick for speeding up alignments with match - E fewer anchors.  Rather than" << endl 
             << "               finding anchors between the read and the genome at every position in the read, " <<endl
             << "               when an anchor is found at position i in a read of length L, the next position " << endl
             << "               in a read to find an anchor is at i+L-E." << endl
             << "               Use this when alignining already assembled contigs." << endl
             << "   -nCandidates n (10)" << endl 
             << "               Keep up to 'n' candidates for the best alignment.  A large value of n will slow mapping" << endl
             << "               because the slower dynamic programming steps are applied to more clusters of anchors" <<endl
             << "               which can be a rate limiting step when reads are very long."<<endl
             << "   -concordant(false)" << endl
             << "               Map all subreads of a zmw (hole) to where the longest full pass subread of the zmw " << endl
             << "               aligned to. This requires to use the region table and hq regions." << endl
             << "               This option only works when reads are in base or pulse h5 format." << endl
             << endl
             << "  Options for Refining Hits." << endl
//             << "   -indelRate i (0.30)" << endl
//             << "               The approximate maximum rate to allow drifting from the diagonal." <<endl << endl
             << "   -sdpTupleSize K (11)" << endl
             << "               Use matches of length K to speed dynamic programming alignments.  This controls" <<endl
             << "               accuracy of assigning gaps in pairwise alignments once a mapping has been found,"<<endl
             << "               rather than mapping sensitivity itself."<<endl
             << "   -scoreMatrix \"score matrix string\" " << endl
             << "               Specify an alternative score matrix for scoring fasta reads.  The matrix is " << endl
             << "               in the format " << endl
             << "                  ACGTN" << endl
             << "                A abcde" << endl
             << "                C fghij" << endl
             << "                G klmno" << endl
             << "                T pqrst" << endl
             << "                N uvwxy" << " . The values a...y should be input as a quoted space separated " << endl
             << "               string: \"a b c ... y\". Lower scores are better, so matches should be less " << endl
             << "               than mismatches e.g. a,g,m,s = -5 (match), mismatch = 6. " << endl
             << "   -affineOpen value (10) " << endl
             << "               Set the penalty for opening an affine alignment." << endl
             << "   -affineExtend a (0)" << endl
             << "               Change affine (extension) gap penalty. Lower value allows more gaps." << endl << endl
             << " Options for overlap/dynamic programming alignments and pairwise overlap for de novo assembly. " << endl
             << "   -useQuality (false)" << endl
             << "               Use substitution/insertion/deletion/merge quality values to score gap and " << endl
             << "               mismatch penalties in pairwise alignments.  Because the insertion and deletion" << endl
             << "               rates are much higher than substitution, this will make many alignments " <<endl
             << "               favor an insertion/deletion over a substitution.  Naive consensus calling methods "<<endl
             << "               will then often miss substitution polymorphisms. This option should be " << endl
             << "               used when calling consensus using the Quiver method.  Furthermore, when " << endl
             << "               not using quality values to score alignments, there will be a lower consensus " << endl
             << "               accuracy in homolymer regions." << endl
             << "   -affineAlign (false)" << endl
             << "               Refine alignment using affine guided align." << endl << endl
             << " Options for filtering reads." << endl
             << "   -minReadLength l(50)" << endl
             << "               Skip reads that have a full length less than l. Subreads may be shorter." << endl 
             << "   -minSubreadLength l(0)" << endl
             << "               Do not align subreads of length less than l." << endl
             << "   -maxScore m (0)" << endl
             << "               Maximum score to output (high is bad, negative good)." << endl << endl
             << " Options for parallel alignment." << endl
             << "   -nproc N (1)" << endl
             << "               Align using N processes.  All large data structures such as the suffix array and " << endl
             << "               tuple count table are shared."<<endl
             << "   -start S (0)" << endl
             << "               Index of the first read to begin aligning. This is useful when multiple instances " << endl
             << "               are running on the same data, for example when on a multi-rack cluster."<<endl
             << "   -stride S (1)" << endl
             << "               Align one read every 'S' reads." << endl << endl
             << " Options for subsampling reads." << endl
             << "   -subsample (0)" << endl
             << "               Proportion of reads to randomly subsample (expressed as a decimal) and align." << endl
             << endl
//             << " Options for dynamic programming alignments. " << endl << endl
//             << "   -ignoreQuality" << endl
//             << "                 Ignore quality values when computing alignments (they still may be used." << endl 
//             << "                 when mapping)." << endl << endl
//             << " -v            Print some verbose information." << endl 
//             << " -V 2          Make verbosity more verbose.  Probably only useful for development." << endl
             << " -h            Print this help file." << endl << endl
             << "To cite BLASR, please use: Chaisson M.J., and Tesler G., Mapping " << endl
             << "single molecule sequencing reads using Basic Local Alignment with " << endl
             << "Successive Refinement (BLASR): Theory and Application, BMC " << endl
             << "Bioinformatics 2012, 13:238." << endl 
             << "Please report any bugs to "
             << "'https://github.com/PacificBiosciences/blasr/issues'." << endl << endl;
  str = helpStream.str();
}

void SetConciseHelp(string &conciseHelp) {
  stringstream strm;
  strm << "blasr - a program to map reads to a genome" << endl
       << " usage: blasr reads genome " << endl
       << " Run with -h for a list of commands " << endl
       << "          -help for verbose discussion of how to run blasr." << endl;
  conciseHelp = strm.str();
}

void PrintDiscussion() {
  cout << "NAME"<<endl;
  cout << "         blasr - Map SMRT Sequences to a reference genome."<< endl << endl;
  cout << "SYNOPSIS" << endl
       << "         blasr reads.fasta genome.fasta " << endl << endl
       << "         blasr reads.fasta genome.fasta -sa genome.fasta.sa" << endl << endl
       << "         blasr reads.bax.h5 genome.fasta [-sa genome.fasta.sa] " << endl << endl
       << "         blasr reads.bax.h5 genome.fasta -sa genome.fasta.sa -maxScore -100 -minMatch 15 ... " << endl << endl
       << "         blasr reads.bax.h5 genome.fasta -sa genome.fasta.sa -nproc 24 -out alignment.out ... " << endl << endl
       << "DESCRIPTION " << endl
       << "  blasr is a read mapping program that maps reads to positions " << endl
       << "  in a genome by clustering short exact matches between the read and" << endl
       << "  the genome, and scoring clusters using alignment. The matches are" << endl
       << "  generated by searching all suffixes of a read against the genome" << endl
       << "  using a suffix array. Global chaining methods are used to score " << endl
       << "  clusters of matches." << endl << endl
       << "  The only required inputs to blasr are a file of reads and a" << endl
       << "  reference genome.  It is exremely useful to have read filtering" << endl
       << "  information, and mapping runtime may decrease substantially when a" << endl
       << "  precomputed suffix array index on the reference sequence is" << endl
       << "  specified." << endl
       << "  " << endl
       << "  Although reads may be input in FASTA format, the recommended input is HDF" << endl
       << "  bax.h5 and plx.h5 files because these contain qualtiy value" << endl
       << "  information that is used in the alignment and produces higher quality" << endl
       << "  variant detection.  " << endl
       << "  " << endl
       << "  Read filtering information is contained in the .bax.h5 input files as" << endl
       << "  well as generated by other post-processing programs with analysis of" << endl
       << "  pulse files and read in from a separate .region.h5 file.  The current" << endl
       << "  set of filters that are applied to reads are high quality region" << endl
       << "  filtering, and adapter filtering.  Regions outside high-quality" << endl
       << "  regions are ignored in mapping.  Reads that contain regions annotated" << endl
       << "  as adapter are split into non-adapter (template) regions, and mapped" << endl
       << "  separately." << endl
       << "  " << endl
       << "  When suffix array index of a genome is not specified, the suffix array is" << endl
       << "  built before producing alignment.   This may be prohibitively slow" << endl
       << "  when the genome is large (e.g. Human).  It is best to precompute the" << endl
       << "  suffix array of a genome using the program sawriter, and then specify" << endl
       << "  the suffix array on the command line using -sa genome.fa.sa." << endl
       << "  " << endl
       << "  The optional parameters are roughly divided into three categories:" << endl
       << "  control over anchoring, alignment scoring, and output. " << endl
       << "  " << endl
       << "  The default anchoring parameters are optimal for small genomes and" << endl
       << "  samples with up to 5% divergence from the reference genome.  The main" << endl
       << "  parameter governing speed and sensitivity is the -minMatch parameter." << endl
       << "  For human genome alignments, a value of 11 or higher is recommended.  " << endl
       << "  Several methods may be used to speed up alignments, at the expense of" << endl
       << "  possibly decreasing sensitivity.  " << endl
       << "  " << endl
//       << "  If the genome is highly repetitive or divergent from the read" << endl
//       << "  sequences, the value of -maxExpand should be increased.  This option" << endl
//       << "  controls how much the search for anchors is expanded past a simple" << endl
//       << "  greedy search.  A value for -maxExpand of 1 is sufficent for" << endl
//       << "  non-repetitive genomes, and values of -maxExpand greater than 5 are" << endl
//       << "  not recommended." << endl
//       << "  " << endl
       << "  Regions that are too repetitive may be ignored during mapping by" << endl
       << "  limiting the number of positions a read maps to with the" << endl
       << "  -maxAnchorsPerPosition option.  Values between 500 and 1000 are effective" << endl
       << "  in the human genome." << endl
       << "  " << endl
       << "  For small genomes such as bacterial genomes or BACs, the default parameters " << endl
       << "  are sufficient for maximal sensitivity and good speed." << endl
       << endl << endl;
}

int CountZero(unsigned char *ptr, int length) {
  int i;
  int nZero = 0;
  for (i = 0; i < length; i++) {
    if (ptr[i] == 0) { ++nZero; }
  }
  return nZero;
}


bool ReadHasMeaningfulQualityValues(FASTQSequence &sequence) {
  if (sequence.qual.Empty() == true) {
    return 0;
  }
  else {
    int q;
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


void StoreRankingStats( vector<T_AlignmentCandidate*> &alignments,
                        VarianceAccumulator<float> &accumPValue, 
                        VarianceAccumulator<float> &accumWeight) {
  int i;
  for (i = 0; i < alignments.size(); i++) {
    alignments[i]->pvalVariance = accumPValue.GetVariance();
    alignments[i]->pvalNStdDev  = accumPValue.GetNStdDev(alignments[i]->clusterScore);
    alignments[i]->weightVariance = accumWeight.GetVariance();
    alignments[i]->weightNStdDev  = accumWeight.GetNStdDev(alignments[i]->clusterWeight);
  }

}



template<typename T_TargetSequence, typename T_QuerySequence, typename TDBSequence>
void AlignIntervals(T_TargetSequence &genome, T_QuerySequence &read, T_QuerySequence &rcRead,
                    WeightedIntervalSet &weightedIntervals,
                    int mutationCostMatrix[][5], 
                    int ins, int del, int sdpTupleSize,
                    int useSeqDB, SequenceIndexDatabase<TDBSequence> &seqDB,
                    vector<T_AlignmentCandidate*> &alignments,
                    MappingParameters &params,
                    int useScoreCutoff, int maxScore,
                    MappingBuffers &mappingBuffers,
                    int procId=0) {
                        
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
    alignment->clusterWeight= (*intvIt).size;
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
    DNALength subreadStart = read.subreadStart;
    DNALength subreadEnd   = read.subreadEnd;
    if ((*intvIt).GetStrandIndex() == Reverse) {
      subreadEnd   = read.MakeRCCoordinate(read.subreadStart) + 1;
      subreadStart = read.MakeRCCoordinate(read.subreadEnd-1);
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
      ((DNASequence)alignment->qAlignedSeq).PrintSeq(cout);
      cout << endl << "aligning reference" << endl;
      ((DNASequence)alignment->tAlignedSeq).PrintSeq(cout);
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

    int subreadLength = forrev[(*intvIt).GetStrandIndex()]->subreadEnd - forrev[(*intvIt).GetStrandIndex()]->subreadStart;
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
        if (alignment->tStrand == 0) {
          matches = (vector<ChainedMatchPos>*) &(*intvIt).matches;
        }
        else  { 
          rcMatches.resize((*intvIt).matches.size());
          //
          // Make the reverse complement of the match list.
          //
          
          // 1. Reverse complement the coordinates.
          for (m = 0; m < (*intvIt).matches.size(); m++) {
            int revCompIndex = rcMatches.size() - m - 1;
            rcMatches[revCompIndex].q = read.MakeRCCoordinate((*intvIt).matches[m].q + (*intvIt).matches[m].l - 1);
            rcMatches[revCompIndex].t = genome.MakeRCCoordinate((*intvIt).matches[m].t + (*intvIt).matches[m].l - 1);
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
        for (m = 0; matches->size() > 0 and m < matches->size() - 1; m++) {
          Block block;
          block.qPos = (*matches)[m].q - alignment->qAlignedSeqPos;
          block.tPos = (*matches)[m].t - alignment->tAlignedSeqPos;
          block.length = (*matches)[m].l;
          
          //
          // Find the lengths of the gaps between anchors.
          //
          int tGap, qGap;
          tGap = (*matches)[m+1].t - ((*matches)[m].t + (*matches)[m].l);
          qGap = (*matches)[m+1].q - ((*matches)[m].q + (*matches)[m].l);
          float gapRatio = (1.0*tGap)/qGap;

          if (tGap > 0 and qGap > 0) {

            DNASequence tSubSeq;
            FASTQSequence qSubSeq;
            DNALength tPos, qPos;
            tPos = block.tPos + block.length;
            qPos = block.qPos + block.length;
            tSubSeq.ReferenceSubstring(alignment->tAlignedSeq, tPos, tGap);
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
                                    params.extendFrontAlignment, params.recurseOver);
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
            alignment->blocks.push_back(block);
            anchorsOnly.blocks.push_back(block);
            alignment->blocks.insert(alignment->blocks.end(),
                                     alignmentInGap.blocks.begin(),
                                     alignmentInGap.blocks.end());
          }
        }

        m = (*matches).size() - 1;
        Block block;
        block.qPos = (*matches)[m].q - alignment->qAlignedSeqPos;
        block.tPos = (*matches)[m].t - alignment->tAlignedSeqPos;

        assert(block.tPos <= alignment->tAlignedSeq.length);
        assert(block.qPos <= alignment->qAlignedSeq.length);

        block.length = (*matches)[m].l;
        alignment->blocks.push_back(block);        
        anchorsOnly.blocks.push_back(block);
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
                //              SMRTDistanceMatrix, params.insertion, params.deletion);
      }
      else {
        alignScore = SDPAlign(alignment->qAlignedSeq, alignment->tAlignedSeq, distScoreFn, 
                              sdpTupleSize, params.sdpIns, params.sdpDel, params.indelRate*3, 
                              *alignment, mappingBuffers, 
                              Local, 
                              params.detailedSDPAlignment, 
                              params.extendFrontAlignment);
        ComputeAlignmentStats(*alignment, alignment->qAlignedSeq.seq, alignment->tAlignedSeq.seq,
                              distScoreFn);
                 //             SMRTDistanceMatrix, params.insertion, params.deletion);
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

          
      DNALength maximumExtendLength = 300;

      if (alignment->blocks.size() > 0 ){
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
            ((DNASequence)genome).CopyAsRC(genomeSuffix, lastAlignedTPos, genomeSuffixLength);
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
        genomeSuffix.Free();

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
            ((DNASequence)genome).MakeRC(genomePrefix, firstAlignedTPos - genomePrefixLength, genomePrefixLength);
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
        }
        else {
          genomePrefix.Free();
        }
      }
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


bool FirstContainsSecond(DNALength aStart, DNALength aEnd, DNALength bStart, DNALength bEnd) {
  return ((bStart > aStart and bEnd <= aEnd) or
          (bStart >= aStart and bEnd < aEnd));
}

template<typename T_Sequence>
bool CheckForSufficientMatch(T_Sequence &read, vector<T_AlignmentCandidate*> &alignmentPtrs, MappingParameters &params) {
  if (alignmentPtrs.size() > 0 and alignmentPtrs[0]->score < params.maxScore) {
    return true;
  }
  else {
    return false;
  }
}

void DeleteAlignments(vector<T_AlignmentCandidate*> &alignmentPtrs, int startIndex=0) {
  int i;
  for (i = startIndex; i < alignmentPtrs.size(); i++ ) {
    delete alignmentPtrs[i];
  }
  alignmentPtrs.resize(0);
}

int RemoveLowQualitySDPAlignments(int readLength, vector<T_AlignmentCandidate*> &alignmentPtrs, MappingParameters &params) {
  // Just a hack.  For now, assume there is at least 1 match per 50 bases.
  int totalBasesMatched = 0;
  int a;
  for (a = 0; a < alignmentPtrs.size(); a++) {
    int b;
    for (b = 0; b < alignmentPtrs[a]->blocks.size(); b++) {
      totalBasesMatched += alignmentPtrs[a]->blocks[b].length;
    }
    int expectedMatches = params.sdpTupleSize/50.0 * readLength;
    if (totalBasesMatched < expectedMatches) {
      delete alignmentPtrs[a];
      alignmentPtrs[a] = NULL;
    }
  }
  int packedAlignmentIndex = 0;
  for (a = 0; a < alignmentPtrs.size(); a++) {
    if (alignmentPtrs[a] != NULL) {
      alignmentPtrs[packedAlignmentIndex] = alignmentPtrs[a];
      packedAlignmentIndex++;
    }
  }
  alignmentPtrs.resize(packedAlignmentIndex);
  return packedAlignmentIndex;
}


template<typename T_Sequence>
int RemoveLowQualityAlignments(T_Sequence &read, vector<T_AlignmentCandidate*> &alignmentPtrs, MappingParameters &params) {
  if (params.verbosity > 0) {
    cout << "checking at least " << alignmentPtrs.size() << " alignments to see if they are accurate." << endl;
  }
  UInt i;
  for (i = 0; i < MIN(params.nCandidates, alignmentPtrs.size()); i++) { 
    if (params.verbosity > 0) {
      cout << "  " << i << " " << alignmentPtrs[i]->score << endl;
    }
    if (alignmentPtrs[i]->blocks.size() == 0 or
        alignmentPtrs[i]->score > params.maxScore) {
      //
      // Since the alignments are sorted according to alignment
      // score, once one of the alignments is too low of a score,
      // all remaining alignments are also too low, and should be
      // removed as well.  Do that all at once.
      //
      if (alignmentPtrs[i]->blocks.size() == 0 and params.verbosity > 0) {
        cout << "Removing empty alignment " << alignmentPtrs[i]->qName << endl;
      }
      if (params.verbosity  > 0) { 
        cout << alignmentPtrs[i]->qName << " alignment " << i << " is too low of a score." << alignmentPtrs[i]->score << endl;
      }
      int deletedIndex = i;
      for (; deletedIndex < alignmentPtrs.size(); deletedIndex++) {
        delete alignmentPtrs[deletedIndex];
        alignmentPtrs[deletedIndex] = NULL;
      }
      alignmentPtrs.erase(i + alignmentPtrs.begin(), alignmentPtrs.end());
      break;
    }
    else {
      if (params.verbosity > 0) {
        cout << "Keeping alignment " << i << " " << alignmentPtrs[i]->qPos << " " << alignmentPtrs[i]->qLength
             << " " << alignmentPtrs[i]->tName << " " << alignmentPtrs[i]->tPos << " " << alignmentPtrs[i]->tLength 
             << " from score: " << alignmentPtrs[i]->score << endl;
      }
    }
  }
  return alignmentPtrs.size();
}

int RemoveOverlappingAlignments(vector<T_AlignmentCandidate*> &alignmentPtrs, MappingParameters &params) {
  vector<unsigned char> alignmentIsContained;
  alignmentIsContained.resize(alignmentPtrs.size());
  std::fill(alignmentIsContained.begin(), alignmentIsContained.end(), false);

  int j;
  int numContained = 0;
  int curNotContained = 0;
    
  if (alignmentPtrs.size() > 0) {
    UInt i;
    for (i = 0; i < alignmentPtrs.size()-1; i++ ){
      T_AlignmentCandidate *aref = alignmentPtrs[i];
      if (aref->pctSimilarity < params.minPctIdentity) {
        continue;
      }
      for (j = i + 1; j < alignmentPtrs.size(); j++ ){
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
            cout << "alignment " << i << " is contained in " << j << endl;
            cout << aref->tAlignedSeqPos << " " <<  alignmentPtrs[j]->tAlignedSeqPos << " "
                 << aref->tAlignedSeqPos + aref->tAlignedSeqLength << " " 
                 << alignmentPtrs[j]->tAlignedSeqPos + alignmentPtrs[j]->tAlignedSeqLength << endl;
          }
        }
        else if (alignmentPtrs[j]->GenomicTBegin() <= aref->GenomicTBegin() and
                 alignmentPtrs[j]->GenomicTEnd()   >= aref->GenomicTEnd() and 
                 alignmentPtrs[i]->tIndex == alignmentPtrs[j]->tIndex) {
          if (params.verbosity >= 2) {
            cout << "ALIGNMENT " << j << " is contained in " << i << endl;
            cout << alignmentPtrs[j]->tAlignedSeqPos << " " <<  aref->tAlignedSeqPos << " "
                 << alignmentPtrs[j]->tAlignedSeqPos + alignmentPtrs[j]->tAlignedSeqLength << " " 
                 <<  aref->tAlignedSeqPos + aref->tAlignedSeqLength << endl;
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
      }
      else {
        alignmentPtrs[curNotContained] = aref;
        ++curNotContained;
      }
    }
    alignmentPtrs.resize(alignmentPtrs.size() - numContained);
  } 
  return alignmentPtrs.size();
}

template<typename T_RefSequence, typename T_Sequence>
void RefineAlignments(vector<T_Sequence*> &bothQueryStrands,
                      T_RefSequence &genome,
                      vector<T_AlignmentCandidate*> &alignmentPtrs, MappingParameters &params, MappingBuffers &mappingBuffers) {

  
  UInt i;
  for (i = 0; i < alignmentPtrs.size(); i++ ) {
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
    

template<typename T_RefSequence, typename T_Sequence>
void PairwiseLocalAlign(T_Sequence &qSeq, T_RefSequence &tSeq, 
                        int k, 
                        MappingParameters &params, T_AlignmentCandidate &alignment, 
                        MappingBuffers &mappingBuffers,
                        AlignmentType alignType=Global) {
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
          //SMRTDistanceMatrix, params.indel, params.indel );

  if (params.scoreType == 1) {
    alignment.score = alignment.sumQVScore;
  }
  
}

template<typename T_RefSequence, typename T_Sequence>
void RefineAlignment(vector<T_Sequence*> &bothQueryStrands,
                     T_RefSequence &genome,
                     T_AlignmentCandidate  &alignmentCandidate, MappingParameters &params,
                     MappingBuffers &mappingBuffers) {


  FASTQSequence qSeq;
  DNASequence   tSeq;
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
    SMRTSequence subread;
    subread.ReferenceSubstring(*bothQueryStrands[0], 
                               bothQueryStrands[0]->subreadStart,
                               (bothQueryStrands[0]->subreadEnd - 
                                bothQueryStrands[0]->subreadStart));

    int drift = ComputeDrift(alignmentCandidate);
    T_AlignmentCandidate refinedAlignment;

    KBandAlign(subread, alignmentCandidate.tAlignedSeq, SMRTDistanceMatrix, 
               params.insertion, params.deletion,
                              drift,
                              mappingBuffers.scoreMat, mappingBuffers.pathMat,
                              refinedAlignment, idsScoreFn, Global);
    refinedAlignment.RemoveEndGaps();
    ComputeAlignmentStats(refinedAlignment, 
                          subread.seq, 
                          alignmentCandidate.tAlignedSeq.seq, 
                          distScoreFn2);
                          //idsScoreFn);
                          //SMRTDistanceMatrix, params.indel, params.indel);
    
    alignmentCandidate.blocks = refinedAlignment.blocks;
    alignmentCandidate.gaps   = refinedAlignment.gaps;
    alignmentCandidate.tPos   = refinedAlignment.tPos;
    alignmentCandidate.qPos   = refinedAlignment.qPos + bothQueryStrands[0]->subreadStart;
    alignmentCandidate.score  = refinedAlignment.score;
  }
  else if (params.useGuidedAlign) {
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
      tSeq.Copy(alignmentCandidate.tAlignedSeq, 
                alignmentCandidate.tPos,
                (alignmentCandidate.blocks[lastBlock].tPos + 
                 alignmentCandidate.blocks[lastBlock].length -
                 alignmentCandidate.blocks[0].tPos));
    
      //      qSeq.ReferenceSubstring(alignmentCandidate.qAlignedSeq,
      qSeq.ReferenceSubstring(*bothQueryStrands[0],
                              alignmentCandidate.qAlignedSeqPos + alignmentCandidate.qPos, 
                              (alignmentCandidate.blocks[lastBlock].qPos +
                               alignmentCandidate.blocks[lastBlock].length));


      if (!params.ignoreQualities && ReadHasMeaningfulQualityValues(alignmentCandidate.qAlignedSeq)) {
        if (params.affineAlign) {
            AffineGuidedAlign(qSeq, tSeq, alignmentCandidate, 
                            idsScoreFn, params.bandSize,
                            mappingBuffers, 
                            refinedAlignment, Global, false);
        }
        else {
            GuidedAlign(qSeq, tSeq, alignmentCandidate, 
                        idsScoreFn, params.guidedAlignBandSize,
                        mappingBuffers, 
                        refinedAlignment, Global, false);
        }
      }
      else {
        if (params.affineAlign) {
            AffineGuidedAlign(qSeq, tSeq, alignmentCandidate, 
                              distScoreFn, params.bandSize,
                            mappingBuffers, 
                            refinedAlignment, Global, false);
        }
        else {
        GuidedAlign(qSeq, tSeq, alignmentCandidate, 
                    distScoreFn, params.guidedAlignBandSize,
                    mappingBuffers,
                    refinedAlignment, Global, false);
        }
      }
      ComputeAlignmentStats(refinedAlignment, 
                            qSeq.seq,
                            tSeq.seq, 
                            distScoreFn2, params.affineAlign);
                            //SMRTDistanceMatrix, params.indel, params.indel);
      //
      // Copy the refine alignment, which may be a subsequence of the
      // alignmentCandidate into the alignment candidate.  
      //

      // First copy the alignment block and gap (the description of
      // the base by base alignment).

      alignmentCandidate.blocks.clear();
      alignmentCandidate.blocks = refinedAlignment.blocks;

      alignmentCandidate.CopyStats(refinedAlignment);

      alignmentCandidate.gaps   = refinedAlignment.gaps;
      alignmentCandidate.score  = refinedAlignment.score;
      alignmentCandidate.nCells = refinedAlignment.nCells;

      // Next copy the information that describes what interval was
      // aligned.  Since the reference sequences of the alignment
      // candidate have been modified, they are reassigned.
      alignmentCandidate.tAlignedSeq.TakeOwnership(tSeq);
      alignmentCandidate.ReassignQSequence(qSeq);
      alignmentCandidate.tAlignedSeqPos    += alignmentCandidate.tPos; 
      alignmentCandidate.qAlignedSeqPos    += alignmentCandidate.qPos;

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
  }
  else {


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

    if (alignmentCandidate.blocks.size() == 0 ){ 
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
                              alignmentCandidate.blocks[lastBlock].qPos + alignmentCandidate.blocks[lastBlock].length);
    }
    else {
      qSeq.ReferenceSubstring(alignmentCandidate.qAlignedSeq, // the subsequence that the alignment points to
                              alignmentCandidate.qPos  + alignmentCandidate.qAlignedSeqPos, 
                              alignmentCandidate.blocks[lastBlock].qPos + alignmentCandidate.blocks[lastBlock].length - alignmentCandidate.blocks[0].qPos);
    }
      
    tSeq.Copy(alignmentCandidate.tAlignedSeq, // the subsequence the alignment points to
              alignmentCandidate.tPos, // ofset into the subsequence
              alignmentCandidate.blocks[lastBlock].tPos + alignmentCandidate.blocks[lastBlock].length - alignmentCandidate.blocks[0].tPos);

    T_AlignmentCandidate refinedAlignment;

    //
    // When the parameter bandSize is 0, set the alignment band size
    // to the drift off the diagonal, plus a little more for wiggle
    // room.  When the parameteris nonzero, use that as a fixed band.
    //
    int k;
    if (params.bandSize == 0) {
      k = abs(drift) * 1.5;
    }
    else {
      k = params.bandSize;
    }
    if (params.verbosity > 0) {
      cout << "drift: " << drift << " qlen: " << alignmentCandidate.qAlignedSeq.length << " tlen: " << alignmentCandidate.tAlignedSeq.length << " k: " << k << endl;
      cout << "aligning in " << k << " * " << alignmentCandidate.tAlignedSeq.length << " " << k * alignmentCandidate.tAlignedSeq.length << endl;
    }
    if (k < 10) {
      k = 10;
    }

    alignmentCandidate.tAlignedSeqPos    += alignmentCandidate.tPos; 
    
    VectorIndex lastSDPBlock = alignmentCandidate.blocks.size() - 1;

    if (alignmentCandidate.blocks.size() > 0) {
      DNALength prevLength =  alignmentCandidate.tAlignedSeqLength -= alignmentCandidate.tPos;
      alignmentCandidate.tAlignedSeqLength = (alignmentCandidate.blocks[lastSDPBlock].tPos 
                                              + alignmentCandidate.blocks[lastSDPBlock].length 
                                              - alignmentCandidate.blocks[0].tPos);
    }
    else {
      alignmentCandidate.tAlignedSeqLength = 0;
    }

    alignmentCandidate.tPos              = 0;
    alignmentCandidate.qAlignedSeqPos    += alignmentCandidate.qPos;

    if (alignmentCandidate.blocks.size() > 0) {
      DNALength prevLength =  alignmentCandidate.qAlignedSeqLength -= alignmentCandidate.qPos; 
      alignmentCandidate.qAlignedSeqLength = (alignmentCandidate.blocks[lastSDPBlock].qPos 
                                              + alignmentCandidate.blocks[lastSDPBlock].length
                                              - alignmentCandidate.blocks[0].qPos);
    }
    else {
      alignmentCandidate.qAlignedSeqLength = 0;
    }
    alignmentCandidate.qPos                = 0;

    alignmentCandidate.blocks.clear();
    alignmentCandidate.tAlignedSeq.Free();
    alignmentCandidate.tAlignedSeq.TakeOwnership(tSeq);
    alignmentCandidate.ReassignQSequence(qSeq);

    if (params.verbosity >= 2) {
      cout << "refining target: " << endl;
      alignmentCandidate.tAlignedSeq.PrintSeq(cout);
      cout << "refining query: " << endl;
      ((DNASequence*)&alignmentCandidate.qAlignedSeq)->PrintSeq(cout);
      cout << endl;
    }
    PairwiseLocalAlign(qSeq, tSeq, k, params, alignmentCandidate, mappingBuffers, Fit);
  }
}

void AssignRefContigLocation(T_AlignmentCandidate &alignment, SequenceIndexDatabase<FASTQSequence> &seqdb, DNASequence &genome) {
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
  }
  else {
    //
    // Flip coordinates into forward strand in order to find the boundaries 
    // of the contig, then reverse them in order to find offset.
    //

    // Find the reverse complement coordinate of the index of the last aligned base.
    assert(alignment.tAlignedSeqLength > 0);
    forwardTPos = genome.MakeRCCoordinate(alignment.tAlignedSeqPos + alignment.tAlignedSeqLength - 1);
    seqDBIndex  = seqdb.SearchForIndex(forwardTPos);

    
    //
    // Find the reverse comlement coordinate of the last base of this
    // sequence.  This would normally be the start of the next contig
    // -1 to get the length, but since an 'N' is added between every
    // pair of sequences, this is -2.
    //
    DNALength reverseTOffset;
    reverseTOffset = genome.MakeRCCoordinate(seqdb.seqStartPos[seqDBIndex+1]-2);
    alignment.tAlignedSeqPos -= reverseTOffset;
  }
}

void AssignRefContigLocations(vector<T_AlignmentCandidate*> &alignmentPtrs, SequenceIndexDatabase<FASTQSequence> &seqdb, DNASequence &genome) {
  
  UInt i;
  for (i = 0; i < alignmentPtrs.size(); i++) {
    T_AlignmentCandidate *aref = alignmentPtrs[i];
    AssignRefContigLocation(*aref, seqdb, genome);
  }
}


template<typename T_RefSequence>
void AssignGenericRefContigName(vector<T_AlignmentCandidate*> &alignmentPtrs, T_RefSequence &genome) {
  UInt i;
  for (i = 0; i < alignmentPtrs.size(); i++) {
    T_AlignmentCandidate *aref = alignmentPtrs[i];
    aref->tName = genome.title;
  }
}


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
             MappingIPC *mapData) {

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
      numKeysMatched   = MapReadToGenome(bwt, read, read.subreadStart, read.subreadEnd, 
                                         mappingBuffers.matchPosList, params.anchorParameters, forwardNumBasesMatched);
      if (!params.forwardOnly) {
        rcNumKeysMatched = MapReadToGenome(bwt, readRC, readRC.subreadStart, readRC.subreadEnd, 
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
                                (DNALength) ((read.subreadEnd - read.subreadStart) * (1 + params.indelRate)), params.nCandidates,
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
                                (DNALength) ((read.subreadEnd - read.subreadStart) * (1 + params.indelRate)), params.nCandidates, 
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
                                (DNALength) ((read.subreadEnd - read.subreadStart) * (1 + params.indelRate)), params.nCandidates,
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
                                (DNALength) ((read.subreadEnd - read.subreadStart) * (1 + params.indelRate)), params.nCandidates, 
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
                                (DNALength) ((read.subreadEnd - read.subreadStart) * (1 + params.indelRate)), params.nCandidates,
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
                                (DNALength) ((read.subreadEnd - read.subreadStart) * (1 + params.indelRate)), params.nCandidates, 
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
                    params.useScoreCutoff, params.maxScore,
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

  if (params.forPicard) {
    int a;
    for (a = 0; a < alignmentPtrs.size(); a++ ) {
      alignmentPtrs[a]->OrderGapsByType();
    }
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
    //    aref->qStrand = aref->readIndex;
    if (aref->tStrand == 0) {
      aref->qName = read.GetName();
    }
    else {
      aref->qName = readRC.GetName();
    }
  }

  AssignRefContigLocations(alignmentPtrs, seqdb, genome);
}


void SumMismatches(SMRTSequence &read,
                   T_AlignmentCandidate &alignment,
                   int fullIntvStart, int fullIntvEnd,
                   int &sum) {
  int alnStart, alnEnd;
  alignment.GetQIntervalOnForwardStrand(alnStart, alnEnd);
  int p;
  sum = 0;
  if (read.substitutionQV.Empty() == false) {
    for (p = fullIntvStart; p < alnStart; p++) {
      sum += read.substitutionQV[p];
    }
    for (p = alnEnd; p < fullIntvEnd; p++) {
      sum += read.substitutionQV[p];
    }
  }
}

int FindMaxLengthAlignment(vector<T_AlignmentCandidate*> alignmentPtrs,
                           int &maxLengthIndex) {
  int i;
  int maxLength = 0;
  maxLengthIndex = -1;

  for (i = 0; i < alignmentPtrs.size(); i++) {
    int qStart, qEnd;
    alignmentPtrs[i]->GetQInterval(qStart, qEnd);
    if (qEnd - qStart > maxLength) {
      maxLengthIndex = i;
      maxLength = qEnd - qStart;
    }
  }
  return (maxLength != -1);
}

bool AlignmentsOverlap(T_AlignmentCandidate &alnA, T_AlignmentCandidate &alnB, float minPercentOverlap) {
  int alnAStart, alnAEnd, alnBStart, alnBEnd;
  bool useForwardStrand=true;
  alnA.GetQInterval(alnAStart, alnAEnd, useForwardStrand);
  alnB.GetQInterval(alnBStart, alnBEnd, useForwardStrand);
  // Look if one alignment encompasses the other
  int ovp = 0;
  if (alnAStart <= alnBStart and alnAEnd >= alnBEnd) {
    ovp = alnBEnd - alnBStart;
  }
  else if (alnBStart <= alnAStart and alnBEnd >= alnAEnd) {
    ovp = alnAEnd - alnAStart;
  }
  else {
    //
    // Look to see if the alignments overlap
    //

    if (alnAEnd >= alnBStart and alnAEnd <= alnBEnd) {
      ovp = alnAEnd - alnBStart;
    }
    else if (alnAStart >= alnBStart and alnAStart <= alnBEnd) {
      ovp = alnBEnd - alnAStart;
    }
  }
  float ovpPercent = (2.0*ovp) / ((alnAEnd - alnAStart) + (alnBEnd - alnBStart));
  // returns true when an overlap is found.
  return ( ovpPercent > minPercentOverlap);

}


void PartitionOverlappingAlignments(vector<T_AlignmentCandidate*> &alignmentPtrs,
                                    vector<set<int> > &partitions,
                                    float minOverlap) {
  if (alignmentPtrs.size() == 0) {
    partitions.clear();
    return;
  }

  set<int>::iterator setIt, setEnd;
  int i, p;
  bool overlapFound = false;
  for (i = 0; i < alignmentPtrs.size(); i++) {
    overlapFound = false;
    for (p = 0; p < partitions.size() and overlapFound == false; p++) {
      setEnd = partitions[p].end();
      for (setIt = partitions[p].begin(); setIt != partitions[p].end() and overlapFound == false; ++setIt) {
        if (AlignmentsOverlap(*alignmentPtrs[i], *alignmentPtrs[*setIt], minOverlap) or
            (alignmentPtrs[i]->QAlignStart() <= alignmentPtrs[*setIt]->QAlignStart()) and
            (alignmentPtrs[i]->QAlignEnd()  > alignmentPtrs[*setIt]->QAlignEnd())) {
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
      partitions.push_back(set<int>());
      partitions[partitions.size()-1].insert(i);
    }
  }
}

void ScaleMapQVByClusterSize(T_AlignmentCandidate &alignment, MappingParameters &params) {
  if (alignment.numSignificantClusters > params.nCandidates) {
    alignment.mapQV = Phred((1-InversePhred(alignment.mapQV))* ((float)params.nCandidates / alignment.numSignificantClusters));
  }
  else if (alignment.numSignificantClusters == 0) {
    alignment.mapQV = 0;
  }
}

void StoreMapQVs(SMRTSequence &read,
                 vector<T_AlignmentCandidate*> &alignmentPtrs, 
                 MappingParameters &params) {
  
  //
  // Only weight alignments for mapqv against eachother if they are overlapping.
  //
  int a;
  vector<set<int> > partitions; // Each set contains alignments that overlap on the read.
  DistanceMatrixScoreFunction<DNASequence, FASTQSequence> distScoreFn;
  distScoreFn.del = params.deletion;
  distScoreFn.ins = params.insertion;
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
  for (a = 0; a < alignmentPtrs.size(); a++) {
    if (params.ignoreQualities == false) {
      alignmentPtrs[a]->probScore = -ComputeAlignmentScore(*alignmentPtrs[a],
                                                          alignmentPtrs[a]->qAlignedSeq,
                                                          alignmentPtrs[a]->tAlignedSeq,
                                                          idsScoreFn) / 10.0;
    }
    else {
      alignmentPtrs[a]->probScore = -ComputeAlignmentScore(*alignmentPtrs[a],
                                                           alignmentPtrs[a]->qAlignedSeq,
                                                           alignmentPtrs[a]->tAlignedSeq,
                                                           distScoreFn) / 10.0;
    }
  }
  PartitionOverlappingAlignments(alignmentPtrs, partitions, params.minFractionToBeConsideredOverlapping);
  
  int p;
  set<int>::iterator partIt, partEnd;
  
  //
  // For each partition, store where on the read it begins, and where
  // it ends. 
  //
  vector<int> partitionBeginPos, partitionEndPos;
  partitionBeginPos.resize(partitions.size());
  partitionEndPos.resize(partitions.size());
  fill(partitionBeginPos.begin(), partitionBeginPos.end(), -1);
  fill(partitionEndPos.begin(), partitionEndPos.end(), -1);
  vector<char> assigned;
  assigned.resize( alignmentPtrs.size());
  fill(assigned.begin(), assigned.end(), false);
                   
  for (p = 0; p < partitions.size(); p++) {
    partEnd = partitions[p].end();
    int alnStart, alnEnd;

    if (partitions[p].size() > 0) {
      partIt = partitions[p].begin();
      alignmentPtrs[*partIt]->GetQInterval(alnStart, alnEnd);
      partitionBeginPos[p] = alnStart; 
      partitionEndPos[p]   = alnEnd;
      ++partIt;
      partEnd = partitions[p].end();
      for (; partIt != partEnd; ++partIt) {
        //  Comment out because all reads are now in the forward strand.
        //  alignmentPtrs[*partIt]->GetQInterval(alnStart, alnEnd, convertToForwardStrand);
        alignmentPtrs[*partIt]->GetQInterval(alnStart, alnEnd);
        if (alnEnd - alnStart > partitionEndPos[p] - partitionBeginPos[p]) {
          partitionBeginPos[p] = alnStart;
          partitionEndPos[p]   = alnEnd;
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

  for (p = 0; p < partitions.size(); p++) {
    partEnd = partitions[p].end();
    int alnStart, alnEnd;
    for (partIt = partitions[p].begin(); partitions[p].size() > 0 and partIt != partEnd; ++partIt) {
      int frontSum = 0, endSum = 0;
      int alnStart, alnEnd;
      int mismatchSum = 0;
      alignmentPtrs[*partIt]->GetQInterval(alnStart, alnEnd, convertToForwardStrand);
      if (alnStart - partitionBeginPos[p] > MAPQV_END_ALIGN_WIGGLE or
          partitionEndPos[p] - alnEnd > MAPQV_END_ALIGN_WIGGLE) {
        SumMismatches(read, *alignmentPtrs[*partIt], partitionBeginPos[p], partitionEndPos[p], mismatchSum);
      }
      //
      // Random sequence can be aligned with about 50% similarity due
      // to optimization, so weight the qv sum 
      //
      alignmentPtrs[*partIt]->probScore += -(mismatchSum) * 0.5;
    }
  }

  //                                           
  // Determine mapqv by summing qvscores in partitions

  float mapQVDenominator = 0;
  for (p = 0; p < partitions.size(); p++) {
    set<int>::iterator nextIt;
    if (partitions[p].size() == 0) {
      continue;
    }
    int index = *partitions[p].begin();

    mapQVDenominator = alignmentPtrs[index]->probScore;

    if (partitions[p].size() > 1) {
      partIt  = partitions[p].begin();
      partEnd = partitions[p].end();
      ++partIt;

      for (; partIt != partEnd; ++partIt) {
        index = *partIt;
        mapQVDenominator = LogSumOfTwo(mapQVDenominator, alignmentPtrs[index]->probScore);
      }
    }
    
    
    for (partIt = partitions[p].begin(); 
         partIt != partitions[p].end(); ++partIt) {
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
      }
      else {
        double log10 = log(10);
        double sub   = alignmentPtrs[*partIt]->probScore - mapQVDenominator;
        double expo = exp(log10*sub);
        double diff = 1.0 - expo;
        int phredValue;
        
        if (expo == 0) {
          phredValue = 0;
        }
        else if (diff == 0) {
          phredValue = MAX_PHRED_SCORE;
        }
        else {
          phredValue = Phred(diff);
        }
        if (phredValue > MAX_PHRED_SCORE) {
          phredValue = MAX_PHRED_SCORE;
        }
        
        alignmentPtrs[*partIt]->mapQV = phredValue;
        assigned[*partIt]=true;
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

//
// The full read is not the subread, and does not have masked off characters.
//
void PrintAlignment(T_AlignmentCandidate &alignment, SMRTSequence &fullRead, MappingParameters &params, AlignmentContext &alignmentContext, ostream &outFile) {
  if (alignment.score > params.maxScore) {
    return;
  }
  if (alignment.pctSimilarity < params.minPctIdentity) {
    return;
  }
  try {
    int lastBlock = alignment.blocks.size() - 1;
    if (params.printFormat == StickPrint) {
      PrintAlignmentStats(alignment, outFile);
      StickPrintAlignment(alignment,
                          (DNASequence&) alignment.qAlignedSeq,
                          (DNASequence&) alignment.tAlignedSeq,
                          outFile,
                          alignment.qAlignedSeqPos, alignment.tAlignedSeqPos);
    }
    else if (params.printFormat == SAM) {
      SAMOutput::PrintAlignment(alignment, fullRead, outFile, alignmentContext, params.clipping);
    }
    else if (params.printFormat == CompareXML) {
      CompareXMLPrintAlignment(alignment,
                               (DNASequence&) alignment.qAlignedSeq, (DNASequence&) alignment.tAlignedSeq,
                               outFile,
                               alignment.qAlignedSeqPos, alignment.tAlignedSeqPos);
    }
    else if (params.printFormat == Vulgar) {
      PrintAlignmentStats(alignment, outFile);
      VulgarPrintAlignment(alignment, outFile);
    }
    else if (params.printFormat == CompareSequencesParsable) {
      PrintCompareSequencesAlignment(alignment, alignment.qAlignedSeq, alignment.tAlignedSeq, outFile);
    }
    else if (params.printFormat == Interval) {
      if (alignment.blocks.size() > 0) {
        IntervalAlignmentPrinter::Print(alignment, outFile);
      }
    }
    else if (params.printFormat == SummaryPrint) {
      if (alignment.blocks.size() > 0) {
        SummaryAlignmentPrinter::Print(alignment, outFile);
      }
    }
  }
  catch (ostream::failure f) {
    cout << "ERROR writing to output file. The output drive may be full, or you  " << endl;
    cout << "may not have proper write permissions." << endl;
    exit(1);
  }
}

vector<T_AlignmentCandidate*>
SelectAlignmentsToPrint(vector<T_AlignmentCandidate*> alignmentPtrs,
                        MappingParameters & params,
                        const int & associatedRandInt) {
  //
  // Select all alignments, unless the parameter placeRepeatsRandomly is set.
  // In this case only one read is selected and it is selected from all
  // equally top scoring hits.
  //
  UInt i;
  int optScore;
  int nOpt = 0;
  int startIndex = 0;
  int endIndex = 0;
 
  if (params.placeRandomly) {
    if (alignmentPtrs.size() > 0) {
      std::sort(alignmentPtrs.begin(), alignmentPtrs.end(), 
                SortAlignmentPointersByScore());
      optScore = alignmentPtrs[0]->score;
      nOpt = 1;
      // First find the minimum score, and count how many times it
      // exists
      for (i = 1; i < alignmentPtrs.size(); i++) { 
        assert(alignmentPtrs[i]->score >= optScore);
        if (alignmentPtrs[i]->score == optScore) {
          nOpt++;
        }
      }

      // Select a random item from all equally top scoring hits:
      // alignmentPtrs[0...nOpt-1].
      startIndex = associatedRandInt % nOpt;
      endIndex = startIndex + 1;
    }
  } else {
    startIndex = 0;
    endIndex   = MIN(params.nBest, alignmentPtrs.size());
  }

  vector<T_AlignmentCandidate*> ret;
  ret.assign(alignmentPtrs.begin() + startIndex, 
             alignmentPtrs.begin() + endIndex);
  return ret;
}

// Print all alignments in vector<T_AlignmentCandidate*> alignmentPtrs
void PrintAlignments(vector<T_AlignmentCandidate*> alignmentPtrs,
                     SMRTSequence &read,
                     MappingParameters &params, ostream &outFile, 
                     AlignmentContext alignmentContext) {
  if (params.nProc > 1) {
#ifdef __APPLE__
    sem_wait(semaphores.writer);
#else
    sem_wait(&semaphores.writer);
#endif
  }
  for (int i = 0; i < alignmentPtrs.size(); i++) { 
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
      continue;
    }
    
    //
    // Configure some of the alignment context before printing.
    //
    if (i > 0 and params.placeRandomly == false) {
      alignmentContext.isPrimary = false;
    }
    else {
      alignmentContext.isPrimary = true;
    }

    if (params.printSAM) {
        DistanceMatrixScoreFunction<DNASequence, FASTASequence> editdistScoreFn(EditDistanceMatrix, 1, 1);
        T_AlignmentCandidate & alignment = *alignmentPtrs[i];
        alignmentContext.editDist = ComputeAlignmentScore(alignment,
            alignment.qAlignedSeq, 
            alignment.tAlignedSeq, 
            editdistScoreFn);
            //EditDistanceMatrix, 1, 1);
    }
    
    PrintAlignment(*alignmentPtrs[i], read, params, alignmentContext, outFile);
  }

  if (params.nProc > 1) {
#ifdef __APPLE__
    sem_post(semaphores.writer);
#else
    sem_post(&semaphores.writer);
#endif
  }

}

template<typename T_Sequence>
bool GetNextReadThroughSemaphore(ReaderAgglomerate &reader, MappingParameters &params, T_Sequence &read, AlignmentContext &context, int & associatedRandInt) {
  //
  // Grab the value of the semaphore for debugging purposes.
  //
  // uncomment when needed
  // int semvalue;
  // if (params.nProc > 1) {
  //  sem_getvalue(&semaphores.reader, &semvalue);
  // }

  //
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
  context.readGroupId = reader.readGroupId;
  
  if (params.nProc > 1) {
#ifdef __APPLE__
    sem_post(semaphores.reader);
#else
    sem_post(&semaphores.reader);
#endif
  }
  return returnValue;
}

void AssignMapQV(vector<T_AlignmentCandidate*> &alignmentPtrs) {
  int i;
  int mapQV = 1;
  if (alignmentPtrs.size() > 1 and alignmentPtrs[0]->score == alignmentPtrs[1]->score) {
    // the top two alignments have the same score, don't consider them as mapped.
    mapQV = 0;
  }
  
  for (i = 0; i < alignmentPtrs.size(); i++) {
    alignmentPtrs[i]->mapQV = mapQV;
  }
}

void PrintAlignmentPtrs(vector <T_AlignmentCandidate*> & alignmentPtrs,
    ostream & out = cout) {
    for(int alignmentIndex = 0; 
        alignmentIndex < alignmentPtrs.size();
        alignmentIndex++) {
        out << "["<< alignmentIndex << "/" 
            << alignmentPtrs.size() << "]" << endl;
        T_AlignmentCandidate *alignment = alignmentPtrs[alignmentIndex];          
        alignment->Print(out);
    }
    out << endl;
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
  newTAlignedSeqPos = UInt((alignment->tAlignedSeqPos > flankSize)?
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
    ((DNASequence) subread).PrintSeq(threadOut);
    threadOut << endl;
    threadOut << "alignedRefSeq: " << endl;
    ((DNASequence) alignedRefSequence).PrintSeq(threadOut);
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


// Given a SMRT sequence and a subread interval, make the subread.
// Input: 
//   smrtRead         - a SMRT sequence 
//   subreadInterval  - a subread interval 
//   params           - mapping parameters
// Output: 
//   subreadSequence - the constructed subread
void MakeSubreadOfInterval(SMRTSequence & subreadSequence,
    SMRTSequence & smrtRead, ReadInterval & subreadInterval, 
    MappingParameters & params) {
    //
    // subreadMapType is a way of limiting the portion of the read
    // that is aligned.  The output is similar, but the
    // computation is slightly different.  The subreadMapType 0
    // was written first, and just creates a hard mask over the
    // regions that are not to be aligned.  The subreadMapType is
    // slightly more formal mode where a new read is pointed to
    // the subread then aligned.
    //
    // subreadMapType of 0 is always used, however eventually it
    // may be faster to go to 1, just 1 isn't tested thoroughly
    // yet.
    //
    // Note, for proper SAM printing, subreadMaptype of 0 is needed.
    //
    int start = subreadInterval.start;
    int end   = subreadInterval.end;
        
    assert(smrtRead.length >= subreadSequence.length);
    if (params.subreadMapType == 0) {
      smrtRead.MakeSubreadAsMasked(subreadSequence, start, end); 
    }
    else if (params.subreadMapType == 1) {
      smrtRead.MakeSubreadAsReference(subreadSequence, start, end); 
    }

    if (!params.preserveReadTitle) {
      smrtRead.SetSubreadTitle(subreadSequence, 
                               subreadSequence.subreadStart,
                               subreadSequence.subreadEnd);
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
                   SMRTSequence & smrtRead) {
  assert(smrtRead.length >= subreadSequence.length);
  // Reverse complement sequence of the subread.
  subreadSequence.MakeRC(subreadSequenceRC);
  // Update start and end positions of subreadSequenceRC in the 
  // coordinate of reverse compelement sequence of the SMRT read.
  subreadSequenceRC.subreadStart = smrtRead.length - subreadSequence.subreadEnd;
  subreadSequenceRC.subreadEnd   = smrtRead.length - subreadSequence.subreadStart;
  subreadSequenceRC.zmwData      = smrtRead.zmwData;
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
void PrintAllReadAlignments(ReadAlignments & allReadAlignments,
                            AlignmentContext & alignmentContext,
                            ostream & outFilePtr,
                            ostream & unalignedFilePtr,
                            MappingParameters & params) {
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

  if (alignmentContext.numProperlyAlignedSubreads == allReadAlignments.subreadAlignments.size()) {
      alignmentContext.allSubreadsProperlyAligned = true;
  }
  alignmentContext.nSubreads = nAlignedSubreads;

  for (subreadIndex = 0; subreadIndex < nAlignedSubreads; subreadIndex++) {
    alignmentContext.subreadIndex = subreadIndex;
    if (subreadIndex < nAlignedSubreads-1 and allReadAlignments.subreadAlignments[subreadIndex+1].size() > 0) {
      alignmentContext.nextSubreadPos = allReadAlignments.subreadAlignments[subreadIndex+1][0]->QAlignStart();
      alignmentContext.nextSubreadDir = allReadAlignments.subreadAlignments[subreadIndex+1][0]->qStrand;
      alignmentContext.rNext = allReadAlignments.subreadAlignments[subreadIndex+1][0]->tName;
      alignmentContext.hasNextSubreadPos = true;
    } else {
      alignmentContext.nextSubreadPos = 0;
      alignmentContext.nextSubreadDir = 0;
      alignmentContext.rNext = "";
      alignmentContext.hasNextSubreadPos = false;
    }
      
    if (allReadAlignments.subreadAlignments[subreadIndex].size() > 0) {
        PrintAlignments(allReadAlignments.subreadAlignments[subreadIndex], 
                        allReadAlignments.subreads[subreadIndex], // the source read
                        // for these alignments
                        params, outFilePtr,//*mapData->outFilePtr,
                        alignmentContext);   
    } else {
      //
      // Print the unaligned sequences.
      //
      if (params.printUnaligned == true) {
        if (params.nProc == 1) {
          //allReadAlignments.subreads[subreadIndex].PrintSeq(*mapData->unalignedFilePtr);
          allReadAlignments.subreads[subreadIndex].PrintSeq(unalignedFilePtr);
        }
        else {
#ifdef __APPLE__
          sem_wait(semaphores.unaligned);
#else
          sem_wait(&semaphores.unaligned);
#endif
          //allReadAlignments.subreads[subreadIndex].PrintSeq(*mapData->unalignedFilePtr);
          allReadAlignments.subreads[subreadIndex].PrintSeq(unalignedFilePtr);
#ifdef __APPLE__
          sem_post(semaphores.unaligned);
#else
          sem_post(&semaphores.unaligned);
#endif
        } // End of nproc > 1.
      } // End of printing  unaligned sequences.
    } // End of finding no alignments for the subread with subreadIndex.
  } // End of printing and processing alignmentContext for each subread.
}

//template<typename T_SuffixArray, typename T_GenomeSequence, typename T_Tuple>
void MapReads(MappingData<T_SuffixArray, T_GenomeSequence, T_Tuple> *mapData) { 
  //
  // Step 1, initialize local pointers to map data 
  // for programming shorthand.
  //
  MappingParameters params = mapData->params;

  DNASuffixArray sarray;
  TupleCountTable<T_GenomeSequence, DNATuple> ct;
  SequenceIndexDatabase<FASTQSequence> seqdb;
  FASTQSequence fastaGenome;
  T_GenomeSequence    genome;
  BWT *bwtPtr;

  mapData->ShallowCopySuffixArray(sarray);
  mapData->ShallowCopyReferenceSequence(genome);
  mapData->ShallowCopySequenceIndexDatabase(seqdb);
  mapData->ShallowCopyTupleCountTable(ct);
  
  bwtPtr = mapData->bwtPtr;
  SeqBoundaryFtr<FASTQSequence> seqBoundary(&seqdb);

  VectorIndex i, j;
  if (params.match != 0) {
  }    

  int numAligned = 0;
  
  SMRTSequence smrtRead, smrtReadRC;
  SMRTSequence unrolledReadRC;
  CCSSequence  ccsRead;
  RegionAnnotation annotation;
  T_Sequence read;
  int readIndex = -1;
  int readsFileIndex;
  bool readIsCCS = false;


  // Print verbose logging to pid.threadid.log for each thread.
  ofstream threadOut;
  if (params.verbosity >= 3) {
    stringstream ss;
    ss << getpid() << "." << pthread_self();
    string threadLogFileName = ss.str() + ".log";
    threadOut.open(threadLogFileName.c_str(), ios::out|ios::app);
  }

  //
  // Reuse the following buffers during alignment.  Since these keep
  // storage contiguous, hopefully this will decrease memory
  // fragmentation.
  //
  MappingBuffers mappingBuffers;
  while (true) {
    //
    // Scan the next read from input.  This may either be a CCS read,
    // or regular read (though this may be aligned in whole, or by
    // subread).
    //
    AlignmentContext alignmentContext;
    // Associate each sequence to read in with a determined random int.
    int associatedRandInt = 0;
    if (mapData->reader->GetFileType() == HDFCCS ||
        mapData->reader->GetFileType() == HDFCCSONLY) { 
      if (GetNextReadThroughSemaphore(*mapData->reader, params, ccsRead, alignmentContext, associatedRandInt) == false) {
        break;
      }
      else {
        if (params.unrollCcs == false) {
          readIsCCS = true;
          smrtRead.Copy(ccsRead);
          ccsRead.SetQVScale(params.qvScaleType);
        }
        else {
          smrtRead.Copy(ccsRead.unrolledRead);
        }
        ++readIndex;
        smrtRead.SetQVScale(params.qvScaleType);
      }
      assert(ccsRead.zmwData.holeNumber == smrtRead.zmwData.holeNumber and
             ccsRead.zmwData.holeNumber == ccsRead.unrolledRead.zmwData.holeNumber);
    } else {
      if (GetNextReadThroughSemaphore(*mapData->reader, params, smrtRead, alignmentContext, associatedRandInt) == false) {
        break;
      }
      else {
        ++readIndex;
        smrtRead.SetQVScale(params.qvScaleType);
      }
    }

    //
    // Only normal (non-CCS) reads should be masked.  Since CCS reads store the raw read, that is masked.
    //
    bool readHasGoodRegion = true;
    if (params.useRegionTable and params.useHQRegionTable) {
      if (readIsCCS) {
        readHasGoodRegion = MaskRead(ccsRead.unrolledRead, ccsRead.unrolledRead.zmwData, *mapData->regionTablePtr);       
      }
      else {
        readHasGoodRegion = MaskRead(smrtRead, smrtRead.zmwData, *mapData->regionTablePtr);
      }
      //
      // Store the high quality start and end of this read for masking purposes when printing.
      //
      int hqStart, hqEnd;
      int score;
      LookupHQRegion(smrtRead.zmwData.holeNumber, *mapData->regionTablePtr, hqStart, hqEnd, score);
      smrtRead.lowQualityPrefix = hqStart;
      smrtRead.lowQualitySuffix = smrtRead.length - hqEnd;
    }
    else {
      smrtRead.lowQualityPrefix = 0;
      smrtRead.lowQualitySuffix = 0;
    }

    //
    // Give the opportunity to align a subset of reads.
    //
    if ((params.maxReadIndex >= 0 and 
         smrtRead.zmwData.holeNumber > params.maxReadIndex)) {
      if (readIsCCS) {
        ccsRead.Free();
      }
      smrtRead.Free();
      break;
    }

    if (readHasGoodRegion == false or 
        (params.readIndex >= 0 and 
         smrtRead.zmwData.holeNumber != params.readIndex) or
        (params.holeNumberRangesStr.size() > 0 and
         not params.holeNumberRanges.contains(smrtRead.zmwData.holeNumber))) {
      //
      // Nothing to do with this read. Skip aligning it entirely.
      //
      if (readIsCCS) {
        ccsRead.Free();
      }
      smrtRead.Free();
      // Stop processing once the specified read index is reached.
      // Eventually this will change to just seek to readIndex, and
      // just align one read anyway.
      if ((params.readIndex >= 0 and 
           smrtRead.zmwData.holeNumber > params.readIndex) or
          (params.holeNumberRangesStr.size() > 0 and 
           smrtRead.zmwData.holeNumber > params.holeNumberRanges.max())){
        break;
      }
      continue;
    }
    

    // 
    // Discard reads that are too small, or not labeled as having any
    // useable/good sequence. 
    //
    if (smrtRead.length < params.minReadLength or readHasGoodRegion == false or 
                          (params.maxReadLength != 0 and 
                           smrtRead.length > params.maxReadLength)) {
      if (readIsCCS) {
        ccsRead.Free();
      }
      smrtRead.Free();
      continue;
    }

    if (smrtRead.qual.Empty() == false and smrtRead.GetAverageQuality() < params.minAvgQual) {
      if (readIsCCS) {
        ccsRead.Free();
      }
      smrtRead.Free();
      continue;
    }

    if (params.verbosity > 1) {
      cout << "aligning read: " << endl;
      smrtRead.PrintSeq(cout);
    }

    smrtRead.MakeRC(smrtReadRC);
    
    if (readIsCCS) {
      ccsRead.unrolledRead.MakeRC(unrolledReadRC);
    }

    //
    // When aligning subreads separately, iterate over each subread, and print the alignments for these.
    //

    ReadAlignments allReadAlignments;
    allReadAlignments.read = smrtRead;

    if (readIsCCS == false and params.mapSubreadsSeparately) {
      // (not readIsCCS and not -noSplitSubreads)
      vector<ReadInterval> subreadIntervals;
      vector<int>          subreadDirections;
      vector<ReadInterval> adapterIntervals; 
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
      }
      else {
        // 
        // Grab the subread & adapter intervals from the entire region table to
        // iterate over.
        //
        CollectSubreadIntervals(smrtRead, mapData->regionTablePtr, subreadIntervals, params.byAdapter);
        CollectAdapterIntervals(smrtRead, mapData->regionTablePtr, adapterIntervals);
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
      int longestSubreadIndex = GetHighQualitySubreadsIntervals(
            subreadIntervals, // a vector of subread intervals.
            subreadDirections, // a vector of subread directions.
            smrtRead.lowQualityPrefix, // hq region start pos.
            smrtRead.length - smrtRead.lowQualitySuffix, // hq end pos.
            params.minSubreadLength); // minimum read length.
      
      // Get index of the (left-most) longest full-pass subread which has
      // both adapters.
      int longestFullSubreadIndex = GetLongestFullSubreadIndex(
          subreadIntervals, adapterIntervals);

      int bestSubreadIndex = (longestFullSubreadIndex >= 0)?
          longestFullSubreadIndex:longestSubreadIndex;

      // Flop all directions if direction of the longest subread is 1.
      if (bestSubreadIndex >= 0 and 
          bestSubreadIndex < subreadDirections.size() and
          subreadDirections[bestSubreadIndex] == 1) {
        UpdateDirections(subreadDirections, true);
      }

      int startIndex = 0;
      int endIndex = subreadIntervals.size();

      if (params.concordant) {
        // Only the longest subread will be aligned in the first round.
        startIndex = max(startIndex, bestSubreadIndex);
        endIndex   = min(endIndex, bestSubreadIndex + 1);
      }

      //
      // Make room for alignments.
      //
      allReadAlignments.Resize(subreadIntervals.size());
      allReadAlignments.alignMode = Subread;

      DNALength intvIndex;
      for (intvIndex = startIndex; intvIndex < endIndex; intvIndex++) {
        SMRTSequence subreadSequence, subreadSequenceRC;
        MakeSubreadOfInterval(subreadSequence, smrtRead, 
                              subreadIntervals[intvIndex], params);
        MakeSubreadRC(subreadSequenceRC, subreadSequence, smrtRead);


        //
        // Store the sequence that is being mapped in case no hits are
        // found, and missing sequences are printed.
        //
        allReadAlignments.SetSequence(intvIndex, subreadSequence);

        vector<T_AlignmentCandidate*> alignmentPtrs;
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
                ct,               // Count table to use word frequencies in the genome to weight matches.
                seqdb,            // Information about the names of
                                  // chromosomes in the genome, and
                                  // where their sequences are in the genome.
                params,           // A huge list of parameters for
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
                mapData           // Some values that are shared
                                  // across threads.
                );

        //
        // No alignments were found, sometimes parameters are
        // specified to try really hard again to find an alignment.
        // This sets some parameters that use a more sensitive search
        // at the cost of time.
        //

        if ((alignmentPtrs.size() == 0 or alignmentPtrs[0]->pctSimilarity < 80) and params.doSensitiveSearch) {
          MappingParameters sensitiveParams = params;
          sensitiveParams.SetForSensitivity();
          MapRead(subreadSequence, subreadSequenceRC, genome, 
                  sarray, *bwtPtr, 
                  seqBoundary, ct, seqdb,
                  sensitiveParams, mapData->metrics, 
                  alignmentPtrs, mappingBuffers, 
                  mapData);
        }

        //
        // Store the mapping quality values.
        //
        if (alignmentPtrs.size() > 0 and 
            alignmentPtrs[0]->score < params.maxScore and 
            params.storeMapQV) {
          StoreMapQVs(subreadSequence, alignmentPtrs, params);
        }

        // 
        // Select alignments for this subread.
        //
        vector<T_AlignmentCandidate*> selectedAlignmentPtrs =
        SelectAlignmentsToPrint(alignmentPtrs, params, 
                                associatedRandInt);
        allReadAlignments.AddAlignmentsForSeq(intvIndex, selectedAlignmentPtrs);

        //
        // Move reference from subreadSequence, which will be freed at
        // the end of this loop to the smrtRead, which exists for the
        // duration of aligning all subread of the smrtRead.
        //
        if (params.subreadMapType == 0) {
          int a;
          for (a = 0; a < alignmentPtrs.size(); a++) {
            if (alignmentPtrs[a]->qStrand == 0) {
              alignmentPtrs[a]->qAlignedSeq.ReferenceSubstring(smrtRead,
                                                               alignmentPtrs[a]->qAlignedSeq.seq - subreadSequence.seq,
                                                               alignmentPtrs[a]->qAlignedSeqLength);
            }
            else {
              alignmentPtrs[a]->qAlignedSeq.ReferenceSubstring(smrtReadRC,
                                                               alignmentPtrs[a]->qAlignedSeq.seq - subreadSequenceRC.seq, 
                                                               alignmentPtrs[a]->qAlignedSeqLength);
            }
          }
          subreadSequence.Free();
          subreadSequenceRC.Free();
        }
      } // End of looping over subread intervals within [startIndex, endIndex).

      if (params.verbosity >= 3) 
          allReadAlignments.Print(threadOut);

      if (params.concordant) {
        allReadAlignments.read = smrtRead;
        allReadAlignments.alignMode = ZmwSubreads;

        if (startIndex >= 0 && 
            startIndex < allReadAlignments.subreadAlignments.size()) {
          vector<T_AlignmentCandidate*> selectedAlignmentPtrs =
              allReadAlignments.CopySubreadAlignments(startIndex);

          for(int alignmentIndex = 0; alignmentIndex < selectedAlignmentPtrs.size(); 
              alignmentIndex++) {
            FlankTAlignedSeq(selectedAlignmentPtrs[alignmentIndex], 
                              seqdb, genome, params.flankSize);
          }

          for (intvIndex = 0; intvIndex < subreadIntervals.size(); intvIndex++) {
            if (intvIndex == startIndex) continue; 
            int passDirection = subreadDirections[intvIndex];
            int passStartBase = subreadIntervals[intvIndex].start;
            int passNumBases  = subreadIntervals[intvIndex].end - passStartBase;
            if (passNumBases <= params.minReadLength) {continue;}

            mapData->metrics.numReads++;
            SMRTSequence subread;
            subread.ReferenceSubstring(smrtRead, passStartBase, passNumBases);
            subread.CopyTitle(smrtRead.title);
            // The unrolled alignment should be relative to the entire read.
            allReadAlignments.SetSequence(intvIndex, smrtRead);

            for (int alnIndex = 0; alnIndex < selectedAlignmentPtrs.size(); alnIndex++) {
              T_AlignmentCandidate * alignment = selectedAlignmentPtrs[alnIndex];
              if (alignment->score > params.maxScore) break;
              AlignSubreadToAlignmentTarget(allReadAlignments, 
                                            subread,
                                            smrtRead,
                                            alignment,
                                            passDirection, 
                                            subreadIntervals[intvIndex],
                                            intvIndex,
                                            params, mappingBuffers, threadOut);
            } // End of aligning this subread to each selected alignment.
          } // End of aligning each subread to where the longest subread aligned to. 
          for(int alignmentIndex = 0; alignmentIndex < selectedAlignmentPtrs.size(); 
              alignmentIndex++) {
            if (selectedAlignmentPtrs[alignmentIndex]) 
              delete selectedAlignmentPtrs[alignmentIndex];
          }
        } // End of if startIndex >= 0 and < subreadAlignments.size()
      } // End of if params.concordant
    } // End of if (readIsCCS == false and params.mapSubreadsSeparately).
    else { // if (readIsCCS or (not readIsCCS and -noSplitSubreads) )
      //
      // The read must be mapped as a whole, even if it contains subreads.
      //
      vector<T_AlignmentCandidate*> alignmentPtrs;
      mapData->metrics.numReads++;
      smrtRead.subreadStart   = 0; smrtRead.subreadEnd   = smrtRead.length;
      smrtReadRC.subreadStart = 0; smrtReadRC.subreadEnd = smrtRead.length;

      MapRead(smrtRead, smrtReadRC, 
              genome, sarray, *bwtPtr, seqBoundary, ct, seqdb, params, mapData->metrics,
              alignmentPtrs, mappingBuffers, mapData);
      
      //
      // Store the mapping quality values.
      //
      if (alignmentPtrs.size() > 0 and 
          alignmentPtrs[0]->score < params.maxScore and 
          params.storeMapQV) {
        StoreMapQVs(smrtRead, alignmentPtrs, params);
      }

      // 
      // Select de novo ccs-reference alignments for subreads to align to.
      //
      vector<T_AlignmentCandidate*> selectedAlignmentPtrs =
      SelectAlignmentsToPrint(alignmentPtrs, params, associatedRandInt);

      //
      // Just one sequence is aligned.  There is one primary hit, and
      // all other are secondary.
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
        }
        else {
          allReadAlignments.alignMode = Fullread;
        }
        allReadAlignments.SetSequence(0, smrtRead);
      }
      else if (readIsCCS) { // if -useccsall or -useccs 
        // Flank alignment candidates to both ends.
        for(int alignmentIndex = 0; alignmentIndex < selectedAlignmentPtrs.size(); 
          alignmentIndex++) {
          FlankTAlignedSeq(selectedAlignmentPtrs[alignmentIndex],
                           seqdb, genome, params.flankSize);
        }

        //
        // Align the ccs subread to where the denovo sequence mapped (explode).
        //
        SMRTSequence readRC;
        
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
        }
        else {
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
          if (passNumBases <= params.minReadLength) { continue; }

          ReadInterval subreadInterval(passStartBase, passStartBase + passNumBases);

          subread.ReferenceSubstring(ccsRead.unrolledRead, passStartBase, passNumBases-1);
          subread.CopyTitle(ccsRead.title);
          // The unrolled alignment should be relative to the entire read.
          allReadAlignments.SetSequence(subreadIndex, ccsRead.unrolledRead);

          int alignmentIndex;
          //
          // Align this subread to all the positions that the de novo
          // sequence has aligned to.
          //
          for (alignmentIndex = 0; alignmentIndex < selectedAlignmentPtrs.size(); alignmentIndex++) {
            T_AlignmentCandidate *alignment = selectedAlignmentPtrs[alignmentIndex];
            if (alignment->score > params.maxScore) break;
            AlignSubreadToAlignmentTarget(allReadAlignments,
                                          subread, ccsRead.unrolledRead, 
                                          alignment,
                                          passDirection,
                                          subreadInterval,
                                          subreadIndex,
                                          params, mappingBuffers, threadOut);
          } // End of aligning this subread to where the de novo ccs has aligned to.
        } // End of alignining all subreads to where the de novo ccs has aligned to.
      } // End of if readIsCCS and !params.useCcsOnly 

    } // End of if not (readIsCCS == false and params.mapSubreadsSeparately) 

    PrintAllReadAlignments(allReadAlignments, alignmentContext,
                           *mapData->outFilePtr, 
                           *mapData->unalignedFilePtr,
                           params);
        
    allReadAlignments.Clear();
    smrtReadRC.Free();
    smrtRead.Free();

    if (readIsCCS) {
      ccsRead.Free();
      unrolledReadRC.Free();
    }
    numAligned++;
    if(numAligned % 100 == 0) {
      mappingBuffers.Reset();
    }
  } // End of while (true).
  if (params.nProc > 1) {
#ifdef __APPLE__
    sem_wait(semaphores.reader);
    sem_post(semaphores.reader);
#else
    sem_wait(&semaphores.reader);
    sem_post(&semaphores.reader);
#endif
  }
  if (params.nProc > 1) {
    pthread_exit(NULL); 
  }
  threadOut.close();
}

float ComputePMatch(float accuracy, int anchorLength) {
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
                                       float &mean, float &variance) {
  mean = 1/pSuccess;
  variance = (1-pSuccess)/ (pow(pSuccess,2));
}

int ComputeExpectedWaitingBases(float mean, float variance, float certainty) {
  float nStdDev;
  assert(FindQNorm(certainty, nStdDev) != false);
  return mean + sqrt(variance) * nStdDev;
}

int main(int argc, char* argv[]) {
  //
  // Configure parameters for refining alignments.
  //
  MappingParameters params;
  ReverseCompressIndex index;
  pid_t parentPID;
  pid_t *pids;
  
  CommandLineParser clp;
  string commandLine;
  string helpString;
  SetHelp(helpString);

  string conciseHelpString;
  SetConciseHelp(conciseHelpString);
  
  stringstream usageSStrm;
  usageSStrm << "   Basic usage: 'blasr reads.{fasta,bax.h5} genome.fasta [-options] " << endl
             << " [option]\tDescription (default_value)." << endl << endl
             << " Input Files." << endl
             << "   reads.fasta is a multi-fasta file of reads.  While any fasta file is valid input, " 
             "it is preferable to use plx.h5 or bax.h5 files because they contain "
             "more rich quality value information." << endl
             << "   reads.bax.h5|reads.plx.h5 Is the native output format in Hierarchical Data Format of "
    "SMRT reads. This is the preferred input to blasr because rich quality"
    "value (insertion,deletion, and substitution quality values) information is "
    "maintained.  The extra quality information improves variant detection and mapping"<<
    "speed." << endl << endl;


  clp.SetHelp(helpString);
  clp.SetConciseHelp(conciseHelpString);
  clp.SetProgramSummary(usageSStrm.str());
  clp.SetProgramName("blasr");

  string versionString;
  GetVersion(versionString);
  clp.SetVersion(versionString);

  //
  // Make the default arguments.
  //
  bool required=true;
  bool optional=false;
  int  trashbinInt;
  float trashbinFloat;
  bool trashbinBool;
  //  clp.RegisterStringListOption("input_files", &params.readsFileNames, "Read files followed by genome.");
  //  clp.RegisterPreviousFlagsAsHidden();
  bool printVerboseHelp = false;
  bool printLongHelp    = false;
  clp.RegisterStringOption("sa", &params.suffixArrayFileName, "");
  clp.RegisterStringOption("ctab", &params.countTableName, "" );
  clp.RegisterStringOption("regionTable", &params.regionTableFileName, "");
  clp.RegisterStringOption("ccsFofn", &params.ccsFofnFileName, "");
  clp.RegisterIntOption("bestn", (int*) &params.nBest, "", CommandLineParser::PositiveInteger);
  clp.RegisterIntOption("limsAlign", &params.limsAlign, "", CommandLineParser::PositiveInteger);
  clp.RegisterFlagOption("printOnlyBest", &params.printOnlyBest, "");
  clp.RegisterFlagOption("outputByThread", &params.outputByThread, "");
  clp.RegisterFlagOption("rbao", &params.refineBetweenAnchorsOnly, "");
  clp.RegisterFlagOption("allowAdjacentIndels", &params.forPicard, "");
  clp.RegisterFlagOption("onegap", &params.separateGaps, "");
  clp.RegisterFlagOption("allowAdjacentIndels", &params.forPicard, "");
  clp.RegisterFlagOption("placeRepeatsRandomly", &params.placeRandomly, "");
  clp.RegisterIntOption("randomSeed", &params.randomSeed, "", CommandLineParser::Integer);
  clp.RegisterFlagOption("extend", &params.extendAlignments, "");
  clp.RegisterIntOption("branchExpand", &params.anchorParameters.branchExpand, "", CommandLineParser::NonNegativeInteger);
  clp.RegisterIntOption("maxExtendDropoff", &params.maxExtendDropoff, "", CommandLineParser::NonNegativeInteger);
  clp.RegisterFlagOption("nucmer", &params.emulateNucmer, "");
  clp.RegisterIntOption("maxExpand", &params.maxExpand, "", CommandLineParser::PositiveInteger);
  clp.RegisterIntOption("minExpand", &params.minExpand, "", CommandLineParser::NonNegativeInteger);
  clp.RegisterStringOption("seqdb",  &params.seqDBName, "");
  clp.RegisterStringOption("anchors",  &params.anchorFileName, "");
  clp.RegisterStringOption("clusters", &params.clusterFileName, "");
  clp.RegisterFlagOption("samplePaths", (bool*) &params.samplePaths, "");
  clp.RegisterFlagOption("noStoreMapQV", &params.storeMapQV, "");
  clp.RegisterFlagOption("nowarp", (bool*) &params.nowarp, "");
  clp.RegisterFlagOption("noRefineAlign", (bool*) &params.refineAlign, "");
  clp.RegisterFlagOption("guidedAlign", (bool*)&params.useGuidedAlign, "");
  clp.RegisterFlagOption("useGuidedAlign", (bool*)&trashbinBool, "");
  clp.RegisterFlagOption("noUseGuidedAlign", (bool*)&params.useGuidedAlign, "");
  clp.RegisterFlagOption("header", (bool*)&params.printHeader, "");
  clp.RegisterIntOption("subreadImplType", &params.subreadMapType, "", CommandLineParser::PositiveInteger);
  clp.RegisterIntOption("bandSize", &params.bandSize, "", CommandLineParser::PositiveInteger);  
  clp.RegisterIntOption("extendBandSize", &params.extendBandSize, "", CommandLineParser::PositiveInteger);  
  clp.RegisterIntOption("guidedAlignBandSize", &params.guidedAlignBandSize, "", CommandLineParser::PositiveInteger);  
  clp.RegisterIntOption("maxAnchorsPerPosition", &params.anchorParameters.maxAnchorsPerPosition, "", CommandLineParser::PositiveInteger);
  clp.RegisterIntOption("stopMappingOnceUnique", (int*) &params.anchorParameters.stopMappingOnceUnique, "", CommandLineParser::NonNegativeInteger);
  clp.RegisterStringOption("out", &params.outFileName, "");
  clp.RegisterIntOption("match", &params.match, "", CommandLineParser::Integer);
  clp.RegisterIntOption("mismatch", &params.mismatch, "", CommandLineParser::Integer);
  clp.RegisterIntOption("minMatch", &params.minMatchLength, "", CommandLineParser::PositiveInteger);
  clp.RegisterIntOption("maxMatch", &params.anchorParameters.maxLCPLength, "", CommandLineParser::NonNegativeInteger);
  clp.RegisterIntOption("maxLCPLength", &params.anchorParameters.maxLCPLength, "", CommandLineParser::NonNegativeInteger);
  clp.RegisterIntOption("indel", &params.indel, "", CommandLineParser::Integer);
  clp.RegisterIntOption("insertion", &params.insertion, "", CommandLineParser::Integer);
  clp.RegisterIntOption("deletion", &params.deletion, "", CommandLineParser::Integer);
  clp.RegisterIntOption("idsIndel", &params.idsIndel, "", CommandLineParser::Integer);
  clp.RegisterIntOption("sdpindel", &params.sdpIndel, "", CommandLineParser::Integer);
  clp.RegisterIntOption("sdpIns", &params.sdpIns, "", CommandLineParser::Integer);
  clp.RegisterIntOption("sdpDel", &params.sdpDel, "", CommandLineParser::Integer);
  clp.RegisterFloatOption("indelRate", &params.indelRate, "", CommandLineParser::NonNegativeFloat);
  clp.RegisterFloatOption("minRatio", &params.minRatio, "", CommandLineParser::NonNegativeFloat); 
  clp.RegisterFloatOption("sdpbypass", &params.sdpBypassThreshold, "", CommandLineParser::NonNegativeFloat);
  clp.RegisterFloatOption("minFrac", &trashbinFloat, "", CommandLineParser::NonNegativeFloat);
  clp.RegisterIntOption("maxScore", &params.maxScore, "", CommandLineParser::Integer);
  clp.RegisterStringOption("bwt", &params.bwtFileName, "");
  clp.RegisterIntOption("m", &params.printFormat, "", CommandLineParser::NonNegativeInteger);
  clp.RegisterFlagOption("sam", &params.printSAM, "");
  clp.RegisterStringOption("clipping", &params.clippingString, "");
  clp.RegisterIntOption("sdpTupleSize", &params.sdpTupleSize, "", CommandLineParser::PositiveInteger);
  clp.RegisterIntOption("pvaltype", &params.pValueType, "", CommandLineParser::NonNegativeInteger);
  clp.RegisterIntOption("start", &params.startRead, "", CommandLineParser::NonNegativeInteger);
  clp.RegisterIntOption("stride", &params.stride, "", CommandLineParser::NonNegativeInteger);
  clp.RegisterFloatOption("subsample", &params.subsample, "", CommandLineParser::PositiveFloat);
  clp.RegisterIntOption("nproc", &params.nProc, "", CommandLineParser::PositiveInteger);
  clp.RegisterFlagOption("sortRefinedAlignments",(bool*) &params.sortRefinedAlignments, "");
  clp.RegisterIntOption("quallc", &params.qualityLowerCaseThreshold, "", CommandLineParser::Integer);
  clp.RegisterFlagOption("v", (bool*) &params.verbosity, "");
  clp.RegisterIntOption("V", &params.verbosity, "Specify a level of verbosity.", CommandLineParser::NonNegativeInteger);
  clp.RegisterIntOption("contextAlignLength", &params.anchorParameters.contextAlignLength, "", CommandLineParser::PositiveInteger);
  clp.RegisterFlagOption("skipLookupTable", &params.anchorParameters.useLookupTable, "");
  clp.RegisterStringOption("metrics", &params.metricsFileName, "");
  clp.RegisterStringOption("lcpBounds", &params.lcpBoundsFileName, "");
  clp.RegisterStringOption("fullMetrics", &params.fullMetricsFileName, "");
  clp.RegisterIntOption("nbranch", &params.anchorParameters.numBranches, "", CommandLineParser::NonNegativeInteger);
  clp.RegisterFlagOption("divideByAdapter", &params.byAdapter, "");
  clp.RegisterFlagOption("useQuality", &params.ignoreQualities, "");
  clp.RegisterFlagOption("noFrontAlign", &params.extendFrontAlignment, "");
  clp.RegisterIntOption("minReadLength", &params.minReadLength, "", CommandLineParser::NonNegativeInteger);
  clp.RegisterIntOption("maxReadLength", &params.maxReadLength, "", CommandLineParser::NonNegativeInteger);
  clp.RegisterIntOption("minSubreadLength", &params.minSubreadLength, "", CommandLineParser::NonNegativeInteger);
  clp.RegisterIntOption("minAvgQual", &params.minAvgQual, "", CommandLineParser::Integer);
  clp.RegisterFlagOption("advanceHalf", &params.advanceHalf, "");
  clp.RegisterIntOption("advanceExactMatches", &params.anchorParameters.advanceExactMatches, "", CommandLineParser::NonNegativeInteger);
  clp.RegisterFlagOption("unrollCcs", &params.unrollCcs, "");
  clp.RegisterFlagOption("useccs", &params.useCcs, "");
  clp.RegisterFlagOption("useccsdenovo", &params.useCcsOnly, "");
  clp.RegisterFlagOption("useccsall", &params.useAllSubreadsInCcs, "");
  clp.RegisterFlagOption("extendDenovoCCSSubreads", &params.extendDenovoCCSSubreads, "");
  clp.RegisterFlagOption("noRefineAlignments", &params.refineAlignments, "");
  clp.RegisterFloatOption("minPctIdentity", &params.minPctIdentity, "", CommandLineParser::NonNegativeFloat);
  clp.RegisterFloatOption("maxPctIdentity", &params.maxPctIdentity, "", CommandLineParser::NonNegativeFloat);
  clp.RegisterIntOption("nCandidates", &params.nCandidates, "", CommandLineParser::NonNegativeInteger);
  clp.RegisterFlagOption("useTemp", (bool*) &params.tempDirectory, "");
  clp.RegisterFlagOption("noSplitSubreads", &params.mapSubreadsSeparately, "");
  clp.RegisterFlagOption("concordant", &params.concordant, "");
  clp.RegisterIntOption("flankSize", &params.flankSize, "", CommandLineParser::NonNegativeInteger);
  clp.RegisterIntOption("subreadMapType", &params.subreadMapType, "", CommandLineParser::NonNegativeInteger);
  clp.RegisterStringOption("titleTable", &params.titleTableName, "");
  clp.RegisterFlagOption("useSensitiveSearch", &params.doSensitiveSearch, "");
  clp.RegisterFlagOption("ignoreRegions", &params.useRegionTable, "");
  clp.RegisterFlagOption("ignoreHQRegions", &params.useHQRegionTable, "");
  clp.RegisterFlagOption("computeAlignProbability", &params.computeAlignProbability, "");
  clp.RegisterStringOption("unaligned", &params.unalignedFileName, "");
  clp.RegisterFlagOption("global", &params.doGlobalAlignment, "");
  clp.RegisterIntOption("globalChainType", &params.globalChainType, "", CommandLineParser::NonNegativeInteger);
  clp.RegisterFlagOption("noPrintSubreadTitle", (bool*) &params.printSubreadTitle, "");
  clp.RegisterIntOption("saLookupTableLength", &params.lookupTableLength, "", CommandLineParser::PositiveInteger);
  clp.RegisterFlagOption("useDetailedSDP", &params.detailedSDPAlignment, "");
  clp.RegisterFlagOption("nouseDetailedSDP", &trashbinBool, "");
  clp.RegisterIntOption("sdpFilterType", &params.sdpFilterType, "", CommandLineParser::NonNegativeInteger);
  clp.RegisterIntOption("scoreType", &params.scoreType, "", CommandLineParser::NonNegativeInteger);
  clp.RegisterFlagOption("help", &params.printDiscussion, "");
  clp.RegisterFlagOption("h", &printVerboseHelp, "");
  clp.RegisterFloatOption("accuracyPrior",    &params.readAccuracyPrior, "", CommandLineParser::NonNegativeFloat);
  clp.RegisterIntOption("readIndex", &params.readIndex, "", CommandLineParser::NonNegativeInteger);
  // holeNumberRangesStr is a string of comma-delimited hole number ranges, such as '1,2,3,10-15'.
  // Blasr only analyzes reads whose hole numbers are in the specified hole number ranges. 
  clp.RegisterStringOption("holeNumbers", &params.holeNumberRangesStr, "");
  clp.RegisterIntOption("maxReadIndex", &params.maxReadIndex, "", CommandLineParser::NonNegativeInteger);
  clp.RegisterIntOption("substitutionPrior",  &params.substitutionPrior, "", CommandLineParser::NonNegativeInteger);
  clp.RegisterIntOption("deletionPrior",  &params.globalDeletionPrior, "", CommandLineParser::NonNegativeInteger);
  clp.RegisterIntOption("recurseOver", &params.recurseOver, "", CommandLineParser::NonNegativeInteger);
  clp.RegisterStringOption("scoreMatrix", &params.scoreMatrixString, "");
  clp.RegisterFlagOption("printDotPlots", &params.printDotPlots, "");
  clp.RegisterFlagOption("preserveReadTitle", &params.preserveReadTitle,"");
  clp.RegisterFlagOption("forwardOnly", &params.forwardOnly,"");
  clp.RegisterFlagOption("affineAlign", &params.affineAlign, "");
  clp.RegisterIntOption("affineOpen", &params.affineOpen, "", CommandLineParser::NonNegativeInteger);
  clp.RegisterIntOption("affineExtend", &params.affineExtend, "", CommandLineParser::NonNegativeInteger);
  clp.RegisterFlagOption("scaleMapQVByNClusters", &params.scaleMapQVByNumSignificantClusters, "", false);
  clp.ParseCommandLine(argc, argv, params.readsFileNames);
  clp.CommandLineToString(argc, argv, commandLine);

  if (printVerboseHelp) {
    cout << helpString << endl;
    exit(0);
  }

  if (params.printDiscussion) {
    PrintDiscussion();
    exit(0);
  }
  if (argc < 3) {
    cout << conciseHelpString;
    exit(1);
  }
  
  int a, b;
  for (a = 0; a < 5; a++ ){
    for (b = 0; b < 5; b++ ){
      if (a != b) {
        SMRTDistanceMatrix[a][b] += params.mismatch;
      }
      else {
        SMRTDistanceMatrix[a][b] += params.match;
      }
    }
  }
  
  if (params.scoreMatrixString != "") {
    if (StringToScoreMatrix(params.scoreMatrixString, SMRTDistanceMatrix) == false) {
      cout << "ERROR. The string " << endl
           << params.scoreMatrixString << endl
           << "is not a valid format.  It should be a quoted, space separated string of " << endl
           << "integer values.  The matrix: " << endl
           << "    A  C  G  T  N" << endl
           << " A  1  2  3  4  5" << endl
           << " C  6  7  8  9 10" << endl
           << " G 11 12 13 14 15" << endl
           << " T 16 17 18 19 20" << endl
           << " N 21 22 23 24 25" << endl
           << " should be specified as \"1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25\"" << endl;
      exit(1);
    }
  }
  
  FileOfFileNames::ExpandFileNameList(params.readsFileNames);

  // The last reads file is the genome
  if (params.readsFileNames.size() > 0) {
    params.genomeFileName = params.readsFileNames[params.readsFileNames.size()-1];
  }

  if (params.printDiscussion == true) {
    PrintDiscussion();
    exit(0);
  }

  cerr << "[INFO] " << GetTimestamp() << " [blasr] started." << endl;
  params.MakeSane();
  params.anchorParameters.verbosity = params.verbosity; 

  //
  // The random number generator is used for subsampling for debugging
  // and testing consensus.
  //
  if (params.useRandomSeed == true) {
    InitializeRandomGenerator(params.randomSeed);
  }
  else {
    InitializeRandomGeneratorWithTime();
  }
  
  //
  // Various aspects of timing are stored here.  However this isn't
  // quite finished.
  //
  MappingMetrics metrics;

  ofstream fullMetricsFile;
  if (params.fullMetricsFileName != "") {
    CrucialOpen(params.fullMetricsFileName, fullMetricsFile, std::ios::out);
    metrics.SetStoreList();
  }

  /*
   * If reading a separate region table, there is a 1-1 correspondence
   * between region table and bas file.
   */
  if (params.readSeparateRegionTable) {
    if (FileOfFileNames::IsFOFN(params.regionTableFileName)) {
      FileOfFileNames::FOFNToList(params.regionTableFileName, params.regionTableFileNames);
    }
    else {
      params.regionTableFileNames.push_back(params.regionTableFileName);
    }
  }

  if (params.regionTableFileNames.size() != 0 and 
      params.regionTableFileNames.size() != params.readsFileNames.size() - 1) {
    cout << "Error, there are not the same number of region table files as input files." << endl;
    exit(1);
  }

  // If reading a separate ccs fofn, there is a 1-1 corresponence 
  // between ccs fofn and base file.
  if (params.readSeparateCcsFofn) {
    if (FileOfFileNames::IsFOFN(params.ccsFofnFileName)) {
      FileOfFileNames::FOFNToList(params.ccsFofnFileName, params.ccsFofnFileNames);
    }
    else {
      params.ccsFofnFileNames.push_back(params.ccsFofnFileName);
    }
  }
  if (params.ccsFofnFileNames.size() != 0 and 
      params.ccsFofnFileNames.size() != params.readsFileNames.size() - 1) {
    cout << "Error, there are not the same number of ccs files as input files." << endl;
    exit(1);
  }
    
  if (params.readsFileNames.size() < 2) {
    cout << "Error, you must provide at least one reads file and a genome file." <<endl;
    exit(1);
  }
  
  //  The input reads files must have file extensions.
  for (int i = 0; i < params.readsFileNames.size()-1; i++) {
    size_t dotPos = params.readsFileNames[i].find_last_of('.');
    if (dotPos == string::npos) {
      cout<<"ERROR, the input reads files must include file extensions."<<endl;
      exit(1);
    }
  }

  // -useQuality can not be used in combination with a fasta input
  if (!params.ignoreQualities) {
    for (int i = 0; i < params.readsFileNames.size()-1; i++) {
      size_t dotPos = params.readsFileNames[i].find_last_of('.');
      assert (dotPos != string::npos); //dotPos should have been checked above
      string suffix = params.readsFileNames[i].substr(dotPos+1);
      if (suffix == "fasta") {
        cout<<"ERROR, you can not use -useQuality option when any of the input reads files are in multi-fasta format."<<endl; 
        exit(1);
      }
    }
  }

  parentPID = getpid();

  SequenceIndexDatabase<FASTASequence> seqdb;
  SeqBoundaryFtr<FASTASequence> seqBoundary(&seqdb);

  //
  // Initialize the sequence index database if it used. If it is not
  // specified, it is initialized by default when reading a multiFASTA
  // file.
  //
  if (params.useSeqDB) {
    ifstream seqdbin;
    CrucialOpen(params.seqDBName, seqdbin);
    seqdb.ReadDatabase(seqdbin);
  }

  //
  // Make sure the reads file exists and can be opened before
  // trying to read any of the larger data structures.
  //
  

  FASTASequence   fastaGenome;
  T_Sequence      genome;
  FASTAReader     genomeReader;

  // 
  // The genome is in normal FASTA, or condensed (lossy homopolymer->unipolymer) 
  // format.  Both may be read in using a FASTA reader.
  //
  if (!genomeReader.Init(params.genomeFileName)) {
    cout << "Could not open genome file " << params.genomeFileName << endl;
    exit(1);
  }

  if (params.printSAM) {
    genomeReader.computeMD5 = true;
  }
  //
  // If no sequence title database is supplied, initialize one when
  // reading in the reference, and consider a seqdb to be present.
  //
  if (!params.useSeqDB) {
    genomeReader.ReadAllSequencesIntoOne(fastaGenome, &seqdb);
    params.useSeqDB = true;
  }
  else {
    genomeReader.ReadAllSequencesIntoOne(fastaGenome);
  }
  genomeReader.Close();
  //
  // The genome may have extra spaces in the fasta name. Get rid of those.
  //
  VectorIndex t;
  for (t = 0; t < fastaGenome.titleLength; t++ ){
    if (fastaGenome.title[t] == ' ') {
      fastaGenome.titleLength = t;
      fastaGenome.title[t] = '\0';
      break;
    }
  }
  genome.seq = fastaGenome.seq;
  genome.length = fastaGenome.length;
  genome.title = fastaGenome.title;
  genome.titleLength = fastaGenome.titleLength;
  genome.ToUpper();


  DNASuffixArray sarray;
  TupleCountTable<T_GenomeSequence, DNATuple> ct;

  int listTupleSize;
  
  ofstream outFile;
  outFile.exceptions(ostream::failbit);
  ofstream unalignedOutFile;
  BWT bwt;
  
  if (params.useBwt) {
    if (bwt.Read(params.bwtFileName) == 0) {
      cout << "ERROR! Could not read the BWT file. " << params.bwtFileName << endl;
      exit(1);
    }
  }
  else {
    if (!params.useSuffixArray) {
      //
      // There was no explicit specification of a suffix
      // array on the command line, so build it on the fly here.
      //
      genome.ToThreeBit();    
      vector<int> alphabet;
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
    }
    else if (params.useSuffixArray) {
      if (sarray.Read(params.suffixArrayFileName)) {
        if (params.minMatchLength != 0) {
          params.listTupleSize = min(8, params.minMatchLength);
        }
        else {
          params.listTupleSize = sarray.lookupPrefixLength;
        }
        if (params.minMatchLength < sarray.lookupPrefixLength) {
          cerr << "WARNING. The value of -minMatch " << params.minMatchLength << " is less than the smallest searched length of " << sarray.lookupPrefixLength << ".  Setting -minMatch to " << sarray.lookupPrefixLength << "." << endl;
          params.minMatchLength = sarray.lookupPrefixLength;
        }
      }
      else {
        cout << "ERROR. " << params.suffixArrayFileName << " is not a valid suffix array. " << endl
             << " Make sure it is generated with the latest version of sawriter." << endl;
        exit(1);
      }
    }
  }
  
  if (params.minMatchLength < sarray.lookupPrefixLength) {
    cerr << "WARNING. The value of -minMatch " << params.minMatchLength << " is less than the smallest searched length of " << sarray.lookupPrefixLength << ".  Setting -minMatch to " << sarray.lookupPrefixLength << "." << endl;
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
  long l;
  TupleMetrics saLookupTupleMetrics;
  if (params.useCountTable) {
    ifstream ctIn;
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
    ofstream titleTableOut;
    CrucialOpen(params.titleTableName, titleTableOut);
    //
    // When using a sequence index database, the title table is simply copied 
    // from the sequencedb. 
    //
    if (params.useSeqDB) {
      titleTable.Copy(seqdb.names, seqdb.nSeqPos-1);
      titleTable.ResetTableToIntegers(seqdb.names, seqdb.nameLengths, seqdb.nSeqPos-1);
    }
    else {
      //
      // No seqdb, so there is just one sequence. Still the user specified a title
      // table, so just the first sequence in the fasta file should be used. 
      //
      titleTable.Copy(&fastaGenome.title, 1);
      titleTable.ResetTableToIntegers(&genome.title, &genome.titleLength, 1);
      fastaGenome.titleLength = strlen(genome.title);
    }
    titleTable.Write(titleTableOut);
  }
  else {
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

  ostream  *outFilePtr = &cout;
  ofstream outFileStrm;
  ofstream unalignedFile;
  ostream *unalignedFilePtr = NULL;
  ofstream metricsOut, lcpBoundsOut;
  ofstream anchorFileStrm;
  ofstream clusterOut, *clusterOutPtr;
 
  if (params.anchorFileName != "") {
    CrucialOpen(params.anchorFileName, anchorFileStrm, std::ios::out);
  }

  if (params.clusterFileName != "") {
    CrucialOpen(params.clusterFileName, clusterOut, std::ios::out);
    clusterOutPtr = &clusterOut;
    clusterOut << "total_size p_value n_anchors read_length align_score read_accuracy anchor_probability min_exp_anchors seq_length" << endl;
  }
  else {
    clusterOutPtr = NULL;
  }

  if (params.outFileName != "") {
    CrucialOpen(params.outFileName, outFileStrm, std::ios::out);
    outFilePtr = &outFileStrm;
  }

  if (params.printHeader) {
    switch(params.printFormat) {
    case(SummaryPrint):
      SummaryAlignmentPrinter::PrintHeader(*outFilePtr);
      break;
    case(Interval):
      IntervalAlignmentPrinter::PrintHeader(*outFilePtr);
      break;
    case(CompareSequencesParsable):
      PrintCompareSequencesHeader(*outFilePtr);
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
    //    lcpBoundsOut << "pos depth width lnwidth" << endl;
  }
  
  //
  // Configure the mapping database.
  //

  MappingData<T_SuffixArray, T_GenomeSequence, T_Tuple> *mapdb = new MappingData<T_SuffixArray, T_GenomeSequence, T_Tuple>[params.nProc];

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
  for (procIndex = 0; procIndex < params.nProc; procIndex++ ){
    pthread_attr_init(&threadAttr[procIndex]);
  }

  //
  // Start the mapping jobs.
  //
  int readsFileIndex = 0;
  if (params.subsample < 1) {
    InitializeRandomGeneratorWithTime();
    reader = new ReaderAgglomerate(params.subsample);
  }
  else {
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

  
  if (params.printSAM) {
    string hdString, sqString, rgString, pgString;
    MakeSAMHDString(hdString);
    *outFilePtr << hdString << endl;
    seqdb.MakeSAMSQString(sqString);
    *outFilePtr << sqString; // this already outputs endl
    for (readsFileIndex = 0; readsFileIndex < params.readsFileNames.size()-1; readsFileIndex++ ) {    
      reader->SetReadFileName(params.readsFileNames[readsFileIndex]);
      reader->Initialize();
      string movieName, movieNameMD5;
      reader->GetMovieName(movieName);
      MakeMD5(movieName, movieNameMD5, 10);
      string chipId;
      ParseChipIdFromMovieName(movieName, chipId);
      *outFilePtr << "@RG\t" << "ID:"<<movieNameMD5<<"\t" << "PU:"<< movieName << "\tSM:"<<chipId << endl;
      reader->Close();
    }
    string commandLineString;
    clp.CommandLineToString(argc, argv, commandLineString);
    MakeSAMPGString(commandLineString, pgString);
    *outFilePtr << pgString << endl;
  }

  
  for (readsFileIndex = 0; readsFileIndex < params.readsFileNames.size()-1; readsFileIndex++ ){ 
    params.readsFileIndex = readsFileIndex;
    //
    // Configure the reader to use the correct read and region
    // file names.
    // 
    reader->SetReadFileName(params.readsFileNames[params.readsFileIndex]);

    //
    // Initialize using already set file names.
    //
    int initReturnValue = reader->Initialize();    

    // Check whether use ccs only.
    if (reader->GetFileType() == HDFCCSONLY) {
       params.useAllSubreadsInCcs = false;
       params.useCcs = params.useCcsOnly = true;
    }

    string changeListIdString;
    reader->hdfBasReader.GetChangeListID(changeListIdString);
    ChangeListID changeListId(changeListIdString);
    params.qvScaleType = DetermineQVScaleFromChangeListID(changeListId);
    if (reader->FileHasZMWInformation() and params.useRegionTable) {
      if (params.readSeparateRegionTable) {
        if (regionTableReader->Initialize(params.regionTableFileNames[params.readsFileIndex]) == 0) {
          cout << "ERROR! Could not read the region table " << params.regionTableFileNames[params.readsFileIndex] <<endl;
          exit(1);
        }
        params.useRegionTable = true;
      }
      else {
        if (reader->HasRegionTable()) {
          if (regionTableReader->Initialize(params.readsFileNames[params.readsFileIndex]) == 0) {
            cout << "ERROR! Could not read the region table " << params.readsFileNames[params.readsFileIndex] <<endl;
            exit(1);
          }
          params.useRegionTable = true;
        }
        else {
          params.useRegionTable = false;
        }
      }
    }
    else {
      params.useRegionTable = false;
    }

    /*
     * Check to see if there is a region table. If there is a separate
     * region table, use that (over the region table in the bas
     * file).  If there is a region table in the bas file, use that,
     * without having to specify a region table on the command line. 
     */
    if (params.useRegionTable) {
      regionTable.Reset();
      regionTableReader->ReadTable(regionTable);
      regionTableReader->Close();
      regionTable.SortTableByHoleNumber();
    }

    //
    // Check to see if there is a separate ccs fofn. If there is a separate
    // ccs fofn, use that over the one in the bas file. 
    //
    //if (params.readSeparateCcsFofn and params.useCcs) {
    //  if (reader->SetCCS(params.ccsFofnFileNames[params.readsFileIndex]) == 0) {
    //    cout << "ERROR! Could not read the ccs file " 
    //         << params.ccsFofnFileNames[params.readsFileIndex] << endl;
    //    exit(1);
    //  }
    // }

    if (reader->GetFileType() != HDFCCS and 
        reader->GetFileType() != HDFBase and
        reader->GetFileType() != HDFPulse and
        params.concordant) {
        cout << "WARNING! Option concordant is only enabled when "
             << "input reads are in PacBio base h5 or pulse h5 format." << endl;
        params.concordant = false;
    }


#ifdef USE_GOOGLE_PROFILER
    char *profileFileName = getenv("CPUPROFILE");
    if (profileFileName != NULL) {
      ProfilerStart(profileFileName);
    }
    else {
      ProfilerStart("google_profile.txt");
    }
#endif

    if (initReturnValue > 0) {
      if (params.nProc == 1) {
        mapdb[0].Initialize(&sarray, &genome, &seqdb, &ct, &index, params, reader, &regionTable, 
                            outFilePtr, unalignedFilePtr, &anchorFileStrm, clusterOutPtr);
        mapdb[0].bwtPtr = &bwt;
        if (params.fullMetricsFileName != "") {
          mapdb[0].metrics.SetStoreList(true);
        }
        if (params.lcpBoundsFileName != "") {
          mapdb[0].lcpBoundsOutPtr = &lcpBoundsOut;
        }
        else {
          mapdb[0].lcpBoundsOutPtr = NULL;
        }

        MapReads(&mapdb[0]);
        metrics.Collect(mapdb[0].metrics);
      }
      else {
        pthread_t *threads = new pthread_t[params.nProc];
        for (procIndex = 0; procIndex < params.nProc; procIndex++ ){ 
          //
          // Initialize thread-specific parameters.
          //
            
          mapdb[procIndex].Initialize(&sarray, &genome, &seqdb, &ct, &index, params, reader, &regionTable, 
                                      outFilePtr, unalignedFilePtr, &anchorFileStrm, clusterOutPtr);
          mapdb[procIndex].bwtPtr      = &bwt;
          if (params.fullMetricsFileName != "") {
            mapdb[procIndex].metrics.SetStoreList(true);
          }
          if (params.lcpBoundsFileName != "") {
            mapdb[procIndex].lcpBoundsOutPtr = &lcpBoundsOut;
          }
          else {
            mapdb[procIndex].lcpBoundsOutPtr = NULL;
          }

          if (params.outputByThread) {
            ofstream *outPtr =new ofstream;
            mapdb[procIndex].outFilePtr = outPtr;
            stringstream outNameStream;
            outNameStream << params.outFileName << "." << procIndex;
            mapdb[procIndex].params.outFileName = outNameStream.str();
            CrucialOpen(mapdb[procIndex].params.outFileName, *outPtr, std::ios::out);
          }
          pthread_create(&threads[procIndex], &threadAttr[procIndex], (void* (*)(void*))MapReads, &mapdb[procIndex]);
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
            delete threads;
            threads = NULL;
        }
      }
    }
    reader->Close();
  }
  
  if (!reader) {delete reader; reader = NULL;}

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
    outFileStrm.close();
  }
  cerr << "[INFO] " << GetTimestamp() << " [blasr] ended." << endl;
  return 0;
}
