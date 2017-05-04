#pragma once
/*
 * ============================================================================
 *
 *       Filename:  RegisterOptions.hpp
 *
 *    Description:
 *
 *        Version:  1.0
 *        Created:  04/29/2015 04:48:26 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Yuan Li (yli), yli@pacificbiosciences.com
 *        Company:  Pacific Biosciences
 *
 * ============================================================================
 */

#include <libconfig.h>
#include <CommandLineParser.hpp>
#include <sstream>

#include "MappingParameters.h"
#include "RegisterFilterOptions.h"
using namespace std;

void RegisterBlasrOptions(CommandLineParser& clp, MappingParameters& params)
{
    int trashbinInt;
    float trashbinFloat;
    bool trashbinBool;
    clp.RegisterStringOption("-sa", &params.suffixArrayFileName, "");
    clp.RegisterStringOption("-ctab", &params.countTableName, "");
    clp.RegisterStringOption("-regionTable", &params.regionTableFileName, "");
    clp.RegisterStringOption("-ccsFofn", &params.ccsFofnFileName, "");
    clp.RegisterIntOption("-bestn", (int*)&params.nBest, "", CommandLineParser::PositiveInteger);
    clp.RegisterIntOption("-limsAlign", &params.limsAlign, "", CommandLineParser::PositiveInteger);
    clp.RegisterFlagOption("-printOnlyBest", &params.printOnlyBest, "");
    clp.RegisterFlagOption("-outputByThread", &params.outputByThread, "");
    clp.RegisterFlagOption("-rbao", &params.refineBetweenAnchorsOnly, "");
    clp.RegisterFlagOption("-onegap", &params.separateGaps, "");
    clp.RegisterFlagOption("-allowAdjacentIndels", &params.allowAdjacentIndels, "", false);
    clp.RegisterFlagOption("-placeRepeatsRandomly", &params.placeRandomly, "");
    clp.RegisterFlagOption("-placeGapConsistently", &params.placeGapConsistently, "");
    clp.RegisterIntOption("-randomSeed", &params.randomSeed, "", CommandLineParser::Integer);
    clp.RegisterFlagOption("-extend", &params.extendAlignments, "");
    clp.RegisterIntOption("-branchExpand", &params.anchorParameters.branchExpand, "",
                          CommandLineParser::NonNegativeInteger);
    clp.RegisterIntOption("-maxExtendDropoff", &params.maxExtendDropoff, "",
                          CommandLineParser::NonNegativeInteger);
    clp.RegisterFlagOption("-nucmer", &params.emulateNucmer, "");
    clp.RegisterIntOption("-maxExpand", &params.maxExpand, "", CommandLineParser::PositiveInteger);
    clp.RegisterIntOption("-minExpand", &params.minExpand, "",
                          CommandLineParser::NonNegativeInteger);
    clp.RegisterStringOption("-seqdb", &params.seqDBName, "");
    clp.RegisterStringOption("-anchors", &params.anchorFileName, "");
    clp.RegisterStringOption("-clusters", &params.clusterFileName, "");
    clp.RegisterFlagOption("-samplePaths", (bool*)&params.samplePaths, "");
    clp.RegisterFlagOption("-noStoreMapQV", &params.storeMapQV, "");
    clp.RegisterFlagOption("-nowarp", (bool*)&params.nowarp, "");
    clp.RegisterFlagOption("-guidedAlign", (bool*)&params.useGuidedAlign, "");
    clp.RegisterFlagOption("-useGuidedAlign", (bool*)&trashbinBool, "");
    clp.RegisterFlagOption("-noUseGuidedAlign", (bool*)&params.useGuidedAlign, "");
    clp.RegisterFlagOption("-header", (bool*)&params.printHeader, "");
    clp.RegisterIntOption("-bandSize", &params.bandSize, "", CommandLineParser::PositiveInteger);
    clp.RegisterIntOption("-extendBandSize", &params.extendBandSize, "",
                          CommandLineParser::PositiveInteger);
    clp.RegisterIntOption("-guidedAlignBandSize", &params.guidedAlignBandSize, "",
                          CommandLineParser::PositiveInteger);
    clp.RegisterIntOption("-maxAnchorsPerPosition",
                          (int*)&params.anchorParameters.maxAnchorsPerPosition, "",
                          CommandLineParser::PositiveInteger);
    clp.RegisterIntOption("-stopMappingOnceUnique",
                          (int*)&params.anchorParameters.stopMappingOnceUnique, "",
                          CommandLineParser::NonNegativeInteger);
    clp.RegisterStringOption("-out", &params.outFileName, "");
    clp.RegisterIntOption("-match", &params.match, "", CommandLineParser::Integer);
    clp.RegisterIntOption("-mismatch", &params.mismatch, "", CommandLineParser::Integer);
    clp.RegisterIntOption("-minMatch", &params.minMatchLength, "",
                          CommandLineParser::PositiveInteger);
    clp.RegisterIntOption("-maxMatch", &params.anchorParameters.maxLCPLength, "",
                          CommandLineParser::NonNegativeInteger);
    clp.RegisterIntOption("-maxLCPLength", &params.anchorParameters.maxLCPLength, "",
                          CommandLineParser::NonNegativeInteger);
    clp.RegisterIntOption("-indel", &params.indel, "", CommandLineParser::Integer);
    clp.RegisterIntOption("-insertion", &params.insertion, "", CommandLineParser::Integer);
    clp.RegisterIntOption("-deletion", &params.deletion, "", CommandLineParser::Integer);
    clp.RegisterIntOption("-idsIndel", &params.idsIndel, "", CommandLineParser::Integer);
    clp.RegisterIntOption("-sdpindel", &params.sdpIndel, "", CommandLineParser::Integer);
    clp.RegisterIntOption("-sdpIns", &params.sdpIns, "", CommandLineParser::Integer);
    clp.RegisterIntOption("-sdpDel", &params.sdpDel, "", CommandLineParser::Integer);
    clp.RegisterFloatOption("-indelRate", &params.indelRate, "",
                            CommandLineParser::NonNegativeFloat);
    clp.RegisterFloatOption("-minRatio", &params.minRatio, "", CommandLineParser::NonNegativeFloat);
    clp.RegisterFloatOption("-sdpbypass", &params.sdpBypassThreshold, "",
                            CommandLineParser::NonNegativeFloat);
    clp.RegisterFloatOption("-minFrac", &trashbinFloat, "", CommandLineParser::NonNegativeFloat);
    clp.RegisterIntOption("-maxScore", &params.maxScore, "", CommandLineParser::Integer);
    clp.RegisterStringOption("-bwt", &params.bwtFileName, "");
    clp.RegisterIntOption("m", &params.printFormat, "", CommandLineParser::NonNegativeInteger);
#ifdef USE_PBBAM
    clp.RegisterFlagOption("-sam", &params.printSAM, "");
    clp.RegisterFlagOption("-bam", &params.printBAM, "");
    // BAM read manipulations
    clp.RegisterFlagOption("-polymerase", &params.polymeraseMode, "", false);
#endif
    clp.RegisterStringOption("-clipping", &params.clippingString, "");
    clp.RegisterIntOption("-sdpTupleSize", &params.sdpTupleSize, "",
                          CommandLineParser::PositiveInteger);
    clp.RegisterIntOption("-pvaltype", &params.pValueType, "",
                          CommandLineParser::NonNegativeInteger);
    clp.RegisterIntOption("-start", &params.startRead, "", CommandLineParser::NonNegativeInteger);
    clp.RegisterIntOption("-stride", &params.stride, "", CommandLineParser::NonNegativeInteger);
    clp.RegisterFloatOption("-subsample", &params.subsample, "", CommandLineParser::PositiveFloat);
    clp.RegisterIntOption("-nproc", &params.nProc, "", CommandLineParser::PositiveInteger);
    clp.RegisterFlagOption("-sortRefinedAlignments", (bool*)&params.sortRefinedAlignments, "");
    clp.RegisterIntOption("-quallc", &params.qualityLowerCaseThreshold, "",
                          CommandLineParser::Integer);
    clp.RegisterFlagOption("v", (bool*)&params.verbosity, "");
    clp.RegisterIntOption("V", &params.verbosity, "Specify a level of verbosity.",
                          CommandLineParser::NonNegativeInteger);
    clp.RegisterIntOption("-contextAlignLength", &params.anchorParameters.contextAlignLength, "",
                          CommandLineParser::PositiveInteger);
    clp.RegisterFlagOption("-skipLookupTable", &params.anchorParameters.useLookupTable, "");
    clp.RegisterStringOption("-metrics", &params.metricsFileName, "");
    clp.RegisterStringOption("-lcpBounds", &params.lcpBoundsFileName, "");
    clp.RegisterStringOption("-fullMetrics", &params.fullMetricsFileName, "");
    clp.RegisterIntOption("-nbranch", &params.anchorParameters.numBranches, "",
                          CommandLineParser::NonNegativeInteger);
    clp.RegisterFlagOption("-divideByAdapter", &params.byAdapter, "");
    clp.RegisterFlagOption("-useQuality", &params.ignoreQualities, "");
    clp.RegisterFlagOption("-noFrontAlign", &params.extendFrontAlignment, "");
    clp.RegisterIntOption("-minReadLength", &params.minReadLength, "",
                          CommandLineParser::NonNegativeInteger);
    clp.RegisterIntOption("-maxReadLength", &params.maxReadLength, "",
                          CommandLineParser::NonNegativeInteger);
    clp.RegisterIntOption("-minSubreadLength", &params.minSubreadLength, "",
                          CommandLineParser::NonNegativeInteger);
    clp.RegisterIntOption("-minRawSubreadScore", &params.minRawSubreadScore, "",
                          CommandLineParser::NonNegativeInteger);
    clp.RegisterIntOption("-minAvgQual", &params.minAvgQual, "", CommandLineParser::Integer);
    clp.RegisterFlagOption("-advanceHalf", &params.advanceHalf, "");
    clp.RegisterIntOption("-advanceExactMatches", &params.anchorParameters.advanceExactMatches, "",
                          CommandLineParser::NonNegativeInteger);
    clp.RegisterFlagOption("-useccs", &params.useCcs, "");
    clp.RegisterFlagOption("-useccsdenovo", &params.useCcsOnly, "");
    clp.RegisterFlagOption("-useccsall", &params.useAllSubreadsInCcs, "");
    clp.RegisterFlagOption("-extendDenovoCCSSubreads", &params.extendDenovoCCSSubreads, "");
    clp.RegisterFlagOption("-noRefineAlignments", &params.refineAlignments, "");
    clp.RegisterFlagOption("-refineConcordantAlignments", &params.refineConcordantAlignments, "");
    clp.RegisterIntOption("-nCandidates", &params.nCandidates, "",
                          CommandLineParser::NonNegativeInteger);
    clp.RegisterFlagOption("-useTemp", (bool*)&params.tempDirectory, "");
    clp.RegisterFlagOption("-noSplitSubreads", &params.mapSubreadsSeparately, "");
    clp.RegisterFlagOption("-concordant", &params.concordant, "");
    // When -concordant is turned on, blasr first selects a subread (e.g., the median length full-pass subread)
    // of a zmw as template, maps the template subread to a reference, then infers directions of all other subreads
    // of the same zmw based on direction of the template, and finally maps all other subreads to the same
    // genomic coordinates as the template. When -concordantAlignBothDirections is turned on, blasr will align
    // all other subreads both forwardly and backwardly, without infering their directions. This is a hidden
    // diagnostic option only useful for analyzing movies which have lots of un-identified or missed adapters such
    // that directions of subreads can not be inferred accurately.
    clp.RegisterFlagOption("-concordantAlignBothDirections", &params.concordantAlignBothDirections,
                           "");
    clp.RegisterIntOption("-flankSize", &params.flankSize, "",
                          CommandLineParser::NonNegativeInteger);
    clp.RegisterStringOption("-titleTable", &params.titleTableName, "");
    clp.RegisterFlagOption("-useSensitiveSearch", &params.doSensitiveSearch, "");
    // ignoreRegions or ignoreHQRegions implies region table must exist (i.e., query is HDF).
    clp.RegisterFlagOption("-ignoreRegions", &params.setIgnoreRegions, "");
    clp.RegisterFlagOption("-ignoreHQRegions", &params.setIgnoreHQRegions, "");
    clp.RegisterFlagOption("-computeAlignProbability", &params.computeAlignProbability, "");
    clp.RegisterStringOption("-unaligned", &params.unalignedFileName, "");
    // Print unaligned reads names only
    clp.RegisterFlagOption("-noPrintUnalignedSeqs", &params.noPrintUnalignedSeqs, "");
    clp.RegisterFlagOption("-global", &params.doGlobalAlignment, "");
    clp.RegisterIntOption("-globalChainType", &params.globalChainType, "",
                          CommandLineParser::NonNegativeInteger);
    clp.RegisterFlagOption("-noPrintSubreadTitle", (bool*)&params.printSubreadTitle, "");
    clp.RegisterIntOption("-saLookupTableLength", &params.lookupTableLength, "",
                          CommandLineParser::PositiveInteger);
    clp.RegisterFlagOption("-useDetailedSDP", &params.detailedSDPAlignment, "");
    clp.RegisterFlagOption("-nouseDetailedSDP", &trashbinBool, "");
    clp.RegisterIntOption("-sdpFilterType", &params.sdpFilterType, "",
                          CommandLineParser::NonNegativeInteger);
    clp.RegisterIntOption("-scoreType", &params.scoreType, "",
                          CommandLineParser::NonNegativeInteger);
    clp.RegisterFlagOption("h", &params.printVerboseHelp, "");
    clp.RegisterFlagOption("-help", &params.printDiscussion, "");
    clp.RegisterFloatOption("-accuracyPrior", &params.readAccuracyPrior, "",
                            CommandLineParser::NonNegativeFloat);
    // holeNumberRangesStr is a string of comma-delimited hole number ranges, such as '1,2,3,10-15'.
    // Blasr only analyzes reads whose hole numbers are in the specified hole number ranges.
    clp.RegisterStringOption("-holeNumbers", &params.holeNumberRangesStr, "");
    clp.RegisterIntOption("-substitutionPrior", &params.substitutionPrior, "",
                          CommandLineParser::NonNegativeInteger);
    clp.RegisterIntOption("-deletionPrior", &params.globalDeletionPrior, "",
                          CommandLineParser::NonNegativeInteger);
    clp.RegisterIntOption("-recurseOver", &params.recurseOver, "",
                          CommandLineParser::NonNegativeInteger);
    clp.RegisterStringOption("-scoreMatrix", &params.scoreMatrixString, "");
    clp.RegisterFlagOption("-printDotPlots", &params.printDotPlots, "");
    clp.RegisterFlagOption("-preserveReadTitle", &params.preserveReadTitle, "");
    clp.RegisterFlagOption("-forwardOnly", &params.forwardOnly, "");
    clp.RegisterFlagOption("-affineAlign", &params.affineAlign, "");
    clp.RegisterIntOption("-affineOpen", &params.affineOpen, "",
                          CommandLineParser::NonNegativeInteger);
    clp.RegisterIntOption("-affineExtend", &params.affineExtend, "",
                          CommandLineParser::NonNegativeInteger);
    clp.RegisterFlagOption("-scaleMapQVByNClusters", &params.scaleMapQVByNumSignificantClusters, "",
                           false);
    clp.RegisterFlagOption("-printSAMQV", &params.printSAMQV, "", false);
    clp.RegisterFlagOption("-cigarUseSeqMatch", &params.cigarUseSeqMatch, "", false);
    clp.RegisterStringListOption("-samQV", &params.samQV, "");
    clp.RegisterFlagOption("-fastMaxInterval", &params.fastMaxInterval, "", false);
    clp.RegisterFlagOption("-aggressiveIntervalCut", &params.aggressiveIntervalCut, "", false);
    clp.RegisterFlagOption("-fastSDP", &params.fastSDP, "", false);
    clp.RegisterStringOption("-concordantTemplate", &params.concordantTemplate, "typicalsubread");

    RegisterFilterOptions(clp, params.minAlnLength, params.minPctSimilarity, params.minPctAccuracy,
                          params.hitPolicyStr, trashbinBool = true, trashbinInt, params.maxScore);
}

const string BlasrHelp(MappingParameters& params)
{
    stringstream helpStream;
    helpStream
        << "   Options for blasr " << endl
        << "   Basic usage: 'blasr reads.{bam|fasta|bax.h5|fofn} genome.fasta [-options] " << endl
        << " option\tDescription (default_value)." << endl
        << endl
        << " Input Files." << endl
        << "   reads.bam   is a PacBio BAM file of reads." << endl
        << "               This is the preferred input to blasr because rich quality" << endl
        << "               value (insertion,deletion, and substitution quality values) information "
           "is "
        << endl
        << "               maintained.  The extra quality information improves variant detection "
           "and mapping"
        << endl
        << "               speed." << endl
        << "   reads.fasta is a multi-fasta file of reads.  While any fasta file is valid input, "
        << endl
        << "   reads.bax.h5|reads.plx.h5 is the old DEPRECATED output format of SMRT reads." << endl
        << "   input.fofn  File of file names accepted." << endl
        << endl
        << "   --sa suffixArrayFile" << endl
        << "               Use the suffix array 'sa' for detecting matches" << endl
        << "               between the reads and the reference.  The suffix" << endl
        << "               array has been prepared by the sawriter program." << endl
        << endl
        << "   --ctab tab " << endl
        << "               A table of tuple counts used to estimate match significance.  This is "
        << endl
        << "               by the program 'printTupleCountTable'.  While it is quick to generate "
           "on "
        << endl
        << "               the fly, if there are many invocations of blasr, it is useful to" << endl
        << "               precompute the ctab." << endl
        << endl
        << "   --regionTable table (DEPRECATED)" << endl
        << "               Read in a read-region table in HDF format for masking portions of reads."
        << endl
        << "               This may be a single table if there is just one input file, " << endl
        << "               or a fofn.  When a region table is specified, any region table inside "
        << endl
        << "               the reads.plx.h5 or reads.bax.h5 files are ignored." << endl
        << endl

        << " Options for modifying reads." << endl
        << "   --noSplitSubreads" << endl
        << "               Do not split subreads at adapters. This is typically only " << endl
        << "               useful when the genome in an unrolled version of a known template, and "
        << endl
        << "               contains template-adapter-reverse_template sequence." << endl
#ifdef USE_PBBAM
        << "               For BAM input it reconstitutes full ZMW reads." << endl
        << endl
        << " Option for modifying BAM reads with --noSplitSubreads." << endl
        << "   --polymerase" << endl
        << "               Instead of reconstituting ZMW reads," << endl
        << "               this option reconstitutes polymerase reads, omitting LQ regions." << endl
        << "               Polymerase reads are aligned, if at least one subread is present."
        << endl
        << endl
#endif
        << "(DEPRECATED) Options for modifying HDF reads." << endl
        << "               There is ancilliary information about substrings of reads " << endl
        << "               that is stored in a 'region table' for each read file.  Because " << endl
        << "               HDF is used, the region table may be part of the .bax.h5 or .plx.h5 "
           "file,"
        << endl
        << "               or a separate file.  A contiguously read substring from the template "
        << endl
        << "               is a subread, and any read may contain multiple subreads. The "
           "boundaries "
        << endl
        << "               of the subreads may be inferred from the region table either directly "
           "or "
        << endl
        << "               by definition of adapter boundaries.  Typically region tables also"
        << endl
        << "               contain information for the location of the high and low quality "
           "regions of"
        << endl
        << "               reads.  Reads produced by spurrious reads from empty ZMWs have a high"
        << endl
        << "               quality start coordinate equal to high quality end, making no usable "
           "read."
        << endl
        << "   --useccs   " << endl
        << "               Align the circular consensus sequence (ccs), then report alignments"
        << endl
        << "               of the ccs subreads to the window that the ccs was mapped to.  Only "
        << endl
        << "               alignments of the subreads are reported." << endl
        << "   --useccsall" << endl
        << "               Similar to -useccs, except all subreads are aligned, rather than just"
        << endl
        << "               the subreads used to call the ccs.  This will include reads that only"
        << endl
        << "               cover part of the template." << endl
        << "   --useccsdenovo" << endl
        << "               Align the circular consensus, and report only the alignment of the ccs"
        << endl
        << "               sequence." << endl
        << "   --ignoreRegions(false)" << endl
        << "               Ignore any information in the region table." << endl
        << "   --ignoreHQRegions (false)Ignore any hq regions in the region table." << endl
        << endl
        << " Alignments To Report." << endl
        << "   --bestn n (10)" << endl
        << "               Report the top 'n' alignments." << endl
        << "   --hitPolicy" << endl
        << "               " << params.hitPolicy.Help(string(15, ' ')) << endl
        << "   --placeRepeatsRandomly (false)" << endl
        << "               DEPRECATED! If true, equivalent to --hitPolicy randombest." << endl
        << "   --placeGapConsistently (false)" << endl
        << "               Place gaps consistently in alignments of a read as alignments " << endl
        << "               of its reverse complementary sequence." << endl
        << "   --randomSeed (0)" << endl
        << "               Seed for random number generator. By default (0), use current time as "
           "seed. "
        << endl
        << "   --noSortRefinedAlignments (false) " << endl
        << "               Once candidate alignments are generated and scored via sparse dynamic "
        << endl
        << "               programming, they are rescored using local alignment that accounts "
        << endl
        << "               for different error profiles." << endl
        << "               Resorting based on the local alignment may change the order the hits "
           "are returned."
        << endl
        << "   --allowAdjacentIndels " << endl
        << "               When specified, adjacent insertion or deletions are allowed. Otherwise, "
           "adjacent "
        << endl
        << "               insertion and deletions are merged into one operation.  Using quality "
           "values "
        << endl
        << "               to guide pairwise alignments may dictate that the higher probability "
           "alignment "
        << endl
        << "               contains adjacent insertions or deletions.  Current tools such as GATK "
           "do not permit"
        << endl
        << "               this and so they are not reported by default." << endl
        << endl
        << " Output Formats and Files" << endl
        << "   --out out (terminal)  " << endl
        << "               Write output to 'out'." << endl
#ifdef USE_PBBAM
        << "   --bam       Write output in PacBio BAM format. This is the preferred output format."
        << endl
        << "               Input query reads must be in PacBio BAM format." << endl
#endif
        << "   --sam       Write output in SAM format. Starting from version 5.2 is no longer "
           "supported"
        << endl
        << "               Use --bam, then translate from .bam to .sam" << endl
        << "   -m t           " << endl
        << "               If not printing SAM, modify the output of the alignment." << endl
        << "                t=" << StickPrint
        << " Print blast like output with |'s connecting matched nucleotides." << endl
        << "                  " << SummaryPrint << " Print only a summary: score and pos." << endl
        << "                  " << CompareXML << " Print in Compare.xml format." << endl
        << "                  " << Vulgar << " Print in vulgar format (DEPRECATED)." << endl
        << "                  " << Interval << " Print a longer tabular version of the alignment."
        << endl
        << "                  " << CompareSequencesParsable
        << " Print in a machine-parsable format that is read by compareSequences.py." << endl
        << "   --header" << endl
        << "               Print a header as the first line of the output file describing the "
           "contents of each column."
        << endl
        << "   --titleTable tab (NULL) " << endl
        << "               Construct a table of reference sequence titles.  The reference "
           "sequences are "
        << endl
        << "               enumerated by row, 0,1,...  The reference index is printed in alignment "
           "results"
        << endl
        << "               rather than the full reference name.  This makes output concise, "
           "particularly when"
        << endl
        << "               very verbose titles exist in reference names." << endl
        << "   --unaligned file" << endl
        << "               Output reads that are not aligned to 'file'" << endl
        << "   --noPrintUnalignedSeqs" << endl
        << "               Must be used together with -unaligned, print unaligned read names only."
        << endl
        << "   --clipping [none|hard|subread|soft] (none)" << endl
        << "               Use no/hard/subread/soft clipping, ONLY for SAM/BAM output." << endl
        << "   --printSAMQV (false)" << endl
        << "               Print quality values to SAM output." << endl
        //             << "   --cigarUseSeqMatch (false)" << endl
        //             << "               CIGAR strings in SAM/BAM output use '=' and 'X' to represent sequence match and mismatch instead of 'M'." << endl << endl
        << " Options for anchoring alignment regions. This will have the greatest effect on speed "
           "and sensitivity."
        << endl
        << "   --minMatch m (12) " << endl
        << "               Minimum seed length.  Higher minMatch will speed up alignment, " << endl
        << "               but decrease sensitivity." << endl
        //             << "   --maxExpand M (1)" << endl
        //             << "               Perform no more than M iterations of searches through the suffix " << endl
        //             << "               array for matches. At each iteration, all matches of length LCPi-M" << endl
        //             << "               are found, where LCPi is the length of the longest common prefix " << endl
        //             << "               between the string at i and anywhere in the genome."<<endl
        //             << "               The number of matches grows as M increases, and can become very large with M > 3." << endl
        << "   --maxMatch l (inf)" << endl
        << "               Stop mapping a read to the genome when the lcp length reaches l.  "
        << endl
        << "               This is useful when the query is part of the reference, for example "
           "when "
        << endl
        << "               constructing pairwise alignments for de novo assembly." << endl
        << "   --maxLCPLength l (inf)" << endl
        << "               The same as -maxMatch." << endl
        << "   --maxAnchorsPerPosition m (10000) " << endl
        << "               Do not add anchors from a position if it matches to more than 'm' "
           "locations in the target."
        << endl
        //             << "   --advanceHalf (false) " << endl
        //             << "               A trick for speeding up alignments at the cost of sensitivity.  If " << endl
        //             << "               a cluster of anchors of size n, (a1,...,an) is found, normally anchors " << endl
        //             << "               (a2,...an) of size n-1 is also clustered to make sure a1 did not decrease the " << endl
        //             << "               cluster score.  When advanceHalf is specified, clustering begins at a_(n/2)."<<endl<< endl
        << "   --advanceExactMatches E (0)" << endl
        << "               Another trick for speeding up alignments with match - E fewer anchors.  "
           "Rather than"
        << endl
        << "               finding anchors between the read and the genome at every position in "
           "the read, "
        << endl
        << "               when an anchor is found at position i in a read of length L, the next "
           "position "
        << endl
        << "               in a read to find an anchor is at i+L-E." << endl
        << "               Use this when alignining already assembled contigs." << endl
        << "   --nCandidates n (10)" << endl
        << "               Keep up to 'n' candidates for the best alignment.  A large value of n "
           "will slow mapping"
        << endl
        << "               because the slower dynamic programming steps are applied to more "
           "clusters of anchors"
        << endl
        << "               which can be a rate limiting step when reads are very long." << endl
        << "   --concordant(false)" << endl
        << "               Map all subreads of a zmw (hole) to where the longest full pass subread "
           "of the zmw "
        << endl
        << "               aligned to. This requires to use the region table and hq regions."
        << endl
        << "               This option only works when reads are in base or pulse h5 format."
        << endl
        << "   --fastMaxInterval(false)" << endl
        << "               Fast search maximum increasing intervals as alignment candidates. The "
           "search "
        << endl
        << "               is not as exhaustive as the default, but is much faster." << endl
        << "   --aggressiveIntervalCut(false)" << endl
        << "               Agreesively filter out non-promising alignment candidates, if there "
        << endl
        << "               exists at least one promising candidate. If this option is turned on, "
        << endl
        << "               Blasr is likely to ignore short alignments of ALU elements." << endl
        << "   --fastSDP(false)" << endl
        << "               Use a fast heuristic algorithm to speed up sparse dynamic programming."
        << endl
        << endl
        << "  Options for Refining Hits." << endl
        //             << "   --indelRate i (0.30)" << endl
        //             << "               The approximate maximum rate to allow drifting from the diagonal." <<endl << endl
        << "   --refineConcordantAlignments(false)" << endl
        << "               Refine concordant alignments. It slightly increases alignment accuracy "
           "at cost of time."
        << endl
        << "   --sdpTupleSize K (11)" << endl
        << "               Use matches of length K to speed dynamic programming alignments.  This "
           "controls"
        << endl
        << "               accuracy of assigning gaps in pairwise alignments once a mapping has "
           "been found,"
        << endl
        << "               rather than mapping sensitivity itself." << endl
        << "   --scoreMatrix \"score matrix string\" " << endl
        << "               Specify an alternative score matrix for scoring fasta reads.  The "
           "matrix is "
        << endl
        << "               in the format " << endl
        << "                  ACGTN" << endl
        << "                A abcde" << endl
        << "                C fghij" << endl
        << "                G klmno" << endl
        << "                T pqrst" << endl
        << "                N uvwxy"
        << " . The values a...y should be input as a quoted space separated " << endl
        << "               string: \"a b c ... y\". Lower scores are better, so matches should be "
           "less "
        << endl
        << "               than mismatches e.g. a,g,m,s = -5 (match), mismatch = 6. " << endl
        << "   --affineOpen value (10) " << endl
        << "               Set the penalty for opening an affine alignment." << endl
        << "   --affineExtend a (0)" << endl
        << "               Change affine (extension) gap penalty. Lower value allows more gaps."
        << endl
        << endl
        << " Options for overlap/dynamic programming alignments and pairwise overlap for de novo "
           "assembly. "
        << endl
        << "   --useQuality (false)" << endl
        << "               Use substitution/insertion/deletion/merge quality values to score gap "
           "and "
        << endl
        << "               mismatch penalties in pairwise alignments.  Because the insertion and "
           "deletion"
        << endl
        << "               rates are much higher than substitution, this will make many alignments "
        << endl
        << "               favor an insertion/deletion over a substitution.  Naive consensus "
           "calling methods "
        << endl
        << "               will then often miss substitution polymorphisms. This option should be "
        << endl
        << "               used when calling consensus using the Quiver method.  Furthermore, when "
        << endl
        << "               not using quality values to score alignments, there will be a lower "
           "consensus "
        << endl
        << "               accuracy in homolymer regions." << endl
        << "   --affineAlign (false)" << endl
        << "               Refine alignment using affine guided align." << endl
        << endl
        << " Options for filtering reads and alignments" << endl
        << "   --minReadLength l(50)" << endl
        << "               Skip reads that have a full length less than l. Subreads may be shorter."
        << endl
        << "   --minSubreadLength l(0)" << endl
        << "               Do not align subreads of length less than l." << endl
        << "   --minRawSubreadScore m(0)" << endl
        << "               Do not align subreads whose quality score in region table is less than "
           "m (quality scores should be in range [0, 1000])."
        << endl
        << "   --maxScore m(-200)" << endl  //params.filterCriteria.scoreCutoff
        << "               Maximum score to output (high is bad, negative good)." << endl
        << "   --minAlnLength" << endl
        << "               " << params.filterCriteria.MinAlnLengthHelp() << endl
        << "   --minPctSimilarity" << endl
        << "               " << params.filterCriteria.MinPctSimilarityHelp() << endl
        << "   --minPctAccuracy" << endl
        << "               " << params.filterCriteria.MinPctAccuracyHelp() << endl
        << endl
        << " Options for parallel alignment." << endl
        << "   --nproc N (1)" << endl
        << "               Align using N processes.  All large data structures such as the suffix "
           "array and "
        << endl
        << "               tuple count table are shared." << endl
        << "   --start S (0)" << endl
        << "               Index of the first read to begin aligning. This is useful when multiple "
           "instances "
        << endl
        << "               are running on the same data, for example when on a multi-rack cluster."
        << endl
        << "   --stride S (1)" << endl
        << "               Align one read every 'S' reads." << endl
        << endl
        << " Options for subsampling reads." << endl
        << "   --subsample (0)" << endl
        << "               Proportion of reads to randomly subsample (expressed as a decimal) and "
           "align."
        << endl
        << "   --holeNumbers LIST " << endl
        << "               When specified, only align reads whose ZMW hole numbers are in LIST."
        << endl
        << "               LIST is a comma-delimited string of ranges, such as '1,2,3,10-13'."
        << endl
        << "               This option only works when reads are in bam, bax.h5 or plx.h5 format."
        << endl
        << endl
        //             << " Options for dynamic programming alignments. " << endl << endl
        //             << "   --ignoreQuality" << endl
        //             << "                 Ignore quality values when computing alignments (they still may be used." << endl
        //             << "                 when mapping)." << endl << endl
        //             << " -v            Print some verbose information." << endl
        //             << " -V 2          Make verbosity more verbose.  Probably only useful for development." << endl
        << " -h            Print this help file." << endl
        << endl
        << "In release v5.1 of BLASR, command-line options will use the " << endl
        << "single dash/double dash convention: " << endl
        << "Character options are preceded by a single dash. (Example: -v) " << endl
        << "Word options are preceded by a double dash. (Example: --verbose) " << endl
        << "Please modify your scripts accordingly when BLASR v5.1 is released. " << endl
        << endl
        << "To cite BLASR, please use: Chaisson M.J., and Tesler G., Mapping " << endl
        << "single molecule sequencing reads using Basic Local Alignment with " << endl
        << "Successive Refinement (BLASR): Theory and Application, BMC " << endl
        << "Bioinformatics 2012, 13:238." << endl
        << "Please report any bugs to "
        << "'https://github.com/PacificBiosciences/blasr/issues'." << endl
        << endl;
    return helpStream.str();
}

const string BlasrConciseHelp(void)
{
    stringstream ss;
    ss << "blasr - a program to map reads to a genome" << endl
       << " usage: blasr reads genome " << endl
       << " Run with -h for a list of commands " << endl
       << "          --help for verbose discussion of how to run blasr." << endl
       << endl
       << "In release v5.1 of BLASR, command-line options will use the " << endl
       << "single dash/double dash convention: " << endl
       << "Character options are preceded by a single dash. (Example: -v) " << endl
       << "Word options are preceded by a double dash. (Example: --verbose) " << endl
       << "Please modify your scripts accordingly when BLASR v5.1 is released. " << endl
       << endl;
    return ss.str();
}

const string BlasrSummaryHelp(void)
{
    stringstream ss;
    ss << "   Basic usage: 'blasr reads.{bam|fasta|bax.h5|fofn} genome.fasta [-options] " << endl
       << " [option]\tDescription (default_value)." << endl
       << endl
       << " Input Files." << endl
       << "   reads.bam is the NEW native output format for SMRT reads."
          "This is the preferred input to blasr because rich quality"
          "value (insertion,deletion, and substitution quality values) information is "
          "maintained.  The extra quality information improves variant detection and mapping"
       << "speed." << endl
       << "   reads.fasta is a multi-fasta file of reads.  While any fasta file is valid input, "
          "it is preferable to use bax.h5 or plx.h5 files because they contain "
          "more rich quality value information."
       << endl
       << "   reads.bax.h5|reads.plx.h5 is the OLD (DEPRECATED) output format of "
          "SMRT reads. "
       << endl
       << "   reads.fofn File of file names accepted." << endl
       << endl;
    return ss.str();
}

const string BlasrDiscussion(void)
{
    stringstream ss;
    ss << "NAME" << endl
       << "         blasr - Map SMRT Sequences to a reference genome." << endl
       << endl
       << "SYNOPSIS" << endl
       << "         blasr reads.bam genome.fasta --bam --out out.bam" << endl
       << endl
       << "         blasr reads.fasta genome.fasta " << endl
       << endl
       << "         blasr reads.fasta genome.fasta --sa genome.fasta.sa" << endl
       << endl
       << "         blasr reads.bax.h5 genome.fasta [--sa genome.fasta.sa] " << endl
       << endl
       << "         blasr reads.bax.h5 genome.fasta --sa genome.fasta.sa --maxScore 100 --minMatch "
          "15 ... "
       << endl
       << endl
       << "         blasr reads.bax.h5 genome.fasta --sa genome.fasta.sa --nproc 24 --out "
          "alignment.out ... "
       << endl
       << endl
       << "DESCRIPTION " << endl
       << "  blasr is a read mapping program that maps reads to positions " << endl
       << "  in a genome by clustering short exact matches between the read and" << endl
       << "  the genome, and scoring clusters using alignment. The matches are" << endl
       << "  generated by searching all suffixes of a read against the genome" << endl
       << "  using a suffix array. Global chaining methods are used to score " << endl
       << "  clusters of matches." << endl
       << endl
       << "  The only required inputs to blasr are a file of reads and a" << endl
       << "  reference genome.  It is exremely useful to have read filtering" << endl
       << "  information, and mapping runtime may decrease substantially when a" << endl
       << "  precomputed suffix array index on the reference sequence is" << endl
       << "  specified." << endl
       << "  " << endl
       << "  Although reads may be input in FASTA format, the recommended input is" << endl
       << "  PacBio BAM files because these contain qualtiy value" << endl
       << "  information that is used in the alignment and produces higher quality" << endl
       << "  variant detection." << endl
       << "  Although alignments can be output in various formats, the recommended " << endl
       << "  output format is PacBio BAM." << endl
       << "  Support to bax.h5 and plx.h5 files will be DEPRECATED." << endl
       << "  Support to region tables for h5 files will be DEPRECATED." << endl
       //<< "  Read filtering information is contained in the .bax.h5 input files as" << endl
       //<< "  well as generated by other post-processing programs with analysis of" << endl
       //<< "  pulse files and read in from a separate .region.h5 file.  The current" << endl
       //<< "  set of filters that are applied to reads are high quality region" << endl
       //<< "  filtering, and adapter filtering.  Regions outside high-quality" << endl
       //<< "  regions are ignored in mapping.  Reads that contain regions annotated" << endl
       //<< "  as adapter are split into non-adapter (template) regions, and mapped" << endl
       //<< "  separately." << endl
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
       << endl
       << endl;
    return ss.str();
}
