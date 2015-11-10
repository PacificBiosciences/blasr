#pragma once

#define REQUIRE_PBBAM_ERROR() \
assert("blasr must be compiled with lib pbbam to perform IO on bam." == 0);

#include <vector>

#include <reads/ReadType.hpp>
#include <utils/FileOfFileNames.hpp>
#include <utils/RangeUtils.hpp>
#include <tuples/TupleMetrics.hpp>
#include <datastructures/anchoring/AnchorParameters.hpp>
#include <qvs/QualityValue.hpp>
#include <format/SAMPrinter.hpp>
#include <algorithms/alignment/AlignmentFormats.hpp>
#include <files/BaseSequenceIO.hpp>
#include <datastructures/alignment/FilterCriteria.hpp>

class MappingParameters {
public:
    //
    // Parameters for global substitution, insertion, and deletion priors.
    //
    float minFractionToBeConsideredOverlapping;
    float indelRate;
    float minRatio;
    int indel;
    int idsIndel;
    int sdpIndel;
    int sdpIns, sdpDel;
    int insertion;
    int deletion;
    int mismatch;
    int sdpTupleSize;
    int match;
    int showAlign;
    bool useScoreCutoff;
    int maxScore;
    int argi;
    int nProc;
    int globalChainType;
    SAMOutput::Clipping clipping;
    string clippingString;
    QVScale qvScaleType;
    vector<string> readsFileNames; // = queryFileNames, genomeFileName
    vector<string> queryFileNames;
    string genomeFileName;
    // Query file type: FASTA/FASTQ/HDF*/PBBAM,
    // Note that mixed query file types is not allowed.
    FileType queryFileType;
    // Query read type, SUBREAD, CCS or UNKNOWN
    // Note that mixed read types is not allowed.
    ReadType::ReadTypeEnum queryReadType;
    vector<string> regionTableFileNames;
    vector<string> ccsFofnFileNames;
    string tupleListName;
    string posTableName;
    string outFileName;
    string suffixArrayFileName;
    string bwtFileName;
    string indexFileName;
    string anchorFileName;
    string clusterFileName;
    VectorIndex nBest;
    int printWindow;
    int doCondense;
    int do4BitComp;
    int cutoff;
    int useSuffixArray;
    int useBwt;
    int useReverseCompressIndex;
    int useTupleList;
    int useSeqDB;
    string seqDBName;
    int useCountTable;
    string countTableName;
    int minMatchLength;
    int listTupleSize;
    int printFormat;
    int maxExpand, minExpand;
    int startRead;
    int stride;
    int pValueType;
    float subsample;
    int sortRefinedAlignments;
    int verbosity;
    bool printSAM;
    bool cigarUseSeqMatch;
    bool printBAM;
    bool storeMapQV;
    bool useRandomSeed;
    int  randomSeed;
    bool placeRandomly;
    bool printHeader;
    bool samplePaths;
    bool warp, nowarp;
    //bool usePrefixLookupTable;
    bool doSensitiveSearch;
    bool emulateNucmer;
    bool refineBetweenAnchorsOnly;
    bool byAdapter;
    bool extendDenovoCCSSubreads;
    TupleMetrics saTupleMetrics;
    TupleMetrics sdpTupleMetrics;
    int lookupTableLength;
    //int branchQualityThreshold;
    int qualityLowerCaseThreshold;
    AnchorParameters anchorParameters;
    int readsFileIndex;
    //int numBranches;
    bool storeMetrics;
    bool ignoreQualities;
    bool extendFrontAlignment;
    bool extendAlignments;
    int  maxExtendDropoff;
    int  minReadLength;
    int  maxReadLength;
    int  minSubreadLength;
    int  minRawSubreadScore;
    int  minAvgQual;
    bool overlap;
    bool advanceHalf;
    int advanceExactMatches;
    float approximateMaxInsertionRate;
    float minPctSimilarity; // [0, 100]
    float minPctAccuracy; // [0, 100]
    bool refineAlignments;
    int nCandidates;
    bool doGlobalAlignment;
    string tempDirectory;
    bool useTitleTable;
    string titleTableName;
    bool readSeparateRegionTable;
    bool readSeparateCcsFofn;
    string regionTableFileName;
    string ccsFofnFileName;
    //float averageMismatchScore;
    bool mapSubreadsSeparately;
    bool concordant;
    int  flankSize;
    bool useRegionTable;
    bool useHQRegionTable;
    bool printUnaligned;
    string unalignedFileName;
    string metricsFileName;
    string lcpBoundsFileName;
    string fullMetricsFileName;
    bool printSubreadTitle;
    bool useCcs;
    bool useAllSubreadsInCcs;
    bool useCcsOnly;
    bool detailedSDPAlignment, nouseDetailedSDPAlignment;
    int  chunkSize;
    int  sdpFilterType;
    bool useGuidedAlign;
    int  guidedAlignBandSize;
    int  bandSize;
    int  extendBandSize;
    bool useQVScore;
    int  scoreType;
    bool printVerboseHelp;
    bool printDiscussion;
    float sdpBypassThreshold;
    bool computeAlignProbability;
    float qvMatchWeight;
    float qvMismatchWeight;
    float qvInsWeight;
    float qvDelWeight;
    float readAccuracyPrior;
    bool  printVersion;
    int   substitutionPrior;
    int   globalDeletionPrior;
    bool  outputByThread;
    int   recurseOver;
    bool  allowAdjacentIndels;
    bool  separateGaps;
    string scoreMatrixString;
    bool  printDotPlots;
    bool  preserveReadTitle;
    bool  forwardOnly;
    bool  printOnlyBest;
    bool  affineAlign;
    int   affineExtend;
    int   affineOpen;
    bool  scaleMapQVByNumSignificantClusters;
    int   limsAlign;
    string holeNumberRangesStr;
    Ranges holeNumberRanges;
    int minAlnLength;
    bool printSAMQV;
    vector<string> samQV;
    SupplementalQVList samQVList;
    bool fastMaxInterval;
    bool aggressiveIntervalCut;
    bool fastSDP;
    string concordantTemplate;
    bool concordantAlignBothDirections;
    FilterCriteria filterCriteria;
    string hitPolicyStr;
    HitPolicy hitPolicy;
    bool enableHiddenPaths;

    void Init() {
        qvMatchWeight = 1.0;
        qvMismatchWeight = 1.0;
        qvInsWeight = 1.0;
        qvDelWeight = 1.0;
        minFractionToBeConsideredOverlapping = 0.75;
        minRatio = 0.25;
        indelRate = 0.3;
        indel    = 5;
        insertion = 4; // asymmetric indel parameters
        deletion  = 5;
        idsIndel = 15;
        sdpIndel = 5;
        sdpIns   = 5;
        sdpDel   = 10;
        sdpTupleSize = 11;
        match = 0;
        mismatch = 0;
        showAlign = 1;
        useScoreCutoff = false;
        maxScore = -200;
        argi = 1;
        nProc = 1;
        readsFileNames.clear();
        queryFileNames.clear();
        genomeFileName = "";
        queryReadType = ReadType::UNKNOWN;
        queryFileType = FileType::None;
        tupleListName = "";
        posTableName = "";
        suffixArrayFileName= "";
        bwtFileName = "";
        indexFileName = "";
        anchorFileName = "";
        outFileName = "";
        nBest = 10;
        nCandidates = 10;
        printWindow = 0;
        doCondense  = 0;
        do4BitComp  = 0;
        pValueType = 0;
        cutoff = 0;
        useSuffixArray = 0;
        useBwt = 0;
        useReverseCompressIndex = 0;
        useTupleList = 0;
        useSeqDB = 0;
        seqDBName = "";
        useCountTable = 0;
        countTableName = "";
        lookupTableLength = 8;
        anchorParameters.minMatchLength = minMatchLength = 12;
        printFormat = SummaryPrint;
        maxExpand = 0;
        minExpand = 0;
        startRead = 0;
        stride    = 1;
        subsample = 1.1;
        listTupleSize = 6;
        sortRefinedAlignments = 1;
        anchorParameters.verbosity = verbosity = 0;
        saTupleMetrics.Initialize(listTupleSize);
        sdpTupleMetrics.Initialize(sdpTupleSize);
        qualityLowerCaseThreshold = 0;
        anchorParameters.branchQualityThreshold = 0;
        readsFileIndex = 0;
        printSAM = false;
        printBAM = false;
        useRandomSeed = false;
        randomSeed = 0;
        placeRandomly = false;
        samplePaths = false;
        nowarp = false;
        storeMapQV = true;
        warp = true;
        extendDenovoCCSSubreads = false;
        storeMetrics = false;
        ignoreQualities = true;
        extendFrontAlignment = false;
        extendAlignments = false;
        maxExtendDropoff = 10;
        minReadLength = 50;
        maxReadLength = 0; // means no max read length
        minSubreadLength = 0;
        minRawSubreadScore = -1; // raw subread score in region table should be in range [0, 1000].
        minAvgQual = 0;
        overlap = false;
        advanceHalf = false;
        refineAlignments = true;
        anchorParameters.advanceExactMatches = advanceExactMatches = 0;
        approximateMaxInsertionRate = 1.30;
        minPctSimilarity = 0;
        minPctAccuracy = 0;
        doGlobalAlignment = false;
        tempDirectory = "";
        useTitleTable = false;
        titleTableName  = "";
        readSeparateRegionTable = false;
        readSeparateCcsFofn = false;
        regionTableFileName = "";
        ccsFofnFileName = "";
        mapSubreadsSeparately=true;
        concordant=false;
        flankSize=40;
        useRegionTable = true;
        useHQRegionTable=true;
        printUnaligned = false;
        unalignedFileName = "";
        globalChainType = 0;
        metricsFileName = "";
        fullMetricsFileName = "";
        doSensitiveSearch = false;
        emulateNucmer = false;
        refineBetweenAnchorsOnly = false;
        printSubreadTitle = true;
        detailedSDPAlignment = true;
        nouseDetailedSDPAlignment = false;
        useCcs     = false;
        useCcsOnly = false;
        useAllSubreadsInCcs = false;
        chunkSize = 10000000;
        sdpFilterType = 0;
        anchorParameters.stopMappingOnceUnique = true;
        useGuidedAlign = true;
        bandSize = 0;
        extendBandSize = 10;
        guidedAlignBandSize = 10;
        useQVScore = false;
        printVerboseHelp = false;
        printDiscussion = false;
        sdpBypassThreshold = 1000000.0;
        scoreType = 0;
        byAdapter = false;
        qvScaleType = PHRED;
        printHeader = false;
        computeAlignProbability = false;
        readAccuracyPrior = 0.85;
        printVersion = false;
        clipping = SAMOutput::none;
        clippingString = "";
        substitutionPrior = 20;
        globalDeletionPrior = 13;
        outputByThread = false;
        recurseOver = 10000;
        allowAdjacentIndels = false;
        separateGaps = false;
        scoreMatrixString = "";
        printDotPlots = false;
        preserveReadTitle = false;
        forwardOnly = false;
        printOnlyBest = false;
        affineAlign = false;
        affineExtend = 0;
        affineOpen = 10;
        scaleMapQVByNumSignificantClusters = false;
        limsAlign = 0;
        holeNumberRangesStr = "";
        minAlnLength = 0;
        printSAMQV = false;
        cigarUseSeqMatch = false;
        samQV.clear();
        samQVList.clear();
        fastMaxInterval = false;
        aggressiveIntervalCut = false;
        fastSDP = false;
        concordantTemplate = "mediansubread"; // typicalsubread or longestsubread
        concordantAlignBothDirections = false;

        hitPolicyStr = "all";
        ResetFilterAndHit();
        enableHiddenPaths = false; //turn off hidden paths.
    }

    MappingParameters()
        : filterCriteria(0, 0, 0, false, Score(0, ScoreSign::NEGATIVE))
        , hitPolicy("all", ScoreSign::NEGATIVE)
    {
        Init();
    }

    void MakeSane() {
        // Expand FOFN
        FileOfFileNames::ExpandFileNameList(readsFileNames);

        // Must have at least a query and a genome
        if (readsFileNames.size() <= 1) {
            cout << "Error, you must provide at least one reads file and a genome file." <<endl;
            exit(1);
        }

        // Separate query reads files and a genome read file
        // The last reads file is the genome
        queryFileNames = readsFileNames;
        queryFileNames.pop_back();
        genomeFileName = readsFileNames.back();

        // Check query file type.
        BaseSequenceIO::DetermineFileTypeByExtension(queryFileNames[0], queryFileType);
        for (size_t i = 1; i < queryFileNames.size(); i++) {
            FileType fileType;
            BaseSequenceIO::DetermineFileTypeByExtension(queryFileNames[i], fileType);
            if (fileType != queryFileType) {
                cout << "ERROR, mixed query file types is not allowed." << endl;
                exit(1);
            }
        }

        // -useQuality can not be used in combination with a fasta input
        if (!ignoreQualities) {
            if (queryFileType == Fasta) {
                cout<<"ERROR, you can not use -useQuality option when any of the input reads files are in multi-fasta format."<<endl;
                exit(1);
            }
        }

        //
        // Fix all logical incompatibilities with parameters.
        //
        if (nowarp) {
            warp = false;
        }

        if (nCandidates < nBest) {
            cerr << "Warning: resetting nCandidates to nBest " << nBest << endl;
            nCandidates = nBest;
        }

        if (placeRandomly and hitPolicyStr != "randombest") {
            cerr << "Warning: placeRepeatsRandomly is deprecated, resetting hit policy to randombest." << endl;
            hitPolicyStr = "randombest";
        }
        if ((hitPolicyStr == "random" or hitPolicyStr == "randombest") and nBest == 1) {
            cerr << "Warning: When attempting to select equivalently scoring reads at random " << endl
                << "the bestn parameter should be greater than one." << endl;
        }

        if (concordant) {
            if (useCcs) {
                concordant = false;
            } else {
                useRegionTable   = true;
                useHQRegionTable = true;
            }
            if (concordantTemplate != "longestsubread" and concordantTemplate != "typicalsubread" and concordantTemplate != "mediansubread") {
                cout << "ERROR, unsupported concordantTemplate: " << concordantTemplate << endl;
                exit(1);
            }
        }

        if (sdpFilterType > 1) {
            cerr << "Warning: using new filter method for SDP alignments.  The parameter is " << endl
                << "either 0 or 1, but " << sdpFilterType << " was specified." << endl;
            sdpFilterType = 1;
        }
        if (sdpFilterType == 0) {
            detailedSDPAlignment = true;
            nouseDetailedSDPAlignment = false;
        }
        if (detailedSDPAlignment == false) {
            sdpFilterType = 1;
        }
        if (useGuidedAlign == true and bandSize == 0) {
            bandSize = 16;
        }
        anchorParameters.minMatchLength = minMatchLength;
        if (suffixArrayFileName != "") {
            useSuffixArray = true;
        }
        if (bwtFileName != "") {
            useBwt = true;
        }
        if (useBwt and useSuffixArray) {
            cout << "ERROR, sa and bwt must be used independently." << endl;
            exit(1);
        }
        if (countTableName != "") {
            useCountTable = true;
        }
        if (metricsFileName != "" or fullMetricsFileName != "") {
            storeMetrics = true;
        }
        if (useCcsOnly) {
            useCcs = true;
        }
        if (useAllSubreadsInCcs == true) {
            useCcs = true;
        }
        if (titleTableName != "") {
            useTitleTable = true;
        }
        if (unalignedFileName != "") {
            printUnaligned = true;
        }
        if (regionTableFileName != "") {
            useRegionTable = true;
            readSeparateRegionTable = true;
        }
        if (ccsFofnFileName != "") {
            readSeparateCcsFofn = true;
        }
        if (nouseDetailedSDPAlignment == true) {
            detailedSDPAlignment = false;
        }
        if (nouseDetailedSDPAlignment == false) {
            detailedSDPAlignment = true;
        }
        if (anchorParameters.maxLCPLength != 0 and anchorParameters.maxLCPLength < anchorParameters.minMatchLength) {
            cerr << "ERROR: maxLCPLength is less than minLCPLength, which will result in no hits." << endl;
        }
        if (subsample < 1 and stride > 1) {
            cout << "ERROR, subsample and stride must be used independently." << endl;
            exit(1);
        }


        if (emulateNucmer) {
            SetEmulateNucmer();
        }

        if (randomSeed != 0) {
            useRandomSeed = true;
        }
        if (printSAM) {
            printFormat = SAM;
        }
        //
        // Parse the clipping.
        //
        if (clippingString == "soft") {
            clipping = SAMOutput::soft;
        }
        else if (clippingString == "hard") {
            clipping = SAMOutput::hard;
        }
        else if (clippingString == "none") {
            clipping = SAMOutput::none;
        }
        else if (clippingString == "subread") {
            clipping = SAMOutput::subread;
        }
        else if (clippingString != "") {
            cout << "ERROR, clipping should either be soft, hard, or none." << endl;
            exit(1);
        }

        if (printBAM) {
#ifndef USE_PBBAM
            REQUIRE_PBBAM_ERROR();
#else
            cigarUseSeqMatch = true; // ALWAYS true for BAM
            printFormat = BAM;
            printSAM = false;
            samQVList.SetDefaultQV();
            printSAMQV = true;
            if (clipping != SAMOutput::soft) {
                // Only support two clipping methods: soft or subread.
                clipping = SAMOutput::subread;
            }
            if (queryFileType != PBBAM and queryFileType != PBDATASET and not enableHiddenPaths) {
                // bax|fasta|fastq -> bam paths are turned off by default
                cout << "ERROR, could not output alignments in BAM unless input reads are in PacBio BAM or DATASET files." << endl;
                exit(1);
            }
            if (outFileName == "") {
                cout << "ERROR, BAM output file must be specified." << endl;
                exit(1);
            }
            if (outputByThread) {
                cout << "ERROR, could not output alignments by threads in BAM format." << endl;
                exit(1);
            }
#endif
        }

        if (limsAlign != 0) {
            mapSubreadsSeparately = false;
            forwardOnly = true;
        }

        if (holeNumberRangesStr.size() > 0) {
            if (not holeNumberRanges.setRanges(holeNumberRangesStr)) {
                cout << "ERROR, could not parse hole number ranges: "
                    << holeNumberRangesStr << "." << endl;
                exit(1);
            }
        }

        if (printSAMQV) {
            if (samQV.size() == 0) {
                samQVList.SetDefaultQV();
            }
            else {
                samQVList.UseQV(samQV);
            }
        }

        if (minRawSubreadScore > 1000) {
            cout << "ERROR, minimum raw subread score should be less than 1000." << endl;
            exit(1);
        }
        if (minRawSubreadScore != -1 and byAdapter) {
            cout << "ERROR, minRawSubreadScore and byAdapter should not be used together." << endl;
            exit(1);
        }
        // Determine query read type
        queryReadType = DetermineQueryReadType();
        // Pass verbosity
        anchorParameters.verbosity = verbosity;

        // Set filter criteria and hit policy
        ResetFilterAndHit();
    }
    void ResetFilterAndHit(void) {
        filterCriteria = FilterCriteria(minAlnLength, minPctSimilarity,
                                        minPctAccuracy, true,
                                        Score(static_cast<float>(maxScore), ScoreSign::NEGATIVE));
        hitPolicy = HitPolicy(hitPolicyStr, ScoreSign::NEGATIVE);
    }

    ReadType::ReadTypeEnum DetermineQueryReadType() {
        if (useCcsOnly or queryFileType == HDFCCSONLY) {
            return ReadType::CCS;
        }
        if (queryFileType == PBBAM) {
            // Read type in BAM may be CCS, SUBREAD, HQREGION or POLYMERASE.
            // Determine it later.
            return ReadType::UNKNOWN;
        }
        if (mapSubreadsSeparately) {
            return ReadType::SUBREAD;
        } else {
            if (useHQRegionTable) {
                return ReadType::HQREGION;
            } else {
                return ReadType::POLYMERASE;
            }
        }
    }

    void SetEmulateNucmer() {
        anchorParameters.stopMappingOnceUnique = true;
        anchorParameters.advanceExactMatches   = 30;
        anchorParameters.maxAnchorsPerPosition = 1;
        sdpBypassThreshold                     = 0.75;
        sdpTupleSize                           = 15;
        anchorParameters.minMatchLength        = 30;
        useGuidedAlign                         = true;
        refineAlignments                       = false;
    }

    void SetForSensitivity() {
        advanceExactMatches = 0;
        anchorParameters.numBranches = 1;
        anchorParameters.maxAnchorsPerPosition = 10000;
    }
};
