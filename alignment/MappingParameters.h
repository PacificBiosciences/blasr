#ifndef MAPPING_PARAMETERS_H_
#define MAPPING_PARAMETERS_H_
#include <vector>

#include "tuples/TupleMetrics.hpp"
#include "anchoring/AnchorParameters.hpp"
#include "qvs/QualityValue.hpp"
#include "format/SAMPrinter.hpp"
#include "algorithms/alignment/AlignmentFormats.h" // not ported
#include "utils/RangeUtils.hpp"

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
    int refineAlign;
    int useScoreCutoff;
    int maxScore;
    int argi;
    int nProc;
    int globalChainType;
    int readIndex;
    SAMOutput::Clipping clipping;
    string clippingString;
    QVScale qvScaleType;
    vector<string> readsFileNames;
    vector<string> regionTableFileNames;
    vector<string> ccsFofnFileNames;
    string tupleListName;
    string posTableName;
    string outFileName;
    string genomeFileName;
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
    int  minAvgQual;
    bool overlap;
    bool advanceHalf;
    int advanceExactMatches;
    float approximateMaxInsertionRate;
    float minPctIdentity;
    float maxPctIdentity;
    bool refineAlignments;
    VectorIndex nCandidates;
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
    bool unrollCcs;
    bool useCcs;
    bool useAllSubreadsInCcs;
    bool useCcsOnly;
    bool detailedSDPAlignment, nouseDetailedSDPAlignment;
    int  chunkSize;
    int  subreadMapType;
    int  sdpFilterType;
    bool useGuidedAlign;
    int  guidedAlignBandSize;
    int  bandSize;
    int  extendBandSize;
    bool useQVScore;
    int  scoreType;
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
    int   maxReadIndex;
    int   recurseOver;
    bool  forPicard;
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
    void Init() {
        readIndex = -1;
        maxReadIndex = -1;
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
        refineAlign = 1;
        useScoreCutoff = 0;
        maxScore = -200;
        argi = 1;
        nProc = 1;
        assert(readsFileNames.empty());
        tupleListName = "";
        genomeFileName = "";
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
        minAvgQual = 0;
        overlap = false;
        advanceHalf = false;
        refineAlignments = true;
        anchorParameters.advanceExactMatches = advanceExactMatches = 0;
        approximateMaxInsertionRate = 1.30;
        minPctIdentity = 0;
        maxPctIdentity = 100.1;
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
        subreadMapType = 0;
        unrollCcs  = false;
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
        forPicard = false;
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
    }

    MappingParameters() {
        Init();
    }

    void MakeSane(){ 
        //
        // Fix all logical incompatibilities with parameters.
        //

        if (nowarp) {
            warp = false;
        }

        if (nCandidates < nBest) {
            cout << "Warning: resetting nCandidates to nBest " << nBest << endl;
            nCandidates = nBest;
        }


        if (placeRandomly and nBest == 1) {
            cout << "Warning: When attempting to select equivalently scoring reads at random " << endl
                << "the bestn parameter should be greater than one." << endl;
        }

        if (concordant) {
            if (useCcs) {
                concordant = false;
            } else {
                useRegionTable   = true;
                useHQRegionTable = true;
            }
        }

        if (sdpFilterType > 1) {
            cout << "Warning: using new filter method for SDP alignments.  The parameter is " << endl
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
        if (maxScore != 0) {
            useScoreCutoff = 1;
        }
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
            cout << "ERROR: maxLCPLength is less than minLCPLength, which will result in no hits." << endl;
        }
        if (subsample < 1 and stride > 1) {
            cout << "ERROR, subsample and stride must be used independently." << endl;
            exit(1);
        }

        if (subreadMapType < 0 or subreadMapType > 1) {
            cout << "Error, subreadImplType must be 0 or 1" << endl;
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
            forPicard = true;
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
        else if (clippingString != "") {
            cout << "ERROR, clipping should either be soft, hard, or none." << endl;
            exit(1);
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


#endif
