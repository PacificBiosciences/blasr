#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

#include "HDFPlsReader.hpp"
#include "amos/AfgBasWriter.hpp"
#include "HDFRegionTableReader.hpp"
#include "reads/RegionTable.hpp"
#include "reads/ReadInterval.hpp"
#include "files/ReaderAgglomerate.hpp"
#include "utils/FileOfFileNames.hpp"
#include "utils/RegionUtils.hpp"
#include "SMRTSequence.hpp"
#include "utils.hpp"

using namespace std;
void PrintUsage() {
        cout << "usage: toAfg input.filetype output.filetype" << endl
                 << "                                         [-minSubreadLength l] " << endl
                 << "                                         [-regionTable regions_file] " << endl
                 << "                                         [-noSplitSubreads]" << endl
                 << "                                         [-useccsdenovo]" << endl
                 << "                                         [-uniformQV QV]" << endl 
        << "Print reads stored in a file (pls|fasta|fastq) as an afg." << endl;
}
int main(int argc, char* argv[]) {

    string inputFileName, outputFileName;

    if (argc < 2) {
        PrintUsage();
        exit(1);
    }
    vector<string> inputFileNames;
    inputFileName = argv[1];
    outputFileName = argv[2];
    int argi = 3;
    RegionTable regionTable;
    string regionsFOFNName = "";
    vector<string> regionFileNames;
    bool splitSubreads = true;
    bool useCCS = false;
    bool useUniformQV = false;
    int uniformQV = 7;
    int minSubreadLength = 1;
    while (argi < argc) {
        if (strcmp(argv[argi], "-regionTable") == 0) {
            regionsFOFNName = argv[++argi];
        }
        else if (strcmp(argv[argi], "-noSplitSubreads") == 0) {
            splitSubreads = false;
        }
        else if (strcmp(argv[argi], "-minSubreadLength") == 0) {
            minSubreadLength = atoi(argv[++argi]);
        }
        else if (strcmp(argv[argi], "-useccsdenovo") == 0) {
            useCCS = true;
        }
        else if (strcmp(argv[argi], "-uniformQV") == 0) {
            useUniformQV = true;
            uniformQV = atoi(argv[++argi]);
        }
        else {
            PrintUsage();
            cout << "ERROR! Option " << argv[argi] << " is not supported." << endl;
        }
        argi++;
    }
         
    if (FileOfFileNames::IsFOFN(inputFileName)) {
        FileOfFileNames::FOFNToList(inputFileName, inputFileNames);
    }
    else {
        inputFileNames.push_back(inputFileName);
    }
    if (regionsFOFNName == "") {
        regionFileNames = inputFileNames;
    }
    else {
        if (FileOfFileNames::IsFOFN(regionsFOFNName)) {
            FileOfFileNames::FOFNToList(regionsFOFNName, regionFileNames);
        }
        else {
            regionFileNames.push_back(regionsFOFNName);
        }
    }

    ofstream fastaOut;
    CrucialOpen(outputFileName, fastaOut);
    int plsFileIndex;
    HDFRegionTableReader hdfRegionReader;
    AfgBasWriter afgWriter;
    if (useUniformQV){
        afgWriter.SetDefaultQuality(uniformQV);
    }

    afgWriter.Initialize(outputFileName);

    for (plsFileIndex = 0; plsFileIndex < inputFileNames.size(); plsFileIndex++) {
        if (splitSubreads) {
            hdfRegionReader.Initialize(regionFileNames[plsFileIndex]);
            hdfRegionReader.ReadTable(regionTable);
            regionTable.SortTableByHoleNumber();
        }
        
        ReaderAgglomerate reader;
        // reader.SkipReadQuality(); // should have been taken care of by *Filter modules
        if (useCCS){
            reader.UseCCS();
        } else {
            reader.IgnoreCCS();
        }
        reader.Initialize(inputFileNames[plsFileIndex]);
        CCSSequence seq; 
        int seqIndex = 0;
        int numRecords = 0;
        vector<ReadInterval> subreadIntervals;
        while (reader.GetNext(seq)){ 
            ++seqIndex;

            if (useUniformQV && seq.qual.data != NULL){
                for (int qvIndex = 0; qvIndex < seq.length; qvIndex++){
                    seq.qual[qvIndex] = uniformQV;
                }
            }

            if (splitSubreads == false) {
                if (seq.length >= minSubreadLength) {
                    afgWriter.Write(seq);
                }
                seq.Free();
                continue;
            }


            DNALength hqReadStart, hqReadEnd;
            int score;
            GetReadTrimCoordinates(seq, seq.zmwData, regionTable, hqReadStart, hqReadEnd, score);
            subreadIntervals.clear(); // clear old, new intervals are appended.
            CollectSubreadIntervals(seq,&regionTable, subreadIntervals);

            if (seq.length == 0 and subreadIntervals.size() > 0) {
                cout << "WARNING! A high quality interval region exists for a read of length 0." <<endl;
                cout << "  The offending ZMW number is " << seq.HoleNumber() << endl;
                seq.Free();
                continue;
            }


            for (int intvIndex = 0; intvIndex < subreadIntervals.size(); intvIndex++) {
                SMRTSequence subreadSequence;
                
                int subreadStart = subreadIntervals[intvIndex].start > hqReadStart ? 
                                   subreadIntervals[intvIndex].start : hqReadStart;
                int subreadEnd   = subreadIntervals[intvIndex].end < hqReadEnd ?
                                   subreadIntervals[intvIndex].end : hqReadEnd;
                int subreadLength = subreadEnd - subreadStart;

                if (subreadLength < minSubreadLength) continue;

                subreadSequence.SubreadStart(subreadStart);
                subreadSequence.SubreadEnd  (subreadEnd);
                subreadSequence.ReferenceSubstring(seq, subreadStart, subreadLength);      

            
                stringstream titleStream;
                titleStream << seq.title << "/" << subreadIntervals[intvIndex].start 
                                         << "_" << subreadIntervals[intvIndex].end;
                subreadSequence.CopyTitle(titleStream.str());
                afgWriter.Write(subreadSequence);
            }
            seq.Free();
        }
        reader.Close();
        hdfRegionReader.Close();
    }
}
