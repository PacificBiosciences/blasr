// Author: Derek Barnett

#ifndef CONVERTERBASE_H
#define CONVERTERBASE_H

#include "IConverter.h"
#include "Settings.h"
#include "HDFBasReader.hpp"
#include <pbbam/BamFile.h>
#include <pbbam/BamHeader.h>
#include <pbbam/BamWriter.h>
#include <pbbam/PbiFile.h>
#include <pbbam/ReadGroupInfo.h>
#include <pbbam/Tag.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <libgen.h>
#include <cstdlib>
#include <climits>
#include <map>
#include <set>
#include <string>
#include <vector>

namespace PacBio {
namespace BAM {

class BamRecordImpl;

} // namespace BAM
} // namespace PacBio

template<typename RecordType = SMRTSequence, typename HdfReader = HDFBasReader>
class ConverterBase : public IConverter
{
public:
    ~ConverterBase(void);

public:
    virtual bool Run(void) final;

protected:
    ConverterBase(Settings& settings);

    virtual bool ConvertFile(HdfReader* reader,
                             PacBio::BAM::BamWriter* writer) =0;

    virtual bool ConvertFile(HdfReader* reader,
                             PacBio::BAM::BamWriter* writer,
                             PacBio::BAM::BamWriter* scrapsWriter) =0;

    virtual bool ConvertRecord(const RecordType& smrtRecord,
                               const int start,
                               const int end,
                               const std::string& rgId,
                               PacBio::BAM::BamRecordImpl* bamRecord);

    virtual bool WriteRecord(const RecordType& smrtRecord,
                             const int recordStart,
                             const int recordEnd,
                             const std::string& readGroupId,
                             PacBio::BAM::BamWriter* writer);

    virtual bool WriteFilteredRecord(const RecordType& smrtRecord,
                                     const int recordStart,
                                     const int recordEnd,
                                     const std::string& readGroupId,
                                     PacBio::BAM::BamWriter* writer);

    virtual bool WriteFilteredRecord(const RecordType& smrtRecord,
                                     const int recordStart,
                                     const int recordEnd,
                                     const std::string& readGroupId,
                                     const uint8_t contextFlags,
                                     PacBio::BAM::BamWriter* writer);

    virtual bool WriteLowQualityRecord(const RecordType& smrtRecord,
                                       const int recordStart,
                                       const int recordEnd,
                                       const std::string& readGroupId,
                                       PacBio::BAM::BamWriter* writer);

    virtual bool WriteAdapterRecord(const RecordType& smrtRecord,
                                    const int recordStart,
                                    const int recordEnd,
                                    const std::string& readGroupId,
                                    PacBio::BAM::BamWriter* writer);

    virtual bool WriteSubreadRecord(const RecordType& smrtRecord,
                                    const int recordStart,
                                    const int recordEnd,
                                    const std::string& readGroupId,
                                    const uint8_t contextFlags,
                                    PacBio::BAM::BamWriter* writer);

    virtual void SetSequenceAndQualities(PacBio::BAM::BamRecordImpl* bamRecord,
                                         const RecordType& smrtRecord,
                                         const int start,
                                         const int length);

    virtual void AddRecordName(PacBio::BAM::BamRecordImpl* bamRecord,
                               const UInt holeNumber,
                               const int start,
                               const int end);

    virtual void AddModeTags(PacBio::BAM::TagCollection* tags,
                             const RecordType& smrtRecord,
                             const int start,
                             const int end);

    virtual HdfReader* InitHdfReader(void);
    virtual void InitReadScores(HdfReader* reader) final;

    virtual bool IsSequencingZmw(const RecordType& record) const final;

    virtual bool LoadChemistryFromMetadataXML(const std::string& baxFn,
                                              const std::string& movieName) final; 

    virtual std::string HeaderReadType(void) const =0;
    virtual std::string ScrapsReadType(void) const =0;
    virtual std::string OutputFileSuffix(void) const =0;
    virtual std::string ScrapsFileSuffix(void) const =0;

    // Settings variable accessors
    virtual std::string MovieName(void);
    virtual std::string ReadGroupId(void);
    virtual std::string ScrapsReadGroupId(void);

protected:
    std::vector<HdfReader*> readers_;
    std::map<HdfReader*, std::string> filenameForReader_;

    std::vector<float> readScores_;
    std::map<UInt, size_t> indexForHoleNumber_; // helper table for read scores (holenumber -> vector index)

    // re-used containers
    PacBio::BAM::BamRecordImpl bamRecord_;
    std::string recordSequence_;
    PacBio::BAM::QualityValues recordDeletionQVs_;
    PacBio::BAM::QualityValues recordInsertionQVs_;
    PacBio::BAM::QualityValues recordMergeQVs_;
    PacBio::BAM::QualityValues recordSubstitutionQVs_;
    std::string recordDeletionTags_;
    std::string recordSubstitutionTags_;
    std::vector<uint16_t> recordRawIPDs_;
    std::vector<uint8_t> recordEncodedIPDs_;
    std::vector<uint16_t> recordRawPulseWidths_;
    std::vector<uint8_t> recordEncodedPulseWidths_;

    // IPD downsampling
    std::vector<uint16_t> framepoints_;
    std::vector<uint8_t> frameToCode_;
    uint16_t maxFramepoint_;

    // store tags
    //
    //     i: signed32, I: unsigned32, C: unsigned8, A: ASCII, Z: string
    //
    // qs:i - 0-based start of query in the polymerase read
    // qe:i - 0-based end of query in the polymerase read
    // zm:i - ZMW hole number
    // np:i - NumPasses (1 for subreads, variable for CCS)
    // rq:i - float in [0.0,1.0] encoding expected accuracy
    // dq:Z - DeletionQV
    // dt:Z - DeletionTag
    // iq:Z - InsertionQV
    // mq:Z - MergeQV
    // sq:Z - SubstitutionQV
    // st:Z - SubstitutionTag
    // ip:B,C *or* B,S - IPD (frames: 8-bit (lossy) or 16-bit (full)
    // pw:B,C *or* B,S - PulseWidth (frames: 8-bit (lossy) or 16-bit (full)
    // sc:A - Scrap-type
    // sz:A - ZMW classification
    //
    // RG:Z - standard SAM/BAM RG tag, contains the corresponding @RG:ID
    //
    static const std::string Tag_zm;
    static const std::string Tag_rq;
    static const std::string Tag_cx;
    static const std::string Tag_sn;
    static const std::string Tag_dq;
    static const std::string Tag_dt;
    static const std::string Tag_iq;
    static const std::string Tag_mq;
    static const std::string Tag_sq;
    static const std::string Tag_st;
    static const std::string Tag_ip;
    static const std::string Tag_pw;
    static const std::string Tag_sc;
    static const std::string Tag_sz;
    static const std::string Tag_RG;

    // store re-used tag values
    //
    // Adapter Tag    = Tag('A', ASCII_CHAR)
    // LowQuality Tag = Tag('L', ASCII_CHAR)
    // Filtered Tag   = Tag('F', ASCII_CHAR)
    // NormalZMW Tag  = Tag('N', ASCII_CHAR)
    //
    static const PacBio::BAM::Tag lowQualityTag_;
    static const PacBio::BAM::Tag adapterTag_;
    static const PacBio::BAM::Tag filteredTag_;
    static const PacBio::BAM::Tag normalZmwTag_;
};

// Static Tag-name initializers
template<typename RecordType, typename HdfReader>
const std::string ConverterBase<RecordType, HdfReader>::Tag_zm = "zm";
template<typename RecordType, typename HdfReader>
const std::string ConverterBase<RecordType, HdfReader>::Tag_rq = "rq";
template<typename RecordType, typename HdfReader>
const std::string ConverterBase<RecordType, HdfReader>::Tag_cx = "cx";
template<typename RecordType, typename HdfReader>
const std::string ConverterBase<RecordType, HdfReader>::Tag_sn = "sn";
template<typename RecordType, typename HdfReader>
const std::string ConverterBase<RecordType, HdfReader>::Tag_dq = "dq";
template<typename RecordType, typename HdfReader>
const std::string ConverterBase<RecordType, HdfReader>::Tag_dt = "dt";
template<typename RecordType, typename HdfReader>
const std::string ConverterBase<RecordType, HdfReader>::Tag_iq = "iq";
template<typename RecordType, typename HdfReader>
const std::string ConverterBase<RecordType, HdfReader>::Tag_mq = "mq";
template<typename RecordType, typename HdfReader>
const std::string ConverterBase<RecordType, HdfReader>::Tag_sq = "sq";
template<typename RecordType, typename HdfReader>
const std::string ConverterBase<RecordType, HdfReader>::Tag_st = "st";
template<typename RecordType, typename HdfReader>
const std::string ConverterBase<RecordType, HdfReader>::Tag_ip = "ip";
template<typename RecordType, typename HdfReader>
const std::string ConverterBase<RecordType, HdfReader>::Tag_pw = "pw";
template<typename RecordType, typename HdfReader>
const std::string ConverterBase<RecordType, HdfReader>::Tag_sc = "sc";
template<typename RecordType, typename HdfReader>
const std::string ConverterBase<RecordType, HdfReader>::Tag_sz = "sz";
template<typename RecordType, typename HdfReader>
const std::string ConverterBase<RecordType, HdfReader>::Tag_RG = "RG";

// Static Tag-Value initializers
template<typename RecordType, typename HdfReader>
const PacBio::BAM::Tag ConverterBase<RecordType, HdfReader>::lowQualityTag_ =
        PacBio::BAM::Tag('L', PacBio::BAM::TagModifier::ASCII_CHAR);

template<typename RecordType, typename HdfReader>
const PacBio::BAM::Tag ConverterBase<RecordType, HdfReader>::adapterTag_ =
        PacBio::BAM::Tag('A', PacBio::BAM::TagModifier::ASCII_CHAR);

template<typename RecordType, typename HdfReader>
const PacBio::BAM::Tag ConverterBase<RecordType, HdfReader>::filteredTag_ =
        PacBio::BAM::Tag('F', PacBio::BAM::TagModifier::ASCII_CHAR);

template<typename RecordType, typename HdfReader>
const PacBio::BAM::Tag ConverterBase<RecordType, HdfReader>::normalZmwTag_ =
        PacBio::BAM::Tag('N', PacBio::BAM::TagModifier::ASCII_CHAR);

// Constructor
template<typename RecordType, typename HdfReader>
ConverterBase<RecordType, HdfReader>::ConverterBase(Settings& settings)
    : IConverter(settings)
{ }

// Destructor
template<typename RecordType, typename HdfReader>
ConverterBase<RecordType, HdfReader>::~ConverterBase(void)
{
    auto end  = readers_.end();
    for (auto iter = readers_.begin(); iter != end; ++iter) {
        HdfReader* r = (*iter);
        if (r) {
            r->Close();
            delete r;
            r = 0;
        }
    }
    readers_.clear();
}

template<typename RecordType, typename HdfReader>
std::string ConverterBase<RecordType, HdfReader>::MovieName(void)
{
    return settings_.movieName;
}

template<typename RecordType, typename HdfReader>
std::string ConverterBase<RecordType, HdfReader>::ReadGroupId(void)
{
    return settings_.readGroupId;
}

template<typename RecordType, typename HdfReader>
std::string ConverterBase<RecordType, HdfReader>::ScrapsReadGroupId(void)
{
    return settings_.scrapsReadGroupId;
}

template<typename RecordType, typename HdfReader>
bool ConverterBase<RecordType, HdfReader>::ConvertRecord(
        const RecordType& smrtRead,
        const int subreadStart,
        const int subreadEnd,
        const std::string& rgId,
        PacBio::BAM::BamRecordImpl* bamRecord)
{
    using namespace PacBio;
    using namespace PacBio::BAM;
    using namespace std;

    // sanity check
    assert(bamRecord);

    const UInt holeNumber   = smrtRead.zmwData.holeNumber;
    const DNALength length = subreadEnd - subreadStart;

    AddRecordName(bamRecord, holeNumber, subreadStart, subreadEnd);

    // store sequence
    // NOTE - qualities are empty (per PacBio BAM spec)
    SetSequenceAndQualities(bamRecord, smrtRead, subreadStart, length);

    // check settings/existence of *QV/*Tag data
    if (settings_.usingDeletionQV && smrtRead.deletionQV.Empty())
    {
        AddErrorMessage("DeletionQV requested but unavailable");
        return false;
    }

    if (settings_.usingInsertionQV && smrtRead.insertionQV.Empty())
    {
        AddErrorMessage("InsertionQV requested but unavailable");
        return false;
    }

    if (settings_.usingMergeQV && smrtRead.mergeQV.Empty())
    {
        AddErrorMessage("MergeQV requested but unavailable");
        return false;
    }

    if (settings_.usingSubstitutionQV && smrtRead.substitutionQV.Empty())
    {
        AddErrorMessage("SubstitutionQV requested but unavailable");
        return false;
    }

    if (settings_.usingDeletionTag && smrtRead.deletionTag == nullptr)
    {
        AddErrorMessage("DeletionTag requested but unavailable");
        return false;
    }

    if (settings_.usingSubstitutionTag && smrtRead.substitutionTag == nullptr)
    {
        AddErrorMessage("SubstitutionTag requested but unavailable");
        return false;
    }

    if (settings_.usingIPD && smrtRead.preBaseFrames == nullptr)
    {
        AddErrorMessage("IPD requested but unavailable");
        return false;
    }

    if (settings_.usingPulseWidth && smrtRead.widthInFrames == nullptr)
    {
        AddErrorMessage("PulseWidth requested but unavailable");
        return false;
    }

    // fetch *QV/*Tag data
    if (settings_.usingDeletionQV) {
        recordDeletionQVs_.assign((uint8_t*)smrtRead.deletionQV.data + subreadStart,
                                  (uint8_t*)smrtRead.deletionQV.data + subreadStart + length);
    }
    if (settings_.usingInsertionQV) {
        recordInsertionQVs_.assign((uint8_t*)smrtRead.insertionQV.data + subreadStart,
                                   (uint8_t*)smrtRead.insertionQV.data + subreadStart + length);
    }
    if (settings_.usingMergeQV) {
        recordMergeQVs_.assign((uint8_t*)smrtRead.mergeQV.data + subreadStart,
                               (uint8_t*)smrtRead.mergeQV.data + subreadStart + length);
    }
    if (settings_.usingSubstitutionQV) {
        recordSubstitutionQVs_.assign((uint8_t*)smrtRead.substitutionQV.data + subreadStart,
                                      (uint8_t*)smrtRead.substitutionQV.data + subreadStart + length);
    }
    if (settings_.usingDeletionTag) {
        recordDeletionTags_.assign((char*)smrtRead.deletionTag + subreadStart,
                                   (char*)smrtRead.deletionTag + subreadStart + length);
    }
    if (settings_.usingSubstitutionTag) {
        recordSubstitutionTags_.assign((char*)smrtRead.substitutionTag + subreadStart,
                                       (char*)smrtRead.substitutionTag + subreadStart + length);
    }

    // fetch IPDs, then maybe encode
    if (settings_.usingIPD) {
        recordRawIPDs_.assign((uint16_t*)smrtRead.preBaseFrames + subreadStart,
                              (uint16_t*)smrtRead.preBaseFrames + subreadStart + length);

        // if not using full data, encode
        if (!settings_.losslessFrames)
            recordEncodedIPDs_ = std::move(Frames::Encode(recordRawIPDs_));
    }

    // fetch PulseWidths, then maybe encode
    if (settings_.usingPulseWidth) {
        recordRawPulseWidths_.assign((uint16_t*)smrtRead.widthInFrames + subreadStart,
                                     (uint16_t*)smrtRead.widthInFrames + subreadStart + length);

        // if not using full data, encode
        if (!settings_.losslessFrames)
            recordEncodedPulseWidths_ = std::move(Frames::Encode(recordRawPulseWidths_));
    }

    TagCollection tags;
    tags[Tag_RG] = rgId;
    tags[Tag_zm] = static_cast<int32_t>(holeNumber);

    // HQRegionSNR, TODO: should I do this in AddModeTags?
    if (HeaderReadType() != "CCS")
    {
        // Stored as 'ACGT' in BAM, no fixed order in SMRTSequence
        vector<float> hqSnr = { smrtRead.HQRegionSnr('A'),
                                smrtRead.HQRegionSnr('C'),   
                                smrtRead.HQRegionSnr('G'),   
                                smrtRead.HQRegionSnr('T')};
        tags[Tag_sn] = hqSnr;
    }

    AddModeTags(&tags, smrtRead, subreadStart, subreadEnd);

    if (!readScores_.empty())
        tags[Tag_rq] = static_cast<float>(readScores_.at(indexForHoleNumber_[holeNumber]));
    else
        tags[Tag_rq] = static_cast<float>(0.0f);

    if (settings_.usingDeletionQV)      tags[Tag_dq] = recordDeletionQVs_.Fastq();
    if (settings_.usingDeletionTag)     tags[Tag_dt] = recordDeletionTags_;
    if (settings_.usingInsertionQV)     tags[Tag_iq] = recordInsertionQVs_.Fastq();
    if (settings_.usingMergeQV)         tags[Tag_mq] = recordMergeQVs_.Fastq();
    if (settings_.usingSubstitutionQV)  tags[Tag_sq] = recordSubstitutionQVs_.Fastq();
    if (settings_.usingSubstitutionTag) tags[Tag_st] = recordSubstitutionTags_;

    if (settings_.usingIPD) {
        if (settings_.losslessFrames)
            tags[Tag_ip] = recordRawIPDs_;
        else
            tags[Tag_ip] = recordEncodedIPDs_;

    }

    if (settings_.usingPulseWidth) {
        if (settings_.losslessFrames)
            tags[Tag_pw] = recordRawPulseWidths_;
        else
            tags[Tag_pw] = recordEncodedPulseWidths_;
    }

    bamRecord->Tags(tags);

    // if we get here, everything should be OK
    return true;
}

template<typename RecordType, typename HdfReader>
bool ConverterBase<RecordType, HdfReader>::WriteRecord(const RecordType& smrtRecord,
                                                       const int recordStart,
                                                       const int recordEnd,
                                                       const std::string& readGroupId,
                                                       PacBio::BAM::BamWriter* writer)
{
    // attempt convert BAX to BAM
    if (!ConvertRecord(smrtRecord,
                       recordStart,
                       recordEnd,
                       readGroupId,
                       &bamRecord_))
    {
        return false;
    }

    // attempt write BAM to file
    try {
        writer->Write(bamRecord_);
    } catch (std::exception&) {
        AddErrorMessage("failed to write BAM record");
        return false;
    }

    return true;
}

template<typename RecordType, typename HdfReader>
bool ConverterBase<RecordType, HdfReader>::WriteFilteredRecord(const RecordType& smrtRecord,
                                                               const int recordStart,
                                                               const int recordEnd,
                                                               const std::string& readGroupId,
                                                               PacBio::BAM::BamWriter* writer)
{
    // attempt convert BAX to BAM
    if (!ConvertRecord(smrtRecord,
                       recordStart,
                       recordEnd,
                       readGroupId,
                       &bamRecord_))
    {
        return false;
    }

    // add scrap tags
    if (!bamRecord_.AddTag(Tag_sz, normalZmwTag_))
    {
        AddErrorMessage("failed to add scrap's zmw classification tag");
        return false;
    }
    if (!bamRecord_.AddTag(Tag_sc, filteredTag_))
    {
        AddErrorMessage("failed to add scrap's filtered tag");
        return false;
    }

    // attempt write BAM to file
    try {
        writer->Write(bamRecord_);
    } catch (std::exception&) {
        AddErrorMessage("failed to write BAM record");
        return false;
    }

    return true;
}

template<typename RecordType, typename HdfReader>
bool ConverterBase<RecordType, HdfReader>::WriteFilteredRecord(const RecordType& smrtRecord,
                                                               const int recordStart,
                                                               const int recordEnd,
                                                               const std::string& readGroupId,
                                                               const uint8_t contextFlags,
                                                               PacBio::BAM::BamWriter* writer)
{
    // attempt convert BAX to BAM
    if (!ConvertRecord(smrtRecord,
                       recordStart,
                       recordEnd,
                       readGroupId,
                       &bamRecord_))
    {
        return false;
    }

    // add scrap tags
    if (!bamRecord_.AddTag(Tag_sz, normalZmwTag_))
    {
        AddErrorMessage("failed to add scrap's zmw classification tag");
        return false;
    }
    if (!bamRecord_.AddTag(Tag_sc, filteredTag_))
    {
        AddErrorMessage("failed to add scrap's filtered tag");
        return false;
    }

    // add context tag
    if (!bamRecord_.AddTag(Tag_cx, contextFlags))
    {
        AddErrorMessage("failed to add context flag tag");
        return false;
    }

    // attempt write BAM to file
    try {
        writer->Write(bamRecord_);
    } catch (std::exception&) {
        AddErrorMessage("failed to write BAM record");
        return false;
    }

    return true;
}

template<typename RecordType, typename HdfReader>
bool ConverterBase<RecordType, HdfReader>::WriteLowQualityRecord(const RecordType& smrtRecord,
                                                                 const int recordStart,
                                                                 const int recordEnd,
                                                                 const std::string& readGroupId,
                                                                 PacBio::BAM::BamWriter* writer)
{
    // attempt convert BAX to BAM
    if (!ConvertRecord(smrtRecord,
                       recordStart,
                       recordEnd,
                       readGroupId,
                       &bamRecord_))
    {
        return false;
    }

    // add scrap tags
    if (!bamRecord_.AddTag(Tag_sz, normalZmwTag_))
    {
        AddErrorMessage("failed to add scrap's zmw classification tag");
        return false;
    }
    if (!bamRecord_.AddTag(Tag_sc, lowQualityTag_))
    {
        AddErrorMessage("failed to add scrap's low-quality region tag");
        return false;
    }

    // attempt write BAM to file
    try {
        writer->Write(bamRecord_);
    } catch (std::exception&) {
        AddErrorMessage("failed to write BAM record");
        return false;
    }

    return true;
}

template<typename RecordType, typename HdfReader>
bool ConverterBase<RecordType, HdfReader>::WriteAdapterRecord(const RecordType& smrtRecord,
                                                              const int recordStart,
                                                              const int recordEnd,
                                                              const std::string& readGroupId,
                                                              PacBio::BAM::BamWriter* writer)
{
    // attempt convert BAX to BAM
    if (!ConvertRecord(smrtRecord,
                       recordStart,
                       recordEnd,
                       readGroupId,
                       &bamRecord_))
    {
        return false;
    }

    // add scrap tags
    if (!bamRecord_.AddTag(Tag_sz, normalZmwTag_))
    {
        AddErrorMessage("failed to add scrap's zmw classification tag");
        return false;
    }
    if (!bamRecord_.AddTag(Tag_sc, adapterTag_))
    {
        AddErrorMessage("failed to add scrap's adapter tag");
        return false;
    }

    // attempt write BAM to file
    try {
        writer->Write(bamRecord_);
    } catch (std::exception&) {
        AddErrorMessage("failed to write BAM record");
        return false;
    }

    return true;
}

template<typename RecordType, typename HdfReader>
bool ConverterBase<RecordType, HdfReader>::WriteSubreadRecord(const RecordType& smrtRecord,
                                                              const int recordStart,
                                                              const int recordEnd,
                                                              const std::string& readGroupId,
                                                              const uint8_t contextFlags,
                                                              PacBio::BAM::BamWriter* writer)
{
    // attempt convert BAX to BAM
    if (!ConvertRecord(smrtRecord,
                       recordStart,
                       recordEnd,
                       readGroupId,
                       &bamRecord_))
    {
        return false;
    }

    // Try to add the additional tag supplied by the caller
    if (!bamRecord_.AddTag(Tag_cx, contextFlags))
    {
        AddErrorMessage("failed to add context flag tag");
        return false;
    }

    // attempt write BAM to file
    try {
        writer->Write(bamRecord_);
    } catch (std::exception&) {
        AddErrorMessage("failed to write BAM record");
        return false;
    }

    return true;
}

template<typename RecordType, typename HdfReader>
void ConverterBase<RecordType, HdfReader>::SetSequenceAndQualities(
        PacBio::BAM::BamRecordImpl* bamRecord,
        const RecordType& smrtRead,
        const int start,
        const int length)
{
    recordSequence_.assign((const char*)smrtRead.seq + start, length);
    bamRecord->SetSequenceAndQualities(recordSequence_);
}

template<typename RecordType, typename HdfReader>
void ConverterBase<RecordType, HdfReader>::AddRecordName(
        PacBio::BAM::BamRecordImpl* bamRecord,
        const UInt holeNumber,
        const int start,
        const int end)
{
    const string name = settings_.movieName + "/"
                      + to_string(holeNumber) + "/"
                      + to_string(start) + "_"
                      + to_string(end);
    bamRecord->Name(name);
}

template<typename RecordType, typename HdfReader>
void ConverterBase<RecordType, HdfReader>::AddModeTags(
        PacBio::BAM::TagCollection* tags,
        const RecordType& smrtRead,
        const int start,
        const int end)
{
    (*tags)["qs"] = start;
    (*tags)["qe"] = end;
    (*tags)["np"] = static_cast<int32_t>(1);
}

template<typename RecordType, typename HdfReader>
HdfReader* ConverterBase<RecordType, HdfReader>::InitHdfReader(void)
{
    HdfReader* reader = new HdfReader;
    reader->IncludeField("Basecall");
    if (HeaderReadType() != "CCS")      reader->IncludeField("HQRegionSNR");
    if (settings_.usingDeletionQV)      reader->IncludeField("DeletionQV");
    if (settings_.usingDeletionTag)     reader->IncludeField("DeletionTag");
    if (settings_.usingInsertionQV)     reader->IncludeField("InsertionQV");
    if (settings_.usingIPD)             reader->IncludeField("PreBaseFrames");
    if (settings_.usingMergeQV)         reader->IncludeField("MergeQV");
    if (settings_.usingPulseWidth)      reader->IncludeField("WidthInFrames");
    if (settings_.usingSubstitutionQV)  reader->IncludeField("SubstitutionQV");
    if (settings_.usingSubstitutionTag) reader->IncludeField("SubstitutionTag");
    return reader;
}

template<typename RecordType, typename HdfReader>
void ConverterBase<RecordType, HdfReader>::InitReadScores(HdfReader* reader)
{
    assert(reader);

    // fetch read scores
    readScores_.clear();
    if (reader->baseCallsGroup.ContainsObject("ZMWMetrics")) {
        HDFGroup zmwMetricsGroup;
        if (zmwMetricsGroup.Initialize(reader->baseCallsGroup.group, "ZMWMetrics")) {
            if (zmwMetricsGroup.ContainsObject("ReadScore")) {
                HDFArray<float> readScoresArray;
                if (readScoresArray.InitializeForReading(zmwMetricsGroup, "ReadScore"))
                    readScoresArray.ReadDataset(readScores_);
            }
        }
    }

    // init holenumber -> index lookup
    indexForHoleNumber_.clear();
    if (!readScores_.empty()) {
        for (size_t i = 0; i < readScores_.size(); ++i) {
            UInt holeNumber;
            reader->zmwReader.GetHoleNumberAt(i, holeNumber);
            indexForHoleNumber_[holeNumber] = i;
        }
    }
}

template<typename RecordType, typename HdfReader>
bool ConverterBase<RecordType, HdfReader>::IsSequencingZmw(const RecordType& record) const
{ return record.zmwData.holeStatus == 0; }

template<typename RecordType, typename HdfReader>
bool ConverterBase<RecordType, HdfReader>::LoadChemistryFromMetadataXML(
        const std::string& baxFn,
        const std::string& movieName)
{
    using boost::property_tree::ptree;

    // get the absolute path of the bax.h5 file and go up 2 directories
    char buf[PATH_MAX + 1];
    char *res = realpath(baxFn.c_str(), buf);

    if (res == nullptr)
        return false;

    // up 1
    res = dirname(res);

    if (res == nullptr)
        return false;

    // up 2
    res = dirname(res);

    if (res == nullptr)
        return false;

    std::string prefix(res);
    std::string path = prefix + '/' + movieName + ".metadata.xml";
    ptree pt;

    try
    {
        read_xml(path, pt);

        bindingKit_        = pt.get<std::string>("Metadata.BindingKit.PartNumber");
        sequencingKit_     = pt.get<std::string>("Metadata.SequencingKit.PartNumber");
        basecallerVersion_ = pt.get<std::string>("Metadata.InstCtrlVer");

        return true;
    }
    catch (...)
    { }

    return false;
}

template<typename RecordType, typename HdfReader>
bool ConverterBase<RecordType, HdfReader>::Run(void)
{
    using namespace PacBio;
    using namespace PacBio::BAM;
    using namespace std;

    set<string> movieNames;

    // initialize input BAX readers
    const auto baxEnd = settings_.inputBaxFilenames.cend();
    for (auto baxIter = settings_.inputBaxFilenames.cbegin(); baxIter != baxEnd; ++baxIter) {
        const string& baxFn = (*baxIter);
        if (baxFn.empty())
            continue;

        HdfReader* reader = InitHdfReader();
        
        // read in mandatory ReadGroupInfo from bax file
        if (reader->Initialize(baxFn) &&
            reader->scanDataReader.fileHasScanData &&
            reader->scanDataReader.initializedRunInfoGroup)
        {
            // FrameRate
            {
                HDFAtom<float> frAtom;
                if (reader->scanDataReader.acqParamsGroup.ContainsAttribute("FrameRate") &&
                    frAtom.Initialize(reader->scanDataReader.acqParamsGroup, "FrameRate"))
                {
                    float localFrameRate;
                    frAtom.Read(localFrameRate);
                    frAtom.dataspace.close();
                    frameRateHz_ = std::to_string(localFrameRate);
                } else {
                    AddErrorMessage("FrameRate is mandatory but unavailable");
                    return false;
                }
            }

            // chemistry triple success flag
            bool success = false;

            // BindingKit
            {
                HDFAtom<std::string> bkAtom;
                if (reader->scanDataReader.runInfoGroup.ContainsAttribute("BindingKit") &&
                    bkAtom.Initialize(reader->scanDataReader.runInfoGroup, "BindingKit"))
                {
                    bkAtom.Read(bindingKit_);
                    bkAtom.dataspace.close();
                } else {
                    goto fallback;
                }
            }

            // SequencingKit
            {
                HDFAtom<std::string> skAtom;
                if (reader->scanDataReader.runInfoGroup.ContainsAttribute("SequencingKit") &&
                    skAtom.Initialize(reader->scanDataReader.runInfoGroup, "SequencingKit"))
                {
                    skAtom.Read(sequencingKit_);
                    skAtom.dataspace.close();
                } else {
                    goto fallback;
                }
            }

            // basecaller ChangeListID
            {
                HDFGroup bcGroup;
                if (reader->pulseDataGroup.ContainsObject("BaseCalls") &&
                    bcGroup.Initialize(reader->pulseDataGroup.group, "BaseCalls"))
                {
                    HDFAtom<std::string> clAtom;
                    if (bcGroup.ContainsAttribute("ChangeListID") &&
                        clAtom.Initialize(bcGroup.group, "ChangeListID"))
                    {
                        clAtom.Read(basecallerVersion_);
                        clAtom.dataspace.close();
                        success = true;
                    }
                    bcGroup.Close();
                }
            }

fallback:
            if (!success && !LoadChemistryFromMetadataXML(baxFn, reader->GetMovieName()))
            {
                AddErrorMessage("BindingKit, SequencingKit, and ChangeListID are mandatory but unavailable");
                return false;
            }

        } else {
            delete reader;
            AddErrorMessage("Failed to properly initialize HDFBasReader");
            return false;
        }

        movieNames.insert(reader->GetMovieName());
        readers_.push_back(reader);
        filenameForReader_[reader] = baxFn;
    }

    if (readers_.empty()) {
        AddErrorMessage("could not open BAX file(s)");
        return false;
    }

    // sanity check that BAX files come from same movie
    if (movieNames.size() != 1) {
        AddErrorMessage("multiple movies detected:");
        for (const auto m : movieNames)
            AddErrorMessage(string("    ")+m);
        return false;
    }
    settings_.movieName = (*movieNames.cbegin());

    // Use the movie name to initialize the ReadGroupId
    settings_.readGroupId = MakeReadGroupId(MovieName(), HeaderReadType());

    // initialize output file(s)
    if (settings_.outputBamPrefix.empty())
        settings_.outputBamPrefix = settings_.movieName;
    settings_.outputBamFilename = settings_.outputBamPrefix + OutputFileSuffix();

    // Separate single-output from dual-output jobs
    if (HeaderReadType() == "SUBREAD" || HeaderReadType() == "HQREGION")
    {
        // setup scram BAM file info
        settings_.scrapsReadGroupId = MakeReadGroupId(MovieName(), ScrapsReadType());
        settings_.scrapsBamFilename = settings_.outputBamPrefix + ScrapsFileSuffix();

        // main conversion of BAX -> BAM records for dual-output jobs
        try {
            BamWriter writer(settings_.outputBamFilename, CreateHeader(HeaderReadType()));
            BamWriter scrapsWriter(settings_.scrapsBamFilename, CreateHeader(ScrapsReadType()));

            for (HdfReader* reader : readers_) {
                assert(reader);
                if (!ConvertFile(reader, &writer, &scrapsWriter))
                    return false;
            }
        } catch (std::exception&) {
            // TODO: get more helpful message here
            AddErrorMessage("failed to convert BAM file");
            return false;
        }

        // make PBI files
        PbiFile::CreateFrom(BamFile{ settings_.outputBamFilename });
        PbiFile::CreateFrom(BamFile{ settings_.scrapsBamFilename });

    } else { 

        assert(settings_.scrapsBamFilename.empty());

        // main conversion of BAX -> BAM records for single-output jobs
        try {
            BamWriter writer(settings_.outputBamFilename, CreateHeader(HeaderReadType()));
          
            for (HdfReader* reader : readers_) {
                assert(reader);
                if (!ConvertFile(reader, &writer))
                    return false;
            }
        } catch (std::exception&) {
            // TODO: get more helpful message here
            AddErrorMessage("failed to convert BAM file");
            return false;
        }

        // make PBI file
        PbiFile::CreateFrom(BamFile{ settings_.outputBamFilename });
    }

    // if we get here, return success
    return true;
}

#endif
