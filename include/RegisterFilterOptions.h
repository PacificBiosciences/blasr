#include "libconfig.h" 
#include "CommandLineParser.hpp"
#include "datastructures/alignment/FilterCriteria.hpp"
#include <sstream>
using namespace std;

/// Register options for filtering alignments.
void RegisterFilterOptions(CommandLineParser & clp, int & minAlnLength,
                           float & minPctSimilarity, float & minPctAccuracy, 
                           string & hitPolicyStr, bool & useScoreCutoff,
                           int & scoreSignInt, int & scoreCutoff) {
    ScoreSign ss = static_cast<ScoreSign>(scoreSignInt);
    Score sc(static_cast<float>(scoreCutoff),  ss);
    FilterCriteria fc(static_cast<DNALength>(minAlnLength), 
                      minPctSimilarity, minPctAccuracy, 
                      useScoreCutoff, sc);

    HitPolicy hp("randombest", ScoreSign::NEGATIVE);

    clp.RegisterIntOption("minAlnLength", &minAlnLength, 
                          fc.MinAlnLengthHelp(),
                          CommandLineParser::PositiveInteger);
    clp.RegisterIntOption("minAlignLength", &minAlnLength, 
                          "Alias of -minAlnLength", 
                          CommandLineParser::PositiveInteger);
    clp.RegisterIntOption("minLength", &minAlnLength, 
                          "Alias of -minAlnLength", 
                          CommandLineParser::PositiveInteger);

    clp.RegisterFloatOption("minPctSimilarity", &minPctSimilarity, 
                            fc.MinPctSimilarityHelp(), 
                            CommandLineParser::PositiveFloat);
    clp.RegisterFloatOption("minPctIdentity", &minPctSimilarity, 
                            "Alias of -minPctSimilarity",
                            CommandLineParser::PositiveFloat);

    clp.RegisterFloatOption("minPctAccuracy", &minPctAccuracy,
                            fc.MinPctAccuracyHelp(),
                            CommandLineParser::PositiveFloat);
    clp.RegisterFloatOption("minAccuracy", &minPctAccuracy,
                            "Alias of -minPctAccuracy",
                            CommandLineParser::PositiveFloat);

    clp.RegisterStringOption("hitPolicy", &hitPolicyStr, hp.Help());
            
    clp.RegisterIntOption("scoreSign", &scoreSignInt, 
                          fc.ScoreSignHelp(),
                          CommandLineParser::Integer);

    clp.RegisterIntOption("scoreCutoff", &scoreCutoff, 
                          fc.ScoreCutoffHelp(),
                          CommandLineParser::Integer);
}
