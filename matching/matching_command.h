#ifndef MATCHING_COMMAND_H
#define MATCHING_COMMAND_H

#include "utility/command_parser.h"
#include <map>
#include <iostream>

enum OptionKeyword {
    QueryGraphFile = 1,     // -q, The query graph file path, compulsive parameter
    DataGraphFile = 2,      // -d, The data graph file path, compulsive parameter
    MaxOutputEmbeddingNum = 3, // -num, The maximum output embedding num
    CSRFilePath = 4                    // -csr, The input csr file path
};

class MatchingCommand : public CommandParser{
private:
	//the pair <actual meaning, command abbreviation>
    std::map<OptionKeyword, std::string> options_key;
    //the pair <actual meaning, command parameter>
    std::map<OptionKeyword, std::string> options_value;

private:
    void processOptions();

public:
    MatchingCommand(int argc, char **argv);

    std::string getDataGraphFilePath() {
        return options_value[OptionKeyword::DataGraphFile];
    }

    std::string getQueryGraphFilePath() {
        return options_value[OptionKeyword::QueryGraphFile];
    }

    std::string getMaximumEmbeddingNum() {
        return options_value[OptionKeyword::MaxOutputEmbeddingNum] == "" ? "MAX" : options_value[OptionKeyword::MaxOutputEmbeddingNum];
    }

    std::string getCSRFilePath() {
        return options_value[OptionKeyword::CSRFilePath] == "" ? "" : options_value[OptionKeyword::CSRFilePath];
    }
};

#endif