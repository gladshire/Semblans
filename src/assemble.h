#pragma once
#include <iostream>
#include <sstream>
#include <boost/filesystem.hpp>
#include <boost/dll.hpp>
#include <algorithm>
#include <cctype>
#include <iomanip>
#include <thread>
#include <chrono>
#include <atomic>
#include <set>
#include <signal.h>
#include "sra.h"
#include "transcript.h"
#include "ini_parse.h"
#include "sra_toolkit.h"
#include "trinity_wrap.h"
#include "star_wrap.h"
#include "print_info.h"
#include "seq.h"
#include "seq_hash.h"

void updateHeaders(std::string fastaFilePath, std::string newPrefix,
                   uintmax_t ram_b);

std::vector<transcript> get_transcript(std::vector<SRA> sras);

bool stringToBool(std::string boolStr);

void makeGroupCheckpoint(std::string cpDir, std::string prefix);

bool groupCheckpointExists(std::string cpDir, std::string prefix);

void makeTransInfoFile(const std::vector<std::pair<std::string, std::string>> & sraRuns, std::string transInfoFileStr);

std::vector<std::string> run_trinity_bulk(std::map<std::string, std::vector<SRA>> sraGroups,
                                          std::string threads, std::string ram_gb,
                                          bool assembSeqsInterest, bool assembSeqsNoInterest,
                                          bool assembAllSeqs,
                                          bool dispOutput, bool retainInterFiles,
                                          std::string logFile, const INI_MAP & cfgIni);
