#pragma once
#include <iostream>
#include <boost/filesystem.hpp>
#include <boost/dll.hpp>
#include <algorithm>
#include <cctype>
#include "sra.h"
#include "transcript.h"
#include "ini_parse.h"
#include "sra_toolkit.h"
#include "trinity_wrap.h"
#include "print_info.h"

std::vector<transcript> get_transcript(std::vector<SRA> sras);

bool stringToBool(std::string boolStr);

void makeGroupCheckpoint(std::string cpDir, std::string prefix);

bool groupCheckpointExists(std::string cpDir, std::string prefix);

void makeTransInfoFile(const std::vector<SRA> & sras, std::string transInfoFileStr);

void run_trinity_bulk(std::map<std::string, std::vector<SRA>> sraGroups,
                      std::string threads, std::string ram_gb,
                      bool dispOutput, bool retainInterFiles,
                      std::string logFile, const INI_MAP & cfgIni);
