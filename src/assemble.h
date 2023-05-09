#pragma once
#include <iostream>
#include <boost/filesystem.hpp>
#include <boost/dll.hpp>
#include "sra.h"
#include "transcript.h"
#include "ini_parse.h"
#include "sra_toolkit.h"
#include "trinity_wrap.h"
#include "print_info.h"

std::vector<transcript> get_transcript(std::vector<SRA> sras);

bool stringToBool(std::string boolStr);

std::vector<transcript> run_trinity_bulk(std::vector<SRA> sras,
                                         std::string threads, std::string ram_gb,
                                         bool mult_sra, bool dispOutput, bool retainInterFiles,
                                         std::string logFile, const INI_MAP & cfgIni);
