#pragma once
#include <iostream>
#include <fstream>
#include <boost/filesystem.hpp>
#include <boost/dll.hpp>
#include "print_info.h"

namespace fs = boost::filesystem;
namespace dl = boost::dll;

#define PATH_PANTHER std::string((dl::program_location().parent_path() / fs::path(std::string("../external/pantherScore/pantherScore2.2.pl").c_str())).c_str())

#define PATH_HMMER std::string((dl::program_location().parent_path() / fs::path(std::string("../external/hmmer/bin/").c_str())).c_str())

#define PANTHER_LIB std::string((dl::program_location().parent_path() / fs::path(std::string("../data/panther_db/famlib/rel/PANTHER18.0_altVersion/hmmscoring/PANTHER18.0").c_str())).c_str())

void pantherScore(std::string transPepIn, std::string outFile, std::string threads,
                  bool dispOutput, std::string logFile);
