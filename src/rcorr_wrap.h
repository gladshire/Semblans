#pragma once
#include <boost/filesystem.hpp>
#include <boost/dll.hpp>
#include <fstream>
#include "sra.h"
#include "print_info.h"

namespace fs = boost::filesystem;
namespace dl = boost::dll;

#define PATH_RCORR std::string((dl::program_location().parent_path() / fs::path(std::string("../external/Rcorrector/run_rcorrector.pl"))).c_str())

void run_rcorr(std::pair<std::string, std::string> sraRun, std::string outDir,
               std::string threads, bool dispOutput, std::string logFile);
