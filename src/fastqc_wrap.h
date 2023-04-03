#pragma once
#include <boost/filesystem.hpp>
#include <boost/dll.hpp>
#include "sra.h"
#include "print_info.h"

namespace fs = boost::filesystem;
namespace dl = boost::dll;

#define PATH_FASTQC std::string((dl::program_location().parent_path() / fs::path(std::string("../external/FastQC/fastqc"))).c_str())

namespace fs = boost::filesystem;

void run_fastqc(std::pair<std::string, std::string> sraRunFiles,
                std::string threads, std::string outFile,
                bool dispOutput, std::string logFile);

void run_fastqc_bulk(std::vector<std::pair<std::string, std::string>> sraRunsInput,
                     std::vector<std::string> outFiles, std::string threads,
                     bool dispOutput, std::string logFile);
