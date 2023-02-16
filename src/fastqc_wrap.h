#pragma once
#include <boost/filesystem.hpp>
#include <boost/dll.hpp>
#include "sra.h"
#include "print_info.h"

namespace fs = boost::filesystem;
namespace dl = boost::dll;

#define PATH_FASTQC std::string((dl::program_location().parent_path() / fs::path(std::string("../external/FastQC/fastqc"))).c_str())

namespace fs = boost::filesystem;

void run_fastqc(SRA sra, std::string threads, std::string outDir);

void run_fastqc_bulk(std::vector<SRA> sras, std::string threads, std::string outDir);
