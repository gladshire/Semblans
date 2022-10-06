#pragma once
#include <boost/filesystem.hpp>
#include "sra.h"
#include "print_info.h"

#define PATH_FASTQC std::string("../lib/FastQC/fastqc")

namespace fs = boost::filesystem;

void run_fastqc(std::vector<SRA> sras, std::string threads);
