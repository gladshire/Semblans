#pragma once
#include <boost/filesystem.hpp>
#include "sra.h"
#include "print_info.h"

#define PATH_RCORR std::string("../lib/Rcorrector-1.0.5/run_rcorrector.pl")

namespace fs = boost::filesystem;

void run_rcorr(std::vector<SRA> sras, std::string threads);
