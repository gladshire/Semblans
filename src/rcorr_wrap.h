#pragma once
#include <boost/filesystem.hpp>
#include <boost/dll.hpp>
#include "sra.h"
#include "print_info.h"

namespace fs = boost::filesystem;
namespace dl = boost::dll;

#define PATH_RCORR std::string((dl::program_location().parent_path() / fs::path(std::string("../external/Rcorrector-1.0.5/run_rcorrector.pl"))).c_str())

void run_rcorr(const std::vector<SRA> & sras, std::string threads);
