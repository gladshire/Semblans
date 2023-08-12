#pragma once
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cstring>
#include <algorithm>
#include <fstream>
#include <boost/filesystem.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/copy.hpp>
#include "print_info.h"
#include "sra_toolkit.h"
#include "sra.h"

namespace io = boost::iostreams;


bool rem_unfix_pe(std::pair<std::string, std::string> sraRunIn,
                  std::pair<std::string, std::string> sraRunOut,
                  uintmax_t ram_b, bool dispOutput, bool compressFiles,
                  std::string logFile);

bool rem_unfix_se(std::string sraRunIn, std::string sraRunOut,
                  uintmax_t ram_b, bool dispOutput, bool compressFiles,
                  std::string logFile);
