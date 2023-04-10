#pragma once
#include <iostream>
#include <fstream>
#include <set>
#include <vector>
#include <boost/filesystem.hpp>
#include <boost/regex.hpp>
#include "sra.h"
#include "transcript.h"
#include "seq.h"
#include "seq_hash.h"
#include "print_info.h"

namespace fs = boost::filesystem;
namespace dl = boost::dll;

std::set<std::string> makeClusterSet(std::ifstream & clustFile);

void filterCorset(std::string transIn, std::string transClust,
                  std::string transLargestClust, std::string transRedund,
                  uintmax_t ram_b, std::string outDir, std::string logFile);
