#pragma once
#include <iostream>
#include <cstring>
#include <algorithm>
#include <fstream>
#include <boost/filesystem.hpp>
#include "print_info.h"
#include "sra_toolkit.h"
#include "sra.h"

void rem_unfix_pe(std::pair<std::string, std::string> sraRunIn,
                  std::pair<std::string, std::string> sraRunOut,
                  uintmax_t ram_b);

void rem_unfix_se(std::string sraRunIn, std::string sraRunOut,
                  uintmax_t ram_b);

void rem_unfix_bulk(std::vector<std::pair<std::string, std::string>> sraRunsIn,
                    std::vector<std::pair<std::string, std::string>> sraRunsOut,
                    std::string ram_gb, std::string logFile);
