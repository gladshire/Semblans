#pragma once
#include <iostream>
#include <cstring>
#include <algorithm>
#include <fstream>
#include <boost/filesystem.hpp>
#include "print_info.h"
#include "sra_toolkit.h"
#include "sra.h"

void rem_unfix_pe(SRA sra, long long int ram_b);

void rem_unfix_se(SRA sra, long long int ram_b);

void rem_unfix_bulk(const std::vector<SRA> & sras, std::string ram_gb, std::string logFile);
