#pragma once
#include <iostream>
#include <fstream>
#include <sstream>
#include <boost/filesystem.hpp>
#include "print_info.h"
#include "sra.h"
#include "thread_pool.h"

void rem_unfix_pe(SRA sra, long long int ram_b);

void rem_unfix_se(SRA sra, long long int ram_b);

void rem_unfix_bulk(std::vector<SRA> sras, std::string threads, std::string ram_gb);
