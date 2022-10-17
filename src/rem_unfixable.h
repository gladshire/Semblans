#pragma once
#include <iostream>
#include <fstream>
#include <boost/filesystem.hpp>
#include "print_info.h"
#include "sra.h"
#include "thread_pool.h"

void rem_unfix_pe(SRA sra);

void rem_unfix_se(SRA sra);

void rem_unfix_bulk(std::vector<SRA> sras, std::string threads);
