#include <iostream>
#include <fstream>
#include <regex>
#include <boost/filesystem.hpp>
#include "sra.h"
#include "thread_pool.h"

void rem_overrep_se(SRA sra);

void rem_overrep_pe(SRA sra);

void rem_overrep_bulk(std::vector<SRA>, std::string threads);
