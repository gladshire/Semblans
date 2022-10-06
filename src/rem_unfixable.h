#include <iostream>
#include <fstream>
#include <thread>
#include <boost/filesystem.hpp>
#include "sra.h"

void rem_unfix_pe(SRA sra);

void rem_unfix_se(SRA sra);

void rem_unfix_bulk(std::vector<SRA> sras, std::string threads);
