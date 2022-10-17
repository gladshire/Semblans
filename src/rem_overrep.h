#pragma once
#include <iostream>
#include <fstream>
#include <sstream>
#include <boost/filesystem.hpp>
#include <boost/regex.hpp>
#include "sra.h"

std::pair<std::vector<std::string>, std::vector<std::string>> get_overrep_seqs_pe(SRA sra);

std::vector<std::string> get_overrep_seqs_se(SRA sra); 
