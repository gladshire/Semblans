#pragma once
#include <iostream>
#include <fstream>
#include <boost/filesystem.hpp>
#include <boost/regex.hpp>
#include "sra.h"
#include "sra_toolkit.h"
#include "print_info.h"

std::pair<std::vector<std::string>, std::vector<std::string>> get_overrep_seqs_pe(SRA sra);

std::vector<std::string> get_overrep_seqs_se(SRA sra); 

void rem_overrep_pe(SRA sra, long long int ram_b,
                    std::pair<std::vector<std::string>, std::vector<std::string>> overrepSeqs);

void rem_overrep_se(SRA sra, long long int ram_b, std::vector<std::string> overrepSeqs);

void rem_overrep_bulk(const std::vector<SRA> & sras, std::string ram_gb);
