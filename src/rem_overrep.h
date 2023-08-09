#pragma once
#include <iostream>
#include <fstream>
#include <boost/filesystem.hpp>
#include <boost/regex.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/copy.hpp>
#include "sra.h"
#include "sra_toolkit.h"
#include "print_info.h"
#include "seq.h"
#include "seq_hash.h"

namespace io = boost::iostreams;


//seqHash * get_overrep_seqs_pe(SRA sra, int & len1, int & len2);
std::pair<std::vector<std::string>, std::vector<std::string>> get_overrep_seqs_pe(SRA sra);

std::vector<std::string> get_overrep_seqs_se(SRA sra); 

bool rem_overrep_pe(std::pair<std::string, std::string> sraRunIn,
                    std::pair<std::string, std::string> sraRunOut,
                    uintmax_t ram_b, bool compressFiles,
                    std::pair<std::vector<std::string>, std::vector<std::string>> overrepSeqs);

bool rem_overrep_se(std::string sraRunIn, std::string sraRunOut,
                    uintmax_t ram_b, bool compressFiles,
                    std::vector<std::string> overrepSeqs);

//void rem_overrep_bulk(std::vector<std::pair<std::string, std::string>> sraRunsIn,
//                      std::string ram_gb, std::string logFile);
