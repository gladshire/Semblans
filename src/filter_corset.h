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

/*
void align_buffer_end(std::ifstream & inFile, char * inFileData,
                      std::streamsize & s);

void fillFastaHash(seqHash fastaHashTable, transcript trans, uintmax_t ram_b);
*/
std::set<std::string> makeClusterSet(std::ifstream & clustFile);

void filterCorset(transcript trans, std::string clusterPath, uintmax_t ram_b,
                  std::string out_dir, std::string logFile);
