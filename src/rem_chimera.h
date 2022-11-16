#pragma once
#include <iostream>
#include <fstream>
#include <boost/filesystem.hpp>
#include <boost/regex.hpp>
#include "sra.h"
#include "transcript.h"
#include "print_info.h"

#define PIDENT_CUTOFF 30
#define LENGTH_CUTOFF 100

namespace fs = boost::filesystem;


void detect_chimera(std::string blastxFile, std::string outDir);

void rem_chimera(transcript transcripts, std::string infoFile);

void rem_chim_bulk(std::vector<transcript> transcriptsV, std::string infoFile);
