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

#define PIDENT_CUTOFF 30
#define LENGTH_CUTOFF 100

namespace fs = boost::filesystem;


std::vector<std::string> getLineVec(std::string line);

double getQcov(std::vector<std::string> currLineVec);

std::vector<std::string> expRange(std::vector<std::string> & hsp1,
                                  std::vector<std::string> & hsp2);

bool isSeparate(std::vector<std::string> hsp1, std::vector<std::string> hsp2);

bool checkBlock(std::vector<std::vector<std::string>> & block, bool multiGene,
                std::ofstream & outFile1, std::ofstream & outFile2);

void detect_chimera(std::string blastxFile, std::string outDir);

std::set<std::string> makeChimeraSet(std::ifstream & chimFile);

void removeChimera(std::string transIn, std::string transOut,
                   std::string infoFilePath, std::string cutFilePath,
                   std::string ram_gb, std::string logFile);
