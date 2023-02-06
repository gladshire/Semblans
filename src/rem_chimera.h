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

void detect_chimera(transcript trans, std::string blastxFile, std::string outDir);

void align_buffer_end(std::ifstream & inFile, char * inFileData, std::streamsize & s);

/*
unsigned long hashFunction(char * key);

void insertHash(std::vector<sequence> * fastaHashTable, uintmax_t lenHashTable,
                std::string header, std::string sequence);

void deleteHash(std::vector<sequence> * fastaHashTable, uintmax_t lenHashTable,
                std::string header);

bool inHashTable(std::vector<sequence> * fastaHashTable, uintmax_t lenHashTable,
                 std::string header);
*/
void fillFastaHash(seqHash fastaHashTable, transcript trans, uintmax_t ram_b);

std::set<std::string> makeChimeraSet(std::ifstream & chimFile);

void removeChimera(transcript trans, std::string infoFilePath, std::string cutFilePath,
                   uintmax_t ram_b, std::string outDir);
