#pragma once
#include <iostream>
#include <fstream>
#include <boost/filesystem.hpp>
#include <boost/regex.hpp>
#include <boost/dll.hpp>
#include "sra.h"
#include "transcript.h"
#include "print_info.h"


namespace fs = boost::filesystem;
namespace dl = boost::dll;

#define PATH_DIAMOND std::string((dl::program_location().parent_path() / fs::path(std::string("../external/diamond-2.1.3/diamond"))).c_str())
#define DMAKEDB (PATH_DIAMOND + " makedb")
#define DBLASTX (PATH_DIAMOND + " blastx")
#define DBLASTP (PATH_DIAMOND + " blastp")

void makeDb(std::string pathProtRef, std::string outDir);

void blastxDiam(transcript transcripts, std::string blastDb,
                std::string threads, std::string outDir);

void blastpDiam(std::string pepFilePath, std::string blastDb,
                std::string threads, std::string outFile);
