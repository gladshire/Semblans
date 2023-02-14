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

#define PATH_BLAST std::string((dl::program_location().parent_path() / fs::path(std::string("../lib/ncbi-blast-2.13.0+-src/c++/ReleaseMT/bin/"))).c_str())
#define MAKEBLASTDB (PATH_BLAST + "makeblastdb")
#define BLASTX (PATH_BLAST + "blastx")
#define BLASTP (PATH_BLAST + "blastp")

// Add wrapper for blastp

void makeBlastDb(std::string pathProfRef, std::string outDir);

void blastx(transcript transcripts, std::string blastDb,
            std::string threads, std::string outDir);

void blastp(std::string pepFilePath, std::string blastDb,
            std::string threads, std::string outFile);
