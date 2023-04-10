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

#define PATH_BLAST std::string((dl::program_location().parent_path() / fs::path(std::string("../external/ncbi-blast-2.13.0+/bin/"))).c_str())
#define MAKEBLASTDB (PATH_BLAST + "makeblastdb")
#define BLASTX (PATH_BLAST + "blastx")
#define BLASTP (PATH_BLAST + "blastp")


void makeBlastDb(std::string pathProtRef, std::string outDir,
                 bool dispOutput, std::string logFile);

void blastx(std::string transIn, std::string blastDb,
            std::string threads, std::string outDir,
            bool dispOutput, std::string logFile);

void blastp(std::string pepFilePath, std::string blastDb,
            std::string threads, std::string outFile,
            bool dispOutput, std::string logFile);
