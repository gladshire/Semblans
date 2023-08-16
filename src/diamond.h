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

#define PATH_DIAMOND std::string((dl::program_location().parent_path() / fs::path(std::string("../external/diamond/diamond"))).c_str())
#define DMAKEDB (PATH_DIAMOND + " makedb")
#define DBLASTX (PATH_DIAMOND + " blastx")
#define DBLASTP (PATH_DIAMOND + " blastp")

void makeBlastDbDiam(std::string pathProtRef, std::string outDir,
                     bool dispOutput, std::string logFile);

void blastxDiam(std::string transIn, std::string blastDb,
                std::string maxEvalue, std::string maxTargetSeqs,
                std::string threads, std::string outDir,
                bool dispOutput, std::string logFile);

void blastpDiam(std::string pepFilePath, std::string blastDb,
                std::string maxEvalue, std::string maxTargetSeqs,
                std::string threads, std::string outFile,
                bool dispOutput, std::string logFile);

