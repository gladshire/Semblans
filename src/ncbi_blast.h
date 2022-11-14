#include <iostream>
#include <fstream>
#include <boost/filesystem.hpp>
#include <boost/regex.hpp>
#include "sra.h"
#include "transcript.h"
#include "print_info.h"


#define PATH_BLAST std::string("../lib/ncbi-blast-2.13.0+-src/c++/ReleaseMT/bin/")
#define MAKEBLASTDB (PATH_BLAST + "makeblastdb")
#define BLASTX (PATH_BLAST + "blastx")

namespace fs = boost::filesystem;


void makeBlastDb(std::string pathProfRef, std::string outDir);

void blastx(transcript transcripts, std::string blastDb,
            std::string threads, std::string outDir);


