#include "sra.h"
#include <boost/filesystem.hpp>

#define PATH_FASTQC "../lib/FastQC/fastqc"

void run_fastqc(std::vector<SRA> sras, std::string threads) {
  std::string outDir(sras[0].get_fastqc_dir().first.parent_path().parent_path().native().c_str());
  std::string outFile;
  std::string inFile1;
  std::string inFile2;
  std::string fastqcFlags = " --extract -t " + threads + " -q -o ";
  for (auto sra : sras) {
    outFile = std::string(outDir + sra.make_file_str());
    system(("mkdir " + outFile).c_str());
    if (!sra.is_paired()) {
      inFile1 = std::string(sra.get_sra_path_raw().first.c_str());
      inFile2 = std::string(sra.get_sra_path_raw().second.c_str());
      system((PATH_FASTQC + fastqcFlags + outFile + " " + inFile1 + " " + inFile2).c_str());
    }
    else {
      inFile1 = std::string(sra.get_sra_path_raw().first.c_str());
      system((PATH_FASTQC + fastqcFlags + outFile + " " + inFile1).c_str());
    }
  } 
}
