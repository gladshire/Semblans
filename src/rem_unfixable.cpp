#include <iostream>
#include <thread>
#include <boost/filesystem.hpp>
#include "sra.h"

void remove_unfixable_reads(std::vector<SRA> sras, std::string threads) {
  std::string outDir(sras[0].get_sra_path_corr().first.parent_path().native().c_str());
  std::string inFile1;
  std::string inFile2;
  std::string outFile1;
  std::string outFile2;

  for (auto sra : sras) {
    if (sra.is_paired()) {
      std::string inFile1(sra.get_sra_path_corr().first.native().c_str());
      std::string inFile2(sra.get_sra_path_corr().second.native().c_str());
      std::string outFile1(sra.get_sra_path_corr().first.stem().native().c_str() +
                           ".fix.fastq");
      std::string outFile2(sra.get_sra_path_corr().second.stem().native().c_str() +
                           ".fix.fastq");
      
    }
    else {

    }
  }
