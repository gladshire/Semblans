#include "rcorr_wrap.h"

void run_rcorr(std::vector<SRA> sras, std::string threads) {
  std::string outDir(sras[0].get_sra_path_corr().first.parent_path().native().c_str());
  std::string inFile1;
  std::string inFile2;
  
  std::string maxCor;
  std::string maxCorK;
  std::string wkProp;
  
  std::string rcorrCmd("perl " + PATH_RCORR + " -t " + threads + " -od " + outDir);
  for (auto sra : sras) {
    std::cout << "\nRunning error correction for:\n" << std::endl;
    summarize_sing_sra(sra);
    if (sra.is_paired()) {
      inFile1 = sra.get_sra_path_raw().first.c_str();
      inFile2 = sra.get_sra_path_raw().second.c_str();
      system((rcorrCmd + " -1 " + inFile1 + " -2 " + inFile2).c_str()); 
    }
    else {
      inFile1 = sra.get_sra_path_raw().first.c_str();
      system((rcorrCmd + " -s " + inFile1).c_str());
    }
  }
} 
