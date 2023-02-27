#include "rcorr_wrap.h"


void run_rcorr(const std::vector<SRA> & sras, std::string threads) {
  std::cout << "\nRunning error correction for:\n" << std::endl;
  summarize_all_sras(sras);

  std::string inFile1;
  std::string inFile2;
  
  std::string maxCor;
  std::string maxCorK;
  std::string wkProp;
 
  int result;
  for (auto sra : sras) {
    std::string outDir(sra.get_sra_path_corr().first.parent_path().c_str());
    std::string rcorrCmd("perl " + PATH_RCORR + " -t " + threads + " -od " + outDir);
    inFile1 = sra.get_sra_path_raw().first.c_str();
    if (fs::exists(sra.get_sra_path_corr().first)) {
      std::cout << "Error-corrected version found for: " << sra.get_accession() << std::endl;
      continue;
    }
    if (sra.is_paired()) {
      inFile2 = sra.get_sra_path_raw().second.c_str();
      result = system((rcorrCmd + " -1 " + inFile1 + " -2 " + inFile2).c_str()); 
    }
    else {
      result = system((rcorrCmd + " -s " + inFile1).c_str());
    }
    if (WIFSIGNALED(result)) {
      std::cout << "Exited with signal " << WTERMSIG(result) << std::endl;
      exit(1);
    }
  }
} 
