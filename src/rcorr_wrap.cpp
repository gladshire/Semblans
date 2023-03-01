#include "rcorr_wrap.h"


void run_rcorr(const std::vector<SRA> & sras, std::string threads,
               bool dispOutput, std::string logFile) {
  std::cout << "\nRunning error correction for:\n" << std::endl;
  summarize_all_sras(sras);

  std::string inFile1;
  std::string inFile2;
  
  std::string maxCor;
  std::string maxCorK;
  std::string wkProp;

  std::string printOut;
  if (dispOutput) {
    printOut = " |& tee -a " + logFile;
  }
  else {
    printOut = " &>> " + logFile;
  }

  int result;
  for (auto sra : sras) {
    std::string outDir(sra.get_sra_path_corr().first.parent_path().c_str());
    inFile1 = sra.get_sra_path_raw().first.c_str();
    if (fs::exists(sra.get_sra_path_corr().first)) {
      std::cout << "Error-corrected version found for: " << std::endl;
      summarize_sing_sra(sra);
      continue;
    }
    std::string rcorrCmd = "perl " + PATH_RCORR + " -t " + threads + " -od " + outDir;
    if (sra.is_paired()) {
      inFile2 = sra.get_sra_path_raw().second.c_str();
      rcorrCmd += " -1 " + inFile1 + " -2 " + inFile2;
      result = system((rcorrCmd + printOut).c_str()); 
    }
    else {
      rcorrCmd += " -s " + inFile1;
      result = system((rcorrCmd + printOut).c_str());
    }
    if (WIFSIGNALED(result)) {
      std::cout << "Exited with signal " << WTERMSIG(result) << std::endl;
      exit(1);
    }
  }
} 
