#include "fastqc_wrap.h"


void run_fastqc(SRA sra, std::string threads, std::string outDir,
                bool dispOutput, std::string logFile) {
  std::string inFile1;
  std::string inFile2;
  std::string outFile;
  if (sra.get_accession() != "") {
    outFile = outDir + "/" + sra.make_file_str();
  }
  else {
    outFile = outDir + "/" + sra.get_file_prefix().first;
  }
  std::string fastqcFlags = " --extract -t " + threads + " -o ";
  int result;
  if (fs::exists(fs::path(outFile.c_str()))) {
    std::cout << "FastQC analysis found for:\n" << std::endl;
    if (sra.is_paired()) {
      std::cout << "  Paired-end run:" << std::endl;
      std::cout << "  " << sra.get_file_prefix().first << std::endl;
      std::cout << "  " << sra.get_file_prefix().second << std::endl;
      std::cout << std::endl; 
    }
    else {
      std::cout << "  Single-end run:" << std::endl;
      std::cout << "  " << sra.get_file_prefix().first << std::endl;
      std::cout << std::endl;
    }
    return;
  }
  system(("mkdir " + outFile).c_str());
  std::string fastqcCmd;
  if (sra.is_paired()) {
    // This is not good. Fix at some point.
    if (outDir == std::string(sra.get_fastqc_dir_2().first.parent_path().parent_path().c_str())) {
      inFile1 = sra.get_sra_path_for_filt().first.c_str();
      inFile2 = sra.get_sra_path_for_filt().second.c_str();
    }
    else {
      inFile1 = sra.get_sra_path_raw().first.c_str();
      inFile2 = sra.get_sra_path_raw().second.c_str();
    }
    fastqcCmd = PATH_FASTQC + fastqcFlags + outFile + " " + inFile1 + " " + inFile2;
    if (dispOutput) {
      fastqcCmd += (" |& tee -a " + logFile);
    }
    else {
      fastqcCmd += (" &>> " + logFile);
    }
    result = system(fastqcCmd.c_str());
  }
  else {
    if (outDir == std::string(sra.get_fastqc_dir_2().first.parent_path().parent_path().c_str())) {
      inFile1 = sra.get_sra_path_for_filt().first.c_str();
    }
    else {
      inFile1 = sra.get_sra_path_raw().first.c_str();
    }
    fastqcCmd = PATH_FASTQC + fastqcFlags + outFile + " " + inFile1;
    if (dispOutput) {
      fastqcCmd += (" |& tee -a " + logFile);
    }
    else {
      fastqcCmd += (" &>> " + logFile);
    }
    result = system((PATH_FASTQC + fastqcFlags + outFile + " " + inFile1).c_str());
  }
  if (WIFSIGNALED(result)) {
    std::cout << "Exited with signal " << WTERMSIG(result) << std::endl;
    exit(1);
  }
}

void run_fastqc_bulk(const std::vector<SRA> & sras, std::string threads, std::string outDir,
                     bool dispOutput, std::string logFile) {
  std::cout << "\nRunning quality analysis for:\n" << std::endl;
  summarize_all_sras(sras);
  for (auto sra : sras) {
    run_fastqc(sra, threads, outDir, dispOutput, logFile);
  }
}
