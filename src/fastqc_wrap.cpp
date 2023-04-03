#include "fastqc_wrap.h"


void run_fastqc(std::pair<std::string, std::string> sraRunFiles,
                std::string threads, std::string outFile,
                bool dispOutput, std::string logFile) {
  std::string inFile1 = sraRunFiles.first;
  std::string inFile2 = sraRunFiles.second;
  bool isPaired;
  if (sraRunFiles.second != "") {
    isPaired = true;
  }
  else {
    isPaired = false;
  }
  std::string fastqcFlags = " --extract -t " + threads + " -o ";
  int result;
  system(("mkdir " + outFile).c_str());
  std::string fastqcCmd;
  if (isPaired) {
    fastqcCmd = PATH_FASTQC + fastqcFlags + outFile + " " + inFile1 + " " + inFile2;
    if (dispOutput) {
      fastqcCmd += (" 2>&1 | tee -a " + logFile);
    }
    else {
      fastqcCmd += (" >>" + logFile + " 2>&1");
    }
    result = system(fastqcCmd.c_str());
  }
  else {
    fastqcCmd = PATH_FASTQC + fastqcFlags + outFile + " " + inFile1;
    if (dispOutput) {
      fastqcCmd += (" 2>&1 | tee -a " + logFile);
    }
    else {
      fastqcCmd += (" >>" + logFile + " 2>&1");
    }
    result = system((PATH_FASTQC + fastqcFlags + outFile + " " + inFile1).c_str());
  }
  if (WIFSIGNALED(result)) {
    std::cout << "Exited with signal " << WTERMSIG(result) << std::endl;
    exit(1);
  }
}

void run_fastqc_bulk(std::vector<std::pair<std::string, std::string>> sraRunsInput,
                     std::vector<std::string> outFiles, std::string threads,
                     bool dispOutput, std::string logFile) {
  //std::cout << "\nRunning quality analysis for:\n" << std::endl;
  //summarize_all_sras(sras, logFile, 2);
  /*for (auto sra : sras) {
    // Check if checkpoint exists
    if (outDir == std::string(sra.get_fastqc_dir_1().first.parent_path().parent_path().c_str())) {
      if (sra.checkpointExists("fastqc1")) {
        logOutput("FastQC analysis found for:", logFile);
        summarize_sing_sra(sra, logFile, 2);
        continue;
      }
    }
    else {
      if (sra.checkpointExists("fastqc2")) {
        logOutput("FastQC analysis found for:", logFile);
        summarize_sing_sra(sra, logFile, 2);
        continue;
      }
    }
    run_fastqc(sra, threads, outDir, dispOutput, logFile);
    // Make checkpoint file
    if (outDir == std::string(sra.get_fastqc_dir_1().first.parent_path().parent_path().c_str())) {
      sra.makeCheckpoint("fastqc1");
    }
    else {
      sra.makeCheckpoint("fastqc2");
    }
  }*/
  for (int i = 0; i < sraRunsInput.size(); i++) {
    run_fastqc(sraRunsInput[i], threads, outFiles[i], dispOutput, logFile);
  }
}
