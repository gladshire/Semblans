#include "fastqc_wrap.h"

// Perform a quality analysis of a given SRA run with FastQC
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
  system(("mkdir " + outFile + " > /dev/null 2>&1").c_str());
  std::string fastqcCmd;
  std::string printOut;
  if (isPaired) {
    fastqcCmd = PATH_FASTQC + fastqcFlags + outFile + " " + inFile1 + " " + inFile2;
  }
  else {
    fastqcCmd = PATH_FASTQC + fastqcFlags + outFile + " " + inFile1;
  }
  if (dispOutput) {
    fastqcCmd += " 2>&1 | tee -a " + logFile;
    logOutput("  Running command: " + fastqcCmd + "\n\n", logFile);
  }
  else {
    fastqcCmd += " >>" + logFile + " 2>&1";
  }
  result = system(fastqcCmd.c_str());
  if (WIFSIGNALED(result)) {
    system("setterm -cursor on");
    std::cout << "Exited with signal " << WTERMSIG(result) << std::endl;
    exit(1);
  }
}

void run_fastqc_bulk(std::vector<std::pair<std::string, std::string>> sraRunsInput,
                     std::vector<std::string> outFiles, std::string threads,
                     bool dispOutput, std::string logFile) {
  for (int i = 0; i < sraRunsInput.size(); i++) {
    run_fastqc(sraRunsInput[i], threads, outFiles[i], dispOutput, logFile);
  }
}
