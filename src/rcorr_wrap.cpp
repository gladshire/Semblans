#include "rcorr_wrap.h"


void run_rcorr(std::vector<std::pair<std::string, std::string>> sraRunFiles,
               std::string outDir, std::string threads, bool dispOutput,
               std::string logFile) {
  /*logOutput("\nRunning error correction for:\n", logFile);
  summarize_all_sras(sras, logFile, 2);
  */
  std::string inFile1;
  std::string inFile2;

  std::string maxCor;
  std::string maxCorK;
  std::string wkProp;

  std::string printOut;
  if (dispOutput) {
    printOut = " 2>&1 | tee -a " + logFile;
  }
  else {
    printOut = " >>" + logFile + " 2>&1";
  }

  bool isPaired;
  int result;
  for (auto sraRun : sraRunFiles) {
    if (sraRun.second != "") {
      isPaired = true;
    }
    else {
      isPaired = false;
    }
    inFile1 = sraRun.first;
    // Check for checkpoint file
    /*if (sra.checkpointExists("corr")) {
      logOutput("Error-corrected version found for: ", logFile);
      summarize_sing_sra(sra, logFile, 2);
      continue;
    }*/
    std::string rcorrCmd = "perl " + PATH_RCORR + " -t " + threads + " -od " + outDir;
    if (isPaired) {
      inFile2 = sraRun.second;
      rcorrCmd += " -1 " + inFile1 + " -2 " + inFile2;
      result = system((rcorrCmd + printOut).c_str()); 
    }
    else {
      rcorrCmd += " -s " + inFile1;
      result = system((rcorrCmd + printOut).c_str());
    }
    if (WIFSIGNALED(result)) {
      logOutput("Exited with signal " + std::to_string(WTERMSIG(result)), logFile);
      exit(1);
    }
    // Create checkpoint
    //sra.makeCheckpoint("corr");
  }
} 
