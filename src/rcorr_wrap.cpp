#include "rcorr_wrap.h"

void run_rcorr(std::pair<std::string, std::string> sraRun, std::string outDir,
               std::string threads, bool dispOutput, std::string logFile) {
  std::string inFile1 = sraRun.first;
  std::string inFile2 = sraRun.second;

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
  if (sraRun.second != "") {
    isPaired = true;
  }
  else {
    isPaired = false;
  }
  std::string rcorrCmd = "perl " + PATH_RCORR + " -t " + threads + " -od " + outDir;
  if (isPaired) {
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
}
