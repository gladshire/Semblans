#include "trimm_wrap.h"


void run_trimmomatic(std::pair<std::string, std::string> sraRunIn,
                     std::pair<std::string, std::string> sraRunOutP,
                     std::pair<std::string, std::string> sraRunOutU,
                     std::string threads, bool dispOutput, std::string logFile) {
  std::string inFile1 = sraRunIn.first;
  std::string inFile2 = sraRunIn.second;
  std::string outFileP1 = sraRunOutP.first;
  std::string outFileP2 = sraRunOutP.second;
  std::string outFileU1 = sraRunOutU.first;
  std::string outFileU2 = sraRunOutU.second;

  std::string printOut;
  if (dispOutput) {
    printOut = " 2>&1 | tee -a " + logFile;
  }
  else {
    printOut = " >>" + logFile + " 2>&1";
  }
  int result;
  bool isPaired;
  if (sraRunIn.second != "") {
    isPaired = true;
  }
  else {
    isPaired = false;
  }
  std::string trimmCmd("java -jar " + PATH_TRIMM);
  std::string trimmFlags("-threads " + threads + " " + "ILLUMINACLIP:" + TRUSEQ_ALL +
                         ":2:30:10 SLIDINGWINDOW:4:5 LEADING:5 TRAILING:5 MINLEN:25");

  if (isPaired) {
    trimmCmd += " PE " + inFile1 + " " + inFile2 + " " + outFileP1 + " " + outFileU1 +
                " " + outFileP2 + " " + outFileU2 + " " + trimmFlags + printOut;
    result = system(trimmCmd.c_str());
  }
  else {
    trimmCmd += " SE " + inFile1 + " " + outFileU1 + " " + trimmFlags + printOut;
    result = system(trimmCmd.c_str());
  }
  if (WIFSIGNALED(result)) {
    system("setterm -cursor on");
    logOutput("Exited with signal " + std::to_string(WTERMSIG(result)), logFile);
    exit(1);
  }
}
