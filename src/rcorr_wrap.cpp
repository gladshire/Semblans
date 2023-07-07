#include "rcorr_wrap.h"

void run_rcorr(std::pair<std::string, std::string> sraRun, std::string outDir,
               std::string threads, bool dispOutput, bool compressFiles,
               std::string logFile) {
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
  std::string rcorrCmd;
  if (isPaired) {
    rcorrCmd = "( perl " + PATH_RCORR + " -t " + threads + " -1 " + inFile1 + " -2 " + inFile2 +
               " -verbose";
  }
  else {
    rcorrCmd = "( perl " + PATH_RCORR + " -t " + threads + " -s " + inFile1 +
               " -verbose";
  }
  if (compressFiles) {
    std::string sraPathL = outDir + "/" + std::string(fs::path(sraRun.first.c_str()).stem().stem().c_str()) + ".cor.fq.gz";
    std::string sraPathR = outDir + "/" + std::string(fs::path(sraRun.second.c_str()).stem().stem().c_str()) + ".cor.fq.gz";
    rcorrCmd += std::string(" -stdout | awk \'{ if ((NR-1) % 8 < 4) {print | \"" +
                            PATH_PIGZ + " --fast -p " + threads + " > " + sraPathL + "\"} " +
                            "else {print | \"" + PATH_PIGZ + " --fast -p " + threads + " > " + sraPathR + "\"} }\' )");
  }
  else {
    rcorrCmd += std::string(" -od " + outDir + ")");
  }
  rcorrCmd += printOut;
  result = system(rcorrCmd.c_str());
  if (WIFSIGNALED(result)) {
    system("setterm -cursor on");
    logOutput("Exited with signal " + std::to_string(WTERMSIG(result)), logFile);
    system("rm ./tmp*");
    exit(1);
  }
}
