#include "rcorr_wrap.h"

// Given an SRA run's sequence files, correct base pair errors with Rcorrector
void run_rcorr(std::pair<std::string, std::string> sraRun, std::string outDir,
               std::string threads, std::string kmerLength, std::string maxCorrK,
               std::string weakProportion, bool dispOutput, bool compressFiles,
               std::string logFile) {
  std::string inFile1 = sraRun.first;
  std::string inFile2 = sraRun.second;

  std::string maxCor;
  std::string maxCorK;
  std::string wkProp;

  std::string printOut;
  
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
               " -k " + kmerLength + " -maxcorK " + maxCorrK + " -wk " + weakProportion +
               " -tmpd " + outDir;
  }
  else {
    rcorrCmd = "( perl " + PATH_RCORR + " -t " + threads + " -s " + inFile1 +
               " -k " + kmerLength + " -maxcorK " + maxCorrK + " -wk " + weakProportion +
               " -tmpd " + outDir;
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
  if (dispOutput) {
    rcorrCmd += " 2>&1 | tee -a " + logFile;
    logOutput("  Running command: " + rcorrCmd + "\n\n", logFile);
  }
  else {
    rcorrCmd += " >>" + logFile + " 2>&1";
  }
  result = system(rcorrCmd.c_str());
  if (WIFSIGNALED(result)) {
    system("setterm -cursor on");
    logOutput("Exited with signal " + std::to_string(WTERMSIG(result)), logFile);
    system("rm ./tmp*");
    exit(1);
  }
}
