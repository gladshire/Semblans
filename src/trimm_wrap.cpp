#include "trimm_wrap.h"

// Given an SRA run's sequence data, trim its reads' adapter sequences
void run_trimmomatic(std::pair<std::string, std::string> sraRunIn,
                     std::pair<std::string, std::string> sraRunOutP,
                     std::pair<std::string, std::string> sraRunOutU,
                     std::string threads, std::string maxSeedMismatch,
                     std::string minScorePaired, std::string minScoreSingle,
                     std::string windowSize, std::string windowMinQuality,
                     std::string minQualityLead, std::string minQualityTrail,
                     std::string minReadLength, std::string numBpCutFront,
                     bool dispOutput, std::string logFile) {
  std::string inFile1 = sraRunIn.first;
  std::string inFile2 = sraRunIn.second;
  std::string outFileP1 = sraRunOutP.first;
  std::string outFileP2 = sraRunOutP.second;
  std::string outFileU1 = sraRunOutU.first;
  std::string outFileU2 = sraRunOutU.second;
  int result;
  bool isPaired;
  
  if (sraRunIn.second != "") {
    isPaired = true;
  }
  else {
    isPaired = false;
  }
  std::string trimmCmd("java -jar " + PATH_TRIMM);
  std::string trimmFlags("-threads " + threads + " " + "ILLUMINACLIP:" +
                         TRUSEQ_ALL + ":" + maxSeedMismatch + ":" + minScorePaired + ":" + minScoreSingle +
                         " SLIDINGWINDOW:" + windowSize + ":" + windowMinQuality +
                         " LEADING:" + minQualityLead + " TRAILING:" + minQualityTrail +
                         " MINLEN:" + minReadLength + " HEADCROP:" + numBpCutFront);

  if (isPaired) {
    trimmCmd = trimmCmd + " PE " + " -phred33 " + inFile1 + " " + inFile2 + " " +
               outFileP1 + " " + outFileU1 + " " + outFileP2 + " " + outFileU2 + " " +
               trimmFlags;
  }
  else {
    trimmCmd = trimmCmd + " SE " + " -phred33 " + inFile1 + " " +
               outFileU1 + " " + trimmFlags;
  }
  if (dispOutput) {
    trimmCmd += " 2>&1 | tee -a " + logFile;
    logOutput("  Running command: " + trimmCmd + "\n\n", logFile);
  }
  else {
    trimmCmd += " >>" + logFile + " 2>&1";
  }
  result = system(trimmCmd.c_str());
  checkExitSignal(result, logFile);
}
