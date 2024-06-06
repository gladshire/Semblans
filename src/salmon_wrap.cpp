#include "salmon_wrap.h"

// Index reads from sequence data using Salmon
void salmon_index(std::string transIn, std::string transIndex, std::string threads,
                  bool dispOutput, std::string logFile) {
  int result;

  std::string salmCmd = PATH_SALMON + " index" + " -t " + transIn +
                        " -i " + transIndex + " -p " + threads;
  if (dispOutput) {
    salmCmd += " 2>&1 | tee -a " + logFile;
    logOutput("  Running command: " + salmCmd + "\n\n", logFile);
  }
  else {
    salmCmd += " >>" + logFile + " 2>&1";
  }
  result = system(salmCmd.c_str());
  checkExitSignal(result, logFile);
}

bool runPaired(std::vector<std::pair<std::string, std::string>> sraRunsIn) {
  int numSingle = 0;
  int numPaired = 0;
  for (auto sraRun : sraRunsIn) {
    if (sraRun.second != "" &&
        sraRun.second != "null") {
      numPaired++;
    }
    else {
      numSingle++;
    }
  }
  if (numPaired >= numSingle) {
    return true;
  }
  else {
    return false;
  }
}

// Quantify reads from sequence data. Assumes all reads are either single or paired-end
void salmon_quant(std::string transIn, std::string transIndex, std::string transQuant,
                  std::vector<std::pair<std::string, std::string>> sraRunsIn,
                  std::string threads, bool dispOutput, std::string logFile) {
  std::string transFilePath(transIn);
  std::string indexFilePath(transIndex);
  std::string quantFilePath(transQuant);
  
  std::string sras1 = "";
  std::string sras2 = "";

  // Determine if more data exists for paired or single runs contituting transcript
  bool morePaired = false;
  morePaired = runPaired(sraRunsIn);
  
  // Construct largest possible list of either paired or single runs
  for (int i = 0; i < sraRunsIn.size(); i++) {
    if (morePaired) {
      if (sraRunsIn[i].second != "" &&
          sraRunsIn[i].second != "null") {
        sras1 += (sraRunsIn[i].first + " ");
        sras2 += (sraRunsIn[i].second + " ");
      }
    }
    else {
      if (sraRunsIn[i].second == "" ||
          sraRunsIn[i].second == "null") {
        sras1 += (sraRunsIn[i].first + " ");
      }
    }
  }

  int result;
  if (morePaired) {
    std::string salmCmd = PATH_SALMON + " quant" + " -i " + indexFilePath + " --dumpEq" +
                          " --writeUnmappedNames " + " --libType" + " A" + " -p " + threads +
                          " -1 " + sras1 + " -2 " + sras2 + " -o " + quantFilePath;
    if (dispOutput) {
      salmCmd += " 2>&1 | tee -a " + logFile;
      logOutput("  Running command: " + salmCmd + "\n\n", logFile);
    }
    else {
      salmCmd += " >>" + logFile + " 2>&1";
    }
    result = system(salmCmd.c_str());
  }
  else {
    std::string salmCmd = PATH_SALMON + " quant" + " -i " + indexFilePath + " --dumpEq" +
                          " --writeUnmappedNames " + " --libType" + " A" + " -p " + threads +
                          " -r " + sras1 + " -o " + quantFilePath;
    if (dispOutput) {
      salmCmd += " 2>&1 | tee -a " + logFile;
      logOutput("  Running command: " + salmCmd + "\n\n", logFile);
    }
    else {
      salmCmd += " >>" + logFile + " 2>&1";
    }
    result = system(salmCmd.c_str());
  }
  replaceChar(logFile, '\r', '\n');
  checkExitSignal(result, logFile);
}
