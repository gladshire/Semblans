#include "salmon_wrap.h"

// TODO: Fix this function. Does not read files in properly
void prepareDecoys(std::vector<std::string> desiredFiles, 
                   std::vector<std::string> decoyFiles,
                   std::string preIndexOutput,
                   std::string decoyFileOutput) {
  std::ifstream currSeqFile;
  std::ofstream decoyFile(decoyFileOutput);
  std::ofstream preIndexFile(preIndexOutput);
  std::string currLine;
  // Iterate through files with sequence targets for mapping
  for (auto currTarget : desiredFiles) {
    currSeqFile.open(currTarget);
    while (std::getline(currSeqFile, currLine)) {
      if (currLine[0] == '>' || currLine[0] == '@') {
        /*for (int i = 0; i < currLine.size(); i++) {
          if (currLine[i] == ' ') {
            currLine[i] = '_';
          }
        }*/
        currLine = currLine.substr(0, currLine.find(" "));
      }
      preIndexFile << currLine << std::endl;
    }
    currSeqFile.close();
  }
  // Iterate through files with decoy sequences
  for (auto currDecoy : decoyFiles) {
    currSeqFile.open(currDecoy);
    while (std::getline(currSeqFile, currLine)) {
      if (currLine[0] == '>' || currLine[0] == '@') {
        /*for (int i = 0; i < currLine.size(); i++) {
          if (currLine[i] == ' ') {
            currLine[i] = '_';
          }
        }*/
        currLine = currLine.substr(0, currLine.find(" "));
        decoyFile << currLine << std::endl;
      }
      preIndexFile << currLine << std::endl;
    }
    currSeqFile.close();
  }
  decoyFile.close();
  preIndexFile.close();
}

// Index reads from sequence data using Salmon
void salmon_index(std::string transIn, std::string transIndex, std::string decoys,
                  std::string threads, bool dispOutput, std::string logFile) {
  std::cout << transIn << std::endl;
  std::cout << decoys << std::endl;

  std::string printOut;
  if (dispOutput) {
    printOut = " 2>&1 | tee -a " + logFile;
  }
  else {
    printOut = " >>" + logFile + " 2>&1";
  }
  int result;

  std::string salm_cmd = PATH_SALMON + " index" + " -t " + transIn +
                         " -i " + transIndex + " --decoys " + decoys +
                         " -p " + threads + printOut;
  result = system(salm_cmd.c_str());
  if (WIFSIGNALED(result)) {
    logOutput("Exited with signal " + std::to_string(WTERMSIG(result)), logFile);
    exit(1);
  }
}

bool runPaired(std::vector<std::pair<std::string, std::string>> sraRunsIn) {
  int numSingle = 0;
  int numPaired = 0;
  for (auto sraRun : sraRunsIn) {
    if (sraRun.second != "") {
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
      if (sraRunsIn[i].second != "") {
        sras1 += (sraRunsIn[i].first + " ");
        sras2 += (sraRunsIn[i].second + " ");
      }
    }
    else {
      if (sraRunsIn[i].second == "") {
        sras1 += (sraRunsIn[i].first + " ");
      }
    }
  }

  std::string printOut;
  if (dispOutput) {
    printOut = " 2>&1 | tee -a " + logFile;
  }
  else {
    printOut = " >>" + logFile + " 2>&1";
  }
  int result;
  
  if (morePaired) {
    std::string salm_cmd = PATH_SALMON + " quant" + " -i " + indexFilePath + " --dumpEq" +
                           " --writeUnmappedNames " + " --libType" + " A" + " -p " + threads +
                           " -1 " + sras1 + " -2 " + sras2 +
                           " --validateMappings" + " -o " + quantFilePath + printOut;
    result = system(salm_cmd.c_str());
  }
  else {
    std::string salm_cmd = PATH_SALMON + " quant" + " -i " + indexFilePath + " --dumpEq" +
                           " --writeUnmappedNames " + " --libType" + " A" + " -p " + threads +
                           " -r " + sras1 + " --validateMappings" + " -o " + quantFilePath +
                           printOut;
    result = system(salm_cmd.c_str());
  }
  if (WIFSIGNALED(result)) {
    system("setterm -cursor on");
    logOutput("Exited with signal " + std::to_string(WTERMSIG(result)), logFile);
    exit(1);
  }
}
