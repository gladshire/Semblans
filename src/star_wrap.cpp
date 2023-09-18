#include "star_wrap.h"

void star_index(std::vector<std::string> fastaFiles, std::string outDir, std::string threads,
                bool dispOutput, std::string logFile) {
  std::string fastasString = "";
  std::string printOut;
  if (dispOutput) {
    printOut = " 2>&1 | tee -a " + logFile;
  }
  else {
    printOut = " >>" + logFile + " 2>&1";
  }
  for (auto fastaFile : fastaFiles) {
    fastasString += (fastaFile + " ");
  }
  int result;
  std::string star_cmd = "( " + PATH_STAR + " --runMode genomeGenerate " +
                         " --genomeFastaFiles " + fastasString +
                         " --genomeDir " + outDir +
                         " --runThreadN " + threads + " )" + printOut;
  result = system(star_cmd.c_str());
  if (WIFSIGNALED(result)) {
    logOutput("Exited with signal " + std::to_string(WTERMSIG(result)), logFile);
    exit(1);
  }
}


void star_map(std::string indexPath, std::string outMap,
              std::pair<std::string, std::string> sraReadsIn,
              std::string maxMapsPerRead, std::string threads,
              bool dispOutput, std::string logFile) {
  std::string sraReadsString;
  std::string printOut;
  if (dispOutput) {
    printOut = " 2>&1 | tee -a " + logFile;
  }
  else {
    printOut = " >>" + logFile + " 2>&1";
  }
  sraReadsString = sraReadsIn.first;
  if (!sraReadsIn.second.empty()) {
    sraReadsString += (" " + sraReadsIn.second);
  }
  int result;
  std::string star_cmd = "( " + PATH_STAR + " --genomeDir " + indexPath +
                         " --readFilesIn " + sraReadsString +
                         " --outFilterMultimapNmax " + maxMapsPerRead +
                         " --outFileNamePrefix " + outMap +
                         " --runThreadN " + threads + " )" + printOut;
  result = system(star_cmd.c_str());
  if (WIFSIGNALED(result)) {
    logOutput("Exited with signal " + std::to_string(WTERMSIG(result)), logFile);
    exit(1);
  }
}
