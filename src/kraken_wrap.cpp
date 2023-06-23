// TODO: Implement post-filter summary

#include "kraken_wrap.h"

std::vector<std::string> get_kraken2_dbs(const INI_MAP &iniFile) {
  std::string krakDbDir = iniFile.at("Kraken2 settings").at("db_directory");
  std::vector<std::string> kraken2Dbs;
  std::string dbPath;
  for (auto db : iniFile.at("Kraken2 filter order")) {
    dbPath = krakDbDir + db.first;
    kraken2Dbs.push_back(dbPath);
  }
  return kraken2Dbs;
}


std::string get_kraken2_conf(const INI_MAP &iniFile) {
  return iniFile.at("Kraken2 settings").at("confidence_threshold");
}


void pre_summary(SRA sra, std::string db, std::string logFile) {
  logOutput("\nFor SRA accession: " + sra.get_accession(), logFile);
}

void run_kraken2(std::pair<std::string, std::string> sraRunIn,
                 std::string sraRunOut, std::string repFile,
                 std::string threads, std::string db, std::string conf_threshold,
                 bool dispOutput, bool compressFiles, std::string logFile) {

  std::string inFile1 = sraRunIn.first;
  std::string inFile2 = sraRunIn.second;
  std::string outDir = std::string(fs::path(sraRunOut).parent_path().c_str());
  std::string krakCmd(PATH_KRAK2 + " --db " + db);
  std::string krakFlags("--threads " + threads + " --confidence " + conf_threshold);
  std::string krakOutput;
  std::string printOut;
  if (compressFiles) {
    krakOutput = "";
  }
  else {
    krakOutput = " --output - --unclassified-out " + sraRunOut;
  }
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
  if (isPaired) {
    krakCmd = krakCmd + " " + krakFlags + " --paired " + krakOutput + " " + inFile1 + " " + inFile2 + " --report " +
              repFile + printOut;
    result = system(krakCmd.c_str());
  }
  else {
    krakCmd = krakCmd + " " + krakFlags + " " + krakOutput + " " + inFile1 + " " + " --report " + repFile + " " + printOut;
    result = system(krakCmd.c_str());
  }
  if (WIFSIGNALED(result)) {
    system("setterm -cursor on");
    logOutput("Exited with signal " + std::to_string(WTERMSIG(result)), logFile);
    exit(1);
  }
}

