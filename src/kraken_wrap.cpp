// TODO: Implement post-filter summary

#include "kraken_wrap.h"

// Obtain the paths to kraken2 databases from user config INI file
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

// Obtain kraken2 confidence threshold from user config INI file
std::string get_kraken2_conf(const INI_MAP &iniFile) {
  return iniFile.at("Kraken2 settings").at("confidence_threshold");
}


void pre_summary(SRA sra, std::string db, std::string logFile) {
  logOutput("\nFor SRA accession: " + sra.get_accession(), logFile);
}

// Given an SRA run, attempt to classify its reads with a given kraken2 database
// Output both the classified and unclassified reads to new files
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
  std::string sraRunOutDir = fs::path(sraRunOut.c_str()).parent_path().c_str();
  fs::path sraRunOutFile(sraRunIn.first.c_str());
  while (!sraRunOutFile.extension().empty()) {
    sraRunOutFile = sraRunOutFile.stem();
  }
  std::string sraRunOutClass(sraRunOutFile.c_str());
  if (*(&sraRunOutClass.back() - 1) == '_') {
    sraRunOutClass.pop_back();
    sraRunOutClass.back() = '#';
  }
  else {
    sraRunOutClass.push_back('#');
  }
  sraRunOutClass = sraRunOutDir + "/" + sraRunOutClass;
  if (compressFiles) {
    krakOutput = "";
  }
  else {
    
    krakOutput = " --output - --unclassified-out " + sraRunOut +
                 " --classified-out " + sraRunOutClass;
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

