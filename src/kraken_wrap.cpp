#include "kraken_wrap.h"

// Obtain the paths to kraken2 databases from user config INI file
std::vector<std::string> get_kraken2_dbs(const INI_MAP &iniFile) {
  std::string krakDbDir = iniFile.at("Kraken2 settings").at("db_directory");
  std::vector<std::string> kraken2Dbs;
  fs::path dbPath;
  if (iniFile.at("Kraken2 filter order").empty()) {
    return kraken2Dbs;
  }
  for (auto db : iniFile.at("Kraken2 filter order")) {
    dbPath = fs::path(krakDbDir.c_str()) / fs::path(db.first.c_str());
    if (!fs::exists(dbPath)) {
      std::cerr << "ERROR: Kraken2 database '" << dbPath.c_str() << "' not found." << std::endl;
      exit(1);
    }
    kraken2Dbs.push_back(dbPath.c_str());
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
                 std::string threads, std::string db, std::string confThreshold,
                 std::string minBaseQuality, std::string minHitGroups,
                 bool keepForeign, bool dispOutput, bool compressFiles,
                 std::string logFile) {

  std::string inFile1 = sraRunIn.first;
  std::string inFile2 = sraRunIn.second;
  std::string outDir = std::string(fs::path(sraRunOut).parent_path().c_str());
  std::string krakCmd(PATH_KRAK2 + " --db " + db);
  std::string krakFlags(" --threads " + threads + " --confidence " + confThreshold +
                        " --minimum-base-quality " + minBaseQuality +
                        " --minimum-hit-groups " + minHitGroups);
  std::string krakOutput;
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
  sraRunOutClass = sraRunOutDir + "/" + sraRunOutClass + "." +
                   std::string(fs::path(db.c_str()).stem().c_str()) + ".fastq";
  if (compressFiles) {
    krakOutput = "";
  }
  else {
    krakOutput = " --output - --unclassified-out " + sraRunOut;
    if (keepForeign) {
      krakOutput += " --classified-out " + sraRunOutClass;
    }
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
    krakCmd = krakCmd + " " + krakFlags + " --paired " + krakOutput + " " + inFile1 + " " + inFile2 +
              " --report " + repFile;
  }
  else {
    krakCmd = krakCmd + " " + krakFlags + " " + krakOutput + " " + inFile1 + " " +
              " --report " + repFile + " ";
  }
  if (dispOutput) {
    krakCmd += " 2>&1 | tee -a " + logFile;
    logOutput("  Running command: " + krakCmd + "\n", logFile);
  }
  else {
    krakCmd += " >>" + logFile + " 2>&1";
  }
  result = system(krakCmd.c_str());
  checkExitSignal(result, logFile);
}
