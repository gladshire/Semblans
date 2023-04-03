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


void run_kraken2(const std::vector<SRA> & sras, std::string threads, std::string db, 
                 std::string conf_threshold, bool selfPass, bool dispOutput,
                 std::string logFile) {
  std::string inFile1;
  std::string inFile2;
  std::string outFile;
  std::string repFile;
  fs::path dbPath(db.c_str());

  std::string krakCmd(PATH_KRAK2 + " --db " + db);
  std::string krakFlags("--threads " + threads + " --unclassified-out");
  std::string printOut;
  if (dispOutput) {
    printOut = " 2>&1 | tee -a " + logFile;
  }
  else {
    printOut = " >>" + logFile + " 2>&1";
  }
  int result;
  for (auto sra : sras) {
    std::string outDir(sra.get_sra_path_for_filt().first.parent_path().c_str());
    repFile = std::string(outDir + "/" + sra.get_file_prefix().first + "." +
                          dbPath.filename().c_str() + ".report");
    // Check for checkpoint
    if (sra.checkpointExists(std::string(dbPath.stem().c_str()) + ".filt")) {
      logOutput("With database: " + db, logFile);
      logOutput("Filtered version found for: ", logFile);
      summarize_sing_sra(sra, logFile, 2);
      continue;
    }
    pre_summary(sra, db, logFile);
    if (selfPass) {
      std::string tmpIn1 = std::string(sra.get_sra_path_for_filt().first.parent_path().c_str()) +
                           "/INPUT1.fq";
      std::rename(sra.get_sra_path_for_filt().first.c_str(), tmpIn1.c_str());
      inFile1 = tmpIn1;
    }
    else {
      inFile1 = sra.get_sra_path_trim_p().first.c_str();
    }
    outFile = sra.get_sra_path_for_filt().first.c_str();
    if (sra.is_paired()) {
      if (selfPass) {
        std::string tmpIn2 = std::string(sra.get_sra_path_for_filt().second.parent_path().c_str()) +
                             "/INPUT2.fq";
        std::rename(sra.get_sra_path_for_filt().second.c_str(), tmpIn2.c_str());
        inFile2 = tmpIn2;
      }
      else {
        inFile2 = sra.get_sra_path_trim_p().second.c_str();
      }
      outFile = std::string(outFile).replace(outFile.length() - 10, 2, "#");
      result = system((krakCmd + " " + krakFlags + " " + outFile + " --paired " + "--output - " +
                       inFile1 + " " + inFile2 + " --confidence " + conf_threshold + " --report " +
                       outDir + "/" + sra.get_file_prefix().first + "." + dbPath.filename().c_str() +
                       ".report" + printOut).c_str());
    }
    else {
      result = system((krakCmd + " " + krakFlags + " " + outFile + "--output - " + inFile1 + " " +
                       " --confidence " + conf_threshold + " --report " + outDir + "/" +
                       sra.get_file_prefix().first + "." + dbPath.filename().c_str() + ".report" +
                       printOut).c_str());
    }
    if (WIFSIGNALED(result)) {
      logOutput("Exited with signal " + std::to_string(WTERMSIG(result)), logFile);
      exit(1);
    }
    // Create checkpoint
    sra.makeCheckpoint(std::string(dbPath.stem().c_str()) + ".filt");
  }
}


void run_kraken2_dbs(const std::vector<SRA> & sras, std::string threads, std::vector<std::string> dbs,
                     std::string conf_threshold, bool dispOutput, std::string logFile) {
  logOutput("\nFiltering foreign reads for:\n", logFile);
  summarize_all_sras(sras, logFile, 2);
  int i = 0;
  for (auto db : dbs) {
    logOutput("\nNow filtering: " + std::string(fs::path(db.c_str()).filename().c_str()),
              logFile);
    if (i == 0) {
      run_kraken2(sras, threads, db, conf_threshold, false, dispOutput, logFile);
    }
    else {
      run_kraken2(sras, threads, db, conf_threshold, true, dispOutput, logFile);
    }
    i++;
  }
}
