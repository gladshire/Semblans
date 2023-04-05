// TODO:
//   - Customize input file locations according to pipeline parameters in config.ini
//   - Alter wrappers to take filenames as input/output rather than SRA objects (modularity)

#include "preprocess.h"

void retrieve_sra_data(const std::vector<SRA> & sras, std::string threads,
                       bool dispOutput, std::string logFile) {
  std::ofstream logStream;
  logStream.open(logFile, std::ios_base::app);
  teedev logger(logStream, std::cout);
  teeStream loggerStream(logger);
  loggerStream << "Retrieving SRA runs for:\n" << std::endl;
  logStream.close();
  summarize_all_sras(sras, logFile, 2);

  prefetch_sra(sras, dispOutput, logFile);
  fasterq_sra(sras, threads, dispOutput, logFile);
}

bool stringToBool(std::string boolStr) {
  bool boolConv;
  for (int i = 0; i < boolStr.length(); i++) {
    boolStr[i] = std::tolower(boolStr[i]);
  }
  boolConv = (boolStr == "true") ? true : false;
  return boolConv;
}

void print_help() {
  winsize w;
  ioctl(STDOUT_FILENO, TIOCGWINSZ, &w);
  std::cout << std::left << std::setw(w.ws_col) << "COMMAND STRUCTURE" << std::endl;
  std::cout << std::left << std::setw(w.ws_col) << "preprocess PATH/TO/CONFIG.INI num_threads RAM_GB" << std::endl;
}


int main(int argc, char * argv[]) {
  if (argc < 4) {
    print_help();
    return 0;
  }
  else {
    INI_MAP cfgIni = make_ini_map(argv[1]);
    std::string threads = argv[2];
    std::string ram_gb = argv[3];
    std::string retainFiles = argv[4];
    std::vector<SRA> sras;
    std::vector<std::string> localDataFiles;
    bool localDataBool = false;
    bool retainInterFiles = stringToBool(argv[4]);
    bool dispOutput = stringToBool(argv[5]);
    std::string logFilePath = cfgIni["General"]["log_file"];
    // Create vector of SRA objects from SRA accessions, using NCBI web API
    sras = get_sras(cfgIni);

    if (!sras.empty()) {
      retrieve_sra_data(sras, threads, dispOutput, logFilePath);
    }
    // Get single/paired filenames of local data
    for (auto fqFileName : cfgIni.at("Local files")) {
      localDataFiles.push_back(fqFileName.first);
    }
    std::pair<std::string, std::string> sraRunsLocal;
    size_t pos;
    for (auto sraRun : localDataFiles) {
      sraRunsLocal.first = "";
      sraRunsLocal.second = "";
      pos = sraRun.find(" ");
      sraRunsLocal.first = sraRun.substr(0, pos);
      if (pos != std::string::npos) {
        sraRun.erase(0, pos + 1);
        pos = sraRun.find(" ");
        sraRunsLocal.second = sraRun.substr(0, pos);
      }
      // Check if files in pair exist. If not, do not add to sras vector
      if (fs::exists(cfgIni["General"]["local_data_directory"] + sraRunsLocal.first) &&
          fs::exists(cfgIni["General"]["local_data_directory"] + sraRunsLocal.second)) {
        sras.push_back(SRA(sraRunsLocal.first, sraRunsLocal.second, cfgIni));
      }
      else {
        if (sraRunsLocal.first != "" &&
            !fs::exists(cfgIni["General"]["local_data_directory"] + sraRunsLocal.first)) {
          std::cout << "ERROR: Local run not found: \"" << sraRunsLocal.first << "\""
                    << std::endl;
        }
        if (sraRunsLocal.second != "" &&
            !fs::exists(cfgIni["General"]["local_data_directory"] + sraRunsLocal.second)) {
          std::cout << "ERROR: Local run not found: \"" << sraRunsLocal.second << "\""
                    << std::endl;
        }
      }
    }
    
    // Check if no SRAs specified
    if (sras.empty()) {
      std::cout << "ERROR: No SRA runs specified. Please check config file" << std::endl;
    }

    // Summarize preprocess task
    logOutput("Paando Preprocess started with following parameters:", logFilePath);
    logOutput("  Config file:     " + std::string(argv[1]), logFilePath);
    logOutput("  Threads (Cores): " + threads, logFilePath);
    logOutput("  Memory (GB):     " + ram_gb, logFilePath);
    logOutput("  SRA runs:\n", logFilePath); 
    summarize_all_sras(sras, logFilePath, 6);
    logOutput("  Retain intermediate files: " + std::string(argv[4]), logFilePath);

    std::string fastqc_dir_1(sras[0].get_fastqc_dir_1().first.parent_path().parent_path().c_str());
    std::string fastqc_dir_2(sras[0].get_fastqc_dir_2().first.parent_path().parent_path().c_str());
    std::vector<std::string> kraken2Dbs = get_kraken2_dbs(cfgIni);
    std::string kraken2_conf = get_kraken2_conf(cfgIni);
    std::pair<std::vector<std::string>, std::vector<std::string>> overrepSeqs;
    make_proj_space(cfgIni);

    // Run initial fastqc on reads
    std::pair<std::string, std::string> currFastqcIn1;
    std::string currFastqcOut1;
    for (auto sra : sras) {
      currFastqcIn1.first = sra.get_sra_path_raw().first.c_str();
      currFastqcIn1.second = sra.get_sra_path_raw().second.c_str();
      currFastqcOut1 = sra.get_fastqc_dir_1().first.parent_path().c_str();
      // Check for checkpoint file
      if (sra.checkpointExists("fastqc1")) {
        logOutput("FastQC analysis found for:", logFilePath);
        summarize_sing_sra(sra, logFilePath, 2);
        continue;
      }
      run_fastqc(currFastqcIn1, threads, currFastqcOut1, dispOutput, logFilePath);
      // Make checkpoint file
      sra.makeCheckpoint("fastqc1");
    }
    
    // Error-correction stage
    std::pair<std::string, std::string> currRcorrIn;
    std::string rcorrOutDir = sras[0].get_sra_path_corr().first.parent_path().c_str();
    for (auto sra : sras) {
      currRcorrIn.first = sra.get_sra_path_raw().first.c_str();
      currRcorrIn.second = sra.get_sra_path_raw().second.c_str();
      
      // Check for checkpoint file
      if (sra.checkpointExists("corr")) {
        logOutput("Error-corrected version found for:", logFilePath);
        summarize_sing_sra(sra, logFilePath, 2);
        continue;
      }
      run_rcorr(currRcorrIn, rcorrOutDir, threads, dispOutput, logFilePath);
      // Make checkpoint file
      sra.makeCheckpoint("corr");
    }

    // Remove reads with unfixable errors
    std::pair<std::string, std::string> currCorrFixIn;
    std::pair<std::string, std::string> currCorrFixOut;
    uintmax_t ram_b = (uintmax_t)stoi(ram_gb) * 1000000000;
    for (auto sra : sras) {
      currCorrFixIn.first = sra.get_sra_path_corr().first.c_str();
      currCorrFixIn.second = sra.get_sra_path_corr().second.c_str();

      currCorrFixOut.first = sra.get_sra_path_corr_fix().first.c_str();
      currCorrFixOut.second = sra.get_sra_path_corr_fix().second.c_str();

      // Check for checkpoint file
      if (sra.checkpointExists("corr.fix")) {
        logOutput("Unfixable error-fixed version found for:", logFilePath);
        summarize_sing_sra(sra, logFilePath, 2);
        continue;
      }
      if (sra.is_paired()) {
        rem_unfix_pe(currCorrFixIn, currCorrFixOut, ram_b); 
      }
      else {
        rem_unfix_se(currCorrFixIn.first, currCorrFixOut.first, ram_b);
      }
      // Make checkpoint file
      sra.makeCheckpoint("corr");
    }
  
    // If not keeping intermediate files, remove error-corrected files
    if (!retainInterFiles) {
      for (auto sra : sras) {
        fs::remove(fs::path(sra.get_sra_path_corr().first.c_str()));
        if (sra.is_paired()) {
          fs::remove(fs::path(sra.get_sra_path_corr().second.c_str()));
        }
      }
    }

    // Adapter sequence trimming stage
    std::pair<std::string, std::string> currTrimIn;
    std::pair<std::string, std::string> currTrimOutP;
    std::pair<std::string, std::string> currTrimOutU;
    for (auto sra : sras) {
      currTrimIn.first = sra.get_sra_path_corr_fix().first.c_str();
      currTrimIn.second = sra.get_sra_path_corr_fix().second.c_str();

      currTrimOutU.first = sra.get_sra_path_trim_u().first.c_str();
      currTrimOutU.second = "";

      currTrimOutP.first = "";
      currTrimOutP.second = "";
      if (sra.is_paired()) {
        currTrimOutU.second = sra.get_sra_path_trim_u().second.c_str();

        currTrimOutP.first = sra.get_sra_path_trim_p().first.c_str();
        currTrimOutP.second = sra.get_sra_path_trim_p().second.c_str();
      }
      // Check for checkpoint file
      if (sra.checkpointExists("trim")) {
        logOutput("Adapter-trimmed version found for:", logFilePath);
        summarize_sing_sra(sra, logFilePath, 2);
        continue;
      }
      run_trimmomatic(currTrimIn, currTrimOutP, currTrimOutU, threads,
                      dispOutput, logFilePath);
      // Make checkpoint file
      sra.makeCheckpoint("trim");
    }
    // If not keeping intermediate files, remove fixed error-corrected files
    if (!retainInterFiles) {
      for (auto sra : sras) {
        fs::remove(fs::path(sra.get_sra_path_corr_fix().first.c_str()));
        if (sra.is_paired()) {
          fs::remove(fs::path(sra.get_sra_path_corr_fix().second.c_str()));
        }
      }
    }
   
    // Run kraken2 to remove foreign reads
    // for all db in dbs
    //   for all sra in sras
    //     if db == dbs.first
    //       run_kraken2(sra.from_last, db)
    //     if db != dbs.first
    //       run_kraken2(sra.from_filt, db)
    //       
    //       
    std::pair<std::string, std::string> currKrakIn;
    std::string krakOutDir;
    std::string currKrakOut;
    std::string currRepStem;
    std::string currRepFile;

    run_kraken2_dbs(sras, threads, kraken2Dbs, kraken2_conf, dispOutput, logFilePath);

    // If not keeping intermediate files, remove trimmomatic output files
    if (!retainInterFiles) {
      for (auto sra : sras) {
        fs::remove(fs::path(sra.get_sra_path_trim_p().first.c_str()));
        if (sra.is_paired()) {
          fs::remove(fs::path(sra.get_sra_path_trim_p().second.c_str()));
        }
      }
    }

    // Run fastqc on all runs
    std::vector<std::pair<std::string, std::string>> fastqcIn2;
    std::pair<std::string, std::string> currFastqcIn2;
    std::vector<std::string> fastqcOut2;
    for (auto sra : sras) {
      currFastqcIn2.first = sra.get_sra_path_for_filt().first.c_str();
      currFastqcIn2.second = sra.get_sra_path_for_filt().second.c_str();
      fastqcIn2.push_back(currFastqcIn1);

      fastqcOut2.push_back(sra.get_fastqc_dir_2().first.parent_path().c_str());
    }
    run_fastqc_bulk(fastqcIn2, fastqcOut2, threads, dispOutput, logFilePath);

    // Remove reads with over-represented sequences
    rem_overrep_bulk(sras, ram_gb, logFilePath);

    // If not keeping intermediate files, remove kraken2 output files
    if (!retainInterFiles) {
      for (auto sra : sras) {
        fs::remove(fs::path(sra.get_sra_path_for_filt().first.c_str()));
        fs::remove(fs::path(fastqc_dir_2.c_str()));
        if (sra.is_paired()) {
          fs::remove(fs::path(sra.get_sra_path_for_filt().second.c_str()));
        }
      }
    }
  }

  return 0;
}
