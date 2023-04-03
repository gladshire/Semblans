// TODO: Allow use of local data
//   - Process input as comma-separated pairs ending with .fastq / .fq


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
    run_fastqc_bulk(sras, threads, fastqc_dir_1, dispOutput, logFilePath);
    
    // Run rcorrector on reads
    run_rcorr(sras, threads, dispOutput, logFilePath); 
    
    // Remove reads with unfixable errors
    rem_unfix_bulk(sras, ram_gb, logFilePath);
    
    // If not keeping intermediate files, remove error-corrected files
    if (!retainInterFiles) {
      for (auto sra : sras) {
        fs::remove(fs::path(sra.get_sra_path_corr().first.c_str()));
        if (sra.is_paired()) {
          fs::remove(fs::path(sra.get_sra_path_corr().second.c_str()));
        }
      }
    }

    // Run trimmomatic to remove adapter seqs
    run_trimmomatic(sras, threads, dispOutput, logFilePath);

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
    run_fastqc_bulk(sras, threads, fastqc_dir_2, dispOutput, logFilePath);

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
