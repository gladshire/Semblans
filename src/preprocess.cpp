// TODO: Allow use of local data
//   - Process input as comma-separated pairs ending with .fastq / .fq


#include "preprocess.h"

void retrieve_sra_data(const std::vector<SRA> & sras, std::string threads) {
  std::cout << "Retrieving SRA runs for:\n" << std::endl;
  summarize_all_sras(sras);

  prefetch_sra(sras);
  fasterq_sra(sras, threads);
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

    // Create vector of SRA objects from SRA accessions, using NCBI web API
    sras = get_sras(cfgIni);

    if (!sras.empty()) {
      retrieve_sra_data(sras, threads);
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

    // Summarize preprocess task
    std::cout << "  PREPROCESS started with following parameters:" << std::endl;
    std::cout << "    Config file:     " << argv[1] << std::endl;
    std::cout << "    Threads (Cores): " << threads << std::endl;
    std::cout << "    Memory (GB):     " << ram_gb << std::endl;
    std::cout << "    SRAs (NCBI):\n      ";
    for (auto sra : sras) {
      if (sra.get_accession() != "") {
        std::cout << sra.get_accession() << std::endl;
      }
    }
    std::cout << "    SRAs (local):\n";
    for (auto sra : sras) {
      if (sra.get_accession() == "") {
        if (sra.is_paired()) {
          std::cout << "      Paired-end run:" << std::endl;
          std::cout << "        " << sra.get_file_prefix().first << std::endl;
          std::cout << "        " << sra.get_file_prefix().second << std::endl;
          std::cout << std::endl;
        }
        else {
          std::cout << "      Single-end run:" << std::endl;
          std::cout << "        " << sra.get_file_prefix().first << std::endl;
          std::cout << std::endl;
        }
      }
    }
    std::cout << "    Retain intermediate files: " << argv[4] << std::endl;

    std::string fastqc_dir_1(sras[0].get_fastqc_dir_1().first.parent_path().parent_path().c_str());
    std::string fastqc_dir_2(sras[0].get_fastqc_dir_2().first.parent_path().parent_path().c_str());
    std::vector<std::string> kraken2Dbs = get_kraken2_dbs(cfgIni);
    std::string kraken2_conf = get_kraken2_conf(cfgIni);
    std::pair<std::vector<std::string>, std::vector<std::string>> overrepSeqs;
    make_proj_space(cfgIni);


    run_fastqc_bulk(sras, threads, fastqc_dir_1);
    run_rcorr(sras, threads); 
    rem_unfix_bulk(sras, ram_gb);
    run_trimmomatic(sras, threads);
    if (!retainInterFiles) {
      system(("rm -rf " + std::string(sras[0].get_sra_path_corr_fix().first.parent_path().c_str())).c_str());
    }
    run_kraken2_dbs(sras, threads, kraken2Dbs, kraken2_conf);
    if (!retainInterFiles) {
      system(("rm -rf " + std::string(sras[0].get_sra_path_trim_p().first.parent_path().c_str())).c_str());
    }
    run_fastqc_bulk(sras, threads, fastqc_dir_2);
    rem_overrep_bulk(sras, ram_gb);
    if (!retainInterFiles) {
      system(("rm -rf " + std::string(sras[0].get_sra_path_for_filt().first.parent_path().c_str())).c_str());
      system(("rm -rf " + fastqc_dir_2).c_str());
    }
  }

  return 0;
}
