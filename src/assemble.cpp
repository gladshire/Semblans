#include "assemble.h"

// TODO:
//   Allow assembly of local short reads
//   Otherwise retrieve from overrepresented folder

std::vector<transcript> get_transcript(std::vector<SRA> sras) {
  std::vector<transcript> transcripts;
  for (auto &sra : sras) {
    transcript currTrans(sra);
    transcripts.push_back(currTrans);
  }
  return transcripts;
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
  std::cout << "\n" << "NAME_OF_PROGRAM" << " - "
            << "A tool for bulk assemblies of de novo transcriptome data" << std::endl;
  std::cout << "\n" << "COMMAND STRUCTURE" << std::endl;
  std::cout << "\n" << "assemble PATH/TO/CONFIG.INI num_threads RAM_GB (--mult)" << std::endl;
}

int main(int argc, char * argv[]) {
  if (argc > 1) {
    std::vector<SRA> sras;
    std::vector<std::string> localDataFiles;
    // Retrieve SRA objects for trinity runs
    INI_MAP cfgIni = make_ini_map(argv[1]);
    std::string logFilePath = cfgIni["General"]["log_file"];

    sras = get_sras(cfgIni);
    
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
    // Get number of threads
    std::string threads = argv[2];
    // Get RAM in GB
    std::string ram_gb = argv[3];
    // Get boolean for multiple sra processing
    std::string mult_sra_str;
    mult_sra_str = argv[4];
    bool mult_sra = stringToBool(argv[4]);
    
    bool dispOutput = stringToBool(argv[5]);
    std::cout << "  ASSEMBLE started with following parameters:" << std::endl;
    std::cout << "    Config file:     " << argv[1] << std::endl;
    std::cout << "    Threads (Cores): " << threads << std::endl;
    std::cout << "    Memory (GB):     " << ram_gb << std::endl;
    std::cout << "    SRAs:" << std::endl;
    summarize_all_sras(sras);
    // Perform assembly with Trinity
    std::vector<transcript> transcriptsSra = run_trinity_bulk(sras, threads, ram_gb, mult_sra,
                                                              dispOutput, logFilePath);
  }
  else {
    print_help();
  }
  return 0;
}
