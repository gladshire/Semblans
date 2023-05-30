#include "assemble.h"

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

void makeTransInfoFile(const std::vector<SRA> & sras, std::string transInfoFileStr) {
  std::ofstream transInfoFile;
  transInfoFile.open(transInfoFileStr);
  for (auto sra : sras) {
    transInfoFile << std::string(sra.get_sra_path_orep_filt().first.c_str());
    if (sra.is_paired()) {
      transInfoFile << " " << std::string(sra.get_sra_path_orep_filt().second.c_str());
    }
    transInfoFile << std::endl;
  }
  transInfoFile.close();
}

std::vector<transcript> run_trinity_bulk(std::vector<SRA> sras,
                                         std::string threads, std::string ram_gb,
                                         bool mult_sra, bool dispOutput, bool retainInterFiles,
                                         std::string logFile, const INI_MAP & cfgIni) {
  logOutput("Starting de-novo assembly for:", logFile);
  summarize_all_sras(sras, logFile, 2);
  INI_MAP_ENTRY cfgPipeline = cfgIni.at("Pipeline");
  std::vector<transcript> sra_transcripts;
  std::string outDir;
  std::vector<std::pair<std::string, std::string>> sraRunsInTrin;
  std::pair<std::string, std::string> currTrinIn;
  std::string currTrinOut;
  for (auto sra : sras) {
    transcript currSraTrans(sra);
    // Check for checkpoint file
    if (sra.checkpointExists("trinity")) {
      logOutput("Assembly checkpoint found for: ", logFile);
      summarize_sing_sra(sra, logFile, 2);
      continue;
    }

    currTrinIn.first = sra.get_sra_path_orep_filt().first.c_str();
    currTrinIn.second = sra.get_sra_path_orep_filt().second.c_str();
    
    if (!ini_get_bool(cfgPipeline.at("remove_overrepresented").c_str(), 0)) {
      if (!ini_get_bool(cfgPipeline.at("filter_foreign_reads").c_str(), 0)) {
        if (!ini_get_bool(cfgPipeline.at("trim_adapter_seqs").c_str(), 0)) {
          if (!ini_get_bool(cfgPipeline.at("error_correction").c_str(), 0)) {
            currTrinIn.first = sra.get_sra_path_raw().first.c_str();
            currTrinIn.second = sra.get_sra_path_raw().second.c_str();
          }
          else {
            currTrinIn.first = sra.get_sra_path_corr_fix().first.c_str();
            currTrinIn.second = sra.get_sra_path_corr_fix().second.c_str();
          }
        }
        else {
          currTrinIn.first = sra.get_sra_path_trim_p().first.c_str();
          currTrinIn.second = sra.get_sra_path_trim_p().second.c_str();
        }
      }
      else {
        currTrinIn.first = sra.get_sra_path_for_filt().first.c_str();
        currTrinIn.second = sra.get_sra_path_for_filt().second.c_str();
      }
    }

    currTrinOut = currSraTrans.get_trans_path_trinity().c_str();
    if (mult_sra) {
      std::vector<std::pair<std::string, std::string>> sraTrinInComb;
      std::vector<SRA> sras_comb = get_sra_to_combine(sras, sra.get_org_name());
      logOutput("Combined assembly chosen using:", logFile);
      summarize_all_sras(sras_comb, logFile, 2);
      for (auto sra : sras_comb) {
        currTrinIn.first = sra.get_sra_path_orep_filt().first.c_str();
        currTrinIn.second = sra.get_sra_path_orep_filt().second.c_str();
    
        if (!ini_get_bool(cfgPipeline.at("remove_overrepresented").c_str(), 0)) {
          if (!ini_get_bool(cfgPipeline.at("filter_foreign_reads").c_str(), 0)) {
            if (!ini_get_bool(cfgPipeline.at("trim_adapter_seqs").c_str(), 0)) {
              if (!ini_get_bool(cfgPipeline.at("error_correction").c_str(), 0)) {
                currTrinIn.first = sra.get_sra_path_raw().first.c_str();
                currTrinIn.second = sra.get_sra_path_raw().second.c_str();
              }
              else {
                currTrinIn.first = sra.get_sra_path_corr_fix().first.c_str();
                currTrinIn.second = sra.get_sra_path_corr_fix().second.c_str();
              }
            }
            else {
              currTrinIn.first = sra.get_sra_path_trim_p().first.c_str();
              currTrinIn.second = sra.get_sra_path_trim_p().second.c_str();
            }
          }
          else {
            currTrinIn.first = sra.get_sra_path_for_filt().first.c_str();
            currTrinIn.second = sra.get_sra_path_for_filt().second.c_str();
          }
        }
        sraTrinInComb.push_back(currTrinIn);
      }
      run_trinity_comb(sraTrinInComb, currTrinOut, threads, ram_gb, dispOutput, logFile);
      // Make file for transcript containing its associated SRAs
      std::string transInfoFileStr(currSraTrans.get_trans_path_trinity().replace_extension(".transInfo").c_str());
      makeTransInfoFile(sras_comb, transInfoFileStr);
    }   
    else {
      logOutput("Single assembly chosen using:", logFile);
      summarize_sing_sra(sra, logFile, 2);
      run_trinity(currTrinIn, currTrinOut, threads, ram_gb, dispOutput, logFile);
      // Make file for transcript containing its associated SRA
      std::string sraInfoFileStr(currSraTrans.get_trans_path_trinity().replace_extension(".transInfo").c_str());
      std::vector<SRA> singSra{sra};
      makeTransInfoFile(singSra, sraInfoFileStr);
    }
    sra_transcripts.push_back(currSraTrans);
  }
  return sra_transcripts;
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

    // Create file space
    make_proj_space(cfgIni, "assemble");

    // Obtain SRAs
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
          logOutput("ERROR: Local run not found: \"" + sraRunsLocal.first + "\"", logFilePath);
        }
        if (sraRunsLocal.second != "" &&
            !fs::exists(cfgIni["General"]["local_data_directory"] + sraRunsLocal.second)) {
          logOutput("ERROR: Local run not found: \"" + sraRunsLocal.second + "\"", logFilePath);
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
    // Get boolean for intermediate file fate
    std::string retainInterFilesStr = argv[5];
    bool retainInterFiles = stringToBool(argv[5]);
    // Get boolean for verbose printing
    bool dispOutput = stringToBool(argv[6]);


    logOutput("Paando Assemble started with following parameters:", logFilePath);
    logOutput("  Config file:     " + std::string(argv[1]), logFilePath);;
    logOutput("  Threads (Cores): " + threads, logFilePath);
    logOutput("  Memory (GB):     " + ram_gb, logFilePath);
    logOutput("  SRA Runs:\n", logFilePath);
    summarize_all_sras(sras, logFilePath, 6);
    // Perform assembly with Trinity
    std::vector<transcript> transcriptsSra = run_trinity_bulk(sras, threads, ram_gb, mult_sra,
                                                              dispOutput, retainInterFiles,
                                                              logFilePath, cfgIni);
  }
  else {
    print_help();
  }
  return 0;
}
