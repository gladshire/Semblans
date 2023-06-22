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

void makeGroupCheckpoint(std::string cpDir, std::string prefix) {
  std::string cpFileName = cpDir + "/" + prefix + ".trinity.ok";
  std::ofstream cpFile;
  cpFile.open(cpFileName);
  cpFile.close();
}

bool groupCheckpointExists(std::string cpDir, std::string prefix) {
  std::string cpFileName = cpDir + "/" + prefix + ".trinity.ok";
  fs::path cpFilePath(cpFileName.c_str());
  if (fs::exists(cpFilePath)) {
    return true;
  }
  else {
    return false;
  }
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


void run_trinity_bulk(std::map<std::string, std::vector<SRA>> sraGroups,
                      std::string threads, std::string ram_gb,
                      bool dispOutput, bool retainInterFiles,
                      std::string logFile, const INI_MAP & cfgIni) {
  logOutput("Starting de-novo assembly (this may take awhile)", logFile);
  INI_MAP_ENTRY cfgPipeline = cfgIni.at("Pipeline");
  INI_MAP_ENTRY assembGroups = cfgIni.at("Assembly groups");
  std::string outDir;
  std::vector<std::pair<std::string, std::string>> sraRunsInTrin;
  std::pair<std::string, std::string> currTrinIn;
  std::string currTrinOut;
  fs::path cpDir;
  fs::path trinDir;
  std::string transInfoFileStr;

  for (auto sraGroup : sraGroups) {
    SRA dummySra = sraGroup.second.at(0);
    cpDir = dummySra.get_fastqc_dir_2().first.parent_path().parent_path().parent_path() / "checkpoints";

    // Check whether checkpoint exists
    if (groupCheckpointExists(std::string(cpDir.c_str()), sraGroup.first)) {
      logOutput("Assembly found for group: " + sraGroup.first, logFile);
      logOutput("  Containing:", logFile);
      summarize_all_sras(sraGroup.second, logFile, 4);
      continue;
    }

    transcript dummyTrans(dummySra);    
    fs::path pathTrinDir = dummyTrans.get_trans_path_trinity().parent_path();

    outDir = std::string(pathTrinDir.c_str());
    currTrinOut = outDir + "/" + sraGroup.first + ".Trinity.fasta";

    for (auto sra : sraGroup.second) {
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
      sraRunsInTrin.push_back(currTrinIn);
    }
    
    if (sraRunsInTrin.size() > 1) {
      // Multi-assembly
      run_trinity_comb(sraRunsInTrin, currTrinOut, threads, ram_gb, dispOutput, logFile);
    }
    else if (sraRunsInTrin.size() == 1) {
      // Single-assembly
      run_trinity(sraRunsInTrin.at(0), currTrinOut, threads, ram_gb, dispOutput, logFile);
    }
    else {
      std::cout << "This should never happen. You fucked up" << std::endl;
    }
    // Make file for transcript containing associated SRAs
    transInfoFileStr = std::string(fs::path(currTrinOut.c_str()).replace_extension(".ti").c_str());
    makeTransInfoFile(sraGroup.second, transInfoFileStr);
    // Create checkpoint for assembly group
    makeGroupCheckpoint(std::string(cpDir.c_str()), sraGroup.first);
  }
}


void print_help() {
  std::cout << "\n" << "NAME_OF_PROGRAM" << " - "
            << "A tool for bulk assemblies of de novo transcriptome data" << std::endl;
  std::cout << "\n" << "COMMAND STRUCTURE" << std::endl;
  std::cout << "\n" << "assemble PATH/TO/CONFIG.INI num_threads RAM_GB (--mult)" << std::endl;
}

int main(int argc, char * argv[]) {
  system("setterm -cursor off");
  if (argc > 1) {
    std::vector<SRA> sras;
    std::vector<std::string> localDataFiles;
    // Get INI config file
    INI_MAP cfgIni = make_ini_map(argv[1]);
    INI_MAP_ENTRY cfgIniGen = cfgIni["General"];

    // Get number of threads
    std::string threads = argv[2];

    // Get RAM in GB
    std::string ram_gb = argv[3];

    // Get boolean for multiple sra processing
    bool mult_sra = stringToBool(argv[4]);

    // Get boolean for intermediate file fate
    bool retainInterFiles = stringToBool(argv[5]);
    
    // Get boolean for verbose printing
    bool dispOutput = stringToBool(argv[6]);
   
    // Get boolean for output file compression
    //bool compressFiles = ini_get_bool(cfgIni["General"]["compress_files"].c_str(), 0);
    bool compressFiles = false;
    // Retrieve SRA objects for trinity runs
    std::string logFilePath((fs::canonical((fs::path(cfgIniGen["output_directory"].c_str()))) /
                             fs::path(cfgIniGen["project_name"].c_str()) /
                             fs::path(cfgIniGen["log_file"].c_str())).c_str());

    // Create file space
    make_proj_space(cfgIni, "assemble");

    // Obtain SRAs
    sras = get_sras(cfgIni, compressFiles);
    
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
        sras.push_back(SRA(sraRunsLocal.first, sraRunsLocal.second, cfgIni, compressFiles));
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

    // Get group specifications for SRAs
    std::map<std::string, std::vector<SRA>> sraGroups;
    INI_MAP_ENTRY assemblyGroups = cfgIni["Assembly groups"];
    std::string currGroupName;
    std::string currIniArrStr;
    std::vector<std::string> iniStrArray;
    std::vector<SRA> currSraGroup;
    // Iterate through user-defined assembly groups in config file
    for (auto assemblyGroup : assemblyGroups) {
      // Get group name and group array string 
      currGroupName = assemblyGroup.first;
      currIniArrStr = assemblyGroup.second;

      // Remove all whitespace from string
      currIniArrStr.erase(remove_if(currIniArrStr.begin(), currIniArrStr.end(), isspace),
                          currIniArrStr.end());
      
      // Tokenize array string into vector of strings
      iniStrArray = getStrArray(currIniArrStr, ",");

      // Match vector's strings with SRAs, filling SRA group vector
      bool sraChosen;
      std::vector<bool> iniEntryFound(iniStrArray.size(), false);
      for (auto sra : sras) {
        sraChosen = false;
        for (int i = 0; i < iniStrArray.size(); i++) {
          // If SRA specified in group, push to group vector
          if (iniStrArray[i] == sra.get_accession() ||
              iniStrArray[i] == sra.get_file_prefix().first ||
              iniStrArray[i] == sra.get_file_prefix().second ||
              iniStrArray[i] == sra.get_file_prefix().first + ".fastq" ||
              iniStrArray[i] == sra.get_file_prefix().second + ".fastq") {
            iniEntryFound[i] = true;
            if (!sraChosen) {
              currSraGroup.push_back(sra);
              sraChosen = true;
            }
          }
        }
        if (!sraChosen) {
          sraGroups.emplace(sra.get_file_prefix().first,
                            std::vector<SRA>{sra});
        }
      }
      for (int i = 0; i < iniEntryFound.size(); i++) {
        if (iniEntryFound[i] == false) {
          logOutput("ERROR:\n  Entry \"" + iniStrArray[i] + "\" in group \"" +
                    currGroupName + "\" not found.\n  Proceeding without it.\n",
                    logFilePath);
        }
      }
      
      // Emplace current SRA group into map
      sraGroups.emplace(currGroupName, currSraGroup);
      currSraGroup.clear();
    }

    logOutput("Semblans Assemble started with following parameters:", logFilePath);
    logOutput("  Config file:     " + std::string(argv[1]), logFilePath);;
    logOutput("  Threads (Cores): " + threads, logFilePath);
    logOutput("  Memory (GB):     " + ram_gb, logFilePath);
    logOutput("  SRA Runs:\n", logFilePath);
    summarize_all_sras(sras, logFilePath, 6);
    // Perform assembly with Trinity
    run_trinity_bulk(sraGroups, threads, ram_gb, dispOutput,
                     retainInterFiles, logFilePath, cfgIni);
  }
  else {
    print_help();
  }
  system("setterm -cursor on");
  return 0;
}
