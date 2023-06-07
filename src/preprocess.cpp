#include "preprocess.h"


std::atomic<bool> procRunning(false);


void progressAnim(int numSpace) {
  const std::string anim[] = {".  ", ".. ", "..."};
  int animIndex = 0;

  while (procRunning) {
    std::cout << "\r";
    for (int i = 0; i < numSpace; i++) {
      std::cout << " ";
    }
    std::cout << anim[animIndex] << std::flush;
    animIndex = (animIndex + 1) % 3;
    std::this_thread::sleep_for(std::chrono::milliseconds(1000));
  }
  std::cout << "\r";
  for (int i = 0; i < numSpace; i++) {
    std::cout << " ";
  }
  std::cout << std::endl;
}

void retrieveSraData(const std::vector<SRA> & sras, std::string threads,
                     bool dispOutput, bool compressOutput,
                     bool retainInterFiles,
                     std::string logFilePath) {
  logOutput("Starting retrieval of raw sequence data", logFilePath);
  // Prefetch raw data
  for (auto sra : sras) {
    // Check for checkpoint file
    if (sra.checkpointExists("sra")) {
      logOutput("  Prefetch data found for:   " + sra.get_accession(), logFilePath);
      continue;
    }
    else {
      logOutput("  Downloading raw data for: ", logFilePath);
      summarize_sing_sra(sra, logFilePath, 4);
     
      procRunning = true;
      std::thread prefProgThread(progressAnim,2);
      prefetch_sra(sra, dispOutput, logFilePath);
      procRunning = false;
      prefProgThread.join();

    }
    // Make checkpoint file
    sra.makeCheckpoint("sra");
  }
  // Dump raw data to FASTQ files
  for (auto sra : sras) {
    // Check for checkpoint file
    if (sra.checkpointExists("dump")) {
      logOutput("  Raw FASTQC dump found for: " + sra.get_accession(), logFilePath);
      continue;
    }
    else {
      logOutput("  Dumping to FASTQ file: ", logFilePath);
      summarize_sing_sra(sra, logFilePath, 4);

      procRunning = true;
      std::thread fqdumpThread(progressAnim,2);
      fasterq_sra(sra, threads, dispOutput, compressOutput, logFilePath);
      procRunning = false;
      fqdumpThread.join();

    }
    // Make checkpoint file
    sra.makeCheckpoint("dump");
    if (!retainInterFiles) {
      fs::remove_all(sra.get_sra_path_raw().first.parent_path() / fs::path(sra.get_accession().c_str()));
    }
  }
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

void fastqcBulk1(const std::vector<SRA> & sras, std::string threads, bool dispOutput,
                 std::string logFilePath) {
  logOutput("Starting initial quality analysis", logFilePath);
  std::pair<std::string, std::string> currFastqcIn1;
  std::string currFastqcOut1;
  for (auto sra : sras) {
    currFastqcIn1.first = sra.get_sra_path_raw().first.c_str();
    currFastqcIn1.second = sra.get_sra_path_raw().second.c_str();
    currFastqcOut1 = sra.get_fastqc_dir_1().first.parent_path().c_str();
    // Check for checkpoint file
    if (sra.checkpointExists("fastqc1")) {
      logOutput("  FastQC analysis found for: " + sra.get_accession(), logFilePath);
      continue;
    }
    logOutput("  Now running quality analysis for:", logFilePath);
    summarize_sing_sra(sra, logFilePath, 4);

    procRunning = true;
    std::thread fqcThread(progressAnim,2);
    run_fastqc(currFastqcIn1, threads, currFastqcOut1, dispOutput, logFilePath);
    procRunning = false;
    fqcThread.join();
    
    // Make checkpoint file
    sra.makeCheckpoint("fastqc1");
  }
}

void fastqcBulk2(const std::vector<SRA> & sras, std::string threads, bool dispOutput,
                 std::string logFilePath, const INI_MAP & cfgIni) {
  logOutput("Starting pre-assembly quality analysis", logFilePath);
  INI_MAP_ENTRY cfgPipeline = cfgIni.at("Pipeline");
  std::pair<std::string, std::string> currFastqcIn;
  std::string currFastqcOut;
  for (auto sra : sras) {
    currFastqcIn.first = sra.get_sra_path_for_filt().first.c_str();
    currFastqcIn.second = sra.get_sra_path_for_filt().second.c_str();
    if (!ini_get_bool(cfgPipeline.at("filter_foreign_reads").c_str(), 0)) {
      if (!ini_get_bool(cfgPipeline.at("trim_adapter_seqs").c_str(), 0)) {
        if (!ini_get_bool(cfgPipeline.at("error_correction").c_str(), 0)) {
          currFastqcIn.first = sra.get_sra_path_raw().first.c_str();
          currFastqcIn.second = sra.get_sra_path_raw().second.c_str();
        }
        else {
          currFastqcIn.first = sra.get_sra_path_corr_fix().first.c_str();
          currFastqcIn.second = sra.get_sra_path_corr_fix().second.c_str();
        }
      }
      else {
        currFastqcIn.first = sra.get_sra_path_trim_p().first.c_str();
        currFastqcIn.second = sra.get_sra_path_trim_p().second.c_str();
      }
    }
    currFastqcOut = sra.get_fastqc_dir_2().first.parent_path().c_str();
    // Check for checkpoint file
    if (sra.checkpointExists("fastqc2")) {
      logOutput("  FastQC analysis found for: " + sra.get_accession(), logFilePath);
      continue;
    }
    logOutput("  Now running quality analysis for:", logFilePath);
    summarize_sing_sra(sra, logFilePath, 4);

    procRunning = true;
    std::thread fqcThread(progressAnim,2);
    run_fastqc(currFastqcIn, threads, currFastqcOut, dispOutput, logFilePath);
    procRunning = false;
    fqcThread.join();

    // Make checkpoint file
    sra.makeCheckpoint("fastqc2");
  }
}


void errorCorrBulk(const std::vector<SRA> & sras, std::string threads,
                   bool dispOutput, bool retainInterFiles, bool compressFiles,
                   std::string logFilePath, const INI_MAP & cfgIni) {
  logOutput("Starting error correction", logFilePath);
  std::pair<std::string, std::string> currRcorrIn;
  std::string rcorrOutDir;
  for (auto sra : sras) {
    rcorrOutDir = sra.get_sra_path_corr().first.parent_path().c_str();
    currRcorrIn.first = sra.get_sra_path_raw().first.c_str();
    currRcorrIn.second = sra.get_sra_path_raw().second.c_str();
     
    // Check for checkpoint file
    if (sra.checkpointExists("corr")) {
      logOutput("  Error-corrected version found for: " + sra.get_accession(), logFilePath);
      continue;
    }
    logOutput("  Now running error correction for:", logFilePath);
    summarize_sing_sra(sra, logFilePath, 4);

    procRunning = true;
    std::thread rcorrThread(progressAnim,2);
    run_rcorr(currRcorrIn, rcorrOutDir, threads, dispOutput, compressFiles, logFilePath);
    procRunning = false;
    rcorrThread.join();

    // Make checkpoint file
    sra.makeCheckpoint("corr");
  }
}

void remUnfixBulk(const std::vector<SRA> & sras, std::string threads, std::string ram_gb,
                  bool dispOutput, bool retainInterFiles, bool compressFiles,
                  std::string logFilePath, const INI_MAP & cfgIni) {
  logOutput("Starting post-correction removal of unfixable reads", logFilePath);
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
      logOutput("  Unfixable error-fixed version found for: " + sra.get_accession(), logFilePath);
      continue;
    }
    logOutput("  Now removing unfixable reads for:", logFilePath);
    summarize_sing_sra(sra, logFilePath, 4);

    procRunning = true;
    std::thread fixThread(progressAnim,2);
    if (sra.is_paired()) {
      rem_unfix_pe(currCorrFixIn, currCorrFixOut, ram_b, compressFiles); 
    }
    else {
      rem_unfix_se(currCorrFixIn.first, currCorrFixOut.first, ram_b);
    }
    procRunning = false;
    fixThread.join();

    // Make checkpoint file
    sra.makeCheckpoint("corr.fix");

    // If deleting intermediates, remove unfixed reads
    if (!retainInterFiles) {
      fs::remove(fs::path(sra.get_sra_path_corr().first.c_str()));
      if (sra.is_paired()) {
        fs::remove(fs::path(sra.get_sra_path_corr().second.c_str()));
      }
    }
  }
}

void trimBulk(const std::vector<SRA> & sras, std::string threads,
              bool dispOutput, bool retainInterFiles,
              std::string logFilePath, const INI_MAP & cfgIni) {
  logOutput("Starting adapter sequence trimming", logFilePath);
  INI_MAP_ENTRY cfgPipeline = cfgIni.at("Pipeline");
  std::pair<std::string, std::string> currTrimIn;
  std::pair<std::string, std::string> currTrimOutP;
  std::pair<std::string, std::string> currTrimOutU;
  for (auto sra : sras) {
    currTrimIn.first = sra.get_sra_path_corr_fix().first.c_str();
    currTrimIn.second = sra.get_sra_path_corr_fix().second.c_str();   

    if (!ini_get_bool(cfgPipeline.at("error_correction").c_str(), 0)) {
      currTrimIn.first = sra.get_sra_path_raw().first.c_str();
      currTrimIn.second = sra.get_sra_path_raw().second.c_str();
    }

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
      logOutput("  Adapter-trimmed version found for: " + sra.get_accession(), logFilePath);
      continue;
    }
    logOutput("  Running adapter sequence trimming for:", logFilePath);
    summarize_sing_sra(sra, logFilePath, 4);

    procRunning = true;
    std::thread trimThread(progressAnim,2);
    run_trimmomatic(currTrimIn, currTrimOutP, currTrimOutU, threads,
                    dispOutput, logFilePath);
    procRunning = false;
    trimThread.join();

    // Make checkpoint file
    sra.makeCheckpoint("trim");

    // If deleting intermediates, remove input files
    if (!retainInterFiles) {
      fs::remove(fs::path(currTrimIn.first));
      if (sra.is_paired()) {
        fs::remove(fs::path(currTrimIn.second));
      }
    }
  }
  // If deleting intermediates, remove directory from previous step
  if (!retainInterFiles) {
    fs::remove_all(fs::path(currTrimIn.first.c_str()).parent_path());
  }
}

void filtForeignBulk(const std::vector<SRA> & sras, std::vector<std::string> krakenDbs,
                     std::string krakenConf, std::string threads,
                     bool dispOutput, bool retainInterFiles,
                     std::string logFilePath, const INI_MAP & cfgIni) {
  logOutput("Starting foreign sequence filtering", logFilePath);
  INI_MAP_ENTRY cfgPipeline = cfgIni.at("Pipeline");
  std::pair<std::string, std::string> firstKrakIn;
  std::pair<std::string, std::string> currKrakIn;
  std::string krakOutDir;
  std::string repFile;
  std::string currKrakOut;
  std::string currRepStem;
  std::string currRepFile;

  
  for (int i = 0; i < krakenDbs.size(); i++) {
    logOutput("\nNow filtering with database: " +
              std::string(fs::path(krakenDbs[i].c_str()).filename().c_str()) + "\n",
              logFilePath);
    for (auto sra : sras) {
      // Check for checkpoint file
      if (sra.checkpointExists(std::string(fs::path(krakenDbs[i]).stem().c_str()) + ".filt")) {
        logOutput("With database: " + krakenDbs[i], logFilePath);
        logOutput("Filtered version found for: " + sra.get_accession(), logFilePath);
        continue;
      }
      logOutput("  Running filter of:", logFilePath);
      summarize_sing_sra(sra, logFilePath, 4);
      krakOutDir = sra.get_sra_path_for_filt().first.parent_path().c_str();
      if (sra.get_accession() == "") {
        repFile = krakOutDir + "/" + sra.get_file_prefix().first + "." +
                  std::string(fs::path(krakenDbs[i]).filename().c_str()) + ".report";
      }
      else {
        repFile = krakOutDir + "/" + sra.make_file_str() + "." +
                  std::string(fs::path(krakenDbs[i]).filename().c_str()) + ".report";
      }
      if (i == 0) {
        currKrakIn.first = sra.get_sra_path_trim_p().first.c_str();
        currKrakIn.second = sra.get_sra_path_trim_p().second.c_str();
        if (!ini_get_bool(cfgPipeline.at("trim_adapter_seqs").c_str(), 0)) {
          if (!ini_get_bool(cfgPipeline.at("error_correction").c_str(), 0)) {
            currKrakIn.first = sra.get_sra_path_raw().first.c_str();
            currKrakIn.second = sra.get_sra_path_raw().second.c_str();
          }
          else {
            currKrakIn.first = sra.get_sra_path_corr_fix().first.c_str();
            currKrakIn.second = sra.get_sra_path_corr_fix().second.c_str();
          }
        }
        firstKrakIn.first = currKrakIn.first;
        firstKrakIn.second = currKrakIn.second;
      }
      else {
        currKrakIn.first = sra.get_sra_path_for_filt().first.c_str();
        currKrakIn.second = sra.get_sra_path_for_filt().second.c_str();
      }
      if (sra.is_paired()) {
        currKrakOut = krakOutDir + "/TMP#.fq";
      }
        else {
        currKrakOut = krakOutDir + "TMP.fq";
      }
      procRunning = true;
      std::thread krakThread(progressAnim,2);
      run_kraken2(currKrakIn, currKrakOut, repFile, threads, krakenDbs[i],
                  krakenConf, dispOutput, logFilePath);
      procRunning = false;
      krakThread.join();
      
      if (sra.is_paired()) {
        std::rename((krakOutDir + "/TMP_1.fq").c_str(),
                    sra.get_sra_path_for_filt().first.c_str());
        std::rename((krakOutDir + "/TMP_2.fq").c_str(),
                    sra.get_sra_path_for_filt().second.c_str());
      }
      else {
        std::rename((krakOutDir + "/TMP.fq").c_str(),
                    sra.get_sra_path_for_filt().first.c_str());
      }
      // Make checkpoint file
      sra.makeCheckpoint(std::string(fs::path(krakenDbs[i]).stem().c_str()) + ".filt");

      // If deleting intermediates, remove Kraken2 input files
      if (i == 0) {
        if (!retainInterFiles) {
          fs::remove(fs::path(firstKrakIn.first));
          if (sra.is_paired()) {
            fs::remove(fs::path(firstKrakIn.second));
          }
        }
      }
    }
    // TODO: Perform summary of kraken2 filter jobs
    //
    // for each sra in sras
    //   currNumFragsFilt <- 0
    //   CurrPerFragsFilt <- 0.0
    //   totNumFragsFilt <- 0
    //   totPerFragsFilt <- 0.0
    //   for each of sra's report files
    //     currDb <- db
    //     open(report file)
    //     currNumFragsFilt <- reportFile.row(1).col(0)
    //     currPerFragsFilt <- reportFile.row(1).col(1)
    //     totNumFragsFilt += currNumFragsFilt
    //     totPerFragsFilt += currPerFragsFilt
    //     PRINT(currDb, currNumFragsFilt, currPerFragsFilt)
    //   PRINT(totNumFragsFilt, totPerFrangsFilt)
  }
  // If deleting intermediates, remove directory from previous step
  if (!retainInterFiles) {
    fs::remove_all(fs::path(firstKrakIn.first.c_str()).parent_path());
  }
}

void remOverrepBulk(const std::vector<SRA> & sras, std::string threads, std::string ram_gb,
                    bool dispOutput, bool retainInterFiles,
                    std::string logFilePath, const INI_MAP & cfgIni) {
  INI_MAP_ENTRY cfgPipeline = cfgIni.at("Pipeline");
  logOutput("Starting removal of overrepresented reads", logFilePath);
  uintmax_t ram_b = (uintmax_t)stoi(ram_gb) * 1000000000;
  std::pair<std::vector<std::string>, std::vector<std::string>> currOrepSeqsPe;
  std::vector<std::string> currOrepSeqsSe;
  std::pair<std::string, std::string> currOrepIn;
  std::pair<std::string, std::string> currOrepOut;
  fs::path fastqcDir;
  for (auto sra : sras) {
    // Check for checkpoint file
    if (sra.checkpointExists("orep.fix")) {
      logOutput("  Overrepresented-filtered version found for: " + sra.get_accession(), logFilePath);
      continue;
    }
    logOutput("  Running removal of overrepresented reads for:", logFilePath);
    summarize_sing_sra(sra, logFilePath, 4);
    currOrepIn.first = sra.get_sra_path_for_filt().first.c_str();
    currOrepIn.second = sra.get_sra_path_for_filt().second.c_str();

    fastqcDir = sra.get_fastqc_dir_2().first.parent_path();
    if (!ini_get_bool(cfgPipeline.at("filter_foreign_reads").c_str(), 0)) {
      if (!ini_get_bool(cfgPipeline.at("trim_adapter_seqs").c_str(), 0)) {
        if (!ini_get_bool(cfgPipeline.at("error_correction").c_str(), 0)) {
          currOrepIn.first = sra.get_sra_path_raw().first.c_str();
          currOrepIn.second = sra.get_sra_path_raw().second.c_str();
        }
        else {
          currOrepIn.first = sra.get_sra_path_corr_fix().first.c_str();
          currOrepIn.second = sra.get_sra_path_corr_fix().second.c_str();
        }
      }
      else {
        currOrepIn.first = sra.get_sra_path_trim_p().first.c_str();
        currOrepIn.second = sra.get_sra_path_trim_p().second.c_str();
      }
    }

    currOrepOut.first = sra.get_sra_path_orep_filt().first.c_str();
    currOrepOut.second = sra.get_sra_path_orep_filt().second.c_str();

    procRunning = true;
    std::thread orepThread(progressAnim,2);
    if (sra.is_paired()) {
      currOrepSeqsPe = get_overrep_seqs_pe(sra);
      rem_overrep_pe(currOrepIn, currOrepOut, ram_b, currOrepSeqsPe);
    }
    else {
      currOrepSeqsSe = get_overrep_seqs_se(sra);
      rem_overrep_se(currOrepIn.first, currOrepOut.first, ram_b, currOrepSeqsSe);
    }
    procRunning = false;
    orepThread.join();

    // Make checkpoint file
    sra.makeCheckpoint("orep.fix");

    // If deleting intermediates, remove input files
    if (!retainInterFiles) {
      fs::remove(currOrepIn.first);
      if (sra.is_paired()) {
        fs::remove(currOrepIn.second);
      }
    }
  }
  // If deleting intermediates, remove directory from previous step
  if (!retainInterFiles) {
    fs::remove_all(fs::path(currOrepIn.first.c_str()).parent_path());
    fs::remove_all(fs::path(fastqcDir));
  }
}

int main(int argc, char * argv[]) {
  system("setterm -cursor off");
  if (argc < 4) {
    print_help();
    return 0;
  }
  else {
    // Obtain contents of .INI configuration file
    INI_MAP cfgIni = make_ini_map(argv[1]);

    // Obtain the number of CPU threads/cores specified
    std::string threads = argv[2];

    // Obtain the amount of RAM/memory specified
    std::string ram_gb = argv[3];

    // Obtain specification for retention of intermediate files in pipeline
    bool retainInterFiles = stringToBool(argv[4]);

    // Obtain specification for verbose printing to terminal
    bool dispOutput = stringToBool(argv[5]);

    // Obtain specification for compression of output files
    bool compressFiles = stringToBool(argv[6]);

    // Obtain path to log file from config file
    std::string logFilePath = cfgIni["General"]["log_file"];
 
    // Make project file structure
    make_proj_space(cfgIni, "preprocess");

    // Declare vector for SRA objects
    std::vector<SRA> sras;

    // Declare vector for local run paths
    std::vector<std::string> localDataFiles;

    // Create vector of SRA objects from SRA accessions, using NCBI web API
    sras = get_sras(cfgIni, compressFiles);

    // Obtain terminal window size for printing purposes
    struct winsize w;
    ioctl(STDOUT_FILENO, TIOCGWINSZ, &w);
    
    if (!sras.empty()) {
      retrieveSraData(sras, threads, dispOutput, compressFiles, retainInterFiles, logFilePath);
      logOutput("Successfully downloaded sequence data", logFilePath);
    }
    // Get single/paired filenames of local data
    for (auto fqFileName : cfgIni.at("Local files")) {
      localDataFiles.push_back(fqFileName.first);
    }
    std::pair<std::string, std::string> sraRunsLocal;
    std::string localDataDir = cfgIni["General"]["local_data_directory"];
    if (localDataDir[localDataDir.size() - 1] != '/') {
      localDataDir += "/";
    }
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
      if (fs::exists(localDataDir + sraRunsLocal.first) &&
          fs::exists(localDataDir + sraRunsLocal.second)) {
        sras.push_back(SRA(sraRunsLocal.first, sraRunsLocal.second, cfgIni, compressFiles));
      }
      else {
        if (sraRunsLocal.first != "" &&
            !fs::exists(localDataDir + sraRunsLocal.first)) {
          std::cout << "ERROR: Local run not found: \"" << sraRunsLocal.first << "\""
                    << std::endl;
        }
        if (sraRunsLocal.second != "" &&
            !fs::exists(localDataDir + sraRunsLocal.second)) {
          std::cout << "ERROR: Local run not found: \"" << sraRunsLocal.second << "\""
                    << std::endl;
        }
      }
    }
    
    // Check if no SRAs specified
    if (sras.empty()) {
      std::cout << "ERROR: No SRA runs specified. Please check config file" << std::endl;
    }

    logOutput("Raw sequence data prepared for pre-assembly processing", logFilePath);
    logOutput("Initiating preprocessing of data ...", logFilePath);

    // Summarize preprocess task
    logOutput("Paando Preprocess started with following parameters:", logFilePath);
    logOutput("  Config file:     " + std::string(argv[1]), logFilePath);
    logOutput("  Threads (Cores): " + threads, logFilePath);
    logOutput("  Memory (GB):     " + ram_gb, logFilePath);
    logOutput("  SRA runs:\n", logFilePath); 
    summarize_all_sras(sras, logFilePath, 6);
    
    std::string retainStr;
    if (retainInterFiles) {
      retainStr = "YES";
    }
    else {
      retainStr = "NO";
    }
    logOutput("  Retain intermediate files: " + retainStr, logFilePath);

    std::string compressStr;
    if (compressFiles) {
      compressStr = "YES";
    }
    else {
      compressStr = "NO";
    }
    logOutput("  Compress all output files: " + compressStr, logFilePath);

    std::string fastqc_dir_1(sras[0].get_fastqc_dir_1().first.parent_path().parent_path().c_str());
    std::string fastqc_dir_2(sras[0].get_fastqc_dir_2().first.parent_path().parent_path().c_str());
    
    INI_MAP_ENTRY cfgPipeline = cfgIni.at("Pipeline");


    // Run initial fastqc on reads
    if (ini_get_bool(cfgPipeline.at("pre_quality_check").c_str(), 0)) {
      fastqcBulk1(sras, threads, dispOutput, logFilePath);
    }

    // Error-correction stage
    if (ini_get_bool(cfgPipeline.at("error_correction").c_str(), 0)) {
      errorCorrBulk(sras, threads, dispOutput, retainInterFiles, compressFiles, logFilePath, cfgIni);
    
      // Remove reads with unfixable errors
      remUnfixBulk(sras, threads, ram_gb, dispOutput, retainInterFiles, compressFiles, logFilePath, cfgIni);
    }

    // Adapter sequence trimming stage
    if (ini_get_bool(cfgPipeline.at("trim_adapter_seqs").c_str(), 0)) {
      trimBulk(sras, threads, dispOutput, retainInterFiles, logFilePath, cfgIni);
    }

    // Run kraken2 to remove foreign reads
    if (ini_get_bool(cfgPipeline.at("filter_foreign_reads").c_str(), 0)) {
      std::vector<std::string> krakenDbs = get_kraken2_dbs(cfgIni);
      std::string krakenConf = get_kraken2_conf(cfgIni);
      filtForeignBulk(sras, krakenDbs, krakenConf, threads,
                      dispOutput, retainInterFiles, logFilePath, cfgIni);
    }

    // Run fastqc on all runs
    fastqcBulk2(sras, threads, dispOutput, logFilePath, cfgIni);

    // Remove reads with over-represented sequences
    if (ini_get_bool(cfgPipeline.at("remove_overrepresented").c_str(), 0)) {
      remOverrepBulk(sras, threads, ram_gb, dispOutput, retainInterFiles, logFilePath, cfgIni);
    }
  }
  system("setterm -cursor on");
  return 0;
}
