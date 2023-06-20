#include "postprocess.h"

std::atomic<bool> procRunning(false);
extern std::vector<std::string> stepDirs;

void progressAnim(int numSpace) {
  const std::string anim[] = {".  ", ".. ", "..."};
  int animIndex = 0;

  while (procRunning) {
    std::cout << "\r";
    for (int i = 0; i < numSpace; i++) {
      std::cout << " ";
    }
    std::cout << anim[animIndex] << std::flush;
    animIndex = (animIndex +1) % 3;
    std::this_thread::sleep_for(std::chrono::milliseconds(1000));
  }
  std::cout << "\r";
  for (int i = 0; i < numSpace; i++) {
    std::cout << " ";
  }
  std::cout << "   " << std::endl;
}

std::vector<transcript> get_transcript(std::vector<SRA> sras) {
  std::vector<transcript> transcripts;
  for (auto & sra : sras) {
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

void blastxDiamBulk(const std::vector<transcript> & transVec, std::string threads,
                    bool dispOutput, std::string logFilePath, const INI_MAP & cfgIni) {
  logOutput("\nStarting BLASTX alignment", logFilePath);
  std::string currTransInDiam;
  std::string currBlastDbName;
  std::string refProt = cfgIni.at("General").at("reference_proteome_path");
  if (!fs::exists(fs::path(refProt.c_str()))) {
    logOutput("ERROR: Reference proteome: " + refProt + " not found", logFilePath);
    logOutput("  Please ensure its path is correctly specified in your config INI file", logFilePath);
    exit(1);
  }
  std::string blastDbDir;
  blastDbDir = cfgIni.at("General").at("output_directory") + "/" +
               cfgIni.at("General").at("project_name") + "/" +
               stepDirs[8] + "/";
  std::string blastDbName;
  blastDbName = std::string(fs::path(refProt.c_str()).stem().c_str());
  for (auto trans : transVec) {
    // Check if BlastX checkpoint exists
    if (trans.checkpointExists("blastx")) {
      logOutput("  BLASTX checkpoint found for: " + trans.get_file_prefix(), logFilePath);
      continue;
    }
    logOutput("  Now running BLASTX alignment for:", logFilePath);
    summarize_sing_trans(trans, logFilePath, 4);
    currTransInDiam = trans.get_trans_path_trinity().c_str();
    makeDb(refProt, blastDbDir, dispOutput, logFilePath);
    // Run BlastX
    currBlastDbName = std::string(fs::path(refProt.c_str()).stem().c_str());
    
    procRunning = true;
    std::thread blastxThread(progressAnim, 2);
    blastxDiam(currTransInDiam, blastDbDir + currBlastDbName, threads,
               blastDbDir, dispOutput, logFilePath);
    procRunning = false;
    blastxThread.join();

    // Create BlastX checkpoint
    trans.makeCheckpoint("blastx");
  }
}

void remChimeraBulk(const std::vector<transcript> & transVec, std::string ram_gb,
                    std::string logFilePath, const INI_MAP & cfgIni) {
  logOutput("\nStarting chimera removal", logFilePath);
  std::string currTransInChim;
  std::string currTransOutChim;
  std::string currBlastx;
  std::string currTransInfo;
  std::string currTransCut;
  std::string chimOutDir;
  for (auto trans : transVec) {
    // Check if chimera removal checkpoint exists 
    if (trans.checkpointExists("chim.fix")) {
      logOutput("  Chimera removal checkpoint found for: " + trans.get_file_prefix(),
                logFilePath);
      continue;
    }
    logOutput("  Now running chimera removal for:", logFilePath);
    summarize_sing_trans(trans, logFilePath, 4);

   
    currTransInChim = trans.get_trans_path_trinity().c_str();
    currTransOutChim = trans.get_trans_path_chimera().c_str();
    currBlastx = trans.get_trans_path_blastx().c_str();
    currTransInfo = trans.get_trans_path_cinfo().c_str();
    currTransCut = trans.get_trans_path_ccut().c_str();
    chimOutDir = trans.get_trans_path_chimera().parent_path().c_str();

    procRunning = true;
    std::thread chimeraThread(progressAnim, 2);
    detect_chimera(currBlastx, chimOutDir);
    removeChimera(currTransInChim, currTransOutChim, currTransInfo, currTransCut, ram_gb,
                  logFilePath);
    procRunning = false;
    chimeraThread.join();

    // Create chimera removal checkpoint
    trans.makeCheckpoint("chim.fix");
  }
}

void salmonBulk(const std::vector<transcript> & transVec, std::string threads,
                bool dispOutput, std::string logFilePath, const INI_MAP & cfgIni) {
  logOutput("\nStarting transcript quantification", logFilePath);
  INI_MAP_ENTRY cfgPipeline = cfgIni.at("Pipeline");
  std::string currTransInSalm;
  std::string currIndex;
  std::string currQuant;
  std::string currTransInfoFileStr;
  std::ifstream currTransInfoFile;
  std::string currLineInfo;
  size_t spacePos;
  size_t nlPos;
  std::vector<std::pair<std::string, std::string>> currSraRunsIn;
  std::pair<std::string, std::string> currSraRun;
  for (auto trans : transVec) {
    currTransInSalm = trans.get_trans_path_chimera().c_str();
    if (!ini_get_bool(cfgPipeline.at("remove_chimera_reads").c_str(), 0)) {
      currTransInSalm = trans.get_trans_path_trinity().c_str();
    }
    currIndex = trans.get_trans_path_index().c_str();
    currQuant = trans.get_trans_path_quant().c_str();
    
    currTransInfoFileStr = std::string(trans.get_trans_path_trinity().replace_extension(".ti").c_str());

    currTransInfoFile.open(currTransInfoFileStr);
    while(getline(currTransInfoFile, currLineInfo)) {
      spacePos = currLineInfo.find(" ");
      nlPos = currLineInfo.find("\n");
      // If line corresponds to paired
      if (spacePos != std::string::npos) {
        // Assign everything up to space as currSraRun.first
        currSraRun.first = currLineInfo.substr(0, spacePos + 1);
        // Assign everything from space to endline as currSraRun.second
        currSraRun.second = currLineInfo.substr(spacePos + 1, nlPos - spacePos);
      }
      else {
        currSraRun.first = currLineInfo;
      }
      currSraRunsIn.push_back(currSraRun);
    }

    // Check if salmon index checkpoint exists
    if (trans.checkpointExists("index")) {
      logOutput("  Indexing checkpoint found for: " + trans.get_file_prefix(), logFilePath);
      continue;
    }
    // Perform salmon index of transcript
    logOutput("  Now mapping reads to transcripts: ", logFilePath);
    summarize_sing_trans(trans, logFilePath, 4);

    procRunning = true;
    std::thread salmIdxThread(progressAnim, 2);
    salmon_index(currTransInSalm, currIndex, threads, dispOutput, logFilePath);
    procRunning = false;
    salmIdxThread.join();

    // Create salmon index checkpoint
    trans.makeCheckpoint("index");

    // Check if salmon quant checkpoint exists
    if (trans.checkpointExists("quant")) {
      logOutput("  Quant checkpoint found for: " + trans.get_file_prefix(), logFilePath);
      continue;
    }
    // Perform salmon quant of transcript
    logOutput("  Now quantifying transcripts: ", logFilePath);
    summarize_sing_trans(trans, logFilePath, 4);

    procRunning = true;
    std::thread salmQntThread(progressAnim, 2);
    salmon_quant(currTransInSalm, currIndex, currQuant, currSraRunsIn, threads,
                 dispOutput, logFilePath);
    procRunning = false;
    salmQntThread.join();

    // Create salmon quant checkpoint
    trans.makeCheckpoint("quant");
  }
}

void corsetBulk(const std::vector<transcript> & transVec, std::string ram_gb,
                bool dispOutput, std::string logFilePath, const INI_MAP & cfgIni) {
  logOutput("\nRemoving redundant transcripts", logFilePath);
  INI_MAP_ENTRY cfgPipeline = cfgIni.at("Pipeline");
  std::string currTransInCors;
  std::string currTransPrefix;
  std::string currEqClassFile;
  std::string currTransClust;
  std::string currTransLargestClust;
  std::string currTransRedund;
  std::string currOutDir;
  uintmax_t ram_b = (uintmax_t)stoi(ram_gb) * 1000000000;
  for (auto trans : transVec) {
    // Check if corset run checkpoint exists
    if (trans.checkpointExists("clust")) {
      logOutput("  Corset checkpoint found for: " + trans.get_file_prefix(), logFilePath);
      continue;
    }
    logOutput("  Now clustering transcripts into genes:", logFilePath);
    summarize_sing_trans(trans, logFilePath, 4);
    currTransInCors = trans.get_trans_path_chimera().c_str();
    if (!ini_get_bool(cfgPipeline.at("remove_chimera_reads").c_str(), 0)) {
      currTransInCors = trans.get_trans_path_trinity().c_str();
    }
    if (trans.get_org_name() == "") {
      currTransPrefix = trans.get_file_prefix();
    }
    else {
      currTransPrefix = trans.make_file_str();
    }
    currEqClassFile = std::string((trans.get_trans_path_quant() /
                                   fs::path("aux_info") /
                                   fs::path("eq_classes.txt.gz")).c_str());
    currTransClust = trans.get_trans_path_clust().c_str();
    currTransLargestClust = trans.get_trans_path_largest().c_str();
    currTransRedund = trans.get_trans_path_redund().c_str();
    currOutDir = trans.get_trans_path_clust().parent_path().c_str();
   
    // Perform corset run
    procRunning = true;
    std::thread clustThread(progressAnim, 2);
    corset_eq_classes(trans.get_file_prefix(), currEqClassFile, currOutDir,
                      dispOutput, logFilePath);
    procRunning = false;
    clustThread.join();

    // Create corset run checkpoint
    trans.makeCheckpoint("clust");

    // Check if corset filter checkpoint exists
    if (trans.checkpointExists("redund.fix")) {
      logOutput("  Corset filtering checkpoint found for: " + trans.get_file_prefix(),
      logFilePath);
      continue;
    }
    // Filter corset output
    logOutput("  Now filtering redundant transcripts:", logFilePath);
    summarize_sing_trans(trans, logFilePath, 4);

    procRunning = true;
    std::thread clustFiltThread(progressAnim, 2);
    filterCorset(currTransInCors, currTransClust, currTransLargestClust, currTransRedund,
                 ram_b, currOutDir, logFilePath);
    procRunning = false;
    clustFiltThread.join();

    // Create corset filtering checkpoint
    trans.makeCheckpoint("redund.fix");
  }
}

void transdecBulk(const std::vector<transcript> & transVec, std::string threads,
                  std::string ram_gb, bool dispOutput, std::string logFilePath,
                  const INI_MAP & cfgIni) {
  logOutput("\nStarting prediction of coding regions", logFilePath);
  INI_MAP_ENTRY cfgPipeline = cfgIni.at("Pipeline");
  std::string currTransInTD;
  std::string currTransCds;
  std::string currTransPep;
  std::string currDb;
  std::string currOutDirTD;
  std::string refProt = cfgIni.at("General").at("reference_proteome_path");
  if (!fs::exists(fs::path(refProt.c_str()))) {
    logOutput("ERROR: Reference proteome: " + refProt + " not found", logFilePath);
    logOutput("  Please ensure its path is correctly specified in your config INI file",
              logFilePath);
    exit(1);
  }
  std::string blastDbDir = cfgIni.at("General").at("output_directory") + "/" +
                           cfgIni.at("General").at("project_name") + "/" +
                           stepDirs[8] + "/";
  std::string blastDbName(fs::path(refProt.c_str()).stem().c_str());
  uintmax_t ram_b = (uintmax_t)stoi(ram_gb) * 1000000000;
  for (auto trans : transVec) {
    // Check if coding region prediction checkpoint exists
    if (trans.checkpointExists("cdr.predict")) {
      logOutput("  Coding region prediction checkpoint found for: " + trans.get_file_prefix(),
                logFilePath);
      continue;
    }
    logOutput("  Now predicting coding regions for:", logFilePath);
    summarize_sing_trans(trans, logFilePath, 4);
    currTransInTD = trans.get_trans_path_largest().c_str();
    if (!ini_get_bool(cfgPipeline.at("cluster_filtering").c_str(), 0)) {
      if (!ini_get_bool(cfgPipeline.at("remove_chimera_reads").c_str(), 0)) {
        currTransInTD = trans.get_trans_path_trinity().c_str();
      }
      else {
        currTransInTD = trans.get_trans_path_chimera().c_str();
      }
    }
    currTransCds = trans.get_trans_path_cds().c_str();
    currTransPep = trans.get_trans_path_prot().c_str();
    
    currDb = blastDbDir + "/" + blastDbName;
    currOutDirTD = trans.get_trans_path_cds().parent_path().c_str();
       
    // Perform transdecoder run to obtain coding sequences / peptides
    procRunning = true;
    std::thread transDecThread(progressAnim, 2);
    run_transdecoder(currTransInTD, currTransCds, currTransPep, threads, ram_b,
                     currDb, currOutDirTD, dispOutput, logFilePath);
    // Create coding region prediction checkpoint
    trans.makeCheckpoint("cdr.predict");
    procRunning = false;
    transDecThread.join();
  }
}

void annotateBulk(const std::vector<transcript> & transVec, std::string threads,
                  std::string ram_gb, bool dispOutput, std::string logFilePath,
                  const INI_MAP & cfgIni) {
  logOutput("\nStarting annotation of transcripts", logFilePath);
  std::string currTransIn;
  std::string currTransPep;
  std::string currTransOut;
  std::string email = "gladshire@gmail.com";
  for (auto trans : transVec) {
    // Check if annotation checkpoint exists 
    if (trans.checkpointExists("annotate")) {
      logOutput("  Annotation checkpoint found for: " + trans.get_file_prefix(), logFilePath);
      continue;
    }
    logOutput("  Now running annotation on:", logFilePath);
    summarize_sing_trans(trans, logFilePath, 4);

    currTransIn = trans.get_trans_path_cds().c_str();
    currTransPep = trans.get_trans_path_prot().c_str();
  
    currTransOut = trans.get_trans_path_annot().c_str();

    // Perform annotation of transcript
    procRunning = true;
    std::thread annotThread(progressAnim, 2);
    annotateTranscript(currTransIn, currTransPep, currTransOut,
                       threads, ram_gb, logFilePath, email);
    procRunning = false;
    annotThread.join();

    // Create annotation checkpoint
    trans.makeCheckpoint("annotate");
  }
}

int main(int argc, char * argv[]) {
  system("setterm -cursor off");
  if (argc > 1) {
    // Get INI config file
    INI_MAP cfgIni = make_ini_map(argv[1]);

    // Get number of threads
    std::string threads = argv[2];
    
    // Get RAM in GB
    std::string ram_gb = argv[3];
    
    // Determine whether to keep intermediate files
    bool retainInterFiles = stringToBool(argv[4]);
    
    // Determine whether to print output of programs
    bool dispOutput = stringToBool(argv[5]);
    
    // Determine whether to compress files
    //bool compressFiles = ini_get_bool(cfgIni["General"]["compress_files"].c_str(), 0);
    bool compressFiles = false;
    
    // Retrieve SRA objects, convert to transcripts
    make_proj_space(cfgIni, "postprocess");
    std::string outputDir = cfgIni["General"]["output_directory"];
    std::string projDir = outputDir + cfgIni["General"]["project_name"] + "/";
    std::string refProt = cfgIni["General"]["reference_proteome_path"];

    // Obtain path to log file from config file
    fs::path logFile(cfgIni["General"]["log_file"].c_str());
    std::string logFilePath;
    if (logFile.filename() == logFile) {
      logFilePath = std::string(fs::canonical((fs::path(cfgIni["General"]["output_directory"].c_str()) /
                                               fs::path(cfgIni["General"]["log_file"].c_str()))).c_str());
    }
    else {
      logFilePath = std::string(fs::canonical((fs::path(cfgIni["General"]["log_file"].c_str()))).c_str());
    }
    
    std::vector<SRA> sras;
    std::vector<std::string> localDataFiles;
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
    if (sras.empty()) {
      std::cout << "ERROR: No SRA runs specified. Please check config file" << std::endl;
    }
    
    std::vector<transcript> transVec;
    fs::path currFile = transcript(sras[0]).get_trans_path_trinity().parent_path();
    fs::directory_iterator fileIter{currFile};
    while (fileIter != fs::directory_iterator{}) {
      // Iterate through files in transcript directory
      // Push transcript to transcript vector
      if (fileIter->path().extension() == ".fasta") {
        transcript currTrans(fileIter->path().c_str(), cfgIni);
        transVec.push_back(currTrans);
      }
      fileIter++;
    }
 

  
    // Summarize program execution parameters
    logOutput("Semblans Postprocess started with following parameters:", logFilePath);
    logOutput("  Config file:     " + std::string(argv[1]), logFilePath);
    logOutput("  Threads (Cores): " + threads, logFilePath);
    logOutput("  Memory (GB):     " + ram_gb, logFilePath);
    logOutput("  Reference Prot:  " + refProt, logFilePath);
    logOutput("  SRA runs:\n", logFilePath);
    summarize_all_sras(sras, logFilePath, 6);

    INI_MAP_ENTRY cfgPipeline = cfgIni.at("Pipeline");
 
    if (ini_get_bool(cfgPipeline.at("remove_chimera_reads").c_str(), 0)) {
      // BlastX alignment of transcript to reference proteome
      blastxDiamBulk(transVec, threads, dispOutput, logFilePath, cfgIni);

      // Detect and remove chimeric transcripts
      remChimeraBulk(transVec, ram_gb, logFilePath, cfgIni);
    }

    if (ini_get_bool(cfgPipeline.at("cluster_filtering").c_str(), 0)) {
      // Perform salmon index of transcripts
      salmonBulk(transVec, threads, dispOutput, logFilePath, cfgIni);
   
      // Perform corset run to cluster transcripts
      corsetBulk(transVec, ram_gb, dispOutput, logFilePath, cfgIni);
    }
    
    // Run transdecoder
    transdecBulk(transVec, threads, ram_gb, dispOutput, logFilePath, cfgIni);
  
    // Annotate transcriptome
    annotateBulk(transVec, threads, ram_gb, dispOutput, logFilePath, cfgIni);
  }
  else {
  
  }
  system("setterm -cursor on");
}
