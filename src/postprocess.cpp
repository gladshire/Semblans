#include "postprocess.h"


extern std::vector<std::string> stepDirs;

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
  logOutput("Starting BLASTX alignment", logFilePath);
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
    logOutput("  Now running BLASTX alignment for:", logFilePath);
    summarize_sing_trans(trans, logFilePath, 4);
    currTransInDiam = trans.get_trans_path_trinity().c_str();
    // Check if BlastX checkpoint exists
    if (trans.checkpointExists("blastx")) {
      logOutput("  BLASTX checkpoint found for: " + trans.get_file_prefix(), logFilePath);
      continue;
    }
    makeDb(refProt, blastDbDir, dispOutput, logFilePath);
    // Run BlastX
    currBlastDbName = std::string(fs::path(refProt.c_str()).stem().c_str());
    blastxDiam(currTransInDiam, blastDbDir + currBlastDbName, threads,
               blastDbDir, dispOutput, logFilePath);
    // Create BlastX checkpoint
    trans.makeCheckpoint("blastx");
  }
}

void remChimeraBulk(const std::vector<transcript> & transVec, std::string ram_gb,
                    std::string logFilePath) {
  logOutput("Starting chimera removal", logFilePath);
  std::string currTransInChim;
  std::string currTransOutChim;
  std::string currBlastx;
  std::string currTransInfo;
  std::string currTransCut;
  std::string chimOutDir;
  for (auto trans : transVec) {
    logOutput("  Now running chimera removal for:", logFilePath);
    summarize_sing_trans(trans, logFilePath, 4);
    // Check if chimera removal checkpoint exists
    if (trans.checkpointExists("chimera")) {
      logOutput("  Chimera removal checkpoint found for: " + trans.get_file_prefix(), logFilePath);
      continue;
    }
    currTransInChim = trans.get_trans_path_trinity().c_str();
    currTransOutChim = trans.get_trans_path_chimera().c_str();
    currBlastx = trans.get_trans_path_blastx().c_str();
    currTransInfo = trans.get_trans_path_cinfo().c_str();
    currTransCut = trans.get_trans_path_ccut().c_str();
    chimOutDir = trans.get_trans_path_chimera().parent_path().c_str();

    detect_chimera(currBlastx, chimOutDir);
    removeChimera(currTransInChim, currTransOutChim, currTransInfo, currTransCut, ram_gb,
                  logFilePath);
    // Create chimera removal checkpoint
    trans.makeCheckpoint("chimera");
  }
}

void salmonBulk(const std::vector<transcript> & transVec, std::string threads,
                bool dispOutput, std::string logFilePath) {
  logOutput("Starting quantification of transcripts", logFilePath);
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
    logOutput("  Now performing quantification for:", logFilePath);
    summarize_sing_trans(trans, logFilePath, 4);
    // Check if clustering checkpoint exists

    currTransInSalm = trans.get_trans_path_chimera().c_str();
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
    salmon_index(currTransInSalm, currIndex, threads, dispOutput, logFilePath);
    // Create salmon index checkpoint
    trans.makeCheckpoint("index");

    // Check if salmon quant checkpoint exists
    if (trans.checkpointExists("quant")) {
      logOutput("  Quant checkpoint found for: " + trans.get_file_prefix(), logFilePath);
      continue;
    }
    // Perform salmon quant of transcript
    salmon_quant(currTransInSalm, currIndex, currQuant, currSraRunsIn, threads,
                 dispOutput, logFilePath);
    // Create salmon quant checkpoint
    trans.makeCheckpoint("quant");
  }
}

void corsetBulk(const std::vector<transcript> & transVec, std::string ram_gb,
                bool dispOutput, std::string logFilePath) {
  logOutput("Removing redundant transcripts for:", logFilePath);
  std::string currTransInCors;
  std::string currTransPrefix;
  std::string currEqClassFile;
  std::string currTransClust;
  std::string currTransLargestClust;
  std::string currTransRedund;
  std::string currOutDir;
  uintmax_t ram_b = (uintmax_t)stoi(ram_gb) * 1000000000;
  for (auto trans : transVec) {
    logOutput("  Now performing cluster-based filtering for:", logFilePath);
    summarize_sing_trans(trans, logFilePath, 4);
    currTransInCors = trans.get_trans_path_chimera().c_str();
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
    // Check if corset run checkpoint exists
    if (trans.checkpointExists("corset")) {
      logOutput("  Corset checkpoint found for: " + trans.get_file_prefix(), logFilePath);
      continue;
    }
    // Perform corset run
    corset_eq_classes(trans.get_file_prefix(), currEqClassFile, currOutDir, dispOutput, logFilePath);
    // Create corset run checkpoint
    trans.makeCheckpoint("corset");

    // Check if corset filter checkpoint exists
    if (trans.checkpointExists("cors_filt")) {
      logOutput("  Corset filtering checkpoint found for: " + trans.get_file_prefix(), logFilePath);
      continue;
    }
    // Filter corset output
    filterCorset(currTransInCors, currTransClust, currTransLargestClust, currTransRedund,
                 ram_b, currOutDir, logFilePath);
    // Create corset filtering checkpoint
    trans.makeCheckpoint("cors_filt");
  }
}

void transdecBulk(const std::vector<transcript> & transVec, std::string threads,
                  std::string ram_gb, bool dispOutput, std::string logFilePath,
                  const INI_MAP & cfgIni) {
  logOutput("Starting generation of coding sequences / peptides", logFilePath);
  std::string currTransInTD;
  std::string currTransCds;
  std::string currTransPep;
  std::string currDb;
  std::string currOutDirTD;
  std::string refProt = cfgIni.at("General").at("reference_proteome_path");
  if (!fs::exists(fs::path(refProt.c_str()))) {
    logOutput("ERROR: Reference proteome: " + refProt + " not found", logFilePath);
    logOutput("  Please ensure its path is correctly specified in your config INI file", logFilePath);
    exit(1);
  }
  std::string blastDbDir = cfgIni.at("General").at("output_directory") + "/" +
                           cfgIni.at("General").at("project_name") + "/" +
                           stepDirs[8] + "/";
  std::string blastDbName(fs::path(refProt.c_str()).stem().c_str());
  uintmax_t ram_b = (uintmax_t)stoi(ram_gb) * 1000000000;
  for (auto trans : transVec) {
    logOutput("  Now building coding sequences / peptides for:", logFilePath);
    summarize_sing_trans(trans, logFilePath, 4);
    currTransInTD = trans.get_trans_path_largest().c_str();
    currTransCds = trans.get_trans_path_cds().c_str();
    currTransPep = trans.get_trans_path_prot().c_str();
    
    currDb = blastDbDir + "/" + blastDbName;
    currOutDirTD = trans.get_trans_path_cds().parent_path().c_str();

    // Check if transdecoder checkpoint exists
    if (trans.checkpointExists("transdecoder")) {
      logOutput("  TransDecoder checkpoint found for: " + trans.get_file_prefix(), logFilePath);
      continue;
    }
    // Perform transdecoder run to obtain coding sequences / peptides
    run_transdecoder(currTransInTD, currTransCds, currTransPep, threads, ram_b,
                     currDb, currOutDirTD, dispOutput, logFilePath);
    // Create transdecoder checkpoint
    trans.makeCheckpoint("transdecoder");
  }
}

void annotateBulk(const std::vector<transcript> & transVec, std::string threads,
                  std::string ram_gb, bool dispOutput, std::string logFilePath,
                  const INI_MAP & cfgIni) {
  logOutput("Starting annotation of transcripts", logFilePath);
  std::string currTransIn;
  std::string currTransPep;
  std::string currTransOut;
  // Check for checkpoint

  std::string email = "gladshire@gmail.com";
  for (auto trans : transVec) {
    logOutput("  Now running annotation on:", logFilePath);
    summarize_sing_trans(trans, logFilePath, 4);

    currTransIn = trans.get_trans_path_cds().c_str();
    currTransPep = trans.get_trans_path_prot().c_str();
  
    currTransOut = trans.get_trans_path_annot().c_str();


    annotateTranscript(currTransIn, currTransPep, currTransOut, threads, ram_gb, logFilePath, email);
  }
}

int main(int argc, char * argv[]) {
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


    // BlastX alignment of transcript to reference proteome
    blastxDiamBulk(transVec, threads, dispOutput, logFilePath, cfgIni);

    // Detect and remove chimeric transcripts
    remChimeraBulk(transVec, ram_gb, logFilePath);

    // Perform salmon index of transcripts
    salmonBulk(transVec, threads, dispOutput, logFilePath);
   
    // Perform corset run to cluster transcripts
    corsetBulk(transVec, ram_gb, dispOutput, logFilePath);
    
    // Run transdecoder
    transdecBulk(transVec, threads, ram_gb, dispOutput, logFilePath, cfgIni);
  
    // Annotate transcriptome
    //annotateBulk(transVec, threads, ram_gb, dispOutput, logFilePath, cfgIni);
  }
  else {
  
  }
}
