// TODO: Implement removal of intermediates for postprocess

#include "postprocess.h"

std::atomic<bool> procRunning(false);
extern std::vector<std::string> stepDirs;


// Given a vector of SRAs, instantiate transcript objects constructed
// off of them
std::vector<transcript> get_transcript(std::vector<SRA> sras) {
  std::vector<transcript> transcripts;
  for (auto & sra : sras) {
    transcript currTrans(sra);
    transcripts.push_back(currTrans);
  }
  return transcripts;
}

// Given a "true" / "false" string, return the appropriate boolean
bool stringToBool(std::string boolStr) {
  bool boolConv;
  for (int i = 0; i < boolStr.length(); i++) {
    boolStr[i] = std::tolower(boolStr[i]);
  }
  boolConv = (boolStr == "true") ? true : false;
  return boolConv;
}

// Given a vector of transcript objects, perform a blastx alignment against the
// user-provided reference proteome
void blastxBulk(const std::vector<transcript> & transVec, std::string threads,
                bool dispOutput, std::string logFilePath,
                std::string refProt, std::string outDir,
                const INI_MAP & cfgIni = {}) {
  logOutput("\nStarting BLASTX alignment\n", logFilePath);
  std::string currTransIn;
  std::string currBlastDbName;
  INI_MAP_ENTRY cfgIniGen;
  INI_MAP_ENTRY cfgIniAlign;
  bool useBlast;
  std::string maxEvalue;
  std::string maxTargetSeqs;
  std::string refProteome;
  std::string blastDbDir;
  std::string blastDbName;
  if (!cfgIni.empty()) {
    cfgIniGen = cfgIni.at("General");
    cfgIniAlign = cfgIni.at("Alignment settings");
    useBlast = ini_get_bool(cfgIniAlign.at("use_blast_instead_of_diamond").c_str(), 0);
    maxEvalue = cfgIniAlign.at("blastx_max_evalue");
    maxTargetSeqs = cfgIniAlign.at("blastx_max_target_seqs");
    refProteome = cfgIniGen.at("reference_proteome_path");
    if (!fs::exists(fs::path(refProt.c_str()))) {
      logOutput("ERROR: Reference proteome: " + refProteome + " not found\n", logFilePath);
      logOutput("  Please ensure its path is correctly specified in your config file\n",
                logFilePath);
      exit(1);
    }
    blastDbDir = cfgIniGen.at("output_directory") + "/" +
                 cfgIniGen.at("project_name") + "/" +
                 stepDirs[8] + "/";
    blastDbName = std::string(fs::path(refProteome.c_str()).stem().c_str());
  }
  else {
    useBlast = false;
    maxEvalue = "0.001";
    maxTargetSeqs = "25";
    refProteome = refProt;
    if (!fs::exists(fs::path(refProt.c_str()))) {
      logOutput("ERROR: Reference proteome: " + refProteome + " not found\n", logFilePath);
      exit(1);
    }
    blastDbDir = outDir + stepDirs[0] + "/";
    blastDbName = std::string(fs::path(refProteome.c_str()).stem().c_str());
  }

  /*
  std::string currTransIn;
  std::string currBlastDbName;
  INI_MAP_ENTRY cfgGen = cfgIni.at("General");
  INI_MAP_ENTRY cfgAlign = cfgIni.at("Alignment settings");
  bool useBlast = ini_get_bool(cfgAlign.at("use_blast_instead_of_diamond").c_str(), 0);
  std::string maxEvalue = cfgAlign.at("blastx_max_evalue");
  std::string maxTargetSeqs = cfgAlign.at("blastx_max_target_seqs");
  std::string refProt = cfgGen.at("reference_proteome_path");
  if (!fs::exists(fs::path(refProt.c_str()))) {
    logOutput("ERROR: Reference proteome: " + refProt + " not found\n", logFilePath);
    logOutput("  Please ensure its path is correctly specified in your config INI file\n",
              logFilePath);
    exit(1);
  }
  std::string blastDbDir;
  blastDbDir = cfgGen.at("output_directory") + "/" +
               cfgGen.at("project_name") + "/" +
               stepDirs[8] + "/";
  std::string blastDbName;
  blastDbName = std::string(fs::path(refProt.c_str()).stem().c_str());
  */
  for (auto trans : transVec) {
    // Check if BlastX checkpoint exists
    if (trans.checkpointExists("blastx")) {
      logOutput("  BLASTX checkpoint found for: " +
                trans.get_file_prefix() + "\n", logFilePath);
      continue;
    }
    logOutput("  Now running BLASTX alignment for:\n", logFilePath);
    summarize_sing_trans(trans, logFilePath, 4);
    currTransIn = trans.get_trans_path_trinity().c_str();
    if (useBlast) {
      makeBlastDb(refProt, blastDbDir, dispOutput, logFilePath);
    }
    else {
      makeBlastDbDiam(refProt, blastDbDir, dispOutput, logFilePath);
    }
    // Run BlastX
    currBlastDbName = std::string(fs::path(refProt.c_str()).stem().c_str());

    if (!dispOutput) {
      procRunning = true;
      std::thread blastxThread(progressAnim, "  ", logFilePath);
      if (useBlast) {
        blastxDiam(currTransIn, blastDbDir + currBlastDbName,
                   maxEvalue, maxTargetSeqs, threads,
                   blastDbDir, dispOutput, logFilePath);
      }
      else {
        blastx(currTransIn, blastDbDir + currBlastDbName,
               maxEvalue, maxTargetSeqs, threads,
               blastDbDir, dispOutput, logFilePath);
      }
      procRunning = false;
      blastxThread.join();
    }
    else {
       if (useBlast) {
         blastx(currTransIn, blastDbDir + currBlastDbName,
                maxEvalue, maxTargetSeqs, threads,
                blastDbDir, dispOutput, logFilePath);
       }
       else {
         blastxDiam(currTransIn, blastDbDir + currBlastDbName,
                    maxEvalue, maxTargetSeqs, threads,
                    blastDbDir, dispOutput, logFilePath);
       }
    }

    // Create BlastX checkpoint
    trans.makeCheckpoint("blastx");
  }
}

// Given a vector of transcript objects, detect and remove chimera transcripts from each's data
void remChimeraBulk(const std::vector<transcript> & transVec, std::string ram_gb,
                    bool retainInterFiles, bool dispOutput, std::string logFilePath,
                    const INI_MAP & cfgIni = {}) {
  logOutput("\nStarting chimera removal\n", logFilePath);
  std::string currTransInChim;
  std::string currTransOutChim;
  std::string currBlastx;
  std::string currTransInfo;
  std::string currTransCut;
  std::string chimOutDir;
  for (auto trans : transVec) {
    // Check if chimera removal checkpoint exists 
    if (trans.checkpointExists("chim.fix")) {
      logOutput("\n  Chimera removal checkpoint found for: " +
                trans.get_file_prefix(), logFilePath);
      continue;
    }
    logOutput("\n  Now running chimera removal for:\n", logFilePath);
    summarize_sing_trans(trans, logFilePath, 4);

   
    currTransInChim = trans.get_trans_path_trinity().c_str();
    currTransOutChim = trans.get_trans_path_chimera().c_str();
    currBlastx = trans.get_trans_path_blastx().c_str();
    currTransInfo = trans.get_trans_path_cinfo().c_str();
    currTransCut = trans.get_trans_path_ccut().c_str();
    chimOutDir = trans.get_trans_path_chimera().parent_path().c_str();

    if (!dispOutput) {
      procRunning = true;
      std::thread chimeraThread(progressAnim, "  ", logFilePath);
      detect_chimera(currBlastx, chimOutDir);
      removeChimera(currTransInChim, currTransOutChim, currTransInfo, currTransCut, ram_gb,
                    logFilePath);
      procRunning = false;
      chimeraThread.join();
    }
    else {
      detect_chimera(currBlastx, chimOutDir);
      removeChimera(currTransInChim, currTransOutChim, currTransInfo, currTransCut, ram_gb,
                    logFilePath);
    }
  
    // Create chimera removal checkpoint
    trans.makeCheckpoint("chim.fix");
  }
}

// Given a vector of transcript objects, create a salmon index for each and then
// quantify each's corresponding reads against its index
void salmonBulk(const std::vector<transcript> & transVec, std::string threads,
                bool retainInterFiles, bool dispOutput, std::string logFilePath,
                std::vector<std::string> readFiles1, std::vector<std::string> readFiles2,
                std::string outDir, const INI_MAP & cfgIni = {}) {
  logOutput("\nStarting transcript quantification\n", logFilePath);
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
  INI_MAP_ENTRY cfgIniPipeline;
  if (!cfgIni.empty()) {
    cfgIniPipeline = cfgIni.at("Pipeline");
  }

  /*
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
  */
  for (auto trans : transVec) {
    currTransInSalm = trans.get_trans_path_chimera().c_str();
    if (!cfgIni.empty()) {
      if (!ini_get_bool(cfgIniPipeline.at("remove_chimera_reads").c_str(), 0)) {
        currTransInSalm = trans.get_trans_path_trinity().c_str();
      }
      currTransInfoFileStr = std::string(trans.get_trans_path_trinity().replace_extension(".ti").c_str());
      currTransInfoFile.open(currTransInfoFileStr);
      while (getline(currTransInfoFile, currLineInfo)) {
        spacePos = currLineInfo.find(" ");
        nlPos = currLineInfo.find("\n");
        if (spacePos != std::string::npos) {
          currSraRun.first = currLineInfo.substr(0, spacePos + 1);
          currSraRun.second = currLineInfo.substr(spacePos + 1, nlPos - spacePos);
        }
        else {
          currSraRun.first = currLineInfo;
        }
        currSraRunsIn.push_back(currSraRun);
      }
      currTransInfoFile.close();
    }
    else {
      for (int i = 0; i < readFiles1.size(); i++) {
        currSraRun.first = readFiles1[i];
        if (!readFiles2.empty()) {
          currSraRun.second = readFiles2[i];
        }
        currSraRunsIn.push_back(currSraRun);
      }
    } 
    currIndex = trans.get_trans_path_index().c_str();
    currQuant = trans.get_trans_path_quant().c_str();
   
    // Check if salmon index checkpoint exists
    if (trans.checkpointExists("index")) {
      logOutput("  Indexing checkpoint found for: " +
                trans.get_file_prefix() + "\n", logFilePath);
    }
    else {
      // Summarize indexing job
      logOutput("  Now producing index for: ", logFilePath);
      summarize_sing_trans(trans, logFilePath, 4);

      // Create index of assembled transcripts
      if (!dispOutput) {
        procRunning = true;
        std::thread salmIdxThread(progressAnim, "  ", logFilePath);
        salmon_index(currTransInSalm, currIndex, threads, dispOutput, logFilePath);
        procRunning = false;
        salmIdxThread.join();
      }
      else {
        salmon_index(currTransInSalm, currIndex, threads, dispOutput, logFilePath);
      }

      // Create index checkpoint
      trans.makeCheckpoint("index");
    }

    // Check if salmon quant checkpoint exists
    if (trans.checkpointExists("quant")) {
      logOutput("  Quant checkpoint found for: " +
                trans.get_file_prefix() + "\n", logFilePath);
    }
    else {
      // Summarize quantification job
      logOutput("  Now quantifying transcripts: \n", logFilePath);
      summarize_sing_trans(trans, logFilePath, 4);

      // Quantify reads mapped to transcript index
      if (!dispOutput) {
        procRunning = true;
        std::thread salmQntThread(progressAnim, "  ", logFilePath);
        salmon_quant(currTransInSalm, currIndex, currQuant, currSraRunsIn, threads,
                     dispOutput, logFilePath);
        procRunning = false;
        salmQntThread.join();
      }
      else {
        salmon_quant(currTransInSalm, currIndex, currQuant, currSraRunsIn, threads,
                     dispOutput, logFilePath);
      }

      // Create salmon quant checkpoint
      trans.makeCheckpoint("quant");
    }
    currSraRunsIn.clear();
  }
}

// Given a vector of transcript objects, determine each's largest cluster / redundant
// transcripts using Corset, and then generate output files for both
void corsetBulk(const std::vector<transcript> & transVec, std::string ram_gb,
                bool retainInterFiles, bool dispOutput, std::string logFilePath,
                const INI_MAP & cfgIni = {}) {
  logOutput("\nRemoving redundant transcripts\n", logFilePath);
  std::string currTransInCors;
  std::string currTransPrefix;
  std::string currEqClassFile;
  std::string currTransClust;
  std::string currTransLargestClust;
  std::string currTransRedund;
  std::string currOutDir;
  INI_MAP_ENTRY cfgIniPipeline;
  uintmax_t ram_b = (uintmax_t)stoi(ram_gb) * 1000000000;
  if (!cfgIni.empty()) {
    cfgIniPipeline = cfgIni.at("Pipeline");
  }
  /*
  INI_MAP_ENTRY cfgPipeline = cfgIni.at("Pipeline");
  std::string currTransInCors;
  std::string currTransPrefix;
  std::string currEqClassFile;
  std::string currTransClust;
  std::string currTransLargestClust;
  std::string currTransRedund;
  std::string currOutDir;
  uintmax_t ram_b = (uintmax_t)stoi(ram_gb) * 1000000000;
  */
  for (auto trans : transVec) {
    // Check if corset run checkpoint exists
    if (trans.checkpointExists("clust")) {
      logOutput("  Corset checkpoint found for: " +
                trans.get_file_prefix() + "\n", logFilePath);
      continue;
    }
    logOutput("  Now clustering transcripts into genes:\n", logFilePath);
    summarize_sing_trans(trans, logFilePath, 4);
    currTransInCors = trans.get_trans_path_chimera().c_str();
    if (!cfgIni.empty()) {
      if (!ini_get_bool(cfgIniPipeline.at("remove_chimera_reads").c_str(), 0)) {
        currTransInCors = trans.get_trans_path_trinity().c_str();
      }
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
    if (!dispOutput) {
      procRunning = true;
      std::thread clustThread(progressAnim, "  ", logFilePath);
      corset_eq_classes(trans.get_file_prefix(), currEqClassFile, currOutDir,
                        dispOutput, logFilePath);
      procRunning = false;
      clustThread.join();
    }
    else {
       corset_eq_classes(trans.get_file_prefix(), currEqClassFile, currOutDir,
                        dispOutput, logFilePath);
    }

    // Create corset run checkpoint
    trans.makeCheckpoint("clust");

    // Check if corset filter checkpoint exists
    if (trans.checkpointExists("redund.fix")) {
      logOutput("\n  Corset filtering checkpoint found for: " +
                trans.get_file_prefix(), logFilePath);
      continue;
    }
    // Filter corset output
    logOutput("\n  Now filtering redundant transcripts:\n", logFilePath);
    summarize_sing_trans(trans, logFilePath, 4);

    if (!dispOutput) {
      procRunning = true;
      std::thread clustFiltThread(progressAnim, "  ", logFilePath);
      filterCorset(currTransInCors, currTransClust, currTransLargestClust, currTransRedund,
                   ram_b, currOutDir, logFilePath);
      procRunning = false;
      clustFiltThread.join();
    }
    else {
      filterCorset(currTransInCors, currTransClust, currTransLargestClust, currTransRedund,
                   ram_b, currOutDir, logFilePath);
    }

    // Create corset filtering checkpoint
    trans.makeCheckpoint("redund.fix");
  }
}

// Given a vector of transcript objects, predict each's coding regions using TransDecoder
void transdecBulk(const std::vector<transcript> & transVec,
                  std::string threads, std::string ram_gb,
                  bool retainInterFiles, bool dispOutput,
                  std::string logFilePath,
                  std::string refProt, std::string outDir,
                  const INI_MAP & cfgIni = {}) {
  logOutput("\nStarting prediction of coding regions\n", logFilePath);
  std::string currTransInTD;
  std::string currTransCds;
  std::string currTransPep;
  std::string currDb;
  std::string currOutDirTD;
  INI_MAP_ENTRY cfgIniGen;
  INI_MAP_ENTRY cfgIniAlign;
  INI_MAP_ENTRY cfgIniPipeline;
  bool useBlast;
  std::string maxEvalue;
  std::string maxTargetSeqs;
  std::string refProteome;
  std::string blastDbDir;
  std::string blastDbName;
  uintmax_t ram_b = (uintmax_t)stoi(ram_gb) * 1000000000;
  if (!cfgIni.empty()) {
    cfgIniGen = cfgIni.at("General");
    cfgIniAlign = cfgIni.at("Alignment settings");
    cfgIniPipeline = cfgIni.at("Pipeline");
    useBlast = ini_get_bool(cfgIniAlign.at("use_blast_instead_of_diamond").c_str(), 0);
    maxEvalue = cfgIniAlign["blastp_max_evalue"];
    maxTargetSeqs = cfgIniAlign["blastp_max_target_seqs"];
    refProteome = cfgIniGen["reference_proteome_path"];
    if (!fs::exists(fs::path(refProteome.c_str()))) {
      logOutput("ERROR: Reference proteome " + refProteome + " not found\n", logFilePath);
      logOutput("  Please ensure its path is correctly specified in your config file\n", logFilePath);
      exit(1);
    }
    blastDbDir = cfgIniGen["output_directory"] + "/" +
                 cfgIniGen["project_name"] + "/" +
                 stepDirs[8] + "/";
    blastDbName = std::string(fs::path(refProteome.c_str()).stem().c_str());
  }
  else {
    useBlast = false;
    maxEvalue = "10";
    maxTargetSeqs = "1";
    refProteome = refProt;
    if (!fs::exists(fs::path(refProteome.c_str()))) {
      logOutput("ERROR: Reference proteome " + refProteome + " not found\n", logFilePath);
      exit(1);
    }
    blastDbDir = outDir + stepDirs[0] + "/";
    blastDbName = std::string(fs::path(refProteome.c_str()).stem().c_str());
  }
  /*
  INI_MAP_ENTRY cfgPipeline = cfgIni.at("Pipeline");
  std::string currTransInTD;
  std::string currTransCds;
  std::string currTransPep;
  std::string currDb;
  std::string currOutDirTD;
  INI_MAP_ENTRY cfgGen = cfgIni.at("General");
  INI_MAP_ENTRY cfgAlign = cfgIni.at("Alignment settings");
  bool useBlast = ini_get_bool(cfgAlign.at("use_blast_instead_of_diamond").c_str(), 0);
  std::string maxEvalue = cfgAlign.at("blastp_max_evalue");
  std::string maxTargetSeqs = cfgAlign.at("blastp_max_target_seqs");
  std::string refProt = cfgGen.at("reference_proteome_path");
  if (!fs::exists(fs::path(refProt.c_str()))) {
    logOutput("ERROR: Reference proteome: " + refProt + " not found\n", logFilePath);
    logOutput("  Please ensure its path is correctly specified in your config INI file\n",
              logFilePath);
    exit(1);
  }
  std::string blastDbDir = cfgIni.at("General").at("output_directory") + "/" +
                           cfgIni.at("General").at("project_name") + "/" +
                           stepDirs[8] + "/";
  std::string blastDbName(fs::path(refProt.c_str()).stem().c_str());
  uintmax_t ram_b = (uintmax_t)stoi(ram_gb) * 1000000000;
  */
  for (auto trans : transVec) {
    // Check if coding region prediction checkpoint exists
    if (trans.checkpointExists("cdr.predict")) {
      logOutput("\n  Coding region prediction checkpoint found for: " +
                trans.get_file_prefix(), logFilePath);
      continue;
    }
    logOutput("\n  Now predicting coding regions for:\n", logFilePath);
    summarize_sing_trans(trans, logFilePath, 4);
    currTransInTD = trans.get_trans_path_largest().c_str();
    if (!cfgIni.empty()) {
      if (!ini_get_bool(cfgIniPipeline.at("cluster_filtering").c_str(), 0)) {
        if (!ini_get_bool(cfgIniPipeline.at("remove_chimera_reads").c_str(), 0)) {
          currTransInTD = trans.get_trans_path_trinity().c_str();
        }
        else {
          currTransInTD = trans.get_trans_path_chimera().c_str();
        }
      }
    }
    else {
      currTransInTD = trans.get_trans_path_chimera().c_str();
    }
    currTransCds = trans.get_trans_path_cds().c_str();
    currTransPep = trans.get_trans_path_prot().c_str();
    
    currDb = blastDbDir + "/" + blastDbName;
    currOutDirTD = trans.get_trans_path_cds().parent_path().c_str();
       
    // Perform transdecoder run to obtain coding sequences / peptides
    if (!dispOutput) {
      procRunning = true;
      std::thread transDecThread(progressAnim, "  ", logFilePath);
      run_transdecoder(currTransInTD, currTransCds, currTransPep, useBlast, maxEvalue,
                       maxTargetSeqs, threads, ram_b, currDb, currOutDirTD, dispOutput,
                       logFilePath);
      procRunning = false;
      transDecThread.join();
    }
    else {
      run_transdecoder(currTransInTD, currTransCds, currTransPep, useBlast, maxEvalue,
                       maxTargetSeqs, threads, ram_b, currDb, currOutDirTD, dispOutput,
                       logFilePath);
    }
    // Create coding region prediction checkpoint
    trans.makeCheckpoint("cdr.predict");
  }
}

// Given a vector of transcript objects, run a HMMER/PANTHER-based annotation for each
void annotateBulk(const std::vector<transcript> & transVec, std::string threads,
                  std::string ram_gb, bool retainInterFiles, bool dispOutput,
                  std::string logFilePath, const INI_MAP & cfgIni = {}) {
  logOutput("\nStarting annotation of transcripts\n", logFilePath);
  std::string currTransIn;
  std::string currTransPep;
  std::string currTransOut;
  for (auto trans : transVec) {
    // Check if annotation checkpoint exists 
    if (trans.checkpointExists("annotate")) {
      logOutput("\n  Annotation checkpoint found for: " +
                trans.get_file_prefix(), logFilePath);
      continue;
    }
    logOutput("\n  Now running annotation on:\n", logFilePath);
    summarize_sing_trans(trans, logFilePath, 4);

    currTransIn = trans.get_trans_path_cds().c_str();
    currTransPep = trans.get_trans_path_prot().c_str();
  
    currTransOut = trans.get_trans_path_annot().c_str();

    // Perform annotation of transcript
    if (!dispOutput) {
      procRunning = true;
      std::thread annotThread(progressAnim, "  ", logFilePath);
      annotateTranscript(currTransIn, currTransPep, currTransOut,
                         threads, ram_gb, dispOutput, logFilePath);
      procRunning = false;
      annotThread.join();
    }
    else {
      annotateTranscript(currTransIn, currTransPep, currTransOut,
                         threads, ram_gb, dispOutput, logFilePath);
    }

    // Create annotation checkpoint
    trans.makeCheckpoint("annotate");
  }
}

std::vector<std::string> getCommaSepStrings(std::string commaSepStrings) {
  std::vector<std::string> stringVec;
  size_t commaInd = commaSepStrings.find(',');
  size_t currPos = 0;
  std::string currStr;
  do {
    currStr = commaSepStrings.substr(currPos, commaInd + 1);
    std::cout << currStr << std::endl;
    stringVec.push_back(currStr);
  } while (commaInd != std::string::npos);
  return stringVec;
}

int main(int argc, char * argv[]) {
  system("setterm -cursor off");
  std::string configPath;
  INI_MAP cfgIni;
  INI_MAP_ENTRY cfgIniGen;
  INI_MAP_ENTRY cfgIniPipeline;
  std::string leftReads;
  std::string rightReads;
  std::string assembly;
  std::string refProt;
  std::string outDir;
  std::string projDir;
  std::string logFilePath;
  std::vector<std::string> readFilesLeft;
  std::vector<std::string> readFilesRight;


  const char * home = std::getenv("HOME");
  std::string threads;
  std::string ram_gb;
  bool retainInterFiles;
  bool dispOutput;
  bool compressFiles = false;

  std::vector<SRA> sras;
  std::vector<std::string> localDataFiles;
  std::pair<std::string, std::string> sraRunsLocal;
  size_t posLocalFiles;

  std::vector<transcript> transVec;
  if (argc > 1) {
    threads = argv[7];
    ram_gb = argv[8];
    retainInterFiles = stringToBool(argv[9]);
    dispOutput = stringToBool(argv[10]);
    // Determine whether sequence files are specified in config file or postprocess call
    configPath = argv[1];
    // Set up postprocess parameters based on command line input
    if (configPath == "null") {
      leftReads = argv[2];
      rightReads = argv[3];
      assembly = argv[4];
      refProt = argv[5];
      outDir = argv[6];
      readFilesLeft = getCommaSepStrings(leftReads);
      readFilesRight = getCommaSepStrings(rightReads);
      logFilePath = "log.txt";
      if (readFilesLeft.size() != readFilesRight.size()) {
        logOutput("ERROR: Number of left/right read files do not match", logFilePath);
        exit(1);
      }
      make_proj_space(outDir, "postprocess");
      transcript transFile(assembly, outDir);
      transVec.push_back(transFile);
      // DO POSTPROCESS STEPS
      
    }
    // Set up postprocess parameters based on user-specified config file
    else {
      // Get INI config file
      cfgIni = make_ini_map(argv[1]);
      cfgIniGen = cfgIni["General"];
      cfgIniPipeline = cfgIni["Pipeline"];
      make_proj_space(cfgIni, "postprocess");
      outDir = cfgIniGen["output_directory"];
      if (outDir[0] == '~') {
        outDir = std::string(home) + outDir.substr(1, outDir.size() - 1);
      }
      projDir = outDir + cfgIniGen["project_name"] + "/";
      refProt = cfgIniGen["reference_proteome_path"];
      if (refProt[0] == '~') {
        refProt = std::string(home) + refProt.substr(1, refProt.size() - 1);
      }
      logFilePath = std::string((fs::canonical(fs::path(cfgIniGen["log_file"].c_str()).parent_path()) /
                                fs::path(cfgIniGen["log_file"].c_str()).filename()).c_str());
      if (logFilePath[0] == '~') {
        logFilePath = std::string(home) + logFilePath.substr(1, logFilePath.size() - 1);
      }
      sras = get_sras(cfgIni, dispOutput, compressFiles, logFilePath);
      for (auto fqFileName : cfgIni.at("Local files")) {
        localDataFiles.push_back(fqFileName.first);
      }
      for (auto sraRun : localDataFiles) {
        sraRunsLocal.first = "";
        sraRunsLocal.second = "";
        posLocalFiles = sraRun.find(' ');
        if (posLocalFiles != std::string::npos) {
          sraRun.erase(0, posLocalFiles + 1);
          posLocalFiles = sraRun.find(' ');
          sraRunsLocal.second = sraRun.substr(0, posLocalFiles);
        }
        if (fs::exists(sraRunsLocal.first) &&
            fs::exists(sraRunsLocal.second)) {
          sras.push_back(SRA(sraRunsLocal.first, sraRunsLocal.second, cfgIni, compressFiles,
                             logFilePath));
        }
        else {
          if (sraRunsLocal.first != "" &&
              !fs::exists(sraRunsLocal.first)) {
            logOutput("ERROR: Local run not found: \"" + sraRunsLocal.first + "\"\n",
                      logFilePath);
          }
          if (sraRunsLocal.second != "" &&
              !fs::exists(sraRunsLocal.second)) {
            logOutput("ERROR: Local run not found: \"" + sraRunsLocal.second + "\"\n",
                      logFilePath);
          }
        }
      }
      if (sras.empty()) {
        logOutput("ERROR: No valid SRA data. Please check the config file.\n",
                  logFilePath);
        exit(1);
      }
      fs::path currFile = transcript(sras[0]).get_trans_path_trinity().parent_path();
      fs::path currPrefix;
      fs::directory_iterator fileIter{currFile};
      while (fileIter != fs::directory_iterator{}) {
        // Iterate through files in Semblans transcript directory
        // Push transcript to transVec transcripts vector
        if (fileIter->path().extension() == ".fasta") {
          currFile = fileIter->path();
          currPrefix = currFile;
          while (!currPrefix.extension().empty()) {
            currPrefix = currPrefix.stem();
          }
          fs::rename(currFile, currFile.parent_path() /
                               fs::path((std::string(currPrefix.c_str()) + ".Trinity.fasta").c_str()));
          transcript currTrans(fileIter->path().c_str(), cfgIni);
          transVec.push_back(currTrans);
        }
        fileIter++;
      }
    }

    // Determine whether to compress files
    //bool compressFiles = ini_get_bool(cfgIni["General"]["compress_files"].c_str(), 0);
    
    // Retrieve SRA objects defined by , convert to transcripts
    /*
    std::vector<SRA> sras;
    std::vector<std::string> localDataFiles;
    sras = get_sras(cfgIni, dispOutput, compressFiles, logFilePath);
    for (auto fqFileName : cfgIni.at("Local files")) {
      localDataFiles.push_back(fqFileName.first);
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
      if (fs::exists(sraRunsLocal.first) &&
          fs::exists(sraRunsLocal.second)) {
        sras.push_back(SRA(sraRunsLocal.first, sraRunsLocal.second, cfgIni, compressFiles,
                           logFilePath));
      }
      else {
        if (sraRunsLocal.first != "" &&
            !fs::exists(sraRunsLocal.first)) {
          logOutput("ERROR: Local run not found: \"" + sraRunsLocal.first + "\"\n",
                    logFilePath);
        }
        if (sraRunsLocal.second != "" &&
            !fs::exists(sraRunsLocal.second)) {
          logOutput("ERROR: Local run not found: \"" + sraRunsLocal.second + "\"\n",
                    logFilePath);
        }
      }
    }
    
    if (sras.empty()) {
      logOutput("ERROR: No valid SRA data. Please check the configuration file.\n",
                logFilePath);
      exit(1);
    }
    
    std::vector<transcript> transVec;
    fs::path currFile = transcript(sras[0]).get_trans_path_trinity().parent_path();
    fs::path currPrefix;
    fs::directory_iterator fileIter{currFile};
    while (fileIter != fs::directory_iterator{}) {
      // Iterate through files in transcript directory
      // Push transcript to transcript vector
      if (fileIter->path().extension() == ".fasta") {
        currFile = fileIter->path();
        currPrefix = currFile;
        while (!currPrefix.extension().empty()) {
          currPrefix = currPrefix.stem();
        }
        fs::rename(currFile, currFile.parent_path() /
                             fs::path((std::string(currPrefix.c_str()) + ".Trinity.fasta").c_str()));
        transcript currTrans(fileIter->path().c_str(), cfgIni);
        transVec.push_back(currTrans);
      }
      fileIter++;
    }
 
    for (auto trans : transVec) {
      std::cout << trans.get_trans_path_trinity().c_str() << std::endl;
    }*/
  
    // Summarize program execution parameters
    logOutput("Semblans Postprocess started with following parameters:\n", logFilePath);
    logOutput("  Config file:     " + std::string(argv[1]) + "\n", logFilePath);
    logOutput("  Threads (Cores): " + threads + "\n", logFilePath);
    logOutput("  Memory (GB):     " + ram_gb + "\n", logFilePath);
    logOutput("  Reference Prot:  " + refProt + "\n", logFilePath);
    logOutput("  SRA runs:\n", logFilePath);
    summarize_all_sras(sras, logFilePath, 6);

    if (configPath != "null") {
      if (ini_get_bool(cfgIniPipeline.at("remove_chimera_reads").c_str(), 0)) {
        // BlastX alignment of transcript to reference proteome
        blastxBulk(transVec, threads, dispOutput, logFilePath, "", "", cfgIni);

        // Detect and remove chimeric transcripts
        remChimeraBulk(transVec, ram_gb, retainInterFiles, dispOutput, logFilePath, cfgIni);
      }
    }
    else {
      blastxBulk(transVec, threads, dispOutput, logFilePath, refProt, outDir);
      remChimeraBulk(transVec, ram_gb, retainInterFiles, dispOutput, logFilePath);
    }

    if (configPath != "null") {
      if (ini_get_bool(cfgIniPipeline.at("cluster_filtering").c_str(), 0)) {
        // Perform salmon index of transcripts
        salmonBulk(transVec, threads, retainInterFiles, dispOutput, logFilePath,
                   {}, {}, "", cfgIni);
   
        // Perform corset run to cluster transcripts
        corsetBulk(transVec, ram_gb, retainInterFiles, dispOutput, logFilePath, cfgIni);
      }
    }
    else {
      salmonBulk(transVec, threads, retainInterFiles, dispOutput, logFilePath,
                 readFilesLeft, readFilesRight, outDir);
      corsetBulk(transVec, ram_gb, retainInterFiles, dispOutput, logFilePath, cfgIni);
    }
    
    // Run transdecoder
    if (configPath != "null") {
      transdecBulk(transVec, threads, ram_gb, retainInterFiles, dispOutput, logFilePath, "", "", cfgIni);
    }
    else {
      transdecBulk(transVec, threads, ram_gb, retainInterFiles, dispOutput, logFilePath, refProt, outDir);
    }

    if (configPath != "null") {
      if (ini_get_bool(cfgIniPipeline.at("annotate_transcripts").c_str(), 0)) {
        // Annotate transcriptome
        annotateBulk(transVec, threads, ram_gb, retainInterFiles, dispOutput, logFilePath, cfgIni);
      }
    }
    logOutput("\nPostprocess finished successfully\n\n", logFilePath);
  }
  else {
  
  }
  system("setterm -cursor on");
}
