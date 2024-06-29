#include "assemble.h"

std::atomic<bool> procRunning(false);

// Update Trinity output transcript headers according to a given prefix string
void updateHeaders(std::string fastaFilePath, std::string newPrefix,
                   uintmax_t ram_b) {

  uintmax_t numBytesFasta = fs::file_size(fastaFilePath.c_str());
  uintmax_t lenHashTable = numBytesFasta / 160;

  seqHash fastaHashTable(lenHashTable, fastaFilePath, ram_b);

  linkedList * hashData = fastaHashTable.getHashData();

  std::string currOldHeader;
  std::string currNewHeader;
  sequence currSeq;
  for (uintmax_t i = 0; i < lenHashTable; i++) {
    if (hashData[i].empty()) {
      continue;
    }
    else {
      hashData[i].updateSeqHead(newPrefix);
    }
  }
  fastaHashTable.dump(fastaFilePath);
}

// Given a FASTA input, divide all SRA reads into those which map to it and those that do not.
// Create two new files for containing the two.
void isolateReads(const std::vector<SRA> & sras, std::string threads,
                  std::string ram_gb, bool dispOutput, const INI_MAP & cfgIni, std::string logFile) {

  logOutput("\nStarting isolation of reads of interest\n", logFile);
  INI_MAP_ENTRY cfgPipeline = cfgIni.at("Pipeline");
  INI_MAP_ENTRY cfgSeqsInt = cfgIni.at("Sequences of interest");
  std::string maxMapsPerRead = cfgIni.at("STAR settings").at("max_maps_per_read");
  std::string fastaIndex;
  std::string fastaMap;
  transcript dummyTrans;
  std::ifstream bamMapFile;
  std::stringstream percentStream;

  // Create string vector of user-defined FASTA files containing sequences to extract
  std::vector<std::string> seqsInterest;
  std::vector<std::string> seqsDecoy;
  for (auto seqFilePath : cfgSeqsInt) {
    seqsInterest.push_back(seqFilePath.first);
  }
  
  uintmax_t ram_b = (uintmax_t)stoi(ram_gb) * 1000000000;
  std::pair<std::string, std::string> currSraIn;
  std::pair<std::string, std::string> currSraUmap;
  std::vector<std::pair<std::string, std::string>> sraRunsIn;
  std::string outDir;
  std::string filePrefix1;
  std::string filePrefix2;

  uintmax_t numBytesReads1;
  uintmax_t lenReadsHash1;
  uintmax_t numBytesReads2;
  uintmax_t lenReadsHash2;
  
  sequence currSeq;
  std::string currHead;
  std::string currSeqData;
  std::string currQual;

  std::string currSeqFilePrefix;
  std::string currLine;
  uintmax_t numReadsSeqInterest = 0;
  std::ifstream currSeqFile;
  // Iterate through sequences of interest, generating a STAR index for each and mapping
  // reads against them
  for (int i = 0; i < seqsInterest.size(); i++) {
    currSeqFile.open(seqsInterest[i]);
    while (std::getline(currSeqFile, currLine)) {
      if (currLine[0] == '>' || currLine[0] == '@') {
        numReadsSeqInterest++;
      }
    }
    logOutput("\n  Now running read extraction using: \"" + seqsInterest[i] + "\"\n", logFile);
    currSeqFilePrefix = std::string(fs::path(seqsInterest[i].c_str()).stem().c_str());
    // Create index of sequence file with STAR
    transcript dummyTrans(seqsInterest[i], cfgIni);
    fastaIndex = std::string(dummyTrans.get_trans_path_trinity().parent_path().c_str()) + "/" +
                 std::string(fs::path(seqsInterest[i].c_str()).stem().c_str()) + "_index";
    // Check if indexing checkpoint exists. If not, generate index
    if (!dummyTrans.checkpointExists(currSeqFilePrefix + ".idx.iso")) {
      fs::create_directory(fastaIndex.c_str());
      if (!dispOutput) {
        procRunning = true;
        std::thread seqIndex(progressAnim, "    Now creating mapping index for \"" + currSeqFilePrefix + "\"", logFile);
        
        star_index(std::vector<std::string>(1, seqsInterest[i]), fastaIndex, threads,
                   dispOutput, logFile);

        procRunning = false;
        seqIndex.join();
      }
      else {
        logOutput("    Now creating mapping index for \"" + currSeqFilePrefix + "\"\n", logFile);

        star_index(std::vector<std::string>(1, seqsInterest[i]), fastaIndex, threads,
                   dispOutput, logFile);
        
      }
      dummyTrans.makeCheckpoint(currSeqFilePrefix + ".idx.iso");
    }
  
    uintmax_t numReads; 
    uintmax_t numMapped;
    int num1k;
    float percentMapped;
    // Iterate through SRA runs, mapping each against the STAR index we just created
    for (auto sra : sras) {
      if (sra.checkpointExists(currSeqFilePrefix + ".iso")) {
        if (!sra.get_accession().empty()) {
          logOutput("\n    Read isolation checkpoint found for: " +
                    sra.get_accession(), logFile);
        }
        else {
          logOutput(std::string("\n    Read isolation checkpoint found for: ") +
                    "\n      " + sra.get_file_prefix().first +
                    "\n      " + sra.get_file_prefix().second,
                    logFile);
        }
        continue;
      }
      logOutput("\n    Extracting reads from:\n", logFile);
      summarize_sing_sra(sra, logFile, 4);
      outDir = sra.get_sra_path_orep_filt().first.parent_path().c_str();

      currSraIn.first = sra.get_sra_path_orep_filt().first.c_str();
      currSraIn.second = sra.get_sra_path_orep_filt().second.c_str();

      // Obtain read input files, based on pipeline preferences
      if (!ini_get_bool(cfgPipeline.at("remove_overrepresented").c_str(), 0)) {
        if (!ini_get_bool(cfgPipeline.at("filter_foreign_reads").c_str(), 0)) {
          if (!ini_get_bool(cfgPipeline.at("trim_adapter_seqs").c_str(), 0)) {
            if (!ini_get_bool(cfgPipeline.at("error_correction").c_str(), 0)) {
              currSraIn.first = sra.get_sra_path_raw().first.c_str();
              currSraIn.second = sra.get_sra_path_raw().second.c_str();
              outDir = sra.get_sra_path_raw().first.parent_path().c_str();
            }
            else {
              currSraIn.first = sra.get_sra_path_corr_fix().first.c_str();
              currSraIn.second = sra.get_sra_path_corr_fix().second.c_str();
              outDir = sra.get_sra_path_corr_fix().first.parent_path().c_str();
            }
          }
          else {
            currSraIn.first = sra.get_sra_path_trim_p().first.c_str();
            currSraIn.second = sra.get_sra_path_trim_p().second.c_str();
            outDir = sra.get_sra_path_trim_p().first.parent_path().c_str();
          }
        }
        else {
          currSraIn.first = sra.get_sra_path_for_filt().first.c_str();
          currSraIn.second = sra.get_sra_path_for_filt().second.c_str();
          outDir = sra.get_sra_path_for_filt().first.parent_path().c_str();
        }
      }
      filePrefix1 = sra.get_file_prefix().first;
      filePrefix2 = sra.get_file_prefix().second;
      
      if (i > 0) {
        currSraIn.first = outDir + "/" + sra.get_file_prefix().first + ".unmapped.fq";
        currSraIn.second = outDir + "/" + sra.get_file_prefix().second + ".unmapped.fq";
      }

      sraRunsIn.push_back(currSraIn);

      // Define name of STAR mapping output file      
      if (!sra.is_paired()) {
        fastaMap = std::string(dummyTrans.get_trans_path_trinity().parent_path().c_str()) + "/" +
                     std::string(fs::path(currSeqFilePrefix.c_str()).stem().c_str()) + "_" + 
                     sra.get_file_prefix().first + "_mapping/";
      }
      else {
        fastaMap = std::string(dummyTrans.get_trans_path_trinity().parent_path().c_str()) + "/" +
                     std::string(fs::path(currSeqFilePrefix.c_str()).stem().c_str()) + "_" +
                     sra.get_file_prefix().first.substr(0, sra.get_file_prefix().first.size() - 1) +
                     "mapping/";
      }

      // Determine whether checkpoint exists for STAR map of current reads against current index.
      // If not, perform a quant run
      if (!sra.checkpointExists(currSeqFilePrefix + ".map.iso")) {
        fs::create_directory(fastaMap.c_str());
        if (!dispOutput) {
          procRunning = true;
          std::thread seqQuant(progressAnim, "      Mapping reads to sequences of interest ",
                               logFile);
          star_map(fastaIndex, fastaMap, currSraIn, maxMapsPerRead, threads, dispOutput, logFile);
          procRunning = false;
          seqQuant.join();
        }
        else {
          logOutput("\n      Mapping reads to sequences of interest\n", logFile);
          star_map(fastaIndex, fastaMap, currSraIn, maxMapsPerRead, threads, dispOutput, logFile);
        }

        // Create checkpoint for STAR mapping of SRA against seqs of interest
        sra.makeCheckpoint(currSeqFilePrefix + ".map.iso");
      }

      // The STAR mapping run has produced a text file 'unmapped_names.txt', which lists all the reads in
      // the SRA run that did NOT map to the index.
      
      // Create sequence hash tables for rapid accession of each read by name
      // Then iterate through all unmapped reads listed in 'unmapped_names.txt', and maintain two sequence hash tables:
      
      //   Unmapped hash table - will contain all reads that did not map to the STAR index (updated after each index)
      //     Mapped hash table - will contain all reads that mapped to the STAR index
      numBytesReads1 = fs::file_size(currSraIn.first.c_str());
      lenReadsHash1 = numBytesReads1 / 160;
      
      procRunning = true;
      std::thread constructHash(progressAnim, "      Creating hash table from forward-ended reads ", logFile);
      seqHash readHashTable(lenReadsHash1, fs::path(currSraIn.first.c_str()), ram_b);
      procRunning = false;
      constructHash.join();
      numReads = readHashTable.getNumItems();
      
      seqHash mappedHash(lenReadsHash1);

      logOutput("      Splitting reads into mapped and unmapped\n", logFile);
      bamMapFile.open(fastaMap + "/Aligned.out.sam");

      // Iterate through headers in unmapped reads file
      // Fill hash tables accordingly
      for (int i = 0; i < numReadsSeqInterest + 3; i++) {
        std::getline(bamMapFile, currLine);
      }
      std::string readName;
      numMapped = 0;
      num1k = 0; 
      while (std::getline(bamMapFile, currLine)) {
        readName = currLine.substr(0, currLine.find("\t"));
        if (!mappedHash.inHashTable(readName)) {
          numMapped++;
          currSeq = readHashTable.getSeq(readName);
          currHead = currSeq.get_header();
          currSeqData = currSeq.get_sequence();
          currQual = currSeq.get_quality();

          readHashTable.deleteHash(readName);
          mappedHash.insertHash(currHead, currSeqData, currQual);
          if (numMapped % 1000 == 0) {
            num1k++;
            std::cout << "\r        Mapped read count: " +
                         std::to_string(1000 * num1k) + " ..." << std::flush;
          }
        }
      }
      percentMapped = (float(numMapped) * 100) / float(numReads);
      logOutput("\r        Mapped read count:   " + std::to_string(numMapped) +
                " (" + getPercent(percentMapped, 2) + "%)", logFile);
      logOutput("\n        Unmapped read count: " + std::to_string(numReads - numMapped) +
                " (" + getPercent(100.0 - percentMapped, 2) + "%)\n", logFile);
      // Dump filled sequence hash tables to new files, containing the mapped
      // and unmapped reads respectively

      procRunning = true;
      std::thread dumpHash(progressAnim, "      Dumping split reads to mapped and unmapped files ", logFile);
      mappedHash.dump(outDir + "/" + filePrefix1 + "." + currSeqFilePrefix + ".mapped.fq");
      readHashTable.dump(outDir + "/" + filePrefix1 + ".unmapped.fq");
      procRunning = false;
      dumpHash.join();
      
      readHashTable.clear();
      mappedHash.clear();
      bamMapFile.close();
      logOutput("      Reads were extracted successfully\n", logFile);

      // If SRA run is paired, perform the same steps as above for the reverse-ended FASTQ file
      if (sra.is_paired()) {
        logOutput("\n", logFile);
        numBytesReads2 = fs::file_size(currSraIn.second.c_str());
        lenReadsHash2 = numBytesReads2 / 160;
        
        procRunning = true;
        std::thread constructHash(progressAnim, "      Creating hash table from reverse-ended reads ", logFile);
        seqHash readHashTable(lenReadsHash2, fs::path(currSraIn.second.c_str()), ram_b);
        procRunning = false;
        constructHash.join();
        numReads = readHashTable.getNumItems();     

        seqHash mappedHash(lenReadsHash2);

        logOutput("      Splitting reads into mapped and unmapped\n", logFile);
        bamMapFile.open(fastaMap + "/Aligned.out.sam");

        // Iterate through headers in unmapped reads file
        // Fill hash tables accordingly
        for (int i = 0; i < numReadsSeqInterest + 3; i++) {
          std::getline(bamMapFile, currLine);
        }
        numMapped = 0;
        num1k = 0;
        while (std::getline(bamMapFile, currLine)) {
          readName = currLine.substr(0, currLine.find("\t"));
          if (!mappedHash.inHashTable(readName)) {
            numMapped++;
            currSeq = readHashTable.getSeq(readName);
            currHead = currSeq.get_header();
            currSeqData = currSeq.get_sequence();
            currQual = currSeq.get_quality();

            readHashTable.deleteHash(currHead);
            mappedHash.insertHash(currHead, currSeqData, currQual);
            if (numMapped % 1000 == 0) {
              num1k++;
              std::cout << "\r        Mapped read count: " +
                           std::to_string(1000 * num1k) + " ..." << std::flush;
            }
          }
        }
        percentMapped = (float(numMapped) * 100) / float(numReads);
        logOutput("\r        Mapped read count:   " + std::to_string(numMapped) +
                  " (" + percentStream.str() + "%)", logFile);
        logOutput("\n        Unmapped read count: " + std::to_string(numReads - numMapped) +
                  " (" + getPercent(100.0 - percentMapped, 2) + "%)\n", logFile);

        // Dump filled sequence hash tables to new files, containing the mapped
        // and unmapped reads respectively

        procRunning = true;
        std::thread dumpHash(progressAnim, "      Dumping split reads to mapped and unmapped files ", logFile);
        mappedHash.dump(outDir + "/" + filePrefix2 + "." + currSeqFilePrefix + ".mapped.fq");
        readHashTable.dump(outDir + "/" + filePrefix2 + ".unmapped.fq");
        procRunning = false;
        dumpHash.join();

        bamMapFile.close();
        logOutput("      Reads were extracted successfully\n", logFile);
      }
      sraRunsIn.clear();
      sra.makeCheckpoint(currSeqFilePrefix + ".iso");
    }
  }
}

// Utility function. Given vector of SRA objects, return vector of transcript
// objects constructed off of them
std::vector<transcript> get_transcript(std::vector<SRA> sras) {
  std::vector<transcript> transcripts;
  for (auto &sra : sras) {
    transcript currTrans(sra);
    transcripts.push_back(currTrans);
  }
  return transcripts;
}

// Utility function. Given a "true" / "false" string, return the corresponding
// boolean
bool stringToBool(std::string boolStr) {
  bool boolConv;
  for (int i = 0; i < boolStr.length(); i++) {
    boolStr[i] = std::tolower(boolStr[i]);
  }
  boolConv = (boolStr == "true") ? true : false;
  return boolConv;
}

// Create a checkpoint for a given assembly group
void makeGroupCheckpoint(std::string cpDir, std::string prefix) {
  std::string cpFileName = cpDir + "/" + prefix + ".trinity.ok";
  std::ofstream cpFile;
  cpFile.open(cpFileName);
  cpFile.close();
}

// Determine whether a group checkpoint file exists
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

// Create a transcript info file containing the read files used in its assembly
void makeTransInfoFile(const std::vector<std::pair<std::string, std::string>> & sraRuns, std::string transInfoFileStr) {
  std::ofstream transInfoFile;
  transInfoFile.open(transInfoFileStr);
  for (auto sra : sraRuns) {
    transInfoFile << sra.first;
    if (sra.second != "") {
      transInfoFile << " " << sra.second;
    }
    transInfoFile << std::endl;
  }
  transInfoFile.close();
}

// Perform a bulk assembly for SRAs and groups of SRAs
std::vector<std::string> run_trinity_bulk(std::map<std::string, std::vector<SRA>> sraGroups,
                                          std::string threads, std::string ram_gb,
                                          bool assembSeqsInterest, bool assembSeqsNoInterest,
                                          bool assembAllSeqs,
                                          bool dispOutput, bool retainInterFiles,
                                          std::string logFile, std::string outDir, std::string outPrefix,
                                          const INI_MAP & cfgIni = {}) {
  logOutput("\nStarting de-novo assembly\n", logFile);
  std::vector<std::string> outFiles;

  // Obtain sections from the config file as map objects
  INI_MAP_ENTRY cfgIniPipeline;
  INI_MAP_ENTRY cfgIniSeqsInt;
  INI_MAP_ENTRY cfgIniAssembGroups;
  std::string inDir;
  std::vector<std::pair<std::string, std::string>> sraRunsInTrin;
  std::vector<std::pair<std::string, std::string>> sraRunsInterest;
  std::vector<std::pair<std::string, std::string>> sraRunsNoInterest;
  std::pair<std::string, std::string> currTrinIn;
  std::string inFilePrefix;
  std::string currTrinOutAll;
  std::string currTrinOutInt;
  std::string currTrinOutNon;
  fs::path cpDir;
  fs::path trinDir;
  std::string trinOutDir;
  std::string transInfoFileStr;
  std::vector<std::string> seqsInterest;
  std::string currFilePrefix1;
  std::string currFilePrefix2;
  if (!cfgIni.empty()) {
    cfgIniPipeline = cfgIni.at("Pipeline");
    cfgIniAssembGroups = cfgIni.at("Assembly groups");
  }

  // Obtain a string vector containing paths to user-defined FASTA sequences for read extraction
  const char * home = std::getenv("HOME");
  std::string currSeqFilePath;
  for (auto seqFilePath : cfgIniSeqsInt) {
    currSeqFilePath = seqFilePath.first;
    if (currSeqFilePath[0] == '~') {
      currSeqFilePath = std::string(home) + currSeqFilePath.substr(1, currSeqFilePath.size() - 1);
    }
    seqsInterest.push_back(currSeqFilePath);
  }
  uintmax_t ram_b = (uintmax_t)stoi(ram_gb) * 1000000000;
  std::string currSeqFilePrefix;

  // Iterate over assembly groups obtained from the "Assembly Groups" section of the config file, which define
  // the SRA runs which should be assembled together, and pass their FASTQ read files into Trinity
  for (auto sraGroup : sraGroups) {
    // Define paths for checkpoint and trinity output directories
    if (!cfgIni.empty()) {
      cpDir = sraGroup.second[0].get_fastqc_dir_2().first.parent_path().parent_path().parent_path().parent_path() / ".checkpoints";
      trinDir = transcript(sraGroup.second[0]).get_trans_path_trinity().parent_path();
      trinOutDir = std::string(trinDir.c_str());
    }
    else {
      cpDir = outDir + "/.checkpoints/";
      trinOutDir = outDir + "/assembly/01-Transcript_assembly/";
    }

    inFilePrefix = sraGroup.first;

    currTrinOutNon = trinOutDir + "/" + inFilePrefix + ".unmapped.Trinity.fasta";
    currTrinOutAll = trinOutDir + "/" + inFilePrefix + ".Trinity.fasta";

    // Iterate over all SRA runs in current assembly group, preparing vectors containing input files for them,
    // as well as their mapped and unmapped read files generated during read isolation/extraction
    for (auto sra : sraGroup.second) {
      if (!cfgIni.empty()) {
        currTrinIn.first = sra.get_sra_path_orep_filt().first.c_str();
        currTrinIn.second = sra.get_sra_path_orep_filt().second.c_str();
        if (!ini_get_bool(cfgIniPipeline.at("remove_overrepresented").c_str(), 0)) {
          if (!ini_get_bool(cfgIniPipeline.at("filter_foreign_reads").c_str(), 0)) {
            if (!ini_get_bool(cfgIniPipeline.at("trim_adapter_seqs").c_str(), 0)) {
              if (!ini_get_bool(cfgIniPipeline.at("error_correction").c_str(), 0)) {
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
      }
      else {
        currTrinIn.first = sra.get_sra_path_raw().first.c_str();
        currTrinIn.second = sra.get_sra_path_raw().second.c_str();
      }
      currFilePrefix1 = fs::path(currTrinIn.first).filename().stem().stem().stem().c_str();
      currFilePrefix2 = fs::path(currTrinIn.second).filename().stem().stem().stem().c_str();
      sraRunsInTrin.push_back(currTrinIn);
      outDir = fs::path(currTrinIn.first).parent_path().c_str();
      inDir = sra.get_sra_path_orep_filt().first.parent_path().c_str();
    }
  
    // Perform assembly for all reads
    if (assembAllSeqs) {
      logOutput("\n  Now assembling all reads", logFile);
      if (!groupCheckpointExists(std::string(cpDir.c_str()), sraGroup.first)) {
        if (sraRunsInTrin.size() > 1) {
          run_trinity_comb(sraRunsInTrin, currTrinOutAll, threads, ram_gb, dispOutput, logFile);
        }
        else {
          run_trinity(sraRunsInTrin.at(0), currTrinOutAll, threads, ram_gb, dispOutput, logFile);
        }
        // Make file for entire assembly containing associated SRAs
        transInfoFileStr = std::string(fs::path(currTrinOutAll.c_str()).replace_extension(".ti").c_str());
        makeTransInfoFile(sraRunsInTrin, transInfoFileStr);

        outFiles.push_back(currTrinOutAll);
        makeGroupCheckpoint(std::string(cpDir.c_str()), sraGroup.first);
        updateHeaders(currTrinOutAll, sraGroup.first, ram_b);
      }
      else {
        logOutput("\n    Global assembly checkpoint found for: " + sraGroup.first, logFile);
      }
      sraRunsInTrin.clear();
    }
  }
  return outFiles;
}

std::vector<std::string> getCommaSepStrings(std::string commaSepStrings) {
  std::vector<std::string> stringVec;
  size_t commaInd;
  size_t currPos = 0;
  std::string currStr;
  do {
    commaInd = commaSepStrings.find(',', currPos);
    currStr = commaSepStrings.substr(currPos, commaInd - currPos);
    stringVec.push_back(currStr);
    currPos = commaInd + 1;
  } while (commaInd != std::string::npos);
  return stringVec;
}


int main(int argc, char * argv[]) {
  system("setterm -cursor off");

  INI_MAP cfgIni;
  INI_MAP_ENTRY cfgIniGen;
  INI_MAP_ENTRY cfgIniSeqsInt;
  INI_MAP_ENTRY cfgIniAssemblyGroups;
  std::string configPath;
  std::string leftReads;
  std::string rightReads;
  std::string outDir;
  std::string outPrefix;
  std::string logFilePath;
  std::vector<std::string> readFilesLeft;
  std::vector<std::string> readFilesRight;

  std::vector<SRA> sras;
  std::vector<std::string> localDataFiles;
  std::vector<std::string> seqsInterest;
  std::vector<std::string> outFiles;
  std::string threads;
  std::string ram_gb;
  bool retainInterFiles;
  bool dispOutput;
  bool entirePipeline;
  bool selectiveAssembly;
  bool assembleInterest;
  bool assembleOthers;
  bool assembleAllSeqs;
  bool compressFiles= false;
 
  if (argc == 11) {
    threads = argv[6];
    ram_gb = argv[7];
    retainInterFiles = stringToBool(argv[8]);
    dispOutput = stringToBool(argv[9]);
    entirePipeline = stringToBool(argv[10]);
    configPath = argv[1];
    selectiveAssembly = false;
    assembleInterest = false;
    assembleOthers = false;
    assembleAllSeqs = true;

    if (configPath == "null") {
      leftReads = argv[2];
      rightReads = argv[3];
      outDir = argv[4];
      outPrefix = argv[5];
      readFilesLeft = getCommaSepStrings(leftReads);
      readFilesRight = getCommaSepStrings(rightReads);
      outDir = std::string((fs::canonical(fs::path(outDir.c_str())).parent_path()).c_str()) + "/" +
               std::string((fs::canonical(fs::path(outDir.c_str())).filename()).c_str()) + "/";      
      logFilePath = outDir + "log.txt";

      make_proj_space(outDir, "assembly");

      if (readFilesLeft.size() != readFilesRight.size()) {
        logOutput("ERROR: Number of left/right read files do not match", logFilePath);
        exit(1);
      }
      for (int i = 0; i < readFilesLeft.size(); i++) {
        if (!fs::exists(readFilesLeft[i].c_str())) {
          logOutput("\nERROR: --left read file '" + readFilesLeft[i] + "' not found\n", logFilePath);
          exit(1);
        }
        if (readFilesRight[i] != "null" && !fs::exists(readFilesRight[i].c_str())) {
          logOutput("\nERROR: --right read file '" + readFilesRight[i] + "' not found\n", logFilePath);
          exit(1);
        }
        sras.push_back(SRA(readFilesLeft[i], readFilesRight[i], outDir, compressFiles, true));
      }
            
    }
    else {
      cfgIni = make_ini_map(argv[1]);
      cfgIniGen = cfgIni["General"];
      logFilePath = std::string((fs::canonical(cfgIniGen["log_file"].c_str()).parent_path() /
                                fs::path(cfgIniGen["log_file"].c_str()).filename()).c_str());
        
      if (!seqsInterest.empty()) {
        for (auto seqFilePath : seqsInterest) {
          if (!fs::exists(fs::path(seqFilePath.c_str()))) {
            logOutput("ERROR:\n  file \"" + seqFilePath + "\" not found.\n  Exiting.\n",
                      logFilePath);
            exit(1);
          }
        }
      }
      else {
        selectiveAssembly = false;
      }
      if (!selectiveAssembly) {
        if (seqsInterest.empty()) {
          //logOutput("\nUser did not define sequences of interest. Skipping selective assembly.\n",
          //          logFilePath);
        }
        assembleInterest = false;
        assembleOthers = false;
      }
    
      // If all assembly options are disabled, exit
      if (!assembleInterest && !assembleOthers && !assembleAllSeqs) {
        logOutput("No assembly option specified. Have a nice day!\n", logFilePath);
        exit(1);
      }
    
      bool compressFiles = false;
    
      // Create file space
      if (entirePipeline) {
        make_proj_space(cfgIni, "all");
      }
      else {
        make_proj_space(cfgIni, "assembly");
      }
      outDir = std::string(fs::canonical(fs::path(cfgIniGen.at("output_directory").c_str()) /
                                         fs::path(cfgIniGen.at("project_name").c_str())).c_str());

      // Obtain SRAs
      if (!dispOutput) {
        procRunning = true;
        std::thread sraRetrieve(progressAnim, "  Obtaining read information from NCBI ", logFilePath);
        sras = get_sras(cfgIni, dispOutput, compressFiles, logFilePath);
        procRunning = false;
        sraRetrieve.join();
      }
      else {
        sras = get_sras(cfgIni, dispOutput, compressFiles, logFilePath);
      }
    
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
        if (fs::exists(sraRunsLocal.first) &&
            fs::exists(sraRunsLocal.second)) {
          sras.push_back(SRA(sraRunsLocal.first, sraRunsLocal.second, cfgIni, compressFiles,
                             logFilePath));
        }
        else {
          if (sraRunsLocal.first != "" &&
              !fs::exists(sraRunsLocal.first)) {
            logOutput("ERROR: Local run not found: \"" + sraRunsLocal.first + "\"", logFilePath);
          }
          if (sraRunsLocal.second != "" &&
              !fs::exists(sraRunsLocal.second)) {
            logOutput("ERROR: Local run not found: \"" + sraRunsLocal.second + "\"", logFilePath);
          }
        }
      }

      // Check if no SRAs specified
      if (sras.empty()) {
        logOutput("ERROR: No valid SRA data.\nPlease check the configuration file.\n",
                  logFilePath);
        exit(1);
      }
    }
    // Get group specifications for SRAs
    std::map<std::string, std::vector<SRA>> sraGroups;
    if (configPath != "null") {
      cfgIniAssemblyGroups = cfgIni["Assembly groups"];
    }
    std::string currGroupName;
    std::string currIniArrStr;
    std::string currSraPrefix1;
    std::string currSraPrefix2;
    size_t numIndex1;
    size_t numIndex2;
    std::vector<std::string> iniStrArray;
    std::vector<SRA> currSraGroup;
    if (cfgIniAssemblyGroups.empty()) {
      for (auto sra : sras) {
        currSraGroup.push_back(sra);
        if (outPrefix == "") {
          outPrefix = sra.get_file_prefix().first + "_" + sra.get_file_prefix().second;
        }
        //sraGroups.emplace(outPrefix, currSraGroup);
        //currSraPrefix1 = sra.get_file_prefix().first;
        //currSraPrefix2 = sra.get_file_prefix().second;
        //sraGroups.emplace(sra.get_file_prefix().first.substr(0, 
        //                  sra.get_file_prefix().first.find_last_of("_")), currSraGroup);
        /*
        numIndex1 = std::string(fs::path(currSraPrefix1.c_str()).stem().c_str()).find("1");
        numIndex2 = std::string(fs::path(currSraPrefix2.c_str()).stem().c_str()).find("2");
        if (numIndex1 == numIndex2) {
          currSraPrefix1.erase(numIndex1 - 1, 2);
          currSraPrefix2.erase(numIndex2 - 1, 2);
        }
        sraGroups.emplace(currSraPrefix1, currSraGroup);
        else {
          sraGroups.emplace(currSraPrefix1, currSraGroup);
        }
        sraGroups.emplace(currSraPrefix.substr(0, currSraPrefix.find("_1", currSraGroup
        */
        //currSraGroup.clear();
      }
      sraGroups.emplace(outPrefix, currSraGroup);
    }
    for (auto assemblyGroup : cfgIniAssemblyGroups) {
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
                    currGroupName + "\" not found.\n  Proceeding without it.\n\n",
                    logFilePath);
        }
      }
      
        // Emplace current SRA group into map
      sraGroups.emplace(currGroupName, currSraGroup);
      currSraGroup.clear();
    }

    logOutput("\nSemblans Assemble started with following parameters:", logFilePath);
    logOutput("\n  Config file:     " + std::string(argv[1]), logFilePath);;
    logOutput("\n  Threads (Cores): " + threads, logFilePath);
    logOutput("\n  Memory (GB):     " + ram_gb, logFilePath);
    logOutput("\n  SRA Runs:\n", logFilePath);
    summarize_all_sras(sras, logFilePath, 6);
    logOutput("\n", logFilePath);
    printBreakLine(logFilePath, 6, 47);
    // Separate out reads of interest
    if (selectiveAssembly) {
      isolateReads(sras, threads, ram_gb, dispOutput,
                   cfgIni, logFilePath);
    }
    // Perform assembly with Trinity
    if (configPath != "null") {
      outFiles = run_trinity_bulk(sraGroups, threads, ram_gb, assembleInterest, 
                                  assembleOthers, assembleAllSeqs, dispOutput,
                                  retainInterFiles, logFilePath, "", "", cfgIni);
    }
    else {
      outFiles = run_trinity_bulk(sraGroups, threads, ram_gb, assembleInterest,
                                  assembleOthers, assembleAllSeqs, dispOutput,
                                  retainInterFiles, logFilePath, outDir, outPrefix);
    }

    logOutput("\nAssemble finished successfully\n\n", logFilePath);
  }
  else {
    logOutput("ERROR: Assemble invoked improperly.\n", logFilePath);
  }

  system("setterm -cursor on");
  return 0;
}
