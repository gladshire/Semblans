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
  std::string fastaIndex;
  std::string fastaMap;
  transcript dummyTrans;
  std::ifstream bamMapFile;

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
    logOutput("  Now running read extraction using: \"" + seqsInterest[i] + "\"\n", logFile);
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
          star_map(fastaIndex, fastaMap, currSraIn, threads, dispOutput, logFile);
          procRunning = false;
          seqQuant.join();
        }
        else {
          logOutput("\n      Mapping reads to sequences of interest\n", logFile);
          star_map(fastaIndex, fastaMap, currSraIn, threads, dispOutput, logFile);
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
      seqHash readHashTable1(lenReadsHash1, fs::path(currSraIn.first.c_str()), ram_b);
      procRunning = false;
      constructHash.join();
      
      seqHash mappedHash1(lenReadsHash1);

      logOutput("      Now splitting reads into mapped and unmapped\n", logFile);
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
        numMapped++;
        readName = currLine.substr(0, currLine.find("\t"));
        if (!mappedHash1.inHashTable(readName)) {
          currSeq = readHashTable1.getSeq(readName);
          currHead = currSeq.get_header();
          currSeqData = currSeq.get_sequence();
          currQual = currSeq.get_quality();

          readHashTable1.deleteHash(readName);
          mappedHash1.insertHash(currHead, currSeqData, currQual);
        }
        if (numMapped % 1000 == 0) {
          num1k++;
          std::cout << "\r        Mapped read count: " +
                       std::to_string(1000 * num1k) + " ..." << std::flush;
        }
      }
      logOutput("\r        Mapped read count: " + std::to_string(numMapped) + "    \n", logFile);
      // Dump filled sequence hash tables to new files, containing the mapped
      // and unmapped reads respectively

      procRunning = true;
      std::thread dumpHash(progressAnim, "      Dumping split reads to mapped and unmapped files ", logFile);
      mappedHash1.dump(outDir + "/" + filePrefix1 + "." + currSeqFilePrefix + ".mapped.fq");
      readHashTable1.dump(outDir + "/" + filePrefix1 + ".unmapped.fq");
      procRunning = false;
      dumpHash.join();
      
      bamMapFile.close();

      // If SRA run is paired, perform the same steps as above for the reverse-ended FASTQ file
      if (sra.is_paired()) {
        logOutput("\n", logFile);
        numBytesReads2 = fs::file_size(currSraIn.second.c_str());
        lenReadsHash2 = numBytesReads2 / 160;
        
        procRunning = true;
        std::thread constructHash(progressAnim, "      Creating hash table from reverse-ended reads ", logFile);
        seqHash readHashTable2(lenReadsHash2, fs::path(currSraIn.second.c_str()), ram_b);
        procRunning = false;
        constructHash.join();
     
        seqHash mappedHash2(lenReadsHash2);

        logOutput("      Now splitting reads into mapped and unmapped\n", logFile);
        bamMapFile.open(fastaMap + "/Aligned.out.sam");

        // Iterate through headers in unmapped reads file
        // Fill hash tables accordingly
        for (int i = 0; i < numReadsSeqInterest + 3; i++) {
          std::getline(bamMapFile, currLine);
        }
        numMapped = 0;
        num1k = 0;
        while (std::getline(bamMapFile, currLine)) {
          numMapped++;
          readName = currLine.substr(0, currLine.find("\t"));
          if (!mappedHash2.inHashTable(readName)) {
            currSeq = readHashTable2.getSeq(readName);
            currHead = currSeq.get_header();
            currSeqData = currSeq.get_sequence();
            currQual = currSeq.get_quality();

            readHashTable2.deleteHash(currHead);
            mappedHash2.insertHash(currHead, currSeqData, currQual);
          }

          // Print number of unmapped reads identified every half million
          if (numMapped % 1000 == 0) {
            num1k++;
            std::cout << "\r        Mapped read count: " +
                         std::to_string(1000 * num1k) + " ..." << std::flush;
          }
        }
        logOutput("\r        Mapped read count: " + std::to_string(numMapped) + "    \n", logFile);

        // Dump filled sequence hash tables to new files, containing the mapped
        // and unmapped reads respectively

        procRunning = true;
        std::thread dumpHash(progressAnim, "      Dumping split reads to mapped and unmapped files ", logFile);
        mappedHash2.dump(outDir + "/" + filePrefix2 + "." + currSeqFilePrefix + ".mapped.fq");
        readHashTable2.dump(outDir + "/" + filePrefix2 + ".unmapped.fq");
        procRunning = false;
        dumpHash.join();

        bamMapFile.close();
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
void run_trinity_bulk(std::map<std::string, std::vector<SRA>> sraGroups,
                      std::string threads, std::string ram_gb,
                      bool assembSeqsInterest, bool assembSeqsNoInterest,
                      bool assembAllSeqs,
                      bool dispOutput, bool retainInterFiles,
                      std::string logFile, const INI_MAP & cfgIni) {
  logOutput("\nStarting de-novo assembly\n", logFile);

  // Obtain sections from the config file as map objects
  INI_MAP_ENTRY cfgPipeline = cfgIni.at("Pipeline");
  INI_MAP_ENTRY cfgSeqsInt = cfgIni.at("Sequences of interest");
  INI_MAP_ENTRY assembGroups = cfgIni.at("Assembly groups");
  std::string outDir;
  std::string inDir;
  std::vector<std::pair<std::string, std::string>> sraRunsInTrin;
  std::vector<std::pair<std::string, std::string>> sraRunsInterest;
  std::vector<std::pair<std::string, std::string>> sraRunsNoInterest;
  std::pair<std::string, std::string> currTrinIn;
  std::string currTrinOutAll;
  std::string currTrinOutInt;
  std::string currTrinOutNon;
  fs::path cpDir;
  fs::path trinDir;
  std::string transInfoFileStr;
  std::vector<std::string> seqsInterest;
  std::string currFilePrefix1;
  std::string currFilePrefix2;
 
  // Obtain a string vector containing paths to user-defined FASTA sequences for read extraction
  const char * home = std::getenv("HOME");
  std::string currSeqFilePath;
  for (auto seqFilePath : cfgSeqsInt) {
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
    cpDir = sraGroup.second[0].get_fastqc_dir_2().first.parent_path().parent_path().parent_path() / ".checkpoints";
    trinDir = transcript(sraGroup.second[0]).get_trans_path_trinity().parent_path();
    outDir = std::string(trinDir.c_str());

    currTrinOutNon = outDir + "/" + sraGroup.first + ".unmapped.Trinity.fasta";
    currTrinOutAll = outDir + "/" + sraGroup.first + ".Trinity.fasta";

    // Iterate over all SRA runs in current assembly group, preparing vectors containing input files for them,
    // as well as their mapped and unmapped read files generated during read isolation/extraction
    for (auto sra : sraGroup.second) {
      currTrinIn.first = sra.get_sra_path_orep_filt().first.c_str();
      currTrinIn.second = sra.get_sra_path_orep_filt().second.c_str();
      //if (assembAllSeqs) {
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
      currFilePrefix1 = fs::path(currTrinIn.first).filename().stem().stem().stem().c_str();
      currFilePrefix2 = fs::path(currTrinIn.second).filename().stem().stem().stem().c_str();
      sraRunsInTrin.push_back(currTrinIn);
      outDir = fs::path(currTrinIn.first).parent_path().c_str();
      //}
      inDir = sra.get_sra_path_orep_filt().first.parent_path().c_str();
  
      // If user preferences specify assembly of unmapped reads, push paths of unmapped read files for current
      // SRA to vector containing the inputs for assembly with Trinity    
      if (assembSeqsNoInterest) {
        currTrinIn.first = outDir + "/" + currFilePrefix1 + ".unmapped.fq";
        if (sra.is_paired()) {
          currTrinIn.second = outDir + "/" + currFilePrefix2 + ".unmapped.fq";
        }
        sraRunsNoInterest.push_back(currTrinIn);
      }

      // Iterate through all user-defined FASTA sequence files, filling vector with paths to mapped input files for
      // assembly with Trinity
      for (auto seqFilePath : seqsInterest) {
        currSeqFilePrefix = std::string(fs::path(seqFilePath.c_str()).stem().c_str());
        currTrinOutInt = outDir + "/" + sraGroup.first + "." + currSeqFilePrefix + ".mapped.Trinity.fasta";

        // Check whether assembly group checkpoint exists
        if (groupCheckpointExists(std::string(cpDir.c_str()), sraGroup.first)) {
          logOutput("Assembly found for: " + sraGroup.first, logFile);
          logOutput("  Containing:", logFile);
          summarize_all_sras(sraGroup.second, logFile, 4);
          continue;
        }

        // If user preferences specify assembly of mapped reads, push paths of mapped read files for current
        // SRA to vector containing the inputs for assembly with Trinity
        if (assembSeqsInterest) {
          currTrinIn.first = outDir + "/" + currFilePrefix1 + "." + currSeqFilePrefix + ".mapped.fq";
          if (sra.is_paired()) {
            currTrinIn.second = outDir + "/" + currFilePrefix2 + "." + currSeqFilePrefix + ".mapped.fq";
          }
          sraRunsInterest.push_back(currTrinIn);
        }

      
        // Perform assembly for reads of interest (those that map to current user-defined FASTA)
        if (assembSeqsInterest) {
          logOutput("  Now assembling reads that map to: \"" + seqFilePath + "\"\n", logFile);
          if (!groupCheckpointExists(std::string(cpDir.c_str()), sraGroup.first + "." +
                                     currSeqFilePrefix + ".mapped")) {
            if (sraRunsInterest.size() > 1) {
              run_trinity_comb(sraRunsInterest, currTrinOutInt, threads, ram_gb, dispOutput, logFile);
            }
            else {
              run_trinity(sraRunsInterest.at(0), currTrinOutInt, threads, ram_gb, dispOutput, logFile);
            }
            // Make file for mapped assembly containing associated SRAs
            transInfoFileStr = std::string(fs::path(currTrinOutInt.c_str()).replace_extension(".ti").c_str());
            makeTransInfoFile(sraRunsInterest, transInfoFileStr);
  
            makeGroupCheckpoint(std::string(cpDir.c_str()), sraGroup.first + "." +
                                currSeqFilePrefix + ".mapped");

            updateHeaders(currTrinOutInt, sraGroup.first + "_" + currSeqFilePrefix + "_mapped", ram_b);
          }
          else {
            logOutput("    Mapped assembly checkpoint found for: " + sraGroup.first + "\n\n", logFile);
          }
          sraRunsInterest.clear();
        }
      }
    }

    // Perform assembly for reads of no interest (those that DO NOT map to current user-defined FASTA)
    if (assembSeqsNoInterest) {
      logOutput("  Now assembling reads that do not map to sequences of interest\n", logFile);
      if (!groupCheckpointExists(std::string(cpDir.c_str()), sraGroup.first + ".unmapped")) {
        if (sraRunsNoInterest.size() > 1) {
          run_trinity_comb(sraRunsNoInterest, currTrinOutNon, threads, ram_gb, dispOutput, logFile);
        }
        else {
          run_trinity(sraRunsNoInterest.at(0), currTrinOutNon, threads, ram_gb, dispOutput, logFile);
        }
        // Make file for unmapped assembly containing associated SRAs
        transInfoFileStr = std::string(fs::path(currTrinOutNon.c_str()).replace_extension(".ti").c_str());
        makeTransInfoFile(sraRunsNoInterest, transInfoFileStr);

        makeGroupCheckpoint(std::string(cpDir.c_str()), sraGroup.first + ".unmapped");
        updateHeaders(currTrinOutNon, sraGroup.first + "_" + currSeqFilePrefix + "_unmapped", ram_b);
      }
      else {
        logOutput("    Unmapped assembly checkpoint found for:\n " + sraGroup.first, logFile);
      }
      sraRunsNoInterest.clear();
    }
    // Perform assembly for all reads
    if (assembAllSeqs) {
      logOutput("  Now assembling all reads", logFile);
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

        makeGroupCheckpoint(std::string(cpDir.c_str()), sraGroup.first);
        updateHeaders(currTrinOutAll, sraGroup.first, ram_b);
      }
      else {
        logOutput("    Global assembly checkpoint found for:\n " + sraGroup.first, logFile);
      }
      sraRunsInTrin.clear();
    }
  }
}


int main(int argc, char * argv[]) {
  system("setterm -cursor off");
  // Get INI config file
  INI_MAP cfgIni = make_ini_map(argv[1]);
  INI_MAP_ENTRY cfgIniGen = cfgIni["General"];

  // Define log file
  std::string logFilePath((fs::canonical((fs::path(cfgIniGen["output_directory"].c_str()))) /
                           fs::path(cfgIniGen["project_name"].c_str()) /
                           fs::path(cfgIniGen["log_file"].c_str())).c_str());
 
  if (argc > 1) {
    std::vector<SRA> sras;
    std::vector<std::string> localDataFiles;

    // Get sequences of interest from config file
    INI_MAP_ENTRY cfgSeqsInt = cfgIni.at("Sequences of interest");
    std::vector<std::string> seqsInterest;
    for (auto seqFilePath : cfgSeqsInt) {
      seqsInterest.push_back(seqFilePath.first);
    }


    // Get number of threads
    std::string threads = argv[2];

    // Get RAM in GB
    std::string ram_gb = argv[3];

    // Get boolean for intermediate file fate
    bool retainInterFiles = stringToBool(argv[4]);
    
    // Get boolean for verbose printing
    bool dispOutput = stringToBool(argv[5]);
  
    // Get sequences of interest from config
    bool selectiveAssembly = true;
    bool assembleInterest = ini_get_bool(cfgIniGen["assemble_seqs_of_interest"].c_str(), 0);
    bool assembleOthers = ini_get_bool(cfgIniGen["assemble_other_seqs"].c_str(), 0);
    bool assembleAllSeqs = ini_get_bool(cfgIniGen["assemble_all_seqs"].c_str(), 0);

        
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
        logOutput("User did not define sequences of interest. Skipping selective assembly.\n",
                  logFilePath);
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
    make_proj_space(cfgIni, "assemble");

    // Obtain SRAs
    if (!dispOutput) {
      procRunning = true;
      std::thread sraRetrieve(progressAnim, "  Obtaining read information from NCBI ", logFilePath);
      sras = get_sras(cfgIni, dispOutput, compressFiles);
      procRunning = false;
      sraRetrieve.join();
    }
    else {
      sras = get_sras(cfgIni, dispOutput, compressFiles);
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

    // Get group specifications for SRAs
    std::map<std::string, std::vector<SRA>> sraGroups;
    INI_MAP_ENTRY assemblyGroups = cfgIni["Assembly groups"];
    std::string currGroupName;
    std::string currIniArrStr;
    std::vector<std::string> iniStrArray;
    std::vector<SRA> currSraGroup;
    // Iterate through user-defined assembly groups in config file
    if (assemblyGroups.empty()) {
      for (auto sra : sras) {
        currSraGroup.push_back(sra);
        sraGroups.emplace(sra.get_file_prefix().first.substr(0, 
                          sra.get_file_prefix().first.find_last_of("_")), currSraGroup);
        currSraGroup.clear();
      }
    }
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
                    currGroupName + "\" not found.\n  Proceeding without it.\n\n",
                    logFilePath);
        }
      }
      
      // Emplace current SRA group into map
      sraGroups.emplace(currGroupName, currSraGroup);
      currSraGroup.clear();
    }

    logOutput("Semblans Assemble started with following parameters:\n", logFilePath);
    logOutput("  Config file:     " + std::string(argv[1]) + "\n", logFilePath);;
    logOutput("  Threads (Cores): " + threads + "\n", logFilePath);
    logOutput("  Memory (GB):     " + ram_gb + "\n", logFilePath);
    logOutput("  SRA Runs:\n\n", logFilePath);
    summarize_all_sras(sras, logFilePath, 6);
    // Separate out reads of interest
    if (selectiveAssembly) {
      isolateReads(sras, threads, ram_gb, dispOutput,
                   cfgIni, logFilePath);
    }
    // Perform assembly with Trinity
    run_trinity_bulk(sraGroups, threads, ram_gb, assembleInterest, 
                     assembleOthers, assembleAllSeqs, dispOutput,
                     retainInterFiles, logFilePath, cfgIni);

    logOutput("\nAssemble finished successfully\n\n", logFilePath);
  }
  else {
    logOutput("ERROR: Assemble invoked improperly.\n", logFilePath);
  }

  system("setterm -cursor on");
  return 0;
}
