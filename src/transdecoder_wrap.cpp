#include "transdecoder_wrap.h"

int MIN_FASTA = 1000;

// Determine if FASTA file meets the minimum requirement for
// number of entries
bool fasta_ok(std::string fastaFile, uintmax_t ram_b) {
  fs::path fastaFilePath(fastaFile.c_str());
  if (!fs::exists(fastaFilePath)) {
    return false;
  }
  uintmax_t numBytesTrans = fs::file_size(fastaFilePath);
  uintmax_t lenHashTable = numBytesTrans / 160;
  uintmax_t fasta_count = 0;
  seqHash fastaHash(lenHashTable, fastaFilePath, ram_b);
  fasta_count = fastaHash.getNumItems();
  if (fasta_count >= MIN_FASTA) {
    return true;
  }
  else {
    return false;
  }
}

// Determine if the blastp output file meets the minimum
// requirement for number of query sequence entries
bool blastpout_ok(std::string blastpFile) {
  fs::path blastpFilePath(blastpFile.c_str());
  if (!fs::exists(blastpFilePath)) {
    return false;
  }
  std::set<std::string> uniqueQuery;
  std::ifstream inFile(blastpFile);
  
  std::string currLine;
  std::string currQuery;
  size_t wsPos = 0;
  while(std::getline(inFile, currLine)) {
    wsPos = currLine.find('\t');
    if (wsPos != std::string::npos) {
      currQuery = currLine.substr(0, wsPos);
      uniqueQuery.insert(currQuery);
    }
  }
  if (uniqueQuery.size() >= MIN_FASTA) {
    return true;
  }
  else {
    return false;
  }
}

// Given a transcripts file, predict the sequence data's coding regions using
// TransDecoder
void run_transdecoder(std::string transIn, std::string transCds, std::string transPep,
                      bool useBlast, std::string maxEvalue, std::string maxTargetSeqs,
                      std::string threads, uintmax_t ram_b,
                      std::string dbPath, std::string outDir,
                      bool dispOutput, std::string logFile) {
  std::string transFilePath(transIn);
  fs::path cdsFilePath(transCds);
  fs::path pepFilePath(transPep);

  std::string transPrefix = fs::path(transIn.c_str()).filename().stem().stem().c_str();

  std::string blastpout = transPrefix + ".blastp.outfmt6";
  std::string printOut;
  if (dispOutput) {
    printOut = " 2>&1 | tee -a " + logFile;
  }
  else {
    printOut = " >>" + logFile + " 2>&1";
  }
  if (fasta_ok(std::string(cdsFilePath.c_str()), ram_b) &&
      fasta_ok(std::string(pepFilePath.c_str()), ram_b)) {
    // Skip TransDecoder
  }
  else {
    std::string transFileName = std::string(fs::path(transFilePath.c_str()).filename().c_str());
    fs::path allpep(std::string(outDir + "/" + transFileName +
                    ".transdecoder_dir/longest_orfs.pep").c_str());
    if (fs::exists(allpep)) {
      //logOutput("Skipping search for long orfs", logFile);
    }
    else {
      // Only operates on un-stranded. Implement stranded later
      std::string tdLongOrfs_cmd = PATH_TRANSD_LONGORFS + " -t " + transFilePath + " -O " +
                                   std::string(allpep.parent_path().c_str()) + printOut;
      int resultLO;

      std::cout << tdLongOrfs_cmd << std::endl;

      resultLO = system(tdLongOrfs_cmd.c_str());
      if (WIFSIGNALED(resultLO)) {
        system("setterm -cursor on");
        logOutput("Exited with signal " + std::to_string(WTERMSIG(resultLO)), logFile);
        exit(1);
      }
      
    }
    if (blastpout_ok(outDir + "/" +  blastpout)) {
      logOutput("Skipping blastp", logFile);
    }
    else {
      if (useBlast) {
        blastp(std::string(allpep.c_str()), dbPath, maxEvalue, maxTargetSeqs, threads,
               outDir + "/" + blastpout, dispOutput, logFile);
      }
      else {
        blastpDiam(std::string(allpep.c_str()), dbPath, maxEvalue, maxTargetSeqs, threads,
                   outDir + "/" + blastpout, dispOutput, logFile);
      }
    }
    if (fasta_ok(std::string(cdsFilePath.c_str()), ram_b) &&
        fasta_ok(std::string(pepFilePath.c_str()), ram_b)) {
      //logOutput("Skip finding final CDS and PEP", logFile);
    }
    else {
      std::string tdPredict_cmd = PATH_TRANSD_PREDICT + " -t " + transFilePath +
                                  " --retain_blastp_hits " + outDir + "/" + blastpout +
                                  " --cpu " + threads + printOut;
      fs::path currDir = fs::current_path();
      fs::current_path(fs::path(outDir.c_str()));
      int resultPD;

      resultPD = system(tdPredict_cmd.c_str());
      if (WIFSIGNALED(resultPD)) {
        system("setterm -cursor on");
        logOutput("Exited with signal " + std::to_string(WTERMSIG(resultPD)), logFile);
        exit(1);
      }
      fs::current_path(currDir);
      std::rename((outDir + "/" + transFileName + ".transdecoder.gff3").c_str(),
                  (outDir + "/" + transPrefix + ".transdecoder.gff3").c_str());
      std::rename((outDir + "/" + transFileName + ".transdecoder.bed").c_str(),
                  (outDir + "/" + transPrefix + ".transdecoder.bed").c_str());
      std::rename((outDir + "/" + transFileName + ".transdecoder.pep").c_str(),
                  (outDir + "/" + transPrefix + ".transdecoder.pep.fasta").c_str());
      std::rename((outDir + "/" + transFileName + ".transdecoder.cds").c_str(),
                  (outDir + "/" + transPrefix + ".transdecoder.cds.fasta").c_str());
    }
  }
}
