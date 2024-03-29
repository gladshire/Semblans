#include "diamond.h"

// Given protein sequence data, use Diamond to create a blast
// database from it
void makeBlastDbDiam(std::string pathProtRef, std::string outDir,
                     bool dispOutput, std::string logFile) {
  fs::path pathProtRefFile(pathProtRef.c_str());
  std::string protRefStr(pathProtRefFile.stem().c_str());
  fs::path outDbFile((outDir + "/" + protRefStr).c_str());
  int result;
  std::string makeDbCmd = DMAKEDB + " --in " + pathProtRef + " --db " +
                          outDbFile.c_str();
  if (dispOutput) {
    makeDbCmd += " 2>&1 | tee -a " + logFile;
    logOutput("  Running command: " + makeDbCmd + "\n\n", logFile);
  }
  else {
    makeDbCmd += " >>" + logFile + " 2>&1";
  }
  result = system(makeDbCmd.c_str());
  checkExitSignal(result, logFile);
}

// Perform sequence translate - protein blastx alignment using Diamond
void blastxDiam(std::string transIn, std::string blastDb,
                std::string maxEvalue, std::string maxTargetSeqs,
                std::string threads, std::string outDir,
                bool dispOutput, std::string logFile) {
  fs::path pathTrans(transIn.c_str());
  std::string pathTransStr(pathTrans.c_str());
  std::string transStr(pathTrans.stem().c_str());

  fs::path outBlastxFile((outDir + "/" + transStr + ".blastx").c_str());
  std::string outBlastxStr(outBlastxFile.c_str());
  int result;
  std::string blastxDiamCmd = DBLASTX + " --db " + blastDb + " --query " + pathTransStr +
                              " --evalue " + maxEvalue + " --max-target-seqs " + maxTargetSeqs +
                              " --outfmt " + "6 qseqid qlen sseqid slen qframe pident nident length mismatch gapopen qstart qend sstart send evalue bitscore" +
                              " --out " + outBlastxStr +
                              " --threads " + threads +
                              " --max-target-seqs 100";
  if (dispOutput) {
    blastxDiamCmd += " 2>&1 | tee -a " + logFile;
    logOutput("  Running command: " + blastxDiamCmd + "\n\n", logFile);
  }
  else {
    blastxDiamCmd += " >>" + logFile + " 2>&1";
  }
  result = system(blastxDiamCmd.c_str());
  checkExitSignal(result, logFile);
}

// Perform protein - protein blastp alignment using Diamond
void blastpDiam(std::string pepFilePath, std::string blastDb,
                std::string maxEvalue, std::string maxTargetSeqs,
                std::string threads, std::string outFile,
                bool dispOutput, std::string logFile) {
  fs::path outBlastpFile(outFile.c_str());
  int result;
  std::string blastpDiamCmd = DBLASTP + " --query " + pepFilePath + " --db " + blastDb +
                               " --evalue " + maxEvalue + " --max-target-seqs " + maxTargetSeqs +
                               " --outfmt 6 " + " --threads " + threads + " --out " + outFile;
  if (dispOutput) {
    blastpDiamCmd += " 2>&1 | tee -a " + logFile;
    logOutput("  Running command: " + blastpDiamCmd + "\n\n", logFile);
  }
  else {
    blastpDiamCmd += " >>" + logFile + " 2>&1";
  }
  result = system(blastpDiamCmd.c_str());
  checkExitSignal(result, logFile);
}
