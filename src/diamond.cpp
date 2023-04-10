#include "diamond.h"

void makeDb(std::string pathProtRef, std::string outDir,
            bool dispOutput, std::string logFile) {
  fs::path pathProtRefFile(pathProtRef.c_str());
  std::string protRefStr(pathProtRefFile.stem().c_str());
  fs::path outDbFile((outDir + "/" + protRefStr).c_str());
  std::string printOut;
  if (dispOutput) {
    printOut = " 2>&1 | tee -a " + logFile;
  }
  else {
    printOut = " >>" + logFile + " 2>&1";
  }
  int result;
  if (fs::exists(outDbFile)) {
    logOutput("BLAST database found for: " + protRefStr, logFile);
    return;
  }
  std::string makeDbCmd = DMAKEDB + " --in " + pathProtRef + " --db " +
                          outDbFile.c_str() + printOut;
  result = system(makeDbCmd.c_str());
  if (WIFSIGNALED(result)) {
    logOutput("Existed with signal " + std::to_string(WTERMSIG(result)), logFile);
    exit(1);
  } 
}

void blastxDiam(std::string transIn, std::string blastDb,
                std::string threads, std::string outDir,
                bool dispOutput, std::string logFile) {
  fs::path pathTrans(transIn.c_str());
  std::string pathTransStr(pathTrans.c_str());
  std::string transStr(pathTrans.stem().c_str());

  fs::path outBlastxFile((outDir + "/" + transStr + ".blastx").c_str());
  std::string outBlastxStr(outBlastxFile.c_str());
  std::string printOut;
  if (dispOutput) {
    printOut = " 2>&1 | tee -a " + logFile;
  }
  else {
    printOut = " >>" + logFile + " 2>&1";
  }
  int result;
  if (fs::exists(outBlastxFile)) {
    logOutput("Diamond BLASTX output found for: " + transStr, logFile);
    return;
  }
  std::string blastxDiamCmd = DBLASTX + " --db " + blastDb + " --query " + pathTransStr +
                              " --evalue " + "0.01" +
                              " --outfmt " + "6 qseqid qlen sseqid slen qframe pident nident length mismatch gapopen qstart qend sstart send evalue bitscore" +
                              " --out " + outBlastxStr +
                              " --threads " + threads +
                              " --max-target-seqs 100" +
                              printOut;
  result = system(blastxDiamCmd.c_str());
  if (WIFSIGNALED(result)) {
    logOutput("Exited with signal " + std::to_string(WTERMSIG(result)), logFile);
    exit(1);
  }
}

void blastpDiam(std::string pepFilePath, std::string blastDb,
                std::string threads, std::string outFile,
                bool dispOutput, std::string logFile) {
  fs::path outBlastpFile(outFile.c_str());
  std::string printOut;
  if (dispOutput) {
    printOut = " 2>&1 | tee -a " + logFile;
  }
  else {
    printOut = " >>" + logFile + " 2>&1";
  }
  int result;
  if (fs::exists(outBlastpFile)) {
    logOutput("Diamond BLASTP output found for: " + pepFilePath, logFile);
    return;
  }
  std::string  blastpDiamCmd = DBLASTP + " --query " + pepFilePath + " --db " + blastDb +
                               " --max-target-seqs 1 --outfmt 6 --evalue 10 " +
                               " --threads " + threads + " --out " + outFile + printOut;
  result = system(blastpDiamCmd.c_str());
  if (WIFSIGNALED(result)) {
    logOutput("Exited with signal " + std::to_string(WTERMSIG(result)), logFile);
    exit(1);
  }
}
